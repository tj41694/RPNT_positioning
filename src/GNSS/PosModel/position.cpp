#include "GNSS/PosModel/position.h"
#include "BaseFunction/basefunction.h"

/* constant --------------------------------------------------------------------------------------- */
#define SQR(x)              ((x)*(x))
#define SQRT(x)             ((x)<0.0?-sqrt(-x):sqrt(x))

#define VAR_CLK_OFF1        3600.0          /* initial variance of pre-estimated receiver clock offset (m^2) */
#define VAR_CLK_DRF1        3600.0          /* initial variance of pre-estimated receiver clock drift (m^2/s^2) */
#define VAR_CLK_OFF2        9E12	        /* initial variance of receiver clock offset (m^2) */
#define VAR_CLK_DRF2        1E6             /* initial variance of receiver clock drift (m^2/s^2) */

#define CJUMP_K1            0.9	            /* clock jump threshold k1 to confirm clock jump (ms) */
#define CJUMP_K2            1E-3            /* clock jump threshold k2 to confirm exact clock jump (ms) */
#define CJUMP_K3            5E-2            /* clock jump threshold k3 to confirm rough clock jump (ms) */
#define CJUMP_UNIT          1E-3            /* clock jump unit (s to ms) */
#define CJUMP_MINSAT        3               /* minimum number of satellites to confirm clock jump */

#define MAXITR              10              /* max number of iteration for point pos */
#define ERR_CBIAS           0.3             /* code bias error std (m) */
#define RES_RATIO           3.0			    /* maximum residual ratio */
#define MAX_POS_VAR         900             /* maximum 3D position formal variance = 30^2 (m^2) */

#define MAX_BIAS_DT         2592000         /* maximum bias time difference = 30 days (s) */

#define ION_RES_THRESHOLD   9.0             /* ionospheric residual threshold (m) */
#define MAX_ION_RES_RATIO   5.0             /* maximum ionospheric residual ratio to average */

const double ION_RES_NORM=SQR(GNSS_LAMBDA_M[iGPS][1])-SQR(GNSS_LAMBDA_M[iGPS][2]);
                                            /* normalization factor of ionosphere residual */
const string SYS_LIST="GREC";			    /* system flag GPS, GLONASS, Galileo, BeiDou */

/* position class --------------------------------------------------------------------------------- */
/* gnss single point position (SPP) class ------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
gnss_single_c::gnss_single_c() {
    /* processing status */
    nsol=ntotal=0;
    stat=vstat=0;
    numsys=ns=0;
    resflag=ratflag=niter=staflag=res_ave=0;
    lsqflag=0;
    posVar=0;
    dsant=0,drant=0;
    satnum=0;
    for (int i=0; i<NSYS; i++) {
        ncode[0][i]=ncode[1][i]=ndoppler[0][i]=ndoppler[1][i]=0;
        Pdop[i]=SNR[i]=multiPath[i]=0;
        preClkOff[i]=0;
    }
    preClkDrf=0;
    preClkVar[0]=VAR_CLK_OFF1; preClkVar[1]=VAR_CLK_DRF1;
    clkRefSys=clkSys1st=-1;
    Pdop[NSYS]=SNR[NSYS]=multiPath[NSYS]=0.0;
    for (int i=0; i<NFREQ; i++) {
        CJ_float[i]=CJ_corr[i]=0;
        CJ_nsat[i]=CJ_ndetect[i]=0;
        CJ_sum[i]=0;
        CJ_con[i]=0;
        CJ_n100[i].assign(100,0);
    }
    n_elAmb=fix_num=0;
    fixNumRatio=fixVarRatio=0;
    dtime=0;
    lastTime=gtime_c();
}
gnss_single_c::~gnss_single_c() {
    rp=NULL; bp=NULL;
    lam=lam1=lam2=NULL;

    usesat.clear(); rovsat.clear(); bassat.clear();

    for (int i=0; i<NFREQ; i++) {
        CJ_n100[i].clear();
        iobs[i].clear(); iobs[i+NFREQ].clear();
    }
    Lobs.clear(); Acoe.clear(); Rvar.clear();
    Rvec.clear(); Xpar.clear(); Xori.clear();
    Rx.clear();   QAA.clear();

    fixIndex.clear();
    ambnum.clear();
    ambElevation[0].clear(); ambElevation[1].clear();
    fix_amb.clear(); sum_var.clear();
    fix_flag.clear();
    fix_Xpar.clear(); fix_Rx.clear();
}
/* Implementation functions ------------------------------------------------------- */
/* correct differenced code bias -----------------------------------------------------
* notice: 
*   1. for BDS-3,   B2a(C5), B2b(C7) and B2a+b(C8) can use the same DCB, e.g. B2-B3
*   2. for BDS-3,   B1A/C(C1) and B2(C2) can use the same DCB, e.g. B1-B3. But the
*                   error can reach to 8 ns (2 - 3 m)
*   3. for Galileo, E5a(E5), E5b(E7) and E5a+b(E8) can use the same DCB, e.g. E5a-E1
* --------------------------------------------------------------------------------- */
void gnss_single_c::correct_DCB() {
    /* post bias solution ok flag */
    bool biasOk=false;
    if ( nav->obsbias.cbts.timediff(obsp->data[0].time)<MAX_BIAS_DT &&
         obsp->data[0].time.timediff(nav->obsbias.cbte)<MAX_BIAS_DT ) {
        biasOk=true;
    }

    for (int i=0; i<obsp->n; i++) {
        gnss_obsd_c *obs=&obsp->data[i];

        /* if have bias product values */
        if ( biasOk && correct_code_bias(obs,nav,opt)>0 ) {
            continue;
        }
        /* or use broadcast tgd instead */
        correct_code_tgd(obs,nav,opt);
    }
}
/* update consecutive number for each satellite ----------------------------------- */
void gnss_single_c::update_sat_nCon() {
    if ( obsp==obsr && obsp->n>0 ) {
        if (obsp->data[0].time.timediff(lastTime)>opt->restime+DTTOL) {
            nsol=0;
        }

        int good=NF;
        for (int i=0; i<obsp->n; i++) {
            /* reset consecutive count number */
            if (obsp->data[i].time.timediff(ssat[obsp->data[i].sat-1].obstime.back())>opt->restime+DTTOL) {
                ssat[obsp->data[i].sat-1].nCon=0;
            }

            good=NF;
            /* test if all required observations are available */
            /*for (int iFreq=0; iFreq<NF; iFreq++) {
                if ( obsp->data[i].P[iFreq]==0.0 ) {
                    good--;
                }
            }
            if (good>0) obsp->data[i].nCon=ssat[obsp->data[i].sat-1].nCon+1;*/
            for (int iFreq=0; iFreq<NF; iFreq++) {
                if ( obsp->data[i].P[iFreq]==0.0 || ( usePhase && obsp->data[i].L[iFreq]==0.0 ) ) {
                    good--;
                    break;
                }
            }
            if (good==NF) obsp->data[i].nCon=ssat[obsp->data[i].sat-1].nCon+1;
        }
        for (int i=0; i<obsp->n; i++) {
            ssat[obsp->data[i].sat-1].nCon=obsp->data[i].nCon;
        }

        lastTime=obsp->data[0].time;
    }
}
/* initialize ephemeris and DCB --------------------------------------------------- */
int gnss_single_c::gnss_ini_EPH_DCB() {
    solp->stat=SOLQ_NONE;

    if (obsp->n<1) { obsp->errorMsg="no observation data\n"; return 0; }

    solp->time=obsp->data[0].time; obsp->errorMsg="";

    /* initialize observation (exc) */
    for (int i=0; i<obsp->n; i++) obsp->data[i].exc=0;

    /* correct differenced code bias */
    correct_DCB();

    /* satellite positons, velocities and clocks */
    if (satfunc->satposclk(obsp,nav)<=4) {
        obsp->errorMsg="no enough available emphemeris!\n ";
        return 0; 
    }

    return 1;
}
/* detect wrong observation with double frequency --------------------------------- */
void gnss_single_c::detect_wrong_obs_2freq() {
    /* detect wrong observation with double frequency */
    vector<double> ion_res;
    vector<int> index_obs;
    int i,ipos,nion=0,nave=0;
    double res_ave;

    for (i=ns=0; i<obsp->n; i++) {
        /* excluded satellite? */
        if (obsp->data[i].exc>0) continue;
        if (exclude_sat(obsp->data[i])) {
            obsp->data[i].exc=1;  continue;
        }
        if (obsp->data[i].P[0]==0.0||obsp->data[i].P[1]==0.0) continue;
        lam=nav->lam[obsp->data[i].sat-1];	/* lambda */
        if (lam[0]==0.0||lam[1]==0.0) continue;
        double norm_ion_res =
            fabs ( (obsp->data[i].P[0]-obsp->data[i].P[1]) / (SQR(lam[0])-SQR(lam[1])) * ION_RES_NORM );
        /* arrange ion_res from small to large */
        for (ipos=0; ipos<ion_res.size(); ipos++) {
            if (ion_res[ipos]>norm_ion_res) break;
        }
        ion_res.insert(ion_res.begin()+ipos,norm_ion_res);
        index_obs.insert(index_obs.begin()+ipos,i);
    }

    /* initialize average residual with the smallest half values */
    nion=ion_res.size();
    if (nion<4) return;
    nave=round(nion/2.0);
    res_ave=ave_vec_n0(ion_res,nave);

    /* detect the rest of observations */
    for (i=nave;i<nion;i++) {
        // if find large ionospheric residual
        if ( ion_res[i]>ION_RES_THRESHOLD && ion_res[i]>MAX_ION_RES_RATIO*res_ave ) {
            obsp->data[ index_obs[i] ].errorMsg="large ionospheric residual!";
            obsp->data[ index_obs[i] ].exc=1;
            continue;
        }
        // if ionospheric residual is OK
        res_ave+=(ion_res[i]-res_ave)/(++nave);
    }
}

/* reset satellite status --------------------------------------------------------- */
void gnss_single_c::resetsat() {
    for (int i=0; i<MAXSAT; i++) {
        for (int j=0; j<NFREQ; j++)
            ssat[i].vsat[j]=ssat[i].snr[j]=ssat[i].multiPath[j]=0;
        ssat[i].multiPath[NFREQ]=0.0;
    }
}
/* carrier-phase bias (fcb) correction -------------------------------------------- */
void gnss_single_c::corr_phase_bias() {
    double lambda;
    int code;

    /* rover observation */
    for (int i=0; i<obsr->n; i++) for (int j=0; j<NFREQ; j++) {
        if (!(code=obsr->data[i].code[j])) continue;
        if ((lambda=nav->lam[obsr->data[i].sat-1][j])==0.0) continue;
        /* correct phase bias (cyc) */
        obsr->data[i].L[j]=nav->ssr[obsr->data[i].sat-1].pbias[code-1]/lambda;
    }
    /* base observation */
    for (int i=0; i<obsb->n; i++) for (int j=0; j<NFREQ; j++) {
        if (!(code=obsb->data[i].code[j])) continue;
        if ((lambda=nav->lam[obsb->data[i].sat-1][j])==0.0) continue;
        /* correct phase bias (cyc) */
        obsb->data[i].L[j]=nav->ssr[obsb->data[i].sat-1].pbias[code-1]/lambda;
    }

}
/* select common satellite of rover and base -------------------------------------- */
void gnss_single_c::selcomsat() {
    usesat.clear(); rovsat.clear(); bassat.clear();

    /* select common satellites */
    for (int i=0,j=0; i<obsr->n&&j<obsb->n; i++,j++) {
        if (obsr->data[i].sat<obsb->data[j].sat) j--;
        else if (obsr->data[i].sat > obsb->data[j].sat) i--;
        else {
            usesat.push_back(obsr->data[i].sat);
            rovsat.push_back(i);
            bassat.push_back(j);
        }
    }
    satnum=usesat.size();
}
/* update observation and ambiguity information for each satellites --------------- */
void gnss_single_c::update_obs_amb_ssat() {
    /* select common satellites for RTK */
    if (opt->mode>=PMODE_DGPS) selcomsat();
    else {
        usesat.clear(); rovsat.clear();
        for (int i=0; i<obsr->n; i++) {
            usesat.push_back(obsr->data[i].sat);
            rovsat.push_back(i);
        }
        satnum=usesat.size();
    }

    /* initialize clock jump detection */
    for (int ifreq=0; ifreq<CJ_nfreq; ifreq++) {
        CJ_float[ifreq]=CJ_corr[ifreq]=0;
        CJ_nsat[ifreq]=0;
        CJ_ndetect[ifreq]=satnum;
    }

    gnss_ssat_c *sss;
    for (int i=0; i<satnum; i++) {
        sss=&ssat[ usesat[i]-1 ];
        lam=nav->lam[ sss->sat-1 ];
        rp=&obsr->data[rovsat[i]]; 
        if (opt->mode>=PMODE_DGPS) bp=&obsb->data[bassat[i]];
        else bp=NULL;

        /* update vectors for current status */
        sss->input_new_obs(rp,bp,CJ_sum);
        /* reset flag according to unsolved time interval */
        sss->test_reset(opt);

        /* detect clock jump */
        for (int ifreq=0; ifreq<CJ_nfreq; ifreq++) {
            if ( sss->reset[ifreq] || sss->amb_obsvb[ifreq][0]!=OBSTYPE_CL || sss->amb_obsvb[ifreq][1]!=OBSTYPE_CL ) {
                CJ_ndetect[ifreq]--;
                continue;
            }
            /* clock jump discriminant value */
            CJ_SF=( sss->P[ifreq].back()-sss->P[ifreq][MAXMEASBUF-1] -
                (sss->L[ifreq].back()-sss->L[ifreq][MAXMEASBUF-1])*lam[ifreq] ) / CLIGHT/CJUMP_UNIT;
            if (CJ_SF>CJUMP_K1) {
                CJ_nsat[ifreq]++;
                CJ_float[ifreq]+=CJ_SF;
            }
        }
    }

    /* determine clock jump for every frequency before update ambiguity */
    clock_jump_everyFreq();

    /* detect cycle slip and update ambiguity statistic */
    for (int i=0; i<satnum; i++) {
        sss=&ssat[ usesat[i]-1 ];
        lam=nav->lam[ sss->sat-1 ];
        rp=&obsr->data[rovsat[i]]; 
        if (opt->mode>=PMODE_DGPS) bp=&obsb->data[bassat[i]];
        else bp=NULL;

        /* update process parameter for current epoch */
        sss->update_pro_par(CJ_corr);
        /* detect cycle slip */
        sss->detect_slip(rp,bp,opt);
        sss->test_reset_LC();
        /* update ambiguity statistic */
        sss->update_amb_statistic(opt);
    }
}
/* test SNR mask ---------------------------------------------------------------------
* test SNR mask
* args   : int    base      I   rover or base-station (0:rover,1:base station)
*          int    freq      I   frequency (0:L1,1:L2,2:L3,...)
*          double el        I   elevation angle (rad)
*          double snr       I   C/N0 (dBHz)
*          snrmask_t *mask  I   SNR mask
* return : status (1:masked,0:unmasked)
* --------------------------------------------------------------------------------- */
int gnss_single_c::test_snr(int rovbas,int freq,double el,double snr,
    const snrmask_t *mask) {
    double minsnr,a;
    int i;

    if (!mask->ena[rovbas]||freq<0||freq>=NFREQ) return 0;

    a=(el*R2D+5.0)/10.0;
    i=(int)floor(a); a-=i;
    if (i<1) minsnr=mask->mask[freq][0];
    else if (i>8) minsnr=mask->mask[freq][8];
    else minsnr=(1.0-a)*mask->mask[freq][i-1]+a*mask->mask[freq][i];

    return snr<minsnr;
}
/* pseudorange measurement error variance ----------------------------------------- */
double gnss_single_c::varerr(const gnss_obsd_c &obs) {
    double fact;
    if (obs.sys!=SYS_BDS) fact=obs.sys==SYS_GLO ? EFACT_GLO : (obs.sys==SYS_SBS ? EFACT_SBS : EFACT_GPS);
    else fact = binary_search(BDS_GEO,BDS_GEO+NBDS_GEO,obs.prn)? EFACT_BDS_G : EFACT_BDS;

    /* elevation factor */
    double sin_elR = obs.azel[1]<opt->con_weight_el? sin(obs.azel[1]) : 1;

    /* IF factor */
    //if (opt->sppiono==IONOOPT_IFLC) fact*=3.0;

    return SQR(fact) * ( SQR(opt->err[0]) + SQR(opt->err[1]/sin_elR) );
}
/* verify excluded satellite ------------------------------------------------------ */
int gnss_single_c::exclude_sat(gnss_obsd_c &obs) {
    /* ephemeris unavailable */
    if (obs.svh!=0) {
        obs.errorMsg="exclude unhealthy satellite!";
        return 1;
    }
    /* BDS */
    if (obs.sys==SYS_BDS){
        /* BDS generation */
        if ( !(opt->bdsGen&BDS_3RD) && obs.prn>=BDS3_MIN ) {
            obs.errorMsg="exclude BDS-3 satellite!";
            return 1;
        }
        else if ( !(opt->bdsGen&BDS_2ND) && obs.prn<BDS3_MIN ) {
            obs.errorMsg="exclude BDS-2 satellite!";
            return 1;
        }
        /* BDS orbit */
        if ( !(opt->bdsOrbit&BDS_ORB_MEO) && binary_search(BDS_MEO,BDS_MEO+NBDS_MEO,obs.prn) ) {
            obs.errorMsg="exclude BDS MEO satellite!";
            return 1;
        }
        else if ( !(opt->bdsOrbit&BDS_ORB_IGSO) && binary_search(BDS_IGSO,BDS_IGSO+NBDS_IGSO,obs.prn) ) {
            obs.errorMsg="exclude BDS IGSO satellite!";
            return 1;
        }
        else if ( !(opt->bdsOrbit&BDS_ORB_GEO) && binary_search(BDS_GEO,BDS_GEO+NBDS_GEO,obs.prn) ) {
            obs.errorMsg="exclude BDS GEO satellite!";
            return 1;
        }
    }

    /* exclude satellite according to option */
    if (opt->exsats[obs.sat-1]==1) {
        obs.errorMsg="exclude satellite by option!";
        return 1;
    }
    /* unselected sat sys */
    if (!(obs.sys&opt->navsys)) {
        obs.errorMsg="exclude satellite of unselected constellation!";
        return 1;
    }

    /* inconsecutive satellite */
    if ( obsp==obsr && nsol+2>2*opt->minConsecutive && obs.nCon<opt->minConsecutive ) {
        obs.errorMsg="exclude inconsecutive satellite!";
        return 1;
    }

    if (obs.sys==SYS_QZS) obs.svh&=0xFE; /* mask QZSS LEX health */

    return 0;
}
/* exclude constellations with only a few available satellites -------------------- */
void gnss_single_c::exc_few_sat() {
    for (int i=0; i<obsp->n; i++) {
        if (obsp->data[i].exc>0||!obsp->data[i].used) continue;
        if (ncd[obsp->data[i].isys]<MINSYSOBS) {
            obsp->data[i].exc=1;
            obsp->data[i].used=0;
            obsp->data[i].errorMsg="exclude constellations with only 1 available satellites!";
        }
    }
}
/* exclude observation with large residual ---------------------------------------- */
void gnss_single_c::exc_largeres(const int last_f) {
    /* exclude large residual observation after the first ajustment */
    if (last_f==1) for (int i=0; i<obsp->n; i++) {
        if (obsp->data[i].exc>0||!obsp->data[i].used) continue;
        if (fabs(obsp->data[i].res[0])>opt->maxres) {
            ncd[obsp->data[i].isys]--;
            obsp->data[i].exc=1;
            obsp->data[i].used=0;
            obsp->data[i].errorMsg="high Pseudorange residual!";
        }
    }
    /* exclude large residual during the process */
    else if (last_f==0) {
        ratflag=0;
        int n_ave=0;
        res_ave=0.0;
        //average residual
        for (int i=0; i<obsp->n; i++) {
            if (obsp->data[i].exc>0||!obsp->data[i].used) continue;
            res_ave+=fabs(obsp->data[i].res[0]);
            n_ave++;
        }
        res_ave/=n_ave;
        //find large residual
        if (res_ave>3.0*opt->maxres) for (int i=0; i<obsp->n; i++) {
            if (obsp->data[i].exc>0||!obsp->data[i].used) continue;
            if (fabs(obsp->data[i].res[0])/res_ave>RES_RATIO) {
                ratflag=1;
                ncd[obsp->data[i].isys]--;
                obsp->data[i].exc=1;
                obsp->data[i].used=0;
                obsp->data[i].errorMsg="large Pseudorange residual ratio!";
            }
        }
    }
}
/* correct and save pre-etimated clock offset and drift --------------------------- */
void gnss_single_c::process_pre_clock() {
    /* clock offset */
    /* update observations with clock offset */
    int iLobs=0;
    for (int i=0; i<obsp->n; i++) {
        if (obsp->data[i].used) {
            Lobs[iLobs]-=preClkOff[obsp->data[i].isys];
            iLobs++;
        }
    }
    /* update solution with clock offset (for PPP) */
    solp->xclk[iGPS]=preClkOff[iGPS]; /* receiver clock bias (s) */
    solp->xclk[iGLO]=preClkOff[iGLO]; /* glo-gps time offset (s) */
    solp->xclk[iGAL]=preClkOff[iGAL]; /* gal-gps time offset (s) */
    solp->xclk[iBDS]=preClkOff[iBDS]; /* bds-gps time offset (s) */

    /* clock drift */
    if ( opt->dynamics ) {
        /* update velocity observations with clock drift */
        for (int i=0; i<vnumL; i++) vLobs[i]-=preClkDrf;
    }
}
/* clock jump detection and correction (for every frequency) ---------------------- */
int gnss_single_c::clock_jump_everyFreq() {
    /* check and repair clock jump for every frequency */
    for (int ifreq=0; ifreq<CJ_nfreq; ifreq++) {
        int id_clock=0;
        //continue if not all satellites was detected "ms" clock jump
        if ( CJ_nsat[ifreq]<CJ_ndetect[ifreq] || CJ_nsat[ifreq]<CJUMP_MINSAT ) continue; 
        /* get average clock jump and its int value */
        CJ_float[ifreq]/=CJ_nsat[ifreq]; id_clock=int(CJ_float[ifreq]);
        double diff_cj=fabs(CJ_float[ifreq]-id_clock);
        /* if it's a exact "ms" clock jump */
        if ( diff_cj<=CJUMP_K2 ) CJ_corr[ifreq]=(CJUMP_UNIT*id_clock*CLIGHT);
        /* if it's a rough "ms" clock jump */
        else if ( diff_cj<=CJUMP_K3 ) CJ_corr[ifreq]=(CJUMP_UNIT*CJ_float[ifreq]*CLIGHT);
        else CJ_corr[ifreq]=0;
        /* update clock jump statistical data */
        CJ_con[ifreq]++;
        CJ_sum[ifreq]+=CJ_corr[ifreq];
        CJ_n100[ifreq].erase(CJ_n100[ifreq].begin());
        CJ_n100[ifreq].push_back(id_clock);
    }

    return 1;
}
/* psendorange with code bias correction and combination -------------------------- */
double gnss_single_c::prange(gnss_obsd_c &obs,const int iter) {
    lam=nav->lam[obs.sat-1];
    double PC,P1,P2,gamma; /* f1^2/f2^2 */
    obs.dcbvar=0.0;

    if (lam[0]==0.0||lam[1]==0.0) return 0.0;
    gamma=SQR(lam[1])/SQR(lam[0]);

    /* test snr mask */
    if (iter>0) {
        if (test_snr(staflag,0,obs.azel[1],obs.SNR[0],&opt->snrmask)) {
            return 0.0;
        }
        if (opt->sppiono==IONOOPT_IFLC) {
            if (test_snr(staflag,1,obs.azel[1],obs.SNR[1],&opt->snrmask)) return 0.0;
        }
    }

    P1=obs.P[0]; P2=obs.P[1];

    if (opt->sppiono==IONOOPT_IFLC) { /* dual-frequency */
        if (P1==0.0||P2==0.0) return 0.0;
        /* iono-free combination */
        PC=(gamma*P1-P2)/(gamma-1.0);
    }
    else { /* single-frequency */
        if (P1==0.0) return 0.0;
        PC=P1;
    }

    obs.dcbvar=SQR(ERR_CBIAS); /* code bias error std added to ovar */

    return PC;
}
/* compute observation\coefficient\covariance matrix ------------------------------ */
int gnss_single_c::code_matrix() {
    numL=0;		    /* number of used observation (size of L) */
    numX=3;         /* number of estimated parameters */
    double blh[3];	/* geodetic position (lat,lon,h) */

    /* reset each matrix (except X=[NSYS+3,1]) */
    Acoe.clear(); Lobs.clear(); Rvec.clear();
    for (int i=0; i<NSYS; i++) ncd[i]=0; //number of code observation of each system
    vector<double> Abuff;

    /* ecef (XYZ) to geodetic position */
    ecef2blh(Xpar.begin(),WGS84,blh);

    for (int i=ns=0; i<obsp->n; i++) {
        double Pc;		/* psudorange with code bias correction */
        double dion=0,ionvar=0;
        double lam_L1;	/* f1 lambda */

        /* reset observation data */
        obsp->data[i].used=0;

        /* excluded satellite? */
        if (obsp->data[i].exc>0) continue;
        if (exclude_sat(obsp->data[i])) {
            obsp->data[i].exc=1;  continue;
        }

        if (!obsp->data[i].sys) {
            obsp->data[i].exc=1;
            obsp->data[i].errorMsg="unknown satellite!";
            continue;
        }

        /* reject duplicated observation data */
        if (i<obsp->n-1&&obsp->data[i].sat==obsp->data[i+1].sat) {
            obsp->data[i].errorMsg="duplicated observation data with next data!";
            i++;
            continue;
        }
        /* geometric distance/azimuth/elevation angle */
        if (
            (obsp->data[i].dist=geodist_GNSS(obsp->data[i].posvel,Xpar.begin(),obsp->data[i].sigvec))<300000.0 ||
            ( satazel(blh,obsp->data[i].sigvec,obsp->data[i].azel) < opt->min_elevation && lsqflag>0 )
            ) {
            obsp->data[i].errorMsg="low elevation!";
            continue;
        }

        /* pseudorange with code bias correction */
        if ((Pc=prange(obsp->data[i],niter))==0.0) {
            obsp->data[i].errorMsg="no avaliable pseudorange!";
            continue;
        }

        /* ionospheric corrections */
        if (resflag&&!sppionf->correction(&obsp->data[i],nav,blh)) {
            obsp->data[i].errorMsg="ionosphere correction error!";
            continue;
        }

        /* GPS-L1 -> L1/B1 */
        if (opt->sppiono!=IONOOPT_IFLC&&(lam_L1=nav->lam[obsp->data[i].sat-1][0])>0.0) {
            double Coe=SQR(lam_L1/GNSS_LAMBDA_M[iGPS][1])*obsp->data[i].ionmap;
            dion=obsp->data[i].dion*Coe;
            ionvar=obsp->data[i].ionvar*Coe*Coe;
        }

        /* tropospheric corrections */
        if (resflag&&!spptrof->correction(&obsp->data[i],blh,0.7)) {
            obsp->data[i].errorMsg="troposphere correction error!";
            continue;
        }

        /* satellite receiver antenna delta */
        dsant=0,drant=0;
        if (resflag) {
            if (opt->posopt[0]>0) recantfunc.recantov(opt,staflag,&obsp->data[i],nav);
            if (opt->sateph==EPHOPT_PREC||opt->sateph==EPHOPT_SSRCOM)
                satantfunc.satantov(opt,Xpar,&obsp->data[i],nav);
            if (opt->sppiono==IONOOPT_IFLC) {
                drant=obsp->data[i].rant[NFREQ];
                dsant=obsp->data[i].sant[NFREQ];
            }
            else { 
                drant=obsp->data[i].rant[0];
                dsant=obsp->data[i].sant[0];
            }
        }
        
        /* pseudorange residual
         * Lnew = Pc - (dist + cdtr - cdts + dtro + dion) */
        Lnew= Pc - (obsp->data[i].dist + preClkOff[obsp->data[i].isys] - CLIGHT*obsp->data[i].dts[0] +
            obsp->data[i].dtro + dion - dsant - drant);
        /* add Lnew to L */
        Lobs.push_back(Lnew);

        /* coefficient matrix */
        for (int j=0; j<3; j++) Acoe.push_back(-obsp->data[i].sigvec[j]); //position
        for (int j=0; j<NSYS; j++) Abuff.push_back( j==obsp->data[i].isys? 1.0 : 0.0 ); //clock offset

        /* error variance
         * Rpc = code + sat + dcb + tro + ion */
        //double Rpc=varerr(obsp->data[i])+obsp->data[i].svar+obsp->data[i].dcbvar+obsp->data[i].trovar+ionvar;
        /* add Rpc to R vector */
        Rvec.push_back(varerr(obsp->data[i]));

        /* update obsp->data[i] */
        obsp->data[i].used=1; obsp->data[i].res[0]=Lobs.back(); obsp->data[i].ovar[0]=Rvec.back();

        /* update observation number */
        ncd[obsp->data[i].isys]++;
        numL++; ns++;
    }

    /* assign available clock offset parameters to Acoe and Xpar */
    for (int i=0; i<NSYS; i++) {
        if (ncd[i]>0) {
            for (int j=numL-1; j>=0; j--) {
                Acoe.insert( Acoe.begin()+(j+1)*numX, Abuff[j*NSYS+i] );
            }
            Xpar.push_back(preClkOff[i]);
            Xest[numX++]=3+i; //add clock offset parameter
        }
    }
    /* update Rx */
    Rx.assign(numX*numX,0);

    return numL;
}
/* compute Doppler observation\coefficient\covariance matrix for velocity --------- */
int gnss_single_c::doppler_matrix() {
    vnumL=0;

    /* reset observation\coefficient\covariance matrix */
    vAcoe.clear(); vLobs.clear(); vRvec.clear();
    for (int i=0; i<NSYS; i++) ndp[i]=0; //number of code observation of each system

    for (int i=0; i<obsp->n; i++) {
        obsp->data[i].vused=0;
        //now only single frequency is available
        lam=nav->lam[obsp->data[i].sat-1];
        if (obsp->data[i].used==0||obsp->data[i].D[0]==0.0||lam[0]==0.0||norm(obsp->data[i].posvel+3,3)<=0.0) {
            continue;
        }

        /* line-of-sight vector in ecef */
        geodist_GNSS(obsp->data[i].posvel,Xpar.begin(),obsp->data[i].sigvec);

        /* satellite velocity relative to receiver in ecef */
        double vel_rs[3];
        for (int j=0;j<3;j++) vel_rs[j]=obsp->data[i].posvel[j+3]-vXpar[j];

        /* range rate with earth rotation correction */
        double rate=dot(vel_rs,obsp->data[i].sigvec,3)+OMGE/CLIGHT*
            ( obsp->data[i].posvel[4]*Xpar[0]+obsp->data[i].posvel[1]*vXpar[0]-
              obsp->data[i].posvel[3]*Xpar[1]-obsp->data[i].posvel[0]*vXpar[1] );

        /* doppler residual */
        Lnew=-lam[0]*obsp->data[i].D[0]-(rate+vXpar[DPLNX-1]-CLIGHT*obsp->data[i].dts[1]);
        vLobs.push_back(Lnew);

        /* add Rpc to R vector */
        vRvec.push_back(SQR(opt->err[5]*lam[0]));

        /* design matrix */
        for (int j=0;j<DPLNX;j++) vAcoe.push_back(j<3 ? -obsp->data[i].sigvec[j] : 1.0);

        obsp->data[i].vused=1;
        vnumL++;
    }

    return vnumL;
}
/* compute observation\coefficient\covariance matrix with pre-position -------- */
int gnss_single_c::code_matrix_prePos() {
    numL=0;		    /* number of used observation (size of L) */
    numX=3;         /* number of estimated parameters */
    double blh[3];	/* geodetic position (lat,lon,h) */

    /* reset each matrix (except X=[NSYS+3,1]) */
    Acoe.clear(); Lobs.clear(); Rvec.clear();
    for (int i=0; i<NSYS; i++) {
        preClkOff[i]=0;
        ncd[i]=0; //number of code observation of each system
    }
    vector<double> Abuff;

    /* ecef (XYZ) to geodetic position */
    ecef2blh(Xpar.begin(),WGS84,blh);

    for (int i=ns=0; i<obsp->n; i++) {
        double Pc;		/* psudorange with code bias correction */
        double dion=0,ionvar=0;
        double lam_L1;	/* f1 lambda */

        /* reset observation data */
        obsp->data[i].used=0;

        /* excluded satellite? */
        if (obsp->data[i].exc>0) continue;
        if (exclude_sat(obsp->data[i])) {
            obsp->data[i].exc=1; continue;
        }

        if (!obsp->data[i].sys) {
            obsp->data[i].exc=1;
            obsp->data[i].errorMsg="unknown satellite!";
            continue;
        }

        /* reject duplicated observation data */
        if (i<obsp->n-1&&obsp->data[i].sat==obsp->data[i+1].sat) {
            obsp->data[i].errorMsg="duplicated observation data with next data!";
            i++;
            continue;
        }
        /* geometric distance/azimuth/elevation angle */
        if (
            (obsp->data[i].dist=geodist_GNSS(obsp->data[i].posvel,Xpar.begin(),obsp->data[i].sigvec))<300000.0
            || ( satazel(blh,obsp->data[i].sigvec,obsp->data[i].azel) < opt->min_elevation )
            ) {
            obsp->data[i].errorMsg="low elevation!";
            continue;
        }

        /* pseudorange with code bias correction */
        if ((Pc=prange(obsp->data[i],1))==0.0) {
            obsp->data[i].errorMsg="no avaliable pseudorange!";
            continue;
        }

        /* ionospheric corrections */
        if (!sppionf->correction(&obsp->data[i],nav,blh)) {
            obsp->data[i].errorMsg="ionosphere correction error!";
            continue;
        }

        /* GPS-L1 -> L1/B1 */
        if (opt->sppiono!=IONOOPT_IFLC&&(lam_L1=nav->lam[obsp->data[i].sat-1][0])>0.0) {
            double Coe=SQR(lam_L1/GNSS_LAMBDA_M[iGPS][1])*obsp->data[i].ionmap;
            dion=obsp->data[i].dion*Coe;
            ionvar=obsp->data[i].ionvar*Coe*Coe;
        }

        /* tropospheric corrections */
        if (!spptrof->correction(&obsp->data[i],blh,0.7)) {
            obsp->data[i].errorMsg="troposphere correction error!";
            continue;
        }

        /* satellite receiver antenna delta */
        dsant=0,drant=0;
        if (opt->posopt[0]>0) recantfunc.recantov(opt,staflag,&obsp->data[i],nav);
        if (opt->sateph==EPHOPT_PREC||opt->sateph==EPHOPT_SSRCOM)
            satantfunc.satantov(opt,Xpar,&obsp->data[i],nav);
        if (opt->sppiono==IONOOPT_IFLC) {
            drant=obsp->data[i].rant[NFREQ];
            dsant=obsp->data[i].sant[NFREQ];
        }
        else { 
            drant=obsp->data[i].rant[0];
            dsant=obsp->data[i].sant[0];
        }

        /* pseudorange residual
        * Lnew = Pc - (dist + cdtr - cdts + dtro + dion) */
        Lnew= Pc - (obsp->data[i].dist - CLIGHT*obsp->data[i].dts[0] +
            obsp->data[i].dtro + dion - dsant - drant);
        /* add Lnew to L */
        Lobs.push_back(Lnew);

        /* coefficient matrix */
        for (int j=0; j<3; j++) Acoe.push_back(-obsp->data[i].sigvec[j]); //position
        for (int j=0; j<NSYS; j++) Abuff.push_back( j==obsp->data[i].isys? 1.0 : 0.0 ); //clock offset
        preClkOff[obsp->data[i].isys]+=Lnew;

        /* error variance
        * Rpc = code + sat + dcb + tro + ion */
        //double Rpc=varerr(obsp->data[i])+obsp->data[i].svar+obsp->data[i].dcbvar+obsp->data[i].trovar+ionvar;
        /* add Rpc to R vector */
        Rvec.push_back(varerr(obsp->data[i]));

        /* update obsp->data[i] */
        obsp->data[i].used=1; obsp->data[i].res[0]=Lobs.back(); obsp->data[i].ovar[0]=Rvec.back();

        /* update observation number */
        ncd[obsp->data[i].isys]++;
        numL++; ns++;
    }
    /* constraint to avoid rank-deficient
    * if no this system */
    for (int i=0; i<NSYS; i++) {
        if (ncd[i]>0) {
            for (int j=numL-1; j>=0; j--) {
                Acoe.insert( Acoe.begin()+(j+1)*numX, Abuff[j*NSYS+i] );
            }
            Xest[numX++]=3+i; //add clock offset parameter
        }
    }

    /* get pre-estimated clock offset */
    if (preClkFlag) for (int i=0; i<NSYS; i++) {
        if (ncd[i]>0) {
            preClkOff[i]/=ncd[i];
        }
    }
    /* initialize Rx */
    Rx.assign(numX*numX,0);

    return numL;
}
/* compute observation\coefficient\covariance matrix with pre-velocity -------- */
int gnss_single_c::doppler_matrix_preVel() {
    vnumL=0;
    preClkDrf=0;

    /* reset observation\coefficient\covariance matrix */
    vAcoe.clear(); vLobs.clear(); vRvec.clear();
    for (int i=0; i<NSYS; i++) ndp[i]=0; //number of code observation of each system

    for (int i=0; i<obsp->n; i++) {
        obsp->data[i].vused=0;
        //now only single frequency is available
        lam=nav->lam[obsp->data[i].sat-1];
        if (obsp->data[i].used==0||obsp->data[i].D[0]==0.0||lam[0]==0.0||norm(obsp->data[i].posvel+3,3)<=0.0) {
            continue;
        }

        /* line-of-sight vector in ecef */
        geodist_GNSS(obsp->data[i].posvel,Xpar.begin(),obsp->data[i].sigvec);

        /* satellite velocity relative to receiver in ecef */
        double vel_rs[3];
        for (int j=0;j<3;j++) vel_rs[j]=obsp->data[i].posvel[j+3]-vXpar[j];

        /* range rate with earth rotation correction */
        double rate=dot(vel_rs,obsp->data[i].sigvec,3)+OMGE/CLIGHT*
            ( obsp->data[i].posvel[4]*Xpar[0]+obsp->data[i].posvel[1]*vXpar[0]-
                obsp->data[i].posvel[3]*Xpar[1]-obsp->data[i].posvel[0]*vXpar[1] );

        /* doppler residual */
        Lnew=-lam[0]*obsp->data[i].D[0]-(rate-CLIGHT*obsp->data[i].dts[1]);
        vLobs.push_back(Lnew);
        preClkDrf+=Lnew;

        /* add Rpc to R vector */
        vRvec.push_back(SQR(opt->err[5]*lam[0]));

        /* design matrix */
        for (int j=0;j<DPLNX;j++) vAcoe.push_back(j<3 ? -obsp->data[i].sigvec[j] : 1.0);

        obsp->data[i].vused=1;
        vnumL++;
    }
    /* get pre-estimated clock drift */
    if ( preClkFlag && vnumL>0 ) {
        preClkDrf/=vnumL;
    }

    return vnumL;
}
/* gnss modeled correction (PCV+ClockJump+PhaseWindUp) ---------------------------- */
void gnss_single_c::model_correction() {
    obsp->used=0;

    for (int i=0; i<obsp->n; i++) {
        gnss_obsd_c *op=&obsp->data[i];
        if (!op->used) continue;

        /* update satellite and receiever antenna phase center correction */
        if (opt->posopt[0]>0) recantfunc.recantov(opt,staflag,op,nav);
        if (opt->sateph==EPHOPT_PREC||opt->sateph==EPHOPT_SSRCOM) {
            satantfunc.satantov(opt,Xpar,op,nav);
        }

        /* calculate phase windup correction */
        if (opt->posopt[1]) {
            if (!satantfunc.phase_windup(op->time,Xpar,op,ssat[op->sat-1].phw[obsp->rcv])) {
                op->used=0; op->exc=1;
                op->errorMsg="phase windup errior!";
                continue;
            }
        }
        obsp->used++;
    }
}
/* update gnss Clock offset ------------------------------------------------------- */
void gnss_single_c::update_clock_offset() {
    for (int i=0; i<obsp->n; i++) {
        gnss_obsd_c *op=&obsp->data[i];
        if (!op->used) continue;
        /* time system and receiver bias offset correction */
        op->dtr=preClkOff[op->isys]/CLIGHT;
    }
}
/* update gnss modeled correction (Clock+PCV+ClockJump+PhaseWindUP) --------------- */
void gnss_single_c::update_model_correction() {
    obsp->used=0;

    for (int i=0; i<obsp->n; i++) {
        gnss_obsd_c *op=&obsp->data[i];
        if (!op->used) continue;
        /* time system and receiver bias offset correction */
        op->dtr=preClkOff[op->isys]/CLIGHT;

        /* update satellite and receiever antenna phase center correction */
        if (opt->posopt[0]>0) recantfunc.recantov(opt,staflag,op,nav);
        if (opt->sateph==EPHOPT_PREC||opt->sateph==EPHOPT_SSRCOM) {
            satantfunc.satantov(opt,Xpar,op,nav);
        }

        /* calculate phase windup correction */
        if (opt->posopt[1]) {
            gtime_c sss=op->time;
            sss.timeadd(-op->dtr);
            if (!satantfunc.phase_windup(sss,Xpar,op,ssat[op->sat-1].phw[obsp->rcv])) {
                op->used=0; op->exc=1;
                op->errorMsg="phase windup errior!";
                continue;
            }
        }

        obsp->used++;
    }
}
/* check position variance -------------------------------------------------------- */
int gnss_single_c::check_pos_var() {
    if (Rx.size()<9) {
        posVar=9999999;
        obsp->errorMsg="parameter number < 3!";
        return 0;
    }
    else {
        posVar = Rx[0] + Rx[1+numX] + Rx[2+2*numX];
        if (posVar>MAX_POS_VAR || posVar<0) {
            string strVar,strNs,strSgm;
            doul2str(0,3," ",SQRT(posVar),strVar);
            int2str(0," ",ns,strNs);
            doul2str(0,3," ",SQRT(adjfunc->sigma0_2),strSgm);
            obsp->errorMsg="SPP 3D position variance is to large: " + strVar +
                           "m! nSat=" + strNs + ", LSQ Sigma_0=" + strSgm;
            return 0;
        }
    }
    return 1;
}
/* calculate pdop value of each satellite system ---------------------------------- */
void gnss_single_c::get_pdop() {
    for (int i=0; i<NSYS+1; i++) Pdop[i]=0.0;
    double pdop_val=0.0;
    if (obsp==obsr&&numL>numX&&Acoe.size()==numL*numX) {
        vector<double> Asys[NSYS+1];
        QAA.assign(9,0.0); Asys[NSYS].clear();
        int usednum=0;
        /* pdop for each single system */
        for (int sys=0; sys<NSYS; sys++) {
            //update Acoe for all systems
            for (int i=0; i<ncd[sys]; i++) for (int j=0; j<3; j++)
                Asys[NSYS].push_back( Acoe[j+(i+usednum)*numX] );
            /* pdop of this system */
            if (ncd[sys]>=3) {
                //Asys of Acoe to calculate pdop
                Asys[sys].assign(ncd[sys]*3,0.0);
                for (int i=0; i<ncd[sys]; i++) for (int j=0; j<3; j++)
                    Asys[sys][j+i*3] = Acoe[j+(i+usednum)*numX];
                //calculate QAA=(Asys'*Asys)^-1
                matmul_vec("TN",3,3,ncd[sys],1.0,Asys[sys],Asys[sys],0.0,QAA);
                if (matinv(QAA,3)==-1) QAA[0]=QAA[4]=QAA[8]=333.3;
                //calculate pdop 
                double pdop_val=QAA[0]+QAA[4]+QAA[8];
                Pdop[sys] = pdop_val>999.0 ? 999.0 : SQRT(pdop_val);
            }
            else if ( ncd[sys]>0 ) Pdop[sys]=999.0;
            usednum+=ncd[sys];
        }
        /* pdop of all system */
        if (usednum>3) {
            //calculate QAA=(Acoe'*Acoe)^-1
            matmul_vec("TN",3,3,usednum,1.0,Asys[NSYS],Asys[NSYS],0.0,QAA);
            if (matinv(QAA,3)==-1) QAA[0]=QAA[4]=QAA[8]=333.3;
            //calculate pdop 
            double pdop_val=QAA[0]+QAA[4]+QAA[8];
            Pdop[NSYS] = pdop_val>999.0 ? 999.0 : SQRT(pdop_val);
        }
    }
}
/* write state information to log_stream ------------------------------------------ */
void gnss_single_c::write_state() {
    if (!log_stream.is_open()) return;

    /* error message */
    if (stat==SOLQ_NONE) {
        /* time */
        obsp->data[0].time.time2str(3);
        log_stream << obsp->data[0].time.sep;
        log_stream << "   " << obsp->errorMsg << "\n";
        for (int i=0; i<obsp->n; i++) {
            if ( obsp->data[i].errorMsg.size()>0 && obsp->data[i].errorMsg.size()<MAXSTRMSG ) 
                log_stream << "    " << ssat[obsp->data[i].sat-1].id << ": "
                << obsp->data[i].errorMsg << "\n";
        }
    }
    else if ( opt->mode==PMODE_SINGLE && obsp==obsr ) {
        /* state message */
        if ( stat==SOLQ_SINGLE || stat==SOLQ_SBAS ) {
            /* time */
            obsp->data[0].time.time2str(3);
            log_stream << obsp->data[0].time.sep;

            /* set float format */
            log_stream.setf(ios::fixed);
            log_stream << setprecision(3);
            /* adjustment time (in the same line of time) */
            log_stream << "   Filter Time:" << setw(5) << dtime;
            /* solution number (in the same line of time) */
            log_stream << "   Solution Number:" << setw(7) << nsol;
            /* standard variance of unit weight in LSQ */
            log_stream << "   LSQ Sigma0: " << setw(9) << SQRT(adjfunc->sigma0_2);
            log_stream << "\n";
            log_stream << setprecision(4);
            /* receiver clock */
            log_stream << "    Clock Error: ";
            for (int i=0; i<NSYS; i++) {
                if (ncd[i]<1) {
                    log_stream << "   " << SYS_LIST[i] << ":" << setw(13) << 0.0 //new Clock
                        << setw(8) << 0.0; //new Clock STD
                }
                else {
                    log_stream << "   " << SYS_LIST[i] << ":" << setw(13) << preClkOff[i] //new Clock
                           << setw(8) << SQRT(solp->vclk[i+i*NSYS]); //new Clock STD
                }
            }
            log_stream << "\n";
            /* satellites number of different system */
            if (opt->logmsg&GNSSLOG_NSAT) {
                int allsat=0;
                log_stream << "    System nSat: ";
                for (int sys=0; sys<NSYS; sys++) {
                    log_stream << "   " << SYS_LIST[sys] << ":" << setw(8) << ncd[sys];
                    allsat+=ncd[sys];
                }
                log_stream << "   All:" << setw(8) << allsat;
                log_stream << "\n";
            }
            /* BODY of state ------------------------- */
            /* PDOP value of different satellite systems */
            if (opt->logmsg&GNSSLOG_PDOP) {
                /* calculate pdop of different satellite system */
                gnss_single_c::get_pdop();
                log_stream << "    PDOP Values: ";
                for (int sys=0; sys<NSYS; sys++)
                    log_stream << "   " << SYS_LIST[sys] << ":" << setw(8) << Pdop[sys];
                if (numL>3) log_stream << "   All:" << setw(8) << Pdop[NSYS];
                log_stream << "\n";
            }
            /* xpar */
            if (opt->logmsg&GNSSLOG_XPAR) {
                log_stream << "    X-parameter: ";
                for (int i=0; i<numX; i++) log_stream << "   " << Xpar[i];
                log_stream << "\n";
            }
            /* satellite parameters (prn, az, el) */
            log_stream << setprecision(4);
            for (int i=0; i<obsp->n; i++) {
                if (!obsp->data[i].used) continue;
                gnss_ssat_c *sss=&ssat[obsp->data[i].sat-1];
                lam=nav->lam[obsp->data[i].sat-1];
                log_stream << "    " << sss->id << ":"; //prn
                log_stream << " A" << setw(8) << obsp->data[i].azel[0]*R2D << "  E" <<
                    setw(7) << obsp->data[i].azel[1]*R2D;//az el

                /* x y z dt and transmisson time*/
                log_stream << "  XYZT:";
                log_stream << setw(15) << obsp->data[i].posvel[0];
                log_stream << setw(15) << obsp->data[i].posvel[1];
                log_stream << setw(15) << obsp->data[i].posvel[2];
                log_stream << setw(11) << obsp->data[i].dts[0]*1E6;
                log_stream << setw(8) << obsp->data[i].time.timediff(obsp->data[i].sigtime);
                /* troposphere */
                log_stream << "  TRO:" << setw(9) << (opt->spptrop==0? 0.0 : obsp->data[i].dtro);
                /* ionosphere delay of L1 (m) */
                log_stream << "  ION:" << setw(10) << (opt->sppiono==0||opt->sppiono==IONOOPT_IFLC? 0.0 :
                    SQR(lam[0]/GNSS_LAMBDA_M[iGPS][1])*obsp->data[i].dion*obsp->data[i].ionmap);
                log_stream;
                /* residual */
                log_stream << "  RES:" << setw(8) << obsp->data[i].res[0];
                /* measurement code */
                if (opt->sppiono==IONOOPT_IFLC) {
                    log_stream << " MEAS: PC-" <<  obsp->data[i].PCode[0] << "-" <<  obsp->data[i].PCode[1];
                }
                else {
                    log_stream << " MEAS: P-" <<  obsp->data[i].PCode[0];
                }
                /* warning message */
                log_stream << "  " << obsp->data[i].warningMsg << "\n";
            }
            /* output excluded satellite messages */
            if (opt->logmsg&GNSSLOG_EXCSAT) {
                for (int i=0; i<obsp->n; i++) {
                    if ( !obsp->data[i].used && obsp->data[i].errorMsg.size()>0 && obsp->data[i].errorMsg.size()<MAXSTRMSG ) 
                        log_stream << "    " << ssat[obsp->data[i].sat-1].id << ": " << obsp->data[i].errorMsg << "\n";
                }
            }
        }
    }
}
/* update satellite sate vector (ssat) -------------------------------------------- */
void gnss_single_c::update_ssat() {
    if (obsp==obsr) {
        for (int i=0; i<MAXSAT; i++) {
            ssat[i].vs=0;
        }
        for (int i=0; i<obsp->n; i++) {
            for (int j=0; j<NFREQ; j++) ssat[obsp->data[i].sat-1].resL[j]=999;
            ssat[obsp->data[i].sat-1].azel[0]=obsp->data[i].azel[0];
            ssat[obsp->data[i].sat-1].azel[1]=obsp->data[i].azel[1];
            ssat[obsp->data[i].sat-1].snr[0]=obsp->data[i].SNR[0];
            if (!obsp->data[i].used) continue;
            ssat[obsp->data[i].sat-1].vs=1;
            ssat[obsp->data[i].sat-1].resP[0]=obsp->data[i].res[0];
        }
    }
}
/* update solution vector (sol) --------------------------------------------------- */
void gnss_single_c::update_sol() {

    solp->type=0;
    solp->nsol=nsol;
    solp->obsTime=solp->time=obsp->data[0].time;
    solp->time.timeadd(-Xpar[3]/CLIGHT);
    if (obsp==obsr) soltime[0]=solp->time;
    solp->NL=numL;
    /* clock */
    solp->xclk[iGPS]=preClkOff[iGPS]; /* receiver clock bias (s) */
    solp->xclk[iGLO]=preClkOff[iGLO]; /* glo-gps time offset (s) */
    solp->xclk[iGAL]=preClkOff[iGAL]; /* gal-gps time offset (s) */
    solp->xclk[iBDS]=preClkOff[iBDS]; /* bds-gps time offset (s) */
    for (int j=0,k=0; j<NSYS; j++) {
        if (ncd[j]>0) {
            solp->vclk[j+j*NSYS]=Rx[ (k+3) + (k+3)*numX ];
            k++;
        }
    }
    /* position */
    for (int j=0; j<3; j++) solp->xpos[j]=Xpar[j];
    for (int j=0; j<3; j++) for (int k=0; k<3; k++)
        solp->vpos[k+j*3]=Rx[k+j*numX];
    solp->ns=(unsigned int)ns;
    solp->age=solp->ratio=0.0;

    solp->stat=opt->sateph==EPHOPT_SBAS ? SOLQ_SBAS : SOLQ_SINGLE;
    /* update X_ALL and Rx_ALL */
    if (opt->mode==PMODE_SINGLE) {
        if (obsp==obsr) for (int i=0; i<3; i++) {
            X_ALL[i]=Xpar[i];
            Rx_ALL[i+i*N_ALL]=Rx[i+i*numX];
            //symmetrize covariance
            for (int j=0; j<i; j++) Rx_ALL[j+i*N_ALL] = Rx_ALL[i+j*N_ALL] = Rx[j+i*numX];
        }
    }
}
/* growth variance of static rover position --------------------------------------- */
double gnss_single_c::growth_rate_static() {
    return tt<3?SQR(0.01*tt):SQR(0.1);
}
/* initialize covariance matrix of amb parameters */
void gnss_single_c::init_RxAMB(unsigned int sat,int freq) {
    X_ALL[NX+(sat-1)*numF+freq]=0;
    for (int i=0; i<N_ALL; i++)
        Rx_ALL[ i + (NX+(sat-1)*numF+freq)*N_ALL ] = Rx_ALL[ NX+(sat-1)*numF+freq + i*N_ALL ] = 0.0;
}
/* initialize covariance matrix of all parameters */
void gnss_single_c::init_RxALL(int ipar) {
    X_ALL[ipar]=0;
    for (int i=0; i<N_ALL; i++) Rx_ALL[ i + ipar*N_ALL ]=Rx_ALL[ ipar + i*N_ALL ]=0.0;
}
/* estimate receiver position ----------------------------------------------------- */
int gnss_single_c::singlepos() {

    /* SPPNX: x,y,z dtr(GPS,GLO,GLA,CMP) */
    Xpar.assign(3,0);
    for (int i=0; i<3; i++) {
        Xpar[i]=solp->xpos[i];
    }
    for (int j=0; j<NSYS; j++) preClkOff[j] = 0 ;

    resflag=ratflag=0; 

    int nnn=0;
    lsqflag=0;
    /* compute receiver position and clock bias */
    for (niter=0; niter<MAXITR; niter++) {
        nnn++;

        /* initialize positioning if large residual ratio detected */
        if (ratflag==1) {
            Xpar.assign(solp->xpos.begin(),solp->xpos.end());
            niter=0;
        }
        else if (Xpar.size()>3) {
            Xpar.erase(Xpar.begin()+3,Xpar.end());
        }

        /* compute observation\coefficient\covariance matrix */
        if (code_matrix()<numX) {
            /* verify if lack of valid observation */
            obsp->errorMsg="lack of valid satellites for SPP!";
            break;
        }

        dtime=0; unsigned int stime=tickget();
        /* adjustment with least square function */
        if ((lsqflag=adjfunc->LSQ_Rvec(Acoe,Lobs,Rvec,Xpar,Rx,numL,numX))==-1) {
            obsp->errorMsg="least square error!";
            break;
        }
        dtime=(int)(tickget()-stime);

        /* update preClkOff using estimated Xpar */
        int iClk=0;
        for (int j=0; j<NSYS; j++) preClkOff[j] = ( ncd[j]>0? Xpar[3+iClk++] : 0 );

        if (lsqflag<100&&niter<MAXITR-1) {
            /* exclude constellations with only a few available satellites */
            exc_few_sat();
            /* exclude observations with large residual */
            exc_largeres(0);
        }
        else {
            /* exclude observation with large residual */
            if (resflag==0) {
                /* exclude constellations with only a few available satellites */
                exc_few_sat();
                exc_largeres(1);
            }
            /* out solution if already exlclude large residual */
            else if (resflag>0) {
                if (check_pos_var()<1) return SOLQ_NONE;

                obsp->used=ns;
                /* update each gnss_obs_c with corrections (Clock+PCO&PCV+ClockJump+PhaseWindUP) */
                update_model_correction();
                /* calculate velocity (only for rover) */
                if ( opt->dynamics && obsp==obsr ) gnss_vel();
                return SOLQ_SINGLE;
            }
            /* restart estimate without large residual observation */
            resflag++; niter=0;
        }
    }

    Lobs.clear(); Acoe.clear();
    Rvec.clear(); Xpar.clear(); Rx.clear();

    return SOLQ_NONE;
}

/* single point position -------------------------------------------------------------
* recnum	: receiver number (1:rover,2:base)
* --------------------------------------------------------------------------------- */
int gnss_single_c::single() {
    warningMsg="";

    /* update consecutive number for each satellite */
    update_sat_nCon();

    /* initialize single position */
    gnss_ini_EPH_DCB();

    /* detect wrong observations with double frequency */
    detect_wrong_obs_2freq();

    /* observation and ambiguity information for each satellites */
    if ( opt->mode>=PMODE_PPP_KINEMA && obsp==obsr ) {
        update_obs_amb_ssat();
    }

    /* estimate receiver position */
    stat=singlepos();

    if (stat) {
        if (opt->mode==PMODE_SINGLE) nsol++; ntotal++;
        /* update satellite state */
        gnss_single_c::update_ssat();
        /* update solution */
        gnss_single_c::update_sol();
    } 
    /* write state information to log_stream */
    gnss_single_c::write_state();

    return stat;
}
/* update solution vectors -----------------------------------------------------------
*  argv :  int      rovbas    I   receiver index (1:rover, 2:base)
* --------------------------------------------------------------------------------- */
void gnss_single_c::init_sol(int rovbas) {
    /* rover station */
    if (rovbas==1) {
        if (sol.back().stat!=SOLQ_NONE) {
            sol.erase(sol.begin()); sol.push_back(gnss_sol_c(this));
        }
        solp=&sol.back();
        obsp=obsr;
        obsp->errorMsg="";
        staflag=0;
        ncd=ncode[0];
        ndp=ndoppler[0];
    }
    /* base station */
    else if (rovbas==2) {
        if (b_sol.back().stat!=SOLQ_NONE) {
            b_sol.erase(b_sol.begin()); b_sol.push_back(gnss_sol_c(this));
        }
        solp=&b_sol.back();
        obsp=obsb;
        obsp->errorMsg="";
        staflag=1;
        ncd=ncode[1];
        ndp=ndoppler[1];
    }
    Xini.assign(SPPNX,1);
    Xest.assign(SPPNX,-1);
    Xest[0]=0; Xest[1]=1; Xest[2]=2;
}
/* preprocess gnss observation ------------------------------------------------ */
int gnss_single_c::preprocess_gnss_obs(const double pos[3], const double vel[3]) {
    warningMsg="";

    init_sol(1);

    /* update consecutive number for each satellite */
    update_sat_nCon();

    /* initialize single position */
    gnss_ini_EPH_DCB();

    /* detect wrong observations with double frequency */
    detect_wrong_obs_2freq();

    /* intiialize position and velocity */
    gnss_get_posvel(pos,vel,NULL,NULL);
    /* assign to Xpar */
    Xpar.assign(3,0);
    vXpar.assign(DPLNX,0);
    for (int i=0; i<3; i++) {
        Xpar[i]=ext_posvel[i];
        vXpar[i]=ext_posvel[i+3];
    }

    /* compute observation\coefficient\covariance matrix */
    if (code_matrix_prePos()<numX) {
        /* verify if lack of valid observation */
        obsp->warnMsg="lack of valid satellites for SPP!";
    }
    /* exculde observations with large residual */
    exc_largeres(0);
    /* exclude constellations with only a few available satellites */
    exc_few_sat();

    return 1;
}
/* get gnss observation equation (SPP) -------------------------------------------- */
int gnss_single_c::gnss_obs_equation() {

    /* position observation */
    /* compute observation\coefficient\covariance matrix */
    if (code_matrix_prePos()<numX) {
        /* verify if lack of valid observation */
        obsp->warnMsg="lack of valid satellites for SPP!";
    }
    /* update Rvar using Rvec */
    Rvar.assign(numL*numL,0.0);
    for (int j=0; j<numL; j++) Rvar[j*numL+j]=Rvec[j];

    /* velocity observation */
    /* velocity observation equation */
    if ( opt->dynamics ) {
        if (doppler_matrix_preVel()<DPLNX) {
            /* verify if lack of valid observation */
            obsp->warnMsg+="lack of valid satellites for velocity!";
        }
        /* update vRvar using vRvec */
        vRvar.assign(vnumL*vnumL,0.0);
        for (int j=0; j<vnumL; j++) vRvar[j*vnumL+j]=vRvec[j];
        /* initialize vRx */
        vRx.assign(DPLNX*DPLNX,0);
    }

    /* correct clock using pre-estimated clock offset and drift */
    if (preClkFlag) {
        process_pre_clock();

        /* update Xpar and Rx for clock offset */
        for (int i=0; i<NSYS; i++) if (ncd[i]>0) Xpar.push_back(preClkOff[i]); //clk offset
        for (int i=3; i<numX; i++) Rx[ i*numX + i ] = VAR_CLK_OFF1; //clk offset noise

        /* update vXpar and vRx for clock drift */
        if ( opt->dynamics ) {
            vRx[DPLNX*DPLNX-1]=VAR_CLK_DRF1;
            vXpar[DPLNX-1]=preClkDrf;
        }
    }
    else {
        /* update Rx */
        for (int i=3; i<numX; i++) {
            Xpar.push_back(0);
            Rx[ i*numX + i ] = VAR_CLK_OFF2; //clk offset noise
        }
        /* update vRx */
        if ( opt->dynamics ) {
            vRx[15]=VAR_CLK_DRF2; //clk offset noise
        }
    }

    return 1;
}
/* update gnss status (SPP) ------------------------------------------------------- */
int gnss_single_c::gnss_update_status() {

    obsp->used=ns;

    nsol++; ntotal++;

    /* update satellite sate vector (ssat) */
    gnss_single_c::update_ssat();
    /* update solution vector (sol) */
    gnss_single_c::update_sol();

    stat=SOLQ_FLOAT;
    /* write state information to log_stream */
    gnss_single_c::write_state();

    return stat;
}
/* base station single-position --------------------------------------------------- */
int gnss_single_c::basepos() {
    init_sol(2);

    /* intialize base position */
    if (opt->baspos!=POSOPT_SINGLE&&opt->mode!=PMODE_MOVEB)
        for (int i=0; i<3; i++) solp->xpos[i]=opt->rb[i];

    return single();
}
/* gnss velocity function ----------------------------------------------------- */
int gnss_single_c::gnss_vel() {

    vstat=VSOLQ_NONE;
    vXpar.assign(DPLNX,0);
    vRx.assign(DPLNX*DPLNX,0);
    lsqflag=0;
    for (int i=0; i<MAXITR; i++) {

        if (doppler_matrix()<DPLNX) {
            /* verify if lack of valid observation */
            obsp->warnMsg="lack of valid satellites for velocity!";
            return 0;
        }

        /* update vRvar using vRvec */
        vRvar.assign(vnumL*vnumL,0.0);
        for (int j=0; j<vnumL; j++) vRvar[j*vnumL+j]=vRvec[j];

        /* adjustment with least square function */
        if ((lsqflag=adjfunc->LSQ(vAcoe,vLobs,vRvar,vXpar,vRx,vnumL,DPLNX))==-1) {
            obsp->warnMsg="least square error for velocity!";
            break;
        }

        if (lsqflag>=100||niter>=MAXITR-1) {
            for (int j=0; j<3; j++) {
                solp->xvel[j]=vXpar[j];
                for (int k=0; k<3; k++) solp->vvel[j*3+k]=vRx[j*DPLNX+k];
            }
            break;
        }
    }
    vstat=VSOLQ_UD;

    return vstat;
}
/* gnss single point position (SPP) function -------------------------------------- */
int gnss_single_c::gnss_pos() {
    init_sol(1);

    return single();
}
