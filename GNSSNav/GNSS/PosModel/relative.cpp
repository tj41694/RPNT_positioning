#include "GNSS/PosModel/position.h"
#include "BaseFunction/basefunction.h"

/* constant --------------------------------------------------------------------------------------- */
#define MINSAT          4                   /* minimum satellite number */
#define SQR(x)          ((x)*(x))
#define SQRT(x)         ((x)<0.0?-sqrt(-x):sqrt(x))

/* inter-system bias */
#define INIT_GRA	        1E-5            /* initial tro gradient */
#define VAR_GRA             2.5E-5          /* initial variance of gradient (m^2) */
#define INIT_HWBIAS	        0               /* initial glo frequency-dependent amb bias */
#define VAR_HWBIAS          0.01            /* initial variance of h/w bias ((m/MHz)^2) */
#define RAT_HWBIAS	        1E-7            /* growth rate of std h/w bias (m/MHz/sqrt(s)) */

/* ambiguity bias */
#define FIXED_AMB_VAR       1E-4            /* fixed ambiguity variance */
#define FIXED_AMB_ADJ       1E-4            /* adjusted variance for fixed ambiguity when have new ambiguity */

/* navigation parameter */
#define VAR_POS             3600.0          /* initial variance of receiver pos (m^2) */
#define VAR_VEL             3600.0			/* initial variance of receiver vel ((m/s)^2) */
#define MAX_ATT_ERR         5               /* maximum attitude error (deg) */

#define VAR_CLK_OFF1        3600.0          /* initial variance of pre-estimated receiver clock offset (m^2) */
#define VAR_CLK_DRF1        3600.0          /* initial variance of pre-estimated receiver clock drift (m^2/s^2) */
#define VAR_CLK_OFF2        9E12	        /* initial variance of receiver clock offset (m^2) */
#define VAR_CLK_DRF2        1E6             /* initial variance of receiver clock drift (m^2/s^2) */

#define TTOL_MOVEB          (1.0+2*DTTOL)

const string SYS_LIST="GREC";			    /* system flag GPS, GLONASS, Galileo, BeiDou */

/* tropospheric parameters */
const string TRO_PAR_Z[]={ "ZTD1:","ZTD2:" };
const string TRO_PAR_ZNE[]={ "ZTD1:","G_N1:","G_E1:","ZTD2:","G_N2:","G_E2:" };
vector<string> TRO_PAR;

/* relative position class ---------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
gnss_relative_c::gnss_relative_c() {
    /* station geographic parameters (ini. in  updatepar()) */
    odt = 0.0;
    for (int i=0; i<3; i++)
        Rxyz[i] = Bxyz[i] = Rblh[i] = Bblh[i] = Rtide[i] = Btide[i] = 0.0;

    baseline = Rbl = 0.0;
    /* reference satellite parameters */
    satnum = numsys = 0;
    for (int i=0; i<NSYS; i++) {
        for (int j=0; j<NFREQ; j++) { L_nsat[i][j] = 0; ambRB1[i][0]=0; }
        irfsat[i] = -1;
        rfsat[0][i] = rfsat[1][i] = 0;
        disRB1[i] = troRB1[i] = ionRB1[i] = 0.0;
    }
    gloRB1[0] = gloRB1[1] = 0.0;
    disRB2 = troRB2 = ionRB2 = gloRB2 = ambRB2 = 0.0;

    rp = bp = NULL;
    for (int sys=0; sys<NSYS; sys++) {
        for (int i=0; i<NFREQ*2; i++)
            nobs[sys][i] = varRB1[sys][i] = 0;
    }

    for (int i=0; i<NFREQ+1; i++) nreset[i]=reset_flag[i]=0;

    iT = iG = iI = iA = 0;
    nA = nI = 0;
    ifact1=ifact2=0.0;
    iniamb=60;
    fix_num = n_Damb = En_Damb = 0;
}
gnss_relative_c::~gnss_relative_c() {
    rp=NULL; bp=NULL;
    lam=lam1=lam2=NULL;

    usesat.clear(); rovsat.clear(); bassat.clear();

    Acoep=NULL;
    for (int sys = 0; sys<NSYS; sys++) Airfsat[sys].clear();
    Asatnum.clear();
    Ais.clear();

    Rxvec.clear();

    ambnum.clear();
    fix_flag.clear();
    fix_amb.clear(); sum_var.clear();
    fix_Xpar.clear(); fix_Rx.clear();
}
/* Implementation functions ------------------------------------------------------- */
/* correct and save pre-etimated clock offset and drift --------------------------- */
void gnss_relative_c::process_pre_clock() {
    /* clock offset */
    /* update solution with clock offset (for PPP) */
    solp->xclk[iGPS]=preClkOff[iGPS]; /* receiver clock bias (s) */
    solp->xclk[iGLO]=preClkOff[iGLO]; /* glo-gps time offset (s) */
    solp->xclk[iGAL]=preClkOff[iGAL]; /* gal-gps time offset (s) */
    solp->xclk[iBDS]=preClkOff[iBDS]; /* bds-gps time offset (s) */

    /* clock drift */
    preClkVar[1]=VAR_CLK_DRF1;
    if ( opt->dynamics ) {
        /* update velocity observations with clock drift */
        for (int i=0; i<vnumL; i++) vLobs[i]-=preClkDrf;
        /* update vXpar and vRx for clock drift */
        vXpar[DPLNX-1]=preClkDrf;
        vRx[DPLNX*DPLNX-1]=preClkVar[1];
    }
}
/* test satellite system ---------------------------------------------------------- */
int gnss_relative_c::test_system(int Sysyem,int System_Num) {
    switch (Sysyem) {
        case SYS_GPS: return System_Num==0;
        case SYS_QZS: return System_Num==0;
        case SYS_SBS: return System_Num==0;
        case SYS_GLO: return System_Num==1;
        case SYS_GAL: return System_Num==2;
        case SYS_BDS: return System_Num==3;
    }
    return 0;
}
/* select common satellites of rover and base ----------------------------------------
 * and get highest elevation satellite of rover ----------------------------------- */
int gnss_relative_c::selcomsat_phase() {
    usesat.clear(); rovsat.clear(); bassat.clear();
    for (int sys=0; sys<NSYS; sys++) for (int freq=0; freq<NFREQ; freq++) {
        L_nsat[sys][freq] = 0;
    }

    /* select common satellites */
    for (int i=0,j=0; i<obsr->n&&j<obsb->n; i++,j++) {
        if (obsr->data[i].sat<obsb->data[j].sat) j--;
        else if (obsr->data[i].sat > obsb->data[j].sat) i--;
        else {
            if ( obsr->data[i].used&&obsb->data[j].used &&
                ncode[0][obsr->data[i].isys]>=MINSYSOBS_RTK &&
                ncode[1][obsr->data[i].isys]>=MINSYSOBS_RTK ) {
                usesat.push_back(obsr->data[i].sat);
                rovsat.push_back(i);
                bassat.push_back(j);
            }
        }
    }
    satnum = usesat.size();

    if (satnum<=0) return 0;

    //initialize reference satellite parameters
    init_refsat();
    /* select reference satellite according to phase observation number and elevation */
    if (opt->mode>PMODE_DGPS) {
        /* get phase observation number of common satellite */
        vector<int> nphase_sat(satnum,0);
        for (int sys=0; sys<NSYS; sys++) { // loop NSYS systems
            irfsat[sys]=-1;
            for (int i=0; i<satnum; i++) { //loop every satellite
                int L_ava[NFREQ] ={ 0 };
                rp = &obsr->data[rovsat[i]]; bp = &obsb->data[bassat[i]];
                if (rp->isys!=sys) continue; //test satellite system
                for (int freq = 0; freq<NF; freq++) { //loop frequency
                    if (rp->L[freq]!=0.0&&rp->P[freq]!=0.0&&bp->L[freq]!=0.0&&bp->P[freq]!=0.0) {
                        nphase_sat[i]++;
                        L_ava[freq] = 1; //available phase obs. (L1,L2,...)
                    }
                }
                if (opt->ionoopt==IONOOPT_IFLC && L_ava[0] && L_ava[1]) L_nsat[sys][0]++;
                if (opt->ionoopt!=IONOOPT_IFLC) {
                    for (int freq = 0; freq<NF; freq++) if (L_ava[freq]) L_nsat[sys][freq]++;
                }
                /* get satellites of rover with the most phase observation and highest elevation  */
                if ( irfsat[sys]<0 || nphase_sat[i]>nphase_sat[irfsat[sys]] ||
                    ( nphase_sat[i]==nphase_sat[irfsat[sys]] &&
                        rp->azel[1]>obsr->data[rovsat[irfsat[sys]]].azel[1] ) ) {
                    irfsat[sys] = i;
                }
            }

            if (irfsat[sys]<0) continue; //continue if no satellite of this system

            rfsat[1][sys] = usesat[irfsat[sys]];
            numsys++;
        }
    }
    else for (int sys=0; sys<NSYS; sys++) {
        for (int i=0; i<satnum; i++) {
            if (obsr->data[rovsat[i]].isys!=sys) continue; //test satellite system
            if (irfsat[sys]<0 || obsr->data[rovsat[i]].azel[1]>obsr->data[rovsat[irfsat[sys]]].azel[1])
                irfsat[sys] = i;
        }
        if (irfsat[sys]>=0) { rfsat[1][sys] = usesat[irfsat[sys]]; numsys++; }
    }

    return satnum;
}
/* initialize reference satellites parameters ------------------------------------- */
void gnss_relative_c::init_refsat() {
    numsys = 0;
    for (int i=0; i<NSYS; i++) {
        irfsat[i] = -1;
        rfsat[1][i] = 0;
        disRB1[i] = troRB1[i] = ionRB1[i] = 0.0;
        for (int j=0; j<NFREQ; j++) {
            ambRB1[i][0]=0;
            varRB1[i][j]=varRB1[i][NFREQ+j]=0;
        }
    }
    gloRB1[0] = gloRB1[1] = 0.0;
}
/* initialize vector and matrix according common satellite ------------------------ */
void gnss_relative_c::init_arrmat(int prePos) {
    Lobs.clear(); Acoe.clear();
    Rvec.clear(); Rvar.clear(); Rx.clear();
    QAA.clear();
    fix_flag.clear();
    fix_amb.clear(); sum_var.clear();
    fixIndex.clear();
    fix_Xpar.clear(); fix_Rx.clear();
    ambnum.assign(numF*satnum,-1);
    ambElevation[0].clear(); ambElevation[1].clear();

    /* number of kinds of parameters */
    /* NP NT NG nI nA */
    if (opt->ionoopt==IONOOPT_CONST) nI=satnum;
    nA = 0;
    numX = NP+NT+NG+nI+nA;
    nTrp = NT/2;
    if (nTrp==3) TRO_PAR.assign(TRO_PAR_ZNE,TRO_PAR_ZNE+6);
    else if (nTrp==1) TRO_PAR.assign(TRO_PAR_Z,TRO_PAR_Z+2);

    /* set start index of kinds of parameters */
    iT=NP; iG=iT+NT; iI=iG+NG; iA=iI+nI;

    /* initialize parameter and its variance vector (no amb) */
    if (prePos==1) Xpar.insert(Xpar.end(),(numX-3),0.0); // if have pre-estimated position
    else Xpar.assign(numX,0.0);
    Rxvec.assign(numX,0.0);

    /* intialize parameter information vector */
    Xini.assign(numX,0); Xest.assign(numX,-1);
    for (int i=0; i<iI; i++) Xest[i]=i;

    /* DD ambiguity fixed rate */
    fixNumRatio=fixVarRatio=fix_num=n_elAmb=0;
}
/* update parameters functions ---------------------------------------------------- */
/* update ambiguity parameters ---------------------------------------------------- */
void gnss_relative_c::updateamb() {
    gnss_ssat_c *sss;
    /* count reset number for each frequnecy */
    for (int f=0; f<NFREQ+1; f++) reset_flag[f]=nreset[f]=0;

    /* LC ambiguity */
    if (opt->ionoopt==IONOOPT_IFLC) {
        /* count reset number */
        for (int i=0; i<satnum; i++) {
            if (ssat[usesat[i]-1].reset[NFREQ]>0) nreset[NFREQ]++;
        }
        //reset_flag[NFREQ] = satnum>7 ? (nreset[NFREQ]>satnum-2) : (nreset[NFREQ]>=satnum);
        //reset_flag[NFREQ] = nreset[NFREQ]>satnum-2;
        reset_flag[NFREQ] = nreset[NFREQ]>=satnum;

        reset_x=reset_flag[NFREQ];
        if (reset_x) nsol=0;

        /* update ambiguity parameter to Xpar and Rxvec */
        for (int sys = 0; sys<NSYS; sys++) {
            if (irfsat[sys]<0) continue; //continue if no observation of this system
            if (L_nsat[sys][0]<MINSYSOBS_RTK) continue; //continue if satellite <2
            for (int i=0; i<satnum; i++) {
                if (obsr->data[rovsat[i]].isys!=sys) continue;
                sss = &ssat[usesat[i]-1];

                /* update ambiguity */
                lam=nav->lam[ sss->sat - 1 ];
                rp=&obsr->data[rovsat[i]]; bp=&obsb->data[bassat[i]];
                /* correct antenna and phase windup for observations and ambiguities */
                sss->correct_obs_amb(NFREQ,rp,bp);
                /* reset ambiguity if this frequency should be reset */
                if (reset_flag[NFREQ]&&sss->reset[NFREQ]<1) sss->reset_amb(NFREQ);
                /* update ambiguity parameters */
                if (sss->update_solved_amb(NFREQ,nsol)>0) {
                    Xpar.push_back(sss->amb[NFREQ]);
                    Rxvec.push_back(sss->ambvar[NFREQ]);
                    ambnum[i] = nA++;
                    //Xest index in all parameters
                    Xest.push_back(NX+(usesat[i]-1)*numF);
                    /* if reset initialize covariance matrix of all parameters */
                    if (sss->fix[NFREQ]==1) {
                        init_RxAMB(sss->sat,0);
                        //reset flag
                        Xini.push_back(1);
                    }
                    else Xini.push_back(0);
                }
            }
        }
    }
    /* uncombined ambiguity */
    else {
        /* count reset number */
        for (int i=0; i<satnum; i++) {
            for (int f=0; f<numF; f++) if (ssat[usesat[i]-1].reset[f]>0) nreset[f]++;
        }
        for (int f=0; f<numF; f++) { //reset flag of each frequency
            //reset_flag[f] = satnum>7 ? (nreset[f]>satnum-2) : (nreset[f]>=satnum);
            //reset_flag[f] = nreset[f]>satnum-2;
            reset_flag[f] = nreset[f]>=satnum;
        }

        reset_x=1;
        for (int f=0; f<numF; f++) reset_x=reset_flag[f]&&reset_x;
        if (reset_x) nsol=0;

        /* update ambiguity parameter to Xpar and Rxvec */
        for (int sys = 0; sys<NSYS; sys++) {
            if (irfsat[sys]<0) continue; //continue if no observation of this system
            for (int freq = 0; freq<numF; freq++) {
                if (L_nsat[sys][freq]<MINSYSOBS_RTK) continue; //continue if satellite <2
                for (int i=0; i<satnum; i++) {
                    if (obsr->data[rovsat[i]].isys!=sys) continue;
                    sss = &ssat[usesat[i]-1];

                    /* update ambiguity */
                    lam=nav->lam[ sss->sat - 1 ];
                    rp=&obsr->data[rovsat[i]]; bp=&obsb->data[bassat[i]];
                    /* correct antenna and phase windup for observations and ambiguities */
                    sss->correct_obs_amb(freq,rp,bp);
                    /* reset ambiguity if this frequency should be reset */
                    if (reset_flag[freq]&&sss->reset[freq]<1) sss->reset_amb(freq);
                    /* update ambiguity parameters */
                    if (sss->update_solved_amb(freq,nsol)>0) {
                        Xpar.push_back(sss->amb[freq]);
                        Rxvec.push_back(sss->ambvar[freq]);
                        /* add a noise if detected reseted ambiguity */
                        if (nreset[freq]&&opt->modear==ARMODE_FIXHOLD) Rxvec.back()+=FIXED_AMB_ADJ;
                        ambnum[i + freq*satnum] = nA++;
                        //Xest index in all parameters
                        Xest.push_back(NX+(usesat[i]-1)*numF+freq);
                        /* if reset initialize covariance matrix of all parameters */
                        if (sss->fix[freq]==1) {
                            init_RxAMB(sss->sat,freq);
                            //reset flag
                            Xini.push_back(1);
                        }
                        else Xini.push_back(0);
                    }
                }
            }
        }
    }
}
/* update geodetic parameters ----------------------------------------------------- */
void gnss_relative_c::updatexyz(int prePos) {
    /* tightly integration */
    if ( prePos==1 ) {
        return; /* the geodetic parameters was updated in "gnss_get_posvel()" */
    }
    /* fix mode */
    if (opt->mode==PMODE_FIXED) {
        for (int i=0; i<3; i++) {
            Xpar[i] = opt->ru[i];
            Rxvec[i] = 1E-8;
            Xini[i]=1;
            init_RxALL(i);
        }
        return;
    }
    /* initialize xyz using SPP result for first epoch */
    //if (norm(sol[MAXSOLBUF-2].xpos.begin(),3) <= 0.0) {
    if (norm(X_ALL.begin(),3) <= 0.0) {
        for (int i=0; i<NP; i++) {
            Xpar[i] = sol.back().xpos[i];
            Rxvec[i] =
                i<3 ? VAR_POS : VAR_VEL;
        }
        return;
    }
    /* static mode */
    if (opt->mode==PMODE_STATIC) {
        for (int i=0; i<3; i++) {
            //Xpar[i] = sol[MAXSOLBUF-2].xpos[i];
            Xpar[i] = X_ALL[i];
            if (ntotal<iniamb) {
                Xini[i]=1;
                init_RxALL(i);
                Rxvec[i] = VAR_POS;
            }
            //else Rxvec[i] = sol[MAXSOLBUF-2].vpos[i+i*3];
            else Rxvec[i] = Rx_ALL[i+i*N_ALL];
        }
        return;
    }
    /* kinematic mode without velocity (only xyz) */
    if ( opt->mode==PMODE_KINEMA || opt->mode==PMODE_MOVEB ) {
        for (int i=0; i<3; i++) {
            Xpar[i] = sol.back().xpos[i];
            Rxvec[i] = VAR_POS;
            Xini[i]=1;
            init_RxALL(i);
        }
        return;
    }
}
/* update troposphere parameters -------------------------------------------------- */
void gnss_relative_c::updatetro() {
    // loop rover and base
    for (int i=0; i<2; i++) {
        /* initialize tro parameters for the first epoch */
        //if ( sol[MAXSOLBUF-2].xtro[i*nTrp]==0.0 || ntotal<iniamb ) {
        if ( X_ALL[iT+i*nTrp]==0.0 || ntotal<iniamb ) {
            const double zazel[] ={ 0.0,PI/2.0 };
            gnss_obsd_c obs;
            double *pos = i==0 ? Rblh : Bblh;
            trofunc->saascorr(&obs,pos,zazel,0.7);
            Xpar[iT+i*nTrp] = obs.dtro;
            Rxvec[iT+i*nTrp]=SQR(opt->std[2]);
            Xini[iT+i*nTrp]=1;
            init_RxALL(iT+i*nTrp);
            /* estimate tro gradient */
            if (nTrp==3) for (int j=1; j<3; j++) {
                Xpar[iT+i*nTrp+j] =INIT_GRA;
                Rxvec[iT+i*nTrp+j]=VAR_GRA;
                Xini[iT+i*nTrp+j]=1;
                init_RxALL(iT+i*nTrp+j);
            }
        }
        /* update tro parameters using last solution */
        else {
            //Xpar[iT+i*nTrp] =sol[MAXSOLBUF-2].xtro[i*nTrp];
            //Rxvec[iT+i*nTrp]=sol[MAXSOLBUF-2].vtro[i*nTrp+i*nTrp*NT]+opt->stdrate[2]*tt;
            Xpar[iT+i*nTrp]  = X_ALL[iT+i*nTrp];
            Rxvec[iT+i*nTrp] = Rx_ALL[ (iT+i*nTrp)*(N_ALL+1) ] + opt->stdrate[2]*tt;
            /* estimate tro gradient */
            if (nTrp==3) for (int j=1; j<3; j++) {
                //Xpar[iT+i*nTrp+j] =sol[MAXSOLBUF-2].xtro[i*nTrp+j];
                //Rxvec[iT+i*nTrp+j]=sol[MAXSOLBUF-2].vtro[i*nTrp+j+(i*nTrp+j)*NT]+opt->stdrate[2]*tt*0.01;
                Xpar[iT+i*nTrp+j]  = X_ALL[iT+i*nTrp+j];
                Rxvec[iT+i*nTrp+j] = Rx_ALL[ (iT+i*nTrp+j)*(N_ALL+1) ] + opt->stdrate[2]*tt*0.01;
            }
        }
    }
}
/* update glonass differenced IFB rate between 2 receivers ------------------------ */
void gnss_relative_c::updateglo() {
    for (int i=0; i<NG; i++) {
        /* initialize glo parameters */
        Xpar[iG+i] = INIT_HWBIAS;
        Rxvec[iG+i] = VAR_HWBIAS;
        Xini[iG+i]=1;
        init_RxALL(iG+i);

        /* update glo parameters */
        /*else {
            Xpar[iG+i] = sol[MAXSOLBUF-2].xglo[i];
            Rxvec[iG+i] = sol[MAXSOLBUF-2].vglo[i + i*NG] + SQR(RAT_HWBIAS)*tt;
        }*/
        /*else {
            Xpar[iG+i] = X_ALL[iG+i];
            Rxvec[iG+i] = Rx_ALL[ (iG+i)*(N_ALL+1) ] + SQR(RAT_HWBIAS)*tt;
        }*/
    }
}
/* update ionosphere parameters --------------------------------------------------- */
void gnss_relative_c::updateion() {
    for (int i=0; i<satnum; i++) {
        ionfunc->correction(&obsr->data[rovsat[i]],nav,Rblh);
        ionfunc->correction(&obsb->data[bassat[i]],nav,Bblh);
        /*double dion=obsb->data[ bassat[i] ].dion - obsr->data[ rovsat[i] ].dion;*/
        /*if (nsol<iniamb) {*/
        Xpar[iI+i]=obsb->data[bassat[i]].dion-obsr->data[rovsat[i]].dion;
        Rxvec[iI+i]=3.0*(obsr->data[rovsat[i]].ionvar+obsb->data[bassat[i]].ionvar);
        int index_I = usesat[i]-1;
        Xest[iI+i]=iI+index_I;
        Xini[iI+i]=1;
        init_RxALL(iI+index_I);
    }
}
/* update parameters covariance matrix -------------------------------------------- */
void gnss_relative_c::updatevar() {
    //get covariance from Rx_ALL
    for (int i=0; i<numX; i++) {
        Rx[i + i*numX] = Rxvec[i]; //variance from Rxvec
        if (Xini[i]>0) continue;
        for (int j=0; j<i; j++) {
            if (Xini[j]>0) continue;
            Rx[j+i*numX] = Rx[i+j*numX] = Rx_ALL[ Xest[j] + Xest[i]*N_ALL ];
        }
    }
}
/* update parameters from previous time to current time --------------------------- */
void gnss_relative_c::updatepar(int prePos) {
    iniamb = opt->iniamb;

    /* reset vector and matrix */
    init_arrmat(prePos);

    /* update ssat observation and ambiguity parameters (not for DGPS mode) */
    if (opt->mode>PMODE_DGPS) updateamb();

    /* update position/velocity */
    updatexyz(prePos);

    /* initialize geodetic position of rover and base */
    for (int i=0; i<3; i++) { Rxyz[i] = Xpar[i]; Bxyz[i] = rb[i]; }
    // tidal correction
    if (opt->tidecorr) {
        gtime_c trov = obsr->data[0].time;
        gtime_c tbas = obsb->data[0].time;
        tidefunc.tidecorr(*(trov.timeadd(-obsr->data[0].dtr)->gpst2utc()),0,Rxyz,Rtide);
        tidefunc.tidecorr(*(tbas.timeadd(-obsb->data[0].dtr)->gpst2utc()),1,Bxyz,Btide);
        for (int i=0; i<3; i++) { Rxyz[i] += Rtide[i]; Bxyz[i] += Btide[i]; }
    }
    // geodetic position of rover and base
    ecef2blh(Rxyz,WGS84,Rblh); ecef2blh(Bxyz,WGS84,Bblh);

    /* update troposphere parameters */
    if ( opt->tropopt==TROPOPT_EST || opt->tropopt==TROPOPT_ESTG ) updatetro();

    /* update glonass differenced IFB rate between 2 receivers */
    if (NG>0 && (opt->navsys&SYS_GLO)) updateglo();

    /* update ionosphere parameters */
    if (opt->ionoopt==IONOOPT_CONST) updateion();

    /* update parameters number and Rx */
    numX += nA;
    fix_flag.assign(nA,0);
    Rx.assign(numX*numX,0.0);
    updatevar();
}
/* time difference between rover and base ----------------------------------------- */
int gnss_relative_c::timediff_rb() {
    /* odt:  time difference between obs. of rover and base (s) */
    /* sol.back().age: time difference bewtween rover and base
           = solution time difference if moving base-line
           = odt                                                */
    odt = obsr->data[0].time.timediff(obsb->data[0].time);
    solp->age = odt;
    if (solp->age > opt->maxtdiff) {
        msg = "time difference between rover and base too large!";
        return 0;
    }
    if (opt->mode==PMODE_MOVEB && solp->age > TTOL_MOVEB) {
        msg = "time sync between rover and base failed!";
        return 0;
    }

    return 1;
}
/* satellite-single-differenced dynamic parameters -----------------------------------
* argv   :  int   isat   number of common satellite
* --------------------------------------------------------------------------------- */
double gnss_relative_c::single_distanc(int isat,int sys) {
    /* R:rover B:base */
    double disR,disB;
    Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;

    for (int i=0; i<3; i++) Acoep->at(i) = -rp->sigvec[i];

    /* distance between satellite and rover/base */
    disR = geodist_GNSS(rp->posvel,Rxyz,rp->sigvec);
    disB = geodist_GNSS(bp->posvel,Bxyz,bp->sigvec);

    satazel(Rblh,rp->sigvec,rp->azel);
    satazel(Bblh,bp->sigvec,bp->azel);

    double sddis = disR - disB;

    return sddis;
}
/* satellite-single-differenced troposphere parameters -------------------------------
* argv   :  int   isat   number of common satellite
* --------------------------------------------------------------------------------- */
double gnss_relative_c::single_troppar(int isat,int sys) {
    Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;

    /* troposphere delay correction and Acoe tro vector */
    double sdtro;
    // rover and isat
    for (int i=0; i<nTrp; i++) trofunc->model_est.Ptro[i] = Xpar[iT + i];
    trofunc->correction(rp,Rblh,0.0);
    for (int i=0; i<nTrp; i++) Acoep->at(iT + i) = trofunc->model_est.Atro[i];
    sdtro = rp->dtro;
    // base and isat
    for (int i=0; i<nTrp; i++) trofunc->model_est.Ptro[i] = Xpar[iT + nTrp + i];
    trofunc->correction(bp,Bblh,0.0);
    for (int i=0; i<nTrp; i++) Acoep->at(iT + nTrp + i) -= trofunc->model_est.Atro[i];
    sdtro -= bp->dtro;

    return sdtro;
}
/* satellite-single-differenced GLO amb-difference parameters ------------------------
* argv   :  int   isat   number of common satellite
*           int   fff      frequency number
* --------------------------------------------------------------------------------- */
double gnss_relative_c::single_gloambd(int isat,int fff,int sys) {
    lam = isat==irfsat[sys]? lam1 : lam2;
    Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;

    /* glo freq-difference parameter and Acoe glo vector */
    double sdglo = 0.0;
    Acoep->at(iG + fff) = CLIGHT/lam[fff]/1E6;
    sdglo = Acoep->at(iG + fff)*Xpar[iG + fff];

    return sdglo;
}
/* satellite-single-differenced ionosphere parameters --------------------------------
* argv   :  int   isat   number of common satellite
* --------------------------------------------------------------------------------- */
double gnss_relative_c::single_ionopar(int isat,int sys) {
    double sdion=0.0;
    /* ionosphere correction with certain model */
    if (opt->ionoopt!=IONOOPT_CONST) {
        // rover and isat
        ionfunc->correction(rp,nav,Rblh);
        sdion = rp->dion;
        // base and isat
        ionfunc->correction(bp,nav,Bblh);
        sdion -= bp->dion;
    }

    /* estimate constrained DDion vertical delay of GPS L1 (m) */
    else {
        Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;
        Acoep->at(iI+isat)=(rp->ionmap+bp->ionmap)/2.0;
        sdion=Acoep->at(iI+isat)*Xpar[iI+isat];
    }

    return sdion;
}
/* satellite-single-differenced ambiguity parameters ---------------------------------
* argv   :  int   isat   number of common satellite
*           int   fff      frequency number
* --------------------------------------------------------------------------------- */
double gnss_relative_c::single_ambtpar(int isat,int fff,int sys) {
    Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;

    double sdamb;
    int namb = ambnum[isat + fff*satnum];
    /* LC ambiguity parameters and Acoe amb vector */
    if (opt->ionoopt==IONOOPT_IFLC) {
        Acoep->at(iA + namb) = 1.0;
        sdamb = Xpar[iA + namb];
    }

    /* freq ambiguity parameters and Acoe amb vector */
    else {
        lam = isat==irfsat[sys]? lam1 : lam2;
        Acoep->at(iA + namb) = lam[fff];
        sdamb = lam[fff]*Xpar[iA + namb];
    }

    return sdamb;
}
/* satellite-single-differenced variance ------------------------------------------ */
double gnss_relative_c::single_variance(int isat,int freq) {
    int fff = freq % numF;
    int sys = ssat[usesat[isat]-1].sys,prn = obsr->data[rovsat[isat]].prn;
    double blerr = baseline*opt->err[4]/1E4,sterr = CLIGHT*opt->sclkstab*odt;
    double fact = 1;
    if (sys!=SYS_BDS) fact *= opt->adjustfunc==ADJUST_HELMERT? EFACT_GPS : (sys==SYS_GLO?
        EFACT_GLO : (sys==SYS_SBS? EFACT_SBS : EFACT_GPS) );
    else fact *= binary_search(BDS_GEO,BDS_GEO+NBDS_GEO,prn)? EFACT_BDS_G : EFACT_BDS;
    /*if (freq < numF) switch (sys){
        case SYS_GPS: sfact=1.000; break;
        case SYS_BDS: sfact=1.455; break;
        case SYS_GLO: sfact=1.433; break;
        case SYS_GAL: sfact=1.403; break;
        default: sfact=1.0;
    }
    else if (freq >= numF) switch (sys){
        case SYS_GPS: sfact=1.983; break;
        case SYS_BDS: sfact=1.642; break;
        case SYS_GLO: sfact=2.037; break;
        case SYS_GAL: sfact=0.957; break;
        default: sfact=1.0;
    }*/

    /* elevation factor */
    double sin_elR = rp->azel[1]<opt->con_weight_el? sin(rp->azel[1]) : 1;

    /* IF factor */
    //if (opt->ionoopt==IONOOPT_IFLC) fact*=3.0;

    /* phase variance */
    if (freq<numF) {
        return SQR(fact) * ( SQR(opt->err[2]) + SQR(opt->err[3]/sin_elR) ) + SQR(blerr) + SQR(sterr);
    }
    /* code variance */
    else return SQR(fact) * ( SQR(opt->err[0]) + SQR(opt->err[1]/sin_elR) ) + SQR(blerr) + SQR(sterr);
}
/* irfsat-single-differenced parameters and variance and lam1 --------------------- */
void gnss_relative_c::irfsat_single() {
    for (int sys = 0; sys<NSYS; sys++) {
        Airfsat[sys].assign(numX,0.0);
        //continue if no observation of this system
        if (irfsat[sys]<0) continue;
        /* irfsat gnss_ssat_c and irfsat lamda */
        gnss_ssat_c *shhh = &ssat[rfsat[1][sys] - 1];
        lam1 = nav->lam[rfsat[1][sys] - 1];

        /* initialize rover & base observation pointer */
        rp = &obsr->data[rovsat[irfsat[sys]]]; bp = &obsb->data[bassat[irfsat[sys]]];
        // distance
        disRB1[sys] = single_distanc(irfsat[sys],sys);
        // troposphere
        troRB1[sys] = single_troppar(irfsat[sys],sys);
        // glo freq-difference
        if (NG>0 && satsys(rfsat[1][sys],NULL)==SYS_GLO)
            for (int i=0; i<numF&&i<NFREQGLO; i++) gloRB1[i] = single_gloambd(irfsat[sys],i,sys);
        // ionosphere
        ionRB1[sys] = single_ionopar(irfsat[sys],sys);
        // ambiguity
        if (opt->mode>PMODE_DGPS) for (int freq = 0; freq<numF; freq++) {
            if (L_nsat[sys][freq]<MINSYSOBS_RTK) continue; //continue if satellite <2
            ambRB1[sys][freq] = 0.0;
            /* if slove ambiguity */
            if (ambnum[irfsat[sys]+freq*satnum]>=0)
                ambRB1[sys][freq] = single_ambtpar(irfsat[sys],freq,sys);
        }

        /* irfsat observation variance */
        for (int freq = 0; freq<numF*2; freq++) {
            varRB1[sys][freq] = single_variance(irfsat[sys],freq);
            shhh->vsat[freq%numF] = 1;
        }
    }
}
/* double-differenced value ----------------------------------------------------------
* argv   :  int   isat   number of common satellite
*           int   freq     frequency number
*           int   sys      satellite system (0:GPS/SBS/QZS,1:GLO,2:GAL,3:BDS)
* process:  compute double-differenced value and its Acoe vector
* --------------------------------------------------------------------------------- */
double gnss_relative_c::double_value(int isat,int freq,int sys) {
    /* observation type */
    int fff = freq % numF;
    ifact1 = (freq<numF?-1.0:1.0)*SQR(lam1[fff]/GNSS_LAMBDA_M[iGPS][1]);
    ifact2 = (freq<numF?-1.0:1.0)*SQR(lam2[fff]/GNSS_LAMBDA_M[iGPS][1]);

    /* single difference of isat */
    // distance
    disRB2 = single_distanc(isat,sys);
    // troposphere
    troRB2 = single_troppar(isat,sys);
    // glo freq-difference
    gloRB2 = 0;
    if (NG>0 && freq<numF&&satsys(usesat[isat],NULL)==SYS_GLO && fff<NFREQGLO) {
        gloRB2 = single_gloambd(isat,fff,sys);
        gloRB2 = gloRB1[fff] - gloRB2;
    }
    // ionosphere
    ionRB2 = single_ionopar(isat,sys);
    // ambiguity
    ambRB2 = 0;
    if (freq<numF) {
        if ( opt->ionoopt!=IONOOPT_IFLC && opt->modear==ARMODE_LCWN ) ambRB2=LC_WN(isat,sys);
        else {
            ambRB2 = single_ambtpar(isat,fff,sys);
            ambRB2 = ambRB1[sys][fff] - ambRB2;
        }
    }

    /* return double-differenced parameters */
    return disRB1[sys] - disRB2 +
        ifact1*ionRB1[sys] - ifact2*ionRB2 +
        troRB1[sys] - troRB2 +
        gloRB2 +
        ambRB2;
}
/* update Rvec vector ------------------------------------------------------------- */
void gnss_relative_c::update_Rvec(int isat,int freq,int sys) {
    /* distance variance */
    Rvec.push_back(varRB1[sys][freq]);
    Rvec.push_back(single_variance(isat,freq));
}
/* update Lobs vector and lam2 ---------------------------------------------------- */
int gnss_relative_c::update_Lobs(int isat,int freq,int sys) {
    double Lo = 0.0;
    int fff = freq % numF;

    /* stanum sat */
    unsigned int sat = usesat[isat];
    gnss_ssat_c *sss = &ssat[sat - 1],*shhh = &ssat[rfsat[1][sys] - 1];

    /* phase observation */
    if (freq<numF) {
        if (ambnum[irfsat[sys] + fff*satnum]>=0 && ambnum[isat + fff*satnum]>=0) {
            /* LC combination */
            if (opt->ionoopt==IONOOPT_IFLC) {
                Lo = shhh->lc12[1] - sss->lc12[1];
            }
            /* Lfreq */
            else {
                Lo = shhh->L[fff].back()*lam1[fff] - sss->L[fff].back()*lam2[fff];
            }
        }
        else return 0;
    }
    /* code observation */
    else {
        double p1 = 0.0,p2 = 0.0;
        /* LC combination */
        if (opt->ionoopt==IONOOPT_IFLC) {
            p1 = shhh->pc12[1];
            p2 = sss->pc12[1];
        }
        /* Lfreq */
        else {
            p1 = shhh->P[fff].back();
            p2 = sss->P[fff].back();
        }
        Lo = p1==0.0 || p2==0.0 ? 0.0 : p1 - p2;
    }

    if (Lo==0.0) return 0;

    Lobs.push_back(Lo);
    return 1;
}
/* update Acoe matrix ------------------------------------------------------------- */
void gnss_relative_c::update_Acoe(int isat,int freq,int sys) {
    /* frequency index */
    int fff = freq % numF;

    /* update Acoe matrix */
    // dynamic pararmeter coefficients
    for (int nx=0; nx<NP; nx++)
        Ais[nx] = Airfsat[sys][nx] - Asatnum[nx];
    // troposphere parameter coefficients
    for (int nx = 0; nx<NT; nx++)
        Ais[iT + nx] = Airfsat[sys][iT + nx] - Asatnum[iT + nx];
    // glo freq-difference parameter coefficients
    if (NG>0 && freq<numF&&sys==1 && fff<NFREQGLO)
        Ais[iG + fff] = Airfsat[sys][iG + fff] - Asatnum[iG + fff];
    // ionosphere parameter coefficients
    if (opt->ionoopt==IONOOPT_CONST) {
        Ais[iI + irfsat[sys]] =  ifact1*Airfsat[sys][iI + irfsat[sys]];
        Ais[iI + isat]      = -ifact2*Asatnum[iI + isat];
    }
    // ambiguity parameter coefficients
    if (opt->mode>PMODE_DGPS&&freq<numF) {
        int ah = ambnum[irfsat[sys] + fff*satnum],as = ambnum[isat + fff*satnum];
        if (ah>=0) Ais[iA + ah] = Airfsat[sys][iA + ah];
        if (as>=0) Ais[iA + as] = -Asatnum[iA + as];
    }
    Acoe.insert(Acoe.end(),Ais.begin(),Ais.end());
}
/* double-differenced LC ambiguity using wide-lane and N1 ----------------- */
double gnss_relative_c::LC_WN(int isat,int sys) {
    gnss_ssat_c *sss=&ssat[usesat[isat]-1],*hhh=&ssat[rfsat[1][sys]-1];
    lam=nav->lam[usesat[isat]-1];
    /* float and integar DD wide-lane and narrow-lane (N1) value */
    double f_wide=hhh->mw12[0]-sss->mw12[0],f_narr=hhh->amb_ave[0]-sss->amb_ave[0];
    double fix_wide,fix_narr;

    /* final wide-lane and narrow-lane (N1) ambiguity */
    double coe_narr=1.0/(1.0/lam[0]+1.0/lam[1]),
        coe_wide=1.0/(lam[1]/lam[0]-1.0)/(1.0/lam[0]+1.0/lam[1]);
    int wide=round(f_wide);
    int flag_wide=fabs(wide-f_wide)<0.25&&
        (hhh->ave_con[NFREQ]>iniamb&&sss->ave_con[NFREQ]>iniamb);
    fix_wide = /*flag_wide? wide : */f_wide;
    fix_narr = f_narr;

    /* DD LC ambiguity */
    double DDLC=coe_wide*fix_wide + coe_narr*fix_narr;

    /* update Xpar,its covariance matrix and coefficient vector */
    int nsss=iA+ambnum[isat],nhhh=iA+ambnum[irfsat[sys]];
    Asatnum[nsss]=1.0;
    /* fixed information */
    if (flag_wide) {
        Xpar[nsss]=Xpar[nhhh]-DDLC;
        fix_flag[nsss-iA]=fix_flag[nhhh-iA]=1; fix_num++;
        Rx[nsss+nsss*numX]=Rx[nhhh+nhhh*numX]=SQR(0.007);
    }

    return DDLC;
}
/* baseline-constraint equation for moving-baseling ------------------------------- */
void gnss_relative_c::base_line() {
    Lobs.push_back(opt->baseline[0] - baseline);
    Acoe.insert(Acoe.end(),numX,0.0);
    for (int i=0; i<3; i++)
        Acoe[i + numL*numX] = (Rxyz[i] - Bxyz[i])/baseline;
    Rbl = SQR(opt->baseline[1]);
    numL++;
}
/* update Rvar matrix ------------------------------------------------------------- */
void gnss_relative_c::update_Rvar() {
    /* initialize Rarr */
    Rvar.assign(numL*numL,0.0);

    for (int sys=0,nt=0; sys<NSYS; sys++) {
        for (int i=opt->mode>PMODE_DGPS? 0:numF; i<numF*2; i++) {
            for (int j=0; j<nobs[sys][i]; j++) {
                for (int k=0; k<nobs[sys][i]; k++) {
                    Rvar[ nt+k + (nt+j)*numL ] = j==k ?
                        Rvec[ 2*(nt+j) ] + Rvec[ 2*(nt+j) + 1 ] :
                        Rvec[ 2*(nt+j) ];
                }
            }
            nt += nobs[sys][i];
        }
    }

    /* baseline variance */
    if (opt->mode==PMODE_MOVEB && opt->baseline[0] > 0.0)
        Rvar[numL*numL - 1] = Rbl;
}
/* double-differenced observation equation -------------------------------------------
* return : number of used observation (size of L)
* --------------------------------------------------------------------------------- */
int gnss_relative_c::double_diff() {
    ns = numL = n_Damb = En_Damb = 0;
    for (int i=0; i<NFREQ*2; i++) {
        iobs[i].assign(satnum,-1);
        for (int j=0; j<NSYS; j++) nobs[j][i]=0;
    }

    /* initialize rover and base position */
    if (norm(Rxyz,3) <= 0.0 || norm(Bxyz,3) <= 0.0) return 0;
    /* initialize base-line parameters */
    baseline = distance(Rxyz,Bxyz,3);

    /* single-difference for reference satellite (irfsat) */
    irfsat_single();

    /* loop of different satellites */
    for (int sys=0; sys<NSYS; sys++) {
        //continue if no observation of this system
        if (irfsat[sys]<0) continue;
        lam1 = nav->lam[rfsat[1][sys] - 1];
        for (int freq = opt->mode>PMODE_DGPS ? 0 : numF; freq<numF*2; freq++) {
            if (opt->mode>PMODE_DGPS&&L_nsat[sys][freq%numF]<MINSYSOBS_RTK) continue; //continue if satellite <2
            for (int isat = 0; isat<satnum; isat++) {
                /* continue if it is different system or reference satellite */
                if (obsr->data[rovsat[isat]].isys!=sys || isat==irfsat[sys])
                    continue;

                /* frequency wavelength */
                lam2 = nav->lam[usesat[isat]-1];

                /* initialize rover & base observation pointer */
                rp = &obsr->data[rovsat[isat]]; bp = &obsb->data[bassat[isat]];

                /* update Lobs vector and lam2 */
                if (!update_Lobs(isat,freq,sys))
                    continue;

                /* initialize Lobs,Rvec,Acoe for new observation */
                /* Rvec.insert(Rvec.end(),2,0.0); initialize in update_Rvec */
                Asatnum.assign(numX,0.0);
                Ais.assign(numX,0.0);

                /* compute double-differenced correction between rover/base and irfsat/isat */
                correction = double_value(isat,freq,sys);

                /* update Lobs.back() */
                Lobs.back() -= correction;

                /* update Aceo matrix */
                update_Acoe(isat,freq,sys);

                /* update Rvec vector */
                update_Rvec(isat,freq,sys);

                /* isat sat */
                unsigned int sat = usesat[isat];
                ssat[sat-1].vsat[freq%numF] = 1;
                iobs[freq][isat]=numL++;
                nobs[sys][freq]++;
                /* ambiguity elevation */
                if (freq<numF) {
                    ambElevation[0].push_back(rp->azel[1]);
                    ambElevation[1].push_back(bp->azel[1]);
                }
            }
            if (freq<numF) n_Damb += L_nsat[sys][freq]-1;
            if (freq<numF&&sys!=0) En_Damb += L_nsat[sys][freq]-1;
        }
    }

    /* add baseline equation if moving-baseline */
    if (opt->mode==PMODE_MOVEB && opt->baseline[0] > 0.0) base_line();

    /* update Rvar using Rvec */
    update_Rvar();
    return 1;
}
/* single to double-difference tansformation matrix (Dx_coe) -------------- */
int gnss_relative_c::single2doul(vector<double> &Dx_coe) {
    int ddanum = 0;
    /* initialize Dx_coe */
    for (int i=0; i<iA; i++) Dx_coe[i + i*numX] = 1.0;

    /* Dx_coe for ambiguity parameters */
    for (int sys = 0; sys<NSYS; sys++) {
        if (irfsat[sys]<0) continue; //continue if no observation of this system
        for (int freq = 0; freq<numF; freq++)
            if (L_nsat[sys][freq]>=MINSYSOBS_RTK) {
                for (int isat = 0; isat<satnum; isat++) {
                    if (obsr->data[rovsat[isat]].isys!=sys ||
                        ambnum[isat + freq*satnum]<0 || isat==irfsat[sys]) continue;
                    //reference satellite amb
                    Dx_coe[iA + ambnum[irfsat[sys] + freq*satnum] + (iA + ddanum)*numX] = 1.0;
                    //other satellite amb
                    Dx_coe[iA + ambnum[isat + freq*satnum] + (iA + ddanum)*numX] = -1.0;
                    ddanum++;
                }
            }
    }

    return ddanum;
}
/* double-differenced ambiguity to single ambiguity ----------------------- */
int gnss_relative_c::ddamb2single(vector<double> &Dx_coe,vector<double> &Damb_Fix,
    vector<double> &R_Damb) {

    /*int nLa = nA + n_Damb; //number of pseudo-observation of ambiguity
    vector<double> Acoe(nLa*nA,0.0),Lamb(nLa,0.0),RLa(nLa*nLa,0.0),Ramb(nA*nA,0.0),da(nA,0.0);

    // Acoe
    for (int i=0; i<nA; i++) Acoe[ i + i*nA ] = 1.0; //single amb. pseudo-obs. Acoe
    for (int i=0; i<n_Damb; i++) for (int j = 0; j<nA; j++)
        Acoe[ j + (i+nA)*nA ] = Dx_coe[ j+iA + (i+iA)*numX ]; //fixed DDamb. pseudo-obs. Acoe
    // Lamb
    for (int i=0; i<nA; i++) Lamb[i] = 0.0; //single amb. pseudo-obs. Lamb
    for (int i=0; i<n_Damb; i++) Lamb[ i + nA ] = -Damb_Fix[i]; //fixed DDamb. pseudo-obs. Lamb
    // RLa
    for (int i=0; i<nA; i++) for (int j = 0; j<nA; j++)
        RLa[ j + i*nLa ] = Rx[ j+iA + (i+iA)*numX ]; //single amb. pseudo-obs. RLa
    for (int i=0; i<n_Damb; i++) RLa[ i+nA + (i+nA)*nLa ] = 0.001; //fixed DDamb. pseudo-obs. RLa

    ambFilter.LSQ(Acoe,Lamb,RLa,da,Ramb,nLa,numX);
    for (int i=0; i<nA; i++) fix_Xpar[i + iA] += da[i];
    */

    int nLa = n_Damb; //number of pseudo-observation of ambiguity
    vector<double> Acoe(nLa*numX,0.0),Lamb(nLa,0.0),RLa(nLa*nLa,0.0),da(numX,0.0);

    // Acoe
    for (int i=0; i<nLa; i++) for (int j = 0; j<nA; j++)
        Acoe[ j+iA + i*numX ] = Dx_coe[ j+iA + (i+iA)*numX ]; //fixed DDamb. pseudo-obs. Acoe
    // Lamb
    for (int i=0; i<nLa; i++) Lamb[i] = -Damb_Fix[i]; //fixed DDamb. pseudo-obs. Lamb
    // RLa
    for (int i=0; i<nLa; i++) { //fixed DDamb. pseudo-obs. RLa
        RLa[ i + i*nLa ] = ( ambElevation[0][i]<opt->minElArHold || ambElevation[1][i]<opt->minElArHold )?
            SQR(R_Damb[ i + i*nLa ]) : FIXED_AMB_VAR;
    }

    /* compute least-square solution */
    int ambGood=ambFilter.KF(Acoe,Lamb,RLa,da,fix_Rx,nLa,numX);

    for (int i=0; i<numX; i++) fix_Xpar[i] += da[i];

    return ambGood;
}
/* get fixed solution ----------------------------------------------------- */
int gnss_relative_c::get_fixsolution() {
    if (n_Damb <= 0) return 0;
    int Dnum = iA + n_Damb;

    /* solve double-differenced ambiguity vectors */
    /* Dx      : double-differenced format Xpar
    *Dx_Coe  : Xpar to Dx Coefficients
    *DR      : Dx_Coe*Rx
    *Damb    : float double-differenced ambiguity
    *fix_amb : fixed double-differenced ambigtuiy
    *Damb_Fix: Damb - fix_amb
    *R_xa    : covariance of ambiguity and other parameters
    *P_Damb  : inverse of R_Damb
    *d_Damb  : P_Damb*Damb_Fix
    *RR      : R_xa*P_Damb */
    vector<double> Dx_Coe(Dnum*numX,0.0),DR(numX*Dnum,0.0),R_Dx(Dnum*Dnum,0.0),Dx(Dnum,1),
        R_Damb(n_Damb*n_Damb,0.0),Damb(n_Damb,0.0),Damb_Fix(n_Damb,0.0),
        R_xa(iA*n_Damb,0.0),P_Damb,d_Damb(n_Damb,0.0),RR(n_Damb*iA,0.0);

    if (single2doul(Dx_Coe)!=n_Damb) return 0; //number of DD ambiguity

    /* initialize of fixed solution */
    fix_amb.assign(n_Damb*2,0.0);

    /* Dx,Damb,DR,R_Dx,R_Damb,R_xa */
    matmul_vec("NN",Dnum,1,numX,1.0,Dx_Coe,Xpar,0,Dx); //Dx
    matmul_vec("NN",Dnum,numX,numX,1.0,Dx_Coe,Rx,0,DR); //DR
    matmul_vec("NT",Dnum,Dnum,numX,1.0,DR,Dx_Coe,0,R_Dx); //R_Dx

    for (int i=0; i<n_Damb; i++) {
        Damb[i] = Dx[iA + i]; //Damb
        for (int j=0; j<n_Damb; j++) R_Damb[j + i*n_Damb] = R_Dx[j+iA + (i+iA)*Dnum]; //R_Damb
        for (int j=0; j<iA; j++) R_xa[j + i*iA] = R_Dx[j + (i+iA)*Dnum]; //R_xa
    }

    /* fixed DD ambiguity solution using lambda/mlambda */
    int flag = lambda.int_amb(Damb,R_Damb,fix_amb,n_Damb,2,sum_var);

    /* if succeed */
    if ( flag==0 ) {
        fixVarRatio = sum_var[0]==0? 1E7 : sum_var[1]/sum_var[0];
        if (fixVarRatio>=opt->thresar[0]) {
            /* initialize of fixed solution */
            fix_Xpar.assign(Xpar.begin(),Xpar.end()); fix_Rx.assign(Rx.begin(),Rx.end());
            for (int i=0; i<n_Damb; i++) Damb_Fix[i] = Damb[i] - fix_amb[i];
            if (opt->modear==ARMODE_FIXHOLD) {
                /* double-differenced ambiguity to single ambiguity */
                if (ddamb2single(Dx_Coe,Damb_Fix,R_Damb)==0) return 1;
            }
            else {
                /* transform float to fixed solution (fix_Xpar = Xpar - R_xa'*P_Damb*Damb_Fix) */
                P_Damb.assign(R_Damb.begin(),R_Damb.end());
                if (matinv(P_Damb,n_Damb)==0) {
                    matmul_vec("NN",n_Damb,1,n_Damb,1.0,P_Damb,Damb_Fix,0,d_Damb);
                    vector<double> dxa(iA,0);
                    matmul_vec("TN",iA,1,n_Damb,-1.0,R_xa,d_Damb,0,dxa);
                    matmul_vec("TN",iA,1,n_Damb,-1.0,R_xa,d_Damb,1.0,fix_Xpar);
                    /* covariance of fixed solution (fix_Rx = Rx - R_xa'*P_Damb*R_xa) */
                    matmul_vec("TN",iA,n_Damb,n_Damb,1.0,R_xa,P_Damb,0,RR);
                    vector<double> fix_RiA(iA*iA,0.0);
                    matmul_vec("NN",iA,iA,n_Damb,-1.0,RR,R_xa,0,fix_RiA);
                    for (int i=0; i<iA; i++) for (int j=0; j<iA; j++) fix_Rx[j+i*numX]+=fix_RiA[j+i*iA];

                    return 1;
                }
            }
        }
    }
    return 0;
}
/* single to double-difference tansformation matrix (considering elevation) ------- */
int gnss_relative_c::single2doul_el(vector<double> &Dx_coe) {
    int ddanum=0;
    n_elAmb=0;
    fixIndex.clear();

    /* Dx_coe for ambiguity parameters */
    for (int sys = 0; sys<NSYS; sys++) {
        if (irfsat[sys]<0) continue; //continue if no observation of this system
        for (int freq = 0; freq<numF; freq++)
            if (L_nsat[sys][freq]>=MINSYSOBS_RTK) {
                for (int isat = 0; isat<satnum; isat++) {
                    if (obsr->data[rovsat[isat]].isys!=sys ||
                        ambnum[isat + freq*satnum]<0 || isat==irfsat[sys]) continue;

                    if (ambElevation[0][ddanum]>opt->minElArFix && ambElevation[1][ddanum]>opt->minElArFix) {
                        fixIndex.push_back(ddanum);
                        //reference satellite amb
                        Dx_coe[iA + ambnum[irfsat[sys] + freq*satnum] + n_elAmb*numX] = 1.0;
                        //other satellite amb
                        Dx_coe[iA + ambnum[isat + freq*satnum] + n_elAmb*numX] = -1.0;
                        n_elAmb++;
                    }
                    ddanum++;
                }
            }
    }

    return ddanum;
}
/* double-differenced ambiguity to single ambiguity (considering elevation) ------- */
int gnss_relative_c::ddamb2single_el(vector<double> &Dx_coe,vector<double> &Damb_Fix,
    vector<double> &R_Damb,vector<double> &Rxa,vector<double> &dxa) {

    int nLa = n_elAmb; //number of pseudo-observation of ambiguity
    vector<double> Lamb(nLa,0.0),RLa(nLa*nLa,0.0);

    // Lamb
    for (int i=0; i<nLa; i++) Lamb[i] = -Damb_Fix[i]; //fixed DDamb. pseudo-obs. Lamb
    // RLa
    for (int i=0; i<nLa; i++) { //fixed DDamb. pseudo-obs. RLa
        RLa[ i + i*nLa ] = ( ambElevation[0][fixIndex[i]]<opt->minElArHold || 
            ambElevation[1][fixIndex[i]]<opt->minElArHold )?
            SQR(R_Damb[ i + i*nLa ]) : FIXED_AMB_VAR;
    }

    /* compute least-square solution */
    int ambGood=ambFilter.KF(Dx_coe,Lamb,RLa,dxa,Rxa,nLa,numX);

    return ambGood;
}
/* get fixed solution (considering elevation) ------------------------------------- */
int gnss_relative_c::get_fixsolution_el() {
    if (n_Damb <= 0) return 0;

    /* solve double-differenced ambiguity with high elevation */
    /* Dx      : double-differenced format Xpar
     * Dx_Coe  : Xpar to Dx Coefficients
     * DR      : Dx_Coe*Rx
     * Damb    : float double-differenced ambiguity
     * fix_amb : fixed double-differenced ambigtuiy
     * Damb_Fix: Damb - fix_amb
     * Rxa     : fixed covariance of ambiguity and other parameters */
    vector<double> Dx_Coe(n_Damb*numX,0.0);

    if ( ambElevation[0].size()!=n_Damb || single2doul_el(Dx_Coe)!=n_Damb || fixIndex.size()!=n_elAmb ) 
        return 0;

    if (n_elAmb<n_Damb)
        n_elAmb=n_elAmb;

    vector<double> DR(n_elAmb*numX,0.0),R_Damb(n_elAmb*n_elAmb,0.0),Damb(n_elAmb,0.0),
        fix_amb(n_elAmb*2,0.0),Damb_Fix(n_elAmb,0.0),dxa(numX,0);

    /* Damb,DR,R_Damb */
    matmul_vec("NN",n_elAmb,1,numX,1.0,Dx_Coe,Xpar,0,Damb); //Damb
    matmul_vec("NN",n_elAmb,numX,numX,1.0,Dx_Coe,Rx,0,DR); //DR
    matmul_vec("NT",n_elAmb,n_elAmb,numX,1.0,DR,Dx_Coe,0,R_Damb); //R_Damb

    /* fixed DD ambiguity solution using lambda/mlambda */
    int flag = lambda.int_amb(Damb,R_Damb,fix_amb,n_elAmb,2,sum_var);

    /* if succeed */
    if ( flag==0 ) {
        fixVarRatio = sum_var[0]==0? 1E7 : sum_var[1]/sum_var[0];
        if (fixVarRatio>=opt->thresar[0]) {
            /* initialize of fixed solution */
            fix_Xpar.assign(Xpar.begin(),Xpar.end()); fix_Rx.assign(Rx.begin(),Rx.end());

            for (int i=0; i<n_elAmb; i++) Damb_Fix[i] = Damb[i] - fix_amb[i];

            /* double-differenced ambiguity to single ambiguity */
            if (ddamb2single_el(Dx_Coe,Damb_Fix,R_Damb,fix_Rx,dxa)==0) {
                for (int i=0; i<numX; i++) fix_Xpar[i] += dxa[i];
                if (opt->modear==ARMODE_CONT) {
                    for (int i=iA; i<numX; i++) {
                        for (int j=iA; j<numX; j++) fix_Rx[j+i*numX] = Rx[j+i*numX];
                    }
                }
                return 1;
            }
        }
    }
    return 0;
}
/* get fixed wide-lane ambiguity -------------------------------------------------- */
int gnss_relative_c::fix_wide_narr(vector<double> &DLC,vector<double> &R_DLC) {
    vector<double> f_wide(n_Damb,0.0),fix_wide(n_Damb,0.0);
    vector<double> f_narr(n_Damb,0.0),fix_narr(2*n_Damb,0.0),narr(n_Damb,0.0); //N1
    vector<int>    lock_flag(n_Damb,0);
    vector<double> coe_LC(n_Damb,0.0),coe_wide(n_Damb,0.0);
    int ndd=0;

    /* fixed double-wide-lane ambiguity */
    for (int sys = 0; sys<NSYS; sys++) {
        if (irfsat[sys]<0) continue; //continue if no observation of this system
        gnss_ssat_c *sss,*hhh=&ssat[rfsat[1][sys] - 1];
        if (L_nsat[sys][0]>=MINSYSOBS_RTK) for (int isat = 0; isat<satnum; isat++) {
            if ( obsr->data[rovsat[isat]].isys!=sys ||
                ambnum[isat]<0 || isat==irfsat[sys]) continue;
            sss=&ssat[usesat[isat]-1];
            double *lami=nav->lam[usesat[isat]-1],*lamh=nav->lam[rfsat[1][sys]-1];
            /* float wide-lane */
            f_wide[ndd] = hhh->mw12[0] - sss->mw12[0];
            f_narr[ndd] = hhh->amb_ave[0] - sss->amb_ave[0];
            /* coe_LC = 1/lam1+1/lam2,coe_wide = 1/(lam2/lam1-1) */
            coe_LC[ndd] = 1.0/lami[0] + 1.0/lami[1];
            coe_wide[ndd] = 1.0/(lami[1]/lami[0]-1.0);
            /* lock count */
            if (sss->ave_con[NFREQ]>iniamb&&hhh->ave_con[NFREQ]>iniamb) lock_flag[ndd]=1;
            /* number of dd ambiguity */
            ndd++;
        }
    }
    /* fixed wide/narrow-lane ambiguity with nearest int value */
    for (int i=0; i<n_Damb; i++) {
        int wide=round(f_wide[i]),narr=round(f_narr[i]);
        fix_wide[i] = fabs(wide-f_wide[i])<0.35&&lock_flag[i]? wide : f_wide[i];
        fix_narr[i] = fabs(narr-f_narr[i])<0.35&&lock_flag[i]? narr : f_narr[i];
    }

    /* fixed LC ambiguity */
    for (int i=0; i<n_Damb; i++)
        fix_amb[i]=1.0/coe_LC[i]*fix_narr[i]+coe_wide[i]/coe_LC[i]*fix_wide[i];

    return 1;
}
/* get LC ambiguity (wide-narrow lane) fixed solution ----------------------------- */
int gnss_relative_c::get_LCfixed() {
    if (n_Damb <= 0) return 0;
    int Dnum = iA + n_Damb; //dd-parameters number
    vector<double> Dx_Coe(Dnum*numX,0.0),fix_LC(n_Damb,0.0);
    /*number of DD ambiguity */
    if (single2doul(Dx_Coe)!=n_Damb) return 0;

    /* solve double-differenced ambiguity vectors */
    /* Dx      : double-differenced format Xpar
    *Dx_Coe  : Xpar to Dx Coefficients
    *DR      : Dx_Coe*Rx
    *DLC     : float double-differenced LC ambiguity
    *amb_fix : fixed double-differenced LC ambigtuiy
    *DLC_Fix : DDLC - amb_fix
    *R_xa    : covariance of ambiguity and other parameters
    *P_DLC   : inverse of R_Damb
    *d_DLC   : P_Damb*Damb_Fix
    *RR      : R_xa*P_Damb */
    vector<double> DR(numX*Dnum,0.0),R_Dx(Dnum*Dnum,0.0),Dx(Dnum,1),
        R_DLC(n_Damb*n_Damb,0.0),DLC(n_Damb,0.0),DLC_Fix(n_Damb,0.0),
        R_xa(iA*n_Damb,0.0),P_DLC,d_DLC(n_Damb,0.0),RR(n_Damb*iA,0.0);

    /* initialize of fixed solution */
    fix_Xpar.assign(Xpar.begin(),Xpar.end()); fix_Rx.assign(Rx.begin(),Rx.end());
    fix_amb.assign(n_Damb,0.0);

    /* Dx,Damb,DR,R_Dx,R_Damb,R_xa */
    matmul_vec("NN",Dnum,1,numX,1.0,Dx_Coe,Xpar,0,Dx); //Dx
    matmul_vec("NN",Dnum,numX,numX,1.0,Dx_Coe,Rx,0,DR); //DR
    matmul_vec("NT",Dnum,Dnum,numX,1.0,DR,Dx_Coe,0,R_Dx); //R_Dx

    for (int i=0; i<n_Damb; i++) {
        DLC[i] = Dx[iA + i]; //Damb
        for (int j=0; j<n_Damb; j++) R_DLC[i+j*n_Damb] = R_Dx[i+iA+(j+iA)*Dnum]; //R_Damb
        for (int j=0; j<iA; j++) R_xa[j+i*iA] = R_Dx[j+(i+iA)*Dnum]; //R_xa
    }

    /* get fixed wide-narrow-lane ambiguity and update LC ambiguity */
    if (fix_wide_narr(DLC,R_DLC)) {
        /* update solution */
        /* transform float to fixed solution (fix_Xpar = Xpar - R_xa'*P_Damb*Damb_Fix) */
        for (int i=0; i<n_Damb; i++) DLC_Fix[i] = DLC[i] - fix_amb[i];
        P_DLC.assign(R_DLC.begin(),R_DLC.end());
        if (matinv(P_DLC,n_Damb)==0) {
            matmul_vec("NN",n_Damb,1,n_Damb,1.0,P_DLC,DLC_Fix,0,d_DLC);
            matmul_vec("TN",iA,1,n_Damb,-1.0,R_xa,d_DLC,1.0,fix_Xpar);
        }
        /* covariance of fixed solution (fix_Rx = Rx - R_xa'*P_Damb*R_xa) */
        matmul_vec("TN",iA,n_Damb,n_Damb,1.0,R_xa,P_DLC,0,RR);
        matmul_vec("NN",iA,iA,n_Damb,-1.0,RR,R_xa,1.0,fix_Rx);

        /* double-differenced ambiguity to single ambiguity */
        ddamb2single(Dx_Coe,DLC_Fix,R_DLC);

        return 1;
    }

    return 0;
}
/* update baseline status --------------------------------------------------------- */
void gnss_relative_c::update_baseline() {
    double rr[3];

    /* update baseline length */
    for (int i=0; i<3; i++) rr[i]=sPar[i]-rb[i];
    baseline=norm(rr,3);

    if ( ( opt->mode==PMODE_MOVEB || opt->baseAtt ) && baseline>1E-3 ) {
        if (baseline<1E-3) {
            for (int i=0; i<3; i++) {
                att_n_b[i]=0;
                att_err[i]=360;
            }
            return;
        }
        double pos[3],enu[3],ecefVar[9],enuVar[9];

        ecef2blh(rb,WGS84,pos);
        ecef2enu(pos,rr,enu);
        /* for test */
        /*enu[0]-=9; enu[1]+=68.4; enu[2]-=2.5;
        baseline=norm(enu,3);*/

        /* attitude */
        double bhor=norm(enu,2);
        att_n_b[2] = -atan2(enu[0],enu[1])*R2D;  //yaw
        att_n_b[0] =  atan2(enu[2],bhor)*R2D;    //pitch
        att_n_b[1] =  0;                         //roll

        /* variance */
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
            ecefVar[j+i*3] = *(sRx + j+i*numX);
        }
        covenu(pos,ecefVar,enuVar);
        // yaw variance
        double errArcY=SQRT(max(enuVar[0],enuVar[4]));
        att_err[2] = errArcY/bhor * R2D;
        if (att_err[2]<0 || att_err[2]>MAX_ATT_ERR ) att_err[2]=360;
        // pitch
        double errArcP=SQRT(max(enuVar[0]+enuVar[4],enuVar[8]));
        att_err[0] = errArcP/baseline * R2D;
        if (att_err[0]<0 || att_err[0]>MAX_ATT_ERR ) att_err[0]=360;
        // roll
        att_err[1]=360;
    }
}
/* calculate pdop value of each satellite system ---------------------------------- */
void gnss_relative_c::get_pdop() {
    for (int i=0; i<NSYS+1; i++) Pdop[i]=0.0;
    double pdop_val=0.0;
    if (Acoe.size()==numL*numX) {
        vector<double> Asys[NSYS+1];
        QAA.assign(9,0.0); Asys[NSYS].assign(0,0.0);
        int allobs=0;
        /* pdop for each single system */
        for (int sys=0,nownum=0; sys<NSYS; sys++) {
            int sysnum=0; int codenum=nobs[sys][numF]+1;
            /* observation number */
            for (int fff=0; fff<numF*2; fff++) sysnum+=nobs[sys][fff];
            /* only available for constellation with >1 satellites */
            if (codenum>1) {
                //update Acoe for all systems
                //reference satellite
                for (int j=0; j<3; j++) Asys[NSYS].push_back( Airfsat[sys][j] );
                //other satellites
                for (int i=1; i<codenum; i++) for (int j=0; j<3; j++)
                    Asys[NSYS].push_back( Airfsat[sys][j] - Acoe[j+(i+nownum)*numX] );
                allobs+=codenum;
            }
            /* pdop of this system */
            if (codenum>=3) {
                //Asys of Acoe to calculate pdop
                Asys[sys].assign(codenum*3,0.0);
                //reference satellite
                for (int j=0; j<3; j++) Asys[sys][j] = Airfsat[sys][j];
                //other satellites
                for (int i=1; i<codenum; i++) for (int j=0; j<3; j++)
                    Asys[sys][j+i*3] = Airfsat[sys][j] - Acoe[j+(i+nownum)*numX];
                //calculate QAA=(Asys'*Asys)^-1
                matmul_vec("TN",3,3,codenum,1.0,Asys[sys],Asys[sys],0,QAA);
                if (matinv(QAA,3)==-1) QAA[0]=QAA[4]=QAA[8]=333.3;
                //calculate pdop 
                double pdop_val=QAA[0]+QAA[4]+QAA[8];
                Pdop[sys] = pdop_val>999.0 ? 999.0 : SQRT(pdop_val);
            }
            else if ( codenum>0 ) Pdop[sys]=999.0;
            nownum+=sysnum;
        }
        /* pdop of all system */
        if (allobs>3) {
            //calculate QAA=(Asys[NSYS]'*Asys[NSYS])^-1
            matmul_vec("TN",3,3,allobs,1.0,Asys[NSYS],Asys[NSYS],0,QAA);
            if (matinv(QAA,3)==-1) QAA[0]=QAA[4]=QAA[8]=333.3;
            //calculate pdop 
            double pdop_val=QAA[0]+QAA[4]+QAA[8];
            Pdop[NSYS] = pdop_val>999.0 ? 999.0 : SQRT(pdop_val);
        }
    }
}
/* average S/N of different systems ----------------------------------------------- */
void gnss_relative_c::ave_SNR() {
    int numall=0; SNR[NSYS]=0.0;
    for (int sys=0; sys<NSYS; sys++) {
        SNR[sys]=0.0;
        int codenum=nobs[sys][numF]+1,sysobs=0;
        if (irfsat[sys]<0||codenum<MINSYSOBS_RTK) continue;
        for (int isat = 0; isat<satnum; isat++) {
            if (obsr->data[rovsat[isat]].isys!=sys) continue;
            SNR[sys] += ssat[usesat[isat]-1].snr[0];
            sysobs++; numall++;
        }
        if (sysobs>0) {
            SNR[NSYS]+=SNR[sys];
            SNR[sys]/=1.0*sysobs;
        }
    }
    if (numall>0) SNR[NSYS]/=numall;
}
/* average multipath of different systems ----------------------------------------- */
void gnss_relative_c::ave_multiPath() {
    int numall=0; multiPath[NSYS]=0.0;
    for (int sys=0; sys<NSYS; sys++) {
        multiPath[sys]=0.0;
        int codenum=nobs[sys][numF]+1,sysobs=0;
        if (irfsat[sys]<0||codenum<MINSYSOBS_RTK) continue;
        for (int isat = 0; isat<satnum; isat++) {
            if (obsr->data[rovsat[isat]].isys!=sys) continue;
            multiPath[sys] += ssat[usesat[isat]-1].multiPath[0];
            sysobs++; numall++;
        }
        if (sysobs>0) {
            multiPath[NSYS]+=multiPath[sys];
            multiPath[sys]/=1.0*sysobs;
        }
    }
    if (numall>0) multiPath[NSYS]/=numall;
}
/* write state information to log_stream ------------------------------------------ */
void gnss_relative_c::write_state() {
    if (!log_stream.is_open()) return;

    if (opt->mode>=PMODE_KINEMA && opt->mode<=PMODE_MOVEB) {
        /* time */
        obsr->data[0].time.time2str(3);
        log_stream << obsr->data[0].time.sep;

        /* error message */
        if ( stat==SOLQ_SINGLE || stat==SOLQ_SBAS ) {
            log_stream << "   ";
            if (msg.size()>0) log_stream << msg << ",   ";
            /* standard variance of unit weight in LSQ */
            log_stream << setprecision(3);
            log_stream << "SPP LSQ Sigma0: " << SQRT(adjfunc->sigma0_2);
            log_stream << "\n";
            for (int i=0; i<obsr->n&&i<MAXOBS; i++) {
                if ( obsr->data[i].errorMsg.size()>0 && obsr->data[i].errorMsg.size()<MAXSTRMSG ) 
                    log_stream << "    " << ssat[obsr->data[i].sat-1].id << ": "
                    << obsr->data[i].errorMsg << "\n";
            }
        }
        /* state message */
        else {
            /* set float format */
            log_stream.setf(ios::fixed);
            log_stream << setprecision(4);

            /* reference satellite (in the same line of time) */
            log_stream << " Reference Sat: ";
            for (int sys = 0; sys<NSYS; sys++)
                if (irfsat[sys]>=0) log_stream << " " << ssat[rfsat[1][sys] - 1].id;
            /* adjustment time (in the same line of time) */
            /* if (opt->adjustfunc!=ADJUST_HELMERT||opt->adjustfunc==ADJUST_HELMERT&&adjfunc->H_niter)*/
            log_stream << "   Filter Time:" << setw(5) << dtime;
            /* solution number (in the same line of time) */
            log_stream << "   Solution Number:" << setw(7) << nsol;
            /* fixed ratio (in the same line of time) */
            if (opt->modear!=ARMODE_OFF) {
                if (opt->modear==ARMODE_LCWN) log_stream << "   Fix Rate: " << setw(8) << fixNumRatio;
                else log_stream << "   Fix Var Ratio: " << setw(8) << fixVarRatio;
            }
            log_stream << "\n";

            /* BODY of state ------------------------- */
            /* tide correction */
            if (opt->tidecorr>0) {
                log_stream << "    TideCorrect:    Rov: ";
                for (int i=0; i<3; i++)
                    log_stream  << setw(9) << Rtide[i];
                log_stream << "   Ref: ";
                for (int i=0; i<3; i++)
                    log_stream  << setw(9) << Btide[i];
                log_stream << "\n";
            }
            /* ionosphere parameters */
            /* troposphere parameters */
            if (NT>0) {
                log_stream << "    Troposphere: ";
                for (int i=0; i<NT; i++) log_stream << "   "  << TRO_PAR[i]
                    << setw(9) << *(sPar + iT+i) << setw(9) << SQRT(*(sRx + (iT+i)*numX+(iT+i)));
                log_stream << "\n";
            }
            /* GLONASS differenced IFB rate paramters */
            if (NG>0) {
                log_stream << "    GLONASS IFB: ";
                for (int i=0; i<NG; i++) log_stream << setw(10) << *(sPar + iG+i);
                log_stream << "\n";
            }
            /* satellites number of different system */
            if (opt->logmsg&GNSSLOG_NSAT) {
                log_stream << "    System nSat: ";
                for (int sys=0; sys<NSYS; sys++) {
                    int nsyso = nobs[sys][numF]>0 ? nobs[sys][numF]+1 : 0;
                    log_stream << "   " << SYS_LIST[sys] << ":" << setw(8) << nsyso;
                }
                log_stream << "   All:" << setw(8) << satnum;
                log_stream << "\n";
            }
            /* PDOP value of different satellite systems */
            if (opt->logmsg&GNSSLOG_PDOP) {
                /* calculate pdop of different satellite system */
                gnss_relative_c::get_pdop();
                log_stream << "    PDOP Values: ";
                for (int sys=0; sys<NSYS; sys++)
                    log_stream << "   " << SYS_LIST[sys] << ":" << setw(8) << Pdop[sys];
                log_stream << "   All:" << setw(8) << Pdop[NSYS];
                log_stream << "\n";
            }
            /* average SNR of different satellite system */
            if (opt->logmsg&GNSSLOG_SNR) {
                ave_SNR();
                log_stream << "    Average SNR: ";
                for (int sys=0; sys<NSYS; sys++) {
                    log_stream << "   " << SYS_LIST[sys] << ":" << setw(8) << SNR[sys];
                }
                log_stream << "   All:" << setw(8) << SNR[NSYS];
                log_stream << "\n";
            }
            /* average multiPath of different satellite system */
            if (opt->logmsg&GNSSLOG_MPATH) {
                ave_multiPath();
                log_stream << "    Multi--Path: ";
                for (int sys=0; sys<NSYS; sys++) {
                    log_stream << "   " << SYS_LIST[sys] << ":" << setw(8) << multiPath[sys];
                }
                log_stream << "   All:" << setw(8) << multiPath[NSYS];
                log_stream << "\n";
            }
            /* helmert component covariance sigma */
            if (opt->adjustfunc==ADJUST_HELMERT&&opt->logmsg&GNSSLOG_HELMERT) {
                log_stream << "    Helmert Sgm: ";
                for (int i=0; i<adjfunc->nsys+1; i++)
                    log_stream << "   " << adjfunc->sys_flag[i] << ":" << setw(8) << adjfunc->sgm2[i];
                log_stream << "   Num:" << adjfunc->H_niter;
                log_stream << "\n";
            }
            /* xpar */
            if (opt->logmsg&GNSSLOG_XPAR) {
                log_stream << "    X-parameter: ";
                for (int i=0; i<numX; i++) log_stream << "   " << *(sPar + i);
                log_stream << "\n";
            }
            /* dx */
            if (opt->logmsg&GNSSLOG_DXPAR) {
                log_stream << "    float-dXpar: ";
                for (int i=0; i<numX; i++) log_stream << "   " << Xpar[i]-Xori[i];
                log_stream << "\n";
            }
            /* satellite parameters (prn,az,el,solv_con,ion,amb) */
            for (int i=0; i<satnum; i++) {
                rp = &obsr->data[rovsat[i]]; bp = &obsb->data[bassat[i]];
                gnss_ssat_c *sss = &ssat[ usesat[i]-1 ];
                log_stream << "    " << sss->id << ":"; //prn
                lam=nav->lam[rp->sat-1];
                /*rover az el */
                log_stream << " A" << setw(9) << rp->azel[0]*R2D << " E" << setw(8) << rp->azel[1]*R2D;
                /* antenna correction */
                log_stream << " DANT:" << setw(9) << rp->rant[0] << setw(9) << rp->rant[1];
                /* x y z dt and transmisson time*/
                log_stream << " XYZT:";
                log_stream << setw(16) << rp->posvel[0];
                log_stream << setw(16) << rp->posvel[1];
                log_stream << setw(16) << rp->posvel[2];
                log_stream << setw(12) << rp->dts[0]*1E6;
                log_stream << setw(9) << rp->time.timediff(obsr->data[i].sigtime);
                /* troposphere */
                if (NT > 0) log_stream << " TRO:" << setw(9) << *(sPar + iT) * rp->m_w;
                else log_stream << " TRO:" << setw(9) << (opt->tropopt==TROPOPT_OFF? 0.0 : rp->dtro);
                /* ionosphere delay of L1 (m) */
                if (nI>0) log_stream << " ION:" << setw(10) << *(sPar + iI+i)
                    << setw(9) << *(sRx + iI+i+(iI+i)*numX);
                else log_stream << " ION:" << setw(10) << (opt->ionoopt==IONOOPT_OFF||opt->ionoopt==IONOOPT_IFLC? 0.0 :
                    SQR(lam[0]/GNSS_LAMBDA_M[iGPS][1])*rp->dion*rp->ionmap);

                /* LC ambiguity */
                if (opt->ionoopt==IONOOPT_IFLC) {
                    int namb = ambnum[i];
                    log_stream << " LC:";
                    log_stream << setw(6) << sss->solv_con[NFREQ];//solv_con
                    //log_stream << setw(12) << Xori[iA+namb]; //orignal ambiguity
                    if (namb>=0) {
                        log_stream << setw(16) << sss->amb[NFREQ]; // new ambiguity
                        log_stream << setw(10) << SQRT(*(sRx + (iA+namb)*numX+(iA+namb)));
                    }
                    else log_stream << setw(16) << 0.0 << setw(10) << 0.0; //amb
                    //log_stream << setw(10) << sss->d_ave[NFREQ+1]; // d_mw12
                    //log_stream << setw(10) << sss->amb_ave[NFREQ]; // d_N1
                    log_stream << " CS:" <<  
                        setw(3) << sss->slip_f[0]<< setw(3) << sss->slip_f[1]; //L1 L2 cycle slip flag

                    //measurement code
                    log_stream << " MEAS: PC-" << setw(3) <<  rp->PCode[0] << "-" << setw(3) <<  rp->PCode[1];
                    log_stream << setw(10) << sss->resP[0];
                    log_stream << " LC-" << setw(3) <<  rp->LCode[0] << "-" << setw(3) <<  rp->LCode[1];
                    log_stream << setw(10) << sss->resL[0];
                }
                /* uncombined ambiguity */
                else {
                    for (int freq = 0; freq<numF; freq++) {
                        int namb = ambnum[i + freq*satnum];
                        log_stream << " L" << freq + 1 << ":";
                        log_stream << setw(6) << sss->solv_con[freq]; //solv_con
                        //log_stream << setw(12) << Xori[iA+namb]; //orignal ambiguity
                        if (namb>=0) {
                            log_stream << setw(16) << sss->amb[freq]; // new ambiguity
                            log_stream << setw(10) << SQRT(*(sRx + (iA+namb)*numX+(iA+namb)));
                        }
                        else log_stream << setw(16) << 0.0 << setw(10) << 0.0; //amb
                        //log_stream << setw(12) << sss->amb_ave[freq]; // d_N1
                        log_stream << " CS:" <<
                            setw(3) << sss->slip_f[freq]; //Ln cycle slip flag

                        //measurement code
                        log_stream << " MEAS: P-" << setw(3) <<  rp->PCode[freq];
                        log_stream << setw(10) << sss->resP[freq];
                        log_stream << " L-" << setw(3) <<  rp->LCode[freq];
                        log_stream << setw(10) << sss->resL[freq];
                    }
                }
                /* warning message */
                log_stream << "  " << rp->warningMsg << "\n";
            }
            /* output excluded satellite messages */
            if (opt->logmsg&GNSSLOG_EXCSAT) {
                for (int i=0; i<obsr->n; i++) {
                    if ( !obsr->data[i].used && obsr->data[i].errorMsg.size()>0 && obsr->data[i].errorMsg.size()<MAXSTRMSG ) 
                        log_stream << "    " << ssat[obsr->data[i].sat-1].id << ": " << obsr->data[i].errorMsg << "\n";
                }
            }
        }
    }
}
/* update satellite sate vector (ssat) -------------------------------------------- */
void gnss_relative_c::update_ssat() {
    for (int i=0; i<MAXSAT; i++) {
        ssat[i].vs=0;
    }
    gnss_ssat_c *sss; int namb;
    for (int isat = 0; isat<satnum; isat++) {
        lam = nav->lam[usesat[isat]-1];
        rp = &obsr->data[rovsat[isat]]; bp = &obsb->data[bassat[isat]];
        sss = &ssat[usesat[isat]-1];
        sss->vs = 1;
        sss->azel[0] = rp->azel[0];  sss->azel[1] = rp->azel[1];
        for (int freq = 0; freq<numF; freq++) {
            sss->snr[freq] = (rp->SNR[freq]+bp->SNR[freq])/2.0;
            /* ionosphere solution */
            if (nI>0) {
                sss->ion_delay=*(sPar + iI+isat);
                sss->ion_var=*(sRx + iI+isat + (iI+isat)*numX);
            }
            /* ambiguity solution and multipath */
            if ((namb = ambnum[isat + freq*satnum])>=0) {
                // LC ambiguity
                if (opt->ionoopt==IONOOPT_IFLC) {
                    // ambiguity
                    sss->amb[NFREQ] = *(sPar+iA+namb);
                    sss->ambvar[NFREQ] = *(sRx + iA+namb + (iA+namb)*numX);
                    if (sss->ave_con[NFREQ]>=iniamb) sss->fix[NFREQ]=2;
                    // update solution indexes
                    sss->solv_con[NFREQ]++;
                    sss->solv_period[NFREQ]+=sss->obstime.back().timediff(sss->solv_time[NFREQ]);
                    sss->solv_time[NFREQ]=sss->obstime.back();
                    // multipath
                    sss->multiPath[NFREQ]=fabs(sss->pc12[1]-sss->lc12[1]+sss->amb[NFREQ]);
                }
                // uncombined ambiguity
                else {
                    // ambiguity
                    sss->amb[freq] = *(sPar+iA+namb)/*Xpar[iA+namb]*/;
                    sss->ambvar[freq] = *(sRx + iA+namb + (iA+namb)*numX);
                    /*if (sss->ambvar[freq]<=0) iniamb=iniamb;*/
                    if ( stat==SOLQ_FIX && opt->modear==ARMODE_FIXHOLD ) {
                        /*if ( (nreset[freq]||sss->fix[freq]!=3) && 
                            rp->azel[1]>opt->minElArHold && bp->azel[1]>opt->minElArHold ) {
                            sss->ambvar[freq]=FIXED_AMB_VAR;
                            for (int i=0; i<numX; i++) //reset covariance matrix
                                *(sRx + i + (iA+namb)*numX) = *(sRx + (iA+namb) + i*numX) = 0.0;
                            sss->fix[freq]=3;
                        }
                        else sss->fix[freq]=2;*/
                        sss->fix[freq]=3;
                    } else if ((sss->ave_con[freq]<iniamb)) {
                        sss->fix[freq]=1;
                    }
                    else sss->fix[freq]=2;
                    // update solution indexes
                    sss->solv_con[freq]++;
                    sss->solv_period[freq]+=sss->obstime.back().timediff(sss->solv_time[freq]);
                    sss->solv_time[freq]=sss->obstime.back();
                    // multipath
                    sss->multiPath[freq]=
                        fabs(sss->P[freq].back()-lam[freq]*(sss->L[freq].back()+sss->amb[freq]));
                }
            }
            /* residual */
            if (iobs[freq][isat]>=0) sss->resL[freq]=adjfunc->V[ iobs[freq][isat] ];
            else sss->resL[freq]=999;
            if (iobs[freq+numF][isat]>=0) sss->resP[freq-numF] = adjfunc->V[ iobs[freq+numF][isat] ];
            else sss->resP[freq]=999;
        }
    }
}
/* update solution vector (sol) --------------------------------------------------- */
void gnss_relative_c::update_sol() {

    soltime[2]=solp->time;
    solp->type = 0;
    solp->NL = numL;
    solp->nsol=nsol;

    /* parameters solution */
    // dynamic parameters
    for (int i=0; i<NP; i++) {
        solp->xpos[i] = *(sPar+i);
        for (int j = 0; j<NP; j++)
            solp->vpos[j + i*3] = *(sRx + j+i*numX);
    }
    // ionosphere parameters
    solp->xion.clear(); solp->vion.clear();
    solp->NI=nI;
    for (int i=0; i<nI; i++) {
        solp->xion.push_back(*(sPar+iI+i));
        for (int j = 0; j<nI; j++)
            solp->vion.push_back(*(sRx + j+iI + (i+iI)*numX));
    }
    // troposphere parameters
    for (int i=0; i<NT; i++) {
        solp->xtro[i] = *(sPar+iT+i);
        for (int j = 0; j<NT; j++)
            solp->vtro[j + i*NT] = *(sRx + j+iT + (i+iT)*numX);
    }
    // GLO freq-diffrerence paramters
    for (int i=0; i<NG; i++) {
        solp->xglo[i] = *(sPar+iG+i);
        for (int j = 0; j<NG; j++)
            solp->vglo[j + i*NG] = *(sRx + j+iG + (i+iG)*numX);
    }
    // ambiguity parameters
    solp->NA = nA;
    solp->xamb.assign(sPar+iA,sPar+iA+nA);
    solp->vamb.assign(nA*nA,0.0);
    for (int i=0; i<nA; i++) for (int j = 0; j<nA; j++)
        solp->vamb[j + i*nA] = *(sRx + j+iA + (i+iA)*numX);

    solp->ns = (unsigned int)satnum;
    if (stat==SOLQ_FIX && opt->modear!=ARMODE_LCWN) solp->ratio=float(fixVarRatio);
    solp->stat = stat;

    /* update Rx_ALL */
    // set Rx_ALL from sRx
    for (int i=0; i<numX; i++) {
        X_ALL[ Xest[i] ] = Xpar[i];
        Rx_ALL[ Xest[i] + Xest[i]*N_ALL] = *(sRx+i+i*numX); //variance from Rxvec
        for (int j=0; j<i; j++) {
            Rx_ALL[ Xest[j] + Xest[i]*N_ALL ] = Rx_ALL[ Xest[i] + Xest[j]*N_ALL ] = *(sRx+j+i*numX);
        }
    }
}

/* preprocess functions ----------------------------------------------------------- */
/* preprocess gnss observation ---------------------------------------------------- */
int gnss_relative_c::preprocess_gnss_obs(const double pos[3],const double vel[3]) {
    warningMsg="";

    init_sol(1);

    /* update consecutive number for each satellite */
    update_sat_nCon();

    /* initialize single position */
    gnss_ini_EPH_DCB();

    /* detect wrong observations with double frequency */
    detect_wrong_obs_2freq();

    /* observation and ambiguity information for each satellites */
    update_obs_amb_ssat();

    /* intiialize position and velocity */
    gnss_get_posvel(pos,vel,NULL,NULL);
    /* assign to Xpar */
    Xpar.assign(3,0);
    vXpar.assign(DPLNX,0);
    for (int i=0; i<3; i++) {
        Xpar[i]=ext_posvel[i];
        vXpar[i]=ext_posvel[i+3];
    }

    /* 1st time to pre-process observations */
    if (code_matrix_prePos()<numX) {
        /* verify if lack of valid observation */
        obsp->warnMsg="lack of valid satellites for SPP!";
    }
    /* exculde observations with large residual */
    exc_largeres(0);
    /* exclude constellations with only a few available satellites */
    exc_few_sat();

    /* 2nd time to pre-process observations */
    if (code_matrix_prePos()<numX) {
        /* verify if lack of valid observation */
        obsp->warnMsg="lack of valid satellites for SPP!";
    }

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

    /* pre-process clock */
    if (preClkFlag) process_pre_clock();
    else {
        preClkVar[0]=VAR_CLK_OFF2;
        preClkVar[1]=VAR_CLK_DRF2;
        /* update vRx */
        if ( opt->dynamics ) {
            vRx[15]=preClkVar[1]; //clk offset noise
        }
    }
    gnss_single_c::update_model_correction();

    return 1;
}
/* get gnss observation equation -------------------------------------------------- */
int gnss_relative_c::gnss_obs_equation() {

    /* carrier phase bias correction */
    if (opt->pppopt.find("-DIS_FCB") != string::npos) corr_phase_bias();

    /* tt : time difference between current and previous (s) */
    tt = soltime[2].time!=0 ? obsr->data[0].time.timediff(soltime[2]) : 0.0;

    if (timediff_rb()==0)
        return 0;

    /* select common satellite and test availability */
    if (selcomsat_phase()<MINSAT) {
        msg = "no enough common satellite for DGNSS!";
        return 0;
    }

    /* reset satellites status */
    resetsat();

    /* test */
    if (nsol==iniamb-1)
        opt->niter=opt->niter;

    /* update parameters */
    updatepar(1);

    /* get double-differenced observation equation */
    double_diff();
    Xori=Xpar;

    return 1;
}
/* update gnss status ------------------------------------------------------------- */
int gnss_relative_c::gnss_update_status() {

    /* point to float solution */
    sPar=Xpar.begin(); sRx=Rx.begin();

    nsol++; ntotal++;

    /* update satellite sate vector (ssat) */
    gnss_relative_c::update_ssat();
    /* update solution vector (sol) */
    gnss_relative_c::update_sol();

    stat=SOLQ_FLOAT;
    /* write state information to log_stream */
    gnss_relative_c::write_state();

    return stat;
}
/* relative position function ----------------------------------------------------- */
int gnss_relative_c::gnss_pos() {
    warningMsg="";

    /* carrier phase bias correction */
    if (opt->pppopt.find("-DIS_FCB")!=string::npos) corr_phase_bias();

    /* base position has already been set */
    /* if no SPP solution for rover yet, run SPP for rover */
    if (fabs(sol.back().obsTime.timediff(obsr->data[0].time))>1E-5) {
        init_sol(1);
        if (single()==SOLQ_NONE) {
            return 0;
        }
    }

    /* base observation? */
    if (obsb->n <= 0) {
        msg = "no base-station observations!";
        gnss_relative_c::write_state();
        return 0;
    }

    /* tt : time difference between current and previous (s) */
    tt = soltime[2].time!=0 ? obsr->data[0].time.timediff(soltime[2]) : 0.0;
    /* time difference between rover and base (s) */
    if (timediff_rb()==0) {
        return 0;
        gnss_relative_c::write_state();
    }

    /* select common satellite and test availability */
    if (selcomsat_phase()<MINSAT) {
        msg = "no enough common satellite for DGNSS!";
        gnss_relative_c::write_state();
        return 0;
    }

    /* reset satellites status */
    resetsat();

    /* test */
    if (nsol==iniamb-1)
        opt->niter=opt->niter;

    /* update parameters */
    updatepar(0);

    /* compute float solution */
    for (int i=0; i<opt->niter; i++) {
        double_diff();
        Xori=Xpar;

        dtime=0; unsigned int stime=tickget();
        /* adjustment with least square function */
        if ((adjfunc->adjustment(Acoe,Lobs,Rvar,Xpar,Rx,numL,numX,numF,nobs))==-1) {
            msg = "filter error!";
            break;
        }
        dtime=(int)(tickget()-stime);
        sum_var.assign(2,0.0);
        stat=SOLQ_FLOAT;
        nsol++; ntotal++;
    }

    /* fix solution */
    if ( stat==SOLQ_FLOAT && opt->ionoopt!=IONOOPT_IFLC && opt->mode!=PMODE_DGPS ) {
        /* uncombined ambiguity */
        if ( opt->modear!=ARMODE_OFF && get_fixsolution_el() ) stat=SOLQ_FIX;
        /* LC ambiguity */
        else if ( nsol>iniamb && opt->modear==ARMODE_LCWN ) {
            stat=SOLQ_FIX;
            fixNumRatio=1.0*fix_num/(1.0*n_Damb);
        }
        else
            stat=SOLQ_FLOAT;
    }

    /* update solution */
    if ( stat==SOLQ_FLOAT || stat==SOLQ_FIX ) {
        if (opt->mode==PMODE_DGPS) stat = SOLQ_DGPS;
        /* final solution iterator */
        if (stat==SOLQ_FIX && opt->modear!=ARMODE_LCWN) {
            sPar=fix_Xpar.begin(); sRx=fix_Rx.begin();
        }
        else {
            sPar=Xpar.begin(); sRx=Rx.begin();
        }

        /* update reference satellite to irfsat */
        for (int i=0; i<NSYS; i++) rfsat[0][i]=rfsat[1][i];

        /* update satellite state */
        gnss_relative_c::update_ssat();
        /* update solution */
        gnss_relative_c::update_sol();

        /* update baseline status */
        update_baseline();
    }

    /* write state information to log_stream */
    gnss_relative_c::write_state();

    return 1;
}
