#include "GNSS/PosModel/position.h"
#include "BaseFunction/basefunction.h"

/* constant --------------------------------------------------------------------------------------- */
#define MINSAT              4               /* minimum satellites number */
#define SQR(x)              ((x)*(x))
#define SQRT(x)             ((x)<0.0?-sqrt(-x):sqrt(x))

#define INIT_GRA	        1E-5			/* initial tro gradient */
#define VAR_GRA             2.5E-5          /* initial variance of gradient (m^2) */
#define INIT_HWBIAS	        0               /* initial glo frequency-dependent amb bias */
#define VAR_HWBIAS          0.01            /* initial variance of h/w bias ((m/MHz)^2) */
#define RAT_HWBIAS	        1E-7            /* growth rate of std h/w bias (m/MHz/sqrt(s)) */

#define VAR_POS             3600.0          /* initial variance of receiver pos (m^2) */
#define VAR_CLK             3600.0          /* initial variance of receiver clock (m^2) */
#define VAR_VEL             3600.0          /* initial variance of receiver vel ((m/s)^2) */

#define VAR_CLK_OFF1        3600.0          /* initial variance of pre-estimated receiver clock offset (m^2) */
#define VAR_CLK_DRF1        3600.0          /* initial variance of pre-estimated receiver clock drift (m^2/s^2) */
#define VAR_CLK_OFF2        9E12	        /* initial variance of receiver clock offset (m^2) */
#define VAR_CLK_DRF2        1E6             /* initial variance of receiver clock drift (m^2/s^2) */

#define TTOL_MOVEB  (1.0+2*DTTOL)

const double freq_glos[]={ FREQ1_GLO,FREQ2_GLO };
const string SYS_LIST="GREC";               /* system flag GPS, GLONASS, Galileo, BeiDou */
const string TRO_PAR[]={ "ZTD1:","G_N1:","G_E1:" };

/* precise point position class ----------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
gnss_ppp_c::gnss_ppp_c() {
    /* station geographic parameters (ini. in  updatepar()) */
    for (int i=0; i<3; i++)
        Rxyz[i] = Rblh[i] = Rtide[i] = 0.0;
    tidevar=0;

    /* centre satellite parameters */
    for (int i=0; i<NSYS; i++) {
        S_nsat[i] = 0;
        clk_sys[i][0]=clk_sys[i][1]=-1;
        irfsat[i] = -1;
        rfsat[1][i] = rfsat[0][i] = 0;
        dist1[i] = dtro1[i] = dion1[i] = 0.0;
        clock1[i]=0;
        for (int j = 0; j<NFREQ; j++) {
            L_nsat[i][j] = nobs[i][j] = 0; 
            ambN1[i][0]=0;
            Lobs1[i][j]=Lobs1[i][NFREQ+j]=0;
            varRB1[i][j]=varRB1[i][NFREQ+j]=0;
        }
    }

    for (int i=0; i<NFREQ+1; i++) nreset[i]=reset_flag[i]=0;

    preClkVar[0]=VAR_CLK;

    satnum=numsys=0;
    iC = iT = iU = iI = iA = 0;
    nC = nU = nI = nA = 0;
    iniamb=60;
    fixVarRatio=0;
    fix_num = 0;
}
gnss_ppp_c::~gnss_ppp_c() {
    rp=NULL; bp=NULL;
    lam=lam1=lam2=NULL;

    usesat.clear(); rovsat.clear();
    for (int i=0; i<NFREQ; i++) {
        CJ_n100[i].clear();
    }
    Rxvec.clear();

    Acoep=NULL;
    for (int sys = 0; sys<NSYS; sys++) Airfsat[sys].clear();
    Asatnum.clear();
    Ais.clear();

    ambnum.clear();
    fix_flag.clear();
    fix_amb.clear(); sum_var.clear();
    fix_Xpar.clear(); fix_Rx.clear();
}

/* Implementation functions ----------------------------------------------------------------------- */
/* correct and save pre-etimated clock offset and drift --------------------------- */
void gnss_ppp_c::process_pre_clock() {
    /* clock offset */
    preClkVar[0]=VAR_CLK_OFF1;
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
int gnss_ppp_c::test_system(int Sysyem,int System_Num) {
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
/* select available satellites of rover ------------------------------------------- */
int gnss_ppp_c::selavasat() {
    for (int sys = 0; sys<NSYS; sys++) {
        S_nsat[sys]=0;
        clk_sys[sys][0]=clk_sys[sys][1];
        clk_sys[sys][1]=-1;
        for (int freq = 0; freq<NFREQ; freq++)
            L_nsat[sys][freq] = 0;
    }
    nC=0;
    satnum = 0; usesat.clear(); rovsat.clear();

    /* select common satellites */
    for (int i=0; i<obsr->n; i++) {
        if ( obsr->data[i].used && ncode[0][obsr->data[i].isys]>=MINSYSOBS ) {
            /* add satellite */
            usesat.push_back(obsr->data[i].sat);
            rovsat.push_back(i);
            S_nsat[obsr->data[i].isys]++;
            /* get phase observations number */
            int L_ava[NFREQ] ={ 0 };
            rp = &obsr->data[i];
            for (int freq = 0; freq<NF; freq++) { //loop frequency
                if ( rp->L[freq]!=0.0 && rp->P[freq]!=0.0 ) L_ava[freq] = 1;
            }
            if (opt->ionoopt==IONOOPT_IFLC && L_ava[0] && L_ava[1]) L_nsat[obsr->data[i].isys][0]++;
            if (opt->ionoopt!=IONOOPT_IFLC) {
                for (int freq = 0; freq<NF; freq++) if (L_ava[freq]) L_nsat[obsr->data[i].isys][freq]++;
            }
        }
    }
    satnum = usesat.size();

    if (satnum <= 0) return 0;

    numsys = 0;
    /* select reference satellite according to phase observation number and elevation */
    if (opt->pppmode==PPP_SINGLE_DIFF) {
        /* initialize reference satellites parameters */
        init_refsat();
        /* get phase observation number of common satellite */
        vector<int> nphase_sat(satnum,0);
        for (int sys = 0; sys<NSYS; sys++) { // loop NSYS systems
            irfsat[sys]=-1;
            for (int i=0; i<satnum; i++) { //loop every satellite
                rp = &obsr->data[rovsat[i]];
                if (rp->isys!=sys) continue; //test satellite system
                for (int freq = 0; freq<NF; freq++) { //loop frequency
                    if ( rp->L[freq]!=0.0 && rp->P[freq]!=0.0 ) nphase_sat[i]++;
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
    else {
        /* get clock information */
        for (int i=0; i<NSYS; i++) {
            if (S_nsat[i]>=MINSYSOBS) {
                numsys++;
                clk_sys[i][1]=nC++;
                if (numsys==1) {
                    clkSys1st=i;
                    if (ntotal<iniamb) clkRefSys=clkSys1st;
                }
            }
            /* estimate pseudo clock offset for reference system that has no observations in this epoch */
            //else if ( clkRefSys==i && ntotal>=iniamb ) {
            else if ( clk_sys[i][0]>=0 && ntotal>=iniamb ) {
                clk_sys[i][1]=nC++;
            }
        }
    }

    return satnum;
}

/* initialize reference satellites parameters ------------------------------------- */
void gnss_ppp_c::init_refsat() {
    for (int i=0; i<NSYS; i++) {
        irfsat[i] = -1;
        rfsat[1][i] = 0;
        dist1[i] = dtro1[i] = dion1[i] = 0.0;
        clock1[i]=0;
        for (int j=0; j<NFREQ; j++) { 
            ambN1[i][0]=0;
            Lobs1[i][j]=Lobs1[i][NFREQ+j]=0;
            varRB1[i][j]=varRB1[i][NFREQ+j]=0;
        }
    }
    gloIFB1[0] = gloIFB1[1] = 0.0;
}
/* initialize vector and matrix according common satellite ------------------------ */
void gnss_ppp_c::init_arrmat(int prePos) {
    Lobs.clear(); Acoe.clear();
    Rvec.clear(); Rvar.clear(); Rx.clear();
    QAA.clear();
    fix_flag.clear();
    fix_amb.clear(); sum_var.clear();
    fix_Xpar.clear(); fix_Rx.clear();
    ambnum.assign(numF*satnum,-1);

    /* number of kinds of parameters */
    /* NP NT nI nU nC nA */
    if (opt->ionoopt==IONOOPT_CONST) nI=satnum;
    else nI=0;
    // if no observations for reference system, use pseudo clock for reference system 
    nC = nC; // initialized in function selavasat()
    nU = nA = 0;
    numX = NP+NT+NG+nI+nC+nA;

    /* set start index of kinds of parameters */
    iT=NP; iG=iT+NT; iI=iG+NG; iC=iI+nI; iA=iC+nC;

    /* initialize parameter and its variance vector (no amb) */
    if (prePos==1) Xpar.insert(Xpar.end(),(numX-3),0.0); // if have pre-estimated position
    else Xpar.assign(numX,0.0);
    Rxvec.assign(numX,0.0);

    /* intialize parameter information vector */
    Xini.assign(numX,0); Xest.assign(numX,-1);
    for (int i=0; i<iI; i++) Xest[i]=i;

    /* ambiguity fixed rate */
    if (opt->modear==ARMODE_FIXHOLD||opt->modear==ARMODE_LCWN) fixVarRatio=fixNumRatio=0;
}

/* update parameters functions ---------------------------------------------------- */
/* update ambiguity parameters ---------------------------------------------------- */
void gnss_ppp_c::updateamb() {
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
            if (S_nsat[sys]<MINSYSOBS) continue; //continue if no observation of this system
            for (int i=0; i<satnum; i++) {
                if (obsr->data[rovsat[i]].isys!=sys) continue;
                sss = &ssat[usesat[i]-1];

                /* update ambiguity */
                lam=nav->lam[ sss->sat - 1 ];
                rp=&obsr->data[rovsat[i]];
                /* correct antenna and phase windup for observations and ambiguities */
                sss->correct_obs_amb(NFREQ,rp,NULL);
                /* reset ambiguity if this frequency should be reset */
                if (reset_flag[NFREQ]&&sss->reset[NFREQ]<1) sss->reset_amb(NFREQ);
                /* update ambiguity parameters */
                if (sss->update_solved_amb(NFREQ,nsol)>0) {
                    Xpar.push_back(sss->amb[NFREQ]);
                    Rxvec.push_back(sss->ambvar[NFREQ]);
                    ambnum[i] = nA++;
                    //Xest index in all parameters
                    Xest.push_back(NX+(usesat[i]-1)*numF);
                    /* if fix==1, initialize covariance matrix of all parameters */
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
            if (S_nsat[sys]<MINSYSOBS) continue; //continue if no observation of this system
            for (int freq = 0; freq<numF; freq++) {
                for (int i=0; i<satnum; i++) {
                    if (obsr->data[rovsat[i]].isys!=sys) continue;
                    sss = &ssat[usesat[i] - 1];

                    /* update ambiguity */
                    lam=nav->lam[ sss->sat - 1 ];
                    rp=&obsr->data[rovsat[i]];
                    /* correct antenna and phase windup for observations and ambiguities */
                    sss->correct_obs_amb(freq,rp,NULL);
                    /* reset ambiguity if this frequency should be reset */
                    if (reset_flag[freq]&&sss->reset[freq]<1) sss->reset_amb(freq);
                    /* update ambiguity parameters */
                    if (sss->update_solved_amb(freq,nsol)>0) {
                        Xpar.push_back(sss->amb[freq]);
                        Rxvec.push_back(sss->ambvar[freq]);
                        /* add a noise if detected reseted ambiguity */
                        if (nreset[freq]&&opt->modear==ARMODE_FIXHOLD) Rxvec.back()=SQR(0.3);
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
void gnss_ppp_c::updatexyz(int prePos) {
    /* tightly integration */
    if ( prePos==1 ) {
        return; /* the geodetic parameters was updated in "gnss_get_posvel()" */
    }

    if ( varFlag&1 ) {
        for (int i=0; i<3; i++) {
            Xpar[i] = ext_posvel[i];
            Rxvec[i] = 3;
            Xini[i]=1;
            init_RxALL(i);
        }
        varFlag=0;
        return;
    }

    /* fixed mode */
    if (opt->mode==PMODE_PPP_FIXED) {
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
    if (opt->mode==PMODE_PPP_STATIC) {
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
    if (opt->mode==PMODE_PPP_KINEMA) {
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
void gnss_ppp_c::updatetro() {
    /* initialize tro parameters for the first epoch */
    //if ( sol[MAXSOLBUF-2].xtro[0]==0.0 || ntotal<iniamb ) {
    if ( X_ALL[iT]==0.0 || ntotal<iniamb ) {
        const double zazel[] ={ 0.0,PI/2.0 };
        gnss_obsd_c obs;
        obs.time=sol.back().time;
        trofunc->saascorr(&obs,Rblh,zazel,0.7);
        Xpar[iT] = obs.dtro;
        //Rxvec[iT]= obs.trovar;
        Rxvec[iT]= SQR(opt->std[2]);
        Xini[iT]=1;
        init_RxALL(iT);
        /* estimate tro gradient */
        if (NT==3) for (int j=1; j<3; j++) {
            Xpar[iT+j] =INIT_GRA;
            Rxvec[iT+j]=VAR_GRA;
            Xini[iT+j]=1;
            init_RxALL(iT+j);
        }
    }
    /* update tro parameters using last solution */
    else {
        //Xpar[iT] =sol[MAXSOLBUF-2].xtro[0];
        //Rxvec[iT]=sol[MAXSOLBUF-2].vtro[0]+opt->stdrate[2]*tt;
        Xpar[iT]  = X_ALL[iT];
        Rxvec[iT] = Rx_ALL[ iT + iT*N_ALL ] + opt->stdrate[2]*tt;
        /* estimate tro gradient */
        if (NT==3) for (int j=1; j<3; j++) {
            //Xpar[iT+j] =sol[MAXSOLBUF-2].xtro[j];
            //Rxvec[iT+j]=sol[MAXSOLBUF-2].vtro[j+j*NT] + opt->stdrate[2]*tt*0.01;
            Xpar[iT+j]  = X_ALL[iT+j];
            Rxvec[iT+j] = Rx_ALL[ (iT+j)*(N_ALL+1) ] + opt->stdrate[2]*tt*0.01;
        }
    }
}
/* update receiver GLONASS IFB parameters --------------------------------- */
void gnss_ppp_c::updateglo() {
    for (int i=0; i<NG; i++) {
        /* initialize receiver GLONASS IFB parameter */
        Xpar[iG+i] = INIT_HWBIAS;
        Rxvec[iG+i] = VAR_HWBIAS;
        Xini[iG+i]=1;
        init_RxALL(iG+i);
    }
}
/* update ionosphere parameters --------------------------------------------------- */
void gnss_ppp_c::updateion() {
    /* loose constraint for ionosphere */
    for (int i=0; i<satnum; i++) {
        ionfunc->correction(&obsr->data[rovsat[i]],nav,Rblh);
        Xpar[iI+i]=obsr->data[rovsat[i]].dion;
        Rxvec[iI+i]=3.0*(obsr->data[rovsat[i]].ionvar);
        int index_I = usesat[i]-1;
        Xest[iI+i]=iI+index_I;
        Xini[iI+i]=1;
        init_RxALL(iI+index_I);
    }
}
/* update clock parameters -------------------------------------------------------- */
void gnss_ppp_c::updateclk() {
    /* estimate pseudo clock offset for reference system that has no observations in this epoch:
     * assume that GPS is the clock reference system, but GPS is not available and Galileo 
     * is the first system to estimate clock in this epoch: assume tGPS0 and tGAL0 are the
     * clock offset of GPS and Galileo in the last epoch, and ptGPS tGAL1 are the pseudo clock
     * offset of GPS and pre-estimated clock offset of Galileo in this epcoh respectively, then:
     * ptGPS = tGAL1 - tGAL0 + tGPS0 */
    if (clkSys1st!=clkRefSys) {
        sol.back().xclk[clkRefSys] = sol.back().xclk[clkSys1st] - X_ALL[iI+NI+clkSys1st];
    }

    /* update clock for reference system */
    Xpar[ iC + clk_sys[clkRefSys][1] ]= sol.back().xclk[clkRefSys];
    Rxvec[ iC + clk_sys[clkRefSys][1] ] = preClkVar[0];
    Xest[ iC + clk_sys[clkRefSys][1] ]=iI+NI+clkRefSys;
    Xini[ iC + clk_sys[clkRefSys][1] ]=1;
    init_RxALL(iI+NI+clkRefSys);

    /* update clock for other systems */
    for (int i=0; i<NSYS; i++) {
        if ( clk_sys[i][1]<0 || i==clkRefSys ) continue;
        //if ( sol[MAXSOLBUF-2].xclk[i]==0 || ntotal<iniamb || opt->ISBmode==ISBMODE_WHITENOISE ) {
        if ( X_ALL[iI+NI+i]==0 || ntotal<iniamb || opt->ISBmode==ISBMODE_WHITENOISE ) {
            Xpar[ iC + clk_sys[i][1] ] = sol.back().xclk[i] - sol.back().xclk[clkRefSys];
            Rxvec[ iC + clk_sys[i][1] ] = preClkVar[0];
            Xini[ iC + clk_sys[i][1] ]=1;
            init_RxALL(iI+NI+i);
        }
        else {
            if ( clk_sys[i][1]<0 || i==clkRefSys ) continue;
            //Xpar[ iC + clk_sys[i][1] ] = sol[MAXSOLBUF-2].xclk[i] - sol[MAXSOLBUF-2].xclk[clkRefSys];
            //Rxvec[ iC + clk_sys[i][1] ] = sol[MAXSOLBUF-2].vclk[i+i*NSYS] + opt->stdrate[3]*tt;
            Xpar[ iC + clk_sys[i][1] ] = X_ALL[iI+NI+i];
            Rxvec[ iC + clk_sys[i][1] ] = Rx_ALL[ (iI+NI+i)*(N_ALL+1) ] + opt->stdrate[3]*tt;
        }
        Xest[ iC + clk_sys[i][1] ]=iI+NI+i;
    }
}
/* update parameters covariance matrix -------------------------------------------- */
void gnss_ppp_c::updatevar() {
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
void gnss_ppp_c::updatepar(int prePos) {
    iniamb = opt->iniamb;

    /* reset vector and matrix */
    init_arrmat(prePos);

    /* update ssat observation and ambiguity parameters */
    updateamb();

    /* update position/velocity */
    updatexyz(prePos);

    /* initialize geodetic position of rover */
    for (int i=0; i<3; i++) Rxyz[i] = Xpar[i];
    // tidal correction
    if (opt->tidecorr) {
        gtime_c trov = obsr->data[0].time;
        tidefunc.tidecorr(*(trov.timeadd(-obsr->data[0].dtr)->gpst2utc()),0,Rxyz,Rtide);
        tidevar=SQR(0.1*norm(Rtide,3));
        for (int i=0; i<3; i++) Rxyz[i] += Rtide[i];
    }
    // geodetic position of rover,base and centre
    ecef2blh(Rxyz,WGS84,Rblh);

    /* update troposphere parameters */
    if ( opt->tropopt==TROPOPT_EST || opt->tropopt==TROPOPT_ESTG ) updatetro();

    /* update glonass IFB parameters */
    if (NG>0 && (opt->navsys&SYS_GLO)) updateglo();

    /* update ionosphere parameters */
    if (opt->ionoopt==IONOOPT_CONST) updateion();

    /* update clock parameters */
    if (nC>0) {
        preClkVar[0]=VAR_CLK;
        updateclk(); // update clock parameters
    }

    /* update parameters number and Rx */
    numX += nA;
    fix_flag.assign(nA,0);
    Rx.assign(numX*numX,0.0);
    updatevar();
}

/* get code/phase observation ------------------------------------------------- */
double gnss_ppp_c::get_codephase_obs(int isat,int freq,int sys) {
    double Lo = 0.0;
    /* frequency wavelength */
    int fff = freq % numF;
    lam = isat==irfsat[sys]? lam1 : lam2;

    /* stanum sat */
    gnss_ssat_c *sss = &ssat[usesat[isat]-1];

    /* phase observation */
    if (freq<numF) {
        if (ambnum[isat + fff*isat]>=0) {
            /* LC combination */
            if (opt->ionoopt==IONOOPT_IFLC) Lo = sss->lc12[1];
            /* Lfreq */
            else Lo = sss->L[fff].back()*lam[fff];
        }
    }
    /* code observation */
    else {
        /* LC combination */
        if (opt->ionoopt==IONOOPT_IFLC) Lo = sss->pc12[1];
        /* Lfreq */
        else Lo = sss->P[fff].back();
    }

    return Lo;
}
/* satellite-undifferenced geodetic parameters ------------------------------------ */
double gnss_ppp_c::model_distanc(int isat,int sys) {
    /* R:rover B:base */
    double disR;
    Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;

    for (int i=0; i<3; i++) Acoep->at(i) = -rp->sigvec[i];

    /* distance between satellite and rover/base */
    disR = geodist_GNSS(rp->posvel,Rxyz,rp->sigvec);

    satazel(Rblh,rp->sigvec,rp->azel);

    return disR;
}
/* satellite-undifferenced troposphere parameters --------------------------------- */
double gnss_ppp_c::model_troppar(int isat,int sys) {
    Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;

    /* troposphere delay correction and Acoep tro vector */
    for (int i=0; i<NT; i++) trofunc->model_est.Ptro[i] = Xpar[iT+i];
    trofunc->correction(rp,Rblh,0.7);
    for (int i=0; i<NT; i++) Acoep->at(iT+i) = trofunc->model_est.Atro[i];

    return rp->dtro;
}
/* receiver GLONASS IFB parameters ---------------------------------------- */
double gnss_ppp_c::single_gloambd(int isat,int fff,int sys) {
    lam = isat==irfsat[sys]? lam1 : lam2;
    Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;
    /* glo freq-difference parameter and Acoep glo vector */
    Acoep->at(iG+fff) = (CLIGHT/lam[fff]-freq_glos[fff])/1E6;

    return Acoep->at(iG+fff)*Xpar[iG + fff];
}
/* satellite-undifferenced ionosphere parameters ---------------------------------- */
double gnss_ppp_c::model_ionopar(int isat,int sys) {
    double sdion=0.0;
    /* ionosphere correction with certain model */
    if (opt->ionoopt!=IONOOPT_CONST) {
        // rover and isat
        ionfunc->correction(rp,nav,Rblh);
        sdion = rp->dion;
    }

    /* estimate constrained DDion vertical delay of GPS L1 (m) */
    else {
        Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;
        Acoep->at(iI+isat)=rp->ionmap;
        sdion=Acoep->at(iI+isat)*Xpar[iI+isat];
    }

    return sdion;
}
/* satellite-undifferenced hardware delay parameters ------------------------------ */
double gnss_ppp_c::model_hdelays(int isat,int fff,int sys) {
    return 0.0;
}
/* satellite-undifferenced clocks parameters -------------------------------------- */
double gnss_ppp_c::model_clocks(int isat,int sys) {
    Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;

    Acoep->at( iC + clk_sys[clkRefSys][1] )=1.0;
    if (sys!=clkRefSys) {
        Acoep->at( iC + clk_sys[sys][1] )=1.0;
        return Xpar[ iC + clk_sys[clkRefSys][1] ] + Xpar[ iC + clk_sys[sys][1] ];
    }
    else return Xpar[ iC + clk_sys[clkRefSys][1] ];
}
/* satellite-undifferenced ambiguity parameters ----------------------------------- */
double gnss_ppp_c::model_ambtpar(int isat,int fff,int sys) {
    Acoep = isat==irfsat[sys]? Airfsat + sys : &Asatnum;
    double sdamb;
    int namb = ambnum[isat + fff*satnum];
    /* LC ambiguity parameters and Acoep amb vector */
    if (opt->ionoopt==IONOOPT_IFLC) {
        Acoep->at(iA+namb) = 1.0;
        sdamb = Xpar[iA+namb];
    }

    /* freq ambiguity parameters and Acoep amb vector */
    else {
        lam = isat==irfsat[sys]? lam1 : lam2;
        Acoep->at(iA+namb) = lam[fff];
        sdamb = lam[fff]*Xpar[iA+namb];
    }

    return sdamb;
}
/* satellite-undifferenced antenna phase center offset & variation model ---------- */
double gnss_ppp_c::model_satantov(int isat,int fff,int sys) {
    return obsr->data[rovsat[isat]].antv[fff];
}
/* satellite-undifferenced variance ----------------------------------------------- */
double gnss_ppp_c::model_variance(int isat,int freq) {
    int fff = freq % numF;
    int sys = ssat[usesat[isat]-1].sys,prn = obsr->data[rovsat[isat]].prn;
    double fact=opt->pppfactor;
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
    //if (opt->ionoopt==IONOOPT_IFLC) fact*=3;

    /* phase variance */
    if (freq<numF) return SQR(fact) * ( SQR(opt->err[2]) + SQR(opt->err[3]/sin_elR) );
    /* code variance */
    else return SQR(fact) * ( SQR(opt->err[0]) + SQR(opt->err[1]/sin_elR) );
}
/* satellite-undifferenced corrections and variance ------------------------------- */
double gnss_ppp_c::error_correction(int isat,int freq,int sys) {
    int fff=freq%numF;
    /* ionosphere factor */
    ifact2 = (freq<numF?-1.0:1.0)*SQR(lam2[fff]/GNSS_LAMBDA_M[iGPS][1]);

    /* undifferenced geodetic parameters */
    dist=model_distanc(isat,sys);

    /* undifferenced troposphere parameters */
    dtro=model_troppar(isat,sys);

    /* receiver GLONASS IFB parameters */
    gloIFB=0;
    if ( NG>0 && freq<numF&&satsys(usesat[isat],NULL)==SYS_GLO && fff<NFREQGLO ) {
        gloIFB=single_gloambd(isat,fff,sys);
    }

    /* satellite-undifferenced ionosphere parameters */
    dion=model_ionopar(isat,sys);

    /* satellite-undifferenced hardware delay parameters */
    dUPD=model_hdelays(isat,fff,sys);

    /* satellite-undifferenced clocks parameters */
    clock = -CLIGHT*rp->dts[0];
    if (nC>0) clock += model_clocks(isat,sys);

    /* satellite-undifferenced ambiguity parameters */
    ambN=0;
    if (freq<numF) ambN=model_ambtpar(isat,fff,sys);

    /* satellite-undifferenced antenna phase center offset & variation correction */
    /* now it is corrected in observation */
    //santv=model_satantov(isat,fff);

    return dist+dtro+gloIFB+ifact2*dion+dUPD+clock+ambN;
}
/* ihsat-single-differenced parameters and variance and lam1 ---------------------- */
void gnss_ppp_c::irfsat_undifference() {
    for (int sys = 0; sys<NSYS; sys++) {
        Airfsat[sys].assign(numX,0.0);
        //continue if no observation of this system
        if (irfsat[sys]<0) continue;

        /* irfsat gnss_ssat_c and irfsat lamda */
        gnss_ssat_c *shhh = &ssat[rfsat[1][sys] - 1];
        lam1 = nav->lam[rfsat[1][sys] - 1];

        /* initialize rover observation pointer */
        rp = &obsr->data[rovsat[irfsat[sys]]]; 

        /* undifferenced geodetic parameters */
        dist1[sys]=model_distanc(irfsat[sys],sys);

        /* undifferenced troposphere parameters */
        dtro1[sys]=model_troppar(irfsat[sys],sys);

        /* satellite-undifferenced ionosphere parameters */
        dion1[sys]=model_ionopar(irfsat[sys],sys);

        /* satellite-undifferenced clock parameters  */
        clock1[sys] = -CLIGHT*rp->dts[0];
        // receiver clocks have been eliminated

        double irfsat_corr1=dist1[sys]+dtro1[sys]+clock1[sys];

        /* Lobs1 */
        for (int freq=0; freq<numF*2; freq++) {
            double Lo = get_codephase_obs(irfsat[sys],freq,sys);
            if (Lo==0) continue;

            int fff = freq % numF;

            /* receiver GLONASS IFB parameters */
            gloIFB=0;
            if ( NG>0 && freq<numF&&satsys(rfsat[1][sys],NULL)==SYS_GLO && fff<NFREQGLO ) {
                gloIFB=gloIFB1[fff]=single_gloambd(irfsat[sys],fff,sys);
            }

            /* satellite-undifferenced ambiguity parameters */
            ambN=0;
            if (freq<numF) ambN=ambN1[sys][fff]=model_ambtpar(irfsat[sys],fff,sys);

            /* ionosphere factor */
            ifact1 = (freq<numF?-1.0:1.0)*SQR(lam1[fff]/GNSS_LAMBDA_M[iGPS][1]);

            Lobs1[sys][freq]=Lo-irfsat_corr1-ifact1*dion1[sys]-gloIFB-ambN;

            /* irfsat observation variance */
            varRB1[sys][freq] = model_variance(irfsat[sys],freq);
            shhh->vsat[freq%numF] = 1;
        }
    }
}

/* update Rvec vector ------------------------------------------------------------- */
void gnss_ppp_c::update_Rvec(int isat,int freq,int sys) {
    if (opt->pppmode==PPP_SINGLE_DIFF) Rvec.push_back(varRB1[sys][freq]);
    Rvec.push_back(model_variance(isat,freq));
}
/* update Lobs vector and lam2 ---------------------------------------------------- */
int gnss_ppp_c::update_Lobs(int isat,int freq,int sys) {
    double Lo = get_codephase_obs(isat,freq,sys);

    if (Lo==0.0) return 0;
    if (opt->pppmode==PPP_SINGLE_DIFF) {
        if (Lobs1[sys][freq]==0.0) return 0;
        Lobs.push_back(Lo-Lobs1[sys][freq]);
    }
    else Lobs.push_back(Lo);
    return 1;
}
/* update Acoe matrix ------------------------------------------------------------- */
void gnss_ppp_c::update_Acoe(int isat,int freq,int sys) {
    if (opt->pppmode==PPP_SINGLE_DIFF) {
        Ais.assign(numX,0.0);
        /* frequency index */
        int fff = freq % numF;

        /* update Acoe matrix */
        // dynamic pararmeter coefficients
        for (int nx=0; nx<NP; nx++)
            Ais[nx] = Asatnum[nx] - Airfsat[sys][nx];

        // troposphere parameter coefficients
        for (int nx = 0; nx<NT; nx++)
            Ais[iT + nx] = Asatnum[iT + nx] - Airfsat[sys][iT + nx];

        // glo freq-difference parameter coefficients
        if (NG>0 && freq<numF&&sys==1 && fff<NFREQGLO)
            Ais[iG + fff] = Asatnum[iG + fff] - Airfsat[sys][iG + fff];

        // ionosphere parameter coefficients
        if (opt->ionoopt==IONOOPT_CONST) {
            Ais[iI + irfsat[sys]] =  -ifact1*Airfsat[sys][iI + irfsat[sys]];
            Ais[iI + isat]      = ifact2*Asatnum[iI + isat];
        }

        // clock parameter coefficients
        // receiver clocks have been eliminated

        // ambiguity parameter coefficients
        if (freq<numF) {
            int ah = ambnum[ irfsat[sys] + fff*satnum ],as = ambnum[ isat + fff*satnum ];
            if (ah>=0) Ais[iA + ah] = -Airfsat[sys][iA + ah];
            if (as>=0) Ais[iA + as] = Asatnum[iA + as];
        }

        Acoe.insert(Acoe.end(),Ais.begin(),Ais.end());
    }
    else Acoe.insert(Acoe.end(),Asatnum.begin(),Asatnum.end());;
}
/* update Rvar matrix ------------------------------------------------------------- */
void gnss_ppp_c::update_Rvar() {
    /* initialize Rarr */
    Rvar.assign(numL*numL,0.0);

    if (opt->pppmode==PPP_SINGLE_DIFF) {
        for (int sys=0,nt=0; sys<NSYS; sys++) {
            for (int i=0; i<numF*2; i++) {
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
    }
    else {
        for (int i=0; i<numL; i++) Rvar[i+i*numL]=Rvec[i];
    }
}

/* LC ambiguity using wide-lane and N1 -------------------------------------------- */
double gnss_ppp_c::LC_WN(int isat,int sys) {
    return 0.0;
}
/* double-differenced observation equation ---------------------------------------- */
int gnss_ppp_c::ppp_obsequaion() {
    ns = numL = 0;
    for (int i=0; i<NFREQ*2; i++) {
        iobs[i].assign(satnum,-1);
        for (int j=0; j<NSYS; j++) nobs[j][i]=0;
    }

    /* initialize rover and base position */
    if (norm(Rxyz,3) <= 0.0) return 0;

    /* undifferenced observation for reference satellite (irfsat) in SD mode */
    if (opt->pppmode==PPP_SINGLE_DIFF) irfsat_undifference();

    /* loop of different satellites */
    for (int sys=0; sys<NSYS; sys++) {
        //continue if no enough observation of this system
        if ( S_nsat[sys]<MINSYSOBS ) continue; 
        //contiune if no reference satellite in SD mode
        if ( opt->pppmode==PPP_SINGLE_DIFF && irfsat[sys]<0 ) continue;
        else lam1 = nav->lam[rfsat[1][sys] - 1];

        for (int freq=0; freq<numF*2; freq++) {
            for (int isat = 0; isat<satnum; isat++) {

                rp = &obsr->data[rovsat[isat]]; // observation data
                lam2 = nav->lam[usesat[isat]-1]; // lambda

                /* continue if not this system or reference satellite */
                if ( rp->isys!=sys || (opt->pppmode==PPP_SINGLE_DIFF && isat==irfsat[sys]) )
                    continue;

                /* update Lobs vector and lam2 */
                if (!update_Lobs(isat,freq,sys))
                    continue;

                /* initialize Lobs,Rvec,Acoe for new observation */
                Asatnum.assign(numX,0.0);

                /* compute double-difference between rover and isat */
                correction = error_correction(isat,freq,sys);

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
            }
        }
    }

    /* update Rvar using Rvec */
    update_Rvar();

    return 1;
}

/* get fixed solution ------------------------------------------------------------- */
int gnss_ppp_c::get_fixsolution() {
    return 0;
}
/* get fixed wide-lane ambiguity -------------------------------------------------- */
int gnss_ppp_c::fix_wide_narr(vector<double> &DLC,vector<double> &R_DLC) {
    return 0;
}
/* get LC ambiguity (wide-narrow lane) fixed solution ----------------------------- */
int gnss_ppp_c::get_LCfixed() {
    return 0;
}

/* calculate pdop value of each satellite system ---------------------------------- */
void gnss_ppp_c::get_pdop() {
    for (int i=0; i<NSYS+1; i++) Pdop[i]=0.0;
    double pdop_val=0.0;
    if (Acoe.size()==numL*numX) {
        vector<double> Asys[NSYS+1];
        QAA.assign(9,0.0); Asys[NSYS].clear();
        int allobs=0;
        /* pdop for each single system */
        for (int sys=0,nownum=0; sys<NSYS; sys++) {
            int sysobs=0; int codenum=nobs[sys][numF];
            /* observation number */
            for (int fff=0; fff<numF*2; fff++) sysobs+=nobs[sys][fff];
            //update Acoe for all systems
            for (int i=0; i<codenum; i++) for (int j=0; j<3; j++)
                Asys[NSYS].push_back( Acoe[j+(i+nownum)*numX] );
            /* pdop of this system */
            if (codenum>=3) {
                //Asys of Acoe to calculate pdop
                Asys[sys].assign(codenum*3,0.0);
                for (int i=0; i<codenum; i++) for (int j=0; j<3; j++)
                    Asys[sys][j+i*3] = Acoe[j+(i+nownum)*numX];
                //calculate QAA=(Asys'*Asys)^-1
                matmul_vec("TN",3,3,codenum,1.0,Asys[sys],Asys[sys],0,QAA);
                if (matinv(QAA,3)==-1) QAA[0]=QAA[4]=QAA[8]=333.3;
                //calculate pdop 
                double pdop_val=QAA[0]+QAA[4]+QAA[8];
                Pdop[sys] = pdop_val>999.0 ? 999.0 : SQRT(pdop_val);
            }
            else if ( codenum>0 ) Pdop[sys]=999.0;
            allobs+=codenum;
            nownum+=sysobs;
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
void gnss_ppp_c::ave_SNR() {
    int numall=0; SNR[NSYS]=0.0;
    for (int sys=0; sys<NSYS; sys++) {
        SNR[sys]=0.0;
        int codenum=nobs[sys][numF],sysobs=0;
        if (codenum<3) continue;
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
void gnss_ppp_c::ave_multiPath() {
    int numall=0; multiPath[NSYS]=0.0;
    for (int sys=0; sys<NSYS; sys++) {
        multiPath[sys]=0.0;
        int codenum=nobs[sys][numF],sysobs=0;
        if (codenum<3) continue;
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
void gnss_ppp_c::write_state() {
    if (!log_stream.is_open()) return;

    if ( opt->mode>=PMODE_PPP_KINEMA && opt->mode<=PMODE_PPP_FIXED ) {
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
            for (int i=0; i<obsr->n; i++) {
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
            if (opt->pppmode==PPP_SINGLE_DIFF) {
                log_stream << " Reference Sat: ";
                for (int sys = 0; sys<NSYS; sys++)
                    if (irfsat[sys]>=0) log_stream << " " << ssat[rfsat[1][sys] - 1].id;
            }
            /* adjustment time (in the same line of time) */
            /* if (opt->adjustfunc!=ADJUST_HELMERT||opt->adjustfunc==ADJUST_HELMERT&&adjfunc->H_niter)*/
            log_stream << "   Filter Time:" << setw(5) << dtime;
            /* solution number (in the same line of time) */
            log_stream << "   Solution Number:" << setw(7) << nsol;
            /* fixed ratio (in the same line of time) */
            if (opt->modear>=ARMODE_FIXHOLD) {
                if (opt->modear==ARMODE_LCWN) log_stream << "   Fix Rate: " << setw(8) << fixNumRatio;
                else log_stream << "   Fix Var Ratio: " << setw(8) << fixVarRatio;
            }
            log_stream << "\n";

            /* BODY of state ------------------------- */
            /* tide correction */
            if (opt->tidecorr>0) {
                log_stream << "    TideCorrect: ";
                for (int i=0; i<3; i++)
                    log_stream  << setw(9) << Rtide[i];
                log_stream << "\n";
            }
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
            /* receiver clock */
            if (nC>0) {
                log_stream << "    Clock Error: ";
                for (int i=0; i<NSYS; i++) {
                    if (clk_sys[i][1]<0) {
                        log_stream << "   " << SYS_LIST[i] << ":" << setw(13) << 0.0 //new Clock
                        << setw(8) << 0.0; //new Clock STD
                    }
                    else {
                        log_stream << "   " << SYS_LIST[i] << ":" << setw(13) << *(sPar + iC+clk_sys[i][1]) //new Clock
                        << setw(8) << SQRT(*( sRx + (iC+clk_sys[i][1])*(numX+1) )); //new Clock STD
                    }
                }
                log_stream << "   Ref: " << SYS_LIST[clkRefSys];
                log_stream << "\n";
            }
            /* satellites number of different system */
            if (opt->logmsg&GNSSLOG_NSAT) {
                log_stream << "    System nSat: ";
                for (int sys=0; sys<NSYS; sys++) {
                    int nsyso = nobs[sys][numF]>0 ? nobs[sys][numF] : 0;
                    log_stream << "   " << SYS_LIST[sys] << ":" << setw(8) << nsyso;
                }
                log_stream << "   All:" << setw(8) << satnum;
                log_stream << "\n";
            }
            /* PDOP value of different satellite systems */
            if (opt->logmsg&GNSSLOG_PDOP) {
                /* calculate pdop of different satellite system */
                gnss_ppp_c::get_pdop();
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
                rp = &obsr->data[rovsat[i]];
                gnss_ssat_c *sss = &ssat[usesat[i]-1];
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
void gnss_ppp_c::update_ssat() {
    for (int i=0; i<MAXSAT; i++) {
        ssat[i].vs = 0;
    }
    gnss_ssat_c *sss; int namb;
    for (int isat = 0; isat<satnum; isat++) {
        lam = nav->lam[usesat[isat]-1];
        rp = &obsr->data[rovsat[isat]];
        sss = &ssat[usesat[isat]-1];
        sss->vs=1;
        sss->azel[0] = rp->azel[0];  sss->azel[1] = rp->azel[1];
        for (int freq = 0; freq<numF; freq++) {
            sss->snr[freq] = rp->SNR[freq];
            /* ionosphere solution */
            if (nI>0) {
                sss->ion_delay=*(sPar + iI+isat);
                sss->ion_var=*(sRx + iI+isat + (iI+isat)*numX);
            }
            /* ambiguity solution and multipath */
            if ((namb = ambnum[isat+freq*satnum])>=0) {
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
                    sss->amb[freq] = *(sPar+iA+namb);
                    sss->ambvar[freq] = *(sRx + iA+namb + (iA+namb)*numX);
                    /*if (sss->ambvar[freq]<=0) iniamb=iniamb;*/
                    if ((sss->ave_con[freq]<iniamb)) { sss->fix[freq]=1; }
                    else if (stat!=SOLQ_FIX) sss->fix[freq]=2;
                    // update solution indexes
                    sss->solv_con[freq]++;
                    sss->solv_period[freq]+=sss->obstime.back().timediff(sss->solv_time[freq]);
                    sss->solv_time[freq]=sss->obstime.back();
                    if (stat==SOLQ_FIX) {
                        if (nreset[freq]||sss->fix[freq]!=3) {
                            sss->ambvar[freq]=SQR(0.03);
                            for (int i=0; i<numX; i++) //reset covariance matrix
                                *(sRx + i + (iA+namb)*numX) = *(sRx + (iA+namb) + i*numX) = 0.0;
                        }
                        sss->fix[freq]=3;
                    }
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
void gnss_ppp_c::update_sol() {

    soltime[1]=solp->time;
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
    // troposphere parameters
    for (int i=0; i<NT; i++) {
        solp->xtro[i] = *(sPar+iT+i);
        for (int j = 0; j<NT; j++)
            solp->vtro[j + i*NT] = *(sRx + j+iT + (i+iT)*numX);
    }
    // GLO IFB paramters
    for (int i=0; i<NG; i++) {
        solp->xglo[i] = *(sPar+iG+i);
        for (int j = 0; j<NG; j++)
            solp->vglo[j + i*NG] = *(sRx + j+iG + (i+iG)*numX);
    }
    // clock parameters
    for (int i=0; i<NSYS; i++) {
        if (clk_sys[i][1]<0) continue;
        if (i==clkRefSys) solp->xclk[i] = *(sPar+iC+clk_sys[i][1]);
        else solp->xclk[i] = *(sPar+iC+clk_sys[clkRefSys][1]) + *(sPar+iC+clk_sys[i][1]);
        for (int j=0; j<NSYS; j++)
            if (clk_sys[j][1]>=0) solp->vclk[j+i*NSYS] = *(sRx + (clk_sys[i][1]+iC)*(numX+1));
    }
    // ionosphere parameters
    solp->xion.clear(); solp->vion.clear();
    solp->NI=nI;
    for (int i=0; i<nI; i++) {
        solp->xion.push_back(*(sPar+iI+i));
        for (int j = 0; j<nI; j++)
            solp->vion.push_back(*(sRx + j+iI + (i+iI)*numX));
    }
    // ambiguity parameters
    solp->NA = nA;
    solp->xamb.assign(sPar+iA,sPar+iA+nA);
    solp->vamb.assign(nA*nA,0.0);
    for (int i=0; i<nA; i++) for (int j = 0; j<nA; j++)
        solp->vamb[j + i*nA] = *(sRx + j+iA + (i+iA)*numX);

    solp->ns = (unsigned int)satnum;
    if (stat==SOLQ_FIX && opt->modear!=ARMODE_LCWN) solp->ratio=float(sum_var[1]/sum_var[0]);
    solp->stat = stat;

    /* update X_ALL and Rx_ALL */
    // set Rx_ALL from sRx
    for (int i=0; i<numX; i++) {
        X_ALL[ Xest[i] ] = *(sPar+i);
        Rx_ALL[ Xest[i] + Xest[i]*N_ALL ] = *(sRx+i+i*numX); //variance from Rxvec
        for (int j=0; j<i; j++) {
            Rx_ALL[ Xest[j] + Xest[i]*N_ALL ] = Rx_ALL[ Xest[i] + Xest[j]*N_ALL ] = *(sRx+j+i*numX);
        }
    }
}

/* preprocess functions ----------------------------------------------------------- */
/* preprocess gnss observation ---------------------------------------------------- */
int gnss_ppp_c::preprocess_gnss_obs(const double pos[3],const double vel[3]) {
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
int gnss_ppp_c::gnss_obs_equation() {

    /* carrier phase bias correction */
    if (opt->pppopt.find("-DIS_FCB") != string::npos) corr_phase_bias();

    /* tt : time difference between current and previous (s) */
    tt = soltime[1].time!=0 ? obsr->data[0].time.timediff(soltime[1]) : 0.0;

    /* select common satellite and test availability */
    if (selavasat()<MINSAT) {
        msg = "no enough satellites for PPP!";
        return 0;
    }

    /* reset satellites status */
    resetsat();

    /* test */
    if (nsol==iniamb-1)
        opt->niter=opt->niter;

    /* update parameters */
    updatepar(1);

    /* get ppp observation equation */
    ppp_obsequaion();
    Xori=Xpar;

    return 1;
}
/* update gnss status ------------------------------------------------------------- */
int gnss_ppp_c::gnss_update_status() {

    /* point to float solution */
    sPar=Xpar.begin(); sRx=Rx.begin();

    nsol++; ntotal++;

    /* update satellite sate vector (ssat) */
    gnss_ppp_c::update_ssat();
    /* update solution vector (sol) */
    gnss_ppp_c::update_sol();

    stat=SOLQ_FLOAT;
    /* write state information to log_stream */
    gnss_ppp_c::write_state();

    return stat;
}
/* precise point position function ------------------------------------------------ */
int gnss_ppp_c::gnss_pos() {
    warningMsg="";

    /* carrier phase bias correction */
    if (opt->pppopt.find("-DIS_FCB") != string::npos) corr_phase_bias();

    /* if no SPP solution yet, run SPP for rover */
    if (obsr->data[0].time.timediff(sol.back().obsTime)>1E-5) {
        init_sol(1);
        if (single()==SOLQ_NONE) {
            return 0;
        }
    }
    else {
        stat=SOLQ_NONE;
        return 0;
    }

    /* tt : time difference between current and previous (s) */
    tt = soltime[1].time!=0 ? obsr->data[0].time.timediff(soltime[1]) : 0.0;

    /* select common satellite and test availability */
    if (selavasat()<MINSAT) {
        msg = "no enough satellites for PPP!";
        gnss_ppp_c::write_state();
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
        ppp_obsequaion();
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
    if ( stat==SOLQ_FLOAT ) {
        /* uncombined ambiguity */
        if ((opt->modear==ARMODE_FIXHOLD) && get_fixsolution()) stat=SOLQ_FIX;
        /* LC ambiguity */
        else if (nsol>iniamb&&opt->modear==ARMODE_LCWN) {
            stat=SOLQ_FIX;
            fixNumRatio=1.0*fix_num/(1.0*nA);
        }
        else
            stat=SOLQ_FLOAT;
    }

    /* update solution */
    if ( stat==SOLQ_FLOAT || stat==SOLQ_FIX ) {
        /* final solution iterator */
        if (stat==SOLQ_FIX && opt->modear!=ARMODE_LCWN) {
            sPar=fix_Xpar.begin(); sRx=fix_Rx.begin();
        }
        else {
            sPar=Xpar.begin(); sRx=Rx.begin();
        }

        /* update reference satellite to irfsat */
        if (opt->pppmode==PPP_SINGLE_DIFF) {
            for (int i=0; i<NSYS; i++) rfsat[0][i]=rfsat[1][i];
        }

        /* update satellite state */
        gnss_ppp_c::update_ssat();
        /* update solution */
        gnss_ppp_c::update_sol();
    }

    /* write state information to log_stream */
    gnss_ppp_c::write_state();

    return 1;
}
