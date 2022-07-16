#ifndef POSITION_H
#define POSITION_H

#include "GNSS/DataClass/data.h"
#include "GNSS/gnss_pro.h"

/* position types --------------------------------------------------------------------------------- */
/* gnss single point position (SPP) class --------------------------------------------------------- */
class gnss_single_c : public gnss_pro_c {
    /* Constructors */
public:
    gnss_single_c();
    virtual ~gnss_single_c();
    /* Implementation functions */
protected:
    /* correct differenced code bias ---------------------------------------------- */
    void correct_DCB();
    /* update consecutive number for each satellite ------------------------------- */
    void update_sat_nCon();
    /* initialize ephemeris and DCB ----------------------------------------------- */
    int gnss_ini_EPH_DCB();
    /* detect wrong observation with double frequency ----------------------------- */
    void detect_wrong_obs_2freq();
    /* reset satellite status ----------------------------------------------------- */
    void resetsat();
    /* carrier-phase bias (fcb) correction ---------------------------------------- */
    void corr_phase_bias();
    /* select common satellite of rover and base ---------------------------------- */
    void selcomsat();
    /* update observation and ambiguity information for each satellites ----------- */
    void update_obs_amb_ssat();

    /* test SNR mask -------------------------------------------------------------- */
    int test_snr(int rovbas,int freq,double el,double snr,
        const snrmask_t *mask);
    /* pseudorange measurement error variance ------------------------------------- */
    double varerr(const gnss_obsd_c &obs);
    /* verify excluded satellite -------------------------------------------------- */
    int exclude_sat(gnss_obsd_c &obs);
    /* exclude constellations with only a few available satellites ---------------- */
    void exc_few_sat();
    /* exclude observation with large residual ------------------------------------ */
    void exc_largeres(const int last_f);
    /* correct and save pre-etimated clock offset and drift ----------------------- */
    virtual void process_pre_clock();

    /* determine and correct clock jump for every frequency ----------------------- */
    int clock_jump_everyFreq();
    /* psendorange with code bias correction and combination ---------------------- */
    double prange(gnss_obsd_c &obs,const int iter);
    /* compute observation\coefficient\covariance matrix -------------------------- */
    int code_matrix();
    /* compute observation\coefficient\covariance matrix for velocity ----- */
    int doppler_matrix();
    /* compute observation\coefficient\covariance matrix with pre-position -------- */
    int code_matrix_prePos();
    /* compute observation\coefficient\covariance matrix with pre-velocity -------- */
    int doppler_matrix_preVel();
    /* gnss modeled correction (PCV+ClockJump+PhaseWindUp) ------------------------ */
    void model_correction();
    /* update gnss Clock offset --------------------------------------------------- */
    void update_clock_offset();
    /* update gnss modeled correction (Clock+PCV+ClockJump+PhaseWindUP) ----------- */
    void update_model_correction();
    /* check position variance ---------------------------------------------------- */
    int check_pos_var();

    /* calculate pdop value of each satellite system ------------------------------ */
    virtual void get_pdop();
    /* write state information to log_stream -------------------------------------- */
    virtual void write_state();
    /* update satellite sate vector (ssat) ---------------------------------------- */
    virtual void update_ssat();
    /* update solution vector (sol) ----------------------------------------------- */
    virtual void update_sol();

    /* growth variance of static rover position ----------------------------------- */
    double growth_rate_static();

    /* initialize covariance matrix of amb parameters ----------------------------- */
    void init_RxAMB(unsigned int sat,int freq);
    /* initialize covariance matrix of all parameters ----------------------------- */
    void init_RxALL(int ipar);

    /* estimate receiver position ------------------------------------------------- */
    int singlepos();
    /* single point position ------------------------------------------------------ */
    int single();
public:
    /* initialize solution vectors ------------------------------------------------ */
    void init_sol(int rovbas);
    /* preprocess gnss observation ------------------------------------------------ */
    virtual int preprocess_gnss_obs(const double pos[3], const double vel[3]);
    /* get gnss observation equation (SPP) ---------------------------------------- */
    virtual int gnss_obs_equation();
    /* update gnss status (SPP) --------------------------------------------------- */
    virtual int gnss_update_status();
    /* gnss velocity function ----------------------------------------------------- */
    int gnss_vel();
    /* base station single-position ----------------------------------------------- */
    int basepos();
    /* virtual gnss position function --------------------------------------------- */
    virtual int gnss_pos();

/* Components */
public:
    /* status */
    int niter;                          /* number of iteration */
    double res_ave;                     /* average residual */
    int dtime;                          /* cost time (ms) */
    gtime_c lastTime;                   /* last observation time */
    int lsqflag;                        /* least-square adjustment flag */
    double posVar;                      /*  */

    /* observation parameters */
    gnss_obsd_c *rp,*bp;                /* obs. pointor of rover and base */
    double *lam,*lam1,*lam2;            /* pointor to lamda vector in nav */

    /* process parameters */
    double Lnew;                        /* new Lobs element */

    /* used common satellite (for ppp and rtk process, set in selcomsat_phase()) */
    vector<unsigned int> usesat;        /* common satellites vector (rover and base) */
    vector<int> rovsat,bassat;          /* common satellites index in rover and base obs */
    int satnum;                         /* common satellites number (rover and base) */
    /* observation parameters (for ppp and rtk process, set in selcomsat_phase()) */
    int L_nsat[NSYS][NFREQ];            /* number of satellites that have available phase obs.
                                         * of NFREQ freq. of NSYS system  */
    int nobs[NSYS][NFREQ*2];            /* number of used NSYS system code/phase observation */
    vector<int> iobs[NFREQ*2];          /* observation index of used satellites */

    /* antenna correction */
    double dsant,drant;                 /* satellites and receiver PCO+PCV */

    /* clock jump detection */
    double CJ_SF;                       /* single clock jump float value */
    double CJ_float[NFREQ];             /* clock jump detected float value in current epoch (ms) */
    int CJ_nsat[NFREQ];                 /* number of satellites detected clock jump in current epoch */
    int CJ_ndetect[NFREQ];              /* number of satellites for clock jump detection in current epoch */
    double CJ_corr[NFREQ];              /* clock jump correction value in current epoch (m) */
    double CJ_sum[NFREQ];               /* cumulative clock jump values since the last reset (m) */
    int CJ_con[NFREQ];                  /* cumulative number of clock jump appearance */
    vector<double> CJ_n100[NFREQ];      /* last 100 clock jumps (ms) */
    double preClkOff[NSYS],preClkDrf;   /* pre-esitmated clock offset and drift */
    double preClkVar[2];                /* pre-estimated clock offset and drift variance */
    int clkRefSys;                      /* clock reference system */
    int clkSys1st;                      /* first system to estimate clock */

    double Pdop[NSYS+1];                /* positioning dilution of precision */
    double SNR[NSYS+1];                 /* average SNR of NSYS GNSS systems */
    double multiPath[NSYS+1];           /* average multipath of NSYS GNSS systems */
    int resflag;                        /* exclude observation with high residual flag */
    int ratflag;                        /* exclude observation with large residual ratio */
    int staflag;                        /* station flag (0:rover, 1:base) */
    int ncode[2][NSYS];                 /* code observation number of each systems [0:rov,1:base] */
    int *ncd;                           /* pointer of ncode[2] */

    /* parameters for estimation of velocity using Doppler */
    int ndoppler[2][NSYS];              /* Doppler observation number of each systems [0:rov,1:base] */
    int *ndp;                           /* pointer of ndoppler[2] */

    /* ambiguity parameters (some initialized in doublediff()) */
    int iniamb;                         /* set initial ambiguity if lock count < iniar */
    int nreset[NFREQ+1];                /* number of reseted ambiguity of each frequency (L1,2,3 and LC) */
    int reset_flag[NFREQ+1];            /* reset flag of each frequency ambiguity (L1,2,3 and LC) */
    int reset_x;                        /* reset flag of all ambiguities */
    int n_elAmb;                        /* number of high elevation ambiguity which will be fixed */
    vector<int> fixIndex;               /* index of high elevation ambiguity which will be fixed */
    vector<int> ambnum;                 /* amb number of one satellite frequency */
    vector<double> ambElevation[2];     /* satellite elevation for fixing ambiguity (0:rover,1:base) */
    vector<double> fix_amb,sum_var;     /* fix_amb: fixed ambiguity vector */
                                        /* sum_var: sum of squared residuals of fix_amb */

    /* fixed solution */
    int fix_num;                        /* number of fixed ambiguity */
    double fixNumRatio,fixVarRatio;     /* fixed number/variance ratio */
    vector<int> fix_flag;               /* fix flag of each DD ambiguity */
    vector<double> fix_Xpar,fix_Rx;     /* fixed parameters and covariance */

    /* messages */
    string warningMsg;                  /* warning messages */
};

/* precise point position class ------------------------------------------------------------------- */
class gnss_ppp_c : public gnss_single_c {
/* Constructors */
public:
    gnss_ppp_c();
    virtual ~gnss_ppp_c();
/* Implementation functions */
protected:
    /* correct and save pre-etimated clock offset and drift ----------------------- */
    virtual void process_pre_clock();
    /* test satellite system ------------------------------------------------------ */
    int test_system(int Sysyem,int System_Num);
    /* select available satellites of rover --------------------------------------- */
    int selavasat();

    /* initialize reference satellites parameters --------------------------------- */
    void init_refsat();
    /* initialize vector and matrix according common satellite ---------------------*/
    void init_arrmat(int prePos);

    /* update parameters functions ------------------------------------------------ */
    /* update ambiguity parameters ------------------------------------------------ */
    void updateamb();
    /* update geodetic parameters ------------------------------------------------- */
    void updatexyz(int prePos);
    /* update troposphere parameters ---------------------------------------------- */
    void updatetro();
    /* update receiver GLONASS IFB parameters ------------------------------------- */
    void updateglo();
    /* update ionosphere parameters ----------------------------------------------- */
    void updateion();
    /* update clock parameters ---------------------------------------------------- */
    void updateclk();
    /* update parameters covariance matrix ---------------------------------------- */
    void updatevar();
    /* update parameters from previous time to current time ----------------------- */
    void updatepar(int prePos);

    /* get code/phase observation ------------------------------------------------- */
    double get_codephase_obs(int isat,int freq,int sys);
    /* undifferenced geodetic parameters ------------------------------------------ */
    double model_distanc(int isat,int sys);
    /* undifferenced troposphere parameters --------------------------------------- */
    double model_troppar(int isat,int sys);
    /* receiver GLONASS IFB parameters ------------------------------------------- */
    double single_gloambd(int isat,int fff,int sys);
    /* satellite-undifferenced ionosphere parameters ------------------------------ */
    double model_ionopar(int isat,int sys);
    /* satellite-undifferenced hardware delay parameters -------------------------- */
    double model_hdelays(int isat,int fff,int sys);
    /* satellite-undifferenced clocks parameters ---------------------------------- */
    double model_clocks(int isat,int sys);
    /* satellite-undifferenced ambiguity parameters ------------------------------- */
    double model_ambtpar(int isat,int fff,int sys);
    /* satellite-undifferenced antenna phase center offset & variation model ------ */
    double model_satantov(int isat,int fff,int sys);
    /* satellite-undifferenced variance ------------------------------------------- */
    double model_variance(int isat,int freq);
    /* satellite-undifferenced corrections and variance --------------------------- */
    double error_correction(int isat,int freq,int sys);
    /* ihsat-single-differenced parameters and variance and lam1 ------------------ */
    void irfsat_undifference();

    /* update Rvec vector --------------------------------------------------------- */
    void update_Rvec(int isat,int freq,int sys);
    /* update Lobs vector and lam2 ------------------------------------------------ */
    int update_Lobs(int isat,int freq,int sys);
    /* update Acoe matrix --------------------------------------------------------- */
    void update_Acoe(int isat,int freq,int sys);
    /* update Rvar matrix --------------------------------------------------------- */
    void update_Rvar();

    /* LC ambiguity using wide-lane and N1 ---------------------------------------- */
    double LC_WN(int isat,int sys);
    /* double-differenced observation equation ------------------------------------ */
    int ppp_obsequaion();

    /* get fixed solution --------------------------------------------------------- */
    int get_fixsolution();
    /* get fixed wide-lane ambiguity ---------------------------------------------- */
    int fix_wide_narr(vector<double> &DLC,vector<double> &R_DLC);
    /* get LC ambiguity (wide-narrow lane) fixed solution ------------------------- */
    int get_LCfixed();

    /* calculate pdop value of each satellite system ------------------------------ */
    virtual void get_pdop();
    /* average S/N of different systems ------------------------------------------- */
    void ave_SNR();
    /* average multipath of different systems ------------------------------------- */
    void ave_multiPath();
    /* write state information to log_stream -------------------------------------- */
    virtual void write_state();
    /* update satellite sate vector (ssat) ---------------------------------------- */
    virtual void update_ssat();
    /* update solution vector (sol) ----------------------------------------------- */
    virtual void update_sol();

    /* preprocess functions ------------------------------------------------------- */

public:
    /* preprocess gnss observation ------------------------------------------------ */
    virtual int preprocess_gnss_obs(const double pos[3], const double vel[3]);
    /* get gnss observation equation (PPP) ---------------------------------------- */
    virtual int gnss_obs_equation();
    /* update gnss status (PPP) --------------------------------------------------- */
    virtual int gnss_update_status();
    /* precise point position function -------------------------------------------- */
    int gnss_pos();

/* Components */
protected:
    /* station geographic parameters (ini. in  updatepar()) */
    double Rxyz[3],Rblh[3];             /* rover position with tidal correction (ecef and geodetic) */
    double Rtide[3],tidevar;            /* tide correction for rover and its variance */


    /* GNSS system information */
    int S_nsat[NSYS];                   /* number of available satellites of each GNSS system */
    int clk_sys[NSYS][2];               /* clock index for GNSS systems (0:last,1:now) */

    /* centre satellite parameters (irfsat) */
    int irfsat[NSYS];                   /* index of reference satellite of NSYS system */
    unsigned int rfsat[2][NSYS];        /* centre satellite of NSYS system (0:last,1:now) */

    /* PPP process parameters */
    double Lobs1[NSYS][NFREQ*2],dist1[NSYS],dtro1[NSYS],gloIFB1[2],dion1[NSYS],
           clock1[NSYS],ambN1[NSYS][NFREQ];
                                        /* chsat undifferenced correction for reference satellite
                                         * ionRB1 represents L1 ion delay */
    double varRB1[NSYS][NFREQ*2];       /* chsat undifferenced variance for reference satellite */

    /* other satellite parameters (isat) */
    double dist,dtro,gloIFB,dion,dUPD,clock,ambN;
    double correction;                  /* double-differenced correction between rover/base and irfsat/isat */

    /* parameters number (ini. in init_arrmat()) */
    int iT,iG,iI,iU,iC,iA;              /* start index of kinds of parameters */
    int nI,nU,nC,nA;                    /* ionosphere, ambiguity parameters (GLO, and all amb) */

    /* ionosphere parameters  */
    double ifact1,ifact2;               /* frequency factor of ionosphere delay */

    /* Rx and Acoe buff */
    vector<double> Rxvec;               /* parameters variance vector */
    vector<double> *Acoep;              /* Acoe pointer to Airfsat[i] or Asatnum */
    vector<double> Airfsat[NSYS];       /* Acoe vector of rfsat single-difference */
    vector<double> Asatnum,Ais;         /* Acoe vector of this satellite */
};

/* relative position class ------------------------------------------------------------------- */
class gnss_relative_c : public gnss_single_c {
/* Constructors */
public:
    gnss_relative_c();
    virtual ~gnss_relative_c();
/* Implementation functions */
protected:
    /* correct and save pre-etimated clock offset and drift ----------------------- */
    virtual void process_pre_clock();
    /* test satellite system ------------------------------------------------------ */
    int test_system(int Sysyem,int System_Num);
    /* select common satellites of rover and base considering phase observation --- */
    int selcomsat_phase();

    /* initialize reference satellites parameters --------------------------------- */
    void init_refsat();
    /* initialize vector and matrix according common satellite -------------------- */
    void init_arrmat(int prePos);

    /* update parameters functions ------------------------------------------------ */
    /* update ambiguity parameters ------------------------------------------------ */
    void updateamb();
    /* update dynamic parameters -------------------------------------------------- */
    void updatexyz(int prePos);
    /* update troposphere parameters ---------------------------------------------- */
    void updatetro();
    /* update receivers GLONASS IFB parameters ------------------------------------ */
    void updateglo();
    /* update ionosphere parameters ----------------------------------------------- */
    void updateion();
    /* update parameters covariance matrix ---------------------------------------- */
    void updatevar();
    /* update parameters from previous time to current time ----------------------- */
    void updatepar(int prePos);

    /* time difference between rover and base ------------------------------------- */
    int timediff_rb();
    /* satellite-single-differenced dynamic parameters ---------------------------- */
    double single_distanc(int isat,int sys);
    /* satellite-single-differenced troposphere parameters ------------------------ */
    double single_troppar(int isat,int sys);
    /* satellite-single-differenced receivers GLONASS IFB parameters -------------- */
    double single_gloambd(int isat,int fff,int sys);
    /* satellite-single-differenced ionosphere parameters ------------------------- */
    double single_ionopar(int isat,int sys);
    /* satellite-single-differenced ambiguity parameters -------------------------- */
    double single_ambtpar(int isat,int fff,int sys);
    /* satellite-single-differenced variance -------------------------------------- */
    double single_variance(int isat,int freq);
    /* ihsat-single-differenced parameters and variance and lam1 ------------------ */
    void irfsat_single();

    /* update Rvec vector --------------------------------------------------------- */
    void update_Rvec(int isat,int freq,int sys);
    /* update Lobs vector and lam2 ------------------------------------------------ */
    int update_Lobs(int isat,int freq,int sys);
    /* update Acoe matrix --------------------------------------------------------- */
    void update_Acoe(int isat,int freq,int sys);
    /* update Rvar matrix --------------------------------------------------------- */
    void update_Rvar();

    /* double-differenced LC ambiguity using wide-lane and N1 --------------------- */
    double LC_WN(int isat,int sys);
    /* baseline-constraint equation for moving-baseling --------------------------- */
    void base_line();
    /* double-differenced value --------------------------------------------------- */
    double double_value(int isat,int freq,int sys);
    /* double-differenced observation equation ------------------------------------ */
    int double_diff();

    /* single to double-difference tansformation matrix (Dx_coe) ------------------ */
    int single2doul(vector<double> &Dx_coe);
    /* double-differenced ambiguity to single ambiguity --------------------------- */
    int ddamb2single(vector<double> &Dx_coe,vector<double> &Damb_Fix,
        vector<double> &R_Damb);
    /* get fixed solution --------------------------------------------------------- */
    int get_fixsolution();
    /* single to double-difference tansformation matrix (considering elevation) --- */
    int single2doul_el(vector<double> &Dx_coe);
    /* double-differenced ambiguity to single ambiguity (considering elevation) --- */
    int ddamb2single_el(vector<double> &Dx_coe,vector<double> &Damb_Fix,
        vector<double> &R_Damb,vector<double> &Rxa,vector<double> &dxa);
    /* get fixed solution (considering elevation) --------------------------------- */
    int get_fixsolution_el();
    /* get fixed wide-lane ambiguity ---------------------------------------------- */
    int fix_wide_narr(vector<double> &DLC,vector<double> &R_DLC);
    /* get LC ambiguity (wide-narrow lane) fixed solution ------------------------- */
    int get_LCfixed();

    /* update baseline status ----------------------------------------------------- */
    void update_baseline();
    /* calculate pdop value of each satellite system ------------------------------ */
    virtual void get_pdop();
    /* average S/N of different systems ------------------------------------------- */
    void ave_SNR();
    /* average multipath of different systems ------------------------------------- */
    void ave_multiPath();
    /* write state information to log_stream -------------------------------------- */
    virtual void write_state();
    /* update satellite sate vector (ssat) ---------------------------------------- */
    virtual void update_ssat();
    /* update solution vector (sol) ----------------------------------------------- */
    virtual void update_sol();
public:
    /* preprocess gnss observation ------------------------------------------------ */
    virtual int preprocess_gnss_obs(const double pos[3], const double vel[3]);
    /* get gnss observation equation (relative) ----------------------------------- */
    virtual int gnss_obs_equation();
    /* update gnss status (relative) ---------------------------------------------- */
    virtual int gnss_update_status();
    /* relative position function ------------------------------------------------- */
    int gnss_pos();
/* Components */
protected:
    /* station geographic parameters (ini. in  updatepar()) */
    double odt;                         /* time difference between obs. of rover and base (s) */
    double Rxyz[3],Bxyz[3];             /* rover and base position with tidal correction (ecef) */
    double Rblh[3],Bblh[3];             /* rover and base position (geodetic) */
    double Rtide[3],Btide[3];           /* tide correction for rover and base */
    double baseline;                    /* base-line length */
    double Rbl;                         /* variance of baseline constraint */

    /* reference satellite parameters (irfsat) */
    unsigned int rfsat[2][NSYS];        /* centre satellite of NSYS system (0:last,1:now) */
    int irfsat[NSYS];                   /* index of reference satellite of NSYS system */
    double disRB1[NSYS],troRB1[NSYS],gloRB1[2],ionRB1[NSYS],ambRB1[NSYS][NFREQ];
                                        /* chsat single-differenced correction for reference satellite
                                         * ionRB1 represents L1 ion delay */
    double varRB1[NSYS][NFREQ*2];       /* chsat single-differenced variance for reference satellite */

    /* other satellite parameters (isat) */
    double disRB2,troRB2,gloRB2,ionRB2,ambRB2;
                                        /* isat single-differenced correction
                                         * ionRB2 represents L1 ion delay */
    double correction;                  /* double-differenced correction between rover/base and irfsat/isat */

    /* parameters (ini. in init_arrmat()) */
    int iT,iG,iI,iA;                    /* start index of kinds of parameters */
    int nI,nA;                          /* ionosphere, ambiguity parameters (GLO, and all amb) */
    int nTrp;                           /* tropospheric parameter number for one station */

    /* Rx and Acoe buff */
    vector<double> Rxvec;               /* parameters variance vector */
    vector<double> *Acoep;              /* Acoe pointer to Airfsat[i] or Asatnum */
    vector<double> Airfsat[NSYS];       /* Acoe vector of rfsat single-difference */
    vector<double> Asatnum,Ais;         /* Acoe vector of atnum single-difference and
                                        * Acoe vector of double-difference */

    /* ionosphere parameters  */
    double ifact1,ifact2;               /* frequency factor of ionosphere delay */

/* fixed solution */
    int n_Damb;                         /* number of all DD ambiguity */
    int En_Damb;                        /* number of all ambiguity (without GLO) */
};
#endif
