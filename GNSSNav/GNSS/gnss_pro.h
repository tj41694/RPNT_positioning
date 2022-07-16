#ifndef GNSS_PRO_H
#define GNSS_PRO_H

#include "Decode/decode.h"
#include "ConfigFile/config.h"
#include "RtkStream/stream.h"
#include "InputFile/inputfile.h"
#include "GNSS/DataClass/data.h"

#include "Adjustment/adjustment.h"
#include "GNSS/AmbModel/ambiguity.h"
#include "GNSS/EphModel/satellite.h"
#include "GNSS/AntModel/antenna.h"
#include "GNSS/TidModel/tide.h"
#include "GNSS/IonModel/ionosphere.h"
#include "GNSS/TroModel/troposphere.h"
#include "GNSS/ParModel/parameter.h"

/* GNSS process class ----------------------------------------------------------------------------- */
class gnss_pro_c {
    /* Constructor */
public:
    gnss_pro_c();
    virtual ~gnss_pro_c();
    /* Implementation functions */
protected:
public:
    /* initialize GNSS process ---------------------------------------------------- */
    void gnss_pro_init();
    /* get receiver position and velocity from external --------------------------- */
    void gnss_get_posvel(const double pos[3],const double vel[3],
                         const double vpos[9],const double vvel[9]);
    /* preprocess gnss observation ------------------------------------------------ */
    virtual int preprocess_gnss_obs(const double pos[3],const double vel[3]);
    /* get gnss observation equation (virtual) ------------------------------------ */
    virtual int gnss_obs_equation();
    /* update gnss status (virtual) ----------------------------------------------- */
    virtual int gnss_update_status();
    /* gnss position function (virtual) ------------------------------------------- */
    virtual int basepos();
    /* gnss velocity function ----------------------------------------------------- */
    virtual int gnss_vel();
    /* gnss position function (virtual) ------------------------------------------- */
    virtual int gnss_pos();
    /* Components */
public:
    /* gnss option */
    gnss_prcopt_c *opt;                 /* GNSS processing options */

    /* processing status */
    gtime_c soltime[3];                 /* last solution time (0:spp, 1:ppp, 2:rtk) */
    vector<gnss_sol_c>  sol;            /* RTK solution */
    vector<gnss_sol_c> b_sol;           /* RTK base solution */
    double rb[6];                       /* base position/velocity (ecef) (m|m/s) */
    double tt;                          /* time difference between current and previous (s) */
    int nfix;                           /* number of continuous fixes of ambiguity */
    vector<gnss_ssat_c> ssat;           /* satellite status */
    int neb;                            /* bytes in error message buffer */
    string errbuf;                      /* error message buffer */
    string msg;                         /* error message */
    int usePhase;                       /* flag of using phase observation (0:no, 1:yes) */
    int nsol,ntotal;                    /* number of continue and total solution */
    int stat,vstat;                     /* positioning and velocity status */
    int ns;                             /* number of valid satellite */
    int numsys;                         /* number of GNSS systems */

    /* position classes */
    gnss_obs_c *obsr,*obsb;             /* observation data (rover, base) */
    gnss_nav_c *nav;                    /* navigation data */
    gnss_sol_c *solp;                   /* gnss_sol_c pointer */
    gnss_obs_c *obsp;                   /* gnss_obs_c pointer */

    /* parameters */
    int NF,numF;                        /* number of positioning frequency and frequency in filter */
    int NP,NV,NT,NG,NC,NU,NI,NA;        /* number of each kind of parameters */
    int NX;                             /* number of all parameters except ambiguity */
    int N_ALL;                          /* number of all parameters with all satellite ambiguity */
    int numX;                           /* current parameters number */
    vector<double> X_ALL;               /* solution of all parameters */
    vector<double> Rx_ALL;              /* covariance of all parameters */
    vector<int> Xest,Xini;              /* estimated index and initialized flag of parameters  */
    int numL;                           /* number of used observation (size of Lobs) */
    int vnumL;                          /* number of doppler observations (size of vLobs) */

    vector<double> Lobs,Acoe,Rvar,Rvec,Xpar,Xori,Rx,QAA;
                                        /* Lobs: observation vector */
                                        /* Acoe: coefficients matrix */
                                        /* Rvar: obs covariance matrix */
                                        /* Rvec: obs variance vector */
                                        /* Xpar: parameters vector */
                                        /* Xori: original parameters vector */
                                        /* Rx  : parameters covariance matrix */
                                        /* QAA : A*A' */

    vector<double> vLobs,vAcoe,vRvar,vRvec,vXpar,vRx;
                                        /* vLobs: observation vector (velocity) */
                                        /* vAcoe: coefficients matrix (velocity) */
                                        /* vRvar: obs covariance matrix (velocity) */
                                        /* vRvec: obs variance vector (velocity) */
                                        /* vXpar: parameters vector (velocity) */
                                        /* vRx  : parameters covariance matrix (velocity) */

    int CJ_nfreq;                       /* number of frequency to detect clock jump */

    /* final solution iterator */
    vector<double>::iterator sPar,sRx;

    double ext_posvel[6];               /* external receiver position and velocity  */
    double ext_vpos[9],ext_vvel[9];     /* external receiver position and velocity covariance  */
    int varFlag;                        /* external receiver position and velocity covariance flag */
    int preClkFlag;                     /* flag to pre-estimate clock offset and drift */
    double att_n_b[3],att_err[3];       /* attitude and error for moving-baseline (deg) */

    /* function classes */
    gnss_parameter_c parafunc;          /* pararmeter number function */
    gnss_satant_c satantfunc;           /* satellite antenna functions (for satfuc) */
    gnss_recant_c recantfunc;           /* receiver antenna functions */
    gnss_satellite_c *satfunc;          /* satellite function class point */
    gnss_tidecorr_c tidefunc;           /* tidal displacement correction functions */
    gnss_ioncorr_c *sppionf;            /* SPP ionosphere delay functions (no estimate) */
    gnss_trocorr_c *spptrof;            /* SPP troposphere delay functions (no estimate) */
    gnss_ioncorr_c *ionfunc;            /* ionosphere delay functions */
    gnss_trocorr_c *trofunc;            /* troposphere delay functions */
    adjfunc_c *adjfunc;                 /* adjustment functions */
    gnss_lambda_c lambda;               /* ambiguity fix function */
    kalmanfilter_c ambFilter;           /* filter for fix and hold ambiguity */

    /* state file */
    string log_path;                    /* log output path */
    fstream log_stream;                 /* log output stream */
};

/* GNSS solution class ---------------------------------------------------------------------------- */
class gnss_sol_c {
/* Constructor */
public:
    gnss_sol_c();
    gnss_sol_c(const gnss_pro_c *gnss_pro);
    ~gnss_sol_c();
/* Implementation functions */
protected:
    /* dynamic covariance to ecef covariance -------------------------------------- */
    void dyc2ecef(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro);
    /* ecef solution -------------------------------------------------------------- */
    void ecef(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro);
    /* ecef position to BLH ------------------------------------------------------- */
    void blh(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro);
    /* ecef position to ENU ------------------------------------------------------- */
    void enu(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro);
    /* ecef position to EMEA ------------------------------------------------------ */
    void nmea(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro);
public:
    /* solution position ---------------------------------------------------------- */
    string forposvel(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro);
    /* solution time -------------------------------------------------------------- */
    string fortime(const gnss_solopt_c *opt);
/* Components */
public:
    gtime_c time;                       /* time (GPST) */
    gtime_c obsTime;                    /* observation time (GPST) */
    unsigned int type;                  /* type (0:xyz-ecef,1:enu-baseline) */
    unsigned int stat;                  /* solution status (SOLQ_???) */
    unsigned int ns;                    /* number of valid satellites */
    float age;                          /* age of differential (s) */
    float ratio;                        /* AR ratio factor for valiation */
    float thres;                        /* AR ratio threshold for valiation */
    int nsol;                           /* solution number */

    /* parameters (only for position mode higher than SPP) */
    int NF,NL;                          /* number of used frequency and observation */
    int NP,NV,NI,NT,NG,NC,NA;           /* number of each kind of parameters */
    vector<double>                      /* parameter vectors */
        xpos,xvel,xtro,xglo,xion,xclk,xamb;
                                        /* dynamic parameters */
                                        /* ionosphere parameters */
                                        /* troposphere parameters */
                                        /* GLO receiver differenced IFB rate (only for relative) */
                                        /* receiver clock parameters (only for ppp) */
                                        /* ambiguity parameters */
    vector<double>
        vpos,vvel,vtro,vglo,vion,vclk,vamb;
                                        /* covariance of parameters */



    /* formated solution (data and string) */
    gtime_c soltime;                    /* solution time */
    string  strtime;

    double refpos[3];                   /* reference position for enu */
    double posvel[6];                   /* position and velocity (m|m/s or deg|deg/s) */
    double posvar[9];                   /* position covariance (m^2 or deg^2) */
    double velvar[9];                   /* velocity covariance (m^2/s^2) */
    double ecefPvar[9],ecefVvar[9];     /* position & velocity covariance in ECEF */
    double lat[3];                      /* latitude (ddd mm ss) */
    double lon[3];                      /* longitude (ddd mm ss) */
    string strpv;
};

/* satellite status class ------------------------------------------------------------------------- */
class gnss_ssat_c {
/* Constructor */
public:
    gnss_ssat_c();
    ~gnss_ssat_c();
    /* Implemetaion functions */
protected:
    /* set coefficients matrix for estimate of gf-slip ---------------------------- */
    void coe_gfslip(const int flag,vector<double> &An,vector<double> &Rn);
    /* detect cycle slip functions ------------------------------------------------ */
    /* detect by LLI -------------------------------------------------------------- */
    void detslip_LLI(const gnss_obsd_c *rov,const gnss_obsd_c *bas);
    /* detect by Code-Phase value ------------------------------------------------- */
    void detslip_PC(const gnss_prcopt_c *opt);
    /* detect by geometry-free combination ---------------------------------------- */
    void detslip_gf(const gnss_prcopt_c *opt);
    /* detect by Melbourne-Wubbena ------------------------------------------------ */
    void detslip_MW(const gnss_prcopt_c *opt);
    /* calculate new ambiguity and initialize statistc data ----------------------- */
    void new_amb_statistic();
public:
    /* initialize vectors with order of polynomial fitting ------------------------ */
    void init_vector(const gnss_prcopt_c *opt);
    /* update carrier wave lengths ------------------------------------------------ */
    void update_lambda(const gnss_nav_c &nav);
    /* update observations for current epoch -------------------------------------- */
    void input_new_obs(const gnss_obsd_c *rov,const gnss_obsd_c *bas,
        const double CJ[NFREQ]);
    /* update process parameter for current epoch --------------------------------- */
    void update_pro_par(const double CJ[NFREQ]);
    /* reset flag according to unsolved time interval ----------------------------- */
    void test_reset(const gnss_prcopt_c *opt);
    /* reset LC ambiguity --------------------------------------------------------- */
    void test_reset_LC();
    /* reset_ambiguity ------------------------------------------------------------ */
    void reset_amb(const int freq);
    /* detect cycle slip ---------------------------------------------------------- */
    void detect_slip(const gnss_obsd_c *rov,const gnss_obsd_c *bas,
        const gnss_prcopt_c *opt);
    /* update ambiguity statistic ------------------------------------------------- */
    void update_amb_statistic(const gnss_prcopt_c *opt);
    /* correct antenna and phase windup for observations and ambiguities ---------- */
    void correct_obs_amb(const int freq,const gnss_obsd_c *rov,const gnss_obsd_c *bas);
    /* update ambiguity parameters to be solved ----------------------------------- */
    int update_solved_amb(const int freq,const int nsol);

/* Components */
public:
    string errorMsg;
    /* state parameters */
    unsigned int sys;                   /* navigation system */
    unsigned int sat;                   /* satellite number */
    string id;                          /* satellite id */
    unsigned int vs;                    /* valid satellite flag single */
    double azel[2];                     /* azimuth/elevation angles {az,el} (rad) */
    double resP[NFREQ];                 /* residuals of pseudorange (m) */
    double resL[NFREQ];                 /* residuals of carrier-phase (m) */
    unsigned int vsat[NFREQ];           /* valid satellite flag */
    double snr[NFREQ];                  /* signal strength (dBHz) */
    int nCon,minCon;                    /* number and minimum of consecutive osbservation */

    /* observation/solution parameters vector (size = gnss_prcopt_c->order) */
    int             nfreq;              /* number of frequencies */
    vector<gtime_c> obstime;            /* observtaion time */
    vector<double>  L[NFREQ];           /* phase observation (cycle) */
    vector<double>  P[NFREQ];           /* code observation (m) */
    vector<double>  D[NFREQ];           /* doppler frequency */
    vector<double>  Lm[NFREQ];          /* phase observation (m) */
    double gf12[2],gf15[2];             /* geometry-free phase L1-L2/L1-L5 (m) (0:last,1:new) */
    double lc12[2],pc12[2];             /* ionosphere-free phase/code L1-L2 (m) (0:last,1:new) */
    double mw12[2],mw15[2];             /* Melbourne-Wubbena L1-L2/L1-L5 (n) (0:average,1:new) */

    /* ionosphere parameters */
    double ion_delay;                   /* vertical ionosphere delay of GPS L1 (m) */
    double ion_var;                     /* variance of ion_var */
    int ion_flag;                       /* ionosphere delay correction flag (0:no,1:yes) */
    vector<double> d_ion;               /* ionosphere delay of each frequency of L1 */

    /* ambiguity parameters */
    double lambda[NFREQ];               /* carrier wave lengths (m) */
    vector<unsigned int> slip[NFREQ];   /* cycle-slip flag vector
                                         * freq(f+1) slip flag
                                         * (bit7-6: rov LLI,   bit5-4: bas LLI,
                                         * bit2   : repaired
                                         * bit1   : half slip, bit0  : slip) */
    int           slip_f[NFREQ];        /* slip detected flag of different strategy
                                         * LLI, poly, gf, MW */
    double        phw[2];               /* phase windup for rover and base (cycle) */
    int           amb_obsvb[NFREQ+1][2];/* ambiguity observability, 0:last,1:new (L1,2,3 and LC) */
    int           amb_solvab[NFREQ+1];  /* ambiguity solvablility */
    double        amb_CP[NFREQ+1];      /* ambiguity calculated by Code-Phase */
    double        amb[NFREQ+1];         /* ambiguity solution vector (L1,2,3 and LC) */
    double        amb_ave[NFREQ+1];     /* average ambiguity vector since reset (L1,2,3 and LC) */
    double        ambvar[NFREQ+1];      /* ambiguity solution variance vector (L1,2,3 and LC) */
    double        init_var[NFREQ+1];    /* initial ambiguity variance (L1,2,3 and LC) */

    int           reset[NFREQ+1];       /* reset ambiguity flag */
    int           fix[NFREQ+1];         /* ambiguity fix flag (1:initial,2:fix,3:hold) for L1,2,3 and LC */
    double        d_ave[NFREQ+2];       /* current correction of average amb (L1,2,3 and LC,mw12) */
    gtime_c       phase_time[NFREQ+1];  /* time of last observed phase (L1,2,3 and LC) */
    gtime_c       solv_time[NFREQ+1];   /* time of last ambiguity solution (L1,2,3 and LC) */
    double        phase_dt[NFREQ+1];    /* time difference between two consecutive observed phases (L1,2,3 and LC) */
    double        lock_period[NFREQ+1]; /* lock time of phase (L1,2,3 and LC) */
    double        solv_period[NFREQ+1]; /* continuously solved time of phase (L1,2,3 and LC) */
    int           lock_con[NFREQ+1];    /* lock count of phase (L1,2,3 and LC) */
    int           cobs_con[NFREQ+1];    /* continuously observed count of phase (L1,2,3 and LC) */
    int           solv_con[NFREQ+1];    /* solved count of ambiguity (L1,2,3 and LC) */
    int           ave_con[NFREQ+1];     /* average count of ambiguity (L1,2,3 and LC) */
    int           MW12_con,MW15_con;    /* averaged count of MW combination */
    unsigned int  half[NFREQ];          /* half-cycle valid flag */

    double        multiPath[NFREQ+1];   /* multipath effects */
};

/* RTK server class ------------------------------------------------------------------------------- */
class gnss_rtksvr_c {
/* Constructor */
public:
    gnss_rtksvr_c();
    ~gnss_rtksvr_c();
/* Implementation functions */
protected:
    /* initialzie function -------------------------------------------------------- */
    /* set sat\rcv antenna information -------------------------------------------- */
    void setpcv(in_gnss_atx_c *atx,int satrcv);
    /* new gnss_pro according to opt ---------------------------------------------- */
    void ini_gnss_pro(gnss_prcopt_c *Prcopt,gnss_filopt_t *Filopt);
    /* initialize stream environment ---------------------------------------------- */
    void strinitcom();
    /* initialize decode format --------------------------------------------------- */
    void inidecode();
    /* initialize stream class ---------------------------------------------------- */
    void inistream();

    /* process function ----------------------------------------------------------- */
    /* sync input streams (if type=STR_FILE) -------------------------------------- */
    void strsync();
    /* write solution header to output stream ------------------------------------- */
    void writesolhead();
    /* update glonass frequency channel number in raw data struct ----------------- */
    void updatefcn();
    /* write solution to each out-stream (stream[3:4])----------------------------- */
    void writesolstr(int index);
public:
    /* initialize observation pointer obsr/obsb (*gnss_pro) ----------------------- */
    int iniobs();
    /* lock/unlock rtk server ----------------------------------------------------- */
    void rtksvrlock();
    void rtksvrunlock();
    /* write solution to output stream -------------------------------------------- */
    void writesol();
    /* input message from stream -------------------------------------------------- */
    /* update rtk server struct --------------------------------------------------- */
    void updatesvr(int ret,int index);
    /* decode receiver raw/rtcm data ---------------------------------------------- */
    int decoderaw(int index);
    /* initialize rtksvr ---------------------------------------------------------- */
    int rtksvrini(all_option_c *option);
    /* start rtksvr --------------------------------------------------------------- */
    int rtksvrstart();
    /* stop rtksvr ---------------------------------------------------------------- */
    void rtksvrstop(char **cmds);

/* Components */
public:
    int state;                          /* server state (0:stop,1:running) */
    int sampling;                       /* observation sampling time */
    int cyctime;                        /* processing cycle (ms) */
    int nmeacycle;                      /* NMEA request cycle (ms) (0:no req) */
    int nmeareq;                        /* NMEA request (0:no,1:nmeapos,2:single sol) */
    double nmeapos[3];                  /* NMEA request position (ecef) (m) */
    int buffsize;                       /* input buffer size (bytes) */
    int format[3];                      /* input format {rov,base,corr} */
    gnss_solopt_c solopt[2];            /* output solution options {sol1,sol2} */
    gnss_prcopt_c *prcopt;              /* GNSS processing options */

    int navsel;                         /* ephemeris select (0:all,1:rover,2:base,3:corr) */
    int nsbs;                           /* number of sbas message */
    int nsol;                           /* number of solution buffer */
    gnss_pro_c *gnss_pro;               /* GNSS process class */
    int nb[3];                          /* bytes in input buffers {rov,base} */
    int nsb[2];                         /* bytes in soulution buffers */
    int npb[3];                         /* bytes in input peek buffers */
    unsigned char *buff[3];             /* input buffers {rov,base,corr} */
    unsigned char *sbuf[2];             /* output buffers {sol1,sol2} */
    unsigned char *pbuf[3];             /* peek buffers {rov,base,corr} */
    unsigned int nmsg[3][10];           /* input message counts */
    decode_data_c *data[3];             /* un-decoded data (raw,rtcm2,rtcm3) for {rov,base,corr} */
    gtime_c ftime[3];                   /* download time {rov,base,corr} */
    string files[3];                    /* download paths {rov,base,corr} */
    gnss_obs_c obs[3];                  /* observation data {rov,base,corr} (give to gnss_pro) */
    gnss_nav_c *nav;                    /* navigation data */
    gnss_sbsmsg_c sbsmsg[MAXSBSMSG];    /* SBAS message buffer */
    int strtype[8];                     /* stream types {rov,base,corr,sol1,sol2,logr,logb,logc} */
    stream_c *stream[8];                /* streams {rov,base,corr,sol1,sol2,logr,logb,logc} */
    stream_c *moni;                     /* monitor stream */
    unsigned int tick;                  /* start tick */
    thread_t thread;                    /* server thread */
    int cputime;                        /* CPU time (ms) for a processing cycle */
    int prcout;                         /* missing observation data count */
    int nave;                           /* number of averaging base pos */
    double rb_ave[3];                   /* averaging base pos */
    string cmds_periodic[3];            /* periodic commands */
    lock_t lock_f;                      /* lock flag */
    int fobs[3];                        /* observation buff number */
    string errorMsg;                    /* error message */
};

/* post-processing server class ------------------------------------------------------------------- */
class gnss_postsvr_c {
    /* Constructor */
public:
    gnss_postsvr_c();
    ~gnss_postsvr_c();
    /* Implementation functions */
protected:
    /* set sat\rcv antenna information -------------------------------------------- */
    void setpcv(in_gnss_atx_c *atx,int satrcv);
    /* new gnss_pro according to opt ---------------------------------------------- */
    void ini_gnss_pro(gnss_prcopt_c *Prcopt,gnss_filopt_t *Filopt);
    /* update station information to prcopt --------------------------------------- */
    void updatesta();

    /* initialize read stream ----------------------------------------------------- */
    int ini_Read_Stream();
    /* rover/base observation synchronization ------------------------------------- */
    int obs_synchron();

    /* write solution header to output stream ------------------------------------- */
    void writesolhead();
    /* write solution to output stream -------------------------------------------- */
    void writesol();
    /* write solution to each solstream ------------------------------------------- */
    void writesolstr(int index);

    /* read observation file from opbsBuff ---------------------------------------- */
    int Read_One_Epoch(int rovbas);
public:
    /* initialize postsvr --------------------------------------------------------- */
    int postsvrini(all_option_c *option);
    /* reset postsvr -------------------------------------------------------------- */
    void postreset();
    /* read Navigation file ------------------------------------------------------- */
    int readNav();
    /* post-position epoch by epoch ----------------------------------------------- */
    int Post_Position_Epoch();
    /* Components */
protected:
    vector<gnss_obs_c> obsBuff[2];      /* observation data buff {0:rover,1:base} */
    int obsBuffIndex[2];                /* observation data buff index {0:rover,1:base} */
    gnss_sta_c rbStaInf[2];             /* rover and base station information {0:rover,1:base} */
    gnss_obs_c *obs[2];                 /* observation data {0:rover,1:base} */
    gnss_obs_c epochObs[2];             /* observation data of one Epoch */
    file_c solstream[2];                /* solution output stream */
    string sat_ant;                     /* satellite antenna file path */
    string rec_ant;                     /* receiver antenna file path */
    string blq_sta;                     /* station ocean-loading tide file path */
    string snx_sta;                     /* station SNX solution file path */
    in_gnss_rnxO_c rover,base;          /* read observation file stream */
    in_gnss_rnxN_c readN;               /* read navigation file stream */
    in_gnss_eph_c readEph;              /* read precise ephemeris stream */
    in_gnss_rnxC_c readC;               /* read precise clock file stream */
public:
    gnss_pro_c *gnss_pro;               /* GNSS process class */
    /* GNSS options */
    gnss_prcopt_c *prcopt;              /* GNSS processing options */
    gnss_pstopt_c *pstopt;              /* GNSS post-procssing options */
    gnss_solopt_c solopt[2];            /* GNSS solution options */
    /* GNSS data */
    gnss_nav_c *nav;                    /* navigation data */
};

/* ROTI process class ----------------------------------------------------------------------------- */
/* single station ROTI class ---------------------------------------------------------------------- */
class gnss_rotipdc_c {
    /* Constructors */
public:
    gnss_rotipdc_c();
    virtual ~gnss_rotipdc_c();
    /* Implementation functions */
protected:
    /* new gnss_pro according to opt ---------------------------------------------- */
    void ini_gnss_pro(gnss_prcopt_c *Prcopt,gnss_filopt_t *Filopt);
    /* update station information to prcopt --------------------------------------- */
    void updatesta();
    /* write head to output file -------------------------------------------------- */
    void write_head_ROTI();

    /* initialize read stream ----------------------------------------------------- */
    int ini_Read_Stream();

    /* verify excluded satellite -------------------------------------------------- */
    int exclude_sat(gnss_obsd_c &obs);
    /* initialize station postion ------------------------------------------------- */
    int initposition();
    /* preprocess of observations ------------------------------------------------- */
    int getSTEC();
    /* write ROTI information to log_stream --------------------------------------- */
    void write_ROTI();
public:
    /* initialize gnss_rotipdc_c ------------------------------------------------------- */
    int rotipdcini(all_option_c *option);
    /* compute ROTI --------------------------------------------------------------- */
    int getROTI();

    /* Components */
protected:
    gnss_pro_c *gnss_pro;               /* GNSS process class */   
    vector<double> STEC[MAXSAT];        /* STEC list of a period for every satellite (TECU) */
    vector<int> CSFlag[MAXSAT];         /* cycle-slip flag of a period for every satellite 
                                        * (0:no slip, 1:slip) */
    vector<gtime_c> STECTime;           /* STEC time list of a period */
    vector<double> dT;                  /* time difference list between two adjacent epochs (s) */
    vector<double> dT2LastEpoch;        /* time difference list between epochs and last epoch (s) */
    vector<double> ROT[MAXSAT];         /* ROT list of a period for every satellite (TECU/min) */
    int numROT[MAXSAT];                 /* ROT number list of a period for every satellite */
    double ROTI[MAXSAT];                /* ROTI of one period (TECU/min) */

    in_gnss_rnxO_c rover;               /* read observation file stream */
    in_gnss_rnxN_c readN;               /* read navigation file stream */
    in_gnss_eph_c readEph;              /* read precise ephemeris stream */
    in_gnss_rnxC_c readC;               /* read precise clock file stream */

    gnss_sta_c staInf;                  /* rover station information */
    gnss_obs_c epochObs;                /* observation data of one Epoch */

    int nEpoch;                         /* epochs number in the current period window */
    int allEpoch;                       /* all epochs number */

public:
    /* GNSS options */
    gnss_prcopt_c *prcopt;              /* GNSS processing options */

    /* GNSS data */
    gnss_nav_c *nav;                    /* navigation data */
    int epochNuminOnePeriod;            /* epoch number of one period */
    int errorCode;                      /* error code */
    string errorMessage;                /* error message */
};

#endif
