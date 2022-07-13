/* Class of Reading Positioning Configure file */

#ifndef CONFIG_H
#define CONFIG_H

#include "gnssnav_lib.h"
#include "BaseFunction/timesys.h"
#include "GNSS/DataClass/data.h"
#include "RtkStream/stream.h"

/* GNSS configuration ----------------------------------------------------------------------------- */
/* struct: option structure ----------------------------------------------------------------------- */
typedef struct {
    const string name;                  /* option name */
    int format;                         /* option format (0:int,1:double,2:string,3:enum) */
    void *var;                          /* pointer to option variable */
    const string comment;               /* option comment/enum labels/unit */
} opt_t;

/* post-processing option class ------------------------------------------------------------------- */
class gnss_pstopt_c {
/* Constructor */
public:
    gnss_pstopt_c();
    ~gnss_pstopt_c();
/* Implementation functions */
protected:
public:
    /* set time parameters -------------------------------------------------------- */
    int set_time_par();
/* Components */
public:
    int predict;                        /* flag of use predicted data */
    double time_inter;                  /* time interval (s) */
    gtime_c time_start,time_end;        /* start and end time */
    vector<gtime_c> time_skip[2];       /* skipped time vector
                                         * [0]: start time 
                                         * [1]: end time */
    int nskip_period;                   /* number of skipped period */
    string rover_obs;                   /* rover observation file */
    string base_obs;                    /* base observation file */
    string nav;                         /* navigation file */
    string prseph;                      /* precise ephemeris file */
    string prsclk;                      /* precise clock file */
    string output[2];                   /* out-put file */
};

/* rtk-processing option class -------------------------------------------------------------------- */
class gnss_rtkopt_c {
/* Constructor */
public:
    gnss_rtkopt_c();
    ~gnss_rtkopt_c();
/* Components */
public:
    int strtype[8];                     /* stream types */
    string strpath[8];                  /* stream paths */
    int strfmt[3];                      /* stream formats */

    int svrcycle;                       /* server cycle (ms) */
    int timeout;                        /* timeout time (ms) */
    int reconnect;                      /* reconnect interval (ms) */
    int nmeacycle;                      /* nmea request cycle (ms) */
    int fswapmargin;                    /* file swap margin (s) */
    int buffsize;                       /* input buffer size (bytes) */
    int navmsgsel;                      /* navigation message select */
    int nmeareq;                        /* nmea request type (0:off,1:lat/lon,2:single) */
    double nmeapos[3];                  /* nmea position (lat/lon) (deg) */
    string proxyaddr;                   /* proxy address (1024) */

    /* unset options */
    string cmds[3];                     /* stream start commands (256) */
    string rropts[3];                   /* receiver and rtcm options (256) */
    stream_c *monitor;                  /* monitor stream */
};

/* processing options class ----------------------------------------------------------------------- */
class gnss_prcopt_c {
/* Constructor */
public:
    gnss_prcopt_c();
    ~gnss_prcopt_c();
/* Implementation functions */
public:
/* Components */
public:
    int mode;                           /* positioning mode (PMODE_???) */
    int pppmode;                        /* PPP mode */
    int soltype;                        /* solution type (0:forward,1:backward,2:combined) */
    int nf;                             /* number of frequencies (1:F1,2:F1+F2,3:F1+F2+F3) */
    int usedF,posF;                     /* number of used (CS detection) and positioning frequencies */
    int freqPriority;                   /* 0:off,1:on, set used frequencies according to priority */
    int navsys;                         /* navigation system */
    double min_elevation;               /* elevation mask angle (rad) */
    double con_weight_el;            /* minimum elevation to const weight for pseudorange (rad) */
    snrmask_t snrmask;                  /* SNR mask */
    int sateph;                         /* satellite ephemeris/clock (EPHOPT_???) */
    int logmsg;                         /* pdop output control */
    int bdsGen;                         /* BDS-2 and BDS-3 satellites control */
    int bdsOrbit;                       /* BDS orbit (1:MEO + 2:IGSO + 4:GEO) */
    int freq_glo[MAX_NF];               /* GLONASS freqeucny priority (e.g.: 1,2,3,4,6) */
    int freq_gal[MAX_NF];               /* Galileo freqeucny priority (e.g.: 1,5,6,7,8) */
    int freq_bds3[MAX_NF];              /* BDS-3 freqeucny priority (e.g.: 1,2,6,7,5,8) */
    int freq_bds2[MAX_NF];              /* BDS-2 freqeucny priority (e.g.: 1,2,6,7,5,8) */

    int modear;                         /* AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold,
                                         * 4:ppp-ar) */
    int order;                          /* order of polynomial fitting (L) to detect cycle slip (3-5) */
    int slipmode;                       /* slip-detect mode (0:obs + 1:CP + 2:geo-free + 4:Melbourne-Wubbena) */
    double slip_std;                    /* m, std variance of repaired cycle slip (1-10) */
    double cs_gf;                       /* cycle-slip detected threshold of geometry-free phase (m) */
    double cs_mw;                       /* cycle-slip detected threshold of Melbourne-Wubbena (m) */
    double sampling;                    /* sec, observation time sampling */
    int minConsecutive;                 /* minimum consecutive number to use satellite when samping <= 15 */
    double restime;                     /* min lock time gap to reset ambiguity (s) */
    int glomodear;                      /* GLONASS AR mode (0:off,1:on,2:auto cal,3:ext cal) */
    int bdsmodear;                      /* BeiDou AR mode (0:off,1:on) */
    int iniamb;                         /* initialize ambiguity if lock count < iniar */
    int maxariter;                      /* max iteration to resolve ambiguity */

    int sppiono;                        /* SPP ionosphere option (no estimate option) */
    int ionoopt;                        /* ionosphere option (IONOOPT_???) */

    int spptrop;                        /* SPP troposphere option (no estimate option) */
    int tropopt;                        /* troposphere option (TROPOPT_???) */

    int dynamics;                       /* dynamics model (0:none,1:velociy) */
    int tidecorr;                       /* earth tide correction (0:off,1:solid,2:solid+otl+pole) */
    int adjustfunc;                     /* adjustment function */
    int adjRobust;                      /* robust mode for adjustment (0:none,1:IGG3-2,2:IGG3-3) */
    int robustRtype;                    /* covariance type in robust function (0:var,1:covar,2:vector) */
    int niter;                          /* number of filter iteration */
    int codesmooth;                     /* code smoothing window size (0:none) */
    int intpref;                        /* interpolate reference obs (for post mission) */
    int sbascorr;                       /* SBAS correction options */
    int sbassatsel;                     /* SBAS satellite selection (0:all) */
    int rovpos;                         /* rover position for fixed mode */
    int baspos;                         /* base position for relative mode */
                                        /* (0:pos in prcopt,  1:average of single pos, */
                                        /*  2:read from file, 3:rinex header, 4:rtcm pos) */
    double pppfactor;                   /* PPP observation noise factor 
                                         * "ppp_noise (m) = pppfactor x obs_noise (m)" */
    double err[6];                      /* measurement error factor */
                                        /* [0]: code error */
                                        /* [1]: elevation dependent code error */
                                        /* [2]: phase error */
                                        /* [3]: elevation dependent phase error */
                                        /* [4]: baseline-length dependent phase error m/10km */
                                        /* [5]: doppler frequency (hz) */

    int ISBmode;                        /* receiver ISB mode (0:whiteNoise,1:randomWalk) */
    double std[3];                      /* initial-state std [0]ambiguity,[1]iono [2]trop */
    double stdrate[4];                  /* growth rate of std [0]ambiguity [1]iono,[2]trop
                                         * [3]ISB */
    double sclkstab;                    /* satellite clock stability (sec/sec) */
    double thresar[8];                  /* AR validation threshold */
    double minElArFix,minElArHold;      /* elevation mask to fix/hold ambiguity (deg) */
    double maxtdiff;                    /* max difference of time (sec) */
    double maxres;                      /* max residual of code in SPP */
    double baseline[2];                 /* baseline length constraint {const,sigma} (m) */
    int baseAtt;                        /* estimate attitude according to solved baseline (0:off,1:on) */

    string name[2];                     /* rover and base name */
    double ru[3];                       /* rover position for fixed mode {x,y,z} (ecef) (m) */
    double rb[3];                       /* base position for relative mode {x,y,z} (ecef) (m) */
    string anttype[2];                  /* antenna types {rover,base} */
    double antdel[2][3];                /* antenna delta {{rov_e,rov_n,rov_u},{ref_e,ref_n,ref_u}} */
    gnss_pcv_c pcvr[2];                 /* receiver antenna parameters {rov,base} */
    unsigned char exsats[MAXSAT];       /* excluded satellites (1:excluded,2:included) */
    int  maxaveep;                      /* max averaging epoches */
    int  initrst;                       /* initialize by restart */
    string rnxopt[2];                   /* rinex options {rover,base} */
    int  posopt[3];                     /* positioning options */
    int  syncsol;                       /* solution sync mode (0:off,1:on) */
    int freqopt;                        /* disable L2-AR */
    string pppopt;                      /* ppp option */
};

/* solution options class ------------------------------------------------------------------------- */
class gnss_solopt_c {
/* Constructors */
public:
    gnss_solopt_c();
    ~gnss_solopt_c();
/* Implementation functions */
public:
    /* write solution header to output stream --------------------------------- */
    int outsolheads(const gnss_prcopt_c *prcopt,unsigned char *buff);
/* Componentss */
public:
    int posf;                           /* solution format (SOLF_???) */
    int refsta;                         /* eun reference station (0:base,1:rover) */
    int times;                          /* time system (TIMES_???) */
    int timef;                          /* time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s) */
    int timeu;                          /* time digits under decimal point */
    int degf;                           /* latitude/longitude format (0:ddd.ddd,1:ddd mm ss) */
    int outhead;                        /* output header (0:no,1:yes) */
    int outopt;                         /* output processing options (0:no,1:yes) */
    int datum;                          /* datum (0:WGS84,1:CGCS2000) */
    int height;                         /* height (0:ellipsoidal,1:geodetic) */
    int geoid;                          /* geoid model (0:EGM96,1:JGD2000) */
    int solstatic;                      /* solution of static mode (0:all,1:single) */
    int sstat;                          /* solution statistics level (0:off,1:states,2:residuals) */
    int trace;                          /* debug trace level (0:off,1-5:debug) */
    double nmeaintv[2];                 /* nmea output interval (s) (<0:no,0:all) */
                                        /* nmeaintv[0]:gprmc,gpgga,nmeaintv[1]:gpgsv */
    
    string sep;                         /* field separator */
    string prog;                        /* program name */
};
/* file options class ----------------------------------------------------------------------------- */
class gnss_filopt_t {
/* Constructor */
public:
    gnss_filopt_t();
    ~gnss_filopt_t();
/* Components */
public:
    string satantp;                     /* satellite antenna parameters file */
    string rcvantp;                     /* receiver antenna parameters file */
    string stapos;                      /* station positions file */
    string geoid;                       /* external geoid data file */
    string iono;                        /* ionosphere data file */
    string ztd;                         /* troposphere ztd file */
    string dcb;                         /* dcb data file */
    string erp;                         /* erp data file */
    string blq;                         /* ocean loading tide blq file */
    string snx;                         /* snx file */
    string tempdir;                     /* ftp/http temporaly directory */
    string geexe;                       /* google earth exec file */
    string solstat;                     /* solution statistics file */
    string test;                        /* debug test file */
};
/* RINEX options class ---------------------------------------------------------------------------- */
class gnss_rnxopt_c {
/* Constructor */
public:
    gnss_rnxopt_c();
    ~gnss_rnxopt_c();
/* Components */
public:
    gtime_c ts,te;                      /* time start/end */
    double tint;                        /* time interval (s) */
    double tunit;                       /* time unit for multiple-session (s) */
    double rnxver;                      /* RINEX version */
    int navsys;                         /* navigation system */
    int obstype;                        /* observation type */
    int freqtype;                       /* frequency type */
    char mask[7][64];                   /* code mask {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    char staid[32];                     /* station id for rinex file name */
    char prog[32];                      /* program */
    char runby[32];                     /* run-by */
    char marker[64];                    /* marker name */
    char markerno[32];                  /* marker number */
    char markertype[32];                /* marker type (ver.3) */
    char name[2][32];                   /* observer/agency */
    char rec[3][32];                    /* receiver #/type/vers */
    char ant[3][32];                    /* antenna #/type */
    double apppos[3];                   /* approx position x/y/z */
    double antdel[3];                   /* antenna delta h/e/n */
    string comment[MAXCOMMENT];         /* comments */
    char rcvopt[256];                   /* receiver dependent options */
    unsigned char exsats[MAXSAT];       /* excluded satellites */
    int scanobs;                        /* scan obs types */
    int outiono;                        /* output iono correction */
    int outtime;                        /* output time system correction */
    int outleaps;                       /* output leap seconds */
    int autopos;                        /* auto approx position */
    int halfcyc;                        /* half cycle correction */
    gtime_c tstart;                     /* first obs time */
    gtime_c tend;                       /* last obs time */
    gtime_c trtcm;                      /* approx log start time for rtcm */
    string tobs[7][MAXOBSTYPE];         /* obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    int nobs[7];                        /* number of obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
};

/* all options ------------------------------------------------------------------------------------ */
class all_option_c {
/* Constructor */
public:
    all_option_c();
    ~all_option_c();
/* Implementation functions */
protected:
    /* reset system options to default ---------------------------------------- */
    void resetsysopts();
    /* string option to enum (int) -------------------------------------------- */
    int str2enum(const string str,const string comment,int *val);
    /* enum (int) to string option -------------------------------------------- */
    int enum2str(string &str,const string comment,int val);
    /* discard space characters at tail --------------------------------------- */
    void chop(string &str);
    /* search option ---------------------------------------------------------- */
    opt_t *searchopt(const string name,const opt_t *opts);
    /* string to option value ------------------------------------------------- */
    int str2opt(opt_t *opt,const string str);
    /* load options ----------------------------------------------------------- */
    int loadopts(const string file,opt_t *opts);
    /* system options buffer to options --------------------------------------- */
    void buff2sysopts();
    /* get system options ----------------------------------------------------- */
    void getsysopts();
public:
    /* GNSS options ----------------------------------------------------------- */
    /* update input and output file options ----------------------------------- */
    void updatefile(const int days);
    /* read post options ------------------------------------------------------ */
    int readpostopt(const string file);
    /* read rtk options ------------------------------------------------------- */
    int readrtkopt(const string file);

    /* Integration options ---------------------------------------------------- */
    /* read all post processing options --------------------------------------- */
    int readpostAll(const string file);
    /* read all rtk processing options ---------------------------------------- */
    int readrtkAll(const string file);

/* Components */
public:
    /* GNSS options */
    gnss_pstopt_c pstopt;
    gnss_rtkopt_c rtkopt;
    gnss_prcopt_c prcopt;
    gnss_solopt_c solopt[2];
    gnss_filopt_t filopt;
    gnss_rnxopt_c rnxopt;
};
#endif
