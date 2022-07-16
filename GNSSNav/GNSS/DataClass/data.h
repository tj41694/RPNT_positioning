/* Data Classes File */

#ifndef DATA_H
#define DATA_H

#include "BaseFunction/timesys.h"

/* class of frequency priority and adoption ------------------------------------------------------- */
class gnss_freq_c {
/* Constructors */
public:
    gnss_freq_c();
    ~gnss_freq_c();
/* Implementation functions */
public:
    /* initialize the GNSS frequency priority with default value ------------------ */
    void init_gnss_freq_priority();
    /* initialize the GNSS used frequency with default value ---------------------- */
    void init_gnss_used_freq();
    /* copy GNSS_USED_FREQ to external matrix ------------------------------------- */
    void copy_used_frep(int ext_freq[MAX_GNSS][NFREQ*2]);

    /* change the GNSS frequency priority with configuration ---------------------- */
    void change_freq_priority(const gnss_prcopt_c *opt);
    /* change the used GNSS frequency according to priority and availability ------ */
    void change_gnss_freq(const int sys,sigind_t *indRov,sigind_t *indBas);
/* Components */
public:
    /* change flag */
    int priority_glo;                   /* GLONASS frequency priority changed flag */
    int priority_gal;                   /* Galileo frequency priority changed flag */
    int priority_bds3;                  /* BDS-3 frequency priority changed flag */
    int priority_bds2;                  /* BDS-2 frequency priority changed flag */

    /* BDS index */
    int BDS_B1_index[2][2];             /* [x][0]:C1(B1A/B1C), [x][1]:C2(B1-2) index for 
                                         * [0][x]:BDS-3 and [1][x]:BDS-2 in GNSS_FREQ_PRI */

    /* gnss frequency priority -------------------------------------------------------
    * the 2nd number of "MAX_NF" values are only set for BDS-2
    * ( GPS + GLO + GAL + BDS + QZS + IRN + LEO ) --------------------------------- */
    int GNSS_FREQ_PRI[MAX_GNSS][MAX_NF*2] = {
        // GPS
        { 1, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // GLONASS
        { 1, 2, 3, 4, 6, 0, 0, 0, 0, 0, 0, 0 },
        // Galileo
        { 1, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0 },
        // BDS
        { 2, 6, 7, 5, 8, 1, 2, 6, 7, 5, 8, 1 }, /* 0:5 for BDS-3 and 6:11 for BDS-2 */
        // QZSS
        { 1, 2, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0 },
        // IRNSS
        { 5, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // SBAS
        { 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    /* used GNSS frequency -----------------------------------------------------------
    * only changed in function copy_used_frep(int ext_freq[MAX_GNSS][NFREQ*2])!!!!!!!!
    * the 2nd number of "MAX_NF" values are only set for BDS-2
    * ( GPS + GLO + GAL + BDS + QZS + IRN + LEO ) --------------------------------- */
    int GNSS_USED_FREQ[MAX_GNSS][NFREQ*2] = {
        // GPS
        { 1, 2, 5, 0, 0, 0 },
        // GLONASS
        { 1, 2, 3, 0, 0, 0 },
        // Galileo
        { 1, 5, 6, 0, 0, 0 },
        // BDS
        { 2, 6, 7, 2, 6, 7 }, /* 0:2 for BDS-3 and 3:5 for BDS-2 */
        // QZSS
        { 1, 2, 5, 0, 0, 0 },
        // IRNSS
        { 5, 9, 0, 0, 0, 0 },
        // SBAS
        { 1, 5, 0, 0, 0, 0 }
    };
};

/* class of one epoch observation data ------------------------------------------------------------ */
class gnss_obsd_c {
/* Constructors */
public:
    gnss_obsd_c();
    ~gnss_obsd_c();
/* Implementation functions */
public:
    /* reset the whole data ------------------------------------------------------- */
    void reset();
    /* reset satellite data ------------------------------------------------------- */
    void satreset();
    /* update signal time use pseudorange------------------------------------------ */
    int sigtime_opsr();
    /* update signal time use satellite clock bias -------------------------------- */
    void sigtime_clk();
    /* update signal time use known position and receiver clock bias -------------- */
    void sigtime_rec(const vector<double> &pos);
/* Components */
public:
    unsigned int rcv;                   /* receiver number 0:rover, 1:base */
    unsigned int sat;                   /* satellite number */
    int prn,sys,isys;                   /* satellite prn and system number/index */
    int nCon;                           /* number of consecutive */
    gtime_c time;                       /* receiver sampling time (GPST) */
    gtime_c sigtime;                    /* satellite send signal time (GPST) */
    double dtr;                         /* receiver clock bias */
    double SNR[NFREQ+NEXOBS];           /* signal strength (BHz) */
    unsigned char LLI[NFREQ+NEXOBS];    /* loss of lock indicator */
    unsigned char code[NFREQ+NEXOBS];   /* code indicator (CODE_???) */
    double L[NFREQ+NEXOBS];             /* observation data carrier-phase (cycle) */
    double P[NFREQ+NEXOBS];             /* observation data pseudorange (m) */
    double D[NFREQ+NEXOBS];             /* observation data doppler frequency (Hz) */
    string LCode[NFREQ+NEXOBS];         /* observation code carrier-phase (cycle) */
    string PCode[NFREQ+NEXOBS];         /* observation code pseudorange (m) */
    string DCode[NFREQ+NEXOBS];         /* observation code doppler frequency (Hz) */
    string SCode[NFREQ+NEXOBS];         /* observation code signal strength (BHz) */
    double DCB[NFREQ+NEXOBS];           /* pseudorange DCB correction */
    double res[NFREQ+NEXOBS];           /* observation residual */
    double ovar[NFREQ+NEXOBS];          /* observation variance */
    double dist;                        /* distance between sat and rec */
    double sigvec[3];                   /* line-of-sight vector */
    double dcbvar;                      /* code bias error variance */
    int used;                           /* observation used flag */
    int vused;                          /* observation used flag for velocity */
    int exc;                            /* observation excluded flag */
    double posvel[6];                   /* satellite position and velocity (ecef) */
    double dts[2];                      /* satellite clock bias and drift */
    double anto[NFREQ][3];              /* satellite antenna phase center offset correction (xyz) */
    double antv[NFREQ];                 /* satellite antenna phase center variance correction */
    double sant[NFREQ+1];               /* satellite antenna phase center correction (L1,2,3 and LC) */
    double azel[2];                     /* satellite azimuth and elevation angles */
    double svar;                        /* satellite position and clock variance */
    int svh;                            /* satellite health flag */
    double ionmap,dion,ionvar;          /* ionosphere delay mapping function,
                                         * correction and variance (GPS L1, m) */
    double ipp[2];                      /* ionospheric pierce point position (lat, lon) */
    double dtro,trovar;                 /* troposphere delay correction and variance */
    double m_h,m_w;                     /* troposphere mapping function */
    double rant[NFREQ+1];               /* receiver antenna phase center correction (L1,2,3 and LC) */
    double phasewp;                     /* phase windup */
    string warningMsg,errorMsg;         /* warning & error message */
};

/* station informtations -------------------------------------------------------------------------- */
class gnss_sta_c {
/* Constructors */
public:
    gnss_sta_c();
    ~gnss_sta_c();
/* Components */
public:
    string name;                        /* marker name */
    string marker;                      /* marker number */
    string antdes;                      /* antenna descriptor */
    string antsno;                      /* antenna serial number */
    string rectype;                     /* receiver type descriptor */
    string recver;                      /* receiver firmware version */
    string recsno;                      /* receiver serial number */
    int antsetup;                       /* antenna setup id */
    int itrf;                           /* ITRF realization year */
    int deltype;                        /* antenna delta type (0:enu,1:xyz) */
    double pos[3];                      /* station position (ecef) (m) */
    double del[3];                      /* antenna position delta (e/n/u or x/y/z) (m) */
    double hgt;                         /* antenna height (m) */
};

/* class of station's information and observation chains ------------------------------------------ */
class gnss_obs_c {
/* Constructors */
public:
    gnss_obs_c();
    ~gnss_obs_c();
    /* Implementaion functions */
public:
    /* reset all data ------------------------------------------------------------- */
    void reset();
/* Components */
public:
    int n;                              /* observation number */
    unsigned char rcv;                  /* receiver number 0:rover, 1:base */
    gnss_sta_c sta;                     /* station parameter type */
    vector<gnss_obsd_c> data;           /* observation chain */
    string errorMsg,warnMsg;            /* error/warning massage */
    int used;                           /* number of used observations */
};

/* BDS/GPS/QZS/GAL broadcast ephemeris class ------------------------------------------------------ */
class gnss_eph_c {
/* Constructors */
public:
    gnss_eph_c();
    ~gnss_eph_c();
/* Components */
public:
    /* parameters refer to the antenna phase center */
    int sat;                            /* satellite number */
    int prn;                            /* satellite PRN */
    int iode,iodc;                      /* IODE,IODC */
    int sva;                            /* SV accuracy (URA index) */
    int svh;                            /* SV health (0:ok) */
    int week;                           /* GPS/QZS: gps week, GAL: galileo week, BDS: BeiDou week */
    int code;                           /* GPS/QZS: code on L2, GAL/CMP: data sources */
    int flag;                           /* GPS/QZS: L2 P data flag, CMP: nav type */
    gtime_c toe,toc,ttr;                /* Toe,Toc,T_trans */
                                        /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,deln,OMGd,idot;
    double crc,crs,cuc,cus,cic,cis;
    double toes;                        /* Toe (s) in week */
    double fit;                         /* fit interval (h) */
    double f0,f1,f2;                    /* SV clock parameters (af0,af1,af2) */
    double tgd[4];                      /* group delay parameters */
                                        /* GPS/QZS:tgd[0]=TGD */
                                        /* GAL    :tgd[0]=BGD1 E5a/E1,tgd[1]=BGD2 E5b/E1 */
                                        /* CMP    :tgd[0]=TGD1 B1/B3, tgd[1]=TGD2 B2/B3 */
    double Arate;                       /* changing rate of A */
    double Adot,ndot;                   /* Adot,ndot for CNAV */
};

/* GLONASS broadcast ephemeris class -------------------------------------------------------------- */
class gnss_geph_c {
/* Constructors */
public:
    gnss_geph_c();
    ~gnss_geph_c();
/* Components */
public:
    int sat;                            /* satellite number */
    int prn;                            /* satellite PRN */
    int iode;                           /* IODE (0-6 bit of tb field) */
    int frq;                            /* satellite frequency number */
    int svh,sva,age;                    /* satellite health, accuracy, age of operation */
    gtime_c toe;                        /* epoch of epherides (gpst) */
    gtime_c tof;                        /* message frame time (gpst) */
    double pos[3];                      /* satellite position (ecef) (m) */
    double vel[3];                      /* satellite velocity (ecef) (m/s) */
    double acc[3];                      /* satellite acceleration (ecef) (m/s^2) */
    double taun,gamn;                   /* SV clock bias (s)/relative freq bias */
    double dtaun;                       /* delay between L1 and L2 (s) */
};

/* SBAS ephemeris class --------------------------------------------------------------------------- */
class gnss_seph_c {
/* Constructors */
public:
    gnss_seph_c();
    ~gnss_seph_c();
/* Components */
public:
    int sat;                            /* satellite number */
    int prn;                            /* satellite PRN */
    gtime_c t0;                         /* reference epoch time (GPST) */
    gtime_c tof;                        /* time of message frame (GPST) */
    int sva;                            /* SV accuracy (URA index) */
    int svh;                            /* SV health (0:ok) */
    double pos[3];                      /* satellite position (m) (ecef) */
    double vel[3];                      /* satellite velocity (m/s) (ecef) */
    double acc[3];                      /* satellite acceleration (m/s^2) (ecef) */
    double af0,af1;                     /* satellite clock-offset/drift (s,s/s) */
};

/* precise ephemeris class ------------------------------------------------------------------------ */
class gnss_peph_c {
/* Constructors */
public:
    gnss_peph_c();
    ~gnss_peph_c();
/* Components */
public:
    gtime_c time;                       /* time (GPST) */
    double pos[MAXSAT][4];              /* satellite position/clock (ecef) (m|s) */
    float  std[MAXSAT][4];              /* satellite position/clock std (m|s) */
    double vel[MAXSAT][4];              /* satellite velocity/clk-rate (m/s|s/s) */
    float  vst[MAXSAT][4];              /* satellite velocity/clk-rate std (m/s|s/s) */
    float  cov[MAXSAT][3];              /* satellite position covariance (m^2) */
    float  vco[MAXSAT][3];              /* satellite velocity covariance (m^2) */
};

/* precise clock class ---------------------------------------------------------------------------- */
class gnss_pclk_c {
/* Constructors */
public:
    gnss_pclk_c();
    ~gnss_pclk_c();
/* Components */
public:
    gtime_c time;                       /* time (GPST) */
    double clk[MAXSAT];                 /* satellite clock (s) */
    float  std[MAXSAT];                 /* satellite clock std (s) */
};

/* almanac class ---------------------------------------------------------------------------------- */
class gnss_alm_c {
/* Constructors */
public:
    gnss_alm_c();
    ~gnss_alm_c();
/* Components */
public:
    int sat;                            /* satellite number */
    int svh;                            /* sv health (0:ok) */
    int svconf;                         /* as and sv config */
    int week;                           /* GPS/QZS: gps week, GAL: galileo week */
    gtime_c toa;                        /* Toa */
                                        /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,OMGd;
    double toas;                        /* Toa (s) in week */
    double f0,f1;                       /* SV clock parameters (af0,af1) */
};

/* TEC grid class --------------------------------------------------------------------------------- */
class gnss_tec_c {
/* Constructors */
public:
    gnss_tec_c();
    ~gnss_tec_c();
/* Components */
public:
    gtime_c time;                       /* epoch time (GPST) */
    int ndata[3];                       /* TEC grid data size {nlat,nlon,nhgt} */
    double rb;                          /* earth radius (km) */
    double lats[3];                     /* latitude start/interval (deg) */
    double lons[3];                     /* longitude start/interval (deg) */
    double hgts[3];                     /* heights start/interval (km) */
    vector<double> data;                /* TEC grid data (tecu) */
    vector<float> rms;                  /* RMS values (tecu) */
};

/* GNSS observation bias class -------------------------------------------------------------------- */
class gnss_bias_c {
    /* Constructors */
public:
    gnss_bias_c();
    ~gnss_bias_c();
    /* Components */
public:
    gtime_c cbts,cbte;                  /* satellite absolute code bias (OSB) start and end time */
    double P1P2[MAXSAT];                /* P1-P2 DCB from ionex file */
    double cbias[MAXSAT][MAXFREQ+1];    /* satellite absolute code bias (OSB) (0:none,1:F1,2:F2,...9:F9) (m) */
    double rbias[2][3];                 /* receiver dcb (0:p1-p2,1:p1-c1,2:p2-c2) (m) */
    double wlbias[MAXSAT];              /* wide-lane bias (cycle) */
    double glo_cpbias[4];               /* glonass code-phase bias {1C,1P,2C,2P} (m) */
};

/* ZTD data class --------------------------------------------------------------------------------- */
class gnss_ztd_c {
/* Constructors */
public:
    gnss_ztd_c();
    ~gnss_ztd_c();
/* Components */
public:
    gtime_c time;                       /* epoch time (GPST) */
    double ztd,rms;                     /* ZTD data and RMS (m) */
};

/* satellite fcb data class ----------------------------------------------------------------------- */
class gnss_fcbd_c {
/* Constructors */
public:
    gnss_fcbd_c();
    ~gnss_fcbd_c();
/* Components */
public:
    gtime_c ts,te;                      /* start/end time (GPST) */
    double bias[MAXSAT][3];             /* fcb value   (cyc) */
    double std[MAXSAT][3];              /* fcb std-dev (cyc) */
};

/* earth rotation parameter data class ------------------------------------------------------------ */
class gnss_erpd_c {
/* Constructors */
public:
    gnss_erpd_c();
    ~gnss_erpd_c();
/* Components */
public:
    double mjd;                         /* mjd (days) */
    double xp,yp;                       /* pole offset (rad) */
    double xpr,ypr;                     /* pole offset rate (rad/day) */
    double ut1_utc;                     /* ut1-utc (s) */
    double lod;                         /* length of day (s/day) */
};

/* earth rotation parameter class ----------------------------------------------------------------- */
class gnss_erp_c {
/* Constructors */
public:
    gnss_erp_c();
    ~gnss_erp_c();
/* Components */
public:
    int n,nmax;                         /* number and max number of data */
    vector<gnss_erpd_c> data;           /* earth rotation parameter data */
};

/* antenna parameter class ------------------------------------------------------------------------ */
class gnss_pcv_c {
/* Constructors */
public:
    gnss_pcv_c();
    ~gnss_pcv_c();
/* Components */
public:
    int sat;                            /* satellite number (0:receiver) */
    int prn;                            /* satellite prn number */
    string type;                        /* antenna type */
    string code;                        /* serial number or satellite code */
    gtime_c ts,te;                      /* valid time start and end */
    double off[NSYS][NFREQ][3];         /* phase center offset e/n/u or x/y/z (m)
                                         * NSYS systems for receiver
                                         * only 0 is used for satellite */
    double zen[3];                      /* zenith (ZEN1 ZEN2 DZEN) */
    int nzen;                           /* nzen = (ZEN2-ZEN1)/DZEN+1 */
    vector<double> var[NSYS][NFREQ];    /* phase center variation (m) */
                                        /* el=90,85,...,0 or nadir=0,1,2,3,... (deg)
                                         * NSYS systems for receiver
                                         * only 0 is used for satellite */
};
/* SBAS fast correction class --------------------------------------------------------------------- */
class gnss_sbsfcorr_c {
/* Constructors */
public:
    gnss_sbsfcorr_c();
    ~gnss_sbsfcorr_c();
/* Components */
public:
    gtime_c t0;                         /* time of applicability (TOF) */
    double prc;                         /* pseudorange correction (PRC) (m) */
    double rrc;                         /* range-rate correction (RRC) (m/s) */
    double dt;                          /* range-rate correction delta-time (s) */
    int iodf;                           /* IODF (issue of date fast corr) */
    short udre;                         /* UDRE+1 */
    short ai;                           /* degradation factor indicator */
};

/* SBAS long term satellite error correction class ------------------------------------------------ */
class gnss_sbslcorr_c {
/* Constructors */
public:
    gnss_sbslcorr_c();
    ~gnss_sbslcorr_c();
/* Components */
public:
    gtime_c t0;                         /* correction time */
    int iode;                           /* IODE (issue of date ephemeris) */
    double dpos[3];                     /* delta position (m) (ecef) */
    double dvel[3];                     /* delta velocity (m/s) (ecef) */
    double daf0,daf1;                   /* delta clock-offset/drift (s,s/s) */
};

/* SBAS satellite correction class ---------------------------------------------------------------- */
class gnss_sbssatp_c {
/* Constructors */
public:
    gnss_sbssatp_c();
    ~gnss_sbssatp_c();
/* Components */
public:
    int sat;                            /* satellite number */
    gnss_sbsfcorr_c fcorr;              /* fast correction */
    gnss_sbslcorr_c lcorr;              /* long term correction */
};

/* SBAS satellite corrections class --------------------------------------------------------------- */
class gnss_sbssat_c {
/* Constructors */
public:
    gnss_sbssat_c();
    ~gnss_sbssat_c();
/* Implementation functions */
protected:
    /* long term correction --------------------------------------------------- */
    int sbslongcorr(gnss_obsd_c *data,double &dclk) const;
    /* fast correction -------------------------------------------------------- */
    int sbsfastcorr(gnss_obsd_c *data,double &prc) const;
public:
    /* sbas satellite ephemeris and clock correction -------------------------- */
    int sbssatcorr(gnss_obsd_c *data) const;
/* Components */
public:
    int iodp;                           /* IODP (issue of date mask) */
    int nsat;                           /* number of satellites */
    int tlat;                           /* system latency (s) */
    gnss_sbssatp_c sat[MAXSAT];         /* satellite correction */
};

/* SBAS ionospheric correction class -------------------------------------------------------------- */
class gnss_sbsigp_c {
/* Constructors */
public:
    gnss_sbsigp_c();
    ~gnss_sbsigp_c();
/* Components */
public:
    gtime_c t0;                         /* correction time */
    short lat,lon;                      /* latitude/longitude (deg) */
    short give;                         /* GIVI+1 */
    float delay;                        /* vertical delay estimate (m) */
};

/* SBAS ionospheric corrections class ------------------------------------------------------------- */
class gnss_sbsion_c {
/* Constructors */
public:
    gnss_sbsion_c();
    ~gnss_sbsion_c();
/* Components */
public:
    int iodi;                           /* IODI (issue of date ionos corr) */
    int nigp;                           /* number of igps */
    gnss_sbsigp_c igp[MAXNIGP];         /* ionospheric correction */
};

/* DGPS/GNSS correction class --------------------------------------------------------------------- */
class gnss_dgps_c {
/* Constructors */
public:
    gnss_dgps_c();
    ~gnss_dgps_c();
/* Components */
public:
    gtime_c t0;                         /* correction time */
    double prc;                         /* pseudorange correction (PRC) (m) */
    double rrc;                         /* range rate correction (RRC) (m/s) */
    int iod;                            /* issue of data (IOD) */
    double udre;                        /* UDRE */
};

/* SSR correction class --------------------------------------------------------------------------- */
class gnss_ssr_c {
/* Constructors */
public:
    gnss_ssr_c();
    ~gnss_ssr_c();
/* Components */
public:
    gtime_c t0[6];                      /* epoch time (GPST) {eph,clk,hrclk,ura,bias,pbias} */
    double udi[6];                      /* SSR update interval (s) */
    int iod[6];                         /* iod ssr {eph,clk,hrclk,ura,bias,pbias} */
    int iode;                           /* issue of data */
    int iodcrc;                         /* issue of data crc for beidou/sbas */
    int ura;                            /* URA indicator */
    int refd;                           /* sat ref datum (0:ITRF,1:regional) */
    double deph[3];                     /* delta orbit {radial,along,cross} (m) */
    double ddeph[3];                    /* dot delta orbit {radial,along,cross} (m/s) */
    double dclk[3];                     /* delta clock {c0,c1,c2} (m,m/s,m/s^2) */
    double hrclk;                       /* high-rate clock corection (m) */
    float  cbias[MAXCODE];              /* code biases (m) */
    double pbias[MAXCODE];              /* phase biases (m) */
    float  stdpb[MAXCODE];              /* std-dev of phase biases (m) */
    double yaw_ang,yaw_rate;            /* yaw angle and yaw rate (deg,deg/s) */
    unsigned char update;               /* update flag (0:no update,1:update) */
};

/* QZSS LEX ephemeris class ----------------------------------------------------------------------- */
class gnss_lexeph_c {
/* Constructors */
public:
    gnss_lexeph_c();
    ~gnss_lexeph_c();
/* Components */
public:
    gtime_c toe;                        /* epoch time (GPST) */
    gtime_c tof;                        /* message frame time (GPST) */
    int sat;                            /* satellite number */
    unsigned char health;               /* signal health (L1,L2,L1C,L5,LEX) */
    unsigned char ura;                  /* URA index */
    double pos[3];                      /* satellite position (m) */
    double vel[3];                      /* satellite velocity (m/s) */
    double acc[3];                      /* satellite acceleration (m/s2) */
    double jerk[3];                     /* satellite jerk (m/s3) */
    double af0,af1;                     /* satellite clock bias and drift (s,s/s) */
    double tgd;                         /* TGD */
    double isc[8];                      /* ISC */
};

/* QZSS LEX ionosphere correction class ----------------------------------------------------------- */
class gnss_lexion_c {
/* Constructors */
public:
    gnss_lexion_c();
    ~gnss_lexion_c();
/* Components */
public:
    gtime_c t0;                         /* epoch time (GPST) */
    double tspan;                       /* valid time span (s) */
    double pos0[2];                     /* reference position {lat,lon} (rad) */
    double coef[3][2];                  /* coefficients lat x lon (3 x 2) */
};

/* stec data class -------------------------------------------------------------------------------- */
class gnss_stec_c {
/* Constructors */
public:
    gnss_stec_c();
    ~gnss_stec_c();
/* Components */
public:
    gtime_c time;                       /* time (GPST) */
    unsigned int sat;                   /* satellite number */
    double ion;                         /* slant ionos delay (m) */
    float std;                          /* std-dev (m) */
    float azel[2];                      /* azimuth/elevation (rad) */
    unsigned char flag;                 /* fix flag */
};

/* trop data class -------------------------------------------------------------------------------- */
class gnss_trop_c {
/* Constructors */
public:
    gnss_trop_c();
    ~gnss_trop_c();
/* Components */
public:
    gtime_c time;                       /* time (GPST) */
    double trp[3];                      /* zenith tropos delay/gradient (m) */
    float std[3];                       /* std-dev (m) */
};

/* ppp corrections class -------------------------------------------------------------------------- */
class gnss_pppcorr_c {
/* Constructors */
public:
    gnss_pppcorr_c();
    ~gnss_pppcorr_c();
/* Components */
public:
    int nsta;                           /* number of stations */
    string stas[MAXSTA];                /* station names */
    double rr[MAXSTA][3];               /* station ecef positions (m) */
    int ns[MAXSTA],nsmax[MAXSTA];       /* number of stec data */
    int nt[MAXSTA],ntmax[MAXSTA];       /* number of trop data */
    gnss_stec_c *stec[MAXSTA];          /* stec data */
    gnss_trop_c *trop[MAXSTA];          /* trop data */
};

/* class of navigation data ----------------------------------------------------------------------- */
class gnss_nav_c {
/* Constructors */
public:
    gnss_nav_c();
    ~gnss_nav_c();
/* Implementation functions */
public:
    /* satellite carrier wave length ---------------------------------------------- */
    double satwavelen(int SatNum,int FrqNum);
    /* update satellite carrier wave lengths -------------------------------------- */
    void update_sat_lambda();
    void update_sat_lambda(vector<gnss_ssat_c> &ssat);

/* Components */
public:
    int n,nmax;                         /* number of broadcast ephemeris */
    int ng,ngmax;                       /* number of glonass ephemeris */
    int ns,nsmax;                       /* number of sbas ephemeris */
    int ne,nemax;                       /* number of precise ephemeris */
    int nc,ncmax;                       /* number of precise clock */
    int na,namax;                       /* number of almanac data */
    int nt,ntmax;                       /* number of tec grid data */
    int nf,nfmax;                       /* number of satellite fcb data */
    vector<gnss_eph_c> eph;             /* GPS/QZS/GAL ephemeris */
    vector<gnss_geph_c> geph;           /* GLONASS ephemeris */
    vector<gnss_seph_c> seph;           /* SBAS ephemeris */
    vector<gnss_peph_c> peph;           /* precise ephemeris */
    vector<gnss_pclk_c> pclk;           /* precise clock */
    vector<gnss_alm_c> alm;             /* almanac data */
    vector<gnss_tec_c> tec;             /* tec grid data */
    vector<gnss_ztd_c> ztd;             /* ZTD data */
    vector<gnss_fcbd_c> fcb;            /* satellite fcb data */
    gnss_erp_c erp;                     /* earth rotation parameters */
    double ocean_par[2][6*11];          /* ocean tide loading parameters {rov,base} */
    double utc_gps[4];                  /* GPS delta-UTC parameters {A0,A1,T,W} */
    double utc_glo[4];                  /* GLONASS UTC GPS time parameters */
    double utc_gal[4];                  /* Galileo UTC GPS time parameters */
    double utc_qzs[4];                  /* QZS UTC GPS time parameters */
    double utc_cmp[4];                  /* BeiDou UTC parameters */
    double utc_irn[4];                  /* IRNSS UTC parameters */
    double utc_sbs[4];                  /* SBAS UTC parameters */
    double ion_gps[8];                  /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_gal[4];                  /* Galileo iono model parameters {ai0,ai1,ai2,0} */
    double ion_qzs[8];                  /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_cmp[8];                  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_irn[8];                  /* IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    int leaps;                          /* leap seconds (s) */
    double lam[MAXSAT][NFREQ];          /* carrier wave lengths (m) */
    char glo_fcn[MAXPRNGLO+1];          /* glonass frequency channel number + 8 */
    gnss_freq_c freq;                   /* gnss frequency priority and adoption */
    gnss_bias_c obsbias;                /* observation bias */
    gnss_pcv_c pcvs[MAXSAT];            /* satellite antenna pcv */
    gnss_sbssat_c sbssat;               /* SBAS satellite corrections */
    gnss_sbsion_c sbsion[MAXBAND+1];    /* SBAS ionosphere corrections */
    gnss_dgps_c dgps[MAXSAT];           /* DGPS corrections */
    gnss_ssr_c ssr[MAXSAT];             /* SSR corrections */
    gnss_lexeph_c lexeph[MAXSAT];       /* LEX ephemeris */
    gnss_lexion_c lexion;               /* LEX ionosphere correction */
    gnss_pppcorr_c pppcorr;             /* ppp corrections */
};

/* SBAS message class ----------------------------------------------------------------------------- */
class gnss_sbsmsg_c {
/* Constructors */
public:
    gnss_sbsmsg_c();
    ~gnss_sbsmsg_c();
/* Implementation functions */
protected:
    /* decode half long term correction (vel code=0) --------------------------- */
    int decode_longcorr0(int p,gnss_sbssat_c *sbssat);
    /* decode half long term correction (vel code=1) --------------------------- */
    int decode_longcorr1(int p,gnss_sbssat_c *sbssat);
    /* decode half long term correction ---------------------------------------- */
    int decode_longcorrh(int p,gnss_sbssat_c *sbssat);
    /* decode type 1: prn masks ------------------------------------------------ */
    int decode_sbstype1(gnss_sbssat_c *sbssat);
    /* decode type 2-5,0: fast corrections ------------------------------------- */
    int decode_sbstype2(gnss_sbssat_c *sbssat);
    /* decode type 6: integrity info ------------------------------------------- */
    int decode_sbstype6(gnss_sbssat_c *sbssat);
    /* decode type 7: fast correction degradation factor ----------------------- */
    int decode_sbstype7(gnss_sbssat_c *sbssat);
    /* decode type 9: geo navigation message ----------------------------------- */
    int decode_sbstype9(gnss_nav_c *nav);
    /* decode type 18: ionospheric grid point masks ---------------------------- */
    int decode_sbstype18(gnss_sbsion_c *sbsion);
    /* decode type 24: mixed fast/long term correction ------------------------- */
    int decode_sbstype24(gnss_sbssat_c *sbssat);
    /* decode type 25: long term satellite error correction -------------------- */
    int decode_sbstype25(gnss_sbssat_c *sbssat);
    /* decode type 26: ionospheric deley corrections --------------------------- */
    int decode_sbstype26(gnss_sbsion_c *sbsion);
public:
    /* update sbas corrections to nav ------------------------------------------ */
    int sbsupdatecorr(gnss_nav_c *nav);
/* Components */
public:
    int week,tow;                       /* receiption time */
    int prn;                            /* SBAS satellite PRN number */
    unsigned char msg[29];              /* SBAS message (226bit) padded by 0 */
};
/* SBAS messages class ---------------------------------------------------------------------------- */
class gnss_ssbs_c {
/* Constructors */
public:
    gnss_ssbs_c();
    ~gnss_ssbs_c();
/* Components */
public:
    int n,nmax;                         /* number of SBAS messages/allocated */
    vector<gnss_sbsmsg_c> msgs;         /* SBAS messages */
};

/* QZSS LEX message class ------------------------------------------------------------------------- */
class gnss_lexmsg_c {
/* Constructors */
public:
    gnss_lexmsg_c();
    ~gnss_lexmsg_c();
/* Implementation functions */
protected:
    /* decode tof and toe field (ref [1] 5.7.2.2.1.1) ----------------------------*/
    int decode_lextof(int i,gtime_c &tof,gtime_c &toe);
    /* decode signal health field (ref [1] 5.7.2.2.1.1) --------------------------*/
    int decode_lexhealth(int i,gtime_c tof,gnss_nav_c *nav);
    /* decode ephemeris and sv clock field (ref [1] 5.7.2.2.1.2) -----------------*/
    int decode_lexeph(int i,gtime_c toe,gnss_nav_c *nav);
    /* decode ionosphere correction field (ref [1] 5.7.2.2.1.3) ------------------*/
    int decode_lexion(int i,gtime_c tof,gnss_nav_c *nav);
    /* convert lex type 12 to rtcm ssr message -----------------------------------*/
    int lex2rtcm(int i,unsigned char *buff);
    /* decode type 10: ephemeris data and clock (ref [1] 5.7.2.2.1,1) ------- */
    int decode_lextype10(gnss_nav_c *nav,gtime_c tof);
    /* decode type 11: ephemeris data and clock (ref [1] 5.7.2.2.1,1) ------- */
    int decode_lextype11(gnss_nav_c *nav,gtime_c tof);
    /* decode type 12: madoca orbit and clock correction -------------------- */
    int decode_lextype12(gnss_nav_c *nav,gtime_c tof);
    /* decode type 20: gsi experiment message (ref [1] 5.7.2.2.2) ----------- */
    int decode_lextype20(gnss_nav_c *nav,gtime_c tof);
public:
    /* update lex corrections ----------------------------------------------- */
    int lexupdatecorr(gnss_nav_c *nav,gtime_c &tof);
/* Components */
public:
    int prn;                            /* satellite PRN number */
    int type;                           /* message type */
    int alert;                          /* alert flag */
    unsigned char stat;                 /* signal tracking status */
    unsigned char snr;                  /* signal C/N0 (0.25 dBHz) */
    unsigned int ttt;                   /* tracking time (ms) */
    unsigned char msg[212];             /* LEX message data part 1695 bits */
};
/* QZSS LEX messages class ------------------------------------------------------------------------ */
class gnss_lex_c {
/* Constructors */
public:
    gnss_lex_c();
    ~gnss_lex_c();
/* Components */
public:
    int n,nmax;                         /* number of LEX messages and allocated */
    vector<gnss_lexmsg_c> msgs;         /* LEX messages */
};
#endif
