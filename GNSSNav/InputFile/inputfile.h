/* Inpute Files */

#ifndef INPUTFILE_H
#define INPUTFILE_H

#include "BaseFunction/timesys.h"
#include "GNSS/DataClass/data.h"
#include "ConfigFile/config.h"

/* GNSS input file -------------------------------------------------------------------------------- */
/* parent class in_gnss_rnx_c ------------------------------------------------------------------------
* O,N,G,H,J,L,C file ------------------------------------------------------------------------------ */
class in_gnss_rnx_c {
/* Constructors */
public:
    in_gnss_rnx_c();
    in_gnss_rnx_c(string file,string option);
    virtual ~in_gnss_rnx_c();
/* Implementation functions */
protected:
    /* set system mask ------------------------------------------------------------ */
    void set_sysmask(const gnss_prcopt_c *prcopt);
public:
    /* test oepn of file stream --------------------------------------------------- */
    int test_open();
    /* close file ----------------------------------------------------------------- */
    void closeF();
    /* base function of read head ------------------------------------------------- */
    int readhead();
    /* vritual function of read head for O,N,G,H,J,L,C ---------------------------- */
    virtual int Head();
    /* virtual function of read body for O,N,G,H,J,L,C ---------------------------- */
    virtual int Body();
/* Components */
protected:
    ifstream inf;                   /* in-stream of RINEX file */
    string opt;                     /* options */
    double ver;                     /* rinex version */
    string type;                    /* type of rinex file */
    int sat_sys;                    /* satellite system */
    int tsys;                       /* time system */
    gtime_c ts;                     /* processing start time */
    gtime_c te;                     /* processing end time */
    double tint;                    /* processing interval */
    int mask;                       /* mask satellite system */
    string buff;                    /* read buff */
public:
    string errorMsg;                /* read rinex file error message */
};

/* derived class in_gnss_rnxO_c ----------------------------------------------------------------------
* O file - Observation ---------------------------------------------------------------------------- */
class in_gnss_rnxO_c : public in_gnss_rnx_c {
/* Constructors */
public:
    in_gnss_rnxO_c();
    in_gnss_rnxO_c(string file,string option,int rcvnum);
    in_gnss_rnxO_c(string file,string option,gtime_c TS,gtime_c TE,int TI,int rcvnum);
    ~in_gnss_rnxO_c();

/* Implementation functions */
protected:
    /* convert rinex obs type ver.2 -> ver.3 -------------------------------------- */
    void convcode(int sys,string str,string &tp);
    /* save slips ----------------------------------------------------------------- */
    void saveslips(unsigned char slips[][NFREQ],gnss_obsd_c &data);
    /* restore slips -------------------------------------------------------------- */
    void restslips(unsigned char slips[][NFREQ],gnss_obsd_c &data);
    /* set signal index ----------------------------------------------------------- */
    void set_index(int syt,string tobss[MAXOBSTYPE],sigind_t &ind,
                   const gnss_freq_c *freq);
    /* decode obs epoch ----------------------------------------------------------- */
    int decode_obsepoch(gtime_c &t,int &flag,vector<int> &sats);
    /* decode obs data ------------------------------------------------------------ */
    int decode_obsdata(gnss_obsd_c &obs);
    /* read "O" file to gnss_obs_c vector ----------------------------------------- */
    int readrnxobsb(int &flag,vector<gnss_obsd_c> &data,const gnss_pstopt_c *pstopt);
    /* inherited functions from class in_gnss_rnx_c */
public:
    /* set used frequency index according to "GNSS_USED_FREQ[MAX_GNSS][NFREQ*2]" -- */
    void set_freq_index(gnss_nav_c *nav);
    /* initialization ------------------------------------------------------------- */
    void ini_ReadO(string file,string option,gtime_c TS,gtime_c TE,int TI,int rcvnum);
    /* read head of "O" file ------------------------------------------------------ */
    virtual int Head(gnss_nav_c *nav,gnss_sta_c *sta);
    /* read boy of "O" file (don't use this function!!!) -------------------------- */
    virtual int Body(gnss_obs_c *obs,const gnss_prcopt_c *prcopt,const gnss_pstopt_c *pstopt);
    /* read content number gnss_obs_c to the buff --------------------------------- */
    int Read2_Obs_Buff(vector <gnss_obs_c> &obs,const gnss_prcopt_c *prcopt,const gnss_pstopt_c *pstopt);
    /* read one epoch body of "O" file -------------------------------------------- */
    int Read_One_Epoch(gnss_obs_c *obs,const gnss_prcopt_c *prcopt,const gnss_pstopt_c *pstopt);

/* Components */
public:
    int rcv;						        /* receiver number */
    int skipflag;
    gtime_c epochFirst,epochEnd;	        /* first and end epoch of observation file */
    string tobs[MAX_GNSS][MAXOBSTYPE]; 	    /* observations' type of all systems */
    sigind_t index[MAX_GNSS]={ {0} };		/* observation's signal index of all systems */
    int lastUsedCode[MAXSAT][NFREQ][MAX_CODE_TYPE];
                                            /* last used observation code for all satellite */

    double val[MAXOBSTYPE];
    unsigned char lli[MAXOBSTYPE];
    int pos[MAXOBSTYPE],p[MAXOBSTYPE];
};

/* derived class in_gnss_rnxN_c -----------------------------------------------------------------------
* N file - Navigation ephemeris -------------------------------------------------------------------- */
class in_gnss_rnxN_c : public in_gnss_rnx_c {
/* Constructors */
public:
    in_gnss_rnxN_c();
    in_gnss_rnxN_c(string file,string option);
    ~in_gnss_rnxN_c();
/* Implementation functions */
protected:
    /* ura value (m) to ura index ------------------------------------------------- */
    int uraindex(double value);

    /* decode glonass ephemeris --------------------------------------------------- */
    int decode_geph(gtime_c toc,int sat,gnss_geph_c *geph,int ndata);
    /* decode geo ephemeris ------------------------------------------------------- */
    int decode_seph(gtime_c toc,int sat,gnss_seph_c *seph);
    /* decode ephemeris ----------------------------------------------------------- */
    int decode_eph(gtime_c toc,int sat,gnss_eph_c *eph);
    /* add data to gnss_nav_c ----------------------------------------------------- */
    void addnav(int sys,gnss_nav_c *nav);
    /* delete data from gnss_nav_c ------------------------------------------------ */
    void delnav(int sys,gnss_nav_c *nav);
    /* test if the new ephemeris data is good to use ------------------------------ */
    int test_good_eph(int sys,gnss_nav_c *nav);
    /* read "O" file to gnss_obs_c vector ----------------------------------------- */
    int readrnxnavb(gnss_nav_c *nav);
public:
    /* initialization ------------------------------------------------------------- */
    void ini_ReadN(string file,string option);
    /* read head of "O" file ------------------------------------------------------ */
    virtual int Head(gnss_nav_c *nav);
    /* read boy of "O" file ------------------------------------------------------- */
    virtual int Body(gnss_nav_c *nav,const gnss_prcopt_c *prcopt);
/* Components */
protected:
    gtime_c toc,last_toc;
    int nvalue,ngood,last_ngood;
    int sat,last_sat,flag;
public:
    double data[64];				/* one satellite data */
};

/* derived class in_gnss_rnxC_c ----------------------------------------------------------------------
* C file - Precise Clock file --------------------------------------------------------------------- */
class in_gnss_rnxC_c : public in_gnss_rnx_c {
/* Constructors */
public:
    in_gnss_rnxC_c();
    in_gnss_rnxC_c(string file,string option);
    ~in_gnss_rnxC_c();
/* Implementation functions */
protected:
    /* read head of "C" file ------------------------------------------------------ */
    virtual int Head(gnss_nav_c *nav);
    /* read boy of "C" file ------------------------------------------------------- */
    virtual int Body(gnss_nav_c *nav,const gnss_prcopt_c *prcopt);
public:
    /* initialization ------------------------------------------------------------- */
    void ini_ReadC(string file,string option);
    /* read precise clocks file --------------------------------------------------- */
    int readclk(gnss_nav_c *nav,const gnss_prcopt_c *opt);
/* Components */
protected:
    int nameLen;                    /* length of Receiver or satellite name */
    int iTime,iClk,iVar;            /* first char place of time, clock and variance */
    int lClk,lVar;                  /* length of one line to clock or variance */
};

/* read earth rotation parameters file (.ERP) ----------------------------------------------------- */
class in_gnss_erp_c {
/* Constructors */
public:
    in_gnss_erp_c();
    in_gnss_erp_c(string file);
    ~in_gnss_erp_c();
/* Implementation functions */
public:
    /* open file ------------------------------------------------------------------ */
    void open_file();
    /* close file ----------------------------------------------------------------- */
    void closeF();
    /* read earth rotation parameters file ---------------------------------------- */
    int readerp(gnss_erp_c *erp);
/* Components */
protected:
    string file_path;				/* erp file path */
    ifstream inf;					/* in-stream of erp file */
};

/* read ocean-loading tide file (.BLQ) ------------------------------------------------------------ */
class in_gnss_blq_c {
/* Constructors */
public:
    in_gnss_blq_c();
    in_gnss_blq_c(string file);
    ~in_gnss_blq_c();
/* Implementation functions */
public:
    /* open file ------------------------------------------------------------------ */
    void open_file();
    /* close file ----------------------------------------------------------------- */
    void closeF();
    /* read earth rotation parameters file ---------------------------------------- */
    int readblq(const string staname,double *ocean_par);
/* Components */
protected:
    string file_path;				/* blq file path */
    ifstream inf;					/* in-stream of blq file */
};

/* read antenna information file (.ATX) ----------------------------------------------------------- */
class in_gnss_atx_c {
/* Constructors */
public:
    in_gnss_atx_c();
    in_gnss_atx_c(string file);
    ~in_gnss_atx_c();
/* Implementation functions */
protected:
    /* convert frequency number to hprtk format (NFREQ) --------------------------- */
    int freq2hprtk(const int prn,const string fstr,int &sys,int &ifreq,const gnss_nav_c *nav);
    /* test satellite system ------------------------------------------------------ */
    int test_sys(const int sys);
public:
    /* open file ------------------------------------------------------------------ */
    void open_file();
    /* close file ----------------------------------------------------------------- */
    void closeF();
    /* read satellite antenna phase center correction file ------------------------ */
    int read_satatx(const gnss_nav_c *nav);
    /* read receiver antenna phase center correction file ------------------------- */
    int read_recatx(const gnss_prcopt_c *opt,const gnss_nav_c *nav);
/* Components */
protected:
    string file_path;				/* atx file path */
    ifstream inf;					/* in-stream of atx file */
    string buff;
public:
    vector<gnss_pcv_c> satpcv;		/* satellite pcv vector */
    gnss_pcv_c recpcv[2];			/* receiver pcv (0:rover, 1:base) */
};

/* read precise ephemeris file -------------------------------------------------------------------- */
class in_gnss_eph_c {
/* Constructors */
public:
    in_gnss_eph_c();
    in_gnss_eph_c(string file,int pred);
    ~in_gnss_eph_c();
/* Implementation functions */
protected:
    /* set system mask ------------------------------------------------------------ */
    void set_sysmask(const gnss_prcopt_c *prcopt);
    /* read precise ephemeris file header ----------------------------------------- */
    int Head();
    /* read precise ephemeris file body ------------------------------------------- */
    int Body(gnss_nav_c *nav);
public:
    /* initialization ------------------------------------------------------------- */
    void ini_readEph(string file,int pred);
    /* open file ------------------------------------------------------------------ */
    void open_file();
    /* close file ----------------------------------------------------------------- */
    void closeF();
    /* read precise ephemeris file ------------------------------------------------ */
    int readsp3(gnss_nav_c *nav,const gnss_prcopt_c *opt);
/* Components */
protected:
    char type;						/* precise ephemeris type */
    string satellite[MAXSAT];		/* satellites flag */
    int ns;							/* number of satellites */
    double bfact[2];				/* base fact for pos/vel and clock bias */
    string tsys;					/* time system */
    string buff;					/* stream buff */
    int mask;						/* mask satellite system */
    string opt;						/* options */
    int pred_flag;					/* flag of use predicted data
                                     * 1: only observed + 2: only predicted + 4: not combined */
public:
    gtime_c ephtime;				/* precise ephemeris time */
    string file_path;				/* sp3/eph file path */
    ifstream inf;					/* in-stream of sp3 file */
};

/* read ionex tec grid file ----------------------------------------------------------------------- */
class in_gnss_ionex_c {
/* Constructors */
public:
    in_gnss_ionex_c();
    in_gnss_ionex_c(string file);
    ~in_gnss_ionex_c();
/* Implementation functions */
protected:
    /* data index (i:lat,j:lon,k:hgt) --------------------------------------------- */
    int dataindex(int i,int j,int k,const int *ndata);

    /* read P1P2 DCB -------------------------------------------------------------- */
    void P1P2DCB(gnss_nav_c *nav);
    /* read head of ionex tec grid file ------------------------------------------- */
    int Head(gnss_nav_c *nav);
    /* add one epoch tec map data to gnss_nav_c ----------------------------------- */
    gnss_tec_c* addtec2nav(gnss_nav_c *nav);
    /* read body of ionex tec grid file ------------------------------------------- */
    int Body(gnss_nav_c *nav);
public:
    /* initialization ------------------------------------------------------------- */
    void ini_readIonex(string file);
    /* open file ------------------------------------------------------------------ */
    void open_file();
    /* close file ----------------------------------------------------------------- */
    void closeF();
    /* read ionex tec grid file --------------------------------------------------- */
    int readIonex(gnss_nav_c *nav);
/* Components */
protected:
    double tec_factor;              /* vtec factor */
    string buff;                    /* stream buff */
    gtime_c iontime;
public:
    double lats[3],lons[3],hgts[3]; /* start/end/interval of lat lon hight */
    double REarth;                  /* radius of Earth */
    double version;                 /* version of ionex file */

    string file_path;               /* ionex tec grid file */
    ifstream inf;                   /* in-stream of ionex file */
};

/* read GNSS code bias file ----------------------------------------------------------------------- */
class in_gnss_bias_c {
/* Constructors */
public:
    in_gnss_bias_c();
    in_gnss_bias_c(string file);
    ~in_gnss_bias_c();
/* Implementation functions */
protected:
    /* read head of GNSS code bias file ------------------------------------------- */
    int Head(gnss_nav_c *nav);
    /* read body of GNSS code bias file ------------------------------------------- */
    int Body(gnss_nav_c *nav);
public:
    /* open file ------------------------------------------------------------------ */
    void open_file();
    /* close file ----------------------------------------------------------------- */
    void closeF();
    /* read GNSS code bias file --------------------------------------------------- */
    int readBias(gnss_nav_c *nav);
/* Components */
protected:
    string buff;					/* stream buff */
    int tsYear,tsDoy,teYear,teDoy,  /* start and end time (int, year:doy:sssss) */
        tsSssss,teSssss; 
    string strY1,strD1,strS1,strY2, /* start and end time (string, year:doy:sssss) */
        strD2,strS2;

public:
    string agency;                  /* agency that created the file */
    gtime_c ts,te;                  /* start and end time */
    int type;                       /* bias type: 0:OSB, 1:DSB, 2:ISB  */
    int mode;                       /* bias mode: 0:relative, 1:absolute */
    int tsys;                       /* time system */

    string file_path;               /* ionex tec grid file */
    ifstream inf;                   /* in-stream of ionex file */
};

/* read GNSS SNX file ----------------------------------------------------------------------------- */
class in_gnss_snx_c {
    /* Constructors */
public:
    in_gnss_snx_c();
    in_gnss_snx_c(string file);
    ~in_gnss_snx_c();
    /* Implementation functions */
protected:
    /* read head of GNSS code bias file ------------------------------------------- */
    int Head();
    /* read body of GNSS code bias file ------------------------------------------- */
    int Body(gnss_prcopt_c *opt,gnss_sta_c sta[2]);
public:
    /* open file ------------------------------------------------------------------ */
    void open_file();
    /* close file ----------------------------------------------------------------- */
    void closeF();
    /* read GNSS code bias file --------------------------------------------------- */
    int readSNX(gnss_prcopt_c *opt,gnss_sta_c sta[2]);
    /* Components */
protected:
    string buff;					/* stream buff */
    int tsYear,tsDoy,teYear,teDoy,  /* start and end time (int, year:doy:sssss) */
        tsSssss,teSssss; 
    string strY1,strD1,strS1,strY2, /* start and end time (string, year:doy:sssss) */
        strD2,strS2;
    string staXYZ[3];               /* station ECEF position solution */
    string staName[2];              /* station 4 character name */

public:
    string agency;                  /* agency that created the file */
    gtime_c ts,te;                  /* start and end time */
    int tsys;                       /* time system */

    string file_path;               /* ionex tec grid file */
    ifstream inf;                   /* in-stream of ionex file */
};

/* read troposphere ZTD file (debug) -------------------------------------------------------------- */
class in_gnss_ztd_c {
public:
    in_gnss_ztd_c();
    in_gnss_ztd_c(int iztd,string file);
    ~in_gnss_ztd_c();
/* Implementation functions */
protected:

public:
    /* open file ------------------------------------------------------------------ */
    void open_file();
    /* close file ----------------------------------------------------------------- */
    void closeF();
    /* read ionex tec grid file --------------------------------------------------- */
    int readZTD(gnss_nav_c *nav);
/* Components */
protected:
    double ztd_factor;				/* vtec factor */
    string buff;					/* stream buff */
    gtime_c ztdtime;
public:
    double ztd;

    string file_path;				/* ionex tec grid file */
    ifstream inf;					/* in-stream of ionex file */
};

#endif