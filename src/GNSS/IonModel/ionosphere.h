#ifndef IONOSPHERE_H
#define IONOSPHERE_H

#include "GNSS/DataClass/data.h"
#include "ConfigFile/config.h"

/* parent ionosphere functions (vertical delay of GPS L1) ----------------------------------------- */
class gnss_ioncorr_c {
/* Constructors */
public:
    gnss_ioncorr_c();
    virtual ~gnss_ioncorr_c();
/* Implementation functions */
protected:
    /* base functions */
public:
    /* ionospheric pierce point position -------------------------------------- */
    double ionppp(const double pos[3],const gnss_obsd_c *obs,const double re,
        const double hion,double *posp);
    /* base ionosphere mapping function --------------------------------------- */
    double ionmapf(const double *pos,const double *azel);
    /* virtual compute ionospere correction (vertical delay of GPS L1) -------------- */
    virtual int correction(gnss_obsd_c *obs,const gnss_nav_c *nav,const double pos[3]);
/* Components */
};

/* subclass ionosphere functions ------------------------------------------------------------------ */
/* ionosphere free combination -------------------------------------------------------------------- */
class gnss_LCion_c : public gnss_ioncorr_c {
/* Constructors */
public:
    gnss_LCion_c();
    ~gnss_LCion_c();
/* Implementation functions */
public:
    /* return 0 delya --------------------------------------------------------- */
    virtual int correction(gnss_obsd_c *obs,const gnss_nav_c *nav,const double pos[3]);
};
/* broadcast ionosphere correction (vertical delay of GPS L1) ------------------------------------- */
class gnss_broadion_c : public gnss_ioncorr_c {
/* Constructors */
public:
    gnss_broadion_c();
    virtual ~gnss_broadion_c();
/* Implementation functions */
protected:
    /* base functions --------------------------------------------------------- */
    /* ionosphere Klobuchar model correction (vertical delay of GPS L1) ------- */
    int klobion(gnss_obsd_c *obs,const double ionpara[8],const double pos[3]);
public:
    /* compute broadcast ionospere correction (vertical delay of GPS L1) ------ */
    virtual int correction(gnss_obsd_c *obs,const gnss_nav_c *nav,const double pos[3]);
};

/* qzss broadcast ionosphere correction (vertical delay of GPS L1) -------------------------------- */
class gnss_qzssion_c : public gnss_broadion_c {
/* Constructors */
public:
    gnss_qzssion_c();
    ~gnss_qzssion_c();
/* Implementation functions */
public:
    /* compute qzss broadcast ionospere correction (vertical delay of GPS L1) - */
    virtual int correction(gnss_obsd_c *obs,const gnss_nav_c *nav,const double pos[3]);
};

/* sbas ionosphere correction (vertical delay of GPS L1) ------------------------------------------ */
class gnss_sbasion_c : public gnss_ioncorr_c {
/* Constructors */
public:
    gnss_sbasion_c();
    ~gnss_sbasion_c();
/* Implementation functions */
protected:
    /* base functions --------------------------------------------------------- */
    /* search igps ------------------------------------------------------------ */
    void searchigp(gtime_c time,const double pos[2],const gnss_sbsion_c *ion,
        const gnss_sbsigp_c **igp,double &x,double &y);
    /* variance of ionosphere correction (give=GIVEI+1) ----------------------- */
    double varicorr(int udre);
public:
    /* compute sbas ionospere correction (vertical delay of GPS L1) ----------- */
    virtual int correction(gnss_obsd_c *obs,const gnss_nav_c *nav,const double pos[3]);
};

/* ionex ionosphere correction (vertical delay of GPS L1) ----------------------------------------- */
class gnss_ionexion_c : public gnss_broadion_c {
/* Constructors */
public:
    gnss_ionexion_c();
    virtual ~gnss_ionexion_c();
/* Implementation functions */
protected:
    /* base functions --------------------------------------------------------- */
    /* data index (i:lat,j:lon,k:hgt) ----------------------------------------- */
    int dataindex(int i,int j,int k,const int *ndata);
    /* interpolate tec grid data ---------------------------------------------- */
    int interptec(const gnss_tec_c *tec,int k,const double *posp,double &value,
        double &rms);
    /* ionosphere delay by tec grid data -------------------------------------- */
    int iondelay(const gnss_tec_c *tec,const double *pos,gnss_obsd_c *obs,
        double *delay,double *var);
public:
    /* compute ionex ionospere correction (vertical delay of GPS L1) ---------- */
    virtual int correction(gnss_obsd_c *obs,const gnss_nav_c *nav,const double pos[3]);
};

/* lex ionosphere correction (vertical delay of GPS L1) ------------------------------------------- */
class gnss_lexioncor_c : public gnss_ioncorr_c {
/* Constructors */
public:
    gnss_lexioncor_c();
    ~gnss_lexioncor_c();
/* Implementation functions */
protected:
    /* base functions --------------------------------------------------------- */
public:
    /* compute lex ionospere correction (vertical delay of GPS L1) ------------ */
    virtual int correction(gnss_obsd_c *obs,const gnss_nav_c *nav,const double pos[3]);
};

/* constrained ionosphere model correction (L1) --------------------------------------------------- */
class gnss_estion_c : public gnss_ionexion_c {
/* Constructors */
public:
    gnss_estion_c();
    ~gnss_estion_c();
/* Implementation functions */
protected:
public:
    /* compute constrained ionosphere model correction (vertical delay of GPS L1) ---- */
    virtual int correction(gnss_obsd_c *obs,const gnss_nav_c *nav,const double pos[3]);
};

#endif
