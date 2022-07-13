#ifndef SATELLTIE_H
#define SATELLTIE_H

#include "GNSS/DataClass/data.h"

#include "GNSS/AntModel/antenna.h"

/* satellite functions ---------------------------------------------------------------------------- */
/* parent class of satellite functions ------------------------------------------------------------ */
class gnss_satellite_c {
/* Constructors */
public:
    gnss_satellite_c();
    virtual ~gnss_satellite_c();
/* Implementation functions */
protected:
    /* virtual satellite clocks function -------------------------------------- */
    virtual int satclk(gnss_obsd_c *data,const gnss_nav_c *nav);
public:
    /* virtual satellite position function ------------------------------------ */
    virtual int satpos(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
    /* compute satellite positions and clocks --------------------------------- */
    virtual int satposclk(gnss_obs_c *obs,const gnss_nav_c *nav);
    /* update satellite positions and clocks with the receiver clock bias ----- */
    void upposclk(gnss_obs_c *obs,const gnss_nav_c *nav,const vector<double> &pos);
/* Components */
public:
    gnss_satant_c *satantfunc;          /* satellite antenna functions (point to gnss_pro_c) */
};

/* subclass of satellite functions ---------------------------------------------------------------- */
/* broadcast ephemeris ---------------------------------------------------------------------------- */
class gnss_broadcast_c : public gnss_satellite_c {
/* Constructors */
public:
    gnss_broadcast_c();
    virtual ~gnss_broadcast_c();
/* Implementation functions */
protected:
    /* broadcast functions ---------------------------------------------------- */
    /* variance by ura ephemeris (ref [1] 20.3.3.3.1.1) ----------------------- */
    double var_uraeph(int ura);
    /* variance by ura ssr (ref [4]) ------------------------------------------ */
    double var_urassr(int ura);
    /* glonass orbit differential equations ----------------------------------- */
    void deq(const double *x,double *xdot,const double *acc);
    /* glonass position and velocity by numerical integration ----------------- */
    void glorbit(double t,double *x,const double *acc);

    /* select GPS/GAL/QZS/CMP ephemeris --------------------------------------- */
    int seleph(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
    /* select GLO ephemeris --------------------------------------------------- */
    int selgeph(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
    /* select SBS ephemeris --------------------------------------------------- */
    int selseph(gnss_obsd_c *data,const gnss_nav_c *nav);
    /* position functions ----------------------------------------------------- */
    /* position from GPS/GAL/QZS/CMP ephemeris -------------------------------- */
    void eph2pos(gnss_obsd_c *data,const gnss_eph_c eph);
    /* position bias from GLO ephemeris --------------------------------------- */
    void geph2pos(gnss_obsd_c *data,const gnss_geph_c geph);
    /* position bias from SBS ephemeris --------------------------------------- */
    void seph2pos(gnss_obsd_c *data,const gnss_seph_c seph);
    /* clock functions -------------------------------------------------------- */
    /* clock bias from GPS/GAL/QZS/CMP ephemeris ------------------------------ */
    void eph2clk(gnss_obsd_c *data,const gnss_eph_c eph);
    /* clock bias from GLO ephemeris ------------------------------------------ */
    void geph2clk(gnss_obsd_c *data,const gnss_geph_c geph);
    /* clock bias from SBS ephemeris ------------------------------------------ */
    void seph2clk(gnss_obsd_c *data,const gnss_seph_c seph);
    /* test availability of navigation file ----------------------------------- */
    int test_nav_availability(gnss_obsd_c *data,const gnss_nav_c *nav,int neph);

    /* broadcast satellite position function ---------------------------------- */
    int broadpos(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
    /* broadcast satellite clocks function ------------------------------------ */
    int satclk(gnss_obsd_c *data,const gnss_nav_c *nav);
public:
    /* broadcast satellite position function called by satposclk() ------------ */
    virtual int satpos(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
};

/* broadcast ephemeris with sbas correction ------------------------------------------------------- */
class gnss_broadsbas_c : public gnss_broadcast_c {
/* Constructors */
public:
    gnss_broadsbas_c();
    ~gnss_broadsbas_c();
/* Implementation functions */
protected:
    /* base functions --------------------------------------------------------------- */
public:
    /* broadcast satellite position function with sbas correction ------------------- */
    int satpos(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
};
/* broadcast ephemeris with ssr_apc correction ---------------------------------------------------- */
class gnss_broadssrapc_c : public gnss_broadcast_c {
/* Constructors */
public:
    gnss_broadssrapc_c();
    ~gnss_broadssrapc_c();
/* Implementation functions */
protected:
    /* base functions --------------------------------------------------------------- */
public:
    /* broadcast satellite position function with ssr_apc correction ---------------- */
    int satpos(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
};
/* broadcast ephemeris with ssr_com correction ---------------------------------------------------- */
class gnss_broadssrcom_c : public gnss_broadcast_c {
/* Constructors */
public:
    gnss_broadssrcom_c();
    ~gnss_broadssrcom_c();
/* Implementation functions */
protected:
    /* base functions --------------------------------------------------------------- */
public:
    /* broadcast satellite position function with ssr_com correction ---------------- */
    int satpos(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
};

/* precise ephemeris ------------------------------------------------------------------------------ */
class gnss_preciseph_c : public gnss_broadcast_c {
/* Constructors */
public:
    gnss_preciseph_c();
    ~gnss_preciseph_c();
/* Implementation functions */
protected:
    /* base precise ephemeris functions --------------------------------------- */
    /* polynomial interpolation by Lagrange's algorithm (time interpolation) -- */
    double interpolLaR(const double *dt,gtime_c *ptime,const double *ppos,
        int n);
    /* precise satellite position of one epoch -------------------------------- */
    int precisepos(gnss_obsd_c *data,const gnss_nav_c *nav);
    /* precise satellite clocks function -------------------------------------- */
    int satclk(gnss_obsd_c *data,const gnss_nav_c *nav);
public:
    /* precise satellite position function ------------------------------------ */
    int satpos(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
    /* compute satellite positions and clocks for test ------------------------ */
    virtual int satposclk(gnss_obs_c *obs,const gnss_nav_c *nav);

/* Components */
protected:
    double scvar;                       /* satellite clock variance */
};

/* qzss lex ephemeris ----------------------------------------------------------------------------- */
class gnss_qzsslex_c : public gnss_broadcast_c {
/* Constructors */
public:
    gnss_qzsslex_c();
    ~gnss_qzsslex_c();
/* Implementation functions */
protected:
    /* base functions --------------------------------------------------------- */
    /* ura value -------------------------------------------------------------- */
    double vareph(int ura);
    /* qzss lex satellite clocks function ------------------------------------- */
    int satclk(gnss_obsd_c *data,const gnss_nav_c *nav);
public:
    /* qzss lex satellite position function ----------------------------------- */
    int satpos(gnss_obsd_c *data,int iode,const gnss_nav_c *nav);
};
#endif