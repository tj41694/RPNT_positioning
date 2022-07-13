#ifndef TIDELOAD_H
#define TIDELOAD_H

#include "ConfigFile/config.h"
#include "GNSS/DataClass/data.h"

/* earth, ocean loading and pole tide function ---------------------------------------------------- */
class gnss_tidecorr_c {
/* Constructor */
public:
    gnss_tidecorr_c();
    ~gnss_tidecorr_c();
    /* Implementaion functions */
protected:
    /*iers mean pole (ref [7] eq.7.25) ---------------------------------------- */
    void iers_mean_pole(double *xp_bar,double *yp_bar);
    /* solar/lunar tides (ref [2] 7) ------------------------------------------ */
    void tide_sl(const double *eu,const double *rp,double GMp,
        const double *pos,double *dr);
    /* compute pole tide (ref [7] eq.7.26) ------------------------------------ */
    void tide_pole();
    /* compute ocean-loading tide (ref [2] 7) --------------------------------- */
    void tide_ocean(int rovbas);
    /* compute solid earth tide (ref [2] 7) ----------------------------------- */
    void tide_solid();
public:
    /* initialize ocean loading parameter pointer ----------------------------- */
    void init_otl(gnss_prcopt_c *opt,gnss_nav_c *nav);
    /* initialize earth rotation parameter pointor ---------------------------- */
    void init_erp(gnss_nav_c *nav);
    /* compute tidal displacement corretion ----------------------------------- */
    void tidecorr(gtime_c time,const int rovbas,const double *xyz,double *dr);
/* Components */
protected:
    gtime_c tut;
    double radius;                      /* station ecef radius */
    double blhpos[2];                   /* station geodetic position {lat,lon(,rad)} */
    double C_e_n[9];                    /* ecef to local coordinates tranformation matrix (3x3) */
    double mode_dr[3][3];               /* dr of different tidal correction {tide,ocean,pole} */
    double denu[2][3];                  /* ocean and pole tide in enu */
    double sun_ecef[3];                 /* sun position in ecef */
    double moon_ecef[3];                /* moon position in ecef */
    double gmst;                        /* greenwich mean sideral time */
    double erpv[5];                     /* earth rotation values{xp,yp,ut1_utc,lod} (rad,rad,s,s/d) */
public:
    int tide_opt;                       /* tide option
                                         * 1: solid earth tide
                                         * 2: ocean tide loading
                                         * 4: pole tide
                                         * 8: elimate permanent deformation */
    double *ocean_par[2];               /* ocean-loading parameters pointor to opt {rov,base} */
    gnss_erp_c *er_par;                 /* earth rotation parameters pointor to opt */
};

#endif
