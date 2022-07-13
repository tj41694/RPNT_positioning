#ifndef ANTENNA_H
#define ANTENNA_H

#include "BaseFunction/timesys.h"
#include "GNSS/DataClass/data.h"
#include "ConfigFile/config.h"

/* satellite antenna phase center offest ---------------------------------------------------------- */
class gnss_satant_c {
/* Constructor */
public:
    gnss_satant_c();
    ~gnss_satant_c();
/* Implementation functions */
protected:
    /* nominal yaw-angle ------------------------------------------------------ */
    double yaw_nominal(double beta,double mu);
    /* yaw-angle of satellite ------------------------------------------------- */
    int yaw_angle(int sat,const char *type,int opt,double beta,double mu,
        double *yaw);
    /* yaw-angle of satellite ------------------------------------------------- */
    int sat_yaw(gtime_c time,const double *rs,double *exs,double *eys);
    /* calculate satellite antenna phase center variation --------------------- */
    void antmodel_s(const gnss_pcv_c *pcv,double nadir,gnss_obsd_c *data);
public:
    /* initialize earth rotation parameter pointor ---------------------------- */
    void init_erp(gnss_nav_c *nav);
    /* satellite antenna phase center offest for one gnss_obsd_c ------------------- */
    void satantoff(gnss_obsd_c *data,const gnss_nav_c *nav);
    /* satellite antenna phase center variation correction -------------------- */
    void satantvar(const vector<double> &XYZ,gnss_obsd_c *data,const gnss_nav_c *nav);
    /* satellite antenna phase center offset and variation correction --------- */
    void satantov(const gnss_prcopt_c *opt,const vector<double> &XYZ,gnss_obsd_c *data,gnss_nav_c *nav);
    /* phase windup correction for one gnss_obsd_c --------------------------------- */
    int phase_windup(gtime_c soltime,const vector<double> &XYZ,gnss_obsd_c *data,double &phw);

/* Components */
protected:
    double *lam;                        /* lambda pointer */

public:
    gnss_erp_c *er_par;                 /* earth rotation parameters pointor to opt */
};

/* receiver antenna phase center offest ----------------------------------------------------------- */
class gnss_recant_c {
/* Constructor */
public:
    gnss_recant_c();
    ~gnss_recant_c();
/* Implementation functions */
protected:
public:
    /* receiver antenna phase center offest and variation correction ---------- */
    void recantov(const gnss_prcopt_c *opt,int rovbas,gnss_obsd_c *data,gnss_nav_c *nav);

/* Components */
protected:
    double *lam;                        /* lambda pointer */
};
#endif