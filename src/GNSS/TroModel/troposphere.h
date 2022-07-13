#ifndef TROPOSPHERE_H
#define TROPOSPHERE_H

#include "GNSS/DataClass/data.h"

/* troposphere model class ------------------------------------------------------------------------ */
class gnss_tromod_c {
    /* Construtors */
public:
    gnss_tromod_c();
    ~gnss_tromod_c();
    /* Implementaion */
public:
/* Components */
public:
    double Ptro[3],Atro[3];              /* troposphere model parameters and its coefficience
                                          * follow the sequency : ZTD, Gradient_North, Gradient_East
                                          * (only used in estimated mode) */
};

/* parent troposphere functions ------------------------------------------------------------------- */
class gnss_trocorr_c {
/* Constructors */
public:
    gnss_trocorr_c();
    virtual ~gnss_trocorr_c();
/* Implementation functions */
protected:
    /* base functions */
    /* triangle functions ----------------------------------------------------- */
    double mapf(double el,double a,double b,double c);
    /* troposphere interpc function ------------------------------------------- */
    double interpc(const double coef[],const double lat);
    /* NMF troposphere mapping function --------------------------------------- */
    double nmftropmapf(const gnss_obsd_c *obs,const double pos[3],double *mapfw);
    /* get meterological parameters ------------------------------------------- */
    void getmet(double lat,double *met);
public:
    /* standard troposphere model (saastamoinen) ------------------------------ */
    int saascorr(gnss_obsd_c *obs,const double pos[3],const double azel[2],
        const double humi);
    /* sbas troposphere model (sbas) ------------------------------------------ */
    int sbascorr(gnss_obsd_c *obs,const double pos[3],const double azel[2]);
    /* virtual troposphere correction ----------------------------------------- */
    virtual int correction(gnss_obsd_c *obs,const double pos[3],const double humi);

/* Components */
public:
    gnss_tromod_c model_est;            /* troposphere model component
                                         * (only used in estimate mode) */
    gnss_nav_c *nav;                    /* gnss navigation data */
};

/* subclass troposphere functions ----------------------------------------------------------------- */
/* saastamoinen model troposphere correction ------------------------------------------------------ */
class gnss_saastro_c : public gnss_trocorr_c {
/* Constructors */
public:
    gnss_saastro_c();
    ~gnss_saastro_c();
/* Implementation functions */
protected:
    /* base functions */
public:
    /* saastamoinen model troposphere correction ------------------------------ */
    int correction(gnss_obsd_c *obs,const double pos[3],const double humi);
};

/* sbas model troposphere correction -------------------------------------------------------------- */
class gnss_sbastro_c : public gnss_trocorr_c {
/* Constructors */
public:
    gnss_sbastro_c();
    ~gnss_sbastro_c();
/* Implementation functions */
protected:
    /* base functions */
public:
    /* sbas model troposphere correction -------------------------------------- */
    int correction(gnss_obsd_c *obs,const double pos[3],const double humi);
};

/* estimated model troposphere correction --------------------------------------------------------- */
class gnss_esttro_c : public gnss_trocorr_c {
/* Constructors */
public:
    gnss_esttro_c();
    ~gnss_esttro_c();
/* Implementation functions */
protected:
    /* base functions */
public:
    /* estimated model troposphere correction ------------------------------------- */
    int correction(gnss_obsd_c *obs,const double pos[3],const double humi);
};

/* ZTD data model troposphere correction ---------------------------------------------------------- */
class gnss_ztdtro_c : public gnss_trocorr_c {
    /* Constructors */
public:
    gnss_ztdtro_c();
    ~gnss_ztdtro_c();
    /* Implementation functions */
protected:
    /* search ZTD data ------------------------------------------------------------ */
    int search_ztd(gnss_obsd_c *obs);
public:
    /* estimated model troposphere correction ------------------------------------- */
    int correction(gnss_obsd_c *obs,const double pos[3],const double humi);
};

#endif
