#ifndef PARAMETER_H
#define PARAMETER_H

#include "ConfigFile/config.h"

/* parameter function ----------------------------------------------------------------------------- */
class gnss_parameter_c {
/* Constructors */
public:
    gnss_parameter_c();
    ~gnss_parameter_c();
    /* Implementaion functions */
public:
    /* frequency number ------------------------------------------------------- */
    int N_Freqency(gnss_prcopt_c *opt);
    /* position parameter number (xyz) ---------------------------------------- */
    int N_Position(gnss_prcopt_c *opt);
    /* velocity parameter number (xyz) ---------------------------------------- */
    int N_Velocity(gnss_prcopt_c *opt);
    /* troposphere parameter number ------------------------------------------- */
    int N_Tro(gnss_prcopt_c *opt);
    /* receiver GLONASS IFB parameters ---------------------------------------- */
    int N_GLOIFB(gnss_prcopt_c *opt);
    /* ionospherer delay (m) in GPS L1 ---------------------------------------- */
    int N_Ion(gnss_prcopt_c *opt);
    /* receiver clock (only for ppp) ------------------------------------------ */
    int N_Clock(gnss_prcopt_c *opt);
    /* receiver uncalibrated phase delay (only for ppp) ----------------------- */
    int N_UPD(gnss_prcopt_c *opt);
    /* ambiguity parameters --------------------------------------------------- */
    int N_Amb(gnss_prcopt_c *opt,const int numF);
};

#endif
