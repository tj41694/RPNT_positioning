#include "GNSS/ParModel/parameter.h"

/* parameter function ----------------------------------------------------------------------------- */
/* Constructor -------------------------------------------------------------------- */
gnss_parameter_c::gnss_parameter_c() {
}
gnss_parameter_c::~gnss_parameter_c() {
}
/* Implementaion functions -------------------------------------------------------- */
/* frequency number --------------------------------------------------------------- */
int gnss_parameter_c::N_Freqency(gnss_prcopt_c *opt) {
    return opt->posF;
}
/* position parameter number (xyz or xyz) ----------------------------------------- */
int gnss_parameter_c::N_Position(gnss_prcopt_c *opt) {
    return 3;
}
/* velocity parameter number (xyz or xyz) ----------------------------------------- */
int gnss_parameter_c::N_Velocity(gnss_prcopt_c *opt) {
    return 3;
}
/* troposphere parameter number --------------------------------------------------- */
int gnss_parameter_c::N_Tro(gnss_prcopt_c *opt) {
    if ( opt->mode==PMODE_SINGLE || opt->tropopt<TROPOPT_EST || opt->tropopt>TROPOPT_ESTG ) return 0;
    return (opt->mode<PMODE_DGPS ? 1 : 2) // PPP or RTK
        * (opt->tropopt==TROPOPT_EST ? 1 : 3);
}
/* receiver GLONASS IFB parameters ------------------------------------------------ */
int gnss_parameter_c::N_GLOIFB(gnss_prcopt_c *opt) {
    return opt->navsys&SYS_GLO&&opt->glomodear>0&&opt->ionoopt!=IONOOPT_IFLC&&opt->mode>PMODE_SINGLE ?
        (NFREQGLO>opt->posF? opt->posF : NFREQGLO) : 0;
}
/* ionospherer delay (m) in GPS L1 ------------------------------------------------ */
int gnss_parameter_c::N_Ion(gnss_prcopt_c *opt) {
    return ( opt->mode>PMODE_SINGLE && opt->ionoopt==IONOOPT_CONST )? MAXSAT : 0;
}
/* only for ppp position ---------------------------------------------------------- */
int gnss_parameter_c::N_Clock(gnss_prcopt_c *opt) {
    return opt->mode<PMODE_DGPS ? NSYS : 0;
}
/* receiver uncalibrated phase delay (only for ppp) ------------------------------- */
int gnss_parameter_c::N_UPD(gnss_prcopt_c *opt) {
    return 0;
}
/* ambiguity parameters ----------------------------------------------------------- */
int gnss_parameter_c::N_Amb(gnss_prcopt_c *opt,const int numF) {
    return opt->mode>PMODE_SINGLE&&opt->mode!=PMODE_DGPS? numF*MAXSAT : 0;
}