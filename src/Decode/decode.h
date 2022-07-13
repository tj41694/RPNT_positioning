#ifndef DECODE_H
#define DECODE_H

#include "gnssnav_lib.h"
#include "GNSS/DataClass/data.h"
#include "BaseFunction/timesys.h"
#include "GNSS/gnss_pro.h"

/* decode data for kinds of formats --------------------------------------------------------------- */
class decode_data_c {
    /* Consstructor */
public:
    decode_data_c();
    virtual ~decode_data_c();
    /* Virtual Implementation functions */
public:
    virtual int decode(unsigned char data);
/* Components */
public:
    gtime_c time;                       /* message time */
    gnss_obs_c obs;                     /* observation data */
    gnss_nav_c nav;                     /* satellite ephemerides */
    gnss_sta_c sta;                     /* station parameters */
    gnss_lexmsg_c lexmsg;               /* LEX message */
    gnss_sbsmsg_c sbsmsg;               /* SBAS message */
    gnss_ssr_c ssr[MAXSAT];             /* output of ssr corrections */
    string opt;                         /* receiver dependent options */
    int ephsat;                         /* sat number of update ephemeris (0:no satellite) */
    int format;                         /* receiver stream format (only for raw_c) */
    gnss_dgps_c *dgps;                  /* output of dgps corrections (only for rtcm2) */
    gnss_rtksvr_c  *Svr;                /* Pointer to RTK server structure
                                         * (when running in that environment otherwise NULL) */
};

#endif
