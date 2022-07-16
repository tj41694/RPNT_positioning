/* Time function head file */

#ifndef TIMESYS_H
#define TIMESYS_H

#include "gnssnav_lib.h"

/* time structure --------------------------------------------------------------------------------- */
class gtime_c {
/* Constructors */
public:
    gtime_c();
    /* initialize with epoch array ------------------------------------------------ */
    gtime_c(const double *epoch);
    gtime_c(int isys);
    ~gtime_c();
public:
/* Implementation function */
    /* string to time ------------------------------------------------------------- */
    int str2time(string s);
    /* ep time to string ---------------------------------------------------------- */
    string time2str(int n);
    /* calender day/time (ep) to time --------------------------------------------- */
    gtime_c *epoch2time(const double *inep);
    /* time to calender day/time (ep) --------------------------------------------- */
    void time2epoch();
    /* gps week time to time ------------------------------------------------------ */
    gtime_c *gpst2time(int week,double sss);
    /* time to gps week tme ------------------------------------------------------- */
    double time2gpst(int *week);
    /* galileo week time to time -------------------------------------------------- */
    gtime_c *gst2time(int week,double sss);
    /* time to galileo week time -------------------------------------------------- */
    double time2gst(int *week);
    /* Beidou week time to time --------------------------------------------------- */
    gtime_c *bdt2time(int week,double sss);
    /* time to Beidou week time --------------------------------------------------- */
    double time2bdt(int *week);
    /* add time ------------------------------------------------------------------- */
    gtime_c *timeadd(double dsec);
    /* difference with other time ------------------------------------------------- */
    double timediff(const gtime_c &t2) const;
    /* get current time in utc ---------------------------------------------------- */
    gtime_c *timeget();
    /* set current time in utc ---------------------------------------------------- */
    void timeset();
    /* gps time to utc ------------------------------------------------------------ */
    gtime_c *gpst2utc();
    /* utc to gps time p----------------------------------------------------------- */
    gtime_c *utc2gpst();
    /* gps time to Beidou time ---------------------------------------------------- */
    gtime_c *gpst2bdt();
    /* Beidou time to gps time ---------------------------------------------------- */
    gtime_c *bdt2gpst();
    /* time to day and sec -------------------------------------------------------- */
    double time2sec(gtime_c &day);
    /* utc to Greenwich mean sidereal time ---------------------------------------- */
    double utc2gmst(double ut1_utc);
    /* day of year to time -------------------------------------------------------- */
    int doy2time(int year,int doy);
    /* time to day of year -------------------------------------------------------- */
    double time2doy();
    /* get time of day ------------------------------------------------------------ */
    double time_of_day();
    /* read leap seconds table ---------------------------------------------------- */
    int read_leaps(const string file);
    /* adjust time considering week handover -------------------------------------- */
    gtime_c *adjweek(gtime_c t0);
    /* adjust time considering week handover -------------------------------------- */
    gtime_c *adjday(gtime_c t0);
    /* screen data by time -------------------------------------------------------- */
    int screent(gtime_c ts,gtime_c te,double tint);
    /* screen time by time interval ----------------------------------------------- */
    int screen_tinter(double tinter,double dttol);
    /* read leap seconds table by text -------------------------------------------- */
    int read_leaps_text(const string file);
    /* read leap seconds table by usno -------------------------------------------- */
    int read_leaps_usno(const string file);
    /* next download time --------------------------------------------------------- */
    gtime_c *nextdltime(const int *topts,int stat);

    /* copy gtime_c --------------------------------------------------------------- */
    gtime_c *copy_gtime(gtime_c t0);
/* Components */
public:
    time_t time;                    /* time (s) expressed by standard time_t */
    double sec;                     /* fraction of second under 1 s */
    double ep[6];                   /* {year,month,day,hour,min,sec} */
    string sep;                     /* ep in string format */
    int doy;                        /* day of year */
    int Gweek,Rweek,Eweek,Cweek;    /* GPS, GLONASS, Galileo, BDS week */
    int Gwwd,Rwwd,Ewwd,Cwwd;        /* GPS, GLONASS, Galileo, BDS days of week */
    double Gws,Rws,Ews,Cws;         /* GPS, GLONASS, Galileo, BDS seconds of week */
    int sys;                        /* time system */
};

#endif