/* Time function source file ---------------------------------------------------------------------- */

#include "BaseFunction/timesys.h"
#include "BaseFunction/basefunction.h"

/* const ------------------------------------------------------------------------------------------ */
/* time reference */
const gtime_c gpst0_(Gt0_ep);
const gtime_c galt0_(Et0_ep);
const gtime_c bdst0_(Ct0_ep);

static double leaps[][MAXLEAPS + 1] ={ /* leap seconds (y,m,d,h,m,s,utc-gpst) */
    { 2017,1,1,0,0,0,-18 },
    { 2015,7,1,0,0,0,-17 },
    { 2012,7,1,0,0,0,-16 },
    { 2009,1,1,0,0,0,-15 },
    { 2006,1,1,0,0,0,-14 },
    { 1999,1,1,0,0,0,-13 },
    { 1997,7,1,0,0,0,-12 },
    { 1996,1,1,0,0,0,-11 },
    { 1994,7,1,0,0,0,-10 },
    { 1993,7,1,0,0,0, -9 },
    { 1992,7,1,0,0,0, -8 },
    { 1991,1,1,0,0,0, -7 },
    { 1990,1,1,0,0,0, -6 },
    { 1988,1,1,0,0,0, -5 },
    { 1985,7,1,0,0,0, -4 },
    { 1983,7,1,0,0,0, -3 },
    { 1982,7,1,0,0,0, -2 },
    { 1981,7,1,0,0,0, -1 },
    { 0   ,0,0,0,0,0,  0 }
};

/* difference with other time ----------------------------------------------------- */
gtime_c::gtime_c() {
    time=0; sec=0.0;
    sys=0; sep="\0";
    for (int i=0; i<6; i++) ep[i]=0.0;
    doy=0;
    Gweek=Rweek=Eweek=Cweek=0;
    Gwwd=Rwwd=Ewwd=Cwwd=0;
    Gws=Rws=Ews=Cws=0;
}
/* initialize with epoch array ---------------------------------------------------- */
gtime_c::gtime_c(const double *epoch) {
    sys=0; sep="\0";
    for (int i=0; i<6; i++) ep[i]=0.0;
    doy=0;
    Gweek=Rweek=Eweek=Cweek=0;
    Gwwd=Rwwd=Ewwd=Cwwd=0;
    Gws=Rws=Ews=Cws=0;
    epoch2time(epoch);
}
gtime_c::gtime_c(int isys) {
    sys=0; sep="\0";
    for (int i=0; i<6; i++) ep[i]=0.0;
    doy=0;
    Gweek=Rweek=Eweek=Cweek=0;
    Gwwd=Rwwd=Ewwd=Cwwd=0;
    Gws=Rws=Ews=Cws=0;
    if (isys==iBDS) epoch2time(Ct0_ep);
    else if (isys==iGAL) epoch2time(Et0_ep);
    else epoch2time(Gt0_ep);
}
gtime_c::~gtime_c() {
}

/* string to time --------------------------------------------------------------------
* convert substring in string to gtime_c struct
* args   : char   *s        I   string ("... yyyy mm dd hh mm ss ...")
----------------------------------------------------------------------------------- */
int gtime_c::str2time(string s) {

    if (sscanf(s.c_str(),"%lf %lf %lf %lf %lf %lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<6)
        return -1;
    if (ep[0]<100) ep[0]+=2000;
    if (ep[0]<=1970||ep[1]==0||ep[2]==0) return -1;

    epoch2time(ep);

    sep=s;

    return 0;
}

/* ep time to string ---------------------------------------------------------------------------------
string sep yyyy/mm/dd hh:mm:ss.ssss...
--------------------------------------------------------------------------------------------------- */
string gtime_c::time2str(int n) {

    if (n<0) n=0; else if (n>12) n=12;
    string str;
    if (1.0-sec<0.5/pow(10.0,n)) { time++; sec=0.0; };
    time2epoch();
    sep=int2str(4,"0",(int)ep[0],str)+"/"+int2str(2,"0",(int)ep[1],str)+"/"+
        int2str(2,"0",(int)ep[2],str)+" "+int2str(2,"0",(int)ep[3],str)+":"+
        int2str(2,"0",(int)ep[4],str)+":"+doul2str(2+n+1,n,"0",ep[5],str);
    return sep;
}

/* calender day/time (ep) to time ----------------------------------------------------------------- */
gtime_c *gtime_c::epoch2time(const double *inep) {
    const int doy[]={ 1,32,60,91,121,152,182,213,244,274,305,335 };

    int days,dsec,year=int(inep[0]),mon=(int)inep[1],day=(int)inep[2];

    if (year<1970||2099<year||mon<1||12<mon) {
        time=0; sec=0; return this;
    }

    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3 ? 1 : 0);
    dsec=(int)floor(inep[5]);
    time=(time_t)days*86400+(time_t)inep[3]*3600+(time_t)inep[4]*60+dsec;
    sec=inep[5]-dsec;

    return this;
}

/* time to calender day/time (ep) --------------------------------------------------------------------
 ep={yyyy,mm,dd,hh,mm,ss.ssss...}
--------------------------------------------------------------------------------------------------- */
void gtime_c::time2epoch() {
    const int mday[]={ /* # of days in a month */
        31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
        31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    int days,dsec,mon,day;

    /* leap year if year%4==0 in 1901-2099 */
    days=(int)(time/86400);
    dsec=(int)(time-(time_t)days*86400);
    for (day=days%1461,mon=0; mon<48; mon++) {
        if (day>=mday[mon]) day-=mday[mon]; else break;
    }
    ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
    ep[3]=dsec/3600; ep[4]=dsec%3600/60; ep[5]=dsec%60+sec;
}

/* gps week time to gtime_c --------------------------------------------------------------------------
* args   : int    week      I   week number in gps time
*          double sss       I   time of week in gps time (s)
--------------------------------------------------------------------------------------------------- */
gtime_c *gtime_c::gpst2time(int week,double sss) {
    epoch2time(Gt0_ep);

    if (sss<-1E9||1E9<sss) sss=0.0;
    time+=(time_t)86400*7*week+(int)sss;
    sec=(long long(sss*1E9)-long long(sss)*1E9)/1.E9;

    Gweek=week;
    Gws=sss;
    Gwwd=int(sss/86400.0);

    return this;
}

/* gtime_c to gps week time --------------------------------------------------------------------------
*  args  : int    *week     IO  week number in gps time (NULL: no output)
--------------------------------------------------------------------------------------------------- */
double gtime_c::time2gpst(int *week) {
    time_t sss=time-gpst0_.time;

    Gweek=(int)(sss/(86400*7));
    Gws=(double)(sss-(double)Gweek*86400*7)+sec;
    Gwwd=(int)((sss-(double)Gweek*86400*7)/86400.0);
    if (week) *week=Gweek;

    return Gws;
}

/* galileo week time to gtime_c ----------------------------------------------------------------------
* args   : int    week      I   week number in gst
*          double sec       I   time of week in gst (s)
--------------------------------------------------------------------------------------------------- */
gtime_c *gtime_c::gst2time(int week,double sss) {
    epoch2time(Et0_ep);

    if (sss<-1E9||1E9<sss) sss=0.0;
    time+=(time_t)86400*7*week+(int)sss;
    sec=(long long(sss*1E9)-long long(sss)*1E9)/1.E9;

    Eweek=week;
    Ews=sss;
    Ewwd=int(sss/86400.0);

    return this;
}

/* gtime_c to galileo week time ----------------------------------------------------------------------
*  args  : int    *week     IO  week number in gst (NULL: no output)
--------------------------------------------------------------------------------------------------- */
double gtime_c::time2gst(int *week) {
    time_t sss=time-galt0_.time;

    Eweek=(int)(sss/(86400*7));
    Ews=(double)(sss-(double)Eweek*86400*7)+sec;
    Ewwd=(int)((sss-(double)Eweek*86400*7)/86400.0);
    if (week) *week=Eweek;

    return Ews;
}

/* Beidou week time to gtime_c -----------------------------------------------------------------------
* args   : int    week      I   week number in bdt
*          double sss       I   time of week in bdt (s)
--------------------------------------------------------------------------------------------------- */
gtime_c *gtime_c::bdt2time(int week,double sss) {
    epoch2time(Ct0_ep);

    if (sss<-1E9||1E9<sss) sss=0.0;
    time+=(time_t)86400*7*week+(int)sss;
    sec=(long long(sss*1E9)-long long(sss)*1E9)/1.E9;

    Cweek=week;
    Cws=sss;
    Cwwd=int(sss/86400.0);

    return this;
}

/* gtime_c to Beidou week time -----------------------------------------------------------------------
* args   : int    *week     IO  week number in bdt (NULL: no output)
--------------------------------------------------------------------------------------------------- */
double gtime_c::time2bdt(int *week) {
    time_t sss=time-bdst0_.time;

    Cweek=(int)(sss/(86400*7));
    Cws=(double)(sss-(double)Cweek*86400*7)+sec;
    Cwwd=(int)((sss-(double)Cweek*86400*7)/86400.0);
    if (week) *week=Cweek;

    return Cws;
}

/* add dsec(s) to gtime_c ------------------------------------------------------------------------- */
gtime_c *gtime_c::timeadd(double dsec) {
    double tt;
    sec+=dsec; tt=floor(sec); time+=(int)tt; sec-=tt;
    return this;
}

/* difference with gtime_c t2 --------------------------------------------------------------------- */
double gtime_c::timediff(const gtime_c &t2)	const
{
    return difftime(time,t2.time)+sec-t2.sec;
}

/* get current time in utc ---------------------------------------------------------------------------
* get current time in utc
* args   : none
* return : current time in utc
*-------------------------------------------------------------------------------------------------- */
static double timeoffset_=0.0;        /* time offset (s) */

gtime_c *gtime_c::timeget() {
#ifdef WIN32
    SYSTEMTIME ts;

    GetSystemTime(&ts); /* utc */
    ep[0]=ts.wYear; ep[1]=ts.wMonth;  ep[2]=ts.wDay;
    ep[3]=ts.wHour; ep[4]=ts.wMinute; ep[5]=ts.wSecond+ts.wMilliseconds*1E-3;
#else
    struct timeval tv;
    struct tm *tt;

    if (!gettimeofday(&tv,NULL)&&(tt=gmtime(&tv.tv_sec))) {
        ep[0]=tt->tm_year+1900; ep[1]=tt->tm_mon+1; ep[2]=tt->tm_mday;
        ep[3]=tt->tm_hour; ep[4]=tt->tm_min; ep[5]=tt->tm_sec+tv.tv_usec*1E-6;
    }
#endif
    epoch2time(ep);

#ifdef CPUTIME_IN_GPST /* cputime operated in gpst */
    gpst2utc();
#endif
    return timeadd(timeoffset_);
}
/* set current time in utc ---------------------------------------------------------------------------
* set current time in utc
* args   : gtime_c          I   current time in utc
* return : none
* notes  : just set time offset between cpu time and current time
*          the time offset is reflected to only timeget()
*          not reentrant
*-------------------------------------------------------------------------------------------------- */
void gtime_c::timeset() {
    gtime_c t0;
    timeoffset_+=timediff(*t0.timeget());
}

/* gpstime to utc ------------------------------------------------------------------------------------
* convert gpstime to utc considering leap seconds
* return : time expressed in utc
* notes  : ignore slight time offset under 100 ns
*-------------------------------------------------------------------------------------------------- */
gtime_c *gtime_c::gpst2utc() {
    gtime_c tu,t0;
    int i;

    for (i=0; leaps[i][0]>0; i++) {
        tu=*this;
        tu.timeadd(leaps[i][6]);
        if (tu.timediff(*t0.epoch2time(leaps[i]))>=0.0) { *this=tu; return this; }
    }
    return this;
}

/* utc to gpstime ------------------------------------------------------------------------------------
* convert utc to gpstime considering leap seconds
* return : time expressed in gpstime
* notes  : ignore slight time offset under 100 ns
*-------------------------------------------------------------------------------------------------- */
gtime_c *gtime_c::utc2gpst() {
    int i;
    gtime_c t0;

    for (i=0; leaps[i][0]>0; i++) {
        if (timediff(*t0.epoch2time(leaps[i]))>=0.0)
            return timeadd(-leaps[i][6]);
    }
    return this;
}

/* gpstime to bdt ------------------------------------------------------------------------------------
* convert gpstime to bdt (beidou navigation satellite system time)
* return : time expressed in bdt
* notes  : ref [8] 3.3, 2006/1/1 00:00 BDT = 2006/1/1 00:00 UTC
*          no leap seconds in BDT
*          ignore slight time offset under 100 ns
*-------------------------------------------------------------------------------------------------- */
gtime_c *gtime_c::gpst2bdt() {
    return timeadd(-14.0);
}
/* bdt to gpstime ------------------------------------------------------------------------------------
* convert bdt (beidou navigation satellite system time) to gpstime
* return : time expressed in gpstime
* notes  : see gpst2bdt()
*-------------------------------------------------------------------------------------------------- */
gtime_c *gtime_c::bdt2gpst() {
    return timeadd(14.0);
}

/* time to day and sec ---------------------------------------------------------------------------- */
double gtime_c::time2sec(gtime_c &day) {
    double sss;
    double ep0[6]={ 0 };
    int i;
    time2epoch();

    sss=ep[3]*3600.0+ep[4]*60.0+ep[5];
    for (i=0; i<3; i++) ep0[i]=ep[i];
    day.epoch2time(ep0);
    return sss;
}

/* utc to gmst ---------------------------------------------------------------------------------------
* convert utc to gmst (Greenwich mean sidereal time)
* args   : gtime_c t        I   time expressed in utc
*          double ut1_utc   I   UT1-UTC (s)
* return : gmst (rad)
*-------------------------------------------------------------------------------------------------- */
double gtime_c::utc2gmst(double ut1_utc) {
    const double ep2000[]={ 2000,1,1,12,0,0 };
    gtime_c tut,tut0,t2000;
    double ut,t1,t2,t3,gmst0,gmst;

    tut=*this; tut.timeadd(ut1_utc);
    ut=tut.time2sec(tut0);
    t1=tut0.timediff(*t2000.epoch2time(ep2000))/86400.0/36525.0;
    t2=t1*t1; t3=t2*t1;
    gmst0=24110.54841+8640184.812866*t1+0.093104*t2-6.2E-6*t3;
    gmst=gmst0+1.002737909350795*ut;

    return fmod(gmst,86400.0)*PI/43200.0; /* 0 <= gmst <= 2*PI */
}

/* day of year to time ---------------------------------------------------------------------------- */
int gtime_c::doy2time(int Year,int Doy) {
    const int doys[]={ 1,32,60,91,121,152,182,213,244,274,305,335 };

    ep[0]=Year; ep[3]=ep[4]=ep[5]=0; doy=Doy;
    /* get month and day */
    for (int i=0; i<12; i++) {
        if (i==11||Doy<(doys[i+1]+(i>0&&Year%4==0?1:0))) {
            ep[1]=i+1;
            ep[2]=Doy-doys[i]-(i>1&&Year%4==0?1:0)+1;
            break;
        }
    }
    epoch2time(ep);
    return 1;
}

/* gtime_c to day of year ----------------------------------------------------------------------------
* convert time to day of year
* return : day of year (days)
*-------------------------------------------------------------------------------------------------- */
double gtime_c::time2doy() {
    double ep0[6]={ 0 };
    gtime_c t0;

    time2epoch();
    ep0[0]=ep[0]; ep0[1]=ep0[2]=1.0; ep0[3]=ep0[4]=ep0[5]=0.0;
    return this->timediff(*t0.epoch2time(ep0))/86400.0+1.0;
}

/* get time of day -----------------------------------------------------------------------------------
* return : time of last day
*-------------------------------------------------------------------------------------------------- */
double gtime_c::time_of_day() {
    time_t sss=time-gpst0_.time;

    double nDay=sss/86400;
    return (double)(sss-nDay*86400)+sec;
}

/* read leap seconds table */
int gtime_c::read_leaps(const string file) {
    int i,n;

    /* read leap seconds table by text or usno */
    if (!(n=read_leaps_text(file))&&!(n=read_leaps_usno(file))) {
        return 0;
    }

    for (i=0; i<7; i++) leaps[n][i]=0.0;

    return 1;
}
/* adjust time considering week handover ---------------------------------- */
gtime_c *gtime_c::adjweek(gtime_c t0) {
    double dt=timediff(t0);
    if (dt < -302400.0) return timeadd(604800.0);
    if (dt >  302400.0) return timeadd(-604800.0);
    return this;
}
/* adjust time considering week handover ---------------------------------- */
gtime_c *gtime_c::adjday(gtime_c t0) {
    double dt=timediff(t0);
    if (dt < -43200.0) return timeadd(86400.0);
    if (dt >  43200.0) return timeadd(-86400.0);
    return this;
}
/* screen by time ------------------------------------------------------------------------------------
* screening by time start, time end, and time interval
* args   :
*		   gtime_c ts    I      time start (ts.time==0:no screening by ts)
*          gtime_c te    I      time end   (te.time==0:no screening by te)
*          double  tint  I      time interval (s) (0.0:no screen by tint)
* return : 1:on condition, 0:not on condition
*-------------------------------------------------------------------------------------------------- */
int gtime_c::screent(gtime_c ts,gtime_c te,double tint) {
    return (tint<=0.0||fmod(time_of_day()+DTTOL,tint)<=DTTOL*2.0)&&
        (ts.time==0||timediff(ts)>=-DTTOL)&&
        (te.time==0||timediff(te)<  DTTOL);
}
/* screen time by time interval ------------------------------------------- */
int gtime_c::screen_tinter(double tinter,double dttol) {
    return (tinter<=0.0||fmod(time_of_day()+dttol,tinter)<=2.0*dttol);
    //return ( tinter<=0.0 || fabs( remainder( remainder(time,tinter)+remainder(sec,tinter), tinter ) )<=dttol );
}
/* read leap seconds table by text -------------------------------------------------------------------
* format : yyyy mm dd hh mm ss ls
--------------------------------------------------------------------------------------------------- */
int gtime_c::read_leaps_text(const string file) {
    ifstream inf;
    string buff;
    int i,n=0,ls,fd;

    inf.open(file,ios::in);
    if (!inf.is_open()) return 0;

    while (getline(inf,buff)&&n<MAXLEAPS) {
        if ((fd=buff.find('#'))!=string::npos) buff.replace(fd,1,"\0");
        if (str2double(buff.substr(0,4),ep[0])==0) continue; /* year */
        for (i=0; i<5; i++)                                  /* mm dd hh mm ss */
            if (str2double(buff.substr(5+3*i,2),ep[i+1])==0) continue;
        if (str2int(buff.substr(20,2),ls)==0) continue;  /* leap second */
        for (i=0; i<6; i++) leaps[n][i]=ep[i];
        leaps[n++][6]=ls;
    }

    inf.close();
    return n;
}
/* read leap seconds table by usno ---------------------------------------------------------------- */
int gtime_c::read_leaps_usno(const string file) {
    static const string months[]={
        "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"
    };
    ifstream inf;
    string buff,month;
    int i,j,y,m,d,n=0;
    double tai_utc,ls[MAXLEAPS][7]={ {0.0} };

    inf.open(file,ios::in);
    if (!inf.is_open()) return 0;

    while (getline(inf,buff)&&n<MAXLEAPS) {
        if (str2int(buff.substr(0,4),y)==0) continue; /* year */
        /* month */
        for (m=0; m<12; m++) if (buff.find(months[m],5)!=string::npos) break;
        if (m>12) continue;
        if (str2int(buff.substr(9,2),d)==0) continue;     /* day */
        /* leap second */
        if ((j=buff.find("TAI-UTC="))==string::npos) continue;
        else str2double(buff.substr(j+8),tai_utc);
        ls[n][0]=y; ls[n][1]=m; ls[n][2]=d; ls[n++][6]=19.0-tai_utc;
    }
    for (i=0; i<n; i++) for (j=0; j<7; j++) leaps[i][j]=ls[n-i-1][j];

    inf.close();
    return n;
}
/* next download time ----------------------------------------------------------------------------- */
gtime_c *gtime_c::nextdltime(const int *topts,int stat) {
    double tow;
    int week,tint;

    /* current time (gpst) */
    timeget()->utc2gpst();
    tow=time2gpst(&week);

    /* next retry time */
    if (stat==0&&topts[3]>0) {
        tow=(floor((tow-topts[2])/topts[3])+1.0)*topts[3]+topts[2];
        return gpst2time(week,tow);
    }

    /* next interval time */
    tint=topts[1]<=0 ? 3600 : topts[1];
    tow=(floor((tow-topts[2])/tint)+1.0)*tint+topts[2];
    gpst2time(week,tow);

    return this;
}

/* copy gtime_c ----------------------------------------------------------- */
gtime_c *gtime_c::copy_gtime(gtime_c t0) {
    time=t0.time;
    sec=t0.sec;
    sep=t0.sep;
    sys=t0.sys;
    for (int i=0; i<6; i++) ep[i]=t0.ep[i];

    return this;
}
