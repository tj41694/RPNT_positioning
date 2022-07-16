#ifndef GNSSNAV_CONST_H
#define GNSSNAV_CONST_H

/* Physical Constants */
#define PI          3.1415926535897932  /* pi */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define SC2RAD      3.1415926535898     /* semi-circle to radian (IS-GPS) */
#define AU          149597870691.0      /* 1 AU (m) */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define AS2R        (D2R/3600.0)        /* arc sec to radian */
#define DH2RS       (D2R/3600.0)        /* deg/hour to rad/s */
#define PPM2NUM     1E-6                /* ppm to number */
#define FT2M        0.3048              /* ft to meter */
#define ABSZERO     -273.15             /* the absolute zero */
#define SSOUND_SEA  1500.0              /* coarse speed of sound in the sea */

#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */
#define MU_ALL      3.986004418E14      /* WGS84 earth gravitational constant (m^3/s^2) */
#define J2          1.082627E-3         /* WGS84 earth's 2nd gravitational constant */
#define J4         -2.37091222E-6       /* WGS84 earth's 4th gravitational constant */
#define J6          6.08347E-9          /* WGS84 earth's 6th gravitational constant */
#define GA22        1.5                 /* a2's 2nd order gravity coefficient      3/2 */
#define GA24       -1.875               /* a2's 4th order gravity coefficient    -15/8 */
#define GA26        2.1875              /* a2's 6th order gravity coefficient    35/16 */
#define GA32       -4.5                 /* a3's 2nd order gravity coefficient     -9/2 */
#define GA34        18.75               /* a3's 4th order gravity coefficient     75/4 */
#define GA36       -45.9375             /* a3's 6th order gravity coefficient  -735/16 */
#define GA44       -21.875              /* a4's 4th order gravity coefficient   -175/8 */
#define GA46        137.8125            /* a4's 6th order gravity coefficient  2205/16 */
#define GA56       -101.0625            /* a6's 6th order gravity coefficient -1617/16 */
#define GB12        3                   /* b1's 2nd order gravity coefficient        3 */
#define GB14       -7.5                 /* b1's 4th order gravity coefficient    -15/2 */
#define GB16        13.125              /* b1's 6th order gravity coefficient    105/8 */
#define GB24        17.5                /* b2's 4th order gravity coefficient     35/2 */
#define GB26       -78.75               /* b2's 6th order gravity coefficient  -945/12 */
#define GB36        86.625              /* b3's 6th order gravity coefficient    693/8 */
#define GRAVITY     9.80665             /* gravity, g (m/s^2) */
#define GRAVITY_0   9.7803253359        /* gravity 0 (for gravity in Nav frame) */
#define M_GRAVITY   9.80665E-3          /* milli gravity, mg ( m/s^2 ) */

#define WGS84       0                   /* WGS84 datum system */
#define CGCS2000    1                   /* CGCS2000 datum system (China) */
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define RP_WGS84    6356752.31425       /* earth polar radius (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */
#define EE_WGS84    0.0818191908425     /* earth eccentricity (WGS84) */
#define RE_CGCS2000 6378137.0           /* earth semimajor axis (CGCS2000) (m) */
#define RP_CGCS2000 6356752.31414       /* earth polar radius (CGCS2000) (m) */
#define FE_CGCS2000 (1.0/298.257222101) /* earth flattening (CGCS2000) */
#define EE_CGCS2000 0.0818191910428     /* earth eccentricity (CGCS2000) */

#define AXIS_NONE   0x00                /* axis: none */
#define AXIS_X      0x01                /* axis: X */
#define AXIS_Y      0x02                /* axis: Y */
#define AXIS_Z      0x04                /* axis: Z */
#define AXIS_TIME   0x10                /* axis: time */

/* GNSS constants */
#define HION        350000.0            /* ionosphere height (m) */
#define MIN_DIST    300000.0            /* minimum distance between satellite and receiver */

#define FREQ1       1.57542E9           /* L1/E1  frequency (Hz) */
#define FREQ2       1.22760E9           /* L2     frequency (Hz) */
#define FREQ5       1.17645E9           /* L5/E5a frequency (Hz) */
#define FREQ6       1.27875E9           /* E6/L6  frequency (Hz) */
#define FREQ7       1.20714E9           /* E5b    frequency (Hz) */
#define FREQ8       1.191795E9          /* E5a+b  frequency (Hz) */
#define FREQ9       2.492028E9          /* S      frequency (Hz) */
#define FREQ1_GLO   1.60200E9           /* GLONASS G1 base  frequency (Hz) */
#define DFRQ1_GLO   0.56250E6           /* GLONASS G1 bias  frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9           /* GLONASS G2 base  frequency (Hz) */
#define DFRQ2_GLO   0.43750E6           /* GLONASS G2 bias  frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9          /* GLONASS G3       frequency (Hz) */
#define FREQC4_GLO  1.600995E9          /* GLONASS G1a CDMA frequency (Hz) */
#define FREQC6_GLO  1.24806E9           /* GLONASS G2a CDMA frequency (Hz) */
#define FREQ1_BDS   1.57542E9           /* BeiDou B1C/B1A C1   frequency for BDS-3 */
#define FREQ2_BDS   1.561098E9          /* BeiDou B1(I) C2     frequency (Hz) */
#define FREQ5_BDS   1.17645E9           /* BeiDou B2a C5       frequency for BDS-3 */
#define FREQ6_BDS   1.26852E9           /* BeiDou B3(I)/B3A C6 frequency (Hz) */
#define FREQ7_BDS   1.20714E9           /* BeiDou B2(I)/B2b C7 frequency (Hz) */
#define FREQ8_BDS   1.191795E9          /* BeiDou B2a+B2b C8   frequency for BDS-3 */

/* BDS constants */
#define BRDC_AREF_MEO   27906100        /* BDS broadcast A reference for MEO (m) */
#define BRDC_AREF_IG    42162200        /* BDS broadcast A reference for IGSO/GEO (m) */
#define BRDC_MSG_B1C    0x01            /* BDS broadcast message type: B1C */
#define BRDC_MSG_B2b    0x02            /* BDS broadcast message type: B2b */
#define BRDC_MSG_B2a    0x04            /* BDS broadcast message type: B2a */
#define BRDC_MSG_B2I    0x08            /* BDS broadcast message type: B2I */
#define BRDC_MSG_B1I    0x10            /* BDS broadcast message type: B1I */
#define BRDC_MSG_B3I    0x20            /* BDS broadcast message type: B3I */

const double Gt0_ep[] ={ 1980,1, 6,0,0,0 };  /* gps time reference */
const double Et0_ep[] ={ 1999,8,22,0,0,0 };  /* galileo system time reference */
const double Ct0_ep[] ={ 2006,1, 1,0,0,0 };  /* beidou time reference */

/* binary values */
const int binary_value[]={ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x100 };

/* rotation constants */
const double OMG_ie_e[3] ={ 0,0,OMGE }; /* earth rotation vector in ECEF */
const double SSM_OMG_e[9] =             /* Skew-symmetric matrix of OMG_ie_e */
        { 0,-OMGE,0,OMGE,0,0,0,0,0 };

/* math constants */
const double I33[9]={ 1,0,0,0,1,0,0,0,1 };
const int AXIS4[4]={ AXIS_X,AXIS_Y,AXIS_Z,AXIS_TIME };

#endif
