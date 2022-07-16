/* Base dependent library file */
#ifndef GNSSNAV_LIB_H
#define GNSSNAV_LIB_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <cstdarg>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>
#include <deque>
#include <cctype>
#include <cstdio>                       /* for file stream read */
#include <cstring>
#include "gnssnav_const.h"

#ifdef WIN32
#define _WINSOCK_DEPRECATED_NO_WARNINGS
//#include <afx.h>
#include <conio.h>
#include <winsock2.h>
#include <windows.h>
#include <Mmsystem.h>
#pragma comment (lib,"winmm.lib") 
#pragma comment(lib,"ws2_32.lib")
#else
#include <pthread.h>                    /* for Unix System thread */ 
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#ifndef __USE_MISC
#define __USE_MISC
#endif
#ifndef CRTSCTS
#define CRTSCTS  020000000000
#endif
#include <errno.h>
#include <termios.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <netdb.h>
#endif

#ifdef WIN_DLL
#define EXPORT __declspec(dllexport)    /* for Windows DLL */
#else
#define EXPORT
#endif

#ifdef WIN32
#define thread_t    HANDLE
#define lock_t      CRITICAL_SECTION
#define initlock(f) InitializeCriticalSection(f)
#define tolock(f)   EnterCriticalSection(f)
#define tounlock(f) LeaveCriticalSection(f)
#define FILEPATHSEP '\\'
#define thrdID      LPDWORD
#else
#define thread_t    pthread_t
#define lock_t      pthread_mutex_t
#define initlock(f) pthread_mutex_init(f,NULL)
#define tolock(f)   pthread_mutex_lock(f)
#define tounlock(f) pthread_mutex_unlock(f)
#define FILEPATHSEP '/'
#endif

/* define name space */
using std::cout;
using std::cin;
using std::setw;
using std::setprecision;

using std::istringstream;
using std::ostringstream;
using std::ifstream;
using std::fstream;
using std::ios;

using std::string;
using std::to_string;

using std::vector;
using std::deque;

/* algorithm */
using std::max;
using std::any_of;
using std::transform;
using std::binary_search;
using std::sort;

#define VER_THIS   "2022.01a"

/* GNSS defined constants ------------------------------------------------------------------------- */
#define SYS_NONE        0x00                /* navigation system: none */
#define SYS_GPS         0x01                /* navigation system: GPS */
#define SYS_GLO         0x02                /* navigation system: GLONASS */
#define SYS_GAL         0x04                /* navigation system: Galileo */
#define SYS_BDS         0x08                /* navigation system: BeiDou */
#define SYS_QZS         0x10                /* navigation system: QZSS */
#define SYS_IRN         0x20                /* navigation system: IRNS */
#define SYS_SBS         0x40                /* navigation system: SBAS */
#define SYS_LEO         0x80                /* navigation system: LEO */
#define SYS_ALL         0xFF                /* navigation system: all */

#define TSYS_GPS        0                   /* time system: GPS time */
#define TSYS_UTC        1                   /* time system: UTC */
#define TSYS_GLO        2                   /* time system: GLONASS time */
#define TSYS_GAL        3                   /* time system: Galileo time */
#define TSYS_QZS        4                   /* time system: QZSS time */
#define TSYS_BDS        5                   /* time system: BeiDou time */
#define TSYS_IRN        6                   /* time system: IRNSS time */

#define EFACT_GPS       1.0                 /* error factor: GPS */
#define EFACT_GLO       1.5                 /* error factor: GLONASS */
#define EFACT_GAL       1.0                 /* error factor: Galileo */
#define EFACT_QZS       1.0                 /* error factor: QZSS */
#define EFACT_BDS       1.0                 /* error factor: BeiDou */
#define EFACT_BDS_G     3.0                 /* error factor: BeiDou GEO */
#define EFACT_IRN       1.5                 /* error factor: IRNSS */
#define EFACT_SBS       3.0                 /* error factor: SBAS */

#define MAXFREQ         9                   /* max NFREQ (1:9) */
#define MAX_NF          6                   /* maximum number of frequencies in one constellation */
#ifndef NFREQ
#define NFREQ           3                   /* number of carrier frequencies */
#endif
#define NFREQGLO        2                   /* number of carrier frequencies of GLONASS */

#ifndef NEXOBS
#define NEXOBS          0                   /* number of extended obs codes */
#endif

#define MINPRNGPS       1                   /* min satellite PRN number of GPS */
#define MAXPRNGPS       32                  /* max satellite PRN number of GPS */
#define NSATGPS         (MAXPRNGPS-MINPRNGPS+1) /* number of GPS satellites */
#define NSYSGPS         1
#define iGPS            0                   /* index of GPS in GNSS process */

#define MINPRNGLO       1                   /* min satellite slot number of GLONASS */
#define MAXPRNGLO       30                  /* max satellite slot number of GLONASS */
#define NSATGLO         (MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO         1
#define iGLO            1                   /* index of GLONASS in GNSS process */

//#ifdef ENAGAL
#define MINPRNGAL       1                   /* min satellite PRN number of Galileo */
#define MAXPRNGAL       40                  /* max satellite PRN number of Galileo */
#define NSATGAL         (MAXPRNGAL-MINPRNGAL+1) /* number of Galileo satellites */
#define NSYSGAL         1
#define iGAL            2                   /* index of Galileo in GNSS process */

//#ifdef ENABDS
#define MINPRNBDS       1                   /* min satellite sat number of BeiDou */
#define MAXPRNBDS       61                  /* max satellite sat number of BeiDou */
#define NSATBDS         (MAXPRNBDS-MINPRNBDS+1) /* number of BeiDou satellites */
#define NSYSBDS         1
#define iBDS            3                   /* index of BDS in GNSS process */
#define BDS_2ND         0x01                /* BDS-2nd generation */
#define BDS_3RD         0x02                /* BDS-3rd generation */
#define BDS_ORB_MEO     0x01                /* BDS MEO orbit */
#define BDS_ORB_IGSO    0x02                /* BDS IGSO orbit */
#define BDS_ORB_GEO     0x04                /* BDS GEO orbit */
#define BDS3_MIN        19                  /* min BDS-3 prn */

#ifdef ENAQZS
#define MINPRNQZS       193                 /* min satellite PRN number of QZSS */
#define MAXPRNQZS       199                 /* max satellite PRN number of QZSS */
#define MINPRNQZS_S     183                 /* min satellite PRN number of QZSS SAIF */
#define MAXPRNQZS_S     189                 /* max satellite PRN number of QZSS SAIF */
#define NSATQZS         (MAXPRNQZS-MINPRNQZS+1) /* number of QZSS satellites */
#define iQZS            4                   /* index of QZSS in GNSS process */
#define NSYSQZS         1
#else
#define MINPRNQZS       0
#define MAXPRNQZS       0
#define MINPRNQZS_S     0
#define MAXPRNQZS_S     0
#define NSATQZS         0
#define iQZS            4                   /* index of QZSS in GNSS process */
#define NSYSQZS         0
#endif

#ifdef ENAIRN
#define MINPRNIRN       1                   /* min satellite sat number of IRNSS */
#define MAXPRNIRN       7                   /* max satellite sat number of IRNSS */
#define NSATIRN         (MAXPRNIRN-MINPRNIRN+1) /* number of IRNSS satellites */
#define iIRN            5                   /* index of IRNSS in GNSS process */
#define NSYSIRN         1
#else
#define MINPRNIRN       0
#define MAXPRNIRN       0
#define NSATIRN         0
#define iIRN            5                   /* index of IRNSS in GNSS process */
#define NSYSIRN         0
#endif

#ifdef ENASBS
#define MINPRNSBS       101                 /* min satellite PRN number of SBAS */
#define MAXPRNSBS       142                 /* max satellite PRN number of SBAS */
#define NSATSBS         (MAXPRNSBS-MINPRNSBS+1) /* number of SBAS satellites */
#define iSBS            6                   /* index of SBAS in GNSS process */
#else
#define MINPRNSBS       0                   /* min satellite PRN number of SBAS */
#define MAXPRNSBS       0                   /* max satellite PRN number of SBAS */
#define NSATSBS         0                   /* number of SBAS satellites */
#define iSBS            6                   /* index of SBAS in GNSS process */
#define NSYSSBS         0
#endif

#ifdef ENALEO
#define MINPRNLEO       1                   /* min satellite sat number of LEO */
#define MAXPRNLEO       10                  /* max satellite sat number of LEO */
#define NSATLEO         (MAXPRNLEO-MINPRNLEO+1) /* number of LEO satellites */
#define NSYSLEO         1
#else
#define MINPRNLEO       0
#define MAXPRNLEO       0
#define NSATLEO         0
#define NSYSLEO         0
#endif

#define MAX_GNSS        7                   /* maximum number of systems */
#define NSYS            (NSYSGPS+NSYSGLO+NSYSGAL+NSYSBDS+NSYSQZS+NSYSIRN+NSYSSBS+NSYSLEO) /* number of systems */

#define MAXSAT          (NSATGPS+NSATGLO+NSATGAL+NSATBDS+NSATQZS+NSATIRN+NSATSBS+NSATLEO) /* number of satellites */
                                            /* max satellite number (1 to MAXSAT) */
#define MAXSTA          255

#ifndef MAXOBS
#define MAXOBS          128                 /* max number of obs in an epoch */
#endif
#define MINSYSOBS       1                   /* min available satellites number of each system in positioning */
#define MINSYSOBS_RTK   2                   /* min available satellites number of each system in RTK positioning */
#define MAXRCV          64                  /* max receiver number (1 to MAXRCV) */
#define MAXOBSTYPE      64                  /* max number of obs type in RINEX */
#define DTTOL           0.005               /* tolerance of time difference (s) */
#define MAXDTOE         7200.0              /* max time difference to GPS Toe (s) */
#define MAXDTOE_QZS     7200.0              /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL     10800.0             /* max time difference to Galileo Toe (s) */
#define MAXDTOE_BDS     21600.0             /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO     1800.0              /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_SBS     360.0               /* max time difference to SBAS Toe (s) */
#define MAXDTOE_S       86400.0             /* max time difference to ephem toe (s) for other */
#define MAXGDOP         300.0               /* max GDOP */

#define INT_SWAP_TRAC   86400.0             /* swap interval of trace file (s) */
#define INT_SWAP_STAT   86400.0             /* swap interval of solution status file (s) */

#define MAXRNAMLEN      20
#define MAXNCON         10                  /* max number of input config file*/
#define MAXEXFILE       1024                /* max number of expanded files */
#define MAXSBSAGEF      30.0                /* max age of SBAS fast correction (s) */
#define MAXSBSAGEL      1800.0              /* max age of SBAS long term corr (s) */
#define MAXSBSURA       8                   /* max URA of SBAS satellite */
#define MAXBAND         10                  /* max SBAS band of IGP */
#define MAXNIGP         201                 /* max number of IGP in SBAS band */
#define MAXNGEO         4                   /* max number of GEO satellites */
#define MAXCOMMENT      10                  /* max number of RINEX comments */
#define MAXSTRPATH      1024                /* max length of stream path */
#define MAXSTRMSG       1024                /* max length of stream message */
#define MAXSTRRTK       8                   /* max number of stream in RTK server */
#define MAXSBSMSG       32                  /* max number of SBAS msg in RTK server */
#define MAXSOLMSG       8191                /* max length of solution message */
#define MAXRAWLEN       4096                /* max length of receiver raw message */
#define MAXERRMSG       4096                /* max length of error/warning message */
#define MAXANT          64                  /* max length of station name/antenna type */
#define MAXMEASBUF      3                   /* number of observation buff (must>2) */
#define MAXSOLBUF       3                   /* max number of solution buffer (must>=2) */
#define MAXOBSBUF       128                 /* max number of observation data buffer */
#define MAXNRPOS        16                  /* max number of reference positions */
#define MAXLEAPS        64                  /* max number of leap seconds table */
#define MAXGISLAYER     32                  /* max number of GIS data layers */
#define MAXRCVCMD       4096                /* max length of receiver commands */

#define RNX2VER         2.10                /* RINEX ver.2 default output version */
#define RNX3VER         3.00                /* RINEX ver.3 default output version */

#define SRNX_LEN        12                  /* name length of RINEX file in short format */
#define SSATCLK_LEN     12                  /* name length of precise orbit and clock file in short format */
#define SDSNX_LEN       12                  /* name length of SNX daily solution file in short format */
#define SDYSNX_LEN      15                  /* name length of SNX daily solution file with year in short format */
#define SWSNX_LEN       11                  /* name length of SNX weekly solution file in short format */
#define SWYSNX_LEN      14                  /* name length of SNX weekly solution file with year in short format */
#define LRNXO_LEN       38                  /* name length of RINEX observation file in long format */
#define LRNXN_LEN       34                  /* name length of RINEX navigation file in long format */
#define LSATCLK_LEN     38                  /* name length of precise orbit and clock file in long format */
#define LSNX_LEN        38                  /* name length of SNX daily/weekly solution file in long format */
#define LBIA_LEN        38                  /* name length of bias file in long format */

#define OBSTYPE_PR      0x01                /* observation type: pseudorange */
#define OBSTYPE_CP      0x02                /* observation type: carrier-phase */
#define OBSTYPE_DOP     0x04                /* observation type: doppler-freq */
#define OBSTYPE_SNR     0x08                /* observation type: SNR */
#define OBSTYPE_ALL     0xFF                /* observation type: all */
#define OBSTYPE_CL      0x03                /* observation type: pseudorange + carrier-phase */

#define FREQTYPE_L1     0x01                /* frequency type: L1/E1 */
#define FREQTYPE_L2     0x02                /* frequency type: L2/B1 */
#define FREQTYPE_L5     0x04                /* frequency type: L5/E5a/L3 */
#define FREQTYPE_L6     0x08                /* frequency type: E6/LEX/B3 */
#define FREQTYPE_L7     0x10                /* frequency type: E5b/B2 */
#define FREQTYPE_L8     0x20                /* frequency type: E5(a+b) */
#define FREQTYPE_L9     0x40                /* frequency type: S */
#define FREQTYPE_ALL    0xFF                /* frequency type: all */

#define FREQ_L1         1                   /* frequency type: L1/E1/B1(A/C) */
#define FREQ_L2         2                   /* frequency type: L2/B1 */
#define FREQ_L3         3                   /* frequency type: G3 */
#define FREQ_L4         4                   /* frequency type: G1a */
#define FREQ_L5         5                   /* frequency type: L5/E5a/B2a */
#define FREQ_L6         6                   /* frequency type: L6/G2a/B3 */
#define FREQ_L7         7                   /* frequency type: E5b/B2b */
#define FREQ_L8         8                   /* frequency type: E5(a+b)/B2(a+b) */
#define FREQ_L9         9                   /* frequency type: S */

#define MAX_CODE_TYPE   4                   /* maximum code type number (C,L,D,S) */
#define CODE_NONE       0                   /* obs code: none or unknown */
#define CODE_L1C        1                   /* obs code: L1C/A,G1C/A,E1C            (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P        2                   /* obs code: L1P,G1P,B1CP               (GPS,GLO,BDS) */
#define CODE_L1W        3                   /* obs code: L1 Z-track                 (GPS) */
#define CODE_L1Y        4                   /* obs code: L1Y                        (GPS) */
#define CODE_L1M        5                   /* obs code: L1M                        (GPS) */
#define CODE_L1N        6                   /* obs code: L1codeless                 (GPS) */
#define CODE_L1S        7                   /* obs code: L1C(D),B1AD                (GPS,BDS,QZS) */
#define CODE_L1L        8                   /* obs code: L1C(P),B1AP                (GPS,BDS,QZS) */
#define CODE_L1D        9                   /* obs code: B1CD (BDS) */
#define CODE_L1A        10                  /* obs code: E1A                        (GAL) */
#define CODE_L1B        11                  /* obs code: E1B,L1Sb                   (GAL,QZS) */
#define CODE_L1X        12                  /* obs code: L1C(D+P),E1B+C,B1CD+P      (GPS,GAL,BDS,QZS) */
#define CODE_L1Z        13                  /* obs code: E1A+B+C,B1AD+P,L1SAIF      (GAL,BDS,QZS) */
#define CODE_L1I        14                  /* (not used) */
#define CODE_L1Q        15                  /* (not used) */
#define CODE_L2C        16                  /* obs code: L2C/A,G1C/A                (GPS,GLO) */
#define CODE_L2D        17                  /* obs code: L2 L1C/A-(P2-P1)           (GPS) */
#define CODE_L2S        18                  /* obs code: L2C(M)                     (GPS,QZS) */
#define CODE_L2L        19                  /* obs code: L2C(L)                     (GPS,QZS) */
#define CODE_L2X        20                  /* obs code: L2C(M+L),B1I+Q             (GPS,BDS,QZS) */
#define CODE_L2P        21                  /* obs code: L2P,G2P                    (GPS,GLO) */
#define CODE_L2W        22                  /* obs code: L2 Z-track                 (GPS) */
#define CODE_L2Y        23                  /* obs code: L2Y                        (GPS) */
#define CODE_L2M        24                  /* obs code: L2M                        (GPS) */
#define CODE_L2N        25                  /* obs code: L2codeless                 (GPS) */
#define CODE_L2I        26                  /* obs code: B1I                        (BDS) */
#define CODE_L2Q        27                  /* obs code: B1Q                        (BDS) */
#define CODE_L3I        28                  /* obs code: G3I                        (GLO) */
#define CODE_L3Q        29                  /* obs code: G3Q                        (GLO) */
#define CODE_L3X        30                  /* obs code: G3I+Q                      (GLO) */
#define CODE_L4A        31                  /* obs code: G4A                        (GLO) */
#define CODE_L4B        32                  /* obs code: G4B                        (GLO) */
#define CODE_L4X        33                  /* obs code: G4X                        (GLO) */
#define CODE_L5I        34                  /* obs code: L5I,E5aI                   (GPS,GAL,QZS,SBS) */
#define CODE_L5Q        35                  /* obs code: L5Q,E5aQ                   (GPS,GAL,QZS,SBS) */
#define CODE_L5X        36                  /* obs code: L5I+Q,E5aI+Q,B2aD+P,L5B+C  (GPS,GAL,BDS,QZS,IRN,SBS) */
#define CODE_L5A        37                  /* obs code: L5A SPS                    (IRN) */
#define CODE_L5B        38                  /* obs code: L5B RS(D)                  (IRN) */
#define CODE_L5C        39                  /* obs code: L5C RS(P)                  (IRN) */
#define CODE_L5D        40                  /* obs code: B2aD,L5S(I)                (BDS,QZS) */
#define CODE_L5P        41                  /* obs code: B2aP,L5S(Q)                (BDS,QZS) */
#define CODE_L5Z        42                  /* obs code: L5S(I+Q)                   (QZS) */
#define CODE_L6A        43                  /* obs code: E6A                        (GAL) */
#define CODE_L6B        44                  /* obs code: E6B                        (GAL) */
#define CODE_L6C        45                  /* obs code: E6C                        (GAL) */
#define CODE_L6X        46                  /* obs code: E6B+C,B3I+Q,L6(D+P)        (GAL,BDS,QZS) */
#define CODE_L6Z        47                  /* obs code: E6A+B+C,B3D+P,L6D+E        (GAL,BDS,QZS) */
#define CODE_L6S        48                  /* obs code: L6D                        (QZS) */
#define CODE_L6L        49                  /* obs code: L6P                        (QZS) */
#define CODE_L6E        50                  /* obs code: L6E                        (QZS) */
#define CODE_L6I        51                  /* obs code: B3I                        (BDS) */
#define CODE_L6Q        52                  /* obs code: B3Q                        (BDS) */
#define CODE_L6D        53                  /* obs code: B3D                        (BDS) */
#define CODE_L6P        54                  /* obs code: B3P                        (BDS) */
#define CODE_L7I        55                  /* obs code: E5bI,B2I                   (GAL,BDS) */
#define CODE_L7Q        56                  /* obs code: E5bQ,B2Q                   (GAL,BDS) */
#define CODE_L7X        57                  /* obs code: E5bI+Q,B2I+Q               (GAL,BDS) */
#define CODE_L7P        58                  /* obs code: E5bI,B2bD                  (GAL,BDS) */
#define CODE_L7D        59                  /* obs code: E5bQ,B2bP                  (GAL,BDS) */
#define CODE_L7Z        60                  /* obs code: E5bI+Q,B2bD+P              (GAL,BDS) */
#define CODE_L8I        61                  /* obs code: E5(a+b)I                   (GAL) */
#define CODE_L8Q        62                  /* obs code: E5(a+b)Q                   (GAL) */
#define CODE_L8X        63                  /* obs code: E5(a+b)I+Q,B2(a+b)D+P      (GAL,BDS) */
#define CODE_L8D        64                  /* obs code: B2(a+b)D                   (BDS) */
#define CODE_L8P        65                  /* obs code: B2(a+b)P                   (BDS) */
#define CODE_L9A        66                  /* obs code: SA SPS                     (IRN) */
#define CODE_L9B        67                  /* obs code: SB RS(D)                   (IRN) */
#define CODE_L9C        68                  /* obs code: SC RS(P)                   (IRN) */
#define CODE_L9X        69                  /* obs code: SB+C                       (IRN) */
#define MAXCODE         69                  /* max number of obs code */

#define PMODE_SINGLE     0              /* positioning mode: single */
#define PMODE_PPP_KINEMA 1              /* positioning mode: PPP-kinemaric */
#define PMODE_PPP_STATIC 2              /* positioning mode: PPP-static */
#define PMODE_PPP_FIXED  3              /* positioning mode: PPP-fixed */
#define PMODE_DGPS       4              /* positioning mode: DGPS/DGNSS */
#define PMODE_KINEMA     5              /* positioning mode: kinematic */
#define PMODE_STATIC     6              /* positioning mode: static */
#define PMODE_FIXED      7              /* positioning mode: fixed */
#define PMODE_MOVEB      8              /* positioning mode: moving-base */

#define PPP_UN_DIFF      0              /* PPP mode: undifferenced */
#define PPP_SINGLE_DIFF  1              /* PPP mode: single-differenced */

#define GNSSLOG_NONE     0x00           /* log message : none */
#define GNSSLOG_NSAT     0x01           /* log message : satellites number */
#define GNSSLOG_PDOP     0x02           /* log message : PDOP values */
#define GNSSLOG_SNR      0x04           /* log message : average SNR */
#define GNSSLOG_MPATH    0x08           /* log message : average multi-path */
#define GNSSLOG_EXCSAT   0x10           /* log message : excluded satellite message */
#define GNSSLOG_HELMERT  0x20           /* log message : helmert filter messages */
#define GNSSLOG_XPAR     0x1000         /* log test message : solved parameter */
#define GNSSLOG_DXPAR    0x2000         /* log test message : parameter correction (dx) */

#define ADJUST_LSA        0             /* adjustment mode: least square adjustment */
#define ADJUST_KALMAN     1             /* adjustment mode: kalman filter */
#define ADJUST_HELMERT    2             /* adjustment mode: kalman filter with helmert */

#define ADJUST_ROBUST_NONE      0       /* adjustment robust mode: none */
#define ADJUST_ROBUST_IGG3_2    1       /* adjustment robust mode: 2-interval IGG3 robust function */
#define ADJUST_ROBUST_IGG3_3    2       /* adjustment robust mode: 2-interval IGG3 robust function */

#define ROBUST_R_DIAGONAL       0       /* covariance type in robust function: only variance (matrix diagonal) */
#define ROBUST_R_ALL            1       /* covariance type in robust function: all covariance (matrix) */
#define ROBUST_R_VECTOR         2       /* covariance type in robust function: variance vector */

#define SOLF_LLH    0                   /* solution format: lat/lon/height */
#define SOLF_XYZ    1                   /* solution format: x/y/z-ecef */
#define SOLF_ENU    2                   /* solution format: e/n/u-baseline */
#define SOLF_NMEA   3                   /* solution format: NMEA-183 */
#define SOLF_STAT   4                   /* solution format: solution status */
#define SOLF_GSIF   5                   /* solution format: GSI F1/F2 */

#define SOLQ_NONE   0                   /* solution status: no solution */
#define SOLQ_FIX    1                   /* solution status: fix */
#define SOLQ_FLOAT  2                   /* solution status: float */
#define SOLQ_SBAS   3                   /* solution status: SBAS */
#define SOLQ_DGPS   4                   /* solution status: DGPS/DGNSS */
#define SOLQ_SINGLE 5                   /* solution status: single */
#define SOLQ_PPP    6                   /* solution status: PPP */
#define SOLQ_DR     7                   /* solution status: dead reconing */
#define MAXSOLQ     7                   /* max number of solution status */

#define VSOLQ_NONE  0                   /* velocity solution: no solution */
#define VSOLQ_UD    1                   /* velocity solution: un-differenced doppler */
#define VSOLQ_SD    2                   /* velocity solution: station-single-differenced doppler */

#define SPPNX       NSYS+3              /* position: maximum number of SPP parameters */
#define DPLNX       4                   /* doppler: number of estimated parameters */

#define TIMES_GPST  0                   /* time system: gps time */
#define TIMES_UTC   1                   /* time system: utc */
#define TIMES_BDT   2                   /* time system: bdt */

#define IONOOPT_OFF    0                /* ionosphere option: correction off */
#define IONOOPT_BRDC   1                /* ionosphere option: broadcast model */
#define IONOOPT_SBAS   2                /* ionosphere option: SBAS model */
#define IONOOPT_IFLC   3                /* ionosphere option: L1/L2 or L1/L5 iono-free LC */
#define IONOOPT_CONST  4                /* ionosphere option: constrained ION model (estimation) */
#define IONOOPT_TEC    5                /* ionosphere option: IONEX TEC model */
#define IONOOPT_QZS    6                /* ionosphere option: QZSS broadcast model */
#define IONOOPT_LEX    7                /* ionosphere option: QZSS LEX ionospehre */
#define IONOOPT_STEC   8                /* ionosphere option: SLANT TEC model */

#define TROPOPT_OFF    0                /* troposphere option: correction off */
#define TROPOPT_SAAS   1                /* troposphere option: Saastamoinen model */
#define TROPOPT_SBAS   2                /* troposphere option: SBAS model */
#define TROPOPT_EST    3                /* troposphere option: ZTD estimation */
#define TROPOPT_ESTG   4                /* troposphere option: ZTD+grad estimation */
#define TROPOPT_ZTD    5                /* troposphere option: ZTD correction */

#define EPHOPT_BRDC    0                /* ephemeris option: broadcast ephemeris */
#define EPHOPT_PREC    1                /* ephemeris option: precise ephemeris */
#define EPHOPT_SBAS    2                /* ephemeris option: broadcast + SBAS */
#define EPHOPT_SSRAPC  3                /* ephemeris option: broadcast + SSR_APC */
#define EPHOPT_SSRCOM  4                /* ephemeris option: broadcast + SSR_COM */
#define EPHOPT_LEX     5                /* ephemeris option: QZSS LEX ephemeris */

#define ISBMODE_WHITENOISE  0           /* receiver ISB mode: white noise */
#define ISBMODE_RANDOMWALK  1           /* receiver ISB mode: random walk */

#define ARMODE_OFF     0                /* AR mode: off */
#define ARMODE_CONT    1                /* AR mode: continuous */
#define ARMODE_INST    2                /* AR mode: instantaneous */
#define ARMODE_FIXHOLD 3                /* AR mode: fix and hold */
#define ARMODE_LCWN    4                /* AR mode: LC with wide/narrow-lane */

#define CSMODE_OBS     0                /* cycle-slip mode: observation LLI */
#define CSMODE_CP      1                /* cycle-slip mode: code - phase */
#define CSMODE_GF      2                /* cycle-slip mode: geo-free combination */
#define CSMODE_MW      4                /* cycle-slip mode: Melbourne-Wubbena combination */

#define SBSOPT_LCORR   1                /* SBAS option: long term correction */
#define SBSOPT_FCORR   2                /* SBAS option: fast correction */
#define SBSOPT_ICORR   4                /* SBAS option: ionosphere correction */
#define SBSOPT_RANGE   8                /* SBAS option: ranging */

#define POSOPT_LLH     0                /* pos option: LLH */
#define POSOPT_XYZ     1                /* pos option: XYZ */
#define POSOPT_SINGLE  2                /* pos option: average of single pos */
#define POSOPT_FILE    3                /* pos option: read from pos file */
#define POSOPT_SNX     4                /* pos option: read from SNX file */
#define POSOPT_RINEX   5                /* pos option: rinex header pos */
#define POSOPT_RTCM    6                /* pos option: rtcm station pos */
#define POSOPT_RAW     7                /* pos option: raw station pos */

#define REFSAT_BASE    0                /* reference station: base */
#define REFSAT_ROVER   1                /* reference station: rover */

#define STR_NONE     0                  /* stream type: none */
#define STR_SERIAL   1                  /* stream type: serial */
#define STR_FILE     2                  /* stream type: file */
#define STR_TCPSVR   3                  /* stream type: TCP server */
#define STR_TCPCLI   4                  /* stream type: TCP client */
#define STR_NTRIPSVR 6                  /* stream type: NTRIP server */
#define STR_NTRIPCLI 7                  /* stream type: NTRIP client */
#define STR_FTP      8                  /* stream type: ftp */
#define STR_HTTP     9                  /* stream type: http */
#define STR_NTRIPC_S 10                 /* stream type: NTRIP caster server */
#define STR_NTRIPC_C 11                 /* stream type: NTRIP caster client */
#define STR_UDPSVR   12                 /* stream type: UDP server */
#define STR_UDPCLI   13                 /* stream type: UDP server */
#define STR_MEMBUF   14                 /* stream type: memory buffer */

#define STRFMT_RTCM2 0                  /* stream format: RTCM 2 */
#define STRFMT_RTCM3 1                  /* stream format: RTCM 3 */
#define STRFMT_OEM4  2                  /* stream format: NovAtel OEMV/4 */
#define STRFMT_OEM3  3                  /* stream format: NovAtel OEM3 */
#define STRFMT_UBX   4                  /* stream format: u-blox LEA-*T */
#define STRFMT_SS2   5                  /* stream format: NovAtel Superstar II */
#define STRFMT_CRES  6                  /* stream format: Hemisphere */
#define STRFMT_STQ   7                  /* stream format: SkyTraq S1315F */
#define STRFMT_GW10  8                  /* stream format: Furuno GW10 */
#define STRFMT_JAVAD 9                  /* stream format: JAVAD GRIL/GREIS */
#define STRFMT_NVS   10                 /* stream format: NVS NVC08C */
#define STRFMT_BINEX 11                 /* stream format: BINEX */
#define STRFMT_RT17  12                 /* stream format: Trimble RT17 */
#define STRFMT_SEPT  13                 /* stream format: Septentrio */
#define STRFMT_CMR   14                 /* stream format: CMR/CMR+ */
#define STRFMT_LEXR  15                 /* stream format: Furuno LPY-10000 */
#define STRFMT_RINEX 16                 /* stream format: RINEX */
#define STRFMT_SP3   17                 /* stream format: SP3 */
#define STRFMT_RNXCLK 18                /* stream format: RINEX CLK */
#define STRFMT_SBAS  19                 /* stream format: SBAS messages */
#define STRFMT_NMEA  20                 /* stream format: NMEA 0183 */
#ifndef EXTLEX
#define MAXRCVFMT    14                 /* max number of receiver format */
#else
#define MAXRCVFMT    15
#endif

#define STR_MODE_R  0x1                 /* stream mode: read */
#define STR_MODE_W  0x2                 /* stream mode: write */
#define STR_MODE_RW 0x3                 /* stream mode: read/write */

#define GEOID_EMBEDDED    0             /* geoid model: embedded geoid */
#define GEOID_EGM96_M150  1             /* geoid model: EGM96 15x15" */
#define GEOID_EGM2008_M25 2             /* geoid model: EGM2008 2.5x2.5" */
#define GEOID_EGM2008_M10 3             /* geoid model: EGM2008 1.0x1.0" */
#define GEOID_GSI2000_M15 4             /* geoid model: GSI geoid 2000 1.0x1.5" */
#define GEOID_RAF09       5             /* geoid model: IGN RAF09 for France 1.5"x2" */

#define COMMENTH    "%"                 /* comment line indicator for solution */
#define MSG_DISCONN "$_DISCONNECT\r\n"  /* disconnect message */
#define PRINT_LEN    4096               /* size of print buff */
#define NOBS_BUFF    2880               /* number of gnss_obs_c buff */

#define DLOPT_FORCE   0x01              /* download option: force download existing */
#define DLOPT_KEEPCMP 0x02              /* download option: keep compressed file */
#define DLOPT_HOLDERR 0x04              /* download option: hold on error file */
#define DLOPT_HOLDLST 0x08              /* download option: hold on listing file */

#define P2_5        0.03125             /* 2^-5 */
#define P2_6        0.015625            /* 2^-6 */
#define P2_11       4.882812500000000E-04 /* 2^-11 */
#define P2_15       3.051757812500000E-05 /* 2^-15 */
#define P2_17       7.629394531250000E-06 /* 2^-17 */
#define P2_19       1.907348632812500E-06 /* 2^-19 */
#define P2_20       9.536743164062500E-07 /* 2^-20 */
#define P2_21       4.768371582031250E-07 /* 2^-21 */
#define P2_23       1.192092895507810E-07 /* 2^-23 */
#define P2_24       5.960464477539063E-08 /* 2^-24 */
#define P2_27       7.450580596923828E-09 /* 2^-27 */
#define P2_29       1.862645149230957E-09 /* 2^-29 */
#define P2_30       9.313225746154785E-10 /* 2^-30 */
#define P2_31       4.656612873077393E-10 /* 2^-31 */
#define P2_32       2.328306436538696E-10 /* 2^-32 */
#define P2_33       1.164153218269348E-10 /* 2^-33 */
#define P2_35       2.910383045673370E-11 /* 2^-35 */
#define P2_38       3.637978807091710E-12 /* 2^-38 */
#define P2_39       1.818989403545856E-12 /* 2^-39 */
#define P2_40       9.094947017729280E-13 /* 2^-40 */
#define P2_43       1.136868377216160E-13 /* 2^-43 */
#define P2_48       3.552713678800501E-15 /* 2^-48 */
#define P2_50       8.881784197001252E-16 /* 2^-50 */
#define P2_55       2.775557561562891E-17 /* 2^-55 */

/* BeiDou satellites GEO/MEO/IGSO prn vector */
#define NBDS_MEO    29                  /* number of BDS MEO satellites */
#define NBDS_IGSO   11                  /* number of BDS IGSO satellites */
#define NBDS_GEO    9                   /* number of BDS GEO satellites */
const int BDS_GEO[]  ={ 1,2,3,4,5,18,59,60,61 };
const int BDS_MEO[]  ={ 11,12,14,19,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,36,37,
                        41,42,43,44,45,46,57,58 };
const int BDS_IGSO[] ={ 6,7,8,9,10,13,15,16,38,39,40 };
const int BDS_3[]    ={ 19,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,36,37,38,39,40,
                        41,42,43,44,45,46,47,59,60,61 };
/* gnss frequency matrix */
const double GNSS_FREQ_Hz[MAX_GNSS][MAXFREQ+1]{
    // GPS
    { 0,     FREQ1,     FREQ2,     0,         0,          FREQ5,     0,          0,         0,         0     },
    // GLONASS
    { 0,     FREQ1_GLO, FREQ2_GLO, FREQ3_GLO, FREQC4_GLO, 0,         FREQC6_GLO, 0,         0,         0     },
    // Galileo
    { 0,     FREQ1,     0,         0,         0,          FREQ5,     FREQ6,      FREQ7,     FREQ8,     0     },
    // BDS
    { 0,     FREQ1_BDS, FREQ2_BDS, 0,         0,          FREQ5_BDS, FREQ6_BDS,  FREQ7_BDS, FREQ8_BDS, 0     },
    // QZSS
    { 0,     FREQ1,     FREQ2,     0,         0,          FREQ5,     FREQ6,      0,         0,         0     },
    // IRNSS
    { 0,     0,         0,         0,         0,          FREQ5,     0,          0,         0,         FREQ9 },
    // SBAS
    { 0,     FREQ1,     FREQ2,     0,         0,          FREQ5,     0,          0,         0,         0     }
};
/* gnss lambda matrix */
const double GNSS_LAMBDA_M[MAX_GNSS][MAXFREQ+1]{
    // GPS
    { 0,     CLIGHT/FREQ1,      CLIGHT/FREQ2,     0,                0,                 CLIGHT/FREQ5,   //0-5
             0,                 0,                0,                0            },                    //6-9
    // GLONASS
    { 0,     CLIGHT/FREQ1_GLO,  CLIGHT/FREQ2_GLO, CLIGHT/FREQ3_GLO, CLIGHT/FREQC4_GLO, 0,              //0-5
             CLIGHT/FREQC6_GLO, 0,                0,                0            },                    //6-9
    // Galileo
    { 0,     CLIGHT/FREQ1,      0,                0,                0,                 CLIGHT/FREQ5,   //0-5
             CLIGHT/FREQ6,      CLIGHT/FREQ7,     CLIGHT/FREQ8,     0            },                    //6-9
    // BDS
    { 0,     CLIGHT/FREQ1_BDS,  CLIGHT/FREQ2_BDS, 0,                0,                 CLIGHT/FREQ5_BDS,//0-5
             CLIGHT/FREQ6_BDS,  CLIGHT/FREQ7_BDS, CLIGHT/FREQ8_BDS, 0            },                     //6-9
    // QZSS
    { 0,     CLIGHT/FREQ1,      CLIGHT/FREQ2,     0,                0,                 CLIGHT/FREQ5,   //0-5
             CLIGHT/FREQ6,      0,                0,                0            },                    //6-9
    // IRNSS
    { 0,     0,                 0,                0,                0,                 CLIGHT/FREQ5,   //0-5
             0,                 0,                0,                CLIGHT/FREQ9 },                    //6-9
    // SBAS
    { 0,     CLIGHT/FREQ1,      CLIGHT/FREQ2,     0,                0,                 CLIGHT/FREQ5,   //0-5
             0,                 0,                0,                0            }                     //6-9
};
/* default gnss frequency priority -------------------------------------------------------------------
* the 2nd number of "MAX_NF" values are only set for BDS-2
* ( GPS + GLO + GAL + BDS + QZS + IRN + LEO ) ----------------------------------------------------- */
const int GNSS_FREQ_PRI_DEF[MAX_GNSS][MAX_NF*2] = {
    // GPS
    { 1, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* don't change GPS frequency */
    // GLONASS
    { 1, 2, 3, 4, 6, 0, 0, 0, 0, 0, 0, 0 },
    // Galileo
    { 1, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0 },
    // BDS
    { 2, 6, 7, 5, 8, 1, 2, 6, 7, 5, 8, 1 }, /* 0:5 for BDS-3 and 6:11 for BDS-2 */
    // QZSS
    { 1, 2, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0 },
    // IRNSS
    { 5, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* don't change IRNSS frequency */
    // SBAS
    { 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  /* don't change SBAS frequency */
};
/* default used GNSS frequency -----------------------------------------------------------------------
* the 2nd number of "MAX_NF" values are only set for BDS-2
* ( GPS + GLO + GAL + BDS + QZS + IRN + LEO ) ----------------------------------------------------- */
const int GNSS_USED_FREQ_DEF[MAX_GNSS][NFREQ*2] = {
    // GPS
    { 1, 2, 5, 0, 0, 0 }, /* don't change GPS frequency */
    // GLONASS
    { 1, 2, 3, 0, 0, 0 },
    // Galileo
    { 1, 5, 6, 0, 0, 0 },
    // BDS
    { 2, 6, 7, 2, 6, 7 }, /* 0:2 for BDS-3 and 3:5 for BDS-2 */
    // QZSS
    { 1, 2, 5, 0, 0, 0 },
    // IRNSS
    { 5, 9, 0, 0, 0, 0 }, /* don't change IRNSS frequency */
    // SBAS
    { 1, 5, 0, 0, 0, 0 }  /* don't change SBAS frequency */
};

typedef void fatalfunc_t(const char *); /* fatal callback function type */

/* shared structures ------------------------------------------------------------------------------ */
/* SNR mask structure ------------------------------------------------------------- */
typedef struct {
    int ena[2];                         /* enable flag {rover,base} */
    double mask[NFREQ][9];              /* mask (dBHz) at 5,10,...85 deg */
} snrmask_t;
/* structure: signal index structure -------------------------------------------------------------- */
typedef struct {                        /* signal index type */
    int n;                              /* number of index */
    string obsCode[MAXOBSTYPE];         /* obs code (C??,L??,D??,S??) */
    int frq[MAXOBSTYPE];                /* signal frequency (1:L1,2:L2,...) */
    int usedFrq[MAXOBSTYPE];            /* used signal frequency (1:F1, 2:F2, 3:F3, else:not used) */
    int BDS2Frq[MAXOBSTYPE];            /* BDS-2 used signal frequency (1:F1, 2:F2, 3:F3, else:not used) */
    int avaiFrq[MAX_NF*2];              /* available frequency */
    int pos[MAXOBSTYPE];                /* signal index in obs data (-1:no) */
    int pri[MAXOBSTYPE];                /* signal priority (15-0) */
    int type[MAXOBSTYPE];               /* type (0:C,1:L,2:D,3:S) */
    int code[MAXOBSTYPE];               /* obs code (CODE_L??) */
    double shift[MAXOBSTYPE];           /* phase shift (cycle) */
} sigind_t;

/* All classes ------------------------------------------------------------------------------------ */
/* basic classes ---------------------------------------------------------------------------------- */
/* mathematic class --------------------------------------------------------------- */
/* 3D vector class ---------------------------------------------------------------- */
class trivector {
/* Constructors */
public:
    trivector() { value[0]=value[1]=value[2]=0; norm=0; }
    trivector(double v1,double v2,double v3) { 
        value[0]=v1; value[1]=v2; value[2]=v3; 
        norm=sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]); 
    }
    ~trivector() {  };
/* Operator Overloading */
public:
    void operator=(const vector<double> &VecSrc) {
        if (VecSrc.size()==3) {
            value[0]=VecSrc[0]; value[1]=VecSrc[1]; value[2]=VecSrc[2];
            norm=sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]); 
        }
    }
    double& operator[](int i) {
        if (i<0||i>2) return value[0];
        return value[i];
    }
    const double& operator[](int i) const {
        if (i<0||i>2) return value[0];
        return value[i];
    }
    double* operator+(int i) {
        if (i<0||i>2) return value;
        return value+i;
    }
    const double* operator+(int i) const {
        if (i<0||i>2) return value;
        return value+i;
    }
/* Implementation */
    double this_norm() { return norm; }
    double get_norm() { return norm=sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]); }
protected:
    double norm;
/* Components */
public:
    double value[3];
};
/* 9D vector class ---------------------------------------------------------------- */
class ninevector {
    /* Constructors */
public:
    ninevector() { value[0]=value[1]=value[2]=0;
                   value[3]=value[4]=value[5]=0;
                   value[6]=value[7]=value[8]=0; }
    ninevector(double v1,double v2,double v3) { value[0]=v1; value[4]=v2; value[8]=v3;
                                                value[1]=value[2]=value[3]=0; 
                                                value[5]=value[6]=value[7]=0; }
    ~ninevector() {  };
    /* Operator Overloading */
public:
    void operator=(const vector<double> &VecSrc) {
        if (VecSrc.size()==9) {
            value[0]=VecSrc[0]; value[1]=VecSrc[1]; value[2]=VecSrc[2];
            value[3]=VecSrc[3]; value[4]=VecSrc[4]; value[5]=VecSrc[5];
            value[6]=VecSrc[6]; value[7]=VecSrc[7]; value[8]=VecSrc[8];
        }
    }
    double& operator[](int i) {
        if (i<0||i>8) return value[0];
        return value[i];
    }
    const double& operator[](int i) const {
        if (i<0||i>8) return value[0];
        return value[i];
    }
    double* operator+(int i) {
        if (i<0||i>8) return value;
        return value+i;
    }
    const double* operator+(int i) const {
        if (i<0||i>8) return value;
        return value+i;
    }
    /* Components */
public:
    double value[9];
};

/* time strcuture class ----------------------------------------------------------- */
class gtime_c;

/* adjustment function class ------------------------------------------------------ */
class adjfunc_c;                        /* parent class of adjustment functions */
class lsadj_c;                          /* least square adjustment function */
class kalmanfilter_c;                   /* kalman filter adjustment function */
class helmert_c;                        /* helmert components covariance estimate
                                         * in kalman filter */

/* configure option class --------------------------------------------------------- */
/* GNSS options */
class gnss_pstopt_c;                    /* post-processing option class */
class gnss_rtkopt_c;                    /* rtk-processing option class */
class gnss_prcopt_c;                    /* GNSS processing options class */
class gnss_solopt_c;                    /* solution options class */
class gnss_filopt_t;                    /* file options class */
class gnss_rnxopt_c;                    /* RINEX options class */
/* all option */
class all_option_c;                     /* all options */

/* stream class ------------------------------------------------------------------- */
class stream_c;                         /* parent class of stream */
class file_c;                           /* stream :: file control class */
class tcpsvr_c;                         /* stream :: tcp server class */
class tcpcli_c;                         /* stream :: tcp cilent class */
class serial_c;                         /* stream :: tcpsvr_c :: serial control class */
class ntrip_c;                          /* stream :: tcpcli_c :: ntrip control class */
class ntripc_c;                         /* stream :: tcpsvr_c :: ntrip caster control class */
class udp_c;                            /* stream :: udp class */
class ftp_c;                            /* stream :: ftp download control class */
class membuf_c;                         /* stream :: memory buffer class */

/* decode data class -------------------------------------------------------------- */
class decode_data_c;                    /* parent class of decode data */
class rtcm_c;                           /* RTCM control struct class */
class rtcm2_c;                          /* RTCM 2 control struct class */
class rtcm3_c;                          /* RTCM 3 control struct class */
class raw_c;                            /* receiver raw data control class */

/* input file class --------------------------------------------------------------- */
class in_gnss_rnx_c;                    /* parent class of read rinxe file */
class in_gnss_rnxO_c;                   /* read rinex file (.O) */
class in_gnss_rnxN_c;                   /* read rinex file (.N,.P) */
class in_gnss_rnxC_c;                   /* read rinex Precise Clock file (.CLK) */
class in_gnss_erp_c;                    /* read earth rotation parameters file (.ERP) */
class in_gnss_blq_c;                    /* read ocean-loading tide file (.BLQ) */
class in_gnss_atx_c;                    /* read antenna information file (.ATX) */
class in_gnss_eph_c;                    /* read precise ephemeris file */
class in_gnss_ionex_c;                  /* read ionex tec grid file */
class in_gnss_bias_c;                   /* read GNSS code bias file */
class in_gnss_snx_c;                    /* read GNSS SNX file */
class in_gnss_ztd_c;                    /* read troposphere ZTD file (debug) */

/* GNSS class ------------------------------------------------------------------------------------- */
/* GNSS processing class ---------------------------------------------------------- */
class gnss_pro_c;                       /* GNSS process class */
class gnss_sol_c;                       /* GNSS solution class */
class gnss_ssat_c;                      /* GNSS satellite status class */
class gnss_rtksvr_c;                    /* GNSS RTK server class */
class gnss_postsvr_c;                   /* GNSS post-processing server class */

/* GNSS data class ---------------------------------------------------------------- */
class gnss_freq_c;                      /* frequency priority and adoption */
class gnss_obsd_c;                      /* one epoch observation data */
class gnss_sta_c;                       /* station informtations */
class gnss_obs_c;                       /* station's information and observation chain */
class gnss_eph_c;                       /* GPS/QZS/GAL broadcast ephemeris class */
class gnss_geph_c;                      /* GLONASS broadcast ephemeris class */
class gnss_seph_c;                      /* SBAS ephemeris class */
class gnss_peph_c;                      /* precise ephemeris class */
class gnss_pclk_c;                      /* precise clock class */
class gnss_alm_c;                       /* almanac class */
class gnss_tec_c;                       /* TEC grid class */
class gnss_bias_c;                      /* GNSS observation bias class */
class gnss_ztd_c;                       /* ZTD data class */
class gnss_fcbd_c;                      /* satellite fcb data class */
class gnss_erpd_c;                      /* earth rotation parameter data class */
class gnss_erp_c;                       /* earth rotation parameter class */
class gnss_pcv_c;                       /* antenna parameter class */
class gnss_sbsfcorr_c;                  /* SBAS fast correction class */
class gnss_sbslcorr_c;                  /* SBAS long term satellite error correction class */
class gnss_sbssatp_c;                   /* SBAS satellite correction class */
class gnss_sbssat_c;                    /* SBAS satellite corrections class */
class gnss_sbsigp_c;                    /* SBAS ionospheric correction class */
class gnss_sbsion_c;                    /* SBAS ionospheric corrections class */
class gnss_dgps_c;                      /* DGPS/GNSS correction class */
class gnss_ssr_c;                       /* SSR correction class */
class gnss_lexeph_c;                    /* QZSS LEX ephemeris class */
class gnss_lexion_c;                    /* QZSS LEX ionosphere correction class */
class gnss_stec_c;                      /* stec data class */
class gnss_trop_c;                      /* trop data class */
class gnss_pppcorr_c;                   /* ppp corrections class */
class gnss_nav_c;                       /* class of navigation data */
class gnss_sbsmsg_c;                    /* SBAS message class */
class gnss_ssbs_c;                      /* SBAS messages class */
class gnss_lexmsg_c;                    /* QZSS LEX message class */
class gnss_lex_c;                       /* QZSS LEX messages class */

/* antenna function class --------------------------------------------------------- */
class gnss_satant_c;                    /* satellite antenna phase center offest */
class gnss_recant_c;                    /* receiver antenna phase center offest */

/* ionosphere correction function class ------------------------------------------- */
class gnss_ioncorr_c;                   /* parent class of ionosphere (GPS L1) */
class gnss_LCion_c;                     /* ionosphere free combination */
class gnss_broadion_c;                  /* broadcast ionosphere correction */
class gnss_qzssion_c;                   /* qzss broadcast ionosphere correction */
class gnss_sbasion_c;                   /* sbas ionosphere correction */
class gnss_ionexion_c;                  /* ionex ionosphere correction */
class gnss_lexioncor_c;                 /* lex ionosphere correction */
class gnss_estion_c;                    /* constrained ionosphere model correction */

/* troposphere correction function class ------------------------------------------ */
class gnss_tromod_c;                    /* troposphere model */
class gnss_trocorr_c;                   /* parent class of troposphere */
class gnss_saastro_c;                   /* saastamoinen model troposphere correction */
class gnss_sbastro_c;                   /* sbas model troposphere correction */
class gnss_esttro_c;                    /* estimated model troposphere correction */
class gnss_ztdtro_c;                    /* ZTD data model troposphere correction */

/* earth tide and ocean load function class --------------------------------------- */
class gnss_tidecorr_c;                  /* all tidal function */

/* satellite function class ------------------------------------------------------- */
class gnss_satellite_c;                 /* parent class of satellite functions */
class gnss_broadcast_c;                 /* broadcast ephemeris */
class gnss_broadsbas_c;                 /* broadcast ephemeris & sbas correction */
class gnss_broadssrapc_c;               /* broadcast ephemeris & ssr_apc correction */
class gnss_broadssrcom_c;               /* broadcast ephemeris & ssr_com correction */
class gnss_preciseph_c;                 /* precise ephemeris */
class gnss_qzsslex_c;                   /* qzss lex ephemeris */

/* parameter function class ------------------------------------------------------- */
class gnss_parameter_c;                 /* parameter function */

/* integar ambiguity strategy class ----------------------------------------------- */
class gnss_intamb_c;                    /* parent class of integar solution */
class gnss_lambda_c;                    /* lambda/mlambda integer LSQ estimation */

/* GNSS position class ------------------------------------------------------------ */
class gnss_single_c;                    /* single point position class */
class gnss_ppp_c;                       /* precise point position class */
class gnss_relative_c;                  /* relative position class */

#endif