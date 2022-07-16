/* Class of Reading Positioning Configure file */
#include "ConfigFile/config.h"
#include "BaseFunction/basefunction.h"

/* constant --------------------------------------------------------------------------------------- */
/* configuration string ----------------------------------------------------------- */
const int sys[]={ SYS_GPS, SYS_GLO, SYS_GAL, SYS_BDS, SYS_QZS, SYS_IRN, SYS_SBS, 0 };
const int iBDSGEN[]={ BDS_2ND, BDS_3RD, 0 };
const int iBDSOBT[]={ BDS_ORB_MEO, BDS_ORB_IGSO, BDS_ORB_GEO, 0 };
const int CSmode[]={ CSMODE_CP, CSMODE_GF, CSMODE_MW, 0 };
const char *s_posmode[]={ "spp","ppp-kinematic","ppp-static","ppp-fixed",
"dgps","rtk-kinematic","rtk-static","rtk-fixed","moving-base","" };
const char *s_pppmod[]={ "un-differenced", "single-differenced" };
const char *s_freq[]={ "single","double","triple","" };
const char *s_soltype[]={ "forward","backward","combined" };
const char *s_sppion[]={ "off","broadcast","sbas","ion-free","ionex-tec","qzs-brdc","qzs-lex","" };
const char *s_ionopt[]={ "off","broadcast","sbas","ion-free","constrain",
                        "ionex tec","qzs","lex","vtec_sf","vtec_ef","gtec","" };
const char *s_spptro[]={ "off","saastamoinen","sbas","ztd-file","" };
const char *s_troopt[]={ "off","saastamoinen","sbas","est-ztd","est-ztd+grad","ztd-file","" };
const char *s_ISBmode[]={ "white noise", "random walk" };
const char *s_sateph[]={ "broadcast","precise","broadcast+sbas","broadcast+ssr apc",
"broadcast+ssr com","qzss lex","" };
const char *s_navsys[]={ "GPS","GLONASS","Galileo","BDS","QZSS","IRNSS","SBAS","" };
const char *s_BDSGEN[]={ "BDS-2","BDS-3","" };
const char *s_BDSOBT[]={ "MEO","IGSO","GEO","" };
const char *s_armode[]={ "off","continuous","instantaneous","fix and hold","" };
const char *s_csmode[]={ "code-phase","geo-free","Melbourne-Wubbena","" };
const char *s_glomode[]={ "off","on","auto calib","external calib","" };
const char *s_adjmode[]={ "Least Square","Kalman Filter","Helmert VCE based Kalman Filter" };
const char *s_adjrobust[]={ "none","IGG3-2","IGG3-3" };
const char *s1[]={ "WGS84","CGCS2000" },*s2[]={ "ellipsoidal","geodetic" };
const char *s3[]={ "GPST","UTC ","BDST " };

/* static option value ---------------------------------------------------------------------------- */
/* system options buffer -----------------------------------------------------*/
static gnss_prcopt_c prcopt_=gnss_prcopt_c();
static gnss_solopt_c solopt1_=gnss_solopt_c(),solopt2_=gnss_solopt_c();
static gnss_filopt_t filopt_=gnss_filopt_t();
static int antpostype_[2];
static double elmask_,min_cw_el_,minElArFix_,minElArHold_;
static double antpos_[2][3];
static string exsats_;
static string snrmask_[NFREQ];
static string freq_glo_,freq_gal_,freq_bds3_,freq_bds2_;

/* constant */
/* common options table --------------------------------------------------------------------------- */
#define SWTOPT          "0:off,1:on"
#define ADJROBUSTOPT    "0:none,1:IGG3-2,2:IGG3-3"
#define ROBUSTVAROPT    "0:var,1:covar,2:vector"
#define SOLOPT          "0:llh,1:xyz,2:enu,3:nmea"
#define POSOPT          "0:llh,1:xyz,2:single,3:posfile,4:snxfile,5:rinexhead,6:rtcm,7:raw"
#define POSOPT2         "0:llh,1:xyz"

/* GNSS options table ----------------------------------------------------------------------------- */
#define MODOPT  "0:spp,1:ppp-kin,2:ppp-sta,3:ppp-fix,4:dgps,5:rtk-kin,6:rtk-sta,7:rtk-fix,8:movingbase"
#define PPPMOD  "0:un-diff,1:single-diff"
#define ADJOPT  "0:lsa,1:kalman,2:helmert"
#define FRQOPT  "1:single,2:double,3:triple"
#define TYPOPT  "0:forward,1:backward,2:combined"
#define SPPION  "0:off,1:brdc,2:sbas,3:ion-free,5:ionex-tec,6:qzs-brdc,7:qzs-lex"
#define IONOPT  "0:off,1:brdc,2:sbas,3:ion-free,4:constrain,5:ionex-tec,6:qzs-brdc,7:qzs-lex,8:stec"
#define SPPTRO  "0:off,1:saas,2:sbas,5:ztd-file"
#define TRPOPT  "0:off,1:saas,2:sbas,3:est-ztd,4:est-ztdgrad,5:ztd-file"
#define EPHOPT  "0:brdc,1:precise,2:brdc+sbas,3:brdc+ssrapc,4:brdc+ssrcom"
#define BDSGEN  "1:BDS-2 + 2:BDS-3"
#define BDSOBT  "1:MEO + 2:IGSO + 4:GEO"
#define NAVOPT  "1:gps + 2:glo + 4:gal + 8:bds + 16:qzs + 32:irn + 64:sbas"
#define SLPOPT  "0:obs+1:CP+2:geo-free+4:Melbourne-Wubbena"
#define GAROPT  "0:off,1:on,2:auto"
#define TSYOPT  "0:gpst,1:utc,2:bdst"
#define TFTOPT  "0:tow,1:hms"
#define DFTOPT  "0:deg,1:dms"
#define DTMOPT  "0:WGS84,1:CGCS2000"
#define HGTOPT  "0:ellipsoidal,1:geodetic"
#define GEOOPT  "0:internal,1:egm96,2:egm08_2.5,3:egm08_1,4:gsi2000"
#define STAOPT  "0:all,1:single"
#define STSOPT  "0:off,1:state,2:residual"
#define ARMOPT  "0:off,1:continuous,2:instantaneous,3:fix-and-hold,4:LC_WN"
#define REFOPT  "0:base,1:rover"
#define TIDEOPT "0:off+1:solid+2:otl+4:pole"
#define PHWOPT  "0:off,1:on,2:precise"
#define LOGOPT  "0:off + 1:nsat + 2:pdop + 4:snr + 8:multi-path + 16:exc-sat + 32:helmert + 4096:xpar + 8192:dx"
#define ISBMOD  "0:whiteNoise,1:randomWalk"
EXPORT opt_t sysopts[]={
    { "gnss-posmode",       3,  (void *)&prcopt_.mode,           MODOPT       },
    { "gnss-pppmode",       3,  (void *)&prcopt_.pppmode,        PPPMOD       },
    { "gnss-frequency",     3,  (void *)&prcopt_.nf,             FRQOPT       },
    { "gnss-freqPriority",  3,  (void *)&prcopt_.freqPriority,   SWTOPT       },
    { "gnss-soltype",       3,  (void *)&prcopt_.soltype,        TYPOPT       },
    { "gnss-elmask",        1,  (void *)&elmask_,                "deg"        },
    { "gnss-con_weight_el", 1,  (void *)&min_cw_el_,             "deg"        },
    { "gnss-snrmask_r",     3,  (void *)&prcopt_.snrmask.ena[0], SWTOPT       },
    { "gnss-snrmask_b",     3,  (void *)&prcopt_.snrmask.ena[1], SWTOPT       },
    { "gnss-snrmask_F1",    2,  (void *)&snrmask_[0],            ""           },
    { "gnss-snrmask_F2",    2,  (void *)&snrmask_[1],            ""           },
    { "gnss-snrmask_F3",    2,  (void *)&snrmask_[2],            ""           },
    { "gnss-dynamics",      3,  (void *)&prcopt_.dynamics,       SWTOPT       },
    { "gnss-tidecorr",      0,  (void *)&prcopt_.tidecorr,       TIDEOPT      },
    { "gnss-sppiono",       3,  (void *)&prcopt_.sppiono,        SPPION       },
    { "gnss-ionoopt",       3,  (void *)&prcopt_.ionoopt,        IONOPT       },
    { "gnss-spptrop",       3,  (void *)&prcopt_.spptrop,        SPPTRO       },
    { "gnss-tropopt",       3,  (void *)&prcopt_.tropopt,        TRPOPT       },
    { "gnss-sateph",        3,  (void *)&prcopt_.sateph,         EPHOPT       },
    { "gnss-posopt1",       3,  (void *)&prcopt_.posopt[0],      SWTOPT       },
    { "gnss-posopt2",       3,  (void *)&prcopt_.posopt[1],      SWTOPT       },
    { "gnss-posopt3",       3,  (void *)&prcopt_.posopt[2],      PHWOPT       },
    { "gnss-exclsats",      2,  (void *)&exsats_,                "prn ..."    },
    { "gnss-bds_gen",       0,  (void *)&prcopt_.bdsGen,         BDSGEN       },
    { "gnss-bds_orbit",     0,  (void *)&prcopt_.bdsOrbit,       BDSOBT       },
    { "gnss-navsys",        0,  (void *)&prcopt_.navsys,         NAVOPT       },
    { "gnss-freq_glo",      2,  (void *)&freq_glo_,              "1,2,3,..."  },
    { "gnss-freq_gal",      2,  (void *)&freq_gal_,              "1,5,6,..."  },
    { "gnss-freq_bds3",     2,  (void *)&freq_bds3_,             "1,2,6,..."  },
    { "gnss-freq_bds2",     2,  (void *)&freq_bds2_,             "1,2,7,..."  },
    { "gnss-logmsg",        0,  (void *)&prcopt_.logmsg,         LOGOPT       },

    { "proc-armode",        3,  (void *)&prcopt_.modear,         ARMOPT       },
    { "proc-order",         0,  (void *)&prcopt_.order,          ""           },
    { "proc-slipmode",      0,  (void *)&prcopt_.slipmode,       SLPOPT       },
    { "proc-slipstd",       1,  (void *)&prcopt_.slip_std,       "m"          },
    { "proc-cs_gf",         1,  (void *)&prcopt_.cs_gf,          "m"          },
    { "proc-cs_mw",         1,  (void *)&prcopt_.cs_mw,          "m"          },
    { "proc-sampling",      1,  (void *)&prcopt_.sampling,       "s"          },
    { "proc-minConse",      0,  (void *)&prcopt_.minConsecutive, ""           },
    { "proc-restime",       1,  (void *)&prcopt_.restime,        "s"          },
    { "proc-gloarmode",     3,  (void *)&prcopt_.glomodear,      GAROPT       },
    { "proc-bdsarmode",     3,  (void *)&prcopt_.bdsmodear,      SWTOPT       },
    { "proc-arthres",       1,  (void *)&prcopt_.thresar[0],     ""           },
    { "proc-arthres1",      1,  (void *)&prcopt_.thresar[1],     ""           },
    { "proc-arthres2",      1,  (void *)&prcopt_.thresar[2],     ""           },
    { "proc-arthres3",      1,  (void *)&prcopt_.thresar[3],     ""           },
    { "proc-arthres4",      1,  (void *)&prcopt_.thresar[4],     ""           },
    { "proc-iniar",         0,  (void *)&prcopt_.iniamb,         "n"          },
    { "proc-minelarfix",    1,  (void *)&minElArFix_,            "deg"        },
    { "proc-minelarhold",   1,  (void *)&minElArHold_,           "deg"        },
    { "proc-maxariter",     0,  (void *)&prcopt_.maxariter,      ""           },
    { "proc-maxage",        1,  (void *)&prcopt_.maxtdiff,       "s"          },
    { "proc-syncsol",       3,  (void *)&prcopt_.syncsol,        SWTOPT       },
    { "proc-maxres",        1,  (void *)&prcopt_.maxres,         ""           },
    { "proc-adjust",        3,  (void *)&prcopt_.adjustfunc,     ADJOPT       },
    { "proc-adjust-robust", 3,  (void *)&prcopt_.adjRobust,      ADJROBUSTOPT },
    { "proc-robust-covar",  3,  (void *)&prcopt_.robustRtype,    ROBUSTVAROPT },
    { "proc-niter",         0,  (void *)&prcopt_.niter,          ""           },
    { "proc-baselen",       1,  (void *)&prcopt_.baseline[0],    "m"          },
    { "proc-basesig",       1,  (void *)&prcopt_.baseline[1],    "m"          },
    { "proc-baseatt",       3,  (void *)&prcopt_.baseAtt,        SWTOPT       },

    { "out1-solformat",     3,  (void *)&solopt1_.posf,          SOLOPT       },
    { "out1-refstation",    3,  (void *)&solopt1_.refsta,        REFOPT       },
    { "out1-outhead",       3,  (void *)&solopt1_.outhead,       SWTOPT       },
    { "out1-outopt",        3,  (void *)&solopt1_.outopt,        SWTOPT       },
    { "out1-timesys",       3,  (void *)&solopt1_.times,         TSYOPT       },
    { "out1-timeform",      3,  (void *)&solopt1_.timef,         TFTOPT       },
    { "out1-timendec",      0,  (void *)&solopt1_.timeu,         ""           },
    { "out1-degform",       3,  (void *)&solopt1_.degf,          DFTOPT       },
    { "out1-datum",         3,  (void *)&solopt1_.datum,         DTMOPT       },
    { "out1-height",        3,  (void *)&solopt1_.height,        HGTOPT       },
    { "out1-geoid",         3,  (void *)&solopt1_.geoid,         GEOOPT       },
    { "out1-solstatic",     3,  (void *)&solopt1_.solstatic,     STAOPT       },
    { "out1-nmeaintv1",     1,  (void *)&solopt1_.nmeaintv[0],   "s"          },
    { "out1-nmeaintv2",     1,  (void *)&solopt1_.nmeaintv[1],   "s"          },
    { "out1-outstat",       3,  (void *)&solopt1_.sstat,         STSOPT       },

    { "out2-solformat",     3,  (void *)&solopt2_.posf,          SOLOPT       },
    { "out2-refstation",    3,  (void *)&solopt2_.refsta,        REFOPT       },
    { "out2-outhead",       3,  (void *)&solopt2_.outhead,       SWTOPT       },
    { "out2-outopt",        3,  (void *)&solopt2_.outopt,        SWTOPT       },
    { "out2-timesys",       3,  (void *)&solopt2_.times,         TSYOPT       },
    { "out2-timeform",      3,  (void *)&solopt2_.timef,         TFTOPT       },
    { "out2-timendec",      0,  (void *)&solopt2_.timeu,         ""           },
    { "out2-degform",       3,  (void *)&solopt2_.degf,          DFTOPT       },
    { "out2-datum",         3,  (void *)&solopt2_.datum,         DTMOPT       },
    { "out2-height",        3,  (void *)&solopt2_.height,        HGTOPT       },
    { "out2-geoid",         3,  (void *)&solopt2_.geoid,         GEOOPT       },
    { "out2-solstatic",     3,  (void *)&solopt2_.solstatic,     STAOPT       },
    { "out2-nmeaintv1",     1,  (void *)&solopt2_.nmeaintv[0],   "s"          },
    { "out2-nmeaintv2",     1,  (void *)&solopt2_.nmeaintv[1],   "s"          },
    { "out2-outstat",       3,  (void *)&solopt2_.sstat,         STSOPT       },

    { "stats-pppfactor",    1,  (void *)&prcopt_.pppfactor,      ""           },
    { "stats-errcode",      1,  (void *)&prcopt_.err[0],         "m"          },
    { "stats-errcodeel",    1,  (void *)&prcopt_.err[1],         "m"          },
    { "stats-errphase",     1,  (void *)&prcopt_.err[2],         "m"          },
    { "stats-errphaseel",   1,  (void *)&prcopt_.err[3],         "m"          },
    { "stats-errphasebl",   1,  (void *)&prcopt_.err[4],         "m/10km"     },
    { "stats-errdoppler",   1,  (void *)&prcopt_.err[5],         "Hz"         },
    { "stats-ISBmode",      3,  (void *)&prcopt_.ISBmode,        ISBMOD       },
    { "stats-stdambs",      1,  (void *)&prcopt_.std[0],         "m"          },
    { "stats-stdiono",      1,  (void *)&prcopt_.std[1],         "m"          },
    { "stats-stdtrop",      1,  (void *)&prcopt_.std[2],         "m"          },
    { "stats-ratambs",      1,  (void *)&prcopt_.stdrate[0],     "m^2/s"      },
    { "stats-rationo",      1,  (void *)&prcopt_.stdrate[1],     "m^2/s"      },
    { "stats-rattrop",      1,  (void *)&prcopt_.stdrate[2],     "m^2/s"      },
    { "stats-ratrISB",      1,  (void *)&prcopt_.stdrate[3],     "m^2/s"      },
    { "stats-clkstab",      1,  (void *)&prcopt_.sclkstab,       "s/s"        },

    { "rover-postype",      3,  (void *)&antpostype_[0],         POSOPT       },
    { "rover-name",         2,  (void *)&prcopt_.name[0],        ""           },
    { "rover-pos1",         1,  (void *)&antpos_[0][0],          "deg|m"      },
    { "rover-pos2",         1,  (void *)&antpos_[0][1],          "deg|m"      },
    { "rover-pos3",         1,  (void *)&antpos_[0][2],          "m|m"        },
    { "rover-anttype",      2,  (void *)&prcopt_.anttype[0],     ""           },
    { "rover-antdele",      1,  (void *)&prcopt_.antdel[0][0],   "m"          },
    { "rover-antdeln",      1,  (void *)&prcopt_.antdel[0][1],   "m"          },
    { "rover-antdelu",      1,  (void *)&prcopt_.antdel[0][2],   "m"          },
    { "base-postype",       3,  (void *)&antpostype_[1],         POSOPT       },
    { "base-name",          2,  (void *)&prcopt_.name[1],        ""           },
    { "base-pos1",          1,  (void *)&antpos_[1][0],          "deg|m"      },
    { "base-pos2",          1,  (void *)&antpos_[1][1],          "deg|m"      },
    { "base-pos3",          1,  (void *)&antpos_[1][2],          "m|m"        },
    { "base-anttype",       2,  (void *)&prcopt_.anttype[1],     ""           },
    { "base-antdele",       1,  (void *)&prcopt_.antdel[1][0],   "m"          },
    { "base-antdeln",       1,  (void *)&prcopt_.antdel[1][1],   "m"          },
    { "base-antdelu",       1,  (void *)&prcopt_.antdel[1][2],   "m"          },
    { "base-maxaveep",      0,  (void *)&prcopt_.maxaveep,       ""           },
    { "base-initrst",       3,  (void *)&prcopt_.initrst,        SWTOPT       },

    { "misc-timeinterp",    3,  (void *)&prcopt_.intpref,        SWTOPT       },
    { "misc-sbasatsel",     0,  (void *)&prcopt_.sbassatsel,     "0:all"      },
    { "misc-rnxopt1",       2,  (void *)&prcopt_.rnxopt[0],      ""           },
    { "misc-rnxopt2",       2,  (void *)&prcopt_.rnxopt[1],      ""           },
    { "misc-pppopt",        2,  (void *)&prcopt_.pppopt,         ""           },

    { "file-satantfile",    2,  (void *)&filopt_.satantp,        ""           },
    { "file-rcvantfile",    2,  (void *)&filopt_.rcvantp,        ""           },
    { "file-staposfile",    2,  (void *)&filopt_.stapos,         ""           },
    { "file-geoidfile",     2,  (void *)&filopt_.geoid,          ""           },
    { "file-ionofile",      2,  (void *)&filopt_.iono,           ""           },
    { "file-ztdfile",       2,  (void *)&filopt_.ztd,            ""           },
    { "file-dcbfile",       2,  (void *)&filopt_.dcb,            ""           },
    { "file-erpfile",       2,  (void *)&filopt_.erp,            ""           },
    { "file-blqfile",       2,  (void *)&filopt_.blq,            ""           },
    { "file-snxfile",       2,  (void *)&filopt_.snx,            ""           },
    { "file-tempdir",       2,  (void *)&filopt_.tempdir,        ""           },
    { "file-geexefile",     2,  (void *)&filopt_.geexe,          ""           },
    { "file-solstatfile",   2,  (void *)&filopt_.solstat,        ""           },
    { "file-testfile",      2,  (void *)&filopt_.test,           ""           },

    { "",0,NULL,"" } /* terminator */
};

/* rtkpro options table --------------------------------------------------------------------------- */
static gnss_rtkopt_c rtkopt_;
#define TIMOPT  "0:gpst,1:utc,2:jst,3:tow"
#define CONOPT  "0:dms,1:deg,2:xyz,3:enu,4:pyl"
#define FLGOPT  "0:off,1:std+2:age/ratio/ns"
#define ISTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,7:ntripcli,8:ftp,9:http"
#define OSTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripsvr"
#define FMTOPT  "0:rtcm2,1:rtcm3,2:oem4,3:oem3,4:ubx,5:ss2,6:hemis,7:skytraq,8:gw10,9:javad,10:nvs,11:binex,12:rt17,15:sp3"
#define NMEOPT  "0:off,1:latlon,2:single"
#define MSGOPT  "0:all,1:rover,2:base,3:corr"
static opt_t rtkopts[]={
    { "inpstr1-type",    3,  (void *)&rtkopt_.strtype[0],         ISTOPT },
    { "inpstr2-type",    3,  (void *)&rtkopt_.strtype[1],         ISTOPT },
    { "inpstr3-type",    3,  (void *)&rtkopt_.strtype[2],         ISTOPT },
    { "inpstr1-path",    2,  (void *)&rtkopt_.strpath[0],         ""     },
    { "inpstr2-path",    2,  (void *)&rtkopt_.strpath[1],         ""     },
    { "inpstr3-path",    2,  (void *)&rtkopt_.strpath[2],         ""     },
    { "inpstr1-format",  3,  (void *)&rtkopt_.strfmt[0],          FMTOPT },
    { "inpstr2-format",  3,  (void *)&rtkopt_.strfmt[1],          FMTOPT },
    { "inpstr3-format",  3,  (void *)&rtkopt_.strfmt[2],          FMTOPT },
    { "inpstr2-nmeareq", 3,  (void *)&rtkopt_.nmeareq,            NMEOPT },
    { "inpstr2-nmealat", 1,  (void *)&rtkopt_.nmeapos[0],         "deg"  },
    { "inpstr2-nmealon", 1,  (void *)&rtkopt_.nmeapos[1],         "deg"  },
    { "inpstr2-nmeahig", 1,  (void *)&rtkopt_.nmeapos[2],         "m"    },
    { "outstr1-type",    3,  (void *)&rtkopt_.strtype[3],         OSTOPT },
    { "outstr2-type",    3,  (void *)&rtkopt_.strtype[4],         OSTOPT },
    { "outstr1-path",    2,  (void *)&rtkopt_.strpath[3],         ""     },
    { "outstr2-path",    2,  (void *)&rtkopt_.strpath[4],         ""     },
    { "logstr1-type",    3,  (void *)&rtkopt_.strtype[5],         OSTOPT },
    { "logstr2-type",    3,  (void *)&rtkopt_.strtype[6],         OSTOPT },
    { "logstr3-type",    3,  (void *)&rtkopt_.strtype[7],         OSTOPT },
    { "logstr1-path",    2,  (void *)&rtkopt_.strpath[5],         ""     },
    { "logstr2-path",    2,  (void *)&rtkopt_.strpath[6],         ""     },
    { "logstr3-path",    2,  (void *)&rtkopt_.strpath[7],         ""     },

    { "misc-svrcycle",   0,  (void *)&rtkopt_.svrcycle,           "ms"   },
    { "misc-timeout",    0,  (void *)&rtkopt_.timeout,            "ms"   },
    { "misc-reconnect",  0,  (void *)&rtkopt_.reconnect,          "ms"   },
    { "misc-nmeacycle",  0,  (void *)&rtkopt_.nmeacycle,          "ms"   },
    { "misc-buffsize",   0,  (void *)&rtkopt_.buffsize,           "bytes"},
    { "misc-navmsgsel",  3,  (void *)&rtkopt_.navmsgsel,          MSGOPT },
    { "misc-proxyaddr",  2,  (void *)&rtkopt_.proxyaddr,          ""     },
    { "misc-fswapmargin",0,  (void *)&rtkopt_.fswapmargin,        "s"    },

    { "",0,NULL,"" }
};

/* postpro options table -------------------------------------------------------------------------- */
static gnss_pstopt_c pstopt_;
string timeskip_="";
#define PREOPT  "1: only observed + 2: only predicted + 4: not combined"
static opt_t pstopts[]={
    { "post-predict",    0,  (void *)&pstopt_.predict,             PREOPT},
    { "post-timeinter",  1,  (void *)&pstopt_.time_inter,             "" },
    { "post-timestart",  2,  (void *)&pstopt_.time_start.sep,         "" },
    { "post-timeend",    2,  (void *)&pstopt_.time_end.sep,           "" },
    { "post-timeskip",   2,  (void *)&timeskip_,                      "" },
    { "post-roverobs",   2,  (void *)&pstopt_.rover_obs,              "" },
    { "post-baseobs",    2,  (void *)&pstopt_.base_obs,               "" },
    { "post-navigation", 2,  (void *)&pstopt_.nav,                    "" },
    { "post-preciseeph", 2,  (void *)&pstopt_.prseph,                 "" },
    { "post-satclock",   2,  (void *)&pstopt_.prsclk,                 "" },
    { "post-output1",    2,  (void *)&pstopt_.output[0],              "" },
    { "post-output2",    2,  (void *)&pstopt_.output[1],              "" },

    { "",0,NULL,"" }
};

/* GNSS configuration ----------------------------------------------------------------------------- */
/* post-processing option class ----------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
* defaults processing options --------------------------------------------------------------------- */
gnss_pstopt_c::gnss_pstopt_c() {
    predict=1;
    time_inter=30;
    time_start=time_end=gtime_c();
    time_skip[0].clear(); time_skip[1].clear();
    nskip_period=0;
}
gnss_pstopt_c::~gnss_pstopt_c() {
}
/* Implementation functions ------------------------------------------------------- */
/* set time parameters ------------------------------------------------------------ */
int gnss_pstopt_c::set_time_par() {
    /* start and end time */
    time_start.str2time(time_start.sep);
    time_end.str2time(time_end.sep);

    /* skipped time */
    nskip_period=0;
    time_skip[0].clear(); time_skip[1].clear();
    int p,pe,pt1e,pt2s,pt2e;
    gtime_c skipts,skipte;
    for (p=0; p<timeskip_.length(); p++) {
        if (timeskip_[p]==' ') continue; //if the first character is space
        pe=p;
        // search the end of one period
        if ((p=timeskip_.substr(p).find(','))==string::npos)
            p=timeskip_.length();
        else p+=pe;
        /* seperate start and end time of one period 
         * yyyy mm dd hh mm ss-yyyy mm dd hh mm ss
         * the start and end time string must be in the format of "yyyy mm dd hh mm ss" */
        if ( p-pe>38 && (pt2s=timeskip_.substr(pe,p-pe).find('-'))!=string::npos ) {
            // find the end of end time string
            for (pt2e=p-1; pt2e-pe>0; ) {
                if ( timeskip_[pt2e]<'0' || timeskip_[pt2e]>'9' ) pt2e--;
                else break;
            }
            pt2e-=pe;
            // if the length of skipped time string satisfies the min value
            if (pt2e>37) {
                // find the end of start time string
                for (pt1e=pt2s-1; pt1e>0; ) {
                    if ( timeskip_[pe+pt1e]<'0' || timeskip_[pe+pt1e]>'9' ) pt1e--;
                    else break;
                }
                if (pt1e>17) {
                    // set start time 
                    skipts.str2time(timeskip_.substr(pe,pt1e+1));
                    // find the start of end time string
                    for (pt2s=pt2s+1; pt2e-pt2s>0; ) {
                        if ( timeskip_[pe+pt2s]<'0' || timeskip_[pe+pt2s]>'9' ) pt2s++;
                        else break;
                    }
                    if (pt2e-pt2s>17) {
                        // set start time
                        skipte.str2time(timeskip_.substr(pe+pt2s,pt2e-pt2s+1));

                        /* if skipte>skipts, assign them to skipped time vector */
                        if (skipte.timediff(skipts)>0) {
                            time_skip[0].push_back(skipts);
                            time_skip[1].push_back(skipte);
                            nskip_period++;
                        }
                    }
                }
            }
        }
    }

    return nskip_period;
}

/* rtk-processing option class -----------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
* defaults processing options --------------------------------------------------------------------- */
gnss_rtkopt_c::gnss_rtkopt_c() {
    for (int i=0; i<7; i++) { strtype[i]=0; strpath[i]="\0"; }
    strfmt[0]=strfmt[1]=1; strfmt[2]=17;
    svrcycle=10; timeout=10000; reconnect=10000; nmeacycle=5000;
    fswapmargin=30; buffsize=32768; navmsgsel=0; nmeareq=0;
    nmeapos[0]=nmeapos[1]=nmeapos[2]=0;
    proxyaddr="\0";

    /* unset options */
    cmds[0]=cmds[1]=cmds[2]="\0";
    rropts[0]=rropts[1]=rropts[2]="\0";
    monitor=NULL;
}
gnss_rtkopt_c::~gnss_rtkopt_c() {
    if (monitor) { delete monitor; monitor=NULL; }
}

/* processing options class --------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
* defaults processing options --------------------------------------------------------------------- */
gnss_prcopt_c::gnss_prcopt_c() {
    int i;
    mode=PMODE_SINGLE;
    pppmode=PPP_UN_DIFF;
    soltype=0;
    nf=2; usedF=2; posF=2;
    freqPriority=0; navsys=SYS_GPS;
    min_elevation=15.0*D2R;
    con_weight_el=90.0*D2R;
    sateph=0; logmsg=0; bdsGen=BDS_3RD; bdsOrbit=BDS_ORB_MEO;
    modear=1; order=0; restime=3;
    sampling = 0.0;
    minConsecutive=3;
    glomodear=1; bdsmodear=1;
    iniamb=300; maxariter=1;
    sppiono=1; spptrop=1;
    ionoopt=0; tropopt=0;
    dynamics=0; tidecorr=0;
    adjustfunc=ADJUST_LSA;
    adjRobust=ADJUST_ROBUST_NONE;
    robustRtype=ROBUST_R_DIAGONAL;
    niter=1; codesmooth=0;
    intpref=0; sbascorr=0;
    sbassatsel=0; rovpos=0;
    baspos=0;
    pppfactor=3;
    err[0]=1; err[1]=1; err[2]=0; err[3]=0.01; err[4]=0.0; err[5]=1.0;
    ISBmode=ISBMODE_WHITENOISE;
    std[0]=30.0; std[1]=1; std[2]=0.3;
    stdrate[0]=0; stdrate[1]=1E-6; stdrate[2]=1E-8; stdrate[3]=0;
    sclkstab=5E-12;
    thresar[0]=3.0; thresar[1]=0.9999; thresar[2]=0.25; thresar[3]=0.1; thresar[4]=0.05;
    thresar[5]=0.0; thresar[6]=0.0; thresar[7]=0.0;
    minElArFix=minElArHold=0.0;
    cs_gf=0.05; cs_mw=2.0; maxtdiff=30.0;
    maxres=30.0;
    ru[0]=ru[1]=ru[2]=rb[0]=rb[1]=rb[2]=0.0;

    for (i=0; i<MAX_NF; i++) {
        freq_glo[i]=freq_gal[i]=freq_bds3[i]=freq_bds2[i]=0;
    }

    for (i=0; i<2; i++) {
        baseline[i]=0.0; anttype[i]="\0"; pcvr[i]=gnss_pcv_c();
        antdel[i][0]=antdel[i][1]=antdel[i][3]=0.0; ;
    }

    baseAtt=0;
    posopt[0]=posopt[1]=posopt[2]=0;
    for (i=0; i<MAXSAT; i++) exsats[i]=0;
    maxaveep=0; initrst=0;
    rnxopt[0]=rnxopt[1]="\0";
    syncsol=0;
    freqopt=0;
    pppopt="\0";
}
gnss_prcopt_c::~gnss_prcopt_c() {
}
/* Implementation functions ------------------------------------------------------- */

/* solution options class ----------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
* defaults solution output options ---------------------------------------------------------------- */
gnss_solopt_c::gnss_solopt_c() {
    posf=SOLF_LLH; times=TIMES_GPST;
    timef=1; timeu=3;
    degf=0; outhead=1;
    outopt=0; datum=0;
    height=0; geoid=0;
    solstatic=0; sstat=0;
    trace=0;
    nmeaintv[0]=0.0; nmeaintv[1]=0.0;
    sep=" "; prog="";
}
gnss_solopt_c::~gnss_solopt_c() {
}
/* Implementation functions ------------------------------------------------------- */
/* write solution header to output stream ----------------------------------------- */
int gnss_solopt_c::outsolheads(const gnss_prcopt_c *prcopt,unsigned char *buff) {
    int i;
    char *p=(char *)buff;
    int t_decimal=timeu<0 ? 0 : (timeu>12 ? 13 : timeu+1);

    if (posf==SOLF_NMEA||posf==SOLF_STAT||posf==SOLF_GSIF) {
        return 0;
    }
    if (outhead) {
        p+=sprintf(p,"%s Program Version       : GNSSNAV ver %s\n",COMMENTH,VER_THIS);
        p+=sprintf(p,"%s Positioning Mode      : %s\n",COMMENTH,s_posmode[prcopt->mode]);
        if ( prcopt->mode>=PMODE_PPP_KINEMA && prcopt->mode<=PMODE_PPP_FIXED ) {
            p+=sprintf(p,"%s PPP Observation Mode  : %s\n",COMMENTH,s_pppmod[prcopt->pppmode]);
        }
        p+=sprintf(p,"%s Frequency             : %s\n",COMMENTH,s_freq[prcopt->nf-1]);
        if (prcopt->mode>PMODE_SINGLE) {
            p+=sprintf(p,"%s Solution Type         : %s\n",COMMENTH,s_soltype[prcopt->soltype]);
        }
        p+=sprintf(p,"%s Elevation Mask        : %.1f deg\n",COMMENTH,prcopt->min_elevation*R2D);
        p+=sprintf(p,"%s Min Const Elevation   : %.1f deg\n",COMMENTH,prcopt->con_weight_el*R2D);
        p+=sprintf(p,"%s Dynamical Parameters  : %s\n",COMMENTH,prcopt->dynamics?"on":"off");
        if (prcopt->mode>PMODE_SINGLE) {
            p+=sprintf(p,"%s Tide Corrections      : %s\n",COMMENTH,prcopt->tidecorr?"on":"off");
        }
        //ionosphere & troposphere mode
        p+=sprintf(p,"%s SPP Ionosphere Mode   : %s\n",COMMENTH,s_sppion[prcopt->sppiono]);
        p+=sprintf(p,"%s SPP Troposphere Mode  : %s\n",COMMENTH,s_spptro[prcopt->spptrop]);
        if (prcopt->mode>PMODE_SINGLE) {
            p+=sprintf(p,"%s Ionosphere Mode       : %s\n",COMMENTH,s_ionopt[prcopt->ionoopt]);
            p+=sprintf(p,"%s Troposphere Mode      : %s\n",COMMENTH,s_troopt[prcopt->tropopt]);
        }
        p+=sprintf(p,"%s Ephemeris Mode        : %s\n",COMMENTH,s_sateph[prcopt->sateph]);
        //navigation system
        p+=sprintf(p,"%s Navigation Systems    :",COMMENTH);
        for (i=0; sys[i]; i++) {
            if (prcopt->navsys&sys[i]) p+=sprintf(p," %s",s_navsys[i]);
        }
        p+=sprintf(p,"\n");
        //BDS configuration
        if (prcopt->navsys&SYS_BDS) {
            p+=sprintf(p,"%s BDS Generation        :",COMMENTH);
            for (i=0; iBDSGEN[i]; i++) {
                if (prcopt->bdsGen&iBDSGEN[i]) p+=sprintf(p," %s",s_BDSGEN[i]);
            }
            p+=sprintf(p,"\n");
            p+=sprintf(p,"%s BDS Orbit Type        :",COMMENTH);
            for (i=0; iBDSOBT[i]; i++) {
                if (prcopt->bdsOrbit&iBDSOBT[i]) p+=sprintf(p," %s",s_BDSOBT[i]);
            }
            p+=sprintf(p,"\n");
        }
        //minimum consecutive number to use satellite
        p+=sprintf(p,"%s Min Consecutive Locks : %d\n",COMMENTH,prcopt->minConsecutive);
        //code residual threshold
        p+=sprintf(p,"%s Code Max Residual     : %.1f m\n",COMMENTH,prcopt->maxres);
        //ambiguity
        if ( prcopt->mode!=PMODE_SINGLE && prcopt->mode!=PMODE_DGPS ) {
            p+=sprintf(p,"%s Ambiguity Mode        : %s\n",COMMENTH,s_armode[prcopt->modear]);
            if (prcopt->navsys&SYS_GLO) {
                p+=sprintf(p,"%s GLONASS IFB Mode      : %s\n",COMMENTH,s_glomode[prcopt->glomodear]);
            }
            if (prcopt->thresar[0]>0.0) {
                p+=sprintf(p,"%s Fix Ambguity Ratio    : %.1f\n",COMMENTH,prcopt->thresar[0]);
            }
            //cycle-slip mode
            p+=sprintf(p,"%s Cycle-Slip Mode       : LLI",COMMENTH);
            for (i=0; CSmode[i]; i++) {
                if (prcopt->slipmode&CSmode[i]) p+=sprintf(p," %s",s_csmode[i]);
            }
            p+=sprintf(p,"\n");
            if (prcopt->slipmode&CSMODE_GF) p+=sprintf(p,"%s CS GF Threshold       : %.3f\n",COMMENTH,prcopt->cs_gf);
            if (prcopt->slipmode&CSMODE_MW) p+=sprintf(p,"%s CS MW Threshold       : %.3f\n",COMMENTH,prcopt->cs_mw);
        }
        if (prcopt->mode>PMODE_SINGLE) {
            //filter mode
            p+=sprintf(p,"%s Adjustment Mode       : %s\n",COMMENTH,s_adjmode[prcopt->adjustfunc]);
        }
        p+=sprintf(p,"%s Adj. Robust Mode      : %s\n",COMMENTH,s_adjrobust[prcopt->adjRobust]);
        //baseline
        if ( prcopt->mode==PMODE_MOVEB && prcopt->baseline[0]>0.0 ) {
            p+=sprintf(p,"%s Baseline Length       : %.4f %.4f m\n",COMMENTH,
                prcopt->baseline[0],prcopt->baseline[1]);
        }
        //observation noise mode
        if ( prcopt->mode>=PMODE_PPP_KINEMA && prcopt->mode<=PMODE_PPP_FIXED ) {
            p+=sprintf(p,"%s PPP Obs Noise Factor  : %.3f\n",COMMENTH,prcopt->pppfactor);
        }
        p+=sprintf(p,"%s Code Constant Noise   : %.3f m\n",COMMENTH,prcopt->err[0]);
        p+=sprintf(p,"%s Code Elevation Noise  : %.3f m\n",COMMENTH,prcopt->err[1]);
        if ( prcopt->mode>PMODE_SINGLE ) {
            p+=sprintf(p,"%s Phase Constant Noise  : %.3f m\n",COMMENTH,prcopt->err[2]);
            p+=sprintf(p,"%s Phase Elevation Noise : %.3f m\n",COMMENTH,prcopt->err[3]);
        }
        if (prcopt->dynamics>0) p+=sprintf(p,"%s Doppler Noise         : %.1f Hz\n",COMMENTH,prcopt->err[5]);
        //ambiguity
        if ( prcopt->mode!=PMODE_SINGLE && prcopt->mode!=PMODE_DGPS ) {
            p+=sprintf(p,"%s Ambiguity Noise       : %.2f m   %.1e m^2/s\n",COMMENTH,prcopt->std[0],prcopt->stdrate[0]);
        }
        //ionosphere
        if ( prcopt->mode>PMODE_SINGLE && prcopt->ionoopt==IONOOPT_CONST ) {
            p+=sprintf(p,"%s Ionosphere Noise      : %.2f m   %.1e m^2/s\n",COMMENTH,prcopt->std[1],prcopt->stdrate[1]);
        }
        //troposhpere
        if ( prcopt->mode>PMODE_SINGLE && (prcopt->tropopt==TROPOPT_EST || prcopt->tropopt==TROPOPT_ESTG) ) {
            p+=sprintf(p,"%s ZTD Noise             : %.2f m   %.1e m^2/s\n",COMMENTH,prcopt->std[2],prcopt->stdrate[2]);
        }
        //receiver clock ISB
        if ( prcopt->mode>=PMODE_PPP_KINEMA && prcopt->mode<=PMODE_PPP_FIXED && prcopt->pppmode==PPP_UN_DIFF ) {
            p+=sprintf(p,"%s Receiver ISB Mode     : %s\n",COMMENTH,s_ISBmode[prcopt->ISBmode]);
            if (prcopt->ISBmode==ISBMODE_RANDOMWALK) {
                p+=sprintf(p,"%s ISB Random Walk       : %.1e m^2/s\n",COMMENTH,prcopt->stdrate[3]);
            }
        }
        //antenna
        for (i=0; i<2; i++) {
            if (i>=1 && (prcopt->mode<PMODE_DGPS)) continue;
            p+=sprintf(p,"%s Antenna Type%d         : %-21s (%7.4f %7.4f %7.4f)\n",COMMENTH,
                i+1,prcopt->anttype[i].c_str(),prcopt->antdel[i][0],prcopt->antdel[i][1],
                prcopt->antdel[i][2]);
        }
        //reference station position
        if ( (PMODE_DGPS<=prcopt->mode&&prcopt->mode<=PMODE_FIXED) || posf==SOLF_ENU ) {
            p+=sprintf(p,"%s Reference Position    :",COMMENTH);
            double refpos[3],pos[3],dmsLon[3],dmsLat[3];
            if (refsta==REFSAT_ROVER) for (i=0; i<3; i++) refpos[i]=prcopt->ru[i];
            else for (i=0; i<3; i++) refpos[i]=prcopt->rb[i];

            if (posf==SOLF_LLH||posf==SOLF_ENU) {
                ecef2blh(refpos,datum,pos);
                if (degf) {
                    deg2dms(pos[0]*R2D,dmsLat,5);
                    deg2dms(pos[1]*R2D,dmsLon,5);
                    p+=sprintf(p,"%4.0f %02.0f %08.5f %3.0f %02.0f %08.5f %10.4f",
                        dmsLon[0],dmsLon[1],dmsLon[2],dmsLat[0],dmsLat[1],
                        dmsLat[2],pos[2]);
                }
                else {
                    p+=sprintf(p,"%14.9f %13.9f %10.4f",pos[1]*R2D,pos[0]*R2D,
                        pos[2]);
                }
                p+=sprintf(p," (%14.4f %14.4f %14.4f)\n",refpos[0],refpos[1],refpos[2]);
            }
            else if (posf==SOLF_XYZ) {
                p+=sprintf(p,"%14.4f %14.4f %14.4f\n",refpos[0],refpos[1],refpos[2]);
            }
        }
        p+=sprintf(p,"%s\n",COMMENTH);

        p+=sprintf(p,"%s (",COMMENTH);
        if (posf==SOLF_XYZ) p+=sprintf(p,"x/y/z-ecef=WGS84");
        else if (posf==SOLF_ENU) p+=sprintf(p,"e/n/u-baseline=WGS84");
        else p+=sprintf(p,"lat/lon/height=%s/%s",s1[datum],s2[height]);
        p+=sprintf(p,",Q=1:fix,2:float,3:sbas,4:dgps,5:spp,ns=# of satellites)\n");
    }
    p+=sprintf(p,"%s  %-*s ",COMMENTH,(timef ? 16 : 8)+t_decimal,s3[times]);

    if (posf==SOLF_LLH) { /* lon/lat/hgt */
        if (degf) {
            p+=sprintf(p,"%16s %16s %12s %3s %5s %3s %8s %8s %8s %8s %8s %8s",
                "longitude(d'\")", "latitude(d'\")", "height(m)",
                 "Q", "nsol", "ns", 
                "sde(m)", "sdn(m)", "sdu(m)", 
                "sden(m)", "sdnu(m)", "sdue(m)");
        }
        else {
            p+=sprintf(p,"%14s %14s %12s %3s %5s %3s %8s %8s %8s %8s %8s %8s",
                "longitude(deg)", "latitude(deg)", "height(m)",
                 "Q", "nsol", "ns", 
                "sde(m)", "sdn(m)", "sdu(m)", 
                "sden(m)", "sdnu(m)", "sdue(m)");
        }
        if (prcopt->dynamics) { /* ve/vn/vu */
            p+=sprintf(p,"%12s %12s %12s %10s %10s %10s %10s %10s %10s",
                "ve(m/s)", "vn(m/s)", "vu(m/s)", 
                "sdve(m/s)", "sdvn(m/s)", "sdu(m/s)", 
                "sdven(m/s)", "sdvnu(m/s)", "sdvue(m/s)");
        }
        p+=sprintf(p,"\n");
    }
    else if (posf==SOLF_XYZ) { /* x/y/z-ecef */
        p+=sprintf(p,"%14s %14s %14s %3s %5s %3s %8s %8s %8s %8s %8s %8s",
            "x-ecef(m)", "y-ecef(m)", "z-ecef(m)", "Q",
             "nsol", "ns", "sdx(m)", 
            "sdy(m)", "sdz(m)", "sdxy(m)", 
            "sdyz(m)", "sdzx(m)");
        if (prcopt->dynamics) { /* vx/vy/vz */
            p+=sprintf(p,"%12s %12s %12s %10s %10s %10s %10s %10s %10s",
                "vx(m/s)", "vy(m/s)", "vz(m/s)", 
                "sdvx(m/s)", "sdvy(m/s)", "sdvz(m/s)", 
                "sdvxy(m/s)", "sdvyz(m/s)", "sdvzx(m/s)");
        }
        p+=sprintf(p,"\n");
    }
    else if (posf==SOLF_ENU) { /* e/n/u-baseline */
        p+=sprintf(p,"%14s %14s %14s %3s %5s %3s %8s %8s %8s %8s %8s %8s",
            "e-baseline(m)", "n-baseline(m)", "u-baseline(m)",
             "Q", "nsol", "ns", 
            "sde(m)", "sdn(m)", "sdu(m)", 
            "sden(m)", "sdnu(m)", "sdue(m)");
        if (prcopt->dynamics) { /* ve/vn/vu */
            p+=sprintf(p,"%12s %12s %12s %10s %10s %10s %10s %10s %10s",
                "ve(m/s)", "vn(m/s)", "vu(m/s)", 
                "sdve(m/s)", "sdvn(m/s)", "sdu(m/s)", 
                "sdven(m/s)", "sdvnu(m/s)", "sdvue(m/s)");
        }
        /* pitch/roll/yaw */
        if ( prcopt->mode==PMODE_MOVEB || ( prcopt->mode>=PMODE_DGPS && prcopt->baseAtt ) ) { 
            p+=sprintf(p,"%11s %11s %11s %10s %10s %10s",
                "pitch(deg)", "roll(deg)", "yaw(deg)", 
                "stdP(deg)", "stdR(deg)", "stdY(deg)" );
        }
        p+=sprintf(p,"\n");
    }
    return p-(char *)buff;
}

/* file options class --------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
gnss_filopt_t::gnss_filopt_t() {
}
gnss_filopt_t::~gnss_filopt_t() {
}

/* RINEX options class -------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
gnss_rnxopt_c::gnss_rnxopt_c() {
}
gnss_rnxopt_c::~gnss_rnxopt_c() {
}

/* all options ---------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
all_option_c::all_option_c() {
    pstopt=gnss_pstopt_c();
    rtkopt=gnss_rtkopt_c();
    prcopt=gnss_prcopt_c();
    solopt[0]=solopt[1]=gnss_solopt_c();
    filopt=gnss_filopt_t();
    rnxopt=gnss_rnxopt_c();
}
all_option_c::~all_option_c() {

}
/* reset static system options to default ----------------------------------------- */
void all_option_c::resetsysopts() {
    int i,j;

    rtkopt_=gnss_rtkopt_c();
    prcopt_=gnss_prcopt_c();
    solopt1_=gnss_solopt_c();
    solopt2_=gnss_solopt_c();
    filopt_=gnss_filopt_t();
    for (i=0; i<2; i++) antpostype_[i]=0;
    elmask_=15;
    min_cw_el_=45;
    minElArFix_=minElArHold_=0.0;
    for (i=0; i<2; i++) for (j=0; j<3; j++) {
        antpos_[i][j]=0.0; snrmask_[j]="";
    }
    exsats_="";
    freq_glo_="";
    freq_gal_="";
    freq_bds3_="";
    freq_bds2_="";
}
/* string option to enum (int) ---------------------------------------------------- */
int all_option_c::str2enum(const string str,const string comment,int *val) {
    size_t strp,nlen=0;
    if ((strp=comment.find(str))==string::npos) return 0;
    if (comment[strp-1]!=':') return 0;
    while (strp>=nlen+2) {
        if (comment[strp-2-nlen]>='0'&&comment[strp-2-nlen]<='9')
            nlen++;
        else break;
    }
    *val=stoi(comment.substr(strp-1-nlen,nlen));
    return 1;
}
/* enum (int) to string option ---------------------------------------------------- */
int all_option_c::enum2str(string &str,const string comment,int val) {
    size_t strp1,strp2;
    string num=to_string(val)+':';
    if ((strp1=comment.find(num))==string::npos) { str=to_string(val); return 0; }
    if ((strp2=comment.find(',',strp1+num.length()))==string::npos&&
        (strp2=comment.find(')',strp1+num.length()))==string::npos) {
        str=comment.substr(strp1+num.length());
        return 1;
    }
    else str=comment.substr(strp1+num.length(),strp2-strp1-num.length());
    return 1;
}
/* discard space characters at tail ----------------------------------------------- */
void all_option_c::chop(string &str) {
    size_t strp;
    if ((strp=str.find('#'))!=string::npos) str.erase(strp); /* '#' means comment */
    while (str.length()>0&&!isgraph((int)str.back())) str.pop_back();
}
/* search option ---------------------------------------------------------------------
* search option record
* args   : string name		I  option name (const)
*          opt_t  *opts     I  options table
*                              (terminated with table[i].name="")
* return : option record (NULL: not found)
-----------------------------------------------------------------------------------*/
opt_t * all_option_c::searchopt(const string name,const opt_t *opts) {
    for (int i=0; opts[i].name!=""; i++)
        if (opts[i].name.find(name)!=string::npos) return (opt_t *)(opts+i);
    return NULL;
}
/* string to option value --------------------------------------------------------- */
int all_option_c::str2opt(opt_t *opt,const string str) {
    switch (opt->format) {
        case 0: *(int	 *)opt->var=stoi(str); break;
        case 1: *(double *)opt->var=stold(str); break;
        case 2: *(string *)opt->var=str; break;
        case 3: return str2enum(str,opt->comment,(int *)opt->var); break;
        default: return 0;
    }
    return 1;
}
/* load options ----------------------------------------------------------------------
* load options from file
* args   : string   file    I  options file
*          opt_t  *opts     IO options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*---------------------------------------------------------------------------------- */
int all_option_c::loadopts(const string file,opt_t *opts) {
    ifstream inf;
    opt_t *opt;
    string buff;
    int n=0;
    size_t strp;

    inf.open(file,ios::in);
    if (!inf.is_open()) return 0;

    while (getline(inf,buff)) {
        n++;
        chop(buff);

        if (buff[0]=='\0') continue;

        if ((strp=buff.find('='))==string::npos) continue;

        string name=buff.substr(0,strp),value=buff.substr(strp+1);
        chop(name);

        if (!(opt=searchopt(name,opts))) continue;

        if (!str2opt(opt,value)) continue;
    }

    inf.close();

    return 1;
}
/* system options buffer to options ----------------------------------------------- */
void all_option_c::buff2sysopts() {
    double pos[3]={ 0 },*rr;
    string buff,id;
    int i,j,sat,p,pe,*ps;

    prcopt_.min_elevation = elmask_    *D2R;
    prcopt_.con_weight_el = min_cw_el_ *D2R;
    prcopt_.minElArFix  = minElArFix_*D2R;
    prcopt_.minElArHold = minElArHold_*D2R;

    /* observation noise */
    if (prcopt_.pppfactor<=0) prcopt_.pppfactor=3;

    /* antenna position */
    for (i=0; i<2; i++) {
        ps=i==0 ? &prcopt_.rovpos : &prcopt_.baspos;
        rr=i==0 ? prcopt_.ru : prcopt_.rb;

        if (antpostype_[i]==0) { /* lat/lon/hgt */
            pos[0]=antpos_[i][0]*D2R;
            pos[1]=antpos_[i][1]*D2R;
            pos[2]=antpos_[i][2];
            blh2ecef(pos,WGS84,rr);
        }
        else if (antpostype_[i]==1) { /* xyz-ecef */
            rr[0]=antpos_[i][0];
            rr[1]=antpos_[i][1];
            rr[2]=antpos_[i][2];
        }
        *ps=antpostype_[i];
    }
    /* excluded satellites */
    for (i=0; i<MAXSAT; i++) prcopt_.exsats[i]=0;
    if (exsats_.length()>2) {
        for (p=0,j=0; p<exsats_.length(); p++) {
            if (exsats_[p]==' ') continue; //if the first character is space
            pe=p;
            if ((p=exsats_.substr(p).find(','))==string::npos)
                p=exsats_.length();
            else p+=pe;
            if (exsats_[pe]=='+') id=exsats_.substr(pe+1,p-pe-1); else id=exsats_.substr(pe,p-pe);
            if (!(sat=satid2no(id))) continue;
            prcopt_.exsats[sat-1]=exsats_[pe]=='+'?2:1;
        }
    }

    /* snrmask */
    for (i=0; i<NFREQ; i++) {
        for (j=0; j<9; j++) prcopt_.snrmask.mask[i][j]=0.0;
        for (p=0,j=0; p<snrmask_[i].length(); p++) {
            if (j>8) break;
            pe=p;
            if ((p=snrmask_[i].substr(p).find(','))==string::npos)
                p=snrmask_[i].length();
            else p+=pe;
            str2double(snrmask_[i].substr(pe,p-pe),prcopt_.snrmask.mask[i][j++]);
        }
    }

    /* frequency priority */
    //GLONASS
    for (i=0; i<MAX_NF; i++) prcopt_.freq_glo[i]=0;
    for (p=0,j=0; p<freq_glo_.length(); p++) {
        if (j>=MAX_NF) break;
        pe=p;
        if ((p=freq_glo_.substr(p).find(','))==string::npos)
            p=freq_glo_.length();
        else p+=pe;
        str2int(freq_glo_.substr(pe,p-pe),prcopt_.freq_glo[j++]);
    }
    //Galileo
    for (i=0; i<MAX_NF; i++) prcopt_.freq_gal[i]=0;
    for (p=0,j=0; p<freq_gal_.length(); p++) {
        if (j>=MAX_NF) break;
        pe=p;
        if ((p=freq_gal_.substr(p).find(','))==string::npos)
            p=freq_gal_.length();
        else p+=pe;
        str2int(freq_gal_.substr(pe,p-pe),prcopt_.freq_gal[j++]);
    }
    //BDS-3
    for (i=0; i<MAX_NF; i++) prcopt_.freq_bds3[i]=0;
    for (p=0,j=0; p<freq_bds3_.length(); p++) {
        if (j>=MAX_NF) break;
        pe=p;
        if ((p=freq_bds3_.substr(p).find(','))==string::npos)
            p=freq_bds3_.length();
        else p+=pe;
        str2int(freq_bds3_.substr(pe,p-pe),prcopt_.freq_bds3[j++]);
    }
    //BDS-2
    for (i=0; i<MAX_NF; i++) prcopt_.freq_bds2[i]=0;
    for (p=0,j=0; p<freq_bds2_.length(); p++) {
        if (j>=MAX_NF) break;
        pe=p;
        if ((p=freq_bds2_.substr(p).find(','))==string::npos)
            p=freq_bds2_.length();
        else p+=pe;
        str2int(freq_bds2_.substr(pe,p-pe),prcopt_.freq_bds2[j++]);
    }

    /* number of frequency (4:L1+L5) */
    if (prcopt_.nf>2) {
        prcopt_.nf=3;
        prcopt_.freqopt=1;
    }
    if ( prcopt_.nf<2 && prcopt_.ionoopt==IONOOPT_IFLC ) {
        prcopt_.usedF=2;
    }
    else prcopt_.usedF=prcopt_.nf;
    if (prcopt_.usedF>1) prcopt_.posF=2; // now at most 2 frequencies are allowed for positioning
    else prcopt_.posF=1;

    /* cycle-slip detection */
    if (prcopt_.posF<2) {
        prcopt_.slipmode|=CSMODE_CP;
    }
}
/* get system options ----------------------------------------------------------------
* get system options
* return : none
* notes  : to load system options, use loadopts() before calling the function
*---------------------------------------------------------------------------------- */
void all_option_c::getsysopts() {

    buff2sysopts();

    /* GNSS options */
    pstopt=pstopt_;
    rtkopt=rtkopt_;
    prcopt=prcopt_;
    solopt[0]=solopt1_;
    solopt[1]=solopt2_;
    filopt=filopt_;
}
/* GNSS options ------------------------------------------------------------------- */
/* update input and output file options ------------------------------------------- */
void all_option_c::updatefile(const int days) {
    size_t pathSep,fileDot;
    int nameLength;
    string day_str,wek_str,wwd_str,yer_str,doy_str,hour_str,min_str;
    int week,wwd,year,doy,hour,min;

    int2str(3,"0",days+1,day_str);

    /* update time to the next day */
    if (pstopt.time_start.time>0&&days) pstopt.time_start.timeadd(days*86400);
    pstopt.time_start.time2str(3);
    if (pstopt.time_end.time>0&&days) pstopt.time_end.timeadd(days*86400);
    pstopt.time_end.time2str(3);

    /* if start time !=0 then change the name of inputted GNSS data file*/
    if (pstopt.time_start.time>0) {
        // convert time to wwwwd and dddd.y
        wwd=int(pstopt.time_start.time2gpst(&week)/86400.0);
        year=int(pstopt.time_start.ep[0]); doy=int(pstopt.time_start.time2doy());
        hour=int(pstopt.time_start.ep[3]); min=int(pstopt.time_start.ep[4]);

        /* transform number to string */
        int2str(4,"0",week,wek_str); int2str(1,"0",wwd,wwd_str);
        int2str(4,"0",year,yer_str); int2str(3,"0",doy,doy_str);
        int2str(3,"0",hour,hour_str); int2str(4,"0",min,min_str);

        /* update gnss files according to updated time */
        // roverobs
        pathSep=pstopt.rover_obs.find_last_of(FILEPATHSEP); fileDot=pstopt.rover_obs.find_last_of(".");
        nameLength = pstopt.rover_obs.length() - (pathSep==string::npos? 0 : pathSep+1);
        if (nameLength==SRNX_LEN && fileDot!=string::npos) {//short format
            pstopt.rover_obs.replace(fileDot-4,3,doy_str); //doy
            pstopt.rover_obs.replace(fileDot+1,2,yer_str.substr(2,2)); //yy
        }
        else if (nameLength==LRNXO_LEN && fileDot!=string::npos) { //long format
            pstopt.rover_obs.replace(fileDot-22,4,yer_str); //year
            pstopt.rover_obs.replace(fileDot-18,3,doy_str); //doy
        }

        // baseobs
        if ((prcopt_.mode>=PMODE_DGPS)) {
            pathSep=pstopt.base_obs.find_last_of(FILEPATHSEP); fileDot=pstopt.base_obs.find_last_of(".");
            nameLength = pstopt.base_obs.length() - (pathSep==string::npos? 0 : pathSep+1);
            if (nameLength==SRNX_LEN && fileDot!=string::npos) {//short format
                pstopt.base_obs.replace(fileDot-4,3,doy_str); //doy
                pstopt.base_obs.replace(fileDot+1,2,yer_str.substr(2,2)); //yy
            }
            else if (nameLength==LRNXO_LEN && fileDot!=string::npos) { //long format
                pstopt.base_obs.replace(fileDot-22,4,yer_str); //year
                pstopt.base_obs.replace(fileDot-18,3,doy_str); //doy
            }
        }
        // navigation
        pathSep=pstopt.nav.find_last_of(FILEPATHSEP); fileDot=pstopt.nav.find_last_of(".");
        nameLength = pstopt.nav.length() - (pathSep==string::npos? 0 : pathSep+1);
        if (nameLength==SRNX_LEN && fileDot!=string::npos) {//short format
            pstopt.nav.replace(fileDot-4,3,doy_str); //doy
            pstopt.nav.replace(fileDot+1,2,yer_str.substr(2,2)); //yy
        }
        else if (nameLength==LRNXN_LEN && fileDot!=string::npos) { //long format
            pstopt.nav.replace(fileDot-18,4,yer_str); //year
            pstopt.nav.replace(fileDot-14,3,doy_str); //doy
        }
        // preciseeph
        if (prcopt_.sateph==EPHOPT_PREC&&pstopt.prseph.length()>0) {
            pathSep=pstopt.prseph.find_last_of(FILEPATHSEP); fileDot=pstopt.prseph.find_last_of(".");
            nameLength = pstopt.prseph.length() - (pathSep==string::npos? 0 : pathSep+1);
            if (nameLength==SSATCLK_LEN && fileDot!=string::npos) {//short format
                pstopt.prseph.replace(fileDot-5,4,wek_str); //wwww
                pstopt.prseph.replace(fileDot-1,1,wwd_str); //wwd
            }
            else if (nameLength==LSATCLK_LEN && fileDot!=string::npos) { //long format
                pstopt.prseph.replace(fileDot-23,4,yer_str); //year
                pstopt.prseph.replace(fileDot-19,3,doy_str); //doy
            }
        }
        // satclock
        if (prcopt_.sateph==EPHOPT_PREC&&pstopt.prsclk.length()>0) {
            pathSep=pstopt.prsclk.find_last_of(FILEPATHSEP); fileDot=pstopt.prsclk.find_last_of(".");
            nameLength = pstopt.prsclk.length() - (pathSep==string::npos? 0 : pathSep+1);
            if (nameLength==SSATCLK_LEN && fileDot!=string::npos) {//short format
                pstopt.prsclk.replace(fileDot-5,4,wek_str); //wwww
                pstopt.prsclk.replace(fileDot-1,1,wwd_str); //wwd
            }
            else if (nameLength==LSATCLK_LEN && fileDot!=string::npos) { //long format
                pstopt.prsclk.replace(fileDot-23,4,yer_str); //year
                pstopt.prsclk.replace(fileDot-19,3,doy_str); //doy
            }
        }
        // DCB file
        if (filopt.dcb.length()>0) {
            pathSep=filopt.dcb.find_last_of(FILEPATHSEP); fileDot=filopt.dcb.find_last_of(".");
            nameLength = filopt.dcb.length() - (pathSep==string::npos? 0 : pathSep+1);
            if (nameLength==LBIA_LEN && fileDot!=string::npos) { //long format
                filopt.dcb.replace(fileDot-23,4,yer_str); //year
                filopt.dcb.replace(fileDot-19,3,doy_str); //doy
            }
        }
        // SNX file
        if (filopt.snx.length()>0) {
            pathSep=filopt.snx.find_last_of(FILEPATHSEP); fileDot=filopt.snx.find_last_of(".");
            nameLength = filopt.snx.length() - (pathSep==string::npos? 0 : pathSep+1);
            if (nameLength==SDSNX_LEN && fileDot!=string::npos) { //daily solution file in short format
                filopt.snx.replace(fileDot-5,4,wek_str); //wwww
                filopt.snx.replace(fileDot-1,1,wwd_str); //wwd
            }
            else if (nameLength==SDYSNX_LEN && fileDot!=string::npos) { //daily solution file with year in short format
                filopt.snx.replace(fileDot-8,2,yer_str.substr(2,2)); //yy
                filopt.snx.replace(fileDot-5,4,wek_str); //wwww
                filopt.snx.replace(fileDot-1,1,wwd_str); //wwd
            }
            else if (nameLength==SWSNX_LEN && fileDot!=string::npos) { //weekly solution file in short format
                filopt.snx.replace(fileDot-4,4,wek_str); //wwww
            }
            else if (nameLength==SWYSNX_LEN && fileDot!=string::npos) { //weekly solution file with year in short format
                filopt.snx.replace(fileDot-7,2,yer_str.substr(2,2)); //yy
                filopt.snx.replace(fileDot-4,4,wek_str); //wwww
            }
            else if (nameLength==LSNX_LEN && fileDot!=string::npos) { //long format
                filopt.snx.replace(fileDot-23,4,yer_str); //year
                filopt.snx.replace(fileDot-19,3,doy_str); //doy
            }
        }
    }
}
/* read post options -------------------------------------------------------------- */
int all_option_c::readpostopt(const string file) {
    /* read system options */
    if (!loadopts(file,sysopts)) return 0;
    /* read post options */
    if (!loadopts(file,pstopts)) return 0;

    /* initialize GNSS time configuration */
    pstopt_.set_time_par();
    prcopt_.sampling=pstopt_.time_inter;
    if (prcopt_.sampling>15) prcopt_.minConsecutive=0;

    getsysopts();

    return 1;
}
/* read rtk options --------------------------------------------------------------- */
int all_option_c::readrtkopt(const string file) {
    /* read system options */
    if (!loadopts(file,sysopts)) return 0;
    /* read rtk options */
    if (!loadopts(file,rtkopts)) return 0;

    if (prcopt_.sampling>15) prcopt_.minConsecutive=0;

    getsysopts();

    return 1;
}

/* Integration options ------------------------------------------------------------ */
/* read all post processing options ----------------------------------------------- */
int all_option_c::readpostAll(const string file) {
    /* read GNSS options */
    if (!loadopts(file,sysopts)) return 0;
    if (!loadopts(file,pstopts)) return 0;

    /* initialize GNSS time configuration */
    pstopt_.set_time_par();
    prcopt_.sampling=pstopt_.time_inter;
    if (prcopt_.sampling>15) prcopt_.minConsecutive=0;

    getsysopts();

    return 1;
}
/* read all rtk processing options ------------------------------------------------ */
int all_option_c::readrtkAll(const string file) {
    /* read GNSS options */
    if (!loadopts(file,sysopts)) return 0;
    if (!loadopts(file,rtkopts)) return 0;

    if (prcopt_.sampling>15) prcopt_.minConsecutive=0;

    getsysopts();

    return 1;
}