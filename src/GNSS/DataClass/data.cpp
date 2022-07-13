/* Data Classes File */
/*
 * 2017-7-25 class of observation data
 *
 *
 *
 *
 *
 */
#include "GNSS/DataClass/data.h"
#include "BaseFunction/basefunction.h"
#include "BaseFunction/timesys.h"
#include "ConfigFile/config.h"
#include "GNSS/gnss_pro.h"

 /* class of frequency priority and adoption ------------------------------------------------------ */
gnss_freq_c::gnss_freq_c() {
    priority_glo=0;
    priority_gal=0;
    priority_bds3=0;
    priority_bds2=0;

    // BDS-3, BDS-2
    BDS_B1_index[0][0]=BDS_B1_index[1][0]=0; // C1(B1A/B1C)
    BDS_B1_index[0][1]=BDS_B1_index[1][1]=1; // C2(B1-2)

    init_gnss_freq_priority();
    init_gnss_used_freq();
}
gnss_freq_c::~gnss_freq_c() {
}
/* Implementation functions ------------------------------------------------------- */
/* initialize the GNSS frequency priority with default value ---------------------- */
void gnss_freq_c::init_gnss_freq_priority() {
    int isys,ifreq;
    for (isys=0; isys<MAX_GNSS; isys++) for (ifreq=0; ifreq<2*MAX_NF; ifreq++) {
        GNSS_FREQ_PRI[isys][ifreq]=GNSS_FREQ_PRI_DEF[isys][ifreq];
    }
}
/* initialize the GNSS used frequency with default value -------------------------- */
void gnss_freq_c::init_gnss_used_freq() {
    int isys,ifreq;
    for (isys=0; isys<MAX_GNSS; isys++) for (ifreq=0; ifreq<2*NFREQ; ifreq++) {
        GNSS_USED_FREQ[isys][ifreq]=GNSS_USED_FREQ_DEF[isys][ifreq];
    }
}
/* copy GNSS_USED_FREQ to external matrix ----------------------------------------- */
void gnss_freq_c::copy_used_frep(int ext_freq[MAX_GNSS][NFREQ*2]) {
    int isys,ifreq;
    for (isys=0; isys<MAX_GNSS; isys++) for (ifreq=0; ifreq<2*NFREQ; ifreq++) {
        ext_freq[isys][ifreq]=GNSS_USED_FREQ[isys][ifreq];
    }
}
/* change the GNSS frequency priority with configuration -------------------------- */
void gnss_freq_c::change_freq_priority(const gnss_prcopt_c *opt) {
    int i,ifreq,nfreq=0;
    int setflag[MAX_NF]={0};

    /* GLONASS */
    //get configuration
    priority_glo=0;
    for (i=0; i<MAX_NF; i++) {
        setflag[i]=0;
        GNSS_FREQ_PRI[iGLO][i]=GNSS_FREQ_PRI_DEF[iGLO][i];
    }
    nfreq=0;
    for (i=0; i<MAX_NF; i++) {
        if (opt->freq_glo[i]<1 || (ifreq=vector_index(opt->freq_glo[i],
            &GNSS_FREQ_PRI_DEF[iGLO][0],MAX_NF))<0 ) {
            continue;
        }
        GNSS_FREQ_PRI[iGLO][nfreq++]=opt->freq_glo[i];
        setflag[ifreq]=1;
    }
    priority_glo=nfreq;
    //set left frequency
    for (i=0,ifreq=nfreq; i<MAX_NF&&ifreq<MAX_NF; i++) {
        if (setflag[i]<1) {
            GNSS_FREQ_PRI[iGLO][ifreq++]=GNSS_FREQ_PRI_DEF[iGLO][i];
        }
    }

    /* Galileo */
    //get configuration
    priority_gal=0;
    for (i=0; i<MAX_NF; i++) {
        setflag[i]=0;
        GNSS_FREQ_PRI[iGAL][i]=GNSS_FREQ_PRI_DEF[iGAL][i];
    }
    nfreq=0;
    for (i=0; i<MAX_NF; i++) {
        if (opt->freq_gal[i]<1 || (ifreq=vector_index(opt->freq_gal[i],
            &GNSS_FREQ_PRI_DEF[iGAL][0],MAX_NF))<0 ) {
            continue;
        }
        GNSS_FREQ_PRI[iGAL][nfreq++]=opt->freq_gal[i];
        setflag[ifreq]=1;
    }
    priority_gal=nfreq;
    //set left frequency
    for (i=0,ifreq=nfreq; i<MAX_NF&&ifreq<MAX_NF; i++) {
        if (setflag[i]<1) {
            GNSS_FREQ_PRI[iGAL][ifreq++]=GNSS_FREQ_PRI_DEF[iGAL][i];
        }
    }

    /* BDS-3 */
    //get configuration
    priority_bds3=0;
    for (i=0; i<MAX_NF; i++) {
        setflag[i]=0;
        GNSS_FREQ_PRI[iBDS][i]=GNSS_FREQ_PRI_DEF[iBDS][i];
    }
    nfreq=0;
    for (i=0; i<MAX_NF; i++) {
        if (opt->freq_bds3[i]<1 || (ifreq=vector_index(opt->freq_bds3[i],
            &GNSS_FREQ_PRI_DEF[iBDS][0],MAX_NF))<0 ) {
            continue;
        }
        GNSS_FREQ_PRI[iBDS][nfreq++]=opt->freq_bds3[i];
        setflag[ifreq]=1;
    }
    priority_bds3=nfreq;
    //set left frequency
    for (i=0,ifreq=nfreq; i<MAX_NF&&ifreq<MAX_NF; i++) {
        if (setflag[i]<1) {
            GNSS_FREQ_PRI[iBDS][ifreq++]=GNSS_FREQ_PRI_DEF[iBDS][i];
        }
    }
    //B1 index
    for (i=0; i<MAX_NF; i++) {
        // C1(B1A/B1C)
        if (GNSS_FREQ_PRI[iBDS][i]==1) BDS_B1_index[0][0]=i;
        // C2(B1-2)
        else if (GNSS_FREQ_PRI[iBDS][i]==2) BDS_B1_index[0][1]=i;
    }

    /* BDS-2 */
    //get configuration
    priority_bds2=0;
    for (i=0; i<MAX_NF; i++) {
        setflag[i]=0;
        GNSS_FREQ_PRI[iBDS][MAX_NF+i]=GNSS_FREQ_PRI_DEF[iBDS][MAX_NF+i];
    }
    nfreq=0;
    for (i=0; i<MAX_NF; i++) {
        if (opt->freq_bds2[i]<1 || (ifreq=vector_index(opt->freq_bds2[i],
            &GNSS_FREQ_PRI_DEF[iBDS][MAX_NF],MAX_NF))<0 ) {
            continue;
        }
        GNSS_FREQ_PRI[iBDS][MAX_NF+nfreq++]=opt->freq_bds2[i];
        setflag[ifreq]=1;
    }
    priority_bds2=nfreq;
    //set left frequency
    for (i=0,ifreq=nfreq; i<MAX_NF&&ifreq<MAX_NF; i++) {
        if (setflag[i]<1) {
            GNSS_FREQ_PRI[iBDS][MAX_NF+ifreq++]=GNSS_FREQ_PRI_DEF[iBDS][MAX_NF+i];
        }
    }
    //B1 index
    for (i=0; i<MAX_NF; i++) {
        // C1(B1A/B1C)
        if (GNSS_FREQ_PRI[iBDS][MAX_NF+i]==1) BDS_B1_index[1][0]=MAX_NF+i;
        // C2(B1-2)
        else if (GNSS_FREQ_PRI[iBDS][MAX_NF+i]==2) BDS_B1_index[1][1]=MAX_NF+i;
    }
}
/* change the used GNSS frequency according to priority and availability -------------
* change the used GNSS frequency "GNSS_USED_FREQ[NSYS][NFREQ+NFREQ]"
* args   : const int    sys      I      GNSS constellation
*          const int   *indRov   IO     rover signal index
*          const int   *indBas   IO     base signal index
* return : (NULL: no output)
* notes  : only change the frequency for Galileo and BDS!!!!!!!!
*---------------------------------------------------------------------------------- */
void gnss_freq_c::change_gnss_freq(const int sys,sigind_t *indRov,sigind_t *indBas) {
    int ifreq,isys;
    switch (sys) {
        case SYS_GPS: isys=iGPS; break;
        case SYS_GLO: isys=iGLO; break;
        case SYS_GAL: isys=iGAL; break;
        case SYS_BDS: isys=iBDS; break;
        case SYS_QZS: isys=iQZS; break;
        case SYS_IRN: isys=iIRN; break;
        case SYS_SBS: isys=iSBS; break;
        default: 	  isys=-1; break;
    }

    if (isys<0) return;

    /* change BDS frequency */
    else if (sys==SYS_BDS) {
        /* loop BDS-3 (i=1) and BDS-2 (i=2) */
        for (int i=0; i<2; i++) {
            /* 2nd, 3rd ... frequency (from C5, C6, C7, C8) */
            int nf=0;
            for (ifreq=0; ifreq<MAX_NF; ifreq++) {
                if ( indRov->avaiFrq[ i*MAX_NF + ifreq ]>0 && 
                    ( indBas==NULL || indBas->avaiFrq[ i*MAX_NF + ifreq ]>0 ) ) {
                    GNSS_USED_FREQ[iBDS][ i*NFREQ + nf++ ] = 
                        GNSS_FREQ_PRI[iBDS][ i*MAX_NF + ifreq ];
                }
                if (nf>=NFREQ) break;
            }
            /* reset unavailable frequency to 0 */
            for (ifreq=nf; ifreq<NFREQ; ifreq++) {
                GNSS_USED_FREQ[iBDS][ i*NFREQ + ifreq ]=0;
            }
        }
    }

    /* change other GNSS frequency */
    else {
        int nf=0;
        for (ifreq=0; ifreq<MAX_NF; ifreq++) {
            if ( indRov->avaiFrq[ifreq]>0 && ( indBas==NULL || indBas->avaiFrq[ifreq]>0 ) ) {
                GNSS_USED_FREQ[isys][nf++] = GNSS_FREQ_PRI[isys][ifreq];
            }
            if (nf>=NFREQ) break;
        }
        /* reset unavailable frequency to 0 */
        for (ifreq=nf; ifreq<NFREQ; ifreq++) {
            GNSS_USED_FREQ[isys][ifreq]=0;
        }
    }
}


 /*class of one epoch observation data ------------------------------------------------------------ */
gnss_obsd_c::gnss_obsd_c() {
    sat=prn=sys=dtr=0;
    nCon=0;
    isys=-1;
    time=sigtime=gtime_c();
    for (int i=0; i<NFREQ+NEXOBS; i++) {
        SNR[i]=LLI[i]=code[i]=L[i]=D[i]=P[i]=DCB[i]=res[i]=ovar[i]=0.0;
        LCode[i]=PCode[i]=DCode[i]=SCode[i]="";
    }
    for (int i=0; i<3; i++) {
        sigvec[i]=0.0;
        for (int j=0; j<2; j++) posvel[i*2+j]=0.0;
    }
    used=vused=0;
    for (int i=0; i<NFREQ; i++) { 
        antv[i]=rant[i]=0.0;
        anto[i][0]=anto[i][1]=anto[i][2]=0.0;
        sant[i]=0;
    }
    rant[NFREQ]=0;
    sant[NFREQ]=0;
    dts[0]=dts[1]=azel[0]=azel[1]=0.0;
    rcv=0;
    dist=dcbvar=exc=sat=svar=svh=ionmap=dion=ionvar=dtro=trovar=m_h=m_w=phasewp=0;
    ipp[0]=ipp[1]=0;
    warningMsg="";
}
gnss_obsd_c::~gnss_obsd_c() {
}
/* Implementation functions ------------------------------------------------------- */
/* reset the whole data ----------------------------------------------------------- */
void gnss_obsd_c::reset() {
    int i;
    dtr=0;
    sat='\0';
    prn=sys=0;
    isys=-1;
    time=gtime_c();
    /* observation data */
    for (i=0; i<NFREQ+NEXOBS; i++) {
        SNR[i]=LLI[i]=code[i]='\0';
        L[i]=P[i]=D[i]=DCB[i]=0.0;
        LCode[i]=PCode[i]=DCode[i]=SCode[i]="";
    }
    for (i=0; i<3; i++) res[i]=ovar[i]=0;
    for (int i=0; i<NFREQ; i++) { 
        antv[i]=rant[i]=0.0;
    }
    rant[NFREQ]=0;
    used=vused=0;
    /* satellite data */
    svar=svh=0;
    for (i=0; i<2; i++) {
        azel[i]=dts[i]=0.0;
        for (int j=0; j<3; j++) posvel[i*3+j]=0.0;
    }
    for (int i=0; i<NFREQ; i++) { 
        anto[i][0]=anto[i][1]=anto[i][2]=0.0;
        sant[i]=0;
    }
    sant[NFREQ]=0;
    /* correction data */
    dist=dcbvar=exc=sat=svar=svh=ionmap=dion=ionvar=dtro=trovar=m_h=m_w=phasewp=0;
    ipp[0]=ipp[1]=0;
    warningMsg="";
    errorMsg="";
}
/* reset satellite data ----------------------------------------------------------- */
void gnss_obsd_c::satreset() {
    /* satellite data */
    svar=svh=0;
    for (int i=0; i<2; i++) {
        azel[i]=dts[i]=0.0;
        for (int j=0; j<3; j++) posvel[i*3+j]=0.0;
    }
    for (int i=0; i<NFREQ; i++) { 
        antv[i]=0.0;
        anto[i][0]=anto[i][1]=anto[i][2]=0.0;
        sant[i]=0;
    }
}
/* updata signal time use pseudorange---------------------------------------------- */
int gnss_obsd_c::sigtime_opsr() {
    double pr;
    int i;

    sigtime = time;
    for (i=0,pr=0.0; i<NFREQ; i++) if ((pr=P[i])!=0.0) break;
    if (i<NFREQ) {
        sigtime.timeadd(-pr/CLIGHT);
        return 1;
    }
    else { errorMsg = "no pseudorange"; return 0; }
}
/* update signal time use satellite clock bias ------------------------------------ */
void gnss_obsd_c::sigtime_clk() {
    sigtime.timeadd(-dts[0]);
}
/* update signal time use known position and receiver clock bias ------------------ */
void gnss_obsd_c::sigtime_rec(const vector<double> &pos) {
    double dis=geodist_GNSS(posvel,pos.begin(),sigvec);
    sigtime=time;
    sigtime.timeadd(-dis/CLIGHT-dtr);
}

/* station informtations -------------------------------------------------------------------------- */
gnss_sta_c::gnss_sta_c() {
    name="        ";
    antsetup=itrf=deltype=hgt=0;
    for (int i=0; i<3; i++) pos[i]=del[i]=0.0;
}
gnss_sta_c::~gnss_sta_c() {
}

/* class of station's information and observation chains ------------------------------------------ */
gnss_obs_c::gnss_obs_c() {
    n=used=0;
    rcv=0;
    sta=gnss_sta_c();
}
gnss_obs_c::~gnss_obs_c() {
    data.clear();
}
void gnss_obs_c::reset() {
    n=0;
    data.clear();
}

/* GPS/QZS/GAL broadcast ephemeris class ---------------------------------------------------------- */
gnss_eph_c::gnss_eph_c() {
    sat=iode=iodc=sva=svh=week=code=flag=0;
    toe=toc=ttr=gtime_c();
    A=e=i0=OMG0=omg=M0=deln=OMGd=idot=0;
    crc=crs=cuc=cus=cic=cis=0;
    toes=fit=f0=f1=f2=Adot=ndot=0;
    tgd[0]=tgd[1]=tgd[2]=tgd[3]=0;
    Arate=0;
}
gnss_eph_c::~gnss_eph_c() {
}

/* GLONASS broadcast ephemeris class -------------------------------------------------------------- */
gnss_geph_c::gnss_geph_c() {
    sat=iode=frq=svh=sva=age=0;
    toe=tof=gtime_c();
    for (int i=0; i<3; i++) pos[i]=vel[i]=acc[i]=0;
    taun=gamn=dtaun=0;
}
gnss_geph_c::~gnss_geph_c() {
}

/* SBAS ephemeris class --------------------------------------------------------------------------- */
gnss_seph_c::gnss_seph_c() {
    sat=sva=svh=0;
    t0=tof=gtime_c();
    af0=af1=0;
    for (int i=0; i<3; i++) pos[i]=vel[i]=acc[i]=0;
}
gnss_seph_c::~gnss_seph_c() {
}

/* precise ephemeris class ------------------------------------------------------------------------ */
gnss_peph_c::gnss_peph_c() {
    time=gtime_c();
    for (int i=0; i<MAXSAT; i++) {
        for (int j=0; j<4; j++)
            pos[i][j]=std[i][j]=vel[i][j]=vst[i][j]=0;
        for (int j=0; j<3; j++)
            cov[i][j]=vco[i][j]=0;
    }
}
gnss_peph_c::~gnss_peph_c() {
}

/* precise clock class ---------------------------------------------------------------------------- */
gnss_pclk_c::gnss_pclk_c() {
    time=gtime_c();
    for (int i=0; i<MAXSAT; i++)
        clk[i]=std[i]=0;
}
gnss_pclk_c::~gnss_pclk_c() {
}

/* almanac class ---------------------------------------------------------------------------------- */
gnss_alm_c::gnss_alm_c() {
    sat=svh=svconf=week=0;
    toa=gtime_c();
    A=e=i0=OMG0=omg=M0=OMGd=0;
    toas=f0=f1=0;
}
gnss_alm_c::~gnss_alm_c() {
}

/* TEC grid class --------------------------------------------------------------------------------- */
gnss_tec_c::gnss_tec_c() {
    time=gtime_c();
    rb=0;
    for (int i=0; i<3; i++)
        lats[i]=lons[i]=hgts[i]=ndata[i]=0;
}
gnss_tec_c::~gnss_tec_c() {
    data.clear(); rms.clear();
}

/* GNSS observation bias class -------------------------------------------------------------------- */
gnss_bias_c::gnss_bias_c() {
    int i,j;

    for (i=0; i<4; i++) {
        /* glonass code-phase bias */
        glo_cpbias[i]=0.0;
    }

    for (i=0; i<MAXSAT; i++) {
        P1P2[i]=0;
        for (j=0; j<=MAXFREQ; j++) cbias[i][j]=0.0;
        wlbias[i]=0.0;
    }
    rbias[0][0]=rbias[0][1]=rbias[0][2]=
        rbias[1][0]=rbias[1][1]=rbias[1][2]=0.0;
}
gnss_bias_c::~gnss_bias_c() {
}

/* ZTD data class --------------------------------------------------------------------------------- */
gnss_ztd_c::gnss_ztd_c() {
    rms=ztd=0;
}
gnss_ztd_c::~gnss_ztd_c() {
}

/* satellite fcb data class ----------------------------------------------------------------------- */
gnss_fcbd_c::gnss_fcbd_c() {
    ts=te=gtime_c();
    for (int i=0; i<MAXSAT; i++)
        bias[i][0]=bias[i][1]=bias[i][2]=
        std[i][0]=std[i][1]=std[i][2]=0;
}
gnss_fcbd_c::~gnss_fcbd_c() {
}

/* earth rotation parameter data class ------------------------------------------------------------ */
gnss_erpd_c::gnss_erpd_c() {
    mjd=xp=yp=xpr=ypr=ut1_utc=lod=0;
}
gnss_erpd_c::~gnss_erpd_c() {
}

/* earth rotation parameter class ----------------------------------------------------------------- */
gnss_erp_c::gnss_erp_c() {
    n=nmax=0;
}
gnss_erp_c::~gnss_erp_c() {
    data.clear();
}

/* antenna parameter class ------------------------------------------------------------------------ */
gnss_pcv_c::gnss_pcv_c() {
    sat=prn=nzen=0;
    ts=te=gtime_c();
    for (int i=0; i<NSYS; i++) {
        for (int f=0; f<NFREQ; f++) {
            off[i][f][0]=off[i][f][1]=off[i][f][2]=0;
            var[i][f].clear();
        }
    }
    for (int i=0; i<3; i++) zen[i]=0;
}
gnss_pcv_c::~gnss_pcv_c() {
    for (int i=0; i<NSYS; i++) {
        for (int f=0; f<NFREQ; f++) {
            var[i][f].clear();
        }
    }
}

/* SBAS fast correction class --------------------------------------------------------------------- */
gnss_sbsfcorr_c::gnss_sbsfcorr_c() {
    t0=gtime_c();
    prc=rrc=dt=iodf=udre=ai=0;
}
gnss_sbsfcorr_c::~gnss_sbsfcorr_c() {
}

/* SBAS long term satellite error correction class ------------------------------------------------ */
gnss_sbslcorr_c::gnss_sbslcorr_c() {
    t0=gtime_c();
    iode=daf0=daf1=0;
    dpos[0]=dpos[1]=dpos[2]=dvel[0]=dvel[1]=dvel[2]=0;
}
gnss_sbslcorr_c::~gnss_sbslcorr_c() {
}

/* SBAS satellite correction class ---------------------------------------------------------------- */
gnss_sbssatp_c::gnss_sbssatp_c() {
    sat=0;
    fcorr=gnss_sbsfcorr_c();
    lcorr=gnss_sbslcorr_c();
}
gnss_sbssatp_c::~gnss_sbssatp_c() {
}

/* SBAS satellite corrections class --------------------------------------------------------------- */
gnss_sbssat_c::gnss_sbssat_c() {
    iodp=0;
    nsat=0;
    tlat=0;
    for (int i=0; i<MAXSAT; i++)
        sat[i]=gnss_sbssatp_c();
}
gnss_sbssat_c::~gnss_sbssat_c() {
}

/* SBAS ionospheric correction class -------------------------------------------------------------- */
gnss_sbsigp_c::gnss_sbsigp_c() {
    t0=gtime_c();
    lat=lon=give=delay=0;
}
gnss_sbsigp_c::~gnss_sbsigp_c() {
}

/* SBAS ionospheric corrections class ------------------------------------------------------------- */
gnss_sbsion_c::gnss_sbsion_c() {
    iodi=0;
    nigp=0;
    for (int i=0; i<MAXNIGP; i++)
        igp[i]=gnss_sbsigp_c();
}
gnss_sbsion_c::~gnss_sbsion_c() {
}

/* DGPS/GNSS correction class --------------------------------------------------------------------- */
gnss_dgps_c::gnss_dgps_c() {
    t0=gtime_c();
    prc=rrc=iod=udre=0;
}
gnss_dgps_c::~gnss_dgps_c() {
}

/* SSR correction class --------------------------------------------------------------------------- */
gnss_ssr_c::gnss_ssr_c() {
    t0[0]=t0[1]=t0[2]=t0[3]=t0[4]=t0[5]=gtime_c();
    update=0;
    iode=iodcrc=ura=refd=hrclk=yaw_ang=yaw_rate=0;
    for (int i=0; i<3; i++) {
        deph[i]=ddeph[i]=dclk[i]=0.0;
        for (int j=0; j<2; j++) udi[i*2+j]=iod[i*2+j]=0;
    }
    for (int i=0; i<MAXCODE; i++) cbias[i]=pbias[i]=stdpb[i]=0.0;
}
gnss_ssr_c::~gnss_ssr_c() {
}

/* QZSS LEX ephemeris class ----------------------------------------------------------------------- */
gnss_lexeph_c::gnss_lexeph_c() {
    toe=tof=gtime_c();
    af0=af1=tgd=sat=0;
    health=ura='\0';
    for (int i=0; i<3; i++)
        pos[i]=vel[i]=acc[i]=jerk[i]=0;
    for (int i=0; i<8; i++)
        isc[i]=0;
}
gnss_lexeph_c::~gnss_lexeph_c() {
}

/* QZSS LEX ionosphere correction class ----------------------------------------------------------- */
gnss_lexion_c::gnss_lexion_c() {
    t0=gtime_c();
    pos0[0]=pos0[1]=tspan=0;
    coef[0][0]=coef[1][0]=coef[2][0]=
        coef[0][1]=coef[1][1]=coef[2][1]=0;
}
gnss_lexion_c::~gnss_lexion_c() {
}

/* stec data class -------------------------------------------------------------------------------- */
gnss_stec_c::gnss_stec_c() {
    time=gtime_c();
    sat=flag='\0';
    ion=std=azel[0]=azel[1]=0;
}
gnss_stec_c::~gnss_stec_c() {
}

/* trop data class -------------------------------------------------------------------------------- */
gnss_trop_c::gnss_trop_c() {
    time=gtime_c();
    trp[0]=trp[1]=trp[2]=
        std[0]=std[1]=std[2]=0;
}
gnss_trop_c::~gnss_trop_c() {
}

/* ppp corrections class -------------------------------------------------------------------------- */
gnss_pppcorr_c::gnss_pppcorr_c() {
    nsta=0;
    for (int i=0; i<MAXSAT; i++) {
        rr[i][0]=rr[i][1]=rr[i][2]=0;
        ns[i]=nsmax[i]=nt[i]=ntmax[i]=0;
        stec[i]=NULL;
        trop[i]=NULL;
    }
}
gnss_pppcorr_c::~gnss_pppcorr_c() {
    for (int i=0; i<MAXSAT; i++) {
        stec[i]=NULL;
        trop[i]=NULL;
    }
}

/* class of navigation data ----------------------------------------------------------------------- */
gnss_nav_c::gnss_nav_c() {
    int i,j;
    /* number of ephemeris */
    n=ng=ns=ne=nc=na=nt=nf=0;
    nmax=ngmax=nsmax=nemax=ncmax=namax=ntmax=nfmax=0;

    for (i=0; i<4; i++) {
        /* utc parameters */
        utc_gps[i]=utc_glo[i]=utc_gal[i]=utc_qzs[i]=utc_cmp[i]=
            utc_irn[i]=utc_sbs[i]=0.0;
        /* model parameters */
        ion_gps[i]=ion_gps[i+1]=0.0;
        ion_gal[i]=0.0;
        ion_qzs[i]=ion_qzs[i+1]=0.0;
        ion_cmp[i]=ion_cmp[i+1]=0.0;
        ion_irn[i]=ion_irn[i+1]=0.0;
    }
    for (i=0; i<66; i++) {
        ocean_par[0][i]=ocean_par[1][0]=0.0;
    }
    leaps=0;
    freq=gnss_freq_c();
    erp=gnss_erp_c();
    sbssat=gnss_sbssat_c();
    lexion=gnss_lexion_c();
    pppcorr=gnss_pppcorr_c();
    obsbias=gnss_bias_c();

    /* MAXSAT array */
    for (i=0; i<MAXSAT; i++) {
        for (j=0; j<NFREQ; j++) lam[i][j]=0.0;
        pcvs[i]=gnss_pcv_c();
        dgps[i]=gnss_dgps_c();
        ssr[i]=gnss_ssr_c();
        lexeph[i]=gnss_lexeph_c();
    }
    for (i=0; i<MAXPRNGLO+1; i++) glo_fcn[i]='\0';
    for (i=0; i<MAXBAND+1; i++) sbsion[i]=gnss_sbsion_c();
}
gnss_nav_c::~gnss_nav_c() {
    eph.clear(); geph.clear(); seph.clear(); peph.clear();
    pclk.clear(); alm.clear(); tec.clear(); ztd.clear(); fcb.clear();
}
/* Implementation functions ------------------------------------------------------- */
/* satellite carrier wave length -------------------------------------------------- */
double gnss_nav_c::satwavelen(int SatNum,int FrqNum) {
    const double freq_glos[]={ 0, FREQ1_GLO,FREQ2_GLO };
    const double dfrq_glos[]={ 0, DFRQ1_GLO,DFRQ2_GLO };
    int prn;
    int i,satSys=satsys(SatNum,&prn);

    /* for GPS */
    if (satSys==SYS_GPS) {
        return GNSS_LAMBDA_M[iGPS][ freq.GNSS_USED_FREQ[iGPS][FrqNum] ];
    }
    /* for GLONASS */
    else if (satSys==SYS_GLO) {
        if ( 0<freq.GNSS_USED_FREQ[iGLO][FrqNum] &&
            freq.GNSS_USED_FREQ[iGLO][FrqNum]<3 ) {
            for (i=0; i<ng; i++) {
                if (geph[i].sat!=SatNum) continue;
                return CLIGHT/(freq_glos[ freq.GNSS_USED_FREQ[iGLO][FrqNum] ]+
                    dfrq_glos[ freq.GNSS_USED_FREQ[iGLO][FrqNum] ]*geph[i].frq);
            }
        }
        else { /* G3,G4,G6 */
            return GNSS_LAMBDA_M[iGLO][ freq.GNSS_USED_FREQ[iGLO][FrqNum] ];
        }
    }
    /* for Galileo */
    else if (satSys==SYS_GAL) {
        return GNSS_LAMBDA_M[iGAL][ freq.GNSS_USED_FREQ[iGAL][FrqNum] ];
    }
    /* for BDS */
    else if (satSys==SYS_BDS) {
        /* BDS-3 */
        if (prn>=BDS3_MIN) return GNSS_LAMBDA_M[iBDS][ freq.GNSS_USED_FREQ[iBDS][FrqNum] ];
        /* BDS-2 */
        else return GNSS_LAMBDA_M[iBDS][ freq.GNSS_USED_FREQ[iBDS][ NFREQ+FrqNum ] ];
    }
    /* for QZSS */
    else if (satSys==SYS_QZS) {
        return GNSS_LAMBDA_M[iQZS][ freq.GNSS_USED_FREQ[iQZS][FrqNum] ];
    }
    /* for IRNSS */
    else if (satSys==SYS_IRN) {
        return GNSS_LAMBDA_M[iIRN][ freq.GNSS_USED_FREQ[iIRN][FrqNum] ];
    }
    /* for SBAS */
    else if (satSys==SYS_SBS) {
        return GNSS_LAMBDA_M[iSBS][ freq.GNSS_USED_FREQ[iSBS][FrqNum] ];
    }

    return 0.0;
}
/* update satellite carrier wave lengths ------------------------------------------ */
void gnss_nav_c::update_sat_lambda() {
    int i,j;
    for (i=0; i<MAXSAT; i++) for (j=0; j<NFREQ; j++) {
        lam[i][j]=satwavelen(i+1,j);
    }
}
void gnss_nav_c::update_sat_lambda(vector<gnss_ssat_c> &ssat) {
    if (ssat.size()<MAXSAT) return;
    for (int i=0; i<MAXSAT; i++) {
        for (int j=0; j<NFREQ; j++) {
            lam[i][j]=satwavelen(i+1,j);
            ssat[i].lambda[j]=lam[i][j];
        }
    }
}

/* SBAS message class ----------------------------------------------------------------------------- */
gnss_sbsmsg_c::gnss_sbsmsg_c() {
    week=tow=prn=0;
    for (int i=0; i<29; i++)
        msg[i]='\0';
}
gnss_sbsmsg_c::~gnss_sbsmsg_c() {
}

/* SBAS messages class ---------------------------------------------------------------------------- */
gnss_ssbs_c::gnss_ssbs_c() {
    n=nmax=0;
}
gnss_ssbs_c::~gnss_ssbs_c() {
    msgs.clear();
}

/* QZSS LEX message class ------------------------------------------------------------------------- */
gnss_lexmsg_c::gnss_lexmsg_c() {
    prn=type=alert=ttt=0;
    stat=snr='\0';
    for (int i=0; i<212; i++)
        msg[i]='\0';
}
gnss_lexmsg_c::~gnss_lexmsg_c() {
}
/* QZSS LEX messages class ------------------------------------------------------------------------ */
gnss_lex_c::gnss_lex_c() {
    n=nmax=0;
}
gnss_lex_c::~gnss_lex_c() {
    msgs.clear();
}
