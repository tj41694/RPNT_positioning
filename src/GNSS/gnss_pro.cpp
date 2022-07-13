#include "GNSS/gnss_pro.h"

#include "BaseFunction/basefunction.h"
#include "GNSS/PosModel/position.h"

#include "Adjustment/adjustment.h"
#include "GNSS/AmbModel/ambiguity.h"
#include "GNSS/EphModel/satellite.h"
#include "GNSS/AntModel/antenna.h"
#include "GNSS/TidModel/tide.h"
#include "GNSS/IonModel/ionosphere.h"
#include "GNSS/TroModel/troposphere.h"
#include "GNSS/ParModel/parameter.h"

/* raw_c */
#include "Decode/rtcm.h"
#include "Decode/raw/binex.h"
#include "Decode/raw/cmr.h"
#include "Decode/raw/crescent.h"
#include "Decode/raw/gw10.h"
#include "Decode/raw/javad.h"
#include "Decode/raw/novatel.h"
#include "Decode/raw/nvs.h"
#include "Decode/raw/septentrio.h"
#include "Decode/raw/skytraq.h"
#include "Decode/raw/superstar2.h"
#include "Decode/raw/trimble17.h"
#include "Decode/raw/ublox.h"

/* constant --------------------------------------------------------------------------------------- */
#define SQR(x)          ((x)<0.0?-(x)*(x):(x)*(x))
#define SQRT(x)         ((x)<0.0?0.0:sqrt(x))
#define MAXION_DT       600.0                       /* max ion-change time (s) */
#define NESTSLIP        999                         /* flag of no estimated cycle slip */
#define MAXITR          10					        /* max number of iteration for point pos */
/* sqrt of covariance ------------------------------------------------------------- */
static double sqvar(double covar) {
    return covar<0.0 ? -sqrt(-covar) : sqrt(covar);
}

/* GNSS process class --------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
gnss_pro_c::gnss_pro_c() {
    /* processing status */
    soltime[0]=soltime[1]=soltime[2]=gtime_c();
    nsol=ntotal=0;
    stat=vstat=0;
    numsys=ns=0;
    /* num */
    neb=nfix=0;
    NF=NP=NV=NT=NG=NI=NC=NU=NA=numF=NX=N_ALL=numX=0;
    CJ_nfreq=2;
    numL=vnumL=0;
    rb[0]=rb[1]=rb[2]=rb[3]=rb[4]=rb[5]=tt=0.0;
    ext_posvel[0]=ext_posvel[1]=ext_posvel[2]=ext_posvel[3]=ext_posvel[4]=ext_posvel[5]=0.0;
    for (int i=0; i<9; i++) {
        ext_vpos[i]=ext_vvel[i]=0;
    }
    att_n_b[0]=att_n_b[1]=att_n_b[2]=0;
    att_err[0]=att_err[1]=att_err[2]=0;
    varFlag=0;
    preClkFlag=0;
    /* class */
    opt=NULL;
    obsr=NULL; obsb=NULL;
    nav=NULL;
    /* string */
    msg = "\0";
    usePhase=0;
    ssat.assign(MAXSAT,gnss_ssat_c());
    /* initialize function classes */
    satantfunc=gnss_satant_c();
    recantfunc=gnss_recant_c();
    parafunc=gnss_parameter_c();
    satfunc=NULL;
    sppionf=NULL;
    spptrof=NULL;
    ionfunc=NULL;
    trofunc=NULL;
    adjfunc=NULL;
}
gnss_pro_c::~gnss_pro_c() {
    opt=NULL; nav=NULL;
    sol.clear(); b_sol.clear();
    if (log_stream.is_open()) log_stream.close();

    Xest.clear(); Xini.clear();
    X_ALL.clear();
    Rx_ALL.clear();
    Lobs.clear(); Acoe.clear(); Rvar.clear();
    Rvec.clear(); Xpar.clear(); Xori.clear();
    Rx.clear();   QAA.clear();
    vLobs.clear(); vAcoe.clear(); vRvar.clear();
    vRvec.clear(); vXpar.clear(); vRx.clear();

    if (sppionf) { delete sppionf; sppionf=NULL; }
    if (spptrof) { delete spptrof; spptrof=NULL; }
    if (satfunc) { delete satfunc; satfunc=NULL; }
    if (ionfunc) { delete ionfunc; ionfunc=NULL; }
    if (trofunc) { delete trofunc; trofunc=NULL; }
    if (adjfunc) { delete adjfunc; adjfunc=NULL; }
}
/* Implementation functions ------------------------------------------------------- */
/* initialize gnss process -------------------------------------------------------- */
void gnss_pro_c::gnss_pro_init() {
    /* initialize parameters */
    NF=parafunc.N_Freqency(opt);
    NP=parafunc.N_Position(opt);
    NV=parafunc.N_Velocity(opt);
    NT=parafunc.N_Tro(opt);
    NG=parafunc.N_GLOIFB(opt);
    NI=parafunc.N_Ion(opt);
    NC=parafunc.N_Clock(opt);
    NX=N_ALL=NP+NT+NC+NI+NG;
    numF = opt->ionoopt==IONOOPT_IFLC ? 1 : NF;
    NA=parafunc.N_Amb(opt,numF);
    N_ALL+=NA;
    if (opt->usedF>2) {
        if(opt->usedF<NFREQ) CJ_nfreq=opt->usedF;
        else CJ_nfreq=NFREQ;
    }
    else CJ_nfreq=2;
    X_ALL.assign(N_ALL,0.0);
    Rx_ALL.assign(N_ALL*N_ALL,0.0);
    /* initialize functions ---------------- */
    /* satellite ephemeris functions */
    switch (opt->sateph) {
        case EPHOPT_BRDC:   satfunc=new gnss_broadcast_c;   break;
        case EPHOPT_PREC:   satfunc=new gnss_preciseph_c;   break;
        case EPHOPT_SBAS:   satfunc=new gnss_broadsbas_c;   break;
        case EPHOPT_SSRAPC: satfunc=new gnss_broadssrapc_c; break;
        case EPHOPT_SSRCOM: satfunc=new gnss_broadssrcom_c; break;
        case EPHOPT_LEX:    satfunc=new gnss_qzsslex_c;     break;
        default: satfunc=new gnss_satellite_c;
    }
    satfunc->satantfunc=&satantfunc;
    /* adjustment functions */
    switch (opt->adjustfunc) {
        case ADJUST_LSA:
             adjfunc=new lsadj_c(opt->adjRobust,opt->robustRtype); break;
        case ADJUST_KALMAN:
             adjfunc=new kalmanfilter_c(opt->adjRobust,opt->robustRtype); break;
        case ADJUST_HELMERT:
             adjfunc=new helmert_c(opt->adjRobust,opt->robustRtype); break;
        default: adjfunc=new kalmanfilter_c(opt->adjRobust,opt->robustRtype);
    }

    /* SPP ionosphere function */
    switch (opt->sppiono) {
        case IONOOPT_IFLC:  sppionf=new gnss_LCion_c;     break;
        case IONOOPT_BRDC:  sppionf=new gnss_broadion_c;  break;
        case IONOOPT_SBAS:  sppionf=new gnss_sbasion_c;   break;
        case IONOOPT_TEC:   sppionf=new gnss_ionexion_c;  break;
        case IONOOPT_QZS:   sppionf=new gnss_qzssion_c;   break;
        case IONOOPT_LEX:   sppionf=new gnss_lexioncor_c; break;
        default:sppionf=new gnss_ioncorr_c;
    }
    /* SPP troposhere function */
    switch (opt->spptrop) {
        case TROPOPT_SAAS: spptrof=new gnss_saastro_c; break;
        case TROPOPT_SBAS: spptrof=new gnss_sbastro_c; break;
        case TROPOPT_ZTD:  spptrof=new gnss_ztdtro_c; break;
        default: spptrof=new gnss_trocorr_c;
    }
    spptrof->nav=nav;
    /* ionosphere function */
    switch (opt->ionoopt) {
        case IONOOPT_IFLC:  ionfunc=new gnss_LCion_c;     break;
        case IONOOPT_BRDC:  ionfunc=new gnss_broadion_c;  break;
        case IONOOPT_SBAS:  ionfunc=new gnss_sbasion_c;   break;
        case IONOOPT_TEC:   ionfunc=new gnss_ionexion_c;  break;
        case IONOOPT_QZS:   ionfunc=new gnss_qzssion_c;   break;
        case IONOOPT_LEX:   ionfunc=new gnss_lexioncor_c; break;
        case IONOOPT_CONST: ionfunc=new gnss_estion_c;  break;
        default:ionfunc=new gnss_ioncorr_c;
    }
    /* troposphere function */
    switch (opt->tropopt) {
        case TROPOPT_SAAS: trofunc=new gnss_saastro_c; break;
        case TROPOPT_SBAS: trofunc=new gnss_sbastro_c; break;
        case TROPOPT_EST:
        case TROPOPT_ESTG: trofunc=new gnss_esttro_c; break;
        case TROPOPT_ZTD:  trofunc=new gnss_ztdtro_c; break;
        default: trofunc=new gnss_trocorr_c;
    }
    trofunc->nav=nav;

    if ( opt->mode==PMODE_SINGLE || opt->mode==PMODE_DGPS ) usePhase=0;
    else usePhase=1;
    /* ambiguity parameters */
    nfix=neb=0;
    /* satellite status vectors [MAXSAT] */
    for (int i=0; i<MAXSAT; i++) {
        ssat[i].init_vector(opt);
        ssat[i].sat=i+1;
        satno2id(i+1,ssat[i].id);
        ssat[i].sys=satsys(i+1,NULL);
    }
    /* solution vectors [MAXSOLBUF] */
    sol.assign(MAXSOLBUF,gnss_sol_c(this));
    b_sol.assign(MAXSOLBUF,gnss_sol_c(this));
}
/* get receiver position and velocity from external ------------------------------- */
void gnss_pro_c::gnss_get_posvel(const double pos[3], const double vel[3],
                                 const double vpos[9],const double vvel[9]) {
    /* save external position and velocity */
    for (int i=0; i<3; i++) {
        ext_posvel[i  ]=pos[i];
        ext_posvel[i+3]=vel[i];
    }
    if ( vpos!=NULL ) {
        for (int i=0; i<3; i++) {
            ext_vpos[i*3+i]=vpos[i*3+i];
            for (int j=0; j<i; j++) {
                ext_vpos[i*3+j]=ext_vpos[j*3+i]=vpos[i*3+j];
            }
        }
        varFlag|=1;
    }
    if ( vvel!=NULL ) {
        for (int i=0; i<3; i++) {
            ext_vvel[i*3+i]=vvel[i*3+i];
            for (int j=0; j<i; j++) {
                ext_vvel[i*3+j]=ext_vvel[j*3+i]=vvel[i*3+j];
            }
        }
        varFlag|=2;
    }
}
/* preprocess gnss observation ---------------------------------------------------- */
int gnss_pro_c::preprocess_gnss_obs(const double pos[3], const double vel[3]) {
    return 0;
}
/* get gnss observation equation (virtual) ---------------------------------------- */
int gnss_pro_c::gnss_obs_equation() {
    return 0;
}
/* update gnss status (virtual) --------------------------------------------------- */
int gnss_pro_c::gnss_update_status() {
    return 0;
}
/* gnss position function (virtual) ----------------------------------------------- */
int gnss_pro_c::basepos() {
    return 0;
}
/* gnss velocity function --------------------------------------------------------- */
int gnss_pro_c::gnss_vel() {
    return 0;
}
/* gnss position function --------------------------------------------------------- */
int gnss_pro_c::gnss_pos() {
    return 0;
}

/* GNSS solution class -------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
gnss_sol_c::gnss_sol_c() {
    nsol=0;
    NF=NL=NI=NT=NA=NG=0;
    NP=NV=3;
    NC=NSYS;
    time=obsTime=soltime=gtime_c();

    xpos.assign(3,0.0); vpos.assign(3*3,0.0);
    xvel.assign(3,0.0); vvel.assign(3*3,0.0);

    for (int i=0; i<6; i++) posvel[i]=0.0;
    for (int i=0; i<9; i++) posvar[i]=velvar[i]=ecefPvar[i]=ecefVvar[i]=0.0;
    for (int i=0; i<3; i++) refpos[i]=lat[i]=lon[i]=0.0;
    age=ratio=thres=0.0;
    type=stat=0;
    ns=0;
}
gnss_sol_c::gnss_sol_c(const gnss_pro_c *gnss_pro) {
    NF=gnss_pro->NF;
    NP=gnss_pro->NP;
    NV=gnss_pro->NV;
    NT=gnss_pro->NT;
    NG=gnss_pro->NG;
    NC=NSYS;

    xpos.assign(3,0.0); vpos.assign(3*3,0.0);
    xvel.assign(3,0.0); vvel.assign(3*3,0.0);
    if (NT) { xtro.assign(NT,0.0); vtro.assign(NT*NT,0.0); }
    if (NG) { xglo.assign(NG,0.0); vglo.assign(NG*NG,0.0); }
    xclk.assign(NSYS,0.0); vclk.assign(NSYS*NSYS,0.0);

    NL=NA=0;

    time=obsTime=soltime=gtime_c();
    for (int i=0; i<6; i++) posvel[i]=0.0;
    for (int i=0; i<9; i++) posvar[i]=velvar[i]=ecefPvar[i]=ecefVvar[i]=0.0;
    for (int i=0; i<3; i++) lat[i]=lon[i]=0.0;
    age=ratio=thres=0.0;
    type=stat='\0';
    ns=0;
}
gnss_sol_c::~gnss_sol_c() {
    xpos.clear(); xvel.clear(); xion.clear(); xtro.clear();
    xamb.clear(); xglo.clear(); xclk.clear();

    vpos.clear(); vvel.clear(); vion.clear(); vtro.clear();
    vamb.clear(); vglo.clear(); vclk.clear();
}
/* Implementation functions ------------------------------------------------------- */
/* dynamic covariance to ecef covariance ------------------------------------------ */
void gnss_sol_c::dyc2ecef(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro) {
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            ecefPvar[j+i*3]=vpos[j+i*3];
    if (gnss_pro->opt->dynamics) {
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                ecefVvar[j+i*3]=vvel[j+i*3];
    }
}
/* ecef solution ------------------------------------------------------------------ */
void gnss_sol_c::ecef(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro) {

    string str;
    strpv="";
    /* position */
    for (int i=0; i<3; i++) strpv+=doul2str(14,4," ",xpos[i],str)+opt->sep;
    /* solution state */
    strpv+=int2str(3," ",stat,str)+opt->sep+int2str(5," ",nsol,str)+opt->sep+
        int2str(3," ",ns,str)+opt->sep;
    /* position variance */
    /* xx yy zz */
    for (int i=0; i<3; i++)
        strpv+=doul2str(8,4," ",SQRT(vpos[i+i*3]),str)+opt->sep;
    /* xy yz zx */
    strpv+=doul2str(8,4," ",sqvar(vpos[1]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(vpos[5]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(vpos[2]),str);

    /* dynamical parameter */
    if (gnss_pro->opt->dynamics) {
        /* velocity */
        for (int i=0; i<3; i++) strpv+=doul2str(12,4," ",xvel[i],str)+opt->sep;
        /* velocity variance */
        /* vxx vyy vzz */
        for (int i=0; i<3; i++)
            strpv+=doul2str(10,4," ",SQRT(vvel[i+i*3]),str)+opt->sep;
        /* vxy vyz vzx */
        strpv+=doul2str(10,4," ",sqvar(vvel[1]),str)+opt->sep;
        strpv+=doul2str(10,4," ",sqvar(vvel[5]),str)+opt->sep;
        strpv+=doul2str(10,4," ",sqvar(vvel[2]),str);
    }
}
/* ecef position to LBH ----------------------------------------------------------- */
void gnss_sol_c::blh(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro) {
    /* compute blh and covariance and save to posvel and posvar */
    ecef2blh(xpos.begin(),opt->datum,posvel);
    dyc2ecef(opt,gnss_pro);
    covenu(posvel,ecefPvar,posvar);
    posvel[0]=posvel[0]*R2D; posvel[1]=posvel[1]*R2D;

    /* geodetic height */
    /*if (opt->height==1) posvel[2]-=;*/

    string str;
    strpv="";
    /* latitude and longitude ddd.ddd or ddd mm ss */
    if (opt->degf) {
        deg2dms(posvel[0],lat,5); //latitude
        deg2dms(posvel[1],lon,5); //longitude
        //longitude
        strpv+=doul2str(4,0," ",lon[0],str)+opt->sep+
            doul2str(2,0,"0",lon[1],str)+opt->sep+
            doul2str(8,5,"0",lon[2],str)+opt->sep;
        //latitude
        strpv+=doul2str(4,0," ",lat[0],str)+opt->sep+
            doul2str(2,0,"0",lat[1],str)+opt->sep+
            doul2str(8,5,"0",lat[2],str)+opt->sep;
    }
    else {
        //longitude latitude
        strpv+=doul2str(14,9," ",posvel[1],str)+opt->sep+
            doul2str(14,9," ",posvel[0],str)+opt->sep;
    }
    /* height */
    strpv+=doul2str(12,4," ",posvel[2],str)+opt->sep;
    /* solution state */
    strpv+=int2str(3," ",stat,str)+opt->sep+int2str(5," ",nsol,str)+opt->sep+
        int2str(3," ",ns,str)+opt->sep;
    /* position covariance */
    strpv+=doul2str(8,4," ",sqvar(posvar[0]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(posvar[4]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(posvar[8]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(posvar[1]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(posvar[5]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(posvar[2]),str);

    /* dynamical parameter */
    if (gnss_pro->opt->dynamics) {
        covenu(posvel,ecefVvar,velvar);
        ecef2enu(posvel,xvel,posvel+3);
        /* velocity */
        for (int i=0; i<3; i++) strpv+=doul2str(12,4," ",posvel[i+3],str)+opt->sep;
        /* velocity variance */
        /* vxx vyy vzz */
        for (int i=0; i<3; i++)
            strpv+=doul2str(10,4," ",SQRT(velvar[i+i*3]),str)+opt->sep;
        /* vxy vyz vzx */
        strpv+=doul2str(10,4," ",sqvar(velvar[1]),str)+opt->sep;
        strpv+=doul2str(10,4," ",sqvar(velvar[5]),str)+opt->sep;
        strpv+=doul2str(10,4," ",sqvar(velvar[2]),str);
    }
}
/* ecef position to ENU ----------------------------------------------------------- */
void gnss_sol_c::enu(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro) {
    /* compute rover enu and covariance and save to posvel and posvar */
    /* pos represents base blh */
    double pos[3],rr[3];
    for (int i=0; i<3; i++) rr[i]=xpos[i]-refpos[i];
    ecef2blh(refpos,opt->datum,pos);
    dyc2ecef(opt,gnss_pro);
    covenu(pos,ecefPvar,posvar);
    ecef2enu(pos,rr,posvel);

    /* write formated solution to strpv */
    string str;
    strpv="";
    /* position */
    for (int i=0; i<3; i++) strpv+=doul2str(14,4," ",posvel[i],str)+opt->sep;
    /* solution state */
    strpv+=int2str(3," ",stat,str)+opt->sep+int2str(5," ",nsol,str)+opt->sep+
        int2str(3," ",ns,str)+opt->sep;
    /* position covariance */
    for (int i=0; i<3; i++) strpv+=doul2str(8,4," ",SQRT(posvar[4*i]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(posvar[1]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(posvar[5]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(posvar[2]),str);

    /* dynamical parameter */
    if (gnss_pro->opt->dynamics) {
        covenu(pos,ecefVvar,velvar);
        ecef2enu(pos,xvel,posvel+3);
        /* velocity */
        for (int i=0; i<3; i++) strpv+=doul2str(12,4," ",posvel[i+3],str)+opt->sep;
        /* velocity variance */
        /* vxx vyy vzz */
        for (int i=0; i<3; i++)
            strpv+=doul2str(10,4," ",SQRT(velvar[i+i*3]),str)+opt->sep;
        /* vxy vyz vzx */
        strpv+=doul2str(10,4," ",sqvar(velvar[1]),str)+opt->sep;
        strpv+=doul2str(10,4," ",sqvar(velvar[5]),str)+opt->sep;
        strpv+=doul2str(10,4," ",sqvar(velvar[2]),str);
    }

    /* attitude */
    if ( gnss_pro->opt->mode==PMODE_MOVEB  || 
        ( gnss_pro->opt->mode>=PMODE_DGPS && gnss_pro->opt->baseAtt ) ) {
        /* attitude */
        for (int i=0; i<3; i++) strpv+=doul2str(11,5," ",gnss_pro->att_n_b[i],str)+opt->sep;
        /* error */
        strpv+=doul2str(10,5," ",gnss_pro->att_err[0],str)+opt->sep;
        strpv+=doul2str(10,5," ",gnss_pro->att_err[1],str)+opt->sep;
        strpv+=doul2str(10,5," ",gnss_pro->att_err[2],str);
    }
}
/* ecef position to EMEA ---------------------------------------------------------- */
void gnss_sol_c::nmea(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro) {
    string str;
    strpv="";
    /* position */
    for (int i=0; i<3; i++) strpv+=doul2str(14,4," ",xpos[i],str)+opt->sep;
    /* solution state */
    strpv+=int2str(3," ",stat,str)+opt->sep+int2str(5," ",nsol,str)+opt->sep+
        int2str(3," ",ns,str)+opt->sep;
    /* position variance */
    /* xx yy zz */
    for (int i=0; i<3; i++)
        strpv+=doul2str(8,4," ",SQRT(vpos[i*3+i]),str)+opt->sep;
    /* xy yz zx */
    strpv+=doul2str(8,4," ",sqvar(vpos[1]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(vpos[5]),str)+opt->sep;
    strpv+=doul2str(8,4," ",sqvar(vpos[2]),str);
}
/* solution position -------------------------------------------------------------- */
string gnss_sol_c::forposvel(const gnss_solopt_c *opt,const gnss_pro_c *gnss_pro) {
    switch (opt->posf) {
        case SOLF_LLH:  blh(opt,gnss_pro);  break;
        case SOLF_ENU:  enu(opt,gnss_pro);  break;
        case SOLF_NMEA: nmea(opt,gnss_pro); break;
        case SOLF_XYZ:
        default:        ecef(opt,gnss_pro);
    }
    return strpv;
}
/* solution time ------------------------------------------------------------------ */
string gnss_sol_c::fortime(const gnss_solopt_c *opt) {
    int timeu=opt->timeu<0 ? 0 : (opt->timeu>12 ? 12 : opt->timeu);
    soltime = time;
    if (opt->times>=TIMES_UTC) soltime.gpst2utc();
    if (opt->times==TIMES_BDT) soltime.gpst2bdt();
    if (opt->timef==1) strtime=soltime.time2str(timeu);
    else {
        int week;
        double gpst=soltime.time2gpst(&week);
        if (86400*7-gpst<0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        string strb;
        strtime=int2str(4,"0",week,strb)+" "+
            doul2str(6+(timeu<=0?0:timeu+1),timeu," ",gpst,strb);
    }
    return strtime;
}

/* satellite status class ----------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
gnss_ssat_c::gnss_ssat_c() {
    for (int i=0; i<NFREQ; i++) {
        resP[i]=resL[i]=0;
        vsat[i]=snr[i]=0;
    }
    id = "";
    sat=sys=vs=0;
    azel[0]=azel[1]=0.0;
    nCon=minCon=0;
    /* ionosphere */
    ion_delay=ion_var=0;
    ion_flag=0;
    /* ambiguity */
    nfreq=2;
    amb_obsvb[0][0]=amb_obsvb[1][0]=amb_obsvb[2][0]=amb_obsvb[NFREQ][0]=0;
    amb_obsvb[0][1]=amb_obsvb[1][1]=amb_obsvb[2][1]=amb_obsvb[NFREQ][1]=0;
    amb_solvab[0]=amb_solvab[1]=amb_solvab[2]=amb_solvab[NFREQ]=0;
    amb_CP[0]=amb_CP[1]=amb_CP[2]=amb_CP[NFREQ]=0;
    amb[0]=amb[1]=amb[2]=amb[NFREQ]=0;
    amb_ave[0]=amb_ave[1]=amb_ave[2]=amb_ave[NFREQ]=0;
    ambvar[0]=ambvar[1]=ambvar[2]=ambvar[NFREQ]=0;
    init_var[0]=init_var[1]=init_var[2]=init_var[NFREQ]=100;
    phase_time[0]=phase_time[1]=phase_time[2]=phase_time[NFREQ]=gtime_c(iGPS);
    solv_time[0]=solv_time[1]=solv_time[2]=solv_time[NFREQ]=gtime_c(iGPS);
    phase_dt[0]=phase_dt[1]=phase_dt[2]=phase_dt[NFREQ]=0;
    lock_period[0]=lock_period[1]=lock_period[2]=lock_period[NFREQ]=0;
    solv_period[0]=solv_period[1]=solv_period[2]=solv_period[NFREQ]=0;
    lock_con[0]=lock_con[1]=lock_con[2]=lock_con[NFREQ]=0;
    cobs_con[0]=cobs_con[1]=cobs_con[2]=cobs_con[NFREQ]=0;
    solv_con[0]=solv_con[1]=solv_con[2]=solv_con[NFREQ]=0;
    ave_con[0]=ave_con[1]=ave_con[2]=ave_con[NFREQ]=0;
    fix[0]=fix[1]=fix[2]=fix[NFREQ]=0;
    reset[0]=reset[1]=reset[2]=reset[NFREQ]=0;
    phw[0]=phw[1]=0.0;
    for (int i=0; i<2; i++) {
        gf12[i]=gf15[i]=0.0;
        lc12[i]=pc12[i]=0.0;
        mw12[i]=mw15[i]=0.0;
    }
    for (int i=0; i<NFREQ; i++) {
        lambda[i]=0;
        slip_f[i]=0;
        half[i]=0;
        d_ave[i]=0.0;
        multiPath[i]=0.0;
    }
    MW12_con=MW15_con=0;
    d_ave[NFREQ]=d_ave[NFREQ+1]=0.0;
    multiPath[NFREQ]=0.0;
}
gnss_ssat_c::~gnss_ssat_c() {
    obstime.clear();
    for (int i=0; i<NFREQ; i++) {
        L[i].clear(); D[i].clear(); P[i].clear();
        Lm[i].clear(); slip[i].clear();
    }
    d_ion.clear();

}
/* Implementation functions ------------------------------------------------------- */
/* set coefficients matrix for estimate of gf-slip -----------------------------------
* argv  : int flag       phase combination type {0:L12+L15,1:L12,2:L15}
* --------------------------------------------------------------------------------- */
void gnss_ssat_c::coe_gfslip(const int flag,vector<double> &An,vector<double> &Rn) {
    /* L12+L15 */
    if (flag==0) {
        //An lam1,-lam2,0; lam1,0,0; 0,lam2,0; lam1,0,-lam5; 0,0,lam5  
        An.assign(5*3,0.0); Rn.assign(5*5,0.0);
        An[0]=An[3]=An[9]=lambda[0],An[1]=-lambda[1]; An[7]=lambda[1];
        An[10]=-lambda[2]; An[14]=lambda[2];
        Rn[0]=Rn[18]=0.01; Rn[6]=Rn[12]=Rn[24]=10.0;
        return;
    }
    else {
        //An lam1,-lam2; lam1,0; lam2,0
        An.assign(3*2,0.0); Rn.assign(3*3,0.0);
        An[0]=An[2]=lambda[0]; An[1]=-lambda[flag]; An[4]=lambda[flag];
        Rn[0]=0.01; Rn[4]=Rn[8]=10.0;
    }
}
/* detect cycle slip functions ---------------------------------------------------- */
/* detect by LLI ------------------------------------------------------------------ */
void gnss_ssat_c::detslip_LLI(const gnss_obsd_c *rov,const gnss_obsd_c *bas) {
    unsigned char LLr[NFREQ]={ 0,0,0 },LLb[NFREQ]={ 0,0,0 };
    for (int i=0; i<nfreq; i++) {
        LLr[i]=( rov? rov->LLI[i] : 0 )&3; LLb[i]=( bas? bas->LLI[i] : 0 )&3;
        half[i]=(LLr[i]&2)||(LLb[i]&2);
        /* slip: 6-7:rover,4-5:base,1:half,0:slip */
        slip[i].back()=(LLr[i]<<6)|(LLb[i]<<4)|(LLr[i]|LLb[i]);
        if (slip[i].back()&3) slip_f[i]|=1;
    }
}
/* detect by detect by Phase-Code value ------------------------------------------- */
void gnss_ssat_c::detslip_PC(const gnss_prcopt_c *opt) {
    /* sigma^2 of new amb_PC and average amb_PC */
    double s2=SQR(15.0),s2_ave,d_amb;
    for (int i=0; i<nfreq; i++) {
        if ( amb_obsvb[i][1]==OBSTYPE_CL && ave_con[i]>0 ) {
            s2_ave=s2/ave_con[i];
            d_amb=fabs(amb_CP[i]-amb_ave[i]);
            if (d_amb>sqrt(s2+s2_ave)) {
                slip[i].back()|=1;
                slip_f[i]|=2;
            }
        }
    }
}
/* detect by geometry-free combination -----------------------------------------------
* if slip checked  gf[n]-gf[n-1] > ion_change_rate * dt + phase_noise
* --------------------------------------------------------------------------------- */
void gnss_ssat_c::detslip_gf(const gnss_prcopt_c *opt) {
    /* if dtime > 3min (ion change to much) return */
    double DN12=0.0,DN15=0.0;

    /* L1L2 */
    if (gf12[0]!=0.0&&gf12[1]!=0.0) {
        DN12 = gf12[1]-gf12[0];
        if (fabs(DN12)>opt->cs_gf) {
            slip[0].back()|=1; slip[1].back()|=1;
            slip_f[0]|=4; slip_f[1]|=4;
        }
    }

    /* L1L5 */
    if (nfreq>2) {
        if (gf15[0]!=0.0&&gf15[1]!=0.0) {
            DN15 = gf15[1]-gf15[0];
            if (fabs(DN15)>opt->cs_gf) {
                slip[0].back()|=1; slip[2].back()|=1;
                slip_f[0]|=4; slip_f[2]|=4;
            }
        }
    }
}
/* detect by Melbourne-Wubbena ---------------------------------------------------- */
void gnss_ssat_c::detslip_MW(const gnss_prcopt_c *opt) {
    /* sigma^2 of new MW and average MW */
    double s2=SQR(opt->cs_mw),s2_ave;
    /* Melbourne-Wubbena ambiguity */
    if (mw12[1]!=0.0&&mw12[0]!=0.0) {    /* mw12 */
        s2_ave=s2/MW12_con;
        if (fabs(mw12[1]-mw12[0])>3.0*sqrt(s2+s2_ave)) {
            slip[0].back()|=1; slip[1].back()|=1;
            slip_f[0]|=8; slip_f[1]|=8;
        }
    }
    if (nfreq>2&&mw15[1]!=0.0&&mw15[0]!=0.0) {    /* mw15 */
        s2_ave=s2/MW15_con;
        if (fabs(mw15[1]-mw15[0])>3.0*sqrt(s2+s2_ave)) {
            slip[0].back()|=1; slip[2].back()|=1;
            slip_f[0]|=8; slip_f[2]|=8;
        }
    }

    return;
}
/* calculate new ambiguity and initialize statistc data --------------------------- */
void gnss_ssat_c::new_amb_statistic() {
    d_ion.erase(d_ion.begin()); d_ion.push_back(0.0); //length for frequency i (m)
    // consider ionosphere effect for ambiguity (ion delay from P1-P2)
    if (ion_flag>0) {
        if (amb_obsvb[0][1]&OBSTYPE_PR&amb_obsvb[1][1]) {
            d_ion.back()=2.0*(P[0].back()-P[1].back())/(SQR(lambda[1]/lambda[0])-1.0);
        }
        else if (amb_obsvb[0][1]&OBSTYPE_PR&amb_obsvb[2][1]) {
            d_ion.back()=2.0*(P[0].back()-P[2].back())/(SQR(lambda[2]/lambda[0])-1.0);
        }
        else if (amb_obsvb[1][1]&OBSTYPE_PR&amb_obsvb[2][1]) {
            d_ion.back()=2.0*(P[1].back()-P[2].back())/(SQR(lambda[2]/lambda[1])-1.0);
        }
    }
    // update uncombined ambiguity
    for (int i=0; i<nfreq; i++) {
        if (amb_obsvb[i][1]==OBSTYPE_CL) {
            if (ion_flag>0) amb_CP[i]=(Lm[i].back()-P[i].back()-SQR(lambda[i]/lambda[0])*d_ion.back())/lambda[i];
            else amb_CP[i]=L[i].back()-P[i].back()/lambda[i];
        }
    }
    // update LC ambiguity
    if (amb_obsvb[NFREQ][1]==OBSTYPE_CL) amb_CP[NFREQ]=lc12[1]-pc12[1];
}
/* public: ------------------------------------------------------------------------ */
/* initialize vectors with order of polynomial fitting ---------------------------- */
void gnss_ssat_c::init_vector(const gnss_prcopt_c *opt) {
    if (opt->usedF>2) nfreq=opt->usedF;
    else nfreq=2;

    obstime.assign(MAXMEASBUF,gtime_c());
    minCon=opt->minConsecutive;

    for (int i=0; i<2; i++) {
        gf12[i]=gf15[i]=0.0;
        lc12[i]=pc12[i]=0.0;
        mw12[i]=mw15[i]=0.0;
    }
    /*nl12[i]=0.0;*/
    for (int i=0; i<NFREQ; i++) {
        L[i].assign(MAXMEASBUF,0.0);
        D[i].assign(MAXMEASBUF,0.0);
        P[i].assign(MAXMEASBUF,0.0);
        Lm[i].assign(MAXMEASBUF,0.0);
        slip[i].assign(MAXMEASBUF,0);
        amb_CP[i]=0.0;
        init_var[i]=SQR(opt->std[0]*5);
    }
    init_var[NFREQ]=SQR(opt->std[0]);
    if (opt->ionoopt!=IONOOPT_OFF) ion_flag=1;
    else ion_flag=0;
    d_ion.assign(MAXMEASBUF,0.0);
}
/* update carrier wave lengths ------------------------------------------------ */
void gnss_ssat_c::update_lambda(const gnss_nav_c &nav) {
    for (int i=0; i<NFREQ; i++) {
        lambda[i]=nav.lam[sat-1][i];
    }
}
/* update observations for current epoch ------------------------------------------ */
void gnss_ssat_c::input_new_obs(const gnss_obsd_c *rov,const gnss_obsd_c *bas,
    const double CJ[NFREQ]) {
    /* update obstime */
    obstime.erase(obstime.begin()); obstime.push_back(rov->time);

    /* observation */
    if (bas) {
        for (int i=0; i<nfreq; i++) {
            P[i].erase(P[i].begin()); P[i].push_back(single_diff(rov->P+i,bas->P+i)); //P
            L[i].erase(L[i].begin()); L[i].push_back(single_diff(rov->L+i,bas->L+i)); //L
            D[i].erase(D[i].begin()); D[i].push_back(single_diff(rov->D+i,bas->D+i)); //D
            if (L[i].back()!=0.0 && CJ!=NULL) L[i].back()+=CJ[i]/lambda[i]; //clock jump correction
        }
    }
    else {
        for (int i=0; i<nfreq; i++) {
            P[i].erase(P[i].begin()); P[i].push_back(rov->P[i]); //P
            L[i].erase(L[i].begin()); L[i].push_back(rov->L[i]); //L
            D[i].erase(D[i].begin()); D[i].push_back(rov->D[i]); //D
            if (L[i].back()!=0.0 && CJ!=NULL) L[i].back()+=CJ[i]/lambda[i]; //clock jump correction
        }
    }
    // uncombined ambiguity
    for (int i=0; i<nfreq; i++) {
        reset[i]=0;
        //initialize cycle slip
        slip_f[i]=0;
        slip[i].erase(slip[i].begin()); slip[i].push_back(0);
        // update observability
        amb_obsvb[i][0]=amb_obsvb[i][1];
        amb_obsvb[i][1]=0;
        if (P[i].back()!=0.0) amb_obsvb[i][1]|=OBSTYPE_PR;
        if (L[i].back()!=0.0) {
            amb_obsvb[i][1]|=OBSTYPE_CP;
            phase_dt[i]=(rov->time.timediff(phase_time[i]));
            phase_time[i]=rov->time;
        }
    }

    // LC ambiguity
    reset[NFREQ]=0;
    // update observability
    amb_obsvb[NFREQ][0]=amb_obsvb[NFREQ][1];
    amb_obsvb[NFREQ][1]=0;
    if (amb_obsvb[0][1]&amb_obsvb[1][1]&OBSTYPE_PR) amb_obsvb[NFREQ][1]|=OBSTYPE_PR;
    if (amb_obsvb[0][1]&amb_obsvb[1][1]&OBSTYPE_CP) {
        amb_obsvb[NFREQ][1]|=OBSTYPE_CP;
        phase_dt[NFREQ]=(rov->time.timediff(phase_time[NFREQ]));
        phase_time[NFREQ]=rov->time;
    }
}
/* update process parameter for current epoch ------------------------------------- */
void gnss_ssat_c::update_pro_par(const double CJ[NFREQ]) {
    /* update phase observation with clock jump */
    for (int i=0; i<nfreq; i++) {
        Lm[i].erase(Lm[i].begin());
        if (L[i].back()!=0.0) {
            if (CJ!=NULL) L[i].back()+=CJ[i]/lambda[i]; //clock jump correction
            Lm[i].push_back( L[i].back()*lambda[i] ); //Lm
        }
        else Lm[i].push_back( 0.0 ); //Lm
    }

    /* gf12, gf15 */
    if (gf12[1]!=0.0) gf12[0]=gf12[1]; gf12[1]=geometry_free(Lm[0].back(),Lm[1].back());
    if (gf15[1]!=0.0) gf15[0]=gf15[1]; gf15[1]=geometry_free(Lm[0].back(),Lm[2].back());

    /* pc12, lc12 */
    pc12[0]=pc12[1]; pc12[1]=iono_free(P[0].back(),P[1].back(),lambda[0],lambda[1]);
    lc12[0]=lc12[1]; lc12[1]=iono_free(Lm[0].back(),Lm[1].back(),lambda[0],lambda[1]);

    /* mw */
    mw12[1]=Mel_Wub(P[0].back(),P[1].back(),Lm[0].back(),Lm[1].back(),lambda[0],lambda[1]); 
    mw15[1]=Mel_Wub(P[0].back(),P[2].back(),Lm[0].back(),Lm[2].back(),lambda[0],lambda[2]);

    /* narrow-lane */
    /*nl12[1]=Narrow(P[0].back(),P[1].back(),Lm[0].back(),Lm[1].back(),lambda[0],lambda[1]);*/

    /* calculate new uncombined ambiguity */
    new_amb_statistic();
}
/* reset flag according to unsolved time interval and last ambiguity solution ----- */
void gnss_ssat_c::test_reset(const gnss_prcopt_c *opt) {
    /* reset uncombined ambiguity if no amb or out of time */
    for (int i=0; i<nfreq; i++) {
        if ( opt->modear==ARMODE_INST || ave_con[i]<1 || phase_dt[i]>opt->restime+DTTOL ) {
            fix[i]=1;
            reset[i]=1;
            amb_ave[i]=0;
            lock_period[i]=0;
            solv_period[i]=0;
            phase_dt[i]=0;
            lock_con[i]=0;
            cobs_con[i]=0;
            solv_con[i]=0;
            ave_con[i]=0;
        }
    }
}
/* reset LC ambiguity ------------------------------------------------------------- */
void gnss_ssat_c::test_reset_LC() {
    /* LC ambiguity */
    if (reset[0]>0||reset[1]>0) {
        fix[NFREQ]=1;
        reset[NFREQ]=1;
        amb_ave[NFREQ]=0;
        lock_period[NFREQ]=0;
        solv_period[NFREQ]=0;
        phase_dt[NFREQ]=0;
        lock_con[NFREQ]=0;
        cobs_con[NFREQ]=0;
        solv_con[NFREQ]=0;
        ave_con[NFREQ]=0;
    }
}
/* reset_ambiguity ---------------------------------------------------------------- */
void gnss_ssat_c::reset_amb(const int freq) {
    //statistic data
    fix[freq]=1;
    reset[freq]=1;
    amb_ave[freq]=0;
    solv_period[freq]=0;
    cobs_con[freq]=0;
    solv_con[freq]=0;
    ave_con[freq]=0;
}
/* detect cycle slip -------------------------------------------------------------- */
void gnss_ssat_c::detect_slip(const gnss_obsd_c *rov,const gnss_obsd_c *bas,
    const gnss_prcopt_c *opt) {
    /* detect by LLI */
    detslip_LLI(rov,bas);
    /* detect by detect by Phase-Code value */
    if (opt->slipmode&CSMODE_CP) detslip_PC(opt);
    /* detect by gf */
    if (opt->slipmode&CSMODE_GF) detslip_gf(opt);
    /* detect by TurboEdit (not used now) */
    if (opt->slipmode&CSMODE_MW) detslip_MW(opt);
    for (int f=0; f<nfreq; f++) if (slip[f].back()&1||half[f]) {
        fix[f]=1; reset[f]=1; amb_ave[f]=0;
        lock_period[f]=0; solv_period[f]=0; phase_dt[f]=0;
        lock_con[f]=0; cobs_con[f]=0; solv_con[f]=0; ave_con[f]=0;
    }
}
/* update ambiguity statistic ----------------------------------------------------- */
void gnss_ssat_c::update_amb_statistic(const gnss_prcopt_c *opt) {
    /* Melbourne-Wubbena ambiguity */
    if (mw12[1]!=0.0) { /* mw12 */
        if (reset[0]||reset[1]) { d_ave[NFREQ+1]=mw12[0]=mw12[1]; MW12_con=1; }
        else {
            d_ave[NFREQ+1]=(mw12[1]-mw12[0])/(++MW12_con);
            mw12[0] += d_ave[NFREQ+1];
        }
    }
    if (mw15[1]!=0.0) { /* mw15 */
        if (reset[0]||reset[2]) { mw15[0]=mw15[1]; MW15_con=1; }
        else mw15[0] += (mw15[1]-mw15[0])/(++MW15_con);
    }

    /* uncombined ambiguity */
    for (int i=0; i<nfreq; i++) {
        if (!(amb_obsvb[i][1]&OBSTYPE_CP)) continue;
        lock_con[i]++;
        cobs_con[i]++;
        lock_period[i]+=phase_dt[i];
        /* update ambiguity and variance */
        if ( amb_obsvb[i][1]&OBSTYPE_PR ) {
            /* average ambiguity */
            d_ave[i]=(amb_CP[i]-amb_ave[i])/(++ave_con[i]);
            amb_ave[i]+=d_ave[i]; //uncombined ambiguity
        }
        if (fix[i]==1) {
            //if (fix[i]==1) ambvar[i]=init_var[i]/(ave_con[i]);
            if (fix[i]==1) ambvar[i]=init_var[i];
        }
        else if (fix[i]==2) ambvar[i]+=opt->stdrate[0]*phase_dt[i];
        else if (fix[i]==3) ambvar[i]+=opt->stdrate[0]*phase_dt[i];
    }

    /* LC ambiguity */
    if (!(amb_obsvb[NFREQ][1]&OBSTYPE_CP)) return;
    lock_con[NFREQ]++;
    cobs_con[NFREQ]++;
    lock_period[NFREQ]+=phase_dt[NFREQ];
    /* update ambiguity and variance */
    if ( amb_obsvb[NFREQ][1]&OBSTYPE_PR ) {
        /* average ambiguity */
        d_ave[NFREQ]=(amb_CP[NFREQ]-amb_ave[NFREQ])/(++ave_con[NFREQ]);
        amb_ave[NFREQ]+=d_ave[NFREQ]; //LC ambiguity
    }
    if (fix[NFREQ]==1) {
        //ambvar[NFREQ]=init_var[NFREQ]/(ave_con[NFREQ]);
        ambvar[NFREQ]=init_var[NFREQ];
    }
    else if (fix[NFREQ]==2) ambvar[NFREQ]+=opt->stdrate[0]*phase_dt[NFREQ];
}
/* correct antenna and phase windup for observations and ambiguities -------------- */
void gnss_ssat_c::correct_obs_amb(const int freq,const gnss_obsd_c *rov,const gnss_obsd_c *bas) {

    if (freq<0||freq>NFREQ) return;
    double dant=0,dphw=0;
    if (bas) {
        dant = rov->rant[freq] + rov->sant[freq] - bas->rant[freq] - bas->sant[freq];
        dphw = phw[0]-phw[1];
    }
    else {
        dant = rov->rant[freq] + rov->sant[freq];
        dphw = phw[0];
    }
    // LC phase
    if (freq==NFREQ) {
        double gamma=SQR(lambda[1])/SQR(lambda[0]);
        if (amb_obsvb[NFREQ][1]&OBSTYPE_PR) pc12[1] += dant; //Code
        if (amb_obsvb[NFREQ][1]&OBSTYPE_CP) lc12[1] += dant - (gamma*dphw*lambda[0]-dphw*lambda[1])/(gamma-1.0); //Phase
        if (amb_obsvb[NFREQ][1]==OBSTYPE_CL) amb_CP[NFREQ] = lc12[1]-pc12[1]; //ambiguity
    }
    // uncombined phase
    else {
        if (amb_obsvb[freq][1]&OBSTYPE_PR) P[freq].back() += dant; //Code
        if (amb_obsvb[freq][1]&OBSTYPE_CP) L[freq].back() += dant/lambda[freq] - dphw; //Phase
        if (amb_obsvb[freq][1]==OBSTYPE_CL) amb_CP[freq] -= dphw; //ambiguity
    }
}
/* update ambiguity parameters ---------------------------------------------------- */
int gnss_ssat_c::update_solved_amb(const int freq,const int nsol) {
    if ( freq<0 || freq>NFREQ ) return 0;

    amb_solvab[freq]=0;

    /* update ambiguity solvability */
    /*if ( amb_obsvb[freq][1]&OBSTYPE_CP && ( solv_con[freq]>0 || amb_obsvb[freq][1]&OBSTYPE_PR ) &&
        ( nsol+2<2*minCon || lock_con[freq]>=minCon ) ) {*/
    if ( amb_obsvb[freq][1]&OBSTYPE_CP && ( solv_con[freq]>0 || amb_obsvb[freq][1]&OBSTYPE_PR ) ) {
        amb_solvab[freq]=1;
        /* update ambiguity and variance */
        if (fix[freq]==1) {
            /* if initialization use average ambiguity */
            //amb[freq]=amb_ave[freq];
            amb[freq]=amb_CP[freq];
        }
    }

    return amb_solvab[freq];
}

/* RTK server class ------------------------------------------------------------------------------- */
/* Constructor -------------------------------------------------------------------- */
gnss_rtksvr_c::gnss_rtksvr_c() {
    int i,j;
    /* num */
    tick=state=sampling=cyctime=nmeacycle=nmeareq=buffsize=navsel=nsbs=nsol=
        cputime=prcout=nave=nsb[0]=nsb[1]=0;
    for (i=0; i<3; i++) {
        rb_ave[i]=nmeapos[i]=nb[i]=npb[i]=fobs[i]=0;
        for (j=0; j<10; j++) nmsg[i][j]=0;
    }
    for (i=0; i<MAXSTRRTK; i++) strtype[i]=0;

    /* classes */
    solopt[0]=solopt[1]=gnss_solopt_c();
    nav=new gnss_nav_c;
    nav->eph.assign(MAXSAT*2,gnss_eph_c());
    nav->geph.assign(NSATGLO*2,gnss_geph_c());
    nav->seph.assign(NSATSBS*2,gnss_seph_c());
    nav->n=MAXSAT*2; nav->ng=NSATGLO*2; nav->ns=NSATSBS*2;

    obs[0]=obs[1]=obs[2]=gnss_obs_c();

    for (i=0; i<MAXSTRRTK; i++) { stream[i]=NULL; }
    for (i=0; i<3; i++) { format[i]=0; buff[i]=pbuf[i]=NULL; data[i]=NULL; }
    sbuf[0]=sbuf[1]=NULL;
    gnss_pro=NULL;
    moni=NULL;

    for (i=0; i<MAXSBSMSG; i++) sbsmsg[i]=gnss_sbsmsg_c();

    initlock(&lock_f);
}
gnss_rtksvr_c::~gnss_rtksvr_c() {
    if (nav) { delete nav; nav=NULL; }
    for (int i=0; i<3; i++) {
        if (buff[i]) { delete[] buff[i]; buff[i]=NULL; }
        if (pbuf[i]) { delete[] pbuf[i]; pbuf[i]=NULL; }
        if (data[i]) { delete data[i]; data[i]=NULL; }
    }
    for (int i=0; i<MAXSTRRTK; i++)
        if (stream[i]) { delete stream[i]; stream[i]=NULL; }
    if (sbuf[0]) { delete[] sbuf[0]; sbuf[0]=NULL; }
    if (sbuf[1]) { delete[] sbuf[1]; sbuf[1]=NULL; }
    if (gnss_pro) { delete gnss_pro; gnss_pro=NULL; }
    moni=NULL;

    prcopt=NULL;
}
/* Implementation functions ------------------------------------------------------- */
/* initialzie function ------------------------------------------------------------ */
/* set sat\rcv antenna information ---------------------------------------------------
* argv   :  int   satrcv   0:satellite,1:receiver
* --------------------------------------------------------------------------------- */
void gnss_rtksvr_c::setpcv(in_gnss_atx_c *atx,int satrcv) {
    /* update satellite antenna to gnss_pro->nav */
    if (satrcv==0) {
        for (int i=0; i<MAXSAT; i++) {
            if (!(satsys(i+1,NULL)&gnss_pro->opt->navsys)) continue;
            for (int np=0; np<atx->satpcv.size(); np++) {
                if (atx->satpcv[np].sat!=i+1) continue;
                if (atx->satpcv[np].te.time>0) continue;
                nav->pcvs[i]=atx->satpcv[np];
                break;
            }
        }
    }
    /* update receiver antenna to gnss_pro->opt */
    else {
        gnss_pro->opt->pcvr[0]=atx->recpcv[0];
        gnss_pro->opt->pcvr[1]=atx->recpcv[1];
    }
}
/* new gnss_pro according to opt -------------------------------------------------- */
void gnss_rtksvr_c::ini_gnss_pro(gnss_prcopt_c *Prcopt,gnss_filopt_t *Filopt) {
    if (gnss_pro) { delete gnss_pro; gnss_pro=NULL; }
    if (Prcopt->mode==PMODE_SINGLE) gnss_pro=new gnss_single_c;
    else if (Prcopt->mode>=PMODE_DGPS) gnss_pro=new gnss_relative_c;
    else  gnss_pro=new gnss_ppp_c;
    gnss_pro->opt=Prcopt;
    /* set base station position */
    for (int i=0; i<6; i++) {
        gnss_pro->rb[i]=i<3 ? Prcopt->rb[i] : 0.0;
    }

    /* set frequency priority according to configuration */
    nav->freq.change_freq_priority(Prcopt);

    /* initialize navigation pointer */
    gnss_pro->nav=nav;
    /* read erp file */
    if (Filopt->erp.size()>10) {
        in_gnss_erp_c erp(Filopt->erp);
        erp.readerp(&nav->erp);
        gnss_pro->tidefunc.init_erp(nav);
        gnss_pro->satantfunc.init_erp(nav);
    }
    /* read blq file */
    if (Filopt->blq.size()>10&&Prcopt->tidecorr&2) {
        in_gnss_blq_c blq(Filopt->blq);
        blq.readblq(Prcopt->name[0],nav->ocean_par[0]);
        if (Prcopt->mode>=PMODE_DGPS)
            blq.readblq(Prcopt->name[1],nav->ocean_par[1]);
        /* tidal displacement functions */
        gnss_pro->tidefunc.init_otl(gnss_pro->opt,nav);
    }
    /* read sat\rec antenna information file */
    if (Filopt->satantp.size()>10) {
        in_gnss_atx_c sat(Filopt->satantp);
        sat.read_satatx(nav);
        setpcv(&sat,0);
    }
    if (Filopt->rcvantp.size()>10) {
        in_gnss_atx_c rcv(Filopt->rcvantp);
        rcv.read_recatx(gnss_pro->opt,nav);
        setpcv(&rcv,1);
    }
    /* open test file */
    if (Filopt->test.length()>5) gnss_pro->log_stream.open(Filopt->test,ios::out);
}
/* initialize stream environment -------------------------------------------------- */
void gnss_rtksvr_c::strinitcom() {
#ifdef WIN32
    WSADATA data;
    WSAStartup(MAKEWORD(2,0),&data);
#endif
}
/* initialize stream class -------------------------------------------------------- */
void gnss_rtksvr_c::inistream() {

    for (int i=0; i<MAXSTRRTK; i++) {
        if (stream[i]) { delete stream[i]; stream[i]=NULL; }
        switch (strtype[i]) {
        case STR_SERIAL:   stream[i]=new serial_c; break;
        case STR_FILE:     stream[i]=new file_c;   break;
        case STR_TCPSVR:   stream[i]=new tcpsvr_c; break;
        case STR_TCPCLI:   stream[i]=new tcpcli_c; break;
        case STR_NTRIPSVR:
        case STR_NTRIPCLI: stream[i]=new ntrip_c;  break;
        case STR_FTP:
        case STR_HTTP:     stream[i]=new ftp_c;    break;
        case STR_NTRIPC_S:
        case STR_NTRIPC_C: stream[i]=new ntripc_c; break;
        case STR_UDPSVR:
        case STR_UDPCLI:   stream[i]=new udp_c;    break;
        case STR_MEMBUF:   stream[i]=new membuf_c; break;
        default: stream[i]=new stream_c;
        }
        stream[i]->Stype=strtype[i];
    }
}
/* initialize decode format ------------------------------------------------------- */
void gnss_rtksvr_c::inidecode() {

    for (int i=0; i<3; i++) {
        if (data[i]) { delete data[i]; data[i]=NULL; }
        switch (format[i]) {
        case STRFMT_RTCM2: data[i]=new rtcm2_c; break;
        case STRFMT_RTCM3: data[i]=new rtcm3_c; break;
        case STRFMT_OEM3:  data[i]=new oem3;   break;
        case STRFMT_OEM4:  data[i]=new oem4;   break;
        case STRFMT_UBX:   data[i]=new ublox;  break;
        case STRFMT_SS2:   data[i]=new ss2;    break;
        case STRFMT_CRES:  data[i]=new cres;   break;
        case STRFMT_STQ:   data[i]=new skyq;   break;
        case STRFMT_GW10:  data[i]=new gw10;   break;
        case STRFMT_JAVAD: data[i]=new javad;  break;
        case STRFMT_NVS:   data[i]=new nvs;    break;
        case STRFMT_BINEX: data[i]=new binex;  break;
        case STRFMT_RT17:  data[i]=new rt17;   break;
        case STRFMT_SEPT:  data[i]=new sbf;    break;
        case STRFMT_LEXR:  data[i]=new decode_data_c;       break;
        case STRFMT_CMR:   data[i]=new cmr; data[i]->Svr=this; break;
        default: data[i]=new decode_data_c;
        }
        data[i]->format=format[i];
    }
}
/* sync input streams (if type=STR_FILE) ------------------------------------------ */
void gnss_rtksvr_c::strsync() {
    if (stream[0]->Stype==STR_FILE&&stream[1]->Stype==STR_FILE)
        stream[0]->strsync(stream+1);
    if (stream[0]->Stype==STR_FILE&&stream[2]->Stype==STR_FILE)
        stream[0]->strsync(stream+2);
}
/* write solution header to output stream ----------------------------------------- */
void gnss_rtksvr_c::writesolhead() {
    unsigned char buff1[8192]={ 0 };
    unsigned char buff2[8192]={ 0 };
    int n;

    /* output 1 */
    if (solopt[0].outopt) {
        n=solopt[0].outsolheads(prcopt,buff1);
        stream[3]->StreamWrite(buff1,n);
    }
    /* output 2 */
    if (solopt[1].outopt) {
        n=solopt[1].outsolheads(prcopt,buff2);
        stream[4]->StreamWrite(buff2,n);
    }
}
/* update glonass frequency channel number in raw data struct --------------------- */
void gnss_rtksvr_c::updatefcn() {
    int i,j,sat,frq;

    for (i=0; i<MAXPRNGLO; i++) {
        sat=satno(SYS_GLO,i+1);

        for (j=0,frq=-999; j<3; j++) {
            if (data[j]->nav.geph[i].sat!=sat) continue;
            frq=data[j]->nav.geph[i].frq;
        }
        if (frq<-7||frq>6) continue;

        for (j=0; j<3; j++) {
            if (data[j]->nav.geph[i].sat==sat) continue;
            data[j]->nav.geph[i].sat=sat;
            data[j]->nav.geph[i].frq=frq;
        }
    }
}
/* write solution to each out-stream (stream[3:4])--------------------------------- */
void gnss_rtksvr_c::writesolstr(int index) {
    unsigned char buff[MAXSOLMSG+1]={ 0 };
    char *p=(char *)buff;

    /* write solution to buff */
    /* [1] write solution time */
    string soltime=gnss_pro->sol.back().fortime(solopt+index);
    p+=sprintf(p,"%s%s",soltime.c_str(),solopt[index].sep.c_str());

    /* [2] position solution */
    string strpv;
    if ( gnss_pro->sol.back().stat==SOLQ_NONE ) return;
    else if (solopt[index].posf==SOLF_ENU) {
        if ( solopt[index].refsta==REFSAT_BASE || prcopt->mode==PMODE_MOVEB ) {
            for (int i=0; i<3; i++) {
                gnss_pro->sol.back().refpos[i]=gnss_pro->rb[i];
            }
        }
        else {
            for (int i=0; i<3; i++) {
                gnss_pro->sol.back().refpos[i]=prcopt->ru[i];
            }
        }
        if ( norm(gnss_pro->sol.back().refpos,3)<1E5 ) return;
    }
    strpv=gnss_pro->sol.back().forposvel(solopt+index,gnss_pro);
    p+=sprintf(p,"%s\n",strpv.c_str());

    /* write solution to stream[index+3] */
    stream[index+3]->StreamWrite(buff,p-(char *)buff);
}

/* initialize observation pointer obsr/obsb (*gnss_pro) --------------------------- */
int gnss_rtksvr_c::iniobs() {
    gnss_pro->obsr=&obs[0];
    gnss_pro->obsb=&obs[1];

    /* test availability of rover observation */
    if (fobs[0]<=0) return 0;
    /* test availability of base observation and time synchronization */
    if ((gnss_pro->opt->mode>=PMODE_DGPS)&&
        (fobs[1]<=0||fabs(gnss_pro->obsr->data[0].time.timediff(gnss_pro->obsb->data[0].time))>=1E-3))
        return 0;
    if (obs[0].data[0].time.timediff(gnss_pro->sol.back().time)<sampling-1E-3) {
        errorMsg="duplicated observation!\n";
        return 0;
    }

    return fobs[0];
}

/* process function --------------------------------------------------------------- */
/* lock/unlock rtk server --------------------------------------------------------- */
void gnss_rtksvr_c::rtksvrlock() {
    tolock(&lock_f);
}
void gnss_rtksvr_c::rtksvrunlock() {
    tounlock(&lock_f);
}
/* write solution to output stream ------------------------------------------------ */
void gnss_rtksvr_c::writesol() {
    writesolstr(0);
    writesolstr(1);
}
/* input message from stream ------------------------------------------------------ */
/* update rtk server struct ------------------------------------------------------- */
void gnss_rtksvr_c::updatesvr(int ret,int index) {
    gnss_eph_c *eph1,*eph2,*eph3;
    gnss_geph_c *geph1,*geph2,*geph3;
    gtime_c tof;
    double pos[3],del[3]={ 0 },dr[3];
    int i,prn,sbssat=gnss_pro->opt->sbassatsel,sys,iode;

    /* observation data */
    if (ret==1) {
        /* initialize obs[index] */
        obs[index].reset();
        if (obs[index].n<MAXOBS) {
            for (i=0; i<data[index]->obs.n; i++) {
                data[index]->obs.data[i].sys=
                    satsys(data[index]->obs.data[i].sat,&data[index]->obs.data[i].prn);
                data[index]->obs.data[i].isys=syscd2num(data[index]->obs.data[i].sys);
                if (gnss_pro->opt->exsats[data[index]->obs.data[i].sat-1]==1||
                    !(data[index]->obs.data[i].sys&gnss_pro->opt->navsys))
                    continue;
                obs[index].data.push_back(data[index]->obs.data[i]);
                obs[index].data.back().rcv=index; //rev flag
            }
            obs[index].n=obs[index].data.size();
            /* arrange observation data */
            sortobs(obs[index]);
        }
        obs[index].rcv=index;
        nmsg[index][0]++;
    }
    /* ephemeris */
    else if (ret==2) {
        if (satsys(data[index]->ephsat,&prn)!=SYS_GLO) {
            if (!navsel||navsel==index+1) {
                eph1=&data[index]->nav.eph[data[index]->ephsat-1];
                eph2=&nav->eph[data[index]->ephsat-1];
                eph3=&nav->eph[data[index]->ephsat-1+MAXSAT];
                if (eph2->ttr.time==0||
                    (eph1->iode!=eph3->iode&&eph1->iode!=eph2->iode)||
                    (eph1->toe.timediff(eph3->toe)!=0.0&&
                        eph1->toe.timediff(eph2->toe)!=0.0)) {
                    *eph3=*eph2;
                    *eph2=*eph1;
                    nav->update_sat_lambda(gnss_pro->ssat);
                }
            }
            nmsg[index][1]++;
        }
        else {
            if (!navsel||navsel==index+1) {
                geph1=&data[index]->nav.geph[prn-1];
                geph2=&nav->geph[prn-1];
                geph3=&nav->geph[prn-1+MAXPRNGLO];
                if (geph2->tof.time==0||
                    (geph1->iode!=geph3->iode&&geph1->iode!=geph2->iode)) {
                    *geph3=*geph2;
                    *geph2=*geph1;
                    nav->update_sat_lambda(gnss_pro->ssat);
                    updatefcn();
                }
            }
            nmsg[index][6]++;
        }
    }
    /* sbas message */
    else if (ret==3) {
        if (sbssat==data[index]->sbsmsg.prn||sbssat==0) {
            if (nsbs<MAXSBSMSG) {
                sbsmsg[nsbs++]=data[index]->sbsmsg;
            }
            else {
                for (i=0; i<MAXSBSMSG-1; i++) sbsmsg[i]=sbsmsg[i+1];
                sbsmsg[i]=data[index]->sbsmsg;
            }
            data[index]->sbsmsg.sbsupdatecorr(nav);
        }
        nmsg[index][3]++;
    }
    /* ion/utc parameters */
    else if (ret==9) {
        if (navsel==0||navsel==index+1) {
            for (i=0; i<8; i++) nav->ion_gps[i]=data[index]->nav.ion_gps[i];
            for (i=0; i<4; i++) nav->utc_gps[i]=data[index]->nav.utc_gps[i];
            for (i=0; i<4; i++) nav->ion_gal[i]=data[index]->nav.ion_gal[i];
            for (i=0; i<4; i++) nav->utc_gal[i]=data[index]->nav.utc_gal[i];
            for (i=0; i<8; i++) nav->ion_qzs[i]=data[index]->nav.ion_qzs[i];
            for (i=0; i<4; i++) nav->utc_qzs[i]=data[index]->nav.utc_qzs[i];
            nav->leaps=data[index]->nav.leaps;
        }
        nmsg[index][2]++;
    }
    /* antenna postion parameters */
    else if (ret==5) {
        if (index==1 && (gnss_pro->opt->baspos==POSOPT_RTCM||gnss_pro->opt->baspos==POSOPT_RAW)) {
            for (i=0; i<3; i++) {
                gnss_pro->rb[i]=data[1]->sta.pos[i];
            }
            /* antenna delta */
            ecef2blh(gnss_pro->rb,WGS84,pos);
            if (data[1]->sta.deltype) { /* xyz */
                del[2]=data[1]->sta.hgt;
                enu2ecef(pos,del,dr);
                for (i=0; i<3; i++) {
                    gnss_pro->rb[i]+=data[1]->sta.del[i]+dr[i];
                }
            }
            else { /* enu */
                enu2ecef(pos,data[1]->sta.del,dr);
                for (i=0; i<3; i++) {
                    gnss_pro->rb[i]+=dr[i];
                }
            }
        }
        nmsg[index][4]++;
    }
    /* dgps correction */
    else if (ret==7) {
        nmsg[index][5]++;
    }
    /* ssr message */
    else if (ret==10) {
        for (i=0; i<MAXSAT; i++) {
            if (!data[index]->ssr[i].update) continue;

            /* check consistency between iods of orbit and clock */
            if (data[index]->ssr[i].iod[0]!=
                data[index]->ssr[i].iod[1]) continue;

            data[index]->ssr[i].update=0;

            iode=data[index]->ssr[i].iode;
            sys=satsys(i+1,&prn);

            /* check corresponding ephemeris exists */
            if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS) {
                if (nav->eph[i].iode!=iode&&
                    nav->eph[i+MAXSAT].iode!=iode) {
                    continue;
                }
            }
            else if (sys==SYS_GLO) {
                if (nav->geph[prn-1].iode!=iode&&
                    nav->geph[prn-1+MAXPRNGLO].iode!=iode) {
                    continue;
                }
            }
            nav->ssr[i]=data[index]->ssr[i];
        }
        nmsg[index][7]++;
    }
    /* lex message */
    else if (ret==31) {
        data[index]->lexmsg.lexupdatecorr(nav,tof);
        nmsg[index][8]++;
    }
    /* error */
    else if (ret==-1) {
        nmsg[index][9]++;
    }
}
/* decode receiver raw/rtcm data -------------------------------------------------- */
int gnss_rtksvr_c::decoderaw(int index) {
    int i,ret;

    /* initialize */
    rtksvrlock();

    for (i=0; i<nb[index]; i++) {

        /* input rtcm/receiver raw data from stream */
        ret=data[index]->decode((unsigned char)buff[index][i]);

        /* update rtk server */
        if (ret>0) updatesvr(ret,index);

        /* observation data received */
        if (ret==1) {
            if (obs[index].n<=MAXOBS) fobs[index]=obs[index].n;
            else { prcout++; fobs[index]=0; }
        }
    }
    nb[index]=0;

    rtksvrunlock();

    return fobs[index];
}

/* initialize rtksvr -------------------------------------------------------------- */
int gnss_rtksvr_c::rtksvrini(all_option_c *option) {
    gtime_c time;
    int i,j,rw;

    if (state) return 0;

    prcopt=&option->prcopt;
    strinitcom();
    sampling=option->prcopt.sampling;
    cyctime=option->rtkopt.svrcycle>1 ? option->rtkopt.svrcycle : 1;
    nmeacycle=option->rtkopt.nmeacycle>1000 ? option->rtkopt.nmeacycle : 1000;
    nmeareq=option->rtkopt.nmeareq;
    for (i=0; i<3; i++) nmeapos[i]=option->rtkopt.nmeapos[i];
    buffsize=option->rtkopt.buffsize>4096 ? option->rtkopt.buffsize : 4096;
    for (i=0; i<3; i++) format[i]=option->rtkopt.strfmt[i];
    for (i=0; i<7; i++) strtype[i]=option->rtkopt.strtype[i];
    navsel=option->rtkopt.navmsgsel;
    nsbs=0;
    nsol=0;
    prcout=0;

    /* initialize gnss_pro */
    ini_gnss_pro(&option->prcopt,&option->filopt);
    gnss_pro->gnss_pro_init();

    if (option->prcopt.initrst) {
        nave=0;
        rb_ave[0]=rb_ave[1]=rb_ave[2]=0.0;
    }

    /* initialize decode and stream format */
    inidecode();
    for (i=0; i<3; i++) { /* input/log streams */
        nb[i]=npb[i]=0;
        if (!(buff[i]=new unsigned char[buffsize])||
            !(pbuf[i]=new unsigned char[buffsize])) {
            return 0;
        }
        for (j=0; j<10; j++) nmsg[i][j]=0;

        /* set receiver and rtcm option */
        data[i]->opt=option->rtkopt.rropts[i];

        /* connect dgps corrections */
        data[i]->dgps=nav->dgps;
    }
    for (i=0; i<2; i++) { /* output peek buffer */
        if (!(sbuf[i]=new unsigned char[buffsize])) {
            return 0;
        }
    }
    /* set solution options */
    solopt[0]=option->solopt[0];
    solopt[1]=option->solopt[1];

    /* update navigation data */
    nav->update_sat_lambda(gnss_pro->ssat);

    /* set monitor stream */
    moni=option->rtkopt.monitor;

    /* initialize streams */
    inistream();
    /* open input streams */
    for (i=0; i<MAXSTRRTK; i++) {
        rw=i<3 ? STR_MODE_R : STR_MODE_W;
        if (option->rtkopt.strtype[i]!=STR_FILE) rw|=STR_MODE_W;
        if (!stream[i]->StreamOpen(option->rtkopt.strpath[i].c_str(),strtype[i],rw)) {
            for (i--; i>=0; i--) stream[i]->StreamClose();
            return 0;
        }
        /* set initial time for rtcm and raw */
        if (i<3) {
            time.timeget()->utc2gpst();
            data[i]->time=
                option->rtkopt.strtype[i]==STR_FILE ? stream[i]->strgettime() : time;
        }
    }
    /* sync input streams (if type=STR_FILE) */
    strsync();

    /* write start commands to input streams */
    for (i=0; i<3; i++) {
        if (option->rtkopt.cmds[i]!="\0") stream[i]->SendCmd(option->rtkopt.cmds[i].c_str());
    }
    /* write solution header to solution streams */
    writesolhead();

    //test
    cout <<"rtksvr initialization is ok!\n";

    return 1;
}

/* thread-start function ---------------------------------------------------------- */
#ifdef WIN32
static DWORD WINAPI rtksvrthread(void *arg)
#else
static void * rtksvrthread(void *arg)
#endif
{
    /* initailize arg to gnss_rtksvr_c */
    gnss_rtksvr_c *svr=(gnss_rtksvr_c *)arg;
    /* compute time and thread run time */
    double cpttime,runtime;
    /* thread-start time and last position-fall time */
    unsigned int startick,lastfall;
    /* position cycle */
    int cycle;
    /* solution time (utc) */
    gtime_c soltime;

    /* initialize svr */
    svr->state=1; svr->tick=tickget();
    lastfall=svr->tick-1000;

    for (cycle=0; svr->state; cycle++) {
        startick=tickget();

        for (int i=0; i<3; i++) {
            /* pointer to buff head and tail */
            unsigned char *bufhead=svr->buff[i]+svr->nb[i],
                *buftail=svr->buff[i]+svr->buffsize;
            int rbufn; //recevied buff number

            /* read receiver raw/rtcm data from input stream */
            if ((rbufn=svr->stream[i]->StreamRead(bufhead,buftail-bufhead))<=0) continue;

            /* write receiver raw/rtcm data to log stream */
            svr->stream[i+5]->StreamWrite(bufhead,rbufn);
            svr->nb[i]+=rbufn;

            /* save peek buffer */
            svr->rtksvrlock();
            rbufn=rbufn<svr->buffsize-svr->npb[i] ? rbufn : svr->buffsize-svr->npb[i];
            memcpy(svr->pbuf[i]+svr->npb[i],bufhead,rbufn);
            svr->npb[i]+=rbufn;
            svr->rtksvrunlock();
        }
        for (int i=0; i<3; i++) {
            /* decode receiver raw/rtcm data */
            svr->decoderaw(i);
        }
        if (svr->iniobs()) {

            /* SPP for base station */
            if ((svr->gnss_pro->opt->mode>=PMODE_DGPS)&&svr->fobs[1]>0) {
                if ((svr->gnss_pro->opt->maxaveep<=0||svr->nave<svr->gnss_pro->opt->maxaveep)&&
                    svr->gnss_pro->basepos()) { //return solution to b_sol
                    svr->nave++;
                    for (int i=0; i<3; i++)
                        svr->rb_ave[i]+=(svr->gnss_pro->b_sol.back().xpos[i]-svr->rb_ave[i])/svr->nave;
                }
                for (int i=0; i<3; i++) {
                    if (svr->gnss_pro->opt->baspos==POSOPT_SINGLE) svr->gnss_pro->rb[i]=svr->rb_ave[i];
                    else if (svr->gnss_pro->opt->mode==PMODE_MOVEB) {
                        svr->gnss_pro->rb[i]=svr->gnss_pro->b_sol.back().xpos[i];
                    }
                }
            }

            /* gnss positioning for rover */
            svr->rtksvrlock();
            svr->gnss_pro->gnss_pos();
            svr->rtksvrunlock();

            /* output solution if sol.stat */
            if (svr->gnss_pro->sol.back().stat!=SOLQ_NONE) {
                /* adjust difference between computer time and UTC time */
                cpttime=(int)(tickget()-startick)/1000.0+DTTOL;
                soltime=svr->gnss_pro->sol.back().time;
                soltime.timeadd(cpttime)->gpst2utc()->timeset();

                /* write solution */
                svr->writesol();
            }
            /* send null solution if no solution (1hz) */
            else if (svr->gnss_pro->sol.back().stat==SOLQ_NONE&&(int)(startick-lastfall)>=1000) {
                svr->writesol();
                lastfall=startick;
            }
        }

        if ((runtime=(int)(tickget()-startick))>0) svr->cputime=runtime;

        /* sleep until next cycle */
        sleepms(svr->cyctime-runtime);
    }
    /* close stream */
    for (int i=0; i<MAXSTRRTK; i++) svr->stream[i]->StreamClose();
#ifdef WIN32
    return 0;
#endif
}

/* start rtksvr ------------------------------------------------------------------- */
int gnss_rtksvr_c::rtksvrstart() {
#ifdef WIN32
    if (!(thread=CreateThread(NULL,0,rtksvrthread,this,0,NULL)))
#else
    if (pthread_create(&thread,NULL,rtksvrthread,this))
#endif
    {
        for (int i=0; i<MAXSTRRTK; i++) stream[i]->StreamClose();
        errorMsg="thread create error\n";
        return 0;
    }
    return 1;
}

/* stop rtksvr -------------------------------------------------------------------- */
void gnss_rtksvr_c::rtksvrstop(char **cmds) {
    /* write stop commands to input streams */
    rtksvrlock();
    for (int i=0; i<3; i++) {
        if (cmds[i]) stream[i]->SendCmd(cmds[i]);
    }
    rtksvrunlock();

    /* stop rtk server */
    state=0;

    /* free rtk server thread */
#ifdef WIN32
    WaitForSingleObject(thread,10000);
    CloseHandle(thread);
#else
    pthread_join(thread,NULL);
#endif
}
