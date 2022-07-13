#include "GNSS/AntModel/antenna.h"
#include "BaseFunction/basefunction.h"

/* Constant --------------------------------------------------------------------------------------- */
#define SQR(x)      ((x)*(x))

/* Functions -------------------------------------------------------------------------------------- */
/* interpolate antenna phase center variation ----------------------------------------
 * check party and decode navigation data word
 * args   : double   ang  I  zenith angle (receiver) or
 *                           nadir angle (satellite) in degree
 *             int   sys  I  GNSS system (0:G,1:R,2:E,3:C) or 0 for satellite
 *             int  freq  I  signal frequency index (0:L1,1:L2,2:L3)
 *           gnss_pcv_c   pcv  I  antenna PCV correction information class
 * return : antenna phase center variation value
 * notes  : see reference [1] 20.3.5.2 user parity algorithm
 * -------------------------------------------------------------------------------- */
static double interpvar(double ang,int sys,int freq,const gnss_pcv_c *pcv)
{
    if (pcv->var[sys][freq].size()<=0) return 0.0;
    double a=(ang-pcv->zen[0])/pcv->zen[2]; /* ang=zen0--zen1 , step=zen2 */
    int i=(int)a;
    if (i<0) return pcv->var[sys][freq][0];
    else if (i>=pcv->nzen-1) return pcv->var[sys][freq][pcv->nzen-1];
    return pcv->var[sys][freq][i]*(1.0-a+i)+pcv->var[sys][freq][i+1]*(a-i);
}

/* satellite antenna phase center offest ---------------------------------------------------------- */
/* Constuctors -------------------------------------------------------------------- */
gnss_satant_c::gnss_satant_c() {
    lam=NULL;
    er_par=NULL;
}
gnss_satant_c::~gnss_satant_c() {
    lam=NULL;
    er_par=NULL;
}
/* Implementation functions ------------------------------------------------------- */
/* nominal yaw-angle ------------------------------------------------------ */
double gnss_satant_c::yaw_nominal(double beta,double mu) {
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite ------------------------------------------------- */
int gnss_satant_c::yaw_angle(int sat,const char *type,int opt,double beta,double mu,
    double *yaw) {
    *yaw=yaw_nominal(beta,mu);
    return 1;
}
/* yaw-angle of satellite ------------------------------------------------- */
int gnss_satant_c::sat_yaw(gtime_c time,const double *rs,double *exs,double *eys) {
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={ 0 };
    int i;

    /* gps time to utc time */
    gtime_c utcTime = time;
    utcTime.gpst2utc();
    /* sun position in ecef */
    if (er_par) geterpv(er_par,utcTime,erpv);
    sunmoonpos(utcTime,erpv,rsun,NULL,NULL);

    /* beta and orbit angle */
    matcpy(ri,rs,6,1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    cross3(ri,ri+3,n);
    cross3(rsun,n,p);
    if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||
        !normv3(p,ep)) return 0;
    beta=PI/2.0-acos(dot(esun,en,3));
    E=acos(dot(es,ep,3));
    mu=PI/2.0+(dot(es,esun,3)<=0?-E:E);
    if (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;

    /* yaw-angle of satellite */
    if ((yaw=yaw_nominal(beta,mu))==0) return 0;

    /* satellite fixed x,y-vector */
    cross3(en,es,ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0; i<3; i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* calculate satellite antenna phase center variation ----------------------------- */
void gnss_satant_c::antmodel_s(const gnss_pcv_c *pcv,double nadir,gnss_obsd_c *data) {
    for (int i=0; i<NFREQ; i++) {
        data->antv[i]=interpvar(nadir*R2D,0,i,pcv);
    }
}
/* initialize earth rotation parameter pointor ------------------------------------ */
void gnss_satant_c::init_erp(gnss_nav_c *nav) {
    er_par=&nav->erp;
}
/* satellite antenna phase center offest for one gnss_obsd_c --------------------------- */
void gnss_satant_c::satantoff(gnss_obsd_c *data,const gnss_nav_c *nav) {
    const double *lam=nav->lam[data->sat-1];
    const gnss_pcv_c *pcv=nav->pcvs+data->sat-1;
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={ 0 };
    double gamma,C1,C2,dant1,dant2;
    int i,j=0,k=1;

    /* gps time to utc time */
    gtime_c utcTime = data->time;
    utcTime.gpst2utc();
    /* sun position in ecef */
    if (er_par) geterpv(er_par,utcTime,erpv);
    sunmoonpos(utcTime,erpv,rsun,NULL,&gmst);

    /* unit vectors of satellite fixed coordinates */
    for (i=0; i<3; i++) r[i]=-data->posvel[i];
    if (!normv3(r,ez)) return;
    for (i=0; i<3; i++) r[i]=rsun[i]-data->posvel[i];
    if (!normv3(r,es)) return;
    cross3(ez,es,r);
    if (!normv3(r,ey)) return;
    cross3(ey,ez,ex);

    //if (NFREQ>=3&&(satsys(data->sat,NULL)&(SYS_GAL|SYS_SBS))) k=2;

    if (NFREQ<2||lam[j]==0.0||lam[k]==0.0) return;
    gamma=SQR(lam[k])/SQR(lam[j]);
    C1=gamma/(gamma-1.0);
    C2=-1.0 /(gamma-1.0);

    /* iono-free LC */
    for (i=0; i<3; i++) {
        /*dant1=pcv->off[0][j][0]*ex[i]+pcv->off[0][j][1]*ey[i]+pcv->off[0][j][2]*ez[i];
        dant2=pcv->off[0][k][0]*ex[i]+pcv->off[0][k][1]*ey[i]+pcv->off[0][k][2]*ez[i];
        data->posvel[i]+=C1*dant1+C2*dant2;*/
        for (int f=0; f<NFREQ; f++) {
            data->anto[f][i]=pcv->off[0][f][0]*ex[i]+pcv->off[0][f][1]*ey[i]+pcv->off[0][f][2]*ez[i];
        }
    }
}
/* satellite antenna phase center variation correction -------------------- */
void gnss_satant_c::satantvar(const vector<double> &XYZ,gnss_obsd_c *data,const gnss_nav_c *nav) {
    const gnss_pcv_c *pcv=nav->pcvs+data->sat-1;

    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;

    for (i=0; i<3; i++) {
        ru[i]=XYZ[i]-data->posvel[i];
        rz[i]=-data->posvel[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return;

    cosa=dot(eu,ez,3);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);

    antmodel_s(pcv,nadir,data);
}
/* satellite antenna phase center offset and variation correction ----------------- */
void gnss_satant_c::satantov(const gnss_prcopt_c *opt,const vector<double> &XYZ,gnss_obsd_c *data,
    gnss_nav_c *nav) {
    /* satellite PCV correction */
    satantvar(XYZ,data,nav);

    /* Note: PCO in xyz was calculated when compute satellite precise positions */
    /* satellite PCO+PCV correction */
    for (int i=0; i<NFREQ; i++) {
        /* compute satellite PCO correction and add satellite PCV correction */
        data->sant[i] = -dot(data->anto[i],data->sigvec,3) - data->antv[i];
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        lam=nav->lam[data->sat-1];
        double gamma=SQR(lam[1])/SQR(lam[0]); // f1^2/f2^2
        data->sant[NFREQ]=(gamma*data->sant[0]-data->sant[1])/(gamma-1.0);
    }
}
/* phase windup correction for one gnss_obsd_c ----------------------------------------- */
int gnss_satant_c::phase_windup(gtime_c soltime,const vector<double> &XYZ,gnss_obsd_c *data,double &phw) {
    double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],C_e_n[9];
    double dr[3],ds[3],drs[3],r[3],blh[3],cosp,ph;
    int i;

    /* satellite yaw attitude model */
    if (!sat_yaw(soltime,data->posvel,exs,eys)) return 0;

    /* unit vector satellite to receiver */
    for (i=0; i<3; i++) r[i]=XYZ[i]-data->posvel[i];
    if (!normv3(r,ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2blh(XYZ.begin(),WGS84,blh);
    blh2Cen(blh,C_e_n);
    exr[0]= C_e_n[3]; exr[1]= C_e_n[4]; exr[2]= C_e_n[5]; /* x = north */
    eyr[0]=-C_e_n[0]; eyr[1]=-C_e_n[1]; eyr[2]=-C_e_n[2]; /* y = west  */

    /* phase windup effect */
    cross3(ek,eys,eks);
    cross3(ek,eyr,ekr);
    for (i=0; i<3; i++) {
        ds[i]=exs[i]-ek[i]*dot(ek,exs,3)-eks[i];
        dr[i]=exr[i]-ek[i]*dot(ek,exr,3)+ekr[i];
    }
    cosp=dot(ds,dr,3)/norm(ds,3)/norm(dr,3);
    if (cosp<-1.0) cosp=-1.0;
    else if (cosp> 1.0) cosp= 1.0;
    ph=acos(cosp)/2.0/PI;
    cross3(ds,dr,drs);
    if (dot(ek,drs,3)<0.0) ph=-ph;

    phw=ph+floor(phw-ph+0.5); /* in cycle */
    return 1;
}

/* receiver antenna phase center offest ----------------------------------------------------------- */
/* Constuctors -------------------------------------------------------------------- */
gnss_recant_c::gnss_recant_c() {
    lam=NULL;
}
gnss_recant_c::~gnss_recant_c() {
    lam=NULL;
}
/* Implementation functions ------------------------------------------------------- */
/* receiver antenna phase center offest and variation correction ------------------ */
void gnss_recant_c::recantov(const gnss_prcopt_c *opt,int rovbas,gnss_obsd_c *data,gnss_nav_c *nav) {
    double e[3],off[3],cosel=cos(data->azel[1]);
    if (data->isys<0) return;

    e[0]=sin(data->azel[0])*cosel;
    e[1]=cos(data->azel[0])*cosel;
    e[2]=sin(data->azel[1]);

    /* antenna center correction */
    for (int i=0; i<NFREQ; i++) {
        //mean phase center and antenna position
        for (int j=0; j<3; j++) off[j] = opt->pcvr[rovbas].off[data->isys][i][j] + opt->antdel[rovbas][j];
        /* phase antenna position correction */
        // include antenna phase center variations
        data->rant[i] = dot(off,e,3) - interpvar(90.0-data->azel[1]*R2D,data->isys,i,&opt->pcvr[rovbas]);
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        lam=nav->lam[data->sat-1];
        double gamma=SQR(lam[1])/SQR(lam[0]); // f1^2/f2^2
        data->rant[NFREQ]=(gamma*data->rant[0]-data->rant[1])/(gamma-1.0);
    }
}