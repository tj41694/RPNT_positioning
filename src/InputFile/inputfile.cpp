/* Inpute Files */

#include "InputFile/inputfile.h"
#include "BaseFunction/basefunction.h"

 /* constant -------------------------------------------------------------------------------------- */
#define SQR(x)      ((x)*(x))

#define NUMSYS          6                   /* number of systems */
#define MAXRNXLEN       (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD      1024                /* max head line position */
#define MINFREQ_GLO     -7                  /* min frequency number glonass */
#define MAXFREQ_GLO     13                  /* max frequency number glonass */
#define NINCOBS         262144              /* inclimental number of obs data */
#define MAXDTIME        60					/* max difference between toc and toe (s) */

#define MIN_BDS_DTGD    1E-12               /* min difference between TGD1 and TGD2 */

const sigind_t csigind={ 0 };

static const string syscodes = "GRECJIS";	/* satellite system codes */
static const string OBSCODES = "CLDS";      /* obs type codes */
static const string frqcodes = "123456789"; /* frequency codes */

static const double ura_eph[] ={            /* ura values (ref [3] 20.3.3.3.1.1) */
    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
    3072.0,6144.0,0.0
};
static const double ura_nominal[] ={        /* ura nominal values */
    2.0,2.8,4.0,5.7,8.0,11.3,16.0,32.0,64.0,128.0,256.0,512.0,1024.0,
    2048.0,4096.0,8192.0
};

/* constant: INS input data ----------------------------------------------------------------------- */
#define N_INS_BUFF          11                  /* number of INS data input buff */
#define N_VER_NOVATEL       24                  /* number of NovAtel version */
#define NOVATEL_TEMP_SIGN   0x8000              /* NOVATEL temperature sign bit */
const string INS_NOVATEL_VERSION[N_VER_NOVATEL]={
    "HG1700-AG58",
    "HG1900",
    "HG1930",
    "HG1700-AG62",
    "HG4930-AN01",
    "SPAN CPT7",
    "IMU-CPT",
    "IMU-KVH1750",
    "IMU-FSAS",
    "LN-200",
    "ISA-100C",
    "uIMU",
    "ADIS16488",
    "IMU-IGM-A1",
    "STIM300",
    "IMU-IGM-S1",
    "G320N",
    "G320N-200Hz",
    "PwrPak7-E1",
    "PwrPak7D-E1",
    "SMART7-S",
    "G370N",
    "PwrPak7-E2",
    "PwrPak7D-E2"
};
/* INS NovAtel data scale factors ( 0: gyroscope (rad), 1: acceleration (m/s) ) */
const double INS_NOVATEL_SFACTOR[N_VER_NOVATEL][2]={
    {        pow(2.0,-33)*100, pow(2.0,-27)*FT2M*100 }, //HG1700-AG58
    {        pow(2.0,-33)*100, pow(2.0,-27)*FT2M*100 }, //HG1900
    {        pow(2.0,-33)*100, pow(2.0,-27)*FT2M*100 }, //HG1930

    {        pow(2.0,-33)*100, pow(2.0,-26)*FT2M*100 }, //HG1700-AG62

    {        pow(2.0,-33)*100,      pow(2.0,-29)*100 }, //HG4930-AN01
    {        pow(2.0,-33)*100,      pow(2.0,-29)*100 }, //SPAN CPT7

    {  0.1/(3600.0*256.0)*100,  0.05/pow(2.0,15)*100 }, //IMU-CPT
    {  0.1/(3600.0*256.0)*200,  0.05/pow(2.0,15)*200 }, //IMU-KVH1750

    { 0.1/pow(2.0,8)*AS2R*200,  0.05/pow(2.0,15)*200 }, //IMU-FSAS

    {        pow(2.0,-19)*200,      pow(2.0,-14)*200 }, //LN-200

    {         pow(1.0,-9)*200,       pow(2.0,-8)*200 }, //ISA-100C
    {         pow(1.0,-9)*200,       pow(2.0,-8)*200 }, //uIMU

    { 720/pow(2.0,31)*D2R*200,   200/pow(2.0,31)*200 }, //ADIS16488
    { 720/pow(2.0,31)*D2R*200,   200/pow(2.0,31)*200 }, //IMU-IGM-A1

    {    pow(2.0,-21)*D2R*125,      pow(2.0,-22)*125 }, //STIM300
    {    pow(2.0,-21)*D2R*125,      pow(2.0,-22)*125 }, //IMU-IGM-S1

    {         0.008/65536*D2R,   0.2/65536*M_GRAVITY }, //G320N
    {         0.008/65536*D2R,   0.2/65536*M_GRAVITY }, //G320N-200Hz

    {         0.008/65536*D2R,   0.2/65536*M_GRAVITY }, //PwrPak7-E1
    {         0.008/65536*D2R,   0.2/65536*M_GRAVITY }, //PwrPak7D-E1
    {         0.008/65536*D2R,   0.2/65536*M_GRAVITY }, //SMART7-S

    {     0.0151515/65536*D2R,   0.4/65536*M_GRAVITY }, //G370N
    {     0.0151515/65536*D2R,   0.4/65536*M_GRAVITY }, //PwrPak7-E2
    {     0.0151515/65536*D2R,   0.4/65536*M_GRAVITY }, //PwrPak7D-E2
};
/* INS NovAtel temperature factors ( temperature = tempFactor[0]*LSB + tempFactor[1] ) */
const double INS_NOVATEL_TFACTOR[N_VER_NOVATEL][2]={
    {           0,          0 }, //HG1700-AG58
    {           0,          0 }, //HG1900
    {           0,          0 }, //HG1930

    {           0,          0 }, //HG1700-AG62

    {           0,          0 }, //HG4930-AN01
    {           0,          0 }, //SPAN CPT7

    {           0,          0 }, //IMU-CPT
    {           1,          0 }, //IMU-KVH1750

    {           0,          0 }, //IMU-FSAS

    {           0,          0 }, //LN-200

    { pow(2.0,-8),          0 }, //ISA-100C
    { pow(2.0,-8),          0 }, //uIMU

    {     0.00565,         25 }, //ADIS16488
    {     0.00565,         25 }, //IMU-IGM-A1

    { pow(2.0,-8),          0 }, //STIM300
    { pow(2.0,-8),          0 }, //IMU-IGM-S1

    {  -0.0037918, 34.9876012 }, //G320N
    {  -0.0037918, 34.9876012 }, //G320N-200Hz

    {  -0.0037918, 34.9876012 }, //PwrPak7-E1
    {  -0.0037918, 34.9876012 }, //PwrPak7D-E1
    {  -0.0037918, 34.9876012 }, //SMART7-S

    {  -0.0037918, 34.9876012 }, //G370N
    {  -0.0037918, 34.9876012 }, //PwrPak7-E2
    {  -0.0037918, 34.9876012 }, //PwrPak7D-E2
};

/* constant: SONAR input data --------------------------------------------------------------------- */
#define N_SONAR_BUFF        11                  /* number of SONAR data input buff */

/* functions of parent class "in_gnss_rnx_c" ---------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
in_gnss_rnx_c::in_gnss_rnx_c() {
    opt="";
    ver=2.10;
    ver=2.10;
    ts=gtime_c();
    te=gtime_c();
    tint=0.0;
    mask=SYS_NONE;
}
in_gnss_rnx_c::in_gnss_rnx_c(string file,string option) {
    inf.open(file,ifstream::in);
    opt=option;
    ver=2.10;
    ts=gtime_c();
    te=gtime_c();
    tint=0.0;
}
in_gnss_rnx_c::~in_gnss_rnx_c() {
    inf.close();
}
/* Implementation functions ------------------------------------------------------- */
/* set system mask ---------------------------------------------------------------- */
void in_gnss_rnx_c::set_sysmask(const gnss_prcopt_c *prcopt) {
    int i,f;
    string sys;

    f=opt.find("-SYS=");
    if (f==string::npos) {
        if (prcopt) mask=prcopt->navsys; return;
        mask=SYS_ALL; return;
    }

    sys=opt.substr(f+5);
    for (i=0; i<(int)(sys.length())&&sys[i]!=' '; i++) {
        switch (sys[i]) {
            case 'G': mask|=SYS_GPS; break;
            case 'R': mask|=SYS_GLO; break;
            case 'E': mask|=SYS_GAL; break;
            case 'J': mask|=SYS_QZS; break;
            case 'C': mask|=SYS_BDS; break;
            case 'I': mask|=SYS_IRN; break;
            case 'S': mask|=SYS_SBS; break;
        }
    }
}
/* test oepn of file stream ------------------------------------------------------- */
int in_gnss_rnx_c::test_open() {
    return inf.is_open()&&!inf.eof();
}
/* close file --------------------------------------------------------------------- */
void in_gnss_rnx_c::closeF() {
    if (inf.is_open()) inf.close();
}
/* base function of read head ----------------------------------------------------- */
int in_gnss_rnx_c::readhead() {

    while (getline(inf,buff).good()) {

        if (buff.find("RINEX VERSION / TYPE")!=string::npos) {
            str2double(buff.substr(0,9),ver);//version
            type = buff.substr(20,1);	//rinex type
            switch (buff[40]) {
                case ' ':
                case 'G': sat_sys=SYS_GPS;  tsys=TSYS_GPS; break;
                case 'R': sat_sys=SYS_GLO;  tsys=TSYS_UTC; break;
                case 'E': sat_sys=SYS_GAL;  tsys=TSYS_GAL; break;
                case 'S': sat_sys=SYS_SBS;  tsys=TSYS_GPS; break;
                case 'J': sat_sys=SYS_QZS;  tsys=TSYS_QZS; break;
                case 'C': sat_sys=SYS_BDS;  tsys=TSYS_BDS; break;
                case 'I': sat_sys=SYS_IRN;  tsys=TSYS_IRN; break;
                case 'M': sat_sys=SYS_NONE; tsys=TSYS_GPS; break; /* mixed */
                default:                                   break;
            }
            break;
        }
    }

    return 1;
}
/* vritual function of read head for O,N,G,H,J,L,C -------------------------------- */
int in_gnss_rnx_c::Head() {
    return 0;
}
/* virtual function of read body for O,N,G,H,J,L,C -------------------------------- */
int in_gnss_rnx_c::Body() {
    return 0;
}

/* Class in_gnss_rnxO_c - read Observation file -----------------------------------------------------
-----------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------- */
in_gnss_rnxO_c::in_gnss_rnxO_c() {
    epochFirst=gtime_c();
    epochEnd=gtime_c();
    for (int i=0; i<MAXSAT; i++) {
        for (int j=0; j<NFREQ; j++) {
            for (int k=0; k<MAX_CODE_TYPE; k++) lastUsedCode[i][j][k]=-1;
        }
    }
}
in_gnss_rnxO_c::in_gnss_rnxO_c(string file,string option,int rcvnum) :
    in_gnss_rnx_c(file,option),rcv(rcvnum) {
    for (int i=0; i<MAXSAT; i++) {
        for (int j=0; j<NFREQ; j++) {
            for (int k=0; k<MAX_CODE_TYPE; k++) lastUsedCode[i][j][k]=-1;
        }
    }
}
in_gnss_rnxO_c::in_gnss_rnxO_c(string file,string option,gtime_c TS,gtime_c TE,int TI,int rcvnum):
    in_gnss_rnx_c(file,option),rcv(rcvnum) {
    ts=TS; te=TE; tint=TI;
    for (int i=0; i<MAXSAT; i++) {
        for (int j=0; j<NFREQ; j++) {
            for (int k=0; k<MAX_CODE_TYPE; k++) lastUsedCode[i][j][k]=-1;
        }
    }
}
in_gnss_rnxO_c::~in_gnss_rnxO_c() {
}

/* Implementation functions ------------------------------------------------------- */
/* convert rinex obs type ver.2 -> ver.3 ------------------------------------------ */
void in_gnss_rnxO_c::convcode(int sys,string str,string &tp) {
    if (!str.compare("P1")) { /* ver.2.11 GPS L1PY,GLO L2P */
        if (sys==SYS_GPS) tp="C1W";
        else if (sys==SYS_GLO) tp="C1P";
    }
    else if (!str.compare("P2")) { /* ver.2.11 GPS L2PY,GLO L2P */
        if (sys==SYS_GPS) tp="C2W";
        else if (sys==SYS_GLO) tp="C2P";
    }
    else if (!str.compare("C1")) { /* ver.2.11 GPS L1C,GLO L1C/A */
        if (ver>=2.12); /* reject C1 for 2.12 */
        else if (sys==SYS_GPS) tp="C1C";
        else if (sys==SYS_GLO) tp="C1C";
        else if (sys==SYS_GAL) tp="C1X"; /* ver.2.12 */
        else if (sys==SYS_QZS) tp="C1C";
        else if (sys==SYS_SBS) tp="C1C";
        else if (sys==SYS_BDS) tp="C2X"; /* ver.2.12 B1 */
    }
    else if (!str.compare("C2")) {
        if (sys==SYS_GPS) {
            if (ver>=2.12) tp="C2W"; /* L2P(Y) */
            else           tp="C2X"; /* L2C */
        }
        else if (sys==SYS_GLO) tp="C2C";
        else if (sys==SYS_QZS) tp="C2X";
        else if (sys==SYS_BDS) tp="C2X"; /* ver.2.12 B1 */
    }
    else if (ver>=2.12&&str[1]=='A') { /* ver.2.12 L1C/A */
        if (sys==SYS_GPS) tp=str.substr(0,1)+"1C";
        else if (sys==SYS_GLO) tp=str.substr(0,1)+"1C";
        else if (sys==SYS_QZS) tp=str.substr(0,1)+"1C";
        else if (sys==SYS_SBS) tp=str.substr(0,1)+"1C";
    }
    else if (ver>=2.12&&str[1]=='B') { /* ver.2.12 GPS L1C */
        if (sys==SYS_GPS) tp=str.substr(0,1)+"1X";
        else if (sys==SYS_QZS) tp=str.substr(0,1)+"1X";
    }
    else if (ver>=2.12&&str[1]=='C') { /* ver.2.12 GPS L2C */
        if (sys==SYS_GPS) tp=str.substr(0,1)+"2X";
        else if (sys==SYS_QZS) tp=str.substr(0,1)+"2X";
    }
    else if (ver>=2.12&&str[1]=='D') { /* ver.2.12 GLO L2C/A */
        if (sys==SYS_GLO) tp=str.substr(0,1)+"2C";
    }
    else if (ver>=2.12&&str[1]=='1') { /* ver.2.12 GPS L1PY,GLO L1P */
        if (sys==SYS_GPS) tp=str.substr(0,1)+"1W";
        else if (sys==SYS_GLO) tp=str.substr(0,1)+"1P";
        else if (sys==SYS_GAL) tp=str.substr(0,1)+"1X"; /* tentative */
        else if (sys==SYS_BDS) tp=str.substr(0,1)+"2X"; /* extension */
    }
    else if (ver<2.12&&str[1]=='1') {
        if (sys==SYS_GPS) tp=str.substr(0,1)+"1C";
        else if (sys==SYS_GLO) tp=str.substr(0,1)+"1C";
        else if (sys==SYS_GAL) tp=str.substr(0,1)+"1X"; /* tentative */
        else if (sys==SYS_QZS) tp=str.substr(0,1)+"1C";
        else if (sys==SYS_SBS) tp=str.substr(0,1)+"1C";
        else if (sys==SYS_BDS) tp=str.substr(0,1)+"2X"; /* extension */
    }
    else if (str[1]=='2') {
        if (sys==SYS_GPS) tp=str.substr(0,1)+"2W";
        else if (sys==SYS_GLO) tp=str.substr(0,1)+"2P";
        else if (sys==SYS_QZS) tp=str.substr(0,1)+"2X";
        else if (sys==SYS_BDS) tp=str.substr(0,1)+"2X"; /* ver.2.12 B1 */
    }
    else if (str[1]=='5') {
        if (sys==SYS_GPS) tp=str.substr(0,1)+"5X";
        else if (sys==SYS_GAL) tp=str.substr(0,1)+"5X";
        else if (sys==SYS_QZS) tp=str.substr(0,1)+"5X";
        else if (sys==SYS_SBS) tp=str.substr(0,1)+"5X";
    }
    else if (str[1]=='6') {
        if (sys==SYS_GAL) tp=str.substr(0,1)+"6X";
        else if (sys==SYS_QZS) tp=str.substr(0,1)+"6X";
        else if (sys==SYS_BDS) tp=str.substr(0,1)+"6X"; /* ver.2.12 B3 */
    }
    else if (str[1]=='7') {
        if (sys==SYS_GAL) tp=str.substr(0,1)+"7X";
        else if (sys==SYS_BDS) tp=str.substr(0,1)+"7X"; /* ver.2.12 B2 */
    }
    else if (str[1]=='8') {
        if (sys==SYS_GAL) tp=str.substr(0,1)+"8X";
    }
}
/* save slips --------------------------------------------------------------------- */
void in_gnss_rnxO_c::saveslips(unsigned char slips[][NFREQ],gnss_obsd_c &data)
{
    int i;
    for (i=0; i<NFREQ; i++) {
        if (data.LLI[i]&1) slips[data.sat-1][i]|=1;
        if (data.LLI[i]&2) slips[data.sat-1][i]|=2;
    }
}
/* restore slips */
void in_gnss_rnxO_c::restslips(unsigned char slips[][NFREQ],gnss_obsd_c &data)
{
    int i;
    for (i=0; i<NFREQ; i++) {
        if (slips[data.sat-1][i]&1) data.LLI[i]|=1;
        slips[data.sat-1][i]=0;
    }
}
/* set signal index --------------------------------------------------------------- */
void in_gnss_rnxO_c::set_index(int syt,string tobss[MAXOBSTYPE],sigind_t &ind,
                               const gnss_freq_c *freq)
{
    size_t p;
    string str,optstr;
    double shift;
    int i,j,k,n,isys;

    switch (syt) {
        case SYS_GPS: isys=iGPS; break;
        case SYS_GLO: isys=iGLO; break;
        case SYS_GAL: isys=iGAL; break;
        case SYS_BDS: isys=iBDS; break;
        case SYS_QZS: isys=iQZS; break;
        case SYS_IRN: isys=iIRN; break;
        case SYS_SBS: isys=iSBS; break;
        default: 	  isys=-1; return;
    }

    for (j=0; j<MAX_NF*2; j++) {
        ind.avaiFrq[j]=0;
    }

    for (i=n=0; tobss[i][0]; i++,n++) {
        ind.obsCode[i]=tobss[i];
        ind.code[i]=obs2code(tobss[i].substr(1),ind.frq+i);
        ind.type[i]=((p=OBSCODES.find(tobss[i][0]))!=string::npos)?(int)p:0;
        ind.pri[i]=getcodepri(syt,ind.code[i],opt);
        ind.usedFrq[i]=ind.BDS2Frq[i]=-1;
        ind.pos[i]=-1;

        if (ind.frq[i]<1) continue;

        /* update available frequency */
        for (j=0; j<MAX_NF; j++) {
            if ( ind.type[i]==0 && ind.frq[i]==freq->GNSS_FREQ_PRI[isys][j] ) ind.avaiFrq[j]=1;
            // BDS-2
            if ( syt==SYS_BDS && ind.type[i]==0 && ind.frq[i]==freq->GNSS_FREQ_PRI[isys][ MAX_NF + j ] )
                ind.avaiFrq[ MAX_NF + j ]=1;
        }
    }

    /* parse phase shift options */
    switch (syt) {
        case SYS_GPS: optstr="-GL%2s=%lf"; break;
        case SYS_GLO: optstr="-RL%2s=%lf"; break;
        case SYS_GAL: optstr="-EL%2s=%lf"; break;
        case SYS_BDS: optstr="-CL%2s=%lf"; break;
        case SYS_QZS: optstr="-JL%2s=%lf"; break;
        case SYS_IRN: optstr="-IL%2s=%lf"; break;
        case SYS_SBS: optstr="-SL%2s=%lf"; break;
    }
    for (p=0; p<opt.length(); p++) {
        if (opt.substr(p,3)==optstr.substr(0,3)) {
            str=opt.substr(p+3,2);
            str2double(opt.substr(p+5),shift);
            for (i=0; i<n; i++) {
                if (code2obs(ind.code[i],NULL).compare(str)!=0) continue;
                ind.shift[i]=shift;
            }
        }
    }

    ind.n=n;
}
/* decode obs epoch --------------------------------------------------------------- */
int in_gnss_rnxO_c::decode_obsepoch(gtime_c &t,int &flag,vector<int> &sats) {
    int i,j,n;
    string satid;

    if (buff.length()<32) return 0;

    if (ver<=2.99) { /* ver.2 */
        str2int(buff.substr(29,3),n);
        if (n<=0) return 0;

        /* epoch flag: 3:new site,4:header info,5:external event */
        str2int(buff.substr(28,1),flag);
        if (3<=flag&&flag<=5) return n;

        if (t.str2time(buff.substr(0,26))!=0) return 0;

        for (i=0,j=32; i<n; i++,j+=3) {
            if (j>=68) {
                if (!getline(inf,buff)) break; /* read next line */
                j=32;
            }
            if (i<MAXOBS) {
                satid=buff.substr(j,3);
                sats[i]=satid2no(satid);
            }
        }
    }
    else { /* ver.3 */
        str2int(buff.substr(32,3),n);
        if (n<=0) return 0;

        str2int(buff.substr(31,1),flag);
        if (3<flag&&flag<=5) return n;

        if (buff[0]!='>'||t.str2time(buff.substr(1,28))!=0)
            return 0;
    }
    return n;
}
/* decode obs data ---------------------------------------------------------------- */
int in_gnss_rnxO_c::decode_obsdata(gnss_obsd_c &obs) {
    sigind_t *ind;
    string satid;
    int i,j,n,m,num,stat=1,k[16],l[16];
    int *pfrq;

    /* get satellite information */
    if (ver>2.99) { /* ver.3 */
        satid=buff.substr(0,3);
        obs.sat=(unsigned int)satid2no(satid);
    }
    obs.sys=satsys(obs.sat,&obs.prn);
    obs.isys=syscd2num(obs.sys);
    if (obs.sat<1) stat=0;
    else if (!(obs.sys&mask)) stat=0;

    if (stat<1) return 0;

    /* signal index pointer */
    switch (obs.sys) {
        case SYS_GLO: ind=index+iGLO; break;
        case SYS_GAL: ind=index+iGAL; break;
        case SYS_BDS: ind=index+iBDS; break;
        case SYS_QZS: ind=index+iQZS; break;
        case SYS_IRN: ind=index+iIRN; break;
        case SYS_SBS: ind=index+iSBS; break;
        default: 	  ind=index+iGPS; break;
    }

    /* initialize val lli and p */
    for (i=0; i<ind->n; i++) {
        val[i]=0;
        lli[i]=0;
        pos[i]=-1;
        p[i]=-1;
    }

    if ( obs.sys==SYS_BDS && obs.prn<BDS3_MIN ) pfrq=ind->BDS2Frq;
    else pfrq=ind->usedFrq;

    /* read all observations */
    for (i=0,j=ver<=2.99?0:3; i<ind->n&&j+13<buff.length(); i++,j+=16) {
        str2double(buff.substr(j,14),val[i]);
        val[i] += ind->shift[i];
        if (j+14<buff.length()) {
            if (buff.substr(j+14,1).compare(" ")==0) lli[i]=0;
            else {
                str2int(buff.substr(j+14,1),num);
                lli[i]=(unsigned char)num&3;
            }
        }
        /* if the last data, read the next line */
        if ( ver<=2.99 && ( (j+17>80&&i+1<ind->n) || 
            (j+16>=buff.length()&&i+1+(int)((80-j-16)/16)<ind->n) ) ) { /* ver.2 */
            if (!getline(inf,buff)) break; /* read next line */
            i+=(int)((80-j-16)/16);
            j=-16;
        }
    }

    /* get available observation position with top priority */
    for (i=0; i<NFREQ; i++) { // loop of frequency
        for (int tp=0; tp<OBSCODES.size(); tp++) { // loop of observation types
            for (j=0,m=-1; j<ind->n; j++) { // loop of observation data
                if ( val[j]!=0 && pfrq[j]==i+1 && 
                    ind->type[j]==tp && ind->pri[j] && 
                    ( m<0 || lastUsedCode[obs.sat-1][i][tp]==j || ind->pri[j]>ind->pri[m] ) ) {
                    m=j;
                    if (lastUsedCode[obs.sat-1][i][tp]==j) break;
                }
            }
            if (m<0) continue;

            for (j=0; j<ind->n; j++) {
                if (ind->type[j]==tp&&ind->code[j]==ind->code[m]) {
                    pos[j]=i;
                    if ( tp==1 && lastUsedCode[obs.sat-1][i][tp]!=j && lastUsedCode[obs.sat-1][i][tp]>-1 ) {
                        lli[j]|=3;
                        obs.warningMsg="Carrier Phase changed from '" + 
                            ind->obsCode[ lastUsedCode[obs.sat-1][i][tp] ] + "' to '" + ind->obsCode[j];
                    }
                    lastUsedCode[obs.sat-1][i][tp]=j;
                }
            }
        }
    }

    for (i=n=m=0; i<ind->n; i++) {
        p[i]=ver<=2.11?pfrq[i]-1:pos[i];

        if (ind->type[i]==0&&p[i]==0) k[n++]=i; /* C1? index */
        if (ind->type[i]==0&&p[i]==1) l[m++]=i; /* C2? index */
    }
    if (ver<=2.11) {
        /* if multiple codes (C1/P1,C2/P2), select higher priority */
        if (n>=2) {
            if (val[k[0]]==0.0&&val[k[1]]==0.0) {
                p[k[0]]=-1; p[k[1]]=-1;
            }
            else if (val[k[0]]!=0.0&&val[k[1]]==0.0) {
                p[k[0]]=0; p[k[1]]=-1;
            }
            else if (val[k[0]]==0.0&&val[k[1]]!=0.0) {
                p[k[0]]=-1; p[k[1]]=0;
            }
            else if (ind->pri[k[1]]>ind->pri[k[0]]) {
                p[k[1]]=0; p[k[0]]=NEXOBS<1?-1:NFREQ;
            }
            else {
                p[k[0]]=0; p[k[1]]=NEXOBS<1?-1:NFREQ;
            }
        }
        if (m>=2) {
            if (val[l[0]]==0.0&&val[l[1]]==0.0) {
                p[l[0]]=-1; p[l[1]]=-1;
            }
            else if (val[l[0]]!=0.0&&val[l[1]]==0.0) {
                p[l[0]]=1; p[l[1]]=-1;
            }
            else if (val[l[0]]==0.0&&val[l[1]]!=0.0) {
                p[l[0]]=-1; p[l[1]]=1;
            }
            else if (ind->pri[l[1]]>ind->pri[l[0]]) {
                p[l[1]]=1; p[l[0]]=NEXOBS<2?-1:NFREQ+1;
            }
            else {
                p[l[0]]=1; p[l[1]]=NEXOBS<2?-1:NFREQ+1;
            }
        }
    }
    /* save obs data */
    for (i=0; i<ind->n; i++) {
        if (p[i]<0||val[i]==0.0) continue;
        switch (ind->type[i]) {
            case 0: obs.P[p[i]]=val[i]; obs.code[p[i]]=ind->code[i]; obs.PCode[p[i]]=ind->obsCode[i]; break;
            case 1: obs.L[p[i]]=val[i]; obs.LLI[p[i]]=lli[i]; obs.LCode[p[i]]=ind->obsCode[i];        break;
            case 2: obs.D[p[i]]=(float)val[i]; obs.DCode[p[i]]=ind->obsCode[i];                       break;
            case 3: obs.SNR[p[i]]=val[i]; obs.SCode[p[i]]=ind->obsCode[i];                            break;
        }
    }

    return 1;
}

/* read "O" file to gnss_obs_c vector -------------------------------------------------- */
int in_gnss_rnxO_c::readrnxobsb(int &flag,vector<gnss_obsd_c> &data,const gnss_pstopt_c *pstopt) {
    gnss_obsd_c obs;
    gtime_c time;
    int i=0,j=0,n=0,nsat=0;
    vector<int> sats(MAXOBS,0);

    /* read record */
    while (getline(inf,buff).good()) {

        /* decode obs epoch */
        if (i==0) {
            skipflag=0;
            if ((nsat=decode_obsepoch(time,flag,sats))<=0) {
                continue;
            }
            else {
                /* utc -> gpst */
                if (tsys==TSYS_UTC) time.utc2gpst();
                time.time2str(3);

                /* screen data by time */
                if (!time.screent(ts,te,tint)) {
                    skipflag=1;
                }

                /* screen data by skipped time */
                if ( pstopt!=NULL && pstopt->nskip_period>0 ) {
                    skipflag=0;
                    for (j=0; j<pstopt->nskip_period; j++) {
                        if (time.screent(pstopt->time_skip[0][j],pstopt->time_skip[1][j],-1)) {
                            skipflag=1;
                            break;
                        }
                    }
                }
            }
        }
        else if (buff.length()<1) continue;
        else if (flag<=2||flag==6) {

            obs.reset();
            obs.time=time; obs.time.sys=satsys(sats[i-1],NULL);
            obs.sat=(unsigned char)sats[i-1];
            obs.rcv=rcv;

            /* decode obs data */
            if (decode_obsdata(obs)&&data.size()<MAXOBS) {
                data.push_back(obs); n++;
            }
        }
        if (++i>nsat) return n;

    }
    return 0;
}
/* set used frequency index according to "GNSS_USED_FREQ[MAX_GNSS][NFREQ*2]" ------ */
void in_gnss_rnxO_c::set_freq_index(gnss_nav_c *nav) {
    int isys=0,imeas=0,ifreq=0,tp=0,j=0,k=0;

    /* loop GNSS */
    for (isys=0; isys<MAX_GNSS; isys++) {
        for (imeas=0; imeas<index[isys].n; imeas++) {
            if (index[isys].frq[imeas]>0) {
                for (ifreq=0; ifreq<NFREQ; ifreq++) {
                    if (index[isys].frq[imeas]==nav->freq.GNSS_USED_FREQ[isys][ifreq])
                        index[isys].usedFrq[imeas]=ifreq+1;
                }
                // BDS-2
                if (isys==iBDS) for (ifreq=0; ifreq<NFREQ; ifreq++){
                    if (index[isys].frq[imeas]==nav->freq.GNSS_USED_FREQ[isys][ NFREQ + ifreq ])
                        index[isys].BDS2Frq[imeas]=ifreq+1;
                }
            }
        }

        /* assign index for highest priority code */
        for (ifreq=0; ifreq<NFREQ; ifreq++) { // loop of frequency
            for (int tp=0; tp<OBSCODES.size(); tp++) { // loop of observation types
                for (j=0,k=-1; j<index[isys].n; j++) { // loop of observation data
                    if (index[isys].usedFrq[j]==ifreq+1&&index[isys].type[j]==tp&&index[isys].pri[j]&&(k<0||index[isys].pri[j]>index[isys].pri[k])) {
                        k=j;
                    }
                }
                if (k<0) continue;

                for (j=0; j<index[isys].n; j++) {
                    if (index[isys].type[j]==tp&&index[isys].code[j]==index[isys].code[k]) index[isys].pos[j]=ifreq;
                }
            }
        }
    }
}
/* initialization ----------------------------------------------------------------- */
void in_gnss_rnxO_c::ini_ReadO(string file,string option,gtime_c TS,gtime_c TE,int TI,int rcvnum) {
    inf.open(file,ifstream::in);
    opt=option;
    ts=TS; te=TE; tint=TI;
    rcv=rcvnum;
}
/* read Observation file head ----------------------------------------------------- */
int in_gnss_rnxO_c::Head(gnss_nav_c *nav,gnss_sta_c *sta) {
    /* chekc inf */
    if (!inf.is_open()) return 0;

    /* default codes for unknown code */
    const string defcodes[]={
        "CWX    ",  /* GPS: L125____ */
        "CC     ",  /* GLO: L12_____ */
        "X XXXX ",  /* GAL: L1_5678_ */
        "CXXX   ",  /* QZS: L1256___ */
        "C X    ",  /* SBS: L1_5____ */
        "XX XX  ",  /* BDS: L1__67__ */
        "  A   A"   /* IRN: L__5___9 */
    };
    int i,j,k,n,nt,prn,fcn,nline=0;
    string str;

    /* base rinex  */
    if (readhead()!=1) { errorMsg="read head error!\n"; return 0; }

    while (getline(inf,buff).good()) {

        if (buff.find("MARKER NAME",60)!=string::npos && sta)
            sta->name=buff.substr(0,60);
        else if (buff.find("MARKER NUMBER",60)!=string::npos && sta)
            sta->marker=buff.substr(0,20);
        else if (buff.find("MARKER TYPE",60)!=string::npos) continue;
        else if (buff.find("OBSERVER / AGENCY",60)!=string::npos) continue;
        else if (buff.find("REC # / TYPE / VERS",60)!=string::npos && sta) {
            sta->recsno=buff.substr(0,20);
            sta->rectype=buff.substr(20,20);
            sta->recver=buff.substr(40,20);
        }
        else if (buff.find("ANT # / TYPE",60)!=string::npos && sta) {
            sta->antsno=buff.substr(0,20);
            sta->antdes=buff.substr(20,20);
        }
        else if (buff.find("APPROX POSITION XYZ",60)!=string::npos && sta)
            for (i=0; i<3; i++) {
                str2double(buff.substr(i*14,14),sta->pos[i]);
            }
        else if (buff.find("ANTENNA: DELTA H/E/N",60)!=string::npos && sta) {
            str2double(buff.substr(0,14),sta->del[2]);  /* h */
            str2double(buff.substr(14,14),sta->del[0]); /* e */
            str2double(buff.substr(28,14),sta->del[1]); /* n */
        }
        else if (buff.find("ANTENNA: DELTA X/Y/Z",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("ANTENNA: PHASECENTER",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("ANTENNA: B.SIGHT XYZ",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("ANTENNA: ZERODIR AZI",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("ANTENNA: ZERODIR XYZ",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("CENTER OF MASS: XYZ",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("SYS / # / OBS TYPES",60)!=string::npos) { /* ver.3 */
            if ((i=syscodes.find(buff[0]))==string::npos) {
                errorMsg="wrong system "+buff.substr(0,1)+"!\n"; return 0;
            }
            str2int(buff.substr(3,3),n);
            for (j=nt=0,k=7; j<n; j++,k+=4) {
                if (k>58) {
                    if (!getline(inf,buff)) break; /* read next line */
                    k=7;
                }
                if (nt<MAXOBSTYPE-1) tobs[i][nt++]=buff.substr(k,3);
            }
            tobs[i][nt]="\0";

            /* change BDS B1 code from 3.02 to 3.03 */
            if (i==iBDS&&ver>2.99&&ver<3.03)
                for (j=0; j<nt; j++) if (tobs[i][j][1]=='1') tobs[i][j].replace(1,1,"2");
            for (j=0; j<nt; j++) {
                if (tobs[i][j][2]) continue;
                if ((k=frqcodes.find(tobs[i][j][1]))==string::npos) continue;
                tobs[i][j].replace(2,1,defcodes[i].substr(k,1));
            }
        }
        else if (buff.find("WAVELENGTH FACT L1/2",60)!=string::npos) continue; /* opt ver.2 */
        else if (buff.find("# / TYPES OF OBSERV",60)!=string::npos) { /* ver.2 */
            str2int(buff.substr(0,6),n);
            for (i=nt=0,j=10; i<n; i++,j+=6) {
                if (j>58) {
                    if (!getline(inf,buff)) break;	/* read next line */
                    j=10;
                }
                if (nt>=MAXOBSTYPE-1) continue;
                if (ver<=2.99) {
                    str=buff.substr(j,2);
                    convcode(SYS_GPS,str,tobs[0][nt]);
                    convcode(SYS_GLO,str,tobs[1][nt]);
                    convcode(SYS_GAL,str,tobs[2][nt]);
                    convcode(SYS_BDS,str,tobs[3][nt]);
                    convcode(SYS_QZS,str,tobs[4][nt]);
                    convcode(SYS_SBS,str,tobs[6][nt]);
                }
                nt++;
            }
            tobs[0][nt]="\0";
        }
        else if (buff.find("SIGNAL STRENGTH UNIT",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("INTERVAL",60)!=string::npos) continue; /* opt */
        else if (buff.find("TIME OF FIRST OBS",60)!=string::npos) {
            epochFirst.str2time(buff.substr(2,41));
            if (!buff.compare(48,3,"GPS")) tsys=TSYS_GPS;
            else if (!buff.compare(48,3,"GLO")) tsys=TSYS_UTC;
            else if (!buff.compare(48,3,"GAL")) tsys=TSYS_GAL;
            else if (!buff.compare(48,3,"QZS")) tsys=TSYS_QZS; /* ver.3.02 */
            else if (!buff.compare(48,3,"BDT")) tsys=TSYS_BDS; /* ver.3.02 */
            else if (!buff.compare(48,3,"IRN")) tsys=TSYS_IRN; /* ver.3.03 */
        }
        else if (buff.find("TIME OF LAST OBS",60)!=string::npos) {
            epochEnd.str2time(buff.substr(2,41));
        }
        else if (buff.find("RCV CLOCK OFFS APPL",60)!=string::npos) continue; /* opt */
        else if (buff.find("SYS / DCBS APPLIED",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("SYS / PCVS APPLIED",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("SYS / SCALE FACTOR",60)!=string::npos) continue; /* opt ver.3 */
        else if (buff.find("SYS / PHASE SHIFTS",60)!=string::npos) continue; /* ver.3.01 */
        else if (buff.find("GLONASS SLOT / FRQ #",60)!=string::npos && nav) /* ver.3.02 */
            for (i=0; i<8; i++) {
                if (buff.compare(8*i+4,1,"R")!=0||!buff.compare(8*i+8,2,"  ")) continue;
                str2int(buff.substr(8*i+5,2),prn);
                str2int(buff.substr(8*i+8,2),fcn);
                if (1<=prn&&prn<=MAXPRNGLO) nav->glo_fcn[prn-1]=fcn+8;
            }
        else if (buff.find("GLONASS COD/PHS/BIS",60)!=string::npos && nav) /* ver.3.02 */
            for (i=0; i<4; i++) {
                if (buff.compare(13*i+1,3,"C1C"))
                    str2double(buff.substr(13*i+5,8),nav->obsbias.glo_cpbias[0]);
                else if (buff.compare(13*i+1,3,"C1P"))
                    str2double(buff.substr(13*i+5,8),nav->obsbias.glo_cpbias[1]);
                else if (buff.compare(13*i+1,3,"C2C"))
                    str2double(buff.substr(13*i+5,8),nav->obsbias.glo_cpbias[2]);
                else if (buff.compare(13*i+1,3,"C2P"))
                    str2double(buff.substr(13*i+5,8),nav->obsbias.glo_cpbias[3]);
            }
        else if (buff.find("LEAP SECONDS",60)!=string::npos && nav) {/* opt */
            str2int(buff.substr(0,6),nav->leaps);
        }
        else if (buff.find("# OF SALTELLITES",60)!=string::npos) continue;/* opt */
        else if (buff.find("PRN / # OF OBS",60)!=string::npos) continue;/* opt */
        else if (buff.find("PGM / RUN BY / DATE",60)!=string::npos) continue;
        else if (buff.find("COMMENT",60)!=string::npos) continue;
        if (buff.find("END OF HEADER",60)!=string::npos)
            break;
        if (++nline>=MAXPOSHEAD && type.compare(" ")==0) return 0; /* no rinex file */
    }

    /* set signal index */
    set_index(SYS_GPS,tobs[0],index[0],&nav->freq);
    set_index(SYS_GLO,tobs[1],index[1],&nav->freq);
    set_index(SYS_GAL,tobs[2],index[2],&nav->freq);
    set_index(SYS_BDS,tobs[3],index[3],&nav->freq);
    set_index(SYS_QZS,tobs[4],index[4],&nav->freq);
    set_index(SYS_IRN,tobs[5],index[5],&nav->freq);
    set_index(SYS_SBS,tobs[6],index[6],&nav->freq);

    return 1;
}
/* read boy of "O" file (don't use this function!!!) ------------------------------ */
int in_gnss_rnxO_c::Body(gnss_obs_c *obs,const gnss_prcopt_c *prcopt,const gnss_pstopt_c *pstopt) {
    vector<gnss_obsd_c> data;
    unsigned char slips[MAXSAT][NFREQ]={ {0} };
    int i,j,n,flag=0;

    if (!obs||!inf.is_open()) return 0;

    /* set system mask */
    set_sysmask(prcopt);

    /* read rinex obs data body */
    while ((n=readrnxobsb(flag,data,pstopt))>=0&&!inf.eof()) {
        if (te.time>0 && data[0].time.timediff(te)>=DTTOL ) return obs->data.size();

        for (i=0; i<n; i++) {
            /* save cycle-slip */
            saveslips(slips,data[i]);
        }

        if (skipflag>0) {
            data.clear();
            continue;
        }

        for (i=0; i<n; i++) {
            /* restore cycle-slip */
            restslips(slips,data[i]);

            obs->rcv=(unsigned char)rcv;

            /* save obs data */
            obs->data.push_back(data[i]);
        }
        data.clear();
    }

    obs->n=obs->data.size();
    closeF();

    return n;

}
/* read content number gnss_obs_c to the buff ------------------------------------------ */
int in_gnss_rnxO_c::Read2_Obs_Buff(vector <gnss_obs_c> &obs,const gnss_prcopt_c *prcopt,
                                   const gnss_pstopt_c *pstopt) {
    if (!inf.is_open()) return 0;

    unsigned char slips[MAXSAT][NFREQ]={ {0} };
    int i,j,n,flag=0;

    /* set system mask */
    set_sysmask(prcopt);

    obs.clear(); obs.push_back(gnss_obs_c()); obs.back().data.clear();
    /* read rinex obs data body */
    while ((n=readrnxobsb(flag,obs.back().data,pstopt))>=0&&!inf.eof()) {
        if ( te.time>0 && obs.back().data.size()>0 && obs.back().data[0].time.timediff(te)>=DTTOL ) {
            return obs.size();
        }

        for (i=0; i<n; i++) {
            /* save cycle-slip */
            saveslips(slips,obs.back().data[i]);
        }

        if (skipflag>0) {
            obs.back().data.clear();
            continue;
        }

        for (i=0; i<n; i++) {
            /* restore cycle-slip */
            restslips(slips,obs.back().data[i]);
        }
        if (n>0) {
            obs.back().n=obs.back().data.size();
            obs.back().rcv=(unsigned char)rcv;
            sortobs(obs.back());
            if (obs.size()<NOBS_BUFF) { 
                obs.push_back(gnss_obs_c()); obs.back().data.clear();
            }
            else break;
        }
    }

    return obs.size();
}
/* read one epoch body of "O" file ------------------------------------------------ */
int in_gnss_rnxO_c::Read_One_Epoch(gnss_obs_c *obs,const gnss_prcopt_c *prcopt,
                                   const gnss_pstopt_c *pstopt) {
    if (!obs||!inf.is_open()) return 0;

    obs->data.clear(); obs->n=0;
    unsigned char slips[MAXSAT][NFREQ]={ {0} };
    int i,j,n,flag=0;

    /* set system mask */
    set_sysmask(prcopt);

    /* read rinex obs data body */
    while ((n=readrnxobsb(flag,obs->data,pstopt))>=0&&!inf.eof()) {
        if ( te.time>0 && obs->data.size()>0 && obs->data[0].time.timediff(te)>=DTTOL ) {
            return n;
        }

        for (i=0; i<n; i++) {
            /* save cycle-slip */
            saveslips(slips,obs->data[i]);
        }

        if (skipflag>0) {
            obs->data.clear();
            continue;
        }

        for (i=0; i<n; i++) {
            /* restore cycle-slip */
            restslips(slips,obs->data[i]);

            obs->rcv=(unsigned char)rcv;
        }
        if (n>0) {
            obs->n=obs->data.size();
            break;
        }
    }

    sortobs(*obs);

    return n;
}

/* derived class in_gnss_rnxN_c ----------------------------------------------------------------------
* N file - Navigation ephemeris ------------------------------------------------------------------- */
/* Constructors ------------------------------------------------------------------- */
in_gnss_rnxN_c::in_gnss_rnxN_c() {
}
in_gnss_rnxN_c::in_gnss_rnxN_c(string file,string option) {
    inf.open(file,ifstream::in);
    opt=option;
    ver=2.10;
    ts=gtime_c();
    te=gtime_c();
    tint=0.0;
}
in_gnss_rnxN_c::~in_gnss_rnxN_c() {
}
/* Implementation functions ------------------------------------------------------- */
/* ura value (m) to ura index ----------------------------------------------------- */
int in_gnss_rnxN_c::uraindex(double value) {
    int i;
    for (i=0; i<15; i++) if (ura_eph[i]>=value) break;
    return i;
}
/* decode glonass ephemeris ------------------------------------------------------- */
int in_gnss_rnxN_c::decode_geph(gtime_c toc,int sat,gnss_geph_c *geph,int ndata) {
    gtime_c tof;
    double tow,tod;
    int week,dow;

    if (satsys(sat,&geph->prn)!=SYS_GLO) {
        return 0;
    }

    geph->sat=sat;

    /* toc rounded by 15 min in utc */
    tow=toc.time2gpst(&week);
    toc.gpst2time(week,floor((tow+450.0)/900.0)*900);
    dow=(int)floor(tow/86400.0);

    /* time of frame in utc */
    tod=ver<=2.99?data[2]:fmod(data[2],86400.0); /* tod (v.2), tow (v.3) in utc */
    tof.gpst2time(week,tod+dow*86400.0);
    tof.adjday(toc);

    geph->toe.copy_gtime(toc)->utc2gpst();   /* toc (gpst) */
    geph->tof.copy_gtime(toc)->utc2gpst();   /* tof (gpst) */

    /* iode = tb (7bit), tb =index of UTC+3H within current day */
    geph->iode=(int)(fmod(tow+10800.0,86400.0)/900.0+0.5);

    geph->taun=-data[0];       /* -taun */
    geph->gamn= data[1];       /* +gamman */

    geph->pos[0]=data[3]*1E3; geph->pos[1]=data[7]*1E3; geph->pos[2]=data[11]*1E3;
    geph->vel[0]=data[4]*1E3; geph->vel[1]=data[8]*1E3; geph->vel[2]=data[12]*1E3;
    geph->acc[0]=data[5]*1E3; geph->acc[1]=data[9]*1E3; geph->acc[2]=data[13]*1E3;

    geph->svh=(int)data[6];
    geph->frq=(int)data[10];
    geph->age=(int)data[14];

    /* some receiver output >128 for minus frequency number */
    if (geph->frq>128) geph->frq-=256;

    if (geph->frq<MINFREQ_GLO||MAXFREQ_GLO<geph->frq) {
        geph->svh=-1;
    }

    /* rinex 3.05 */
    if (ver>=3.05) {
        if ( ndata>15 && fabs(data[15])< 1 ) geph->dtaun=data[15];
    }

    return 1;
}
/* decode geo ephemeris ----------------------------------------------------------- */
int in_gnss_rnxN_c::decode_seph(gtime_c toc,int sat,gnss_seph_c *seph) {
    int week;

    if (satsys(sat,&seph->prn)!=SYS_SBS) {
        return 0;
    }

    seph->sat=sat;
    seph->t0 =toc;

    toc.time2gpst(&week);
    seph->tof.gpst2time(week,data[2])->adjweek(toc);

    seph->af0=data[0];
    seph->af1=data[1];

    seph->pos[0]=data[3]*1E3; seph->pos[1]=data[7]*1E3; seph->pos[2]=data[11]*1E3;
    seph->vel[0]=data[4]*1E3; seph->vel[1]=data[8]*1E3; seph->vel[2]=data[12]*1E3;
    seph->acc[0]=data[5]*1E3; seph->acc[1]=data[9]*1E3; seph->acc[2]=data[13]*1E3;

    seph->svh=(int)data[6];
    seph->sva=uraindex(data[10]);

    return 1;
}
/* decode ephemeris --------------------------------------------------------------- */
int in_gnss_rnxN_c::decode_eph(gtime_c toc,int sat,gnss_eph_c *eph) {
    int sys;

    sys=satsys(sat,&eph->prn);

    if (!(sys&(SYS_GPS|SYS_GAL|SYS_QZS|SYS_BDS|SYS_IRN))) {
        return 0;
    }

    eph->sat=sat;
    eph->toc=toc;

    eph->f0=data[0];
    eph->f1=data[1];
    eph->f2=data[2];

    eph->e=data[8]; eph->i0  =data[15]; eph->OMG0=data[13];
    eph->omg =data[17]; eph->M0 =data[6]; eph->deln=data[5]; eph->OMGd=data[18];
    eph->idot=data[19]; eph->crc=data[16]; eph->crs =data[4]; eph->cuc =data[7];
    eph->cus =data[9]; eph->cic=data[12]; eph->cis =data[14];

    if (sys==SYS_GPS||sys==SYS_QZS) {
        eph->A=SQR(data[10]);
        eph->iode=(int)data[3];       /* IODE */
        eph->iodc=(int)data[26];      /* IODC */
        eph->toes=     data[11];      /* toe (s) in gps week */
        eph->week=(int)data[21];      /* gps week */
        eph->toe.gpst2time(eph->week,data[11])->adjweek(toc);
        eph->ttr.gpst2time(eph->week,data[27])->adjweek(toc);

        eph->code=(int)data[20];      /* GPS: codes on L2 ch */
        eph->svh =(int)data[24];      /* sv health */
        eph->sva=uraindex(data[23]);  /* ura (m->index) */
        eph->flag=(int)data[22];      /* GPS: L2 P data flag */

        eph->tgd[0]=   data[25];      /* TGD */
        if (sys==SYS_GPS) {
            eph->fit=data[28];        /* fit interval (h) */
        }
        else {
            eph->fit=data[28]==0.0?1.0:2.0; /* fit interval (0:1h,1:>2h) */
        }
    }
    else if (sys==SYS_GAL) { /* GAL v.3.05 */
        eph->A=SQR(data[10]);
        eph->iode=(int)data[3];       /* IODnav */
        eph->toes=     data[11];      /* toe (s) in galileo week */
        eph->week=(int)data[21];      /* gal week = gps week */
        eph->toe.gpst2time(eph->week,data[11])->adjweek(toc);
        eph->ttr.gpst2time(eph->week,data[27])->adjweek(toc);

        eph->code=(int)data[20];      /* data sources */
                                      /* bit 0 set: I/NAV E1-B */
                                      /* bit 1 set: F/NAV E5a-I */
                                      /* bit 2 set: F/NAV E5b-I */
                                      /* bit 8 set: af0-af2 toc are for E5a.E1 */
                                      /* bit 9 set: af0-af2 toc are for E5b.E1 */
        eph->svh =(int)data[24];      /* sv health */
                                      /* bit     0: E1B DVS */
                                      /* bit   1-2: E1B HS */
                                      /* bit     3: E5a DVS */
                                      /* bit   4-5: E5a HS */
                                      /* bit     6: E5b DVS */
                                      /* bit   7-8: E5b HS */
        eph->sva =uraindex(data[23]); /* ura (m->index) */

        eph->tgd[0]=data[25];         /* BGD E5a/E1 */
        eph->tgd[1]=data[26];         /* BGD E5b/E1 */
    }
    else if (sys==SYS_BDS) { /* BeiDou v.3.05 */
        eph->toc.copy_gtime(eph->toc)->bdt2gpst();         /* bdt -> gpst */
        eph->iode=(int)data[3];                            /* AODE */
        eph->iodc=(int)data[28];                           /* AODC */
        eph->toes=     data[11];                           /* toe (s) in bdt week */
        eph->week=(int)data[21];                           /* bdt week */
        eph->code=(int)data[20];                           /* data sources */
                                                           /* bit 0 from B1C B-CNAV1 */
                                                           /* bit 1 from B2b B-CNAV1 */
                                                           /* bit 2 from B2a B-CNAV2 */
                                                           /* bit 3 from B2I broadcast */
                                                           /* bit 4 from B1I broadcast */
                                                           /* bit 5 from B3I broadcast */
                                                           /* If it is 0, source is B1I broadcast */
        eph->toe.bdt2time(eph->week,data[11])->bdt2gpst(); /* bdt -> gpst */
        eph->ttr.bdt2time(eph->week,data[27])->bdt2gpst(); /* bdt -> gpst */
        eph->toe.adjweek(toc);
        eph->ttr.adjweek(toc);

        eph->svh =(int)data[24];      /* satH1 */
        eph->sva=uraindex(data[23]);  /* ura (m->index) */

        //data message source
        if ( eph->code&BRDC_MSG_B1C || eph->code&BRDC_MSG_B2b || eph->code&BRDC_MSG_B2a ) {
            eph->Arate=data[22];
            if (binary_search(BDS_MEO,BDS_MEO+NBDS_MEO,eph->prn)) eph->A=data[10]+BRDC_AREF_MEO;
            else eph->A=data[10]+BRDC_AREF_IG;
            if (eph->code&BRDC_MSG_B1C) eph->tgd[0]=data[25]; /* TGD1 B1/B3 */
            else eph->tgd[1]=data[26]; /* TGD2 B2/B3 */
        }
        else {
            eph->A=SQR(data[10]);
            eph->tgd[0]=data[25]; /* TGD1 B1/B3 */
            if (fabs(data[25]-data[26])>MIN_BDS_DTGD) eph->tgd[1]=data[26]; /* TGD2 B2/B3 */
        }
    }
    else if (sys==SYS_IRN) { /* IRNSS v.3.03 */
        eph->A=SQR(data[10]);
        eph->iode=(int)data[3];       /* IODEC */
        eph->toes=     data[11];      /* toe (s) in irnss week */
        eph->week=(int)data[21];      /* irnss week */
        eph->toe.gpst2time(eph->week,data[11])->adjweek(toc);
        eph->ttr.gpst2time(eph->week,data[27])->adjweek(toc);
        eph->svh =(int)data[24];      /* sv health */
        eph->sva=uraindex(data[23]);  /* ura (m->index) */
        eph->tgd[0]=   data[25];      /* TGD */
    }

    if (eph->iode<0||1023<eph->iode) {
        eph->svh=-1;
    }
    if (eph->iodc<0||1023<eph->iodc) {
        eph->svh=-1;
    }
    eph->toc.time2str(3);
    eph->toe.time2str(3);

    return 1;
}
/* add data to gnss_nav_c --------------------------------------------------------- */
void in_gnss_rnxN_c::addnav(int sys,gnss_nav_c *nav) {
    if (sys==SYS_GLO) {
        nav->ngmax=2*++nav->ng;
        nav->geph.push_back(gnss_geph_c());
    }
    else if (sys==SYS_SBS) {
        nav->nsmax=2*++nav->ns;
        nav->seph.push_back(gnss_seph_c());
    }
    else if (sys==SYS_BDS||sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS) {
        nav->nmax=2*++nav->n;
        nav->eph.push_back(gnss_eph_c());
    }
}
/* delete data from gnss_nav_c ---------------------------------------------------- */
void in_gnss_rnxN_c::delnav(int sys,gnss_nav_c *nav) {
    if (sys==SYS_GLO) {
        nav->ngmax=2*--nav->ng;
        nav->geph.pop_back();
    }
    else if (sys==SYS_SBS) {
        nav->nsmax=2*--nav->ns;
        nav->seph.pop_back();
    }
    else if (sys==SYS_BDS||sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS) {
        nav->nmax=2*--nav->n;
        nav->eph.pop_back();
    }
}
/* test if the new ephemeris data is good to use ---------------------------------- */
int in_gnss_rnxN_c::test_good_eph(int sys,gnss_nav_c *nav) {
    if (last_sat==sat) {
        if (fabs(last_toc.timediff(toc))<MAXDTIME) {
            // use ephemeris with more available values
            if (last_ngood<ngood) {
                delnav(sys,nav);
                flag=1;
            }
            else flag=0;
        }
        else flag=1;
    }
    else flag=1;
    return flag;
}
/* read "O" file to gnss_obs_c vector --------------------------------------------- */
int in_gnss_rnxN_c::readrnxnavb(gnss_nav_c *nav) {
    toc=gtime_c(),last_toc=gtime_c();
    nvalue=0,ngood=0,last_ngood=0;
    int j,prn,sp=3,sys=sat_sys;
    sat=last_sat=0;
    flag=1;
    string id;

    while (getline(inf,buff).good()) {
        if (buff.compare(0,3,"   ")!=0) nvalue=ngood=0;
        /* first line */
        if (nvalue==0) {
            /* decode satellite field */
            if (ver>=3.0||sat_sys==SYS_GAL||sat_sys==SYS_QZS||sat_sys==SYS_NONE) { /* ver.3 or GAL/QZS */
                id=buff.substr(0,3);
                sat=satid2no(id);
                sp=4;
                if (ver>=3.0) sys=satsys(sat,NULL);
            }
            else {
                str2int(buff.substr(0,2),prn);

                if (sys==SYS_SBS) sat=satno(SYS_SBS,prn+100);

                else if (sys==SYS_GLO) sat=satno(SYS_GLO,prn);

                //else if (93<=prn&&prn<=97) { /* extension */
                //    sat=satno(SYS_QZS,prn+100);
                //}
                else sat=satno(SYS_GPS,prn);
            }
            /* decode toc field */
            if (toc.str2time(buff.substr(sp,19))) flag=0;
            else flag=1;
            /* decode data fields */
            for (j=0; j<3; j++) {
                str2double(buff.substr(sp+19*(j+1),19),data[nvalue++]);
                if (data[nvalue-1]!=0) ngood++;
            }
        }
        /* next line */
        else if (flag==1) {
            /* decode data fields */
            for (j=0; j<4; j++) {
                if (sp+19*(j+1)<=buff.size()) str2double(buff.substr(sp+19*j,19),data[nvalue++]);
                else data[nvalue++]=0.0;
                if (data[nvalue-1]!=0) ngood++;
            }
            /* decode ephemeris */
            if (sys==SYS_GLO&&nvalue>=15) {
                if (!(mask&sys)||!test_good_eph(sys,nav)) continue;
                addnav(sys,nav);
                decode_geph(toc,sat,&nav->geph.back(),nvalue);
                last_sat=sat; last_toc=toc;
                last_ngood=ngood;
                continue;
            }
            else if (sys==SYS_SBS&&nvalue>=15) {
                if (!(mask&sys)||!test_good_eph(sys,nav)) continue;
                addnav(sys,nav);
                decode_seph(toc,sat,&nav->seph.back());
                last_sat=sat; last_toc=toc;
                last_ngood=ngood;
                continue;
            }
            else if ( ( nvalue>=29 && ver>=3.0 ) || nvalue>=31 ) {
                if (!(mask&sys)||!test_good_eph(sys,nav)) continue;
                addnav(sys,nav);
                decode_eph(toc,sat,&nav->eph.back());
                last_sat=sat; last_toc=toc;
                last_ngood=ngood;
                continue;
            }
        }
    }
    return 1;
}
/* initialization ----------------------------------------------------------------- */
void in_gnss_rnxN_c::ini_ReadN(string file,string option) {
    inf.open(file,ifstream::in);
    opt=option;
}
/* read head of "N" file ---------------------------------------------------------- */
int in_gnss_rnxN_c::Head(gnss_nav_c *nav) {
    /* chekc inf */
    if (!inf.is_open()) return 0;

    int nline=0,i,j;
    string str;

    /* base rinex  */
    if (readhead()!=1) return 0;

    while (getline(inf,buff).good()) {
        if (buff.find("ION ALPHA",60)!=string::npos) { /* opt ver.2 */
            if (nav) {
                for (i=0,j=2; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_gps[i]);
            }
        }
        else if (buff.find("ION BETA",60)!=string::npos) { /* opt ver.2 */
            if (nav) {
                for (i=0,j=2; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_gps[i+4]);
            }
        }
        else if (buff.find("DELTA-UTC: A0,A1,T,W",60)!=string::npos) { /* opt ver.2 */
            if (nav) {
                for (i=0,j=3; i<2; i++,j+=19) str2double(buff.substr(j,12),nav->utc_gps[i]);
                for (; i<4; i++,j+=9) str2double(buff.substr(j,9),nav->utc_gps[i]);
            }
        }
        else if (buff.find("IONOSPHERIC CORR",60)!=string::npos) { /* opt ver.3 */
            if (nav) {
                if (buff.compare(0,4,"GPSA")==0) {
                    for (i=0,j=5; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_gps[i]);
                }
                else if (buff.compare(0,4,"GPSB")==0) {
                    for (i=0,j=5; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_gps[i+4]);
                }
                else if (buff.compare(0,3,"GAL")==0) {
                    for (i=0,j=5; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_gal[i]);
                }
                else if (buff.compare(0,4,"QZSA")==0) { /* v.3.02 */
                    for (i=0,j=5; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_qzs[i]);
                }
                else if (buff.compare(0,4,"QZSB")==0) { /* v.3.02 */
                    for (i=0,j=5; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_qzs[i+4]);
                }
                else if (buff.compare(0,4,"BDSA")==0) { /* v.3.02 */
                    for (i=0,j=5; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_cmp[i]);
                }
                else if (buff.compare(0,4,"BDSB")==0) { /* v.3.02 */
                    for (i=0,j=5; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_cmp[i+4]);
                }
                else if (buff.compare(0,4,"IRNA")==0) { /* v.3.03 */
                    for (i=0,j=5; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_irn[i]);
                }
                else if (buff.compare(0,4,"IRNB")==0) { /* v.3.03 */
                    for (i=0,j=5; i<4; i++,j+=12) str2double(buff.substr(j,12),nav->ion_irn[i+4]);
                }
            }
        }
        else if (buff.find("TIME SYSTEM CORR",60)!=string::npos) { /* opt ver.3 */
            if (nav) {
                if (buff.compare(0,4,"GPUT")==0) {
                    str2double(buff.substr(5,17),nav->utc_gps[0]);
                    str2double(buff.substr(22,16),nav->utc_gps[1]);
                    str2double(buff.substr(38,7),nav->utc_gps[2]);
                    str2double(buff.substr(45,5),nav->utc_gps[3]);
                }
                else if (buff.compare(0,4,"GLUT")==0) {
                    str2double(buff.substr(5,17),nav->utc_glo[0]);
                    str2double(buff.substr(22,16),nav->utc_glo[1]);
                    str2double(buff.substr(38,7),nav->utc_glo[2]);
                    str2double(buff.substr(45,5),nav->utc_glo[3]);
                }
                else if (buff.compare(0,4,"GAUT")==0) { /* v.3.02 */
                    str2double(buff.substr(5,17),nav->utc_gal[0]);
                    str2double(buff.substr(22,16),nav->utc_gal[1]);
                    str2double(buff.substr(38,7),nav->utc_gal[2]);
                    str2double(buff.substr(45,5),nav->utc_gal[3]);
                }
                else if (buff.compare(0,4,"QZUT")==0) { /* v.3.02 */
                    str2double(buff.substr(5,17),nav->utc_qzs[0]);
                    str2double(buff.substr(22,16),nav->utc_qzs[1]);
                    str2double(buff.substr(38,7),nav->utc_qzs[2]);
                    str2double(buff.substr(45,5),nav->utc_qzs[3]);
                }
                else if (buff.compare(0,4,"BDUT")==0) { /* v.3.02 */
                    str2double(buff.substr(5,17),nav->utc_cmp[0]);
                    str2double(buff.substr(22,16),nav->utc_cmp[1]);
                    str2double(buff.substr(38,7),nav->utc_cmp[2]);
                    str2double(buff.substr(45,5),nav->utc_cmp[3]);
                }
                else if (buff.compare(0,4,"SBUT")==0) { /* v.3.02 */
                    str2double(buff.substr(5,17),nav->utc_cmp[0]);
                    str2double(buff.substr(22,16),nav->utc_cmp[1]);
                    str2double(buff.substr(38,7),nav->utc_cmp[2]);
                    str2double(buff.substr(45,5),nav->utc_cmp[3]);
                }
                else if (buff.compare(0,4,"IRUT")==0) { /* v.3.03 */
                    str2double(buff.substr(5,17),nav->utc_irn[0]);
                    str2double(buff.substr(22,16),nav->utc_irn[1]);
                    str2double(buff.substr(38,7),nav->utc_irn[2]);
                    str2double(buff.substr(45,5),nav->utc_irn[3]);
                }
            }
        }
        else if (buff.find("LEAP SECONDS",60)!=string::npos) { /* opt */
            if (nav) str2int(buff.substr(0,6),nav->leaps);
        }
        if (buff.find("END OF HEADER",60)!=string::npos) return 1;
        if (++nline>=MAXPOSHEAD && type.compare(" ")==0) break; /* no rinex file */
    }
    return 0;
}
/* read boy of "N" file ----------------------------------------------------------- */
int in_gnss_rnxN_c::Body(gnss_nav_c *nav,const gnss_prcopt_c *prcopt) {
    if (!nav) return 0;

    set_sysmask(prcopt);

    readrnxnavb(nav);

    closeF();
    return nav->n>0||nav->ng>0||nav->ns>0;
}

/* derived class in_gnss_rnxC_c ----------------------------------------------------------------------
* C file - Precise Clock file --------------------------------------------------------------------- */
in_gnss_rnxC_c::in_gnss_rnxC_c() {
    ver=3.02;
    nameLen=4;
    iTime=4+nameLen;
    iClk=36+nameLen;
    iVar=56+nameLen;
    lClk=55+nameLen;
    lVar=75+nameLen;
}
in_gnss_rnxC_c::in_gnss_rnxC_c(string file,string option) {
    inf.open(file,ifstream::in);
    opt=option;
    ver=3.02;
    nameLen=4;
    iTime=4+nameLen;
    iClk=36+nameLen;
    iVar=56+nameLen;
    lClk=55+nameLen;
    lVar=75+nameLen;
}
in_gnss_rnxC_c::~in_gnss_rnxC_c() {
}
/* Implementation functions ------------------------------------------------------- */
/* initialization ----------------------------------------------------------------- */
void in_gnss_rnxC_c::ini_ReadC(string file,string option) {
    inf.open(file,ifstream::in);
    opt=option;
}
/* read head of "C" file ---------------------------------------------------------- */
int in_gnss_rnxC_c::Head(gnss_nav_c *nav) {
    /* chekc inf */
    if (!inf.is_open()) return 0;

    int nline=0;
    string str;

    /* base rinex  */
    if (readhead()!=1) return 0;

    /* length of name */
    if (ver>3.03) nameLen=9;
    else nameLen=4;
    iTime=4+nameLen;
    iClk=36+nameLen;
    iVar=56+nameLen;
    lClk=55+nameLen;
    lVar=75+nameLen;

    /* time system */
    while (getline(inf,buff).good()) {
        if (buff.find("TIME SYSTEM ID")!=string::npos) {
            if (!buff.compare(3,3,"GPS")) tsys=TSYS_GPS;
            else if (!buff.compare(3,3,"GLO")) tsys=TSYS_UTC;
            else if (!buff.compare(3,3,"GAL")) tsys=TSYS_GAL;
            else if (!buff.compare(3,3,"QZS")) tsys=TSYS_QZS; /* ver.3.02 */
            else if (!buff.compare(3,3,"BDT")) tsys=TSYS_BDS; /* ver.3.02 */
            else if (!buff.compare(3,3,"IRN")) tsys=TSYS_IRN; /* ver.3.03 */
            continue;
        }
        else if (buff.find("END OF HEADER")!=string::npos) return 1;
        if (++nline>=MAXPOSHEAD && type.compare(" ")==0) break; /* no rinex file */
    }
    return 0;
}
/* read boy of "C" file ----------------------------------------------------------- */
int in_gnss_rnxC_c::Body(gnss_nav_c *nav,const gnss_prcopt_c *prcopt) {
    gtime_c Ctime=gtime_c();
    gnss_pclk_c pclk;
    int sat;
    string satid="";

    if (!nav) return 0;

    /* set system mask */
    set_sysmask(prcopt);

    nav->pclk.clear(); nav->nc=nav->ncmax=0;

    while (getline(inf,buff).good()) {

        if ( buff.length()>=lClk && buff.compare(0,2,"AS")==0 ) {
            if (Ctime.str2time(buff.substr(iTime,26))!=0) {
                continue;
            }
            satid=buff.substr(3,nameLen);

            /* only read AS (satellite clock) record */
            if (!(sat=satid2no(satid))) continue;

            if (!(satsys(sat,NULL)&mask)) continue;

            double data[2]={0};
            if ( buff.length()>=lClk && !str2double(buff.substr(iClk,19),data[0]) ) continue;
            if (buff.length()>=lVar) str2double(buff.substr(iVar,19),data[1]);

            if (nav->nc<=0||fabs(Ctime.timediff(nav->pclk[nav->nc-1].time))>1E-9) {
                nav->nc++;
                nav->pclk.push_back(pclk);
                nav->pclk[nav->nc-1].time = Ctime;
            }
            nav->pclk[nav->nc-1].clk[sat-1]=data[0];
            nav->pclk[nav->nc-1].std[sat-1]=data[1];
        }
        else continue;

    }
    nav->ncmax=nav->nc*2;
    return nav->nc>0;
}
/* read precise clocks file ------------------------------------------------------- */
int in_gnss_rnxC_c::readclk(gnss_nav_c *nav,const gnss_prcopt_c *opt) {
    Head(nav);
    Body(nav,opt);

    closeF();

    return 1;
}

/* read earth rotation parameters file ------------------------------------------------------------ */
/* Constructors ------------------------------------------------------------------- */
in_gnss_erp_c::in_gnss_erp_c() {
}
in_gnss_erp_c::in_gnss_erp_c(string file) {
    inf.open(file,ifstream::in);
    file_path=file;
}
in_gnss_erp_c::~in_gnss_erp_c() {
    closeF();
}
/* Implementation functions ------------------------------------------------------- */
/* open file ---------------------------------------------------------------------- */
void in_gnss_erp_c::open_file() {
    inf.open(file_path);
}
/* close file --------------------------------------------------------------------- */
void in_gnss_erp_c::closeF() {
    if (inf.is_open()) inf.close();
}
/* read earth rotation parameters file -------------------------------------------- */
int in_gnss_erp_c::readerp(gnss_erp_c *erp) {
    string buff;
    double v[14]={ 0 };

    if (!inf.is_open()) return 0;

    /* initialize erp->data */
    erp->data.clear();

    while (getline(inf,buff).good()) {
        if (sscanf(buff.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10,v+11,v+12,v+13)<14)
            continue;
        erp->data.push_back(gnss_erpd_c());
        erp->data.back().mjd=v[0];
        erp->data.back().xp=v[1]*1E-6*AS2R;
        erp->data.back().yp=v[2]*1E-6*AS2R;
        erp->data.back().ut1_utc=v[3]*1E-7;
        erp->data.back().lod=v[4]*1E-7;
        erp->data.back().xpr=v[12]*1E-6*AS2R;
        erp->data.back().ypr=v[13]*1E-6*AS2R;
    }
    erp->n=erp->data.size();
    erp->nmax=erp->n+128;

    inf.close();
    return 1;
}

/* read ocean-loading tide (.BLQ) file ------------------------------------------------------------ */
/* Constructors ------------------------------------------------------------------- */
in_gnss_blq_c::in_gnss_blq_c() {
}
in_gnss_blq_c::in_gnss_blq_c(string file) {
    inf.open(file,ifstream::in);
    file_path=file;
}
in_gnss_blq_c::~in_gnss_blq_c() {
    closeF();
}
/* Implementation functions ------------------------------------------------------- */
/* open file ---------------------------------------------------------------------- */
void in_gnss_blq_c::open_file() {
    inf.open(file_path);
}
/* close file --------------------------------------------------------------------- */
void in_gnss_blq_c::closeF() {
    if (inf.is_open()) inf.close();
}
/* read earth rotation parameters file -------------------------------------------- */
int in_gnss_blq_c::readblq(const string staname,double *ocean_par) {
    string buff,NAME="    ";
    if (!inf.is_open()) inf.open(file_path,ifstream::in);

    /* 4-character station name */
    transform(staname.begin(),staname.begin()+4,NAME.begin(),::toupper);

    while (getline(inf,buff).good()) {
        if (buff.size()<2||buff.compare(0,2,"$$")==0) continue;
        /* read blq value if station name is right */
        if (buff.compare(2,4,NAME)==0) {
            double v[11]={ 0 };
            int n=0;
            while (getline(inf,buff).good()) {
                if (buff.compare(0,2,"$$")==0) continue;
                if (sscanf(buff.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10)<11) continue;
                for (int i=0; i<11; i++) ocean_par[n+i*6]=v[i];
                if (++n==6) { inf.close(); return 1; }
            }
            inf.close(); return 0;
        }
    }

    inf.close();
    return 1;
}

/* read antenna information file (.ATX) ----------------------------------------------------------- */
/* Constructors ------------------------------------------------------------------- */
in_gnss_atx_c::in_gnss_atx_c() {
}
in_gnss_atx_c::in_gnss_atx_c(string file) {
    inf.open(file,ifstream::in);
    file_path=file;
}
in_gnss_atx_c::~in_gnss_atx_c() {
    satpcv.clear();
    inf.close();
}
/* Implementation functions ------------------------------------------------------- */
/* convert frequency number to hprtk format (NFREQ) ------------------------------- */
int in_gnss_atx_c::freq2hprtk(const int prn,const string fstr,int &sys,int &ifreq,
    const gnss_nav_c *nav) {

    int i=0;

    switch (fstr[0]) {
        case 'G': sys=iGPS; break; //GPS
        case 'R': sys=iGLO; break; //GLONASS
        case 'E': sys=iGAL; break; //Galileo
        case 'C': sys=iBDS; break; //Beidou
        default: sys=-1;
    }
    /* frequency number in antenna file */
    str2int(fstr.substr(1,2),ifreq);

    if ( ifreq<1 || sys<0 ) return 0;

    /* convert ifreq to hprtk format for different systems */
    /* Beidou */
    if (sys==iBDS) {
        /* BDS-3 */
        if (prn>=BDS3_MIN) {
            for (i=0; i<NFREQ; i++) {
                if (ifreq==nav->freq.GNSS_USED_FREQ[iBDS][i]) {
                    ifreq=i+1;
                    return ifreq;
                }
            }
            /* use C07 PCO/PCV for C05 and C08 */
            if ( ifreq==7 ) {
                if ( nav->freq.GNSS_USED_FREQ[iBDS][1]==5 || nav->freq.GNSS_USED_FREQ[iBDS][1]==8 ) {
                    ifreq=2;
                    return ifreq;
                }
                if ( nav->freq.GNSS_USED_FREQ[iBDS][2]==5 || nav->freq.GNSS_USED_FREQ[iBDS][2]==8 ) {
                    ifreq=3;
                    return ifreq;
                }
            }
        }
        /* BDS-2 */
        else {
            for (i=0; i<NFREQ; i++) {
                if (ifreq==nav->freq.GNSS_USED_FREQ[iBDS][NFREQ+i]) {
                    ifreq=i+1;
                    return ifreq;
                }
            }
        }
    }
    /* other GNSS */
    else {
        for (i=0; i<NFREQ; i++) {
            if (ifreq==nav->freq.GNSS_USED_FREQ[sys][i]) {
                ifreq=i+1;
                return ifreq;
            }
        }
    }

    ifreq=0;
    return ifreq;
}
/* test satellite system ---------------------------------------------------------- */
int in_gnss_atx_c::test_sys(const int sys) {
    return sys==iGPS||sys==iGLO||sys==iGAL||sys==iBDS;
}
/* open file ---------------------------------------------------------------------- */
void in_gnss_atx_c::open_file() {
    inf.open(file_path,ifstream::in);
}
/* close file --------------------------------------------------------------------- */
void in_gnss_atx_c::closeF() {
    if (inf.is_open()) inf.close();
}
/* read satellite antenna phase center correction file ---------------------------- */
int in_gnss_atx_c::read_satatx(const gnss_nav_c *nav) {
    int stat=0;
    int sys=0,freq=0;

    if (!inf.is_open()) inf.open(file_path,ifstream::in);

    /* read atx information to pcv vector */
    while (getline(inf,buff).good()) {
        if (buff.length()<60||buff.find("COMMENT",60)!=string::npos) continue;

        /* start of one antenna information */
        if (buff.find("START OF ANTENNA",60)!=string::npos) {
            /* check if this is a satellite antenna */
            if (getline(inf,buff).good()) {
                if (buff.find("TYPE / SERIAL NO",60)!=string::npos) {
                    if (buff.compare(23,8,"        ")==0) {
                        satpcv.push_back(gnss_pcv_c());
                        satpcv.back().type=buff.substr(0,20);
                        satpcv.back().code=buff.substr(20,20);
                        /* satellite */
                        satpcv.back().sat=satid2no(satpcv.back().code);
                        satsys(satpcv.back().sat,&satpcv.back().prn);
                        stat=1;
                    }
                }
            }
        }
        /* end of one antenna information */
        if (buff.find("END OF ANTENNA",60)!=string::npos) {
            for (int f=0; f<NFREQ; f++) {
                /* use "f-1" offest if no offest for this frequency */
                if ( f>0 && norm(satpcv.back().off[0][f],3)<=0.0 )
                    matcpy(satpcv.back().off[0][f],satpcv.back().off[0][f-1],3,1);
                /* use "f-1" variantion if no variations for this frequency */
                if ( satpcv.back().var[0][f].size()<=0 &&satpcv.back().nzen>0 ) {
                    if (f==0) satpcv.back().var[0][0].assign(satpcv.back().nzen,0.0);
                    else satpcv.back().var[0][f]=satpcv.back().var[0][f-1];
                }
            }
            stat=0;
        }
        if (!stat) continue;

        else if (buff.find("ZEN1 / ZEN2 / DZEN",60)!=string::npos) {
            str2double(buff.substr(2,6),satpcv.back().zen[0]);
            str2double(buff.substr(8,6),satpcv.back().zen[1]);
            str2double(buff.substr(14,6),satpcv.back().zen[2]);
            satpcv.back().nzen =
                (int)((satpcv.back().zen[1]-satpcv.back().zen[0])/satpcv.back().zen[2])+1;
        }
        else if (buff.find("VALID FROM",60)!=string::npos) {
            if (!satpcv.back().ts.str2time(buff.substr(0,43))) continue;
        }
        else if (buff.find("VALID UNTIL",60)!=string::npos) {
            if (!satpcv.back().te.str2time(buff.substr(0,43))) continue;
        }
        else if (buff.find("START OF FREQUENCY",60)!=string::npos) {
            freq2hprtk(satpcv.back().prn,buff.substr(3,3),sys,freq,nav);
            if (freq<1||NFREQ<freq) continue;
        }
        else if (buff.find("END OF FREQUENCY",60)!=string::npos) {
            freq=0;
        }
        /* phase center offest values (X Y Z) */
        else if (buff.find("NORTH / EAST / UP",60)!=string::npos) {
            double neu[3]={ 0.0 };
            if (freq<1||NFREQ<freq) continue;
            if (sscanf(buff.c_str(),"%lf %lf %lf",neu,neu+1,neu+2)<3) continue;
            satpcv.back().off[0][freq-1][0]=1E-3*neu[0]; /* x */
            satpcv.back().off[0][freq-1][1]=1E-3*neu[1]; /* y */
            satpcv.back().off[0][freq-1][2]=1E-3*neu[2]; /* z */
        }
        /* phase center variantion values */
        else if (buff.find("NOAZI")!=string::npos) {
            if (freq<1||NFREQ<freq) continue;
            for (int i=0; i<satpcv.back().nzen; i++) {
                satpcv.back().var[0][freq-1].push_back(0.0);
                if (buff.length()>=8+i*8+8)
                    str2double(buff.substr(8+i*8,8),satpcv.back().var[0][freq-1][i]);
                satpcv.back().var[0][freq-1][i]*=1E-3;
            }
        }
    }

    closeF();
    return 1;
}
/* read receiver antenna phase center correction file ----------------------------- */
int in_gnss_atx_c::read_recatx(const gnss_prcopt_c *opt,const gnss_nav_c *nav) {
    int recbas=0,recnum=0;
    int sys=0,freq=0;

    if (!inf.is_open()) inf.open(file_path,ifstream::in);

    /* read atx information to recevier/base pcv vector */
    while (getline(inf,buff)) {
        if (buff.length()<60||buff.find("COMMENT",60)!=string::npos) continue;

        /* start of one antenna information */
        if (buff.find("START OF ANTENNA",60)!=string::npos) {
            /* check if rover or base antenna */
            if (getline(inf,buff).good()) {
                if (buff.find("TYPE / SERIAL NO",60)!=string::npos) {
                    if (opt->anttype[0].compare(0,20,buff.substr(0,20))==0)  recbas=1;
                    else if (opt->anttype[1].compare(0,20,buff.substr(0,20))==0)  recbas=2;
                    else recbas=0;
                    /* antnenna type and number */
                    if (recbas) {
                        recpcv[recbas-1].type=buff.substr(0,20);
                        recpcv[recbas-1].code=buff.substr(20,20);
                        recnum++;
                    }
                }
            }
        }
        /* end of one antenna information */
        if (buff.find("END OF ANTENNA",60)!=string::npos) {
            if (recbas>0) {
                /* for GPS */
                for (int f=0; f<NFREQ; f++) {
                    /* use "f-1" offest if no offest for this frequency */
                    if ( f>0 && norm(recpcv[recbas-1].off[0][f],3)<=0.0 )
                        matcpy(recpcv[recbas-1].off[0][f],recpcv[recbas-1].off[0][f-1],3,1);
                    /* use "f-1" variantion if no variations for this frequency */
                    if ( recpcv[recbas-1].var[0][f].size()<=0 &&recpcv[recbas-1].nzen>0 ) {
                        if (f==0) recpcv[recbas-1].var[0][0].assign(recpcv[recbas-1].nzen,0.0);
                        else recpcv[recbas-1].var[0][f]=recpcv[recbas-1].var[0][f-1];
                    }
                }
                /* for other GNSS */
                for (int j=1; j<NSYS; j++) {
                    /* use GPS antenna correction if some system correction does not exists */
                    for (int f=0; f<NFREQ; f++) {
                        /* use "f-1" offest if no offest for this frequency (only GPS) */
                        if ( norm(recpcv[recbas-1].off[j][f],3)<=0.0 )
                            matcpy(recpcv[recbas-1].off[j][f],recpcv[recbas-1].off[0][f],3,1);
                        /* use "f-1" variantion if no variations for this frequency (only GPS) */
                        if ( recpcv[recbas-1].var[j][f].size()<=0 &&recpcv[recbas-1].nzen>0 ) {
                            recpcv[recbas-1].var[j][f]=recpcv[recbas-1].var[0][f];
                        }
                    }
                }
            }
            recbas=0;
        }
        if (!recbas) continue;

        else if (buff.find("ZEN1 / ZEN2 / DZEN",60)!=string::npos) {
            str2double(buff.substr(2,6),recpcv[recbas-1].zen[0]);
            str2double(buff.substr(8,6),recpcv[recbas-1].zen[1]);
            str2double(buff.substr(14,6),recpcv[recbas-1].zen[2]);
            recpcv[recbas-1].nzen =
                (int)((recpcv[recbas-1].zen[1]-recpcv[recbas-1].zen[0])/recpcv[recbas-1].zen[2])+1;
        }
        else if (buff.find("VALID FROM",60)!=string::npos) {
            if (!recpcv[recbas-1].ts.str2time(buff.substr(0,43))) continue;
        }
        else if (buff.find("VALID UNTIL",60)!=string::npos) {
            if (!recpcv[recbas-1].te.str2time(buff.substr(0,43))) continue;
        }
        else if (buff.find("START OF FREQUENCY",60)!=string::npos) {
            freq2hprtk(0,buff.substr(3,3),sys,freq,nav);
            if (freq<1||NFREQ<freq||!test_sys(sys)) continue;
        }
        else if (buff.find("END OF FREQUENCY",60)!=string::npos) {
            freq=0;
        }
        /* phase center offest values (E N U) */
        else if (buff.find("NORTH / EAST / UP",60)!=string::npos) {
            double neu[3]={ 0.0 };
            if (freq<1||NFREQ<freq||!test_sys(sys)) continue;
            if (sscanf(buff.c_str(),"%lf %lf %lf",neu,neu+1,neu+2)<3) continue;
            recpcv[recbas-1].off[sys][freq-1][0]=1E-3*neu[1]; /* e */
            recpcv[recbas-1].off[sys][freq-1][1]=1E-3*neu[0]; /* n */
            recpcv[recbas-1].off[sys][freq-1][2]=1E-3*neu[2]; /* u */
        }
        /* phase center variantion values */
        else if (buff.find("NOAZI")!=string::npos) {
            if (freq<1||NFREQ<freq) continue;
            for (int i=0; i<recpcv[recbas-1].nzen; i++) {
                recpcv[recbas-1].var[sys][freq-1].push_back(0.0);
                if (buff.length()>=8+i*8+8)
                    str2double(buff.substr(8+i*8,8),recpcv[recbas-1].var[sys][freq-1][i]);
                recpcv[recbas-1].var[sys][freq-1][i]*=1E-3;
            }
        }
    }

    if (opt->anttype[0].compare(0,20,opt->anttype[1])==0) recpcv[1]=recpcv[0];

    closeF();
    return 1;
}

/* read precise ephemeris file -------------------------------------------------------------------- */
/* Constructors ------------------------------------------------------------------- */
in_gnss_eph_c::in_gnss_eph_c() {
    type=' ';
    ephtime=gtime_c();
    for (int i=0; i<MAXSAT; i++) satellite[i]=" ";
    bfact[0]=bfact[1]=0;
    tsys="GPS";
    ns=0;
    mask=SYS_ALL;
    opt="";
    pred_flag=1;
}
in_gnss_eph_c::in_gnss_eph_c(string file,int pred) {
    type=' ';
    ephtime=gtime_c();
    for (int i=0; i<MAXSAT; i++) satellite[i]=" ";
    bfact[0]=bfact[1]=0;
    tsys="GPS";
    ns=0;
    mask=SYS_ALL;
    opt="";
    pred_flag=pred;
    inf.open(file,ifstream::in);
    file_path=file;
}
in_gnss_eph_c::~in_gnss_eph_c() {
    closeF();
}
/* Implementation functions ------------------------------------------------------- */
void in_gnss_eph_c::set_sysmask(const gnss_prcopt_c *prcopt) {
    int i,f;
    string sys;

    f=opt.find("-SYS=");
    if (f==string::npos) {
        if (prcopt) mask=prcopt->navsys; return;
        mask=SYS_ALL; return;
    }

    sys=opt.substr(f+5);
    for (i=0; i<(int)(sys.length())&&sys[i]!=' '; i++) {
        switch (sys[i]) {
        case 'G': mask|=SYS_GPS; break;
        case 'R': mask|=SYS_GLO; break;
        case 'E': mask|=SYS_GAL; break;
        case 'J': mask|=SYS_QZS; break;
        case 'C': mask|=SYS_BDS; break;
        case 'I': mask|=SYS_IRN; break;
        case 'S': mask|=SYS_SBS; break;
        }
    }
}
/* read precise ephemeris file header --------------------------------------------- */
int in_gnss_eph_c::Head() {
    if (!inf.is_open()) return 0;
    int column=0,nsat=0;
    /* name of last block (# ,##,+ ,++,%c,%f,%i,/*) */
    string last_block;

    while (getline(inf,buff).good()) {
        /* stop first line of body */
        if (!buff.substr(0,2).compare("%i")||column>30) break;
        /* 1st line of head */
        if (column==0) {
            type=buff[2];
            ephtime.str2time(buff.substr(3,28));
        }
        /* satellite list */
        else if (!buff.substr(0,2).compare("+ ")) {
            /* first line of "+  " */
            if (!last_block.compare("##")||column==2) str2int(buff.substr(2,4),ns);
            for (int i=0; i<17&&nsat<ns; i++) {
                satellite[nsat++]=buff.substr(9+i*3,3);
            }
        }
        else if (!buff.substr(0,2).compare("%c")&&!last_block.compare("++"))
            tsys=buff.substr(9,3);
        else if (!buff.substr(0,2).compare("%f")&&!last_block.compare("%c")) {
            str2double(buff.substr(3,10),bfact[0]);
            str2double(buff.substr(14,12),bfact[1]);
        }

        last_block=buff.substr(0,2);
        column++;
    }
    return 1;
}
/* read precise ephemeris file body ----------------------------------------------- */
int in_gnss_eph_c::Body(gnss_nav_c *nav) {
    if (!inf.is_open()) return 0;
    gnss_peph_c peph;
    gnss_pclk_c pclk;
    gtime_c time;
    double val,std=0.0,base;
    int sat,sys,prn,n=ns*(type=='P'?1:2),pred_o,pred_c,v;

    nav->ne=nav->nemax=0;

    while (getline(inf,buff).good()) {
        if (!buff.compare("EOF")) break; //end of file

        if (buff[0]!='*'||time.str2time(buff.substr(3,28))) continue;

        if (!tsys.compare("UTC")) time.utc2gpst(); /* utc->gpst */

        peph.time = pclk.time = time;

        for (int i=0; i<MAXSAT; i++) {
            for (int j=0; j<4; j++) {
                peph.pos[i][j]=0.0;
                peph.std[i][j]=0.0f;
                peph.vel[i][j]=0.0;
                peph.vst[i][j]=0.0f;
            }
            for (int j=0; j<3; j++) {
                peph.cov[i][j]=0.0f;
                peph.vco[i][j]=0.0f;
            }
            pclk.clk[i]=pclk.std[i]=0.0;
        }
        for (int i=pred_o=pred_c=v=0; i<n&&getline(inf,buff); i++) {

            if (buff.length()<4||(buff[0]!='P'&&buff[0]!='V')) continue;

            sys=buff[1]==' '?SYS_GPS:code2sys(buff[1]);
            str2int(buff.substr(2,2),prn);
            if (sys==SYS_SBS) prn+=100;
            else if (sys==SYS_QZS) prn+=192; /* extension to sp3-c */

            if (!(sys&mask)) continue;
            if (!(sat=satno(sys,prn))) continue;

            if (buff[0]=='P') {
                pred_c=buff.length()>=76&&buff[75]=='P';
                pred_o=buff.length()>=80&&buff[79]=='P';
            }
            for (int j=0; j<4; j++) {
                std=0.0;

                /* read option for predicted value */
                if (j< 3&&(pred_flag&1)&& pred_o) continue;
                if (j< 3&&(pred_flag&2)&&!pred_o) continue;
                if (j==3&&(pred_flag&1)&& pred_c) continue;
                if (j==3&&(pred_flag&2)&&!pred_c) continue;

                str2double(buff.substr(4+j*14,14),val);
                if (buff.length()>=80) str2double(buff.substr(61+j*3,j<3?2:3),std);

                if (buff[0]=='P') { /* position */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
                        v=1; /* valid epoch */
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
                    }
                }
                else if (v) { /* velocity */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.vel[sat-1][j]=val*(j<3?0.1:1E-10);
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.vst[sat-1][j]=(float)(pow(base,std)*(j<3?1E-7:1E-16));
                    }
                }
                /* precise clock */
                if (j==3&&v) {
                    pclk.clk[sat-1]=peph.pos[sat-1][j];
                    pclk.std[sat-1]=peph.std[sat-1][j];
                }
            }
        }
        if (v) {
            nav->peph.push_back(peph);
            nav->ne++; nav->nemax=2*nav->ne;
            nav->pclk.push_back(pclk);
            nav->nc++; nav->ncmax=2*nav->nc;
        }
    }

    return 1;
}
/* initialization --------------------------------------------------------- */
void in_gnss_eph_c::ini_readEph(string file,int pred) {
    type=' ';
    ephtime=gtime_c();
    for (int i=0; i<MAXSAT; i++) satellite[i]=" ";
    bfact[0]=bfact[1]=0;
    tsys="GPS";
    ns=0;
    mask=SYS_ALL;
    opt="";
    pred_flag=pred;
    inf.open(file,ifstream::in);
    file_path=file;
}
/* open file ---------------------------------------------------------------------- */
void in_gnss_eph_c::open_file() {
    inf.open(file_path,ifstream::in);
}
/* close file --------------------------------------------------------------------- */
void in_gnss_eph_c::closeF() {
    if (inf.is_open()) inf.close();
}
/* read precise ephemeris file ---------------------------------------------------- */
int in_gnss_eph_c::readsp3(gnss_nav_c *nav,const gnss_prcopt_c *opt) {
    set_sysmask(opt);

    Head();
    Body(nav);

    closeF();

    return 1;
}
/* read ionex tec grid file ----------------------------------------------------------------------- */
/* Constructors ------------------------------------------------------------------- */
in_gnss_ionex_c::in_gnss_ionex_c() {
    hgts[0]=hgts[1]=450; hgts[2]=0.0;
    lats[0]=87.5; lats[1]=-87.5; lats[2]=-2.5;
    lons[0]=-180.0; lons[1]=180.0; lons[2]=5.0;
    REarth=6371.0;
    tec_factor=-1.0;
}
in_gnss_ionex_c::in_gnss_ionex_c(string file) {
    hgts[0]=hgts[1]=450; hgts[2]=0.0;
    lats[0]=87.5; lats[1]=-87.5; lats[2]=-2.5;
    lons[0]=-180.0; lons[1]=180.0; lons[2]=5.0;
    REarth=6371.0;
    tec_factor=-1.0;
    /* file and stream */
    file_path=file;
    inf.open(file,ifstream::in);
}
in_gnss_ionex_c::~in_gnss_ionex_c() {
    closeF();
}
/* Implementation functions ------------------------------------------------------- */
/* data index (i:lat,j:lon,k:hgt) ------------------------------------------------- */
int in_gnss_ionex_c::dataindex(int i,int j,int k,const int *ndata) {
    if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
    return i+ndata[0]*(j+ndata[1]*k);
}
/* read P1P2 DCB ------------------------------------------------------------------ */
void in_gnss_ionex_c::P1P2DCB(gnss_nav_c *nav) {
    int sat;

    for (int i=0; i<MAXSAT; i++) nav->obsbias.P1P2[i]=0.0;

    while (getline(inf,buff).good()) {
        if (buff.length()<60) continue;

        if (buff.find("PRN / BIAS / RMS")!=string::npos) {

            if (!(sat=satid2no(buff.substr(3,3)))) {
                continue;
            }
            str2double(buff.substr(6,10),nav->obsbias.P1P2[sat-1]);
        }
        else if (buff.find("END OF AUX DATA")!=string::npos) break;
    }
}
/* read head of ionex tec grid file ----------------------------------------------- */
int in_gnss_ionex_c::Head(gnss_nav_c *nav) {
    if (!inf.is_open()) return 0;
    while (getline(inf,buff).good()) {
        if (buff.length()<60) continue;

        if (buff.find("IONEX VERSION / TYPE")!=string::npos) {
            if (buff[20]=='I') str2double(buff.substr(0,8),version);
        }
        else if (buff.find("BASE RADIUS")!=string::npos) {
            str2double(buff.substr(0,8),REarth);
        }
        else if (buff.find("HGT1 / HGT2 / DHGT")!=string::npos) {
            str2double(buff.substr(2,6),hgts[0]);
            str2double(buff.substr(8,6),hgts[1]);
            str2double(buff.substr(14,6),hgts[2]);
        }
        else if (buff.find("LAT1 / LAT2 / DLAT")!=string::npos) {
            str2double(buff.substr(2,6),lats[0]);
            str2double(buff.substr(8,6),lats[1]);
            str2double(buff.substr(14,6),lats[2]);
        }
        else if (buff.find("LON1 / LON2 / DLON")!=string::npos) {
            str2double(buff.substr(2,6),lons[0]);
            str2double(buff.substr(8,6),lons[1]);
            str2double(buff.substr(14,6),lons[2]);
        }
        else if (buff.find("EXPONENT")!=string::npos) {
            str2double(buff.substr(0,6),tec_factor);
        }
        else if (buff.find("START OF AUX DATA")!=string::npos&&
            buff.find("DIFFERENTIAL CODE BIASES")!=string::npos) {
            P1P2DCB(nav);
        }
        else if (buff.find("END OF HEADER")!=string::npos) {
            return version;
        }
    }
    return 0;
}
/* add one epoch tec map data to gnss_nav_c ------------------------------------ */
gnss_tec_c*  in_gnss_ionex_c::addtec2nav(gnss_nav_c *nav) {
    int ndata[3];

    ndata[0]=nitem(lats);
    ndata[1]=nitem(lons);
    ndata[2]=nitem(hgts);
    if (ndata[0]<=1||ndata[1]<=1||ndata[2]<=0) return NULL;

    nav->tec.push_back(gnss_tec_c());
    if (nav->nt>=nav->ntmax) {
        nav->ntmax+=256;
    }
    nav->tec.back().rb=REarth;
    for (int i=0; i<3; i++) {
        nav->tec.back().ndata[i]=ndata[i];
        nav->tec.back().lats[i]=lats[i];
        nav->tec.back().lons[i]=lons[i];
        nav->tec.back().hgts[i]=hgts[i];
    }
    int n=ndata[0]*ndata[1]*ndata[2];

    nav->tec.back().data.assign(n,0.0);
    nav->tec.back().rms.assign(n,0.0);

    for (int i=0; i<n; i++) {
        nav->tec.back().data[i]=0.0;
        nav->tec.back().rms[i]=0.0f;
    }
    nav->nt++;
    return &nav->tec.back();
}
/* read body of ionex tec grid file ----------------------------------------------- */
int in_gnss_ionex_c::Body(gnss_nav_c *nav) {
    if (!inf.is_open()) return 0;

    gnss_tec_c *p=NULL;
    int type=0;

    while (getline(inf,buff).good()) {
        if (buff.length()<60) continue;

        if (buff.find("START OF TEC MAP")!=string::npos) {
            if (p=(addtec2nav(nav))) type=1;
        }
        else if (buff.find("END OF TEC MAP")!=string::npos) {
            type=0;
            p=NULL;
        }
        else if (buff.find("START OF RMS MAP")!=string::npos) {
            type=2;
            p=NULL;
        }
        else if (buff.find("END OF RMS MAP")!=string::npos) {
            type=0;
            p=NULL;
        }
        else if (buff.find("EPOCH OF CURRENT MAP")!=string::npos) {
            if (iontime.str2time(buff.substr(0,36))) continue;

            if (type==2) {
                for (int i=nav->nt-1; i>=0; i--) {
                    if (fabs(iontime.timediff(nav->tec[i].time))>=1.0) continue;
                    p=&nav->tec[i];
                    break;
                }
            }
            else if (p) p->time=iontime;
        }
        else if (buff.find("LAT/LON1/LON2/DLON/H")!=string::npos&&p) {
            double lat,lon[3],hgt,x;
            str2double(buff.substr(2,6),lat);
            str2double(buff.substr(8,6),lon[0]);
            str2double(buff.substr(14,6),lon[1]);
            str2double(buff.substr(20,6),lon[2]);
            str2double(buff.substr(26,6),hgt);

            int i=getindex(lat,p->lats);
            int k=getindex(hgt,p->hgts);
            int n=nitem(lon);

            for (int m=0; m<n; m++) {
                int index;

                if (m%16==0&&!getline(inf,buff)) break;

                int j=getindex(lon[0]+lon[2]*m,p->lons);
                if ((index=dataindex(i,j,k,p->ndata))<0) continue;

                str2double(buff.substr(m%16*5,5),x);
                if (fabs(x-9999.0)<=1E-6) continue;

                if (type==1) p->data[index]=x*pow(10.0,tec_factor);
                else p->rms[index]=(float)(x*pow(10.0,tec_factor));
            }
        }
    }
    return 1;
}
/* initialization ----------------------------------------------------------------- */
void in_gnss_ionex_c::ini_readIonex(string file) {
    hgts[0]=hgts[1]=450; hgts[2]=0.0;
    lats[0]=87.5; lats[1]=-87.5; lats[2]=-2.5;
    lons[0]=-180.0; lons[1]=180.0; lons[2]=5.0;
    REarth=6371.0;
    tec_factor=-1.0;
    /* file and stream */
    file_path=file;
    inf.open(file,ifstream::in);
}
/* open file ---------------------------------------------------------------------- */
void in_gnss_ionex_c::open_file() {
    inf.open(file_path,ifstream::in);
}
/* close file --------------------------------------------------------------------- */
void in_gnss_ionex_c::closeF() {
    if (inf.is_open()) inf.close();
}
/* read ionex tec grid file ------------------------------------------------------- */
int in_gnss_ionex_c::readIonex(gnss_nav_c *nav) {
    if (Head(nav)<1) { closeF(); return 0; }
    if (Body(nav)<1) { closeF(); return 0; }

    closeF();
    return 1;
}

/* read GNSS code bias file ----------------------------------------------------------------------- */
/* Constructors ------------------------------------------------------------------- */
in_gnss_bias_c::in_gnss_bias_c() {
    tsys=TSYS_GPS;
    type=mode=0;
    tsYear=tsDoy=tsSssss=0;
    teYear=teDoy=teSssss=0;
}
in_gnss_bias_c::in_gnss_bias_c(string file) {
    tsys=TSYS_GPS;
    type=mode=0;
    tsYear=tsDoy=tsSssss=0;
    teYear=teDoy=teSssss=0;
    /* file and stream */
    file_path=file;
    inf.open(file,ifstream::in);
}
in_gnss_bias_c::~in_gnss_bias_c() {
    closeF();
}
/* Implementation functions ------------------------------------------------------- */
/* read head of GNSS code bias file ----------------------------------------------- */
int in_gnss_bias_c::Head(gnss_nav_c *nav) {
    /* chekc inf */
    if (!inf.is_open()) return 0;

    mode=0;
    while (getline(inf,buff).good()) {
        if (buff.length()<80) continue;
        /* available time */
        if (buff.compare(0,5,"%=BIA")==0) {
            agency=buff.substr(30,3);
            strY1=buff.substr(34,4); str2int(strY1,tsYear);
            strD1=buff.substr(39,3); str2int(strD1,tsDoy);
            strS1=buff.substr(43,5); str2int(strS1,tsSssss);
            strY2=buff.substr(49,4); str2int(strY2,teYear);
            strD2=buff.substr(54,3); str2int(strD2,teDoy);
            strS2=buff.substr(58,5); str2int(strS2,teSssss);
            if (tsYear<1980||teYear<1980||tsDoy<0||teDoy<0||tsDoy>366||teDoy>366||
                tsSssss<0||teSssss<0||tsSssss>86400||teSssss>86400) {
                return 0;
            }
            ts.doy2time(tsYear,tsDoy);
            ts.timeadd(tsSssss);
            te.doy2time(teYear,teDoy);
            te.timeadd(teSssss);
            nav->obsbias.cbts=ts;
            nav->obsbias.cbte=te;
        }
        else if (buff.compare(1,9,"BIAS_MODE")==0) {
            if (buff.compare(41,8,"ABSOLUTE")==0) mode=1;
            else if (buff.compare(41,8,"RELATIVE")==0) mode=0;
        }
        else if (buff.compare(1,11,"TIME_SYSTEM")==0) {
            switch (buff[41]) {
                case ' ':
                case 'M': /* mixed */
                case 'G': tsys=TSYS_GPS; break;
                case 'R': tsys=TSYS_UTC; break;
                case 'E': tsys=TSYS_GAL; break;
                case 'C': tsys=TSYS_BDS; break;
                case 'J': tsys=TSYS_QZS; break;
                case 'I': tsys=TSYS_IRN; break;
                case 'S': tsys=TSYS_GPS; break;
                default:                 break;
            }
        }
        else if (buff.compare(0,14,"*BIAS SVN_ PRN")==0) {
            return 1;
        }
    }

    /* now only absolute mode is used in the software */
    return mode;
}
/* read body of GNSS code bias file ----------------------------------------------- */
int in_gnss_bias_c::Body(gnss_nav_c *nav) {
    /* chekc inf */
    if (!inf.is_open()) return 0;

    for (int i=0; i<MAXSAT; i++) for (int j=0; j<=MAXFREQ; j++) nav->obsbias.cbias[i][j]=0;

    int sat=0,freq=0,nsat=0;
    double bias=0;
    while (getline(inf,buff).good()) {
        if (buff.length()<91) continue;
        /* now only OSB is used in the software */
        if (buff.compare(1,3,"OSB")==0) {
            sat=0; freq=0;
            bias=0;
            sat=satid2no(buff.substr(11,3));
            //Code
            if ( sat<1 || buff[25]!='C' ) continue;
            //freqeunc
            str2int(buff.substr(26,1),freq);
            if ( freq<1 || freq>9 || nav->obsbias.cbias[sat-1][freq]!=0 ) continue;
            //time
            if ( buff.substr(35,4).compare(strY1)!=0 || buff.substr(40,3).compare(strD1)!=0 || 
                 buff.substr(44,5).compare(strS1)!=0 || buff.substr(50,4).compare(strY2)!=0 ||
                 buff.substr(55,3).compare(strD2)!=0 || buff.substr(59,5).compare(strS2)!=0) 
                continue;
            str2double(buff.substr(70,21),bias);
            nav->obsbias.cbias[sat-1][freq]=bias*1E-9; //ns to s
            nsat++;
        }
    }
    return nsat;
}
/* open file ---------------------------------------------------------------------- */
void in_gnss_bias_c::open_file() {
    inf.open(file_path,ifstream::in);
}
/* close file --------------------------------------------------------------------- */
void in_gnss_bias_c::closeF() {
    if (inf.is_open()) inf.close();
}
/* read GNSS code bias file --------------------------------------------------- */
int in_gnss_bias_c::readBias(gnss_nav_c *nav) {
    if (Head(nav)<1) { closeF(); return 0; }
    if (Body(nav)<1) { closeF(); return 0; }

    closeF();
    return 1;
}

/* read GNSS SNX file ----------------------------------------------------------------------------- */
/* Constructors ------------------------------------------------------------------- */
in_gnss_snx_c::in_gnss_snx_c() {
    tsys=TSYS_GPS;
    tsYear=tsDoy=tsSssss=0;
    teYear=teDoy=teSssss=0;
    staXYZ[0]=staXYZ[1]=staXYZ[2]="";
    staName[0]=staName[1]="";
}
in_gnss_snx_c::in_gnss_snx_c(string file) {
    tsys=TSYS_GPS;
    tsYear=tsDoy=tsSssss=0;
    teYear=teDoy=teSssss=0;
    /* file and stream */
    file_path=file;
    inf.open(file,ifstream::in);
    staXYZ[0]=staXYZ[1]=staXYZ[2]="";
    staName[0]=staName[1]="";
}
in_gnss_snx_c::~in_gnss_snx_c() {
    closeF();
}
/* Implementation functions ------------------------------------------------------- */
/* read head of GNSS code bias file ----------------------------------------------- */
int in_gnss_snx_c::Head() {
    /* chekc inf */
    if (!inf.is_open()) return 0;

    if ( getline(inf,buff).good() && buff.length()>56 ) {
        /* header line */
        if (buff.compare(0,5,"%=SNX")==0) {
            agency=buff.substr(28,3);
            strY1=buff.substr(32,4); strY1="20"+strY1; str2int(strY1,tsYear);
            strD1=buff.substr(35,3); str2int(strD1,tsDoy);
            strS1=buff.substr(39,5); str2int(strS1,tsSssss);
            strY2=buff.substr(45,4); strY2="20"+strY2; str2int(strY2,teYear);
            strD2=buff.substr(48,3); str2int(strD2,teDoy);
            strS2=buff.substr(52,5); str2int(strS2,teSssss);
            if (tsYear<1980||teYear<1980||tsDoy<0||teDoy<0||tsDoy>366||teDoy>366||
                tsSssss<0||teSssss<0||tsSssss>86400||teSssss>86400) {
                return 0;
            }
            ts.doy2time(tsYear,tsDoy);
            ts.timeadd(tsSssss);
            te.doy2time(teYear,teDoy);
            te.timeadd(teSssss);

            return 1;
        }
    }

    return 0;
}
/* read body of GNSS code bias file ----------------------------------------------- */
int in_gnss_snx_c::Body(gnss_prcopt_c *opt,gnss_sta_c sta[2]) {
    /* chekc inf */
    if (!inf.is_open()) return 0;

    int solFlag=0,nxyz=0;
    staXYZ[0]=staXYZ[1]=staXYZ[2]="";
    staName[0]="0000"; staName[1]="0000";
    std::transform(sta[0].name.begin(),sta[0].name.begin()+4,staName[0].begin(),::toupper);
    std::transform(sta[1].name.begin(),sta[1].name.begin()+4,staName[1].begin(),::toupper);
    while (getline(inf,buff).good()) {
        if ( buff.length()>17 && buff.compare(0,18,"+SOLUTION/ESTIMATE")==0 ) {
            solFlag=1;
        }
        else if ( solFlag>0 && buff.length()>67 && buff.compare(14,4,staName[0],0,4)==0 ) {
            if ( buff.compare(7,4,"STAX")==0 && opt->rovpos==POSOPT_SNX ) {
                staXYZ[0]=buff.substr(47,21);
                str2double(staXYZ[0],opt->ru[0]); // f_ib_b Z
                nxyz++;
            }
            if ( buff.compare(7,4,"STAY")==0 && opt->rovpos==POSOPT_SNX ) {
                staXYZ[1]=buff.substr(47,21);
                str2double(staXYZ[1],opt->ru[1]); // f_ib_b Z
                nxyz++;
            }
            if ( buff.compare(7,4,"STAZ")==0 && opt->rovpos==POSOPT_SNX ) {
                staXYZ[2]=buff.substr(47,21);
                str2double(staXYZ[2],opt->ru[2]); // f_ib_b Z
                nxyz++;
            }
        }
        else if ( solFlag>0 && buff.length()>67 && buff.compare(14,4,staName[1],0,4)==0 ) {
            if ( buff.compare(7,4,"STAX")==0 && opt->baspos==POSOPT_SNX ) {
                staXYZ[0]=buff.substr(47,21);
                str2double(staXYZ[0],opt->rb[0]); // f_ib_b Z
                nxyz++;
            }
            if ( buff.compare(7,4,"STAY")==0 && opt->baspos==POSOPT_SNX ) {
                staXYZ[1]=buff.substr(47,21);
                str2double(staXYZ[1],opt->rb[1]); // f_ib_b Z
                nxyz++;
            }
            if ( buff.compare(7,4,"STAZ")==0 && opt->baspos==POSOPT_SNX ) {
                staXYZ[2]=buff.substr(47,21);
                str2double(staXYZ[2],opt->rb[2]); // f_ib_b Z
                nxyz++;
            }
        }
        else if ( buff.length()>17 && buff.compare(0,18,"-SOLUTION/ESTIMATE")==0 ) {
            solFlag=0;
            break;
        }
        else if ( buff.length()>6 && buff.compare(0,7,"%ENDSNX")==0 ) {
            solFlag=0;
            break;
        }
    }
    return nxyz;
}
/* open file ---------------------------------------------------------------------- */
void in_gnss_snx_c::open_file() {
    inf.open(file_path,ifstream::in);
}
/* close file --------------------------------------------------------------------- */
void in_gnss_snx_c::closeF() {
    if (inf.is_open()) inf.close();
}
/* read GNSS code bias file --------------------------------------------------- */
int in_gnss_snx_c::readSNX(gnss_prcopt_c *opt,gnss_sta_c sta[2]) {
    if (Head()<1) { closeF(); return 0; }
    if (Body(opt,sta)<1) { closeF(); return 0; }

    closeF();
    return 1;
}

/* read troposphere ZTD file (debug) -------------------------------------------------------------- */
/* Constructors ------------------------------------------------------------------- */
in_gnss_ztd_c::in_gnss_ztd_c() {
    ztd_factor=0;
}
in_gnss_ztd_c::in_gnss_ztd_c(int iztd,string file) {
    ztd_factor=iztd;
    /* file and stream */
    file_path=file;
    inf.open(file,ifstream::in);
}
in_gnss_ztd_c::~in_gnss_ztd_c() {
    closeF();
}
/* Implementation functions ------------------------------------------------------- */
/* open file ---------------------------------------------------------------------- */
void in_gnss_ztd_c::open_file() {
    inf.open(file_path,ifstream::in);
}
/* close file --------------------------------------------------------------------- */
void in_gnss_ztd_c::closeF() {
    if (inf.is_open()) inf.close();
}
/* read ionex tec grid file ------------------------------------------------------- */
int in_gnss_ztd_c::readZTD(gnss_nav_c *nav) {
    nav->ztd.clear();
    if ( !inf.is_open() || ztd_factor<1 ) return 0;

    while (getline(inf,buff).good()) {
        //get ztd time
        if ( buff.size()>29 && buff[0]=='>' && ztdtime.str2time(buff.substr(2,28))==0 ) {
            //get ztd
            if ( getline(inf,buff) && buff.size()>328 ) {
                bool ztdOk=true;
                string strValue;
                buff=buff.substr(312);
                strValue=strtok((char *)buff.c_str(),",");
                for (int i=1; i<ztd_factor; i++) {
                    strValue=strtok(NULL,",");
                }
                if (strValue.size()<1) {
                    ztdOk=false;
                    continue;
                }
                if (ztdOk==true) {
                    ztd=-1;
                    str2double(strValue,ztd);
                    if (ztd>=0 && ztd<9999) {
                        if (ztd_factor==3) ztd*=1E-3;
                        nav->ztd.push_back(gnss_ztd_c());
                        nav->ztd.back().time=ztdtime;
                        nav->ztd.back().ztd=ztd;
                    }
                }
            }
        }
    }
    return 1;
}