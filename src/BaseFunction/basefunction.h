/* Base functions for calculation -----------------------------------------------------------------
* options : -DLAPACK   use LAPACK/BLAS
*           -DMKL      use Intel MKL
*           -DWIN32    use WIN32 API
*           -DNOCALLOC no use calloc for zero matrix
*           -DIERS_MODEL use GMF instead of NMF
*           -DDLL      built for shared library
*           -DCPUTIME_IN_GPST cputime operated in gpst
*
* references :
*     [1] IS-GPS-200D, Navstar GPS Space Segment/Navigation User Interfaces,
*         7 March, 2006
*     [2] RTCA/DO-229C, Minimum operational performanc standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         RTCA inc, November 28, 2001
*     [3] M.Rothacher, R.Schmid, ANTEX: The Antenna Exchange Format Version 1.4,
*         15 September, 2010
*     [4] A.Gelb ed., Applied Optimal Estimation, The M.I.T Press, 1974
*     [5] A.E.Niell, Global mapping functions for the atmosphere delay at radio
*         wavelengths, Jounal of geophysical research, 1996
*     [6] W.Gurtner and L.Estey, RINEX The Receiver Independent Exchange Format
*         Version 3.00, November 28, 2007
*     [7] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*         May 2009
*     [8] China Satellite Navigation Office, BeiDou navigation satellite system
*         signal in space interface control document, open service signal B1I
*         (version 1.0), Dec 2012
*     [9] J.Boehm, A.Niell, P.Tregoning and H.Shuh, Global Mapping Function
*         (GMF): A new empirical mapping function base on numerical weather
*         model data, Geophysical Research Letters, 33, L07304, 2006
*     [10] GLONASS/GPS/Galileo/Compass/SBAS NV08C receiver series BINR interface
*         protocol specification ver.1.3, August, 2012
*
--------------------------------------------------------------------------------------------------- */
#ifndef BASEFUNCTION_H
#define BASEFUNCTION_H

#include "BaseFunction/timesys.h"
#include "GNSS/DataClass/data.h"

/* constants -------------------------------------------------------------------------------------- */
#define POLYCRC32   0xEDB88320u /* CRC32 polynomial */
#define POLYCRC24Q  0x1864CFBu  /* CRC24Q polynomial */
/* chi-sqr(n) (alpha=0.001) ------------------------------------------------------- */
const double ChiSqures[100] ={
    10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
    31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
    46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
    61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
    74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
    88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
    101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
    113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
    126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
    138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};
/* stream format strings ---------------------------------------------------------- */
const string FormatStrs[32] ={
    "RTCM 2",                   /*  0 */
    "RTCM 3",                   /*  1 */
    "NovAtel OEM6",             /*  2 */
    "NovAtel OEM3",             /*  3 */
    "u-blox",                   /*  4 */
    "Superstar II",             /*  5 */
    "Hemisphere",               /*  6 */
    "SkyTraq",                  /*  7 */
    "GW10",                     /*  8 */
    "Javad",                    /*  9 */
    "NVS BINR",                 /* 10 */
    "BINEX",                    /* 11 */
    "Trimble RT17",             /* 12 */
    "Septentrio",               /* 13 */
    "CMR/CMR+",                 /* 14 */
    "LEX Receiver",             /* 15 */
    "RINEX",                    /* 16 */
    "SP3",                      /* 17 */
    "RINEX CLK",                /* 18 */
    "SBAS",                     /* 19 */
    "NMEA 0183",                /* 20 */
    ""
};
/* functions -------------------------------------------------------------------------------------- */
/* decoded binary data to ASCII data -------------------------------------------------------------- */
/* extract unsigned/signed bits --------------------------------------------------- */
unsigned int getbitu(const unsigned char *ChBuff,int BitPos,int BitLen);
int getbits(const unsigned char *ChBuff,int BitPos,int BitLen);
/* set unsigned/signed bits ------------------------------------------------------- */
void setbitu(unsigned char *ChBuff,int BitPos,int BitLen,unsigned int ByteData);
/* crc-24q parity ----------------------------------------------------------------- */
unsigned int rtk_crc24q(const unsigned char *ChBuff,int BitLen);
/* crc-16 parity ------------------------------------------------------------------ */
unsigned short rtk_crc16(const unsigned char *ChBuff,int BitLen);
/* decode navigation data word ---------------------------------------------------- */
int decode_word(unsigned int NavWord,unsigned char *ChData);

/* vector and matrix functions -------------------------------------------------------------------- */
/* get element index in the vector ---------------------------------------------------
* args   : double *VecA          I   matrix to be processed
*          int    SizeA          I   size of VecA
* --------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
int vector_index(Iter1 Value,Iter2 VecA,const int SizeA) {
    for (int i=0; i<SizeA; i++) {
        if (Value==VecA[i]) return i;
    }
    return -1;
}
/* initialize identity matrix --------------------------------------------------------
* args   : double *VecA          I   matrix to be processed
*          int    dimA           I   dimension of VecA
* --------------------------------------------------------------------------------- */
template <typename Iter>
void eyemat(Iter VecA,const int dimA) {
    for (int i=0; i<dimA; i++)
        for (int j=0; j<dimA; j++)
            VecA[j+i*dimA] = i==j ? 1.0 : 0.0;
}
/* initialize identity matrix for zero matrix ----------------------------------------
* args   : double *VecA          I   matrix to be processed
*          int    dimA           I   dimension of VecA
* --------------------------------------------------------------------------------- */
template <typename Iter>
void eyemat_zero(Iter VecA,const int dimA) {
    for (int i=0; i<dimA; i++)
        VecA[i+i*dimA] = 1.0;
}
/* initialize zero matrix ------------------------------------------------------------
* args   : double *VecA          I   matrix to be processed
*          int    dimA           I   dimension of VecA
* --------------------------------------------------------------------------------- */
template <typename Iter>
void zeromat(Iter VecA,const int dimA) {
    for (int i=0; i< dimA; i++) VecA[i] = 0.0;
}
/* sum of vector and matrix ----------------------------------------------------------
* args   : double *VecA          I   vector a (n x 1)
*          int    SizeA          I   size of vector a
* return : sum of VecA
* --------------------------------------------------------------------------------- */
template <typename Iter>
double sum(const Iter VecA,const int SizeA) {
    double S=0.0;
    for (int i=0; i<SizeA; i++) {
        S+=VecA[i];
    }
    return S;
}
/* average of vector -------------------------------------------------------------- */
template <typename Iter>
double ave_vec_n0(const Iter VecA,const int SizeA) {
    double S=0.0;
    int nA=0;
    for (int i=0; i<SizeA; i++) {
        if (VecA[i]!=0.0) { 
            S+=VecA[i];
            nA++;
        }
    }
    return S/nA;
}
/* squares average of vector ------------------------------------------------------ */
template <typename Iter>
double sqr_ave_n0(const Iter VecA,const int SizeA) {
    double S=0.0;
    int nA=0;
    for (int i=0; i<SizeA; i++) {
        if (VecA[i]!=0.0) { 
            S+=VecA[i]*VecA[i];
            nA++;
        }
    }
    return S/nA;
}
/* std of vector ------------------------------------------------------------------ */
template <typename Iter>
double std_vec_n0(const Iter VecA,const int SizeA) {
    double ave=ave_vec_n0(VecA,SizeA),std=0;
    int nA=0;
    for (int i=0; i<SizeA; i++) {
        if (VecA[i]!=0.0) { 
            std+=(VecA[i]-ave)*(VecA[i]-ave);
            nA++;
        }
    }
    return nA>1? sqrt(std/nA) : 0;
}
/* rms of vector ------------------------------------------------------------------ */
template <typename Iter>
double rms_vec_n0(const Iter VecA,const int SizeA) {
    double ave=ave_vec_n0(VecA,SizeA),std=0;
    int nA=0;
    for (int i=0; i<SizeA; i++) {
        if (VecA[i]!=0.0) { 
            std+=(VecA[i]-ave)*(VecA[i]-ave);
            nA++;
        }
    }
    return nA>1? sqrt(std/(nA-1)) : 0;
}
/* update average (update NumVal) ------------------------------------------------- */
double update_ave_updateN(double &Ave,int &NumVal,const double NewVal);
/* update average with new NumVal ------------------------------------------------- */
double update_ave(double &Ave,const int NumVal,const double NewVal);
/* remove one value from average (update NumVal) ---------------------------------- */
double remove_from_ave_updateN(double &Ave,int &NumVal,const double OneVal);
/* remove one value from average with new NumVal ---------------------------------- */
double remove_from_ave(double &Ave,const int NumVal,const double OneVal);
/* update variance (update NumVal) ------------------------------------------------ */
double update_var_updateN(double &Var,int &NumVal,const double Ave,const double NewVal);
/* update variance with new NumVal ------------------------------------------------ */
double update_var(double &Var,const int NumVal,const double Ave,const double NewVal);
/* update the average for a moving window ----------------------------------------- */
double update_ave_window(double &Ave,const int NumVal,const double NewVal,const double OldVal);
/* update the variance for a moving window ---------------------------------------- */
double update_var_window(double &Var,const int NumVal,const double NewAve,const double NewVal,const double OldVal);
/* inner product ---------------------------------------------------------------------
* inner product of vectors
* args   : double *VecA,*VecB     I   vector a,b (n x 1)
*          int    SizeVec         I   size of vector a,b
* return : VecA'*VecB
* --------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
double dot(const Iter1 VecA,const Iter2 VecB,int SizeVec) {
    double dInn=0.0;

    while (--SizeVec>=0) dInn+=VecA[SizeVec]*VecB[SizeVec];
    return dInn;
}
/* euclid norm -----------------------------------------------------------------------
* euclid norm of vector
* args   : double *VecA        I   vector a (n x 1)
*          int    SizeVec      I   size of vector a
* return : | VecA |
*---------------------------------------------------------------------------------- */
template <typename Iter>
double norm(const Iter VecA,int SizeVec) {
    return sqrt(dot(VecA,VecA,SizeVec));
}
/* geometrical distance of 2 vectors -------------------------------------------------
* args   : double *VecA,*VecB  I   vector A,B (SizeAB * 1)
         : int    SizeAB       I   size of vector A,B
* return : distance between VecA and VecB
* --------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
double distance(const Iter1 VecA,const Iter2 VecB,int SizeAB) {
    vector<double> vecAB(SizeAB,0.0);
    for (int i=0; i<SizeAB; i++)
        vecAB[i] = VecA[i] - VecB[i];
    return norm(vecAB.begin(),SizeAB);
}
/* normalize 3d vector ------------------------------------------------------------ */
template <typename Iter1,typename Iter2>
int normv3(const Iter1 VecA,Iter2 VecB) {
    double r;
    if ((r=norm(VecA,3))<=0.0) return 0;
    VecB[0]=VecA[0]/r;
    VecB[1]=VecA[1]/r;
    VecB[2]=VecA[2]/r;
    return 1;
}
/* outer product of 3d vectors ---------------------------------------------------- */
template <typename Iter1,typename Iter2,typename Iter3>
void cross3(const Iter1 VecA,const Iter2 VecB,Iter3 VecC) {
    VecC[0]= VecA[1]*VecB[2] - VecA[2]*VecB[1];
    VecC[1]= VecA[2]*VecB[0] - VecA[0]*VecB[2];
    VecC[2]= VecA[0]*VecB[1] - VecA[1]*VecB[0];
}
/* trace of matrix ---------------------------------------------------------------- */
template <typename Iter>
double matrace(const Iter VecA,const int SizeA) {
    double trace=0.0;
    for (int i=0; i<SizeA; i++) {
        trace+=VecA[i+i*SizeA];
    }
    return trace;
}
/* small Euler angles to Skew-symmetric matrix --------------------------------------
* compute small Euler angles to the corresponding Skew-symmetric matrix
* args   : double *sEuler      I   small Euler angles (rad)
*          double *SSMat       O   corresponding Skew-symmetric matrix
* return : none
* note   : matirix stored by row-major order
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
void sEuler2SSM(const Iter1 sEuler,Iter2 SSMat) {
    SSMat[0]= 0;         SSMat[1]=-sEuler[2]; SSMat[2]= sEuler[1];
    SSMat[3]= sEuler[2]; SSMat[4]= 0;         SSMat[5]=-sEuler[0];
    SSMat[6]=-sEuler[1]; SSMat[7]= sEuler[0]; SSMat[8]= 0;
}
/* copy matrix -----------------------------------------------------------------------
* copy MatSrc to MatTrg
* args   : double *MatTrg       O   target matrix
*          double *MatSrc       I   source matrix
*          int     n            I   row number
*          int     m            I   column number
* return : none
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
void matcpy(Iter1 MatTrg,const Iter2 MatSrc,int n,int m) {
    for (int i=0; i<n*m; i++) MatTrg[i]= MatSrc[i];
}
/* math (vector, matrix) functions ---------------------------------------------------------------- */
/* initialize matrix with source matrix and a ratio ----------------------------------
* initialize MatTrg with MatSrc and Ratio: MatTrg = Ratio*MatSrc
* args   : double *MatTrg       O   target matrix
*          double *MatSrc       I   source matrix
*          int     n            I   row number
*          int     m            I   column number
*          double  Ratio        I   ratio
* return : none
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
void matcpy(Iter1 MatTrg,const Iter2 MatSrc,int n,int m,double Ratio) {
    for (int i=0; i<n*m; i++) MatTrg[i]= Ratio*MatSrc[i];
}
/* transpose matrix ------------------------------------------------------------------
* transpose MatTrg
* args   : double *MatTrg       O   target matrix
*          int     n            I   row number
*          int     m            I   column number
* return : none
*---------------------------------------------------------------------------------- */
template <typename Iter1>
void mattrans(Iter1 MatTrg,int n,int m) {
    vector <double> mBuff(n*m,0.0);
    for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
        mBuff[j*n+i]=MatTrg[i*m+j];
    }
    for (int i=0; i<n*m; i++) { 
        MatTrg[i]=mBuff[i]; 
    }
    mBuff.clear();
}
/* multiply matrix -------------------------------------------------------------------
*   1 NN : A(sizeA,sizeAB) * B(sizeAB,sizeB)
*   2 NT : A(sizeA,sizeAB) * B(sizeB,sizeAB)
*   3 TN : A(sizeAB,sizeA) * B(sizeAB,sizeB)
*   4 TT : A(sizeAB,sizeA) * B(sizeB,sizeAB)
* result : MatC(sizeB,sizeA)
* comment: matrix multiplication (contrary to the normal linear algebra)
* --------------------------------------------------------------------------------- */
void matmul_pnt(const string TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
    const double *MatA,const double *MatB,double CoeC,double *MatC);
void matmul_vec(const string TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
    const vector<double> &MatA,const vector<double> &MatB,double CoeC,vector<double> &MatC);
void matmul_vec(const string TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
    const vector<double> &MatA,const vector<double> &MatB,double CoeC,vector<double> &MatC,
    const int StartA,const int StartB,const int StartC);
template <typename Iter1,typename Iter2,typename Iter3>
void matmul(const string TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
    const Iter1 MatA,const Iter2 MatB,double CoeC,Iter3 MatC) {

    if (TraFlag.size()<2) return;
    double LineAB;
    int i,j,x,f=TraFlag[0]=='N' ? (TraFlag[1]=='N' ? 1 : 2) : (TraFlag[1]=='N' ? 3 : 4);

    //transpose MatC
    if ( TraFlag.size()>2 && TraFlag[2]=='T' ) {
        for (i=0; i<SizeA; i++) for (j=0; j<SizeB; j++) {
            LineAB=0.0;
            switch (f) {
                case 1: for (x=0; x<SizeAB; x++) LineAB+= MatA[x+i*SizeAB] *MatB[j+x*SizeB]; break;
                case 2: for (x=0; x<SizeAB; x++) LineAB+= MatA[x+i*SizeAB] *MatB[x+j*SizeAB]; break;
                case 3: for (x=0; x<SizeAB; x++) LineAB+= MatA[i+x*SizeA]  *MatB[j+x*SizeB]; break;
                case 4: for (x=0; x<SizeAB; x++) LineAB+= MatA[i+x*SizeA]  *MatB[x+j*SizeAB]; break;
            }
            if (CoeC==0.0) MatC[i+j*SizeA] = CoeAB*LineAB;
            else MatC[i+j*SizeA] = CoeAB*LineAB + CoeC*MatC[i+j*SizeA];
        }
    }
    else {
        for (i=0; i<SizeA; i++) for (j=0; j<SizeB; j++) {
            LineAB=0.0;
            switch (f) {
                case 1: for (x=0; x<SizeAB; x++) LineAB+= MatA[x+i*SizeAB] *MatB[j+x*SizeB]; break;
                case 2: for (x=0; x<SizeAB; x++) LineAB+= MatA[x+i*SizeAB] *MatB[x+j*SizeAB]; break;
                case 3: for (x=0; x<SizeAB; x++) LineAB+= MatA[i+x*SizeA]  *MatB[j+x*SizeB]; break;
                case 4: for (x=0; x<SizeAB; x++) LineAB+= MatA[i+x*SizeA]  *MatB[x+j*SizeAB]; break;
            }
            if (CoeC==0.0) MatC[i*SizeB+j] = CoeAB*LineAB;
            else MatC[i*SizeB+j] = CoeAB*LineAB + CoeC*MatC[i*SizeB+j];
        }
    }
}
/* (static) LU decomposition ------------------------------------------------------ */
int LUdcmp(vector<double> &MatSrc,int Size,vector<int> &index);
/* LU back-substitution ----------------------------------------------------------- */
void LUbksb(vector<double> &MatB,int Size,vector<int> &index,vector<double> &MatA,int startA);
/* inverse of matrix -------------------------------------------------------------- */
int matinv(vector<double> &MatSrc,int Size);
/* polynomial coefficients estimate ----------------------------------------------- */
int polyest(const vector<double> &A,const vector<double> &L,vector<double> &Coe,
    const int numL,const int numX,double *sigma);
/* solve linear equation ---------------------------------------------------------- */
int solve_line(const string TransFlag,const vector<double> &A,const vector<double> &Y,
    vector<double> &X,int Xnum,int SolNum);
/* get index ---------------------------------------------------------------------- */
int getindex(double value,const double *range);
/* get number of items ------------------------------------------------------------ */
int nitem(const double *range);

/* data format transfer functions ----------------------------------------------------------------- */
/* replace keywords in file path -------------------------------------------------- */
int reppath(const string SrcPath,string &DstPath,gtime_c GpsTime,
    const string RovID,const string BaseID);
/* convert double number to string ------------------------------------------------ */
string doul2str(int StrLen,int DecLen,const string StrFiller,const double SrcNum,string &DstStr);
/* convert int number to string --------------------------------------------------- */
string int2str(int StrLen,const string StrFiller,const int SrcNum,string &Dststr);
/* convert string to double number ------------------------------------------------ */
int str2double(const string SrcStr,double &DstNum);
/* convert string to int number --------------------------------------------------- */
int str2int(const string SrcStr,int &DstNum);
/* convert between matrix and vector, Iter1 to Iter2 ------------------------------ */
template <typename Iter1,typename Iter2>
int vecarr(const Iter1 VecArr1,Iter2 VecArr2,int SizeVA) {
    while (--SizeVA>=0) VecArr1[SizeVA] = VecArr2[SizeVA];
    return 1;
}

/* time and position transfer functions ----------------------------------------------------------- */
/* get tick time ------------------------------------------------------------------ */
unsigned int tickget();
/* sleep ms ----------------------------------------------------------------------- */
void sleepms(int ms);
/* adjust gps week number --------------------------------------------------------- */
int adjgpsweek(int UnWeek);
/* convert degree to deg-min-sec -----------------------------------------------------
* convert degree to degree-minute-second
* args   : double  Degree       I   degree
*          double *DMS          O   degree-minute-second {deg,min,sec}
*          int     NumDec       I   number of decimals of second
* return : none
*---------------------------------------------------------------------------------- */
template <typename Iter1>
void deg2dms(const double Degree,Iter1 DMS,int NumDec)
{
    double sign=Degree<0.0 ? -1.0 : 1.0,a=fabs(Degree);
    double unit=pow(0.1,NumDec);
    DMS[0]=floor(a); a=(a-DMS[0])*60.0;
    DMS[1]=floor(a); a=(a-DMS[1])*60.0;
    DMS[2]=floor(a/unit+0.5)*unit;
    if (DMS[2]>=60.0) {
        DMS[2]=0.0;
        DMS[1]+=1.0;
        if (DMS[1]>=60.0) {
            DMS[1]=0.0;
            DMS[0]+=1.0;
        }
    }
    DMS[0]*=sign;
}
/* convert deg-min-sec to degree -----------------------------------------------------
* convert degree-minute-second to degree
* args   : double *DMS      I   degree-minute-second {deg,min,sec}
* return : degree
*---------------------------------------------------------------------------------- */
template <typename Iter1>
double dms2deg(const Iter1 DMS)
{
    double sign=DMS[0]<0.0 ? -1.0 : 1.0;
    return sign*(fabs(DMS[0])+DMS[1]/60.0+DMS[2]/3600.0);
}

/* geography and astronomy functions -------------------------------------------------------------- */
/* convert vector between NED and ENU in navigation frame ------------------------- */
void cvt_NED_ENU(double VecTrg[3]);
/* convert vector between FRD and RFD in body frame ------------------------------- */
void cvt_FRD_RFU(double VecTrg[3]);
/* convert the sign of yaw in attitude angle -------------------------------------- */
void change_yaw(double att_n_b[3]);
/* using angle increment to update coordinate transformation matrix --------------- */
void angIncrem2Cno(const double AttIncrem[3],double C_new_old[9]);
/* using angle increment to compute average coordinate transformation matrix ------ */
void angIncrem2Cave(const double AttIncrem[3],const double ERAIncrem[3],
    const double C_old[9],double ave_C_new[9]);
/* calulates the meridian and transverse radii of curvature ----------------------- */
void radiiCurvature(const double Latitude,int DatumSys,double &R_E,double &R_N);
/* calulates angular rate in Nav frame -------------------------------------------- */
void angularRateNav(const double BlhPos[3],int DatumSys,const double v_eb_n[3],
    double omg_en_n[3],double *R_EN);
/* normalize coordinate transformation matrix ----------------------------------------
* compute coordinate transformation matrix (CTM) from the corresponding Euler angles
* args   : double *CTranMat    O   coordinate transformation matrix
* return : none
* note   : matirix stored by row-major order
*          for C_b_n, yaw is postive from north to west (anti-clockwise)
*---------------------------------------------------------------------------------- */
void normCTM(double CTranMat[9]);
/* Euler angles to the coordinate transformation matrix ------------------------------
* compute coordinate transformation matrix (CTM) from the corresponding Euler angles
* args   : double *Euler       I   Euler angles {pitch, roll, yaw} (rad)
*          double *CTranMat    O   corresponding coordinate transformation matrix
* return : none
* note   : matirix stored by row-major order
*          for C_b_n, yaw is postive from north to west (anti-clockwise)
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
void Euler2CTM(const Iter1 Euler,Iter2 CTranMat) {
    double sP=sin(Euler[0]), cP=cos(Euler[0]), sR=sin(Euler[1]), cR=cos(Euler[1]),
            sY=sin(Euler[2]), cY=cos(Euler[2]);

    CTranMat[0]= cY*cR-sY*sP*sR; CTranMat[1]=-sY*cP; CTranMat[2]= cY*sR+sY*sP*cR;
    CTranMat[3]= sY*cR+cY*sP*sR; CTranMat[4]= cY*cP; CTranMat[5]= sY*sR-cY*sP*cR;
    CTranMat[6]=-cP*sR;          CTranMat[7]= sP;    CTranMat[8]= cP*cR;
}
/* coordinate transformation matrix to the Euler angles ------------------------------
* compute Euler angles from the corresponding coordinate transformation matrix (CTM)
* args   : double *CTranMat    I   corresponding coordinate transformation matrix
*          double *Euler       O   Euler angles {pitch, roll, yaw} (rad)
* return : none
* note   : matirix stored by row-major order
*          for C_b_n, yaw is postive from north to west (anti-clockwise)
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
void CTM2Euler(const Iter1 CTranMat,Iter2 Euler) {
    if (fabs(CTranMat[7])<=0.999999) {
        Euler[0]= asin(CTranMat[7]); //pitch
        Euler[1]=-atan2(CTranMat[6],CTranMat[8]); //roll
        Euler[2]=-atan2(CTranMat[1],CTranMat[4]); //yaw
    }
    else {
        Euler[0]= asin(CTranMat[7]); //pitch
        Euler[1]= atan2(CTranMat[2],CTranMat[0]); //roll
        Euler[2]= 0; //yaw
    }
}
/* ecef (BLH) to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transformation matrix
* args   : double *BlhPos      I   geodetic position {lat,lon} (rad)
*          double *CTranMat    O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by row-major order
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
void blh2Cen(const Iter1 BlhPos,Iter2 CTranMat) {
    double sinLat=sin(BlhPos[0]),cosLat=cos(BlhPos[0]),sinLon=sin(BlhPos[1]),cosLon=cos(BlhPos[1]);

    CTranMat[0]=-sinLon;        CTranMat[1]= cosLon;        CTranMat[2]=0.0;
    CTranMat[3]=-sinLat*cosLon; CTranMat[4]=-sinLat*sinLon; CTranMat[5]=cosLat;
    CTranMat[6]= cosLat*cosLon; CTranMat[7]= cosLat*sinLon; CTranMat[8]=sinLat;
}
/* transform ecef to geodetic postion ------------------------------------------------
* transform ecef position to geodetic position
* args   : double *XyzPos      I   ecef position {x,y,z} (m)
*          int     DatumSys    I   datum {0:WGS84,1:CGCS2000}
*          double *BlhPos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
Iter1/Iter2 is a double pointer or a vector<double> iterator
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
void ecef2blh(const Iter1 XyzPos,int DatumSys,Iter2 BlhPos)
{
    double re,fe;
    if (DatumSys==CGCS2000) { re=RE_CGCS2000; fe=FE_CGCS2000; }
    else { re=RE_WGS84; fe=FE_WGS84; }
    double e2=fe*(2.0-fe),r2=dot(XyzPos,XyzPos,2),z,zk,v=re,sinp;

    for (z=XyzPos[2],zk=0.0; fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=re/sqrt(1.0-e2*sinp*sinp);
        z=XyzPos[2]+v*e2*sinp;
    }
    BlhPos[0]=r2>1E-12 ? atan(z/sqrt(r2)) : (XyzPos[2]>0.0 ? PI/2.0 : -PI/2.0);
    BlhPos[1]=r2>1E-12 ? atan2(XyzPos[1],XyzPos[0]) : 0.0;
    BlhPos[2]=sqrt(r2+z*z)-v;
}
/* ecef (XYZ) to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transformation matrix
* args   : double *XyzPos      I   geodetic position {x,y,z} (m)
*          int     DatumSys    I   datum {0:WGS84,1:CGCS2000}
*          double *CTranMat    O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by row-major order
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
void xyz2Cen(const Iter1 XyzPos,int DatumSys,Iter2 CTranMat) {
    double BlhPos[3];
    ecef2blh(XyzPos,DatumSys,BlhPos);

    double sinLat=sin(BlhPos[0]),cosLat=cos(BlhPos[0]),sinLon=sin(BlhPos[1]),cosLon=cos(BlhPos[1]);

    CTranMat[0]=-sinLon;        CTranMat[1]= cosLon;        CTranMat[2]=0.0;
    CTranMat[3]=-sinLat*cosLon; CTranMat[4]=-sinLat*sinLon; CTranMat[5]=cosLat;
    CTranMat[6]= cosLat*cosLon; CTranMat[7]= cosLat*sinLon; CTranMat[8]=sinLat;
}
/* transform ecef vector to local tangental coordinate ---------------------------- */
template <typename Iter1,typename Iter2,typename Iter3>
void ecef2enu(const Iter1 BlhPos,const Iter2 SightVec,Iter3 EnuPos) {
    double C_e_n[9];

    blh2Cen(BlhPos,C_e_n);
    matmul("NN",3,1,3,1.0,C_e_n,SightVec,0.0,EnuPos);
}
/* calculate the difference in ecef and transform to local tangental coordinate ------
* dxyz = XyzPos2 - XyzPos1, EnuPos = C_e_n * dxyz */
template <typename Iter1,typename Iter2,typename Iter3>
void ecef2denu(const Iter1 XyzPos1,const Iter2 XyzPos2,int DatumSys,Iter3 EnuPos) {
    double C_e_n[9],dxyz[3];

    dxyz[0]=XyzPos2[0]-XyzPos1[0];
    dxyz[1]=XyzPos2[1]-XyzPos1[1];
    dxyz[2]=XyzPos2[2]-XyzPos1[2];
    xyz2Cen(XyzPos1,DatumSys,C_e_n);
    matmul("NN",3,1,3,1.0,C_e_n,dxyz,0.0,EnuPos);
}
template <typename Iter1>
void sssssssssss(Iter1 &XyzPos1) {
    double x=10,y=XyzPos1[0];
    x=y;
}
/* transform local vector to ecef coordinate -----------------------------------------
* transform local tangental coordinate vector to ecef
* args   : double *BlhPos      I   geodetic position {lat,lon} (rad)
*          double *EnuPos      I   vector in local tangental coordinate {e,n,u}
*          double *XyzVec      O   difference vector in ecef coordinate {x,y,z}
* return : none
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2,typename Iter3>
void enu2ecef(const Iter1 BlhPos,const Iter2 EnuPos,Iter3 XyzVec) {
    double C_e_n[9];

    blh2Cen(BlhPos,C_e_n);
    matmul("TN",3,1,3,1.0,C_e_n,EnuPos,0.0,XyzVec);
}
/* transform covariance to local tangental coordinate --------------------------------
* transform ecef covariance to local tangental coordinate
* args   : double *BlhPos      I   geodetic position {lat,lon} (rad)
*          double *XyzVar      I   covariance in ecef coordinate
*          double *BlhVar      O   covariance in local tangental coordinate
* return : none
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2,typename Iter3>
void covenu(const Iter1 BlhPos,const Iter2 XyzVar,Iter3 BlhVar)
{
    double C_e_n[9],EP[9];

    blh2Cen(BlhPos,C_e_n);
    matmul("NN",3,3,3,1.0,C_e_n,XyzVar,0.0,EP);
    matmul("NT",3,3,3,1.0,EP,C_e_n,0.0,BlhVar);
}
/* transform covariance to local tangental coordinate --------------------------------
* transform ecef covariance to local tangental coordinate
* args   : double *BlhPos      I   geodetic position {lat,lon} (rad)
*          double *BlhVar      I   covariance in local tangental coordinate
*          double *XyzVar      O   covariance in xyz coordinate
* return : none
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2,typename Iter3>
void covecef(const Iter1 BlhPos,const Iter2 BlhVar,Iter3 XyzVar)
{
    double C_e_n[9],EQ[9];

    blh2Cen(BlhPos,C_e_n);
    matmul("TN",3,3,3,1.0,C_e_n,BlhVar,0.0,EQ);
    matmul("NN",3,3,3,1.0,EQ,C_e_n,0.0,XyzVar);
}
/* transform geodetic to ecef position -----------------------------------------------
* transform geodetic position to ecef position
* args   : double *BlhPos      I   geodetic position {lat,lon,h} (rad,m)
*          int     DatumSys    I   datum {0:WGS84,1:CGCS2000}
*          double *XyzPos      O   ecef position {x,y,z} (m)
* return : none
* notes  : WGS84, ellipsoidal height
Iter1/Iter2 is a double pointer or a vector<double> iterator
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
void blh2ecef(const Iter1 BlhPos,int DatumSys,Iter2 XyzPos)
{
    double re,fe;
    if (DatumSys==CGCS2000) { re=RE_CGCS2000; fe=FE_CGCS2000; }
    else { re=RE_WGS84; fe=FE_WGS84; }
    double sinLat=sin(BlhPos[0]),cosLat=cos(BlhPos[0]),sinLon=sin(BlhPos[1]),cosLon=cos(BlhPos[1]);
    double e2=fe*(2.0-fe),v=re/sqrt(1.0-e2*sinLat*sinLat);

    XyzPos[0]=(v+BlhPos[2])*cosLat*cosLon;
    XyzPos[1]=(v+BlhPos[2])*cosLat*sinLon;
    XyzPos[2]=(v*(1.0-e2)+BlhPos[2])*sinLat;
}
/* geometric distance on the earth ------------------------------------------------ */
template <typename Iter1,typename Iter2,typename Iter3>
double geodist_earth(const Iter1 Pos1,const Iter2 Pos2,Iter3 SightVec) {
    double r;

    /* original satellite position */
    SightVec[0]=Pos1[0]-Pos2[0];
    SightVec[1]=Pos1[1]-Pos2[1];
    SightVec[2]=Pos1[2]-Pos2[2];
    r=norm(SightVec,3);

    SightVec[0]/=r; SightVec[1]/=r; SightVec[2]/=r;

    return r;
}
/* geometric distance for GNSS ---------------------------------------------------- */
template <typename Iter1,typename Iter2,typename Iter3>
double geodist_GNSS(const Iter1 Pos1,const Iter2 Pos2,Iter3 SightVec) {
    double r;

    /* original satellite position */
    SightVec[0]=Pos1[0]-Pos2[0];
    SightVec[1]=Pos1[1]-Pos2[1];
    SightVec[2]=Pos1[2]-Pos2[2];
    r=norm(SightVec,3);

    double dr=OMGE*(Pos1[0]*Pos2[1] - Pos1[1]*Pos2[0])/CLIGHT;

    SightVec[0]/=r; SightVec[1]/=r; SightVec[2]/=r;

    return r+dr;
}
/* satellite azimuth/elevation angle ---------------------------------------------- */
template <typename Iter1,typename Iter2,typename Iter3>
double satazel(const Iter1 BlhPos,const Iter2 SightVec,Iter3 AzEl) {
    double az=0.0,el=PI/2.0,enu[3];

    if ((BlhPos[2])>-RE_WGS84) {
        ecef2enu(BlhPos,SightVec,enu);
        az=dot(enu,enu,2)<1E-12 ? 0.0 : atan2(enu[0],enu[1]);
        if (az<0.0) az+=2*PI;
        el=asin(enu[2]);
    }
    AzEl[0]=az; AzEl[1]=el;
    return el;
}
/* astronomical arguments: f={l,l',F,D,OMG} (rad) --------------------------------- */
void ast_args(double t,double *f);
/* eci to ecef transformation matrix ---------------------------------------------- */
void eci2ecef(gtime_c tutc,const double *erpv,double *U,double *gmst);
/* sun and moon position ---------------------------------------------------------- */
void sunmoonpos(gtime_c UT1Time,double *ERPValue,double *SunPos,
    double *MoonPos,double *gmst);
/* get earth rotation parameter values -------------------------------------------- */
int geterpv(const gnss_erp_c *Earth_Par,gtime_c GPS_Time,double *Earth_Value);
/* get 6 order earth gravity in Nav --------------------------------------------------
* get 6 order earth gravity in Nav
* args   : double *BlhPos      I   geodetic position {lat,lon,h} (rad,m)
*          double *G_Nav       O   gravity in ecef  (m/s^2)
* return : 1:ok 0:error position
Iter1/Iter2 is a double pointer or a vector<double> iterator
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
int g6_in_nav(const Iter1 BlhPos,Iter2 G_Nav) {
    double sinL2 = sin(BlhPos[0])*sin(BlhPos[0]);
    double gg0 = GRAVITY_0 * (1 + 0.001931853 * sinL2) / sqrt(1 - EE_WGS84*EE_WGS84 * sinL2);
    double RE2=RE_WGS84*RE_WGS84;

    G_Nav[0] = 0;
    G_Nav[1] = -8.08E-9 * BlhPos[2] * sin(2 * BlhPos[0]);
    G_Nav[2] = -gg0 * (1 - (2 / RE_WGS84) * (1 + FE_WGS84 * (1 - 2 * sinL2) +
        (OMGE*OMGE * RE2 * RP_WGS84 / MU_ALL)) * BlhPos[2] + (3 * BlhPos[2]*BlhPos[2] / RE2));

    return 1;
}
/* get 6 order earth gravity in ECEF -------------------------------------------------
* get 6 order earth gravity in ECEF
* args   : double *XyzPos      I   ecef position {x,y,z} (m)
*          double *G_ECEF      O   gravity in ecef  (m/s^2)
* return : 1:ok 0:error position
Iter1/Iter2 is a double pointer or a vector<double> iterator
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
int g6_in_ecef(const Iter1 XyzPos,Iter2 G_ECEF) {
    double r=norm(XyzPos,3),r2=r*r,OMG2=OMGE*OMGE; //radius
    if (r-RE_WGS84<-1E5) return 0;
    double p=RE_WGS84/r,p2=p*p,p4=p2*p2,p6=p4*p2,
        t=XyzPos[2]/r,t2=t*t,t4=t2*t2,t6=t4*t2; //t=z/r
    //basic coefficients
    double a1r=-MU_ALL/r2/r,a2=1+GA22*J2*p2+GA24*J4*p4+GA26*J6*p6,
        a3=GA32*J2*p2+GA34*J4*p4+GA36*J6*p6,a4=GA44*J4*p4+GA46*J6*p6,a5=GA56*J6*p6,
        b1=GB12*J2*p2+GB14*J4*p4+GB16*J6*p6,b2=GB24*J4*p4+GB26*J6*p6,b3=GB36*J6*p6;
    //gravity coefficients
    double c1=a2,c2=a3-b1,c3=a4-b2,c4=a5-b3,
        d1=a2+b1,d2=c2+b2,d3=c3+b3,d4=c4,
        cXY=c1+c2*t2+c3*t4+c4*t6,cZ=d1+d2*t2+d3*t4+d4*t6;

    G_ECEF[0]=(a1r*cXY+OMG2)*XyzPos[0]; //in X
    G_ECEF[1]=(a1r*cXY+OMG2)*XyzPos[1]; //in Y
    G_ECEF[2]=(a1r*cZ)*XyzPos[2]; //in Z
    return 1;
}
/* get earth gravity error coefficients in ECEF --------------------------------------
* get earth gravity error coefficients in ECEF
* args   : double *XyzPos      I   ecef position {x,y,z} (m)
*          double *Gcoe_ECEF   O   gravity in ecef [3*3] (s^-2) 
* return : 1:ok 0:error position
Iter1/Iter2 is a double pointer or a vector<double> iterator
*---------------------------------------------------------------------------------- */
template <typename Iter1,typename Iter2>
int gcoe_in_ecef(const Iter1 XyzPos,Iter2 Gcoe_ECEF) {
    double r=norm(XyzPos,3),r2=r*r,kMr3=MU_ALL/r2/r,OMG2=OMGE*OMGE; //radius
    if (r-RE_WGS84<-1E5) return 0;
    double x2=XyzPos[0]*XyzPos[0],y2=XyzPos[1]*XyzPos[1],z2=XyzPos[2]*XyzPos[2],
        xy=XyzPos[0]*XyzPos[1],yz=XyzPos[1]*XyzPos[2],zx=XyzPos[2]*XyzPos[0];
    Gcoe_ECEF[0]=kMr3*(3*x2/r2-1)+OMG2; //Gcoe_ECEF[0,0]
    Gcoe_ECEF[1]=kMr3*(3*xy/r2); //Gcoe_ECEF[0,1]
    Gcoe_ECEF[2]=kMr3*(3*zx/r2); //Gcoe_ECEF[0,2]
    Gcoe_ECEF[3]=kMr3*(3*xy/r2); //Gcoe_ECEF[1,0]
    Gcoe_ECEF[4]=kMr3*(3*y2/r2-1)+OMG2; //Gcoe_ECEF[1,1]
    Gcoe_ECEF[5]=kMr3*(3*yz/r2); //Gcoe_ECEF[1,2]
    Gcoe_ECEF[6]=kMr3*(3*zx/r2); //Gcoe_ECEF[2,0]
    Gcoe_ECEF[7]=kMr3*(3*yz/r2); //Gcoe_ECEF[2,1]
    Gcoe_ECEF[8]=kMr3*(3*z2/r2-1); //Gcoe_ECEF[2,2]
    return 0;
}

/* satellite data functions ----------------------------------------------------------------------- */
/* satellite system+prn/slot number to satellite number --------------------------- */
int satno(int StaSys,int PrnNum);
/* satellite number to satellite system ------------------------------------------- */
int satsys(int SatNum,int *PrnNum);
/* system code to system number (0 -- NSYS-1) ------------------------------------- */
int syscd2num(int sys);
/* satellite id to satellite number ----------------------------------------------- */
int satid2no(string SatID);
/* satellite number to satellite id ----------------------------------------------- */
int satno2id(int sat,string &SatID);
/* used satellite system to string format flag (GCRE) ----------------------------- */
string satsys_flag(const int navsys);

/* observation and code transfer functions -------------------------------------------------------- */
/* obs type string to obs code ---------------------------------------------------- */
int obs2code(string ObsCode,int *ObsFre);
/* obs code to obs code string ---------------------------------------------------- */
string code2obs(unsigned char ObsCode,int *ObsFre);
/* satellite code to satellite system --------------------------------------------- */
int code2sys(char SysCode);
/* get code priority -------------------------------------------------------------- */
int getcodepri(int SatSys,unsigned char ObsCode,string CodeOpt);
/* add fatal callback function ---------------------------------------------------- */
void add_fatal(fatalfunc_t *func);

/* GPS data functions ----------------------------------------------------------------------------- */
/* arrange observation data ------------------------------------------------------- */
int sortobs(gnss_obs_c &SrcObs);
/* station-cross single-difference observation ------------------------------------ */
double single_diff(const double *ObsRov,const double *ObsBas);
/* compute ionosphere-free combination -------------------------------------------- */
double iono_free(const double &ObsF1,const double &ObsF2,const double &lam1,const double &lam2);
/* compute geometry-free combination ---------------------------------------------- */
double geometry_free(const double &ObsF1,const double &ObsF2);
/* compute Melbourne-Wubbena combination ------------------------------------------ */
double Mel_Wub(const double &ObsP1,const double &ObsP2,const double &ObsL1,const double &ObsL2,
    const double &lam1,const double &lam2);
/* compute narrow-lane ambiguity -------------------------------------------------- */
double Narrow(const double &ObsP1,const double &ObsP2,const double &ObsL1,const double &ObsL2,
    const double &lam1,const double &lam2);
/* compute ambiguity combination -------------------------------------------------- */
double amb_cmb(const double Lr,const double P);
/* compute single time-cross difference ------------------------------------------- */
double single_time(const double obs_t1,const double obs_t2);
/* get tgd parameter (m) ---------------------------------------------------------- */
int gettgd(gnss_obsd_c *obs,const gnss_nav_c *nav,double tgd[2]);
/* correct code observation with obsolute bias ------------------------------------ */
int correct_code_bias(gnss_obsd_c *obs,const gnss_nav_c *nav,const gnss_prcopt_c *opt);
/* correct code observation with navigation tgd ----------------------------------- */
int correct_code_tgd(gnss_obsd_c *obs,const gnss_nav_c *nav,const gnss_prcopt_c *opt);

/* system functions ------------------------------------------------------------------------------- */
/* execute command ---------------------------------------------------------------- */
int execcmd(const string StrCmd);
/* create directory --------------------------------------------------------------- */
void createdir(const string StrPath);
/* uncompress file ---------------------------------------------------------------- */
int rtk_uncompress(const string SrcFile,string UncFile);

#endif
