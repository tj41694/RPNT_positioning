#include "Adjustment/adjustment.h"
#include "BaseFunction/timesys.h"
#include "BaseFunction/basefunction.h"

/* const ------------------------------------------------------------------------------------------ */
#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?sqrt(fabs(x)):sqrt(x))

/* robust constants ------------------------------------------------------------------------------- */
#define ROBUST_MAX_ITERATION    3       /* robust maximum iteration number */
#define ROBUST_RES_FACTOR       1.483   /* IGG3 robust residual factor */
#define ROBUST_CUTOFF           10000.0 /* cutoff of expansion factor in robust function */
#define ROBUST_CUTOFF_INVERSE   0.0001  /* inverse of ROBUST_CUTOFF */
#define ROBUST_CUTOFF_SQRT      100.0   /* sqrt(ROBUST_CUTOFF) */
#define ROBUST_S2_K0            1.5     /* default 2-segment IGG3 robust threshold k0 */
#define ROBUST_S3_K0            1.5     /* default 3-segment IGG3 robust threshold k0 */
#define ROBUST_S3_K1            6.0     /* default 3-segment IGG3 robust threshold k1 */

/* least-square constants ------------------------------------------------------------------------- */
#define LSQ_FINISH_DX           1E-4    /* maximum dx to finish least-square adjustment */
#define LSQ_CONVERGE_DX         100     /* maximum dx to determine least-square converged */

/* Helmert constants ------------------------------------------------------------------------------ */
const int    MAX_ITER=3;				/* maximum iteration number of Helmert algorithm */
const double RXL_THRES=0.3;				/* Rx sigma lower threshold */
const double RXH_THRES=3.0;				/* Rx sigma higher threshold */
const double RLOW_THRES=0.01;			/* R sigma lower threshold */
const double RHIG_THRES=100.0;			/* R sigma higher threshold */
const double MAX_DSGM=0.2;				/* maximum difference of sgm with 1.0 */
const double SGM_THRES=0.07;			/* sigma threshold to stop */
const double BIG_SGM=10.0;				/* replace with this value if sgm is negative */
const string SYS_LIST="GREC";			/* system list of GPS, GLONASS, Galileo, BeiDou */
const string PC_LIST="PC";				/* observation type list of P:phase, C:code */

/* parent adjustment functions -------------------------------------------------------------------- */
/* Constructor -------------------------------------------------------------------- */
adjfunc_c::adjfunc_c() {
    flag=0;
    sigma0_2=0;
    sig0=0;
    H_niter=0; H_flag=1;
    iteraRobust=0;
    robustFcp[0][0]=&adjfunc_c::robust_empty;
    robustFcp[0][1]=&adjfunc_c::robust_yang_IGG3_2s_Rdiag;
    robustFcp[0][2]=&adjfunc_c::robust_yang_IGG3_3s_Rdiag;
    robustFcp[1][0]=&adjfunc_c::robust_empty;
    robustFcp[1][1]=&adjfunc_c::robust_yang_IGG3_2s;
    robustFcp[1][2]=&adjfunc_c::robust_yang_IGG3_3s;
    robustFcp[2][0]=&adjfunc_c::robust_empty;
    robustFcp[2][1]=&adjfunc_c::robust_yang_IGG3_2s_Rvec;
    robustFcp[2][2]=&adjfunc_c::robust_yang_IGG3_3s_Rvec;
    robustMode=ADJUST_ROBUST_NONE;
    robustRtype=ROBUST_R_DIAGONAL;
    adjustFcp=&adjfunc_c::LSQ;
    robust2sK0=ROBUST_S2_K0;
    robust3sK0=ROBUST_S3_K0;
    robust3sK1=ROBUST_S3_K1;
    robust3sK1_1=1/ROBUST_S3_K1;
}adjfunc_c::adjfunc_c(int Robust,int RobustCovar) {
    flag=0;
    sigma0_2=0;
    sig0=0;
    H_niter=0; H_flag=1;
    iteraRobust=0;
    robustFcp[0][0]=&adjfunc_c::robust_empty;
    robustFcp[0][1]=&adjfunc_c::robust_yang_IGG3_2s_Rdiag;
    robustFcp[0][2]=&adjfunc_c::robust_yang_IGG3_3s_Rdiag;
    robustFcp[1][0]=&adjfunc_c::robust_empty;
    robustFcp[1][1]=&adjfunc_c::robust_yang_IGG3_2s;
    robustFcp[1][2]=&adjfunc_c::robust_yang_IGG3_3s;
    robustFcp[2][0]=&adjfunc_c::robust_empty;
    robustFcp[2][1]=&adjfunc_c::robust_yang_IGG3_2s_Rvec;
    robustFcp[2][2]=&adjfunc_c::robust_yang_IGG3_3s_Rvec;
    robustMode=Robust;
    robustRtype=RobustCovar;
    if (robustMode>ADJUST_ROBUST_NONE) {
        if (robustRtype==ROBUST_R_VECTOR) adjustFcp=&adjfunc_c::robust_LSQ_Rvec;
        else adjustFcp=&adjfunc_c::robust_LSQ;
    }
    else {
        if (robustRtype==ROBUST_R_VECTOR) adjustFcp=&adjfunc_c::LSQ_Rvec;
        else adjustFcp=&adjfunc_c::LSQ;
    }
    robust2sK0=ROBUST_S2_K0;
    robust3sK0=ROBUST_S3_K0;
    robust3sK1=ROBUST_S3_K1;
    robust3sK1_1=1/ROBUST_S3_K1;
}
adjfunc_c::~adjfunc_c() {
    P.clear(); AP.clear(); APL.clear(); LP.clear();
    RxA.clear(); Q.clear(); E.clear(); ERx.clear(); KR.clear();
    helmertR.clear(); newRx.clear(); buffRx.clear();
    sgm.clear(); fsgm.clear();
    sgm_all.clear();

    V.clear(); dX.clear(); K.clear(); QQQ.clear();
    normV.clear(); sortV.clear();
    robustFactor.clear();
    sgm2.clear(); pre_fsgm.clear();
    sys_flag.clear(); neg_flag.clear();
}
/* Implementaion functions -------------------------------------------------------- */
/* least-square functions --------------------------------------------------------- */
/* least-square: estimate dx ------------------------------------------------------ */
int adjfunc_c::LSQ_est_dx(const vector<double> &A,const vector<double> &L,vector<double> &R,
    int numL,int numX) {

    matcpy(P.begin(),R.begin(),numL,numL);//weight matrix P
    if (matinv(P,numL)==0) { //P=R^-1
        /* compute AP,APA */
        matmul_vec("TN",numX,numL,numL,1.0,A,P,0.0,AP); //AP
        matmul_vec("NN",numX,1,numL,1.0,AP,L,0.0,APL);  //APL
        /* cofactor matrix Q */
        matmul_vec("NN",numX,numX,numL,1.0,AP,A,0.0,QQQ);
        if (matinv(QQQ,numX)==0) {
            /* compute dX */
            matmul_vec("NN",numX,1,numX,1.0,QQQ,APL,0.0,dX);  //dX
            /* residual of L */
            V.assign(L.begin(),L.end());
            matmul_vec("NN",numL,1,numX,-1.0,A,dX,1.0,V);
            return flag=0;
        }
        else return flag=-1;
    }
    else return flag=-1;
}
/* least-square: estimate dx for variance vector Rvec ----------------------------- */
int adjfunc_c::LSQ_est_dx_Rvec(const vector<double> &A,const vector<double> &L,
    vector<double> &Rvec,int numL,int numX) {
    /* get AP */
    for (int i=0; i<numX; i++) for (int j=0; j<numL; j++) {
        AP[i*numL+j] = A[j*numX+i]/Rvec[j];
    }
    /* compute APA */
    matmul_vec("NN",numX,1,numL,1.0,AP,L,0.0,APL); //APL
    /* cofactor matrix Q */
    matmul_vec("NN",numX,numX,numL,1.0,AP,A,0.0,QQQ);
    if (matinv(QQQ,numX)==0) {
        /* compute dX */
        matmul_vec("NN",numX,1,numX,1.0,QQQ,APL,0.0,dX); //dX
        /* residual of L */
        V.assign(L.begin(),L.end());
        matmul_vec("NN",numL,1,numX,-1.0,A,dX,1.0,V);
        return flag=0;
    }
    else return flag=-1;
}
/* least-squre: update x and Rx --------------------------------------------------- */
int adjfunc_c::LSQ_update_xRx(const vector<double> &L,vector<double> &X,
    vector<double> &Rx,int numL,int numX) {
    /* update X and Rx */
    for (int i=0; i<numX; i++) X[i]+=dX[i];
    /* compute Rx */
    LP.assign(numL,1);
    if (numL>numX) {
        matmul_vec("TN",1,numL,numL,1.0,L,P,0.0,LP); //LP
        matmul("NN",1,1,numL,1.0,LP.begin(),L.begin(),0.0,&VPV); //LPL to VPV
        matmul("TN",1,1,numX,-1.0,APL.begin(),dX.begin(),1.0,&VPV); //VPV=LPL-(APL)X
        //coe=VPV[0];
        sigma0_2=VPV/(numL-numX);
    }
    Rx.assign(QQQ.begin(),QQQ.end());
    for (size_t i=0; i<Rx.size(); i++) Rx[i]=Rx[i];

    double normDx=norm(dX.begin(),numX);
    if (normDx>LSQ_CONVERGE_DX) return flag=0;
    else if (normDx>LSQ_FINISH_DX) return flag=1;
    else return flag=100;
}
/* least-squre: update x and Rx for variance vector Rvec -------------------------- */
int adjfunc_c::LSQ_update_xRx_Rvec(const vector<double> &L,vector<double> &Rvec,
    vector<double> &X,vector<double> &Rx,int numL,int numX) {
    /* update X and Rx */
    for (int i=0; i<numX; i++) X[i]+=dX[i];
    /* compute Rx */
    if (numL>numX) {
        VPV=0;
        for (int i=0; i<numX; i++) VPV+=L[i]*L[i]/Rvec[i]; //VPV
        matmul("TN",1,1,numX,-1.0,APL.begin(),dX.begin(),1.0,&VPV); //VPV=LPL-(APL)X
        //coe=VPV[0];
        sigma0_2=VPV/(numL-numX);
    }
    Rx.assign(QQQ.begin(),QQQ.end());
    for (int i=0; i<Rx.size(); i++) Rx[i]=Rx[i];

    double normDx=norm(dX.begin(),numX);
    if (normDx>LSQ_CONVERGE_DX) return flag=0;
    else if (normDx>LSQ_FINISH_DX) return flag=1;
    else return flag=100;
}
/* kalman filter functions -------------------------------------------------------- */
/* kalman filter : estimate dx ---------------------------------------------------- */
int adjfunc_c::KF_est_dx(const vector<double> &A,const vector<double> &L,vector<double> &R,
    vector<double> &Rx,int numL,int numX) {

    matcpy(Q.begin(),R.begin(),numL,numL);
    matmul_vec("NT",numX,numL,numX,1.0,Rx,A,0.0,RxA); // Rx*A'
    matmul_vec("NN",numL,numL,numX,1.0,A,RxA,1.0,Q); //Q=A*RxA'+R
    if (matinv(Q,numL)==0) {
        matmul_vec("NN",numX,numL,numL,1.0,RxA,Q,0.0,K); // K=RxA'*Q^-1
        matmul_vec("NN",numX,1,numL,1.0,K,L,0.0,dX); // dX=K*L
        /* residual of L */
        V.assign(L.begin(),L.end());
        matmul_vec("NN",numL,1,numX,-1.0,A,dX,1.0,V);
        return flag=0;
    }
    return flag=-1;
}
/* kalman filter: estimate dx for variance vector Rvec ---------------------------- */
int adjfunc_c::KF_est_dx_Rvec(const vector<double> &A,const vector<double> &L,
    vector<double> &Rvec,vector<double> &Rx,int numL,int numX) {

    matmul_vec("NT",numX,numL,numX,1.0,Rx,A,0.0,RxA); //Rx*A'
    matmul_vec("NN",numL,numL,numX,1.0,A,RxA,0.0,Q); //Q=A*RxA'
    for (int i=0; i<numL; i++) Q[i*numL+i]+=Rvec[i]; //Q=A*RxA'+R
    if (matinv(Q,numL)==0) {
        matmul_vec("NN",numX,numL,numL,1.0,RxA,Q,0.0,K); //K=RxA'*Q^-1
        matmul_vec("NN",numX,1,numL,1.0,K,L,0.0,dX); //dX=K*L
        /* residual of L */
        V.assign(L.begin(),L.end());
        matmul_vec("NN",numL,1,numX,-1.0,A,dX,1.0,V);
        return flag=0;
    }
    return flag=-1;
}
/* kalman filter: update x and Rx ------------------------------------------------- */
int adjfunc_c::KF_update_xRx(const vector<double> &A,const vector<double> &L,vector<double> &R,
    vector<double> &X,vector<double> &Rx,int numL,int numX) {

    matmul_vec("NN",numX,numX,numL,-1.0,K,A,1.0,E); //E=E-K*A
    //matmul_vec("NN",numX,numX,numX,1.0,E,Rx,0.0,QQQ); //QQQ=E*Rx
    //equivalent QQQ = E*RxE'+ KRK'
    matmul_vec("NN",numX,numX,numX,1.0,E,Rx,0.0,ERx); //ERx=E*Rx
    matmul_vec("NT",numX,numX,numX,1.0,ERx,E,0.0,QQQ); //QQQ=E*RxE'
    matmul_vec("NN",numX,numL,numL,1.0,K,R,0.0,KR); //KR=K*R
    matmul_vec("NT",numX,numX,numL,1.0,KR,K,1.0,QQQ); //QQQ=E*RxE'+KRK'
    /* update X */
    for (int i=0; i<numX; i++) X[i]+=dX[i];
    /* update Rx */
    matcpy(Rx.begin(),QQQ.begin(),numX,numX);

    return flag;
}
/* kalman filter: update x and Rx for variance vector Rvec ------------------------ */
int adjfunc_c::KF_update_xRx_Rvec(const vector<double> &A,const vector<double> &L,
    vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX){

    matmul_vec("NN",numX,numX,numL,-1.0,K,A,1.0,E); //E=E-K*A
    //matmul_vec("NN",numX,numX,numX,1.0,E,Rx,0.0,QQQ); //QQQ=E*Rx
    //equivalent QQQ = E*RxE'+ KRK'
    matmul_vec("NN",numX,numX,numX,1.0,E,Rx,0.0,ERx); //ERx=E*Rx
    matmul_vec("NT",numX,numX,numX,1.0,ERx,E,0.0,QQQ); //QQQ=E*RxE'
    for (int i=0; i<numX; i++) for (int j=0; j<numL; j++) { //KR=K*R
        KR[i*numL+j] = K[i*numL+j]*Rvec[j];
    }
    matmul_vec("NT",numX,numX,numL,1.0,KR,K,1.0,QQQ); //QQQ=E*RxE'+KRK'
    /* update X */
    for (int i=0; i<numX; i++) X[i]+=dX[i];
    /* update Rx */
    matcpy(Rx.begin(),QQQ.begin(),numX,numX);

    return flag;
}
/* Helmert component covariance estimate for multi-GNSS observation (tight) ----------
* formula : S * sgm2 = W
* base matrix:
* P   = R^-1
* Px  = Qx = Rx^-1
* for the ith system :
* Sii = ni - 2tr(Qi*QQQ^-1) + tr(Qi*QQQ^-1*Qi*QQQ^-1)
* Sij = tr(Qi*QQQ^-1*Qj*QQQ^-1) (j!=i)
* Wi  = Vi*Pi*Vi
* --------------------------------------------------------------------------------- */
int adjfunc_c::helmert_tight(const vector<double> &A,int numL,int numX,int numF,
    const int nobs[NSYS][NFREQ*2]) {
    /* initialization */
    vector<double> PP(helmertR.begin(),helmertR.end()),Px(newRx.begin(),newRx.end()),
        Qi[NSYS+1],Pi,VPi,Ai,APi,S,W;
    /* compute P */
    if (matinv(PP,numL)==-1||matinv(Px,numX)==-1) return -1;

    sys_flag.clear();
    for (int i=0; i<NSYS+1; i++) sys_n[i]=0;
    /* loop of satellite systems to get Qi and W */
    for (int sys=nsys=0,nnow=0; sys<NSYS; sys++) {
        for (int fff=0; fff<numF*2; fff++) sys_n[nsys]+=nobs[sys][fff];
        if (sys_n[nsys]==0) continue;

        /* initialize vectors */
        VPi.assign(sys_n[nsys],0.0); Pi.assign(sys_n[nsys]*sys_n[nsys],0.0);
        APi.assign(numX*sys_n[nsys],0.0);
        Ai.assign(A.begin()+nnow*numX,A.begin()+(nnow+sys_n[nsys])*numX); //Ai
        Qi[nsys].assign(numX*numX,0.0); W.push_back(0.0);

        for (int i=0; i<sys_n[nsys]; i++) for (int j=0; j<sys_n[nsys]; j++)
            Pi[j + i*sys_n[nsys]]=PP[j+nnow + (i+nnow)*numL]; // Pi
        matmul_vec("TN",1,sys_n[nsys],sys_n[nsys],1.0,V,Pi,0.0,VPi,nnow,0,0); //VPi
        matmul_vec("NN",1,1,sys_n[nsys],1.0,VPi,V,0.0,W,0,nnow,nsys); //W[nsys]

        matmul_vec("TN",numX,sys_n[nsys],sys_n[nsys],1.0,Ai,Pi,0.0,APi); //APi
        matmul_vec("NN",numX,numX,sys_n[nsys],1.0,APi,Ai,0.0,Qi[nsys]); //Qi[nsys]

        /* update system information */
        sys_flag.push_back("A"+SYS_LIST.substr(sys,1)); //A+(GREC)
        nnow+=sys_n[nsys]; nsys++;
    }
    if (nsys<0) return -1;
    /* Qi and W of X parameter */
    sys_n[nsys]=numX; sys_flag.push_back("X");
    int nsys_all=nsys+1; vector<double> dXP(numX,0.0);
    Qi[nsys].assign(Px.begin(),Px.end()); //Qx
    W.push_back(0.0);
    matmul_vec("TN",1,numX,numX,1.0,dX,Px,0.0,dXP); //dXP
    matmul_vec("NN",1,1,numX,1.0,dXP,dX,0.0,W,0,0,nsys); //Wx

    /* initilaize sgm2 and S */
    sgm.assign(nsys_all,0.0); S.assign(nsys_all*nsys_all,0.0);

    /* S */
    vector<double> QiQ(numX*numX,0.0),QjQ(numX*numX,0.0),QQQQ(numX*numX,0.0);
    for (int i=0; i<nsys_all; i++) {
        matmul_vec("NN",numX,numX,numX,1.0,QQQ,Qi[i],0.0,QiQ); //QiQ
        for (int j=0; j<nsys_all; j++) {
            matmul_vec("NN",numX,numX,numX,1.0,QQQ,Qi[j],0.0,QjQ); //QjQ
            matmul_vec("NN",numX,numX,numX,1.0,QiQ,QjQ,0.0,QQQQ); //QQQQ
            S[j+i*nsys_all]+=matrace(QQQQ.begin(),numX);
            //this system
            if (j==i) S[j+i*nsys_all]+=sys_n[i]-2*matrace(QiQ,numX);
        }
    }
    if (solve_line("NNT",S,W,sgm,nsys_all,1)) return -1;

    /* normalize and save sgm */
    if (H_niter==0) {
        sgm2.assign(nsys_all,1.0); fsgm.assign(nsys_all,1.0); neg_flag.assign(nsys_all,0);
    }
    for (int i=0; i<nsys_all; i++) {
        fsgm[i]=fabs(sgm[i]/sgm[0]);
        sgm2[i]*=fsgm[i];
    }

    /* update helmertR and newRx */
    for (int sys=0,nnow=0; sys<nsys; sys++) { //helmertR
        for (int i=0; i<sys_n[sys]; i++) for (int j=0; j<sys_n[sys]; j++)
            helmertR[j+nnow + (i+nnow)*numL]*=fsgm[sys];
        nnow+=sys_n[sys];
    }
    for (int i=0; i<numX; i++) for (int j=0; j<numX; j++) newRx[j+i*numX]*=fsgm[nsys]; //newRx

    /* continue helmert algorithm? */
    for (int i=0; i<nsys_all; i++) if (fabs(fsgm[i]-1.0)>SGM_THRES) { H_flag=1; break; }
    else H_flag=0;

    PP.clear(); Px.clear();
    Pi.clear(); VPi.clear(); Ai.clear(); APi.clear();
    S.clear(); W.clear();
    for (int i=0; i<NSYS+1; i++) Qi[i].clear();

    return 0;
}
/* Helmert component covariance estimate for multi-GNSS observation (simple) ---------
* formula : S * sgm2 = W
* base matrix:
* P   = R^-1
* Px  = Qx = Rx^-1
* for the ith system :
* Sii = ni - tr(Qi*QQQ^-1)
* Sij = 0 (j!=i)
* Wi  = Vi*Pi*Vi
* --------------------------------------------------------------------------------- */
int adjfunc_c::helmert_simple(const vector<double> &A,int numL,int numX,int numF,
    const int nobs[NSYS][NFREQ*2]) {
    /* initialization */
    vector<double> PP(helmertR.begin(),helmertR.end()),Px(newRx.begin(),newRx.end()),
        Qi[NSYS+1],Pi,VPi,Ai,APi,S,W;

    /* compute P */
    if (matinv(PP,numL)==-1||matinv(Px,numX)==-1) return -1;

    sys_flag.clear();
    for (int i=0; i<NSYS+1; i++) sys_n[i]=0;
    /* loop of satellite systems to get Qi and W */
    for (int sys=nsys=0,nnow=0; sys<NSYS; sys++) {
        for (int fff=0; fff<numF*2; fff++) sys_n[nsys]+=nobs[sys][fff];
        if (sys_n[nsys]==0) continue;

        /* initialize vectors */
        VPi.assign(sys_n[nsys],0.0); Pi.assign(sys_n[nsys]*sys_n[nsys],0.0);
        APi.assign(numX*sys_n[nsys],0.0);
        Ai.assign(A.begin()+nnow*numX,A.begin()+(nnow+sys_n[nsys])*numX); //Ai
        Qi[nsys].assign(numX*numX,0.0); W.push_back(0.0);

        for (int i=0; i<sys_n[nsys]; i++) for (int j=0; j<sys_n[nsys]; j++)
            Pi[j + i*sys_n[nsys]]=PP[j+nnow + (i+nnow)*numL]; // Pi
        matmul_vec("TN",1,sys_n[nsys],sys_n[nsys],1.0,V,Pi,0.0,VPi,nnow,0,0); //VPi
        matmul_vec("NN",1,1,sys_n[nsys],1.0,VPi,V,0.0,W,0,nnow,nsys); //W[nsys]

        matmul_vec("TN",numX,sys_n[nsys],sys_n[nsys],1.0,Ai,Pi,0.0,APi); //APi
        matmul_vec("NN",numX,numX,sys_n[nsys],1.0,APi,Ai,0.0,Qi[nsys]); //Qi[nsys]

        /* update system information */
        sys_flag.push_back("A"+SYS_LIST.substr(sys,1)); //A+(GREC)
        nnow+=sys_n[nsys]; nsys++;
    }
    if (nsys<0) return -1;

    /* Qi and W of X parameter */
    sys_n[nsys]=numX; sys_flag.push_back("X");
    int nsys_all=nsys+1; vector<double> dXP(numX,0.0);
    Qi[nsys].assign(Px.begin(),Px.end()); //Qx
    W.push_back(0.0);
    matmul_vec("TN",1,numX,numX,1.0,dX,Px,0.0,dXP); //dXP
    matmul_vec("NN",1,1,numX,1.0,dXP,dX,0.0,W,0,0,nsys); //Wx

    /* initilaize sgm2 and S */
    sgm.assign(nsys_all,0.0); S.assign(nsys_all,0.0);

    /* S */
    vector<double> QiQ(numX*numX,0.0),QQQQ(numX*numX,0.0);
    for (int i=0; i<nsys_all; i++) {
        matmul_vec("NN",numX,numX,numX,1.0,QQQ,Qi[i],0.0,QiQ); //QiQ
        S[i]=sys_n[i]-matrace(QiQ,numX);
        sgm[i]=W[i]/S[i];
    }

    /* normalize and save sgm */
    if (H_niter==0) {
        sgm2.assign(nsys_all,1.0); fsgm.assign(nsys_all,1.0);
    }
    for (int i=0; i<nsys; i++) {
        fsgm[i]=sgm[i]/sgm[0];
        sgm2[i]*=fsgm[i];
    }

    /* update helmertR and newRx */
    for (int sys=0,nnow=0; sys<nsys; sys++) { //helmertR
        for (int i=0; i<sys_n[sys]; i++) for (int j=0; j<sys_n[sys]; j++)
            helmertR[j+nnow + (i+nnow)*numL]*=fsgm[sys];
        nnow+=sys_n[sys];
    }
    for (int i=0; i<numX; i++) for (int j=0; j<numX; j++) newRx[j+i*numX]*=fsgm[nsys]; //newRx

    /* continue helmert algorithm? */
    for (int i=0; i<nsys; i++) if (fabs(fsgm[i]-1.0)>SGM_THRES) { H_flag=1; break; }
    else H_flag=0;

    PP.clear(); Px.clear();
    Pi.clear(); VPi.clear(); Ai.clear(); APi.clear();
    S.clear(); W.clear();
    for (int i=0; i<NSYS+1; i++) Qi[i].clear();

    return 0;
}
/* HVCE for multi-GNSS observation (phase and code) ------------------------------- */
int adjfunc_c::helmert_pc(const vector<double> &A,int numL,int numX,int numF,
    const int nobs[NSYS][NFREQ*2]) {
    /* initialization */
    vector<double> PP(helmertR.begin(),helmertR.end()),Px(newRx.begin(),newRx.end()),
        Qi[NSYS*2+1],Pi,VPi,Ai,APi,S,W;

    /* compute P */
    if (matinv(PP,numL)==-1||matinv(Px,numX)==-1) return -1;

    /* loop of satellite systems to get Qi and W */
    int sysobs[NSYS*2+1]={ 0 };
    int C0[NSYS*2+1]={ 0 };
    sys_flag.clear();
    for (int sys=nsys=0,nnow=0; sys<NSYS; sys++) {
        for (int pc=0; pc<2; pc++) {
            for (int fff=0; fff<numF; fff++) sysobs[nsys]+=nobs[sys][pc*numF+fff];
            if (sysobs[nsys]<=0) continue;
            C0[nsys]=pc;

            /* initialize vectors */
            VPi.assign(sysobs[nsys],0.0); Pi.assign(sysobs[nsys]*sysobs[nsys],0.0);
            APi.assign(numX*sysobs[nsys],0.0);
            Ai.assign(A.begin()+nnow*numX,A.begin()+(nnow+sysobs[nsys])*numX); //Ai
            Qi[nsys].assign(numX*numX,0.0); W.push_back(0.0);

            for (int i=0; i<sysobs[nsys]; i++) for (int j=0; j<sysobs[nsys]; j++)
                Pi[j + i*sysobs[nsys]]=PP[j+nnow + (i+nnow)*numL]; // Pi
            matmul_vec("TN",1,sysobs[nsys],sysobs[nsys],1.0,V,Pi,0.0,VPi,nnow,0,0); //VPi
            matmul_vec("NN",1,1,sysobs[nsys],1.0,VPi,V,0.0,W,0,nnow,nsys); //W[nsys]

            matmul_vec("TN",numX,sysobs[nsys],sysobs[nsys],1.0,Ai,Pi,0.0,APi); //APi
            matmul_vec("NN",numX,numX,sysobs[nsys],1.0,APi,Ai,0.0,Qi[nsys]); //Qi[nsys]

            /* update system information */
            sys_flag.push_back(PC_LIST[pc]+SYS_LIST.substr(sys,1)); //(PC)+(GREC)
            nnow+=sysobs[nsys]; nsys++;
        }
    }
    if (nsys<0) return -1;
    /* Qi and W of X parameter */
    sysobs[nsys]=numX; sys_flag.push_back("X"); C0[nsys]=1;
    int nsys_all=nsys+1; vector<double> dXP(numX,0.0);
    Qi[nsys].assign(Px.begin(),Px.end()); //Qx
    W.push_back(0.0);
    matmul_vec("TN",1,numX,numX,1.0,dX,Px,0.0,dXP); //dXP
    matmul_vec("NN",1,1,numX,1.0,dXP,dX,0.0,W,0,0,nsys); //Wx

    /* initilaize sgm2 and S */
    sgm.assign(nsys_all,0.0); S.assign(nsys_all,0.0);

    /* S */
    vector<double> QiQ(numX*numX,0.0),QQQQ(numX*numX,0.0);
    for (int i=0; i<nsys_all; i++) {
        for (int j=0; j<numX; j++) for (int k=0; k<numX; k++)
            QQQQ[k+j*numX]+=Qi[i][k+j*numX];
    }
    matinv(QQQQ,numX);
    for (int i=0; i<nsys_all; i++) {
        matmul_vec("NN",numX,numX,numX,1.0,QQQQ,Qi[i],0.0,QiQ); //QiQ
        S[i]=sysobs[i]-matrace(QiQ,numX);
        sgm[i]=W[i]/S[i];
    }

    /* normalize and save sgm */
    if (H_niter==0) {
        sgm2.assign(nsys_all,1.0); fsgm.assign(nsys_all,1.0);
    }
    /* end helmert if not convergent */
    for (int i=0; i<nsys_all; i++) {
        fsgm[i]=sgm[i]/sgm[0];
        if (i==nsys) break;
        if (H_niter&&fabs(pre_fsgm[i]-1.0)>MAX_DSGM&&fabs(fsgm[i]-1.0)>MAX_DSGM) {
            H_flag=0; H_niter=0;
            fsgm.assign(pre_fsgm.begin(),pre_fsgm.end());
            return -1;
        }
        /* stop if sgm is weird */
        if (fabs(fsgm[i])<=RLOW_THRES || fabs(fsgm[i])>=RHIG_THRES) {
            H_niter=0; H_flag=0;
            sgm2.assign(nsys_all,0.0); fsgm.assign(nsys_all,0.0);
            return -1;
        }
    }
    pre_fsgm.assign(fsgm.begin(),fsgm.end());
    /* update sgm2 */
    for (int i=0; i<nsys_all; i++) sgm2[i]*=fsgm[i];


    /* update helmertR and newRx */
    for (int sys=0,nnow=0; sys<nsys; sys++) { //helmertR
        for (int i=0; i<sysobs[sys]; i++) for (int j=0; j<sysobs[sys]; j++)
            helmertR[j+nnow + (i+nnow)*numL]*=fsgm[sys];
        nnow+=sysobs[sys];
    }
    //for (int i=0; i<numX; i++) for (int j=0; j<numX; j++) newRx[j+i*numX]*=fsgm[nsys]; //newRx

    /* continue helmert algorithm? */
    for (int i=0; i<nsys; i++) {
        /* continue if can be convergent */
        if (fabs(fsgm[i]-1.0)>SGM_THRES) { H_flag=1; break; }
        /* end loop if already converged */
        else H_flag=0;
    }

    PP.clear(); Px.clear();
    Pi.clear(); VPi.clear(); Ai.clear(); APi.clear();
    S.clear(); W.clear(); QQQQ.clear();
    for (int i=0; i<NSYS*2+1; i++) Qi[i].clear();

    return 0;
}
/* Helmert component covariance estimate for multi-GNSS observation (phase) ------- */
int adjfunc_c::helmert_phase(const vector<double> &A,int numL,int numX,int numF,
    const int nobs[NSYS][NFREQ*2]) {
    /* initialization */
    vector<double> PP(helmertR.begin(),helmertR.end()),Px(newRx.begin(),newRx.end()),
        Qi[NSYS+1],Pi,VPi,Ai,APi,S,W;

    /* compute P */
    if (matinv(PP,numL)==-1||matinv(Px,numX)==-1) return -1;

    for (int i=0; i<NSYS+1; i++) sys_n[i]=0;
    /* loop of satellite systems to get Qi and W */
    int sysobs[NSYS]={ 0 };
    sys_flag.clear();
    for (int sys=nsys=0,nnow=0; sys<NSYS; sys++) {
        for (int fff=0; fff<numF*2; fff++) {
            sysobs[nsys]+=nobs[sys][fff];
            if (fff<numF) sys_n[nsys]+=nobs[sys][fff];
        }
        if (sys_n[nsys]<=0) continue;

        /* initialize vectors */
        VPi.assign(sys_n[nsys],0.0); Pi.assign(sys_n[nsys]*sys_n[nsys],0.0);
        APi.assign(numX*sys_n[nsys],0.0);
        Ai.assign(A.begin()+nnow*numX,A.begin()+(nnow+sys_n[nsys])*numX); //Ai
        Qi[nsys].assign(numX*numX,0.0); W.push_back(0.0);

        for (int i=0; i<sys_n[nsys]; i++) for (int j=0; j<sys_n[nsys]; j++)
            Pi[j + i*sys_n[nsys]]=PP[j+nnow + (i+nnow)*numL]; // Pi
        matmul_vec("TN",1,sys_n[nsys],sys_n[nsys],1.0,V,Pi,0.0,VPi,nnow,0,0); //VPi
        matmul_vec("NN",1,1,sys_n[nsys],1.0,VPi,V,0.0,W,0,nnow,nsys); //W[nsys]

        matmul_vec("TN",numX,sys_n[nsys],sys_n[nsys],1.0,Ai,Pi,0.0,APi); //APi
        matmul_vec("NN",numX,numX,sys_n[nsys],1.0,APi,Ai,0.0,Qi[nsys]); //Qi[nsys]

        /* update system information */
        sys_flag.push_back("P"+SYS_LIST.substr(sys,1)); //P+(GREC)
        nnow+=sysobs[nsys]; nsys++;
    }
    if (nsys<0) return -1;
    /* Qi and W of X parameter */
    sys_n[nsys]=numX; sys_flag.push_back("X");
    int nsys_all=nsys+1; vector<double> dXP(numX,0.0);
    Qi[nsys].assign(Px.begin(),Px.end()); //Qx
    W.push_back(0.0);
    matmul_vec("TN",1,numX,numX,1.0,dX,Px,0.0,dXP); //dXP
    matmul_vec("NN",1,1,numX,1.0,dXP,dX,0.0,W,0,0,nsys); //Wx

    /* initilaize sgm2 and S */
    sgm.assign(nsys_all,0.0); S.assign(nsys_all,0.0);

    /* S */
    vector<double> QiQ(numX*numX,0.0),QQQQ(numX*numX,0.0);
    for (int i=0; i<nsys_all; i++) {
        matmul_vec("NN",numX,numX,numX,1.0,QQQ,Qi[i],0.0,QiQ); //QiQ
        S[i]=sys_n[i]-matrace(QiQ,numX);
        sgm[i]=W[i]/S[i];
    }

    /* normalize and save sgm */
    if (H_niter==0) {
        sgm2.assign(nsys_all,1.0); fsgm.assign(nsys_all,1.0);
    }
    for (int i=0; i<nsys_all; i++) {
        fsgm[i]=sgm[i]/sgm[0];
        sgm2[i]*=fsgm[i];
    }

    /* update helmertR and newRx */
    for (int sys=0,nnow=0; sys<nsys; sys++) { //helmertR
        for (int i=0; i<sys_n[sys]; i++) for (int j=0; j<sys_n[sys]; j++)
            helmertR[j+nnow + (i+nnow)*numL]*=fsgm[sys];
        nnow+=sysobs[sys];
    }
    for (int i=0; i<numX; i++) for (int j=0; j<numX; j++) newRx[j+i*numX]*=fsgm[nsys]; //newRx

    /* continue helmert algorithm? */
    for (int i=0; i<nsys_all; i++) if (fabs(fsgm[i]-1.0)>SGM_THRES) { H_flag=1; break; }
    else H_flag=0;

    PP.clear(); Px.clear();
    Pi.clear(); VPi.clear(); Ai.clear(); APi.clear();
    S.clear(); W.clear();
    for (int i=0; i<NSYS+1; i++) Qi[i].clear();

    return 0;
}
/* Helmert component covariance estimate for multi-GNSS observation (code) -------- */
int adjfunc_c::helmert_code(const vector<double> &A,int numL,int numX,int numF,
    const int nobs[NSYS][NFREQ*2]) {
    /* initialization */
    vector<double> PP(helmertR.begin(),helmertR.end()),Px(newRx.begin(),newRx.end()),
        Qi[NSYS+1],Pi,VPi,Ai,APi,S,W;

    /* compute P */
    if (matinv(PP,numL)==-1||matinv(Px,numX)==-1) return -1;

    for (int i=0; i<NSYS+1; i++) sys_n[i]=0;
    /* loop of satellite systems to get Qi and W */
    int sysobs[NSYS]={ 0 },nphase[NSYS]={ 0 };
    sys_flag.clear();
    for (int sys=nsys=0,nnow=0; sys<NSYS; sys++) {
        for (int fff=0; fff<numF*2; fff++) {
            sysobs[nsys]+=nobs[sys][fff];
            if (fff>=numF) sys_n[nsys] +=nobs[sys][fff];
            if (fff< numF) nphase[nsys]+=nobs[sys][fff];
        }
        if (sys_n[nsys]<=0) continue;

        /* initialize vectors */
        VPi.assign(sys_n[nsys],0.0); Pi.assign(sys_n[nsys]*sys_n[nsys],0.0);
        APi.assign(numX*sys_n[nsys],0.0);
        Ai.assign(A.begin()+(nnow+nphase[nsys])*numX,A.begin()+(nnow+nphase[nsys]+sys_n[nsys])*numX); //Ai
        Qi[nsys].assign(numX*numX,0.0); W.push_back(0.0);

        for (int i=0; i<sys_n[nsys]; i++) for (int j=0; j<sys_n[nsys]; j++)
            Pi[j + i*sys_n[nsys]]=PP[j+nnow+nphase[nsys] + (i+nnow+nphase[nsys])*numL]; // Pi
        matmul_vec("TN",1,sys_n[nsys],sys_n[nsys],1.0,V,Pi,0.0,VPi,nnow+nphase[nsys],0,0); //VPi
        matmul_vec("NN",1,1,sys_n[nsys],1.0,VPi,V,0.0,W,0,nnow+nphase[nsys],nsys); //W[nsys]

        matmul_vec("TN",numX,sys_n[nsys],sys_n[nsys],1.0,Ai,Pi,0.0,APi); //APi
        matmul_vec("NN",numX,numX,sys_n[nsys],1.0,APi,Ai,0.0,Qi[nsys]); //Qi[nsys]

        /* update system information */
        sys_flag.push_back("C"+SYS_LIST.substr(sys,1)); //C+(GREC)
        nnow+=sysobs[nsys]; nsys++;
    }
    if (nsys<0) return -1;
    /* Qi and W of X parameter */
    sys_n[nsys]=numX; sys_flag.push_back("X");
    int nsys_all=nsys+1; vector<double> dXP(numX,0.0);
    Qi[nsys].assign(Px.begin(),Px.end()); //Qx
    W.push_back(0.0);
    matmul_vec("TN",1,numX,numX,1.0,dX,Px,0.0,dXP); //dXP
    matmul_vec("NN",1,1,numX,1.0,dXP,dX,0.0,W,0,0,nsys); //Wx

    /* initilaize sgm2 and S */
    sgm.assign(nsys_all,0.0); S.assign(nsys_all,0.0);

    /* S */
    vector<double> QiQ(numX*numX,0.0),QQQQ(numX*numX,0.0);
    for (int i=0; i<nsys_all; i++) {
        matmul_vec("NN",numX,numX,numX,1.0,QQQ,Qi[i],0.0,QiQ); //QiQ
        S[i]=sys_n[i]-matrace(QiQ,numX);
        sgm[i]=W[i]/S[i];
    }

    /* normalize and save sgm */
    if (H_niter==0) {
        sgm2.assign(nsys_all,1.0); fsgm.assign(nsys_all,1.0);
    }
    for (int i=0; i<nsys; i++) {
        fsgm[i]=sgm[i]/sgm[0];
        sgm2[i]*=fsgm[i];
    }

    /* update helmertR and newRx */
    for (int sys=0,nnow=0; sys<nsys; sys++) { //helmertR
        for (int i=0; i<sys_n[sys]; i++) for (int j=0; j<sys_n[sys]; j++)
            helmertR[j+nnow+nphase[sys] + (i+nnow+nphase[sys])*numL]*=fsgm[sys];
        nnow+=sysobs[sys];
    }
    for (int i=0; i<numX; i++) for (int j=0; j<numX; j++) newRx[j+i*numX]*=fsgm[nsys]; //newRx

    /* continue helmert algorithm? */
    for (int i=0; i<nsys; i++) if (fabs(fsgm[i]-1.0)>SGM_THRES) { H_flag=1; break; }
    else H_flag=0;

    PP.clear(); Px.clear();
    Pi.clear(); VPi.clear(); Ai.clear(); APi.clear();
    S.clear(); W.clear();
    for (int i=0; i<NSYS+1; i++) Qi[i].clear();

    return 0;
}
/* robust functions --------------------------------------------------------------- */
/* robust empty function ---------------------------------------------------------- */
int adjfunc_c::robust_empty(int numL,int numX,vector<double> &R) {
    return largeResFlag=0;
}
/* Yang IGG3 robust algorithm with 2 intervals --------------------------------------
*      sV(i) = V(i)/sqrt(R(i)) (standard residual)
* factor(ii) = 1            , if sV(i) <= k0
* factor(ii) = sV(i)/k0     , if sV(i) >  k0
* factor(ij) = sqrt(factor(ii)*factor(jj))
* -------------------------------------------------------------------------------- */
int adjfunc_c::robust_yang_IGG3_2s(int numL,int numX,vector<double> &R) {
    int i;
    int halfL=numL/2;
    normV.assign(numL,0);
    robustFactor.assign(numL,0); largeResFlag=0;

    /* get sigma0 */
    sig0=0;
    for (i=0; i<numL; i++) normV[i]=fabs(V[i]/SQRT(R[i+i*numL]));
    sortV.assign(normV.begin(),normV.end());
    sort(sortV.begin(),sortV.end());
    sig0=ROBUST_RES_FACTOR*sortV[halfL];

    /* loop to get factors and update R */
    for (i=0; i<numL; i++) {
        robustFactor[i]=normV[i]/sig0/robust2sK0;
        /* large residual detected */
        if (robustFactor[i]>ROBUST_CUTOFF) {
            robustFactor[i]=ROBUST_CUTOFF;
            largeResFlag++;
            /* update covariance */
            R[i+i*numL]*=robustFactor[i];
            for (int j=0; j<i; j++) R[j+i*numL]=(R[i+j*numL]*=ROBUST_CUTOFF_SQRT);
            for (int j=i+1; j<numL; j++) R[j+i*numL]=(R[i+j*numL]*=ROBUST_CUTOFF_SQRT);
        }
        else if (robustFactor[i]>1.0) {
            buffFactor=sqrt(robustFactor[i]);
            largeResFlag++;
            /* update covariance */
            R[i+i*numL]*=robustFactor[i];
            for (int j=0; j<i; j++) R[j+i*numL]=(R[i+j*numL]*=buffFactor);
            for (int j=i+1; j<numL; j++) R[j+i*numL]=(R[i+j*numL]*=buffFactor);
        }
        else robustFactor[i]=1.0;
    }
    return largeResFlag;
}
/* Yang IGG3 robust algorithm with 3 intervals --------------------------------------
*      sV(i) = V(i)/sqrt(R(i)) (standard residual)
* factor(ii) = 1                               , if sV(i) <= k0
* factor(ii) = sV(i)/k0*((k1-k0)/(k1-sV(i))^2  , if sV(i) >  k0  and  sV(i) <  k1
*              (cut off in 10000)
* factor(ii) = 10000                           , if sV(i) >= k1
* factor(ij) = sqrt(factor(ii)*factor(jj))
* -------------------------------------------------------------------------------- */
int adjfunc_c::robust_yang_IGG3_3s(int numL,int numX,vector<double> &R) {
    int i;
    double factor;
    int halfL=numL/2;
    normV.assign(numL,0);
    robustFactor.assign(numL,0); largeResFlag=0;

    /* get sigma0 */
    sig0=0;
    for (i=0; i<numL; i++) normV[i]=fabs(V[i]/SQRT(R[i+i*numL]));
    sortV.assign(normV.begin(),normV.end());
    sort(sortV.begin(),sortV.end());
    sig0=ROBUST_RES_FACTOR*sortV[halfL];

    /* loop to get factors and update R */
    for (i=0; i<numL; i++) {
        buffFactor=normV[i]/sig0;
        /* large residual detected */
        if (buffFactor>robust3sK0) {
            if (buffFactor>=robust3sK1) { //cutoff at k1
                robustFactor[i]=ROBUST_CUTOFF;
                largeResFlag++;
                /* update covariance */
                R[i+i*numL]*=robustFactor[i];
                for (int j=0; j<i; j++) R[j+i*numL]=(R[i+j*numL]*=ROBUST_CUTOFF_SQRT);
                for (int j=i+1; j<numL; j++) R[j+i*numL]=(R[i+j*numL]*=ROBUST_CUTOFF_SQRT);
            }
            else {
                robust3sFactor=robust3sK0*SQR((robust3sK1-buffFactor)/(robust3sK1-robust3sK0))/buffFactor;
                /* sV>k1 or factor>cutoff */
                if ( robust3sFactor<=ROBUST_CUTOFF_INVERSE ) {
                    robustFactor[i]=ROBUST_CUTOFF;
                    largeResFlag++;
                    /* update covariance */
                    R[i+i*numL]*=robustFactor[i];
                    for (int j=0; j<i; j++) R[j+i*numL]=(R[i+j*numL]*=ROBUST_CUTOFF_SQRT);
                    for (int j=i+1; j<numL; j++) R[j+i*numL]=(R[i+j*numL]*=ROBUST_CUTOFF_SQRT);
                }
                else {
                    robustFactor[i]=1.0/SQRT(robust3sFactor);
                    largeResFlag++;
                    /* update covariance */
                    for (int j=0; j<i; j++) R[j+i*numL]=(R[i+j*numL]*=robustFactor[i]);
                    for (int j=i+1; j<numL; j++) R[j+i*numL]=(R[i+j*numL]*=robustFactor[i]);
                    robustFactor[i]=1.0/robust3sFactor;
                    R[i+i*numL]*=robustFactor[i];
                }
            }
        }
        else robustFactor[i]=1.0;
    }
    return largeResFlag;
}
/* Yang IGG3 robust algorithm with 2 intervals only for R diagonal ---------------- */
int adjfunc_c::robust_yang_IGG3_2s_Rdiag(int numL,int numX,vector<double> &R) {
    int i;
    int halfL=numL/2;
    normV.assign(numL,0);
    robustFactor.assign(numL,0); largeResFlag=0;

    /* get sigma0 */
    sig0=0;
    for (i=0; i<numL; i++) normV[i]=fabs(V[i]/SQRT(R[i*numL+i]));
    sortV.assign(normV.begin(),normV.end());
    sort(sortV.begin(),sortV.end());
    sig0=ROBUST_RES_FACTOR*sortV[halfL];

    /* loop to get factors and update Rvec */
    for (i=0; i<numL; i++) {
        robustFactor[i]=normV[i]/sig0/robust2sK0;
        /* large residual detected */
        if (robustFactor[i]>ROBUST_CUTOFF) {
            robustFactor[i]=ROBUST_CUTOFF;
            largeResFlag++;
            /* update covariance */
            R[i*numL+i]*=ROBUST_CUTOFF;
        }
        else if (robustFactor[i]>1.0) {
            largeResFlag++;
            /* update covariance */
            R[i*numL+i]*=robustFactor[i];
        }
        else robustFactor[i]=1.0;
    }
    return largeResFlag;
}
/* Yang IGG3 robust algorithm with 3 intervals only for R diagonal ---------------- */
int adjfunc_c::robust_yang_IGG3_3s_Rdiag(int numL,int numX,vector<double> &R) {
    int i;
    double factor;
    int halfL=numL/2;
    normV.assign(numL,0);
    robustFactor.assign(numL,0); largeResFlag=0;

    /* get sigma0 */
    sig0=0;
    for (i=0; i<numL; i++) normV[i]=fabs(V[i]/SQRT(R[i*numL+i]));
    sortV.assign(normV.begin(),normV.end());
    sort(sortV.begin(),sortV.end());
    sig0=ROBUST_RES_FACTOR*sortV[halfL];

    /* loop to get factors and update R */
    for (i=0; i<numL; i++) {
        buffFactor=normV[i]/sig0;
        /* large residual detected */
        if (buffFactor>robust3sK0) {
            if (buffFactor>=robust3sK1) { //cutoff at k1
                robustFactor[i]=ROBUST_CUTOFF;
                largeResFlag++;
                R[i*numL+i]*=ROBUST_CUTOFF;
            }
            else {
                robust3sFactor=robust3sK0*SQR((robust3sK1-buffFactor)/(robust3sK1-robust3sK0))/buffFactor;
                /* sV>k1 or factor>cutoff */
                if ( robust3sFactor<=ROBUST_CUTOFF_INVERSE ) {
                    robustFactor[i]=ROBUST_CUTOFF;
                    largeResFlag++;
                    R[i*numL+i]*=ROBUST_CUTOFF;
                }
                else {
                    robustFactor[i]=1.0/robust3sFactor;
                    largeResFlag++;
                    /* update covariance */
                    R[i*numL+i]*=robustFactor[i];
                }
            }
        }
        else robustFactor[i]=1.0;
    }
    return largeResFlag;
}
/* Yang IGG3 robust algorithm with 2 intervals for variance vector Rvec ----------- */
int adjfunc_c::robust_yang_IGG3_2s_Rvec(int numL,int numX,vector<double> &Rvec) {
    int i;
    int halfL=numL/2;
    normV.assign(numL,0);
    robustFactor.assign(numL,0); largeResFlag=0;

    /* get sigma0 */
    sig0=0;
    for (i=0; i<numL; i++) normV[i]=fabs(V[i]/SQRT(Rvec[i]));
    sortV.assign(normV.begin(),normV.end());
    sort(sortV.begin(),sortV.end());
    sig0=ROBUST_RES_FACTOR*sortV[halfL];

    /* loop to get factors and update Rvec */
    for (i=0; i<numL; i++) {
        robustFactor[i]=normV[i]/sig0/robust2sK0;
        /* large residual detected */
        if (robustFactor[i]>ROBUST_CUTOFF) {
            robustFactor[i]=ROBUST_CUTOFF;
            largeResFlag++;
            /* update covariance */
            Rvec[i]*=ROBUST_CUTOFF;
        }
        else if (robustFactor[i]>1.0) {
            largeResFlag++;
            /* update covariance */
            Rvec[i]*=robustFactor[i];
        }
        else robustFactor[i]=1.0;
    }
    return largeResFlag;
}
/* Yang IGG3 robust algorithm with 3 intervals for variance vector Rvec ----------- */
int adjfunc_c::robust_yang_IGG3_3s_Rvec(int numL,int numX,vector<double> &Rvec) {
    int i;
    double factor;
    int halfL=numL/2;
    normV.assign(numL,0);
    robustFactor.assign(numL,0); largeResFlag=0;

    /* get sigma0 */
    sig0=0;
    for (i=0; i<numL; i++) normV[i]=fabs(V[i]/SQRT(Rvec[i]));
    sortV.assign(normV.begin(),normV.end());
    sort(sortV.begin(),sortV.end());
    sig0=ROBUST_RES_FACTOR*sortV[halfL];

    /* loop to get factors and update R */
    for (i=0; i<numL; i++) {
        buffFactor=normV[i]/sig0;
        /* large residual detected */
        if (buffFactor>robust3sK0) {
            if (buffFactor>=robust3sK1) { //cutoff at k1
                robustFactor[i]=ROBUST_CUTOFF;
                largeResFlag++;
                Rvec[i]*=ROBUST_CUTOFF;
            }
            else {
                robust3sFactor=robust3sK0*SQR((robust3sK1-buffFactor)/(robust3sK1-robust3sK0))/buffFactor;
                /* sV>k1 or factor>cutoff */
                if ( robust3sFactor<=ROBUST_CUTOFF_INVERSE ) {
                    robustFactor[i]=ROBUST_CUTOFF;
                    largeResFlag++;
                    Rvec[i]*=ROBUST_CUTOFF;
                }
                else {
                    robustFactor[i]=1.0/robust3sFactor;
                    largeResFlag++;
                    /* update covariance */
                    Rvec[i]*=robustFactor[i];
                }
            }
        }
        else robustFactor[i]=1.0;
    }
    return largeResFlag;
}
/* Yang IGG3 robust algorithm with 2 intervals for each constellation (GNSS) ---------
 *      sV(i) = V(i)/sqrt(R(i)) (standard residual)
 * factor(ii) = 1            , if sV(i) <= k0
 * factor(ii) = sV(i)/k0     , if sV(i) >  k0
 * factor(ij) = sqrt(factor(ii)*factor(jj))
 * -------------------------------------------------------------------------------- */
int adjfunc_c::robust_yang_IGG3_2s_GNSS(int numL,int numX,int numF,vector<double> &R,
    const int nobs[NSYS][NFREQ*2]) {

    int halfL=numL/2;
    normV.assign(numL,0);
    robustFactor.assign(numL,0); largeResFlag=0;

    /* get sigma0 */
    sig0=0;
    for (int i=0; i<numL; i++) normV[i] = fabs(V[i]/SQRT(R[i+i*numL]));
    sortV.assign(normV.begin(),normV.end());
    sort(sortV.begin(),sortV.end());
    sig0=ROBUST_RES_FACTOR*sortV[halfL];

    for (int sys=0,nnow=0; sys<NSYS; sys++) {
        for (int ipc=0; ipc<2; ipc++) {
            int numobs=0;
            for (int fff=0; fff<numF; fff++) numobs+=nobs[sys][ipc*numF+fff];
            if (numobs<=0) continue;
            for (int i=0; i<numobs; i++) {
                robustFactor[i+nnow]=normV[i+nnow]/sig0/robust2sK0;
                /* large residual detected */
                if (robustFactor[i+nnow]>ROBUST_CUTOFF) {
                    robustFactor[i+nnow]=ROBUST_CUTOFF;
                    largeResFlag++;
                    /* update covariance */
                    for (int j=0; j<numobs; j++) {
                        R[j+nnow+(i+nnow)*numL]*=ROBUST_CUTOFF_SQRT;
                        R[i+nnow+(j+nnow)*numL]*=ROBUST_CUTOFF_SQRT;
                    }
                }
                else if (robustFactor[i+nnow]>1.0) {
                    buffFactor=sqrt(robustFactor[i+nnow]);
                    largeResFlag++;
                    /* update covariance */
                    for (int j=0; j<numobs; j++) {
                        R[j+nnow+(i+nnow)*numL]*=buffFactor;
                        R[i+nnow+(j+nnow)*numL]*=buffFactor;
                    }
                }
                else robustFactor[i+nnow]=1.0;
            }
            nnow+=numobs;
        }
    }
    return largeResFlag;
}
/* Yang IGG3 robust algorithm with 3 intervals for each constellation (GNSS) ---------
*      sV(i) = V(i)/sqrt(R(i)) (standard residual)
* factor(ii) = 1                               , if sV(i) <= k0
* factor(ii) = sV(i)/k0*((k1-k0)/(k1-sV(i))^2  , if sV(i) >  k0  and  sV(i) <  k1
*              (cut off in 10000)
* factor(ii) = 10000                           , if sV(i) >= k1
* factor(ij) = sqrt(factor(ii)*factor(jj))
* -------------------------------------------------------------------------------- */
int adjfunc_c::robust_yang_IGG3_3s_GNSS(int numL,int numX,int numF,vector<double> &R,
    const int nobs[NSYS][NFREQ*2]) {

    int halfL=numL/2;
    normV.assign(numL,0);
    robustFactor.assign(numL,0); largeResFlag=0;

    /* get sigma0 */
    sig0=0;
    for (int i=0; i<numL; i++) normV[i] = fabs(V[i]/SQRT(R[i+i*numL]));
    sortV.assign(normV.begin(),normV.end());
    sort(sortV.begin(),sortV.end());
    sig0=ROBUST_RES_FACTOR*sortV[halfL];

    /* loop to get factors */
    for (int sys=0,nnow=0; sys<NSYS; sys++) {
        for (int ipc=0; ipc<2; ipc++) {
            int numobs=0;
            for (int fff=0; fff<numF; fff++) numobs+=nobs[sys][ipc*numF+fff];
            if (numobs<=0) continue;
            for (int i=0; i<numobs; i++) {
                buffFactor=normV[i+nnow]/sig0;
                /* large residual detected */
                if (buffFactor>robust3sK0) {
                    if (buffFactor>=robust3sK1) { //cutoff at k1
                        robustFactor[i+nnow]=ROBUST_CUTOFF;
                        largeResFlag++;
                        for (int j=0; j<numobs; j++) {
                            R[j+nnow+(i+nnow)*numL]*=ROBUST_CUTOFF_SQRT;
                            R[i+nnow+(j+nnow)*numL]*=ROBUST_CUTOFF_SQRT;
                        }
                    }
                    else {
                        robust3sFactor=robust3sK0*SQR((robust3sK1-buffFactor)/(robust3sK1-robust3sK0))/buffFactor;
                        /* sV>k1 or factor>cutoff */
                        if ( robust3sFactor<=ROBUST_CUTOFF_INVERSE ) {
                            robustFactor[i+nnow]=ROBUST_CUTOFF;
                            largeResFlag++;
                            for (int j=0; j<numobs; j++) {
                                R[j+nnow+(i+nnow)*numL]*=ROBUST_CUTOFF_SQRT;
                                R[i+nnow+(j+nnow)*numL]*=ROBUST_CUTOFF_SQRT;
                            }
                        }
                        else {
                            robustFactor[i+nnow]=1.0/SQRT(robust3sFactor);
                            largeResFlag++;
                            /* update covariance */
                            for (int j=0; j<numobs; j++) {
                                R[j+nnow+(i+nnow)*numL]*=robustFactor[i+nnow];
                                R[i+nnow+(j+nnow)*numL]*=robustFactor[i+nnow];
                            }
                            robustFactor[i+nnow]=1.0/robust3sFactor;
                        }
                    }
                }
                else robustFactor[i+nnow]=1.0;
            }
            nnow+=numobs;
        }
    }
    return largeResFlag;
}
/* least-square adjustment function --------------------------------------------------
/* estimate parameter X array and covariance Rx:
*   P=R^-1
*   QQQ=(A'*P*A)^-1, APL=A'*P*L
*   xp=x+QQQ*APL, Rx(new)=QQQ
* --------------------------------------------------------------------------------- */
int adjfunc_c::LSQ(const vector<double> &A,const vector<double> &L,
    vector<double> &R,vector<double> &X,vector<double> &Rx,int numL,int numX) {

    sigma0_2=1E9;
    if (numL<numX) return flag=-1;
    /* initialize matrixes */
    AP.assign(A.size(),0); //numX x numL
    APL.assign(numX,0); dX.assign(numX,0.0); 
    QQQ.assign(Rx.size(),0.0); //numX x numX
    P.assign(R.begin(),R.end()); //weight matrix P

    /* least square process */
    if (LSQ_est_dx(A,L,R,numL,numX)>=0) { //estimate dx
        LSQ_update_xRx(L,X,Rx,numL,numX); //update x and Rx
    }

    P.clear(); AP.clear(); APL.clear(); LP.clear();

    return flag;
}
/* least-square adjustment function for variance vector Rvec -------------------------
/* estimate parameter X array and covariance Rx:
*   !!!Rvec consists of the diagonal elements of obs covariance matrix R!!!
*   A1=A/R,L1=L/R
*   QQQ=(A1'*A1)^-1, APL=A1'*L1
*   xp=x+QQQ*APL, Rx(new)=QQQ
* --------------------------------------------------------------------------------- */
int adjfunc_c::LSQ_Rvec(const vector<double> &A,const vector<double> &L,
    vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX) {

    sigma0_2=1E9;
    if (numL<numX) return flag=-1;
    /* initialize matrixes */
    AP.assign(A.size(),0); //numX x numL
    APL.assign(numX,0); dX.assign(numX,0.0);
    QQQ.assign(Rx.size(),0.0); //numX x numX

    /* least square process */
    if (LSQ_est_dx_Rvec(A,L,Rvec,numL,numX)>=0) { //estimate dx
        LSQ_update_xRx_Rvec(L,Rvec,X,Rx,numL,numX); //update x and Rx
    }

    AP.clear(); APL.clear();

    return flag;
}
/* kalman filter ---------------------------------------------------------------------
/* estimate parameter X array and covariance Rx
*   kalman filter state update as follows:
*   Q=A*Rx*A'+R, K=Rx*A'*Q^-1
*   xp=x+K*L, Rx(new)=(E-K*A)*Rx
* --------------------------------------------------------------------------------- */
int adjfunc_c::KF(const vector<double> &A,const vector<double> &L,
    vector<double> &R,vector<double> &X,vector<double> &Rx,int numL,int numX) {

    RxA.assign(A.size(),0.0); //numX x numL
    Q.assign(R.size(),0.0); //numL x numL
    E.assign(Rx.size(),0.0); eyemat_zero(E.begin(),numX);//numX x numX
    ERx.assign(Rx.size(),0.0); //numX x numX
    KR.assign(A.size(),0.0); //numX x numL
    dX.assign(numX,0.0); K.assign(A.size(),0.0); //numX x numL
    QQQ.assign(Rx.size(),0.0); //numX x numX

    /* kalman filter process */
    if (KF_est_dx(A,L,R,Rx,numL,numX)>=0) { //estimate dx
        KF_update_xRx(A,L,R,X,Rx,numL,numX); //update x and Rx
    }

    RxA.clear(); Q.clear(); E.clear(); ERx.clear(); KR.clear();

    return flag;
}
/* kalman filter for variance vector Rvec ----------------------------------------- */
int adjfunc_c::KF_Rvec(const vector<double> &A,const vector<double> &L,
    vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX){

    RxA.assign(A.size(),0.0); //numX x numL
    Q.assign(numL*numL,0.0); //numL x numL
    E.assign(Rx.size(),0.0); eyemat_zero(E.begin(),numX);//numX x numX
    ERx.assign(Rx.size(),0.0); //numX x numX
    KR.assign(A.size(),0.0); //numX x numL
    dX.assign(numX,0.0); K.assign(A.size(),0.0); //numX x numL
    QQQ.assign(Rx.size(),0.0); //numX x numX

    /* kalman filter process */
    if (KF_est_dx_Rvec(A,L,Rvec,Rx,numL,numX)>=0) { //estimate dx
        KF_update_xRx_Rvec(A,L,Rvec,X,Rx,numL,numX); //update x and Rx
    }

    RxA.clear(); Q.clear(); E.clear(); ERx.clear(); KR.clear();

    return flag;
}
/* robust least-square adjustment function ---------------------------------------- */
int adjfunc_c::robust_LSQ(const vector<double> &A,const vector<double> &L,
    vector<double> &R,vector<double> &X,vector<double> &Rx,int numL,int numX) {

    sigma0_2=1E9;
    if (numL<numX) return flag=-1;
    /* initialize matrixes */
    AP.assign(A.size(),0); //numX x numL
    APL.assign(numX,0); dX.assign(numX,0.0); 
    QQQ.assign(Rx.size(),0.0); //numX x numX
    P.assign(R.begin(),R.end()); //weight matrix P

    /* least square and robust process */
    if (LSQ_est_dx(A,L,R,numL,numX)>=0) {
        (this->*robustFcp[robustRtype][robustMode])(numL,numX,R);
        //update dx, x and Rx
        if (LSQ_est_dx(A,L,R,numL,numX)>=0) LSQ_update_xRx(L,X,Rx,numL,numX);
    }

    P.clear(); AP.clear(); APL.clear(); LP.clear();

    return flag;
}
/* robust least-square adjustment function for variance vector Rvec --------------- */
int adjfunc_c::robust_LSQ_Rvec(const vector<double> &A,const vector<double> &L,
    vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX) {

    sigma0_2=1E9;
    if (numL<numX) return flag=-1;
    /* initialize matrixes */
    AP.assign(A.size(),0); //numX x numL
    APL.assign(numX,0); dX.assign(numX,0.0);
    QQQ.assign(Rx.size(),0.0); //numX x numX

    /* least square and robust process */
    if (LSQ_est_dx_Rvec(A,L,Rvec,numL,numX)>=0) {
        (this->*robustFcp[2][robustMode])(numL,numX,Rvec);
        //update dx, x and Rx
        if (LSQ_est_dx_Rvec(A,L,Rvec,numL,numX)>=0) LSQ_update_xRx_Rvec(L,Rvec,X,Rx,numL,numX);
    }

    AP.clear(); APL.clear();

    return flag;
}
/* robust kalman filter ----------------------------------------------------------- */
int adjfunc_c::robust_KF(const vector<double> &A,const vector<double> &L,
    vector<double> &R,vector<double> &X,vector<double> &Rx,int numL,int numX) {

    RxA.assign(A.size(),0.0); //numX x numL
    Q.assign(R.size(),0.0); //numL x numL
    E.assign(Rx.size(),0.0); eyemat_zero(E.begin(),numX);//numX x numX
    ERx.assign(Rx.size(),0.0); //numX x numX
    KR.assign(A.size(),0.0); //numX x numL
    dX.assign(numX,0.0); K.assign(A.size(),0.0); //numX x numL
    QQQ.assign(Rx.size(),0.0); //numX x numX

    /* kalman filter and robust process */
    if (KF_est_dx(A,L,R,Rx,numL,numX)>=0) {
        (this->*robustFcp[robustRtype][robustMode])(numL,numX,R);
        //update dx, x and Rx
        if (KF_est_dx(A,L,R,Rx,numL,numX)>=0) KF_update_xRx(A,L,R,X,Rx,numL,numX);
    }

    RxA.clear(); Q.clear(); E.clear();
    ERx.clear(); KR.clear();

    return flag;
}
/* robust kalman filter for variance vector Rvec ---------------------------------- */
int adjfunc_c::robust_KF_Rvec(const vector<double> &A,const vector<double> &L,
    vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX){

    RxA.assign(A.size(),0.0); //numX x numL
    Q.assign(numL*numL,0.0); //numL x numL
    E.assign(Rx.size(),0.0); eyemat_zero(E.begin(),numX);//numX x numX
    ERx.assign(Rx.size(),0.0); //numX x numX
    KR.assign(A.size(),0.0); //numX x numL
    dX.assign(numX,0.0); K.assign(A.size(),0.0); //numX x numL
    QQQ.assign(Rx.size(),0.0); //numX x numX

    /* kalman filter and robust process */
    if (KF_est_dx_Rvec(A,L,Rvec,Rx,numL,numX)>=0) {
        (this->*robustFcp[2][robustMode])(numL,numX,Rvec);
        //update dx, x and Rx
        if (KF_est_dx_Rvec(A,L,Rvec,Rx,numL,numX)>=0) KF_update_xRx_Rvec(A,L,Rvec,X,Rx,numL,numX);
    }

    RxA.clear(); Q.clear(); E.clear(); ERx.clear(); KR.clear();

    return flag;
}
/* estimate parameter X array and covariance Rx (virtual)-------------------------- */
int adjfunc_c::adjustment(const vector<double> &A,const vector<double> &L,
    vector<double> &R,vector<double> &X,vector<double> &Rx,
    int numL,int numX,int numF,const int nobs[NSYS][NFREQ*2]) {
    return -1;
}

/* subclass adjustment functions ------------------------------------------------------------------ */
/* least square adjustment function --------------------------------------------------------------- */
/* Constructor -------------------------------------------------------------------- */
lsadj_c::lsadj_c() {
    adjustFcp=&adjfunc_c::LSQ;
}
lsadj_c::lsadj_c(int Robust,int RobustCovar) : adjfunc_c(Robust,RobustCovar) {
}
lsadj_c::~lsadj_c() {
}
/* Implementaion functions -------------------------------------------------------- */
/* estimate parameter X array and covariance Rx --------------------------------------
*   nobs[NSYS][NFREQ*2] not used
* --------------------------------------------------------------------------------- */
int lsadj_c::adjustment(const vector<double> &A,const vector<double> &L,
    vector<double> &R,vector<double> &X,vector<double> &Rx,
    const int numL,const int numX,int numF,const int nobs[NSYS][NFREQ*2]) {

    return (this->*adjustFcp)(A,L,R,X,Rx,numL,numX);
}

/* kalman filter adjustment function -------------------------------------------------------------- */
/* Constructor -------------------------------------------------------------------- */
kalmanfilter_c::kalmanfilter_c() {
}
kalmanfilter_c::kalmanfilter_c(int Robust,int RobustCovar) : adjfunc_c(Robust,RobustCovar) {
    if (robustMode>ADJUST_ROBUST_NONE) {
        if (robustRtype==ROBUST_R_VECTOR) adjustFcp=&adjfunc_c::robust_KF_Rvec;
        else adjustFcp=&adjfunc_c::robust_KF;
    }
    else {
        if (robustRtype==ROBUST_R_VECTOR) adjustFcp=&adjfunc_c::KF_Rvec;
        else adjustFcp=&adjfunc_c::KF;
    }
}
kalmanfilter_c::~kalmanfilter_c() {
}
int kalmanfilter_c::adjustment(const vector<double> &A,const vector<double> &L,
    vector<double> &R,vector<double> &X,vector<double> &Rx,
    int numL,int numX,int numF,const int nobs[NSYS][NFREQ*2]) {

    return (this->*adjustFcp)(A,L,R,X,Rx,numL,numX);
}

/* helmert components covariance estimate for kalman filter adjustment function ------------------- */
/* estimate parameter X array and covariance Rx ----------------------------------- */
helmert_c::helmert_c() {
}
helmert_c::helmert_c(int Robust,int RobustCovar) {
}
helmert_c::~helmert_c() {
    helmertR.clear(); newRx.clear();
    buffRx.clear();
}
/* Implementaion functions -------------------------------------------------------- */
/* for test reset covariance matrix of observations ------------------------------- */
void helmert_c::reset_var(vector<double> &R,int numL,int numF,const int nobs[NSYS][NFREQ*2]) {
    int nnow=0;
    for (int i=0; i<NSYS; i++) for (int j=0; j<numF*2; j++) {
        if (nobs[i][j]<=0) continue;
        if (j>=numF) for (int k1=0; k1<nobs[i][j]; k1++) for (int k2=0; k2<nobs[i][j]; k2++)
            R[nnow+k2+(nnow+k1)*numL]*=3.0;
        nnow+=nobs[i][j];
    }
}
int helmert_c::adjustment(const vector<double> &A,const vector<double> &L,
    vector<double> &R,vector<double> &X,vector<double> &Rx,
    int numL,int numX,int numF,const int nobs[NSYS][NFREQ*2]) {

    helmertR.assign(R.begin(),R.end()); newRx.assign(Rx.begin(),Rx.end());
    vector<double> X_ori(X);
    H_flag=1;

    for (H_niter=0; H_flag&&H_niter<MAX_ITER; H_niter++) {
        X_ori.assign(X.begin(),X.end()); buffRx.assign(newRx.begin(),newRx.end());
        if (kalmanfilter_c::robust_KF(A,L,helmertR,X_ori,buffRx,numL,numX)) {
            //reset_var(R,numL,numF,nobs);
            return kalmanfilter_c::robust_KF(A,L,R,X,Rx,numL,numX);
        }
        if (helmert_pc(A,numL,numX,numF,nobs)) {
            //reset_var(R,numL,numF,nobs);
            return kalmanfilter_c::robust_KF(A,L,R,X,Rx,numL,numX);
        }
    }

    /* finial solution */
    Rx.assign(newRx.begin(),newRx.end());
    return kalmanfilter_c::robust_KF(A,L,helmertR,X,Rx,numL,numX);
}
