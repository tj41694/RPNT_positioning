#ifndef ADJUSTMENT_H
#define ADJUSTMENT_H

#include "gnssnav_lib.h"

/* parent adjustment functions -------------------------------------------------------------------- */
class adjfunc_c {
/* Constructor */
public:
    adjfunc_c();
    adjfunc_c(int Robust,int RobustCovar);
    virtual ~adjfunc_c();
    /* Implementaion functions */
protected:
    /* least-square functions ----------------------------------------------------- */
    /* least-square: estimate dx -------------------------------------------------- */
    int LSQ_est_dx(const vector<double> &A,const vector<double> &L,vector<double> &R,
        int numL,int numX);
    /* least-squre: estimate dx for variance vector Rvec -------------------------- */
    int LSQ_est_dx_Rvec(const vector<double> &A,const vector<double> &L,
        vector<double> &Rvec,int numL,int numX);
    /* least-squre: update x and Rx ----------------------------------------------- */
    int LSQ_update_xRx(const vector<double> &L,vector<double> &X,
        vector<double> &Rx,int numL,int numX);
    /* least-squre: update x and Rx for variance vector Rvec ---------------------- */
    int LSQ_update_xRx_Rvec(const vector<double> &L,vector<double> &Rvec,
        vector<double> &X,vector<double> &Rx,int numL,int numX);

    /* kalman filter functions ---------------------------------------------------- */
    /* kalman filter: estimate dx-------------------------------------------------- */
    int KF_est_dx(const vector<double> &A,const vector<double> &L,vector<double> &R,
        vector<double> &Rx,int numL,int numX);
    /* kalman filter: estimate dx for variance vector Rvec ------------------------ */
    int KF_est_dx_Rvec(const vector<double> &A,const vector<double> &L,
        vector<double> &Rvec,vector<double> &Rx,int numL,int numX);
    /* kalman filter: update x and Rx --------------------------------------------- */
    int KF_update_xRx(const vector<double> &A,const vector<double> &L,vector<double> &R,
        vector<double> &X,vector<double> &Rx,int numL,int numX);
    /* kalman filter: update x and Rx for variance vector Rvec -------------------- */
    int KF_update_xRx_Rvec(const vector<double> &A,const vector<double> &L,
        vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX);

    /* HVCE for multi-GNSS observation (tight) ------------------------------------ */
    int helmert_tight(const vector<double> &A,int numL,int numX,int numF,
        const int nobs[NSYS][NFREQ*2]);
    /* HVCE for multi-GNSS observation (simple) ----------------------------------- */
    int helmert_simple(const vector<double> &A,int numL,int numX,int numF,
        const int nobs[NSYS][NFREQ*2]);
    /* HVCE for multi-GNSS observation (phase and code) --------------------------- */
    int helmert_pc(const vector<double> &A,int numL,int numX,int numF,
        const int nobs[NSYS][NFREQ*2]);
    /* HVCE for multi-GNSS observation (phase) ------------------------------------ */
    int helmert_phase(const vector<double> &A,int numL,int numX,int numF,
        const int nobs[NSYS][NFREQ*2]);
    /* HVCE for multi-GNSS observation (code) ------------------------------------- */
    int helmert_code(const vector<double> &A,int numL,int numX,int numF,
        const int nobs[NSYS][NFREQ*2]);

    /* robust functions ----------------------------------------------------------- */
    /* robust empty function ------------------------------------------------------ */
    int robust_empty(int numL,int numX,vector<double> &R);
    /* Yang IGG3 robust algorithm with 2 intervals -------------------------------- */
    int robust_yang_IGG3_2s(int numL,int numX,vector<double> &R);
    /* Yang IGG3 robust algorithm with 3 intervals -------------------------------- */
    int robust_yang_IGG3_3s(int numL,int numX,vector<double> &R);
    /* Yang IGG3 robust algorithm with 2 intervals only for R diagonal ------------ */
    int robust_yang_IGG3_2s_Rdiag(int numL,int numX,vector<double> &R);
    /* Yang IGG3 robust algorithm with 3 intervals only for R diagonal ------------ */
    int robust_yang_IGG3_3s_Rdiag(int numL,int numX,vector<double> &R);
    /* Yang IGG3 robust algorithm with 2 intervals for variance vector Rvec ------- */
    int robust_yang_IGG3_2s_Rvec(int numL,int numX,vector<double> &Rvec);
    /* Yang IGG3 robust algorithm with 3 intervals for variance vector Rvec ------- */
    int robust_yang_IGG3_3s_Rvec(int numL,int numX,vector<double> &Rvec);
    /* Yang IGG3 robust algorithm with 2 intervals for each constellation (GNSS) -- */
    int robust_yang_IGG3_2s_GNSS(int numL,int numX,int numF,vector<double> &R,
        const int nobs[NSYS][NFREQ*2]);
    /* Yang IGG3 robust algorithm with 3 intervals for each constellation (GNSS) -- */
    int robust_yang_IGG3_3s_GNSS(int numL,int numX,int numF,vector<double> &R,
        const int nobs[NSYS][NFREQ*2]);
public:
    /* least-square adjustment function ------------------------------------------- */
    int LSQ(const vector<double> &A,const vector<double> &L,vector<double> &R,
        vector<double> &X,vector<double> &Rx,int numL,int numX);
    /* least-square adjustment function for variance vector Rvec ------------------ */
    int LSQ_Rvec(const vector<double> &A,const vector<double> &L,
        vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX);
    /* kalman filter -------------------------------------------------------------- */
    int KF(const vector<double> &A,const vector<double> &L,vector<double> &R,
        vector<double> &X,vector<double> &Rx,int numL,int numX);
    /* kalman filter for variance vector Rvec ------------------------------------- */
    int KF_Rvec(const vector<double> &A,const vector<double> &L,
        vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX);
    /* robust least-square adjustment function ------------------------------------ */
    int robust_LSQ(const vector<double> &A,const vector<double> &L,vector<double> &R,
        vector<double> &X,vector<double> &Rx,int numL,int numX);
    /* robust least-square adjustment function for variance vector Rvec ----------- */
    int robust_LSQ_Rvec(const vector<double> &A,const vector<double> &L,
        vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX);
    /* robust kalman filter ------------------------------------------------------- */
    int robust_KF(const vector<double> &A,const vector<double> &L,vector<double> &R,
        vector<double> &X,vector<double> &Rx,int numL,int numX);
    /* robust kalman filter for variance vector Rvec ------------------------------ */
    int robust_KF_Rvec(const vector<double> &A,const vector<double> &L,
        vector<double> &Rvec,vector<double> &X,vector<double> &Rx,int numL,int numX);
    /* estimate parameter X array and covariance Rx ------------------------------- */
    virtual int adjustment(const vector<double> &A,const vector<double> &L,
        vector<double> &R,vector<double> &X,vector<double> &Rx,
        int numL,int numX,int numF,const int nobs[NSYS][NFREQ*2]);

/* Components */
protected:
    /* adjustment function pointer ------------------------------------------------ */
    /* robust function pointer */
    typedef int (adjfunc_c::*robust_fcp) (int,int,vector<double>&);
    /* adjustment function pointer */
    typedef int (adjfunc_c::*adjust_fcp) (const vector<double>&,const vector<double>&,
        vector<double>&,vector<double>&,vector<double>&,int,int);

    /* normal parameters */
    /* least square process parameters */
    vector<double> P,AP,APL,LP;         /* R^-1, A'*P, A'*P*L, L*P */
    /* kalman filter process parameters */
    vector<double> RxA,Q,E,ERx,KR;      /* Rx*A', Q=A*RxA'+R, E=I, E*Rx, K*R */
    /* Helmert component covariance estimate for multi-GNSS observation */
    vector<double> helmertR;            /* covariance matrix R in helmert process */
    vector<double> newRx,buffRx;        /* corrected and buff covariance matrix Rx */
    int sys_n[NSYS+1];                  /* observation number of different systems */
    int H_flag;                         /* flag of continue helmert algorithm */
    vector<double> sgm,fsgm;            /* original and normalized sigma^2 in one iteration */
    vector<double> sgm_all;             /* all sigma^2 after every iteration */
    /* for robust */
    double buffFactor;                  /* buff factor (robustFactor^1/2) */
    double robust2sK0;                  /* 2-interval IGG3 robust threshold k0 */
    double robust3sK0;                  /* 3-interval IGG3 robust threshold k0 */
    double robust3sK1;                  /* 3-interval IGG3 robust threshold k1 */
    double robust3sK1_1;                /* inverse of 3-interval IGG3 robust threshold k1 (1/k1) */
    double robust3sFactor;              /* 3-interval IGG3 robust factor */

public:
    robust_fcp robustFcp[3][3];         /* robust function pointer (1st-dim: type of R, 2nd-dim: robust mode)
                                         * type of R = 0:var,1:covar,2:vector
                                         * robust mode = 0:none,1:IGG3-2,2:IGG3-3 */
    adjust_fcp adjustFcp;               /* default adjustment function pointer */

    vector<double> V;                   /* residual vector */
    vector<double> dX;                  /* parameter correction value */
    vector<double> K;                   /* gain matrix of Kalman filter */
    vector<double> QQQ;                 /* covariance matrix of estimated parameters */
    double VPV;                         /* sum of weighted observation residual */
    int nsys;                           /* number of systems */
    int H_niter;                        /* Helmert iteration number */
    int largeResFlag;                   /* large residual flag */
    double sigma0_2;                    /* standard variance of unit weight in LSQ */
    int iteraRobust;                    /* robust iteration number */
    int robustMode;                     /* robust mode (0:none,1:IGG3-2,2:IGG3-3) */
    int robustRtype;                    /* covariance type in robust function (0:var,1:covar,2:vector) */
    double sig0;                        /* normalized variance of unit weight for robust */
    vector<double> normV,sortV;         /* normalized residual vector by dividing sqrt(R) */
    vector<double> robustFactor;        /* large residual factor for observations */
    vector<double> sgm2;                /* final unit weight covariance (sigma^2) */
    vector<double> pre_fsgm;            /* previous fsgm */
    vector<string> sys_flag;            /* system flag */
    vector<int> neg_flag;               /* flag of negative (sigma^2) */
    vector<double> precentage;          /* average precentage after each iteration */
    int flag;                           /* adjustment flag */
    string warnMsg,errorMsg;            /* warning and error messages */
};

/* subclass adjustment functions ------------------------------------------------------------------ */
/* least square adjustment function --------------------------------------------------------------- */
class lsadj_c : public adjfunc_c {
/* Constructor */
public:
    lsadj_c();
    lsadj_c(int Robust,int RobustCovar);
    ~lsadj_c();
    /* Implementaion functions */
public:
    /* estimate parameter X array and covariance Rx --------------------------- */
    int adjustment(const vector<double> &A,const vector<double> &L,
        vector<double> &R,vector<double> &X,vector<double> &Rx,
        int numL,int numX,int numF,const int nobs[NSYS][NFREQ*2]);
};

/* kalman filter adjustment function -------------------------------------------------------------- */
class kalmanfilter_c : public adjfunc_c {
/* Constructor */
public:
    kalmanfilter_c();
    kalmanfilter_c(int Robust,int RobustCovar);
    ~kalmanfilter_c();
    /* Implementaion functions */
        /* base functions */
public:
    /* estimate parameter X array and covariance Rx --------------------------- */
    virtual int adjustment(const vector<double> &A,const vector<double> &L,
        vector<double> &R,vector<double> &X,vector<double> &Rx,
        int numL,int numX,int numF,const int nobs[NSYS][NFREQ*2]);
};

/* helmert variance component estimate for kalman filter adjustment function ---------------------- */
class helmert_c : public kalmanfilter_c {
/* Constructor */
public:
    helmert_c();
    helmert_c(int Robust,int RobustCovar);
    ~helmert_c();
    /* Implementaion functions */
    /* for test reset covariance matrix of observations */
    void reset_var(vector<double> &R,int numL,int numF,const int nobs[NSYS][NFREQ*2]);
    /* base functions */
public:
    /* estimate parameter X array and covariance Rx --------------------------- */
    int adjustment(const vector<double> &A,const vector<double> &L,
        vector<double> &R,vector<double> &X,vector<double> &Rx,
        int numL,int numX,int numF,const int nobs[NSYS][NFREQ*2]);
};
#endif
