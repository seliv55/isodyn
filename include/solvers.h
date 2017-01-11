//---------------------------------------------------------------------------
#ifndef mainH
extern "C" void d02mwf_(int& NEQ,int& MAXORD,char& JCEVAL,double& HMAX,double& H0,int& ITOL,int ICOM[],int& LICOM,double COM[],int& LCOM,int& IFAIL);
extern "C" void d02nef_(int& NEQ,double& T,double& TOUT,double Y[],double YDOT[],double *RTOL,double *ATOL,int& ITASK,
void (*RES)(int& NEQ,double& T,double Y[],double YDOT[],double R[],int& IRES,int *IUSER,double *RUSER),
void (*JAC)(int NEQ,double T,double Y[],double YDOT[],double **PD,double& CJ,int *IUSER,double *RUSER),
int ICOM[],double COM[],int& LCOM,int *IUSER,double *RUSER,int& IFAIL);

extern "C" void resd02(int& NEQ,double& T,double Y[],double YDOT[],double R[],int& IRES,int *IUSER,double *RUSER);
extern void jacd02(int NEQ,double T,double Y[],double YDOT[],double **PD,double& CJ,int *IUSER,double *RUSER);

extern "C" void c05qbf_(void (*fcn)(const int& n,double x[],double fvec[],int* iuser,double* ruser,int& iflag),
  const int& N,double X[],double FVEC[],double& XTOL,int* IUSER,double* RUSER,int& IFAIL);
extern void fcn(const int& n,double x[],double fvec[],int* iuser,double* ruser,int& iflag);

extern "C" void hybrd1_(void (*fcn1)(const int& n,double x[],double fvec[],int& iflag),
  const int& N,double X[],double FVEC[],double& XTOL,int& info,double wa[],int& lwa);
extern void fcn1(const int& n,double x[],double fvec[],int& iflag);

extern "C" {
void __hompack90_MOD_fixpqf(const int* N,double* Y,int* IFLAG,double* ARCRE,double* ARCAE,double* ANSRE,double* ANSAE,int* TRACE,double* A,double* SSPAR,int* NFE,double* ARCLEN);

void __hompack90_MOD_fixpnf(const int* N,double* Y,int* IFLAG,double* ARCRE,double* ARCAE,double* ANSRE,double* ANSAE,int* TRACE,double* A,double* SSPAR,int* NFE,double* ARCLEN,bool* POLY_SWITCH);

void __hompack90_MOD_fixpdf(const int* N,double* Y,int* IFLAG,double* ARCTOL,double *EPS,int* TRACE,double* A,int* NDIMA,int* NFE,double* ARCLEN);

 void f_(double x[],double fvec[]);
 void fjac_(double x[],double fvec[], int& k);
 void rho_();
 void rhoa_();
 void rhojac_();
 void rhojs_();
 void fjacs_();
 }
//      SUBROUTINE FIXPNF(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,SSPAR,NFE,ARCLEN,POLY_SWITCH)

extern "C" void ddassl_(
void (*res)(const double& time, double y[], const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]),
const int& noOfEquations,
double& currentTime,
double initialY[],
double initialYPrime[],
const double& finalTime,
int info[15],
double &relativeTolerance,
double &absoluteTolerance,
int& outputStatusFlag,
double dWorkArray[],
const int& lengthOfDWork,
int iWorkArray[],
const int& lengthOfIWork,
const double rParArray[],
const int iParArray[]
,void (*jacobian)(const double& time,  double y[], const double yprime[], double **PD, double& CJ,const double rPar[],const int iPar[]));
 extern void res(const double& T, double y[],const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]);
    extern void jac(const double& time,  double y[], const double yprime[], double **PD, double& CJ,const double rPar[],const  int iPar[]);
    extern void isores(const double& T, double y[],const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]);

extern void tsolve(const double tmax);
extern void tisolve(const double tmax);
extern void ddsolve(const double tmax);
extern void d02nefsv(const double tmax);
extern void derivs(const DP x, Vec_IO_DP &y, Vec_O_DP &dydx);
extern void integ(const double tmax);
extern void derivsl(const DP x, Vec_IO_DP &y, Vec_O_DP &dydx);
extern void integrbs(const double tmax);
//extern void  NR::jacobn_s(const DP x, Vec_I_DP &y, Vec_O_DP &dfdx, Mat_O_DP &dfdy);
extern double setgrad(double xi0,const double tmax); 
extern double fitonerand(const double a,double tmax);
extern void fit(const double tmax);
extern void coord(const double tmax);
extern void grad(double fact, double tmax) ;
extern void reshit(const double tmax);
extern double getmax();
extern DP dxsav;   // defining declarations
extern int kmax,kount;
extern Vec_DP *xp_p;
extern Mat_DP *yp_p;
extern Vec_INT *ija_p;
extern Vec_DP *sa_p;
extern int nrhs;   // counts function evaluations
extern Vec_DP *x_p;
extern Mat_DP *d_p;
extern Ldistr horse;
#define mainH
#endif
