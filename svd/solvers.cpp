//---------------------------------------------------------------------------
#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
#include "StiffIntegratorT.h"
#include "NonStiffIntegratorT.h"
#include "solvers.h"
//---------------------------------------------------------------------------
using namespace std;
   extern string foc, kin, kinflx, kinc;
DP dxsav;  
int kmax,kount;
Vec_DP *xp_p;
Mat_DP *yp_p;
Vec_INT *ija_p;
Vec_DP *sa_p;
int nrhs;   // counts function evaluations
Vec_DP *x_p;
Mat_DP *d_p;
double Vt;
string kin0;
void res(const double& T, double y[],const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]){
       double f[numx];
	Problem.f(y, f);
	for(unsigned i=0;i<numx;i++) delta[i]=f[i]-yprime[i];
}

void jac(const double& time,  double y[], const double yPrime[], double **PD, double& CJ, const double rPar[], const int iPar[]){}

void ddsolve(const double tout) {
        int i,info[15],idid=0,lrw=1700,liw=1700,iwork[1700],ipar[2],ires=0;
double t=0.,y[numx],yprime[numx],rtol=0.005,atol=1.0e-7, rwork[1700], rpar[2];
        Problem.init();
  for(i=0;i<numx;i++) y[i]=xx[i];
  Problem.f(xx,yprime);
        for(i=0;i<15;i++) info[i]=0;
//        info[1]=1;
//        info[2]=1;
//        info[4]=1;
//        info[10]=1;
	//	printf("infohere ");
	//	printf("%i %i %i %i\n",info[0], info[6], info[10],info[14]);
	//printf("herenow \n");

ddassl_(res,numx,t,xx,yprime,tout,info,rtol,atol,idid,rwork,lrw,iwork, liw,  rpar, ipar, jac);
         Problem.fin(xx);
}

//void d02nefsv(const double tout) {
//int i,info[15],itask,lcom=413990,licom=990,icom[990],iuser[2],ires=0, maxord=5, itol=0, ifail=0;
//double t=0., yprime[Nn], rtol[Nn], atol[Nn], com[413990], ruser[2], h0=0.1e-6, hmax=0.35, to=tout;
//  rtol[0]=0.0032; atol[0]=1.0e-9;
//  char *jceval="Numeric";
//d02mwf_(Nn,maxord,jceval[0],hmax,h0,itol,icom,licom,com,lcom,ifail);
//       horse.ssc(ystart);
//       horse.distr(ystart, yprime);
//d02nef_(Nn,t,to,ystart,yprime,rtol,atol,itask,resd02,jacd02,icom,com,lcom,iuser, ruser,ifail );
//    dlt=horse.dilut();
//	for(unsigned i=0;i<Nn;i++) ystart[i]=y[i];
//}
void resd02(int& neq,double& t,double y[],double ydot[],double r[],int& ires,int *iuser,double *ruser){
        double f[horse.getNn()];
	horse.distr(y, f);
	for(unsigned i=0;i<horse.getNn();i++) r[i]=f[i]-ydot[i];
}
void jacd02(int NEQ,double t,double y[],double YDOT[],double **PD,double& cj,int *iuser,double *ruser){};

void isores(const double& T, double y[],const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]){
       double f[horse.getNn()];
	horse.distr(y, f);
	for(unsigned i=0;i<horse.getNn();i++) delta[i]=f[i]-yprime[i];
}

double Ldistr::ddisolve() {
        int i,info[15],idid=0,lrw=800000,liw=1190,iwork[1190],ipar[2],ires=0,ikin=55;
double t=0.,xi=0.,rtol=0.0075,atol=1.0e-9,h0=0.1e-5, hmax=25.2, rpar[2], *rwork=new double[800000];
        for(i=0;i<15;i++)   info[i]=0;
      info[6]=1; rwork[1]=hmax;//set max step
//      info[7]=1;  rwork[2]=h0; 
//       info[9]=1;//non-negativity constrain
       info[10]=1;//set initial yprime
    double *ys=new double[888],*yprime=new double[888];
        ssc(ys);
          ostringstream foc1, kinet, kicon, flxkin;
    foc1.precision(4); kinet.precision(4);
 massfr();
    showdescr(kinet,expm0);  show(kinet,0);
    showdescr(kicon,expcon);  showcon(kicon,0);
         double *potok=new double [getntime()*nflx]; double *ppotok[getntime()];
         for(int i=0;i<getntime();i++) ppotok[i]=&potok[i*nflx];
 double tout;
    sklad(0);
  for(int j=0;j<nflx;j++) ppotok[0][j]=flx[j]*1000.*dt;
        for(int i=1;i<(getntime());i++) {double tm=(tex[i]-t)/(double)ikin;
        for(int k=0;k<ikin;k++){ tout=t+tm;
        distr(ys, yprime);
ddassl_(isores,getNn(),t,ys,yprime,tout,info,rtol,atol,idid,rwork,lrw,iwork, liw,  rpar, ipar, jac);
cout<<"y[0]="<<ys[35]<<endl;
 if(idid<0) {  throw("dassl problem"); }
massfr();
   show(kinet,tout);  showcon(kicon,tout);
    t=tout; cout<<"t="<<t<<endl;
    }
xi += xits(i);
xi += xicon(i);  sklad(i);
      for(int j=0;j<nflx;j++) ppotok[i][j]=flx[j]*1000.*dt;
    }
wrikin(foc1,getntime());
wricon(foc1,getntime());
  for(int j=0;j<nflx;j++) { flxkin<<Problem.fid[j]<<" "; for(int i=0;i<getntime();i++) flxkin<<ppotok[i][j]<<" "; flxkin<<"\n";}
foc=foc1.str(); kin=kinet.str(); kinflx=flxkin.str(); kinc=kicon.str();
            delete[] potok;
             delete[] rwork; delete[] ys; delete[] yprime;
return xi;}

/**/
void insT::Function(double x, double *y, double *f)
{
	horse.distr(y, f);
}
void isosT::Function(double x, double *y, double *f)
{
	horse.distr(y, f);
}
void isT::Function(double x, double *y, double *f)
{
	Problem.f(y, f);
}
void Jacobian(double x, double *y, double **J){}
void isosT::Jacobian(double x, double *y, double **dfdy){
	int i,j;
//	int Nn=y.size();
        double dydx0[horse.getNn()];
        double dydx1[horse.getNn()];
	double dy,aa;
		Function(x,y,dydx0);
	for ( i=0;i<horse.getNn();i++) { aa=y[i];
	   if(aa>0.0001) dy=aa*0.01; else dy=0.000001;
		y[i] += dy;
		Function(x,y,dydx1);
		y[i] = aa;
		for ( j=0;j<horse.getNn();j++) {
		dfdy[j][i] = (dydx1[j] - dydx0[j]) / dy;
		}
	}
}

void Mass(double **M){} // Mass
void tsolve(const double tmax){
        Problem.init();
	// dimension of problem
//	double yy[numx];
//  for(int i=0;i<numx;i++) yy[i]=xx[i];
	// initial value for x
	 int kmax=25; 
   double xbeg(0.0), dx = tmax/((double)kmax), xend=dx;
	// rtoler and atoler are scalars
	int itoler(0);
	// relative tolerance
	double *rtoler = new double(1.0e-3);
	// absolute tolerance
	double *atoler = new double(1.0e-5);
	// use SolutionOutput routine
	const int iout(0);
	// initial step size
	double hinit(0.01);
	// analytical Jacobian function provided
	 int ijac(0);
	// number of non-zero rows below main diagonal of Jacobian
	int mljac(numx);
	// number of non-zero rows above main diagonal of Jacobian
	int mujac(numx);
	// Mass matrix routine is identity
	const int imas(0);
	int mlmas(numx);
	int mumas(0);
	
	// Use default values (see header files) for these parameters:
	double hmax(0.0);
	int nmax(0);
	double uround(0.0), safe(0.), facl(0.0), facr(0.0);
	int nit(0);
	bool startn(false);
	int nind1(0), nind2(0), nind3(0), npred(0), m1(0), m2(0);
	bool hess(false);
	double fnewt(0.0), quot1(0.0), quot2(0.0), thet(0.0);
	
ostringstream skin; 
	skin<<"t ";
          for (int j=0;j<numx;j++) skin<<Parray::namex[j]<<" "; skin<<"\n0 ";
          for (int j=0;j<numx;j++) skin<<xx[j]<<" "; skin<<endl;
for(int i=0;i<kmax;i++){
	isT stiffT(numx, xx, xbeg, xend, dx, itoler, rtoler, atoler,
		iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,
		mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred,
		m1, m2, hess, fnewt, quot1, quot2, thet);
	stiffT.Integrate();
        skin<<xend<<" ";
	for(int i=0;i<numx; i++) skin<<xx[i]<<" "; skin<<endl;
	xbeg=xend; xend += dx;
	}	
	kin0=skin.str();
         Problem.fin(xx);// cout << pyri <<endl;
/*     for(int i=0;i<nflx;i++) {
         cout<<i<<") "<<Problem.namef[i]<<"= "<<flx[i]*1000.*dt<<endl;
         if(!(i-aldf)) cout<<"aldf="<<flx[aldf]*1000.*dt<<"; aldr="<<flx[aldf+1]*1000.*dt<<endl;
         if(!(i-cs0)) cout<<"cs0="<<flx[cs0]*1000.*dt<<endl;
         if(!(i-coain)) cout<<"coain="<<flx[coain]*1000.*dt<<endl;
         }*/
	delete rtoler;
	delete atoler;
}
void tisolve(const double tmax){
	double y[horse.getNn()];
        horse.ssc(y);
	double xbeg(0.0), xend = tmax, dx(10.0);
	// rtoler and atoler are scalars
	int itoler(0);
	// relative tolerance
	double *rtoler = new double(1.0e-7);
	// absolute tolerance
	double *atoler = new double(1.0e-7);
	// use SolutionOutput routine
	const int iout(1);
	// initial step size
	double hinit(0.0);
	// analytical Jacobian function provided
	 int ijac(0);
	// number of non-zero rows below main diagonal of Jacobian
	int mljac(0);
	// number of non-zero rows above main diagonal of Jacobian
	int mujac(0);
	// Mass matrix routine is identity
	const int imas(0);
	int mlmas(nmet);
	int mumas(0);
	// Use default values (see header files) for these parameters:
	double hmax(0.0);
	int nmax(1200000);
	double uround(0.0), safe(0.), facl(0.0), facr(0.0);
	int nit(0);
	bool startn(false);
	int nind1(0), nind2(0), nind3(0), npred(0), m1(0), m2(0);
	bool hess(false);
	double fnewt(0.0), quot1(0.0), quot2(0.0), thet(0.1);
	
	double beta = 0.03;
	int nstiff = -1;
	int nrdens = horse.getNn();
	unsigned *icont = NULL;
//   cout<<"t0="<<t0<<"; t1="<<t1<<endl;
/*	insT nonstiffT(Nn, y, xbeg, xend, dx, nrdens, itoler, rtoler, atoler, iout, 
	  hinit, hmax, nmax, uround, safe, facl, facr, beta, nstiff, icont);
	nonstiffT.Integrate();
*/	
	isosT stiso(horse.getNn(), y, xbeg, xend, dx, itoler, rtoler, atoler,
		iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,
		mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred,
		m1, m2, hess, fnewt, quot1, quot2, thet);
	stiso.Integrate();
	//for(unsigned i=0;i<Nn;i++) ystart[i]=y[i];

	delete rtoler;
	delete atoler;
}
   
/*

void KinTot::f( REAL t, const vector& x, vector& der ) const{
	const double *py=&x(1); double *pdydt=&der(1);
	Problem.f(py, pdydt);
}
void KinTot::solve(const double tmax){
        vector InitialValue(numx);
        Problem.init();
	for (int i=0;i<numx;i++) InitialValue(i+1)=xx[i];
        IVPParameterBase Parameter(0.,tmax,InitialValue);
        BDF Solver;
        Solver.AssignExpIVP( this, &Parameter );
//	Solver.SetAbsoluteTolerance(0.002);
//	Solver.SetRelativeTolerance(0.002);
//	Solver.SetInitialStepsize( 0.00000001 );
//	Solver.SetMaxStepsizeFactor(0.00015);
	Solver.SetMaxSteps( 3500000 );
        Solver.Integrate();
         Problem.fin();// cout << pyri <<endl;
}

void IsoSlv::f( REAL t, const vector& x, vector& der ) const{
	const DP *py=&x(1); DP *pdydt=&der(1);
	horse.distr(py, pdydt);
}
void Ldistr::integrbs(const double tmax){
        const int Nn=getN();
        IsoSlv marca(Nn);
        vector InitialValue(Nn);
	double* pin=&InitialValue(1);
        ssc(pin);
        IVPParameterBase Parameter(0.,tmax,InitialValue);
        BDF Solver;
    //AdamsMoulton Solver;
    //    DormandPrince853 Solver;
   // DormandPrince5 Solver;
   //   CalvoMontijanoRandez6 Solver;
        Solver.AssignExpIVP( &marca, &Parameter );
	Solver.SetAbsoluteTolerance(0.05);
	Solver.SetRelativeTolerance(0.05);
	Solver.SetInitialStepsize( 0.001 );
	Solver.SetMaxStepsizeFactor(0.003);
	Solver.SetMaxSteps( 350000 );
        Solver.Integrate();
}          

void ode( const DP t,  Vec_I_DP &x, Vec_O_DP &der ) {
	const double *py=&x[0]; double *pdydt=&der[0];
	Problem.f(py, pdydt);
}
void reshit(const double tmax){
        DP eps=1.0e-6,h1=0.0000001,hmin=0.0,x1=0.0;
        Vec_DP ys(numx);
        for (int i=0;i<numx;i++) ys[i]=xx[i];
        nrhs=0;
        int nbad,nok;
        Problem.init();
    NR::odeint(ys,x1,tmax,eps,h1,hmin,nok,nbad,ode,NR::rkqs);
         Problem.fin();// cout << pyri <<endl;
}*/
void  NR::jacobn_s(const DP x, Vec_IO_DP &y, Vec_O_DP &dfdx, Mat_O_DP &dfdy)
{
	int i,j;
//	int Nn=y.size();
	for (i=0;i<horse.getNn();i++) dfdx[i]=0.0; 
        Vec_DP dydx0(horse.getNn());
        Vec_DP dydx1(horse.getNn());
	double dy,aa;
		derivsl(x,y,dydx0);
	for ( i=0;i<horse.getNn();i++) { aa=y[i]; 
		if (y[i]>0.0001) dy= y[i]*0.01;
		  else dy= 0.00001;
		  y[i] += dy;
		derivsl(x,y,dydx1);
		y[i] = aa;
		for ( j=0;j<horse.getNn();j++) {
		dfdy[j][i] = (dydx1[j] - dydx0[j]) / dy;
		}
	}
}

