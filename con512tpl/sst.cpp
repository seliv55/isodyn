#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
#include "solvers.h"
using namespace std;

void fcn(const int& numx,double x[],double fvec[],int* iuser,double* ruser,int& iflag){
	Problem.f(x, fvec);
}
void fcn1(const int& numx,double x[],double fvec[],int& iflag){
	Problem.f(x, fvec);
}

 void rho_(){};
 void rhoa_(){};
 void rhojac_(){};
 void rhojs_(){};
 void fjacs_(){};

void f_(double x[],double fvec[]){
	Problem.f(x, fvec);
}
void fjac_(double x[],double y[],int& k){
  double dx,aa, dydt[numx], dydt1[numx];
	Problem.f(x, dydt);
	 aa=x[k];
      if(aa>0.0001) dx=aa*0.01; else dx=0.000001;
	x[k] += dx;
	 Problem.f(x,dydt1);
	  x[k] = aa;
	for (int j=0;j<numx;j++) y[j] = (dydt1[j]-dydt[j]) / dx;
}

inline void homology(){
  int TRACE=0,iflg(-1),NFE,NDIMA; bool POLY_SWITCH;
   double  ARCTOL=1e-10,EPS=0.5e-11,ARCRE=0.5e-9, ARCAE=0.5e-9, ANSRE=1.0e-10, ANSAE=1.0e-10, A[nrea],ARCLEN;
    double  SSP4[4]={1e-7,1e-3,0.,0.},SSP8[8]={0.,0.,0.,1e-7,1e-3,0.,0.,0.},Y[numx+1];
     for(int i=0;i<numx;i++) Y[i+1]=xx[i];
   __hompack90_MOD_fixpdf(&numx,Y,&iflg,&ARCTOL,&EPS,&TRACE,A,&NDIMA, &NFE, &ARCLEN);
   __hompack90_MOD_fixpqf(&numx,Y,&iflg,&ARCRE,&ARCAE,&ANSRE,&ANSAE,&TRACE,A,SSP4, &NFE, &ARCLEN);
//   __hompack90_MOD_fixpnf(&numx,Y,&iflg,&ARCRE,&ARCAE,&ANSRE,&ANSAE,&TRACE,A,SSP8, &NFE, &ARCLEN,&POLY_SWITCH);
    for(int i=0;i<numx;i++) xx[i]=Y[i+1];
}

void Fit::setsst(){
 int ifail(0),lwa=numx*(3*numx+14)/2+1,iuser[numx];
  double xtol(1e-11), ruser[numx], wa[lwa], der[numx],ma(0.),mi(0.);
//     homology();
//   c05qbf_(fcn,numx,xx,der,xtol,iuser,ruser,ifail);
   hybrd1_(fcn1,numx,xx,der,xtol,ifail,wa,lwa);
    f(xx,der); 
 for(int i=0;i<numx;i++) {if(der[i]>ma) ma=der[i]; else if(der[i]<mi) mi=der[i];}
          cout<<"deriv: max="<<ma<<"; min="<<mi<<endl;
  }

double Fit::jacobian(double *y){
  double dy,aa, dydt[numx], bol;
    double *dx= new double[numx*numx], *dfdy[numx];
    		for(int i=0;i<numx;i++) dfdy[i]=&dx[i*numx]; 
  setsst();
    for (int i=0;i<numx;i++) { aa=y[i];
      if(aa>0.0001) dy=aa*0.01; else dy=0.000001;
	y[i] += dy;
	 f(y,dydt);
	  y[i] = aa;
	for (int j=0;j<numx;j++) dfdy[j][i] = (dydt[j]) / dy;
	}
	bol=eigen(dx);
		delete[] dx;
 return bol;}

double Fit::eigen(double dx[]){
 Mat_DP a(dx,numx,numx);
  complex<double> *sobst=new complex<double> [numx];
   Vec_CPLX_DP wri(sobst,numx);
     NR::balanc(a);
       NR::elmhes(a);
        NR::hqr(a,wri); double bol=wri[0].real();
        for (int i=0;i<numx;i++){if(wri[i].real()>bol) bol=wri[i].real();}
          cout << "max real eigenvalue: " << bol<<endl;
 return bol;}

void Fit::chpar(int ipar,int ovar,double factor,int iter)  {
  double a = rea[ipar].v();
   for(int i=0;i<iter;i++) {a*=factor;
    rea[ipar].setVm(a);
     jacobian(xx);
      cout<<"x="<<xx[ovar]<<"; par="<<rea[ipar].v()<<endl;
       }
   }

