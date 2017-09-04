//---------------------------------------------------------------------------
#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
#include "solvers.h"
//---------------------------------------------------------------------------
using namespace std;
    string foc, kin, kinflx, kinc;
void derivsl(const DP x, Vec_IO_DP &y, Vec_O_DP &dydx){
	double *py=&y[0]; 
	DP *pdydt=&dydx[0];
        nrhs++; 
        Vt=Vi*exp(mu*x);
	horse.distr(py, pdydt);
}

double Ldistr::integrbs(){
  DP eps=1.0e-6,h1=0.00001,hmin=1.0e-11,x1=0.0, tm,xi=0., xii;
   const int KMAX(2),icol(nflx);  Vec_DP yy(Nn);  DP *pyinit = &yy[0];
     ssc(pyinit);
      xp_p=new Vec_DP(KMAX); yp_p=new Mat_DP(Nn,KMAX);
        Vec_DP &xp=*xp_p;  Mat_DP &yp=*yp_p;
        double *potok=new double [ntime*icol]; double *ppotok[ntime];
         for(int i=0;i<ntime;i++) ppotok[i]=&potok[i*icol];
    int nbad,nok; nrhs=0; kmax=KMAX;
          ostringstream foc1, kinet, kicon, flxkin;
    foc1.precision(4); kinet.precision(4);
 massfr();
    showdescr(kinet,expm0);  show(kinet,0);
    showdescr(kicon,expcon);  showcon(kicon,0);
 double xfin; sklad(0);
  for(int j=0;j<nflx;j++) ppotok[0][j]=flx[j]*1000.*dt;
        for(int i=1;i<(ntime);i++) {tm=(tex[i]-x1)/10.;
        dxsav = tm/((double)(KMAX-1));
        for(int k=0;k<10;k++){ xfin=x1+tm;
    NR::odeint(yy,x1,xfin,eps,h1,hmin,nok,nbad,derivsl,NR::rkqs);
massfr();   show(kinet,xfin);  showcon(kicon,xfin);
    x1=xfin;
    }
xi += xits(i);
xi += xicon(i);  sklad(i);
     x1=tex[i];  for(int j=0;j<nflx;j++) ppotok[i][j]=flx[j]*1000.*dt;
    }
wrikin(foc1,ntime);
wricon(foc1,ntime);
  for(int j=0;j<nflx;j++) { flxkin<<Problem.fid[j]<<" "; for(int i=0;i<ntime;i++) flxkin<<ppotok[i][j]<<" "; flxkin<<"\n";}
foc=foc1.str(); kin=kinet.str(); kinflx=flxkin.str(); kinc=kicon.str();
        delete yp_p;
        delete xp_p;
        delete[] potok;
 return xi;}

      void Ldistr::setfige()   {
  ostringstream fige1,fige2;
//      readExp(fex2); wriconex(fige2);
//        ofstream fo("excon"); fo<<fige2.str(); fo.close();      
//      readExp(fex1); wrim0ex(fige1);
//          fo.open("exm0"); fo<<fige1.str(); fo.close();
      }
//Email: demin@insysbio.ru Tel: +7 910 4449284
