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
        Vt=Vi;//*exp(mu*x);
	horse.distr(py, pdydt);
}

double Ldistr::integrbs(){
  DP eps=1.0e-6,h1=0.00001,hmin=1.0e-11,t=0.0, tm,xi=0., xii;
     Vec_DP yy(Nn); DP *pyinit = &yy[0];
   const int KMAX(2),icol(nflx);
      xp_p=new Vec_DP(KMAX); yp_p=new Mat_DP(Nn,KMAX);
        Vec_DP &xp=*xp_p;  Mat_DP &yp=*yp_p;
        double *potok=new double [ntime*icol]; double *ppotok[ntime];
         for(int i=0;i<ntime;i++) ppotok[i]=&potok[i*icol];
    int nbad,nok; nrhs=0; kmax=KMAX;
          ostringstream foc1, kinet, kicon, flxkin;
    foc1.precision(4); kinet.precision(4);
    foc1<<"* Metabolite t: ";
    for(int i=0;i<(ntime);i++) foc1<<setw(10)<<tex[i]; foc1<<setw(12)<<" : exper <-> xi²";
   ssc(pyinit);
    massfr();
     showdescr(kinet,expm0);  show(kinet,0);
     showdescr(kicon,expcon);  showcon(kicon,0);
 double tout; sklad(0);
  for(int j=0;j<nflx;j++) ppotok[0][j]=flx[j]*1000.*flx[rdt];
         setiso(pyinit);
        for(int i=1;i<(ntime);i++) {tm=(tex[i]-t)/10.;
        dxsav = tm/((double)(KMAX-1));
        for(int k=0;k<10;k++){ tout=t+tm;
    NR::odeint(yy,t,tout,eps,h1,hmin,nok,nbad,derivsl,NR::rkqs);
  massfr();   show(kinet,tout);
    showcon(kicon,tout); 
    t=tout;
    } sklad(i);
xi += xits(i);
xi += xicon(i); 
      for(int j=0;j<nflx;j++) ppotok[i][j]=flx[j]*1000.*flx[rdt];
    }
wrikin(foc1,ntime);
wricon(foc1,ntime);
  for(int j=0;j<nflx;j++) { flxkin<<Problem.fid[j]<<" "; for(int i=0;i<ntime;i++) flxkin<<ppotok[i][j]<<" "; flxkin<<'\n';}
foc=foc1.str(); kin=kinet.str(); kinflx=flxkin.str(); kinc=kicon.str();// cout<<kin<<endl;
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
