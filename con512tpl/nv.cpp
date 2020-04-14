#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
#include "solvers.h"
#include "analis.h"
using namespace std;
const int hk=0, pyrclac=hk+1, lacin=pyrclac+1, laccpyr=lacin+1, pyrclacc=laccpyr+1, pdh=pyrclacc+1, citakg=pdh+1, akgsuc=citakg+1, sucmal=akgsuc+1, pc=sucmal+1, malicm=pc+1, oacd=malicm+1, akgdcm=oacd+1, liase=akgdcm+1, akgcit1=liase+1, gln_in=akgcit1+1, ala_o=gln_in+1, cs0=ala_o+1, D=cs0+1, rdt=D+1, nrea=rdt+1;

const int nflx=nrea;
const int npyrc=0, npyr=npyrc+1, ncoa=npyr+1, noac=ncoa+1, ncit=noac+1, nakg=ncit+1, nakgc=nakg+1, nsuc=nakgc+1, nmal=nsuc+1, nlacc=nmal+1, ngl=nlacc+1, nlac=ngl+1, ngln=nlac+1, nmet=ngln+1;

const int numx=ngl;

Metab Ldistr::pyrc(3,"pyr"), Ldistr::pyr(3,"Pyr"), Ldistr::coa(2,"CoA"), Ldistr::oac(4,"oac"), Ldistr::cit(6,"Cit"), Ldistr::akg(5,"aKg"), Ldistr::akgc(5,"akgc"), Ldistr::suc(4,"suc"), Ldistr::mal(4,"Mal"), Ldistr::lacc(3,"lacc"), Ldistr::gl(3,"Gluc"), Ldistr::lac(3,"Lac"), Ldistr::gln(5,"Glutamin");
	Fit Problem;
	const double thft(1.);
	double xx[nmet],flx[nflx],fluxes[nflx];
	double xinit1[nmet],xinit2[nmet];
	string Parray::fid[nflx],Parray::fname[nflx],Parray::fschem[nflx], Parray::namex[nmet];
	Reapar Parray::rea[nrea];
	double Analis::nv1[nrea], Analis::nv2[nrea];
Metab *Ldistr::met[13];

 void Ldistr::setmet(){met[0]=&pyrc; met[1]=&pyr; met[2]=&coa; met[3]=&oac; met[4]=&cit; met[5]=&akg; met[6]=&akgc; met[7]=&suc; met[8]=&mal; met[9]=&lacc; met[10]=&gl; met[11]=&lac; met[12]=&gln; 
 lmet=13;  }
 void Ldistr::setcon(){met[0]->setconc(xx[npyrc]); met[1]->setconc(xx[npyr]); met[2]->setconc(xx[ncoa]); met[3]->setconc(xx[noac]); met[4]->setconc(xx[ncit]); met[5]->setconc(xx[nakg]); met[6]->setconc(xx[nakgc]); met[7]->setconc(xx[nsuc]); met[8]->setconc(xx[nmal]); met[9]->setconc(xx[nlacc]); met[10]->setconc(xx[ngl]); met[11]->setconc(xx[nlac]); met[12]->setconc(xx[ngln]);  }
void Fit::f(const double *y,double *dydx) {
	for(int i=0;i<nmet;i++) dydx[i]=0.;
	for(int i=0;i<nflx;i++) flx[i]=0.;
flx[hk]= rea[hk].v(); 	dydx[npyrc] += flx[hk];  
flx[pyrclac]= rea[pyrclac].v(y[npyrc]); 	dydx[npyrc] -= flx[pyrclac];  
flx[lacin]= rea[lacin].v(); 	dydx[nlacc] += flx[lacin];  
flx[laccpyr]= rea[laccpyr].v(y[nlacc]); 	dydx[nlacc] -= flx[laccpyr];  dydx[npyr] += flx[laccpyr];  
flx[pyrclacc]= rea[pyrclacc].v(y[npyrc]); 	dydx[npyrc] -= flx[pyrclacc];  dydx[nlacc] += flx[pyrclacc];  
flx[pdh]= rea[pdh].v(y[npyr]); 	dydx[npyr] -= flx[pdh];  dydx[ncoa] += flx[pdh];  
flx[citakg]= rea[D].v()*rea[citakg].v(y[ncit]); 	dydx[ncit] -= flx[citakg];  dydx[nakg] += flx[citakg];  
flx[akgsuc]= rea[D].v()*rea[akgsuc].v(y[nakg]); 	dydx[nakg] -= flx[akgsuc];  dydx[nsuc] += flx[akgsuc];  
flx[sucmal]= rea[D].v()*rea[sucmal].v(y[nsuc]); 	dydx[nsuc] -= flx[sucmal];  dydx[nmal] += flx[sucmal];  
flx[pc]= rea[pc].v(y[npyr]); 	dydx[npyr] -= flx[pc];  dydx[nmal] += flx[pc];  
flx[malicm]= rea[malicm].v(y[nmal]); 	dydx[nmal] -= flx[malicm];  dydx[npyr] += flx[malicm];  
flx[oacd]= rea[oacd].v(y[noac]); 	dydx[noac] -= flx[oacd];  dydx[nmal] += flx[oacd];  
flx[akgdcm]= rea[akgdcm].v(y[nakgc]); 	dydx[nakgc] -= flx[akgdcm];  dydx[nakg] += flx[akgdcm];  
flx[liase]= rea[liase].v(y[ncit]); 	dydx[ncit] -= flx[liase];  dydx[noac] += flx[liase];  
flx[akgcit1]= rea[akgcit1].v(y[nakgc]); 	dydx[nakgc] -= flx[akgcit1];  dydx[ncit] += flx[akgcit1];  
flx[gln_in]= rea[gln_in].v(); 	dydx[nakgc] += flx[gln_in];  
flx[ala_o]= rea[ala_o].v(y[npyrc]); 	dydx[npyrc] -= flx[ala_o];  
flx[cs0]= rea[D].v()*rea[cs0].v(y[nmal], y[ncoa]); 	dydx[nmal] -= flx[cs0];  dydx[ncoa] -= flx[cs0];  dydx[ncit] += flx[cs0];  
flx[D]= rea[D].v(); 	
flx[rdt]= rea[rdt].v(); 	
for(int i=0;i<nmet;i++) dydx[i]*=(flx[rdt]/Vi);
}

void Fit::ff(const double *y,double *dydx) {
	f(y,dydx);
	dydx[ngl] = (-flx[hk])*flx[rdt];
	dydx[nlac] = (-flx[lacin])*flx[rdt];
	dydx[ngln] = (-flx[gln_in])*flx[rdt];
}

void Parray::init(){ft3=10.; fh6=7.;}
void Parray::fin(double y[]){
	
flfor(y);
}
void Parray::flfor(double *y){
for(int i=0;i<nflx;i++) fluxes[i] = flx[i] * flx[rdt]/Vi;
fluxes[pyrclac] /= y[npyrc];
fluxes[laccpyr] /= y[nlacc];
fluxes[pyrclacc] /= y[npyrc];
fluxes[pdh] /= y[npyr];
fluxes[citakg] /= y[ncit];
fluxes[akgsuc] /= y[nakg];
fluxes[sucmal] /= y[nsuc];
fluxes[pc] /= y[npyr];
fluxes[malicm] /= y[nmal];
fluxes[oacd] /= y[noac];
fluxes[akgdcm] /= y[nakgc];
fluxes[liase] /= y[ncit];
fluxes[akgcit1] /= y[nakgc];
fluxes[ala_o] /= y[npyrc];
fluxes[cs0] /= y[nmal]*y[ncoa];
}
