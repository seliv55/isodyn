#include <iostream>
#include "nums.hh"
#include "modlab.h"
using namespace std;
void Ldistr::ssc(double *pyinit) {setiso(pyinit);
	 for(int i=0;i<numx;i++) pyinit[i]=xx[i];
for(int i=0;i<(lmet);i++) {met[i]->set0(met[i]->getconc()[0].mean);}// cout<<met[i]->getdescr()<<"="<<(met[i]->getconc())->mean<<"\n";
for(int i=0;i<(lmetb);i++) metb[i]->set0(xx[i+lmet]);
for(int i=0;i<(lmetk);i++) metk[i]->set0(xx[i+lmet+lmetb]);
met[itrac]->set0(met[itrac]->getconc()[0].mean);
met[itrac]->iso[markis]=met[itrac]->getconc()[0].mean*marfrac;
met[itrac]->iso[0]=met[itrac]->getconc()[0].mean*(1-marfrac);
 cout<<"Tracer="<<met[itrac]->getdescr()<<" fract="<<marfrac<<" marked="<<met[itrac]->iso[markis]<<" m0="<<met[itrac]->iso[0]<<endl;
}

