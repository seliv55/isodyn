#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "modlab.h"
using namespace std;
void Ldistr::ssc(double *pyinit) { setiso(pyinit);
	 for(int i=0;i<nmet;i++) pyinit[i]=xx[i];
for(int i=0;i<(lmet);i++) {met[i]->set0(met[i]->getconc()[0].mean);met[i]->sumt();}// cout<<met[i]->getdescr()<<"="<<(met[i]->getconc())->mean<<"\n";
for(int i=0;i<(lmetk);i++) metk[i]->set0(xx[i+lmet+lmetb]);
met[itrac]->set0(0.);
met[itrac]->iso[markis]=met[itrac]->getconc()[0].mean*marfrac;
met[itrac]->iso[0]=met[itrac]->getconc()[0].mean*(1-marfrac);
massfr();
sklad(0);
 cout<<"Tracer="<<met[itrac]->getdescr()<<" con="<<met[itrac]->getconc()[0].mean<<" marked="<<met[itrac]->iso[markis]/met[itrac]->getconc()[0].mean<<" m0con="<<met[itrac]->iso[0]<<" m0="<<gl.exper[0].shkin(0)<<" sumiso="<<met[itrac]->sumt()<<" calc[0]="<<gl.exper[0].getcalc()[0]/gl.exper[0].getcalc()[gl.exper[0].getmi()+1]<<" calc[mi+1]="<<gl.exper[0].getcalc()[gl.exper[0].getmi()+1]<<" mi="<<gl.exper[0].getmi()<<endl;
}

