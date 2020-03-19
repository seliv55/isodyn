#include <iostream>
#include "nums.hh"
#include "modlab.h"
using namespace std;
void Ldistr::ssc(double *pyinit) { setiso(pyinit);
	 for(int i=0;i<nmet;i++) pyinit[i]=xx[i];
  for(int i=0;i<(lmet);i++) {met[i]->set0(met[i]->getconc()[0].mean); met[i]->sumt();}
  met[itrac]->set0(0.);
  met[itrac]->iso[markis]=met[itrac]->getconc()[0].mean*marfrac;
  met[itrac]->iso[0]=met[itrac]->getconc()[0].mean*(1-marfrac);
  }

