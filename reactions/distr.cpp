#include <iostream>
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
using namespace std;
void Ldistr::distr(double *py,double *pdydt) {
	double NOL=0.;
	setiso(py); setdiso(pdydt);
	for (int i=0;i<Nn;i++) pdydt[i]=0.;
	Problem.ff(py,pdydt);
	Problem.fin(py);/**/
gl.input(pyrc,fluxes[hk]/gl.sumt());
pyrc.output(fluxes[pyrclac]);
lac.input(lacc,fluxes[lacin]);
lacc.input(pyr,fluxes[laccpyr]);
pyrc.input(lacc,fluxes[pyrclacc]);
pyrc.input(pyr,fluxes[pyrdcm]);
pyr.cutfirst(coa,fluxes[pdh]);
cit.icdh(akg,fluxes[citakg]);
akg.decarb(suc,fluxes[akgsuc]);
suc.input(mal,fluxes[sucmal]);
pyr.carb(mal,fluxes[pc]);
mal.decarb(pyr,fluxes[malicm]);
oac.input(mal,fluxes[oacd]);
akgc.input(akg,fluxes[akgdcm]);
cit.split(oac,coa,fluxes[liase]);
akgc.icdhr(cit,fluxes[akgcit1]);
gln.input(akgc,fluxes[gln_in]);
pyrc.output(fluxes[ala_o]);
cit.condence(mal,coa,fluxes[cs0]);
gl.volume(Vt);
	lac.volume(Vt);
	gln.volume(Vt);
	symm(mal.getisot());
//xx[npyrc]=pyrc.sumt(); xx[npyr]=pyr.sumt(); xx[ncoa]=coa.sumt(); xx[noac]=oac.sumt(); xx[ncit]=cit.sumt(); xx[nakg]=akg.sumt(); xx[nakgc]=akgc.sumt(); xx[nsuc]=suc.sumt(); xx[nmal]=mal.sumt(); xx[nlacc]=lacc.sumt(); xx[ngl]=gl.sumt(); xx[nlac]=lac.sumt(); xx[ngln]=gln.sumt(); 
//xx[npyrc]=py[npyrc]; xx[npyr]=py[npyr]; xx[ncoa]=py[ncoa]; xx[noac]=py[noac]; xx[ncit]=py[ncit]; xx[nakg]=py[nakg]; xx[nakgc]=py[nakgc]; xx[nsuc]=py[nsuc]; xx[nmal]=py[nmal]; xx[nlacc]=py[nlacc]; xx[ngl]=py[ngl]; xx[nlac]=py[nlac]; xx[ngln]=py[ngln]; 
}
