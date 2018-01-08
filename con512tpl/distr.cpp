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
	h6.sett(); p5.sett(); s7.sett(); t3.sett(); e4.sett();
	Problem.f(py,pdydt);
	Problem.ff(py,pdydt);
	Problem.fin(py);/**/
	double xthf=thft-cthf.sumt();
gl.input(h6,fluxes[hk]/gl.sumt());
h6.input(fbp,fluxes[pfk],fluxes[fbpase]);
t3.input(pep,fluxes[t3pep],fluxes[pept3]);
pep.input(pyr,fluxes[pk]);
pyr.input(lac,fluxes[pyrlac],fluxes[lacpyr]);
pyr.input(pyrm,fluxes[pyrdcm],fluxes[pyrdmc]);
pyrm.cutfirst(coa,fluxes[pdh]);
cit.icdh(akg,fluxes[citakg]);
akg.decarb(fum,fluxes[akgfum]);
fum.input(mal,fluxes[fumal],fluxes[malfum]);
mal.input(oa,fluxes[maloa],fluxes[oamal]);
pyrm.carb(oa,fluxes[pc]);
mal.decarb(pyrm,fluxes[malicm]);
oac.decarb(pyr,fluxes[malicc]);
h6.cutfirst(p5,fluxes[ppp]);
oac.input(mal,fluxes[oacd],fluxes[mald]);
cit.input(citc,fluxes[citdmc],fluxes[citdcm]);
akg.input(akgc,fluxes[akgdmc],fluxes[akgdcm]);
citc.split(coac,oac,fluxes[coaout]);
citc.icdh(akgc,fluxes[citakg1]);
akgc.icdhr(citc,fluxes[akgcit1]);
gln.input(akgc,fluxes[gln_in],fluxes[gln_out]);
glu.input(akgc,fluxes[gluin],fluxes[gluout]);
t3.input(ser,fluxes[t3ser]);
ser.input(pyrm,fluxes[serpyr]);
oac.input(asp,fluxes[asp_o],fluxes[asp_i]);
pyr.input(ala,fluxes[ala_o],fluxes[ala_i]);
p5.input(rna,fluxes[r5_o],fluxes[r5_i]);
pyrm.diso[0]+=fluxes[cystin];
pro.input(akgc,fluxes[proin],fluxes[proout]);
akg.diso[0]+=fluxes[kgin];
coa.diso[0]+=fluxes[coain];
gln.output(fluxes[gln_pr]);
ser.output(fluxes[ser_pr]);
asp.output(fluxes[asp_pr]);
ala.output(fluxes[ala_pr]);
pro.output(fluxes[pro_pr]);
ala.diso[0]+=fluxes[trpala];
ser.split(gly,cthf,fluxes[sergly]*xthf);
ser.condence(gly,cthf,fluxes[glyser]);
cit.condence(oa,coa,fluxes[cs0]);
fbp.splinverse(t3,t3,fluxes[rald],fluxes[rald+1]);
//h6.split(t3,dhe,(fluxes[tafl]+fluxes[tafl+2])/h6.sumt());
//s7.condence(e4,dhe,(fluxes[tafl]+fluxes[tafl+3])/e4.sumt()/dhe.sumt());
//s7.split(e4,dhe,(fluxes[tafl+1]+fluxes[tafl+3])/s7.sumt());
//h6.condence(t3,dhe,(fluxes[tafl+1]+fluxes[tafl+2])/t3.sumt()/dhe.sumt());

p5.split(t3,gae,(fluxes[tkfl]+fluxes[p5f6]+fluxes[p5g3i])/p5.sumt());
s7.condence(p5,gae,(fluxes[tkfl]+fluxes[f6s7]+fluxes[s7p5i])/p5.sumt()/gae.sumt());
s7.split(p5,gae,(fluxes[s7p5]+fluxes[s7f6]+fluxes[s7p5i])/s7.sumt());
p5.condence(t3,gae,(fluxes[s7p5]+fluxes[f6p5]+fluxes[p5g3i])/t3.sumt()/gae.sumt());
h6.split(e4,gae,(fluxes[f6p5]+fluxes[f6s7]+fluxes[f6e4i])/h6.sumt());
h6.condence(e4,gae,(fluxes[s7f6]+fluxes[p5f6]+fluxes[f6e4i])/e4.sumt()/gae.sumt());
//spInvsl(fbp.iso,fbp.diso,t3.iso,t3.diso,fluxes[rald+2],t3.sumt());
	h6.tka(t3,e4,s7,fluxes[tafl]);
	s7.tka(e4,t3,h6,fluxes[tafl+1]);
	h6.invista(t3,fluxes[tafl+2]);
	s7.invista(e4,fluxes[tafl+3]);
//p5.tkk(t3,p5,s7,fluxes[tkfl]);
//s7.tkk(p5,t3,p5,fluxes[tkfl+1]);
//h6.tkk(e4,t3,p5,fluxes[tkfl+2]);
//p5.tkk(t3,e4,h6,fluxes[tkfl+3]);
//h6.tkk(e4,p5,s7,fluxes[tkfl+4]);
//s7.tkk(p5,e4,h6,fluxes[tkfl+5]);
//p5.invistk(t3,fluxes[tkfl+6]);
//h6.invistk(e4,fluxes[tkfl+7]);
//s7.invistk(p5,fluxes[tkfl+8]);
	gl.volume(Vt);
	lac.volume(Vt);
	glu.volume(Vt);
	gln.volume(Vt);
	ala.volume(Vt);
	asp.volume(Vt);
	ser.volume(Vt);
	gly.volume(Vt);
	pro.volume(Vt);
	rna.volume(Vt);
	symm(fum.getisot());
//xx[nh6]=h6.sumt(); xx[nfbp]=fbp.sumt(); xx[nt3]=t3.sumt(); xx[npep]=pep.sumt(); xx[npyr]=pyr.sumt(); xx[npyrm]=pyrm.sumt(); xx[ncoa]=coa.sumt(); xx[ncoac]=coac.sumt(); xx[nagl]=agl.sumt(); xx[noa]=oa.sumt(); xx[noac]=oac.sumt(); xx[ncit]=cit.sumt(); xx[ncitc]=citc.sumt(); xx[nakg]=akg.sumt(); xx[nakgc]=akgc.sumt(); xx[nfum]=fum.sumt(); xx[nmal]=mal.sumt(); xx[np5]=p5.sumt(); xx[ne4]=e4.sumt(); xx[ns7]=s7.sumt(); xx[ncthf]=cthf.sumt(); xx[ndhe]=dhe.sumt(); xx[ngae]=gae.sumt(); xx[ngl]=gl.sumt(); xx[nlac]=lac.sumt(); xx[nglu]=glu.sumt(); xx[ngln]=gln.sumt(); xx[nala]=ala.sumt(); xx[nasp]=asp.sumt(); xx[nser]=ser.sumt(); xx[ngly]=gly.sumt(); xx[npro]=pro.sumt(); xx[nrna]=rna.sumt(); 
xx[nh6]=py[nh6]; xx[nfbp]=py[nfbp]; xx[nt3]=py[nt3]; xx[npep]=py[npep]; xx[npyr]=py[npyr]; xx[npyrm]=py[npyrm]; xx[ncoa]=py[ncoa]; xx[ncoac]=py[ncoac]; xx[nagl]=py[nagl]; xx[noa]=py[noa]; xx[noac]=py[noac]; xx[ncit]=py[ncit]; xx[ncitc]=py[ncitc]; xx[nakg]=py[nakg]; xx[nakgc]=py[nakgc]; xx[nfum]=py[nfum]; xx[nmal]=py[nmal]; xx[np5]=py[np5]; xx[ne4]=py[ne4]; xx[ns7]=py[ns7]; xx[ncthf]=py[ncthf]; xx[ndhe]=py[ndhe]; xx[ngae]=py[ngae]; xx[ngl]=py[ngl]; xx[nlac]=py[nlac]; xx[nglu]=py[nglu]; xx[ngln]=py[ngln]; xx[nala]=py[nala]; xx[nasp]=py[nasp]; xx[nser]=py[nser]; xx[ngly]=py[ngly]; xx[npro]=py[npro]; xx[nrna]=py[nrna]; 
}
