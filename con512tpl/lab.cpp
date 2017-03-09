#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "modlab.h"
using namespace std;

int Ldistr::getN() {
gl.ny=nmet;
glycog.ny=gl.ny+gl.getlen();
lac.ny=glycog.ny+glycog.getlen();
glu.ny=lac.ny+lac.getlen();
glu25.ny=glu.ny;
gln.ny=glu.ny+glu.getlen();
ala.ny=gln.ny+gln.getlen();
ser.ny=ala.ny+ala.getlen();
agl.ny=ser.ny+ser.getlen();
rna.ny=agl.ny+agl.getlen();
asp.ny=rna.ny+rna.getlen();
pro.ny=asp.ny+asp.getlen();
h6.ny=pro.ny+pro.getlen();
fbp.ny=h6.ny+h6.getlen();
t3.ny=fbp.ny+fbp.getlen();
s7.ny=t3.ny+t3.getlen();
pep.ny=s7.ny+s7.getlen();
pyr.ny=pep.ny+pep.getlen();
pyrm.ny=pyr.ny+pyr.getlen();
coa.ny=pyrm.ny+pyrm.getlen();
coac.ny=coa.ny+coa.getlen();
gly.ny=coac.ny+coac.getlen();
cthf.ny=gly.ny+gly.getlen();
oa.ny=cthf.ny+cthf.getlen();
oac.ny=oa.ny+oa.getlen();
cit.ny=oac.ny+oac.getlen();
citc.ny=cit.ny+cit.getlen();
akg.ny=citc.ny+citc.getlen();
akgc.ny=akg.ny+akg.getlen();
fum.ny=akgc.ny+akgc.getlen();
mal.ny=fum.ny+fum.getlen();
p5.ny=mal.ny+mal.getlen();
e4.ny=p5.ny+p5.getlen();
return (e4.ny+e4.getlen());
}

void Ldistr::setdiso(double *pyinit) {
gl.diso= &pyinit[gl.ny];
glycog.diso= &pyinit[glycog.ny];
lac.diso= &pyinit[lac.ny];
glu.diso= &pyinit[glu.ny];
glu25.diso= &pyinit[glu.ny];
gln.diso= &pyinit[gln.ny];
ala.diso= &pyinit[ala.ny];
ser.diso= &pyinit[ser.ny];
agl.diso= &pyinit[agl.ny];
rna.diso= &pyinit[rna.ny];
asp.diso= &pyinit[asp.ny];
pro.diso= &pyinit[pro.ny];
h6.diso= &pyinit[h6.ny];
fbp.diso= &pyinit[fbp.ny];
t3.diso= &pyinit[t3.ny];
s7.diso= &pyinit[s7.ny];
pep.diso= &pyinit[pep.ny];
pyr.diso= &pyinit[pyr.ny];
pyrm.diso= &pyinit[pyrm.ny];
coa.diso= &pyinit[coa.ny];
coac.diso= &pyinit[coac.ny];
gly.diso= &pyinit[gly.ny];
cthf.diso= &pyinit[cthf.ny];
oa.diso= &pyinit[oa.ny];
oac.diso= &pyinit[oac.ny];
cit.diso= &pyinit[cit.ny];
citc.diso= &pyinit[citc.ny];
akg.diso= &pyinit[akg.ny];
akgc.diso= &pyinit[akgc.ny];
fum.diso= &pyinit[fum.ny];
mal.diso= &pyinit[mal.ny];
p5.diso= &pyinit[p5.ny];
e4.diso= &pyinit[e4.ny];
}

void Ldistr::setiso(double *pyinit) {
gl.iso= &pyinit[gl.ny];
glycog.iso= &pyinit[glycog.ny];
lac.iso= &pyinit[lac.ny];
glu.iso= &pyinit[glu.ny];
glu25.iso= &pyinit[glu.ny];
gln.iso= &pyinit[gln.ny];
ala.iso= &pyinit[ala.ny];
ser.iso= &pyinit[ser.ny];
agl.iso= &pyinit[agl.ny];
rna.iso= &pyinit[rna.ny];
asp.iso= &pyinit[asp.ny];
pro.iso= &pyinit[pro.ny];
h6.iso= &pyinit[h6.ny];
fbp.iso= &pyinit[fbp.ny];
t3.iso= &pyinit[t3.ny];
s7.iso= &pyinit[s7.ny];
pep.iso= &pyinit[pep.ny];
pyr.iso= &pyinit[pyr.ny];
pyrm.iso= &pyinit[pyrm.ny];
coa.iso= &pyinit[coa.ny];
coac.iso= &pyinit[coac.ny];
gly.iso= &pyinit[gly.ny];
cthf.iso= &pyinit[cthf.ny];
oa.iso= &pyinit[oa.ny];
oac.iso= &pyinit[oac.ny];
cit.iso= &pyinit[cit.ny];
citc.iso= &pyinit[citc.ny];
akg.iso= &pyinit[akg.ny];
akgc.iso= &pyinit[akgc.ny];
fum.iso= &pyinit[fum.ny];
mal.iso= &pyinit[mal.ny];
p5.iso= &pyinit[p5.ny];
e4.iso= &pyinit[e4.ny];
}

void Ldistr::massfr() {
gl.percent(); glycog.percent(); lac.percent(); glu.percent(1,1); glu25.percent(1,0); gln.percent(); ala.percent(); ser.percent(); agl.percent(); rna.percent(); asp.percent(); pro.percent(); pyr.percent(); pyrm.percent(); oa.percent(); oac.percent(); cit.percent(); citc.percent(); akg.percent(); fum.percent(); mal.percent(); coac.percent(); gly.percent(); cthf.percent(); }

double Ldistr::xits(int its) {int itp=its-1; double xi=0;
for(int i=0;i<lmet;i++) xi += met[i]-> chisq(its,met[i]->getmi());
return xi;}

void Ldistr::wriconex(ostringstream& fo) {
//t=1;  glc=t+1; glsd=glc+1;  lacc=glsd+1; lacsd=lacc+1;
  data *b;
 fo<<"time "; for(int i=0;i<lmet;i++) fo<<met[i]->getdescr()<<" sd "; fo<<endl;
 for(int its=0;its<ntime;its++) { fo<<tex[its];
   for(int i=0;i<lmet;i++){ b=met[i]->getconc();   fo<<" "<<b[its].mean<<" "<<b[its].sd; }
      fo<<"\n";         }
}

void Ldistr::wrim0ex(ostringstream& fo) {
//t=1;  lac=t+1; lacsd=lac+1; glu=lacsd+1; glusd=glu+1; ala=glusd+1; alasd=ala=+1; gly=alasd+1; glysd=gly+1; ser=glysd+1; sersd=ser+1; pro=sersd+1; prosd=pro+1; cit=prosd+1; citsd=cit+1; agl=citsd+1; aglsd=agl+1; asp=aglsd+1; aspsd=asp+1; fum=aspsd+1; fumsd=fum+1; mal=fumsd+1; malsd=mal+1; coa=malsd+1; coasd=coa+1;
  data *b;
fo<<"time "; for(int i=0;i<lmet;i++) fo<<met[i]->getdescr()<<" sd "; fo<<endl;
   for(int its=0;its<ntime;its++) {  fo<<tex[its];
       for(int i=0;i<lmet;i++){ b=met[i]->getexper(its);   fo<<" "<<b[0].mean<<" "<<b[0].sd; }
      fo<<"\n";         }
       }
void Ldistr::show(ostringstream& fo,double xfin) { fo<<setw(5)<<xfin;
for(int i=0;i<lmet;i++) met[i]-> showm0(fo);}

//cglct=t+1; clact=cglct+1; lact=clact+1; glut=lact+1; alat=glut+1; glyt=alat+1; sert=glyt+1; prot=sert+1; citmt=prot+1; citct=citmt+1; citsum=citct+1; aglt=citsum+1; aspt=aglt+1; fumt=aspt+1; malt=fumt+1; oact=malt+1; malsum=oact+1;
//void Ldistr::showmi(ostringstream& fo,int nt) { 
//  fo<<"Ser ";  ser.showmi(fo,nt);}


double Ldistr::label() {
return (h6.sumt()+fbp.sumt()+t3.sumt()+s7.sumt()+pep.sumt()+pyr.sumt()+pyrm.sumt()+coa.sumt()+oa.sumt()+oac.sumt()+cit.sumt()+citc.sumt()+akg.sumt()+akgc.sumt()+fum.sumt()+mal.sumt()+p5.sumt()+e4.sumt());}

void Ldistr::flback(){ cout<<" ∑isotopomers-variable:"<<"\n";
cout<<" h6="<<h6.sumt()<<"-"<<xx[nh6]<<" fbp="<<fbp.sumt()<<"-"<<xx[nfbp]<<" t3="<<t3.sumt()<<"-"<<xx[nt3]<<" s7="<<s7.sumt()<<"-"<<xx[ns7]<<" pep="<<pep.sumt()<<"-"<<xx[npep]<<" pyr="<<pyr.sumt()<<"-"<<xx[npyr]<<" pyrm="<<pyrm.sumt()<<"-"<<xx[npyrm]<<" coa="<<coa.sumt()<<"-"<<xx[ncoa]<<" oa="<<oa.sumt()<<"-"<<xx[noa]<<" oac="<<oac.sumt()<<"-"<<xx[noac]<<" cit="<<cit.sumt()<<"-"<<xx[ncit]<<" citc="<<citc.sumt()<<"-"<<xx[ncitc]<<" akg="<<akg.sumt()<<"-"<<xx[nakg]<<" akgc="<<akgc.sumt()<<"-"<<xx[nakgc]<<" fum="<<fum.sumt()<<"-"<<xx[nfum]<<" mal="<<mal.sumt()<<"-"<<xx[nmal]<<" p5="<<p5.sumt()<<"-"<<xx[np5]<<" e4="<<e4.sumt()<<"-"<<xx[ne4]<<" gly="<<gly.sumt()<<"-"<<xx[ngly]<<" thf="<<cthf.sumt()<<"-"<<xx[ncthf]<<" atp="<<xx[natp]<<" nad="<<xx[nnad]<<"\n";}

