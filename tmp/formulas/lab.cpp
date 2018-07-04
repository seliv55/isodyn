#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "modlab.h"
using namespace std;

int Ldistr::getN() {
h6.ny=nmet;
fbp.ny=h6.ny+h6.getlen();
t3.ny=fbp.ny+fbp.getlen();
pep.ny=t3.ny+t3.getlen();
pyr.ny=pep.ny+pep.getlen();
pyrm.ny=pyr.ny+pyr.getlen();
coa.ny=pyrm.ny+pyrm.getlen();
coac.ny=coa.ny+coa.getlen();
agl.ny=coac.ny+coac.getlen();
oa.ny=agl.ny+agl.getlen();
oac.ny=oa.ny+oa.getlen();
cit.ny=oac.ny+oac.getlen();
citc.ny=cit.ny+cit.getlen();
akg.ny=citc.ny+citc.getlen();
akgc.ny=akg.ny+akg.getlen();
fum.ny=akgc.ny+akgc.getlen();
mal.ny=fum.ny+fum.getlen();
p5.ny=mal.ny+mal.getlen();
e4.ny=p5.ny+p5.getlen();
s7.ny=e4.ny+e4.getlen();
cthf.ny=s7.ny+s7.getlen();
gl.ny=cthf.ny+cthf.getlen();
lac.ny=gl.ny+gl.getlen();
glu.ny=lac.ny+lac.getlen();
gln.ny=glu.ny+glu.getlen();
ala.ny=gln.ny+gln.getlen();
asp.ny=ala.ny+ala.getlen();
ser.ny=asp.ny+asp.getlen();
gly.ny=ser.ny+ser.getlen();
pro.ny=gly.ny+gly.getlen();
rna.ny=pro.ny+pro.getlen();
return (rna.ny+rna.getlen());
}

void Ldistr::setdiso(double *pyinit) {
h6.diso= &pyinit[h6.ny];
fbp.diso= &pyinit[fbp.ny];
t3.diso= &pyinit[t3.ny];
pep.diso= &pyinit[pep.ny];
pyr.diso= &pyinit[pyr.ny];
pyrm.diso= &pyinit[pyrm.ny];
coa.diso= &pyinit[coa.ny];
coac.diso= &pyinit[coac.ny];
agl.diso= &pyinit[agl.ny];
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
s7.diso= &pyinit[s7.ny];
cthf.diso= &pyinit[cthf.ny];
gl.diso= &pyinit[gl.ny];
lac.diso= &pyinit[lac.ny];
glu.diso= &pyinit[glu.ny];
gln.diso= &pyinit[gln.ny];
ala.diso= &pyinit[ala.ny];
asp.diso= &pyinit[asp.ny];
ser.diso= &pyinit[ser.ny];
gly.diso= &pyinit[gly.ny];
pro.diso= &pyinit[pro.ny];
rna.diso= &pyinit[rna.ny];
}

void Ldistr::setiso(double *pyinit) {
h6.iso= &pyinit[h6.ny];
fbp.iso= &pyinit[fbp.ny];
t3.iso= &pyinit[t3.ny];
pep.iso= &pyinit[pep.ny];
pyr.iso= &pyinit[pyr.ny];
pyrm.iso= &pyinit[pyrm.ny];
coa.iso= &pyinit[coa.ny];
coac.iso= &pyinit[coac.ny];
agl.iso= &pyinit[agl.ny];
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
s7.iso= &pyinit[s7.ny];
cthf.iso= &pyinit[cthf.ny];
gl.iso= &pyinit[gl.ny];
lac.iso= &pyinit[lac.ny];
glu.iso= &pyinit[glu.ny];
gln.iso= &pyinit[gln.ny];
ala.iso= &pyinit[ala.ny];
asp.iso= &pyinit[asp.ny];
ser.iso= &pyinit[ser.ny];
gly.iso= &pyinit[gly.ny];
pro.iso= &pyinit[pro.ny];
rna.iso= &pyinit[rna.ny];
}

void Ldistr::massfr() {
pyr.percent(); agl.percent(); oac.percent(); cit.percent(); akg.percent(); fum.percent(); mal.percent(); gl.percent(); lac.percent(); gln.percent(); }

double Ldistr::xits(int its) {int itp=its-1; double xi=0;
xi += pyr.chisq(its);
xi += agl.chisq(its);
xi += oac.chisq(its);
xi += cit.chisq(its);
xi += akg.chisq(its);
xi += fum.chisq(its);
xi += mal.chisq(its);
xi += gl.chicon(its);
xi += gl.chisq(its);
xi += lac.chicon(its);
xi += gln.chisq(its);
return xi;}

void Ldistr::wriex(ostringstream& fo) {
 for(int its=0;its<ntime;its++) { fo<<tex[its];
  fo<<" "<<pyr.getexper(its)[0].mean<<" "<<pyr.getexper(its)[0].sd;
  fo<<" "<<agl.getexper(its)[0].mean<<" "<<agl.getexper(its)[0].sd;
  fo<<" "<<oac.getexper(its)[0].mean<<" "<<oac.getexper(its)[0].sd;
  fo<<" "<<cit.getexper(its)[0].mean<<" "<<cit.getexper(its)[0].sd;
  fo<<" "<<akg.getexper(its)[0].mean<<" "<<akg.getexper(its)[0].sd;
  fo<<" "<<fum.getexper(its)[0].mean<<" "<<fum.getexper(its)[0].sd;
  fo<<" "<<mal.getexper(its)[0].mean<<" "<<mal.getexper(its)[0].sd;
  fo<<" "<<gl.getconc()[its].mean<<" "<<gl.getconc()[its].sd;
  fo<<" "<<gl.getexper(its)[0].mean<<" "<<gl.getexper(its)[0].sd;
  fo<<" "<<lac.getconc()[its].mean<<" "<<lac.getconc()[its].sd;
  fo<<" "<<gln.getexper(its)[0].mean<<" "<<gln.getexper(its)[0].sd;fo<<'\n';}}

void Ldistr::show(ostringstream& fo,double xfin) { static int itit=0;
if(itit==0) fo<<"time pyr_m_0 agl_m_0 oac_m_0 cit_m_0 akg_m_0 fum_m_0 mal_m_0 gl_c gl_m_0 lac_c gln_m_0 \n"; itit++;
fo<<setw(5)<<xfin;pyr.showm0(fo); agl.showm0(fo); oac.showm0(fo); cit.showm0(fo); akg.showm0(fo); fum.showm0(fo); mal.showm0(fo); gl.showcon(fo); gl.showm0(fo); lac.showcon(fo); gln.showm0(fo); fo<<'\n';}

double Ldistr::label() {
return (h6.sumt()+fbp.sumt()+t3.sumt()+pep.sumt()+pyr.sumt()+pyrm.sumt()+coa.sumt()+coac.sumt()+agl.sumt()+oa.sumt()+oac.sumt()+cit.sumt()+citc.sumt()+akg.sumt()+akgc.sumt()+fum.sumt()+mal.sumt()+p5.sumt()+e4.sumt()+s7.sumt()+cthf.sumt());}

void Ldistr::flback(){ cout<<" ∑isotopomers-variable:"<<"\n";
cout<<" nh6="<<h6.sumt()<<"-"<<xx[nh6]<<" nfbp="<<fbp.sumt()<<"-"<<xx[nfbp]<<" nt3="<<t3.sumt()<<"-"<<xx[nt3]<<" npep="<<pep.sumt()<<"-"<<xx[npep]<<" npyr="<<pyr.sumt()<<"-"<<xx[npyr]<<" npyrm="<<pyrm.sumt()<<"-"<<xx[npyrm]<<" ncoa="<<coa.sumt()<<"-"<<xx[ncoa]<<" ncoac="<<coac.sumt()<<"-"<<xx[ncoac]<<" nagl="<<agl.sumt()<<"-"<<xx[nagl]<<" noa="<<oa.sumt()<<"-"<<xx[noa]<<" noac="<<oac.sumt()<<"-"<<xx[noac]<<" ncit="<<cit.sumt()<<"-"<<xx[ncit]<<" ncitc="<<citc.sumt()<<"-"<<xx[ncitc]<<" nakg="<<akg.sumt()<<"-"<<xx[nakg]<<" nakgc="<<akgc.sumt()<<"-"<<xx[nakgc]<<" nfum="<<fum.sumt()<<"-"<<xx[nfum]<<" nmal="<<mal.sumt()<<"-"<<xx[nmal]<<" np5="<<p5.sumt()<<"-"<<xx[np5]<<" ne4="<<e4.sumt()<<"-"<<xx[ne4]<<" ns7="<<s7.sumt()<<"-"<<xx[ns7]<<" ncthf="<<cthf.sumt()<<"-"<<xx[ncthf]<<" ngl="<<gl.sumt()<<"-"<<xx[ngl]<<" nlac="<<lac.sumt()<<"-"<<xx[nlac]<<" nglu="<<glu.sumt()<<"-"<<xx[nglu]<<" ngln="<<gln.sumt()<<"-"<<xx[ngln]<<" nala="<<ala.sumt()<<"-"<<xx[nala]<<" nasp="<<asp.sumt()<<"-"<<xx[nasp]<<" nser="<<ser.sumt()<<"-"<<xx[nser]<<" ngly="<<gly.sumt()<<"-"<<xx[ngly]<<" npro="<<pro.sumt()<<"-"<<xx[npro]<<" nrna="<<rna.sumt()<<"-"<<xx[nrna]<<"\n";}

