#include <iostream>
#include <cmath>
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
using namespace std;
double Vi, xribi, xasp=1.,mu;
double Ldistr::readExp (string fn) {
string aaa;  ifstream fi(fn.c_str()); double Ti,ts1;  mu=0.; fi>> dt; getline(fi,aaa);
tex[0]=0.;  for(ntime=1;;ntime++) {fi>>tex[ntime]; tex[ntime] *= 60.; if(tex[ntime]<0) break;}
fi>> Vi; getline(fi,aaa); double Nc[ntime]; for(int i=0;i<ntime;i++) fi>>Nc[i]; getline(fi,aaa);
for(int i=1;i<ntime;i++) mu += log(Nc[i]/Nc[0])/tex[i];
 mu /= (double)(ntime-1); getline(fi,aaa);
h6.setconc(xx[nh6]);
fbp.setconc(xx[nfbp]);
t3.setconc(xx[nt3]);
pep.setconc(xx[npep]);
pyr.setconc(xx[npyr]);
pyrm.setconc(xx[npyrm]);
coa.setconc(xx[ncoa]);
coac.setconc(xx[ncoac]);
agl.setconc(xx[nagl]);
oa.setconc(xx[noa]);
oac.setconc(xx[noac]);
cit.setconc(xx[ncit]);
citc.setconc(xx[ncitc]);
akg.setconc(xx[nakg]);
akgc.setconc(xx[nakgc]);
fum.setconc(xx[nfum]);
mal.setconc(xx[nmal]);
p5.setconc(xx[np5]);
e4.setconc(xx[ne4]);
s7.setconc(xx[ns7]);
cthf.setconc(xx[ncthf]);
gl.setconc(xx[ngl]);
lac.setconc(xx[nlac]);
glu.setconc(xx[nglu]);
gln.setconc(xx[ngln]);
ala.setconc(xx[nala]);
asp.setconc(xx[nasp]);
ser.setconc(xx[nser]);
gly.setconc(xx[ngly]);
pro.setconc(xx[npro]);
rna.setconc(xx[nrna]);
	gl.readc(fi,ntime);
	lac.readc(fi,ntime);
pyr.read(fi,ntime);
agl.read(fi,ntime);
oac.read(fi,ntime);
cit.read(fi,ntime);
akg.read(fi,ntime);
fum.read(fi,ntime);
mal.read(fi,ntime);
gl.read(fi,ntime);
gln.read(fi,ntime);
fi.close();
return ts1;}
