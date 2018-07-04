#include <iostream>
#include "nums.hh"
#include "modlab.h"
using namespace std;
void Ldistr::ssc(double *pyinit) {setiso(pyinit);
for(int i=0;i<nmet;i++) pyinit[i]=xx[i];
h6.set0(h6.getconc()[0].mean);
fbp.set0(fbp.getconc()[0].mean);
t3.set0(t3.getconc()[0].mean);
pep.set0(pep.getconc()[0].mean);
pyr.set0(pyr.getconc()[0].mean);
pyrm.set0(pyrm.getconc()[0].mean);
coa.set0(coa.getconc()[0].mean);
coac.set0(coac.getconc()[0].mean);
agl.set0(agl.getconc()[0].mean);
oa.set0(oa.getconc()[0].mean);
oac.set0(oac.getconc()[0].mean);
cit.set0(cit.getconc()[0].mean);
citc.set0(citc.getconc()[0].mean);
akg.set0(akg.getconc()[0].mean);
akgc.set0(akgc.getconc()[0].mean);
fum.set0(fum.getconc()[0].mean);
mal.set0(mal.getconc()[0].mean);
p5.set0(p5.getconc()[0].mean);
e4.set0(e4.getconc()[0].mean);
s7.set0(s7.getconc()[0].mean);
cthf.set0(cthf.getconc()[0].mean);
gl.set0(gl.getconc()[0].mean);
lac.set0(lac.getconc()[0].mean);
glu.set0(glu.getconc()[0].mean);
gln.set0(gln.getconc()[0].mean);
ala.set0(ala.getconc()[0].mean);
asp.set0(asp.getconc()[0].mean);
ser.set0(ser.getconc()[0].mean);
gly.set0(gly.getconc()[0].mean);
pro.set0(pro.getconc()[0].mean);
rna.set0(rna.getconc()[0].mean);
gl.set0(gl.getconc()[0].mean*gl.getexper(0)[0].mean);
gl.iso[48]=gl.getconc()[0].mean*gl.getexper(0)[2].mean;
gl.iso[63]=gl.getconc()[0].mean*gl.getexper(0)[6].mean;
}

