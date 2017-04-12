#include <iostream>
#include "nums.hh"
#include "modlab.h"
using namespace std;
void Ldistr::ssc(double *pyinit) {setiso(pyinit);
	 for(int i=0;i<nmet;i++) pyinit[i]=xx[i];
gl.set0(gl.getconc()[0].mean);
glycog.set0(glycog.getconc()[0].mean);
lac.set0(lac.getconc()[0].mean);
glu.set0(glu.getconc()[0].mean);
glu25.set0(glu.getconc()[0].mean);
gln.set0(gln.getconc()[0].mean);
ala.set0(xx[nala]);
coac.set0(coac.getconc()[0].mean);
gly.set0(gly.getconc()[0].mean);
ser.set0(xx[nser]);
rna.set0(ser.getconc()[0].mean);
agl.set0(1.1);
asp.set0(asp.getconc()[0].mean);
pro.set0(pro.getconc()[0].mean);
gln.set0(gln.getconc()[0].mean*gln.getexper(0)[0].mean);
gln.iso[31]=gln.getconc()[0].mean*gln.getexper(0)[5].mean;
gl.set0(gl.getconc()[0].mean*gl.getexper(0)[0].mean); 
gl.iso[1] =gl.getconc()[0].mean*gl.getexper(0)[1].mean;
gl.iso[48]=gl.getconc()[0].mean*gl.getexper(0)[2].mean;
gl.iso[49]=gl.getconc()[0].mean*gl.getexper(0)[3].mean/4.;
gl.iso[50]=gl.iso[49];
gl.iso[52]=gl.iso[49];
gl.iso[56]=gl.iso[49];
gl.iso[51]=gl.getconc()[0].mean*gl.getexper(0)[4].mean/6.;
gl.iso[53]=gl.iso[51];
gl.iso[54]=gl.iso[51];
gl.iso[57]=gl.iso[51];
gl.iso[58]=gl.iso[51];
gl.iso[60]=gl.iso[51];
gl.iso[63]=gl.getconc()[0].mean*gl.getexper(0)[6].mean;
h6.set0(xx[nh6]);
fbp.set0(xx[nfbp]);
t3.set0(xx[nt3]);
s7.set0(xx[ns7]);
pep.set0(xx[npep]);
pyr.set0(xx[npyr]);
pyrm.set0(xx[npyrm]);
coa.set0(xx[ncoa]);
cthf.set0(xx[ncthf]);
oa.set0(xx[noa]);
oac.set0(xx[noac]);
cit.set0(xx[ncit]);
citc.set0(xx[ncitc]);
akg.set0(xx[nakg]);
akgc.set0(xx[nakgc]);
fum.set0(xx[nfum]);
mal.set0(xx[nmal]);
p5.set0(xx[np5]);
e4.set0(xx[ne4]);
}

