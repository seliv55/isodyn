#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
#include "solvers.h"
#include "analis.h"
using namespace std;
const int hk=0, pfk=hk+1, fbpase=pfk+1, t3pep=fbpase+1, pept3=t3pep+1, pk=pept3+1, pyrlac=pk+1, lacpyr=pyrlac+1, pyrdcm=lacpyr+1, pyrdmc=pyrdcm+1, pdh=pyrdmc+1, citakg=pdh+1, akgfum=citakg+1, fumal=akgfum+1, malfum=fumal+1, maloa=malfum+1, oamal=maloa+1, pc=oamal+1, malicm=pc+1, malicc=malicm+1, ppp=malicc+1, oacd=ppp+1, mald=oacd+1, citdmc=mald+1, citdcm=citdmc+1, akgdmc=citdcm+1, akgdcm=akgdmc+1, coaout=akgdcm+1, citakg1=coaout+1, akgcit1=citakg1+1, gln_in=akgcit1+1, gln_out=gln_in+1, gluin=gln_out+1, gluout=gluin+1, t3ser=gluout+1, serpyr=t3ser+1, asp_o=serpyr+1, asp_i=asp_o+1, ala_o=asp_i+1, ala_i=ala_o+1, r5_o=ala_i+1, r5_i=r5_o+1, cystin=r5_i+1, proin=cystin+1, proout=proin+1, kgin=proout+1, coain=kgin+1, gln_pr=coain+1, ser_pr=gln_pr+1, asp_pr=ser_pr+1, ala_pr=asp_pr+1, pro_pr=ala_pr+1, trpala=pro_pr+1, mthf=trpala+1, thf=mthf+1, sergly=thf+1, glyser=sergly+1, cs0=glyser+1, D=cs0+1, atpase=D+1, resp=atpase+1, rald=resp+1, rta=rald+1, rtk=rta+1, nrea=rtk+1;

const int aldfl=rald, aldrev=aldfl+1, aldfli=aldrev+1, aldi1=aldfli+1, tafl=aldi1+1, s7f6a=tafl+1, f6g3a=s7f6a+1, s7e4a=f6g3a+1, tkfl=s7e4a+1, s7p5=tkfl+1, f6p5=s7p5+1, p5f6=f6p5+1, f6s7=p5f6+1, s7f6=f6s7+1, p5g3i=s7f6+1, f6e4i=p5g3i+1, s7p5i=f6e4i+1, nflx=s7p5i+1;

const int nfbp=0, nt3=nfbp+1, npep=nt3+1, npyr=npep+1, npyrm=npyr+1, ncoa=npyrm+1, ncoac=ncoa+1, nagl=ncoac+1, noa=nagl+1, noac=noa+1, ncit=noac+1, ncitc=ncit+1, nakg=ncitc+1, nakgc=nakg+1, nfum=nakgc+1, nmal=nfum+1, ne4=nmal+1, ns7=ne4+1, nh6=ns7+1, np5=nh6+1, ncthf=np5+1, ngae=ncthf+1, ndhe=ngae+1, n_atp=ndhe+1, n_nad=n_atp+1, ngl=n_nad+1, nlac=ngl+1, nglu=nlac+1, ngln=nglu+1, nala=ngln+1, nasp=nala+1, nser=nasp+1, ngly=nser+1, npro=ngly+1, nrna=npro+1, nmet=nrna+1;

const int numx=ngl;

Metab_data Ldistr::fbp(6,"fbp"), Ldistr::t3(3,"t3"), Ldistr::pep(3,"pep"), Ldistr::pyr(3,"pyr"), Ldistr::pyrm(3,"Pyr"), Ldistr::coa(2,"CoA"), Ldistr::coac(2,"coac"), Ldistr::agl(3,"Glycerol"), Ldistr::oa(4,"Oaa"), Ldistr::oac(4,"oac"), Ldistr::cit(6,"Cit"), Ldistr::citc(6,"citc"), Ldistr::akg(5,"aKg"), Ldistr::akgc(5,"akgc"), Ldistr::fum(4,"Fum"), Ldistr::mal(4,"Mal"), Ldistr::e4(4,"ne4"), Ldistr::gl(6,"Gluc"), Ldistr::lac(3,"Lac"), Ldistr::glu(5,"Glutamate2-5"), Ldistr::gln(5,"Glutamin"), Ldistr::ala(3,"Ala"), Ldistr::asp(4,"Asp"), Ldistr::ser(3,"Ser"), Ldistr::gly(2,"Gly"), Ldistr::pro(5,"Pro"), Ldistr::rna(5,"Rib");
ketose Ldistr::s7(7,"s7"), Ldistr::h6(6,"h6"), Ldistr::p5(5,"rib");
Metab Ldistr::cthf(1,"cthf"), Ldistr::gae(2,"gae"), Ldistr::dhe(3,"dhe");
	Fit Problem;
	const double thft(1.);
	double dt,xx[nmet],flx[nflx],fluxes[nflx];
	double xinit1[nmet],xinit2[nmet];
	string Parray::fid[nflx],Parray::fname[nflx],Parray::fschem[nflx], Parray::namex[numx];
	Reapar Parray::rea[nrea];
	double Analis::nv1[nrea], Analis::nv2[nrea];
Metab_data *Ldistr::met[27];
Metab *Ldistr::metb[3];
 ketose *Ldistr::metk[3];

 void Ldistr::setmet(){met[0]=&fbp; met[1]=&t3; met[2]=&pep; met[3]=&pyr; met[4]=&pyrm; met[5]=&coa; met[6]=&coac; met[7]=&agl; met[8]=&oa; met[9]=&oac; met[10]=&cit; met[11]=&citc; met[12]=&akg; met[13]=&akgc; met[14]=&fum; met[15]=&mal; met[16]=&e4; met[17]=&gl; met[18]=&lac; met[19]=&glu; met[20]=&gln; met[21]=&ala; met[22]=&asp; met[23]=&ser; met[24]=&gly; met[25]=&pro; met[26]=&rna; metb[0]=&cthf; metb[1]=&gae; metb[2]=&dhe; metk[0]=&s7; metk[1]=&h6; metk[2]=&p5; 
 lmet=27; lmetb=3; lmetk=3;  }
 void Ldistr::setcon(){met[0]->setconc(xx[nfbp]); met[1]->setconc(xx[nt3]); met[2]->setconc(xx[npep]); met[3]->setconc(xx[npyr]); met[4]->setconc(xx[npyrm]); met[5]->setconc(xx[ncoa]); met[6]->setconc(xx[ncoac]); met[7]->setconc(xx[nagl]); met[8]->setconc(xx[noa]); met[9]->setconc(xx[noac]); met[10]->setconc(xx[ncit]); met[11]->setconc(xx[ncitc]); met[12]->setconc(xx[nakg]); met[13]->setconc(xx[nakgc]); met[14]->setconc(xx[nfum]); met[15]->setconc(xx[nmal]); met[16]->setconc(xx[ne4]); met[17]->setconc(xx[ngl]); met[18]->setconc(xx[nlac]); met[19]->setconc(xx[nglu]); met[20]->setconc(xx[ngln]); met[21]->setconc(xx[nala]); met[22]->setconc(xx[nasp]); met[23]->setconc(xx[nser]); met[24]->setconc(xx[ngly]); met[25]->setconc(xx[npro]); met[26]->setconc(xx[nrna]); metb[0]->setconc(xx[ncthf]); metb[1]->setconc(xx[ngae]); metb[2]->setconc(xx[ndhe]); metk[0]->setconc(xx[ns7]); metk[1]->setconc(xx[nh6]); metk[2]->setconc(xx[np5]);  }
void Fit::f(const double *y,double *dydx) {
	for(int i=0;i<numx;i++) dydx[i]=0.;
	 for(int i=0;i<nflx;i++) flx[i]=0.;
	 double amp = -(sqrt(4.*xx[n_atp]*tan-3.*xx[n_atp]*xx[n_atp])-2.*tan+xx[n_atp])/2.;
	 double a_dp = (sqrt(xx[n_atp])*sqrt(4.*tan-3.*xx[n_atp])-xx[n_atp])/2.;
	 double h_nad = tnad-xx[n_nad];
	 xthf=thft-xx[ncthf];
flx[hk]= rea[hk].v( y[n_atp]);                     dydx[n_atp] -= flx[hk];  dydx[nh6] += flx[hk];  
flx[pfk]= rea[pfk].v( y[nh6], y[n_atp]);           dydx[nh6] -= flx[pfk];  dydx[n_atp] -= flx[pfk];  dydx[nfbp] += flx[pfk];  
flx[fbpase]= rea[fbpase].v( y[nfbp]);              dydx[nfbp] -= flx[fbpase];  dydx[nh6] += flx[fbpase];  
flx[t3pep]= rea[t3pep].v( y[nt3], y[n_nad], a_dp); dydx[nt3] -= flx[t3pep];  dydx[n_nad] -= flx[t3pep];  dydx[n_atp] += flx[t3pep];  dydx[npep] += flx[t3pep];  
flx[pept3]= rea[pept3].v( y[npep], y[n_atp], h_nad); dydx[npep] -= flx[pept3];  dydx[n_atp] -= flx[pept3];  dydx[n_nad] += flx[pept3];  dydx[nt3] += flx[pept3];  
flx[pk]= rea[pk].v( y[npep], a_dp);                dydx[npep] -= flx[pk];  dydx[n_atp] += flx[pk];  dydx[npyr] += flx[pk];  
flx[pyrlac]= rea[pyrlac].v( y[npyr], h_nad);       dydx[npyr] -= flx[pyrlac];  dydx[n_nad] += flx[pyrlac];  
flx[lacpyr]= rea[lacpyr].v( y[n_nad], xx[nlac]);   dydx[n_nad] -= flx[lacpyr];  dydx[npyr] += flx[lacpyr];  
flx[pyrdcm]= rea[pyrdcm].v( y[npyr]);              dydx[npyr] -= flx[pyrdcm];  dydx[npyrm] += flx[pyrdcm];  
flx[pyrdmc]= rea[pyrdmc].v( y[npyrm]);             dydx[npyrm] -= flx[pyrdmc];  dydx[npyr] += flx[pyrdmc];  
flx[pdh]= rea[pdh].v( y[npyrm], y[n_nad]);         dydx[npyrm] -= flx[pdh];  dydx[n_nad] -= flx[pdh];  dydx[ncoa] += flx[pdh];  
flx[citakg]= rea[D].v()*rea[citakg].v( y[ncit], y[n_nad]); dydx[ncit] -= flx[citakg];  dydx[n_nad] -= flx[citakg];  dydx[nakg] += flx[citakg];  
flx[akgfum]= rea[D].v()*rea[akgfum].v( y[nakg], y[n_nad], a_dp); dydx[nakg] -= flx[akgfum];  dydx[n_nad] -= flx[akgfum];  dydx[n_atp] += flx[akgfum];  dydx[nfum] += flx[akgfum];  
flx[fumal]= rea[D].v()*rea[fumal].v( y[nfum]);     dydx[nfum] -= flx[fumal];  dydx[nmal] += flx[fumal];  
flx[malfum]= rea[D].v()*rea[malfum].v( y[nmal]);   dydx[nmal] -= flx[malfum];  dydx[nfum] += flx[malfum];  
flx[maloa]= rea[maloa].v( y[nmal], y[n_nad]);      dydx[nmal] -= flx[maloa];  dydx[n_nad] -= flx[maloa];  dydx[noa] += flx[maloa];  
flx[oamal]= rea[oamal].v( y[noa], h_nad);          dydx[noa] -= flx[oamal];  dydx[n_nad] += flx[oamal];  dydx[nmal] += flx[oamal];  
flx[pc]= rea[pc].v( y[npyrm], y[n_atp]);           dydx[npyrm] -= flx[pc];  dydx[n_atp] -= flx[pc];  dydx[noa] += flx[pc];  
flx[malicm]= rea[malicm].v( y[nmal], y[n_nad]);    dydx[nmal] -= flx[malicm];  dydx[n_nad] -= flx[malicm];  dydx[npyrm] += flx[malicm];  
flx[malicc]= rea[malicc].v( y[noac], y[n_nad]);    dydx[noac] -= flx[malicc];  dydx[n_nad] -= flx[malicc];  dydx[npyr] += flx[malicc];  
flx[ppp]= rea[ppp].v( y[nh6]);                     dydx[nh6] -= flx[ppp];  dydx[np5] += flx[ppp];  
flx[oacd]= rea[oacd].v( y[noac]);                  dydx[noac] -= flx[oacd];  dydx[nmal] += flx[oacd];  
flx[mald]= rea[mald].v( y[nmal]);                  dydx[nmal] -= flx[mald];  dydx[noac] += flx[mald];  
flx[citdmc]= rea[citdmc].v( y[ncit]);              dydx[ncit] -= flx[citdmc];  dydx[ncitc] += flx[citdmc];  
flx[citdcm]= rea[citdcm].v( y[ncitc]);             dydx[ncitc] -= flx[citdcm];  dydx[ncit] += flx[citdcm];  
flx[akgdmc]= rea[akgdmc].v( y[nakg]);              dydx[nakg] -= flx[akgdmc];  dydx[nakgc] += flx[akgdmc];  
flx[akgdcm]= rea[akgdcm].v( y[nakgc]);             dydx[nakgc] -= flx[akgdcm];  dydx[nakg] += flx[akgdcm];  
flx[coaout]= rea[coaout].v( y[ncitc], y[n_atp]);   dydx[ncitc] -= flx[coaout];  dydx[n_atp] -= flx[coaout];  dydx[noac] += flx[coaout];  
flx[citakg1]= rea[citakg1].v( y[ncitc], y[n_nad]); dydx[ncitc] -= flx[citakg1];  dydx[n_nad] -= flx[citakg1];  dydx[nakgc] += flx[citakg1];  
flx[akgcit1]= rea[akgcit1].v( y[nakgc], h_nad);    dydx[nakgc] -= flx[akgcit1];  dydx[n_nad] += flx[akgcit1];  dydx[ncitc] += flx[akgcit1];  
flx[gln_in]= rea[gln_in].v( xx[ngln]);             dydx[nakgc] += flx[gln_in];  
flx[gln_out]= rea[gln_out].v( y[nakgc]);           dydx[nakgc] -= flx[gln_out];  
flx[gluin]= rea[gluin].v( xx[nglu]);               dydx[nakgc] += flx[gluin];  
flx[gluout]= rea[gluout].v( y[nakgc]);             dydx[nakgc] -= flx[gluout];  
flx[t3ser]= rea[t3ser].v( y[nt3]);                 dydx[nt3] -= flx[t3ser];  
flx[serpyr]= rea[serpyr].v( xx[nser]);             dydx[npyrm] += flx[serpyr];  
flx[asp_o]= rea[asp_o].v( y[noac]);                dydx[noac] -= flx[asp_o];  
flx[asp_i]= rea[asp_i].v();                        dydx[noac] += flx[asp_i];  
flx[ala_o]= rea[ala_o].v( y[npyr]);                dydx[npyr] -= flx[ala_o];  
flx[ala_i]= rea[ala_i].v();                        dydx[npyr] += flx[ala_i];  
flx[r5_o]= rea[r5_o].v( y[np5]);                   dydx[np5] -= flx[r5_o];  
flx[r5_i]= rea[r5_i].v();                          dydx[np5] += flx[r5_i];  
flx[cystin]= rea[cystin].v();                      dydx[npyrm] += flx[cystin];  
flx[proin]= rea[proin].v();                        dydx[nakgc] += flx[proin];  
flx[proout]= rea[proout].v( y[nakgc]);             dydx[nakgc] -= flx[proout];  
flx[kgin]= rea[kgin].v();                          dydx[nakg] += flx[kgin];  
flx[coain]= rea[coain].v();                        dydx[ncoa] += flx[coain];  
flx[gln_pr]= rea[gln_pr].v( xx[ngln]);             
flx[ser_pr]= rea[ser_pr].v( xx[nser]);             
flx[asp_pr]= rea[asp_pr].v( xx[nasp]);             
flx[ala_pr]= rea[ala_pr].v( xx[nala]);             
flx[pro_pr]= rea[pro_pr].v( xx[npro]);             
flx[trpala]= rea[trpala].v();                      
flx[mthf]= rea[mthf].v( y[ncthf]);                 dydx[ncthf] -= flx[mthf];  
flx[thf]= rea[thf].v( y[ncthf]);                   dydx[ncthf] -= flx[thf];  
flx[sergly]= rea[sergly].v( xx[nser], xthf);       
flx[glyser]= rea[glyser].v( xx[ngly], y[ncthf]);   dydx[ncthf] -= flx[glyser];  
flx[cs0]= rea[D].v()*rea[cs0].v( y[noa], y[ncoa]); dydx[noa] -= flx[cs0];  dydx[ncoa] -= flx[cs0];  dydx[ncit] += flx[cs0];  
flx[D]= rea[D].v();                                
flx[atpase]= rea[atpase].v( y[n_atp]);             dydx[n_atp] -= flx[atpase];  
flx[resp]= rea[resp].v( h_nad, a_dp);              dydx[n_nad] += flx[resp];  dydx[n_atp] += flx[resp];  dydx[n_atp] += flx[resp];  dydx[n_atp] += flx[resp];  
aldolase.st1fl(&flx[aldfl], y[nfbp], y[nt3]);      dydx[nfbp] -= flx[aldfl];    dydx[nt3] += 2.*flx[aldfl];
ta.st1fl(&flx[tafl], y[nh6]/fh6, y[nt3]/ft3, y[ne4], y[ns7]);
				dydx[nh6] -= flx[tafl];	 dydx[nt3] += flx[tafl];
				dydx[ne4] -= flx[tafl];	 dydx[ns7] += flx[tafl];
tk.st1fl(&flx[tkfl], y[nt3]/ft3, y[np5], y[ne4], y[nh6]/fh6, y[np5], y[ns7]);
				dydx[np5] -= flx[tkfl];	dydx[nt3] += flx[tkfl];
				dydx[ns7] -= flx[tkfl+1];	dydx[np5] += flx[tkfl+1];
				dydx[nh6] -= flx[tkfl+2];	dydx[ne4] += flx[tkfl+2];
for(int i=0;i<numx;i++) dydx[i]*=(dt/Vi);
}

void Fit::ff(const double *y,double *dydx) {
	dydx[ngl] = (- flx[hk])*dt;
	dydx[nlac] = (+flx[pyrlac]- flx[lacpyr])*dt;
	dydx[nglu] = (+flx[gluout]- flx[gluin])*dt;
	dydx[ngln] = (+flx[gln_out]- flx[gln_in]- flx[gln_pr])*dt;
	dydx[nala] = (+flx[ala_o]- flx[ala_i]+flx[trpala]- flx[ala_pr])*dt;
	dydx[nasp] = (+flx[asp_o]- flx[asp_i]- flx[asp_pr])*dt;
	dydx[nser] = (+flx[t3ser]- flx[serpyr]- flx[ser_pr]+flx[glyser]- flx[sergly])*dt;
	dydx[ngly] = (+flx[sergly]- flx[glyser])*dt;
	dydx[npro] = (+flx[proout]- flx[proin]- flx[pro_pr])*dt;
	dydx[nrna] = (+flx[r5_o]- flx[r5_i])*dt;
}

void Parray::init(){ft3=10.; fh6=7.; 
	tk.setk(rea[rtk].getpar());
	ta.setk(rea[rta].getpar());
	aldolase.setk(rea[rald].getpar());}
void Parray::fin(double y[]){
	aldolase.st1fl(&flx[aldfl], y[nfbp], y[nt3]);
	tk.st1fl(&flx[tkfl], y[nt3]/ft3, y[np5], y[ne4], y[nh6]/fh6, y[np5], y[ns7]);
	ta.st1fl(&flx[tafl], y[nh6]/fh6, y[nt3]/ft3, y[ne4], y[ns7]);
	aldolase.st2fl(&flx[aldfl], y[nfbp], y[nt3]);
	tk.st2fl(&flx[tkfl], y[nt3]/ft3, y[np5], y[ne4], y[nh6]/fh6, y[np5], y[ns7]);
	ta.st2fl(&flx[tafl], y[nh6]/fh6, y[nt3]/ft3, y[ne4], y[ns7]);
flfor(y);
}
void Parray::flfor(double *y){
for(int i=0;i<nflx;i++) fluxes[i] = flx[i] * dt/Vi;
fluxes[pfk] /= y[nh6];
fluxes[fbpase] /= y[nfbp];
fluxes[t3pep] /= y[nt3];
fluxes[pept3] /= y[npep];
fluxes[pk] /= y[npep];
fluxes[pyrlac] /= y[npyr];
fluxes[lacpyr] /= xx[nlac];
fluxes[pyrdcm] /= y[npyr];
fluxes[pyrdmc] /= y[npyrm];
fluxes[pdh] /= y[npyrm];
fluxes[citakg] /= y[ncit];
fluxes[akgfum] /= y[nakg];
fluxes[fumal] /= y[nfum];
fluxes[malfum] /= y[nmal];
fluxes[maloa] /= y[nmal];
fluxes[oamal] /= y[noa];
fluxes[pc] /= y[npyrm];
fluxes[malicm] /= y[nmal];
fluxes[malicc] /= y[noac];
fluxes[ppp] /= y[nh6];
fluxes[oacd] /= y[noac];
fluxes[mald] /= y[nmal];
fluxes[citdmc] /= y[ncit];
fluxes[citdcm] /= y[ncitc];
fluxes[akgdmc] /= y[nakg];
fluxes[akgdcm] /= y[nakgc];
fluxes[coaout] /= y[ncitc];
fluxes[citakg1] /= y[ncitc];
fluxes[akgcit1] /= y[nakgc];
fluxes[gln_in] /= xx[ngln];
fluxes[gln_out] /= y[nakgc];
fluxes[gluin] /= xx[nglu];
fluxes[gluout] /= y[nakgc];
fluxes[t3ser] /= y[nt3];
fluxes[serpyr] /= xx[nser];
fluxes[asp_o] /= y[noac];
fluxes[ala_o] /= y[npyr];
fluxes[r5_o] /= y[np5];
fluxes[proout] /= y[nakgc];
fluxes[gln_pr] /= xx[ngln];
fluxes[ser_pr] /= xx[nser];
fluxes[asp_pr] /= xx[nasp];
fluxes[ala_pr] /= xx[nala];
fluxes[pro_pr] /= xx[npro];
fluxes[mthf] /= y[ncthf];
fluxes[thf] /= y[ncthf];
fluxes[sergly] /= xx[nser]*xthf;
fluxes[glyser] /= xx[ngly]*y[ncthf];
fluxes[cs0] /= y[noa]*y[ncoa];
fluxes[rald] /= xx[nfbp];
fluxes[rald+1] /= (xx[nt3]*xx[nt3]);
fluxes[rald+2] /= xx[nfbp];
}

