#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
#include "solvers.h"
#include "analis.h"
using namespace std;
const int hk=0, pfk=hk+1, fbpase=pfk+1, t3pep=fbpase+1, pept3=t3pep+1, pk=pept3+1, pyrlac=pk+1, lacpyr=pyrlac+1, pyrdcm=lacpyr+1, pyrdmc=pyrdcm+1, pdh=pyrdmc+1, citakg=pdh+1, akgfum=citakg+1, fumal=akgfum+1, malfum=fumal+1, fumout=malfum+1, maloa=fumout+1, oamal=maloa+1, pc=oamal+1, malicm=pc+1, malicc=malicm+1, pepck=malicc+1, ppp=pepck+1, oacd=ppp+1, mald=oacd+1, citdmc=mald+1, citdcm=citdmc+1, akgdmc=citdcm+1, akgdcm=akgdmc+1, coaout=akgdcm+1, coar=coaout+1, citakg1=coar+1, akgcit1=citakg1+1, gln_in=akgcit1+1, gln_out=gln_in+1, gln_pr=gln_out+1, gluin=gln_pr+1, gluout=gluin+1, t3ser=gluout+1, serpyr=t3ser+1, ser_pr=serpyr+1, sergly=ser_pr+1, glyser=sergly+1, thf=glyser+1, mthf=thf+1, asp_o=mthf+1, asp_i=asp_o+1, asp_pr=asp_i+1, ala_o=asp_pr+1, ala_i=ala_o+1, trpala=ala_i+1, ala_pr=trpala+1, r5_o=ala_pr+1, r5_i=r5_o+1, cystin=r5_i+1, proin=cystin+1, proout=proin+1, pro_pr=proout+1, kgin=pro_pr+1, coain=kgin+1, aglin=coain+1, aglout=aglin+1, glycogin=aglout+1, glycogout=glycogin+1, D=glycogout+1, cs0=D+1, resp=cs0+1, atpase=resp+1, rald=atpase+1, rta=rald+1, rtk=rta+1, nrea=rtk+1, 
aldf=atpase+1, aldrev=aldf+1, aldfli=aldrev+1, aldi1=aldfli+1, tafl=aldi1+1, s7f6a=tafl+1, f6g3a=s7f6a+1, s7e4a=f6g3a+1, tkfl=s7e4a+1, s7p5=tkfl+1, f6p5=s7p5+1, p5f6=f6p5+1, f6s7=p5f6+1, s7f6=f6s7+1, p5g3i=s7f6+1, f6e4i=p5g3i+1, s7p5i=f6e4i+1, nflx =s7p5i+1;

const int nh6=0, nfbp=nh6+1, nt3=nfbp+1, npep=nt3+1, npyr=npep+1, npyrm=npyr+1, ncoa=npyrm+1, noa=ncoa+1, noac=noa+1, ncit=noac+1, ncitc=ncit+1, nakg=ncitc+1, nakgc=nakg+1, nfum=nakgc+1, nmal=nfum+1, np5=nmal+1, ne4=np5+1, ns7=ne4+1, natp=ns7+1, nnad=natp+1, ngl=nnad+1, nlac=ngl+1, nglu=nlac+1, ngln=nglu+1, nala=ngln+1, nasp=nala+1, nser=nasp+1, npro=nser+1, nrna=npro+1, ngly=nrna+1, ncthf=ngly+1, ncoac=ncthf+1, nagl=ncoac+1, nglycog=nagl+1, nmet=nglycog+1, numx=nnad+1;

   Fit Problem;
   const double thft(1.);
   double dt, xx[nmet],flx[nflx],fluxes[nflx];
   double xm0,xinit1[nmet],xinit2[nmet];
   int Parray::par[nrea];
   string Parray::fid[nflx],Parray::fname[nflx],Parray::fschem[nflx], Parray::namex[numx];
   Reapar Parray::rea[nrea];
   double Analis::nv1[nrea], Analis::nv2[nrea];
void Parray::rnames(ifstream& fi){
   for (int i=0;i<nflx;i++)
        fi>>fid[i]>>fid[i]>>fname[i]>>fschem[i];
        }
   
void Fit::ff(const double *y,double *dydx) { for(int i=0;i<numx;i++) dydx[i]=0.;
dydx[ngl] = - flx[hk];                  //glucose
dydx[nlac] = flx[pyrlac]-flx[lacpyr];   //lactate
dydx[nglu] = flx[gluout]-flx[gluin];    //glutamate
dydx[ngln] = flx[gln_out]-flx[gln_in]-flx[gln_pr];//glutamine
dydx[nala] = flx[ala_o]-flx[ala_i]+flx[trpala]-flx[ala_pr];//alanine
dydx[nasp] = flx[asp_o]-flx[asp_i]-flx[asp_pr];//aspartate
dydx[nser] = flx[t3ser]-flx[serpyr]-flx[ser_pr]+flx[glyser]-flx[sergly];//serine
dydx[ngly] = flx[sergly]-flx[glyser];   //glycine
dydx[ncthf] = flx[sergly]-flx[glyser]+flx[mthf]-flx[thf];//tetrahydrofolate
dydx[npro] = flx[proout]-flx[proin]-flx[pro_pr];//proline
dydx[nrna] = flx[r5_o]-flx[r5_i];       //ribose
dydx[ncoac] = flx[coaout] - flx[coar];  //accoa
dydx[nagl] = flx[aglin] - flx[aglout];  //a-glycerate
dydx[nglycog] = flx[glycogin] - flx[glycogout];//glycogen
}
void Fit::f(const double *y,double *dydx) { for(int i=0;i<numx;i++) dydx[i]=0.;
                                            for(int i=0;i<nflx;i++) flx[i]=0.;
 double amp = -(sqrt(4.*xx[natp]*tan-3.*xx[natp]*xx[natp])-2.*tan+xx[natp])/2.;
   double adp = (sqrt(xx[natp])*sqrt(4.*tan-3.*xx[natp])-xx[natp])/2.;
     double nadh = tnad-xx[nnad];
      xthf=thft-xx[ncthf];
 flx[hk]= rea[hk].v();//*y[natp]/(nv[khk]+y[natp]+y[nh6]/nv[kihk]);
   dydx[natp]-= flx[hk]; dydx[nh6] += flx[hk];
     double yh6=pow(y[nh6],1.5);
 flx[pfk]= rea[pfk].v(yh6,amp);
 flx[fbpase]= rea[fbpase].ving(y[nfbp],amp,y[nh6]);
                                  dydx[nh6] -= flx[pfk];    dydx[nfbp] += flx[pfk];
                                  dydx[natp] -= flx[pfk]; 
                                  dydx[nfbp] -= flx[fbpase]; dydx[nh6] += flx[fbpase];
 flx[glycogin]=rea[glycogin].v(y[nh6]); dydx[nh6] += flx[glycogin];
 flx[glycogout]=rea[glycogout].v();     dydx[nh6] -= flx[glycogout];
aldolase.st1fl(&flx[aldf], y[nfbp], y[nt3]); dydx[nfbp] -= flx[aldf]; dydx[nt3] += 2.*flx[aldf];
 flx[t3pep]= rea[t3pep].v(y[nt3],adp,y[nnad]); 
 flx[pept3]= rea[pept3].v(y[npep],y[natp],nadh);
                                  dydx[nt3] -= flx[t3pep];       dydx[npep] += flx[t3pep];
                                  dydx[nnad]-=  flx[t3pep];      dydx[natp] += flx[t3pep];
                                  dydx[npep] -= flx[pept3];      dydx[nt3] += flx[pept3];
 flx[pk]=rea[pk].v(y[npep],adp);  
                                  dydx[npep] -= flx[pk];         dydx[npyr] += flx[pk];
                                                                 dydx[natp] += flx[pk];
 flx[pyrlac]= rea[pyrlac].v(y[npyr],nadh);
 flx[lacpyr]= rea[lacpyr].v(xx[nlac],y[nnad]);
                              dydx[npyr] -= flx[pyrlac];     dydx[nnad] += flx[pyrlac];
                              dydx[nnad] -= flx[lacpyr];     dydx[npyr] += flx[lacpyr];
flx[atpase]=rea[atpase].v(y[natp]);
                                          dydx[natp] += (flx[resp]*3.+flx[akgfum]*2.-flx[atpase]);
flx[pyrdcm]= rea[pyrdcm].v(y[npyr]); dydx[npyr] -= flx[pyrdcm]; dydx[npyrm] += flx[pyrdcm];
 flx[pyrdmc]= rea[pyrdmc].v(y[npyrm]); dydx[npyrm] -= flx[pyrdmc]; dydx[npyr] += flx[pyrdmc];
flx[ppp]= rea[ppp].v(y[nh6]);            dydx[nh6] -= flx[ppp];     dydx[np5] += flx[ppp];
 flx[pdh]= rea[pdh].v(y[npyrm],y[nnad]);
                              dydx[npyrm] -= flx[pdh];       dydx[ncoa] += flx[pdh];
                              dydx[nnad] -=  flx[pdh];
 flx[D]= rea[D].v();
 flx[citakg]= rea[D].v()*rea[citakg].v(y[ncit],y[nnad]);
                              dydx[ncit] -= flx[citakg];    dydx[nakg] += flx[citakg];
                              dydx[nnad] -=  flx[citakg];
 flx[akgfum]= rea[D].v()*rea[akgfum].v(y[nakg],y[nnad],adp);
                                dydx[nakg] -= flx[akgfum];    dydx[nfum] += flx[akgfum];
                                dydx[nnad] -= flx[akgfum];    dydx[natp] += flx[akgfum];
 flx[fumal]= rea[D].v()*rea[fumal].v(y[nfum]);   dydx[nfum] -= flx[fumal]; dydx[nmal] += flx[fumal];
 flx[malfum]= rea[D].v()*rea[malfum].v(y[nmal]); dydx[nmal] -= flx[malfum]; dydx[nfum] += flx[malfum];
 flx[fumout]= rea[fumout].v(y[nfum]);            dydx[nfum] -= flx[fumout];
 flx[maloa]= rea[maloa].v(y[nmal], y[nnad]);
 flx[oamal]= rea[oamal].v(y[noa], nadh);
                                dydx[nmal] -= flx[maloa];      dydx[noa] += flx[maloa];
                                dydx[noa] -= flx[oamal];       dydx[nmal] += flx[oamal];
                                dydx[nnad]-=  flx[maloa];      dydx[nnad] += flx[oamal];
flx[cs0]= rea[D].v()*rea[cs0].v(y[ncoa], y[noa]);
				dydx[ncoa] -= flx[cs0];	dydx[noa] -= flx[cs0];	dydx[ncit] += flx[cs0];
 flx[pc]= rea[pc].v(y[npyrm], y[natp]);
                                dydx[npyrm] -= flx[pc];        dydx[noa] += flx[pc];
                                dydx[natp] -= flx[pc];
 flx[malicm]= rea[malicm].v(y[nmal], y[nnad]);
                                dydx[nmal] -= flx[malicm];     dydx[npyrm] += flx[malicm];
                                dydx[nnad]-=  flx[malicm];
 flx[pepck]= rea[pepck].v(y[noa], y[natp]);
                                dydx[noa] -= flx[pepck];       dydx[npep] += flx[pepck];
                                dydx[natp] -= flx[pepck];
 flx[malicc]= rea[malicc].v(y[noac], y[nnad]);
                                dydx[noac] -= flx[malicc]; dydx[npyr] += flx[malicc]; dydx[nnad] -= flx[malicc];
flx[oacd]= rea[oacd].v(y[noac]);           dydx[noac] -= flx[oacd];       dydx[nmal] += flx[oacd];
flx[mald]= rea[mald].v( y[nmal]);           dydx[nmal] -= flx[mald];       dydx[noac] += flx[mald];
flx[citdmc]= rea[citdmc].v(y[ncit]);     dydx[ncit] -= flx[citdmc];     dydx[ncitc] += flx[citdmc];
flx[citdcm]= rea[citdcm].v( y[ncitc]);    dydx[ncitc] -= flx[citdcm];    dydx[ncit] += flx[citdcm];
flx[akgdmc]= rea[akgdmc].v( y[nakg]);   dydx[nakg] -= flx[akgdmc];     dydx[nakgc] += flx[akgdmc];
flx[akgdcm]= rea[akgdcm].v( y[nakgc]);    dydx[nakgc] -= flx[akgdcm];    dydx[nakg] += flx[akgdcm];
flx[coaout]= rea[coaout].v(y[ncitc], y[natp]);
                                dydx[ncitc] -= flx[coaout];    dydx[noac] += flx[coaout];
                                dydx[natp] -= flx[coaout];
flx[coar]= rea[coar].v(xx[ncoac]);
flx[citakg1]= rea[citakg1].v(y[ncitc]); 
flx[akgcit1]= rea[akgcit1].v(y[nakgc], y[nnad]); 
                                dydx[ncitc] -= flx[citakg1];   dydx[nakgc] += flx[citakg1];
                                dydx[nakgc] -= flx[akgcit1];   dydx[ncitc] += flx[akgcit1];
                                dydx[nnad] -=  flx[citakg1];   dydx[nnad] += flx[akgcit1];
flx[gln_in]= rea[gln_in].vin(xx[ngln], y[nakgc]);     dydx[nakgc] += flx[gln_in];
flx[gln_out]= rea[gln_out].v( y[nakgc]); dydx[nakgc] -= flx[gln_out];   
flx[gln_pr]= rea[gln_pr].v(xx[ngln]);
flx[gluin]= rea[gluin].v(xx[nglu]);                        dydx[nakgc] += flx[gluin];
flx[gluout]= rea[gluout].v(y[nakgc]); dydx[nakgc] -= flx[gluout];    
flx[t3ser]= rea[t3ser].v( y[nt3]);  dydx[nt3] -= flx[t3ser];
flx[serpyr]= rea[serpyr].v(xx[nser]);                  dydx[npyrm] += flx[serpyr];
flx[ser_pr]= rea[ser_pr].v(xx[nser]);
flx[sergly]= rea[sergly].v(xx[nser],xthf);
flx[glyser]= rea[glyser].v(xx[ngly],xx[ncthf]); 
flx[mthf]= rea[mthf].v(xthf);
flx[thf]= rea[thf].v(xx[ncthf]);
flx[asp_o]= rea[asp_o].v( y[noac]);  dydx[noac] -= flx[asp_o];      
flx[asp_i]= rea[asp_i].v();                  dydx[nfum] += flx[asp_i];
flx[asp_pr]= rea[asp_pr].v(xx[nasp]);/**/
flx[ala_o]= rea[ala_o].v( y[nt3]); dydx[nt3] -= flx[ala_o];      
flx[ala_i]= rea[ala_i].v(xx[nala]);                       dydx[npyr] += flx[ala_i];
flx[trpala]= rea[trpala].v();
flx[ala_pr]= rea[ala_pr].v(xx[nala]);
flx[r5_o]= rea[r5_o].v( y[np5]);   dydx[np5] -= flx[r5_o];        
flx[r5_i]= rea[r5_i].v(y[nrna]);                             dydx[np5] += flx[r5_i];
flx[cystin]= rea[cystin].v();                                       dydx[npyr] += flx[cystin];
flx[proin]= rea[proin].v(xx[npro]);                   dydx[nakgc] += flx[proin];
flx[proout]=rea[proout].v(y[nakgc]);    dydx[nakgc] -= flx[proout];    
flx[pro_pr]= rea[pro_pr].v(xx[npro]);
flx[kgin]= rea[kgin].v();                                  dydx[nakg] += flx[kgin];/**/
flx[coain]= rea[coain].v();                                       dydx[ncoa] += flx[coain]; 
flx[aglin]= rea[aglin].v(xx[nagl]);     dydx[nt3] +=flx[aglin];
flx[aglout]= rea[aglout].v(y[nt3]);        dydx[nt3] -=flx[aglout];
flx[resp]=rea[resp].v(nadh,adp);
                                          dydx[nnad] += flx[resp];
ta.st1fl(&flx[tafl], y[nh6]/fh6, y[nt3]/ft3, y[ne4], y[ns7]);
				dydx[nh6] -= flx[tafl];		 dydx[nt3] += flx[tafl];
				dydx[ne4] -= flx[tafl];		 dydx[ns7] += flx[tafl];
tk.st1fl(&flx[tkfl], y[nt3]/ft3, y[np5], y[ne4], y[nh6]/fh6, y[np5], y[ns7]);
				dydx[np5] -= flx[tkfl];		dydx[nt3] += flx[tkfl];
				dydx[ns7] -= flx[tkfl+1];	dydx[np5] += flx[tkfl+1];
				dydx[nh6] -= flx[tkfl+2];	dydx[ne4] += flx[tkfl+2];
//					dydx[nnad]=0.;dydx[natp]=0.;
}

void Parray::init(){ft3=10.; fh6=7.; 
	tk.setk(rea[rtk].getpar() );
	ta.setk(rea[rta].getpar());
	aldolase.setk(rea[rald].getpar());}
void Parray::fin(double y[]){
//cout<<"serfl="<<flx[serpyr]<<"; vm="<<nv[Vserpyr]<<"; ser0="<<y[nser]<<"; "<<Problem.MM(nv[Vserpyr],0.02,y[nser])<<endl;
	aldolase.st1fl(&flx[aldf], y[nfbp], y[nt3]);
	tk.st1fl(&flx[tkfl], y[nt3]/ft3, y[np5], y[ne4], y[nh6]/fh6, y[np5], y[ns7]);
	ta.st1fl(&flx[tafl], y[nh6]/fh6, y[nt3]/ft3, y[ne4], y[ns7]);
	aldolase.st2fl(&flx[aldf], y[nfbp], y[nt3]);
	tk.st2fl(&flx[tkfl], y[nt3]/ft3, y[np5], y[ne4], y[nh6]/fh6, y[np5], y[ns7]);
	ta.st2fl(&flx[tafl], y[nh6]/fh6, y[nt3]/ft3, y[ne4], y[ns7]);
flfor(y);
}
void Parray::flfor(double y[]){
for(int i=0;i<nflx;i++) fluxes[i] = flx[i] * dt/Vi;
fluxes[pfk] /= y[nh6];
fluxes[fbpase] /= y[nfbp];
fluxes[t3pep] /= y[nt3];
fluxes[pept3] /= y[npep];
fluxes[pk] /= y[npep];
fluxes[pyrlac] /= y[npyr];
fluxes[pyrdcm] /= y[npyr];
fluxes[pyrdmc] /= y[npyrm];
fluxes[pdh] /= y[npyrm];
fluxes[citakg] /= y[ncit];
fluxes[akgfum] /= y[nakg];
fluxes[fumal] /= y[nfum];
fluxes[malfum] /= y[nmal];
fluxes[fumout] /= y[nfum];
fluxes[maloa] /= y[nmal];
fluxes[oamal] /= y[noa];
fluxes[pc] /= y[npyrm]; 
fluxes[malicm] /= y[nmal];
fluxes[malicc] /= y[noac];
fluxes[pepck] /= y[noa];
fluxes[ppp] /= y[nh6];
fluxes[oacd] /= y[noac];
fluxes[mald] /= y[nmal];
fluxes[citdmc] /= y[ncit];
fluxes[citdcm] /= y[ncitc];
fluxes[akgdmc] /= y[nakg];
fluxes[akgdcm] /= y[nakgc];
fluxes[coaout] /= y[ncitc];
fluxes[coar] /= xx[ncoac];
fluxes[citakg1] /= y[ncitc];
fluxes[akgcit1] /= y[nakgc];
fluxes[gln_out] /= y[nakgc];
fluxes[gln_in]  /= xx[ngln];
fluxes[gln_pr]  /= xx[ngln];
fluxes[gluin]  /= xx[nglu];
fluxes[gluout] /= y[nakgc];
fluxes[t3ser] /= y[nt3];
fluxes[serpyr]  /= xx[nser];
fluxes[ser_pr]  /= xx[nser];
fluxes[sergly]  /= (xthf*xx[nser]);
fluxes[glyser]  /= (xx[ncthf]*xx[ngly]);
fluxes[mthf] /= xthf;
fluxes[thf] /= xx[ncthf];
fluxes[asp_o] /= y[noac];
fluxes[asp_pr]  /= xx[nasp];
fluxes[ala_o] /= y[nt3];
fluxes[ala_i]  /= xx[nala];
fluxes[ala_pr]  /= xx[nala];
fluxes[r5_o] /= y[np5];
fluxes[proout] /= y[nakgc];
fluxes[proin]  /= xx[npro];
fluxes[pro_pr]  /= xx[npro];
fluxes[aglin]  /= xx[nagl];
fluxes[aglout]  /= y[nt3];
 fluxes[glycogin] /= y[nh6];
 fluxes[glycogout] /= xx[nglycog];
fluxes[cs0] /= (y[ncoa]*y[noa]);
fluxes[aldf] /= y[nfbp];
fluxes[aldf+1] /= (y[nt3]*y[nt3]);
fluxes[aldf+2] /= y[nfbp];
}

void Fit::cont(const int chpar,const double pmin,const double pint){
// get a set of solutions increasing or decreasing a parameter in an interval
 const int nss(50);
   double dif,fact(1.0),f1(1.02),mx(0.),ss[2][nss],df[33], pdf[33], a=rea[chpar].v();
   const double dp=pint/(double)nss;
   ostringstream fn; ofstream fo("00000"); int j(0);
   for(int k=0;k<2;k++){
    for(int ii=0;ii<=nss;ii++){
     if(k) rea[chpar].setVm(pmin+(pint-ii*dp));
      else rea[chpar].setVm(pmin+ii*dp);
  try{ tsolve(18000.); 
	} catch( char const* str ){cout << "exception: "<< str <<endl;}
   cout<<rea[chpar].getpar()<<": atp="<<xx[natp]<<"; nad="<<xx[nnad]<<"; pyr="<<xx[npyr]<<endl; 
   fo<<rea[chpar].getpar()<<" "<<xx[natp]<<" "<<xx[nnad]<<" "<<xx[npyr]<<endl; 
      } }
   rea[chpar].setVm(a); }

