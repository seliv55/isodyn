//---------------------------------------------------------------------------
//#include <fstream>
#include <iostream>
#pragma hdrstop
//#include "nr.h"
#include "nums.hh"
#include "tk.hh"

//---------------------------------------------------------------------------
using namespace std;
//#pragma package(smart_init)
//         enum {nh6=1,nt3,npep=5,npyr,ncoa,noaa,ncit,nglgn,np5=14,ne4,ns7,nmal};
void TK::setk( double *pv) {
     const double kdX5 = 0.055*pv[1],kdG3 = 0.0344*pv[2],kdR5 = 0.27*pv[3],
           kdF6 = 0.9*pv[4], diss = 5000, kefg = 29.7, kexr = 0.48;
           e0=0.0025*pv[0];		//e0tk
		k0[1] = diss;	k[1]=k0[1]/kdX5;
		k[2] = 5300.;	k0[2]=5300.;
		k[3] = diss;	k0[3]=k[3]/kdG3;
		k0[4] = diss;	k[4] = k0[4]/kdR5;
		k[5] = 5300.;	k0[5]= 1770.;
		k[6] = diss;	k0[6] = k[6]*k[2]*k[5]*kdG3*kexr/k0[2]/k0[5]/kdR5/kdX5;
		k0[7] = diss;	k[7] = k0[7]/kdF6;
		k[8] = 5300.*pv[6];	k0[8] = 5300.*pv[5];
		k[9] = diss;	k0[9] = k[9]*k[8]*k0[2]*kdX5*kefg/k0[8]/k[2]/kdF6/kdG3;
const double &k1=k[1],&k2=k[2],&k3=k[3],&k4=k[4],&k5=k[5],&k6=k[6],&k7=k[7],&k8=k[8],&k9=k[9],
&k01=k0[1],&k02=k0[2],&k03=k0[3],&k04=k0[4],&k05=k0[5],&k06=k0[6],&k07=k0[7],&k08=k0[8],&k09=k0[9];
// puts together the constants, making easier the flux calculations during the ODE solving
const double tmp1= (k8*k9 + k07*(k08 + k9)), tmp2= (k5*k6 + k04*(k05 + k6)), tmp3= (k2*k3 + k01*(k02 + k3)),
         tmp4= (k01 + k2), tmp5= (k02 + k3), tmp6= (k05 + k6), tmp7= (k04 + k5), tmp8= (k08 + k9), tmp9= (k07 + k8);
      		kDe[0].a[0] = k01*k02*k03*tmp2*tmp1;
		kDe[0].a[1] = k07*k08*k09*tmp3*tmp2;
		kDe[0].a[2] = tmp3*k4*k5*k6*tmp1;
//De=g3*kDe[0] + ce4*kDe[1] + kDe[2]*r5i;
	kDe[1].a[0] = k02*k03*k7*k8*k9*tmp2;
	kDe[1].a[1] = k07*k08*k09*k1*tmp5*tmp2;
	kDe[1].a[2] = k1*tmp5*k4*k5*k6*tmp1;
	kDe[1].a[3] = k02*k03*k04*k05*k06*tmp1;
	kDe[1].a[4] = k02*k03*k1*tmp2*tmp1;
//Ded1=f6i*g3i*kDed1[0] + ce4*x5*kDed1[1] + kDed1[2]*r5i*x5 + g3*cs7*kDed1[3] + g3*x5*kDed1[4]);
   		kDe[2].a[0] = tmp2*k07*k08*k09*k1*k2;
		kDe[2].a[1] = tmp2*k03*tmp4*k7*k8*k9;
		kDe[2].a[2] = k1*k2*k4*k5*k6*tmp1;
		kDe[2].a[3] = k03*k04*k05*k06*tmp4*tmp1;
		kDe[2].a[4] = k03*k1*k2*tmp2*tmp1;
//Dega1=(ce4*x5*a0 + f6i*g3*a1) + (r5i*x5*a2 + g3*(cs7*a3 + x5*a4));;
	kDe[3].a[0] = tmp3*tmp2*k7*k8*k9;
	kDe[3].a[1] = k04*k05*k06*tmp3*tmp1;
	kDe[3].a[2] = k1*k2*k3*tmp2*tmp1;
//Deg=f6i*kDeg[0] + cs7*kDeg[1] + x5*kDeg[2];
		kDe[4].a[0] = tmp3*k05*k06*k07*k08*k09;
		kDe[4].a[1] = tmp3*k4*tmp6*k7*k8*k9;
		kDe[4].a[2] = k01*k02*k03*k05*k06*tmp1;
		kDe[4].a[3] = k4*k05*k06*tmp3*tmp1;
		kDe[4].a[4] = k4*k1*k2*k3*tmp6*tmp1;
//Dega2=(ce4*cs7*a0 + f6i*r5i*a1) + (g3*cs7*a2 + r5i*(cs7*a3 + x5*a4));
	kDe[5].a[0] = tmp3*k06*k07*k08*k09*tmp7;
	kDe[5].a[1] = tmp3*k4*k5*k7*k8*k9;
	kDe[5].a[2] = k1*k2*k3*k4*k5*tmp1;
	kDe[5].a[3] = k06*tmp3*k4*k5*tmp1;
	kDe[5].a[4] = k01*k02*k03*k06*tmp7*tmp1;
//Ded2=(ce4*cs7*a0 + f6i*r5i*a1) + (r5i*(x5*a2 + cs7*a3) + g3*cs7*a4);
       	kDe[6].a[0] = tmp2*k08*k09*k1*k2*k3;
	kDe[6].a[1] = tmp2*k01*k02*k03*k7*tmp8;
	kDe[6].a[2] = k08*k09*tmp3*k04*k05*k06;
	kDe[6].a[3] = k08*k09*tmp3*tmp2*k7;
	kDe[6].a[4] = tmp3*k4*k5*k6*k7*tmp8;
//Ded3=(ae4*ax5*a0 + af6*ag3*a1) + (ae4*(as7*a2 + af6*a3) + af6*ar5*a4);
		kDe[7].a[0] = tmp2*k01*k02*k03*k7*k8;
		kDe[7].a[1] = tmp2*k09*k1*k2*k3*tmp9;
		kDe[7].a[2] = tmp3*k4*k5*k6*k7*k8;
		kDe[7].a[3] = tmp3*k09*tmp2*k7*k8;
		kDe[7].a[4] = tmp3*k09*k04*k05*k06*tmp9;
//Dega3=(af6*ag3* + ae4*ax5*) + (af6*ar5* + ae4*(af6* + as7*));
}
void TK::Graf(double* const Dsp, double& g3i,double& r5i,double& f6i,double& x5i,double& ce4,double& cs7) const{
//defines the trees for enzyme species
Dsp[0]=g3i*kDe[0].a[0]+ ce4*kDe[0].a[1] + kDe[0].a[2]*r5i;//e
Dsp[1]=f6i*g3i*kDe[1].a[0] + (ce4*kDe[1].a[1] + kDe[1].a[2]*r5i)*x5i + g3i*(cs7*kDe[1].a[3] + x5i*kDe[1].a[4]);//ex
Dsp[2]=ce4*x5i*kDe[2].a[0] + f6i*g3i*kDe[2].a[1] + r5i*x5i*kDe[2].a[2] + g3i*cs7*kDe[2].a[3] + g3i*x5i*kDe[2].a[4];//egg
Dsp[3]=f6i*kDe[3].a[0] + cs7*kDe[3].a[1] + x5i*kDe[3].a[2];//eg
Dsp[4]=ce4*cs7*kDe[4].a[0] + f6i*r5i*kDe[4].a[1] + g3i*cs7*kDe[4].a[2] + r5i*cs7*kDe[4].a[3] + x5i*r5i*kDe[4].a[4];//egr
Dsp[5]=ce4*cs7*kDe[5].a[0] + f6i*r5i*kDe[5].a[1] + x5i*r5i*kDe[5].a[2] + r5i*cs7*kDe[5].a[3] + kDe[5].a[4]*g3i*cs7;//es
Dsp[6]=Dsp[0];//e
Dsp[7]=ce4*x5i*kDe[6].a[0] + f6i*g3i*kDe[6].a[1] + kDe[6].a[2]*ce4*cs7 + ce4*f6i*kDe[6].a[3] + f6i*r5i*kDe[6].a[4];//ef
Dsp[8]=f6i*g3i*kDe[7].a[0] + kDe[7].a[1]*ce4*x5i + f6i*r5i*kDe[7].a[2] + ce4*f6i*kDe[7].a[3] + ce4*cs7*kDe[7].a[4];//ege
	Dsp[9]=Dsp[3];
	double sum=0; for(int i=1;i<9;i++) sum += Dsp[i];
	Dsp[10] = sum;
}

void TK::st1fl(double fl[], double& x5i,double& cs7,double& f6i,double& g3i,double& r5i,double& ce4) {
// calculates the elementary reaction rates and net fluxes
// for all the three TK reactions during the ODE solving.
	double Dsp[11];
	Graf(Dsp, g3i,r5i,f6i,x5i,ce4,cs7);
	for(int i=1;i<10;i++) {v[i]=Dsp[i-1]*e0*k[i]/Dsp[10];
		v0[i]=Dsp[i]*e0*k0[i]/Dsp[10];}
	fl[0] = v[2] - v0[2];	//flx13,x5p->g3p
	fl[1] = v0[5] - v[5];	//flx14,s7p->r5p
	fl[2] = v[8] - v0[8];	//flx15,f6p->e4p
}
void TK::st2fl( double* const fl, double& x5i,double& cs7,double& f6i,double& g3i,double& r5i,double& ce4) {
// calculates the isotope exchange fluxes for all the three TK reactions to be
//used for isotopomer computing.
	v[1]*=x5i; v0[3]*=g3i; v[4]*=r5i; v0[6]*=cs7; v[7]*=f6i; v0[9]*=ce4;
	double tmp[7],fract;
	tmp[0]=(v0[9] + v[8]);
	tmp[1]=(v0[6] + v[5]);
	tmp[2]=v0[8]*v0[9];
	tmp[3]=v0[5]*v0[6];
	tmp[4]=(v0[2]*v0[3] + v[1]*(v0[3] + v[2]));
	tmp[5]=(tmp[2] + v[7]*tmp[0]);
	tmp[6]=(tmp[3] + v[4]*tmp[1]);
	double denom=v0[2]*v0[3]*(v0[4]*tmp[3]*(tmp[2]+v[7]*tmp[0])+(tmp[3]+v[4]*tmp[1])*v[7]*v[8]*v[9]) + 
		v[1]*(v0[3]*(v0[4]*tmp[3]*(tmp[2] + v[7]*tmp[0]) + (tmp[3] + v[4]*tmp[1])*v[7]*v[8]*v[9]) + 
		v[2]*(v0[4]*tmp[3]*(tmp[2] + v[7]*tmp[0]) + (tmp[3] + v[4]*tmp[1])*(tmp[2]*v[3] + 
		v[7]*(tmp[0]*v[3] + v[8]*v[9]))));
fract = v[1]*v[2]*v[3]*v[4]*v[5]*tmp[5]/denom;		fl[0] = v[6]*fract;	//x5i->cs7,flx32
fract = v0[2]*v0[3]*v0[4]*v0[5]*v0[6]*tmp[5]/denom;	fl[1] = v0[1]*fract;	//cs7->x5i,flx33
fract = v0[2]*v0[3]*tmp[6]*v[7]*v[8]*v[9]/denom;	fl[2] = v0[1]*fract;	//f6i->x5i,flx34
fract = tmp[2]*v[1]*v[2]*v[3]*tmp[6]/denom;		fl[3] = v0[7]*fract;    //x5i->f6i,flx35
fract = tmp[4]*v[4]*v[5]*v[7]*v[8]*v[9]/denom;		fl[4] = v[6]*fract;	//f6i->cs7,flx36
fract = v0[4]*v0[5]*v0[6]*v0[8]*v0[9]*tmp[4]/denom;	fl[5] = v0[7]*fract;	//cs7->f6i,flx37
		denom = v0[2]*v0[3] + v0[3]*v[1] + v[1]*v[2];
fract = v[1]*v[2]/denom;      fl[6] = v[3]*fract-fl[0]-fl[3];	//x5i<->g3i,flx38
		denom=v0[8]*v0[9] + v0[9]*v[7] + v[7]*v[8];
fract = v[7]*v[8]/denom;      fl[7]=v[9]*fract-fl[2]-fl[4];	//f6i<->ce4,flx39
		denom=v0[5]*v0[6] + v0[6]*v[4] + v[4]*v[5];
fract = v[4]*v[5]/denom;      fl[8]=v[6]*fract-fl[0]-fl[4];	//r5i<->cs7,flx40
}
//      Doub3 TA::kDe[6];
//      double TA::k[7],TA::k0[7],TA::v[7],TA::v0[7];
//      double TA::e0;
void TA::setk(double *pv) {
     const double kefe = 0.37, diss = 5000., kdEA = 0.032*pv[1], kdFA = 0.72*pv[2],
           kdGA = 0.07*pv[3];
	   e0=0.0025*pv[0];	//e0ta
     k0[1] = diss; k[1]=k0[1]/kdFA;   k[2] = 5300.*pv[4];	k0[2]=5300.;
     k[3] = diss;  k0[3]=k[3]/kdGA;   k0[4] = diss;             k[4] = k0[4]/kdEA;
     k[5] = 5300.;     k0[5]= 1770.;  k[6] = diss;
     k0[6] = k[6]*k[2]*kdGA*k[5]*kefe/kdFA/k0[2]/kdEA/k0[5];
	kDe[0].a[0] = (k0[1]*k0[2]*k0[3]*k0[4]*k0[5] + k0[1]*k0[2]*k0[3]*k0[4]*k[6] + k0[1]*k0[2]*k0[3]*k[5]*k[6]);
	kDe[0].a[1] = (k0[1]*k0[2]*k[4]*k[5]*k[6] + k0[1]*k[3]*k[4]*k[5]*k[6] + k[2]*k[3]*k[4]*k[5]*k[6]);
	kDe[1].a[0] = (k0[2]*k0[3]*k0[4]*k0[5]*k0[6]);
	kDe[1].a[1] = (k0[2]*k0[3]*k0[4]*k0[5]*k[1] + k0[2]*k0[3]*k0[4]*k[1]*k[6] + k0[2]*k0[3]*k[1]*k[5]*k[6]);
	kDe[1].a[2] = (k0[2]*k[1]*k[4]*k[5]*k[6] + k[1]*k[3]*k[4]*k[5]*k[6]);
		kDe[2].a[0] = (k0[1]*k0[3]*k0[4]*k0[5]*k0[6] + k0[3]*k0[4]*k0[5]*k0[6]*k[2]);
		kDe[2].a[1] =(k0[3]*k0[4]*k0[5]*k[1]*k[2] + k0[3]*k0[4]*k[1]*k[2]*k[6] + k0[3]*k[1]*k[2]*k[5]*k[6]);
		kDe[2].a[2] =(k[1]*k[2]*k[4]*k[5]*k[6]);
	kDe[3].a[0] = (k0[1]*k0[2]*k0[4]*k0[5]*k0[6] + k0[1]*k0[4]*k0[5]*k0[6]*k[3] + k0[4]*k0[5]*k0[6]*k[2]*k[3]);
	kDe[3].a[1] = (k0[4]*k0[5]*k[1]*k[2]*k[3] + k0[4]*k[1]*k[2]*k[3]*k[6] + k[1]*k[2]*k[3]*k[5]*k[6]);
		kDe[4].a[0] = (k0[1]*k0[2]*k0[3]*k0[5]*k0[6]);
		kDe[4].a[1] = (k0[1]*k0[2]*k0[5]*k0[6]*k[4] + k0[1]*k0[5]*k0[6]*k[3]*k[4] + k0[5]*k0[6]*k[2]*k[3]*k[4]);
		kDe[4].a[2] = (k0[5]*k[1]*k[2]*k[3]*k[4] + k[1]*k[2]*k[3]*k[4]*k[6]);
	kDe[5].a[0] = (k0[1]*k0[2]*k0[3]*k0[4]*k0[6] + k0[1]*k0[2]*k0[3]*k0[6]*k[5]);
	kDe[5].a[1] = (k0[1]*k0[2]*k0[6]*k[4]*k[5] + k0[1]*k0[6]*k[3]*k[4]*k[5] + k0[6]*k[2]*k[3]*k[4]*k[5]);
	kDe[5].a[2] = (k[1]*k[2]*k[3]*k[4]*k[5]);
}
void TA::Graf(double* const Dsp, double& g3i,double& f6i,double& ce4,double& cs7) const{
//defines the trees for enzyme species
	Dsp[0] = g3i*kDe[0].a[0] + ce4*kDe[0].a[1];
	Dsp[1] = g3i*cs7*kDe[1].a[0] + g3i*f6i*kDe[1].a[1] + ce4*f6i*kDe[1].a[2];//ef
	Dsp[2] = g3i*cs7*kDe[2].a[0] + g3i*f6i*kDe[2].a[1] + ce4*f6i*kDe[2].a[2];//egg
	Dsp[3] = cs7*kDe[3].a[0] + f6i*kDe[3].a[1];//eg
	Dsp[4] = g3i*cs7*kDe[4].a[0] + ce4*cs7*kDe[4].a[1] + ce4*f6i*kDe[4].a[2];//ege
	Dsp[5] = g3i*cs7*kDe[5].a[0] + ce4*cs7*kDe[5].a[1] + ce4*f6i*kDe[5].a[2];//es
	Dsp[6]=Dsp[0];//e
	double sum=0; for(int i=0;i<6;i++) sum += Dsp[i];
	Dsp[7] = sum;
}
void TA::st1fl(double fl[], double& f6i,double& ce4, double& g3i,double& cs7) {
// calculates the elementary reaction rates and net flux
// for the TA reaction during the ODE solving.
	double Dsp[9];
	Graf(Dsp, g3i,f6i,ce4,cs7);
	for(int i=1;i<7;i++) {v[i]=Dsp[i-1]*e0*k[i]/Dsp[7];
		v0[i]=Dsp[i]*e0*k0[i]/Dsp[7];}
	fl[0] = (v[2] - v0[2]);	//flx16,f6p->g3p,e4p->s7p
}
void TA::st2fl(double* const fl,double& f6i,double& ce4, double& g3i,double& cs7) {
// calculates the elementary reaction rates and net flux
// for the TA reaction during the ODE solving.
	v[1]*=f6i; v0[3]*=g3i; v[4]*=ce4; v0[6]*=cs7;
	double Denom = v0[2]*v0[3]*v0[4]*v0[5]*v0[6] + v0[3]*v0[4]*v0[5]*v0[6]*v[1] + v0[4]*v0[5]*v0[6]*v[1]*v[2] + 
		v0[5]*v0[6]*v[1]*v[2]*v[3] + v0[6]*v[1]*v[2]*v[3]*v[4] + v[1]*v[2]*v[3]*v[4]*v[5]; 
        fl[0] = (v[1]*v[2]*v[3]*v[4]*v[5]*v[6])/Denom;		//f6i->cs7,flx41
	fl[1] = (v0[1]*v0[2]*v0[3]*v0[4]*v0[5]*v0[6])/Denom;    //cs7->f6i,flx42
		Denom = v0[2]*v0[3] + v0[3]*v[1] + v[1]*v[2];
	fl[2] = (v[1]*v[2]*v[3])/Denom - fl[0];				//f6i->g3i,flx43
		Denom = v0[5]*v0[6] + v0[6]*v[4] + v[4]*v[5]; 
	fl[3] = (v0[4]*v0[5]*v0[6])/Denom - fl[1];			//cs7->ce4,flx44
}
void Aldolase::setk(double *pv) {
e0=pv[0]; 
k[1]=pv[1]; 
k[3] = 0.0014*k[1]; 
k[2] = 0.0014*k[1]; 
k0[1] = 0.0014*k[1];
k0[2] = 0.376051*k[1]; 
k0[3] = 0.0376051*k[1]; 
	de[0]= k[2]*k[3]+k[3]*k0[1]; de[1]=k0[1]*k0[2];
	defbp[0] =k[1]*k[3]; defbp[1] =k0[2]*k0[3]; defbp[2] =k[1]*k0[2];
	dedh[0] = k[1]*k[2]; dedh[1] =k0[1]*k0[3]+k[2]*k0[3];
}
void Aldolase::st1fl(double fl[],double& cfbp, double& dhap) {
	De = de[0]+de[1]*dhap;
	Defbp = defbp[0]*cfbp+defbp[1]*dhap*dhap+defbp[2]*cfbp*dhap;
	Dedhap = dedh[0]*cfbp+dedh[1]*dhap;
	Denom = De+Defbp+Dedhap;
	v[1] = (De*e0*k[1]*cfbp)/Denom;
	v0[1] = (Defbp*e0*k0[1])/Denom;
	fl[0] = (v[1] - v0[1]);	
}
void Aldolase::st2fl(double* const fl,double& cfbp, double& dhap) {
	v[1] = (De*e0*k[1]*cfbp)/Denom;
	v0[1] = (Defbp*e0*k0[1])/Denom;
	v[2] = (Defbp*e0*k[2])/Denom;
	v0[2] = (Dedhap*e0*k0[2]*dhap)/Denom;
	v[3] = (Dedhap*e0*k[3])/Denom;
	v0[3] = (De*e0*k0[3]*dhap)/Denom;
	dedh[2] = (v[1]*v[2]+v[1]*v0[3]+v0[2]*v0[3]);
	fl[0] = v[1]*v[2]*v[3]/dedh[2];
	fl[1] = v0[1]*v0[2]*v0[3]/dedh[2];
	fl[2] =v[1]*v[2]/(v[1]+v0[2])-fl[0];
	fl[3] = v0[1]*v0[2]/(v[1]+v0[2])-fl[1];
}/**/


