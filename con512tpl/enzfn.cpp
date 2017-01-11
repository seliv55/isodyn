#include <iostream>
#include "nums.hh"
#include "modlab.h"


//---------------------------------------------------------------------------
using namespace std;
void Ldistr::split(double *h6,double *dh6,double *t1,double *dt1,const double vf,const double vr) {
	int k,k1,i; double x;
	for(i=0;i<64;i++) {
		int ifr = (i >> 3);
		if (ifr==1) k1=4;
		else if (ifr==3) k1=6;
		else if (ifr==4) k1=1;
		else if (ifr==6) k1=3;
		else k1=ifr;
		k = (i & 7);
            x = vf*(h6[i]) - vr*(t1[k])*(t1[k1]);
		dh6[i] -= x; dt1[k] += x;
		dt1[k1] += x;
	}
}
void Ldistr::spInvsl(double *h6,double *dh6,double *t1,double *dt1,const double vf,const double tsum) {
	int k,k1,i,ifr,j; double x,xb,xt;
	for(i=0;i<64;i++) {
		ifr = (i >> 3); k = (i&7);
		if (ifr==1) k1=4;
		else if (ifr==3) k1=6;
		else if (ifr==4) k1=1;
		else if (ifr==6) k1=3;
		else k1=ifr;
            x = vf*(h6[i]);
		dh6[i] -= x; 
		dt1[k1] += x;
		xb = x/tsum;
		for (j=0;j<8;j++) {
			xt = xb*t1[j];
		if (j==1) k1=4;
		else if (j==3) k1=6;
		else if (j==4) k1=1;
		else if (j==6) k1=3;
		else k1=j;
		ifr = ((k1<<3)|k);
		dt1[j] -= xt;
		dh6[ifr] += xt;
		}
	}
}
void Ldistr::csyn(double *coai,double *dcoai,double *oaai,double *doaai,double *dciti,const double v) {
	int i1,i,j,ic; double x,dx;
	for(i=0;i<4;i++) {
		i1=(i<<4); x=coai[i]*v;
                j=0;
		while(j<16) {
	 		dx=x*oaai[j];
			dcoai[i] -= dx; doaai[j] -= dx;
			ic=(j|i1); dciti[ic] += dx;j++;
		}
	}
}
void Ldistr::symm (double *s) {
	int io[6]={1, 2, 3,  5,  7,  11};
	int oi[6]={8, 4, 12, 10, 14, 13};
	for(int i=0;i<6;i++) {
		double x=(s[io[i]]+s[oi[i]])/2.;
		s[io[i]]=x;	s[oi[i]]=x;
	}
}
double Ldistr::ser_gly (double v,double xthf) {
 double x,sum(0.); int i1, i2;
  for(int i=0;i<8;i++) { i1=(i>>1); i2=(i&1);
     x=v*ser.iso[i]*xthf;
       ser.diso[i] -=x; gly.diso[i1]+=x; cthf.diso[i2]+=x; sum+=x; }
    return sum;
}
double Ldistr::gly_ser (const double v) {
 double x,x1,sum(0.); int i1;
 for(int j=0;j<2;j++){x1=v*cthf.iso[j]; 
 for(int i=0;i<4;i++) { i1=((i<<1)|j);
     x=gly.iso[i]*x1;// cout<<"x="<<x<<endl;
          gly.diso[i] -=x; cthf.diso[j]-=x; ser.diso[i1]+=x; sum+=x; }}
    return sum;
}

