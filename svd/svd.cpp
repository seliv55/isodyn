#include <iostream>
#include "nr.h"
//#include "bdf.hpp"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
#include "analis.h"
//---------------------------------------------------------------------------
//u=aa*ai
void Analis::matmult(Mat_I_DP &aa, Mat_I_DP &ai,Mat_O_DP &u){
	int k,l,j,n=aa.nrows();
          for (k=0;k<n;k++) {
            for (l=0;l<n;l++) {
              u[k][l]=0.0;
              for (j=0;j<n;j++)
                u[k][l] += (aa[k][j]*ai[j][l]);
            }
          }
}
void Analis::revsvd(Mat_I_DP &aa,Vec_I_DP &w, Mat_I_DP &ai,Mat_O_DP &u){
	int k,l,j,n=aa.nrows();
          for (k=0;k<n;k++) {
            for (l=0;l<n;l++) {
              u[k][l]=0.0;
              for (j=0;j<n;j++) 
                u[k][l] += (aa[k][j]*ai[l][j]*w[j]);
            }
          }
}
inline void matdisp(Mat_I_DP &ai){
	int k,l,n=ai.nrows();
          for (k=0;k<n;k++) {
            for (l=0;l<n;l++) cout << setw(12) << ai[k][l];
            cout << endl;
          }
}
void Analis::hessian(int ndim,int np,int nex,double xi0,double st[], Mat_DP& aa,Vec_DP& b,int par1[],int nt) {
     double fact(1.03),xi1,a,da;
     int i,j,k,l,nex1;
          Mat_DP dex(ndim,nex);
	Problem.restoreVm(nrea,nv1);//gets nv
	for(int i=0;i<numx;i++) xx[i]=xinit1[i];//gets xx
     for (i=np;i<ndim;i++) {
        a=Problem.rea[par1[i]].v();
        Problem.rea[par1[i]].setVm(a*fact);
        da= a*(fact - 1.) ;
       xi1= get<0>(solve());  nex1=horse.stor(&dex[i][0],nt); horse.diff(da,&dex[i][0],st,nex);
        
        b[i] = -(xi1-xi0)/da/2.;
	Problem.restoreVm(nrea,nv1);//gets nv
	for(int i=0;i<numx;i++) xx[i]=xinit1[i];//gets xx
             }
        for (j=np;j<ndim;j++)
          for (k=0;k<=j;k++) {
        aa[j][k]=0.;
     for (i=0;i<nex;i++) {aa[j][k] += dex[j][i]*dex[k][i]; aa[k][j] = aa[j][k];}
     }
}
void Analis::grad(int nt) {
     int n=Problem.getparsize(), m(1), nex=horse.getmicon();
	int par1[n], k, l, np(0), ndim(8);
     double stcalc[nex];
     ofstream fi1("hes1.csv");
     Mat_DP A(n,n);
	Problem.storeVms(nrea,nv1);//saves nv
	for(int i=0;i<numx;i++) xinit1[i]=xx[i];//saves xx
          double xi0= get<0>(solve()); nex=horse.stor(stcalc,nt);
              cout<<"points: "<<nex<<endl;
   while(n>=ndim) {
Mat_DP aa(ndim,ndim),ai(ndim,ndim),u(ndim,ndim),v(ndim,ndim),uu(ndim,ndim),uuu(ndim,ndim);
          Vec_DP w(ndim),z(ndim),b(ndim);
            for (k=0;k<ndim;k++) for (l=0;l<=k;l++) aa[k][l]=aa[l][k]=A[k][l];
        hessian(ndim,np,nex,xi0,stcalc,aa,b,par1,nt);
        ai=aa; 
//          NR::gaussj(ai,x);
        NR::svdcmp(ai,w,v);
	revsvd(ai,w,v,uu);
        for (k=0;k<ndim;k++) z[k] = 1.0/w[k];
	revsvd(v,z,ai,u);
	matmult(u,uu, uuu);
        int flag(0);
   for (k=0;k<ndim;k++) {cout << setw(12) << uuu[k][k];
       if((uuu[k][k]>0.999)&&(uuu[k][k]<1.001));  else flag++;
       }   cout << endl;
       if(flag){ for (k=ndim;k<n;k++) par1[k-1] = par1[k]; n--;}
       else {for (k=0;k<ndim;k++) cout << par1[k]<<"; "; 
            for (k=np;k<ndim;k++) for (l=0;l<ndim;l++) A[k][l]=aa[k][l];
            cout <<endl; np=0;
      for (k=0;k<ndim;k++) {for (l=0;l<ndim;l++) fi1<<aa[k][l]<<" "; fi1<<endl;}
            ndim++;}
      }
      ndim--; fi1.close();/* fi1.open("hes2.csv");
      for (k=0;k<ndim;k++) cout << par1[k]<<"; ";
            cout << endl<<" ndim="<<ndim<<"; n="<<n <<"; np="<<np<< endl;
Mat_DP aa(ndim,ndim),v(ndim,ndim),u(ndim,ndim),uu(ndim,ndim),uuu(ndim,ndim);
          Vec_DP w(ndim),z(ndim),b(ndim);
        hessian(ndim,0,aa,b,par1,tmax);
      for (k=0;k<ndim;k++) {for (l=0;l<ndim;l++) fi1<<aa[k][l]<<" "; fi1<<endl;}
        NR::svdcmp(aa,w,v);
	revsvd(aa,w,v,uu);
        for (k=0;k<ndim;k++) z[k] = 1.0/w[k];
	revsvd(v,z,aa,u);
	matmult(u,uu, uuu);
        cout << endl << "diag c=a^(-1):" << endl;
   for (k=0;k<ndim;k++) cout << setw(12) << (sqrt(u[k][k])/Problem.getVal(par1[k]));//diag a^(-1) 
        cout << endl << "diag a*a^(-1):" << endl;
   for (k=0;k<ndim;k++) {cout << setw(12) << uuu[k][k];
       }   cout << endl;
        DP wmax=0.,wmin;
          cout << "diag W:" << endl;
            for (l=0;l<ndim;l++) {if(wmax<w[l]) wmax=w[l]; cout << setw(12) << w[l];}//W  
            cout << endl;
          wmin=wmax*(1.0e-11);
          for (k=0;k<ndim;k++){
            if (w[k] < wmin) w[k]=0.0; cout << setw(12) << w[k];
            }//W 
            cout << endl;
            NR::svbksb(aa,w,v,b,z);
            cout << " solution vector is:" << endl;//solution and substitution
     for (k=0;k<ndim;k++) {
     cout << setw(12) << z[k]; Problem.setVal(par1[k], (Problem.getVal(par1[k])+z[k]));
     }
     cout << endl;
     cout <<"chi square="<< solve(tmax)<<endl;*/
}

