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
inline void diff(double da,double dd[],double st[],int nex){for(int i=0;i<nex;i++) {dd[i] -= st[i]; dd[i]  /= da;}}
 
void Analis::hessian(vector<int> par1,int ndim,int np,int nex,double xi0,double st[], Mat_DP& aa, Mat_DP& dex,Vec_DP& b) {
     double fact(1.03),xi1,a,da;
     int i,j,k,l,nex1;
     for (i=np;i<ndim;i++) {
        a=Problem.rea[par1[i]].v();
        Problem.rea[par1[i]].setVm(a*fact);
        da= a*(fact - 1.) ;
       xi1= get<0>(solve());  nex1=horse.stor(&dex[i][0]); diff(da,&dex[i][0],st,nex);
        
        b[i] = -(xi1-xi0)/da/2.;
        Problem.rea[par1[i]].setVm(a);
             }
        for (j=np;j<ndim;j++){
          for (k=0;k<=j;k++) {
        aa[j][k]=0.;
     for (i=0;i<nex;i++) aa[j][k] += dex[j][i]*dex[k][i]; aa[k][j] = aa[j][k]; cout<<aa[j][k]<<"  "; } cout<<endl;
        }
     }

void Analis::grad() {
     int n=Problem.getparsize(), m(1), nex=horse.getmicon();
	int k, l, np(0), ndim(8);
      vector<int> par1=Problem.getFitPar();
     ofstream fi1("hes1.csv");
     Mat_DP A(n/2,n/2), D(n/2,nex);
	Problem.storeVms(nrea,nv1);//saves nv
          double xi0= get<0>(solve());
     double stcalc[nex]; nex=horse.stor(stcalc);
      for(int i=0;i<nex;i++) cout<<stcalc[i]<<"  ";  cout<<":stcalc"<<endl;
              cout<<"points: "<<nex<<endl;
   while(n>=ndim) { cout<<"par1: "; for(int i=0;i<ndim;i++) cout<<par1[i]<<" "; cout<<'\n';
Mat_DP aa(ndim,ndim),dd(ndim,nex),ai(ndim,ndim),u(ndim,ndim),v(ndim,ndim),uu(ndim,ndim),uuu(ndim,ndim);
          Vec_DP w(ndim),z(ndim),b(ndim);
            for (k=0;k<ndim;k++){
             for (l=0;l<=k;l++) aa[k][l]=aa[l][k]=A[k][l];
             for (l=0;l<nex;l++) dd[k][l]=D[k][l];
            }
        hessian(par1,ndim,np,nex,xi0,stcalc,aa,dd,b); 
        ai=aa;
            for (k=np;k<ndim;k++){
             for (l=0;l<ndim;l++) A[k][l]=aa[k][l];
             for (l=0;l<nex;l++) D[k][l]=dd[k][l];
            }
            
//          NR::gaussj(ai,x);
        NR::svdcmp(ai,w,v);
	revsvd(ai,w,v,uu);
        for (k=0;k<ndim;k++) z[k] = 1.0/w[k];
	revsvd(v,z,ai,u);
	matmult(u,uu, uuu);
        int flag(0);
   for (k=0;k<ndim;k++) {cout << setw(12) << uuu[k][k];
       if((uuu[k][k]<0.98)||(uuu[k][k]>1.02)) flag++;
       }   cout << endl;
       if(flag){ for (k=ndim;k<n;k++) par1[k-1] = par1[k]; n--;}
       else {for (k=0;k<ndim;k++) cout << par1[k]<<"; ";  cout <<endl;
      for (k=0;k<ndim;k++) {for (l=0;l<ndim;l++) fi1<<aa[k][l]<<" "; fi1<<endl;}
            ndim++;} np=ndim-1;
          while(Problem.rea[par1[ndim-1]].v()<1e-7) { for (k=ndim;k<n;k++) par1[k-1] = par1[k]; n--;}
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

