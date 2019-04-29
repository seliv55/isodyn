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


void Analis::sensitivity(double *grd,bool *gsign,vector<int>& par1,const double fac,vector<int>& sens){
  int ngrd=par1.size(); tuple<double,double,time_t> sol; double x01,pnew;
   cout<<"\nSensitivity: x00="<<x00<<'\n'; trs=1.05;
   cout<<"par1: "; for(int i=0;i<par1.size();i++) cout<<par1[i]<<" "; cout<<'\n';
  for(int i=0;i<ngrd;i++)
   if(nv1[par1[i]]>1e-7){
             pnew=nv1[par1[i]]*fac; cout<<par1[i];
             Problem.rea[par1[i]].setVm(pnew);
       try{ sol =solve();
            x01 = get<0>(sol);
            grd[i]=(x01-x00)/(fac-1.);
//            grd[i]=x01*suxx/x00/xmin;
         if(grd[i]<=0.) {
              gsign[i]=false;
            }
         else { pnew=nv1[par1[i]]/fac; Problem.rea[par1[i]].setVm(pnew);
                sol =solve();
                gsign[i]=true;
                x01 = get<0>(sol);
            grd[i]=(x01-x00)/(fac-1.);
//                grd[i]=x01*suxx/x00/xmin;
            } string saa;
            if(grd[i]<(-trs)) {saa="\t-"; sens.push_back(par1[i]); }
            else if(grd[i]>trs) {saa="\t+"; sens.push_back(par1[i]); }
            else if((grd[i]<0)&&(suxx<xmin)&&(tf<tmin)) {saa="\t*"; setx00(x01,tf,suxx); nv1[par1[i]]=pnew;Problem.write(sol,ifn); }
            else saa="\t";
          cout<<saa<<"xi="<<x01<<" sux="<<suxx<<" t="<<tf<<saa<<"GRAD="<<grd[i]<<"\t"<<gsign[i]<< "\tpnew="<<pnew<<'\n';
       }catch( char const* str ){cout << "Analis::sensitivity: "<< str <<endl; Problem.restoreVm(nrea,nv2);
                 for(int i=0;i<numx;i++) xinit1[i]=xinit2[i]; }
          Problem.rea[par1[i]].setVm(nv1[par1[i]]);
  }
  else grd[i]=0.;
}
double Analis::checkgroup(bool* gsign, double* grd,int imi,vector<int>& vpar,vector<int>& vind,double ff){
        tuple<double,double,time_t> sol; double n_v,xi;
           cout<<"\n check group, grd_max=:"<<grd[imi]<<'\n';
    for(int i=0;i<vpar.size();i++) { double fac=(ff-1)*grd[vind[i]]/grd[imi] +1.;
       n_v=(gsign[vind[i]]) ? nv1[vpar[i]]/fac : nv1[vpar[i]]*fac;
       Problem.rea[vpar[i]].setVm(n_v);
    try{ 
      sol =solve(); xi=get<0>(sol);
      cout<<"  par["<<vpar[i]<<"]="<<nv1[vpar[i]]<<" - "<<n_v<<"\txi="<<xi<<'\n';
      if((xi/x00<0.999)&&(suxx/xmin<1.1)&&(tf/tmin<1.3)) {Problem.write(sol,ifn); nv1[vpar[i]]=n_v; setx00(xi,tf,suxx);}
       else {Problem.rea[vpar[i]].setVm(nv1[vpar[i]]);}
    }catch( char const* str ){cout << "Analis::checkgroup: "<< str <<endl; Problem.restoreVm(nrea,nv2);
                 for(int i=0;i<numx;i++) xinit1[i]=xinit2[i];}
  }  
 return x00;
}
int Analis::select_pars(double *grd,vector<int>&par1,vector<int>& vpar,vector<int>& vind){
  int imi(-1),ngrd=par1.size(),nvpar(0); double gmin(1.),thr(-2.5); cout<<"Selected parameters: ";
  for(int i=0;i<ngrd;i++) if(grd[i]<thr) {
     vpar.push_back(par1[i]); vind.push_back(i); cout<<par1[i]<<" ";
     if(grd[i]<gmin) {gmin=grd[i]; imi=i;}
    }
  if((nvpar=vpar.size())>0){cout<<'\n'<<"par[imin]="<<par1[imi]<<" grd="<<gmin<<" v="<<Problem.rea[par1[imi]].v()<<'\n';
  for(int i=0;i<nvpar;i++) { cout<<'\n'<<vpar[i]<<" "<<grd[vind[i]];}
    }
return imi;}

void Analis::grdesc(const double lambda) { 
  cout<<"Gradient descent\n";
  int nex=horse.getmicon();  Problem.storeVms(nrea,nv1);
  vector<int> par1=Problem.getFitPar(); 
//  tuple<double,double,time_t> sol=solve(); x00 = get<0>(sol); xmin=suxx; tmin=tf;
//  steps of reduction of sensitive parameters:
  for(;;){  int ngrd=par1.size(); double grd[ngrd]; bool gsign[ngrd];
    vector<int> vpar,vind,sens;
//  check sensitivity:
    sensitivity(grd,gsign,par1,lambda,sens);
//  select pars of max sensitivity
    int imi=select_pars(grd,par1,vpar,vind);
   if(imi<0) return;
   par1=sens;
//  checking
    chimin=x00;
    checkgroup(gsign, grd, imi,vpar,vind,lambda);
   if(x00>(chimin-0.5)) return;
    
   for(;;){ chimin=x00;
    checkgroup(gsign, grd, imi,vpar,vind,lambda);
   if(x00>(chimin-1.5)) break;}
  }
}
 
double Analis::checkmax(bool bsig,double& oldv,const double& f1,int& parm){
 cout<<"\tcheck max:"; tuple<double,double,time_t> sol;
  double xi, newv1;
//  initial check of f1
     newv1= ((bsig) ? oldv/f1 : oldv*f1);
     Problem.rea[parm].setVm(newv1); 
     try{ sol =solve(); xi=get<0>(sol);  cout<<" xi="<<xi<<'\n';
     if((xi/x00<0.999)&&(suxx/xmin<1.1)&&(tf/tmin<1.3)) { Problem.write(sol,ifn); oldv=newv1; x00=xi;}
         else {Problem.rea[parm].setVm(oldv); }
     }catch( char const* str ){cout << "Analis::checkmax: "<< str <<endl; Problem.restoreVm(nrea,nv2);
                 for(int i=0;i<numx;i++) xinit1[i]=xinit2[i];}
 return x00;}

//double Analis::grdes(bool* gsign, double* grd,int imi,vector<int>& vpar,vector<int>& vind,double ff,double* oval){
//        tuple<double,double,time_t> sol; double n_v,xi;
//           cout<<" check grad:"<<'\n'; Problem.storeVms(nrea,nv1);
//    for(int i=0;i<vpar.size();i++) { double fac=(ff-1)*grd[vind[i]]/grd[imi] +1.;
//       n_v=(gsign[vind[i]]) ? oval[vind[i]]/fac : oval[vind[i]]*fac;
//       Problem.rea[vpar[i]].setVm(n_v);
//  }  
//    try{ 
//      sol =solve(); xi=get<0>(sol);
//    }catch( char const* str ){cout << "Analis::grdes: "<< str <<endl; Problem.restoreVm(nrea,nv2);
//                 for(int i=0;i<numx;i++) xinit1[i]=xinit2[i];}
//         if((xi/x00<0.999)&&(suxx/xmin<1.1)&&(tf/tmin<1.3)) {Problem.write(sol,ifn); x00=xi;}
//       else {Problem.restoreVm(nrea,nv1);}
// return x00;
//}

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

