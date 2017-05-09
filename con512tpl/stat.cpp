#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
using namespace std;

/*double Parray::sumx() {
 return (xx[nh6]+xx[nt3]+xx[nfbp]+xx[npep]+xx[npyr]+xx[noa]+ xx[noa1]+ xx[nmal]+ xx[nakg]+ xx[nakg1]+ xx[ncit]+ xx[ncit1]+ xx[ncoa]+ xx[np5]+ xx[ne4]+ xx[ns7]);}*/
 inline double positive(double a){return sqrt(a*a);};
 
double Fit::dermax(){    double dx[numx];     f(xx,dx);
      double ma=positive(dx[0]); for(int i=1;i<numx;i++){
            double xpos=positive(dx[i]); if(xpos>ma) ma=xpos;}
                                             return ma;}
                                             
double Fit::read(int &t,double &c, string fn){
 int i; string aaa;  ifstream fi(fn.c_str()); string ppp;
 static int flg=0;
  if(!fi.good()) return (1.e12);
    for (i=0;i<nrea;i++) rea[i].read(fi,flg); 
	for(i=1;;i++) {fi>>par[i]; if(par[i]<0) {par[0]=i; break;}}
	for (i=0;i<numx;i++) {fi>>aaa>>namex[i]>>xx[i];}
        for (i=0;i<nflx;i++) fi>>aaa>>namef[i]>>flx[i]; if(flx[pfk]<1e-7) cout<<fn<<"!!!: hk=0"<<endl; 
		fi >> xi >> t >> c;// if((flx[0]<0.1)||(flx[0]>0.2)) xi += 100.;
	fi.close(); flg++; //cout<<fn<<"; xi="<<xi<<endl;
return xi;}

void Fit::write (time_t tf, int& ifn,const double xi0,const double xm,bool flg) const {
         int i; stringstream fn;
            if(flg) ifn=setnumofi();  if((ifn>fnfin)&&(ifn<(fnfin+5))) throw(invalid_argument("limit number of iterations reached"));
    fn<<outdir<<ifn; //sprintf(fn,"%i",ifn);
  ofstream fi(fn.str().c_str());
    for (i=0;i<nrea;i++) rea[i].write(fi,i); 
	for (i=1; i<par[0]; i++) fi << par[i]<<" "; fi << "-1" << endl;
  for (i=0;i<numx;i++) fi<<i<<") "<<setw(9)<<left<<namex[i]<<" "<<xx[i]<<endl;
	for (i=0;i<nflx;i++) fi<<i<<") "<<setw(9)<<left<<namef[i]<< " "<<(flx[i]*1000.*dt)<<endl;
		fi << xi0 <<endl;
		fi << tf <<endl;
		fi << xm <<endl;
		fi << dif <<endl;
	fi.close();
cout <<"File saved: "<<ifn << ": xi=" << xi0 << "; xm=" << xm <<"; time="<<((float)tf/CLOCKS_PER_SEC)<<"; dif="<<dif<<endl;	
}

void Fit::wstorefl (const char fn1[],int numpar,const double** m,string name[]) {
        ofstream fi(fn1);
   fi << " 95% +- 99% +- bestfit fimin fimax xmin xmax" << endl;
//	        for (int j=0;j<parsets;j++) fi<<" "<<j; fi <<endl;
	        for (int i=0;i<numpar;i++) {
			fi<<name[i]<<" ";
//			cout<<setw(11)<<name[i]<<" ";
	double mn=m[0][i], mx=m[0][i]; int imn=0, imx=0;
	        for (int j=0;j<i95;j++) {
			if(m[j][i]>mx) {mx=m[j][i]; imx=j;}
			else if((m[j][i]<mn)&&(m[j][i]>1e-13)) {mn=m[j][i]; imn=j;}
		 }   fi  << mn <<" "<< mx <<" ";
		 
	        for (int j=i95;j<i99;j++) {
			if(m[j][i]>mx) {mx=m[j][i]; imx=j;}
			else if((m[j][i]<mn)&&(m[j][i]>1e-13)) {mn=m[j][i]; imn=j;}
 }   fi  << mn <<" "<< mx <<" "<< m[0][i]<<" "<<imn<<" "<<imx <<" ";
   fi <<(m[imn][nflx]-m[0][nflx])<<" "<<(m[imx][nflx]-m[0][nflx])<<endl;
//   cout << (mx+mn)/2. <<" "<< (mx-mn)/2. <<" "<< m[0][i] <<endl;
                }
   fi<<"chi         "<<m[0][numpar]; fi <<endl;
   fi<<"time         "<<m[0][numpar+1]; fi <<endl;
   fi<<"conc         "<<m[0][numpar+2]; fi <<endl;

	fi.close();
}
void Fit::readst( int* b){
          int i,tm; char fn[15]; string aaa;
          static int flg=0;
        double mpar[i99][nrea+3], mfl[i99][nflx+3];
        const double *pmp[i99],*pmf[i99];
	ifstream fi;
 for (int iset=0;iset<i99;iset++) {
	pmp[iset]=&mpar[iset][0]; pmf[iset]=&mfl[iset][0];
	 stringstream fn; fn<<outdir<<(iset+1);
          fi.open(fn.str().c_str());
	for (i=0;i<nrea;i++) {fi>>aaa>>aaa>>aaa>>aaa>>mpar[iset][i]; getline(fi,aaa);}
	getline(fi,aaa);
	for (i=0;i<numx;i++) fi>>aaa>>namex[i]>>xx[i];
        for (i=0;i<nflx;i++) fi>>aaa>>namef[i]>>mfl[iset][i];
/*      mfl[iset][mdh_net]=mfl[iset][maloa]-mfl[iset][oamal];
        mfl[iset][aldfl+2]=mfl[iset][aldfl]-mfl[iset][aldfl+1];
        mfl[iset][ckg_net]=mfl[iset][citakg1]-mfl[iset][idhr];
        mfl[iset][kgd_net]=mfl[iset][akgdf]-mfl[iset][akgdr];
        mfl[iset][lac_net]=mfl[iset][laci_o]-mfl[iset][laco_i];
        mfl[iset][glu_net]=mfl[iset][glu_out]-mfl[iset][glu_in];
mfl[iset][nadhf]=mfl[iset][pdh]+mfl[iset][akgfum]+(mfl[iset][citakg]+mfl[iset][r5in])*0.5;*/
		fi >> mfl[iset][nflx]>>mfl[iset][nflx+1]>>mfl[iset][nflx+2];
		mpar[iset][nrea]= mfl[iset][nflx];
                mpar[iset][nrea+1]= mfl[iset][nflx+1];
                mpar[iset][nrea+2]= mfl[iset][nflx+2];
	fi.close();
 } cout<<namef[0]<<endl;
//    cout<<"Parameters:"<<setw(11)<<"mean"<<setw(11)<<"SD"<<setw(11)<<"SE"<<endl;
         wstorefl("statpar", nrea, pmp,namef);
//    cout<<"\nFluxes:    "<<setw(11)<<"mean"<<setw(11)<<"SD"<<setw(11)<<"SE"<<endl;
         wstorefl("statfl", nflx, pmf,namef);
}
void Fit::stat(const int NP ){
        Vec_DP a(NP), conc(NP); Vec_INT b(NP),t(NP);
        int* pb=&b[0];
        int i,sys;
 for ( i=1;i<=NP;i++) {
	  stringstream fn; fn<<outdir<<i; cout<<fn.str().c_str()<<endl;
       a[i-1] = read(t[i-1],conc[i-1],fn.str().c_str());
       b[i-1] = i;
 }
        NR::sort2a(a,b);
        cout << endl << "After sorting, Parameter file (xi2) are:" << endl;
	cout.precision(4);
	for (i=0;i<NP;i++) {
		double df=(a[i]-a[0]);
		if (df<1.*1.0) {i68++; i90++; i95++; i99++;}
		else if (df<1.*2.71) {i90++; i95++; i99++;}
		else if (df<1.*4.0) { i95++; i99++;}
		else if (df<1.0*6.63)  i99++;
		else break;
	}
	for(i=0;i<i99;i++){
cout<<setw(3)<<b[i]<<" ("<<a[i]<<"; "<<conc[i]<<"; "<<t[i]<<") "; 
if ((i%3)==0) cout<<endl;
	    }
	cout<<endl;
readst( pb);
int chr=NP;//i99;
 for (i=0;i<chr;i++) {
	  stringstream fn; fn<<"mv "<<outdir<<b[i]<<" "<<outdir<<(i+1)<<"a";
		sys=system(fn.str().c_str());
	}
 for (i=0;i<chr;i++) {
	  stringstream fn; fn<<"mv "<<outdir<<(i+1)<<"a "<<outdir<<(i+1);
		sys=system(fn.str().c_str());
	}
	cout<<"selected " <<i99<<" from "<<NP<<endl;
}

