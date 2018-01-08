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
                                             
tuple<double,double,time_t> Fit::read(string fn){
 int i, ip; string aaa;  ifstream fi(fn.c_str()); string ppp;
 static int flg=0; tuple<double,double,time_t> sol(1.e12,0.,0);
  if(!fi.good()) return sol;
    for (i=0;i<nrea;i++) rea[i].read(fi,flg);// cout<<"Nreact="<<nrea<<"; last="<<rea[nrea-1].v()<<endl;
	for(i=0;;i++) {fi>>ip; if(ip<0) break; par.push_back(ip);} getline(fi,spar);
	for (i=0;i<nmet;i++) {fi>>aaa>>namex[i]>>xx[i];} //cout<<aaa<<" "<<namex[nmet-1]<<" "<<xx[nmet-1]<<endl;
        for (i=0;i<nflx;i++) fi>>aaa>>fid[i]>>flx[i]; if(flx[pfk]<1e-7) cout<<fn<<"!!!: hk=0"<<endl; 
		fi >> get<0>(sol) >> get<2>(sol) >> get<1>(sol);// if((flx[0]<0.1)||(flx[0]>0.2)) xi += 100.;
	fi.close(); flg++; //cout<<fn<<"; xi="<<xi<<endl;
return sol;}

void Fit::write (tuple<double,double,time_t> sol, int& ifn,bool flg) const {
         int i; stringstream fn;
            if(flg) ifn=setnumofi();  if((ifn>fnfin)&&(ifn<(fnfin+5))) throw(invalid_argument("limit number of iterations reached"));
    fn<<outdir<<ifn; //sprintf(fn,"%i",ifn);
  ofstream fi(fn.str().c_str());
    for (i=0;i<nrea;i++) rea[i].write(fi,i); 
	for (i=0;i<par.size();i++) fi << par[i]<<" "; fi << "-1"<<spar<<"\n";
  for (i=0;i<nmet;i++) fi<<i<<") "<<setw(9)<<left<<namex[i]<<" "<<xx[i]<<"\n";
	for (i=0;i<nflx;i++) fi<<i<<") "<<setw(9)<<left<<fid[i]<< " "<<(flx[i]*1000.*dt)<<"\n";
		fi << get<0>(sol) <<"\n";
		fi << get<2>(sol) <<"\n";
		fi << get<1>(sol) <<endl;
	fi.close();
cout <<"File saved: "<<ifn << ": xi=" << get<0>(sol) << "; xm=" << get<1>(sol) <<"; time="<<((float)get<2>(sol)/CLOCKS_PER_SEC)<<endl;	
}
 
void Parray::rnames(ifstream& fi){
   for (int i=0;i<nflx;i++)
        fi>>fid[i]>>fid[i]>>fname[i]>>fschem[i];
        }

void Fit::wstorefl (int numpar,const double** m,string name[]) {
     stringstream finame; finame<<outdir<<"names"; ifstream fii(finame.str().c_str());
     rnames(fii);
        ofstream fi(flmain.c_str());
   fi << "Confidence_level: 0.99\n Reaction_id Lower_bound Upper_bound name scheme\n";
//	        for (int j=0;j<parsets;j++) fi<<" "<<j; fi <<endl;
	        for (int i=0;i<numpar;i++) {
			fi<<fid[i]<<" ";
//			cout<<setw(11)<<name[i]<<" ";
	double mn=m[0][i], mx=m[0][i]; int imn=0, imx=0;
//	        for (int j=0;j<i95;j++) {
//			if(m[j][i]>mx) {mx=m[j][i]; imx=j;}
//			else if((m[j][i]<mn)&&(m[j][i]>1e-13)) {mn=m[j][i]; imn=j;}
//		 }   fi  << mn <<" "<< mx <<" ";
		 
	        for (int j=0;j<i99;j++) {
			if(m[j][i]>mx) {mx=m[j][i]; imx=j;}
			else if((m[j][i]<mn)&&(m[j][i]>1e-13)) {mn=m[j][i]; imn=j;}
 }   fi  << mn <<" "<< mx <<" "<<fname[i]<<" "<<fschem[i]<<"\n";
//   fi <<(m[imn][nflx]-m[0][nflx])<<" "<<(m[imx][nflx]-m[0][nflx])<<endl;
//   cout << (mx+mn)/2. <<" "<< (mx-mn)/2. <<" "<< m[0][i] <<endl;
                }
   fi<<"chi "<<m[0][numpar] <<"\n";
   fi<<"time "<<m[0][numpar+1] <<"\n";

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
        for (i=0;i<nflx;i++) fi>>aaa>>fid[i]>>mfl[iset][i];
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
 } cout<<fid[0]<<endl;
//    cout<<"\nFluxes:    "<<setw(11)<<"mean"<<setw(11)<<"SD"<<setw(11)<<"SE"<<endl;
         wstorefl(nflx, pmf,fid);
}
void Fit::stat(const int NP ){
        Vec_DP a(NP), conc(NP), *ac; Vec_INT b(NP),t(NP);
        int i,sys;
 for ( i=1;i<=NP;i++) {
	  stringstream fn; fn<<outdir<<i; cout<<fn.str().c_str()<<endl;
       tie(a[i-1],conc[i-1],t[i-1]) = read(fn.str().c_str());
//       a[i-1] = read(t[i-1],conc[i-1],fn.str().c_str());
       b[i-1] = i;
 }
        ac=&a;
//        ac=&conc;
        NR::sort2a(*ac,b);
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
readst( &b[0]);
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

