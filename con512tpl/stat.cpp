#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
using namespace std;

 inline double positive(double a){return sqrt(a*a);};
 
double Fit::dermax(){    double dx[numx];     f(xx,dx);
      double ma=positive(dx[0]); for(int i=1;i<numx;i++){
            double xpos=positive(dx[i]); if(xpos>ma) ma=xpos;}
                                             return ma;}
                                             
tuple<double,double,time_t> Fit::read(string fn){ 
 int i, ip; string aaa;  ifstream fi(fn.c_str());
 static int flg=0; tuple<double,double,time_t> sol(1.e12,0.,0);
  if(!fi.good()) {cout<<" bad file! ";return sol;}
    for (i=0;i<nrea;i++) rea[i].read(fi,flg); //cout<<"Nreact="<<nrea<<"; last="<<rea[nrea-1].v()<<'\n';
        par.clear();
	for(i=0;;i++) {fi>>ip; if(ip<0) break; par.push_back(ip);} getline(fi,spar);
	for (i=0;i<nmet;i++) {fi>>aaa>>namex[i]>>xx[i];} //cout<<nmet<<" "<<namex[nmet-1]<<" "<<xx[nmet-1]<<endl;
		fi >>aaa>> get<0>(sol)>>aaa >> get<2>(sol) >>aaa>> get<1>(sol);// cout<<aaa<<get<1>(sol)<<endl;
//        for (i=0;i<nflx;i++) fi>>aaa>>fid[i]>>flx[i];//cout<<nflx<<" "<<fid[nflx-1]<<" "<<flx[nflx-1]<<endl;
	fi.close(); flg++; //cout<<fn<<"; xi="<<xi<<endl;
return sol;}

void Fit::write (tuple<double,double,time_t> sol, int& ifn,bool flg)  {
         int i; stringstream fn;
            if(flg) ifn=setnumofi();
            if((ifn>fnfin)&&(ifn<(fnfin+5))) throw(invalid_argument("limit number of iterations reached"));
    fn<<outdir<<'/'<<ifn; //sprintf(fn,"%i",ifn);
  ofstream fi(fn.str().c_str());
    for (i=0;i<nrea;i++) rea[i].write(fi,i); 
    for (i=0;i<par.size();i++) fi << par[i]<<" "; fi << "-1"<<spar<<"\n";
    for(int i=0;i<nmet;i++) xx[i]=xinit1[i];                             //take initial values
 for (i=0;i<nmet;i++) fi<<i<<") "<<setw(9)<<left<<namex[i]<<" "<<xx[i]<<"\n";
		fi<<"χ= " << get<0>(sol) <<"\n";
		fi<<"t= " << get<2>(sol) <<"\n";
		fi<<"C= " << get<1>(sol) <<endl;
	for (i=0;i<rdt;i++) fi<<i<<") "<<setw(9)<<left<<rea[i].getname()<< " "<<(flx[i]*1000.*flx[rdt])<<"\n";
	fi.close();
cout <<"File saved: "<<ifn << ": xi=" << get<0>(sol) << "; xm=" << get<1>(sol) <<"; time="<<((float)get<2>(sol)/CLOCKS_PER_SEC)<<endl;	
}
 
void Parray::rnames(ifstream& fi){
   for (int i=0;i<nflx;i++)
        fi>>fid[i]>>fid[i]>>fname[i]>>fschem[i];
        }
void Fit::wstorefl (int numpar,const double** m,string fid[]) {
        ofstream fi(flmain.c_str());
   fi << "Confidence_level: 0.99\n Reaction_id Lower_bound Upper_bound\n";
//	        for (int j=0;j<parsets;j++) fi<<" "<<j; fi <<endl;
	 for (int i=0;i<numpar;i++) {	fi<<fid[i]<<" ";
	   double mn=m[0][i], mx=m[0][i]; int imn=0, imx=0;
	   for (int j=0;j<i99;j++) {
	    if(m[j][i]>mx) {mx=m[j][i]; imx=j;}
	    else if((m[j][i]<mn)&&(m[j][i]>1e-13)) {mn=m[j][i]; imn=j;}
                                   }
           fi  << mn <<" "<< mx<<"\n"; }
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
	 stringstream fn; fn<<outdir<<'/'<<(iset+1);
          fi.open(fn.str().c_str());
	for (;;) {getline(fi,aaa); if(aaa.find("C= ")+1) break;}
	
        for (i=0;i<nflx-2;i++) {fi>>aaa>>fid[i]>>mfl[iset][i];}
                     if(mfl[iset][50]>0) cout<<" set="<<iset<<" fl="<<mfl[iset][50]<<'\n';
		fi >> mfl[iset][nflx]>>mfl[iset][nflx+1]>>mfl[iset][nflx+2];
	fi.close();
 }
         wstorefl(nflx-2, pmf,fid);
}
void Fit::stat(const int NP ){ //reorders parameters files ascending with respect to chi square
        Vec_DP xi(NP), conc(NP), *ac; Vec_INT b(NP),t(NP);
        int i,sys; cout<<"NP= "<<NP<<'\n';
 for ( i=1;i<=NP;i++) {
	  stringstream fn; fn<<outdir<<'/'<<i; 
       tie(xi[i-1],conc[i-1],t[i-1]) = read(fn.str().c_str());
//       xi[i-1] = read(t[i-1],conc[i-1],fn.str().c_str());
       b[i-1] = i;
 }
        ac=&xi;
//        ac=&conc;
        NR::sort2a(*ac,b);
        cout << endl << "After sorting by xi2, Parameter files are:" << '\n';
	cout.precision(4);
	for (i=0;i<NP;i++) {
		double df=(xi[i]-xi[0]);
		if (df<1.*1.0) {i68++; i90++; i95++; i99++;}
		else if (df<1.*2.71) {i90++; i95++; i99++;}
		else if (df<1.*4.0) { i95++; i99++;}
		else if (df<1.0*6.63)  i99++;
		else break;
	}
	for(i=0;i<i99;i++){
cout<<setw(3)<<b[i]<<" ("<<xi[i]<<"; "<<conc[i]<<"; "<<t[i]<<") "; 
if ((i%3)==0) cout<<'\n';
	    }
	cout<<endl;
readst( &b[0]);
int chr=NP;//i99;
 for (i=0;i<chr;i++) {
	  stringstream fn; fn<<"mv "<<outdir<<'/'<<b[i]<<" "<<outdir<<'/'<<(i+1)<<"a";
		sys=system(fn.str().c_str());
	}
 for (i=0;i<chr;i++) {
	  stringstream fn; fn<<"mv "<<outdir<<'/'<<(i+1)<<"a "<<outdir<<'/'<<(i+1);
		sys=system(fn.str().c_str());
	}
	cout<<"selected " <<i99<<" from "<<NP<<endl;
}

