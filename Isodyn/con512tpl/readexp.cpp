#include <iostream>
#include <cmath>
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
using namespace std;
double Vi, xribi, xasp=1.,mu;
double Ldistr::readExp (string fn) {
string aaa;  ifstream fi(fn.c_str()); double Ti,ts1;  mu=0.; dt=0.6;  Vi=0.014; tex[0]=0.;  
      for(ntime=1;;ntime++) {fi>>tex[ntime]; tex[ntime] *= 60.; if(tex[ntime]<0) break;}//time points
 double Nc[ntime]; for(int i=0;i<ntime;i++) fi>>Nc[i]; getline(fi,aaa);//cells number
 
for(int i=1;i<ntime;i++) mu += log(Nc[i]/Nc[0])/tex[i];
 mu /= ((double)(ntime-1)*1.0); getline(fi,aaa);
	 xx[ngl]=10.; gl.setconc(xx[ngl]); gl.setex0();
	 xx[nlac]=0.001; lac.setconc(xx[nlac]);
	 xx[nglu]=0.1; glu.setconc(xx[nglu]); glu25.setconc(xx[nglu]);
	 xx[ngln]=4.9; gln.setconc(xx[ngln]); gln.setex0();
	 xx[nala]=0.4; ala.setconc(xx[nala]);
	 xx[ngly]=0.4; gly.setconc(xx[ngly]);
	 xx[nser]=0.4; ser.setconc(xx[nser]);
	 xx[nasp]=0.0001; asp.setconc(xx[nasp]);
	 xx[npro]=0.0001; pro.setconc(xx[npro]);
	 xx[ncoac]=0.007; coac.setconc(xx[ncoac]);
	 xx[nrna]=0.007; coac.setconc(xx[nrna]);

	 xx[nagl]=1.0;
	 xx[ncthf]=0.5;	 lmet=sizeof(met)/sizeof(*met);
	 bool a=0;
    while(!fi.eof()){
      fi>>aaa; if(fi.eof()) break; cout<<"descr: "<<aaa<<endl;
   int k(0); for(int i=0;i<lmet;i++){
      int imatch=aaa.find(met[i]->getdescr());
        if(imatch != (int)std::string::npos){met[i]->setex0();
         cout<<i<<": "<<met[i]->getdescr()<<endl;
          met[i]->read(fi,1); k++; break;}
                             }
           if(k==0){getline(fi,aaa);getline(fi,aaa);cout<<"k="<<k<<" aaa="<<aaa<<endl;}
    }
     gl.setrav(0,1);
return ts1;}
