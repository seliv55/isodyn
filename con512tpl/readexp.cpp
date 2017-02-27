#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
using namespace std;
double Vi, xribi, xasp=1.,mu;

double Ldistr::readExp (char fn[]) {
  ifstream fi(fn); double Ti,ts1;  mu=0.; dt=0.6;  Vi=0.014; 
   vector<Iso> result;
     rcsv(fi, result );  int lres=result.size();

 double Nc[ntime]; for(int i=0;i<ntime;i++) Nc[i]=(double)i+1.;//cells number
 
for(int i=1;i<ntime;i++) mu += log(Nc[i]/Nc[0])/tex[i];
 mu /= ((double)(ntime-1)*1.0);
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
	 
   for(int j=0;j<lres;j++){
       string* s=result[j].getname();
       int niso=result[j].getniso();
       double* mid=result[j].getmid();
       double* std=result[j].getsd();
    int k(0);  for(int i=0;i<lmet;i++){
      int imatch=(*s).find(met[i]->getdescr());
        if(imatch != (int)std::string::npos){ met[i]->setex0();
         cout<<i<<": "<<met[i]->getdescr()<<" niso="<<niso<<" m0="<<mid[0]<<" sd0="<<std[0]<<endl;
         met[i]->sex(niso,mid,std,1); k++; break; }
                             }
   if(k==0){cout<<"k="<<k<<" no metabolite match?!?!?!\n";}
   }
//    while(!fi.eof()){
//      fi>>aaa; if(fi.eof()) break; cout<<"descr: "<<aaa<<endl;
//   int k=0; for(int i=0;i<lmet;i++){
//      int imatch=aaa.find(met[i]->getdescr());
//        if(imatch != (int)std::string::npos){met[i]->setex0();
//         cout<<i<<": "<<met[i]->getdescr()<<endl;
//          met[i]->read(fi,1); k++; break;}
//                             }
//           
//    }
     gl.setrav(0,1);
   for(int i=0;i<lres;i++) {result[i].delmid(); delete[] result[i].sd;} result.clear();
  return ts1;}

vector<string> Ldistr::spli(stringstream& test,char a){
    vector<string> seglist;    string segment;
  while(getline(test, segment, a)) seglist.push_back(segment); return seglist;
}

const int trac=0, lab=trac+1, abund=lab+1, inj=abund+1, repl=inj+1, etime=repl+1, emet=etime+1, efrg=emet+1, formula=efrg+1, intens=formula+1, isotopol=intens+1, conc=isotopol+1, nepar=conc+1;

void Ldistr::defcol(int nucol[],vector<string> vstr){
     int len=vstr.size(); string a; int pos;
  for(int i=0;i<len;i++){
     if(vstr[i].find("tracer")!=string::npos) nucol[trac]=i;
     else if(vstr[i].find("label")!=string::npos) nucol[lab]=i;
     else if(((pos=vstr[i].find("abundance"))!=string::npos)&&(pos<3)) nucol[abund]=i;
     else if(vstr[i].find("inject")!=string::npos) nucol[inj]=i;
     else if(vstr[i].find("Replicate")!=string::npos) nucol[repl]=i;
     else if(vstr[i].find("time")!=string::npos) nucol[etime]=i;
     else if(vstr[i].find("Metabolite")!=string::npos) nucol[emet]=i;
     else if(vstr[i].find("atomic")!=string::npos) nucol[efrg]=i;
     else if(vstr[i].find("formula")!=string::npos) nucol[formula]=i;
     else if(vstr[i].find("intens")!=string::npos) nucol[intens]=i;
     else if((vstr[i].find("isotopol")!=string::npos)&&(vstr[i].length()<15)) nucol[isotopol]=i;
     else if(vstr[i].find("concentration")!=string::npos) nucol[conc]=i;
  }
}

void Ldistr::rcsv(ifstream& fi,vector<Iso>& result ){
   string aaa;
   int cols[nepar];
    vector<string> titl, strok;  result.clear();
//   ifstream fi("midcorout.csv");
//     ifstream fi("exchangeformat.csv");
// List of column Titles
    getline(fi,aaa);
     stringstream test(aaa); 
      titl=spli(test,',');
   defcol(cols,titl);

//get rows
         int nstrok(0);
    while(getline(fi,aaa)) {strok.push_back(aaa); nstrok++;}
  
//  string Matrix of data for the analysis; segline[nstrok][columns #]
 vector<string> segline[nstrok]; int len;
   for(int j=0;j<nstrok;j++){

// split rows by "
   test.clear(); test.str(strok[j]);
      segline[j]=spli(test,'"'); len=segline[j].size();

// split rows by ,      
   if(len==1) {test.clear(); test.str(strok[j]);
          segline[j].clear(); segline[j]=spli(test,',');  len=segline[j].size();}//}

//erase wrong columns   
   for(int i=0;i<len;i++)
        if(segline[j][i]==",") {segline[j].erase(segline[j].begin()+i); len--; i--;}
   if(segline[j][0].length()==0) {segline[j].erase(segline[j].begin()); len--;}

//split joined columns       
   for(int i=0;i<len;i++)
    if((segline[j][i].length()>1)&&(segline[j][i].at(0)==',')) {
       stringstream ttt(segline[j][i]);     segline[j].erase(segline[j].begin()+i);
         vector<string> vstav; vstav=spli(ttt,',');
         if(vstav[0].length()<2) vstav.erase(vstav.begin());
         segline[j].insert(segline[j].begin()+i,vstav.begin(),vstav.end());
          }
          
   }
    cout<<"titles - data: "<< titl.size()<<" - "<<segline[5].size()<<"\n";
    
    tex.clear(); tex.push_back(0.);  
  for(int i=0;i<nstrok;i++) if(segline[0][cols[etime]]!="") {
   string s=segline[i][cols[etime]]; int j(0); ntime=tex.size(); double d=strtod(s.c_str(), NULL)*60.;
    for(j=0;j<ntime;j++) if(d==tex[j]) break; if(j==ntime) tex.push_back(d); }
     ntime=tex.size();  cout<<"incubation times (h): "; for(int i=0;i<ntime;i++) cout<<tex[i]/60.<<" "; cout<<"\n";
   
//convert MID from string to double
   double dis[nstrok];
   for(int i=0;i<nstrok;i++) {dis[i]=0.; if(segline[i][cols[conc]].length()>2) 
     dis[i]=strtod(segline[i][cols[conc]].c_str(), NULL);}
         
    int iro=0, lef=segline[0][cols[efrg]].at(4)-segline[0][cols[efrg]].at(1)+2;
     aaa=segline[0][cols[emet]]; 
     string chinj=segline[0][cols[inj]], sabun=segline[0][cols[abund]];
   
// orden mid for each injection
   vector<Iso> liso;
   
while(iro<(nstrok-1)){

   while((segline[iro][cols[emet]]==aaa)&&(segline[iro][cols[abund]]==sabun)){
//   for(int j=0;j<5;j++){
    Iso iso(lef); iso.setname(segline[iro][cols[emet]]);  int iiso(0); 
    while(segline[iro][cols[inj]]==chinj){
     if((segline[iro][cols[conc]].length()>2)&&(segline[iro][cols[abund]]==sabun)){
        iso.setmid(iiso,dis[iro]); iiso++; }
         iro++; }
    chinj=segline[iro][cols[inj]]; liso.push_back(iso); 
   }
   
   int lisi(liso.size());// for(int i=0;i<lisi;i++) liso[i].showmid();
// statistics grouping all injection
    Iso iso(lef); iso.setname(segline[iro-2][cols[emet]]);
    iso.calmesd(liso); iso.showmid("_mean");  iso.showsd("_sd"); result.push_back(iso);
   for(int i=0;i<liso.size();i++) liso[i].delmid(); liso.clear();
   
   while(segline[iro][cols[abund]]=="") if(iro<nstrok) iro++;
   if(iro<nstrok){
     lef=segline[iro][cols[efrg]].at(4)-segline[iro][cols[efrg]].at(1)+2;
     aaa=segline[iro][cols[emet]]; 
     chinj=segline[iro][cols[inj]];
  }  }   }

