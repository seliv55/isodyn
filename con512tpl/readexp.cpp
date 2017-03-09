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
     Tracer l13c=rcsv(fi, result );  int lres=result.size();
     cout<<" prueba: "<<l13c.fract<<endl;

 double Nc[ntime]; for(int i=0;i<ntime;i++) Nc[i]=(double)i+1.;//cells number
 
for(int i=1;i<ntime;i++) mu += log(Nc[i]/Nc[0])/tex[i];
 mu /= ((double)(ntime-1)*1.0);
	 xx[ngl]=10.; gl.setconc(xx[ngl]); gl.setex0();
	 xx[nglycog]=3.; glycog.setconc(xx[nglycog]);// glycog.setex0();
	 xx[nlac]=0.001; lac.setconc(xx[nlac]);
	 xx[nglu]=4.1; glu.setconc(xx[nglu]); glu25.setconc(xx[nglu]);
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
	 
    int k(0);
       string* s=l13c.getname(); cout<<"trac name="<<*s<<endl;
       int niso=l13c.getniso();
       data* mid=l13c.getmid();
      for(int i=0;i<lmet;i++){
      int imatch=(*s).find(met[i]->getdescr());
        if(imatch != (int)std::string::npos){ met[i]->setex0();
         cout<<i<<": "<<met[i]->getdescr()<<" niso="<<niso<<" m0="<<mid[0].mean<<" sd0="<<mid[0].sd<<endl;
         met[i]->sex(niso,mid,0); met[i]->setrav(1,0); k++; break; }
                             }
	 
   for(int j=0;j<lres;j++){
        s=result[j].getname();
        niso=result[j].getniso();
        mid=result[j].getmid();
    k=0;
      for(int i=0;i<lmet;i++){
      int imatch=(*s).find(met[i]->getdescr());
        if(imatch != (int)std::string::npos){ met[i]->setex0();
         cout<<i<<": "<<met[i]->getdescr()<<" niso="<<niso<<" m0="<<result[j].getmid()[0].mean<<" sd0="<<result[j].getmid()[0].sd<<endl;
         met[i]->sex(niso,result[j].getmid(),1); k++; break; }
                             }
   if(k==0){cout<<*s<<" k="<<k<<" no metabolite match?!?!?!\n";}
   }
    
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

int Ldistr::c13pos(string& s,int& nc,int& nlab){
  stringstream ttt(s); int ibin(0),ilab; nlab=0;
      vector<string> isopos; isopos=spli(ttt,','); nc=isopos.size();
       for(int i=0;i<nc;i++){ ilab=((int)isopos[i].at(0)-(int)'0');
        nlab+=ilab; ibin+=(ilab<<(nc-1-i));
       }
  return ibin; }
      
Tracer Ldistr::rcsv(ifstream& fi,vector<Iso>& result ){
   string aaa;
   int cols[nepar];
    vector<string> titl, strok;  result.clear();
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
// set tracer
  int nc, ibin,nlab;
  if(segline[0][cols[trac]]!="")  ibin=c13pos(segline[0][cols[lab]],nc,nlab);
       cout<<"c13pos: "; cout<<ibin<<" "<<nc<<" "<<nlab<<"\n";
   Tracer labmet(nc,segline[0][cols[trac]],ibin,strtod(segline[0][cols[abund]].c_str(), NULL),nlab);
          labmet.setrac(); cout<<"trac name="<<*labmet.getname()<<endl;
//  tex: time points of measures during incubation
    tex.clear(); tex.push_back(0.);
      for(int i=0;i<nstrok;i++) if(segline[i][cols[trac]]!="") { int j;
   double d=strtod(segline[i][cols[etime]].c_str(), NULL)*60.;
    for(j=0;j<tex.size();j++) if(d==tex[j]) break;
     if(j==tex.size()) tex.push_back(d); }     ntime=tex.size();
       cout<<"incubation times (h): "; for(int i=0;i<ntime;i++) cout<<tex[i]/60.<<" "; cout<<"\n";
   
//convert MID from string to double in whole column
   double dis[nstrok];
   for(int i=0;i<nstrok;i++) {dis[i]=0.; if(segline[i][cols[conc]].length()>2) 
     dis[i]=strtod(segline[i][cols[conc]].c_str(), NULL);}
         
    int iro=0, // row #
     lef=segline[0][cols[efrg]].at(4)-segline[0][cols[efrg]].at(1)+2; // fragment length
     aaa=segline[0][cols[emet]]; //metabolite name
    string chinj=segline[0][cols[inj]],//injection #
     strac=segline[0][cols[trac]]; //13C tracer
   
   vector<Iso> liso;
   
while(iro<(nstrok-1)){
// ordenate mid for each injection
   while((segline[iro][cols[emet]]==aaa)&&(segline[iro][cols[trac]]==strac)){
    Iso iso(lef,segline[iro][cols[emet]]);  int iiso(0); 
    while(segline[iro][cols[inj]]==chinj){
     if((segline[iro][cols[conc]].length()>2)){
        iso.setmid(iiso,dis[iro]); iiso++; }
         iro++; }
    chinj=segline[iro][cols[inj]]; liso.push_back(iso); 
   }
   
// statistics grouping all injection for each metabolite
    Iso iso(lef,segline[iro-2][cols[emet]]);
    iso.calmesd(liso); iso.showmid("_mean");  iso.showsd("_sd"); result.push_back(iso);
   for(int i=0;i<liso.size();i++) liso[i].delmid(); liso.clear();
   
   while(segline[iro][cols[trac]]=="") if(iro<nstrok) iro++;
   if(iro<nstrok){
     lef=segline[iro][cols[efrg]].at(4)-segline[iro][cols[efrg]].at(1)+2;
     aaa=segline[iro][cols[emet]]; 
     chinj=segline[iro][cols[inj]];
  }  }  
  return labmet; }

