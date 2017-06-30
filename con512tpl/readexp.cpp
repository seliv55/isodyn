#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <unordered_set>
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
using namespace std;
double Vi, xribi, xasp=1.,mu;

double Ldistr::readExp (char fn[]) {
  ifstream fi(fn); double Ti,ts1;  mu=0.; dt=0.6;  Vi=0.014; 
     Tracer l13c=rcsv(fi, result );  int lres=result.size();
     cout<<" 13C%: "<<*l13c.getname()<<", fraction: "<<l13c.fract<<endl;

 double Nc[ntime]; for(int i=0;i<ntime;i++) Nc[i]=(double)i+1.;//cells number
 
for(int i=1;i<ntime;i++) mu += log(Nc[i]/Nc[0])/tex[i];
 mu /= ((double)(ntime-1)*1.0);
           for(int i=0;i<nfbp;i++)  met[i]->setconc(xx[i]);
	 gl.setex0(); gln.setex0();
	  glu.setconc(xx[nglu]); glu25.setconc(xx[nglu]);
	  

	 xx[nagl]=1.0;
	 xx[ncthf]=0.5;	 lmet=nfbp;//sizeof(met)/sizeof(*met); 
	 
       l13c.setmid(markis,marfrac*100.);  l13c.setmid(0,100*(1-marfrac)); 
        itrac=findmet(l13c.getname(),l13c.getniso(),l13c.getmid());
	 
   for(int j=0;j<lres;j++){ findmet(result[j].getname(),result[j].getniso(),result[j].getmid());   }
    
  return ts1;}

vector<string> Ldistr::spli(stringstream& test,char a){
    vector<string> seglist;    string segment;
  while(getline(test, segment, a)) seglist.push_back(segment); return seglist;
}

const int trac=0, lab=trac+1, abund=lab+1, injec=abund+1, etime=injec+1, emet=etime+1, efrg=emet+1, formula=efrg+1, intens=formula+1, isotopol=intens+1, conc=isotopol+1, nepar=conc+1;

void Ldistr::defcol(int nucol[],vector<string> vstr){
     int len=vstr.size(); string a; int pos;
  for(int i=0;i<len;i++){
     if(vstr[i].find("tracer")!=string::npos) nucol[trac]=i;
     else if(vstr[i].find("label")!=string::npos) nucol[lab]=i;
     else if(((pos=vstr[i].find("abundance"))!=string::npos)&&(pos<3)) nucol[abund]=i;
     else if(vstr[i].find("injection")!=string::npos) nucol[injec]=i;
     else if(vstr[i].find("time")!=string::npos) nucol[etime]=i;
     else if(vstr[i].find("Metabolite")!=string::npos) nucol[emet]=i;
     else if(vstr[i].find("atomic positions")!=string::npos) nucol[efrg]=i;
     else if(vstr[i].find("formula")!=string::npos) nucol[formula]=i;
     else if(vstr[i].find("intens")!=string::npos) nucol[intens]=i;
     else if(vstr[i]=="\"isotopologue\"") {nucol[isotopol]=i;} //cout<<"isotopol: "<<nucol[isotopol]<<endl;} 
     else if(vstr[i].find("gue abund")!=string::npos) nucol[conc]=i; 
  }
}

int Ldistr::findmet(string *s,int niso,data* mid) { int k(0);
      for(int i=0;i<=lmet;i++){
      size_t imatch=(*s).find(met[i]->getdescr()); 
        if(imatch +1){ met[i]->setex0(); k=i;
         cout<<i<<": "<<met[i]->getdescr()<<" niso="<<niso<<" m0="<<mid[0].mean<<" sd0="<<mid[0].sd<<endl;
         met[i]->sex(niso,mid,1); expm0.push_back(met[i]); break;}
                             }
   if(k==0){cout<<*s<<" k="<<k<<" no metabolite match?!?!?!\n";}
     return k;}

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
    vector<string> titl, strok, substrok;//  result.clear();
//** 1. CONVERSION OF DATA FROM STRINGS
// List of column Titles
    getline(fi,aaa);
     stringstream test(aaa); 
      titl=spli(test,',');
   defcol(cols,titl);
//get rows
    while(getline(fi,aaa)) strok.push_back(aaa);
    
         int nstrok(strok.size());
   unordered_set<string> metka;
    for(int i=0;i<nstrok;i++){ size_t pos=strok[i].find("C13]-"); if(pos!=string::npos){
     string bbb=strok[i].substr(pos,10); metka.emplace(bbb); }       }
    for(const string& it: metka) cout<<it<<"\n";
    
    for(int i=0;i<nstrok;i++) if(strok[i].find(*(metka.begin()))!=string::npos) substrok.push_back(strok[i]);
    nstrok=substrok.size();
    
//  string Matrix of data for the analysis; segline[nstrok][columns #]
 vector<string> segline[nstrok]; int len;
   for(int j=0;j<nstrok;j++){

// split rows by "
   test.clear(); test.str(substrok[j]);
      segline[j]=spli(test,'"'); len=segline[j].size();

// split rows by ,      
   if(len==1) {test.clear(); test.str(substrok[j]);
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
//** 2. GETTING DATA IN APPROPRIATE FORMAT

// set tracer
  int nc, nlab;
  if(segline[0][cols[trac]]!="") { markis=c13pos(segline[0][cols[lab]],nc,nlab);
     marfrac=strtod(segline[0][cols[abund]].c_str(), NULL)/100.;
      }
       cout<<"Tracer="<<segline[0][cols[trac]]<<"; Isotopomers="<<markis<<"; Carbons="<<nc<<"; Labels="<<nlab<<"; Fraction="<<marfrac<<"\n";
   Tracer labmet(nc,0.,segline[0][cols[trac]],markis,marfrac,nlab);
          labmet.setrac();
//  tex: time points of measures during incubation
    tex.clear(); tex.push_back(0.);
   double dis[nstrok];
//convert MID from string to double in whole column
    int iro=0, // row #
     lef=(int)segline[0][cols[efrg]].at(4)-(int)segline[0][cols[efrg]].at(1)+2; // number of isotopomers in fragment
    string metnm=segline[0][cols[emet]], //metabolite name
           chinj=segline[iro][cols[isotopol]],//injection #
           strac=segline[0][cols[trac]]; //13C tracer
   cout<<"isotopologue="<<segline[iro][cols[isotopol]]<<endl;
   vector<Iso> liso;   int iiso(0); Iso *iso;
   
    
      for(int i=0;i<nstrok;i++){
       if(segline[i][cols[trac]]!="") { int j;
   double d=strtod(segline[i][cols[etime]].c_str(), NULL)*60.;
    for(j=0;j<tex.size();j++) if(d==tex[j]) break;
     if(j==tex.size()) tex.push_back(d);}
     
      dis[i]=0.; if(segline[i][cols[conc]].length()>2) 
     dis[i]=strtod(segline[i][cols[conc]].c_str(), NULL);
     }
          ntime=tex.size();
       cout<<"incubation times (h): "; for(int i=0;i<ntime;i++) cout<<tex[i]/60.<<" "; cout<<"\n";
   
while(iro<(nstrok-1)){
   while(segline[iro][cols[emet]]==metnm){ // extract mid for each metabolite from all injections
         if(iiso==0) {iso=new Iso(lef,strtod(segline[iro][cols[etime]].c_str(), NULL)*60.,segline[iro][cols[emet]]);
          if(segline[iro][cols[isotopol]]==chinj) iro++;}
    while((segline[iro][cols[isotopol]]!=chinj)&&(iiso<lef)){  // ordenate mid for each injection
                   iso->setmid(iiso,dis[iro]); iiso++; iro++;}
                   if(iiso==lef){ liso.push_back(*iso);}
                   
         if(segline[iro][cols[isotopol]]!=chinj) {iro++; iiso++;}
         else { iiso=0; iro++;}
    if(iro>=nstrok)  {liso.push_back(*iso); break;}
   }
// statistics grouping all injection for each metabolite
     for(int i=0;i<liso.size();i++) liso[i].showmid();
    iso=new Iso(lef,0.,segline[iro-2][cols[emet]]);
    iso->calmesd(liso); iso->showmid("_mean");  iso->showsd("_sd"); cout<<"\n";
    if((*labmet.getname()).find(*iso->getname())+1) {marfrac=iso->getmid()[nlab].mean*0.01; }
     result.push_back(*iso);
    liso.clear();
//   while(segline[iro][cols[trac]]=="") if(iro<nstrok) iro++;
   if(iro<(nstrok-1)){ iiso=0;
     lef=(int)segline[iro][cols[efrg]].at(4)-(int)segline[iro][cols[efrg]].at(1)+2;
     metnm=segline[iro][cols[emet]];
  }    }
  return labmet; }

