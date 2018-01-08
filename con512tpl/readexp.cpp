#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
using namespace std;
double Vi, xribi, xasp=1.,mu;

double Ldistr::readExp (char fn[],int ntr) {
  ifstream fi(fn); double Ti,ts1;  mu=0.;  Vi=0.014; 
     Tracer l13c=rcsv(fi, result,ntr );  int lres=result.size();
     cout<<"readExp: 13C-substrate: "<<*l13c.getname()<<", fraction: "<<l13c.fract<<endl;

 double Nc[ntime]; for(int i=0;i<ntime;i++) Nc[i]=(double)i+1.;//**??** cells number
 
for(int i=1;i<ntime;i++) mu += log(Nc[i]/Nc[0])/tex[i];
 mu /= ((double)(ntime-1)*1.0);
           setcon();
	  glu.setconc(xx[nglu]); glu25.setconc(xx[nglu]);

	 
       l13c.setmid(markis,marfrac);  l13c.setmid(0,(1-marfrac)); 
        itrac=findmet(l13c); cout<<"readExp: tracer="<<met[itrac]->getdescr()<<'\n';
	 
   for(int j=0;j<lres;j++){ findmet(result[j]);   }
//	 gl.setex0(); gln.setex0();
//    for(int j=0;j<expm0.size();j++) cout<<expm0[j]->getdescr()<<endl;
  return ts1;}

vector<string> Ldistr::spli(stringstream& test,char a){
    vector<string> seglist;    string segment;
  while(getline(test, segment, a)) seglist.push_back(segment); return seglist;
}

const int trac=0, lab=trac+1, abund=lab+1, injec=abund+1, etime=injec+1, emet=etime+1, efrg=emet+1, formula=efrg+1, intens=formula+1, isotopol=intens+1, conc=isotopol+1, nepar=conc+1;

void Ldistr::defcol(int nucol[],vector<string> vstr){
     int len=vstr.size(); string a; int pos;
  for(int i=0;i<len;i++){
     if(vstr[i].find("tracer")+1) nucol[trac]=i;
     else if(vstr[i].find("label")+1) nucol[lab]=i;
     else if(((pos=vstr[i].find("abundance"))+1)&&(pos<3)) nucol[abund]=i;
     else if(vstr[i].find("injection")+1) nucol[injec]=i;
     else if(vstr[i].find("time")+1) nucol[etime]=i;
     else if(vstr[i].find("Metabolite")+1) nucol[emet]=i;
     else if(vstr[i].find("atomic positions")+1) nucol[efrg]=i;
     else if(vstr[i].find("formula")+1) nucol[formula]=i;
     else if(vstr[i].find("intens")+1) nucol[intens]=i;
     else if(vstr[i]=="\"isotopologue\"") {nucol[isotopol]=i;} //cout<<"isotopol: "<<nucol[isotopol]<<endl;
     else if(vstr[i].find("gue abund")+1) nucol[conc]=i; 
  }
}

int Ldistr::findmet(Iso& iso) { int k(-1); 
//find met[i] for which there are experimental data presented in the unit iso
     int j(0); while(iso.gett()>(tex[j]+1e-7)) j++;
   for(int i=0;i<lmet;i++)  if(iso.getname()->find(met[i]->getdescr())+1){ k=i;
     for(int iex=0;iex<expm0.size();iex++) if(met[i]->getdescr()==expm0[iex]->getdescr()){k=-2; break;}//repetitions
          if(k>=0) {met[i]->setexper(ntime);  expm0.push_back(met[i]);} //if met[i] was not present in expm0, add met[i]
       met[i]->sex(iso.getniso(),iso.getmid(),j); //set experimental mid for a given met & time
           cout<<"findmet "<<met[i]->getdescr()<<" m0="<<iso.getmid()[0].mean<<" t="<<iso.gett()<<endl; break;}
           
     if(k==-1) for(int i=0;i<lmetk;i++)  if(iso.getname()->find(metk[i]->getdescr())+1){ k=i;
     for(int iex=0;iex<kexpm0.size();iex++) if(metk[i]->getdescr()==kexpm0[iex]->getdescr()){k=-2; break;}//repetitions
         if(k>=0) {metk[i]->setexper(ntime);  kexpm0.push_back(metk[i]);}//if met[i] was not present in kexpm0, add met[i]
       metk[i]->sex(iso.getniso(),iso.getmid(),j); //set experimental mid for a given met & time
           cout<<"findmet "<<metk[i]->getdescr()<<" m0="<<iso.getmid()[0].mean<<" t="<<iso.gett()<<endl; break;}
     if(k==-1)cout<<"findmet "<<(*iso.getname())<<" no metabolite match?!?!?!\n";
     return k;}

int Ldistr::c13pos(string& s,int& nc,int& nlab){
  stringstream ttt(s); int ibin(0),ilab; nlab=0;
      vector<string> isopos; isopos=spli(ttt,','); nc=isopos.size();
       for(int i=0;i<nc;i++){ ilab=((int)isopos[i].at(0)-(int)'0');
        nlab+=ilab; ibin+=(ilab<<(nc-1-i));
       }
  return ibin; }
      
set<string> Ldistr::findopt(string a, vector<string> strok){
         int nstrok(strok.size());
   set<string> metka;               // set labeled substrates
    for(int i=0;i<nstrok;i++){ size_t pos=strok[i].find("C13]-");
     if(pos+1){ metka.insert(strok[i].substr(pos,17)); }       }
    for(set<string>::iterator it=metka.begin(); it!=metka.end(); it++) cout<<"findopt "<<*it<<'\n';
  return metka;
}
 void Ldistr::splitstrings(vector<string> segline[],int nstrok,vector<string> substrok){
   int len;
   for(int j=0;j<nstrok;j++){

// split rows by "
  stringstream test; test.str(substrok[j]); 
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
 }
    
Tracer Ldistr::rcsv(ifstream& fi,vector<Iso>& result,int mar ){
   string aaa;
   int cols[nepar];
    vector<string> titl, strok;//  result.clear();
//** 1. CONVERSION OF DATA FROM STRINGS
// List of column Titles
    getline(fi,aaa);           // read first line
     stringstream test(aaa); 
      titl=spli(test,',');    // names and positions of columns
   defcol(cols,titl);
//get rows
    while(getline(fi,aaa)) strok.push_back(aaa);
    int nstrok=strok.size();
    
     set<string> metka=findopt("C13]-",strok);
    
//                              chosing strings corresponding to the labeled substrate
    vector<string> substrok;//  result.clear();
     set<string>::iterator itm=metka.begin(); for(int i=0;i<mar;i++) itm++;
    for(int i=0;i<nstrok;i++) if(strok[i].find(*itm)+1) substrok.push_back(strok[i]);
    nstrok=substrok.size(); cout<<"rcsv: "<<nstrok<<" strok\n";
    
//  string Matrix of data for the analysis; segline[nstrok][columns #]
  vector<string> segline[nstrok];
  splitstrings(segline, nstrok,substrok);
    cout<<"rcsv: titles - data: "<< titl.size()<<" - "<<segline[5].size()<<"\n";
//** 2. GETTING DATA IN APPROPRIATE FORMAT

// set tracer
  int nc, nlab;
  if(segline[0][cols[trac]]!="") { markis=c13pos(segline[0][cols[lab]],nc,nlab);
     marfrac=stod(segline[0][cols[abund]])/100.;
      }
   Tracer labmet(nc,0.,segline[0][cols[trac]],markis,marfrac,nlab);
          labmet.setrac(); 
   cout<<"rcsv: Tracer="<<*labmet.getname()<<"; Isotopomer="<<markis<<"; Carbons="<<nc;
   cout<<"; Labels="<<nlab<<"; Fraction="<<marfrac<<'\n';
          
//  tex: time points of measures during incubation
    tex.clear();
     set<string> stex; stex.insert("0");
    for(int i=0;i<nstrok;i++) stex.insert(segline[i][cols[etime]]); int itex=0;
    cout<<"rcsv: tex="; 
   for(set<string>::iterator it=stex.begin();it!=stex.end();it++){
     tex.push_back(stod(*it)*60);cout<<tex[itex]<<" ";++itex;} cout << '\n';
          ntime=tex.size();
     
// Take values of some parameters
 int lef=(int)segline[0][cols[efrg]].at(4)-(int)segline[0][cols[efrg]].at(1)+2; // number of isotopomers in fragment
    string metnm=segline[0][cols[emet]], //metabolite name
           chinj=segline[0][cols[isotopol]];//injection #
//           strac=segline[0][cols[trac]]; //13C tracer
//convert MID from string to double in whole column
    double dis[nstrok];
      for(int i=0;i<nstrok;i++){ dis[i]=0.;
       if(segline[i][cols[conc]].length()>2)  dis[i]=stod(segline[i][cols[conc]])*0.01;
     }
   
   vector<Iso> liso,liso1;  Iso *iso;  int iiso(0),iro(0);
while(iro<(nstrok)){ cout<<"rcsv: name="<<metnm<<" cinj="<<chinj<<" "<<segline[iro][cols[isotopol]]<<'\n';
   while(segline[iro][cols[emet]]==metnm){ // take data for metabolite metnm
         if(iiso==0) {iso=new Iso(lef,stod(segline[iro][cols[etime]])*60.,segline[iro][cols[emet]]);
          if(segline[iro][cols[isotopol]]==chinj) iro++;}
    while((segline[iro][cols[isotopol]]!=chinj)&&(iiso<lef)){  // ordenate mid for each injection
                   iso->setmid(iiso,dis[iro]); iiso++; iro++;}
                   if(iiso==lef) { liso.push_back(*iso);}
                   
                   
         if(segline[iro][cols[isotopol]]!=chinj) {iro++; iiso++;}
         else { iiso=0; iro++;}
    if(iro>=nstrok)  {liso.push_back(*iso); break;}
   }
//   liso.push_back(labmet);
// statistics grouping all injection for each metabolite
     for(int j=0;j<tex.size();j++){
   for(int i=0;i<liso.size();i++) if(tex[j]==liso[i].gett()) liso1.push_back(liso[i]);
//     for(int i=0;i<liso1.size();i++) liso1[i].showmid();
    if(liso1.size()){ iso=new Iso(lef,tex[j],segline[iro-2][cols[emet]]);
    iso->calmesd(liso1);
     iso->showmid("_mean");  iso->showsd("_sd");
     result.push_back(*iso); 
    liso1.clear();    
//    correct labeled fraction
    if(iso->gett()<1e-7){
    if((*labmet.getname()).find(*iso->getname())+1){
     cout<<"**label in "<<*iso->getname()<<" t="<<iso->gett()<<endl;
    marfrac=iso->getmid()[nlab].mean;    }
                 }
        delete iso; cout<<'\n';   }
     }
         liso.clear();
   if(iro<(nstrok-1)){ iiso=0;
     lef=(int)segline[iro][cols[efrg]].at(4)-(int)segline[iro][cols[efrg]].at(1)+2;
     metnm=segline[iro][cols[emet]];  }
         }
  return labmet; }

