#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "modlab.h"
using namespace std;

int Ldistr::getN() {
met[0]->ny=numx;
for(int i=1;i<lmet;i++) met[i]->ny=met[i-1]->ny+met[i-1]->getlen();
  glu25.ny=glu.ny;
metb[0]->ny= met[lmet-1]->ny+met[lmet-1]->getlen();
for(int i=1;i<lmetb;i++) metb[i]->ny=metb[i-1]->ny+metb[i-1]->getlen();
metk[0]->ny= metb[lmetb-1]->ny+metb[lmetb-1]->getlen();
for(int i=1;i<lmetk;i++) metk[i]->ny=metk[i-1]->ny+metk[i-1]->getlen();
return (Nn=metk[lmetk-1]->ny+metk[lmetk-1]->getlen());
}

void Ldistr::setdiso(double *pyinit) {
for(int i=0;i<=lmet;i++) met[i]->diso= &pyinit[met[i]->ny]; 
for(int i=0;i<lmetb;i++) metb[i]->diso= &pyinit[metb[i]->ny];
for(int i=0;i<lmetk;i++) metk[i]->diso= &pyinit[metk[i]->ny];
}

void Ldistr::setiso(double *pyinit) {
for(int i=0;i<=(lmet);i++) met[i]->iso= &pyinit[met[i]->ny]; 
for(int i=0;i<(lmetb);i++) metb[i]->iso= &pyinit[metb[i]->ny];
for(int i=0;i<(lmetk);i++) metk[i]->iso= &pyinit[metk[i]->ny];
}

void Ldistr::sklad(int itime){
 for(int i=0;i<expm0.size();i++) expm0[i]->skladm0(itime);
 for(int i=0;i<expcon.size();i++) expcon[i]->skladc(itime);
 }

void Ldistr::massfr() {
for(int i=0;i<expm0.size();i++) expm0[i]->percent(); glu.percent(1,1); glu25.percent(1,0);
}

double Ldistr::xits(int its) {int itp=its-1; double xi=0;
for(int i=0;i<expm0.size();i++) xi += expm0[i]-> chisq(its,expm0[i]->getmi());
return xi;}

void Ldistr::wriconex(ostringstream& fo) {
//t=1;  glc=t+1; glsd=glc+1;  lacc=glsd+1; lacsd=lacc+1;
  data *b;
 fo<<"time "; for(int i=0;i<expcon.size();i++) fo<<expcon[i]->getdescr()<<" sd "; fo<<endl;
 for(int its=0;its<ntime;its++) { fo<<tex[its];
   for(int i=0;i<expcon.size();i++){ b=expcon[i]->getconc();   fo<<" "<<b[its].mean<<" "<<b[its].sd; }
      fo<<"\n";         }
}

void Ldistr::wrim0ex(ostringstream& fo) {
//t=1;  lac=t+1; lacsd=lac+1; glu=lacsd+1; glusd=glu+1; ala=glusd+1; alasd=ala=+1; gly=alasd+1; glysd=gly+1; ser=glysd+1; sersd=ser+1; pro=sersd+1; prosd=pro+1; cit=prosd+1; citsd=cit+1; agl=citsd+1; aglsd=agl+1; asp=aglsd+1; aspsd=asp+1; fum=aspsd+1; fumsd=fum+1; mal=fumsd+1; malsd=mal+1; coa=malsd+1; coasd=coa+1;
  data *b;
fo<<"time "; for(int i=0;i<expm0.size();i++) fo<<expm0[i]->getdescr()<<" sd "; fo<<endl;
   for(int its=0;its<ntime;its++) {  fo<<tex[its];
       for(int i=0;i<expm0.size();i++){ b=expm0[i]->getexper(its);   fo<<" "<<b[0].mean<<" "<<b[0].sd; }
      fo<<"\n";         }
       }
void Ldistr::show(ostringstream& fo,double xfin) { fo<<setw(5)<<xfin;
for(int i=0;i<expm0.size();i++) expm0[i]-> showm0(fo);}

double Ldistr::consum() { double sum(0.);
 for(int i=0;i<=lmet;i++) { int j;
   for(j=0;j<expcon.size();j++) if(met[i]->getdescr()==expcon[j]->getdescr()) break;
     if(j==expcon.size()) {sum += met[i]->sumt(); cout<<met[i]->getdescr()<<"  ";}}
 for(int i=0;i<lmetb;i++) {sum += metb[i]->sumt(); cout<<metb[i]->getdescr()<<"  ";}
 for(int i=0;i<lmetk;i++) {sum += metk[i]->sumt(); cout<<metk[i]->getdescr()<<"  ";}cout<<endl;
return sum;}


