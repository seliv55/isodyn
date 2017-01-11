#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

//---------------------------------------------------------------------------
using namespace std;

int main ( int argc, char *argv[] ) {
   ifstream fi("model");  int kfl, knkin, intra, iextra, kiso;
     string aaa, stri, siso;
//reding the scheme of the model:
getline(fi,aaa);
getline(fi,aaa);
getline(fi,aaa);
       fi>>aaa>>kfl; getline(fi,aaa); getline(fi,stri);
       
     string flname[kfl+18], subst[kfl], prod[kfl], eq[kfl], isos[kfl], isop[kfl], isom[kfl], revfl[kfl];
     int i=0;
   for(;;) {fi>>aaa; if(aaa[0]=='/') break; else {getline(fi,flname[i]); i++;}}
   ofstream fo("model1");
   fo<<"fluxes  "<<i<<"\n"<<stri<<"\n";
    for(int j=0;j<i;j++) fo<<j<<" "<<flname[j]<<"\n";
    return 0;}
   

