//---------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

//---------------------------------------------------------------------------
using namespace std;
string partk="e0tk", parta="e0ta", parald="e0ald", ssym="fum";
string setk1, seta1, setald1;

inline string mm(int &ii,string par[],string& rpar, string fln, string sub, string prd, bool ff=0){ 
   string sfl="flx["+fln+"]"; par[ii]="V"+fln; rpar+="rea["+fln+"].v="+par[ii]+"; ";
    string foo =sfl+"= ";
     if(ff) foo += "nv[D]*";
      foo +="MM(nv["+par[ii]+"], ";      ii++;
       par[ii]="K"+fln; rpar+="rea["+fln+"].km="+par[ii]+"; ";
        foo +="nv["+par[ii]+"], xx["+sub+"]);"; ii++;   while(foo.length()<57) foo +=" ";
         foo +="dydx["+sub+"] -= "+sfl+";";             while(foo.length()<88) foo +=" ";
          if(prd!="0") foo +="dydx["+prd+"] += "+sfl+";"; foo +="\n";
return foo;}

inline string con(int &ii,string par[],string& rpar, string fln, string prd){ 
   string sfl="flx["+fln+"]"; par[ii]="V"+fln;
     rpar+="rea["+fln+"].v="+par[ii]+"; ";
      string foo =sfl+"= nv["+par[ii]+"];";    ii++;  while(foo.length()<88) foo +=" ";
          foo +="dydx["+prd+"] += "+sfl+";"; foo +="\n";
return foo;}

inline string cs(int &ii,string par[], string csfl, string s1, string s2,string prd){ 
    string sfl="flx["+csfl+"]"; par[ii]="V"+csfl; 
      string foo =sfl+"= nv[D]*MM(nv["+par[ii]+"], ";      ii++;
       par[ii]="K"+s1.substr(1);
        foo +="nv["+par[ii]+"], xx["+s1+"])*MM(1.,nv["; ii++; 
       par[ii]="K"+s2.substr(1);
        foo +=par[ii]+"], xx["+s2+"]);"; ii++; 
         foo +="\n\t\t\t\tdydx["+s1+"] -= "+sfl+";";  
         foo +="\tdydx["+s2+"] -= "+sfl+";";   
           foo +="\tdydx["+prd+"] += "+sfl+";\n";
return foo;}

inline string ald(int &ii,string par[],string aldfl,string alds,string aldp){ 
    string sfl="flx["+aldfl+"]"; par[ii]=parald; ii++; par[ii]="Kald"; ii++;
      setald1 ="aldolase.st1fl(&"+sfl+", xx[n"+alds+"], xx[n"+aldp+"]);"; string foo=setald1;
      while(foo.length()<57) foo +=" ";  foo +="dydx[n"+alds+"] -= "+sfl+";";
      while(foo.length()<88) foo +=" ";  foo +="dydx[n"+aldp+"] += 2.*"+sfl+";\n";     
return foo;}

inline string ta(int &ii,string par[],string tafl,string tak1,string taa1,string taa2,string tak2){ 
    string sfl="flx["+tafl+"]"; par[ii]=parta; ii++;
    par[ii]="KdEa"; ii++; par[ii]="KdFa"; ii++; par[ii]="KdGa"; ii++; par[ii]="Kcfg"; ii++;
      seta1 ="ta.st1fl(&"+sfl+", xx[n"+tak1+"]/fh6, xx[n"+taa1+"]/ft3, xx[n"+taa2+"], xx[n"+tak2+"]);"; string foo=seta1;
      foo +="\n\t\t\t\t\tdydx[n"+tak1+"] -= "+sfl+";\t\t dydx[n"+taa1+"] += "+sfl+";";
      foo +="\n\t\t\t\t\tdydx[n"+taa2+"] -= "+sfl+";\t\t dydx[n"+tak2+"] += "+sfl+";\n";    
return foo;}

inline string tk(int &ii,string par[],string tkfl,string tka1,string tkk1,string tka2,string tkk2,string tka3,string tkk3){ 
    string sfl="flx["+tkfl+"]", sfl1="flx["+tkfl+"+1]", sfl2="flx["+tkfl+"+2]"; 
    par[ii]=partk; ii++; par[ii]="Kdx5"; ii++; par[ii]="Kdg3"; ii++;
    par[ii]="Kdr5"; ii++; par[ii]="Kdf6"; ii++; par[ii]="Kcef"; ii++; par[ii]="Kcfe"; ii++;
  setk1="tk.st1fl(&"+sfl+", xx[n"+tka1+"]/ft3, xx[n"+tkk1+"], xx[n"+tka2+"], xx[n"+tkk2+"]/fh6, xx[n"+tka3+"], xx[n"+tkk3+"]);";
     string foo=setk1+"\n\t\t\t\t\tdydx[n"+tkk1+"] -= "+sfl+";\t\tdydx[n"+tka1+"] += "+sfl+";\n"; 
      foo +="\t\t\t\t\tdydx[n"+tkk3+"] -= "+sfl1+";\t\tdydx[n"+tka3+"] += "+sfl1+";\n";    
      foo +="\t\t\t\t\tdydx[n"+tkk2+"] -= "+sfl2+";\t\tdydx[n"+tka2+"] += "+sfl2+";\n";    
return foo;}

inline string listint(string sname[],int inum,string sfin){
  string aaa="const int "+sname[0]+"=0, ";
   for(int i=1;i<inum;i++)   aaa+=sname[i]+"="+sname[i-1]+"+1, ";
     aaa+=sfin+"="+sname[inum-1]+"+1;\n\n";
 return aaa;}
inline string listname(string sname[],int inum,string sfin){
  string aaa="extern const int ";
   for(int i=0;i<inum;i++)   aaa+=sname[i]+", ";   aaa+=sfin+";\n\n";
 return aaa;}
  inline string setvals(string snames[],int nelem, string sval){ string aaa="";
    for(int i=0;i<nelem;i++) {ostringstream snum; snum<<i; aaa+=snum.str()+") "+snames[i]+"\t"+sval+"\n";}
      return aaa;}

int main ( int argc, char *argv[] ) {
   ifstream fi("model");  int kfl, knkin, intra, iextra, kiso;
     string aaa, stri, siso;
// variables of kinetic model   
       fi>>aaa>>knkin;   string newkin[knkin], rst[knkin],rr[knkin]; cout<<aaa<<endl;  getline(fi,aaa); getline(fi,aaa);
   for(int i=0;i<knkin;i++) {fi>> newkin[i]; rr[i]="r"+newkin[i]+"[]";  rst[i]="const int "+rr[i]+"={";}
//reding the scheme of the model:
       fi>>aaa>>kfl; getline(fi,aaa); getline(fi,aaa);
       
     string flname[kfl+18], subst[kfl], prod[kfl], eq[kfl], isos[kfl], isop[kfl], isom[kfl], revfl[kfl];
   for(int i=0;i<kfl;i++){ fi>>aaa>> flname[i]>> subst[i]>>prod[i]>>eq[i]>>isos[i]>>isop[i]>>isom[i]>>revfl[i];
     bool err=true; for(int j=0;j<knkin;j++) if(subst[i]==newkin[j]) {rst[j]+="-"+flname[i]+","; err=false;}
                             else if(prod[i]==newkin[j]) {rst[j]+=flname[i]+","; err=false;}
                             else if(newkin[j]=="ne4") err=false;
                              else if(newkin[j]=="ns7") err=false;
         if(err) throw flname[i]+": subst/prod not declared"; 
         } string sss="const int *Parray::react[numx]; React Parray::rea[nflx];\n";
         for(int j=0;j<knkin;j++) {rst[j]+="99};\n"; sss+=rst[j];} sss+="void Parray::setreact(){\n";
         for(int j=0;j<knkin;j++) {sss+="react["+newkin[j]+"]=r"+newkin[j]+";\n";}
         sss+="};\n";
   
     string csfl[2], css1[2], css2[2],csp[2],ics1[2],ics2[2],icsp[2],csm[2]; fi>>aaa;
   for(int i=0;i<2;i++) fi>>csfl[i]>>css1[i]>>css2[i]>>csp[i]>>ics1[i]>>ics2[i]>>icsp[i]>>csm[i];//citrate synthase
       cout<<"csfl[0]="<<csfl[0]<<endl;
     string aldfl, alds, aldp;
   fi>>aaa>>aldfl>>alds>>aldp;//aldolase
    cout<<"aldfl="<<aldfl<<endl;
     string tafl, tak1,taa1,taa2,tak2;
   fi>>aaa>>tafl>>tak1>>taa1>>taa2>>tak2;//transaldolase
   
     string tkfl, tka1,tkk1,tka2,tkk2,tka3,tkk3;
   fi>>aaa>>tkfl>>tka1>>tkk1>>tka2>>tkk2>>tka3>>tkk3;//transketoolase
   
          getline(fi,aaa); getline(fi,aaa);
       fi>>aaa>>iextra>>aaa>>intra; kiso=intra+iextra; int ncarb[kiso]; cout<<"kiso="<<kiso<<endl;
       
     string stype[kiso],  newiso[kiso], edata[kiso], xi[kiso];// read iso: type,name, edata(conc,lab,0), chi(yes,no)
   for(int i=0;i<kiso;i++) {fi>>aaa>>stype[i]>>newiso[i]>>edata[i]>>xi[i]; ncarb[i]=atoi(stype[i].substr(stype[i].size()-2,1).c_str());}
         fi.close();
         
//setting ODE model, conversion of fluxes for iso model (flfor):
       int ii(0);
  string par[2*kfl+17], fo="void Fit::f(const double *y,double *dydx) {\nfor(int i=0;i<numx;i++) dydx[i]=0.;\n";
  string flfor="void Parray::flfor(){\nfor(int i=0;i<nflx;i++) fluxes[i] = flx[i] * dt/Vi;\n";
  string rpar="void Parray::reapar(){\n";
   for(int i=0;i<kfl;i++) {
     if(eq[i]=="mm") {  fo += mm(ii,par,rpar,flname[i],subst[i],prod[i]);
                flfor +="fluxes["+flname[i]+"] /= xx["+ subst[i]+ "];\n"; }
      else if(eq[i]=="*mm") {  fo += mm(ii,par,rpar,flname[i],subst[i],prod[i],1);//flux multiplied by D
                flfor +="fluxes["+flname[i]+"] /= xx["+ subst[i]+ "];\n"; }
       else if(eq[i]=="con") { fo += con(ii,par,rpar,flname[i],prod[i]);}
       else if(eq[i]=="con1") { fo += con(ii,par,rpar,flname[i],prod[i]);
                 if(isos[i]!="0") flfor +="fluxes["+flname[i]+"]  /= horse.get"+ isos[i]+ "();\n";
                  else  flfor +="fluxes["+flname[i]+"]  /= horse.get"+ isop[i-1]+ "();\n"; }
       else if(eq[i]=="out") {par[ii]="V"+flname[i]; fo +="flx["+flname[i]+"]=nv["+par[ii]+"];\n";
         rpar+="rea["+flname[i]+"].v="+par[ii]+"; "; ii++;
          flfor +="fluxes["+flname[i]+"]  /= horse.get"+ isos[i]+ "();\n";  }
       else {par[ii]="V"+flname[i]; fo +="flx["+flname[i]+"]=nv["+par[ii]+"];\n";
         rpar+="rea["+flname[i]+"].v="+par[ii]+"; "; ii++;}
         }
           par[ii]="D"; ii++; 
           rpar+="};\n";
      fo += cs(ii,par,csfl[0],css1[0],css2[0],csp[0]); flfor +="fluxes["+csfl[0]+"] /= (xx["+ css1[0]+ "]*xx["+ css2[0]+ "]);\n";
                                          int ifl=kfl; flname[ifl]=csfl[0]; ifl++;
//      fo += cs(ii,par,csfl[1],css1[1],css2[1],csp[1]); flfor +="fluxes["+csfl[1]+"] /= (xx["+ css1[1]+ "]*xx["+ css2[1]+ "]);\n";
//                                          flname[ifl]=csfl[1]; ifl++;
      fo += ald(ii,par,aldfl, alds, aldp); flname[ifl]=aldfl; ifl++;  flname[ifl]="aldrev"; ifl++;
                                          flname[ifl]="aldfli"; ifl++;  flname[ifl]="aldi1"; ifl++;
                                          flfor +="fluxes["+aldfl+"] /= xx[n"+ alds+ "];\n";
                                          flfor +="fluxes["+aldfl+"+1] /= (xx[n"+ aldp+ "]*xx[n"+ aldp+ "]);\n";
                                          flfor +="fluxes["+aldfl+"+2] /= xx[n"+ alds+ "];\n";
      fo += ta(ii,par,tafl, tak1,taa1,taa2,tak2); flname[ifl]=tafl; ifl++; flname[ifl]="s7f6a";ifl++;
                                                   flname[ifl]="f6g3a"; ifl++; flname[ifl]="s7e4a"; ifl++;
      fo += tk(ii,par,tkfl, tka1,tkk1,tka2,tkk2,tka3,tkk3); flname[ifl]=tkfl; ifl++;
                 flname[ifl]="s7p5"; ifl++; flname[ifl]="f6p5"; ifl++; flname[ifl]="p5f6"; ifl++; flname[ifl]="f6s7"; ifl++;
                 flname[ifl]="s7f6"; ifl++; flname[ifl]="p5g3i"; ifl++; flname[ifl]="f6e4i"; ifl++; flname[ifl]="s7p5i"; ifl++;
               fo+="}\n\n"; flfor+="}\n\n";
               
//declaration of the list of fluxes and variables for "nums.hh"
  string hlfl=listname(flname,ifl,"nflx");   hlfl+=listname(newkin,knkin,"numx"); hlfl+=listname(rr,knkin,"nNV");
     hlfl+="extern const int Nn,e0ald, D, Vcs0;\n";
     hlfl+="extern double dt,dif,dif0, suxx, Vi,Vt, mu, xx[], fluxes[], flx[],ystart[];\n";
     hlfl+="extern double nv1[],nv2[],xinit1[],xinit2[];\nextern time_t ts,tf,tcal;\n";
     hlfl+="extern std::string foc, kin, fex1, fex2;\nextern int ifn;\n";
                                 ofstream fout("../include/nums.hh"); fout<<hlfl; fout.close();
                                 
 //file with parameters values
               aaa=setvals(par,ii, "0.01")+"1 3 5 7 9 11\n"+setvals(newkin,knkin, "0.01")+setvals(flname,ifl, "0.001");
 fout.open("parameters"); fout<<aaa; fout.close();
    
 aaa="#include <iostream>\n#include \"nr.h\"\n#include \"nums.hh\"\n#include \"tk.hh\"\n#include \"nv.hh\"\n#include \"modlab.h\"\n";
 aaa += "using namespace std;\n";
//definition of the list of fluxes, parameters and variables for "nv.cpp"
 aaa +=listint(flname,ifl,"nflx") + listint(par,ii,"nNV") + listint(newkin,knkin,"numx"); aaa+=("Fit Problem(nNV);\n"+sss);
 aaa +="double dt,xx[numx],flx[nflx],fluxes[nflx];\ndouble xm0,nv1[nNV],nv2[nNV],xinit1[numx],xinit2[numx];\n";
  aaa +=" int Parray::par[nNV]; string Parray::namepar[nNV], Parray::namef[nflx], Parray::namex[numx];\n"+rpar+fo;
  aaa +="void Parray::init(){ft3=10.; fh6=7.;\n\ttk.setk(&nv["+partk+"]);\n\tta.setk(&nv["+parta+"]);\n\taldolase.setk(&nv["+parald+"]);}\n";
  aaa +="void Parray::fin(){\n\t"+setald1+"\n\t"+setk1+"\n\t"+seta1+"\n\t";
  setald1[11]='2'; setk1[5]='2'; seta1[5]='2';
                                    aaa +=setald1+"\n\t"+setk1+"\n\t"+seta1+"\nflfor();\n}\n"+flfor;
  fout.open("../con512tpl/nv.cpp"); fout<<aaa;   fout.close();

 //distr.cpp - begin   
   string sdist="#include <iostream>\n#include \"nums.hh\"\n#include \"tk.hh\"\n#include \"nv.hh\"\n#include \"modlab.h\"\n";
 sdist += "using namespace std;\nvoid Ldistr::distr(double *py,double *pdydt) {\n\tdouble NOL=0.;\n\tsetiso(py); setdiso(pdydt);\n\t";
 sdist += "for (int i=0;i<Nn;i++) pdydt[i]=0.;\n\t";
  for(int i=iextra;i<kiso;i++) if(stype[i].substr(0,3)=="ald") sdist+=newiso[i]+".sett(); ";
                            else if(stype[i].substr(0,3)=="ket") sdist+=newiso[i]+".sett(); "; sdist+="\n";
     for(int i=0;i<kfl;i++) { if(isos[i]!="0"){
            if(isos[i][0]=='/') sdist +=isos[i].substr(1)+isom[i]+isop[i]+",fluxes["+flname[i]+"]/"+isos[i].substr(1)+".sumt()";
               else if(isop[i]!="0") sdist +=isos[i]+isom[i]+isop[i]+",fluxes["+flname[i]+"]";
           if(revfl[i]!="0") {if(revfl[i]!="NOL") sdist += ", fluxes["+revfl[i]+"]"; else sdist += ", 0.";}
           else if(isop[i]=="0") sdist +=isos[i]+isom[i]+"fluxes["+flname[i]+"]";
                               sdist +=");\n";}
                  else if(isop[i]!="0") sdist += isop[i]+".diso[0]+=fluxes["+flname[i]+"];\n";}
    sdist +="csyn("+ics1[0]+".iso,"+ics1[0]+".diso,"+ics2[0]+".iso,"+ics2[0]+".diso,"+icsp[0]+".diso,fluxes["+csfl[0]+"]);\n";
//    sdist +="csyn("+ics1[1]+".iso,"+ics1[1]+".diso,"+ics2[1]+".iso,"+ics2[1]+".diso,"+icsp[1]+".diso,fluxes["+csfl[1]+"]);\n";
    sdist +="split("+alds+".iso,"+alds+".diso,"+aldp+".iso,"+aldp+".diso,fluxes["+aldfl+"],fluxes["+aldfl+"+1]);\n";
    sdist +="spInvsl("+alds+".iso,"+alds+".diso,"+aldp+".iso,"+aldp+".diso,fluxes["+aldfl+"+2],"+aldp+".sumt());\n\t";

    sdist +=tak1+".tka<"+taa2.substr(1)+">("+taa1+","+taa2+","+tak2+",fluxes["+tafl+"]);\n\t";
    sdist +=tak2+".tka<"+taa1.substr(1)+">("+taa2+","+taa1+","+tak1+",fluxes["+tafl+"+1]);\n\t";
    sdist +=tak1+".invista("+taa1+",fluxes["+tafl+"+2]);\n\t";
    sdist +=tak2+".invista("+taa2+",fluxes["+tafl+"+3]);\n";
    
    sdist +=tkk1+".tkk<"+tka3.substr(1)+">("+tka1+","+tka3+","+tkk3+",fluxes["+tkfl+"]);\n";
    sdist +=tkk3+".tkk<"+tka1.substr(1)+">("+tka3+","+tka1+","+tkk1+",fluxes["+tkfl+"+1]);\n";
    sdist +=tkk2+".tkk<"+tka1.substr(1)+">("+tka2+","+tka1+","+tkk1+",fluxes["+tkfl+"+2]);\n";
    sdist +=tkk1+".tkk<"+tka2.substr(1)+">("+tka1+","+tka2+","+tkk2+",fluxes["+tkfl+"+3]);\n";
    sdist +=tkk2+".tkk<"+tka3.substr(1)+">("+tka2+","+tka3+","+tkk3+",fluxes["+tkfl+"+4]);\n";
    sdist +=tkk3+".tkk<"+tka2.substr(1)+">("+tka3+","+tka2+","+tkk2+",fluxes["+tkfl+"+5]);\n";
    sdist +=tkk1+".invistk("+tka1+",fluxes["+tkfl+"+6]);\n";
    sdist +=tkk2+".invistk("+tka2+",fluxes["+tkfl+"+7]);\n";
    sdist +=tkk3+".invistk("+tka3+",fluxes["+tkfl+"+8]);\n\t";
    for(int i=0;i<iextra;i++) sdist +=newiso[i]+".volume(Vt);\n\t";
                     sdist +="symm("+ssym+".getisot());\n}\n";
 fout.open("../con512tpl/distr.cpp"); fout<<sdist;  fout.close();
 
 //readexp.cpp -begin                                     
 aaa="#include <iostream>\n#include <cmath>\n#include \"nums.hh\"\n#include \"tk.hh\"\n#include \"nv.hh\"\n#include \"modlab.h\"\n";
 aaa += "using namespace std;\ndouble Vi, xribi, xasp=1.,mu;\ndouble Ldistr::readExp (string fn) {\n";
 aaa +="string aaa;  ifstream fi(fn.c_str()); double Ti,ts1;  mu=0.; fi>> dt; getline(fi,aaa);\n";
 aaa +="tex[0]=0.;  for(ntime=1;;ntime++) {fi>>tex[ntime]; tex[ntime] *= 60.; if(tex[ntime]<0) break;}\n";
 aaa +="fi>> Vi; getline(fi,aaa); double Nc[ntime]; for(int i=0;i<ntime;i++) fi>>Nc[i]; getline(fi,aaa);\n";
 aaa +="for(int i=1;i<ntime;i++) mu += log(Nc[i]/Nc[0])/tex[i]; mu/=1.;\n mu /= (double)(ntime-1); getline(fi,aaa);\n";
 for(int i=0;i<kiso;i++) if(edata[i][0]=='c') aaa +="\t"+newiso[i]+".readc(fi,c"+newiso[i]+",ntime);\n";
 for(int i=0;i<kiso;i++) if(edata[i][1]=='l') {int j=edata[i].length();
       if(j==2) aaa +=newiso[i]+".read(fi,f"+newiso[i]+",ntime,e"+newiso[i]+");\n";
    else {if(j>2) {sdist=newiso[i]+edata[i].substr(2,2); aaa +=newiso[i]+".read(fi,f"+sdist+",ntime,e"+sdist+");\n";}
       if(j>4) {sdist=newiso[i]+edata[i].substr(4,2); aaa +=newiso[i]+".read(fi,f"+sdist+",ntime,e"+sdist+");\n";}}
                      }
                aaa +="fi.close();\nreturn ts1;}\n";
 fout.open("../con512tpl/readexp.cpp"); fout<<aaa; fout.close();

//functions for "lab.cpp":
   aaa="#include <iostream>\n#include \"nr.h\"\n#include \"nums.hh\"\n#include \"modlab.h\"\n";
   aaa+="using namespace std;\n\nint Ldistr::getN() {\n"+newiso[0]+".ny=0;\n";
    for(int i=1;i<kiso;i++) aaa+=newiso[i]+".ny="+newiso[i-1]+".ny+"+newiso[i-1]+".getlen();\n";
                  aaa+="return ("+newiso[kiso-1]+".ny+"+newiso[kiso-1]+".getlen());\n}\n\n";

      aaa+="void Ldistr::setdiso(double *pyinit) {\n"+newiso[0]+".diso= &pyinit["+newiso[0]+".ny];\n";
    for(int i=1;i<kiso;i++) aaa +=newiso[i]+".diso= &pyinit["+newiso[i]+".ny];\n"; aaa+="}\n\n";

     aaa+="void Ldistr::setiso(double *pyinit) {\n"+newiso[0]+".iso= &pyinit["+newiso[0]+".ny];\n";
        for(int i=1;i<kiso;i++) aaa +=newiso[i]+".iso= &pyinit["+newiso[i]+".ny];\n";   aaa+="}\n\n";

  siso =aaa+"void Ldistr::massfr() {\n";
   for(int i=0;i<kiso;i++) if((edata[i][0]=='c')||(edata[i][0]=='l')) siso +=newiso[i]+".percent(); ";
      siso +="}\n\n";

  siso +="double Ldistr::xits(int its) {int itp=its-1; double xi=0;\n";
 for(int i=0;i<kiso;i++) {if(edata[i][0]=='c') {
     string sxic="xic"+newiso[i]+"[itp]";
      siso +=sxic+"="+newiso[i]+".chicon(its,c"+newiso[i]+"); xi +="+sxic+"*"+sxic+";\n";}
   if(edata[i].size()>2) sdist=newiso[i]+edata[i].substr(2,2); else sdist=newiso[i];
  {if(xi[i].size()>2)
   siso +="xi +=(xi"+sdist+"[itp]="+newiso[i]+".chisqsum("+xi[i]+",its,e"+sdist+","+stype[i].substr(stype[i].size()-2,1)+"));\n"; 
  else if(xi[i]=="1")
   siso +="xi +=(xi"+sdist+"[itp]="+newiso[i]+".chisq(its,e"+sdist+","+stype[i].substr(stype[i].size()-2,1)+"));\n";} }
     siso += "return xi;}\n\n";
     
   siso +="void Ldistr::wrifig(ostringstream& fo) {\n//t=1; "; string bbb="t";
    aaa="for(int its=0;its<ntime;its++) { fo<<tex[its];\n";
    for(int i=0;i<kiso;i++) {
      if(edata[i][0]=='c'){ siso +=" "+newiso[i]+"c="+bbb+"+1; "+newiso[i]+"sd="+newiso[i]+"c+1; ";bbb=newiso[i]+"sd";
          aaa+="fo<<\" \"<<c"+newiso[i]+"[its].mean<<\" \"<<c"+newiso[i]+"[its].sd;\n";}
        if(xi[i]=="1"){        siso +=" "+newiso[i]+"0="+bbb+"+1; "+newiso[i]+"0sd="+newiso[i]+"0+1; ";bbb=newiso[i]+"0sd";
          aaa+="fo<<\" \"<<e"+newiso[i]+"[its][0].mean<<\" \"<<e"+newiso[i]+"[its][0].sd;\n";}
        else if (edata[i].size()>2){ sdist=newiso[i]+edata[i].substr(2,2);  
           siso +=" "+sdist+"0="+bbb+"+1; "+sdist+"0sd="+sdist+"0+1; ";bbb=sdist+"0sd";
          aaa+="fo<<\" \"<<e"+sdist+"[its][0].mean<<\" \"<<e"+sdist+"[its][0].sd;\n";}
         else if(xi[i].size()>2){ siso +=" "+newiso[i]+"0="+bbb+"+1; "+newiso[i]+"0sd="+newiso[i]+"0+1; ";bbb=newiso[i]+"0sd";
               aaa+="fo<<\" \"<<e"+newiso[i]+"[its][0].mean<<\" \"<<e"+newiso[i]+"[its][0].sd;\n";}
              }
               aaa.erase(aaa.end()-2,aaa.end()); aaa+="<<\"\\n\";}\n}\n";
        siso+=";\n"+aaa; aaa="";

     siso+="void Ldistr::show(ostringstream& fo,double xfin) { fo<<setw(5)<<xfin;\n";
      aaa="//t=1; "; bbb="t";
    for(int i=0;i<kiso;i++) {
     if(edata[i][0]=='c') {siso +=newiso[i]+".showcon(fo); ";aaa+="c"+newiso[i]+"="+bbb+"+1; "; bbb="c"+newiso[i];}
      if(xi[i]=="1") { siso+=newiso[i]+".showm0(fo); ";aaa+=newiso[i]+"="+bbb+"+1; "; bbb=newiso[i];}
      else if(xi[i].size()>2){siso+=newiso[i]+".showsum("+xi[i]+",fo); ";aaa+="s"+newiso[i]+"="+bbb+"+1; ";bbb="s"+newiso[i];}
      }
             siso+="fo<<\"\\n\";}\n"+aaa+"\n"; aaa="";
             
     siso+="void Ldistr::showex(ostringstream& fo,int its) { fo.precision(4);\n fo<<setw(5)<<\"Data \"";
    for(int i=0;i<kiso;i++) {if(edata[i][0]=='c') {siso +="<<setw(8)<< c"+newiso[i]+"[its].mean";
                                                aaa +="<<setw(8)<< xic"+newiso[i]+"[itx]";}
                             if(xi[i]=="1") {siso +="<<setw(8)<< e"+newiso[i]+"[its][0].mean";
                                        aaa +="<<setw(8)<< xi"+newiso[i]+"[itx]";}
              if (edata[i].size()>2){ sdist=newiso[i]+edata[i].substr(2,2);  siso+="<<setw(8)<< e"+sdist+"[its][0].mean";
                                     aaa +="<<setw(8)<< xi"+sdist+"[itx]";}
              else if(xi[i].size()>2) {siso += "<<setw(8)<< e"+newiso[i]+"[its][0].mean";
                                     aaa +="<<setw(8)<< xi"+newiso[i]+"[itx]";}
              }
                 siso+="<<\"\\n\\n\";\nint itx=its-1;\nfo<<\"chi: \""+aaa+"<<\"\\n\";}\n\n";
                             
    siso+="double Ldistr::label() {\nreturn (";
    for(int i=iextra;i<kiso;i++) siso+=newiso[i]+".sumt()+";
       siso.erase(siso.end()-1); siso+=");}\n\n";
    siso+="void Ldistr::flback(){ cout<<\" ∑isotopomers-variable:\"<<\"\\n\";\ncout";
    for(int i=iextra;i<kiso;i++) siso+="<<\" "+newiso[i]+"=\"<<"+newiso[i]+".sumt()<<\"-\"<<xx[n"+newiso[i]+"]";
                           siso+="<<\"\\n\";}\n\n";
 fout.open("../con512tpl/lab.cpp"); fout<<siso; fout.close();/**/
     
    siso="#include <iostream>\n#include \"nums.hh\"\n#include \"modlab.h\"\nusing namespace std;\n";
                                        siso+="void Ldistr::ssc(double *pyinit) {setiso(pyinit);\n";
    for(int i=0;i<iextra;i++) if(edata[i][0]=='c') siso+=newiso[i]+".set0(c"+newiso[i]+"[0].mean);\n";
                              else siso+=newiso[i]+".set0(0.1);\n";
        siso+="gln.iso[0]=cgln[0].mean*egln[0][0].mean;\ngln.iso[31]=cgln[0].mean*egln[0][5].mean;\n";
        siso+="gl.iso[0]=cgl[0].mean*egl[0][0].mean;\ngl.iso[48]=cgl[0].mean*egl[0][2].mean;\n";
        siso+="gl.iso[49]=cgl[0].mean*egl[0][3].mean/4.;\ngl.iso[50]=gl.iso[49];\ngl.iso[52]=gl.iso[49];\n";
        siso+="gl.iso[56]=gl.iso[49];\ngl.iso[51]=cgl[0].mean*egl[0][4].mean/6.;\ngl.iso[53]=gl.iso[51];\n";
        siso+="gl.iso[54]=gl.iso[51];\ngl.iso[57]=gl.iso[51];\ngl.iso[58]=gl.iso[51];\ngl.iso[60]=gl.iso[51];\n";
     for(int i=iextra;i<kiso;i++) siso+=newiso[i]+".set0(xx[n"+newiso[i]+"]);\n";
                  siso+="}\n\n";
 fout.open("../con512tpl/ssc.cpp");    fout<<siso; fout.close();/**/

  
        stri= "";   for(int i=0;i<kiso;i++) stri += stype[i]+" "+newiso[i]+";\n";
       aaa="data "; for(int i=0;i<iextra;i++) if(edata[i][0]=='c') aaa +="c"+newiso[i]+"[tt],";
 for(int i=0;i<kiso;i++) if(edata[i][1]=='l') {int j=edata[i].size(); 
      int a=1+ncarb[i]; ostringstream ao; ao<<a;
       if(j==2) aaa +="e"+newiso[i]+"[tt]["+ao.str()+"],";
        else {if(j>2) {sdist=newiso[i]+edata[i].substr(2,2); aaa +="e"+sdist+"[tt]["+ao.str()+"],";}
             if(j>4) {sdist=newiso[i]+edata[i].substr(4,2); aaa +="e"+sdist+"[tt]["+ao.str()+"],";}}
                      }
                      aaa[aaa.size()-1]=';';   aaa +="\ndouble ";
    for(int i=0;i<iextra;i++) if(edata[i][0]=='c') aaa +="xic"+newiso[i]+"[tt],";
 for(int i=0;i<kiso;i++) if(edata[i][1]=='l') {int j=edata[i].size(); 
       if(j==2) aaa +="xi"+newiso[i]+"[tt],";
        else {if(j>2) {sdist=newiso[i]+edata[i].substr(2,2); aaa +="xi"+sdist+"[tt],";}}
                      }
                      aaa[aaa.size()-1]=';';   aaa +="\nint ";
 for(int i=0;i<kiso;i++) if(edata[i][1]=='l') {int j=edata[i].size(); 
       if(j==2) aaa +="f"+newiso[i]+",";
        else {if(j>2) {sdist=newiso[i]+edata[i].substr(2,2); aaa +="f"+sdist+",";}}
                      }
                      aaa[aaa.size()-1]=';';   aaa +="\n";
 fout.open("output"); fout<<stri<<aaa;  fout.close();
        return 0;}

