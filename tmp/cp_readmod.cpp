//---------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdlib>
#include <vector>

//---------------------------------------------------------------------------
using namespace std;
string partk="e0tk", parta="e0ta", ssym="fum";
string setk1, seta1, setald1;

inline string mm(string fln, string *sub, int nsub, ostringstream& parsets, bool ff=0){ //cout<<sub[i]<<'\n';
    string sfl="flx["+fln+"]", sub1="",ssub="", sdx="", foo =sfl+"= ";
     if(ff) foo += "rea[D].v()*"; int typv(0);
 for(int i=0;i<nsub;i++){ sub1=sub[i];
  if(sub1.at(0)=='-') if(sub1.size()>1){if(sub1.at(1)=='n') ssub +=", y["+sub1.erase(0,1)+"]";
                       else ssub +=", "+sub1.erase(0,1); typv++;}}
      parsets<<fln<<" "<<typv+1<<" v0= 0.01 ";
       for(int i=0;i<typv;i++) parsets<<"K"<<i<<"= 0.01 "; parsets<<'\n';
    ssub=ssub.erase(0,1);
      foo +="rea["+fln+"].v("+ssub+"); "; ssub="";
 for(int i=0;i<nsub;i++){ sub1=sub[i];
  if(sub1.at(0)=='-') {if(sub1.size()>1){if(sub1.at(1)=='n') ssub +="dydx["+sub1.erase(0,1)+"] -= "+sfl+";  ";}}
   else ssub +="dydx["+sub1+"] += "+sfl+";  ";
  }
         while(foo.length()<51) foo +=" "; foo +=ssub+'\n'; 
return foo;}

inline string rinput(ifstream& fi, string flname){string ssub,spr,flrev, rdist,sf="";
 fi>>ssub>>spr>>flrev; if(ssub.at(0)=='/') {sf=ssub+".sumt()"; ssub=ssub.substr(1); }
  rdist=ssub+".input("+ spr+",fluxes["+flname+"]"+sf;//read product
 if(flrev.length()-1) rdist += ",fluxes["+flrev+"]";//read reverse flux
   rdist += ");\n";
         return rdist;}

inline string irrev(ifstream& fi, string flname){string rtype, smet, rdist;
 fi>>rtype>>smet; rdist=smet+"."+rtype+"(";//read substrate
 fi>>smet; rdist += smet+",fluxes["+flname+"]);\n";//read product
         return rdist;}

inline string r3met(ifstream& fi){string rtype,smet0,smet1, smet2, flrev,flname, rdist;
 fi>>rtype>>smet0>>smet1>>smet2>>flname>>flrev;//read substrates
 cout << flname<<endl;
  rdist=smet0+"."+rtype+"(" + smet1+","+smet2+","+flname;
  if(flrev.length()-1) rdist +=",fluxes["+flrev+"]";
 rdist += ");\n";//read product
         return rdist;}

inline string rout(ifstream& fi, string flname){string rtype, smet, rdist;
 fi>>smet; rdist=smet+".output(fluxes["+flname+"]);\n";//read substrate
         return rdist;}

inline string rm0in(ifstream& fi, string flname){string rtype, smet, rdist;
 fi>>smet; rdist=smet+".diso[0]+=fluxes["+flname+"];\n";//read substrate
         return rdist;}

inline string ald(ostringstream& parset,string aldfl,string alds,string aldp){ 
    string sfl="flx["+aldfl+"]"; parset<<"ald   2  V= 0.00693115  K1= 14870\n";
          setald1 ="aldolase.st1fl(&"+sfl+", y[n"+alds+"], y[n"+aldp+"]);"; string foo=setald1;
      while(foo.length()<51) foo +=" ";  foo +="dydx[n"+alds+"] -= "+sfl+";";
      while(foo.length()<80) foo +=" ";  foo +="dydx[n"+aldp+"] += 2.*"+sfl+";\n";     
return foo;}

inline string ta(ostringstream& parset,string tafl,string tak1,string taa1,string taa2,string tak2){ 
    string sfl="flx["+tafl+"]"; parset<<"ta   5  v0= 0.000113773  K1= 1.45  K2= 1.6456  K3= 0.010418  K4= 0.544787\n";
      seta1 ="ta.st1fl(&"+sfl+", y[n"+tak1+"]/fh6, y[n"+taa1+"]/ft3, y[n"+taa2+"], y[n"+tak2+"]);"; string foo=seta1;
      foo +="\n\t\t\t\tdydx[n"+tak1+"] -= "+sfl+";\t dydx[n"+taa1+"] += "+sfl+";";
      foo +="\n\t\t\t\tdydx[n"+taa2+"] -= "+sfl+";\t dydx[n"+tak2+"] += "+sfl+";\n";    
return foo;}

inline string tk(ostringstream& parset,string tkfl,string tka1,string tkk1,string tka2,string tkk2,string tka3,string tkk3){ 
    string sfl="flx["+tkfl+"]", sfl1="flx["+tkfl+"+1]", sfl2="flx["+tkfl+"+2]"; 
    parset<<"tk   7  V= 1.29301e-07  K1= 0.03  K2= 2  K3= 3  K4= 1.63335  K5= 3  K6= 11\n";
  setk1="tk.st1fl(&"+sfl+", y[n"+tka1+"]/ft3, y[n"+tkk1+"], y[n"+tka2+"], y[n"+tkk2+"]/fh6, y[n"+tka3+"], y[n"+tkk3+"]);";
     string foo=setk1+"\n\t\t\t\tdydx[n"+tkk1+"] -= "+sfl+";\tdydx[n"+tka1+"] += "+sfl+";\n"; 
      foo +="\t\t\t\tdydx[n"+tkk3+"] -= "+sfl1+";\tdydx[n"+tka3+"] += "+sfl1+";\n";    
      foo +="\t\t\t\tdydx[n"+tkk2+"] -= "+sfl2+";\tdydx[n"+tka2+"] += "+sfl2+";\n";    
return foo;}

inline string listint(vector<string> sname,string sfin,string sin="0"){
  string aaa="const int "+sname[0]+"="+sin+", "; int inum=sname.size();
   for(int i=1;i<inum;i++)   aaa+=sname[i]+"="+sname[i-1]+"+1, ";
     aaa+=sfin+"="+sname[inum-1]+"+1;\n\n";
 return aaa;}
inline string listname(vector<string> sname, string sfin){
  string aaa="extern const int "; int inum=sname.size();
   for(int i=0;i<inum;i++)   aaa+=sname[i]+", ";   aaa+=sfin+";\n\n";
 return aaa;}
  inline string setvals(vector<string> snames, string sval){ string aaa=""; int nelem=snames.size();
    for(int i=0;i<nelem;i++) {ostringstream snum; snum<<i; aaa+=snum.str()+") "+snames[i]+"\t"+sval+"\n";}
      return aaa;}

int main ( int argc, char *argv[] ) {
   ifstream fi("model");  int knkin;
     string aaa="", stri, siso;
     int lmet=0, lmetb=0, lmetk=0;
     ostringstream smeta,smetb,smetk;
     ostringstream scona,sconb,sconk;
// variables of kinetic model   
       vector<string> stype, mname, edata, xi;// metabolite type, name, edata(conc,lab,0), chi(yes,no)
         getline(fi,aaa);
   for(int i=0; ;i++) {
    fi>> aaa; if(aaa=="fin") break;
     stype.push_back(aaa);//type of metabolite
    fi>> aaa; mname.push_back(aaa);//names/numbrs of metabolites met[i]->setconc(
     if(stype[i].find("Metab_d")+1) { smeta <<"met["<<lmet<<"]=&"<<mname[i].substr(1)<<"; ";//setmet()
         scona <<"met["<<lmet<<"]->setconc(xx["<<mname[i]<<"]); "; lmet++;}//setcon(),nv.cpp
     if(stype[i].find("Metab<")+1) { smetb <<"metb["<<lmetb<<"]=&"<<mname[i].substr(1)<<"; ";
      sconb <<"metb["<<lmetb<<"]->setconc(xx["<<mname[i]<<"]); "; lmetb++;}
     if(stype[i].find("keto")+1) { smetk <<"metk["<<lmetk<<"]=&"<<mname[i].substr(1)<<"; ";
     sconk <<"metk["<<lmetk<<"]->setconc(xx["<<mname[i]<<"]); "; lmetk++;}
              if(aaa=="numx") knkin=i;
    fi>> aaa; edata.push_back(aaa);//names/numbrs of metabolites
    fi>> aaa; xi.push_back(aaa);//names/numbrs of metabolites
   }
     int nmet=stype.size();
    ostringstream parsets;
    smetk<<"\n lmet="<<lmet<<"; lmetb="<<lmetb<<"; lmetk="<<lmetk<<"; ";
   cout<<"variables: "<<nmet<<endl;
//reding the scheme of the model:
       
     vector<string> flname, fladd;
     int numsub, eq;
//     nv.cpp:
  string f_ff="void Fit::f(const double *y,double *dydx) {\n\tfor(int i=0;i<numx;i++) dydx[i]=0.;\n\t for(int i=0;i<nflx;i++) flx[i]=0.;\n\t double amp = -(sqrt(4.*xx[n_atp]*tan-3.*xx[n_atp]*xx[n_atp])-2.*tan+xx[n_atp])/2.;\n\t double a_dp = (sqrt(xx[n_atp])*sqrt(4.*tan-3.*xx[n_atp])-xx[n_atp])/2.;\n\t double h_nad = tnad-xx[n_nad];\n\t xthf=thft-xx[ncthf];\n";// functions f & ff
  string flfor="void Parray::flfor(double *y){\nfor(int i=0;i<nflx;i++) fluxes[i] = flx[i] * dt/Vi;\n";// function flfor
       string stdist;
    
// reactions (fluxes)
   for(int i=0; ;i++){ fi>>aaa; if(aaa.substr(0,3)=="fin") break;
                       flname.push_back(aaa);
                       fi>>numsub; string subst[numsub]; //substrates
           string s_fl ="fluxes["+aaa+"] /= "; string sy=""; int is(0);
     for(int i1=0;i1<numsub;i1++){  fi>>subst[i1];                            //read substrates
      if(subst[i1].at(0)=='-')if(subst[i1].size()>2)if(subst[i1].at(2)!='_') {
           is++;  if(is>1) s_fl +="*";
           sy=subst[i1];
           sy.erase(0,1);
           if(sy.at(0)=='n') s_fl +="y[";
           s_fl +=sy;
           if(sy.at(0)=='n') s_fl +="]"; } }
               if(sy.size()>2) flfor +=s_fl+";\n";
               fi>>eq; parsets<<i<<" ";
       f_ff += mm(aaa,subst,numsub,parsets,eq); // function f
    fi>>aaa; if(aaa=="input"){aaa=rinput( fi, flname[i]); stdist+=aaa; cout<<aaa<<endl;}
      else if(aaa=="irr"){aaa=irrev( fi, flname[i]); stdist+=aaa; cout<<aaa<<endl;}
      else if(aaa=="output"){aaa=rout( fi, flname[i]); stdist+=aaa; cout<<aaa<<endl;}
      else if(aaa=="m0inp"){aaa=rm0in( fi, flname[i]); stdist+=aaa; cout<<aaa<<endl;}
      else if(aaa=="r3met"){aaa=r3met( fi); stdist+=aaa; cout<<aaa<<endl;}
      else getline(fi,aaa);
         }
      int kfl=flname.size();  cout<<"fluxes: "<<kfl<<'\n';

     string aldfl, alds, aldp; int nr;
   fi>>aaa>>aldfl>>alds>>aldp;//aldolase
   flname.push_back(aldfl);
   parsets<<kfl<<" ";
      f_ff += ald(parsets,"aldfl", alds, aldp);
       stdist+=r3met( fi);
           fladd.push_back("aldfl"); fladd.push_back("aldrev");
           fladd.push_back("aldfli");  fladd.push_back("aldi1");
                 flfor +="fluxes["+aldfl+"] /= xx[n"+ alds+ "];\n";
                 flfor +="fluxes["+aldfl+"+1] /= (xx[n"+ aldp+ "]*xx[n"+ aldp+ "]);\n";
                 flfor +="fluxes["+aldfl+"+2] /= xx[n"+ alds+ "];\n";
     string tafl, tak1,taa1,taa2,tak2;
   fi>>aaa>>nr>>tafl>>tak1>>taa1>>taa2>>tak2;//transaldolase
   flname.push_back(tafl);
   parsets<<kfl+1<<" ";
      f_ff += ta(parsets,"tafl", tak1,taa1,taa2,tak2);
       for(int i=0;i<nr;i++) stdist+=("//"+r3met( fi));
       fladd.push_back("tafl"); fladd.push_back("s7f6a");
       fladd.push_back("f6g3a"); fladd.push_back("s7e4a");
   
     string tkfl, tka1,tkk1,tka2,tkk2,tka3,tkk3;
   fi>>aaa>>nr>>tkfl>>tka1>>tkk1>>tka2>>tkk2>>tka3>>tkk3;//transketoolase
   flname.push_back(tkfl);
   parsets<<kfl+2<<" ";
      f_ff += tk(parsets,"tkfl", tka1,tkk1,tka2,tkk2,tka3,tkk3); stdist+='\n';
       for(int i=0;i<nr;i++) stdist+=r3met( fi);
       fladd.push_back("tkfl"); fladd.push_back("s7p5");
       fladd.push_back("f6p5"); fladd.push_back("p5f6");
       fladd.push_back("f6s7"); fladd.push_back("s7f6");
       fladd.push_back("p5g3i"); fladd.push_back("f6e4i"); fladd.push_back("s7p5i");
               f_ff+="for(int i=0;i<numx;i++) dydx[i]*=(dt/Vi);\n}\n\n"; flfor+="}\n\n";
   
          fi>>aaa;
          f_ff+="void Fit::ff(const double *y,double *dydx) {\n";
          vector<string> smet; int numfl;
       for(int i=0; ;i++) {        fi>>aaa;
         if(aaa.substr(0,3)=="fin") break; smet.push_back(aaa); f_ff+="\tdydx["+aaa+"] = (";
                                  fi>>numfl;
        for(int j=0;j<numfl;j++) {fi>>aaa;
        if(aaa.at(0)=='-') {f_ff+="- "; aaa.erase(0,1);} else f_ff+="+"; f_ff+="flx["+aaa+"]";}
         f_ff+=")*dt;\n";
       }  f_ff+="}\n\n";
                  
         fi.close();
//     cout<<f_ff<<endl;
         
//setting ODE model, conversion of fluxes for iso model (flfor):
      
             
//declaration of the list of fluxes and variables for "nums.hh"
  ofstream numshh("../include/nums.hh");
  string hlfl=listname(flname,"nrea")+listname(fladd,"nflx") + listname(mname,"nmet");//+listname(smet,"nmet");
     hlfl+="\nextern const double thft;\n";
     hlfl+="extern double mader, dt,dif,dif0, suxx, Vi,Vt, mu, xx[], fluxes[], flx[];\n";
     hlfl+="extern double nv1[],nv2[],xinit1[],xinit2[];\nextern time_t ts,tf,tcal;\n";
     hlfl+="extern std::string foc, kin, kinc;\nextern char *fex1, *fex2;\nextern int ifn;\n";
          numshh<<hlfl; numshh.close(); hlfl=""; //writing "nums.hh"
 //file with parameters values
               parsets<<"1 3 5 7 9 11 -1\n"+setvals(mname, "0.01")+setvals(flname, "0.001");
      ofstream param("parameters");
  param<<parsets.str(); param.close(); parsets.str("");
  
  parsets<<"Metab_data *Ldistr::met["<<lmet<<"];\nMetab *Ldistr::metb["<<lmetb<<"];\n ketose *Ldistr::metk["<<lmetk<<"];\n";
  
    
string nvcpp="#include <iostream>\n#include \"nr.h\"\n#include \"nums.hh\"\n#include \"tk.hh\"\n#include \"nv.hh\"\n#include \"modlab.h\"\n#include \"solvers.h\"\n#include \"analis.h\"\n";
 nvcpp += "using namespace std;\n";
//definition of the list of fluxes, parameters and variables for "nv.cpp"
 nvcpp +=listint(flname,"nrea")+listint(fladd,"nflx","rald") + listint(mname,"nmet");//+listint(smet,"nmet","numx");
  nvcpp+=("\tFit Problem;\n\tconst double thft(1.);\n\tdouble dt,xx[nmet],flx[nflx],fluxes[nflx];\n\tdouble xinit1[nmet],xinit2[nmet];\n\t");
  nvcpp +="string Parray::fid[nflx],Parray::fname[nflx],Parray::fschem[nflx], Parray::namex[numx];\n\tReapar Parray::rea[nrea];\n\tdouble Analis::nv1[nrea], Analis::nv2[nrea];\n"+parsets.str()+"\n void Ldistr::setmet(){";
   parsets.str("");
  nvcpp +=smeta.str()+smetb.str()+smetk.str() + " }\n void Ldistr::setcon(){"+
    scona.str()+sconb.str()+sconk.str()+ " }\n"+ f_ff;
  nvcpp +="void Parray::init(){ft3=10.; fh6=7.;\n\ttk.setk(rea[rtk].getpar());\n\tta.setk(rea[rta].getpar());\n\taldolase.setk(rea[rald].getpar());}\n";
  nvcpp +="void Parray::fin(double y[]){\n\t"+setald1+"\n\t"+setk1+"\n\t"+seta1+"\n\t";
  setald1[11]='2'; setk1[5]='2'; seta1[5]='2';
                                    nvcpp +=setald1+"\n\t"+setk1+"\n\t"+seta1+"\nflfor(y);\n}\n"+flfor;
  ofstream fout("../con512tpl/nv.cpp"); fout<<nvcpp;   fout.close(); nvcpp="";

 //distr.cpp - begin   
   string sdist="#include <iostream>\n#include \"nums.hh\"\n#include \"tk.hh\"\n#include \"nv.hh\"\n#include \"modlab.h\"\n";
 sdist += "using namespace std;\nvoid Ldistr::distr(double *py,double *pdydt) {\n\tdouble NOL=0.;\n\tsetiso(py); setdiso(pdydt);\n\t";
 sdist += "for (int i=0;i<Nn;i++) pdydt[i]=0.;\n\t";
  for(int i=0;i<nmet;i++) if((stype[i].size()>3)&&(stype[i].substr(0,3)=="ket")) sdist+=mname[i].substr(1)+".sett(); "; sdist+="t3.sett(); e4.sett();\n"; 
//     cout<<mname[i]<<endl;
   sdist += "\tProblem.f(py,pdydt);\n\tProblem.ff(py,pdydt);\n\tProblem.fin(py);/**/\n\tdouble xthf=thft-cthf.sumt();\n";
    sdist +=stdist;
//    sdist +="csyn("+ics1[1]+".iso,"+ics1[1]+".diso,"+ics2[1]+".iso,"+ics2[1]+".diso,"+icsp[1]+".diso,fluxes["+csfl[1]+"]);\n";
//    sdist +="split("+alds+".iso,"+alds+".diso,"+aldp+".iso,"+aldp+".diso,fluxes["+aldfl+"],fluxes["+aldfl+"+1]);\n";
    sdist +="//spInvsl("+alds+".iso,"+alds+".diso,"+aldp+".iso,"+aldp+".diso,fluxes["+aldfl+"+2],"+aldp+".sumt());\n\t";

    sdist +=tak1+".tka("+taa1+","+taa2+","+tak2+",fluxes[tafl]);\n\t";
    sdist +=tak2+".tka("+taa2+","+taa1+","+tak1+",fluxes[tafl+1]);\n\t";
    sdist +=tak1+".invista("+taa1+",fluxes[tafl+2]);\n\t";
    sdist +=tak2+".invista("+taa2+",fluxes[tafl+3]);\n//";
    
    sdist +=tkk1+".tkk("+tka1+","+tka3+","+tkk3+",fluxes[tkfl]);\n//";
    sdist +=tkk3+".tkk("+tka3+","+tka1+","+tkk1+",fluxes[tkfl+1]);\n//";
    sdist +=tkk2+".tkk("+tka2+","+tka1+","+tkk1+",fluxes[tkfl+2]);\n//";
    sdist +=tkk1+".tkk("+tka1+","+tka2+","+tkk2+",fluxes[tkfl+3]);\n//";
    sdist +=tkk2+".tkk("+tka2+","+tka3+","+tkk3+",fluxes[tkfl+4]);\n//";
    sdist +=tkk3+".tkk("+tka3+","+tka2+","+tkk2+",fluxes[tkfl+5]);\n//";
    sdist +=tkk1+".invistk("+tka1+",fluxes[tkfl+6]);\n//";
    sdist +=tkk2+".invistk("+tka2+",fluxes[tkfl+7]);\n//";
    sdist +=tkk3+".invistk("+tka3+",fluxes[tkfl+8]);\n\t";
    for(int i=knkin+1;i<nmet;i++) sdist +=mname[i].substr(1)+".volume(Vt);\n\t";
                     sdist +="symm("+ssym+".getisot());\n//";
    for(int i=0;i<nmet;i++) if(stype[i]!="0") sdist += "xx["+mname[i]+"]="+mname[i].substr(1)+".sumt(); ";  sdist += "\n";
    for(int i=0;i<nmet;i++) if(stype[i]!="0") sdist += "xx["+mname[i]+"]=py["+mname[i]+"]; "; 
    sdist += "\n}\n";
 fout.open("../con512tpl/distr.cpp"); fout<<sdist;  fout.close(); sdist="";
 
 stri= "Metab_data"; for(int i=0;i<lmet;i++) if(stype[i]!="0") stri += " "+mname[i].substr(1)+";\n";
 stri += "void sklad(int itime){\n";
 for(int i=0;i<nmet;i++) {if(edata[i][0]=='c') stri += mname[i].substr(1)+".skladc(itime); ";
                          if(edata[i][1]=='l') stri += mname[i].substr(1)+".skladm0(itime); ";}
                          stri += "}\n";
 stri += "void wrikin(std::ostringstream& so, int nt){\n";
 for(int i=0;i<nmet;i++) {if(edata[i][0]=='c') stri += mname[i].substr(1)+".wrikinc(\""+mname[i]+"\",so,nt); ";
                          if(edata[i][1]=='l') stri += mname[i].substr(1)+".wrikinm0(\""+mname[i]+"\",so,nt); ";}
                          stri += "}\n";
 fout.open("output"); fout<<stri;  fout.close();
        return 0;}

