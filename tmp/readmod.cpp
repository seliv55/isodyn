//---------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdlib>
#include <vector>
#include <tuple> 

//---------------------------------------------------------------------------
using namespace std;
string partk="e0tk", parta="e0ta", ssym="mal";
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
   ifstream fi("model");  int  nxin(0); bool bx(0);
     string aaa="", siso,snx;
     int lmet=0, lmetk=0;
     ostringstream smeta,smetb,smetk;
     ostringstream scona,sconb,sconk;
     ostringstream sdata,sket,sMe;
// variables of kinetic model   
       vector<tuple<int, string, string>> mdat, mdatk;
       vector<string> mname;
         getline(fi,aaa);
   for(int i=0; ;i++) {//names of metabolite
    fi>> aaa; if(aaa=="fin") break;
         else if(aaa.find("intern_")+1) {nxin=i; bx=1;}
         else {tuple<int, string, string> namc; int nc; string sss, snm;
          fi>>nc>>sss>>snm; mname.push_back(sss); namc=make_tuple(nc,sss,snm);
               if(bx) {snx="const int numx="+sss+";\n\n"; bx=0;}
           if(aaa.find("Metab")+1) { mdat.push_back(namc);
              smeta <<"met["<<lmet<<"]=&"<<get<1>(namc).substr(1)<<"; ";
              scona <<"met["<<lmet<<"]->setconc(xx["<<get<1>(namc)<<"]); ";
              sdata<<" Ldistr::"<<sss.substr(1)<<"("<<nc<<","<<snm<<"),";
                  lmet++; }
           else if(aaa.find("keto")+1) { mdatk.push_back(namc);
              smetk <<"metk["<<lmetk<<"]=&"<<get<1>(namc).substr(1)<<"; ";
              sconk <<"metk["<<lmetk<<"]->setconc(xx["<<get<1>(namc)<<"]); ";
              sket<<" Ldistr::"<<sss.substr(1)<<"("<<nc<<","<<snm<<"),";
                  lmetk++; }
                }
              }//names of metabolite
   
     int nmet=mname.size();
    ostringstream parsets;
    lmet=mdat.size();   lmetk=mdatk.size(); 
    smetk<<"\n lmet="<<lmet<<"; lmetk="<<lmetk<<"; ";
   cout<<"variables: "<<nmet<<endl;
//reding the scheme of the model:
       
     vector<string> flname, fladd, bal;
     for(int i=0; ;i++){
      getline(fi,aaa); if(aaa.find("end_bal")+1) break;
      bal.push_back(aaa);
      cout<<bal[i]<<'\n';
     }
     int numsub, eq;
//     nv.cpp:
  string f_ff="void Fit::f(const double *y,double *dydx) {\n\tfor(int i=0;i<nmet;i++) dydx[i]=0.;\n\tfor(int i=0;i<nflx;i++) flx[i]=0.;\n";// functions f & ff
  for(int i=1;i<bal.size();i++) f_ff +="\t double " + bal[i] +"\n";
   bal.clear();
  string flfor="void Parray::flfor(double *y){\nfor(int i=0;i<nflx;i++) fluxes[i] = flx[i] * flx[rdt]/Vi;\n";
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
    fi>>aaa; if(aaa=="input"){aaa=rinput( fi, flname[i]); stdist+=aaa; }
      else if(aaa=="irr"){aaa=irrev( fi, flname[i]); stdist+=aaa; }
      else if(aaa=="output"){aaa=rout( fi, flname[i]); stdist+=aaa; }
      else if(aaa=="m0inp"){aaa=rm0in( fi, flname[i]); stdist+=aaa; }
      else if(aaa=="r3met"){aaa=r3met( fi); stdist+=aaa; }
      else getline(fi,aaa);
         }
      int kfl=flname.size();  cout<<"fluxes: "<<kfl<<'\n';
      fi>>aaa; 
     string tafl, tak1,taa1,taa2,tak2;
   
     string tkfl, tka1,tkk1,tka2,tkk2,tka3,tkk3;
      cout<<"parsets="<<parsets.str()<<'\n';
               f_ff+="for(int i=0;i<nmet;i++) dydx[i]*=(flx[rdt]/Vi);\n}\n\n"; flfor+="}\n\n";
   
//          fi>>aaa;
          f_ff+="void Fit::ff(const double *y,double *dydx) {\n\tf(y,dydx);\n";
          vector<string> smet; int numfl;
       for(int i=0; ;i++) {        fi>>aaa;
cout<<" *** pass ***"<<aaa<<endl;
         if(aaa.substr(0,3)=="fin") break;
         smet.push_back(aaa); f_ff+="\tdydx["+aaa+"] = (";
                                  fi>>numfl;
        for(int j=0;j<numfl;j++) {fi>>aaa;
        if(aaa.at(0)=='-') {f_ff+="- "; aaa.erase(0,1);}
        else f_ff+="+";
        f_ff+="flx["+aaa+"]";}
         f_ff+=")*flx[rdt];\n";
       }  f_ff+="}\n\n";
                  
         fi.close();
         
//setting ODE model, conversion of fluxes for iso model (flfor):
//declaration of the list of fluxes and variables for "nums.hh"
  ofstream numshh("../include/nums.hh");
  string hlfl=listname(flname,"nrea")+listname(fladd,"nflx") + listname(mname,"numx, nmet");//+listname(smet,"nmet");
     hlfl+="\nextern const double thft;\n";
     hlfl+="extern double mader, dif,dif0, suxx, Vi,Vt, mu, xx[], fluxes[], flx[];\n";
     hlfl+="extern double nv1[],nv2[],xinit1[],xinit2[];\nextern time_t ts,tf,tcal;\n";
     hlfl+="extern std::string foc, kin, kinc;\nextern char *fex1, *fex2;\nextern int ifn;\n";
          numshh<<hlfl; numshh.close(); hlfl=""; //writing "nums.hh"
 //file with parameters values
               parsets<<"1 3 5 7 9 11 -1\n"+setvals(mname, "0.01")+setvals(flname, "0.001");
      ofstream param("parameters");
  param<<parsets.str(); param.close(); parsets.str("");
  
  parsets<<"Metab *Ldistr::met["<<lmet<<"];\n";
  if(lmetk>0)  parsets<<" ketose *Ldistr::metk["<<lmetk<<"];\n";
  
    
string nvcpp="#include <iostream>\n#include \"nr.h\"\n#include \"nums.hh\"\n#include \"tk.hh\"\n#include \"nv.hh\"\n#include \"modlab.h\"\n#include \"solvers.h\"\n#include \"analis.h\"\n";
 nvcpp += "using namespace std;\n";
//definition of the list of fluxes, parameters and variables for "nv.cpp"
 nvcpp +=listint(flname,"nrea")+"const int nflx=nrea;\n" + listint(mname,"nmet")+snx;//+listint(smet,"nmet","numx")+listint(fladd,"nflx","rald");
string aa="Metab"+sdata.str(); aa.erase(aa.end()-1); 
if(lmetk>0)  {aa+=";\nketose"+sket.str(); aa.erase(aa.end()-1);}
 nvcpp +=aa+";\n";
  nvcpp+=("\tFit Problem;\n\tconst double thft(1.);\n\tdouble xx[nmet],flx[nflx],fluxes[nflx];\n\tdouble xinit1[nmet],xinit2[nmet];\n\t");
  nvcpp +="string Parray::fid[nflx],Parray::fname[nflx],Parray::fschem[nflx], Parray::namex[nmet];\n\tReapar Parray::rea[nrea];\n\tdouble Analis::nv1[nrea], Analis::nv2[nrea];\n"+parsets.str()+"\n void Ldistr::setmet(){";
   parsets.str("");
  nvcpp +=smeta.str()+smetb.str()+smetk.str() + " }\n void Ldistr::setcon(){"+
    scona.str()+sconb.str()+sconk.str()+ " }\n"+ f_ff;
  nvcpp +="void Parray::init(){ft3=10.; fh6=7.;}\n";
  nvcpp +="void Parray::fin(double y[]){\n\t"+setald1+"\n\t"+setk1+"\n\t"+seta1+"\n\t";
  setald1[11]='2'; setk1[5]='2'; seta1[5]='2';
nvcpp +=setald1+"\n\t"+setk1+"\n\t"+seta1+"\nflfor(y);\n}\n"+flfor;
  ofstream fout("../con512tpl/nv.cpp"); fout<<nvcpp;   fout.close(); nvcpp="";

 //distr.cpp - begin   
   string sdist="#include <iostream>\n#include \"nums.hh\"\n#include \"tk.hh\"\n#include \"nv.hh\"\n#include \"modlab.h\"\n";
 sdist += "using namespace std;\nvoid Ldistr::distr(double *py,double *pdydt) {\n\tdouble NOL=0.;\n\tsetiso(py); setdiso(pdydt);\n\t";
 sdist += "for (int i=0;i<Nn;i++) pdydt[i]=0.;\n\t";
  for(int i=0;i<mdatk.size();i++) sdist+=get<1>(mdatk[i]).substr(1)+".sett(); "; sdist+="\n"; 
   sdist += "\tProblem.ff(py,pdydt);\n\tProblem.fin(py);/**/\n";
    sdist +=stdist;
    for(int i= nxin;i<nmet;i++) sdist +=mname[i].substr(1)+".volume(Vt);\n\t";
                     sdist +="symm("+ssym+".getisot());\n//";
    for(int i=0;i<mdat.size();i++) sdist += "xx["+get<1>(mdat[i])+"]="+get<1>(mdat[i]).substr(1)+".sumt(); ";  sdist += "\n//";
    for(int i=0;i<mdatk.size();i++) sdist += "xx["+get<1>(mdatk[i])+"]="+get<1>(mdatk[i]).substr(1)+".sumt(); ";  sdist += "\n//";
    
    for(int i=0;i<nmet;i++) sdist += "xx["+mname[i]+"]=py["+mname[i]+"]; "; 
    sdist += "\n}\n";
 fout.open("../con512tpl/distr.cpp"); fout<<sdist;  fout.close(); sdist="";
 
 string stri= "Metab "+get<1>(mdat[0]).substr(1);
 ostringstream stcons;
  stcons<< (get<1>(mdat[0]).substr(1)+"(")<<get<0>(mdat[0])<<(", "+get<2>(mdat[0])+")");
  for(int i=1;i<mdat.size();i++){stri += ", "+get<1>(mdat[i]).substr(1);
  stcons<<(", "+get<1>(mdat[i]).substr(1)+"(")<<get<0>(mdat[i])<<(", "+get<2>(mdat[i])+")");
  }   stri += ";\n";
  
 if(lmetk>0){
  stri += "ketose "+get<1>(mdatk[0]).substr(1);
  for(int i=1;i<mdatk.size();i++) stri += ", "+get<1>(mdatk[i]).substr(1);  stri += ";\n";
  for(int i=0;i<mdatk.size();i++)
  stcons<<(", "+get<1>(mdatk[i]).substr(1)+"(")<<get<0>(mdatk[i])<<(", "+get<2>(mdatk[i])+")");
 }  
 fout.open("output"); fout<<stri<<stcons.str();  fout.close();
        return 0;}

