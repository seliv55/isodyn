// clearing vectors
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
inline string mm(string fln, vector<string> sub, vector<string> prd, ostringstream& parsets, bool ff=0){ //cout<<sub[i]<<'\n';
    string sfl="flx["+fln+"]", sub1="",ssub="", sdx="", foo =sfl+"= ";
     if(ff) foo += "rea[D].v()*"; int typv(0);
 for(unsigned i=0;i<sub.size();i++){ sub1=sub[i];
   if(sub1.size()>1){if(sub1.at(0)=='n') ssub +=", y["+sub1+"]";
                       else ssub +=", "+sub1; typv++;}}
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

int main ()
{ ifstream fi("train"); string aaa,bbb;
    ostringstream parsets;
     fi>>aaa; //if(aaa.substr(0,3)=="fin") break;
     vector<string> flname, subst, prod; //substrates
     flname.push_back(aaa);
     string s_fl ="fluxes["+aaa+"] /= "; string sy=""; int is(0);
      int i1=0;
     while(bbb!="->"){
       fi>>bbb;  if(bbb=="->") break;              //read substrates
       subst.push_back(bbb); 
       if(subst[i1].size()>2) {
           is++;  if(is>1) s_fl +="*";
           sy=subst[i1];
           if(sy.at(0)=='n') s_fl +="y[";
           s_fl +=sy;
           if(sy.at(0)=='n') s_fl +="]"; }
           i1++;}
cout<<bbb<<'\n'; i1=0;
     while(bbb!="$"){ 
       fi>>bbb;  if(bbb=="$") break;      //read substrates
       if(bbb != "+") {
       prod.push_back(bbb); 
       if(prod[i1].size()>2) {
           is++;  if(is>1) s_fl +="*";
           sy=prod[i1];
           if(sy.at(0)=='n') s_fl +="y[";
           s_fl +=sy;
           if(sy.at(0)=='n') s_fl +="]"; }
           cout<<prod[i1]<<'\n';
       i1++;}
           } cout<<s_fl<<'\n'; i=0
               fi>>eq; parsets<<i<<" ";
       f_ff += mm(aaa,subst,prod,parsets,eq); // function f
      
   std::vector<std::string> myvector;
  std::cout <<myvector.size()<<'\n';
  myvector.push_back ("asd");
  std::cout <<myvector[0]<<'\n';
  myvector.push_back ("ert");
  myvector.push_back ("xcv");

  std::cout << "myvector contains:";
  for (unsigned i=0; i<myvector.size(); i++)
    std::cout << ' ' << myvector[i];
  std::cout << '\n';

  myvector.clear();
  myvector.push_back ("abcd");
  myvector.push_back ("efgh");

  std::cout << "myvector contains:";
  for (unsigned i=0; i<myvector.size(); i++)
    std::cout << ' ' << myvector[i];
  std::cout << '\n';

  return 0;
}
