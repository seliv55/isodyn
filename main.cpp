//---------------------------------------------------------------------------
#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
#include "analis.h"
#include "solvers.h"
//---------------------------------------------------------------------------
//#pragma package(smart_init)
using namespace std;
 Ldistr horse;
  int  ifn=0; 
   double dlt, dif, suxx, mader;
//   const double flf=60.*1000.*5./7.;
    time_t ts,tf; 
  extern string foc, kin0, kin, kinflx;
    
tuple<double,double,time_t> solve(){
//solving the kinetic model:
   ts=clock(); double xi;
 for(int i=0;i<nmet;i++) xx[i]=xinit1[i]; //take initial values
//    tsolve(1400.); //mader=Problem.dermax();
    ddsolve(1400.); //mader=Problem.dermax();
//    Problem.shownx(nmet,xx);
//   xi=horse.ddisolve();
    // Problem.jacobian(xx);  tsolve(37000.); for(int i=0;i<numx;i++) xinit1[i]=xx[i];
      xi=horse.integrbs();
  tf=clock()-ts;
   suxx=horse.consum(); // suxx += dif;
//   for(int i=0;i<nmet;i++) xx[i]=xinit1[i]; //set iv for a case of saving them
   return make_tuple(xi,suxx,tf);}
   
inline void chekxi(char *efi){
   stringstream fn; int iin, ifin;  double xi0;
    cout<<"first file: "; cin>>iin;
    cout<<"\nlast file: "; cin>>ifin; 
     for(int i=iin;i<ifin;i++) {fn<<*Problem.getodir()<<i; cout<<i<<": "<<endl;
      Problem.read(fn.str().c_str()); cout<<fn.str().c_str()<<endl;
      horse.readExp(efi); for(int ii=0;ii<numx;ii++) cout<<(xinit1[ii]=xx[ii])<<endl;
        tuple<double,double,time_t> sol=solve();            
   Problem.write(sol,i,0);      } /**/
  }

    tuple<double,double,time_t> sol0;
    
string Ldistr::read_con(ifstream& fi){
// read par-file, output dir, flux conf. inter: main and to compare
 string aaa,gpl="gnuplot  -e \"con=3;m0=4\" xplt.p"; 
     int isu; 
      fi>>aaa;
      getline(fi,aaa);
      stringstream onel(aaa);  int ii(0);
      while(getline(onel,aaa,' ')) if(aaa.length()){double ta=stod(aaa); texcon.push_back(ta); ii++; cout<<ta<<" ";}
       ntime=ii; 
      while(!fi.eof()){ onel.str(std::string());
       fi>>aaa; if(fi.eof()) break; int k;
        for(k=0;k<=lmet;k++) if(aaa.find(met[k]->getdescr())+1) {
           for(int j=0;j<ntime;j++) {
       double cc; fi>>cc; met[k]->setconc(cc,j);
       cout<<"t="<<texcon[j]<<" "<<met[k]->getdescr()<<" "<<met[k]->getconc()[j].mean<<'\n';
           } getline(fi,aaa);
//       xx[k]=met[k]->getconc()[0].mean;
             expcon.push_back(met[k]);  break;}
             if(k==lmet) {for(int j=0;j<ntime;j++) fi>>aaa; getline(fi,aaa);}
                }
      return gpl;}
      
void Ldistr::setflcon(){ int i, ifound;//set flag if only concentration, not m0, is available
       for(int j=0;j<expcon.size();j++){
//       cout<<expcon[j]->getdescr()<<" "<<expcon[j]->getconc()[0].mean<<'\n';
        for(i=0;i<expm0.size();i++){
       ifound=expcon[j]->getdescr().find(expm0[i]->getdescr())+1; if(ifound) break; }
       if(!ifound) {expcon[j]->rizeflcon(); //cout<<expcon[j]->getdescr()<<" flcon=1\n";
       }
//         else cout<<expcon[j]->getdescr()<<" flcon=0\n";
       }
}
       
int main( int argc, char *argv[] ){
   cout<<"Nn="<<horse.getN()<<endl;
   double tmp,xi0; 
   int itmp; bool check;
         
   for(int i=0;i<argc;i++) cout<<argv[i]<<" "; cout<<argc<<":args\n"; 
       
     string sfpar(argv[3]), sout(argv[4]), sflmain(argv[5]), sflcomp(argv[6]), name;
      ifstream fi(argv[2]);  int ntr(0),Nfi=atoi(argv[7]);
     Problem.read(sfpar);    //read parameters
     Problem.setodir(sout,sflmain,sflcomp); //set output directory

      if(*Problem.getodir()=="glut") ntr=1; 
      ifn=Problem.setnumofi(); //number of parameter files
//order parameter files by increasing of χ2:
  if(argv[8][0]=='S')  { cout<<"statistics, "<<argv[4]<<endl; Problem.stat(ifn-1); return 0; }
  else if ((argc>3)&&(argv[3][0]=='x')) {chekxi(argv[1]); return 0; } // check χ2
   else{
     cout.precision(3);
        cout<<"** readexp**"<<'\n';
        horse.readExp(argv[1],ntr);            // read experimental data 
        cout<<"** readexp**"<<'\n';
     string gpl=horse.read_con(fi);  // read concentrations
        horse.setflcon();
      int m0len=horse.wrim0ex("exm0");                // set experimental data for figure
      int conlen=horse.wriconex("excon");                // set experimental data for figure
     for(int i=0;i<nmet;i++) {xinit1[i]=xx[i]; xinit2[i]=xx[i];}//copy initial values
	try{ cout<<'\n';  ts=clock(); 
       Problem.shownx(nmet,xx);    cout<<'\n';      // print concentrations on screen
  sol0=solve();//horse.integrbs();                 //solve ODEs for isotopomers
       Problem.shownx(nmet,xx);  // print concentrations on screen
       xi0=get<0>(sol0);
     cout<<" Σxi²="<<xi0 <<"\n";        //final results
     cout<<foc; 
        ofstream kkin("kinGlc"); kkin<<kin<<endl; kkin.close();
         kkin.open("kincon"); kkin<<kinc<<endl; kkin.close();
	} catch( char const* str ){cout << "exception: "<< str <<endl;}
        int ifn0=250002; Problem.write(sol0,ifn0,0);
		Analis analis;
		analis.setx00(get<0>(sol0),tf,suxx);
		Problem.storeVms(nrea,analis.nv1);
		Problem.storeVms(nrea,analis.nv2);
//	Problem.cont(121,0.0001, 0.015);
//chekxi(1,33);
          int sys=system(gpl.c_str());//gnuplot -e 'var=value' script.gp
		srand(time(NULL));
   Problem.setfnfin(ifn+Nfi);
     try{  switch(argv[8][0]){
   case 'N':
           analis.grad();  break;
   case 'F':
           for(int i=0;i<5;i++){analis.coord(0.2,1.07,sfpar);
             analis.grdesc(1.05); 
           } break;
   case 'G':
           analis.descent(1.05);
//           analis.grdesc(1.05);
           break;
   case 'C':
           analis.confidence(1.15,1.07); break;
   default: 
        cout << "single simulation finished\n" ; 
            }
  } catch(const invalid_argument&) { cout<<Nfi<<" files saved!\n"; }
  }
return (int)xi0;}

