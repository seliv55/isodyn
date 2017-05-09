//---------------------------------------------------------------------------
#include <iostream>
#include <string>
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
  int  ifn=0; const int Nn= horse.getN(), neq=870;
   double dlt, dif, suxx, mader;
//   const double flf=60.*1000.*5./7.;
    time_t ts,tf; char *fex1, *fex2;
  extern string foc, kin0, kin, kinflx;
    
tuple<double,double,time_t> solve(){
//solving the kinetic model:
   ts=clock(); 
 for(int i=0;i<nmet;i++) xx[i]=xinit1[i]; //take initial values
    tsolve(37000.); mader=Problem.dermax();
//  for(int i=0;i<10;i++)   ddsolve(3700.);
    // Problem.jacobian(xx);  tsolve(37000.); for(int i=0;i<numx;i++) xinit1[i]=xx[i];
      double xi=horse.integrbs();
//horse.readExp(fex2);for(int i=0;i<numx;i++) xx[i]=xinit1[i];
//      tsolve(37000.); 
//      xi += horse.integrbs();
  tf=clock()-ts;
   suxx=horse.consum(); // suxx += dif;
   for(int i=0;i<nmet;i++) xx[i]=xinit1[i]; //set iv for a case of saving them
   return make_tuple(xi,suxx,tf);}
   
inline void chekxi(char *efi){
   stringstream fn; int iin, ifin;  double xi0;
    cout<<"first file: "; cin>>iin;
    cout<<"\nlast file: "; cin>>ifin; 
     for(int i=iin;i<ifin;i++) {fn<<*Problem.getodir()<<i; cout<<i<<": "<<endl;
      Problem.read(fn.str().c_str()); cout<<fn.str().c_str()<<endl;
      horse.readExp(efi); for(int ii=0;ii<nmet;ii++) cout<<(xinit1[ii]=xx[ii])<<endl;
        tuple<double,double,time_t> sol=solve();            
   Problem.write(sol,i,0);      } /**/
  }

    tuple<double,double,time_t> sol0;
    
//void Ldistr::read_con(ifstream& fi, string& arg1){int kmet;
//       string cell=arg1.substr(arg1.find_last_of('/')+1); cout<<"cell="<<cell<<endl;
//      string aaa; int isu, ntime(2);
//      fi>>isu>>aaa>>ntime>>aaa; double conc;
//      while(!fi.eof()){getline(fi,aaa);  if((aaa.find(cell)+1)) break;}
//      if(!fi.eof()){
//      for(int i=0;i<ntime;i++) {double ta; fi>>ta>>aaa; texcon.push_back(ta); cout<<"tex="<<texcon[i]<<endl;}
//      for(int i=0;i<isu;i++){
//       fi>>aaa; cout<<aaa<<"  "; for(int k=0;k<10;k++)
//              if(aaa.find(met[k]->getdescr())+1) {kmet=k;cout<<met[k]->getdescr()<<endl; break;}
//        for(int j=0;j<ntime;j++) {fi>>conc; met[kmet]->setconc(conc,j);}
//      }
//   }   }

int main( int argc, char *argv[] ){
   double tmp,xi0;ofstream kkin("kinxx"); 
//   int a[][3]={{1,2,3},{4,5,6}}}; cout<<"a11="<<a[1][1]<<endl;
   int itmp; bool check;
         fex1=argv[1]; fex2=fex1;
         
//     string arg1(argv[1]),name; ifstream fi("xconc");
//     horse.read_con(fi,arg1);
         
    string indir= Problem.setodir(argv[2]); //set output directory
      ifn=Problem.setnumofi(); //number of parameter files
     
  if((argc>3)&&(argv[3][0]=='s'))  { cout<<argv[3][0]<<endl; Problem.stat(ifn-1); return 0; } //order parameter files by increasing of χ2
  else if ((argc>3)&&(argv[3][0]=='x')) {chekxi(argv[1]); return 0; } // check χ2
   else{
     cout.precision(3);
     sol0=Problem.read(argv[2]);    //read parameters
//      horse.setfige();                // set experimental data for figure
        horse.readExp(argv[1]);            // read experimental data 
     for(int i=0;i<nmet;i++) {xinit1[i]=xx[i]; xinit2[i]=xx[i];}//copy initial values
	try{   ts=clock();
    tsolve(37000.);                     //solve ODEs for total concentrations
//     mader=Problem.dermax(); Problem.jacobian(xx);  tsolve(37000.);
//    tsolve(35000.); Problem.jacobian(xx); 
//         cout<<"aldf="<<flx[aldf]<<"; aldr="<<flx[aldf+1]<<endl;
//       Problem.chpar(28,noa,0.99,27);
//    tsolve(35000.);  for(int i=0;i<numx;i++)  xinit1[i]=xx[i];
     Problem.shownx(numx,xx);         // print concentrations on screen
     kkin<<kin0;  kkin.close();         //save concentration dynamics to "kinxx"
     int sys=system("plkin.p "+indir+"sconc.png"); //gnuplot -e 'var=value' script.gp
//  xi0=horse.integrbs();                          //solve ODEs for isotopomers
     sol0=solve();//horse.integrbs();              //solve ODEs for isotopomers
     Problem.shownx(numx,xx);         // print concentrations on screen
     cout<<" Σxi²="<<get<0>(sol0) <<"\n";        //final results
     cout<<setw(9)<<"*"<<"*Metab   *   init    Final : "<<setw(17)<<"exper -> xi²\n" <<foc; 
//        kkin.open("kinGlc"); kkin<<kin<<endl; kkin.close();
//        kkin.open("kinflx"); kkin<<kinflx<<endl; kkin.close();
        
//        horse.readExp(fex2);
//      for(int i=0;i<numx;i++) xx[i]=xinit1[i];
//      tsolve(37000.);
//   xi0+=horse.integrbs();
//  xi0+=horse.ddisolve();
//       tf=clock()-ts;
//cout<<comm<<endl;
	} catch( char const* str ){cout << "exception: "<< str <<endl;}
//     cout<<foc; kkin.open("kinetics"); kkin<<kin<<endl; kkin.close();
//      suxx=horse.label(); 
		Analis analis;
		analis.setx00(xi0,tf,suxx);
//           for(int i=0;i<nmet;i++) xx[i]=xinit1[i];
                int ifn0=250002; Problem.write(sol0,ifn0,0);
//	Problem.cont(121,0.0001, 0.015);
//chekxi(1,33);
//          int sys=system("gnuplot xplt.p");//gnuplot -e 'var=value' script.gp
		srand(time(NULL));
 if (argc>3) //if(argv[3][0]=='g') analis.grad(1000); else 
{  stringstream stx(argv[3]); int ia; stx>>ia; cout<<ia<<endl; Problem.setfnfin(ifn+ia);
   
// 	cout<<"parameter set="<<ia<<endl;
           try{ analis.confidence(1.15,1.07);} catch(const invalid_argument&){cout<<ia<<" files saved!\n"; return 0;}
	         Problem.stat(ifn-1);
                  sol0=Problem.read("1");
                   for(int i=0;i<numx;i++) xinit1[i]=xx[i];
                  }
//               analis.sensitiv(tmax);
//               analis.swarm(tmax,111);
//               analis.grad(tmax);
}
return 0;}

