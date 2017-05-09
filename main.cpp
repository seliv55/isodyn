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
  int  ifn=0; const int Nn= horse.getN(), neq=870;
   double dlt, dif, suxx, mader;
//   const double flf=60.*1000.*5./7.;
    time_t ts,tf; char *fex1, *fex2;
  extern string foc, kin0, kin, kinflx;
    
double solve(){
//solving the kinetic model:
   ts=clock(); 
 for(int i=0;i<nmet;i++) xx[i]=xinit1[i];
    tsolve(37000.); mader=Problem.dermax();
    // Problem.jacobian(xx);  tsolve(37000.); for(int i=0;i<numx;i++) xinit1[i]=xx[i];
      double xi=horse.integrbs();
//horse.readExp(fex2);for(int i=0;i<numx;i++) xx[i]=xinit1[i];
//      tsolve(37000.); 
//      xi += horse.integrbs();
  tf=clock()-ts;for(int i=0;i<nmet;i++) xx[i]=xinit1[i];
   suxx=horse.consum(); // suxx += dif;
    dif = suxx;
   return (xi);}
   
inline void chekxi(char *efi){
   stringstream fn; int itmp, iin, ifin;  double tmp,xi0;
    cout<<"first file: "; cin>>iin;
    cout<<"\nlast file: "; cin>>ifin; 
     for(int i=iin;i<ifin;i++) {fn<<*Problem.getodir()<<i; cout<<i<<": "<<endl;
      Problem.read(itmp,tmp,fn.str().c_str()); cout<<fn.str().c_str()<<endl;
      horse.readExp(efi); for(int ii=0;ii<nmet;ii++) cout<<(xinit1[ii]=xx[ii])<<endl;
         xi0=solve();            
   Problem.write(tf,i,xi0,suxx,0);      } /**/
  }

int main( int argc, char *argv[] ){
   double tmp,xi0;ofstream kkin("kinxx"); 
//   int a[][3]={{1,2,3},{4,5,6}}}; cout<<"a11="<<a[1][1]<<endl;
   int itmp; bool check;
         fex1=argv[1]; fex2=fex1;
    Problem.setodir(argv[2]);  ifn=Problem.setnumofi(); 
     
  if((argc>3)&&(argv[3][0]=='s'))  { cout<<argv[3][0]<<endl; Problem.stat(ifn-1); return 0; } //order parameter files by increasing of χ2
  else if ((argc>3)&&(argv[3][0]=='x')) {chekxi(argv[1]); return 0; } // check χ2
   else{
     cout.precision(3);
     xi0=Problem.read(itmp,tmp,argv[2]);    //read parameters
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
        int sys=system("plkin.p");//gnuplot -e 'var=value' script.gp
  xi0=solve();//horse.integrbs();                 //solve ODEs for isotopomers
//	horse.flback();                 //print Σisotopomers on screen
     cout<<" Σxi²="<<xi0 <<"\n";        //final results
     cout<<setw(9)<<"* "<<"*Metab  *  init  Final : exper -> xi²"<<"\n" <<foc; 
//        kkin.open("kinGlc"); kkin<<kin<<endl; kkin.close();
//        kkin.open("kinflx"); kkin<<kinflx<<endl; kkin.close();
        
//        horse.readExp(fex2);
//      for(int i=0;i<numx;i++) xx[i]=xinit1[i];
//      tsolve(37000.);
//   xi0+=horse.integrbs();
//  xi0+=horse.ddisolve();
       tf=clock()-ts;
	} catch( char const* str ){cout << "exception: "<< str <<endl;}
//     cout<<foc; kkin.open("kinetics"); kkin<<kin<<endl; kkin.close();
//      suxx=horse.label(); 
		Analis analis;
		analis.setx00(xi0,tf,suxx);
//           for(int i=0;i<nmet;i++) xx[i]=xinit1[i];
                int ifn0=250002; Problem.write(tf,ifn0,xi0,suxx,0);
//	Problem.cont(121,0.0001, 0.015);
//chekxi(1,33);
//          int sys=system("gnuplot xplt.p");//gnuplot -e 'var=value' script.gp
		srand(time(NULL));
 if (argc>3) {  stringstream stx(argv[3]); int ia; stx>>ia; cout<<ia<<endl; Problem.setfnfin(ifn+ia);
   
// 	cout<<"parameter set="<<ia<<endl;
               analis.confidence(1.15,1.07);
	         Problem.stat(ifn-1);
                  xi0=Problem.read(itmp,tmp,"1");
                   for(int i=0;i<numx;i++) xinit1[i]=xx[i];
                  }
//               analis.sensitiv(tmax);
//               analis.swarm(tmax,111);
//               analis.grad(tmax);
}}

