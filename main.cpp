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
  int  ifn=0;  int Nn;
   double dlt, dif, suxx, mader;
//   const double flf=60.*1000.*5./7.;
    time_t ts,tf; char *fex1, *fex2;
  extern string foc, kin0, kin, kinflx;
    
tuple<double,double,time_t> solve(){
//solving the kinetic model:
   ts=clock(); 
 for(int i=0;i<numx;i++) xx[i]=xinit1[i]; //take initial values
    tsolve(37000.); mader=Problem.dermax();
//  for(int i=0;i<10;i++)   ddsolve(3700.);
    // Problem.jacobian(xx);  tsolve(37000.); for(int i=0;i<numx;i++) xinit1[i]=xx[i];
      double xi=horse.integrbs();
//horse.readExp(fex2);for(int i=0;i<numx;i++) xx[i]=xinit1[i];
//      tsolve(37000.); 
//      xi += horse.integrbs();
  tf=clock()-ts;
   suxx=horse.consum(); // suxx += dif;
   for(int i=0;i<numx;i++) xx[i]=xinit1[i]; //set iv for a case of saving them
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
    
string Ldistr::read_con(ifstream& fi, string& arg1,int& nfi){
// read par-file, output dir, flux conf. inter: main and to compare
 string aaa, spar, sout, sflmain, sflcomp, gpl; 
 fi>>aaa>>spar>>aaa>>sout>>aaa>>nfi>>aaa>>sflmain>>aaa>>sflcomp;
 getline(fi,aaa); 
 getline(fi,gpl); cout<<gpl<<endl;
     Problem.read(spar.c_str());    //read parameters
// extract cell_type/conditions, numbers: of substrates, time points
   string cell=arg1.substr(arg1.find_last_of('/')+1); cout<<"cell="<<cell<<endl;
     Problem.setodir(sout,sflmain,sflcomp); //set output directory
     int isu;
      fi>>isu>>aaa>>ntime>>aaa; 
// localize cell_type/conditions
      while(!fi.eof()){getline(fi,aaa);  if((aaa.find(cell)+1)) break;}
      if(!fi.eof()){
      for(int i=0;i<ntime;i++) {double ta; fi>>ta>>aaa; texcon.push_back(ta); }
      for(int i=0;i<isu;i++){
       fi>>aaa; int k;
        for(k=0;k<=lmet;k++) if(aaa.find(met[k]->getdescr())+1) {
           for(int j=0;j<ntime;j++) {
       double cc; fi>>cc; met[k]->setconc(cc,j); xx[k]=met[k]->getconc()[0].mean;}
             expcon.push_back(met[k]);  break;}
             if(k==lmet) for(int j=0;j<ntime;j++) fi>>aaa;
                }
   }   return gpl;}
void Ldistr::setflcon(){ int i, ifound;
       for(int j=0;j<expcon.size();j++){
        for(i=0;i<expm0.size();i++){
       ifound=expcon[j]->getdescr().find(expm0[i]->getdescr())+1; if(ifound) break; }
       if(!ifound) {expcon[j]->rizeflcon(); cout<<expcon[j]->getdescr()<<" flcon=1\n";}
         else cout<<expcon[j]->getdescr()<<" flcon=0\n";
       }
}

int main( int argc, char *argv[] ){
// run: ./isodyn.out  experim_MID-lile concentr_and_param_info_file [sx_NumOfFileSavedInFitting]
// 's'-statistics; 'x'- check χ²
 cout<<"Nn="<<horse.getN()<<endl;
   double tmp,xi0;ofstream kkin("kinxx"); 
//   int a[][3]={{1,2,3},{4,5,6}}}; cout<<"a11="<<a[1][1]<<endl;
   int itmp; bool check;
         fex1=argv[1]; fex2=fex1;
         
         
     string arg1(argv[1]),name; ifstream fi(argv[2]);  int ntr(0),Nfi;
     string gpl=horse.read_con(fi,arg1,Nfi); if(*Problem.getodir()=="glut/") ntr=1; 
      ifn=Problem.setnumofi(); //number of parameter files
      cout<<ntr<<"=ntr "<<ifn<<"=ifn\n";
     
  if((argc>3)&&(argv[3][0]=='s'))  { cout<<argv[3][0]<<endl; Problem.stat(ifn-1); return 0; } //order parameter files by increasing of χ2
  else if ((argc>3)&&(argv[3][0]=='x')) {chekxi(argv[1]); return 0; } // check χ2
   else{
     cout.precision(3);
//     sol0=Problem.read(argv[2]);    //read parameters
        horse.readExp(argv[1],ntr);            // read experimental data 
        horse.setflcon();
      int m0len=horse.wrim0ex("exm0");                // set experimental data for figure
      int conlen=horse.wriconex("excon");                // set experimental data for figure
     for(int i=0;i<numx;i++) {xinit1[i]=xx[i]; xinit2[i]=xx[i];}//copy initial values
	try{   ts=clock();
    tsolve(37000.);                     //solve ODEs for total concentrations
//     mader=Problem.dermax(); Problem.jacobian(xx);  tsolve(37000.);
//    tsolve(35000.); Problem.jacobian(xx); 
//         cout<<"aldf="<<flx[aldf]<<"; aldr="<<flx[aldf+1]<<endl;
//       Problem.chpar(28,noa,0.99,27);
//    tsolve(35000.);  for(int i=0;i<numx;i++)  xinit1[i]=xx[i];
       Problem.shownx(numx,xx);         // print concentrations on screen
     kkin<<kin0;  kkin.close();         //save concentration dynamics to "kinxx"
        int sys=system("gnuplot plkin.p");//gnuplot -e 'var=value' script.gp
//  xi0=horse.integrbs();                 //solve ODEs for isotopomers
  sol0=solve();//horse.integrbs();                 //solve ODEs for isotopomers
       Problem.shownx(numx,xx);         // print concentrations on screen
     cout<<" Σxi²="<<get<0>(sol0) <<"\n";        //final results
     cout<<setw(9)<<"*"<<"*Metab   *   init    Final : "<<setw(17)<<"exper -> xi²\n" <<foc; 
        kkin.open("kinGlc"); kkin<<kin<<endl; kkin.close();
         kkin.open("kincon"); kkin<<kinc<<endl; kkin.close();
//        kkin.open("kinflx"); kkin<<kinflx<<endl; kkin.close();
        
//        horse.readExp(fex2);
//      for(int i=0;i<numx;i++) xx[i]=xinit1[i];
//      tsolve(37000.);
//   xi0+=horse.integrbs();
//  xi0+=horse.ddisolve();
//       tf=clock()-ts;
	} catch( char const* str ){cout << "exception: "<< str <<endl;}
//     cout<<foc; kkin.open("kinetics"); kkin<<kin<<endl; kkin.close();
//      suxx=horse.label(); 
		Analis analis;
		analis.setx00(xi0,tf,suxx);
//           for(int i=0;i<nmet;i++) xx[i]=xinit1[i];
                int ifn0=250002; Problem.write(sol0,ifn0,0);
//	Problem.cont(121,0.0001, 0.015);
//chekxi(1,33);
          int sys=system(gpl.c_str());//gnuplot -e 'var=value' script.gp
		srand(time(NULL));
                cout<<"niso="<<horse.getmicon()<<'\n';
               analis.grad(1);
 if (argc>3) {  Problem.setfnfin(ifn+Nfi);
  if(argv[3][0]=='f') analis.coord(0.03,1.07); else {
     try{ analis.confidence(1.15,1.07);} catch(const invalid_argument&){cout<<Nfi<<" files saved!\n"; return 0;}
	         Problem.stat(ifn-1);
                  sol0=Problem.read("1");
                   for(int i=0;i<numx;i++) xinit1[i]=xx[i];
                  }
 }
//               analis.sensitiv(tmax);
//               analis.swarm(tmax,111);
}
return 0;}

