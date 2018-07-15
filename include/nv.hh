//---------------------------------------------------------------------------
#ifndef NVH
#define NVH
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <tuple>
class Reapar{
   std::string name, *spar;
    int npar; double *par;
  public:
	double v(){return par[0];}
	double v(double x){return par[0]*x/(par[1]+x);}
	double v(double x,double y){return par[0]*x/(par[1]+x)*y/(par[2]+y);}
	double v(double x,double y,double z){return par[0]*x/(par[1]+x)*y/(par[2]+y)*z/(par[3]+z);}
	double ving(double x,double y,double z){return par[0]*x/(par[1]+x+y/par[2]+z);}
	double vin(double x,double y){return par[0]*x/(par[1]+x)*(par[2]/y+par[3]);}
	
  void read(std::ifstream& fi,int& flg){fi>>name>>name>>npar;
   if(!flg) {par=new double[npar]; spar=new std::string[npar];}
     for(int i=0;i<npar;i++) fi>>this->spar[i]>>this->par[i]; }
     
 void write(std::ofstream& fo,int ir){fo<<ir<<std::setw(9)<<name<<std::setw(4)<<npar;
  for(int i=0;i<npar;i++) fo<<"  "<<this->spar[i]<<" "<<this->par[i]; fo<<"\n";}
  
	void setVm(double a,int np=0) {par[np] = a;}
	double chanVm(double ff,int np=0) {double old=par[np]; par[np] *= ff; return old;}
	double *getpar(){return par;}
	std::string& getname(){return name;}
	int genpar(){return npar;}
		Reapar() {}
	virtual ~Reapar(){delete[] par; delete[] spar; }
};

class Parray {
protected:
	double ft3, fh6, xi, tan, tnad, xthf;
	std::vector<int> par;
        TK tk;
        TA ta;
        Aldolase aldolase;
	void flfor(double y[]);
public:
// static std::string namef[], namex[];
 static std::string fid[],fname[],fschem[], namex[];
	static Reapar rea[];
       double sumx();
       void result() const;
	void storeVms(int len,double a[]){for(int i=0;i<len;i++) a[i]=rea[i].v();};
	void restoreVm(int len,double a[]) {for(int i=0;i<len;i++) rea[i].setVm(a[i]);}
	 int* getpar(){return &par[0];}
      bool belong(int ip) const;
        void init();
        void fin(double y[]);
        void rnames(std::ifstream& fi);
  Parray( ): tan(1.1),tnad(1.){ }
  virtual ~Parray(){  }
};
class Fit: public Parray{
	int i68, i90, i95, i99, fnfin;
	std::string outdir, flmain,flcomp, spar;
    public:	
      void stat(const int NP );
      std::tuple<double,double,time_t> read(std::string fn);
      void readst( int* );
      void setfnfin(int ia){fnfin=ia;}
      void wstorefl (int numpar, const double** m,std::string n[]);
      void write (std::tuple<double,double,time_t> sol, int& ifn,bool flg=true) ;
      void perturb(const double f1);
      void cont(const int,const  double,const double);
      void fitc(double dc,double dm,int iin,int iout);
      std::vector<int> getFitPar(){return par;}
      void setsst();
      double jacobian(double *y);
      double eigen(double dx[]);
      void chpar(int ipar,int ovar,double factor,int);
      void f(const double *y,double *dydx);
      void ff(const double *y,double *dydx);
      int getparsize(){for(unsigned i=0;i<par.size();i++) std::cout<<par[i]<<" "; std::cout<<"\n"; return par.size();}
 
    void setodir(std::string outd,std::string flm,std::string flc){ outdir=outd; flmain=flm; flcomp=flc;}
         
    std::string* getodir(){return &outdir;}
    std::string* getflm(){return &flmain;}
    std::string* getflc(){return &flcomp;}
         
    int setnumofi()const{ for(int i=1;;i++) {std::stringstream fn;
       fn<<outdir<<i;   std::ifstream checkfi(fn.str().c_str()); checkfi.close(); 
	 if(!checkfi.good())  return (i);}    }

     void getListFit(int list[]) const {  
	   for (int i=0;i<par.size();i++) list[i] = par[i]; }
	   
        void shownx(int nx,double dx[]){for(int i=0;i<nx;i++) std::cout<<namex[i]<<":"<<dx[i]<<"; ";
                                             std::cout<<std::endl;}
        double dermax();
	  Fit( ): Parray() {}
	   ~Fit(){}
};
extern Fit Problem;
//extern double fluxes[],dt,xx[];
//---------------------------------------------------------------------------
#endif


