//---------------------------------------------------------------------------
#ifndef NVH
#define NVH
#include <fstream>
#include <sstream>
#include <iomanip>
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
  
	void setVm(double a) {par[0] = a;}
	double chanVm(double ff) {double old=par[0]; par[0] *= ff; return old;}
	double *getpar(){return par;}
	std::string& getname(){return name;}
		Reapar() {}
	virtual ~Reapar(){delete[] par; delete[] spar; }
};

class Parray {
protected:
	double ft3, fh6, xi, tan, tnad, xthf;
	 static int par[];
        TK tk;
        TA ta;
        Aldolase aldolase;
	void flfor(double y[]);
public:
 static std::string namef[], namex[];
	static Reapar rea[];
       double sumx();
       void result() const;
	void storeVms(int len,double a[]){for(int i=0;i<len;i++) a[i]=rea[i].v();};
	void restoreVm(int len,double a[]) {for(int i=0;i<len;i++) rea[i].setVm(a[i]);}
	static int* getpar(){return &par[0];}
      bool belong(int ip) const;
        void init();
        void fin(double y[]);
  Parray( ): tan(1.1),tnad(1.){ }
  virtual ~Parray(){  }
};
class Fit: public Parray{
	int i68, i90, i95, i99, fnfin;
	std::string outdir;
    public:	
      void stat(const int NP );
      double read(int &t,double &c, std::string fn);
      void readst( int* );
      void setfnfin(int ia){fnfin=ia;}
      void wstorefl (const char fn1[], int numpar, const double** m,std::string n[]);
      void write (time_t tf, int& fn,const  double xi0, const double xm,bool flg=true) const ;
      void perturb(const double f1);
      void cont(const int,const  double,const double);
      void setpar(int p[]){ for(int i=1;;i++) {par[i]=p[i-1]; if(par[i]<0) {par[0]=i; break;}}}
      void fitc(double dc,double dm,int iin,int iout);
      int* getFitPar(){return &par[0];}
      void setsst();
      double jacobian(double *y);
      double eigen(double dx[]);
      void chpar(int ipar,int ovar,double factor,int);
      void f(const double *y,double *dydx);
      void ff(const double *y,double *dydx);
 
    void setodir(char *filo){ std::string infi(filo);
         int pos=infi.find_last_of('/');   outdir=infi.substr(0,pos+1); std::cout<<outdir<<std::endl; }
    std::string* getodir(){return &outdir;}
         
    int setnumofi()const{ for(int i=1;;i++) {std::stringstream fn;
       fn<<outdir<<i;   std::ifstream checkfi(fn.str().c_str()); checkfi.close(); 
	 if(!checkfi.good())  return (i);}    }

     int getListFit(int list[]) const {      int k=1;
	   while (par[k]+1) {list[k-1] = par[k]; k++;} return (k-1);}
	   
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


