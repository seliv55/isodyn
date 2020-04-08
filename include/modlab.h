//---------------------------------------------------------------------------

#ifndef labH
#define labH
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <set>
#include <cmath>
const int tt=5;

class data {
public:
	 double mean, sd;
        data(double a=0.):mean(a),sd(a){}
        ~data(){}
};

class Iso {
 protected:
  int niso; double ttime;  std::string name; data *mid;
 public:
    void showmid(std::string s=""){for(int i=0;i<niso;i++) std::cout<<mid[i].mean<<" "; std::cout<<name+s<<" time="<<ttime<<'\n';}
    void showsd(std::string s=""){for(int i=0;i<niso;i++) std::cout<<mid[i].sd<<" "; std::cout<<name+s<<'\n';}
    void calmesd(std::vector<Iso>& linj){
       int len=linj.size();
         for(int ic=0;ic<niso;ic++) { mid[ic].mean=0.; mid[ic].sd=0.;
          for(int i=0;i<len;i++) mid[ic].mean+=linj[i].mid[ic].mean; 
          mid[ic].mean/=(double)len; }
         if(len>1) for(int ic=0;ic<niso;ic++){ 
        for(int i=0;i<len;i++){ double a=mid[ic].mean-linj[i].mid[ic].mean; mid[ic].sd+=a*a; }
        mid[ic].sd/=(double)(len-1); mid[ic].sd=sqrt(mid[ic].sd); }
           }
    void setname(std::string s){name=s;}
    void setmid(int mi, double d){mid[mi].mean=d;}
    std::string* getname(){return &name;}
    int getniso(){return niso;}
    double gett(){return ttime;}
    data* getmid(){return mid;}
    
  Iso(int n, double tti, std::string s){niso=n; ttime=tti; name=s;
   mid=new data[n]; mid[0].mean=1.; mid[0].sd=0.1;
   for(int i=1;i<niso;i++){ mid[i].mean=0.; mid[i].sd=0.1; }}
//  Iso(const Iso &obj){mid=new data[niso];
//        for(int i=0;i<niso;i++){ mid[i].mean=obj.mid[i].mean; mid[i].sd=obj.mid[i].sd; }}
  virtual ~Iso(){delete[] mid;}
};

class Tracer:public Iso {
 protected:
  int nmark,nmi;
 public:
  double fract;
  void setrac(){mid[0].mean=(100.-fract); mid[nmi].mean=fract;}
 Tracer(int n, double tti,std::string s,int nc,double abun,int nlab):Iso(n,tti,s){nmark=nc; fract=abun;
   nmi=nlab;}
 ~Tracer(){}
};

class Exper{
        int mi, frbeg;
        data *exper[tt];
        double *calc,*calstor[tt],xi[tt], xicon[tt],kinm0[tt],time[tt];
        std::string descr;
 public:
void setex0(){
     exper[0][0].mean=1.; exper[0][0].sd=0.01;
      for(int i=1;i<=mi;i++){ exper[0][i].mean=0.; exper[0][i].sd=0.01;  }
}

void setime(double tt){}

void setrav(int nt1,int nt2){
      for(int i=0;i<=mi;i++){ exper[nt1][i].mean=exper[nt2][i].mean; exper[nt1][i].sd=exper[nt2][i].sd;  }
}

void rex(std::ifstream& ifi, int nt){// reads experimental data
      std::string aaa; double sdd;
      ifi>>aaa;
/*      std::cout<<" rex "<<aaa<<" ";*/
    for(int i=0;i<=mi;i++){ ifi>>exper[nt][i].mean;
     std::cout<<exper[nt][i].mean<<" ";
     }  std::cout<<'\n';
      getline(ifi,aaa);
      ifi>>aaa;
/*      std::cout<<"\n rex "<<aaa<<" ";*/
       for(int i=0;i<=mi;i++){ifi>>sdd; if(sdd<0.01)sdd=0.01; exper[nt][i].sd=sdd;
/*        std::cout<<exper[nt][i].sd<<" ";*/
        } //std::cout<<'\n';
/*       getline(ifi,aaa);*/
}

double percent(const double *iso,int len){
		double tot=0.; int nc,sum;
		for(int i=0;i<mi+2;i++) calc[i]=0.;
   for(nc=0;nc<len;nc++) { int i=((nc>>frbeg)&((1<<mi)-1)); sum=0;
	while(i>0){sum += (i & 1);i = (i >> 1);}     calc[sum] += iso[nc];
		          }//end for
	for(int i=0;i<=mi;i++) tot += calc[i];
//	for(int i=0;i<=mi;i++) calc[i] /= tot;
return calc[mi+1]=tot;}

int stormi(double dis[],int nt) {
     for(int j=1;j<nt;j++) for(int i=0;i<mi;i++) dis[i+(j-1)*mi]=this->calstor[j][i];
return mi*(nt-1);}
 
data * getexper(int nt){return &exper[nt][0];}

std::string& getdescr(){return descr;} 

void skladm0(int itime){ kinm0[itime]=calc[0]/calc[mi+1];
   std::cout<<descr<<" time="<<itime<<" conc="<<calc[mi+1]<<'\n';
     for(int i=0;i<=mi;i++) std::cout<<calc[i]/calc[mi+1]<<" ";
                        std::cout<<'\n'; 
                        }

int getmi(){ return mi;}
 
 void showsum(double calc2[],std::ostringstream& fo){
   fo<<" "<<std::setw(7)<<(this->calc[0]+calc2[0])/(this->calc[mi+1]+calc2[mi+1]);
 }
 void showm0(std::ostringstream& fo){
   fo<<" "<<std::setw(7)<<this->calc[0]/this->calc[mi+1];
 }

double chisq(int nt) {  xi[nt]=0.; double a; 
  for(int k=0;k<mi;k++) if(exper[nt][k].sd>5.e-3){
   calstor[nt][k]=calc[k]/calc[mi+1];
    a=(exper[nt][k].mean - calstor[nt][k])/(exper[nt][k].sd); xi[nt] += a*a;}
       else calstor[nt][k]=0.;
         return xi[nt];}
  
double chisqsum(double s2[],int nt,int nl) { xi[nt]=0; 
	double a, b=this->calc[mi+1]+s2[mi+1];
  if(exper[nt][0].sd>1.e-5) for(int k=0;k<nl;k++) {
     a=(exper[nt][k].mean-(this->calc[k]+s2[k])/b)/(exper[nt][k].sd);
      xi[nt] += a*a; }
  return xi[nt];}

void wrikinm0(std::ostringstream& so, int nt) {
     so.precision(3);
      so<<"\n"<<std::setw(12)<<descr<<"_m0: ";
    for(int i=0;i<nt;i++)  so<<std::setw(9)<<this->kinm0[i];
    so<<" : "<<std::setw(7)<<exper[nt-1][0].mean<<" -> "<<xi[nt-1];
 }
 void shexper(int nt){ std::cout<<getdescr()<<" ";
   for(int j=0;j<nt;j++){
     for(int i=0;i<=mi;i++) std::cout<<exper[j][i].mean<<" "; 
     std::cout<<'\n';} std::cout<<'\n';
}
double* getcalc(){return &calc[0];}
double shkin(int itp){return this->kinm0[itp];}

double dilut( double& dlt,const double *iso,const int len,const int nt) {
	     int i; double xi(0.0),a(0.0);
	   this->percent(iso,len);
	dlt = (exper[nt][0].mean-this->calc[0])/(1.-exper[nt][0].mean);
		if(dlt>0) {
		this->calc[0] += dlt;
		for(i=0;i<=mi;i++) this->calc[i] /= (1.+dlt);
		this->calc[mi+1] *= (1.+dlt);
		}
		for(i=0;i<(mi+2);i++) {    //          for(int i=1;i<repl;i++)
			a = (exper[nt][i].mean-this->calc[i])/(exper[nt][i].sd);
			xi += a*a;
		}
		return xi;
}

   Exper(int n,int beg,std::string simya=""): mi(n), frbeg(beg), descr(simya) { for(int i=0;i<tt;i++) {
      exper[i]=new data[mi+1]; calstor[i]=new double[mi+1];} exper[0][0].mean=1.; calc=new double[mi+2]; } //
   Exper(const Exper &obj) = default; //{//
//   for(int i=0;i<tt;i++) {      exper[i]=new data[mi+1]; calstor[i]=new double[mi+1];
//   for(int j=0;j<mi+2;j++){exper[i][j]=obj.exper[i][j]; calstor[i][j]=obj.calstor[i][j]; } }
//      exper[0][0].mean=1.; calc=new double[mi+2]; for(int j=0;j<mi+3;j++) calc[j]=obj.calc[j];} //
   virtual ~Exper(){};//for(int i=0;i<tt;i++) {delete[] exper[i]; delete[] calstor[i];} delete[] calc;}
};

class Metab{
 protected:
      const int N, len;
      bool flag;
      double fc,tcon, xicon[tt], kinc[tt];
  std::string descr;
        data conc[tt];
 public:
  std::vector<Exper>  exper;
  bool flcon;
	double *iso, *diso;
	int ny;

void wrikin(std::ostringstream& so, int nt){
   for(unsigned int i=0;i<exper.size();i++) exper[i].wrikinm0(so,nt);
 }
 
void percent(){ for(unsigned int j=0;j<this->exper.size();j++) tcon=this->exper[j].percent(iso,len);}

void readc(std::ifstream& fi,  int nt){ std::string aaa;
     fi>>aaa;
     for(int i=0;i<nt;i++) fi>>conc[i].mean;
     for(int i=0;i<nt;i++) {fi>>conc[i].sd; if(conc[i].sd<0.01) conc[i].sd=0.01;}
}

 void showcon(std::ostringstream& fo){ fo<<" "<<std::setw(7)<<this->sumt();}

void showm0(std::ostringstream& fo,double xfin) {
   for(unsigned int j=0;j<exper.size();j++) exper[j].showm0(fo);
}

void skladc(int itime){ kinc[itime]=tcon; }

void skladm0(int itime){ for(unsigned int j=0;j<(exper.size());j++) exper[j].skladm0(itime);}

void rizeflcon(){flcon=1;}
  
data * getconc(){return conc;}
 
void setconc(double a, int nt=0){ conc[nt].mean=a; conc[nt].sd=a*0.1; }
 
double xits(int its) { double xi=0;
/*for(unsigned int j=0;j<exper.size();j++)*/
 xi += exper[0].chisq(its);
return xi;}

double chicon(int nt) {
	double a;
     a=(conc[nt].mean - this->tcon)/(this->conc[nt].sd);
  return (xicon[nt]=a*a);}

void wrikinc(std::ostringstream& so, int nt) {
  so.precision(3);
    so<<"\n"<<std::setw(12)<<descr<<"_c: ";
    for(int i=0;i<nt;i++) so<<std::setw(9)<<this->kinc[i];
    so<<"** "<<std::setw(7)<<conc[nt-1].mean<<" -> "<<xicon[nt-1];
 }

int getN(){return N;}

std::string getdescr(){return descr;}
 
double sumt(){	double sum=0.; int i;
		for (i=0;i<len;i++) sum += iso[i];
return tcon=sum;	}

double input(Metab& s2,const double& vi,const double vo=0.){
		double x, sum=0.;
	for(int i=0;i<len;i++) {
		x=(vi*(iso[i])-vo*(s2.iso[i]));
		diso[i] -= x; s2.diso[i] += x; sum += x;		}
return sum;}

void output(const double& vi){	double x;
	for(int i=0;i<len;i++) {
		x=vi*iso[i];	diso[i] -= x;}
}
 
void set0(const double& s0){	iso[0]=s0;
		for(int i=1;i<len;i++) iso[i]=0.;
}

int getlen() const {return len;	}

double* getisot() {return iso;	}

void coriso(const double xkin,const double xiso){
 double ff=xkin/xiso;
  for(int i=0;i<len;i++) iso[i] *= ff;
}

double cutfirst(Metab& p,const double& v) {
	double dx, sum=0.;
	int i1;
	int lenp=((len>>1)-1);
	for(int i=0;i<this->len;i++) {i1 = (i & lenp);
		dx = v*(this->iso[i]);
	this->diso[i] -= dx; p.diso[i1] += dx; sum += dx;}
return sum;}
	

void volume(double factor) {for (int i=0; i<len;i++) diso[i] *= factor;	}

void carb (Metab& pr,const double& v){
	int ipr; double dx;
	for(int i=0;i<this->len;i++) {
		dx = v*(this->iso[i]);
		this->diso[i] -= dx;
	ipr=(i<<1);	pr.diso[ipr] += dx;}
}

void decarb (Metab& pr,const double& v) {
	int i,ipr; double dx;
	for(i=0;i<this->len;i++) {
		dx = v*(this->iso[i]);
		this->diso[i] -= dx;
		ipr=(i>>1);	pr.diso[ipr] += dx;}
}

void icdh (Metab& pr,const double& v) {
	int i,ipr; double dx;
	for(i=0;i<this->len;i++) {
		dx = v*(this->iso[i]);
		this->diso[i] -= dx;
	ipr=((i&7)|((i>>1)&24)); pr.diso[ipr] += dx;}
}

void icdhr (Metab& pr,const double& v) {
	int i,ipr; double dx;
	for(i=0;i<this->len;i++) {
		dx = v*(this->iso[i]);
		this->diso[i] -= dx;
	ipr=((i&7)|((i&24)<<1)); pr.diso[ipr] += dx;}
}

double liase (Metab& pr,Metab& ac,const double& v) {
		int i,ipr; double dx, sum=0.;
		for(i=0;i<this->len;i++) {
			dx = v*(this->iso[i]);
			this->diso[i] -= dx;
		ipr=(i&15); pr.diso[ipr] += dx;
		ipr=(i>>4); ac.diso[ipr] += dx; sum += dx;		}
return sum;}

double getcon() const {return tcon;}
        
void split(Metab& pr1, Metab& pr2,const double v){
  int l1=pr1.getlen()-1, n1=pr1.getN();
   for(int i=0;i<this->len;i++) {
    int ip1=(i&l1), ip2=(i>>n1);
    double x=this->iso[i]*v;
    this->diso[i] -= x;
    pr1.diso[ip1] += x;// pr2.diso[ip2] += x;
       }
}
 
void splinverse(Metab& pr1, Metab& pr2,const double vf,const double vr){
  int l1=pr1.getlen()-1, n1=pr1.getN(), n2=pr2.getN();
   for(int i=0;i<this->len;i++) {
    int ip1(i&l1), ip(i>>n1),ip2(0);
    for(int i2=0;i2<n2;i2++) if(ip&(1<<i2)) ip2=(ip2|(1<<(n2-i2-1)));
    double x=this->iso[i]*vf-pr1.iso[ip1]*pr2.iso[ip2]*vr;
    this->diso[i] -= x;
    pr1.diso[ip1] += x; pr2.diso[ip2] += x;
        }
}
 
void condence(Metab& s1, Metab& s2,const double v){
  int n1=s2.getN(),ip,ip1; double x, dx;
   for(int i=0;i<s1.getlen();i++) { x=s1.iso[i]*v; ip=(i<<n1); //std::cout<<"i="<<i<<" ip="<<ip;
        for(int j=0;j<s2.getlen();j++){dx=x*s2.iso[j]; ip1=(ip|j); //std::cout<<"j="<<j<<" ip1="<<ip1<<std::endl;
//     std::cout<<"ip="<<ip1<<" i="<<i<<" j="<<j<<std::endl;
    this->diso[ip1] += dx;
    s2.diso[j] -= dx; s1.diso[i] -= dx; }}
}

	Metab(int n,std::string simya=""): N(n), len(1<<n), descr(simya) { } //
	virtual ~Metab(){ }
};

class ketose:public Metab {
 public:
        double ta123[8],tk12[4];
         void sett(){
		int i;
		for (i=0;i<8;i++) ta123[i]=0.;
		for (i=0;i<4;i++) tk12[i]=0.;
		int lena=(N-3);
		for (i=0;i<this->len;i++) ta123[i>>lena] += this->iso[i];
		for (i=0;i<8;i++) tk12[i>>1] += ta123[i];
		this->tcon=0; for (i=0;i<4;i++) this->tcon += tk12[i];
//		this->tcon= this->sumt();
		for (i=0;i<8;i++) ta123[i] /= this->tcon;
		for (i=0;i<4;i++) tk12[i] /= this->tcon;
	}
 void invistk(Metab& ald,const double v) {   
  int i,ii,j, ik,ib,lena=this->len>>2;
   double x,vk, va,sumx=0;
    const double *pak;
       pak=&this->tk12[0]; ib=4;
	vk=v/this->getcon();
   for(j=0;j<ib;j++){
     for(i=0;i<lena;i++) {ii=(i + j*lena);
 	x=this->iso[ii]*vk;
          sumx +=x;
         this->diso[ii] -= x;
          ald.diso[i] += x; }	
    }
     va= sumx/ald.getcon();
     for(j=0;j<lena;j++) {
          x= va*(ald.iso[j]);
            ald.diso[j] -= x;
        for (i=0;i<ib;i++) {
             ik=(j + i*lena);
                this->diso[ik] += x*pak[i];}
       }
 }
 void invista(Metab& ald,const double v) {   
  int i,ii,j, ik,ib,lena=this->len>>3;
   double x,vk, va,sumx=0;
    const double *pak;
       pak=&this->ta123[0]; ib=8;
	vk=v/this->getcon();
   for(j=0;j<ib;j++){
     for(i=0;i<lena;i++) {ii=(i + j*lena);
 	x=this->iso[ii]*vk;
          sumx +=x;
         this->diso[ii] -= x;
          ald.diso[i] += x; }	
    }
     va= sumx/ald.getcon();
     for(j=0;j<lena;j++) {
          x= va*(ald.iso[j]);
            ald.diso[j] -= x;
        for (i=0;i<ib;i++) {
             ik=(j + i*lena);
                this->diso[ik] += x*pak[i];}
       }
 }
 void tkk(Metab& ald1, Metab& ald2, ketose& ket2,const double v) {
     int i,ii,j,ik,ib,lena1=this->len>>2, lena2=ald2.getlen();
     double x,vk, va,sumx=0.;
     const double *pak;
   pak=&this->tk12[0]; ib=4; 
     vk=v/this->getcon();//kss; //vk /= (1+vk);
     ik=this->len; j=lena1-1;
         for (i=0;i<ik;i++) {ii=(i&j);
             x=this->iso[i]*vk;  
             sumx +=x;
             this->diso[i] -= x;   
             ald1.diso[ii] += x;
         }
    va= sumx/ald2.getcon();//suma;
    for(j=0;j<lena2;j++) {x= va*(ald2.iso[j]); ald2.diso[j] -= x;
        for (i=0;i<ib;i++) {ik=(j + i*lena2);
            ket2.diso[ik] += x*pak[i];    
            }
    }
 }

 void tka(Metab& ald1, Metab& ald2, ketose& ket2,const double v) {
     int i,ii,j,ik,ib,lena1=this->len>>3, lena2=ald2.getlen();
     double x,vk, va,sumx=0.;
     const double *pak;
        pak=&this->ta123[0]; ib=8; 
     vk=v/this->getcon();
     ik=this->len; j=lena1-1;
         for (i=0;i<ik;i++) {ii=(i&j);
             x=this->iso[i]*vk;  
             sumx +=x;
             this->diso[i] -= x;   
             ald1.diso[ii] += x;
         }
    va= sumx/ald2.getcon();
    for(j=0;j<lena2;j++) {x= va*(ald2.iso[j]); ald2.diso[j] -= x;
        for (i=0;i<ib;i++) {ik=(j + i*lena2);
            ket2.diso[ik] += x*pak[i];    
            }
    }
 }
        ketose(int n,std::string simya=""):Metab(n,simya){}
        ~ketose(){}
        };
        
class Ldistr {
int ntime,lmet, lmetk, itrac,markis,Nn;
static Metab pyrc, pyr, coa, oac, cit, akg, akgc, suc, mal, lacc, gl, lac, gln;
  double lacout,coaefl,tca,fpdh,marfrac,dsum;
  std::vector<Iso> result;
      void symm (double *s);
 void csyn(double *coai,double *dcoai,double *oai,double *doai,double *dciti,const double v);
 void split(double *h6,double *dh6,double *t1,double *dt1,const double vf,const double vr);
 void spInvsl(double *h6,double *dh6,double *t1,double *dt1,const double vf,const double tsum);
 void ast (double vf,double vr);


 void fitmet(double& xic,int iout,int iin,double fac);
 double fitm(double& xic,int ipar,double fac);
 double ser_gly (double v, double xthf);
 double gly_ser (double v);
  std::set<std::string> findopt(std::string a, std::vector<std::string> strok);
 void splitstrings(std::vector<std::string> segline[],int nstrok,std::vector<std::string> substrok);
 public:
 void wrikin(std::ostringstream& so, int nt){
   for(unsigned int i=0;i<expm0.size();i++) expm0[i]->wrikin(so,nt);
 }
 void shexper( int nt,int nexp=0){
   for(unsigned int i=0;i<expm0.size();i++) expm0[i]->exper[nexp].shexper(nt);
 }
 
 
 void wricon(std::ostringstream& so, int nt,int nexp=0){ so<<"\n";
   for(unsigned int i=0;i<expcon.size();i++) expcon[i]->wrikinc(so,nt);
   so<<"\n";
 }

  std::vector<double> tex, texcon;
  std::vector<int> vnn;
	void shocon();
       void setdiso(double *pyinit);
       void setiso(double *pyinit);
   void sklad(int itime);
   int wrim0ex(std::string );
   int wriconex(std::string fn);
   static Metab *met[];
   static ketose *metk[];
      std::vector<Metab*> expcon, expm0;
      std::vector<ketose*> kexpm0;
       void setmet();
       void setcon();
    void setfige();
    std::string read_con(std::ifstream& fi);
//       double setLacInit(double fact){return lac.setInCon(fact);}
       void read (char *finame, int ndat);
       void flback();
       int getN();
       int getNn(){return (Nn);}
       int getntime(){return (ntime);}
       void ssc(double *);
       void distr(double *py,double *pdydt);
       double consum();
      void readExp(char *fn,int ntr=0);
      int sexm0(std::string nm);
      void rmid(std::ifstream& ifi);
      std::vector<std::string> spli(std::stringstream& test,char a);
      int c13pos(std::string& s,int& nc,int& nlab);
      void defcol(int nucol[],std::vector<std::string> vstr);
      int findmet(Iso&);
      Tracer rcsv(std::ifstream& fi,std::vector<Iso>& result,int mar=0 );
        void show(std::ostringstream& fo,double xfin);
        void showcon(std::ostringstream& fo,double xfin);
        void showdescr(std::ostringstream& fo,std::vector<Metab*>);
        double xits(int its);
        double xicon(int its);
        void massfr();
        double integrbs();
	double ddisolve();
        int stor(double dist[]);
        void setflcon();
        int getmicon();
	Ldistr() {setmet();  }
	~Ldistr(void) { result.clear();}
};

extern Ldistr horse;

//---------------------------------------------------------------------------
#endif
 
