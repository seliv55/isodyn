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
          for(int i=0;i<len;i++) mid[ic].mean+=linj[i].mid[ic].mean; mid[ic].mean/=(double)len; }
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
  virtual ~Iso(){}
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

class Metab{
 protected:
      const int N, len;
      bool flag;
      double fc,tcon, *calc, kinc[tt],kinm0[tt];
      data conc[tt];
  std::string descr;
 public:
	double *iso, *diso;
	int ny;
 double* getcalc(){return &calc[0];}
 int getN(){return N;}
 std::string getdescr(){return descr;}
	
 double percent(const int left=0,const int riht=0){
		double tot=0.; int nc,sum1;
		for(int i=0;i<N+2;i++) calc[i]=0.;
   for(nc=0;nc<len;nc++) {
    int i=((nc>>riht)&((len>>left>>riht)-1)); sum1=0;
	while(i>0){sum1 += (i & 1);i = (i >> 1);}
			calc[sum1] += iso[nc];
		}//end for
	for(int i=0;i<=N;i++) tot += calc[i];
//	for(int i=0;i<=N;i++) calc[i] /= tot;
		return calc[N+1]=tot;
	}
	
 void setconc(double a, int nt=0){ conc[nt].mean=a; conc[nt].sd=a*0.1; }
 
 void skladc(int itime){  kinc[itime]=calc[N+1]; }
 
 void skladm0(int itime){ kinm0[itime]=calc[0]/calc[N+1]; }
 
 double sumt(){
		double sum=0.; unsigned i;
		for (i=0;i<len;i++) sum += iso[i];
		return calc[N+1]=sum;
	}
 double input(Metab& s2,const double& vi,const double& vo=0.){
		double x, sum=0.;
		for(int i=0;i<len;i++) {
			x=(vi*(iso[i])-vo*(s2.iso[i]));
			diso[i] -= x; s2.diso[i] += x; sum += x;
		}
	return sum;}
 void output(const double& vi){
		double x, sum=0.;
		for(int i=0;i<len;i++) {
			x=vi*iso[i];	diso[i] -= x;}
	}
 
 void set0(const double& s0){
		iso[0]=s0;
		for(int i=1;i<len;i++) iso[i]=0.;
	}
 int getlen() const {return len;	}
 double* getisot() {return iso;	}

 void showsum(Metab& s2,std::ostringstream& fo){
   fo<<" "<<std::setw(7)<<(this->calc[0]+s2.calc[0])/(this->calc[N+1]+s2.calc[N+1]);
 }
 void showm0(std::ostringstream& fo){ double m0=this->calc[0]/this->calc[N+1];
   fo<<" "<<std::setw(7)<<m0;
 }

 void showcon(std::ostringstream& fo){
  fo<<" "<<std::setw(7)<<this->calc[N+1];
 }

 void coriso(const double xkin,const double xiso){
 double ff=xkin/xiso;
  for(int i=0;i<len;i++) iso[i] *= ff;
  }

	void wcalc(std::ofstream& fi) {
	unsigned k;
	percent();
		for (k=0;k<(N+2);k++) 
       if((N+1)-k) fi<<"m"<<k<<" "<<calc[k]/calc[N+1]<<std::endl;
         else fi<<"conc: "<<calc[k]/calc[N+1]<<std::endl; 
			}
 double cutfirst(Metab& p,const double& v) {
		double dx, sum=0.;
		int i=0,i1;
		int lenp=((len>>1)-1);
		for(int i=0;i<this->len;i++) {i1 = (i & lenp);
			dx = v*(this->iso[i]);
		this->diso[i] -= dx; p.diso[i1] += dx; sum += dx;
		}
	return sum;}
	

	void volume(double factor) {
	for (int i=0; i<len;i++) diso[i] *= factor;
	}
	void carb (Metab& pr,const double& v){
		int ipr; double dx;
		for(int i=0;i<this->len;i++) {
			dx = v*(this->iso[i]);
			this->diso[i] -= dx;
		ipr=(i<<1);		pr.diso[ipr] += dx;
		}
	}
	void decarb (Metab& pr,const double& v) {
		int i,ipr; double dx;
		for(i=0;i<this->len;i++) {
			dx = v*(this->iso[i]);
			this->diso[i] -= dx;
		ipr=(i>>1);		pr.diso[ipr] += dx;
		}
	}
	void icdh (Metab& pr,const double& v) {
		int i,ipr; double dx;
		for(i=0;i<this->len;i++) {
			dx = v*(this->iso[i]);
			this->diso[i] -= dx;
		ipr=((i&7)|((i>>1)&24)); pr.diso[ipr] += dx;
		}
	}
	void icdhr (Metab& pr,const double& v) {
		int i,ipr; double dx;
		for(i=0;i<this->len;i++) {
			dx = v*(this->iso[i]);
			this->diso[i] -= dx;
		ipr=((i&7)|((i&24)<<1)); pr.diso[ipr] += dx;
		}
	}
	double liase (Metab& pr,Metab& ac,const double& v) {
		int i,ipr; double dx, sum=0.;
		for(i=0;i<this->len;i++) {
			dx = v*(this->iso[i]);
			this->diso[i] -= dx;
		ipr=(i&15); pr.diso[ipr] += dx;
		ipr=(i>>4); ac.diso[ipr] += dx; sum += dx;
		}
	return sum;}
       void sett(){ tcon= this->sumt();}
       
        double getcon() const {return tcon;}
        
 void split(Metab& pr1, Metab& pr2,const double v){
  int l1=pr1.getlen()-1, n1=pr1.getN();
   for(int i=0;i<this->len;i++) {
    int ip1=(i&l1), ip2=(i>>n1);
    double x=this->iso[i]*v;
    this->diso[i] -= x;
    pr1.diso[ip1] += x; pr2.diso[ip2] += x;
    }
 }
 
 void splinverse(Metab& pr1, Metab& pr2,const double vf,const double vr){
  int l1=pr1.getlen()-1, n1=pr1.getN(), n2=pr2.getN();
   for(int i=0;i<this->len;i++) {
    int ip1(i&l1), ip(i>>n1),ip2(0);
    for(int i=0;i<n2;i++) if(ip&(1<<i)) ip2=(ip2|(1<<(n2-i-1)));
    double x=this->iso[i]*vf-pr1.iso[ip1]*pr2.iso[ip2]*vr;
    this->diso[i] -= x;
    pr1.diso[ip1] += x; pr2.diso[ip2] += x;
    }
 }
 
 void condence(Metab& s1, Metab& s2,const double v){
  int n1=s1.getN(),ip,ip1; double x, dx;
   for(int i=0;i<s2.getlen();i++) { x=s2.iso[i]*v; ip=(i<<n1); 
        for(int j=0;j<s1.getlen();j++){dx=x*s1.iso[j]; ip1=(ip|j);
//     std::cout<<"ip="<<ip1<<" i="<<i<<" j="<<j<<std::endl;
    this->diso[ip1] += dx;
    s1.diso[j] -= dx; s2.diso[i] -= dx; }}
 }
        
	Metab(int n,std::string simya=""): N(n), len(1<<n), descr(simya) { calc=new double[N+2];} //
	virtual ~Metab(){ delete[] calc;}
};

class Metab_data:public Metab {
  data  *exper[tt]; double *calstor[tt];
  double xicon[tt], xi[tt];
  int mi;
	public:
  bool flcon;
  void setex0(){
     exper[0][0].mean=1.; exper[0][0].sd=0.01;
      for(int i=1;i<=mi;i++){ exper[0][i].mean=0.; exper[0][i].sd=0.01;  }
                         }
   
 void setrav(int nt1,int nt2){
      for(int i=0;i<=N;i++){ exper[nt1][i].mean=exper[nt2][i].mean; exper[nt1][i].sd=exper[nt2][i].sd;  }
                         }
   
 void sex(int niso, data mid[],int nt){ mi=niso;//exper[nt]=mid;
   for(int i=0;i<niso;i++) exper[nt][i].mean=mid[i].mean;
   for(int i=0;i<niso;i++) {exper[nt][i].sd=mid[i].sd; if(exper[nt][i].sd<0.01) exper[nt][i].sd=0.01;}
//   for(int i=0;i<niso;i++) { if(mid[i].sd<0.01) mid[i].sd=0.01;}
	}
	
 void readc(std::ifstream& fi,  int nt){ std::string aaa;
     fi>>aaa;
     for(int i=0;i<nt;i++) fi>>conc[i].mean;
     for(int i=0;i<nt;i++) {fi>>conc[i].sd; if(conc[i].sd<0.01) conc[i].sd=0.01;}
	}

 data * getconc(){return conc;}
 
 int getmi(){ return mi;}
 
 int storc(double dis[],int nt) {if(conc[nt].mean>1e-5) dis[0]=this->calc[N+1]/conc[nt].mean; else dis[0]=0.; return 1;}
 
 int stormi(double dis[],int nt) {
     for(int j=1;j<nt;j++) for(int i=0;i<mi;i++) dis[i+(j-1)*mi]=this->calstor[j][i];
      return mi*(nt-1);}
 
 void rizeflcon(){flcon=1;}
  
 data * getexper(int nt){return &exper[nt][0];}

 void showmi(std::ostringstream& fo,int nt){
   for (int k=0;k<=N;k++) fo<<" "<<std::setw(7)<<(this->calc[k]/this->calc[N+1]); fo<<"\n";
   fo<<"Data";for (int k=0;k<=N;k++) fo<<" "<<std::setw(7)<<exper[nt][k].mean; fo<<"\n";
 }

 double chisq(int nt) {     xi[nt]=0.; double a; 
  for(int k=0;k<mi;k++) if(exper[nt][k].sd>5.e-3){
   calstor[nt][k]=calc[k]/calc[N+1];
    a=(exper[nt][k].mean - calstor[nt][k])/(exper[nt][k].sd); xi[nt] += a*a;}
       else calstor[nt][k]=0.;
         return xi[nt];}
  
 double chisqsum(Metab& s2,int nt,int nl) {
	double a, b=this->calc[N+1]+s2.getcalc()[N+1];
	 xi[nt]=0; 
  if(exper[nt][0].sd>1.e-5) for(int k=0;k<nl;k++) {
     a=(exper[nt][k].mean-(this->calc[k]+s2.getcalc()[k])/b)/(exper[nt][k].sd);
      xi[nt] += a*a; }
  return xi[nt];}

 double chicon(int nt) {
	double a;
     a=(conc[nt].mean-this->calc[N+1])/(conc[nt].sd);
  return (xicon[nt]=a*a);}

 void wrikinc(std::ostringstream& so, int nt) {
  so.precision(3);
    so<<"\n"<<std::setw(12)<<descr<<"_c: ";
    for(int i=0;i<nt;i++) so<<std::setw(9)<<this->kinc[i];
    so<<"** "<<std::setw(7)<<conc[nt-1].mean<<" -> "<<xicon[nt-1];
 }
 void wrikinm0(std::ostringstream& so, int nt) {
     if(xi[nt-1]>0){ so.precision(3);
      so<<"\n"<<std::setw(12)<<descr<<"_m0: ";
    for(int i=0;i<nt;i++)  so<<std::setw(9)<<this->kinm0[i];
    so<<" : "<<std::setw(7)<<exper[nt-1][0].mean<<" -> "<<xi[nt-1]; }
 }
void setexper(int nt){for(int i=0;i<nt;i++) {exper[i]=new data[N+1]; calstor[i]=new double[N+1];} exper[0][0].mean=1.; }
void delexper(int nt){for(int i=0;i<nt;i++) {delete[] exper[i]; delete[] calstor[i];}}
  Metab_data(int n,std::string simya=" "):Metab(n,simya){flcon=0;}
  ~Metab_data(){}
};

class ketose:public Metab_data {
 public:
        double ta123[8],tk12[4];
 double dilut(double& dlt,data exper[],const int riht=0,const int left=0) {
	     unsigned i; double xi(0.0),a(0.0);
	   this->percent( riht, left);
	dlt = (exper[0].mean-this->calc[0])/(1.-exper[0].mean);
		if(dlt>0) {
		this->calc[0] += dlt;
		for(i=0;i<=N;i++) this->calc[i] /= (1.+dlt);this->calc[N+1] *= (1.+dlt);
		}
		for(i=0;i<(N+2);i++) {    //          for(int i=1;i<repl;i++)
			a = (exper[i].mean-this->calc[i])/(exper[i].sd);
			xi += a*a;
		}
		return xi;
	}
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
        ketose(int n,std::string simya=""):Metab_data(n,simya){}
        ~ketose(){}
        };
        
 class twofrag:public Metab_data{
   double *calc1;
   data *exper1[tt];
  int mi1;
 public:
   
 void setex10(){
     exper1[0][0].mean=1.; exper1[0][0].sd=0.01;
      for(int i=1;i<=N;i++){ exper1[0][i].mean=0.; exper1[0][i].sd=0.01;  }
                         }
                         
 void read1(std::ifstream& fi,int nt){
	unsigned i; std::string aaa;//i: isotopomer; j: time
       fi >> mi1;
   for(i=0;i<=mi1;i++) { fi>>exper1[nt][i].mean;} fi>>aaa; fi>>aaa;
   for(i=0;i<=mi1;i++) {fi>>exper1[nt][i].sd; if(exper1[nt][i].sd<0.01) exper1[nt][i].sd=0.01;}
	}

        twofrag(int n,std::string simya=""):Metab_data(n,simya){}
        ~twofrag(){}
};

class Ldistr {
	int ntime,lmet,lmetb, lmetk, itrac,markis,Nn;
  Metab_data gl, lac, glu, gln, rna, glycog, pro, asp, ala, ser, agl, pyrm, coa, coac, gly, oa, oac, cit, akg, fum, mal, glu25, fbp, t3, pep, pyr, citc, akgc, e4;
  Metab cthf,dhe,gae;
  ketose h6, s7, p5;
  double lacout,coaefl,tca,fpdh,marfrac;
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
   for(int i=0;i<expm0.size();i++) expm0[i]->wrikinm0(so,nt);
 }
 
 void wricon(std::ostringstream& so, int nt){ so<<"\n";
   for(int i=0;i<expcon.size();i++) expcon[i]->wrikinc(so,nt); so<<"\n";
 }

  std::vector<double> tex, texcon;
       void setdiso(double *pyinit);
       void setiso(double *pyinit);
   void sklad(int itime);
   int wrim0ex(std::string );
   int wriconex(std::string fn);
      static Metab_data *met[];
      static Metab *metb[];
      static ketose *metk[];
      std::vector<Metab_data*> expcon, expm0;
      std::vector<ketose*> kexpm0;
       void setmet();
       void setcon();
    void setfige();
    std::string read_con(std::ifstream& fi, std::string& arg1,int& nfi);
//       double setLacInit(double fact){return lac.setInCon(fact);}
       void read (char *finame, int ndat);
       void flback();
       int getN();
       int getNn(){return (Nn);}
       int getntime(){return (ntime);}
       void ssc(double *);
       void distr(double *py,double *pdydt);
       double consum();
       void showmi(std::ostringstream& fo,int nt);
      double readExp(char *fn,int ntr=0);
      std::vector<std::string> spli(std::stringstream& test,char a);
      int c13pos(std::string& s,int& nc,int& nlab);
      void defcol(int nucol[],std::vector<std::string> vstr);
      int findmet(Iso&);
      Tracer rcsv(std::ifstream& fi,std::vector<Iso>& result,int mar=0 );
        void show(std::ostringstream& fo,double xfin);
        void showcon(std::ostringstream& fo,double xfin);
        void showdescr(std::ostringstream& fo,std::vector<Metab_data*>);
        double xits(int its);
        double xicon(int its);
        void massfr();
        double integrbs();
	double ddisolve();
        int stor(double dist[]);
        void setflcon();
        int getmicon();
	Ldistr(): gl(6,"Gluc"), lac(3,"Lac"), glu(5,"Glutamate2-5"), gln(5,"Glutamin"), rna(5,"Rib"), glycog(6,"Glycog"), pro(5,"Pro"), asp(4,"Asp"), ala(3,"Ala"), ser(3,"Ser"), agl(3,"Glycerol"), pyrm(3,"Pyr"), coa(2,"CoA"), coac(2,"coac"), gly(2,"Gly"),  oa(4,"Oaa"), oac(4,"oac"), cit(6,"Cit"), akg(5,"aKg"), fum(4,"Fum"), mal(4,"Mal"), glu25(5,"glutamate2-5"),
	 fbp(6), t3(3), pep(3), pyr(3), cthf(1),dhe(3),gae(2), citc(6), akgc(5), e4(4),
	  h6(6), s7(7), p5(5,"rib")  {setmet(); getN(); }
	~Ldistr(void) { result.clear(); for(int i=0;i<expm0.size();i++) expm0[i]->delexper(ntime);}
};

extern Ldistr horse;

//---------------------------------------------------------------------------
#endif
 
