#ifndef ANALISH
#define ANALISH
extern double ystart[];
class Analis {
double chimin, xmin, tmin, dif1, xm, x00;
static double nv1[],nv2[];
double descent(double, int ip=(-1)); 
 void rconfint(ifstream& fi, double ami[],double ama[]);
 void stepdown(double oldp, int& i,double fdes);
double dermax();
double xmax();
double minx(){double mx=xx[0]; for(int i=1;i<numx;i++) if(xx[i]<mx)mx=xx[i]; return mx;}
void revsvd(Mat_I_DP &aa,Vec_I_DP &w, Mat_I_DP &ai,Mat_O_DP &u);
void matmult(Mat_I_DP &aa, Mat_I_DP &ai,Mat_O_DP &u);
void hessian(int ndim,int np,int nex,double xi0,double st[], Mat_DP& aa,Vec_DP& b,int par1[],double tmax);
void cross(double p1[],double p2[]);
void mutate(int npf,int par[]);
public:
double getmax();
void setx00(const double x){x00=x;}
void fitm0(const double tmax);
void sensitiv(const double tmax);
void coord(const double f1, double);
void confidence(double,double);
void grad(double tmax);
void swarm(const double tmax,const int ngen);
void genetic(const double tmax,const int ngen);
Analis(){}
~Analis(){}
};
extern double solve();
#endif

