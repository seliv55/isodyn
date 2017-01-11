//---------------------------------------------------------------------------

#ifndef TKH
#define TKH
//#include "par.h"
class Doub3 {
public:
       double a[3];
       Doub3(){a[0]=0.;a[1]=0.;a[2]=0.;}
       ~Doub3(){}
};

class Doub5 {
public:
       double a[5];
       Doub5(void){a[0]=0.;a[1]=0.;a[2]=0.;a[3]=0.;a[4]=0.;}
       ~Doub5(){}
};
class TK{
       Doub5 kDe[8];
       double k[11],k0[11],v[11],v0[11];
       double e0;
public:
       void setk(double* );
       void Graf(double* const Dsp, double& g3i,double& r5i,double& f6i,double& x5i,double& ce4,double& cs7) const;
       void st1fl(double fl[], double& x5i,double& cs7,double& f6i,double& g3i,double& r5i,double& ce4) ;
       void st2fl(double fl[], double& x5i,double& cs7,double& f6i,double& g3i,double& r5i,double& ce4);
       TK(){}
       ~TK(){}
};
class TA{
       Doub3 kDe[7];
       double k[9],k0[9],v[9],v0[9];
       double e0;
public:
       void setk( double* );
       void Graf(double* const Dsp, double& g3i,double& f6i,double& ce4,double& cs7) const;
       void st1fl(double fl[], double& f6i,double& ce4, double& g3i,double& cs7) ;
       void st2fl(double fl[], double& f6i,double& ce4, double& g3i,double& cs7);
       TA(){}
       ~TA(){}
};
class Aldolase {
	double De,Defbp,Dedhap,Denom;
	double de[2],defbp[3],dedh[3];
       double k[4],k0[4],v[4],v0[4];
       double e0;
public:
       void setk( double* );
       void st1fl(double* const fl,double& g3i, double& dhap) ;
       void st2fl(double* const fl,double& g3i, double& dhap);
       Aldolase(){}
       ~Aldolase(){}
};
//---------------------------------------------------------------------------
#endif
