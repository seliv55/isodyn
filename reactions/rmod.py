def mm( x):
   subend=x.index('->'); prend=x.index('$')
   sfl="flx["+x[0]+"]"; foo =sfl+"= "
   sub1="";ssub=""; sdx=""; parsets=''
   if int(x[prend+1]):
    foo += "rea[D].v()*"
   parsets= x[0]+' '+str(subend)+' v0 '
   if subend>1:
    if x[1][0]=='n':
       ssub += 'y[' +x[1]+']'
       parsets += ' k1 '
#    else: 
#       ssub += x[1]
#       parsets += ' k1 '
    if subend>2:
       ssub += ', y[' +x[2]+']'
       parsets += ' k2 '
   parsets+='\n'
   foo +="rea["+x[0]+"].v("+ssub+")"
   ssub=""
   for s in x[1:subend]:
     if s[0]=='n': ssub +="dydx["+s+"] -= "+sfl+";  "
   for s in x[(subend+1):prend]:
     if s[0]=='n': ssub +="dydx["+s+"] += "+sfl+";  "
   foo +='; \t'+ssub+'\n';
   return foo, parsets

def rinput(x):
   sf=''
   ssub=x[1][1:]; spr=x[3][1:]; #lrev=x[prend+5]
   if ssub[0]=='/': sf=ssub+".sumt()"; ssub=ssub[1:]
   rdist=ssub+".input("+ spr+",fluxes["+x[0]+"]"+sf
#   if len(flrev)-1: rdist += ",fluxes["+flrev+"]"
   rdist += ");\n"
   return rdist

def r3met(x,a):
   if a[1].find('+')>0:
     rtype='split'
     rdist=x[1][1:]+"."+rtype+"(" + x[3][1:]+","+x[4][1:]+",fluxes["+x[0]+"]"
   if a[0].find('+')>0:
     rtype='condence'
     rdist=x[4][1:]+"."+rtype+"(" + x[2][1:]+","+x[1][1:]+",fluxes["+x[0]+"]"
#   if len(flrev)-1: rdist +=",fluxes["+flrev+"]"
   rdist += ");\n"
   return rdist

def rm0in(x):
   prend=x.index('$')
   smet=x[prend+3]
   rdist=smet+".diso[0]+=fluxes["+x[0]+"];\n"
   return rdist

def rout(x):
   smet=x[1][1:]
   rdist=smet+".output(fluxes["+x[0]+"]);\n"
   return rdist

def irrev(x,a):
   if a[0].find(a[1])>0:
     rtype='cutfirst'
   elif a[0].find(a[1])==0:
     rtype='decarb'
   elif a[1].find(a[0])==0:
     rtype='carb';
   elif a[0].find(a[1][2:])==3:
     rtype='icdh';
   elif a[0].find(a[1][3:])==2:
     rtype='icdhr';
   smet=x[1][1:]; pmet=x[3][1:]
   rdist=smet+"."+rtype+"("+pmet+",fluxes["+x[0]+"]);\n"
   return rdist

def dextern(x):
  a="\tdydx["+x[0]+"] = ("
  for aa in x[1:]: a += aa[0]+'flx['+aa[1:]+']'
  a += ")*flx[rdt];\n"
  return a

def listname(sfirst, sname, sfin=''):
  aaa=sfirst
  for x in sname:  aaa+=x+", "
  if sfin=='': aaa=aaa[0:(len(aaa)-2)]
  aaa+=sfin+";\n\n";
  return aaa

def setvals(snames, sval):
  aaa=""
  for  i in range(len(snames)):
    aaa+=str(i)+") "+snames[i]+"\t"+sval+"\n"
  return aaa

def listint( sname, sfin, sin="0"):
  aaa="const int "+sname[0]+"="+sin+", "
  for i in range(1,len(sname)):
     aaa+=sname[i]+"="+sname[i-1]+"+1, " 
  aaa+=sfin+"="+sname[len(sname)-1]+"+1;\n\n"
  return aaa

lines_list = open('model').read().splitlines()
nvar=lines_list.index('fin')
vnut=lines_list.index('intern_var')
inmet=lines_list[1:vnut]
splinmet=[x.split() for x in inmet]
mnami=[x[2] for x in splinmet]
oumet=lines_list[(vnut+1):nvar]
sploumet=[x.split() for x in oumet]
mnamo=[x[2] for x in sploumet]
mname=mnami+mnamo
snx="const int numx="+mnamo[0]+";\n\n"
smeta=""; scona=""; sdata=""; lmet=0;
for x in splinmet:
   smeta += 'met['+str(lmet)+']=&'+x[2][1:]+'; '
   scona += "met["+str(lmet)+"]->setconc(xx["+x[2]+"]); ";
   sdata += " Ldistr::"+x[2][1:]+"("+x[1]+","+x[3]+"),";
   lmet+=1

for x in sploumet:
   smeta += 'met['+str(lmet)+']=&'+x[2][1:]+'; '
   scona += "met["+str(lmet)+"]->setconc(xx["+x[2]+"]); ";
   sdata += " Ldistr::"+x[2][1:]+"("+x[1]+","+x[3]+"),";
   lmet+=1

smeta += "\n lmet="+str(lmet)+"; "
f_ff="void Fit::f(const double *y,double *dydx) {\n\tfor(int i=0;i<nmet;i++) dydx[i]=0.;\n\tfor(int i=0;i<nflx;i++) flx[i]=0.;\n";
flfor="void Parray::flfor(double *y){\nfor(int i=0;i<nflx;i++) fluxes[i] = flx[i] * flx[rdt]/Vi;\n"
rstar=lines_list.index('reactions:')
rend=lines_list.index('rend')
react=lines_list[(rstar+1):rend]
splir=[x.split() for x in react]
flname=[x[0] for x in splir]
a1=''; a2=''; pars=''; s_fl=''; i=0
stdist=""
for x in splir:
   if x[1][0] != '/':
     s_fl ="fluxes["+x[0]+"] /= "
     if (x[1][0]=='n') & (x.index('->') == 2):
       s_fl += 'y[' +x[1]+']'
       flfor +=s_fl+";\n"
     if (x[1][0]=='y') & (x.index('->') == 2):
       s_fl += x[1]
       flfor +=s_fl+";\n"
     if x.index('->') == 3:
       s_fl += 'y['+x[1]+']*y['+x[2]+']'
       flfor +=s_fl+";\n"
   a1,a2 = mm(x)
   f_ff += a1
   pars += str(i)+') '+ a2
   i+=1
   prend=x.index('$')
   atr=x[prend+2].split('>')
   if len(atr)==2:
     if atr[0]==atr[1]:
       aaa=rinput( x)
       stdist += aaa
     elif atr[1]=='':
       aaa=rout(x)
       stdist+=aaa
     elif prend<5:
       aaa=irrev(x,atr)
       stdist += aaa
     elif prend==5:
       aaa=r3met(x,atr)
       stdist+=aaa

flfor += '}\n'
kfl=len(react)
a=pars.split('\n'); pars=''
a.remove('')
fipar=open('../output/1').read().splitlines()
rind=[fipar.index(y) for y in fipar if 'rdt' in y][0]
reapar=fipar[0:(rind+1)]
i=0
for x in splir:
  ind=[reapar.index(y) for y in reapar if x[0] in y]
  if len(ind)>0:
     a[i]=fipar[ind[0]]
  i+=1

for x in a:
  pars += x+'\n'

f_ff+="for(int i=0;i<nmet;i++) dydx[i]*=(flx[rdt]/Vi);\n}\n\n"
f_ff+="void Fit::ff(const double *y,double *dydx) {\n\tf(y,dydx);\n";
fin=lines_list[rend:].index('fin')
splir=[x.split() for x in lines_list[(rend+2):(rend+fin)]]
for x in splir: f_ff += dextern(x)

f_ff+="}\n\n"
hlfl=listname('extern const int ',flname,"nrea, nflx")
hlfl+=listname('extern const int ',mname,"numx, nmet");
hlfl+="\nextern const double thft;\n";
hlfl+="extern double mader, dif,dif0, suxx, Vi,Vt, mu, xx[], fluxes[], flx[];\n";
hlfl+="extern double nv1[],nv2[],xinit1[],xinit2[];\nextern time_t ts,tf,tcal;\n";
hlfl+="extern std::string foc, kin, kinc;\nextern char *fex1, *fex2;\nextern int ifn;\n";
numshh=open("../include/nums.hh",'w')
numshh.write(hlfl); numshh.close(); hlfl=""

metind=[fipar.index(y) for y in fipar if 't=' in y][0]-1
metpar=fipar[(rind+2):metind]
meti=setvals(mname, "0.01").split('\n')
meti.remove('')
i=0
for x in mname:
  ind=[metpar.index(y) for y in metpar if x in y]
  if len(ind)>0:
     meti[i]=metpar[ind[0]]
  i+=1

ini=''
for x in meti:
  ini += x+'\n'

pars += fipar[rind+1]+"\n"+ini+setvals(flname, "0.001")
param=open("parameters",'w')
param.write(pars); param.close(); pars=""

parsets="Metab *Ldistr::met["+str(lmet)+"];\n";

nvcpp="#include <iostream>\n#include \"nr.h\"\n#include \"nums.hh\"\n#include \"tk.hh\"\n#include \"nv.hh\"\n#include \"modlab.h\"\n#include \"solvers.h\"\n#include \"analis.h\"\n";
nvcpp += "using namespace std;\n"
#definition of the list of fluxes, parameters and variables for "nv.cpp"
nvcpp +=listint(flname,"nrea")+"const int nflx=nrea;\n" + listint(mname,"nmet")+snx
aa="Metab"+sdata[0:(len(sdata)-1)]
nvcpp +=aa+";\n";
nvcpp+=("\tFit Problem;\n\tconst double thft(1.);\n\tdouble xx[nmet],flx[nflx],fluxes[nflx];\n\tdouble xinit1[nmet],xinit2[nmet];\n\t");
nvcpp +="string Parray::fid[nflx],Parray::fname[nflx],Parray::fschem[nflx], Parray::namex[nmet];\n\tReapar Parray::rea[nrea];\n\tdouble Analis::nv1[nrea], Analis::nv2[nrea];\n"+parsets+"\n void Ldistr::setmet(){";
parsets=""
nvcpp +=smeta + " }\n void Ldistr::setcon(){"+ scona+ " }\n"+ f_ff;
nvcpp +="void Parray::init(){ft3=10.; fh6=7.;}\n";
nvcpp +="void Parray::fin(double y[]){\n\t"
nvcpp +="\nflfor(y);\n}\n"+flfor;
fout=open("../con512tpl/nv.cpp",'w')
fout.write(nvcpp)
fout.close(); nvcpp="";

sdist="#include <iostream>\n#include \"nums.hh\"\n#include \"tk.hh\"\n#include \"nv.hh\"\n#include \"modlab.h\"\n"
sdist += "using namespace std;\nvoid Ldistr::distr(double *py,double *pdydt) {\n\tdouble NOL=0.;\n\tsetiso(py); setdiso(pdydt);\n\t"
sdist += "for (int i=0;i<Nn;i++) pdydt[i]=0.;\n\t"
sdist += "Problem.ff(py,pdydt);\n\tProblem.fin(py);/**/\n";
sdist +=stdist;
for x in sploumet: sdist += x[2][1:]+".volume(Vt);\n\t";

sdist +="symm(mal.getisot());\n";
for x in mname: sdist += "py["+x+"]="+x[1:]+".sumt(); "

sdist += "\n//";
for x in mname: sdist += "xx["+x+"]=py["+x+"]; ";

sdist += "\n}\n";
fout=open("../con512tpl/distr.cpp",'w'); fout.write(sdist);  fout.close()

mn=[x[1:] for x in mname]
stri= listname("static Metab ", mn)
stri= stri[0:(len(stri)-2)]
fout=open("output",'w'); fout.write(stri);  fout.close()

lines_list = open('../include/modlab.h').read().splitlines()
nvar=lines_list.index('class Ldistr {')
lines_list[nvar+2]=stri
aa=''
for x in lines_list:
  aa+=x+'\n' 

fout=open("../include/modlab.h",'w'); fout.write(aa);  fout.close()

