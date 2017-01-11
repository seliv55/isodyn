set t pngcairo enhanced dashed font "arial,12" size 1000,500; set output "figs/con.png";
#set t svg enhanced dashed font "arial,12" size 1000,500; set output "figs/con.svg";
set multiplot layout 2,5
set border 3;
set xtics 10 nomirror;
set ytics nomirror;
fkey=0.985;
set autoscale;

# * Numbers
t=1;  glc=t+1; glcsd=glc+1;  lacc=glcsd+1; laccsd=lacc+1;
lac=t+1; lacsd=lac+1; glu=lacsd+1; glusd=glu+1; ala=glusd+1; alasd=ala+1; gly=alasd+1; glysd=gly+1; ser=glysd+1; sersd=ser+1; pro=sersd+1; prosd=pro+1; cit=prosd+1; citsd=cit+1; agl=citsd+1; aglsd=agl+1; asp=aglsd+1; aspsd=asp+1; fum=aspsd+1; fumsd=fum+1; mal=fumsd+1; malsd=mal+1; coa=malsd+1; coasd=coa+1;

cglct=t+1; clact=cglct+1; lact=clact+1; pyrt=lact+1; pyrmt=pyrt+1; pyrsum=pyrmt+1; glut=pyrsum+1; alat=glut+1; glyt=alat+1; sert=glyt+1; prot=sert+1; citmt=prot+1; citct=citmt+1; citsum=citct+1; aglt=citsum+1; aspt=aglt+1; fumt=aspt+1; malt=fumt+1; oact=malt+1; malsum=oact+1;


set style line 1 lt 1 lw 1.5 lc rgb "#204a87"
set style line 2 lt 2 lw 1.5 lc rgb "#204a87"
set style line 3 lt 3 lw 1.5 lc rgb "#204a87"
set style line 4 lt 1 lw 1.5 lc rgb "#f57900"

# * Glucose
set label 1 "A"  at graph 0.15,0.975;
set ylabel "Concentration (mM)" offset 2,-1;
set ytics 2;
set key at graph 0.35,0.96 left Left  spacing 1.0 samplen 1 font "arial,11";
set xrange [0 to 25];
plot  [*:*] [0:*] 'excon' u ($1)/60:glc:glcsd w e pt 6 lc rgb "#204a87" not,'kinGlc' u ($1)/60:cglct w l ls 1 t "Glc",\
 'excon' u ($1)/60:lacc:laccsd w e pt 6 lc rgb "#f57900" not,'kinGlc' u ($1)/60:clact w l ls 4 t "Lac";

# * malate, akg
set ylabel "m0: malate" offset 2,0;
set key at graph 0.3,0.85
set label 1 "B";
set ytics 0.1;
plot 'exm0' u ($1)/60:mal:malsd w e pt 6 lc rgb "#204a87" not,\
 'kinGlc' u ($1)/60:malt w l ls 3 t "malm", 'kinGlc' u ($1)/60:oact w l ls 2 t "malc", 'kinGlc' u ($1)/60:malsum w l ls 1 t "sum";

# * citrate
set label 1 "C";
set ylabel "Citrate" offset 2,0;
plot 'exm0' u ($1)/60:cit:citsd w e pt 6 lc rgb "#204a87" not, 'kinGlc' u ($1)/60:citmt w l ls 3 t "citm", 'kinGlc' u ($1)/60:citct w l ls 2 t "citc", 'kinGlc' u ($1)/60:citsum w l ls 1 t "sum";

# * Lactate
set label 1 "D";
set ylabel "Lactate" offset 2,0;
plot 'exm0' u ($1)/60:lac:lacsd w e pt 6 lc rgb "#204a87" not,\
 'kinGlc' u ($1)/60:lact w l ls 3 t "lac", 'kinGlc' u ($1)/60:pyrt w l ls 1 t "pyr_g", 'kinGlc' u ($1)/60:pyrmt w l ls 2 t "pyr", 'kinGlc' u ($1)/60:pyrsum w l ls 4 t "sum",;

# * Glutamate
set label 1 "E";
set ylabel "Glutamate" offset 2,0;
plot 'exm0' u ($1)/60:glu:glusd w e pt 6 lc rgb "#204a87" t "glu"\
, 'kinGlc' u ($1)/60:glut w l ls 4 t "m0";

set xlabel "Time (h)" offset 0,0.5;

# * Alanine
set label 1 "F";
set ylabel "Alanine" offset 2,0;
plot 'exm0' u ($1)/60:ala:alasd w e pt 6 lc rgb "#204a87" t "ala"\
, 'kinGlc' u ($1)/60:alat w l ls 4 t "m0";

# * Serine
set label 1 "G";
set ylabel "Serine" offset 2,0;
plot 'exm0' u ($1)/60:ser:sersd w e pt 6 lc rgb "#204a87" t "ser"\
, 'kinGlc' u ($1)/60:sert w l ls 4 t "m0";

# * Proline
set label 1 "H";
set ylabel "Proline" offset 2,0;
plot 'exm0' u ($1)/60:pro:prosd w e pt 6 lc rgb "#204a87" t "pro"\
, 'kinGlc' u ($1)/60:prot w l ls 4 t "m0";

# * Glycerol
set label 1 "I";
set ylabel "Glycerol" offset 2,0;
plot 'exm0' u ($1)/60:agl:aglsd w e pt 6 lc rgb "#204a87" t "agl"\
, 'kinGlc' u ($1)/60:aglt w l ls 4 t "m0";

# * Fumarate
set label 1 "J";
set ylabel "Fumarate" offset 2,0;
plot 'exm0' u ($1)/60:fum:fumsd w e pt 6 lc rgb "#204a87" t "fum"\
, 'kinGlc' u ($1)/60:fumt w l ls 4 t "m0", ;

unset multiplot;

