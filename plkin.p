set t pngcairo enhanced dashed font "arial,12" size 1000,500;
set output filename;
set multiplot layout 2,5
set border 3;
set xrange [0:500]
set xtics 200 nomirror;
set ytics nomirror;
fkey=0.985;

# * Numbers
t=1;  nh6=t+1; nfbp=nh6+1; nt3=nfbp+1; npep=nt3+1; npyr=npep+1; npyrm=npyr+1; ncoa=npyrm+1; noa=ncoa+1; noac=noa+1; ncit=noac+1; ncitc=ncit+1; nakg=ncitc+1; nakgc=nakg+1; nfum=nakgc+1; nmal=nfum+1; np5=nmal+1; ne4=np5+1; ns7=ne4+1; natp=ns7+1; nnad=natp+1;

set style line 1 lt 1 lw 1.5 lc rgb "#f57900"
set style line 2 lt 1 lw 1.5 lc rgb "#4e9a06"

# * Glucose
set label 1 "A"  at graph 0.15,0.95;
set yrange [*:*]
set key at graph 0.6,0.3 left Right  spacing 1.0 samplen 1 font "arial,11";
plot fn2 u ($1)/60:nfbp w l ls 1 t "fbp", fn2 u ($1)/60:npyr w l ls 2 t "pyr"

set label 1 "B";
plot fn2 u ($1)/60:natp w l ls 1 t "atp", fn2 u ($1)/60:nnad w l ls 2 t "nad"

set label 1 "C";
plot fn2 u ($1)/60:ncoa w l ls 1 t "accoa"

set label 1 "D";
plot fn2 u ($1)/60:ncit w l ls 1 t "cit"

set label 1 "E";
plot fn2 u ($1)/60:nakg w l ls 1 t "akg"

set label 1 "F";
plot fn2 u ($1)/60:nmal w l ls 1 t "mal"

set label 1 "G";
plot fn2 u ($1)/60:npep w l ls 1 t "pep"

set label 1 "H";
plot fn2 u ($1)/60:ne4 w l ls 1 t "e4"

set label 1 "I";
plot fn2 u ($1)/60:noa w l ls 1 t "oa"

set label 1 "J";
plot fn2 u ($1)/60:npyrm w l ls 1 t "pyrm"

unset multiplot;

