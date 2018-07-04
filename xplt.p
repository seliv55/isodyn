set t pngcairo enhanced dashed font "arial,12" size 1000,300; set output "figs/m0.png";
#set t svg enhanced dashed font "arial,12" size 1000,300; set output "figs/con.svg";
set multiplot layout 1,m0
set border 3;
set xrange [0:1500]
set xtics 500 nomirror;
set ytics nomirror;
fkey=0.985;
set autoscale;
set key autotitle columnhead
set xlabel 'min' offset 0,0.5

set style line 1 lt 1 lw 1.5 lc rgb "#204a87"
set style line 2 lt 2 lw 1.5 lc rgb "#204a87"
set style line 3 lt 3 lw 1.5 lc rgb "#204a87"
set style line 4 lt 1 lw 1.5 lc rgb "#f57900"
set ylabel 'm0' offset 2,0
do for [n=1:m0]{
plot 'exm0' u 1:(column(n*2)) w p pt 3 t columnhead(n*2), 'kinGlc' u 1:(column(n+1)) w l
  }
unset multiplot;
set output "figs/con.png";
set multiplot layout 1,con
set ylabel 'concentration' offset 2,0
do for [n=1:con]{
plot 'excon' u 1:(column(n*2)) w p pt 3 t columnhead(n*2), 'kincon' u 1:(column(n+1)) w l
  }
unset multiplot;

