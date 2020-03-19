set t pngcairo enhanced dashed font "arial,12" size 1000,900; set output "figs/m0.png";
#set t svg enhanced dashed font "arial,12" size 1000,600; set output "figs/con.svg";
set multiplot layout 3,m0
set border 3;
#set xrange [0:1500]
set yrange [0:100]
set xtics 500 nomirror;
set ytics nomirror;
fkey=0.985;
#set autoscale;
set key autotitle columnhead
set xlabel 'time(min)' offset 0,0.5
set key right above
set style line 1 lt 1 lw 1.5 lc rgb "#204a87"
set style line 2 lt 2 lw 1.5 lc rgb "#204a87"
set style line 3 lt 3 lw 1.5 lc rgb "#204a87"
set style line 4 lt 1 lw 1.5 lc rgb "#f57900"
set ylabel '% unlabeled' offset 2,0
do for [n=1:(3*m0)]{
plot 'exm0' u 1:(100*column(n*2)):(100*column(n*2+1)) w err pt 2 t columnhead(n*2), 'kinGlc' u 1:(100*column(n+1)) w l ls 4
  }
unset multiplot;
set t pngcairo enhanced dashed font "arial,12" size 700,300;
set output "figs/con.png";
set multiplot layout 1,con
set autoscale;
set ylabel 'concentration(mM)' offset 2,0
do for [n=1:con]{
plot 'excon' u 1:(column(n*2)):(column(n*2+1)) w err pt 2 t columnhead(n*2), 'kincon' u 1:(column(n+1)) w l
  }
unset multiplot;
# w p pt 3 t columnhead(n*2), 'kinGlc' u 1:(100*column(n+1)) w l
