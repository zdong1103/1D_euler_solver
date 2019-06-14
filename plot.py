# Gnuplot file

set terminal qt size 1000,300 enhanced persist font "Times New Roman,16"
set multiplot layout 1,3
set linestyle 1 dt 1 lc -1 lw 1 # solid black line
set linestyle 2 dt 2 lc -1 lw 1 # - - black line
set linestyle 3 dt 4 lc -1 lw 1 # -.-. black line
set linestyle 4 dt 1 lc 1 lw 1 # solid red line

set ylabel offset 0.5
set xlabel offset 0
set grid
set format y "%.4g"

set ylabel "rho"
set yrange [:1.02]
set xlabel "x"
plot "./output/out.00048.dat" u 2:3 w lines ls 1 notitle

set ylabel "U"
unset yrange
plot "./output/out.00048.dat" u 2:4 w lines ls 1 notitle

set ylabel "p"
set yrange [:1.02]
plot "./output/out.00048.dat" u 2:5 w lines ls 1 notitle
