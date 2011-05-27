#!/usr/bin/gnuplot

set terminal tikz standalone color solid size 9cm,8cm
set output "polarplot.tex"

# Settingup stuff 
set polar
set angles degrees
set size ratio -1
set yrange [ 0 :.016]
set xrange [-.012 :.012]
set grid polar 30 lc rgb "gray" lt 1

# Removing unwanted stuff
unset key
unset xtics
unset ytics
unset mxtics
unset mytics
unset x2tics
unset y2tics

# Setting tick labels around the figure
set label 1 "-60$^\\circ$"  at -0.0104*cos(30), 0.0104*sin(30) rotate by  60 center 
set label 2 "60$^\\circ$"   at  0.0104*cos(30), 0.0104*sin(30) rotate by -60 center
set label 3 "-30$^\\circ$"  at -0.0154*cos(60), 0.0154*sin(60) rotate by  30 center
set label 4 "30$^\\circ$"   at  0.0154*cos(60), 0.0154*sin(60) rotate by -30 center
set label 5 "0$^\\circ$"    at  0.0154*cos(90), 0.0154*sin(90) center

set label 6 "0.005" at 0.0053*sin(15), 0.0053*cos(15) rotate by -15 center
set label 7 "0.010" at 0.0103*sin(15), 0.0103*cos(15) rotate by -15 center
set label 8 "0.015" at 0.0153*sin(15), 0.0153*cos(15) rotate by -15 center

# Creating linestyles
set style line 1 lt 1 lw 1 lc rgb "blue"
set style line 2 lt 2 lw 1 lc rgb "red"

# Setting arrow for incident light 
set arrow from .015*cos(90) ,.015*sin(90)  to .013*cos(90)  ,.013*sin(90)  filled ls 1
set arrow from .015*cos(130),.015*sin(130) to .013*cos(130) ,.013*sin(130) filled ls 2

# Plotting the totoal mean diffusive scattering intensity
plot 'res0.dat' using (-$1+90):($2+$3) with lines ls 1,\
	 'res40.dat' using (-$1+90):($2+$3) with lines ls 2
