#!/usr/bin/gnuplot --persist

clear
reset

set terminal tikz standalone color solid size 7.5cm, 5cm
set output "res40.tex"

# Some setup
unset key
set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0
set size 1,1
set origin 0,0
set multiplot

# Main plot
set size 1,1
set origin 0,0

set xrange [-90:90]
set yrange [0:0.015]
set xtics 30
set mxtics 3
set ytics 0.005
set mytics 5

set label 1 "1-1"   at 0, 0.0087 	right 
set label 2 "2-2"   at 0, 0.0025 	right 
set label 3 "-3-1"	at 0, 0.0008 	right 
set label 4 "total"	at 0, 0.0102 	right

set xlabel "$\\theta_\\mathrm{s}\\, [\^\\circ]$"
set ylabel "$I(\\theta_\\mathrm{s}|\\theta_0)\\, [\\mathrm{rad}\^{-1}]$"
plot "res40final.dat" using 1:2 			with lines title "1-1" 		lt rgb "red", \
	 "res40final.dat" using 1:3 			with lines title "2-2" 		lt rgb "green", \
	 "res40final.dat" using 1:4 			with lines title "2-2" 		lt rgb "orange", \
	 "res40final.dat" using 1:($2+$3-$4) 	with lines title "total" 	lt rgb "blue"

# Inset
unset xlabel
unset ylabel
unset label 1
unset label 2
unset label 3
unset label 4
set size 0.2,0.25
set origin 0.75,0.70

set xrange [-42:-38]
set yrange [.008:.0095]
set xtics 2
set mxtics 2
set ytics 0.001
set mytics 2

plot "res40zoom.dat" using 1:($2+$3-$4) with lines lt rgb "blue"

unset multiplot
