#!/usr/bin/gnuplot --persist

set terminal tikz standalone color solid size 10cm,6cm
set output "res0.tex"
set xrange [-90:90]

unset key
set label "1-1"   at 5, 0.0075 
set label "2-2"   at 5, 0.003
set label "total" at 5, 0.01

set xlabel "$\\theta_0\\, [\^\\circ]$"
set ylabel "$I(\\theta_s|\\theta_0)\\, [rad\^{-1}]$"
plot "res0.dat" using 1:2 with lines title       "1-1", \
	 "res0.dat" using 1:3 with lines title       "2-2", \
	 "res0.dat" using 1:($2+$3) with lines title "total"

