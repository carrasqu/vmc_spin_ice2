set terminal postscript color enhanced "Helvetica,25"
set output 'cv.eps'

set multiplot
set origin 0.0,0.0
set size 0.56,1
unset key
set xlabel 'T'
set ylabel 'Energy'
set xtics 0.5

p 'L_1/data_1.dat' u 1:2 w lp lc 0 title 'L=1', 'L_1/data_1.dat' u 1:2:3 w e lc 0  notitle, 'L_2/data_2.dat' u 1:2 w lp lc 1 title 'L=2', 'L_2/data_2.dat' u 1:2:3 w e lc 1  notitle ,'L_3/data_3.dat' u 1:2 w lp lc 3 title 'L=3', 'L_3/data_3.dat' u 1:2:3 w e lc 3  notitle,'L_4/data_4.dat' u 1:2 w lp lc 4 title 'L=4', 'L_4/data_4.dat' u 1:2:3 w e lc 4  notitle, 'en.dat' u 1:($2+4) w lp lc rgb "#FF7F50" title 'ED L=1' 


set origin 0.5,0.0
set size 0.56,1
set key
set xlabel 'T'
set ylabel 'specific heat'
p 'L_1/data_1.dat' u 1:4 w lp lc 0 title 'L=1', 'L_1/data_1.dat' u 1:4:5 w e  lc 0 notitle, 'L_2/data_2.dat' u 1:4 w lp lc 1 title 'L=2', 'L_2/data_2.dat' u 1:4:5 w e  lc 1 notitle, 'L_3/data_3.dat' u 1:4 w lp lc 3 title 'L=3', 'L_3/data_3.dat' u 1:4:5 w e  lc 3 notitle,'L_4/data_4.dat' u 1:4 w lp lc 4 title 'L=4', 'L_4/data_4.dat' u 1:4:5 w e  lc 4 notitle, 'cv.dat' u 1:2 w lp lc rgb "#FF7F50" title 'ED L=1'
