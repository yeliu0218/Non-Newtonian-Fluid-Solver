#
# run gnuplot in the directory containing the result files
# > gnuplot visu/error.plot
#

out = 0

if( out==0 ) set terminal x11 ;
if( out==0 ) set size square 1.0 ;

if( out==1 ) set terminal postscript eps colour "Helvetica" ;
if( out==1 ) set output "error.eps" ;
if( out==1 ) set size square 0.8 ;

if( out==2 ) set terminal svg fname "Bitstream Vera Sans Mono" fsize 18 ;

set grid

set mxtics 10
set mytics 10

#set xrange [-1.0:1.0]
#set yrange [-0.2:1.4]

#--------------------------------------
# error norms
#--------------------------------------

set xlabel "h"
set ylabel "error"
set key top left

set logscale x
set logscale y

# manually enter the L2 solution norm
sol_L2_norm = 1.0
# manually enter the H1 solution norm
sol_H1_norm = 1.0 # 

plot "error_norms.txt" u 3:($6/sol_L2_norm) title "L2 error norm" w lp lw 2 ,\
     "error_norms.txt" u 3:($8/sol_H1_norm) title "H1 error norm" w lp lw 2

if( out==0 ) pause mouse
