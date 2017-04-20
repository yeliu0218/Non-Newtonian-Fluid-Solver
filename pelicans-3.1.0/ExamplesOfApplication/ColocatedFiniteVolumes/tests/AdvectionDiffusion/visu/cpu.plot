#
# run gnuplot in the directory containing the result files
# > gnuplot visu/cpu.plot
#

out = 0

if( out==0 ) set terminal x11 ;
if( out==0 ) set size square 1.0 ;

if( out==1 ) set terminal postscript eps colour "Helvetica" ;
if( out==1 ) set output "cutline.eps" ;
if( out==1 ) set size square 0.8 ;

if( out==2 ) set terminal svg fname "Bitstream Vera Sans Mono" fsize 18 ;

set mxtics 2 ### minor tics marks along the x axis
set mytics 2 ### minor tics marks along the x axis

set xlabel "x"
set ylabel "u"

set key top right ### legend on top right

plot "cutline_00002.txt" using 3:5 title "analytic" with lines \
                                                    linewidth 2, \
     "cutline_00002.txt" using 3:4 title "computed" with lines \
                                                    linewidth 2

if( out==0 ) pause mouse
