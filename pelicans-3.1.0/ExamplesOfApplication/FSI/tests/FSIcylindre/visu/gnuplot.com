plot "reference.txt"  using 1:($2+2) with linespoints \
      pointtype 1 title 'exact', \
     "numeric.txt" using 1:($2+2) with points \
      pointtype 9 pointsize 1.25 title 'numeric.txt'
set xlabel 't [s]'
set ylabel 'FS interface displacement [m]'