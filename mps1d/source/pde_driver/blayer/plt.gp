set xrange [0:1]
plot 'gmres_pdesoln.dat' u 1:2 w lp pt 6 pi 5 ps 2 title "gmres mgrid",\
     'exactsoln.dat' u 1:2 w lp pt 4 pi 7 ps 2 title "exact"
