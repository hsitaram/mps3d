set xrange [0:1]
plot 'gmres_pdesoln.dat' u 1:2 w lp pt 4 pi 250 ps 2 title "gmres mgrid",\
(x**2+x)*exp(-1) title "exact"
