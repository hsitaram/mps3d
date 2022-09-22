set terminal x11
plot 'gmres_pdesoln.dat' u 1:2 w lp pt 4 pi 100 ps 2 title "gmres mgrid",\
'mgrid_pdesoln.dat' u 1:2 w lp pt 5 pi 200 ps 2 title "mgrid",\
'gs_pdesoln.dat' u 1:2 w lp lw 3 pt 5 pi 300 ps 2 title "gauss seidel",\
'gs_pdesoln.dat' u 1:($1-$1**2) w lp pt 7 pi 300 ps 2 title "exact"
