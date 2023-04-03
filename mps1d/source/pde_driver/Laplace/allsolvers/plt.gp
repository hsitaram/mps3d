set key left
plot 'gmres_pdesoln.dat' u 1:2 w lp pt 4 pi 250 ps 2 title "gmres mgrid",\
'mgrid_pdesoln.dat' u 1:2 w lp pt 5 ps 2 pi 200 title "mgrid",\
'gs_pdesoln.dat' u 1:2 w lp lw 3 pt 5 ps 2 pi 205 title "gauss seidel",\
'gs_pdesoln.dat' u 1:1 w l pt 7 ps 2 title "exact"
