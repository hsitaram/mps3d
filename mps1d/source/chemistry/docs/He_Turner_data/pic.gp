set terminal pdf
set output "reactrates.pdf"
set termoption font "Helvetica,20"
set key spacing 1.2
#set lmargin 8
#set rmargin 3
#set bmargin 4
set xlabel "Mean energy (eV)"
set ylabel "Rate coefficient (m^3/s)"
set xrange [1.0:100]
set yrange [1e-30:1e-1]
set log
plot 'fitting_24.58eV.dat' u 1:2 w lp lc 7 pt 6 pi 5 lw 2 lt 1 title "He + e -> He+ + 2e",\
     'fitting_24.58eV.dat' u 1:(70e-15*exp(-74/$1)) w l lc 8 lw 2 lt 9 title "",\
     'fitting_19.00eV.dat' u 1:2 w lp lc 8 lw 2 pt 6 pi 5 lt 1 title "He + e -> He* + e",\
     'fitting_19.00eV.dat' u 1:(10e-15*exp(-48/$1)) w l lc 8 lw 2 lt 9 title "",\
     'fitting_19.80eV.dat' u 1:2 w lp lc 9 lw 2 pt 6 pi 5 lt 1 title "He + e -> He**+ + e",\
     'fitting_19.80eV.dat' u 1:(10e-15*exp(-52/$1)) w l lc 8 lw 2 lt 9 title "Fit"
