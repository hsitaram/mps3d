set terminal pdf
set output "mob.pdf"
set termoption font "Helvetica,20"
set key spacing 1.2
#set lmargin 8
#set rmargin 3
#set bmargin 4
set xlabel "Mean energy (eV)"
set ylabel "Mobility*N (SI units)"
set log
plot 'fitting_Mobility.dat' u 1:2 w lp lc 7 lw 2 lt 1 title "BOLSIG output",\
     'fitting_Mobility.dat' u 1:(exp(55.0+0.3942*log($1)+2.134/$1-0.6433/($1**2)+0.07112/($1**3))) w l lw 2 title "fit"
