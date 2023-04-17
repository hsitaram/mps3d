<<<<<<< HEAD
CHEMFILE=TST_He_N2_chem
=======
CHEMFILE=He_chem
>>>>>>> mps3d_origin/main
FILE=../fundconstants.f90 ${CHEMFILE}.f90
OBJFILE=${CHEMFILE}.o fundconstants.o

CHEMDRIVER: 
	gfortran -c ${FILE}
	gfortran ${OBJFILE} chem_driver.f90 -o chemdriver
clean:
	rm -rf *.o *.mod chemdriver
