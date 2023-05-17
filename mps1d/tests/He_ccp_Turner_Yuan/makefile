CHEMFILEPREFIX=TST_He_ccp_Turner_Yuan_chem
CHEMOBJFILE = ${CHEMFILEPREFIX}.o
SOURCEDIR = ../source
CHEMDIR  = ${SOURCEDIR}/chemistry

OBJFILES = ${SOURCEDIR}/fundconstants.o ${CHEMDIR}/${CHEMOBJFILE} ${SOURCEDIR}/vectorfunctions.o\
	   ${SOURCEDIR}/gmres.o  ${SOURCEDIR}/convdiff.o ${SOURCEDIR}/solver_zdp.o
MAINFILE = ${SOURCEDIR}/main.f90

FORT = gfortran
#FORTFLAGS = -Wunused -fbounds-check
FORTFLAGS = -O3
#FORTFLAGS = 

%.o: %.f90
	${FORT} -c ${FORTFLAGS} $< -o $@

fortranconvdiff: ${OBJFILES}
	${FORT} ${FORTFLAGS} ${OBJFILES} ${MAINFILE} -o plsolver_${CHEMFILEPREFIX}.exe

clean:
	rm -rf *.o ${SOURCEDIR}/*.o ${CHEMDIR}/*.o *.mod *.exe *.dat
