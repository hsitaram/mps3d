SOURCEDIR = ../../../

OBJFILES = ${SOURCEDIR}/fundconstants.o ${SOURCEDIR}/vectorfunctions.o ${SOURCEDIR}/gmres.o ${SOURCEDIR}/convdiff.o
MAINFILE = laplace.f90

COMPTYPE=GNU
#COMPTYPE=INTEL
ifeq (${COMPTYPE},GNU)
	FORT = gfortran
else
	FORT = ifort
endif
#FORTFLAGS = -Wunused -fbounds-check
FORTFLAGS = -O2
#FORTFLAGS = 

%.o: %.f90
	${FORT} -c ${FORTFLAGS} $< -o $@

fortranconvdiff: ${OBJFILES}
	${FORT} ${FORTFLAGS} ${OBJFILES} ${MAINFILE} -o gpde.exe

clean:
	rm -rf *.o ${SOURCEDIR}/*.o *.mod *.exe *.dat
