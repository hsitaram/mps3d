#CHEMFILEPREFIX = H2_O2_chem
CHEMFILEPREFIX = ion_elec_chem
CHEMFILE = ${CHEMFILEPREFIX}.f90
CHEMOBJFILE = ${CHEMFILEPREFIX}.o

MGRIDLIBDIR = ../mgridlib
SOLVERDIR = ../plasma_pde
CHEMDIR  = ${SOLVERDIR}/chemistry

MGRIDLIBFILES = ${MGRIDLIBDIR}/globalvars.o ${MGRIDLIBDIR}/vectorfunctions.o \
		${MGRIDLIBDIR}/par_decompose.o ${MGRIDLIBDIR}/discretepde.o \
		${MGRIDLIBDIR}/cartwrite.o ${MGRIDLIBDIR}/linearsolvers.o ${MGRIDLIBDIR}/gmres.o \
		${MGRIDLIBDIR}/mgridsteps.o 

SOLVERFILES   = ${SOLVERDIR}/fundconstants.o \
		${CHEMDIR}/${CHEMOBJFILE} ${SOLVERDIR}/poissoneq.o \
		${SOLVERDIR}/edeneq.o \
		${SOLVERDIR}/ionseq.o \
		${SOLVERDIR}/neutralseq.o \
		${SOLVERDIR}/eenergyeq.o \
	       	${SOLVERDIR}/inputs.o \
		${SOLVERDIR}/solve_manager.o \

OBJFILES = globalvars.o vectorfunctions.o inputs.o par_decompose.o \
	   discretepde.o cartwrite.o linearsolvers.o gmres.o\
	   mgridsteps.o poissoneq.o \
	   edeneq.o ionseq.o neutralseq.o eenergyeq.o \
	   solve_manager.o ${CHEMDIR}/${CHEMOBJFILE}

FORT = mpif90
FFLAGS = -O2
#FFLAGS = -g -check bounds

EXECNAME = mps3d_${CHEMFILEPREFIX}.exe

%.o: %.f90
	${FORT} -c ${FFLAGS} $< -o $@

mgridsolve: ${MGRIDLIBFILES} ${SOLVERFILES}
	mv ${MGRIDLIBDIR}/*.o .
	mv ${SOLVERDIR}/*.o .
	${FORT} ${FFLAGS} ${SOLVERDIR}/main.f90 ${OBJFILES} -o ${EXECNAME}

clean:
	rm -rf ${MGRIDLIBDIR}/*.o ${SOLVERDIR}/*.o ${CHEMDIR}/*.o *.o *.mod *.exe mgridsolve
