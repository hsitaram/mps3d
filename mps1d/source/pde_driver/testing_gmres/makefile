FILES=vectorfunctions.o gmres.o phys_model.o

FORT=gfortran

%.o: %.f90
	${FORT} -c $< -o $@
gmrestest:${FILES}
	${FORT} ${FILES} main.f90 -o gmrestest
clean:
	rm -rf *.o *.mod ./gmrestest
