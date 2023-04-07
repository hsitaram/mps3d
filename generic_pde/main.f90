program mgrid

	use solve_manager
	implicit none

	integer :: ierr

	call MPI_Init(ierr)
	call solversetup()

#ifdef LAPLACE_SOLVE
        print *,"in laplace solve"
	call timestepping(0.1d0,0.1d0)
#else
        print *,"not in laplace solve"
	call timestepping(0.1d0,1.d0)
#endif

	call MPI_Finalize(ierr)

end program mgrid
