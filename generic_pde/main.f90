program mgrid

	use solve_manager
	implicit none

	integer :: ierr

	call MPI_Init(ierr)
	call solversetup()

	call timestepping(0.1d0,0.1d0)
	call MPI_Finalize(ierr)

end program mgrid
