program mgrid

	use solve_manager
	implicit none

	integer :: ierr

	call MPI_Init(ierr)
	call solversetup()

	call timestepping()
	call MPI_Finalize(ierr)

end program mgrid
