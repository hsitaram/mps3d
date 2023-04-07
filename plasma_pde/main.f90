program mgrid

	use solve_manager
	implicit none

	integer :: ierr
	real*8  :: time_stepping_wall_time

	call MPI_Init(ierr)
	call solversetup()
	
	if(g_myproc .eq. g_rootproc) time_stepping_wall_time = -MPI_WTime()
	call timestepping()
	if(g_myproc .eq. g_rootproc) time_stepping_wall_time = &
	time_stepping_wall_time + MPI_WTime()

	if(g_myproc .eq. g_rootproc) write(*,*) "Wall clock time:",&
        time_stepping_wall_time

	call MPI_Finalize(ierr)

end program mgrid
