program main

      use plasma_solver
      implicit none

      real :: time_start, time_end

      call init()
      call cpu_time(time_start)
      call timestepping()
      call cpu_time(time_end)
      write(*,*) "Time taken by timestepping = ", time_end - time_start 

end program main
