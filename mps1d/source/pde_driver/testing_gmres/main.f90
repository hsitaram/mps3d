!===============================================
program testgmres
	use phys_model
        use solvergmres	
	implicit none
	 integer, parameter :: NP=64
	 real*8 :: b(NP),x0(NP),x(NP)
	 integer :: maxkspdim
	 integer :: nrestarts
	 real*8  :: tol,init_res
	 logical :: success

	  maxkspdim=64
	  nrestarts=1
	  tol=1e-30
	 

	 b=0.d0
	 x0=1.d0
	 x=1.50

	b(1) =1.d0
	b(NP)=2.d0

	call performgmres(b,x0,x,maxkspdim,NP,nrestarts,&
		poisson_AX,no_precond,tol,success,.true.,init_res)

	call printsoln(x,NP)

	
end program testgmres
!===============================================
subroutine printsoln(x,n)
	
	implicit none
	integer,intent(in) :: n
	real*8, intent(in) :: x(n)

	integer :: i,fptr

	fptr=666
	open(file='soln.dat',unit=fptr)

	do i=1,n
		write(fptr,*)i,x(i)
	enddo

	close(fptr)

end subroutine printsoln
!===============================================
