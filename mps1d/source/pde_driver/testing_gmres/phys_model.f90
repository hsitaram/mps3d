module phys_model
	implicit none
	contains

!===============================================
subroutine poisson_AX(AX,X,n)

	integer,intent(in) :: n
	real*8,intent(in)  :: X(n)
	real*8,intent(out) :: AX(n)

	integer :: i

	AX(1) = X(1)
	AX(n) = X(n)

	do i=2,n-1
		AX(i)=X(i-1)+X(i+1)-2*X(i)
	enddo
	
end subroutine poisson_AX
!===============================================
subroutine no_precond(MinvX,X,n)

	integer,intent(in) :: n
	real*8,intent(in)  :: X(n)
	real*8,intent(out) :: MinvX(n)

	integer :: i

	MinvX = X
	
end subroutine no_precond
!===============================================
subroutine gj_precond(MinvX,X,n)

	integer,intent(in) :: n
	real*8,intent(in)  :: X(n)
	real*8,intent(out) :: MinvX(n)

	integer :: i
	
	MinvX(1) = X(1)

	do i=2,n-1
		MinvX(i) = X(i)/(-2.0)
	enddo
	
	MinvX(n) = X(n)
	
end subroutine gj_precond
!===============================================
subroutine gs_precond(MinvX,X,n)

	integer,intent(in) :: n
	real*8,intent(in)  :: X(n)
	real*8,intent(out) :: MinvX(n)

	integer :: i,niter,it

	niter = 1	
	MinvX    = 0.d0	
	MinvX(1) = X(1)
	MinvX(n) = X(n)

	do it=1,niter
	do i=2,n-1
		MinvX(i) = (X(i)-MinvX(i-1)-MinvX(i+1))/(-2.0)
	enddo
	enddo
	
	
end subroutine gs_precond
!===============================================
end module phys_model
