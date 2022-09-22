program linsolve_driver

      use linearsolvers
      implicit none

      integer, parameter :: n=50
      real*8,allocatable :: mat(:,:)
      real*8,allocatable :: b(:)
      real*8,allocatable :: x(:)
      logical :: success
      integer :: i


      allocate(mat(n,n))
      allocate(b(n))
      allocate(x(n))

      mat  = 0.d0
      b(1) = 1.d0
      b(n) = 2.d0
      x=0.d0

      do i=1,n
      	mat(i,i)=2.d0

	if(i .gt. 1) then
		mat(i,i-1)=-1.d0
	endif

	if(i .lt. n) then
		mat(i,i+1)=-1.d0
	endif
      enddo

      !b(3) = 1.d0
      !b(n) = 2.d0
      !mat(1,1) =  0.d0
      !mat(2,1) = -1.d0
      !mat(3,1) =  2.d0
      !mat(4,1) =  1.d0

      !mat(1,2) = -1.d0
      !mat(3,2) = -1.d0

      !mat(1,3) =  1.d0
      !mat(2,3) = -1.d0
      !mat(3,3) = -1.d0
      !mat(4,3) = -1.d0

      !mat(1,4) = -1.d0


      !print *,mat(1,:)
      !print *,mat(2,:)
      !print *,mat(3,:)
      !print *,mat(4,:)
      !print *,"============="

      call sparse_directsolve(mat,b,n,x,success)
      
      !print *,mat(1,:)
      !print *,mat(2,:)
      !print *,mat(3,:)
      !print *,mat(4,:)
      !print *,"============="
      
      print *,"x:",x

      if(success .eqv. .false.) print *,"failed"

end program linsolve_driver
