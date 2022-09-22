module linearsolvers
      implicit none
	     real*8, parameter  :: NEARZERO   = 0.000001
      contains

!================================================================================
	subroutine sparse_directsolve(mat,b,n,x,success)

		!can be written better to use column matrix elements
		!rather than rows.

      		integer, intent(in)    :: n
	     	real*8, intent(inout)  :: mat(n,n)
		real*8, intent(inout)  :: b(n)
		real*8, intent(out)    :: x(n)
		logical,intent(out)    :: success

		integer :: i,j,k
		real*8 :: fac
		real*8 :: pivot_temp(n)
		real*8 :: bpivot_temp
		real*8 :: diag,sumoffdiag

		logical :: nonzeroexists

		success=.false.

		!loop over diagonal values
		do i=1,n-1

			success=.true.

			if(abs(mat(i,i)) .lt. NEARZERO) then

				nonzeroexists=.false.
				do j=i+1,n

					if(abs(mat(j,i)) .gt. NEARZERO) then
						nonzeroexists=.true.
						exit
					endif

				enddo

				if(nonzeroexists .eqv. .true.) then
					
					pivot_temp = mat(i,:)
					mat(i,:)   = mat(j,:)
					mat(j,:)   = pivot_temp

					bpivot_temp = b(i)
					b(i)        = b(j)
					b(j)        = bpivot_temp
				else
					!matrix determinant is 0
					success=.false.
					exit
				endif

			endif

			do j=i+1,n
				
				if(abs(mat(j,i)) .gt. NEARZERO) then

					fac = mat(j,i)/mat(i,i)
						
					do k=1,n
						if((abs(mat(j,k)) .gt. NEARZERO) .or. (abs(mat(i,k)) .gt. NEARZERO)) then
							mat(j,k)=mat(j,k)-fac*mat(i,k)
						endif
					enddo	
					
					b(j)=b(j) - fac*b(i)
				endif
			enddo
		enddo	

		if(n .eq. 1) then
			success=.true.
		endif
		!solve for x
		if(success .eqv. .true.) then

			do i=n,1,-1
				
				diag = mat(i,i)
				sumoffdiag = 0.d0

				do j=i+1,n
					sumoffdiag = sumoffdiag + mat(i,j)*x(j)
				enddo

				x(i) = (b(i)-sumoffdiag)/diag

			enddo

		endif
				

	end subroutine sparse_directsolve
!================================================================================
      
end module linearsolvers
