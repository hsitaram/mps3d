module vectorfunctions
      use fundconstants
      implicit none
      contains
  !===================================================
  subroutine findnorm(norm,v1,n)

    integer,intent(in) :: n
    real*8,intent(in) :: v1(n)
    real*8, intent(out) :: norm

    integer :: i

    norm = ZERO

    do i=1,n
       norm=norm+v1(i)*v1(i)
    enddo

    norm=sqrt(norm)

  end subroutine findnorm
  !===================================================
  subroutine innerproduct(v1dotv2,v1,v2,n)

    integer, intent(in) :: n
    real*8, intent(in) :: v1(n),v2(n)
    real*8, intent(out) :: v1dotv2

    integer :: i

    v1dotv2 = ZERO

    do i=1,n
       v1dotv2=v1dotv2+v1(i)*v2(i)
    enddo

  end subroutine innerproduct
  !===================================================

end module vectorfunctions
