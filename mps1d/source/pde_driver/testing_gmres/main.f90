!===============================================
program testgmres
    use phys_model
    use solvergmres	
    implicit none
    integer, parameter :: NP=10
    real*8 :: b(NP),x0(NP),x(NP)
    integer :: maxkspdim
    integer :: nrestarts
    real*8  :: tol,init_res
    logical :: success

    maxkspdim=3
    nrestarts=5
    tol=1e-10


    b=0.d0
    x0=0.d0
    x=0.d0

    b(1) =1.d0
    b(NP)=2.d0

    print *,"gmres with no preconditioning"
    call performgmres(b,x0,x,maxkspdim,NP,nrestarts,&
        poisson_AX,no_precond,tol,success,.true.,init_res)
    
    b=0.d0
    x0=0.d0
    x=0.d0

    b(1) =1.d0
    b(NP)=2.d0
    
    print *,"gmres with Gauss Seidel preconditioning"
    call performgmres(b,x0,x,maxkspdim,NP,nrestarts,&
        poisson_AX,gs_precond,tol,success,.true.,init_res)
    
    b=0.d0
    x0=0.d0
    x=0.d0

    b(1) =1.d0
    b(NP)=2.d0
    
    print *,"gmres with Gauss Jacobi preconditioning"
    call performgmres(b,x0,x,maxkspdim,NP,nrestarts,&
        poisson_AX,gj_precond,tol,success,.true.,init_res)

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
