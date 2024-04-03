program mms
    use convdiff
    use solvergmres_module
    implicit none

    integer :: np

    real*8  :: xmin,xmax
    real*8  :: initial_res
    logical :: printflag
    integer :: i,tstep
    logical :: dirc_bc_flags(2),flux_bc_flags(2)
    real*8  :: dircvals(2),fluxvals(2)
    real*8  :: timederivfactor


    real*8  :: dt,dx,gmres_tol,mgrid_tol,norm,x
    integer :: itmax_restart,maxkspdim
    integer :: ngs_it,ngs_inner_it,nvcycles
    logical :: success
    real*8 :: t,fintime
    real*8  :: l2err,exactsoln

    real*8,allocatable  :: dcoeff(:),vel(:),reac(:),source(:)
    real*8,allocatable  :: b(:),phi(:),phiold(:)
    real*8,allocatable  :: res(:),AX(:)
    integer :: fluxscheme

    print *,"enter no: of points (2^n+1)" 
    read(*,*) np
    fluxscheme=3 !central

    allocate(dcoeff(np))
    allocate(vel(np))
    allocate(reac(np))
    allocate(source(np))
    allocate(b(np))
    allocate(phi(np))
    allocate(phiold(np))
    allocate(res(np))
    allocate(AX(np))


    dcoeff = ZERO
    vel    = ZERO
    reac   = ZERO
    source = ZERO
    phi    = ZERO
    b      = ZERO
    res    = ZERO
    AX     = ZERO

    xmin   = ZERO
    xmax   = ONE
    dx     = (xmax-xmin)/(np-1)
    dt      = 0.1

    !==================================================
    !general PDE is of the form
    ! du/dt + d/dx (cu) = d/dx (k du/dx) + ru + s
    !==================================================

    !==================================================
    !PDE being solved
    !d/dx (-u)=d/dx( 0.5 (du/dx)) + 2 u - x^2
    !exact solution 0.5(x^2-x)

    !boundary conditions
    !dirichlet u=0 at  x=0
    !dirichlet u=0 at  x=1
    !flux f=0.25 at  x=0
    !flux f=-0.25 at  x=1

    timederivfactor = ZERO
    dirc_bc_flags(1) = .true.
    dirc_bc_flags(2) = .true.
    flux_bc_flags(1) = .false.
    flux_bc_flags(2) = .false.

    dircvals(1) = ZERO
    dircvals(2) = ZERO
    fluxvals(1) = 0.25
    fluxvals(2) = -0.25

    do i=1,np
        x = xmin + (i-1)*dx
        vel(i)    =  -1.d0
        dcoeff(i) =  0.5
        reac(i)   =  2.d0
        source(i) =  -x**2
    enddo
    !==================================================

    phi = ZERO
    phiold = phi
    gmres_tol = 1e-9
    itmax_restart = 10 
    maxkspdim     = np
    printflag     = .true.

    call findrhs(b,phiold,timederivfactor,source,dirc_bc_flags, &
        flux_bc_flags,dircvals,fluxvals,dx,dt,np)

    call performgmres(b,phiold,phi,timederivfactor,&
        vel,dcoeff,reac,dirc_bc_flags,&
        flux_bc_flags,dircvals,fluxvals,dx,dt,&
        maxkspdim,np,itmax_restart,findAX,mgridprecond,&
        gmres_tol,success,printflag,initial_res,fluxscheme)


    call writesoln("gmres_pdesoln.dat",phi,np,xmin,dx)
    print *,"==================================="    

    l2err=0.d0    
    do i=1,np
        x = xmin + (i-1)*dx
        exactsoln=0.5d0*(x**2-x)
        l2err = l2err + (phi(i)-exactsoln)**2
    enddo

    l2err=sqrt(l2err/np)

    print *,dx,dt,l2err

    deallocate(dcoeff)
    deallocate(vel)
    deallocate(reac)
    deallocate(source)
    deallocate(b)
    deallocate(phi)
    deallocate(phiold)
    deallocate(res)
    deallocate(AX)

end program mms

subroutine writesoln(fname,soln,np,xmin,dx)

    implicit none
    character(LEN=*), intent(in) :: fname
    integer, intent(in) :: np
    real*8, intent(in)  :: xmin,dx
    real*8, intent(in) :: soln(np)

    integer :: fnum,i

    fnum=13

    open(unit=fnum,file=trim(fname))

    do i=1,np
        write(fnum,'(F20.8,F20.8)') xmin+(i-1)*dx,soln(i)
    enddo

    close(fnum)

end subroutine writesoln
