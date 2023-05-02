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
    integer :: ntsteps
    real*8  :: l2err,exactsoln
    
    real*8,allocatable  :: dcoeff(:),vel(:),reac(:),source(:)
    real*8,allocatable  :: b(:),phi(:),phiold(:)
    real*8,allocatable  :: res(:),AX(:)

    print *,"enter no: of points (2^n+1), no: of timesteps"    
    read(*,*) np,ntsteps

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
    fintime = 1.0
    t       = 0.d0
    dt      = fintime/ntsteps

    !==================================================
    !general PDE is of the form
    ! du/dt + d/dx (cu) = d/dx (k du/dx) + ru + s
    !==================================================

    !==================================================
    !PDE being solved
    !du/dt + d/dx (x u)=d/dx( (x^2+1) du/dx) - u - exp(-t) (3x^2+2)
    !exact solution (x^2+x)*exp(-t)

    !boundary conditions
    !dirichlet u=0 at  x=0
    !dirichlet u=2 exp(-t) at  x=1

    timederivfactor = ONE
    dirc_bc_flags(1) = .true.
    dirc_bc_flags(2) = .true.
    flux_bc_flags(1) = .false.
    flux_bc_flags(2) = .false.

    dircvals(1) = ZERO
    dircvals(2) = 2.d0*exp(-t)
    fluxvals(1) = ZERO
    fluxvals(2) = ZERO

    do i=1,np
        x = xmin + (i-1)*dx
        vel(i)    =  x
        dcoeff(i) =  x**2+ONE
        reac(i)   =  -1.d0
        source(i) =  -(3.d0*x**2+2.0)*exp(-t)
    enddo
    !==================================================
    
    phi = ZERO
    phiold = phi
    gmres_tol = 1e-9
    itmax_restart = 10 
    maxkspdim     = 10
    printflag     = .true.

    do tstep=1,ntsteps

        print *,"==================================="    
        do i=1,np
            x = xmin + (i-1)*dx
            source(i) =-(3.d0*x**2+2.0)*exp(-t)
        enddo

        dircvals(1) = ZERO
        dircvals(2) = 2.d0*exp(-t)

        call findrhs(b,phiold,timederivfactor,source,dirc_bc_flags, &
            flux_bc_flags,dircvals,fluxvals,dx,dt,np)

        call performgmres(b,phiold,phi,timederivfactor,&
            vel,dcoeff,reac,dirc_bc_flags,&
            flux_bc_flags,dircvals,fluxvals,dx,dt,&
            maxkspdim,np,itmax_restart,findAX,mgridprecond,&
            gmres_tol,success,printflag,initial_res)

        t=t+dt
        phiold=phi
        print *,"time,dt:",t,dt,"================="  

    enddo

    call writesoln("gmres_pdesoln.dat",phi,np,xmin,dx)
    print *,"==================================="    
    
    l2err=0.d0    
    do i=1,np
            x = xmin + (i-1)*dx
            exactsoln=(x**2+x)*exp(-fintime)
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
        write(fnum,'(F12.8,F12.8)') xmin+(i-1)*dx,soln(i)
    enddo

    close(fnum)

end subroutine writesoln
