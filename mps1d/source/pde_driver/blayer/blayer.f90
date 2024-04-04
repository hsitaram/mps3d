program blayer
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
    real*8  :: l2err
    
    real*8, allocatable :: exactsoln(:)

    real*8,allocatable  :: dcoeff(:),vel(:),reac(:),source(:)
    real*8,allocatable  :: b(:),phi(:),phiold(:)
    real*8,allocatable  :: res(:),AX(:)
    integer :: fluxscheme
    real*8 :: const_vel,const_dcoeff

    print *,"enter no: of points (2^n+1) (e.g. 65)" 
    read(*,*) np
    print *,"enter fluxscheme (1-upwind,2-WAF,3-central,4-SG)"
    read(*,*) fluxscheme 
    print *,"enter vel, dcoeff: (e.g v=-1,dcoeff=0.1) note: do try negative velocities!!"
    read(*,*) const_vel,const_dcoeff


    allocate(dcoeff(np))
    allocate(vel(np))
    allocate(reac(np))
    allocate(source(np))
    allocate(b(np))
    allocate(phi(np))
    allocate(phiold(np))
    allocate(res(np))
    allocate(AX(np))
    allocate(exactsoln(np))


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
    !d/dx (c u)=d/dx( D (du/dx))
    !exact solution (exp(cx/D)-1)/(exp(c/D)-1)

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

    dircvals(1) = 1.d0
    dircvals(2) = 0.d0
    fluxvals(1) = 0.d0
    fluxvals(2) = 0.d0

    do i=1,np
        x = xmin + (i-1)*dx
        vel(i)    =  const_vel
        dcoeff(i) =  const_dcoeff
        reac(i)   =  0.d0
        source(i) =  0.d0
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
        exactsoln(i)=(exp(const_vel*x/const_dcoeff)-exp(const_vel/const_dcoeff))/(1.d0-exp(const_vel/const_dcoeff))
        l2err = l2err + (phi(i)-exactsoln(i))**2
    enddo
    
    call writesoln("exactsoln.dat",exactsoln,np,xmin,dx)
    print *,"==================================="    

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

end program blayer

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
