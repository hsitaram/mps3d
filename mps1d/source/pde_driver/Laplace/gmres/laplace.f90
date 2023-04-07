program generic_pde
    use convdiff
    use solvergmres_module
    implicit none

    integer, parameter :: np = 5

    real*8  :: xmin,xmax
    real*8  :: initial_res,final_res
    logical :: printflag
    integer :: i,j
    logical :: dirc_bc_flags(2),flux_bc_flags(2)
    real*8  :: dircvals(2),fluxvals(2)
    real*8  :: timederivfactor

    real*8  :: dcoeff(np),vel(np),reac(np),source(np)
    real*8  :: b(np),phi(np),phiold(np)
    real*8  :: res(np),AX(np)

    real*8  :: dt,dx,gmres_tol,mgrid_tol,norm,x
    integer :: itmax_restart,maxkspdim
    integer :: ngs_it,ngs_inner_it,nvcycles
    logical :: success

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
    dt     = 0.1

    !==================================================
    
    !PDE being solved
    !d2u/dx2 = 0
    !exact solution x

    !here c=0, k=1, r=0, s=0

    !boundary conditions
    !dirichlet u=0 at  x=0
    !dirichlet u=1 at  x=1
    !or
    !flux in positive x direction f=(cu-kdu/dx)=-1 at x=1

    timederivfactor = ZERO
    dirc_bc_flags(1) = .true.
    dirc_bc_flags(2) = .true.
    flux_bc_flags(1) = .false.
    flux_bc_flags(2) = .false.

    dircvals(1) = ONE
    dircvals(2) = TWO
    fluxvals(1) = ZERO
    fluxvals(2) = ZERO

    do i=1,np
        x = xmin + (i-1)*dx
        vel(i)    =  ZERO
        dcoeff(i) =  ONE
        reac(i)   =  ZERO
        source(i) =   ZERO
    enddo
    !==================================================

    phi = ZERO
    phiold = phi
    gmres_tol = 1e-12
    itmax_restart = 1 
    maxkspdim     = 2
    printflag     = .true.

    call findrhs(b,phiold,timederivfactor,source,dirc_bc_flags, &
        flux_bc_flags,dircvals,fluxvals,dx,dt,np)
    
    print *,"b:",b

    call performgmres(b,phiold,phi,timederivfactor,&
        vel,dcoeff,reac,dirc_bc_flags,&
        flux_bc_flags,dircvals,fluxvals,dx,dt,&
        maxkspdim,np,itmax_restart,findAX,noprecond,&
        gmres_tol,success,printflag,initial_res)

    print *,"phi:",phi

end program generic_pde
