program generic_pde
    use convdiff
    use solvergmres_module
    implicit none

    integer, parameter :: np = 4097

    real*8  :: xmin,xmax
    real*8  :: residual
    logical :: printflag
    integer :: i,j
    logical :: dirc_bc_flags(2),flux_bc_flags(2)
    real*8  :: dircvals(2),fluxvals(2)
    real*8  :: timederivfactor

    real*8  :: dcoeff(np),vel(np),reac(np),source(np)
    real*8  :: b(np),phi(np),phiold(np)
    real*8  :: res(np),AX(np)
    
    real*8  :: dt,dx,gmres_tol,norm,x
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
    
    !general PDE is of the form
    ! du/dt + d/dx (cu) = d/dx (k du/dx) + ru + s

    !PDE being solved
    !d2u/dx2 + u + x^2 - x + 2 = 0
    !exact solution x-x^2

    !here c=0, k=1, r=1, s=x^2-x+2

    !boundary conditions
    !dirichlet u=0 at  x=0

    !flux in positive x direction f=(cu-kdu/dx)=1 at x=1
    
    timederivfactor = ZERO
      	
    do i=1,np
            x = xmin + (i-1)*dx
	    vel(i)    =  ZERO
            dcoeff(i) =  ONE
      	    !reac(i)   =  ONE
      	    reac(i)   =  ZERO
	    !source(i) =  x**2 - x + TWO
	    source(i) =  ZERO
            phi(i)    =  ONE - x
    enddo

    dirc_bc_flags(1) = .true.
    !dirc_bc_flags(2) = .false.
    dirc_bc_flags(2) = .true.
    flux_bc_flags(1) = .false.
    !flux_bc_flags(2) = .true.
    flux_bc_flags(2) = .false.

    dircvals(1) = ZERO
    dircvals(2) = ZERO
    fluxvals(1) = ZERO
    fluxvals(2) = ONE

    print *,"gauss seidel======================"   
    phiold = phi
    ngs_it = 50 
    ngs_inner_it = 1000

    call findrhs(b,phiold,timederivfactor,source,dirc_bc_flags, &
			flux_bc_flags,dircvals,fluxvals,dx,dt,np)
    do i=1,ngs_it
		call gauss_seidel_smoothing(res,b,phi,timederivfactor,vel,dcoeff,reac,&
				dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,np,ngs_inner_it,NEARZERO)

                call findnorm(norm,res,np)
                print *,"iteration:",i,"residual norm=",norm
    enddo
    call writesoln("gs_pdesoln.dat",phi,np,xmin,dx)
    print *,"==================================="    
    
    print *,"multigrid=========================="    
    phiold   = phi
    nvcycles = 40

    call findrhs(b,phiold,timederivfactor,source,dirc_bc_flags, &
			flux_bc_flags,dircvals,fluxvals,dx,dt,np)
    do i=1,nvcycles
		call dovcycle(phi,b,timederivfactor,vel,dcoeff,reac,&
		dirc_bc_flags,flux_bc_flags,dircvals,&
			fluxvals,dx,dt,np)

                call findAX(AX,phi,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,np)

                res = b-AX
                call findnorm(norm,res,np)
                print *,"residual norm=",norm

    enddo
    call writesoln("mgrid_pdesoln.dat",phi,np,xmin,dx)
    print *,"==================================="    


    print *,"GMRES-mgrid========================"    
    phiold = phi
    gmres_tol = NEARZERO
    itmax_restart = 50
    maxkspdim     = 5
    printflag     = .true.

    call findrhs(b,phiold,timederivfactor,source,dirc_bc_flags, &
			flux_bc_flags,dircvals,fluxvals,dx,dt,np)

    call performgmres(b,phiold,phi,timederivfactor,&
    			vel,dcoeff,reac,dirc_bc_flags,&
    			flux_bc_flags,dircvals,fluxvals,dx,dt,&
    			maxkspdim,np,itmax_restart,findAX,mgridprecond,&
    			gmres_tol,success,printflag,residual)

    call writesoln("gmres_pdesoln.dat",phi,np,xmin,dx)
    print *,"==================================="    

end program generic_pde

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
        write(fnum,'(F10.5,F10.5)') xmin+(i-1)*dx,soln(i)
    enddo

    close(fnum)

end subroutine writesoln
