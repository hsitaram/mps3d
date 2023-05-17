module convdiff
use vectorfunctions
use solvergmres_module
implicit none

integer,parameter :: minsize=3

contains
!===============================================================
subroutine findupwindconvflux(cL,cR,uL,uR,flux)

	real*8,intent(in) :: cL,cR
	real*8,intent(in) :: uL,uR
	real*8,intent(out) :: flux

	real*8 :: c_half

	c_half=0.5*(cL+cR)
	
	if(c_half .ge. 0.d0) then
	    flux=c_half*uL
	else
            flux=c_half*uR
	endif
	     
end subroutine findupwindconvflux
!===============================================================
subroutine finddiffusiveflux(DL,DR,uL,uR,dx,flux)

	real*8,intent(in) :: DL,DR
	real*8,intent(in) :: uL,uR
	real*8,intent(in) :: dx
	real*8,intent(out) :: flux

	real*8 :: D_half

	D_half = HALF*(DL+DR)
	flux = D_half*(uR-uL)/dx

end subroutine finddiffusiveflux
!===============================================================
subroutine findAX(AX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,n)

      integer, intent(in) :: n
      real*8, intent(inout) :: AX(n)
      real*8, intent(in)  :: X(n)

      real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
      real*8,intent(in)  :: dx,dt
      logical,intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
      real*8, intent(in) :: dircvals(2),fluxvals(2)
      real*8, intent(in) :: timederivfactor
      
      integer :: i
      real*8 :: flux
      real*8 :: dx2

      dx2 = dx*dx

      AX(:) = timederivfactor*X(:)/dt
      
      !convection and diffusion terms
      do i=1,n-1
      	call findupwindconvflux(vel(i),vel(i+1),X(i),X(i+1),flux)
	AX(i)   = AX(i)   + flux/dx
	AX(i+1) = AX(i+1) - flux/dx
	call finddiffusiveflux(dcoeff(i),dcoeff(i+1),X(i),X(i+1),dx,flux)
	AX(i)   = AX(i)   - flux/dx
	AX(i+1) = AX(i+1) + flux/dx
      enddo

      !print *,"AX inside findAX:",AX

      !reaction term
      do i=1,n
      	AX(i) = AX(i) - reac(i)*X(i)
      enddo

      !boundary conditions
      if(dirc_bc_flags(1) .eqv. .true.) then
	      AX(1)=X(1)
      endif
      if(dirc_bc_flags(2) .eqv. .true.) then
	     AX(n)=X(n)
      endif

end subroutine findAX
!===============================================================
subroutine findrhs(b,xold,timederivfactor,&
		source,dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n)

      integer, intent(in) :: n
      real*8, intent(inout) :: b(n)
      real*8, intent(in) :: xold(n),source(n)
      logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
      real*8, intent(in) :: dircvals(2),fluxvals(2)
      real*8, intent(in) :: dx,dt
      real*8, intent(in) :: timederivfactor
      integer :: i

      do i=1,n
      	b(i) = timederivfactor*xold(i)/dt + source(i)
      enddo	

      if(dirc_bc_flags(1) .eqv. .true.) then
	      b(1) = dircvals(1)
      endif
      if(dirc_bc_flags(2) .eqv. .true.) then
	      b(n) = dircvals(2)
      endif
      if(flux_bc_flags(1) .eqv. .true.) then
	      b(1) = b(1) + fluxvals(1)/(dx)
      endif
      if(flux_bc_flags(2) .eqv. .true.) then
	      b(n) = b(n) - fluxvals(2)/(dx)
      endif

end subroutine findrhs
!===============================================================
subroutine noprecond(MinvX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
				flux_bc_flags,dircvals,fluxvals,&
				dx,dt,n)
      integer, intent(in) :: n
      real*8, intent(inout) :: MinvX(n)
      real*8, intent(in)  :: X(n)
	
      logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
      real*8, intent(in) :: dircvals(2),fluxvals(2)
      real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
      real*8,intent(in)  :: dx,dt
      real*8,intent(in)  :: timederivfactor

      MinvX = X

end subroutine noprecond
!===============================================================
subroutine gauss_seidel_smoothing(res,b,X,timederivfactor,vel,dcoeff,reac,&
		dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n,maxiter,tol)

	integer, intent(in)   :: n
	real*8, intent(in)    :: b(n)
	real*8, intent(inout) :: X(n)
	real*8, intent(inout) :: res(n)
        real*8, intent(in)    :: tol
	integer, intent(in) :: maxiter

        logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in) :: dircvals(2),fluxvals(2)
        real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
        real*8,intent(in)  :: dx,dt
	real*8,intent(in)  :: timederivfactor

	integer :: i,it

	real*8 :: diag
	real*8 :: cL,cR,chalf
	real*8 :: dL,dR,dhalf
	real*8 :: offdiag
	real*8 :: AX(n)
	real*8 :: dx2
        real*8 :: resnorm

	dx2 = dx*dx
        AX  = ZERO

	!gauss seidel iterations
	do it=1,maxiter
	     
	     do i=1,n

	        diag = -reac(i) + timederivfactor*ONE/dt

		offdiag=ZERO

		!right face
		if(i .lt. n) then
			!convection term
			cL = vel(i)
			cR = vel(i+1)
			chalf = HALF*(cL+cR)
			if(chalf .ge. 0) then
				diag = diag + chalf/dx
			else
				offdiag = offdiag + chalf*X(i+1)/dx
			endif
		
			!diffusion term
			dL = dcoeff(i)
			dR = dcoeff(i+1)
			dhalf = HALF*(dL + dR)
		
			diag = diag + dhalf/dx2
			offdiag = offdiag - X(i+1)*dhalf/dx2
		endif
	     	
	     	!left face
		if(i .gt. 1) then

			!convection term
			cL = vel(i-1)
			cR = vel(i)
			chalf = HALF*(cL+cR)
			if(chalf .ge. 0) then
				offdiag = offdiag - chalf*X(i-1)/dx
			else
				diag = diag - chalf/dx
			endif
		
			!diffusion term
			dL = dcoeff(i-1)
			dR = dcoeff(i)
			dhalf = HALF*(dL+dR)
		
			diag = diag + dhalf/dx2
			offdiag = offdiag - X(i-1)*dhalf/dx2
		endif

		if(i .eq. 1) then
			if(dirc_bc_flags(1) .eqv. .true.) then
				diag    = ONE
				offdiag = ZERO
			endif
		endif	

		if(i .eq. n) then
			if(dirc_bc_flags(2) .eqv. .true.) then
				diag    = ONE
				offdiag = ZERO
			endif
		endif

		X(i) = (b(i)-offdiag)/diag

		enddo
	
                call  findAX(AX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,n)

	        res = b - AX

                call findnorm(resnorm,res,n)

                if(resnorm .lt. tol) then
                   !print *,"it:",it
                   exit
                endif 

	enddo

	call  findAX(AX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,n)

	res = b - AX 

end subroutine gauss_seidel_smoothing
!===============================================================
subroutine restriction(Xh,X2h,n2h) 

	integer, intent(in) :: n2h
	real*8, intent(inout) :: Xh(2*n2h-1)
	real*8, intent(inout) :: X2h(n2h)

	integer :: i

	X2h(1)     = Xh(1)
	X2h(n2h) = Xh(2*n2h-1)

	do i=2,n2h-1
	   X2h(i)=ONEFOURTH*(Xh(2*i-2) + TWO*Xh(2*i-1) + Xh(2*i))
	enddo


end subroutine restriction
!===============================================================
subroutine prolong(Xh,X2h,n2h)

	integer, intent(in) :: n2h
	real*8, intent(inout) :: Xh(2*n2h-1)
	real*8, intent(inout) :: X2h(n2h)

	integer :: i

	do i=1,n2h-1
	   Xh(2*i-1) = X2h(i)
	   Xh(2*i)   = HALF*(X2h(i) + X2h(i+1))
	enddo

	Xh(2*n2h-1) = X2h(n2h)
	
end subroutine prolong
!===============================================================
recursive subroutine dovcycle(X,b,timederivfactor,vel,dcoeff,reac,&
		dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n)

	integer, intent(in)   :: n
	real*8, intent(inout) :: X(n)
	real*8, intent(inout) :: b(n)
        
	logical, intent(in)  :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in)   :: dircvals(2),fluxvals(2)
        real*8, intent(in)   :: vel(n),dcoeff(n),reac(n)
        real*8 ,intent(in)   :: dx,dt
	real*8 ,intent(in)   :: timederivfactor

	real*8  :: resh(n)
	real*8 :: eh(n)
	real*8 :: velh(n)
	real*8 :: dcoeffh(n)
	real*8 :: reach(n)
	
	real*8 :: res2h(n/2+1)
	real*8 :: e2h(n/2+1)
	real*8 :: vel2h(n/2+1)
	real*8 :: dcoeff2h(n/2+1)
	real*8 :: reac2h(n/2+1)
        real*8 :: Xmin(minsize)
        real*8 :: AX(minsize)
        logical :: success
        logical :: printflag
        real*8 :: initial_res

        success = .false.
        printflag = .false.


	velh    = vel
	dcoeffh = dcoeff
	reach   = reac

        resh     = ZERO
        eh       = ZERO

        res2h    = ZERO
        e2h      = ZERO
        vel2h    = ZERO
        dcoeff2h = ZERO
        reac2h   = ZERO
        Xmin     = ZERO

        if(n .le. minsize) then
            !call gauss_seidel_smoothing(resh,b,X,timederivfactor,velh,dcoeffh,reach,&
            !		dirc_bc_flags,flux_bc_flags,dircvals,&
            !		fluxvals,dx,dt,n,10000,NEARZERO)

            call performgmres(b,Xmin,X,timederivfactor,&
                velh,dcoeffh,reach,dirc_bc_flags,&
                flux_bc_flags,dircvals,fluxvals,dx,dt,&
                n,n,1,findAX,noprecond,&
                1.d-12,success,printflag,initial_res)

            AX = ZERO

            call  findAX(AX,X,timederivfactor,velh,dcoeffh,reach,dirc_bc_flags,&
                flux_bc_flags,dircvals,fluxvals,dx,dt,n)

            resh = b - AX 

        else

            !initial smoothing
            call gauss_seidel_smoothing(resh,b,X,timederivfactor,velh,dcoeffh,reach,&
                dirc_bc_flags,flux_bc_flags,dircvals,&
                fluxvals,dx,dt,n,1,NEARZERO)

            !restriction of residual from fine to coarse grid
            call restriction(resh,res2h,n/2+1)

            call restriction(reach,reac2h,n/2+1)
            call restriction(dcoeffh,dcoeff2h,n/2+1)
            call restriction(velh,vel2h,n/2+1)

            call dovcycle(e2h,res2h,timederivfactor,vel2h,dcoeff2h,reac2h,&
                dirc_bc_flags,flux_bc_flags,dircvals,&
                fluxvals,TWO*dx,dt,n/2+1)

            !prolong error from coarse to fine grid
            call prolong(eh,e2h,n/2+1)

            !update
            X(:) = X(:) + eh(:)

            !post smooth
            call gauss_seidel_smoothing(resh,b,X,timederivfactor,velh,dcoeffh,reach,&
                dirc_bc_flags,flux_bc_flags,dircvals,&
                fluxvals,dx,dt,n,1,NEARZERO)

        endif

    end subroutine dovcycle
    !===============================================================
    subroutine mgridprecond(MinvX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
            flux_bc_flags,dircvals,fluxvals,&
            dx,dt,n)

        integer, intent(in)    :: n
        real*8, intent(inout)    :: MinvX(n)
        real*8, intent(inout)  :: X(n)

        logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in) :: dircvals(2),fluxvals(2)
        real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
        real*8,intent(in)  :: dx,dt
        real*8,intent(in)  :: timederivfactor

        integer :: nvcycles,i
        !real*8 :: AX(n),res(n)
        !real*8 :: norm
        !integer :: j

        !norm=0.d0
        nvcycles=1

        do i=1,nvcycles
            call dovcycle(MinvX,X,timederivfactor,vel,dcoeff,reac,&
                dirc_bc_flags,flux_bc_flags,dircvals,&
                fluxvals,dx,dt,n)
        enddo


    end subroutine mgridprecond
    !===============================================================

end module convdiff
