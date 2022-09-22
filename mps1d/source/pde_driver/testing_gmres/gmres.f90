module solvergmres

  use vectorfunctions
  implicit none

contains

  !===================================================
  subroutine arnoldialgorithm(j,maxkspdim,n,Hmat,&
                              kspvecs,findAX,precond,&
			      lucky,luckykspdim,nanflag)

    integer,intent(in)    :: j,n,maxkspdim
    real*8, intent(inout) :: Hmat(maxkspdim+1,maxkspdim)
    real*8, intent(inout) :: kspvecs(n,maxkspdim+1)
    logical,intent(inout) :: lucky
    logical,intent(inout) :: nanflag
    integer,intent(inout) :: luckykspdim

    external :: findAX
    external :: precond

    real*8 :: Avj(n)
    real*8 :: MinvAvj(n)
    real*8 :: wj(n)
    real*8 :: vi(n),vj(n)
    
    integer :: i


    Avj     = ZERO
    MinvAvj = ZERO
    wj      = ZERO
    vi      = ZERO
    vj      = ZERO

    lucky   = .false.
    nanflag = .false.
		  
    vj = kspvecs(:,j)
    call findAX(Avj,vj,n)
    call precond(MinvAvj,Avj,n)
    do i=1,j
	vi = kspvecs(:,i)	
	call innerproduct(Hmat(i,j),MinvAvj,vi,n)
    enddo

    wj = MinvAvj

    do i=1,j
	vi = kspvecs(:,i)
	wj=wj - Hmat(i,j)*vi
    enddo

    call findnorm(Hmat(j+1,j),wj,n)

    if(Hmat(j+1,j) > VERYSMALL) then
	kspvecs(:,j+1)=wj(:)/Hmat(j+1,j)

    else if(isit_nan(Hmat(j+1,j)) .eqv. .false.) then
		print *,"hmat:",Hmat(j+1,j)
		lucky=.true.
                luckykspdim=j;
    else
        print *,"NaN or Inf detected:",Hmat(j+1,j)
        nanflag=.true.;
    endif

  end subroutine arnoldialgorithm
  !===================================================
  subroutine leastsquaresminimize(y,Hmat,m,maxkspdim,beta)

		  integer,intent(in)   :: m
                  integer,intent(in)   :: maxkspdim
		  real*8,intent(inout) :: Hmat(maxkspdim+1,maxkspdim)
		  real*8,intent(out)   :: y(maxkspdim)
		  real*8,intent(in)    :: beta

		  real*8 :: c,s,h_up,h_down,dtr
		  real*8 :: val1,val2
		  real*8 :: beta_e1(m+1)
                  real*8 :: Hmatcopy(maxkspdim+1,maxkspdim)

		  integer :: i,j

                  Hmatcopy   = Hmat
		  beta_e1    = ZERO
		  beta_e1(1) = beta;

		  do i=1,m

		  	h_up   = Hmatcopy(i,i)
			h_down = Hmatcopy(i+1,i)

			dtr = sqrt(h_up*h_up+h_down*h_down)

			c = h_up/dtr
			s = h_down/dtr

			do j=1,m
				
				h_up   = Hmatcopy(i,j)
				h_down = Hmatcopy(i+1,j)

				Hmatcopy(i,j)   =  c*h_up  + s*h_down
				Hmatcopy(i+1,j) = -s*h_up  + c*h_down 

			enddo

			val1 =  c*beta_e1(i)  + s*beta_e1(i+1); 
			val2 = -s*beta_e1(i)  + c*beta_e1(i+1);

			beta_e1(i)   = val1
			beta_e1(i+1) = val2	
		enddo


		y(m) = beta_e1(m)/Hmatcopy(m,m)

		do i=m-1,1,-1

			y(i)=beta_e1(i)

			do j=i+1,m
				y(i)=y(i)-Hmatcopy(i,j)*y(j)
			enddo

			y(i) = y(i)/Hmatcopy(i,i)

		enddo

  end subroutine leastsquaresminimize
  !===================================================
  subroutine performgmres(b,x0,x,maxkspdim,n,nrestarts,findAX,&
			precond,tol,success,printflag,initial_res)

	external :: findAX, precond

	integer,intent(in)    :: maxkspdim,n,nrestarts
	real*8, intent(in)    :: b(n)
        real*8, intent(inout) :: x0(n)
	real*8, intent(inout) :: x(n)
        real*8, intent(in)    :: tol
        logical, intent(out)  :: success
        logical, intent(in)   :: printflag

        real*8, intent(out)   :: initial_res

	integer :: i,j

	real*8 :: r0(n)
	real*8 :: Minvr(n)
	real*8 :: Ax0(n),Ax(n)
	real*8 :: r(n),v1(n)
	real*8 :: beta,residnorm
	real*8 :: y(maxkspdim)

	real*8 :: kspvectors(n,maxkspdim+1)
	real*8 :: Hmat(maxkspdim+1,maxkspdim)

	real*8 :: eps

	logical :: lucky
        logical :: nanflag
        integer :: luckykspdim
        integer :: kspdim,optkspdim

	Hmat       = ZERO
	kspvectors = ZERO
        r0         = ZERO
        Minvr      = ZERO
        Ax0        = ZERO
        Ax         = ZERO
        r          = ZERO
        v1         = ZERO
        y          = ZERO
	eps        = tol

	call findAX(Ax0,x0,n)
	r0 = b-Ax0

	call findnorm(initial_res,r0,n)

	call precond(Minvr,r0,n)
	r  = Minvr
	x  = x0


	lucky      = .false.
        nanflag    = .false.
        success    = .true.

	call findnorm(beta,r,n)
        i = 0
        if(printflag .eqv. .true.) print *,"restart iteration:",i,"residual norm:",beta

	do i=1,nrestarts
		  	
	    call findnorm(beta,r,n)
            x0 = x
	    v1 = r/beta
            kspvectors(:,1) = v1(:)

            do kspdim=1,maxkspdim

                x = x0

  		call arnoldialgorithm(kspdim,maxkspdim,n,&
			Hmat,kspvectors,findAX,precond,&
			lucky,luckykspdim,nanflag)

                if(nanflag .eqv. .true.) then
                    success=.false.
                    exit
                endif
  		
                call leastsquaresminimize(y,Hmat,kspdim,maxkspdim,beta)

		do j=1,kspdim
			x=x+y(j)*kspvectors(:,j)
		enddo

		call findAX(Ax,x,n)
		r = b - Ax
		call precond(Minvr,r,n)
		r = Minvr
	        call findnorm(residnorm,r,n)

                print *,"residnorm, kspdim:",residnorm,kspdim

		if(lucky .eqv. .true.) then
                    if(printflag .eqv. .true.) print *,"lucky condition"
                    exit
		endif

		if((residnorm .le. eps) .and. (kspdim .gt. 1)) then
			exit
		endif

	    enddo
	    optkspdim=merge(maxkspdim,kspdim,(kspdim .gt. maxkspdim))

            if((residnorm .le. eps) .or. (lucky .eqv. .true.) .or. (success .eqv. .false.)) then
                if(printflag .eqv. .true.) print *,"restart iteration:",i,"residual norm:",residnorm,"ksp size:",optkspdim
                exit
            endif
        
            if(printflag .eqv. .true.) print *,"restart iteration:",i,"residual norm:",residnorm,"ksp size:",optkspdim

        enddo
        
				
  end subroutine performgmres
  !===================================================
  function isit_nan(x) result(isnan)

          real*8 :: x
          logical :: isnan

          isnan = .false.
          if(x .ne. x) then
              isnan = .true.
          endif

          if((x/x .ne. 1.0) .and. (x .ne. 0)) then
              isnan = .true.
          endif

  end function
  !===================================================

end module solvergmres
