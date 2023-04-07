module solvergmres_module

  use vectorfunctions
  implicit none

contains

  !===================================================
  subroutine arnoldialgorithm(j,maxkspdim,&
			      timederivflag,dt,lx,ly,lz,&
                              vel,dcoeff,reac,&
			      bc_codes,bcvals,&
                              lrank,rrank,brank,trank,krank,frank,&
                              llenx,lleny,llenz,Hmat,&
                              kspvecs,findAX,precond,&
			      lucky,luckykspdim,nanflag)

    logical, intent(in)   :: timederivflag
    integer,intent(in)    :: j,maxkspdim,lx,ly,lz
    character(LEN=4),intent(in)         :: bc_codes(NFACES)
    type(boundarycondition), intent(inout) :: bcvals

    real*8, intent(inout) :: Hmat(maxkspdim+1,maxkspdim)
    real*8, intent(inout) :: kspvecs(-g_nglayers+1:lx+g_nglayers,&
                                     -g_nglayers+1:ly+g_nglayers,&
                                     -g_nglayers+1:lz+g_nglayers,maxkspdim+1)

    real*8, intent(in)    :: vel(-g_nglayers+1:lx+g_nglayers,&
                                 -g_nglayers+1:ly+g_nglayers,&
                                 -g_nglayers+1:lz+g_nglayers)

    real*8, intent(in)    :: dcoeff(-g_nglayers+1:lx+g_nglayers,&
                                 -g_nglayers+1:ly+g_nglayers,&
                                 -g_nglayers+1:lz+g_nglayers)

    real*8, intent(in)    :: reac(-g_nglayers+1:lx+g_nglayers,&
                                 -g_nglayers+1:ly+g_nglayers,&
                                 -g_nglayers+1:lz+g_nglayers)

    real*8, intent(in)    :: dt,llenx,lleny,llenz
    integer,intent(in)    :: lrank,rrank,brank,trank,krank,frank
    logical,intent(inout) :: lucky
    logical,intent(inout) :: nanflag
    integer,intent(inout) :: luckykspdim

    external :: findAX
    external :: precond

    integer :: n, n_with_ghst

    real*8 :: Avj(-g_nglayers+1:lx+g_nglayers,&
                 -g_nglayers+1:ly+g_nglayers,&
                 -g_nglayers+1:lz+g_nglayers)

    real*8 :: MinvAvj(-g_nglayers+1:lx+g_nglayers,&
                 -g_nglayers+1:ly+g_nglayers,&
                 -g_nglayers+1:lz+g_nglayers)

    real*8 :: wj(-g_nglayers+1:lx+g_nglayers,&
                 -g_nglayers+1:ly+g_nglayers,&
                 -g_nglayers+1:lz+g_nglayers)

    real*8 :: vi(-g_nglayers+1:lx+g_nglayers,&
                 -g_nglayers+1:ly+g_nglayers,&
                 -g_nglayers+1:lz+g_nglayers)

    real*8 :: vj(-g_nglayers+1:lx+g_nglayers,&
                 -g_nglayers+1:ly+g_nglayers,&
                 -g_nglayers+1:lz+g_nglayers)
    
    integer :: i
    real*8 :: diag(lx,ly,lz) 
    real*8 :: offdiag(lx,ly,lz)


    n_with_ghst = (lx+2*g_nglayers)*(ly+2*g_nglayers)*(lz+2*g_nglayers)
    n           = lx*ly*lz

    Avj     = ZERO
    MinvAvj = ZERO
    wj      = ZERO
    vi      = ZERO
    vj      = ZERO

    lucky   = .false.
    nanflag = .false.
		  
    vj = kspvecs(:,:,:,j)

    call findAX(vj,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,&
		            bc_codes,bcvals,Avj(1:lx,1:ly,1:lz),diag,offdiag,&
                            lrank,rrank,brank,trank,krank,frank,&
                            llenx,lleny,llenz)

    call precond(MinvAvj,Avj,timederivflag,dt,lx,ly,lz,vel,&
                    dcoeff,reac,bc_codes,bcvals,&
                    lrank,rrank,brank,trank,krank,frank,&
                    llenx,lleny,llenz,1)

    do i=1,j
	vi = kspvecs(:,:,:,i)	
	call innerproduct(Hmat(i,j),MinvAvj(1:lx,1:ly,1:lz),vi(1:lx,1:ly,1:lz),n)
    enddo

    wj = MinvAvj

    do i=1,j
	vi = kspvecs(:,:,:,i)
	wj=wj - Hmat(i,j)*vi
    enddo

    call findnorm(Hmat(j+1,j),wj,n_with_ghst)

    if(Hmat(j+1,j) > NEARZERO) then
	kspvecs(:,:,:,j+1)=wj(:,:,:)/Hmat(j+1,j)

    else if(isit_nan(Hmat(j+1,j)) .eqv. .false.) then
		print *,"Hmat:",Hmat(j+1,j)
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
  subroutine performgmres(b,x0,x,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,&
				  bc_codes,bcvals,&
				  lrank,rrank,brank,trank,krank,frank,&
				  llenx,lleny,llenz,&
				  maxkspdim,nrestarts,findAX,&
				  precond,tol,success,printflag,initial_res,final_res)

	external :: findAX, precond

        logical, intent(in)   :: timederivflag
	integer,intent(in)    :: maxkspdim,nrestarts
        integer,intent(in)    :: lx,ly,lz
        real*8, intent(in)    :: dt,llenx,lleny,llenz
        integer,intent(in)    :: lrank,rrank,brank,trank,krank,frank
	character(LEN=4),intent(in)         :: bc_codes(NFACES)
	type(boundarycondition), intent(inout) :: bcvals

        
	real*8, intent(in)    :: b(lx,ly,lz)

        real*8, intent(inout) :: x0(-g_nglayers+1:lx+g_nglayers,&
                                    -g_nglayers+1:ly+g_nglayers,&
                                    -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout) :: x(-g_nglayers+1:lx+g_nglayers,&
                                   -g_nglayers+1:ly+g_nglayers,&
                                   -g_nglayers+1:lz+g_nglayers)

	real*8, intent(in)    :: vel(-g_nglayers+1:lx+g_nglayers,&
                                     -g_nglayers+1:ly+g_nglayers,&
                                     -g_nglayers+1:lz+g_nglayers,NDIM)

	real*8, intent(in)    :: dcoeff(-g_nglayers+1:lx+g_nglayers,&
                                     -g_nglayers+1:ly+g_nglayers,&
                                     -g_nglayers+1:lz+g_nglayers)

	real*8, intent(in)    :: reac(-g_nglayers+1:lx+g_nglayers,&
                                      -g_nglayers+1:ly+g_nglayers,&
                                      -g_nglayers+1:lz+g_nglayers)

        real*8, intent(in)    :: tol
        logical, intent(out)  :: success
        logical, intent(in)   :: printflag
        real*8, intent(out)   :: initial_res
        real*8, intent(out)   :: final_res

	integer :: i,j

	real*8                 :: r0(-g_nglayers+1:lx+g_nglayers,&
                                     -g_nglayers+1:ly+g_nglayers,&
                                     -g_nglayers+1:lz+g_nglayers)

	real*8                 :: Minvr(-g_nglayers+1:lx+g_nglayers,&
                                        -g_nglayers+1:ly+g_nglayers,&
                                        -g_nglayers+1:lz+g_nglayers)

	real*8                 :: Ax0(-g_nglayers+1:lx+g_nglayers,&
                                      -g_nglayers+1:ly+g_nglayers,&
                                      -g_nglayers+1:lz+g_nglayers)

	real*8                 :: Ax(-g_nglayers+1:lx+g_nglayers,&
                                      -g_nglayers+1:ly+g_nglayers,&
                                      -g_nglayers+1:lz+g_nglayers)

	real*8                :: v1(-g_nglayers+1:lx+g_nglayers,&
                                    -g_nglayers+1:ly+g_nglayers,&
                                    -g_nglayers+1:lz+g_nglayers)

	real*8                :: r(-g_nglayers+1:lx+g_nglayers,&
                                   -g_nglayers+1:ly+g_nglayers,&
                                   -g_nglayers+1:lz+g_nglayers)
	
        real*8 :: diag(lx,ly,lz)
	real*8 :: offdiag(lx,ly,lz)
	real*8 :: beta,residnorm
	real*8 :: y(maxkspdim)

	real*8 :: kspvectors(-g_nglayers+1:lx+g_nglayers,&
                             -g_nglayers+1:ly+g_nglayers,&
                             -g_nglayers+1:lz+g_nglayers,&
                             maxkspdim+1)

	real*8 :: Hmat(maxkspdim+1,maxkspdim)

	real*8 :: eps

	logical :: lucky
        logical :: nanflag
        integer :: luckykspdim
        integer :: kspdim,n,n_with_ghst

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

        n           = lx*ly*lz
        n_with_ghst = (lx+2*g_nglayers)*(ly+2*g_nglayers)*(lz+2*g_nglayers)

	call findAX(x0,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,&
		    bc_codes,bcvals,Ax0(1:lx,1:ly,1:lz),diag,offdiag,&
                    lrank,rrank,brank,trank,krank,frank,&
                    llenx,lleny,llenz)
        
	r0(1:lx,1:ly,1:lz) = b(1:lx,1:ly,1:lz)-Ax0(1:lx,1:ly,1:lz)

	call findnorm(initial_res,r0,n_with_ghst)

	call precond(Minvr,r0,timederivflag,dt,lx,ly,lz,vel,&
                    dcoeff,reac,bc_codes,bcvals,&
                    lrank,rrank,brank,trank,krank,frank,&
                    llenx,lleny,llenz,1)

	r = Minvr
	x  = x0

	lucky      = .false.
        nanflag    = .false.
        success    = .true.

	call findnorm(beta,r,n_with_ghst)
        i = 0
        if(printflag .eqv. .true.) print *,"restart iteration:",i,"residual norm:",beta

	do i=1,nrestarts
		  	
	    call findnorm(beta,r,n_with_ghst)
	    if(abs(beta) .lt. tol) exit
            x0 = x
	    v1 = r/beta
            kspvectors(:,:,:,1) = v1(:,:,:)

            do kspdim=1,maxkspdim
		
                x = x0
  		call  arnoldialgorithm(kspdim,maxkspdim,&
			      timederivflag,dt,lx,ly,lz,&
                              vel,dcoeff,reac,&
			      bc_codes,bcvals,&
                              lrank,rrank,brank,trank,krank,frank,&
                              llenx,lleny,llenz,Hmat,&
                              kspvectors,findAX,precond,&
			      lucky,luckykspdim,nanflag)

                if(nanflag .eqv. .true.) then
                    success=.false.
                    exit
                endif

                call leastsquaresminimize(y,Hmat,kspdim,maxkspdim,beta)

		do j=1,kspdim
			x=x+y(j)*kspvectors(:,:,:,j)
		enddo
		
	        call findAX(x,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,&
		            bc_codes,bcvals,Ax(1:lx,1:ly,1:lz),diag,offdiag,&
                            lrank,rrank,brank,trank,krank,frank,&
                            llenx,lleny,llenz)

		r(1:lx,1:ly,1:lz) = b(1:lx,1:ly,1:lz) - Ax(1:lx,1:ly,1:lz)

		call precond(Minvr,r,timederivflag,dt,lx,ly,lz,vel,&
                    dcoeff,reac,bc_codes,bcvals,&
                    lrank,rrank,brank,trank,krank,frank,&
                    llenx,lleny,llenz,1)

		r = Minvr
	        call findnorm(residnorm,r,n_with_ghst)

		if(lucky .eqv. .true.) then
                    if(printflag .eqv. .true.) print *,"lucky condition",residnorm
                    exit
		endif

		if((residnorm .le. eps) .and. (kspdim .gt. 1)) then
			exit
		endif

	    enddo

            if((residnorm .le. eps) .or. (lucky .eqv. .true.) .or. (success .eqv. .false.)) then
                if(printflag .eqv. .true.) print *,"restart iteration:",i,"residual norm:",residnorm
		final_res = residnorm
                exit
            endif
        
            if(printflag .eqv. .true.) print *,"restart iteration:",i,"residual norm:",residnorm

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

end module solvergmres_module
