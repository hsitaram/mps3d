module solvergmres_module

    use vectorfunctions
    implicit none

contains

    !===================================================
    subroutine arnoldialgorithm(j,maxkspdim,n,&
            timederivfactor,vel,dcoeff,&
            reac,dirc_bc_flags,flux_bc_flags,&
            dircvals,fluxvals,dx,dt,Hmat,&
            kspvecs,findAX,precond,&
            lucky,luckykspdim,nanflag,fluxscheme)

        integer,intent(in)    :: j,n,maxkspdim
        real*8, intent(inout) :: Hmat(maxkspdim+1,maxkspdim)
        real*8, intent(inout) :: kspvecs(n,maxkspdim+1)
        logical,intent(inout) :: lucky
        logical,intent(inout) :: nanflag
        integer,intent(inout) :: luckykspdim
        real*8, intent(in)    :: vel(n),dcoeff(n),reac(n)
        real*8, intent(in)    :: dx,dt
        logical,intent(in)    :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in)    :: dircvals(2),fluxvals(2)
        real*8, intent(in)    :: timederivfactor
        integer, intent(in)   :: fluxscheme

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

        call findAX(Avj,vj,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
            flux_bc_flags,dircvals,fluxvals,dx,dt,n,fluxscheme)

        call precond(MinvAvj,Avj,timederivfactor,vel,dcoeff,reac,&
            dirc_bc_flags,&
            flux_bc_flags,dircvals,fluxvals,&
            dx,dt,n,fluxscheme)

        !Gram Schmidt orthogonalization
        wj = MinvAvj
        do i=1,j
            vi = kspvecs(:,i)	
            call innerproduct(Hmat(i,j),MinvAvj,vi,n)
            wj=wj - Hmat(i,j)*vi
        enddo

        call findnorm(Hmat(j+1,j),wj,n)

        if(Hmat(j+1,j) > VERYSMALL) then
            kspvecs(:,j+1)=wj(:)/Hmat(j+1,j)

        else if(isit_nan(Hmat(j+1,j)) .eqv. .false.) then
            lucky=.true.
            luckykspdim=j;
        else
            print *,"NaN or Inf detected:",Hmat(j+1,j)
            nanflag=.true.
        endif

    end subroutine arnoldialgorithm
    !===================================================
    subroutine uppertrisolve(Hmat,maxkspdim,beta_e1,k,y)

        integer,intent(in)   :: k
        integer,intent(in)   :: maxkspdim
        real*8,intent(inout) :: Hmat(maxkspdim+1,maxkspdim)
        real*8,intent(out)   :: y(maxkspdim)
        real*8,intent(inout) :: beta_e1(maxkspdim)

        integer :: i,j

        do i=k,1,-1
            y(i)=beta_e1(i)
            do j=i+1,k 
                y(i) = y(i) - Hmat(i,j) * y(j)
            enddo
            y(i) = y(i) / Hmat(i,i)
        enddo

    end subroutine uppertrisolve
    !===================================================
    subroutine givens_rotate(r, s, cosine, sine)

        real*8, intent(in) :: r,s
        real*8, intent(out) :: cosine, sine

        if (r .eq. 0) then
            cosine = 0.d0
            sine = 1.d0
        else
            cosine = r/sqrt(r**2 + s**2)
            sine = s/sqrt(r**2 + s**2)
        endif

    end subroutine givens_rotate
    !===================================================
    subroutine triangularize(Hessmat,maxkspdim,cos_arr,sin_arr,beta_e1,k,error)

        integer,intent(in)   :: k
        integer,intent(in)   :: maxkspdim
        real*8,intent(inout) :: Hessmat(maxkspdim+1,maxkspdim)
        real*8,intent(inout) :: cos_arr(maxkspdim)
        real*8,intent(inout) :: sin_arr(maxkspdim)
        real*8,intent(inout) :: beta_e1(maxkspdim+1)
        real*8,intent(out)   :: error

        real*8 :: temp, cosine, sine
        integer :: i

        !apply for ith column
        do i=1,k-1
            temp = cos_arr(i) * Hessmat(i,k) + sin_arr(i) * Hessmat(i + 1,k);
            Hessmat(i+1,k) = -sin_arr(i) * Hessmat(i,k) + cos_arr(i) * Hessmat(i + 1,k);
            Hessmat(i,k) = temp
        enddo

        !update the next sin cos values for rotation
        call givens_rotate(Hessmat(k,k), Hessmat(k + 1,k), cosine, sine)

        !eliminate H(i + 1, i)
        Hessmat(k,k) = cosine * Hessmat(k,k) + sine * Hessmat(k+1,k)
        Hessmat(k+1,k) = 0.d0

        cos_arr(k) = cosine;
        sin_arr(k) = sine;

        beta_e1(k+1) = -sine * beta_e1(k);
        beta_e1(k)   = cosine * beta_e1(k);


        error = abs(beta_e1(k + 1))

    end subroutine triangularize
    !===================================================
    subroutine printvec(vec,n)

        integer, intent(in) :: n
        real*8, intent(in) :: vec(n)
        
        integer :: i
                
        do i=1,n
            write(*,'(F10.4)') vec(i)
        enddo

    end subroutine printvec
    !===================================================
    subroutine printmat(mat,m,n)
        
        integer, intent(in) :: m,n
        real*8, intent(in) :: mat(m,n)
        
        integer :: i,j

        do, i=1,m
            write(*,'(100g15.4)') ( mat(i,j), j=1,n )
        enddo


    end subroutine printmat
    !===================================================
    subroutine performgmres(b,x0,x,timederivfactor,vel,&
            dcoeff,reac,&
            dirc_bc_flags,flux_bc_flags,&
            dircvals,fluxvals,dx,dt,maxkspdim,&
            n,nrestarts,findAX,&
            precond,tol,success,printflag,initial_res,fluxscheme)

        external :: findAX, precond

        integer,intent(in)    :: maxkspdim,n,nrestarts
        real*8, intent(in)    :: b(n)
        real*8, intent(inout) :: x0(n)
        real*8, intent(inout) :: x(n)
        real*8, intent(in)    :: tol
        logical, intent(out)  :: success
        logical, intent(in)   :: printflag
        integer, intent(in)   :: fluxscheme

        real*8, intent(out)   :: initial_res
        real*8, intent(in)    :: vel(n),dcoeff(n),reac(n)
        real*8, intent(in)    :: dx,dt
        logical,intent(in)    :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in)    :: dircvals(2),fluxvals(2)
        real*8, intent(in)    :: timederivfactor

        integer :: i,j

        real*8 :: r0(n) !initial residual
        real*8 :: Minvr(n) !preconditioned residual
        real*8 :: Minvb(n) !preconditioned rhs
        real*8 :: Ax0(n),Ax(n) 
        real*8 :: r(n),v1(n) !residual and first ksp vector

        real*8 :: beta,residnorm
        real*8 :: y(maxkspdim)
        real*8 :: cos_arr(maxkspdim)
        real*8 :: sin_arr(maxkspdim)
        real*8 :: beta_e1(maxkspdim+1)

        real*8 :: kspvectors(n,maxkspdim+1)
        real*8 :: Hmat(maxkspdim+1,maxkspdim)

        real*8 :: eps
        real*8 :: b_norm

        logical :: lucky
        logical :: nanflag
        integer :: luckykspdim
        integer :: kspdim, optkspdim

        eps        = tol
        x = x0
        
        Ax0        = ZERO
        call findAX(Ax0,x0,timederivfactor,vel,dcoeff,reac,&
            dirc_bc_flags,flux_bc_flags,&
            dircvals,fluxvals,dx,dt,n,fluxscheme)

        !initial residual
        r0 = b-Ax0
        call findnorm(initial_res,r0,n)


        lucky      = .false.
        nanflag    = .false.
        success    = .true.

        do i=1,nrestarts

            Hmat       = ZERO
            kspvectors = ZERO
            r0         = ZERO
            Minvr      = ZERO
            Ax0        = ZERO
            Ax         = ZERO
            r          = ZERO
            v1         = ZERO
            y          = ZERO
            cos_arr    = ZERO
            sin_arr    = ZERO
            beta_e1    = ZERO
            Minvb      = ZERO

            call precond(Minvb,b,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
                flux_bc_flags,dircvals,fluxvals,dx,dt,n,fluxscheme)
            call findnorm(b_norm,Minvb,n)

            if(b_norm .eq. 0.d0) then
                print *,"rhs is zero"
                x=0.d0
                exit
            endif

            x0 = x

            !find Ax0
            call findAX(Ax0,x0,timederivfactor,vel,dcoeff,reac,&
                dirc_bc_flags,flux_bc_flags,&
                dircvals,fluxvals,dx,dt,n,fluxscheme)

            !initial residual
            r0 = b-Ax0

            !precondition residual
            call precond(Minvr,r0,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
                flux_bc_flags,dircvals,fluxvals,dx,dt,n,fluxscheme)

            !assign residual to be the preconditioned residual
            r  = Minvr

            !find norm of r and assign to beta
            call findnorm(beta,r,n)

            !first ksp vector
            v1 = r/beta
            kspvectors(:,1) = v1(:)
            beta_e1(1)      = beta

            if(printflag .eqv. .true.) print *,"restart iteration:",i,&
                "normalized residual norm:",beta/b_norm

            optkspdim=maxkspdim
            do kspdim=1,maxkspdim

                x = x0

                call arnoldialgorithm(kspdim,maxkspdim,n,&
                    timederivfactor,vel,dcoeff,&
                    reac,dirc_bc_flags,flux_bc_flags,&
                    dircvals,fluxvals,dx,dt,Hmat,kspvectors,findAX,precond,&
                    lucky,luckykspdim,nanflag,fluxscheme)

                !call printmat(kspvectors,n,maxkspdim+1)

                if(nanflag .eqv. .true.) then
                    call abort()
                    success=.false.
                    exit
                endif

                !print *,"before triangularize"
                !call printmat(Hmat,maxkspdim+1,maxkspdim)
                
                call triangularize(Hmat,maxkspdim,cos_arr,sin_arr,beta_e1,kspdim,residnorm)
                
                !print *,"after triangularize"
                !call printmat(Hmat,maxkspdim+1,maxkspdim)

                if(printflag .eqv. .true.) print *,"normalized residnorm, kspdim:",residnorm/b_norm,kspdim

                if(lucky .eqv. .true.) then
                    optkspdim=kspdim
                    if(printflag .eqv. .true.) print *,"lucky condition"
                    exit
                endif

                if(residnorm/b_norm .le. eps) then
                    optkspdim=kspdim
                    exit
                endif

            enddo
            call uppertrisolve(Hmat,maxkspdim,beta_e1,optkspdim,y)

            do j=1,optkspdim
                x=x+y(j)*kspvectors(:,j)
            enddo

            !print *,"y:",y

            if((residnorm/b_norm .le. eps) .or. (lucky .eqv. .true.) .or. (success .eqv. .false.)) then
                if(printflag .eqv. .true.) print *,"after restart iteration:",i,"normalized residual norm < tol:",&
                    residnorm/b_norm,tol
                exit
            endif

        enddo
        
        !Ax        = ZERO
        !call findAX(Ax,x,timederivfactor,vel,dcoeff,reac,&
        !    dirc_bc_flags,flux_bc_flags,&
        !    dircvals,fluxvals,dx,dt,n)

        !initial residual
        !r = b-Ax
        !call findnorm(final_res,r,n)


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
