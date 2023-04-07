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
            !print *,"hmat:",Hmat(j+1,j)
            lucky=.true.
            luckykspdim=j;
        else
            print *,"NaN or Inf detected:",Hmat(j+1,j)
            nanflag=.true.;
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
        real*8,intent(inout) :: beta_e1(maxkspdim)
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
        real*8 :: Minvb(n)
        real*8 :: Ax0(n),Ax(n)
        real*8 :: r(n),v1(n)
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
        cos_arr    = ZERO
        sin_arr    = ZERO
        beta_e1    = ZERO

        call precond(Minvb,b,n)
        call findnorm(b_norm,Minvb,n)

        x  = x0

        lucky      = .false.
        nanflag    = .false.
        success    = .true.


        do i=1,nrestarts

            x0 = x
            call findAX(Ax0,x0,n)
            r0 = b-Ax0
            call findnorm(initial_res,r0,n)
            call precond(Minvr,r0,n)
            
            r  = Minvr
            call findnorm(beta,r,n)
            
            v1 = r/beta
            kspvectors(:,1) = v1(:)
            beta_e1(1)      = beta
            cos_arr         = ZERO
            sin_arr         = ZERO
        
            if(printflag .eqv. .true.) print *,"restart iteration:",i,&
                "normalized residual norm:",beta/b_norm

            optkspdim=maxkspdim
            do kspdim=1,maxkspdim

                x = x0

                call arnoldialgorithm(kspdim,maxkspdim,n,&
                    Hmat,kspvectors,findAX,precond,&
                    lucky,luckykspdim,nanflag)
            
                if(nanflag .eqv. .true.) then
                    success=.false.
                    exit
                endif
                
                call triangularize(Hmat,maxkspdim,cos_arr,sin_arr,beta_e1,kspdim,residnorm)
                
                print *,"normalized residnorm, kspdim:",residnorm/b_norm,kspdim

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
