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

        !Gram Schmidt orthogonalization
        wj = MinvAvj
        do i=1,j
            vi = kspvecs(:,:,:,i)	
            call innerproduct(Hmat(i,j),MinvAvj(1:lx,1:ly,1:lz),vi(1:lx,1:ly,1:lz),n)
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

        real*8                 :: Minvb(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8                 :: Ax0(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8                 :: Ax(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8                :: r(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8                :: v1(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8 :: diag(lx,ly,lz)
        real*8 :: offdiag(lx,ly,lz)
        real*8 :: beta,residnorm
        real*8 :: y(maxkspdim)
        real*8 :: cos_arr(maxkspdim)
        real*8 :: sin_arr(maxkspdim)
        real*8 :: beta_e1(maxkspdim+1)

        real*8 :: kspvectors(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers,&
            maxkspdim+1)

        real*8 :: Hmat(maxkspdim+1,maxkspdim)

        real*8 :: eps
        real*8 :: b_norm

        logical :: lucky
        logical :: nanflag
        integer :: luckykspdim
        integer :: kspdim,optkspdim
        integer :: n,n_with_ghst

        eps        = tol

        n           = lx*ly*lz
        n_with_ghst = (lx+2*g_nglayers)*(ly+2*g_nglayers)*(lz+2*g_nglayers)


        call findnorm(initial_res,r0,n_with_ghst)

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

            call precond(Minvb,b,timederivflag,dt,lx,ly,lz,vel,&
                dcoeff,reac,bc_codes,bcvals,&
                lrank,rrank,brank,trank,krank,frank,&
                llenx,lleny,llenz,1)

            call findnorm(b_norm,Minvb,n_with_ghst)

            if(b_norm .eq. 0) then
                print *,"rhs is zero"
                x=0.d0
                exit
            endif

            x0 = x

            !find AX0
            call findAX(x0,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,&
                bc_codes,bcvals,Ax0(1:lx,1:ly,1:lz),diag,offdiag,&
                lrank,rrank,brank,trank,krank,frank,&
                llenx,lleny,llenz)

            !initial residual
            r0(1:lx,1:ly,1:lz) = b(1:lx,1:ly,1:lz)-Ax0(1:lx,1:ly,1:lz)

            call precond(Minvr,r0,timederivflag,dt,lx,ly,lz,vel,&
                dcoeff,reac,bc_codes,bcvals,&
                lrank,rrank,brank,trank,krank,frank,&
                llenx,lleny,llenz,1)

            !assign residual to be the preconditioned residual
            r = Minvr

            !find norm of r and assign to beta
            call findnorm(beta,r,n_with_ghst)

            !first ksp vector
            v1 = r/beta
            kspvectors(:,:,:,1) = v1(:,:,:)
            beta_e1(1)      = beta



            call findnorm(beta,r,n_with_ghst)
            if(abs(beta) .lt. tol) exit
            x0 = x
            v1 = r/beta
            kspvectors(:,:,:,1) = v1(:,:,:)

            if(printflag .eqv. .true.) print *,"restart iteration:",i,&
                "normalized residual norm:",beta/b_norm

            optkspdim=maxkspdim
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

                call triangularize(Hmat,maxkspdim,cos_arr,sin_arr,beta_e1,kspdim,residnorm)

                if(lucky .eqv. .true.) then
                    optkspdim=kspdim
                    if(printflag .eqv. .true.) print *,"lucky condition"
                    exit
                endif

            enddo
            call uppertrisolve(Hmat,maxkspdim,beta_e1,optkspdim,y)

            do j=1,optkspdim
                x=x+y(j)*kspvectors(:,:,:,j)
            enddo

            if((residnorm/b_norm .le. eps) .or. (lucky .eqv. .true.) .or. (success .eqv. .false.)) then
                if(printflag .eqv. .true.) print *,"after restart iteration:",i,"normalized residual norm < tol:",&
                    residnorm/b_norm,tol
                exit
            endif



        enddo

        Ax=ZERO        
        call findAX(x,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,&
            bc_codes,bcvals,Ax(1:lx,1:ly,1:lz),diag,offdiag,&
            lrank,rrank,brank,trank,krank,frank,&
            llenx,lleny,llenz)

        r(1:lx,1:ly,1:lz) = b(1:lx,1:ly,1:lz) - Ax(1:lx,1:ly,1:lz)
        call findnorm(beta,r,n_with_ghst)

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
