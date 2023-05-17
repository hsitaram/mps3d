module discretepde

    use globalvars
    use par_decompose

    implicit none

    !TODO: simplify boundary if loops
    !TODO: add soln/dt in find_b instead of application side

contains
    !=========================================================================================================
    subroutine find_redblack_indices(lx,ly,lz,redindices,blackindices)

        integer,intent(in)  :: lx,ly,lz	
        integer,allocatable :: redindices(:)
        integer,allocatable :: blackindices(:)

        integer :: i,j,k
        integer :: gl_i,gl_j,gl_k 
        integer :: bindex,rindex,onedindex
        integer :: redsize,blacksize

        redsize   = lx*ly*lz/2 + 1
        blacksize = lx*ly*lz/2 + 1

        !addition of 1 is for safety, if there are odd number
        !of cells, reds or blacks may be more than the other
        !we also dont which one will be more
        !it will be different on different processors.
        !It is corrected in the smoothing routine
        allocate(redindices(redsize))
        allocate(blackindices(blacksize))	

        redindices   = DUMMYVAL
        blackindices = DUMMYVAL

        bindex=1
        rindex=1
        do k=1,lz
            do j=1,ly
                do i=1,lx

                    call local_to_global_index(i,j,k,&
                        g_myproc,gl_i,gl_j,gl_k)

                    call find1dindex(i,j,k,lx,ly,lz,onedindex)

                    if(mod((gl_i+gl_j+gl_k),2) .eq. 0) then

                        blackindices(bindex) = onedindex
                        bindex=bindex+1

                    else
                        redindices(rindex) = onedindex
                        rindex=rindex+1
                    endif
                enddo
            enddo
        enddo


    end subroutine find_redblack_indices
    !==================================================================================================================
    subroutine perform_gjacobi_smoothing(solnvar,timederivflag,dt,lx,ly,lz,applybcflag,mgridsource,res,&
            vel,dcoeff,reac,source,&
            bc_codes,bcvals,&
            lrank,rrank,brank,trank,krank,frank,&
            llenx,lleny,llenz,niter)

        real*8,  intent(in)    :: dt
        integer, intent(in)    :: lx,ly,lz
        integer, intent(in)    :: niter

        real*8, intent(inout)  :: solnvar(-g_nglayers+1:lx+g_nglayers, &
            -g_nglayers+1:ly+g_nglayers, &
            -g_nglayers+1:lz+g_nglayers)

        integer, intent(in)   :: lrank,rrank,brank,trank,krank,frank
        real*8, intent(inout)    :: res(lx,ly,lz)

        real*8,intent(in)     :: llenx,lleny,llenz
        real*8,intent(in)     :: mgridsource(lx,ly,lz)
        logical,intent(in)    :: timederivflag
        logical,intent(in)    :: applybcflag

        character(LEN=4),intent(in)         :: bc_codes(NFACES)
        type(boundarycondition), intent(inout) :: bcvals

        real*8, intent(inout)    :: vel(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers,NDIM)

        real*8, intent(inout)    :: dcoeff(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout)    :: reac(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout)    :: source(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8   :: AX(lx,ly,lz)
        real*8   :: diag(lx,ly,lz)
        real*8   :: offdiag(lx,ly,lz)
        real*8   :: b(lx,ly,lz)

        integer :: it

        call find_bvec(solnvar,dt,lx,ly,lz,mgridsource,vel,dcoeff,source,bc_codes,bcvals,&
            applybcflag,b,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)

        do it=1,niter

            call find_AX(solnvar,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,bc_codes,bcvals,&
                AX,diag,offdiag,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)


            solnvar(1:lx,1:ly,1:lz) = (b - offdiag)/diag

        enddo

        call find_AX(solnvar,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,bc_codes,bcvals,&
            AX,diag,offdiag,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)
        res = b - AX

    end subroutine perform_gjacobi_smoothing
    !============================================================================================================
    subroutine perform_gseidel_smoothing(solnvar,timederivflag,dt,lx,ly,lz,applybcflag,mgridsource,res,&
            vel,dcoeff,reac,source,&
            bc_codes,bcvals,&
            lrank,rrank,brank,trank,krank,frank,&
            llenx,lleny,llenz,niter)

        real*8, intent(in)    :: dt
        integer, intent(in)   :: lx,ly,lz
        integer, intent(in)   :: niter

        integer, intent(in)   :: lrank,rrank,brank,trank,krank,frank
        real*8, intent(in)    :: llenx,lleny,llenz

        real*8, intent(inout) :: solnvar(-g_nglayers+1:lx+g_nglayers, &
            -g_nglayers+1:ly+g_nglayers, &
            -g_nglayers+1:lz+g_nglayers)

        logical, intent(in)   :: timederivflag
        logical, intent(in)   :: applybcflag

        real*8, intent(inout) :: res(lx,ly,lz)
        real*8,intent(in)     :: mgridsource(lx,ly,lz)

        real*8, intent(inout)    :: vel(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers,NDIM)

        real*8, intent(inout)    :: dcoeff(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout)    :: reac(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout)    :: source(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)


        character(LEN=4),intent(in)         :: bc_codes(NFACES)
        type(boundarycondition), intent(inout) :: bcvals

        integer :: ind,onedindex
        integer,allocatable :: redindices(:),blackindices(:)
        integer :: i,j,k,it
        integer :: redsize,blacksize

        real*8   :: AX(lx,ly,lz)
        real*8   :: diag(lx,ly,lz)
        real*8   :: offdiag(lx,ly,lz)
        real*8   :: b(lx,ly,lz)

        call find_redblack_indices(lx,ly,lz,redindices,blackindices)

        redsize=size(redindices)
        if(redindices(redsize) .eq. DUMMYVAL) then
            redsize=redsize-1
        endif

        blacksize=size(blackindices)
        if(blackindices(blacksize) .eq. DUMMYVAL) then
            blacksize=blacksize-1
        endif

        call find_bvec(solnvar,dt,lx,ly,lz,mgridsource,vel,dcoeff,source,bc_codes,bcvals,&
            applybcflag,b,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)

        !print *,g_myproc,"found bvec"

        do it=1,niter

            call find_AX(solnvar,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,bc_codes,bcvals,&
                AX,diag,offdiag,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)

            do ind=1,redsize

                onedindex=redindices(ind)
                call find3dindex(onedindex,lx,ly,lz,i,j,k)

                solnvar(i,j,k) = (b(i,j,k)-offdiag(i,j,k))/diag(i,j,k)
            enddo



            call find_AX(solnvar,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,bc_codes,bcvals,&
                AX,diag,offdiag,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)

            !update black indices
            do ind=1,blacksize

                onedindex=blackindices(ind)
                call find3dindex(onedindex,lx,ly,lz,i,j,k)

                solnvar(i,j,k) = (b(i,j,k)-offdiag(i,j,k))/diag(i,j,k)
            enddo

        enddo

        call find_AX(solnvar,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,bc_codes,bcvals,&
            AX,diag,offdiag,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)
        res = b - AX


    end subroutine perform_gseidel_smoothing
    !============================================================================================================
    subroutine gseidel_precond(MinvX,X,timederivflag,dt,lx,ly,lz,&
            vel,dcoeff,reac,bc_codes,bcvals,&
            lrank,rrank,brank,trank,krank,frank,&
            llenx,lleny,llenz,niter)

        real*8, intent(in)    :: dt
        integer, intent(in)   :: lx,ly,lz
        integer, intent(in)   :: niter

        integer, intent(in)   :: lrank,rrank,brank,trank,krank,frank
        real*8, intent(in)    :: llenx,lleny,llenz

        real*8, intent(inout) :: X(-g_nglayers+1:lx+g_nglayers, &
            -g_nglayers+1:ly+g_nglayers, &
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout) :: MinvX(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        logical, intent(in)   :: timederivflag

        real*8, intent(inout)    :: vel(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers,NDIM)

        real*8, intent(inout)    :: dcoeff(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout)    :: reac(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)


        character(LEN=4),intent(in)         :: bc_codes(NFACES)
        type(boundarycondition), intent(inout) :: bcvals

        integer :: ind,onedindex
        integer,allocatable :: redindices(:),blackindices(:)
        integer :: i,j,k,it
        integer :: redsize,blacksize

        real*8   :: AX(lx,ly,lz)
        real*8   :: diag(lx,ly,lz)
        real*8   :: offdiag(lx,ly,lz)
        real*8   :: b(lx,ly,lz)

        call find_redblack_indices(lx,ly,lz,redindices,blackindices)

        redsize=size(redindices)
        if(redindices(redsize) .eq. DUMMYVAL) then
            redsize=redsize-1
        endif

        blacksize=size(blackindices)
        if(blackindices(blacksize) .eq. DUMMYVAL) then
            blacksize=blacksize-1
        endif

        MinvX = ZERO		
        
        if(lx*ly*lz .eq. 1) then

            call find_AX(MinvX,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,bc_codes,bcvals,&
                AX,diag,offdiag,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)

            MinvX(1:lx,1:ly,1:lz)=X(1:lx,1:ly,1:lz)/diag(1:lx,1:ly,1:lz)
        
        else
            do it=1,niter

                call find_AX(MinvX,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,bc_codes,bcvals,&
                    AX,diag,offdiag,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)

                do ind=1,redsize

                    onedindex=redindices(ind)
                    call find3dindex(onedindex,lx,ly,lz,i,j,k)

                    MinvX(i,j,k) = (X(i,j,k)-offdiag(i,j,k))/diag(i,j,k)
                enddo


                call find_AX(MinvX,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,bc_codes,bcvals,&
                    AX,diag,offdiag,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)

                !update black indices
                do ind=1,blacksize

                    onedindex=blackindices(ind)
                    call find3dindex(onedindex,lx,ly,lz,i,j,k)

                    MinvX(i,j,k) = (X(i,j,k)-offdiag(i,j,k))/diag(i,j,k)
                enddo

            enddo
        endif

    end subroutine gseidel_precond
    !============================================================================================================
    subroutine no_precond(MinvX,X,timederivflag,dt,lx,ly,lz,&
            vel,dcoeff,reac,bc_codes,bcvals,&
            lrank,rrank,brank,trank,krank,frank,&
            llenx,lleny,llenz,niter)

        real*8, intent(in)    :: dt
        integer, intent(in)   :: lx,ly,lz
        integer, intent(in)   :: niter

        integer, intent(in)   :: lrank,rrank,brank,trank,krank,frank
        real*8, intent(in)    :: llenx,lleny,llenz

        real*8, intent(inout) :: X(-g_nglayers+1:lx+g_nglayers, &
            -g_nglayers+1:ly+g_nglayers, &
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout) :: MinvX(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        logical, intent(in)   :: timederivflag

        real*8, intent(inout)    :: vel(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers,NDIM)

        real*8, intent(inout)    :: dcoeff(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout)    :: reac(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)


        character(LEN=4),intent(in)         :: bc_codes(NFACES)
        type(boundarycondition), intent(inout) :: bcvals


        MinvX = X

    end subroutine no_precond
    !============================================================================================================
    subroutine find_bvec(solnvar,dt,lx,ly,lz,mgridsource,vel,dcoeff,source,bc_codes,bcvals,&
            applybcflag,b,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)

        real*8, intent(in)    :: dt
        integer, intent(in)   :: lx,ly,lz
        real*8, intent(in)    :: solnvar(-g_nglayers+1:lx+g_nglayers, &
            -g_nglayers+1:ly+g_nglayers, &
            -g_nglayers+1:lz+g_nglayers)

        integer, intent(in)   :: lrank,rrank,brank,trank,krank,frank

        real*8, intent(in)    :: llenx,lleny,llenz

        real*8,intent(in)     :: mgridsource(lx,ly,lz)

        real*8, intent(inout)    :: vel(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers,NDIM)

        real*8, intent(inout)    :: dcoeff(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout)    :: source(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        logical,intent(in)    :: applybcflag

        character(LEN=4), intent(in) :: bc_codes(NFACES)
        type(boundarycondition), intent(inout) :: bcvals

        real*8,intent(out)    :: b(lx,ly,lz)

        integer :: i,j,k
        real*8 :: dx,dy,dz
        real*8 :: vel_b

        dx = llenx/lx
        dy = lleny/ly
        dz = llenz/lz

        b=0.d0

        if(applybcflag .eqv. .true.) then

            if(lrank .lt. 0) then

                i = 1
                if(bc_codes(LEFT) .eq. 'DIRC') then

                    do k=1,lz
                        do j=1,ly

                            vel_b = vel(i,j,k,XDIR)

                            b(i,j,k)    = b(i,j,k) + TWO * dcoeff(i,j,k) * bcvals%left(j,k)/(dx**2)

                            !upwind
                            if(vel_b .gt. ZERO) then
                                b(i,j,k) = b(i,j,k) + vel_b * bcvals%left(j,k)/dx
                            endif	

                        enddo
                    enddo

                else if(bc_codes(LEFT) .eq. "FLUX") then

                    do k=1,lz
                        do j=1,ly
                            b(i,j,k) = b(i,j,k) - bcvals%left(j,k)/dx
                        enddo
                    enddo
                endif
            endif

            if(rrank .lt. 0) then

                i = lx
                if(bc_codes(RIGHT) .eq. 'DIRC') then

                    do k=1,lz
                        do j=1,ly

                            vel_b = vel(i,j,k,XDIR)

                            b(i,j,k)    = b(i,j,k) + TWO * dcoeff(i,j,k) * bcvals%right(j,k)/(dx**2)

                            !upwind
                            if(vel_b .lt. ZERO) then
                                b(i,j,k) = b(i,j,k) - vel_b * bcvals%right(j,k)/dx
                            endif	

                        enddo
                    enddo

                else if(bc_codes(RIGHT) .eq. "FLUX") then

                    do k=1,lz
                        do j=1,ly
                            b(i,j,k) = b(i,j,k) - bcvals%right(j,k)/dx
                        enddo
                    enddo
                endif
            endif

            if(brank .lt. 0) then

                j = 1
                if(bc_codes(BOTTOM) .eq. 'DIRC') then

                    do k=1,lz
                        do i=1,lx

                            vel_b = vel(i,j,k,YDIR)

                            b(i,j,k)    = b(i,j,k) + TWO * dcoeff(i,j,k) * bcvals%bottom(i,k)/(dy**2)

                            !upwind
                            if(vel_b .gt. ZERO) then
                                b(i,j,k) = b(i,j,k) + vel_b * bcvals%bottom(i,k)/dy
                            endif	

                        enddo
                    enddo

                else if(bc_codes(BOTTOM) .eq. "FLUX") then

                    do k=1,lz
                        do i=1,lx
                            b(i,j,k) = b(i,j,k) - bcvals%bottom(i,k)/dy
                        enddo
                    enddo
                endif
            endif

            if(trank .lt. 0) then

                j = ly
                if(bc_codes(TOP) .eq. 'DIRC') then

                    do k=1,lz
                        do i=1,lx

                            vel_b = vel(i,j,k,YDIR)

                            b(i,j,k)    = b(i,j,k) + TWO * dcoeff(i,j,k) * bcvals%top(i,k)/(dy**2)

                            !upwind
                            if(vel_b .lt. ZERO) then
                                b(i,j,k) = b(i,j,k) - vel_b * bcvals%top(i,k)/dy
                            endif	

                        enddo
                    enddo

                else if(bc_codes(TOP) .eq. "FLUX") then

                    do k=1,lz
                        do i=1,lx
                            b(i,j,k) = b(i,j,k) - bcvals%top(i,k)/dy
                        enddo
                    enddo
                endif
            endif

            if(krank .lt. 0) then

                k = 1
                if(bc_codes(BACK) .eq. 'DIRC') then

                    do j=1,ly
                        do i=1,lx

                            vel_b = vel(i,j,k,ZDIR)

                            b(i,j,k)    = b(i,j,k) + TWO * dcoeff(i,j,k) * bcvals%back(i,j)/(dz**2)

                            !upwind
                            if(vel_b .gt. ZERO) then
                                b(i,j,k) = b(i,j,k) + vel_b * bcvals%back(i,j)/dz
                            endif	

                        enddo
                    enddo

                else if(bc_codes(BACK) .eq. "FLUX") then

                    do j=1,ly
                        do i=1,lx
                            b(i,j,k) = b(i,j,k) - bcvals%back(i,j)/dz
                        enddo
                    enddo
                endif
            endif

            if(frank .lt. 0) then

                k = lz
                if(bc_codes(FRONT) .eq. 'DIRC') then

                    do j=1,ly
                        do i=1,lx

                            vel_b = vel(i,j,k,ZDIR)

                            b(i,j,k)    = b(i,j,k) + TWO * dcoeff(i,j,k) * bcvals%front(i,j)/(dz**2)

                            !upwind
                            if(vel_b .lt. ZERO) then
                                b(i,j,k) = b(i,j,k) - vel_b * bcvals%front(i,j)/dz
                            endif	

                        enddo
                    enddo

                else if(bc_codes(FRONT) .eq. "FLUX") then

                    do j=1,ly
                        do i=1,lx
                            b(i,j,k) = b(i,j,k) - bcvals%front(i,j)/dz
                        enddo
                    enddo
                endif
            endif

        endif

        !add source term
        b = b + mgridsource

        !only for step 1 of v cycle this is true
        if(applybcflag .eqv. .true.) then
            b = b + source(1:lx,1:ly,1:lz)
        endif

    end subroutine find_bvec
    !============================================================================================================
    subroutine findcellcontributions(soln_l,soln_r,k_l,k_r,vel_l,vel_r,&
            a_by_v,d_lr,diag_l,diag_r,offdiag_l,offdiag_r)

        !Note: left to right is always in positive x,y or z directions		
        !vel_l and vel_r are components of velocities along x,y or z
        !directions

        !signs are assumed to be
        !du/dt + del.(vel u) = del.(k del u) + reac*u + source

        real*8, intent(in) :: soln_l,soln_r
        real*8, intent(in) :: k_l,k_r
        real*8, intent(in) :: vel_l,vel_r
        real*8, intent(in) :: a_by_v,d_lr

        real*8, intent(out) :: diag_l,diag_r
        real*8, intent(out) :: offdiag_l,offdiag_r

        real*8 :: k_half,vel_half

        diag_l=0.d0
        diag_r=0.d0

        offdiag_l=0.d0
        offdiag_r=0.d0

        !diffusion contribution
        k_half = HALF*(k_l + k_r)
        diag_l = diag_l + (k_half/d_lr)*a_by_v
        diag_r = diag_r + (k_half/d_lr)*a_by_v

        offdiag_l = offdiag_l - soln_r*(k_half/d_lr)*a_by_v
        offdiag_r = offdiag_r - soln_l*(k_half/d_lr)*a_by_v


        !convection contribution (UPWIND)
        ! Note the sign convention
        ! Flux is positive if going out of cell
        ! Flux is negative if going into cell
        ! Therefore, the signs are decided based on sign of
        ! vel_half and whether the flux is in or out
        vel_half = HALF*(vel_l + vel_r)
        if(vel_half .gt. ZERO) then
            diag_l    = diag_l    + vel_half*a_by_v
            offdiag_r = offdiag_r - vel_half*soln_l*a_by_v 
        else
            diag_r    = diag_r    - vel_half*a_by_v
            offdiag_l = offdiag_l + vel_half*soln_r*a_by_v
        endif

    end subroutine findcellcontributions
    !============================================================================================================
    subroutine find_AX(solnvar,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,bc_codes,bcvals,&
            AX,diag,offdiag,lrank,rrank,brank,trank,krank,frank,llenx,lleny,llenz)

        real*8, intent(in)    :: dt
        integer, intent(in)   :: lx,ly,lz
        type(boundarycondition), intent(inout) :: bcvals
        real*8, intent(inout) :: solnvar(-g_nglayers+1:lx+g_nglayers, &
            -g_nglayers+1:ly+g_nglayers, &
            -g_nglayers+1:lz+g_nglayers)

        logical, intent(in)   :: timederivflag

        character(LEN=4),intent(in) :: bc_codes(NFACES)

        real*8, intent(inout)    :: vel(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers,NDIM)

        real*8, intent(inout)    :: dcoeff(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        real*8, intent(inout)    :: reac(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,&
            -g_nglayers+1:lz+g_nglayers)

        integer, intent(in)   :: lrank,rrank,brank,trank,krank,frank
        real*8,intent(in)     :: llenx,lleny,llenz

        real*8,intent(out)    :: AX(lx,ly,lz)
        real*8,intent(out)    :: diag(lx,ly,lz)
        real*8,intent(out)    :: offdiag(lx,ly,lz)

        real*8 :: idmat(lx,ly,lz)

        integer :: imin,imax
        integer :: jmin,jmax
        integer :: kmin,kmax

        integer :: i,j,k

        integer :: lcell,rcell
        real*8 :: dx,dy,dz

        real*8 :: diag_l,diag_r
        real*8 :: offdiag_l,offdiag_r
        real*8 :: soln_l,soln_r

        real*8 :: dcoeff_l,dcoeff_r
        real*8 :: vel_l,vel_r
        real*8 :: dcoeff_b,vel_b

        dx = llenx/lx
        dy = lleny/ly
        dz = llenz/lz

        AX   = ZERO
        diag = ZERO
        offdiag = ZERO
        idmat = 1.d0

        !finite volume form without t derivative
        ! \sum( v phi - D \grad phi ) . n A/V  - R = S

        !exchange halo data
        if((lrank+rrank+brank+trank+krank+frank) .gt. 6*MPI_PROC_NULL) then

            call exchangehalodata(solnvar,lx,ly,lz)

            call exchangehalodata(vel(:,:,:,XDIR),lx,ly,lz)
            call exchangehalodata(vel(:,:,:,YDIR),lx,ly,lz)
            call exchangehalodata(vel(:,:,:,ZDIR),lx,ly,lz)

            call exchangehalodata(dcoeff,lx,ly,lz)
            call exchangehalodata(reac,lx,ly,lz)
        endif

        !Note that ghost values at physical boundaries are zero
        !============================================================================================
        !X direction
        !============================================================================================
        !loop over all interior x faces
        imin = 0
        imax = lx
        !only set min and max for physical boundaries
        !processor boundaries are still interior
        if(lrank .lt. 0) imin = 1
        if(rrank .lt. 0) imax = lx-1

        do k=1,lz
            do j=1,ly
                do i=imin,imax

                    !flux for simple laplace solve
                    lcell = i
                    rcell = i+1

                    soln_l = solnvar(lcell,j,k)
                    soln_r = solnvar(rcell,j,k)

                    dcoeff_l = dcoeff(lcell,j,k)
                    dcoeff_r = dcoeff(rcell,j,k)

                    vel_l = vel(lcell,j,k,XDIR)
                    vel_r = vel(rcell,j,k,XDIR)

                    call findcellcontributions(soln_l,soln_r,dcoeff_l,dcoeff_r,vel_l,vel_r,1/dx,dx,&
                        diag_l,diag_r,offdiag_l,offdiag_r)

                    if(lcell .gt. 0)  then

                        diag(lcell,j,k)    =    diag(lcell,j,k)  + diag_l
                        offdiag(lcell,j,k) = offdiag(lcell,j,k)  + offdiag_l

                    endif

                    if(rcell .le. lx) then

                        diag(rcell,j,k)    =    diag(rcell,j,k) + diag_r
                        offdiag(rcell,j,k) = offdiag(rcell,j,k) + offdiag_r

                    endif
                enddo
            enddo
        enddo

        if(lrank .lt. 0) then

            i = 1
            if(bc_codes(LEFT) .eq. 'DIRC') then

                do k=1,lz
                    do j=1,ly

                        vel_b    = vel(i,j,k,XDIR)
                        dcoeff_b = dcoeff(i,j,k)
                        !dcoeff_b = 1.d0

                        diag(i,j,k) =    diag(i,j,k) + TWO * dcoeff_b/(dx**2)

                        !upwind
                        if(vel_b .lt. ZERO) then
                            diag(i,j,k) = diag(i,j,k) - vel_b/dx
                        endif	
                    enddo
                enddo

            else if(bc_codes(LEFT) .eq. 'ZGRD') then

                do k=1,lz
                    do j=1,ly

                        vel_b    = vel(i,j,k,XDIR)
                        !no diffusive flux, only convective	
                        diag(i,j,k) = diag(i,j,k) - vel_b/dx
                    enddo
                enddo
            endif

        endif

        if(rrank .lt. 0) then

            i = lx
            if(bc_codes(RIGHT) .eq. 'DIRC') then

                do k=1,lz
                    do j=1,ly

                        vel_b    = vel(i,j,k,XDIR)
                        dcoeff_b = dcoeff(i,j,k)
                        !dcoeff_b = 1.d0

                        diag(i,j,k) =    diag(i,j,k) + TWO * dcoeff_b/(dx**2)

                        !upwind
                        if(vel_b .gt. ZERO) then
                            diag(i,j,k) = diag(i,j,k) + vel_b/dx
                        endif	
                    enddo
                enddo

            else if(bc_codes(RIGHT) .eq. 'ZGRD') then

                do k=1,lz
                    do j=1,ly

                        vel_b    = vel(i,j,k,XDIR)
                        !no diffusive flux, only convective	
                        diag(i,j,k) = diag(i,j,k) + vel_b/dx
                    enddo
                enddo
            endif
        endif
        !============================================================================================

        !============================================================================================
        !Y direction
        !============================================================================================
        jmin = 0
        jmax = ly
        !only set min and max for physical boundaries
        !processor boundaries are still interior
        if(brank .lt. 0) jmin = 1
        if(trank .lt. 0) jmax = ly-1

        do k=1,lz

            do j=jmin,jmax

                !flux for simple laplace solve
                lcell = j
                rcell = j+1

                do i=1,lx

                    soln_l = solnvar(i,lcell,k)
                    soln_r = solnvar(i,rcell,k)

                    dcoeff_l = dcoeff(i,lcell,k)
                    dcoeff_r = dcoeff(i,rcell,k)

                    vel_l = vel(i,lcell,k,YDIR)
                    vel_r = vel(i,rcell,k,YDIR)

                    call findcellcontributions(soln_l,soln_r,dcoeff_l,dcoeff_r,vel_l,vel_r,1/dy,dy,&
                        diag_l,diag_r,offdiag_l,offdiag_r)


                    if(lcell .gt. 0)  then

                        diag(i,lcell,k) =    diag(i,lcell,k) + diag_l
                        offdiag(i,lcell,k) = offdiag(i,lcell,k) + offdiag_l

                    endif

                    if(rcell .le. ly) then

                        diag(i,rcell,k)    = diag(i,rcell,k) + diag_r
                        offdiag(i,rcell,k) = offdiag(i,rcell,k) + offdiag_r

                    endif
                enddo
            enddo

        enddo

        if(brank .lt. 0) then

            j = 1
            if(bc_codes(BOTTOM) .eq. 'DIRC') then

                do k=1,lz
                    do i=1,lx

                        vel_b    = vel(i,j,k,YDIR)
                        dcoeff_b = dcoeff(i,j,k)
                        !dcoeff_b = 1.d0

                        diag(i,j,k) = diag(i,j,k) + TWO * dcoeff_b/(dy**2)

                        !upwind
                        if(vel_b .lt. ZERO) then
                            diag(i,j,k) = diag(i,j,k) - vel_b/dy
                        endif	
                    enddo
                enddo

            else if(bc_codes(BOTTOM) .eq. 'ZGRD') then

                do k=1,lz
                    do i=1,lx

                        vel_b    = vel(i,j,k,YDIR)
                        !no diffusive flux, only convective	
                        diag(i,j,k) = diag(i,j,k) - vel_b/dy
                    enddo
                enddo
            endif
        endif

        if(trank .lt. 0) then

            j = ly
            if(bc_codes(TOP) .eq. 'DIRC') then

                do k=1,lz
                    do i=1,lx

                        vel_b    = vel(i,j,k,YDIR)
                        dcoeff_b = dcoeff(i,j,k)
                        !dcoeff_b = 1.d0

                        diag(i,j,k) = diag(i,j,k) + TWO * dcoeff_b/(dy**2)

                        !upwind
                        if(vel_b .gt. ZERO) then
                            diag(i,j,k) = diag(i,j,k) + vel_b/dy
                        endif	
                    enddo
                enddo

            else if(bc_codes(TOP) .eq. 'ZGRD') then

                do k=1,lz
                    do i=1,lx

                        vel_b    = vel(i,j,k,YDIR)
                        !no diffusive flux, only convective	
                        diag(i,j,k) = diag(i,j,k) + vel_b/dy
                    enddo
                enddo
            endif
        endif
        !============================================================================================
        !print *,"g_myproc after y:",g_myproc,diag

        !============================================================================================
        !Z direction
        !============================================================================================
        kmin = 0
        kmax = lz
        !only set min and max for physical boundaries
        !processor boundaries are still interior
        if(krank .lt. 0) kmin = 1
        if(frank .lt. 0) kmax = lz-1

        do k=kmin,kmax

            !flux for simple laplace solve
            lcell = k
            rcell = k+1

            do j=1,ly
                do i=1,lx
                    soln_l = solnvar(i,j,lcell)
                    soln_r = solnvar(i,j,rcell)

                    dcoeff_l = dcoeff(i,j,lcell)
                    dcoeff_r = dcoeff(i,j,rcell)

                    vel_l = vel(i,j,lcell,ZDIR)
                    vel_r = vel(i,j,rcell,ZDIR)

                    call findcellcontributions(soln_l,soln_r,dcoeff_l,dcoeff_r,vel_l,vel_r,1/dz,dz,&
                        diag_l,diag_r,offdiag_l,offdiag_r)

                    if(lcell .gt. 0)  then

                        diag(i,j,lcell) =    diag(i,j,lcell) + diag_l
                        offdiag(i,j,lcell) = offdiag(i,j,lcell) + offdiag_l

                    endif

                    if(rcell .le. lz) then

                        diag(i,j,rcell)    = diag(i,j,rcell) + diag_r
                        offdiag(i,j,rcell) = offdiag(i,j,rcell) + offdiag_r

                    endif
                enddo
            enddo
        enddo


        if(krank .lt. 0) then

            k = 1
            if(bc_codes(BACK) .eq. 'DIRC') then

                do j=1,ly
                    do i=1,lx

                        vel_b    = vel(i,j,k,ZDIR)
                        dcoeff_b = dcoeff(i,j,k)
                        !dcoeff_b = 1.d0

                        diag(i,j,k) = diag(i,j,k) + TWO*dcoeff_b/(dz**2)

                        !upwind
                        if(vel_b .lt. ZERO) then
                            diag(i,j,k) = diag(i,j,k) - vel_b/dz
                        endif	
                    enddo
                enddo

            else if(bc_codes(BACK) .eq. 'ZGRD') then

                do j=1,ly
                    do i=1,lx
                        vel_b    = vel(i,j,k,ZDIR)
                        !no diffusive flux, only convective	
                        diag(i,j,k) = diag(i,j,k) - vel_b/dz
                    enddo
                enddo
            endif
        endif

        if(frank .lt. 0) then

            k = lz
            if(bc_codes(FRONT) .eq. 'DIRC') then

                do j=1,ly
                    do i=1,lx
                        vel_b    = vel(i,j,k,ZDIR)
                        dcoeff_b = dcoeff(i,j,k)
                        !dcoeff_b = 1.d0

                        diag(i,j,k) = diag(i,j,k) + TWO*dcoeff_b/(dz**2)

                        !upwind
                        if(vel_b .gt. ZERO) then
                            diag(i,j,k) = diag(i,j,k) + vel_b/dz
                        endif	
                    enddo
                enddo

            else if(bc_codes(FRONT) .eq. 'ZGRD') then

                do j=1,ly
                    do i=1,lx
                        vel_b    = vel(i,j,k,ZDIR)
                        !no diffusive flux, only convective	
                        diag(i,j,k) = diag(i,j,k) + vel_b/dz
                    enddo
                enddo
            endif
        endif

        !reaction contribution
        diag = diag - reac(1:lx,1:ly,1:lz)

        !time derivative contribution
        if(timederivflag .eqv. .true.) diag = diag + idmat/dt

        !============================================================================================
        !now find AX and update b with source term
        do k=1,lz
            do j=1,ly
                do i=1,lx
                    AX(i,j,k) = diag(i,j,k)*solnvar(i,j,k)+offdiag(i,j,k)
                enddo
            enddo
        enddo
        !============================================================================================

    end subroutine find_AX
    !============================================================================================================
    subroutine update_boundary_ghostvalues(solnvar,bc_codes,bcvals,lx,ly,lz)

        character(LEN=4),intent(in) :: bc_codes(NFACES)
        type(boundarycondition), intent(in) :: bcvals

        integer, intent(in)   :: lx,ly,lz
        real*8, intent(inout) :: solnvar(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,-g_nglayers+1:lz+g_nglayers)

        integer :: i,j,k
        real*8 :: x,y,z

        !for dirichlet boundaries, a second order value is assigned to 
        !ghost cells.
        !  --g--|(d)---c---|
        ! phi_d = 0.5*(phi_c + phi_g)

        !g_nglayers is 1
        !for higher order cases we have to do something else

        call exchangehalodata(solnvar,lx,ly,lz)


        !==================Left side=================================
        if(g_lrank .lt. 0) then

            if(bc_codes(LEFT) .eq. 'DIRC') then

                do k=1,lz
                    do j=1,ly
                        solnvar(0,j,k) = TWO*bcvals%left(j,k) - solnvar(1,j,k)
                    enddo
                enddo
            else
                do k=1,lz
                    do j=1,ly
                        solnvar(0,j,k) = solnvar(1,j,k)
                    enddo
                enddo
            endif
        endif
        !=============================================================


        !==================Right side=================================
        if(g_rrank .lt. 0) then

            if(bc_codes(RIGHT) .eq. 'DIRC') then
                do k=1,lz
                    do j=1,ly
                        solnvar(lx+g_nglayers,j,k) = TWO*bcvals%right(j,k) - solnvar(lx,j,k)
                    enddo
                enddo
            else
                do k=1,lz
                    do j=1,ly
                        solnvar(lx+g_nglayers,j,k) = solnvar(lx,j,k)
                    enddo
                enddo
            endif
        endif
        !==============================================================


        !==================Bottom side=================================
        if(g_brank .lt. 0) then

            if(bc_codes(BOTTOM) .eq. 'DIRC') then
                do k=1,lz
                    do i=1,lx
                        solnvar(i,0,k) = TWO*bcvals%bottom(i,k) - solnvar(i,1,k)
                    enddo
                enddo
            else
                do k=1,lz
                    do i=1,lx
                        solnvar(i,0,k) = solnvar(i,1,k)
                    enddo
                enddo
            endif
        endif
        !==============================================================


        !==================Top side====================================
        if(g_trank .lt. 0) then

            if(bc_codes(TOP) .eq. 'DIRC') then
                do k=1,lz
                    do i=1,lx
                        solnvar(i,ly+g_nglayers,k) = TWO*bcvals%top(i,k) - solnvar(i,ly,k)
                    enddo
                enddo

            else
                do k=1,lz
                    do i=1,lx
                        solnvar(i,ly+g_nglayers,k) = solnvar(i,ly,k)
                    enddo
                enddo
            endif
        endif
        !==============================================================


        !==================Back side====================================
        if(g_krank .lt. 0) then
            if(bc_codes(BACK) .eq. 'DIRC') then

                do j=1,ly
                    do i=1,lx
                        solnvar(i,j,0) = TWO*bcvals%back(i,j) - solnvar(i,j,1)
                    enddo
                enddo
            else
                do j=1,ly
                    do i=1,lx
                        solnvar(i,j,0) = solnvar(i,j,1)
                    enddo
                enddo

            endif
        endif
        !==============================================================


        !================Front side====================================
        if(g_frank .lt. 0) then
            if(bc_codes(FRONT) .eq. 'DIRC') then

                do j=1,ly
                    do i=1,lx
                        solnvar(i,j,lz+g_nglayers) = TWO*bcvals%front(i,j) - solnvar(i,j,lz)
                    enddo
                enddo
            else
                do j=1,ly
                    do i=1,lx
                        solnvar(i,j,lz+g_nglayers) = solnvar(i,j,lz)
                    enddo
                enddo
            endif
        endif
        !==============================================================


    end subroutine update_boundary_ghostvalues
    !=================================================================================================================
    subroutine apply_null_bcs(solnvar,lx,ly,lz)

        integer, intent(in)   :: lx,ly,lz
        real*8, intent(inout) :: solnvar(-g_nglayers+1:lx+g_nglayers,&
            -g_nglayers+1:ly+g_nglayers,-g_nglayers+1:lz+g_nglayers)

        integer :: i,j,k

        if(g_lrank .lt. 0) then
            do k=1,lz
                do j=1,ly
                    solnvar(0,j,k) = ZERO
                enddo
            enddo
        endif

        if(g_rrank .lt. 0) then
            do k=1,lz
                do j=1,ly
                    solnvar(lx+g_nglayers,j,k) = ZERO
                enddo
            enddo
        endif

        if(g_brank .lt. 0) then
            do k=1,lz
                do i=1,lx
                    solnvar(i,0,k) = ZERO
                enddo
            enddo
        endif

        if(g_trank .lt. 0) then
            do k=1,lz
                do i=1,lx
                    solnvar(i,ly+g_nglayers,k) = ZERO
                enddo
            enddo
        endif

        if(g_krank .lt. 0) then
            do j=1,ly
                do i=1,lx
                    solnvar(i,j,0) = ZERO
                enddo
            enddo
        endif

        if(g_frank .lt. 0) then
            do j=1,ly
                do i=1,lx
                    solnvar(i,j,lz+g_nglayers) = ZERO
                enddo
            enddo
        endif

    end subroutine apply_null_bcs
    !=================================================================================================================
end module discretepde
