module pdedefs

    use globalvars
    use mgridsteps

    implicit none

    real*8, allocatable :: pdesoln(:,:,:)
    character(LEN=20)   :: pdescalarname

    real*8, allocatable, private    :: dcoeff(:,:,:)
    real*8, allocatable, private    :: vel(:,:,:,:)
    real*8, allocatable, private    :: source(:,:,:)
    real*8, allocatable, private    :: reac(:,:,:)
    character(LEN=4),private        :: pde_bc_codes(NFACES)
    type(boundarycondition),private :: pde_bcvals

    integer :: pde_laplsolver,pde_nvcycles

    logical,private :: timederivflag

contains
    !========================================================================
    subroutine pde_initialize(laplsolveflag,nvcycles,&
            bc_codes,bc_params)

        integer, intent(in) :: laplsolveflag
        integer, intent(in) :: nvcycles
        character(LEN=*)  :: bc_codes(NFACES)
        real*8,intent(in) :: bc_params(NFACES)
        integer :: i

        pde_laplsolver = laplsolveflag
        pde_nvcycles = nvcycles

        allocate(pdesoln(g_lx+2*g_nglayers,&
            g_ly+2*g_nglayers,&
            g_lz+2*g_nglayers))

        allocate(vel(g_lx+2*g_nglayers,&
            g_ly+2*g_nglayers,&
            g_lz+2*g_nglayers,NDIM))

        allocate(dcoeff(g_lx+2*g_nglayers,&
            g_ly+2*g_nglayers,&
            g_lz+2*g_nglayers))

        allocate(reac(g_lx+2*g_nglayers,&
            g_ly+2*g_nglayers,&
            g_lz+2*g_nglayers))

        allocate(source(g_lx+2*g_nglayers,&
            g_ly+2*g_nglayers,&
            g_lz+2*g_nglayers))


        pdesoln = ZERO
        vel     = ZERO
        dcoeff  = ZERO
        reac    = ZERO
        source  = ZERO

        timederivflag = .true.
        if(laplsolveflag .eq. 1) then
            timederivflag = .false.
        endif

        call pde_init()
        call pde_update_bcs(bc_codes,bc_params)

    end subroutine pde_initialize
    !=========================================================================
    subroutine pde_init()

        real*8  :: x,y,z
        integer :: i,j,k

        pdescalarname="Phi"
        pdescalarname=trim(pdescalarname)

        pdesoln=ZERO

        do k=2,g_lz+1
            do j=2,g_ly+1
                do i=2,g_lx+1

                    x = g_offx + (i-2)*g_dx + HALF*g_dx
                    y = g_offy + (j-2)*g_dy + HALF*g_dy
                    z = g_offz + (k-2)*g_dz + HALF*g_dz

                    pdesoln(i,j,k) = 0.5*(x**2-x)
                    !pdesoln(i,j,k) = 0.d0

                enddo
            enddo
        enddo



    end subroutine pde_init
    !=================================================================
    subroutine pde_update_bcs(bc_codes,bc_params)

        character(LEN=*)  :: bc_codes(NFACES)
        real*8,intent(in) :: bc_params(NFACES)
        integer :: i,j,k

        call allocate_bcvalues(pde_bcvals,g_lx,g_ly,g_lz)
        pde_bc_codes = bc_codes

        if(g_lrank .lt. 0) then

            i = 1
            if(bc_codes(LEFT) .eq. 'DIRC') then

                do k=1,g_lz
                    do j=1,g_ly
                        pde_bcvals%left(j,k)=bc_params(LEFT)
                    enddo
                enddo
            endif
            if(bc_codes(LEFT) .eq. 'FLUX') then

                do k=1,g_lz
                    do j=1,g_ly
                        pde_bcvals%left(j,k)=bc_params(LEFT)
                    enddo
                enddo
            endif
        endif

        if(g_rrank .lt. 0) then

            i = g_lx
            if(bc_codes(RIGHT) .eq. 'DIRC') then

                do k=1,g_lz
                    do j=1,g_ly
                        pde_bcvals%right(j,k)=bc_params(RIGHT)
                    enddo
                enddo
            endif
            if(bc_codes(RIGHT) .eq. 'FLUX') then

                do k=1,g_lz
                    do j=1,g_ly
                        pde_bcvals%right(j,k)=bc_params(RIGHT)
                    enddo
                enddo
            endif
        endif

        if(g_brank .lt. 0) then

            j = 1
            if(bc_codes(BOTTOM) .eq. 'DIRC') then

                do k=1,g_lz
                    do i=1,g_lx
                        pde_bcvals%bottom(i,k)=bc_params(BOTTOM)
                    enddo
                enddo
            endif
            if(bc_codes(BOTTOM) .eq. 'FLUX') then

                do k=1,g_lz
                    do i=1,g_lx
                        pde_bcvals%bottom(i,k)=bc_params(BOTTOM)
                    enddo
                enddo
            endif
        endif

        if(g_trank .lt. 0) then

            j = g_ly
            if(bc_codes(TOP) .eq. 'DIRC') then

                do k=1,g_lz
                    do i=1,g_lx
                        pde_bcvals%top(i,k)=bc_params(TOP)
                    enddo
                enddo
            endif

            if(bc_codes(TOP) .eq. 'FLUX') then

                do k=1,g_lz
                    do i=1,g_lx
                        pde_bcvals%top(i,k)=bc_params(TOP)
                    enddo
                enddo
            endif
        endif

        if(g_krank .lt. 0) then

            k = 1
            if(bc_codes(BACK) .eq. 'DIRC') then

                do j=1,g_ly
                    do i=1,g_lx
                        pde_bcvals%back(i,j)=bc_params(BACK)
                    enddo
                enddo
            endif

            if(bc_codes(BACK) .eq. 'FLUX') then

                do j=1,g_ly
                    do i=1,g_lx
                        pde_bcvals%back(i,j)=bc_params(BACK)
                    enddo
                enddo
            endif
        endif

        if(g_frank .lt. 0) then

            k = g_lz
            if(bc_codes(FRONT) .eq. 'DIRC') then

                do j=1,g_ly
                    do i=1,g_lx
                        pde_bcvals%front(i,j)=bc_params(FRONT)
                    enddo
                enddo
            endif

            if(bc_codes(FRONT) .eq. 'FLUX') then

                do j=1,g_ly
                    do i=1,g_lx
                        pde_bcvals%front(i,j)=bc_params(FRONT)
                    enddo
                enddo
            endif
        endif


    end subroutine pde_update_bcs
    !=================================================================
    subroutine pde_update_transport(dt)

        real*8, intent(in) :: dt
        real*8  :: x,y,z
        integer :: i,j,k


        do k=2,g_lz+1
            do j=2,g_ly+1
                do i=2,g_lx+1

                    x = g_offx + (i-2)*g_dx + HALF*g_dx
                    y = g_offy + (j-2)*g_dy + HALF*g_dy
                    z = g_offz + (k-2)*g_dz + HALF*g_dz

                    if(pde_laplsolver .eq. 1) then
                        dcoeff(i,j,k) = 1.d0
                    else
                        vel(i,j,k,XDIR) =  -1.d0	
                        vel(i,j,k,YDIR) =  0.d0	
                        vel(i,j,k,ZDIR) =  0.d0

                        dcoeff(i,j,k)   = 0.5d0
                        reac(i,j,k)     = 2.d0
                        source(i,j,k)   = -x**2
                    endif

                enddo
            enddo
        enddo

        if(timederivflag .eqv. .true.) source = source + pdesoln/dt


    end subroutine pde_update_transport
    !=================================================================
    subroutine pde_solve(dt)

        real*8, intent(in) :: dt

        integer :: i,ii,jj,kk
        integer :: nvcycles
        real*8 :: resnorm
        real*8,allocatable:: sterm(:,:,:)
        real*8,allocatable:: res(:,:,:)
        real*8 :: solve_time

        real*8 :: err_tol
        real*8 :: solnerr,x

        err_tol=ERRTOLR*0.001
        nvcycles=pde_nvcycles

        allocate(sterm(g_lx,g_ly,g_lz))
        allocate(res(g_lx,g_ly,g_lz))

        if(g_myproc .eq. g_rootproc) solve_time=-MPI_Wtime()

        call pde_update_transport(dt)

        do i=1,nvcycles

            sterm = 0.d0

            call perform_vcycle(pdesoln,timederivflag,dt,&
                g_lx,g_ly,g_lz,vel,dcoeff,reac,source,&
                pde_bc_codes,pde_bcvals,sterm,resnorm)

            !Sanity check; see if AX and b are found correctly=======================

            !call perform_gjacobi_smoothing(pdesoln,timederivflag,dt,g_lx,g_ly,g_lz,.true.,sterm,res,&
            !			vel,dcoeff,reac,source,&
            !			pde_bc_codes,pde_bcvals,&
            !			g_lrank,g_rrank,g_brank,g_trank,g_krank,g_frank,&
            !			g_llenx,g_lleny,g_llenz,1)

            !call compute_norm(res,g_lx,g_ly,g_lz,resnorm)

            !Sanity check; see if AX and b are found correctly=======================

            solnerr=0.d0
            do kk=2,g_lz+1
                do jj=2,g_ly+1
                    do ii=2,g_lx+1
                    
                        x = g_offx + (i-2)*g_dx + HALF*g_dx

                        if(pde_laplsolver .eq. 1) then
                                solnerr=solnerr+(pdesoln(ii,jj,kk)-1.d0)**2
                        else
                                solnerr=solnerr+(pdesoln(ii,jj,kk)-0.5d0*(x**2-x))**2
                        endif

                    enddo
                enddo
            enddo

            solnerr=solnerr/(g_lx*g_ly*g_lz)

            if(g_myproc .eq. g_rootproc) print *,"it:",i,resnorm,sqrt(solnerr)

            if(resnorm .le. err_tol) exit

        enddo

        if(g_myproc .eq. g_rootproc) solve_time=solve_time+MPI_Wtime()
        if(g_myproc .eq. g_rootproc) print *,"solve time:",solve_time

    end subroutine pde_solve
    !=================================================================
end module pdedefs
