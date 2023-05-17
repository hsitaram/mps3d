module solve_manager

    use globalvars
    use cartwrite
    use pdedefs
    use par_decompose
    use inputs 
    implicit none

    integer, private :: num_eq
    real*8,allocatable,private :: all_solnvars(:,:,:,:)
    character(LEN=MAXSTRSIZE),allocatable,private :: all_scalarnames(:)
    character(LEN=4),private :: solver_bc_codes(NFACES)
    real*8, private :: solver_bcparams(NFACES)
    real*8, private :: bcparams(NFACES)
    integer,private :: solve_nscalars
    integer :: laplsolveflag
    integer :: nvcycles
    real*8  :: dt,tfinal

contains
    !===================================================================================
    subroutine solversetup()


        call readinputfile("INPUTS",solver_bc_codes,solver_bcparams)
        call readsolverinpfile("SOLVER")
        call domaindecompose()

        num_eq = 1

        call pde_initialize(laplsolveflag,nvcycles,&
            solver_bc_codes,solver_bcparams)

        allocate(all_solnvars(g_lx+2*g_nglayers,&
            g_ly+2*g_nglayers,&
            g_lz+2*g_nglayers,num_eq))

        allocate(all_scalarnames(num_eq))


    end subroutine solversetup
    !===================================================================================
    subroutine readsolverinpfile(solverfilename)

        integer :: fnum
        character(LEN=*), intent(in) :: solverfilename
        character(LEN=MAXSTRSIZE) :: temp

        fnum=14

        open(unit=fnum,file=solverfilename)
        read(fnum,*) temp
        read(fnum,*) temp
        read(fnum,*) temp

        read(fnum,*) temp,laplsolveflag
        read(fnum,*) temp,dt
        read(fnum,*) temp,tfinal
        read(fnum,*) temp,nvcycles

        close(fnum)

    end subroutine readsolverinpfile
    !============================================================
    subroutine writeoutputfile(fname)

        character(LEN=*),intent(in) :: fname

        all_solnvars(:,:,:,num_eq)    = pdesoln
        all_scalarnames(num_eq)       = pdescalarname
        call writeoutput(fname,all_solnvars,num_eq,all_scalarnames,&
            g_lx,g_ly,g_lz)

    end subroutine writeoutputfile
    !===================================================================================
    subroutine timestepping()

        real*8 :: t

        t=ZERO
        call writeoutputfile('output0')

        do while(t .lt. tfinal)

            call pde_solve(dt)
            t=t+dt
            if(g_myproc .eq. g_rootproc) print *,"=====time:",t,"=========="

        enddo

        call writeoutputfile('output1')

    end subroutine timestepping
    !===================================================================================

end module solve_manager
