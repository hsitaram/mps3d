module edeneq

      use par_decompose
      use globalvars
      use mgridsteps
      use chem_module

      implicit none

      real*8, private     :: nemin
      real*8, allocatable :: edensoln(:,:,:)
      character(LEN=20)   :: edenname
      
      real*8, allocatable, private    :: dcoeff(:,:,:)
      real*8, allocatable, private    :: vel(:,:,:,:)
      real*8, allocatable, private    :: source(:,:,:)
      real*8, allocatable, private    :: reac(:,:,:)
      character(LEN=4),private        :: eden_bc_codes(NFACES)
      type(boundarycondition),private :: eden_bcvals

      real*8, private :: lscale,tscale 
      real*8, private :: edenscale

      real*8, private :: vel_scaling
      real*8, private :: dcoeff_scaling
      real*8, private :: reac_scaling
      real*8, private :: source_scaling

      logical, private :: timederivflag
      

      contains
!========================================================================
	subroutine eden_initialize(eden_initial,bc_codes,bc_params,&
				length_scale,time_scale,var_scale)

	      real*8, intent(in) :: eden_initial(g_lx+2*g_nglayers,&
			      		 	 g_ly+2*g_nglayers,&
						 g_lz+2*g_nglayers)
	      
	      character(LEN=*), intent(in)  :: bc_codes(NFACES)
	      real*8,intent(in) :: bc_params(NFACES)
	      real*8, intent(in) :: length_scale,time_scale
	      real*8, intent(in) :: var_scale
	      
	      allocate(edensoln(g_lx+2*g_nglayers,&
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

	      nemin = 1.0d11

	      edenscale = var_scale
	      lscale   = length_scale
	      tscale   = time_scale

	      nemin = nemin/edenscale

	      edensoln = ZERO
	      vel     = ZERO
	      dcoeff  = ZERO
	      reac    = ZERO
	      source  = ZERO

	      !These are factors multiplied (not divided) with
	      !transport coefficients
	      vel_scaling    = tscale/lscale
	      dcoeff_scaling = tscale/(lscale**2)
	      reac_scaling   = tscale
	      source_scaling = tscale/edenscale 

	      timederivflag = .true.
	
	      call eden_init(eden_initial)
	      call eden_set_bcs(bc_codes,bc_params)

	end subroutine eden_initialize
!=========================================================================
      subroutine eden_init(eden_initial)
		
	        real*8, intent(in) :: eden_initial(g_lx+2*g_nglayers,&
		   	      		 	   g_ly+2*g_nglayers,&
  						   g_lz+2*g_nglayers)
		real*8  :: x,y,z
	        integer :: i,j,k

		edenname="Electron_density (#/m3)"
		edenname=trim(edenname)

		edensoln=ZERO

		do k=2,g_lz+1
			do j=2,g_ly+1
				do i=2,g_lx+1
					
					x = (g_offx + (i-2)*g_dx + HALF*g_dx)*lscale
					y = (g_offy + (j-2)*g_dy + HALF*g_dy)*lscale
					z = (g_offz + (k-2)*g_dz + HALF*g_dz)*lscale
					
					!g_solnvars(i,j,k,1) = x**2+y**2+z**2
					!edensoln(i,j,k) = 0.5*(x**2-x)
					edensoln(i,j,k) = 0.d0

				enddo
			enddo
		enddo

		edensoln = eden_initial/edenscale

      end subroutine eden_init
!=======================================================================================
      subroutine eden_set_bcs(bc_codes,bc_params)

		character(LEN=*)  :: bc_codes(NFACES)
		real*8,intent(in) :: bc_params(NFACES)

		integer :: i,j,k

	  	call allocate_bcvalues(eden_bcvals,g_lx,g_ly,g_lz)
		eden_bc_codes = bc_codes

		!Set dirichlet conditions from boundary
		!Flux/zgrd conditions dont require boundary terms
		if(g_lrank .lt. 0) then
	
			i = 1
			do k=1,g_lz
				do j=1,g_ly
					eden_bcvals%left(j,k)=bc_params(LEFT)/edenscale
				enddo
			enddo

		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			do k=1,g_lz
				do j=1,g_ly
					eden_bcvals%right(j,k)=bc_params(RIGHT)/edenscale
				enddo
			enddo
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			do k=1,g_lz
				do i=1,g_lx
					eden_bcvals%bottom(i,k)=bc_params(BOTTOM)/edenscale
				enddo
			enddo
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			do k=1,g_lz
				do i=1,g_lx
					eden_bcvals%top(i,k)=bc_params(TOP)/edenscale
				enddo
			enddo
		endif

		if(g_krank .lt. 0) then

			k = 1
			do j=1,g_ly
				do i=1,g_lx
					eden_bcvals%back(i,j)=bc_params(BACK)/edenscale
				enddo
			enddo
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			do j=1,g_ly
				do i=1,g_lx
					eden_bcvals%front(i,j)=bc_params(FRONT)/edenscale
				enddo
			enddo
		endif


      end subroutine eden_set_bcs
!===========================================================================================
      subroutine eden_update_bcs(elecden,electemp)

		real*8, intent(in) :: electemp(-g_nglayers+1:g_lx+g_nglayers,&
					       -g_nglayers+1:g_ly+g_nglayers,&
					       -g_nglayers+1:g_lz+g_nglayers)
		!elecden is scaled electron density
		real*8, intent(in) :: elecden (-g_nglayers+1:g_lx+g_nglayers,&
					       -g_nglayers+1:g_ly+g_nglayers,&
					       -g_nglayers+1:g_lz+g_nglayers)

		real*8 :: cbar,const_1
		integer :: i,j,k

		const_1 = 8.d0*k_BMAN/PI/M_ELEC

		if(g_lrank .lt. 0) then
	
			i = 1
			if(eden_bc_codes(LEFT) .eq. 'DIRC') then

				do k=1,g_lz
					do j=1,g_ly
						eden_bcvals%left(j,k)=eden_bcvals%left(j,k)
					enddo
				enddo

			else if(eden_bc_codes(LEFT) .eq. 'FLUX') then
				
				do k=1,g_lz
					do j=1,g_ly
					
						cbar = sqrt(const_1*electemp(i,j,k))	
						eden_bcvals%left(j,k) = 0.25*elecden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif

		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			if(eden_bc_codes(RIGHT) .eq. 'DIRC') then
	
				do k=1,g_lz
					do j=1,g_ly
						eden_bcvals%right(j,k)=eden_bcvals%right(j,k)
					enddo
				enddo

			else if(eden_bc_codes(RIGHT) .eq. 'FLUX') then
				
				do k=1,g_lz
					do j=1,g_ly
						
						cbar = sqrt(const_1*electemp(i,j,k))	
						eden_bcvals%right(j,k) = 0.25*elecden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			if(eden_bc_codes(BOTTOM) .eq. 'DIRC') then

				do k=1,g_lz
					do i=1,g_lx
						eden_bcvals%bottom(i,k)=eden_bcvals%bottom(i,k)
					enddo
				enddo

			else if(eden_bc_codes(BOTTOM) .eq. 'FLUX') then
				
				do k=1,g_lz
					do i=1,g_lx
						
						cbar = sqrt(const_1*electemp(i,j,k))	
						eden_bcvals%bottom(i,k) = 0.25*elecden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			if(eden_bc_codes(TOP) .eq. 'DIRC') then

				do k=1,g_lz
					do i=1,g_lx
						eden_bcvals%top(i,k)=eden_bcvals%top(i,k)
					enddo
				enddo
			
			else if(eden_bc_codes(TOP) .eq. 'FLUX') then
				
				do k=1,g_lz
					do i=1,g_lx
						
						cbar = sqrt(const_1*electemp(i,j,k))	
						eden_bcvals%top(i,k) = 0.25*elecden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif

		if(g_krank .lt. 0) then

			k = 1
			if(eden_bc_codes(BACK) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						eden_bcvals%back(i,j)=eden_bcvals%back(i,j)
					enddo
				enddo
			
			else if(eden_bc_codes(BACK) .eq. 'FLUX') then
				
				do j=1,g_ly
					do i=1,g_lx
						
						cbar = sqrt(const_1*electemp(i,j,k))	
						eden_bcvals%back(i,j) = 0.25*elecden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			if(eden_bc_codes(FRONT) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						eden_bcvals%front(i,j)=eden_bcvals%front(i,j)
					enddo
				enddo
			
			else if(eden_bc_codes(FRONT) .eq. 'FLUX') then
				
				do j=1,g_ly
					do i=1,g_lx
						
						cbar = sqrt(const_1*electemp(i,j,k))	
						eden_bcvals%front(i,j) = 0.25*elecden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif


      end subroutine eden_update_bcs
!===========================================================================================
      subroutine eden_update_transport(dt,elecfield,electemp,elecprod,numden,Tg,Pg)
	
		real*8, intent(in) :: dt

		real*8, intent(in) :: elecfield(-g_nglayers+1:g_lx+g_nglayers,&
						-g_nglayers+1:g_ly+g_nglayers,&
						-g_nglayers+1:g_lz+g_nglayers,NDIM)

		real*8, intent(in) :: electemp(-g_nglayers+1:g_lx+g_nglayers,&
					       -g_nglayers+1:g_ly+g_nglayers,&
					       -g_nglayers+1:g_lz+g_nglayers)
	
		real*8, intent(inout) :: elecprod(-g_nglayers+1:g_lx+g_nglayers,&
					          -g_nglayers+1:g_ly+g_nglayers,&
					          -g_nglayers+1:g_lz+g_nglayers)

		real*8, intent(in) :: numden(-g_nglayers+1:g_lx+g_nglayers,&
					     -g_nglayers+1:g_ly+g_nglayers,&
					     -g_nglayers+1:g_lz+g_nglayers,nspecies)
		
		real*8, intent(in) :: Tg,Pg
	
		real*8 :: efield_mag
		real*8 :: mu_e


		real*8  :: x,y,z
	        integer :: i,j,k,nsp

		real*8 :: specarray(nspecies)
		real*8 :: specprod

		vel = ZERO
		dcoeff = ZERO
		reac = ZERO
		source = ZERO

		do k=1,g_lz
			do j=1,g_ly
				do i=1,g_lx

					do nsp=1,nspecies
						specarray(nsp)=numden(i,j,k,nsp)
					enddo

					efield_mag = elecfield(i,j,k,XDIR)**2
					efield_mag = efield_mag + elecfield(i,j,k,YDIR)**2
					efield_mag = efield_mag + elecfield(i,j,k,ZDIR)**2

					efield_mag = sqrt(efield_mag)

					mu_e      =   getspecmobility(especnum,specarray,efield_mag,&
								electemp(i,j,k),Tg,Pg)
					
					vel(i+1,j+1,k+1,XDIR) = mu_e*elecfield(i,j,k,XDIR)
					vel(i+1,j+1,k+1,YDIR) = mu_e*elecfield(i,j,k,YDIR)
					vel(i+1,j+1,k+1,ZDIR) = mu_e*elecfield(i,j,k,ZDIR)
		
					dcoeff(i+1,j+1,k+1)   = getspecdcoeff(especnum,specarray,efield_mag,&
									electemp(i,j,k),Tg,Pg)
					
					call getspecproduction(especnum,electemp(i,j,k),Tg,&
								specarray,specprod,efield_mag)

					source(i+1,j+1,k+1)   =   source(i+1,j+1,k+1) + specprod

				enddo
			enddo
		enddo
		
		!vel = ZERO
		!dcoeff = ZERO
		!reac = ZERO
		!source = ZERO
		elecprod(1:g_lx, 1:g_ly, 1:g_lz) = source(2:g_lx+1, 2:g_ly+1, 2:g_lz+1)

		dcoeff   = dcoeff * dcoeff_scaling
		vel      = vel * vel_scaling
		reac     = reac * reac_scaling
		source   = source * source_scaling
		
		if(timederivflag .eqv. .true.) source = source + edensoln/dt


      end subroutine eden_update_transport
!=========================================================================================
      subroutine eden_solve(time,dt,elecfield,electemp,numden,elecprod,Tg,Pg,printflag)

	real*8, intent(in) :: time
	real*8, intent(in) :: dt,Tg,Pg
	
	logical, intent(in) :: printflag
	real*8, intent(in) :: elecfield(-g_nglayers+1:g_lx+g_nglayers,&
					-g_nglayers+1:g_ly+g_nglayers,&
					-g_nglayers+1:g_lz+g_nglayers,NDIM)

	real*8, intent(in) :: electemp(-g_nglayers+1:g_lx+g_nglayers,&
				       -g_nglayers+1:g_ly+g_nglayers,&
				       -g_nglayers+1:g_lz+g_nglayers)
	

	real*8, intent(in) :: numden(-g_nglayers+1:g_lx+g_nglayers,&
				     -g_nglayers+1:g_ly+g_nglayers,&
				     -g_nglayers+1:g_lz+g_nglayers,nspecies)
	
	real*8, intent(inout) :: elecprod(-g_nglayers+1:g_lx+g_nglayers,&
				          -g_nglayers+1:g_ly+g_nglayers,&
				          -g_nglayers+1:g_lz+g_nglayers)

	integer :: i,j,k
	integer :: nvcycles
	real*8 :: resnorm
	real*8,allocatable:: sterm(:,:,:)
	real*8,allocatable:: res(:,:,:)
	real*8 :: solve_time

	real*8 :: err_tol
	real*8 :: dt_s
	
	dt_s = dt/tscale
	
	err_tol=ERRTOLR
	nvcycles=20
	
	allocate(sterm(g_lx,g_ly,g_lz))
	allocate(res(g_lx,g_ly,g_lz))
	      
	if(g_myproc .eq. g_rootproc) solve_time=-MPI_Wtime()

	call eden_update_transport(dt_s,elecfield,electemp,elecprod,numden,Tg,Pg)
	call eden_update_bcs(edensoln,electemp)

	do i=1,nvcycles

		sterm = 0.d0
		
		call perform_vcycle(edensoln,timederivflag,dt_s,g_lx,g_ly,g_lz,vel,dcoeff,reac,source,&
				eden_bc_codes,eden_bcvals,sterm,resnorm)
	
		if((printflag .eqv. .true.) .and. (g_myproc .eq. g_rootproc)) print *,"it:",i,resnorm

		if(resnorm .le. err_tol) exit
	
	enddo

	if(g_myproc .eq. g_rootproc) solve_time=solve_time+MPI_Wtime()
	if((g_myproc .eq. g_rootproc) .and. (printflag .eqv. .true.)) print *,"solve time:",solve_time


	!floor electron density
	do k=1,g_lz
		do j=1,g_ly
			do i=1,g_lx
		       		
				if(edensoln(i+1,j+1,k+1) .lt. nemin) then
					edensoln(i+1,j+1,k+1)=nemin
				endif
			enddo
		enddo
	enddo	

	call exchangehalodata(edensoln,g_lx,g_ly,g_lz)
	call update_boundary_ghostvalues(edensoln,eden_bc_codes,&
				eden_bcvals,g_lx,g_ly,g_lz)

      end subroutine eden_solve
!=================================================================
end module edeneq
