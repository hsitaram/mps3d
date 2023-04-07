module neutralseq

      use globalvars
      use mgridsteps
      use chem_module

      implicit none

      real*8, allocatable             :: neutralsoln(:,:,:,:)
      character(LEN=20),allocatable   :: neutralnames(:)
      
      real*8, allocatable, private    :: dcoeff(:,:,:)
      real*8, allocatable, private    :: vel(:,:,:,:)
      real*8, allocatable, private    :: source(:,:,:)
      real*8, allocatable, private    :: reac(:,:,:)
      character(LEN=4),private        :: neutrals_bc_codes(NFACES)
      type(boundarycondition),private :: neutrals_bcvals

      real*8, private :: lscale,tscale 
      real*8, private :: neutralscale

      real*8, private :: vel_scaling
      real*8, private :: dcoeff_scaling
      real*8, private :: reac_scaling
      real*8, private :: source_scaling

      logical, private :: timederivflag
      

      contains
!========================================================================
	subroutine neutrals_initialize(neutrals_initial,bc_codes,bc_params,&
				length_scale,time_scale,var_scale)

	      real*8, intent(in) :: neutrals_initial(g_lx+2*g_nglayers,&
			      		     	     g_ly+2*g_nglayers,&
						     g_lz+2*g_nglayers,&
						     solved_neutrals_num)
	
	      character(LEN=*), intent(in)  :: bc_codes(NFACES)
	      real*8,intent(in) :: bc_params(NFACES)
	      real*8, intent(in) :: length_scale,time_scale
	      real*8, intent(in) :: var_scale
	      
	      allocate(neutralsoln(g_lx+2*g_nglayers,&
			       g_ly+2*g_nglayers,&
			       g_lz+2*g_nglayers,solved_neutrals_num))
	      
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


	      neutralscale = var_scale
	      lscale   = length_scale
	      tscale   = time_scale

	      neutralsoln = ZERO
	      vel     = ZERO
	      dcoeff  = ZERO
	      reac    = ZERO
	      source  = ZERO

	      !These are factors multiplied (not divided) with
	      !transport coefficients
	      vel_scaling    = tscale/lscale
	      dcoeff_scaling = tscale/(lscale**2)
	      reac_scaling   = tscale
	      source_scaling = tscale/neutralscale 

	      timederivflag = .true.
	
	      call neutrals_init(neutrals_initial)
	      call neutrals_set_bcs(bc_codes,bc_params)

	end subroutine neutrals_initialize
!=========================================================================
      subroutine neutrals_init(neutrals_initial)
		
	        real*8, intent(in) :: neutrals_initial(g_lx+2*g_nglayers,&
		  	      		     	       g_ly+2*g_nglayers,&
						       g_lz+2*g_nglayers,&
						       solved_neutrals_num)

		real*8  :: x,y,z
	        integer :: i,j,k,nsp

		allocate(neutralnames(solved_neutrals_num))

		do nsp=1,solved_neutrals_num
			neutralnames(nsp) = specnames(nsp-1+neutralspecmin)//'_Density (#/m3)'
			neutralnames(nsp) = trim(neutralnames(nsp))
		enddo

		neutralsoln=ZERO

		do k=2,g_lz+1
			do j=2,g_ly+1
				do i=2,g_lx+1
					
					x = (g_offx + (i-2)*g_dx + HALF*g_dx)*lscale
					y = (g_offy + (j-2)*g_dy + HALF*g_dy)*lscale
					z = (g_offz + (k-2)*g_dz + HALF*g_dz)*lscale
					
					!g_solnvars(i,j,k,1) = x**2+y**2+z**2
					!ionsoln(i,j,k) = 0.5*(x**2-x)
					!ionsoln(i,j,k) = 0.d0

				enddo
			enddo
		enddo

		neutralsoln = neutrals_initial/neutralscale


      end subroutine neutrals_init
!=======================================================================================
      subroutine neutrals_set_bcs(bc_codes,bc_params)

		character(LEN=*)  :: bc_codes(NFACES)
		real*8,intent(in) :: bc_params(NFACES)

		integer :: i,j,k

	  	call allocate_bcvalues(neutrals_bcvals,g_lx,g_ly,g_lz)
		neutrals_bc_codes = bc_codes

		!Set dirichlet conditions from boundary
		!Flux/zgrd conditions dont require boundary terms
		if(g_lrank .lt. 0) then
	
			i = 1
			do k=1,g_lz
				do j=1,g_ly
					neutrals_bcvals%left(j,k)=bc_params(LEFT)/neutralscale
				enddo
			enddo

		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			do k=1,g_lz
				do j=1,g_ly
					neutrals_bcvals%right(j,k)=bc_params(RIGHT)/neutralscale
				enddo
			enddo
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			do k=1,g_lz
				do i=1,g_lx
					neutrals_bcvals%bottom(i,k)=bc_params(BOTTOM)/neutralscale
				enddo
			enddo
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			do k=1,g_lz
				do i=1,g_lx
					neutrals_bcvals%top(i,k)=bc_params(TOP)/neutralscale
				enddo
			enddo
		endif

		if(g_krank .lt. 0) then

			k = 1
			do j=1,g_ly
				do i=1,g_lx
					neutrals_bcvals%back(i,j)=bc_params(BACK)/neutralscale
				enddo
			enddo
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			do j=1,g_ly
				do i=1,g_lx
					neutrals_bcvals%front(i,j)=bc_params(FRONT)/neutralscale
				enddo
			enddo
		endif


      end subroutine neutrals_set_bcs
!===========================================================================================
      subroutine neutrals_update_bcs(neutralnum,neutralden,Tg,Pg)

		real*8,  intent(in) :: Tg,Pg
		integer, intent(in) :: neutralnum

		!neutralden is scaled (non-dimesional)
		real*8, intent(in) :: neutralden(-g_nglayers+1:g_lx+g_nglayers,&
						 -g_nglayers+1:g_ly+g_nglayers,&
						 -g_nglayers+1:g_lz+g_nglayers)

		integer :: i,j,k,neutralspecnum,nsp
		
		real*8  :: cbar,const_1
		real*8  :: specarray(nspecies)
		real*8  :: efield_mag

		neutralspecnum = neutralnum - 1 + neutralspecmin

		const_1 = 8.d0*k_BMAN/PI/molmass(neutralspecnum)
		cbar    = sqrt(const_1*Tg)

		if(g_lrank .lt. 0) then
	
			i = 1
			if(neutrals_bc_codes(LEFT) .eq. 'DIRC') then

				do k=1,g_lz
					do j=1,g_ly
						neutrals_bcvals%left(j,k)=neutrals_bcvals%left(j,k)
					enddo
				enddo

			else if(neutrals_bc_codes(LEFT) .eq. 'FLUX') then
				
				
				do k=1,g_lz
					do j=1,g_ly
						neutrals_bcvals%left(j,k)  = 0.25*neutralden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif

		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			if(neutrals_bc_codes(RIGHT) .eq. 'DIRC') then
	
				do k=1,g_lz
					do j=1,g_ly
						neutrals_bcvals%right(j,k)=neutrals_bcvals%right(j,k)
					enddo
				enddo

			else if(neutrals_bc_codes(RIGHT) .eq. 'FLUX') then
				
	
				do k=1,g_lz
					do j=1,g_ly
						neutrals_bcvals%right(j,k)  = 0.25*neutralden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			if(neutrals_bc_codes(BOTTOM) .eq. 'DIRC') then

				do k=1,g_lz
					do i=1,g_lx
						neutrals_bcvals%bottom(i,k)=neutrals_bcvals%bottom(i,k)
					enddo
				enddo

			else if(neutrals_bc_codes(BOTTOM) .eq. 'FLUX') then
				
				do k=1,g_lz
					do i=1,g_lx

						neutrals_bcvals%bottom(i,k)  = 0.25*neutralden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			if(neutrals_bc_codes(TOP) .eq. 'DIRC') then

				do k=1,g_lz
					do i=1,g_lx
						neutrals_bcvals%top(i,k)=neutrals_bcvals%top(i,k)
					enddo
				enddo
			
			else if(neutrals_bc_codes(TOP) .eq. 'FLUX') then
				
				do k=1,g_lz
					do i=1,g_lx
						neutrals_bcvals%top(i,k)  = 0.25*neutralden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif

		if(g_krank .lt. 0) then

			k = 1
			if(neutrals_bc_codes(BACK) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						neutrals_bcvals%back(i,j)=neutrals_bcvals%back(i,j)
					enddo
				enddo
			
			else if(neutrals_bc_codes(BACK) .eq. 'FLUX') then
				
				do j=1,g_ly
					do i=1,g_lx
						neutrals_bcvals%back(i,j)  = 0.25*neutralden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			if(neutrals_bc_codes(FRONT) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						neutrals_bcvals%front(i,j)=neutrals_bcvals%front(i,j)
					enddo
				enddo
			
			else if(neutrals_bc_codes(FRONT) .eq. 'FLUX') then
				
				do j=1,g_ly
					do i=1,g_lx
						neutrals_bcvals%front(i,j)  = 0.25*neutralden(i,j,k)*cbar*vel_scaling
					enddo
				enddo
			endif
		endif


      end subroutine neutrals_update_bcs
!===========================================================================================
      subroutine neutrals_update_transport(neutralnum,dt,elecfield,electemp,neutralprod,numden,Tg,Pg)
	
		integer, intent(in) :: neutralnum
		real*8, intent(in) :: dt

		real*8, intent(in) :: elecfield(-g_nglayers+1:g_lx+g_nglayers,&
						-g_nglayers+1:g_ly+g_nglayers,&
						-g_nglayers+1:g_lz+g_nglayers,NDIM)

		real*8, intent(in) :: electemp(-g_nglayers+1:g_lx+g_nglayers,&
					       -g_nglayers+1:g_ly+g_nglayers,&
					       -g_nglayers+1:g_lz+g_nglayers)

		real*8, intent(in) :: numden(-g_nglayers+1:g_lx+g_nglayers,&
					     -g_nglayers+1:g_ly+g_nglayers,&
					     -g_nglayers+1:g_lz+g_nglayers,nspecies)
	
		real*8, intent(inout) :: neutralprod(-g_nglayers+1:g_lx+g_nglayers,&
				   	             -g_nglayers+1:g_ly+g_nglayers,&
					             -g_nglayers+1:g_lz+g_nglayers)
		
		real*8, intent(in) :: Tg,Pg

		real*8 :: efield_mag
		real*8 :: mu_i

		real*8  :: x,y,z
	        integer :: i,j,k,nsp
		integer :: neutralspecnum

		real*8 :: specarray(nspecies)
		real*8 :: specprod

		vel = ZERO
		dcoeff = ZERO
		reac = ZERO
		source = ZERO

		neutralspecnum = neutralnum - 1 + neutralspecmin

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

					dcoeff(i+1,j+1,k+1)   = getspecdcoeff(neutralspecnum,specarray,efield_mag,&
									electemp(i,j,k),Tg,Pg)
					
					call getspecproduction(neutralspecnum,electemp(i,j,k),Tg,&
								specarray,specprod,efield_mag)

					source(i+1,j+1,k+1)   =   source(i+1,j+1,k+1) + specprod

				enddo
			enddo
		enddo

		neutralprod(1:g_lx, 1:g_ly, 1:g_lz) = source(2:g_lx+1, 2:g_ly+1, 2:g_lz+1)
		
		dcoeff   = dcoeff * dcoeff_scaling
		vel      = vel * vel_scaling
		reac     = reac * reac_scaling
		source   = source * source_scaling
		
		if(timederivflag .eqv. .true.) source = source + neutralsoln(:,:,:,neutralnum)/dt
		

      end subroutine neutrals_update_transport
!=========================================================================================
      subroutine neutrals_solve(time,dt,elecfield,electemp,numden,neutralprod,Tg,Pg,printflag)

	real*8, intent(in) :: time,dt,Tg,Pg
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
	
	real*8, intent(inout) :: neutralprod(-g_nglayers+1:g_lx+g_nglayers,&
			   	             -g_nglayers+1:g_ly+g_nglayers,&
				             -g_nglayers+1:g_lz+g_nglayers,solved_neutrals_num)

	integer :: i,nsp
	integer :: nvcycles
	real*8 :: resnorm
	real*8,allocatable:: sterm(:,:,:)
	real*8 :: solve_time

	real*8 :: err_tol,dt_s
	
	err_tol=ERRTOLR
	nvcycles=20

	dt_s = dt/tscale
	
	allocate(sterm(g_lx,g_ly,g_lz))
	      

	do nsp=1,solved_neutrals_num

		if(g_myproc .eq. g_rootproc) solve_time=-MPI_Wtime()
	
		if(printflag .eqv. .true.) call printline("solving "//trim(specnames(neutralspecmin+nsp-1))//" Density-------------")
		call neutrals_update_transport(nsp,dt_s,elecfield,electemp,neutralprod(:,:,:,nsp),numden,Tg,Pg)
		call neutrals_update_bcs(nsp,neutralsoln(:,:,:,nsp),Tg,Pg)
	
	
		do i=1,nvcycles

			sterm = 0.d0
		
			call perform_vcycle(neutralsoln(:,:,:,nsp),timederivflag,dt_s,&
					g_lx,g_ly,g_lz,vel,dcoeff,reac,source,&
					neutrals_bc_codes,neutrals_bcvals,sterm,resnorm)

			if((printflag .eqv. .true.) .and. (g_myproc .eq. g_rootproc)) print *,"it:",i,resnorm

			if(resnorm .le. err_tol) exit
	
		enddo
		
		call exchangehalodata(neutralsoln(:,:,:,nsp),g_lx,g_ly,g_lz)
		call update_boundary_ghostvalues(neutralsoln(:,:,:,nsp),neutrals_bc_codes,&
				neutrals_bcvals,g_lx,g_ly,g_lz)

		if(g_myproc .eq. g_rootproc) solve_time=solve_time+MPI_Wtime()
		if((g_myproc .eq. g_rootproc) .and. (printflag .eqv. .true.))  print *,"solve time:",solve_time
		
		if(printflag .eqv. .true.) call printline("------------------------------------------------------------------")
	enddo

      end subroutine neutrals_solve
!=================================================================
end module neutralseq
