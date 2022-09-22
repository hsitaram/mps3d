module ionseq

      use par_decompose
      use globalvars
      use mgridsteps
      use chem_module

      implicit none

      real*8, allocatable                     :: ionsoln(:,:,:,:)
      character(LEN=MAXSTRSIZE),allocatable   :: ionnames(:)
      
      real*8, allocatable, private    :: dcoeff(:,:,:)
      real*8, allocatable, private    :: vel(:,:,:,:)
      real*8, allocatable, private    :: source(:,:,:)
      real*8, allocatable, private    :: reac(:,:,:)
      character(LEN=4),private        :: ions_bc_codes(NFACES)
      type(boundarycondition),private :: ions_bcvals

      real*8, private :: lscale,tscale 
      real*8, private :: ionscale

      real*8, private :: vel_scaling
      real*8, private :: dcoeff_scaling
      real*8, private :: reac_scaling
      real*8, private :: source_scaling

      logical, private :: timederivflag
      

      contains
!========================================================================
	subroutine ions_initialize(ions_initial,bc_codes,bc_params,&
				length_scale,time_scale,var_scale)

	      real*8, intent(in) :: ions_initial(g_lx+2*g_nglayers,&
			      		 	 g_ly+2*g_nglayers,&
						 g_lz+2*g_nglayers,&
						 solved_ions_num)
	      
	      character(LEN=*), intent(in)  :: bc_codes(NFACES)
	      real*8,intent(in) :: bc_params(NFACES)
	      real*8, intent(in) :: length_scale,time_scale
	      real*8, intent(in) :: var_scale
	      
	      allocate(ionsoln(g_lx+2*g_nglayers,&
			       g_ly+2*g_nglayers,&
			       g_lz+2*g_nglayers,solved_ions_num))
	      
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



	      ionscale = var_scale
	      lscale   = length_scale
	      tscale   = time_scale

	      ionsoln = ZERO
	      vel     = ZERO
	      dcoeff  = ZERO
	      reac    = ZERO
	      source  = ZERO

	      !These are factors multiplied (not divided) with
	      !transport coefficients
	      vel_scaling    = tscale/lscale
	      dcoeff_scaling = tscale/(lscale**2)
	      reac_scaling   = tscale
	      source_scaling = tscale/ionscale 

	      timederivflag = .true.
	
	      call ions_init(ions_initial)
	      call ions_set_bcs(bc_codes,bc_params)

	end subroutine ions_initialize
!=========================================================================
      subroutine ions_init(ions_initial)
		
	        real*8, intent(in) :: ions_initial(g_lx+2*g_nglayers,&
		  	      		 	   g_ly+2*g_nglayers,&
						   g_lz+2*g_nglayers,&
						   solved_ions_num)
		real*8  :: x,y,z
	        integer :: i,j,k,nsp
	      
		allocate(ionnames(solved_ions_num))

		do nsp=1,solved_ions_num
			ionnames(nsp)=specnames(nsp-1+ionspecmin)//'_Density (#/m3)'
			ionnames(nsp)=trim(ionnames(nsp))
		enddo

		ionsoln=ZERO

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

		ionsoln = ions_initial/ionscale

      end subroutine ions_init
!=======================================================================================
      subroutine ions_set_bcs(bc_codes,bc_params)

		character(LEN=*)  :: bc_codes(NFACES)
		real*8,intent(in) :: bc_params(NFACES)

		integer :: i,j,k

	  	call allocate_bcvalues(ions_bcvals,g_lx,g_ly,g_lz)
		ions_bc_codes = bc_codes

		!Set dirichlet conditions from boundary
		!Flux/zgrd conditions dont require boundary terms
		if(g_lrank .lt. 0) then
	
			i = 1
			do k=1,g_lz
				do j=1,g_ly
					ions_bcvals%left(j,k)=bc_params(LEFT)/ionscale
				enddo
			enddo

		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			do k=1,g_lz
				do j=1,g_ly
					ions_bcvals%right(j,k)=bc_params(RIGHT)/ionscale
				enddo
			enddo
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			do k=1,g_lz
				do i=1,g_lx
					ions_bcvals%bottom(i,k)=bc_params(BOTTOM)/ionscale
				enddo
			enddo
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			do k=1,g_lz
				do i=1,g_lx
					ions_bcvals%top(i,k)=bc_params(TOP)/ionscale
				enddo
			enddo
		endif

		if(g_krank .lt. 0) then
			
			k = 1
			do j=1,g_ly
				do i=1,g_lx
					ions_bcvals%back(i,j)=bc_params(BACK)/ionscale
				enddo
			enddo
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			do j=1,g_ly
				do i=1,g_lx
					ions_bcvals%front(i,j)=bc_params(FRONT)/ionscale
				enddo
			enddo
		endif


      end subroutine ions_set_bcs
!===========================================================================================
      subroutine ions_update_bcs(ionnum,ionden,elecfield,numden,electemp,Tg,Pg)

		real*8, intent(in) :: Tg,Pg
		integer, intent(in) :: ionnum

		real*8, intent(in) :: elecfield(-g_nglayers+1:g_lx+g_nglayers,&
						-g_nglayers+1:g_ly+g_nglayers,&
						-g_nglayers+1:g_lz+g_nglayers,NDIM)

		real*8, intent(in) :: electemp(-g_nglayers+1:g_lx+g_nglayers,&
					       -g_nglayers+1:g_ly+g_nglayers,&
					       -g_nglayers+1:g_lz+g_nglayers)
		
		!ionden is scaled (non-dimensional)
		real*8, intent(in) :: ionden(-g_nglayers+1:g_lx+g_nglayers,&
					     -g_nglayers+1:g_ly+g_nglayers,&
					     -g_nglayers+1:g_lz+g_nglayers)

		real*8, intent(in) :: numden(-g_nglayers+1:g_lx+g_nglayers,&
					     -g_nglayers+1:g_ly+g_nglayers,&
					     -g_nglayers+1:g_lz+g_nglayers,nspecies)

		integer :: i,j,k,ionspecnum,nsp
		
		real*8  :: mu_i
		real*8  :: specarray(nspecies)
		real*8  :: efield_mag

		ionspecnum = ionnum - 1 + ionspecmin

		if(g_lrank .lt. 0) then
	
			i = 1
			if(ions_bc_codes(LEFT) .eq. 'DIRC') then

				do k=1,g_lz
					do j=1,g_ly
						ions_bcvals%left(j,k)=ions_bcvals%left(j,k)/ionscale
					enddo
				enddo

			else if(ions_bc_codes(LEFT) .eq. 'FLUX') then
				
				
				do k=1,g_lz
					do j=1,g_ly
				
						do nsp=1,nspecies
							specarray(nsp)=numden(i,j,k,nsp)
						enddo
				
						efield_mag = elecfield(i,j,k,XDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,YDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,ZDIR)**2
						efield_mag = sqrt(efield_mag)
	
						mu_i = getspecmobility(ionspecnum,specarray,efield_mag,&
								electemp(i,j,k),Tg,Pg)

						if(mu_i*elecfield(i,j,k,XDIR) .lt. 0) then
							!outward flux has to be positive
							ions_bcvals%left(j,k) = -mu_i*ionden(i,j,k)*elecfield(i,j,k,XDIR)*vel_scaling
						else
							ions_bcvals%left(j,k)  = ZERO
						endif
					enddo
				enddo
			endif

		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			if(ions_bc_codes(RIGHT) .eq. 'DIRC') then
	
				do k=1,g_lz
					do j=1,g_ly
						ions_bcvals%right(j,k)=ions_bcvals%right(j,k)/ionscale
					enddo
				enddo

			else if(ions_bc_codes(RIGHT) .eq. 'FLUX') then
				
	
				do k=1,g_lz
					do j=1,g_ly
				
						do nsp=1,nspecies
							specarray(nsp)=numden(i,j,k,nsp)
						enddo
				
						efield_mag = elecfield(i,j,k,XDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,YDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,ZDIR)**2
						efield_mag = sqrt(efield_mag)
						
						mu_i = getspecmobility(ionspecnum,specarray,efield_mag,&
								electemp(i,j,k),Tg,Pg)

						if(mu_i*elecfield(i,j,k,XDIR) .gt. 0) then
							!outward flux has to be positive
							ions_bcvals%right(j,k) = mu_i*ionden(i,j,k)*elecfield(i,j,k,XDIR)*vel_scaling
						else
							ions_bcvals%right(j,k)  = ZERO
						endif
					enddo
				enddo
			endif
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			if(ions_bc_codes(BOTTOM) .eq. 'DIRC') then

				do k=1,g_lz
					do i=1,g_lx
						ions_bcvals%bottom(i,k)=ions_bcvals%bottom(i,k)/ionscale
					enddo
				enddo

			else if(ions_bc_codes(BOTTOM) .eq. 'FLUX') then
				
				do k=1,g_lz
					do i=1,g_lx

						do nsp=1,nspecies
							specarray(nsp)=numden(i,j,k,nsp)
						enddo
				
						efield_mag = elecfield(i,j,k,XDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,YDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,ZDIR)**2
						efield_mag = sqrt(efield_mag)
						
						mu_i = getspecmobility(ionspecnum,specarray,efield_mag,&
								electemp(i,j,k),Tg,Pg)

						if(mu_i*elecfield(i,j,k,YDIR) .lt. 0) then
							!outward flux has to be positive
							ions_bcvals%bottom(i,k) = -mu_i*ionden(i,j,k)*elecfield(i,j,k,YDIR)*vel_scaling
						else
							ions_bcvals%bottom(i,k)  = ZERO
						endif
					enddo
				enddo
			endif
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			if(ions_bc_codes(TOP) .eq. 'DIRC') then

				do k=1,g_lz
					do i=1,g_lx
						ions_bcvals%top(i,k)=ions_bcvals%top(i,k)/ionscale
					enddo
				enddo
			
			else if(ions_bc_codes(TOP) .eq. 'FLUX') then
				
				do k=1,g_lz
					do i=1,g_lx

						do nsp=1,nspecies
							specarray(nsp)=numden(i,j,k,nsp)
						enddo
				
						efield_mag = elecfield(i,j,k,XDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,YDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,ZDIR)**2
						efield_mag = sqrt(efield_mag)
						
						mu_i = getspecmobility(ionspecnum,specarray,efield_mag,&
								electemp(i,j,k),Tg,Pg)

						if(mu_i*elecfield(i,j,k,YDIR) .gt. 0) then
							!outward flux has to be positive
							ions_bcvals%top(i,k) = mu_i*ionden(i,j,k)*elecfield(i,j,k,YDIR)*vel_scaling
						else
							ions_bcvals%top(i,k)  = ZERO
						endif
						
					enddo
				enddo
			endif
		endif

		if(g_krank .lt. 0) then

			k = 1
			if(ions_bc_codes(BACK) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						ions_bcvals%back(i,j)=ions_bcvals%back(i,j)/ionscale
					enddo
				enddo
			
			else if(ions_bc_codes(BACK) .eq. 'FLUX') then
				
				do j=1,g_ly
					do i=1,g_lx
						
						do nsp=1,nspecies
							specarray(nsp)=numden(i,j,k,nsp)
						enddo
				
						efield_mag = elecfield(i,j,k,XDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,YDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,ZDIR)**2
						efield_mag = sqrt(efield_mag)
						
						mu_i = getspecmobility(ionspecnum,specarray,efield_mag,&
								electemp(i,j,k),Tg,Pg)

						if(mu_i*elecfield(i,j,k,ZDIR) .lt. 0) then
							!outward flux has to be positive
							ions_bcvals%back(i,j) = -mu_i*ionden(i,j,k)*elecfield(i,j,k,ZDIR)*vel_scaling
						else
							ions_bcvals%back(i,j)  = ZERO
						endif
						
					enddo
				enddo
			endif
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			if(ions_bc_codes(FRONT) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						ions_bcvals%front(i,j)=ions_bcvals%front(i,j)/ionscale
					enddo
				enddo
			
			else if(ions_bc_codes(FRONT) .eq. 'FLUX') then
				
				do j=1,g_ly
					do i=1,g_lx
						
						do nsp=1,nspecies
							specarray(nsp)=numden(i,j,k,nsp)
						enddo
				
						efield_mag = elecfield(i,j,k,XDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,YDIR)**2
						efield_mag = efield_mag + elecfield(i,j,k,ZDIR)**2
						efield_mag = sqrt(efield_mag)
						
						mu_i = getspecmobility(ionspecnum,specarray,efield_mag,&
								electemp(i,j,k),Tg,Pg)

						if(mu_i*elecfield(i,j,k,ZDIR) .gt. 0) then
							!outward flux has to be positive
							ions_bcvals%front(i,j) = mu_i*ionden(i,j,k)*elecfield(i,j,k,ZDIR)*vel_scaling
						else
							ions_bcvals%front(i,j)  = ZERO
						endif
					enddo
				enddo
			endif
		endif


      end subroutine ions_update_bcs
!===========================================================================================
      subroutine ions_update_transport(ionnum,dt,elecfield,electemp,numden,Tg,Pg)
	
		integer, intent(in) :: ionnum
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
		
		real*8, intent(in) :: Tg,Pg

		real*8 :: efield_mag
		real*8 :: mu_i


		real*8  :: x,y,z
	        integer :: i,j,k,nsp
		integer :: ionspecnum

		real*8 :: specarray(nspecies)
		real*8 :: specprod

		vel = ZERO
		dcoeff = ZERO
		reac = ZERO
		source = ZERO

		ionspecnum = ionnum - 1 + ionspecmin

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

					mu_i      =   getspecmobility(ionspecnum,specarray,efield_mag,&
								electemp(i,j,k),Tg,Pg)

					
					vel(i+1,j+1,k+1,XDIR) = mu_i*elecfield(i,j,k,XDIR)
					vel(i+1,j+1,k+1,YDIR) = mu_i*elecfield(i,j,k,YDIR)
					vel(i+1,j+1,k+1,ZDIR) = mu_i*elecfield(i,j,k,ZDIR)
		
					dcoeff(i+1,j+1,k+1)   = getspecdcoeff(ionspecnum,specarray,efield_mag,&
									electemp(i,j,k),Tg,Pg)
					
					
					call getspecproduction(ionspecnum,electemp(i,j,k),Tg,&
								specarray,specprod,efield_mag)

					source(i+1,j+1,k+1)   =   source(i+1,j+1,k+1) + specprod

				enddo
			enddo
		enddo
		
		dcoeff   = dcoeff * dcoeff_scaling
		vel      = vel * vel_scaling
		reac     = reac * reac_scaling
		source   = source * source_scaling
		
		if(timederivflag .eqv. .true.) source = source + ionsoln(:,:,:,ionnum)/dt
		

      end subroutine ions_update_transport
!=========================================================================================
      subroutine ions_solve(time,dt,elecfield,electemp,numden,Tg,Pg,printflag)

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

	integer :: i,nsp
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
	      

	do nsp=1,solved_ions_num

		if(g_myproc .eq. g_rootproc) solve_time=-MPI_Wtime()
	
		if(printflag .eqv. .true.) call printline("solving "//trim(specnames(ionspecmin+nsp-1))//" Density-------------")
		call ions_update_transport(nsp,dt_s,elecfield,electemp,numden,Tg,Pg)
		call ions_update_bcs(nsp,ionsoln(:,:,:,nsp),elecfield,numden,electemp,Tg,Pg)
	
		do i=1,nvcycles

			sterm = 0.d0
		
			call perform_vcycle(ionsoln(:,:,:,nsp),timederivflag,dt_s,&
					g_lx,g_ly,g_lz,vel,dcoeff,reac,source,&
					ions_bc_codes,ions_bcvals,sterm,resnorm)

			if((g_myproc .eq. g_rootproc) .and. (printflag .eqv. .true.)) print *,"it:",i,resnorm

			if(resnorm .le. err_tol) exit
	
		enddo

		if(g_myproc .eq. g_rootproc) solve_time=solve_time+MPI_Wtime()
		if((g_myproc .eq. g_rootproc) .and. (printflag .eqv. .true.)) print *,"solve time:",solve_time

		call exchangehalodata(ionsoln(:,:,:,nsp),g_lx,g_ly,g_lz)
		call update_boundary_ghostvalues(ionsoln(:,:,:,nsp),ions_bc_codes,&
				ions_bcvals,g_lx,g_ly,g_lz)

		if(printflag .eqv. .true.) call printline("------------------------------------------------------------------")
	enddo


      end subroutine ions_solve
!=================================================================
end module ionseq
