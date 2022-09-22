module poissoneq

      use globalvars
      use mgridsteps
      use discretepde
      use chem_module

      implicit none

      real*8, allocatable :: poissonsoln(:,:,:)
      character(LEN=20)   :: poissonscalarname
      
      real*8, allocatable, private    :: dcoeff(:,:,:)
      real*8, allocatable, private    :: vel(:,:,:,:)
      real*8, allocatable, private    :: source(:,:,:)
      real*8, allocatable, private    :: reac(:,:,:)
      character(LEN=4),private        :: poisson_bc_codes(NFACES)
      type(boundarycondition),private :: poisson_bcvals

      real*8, private :: lscale,tscale 
      real*8, private :: potscale

      real*8, private :: vel_scaling
      real*8, private :: dcoeff_scaling
      real*8, private :: reac_scaling
      real*8, private :: source_scaling

      logical, private :: timederivflag
      logical, private :: gradient_limiter_flag
      
      !this is a cheat to do space and time dep.
      !boundary conditions
      real*8, private :: poisson_bcparams(NFACES)

      contains
!========================================================================
	subroutine poisson_initialize(pot_initial,bc_codes,bc_params,&
				length_scale,time_scale,var_scale)

	      real*8, intent(in) :: pot_initial(g_lx+2*g_nglayers,&
			      			g_ly+2*g_nglayers,&
						g_lz+2*g_nglayers)

	      character(LEN=*), intent(in)  :: bc_codes(NFACES)
	      real*8,intent(in) :: bc_params(NFACES)
	      real*8, intent(in) :: length_scale,time_scale
	      real*8, intent(in) :: var_scale
	      
	      allocate(poissonsoln(g_lx+2*g_nglayers,&
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


	      potscale = var_scale
	      lscale   = length_scale
	      tscale   = time_scale

	      poissonsoln = ZERO
	      vel     = ZERO
	      dcoeff  = ZERO
	      reac    = ZERO
	      source  = ZERO

	      !These are factors multiplied (not divided) with
	      !transport coefficients
	      vel_scaling    = tscale/lscale
	      dcoeff_scaling = tscale/(lscale**2)
	      reac_scaling   = tscale
	      source_scaling = tscale/potscale 

	      timederivflag = .false.
	      !gradient_limiter_flag = .true.
	      gradient_limiter_flag = .false.
	
	      call poisson_init(pot_initial)
	      call poisson_set_bcs(bc_codes,bc_params)

	end subroutine poisson_initialize
!=========================================================================
      subroutine poisson_init(pot_initial)
		
	        real*8, intent(in) :: pot_initial(g_lx+2*g_nglayers,&
		    	      			  g_ly+2*g_nglayers,&
						  g_lz+2*g_nglayers)
		real*8  :: x,y,z
	        integer :: i,j,k

		poissonscalarname="Phi"
		poissonscalarname=trim(poissonscalarname)

		poissonsoln=ZERO

		do k=2,g_lz+1
			do j=2,g_ly+1
				do i=2,g_lx+1
					
					x = (g_offx + (i-2)*g_dx + HALF*g_dx)*lscale
					y = (g_offy + (j-2)*g_dy + HALF*g_dy)*lscale
					z = (g_offz + (k-2)*g_dz + HALF*g_dz)*lscale
					
					!g_solnvars(i,j,k,1) = x**2+y**2+z**2
					!poissonsoln(i,j,k) = 0.5*(x**2-x)
					poissonsoln(i,j,k) = 0.d0

				enddo
			enddo
		enddo

		poissonsoln = pot_initial/potscale

      end subroutine poisson_init
!=======================================================================================
      subroutine poisson_set_bcs(bc_codes,bc_params)

		character(LEN=*)  :: bc_codes(NFACES)
		real*8,intent(in) :: bc_params(NFACES)

		integer :: i,j,k
		
	  	call allocate_bcvalues(poisson_bcvals,g_lx,g_ly,g_lz)
		poisson_bc_codes = bc_codes
		poisson_bcparams = bc_params

		if(g_lrank .lt. 0) then
	
			i = 1
			do k=1,g_lz
				do j=1,g_ly
					poisson_bcvals%left(j,k)=bc_params(LEFT)/potscale
				enddo
			enddo
		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			do k=1,g_lz
				do j=1,g_ly
					poisson_bcvals%right(j,k)=bc_params(RIGHT)/potscale
				enddo
			enddo
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			do k=1,g_lz
				do i=1,g_lx
					poisson_bcvals%bottom(i,k)=bc_params(BOTTOM)/potscale
				enddo
			enddo
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			do k=1,g_lz
				do i=1,g_lx
					poisson_bcvals%top(i,k)=bc_params(TOP)/potscale
				enddo
			enddo
		endif

		if(g_krank .lt. 0) then

			k = 1
			do j=1,g_ly
				do i=1,g_lx
					poisson_bcvals%back(i,j)=bc_params(BACK)/potscale
				enddo
			enddo
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			do j=1,g_ly
				do i=1,g_lx
					poisson_bcvals%front(i,j)=bc_params(FRONT)/potscale
				enddo
			enddo
		endif


      end subroutine poisson_set_bcs
!===========================================================================================
      subroutine poisson_update_bcs(time)

		!actual time in seconds
		real*8, intent(in) :: time
		integer :: i,j,k
		real*8 :: x,y,z
		real*8 :: x_cent,y_cent,z_cent,rad
		real*8 :: func
		real*8 :: dist2

		x_cent =  g_prob_specific_params(1)
		y_cent =  g_prob_specific_params(2)
		z_cent =  g_prob_specific_params(3)
		rad    =  g_prob_specific_params(4)
		

		if(g_lrank .lt. 0) then
	
			i = 1
			if(poisson_bc_codes(LEFT) .eq. 'DIRC') then

				do k=1,g_lz
					do j=1,g_ly
						poisson_bcvals%left(j,k)=poisson_bcvals%left(j,k)
					enddo
				enddo
			endif
		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			if(poisson_bc_codes(RIGHT) .eq. 'DIRC') then
	
				do k=1,g_lz
					do j=1,g_ly
						poisson_bcvals%right(j,k)=poisson_bcvals%right(j,k)
					enddo
				enddo
			endif
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			y = (g_offy + (j-1)*g_dy + HALF*g_dy)*lscale
			
			if(poisson_bc_codes(BOTTOM) .eq. 'DIRC') then

				do k=1,g_lz
					do i=1,g_lx
						poisson_bcvals%bottom(i,k)=poisson_bcvals%bottom(i,k)
					enddo
				enddo
			endif
		
			if(poisson_bc_codes(BOTTOM) .eq. 'DZGD') then

				do k=1,g_lz
					z = (g_offz + (k-1)*g_dz + HALF*g_dz)*lscale
					
					do i=1,g_lx

						x = (g_offx + (i-1)*g_dx + HALF*g_dx)*lscale
						dist2 = ((x-x_cent)**2+(z-z_cent)**2)

						if(dist2 .lt. (rad**2)) then

							func = exp(-dist2/rad**2)
							poisson_bcvals%bottom(i,k)=poisson_bcparams(BOTTOM)*func/potscale
							!poisson_bcvals%bottom(i,k)=poisson_bcparams(BOTTOM)/potscale

						else
							poisson_bcvals%bottom(i,k)=JUNKVAL
					 	endif	
					enddo
				enddo

			endif
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			!y = (g_offy + (j-1)*g_dy + HALF*g_dy)*lscale

			if(poisson_bc_codes(TOP) .eq. 'DIRC') then

				do k=1,g_lz
					!z = (g_offz + (k-1)*g_dz + HALF*g_dz)*lscale
					
					do i=1,g_lx
						!x = (g_offx + (i-1)*g_dx + HALF*g_dx)*lscale

						!func = exp(-((x-x_cent)**2+(z-z_cent)**2)/rad**2)
						poisson_bcvals%top(i,k)=poisson_bcvals%top(i,k)
						!poisson_bcvals%top(i,k)=1000.d0*func/potscale

					enddo
				enddo
			endif
		endif

		if(g_krank .lt. 0) then

			k = 1
			if(poisson_bc_codes(BACK) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						poisson_bcvals%back(i,j)=poisson_bcvals%back(i,j)
					enddo
				enddo
			endif
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			if(poisson_bc_codes(FRONT) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						poisson_bcvals%front(i,j)=poisson_bcvals%front(i,j)
					enddo
				enddo
			endif
		endif


      end subroutine poisson_update_bcs
!===========================================================================================
      subroutine poisson_update_transport(dt,elecden,iondens)
	
		real*8, intent(in) :: dt
		real*8  :: x,y,z
	        integer :: i,j,k
	
		real*8, intent(in) :: elecden(-g_nglayers+1:g_lx+g_nglayers,&
					      -g_nglayers+1:g_ly+g_nglayers,&
					      -g_nglayers+1:g_lz+g_nglayers)

		real*8, intent(in) :: iondens(-g_nglayers+1:g_lx+g_nglayers,&
					      -g_nglayers+1:g_ly+g_nglayers,&
					      -g_nglayers+1:g_lz+g_nglayers,solved_ions_num)

		integer :: nsp
		real*8 :: total_ionden
		real*8 :: eden

		vel      = ZERO
		dcoeff   = ZERO
		reac	 = ZERO
		source   = ZERO

		do k=1,g_lz
			do j=1,g_ly
				do i=1,g_lx

					eden         = elecden(i,j,k)
					total_ionden = ZERO

					do nsp=1,solved_ions_num
				       		total_ionden = total_ionden+&
								spec_charge(ionspecmin+nsp-1)*&
								iondens(i,j,k,nsp)
					enddo

					!source   = eden*(1.d0 - total_ionden/eden) 
					source(i+1,j+1,k+1)   = (total_ionden-eden) 
				enddo
			enddo
		enddo

		source = source*ECHARGE/EPS0
		
		dcoeff   = 1.d0 * dcoeff_scaling
		source   = source * source_scaling

		if(timederivflag .eqv. .true.) source = source + poissonsoln/dt
		

      end subroutine poisson_update_transport
!=========================================================================================
      subroutine poisson_solve(time,dt,elecden,iondens,printflag)

	!dt and time in seconds, not scaled
	real*8, intent(in) :: dt
	real*8, intent(in) :: time
	
	logical,intent(in) :: printflag

	real*8, intent(in) :: elecden(-g_nglayers+1:g_lx+g_nglayers,&
				      -g_nglayers+1:g_ly+g_nglayers,&
				      -g_nglayers+1:g_lz+g_nglayers)

	real*8, intent(in) :: iondens(-g_nglayers+1:g_lx+g_nglayers,&
				      -g_nglayers+1:g_ly+g_nglayers,&
				      -g_nglayers+1:g_lz+g_nglayers,solved_ions_num)

	integer :: i
	integer :: nvcycles
	real*8 :: resnorm
	real*8,allocatable:: sterm(:,:,:)
	real*8,allocatable:: res(:,:,:)
	real*8 :: solve_time
	real*8 :: dt_s
	real*8 :: needle_voltage

	real*8 :: err_tol
	
	err_tol=ERRTOLR
	nvcycles=20

	dt_s = dt/tscale
	
	allocate(sterm(g_lx,g_ly,g_lz))
	allocate(res(g_lx,g_ly,g_lz))
	      
	if(g_myproc .eq. g_rootproc) solve_time=-MPI_Wtime()

	call poisson_update_transport(dt_s,elecden,iondens)
	call poisson_update_bcs(time)
	
	do i=1,nvcycles

		sterm = 0.d0
		
		call perform_vcycle(poissonsoln,timederivflag,dt_s,g_lx,g_ly,g_lz,vel,dcoeff,reac,source,&
				poisson_bc_codes,poisson_bcvals,sterm,resnorm)
	
		!Sanity check; see if AX and b are found correctly=======================
		
		!call perform_gjacobi_smoothing(poissonsoln,dt,g_lx,g_ly,g_lz,.true.,sterm,res,&
		!			vel,dcoeff,reac,source,&
		!			poisson_bc_codes,poisson_bcvals,&
		!			g_lrank,g_rrank,g_brank,g_trank,g_krank,g_frank,&
		!			g_llenx,g_lleny,g_llenz,1)
		
		!call compute_norm(res,g_lx,g_ly,g_lz,resnorm)
		
		!Sanity check; see if AX and b are found correctly=======================

		if((printflag .eqv. .true.) .and. (g_myproc .eq. g_rootproc)) print *,"it:",i,resnorm

		if(resnorm .le. err_tol) exit
	
	enddo

	if(g_myproc .eq. g_rootproc) solve_time=solve_time+MPI_Wtime()
	if((g_myproc .eq. g_rootproc) .and. (printflag .eqv. .true.))  print *,"solve time:",solve_time

	call exchangehalodata(poissonsoln,g_lx,g_ly,g_lz)
	call update_boundary_ghostvalues(poissonsoln,poisson_bc_codes,&
				poisson_bcvals,g_lx,g_ly,g_lz)
	!correct potential with charge simulation technique (Positive and
	!negative streamers in ambient air: modelling evolution and velocities
	!A Luque, V Ratushnaya, U Ebert - J Phys D:Appl Phys, 2008 
        if(g_prob_specific_params(5) .ne. ZERO) then
		
		if(g_prob_specific_params(12) .ne. ZERO) then
			needle_voltage = timedependentvoltage(time)
		else
			needle_voltage = g_prob_specific_params(5)
		endif

		if((printflag .eqv. .true.) .and. (g_myproc .eq. g_rootproc)) print *,"needle_voltage:",needle_voltage
		
		call correctpotential(poissonsoln,needle_voltage,g_prob_specific_params(6),&
				g_prob_specific_params(7),g_prob_specific_params(8),g_prob_specific_params(9),&
				g_prob_specific_params(10),g_prob_specific_params(11))
	endif
	

      end subroutine poisson_solve
!=================================================================
      subroutine correctpotential(potential,phi_needle,needlex,needley,needlez,px,py,pz)
		
		real*8, intent(in) :: phi_needle,needlex,needley,needlez
		real*8, intent(in) :: px,py,pz
		real*8, intent(inout) :: potential(-g_nglayers+1:g_lx+g_nglayers,&
						   -g_nglayers+1:g_ly+g_nglayers,&
						   -g_nglayers+1:g_lz+g_nglayers)

		real*8 :: q

		integer :: enclosing_proc
		integer :: proc_x,proc_y,proc_z

		real*8 :: pot_value

		real*8 :: r_needle,r
		integer :: i,j,k,ierr
		real*8 :: x,y,z
		
		

		proc_x = floor((px/lscale - g_xorigin)/g_llenx)+1
		proc_y = floor((py/lscale - g_yorigin)/g_lleny)+1
		proc_z = floor((pz/lscale - g_zorigin)/g_llenz)+1

		call find1dindex(proc_x,proc_y,proc_z,g_px,g_py,g_pz,enclosing_proc)

		enclosing_proc = enclosing_proc-1
		if((enclosing_proc .lt. 0) .or. (enclosing_proc .gt. (g_nprocs-1)) )  then
			call printline("enclosing processor not found in correct potential routine")
			call MPI_Abort(g_comm,666,ierr)
		endif
		
		!z = (g_offz + (k-1)*g_dz + HALF*g_dz)*lscale

		if(g_myproc .eq. enclosing_proc) then
			
			i = (px/lscale - g_offx)/g_dx + 1
			j = (py/lscale - g_offy)/g_dy + 1
			k = (pz/lscale - g_offz)/g_dz + 1
			
			pot_value = potential(i,j,k)*potscale

			r_needle =sqrt((needlex-px)**2+(needley-py)**2+(needlez-pz)**2)
			
			q = (phi_needle-pot_value)*r_needle

		endif
		call MPI_Bcast(q,1,MPI_REAL8,enclosing_proc,g_comm,ierr)
		call MPI_Barrier(g_comm,ierr)
		!print *,"q=",q
		do k=1,g_lz
			z = (g_offz + (k-1)*g_dz + HALF*g_dz)*lscale
			
			do j=1,g_ly
				y = (g_offy + (j-1)*g_dy + HALF*g_dy)*lscale

				do i=1,g_lx
					x = (g_offx + (i-1)*g_dx + HALF*g_dx)*lscale
				
					r = sqrt((needlex-x)**2+(needley-y)**2+(needlez-z)**2)	
					potential(i,j,k)=potential(i,j,k)+(q/r)/potscale

			        enddo
			enddo
		enddo

		
      end subroutine correctpotential
!=================================================================
      subroutine minmodlimiter(ntr,dtr,phi)
		
		real*8,intent(in) :: ntr,dtr
		real*8,intent(out) :: phi

		real*8 :: r

		if(dtr .ne. 0) then
			r   = ntr/dtr
			phi = max(ZERO,min(1.d0,r))	
		else
			phi = 1.d0
		endif

      end subroutine minmodlimiter
!=================================================================
      subroutine update_efield(efield)

		real*8, intent(out) :: efield(-g_nglayers+1:g_lx+g_nglayers,&
					      -g_nglayers+1:g_ly+g_nglayers,&
					      -g_nglayers+1:g_lz+g_nglayers, NDIM)
      
		character(LEN=4)                :: efield_bc_codes(NFACES)
      		type(boundarycondition)         :: efield_bcvals

		integer :: i,j,k
		integer :: soln_i,soln_j,soln_k
		real*8  :: phi,ntr,dtr
		real*8  :: val_gt_one

		val_gt_one = TWO 
		!Electric field at physical boundaries are copied 
		!to ghost cell
		efield_bc_codes = 'ZGRD'
	  	call allocate_bcvalues(efield_bcvals,g_lx,g_ly,g_lz)
		

		call exchangehalodata(poissonsoln,g_lx,g_ly,g_lz)
		call update_boundary_ghostvalues(poissonsoln,poisson_bc_codes,&
				poisson_bcvals,g_lx,g_ly,g_lz)

		do k=1,g_lz
			soln_k=k+1
			
			do j=1,g_ly
				soln_j = j+1
			
				do i=1,g_lx
					soln_i = i+1

					efield(i,j,k,XDIR) = 0.5*(poissonsoln(soln_i+1,soln_j,soln_k) - &
							     poissonsoln(soln_i-1,soln_j,soln_k))/g_dx
					
					efield(i,j,k,YDIR) = 0.5*(poissonsoln(soln_i,soln_j+1,soln_k) - &
							     poissonsoln(soln_i,soln_j-1,soln_k))/g_dy
					
				        efield(i,j,k,ZDIR) = 0.5*(poissonsoln(soln_i,soln_j,soln_k+1) - &
							     poissonsoln(soln_i,soln_j,soln_k-1))/g_dz

				enddo
			enddo
		enddo


		if(gradient_limiter_flag .eqv. .true.) then
			do k=1,g_lz
				soln_k=k+1
			
				do j=1,g_ly
					soln_j = j+1
			
					do i=1,g_lx
						soln_i = i+1

						!X_direction============================================
						ntr = (poissonsoln(soln_i,soln_j,soln_k) - &
								     poissonsoln(soln_i-1,soln_j,soln_k))

						dtr = (poissonsoln(soln_i+1,soln_j,soln_k) - &
								poissonsoln(soln_i,soln_j,soln_k))
				
						call minmodlimiter(ntr,dtr,phi)	

						efield(i,j,k,XDIR) = phi*efield(i,j,k,XDIR)			
						!=======================================================	
					
						!Y_direction============================================
						ntr = (poissonsoln(soln_i,soln_j,soln_k) - &
								     poissonsoln(soln_i,soln_j-1,soln_k))

						dtr = (poissonsoln(soln_i,soln_j+1,soln_k) - &
								poissonsoln(soln_i,soln_j,soln_k))
				
						call minmodlimiter(ntr,dtr,phi)	

						efield(i,j,k,YDIR) = phi*efield(i,j,k,YDIR)				
						!=======================================================	
						
						!Z_direction============================================
						ntr = (poissonsoln(soln_i,soln_j,soln_k) - &
								     poissonsoln(soln_i,soln_j,soln_k-1))

						dtr = (poissonsoln(soln_i,soln_j,soln_k+1) - &
								poissonsoln(soln_i,soln_j,soln_k))
				
						call minmodlimiter(ntr,dtr,phi)	

						efield(i,j,k,ZDIR) = phi*efield(i,j,k,ZDIR)				
						!=======================================================	
					enddo
				enddo
			enddo
		endif

		! scaling to dimensional value (V/m)
		! efield is negative gradient of potential
		efield = -efield * (potscale/lscale)
		
		call exchangehalodata(efield(:,:,:,XDIR),g_lx,g_ly,g_lz)
		call exchangehalodata(efield(:,:,:,YDIR),g_lx,g_ly,g_lz)
		call exchangehalodata(efield(:,:,:,ZDIR),g_lx,g_ly,g_lz)
		
		call update_boundary_ghostvalues(efield(:,:,:,XDIR),efield_bc_codes,&
				efield_bcvals,g_lx,g_ly,g_lz)

		call update_boundary_ghostvalues(efield(:,:,:,YDIR),efield_bc_codes,&
				efield_bcvals,g_lx,g_ly,g_lz)
		
		call update_boundary_ghostvalues(efield(:,:,:,ZDIR),efield_bc_codes,&
				efield_bcvals,g_lx,g_ly,g_lz)

      end subroutine update_efield
!=================================================================
      function timedependentvoltage(time) result(voltage)

	real*8, intent(in) :: time !in seconds
	real*8 :: voltage

	real*8 :: rise_time,flat_time,fall_time
	real*8 :: V0

	real*8 :: t1,t2,t3

	rise_time = g_prob_specific_params(13)
	flat_time = g_prob_specific_params(14)
	fall_time = g_prob_specific_params(15)

	t1 = rise_time
	t2 = rise_time + flat_time
	t3 = rise_time + flat_time + fall_time

	V0 = g_prob_specific_params(5)

	if(time .lt. t1) then
		voltage = V0*time/rise_time
	else if((time .ge. t1) .and. (time .lt. t2)) then
		voltage = V0
	else if((time .ge. t2) .and. (time .lt. t3)) then
		voltage = V0*(1.d0 - (time-t2)/fall_time)
	else
		voltage = ZERO
	endif
		 
      end function timedependentvoltage
!=================================================================
end module poissoneq
