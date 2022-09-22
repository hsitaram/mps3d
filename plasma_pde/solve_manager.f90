module solve_manager

      use chem_module	
      use globalvars
      use cartwrite
      use poissoneq
      use edeneq
      use ionseq
      use neutralseq
      use eenergyeq
      use par_decompose
      use inputs
      use mgridsteps
      implicit none

      !public variables
      integer, private :: eden_eqnum
      integer, private :: ions_eqnum
      integer, private :: neutrals_eqnum
      integer, private :: pot_eqnum
      integer, private :: etemp_eqnum

      integer, private :: edensoln_num
      integer, private :: potsoln_num
      integer, private :: elecfield_num
      integer, private :: etempsoln_num
      integer, private :: eenergysoln_num
      integer, private :: eprod_num
      
      integer, private :: num_eq
      integer, private :: allsolnvars_num

      real*8, private  :: timestep
      real*8, private  :: finaltime
      integer, private :: maxtimesteps
      integer, private :: outputstep
      integer, private :: stdoutstep

      real*8, private  :: bgpressure
      real*8, private  :: bgtemperature

      real*8, private  :: lengthscale
      real*8, private  :: timescale
      real*8, private  :: elecdenscale
      real*8, private  :: iondenscale
      real*8, private  :: specdenscale
      real*8, private  :: potscale
      real*8, private  :: etempscale
      real*8, private  :: eenergyscale

      real*8, private  :: sec_elec_emission

      real*8, private,allocatable :: spec_init(:)
      real*8, private :: pot_init
      real*8, private :: electemp_init

      
      real*8, allocatable,private :: all_solnvars(:,:,:,:)
      character(LEN=MAXSTRSIZE),allocatable,private :: all_scalarnames(:)
      character(LEN=4),allocatable,private :: solver_bc_codes(:,:)
      real*8, allocatable, private :: solver_bcparams(:,:)
      real*8, allocatable, private :: avg_neutralprod(:,:,:,:)
      logical, private :: printproductionflag;

    	contains
!===================================================================================
	subroutine readsolverinputfile(solverfilename)

	      character(LEN=*), intent(in) :: solverfilename

	      character(LEN=MAXSTRSIZE) :: temp
	      character(LEN=MAXSTRSIZE) :: eden_bc_code
	      character(LEN=MAXSTRSIZE) :: spec_bc_code
	      character(LEN=MAXSTRSIZE) :: pot_bc_code
	      character(LEN=MAXSTRSIZE) :: etemp_bc_code
	      
	      real*8 :: eden_bcval,spec_bcval,pot_bcval,etemp_bcval
	      integer :: fnum,i,nsp
	      
	      fnum=17

	      open(unit=fnum,file=solverfilename)
		
	      !TIME INFO===========================================
	      read(fnum,*) temp
	      read(fnum,*) temp
	      read(fnum,*) temp
	      
	      read(fnum,*) temp,timestep
	      read(fnum,*) temp,finaltime
      	      read(fnum,*) temp,maxtimesteps
      	      read(fnum,*) temp,outputstep
      	      read(fnum,*) temp,stdoutstep
	      !====================================================
	      !GAS PRESSURE AND TEMPERTURE=========================
	      read(fnum,*) temp
	      read(fnum,*) temp
	      read(fnum,*) temp
	      
	      read(fnum,*) temp,bgpressure
	      read(fnum,*) temp,bgtemperature
	      !====================================================
	      
	      !SCALES===============================================
	      read(fnum,*) temp
	      read(fnum,*) temp
	      read(fnum,*) temp
	    
	      read(fnum,*) temp,lengthscale
	      read(fnum,*) temp,timescale 
	      read(fnum,*) temp,elecdenscale
	      read(fnum,*) temp,iondenscale
	      read(fnum,*) temp,specdenscale
	      read(fnum,*) temp,potscale
	      read(fnum,*) temp,etempscale
	      !=====================================================
	      
	      !INITIAL CONDITIONS==================================
	      read(fnum,*) temp
	      read(fnum,*) temp
	      read(fnum,*) temp
	      
	      !number densities
	      do i=1,nspecies
			read(fnum,*) temp,spec_init(i)
	      enddo

	      !potential and electron temperature
	      read(fnum,*) temp,pot_init
	      read(fnum,*),temp,electemp_init
	      !=====================================================


	      !BOUNDARY CONDITIONS==================================
	      read(fnum,*) temp
	      read(fnum,*) temp
	      read(fnum,*) temp
	      
	      do i=1,NFACES
	      	
	      		read(fnum,*) temp,eden_bc_code,spec_bc_code,&
					pot_bc_code,etemp_bc_code,&
				  	eden_bcval,spec_bcval,&
					pot_bcval,etemp_bcval

			!print *,eden_bcval,spec_bcval,pot_bcval,etemp_bcval

		  	solver_bc_codes(i,eden_eqnum) = eden_bc_code
		 	solver_bcparams(i,eden_eqnum) = eden_bcval

		  	solver_bc_codes(i,ions_eqnum) = spec_bc_code
			solver_bcparams(i,ions_eqnum) = spec_bcval
		  	
			solver_bc_codes(i,neutrals_eqnum) = spec_bc_code
			solver_bcparams(i,neutrals_eqnum) = spec_bcval

			solver_bc_codes(i,pot_eqnum)  =  pot_bc_code
		 	solver_bcparams(i,pot_eqnum)  = pot_bcval

			solver_bc_codes(i,etemp_eqnum) = etemp_bc_code
			solver_bcparams(i,etemp_eqnum) = etemp_bcval
	      enddo
	      !====================================================

	      !PROBLEM SPECIFIC PARAMETERS=========================
	      read(fnum,*) temp
	      read(fnum,*) temp
	      read(fnum,*) temp

	      read(fnum,*) g_prob_specific_params(:)
	      !====================================================

	      close(fnum)

	end subroutine readsolverinputfile
!===================================================================================
	subroutine scalelengths()

		g_xlen = g_xlen/lengthscale
		g_ylen = g_ylen/lengthscale
		g_zlen = g_zlen/lengthscale

		g_xorigin = g_xorigin/lengthscale
		g_yorigin = g_yorigin/lengthscale
		g_zorigin = g_zorigin/lengthscale

		g_dx = g_dx/lengthscale
		g_dy = g_dy/lengthscale
		g_dz = g_dz/lengthscale

		g_offx = g_offx/lengthscale
		g_offy = g_offy/lengthscale
		g_offz = g_offz/lengthscale

		g_llenx = g_llenx/lengthscale
		g_lleny = g_lleny/lengthscale
		g_llenz = g_llenz/lengthscale

	end subroutine scalelengths
!===================================================================================
        subroutine solversetup()

		integer :: i
		real*8,allocatable  :: eenergy_init(:,:,:)

		call initializechemistry()
		
		! all species, potential and electron temperature
		num_eq = 1 !for electrons
		num_eq = num_eq + 1 !for ions
		num_eq = num_eq + 1 !for neutrals
		num_eq = num_eq + 1 + 1 !potential and electron temperature

		allocate(solver_bc_codes(NFACES,num_eq))
		allocate(solver_bcparams(NFACES,num_eq))
		allocate(spec_init(nspecies))

	        !eq nums are for boundary conditions and
		!management here	
		eden_eqnum     = 1
		ions_eqnum     = eden_eqnum + 1
		neutrals_eqnum = ions_eqnum + 1
		pot_eqnum      = neutrals_eqnum + 1
		etemp_eqnum    = pot_eqnum  + 1

		!soln nums are for passing solutions between
		!equations
		allsolnvars_num = nspecies + 6
		edensoln_num    = especnum
	        potsoln_num     = nspecies + 1
		elecfield_num   = nspecies + 2
		etempsoln_num   = nspecies + 5
		eenergysoln_num = nspecies + 6
		eprod_num       = nspecies + 7

		!track electron and neutral species production
		allsolnvars_num = allsolnvars_num + 1
		allsolnvars_num = allsolnvars_num + solved_neutrals_num 

		!Read inputs====================================================
		call readgeominputfile('GEOM_INPUTS')
		call readsolverinputfile('PLASMA_INPUTS')
		!===============================================================
		
		!Parallel decomposition=========================================
		call domaindecompose()
		!===============================================================

		!===============================================================
		!scale lengths
		call scalelengths()
		!===============================================================

		!===============================================================
		!all solnvars include all species,potential, electric field
		!and electron temperature
		allocate(all_solnvars(g_lx+2*g_nglayers,&
				      g_ly+2*g_nglayers,&
				      g_lz+2*g_nglayers,allsolnvars_num))

		allocate(all_scalarnames(allsolnvars_num))
		!===============================================================

		!For H2 O2 ignition runs========================================
		if(g_prob_specific_params(16) .gt. 0) then
			printproductionflag=.true.
		else
			printproductionflag=.false.
		endif

		if(solved_neutrals_num .gt. 0 .and. printproductionflag .eqv. .true.) then
			allocate(avg_neutralprod(g_lx,g_ly,g_lz,solved_neutrals_num))
			avg_neutralprod = ZERO
		endif
		!===============================================================

		!Setting names for variables====================================
		all_scalarnames ='DEFAULT'
		do i=1,nspecies
			all_scalarnames(i) = trim(specnames(i))//'_density (#/m3)'
		enddo

		all_scalarnames(potsoln_num)    = 'Potential (V)'

		all_scalarnames(elecfield_num-1+XDIR) = 'Electric_field_X (V/m)'
		all_scalarnames(elecfield_num-1+YDIR) = 'Electric_field_Y (V/m)'
		all_scalarnames(elecfield_num-1+ZDIR) = 'Electric_field_Z (V/m)'
		
		all_scalarnames(etempsoln_num)    = 'Electron_Temp (K)'
		all_scalarnames(eenergysoln_num)  = 'Electron_energy (J/m3)'
		all_scalarnames(eprod_num)        = 'E_prod (#/m3/s)'

		do i=1,solved_neutrals_num
			all_scalarnames(eprod_num+i) = trim(specnames(neutralspecmin+i-1))&
								//'_prod (#/m3/s)'
		enddo
		!================================================================

		!setting initial conditions======================================
		all_solnvars    = ZERO
		do i=1,nspecies
			all_solnvars(:,:,:,i) = spec_init(i)
		enddo
	
		all_solnvars(:,:,:,potsoln_num)   = pot_init 

		all_solnvars(:,:,:,elecfield_num-1+XDIR) = ZERO
		all_solnvars(:,:,:,elecfield_num-1+YDIR) = ZERO
		all_solnvars(:,:,:,elecfield_num-1+ZDIR) = ZERO
		
		all_solnvars(:,:,:,etempsoln_num) = electemp_init
	        allocate(eenergy_init(g_lx+2*g_nglayers,&
				    g_ly+2*g_nglayers,&
				    g_lz+2*g_nglayers))

	         eenergyscale = elecdenscale*k_BMAN*etempscale
	         eenergy_init = 1.5*all_solnvars(:,:,:,edensoln_num)*&
			      k_BMAN*all_solnvars(:,:,:,etempsoln_num)

	        all_solnvars(:,:,:,eenergysoln_num) = eenergy_init
	       !================================================================	
		
	       call poisson_initialize(all_solnvars(:,:,:,potsoln_num),&
			               solver_bc_codes(:,pot_eqnum),&
				       solver_bcparams(:,pot_eqnum),&
				       lengthscale,timescale,potscale)
	       
	       call eden_initialize(all_solnvars(:,:,:,edensoln_num),&
			               solver_bc_codes(:,eden_eqnum),&
				       solver_bcparams(:,eden_eqnum),&
				       lengthscale,timescale,elecdenscale)
	       
	       call ions_initialize(all_solnvars(:,:,:,ionspecmin:ionspecmax),&
			               solver_bc_codes(:,ions_eqnum),&
				       solver_bcparams(:,ions_eqnum),&
				       lengthscale,timescale,iondenscale)

	      if(solved_neutrals_num .gt. 0) then

		      call neutrals_initialize(all_solnvars(:,:,:,neutralspecmin:neutralspecmax),&
				      		solver_bc_codes(:,neutrals_eqnum),&
						solver_bcparams(:,neutrals_eqnum),&
						lengthscale,timescale,specdenscale)
	      endif


	      call eenergy_initialize(eenergy_init,&
			      	      solver_bc_codes(:,etemp_eqnum),&
				      solver_bcparams(:,etemp_eqnum),&
				      lengthscale,timescale,eenergyscale)

	end subroutine solversetup
!===================================================================================
	subroutine writeoutputfile(fname)
	
		character(LEN=*), intent(in) :: fname

		call writeoutput(fname,all_solnvars,allsolnvars_num,all_scalarnames,&
				g_lx,g_ly,g_lz)

	end subroutine writeoutputfile
!===================================================================================
	subroutine timestepping()

		real*8 :: t
		integer :: timestepnum
		integer :: outputnum
		character(LEN=7) :: outputnumstr
		integer :: nsp

		logical :: printflag

		t=ZERO
		timestepnum=0
		call writeoutputfile('output_initial')

		do while((timestepnum .lt. maxtimesteps) .and. (t .lt. finaltime))

			!---------------------------------------------------------------------------------------------
			printflag = .false.
			if(mod(timestepnum,stdoutstep) .eq. 0) printflag=.true.
			
			if(mod(timestepnum,outputstep) .eq. 0) then

				outputnum=timestepnum/outputstep
				write(outputnumstr,'(I7.7)') outputnum
				call writeoutputfile('output_'//outputnumstr)
				
				if(printproductionflag) then

					avg_neutralprod = avg_neutralprod + &
						all_solnvars(:,:,:,eprod_num+1:eprod_num+solved_neutrals_num)

				endif

			endif
			!---------------------------------------------------------------------------------------------


			!--------------------------------------------------------------------------------------------			
			if(printflag .eqv. .true.) then
				call printline("Potential solve==================================================================")
			endif

			call poisson_solve(t,timestep,all_solnvars(:,:,:,especnum),&
					all_solnvars(:,:,:,ionspecmin:ionspecmax),printflag)
			
			all_solnvars(:,:,:,potsoln_num) = poissonsoln*potscale
			call update_efield(all_solnvars(:,:,:,elecfield_num:elecfield_num+NDIM-1))
			
			if(printflag .eqv. .true.) then
				call printline("=================================================================================")
			endif
			!--------------------------------------------------------------------------------------------			
			
			!--------------------------------------------------------------------------------------------			
			if(printflag .eqv. .true.) then
				call printline("Electron energy solve==================================================================")
			endif

			call eenergy_solve(t,timestep,all_solnvars(:,:,:,potsoln_num),all_solnvars(:,:,:,elecfield_num:elecfield_num+NDIM-1),&
					   all_solnvars(:,:,:,etempsoln_num),all_solnvars(:,:,:,1:nspecies),&
				  	   bgtemperature,bgpressure,printflag)
			
			call update_electrontemp(all_solnvars(:,:,:,etempsoln_num),all_solnvars(:,:,:,especnum))
			all_solnvars(:,:,:,eenergysoln_num) = eenergysoln*eenergyscale
			
			if(printflag .eqv. .true.) then
				call printline("=================================================================================")
			endif
			!--------------------------------------------------------------------------------------------			



			!--------------------------------------------------------------------------------------------			
			if(printflag .eqv. .true.) then
				 call printline("Electron solve===================================================================")
			endif
			
			call eden_solve(t,timestep,all_solnvars(:,:,:,elecfield_num:elecfield_num+NDIM-1),&
					all_solnvars(:,:,:,etempsoln_num),all_solnvars(:,:,:,1:nspecies),&
					all_solnvars(:,:,:,eprod_num),bgtemperature,bgpressure,printflag)
			all_solnvars(:,:,:,edensoln_num) = edensoln*elecdenscale
			
			if(printflag .eqv. .true.) then
				call printline("=================================================================================")
			endif
			!--------------------------------------------------------------------------------------------			


			!--------------------------------------------------------------------------------------------			
			if(printflag .eqv. .true.) then
				call printline("Ions solve=======================================================================")
			endif
			call ions_solve(t,timestep,all_solnvars(:,:,:,elecfield_num:elecfield_num+NDIM-1),&
					all_solnvars(:,:,:,etempsoln_num),all_solnvars(:,:,:,1:nspecies),&
					bgtemperature,bgpressure,printflag)
			all_solnvars(:,:,:,ionspecmin:ionspecmax) = ionsoln*iondenscale
			
			if(printflag .eqv. .true.) then
				 call printline("=================================================================================")
			endif
			!--------------------------------------------------------------------------------------------			


			!--------------------------------------------------------------------------------------------			
			if(solved_neutrals_num .gt. 0) then
				if(printflag .eqv. .true.) then
					 call printline("neutral solve===========================================================")
				endif
				
				call neutrals_solve(t,timestep,all_solnvars(:,:,:,elecfield_num:elecfield_num+NDIM-1),&
						    all_solnvars(:,:,:,etempsoln_num),all_solnvars(:,:,:,1:nspecies),&
						    all_solnvars(:,:,:,eprod_num+1:eprod_num+solved_neutrals_num),&
						    bgtemperature,bgpressure,printflag)

				all_solnvars(:,:,:,neutralspecmin:neutralspecmax)=neutralsoln*specdenscale
				
				if(printflag .eqv. .true.) then
					call printline("========================================================================")
				endif
			endif
			!--------------------------------------------------------------------------------------------			
			

			t=t+timestep
			timestepnum=timestepnum+1
			if((printflag .eqv. .true.) .and. (g_myproc .eq. g_rootproc)) then
				print *,"******************time:",t,"***********************************************"
			endif

		enddo

		avg_neutralprod = avg_neutralprod*outputstep/timestepnum

		if(printproductionflag .eqv. .true.) then
			call write_avg_neutralprod()
		endif

		call writeoutputfile('output_final')

	end subroutine timestepping
!===================================================================================
	subroutine write_avg_neutralprod()

		integer :: nrestricts
		type(flexarray3d),allocatable :: restrictedarrays(:,:)
		integer :: i,nsp
		integer :: l,m,n
		integer :: finalextents(NDIM),finalextents_global(NDIM)
		real*8 :: boxlo(NDIM),boxhi(NDIM)

		real*8, allocatable :: globalrestrictedarray(:,:,:,:)

		nrestricts = g_prob_specific_params(16)

		!flex arrays for each restricted level
		allocate(restrictedarrays(nrestricts+1,solved_neutrals_num))

		!global array for printing
		allocate(globalrestrictedarray(g_px*g_lx/(2**nrestricts),&
					       g_py*g_ly/(2**nrestricts),&
					       g_pz*g_lz/(2**nrestricts),&
					       solved_neutrals_num))
		
		!last restricted level
		finalextents(XDIR) = g_lx/(2**nrestricts)
		finalextents(YDIR) = g_ly/(2**nrestricts)
		finalextents(ZDIR) = g_lz/(2**nrestricts)

		!global number of cells
		finalextents_global(XDIR) = finalextents(XDIR)*g_px
		finalextents_global(YDIR) = finalextents(YDIR)*g_py
		finalextents_global(ZDIR) = finalextents(ZDIR)*g_pz

		!corners of the box
		boxlo(XDIR) = g_xorigin*lengthscale 
		boxlo(YDIR) = g_yorigin*lengthscale 
		boxlo(ZDIR) = g_zorigin*lengthscale 
		
		boxhi(XDIR) = (g_xorigin + g_xlen)*lengthscale 
		boxhi(YDIR) = (g_yorigin + g_ylen)*lengthscale 
		boxhi(ZDIR) = (g_zorigin + g_zlen)*lengthscale 

		!allocate the arrays in restricted levels
		do nsp=1,solved_neutrals_num

			do i=1,nrestricts+1
	
				l = g_lx/(2**(i-1))
				m = g_ly/(2**(i-1))
				n = g_lz/(2**(i-1))

				allocate(restrictedarrays(i,nsp)%arr(l,m,n))
			enddo

		enddo

		!initialize the top most level
		do nsp=1,solved_neutrals_num
			restrictedarrays(1,nsp)%arr = avg_neutralprod(:,:,:,nsp)
		enddo

		!perform restrictions
		do nsp=1,solved_neutrals_num
			
			do i=1,nrestricts

				l = g_lx/(2**(i-1))
				m = g_ly/(2**(i-1))
				n = g_lz/(2**(i-1))

				call restriction(restrictedarrays(i,nsp)%arr,l,m,n,restrictedarrays(i+1,nsp)%arr)
				
			enddo
		enddo

				
		do nsp=1,solved_neutrals_num
			
			call gatherterms(restrictedarrays(nrestricts+1,nsp)%arr,finalextents(XDIR),&
					finalextents(YDIR),finalextents(ZDIR),globalrestrictedarray(:,:,:,nsp))

		enddo

		if(g_myproc .eq. g_rootproc) then	
			call writeprodfile(globalrestrictedarray,finalextents_global,boxlo,boxhi)
		endif

		deallocate(restrictedarrays)
		deallocate(globalrestrictedarray)
				

	end subroutine write_avg_neutralprod
!===================================================================================
subroutine writeprodfile(prod,np,boxlo,boxhi)
	
	integer, intent(in) :: np(NDIM)
	real*8, intent(in)  :: boxlo(NDIM),boxhi(NDIM)
	real*8, intent(in)  :: prod(np(XDIR),np(YDIR),np(ZDIR),solved_neutrals_num)

	integer :: fnum,nsp

	fnum=13
	
	open(unit=fnum,file='prodterms.ho')

	write(fnum, *) boxlo(XDIR),boxhi(XDIR)
	write(fnum, *) boxlo(YDIR),boxhi(YDIR)
	write(fnum, *) boxlo(ZDIR),boxhi(ZDIR)

	write(fnum,*) np(XDIR),np(YDIR),np(ZDIR)

	do nsp=1,solved_neutrals_num
		write(fnum,*) prod(:,:,:,nsp)
	enddo

	close(fnum)	

end subroutine writeprodfile
!===================================================================================

end module solve_manager
