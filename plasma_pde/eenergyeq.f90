module eenergyeq

      use globalvars
      use mgridsteps
      use chem_module

      implicit none

      real*8, allocatable :: eenergysoln(:,:,:)
      character(LEN=20)   :: eenergyname
      
      real*8, private                 :: etempmin
      real*8, allocatable, private    :: dcoeff(:,:,:)
      real*8, allocatable, private    :: vel(:,:,:,:)
      real*8, allocatable, private    :: source(:,:,:)
      real*8, allocatable, private    :: reac(:,:,:)
      character(LEN=4),private        :: eenergy_bc_codes(NFACES)
      type(boundarycondition),private :: eenergy_bcvals

      real*8, private :: lscale,tscale 
      real*8, private :: eenergyscale

      real*8, private :: vel_scaling
      real*8, private :: dcoeff_scaling
      real*8, private :: reac_scaling
      real*8, private :: source_scaling

      logical, private :: timederivflag

      !this is a cheat, the solve manager only gets
      !electron temperature bcs
      !we need current electron density to compute energy
      real*8, private :: eenergy_bcparams(NFACES)
      

      contains
!========================================================================
	subroutine eenergy_initialize(eenergy_initial,bc_codes,bc_params,&
				length_scale,time_scale,var_scale)

	      real*8, intent(in) :: eenergy_initial(g_lx+2*g_nglayers,&
			      		    	    g_ly+2*g_nglayers,&
						    g_lz+2*g_nglayers)
	      
	      character(LEN=*), intent(in)  :: bc_codes(NFACES)
	      real*8,intent(in) :: bc_params(NFACES)
	      real*8, intent(in) :: length_scale,time_scale
	      real*8, intent(in) :: var_scale
	      
	      allocate(eenergysoln(g_lx+2*g_nglayers,&
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

	      etempmin = 0.1*eVtoK

	      eenergyscale = var_scale
	      lscale   = length_scale
	      tscale   = time_scale

	      eenergysoln = ZERO
	      vel     = ZERO
	      dcoeff  = ZERO
	      reac    = ZERO
	      source  = ZERO

	      !These are factors multiplied (not divided) with
	      !transport coefficients
	      vel_scaling    = tscale/lscale
	      dcoeff_scaling = tscale/(lscale**2)
	      reac_scaling   = tscale
	      source_scaling = tscale/eenergyscale 

	      timederivflag = .true.
	
	      call eenergy_init(eenergy_initial)
	      call eenergy_set_bcs(bc_codes,bc_params)

	end subroutine eenergy_initialize
!=========================================================================
      subroutine eenergy_init(eenergy_initial)
		
	        real*8, intent(in) :: eenergy_initial(g_lx+2*g_nglayers,&
			      		    	    g_ly+2*g_nglayers,&
						    g_lz+2*g_nglayers)
		real*8  :: x,y,z
	        integer :: i,j,k

		eenergyname="Electron_energy (#/m3)"
		eenergyname=trim(eenergyname)

		eenergysoln=ZERO

		do k=2,g_lz+1
			do j=2,g_ly+1
				do i=2,g_lx+1
					
					x = (g_offx + (i-2)*g_dx + HALF*g_dx)*lscale
					y = (g_offy + (j-2)*g_dy + HALF*g_dy)*lscale
					z = (g_offz + (k-2)*g_dz + HALF*g_dz)*lscale
					
					!g_solnvars(i,j,k,1) = x**2+y**2+z**2
					!eenergysoln(i,j,k) = 0.5*(x**2-x)
					eenergysoln(i,j,k) = 0.d0

				enddo
			enddo
		enddo

		eenergysoln = eenergy_initial/eenergyscale

      end subroutine eenergy_init
!=======================================================================================
      subroutine eenergy_set_bcs(bc_codes,bc_params)

		character(LEN=*)  :: bc_codes(NFACES)
		real*8,intent(in) :: bc_params(NFACES)

		integer :: i,j,k

	  	call allocate_bcvalues(eenergy_bcvals,g_lx,g_ly,g_lz)
		eenergy_bc_codes = bc_codes

		!copying bc_params
		eenergy_bcparams = bc_params

		!Set dirichlet conditions from boundary
		!Flux/zgrd conditions dont require boundary terms
		if(g_lrank .lt. 0) then
	
			i = 1
			do k=1,g_lz
				do j=1,g_ly
					eenergy_bcvals%left(j,k)=bc_params(LEFT)/eenergyscale
				enddo
			enddo

		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			do k=1,g_lz
				do j=1,g_ly
					eenergy_bcvals%right(j,k)=bc_params(RIGHT)/eenergyscale
				enddo
			enddo
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			do k=1,g_lz
				do i=1,g_lx
					eenergy_bcvals%bottom(i,k)=bc_params(BOTTOM)/eenergyscale
				enddo
			enddo
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			do k=1,g_lz
				do i=1,g_lx
					eenergy_bcvals%top(i,k)=bc_params(TOP)/eenergyscale
				enddo
			enddo
		endif

		if(g_krank .lt. 0) then

			k = 1
			do j=1,g_ly
				do i=1,g_lx
					eenergy_bcvals%back(i,j)=bc_params(BACK)/eenergyscale
				enddo
			enddo
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			do j=1,g_ly
				do i=1,g_lx
					eenergy_bcvals%front(i,j)=bc_params(FRONT)/eenergyscale
				enddo
			enddo
		endif


      end subroutine eenergy_set_bcs
!===========================================================================================
      subroutine eenergy_update_bcs(elecenergy,elecden)

		real*8, intent(in) :: elecden (-g_nglayers+1:g_lx+g_nglayers,&
					       -g_nglayers+1:g_ly+g_nglayers,&
					       -g_nglayers+1:g_lz+g_nglayers)
		
	        !elecenergy is scaled electron energy
		real*8, intent(in) :: elecenergy(-g_nglayers+1:g_lx+g_nglayers,&
						 -g_nglayers+1:g_ly+g_nglayers,&
						 -g_nglayers+1:g_lz+g_nglayers)


		real*8 :: cbar,const_1,etemp,flux
		real*8 :: three_by_two,four_by_three
		integer :: i,j,k

		three_by_two = 1.5d0
		four_by_three = 4.d0/3.d0
		const_1 = 8.d0*k_BMAN/PI/M_ELEC

		if(g_lrank .lt. 0) then
	
			i = 1
			if(eenergy_bc_codes(LEFT) .eq. 'DIRC') then

				do k=1,g_lz
					do j=1,g_ly
						eenergy_bcvals%left(j,k)=three_by_two*k_BMAN*&
								elecden(i,j,k)*eenergy_bcparams(LEFT)/eenergyscale
					enddo
				enddo

			else if(eenergy_bc_codes(LEFT) .eq. 'FLUX') then
				
				do k=1,g_lz
					do j=1,g_ly
						
						etemp = elecenergy(i,j,k)*eenergyscale/(three_by_two*elecden(i,j,k)*k_BMAN)	
						cbar = sqrt(const_1*etemp)	
						flux = four_by_three*elecenergy(i,j,k)*cbar

						eenergy_bcvals%left(j,k) = flux*vel_scaling
					enddo
				enddo
			endif

		endif
		
		if(g_rrank .lt. 0) then

			i = g_lx
			if(eenergy_bc_codes(RIGHT) .eq. 'DIRC') then
	
				do k=1,g_lz
					do j=1,g_ly
						eenergy_bcvals%right(j,k)=three_by_two*k_BMAN*&
								elecden(i,j,k)*eenergy_bcparams(RIGHT)/eenergyscale
					enddo
				enddo

			else if(eenergy_bc_codes(RIGHT) .eq. 'FLUX') then
				
				do k=1,g_lz
					do j=1,g_ly
						
						etemp = elecenergy(i,j,k)*eenergyscale/(three_by_two*elecden(i,j,k)*k_BMAN)	
						cbar = sqrt(const_1*etemp)	
						flux = four_by_three*elecenergy(i,j,k)*cbar
						
						eenergy_bcvals%right(j,k) = flux*vel_scaling
					enddo
				enddo
			endif
		endif
			
		if(g_brank .lt. 0) then

			j = 1
			if(eenergy_bc_codes(BOTTOM) .eq. 'DIRC') then

				do k=1,g_lz
					do i=1,g_lx
						eenergy_bcvals%bottom(i,k)=three_by_two*k_BMAN*&
								elecden(i,j,k)*eenergy_bcparams(BOTTOM)/eenergyscale
					enddo
				enddo

			else if(eenergy_bc_codes(BOTTOM) .eq. 'FLUX') then
				
				do k=1,g_lz
					do i=1,g_lx
						
						etemp = elecenergy(i,j,k)*eenergyscale/(three_by_two*elecden(i,j,k)*k_BMAN)	
						cbar = sqrt(const_1*etemp)	
						flux = four_by_three*elecenergy(i,j,k)*cbar
						
						eenergy_bcvals%bottom(i,k) = flux*vel_scaling
					enddo
				enddo
			endif
		endif
		
		if(g_trank .lt. 0) then

			j = g_ly
			if(eenergy_bc_codes(TOP) .eq. 'DIRC') then

				do k=1,g_lz
					do i=1,g_lx
						eenergy_bcvals%top(i,k)=three_by_two*k_BMAN*&
								elecden(i,j,k)*eenergy_bcparams(TOP)/eenergyscale
					enddo
				enddo
			
			else if(eenergy_bc_codes(TOP) .eq. 'FLUX') then
				
				do k=1,g_lz
					do i=1,g_lx
						
						etemp = elecenergy(i,j,k)*eenergyscale/(three_by_two*elecden(i,j,k)*k_BMAN)	
						cbar = sqrt(const_1*etemp)	
						flux = four_by_three*elecenergy(i,j,k)*cbar
						
						eenergy_bcvals%top(i,k) = flux*vel_scaling
					enddo
				enddo
			endif
		endif

		if(g_krank .lt. 0) then

			k = 1
			if(eenergy_bc_codes(BACK) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						eenergy_bcvals%back(i,j)=three_by_two*k_BMAN*&
								elecden(i,j,k)*eenergy_bcparams(BACK)/eenergyscale
					enddo
				enddo
			
			else if(eenergy_bc_codes(BACK) .eq. 'FLUX') then
				
				do j=1,g_ly
					do i=1,g_lx
						
						etemp = elecenergy(i,j,k)*eenergyscale/(three_by_two*elecden(i,j,k)*k_BMAN)	
						cbar = sqrt(const_1*etemp)	
						flux = four_by_three*elecenergy(i,j,k)*cbar
						
						eenergy_bcvals%back(i,j) = flux*vel_scaling
					enddo
				enddo
			endif
		endif
		
		if(g_frank .lt. 0) then

			k = g_lz
			if(eenergy_bc_codes(FRONT) .eq. 'DIRC') then

				do j=1,g_ly
					do i=1,g_lx
						eenergy_bcvals%front(i,j)=three_by_two*k_BMAN*&
								elecden(i,j,k)*eenergy_bcparams(FRONT)/eenergyscale
					enddo
				enddo
			
			else if(eenergy_bc_codes(FRONT) .eq. 'FLUX') then
				
				do j=1,g_ly
					do i=1,g_lx
						
						etemp = elecenergy(i,j,k)*eenergyscale/(three_by_two*elecden(i,j,k)*k_BMAN)	
						cbar = sqrt(const_1*etemp)	
						flux = four_by_three*elecenergy(i,j,k)*cbar
						
						eenergy_bcvals%front(i,j) = flux*vel_scaling
					enddo
				enddo
			endif
		endif


      end subroutine eenergy_update_bcs
!===========================================================================================
      subroutine jouleheatingterm(elecfield,electemp,numden,Tg,Pg,jheat)
      
		real*8, intent(in) :: elecfield(-g_nglayers+1:g_lx+g_nglayers,&
						-g_nglayers+1:g_ly+g_nglayers,&
						-g_nglayers+1:g_lz+g_nglayers,NDIM)

		real*8, intent(in) :: numden(-g_nglayers+1:g_lx+g_nglayers,&
					     -g_nglayers+1:g_ly+g_nglayers,&
					     -g_nglayers+1:g_lz+g_nglayers,nspecies)
		
		real*8, intent(in) ::  electemp(-g_nglayers+1:g_lx+g_nglayers,&
						-g_nglayers+1:g_ly+g_nglayers,&
				  		-g_nglayers+1:g_lz+g_nglayers)
		real*8, intent(in) ::  Tg,Pg
		real*8, intent(out) :: jheat(g_lx,g_ly,g_lz)

		integer :: i,j,k,nsp

		real*8 :: elecden(-g_nglayers+1:g_lx+g_nglayers,&
				  -g_nglayers+1:g_ly+g_nglayers,&
				  -g_nglayers+1:g_lz+g_nglayers)

		real*8 :: jheatx,jheaty,jheatz
		real*8 :: mu_e,D_e,efield_mag

		real*8 :: etemp,three_by_two,specarray(nspecies)

		three_by_two = 1.5d0
		elecden = numden(:,:,:,especnum)

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
					
					etemp = electemp(i,j,k)
					
					mu_e      =   getspecmobility(especnum,specarray,efield_mag,&
								etemp,Tg,Pg)
					
					D_e       =    getspecdcoeff(especnum,specarray,efield_mag,&
									etemp,Tg,Pg)

					call findjheating_eachdir(i,j,k,g_lx,elecden,elecfield,mu_e,D_e,XDIR,jheatx)
					call findjheating_eachdir(i,j,k,g_ly,elecden,elecfield,mu_e,D_e,YDIR,jheaty)
					call findjheating_eachdir(i,j,k,g_lz,elecden,elecfield,mu_e,D_e,ZDIR,jheatz)

					jheat(i,j,k) = jheatx+jheaty+jheatz
					
				enddo
			enddo
		enddo
		
      end subroutine jouleheatingterm 
!===========================================================================================
      subroutine findjheating_eachdir(ii,jj,kk,np,elecden,elecfield,mu_e,D_e,dir,jheat)

      	        integer, intent(in) :: ii,jj,kk
		integer, intent(in) :: np
		integer, intent(in) :: dir
		
		real*8, intent(in) :: elecden(-g_nglayers+1:g_lx+g_nglayers,&
					      -g_nglayers+1:g_ly+g_nglayers,&
					      -g_nglayers+1:g_lz+g_nglayers)

		real*8, intent(in) :: elecfield(-g_nglayers+1:g_lx+g_nglayers,&
						-g_nglayers+1:g_ly+g_nglayers,&
						-g_nglayers+1:g_lz+g_nglayers,NDIM)

		real*8, intent(in) :: mu_e,D_e

		real*8, intent(out) :: jheat

		real*8,allocatable :: ne(:)
		real*8,allocatable :: efield(:)

		integer :: i
		real*8 :: grad_ne
		real*8 :: Gamma_left,Gamma_right,efield_mag
		real*8 :: elec_field_left,elec_field_right,e_den
		real*8 :: dx
		
		Gamma_left       = ZERO
		Gamma_right      = ZERO
		elec_field_left  = ZERO
		elec_field_right = ZERO

		if(dir .eq. XDIR) then

			allocate(ne(g_lx+2*g_nglayers))
			allocate(efield(g_lx+2*g_nglayers))

			ne = elecden(:,jj,kk)
			efield = elecfield(:,jj,kk,XDIR)

			i = ii+1
			dx = g_dx*lscale

		else if(dir .eq. YDIR) then

			allocate(ne(g_ly+2*g_nglayers))
			allocate(efield(g_ly+2*g_nglayers))

			ne = elecden(ii,:,kk)
			efield = elecfield(ii,:,kk,YDIR)

			i = jj+1
			dx = g_dy*lscale

		else
			allocate(ne(g_lz+2*g_nglayers))
			allocate(efield(g_lz+2*g_nglayers))

			ne = elecden(ii,jj,:)
			efield = elecfield(ii,jj,:,ZDIR)

			i = kk+1
			dx = g_dz*lscale
		endif

		
		!Left side
		grad_ne =  (ne(i)-ne(i-1))/dx
		elec_field_left  = (efield(i)+efield(i-1))*0.5
	
		if(mu_e*elec_field_left .gt. 0) then
			e_den = ne(i-1)
		else
			e_den = ne(i)
		endif
		Gamma_left = (mu_e*e_den*elec_field_left - D_e*grad_ne)
		
	
		!Right side		
		grad_ne =  (ne(i+1)-ne(i))/dx
		elec_field_right  = (efield(i)+efield(i+1))*0.5
			
		if(mu_e*elec_field_right .gt. 0) then
			e_den = ne(i)
		else
			e_den = ne(i+1)
		endif
		Gamma_right = (mu_e*e_den*elec_field_right - D_e*grad_ne)

		
	
		jheat = -ECHARGE*0.5*(elec_field_left+elec_field_right)*&
					0.5*(Gamma_left+Gamma_right)

		if(jheat .lt. ZERO) then
			jheat=ZERO
		endif


      end subroutine findjheating_eachdir
!===========================================================================================
      subroutine elastic_colterm(elecfield,electemp,numden,Tg,Pg,elastic_col)

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
		real*8, intent(out) :: elastic_col(g_lx,g_ly,g_lz)

		integer :: i,j,k,nsp

		real*8 :: elecden(-g_nglayers+1:g_lx+g_nglayers,&
				  -g_nglayers+1:g_ly+g_nglayers,&
				  -g_nglayers+1:g_lz+g_nglayers)

		real*8 :: etemp,efield_mag
	
		real*8 :: nu
		real*8 :: m_e
		real*8 :: m_bg
		real*8 :: specarray(nspecies)

		m_e  = molmass(especnum)
		m_bg = molmass(bgspecnum)

		elecden = numden(:,:,:,especnum)
	
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

				etemp = electemp(i,j,k)

				nu = getelectroncollisionfrequency(specarray,efield_mag,etemp,Tg,Pg)

				elastic_col(i,j,k)= 3.d0*k_BMAN*elecden(i,j,k)*(etemp-Tg)*(m_e/m_bg)*nu
			enddo
		    enddo
		enddo
	
end subroutine elastic_colterm
!=============================================================================
subroutine inelastic_colterm(elecfield,electemp,numden,Tg,Pg,inelastic_col)
	
		real*8, intent(in) :: elecfield(-g_nglayers+1:g_lx+g_nglayers,&
						-g_nglayers+1:g_ly+g_nglayers,&
						-g_nglayers+1:g_lz+g_nglayers,NDIM)

		real*8, intent(in) :: numden(-g_nglayers+1:g_lx+g_nglayers,&
					     -g_nglayers+1:g_ly+g_nglayers,&
					     -g_nglayers+1:g_lz+g_nglayers,nspecies)
		
		real*8, intent(in) :: electemp(-g_nglayers+1:g_lx+g_nglayers,&
					       -g_nglayers+1:g_ly+g_nglayers,&
					       -g_nglayers+1:g_lz+g_nglayers)
		real*8, intent(in) :: Tg,Pg
		real*8, intent(out) :: inelastic_col(g_lx,g_ly,g_lz)

		integer :: i,j,k,nsp

		real*8 :: elecden(-g_nglayers+1:g_lx+g_nglayers,&
				  -g_nglayers+1:g_ly+g_nglayers,&
				  -g_nglayers+1:g_lz+g_nglayers)

		real*8 :: etemp,three_by_two,efield_mag
		real*8 :: specarray(nspecies)
		real*8 :: inelterm

		three_by_two = 1.5d0
		elecden = numden(:,:,:,especnum)

		do k=1,g_lz
			do j=1,g_ly
				do i=1,g_lx
	
					do nsp=1,nspecies
						specarray(nsp)=numden(i,j,k,nsp)
					enddo

					etemp = electemp(i,j,k)
				
					efield_mag = elecfield(i,j,k,XDIR)**2
					efield_mag = efield_mag + elecfield(i,j,k,YDIR)**2
					efield_mag = efield_mag + elecfield(i,j,k,ZDIR)**2

					efield_mag = sqrt(efield_mag)
	  
	  				call getelectroninelasticterm(etemp,Tg,&
			  			specarray,inelterm,efield_mag)

					inelastic_col(i,j,k) = inelterm*ECHARGE !convert eV to joules

	 		 	enddo
			enddo
		 enddo	

end subroutine inelastic_colterm
!=============================================================================
      subroutine eenergy_update_transport(dt,elecpot,elecfield,electemp,numden,Tg,Pg)
	
		real*8, intent(in) :: dt

		real*8, intent(in) :: elecfield(-g_nglayers+1:g_lx+g_nglayers,&
						-g_nglayers+1:g_ly+g_nglayers,&
						-g_nglayers+1:g_lz+g_nglayers,NDIM)
		
		real*8, intent(in) :: elecpot(-g_nglayers+1:g_lx+g_nglayers,&
					      -g_nglayers+1:g_ly+g_nglayers,&
					      -g_nglayers+1:g_lz+g_nglayers)

		real*8, intent(in) :: numden(-g_nglayers+1:g_lx+g_nglayers,&
					     -g_nglayers+1:g_ly+g_nglayers,&
					     -g_nglayers+1:g_lz+g_nglayers,nspecies)
		
		!this is scaled (non-dimensional)
		real*8, intent(in) :: electemp(-g_nglayers+1:g_lx+g_nglayers,&
				               -g_nglayers+1:g_ly+g_nglayers,&
					       -g_nglayers+1:g_lz+g_nglayers)
		
		real*8, intent(in) :: Tg,Pg

		real*8 :: elecden(-g_nglayers+1:g_lx+g_nglayers,&
				  -g_nglayers+1:g_ly+g_nglayers,&
				  -g_nglayers+1:g_lz+g_nglayers)

		real*8 :: efield_mag
		real*8 :: mu_e

	        integer :: i,j,k,nsp

		real*8 :: specarray(nspecies)
		real*8 :: specprod
		real*8 :: five_by_three
		real*8 :: three_by_two
		real*8 :: etemp

		real*8 :: jheatingterm(g_lx,g_ly,g_lz)
		real*8 :: elasticterm(g_lx,g_ly,g_lz)
		real*8 :: inelasticterm(g_lx,g_ly,g_lz)

		three_by_two  = 1.5d0
		five_by_three = 5.d0/3.d0

		elecden = numden(:,:,:,especnum)

		vel    = ZERO
		dcoeff = ZERO
		reac   = ZERO
		source = ZERO

      		call jouleheatingterm(elecfield,electemp,numden,Tg,Pg,jheatingterm)
		call inelastic_colterm(elecfield,electemp,numden,Tg,Pg,inelasticterm)
      		call elastic_colterm(elecfield,electemp,numden,Tg,Pg,elasticterm)

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

					etemp = electemp(i,j,k)

					mu_e      =   getspecmobility(especnum,specarray,efield_mag,&
								etemp,Tg,Pg)

					
					vel(i+1,j+1,k+1,XDIR) = five_by_three * mu_e*elecfield(i,j,k,XDIR)
					vel(i+1,j+1,k+1,YDIR) = five_by_three * mu_e*elecfield(i,j,k,YDIR)
					vel(i+1,j+1,k+1,ZDIR) = five_by_three * mu_e*elecfield(i,j,k,ZDIR)
		
					dcoeff(i+1,j+1,k+1)   = five_by_three * getspecdcoeff(especnum,&
									specarray,efield_mag,&
									etemp,Tg,Pg)
					

				enddo
			enddo
		enddo

		source(2:g_lx+1,2:g_ly+1,2:g_lz+1) = jheatingterm - elasticterm + inelasticterm

		dcoeff   = dcoeff * dcoeff_scaling
		vel      = vel    * vel_scaling
		reac     = reac   * reac_scaling
		source   = source * source_scaling
		
		if(timederivflag .eqv. .true.) source = source + eenergysoln/dt
		

      end subroutine eenergy_update_transport
!=========================================================================================
      subroutine eenergy_solve(time,dt,elecpot,elecfield,electemp,numden,Tg,Pg,printflag)

	real*8, intent(in) :: time,dt,Tg,Pg
	logical, intent(in) :: printflag

	!passing electron temperature to be consistent with 
	!other equations. Also advantageous because I dont have
	!to calculate electron temperature everytime with ee = 3/2 ne kB T
	real*8, intent(in) :: electemp(-g_nglayers+1:g_lx+g_nglayers,&
				       -g_nglayers+1:g_ly+g_nglayers,&
				       -g_nglayers+1:g_lz+g_nglayers)

	real*8, intent(in) :: elecfield(-g_nglayers+1:g_lx+g_nglayers,&
					-g_nglayers+1:g_ly+g_nglayers,&
					-g_nglayers+1:g_lz+g_nglayers,NDIM)
	
	real*8, intent(in) :: elecpot(-g_nglayers+1:g_lx+g_nglayers,&
				      -g_nglayers+1:g_ly+g_nglayers,&
				      -g_nglayers+1:g_lz+g_nglayers)

	real*8, intent(in) :: numden(-g_nglayers+1:g_lx+g_nglayers,&
				     -g_nglayers+1:g_ly+g_nglayers,&
				     -g_nglayers+1:g_lz+g_nglayers,nspecies)

	integer :: i
	integer :: nvcycles
	real*8 :: resnorm
	real*8,allocatable:: sterm(:,:,:)
	real*8 :: solve_time

	real*8 :: err_tol,dt_s
	
	err_tol=ERRTOLR
	nvcycles=20

	dt_s = dt/tscale
	
	allocate(sterm(g_lx,g_ly,g_lz))
	      
	if(g_myproc .eq. g_rootproc) solve_time=-MPI_Wtime()

	call eenergy_update_transport(dt_s,elecpot,elecfield,electemp,numden,Tg,Pg)
	call eenergy_update_bcs(eenergysoln,numden(:,:,:,especnum))
	
	do i=1,nvcycles

		sterm = 0.d0
		
		call perform_vcycle(eenergysoln,timederivflag,dt_s,g_lx,g_ly,g_lz,vel,dcoeff,reac,source,&
				eenergy_bc_codes,eenergy_bcvals,sterm,resnorm)
	
		if((printflag .eqv. .true.) .and. (g_myproc .eq. g_rootproc)) print *,"it:",i,resnorm

		if(resnorm .le. err_tol) exit
	
	enddo
	
	call exchangehalodata(eenergysoln,g_lx,g_ly,g_lz)
	call update_boundary_ghostvalues(eenergysoln,eenergy_bc_codes,&
				eenergy_bcvals,g_lx,g_ly,g_lz)

	if(g_myproc .eq. g_rootproc) solve_time=solve_time+MPI_Wtime()
	if((printflag .eqv. .true.) .and. (g_myproc .eq. g_rootproc)) print *,"solve time:",solve_time

      end subroutine eenergy_solve
!=================================================================
      subroutine update_electrontemp(electemp,elecden)

	    real*8, intent(out) :: electemp(-g_nglayers+1:g_lx+g_nglayers,&
			    		    -g_nglayers+1:g_ly+g_nglayers,&
					    -g_nglayers+1:g_lz+g_nglayers)
	    
	    real*8, intent(in) :: elecden(-g_nglayers+1:g_lx+g_nglayers,&
			    		  -g_nglayers+1:g_ly+g_nglayers,&
					  -g_nglayers+1:g_lz+g_nglayers)

	    real*8 :: three_by_two
	    integer :: i,j,k

	    three_by_two=1.5


	    do k=1,g_lz
	    	do j=1,g_ly
		    do i=1,g_lx

		    	    electemp(i,j,k)=eenergysoln(i+1,j+1,k+1)*eenergyscale&
					    /(three_by_two*elecden(i,j,k)*k_BMAN)

			    if(electemp(i,j,k) .lt. etempmin) then
				   
				    electemp(i,j,k)    = etempmin
				    eenergysoln(i,j,k) = three_by_two*elecden(i,j,k)*&
						    	k_BMAN*etempmin/eenergyscale

			    endif

		    enddo
		enddo
	   enddo

      end subroutine update_electrontemp
!=================================================================
end module eenergyeq
