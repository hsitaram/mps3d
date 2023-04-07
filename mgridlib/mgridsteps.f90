module mgridsteps

      use globalvars
      use discretepde
      use linearsolvers
      use solvergmres_module
      implicit none
      contains

!========================================================================================================================
subroutine restriction(phi_fine,lx,ly,lz,phi_coarse)

	integer, intent(in) :: lx,ly,lz
	real*8, intent(in)  :: phi_fine(lx,ly,lz)
      	real*8,  intent(out) :: phi_coarse(lx/2,ly/2,lz/2)

	integer :: lcx,lcy,lcz
	integer :: i,j,k
	real*8 :: one_eighth

	lcx=lx/2
	lcy=ly/2
	lcz=lz/2

	one_eighth=1.d0/8.d0
	phi_coarse=0.d0

		
	do k=1,lcz

		do j=1,lcy
	
			do i=1,lcx

				phi_coarse(i,j,k) = phi_coarse(i,j,k) + phi_fine(2*i-1,2*j-1,2*k-1) + &
									      phi_fine(2*i-0,2*j-1,2*k-1)
								      
				phi_coarse(i,j,k) = phi_coarse(i,j,k) + phi_fine(2*i-1,2*j-0,2*k-1) + &
									      phi_fine(2*i-0,2*j-0,2*k-1)

				phi_coarse(i,j,k) = phi_coarse(i,j,k) + phi_fine(2*i-1,2*j-1,2*k-0) + &
									      phi_fine(2*i-0,2*j-1,2*k-0)

				phi_coarse(i,j,k) = phi_coarse(i,j,k) + phi_fine(2*i-1,2*j-0,2*k-0) + &
									      phi_fine(2*i-0,2*j-0,2*k-0)

				phi_coarse(i,j,k) = one_eighth*phi_coarse(i,j,k)

			enddo
		enddo
	enddo

end subroutine restriction
!========================================================================================================================
subroutine restrict_boundary(bcvals_fine,lx,ly,lz,bcvals_coarse)

	type(boundarycondition), intent(inout) :: bcvals_fine
	integer, intent(in) :: lx,ly,lz
	type(boundarycondition), intent(inout) :: bcvals_coarse
	
	real*8  :: one_fourth
	integer :: lcx,lcy,lcz
	integer :: i,j,k

	lcx=lx/2
	lcy=ly/2
	lcz=lz/2

	one_fourth = 0.25

	!LEFT and RIGHT
	do k=1,lcz
       	
		do j=1,lcy
			
			bcvals_coarse%left(j,k) = one_fourth*&
						 (bcvals_fine%left(2*j-1,2*k-1) + &
						  bcvals_fine%left(2*j-0,2*k-1) + &
						  bcvals_fine%left(2*j-1,2*k-0) + &
						  bcvals_fine%left(2*j-0,2*k-0))
			
			bcvals_coarse%right(j,k) = one_fourth*&
						 (bcvals_fine%right(2*j-1,2*k-1) + &
						  bcvals_fine%right(2*j-0,2*k-1) + &
						  bcvals_fine%right(2*j-1,2*k-0) + &
						  bcvals_fine%right(2*j-0,2*k-0))
		enddo
	enddo
	
	!BOTTOM and TOP
	do k=1,lcz
       	
		do i=1,lcx
			
			bcvals_coarse%bottom(i,k) = one_fourth*&
						 (bcvals_fine%bottom(2*i-1,2*k-1) + &
						  bcvals_fine%bottom(2*i-0,2*k-1) + &
						  bcvals_fine%bottom(2*i-1,2*k-0) + &
						  bcvals_fine%bottom(2*i-0,2*k-0))
			
			bcvals_coarse%top(i,k)  = one_fourth*&
						 (bcvals_fine%top(2*i-1,2*k-1) + &
						  bcvals_fine%top(2*i-0,2*k-1) + &
						  bcvals_fine%top(2*i-1,2*k-0) + &
						  bcvals_fine%top(2*i-0,2*k-0))
		enddo
	enddo
	
	!BACK and FRONT
	do j=1,lcy
       	
		do i=1,lcx
			
			bcvals_coarse%back(i,j) = one_fourth*&
						 (bcvals_fine%back(2*i-1,2*j-1) + &
						  bcvals_fine%back(2*i-0,2*j-1) + &
						  bcvals_fine%back(2*i-1,2*j-0) + &
						  bcvals_fine%back(2*i-0,2*j-0))
			
			bcvals_coarse%front(i,j) = one_fourth*&
						 (bcvals_fine%front(2*i-1,2*j-1) + &
						  bcvals_fine%front(2*i-0,2*j-1) + &
						  bcvals_fine%front(2*i-1,2*j-0) + &
						  bcvals_fine%front(2*i-0,2*j-0))
		enddo
	enddo

end subroutine restrict_boundary
!========================================================================================================================
subroutine prolongation0(phi_coarse,lx,ly,lz,phi_fine)

	integer, intent(in) :: lx,ly,lz
	real*8, intent(in)  :: phi_coarse(lx,ly,lz)

	real*8, intent(out)::  phi_fine(2*lx,2*ly,2*lz)

	integer :: lfx,lfy,lfz
	integer :: i,j,k
	integer :: ic,jc,kc
	
	lfx = 2*lx
	lfy = 2*ly
	lfz = 2*lz
		
	!print *,"g_myproc:",g_myproc,lfx,lfy,lfz

	do k=1,lfz

		if(mod(k,2) .eq. 0) kc=k/2
		if(mod(k,2) .eq. 1) kc=(k+1)/2	

		do j=1,lfy
				
			if(mod(j,2) .eq. 0) jc=j/2
			if(mod(j,2) .eq. 1) jc=(j+1)/2	

			do i=1,lfx
			
				if(mod(i,2) .eq. 0) ic=i/2
				if(mod(i,2) .eq. 1) ic=(i+1)/2	
				
				phi_fine(i,j,k) = phi_coarse(ic,jc,kc)
			
			enddo
		enddo
	enddo

end subroutine prolongation0
!=============================================================================================================================
subroutine prolongation1(phi_coarse,lx,ly,lz,phi_fine)

	integer, intent(in) :: lx,ly,lz
	real*8, intent(in)  :: phi_coarse(-g_nglayers+1:lx+g_nglayers, &
					  -g_nglayers+1:ly+g_nglayers, &
					  -g_nglayers+1:lz+g_nglayers)

	real*8, intent(out)::  phi_fine(-g_nglayers+1:2*lx+g_nglayers, &
					-g_nglayers+1:2*ly+g_nglayers, &
					-g_nglayers+1:2*lz+g_nglayers)

	integer :: lfx,lfy,lfz
	integer :: i,j,k          !fine index
	integer :: ic,jc,kc       !coarse index
	
	real*8 :: dphi_x,dphi_y,dphi_z

	lfx = 2*lx
	lfy = 2*ly
	lfz = 2*lz
		
	!print *,"g_myproc:",g_myproc,lfx,lfy,lfz

	do kc=1,lz
		do jc=1,ly
			do ic=1,lx
				
				dphi_x = 0.5*( phi_coarse(ic+1,jc,kc) - phi_coarse(ic-1,jc,kc) )
				dphi_y = 0.5*( phi_coarse(ic,jc+1,kc) - phi_coarse(ic,jc-1,kc) )
				dphi_z = 0.5*( phi_coarse(ic,jc,kc+1) - phi_coarse(ic,jc,kc-1) )

				phi_fine(2*ic-1, 2*jc-1, 2*kc-1) = phi_coarse(ic,jc,kc) + 0.25*(-dphi_x -dphi_y -dphi_z)
				phi_fine(2*ic-0, 2*jc-1, 2*kc-1) = phi_coarse(ic,jc,kc) + 0.25*(+dphi_x -dphi_y -dphi_z)
				phi_fine(2*ic-1, 2*jc-0, 2*kc-1) = phi_coarse(ic,jc,kc) + 0.25*(-dphi_x +dphi_y -dphi_z)
				phi_fine(2*ic-0, 2*jc-0, 2*kc-1) = phi_coarse(ic,jc,kc) + 0.25*(+dphi_x +dphi_y -dphi_z)
				
				phi_fine(2*ic-1, 2*jc-1, 2*kc-0) = phi_coarse(ic,jc,kc) + 0.25*(-dphi_x -dphi_y +dphi_z)
				phi_fine(2*ic-0, 2*jc-1, 2*kc-0) = phi_coarse(ic,jc,kc) + 0.25*(+dphi_x -dphi_y +dphi_z)
				phi_fine(2*ic-1, 2*jc-0, 2*kc-0) = phi_coarse(ic,jc,kc) + 0.25*(-dphi_x +dphi_y +dphi_z)
				phi_fine(2*ic-0, 2*jc-0, 2*kc-0) = phi_coarse(ic,jc,kc) + 0.25*(+dphi_x +dphi_y +dphi_z)
				

			enddo
		enddo
	enddo

end subroutine prolongation1
!=============================================================================================================================
subroutine bottomsolve(solnvar,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,source,bc_codes,serial_bcvals,&
				applybcflag,mgridsterm,resnorm)

	!need to obtain serial_bcvals
	!not necessary for now, as bcflags are set to false.
	real*8, intent(in)  :: dt
	integer, intent(in) :: lx,ly,lz
	real*8, intent(inout) :: solnvar(-g_nglayers + 1:lx + g_nglayers, &
  				         -g_nglayers + 1:ly + g_nglayers, &
				         -g_nglayers + 1:lz + g_nglayers)
        
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

	real*8, intent(in)    :: mgridsterm(lx,ly,lz)
	
	character(LEN=4),intent(in)         :: bc_codes(NFACES)
	type(boundarycondition), intent(inout) :: serial_bcvals

	logical, intent(in) :: timederivflag
	logical, intent(in) :: applybcflag

	real*8, intent(out) :: resnorm

	integer :: nx,ny,nz
	integer :: ierr,onedindex
	integer :: offx,offy,offz
	integer :: lrank,rrank,brank,trank,krank,frank

	integer :: i,j,k
	integer :: itmax
	real*8  :: serial_llenx,serial_lleny,serial_llenz

	real*8  :: serialsoln(-g_nglayers + 1:lx*g_px + g_nglayers, &
			      -g_nglayers + 1:ly*g_py + g_nglayers, &
			      -g_nglayers + 1:lz*g_pz + g_nglayers)

	real*8  :: serialsoln0(-g_nglayers + 1:lx*g_px + g_nglayers, &
			      -g_nglayers + 1:ly*g_py + g_nglayers, &
			      -g_nglayers + 1:lz*g_pz + g_nglayers)

	real*8  :: serial_mgridsterm(lx*g_px,ly*g_py,lz*g_pz)
	
	real*8  :: serial_vel(-g_nglayers+1:lx*g_px+g_nglayers,&
			      -g_nglayers+1:ly*g_py+g_nglayers,&
			      -g_nglayers+1:lz*g_pz+g_nglayers,NDIM)

	real*8  :: serial_dcoeff(-g_nglayers+1:lx*g_px+g_nglayers,&
			         -g_nglayers+1:ly*g_py+g_nglayers,&
				 -g_nglayers+1:lz*g_pz+g_nglayers)

	real*8  :: serial_reac(-g_nglayers+1:lx*g_px+g_nglayers,&
			       -g_nglayers+1:ly*g_py+g_nglayers,&
			       -g_nglayers+1:lz*g_pz+g_nglayers)

	real*8  :: serial_source(-g_nglayers+1:lx*g_px+g_nglayers,&
				 -g_nglayers+1:ly*g_py+g_nglayers,&
				 -g_nglayers+1:lz*g_pz+g_nglayers)

	real*8 :: res(lx*g_px,ly*g_py,lz*g_pz)
	real*8 :: b(lx*g_px,ly*g_py,lz*g_pz)
	real*8 :: err_tol
	logical :: err_achieved

	integer :: pindx,pindy,pindz,proc
	integer :: maxkspdim,nrestarts
	logical :: gmresprintflag
	real*8  :: initial_res_gmres

	nx = lx*g_px
	ny = ly*g_py
	nz = lz*g_pz
	
	itmax=20*nx*ny*nz
	!itmax=835
	err_tol=ERRTOLR*0.01

	serialsoln    = ZERO
	res           = ZERO
	serial_vel    = ZERO
	serial_dcoeff = ZERO
	serial_reac   = ZERO
	serial_source = ZERO

	call gatherterms(solnvar(1:lx,1:ly,1:lz), lx,ly,lz,serialsoln(1:lx*g_px,1:ly*g_py,1:lz*g_pz))

	call gatherterms(vel(1:lx,1:ly,1:lz,XDIR),lx,ly,lz,serial_vel(1:lx*g_px,1:ly*g_py,1:lz*g_pz,XDIR))
	call gatherterms(vel(1:lx,1:ly,1:lz,YDIR),lx,ly,lz,serial_vel(1:lx*g_px,1:ly*g_py,1:lz*g_pz,YDIR))
	call gatherterms(vel(1:lx,1:ly,1:lz,ZDIR),lx,ly,lz,serial_vel(1:lx*g_px,1:ly*g_py,1:lz*g_pz,ZDIR))
	
	call gatherterms(dcoeff(1:lx,1:ly,1:lz),lx,ly,lz,serial_dcoeff(1:lx*g_px,1:ly*g_py,1:lz*g_pz))

	call gatherterms(reac  (1:lx,1:ly,1:lz),lx,ly,lz,serial_reac  (1:lx*g_px,1:ly*g_py,1:lz*g_pz))
	
	call gatherterms(source(1:lx,1:ly,1:lz),lx,ly,lz,serial_source(1:lx*g_px,1:ly*g_py,1:lz*g_pz))
	
	call gatherterms(mgridsterm,lx,ly,lz,serial_mgridsterm)

	lrank=MPI_PROC_NULL
	rrank=MPI_PROC_NULL
	brank=MPI_PROC_NULL
	trank=MPI_PROC_NULL
	krank=MPI_PROC_NULL
	frank=MPI_PROC_NULL
	serial_llenx = g_llenx*g_px
	serial_lleny = g_lleny*g_py
	serial_llenz = g_llenz*g_pz

	if(g_myproc .eq. g_rootproc) then

		err_achieved=.false.
		serialsoln0 = serialsoln
		maxkspdim = nx*ny*nz/2
		if(maxkspdim .lt. 2) maxkspdim=6
		if(maxkspdim .gt. 10) maxkspdim=10
		nrestarts = 10
		gmresprintflag = .false.
		
		call find_bvec(serialsoln,dt,nx,ny,nz,serial_mgridsterm,serial_vel,serial_dcoeff,serial_source,&
					bc_codes,serial_bcvals,applybcflag,b,&
					lrank,rrank,brank,trank,krank,frank,serial_llenx,serial_lleny,serial_llenz)

  		call performgmres(b,serialsoln0,serialsoln,timederivflag,dt,nx,ny,nz,serial_vel,serial_dcoeff,&
				  serial_reac,bc_codes,serial_bcvals,&
				  lrank,rrank,brank,trank,krank,frank,&
				  serial_llenx,serial_lleny,serial_llenz,&
				  maxkspdim,nrestarts,find_AX,&
				  gseidel_precond,err_tol,err_achieved,gmresprintflag,initial_res_gmres,resnorm)
			
		
		if(err_achieved .eqv. .false.) then
			print *,"hit max iterations",itmax,initial_res_gmres,resnorm
			call MPI_Abort(g_comm,666,ierr)
		endif

	endif
	call MPI_Barrier(g_comm,ierr)

	call MPI_Bcast(serialsoln,(nx+2*g_nglayers)*(ny+2*g_nglayers)*(nz+2*g_nglayers),&
			MPI_REAL8,g_rootproc,g_comm,ierr)


	offx = (g_pindx-1)*lx
	offy = (g_pindy-1)*ly
	offz = (g_pindz-1)*lz

	do k=1,lz
		do j=1,ly
			do i=1,lx
				solnvar(i,j,k)=serialsoln(i+offx,j+offy,k+offz)
			enddo
		enddo
	enddo


end subroutine bottomsolve
!=============================================================================================================================
recursive subroutine  perform_vcycle(solnvar,timederivflag,dt,lx,ly,lz,vel,dcoeff,reac,source,bc_codes,bcvals,mgridsterm,resnorm)

	real*8, intent(in)     :: dt
	integer, intent(in)    :: lx,ly,lz

	logical, intent(in)    :: timederivflag
	
	character(LEN=4),intent(in)            :: bc_codes(NFACES)
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

	real*8,  intent(inout) :: solnvar (-g_nglayers+1:lx+g_nglayers, &
				           -g_nglayers+1:ly+g_nglayers, &
				           -g_nglayers+1:lz+g_nglayers)

	real*8, intent(inout)  :: mgridsterm(lx,ly,lz)
	
	real*8, intent(out) :: resnorm

	!local variables
	integer :: nsmoothsteps,smoothsteps
	logical :: applybcflag
	real*8  :: resnorm_2h
	
	real*8 :: res(lx,ly,lz)
	
        real*8 :: restrictres(lx/2,ly/2,lz/2)

	real*8 :: restrictsolnold(-g_nglayers+1:lx/2+g_nglayers, &
				  -g_nglayers+1:ly/2+g_nglayers, &
				  -g_nglayers+1:lz/2+g_nglayers)

	real*8 :: restrictsolnnew(-g_nglayers+1:lx/2+g_nglayers, &
				  -g_nglayers+1:ly/2+g_nglayers, &
				  -g_nglayers+1:lz/2+g_nglayers)
	
	real*8 :: error_2h       (-g_nglayers+1:lx/2+g_nglayers, &
				  -g_nglayers+1:ly/2+g_nglayers, &
				  -g_nglayers+1:lz/2+g_nglayers)

	real*8 :: error_h        (-g_nglayers+1:lx+g_nglayers, &
				  -g_nglayers+1:ly+g_nglayers, &
				  -g_nglayers+1:lz+g_nglayers)


	real*8 :: restrictsterm (lx/2,ly/2,lz/2)
	
	!coarse grid transport properties
	real*8   :: vel_c   (-g_nglayers+1:lx/2+g_nglayers, &
			     -g_nglayers+1:ly/2+g_nglayers, &
	     		     -g_nglayers+1:lz/2+g_nglayers,NDIM)

	real*8   :: dcoeff_c(-g_nglayers+1:lx/2+g_nglayers,&
			     -g_nglayers+1:ly/2+g_nglayers,&
			     -g_nglayers+1:lz/2+g_nglayers)

	real*8   :: reac_c  (-g_nglayers+1:lx/2+g_nglayers,&
			     -g_nglayers+1:ly/2+g_nglayers,&
			     -g_nglayers+1:lz/2+g_nglayers)

	real*8   :: source_c(-g_nglayers+1:lx/2+g_nglayers,&
			     -g_nglayers+1:ly/2+g_nglayers,&
			     -g_nglayers+1:lz/2+g_nglayers)

	type(boundarycondition) :: bcvals_c
	type(boundarycondition) :: bottomsolve_bcvals
	
	real*8 :: dummy1(lx/2,ly/2,lz/2)
	real*8 :: dummy2(lx/2,ly/2,lz/2)
	real*8 :: AX(lx/2,ly/2,lz/2)
	logical :: printstuff

	!initialize
	nsmoothsteps=2
	res             = ZERO
	restrictres     = ZERO
	restrictsolnold = ZERO
	restrictsolnnew = ZERO

	error_2h        = ZERO
	error_h         = ZERO
	restrictsterm   = ZERO

	dummy1          = ZERO
	dummy2          = ZERO
	AX              = ZERO

	vel_c           = ZERO
	dcoeff_c        = ZERO
	reac_c          = ZERO
	source_c        = ZERO

	printstuff=.false.

	!steps

	!1) perform smoothing
	!2) restrict fine solution (u_2h)
	!3) restrict fine residual (r_2h)
	!4) solve coarse grid problem A_2h (v) = r_2h + A_2h (u_2h)
	!5) find coarse grid error e_2h = v-u_2h
	!6) prolong error to fine grid e_h
	!7) correct fine grid soluion u_h = u_h + e_h
	!8) post smooth


	if((lx .eq. g_lx) .and. (ly .eq. g_ly) .and. (lz .eq. g_lz)) then
		applybcflag=.true.
	else
		applybcflag=.false.
	endif
	
	if((lx .le. LMIN) .or. (ly .le. LMIN) .or. (lz .le. LMIN)) then

		call allocate_bcvalues(bottomsolve_bcvals,lx*g_px,ly*g_py,lz*g_pz)
	
		call bottomsolve(solnvar,timederivflag,dt,lx,ly,lz,&
				vel,dcoeff,reac,source,bc_codes,bottomsolve_bcvals,&
				applybcflag,mgridsterm,resnorm)
	else

		call allocate_bcvalues(bcvals_c,lx/2,ly/2,lz/2)
	
		call perform_gseidel_smoothing(solnvar,timederivflag,dt,lx,ly,lz,applybcflag,&
			mgridsterm,res,vel,dcoeff,reac,source,&
			bc_codes,bcvals,&
			g_lrank,g_rrank,g_brank,g_trank,g_krank,g_frank,&
			g_llenx,g_lleny,g_llenz,nsmoothsteps)
		
		call restriction(solnvar(1:lx,1:ly,1:lz),lx,ly,lz,restrictsolnold(1:lx/2,1:ly/2,1:lz/2))
		
		call restriction(res,lx,ly,lz,restrictres)
		
		call restriction(vel(1:lx,1:ly,1:lz,XDIR),lx,ly,lz,vel_c(1:lx/2,1:ly/2,1:lz/2,XDIR))
		call restriction(vel(1:lx,1:ly,1:lz,YDIR),lx,ly,lz,vel_c(1:lx/2,1:ly/2,1:lz/2,YDIR))
		call restriction(vel(1:lx,1:ly,1:lz,ZDIR),lx,ly,lz,vel_c(1:lx/2,1:ly/2,1:lz/2,ZDIR))
		
		call restriction(dcoeff(1:lx,1:ly,1:lz),lx,ly,lz,dcoeff_c(1:lx/2,1:ly/2,1:lz/2))
		call restriction(reac  (1:lx,1:ly,1:lz),lx,ly,lz,reac_c  (1:lx/2,1:ly/2,1:lz/2))
		call restriction(source(1:lx,1:ly,1:lz),lx,ly,lz,source_c(1:lx/2,1:ly/2,1:lz/2))
		
		!need to find restricted boundary values (may not be necessary!)
		call restrict_boundary(bcvals,lx,ly,lz,bcvals_c)	

		call find_AX(restrictsolnold,timederivflag,dt,lx/2,ly/2,lz/2,vel_c,dcoeff_c,reac_c,bc_codes,bcvals,&
				AX,dummy1,dummy2,&
				g_lrank,g_rrank,g_brank,g_trank,g_krank,&
				g_frank,g_llenx,g_lleny,g_llenz)

		restrictsterm = AX + restrictres
	
		call update_boundary_ghostvalues(restrictsolnold,bc_codes,bcvals_c,lx/2,ly/2,lz/2)
		
		restrictsolnnew = restrictsolnold
		call perform_vcycle(restrictsolnnew,timederivflag,dt,&
				lx/2,ly/2,lz/2,vel_c,dcoeff_c,reac_c,source_c,&
				bc_codes,bcvals_c,&
				restrictsterm,resnorm_2h)
		
		call update_boundary_ghostvalues(restrictsolnnew,bc_codes,bcvals_c,lx/2,ly/2,lz/2)

		error_2h = restrictsolnnew - restrictsolnold

		!call prolongation0(error_2h(1:lx/2,1:ly/2,1:lz/2),lx/2,ly/2,lz/2,error_h(1:lx,1:ly,1:lz))
		call prolongation1(error_2h,lx/2,ly/2,lz/2,error_h)
		
		solnvar = solnvar + error_h
		
		call perform_gseidel_smoothing(solnvar,timederivflag,dt,lx,ly,lz,applybcflag,&
			mgridsterm,res,vel,dcoeff,reac,source,&
			bc_codes,bcvals,&
			g_lrank,g_rrank,g_brank,g_trank,g_krank,g_frank,&
			g_llenx,g_lleny,g_llenz,nsmoothsteps)
		
		call compute_norm(res,lx,ly,lz,resnorm)

	endif

end subroutine perform_vcycle
!=============================================================================================================================

end module mgridsteps
