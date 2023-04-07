module par_decompose

	use globalvars	
	implicit none
contains
!=====================================================================
	subroutine domaindecompose()

		integer :: ierr
		integer :: maxprocs

		g_comm = MPI_COMM_WORLD
		g_nprocs = g_px*g_py*g_pz
		g_rootproc = 0
		g_lx = g_Nx/g_px;
		g_ly = g_Ny/g_py;
		g_lz = g_Nz/g_pz;
		g_nglayers = 1;
	      
		g_dx = g_xlen/g_Nx
	      	g_dy = g_ylen/g_Ny
	      	g_dz = g_zlen/g_Nz


		call MPI_Comm_rank(g_comm,g_myproc,ierr)
		call MPI_Comm_size(g_comm,maxprocs,ierr)

		if(maxprocs .ne. g_nprocs) then
			print *,"processor count does not match MPI call"
			call MPI_Abort(g_comm,ierr)
		endif
			
		call find3dindex(g_myproc+1,g_px,g_py,g_pz,g_pindx, &
				g_pindy,g_pindz)

		g_offx = g_xorigin + (g_pindx-1)*g_lx*g_dx
		g_offy = g_yorigin + (g_pindy-1)*g_ly*g_dy
		g_offz = g_zorigin + (g_pindz-1)*g_lz*g_dz

		g_llenx = g_lx*g_dx
		g_lleny = g_ly*g_dy
		g_llenz = g_lz*g_dz

                g_rrank=g_myproc + 1
		g_lrank=g_myproc - 1
                
		g_trank=g_myproc + g_px
                g_brank=g_myproc - g_px
                
                g_frank=g_myproc + g_px*g_py
		g_krank=g_myproc - g_px*g_py

    		g_nvalidranks=6

    		if(g_pindx==1) then
       			g_lrank=MPI_PROC_NULL
       			g_nvalidranks=g_nvalidranks-1
    		endif
    		if(g_pindx==g_px) then
       			g_rrank=MPI_PROC_NULL
       			g_nvalidranks=g_nvalidranks-1
    		endif
    		
		if(g_pindy==1) then
       			g_brank=MPI_PROC_NULL
       			g_nvalidranks=g_nvalidranks-1
    		endif
    		if(g_pindy==g_py) then
       			g_trank=MPI_PROC_NULL
       			g_nvalidranks=g_nvalidranks-1
    		endif

    		if(g_pindz==1) then
       			g_krank=MPI_PROC_NULL
       			g_nvalidranks=g_nvalidranks-1
    		endif
    		if(g_pindz==g_pz) then
       			g_frank=MPI_PROC_NULL
       			g_nvalidranks=g_nvalidranks-1
    		endif

		print *,"myproc:",g_myproc,"l,r,b,t,f,k:",g_lrank,&
				g_rrank,g_brank,g_trank,g_frank,g_krank

	end subroutine domaindecompose
!=======================================================================
	subroutine postsend(sendarr,dest,req)

		real*8, intent(in) :: sendarr(:)
		integer, intent(in) :: dest
		integer, intent(in) :: req
		
		integer :: tag
		integer :: ierr

		tag = g_myproc*10+dest+100

		!print *,"g_myproc:",g_myproc,"tag:",tag

		call MPI_Isend(sendarr,size(sendarr),MPI_REAL8,dest,tag,g_comm,req,ierr)


	end subroutine postsend
!=======================================================================
	subroutine postrecv(recvarr,src,req)

		real*8, intent(in) :: recvarr(:)
		integer, intent(in) :: src
		integer, intent(in) :: req
		
		integer :: tag
		integer :: ierr

		tag = src*10+g_myproc+100
		
		!print *,"g_myproc:",g_myproc,"tag:",tag

		call MPI_Irecv(recvarr,size(recvarr),MPI_REAL8,src,tag,g_comm,req,ierr)

	end subroutine postrecv
!=======================================================================
	subroutine exchangehalodata(solnvar,lx,ly,lz)
	
		integer, intent(in)   :: lx,ly,lz
		real*8, intent(inout) :: solnvar(-g_nglayers+1:lx+g_nglayers, &
						 -g_nglayers+1:ly+g_nglayers, &
						 -g_nglayers+1:lz+g_nglayers)
		

		real*8,allocatable :: recvbuf_l(:),recvbuf_r(:)
		real*8,allocatable :: recvbuf_b(:),recvbuf_t(:)
		real*8,allocatable :: recvbuf_k(:),recvbuf_f(:)

		real*8,allocatable :: sendbuf_l(:),sendbuf_r(:)
		real*8,allocatable :: sendbuf_b(:),sendbuf_t(:)
		real*8,allocatable :: sendbuf_k(:),sendbuf_f(:)

		integer :: pos,i,j,k
		integer :: req(2*NFACES),ierr
		integer :: stat(MPI_STATUS_SIZE,2*NFACES)


		allocate(recvbuf_l(ly*lz*g_nglayers))
		allocate(recvbuf_r(ly*lz*g_nglayers))

		allocate(recvbuf_b(lx*lz*g_nglayers))
		allocate(recvbuf_t(lx*lz*g_nglayers))
		
		allocate(recvbuf_k(lx*ly*g_nglayers))
		allocate(recvbuf_f(lx*ly*g_nglayers))
	


		allocate(sendbuf_l(ly*lz*g_nglayers))
		allocate(sendbuf_r(ly*lz*g_nglayers))
		
		allocate(sendbuf_b(lx*lz*g_nglayers))
		allocate(sendbuf_t(lx*lz*g_nglayers))
		
		allocate(sendbuf_k(lx*ly*g_nglayers))
		allocate(sendbuf_f(lx*ly*g_nglayers))

		recvbuf_l = 0.d0
		recvbuf_r = 0.d0
		recvbuf_b = 0.d0
		recvbuf_t = 0.d0
		recvbuf_k = 0.d0
		recvbuf_f = 0.d0

		sendbuf_l = 0.d0
		sendbuf_r = 0.d0
		sendbuf_b = 0.d0
		sendbuf_t = 0.d0
		sendbuf_k = 0.d0
		sendbuf_f = 0.d0
		
		pos=1
		do k=1,lz
			do j=1,ly
				sendbuf_l(pos)=solnvar(1,j,k)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do k=1,lz
			do j=1,ly
				sendbuf_r(pos)=solnvar(lx,j,k)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do k=1,lz
			do i=1,lx
				sendbuf_b(pos)=solnvar(i,1,k)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do k=1,lz
			do i=1,lx
				sendbuf_t(pos)=solnvar(i,ly,k)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do j=1,ly
			do i=1,lx
				sendbuf_k(pos)=solnvar(i,j,1)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do j=1,ly
			do i=1,lx
				sendbuf_f(pos)=solnvar(i,j,lz)
				pos=pos+1
			enddo
		enddo
		
		call postrecv(recvbuf_l,g_lrank,req(LEFT))	
		call postrecv(recvbuf_r,g_rrank,req(RIGHT))	
		
		call postrecv(recvbuf_b,g_brank,req(BOTTOM))	
		call postrecv(recvbuf_t,g_trank,req(TOP))	
		
		call postrecv(recvbuf_k,g_krank,req(BACK))	
		call postrecv(recvbuf_f,g_frank,req(FRONT))	
		


		call postsend(sendbuf_l,g_lrank,req(LEFT+NFACES))	
		call postsend(sendbuf_r,g_rrank,req(RIGHT+NFACES))	
		
		call postsend(sendbuf_b,g_brank,req(BOTTOM+NFACES))	
		call postsend(sendbuf_t,g_trank,req(TOP+NFACES))	
		
		call postsend(sendbuf_k,g_krank,req(BACK+NFACES))	
		call postsend(sendbuf_f,g_frank,req(FRONT+NFACES))

		call MPI_Waitall(2*NFACES,req,stat,ierr)

		!need to update ghost nodes
		pos=1
		do k=1,lz
			do j=1,ly
				solnvar(0,j,k) = recvbuf_l(pos)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do k=1,lz
			do j=1,ly
				solnvar(lx+g_nglayers,j,k) = recvbuf_r(pos)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do k=1,lz
			do i=1,lx
				solnvar(i,0,k) = recvbuf_b(pos)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do k=1,lz
			do i=1,lx
				solnvar(i,ly+g_nglayers,k) = recvbuf_t(pos)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do j=1,ly
			do i=1,lx
				solnvar(i,j,0) = recvbuf_k(pos)
				pos=pos+1
			enddo
		enddo
		
		pos=1
		do j=1,ly
			do i=1,lx
				solnvar(i,j,lz+g_nglayers) = recvbuf_f(pos)
				pos=pos+1
			enddo
		enddo

	end subroutine exchangehalodata
!=================================================================================================================
	subroutine gatherterms(term,lx,ly,lz,gatheredterm)
		
		integer, intent(in) :: lx,ly,lz
		real*8, intent(in)  :: term(lx,ly,lz)
		real*8, intent(out) :: gatheredterm(lx*g_px,ly*g_py,lz*g_pz)
	
		real*8  :: recvbuf(lx,ly,lz,g_nprocs)
		integer :: i,j,k,ierr
		integer :: offx,offy,offz
		integer :: pindx,pindy,pindz,proc

		gatheredterm=ZERO

		call MPI_gather(term,lx*ly*lz,MPI_REAL8,recvbuf,lx*ly*lz,&
			MPI_REAL8,g_rootproc,g_comm,ierr)

		call MPI_Barrier(g_comm,ierr)

		if(g_myproc .eq. g_rootproc) then

			do proc=1,g_nprocs
			
				call find3dindex(proc,g_px,g_py,g_pz,pindx,pindy,pindz)
			
				offx = (pindx-1)*lx
				offy = (pindy-1)*ly
				offz = (pindz-1)*lz

				do k=1,lz
					do j=1,ly
						do i=1,lx
							gatheredterm(i+offx,j+offy,k+offz)=recvbuf(i,j,k,proc)
						enddo
					enddo
				enddo
			enddo
		endif

	end subroutine gatherterms
!===============================================================================
	
end module par_decompose
