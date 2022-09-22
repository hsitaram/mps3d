module globalvars

implicit none
include "mpif.h"

!all global variable names are g_
!distinguish from local variables
	integer, parameter :: MAXSTRSIZE = 300
	integer, parameter :: DOUBLESIZE = 8
	integer, parameter :: INTSIZE    = 4
	integer, parameter :: NFACES     = 6
	integer, parameter :: LEFT       = 1
	integer, parameter :: RIGHT      = 2
	integer, parameter :: BOTTOM     = 3
	integer, parameter :: TOP        = 4
	integer, parameter :: BACK       = 5
	integer, parameter :: FRONT      = 6
	integer, parameter :: DUMMYVAL   = -300000
	real*8, parameter  :: NEARZERO   = 1.d-25
	real*8, parameter  :: ZERO 	 = 0.d0
	real*8, parameter  :: TWO        = 2.d0
	real*8, parameter  :: HALF	 = 0.5d0
	integer, parameter :: LMIN       = 1
	integer, parameter :: NDIM       = 3
	integer, parameter :: XDIR       = 1
	integer, parameter :: YDIR       = 2
	integer, parameter :: ZDIR       = 3
	real*8, parameter  :: JUNKVAL    = -6666666.d0
	real*8, parameter  :: ERRTOLR    = 1.d-11 
	


	integer :: g_Nx,g_Ny,g_Nz !number of cells along each dir
	integer :: g_px,g_py,g_pz !number of procs along each dir
	integer :: g_lx,g_ly,g_lz !number of cells per proc
	integer :: g_pindx,g_pindy,g_pindz !x,y,z coordinates of procs (starting from 1)
	integer :: g_nglayers
	integer :: g_nprocs
	integer :: g_rootproc
	
	integer :: g_lrank,g_rrank !along x
	integer :: g_trank,g_brank !along y
	integer :: g_frank,g_krank !along z
	integer :: g_myproc
	integer :: g_comm
	integer :: g_nvalidranks
	
	real*8 ::  g_xlen,g_ylen,g_zlen
	real*8 ::  g_xorigin,g_yorigin,g_zorigin
	real*8 ::  g_dx,g_dy,g_dz
	real*8 ::  g_offx,g_offy,g_offz
	real*8 ::  g_llenx,g_lleny,g_llenz !local lengths for each processor

	!some parameters from input file, specific to certain problems
	real*8 :: g_prob_specific_params(20)	
	
	type boundarycondition

		real*8,allocatable :: left(:,:)
		real*8,allocatable :: right(:,:)
		real*8,allocatable :: bottom(:,:)
		real*8,allocatable :: top(:,:)
		real*8,allocatable :: back(:,:)
		real*8,allocatable :: front(:,:)

	end type boundarycondition

	type flexarray3d
		real*8,allocatable :: arr(:,:,:)
	end type flexarray3d
contains
!==================================================================================
	subroutine find3dindex(ind,npx,npy,npz,i,j,k)

		integer, intent(in) :: ind,npx,npy,npz
		integer, intent(out) :: i,j,k

		integer :: cindex !index in c/c++

		cindex = ind-1

		!fortran index
		i = mod(cindex,npx)+1
		j = mod(cindex/npx,npy)+1
		k = cindex/(npx*npy)+1

	end subroutine find3dindex
!=================================================================================
	subroutine find1dindex(i,j,k,npx,npy,npz,ind)

		integer,intent(in)  :: i,j,k,npx,npy,npz
		integer,intent(out) :: ind

		!fortran index
		ind=(k-1)*npx*npy+(j-1)*npx+(i-1)+1

	end subroutine find1dindex
!=================================================================================
	subroutine global_to_local_index(indx,indy,indz,lindx,lindy,lindz,proc)

		integer, intent(in) :: indx,indy,indz
		integer, intent(out) :: lindx,lindy,lindz,proc

		integer :: pindx,pindy,pindz

		pindx = ((indx-1)/g_lx)+1
		pindy = ((indy-1)/g_ly)+1
		pindz = ((indz-1)/g_lz)+1
		
		call find1dindex(pindx,pindy,pindz,g_px,g_py,g_pz,proc)

		!proc count is c index
		proc=proc-1

		lindx=indx-(pindx-1)*g_lx
		lindy=indy-(pindy-1)*g_ly
		lindz=indz-(pindz-1)*g_lz

	end subroutine global_to_local_index
!===============================================================================
	subroutine local_to_global_index(indx,indy,indz,proc,gindx,gindy,gindz)

		integer, intent(in)  :: indx,indy,indz,proc
		integer, intent(out) :: gindx,gindy,gindz

		integer :: pindx,pindy,pindz

		call find3dindex(proc+1,g_px,g_py,g_pz,pindx,pindy,pindz)

		gindx = (pindx-1)*g_lx + indx
		gindy = (pindy-1)*g_ly + indy
		gindz = (pindz-1)*g_lz + indz

	end subroutine local_to_global_index
!===============================================================================
	subroutine compute_norm(var,lx,ly,lz,norm)
		
		integer, intent(in) :: lx,ly,lz
		real*8, intent(in)  :: var(lx,ly,lz)
		real*8, intent(out) :: norm

		real*8  :: norm_local
		real*8  :: norm_global

		integer :: ierr

		norm_global = ZERO
		norm_local  = ZERO

		call compute_normsq_local(var,lx,ly,lz,norm_local)

		call MPI_Allreduce(norm_local,norm_global,1,MPI_REAL8,MPI_SUM,g_comm,ierr)

		norm_global = norm_global/(lx*ly*lz*g_nprocs)

		norm = sqrt(norm_global)


	end subroutine compute_norm
!===============================================================================
	subroutine compute_normsq_local(var,lx,ly,lz,norm)
		
		integer, intent(in) :: lx,ly,lz
		real*8, intent(in)  :: var(lx,ly,lz)
		real*8, intent(out) :: norm

		integer :: i,j,k

		norm  = ZERO

		do k=1,lz
			do j=1,ly
				do i=1,lx
					norm = norm + var(i,j,k)*var(i,j,k)
				enddo
			enddo
		enddo
	
	end subroutine compute_normsq_local
!===============================================================================
	subroutine print3darray(var,minx,miny,minz,maxx,maxy,maxz)

		integer,intent(in) :: minx,miny,minz,maxx,maxy,maxz		
		real*8,intent(in) :: var(minx:maxx,miny:maxy,minz:maxz)


		integer :: j,k

		print *,"*****************************************"
		print *,"RANK:",g_myproc
		print *
		do k=minz,maxz
			do j=miny,maxy
				print *,k,j,var(minx:maxx,j,k)
			enddo
		enddo
		print *,"*****************************************"

	end subroutine print3darray
!===============================================================================
	subroutine allocate_bcvalues(bcvals,lx,ly,lz)

		type(boundarycondition), intent(inout) :: bcvals
		integer,intent(in) :: lx,ly,lz

		allocate(bcvals%left(ly,lz))
		allocate(bcvals%right(ly,lz))
		
		allocate(bcvals%bottom(lx,lz))
		allocate(bcvals%top(lx,lz))

		allocate(bcvals%back(lx,ly))
		allocate(bcvals%front(lx,ly))

		bcvals%left    = ZERO
		bcvals%right   = ZERO
		bcvals%bottom  = ZERO
		bcvals%top     = ZERO
		bcvals%back    = ZERO
		bcvals%front   = ZERO

	end subroutine allocate_bcvalues
!===============================================================================
	subroutine printline(str)
	
	character(LEN=*), intent(in) :: str
	if(g_myproc .eq. g_rootproc) print *,str

	end subroutine printline
!===============================================================================

end module globalvars
