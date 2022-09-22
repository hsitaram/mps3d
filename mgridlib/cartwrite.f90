module cartwrite

use globalvars
implicit none

contains
!===============================================================================================
subroutine Writeserialvtrfile(filename,solnvars,nscalars,scalarnames,&
				lx,ly,lz,lenx,leny,lenz,offx,offy,offz)

character(LEN=*), intent(in) :: filename
integer, intent(in) :: lx,ly,lz
integer, intent(in) :: nscalars
character(LEN=*), intent(in) :: scalarnames(nscalars)

real*8, intent(in) :: solnvars(-g_nglayers+1:lx+g_nglayers,&
	-g_nglayers+1:ly+g_nglayers,-g_nglayers+1:lz+g_nglayers,nscalars)
real*8, intent(in) :: lenx,leny,lenz
real*8, intent(in) :: offx,offy,offz

character(LEN=MAXSTRSIZE) :: fname
character(LEN=MAXSTRSIZE) :: str
character(LEN=8) :: offset
character(LEN=8) :: xmin_str
character(LEN=8) :: ymin_str
character(LEN=8) :: zmin_str
character(LEN=8) :: xmax_str
character(LEN=8) :: ymax_str
character(LEN=8) :: zmax_str
character(LEN=MAXSTRSIZE) :: global_extent
character(LEN=MAXSTRSIZE) :: local_extent
character :: slashn

integer :: ierr
integer :: filenum
integer :: i,j,k
integer :: byte_offset
integer :: scalars_size

real*8,allocatable :: local_xc(:),local_yc(:),local_zc(:)
real*8 :: dx,dy,dz
integer :: Nx,Ny,Nz

filenum = 9
slashn  = char(10)

fname = trim(filename)//'.vtr'

allocate(local_xc(lx+1))
allocate(local_yc(ly+1))
allocate(local_zc(lz+1))
	
dx = lenx/lx
dy = leny/ly
dz = lenz/lz

Nx = lx
Ny = ly
Nz = lz

do i=1,lx+1
	local_xc(i) = offx + (i-1)*dx
enddo

do i=1,ly+1
	local_yc(i) = offy + (i-1)*dy
enddo

do i=1,lz+1
	local_zc(i) = offz + (i-1)*dz
enddo

scalars_size   = lx*ly*lz * DOUBLESIZE

write(xmin_str(1:8),'(i8)') 0
write(ymin_str(1:8),'(i8)') 0
write(zmin_str(1:8),'(i8)') 0
write(xmax_str(1:8),'(i8)') Nx
write(ymax_str(1:8),'(i8)') Ny
write(zmax_str(1:8),'(i8)') Nz

global_extent = trim(xmin_str)//' '//trim(xmax_str)//' '//trim(ymin_str)//' '&
		//trim(ymax_str)//trim(zmin_str)//' '//trim(zmax_str)

write(xmin_str(1:8),'(i8)') 0
write(ymin_str(1:8),'(i8)') 0
write(zmin_str(1:8),'(i8)') 0
write(xmax_str(1:8),'(i8)') lx
write(ymax_str(1:8),'(i8)') ly
write(zmax_str(1:8),'(i8)') lz

local_extent = trim(xmin_str)//' '//trim(xmax_str)//' '//trim(ymin_str)//' '&
		//trim(ymax_str)//trim(zmin_str)//' '//trim(zmax_str)

byte_offset = 0                             
open(unit=filenum,file=fname,form='unformatted',access='stream')

str = '<?xml version="1.0"?>'//slashn                                                                                                   
write(filenum) trim(str)

str = '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'//slashn		                               
write(filenum) trim(str)

str = '  <RectilinearGrid WholeExtent="'//trim(global_extent)//'">'//slashn                                       
write(filenum) trim(str)

str = '    <Piece Extent="'//trim(local_extent)//'">'//slashn
write(filenum) trim(str)

str = '      <PointData>  </PointData>'//slashn
write(filenum) trim(str)

str = '      <CellData> '//slashn                                                       
write(filenum) trim(str)

do i=1,nscalars

	write(offset(1:8),'(i8)') byte_offset
	str = '         <DataArray type="Float64" Name="'//trim(scalarnames(i))//&
		'" format="appended" offset="'//offset//'"           />'//slashn
	write(filenum) trim(str)
	byte_offset = byte_offset + intsize + scalars_size

enddo

str = '      </CellData>'//slashn            
write(filenum) trim(str)

str = '      <Coordinates>'//slashn                                                                                   
write(filenum) trim(str)

write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="X"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + intsize + (lx+1)*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="Y"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + intsize + (ly+1)*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="Z"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

str = '      </Coordinates>'//slashn       
write(filenum) trim(str)

str = '    </Piece>'//slashn                                                                                                  
write(filenum) trim(str)

str = '  </RectilinearGrid>'//slashn                                                                                                   
write(filenum) trim(str)

str = '  <AppendedData encoding="raw">'//slashn                                                                                         
write(filenum) trim(str)

str = '_'                                                                                                                           
write(filenum) trim(str)

do i=1,nscalars
write(filenum) scalars_size  ,  solnvars(1:lx,1:ly,1:lz,i)
enddo

write(filenum) (lx+1)*DOUBLESIZE  , local_xc(:)
write(filenum) (ly+1)*DOUBLESIZE  , local_yc(:)
write(filenum) (lz+1)*DOUBLESIZE  , local_zc(:)

str = slashn//'  </AppendedData>'//slashn                                                                             
write(filenum) trim(str)

str = '</VTKFile>'//slashn                                                                                                             
write(filenum) trim(str)

close(filenum)

end subroutine Writeserialvtrfile
!===============================================================================================
subroutine Writevtrfile(dirname,solnvars,nscalars,scalarnames,lx,ly,lz)

character(LEN=*), intent(in) :: dirname
integer, intent(in) :: lx,ly,lz
integer, intent(in) :: nscalars
character(LEN=*), intent(in) :: scalarnames(nscalars)
real*8, intent(in) :: solnvars(-g_nglayers+1:lx+g_nglayers,&
	-g_nglayers+1:ly+g_nglayers,-g_nglayers+1:lz+g_nglayers,nscalars)

character(LEN=MAXSTRSIZE) :: fname
character(LEN=MAXSTRSIZE) :: str
character(LEN=8) :: offset
character(LEN=8) :: xmin_str
character(LEN=8) :: ymin_str
character(LEN=8) :: zmin_str
character(LEN=8) :: xmax_str
character(LEN=8) :: ymax_str
character(LEN=8) :: zmax_str
character(LEN=8) :: proc_str
character(LEN=MAXSTRSIZE) :: global_extent
character(LEN=MAXSTRSIZE) :: local_extent
character :: slashn

integer :: ierr
integer :: filenum
integer :: i,j,k
integer :: byte_offset
integer :: scalars_size

real*8,allocatable :: local_xc(:),local_yc(:),local_zc(:)
real*8 :: dx,dy,dz

filenum = 9
slashn  = char(10)

write(proc_str,'(i8.8)') g_myproc
fname = dirname//'/file_'//trim(proc_str)//'.vtr'

allocate(local_xc(lx+1))
allocate(local_yc(ly+1))
allocate(local_zc(lz+1))
	
dx = g_llenx/lx
dy = g_lleny/ly
dz = g_llenz/lz

do i=1,lx+1
	local_xc(i) = g_offx + (i-1)*dx
enddo

do i=1,ly+1
	local_yc(i) = g_offy + (i-1)*dy
enddo

do i=1,lz+1
	local_zc(i) = g_offz + (i-1)*dz
enddo

scalars_size   = lx*ly*lz * DOUBLESIZE
	
write(xmin_str(1:8),'(i8)') lx*(g_pindx-1)
write(ymin_str(1:8),'(i8)') ly*(g_pindy-1)
write(zmin_str(1:8),'(i8)') lz*(g_pindz-1)
write(xmax_str(1:8),'(i8)') lx*g_pindx
write(ymax_str(1:8),'(i8)') ly*g_pindy
write(zmax_str(1:8),'(i8)') lz*g_pindz

local_extent = trim(xmin_str)//' '//trim(xmax_str)//' '//trim(ymin_str)//' '&
		//trim(ymax_str)//trim(zmin_str)//' '//trim(zmax_str)

byte_offset = 0                             
open(unit=filenum,file=fname,form='unformatted',access='stream')

str = '<?xml version="1.0"?>'//slashn                                                                                                   
write(filenum) trim(str)

str = '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'//slashn		                               
write(filenum) trim(str)

str = '  <RectilinearGrid WholeExtent="'//trim(local_extent)//'">'//slashn                                       
write(filenum) trim(str)

str = '    <Piece Extent="'//trim(local_extent)//'">'//slashn
write(filenum) trim(str)

str = '      <PointData>  </PointData>'//slashn
write(filenum) trim(str)

str = '      <CellData> '//slashn                                                       
write(filenum) trim(str)

do i=1,nscalars

	write(offset(1:8),'(i8)') byte_offset
	str = '         <DataArray type="Float64" Name="'//trim(scalarnames(i))//&
		'" format="appended" offset="'//offset//'"           />'//slashn
	write(filenum) trim(str)
	byte_offset = byte_offset + intsize + scalars_size

enddo

str = '      </CellData>'//slashn            
write(filenum) trim(str)

str = '      <Coordinates>'//slashn                                                                                   
write(filenum) trim(str)

write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="X"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + intsize + (lx+1)*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="Y"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + intsize + (ly+1)*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="Z"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

str = '      </Coordinates>'//slashn       
write(filenum) trim(str)

str = '    </Piece>'//slashn                                                                                                  
write(filenum) trim(str)

str = '  </RectilinearGrid>'//slashn                                                                                                   
write(filenum) trim(str)

str = '  <AppendedData encoding="raw">'//slashn                                                                                         
write(filenum) trim(str)

str = '_'                                                                                                                           
write(filenum) trim(str)

do i=1,nscalars
write(filenum) scalars_size  ,  solnvars(1:lx,1:ly,1:lz,i)
enddo

write(filenum) (lx+1)*DOUBLESIZE  , local_xc(:)
write(filenum) (ly+1)*DOUBLESIZE  , local_yc(:)
write(filenum) (lz+1)*DOUBLESIZE  , local_zc(:)

str = slashn//'  </AppendedData>'//slashn                                                                             
write(filenum) trim(str)

str = '</VTKFile>'//slashn                                                                                                             
write(filenum) trim(str)

close(filenum)

end subroutine Writevtrfile
!===============================================================================================
subroutine Writepvtrfile(dirname,nscalars,scalarnames,lx,ly,lz)

	character(LEN=*), intent(in) :: dirname
	integer, intent(in) :: lx,ly,lz
	integer,intent(in)  :: nscalars
	character(LEN=*) :: scalarnames(nscalars)

	character(LEN=MAXSTRSIZE) :: fname
	character(LEN=MAXSTRSIZE) :: vtrfname
	character(LEN=MAXSTRSIZE*5) :: str
	character(LEN=8) :: xmin_str
	character(LEN=8) :: ymin_str
	character(LEN=8) :: zmin_str
	character(LEN=8) :: xmax_str
	character(LEN=8) :: ymax_str
	character(LEN=8) :: zmax_str
	character(LEN=8) :: proc_str
	character(LEN=MAXSTRSIZE) :: global_extent
	character(LEN=MAXSTRSIZE) :: local_extent
	character :: slashn

	integer :: ierr
	integer :: filenum
	integer :: i
	integer :: pind_x,pind_y,pind_z
	integer :: Nx,Ny,Nz

	Nx = lx*g_px
	Ny = ly*g_py
	Nz = lz*g_pz

	call system('mkdir '//dirname)

	filenum = 9
	slashn  = char(10)
	fname = dirname//'.pvtr'

	write(xmin_str(1:8),'(i8)') 0
	write(ymin_str(1:8),'(i8)') 0
	write(zmin_str(1:8),'(i8)') 0
	write(xmax_str(1:8),'(i8)') Nx
	write(ymax_str(1:8),'(i8)') Ny
	write(zmax_str(1:8),'(i8)') Nz

	global_extent = trim(xmin_str)//' '//trim(xmax_str)//' '//trim(ymin_str)//' '&
		//trim(ymax_str)//trim(zmin_str)//' '//trim(zmax_str)

	open(unit=filenum,file=fname,form='unformatted',access='stream')

	str = '<?xml version="1.0"?>'//slashn                                                                                                   
	write(filenum) trim(str)

	str = '<VTKFile type="PRectilinearGrid" version="0.1" byte_order="LittleEndian">'//slashn		                               
	write(filenum) trim(str)

	str = '<PRectilinearGrid WholeExtent="'//trim(global_extent)//'" GhostLevel="0">'//slashn                                       
	write(filenum) trim(str)


	str = '<PPointData>  </PPointData>'//slashn
	write(filenum) trim(str)

	str = '<PCellData> '//slashn                                                       
	write(filenum) trim(str)

	do i=1,nscalars

		str = '<PDataArray type="Float64" Name="'//trim(scalarnames(i))//'"/>'//slashn
		write(filenum) trim(str)

	enddo

	str = '</PCellData>'//slashn            
	write(filenum) trim(str)

	str = '<PCoordinates>'//slashn                                                                                   
	write(filenum) trim(str)

	str = '<PDataArray type="Float64" Name="X"/>'//slashn
	write(filenum) trim(str)
	str = '<PDataArray type="Float64" Name="Y"/>'//slashn
	write(filenum) trim(str)
	str = '<PDataArray type="Float64" Name="Z"/>'//slashn
	write(filenum) trim(str)

	str = '</PCoordinates>'//slashn       
	write(filenum) trim(str)


	do i=1,g_nprocs
		
	    call find3dindex(i,g_px,g_py,g_pz,pind_x,pind_y,pind_z)
	    
	    write(xmin_str(1:8),'(i8)') lx*(pind_x-1)
            write(ymin_str(1:8),'(i8)') ly*(pind_y-1)
            write(zmin_str(1:8),'(i8)') lz*(pind_z-1)
            write(xmax_str(1:8),'(i8)') lx*pind_x
            write(ymax_str(1:8),'(i8)') ly*pind_y
            write(zmax_str(1:8),'(i8)') lz*pind_z

	    local_extent = trim(xmin_str)//' '//trim(xmax_str)//&
			    ' '//trim(ymin_str)//' '//trim(ymax_str)//&
			    trim(zmin_str)//' '//trim(zmax_str)

	    write(proc_str,'(i8.8)') i-1
            vtrfname = dirname//'/file_'//trim(proc_str)//'.vtr'


	    str = '<Piece Extent="'//trim(local_extent)//&
			    '" Source="'//trim(vtrfname)//'"/>'//slashn    
	    !print *,str
	    write(filenum) trim(str)

	 enddo

	str = '</PRectilinearGrid>'//slashn                                                                 
	write(filenum) trim(str)

	str = '</VTKFile>'//slashn                                           
	write(filenum) trim(str)

	close(filenum)

end subroutine Writepvtrfile
!===============================================================================================
subroutine writeoutput(dirname,solnvars,nscalars,scalarnames,lx,ly,lz)

	character(LEN=*), intent(in) :: dirname
	integer, intent(in) :: lx,ly,lz
	integer, intent(in) :: nscalars
	
	character(LEN=*), intent(in) :: scalarnames(nscalars)
	real*8, intent(in) :: solnvars(-g_nglayers+1:lx+g_nglayers, &
				       -g_nglayers+1:ly+g_nglayers, &
				       -g_nglayers+1:lz+g_nglayers, &
				        nscalars)

	integer :: ierr
	integer :: i,j,k


	if(g_myproc .eq. g_rootproc) then
		call Writepvtrfile(dirname,nscalars,scalarnames,lx,ly,lz)
	endif
	call MPI_Barrier(g_comm,ierr)
	
        call Writevtrfile(dirname,solnvars,nscalars,scalarnames,lx,ly,lz) 	

end subroutine writeoutput
!===============================================================================================

end module cartwrite

