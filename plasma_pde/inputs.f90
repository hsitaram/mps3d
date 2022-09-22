module inputs
	use globalvars
implicit none

contains
!============================================================
      subroutine readgeominputfile(inputfilename)
	
	      character(LEN=*), intent(in) :: inputfilename
	      
	      character(LEN=MAXSTRSIZE) :: temp
	      integer :: fnum

	      fnum=14

	      open(unit=fnum,file=inputfilename)

	      !GEOMETRY==========================================
	      read(fnum,*) temp
	      read(fnum,*) temp
	      read(fnum,*) temp
	      
	      read(fnum,*) temp,g_xorigin
	      read(fnum,*) temp,g_yorigin
	      read(fnum,*) temp,g_zorigin
	      read(fnum,*) temp,g_xlen
	      read(fnum,*) temp,g_ylen
	      read(fnum,*) temp,g_zlen
	      !====================================================
	      
	      !GRID AND PROCESSORS================================
	      read(fnum,*) temp
	      read(fnum,*) temp
	      read(fnum,*) temp
	      
	      read(fnum,*) temp,g_Nx
	      read(fnum,*) temp,g_Ny
	      read(fnum,*) temp,g_Nz
	      read(fnum,*) temp,g_px
	      read(fnum,*) temp,g_py
	      read(fnum,*) temp,g_pz
	      read(fnum,*) temp
	      !====================================================

	      close(fnum)

      end subroutine readgeominputfile
!============================================================

end module inputs
