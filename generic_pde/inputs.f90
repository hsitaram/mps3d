module inputs
use globalvars
implicit none

contains
!============================================================
      subroutine readinputfile(inputfilename,bc_codes,bcparams)
	
	      integer :: fnum
	      character(LEN=*), intent(in) :: inputfilename
	      character(LEN=MAXSTRSIZE) :: temp
	      
	      character(LEN=4),intent(out) :: bc_codes(NFACES)
	      real*8, intent(out) :: bcparams(NFACES)

	      fnum=14

	      open(unit=fnum,file=inputfilename)
	      read(fnum,*) temp
	      read(fnum,*) temp
	      read(fnum,*) temp
	      
	      read(fnum,*) temp,g_xorigin
	      read(fnum,*) temp,g_yorigin
	      read(fnum,*) temp,g_zorigin
	      read(fnum,*) temp,g_xlen
	      read(fnum,*) temp,g_ylen
	      read(fnum,*) temp,g_zlen
	      
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
	      read(fnum,*) temp
	      read(fnum,*) temp
	      
	      read(fnum,*) temp,bc_codes(LEFT) ,  bcparams(LEFT)
	      read(fnum,*) temp,bc_codes(RIGHT),  bcparams(RIGHT)
	      read(fnum,*) temp,bc_codes(BOTTOM), bcparams(BOTTOM)
	      read(fnum,*) temp,bc_codes(TOP),    bcparams(TOP)
	      read(fnum,*) temp,bc_codes(BACK),   bcparams(BACK)
	      read(fnum,*) temp,bc_codes(FRONT),  bcparams(FRONT)

	      g_dx = g_xlen/g_Nx
	      g_dy = g_ylen/g_Ny
	      g_dz = g_zlen/g_Nz

      end subroutine readinputfile
!============================================================

end module inputs
