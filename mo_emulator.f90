!#ifdef __xlC__
!@PROCESS HOT
!@PROCESS XLF90(NOSIGNEDZERO)
!#else
!#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
!#define SWDIV_NOCHK(a,b) ((a)/(b))
!#endif

MODULE mo_emulator

  !USE mo_kind,                 ONLY : wp
  !USE mo_physical_constants,   ONLY : vtmpc1, cpd, grav
  USE m_gp
  USE m_gp_example

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: emulator

CONTAINS

  SUBROUTINE train()
	!Subroutine used for training the emulator
	!read input variables from netcdf file
	file="emulator_input.nc"
        iret = nf_open(file, NF_NOWRITE, ncid)
        IF (iret /= NF_NOERR) STOP 'Emulator NetCDF File not opened'
	

  END SUBROUTINE train

  SUBROUTINE emulator()              
	!This function is called from echam physics, should return emulated fields, e.g cloud cover
	TYPE (SparseGP) gp2
  	TYPE (DenseGP)  gpDense2

  	write(*,*) 'Emulator called'
  

  END SUBROUTINE emulator
	

END MODULE mo_emulator

