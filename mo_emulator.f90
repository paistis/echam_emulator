#ifdef __xlC__
@PROCESS HOT
@PROCESS XLF90(NOSIGNEDZERO)
#else
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#define SWDIV_NOCHK(a,b) ((a)/(b))
#endif

MODULE mo_emulator

  USE mo_kind,                 ONLY : wp
  USE mo_physical_constants,   ONLY : vtmpc1, cpd, grav
  USE m_gp
  USE m_gp_example

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: emulator

CONTAINS

  SUBROUTINE train()
  END SUBROUTINE train

  SUBROUTINE emulator()              

	TYPE (SparseGP) gp2
  	TYPE (DenseGP)  gpDense2

  	write(*,*) 'Emulator called'
  

  END SUBROUTINE emulator
	

END MODULE mo_emulator

