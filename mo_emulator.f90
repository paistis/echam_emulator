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
  type(SparseGP), allocatable :: gp_emu
  PUBLIC :: emulator,train

CONTAINS

  SUBROUTINE calcInversions(klev,iprofile,pblh,inv,mean)
    INTEGER, INTENT(IN)    :: klev,pblh
    real(wp), INTENT(IN)   :: iprofile,pblh
    real(wp), INTENT(OUT)  :: inv,mean
    WRITE(*,*) 'Calculate emulator inversion'
    iprofile= iprofile(pblh:30)
    profile_min=min(iprofile)
    profile_max=max(iprofile)
    inv=profile_max - profile_min
    mean = sum(profile)/size(profile)
  END SUBROUTINE calcInversions

  SUBROUTINE calcAVG(klev, iprofile)
    INTEGER, INTENT(IN)    :: klev
    real(wp), INTENT(IN)   :: iprofile 
     WRITE(*,*) 'Calculate mean'
  END SUBROUTINE calcAVG

  SUBROUTINE calcEmuInputs(klev,pnun,tpot_emu,pbl_emu,clw_emu,sphum_emu,emuInputVec)
     INTEGER, INTENT(IN)    :: klev,pnum
     REAL(wp), INTENT(IN)   :: tpot_emu(klev),pbl_emu(klev),clw_emu(klev),sphum_emu(klev)
     REAL(wp), INTENT(OUT)  :: emuInputVec(pnum)
     WRITE(*,*) 'Calculate emulator inputs'
  END SUBROUTINE calcEmuInputs

  SUBROUTINE train()
     allocate(gp_emu, source=SparseGP('gp_test.out'))
     WRITE(*,*) 'Trainign emulator'
  END SUBROUTINE train

  SUBROUTINE emulator(kproma, kbdim, ktdia, klev, klevp1, &
                      land_fraction, philon, philat, pfull, phalf, tfull, tpotfull, &
		      pbl,pxlm1,pqm1, &
                      lmask_emul, paclc_emul, pxl_emul, pacdnc_emul)          

! Input for the emulator (to be expanded ...)

    INTEGER, INTENT(IN)    :: kbdim, klevp1, klev, kproma, ktdia

    REAL(wp), INTENT(in) :: &
      land_fraction(kbdim), & ! Land fraction
      philon(kbdim),        & ! Longitude (degrees)
      philat(kbdim),        & ! Latitude (degrees)
      pfull(kbdim,klev),    & ! Full-level pressure (Pa)  
      phalf(kbdim,klevp1),  & ! Half-level pressure (Pa)  
      tfull(kbdim,klev),    & ! Full-level temperature (K)  
      tpotfull(kbdim,klev), & ! Full-level potential temperature (K)
      pbl(kbdim),           & ! planetary boundary layer heigh
      pxlm1(kbdim,klev),    & ! cloud liquid water content
      pqm1(kbdim,klev)       ! specific humidity

! Output from the emulator

    LOGICAL, INTENT(out) :: &
      lmask_emul(kbdim,klev)   ! Indicates grid cells where the emulator is used

    REAL(wp), INTENT(out) :: &
      paclc_emul(kbdim,klev),& ! Cloud fraction from the emulator 
      pxl_emul(kbdim,klev),  & ! Cloud liquid water [kg/kg] from the emulator  
      pacdnc_emul(kbdim,klev)  ! Cloud droplet number concentration [1/m3] 
                               ! from the emulator
! Local variables   
    INTEGER :: jl,jk
    LOGICAL :: lstratocum(kbdim)
    INTEGER :: pnum = 5            ! number of parameter used in emulator
    REAL(wp) :: emuInputVec(pnun)
! Initialize the output
    lmask_emul(1:kproma,:) = .FALSE.
    paclc_emul(1:kproma,:) = 0.0_wp
    pxl_emul(1:kproma,:) = 0.0_wp
    pacdnc_emul(1:kproma,:) = 0.0_wp

!calculate emultor inputvalues
    
     

! For testing, set LMASK_EMUL and PACLC_EMUL based on
! - surface type (longitude, latitude)
! - pressure
! - temperature?

    DO jl=1,kproma
      lstratocum(jl)=.FALSE.
! Peruvian stratocumulus region (70-105 W, 5-35 S)
      IF (((philon(jl)>255.).AND.(philon(jl)<290.))&
           .AND.((philat(jl)>-35.).AND.(philat(jl)<-5.))) lstratocum(jl)=.TRUE.
! Californian stratocumulus region (110-130 W, 20-40 N)
      IF (((philon(jl)>230.).AND.(philon(jl)<250.))&
           .AND.((philat(jl)>20.).AND.(philat(jl)<40.))) lstratocum(jl)=.TRUE.
! Namibian stratocumulus region (10W-10 E, 5-30 S)
      IF (((philon(jl)>350.).OR.(philon(jl)<10.))&
           .AND.((philat(jl)>-30.).AND.(philat(jl)<-5.))) lstratocum(jl)=.TRUE.

      IF (land_fraction(jl) > 0.5) lstratocum(jl)=.FALSE.
       
      IF (lstratocum(jl)) THEN
        DO jk=ktdia,klev  
! Use the emulator in "stratocumulus regions" for levels below 800 hPa, 
! and only for above-zero temperatures 
          IF ((pfull(jl,jk) > 80000.).AND.(tfull(jl,jk)>273.15)) THEN
            lmask_emul(jl,jk) = .TRUE.
! Set the cloud fraction to some arbitrary value
	    CALL calcEmuInputs(klev,pnun,tpot_emu,pbl_emu,clw_emu,sphum_emu,emuInputVec)
	    !if in correct region calc emulator input values
            paclc_emul(jl,jk) = gp_emu%predict(emuInputVec,pnum)
          END IF
        ENDDO 
      END IF 

    ENDDO

!   write(*,*) 'EMULATOR test'
  
  END SUBROUTINE emulator
	
END MODULE mo_emulator
