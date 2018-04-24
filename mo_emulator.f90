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
  USE mo_mpi, ONLY: p_parallel,p_parallel_io,p_io,p_bcast

  IMPLICIT NONE
 

  PRIVATE
  CLASS(BaseGP), allocatable :: gp_emu
  PUBLIC :: emulator,train,initEmulatorStream

CONTAINS

  SUBROUTINE initEmulatorStream
	write(*,*) 'Init emulator streams'
  END SUBROUTINE initEmulatorStream

  SUBROUTINE calcInversions(klev,iprofile,pblh,inv,mean)
    INTEGER, INTENT(IN)    :: klev,pblh
    real(wp), INTENT(IN)   :: iprofile(klev)
    real(wp), INTENT(OUT)  :: inv,mean
    real(wp)   		   :: iprofile_tmp(klev),profile_min,profile_max
    WRITE(*,*) 'Calculate emulator inversion'
    iprofile_tmp= iprofile(pblh:30)
    profile_min=minval(iprofile_tmp)
    profile_max=maxval(iprofile_tmp)
    inv=profile_max - profile_min
    mean = sum(iprofile_tmp)/size(iprofile_tmp)
  END SUBROUTINE calcInversions


  SUBROUTINE calcEmuInputs(klev,pnum,tpot_emu,pbl_emu,clw_emu,sphum_emu,emuInputVec)
     INTEGER, INTENT(IN)    :: klev,pnum
     REAL(wp), INTENT(IN)   :: tpot_emu(klev),clw_emu(klev),sphum_emu(klev),pbl_emu
     REAL(wp), INTENT(OUT)  :: emuInputVec(pnum)
     REAL(wp)		    :: tpot_mean,tpot_inv,h2o_mean,h2o_inv,h2o(klev)
     h2o = clw_emu+sphum_emu
     CALL calcInversions(klev,tpot_emu,int(pbl_emu),tpot_inv,tpot_mean)
     CALL calcInversions(klev,h2o,int(pbl_emu),h2o_inv,h2o_mean)
     emuInputVec(1)=h2o_inv
     emuInputVec(2)=tpot_inv
     emuInputVec(3)=h2o(int(pbl_emu))
     emuInputVec(4)=tpot_emu(int(pbl_emu))
     emuInputVec(5)=1000000
     WRITE(*,*) 'Calculate emulator inputs'
  END SUBROUTINE calcEmuInputs

  SUBROUTINE train()
  	character(len=100)          :: label
  	character(len=*), parameter :: filename = 'gp_test.out'
  	integer u

  	open(newunit=u, file=filename)
  	read (u,'(A)') label
  	close(u)
  
  	select case (label)
  		case('SparseGP')
     			allocate(gp_emu, source = SparseGP(filename))
  		case('DenseGP')
     			allocate(gp, source = DenseGP(filename))
  		case default
     			print *, "Incompatible data file (unrecognised GP type '", label, "')"
     			stop 1
		end select
     	IF (p_parallel) THEN ! In parallel mode        
       		CALL p_bcast(gp_emu,p_io)
     !WRITE(*,*) 'Trainign emulator'
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
    REAL(wp) :: emuInputVec(5)
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
	    CALL calcEmuInputs(klev,pnum,tpotfull(jl,jk),pbl(jk),pxlm1(jl,jk),pqm1(jl,jk),emuInputVec)
	    !if in correct region calc emulator input values
            paclc_emul(jl,jk) = gp_emu%predict(emuInputVec,pnum)
          END IF
        ENDDO 
      END IF 

    ENDDO

!   write(*,*) 'EMULATOR test'
  
  END SUBROUTINE emulator
	
END MODULE mo_emulator
