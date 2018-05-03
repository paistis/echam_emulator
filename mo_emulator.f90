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
  ! type definition
  TYPE t_emulator
    REAL(wp), POINTER          :: tpot_inv(:,:)       ! potential temperature inversion strengt [m]
    REAL(wp), POINTER          :: tpot_pbl(:,:)      ! potential temperature at pbl [m]
    REAL(wp), POINTER          :: h2o_inv(:,:)       ! h2o inversion
    REAL(wp), POINTER          :: h2o_pbl(:,:)    ! h2o mixin ration at pbl[kg]
    REAL(wp), POINTER          :: pbl_num(:,:)     ! particle number at pbl[m3]
    REAL(wp), POINTER          :: emu_cf(:,:)  ! emulated cloud fraction[m]
  END TYPE t_emulator

  ! emulator_stream
  ! emulatorvars is used to check namelist input for valid names
  INTEGER, PARAMETER           :: nemulatorvars=6! s.stadtler
  CHARACTER(LEN=32)            :: emulatorvars(1:nemulatorvars)= &
                                (/'tpot_inv            ', &
                                  'tpot_pbl          ', &
                                  'h2o_inv            ', &
                                  'h2o_pbl         ', &  
                                  'pbl_num          ', &  
                                  'emu_cf       ' /)


  ! variable pointers
  TYPE(t_emulator)               :: emulators
CONTAINS

  SUBROUTINE initEmulatorStream
    
    USE mo_control,             ONLY: nlev
    USE mo_submodel_streams,    ONLY: emulator_lpost, emulator_tinterval, emulatornam
    USE mo_string_utls,         ONLY: st1_in_st2_proof
    USE mo_util_string,         ONLY: tolower
    USE mo_exception,           ONLY: finish
    USE mo_memory_base,         ONLY: t_stream, new_stream, &
                                      default_stream_setting, &
                                      add_stream_element, &
                                      AUTO, BELOWSUR
    ! local variables
    INTEGER, PARAMETER             :: ndefault = 6
  CHARACTER(LEN=32)            :: defnam(1:nemulatorvars)= &
                                (/'tpot_inv            ', &
                                  'tpot_pbl          ', &
                                  'h2o_inv            ', &
                                  'h2o_pbl         ', &  
                                  'pbl_num          ', &  
                                  'emu_cf       ' /)
 
    TYPE (t_stream), POINTER       :: emulator_stream
    INTEGER                        :: ierr
    LOGICAL                        :: lpost
    !-- handle ALL and DEFAULT options
    IF (TRIM(tolower(emulatornam(1))) == 'all')     emulatornam(1:nemulatorvars) = emulatorvars(:)
    IF (TRIM(tolower(emulatornam(1))) == 'default') emulatornam(1:ndefault) = defnam(:)

    !-- check that all variable names from namelist are valid
    IF (.NOT. st1_in_st2_proof( emulatornam, emulatorvars, ierr=ierr) ) THEN
      IF (ierr > 0) CALL finish ( 'ini_emulator_stream', 'variable '// &
                                  emulatornam(ierr)//' does not exist in emulator stream' )
    END IF

    !-- open new stream
    !SF #383: registering emulator to rerun storage
    CALL new_stream (emulator_stream,'emulator',lpost=emulator_lpost,lrerun=.TRUE., &
         interval=emulator_tinterval)
    CALL default_stream_setting (emulator_stream, lrerun = .TRUE., &
         contnorest = .TRUE., table = 199, &
         laccu = .false., code = AUTO)

    !-- add individual variables to stream
    lpost = st1_in_st2_proof( 'tpot_inv', emulatornam)
    CALL add_stream_element (emulator_stream, 'tpot_inv', emulators%tpot_inv, &
         longname = 'Potential temperature inversion', &
         units = 'K', lpost = lpost)
    
    !lpost = st1_in_st2_proof( 'tpot_pbl', emulatornam)
    !CALL add_stream_element (emulator_stream, 'tpot_pbl', emulators%tpot_pbl, &                                            &
    !     longname = 'potential temperatuer at boundary layer', &
    !     units = 'K', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'h2o_inv', emulatornam)
    CALL add_stream_element (emulator_stream, 'h2o_inv', emulators%h2o_inv, &
         longname = 'h2o inversion', &
         units = 'kg', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'h2o_pbl', emulatornam)
    CALL add_stream_element (emulator_stream, 'h2o_pbl', emulators%h2o_pbl, &
         longname = 'h2o in plb', &
         units = 'm3', lpost = lpost)

    lpost = st1_in_st2_proof( 'pbl_num', emulatornam)
    CALL add_stream_element (emulator_stream, 'pbl_num', emulators%pbl_num, &
         longname = 'number of particles in pbl', &
         units = 'kg kg-1', lpost=lpost)

    lpost = st1_in_st2_proof( 'emu_cf', emulatornam)
    CALL add_stream_element (emulator_stream, 'emu_cf', emulators%emu_cf, &
         longname = 'cloud fraction', &
         units = '1', lpost=lpost)
  END SUBROUTINE initEmulatorStream

  SUBROUTINE calcInversions(klev,iprofile,pblh,inv,mean)
    INTEGER, INTENT(IN)    :: klev,pblh
    real(wp), INTENT(IN)   :: iprofile(klev)
    real(wp), INTENT(OUT)  :: inv,mean
    real(wp)   		   :: iprofile_tmp(klev),profile_min,profile_max
    !WRITE(*,*) 'Calculate emulator inversion'
    iprofile_tmp= iprofile(pblh:klev)
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
     !WRITE(*,*) 'Calculate emulator inputs'
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
     			allocate(gp_emu, source = DenseGP(filename))
  		case default
     			print *, "Incompatible data file (unrecognised GP type '", label, "')"
     			stop 1
		end select
     	!IF (p_parallel) THEN ! In parallel mode        
       	!	CALL p_bcast(gp_emu,p_io)
	!END IF
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
    emulators%emu_cf(:,:) = 00.0_wp
    emulators%tpot_inv(:,:) =0.0_wp
    emulators%tpot_pbl(:,:) =0.0_wp
    emulators%h2o_pbl(:,:) =  0.0_wp
    emulators%h2o_inv(:,:) = 0.0_wp
    emulators%pbl_num(:,:) =  0.0_wp
!calculate emultor inputvalues
    
     !call train()

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
	    emuInputVec(:) = 0.0_wp
	    CALL calcEmuInputs(klev,pnum,tpotfull(jl,:),pbl(jk),pxlm1(jl,:),pqm1(jl,:),emuInputVec)
	    !if in correct region calc emulator input values
            paclc_emul(jl,jk) = gp_emu%predict(emuInputVec,pnum)
	    !emulators%emu_cf(jl,jk) = paclc_emul(jl,jk)
	    !emulators%tpot_inv(jl,jk) = emuInputVec(2)
            emulators%tpot_pbl(jl,jk) = 0.0_wp ! emuInputVec(4)
            !emulators%h2o_pbl(jl,jk) =  emuInputVec(3)
            !emulators%h2o_inv(jl,jk) = emuInputVec(1)
            !emulators%pbl_num(jl,jk) =  emuInputVec(5)
          END IF
        ENDDO 
      END IF 

    ENDDO

!   write(*,*) 'EMULATOR test'
  
  END SUBROUTINE emulator
	
END MODULE mo_emulator
