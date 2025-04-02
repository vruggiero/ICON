!
!+ ASCAT soil moisture operator
!
MODULE mo_soil_obs
!
! Description:
!   ASCAT soil moisture operator
!
! Current Maintainer: Harald Anlauf (DWD), Valerio Cardinali (Comet)
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de, valerio.cardinali@aeronautica.difesa.it
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1           15-02-2021 Valerio Cardinali (Comet)
!              implementation of H operator for soil moisture
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin     DWD    2018  original code
! Valerio Cardinali  Comet  2021  implementation of H operator
!==============================================================================

  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_exception,  only: finish           ! abort routine
  use mo_kind,       only: wp, sp           ! kind parameters
  use mo_time,       only: cyyyymmddhhmmss  ! convert time to string
  use mo_mpi_dace,   only: dace,           &! MPI group info
                           p_bcast
  use mo_soil,       only: ind_wso,        &! index from soil moisture
                           n_st,           &! number of soil types
                           cporv,          &! pore volume    (fraction of volume)
!                          cfcap,          &! field capacity (fraction of volume)
!                          cpwp,           &! plant wilting point (fraction of volume)
                           cadp,           &! air dryness point (fraction of volume)
                           ST_ICE,ST_ROCK, &! ice and rock soiltypes
                           ST_SEAWATER, ST_SEAICE  ! seawater and seaice soiltypes
  use mo_namelist,   only: position_nml,   &! routine to position nml group
                           nnml,           &! namelist fortran unit number
                           POSITIONED       ! position_nml: OK    return flag
  use mo_wmo_tables, only: WMO6_LATLON,    &! regular grid types
                           WMO6_ROTLL,     &!
                           WMO6_GAUSSIAN,  &!
                           DWD6_ICON        ! ICON grid type
  !-----------------------------
  ! access observation data type
  !-----------------------------
  use mo_t_datum,    only: t_datum,        &! report body derived type
                           rvind,          &! missing value indicator (real)
                           SRC_DER          ! derived quantity flag value
  use mo_t_obs,      only: t_obs,          &! observation derived type
                           t_spot,         &! observation header derived type
                           t_head,         &! observation header derived type
                           t_mcol,         &! model column descriptor
!                          t_mcols,        &! set of model column descriptors
                           source,         &! list of report source files
                           set_xuv,        &! set unit vectors, zenith angle
                           new_spot,       &! reserve memory
                           new_obs,        &! reserve memory
                           new_int,        &! reserve memory for interp.space
!                          set_int_insitu, &! set interpolation space
!                          set_vqc_insitu, &! subroutine to set VQC bounds
!                          shrink_report,  &! remove passive observations
!                          ITY_ICOL,       &! horizontally interpolated input
                           ITY_MCOLS,      &! interpolation type: column
                           CHR_LIN,        &! H is linear
                           TSK_INIT,       &! task: initialisation
                           TSK_READ,       &!       read observations
                           TSK_SET_CHR,    &!       set characteristics
                           TSK_SETUP_COLS, &!       set model columns
                           TSK_SET_CHR,    &!       set characteristics
                           TSK_SETUP_FUL0, &!       set interpolation space
                           TSK_SETUP_FULL, &!       set PSAS-space
                           TSK_Y,          &!       run forward operator
                           TSK_R,          &!       set observational error
                           TSK_SHRINK,     &!       release unused observations
                           SOIL             ! f90 module number
  use mo_obs_set,    only: t_obs_block      ! observation derived type
  use mo_obs_set,    only: t_obs_block      ! observation derived type
  use mo_obs_tables, only: rept_use,       &! report type usage table
                           decr_rpt_use,   &! degrade status of report
                           idb_dbk,        &! index in table rept_stat
                           check_report_1, &! basic checks on reports
                           check_report_0   ! basic checks on reports
  use mo_fdbk_tables,only: OT_SOIL,        &! soil observations
                           OC_ASCWS,       &! ascat soil moisture observations codetype
                           VN_NSOILM,      &! normalized soil moisture code
                           VN_SOILM,       &! volumetric soil moisture code
                           VN_DEPTH         ! depth below surface flag
  use mo_t_use,      only: t_use,          &! data type to hold state
                           use_0,          &! default values of type use
                           decr_use,       &! decrease the state of a report
!                          STAT_ACTIVE,    &! flag for active  observation
                           STAT_PASSIVE,   &! flag for passive observation
                           STAT_DISMISS,   &! flag for dismissed observation
                           STAT_REJECTED,  &! flag for rejected  observation
                           CHK_NOTUSED,    &! flag for unused or redundant data
                           CHK_INSDAT,     &! flag for insufficient data
                           CHK_DATASET,    &! flag for dataset specific flags
!                          CHK_BIASCOR,    &! flag for bias correction not appl.
                           CHK_OPERATOR,   &! obs. operator not applicable
                           CHK_SURF,       &! surface type flag
                           CHK_RULE,       &! rule flag value
                           CHK_FG,         &! first-guess check flag
                           CHK_DOMAIN       ! flag for out of domain
  !-----------------------------
  ! access atmospheric data type
  !-----------------------------
  use mo_atm_grid,   only: t_grid           ! grid derived data type
  use mo_atm_state,  only: t_atm            ! atm. state data type
  use mo_t_col,      only: t_cols,         &! model columns data type
                           COL_SOIL         ! specification of fields
  !------------------------------------
  ! Interface to interpolation routines
  !------------------------------------
  use mo_grid_intpol,only: idx_init         ! get grid-indices & interp.coeff.
  !------------------------------------
  ! access vector and matrix data types
  !------------------------------------
  use mo_dec_matrix, only: t_vector_segm    ! vector segment data type
  !--------------------------------------
  ! obs data header info read from NetCDF
  !--------------------------------------
  use mo_head_netcdf,only:                &!
                           istidn,        &! here: satellite identifier
!                          imissing,      &! NetCDF _FillValue for integer
                           rmissing,      &! NetCDF _FillValue for reals
                           s2ikz,         &! DWD-internal classifier
                           s1cent,        &! data centre
                           stime,         &! header observation time (section1)
                           db_time,       &! data bank time
                           s1cents,       &! data sub centre
                           mlah,          &! latitude
                           mloh,          &! longitude
                           obs_time,      &! body observation time
!                          ystidn          ! any type of station identifier as variable
                           get_int,       &! read integer from NetCDF file
                           get_real        ! read real    from NetCDF file
  implicit none

!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: process_soil               ! general purpose SOIL processing routine
  public :: read_soil_netcdf           ! read SOIL observation from netCDF file
  public :: read_nml_observations_soil ! read namelist /OBSERVATIONS_SOIL/
  public :: t_cdf_soil                 ! derived type to describe CDFs

!------------------------------------------------------------------------------
  !=========================
  ! Derived type definitions
  !=========================
  !-------------------------------------------------------------------
  ! Parameters for 13 percentile values: 0,5,10,20,30,...,80,90,95,100
  ! c.f. Y.Y. Liu et al., Hydrol. Earth Syst. Sci., 15, 425-436 (2011)
  !-------------------------------------------------------------------
  integer, parameter :: n_cdf_points = 13     ! number of CDF intervals
  type t_cdf_soil
     real(wp) :: perc     (1:n_cdf_points)    ! percentile
     real(wp) :: ascat    (1:n_cdf_points)    ! ascat sm obs (in %)
     real(wp) :: slope    (1:n_cdf_points)
     real(wp) :: intercept(1:n_cdf_points)
  end type t_cdf_soil

  !================================================
  ! Variables in Namelist Group /OBSERVATIONS_SOIL/
  !================================================
  logical          :: l_qc              = .false. ! flag to activate quality check controls
  logical          :: l_compl           = .false. ! check on topographic complexity
  logical          :: l_err             = .false. ! check on estimated soil moisture error
  logical          :: l_wet             = .false. ! check on inundation and wetland fraction
  logical          :: l_flag            = .true.  ! check on processing flag
  logical          :: l_stdv            = .false. ! check to discard obs if |o - fg| is greater than a value
  logical          :: flag_monitoring   = .false. ! flag to use MEC for wso monitoring (true) or for real use
  logical          :: flag_CDF_formula  = .false.
  logical          :: flag_NORM_formula = .true.
  logical          :: flag_soiltype1    = .false.
  logical          :: flag_soiltype9    = .true.
  logical          :: flag_soiltype16   = .false.
  logical          :: flag_lsm          = .true.
  logical          :: l_errinflation    = .true.  ! inflate wsoil obs error by a factor of 2
  real(wp)         :: stdv_inflation    =  4._wp  ! default value for stdv_inflation
  real(wp)         :: stdv_oi (n_st)    = -1._wp  ! stdv. for obs.increment check
  type(t_cdf_soil) :: cdf_soil(n_st)

  namelist /OBSERVATIONS_SOIL/ cdf_soil, l_qc, l_compl, l_err, l_wet, l_flag,     &
                               l_stdv, flag_CDF_formula, l_errinflation,          &
                               flag_NORM_formula, flag_soiltype9, flag_soiltype1, &
                               flag_soiltype16, flag_lsm, flag_monitoring,        &
                               stdv_inflation, stdv_oi

  !=================
  ! module variables
  !=================

!==============================================================================
contains
!==============================================================================

  subroutine read_nml_observations_soil
    !========================================
    ! read namelist group /OBSERVATIONS_SOIL/
    !========================================
    integer :: ierr, j
    logical :: first  = .true.   ! read namelist only once

    if (.not. first) return
    first             = .false.
    !-----------------------------------------------------------------------
    ! set defaults (should be consistent with default initialisations above)
    !-----------------------------------------------------------------------
    flag_CDF_formula  = .false.
    flag_NORM_formula = .true.
    l_qc              = .false.     ! flag to activate quality check controls
    l_compl           = .false.     ! check on topographic complexity
    l_err             = .false.     ! check on estimated soil moisture error
    l_wet             = .false.     ! check on inundation and wetland fraction
    l_flag            = .true.      ! check on processing flag
    l_stdv            = .false.     ! check to discard obs if |o - fg| is greater than a value
    flag_monitoring   = .false.     ! flag to use MEC for wso monitoring (true) or for real use
    flag_soiltype1    = .false.
    flag_soiltype9    = .true.
    flag_soiltype16   = .false.
    flag_lsm          = .true.
    l_errinflation    = .true.
    stdv_inflation    = 4._wp
    stdv_oi           = -1._wp
    do j=1,n_st
       cdf_soil(j)%perc     (1:n_cdf_points) = 0._wp
       cdf_soil(j)%ascat    (1:n_cdf_points) = 0._wp
       cdf_soil(j)%slope    (1:n_cdf_points) = 0._wp
       cdf_soil(j)%intercept(1:n_cdf_points) = 0._wp
    enddo

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('OBSERVATIONS_SOIL', status=ierr)
    end if
    call p_bcast (ierr, dace% pio)
    if (ierr /= POSITIONED) then
      if (dace% lpio) then
        write (6,'(a)')   repeat('-',79)
        write (6,'(a)') ' Namelist /OBSERVATIONS_SOIL/ not present'
      end if
      return
    end if
    if (dace% lpio) then
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=OBSERVATIONS_SOIL, iostat=ierr)
        if (ierr/=0) call finish ('read_nml_observations_soil',          &
                                  'ERROR in namelist /OBSERVATIONS_SOIL/')
#else
        read (nnml ,nml=OBSERVATIONS_SOIL)
#endif
      end select

      !---------
      ! printout
      !---------
      write(6,'(a)') repeat('_',79)
      write(6,'()')
      write(6,'(a)')      '  Soil Observations namelist /OBSERVATIONS_SOIL/ :'
      write(6,'()')
      write(6,'(a,l1)')'    l_qc      = ', l_qc
      write(6,'(a,l1)')'    l_compl   = ', l_compl
      write(6,'(a,l1)')'    l_err     = ', l_err
      write(6,'(a,l1)')'    l_wet     = ', l_wet
      write(6,'(a,l1)')'    l_flag    = ', l_flag
      write(6,'(a,l1)')'    l_stdv    = ', l_stdv
      write(6,'(a,l1)')'    l_errinflation       = ', l_errinflation
      write(6,'(a,l1)')'    flag_CDF_formula     = ', flag_CDF_formula
      write(6,'(a,l1)')'    flag_NORM_formula    = ', flag_NORM_formula
      write(6,'(a,l1)')'    flag_monitoring      = ', flag_monitoring
      write(6,'(a,l1)')'    flag_soiltype1       = ', flag_soiltype1
      write(6,'(a,l1)')'    flag_soiltype9       = ', flag_soiltype9
      write(6,'(a,l1)')'    flag_soiltype16      = ', flag_soiltype16
      write(6,'(a,l1)')'    flag_lsm             = ', flag_lsm
      write(6,'(a,f10.5)')'    stdv_inflation       = ', stdv_inflation
      do j=1,n_st
       select case (j)
       case default
        write(6,'(a,i2,a,f10.5)')'    stdv_oi(',j,')          = ', stdv_oi(j)
       case (ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE)
        ! skip
       end select
      enddo
      do j=1,n_st
       select case (j)
       case default
        write(6,'(a,i2,a,13f10.5)')'  cdf_soil(',j,')%perc      = ', cdf_soil(j)%perc
        write(6,'(a,i2,a,13f10.5)')'  cdf_soil(',j,')%ascat     = ', cdf_soil(j)%ascat
        write(6,'(a,i2,a,13f10.5)')'  cdf_soil(',j,')%slope     = ', cdf_soil(j)%slope
        write(6,'(a,i2,a,13f10.5)')'  cdf_soil(',j,')%intercept = ', cdf_soil(j)%intercept
       case (ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE)
        ! skip
       end select
      enddo
    end if

    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    do j=1,n_st
       call p_bcast (cdf_soil(j)% perc     ,dace% pio)
       call p_bcast (cdf_soil(j)% ascat    ,dace% pio)
       call p_bcast (cdf_soil(j)% slope    ,dace% pio)
       call p_bcast (cdf_soil(j)% intercept,dace% pio)
    enddo
    call p_bcast  (l_qc                    ,dace% pio)
    call p_bcast  (l_compl                 ,dace% pio)
    call p_bcast  (l_err                   ,dace% pio)
    call p_bcast  (l_wet                   ,dace% pio)
    call p_bcast  (l_flag                  ,dace% pio)
    call p_bcast  (l_stdv                  ,dace% pio)
    call p_bcast  (flag_CDF_formula        ,dace% pio)
    call p_bcast  (flag_NORM_formula       ,dace% pio)
    call p_bcast  (flag_monitoring         ,dace% pio)
    call p_bcast  (flag_soiltype1          ,dace% pio)
    call p_bcast  (flag_soiltype9          ,dace% pio)
    call p_bcast  (flag_soiltype16         ,dace% pio)
    call p_bcast  (flag_lsm                ,dace% pio)
    call p_bcast  (l_errinflation          ,dace% pio)
    call p_bcast  (stdv_inflation          ,dace% pio)
    call p_bcast  (stdv_oi                 ,dace% pio)

    if (count ([flag_soiltype1, flag_soiltype9, flag_soiltype16]) /= 1) then
       call finish ("read_nml_observations_soil", &
                    "need one out of flag_soiltype1, flag_soiltype9, flag_soiltype16")
    end if
    !-------------------------
    ! Checks on cdf_soil input
    !-------------------------
    if (flag_CDF_formula) then
       do j=1,n_st
          select case (j)
          case default
             ! check for valid settings (TODO)
             if (all (cdf_soil(j)%ascat(1:n_cdf_points  ) == 0._wp) .or. &
                 any (cdf_soil(j)%ascat(1:n_cdf_points-1) >=             &
                      cdf_soil(j)%ascat(2:n_cdf_points  )         )      ) then
                call finish ("read_nml_observations_soil", "bad cdf_soil%ascat")
             end if
          case (ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE)
             ! skip
          end select
       enddo
    end if

  end subroutine read_nml_observations_soil

!------------------------------------------------------------------------------
  subroutine process_soil (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
                            state)
  integer            ,intent(in)             :: task    ! what to do
  type(t_spot)       ,intent(inout),optional :: spot    ! SPOT observations
  type(t_obs_block)  ,intent(inout),optional :: obs     ! observation data type
  type(t_atm)        ,intent(in)             :: atm     ! atmospheric state
  type(t_cols)       ,intent(in)   ,optional :: cols    ! model columns
  type(t_vector_segm),intent(in)   ,optional :: xi      ! interpolated values
  type(t_vector_segm),intent(inout),optional :: y       ! observed quantity
  real(wp)           ,intent(inout),optional :: Jo      ! obs. cost funct. Jo
  type(t_atm)        ,intent(inout),optional :: Jo_atm  ! gradient:d Jo/d atm
  integer            ,intent(in)   ,optional :: state   ! status flag

    !================
    ! local variables
    !================
    integer  :: tsk      ! task (local copy)
    integer  :: i, n, j  ! observation index and len
    integer  :: m, mm    ! loop index
    integer  :: k        ! index for sparse matrix representation
    integer  :: index_2  ! index useful in calculations
    integer  :: ic, ic2  ! indices for temporal interpolation
    integer  :: soiltyp  ! soil type
    integer  :: st(n_st) ! soil type count
    integer  :: order    ! interpolation order (determines no. gridpoints)
    integer  :: np       ! number of points
    integer  :: ij(2)    ! horizontal index subset
    real(wp) :: lsm      ! land fraction (land sea mask)
    real(wp) :: w1, w2   ! temporal interpolation weights
    real(wp) :: w_h      ! horizontal integration weight
    real(wp) :: sumw     ! weight sum
    real(wp) :: w_so(2)  ! soil moisture
    real(wp) :: wso_eff  ! effective soil moisture for layer 1 (1cm)
    real(wp) :: ssm_raw  ! original surface soil moisture data
    real(wp) :: ssm_tr   ! CDF transformed surface soil moisture data
    real(wp) :: tmp      ! temporary
    type(t_grid), pointer :: grid
    type(t_mcol), pointer :: mc(:)
    integer, parameter    :: mp = size (spot% col% h% imc,dim=1)
    real(wp)              :: d(mp)      ! Auxiliary array for distance estimate
    integer               :: soiltyp_m(mp)
    real(wp), allocatable :: ws_obs(:)     ! obs values in m3/m3 in the 9 or 16 points closest to the observations
    !======================
    ! executable statements
    !======================
    tsk = task
    if(tsk==0) return

    !==================================================
    ! skip tasks not required for this observation type
    !==================================================
    tsk = task
    tsk = iand (task, not (&
          TSK_READ         &! BUFR file is read in module mo_obs
        + TSK_SHRINK       &! release unused obs. in report
    ))
    if (tsk == 0) return

    !===============
    ! initialisation
    !===============
    if (iand (TSK_INIT,tsk) /= 0) then
      call read_nml_observations_soil   ! read namelist /OBSERVATIONS_SOIL/
      tsk = tsk - TSK_INIT
      if (tsk == 0) return
    endif

    grid => atm% grid
    if (.not. associated (grid% soiltyp))                    &
         call finish ("process_soil","soiltyp not associated")

    !=========================================
    ! TSK_SETUP_COLS: specify required columns
    !=========================================
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then
      select case (grid% gridtype)
      case (WMO6_LATLON, WMO6_ROTLL, WMO6_GAUSSIAN)
         if (flag_soiltype1) then
            order = 1   ! nearest neighbor
         !else if () then
         !   order = 2   ! 4-point stencil
         else if (flag_soiltype9) then
            order = 4    ! use 16-point stencil as template and adjust below
         else
            order = 4
         end if
      case (DWD6_ICON)
         if (flag_soiltype1) then
            order = 1   ! nearest neighbor
         else if (flag_soiltype9) then
            order = 2   ! 3 closest gridpoints
         else
            order = 4   ! 6 surrounding gridpoints
         end if
      case default
         write(0,*) "process_soil: unsupported gridtype =", grid% gridtype
         call finish ("process_soil","unsupported gridtype")
      end select
      !---------------------------------------------------------------
      ! We (ab)use the setup for horizontal interpolation to determine
      ! the gridpoints contributing to the stencil and adjust the
      ! interpolation weights later (e.g. cell area weighting).
      !---------------------------------------------------------------
      call idx_init (       &
             spot% col% c,  &! <-  column descriptor
             spot% col% h,  &!  -> interpolation coefficients
             obs% o% mc,    &! <-> model column descriptors
             COL_SOIL,      &! <-  fields required
             0,             &! <-  tracers required
             grid,          &! <-  model grid
             spot% i_time,  &! <-  time slot
             spot% w_time,  &! <-  time interpolation weight
       order=order          )! <~  nearest neighbour flag

      try: do
         if (spot% col% h% imc(1,1) == 0) then
            call decr_rpt_use (spot, CHK_DOMAIN, &
                               use=min(STAT_DISMISS, rept_use(OT_SOIL)% use(CHK_DOMAIN)))
            exit try
         end if

         mc => obs% o% mc% c

         select case (grid% gridtype)
         case (WMO6_LATLON, WMO6_ROTLL, WMO6_GAUSSIAN)
            !---------------------------------------------------------
            ! Remove from stencil gridpoints in same line (row/column)
            ! as most distant point.  Weights will be adjusted later.
            !---------------------------------------------------------
            if (flag_soiltype9) then
               np = 0
               do i = 1, mp
                  ic = spot% col% h% imc(i,1)
                  if (ic == 0) exit
                  np = np + 1
                  ij = mc(ic)% ijdtp(1:2)
                  d(i) = sum ((grid% xnglob(ij(1),ij(2),1:3,1) - spot% col% c% x)**2)
               end do
               i  = maxloc (d(1:np), dim=1)
               ic = spot% col% h% imc(i,1)
               ij = mc(ic)% ijdtp(1:2)
!              if (spot% col% h% ijdp(1) /= ij(1)) &
!                   write(0,*) "### spot% col% h% ijdp=",spot% col% h% ijdp(1)
               np = 0
               do i = 1, mp
                  ic = spot% col% h% imc(i,1)
                  if (ic == 0) exit
                  if (mc(ic)% ijdtp(1) == ij(1)) cycle
                  if (mc(ic)% ijdtp(2) == ij(2)) cycle
                  np = np + 1
                  if (np < i) then
                     spot% col% h% imc(np,:) = spot% col% h% imc(i,:)
                  end if
               end do
               spot% col% h% imc(np+1,1) = 0
               if (np /= 9) then
                  write(0,*) "### process_soil: np/=9, np =", np
                  write(0,*) "### lon,lat=",spot%col%c% dlon, spot%col%c% dlat
                  call decr_rpt_use (spot, CHK_DOMAIN, STAT_DISMISS)
                  exit try
               end if
            end if
         end select

         !-------------------
         ! check the np value
         !-------------------
         np = 0
         do i = 1, size (spot% col% h% imc,1)
            ic = spot% col% h% imc(i,1)
            if (ic == 0) exit
            np = np + 1
         enddo

         !--------------------------------------------------------------
         ! temporarily set non-zero equal weights (TODO: area weighting)
         !--------------------------------------------------------------
         if (np > 1) then
            spot% col% h% w(1:np) = 1._wp / np
         end if
         exit try
      end do try
      tsk = tsk - TSK_SETUP_COLS
      if (tsk == 0) return
    endif

    !=============================================
    ! TSK_SET_CHR: set observation characteristics
    !=============================================
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      spot% int_type  = ITY_MCOLS  ! expect single model column
      spot% cost      = 1._wp      ! operator is cheap
      spot% nr        = spot% o% n ! diagonal R
      spot% char = CHR_LIN         ! this is a linear operator
      tsk = tsk - TSK_SET_CHR
      if (tsk == 0) return
    endif

    !===========================================================
    ! tsk == TSK_SETUP_FUL0:
    ! interpolation space
    ! (currently not used. required for 3dvar, not verification)
    !===========================================================
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
      call new_int (obs% o, spot, 0)
      tsk = tsk - TSK_SETUP_FUL0
      if (tsk == 0) return
    endif

    !===========================================
    ! tsk == TSK_SETUP_FULL:
    ! setup description of PSAS-space
    ! observed values were set up while reading.
    !===========================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      tsk = tsk - TSK_SETUP_FULL
      if (tsk == 0) return
    endif

    !=========================================
    ! run operator or set up H (Jacobi-matrix)
    !=========================================
    if (iand (TSK_Y,tsk) /= 0) then
      if (spot% pe_eval == dace% pe) then
        !-------------------------------
        ! get indices, check consistence
        !-------------------------------
        i = spot%o% i + 1
        n = spot%o% n
        if (n >  2) call finish ('process_soil','o%n > 2')
        if (obs% o% varno(i) /= VN_NSOILM)        &
          call finish ('process_soil','not NSOILM')
        if (n == 2) then
          if (obs% o% varno(i+1) /= VN_SOILM)    &
          call finish ('process_soil','not SOILM')
        end if
        obs% o% body(i:i+n-1)% plev = spot% ps_bg
        !--------------------------------
        ! checks for land, not stone etc.
        ! adjust weights for averaging.
        !--------------------------------
        st = 0
        np = 0
        sumw = 0._wp
        mc => obs% o% mc% c
        soiltyp_m = 0
        do j = 1, size (spot% col% h% imc,1)
           ic = spot% col% h% imc(j,1)
           if (ic == 0) exit
           np = np + 1
           soiltyp = grid% soiltyp (mc(ic)% ijdtp(1),    &
                                    mc(ic)% ijdtp(2), 1, &
                                    mc(ic)% ijdtp(3)     )
           lsm     = grid% lsm     (mc(ic)% ijdtp(1),    &
                                    mc(ic)% ijdtp(2), 1, &
                                    mc(ic)% ijdtp(3)     )
           soiltyp_m(np) = soiltyp
           select case (soiltyp)
           case default
              st(soiltyp) = st(soiltyp) + 1
           case (ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE)
              call decr_rpt_use (spot, CHK_SURF, comment="soiltype")
              if (spot% use% state <= STAT_DISMISS) exit
           !case (:0, n_st+1:)          ! must never happen...
           !   call decr_rpt_use (spot, CHK_SURF, STAT_DISMISS)
           !   exit try
           end select
           if (flag_lsm .and. lsm < 0.5_wp) then
              call decr_rpt_use (spot, CHK_SURF, comment="lsm")
              if (spot% use% state <= STAT_DISMISS) exit
           end if
        end do
print *, 'spot e soiltyp_m: ', i,soiltyp_m(:)
        !----------------------------------------------
        ! determine representative or dominant soiltype
        !----------------------------------------------
        !if (all (st == 0)) then
        !   call decr_rpt_use (spot, CHK_SURF, STAT_DISMISS, comment="soiltype")
        !   exit try
        !end if
        soiltyp = maxloc (st, dim=1)
        spot% soiltype = soiltyp
        !------------------------------------------------------------
        ! derive equivalent volumetric soil moisture from observation
        !------------------------------------------------------------
        if (n == 2) then
           if (flag_CDF_formula) then
    if (flag_soiltype1) then
             soiltyp = spot% soiltype
             if (soiltyp<=ST_ROCK .or. soiltyp>=ST_SEAWATER) then
              obs% o% body(i+1)% o  = rvind
              obs% o% body(i+1)% eo = rvind
             else
              ssm_raw = obs% o% body(i)% o        ! observed SMI (0..1)
              tmp     = ssm_raw * 100._wp         ! convert to %
              if (tmp >= cdf_soil(soiltyp)% ascat(      1     ) .and. &
                  tmp <= cdf_soil(soiltyp)% ascat(n_cdf_points)       ) then
                 do m = 2, n_cdf_points
                    if (tmp <= cdf_soil(soiltyp)% ascat(m)) exit
                 end do
                 m = m - 1
                 tmp    =           cdf_soil(soiltyp)%slope(m)    ! Jacobian
                 ssm_tr = ssm_raw * cdf_soil(soiltyp)%slope(m) + &! vol. soil
                                    cdf_soil(soiltyp)%intercept(m)! moisture
                 obs% o% body(i+1)% o  = ssm_tr
                 obs% o% body(i+1)% eo = obs% o% body(i)% eo * ssm_tr
              else
                 !------------------------------------
                 ! Out of domain of CDF transformation
                 !------------------------------------
                 call decr_use (obs% o% body(i  )% use, &
                                check = CHK_FG,         &
                                lflag = .true.          )
                 call decr_use (obs% o% body(i+1)% use, &
                                state = STAT_DISMISS,   &
                                check = CHK_OPERATOR,   &
                                lflag = .true.          )
                 obs% o% body(i+1)% o = rvind
              end if
     end if
    else ! not flag_soiltype1
     allocate (ws_obs(np))
     ws_obs(:) = 0._wp
     index_2 = 0
     ssm_raw = obs% o% body(i)% o        ! observed SMI (0..1)
     tmp     = ssm_raw * 100._wp
     do m=1,np
      if (soiltyp_m(m)>ST_ROCK .and. soiltyp_m(m)<ST_SEAWATER) then
       if (tmp >= cdf_soil(soiltyp_m(m))% ascat(      1     ) .and. &
           tmp <= cdf_soil(soiltyp_m(m))% ascat(n_cdf_points)       ) then
                 do mm = 2, n_cdf_points
                    if (tmp <= cdf_soil(soiltyp_m(m))% ascat(mm)) exit
                 end do
                 mm = mm - 1
                 tmp    =           cdf_soil(soiltyp_m(m))%slope(mm)    ! Jacobian
                 ssm_tr = ssm_raw * cdf_soil(soiltyp_m(m))%slope(mm) + &! vol. soil
                                    cdf_soil(soiltyp_m(m))%intercept(mm)! moisture
                 ws_obs(m)  = ssm_tr
                 index_2 = index_2 + 1
               else
                 !------------------------------------
                 ! Out of domain of CDF transformation
                 !------------------------------------
                 call decr_use (obs% o% body(i  )% use, &
                                check = CHK_FG,         &
                                lflag = .true.          )
                 call decr_use (obs% o% body(i+1)% use, &
                                state = STAT_DISMISS,   &
                                check = CHK_OPERATOR,   &
                                lflag = .true.          )
                 ws_obs(m) = rvind
               end if
      end if
     end do
     if (index_2 /= 0) then
              obs% o% body(i+1)% o = sum(ws_obs)/index_2
              obs% o% body(i+1)% eo = obs% o% body(i)% eo * obs% o% body(i+1)% o
     else
              obs% o% body(i+1)% o  = rvind
              obs% o% body(i+1)% eo = rvind
     endif
     deallocate(ws_obs)
    end if ! not flag_soiltype1

           else ! if (flag_NORM_formula) then
    if (flag_soiltype1) then
              soiltyp = spot% soiltype
      if (soiltyp<=ST_ROCK .or. soiltyp>=ST_SEAWATER) then
               obs% o% body(i+1)% o  = rvind
               obs% o% body(i+1)% eo = rvind
      else
               ssm_raw = obs% o% body(i)% o                ! observed SMI (0..1)
               tmp     = Cporv(soiltyp) - Cadp(soiltyp)    ! Jacobian
               ssm_tr  = ssm_raw * tmp  + Cadp(soiltyp)
               obs% o% body(i+1)% o  = ssm_tr
               obs% o% body(i+1)% eo = obs% o% body(i)% eo * ssm_tr
      endif
    else ! not flag_soiltype1
      allocate (ws_obs(np))
      ws_obs(:) = 0._wp
      index_2 = 0
      ssm_raw = obs% o% body(i)% o        ! observed SMI (0..1)
      tmp     = ssm_raw * 100._wp
      print *, 'nppp: ', np, soiltyp_m(:)
      do m=1,np
       if (soiltyp_m(m)>ST_ROCK .and. soiltyp_m(m)<ST_SEAWATER) then
                 do mm = 2, n_cdf_points
                    if (tmp <= cdf_soil(soiltyp_m(m))% ascat(mm)) exit
                 end do
                 mm = mm - 1
                 tmp    =         Cporv(soiltyp_m(m)) - Cadp(soiltyp_m(m))    ! Jacobian
                 ssm_tr = ssm_raw * tmp  + Cadp(soiltyp_m(m))
                 ws_obs(m)  = ssm_tr
                 index_2 = index_2 + 1
               end if
      end do
      if (index_2 /= 0) then
               obs% o% body(i+1)% o = sum(ws_obs)/index_2
               obs% o% body(i+1)% eo = obs% o% body(i)% eo * obs% o% body(i+1)% o
      else
               obs% o% body(i+1)% o  = rvind
               obs% o% body(i+1)% eo = rvind
      endif
      deallocate(ws_obs)
    end if  ! flag_soiltyp1
           end if  !flag_CDF_formula
        end if
        !------------------
        ! evaluate operator
        !------------------
        sumw = 0._wp
        w_so = 0._wp
        np   = 0
        do m = 1, size (spot% col% h% imc,1)
           ic = spot% col% h% imc(m,1)
           if (ic == 0) exit
           !---------------------------
           ! all or only same soiltype?
           !---------------------------
           !soiltyp = grid% soiltyp (cols% col(ic)% i,&
           !                         cols% col(ic)% j,&
           !                      1, cols% col(ic)% l )
           !if (soiltyp /= spot% soiltype) cycle
           np   = np + 1
           w_h  = 1._wp
           sumw = sumw + w_h
           if (spot% w_time == 0._wp) then
             !--------------------------
             ! no temporal interpolation
             !--------------------------
             w_so = w_so + w_h * cols% col(ic)% w_so(1:2)
           else
             !-----------------------
             ! temporal interpolation
             !-----------------------
             ic2  = spot% col% h% imc(m,2)
             w2   = spot% w_time
             w1   = 1._wp - w2
             w_so = w_so + w_h * ( w1 * cols% col(ic )% w_so(1:2) &
                                 + w2 * cols% col(ic2)% w_so(1:2) )
           endif
        end do
        if (np > 0) w_so = w_so / sumw
        !---------------------------------------------------------------
        ! derive equivalent soil moisture (kg/m^2) for layer 1 (top 1cm)
        ! by effectively averaging over layers 1 (0-1cm) and 2 (1-3cm):
        !---------------------------------------------------------------
        wso_eff = (w_so(1) + 0.5_wp * w_so(2)) * 0.5_wp   ! avg. over top 2cm
        !-------------------------------------------
        ! soil moisture index
        ! mode:  1: use adp,porv; 2: use pwp,fcap
        !        ind_wso (wso, soiltyp, level, mode)
        ! derive normalized soil moisture
        !-------------------------------------------
        if (spot% soiltype > ST_ROCK .and. spot% soiltype < ST_SEAWATER) then
           y%x (i) = ind_wso (wso_eff, spot% soiltype, 1, 1)
        else
           y%x (i) = rvind
        end if
        if (n == 2) then
           y%x(i+1) = wso_eff / (1000*0.01_wp) ! convert kg/m^2/0.01m -> m^3/m^3
        end if

        if (obs% o% body(i+1)% o == rvind) y% x (i+1) = rvind

        !-----------------------------------------------
        ! Gross check based on obs increments |obs - fg|
        !-----------------------------------------------
        if (l_stdv) then
           tmp = stdv_inflation * stdv_oi(spot% soiltype)
           if (tmp > 0._wp) then
              if (n == 2) then
                 if (abs (obs%o%body(i+1)% o - y%x(i+1)) >= tmp) then
                    call decr_use (obs%o%body(i  )% use, STAT_PASSIVE, check=CHK_FG)
                    call decr_use (obs%o%body(i+1)% use, STAT_PASSIVE, check=CHK_FG)
                 end if
              else
                 if (abs (obs%o%body(i  )% o - y%x(i  )) >= tmp) then
                    call decr_use (obs%o%body(i  )% use, STAT_PASSIVE, check=CHK_FG)
                 end if
              end if
           end if
        end if
      endif
      tsk = tsk - TSK_Y
      if (tsk == 0) return
    endif

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    !===============================================
    if (iand (TSK_R,tsk) /= 0) then
     if (dace% pe == obs% o% pe) then
      n = spot% o% n
      k = obs%  R% ia (spot% o% i+1)
      i = spot% o% i+1
      obs% R% ia (i) = k
      select case (obs% o% varno (spot%o% i+1))
      case (VN_NSOILM)
        obs% R% packed(k) = 0.5_wp
      case default
        call finish('process_soil (TSK_R)','invalid observation type')
      end select
      obs% R% ja (k) = i
      k = k + 1
      obs% R% ia (spot% o% i + n + 1) = k
     endif
     tsk = tsk - TSK_R
     if (tsk == 0) return
    endif


!   .....

    !===========
    ! left tasks
    !===========
    if (tsk /= 0) then
      if (dace% lpio) write (6,*) 'process_soil: unknown task',tsk
      call finish ('process_soil','unknown task')
    endif

  end subroutine process_soil

!==============================================================================

  subroutine read_soil_netcdf (ifile, i_source, obs, head, lkeep, nkeep)
  integer       ,intent(in)           :: ifile     ! Number of netCDF file read
  integer       ,intent(inout)        :: i_source  ! number of records in source-file
  type (t_obs)  ,intent(inout)        :: obs       ! observations data type to set
  type (t_head) ,intent(in)           :: head      ! header data already encoded
  logical       ,intent(out)          :: lkeep     ! accept observation ?
  integer       ,intent(out)          :: nkeep     ! number of accepted obsvs.

    !================
    ! local variables
    !================
    type (t_use)          :: use        ! status variable
    type (t_head)         :: hd         ! report header
    type (t_spot)         :: spt0, spti ! report meta data
    type (t_spot) ,save   :: empty      !
    integer               :: len_report ! number of reports from file
    integer               :: n          ! number of reports to process
    integer               :: is         ! report index variable
    integer               :: entry      ! total report count
    logical               :: lk         ! flag to keep this report
    real(sp) ,allocatable :: nssm  (:)  ! surface soil moisture
    integer               :: subtyp     ! observation subtype (Datenbankkennz.)
    integer               :: subtyp_idx ! observation subtype index
    !===========================================
    !list of local variables for quality control
    !===========================================
    integer  ,allocatable :: nsmpf (:)  ! sm processing flag
    integer  ,allocatable :: nsmcf (:)  ! sm correction flag
    real(sp) ,allocatable :: nsmq  (:)  ! sm quality (%)
    real(sp) ,allocatable :: neesm (:)  ! estimated error in surface soil moisture (%)
    integer  ,allocatable :: nsnco (:)  ! snow cover (%)
    real(sp) ,allocatable :: nflsf (:)  ! frozen land surface fraction (%)
    real(sp) ,allocatable :: niawf (:)  ! inundation and wetland fraction (%)
    real(sp) ,allocatable :: ntoco (:)  ! topographic complexity (%)
    integer  ,allocatable :: nctcn (:)  ! cross-track cell number
    integer  ,allocatable :: msain (:)  ! Satellite instrument

    lkeep = .false.
    nkeep = 0

    !------------------------
    ! get dimension of fields
    !------------------------
    if (.not. allocated (istidn)) call finish ('read_soil_netcdf','istidn')
    len_report = size (istidn)
    n          = min (len_report, rept_use(OT_SOIL)% max_proc)

    !----------------
    ! allocate arrays
    !----------------
    allocate (nssm (len_report))
    allocate (nsmq (len_report))
    allocate (neesm(len_report))
    allocate (nsnco(len_report))
    allocate (nflsf(len_report))
    allocate (niawf(len_report))
    allocate (ntoco(len_report))
    allocate (nsmcf(len_report))
    allocate (nsmpf(len_report))
    allocate (nctcn(len_report))
    allocate (msain(len_report))

    !--------------------------------
    ! read variables from NetCDF file
    !--------------------------------
    call get_real (nssm, 'NSSM' , -999._sp)
    call get_real (nsmq, 'NSMQ' , -999._sp)
    call get_real (neesm,'NEESM', -999._sp)
    call get_int  (nsnco,'NSNCO', -1      )
    call get_real (nflsf,'NFLSF', -999._sp)
    call get_real (niawf,'NIAWF', -999._sp)
    call get_real (ntoco,'NTOCO', -999._sp)
    call get_int  (nsmcf,'NSMCF', -1      )
    call get_int  (nsmpf,'NSMPF', -1      )
    call get_int  (nctcn,'NCTCN', -1      )
    call get_int  (msain,'MSAIN', -1      )

    if (all (nssm == -999._sp)) return

    !-------------------------------
    ! preset total number of reports
    !-------------------------------
    entry = sum (source(1:ifile-1)% entries)

    !------------------
    ! loop over reports
    !------------------
    do is = 1, n

      entry   = entry    + 1

      !-----------------------
      ! initialize usage flags
      !-----------------------
      use = use_0

      !--------------------
      ! define head section
      !--------------------
      hd            = head
      hd% modtype   = SOIL
      hd% obstype   = OT_SOIL
      hd% codetype  = OC_ASCWS
      hd% time      = stime   (is)
      hd% db_time   = db_time (is)
      hd% source    = ifile
      hd% record    = is
      hd% id        = entry
      hd% center    = s1cent  (is)
      hd% subcenter = s1cents (is)
      subtyp        = s2ikz   (is)
      if (subtyp > 0) then
         subtyp_idx = idb_dbk (subtyp, OT_SOIL)
         hd% idbk   = subtyp_idx
      end if

!     if (is < 100) then
      if (.FALSE. ) then
        write (6,'()')
        write (6,'( 8(a16, i8,/),   &
                  & 2(a16, a ,/),   &
                  & 6(a16, i8,/) )' )                      &
          'pe='         ,dace% pe,                         &
          'head is='    ,is,                               &
          'obstype='    , hd% obstype  ,                   &
          'dbkz='       , hd% dbkz     ,                   &
          'modtype='    , hd% modtype  ,                   &
          'buf_type='   , hd% buf_type ,                   &
          'buf_subtype=', hd% buf_subtype,                 &
          'codetype='   , hd% codetype ,                   &
          'time='       , cyyyymmddhhmmss (hd% time)   ,   &
          'db_time='    , cyyyymmddhhmmss (hd% db_time),   &
          'dbk='        , hd% idbk     ,                   &
          'source='     , hd% source   ,                   &
          'record='     , hd% record   ,                   &
          'id='         , hd% id       ,                   &
          'center='     , hd% center   ,                   &
          'subcenter='  , hd% subcenter
      endif

      !--------------------------------------------
      ! perform simple generic check on report type
      !--------------------------------------------
      call check_report_0 (use, hd, 1)
      if (use% state <= STAT_DISMISS) cycle

      !------------------
      ! create new report
      !------------------
      spt0             = empty
      spt0% use        = use
      spt0% hd         = hd
      spti             = spt0

      !--------------------------
      ! check for sufficient data
      !--------------------------
      if (nssm  (is) == -999._sp .or. &
          mlah  (is) == rmissing .or. &
          mloh  (is) == rmissing      ) then
        call decr_rpt_use (spti, CHK_INSDAT, comment='read_soil_netcdf')
        cycle
      endif

      !------------------------------
      ! process report header entries
      !------------------------------
      spti% col% c% dlat = mlah     (is)
      spti% col% c% dlon = mloh     (is)
      spti% sttyp        = msain    (is)
!     spti% z            =
      spti% actual_time  = obs_time (is)
      spti% ident        = istidn   (is)
      spti% col% nlev    = 1
      call set_xuv (spti)
      !------------------------
      ! set center / processing
      !------------------------
      spti% statid = 'ASCAT-SOIL'

!     if (is < 100) then
      if (.FALSE. ) then
        write (6,'()')
        write (6,'(   a20, i6  ,a, /, &
                  & 2(a20,f8.3 ,   /),&
                  &   a20, a   ,   / ,&
                  &   a20, a   ,   / ,&
                  & 6(a20, i5      /))' )                    &
              'pe=',dace% pe,'  spti ',                      &
              'spti% col% c% dlat = ', spti% col% c% dlat ,  &
              'spti% col% c% dlon = ', spti% col% c% dlon ,  &
              'spti% actual_time  = ', cyyyymmddhhmmss (spti% actual_time),  &
              'spti% statid       = ', spti% statid       ,  &
              'spti% ident        = ', spti% ident        ,  &
              'spti% col% nlev    = ', spti% col% nlev
      endif

      !----------------
      ! standard checks
      !----------------
      call check_report_1 (spti)
      lk = spti% use% state > STAT_DISMISS
      if (lk) then
        call check_store_soil (spti, obs, lk, nssm(is), nsmpf(is), nsmcf(is), &
                               nsmq(is), neesm(is), nsnco(is), nflsf(is),     &
                               niawf(is), ntoco(is), nctcn(is)                )
        if (lk) nkeep = nkeep + 1
      endif
    end do

    write(6,*) '  read_soil_netcdf: keep',nkeep,' reports out of',n

  end subroutine read_soil_netcdf

!==============================================================================

  subroutine check_store_soil (spot, obs, lkeep, nssm,nsmpf, nsmcf, nsmq, &
                               neesm, nsnco, nflsf, niawf, ntoco, nctcn   )
  type(t_spot)    ,intent(inout)  :: spot    ! report derived type
  type(t_obs)     ,intent(inout)  :: obs     ! observations derived type
  logical         ,intent(out)    :: lkeep   ! keep or reject observation
  real(sp)        ,intent(in)     :: nssm    ! surface soil moisture
  integer         ,intent(in)     :: nsmpf   ! sm processing flag
  integer         ,intent(in)     :: nsmcf   ! sm correction flag
  real(sp)        ,intent(in)     :: nsmq    ! sm quality
  real(sp)        ,intent(in)     :: neesm   ! estimated error in surface sm
  integer         ,intent(in)     :: nsnco   ! snow cover
  real(sp)        ,intent(in)     :: nflsf   ! frozen land surface fraction
  real(sp)        ,intent(in)     :: niawf   ! inundation and wetland fraction
  real(sp)        ,intent(in)     :: ntoco   ! topographic complexity
  integer         ,intent(in)     :: nctcn   ! cross track cell number

    type(t_spot),pointer :: spt          ! temporary
    type(t_datum)        :: bod          ! body derived type
    integer              :: id           ! observation id
    integer              :: i1, in       ! index range
    integer              :: no           ! number of observations in the report

    lkeep = .true.

!   !----------------------
!   ! evaluate quality flag
!   !----------------------
!   call decr_rpt_use (spot, CHK_NOTUSED, STAT_PASSIVE, &
!                                         comment='....')
!   !------------------------------
!   ! exit if no valid data present
!   !------------------------------
!   if (.........) then
!     lkeep = .false.
!     call decr_rpt_use (spot, CHK_INSDAT, comment='....')
!     return
!   endif
!   !--------------------------
!   ! report selection (filter)
!   !--------------------------
!  if (...) call decr_rpt_use (spot, CHK_BLACKLIST)

    !==========================
    !ASCAT-SOIL QUALITY CONTROL
    !==========================

    if (l_flag) then
       ! BUFR Code Table 0 40 006 - Soil moisture processing flag
       ! Bit No.
       !    1   | Not soil
       !    2   | Sensitivity to soil moisture below limit
       !    3   | Azimuthal noise above limit
       !    4   | Backscatter Fore-Aft beam out of range
       !    5   | Slope Mid-Fore beam out of range
       !    6   | Slope Mid-Aft beam out of range
       !    7   | Soil moisture below -20%
       !    8   | Soil moisture above 120%
       !  9-15  | Reserved
       ! All 16 | Missing value
       if (nsmpf /= 0) call decr_rpt_use (spot, CHK_DATASET,            &
                                          comment='l_flag quality check')
    endif

    if (l_qc) then
       if (l_compl) then
          if (ntoco>20.0) call decr_rpt_use (spot, CHK_RULE, STAT_REJECTED, &
                                             comment='l_compl quality check')
       endif
       if (l_err) then
          if (neesm> 8.0) call decr_rpt_use (spot, CHK_RULE, STAT_REJECTED,&
                                             comment='l_err quality check' )
       endif
       if (l_wet) then
          if (niawf>15.0) call decr_rpt_use (spot, CHK_RULE, STAT_REJECTED,&
                                             comment='l_wet quality check' )
       endif
    endif   !l_qc

    if (spot% use% state > STAT_DISMISS) then
      !--------------------------------------------------
      ! new report header entry in DACE  observation list
      !--------------------------------------------------
      call new_spot (obs,1, set_id=.true.)
      spt => obs% spot (obs% n_spot)
      id      = spt% id
      spt     = spot
      spt% id = id
      if (flag_monitoring) then
         no   = 1
      else
         no   = 2
      end if
      !------------------
      ! fill in body info
      !------------------
      bod % use % state = spot% use% state
      bod % use % check = spot% use% check ! some checks are already applied...
      bod % mn          = 'nsoilm'
      call new_obs (obs, no, spot=spt)
      i1 = spt% o% i+1
      in = spt% o% i + spt% o% n
      obs % varno (i1 : in)          = VN_NSOILM       ! varno code
      obs %  olev (i1 : in)          = 0.01_wp         ! layer 0-2 cm below surface
      obs %  body (i1 : in)          = bod
      obs %  body (i1 : in)% lev_typ = VN_DEPTH        ! level type: below surface
      obs %  body (i1 : in)% o       = nssm  * 0.01_sp ! surface soil moisture
      obs %  body (i1 : in)% ac      = neesm * 0.01_sp ! estimated error in surface SM
      obs %  body (i1 : in)% eo      = neesm * 0.01_sp ! observation error
      if (l_errinflation) then
         obs% body(i1 : in)% eo      = obs% body(i1:in)% ac * 2.
      end if
!     obs %  body (i1 : in)% pcc     =                 ! << quality measure

      ! Derived normalized soil moisture will be filled in later...
      if (no == 2) then
         obs% varno (in)             = VN_SOILM        ! varno code
         obs% body  (in)% src        = SRC_DER
         obs% body  (in)% o          = rvind
         obs% body  (in)% ac         = -1._sp
         obs% body  (in)% eo         = -1._sp
         call decr_use (obs% body(i1)% use, STAT_PASSIVE, check=CHK_NOTUSED)
      end if

    else
      call decr_rpt_use (spot, CHK_INSDAT, STAT_DISMISS)
      lkeep = .false.
    endif

  end subroutine check_store_soil

!==============================================================================
end module mo_soil_obs
