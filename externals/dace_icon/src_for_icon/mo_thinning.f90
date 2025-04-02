!
!+ Observation thinning routines
!
MODULE mo_thinning
!
! Description:
!    Observation thinning routines.
!
! Current Maintainer: DWD, Robin Faulwetter, Harald Anlauf
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Hendrik Reich
!  Changes for rotated grids
! V1_7         2009/08/24 Andreas Rhodin
!  Change comments
! V1_8         2009/12/09 Harald Anlauf
!  check_domain: add check for observation above model top
! V1_9         2010/04/20 Andreas Rhodin
!  option to exclude observations near boundaries in out of domain check.
!  new namelist parameter: state, pref_retrv, pref_satids, pref_center
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  implement new parameters 'instr', 'pass'; extend checks to height_t,height_b
! V1_15        2011/12/06 Andreas Rhodin
!  cleanup
! V1_17        2011/12/21 Andreas Rhodin
!  increase ncodes (number of dbkz, satids etc) from 10 to 20
!  fix initialization of def_thin
! V1_19        2012-04-16 Andreas Rhodin
!  new namelist variable 'keep': keep report if state >= 'keep'
! V1_20        2012-06-18 Andreas Rhodin
!  cleanup: remove unused variables
! V1_22        2013-02-13 Andreas Rhodin
!  option to perform thinning before first guess check
! V1_26        2013/06/27 Harald Anlauf
!  Changes for GRIB2/GRIB_API/ICON
! V1_27        2013-11-08 Robin Faulwetter
!  Fixed undesired behaviour of thinning
! V1_29        2014/04/02 Andreas Rhodin
!  consistently use n_ot instead of n_obstype
! V1_42        2015-06-08 Harald Anlauf
!  implement non-global ICON grid
! V1_43        2015-08-19 Andreas Rhodin
!  skip out of domain check and TSK_SHRINK for COSMO-MEC
! V1_47        2016-06-06 Andreas Rhodin
!  apply out-of-domain-check to COSMO observation operators (for MEC)
! V1_48        2016-10-06 Robin Faulwetter
!  Implemented level-based thinning
! V1_50        2017-01-09 Robin Faulwetter
!  Restructured/unified cloud detection for radiances.
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! A.Rhodin, DWD 2003-2008
!===================================================================

  use mo_kind,       only: wp              ! working precision kind parameter
  use mo_namelist,   only: position_nml,  &! routine to position nml group
                           nnml,          &! namelist fortran unit number
                           POSITIONED      ! position_nml: OK    return flag
  use mo_t_datum,    only: rvind           ! invalid value
  use mo_t_obs,      only: t_obs,         &! observation data type
                           t_spot,        &! report meta data type
!                          COSMO,         &! COSMO module type
                           t_ilev,        &! information on levels to be thinned.
                           empty_ilev,    &!
                           operator(==),  &! comparison of t_ilev
                           match_ilev      ! check, whether data matches with t_thlev info
  use mo_mpi_dace,   only: dace,          &! MPI group info
                           p_alltoall,    &! generic MPI_alltoall(V) routine
                           p_bcast,       &! generic broadcast routine
                           p_sum           ! generic sum routine
  use mo_exception,  only: finish,        &! abort on error condition
                           message         ! print a warning
  use mo_t_use,      only: decr_use,      &! decrease the state of a datum
                           CHK_THIN,      &!
                           CHK_DOMAIN,    &!
                           CHK_AREA,      &!
                           CHK_HEIGHT,    &!
                           STAT_DISMISS,  &!
                           STAT_DEFAULT,  &!
                           STAT_INVALID,  &!
                           STAT_OBS_ONLY, &!
                           STAT_REJECTED, &!
                           STAT_PAS_REJ,  &!
!                          STAT_PASSIVE,  &!
                           stat_mnem,     &! mnemonics for states
                           stat_key        ! derive key from mnemonic
  use mo_time,       only: operator(-),   &! calculate time difference
                           hours           ! derive hours  from time data-type
  use mo_obs_tables, only: decr_rpt_use,  &! change use-flags of report
                           rept_use,      &! report type usage table
                           write_pending, &! write pending dismissed reports
                           obstyp          ! table of observation type mnem.
  use mo_atm_grid,   only: t_grid,        &! type definition for grid inform.
                           construct,     &! constructor routine
                           destruct,      &! destructor routine
!                          phirot2phi,    &! trafo rot->geo (phi)
                           phi2phirot,    &! trafo geo->rot (phi)
!                          rlarot2rla,    &! trafo rot->geo (lam)
                           rla2rlarot,    &! trafo geo->rot (lam)
                           construct_local ! local (rotated) lat-lon grid
! use mo_icon_grid,  only: dist_to_bound   ! determine distance to boundary
  use mo_wmo_tables, only: DWD6_ICON,     &! ICON global unstructured grid
                           DWD6_ICOSAHEDRON! Icosahedral based triangular grid
  use mo_satid,      only: id_mnem=>satid  ! derive satid from mnemonic
! use coordinates,   only: cartesian       !+++ required for IBM xlf V12.1
  use mo_grid_intpol,only: grid_indices,  &! derive grid indices of neighbors
                           mp              ! max. number of interpolation points
  use mo_run_params, only: ana_time,      &! analysis time
                           nproc1,nproc2   ! partitioning information
  use mo_system,     only: flush           ! Flush I/O buffer
  use mo_fdbk_tables,only: n_ot,          &! no of different report types
                           OT_AIREP,      &! AIREP observation type id
!                          OT_SYNOP,      &! SYNOP observation type id
!                          OT_SATEM,      &! SATEM observation type id
                           OT_TEMP,       &! TEMP  observation type id
                           OT_PILOT,      &! PILOT observation type id
                           OT_RAD,        &! RAD   observation type id
                           OT_RADAR,      &! RADAR observation type id
!                          VN_P,          &! pressure             code
                           VN_HEIGHT,     &! height               code
                           VN_HOSAG,      &! height above ground  code
                           VN_FLEV         ! nominal flight level code
  use mo_rad,        only: m_chan,        &! max. number of channels
                           n_set,         &! number of radiance datasets
                           set_indx,      &! get index of rad_set
                           rad_set,       &! rad_set array
                           n_set           ! number of radiance datasets
! use mo_instrid,    only: instr_rttov     ! RTTOV Id to WMO Id
  use mo_fg_checks,  only: check_suff      !
  use mo_dec_matrix, only: t_vector        ! vector data type
  use mo_tovs,       only: calc_tovs_quality!calc. quality of TOVS obs.

  implicit none

  private
  public :: thinning      ! new thinning routine
  public :: read_nml_thin ! read namelist /THINNING/
  public :: check_domain  ! check for out of domain
  !======================
  ! data type definitions
  !======================
  integer, parameter :: ncrit  =  8
  integer, parameter :: nterm  =  7
  integer, parameter :: ncodes = 20
  integer, parameter :: mrules = 100

  type t_crit
    character(len=8) :: name
    real(wp)         :: weight
    real(wp)         :: bound
  end type t_crit

  type t_thin
    character(len=40):: comment               ! comment
    integer          :: obstype               ! CMA observation type
    integer          :: codetype     (ncodes) ! CMA code type
    integer          :: dbkz         (ncodes) ! DWD Datenbankkennzahl
    integer          :: satid        (ncodes) ! satellite ID
    integer          :: instr        (ncodes) ! grid used for mapping
    integer          :: pref_satids(2,ncodes) ! preferred satids
    integer          :: pref_center(2,ncodes) ! preferred centers
    integer          :: pref_retrv (2,ncodes) ! preferred retrieval type
    integer          :: pref_grids (2,ncodes) ! preferred grids (i.e. assimilated instruments for TOVS)
    integer          :: ni                    ! GME resolution parameter
    integer          :: nlev                  ! number of levels
    real(wp)         :: dlev                  ! level increment (Pa)
    real(wp)         :: d_km                  ! horizontal resolution (km)
    real(wp)         :: plevel            (2) ! vertical pressure range (Pa)
    type (t_crit)    :: crit (ncrit)          ! criteria
    integer          :: rule (ncrit,nterm)    ! rules
    integer          :: state                 ! state for thinned observations
    integer          :: keep                  ! state to keep in any case
    integer          :: pass                  ! 0:apply before fg check, 1:after
    type(t_ilev)     :: levels       (m_chan) ! levels to be thinned
    type(t_ilev)     :: bands            (10) ! bands to be thinned
    integer          :: quality_mode          ! mode for quality calculation for thinning on selected levels
  end type t_thin

  integer            :: i                     ! Dummy index

  type (t_thin),save :: def_thin =                         &!
        t_thin (          '',                              &! comment
                          0,                               &! obstype
                       (/-1, (/( 0,i=2,ncodes)/) /),       &! codetype
                       (/-1, (/(-9,i=2,ncodes)/) /),       &! dbkz
                       (/-1, (/( 0,i=2,ncodes)/) /),       &! satid
                       (/-1, (/(-9,i=2,ncodes)/) /),       &! instr
                         -1,                               &! pref_satids
                         -1,                               &! pref_center
                         -1,                               &! pref_retrv
                         -1,                               &! pref_grids
                          0,                               &! ni
                          1,                               &! nlev
                          0.,                              &! dlev
                          0.,                              &! d_km
                          [ rvind, rvind],                 &! plevel range
                 (/t_crit('location',   1.,    0. ),       &! crit
                   t_crit('time    ',   1.,   -3. ),       &
                   t_crit('sequence',   1.,    0. ),       &
                   t_crit('data    ',   1.,    0. ),       &
                   t_crit('quality ',   1.,   -1. ),       &
                   t_crit('vertical',   1.,    0. ),       &
                   t_crit('status'  ,   1.,    0. ),       &
                   t_crit('pref    ', 100., -100. )/),     &
                  reshape((/0,0,0,0,0,0,1,0,               &! rules: status
                            0,0,0,0,0,0,0,1,               &!        preference
                            0,1,0,0,0,0,0,0,               &!        time
                            0,0,0,0,1,0,0,0,               &!        quality
                            1,0,0,0,0,1,0,0,               &!        distance
                            0,0,0,1,0,0,0,0,               &!        data
                            0,0,1,0,0,0,0,0/),             &! .. sequence
                          (/ncrit,nterm/)),                &!
                            STAT_DEFAULT,                  &! state
                            STAT_INVALID,                  &! keep
                          1,                               &! pass
                          (/(empty_ilev,i=1,m_chan)/),     &! levels
                          (/(empty_ilev,i=1,10    )/),     &! bands
                         -1                               ) ! quality_mode

  integer               :: nrules = 0
  type (t_thin) ,target :: thin_rules (mrules)
  save                  :: thin_rules

  type t_key
    integer  :: pe
    integer  :: box
    integer  :: spot
    real(wp) :: value (nterm)
    integer  :: n_lev
  end type t_key

  type t_keyx
    integer      :: idx(4)
    integer      :: k
    type (t_key) :: key
  end type t_keyx

!+++++++++++++++++++++++++++++++++++++++++++++++
! dont use generic interface for p_alltoall_keyx
! (workaround for INTEL ifort compiler)
!
! interface p_alltoall
!   module procedure p_alltoall_keyx
! end interface p_alltoall
!+++++++++++++++++++++++++++++++++++++++++++++++

!==============================================================================
contains
!==============================================================================
  subroutine print_thin_rules
  !---------------------------
  ! print thinning rules
  ! set by namelist /THINNING/
  !---------------------------
    integer :: i
    if (dace% lpio) then
      write(6,'(a)') repeat ('-',79)
      write(6,'(a)')
      write(6,'(a)') '  thinning rules:'
      do i=1,nrules
        call print_thin_rule (thin_rules(i), ruleno=i)
      end do
      write(6,'(a)')
    end if
  end subroutine print_thin_rules
!------------------------------------------------------------------------------
  subroutine print_thin_rule (thin, ruleno)
    type (t_thin) ,intent(in)           :: thin
    integer       ,intent(in), optional :: ruleno
    integer :: j, k
    real    :: d_gl
    write(6,'(a)')
    if (present (ruleno)) then
      write(6,'(i3,1x,a)')              ruleno,            trim   (thin% comment)
    else
      write(6,'(4x,a)')                                    trim   (thin% comment)
    end if
    write(6,'(a)')
    write(6,'(a,a)')               '    obstype      = ',  obstyp (thin% obstype)% name
    j = count (thin% codetype >= 0)
    if (thin% codetype(1) == -1) j = 0
    if (j > 0) &
      write(6,'(a,99i5)')          '    codetype     = ',          thin% codetype(:j)
    j = count (thin% dbkz >= 0)
    if (thin% dbkz(1) == -1) j = 0
    if (j > 0) &
      write(6,'(a,99i5)')          '    dbkz         = ',          thin% dbkz(:j)
    j = count (thin% satid > 0)
    if (thin% satid(1) == -1) j = 0
    if (j > 0) &
      write(6,'(a,99i5)')          '    satid        = ',          thin% satid(:j)
    j = count (thin% instr >= 0)
    if (thin% instr(1) == -1) j = 0
    if (j > 0) &
      write(6,'(a,99i5)')          '    instr        = ',          thin% instr(:j)
    j = count (thin% pref_satids >= 0)
    if (j > 0) then
      write(6,'(a,99i5,/)')        '    pref_satids  = ',          thin% pref_satids(1,:j)
      write(6,'(a,99i5,/)')        '    pref         = ',          thin% pref_satids(2,:j)
    endif
    j = count (thin% pref_center >= 0)
    if (j > 0) then
      write(6,'(a,99i5,/)')        '    pref_center  = ',          thin% pref_center(1,:j)
      write(6,'(a,99i5,/)')        '    pref         = ',          thin% pref_center(2,:j)
    endif
    j = count (thin% pref_retrv >= 0)
    if (j > 0) then
      write(6,'(a,99i5,/)')        '    pref_retrv   = ',          thin% pref_retrv (1,:j)
      write(6,'(a,99i5,/)')        '    pref         = ',          thin% pref_retrv (2,:j)
    endif
    j = count (thin% pref_grids >= 0)
    if (j > 0) then
      write(6,'(a,99i5,/)')        '    pref_grids   = ',          thin% pref_grids(1,:j)
      write(6,'(a,99i5,/)')        '    pref         = ',          thin% pref_grids(2,:j)
    endif
    d_gl = 0; if (thin% ni > 0) d_gl = 7054./thin% ni
    write(6,'(a,i5,f6.0,a)')       '    ni           = ',          thin% ni, d_gl, ' km'
    write(6,'(a,f7.1,a)')          '    d_km         = ',          thin% d_km, ' km'
    write(6,'(a,i5)')              '    nlev         = ',          thin% nlev
    write(6,'(a,i5)')              '    dlev(hPa)    = ',    nint (thin% dlev/100)
    if (thin% plevel(1) /= rvind) &
      write(6,'(a,2f7.1)')         '    plevel(hPa)  = ',          thin% plevel/100
    write(6,'(a,a)')               '    state        = ',stat_mnem(thin% state)
    write(6,'(a,a)')               '    keep         = ',stat_mnem(thin% keep)
    write(6,'(a,i5)')              '    pass         = ',          thin% pass
    j = count (thin% levels(:)%value > 0)
    if (j > 0) then
      if (any (thin% levels(:)%instr > 0)) then
        write(6,'(a,*(5(3x,I6,1x,"(instr ",I2.2,")"),/,19x))') &
                                   '    levels       = ',pack(     thin% levels,           &
                                                              mask=thin%levels(:)%value > 0)
      else
        write(6,'(a,*(20(1x,I6,:),/,19x))') &
                                   '    levels       = ',pack(     thin% levels(:)%value,  &
                                                              mask=thin%levels(:)%value > 0)
      end if
    end if
    j = count (thin% bands (:)%value > 0)
    if (j > 0) then
      if (any (thin% bands (:)%instr > 0)) then
        write(6,'(a,*(3x,I6,1x,"(instr ",I2.2,")"))') &
                                   '    bands        = ',pack(     thin% bands ,          &
                                                              mask=thin%bands(:)%value > 0)
      else
        write(6,'(a,*(1x,I6))')    '    bands        = ',pack(     thin% bands (:)%value, &
                                                              mask=thin%bands(:)%value > 0)
      end if
    end if
    if (thin% quality_mode >= 0) then
      write(6,'(a,i5)')            '    quality_mode = ',          thin% quality_mode
    end if
    write(6,'(a)')
    write(6,'(a)')                 '    criterium   weight     bound'
    do j=1,ncrit
      write(6,'(a,a,2f10.2)')'    ',thin% crit(j)% name,   &
                                    thin% crit(j)% weight, &
                                    thin% crit(j)% bound
    end do
    write(6,'(a)')
    do k=1,nterm
      if (any (thin% rule (:,k) /= 0)) then
        write(6,'(a,i2,a)',advance='no')'    rule ',k,' : '
        do j=1,ncrit
          if (thin% rule (j,k) /= 0) &
            write(6,'(1x,a)',advance='no') thin% crit(j)% name
        end do
        write(6,'()')
      endif
    end do
  end subroutine print_thin_rule
!------------------------------------------------------------------------------
  subroutine read_nml_thin

    !--------------------------------------------
    ! namelist variables declaration
    ! corresponding to components of type t_thin
    !--------------------------------------------
    character(len=8) :: obstype           ! CMA obstype
    integer          :: codetype (ncodes) ! CMA codetypes
    integer          :: dbkz     (ncodes) ! DWD Datenbankkennzahl
    integer          :: satid    (ncodes) ! satellite ID (integer)
    character(len=8) :: satids   (ncodes) ! satellite (string, overwrite satid)
    integer          :: instr    (ncodes) ! instrument used for mapping
    integer          :: ni                ! GME resolution parameter
    real(wp)         :: dlev              ! level increment (hPa)
    real(wp)         :: d_km              ! horizontal resolution (km)
    real(wp)         :: plevel   (2)      ! plevel range (hPa)
    character(len=40):: comment           ! any comment
    type (t_crit)    :: crit1             ! criterium 1
    type (t_crit)    :: crit2             ! criterium 2
    type (t_crit)    :: crit3             ! criterium 3
    type (t_crit)    :: crit4             ! criterium 4
    type (t_crit)    :: crit5             ! criterium 5
    type (t_crit)    :: crit6             ! criterium 6
    type (t_crit)    :: crit7             ! criterium 7
    type (t_crit)    :: crit8             ! criterium 8
    character(len=8) :: rule1 (ncrit)          ! rule 1
    character(len=8) :: rule2 (ncrit)          ! rule 2
    character(len=8) :: rule3 (ncrit)          ! rule 3
    character(len=8) :: rule4 (ncrit)          ! rule 4
    character(len=8) :: rule5 (ncrit)          ! rule 5
    character(len=8) :: rule6 (ncrit)          ! rule 6
    character(len=8) :: rule7 (ncrit)          ! rule 7
    integer          :: pref_satids (2,ncodes) ! satid,          preference
    integer          :: pref_center (2,ncodes) ! center,         preference
    integer          :: pref_retrv  (2,ncodes) ! retrieval type, preference
    integer          :: pref_grids  (2,ncodes) ! grid,           preference
    character(len=8) :: state                  ! state for thinned observations
    character(len=8) :: keep                   ! state to keep in any case
    integer          :: pass                   ! 0:apply before fg check, 1:after
    type(t_ilev)     :: levels(m_chan)         ! levels to be thinned
    type(t_ilev)     :: bands(10)              ! bands to be thinned (e.g. IASI band)
    integer          :: quality_mode           ! mode for quality calculation for thinning on selected levels

    namelist /THINNING/                                                       &
                    obstype, codetype, dbkz, satids, instr, ni, dlev, comment,&
                    crit1, crit2, crit3, crit4, crit5, crit6, crit7, crit8,   &
                    rule1, rule2, rule3, rule4, rule5, rule6, rule7,          &
                    pref_satids, pref_center, pref_retrv, pref_grids, state,  &
                    keep, pass, levels, bands, quality_mode, d_km, plevel
    !----------------------
    ! other local variables
    !----------------------
    logical ,save :: first = .true.
    integer       :: ierr
    integer       :: obst
    integer       :: i, l
#if defined(__GNUC__) && !defined(__G95__)
    character(len=80) :: iomsg          ! Fortran 2003 feature
    integer :: ios
#endif
#if defined(__ibm__)
    integer :: ios
#endif
    !----------------------------------------------
    ! include interfaces for external bcast routine
    !----------------------------------------------
    interface
      subroutine p_bcast_derivedtype (buffer, count, source, comm)
        import :: t_crit
        type(t_crit)               :: buffer    ! variable to bcast
        integer      ,intent(in)   :: count     ! len(byte) of variable
        integer      ,intent(in)   :: source    ! source processor index
        integer      ,intent(in)   :: comm      ! communicator
      end subroutine p_bcast_derivedtype
    end interface

    !------------------------------------------------
    ! read set of namelist group /THINNING/ only once
    !------------------------------------------------
    if (.not. first) return
    !------------------------------------
    ! repeatedly read namelist /THINNING/
    !------------------------------------
    do
      !----------------------------------------
      ! preset namelist variables with defaults
      !----------------------------------------
      comment      =  ''
      obstype      =  ''
      codetype     =  0; codetype(1) = -1
      dbkz         = -9; dbkz    (1) = -1
      satid        =  0; satid   (1) = -1
      satids       =  ''
      instr        = -9; instr   (1) = -1
      ni           =  0
      dlev         =  0._wp
      d_km         =  0._wp
      plevel       =  rvind
      crit1        =  t_crit (' ',0.,0.)
      crit2        =  t_crit (' ',0.,0.)
      crit3        =  t_crit (' ',0.,0.)
      crit4        =  t_crit (' ',0.,0.)
      crit5        =  t_crit (' ',0.,0.)
      crit6        =  t_crit (' ',0.,0.)
      crit7        =  t_crit (' ',0.,0.)
      crit8        =  t_crit (' ',0.,0.)
      rule1        =  ''
      rule2        =  ''
      rule3        =  ''
      rule4        =  ''
      rule5        =  ''
      rule6        =  ''
      rule7        =  ''
      pref_satids  =  -1
      pref_center  =  -1
      pref_retrv   =  -1
      pref_grids   =  -1
      state        =  ''
      keep         =  ''
      pass         =   1
      levels       =  t_ilev(-1,-1)
      bands        =  t_ilev(-1,-1)
      quality_mode =  -1
      !--------------
      ! read namelist
      !--------------
      if (dace% lpio) then
        call position_nml ('THINNING' ,lrewind=first ,status=ierr)
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=THINNING, iostat=ios)
          if (ios/=0) &
            call finish ('read_nml_thin','ERROR in namelist /THINNING/')
#elif defined(__GNUC__) && !defined(__G95__)
          ! gfortran: use Fortran 2003 feature iomsg
          iomsg = ""
          read (nnml, nml=THINNING, iostat=ios, iomsg=iomsg)
          if (ios/=0) then
             write (0,*) "iostat =", ios
             write (0,*) "iomsg  = ", trim (iomsg)
             call finish ('read_nml_thin','ERROR in namelist /THINNING/')
          end if
#else
          read (nnml ,nml=THINNING)
#endif
        end select
      endif
      first = .false.
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit
      !-------
      ! checks
      !-------
      if (dace% lpio) then
        if (any(codetype < 0)) then
          codetype (1 ) = -1
          codetype (2:) =  0
        endif
        if (any(dbkz == -1)) then
          dbkz (1 ) = -1
          dbkz (2:) = -9
        endif
        if (any(satids/='')) then
          satid = id_mnem (satids)
          do i=1,size(satid)
          if (satid(i)==0) &
            call finish('read_nml_thin','invalid satid '//satids(i))
          end do
          satid = max (satid, 0)
        endif
        if (any(satid < 0)) then
          satid (1 ) = -1
          satid (2:) =  0
        endif
        if (any(instr == -1)) then
          instr (1 ) = -1
          instr (2:) = -9
        endif
        obst = 0
        do i = 1, n_ot
          if (obstype == obstyp(i)% name) obst = i
        end do
        if (obst == 0) call finish ('namelist /THINNING/',      &
                                    'invalid obstype: '//obstype)
        dlev = dlev * 100._wp ! hPa -> Pa
        where (plevel /= rvind) plevel = plevel * 100._wp ! hPa -> Pa
      endif
      !-----------------------------------
      ! broadcast, store rule in structure
      !-----------------------------------
      l = size (transfer (crit1 ,(/' '/)))
      call p_bcast             (obst           ,dace% pio)
      call p_bcast             (codetype       ,dace% pio)
      call p_bcast             (dbkz           ,dace% pio)
      call p_bcast             (satid          ,dace% pio)
      call p_bcast             (instr          ,dace% pio)
      call p_bcast             (ni             ,dace% pio)
      call p_bcast             (dlev           ,dace% pio)
      call p_bcast             (d_km           ,dace% pio)
      call p_bcast             (plevel         ,dace% pio)
      call p_bcast             (comment        ,dace% pio)
      call p_bcast_derivedtype (crit1, l       ,dace% pio, dace% comm)
      call p_bcast_derivedtype (crit2, l       ,dace% pio, dace% comm)
      call p_bcast_derivedtype (crit3, l       ,dace% pio, dace% comm)
      call p_bcast_derivedtype (crit4, l       ,dace% pio, dace% comm)
      call p_bcast_derivedtype (crit5, l       ,dace% pio, dace% comm)
      call p_bcast_derivedtype (crit6, l       ,dace% pio, dace% comm)
      call p_bcast_derivedtype (crit7, l       ,dace% pio, dace% comm)
      call p_bcast_derivedtype (crit8, l       ,dace% pio, dace% comm)
      call p_bcast             (rule1          ,dace% pio)
      call p_bcast             (rule2          ,dace% pio)
      call p_bcast             (rule3          ,dace% pio)
      call p_bcast             (rule4          ,dace% pio)
      call p_bcast             (rule5          ,dace% pio)
      call p_bcast             (rule6          ,dace% pio)
      call p_bcast             (rule7          ,dace% pio)
      call p_bcast             (pref_satids    ,dace% pio)
      call p_bcast             (pref_center    ,dace% pio)
      call p_bcast             (pref_retrv     ,dace% pio)
      call p_bcast             (pref_grids     ,dace% pio)
      call p_bcast             (state          ,dace% pio)
      call p_bcast             (keep           ,dace% pio)
      call p_bcast             (pass           ,dace% pio)
      call p_bcast             (levels(:)%value,dace% pio)
      call p_bcast             (levels(:)%instr,dace% pio)
      call p_bcast             (bands(:)%value ,dace% pio)
      call p_bcast             (bands(:)%instr ,dace% pio)
      call p_bcast             (quality_mode   ,dace% pio)

      call new_thin_rule (                                                    &
                       obst, codetype, dbkz, satid, instr, ni, dlev, comment, &
                       crit1, crit2, crit3, crit4, crit5, crit6, crit7, crit8,&
                       rule1, rule2, rule3, rule4, rule5, rule6, rule7,       &
                       pref_satids, pref_center, pref_retrv, pref_grids,      &
                       state, keep, pass, levels, bands, quality_mode, d_km,  &
                       plevel)
    end do
    call print_thin_rules

  end subroutine read_nml_thin
!------------------------------------------------------------------------------
  subroutine new_thin_rule (                                                  &
                     obstype, codetype, dbkz, satid, instr, ni, dlev, comment,&
                     crit1, crit2, crit3, crit4, crit5, crit6, crit7, crit8,  &
                     rule1, rule2, rule3, rule4, rule5, rule6, rule7,         &
                     pref_satids, pref_center, pref_retrv, pref_grids,        &
                     state, keep, pass, levels, bands, quality_mode, d_km,    &
                     plevel)
  !-------------------------------
  ! set up a new rule for thinning
  !-------------------------------
  integer         ,intent(in)           :: obstype      ! CMA observation type
  integer         ,intent(in)           :: codetype (:) ! CMA code types
  integer         ,intent(in)           :: dbkz     (:) ! DWD Datenbankkennzahl
  integer         ,intent(in)           :: satid    (:) ! satellite ID
  integer         ,intent(in)           :: instr    (:) ! mapping instrument
  integer         ,intent(in)           :: ni           ! resolution parameter
  real(wp)        ,intent(in) ,optional :: dlev         ! level increment (Pa)
  character(len=*),intent(in) ,optional :: comment      ! any comment
  type(t_crit)    ,intent(in) ,optional :: crit1        ! criterium 1
  type(t_crit)    ,intent(in) ,optional :: crit2        ! criterium 2
  type(t_crit)    ,intent(in) ,optional :: crit3        ! criterium 3
  type(t_crit)    ,intent(in) ,optional :: crit4        ! criterium 4
  type(t_crit)    ,intent(in) ,optional :: crit5        ! criterium 5
  type(t_crit)    ,intent(in) ,optional :: crit6        ! criterium 6
  type(t_crit)    ,intent(in) ,optional :: crit7        ! criterium 7
  type(t_crit)    ,intent(in) ,optional :: crit8        ! criterium 8
  character(len=*),intent(in) ,optional :: rule1(ncrit)     ! rule 1
  character(len=*),intent(in) ,optional :: rule2(ncrit)     ! rule 2
  character(len=*),intent(in) ,optional :: rule3(ncrit)     ! rule 3
  character(len=*),intent(in) ,optional :: rule4(ncrit)     ! rule 4
  character(len=*),intent(in) ,optional :: rule5(ncrit)     ! rule 5
  character(len=*),intent(in) ,optional :: rule6(ncrit)     ! rule 6
  character(len=*),intent(in) ,optional :: rule7(ncrit)     ! rule 7
  integer         ,intent(in) ,optional :: pref_satids(:,:) ! preferred satids
  integer         ,intent(in) ,optional :: pref_center(:,:) ! preferred center
  integer         ,intent(in) ,optional :: pref_retrv (:,:) !  " retrieval type
  integer         ,intent(in) ,optional :: pref_grids (:,:) ! preferred grids
  character(len=*),intent(in) ,optional :: state            ! status
  character(len=*),intent(in) ,optional :: keep             ! status to keep
  integer         ,intent(in) ,optional :: pass             ! pass
  type(t_ilev)    ,intent(in) ,optional :: levels(:)        ! levels
  type(t_ilev)    ,intent(in) ,optional :: bands (:)        ! bands
  integer         ,intent(in) ,optional :: quality_mode     ! mode for quality calculation for thinning on selected levels
  real(wp)        ,intent(in) ,optional :: d_km         ! horizontal resolution (km)
  real(wp)        ,intent(in) ,optional :: plevel(:)    ! level range (Pa)

    !----------------
    ! local variables
    !----------------
    type (t_thin)    :: thin
    type (t_crit)    :: crit  (ncrit)
    character(len=8) :: rules (ncrit,nterm)
    integer          :: i, j

    !----------------------------
    ! process mandatory arguments
    !----------------------------
    thin           = def_thin
    thin% obstype  = obstype
    thin% codetype = codetype
    thin% dbkz     = dbkz
    thin% satid    = satid
    thin% instr    = instr
    thin% ni       = ni

    !---------------------------
    ! process optional arguments
    !---------------------------
    if (present (comment)) thin% comment = comment
    if (present (dlev   )) thin% dlev    = dlev
    if (thin% dlev > 0._wp) thin% nlev = max(1, nint(105000._wp/thin% dlev))
    if (present (d_km   )) thin% d_km    = d_km
    if (present (pref_satids)) &
      thin% pref_satids (:,1:size(pref_satids,2)) = pref_satids
    if (present (pref_center)) &
      thin% pref_center (:,1:size(pref_center,2)) = pref_center
    if (present (pref_retrv )) &
      thin% pref_retrv  (:,1:size(pref_retrv ,2)) = pref_retrv
    if (present (pref_grids)) &
      thin% pref_grids  (:,1:size(pref_grids ,2)) = pref_grids
    thin% state = rept_use (obstype)% use (CHK_THIN)
    if (pass == 0) thin% state = min (thin% state, STAT_DISMISS)
    if (present(state)) then
      if (state /= '')  thin% state = stat_key (state)
      if (thin% state < 0) &
        call finish('namelist /THINNING/','invalid state: '//state)
    endif
    if (present(keep)) then
      if (keep /= '')  thin% keep = stat_key (keep)
      if (thin% keep < 0) &
        call finish('namelist /THINNING/','invalid keep: '//keep)
    endif
    if (present(pass))  thin% pass  = pass
    crit(:)% name = ''
    if (present (crit1)) crit(1) = crit1
    if (present (crit2)) crit(2) = crit2
    if (present (crit3)) crit(3) = crit3
    if (present (crit4)) crit(4) = crit4
    if (present (crit5)) crit(5) = crit5
    if (present (crit6)) crit(6) = crit6
    if (present (crit7)) crit(7) = crit7
    if (present (crit8)) crit(8) = crit8
    do i=1,ncrit
      where  (crit(i)% name == thin%crit% name)
        thin%crit% bound  = crit(i)% bound
        thin%crit% weight = crit(i)% weight
      endwhere
      if (all(crit(i)% name /= thin%crit% name .and.                     &
              crit(i)% name /= ''                    ))                  &
        call finish ('namelist /THINNING/','invalid name: '//crit(i)% name)
    end do

    rules = ''
    if (present (rule1)) rules(:,1) = rule1
    if (present (rule2)) rules(:,2) = rule2
    if (present (rule3)) rules(:,3) = rule3
    if (present (rule4)) rules(:,4) = rule4
    if (present (rule5)) rules(:,5) = rule5
    if (present (rule6)) rules(:,6) = rule6
    if (present (rule7)) rules(:,7) = rule7

    do i=1,nterm
      if (rules(1,i) /= '') then
        thin% rule (:,i) = 0
        do j=1,ncrit
          if (any(rules(:,i) == thin% crit(j)% name)) thin% rule (j,i) = 1
          if (all(rules(:,i) /= thin% crit(j)% name .and.             &
                  rules(:,i) /= ''                  .and.             &
                  rules(:,i) /= 'none'))                              &
            call finish ('namelist /THINNING/','invalid name in rules')
        end do
      end if
    end do

    if (present(levels)) then
      j = 0
      do i = 1, size(levels)
        if (levels(i)%value > 0) then
          j = j + 1
          if (j > size(thin%levels)) then
            call finish ('read_nml_thin','increase size od thin%levels')
          end if
          thin%levels(j) = levels(i)
        end if
      end do
    end if
    if (present(bands)) then
      j = 0
      do i = 1, size(bands)
        if (bands(i)%value > 0) then
          j = j + 1
          if (j > size(thin%bands)) then
            call finish ('read_nml_thin','increase size od thin%bands')
          end if
          thin%bands(j) = bands(i)
        end if
      end do
    end if
    if (present(quality_mode)) then
      thin% quality_mode = quality_mode
    end if

    if (present (plevel)) then
       if (size (plevel) /= 2) call finish ('read_nml_thin','bad size of plevel')
       thin% plevel = plevel
    end if

    !-----------------------------------------
    ! issue warning if rule is already present
    !-----------------------------------------
    if (dace% lpio) then
      do i = 1, nrules
        if (     thin_rules(i)% obstype  == thin% obstype .and. &
            all (thin_rules(i)% codetype == thin% codetype.and. &
                 thin_rules(i)% dbkz     == thin% dbkz    .and. &
                 thin_rules(i)% satid    == thin% satid   .and. &
                 thin_rules(i)% instr    == thin% instr   .and. &
                 thin_rules(i)% state    == thin% state   .and. &
                 thin_rules(i)% ni       == thin% ni     ).and. &
            all (thin_rules(i)% plevel   == thin% plevel ).and. &
            all (thin_rules(i)% levels   == thin% levels ).and. &
            all (thin_rules(i)% bands    == thin% bands  )) then
          write (6,*) 'WARNING: multiple thinning rules given !'
          write (0,*) 'WARNING: multiple thinning rules given !'
!         thin_rules(i) = thin
!         return
        endif
      end do
    endif
    !--------------------------------------------
    ! issue warning if both ni and d_km are given
    !--------------------------------------------
    if (dace% lpio .and. ni > 0 .and. d_km > 0) then
      write (6,*) 'WARNING: ni and d_km specified, ni takes precedence: ', &
           trim (thin% comment)
      write (0,*) 'WARNING: ni and d_km specified, ni takes precedence: ', &
           trim (thin% comment)
    end if
    !-----------------------------
    ! else insert new rule in list
    !-----------------------------
    nrules = nrules + 1
    if (nrules > mrules) call finish ('read_nml_thin','increase mrules')
    thin_rules(nrules) = thin

  end subroutine new_thin_rule
!------------------------------------------------------------------------------
  function cmp_keys (key1, key2, rule) result (use2)
  type (t_key), intent(in) :: key1
  type (t_key), intent(in) :: key2
  integer,      intent(in) :: rule
  logical                  :: use2
    integer :: i
    do i = 1, nterm
      if (key1% value(i) == key2% value(i)) cycle
      use2 = key2% value(i) > key1% value(i)
      return
    end do
    write(0,*)
    write(0,*) "cmp_keys: rule:", rule
    write(0,*) "cmp_keys: spots", key1% spot, key2% spot
    do i = 1, nterm
      write(0,*) 'cmp_keys: values =',key1% value(i) !, key2% value(i)
    end do
!   call finish  ('cmp_keys','equal scores')
    call message ('cmp_keys','equal scores')
    write(0,*)
    call flush (0)
    use2 = .false.
  end function cmp_keys
!------------------------------------------------------------------------------
  subroutine thinning (obs, fg_grid, pass, e_bg, e_o)
  type (t_obs)   ,intent(inout)        :: obs (:) ! observations
  type (t_grid)  ,intent(in)           :: fg_grid ! first guess grid
  integer        ,intent(in)           :: pass    ! 0/1 before/after fg check
  type(t_vector) ,intent(in), optional :: e_bg    ! background error
  type(t_vector) ,intent(in), optional :: e_o     ! observation error

  !-----------------
  ! thinning routine
  !-----------------
    target                  :: obs
    integer                 :: ir       ! rule index
    integer                 :: ib       ! box  index
    integer                 :: is       ! spot index
    type(t_grid)            :: grid     ! hexagonal grid descriptor
    type(t_thin) ,pointer   :: th       ! pointer to actual thinning rule
    type(t_obs)  ,pointer   :: oi       ! pointer to observation box
    type(t_obs)  ,pointer   :: oi_      ! pointer to observation box
    type(t_spot) ,pointer   :: si       ! pointer to report
    integer                 :: np       ! number of processed observations
    integer                 :: nr       ! number of rejected  observations
    integer                 :: na       ! number of currently active obsv.
    real(wp)                :: dp       ! normalised pressure difference
    type (t_key)            :: key2     ! sorting (thinning) key
    type (t_key) ,pointer   :: key1     ! sorting (thinning) key
    integer ,allocatable    :: ikey(:,:,:,:) ! key field (on grid)
    integer                 :: ik       ! key index
    real(wp)                :: v(ncrit) ! sorting (thinning) quantity
    logical                 :: lrej     ! flag for rejection
    integer                 :: i,j,k,d,l! indices
    integer                 :: iset,n,nc! used for band to level transformation
    integer                 :: idx(mp,4)! grid indices of neighbors (pts,indx)
!   integer                 :: ml(1)    ! maxloc
    integer                 :: nn       ! number of neighbours
    real(wp)                :: w(mp)    ! weights
    real(wp)                :: w1,wd,wdm! used for gridpoint choice
    integer                 :: i1       ! index of gridpoint with max weight
    integer                 :: idx1(4)  ! index of nearest gridpoint
    integer                 :: statold  ! status of report before rejection
    integer                 :: statrej  ! status of report after  rejection
    integer                 :: n_lev    ! number of levels to be thinned in report
    integer ,allocatable    :: ilev(:)  ! levels indices to be thinned in report
    real(wp), parameter     :: tol = 1.E-12 ! tolerance for equal gridpoint weights
    character(10)           :: cthin    ! comment prefix for thinning rule
    !------------------------------
    ! arrays used for communication
    !------------------------------
    type (t_keyx) ,pointer  :: keyx (:) ! keys in observation space, unsorted
    type (t_keyx) ,pointer  :: keyo (:) ! keys in observation space, sorted
    type (t_keyx) ,pointer  :: tmp  (:) ! temporary
    type (t_keyx) ,pointer  :: keym (:) ! keys in model space, sorted
    integer                 :: snd_cnt (0:dace% npe-1) ! no. keys sent
    integer                 :: snd_dsp (0:dace% npe-1) ! start indices sent
    integer                 :: rcv_cnt (0:dace% npe-1) ! no. keys received
    integer                 ::     dsp (0:dace% npe-1) ! key indices
    logical       ,pointer  :: rejm (:) ! mark rejected keys in modelspace
    logical       ,pointer  :: rejo (:) !      rejected keys in obsspace
    integer                 :: pe       ! processor element index
    integer                 :: nam      ! keys accepted, in model space
    integer                 :: nlevs

    allocate (keyx(1000))
    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write (6,'(a)')    repeat('-',79)
      write (6,'(a)')
      write (6,'(a,i2)') '  thinning: pass =',pass
    endif
    !-------------------------
    ! loop over thinning rules
    !-------------------------
    do ir = 1, nrules
      th => thin_rules (ir)
      if (th% pass /= pass) cycle
      if (th% ni == 0 .and. fg_grid% gridtype == DWD6_ICOSAHEDRON) &
          th% ni = fg_grid% ni
      if (fg_grid% global) then
         if (th% ni <= 0  ) cycle
      end if
      !-----------------------
      ! set up data structures
      !-----------------------
      if (fg_grid% global .or. th% ni > 0) then
         call construct (grid, gridtype=DWD6_ICOSAHEDRON, ni=th% ni,  &
                         nproc1=nproc1, nproc2=nproc2, comm=dace% comm)
      else
         call construct_local (grid, fg_grid, d_km=th% d_km)
      end if
      allocate (ikey (grid% lbg(1):grid% ubg(1), &
                      grid% lbg(2):grid% ubg(2), th% nlev, grid% nd))
      ikey    = 0
      statrej = th% state
      nlevs   = count (th% levels(:)% value > 0)

      if (th% obstype == OT_RAD) then
        if (any(th% bands(:)% value > 0)) then
          ! Convert band numbers into channel numbers
          l = count(th% levels(:)% value > 0)
          if (l > 0) th% levels(1:l) = pack(th% levels(:), mask=th% levels(:)% value > 0)
          do i = 1, size(th% bands)
            if (th% bands(i)% value > 0) then
              if (th% instr(1) < 0) then
                call finish ('namelist /THINNING/', '"th% bands" require "instr"')
              end if
              if (th% bands(i)% instr >= 0) then
                k = th% bands(i)% instr
              else
                k = th%instr(1)
              end if
              if (th% satid(1) >= 0) then
                iset = set_indx(rad_set(1:n_set), satid=th%satid(1), grid=k)
              else
                iset = set_indx(rad_set(1:n_set), grid=k)
              end if
              if (iset > 0) then
                if (.not.associated(rad_set(iset)%band)) then
                  call finish ('namelist /THINNING/', 'no band information available.')
                end if
                nc = rad_set(iset)%n_chan
                n  = count(rad_set(iset)%band(1:nc) == th% bands(i)%value)
                if (l + n > size(th% levels)) then
                  call finish ('namelist /THINNING/', 'not enough th% levels')
                end if
                th% levels(l+1:l+n)%value = pack(rad_set(iset)%chan(1:nc), mask=rad_set(iset)%band(1:nc)==th% bands(i)%value)
                if (th% bands(i)% instr >= 0) th% levels(l+1:l+n)% instr = th% bands(i)% instr
                l = l + n
              end if
            end if
          end do
          th% bands(:) = t_ilev(-1,-1)
        end if
        nlevs = count(th%levels(:)%value > 0)
        call calc_tovs_quality(1, e_bg=e_bg, e_o=e_o, instrs=th%instr(:), satids=th%satid(:), &
                               tl=th%levels(1:nlevs), obs=obs, qmode=th%quality_mode)
      else
        nlevs = count(th%levels(:)%value > 0)
      end if

      !---------
      ! printout
      !---------
      if (dace% lpio) then
        call print_thin_rule (th, ruleno=ir)
      endif
      write (cthin,'("thin: ",i0)') ir
      !=====================
      ! part I: local checks
      !=====================
      !-----------------------
      ! loop over observations
      !-----------------------
      np = 0
      nr = 0
      na = 0
      do ib = 1, size(obs)
        if (obs(ib)% pe /= dace% pe) cycle
        oi => obs(ib)
        do is = 1, oi% n_spot
          si  => oi% spot (is)
          oi_ => oi
          !--------------------------------------------
          ! filter observations applicable to this rule
          !--------------------------------------------
          if (    th% obstype     /= si% hd% obstype  ) cycle
          if (    th% codetype(1) /= -1          .and.&
              all(th% codetype    /= si% hd% codetype)) cycle
          if (    th% dbkz    (1) /= -1          .and.&
              all(th% dbkz        /= si% hd% dbkz    )) cycle
          if (    th% satid   (1) /= -1          .and.&
              all(th% satid       /= si% hd% satid   )) cycle
          if (    th% instr   (1) /= -1          .and.&
              all(th% instr       /= si% hd% grid_id )) cycle
          if (        statrej     >= si% use% state   ) cycle
          if (    th% plevel  (1) /= rvind       .and.&
                  si% ps          /= rvind       .and.&
                 (th% plevel  (1) <  si% ps .or. &
                  th% plevel  (2) >= si% ps      )    ) cycle
          np = np + 1
          !------------------------------------------------------------------------
          ! check whether the whole report or only selected levels shall be thinned
          !------------------------------------------------------------------------
          n_lev = 0
          if (nlevs > 0) then
            allocate(ilev(si%o%n))
            call match_ilev(th%levels(1:nlevs), n_lev, si=si, oi=oi_, ind=ilev)
          end if
          !---------------------------
          ! find neighbour grid-points
          !---------------------------
          call Grid_Indices   &
            (si% col% c% dlon,& ! <-- geodetic longitude
             si% col% c% dlat,& ! <-- geodetic latitude
             grid,            & ! <-- grid data type
             idx,             & ! --> Grid point indices [Point, index]
             w,               & ! --> Weight
             nn)                ! --> number of points returned
!         if (nn/=3) call finish('thinning','nn /= 3')
          if (nn < 3 .or. nn > 4) call finish('thinning','nn /= 3, 4')
          !-------------------
          ! Choose a gridpoint
          !-------------------
          ! In order to obtain portable results, it is important to make sure that
          ! the choice is not affected by numeric noise. Therefore, we do not rely on
          ! tiny differences (<tol) in w. Weights with differences < tol are considered
          ! "equal". In case the two biggest weights are "equal", we have to select one
          ! of these gridpoints. For a lat/lon grid, we select the first grid point.
          ! For the GME grid we have to take into account a potential cyclic shift of
          ! the gridpoints. We should select the same grid point for different cyclic
          ! shifts. We arbitrarily select the first gridpoint in the pair of neighbored
          ! gridpoints with biggest weights.
#if 1  /* Implementation 1 */
          w1  = w(1)
          i1  = 1
          wdm = 0._wp
          do i = 2, nn
            wd = w(i) - w1
            wdm = max(wdm, abs(wd))
            if (wd > tol) then
              w1  = w(i)
              i1  = i
            elseif (abs(wd) <= tol) then
              ! Weights are "equal"
              if (nn == 4) then
                ! Lat/lon grid. Do nothing, keep gridpoint with smallest index i.
              else
                ! GME grid
                if (i1 == i-1) then
                  ! Keep choice of previous (neighbored) gridpoint
                else
                  w1 = w(i)
                  i1 = i
                end if
              end if
            end if
          end do
          i = i1
#else  /* Alternative implementation of the basic idea.  */
          w1  = maxval (w(1:nn))        ! Actual max. weight
          wdm = minval (w(1:nn))        ! Actual min. weight
          wdm = w1 - wdm                ! Difference max. minus min. weight
          wd  = w1 - tol                ! Threshold weight
          select case (nn)
          case default                  ! Lat/Lon grid
             do i = 1, 4
                if (w(i) > wd) exit     ! Find first weight exceeding treshold
             end do
          case (3)                      ! GME grid
             do i = 2, 3
                if (w(i) > wd) exit     ! Find first weight exceeding treshold
             end do
             select case (i)
             case (2)
                i = 1   ! Points 1 & 2 have close to equal weights, use point 1
             case (3)
                i = 3   ! Points 1 & 3 have close to equal weights, use point 3
             case default
                i = 1   ! Point 1 has definitely largest weight, keep
             end select
          end select
#endif
          ! Very unlikely case: all weights are "equal"
          if (wdm <= tol) then
            i = maxloc(idx(1:nn,1)*1000000._wp + idx(1:nn,2)*1000._wp + idx(1:nn,3), 1)
          end if
          !ml = maxloc (w(1:nn))
          !i = ml(1)

          !--------------------------------------------
          ! move nearest grid-point into first position
          !--------------------------------------------
          if (i/=1) then
            w1    = w (i)
            w (i) = w (1)
            w (1) = w1
            idx1      = idx (i,:)
            idx (i,:) = idx (1,:)
            idx (1,:) = idx1
          endif
          !--------------------
          ! find vertical index
          !--------------------
          k  = 1
          dp = 1._wp
          if (th% nlev > 1) then
            dp = si% ps / th% dlev
            k  = nint (dp)
            dp = 1._wp - 2._wp * abs(dp - k)
            k  = min (k+1, th% nlev)
          endif
          !----------------------------------------
          ! derive quantities required for thinning
          !----------------------------------------
          v(1) = max(w(1)-maxval(w(2:nn)),0._wp)            ! location
          v(2) = - abs (hours (si% actual_time - ana_time)) ! time
!         v(3) = si% hd% id                                 ! sequence
!         v(3) = si% id                                     ! sequence
          if (si% hd% id > 0) then
            v(3) = si% hd% id + 1.e10_wp * si% hd% source   ! sequence
          else
            v(3) = si% id     + 1.e10_wp * si% hd% source   ! sequence
          end if
          v(4) = si% o% n                                   ! data
          v(5) = si% pcc / 100._wp                          ! quality
          v(6) = dp                                         ! vertical
          v(7) = si% use% state                             ! status
          if (n_lev > 0) then
            v(4) = n_lev
            v(7) = 0
            do i = 1, n_lev
              j = ilev(i)
              if (oi_% body(j)% use% state /= STAT_REJECTED) then
                v(7) = max(int(v(7)), int(oi_% body(j)% use% state))
              else
                v(7) = max(int(v(7)), STAT_PAS_REJ)
              end if
            end do
            select case(si% hd% obstype)
            case(OT_RAD)
              call calc_tovs_quality(2, e_bg=e_bg, e_o=e_o, instrs=th%instr(:), satids=th%satid(:), &
                               tl=th%levels(1:nlevs), oi=oi_, si=si, qual=v(5), qmode=th%quality_mode)
              oi_% body(ilev(1:n_lev))% pcc = v(5)
            case default
              if (n_lev > 0) then
                v(5) = sum(oi_% body(ilev(1:n_lev))% pcc) / (1.*n_lev)
              else
                v(5) = 0._wp
              end if
            end select
            v(5) = v(5) / 100._wp
          end if
          if (nlevs > 0) deallocate(ilev)

          ! REJECTED observations should not be preferred over PASSIVE observations
          if (v(7) == STAT_REJECTED) v(7) = STAT_PAS_REJ
          v(8) = 0._wp                                      ! preference
          !--------------------------------------
          ! set preferences for satids or centers
          !--------------------------------------
          do l = 1, ncodes
            if (th% pref_satids(1,l) < 0) exit
            if (th% pref_satids(1,l) == si% hd% satid) &
              v(8) = v(8) + th% pref_satids(2,l)
          end do
          do l = 1, ncodes
            if (th% pref_center(1,l) < 0) exit
            if (th% pref_center(1,l) == si% hd% center) &
              v(8) = v(8) + th% pref_center(2,l)
          end do
          do l = 1, ncodes
            if (th% pref_retrv (1,l) < 0) exit
            if (th% pref_retrv (1,l) == si% stret) &
              v(8) = v(8) + th% pref_retrv (2,l)
          end do
          do l = 1, ncodes
            if (th% pref_grids(1,l) < 0) exit
            if (th% pref_grids(1,l) == si% hd% grid_id) &
              v(8) = v(8) + th% pref_grids(2,l)
          end do
          !-----------------
          ! check for bounds
          !-----------------
          lrej = any (v(1:ncrit) <  th% crit(1:ncrit)% bound)
          if (lrej) then
            statold = si% use% state
!            call decr_rpt_use (si, CHK_THIN, use=th%state, comment='bound')
            call decr_use_thin(comment='bound')
            if (statold /= si% use% state) nr = nr + 1
            cycle
          endif
          !-----------
          ! derive key
          !-----------
          key2% pe    = dace% pe
          key2% box   = ib
          key2% spot  = is
          key2% n_lev = n_lev
          key2% value = 0._wp
          do l = 1, nterm
            key2% value(l) = sum (th% rule(:,l) * th% crit(:)% weight * v(:))
          end do
          !---------------
          ! derive indices
          !---------------
          i    =  idx (1,1) ! horizontal index 1
          j    =  idx (1,2) ! horizontal index 2
          d    =  idx (1,3) ! diamond    index
          ik   =  ikey (i,j,k,d)
          !-----------------------------------------
          ! first obs at this gridpoint, just insert
          !-----------------------------------------
          if (ik == 0) then
            if (na >= size(keyx)) then
              tmp => keyx
              allocate (keyx(2*na))
              keyx(:na) = tmp(:na)
              deallocate (tmp)
            endif
            na = na + 1
            ikey (i,j,k,d) = na
            keyx(na)% key  = key2
            keyx(na)% idx  = idx (1,:)
            keyx(na)% k    = k
            cycle
          endif
          !-------------
          ! compare keys
          !-------------
          key1 => keyx(ik)% key
          if (cmp_keys (key1, key2, ir)) then             ! if (use key 2)
            si   => obs(key1% box)% spot (key1% spot)     ! reject  key 1
            oi_  => obs(key1% box)
            n_lev=  key1%n_lev
            key1 =  key2                                  ! keep key 2
          endif
          if (si% use% state >= th% keep) cycle
!          call decr_rpt_use (si, CHK_THIN, use=th%state, comment='key')
          call decr_use_thin(comment='key')
        end do
      end do
      !=========================
      ! Part II: parallel checks
      !=========================
      if (dace% npe > 1) then
        !--------------------------------
        ! reorder keys for sending to PEs
        !--------------------------------
        snd_cnt = 0
        do i=1,na
          pe = keyx(i)% idx(4)
          snd_cnt (pe) = snd_cnt (pe) + 1
        end do
        snd_dsp (0) = 1
        do i = 1, dace% npe-1
          snd_dsp (i) = snd_dsp (i-1) + snd_cnt (i-1)
        end do
        dsp = snd_dsp
        allocate (keyo(na))
        allocate (rejo(na))
        do i=1,na
          pe = keyx(i)% idx(4)
          ik = dsp (pe)
          keyo (ik) = keyx(i)
          dsp (pe) = ik + 1
        end do
        !--------------------
        ! send / receive data
        !--------------------
        call p_alltoall (snd_cnt, rcv_cnt)
        nam = sum (rcv_cnt)
        allocate (keym (nam))
        allocate (rejm (nam))
        rejm = .false.
        call p_alltoall_keyx (keyo, keym, sendcounts=snd_cnt, &
                                          recvcounts=rcv_cnt)
        !-------------
        ! compare keys
        !-------------
        ikey = 0
        do is = 1, nam
          !---------------
          ! derive indices
          !---------------
          i    =  keym(is)% idx (1) ! horizontal index 1
          j    =  keym(is)% idx (2) ! horizontal index 2
          d    =  keym(is)% idx (3) ! diamond    index
          k    =  keym(is)% k
          key2 =  keym(is)% key
          ik   =  ikey (i,j,k,d)
          !-----------------------------------------
          ! first obs at this gridpoint, just insert
          !-----------------------------------------
          if (ik == 0) then
            ikey (i,j,k,d) = is
            cycle
          endif
          !-------------
          ! compare keys
          !-------------
          key1 => keym(ik)% key
          if (cmp_keys (key1, key2, ir)) then
            rejm(ik) = .true.
            ikey (i,j,k,d) = is
          else
            rejm(is) = .true.
          endif
        end do
        !------------------
        ! send results back
        !------------------
        call p_alltoall (rejm, rejo, sendcounts=rcv_cnt, recvcounts=snd_cnt)
        !---------------------
        ! dismiss observations
        !---------------------
        do i= 1, na
          if (rejo(i)) then
            si    => obs(keyo(i)% key% box)% spot (keyo(i)% key% spot)
            oi_   => obs(keyo(i)% key% box)
            n_lev =  keyo(i)% key% n_lev
            if (si% use% state >= th% keep) cycle
!            call decr_rpt_use (si, CHK_THIN, use=th%state, comment='key')
            call decr_use_thin(comment='key')
          endif
        end do
        !--------
        ! cleanup
        !--------
        deallocate (keyo)
        deallocate (keym)
        deallocate (rejo)
        deallocate (rejm)
      endif
      !---------
      ! printout
      !---------
      np = p_sum (np)
      nr = p_sum (nr)
      if (dace% lpio) then
        write (6,'()')
        write (6,'(a,i8)') '  processed =', np
        write (6,'(a,i8)') '  rejected  =', nr
        write (6,'(a,i8)') '  accepted  =', np - nr
      endif
      !---------------------
      ! free data structures
      !---------------------
      call destruct (grid)
      deallocate (ikey)
      call write_pending

      if (nlevs > 0) then
         call check_suff(obs)
      end if
    end do
    deallocate (keyx)

  contains

    subroutine decr_use_thin(comment)
      character(len=*), intent(in)    :: comment

      integer, allocatable :: ind(:)
      integer              :: i, j
      logical              :: l_rej

      if (n_lev > 0) then
        allocate(ind(si%o%n))
        call match_ilev(th%levels(1:nlevs), n_lev, si=si, oi=oi_, ind=ind)
        l_rej = .false.
        do i = 1, n_lev
          j = ind(i)
          statold = oi_% body(j)% use% state
          call decr_use (oi_% body(j)% use, th%state, CHK_THIN)
          if (statold /= oi_% body(j)% use% state) l_rej = .true.
        end do
        if (l_rej) nr = nr + 1
      else
        statold = si% use% state
        call decr_rpt_use (si, CHK_THIN, use=th%state,           &
                           comment=trim (cthin) // " " // comment)
        if (statold /= si% use% state) nr = nr + 1
      end if

    end subroutine decr_use_thin

  end subroutine thinning
!==============================================================================
  subroutine check_domain (obs, grid, local, horizontal, lcheck_ptopf)
  type(t_obs)  ,intent(inout)        :: obs          ! observations
  type(t_grid) ,intent(in)           :: grid         ! grid (to derive domain)
  logical      ,intent(in) ,optional :: local        ! no parallel mode ?
  logical      ,intent(in) ,optional :: horizontal   ! check horizontal domain?
  logical      ,intent(in) ,optional :: lcheck_ptopf ! discard obs for plev < ptopf?
  !----------------------------------
  ! check for out of domain condition
  !----------------------------------

    integer               :: is                     ! spot index
    real(wp)              :: lonw, lone, lats, latn ! domain boundaries
    real(wp)              :: dlon, dlat             ! longitude, latitude
    real(wp)              :: dist                   ! distance to boundary (degree)
    integer               :: c_ctrl                 ! distance to boundary (grid-points)
    type(t_spot) ,pointer :: s                      ! spot pointer
    integer               :: np, na, ne, nd, nc     ! processed/accepted
    logical               :: global                 ! global grid?
    logical               :: lhor                   ! do horizontal check
    integer               :: i1, in, il             ! level indices
    integer               :: ng, nt, nb, nr         ! above top, bottom; outside
    real(wp)              :: excl_bnd               ! exclude obs.at boundaries
    real(wp)              :: height_t               ! top    height [hPa]
    real(wp)              :: height_b               ! bottom height [hPa]
    real(wp)              :: height_top             ! top    height [m]
    real(wp)              :: height_bot             ! bottom height [m]
    logical               :: l                      ! no parallel mode
    integer               :: stat_height            ! status for height check
    integer               :: stat_domain            ! status for above domain
    integer               :: idx(3,4)               ! Grid indices
    real(wp)              :: w  (3)                 ! Weights
    integer               :: npts                   ! no. of points
    logical               :: llcheck_ptopf          ! by default discard plev<ptopf
    integer               :: d_b(3)                 ! ICON-LAM: boundary distance
    real(wp)              :: pobs                   ! "pressure" of observation
    real(wp)              :: zobs                   ! "height"   of observation
    logical               :: lpobs

    global =  (grid% gridtype == DWD6_ICOSAHEDRON)                     &
         .or. (grid% gridtype == DWD6_ICON         .and. grid% global) &
         .or. (grid% cyc_x .and. grid% poly)
    if (present (horizontal)) global = global .and. horizontal
    llcheck_ptopf = .true.; if (present (lcheck_ptopf)) llcheck_ptopf = lcheck_ptopf
    lhor = .not.global    ; if (present (horizontal))   lhor          = lhor .and. horizontal
    l    = .false.        ; if (present (local))        l             = local

    !-------------------------
    ! derive domain boundaries
    !-------------------------
    lonw = -huge(1._wp)
    lone =  huge(1._wp)
    lats = -huge(1._wp)
    latn =  huge(1._wp)
    if (lhor) then
       if (grid% gridtype == DWD6_ICON) then
          lonw = minval (grid% dlon)
          lone = maxval (grid% dlon)
          lats = minval (grid% dlat)
          latn = maxval (grid% dlat)
          !------------------------------------------
          ! count distance of grid-points to boundary
          !------------------------------------------
#ifdef _CRAYFTN
#warning "Code disabled due to Cray bug with INTENT (Case #214995)"
#endif
!         if (any (rept_use(:)% excl_bnd > 0._wp) .and.                &
!              .not. allocated (grid% icongrid% patch% cells% c_ctrl)) &
!           call dist_to_bound (grid% icongrid% patch, 50)
       else
          if (.not. grid% cyc_x) then
             lonw = grid% dlon(1)
             lone = grid% dlon(grid% nx)
          endif
          if (.not. grid% poly) then
             lats = grid% dlat(1)
             latn = grid% dlat(grid% ny)
          endif
       end if
    end if
    if (lonw >= lone) then
       write(0,*) "lonw,lone=", lonw,lone
       call finish ('check_domain','lonw >= lone')
    end if
    if (lats >= latn) then
       write(0,*) "lats,latn=", lats,latn
       call finish ('check_domain','lats >= latn')
    end if

    !-----------------------
    ! loop over spots, check
    !-----------------------
    np = 0
    na = 0
    ng = 0
    nt = 0
    nb = 0
    ne = 0
    nd = 0
    nc = 0
    nr = 0
    do is = 1, obs% n_spot
      np = np + 1
      s => obs% spot(is)
      !---------------------------------------------------
      ! Height check (conventional observation types only)
      !---------------------------------------------------
!     if(s% hd% modtype==COSMO) then
!       nc = nc + 1
!       cycle
!     endif
      select case (s% hd% obstype)
!     case (OT_SYNOP:OT_SATEM)
      case (OT_TEMP,OT_PILOT,OT_AIREP)
        height_t    = rept_use(s% hd% obstype)% height_t * 100._wp ! hPa -> Pa
        height_b    = rept_use(s% hd% obstype)% height_b * 100._wp ! hPa -> Pa
        height_top  = rept_use(s% hd% obstype)% height_top         ! m
        height_bot  = rept_use(s% hd% obstype)% height_bot         ! m
        stat_height = rept_use(s% hd% obstype)% use (CHK_HEIGHT)
        stat_domain = rept_use(s% hd% obstype)% use (CHK_DOMAIN)
        stat_domain = max (stat_domain, STAT_DISMISS)
        i1 = s% o% i + 1
        in = s% o% i + s% o% n
        do il = i1, in
          !---------------------------------------
          ! Prefer using plev when it has been set (for now)
          !---------------------------------------
          if (obs% body(il)% plev > 0.) then
             pobs  = obs% body(il)% plev
             lpobs = .true.
          else
             lpobs = .false.
             select case (obs% body(il)% lev_typ)
             case (VN_HEIGHT, VN_FLEV)
                zobs  = obs% olev(il)
             case (VN_HOSAG)
                zobs  = obs% olev(il) + s% z
             case default                      ! includes VN_P and -1 (~ VN_P)
                pobs  = obs% body(il)% plev
                lpobs = .true.
             end select
          end if

          if (lpobs) then
             if (pobs < grid% ptopf .and. llcheck_ptopf) then
               ng = ng + 1
               call decr_use (obs% body(il)% use, stat_domain , CHK_DOMAIN)
             else if (pobs < height_t) then
               nt = nt + 1
               call decr_use (obs% body(il)% use, stat_height,  CHK_HEIGHT)
             else if (pobs > height_b) then
               nb = nb + 1
               call decr_use (obs% body(il)% use, stat_height,  CHK_HEIGHT)
             end if
          else
             if (zobs > grid% htopf .and. llcheck_ptopf) then
               ng = ng + 1
               call decr_use (obs% body(il)% use, stat_domain , CHK_DOMAIN)
             else if (zobs > height_top) then
               nt = nt + 1
               call decr_use (obs% body(il)% use, stat_height,  CHK_HEIGHT)
             else if (zobs < height_bot) then
               nb = nb + 1
               call decr_use (obs% body(il)% use, stat_height,  CHK_HEIGHT)
             end if
          end if
        end do
      end select
      !------------------------------------------------
      ! Remote sensing observations: the actual station
      ! may be located outside the domain (e.g. radar)
      !------------------------------------------------
      if (lhor .and. grid% gridtype == DWD6_ICON) then
         select case (s% hd% obstype)
         case (OT_RADAR)
            stat_domain = rept_use(s% hd% obstype)% use (CHK_DOMAIN)
            stat_domain = min (stat_domain, STAT_OBS_ONLY)
            i1 = s% o% i + 1
            in = s% o% i + s% o% n
            do il = i1, in
               if (obs% body(il)% use% state <= stat_domain) cycle
               dlon = obs% body(il)% lon
               dlat = obs% body(il)% lat
               !---------------------------------------------
               ! Locate observation within the irregular grid
               !---------------------------------------------
               call Grid_Indices      &
                    (dlon,            & ! <-- geodetic longitude
                     dlat,            & ! <-- geodetic latitude
                     grid,            & ! <-- grid data type
                     idx,             & ! --> Grid point indices [Point, index]
                     w,               & ! --> Weights
                     npts,            & ! --> number of points returned
                     lunique = .false.) ! <-- set unique index
               if (npts == 0) then
                  !-----------------------------
                  ! observation is out of domain
                  !-----------------------------
                  nr = nr + 1
                  call decr_use (obs% body(il)% use, stat_domain, CHK_DOMAIN)
                  cycle
               end if
            end do
         end select
      end if
      !---------------------------------
      ! Area check (not for global grid)
      !---------------------------------
      if (lhor) then
         excl_bnd = rept_use(s% hd% obstype)% excl_bnd
         !-----------------------------------
         ! special treatment for rotated grid
         !-----------------------------------
         if (grid% rot) then
            dlon = rla2rlarot(s% col% c% dlat, &
                              s% col% c% dlon, &
                              grid% dlatr,     &
                              grid% dlonr,     &
                              0._wp            )
            dlat = phi2phirot(s% col% c% dlat, &
                              s% col% c% dlon, &
                              grid% dlatr,     &
                              grid% dlonr      )
         else
            dlon = s% col% c% dlon
            dlat = s% col% c% dlat
         end if

         if (dlon < lonw      )  dlon = dlon + 360._wp
         if (dlon < lonw .or. &
             dlon > lone .or. &
             dlat < lats .or. &
             dlat > latn      )  then
            nd = nd + 1
            call decr_rpt_use (s, CHK_DOMAIN)
            cycle
         end if

         select case (grid% gridtype)
         case default
            if (dlon < lonw + excl_bnd .or. &
                dlon > lone - excl_bnd .or. &
                dlat < lats + excl_bnd .or. &
                dlat > latn - excl_bnd      )  then
               ne = ne + 1
               call decr_rpt_use (s, CHK_AREA)
               cycle
            end if
         case (DWD6_ICON)
            !------------------------------------------------
            ! for ICON-LAM / ICON-NEST locate the observation
            ! on the irregular grid
            !------------------------------------------------
            call Grid_Indices      &
                 (dlon,            & ! <-- geodetic longitude
                  dlat,            & ! <-- geodetic latitude
                  grid,            & ! <-- grid data type
                  idx,             & ! --> Grid point indices [Point, index]
                  w,               & ! --> Weights
                  npts,            & ! --> number of points returned
                  lunique = .false.) ! <-- set unique index
            if (npts == 0) then
               !-----------------------------
               ! observation is out of domain
               !-----------------------------
               nd = nd + 1
               call decr_rpt_use (s, CHK_DOMAIN)
               cycle
            else
               !----------------------------------------------
               ! check if observation is too close to boundary
               !----------------------------------------------
              if (allocated (grid% icongrid% patch% cells% c_ctrl)) then
                d_b(1) = grid% icongrid% patch% cells% c_ctrl(idx(1,1),idx(1,2))
                d_b(2) = grid% icongrid% patch% cells% c_ctrl(idx(2,1),idx(2,2))
                d_b(3) = grid% icongrid% patch% cells% c_ctrl(idx(3,1),idx(3,2))
                c_ctrl = minval (d_b)
!if (c_ctrl > 0 .and. c_ctrl <= 3) then
!write(0,*) dace%pe, 'check_domain: ### ',s%statid, real(dlat), real(dlon), d_b
!end if
                if (c_ctrl > 0) then
                  dist  = c_ctrl * grid% d_deg * sqrt(2._wp)
                  if (dist < excl_bnd) then
                     ne = ne + 1
                     call decr_rpt_use (s, CHK_AREA)
                     cycle
                  end if
                endif
              endif
            endif
         end select
      end if
      na = na + 1
    end do

    !-------------
    ! print report
    !-------------
    if (.not. l) then
      np = p_sum (np)
      na = p_sum (na)
      ng = p_sum (ng)
      nt = p_sum (nt)
      nb = p_sum (nb)
      ne = p_sum (ne)
      nd = p_sum (nd)
!     nc = p_sum (nc)
      nr = p_sum (nr)
    endif
    if (l .or. dace% lpio) then
      write(6,'( )')
      write(6,'(a)')        '  check for out-of-domain condition'
      write(6,'(a,l1)')     '    horizontal check = ',lhor
      write(6,'( )')
      if (lhor) then
      write(6,'(a,2f10.2)') '    lat (rotated coordinates)  = ',lats,latn
      write(6,'(a,2f10.2)') '    lon (rotated coordinates)  = ',lonw,lone
      end if
      write(6,'(a,f10.2)')  '    model top [hPa]            = ',grid% ptopf/100._wp
      write(6,'( )')
      write(6,'(a,i7)')     '    reports processed          = ',np
!     write(6,'(a,i7)')     '    reports skipped (COSMO)    = ',nc
      write(6,'(a,i7)')     '    reports accepted           = ',na
      write(6,'(a,i7)')     '    reports out of domain      = ',nd
      write(6,'(a,i7)')     '    reports excluded at bounds = ',ne
      write(6,'( )')
      write(6,'(a,l1)')     '    lcheck_ptopf               = ',llcheck_ptopf
      write(6,'(a,i7)')     '    levels above grid          = ',ng
      write(6,'(a,i7)')     '    levels above top           = ',nt
      write(6,'(a,i7)')     '    levels below bot           = ',nb
      write(6,'(a,i7)')     '    out of domain: radar obs.  = ',nr
      write(6,'( )')
      write(6,'(a)')        repeat('-',79)
    endif

    !--------------------------------------
    ! abort in case of no reports in domain
    !--------------------------------------
!   if (na==0) call finish('check_domain','no reports accepted')
    if (np - nc > 0 .and. na == 0) then
       if (dace% lpio) call message ('check_domain','no reports accepted')
    end if

  end subroutine check_domain
!==============================================================================
#define DERIVED type(t_keyx)
#undef MPI_TYPE
#define p_alltoall_DERIVED p_alltoall_keyx
#include "p_alltoall_derived.incf"
#undef DERIVED
!==============================================================================
end module mo_thinning
