!
!+ set up complex rules for observation data usage
!
MODULE mo_obs_rules
!
! Description:
!   This module allows to build up rules for observation operators
!   in dependence of observation type (3D-Var-type, BUFR-type, BUFR-subtype,
!   databank-'Kennzahl'), location (region defined by longitude/latitude,
!   height defined by pressure), parameter (temperature, wind, ...)
!
!   For each observation type default values are coded here and may be
!   overwritten using the namelist group /RULES/.
!
! Current Maintainer: DWD, Harald Anlauf, Alexander Cress
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Harald Anlauf
!  #ifdef USE_SEQUENCE
! V1_7         2009/08/24 Andreas Rhodin
!  new option for artificial data generation PRC_ADD_1: (add 1 to observation)
! V1_8         2009/12/09 Andreas Rhodin
!  add PRC_ADD_R, PRC_SUB_R (add,subtract obs.err. to observation)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  new rules:fr_land,_sea,phase,zenith_angle,zlim,upper bound on observed values
! V1_17        2011/12/21 Andreas Rhodin
!  increase number of satids in namelist /rules/ from 20 to 30
! V1_20        2012-06-18 Harald Anlauf
!  get_rule: properly handle asurf/bsurf
! V1_22        2013-02-13 Andreas Rhodin
!  make 'iu1' public
! V1_23        2013-03-26 Andreas Rhodin
!  namelist /RULES/ add parameters set%v_loc, %h_loc, %ekf_pass for LETKF
! V1_27        2013-11-08 Andreas Rhodin
!  increase number of dbkz entries in /rules/ to 50
! V1_28        2014/02/26 Robin Faulwetter
!  Introduced tsurf in /RULES/.
!  Bugfixes in MW cloud check, /RULES/ usage for radiances, bias correction.
!  Introduced grtim_*_POL.
!  obs - fg based cloud checks for AMSUA and ATMS.
!  New features for diagnosing the most important bias correction predictors
! V1_31        2014-08-21 Andreas Rhodin
!  /RULES/: allow for zlim(2) < zlim(1) for rules valid outside this range
! V1_35        2014-11-07 Harald Anlauf
!  Change iu1 to actually integer(i1); adaptions for specific compilers
! V1_37        2014-12-23 Robin Faulwetter
!  Option for rules that shall only be applied if a instrument is missing
! V1_47        2016-06-06 Andreas Rhodin
!  Add obstype as selection criterium in namelist /rules/.
!  Disentangle /RULES/ default settings for ICON/COSMO/MEC.
! V1_48        2016-10-06 Andreas Rhodin
!  minor cosmetic change
! V1_50        2017-01-09 Robin Faulwetter
!  option for rules: be applied if certain levels/channels are missing
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2003-2008  original source
! Gerhard Paul    DWD  2008       SKY data base
! Harald Anlauf   DWD  2008       optimizations for SX8
!==============================================================================
  !=============
  ! modules used
  !=============
  use mo_kind,       only: wp, sp, i1, i2    ! kind parameters
  use mo_exception,  only: finish            ! abort routine
  use mo_namelist,   only: position_nml,    &! routine to position nml group
                           nnml,            &! namelist fortran unit number
                           POSITIONED        ! position_nml: OK    return flag
  use mo_mpi_dace,   only: dace,            &! DACE communication info
                           p_bcast           ! broadcast routine
  use mo_t_use,      only: stat_key          ! id for 'states'
  use mo_run_params, only: ana_time          ! analysis time
  use mo_time,       only: ihh               ! derive hh (hours) from time
  use mo_satid,      only: satid             ! derive satellite id from name
  use mo_physics,    only: t0c               ! 273.15 K
  use mo_fdbk_tables,only: init_fdbk_tables,&! initialise tables
              tab_varno => varno             ! variable number table
  use mo_t_table,    only: value_name,      &! derive table entry from name
                           name_value        ! find name of table entry
  implicit none

  !================
  ! public entities
  !================
  private
  public :: read_nml_rules ! read namelist /RULES/
  public :: print_rules    ! print the rules
  public :: new_rule       ! set a rules
  public :: get_rule       ! get a rule
  public :: t_set          ! data type definition
  public :: empty_set      ! empty data type
  public :: nc             ! number of 'channels' accounted for
  public :: rud            ! undefined real    value
  public :: iud, iu1       ! undefined integer values
  public :: PRC_NONE, PRC_ADD_1, PRC_ADDFG, PRC_ADDOBSER, PRC_ADDOBSE2, &
            PRC_OBS_FG, PRC_ADD_R, PRC_SUB_R
  !--------------------------------------
  ! Instrument "level" information
  !--------------------------------------
  public :: t_ilev
  public :: empty_ilev
  public :: operator(==)
  !---------------------
  ! constants definition
  !---------------------
  integer    ,parameter :: iud    = -1000000000     ! undefined value
  real(sp)   ,parameter :: rud    = -1000000000._sp ! undefined value
  integer(i1),parameter :: iu1    = -127            ! undefined for integer(i1)
  integer(i2),parameter :: iu2    = -32767          ! undefined for integer(i2)
  integer    ,parameter :: nkz    = 60              ! number of 'Datenbankkennzif.'
! integer    ,parameter :: nc     = 50              ! number of channels, constit.
  integer    ,parameter :: nc     =  1              !+++ currently not used +++
  integer    ,parameter :: nsat   = 30              ! number of Satellite ids
  integer    ,parameter :: nv     = 10              ! number of varnos
  integer    ,parameter :: ninstr = 16              ! number of instrument ids
  integer    ,parameter :: ninstm = 5               ! number of missing instrument ids
  integer    ,parameter :: nlev   = 5               ! number of channels
  !----------------
  ! processing flag
  !----------------
  integer(i1) ,parameter :: PRC_NONE     =  0 ! none
  integer(i1) ,parameter :: PRC_ADD_1    =  1 ! add 1 to observation
  integer(i1) ,parameter :: PRC_ADDFG    =  2 ! add first guess (obs.is increm.)
  integer(i1) ,parameter :: PRC_ADDOBSER =  4 ! add random observational error
  integer(i1) ,parameter :: PRC_ADDOBSE2 =  8 ! add random obs.err. / sqrt(2)
  integer(i1) ,parameter :: PRC_OBS_FG   = 16 ! set observation to fg
  integer(i1) ,parameter :: PRC_ADD_R    = 32 ! add obs.err. to observation
  integer(i1) ,parameter :: PRC_SUB_R    = 64 ! subtract obs.err. from obs.

  !======================
  ! Data type definitions
  !======================
  !------------------------------------
  ! settings specific for one parameter
  !------------------------------------
  type t_set
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer     :: use       = iud   ! how to use this parameter
    real(sp)    :: sgm_oc(4) = rud   ! observation error check bound
                                     ! reject:  |o_err|        > sgm_oc(1)
                                     !          |o_err/o|      > sgm_oc(2)
                                     !          |o_err/fg_err| > sgm_oc(3)
                                     ! enlarge: |o_err/fg_err| < sgm_oc(4)
    real(sp)    :: sgm_fg(3) = rud   ! first guess check bound
    real(sp)    :: bnd_fg(2) = rud   ! bounds on first guess values
    real(sp)    :: bnd_obs   = rud   ! upper bound on observed values
    real(sp)    :: sgm_vq    = rud   ! variational quality control bound
    integer     :: frm_vq    = iud   ! var. quality control formulation
    real(sp)    :: v_loc     = 0._sp ! specific vertical   localisation scale
    real(sp)    :: h_loc     = 0._sp ! specific horizontal localisation scale
    integer(i2) :: m_rej     = iu2   ! for variational quality control
    integer(i1) :: prc       = iu1   ! processing flag
    integer(i1) :: ekf_pass  = 0     ! LETKF pass
  end type t_set

  !------------------------------------------------------------
  ! Instrument "level" information, e.g. channel identification
  !------------------------------------------------------------
  type t_ilev
    integer          :: value  = -1  ! level, e.g. channel
    integer          :: instr  = -1  ! instrument
  end type t_ilev
  type (t_ilev),parameter :: empty_ilev  = t_ilev(-1, -1)

  !--------------------------
  ! data type to hold 'rules'
  !--------------------------
  type t_rules
    character(len=64) :: comment = ''
    integer           :: modtype       = iud        ! DACE module type
    integer           :: obstype       = iud        ! observation type
    integer           :: codetype      = iud        ! code type
    integer           :: bf_type       = iud        ! BUFR-message-type
    integer           :: bf_subt       = iud        ! BUFR-message-subtype
    integer           :: center        = iud        ! Processing center
    integer           :: db_kz (nkz)   = iud        ! Databank-'Kennzahl'
    integer           :: satid(nsat)   = iud        ! Satellite id
    integer           :: varno(nv)     = iud        ! varno integer
    logical           :: l_stname      = .false.    ! station name check active
    character(len=10) :: stname        = ''         ! station name

    real(wp)    :: lat       (2)       =  (/ rud, -rud /) ! Latitude  bounds
    real(wp)    :: lon       (2)       =  (/ rud, -rud /) ! Longitude bounds
    real(wp)    :: plim      (2)       =  (/ rud, -rud /) ! Pressure  bounds [Pa]
    real(wp)    :: zlim      (2)       =  (/ rud, -rud /) ! Height    bounds [m]
    real(wp)    :: mdzlim    (2)       =  (/ rud, -rud /) ! Model orography b.[m]
    real(wp)    :: sol_zenith(2)       =  (/ rud, -rud /) ! Solar zenith angle [deg]
    real(wp)    :: sat_zenith(2)       =  (/ rud, -rud /) ! Satellite zenith angle [deg]
    real(wp)    :: bsurf               =     rud          ! Pressure below surface
    real(wp)    :: asurf               =     rud          ! Pressure above surface
    real(wp)    :: fr_land             =     rud          ! Land fraction required
    real(wp)    :: fr_sea              =     rud          ! Sea  fraction required
    real(wp)    :: tsurf     (2)       =  (/ rud, -rud /) ! Surface temp. bounds
    integer     :: phase     (2)       =  (/ iud,  iud /) ! Phase or flag range
    integer     :: instrs     (ninstr) =  iud             ! Instruments
    integer     :: instrs_miss(ninstm) =  iud             ! Instruments that are missing
    type(t_ilev):: levs_miss(nlev)     = empty_ilev       ! Levels that are missing
                                                          ! (defined by instr and chan number)
    logical     :: xlonlat             =  .false.         ! Exclude area
    logical     :: rlon                =  .false.         ! Restrict longitudes
    logical     :: rlat                =  .false.         ! Restrict latitudes
    logical     :: lhour               =  .true.          ! Rule currently active?

    integer     :: use                 = iud   ! how to use (USE_NOT, _ASS, _MON)
    integer     :: verb                = iud   ! verbosity flag (diagnostic printout)
    integer     :: msl                 = iud   ! reduce pressure to mean sea level

    type(t_set) :: t                     ! temperature specific settings
    type(t_set) :: gp                    ! geopotential
    type(t_set) :: uv                    ! wind
    type(t_set) :: q                     ! humidity
    type(t_set) :: p                     ! pressure
    type(t_set) :: o                     ! other
    type(t_set) :: c (nc)                ! channels, constituents ...
  end type t_rules

!------------------------------------------------------------------------------
  interface operator (==)
    module procedure ilev_eq
  end interface
!------------------------------------------------------------------------------

  !=================
  ! module variables
  !=================
  type (t_rules) ,save :: ru (100)     ! storage for rules
  integer        ,save :: n_ru = 0     ! number of rules

  type (t_rules) ,save :: empty        ! 'empty' rule

  type (t_set) ,parameter :: &
    empty_set = t_set (iud,rud,rud,rud,rud,rud,iud,0._sp,0._sp,iu2,iu1,0)

contains
!------------------------------------------------------------------------------
  subroutine print_rules

   integer          :: i, nz, ns, ni, nim, ncm, j
   integer          :: nvv
!  integer          :: j
!  character(len=4) :: c

   if (.not. dace% lpio) return

   write (6,'(a)') repeat ('-',79)
   write (6,'()')
   write (6,'(a)') ' Rules for observation processing'

   do i=1,n_ru

     nz  = count (ru (i)% db_kz       /= iud)
     ns  = count (ru (i)% satid       /= iud)
     nvv = count (ru (i)% varno       /= iud)
     ni  = count (ru (i)% instrs      /= iud)
     nim = count (ru (i)% instrs_miss /= iud)
     ncm = count (ru (i)% levs_miss(:)% value > 0 .and. ru (i)% levs_miss(:)% instr > 0)

     write (6,'()')
     write (6,'(i3,1x,2a)') i,'Filter: ',trim(ru (i)% comment)
     write (6,'()')
     write (6,'(4x,a,i5)')         'DACE module type     = ',ru (i)% modtype
     write (6,'(4x,a,i5)')         'observation type     = ',ru (i)% obstype
     write (6,'(4x,a,i5)')         'code type            = ',ru (i)% codetype
     write (6,'(4x,a,i5)')         'BUFR-message-type    = ',ru (i)% bf_type
     write (6,'(4x,a,i5)')         'BUFR-message-subtype = ',ru (i)% bf_subt
     write (6,'(4x,a,i5)')         'Processing center    = ',ru (i)% center
     write (6,'(4x,a,(10i6):/(26x,10i6))')&
                                   'Databank-Kennzahl    =' ,ru (i)% db_kz(:nz)
     write (6,'(4x,a,(10i6):/(26x,10i6))')&
                                   'Satid                =' ,ru (i)% satid(:ns)
     if (ru (i)% varno(1) /= iud) &
      write (6,'(4x,a,*(1x,a12))') &
            'Varno                =', name_value (tab_varno, ru (i)% varno(:nvv))
     if (ru (i)% stname /= '') &
      write(6,'(4x,a,a)')          'Station Name         = ',ru (i)% stname
     write (6,'()')
     write (6,'(4x,a)')            'Optional Filter:'
     write (6,'(4x,a,2(1x,f10.2))')'Latitude  bounds     = ',ru (i)% lat
     write (6,'(4x,a,2(1x,f10.2))')'Longitude bounds     = ',ru (i)% lon
     write (6,'(4x,a,2(1x,l1))')   'Exclude   area       = ',ru (i)% xlonlat
     write (6,'(4x,a,2(1x,l1))')   'Rule active          = ',ru (i)% lhour
     write (6,'(4x,a,2(1x,f10.2))')'Pressure  bounds     = ',ru (i)% plim
     write (6,'(4x,a,2(1x,f10.2))')'Height    bounds [m] = ',ru (i)% zlim
     write (6,'(4x,a,2(1x,f10.2))')'Model orography b.[m]= ',ru (i)% mdzlim
     write (6,'(4x,a,2(1x,f10.2))')'Solar zenith angle   = ',ru (i)% sol_zenith
     write (6,'(4x,a,2(1x,f10.2))')'Sat. zenith angle    = ',ru (i)% sat_zenith
     write (6,'(4x,a,  1x,f10.2 )')'Pressure  below surf.= ',ru (i)% bsurf
     write (6,'(4x,a,  1x,f10.2 )')'Pressure  above surf.= ',ru (i)% asurf
     write (6,'(4x,a,  1x,f10.2 )')'Land fraction        = ',ru (i)% fr_land
     write (6,'(4x,a,  1x,f10.2 )')'Sea  fraction        = ',ru (i)% fr_sea
     write (6,'(4x,a,2(1x,f10.2))')'Surface temp.        = ',ru (i)% tsurf
     write (6,'(4x,a,2(1x,i10))')  'Phase or flag range  = ',ru (i)% phase
     write (6,'(4x,a,(10i5):/(27x,10i5))')&
                                   'Instruments          = ',ru (i)% instrs     (:ni )
     write (6,'(4x,a,(10i5):/(27x,10i5))')&
                                   'Instruments missing  = ',ru (i)% instrs_miss(:nim)
     if (ncm > 0) then
       write (6,'(4x,a,(2i5))')    'Levels missing       = ',ru (i)% levs_miss(1)
       do j = 2, ncm
         write (6,'(27x,(2i5))')                             ru (i)% levs_miss(j)
       end do
     else
       write (6,'(4x,a)')          'Levels missing       = '
     end if

     write (6,'()')
     write (6,'(4x,a)')            'Flags:'
     write (6,'(4x,a,i5)')         'Use flag             = ',ru (i)% use
     write (6,'(4x,a,i5)')         'Verbosity flag       = ',ru (i)% verb
     write (6,'()')
     write (6,'(4x,a)')            'Parameter specific:'
     write (6,'(1x,6x,2(1x,a5),11(1x,a10),1x,a6,1x,a3,2(1x,a10),1x,a6)')   &
       'use','m_rej', 'sgm_o_abs','sgm_o_rel','sgm_o_fgh','sgm_o_fgl',     &
       'bnd_fg_lw','bnd_fg_up','bnd_obs', 'sgm_fg_p',          &
       'sgm_fg_a','sgm_fg_eo','sgm_vq','frm_vq','prc','v_loc','h_loc','ekfpas'
     call pr_spec (ru (i)% t ,'t')
     call pr_spec (ru (i)% gp,'gp')
     call pr_spec (ru (i)% uv,'uv')
     call pr_spec (ru (i)% q ,'q')
     call pr_spec (ru (i)% p ,'p')
     call pr_spec (ru (i)% o ,'o')
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! comment loop on zero size array for IBM/xlf12.1.0.x .
! IBM/xlf12.1.0.x: The value of the DO-loop increment should be negative when
!                  initial value is greater than the terminal value.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     do j = 1, nc
!       write (c,'("c",i3.3)') j
!       call pr_spec (ru (i)% c(j), c)
!     end do
   end do
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine pr_spec (s,c)
   type (t_set)    , intent(in) :: s
   character(len=*), intent(in) :: c
     if (s% use      /= iud   .or. &
     any(s% sgm_oc   /= rud)  .or. &
     any(s% sgm_fg   /= rud)  .or. &
     any(s% bnd_fg   /= rud)  .or. &
         s% bnd_obs  /= rud   .or. &
         s% sgm_vq   /= rud   .or. &
         s% frm_vq   /= iud   .or. &
         s% v_loc    /= 0._sp .or. &
         s% h_loc    /= 0._sp .or. &
         s% m_rej    /= iu2   .or. &
         s% ekf_pass /= 0     .or. &
         s% prc      /= iu1        ) then
!#if defined (__sun)
!       !------------------------------------------------------------
!       ! Temporary work around for FP exception with formatted write
!       !------------------------------------------------------------
!       write (6,'(4x,a6,2(1x,i5),12(1x,es10.3),1x,i6,i4)') c, s% use, s% m_rej, &
!         s% sgm_oc, s% bnd_fg, s% bnd_obs, s% sgm_fg, s% sgm_vq, s% frm_vq, s% prc
!#else
       write (6,'(1x,a6,2(1x,i5),11(1x,f10.3),1x,i6,1x,i3,2(1x,es10.3),1x,i6)') &
         c, s% use, s% m_rej, s% sgm_oc, s% bnd_fg, s% bnd_obs, s% sgm_fg,      &
         s% sgm_vq, s% frm_vq, s% prc, s%v_loc, s%h_loc, s%ekf_pass
!#endif
     endif
   end subroutine pr_spec
  end subroutine print_rules
!------------------------------------------------------------------------------
  subroutine new_rule (comment,                                               &
                       type, obstype, codetype, bf_type,bf_subt,db_kz, center,&
                       stname, satids, lat, lon, plim, zlim, mdzlim, varnos,  &
                       sol_zenith, sat_zenith, bsurf,  asurf, hours, fr_land, &
                       fr_sea, tsurf, phase, instrs, instrs_miss, levs_miss,  &
                       use, verb, msl, t, gp, uv, q, p, o, c                  )

  character(len=*) ,intent(in) ,optional :: comment
  integer          ,intent(in) ,optional :: type
  integer          ,intent(in) ,optional :: obstype
  integer          ,intent(in) ,optional :: codetype
  integer          ,intent(in) ,optional :: bf_type
  integer          ,intent(in) ,optional :: bf_subt
  integer          ,intent(in) ,optional :: db_kz (:)
  integer          ,intent(in) ,optional :: center
  character(len=*) ,intent(in) ,optional :: stname
  character(len=*) ,intent(in) ,optional :: satids(:)
  character(len=*) ,intent(in) ,optional :: varnos(:)

  real(wp)         ,intent(in) ,optional :: lat        (2)
  real(wp)         ,intent(in) ,optional :: lon        (2)
  real(wp)         ,intent(in) ,optional :: plim       (2)
  real(wp)         ,intent(in) ,optional :: zlim       (2)
  real(wp)         ,intent(in) ,optional :: mdzlim     (2)
  real(wp)         ,intent(in) ,optional :: sol_zenith (2)
  real(wp)         ,intent(in) ,optional :: sat_zenith (2)
  real(wp)         ,intent(in) ,optional :: bsurf
  real(wp)         ,intent(in) ,optional :: asurf
  real(wp)         ,intent(in) ,optional :: fr_land
  real(wp)         ,intent(in) ,optional :: fr_sea
  real(wp)         ,intent(in) ,optional :: tsurf      (2)
  integer          ,intent(in) ,optional :: phase      (2)
  integer          ,intent(in) ,optional :: instrs     (:)
  integer          ,intent(in) ,optional :: instrs_miss(:)
  type(t_ilev)     ,intent(in) ,optional :: levs_miss  (:)
  integer          ,intent(in) ,optional :: hours(:)

  integer          ,intent(in) ,optional :: use
  integer          ,intent(in) ,optional :: verb
  integer          ,intent(in) ,optional :: msl
  type(t_set)      ,intent(in) ,optional :: t
  type(t_set)      ,intent(in) ,optional :: gp
  type(t_set)      ,intent(in) ,optional :: uv
  type(t_set)      ,intent(in) ,optional :: q
  type(t_set)      ,intent(in) ,optional :: p
  type(t_set)      ,intent(in) ,optional :: o
  type(t_set)      ,intent(in) ,optional :: c (:)

    integer :: n, l

    if (present (db_kz)) then
      if (size(db_kz) > size(ru(1)% db_kz)) then
        if (present (comment)) write (0,*) "new_rule: ", trim (comment)
        write (0,*) "new_rule: size(db_kz), nkz =", size(db_kz), nkz
        call finish ('new_rule','size(db_kz) > nkz')
      end if
    endif

    n_ru = n_ru + 1
    if (n_ru > size (ru)) &
      call finish ('new_rule','max number of rules exhausted')
    n = n_ru
    ru(n) = empty
    if (present (comment))  ru(n)% comment  = comment
    if (present (type   ))  ru(n)% modtype  = type
    if (present (obstype))  ru(n)% modtype  = obstype
    if (present (codetype)) ru(n)% codetype = codetype
    if (present (bf_type))  ru(n)% bf_type  = bf_type
    if (present (bf_subt))  ru(n)% bf_subt  = bf_subt
    if (present (db_kz  ))  ru(n)% db_kz (1:size(db_kz)) = db_kz
    if (present (center ))  ru(n)% center   = center
    if (present (stname ))  ru(n)% stname   = stname
    if (present (stname ))  ru(n)% l_stname = stname /= "" ! usable station id?
    if (present (satids)) then
       if (size (satids) > nsat) &
            call finish ('new_rule','max number of satids exhausted')
       ru(n)% satid(1:size(satids)) = satid (satids)
       if (any (ru(n)% satid == 0)) then
          l = minloc (abs (ru(n)% satid), dim=1)
          call finish ('new_rule','invalid satellite name: '//satids(l))
       end if
    end if
    if (present (varnos)) then
       if (size (varnos) > nv) &
         call finish ('new_rule','max number of varnos exhausted')
        do l = 1, nv
          if (varnos(l) == '') exit
          ru(n)% varno(l) = value_name (tab_varno, varnos(l))
          if (ru(n)% varno(l) == -999) &
          call finish ('new_rule','invalid varno: '//varnos(l))
        end do
    end if
    if (present (lat       )) ru(n)% lat        = lat
    if (present (lon       )) ru(n)% lon        = lon
    if (present (plim      )) ru(n)% plim       = plim
    if (present (zlim      )) ru(n)% zlim       = zlim
    if (present (mdzlim    )) ru(n)% mdzlim     = mdzlim
    if (present (sol_zenith)) ru(n)% sol_zenith = sol_zenith
    if (present (sat_zenith)) ru(n)% sat_zenith = sat_zenith
    if (present (bsurf     )) ru(n)% bsurf      = bsurf
    if (present (asurf     )) ru(n)% asurf      = asurf
    if (present (fr_land   )) ru(n)% fr_land    = fr_land
    if (present (fr_sea    )) ru(n)% fr_sea     = fr_sea
    if (present (tsurf     )) ru(n)% tsurf      = tsurf
    if (present (phase     )) ru(n)% phase      = phase
    if (present (instrs    )) then
      l = size(instrs)
      if (l > size(ru(n)% instrs)) call finish('new rule', 'too many instrs')
      ru(n)% instrs(1:l) = instrs(1:l)
    end if
    if (present (instrs_miss)) then
      l = size(instrs_miss)
      if (l > size(ru(n)% instrs_miss)) call finish('new rule', 'too many instrs_miss')
      ru(n)% instrs_miss(1:l) = instrs_miss(1:l)
    end if
    if (present (levs_miss)) then
      l = size(levs_miss)
      if (l > size(ru(n)% levs_miss)) call finish('new rule', 'too many levs_miss')
      ru(n)% levs_miss(1:l) = levs_miss(1:l)
    end if
    if (present (hours     )) then
       if (any (hours>=0)) ru(n)% lhour   = any (hours == ihh (ana_time))
    end if
    if (present (use    )) ru(n)% use     = use
    if (present (verb   )) ru(n)% verb    = verb
    if (present (msl    )) ru(n)% msl     = msl
    if (present (t      )) ru(n)% t       = t
    if (present (gp     )) ru(n)% gp      = gp
    if (present (uv     )) ru(n)% uv      = uv
    if (present (q      )) ru(n)% q       = q
    if (present (p      )) ru(n)% p       = p
    if (present (o      )) ru(n)% o       = o
    if (present (c      )) ru(n)% c(1:size(c)) = c

  end subroutine new_rule
!------------------------------------------------------------------------------
  subroutine get_rule (type, bf_type, bf_subt, db_kz, stname, obstype,       &
                       codetype, satid, center, varno, lat, lon, pobs, zobs, &
                       sol_zenith, sat_zenith, psurf, zsurf, fr_land, tsurf, &
                       phase, instr, instrs, levs, use, rule_num, verb, msl, &
                       t, gp, uv, q, p, o, c)
  integer         ,intent(in)            :: type
  integer         ,intent(in)            :: bf_type
  integer         ,intent(in)            :: bf_subt
  integer         ,intent(in)            :: db_kz
  character(len=*),intent(in)            :: stname
  integer         ,intent(in)            :: obstype
  integer         ,intent(in)            :: codetype
  integer         ,intent(in)  ,optional :: satid
  integer         ,intent(in)  ,optional :: center
  integer         ,intent(in)  ,optional :: varno
  real(wp)        ,intent(in)  ,optional :: lat
  real(wp)        ,intent(in)  ,optional :: lon
  real(wp)        ,intent(in)  ,optional :: pobs
  real(wp)        ,intent(in)  ,optional :: zobs
  real(wp)        ,intent(in)  ,optional :: sol_zenith
  real(wp)        ,intent(in)  ,optional :: sat_zenith
  real(wp)        ,intent(in)  ,optional :: psurf
  real(wp)        ,intent(in)  ,optional :: zsurf
  real(wp)        ,intent(in)  ,optional :: fr_land
  real(wp)        ,intent(in)  ,optional :: tsurf
  integer         ,intent(in)  ,optional :: phase
  integer         ,intent(in)  ,optional :: instr
  integer         ,intent(in)  ,optional :: instrs(:)
  type(t_ilev)    ,intent(in)  ,optional :: levs(:)
  integer         ,intent(out) ,optional :: use
  integer         ,intent(out) ,optional :: rule_num    ! rule no. setting 'use'
  integer         ,intent(out) ,optional :: verb
  integer         ,intent(out) ,optional :: msl
  type(t_set)     ,intent(out) ,optional :: t
  type(t_set)     ,intent(out) ,optional :: gp
  type(t_set)     ,intent(out) ,optional :: uv
  type(t_set)     ,intent(out) ,optional :: q
  type(t_set)     ,intent(out) ,optional :: p
  type(t_set)     ,intent(out) ,optional :: o
  type(t_set)     ,intent(out) ,optional :: c (:)

    logical  :: m    (n_ru)
    integer  :: n
    integer  :: i, j, k
    real(wp) :: llon (n_ru)
    logical  :: l_avail

    n = n_ru

    !-------------------------------------------
    ! Currently active rules (assimilation date)
    !-------------------------------------------
    m = ru(:n)% lhour

    !------------
    ! check types
    !------------
    where (ru(:n)% modtype  /= iud .and. ru(:n)% modtype  /= type)     m = .false.
    where (ru(:n)% bf_type  /= iud .and. ru(:n)% bf_type  /= bf_type)  m = .false.
    where (ru(:n)% bf_subt  /= iud .and. ru(:n)% bf_subt  /= bf_subt)  m = .false.
    where (ru(:n)% obstype  /= iud .and. ru(:n)% obstype  /= obstype)  m = .false.
    where (ru(:n)% codetype /= iud .and. ru(:n)% codetype /= codetype) m = .false.

    if (db_kz /= iud) then
l1:   do i=1,n
        if (.not.m(i) .or. ru(i)% db_kz(1) == iud)  cycle l1
l2:     do j = 1, nkz
          if (ru(i)% db_kz(j) == db_kz)             cycle l1
          if (ru(i)% db_kz(j) == iud)               exit  l2
        end do l2
        m(i) = .false.
      end do l1
    end if

    if (present (center)) then
      if (center /= iud .and. center /= -1) then
        where (ru(:n)% center /= iud .and. ru(:n)% center /= center) m = .false.
      end if
    else
      where   (ru(:n)% center /= iud)                                m = .false.
    end if

    if (present (satid)) then
       if (satid /= iud .and. satid /= -1) then
l3:       do i=1,n
             if (.not.m(i) .or. ru(i)% satid(1) == iud)  cycle l3
l4:          do j = 1, nsat
                if (ru(i)% satid(j) == satid)            cycle l3
                if (ru(i)% satid(j) == iud  )            exit  l4
             end do l4
             m(i) = .false.
          end do l3
       end if
    else
       !--------------------------------------------------------
       ! Special treatment of rules using satids instead of dbkz
       !--------------------------------------------------------
       where (ru(:n)% db_kz(1) == iud .and. &
              ru(:n)% satid(1) /= iud) m = .false.
    end if


    !------------
    ! check varno
    !------------
    if (present (varno)) then
      if (varno /= iud) then
l5:     do i=1,n
          if (.not.m(i) .or. ru(i)% varno(1) == iud)  cycle l5
l6:       do j=1,nv
            if (ru(i)% varno(j) == varno) cycle l5
            if (ru(i)% varno(j) == iud  ) exit  l6
          end do l6
          m(i) = .false.
        end do l5
      end if
    end if
    ! Moved from above:
    !where (ru(:n)% stname /= '' .and. ru(:n)% stname /= stname)  m = .false.
    !
    ! Optimization for NEC SX: mask slow string comparisons explicitly
    do i=1,n
      if (m(i) .and. ru(i)% l_stname) then
        if (ru(i)% stname /= stname) m(i) = .false.
      end if
    end do

    !---------------
    ! check location
    !---------------
    if (present(lon) .and. present(lat)) then
      llon = lon
      where (llon(:n) < ru(:n)% lon(1)) llon(:n) = llon(:n) + 360._wp
      where (ru(:n)% xlonlat)
        where (ru(:n)% lon(2) >= llon(:n) .and. &
               ru(:n)% lat(1) <= lat      .and. &
               ru(:n)% lat(2) >= lat) m = .false.
      endwhere
      where (ru(:n)% rlon)
        where (ru(:n)% lon(2) <  llon(:n)) m = .false.
      endwhere
      where (ru(:n)% rlat)
        where (ru(:n)% lat(1) > lat .or. ru(:n)% lat(2) < lat) m = .false.
      endwhere
    else
      where (ru(:n)% xlonlat .or. ru(:n)% rlon .or. ru(:n)% rlat) m = .false.
    endif

    if (present(pobs)) then
      where ((ru(:n)% plim(1) /= empty% plim(1) .and. &
              ru(:n)% plim(1) >         pobs   ) .or. &
             (ru(:n)% plim(2) /= empty% plim(2) .and. &
              ru(:n)% plim(2) <         pobs   )      ) m = .false.
    else
      where (ru(:n)% plim(1) /= empty% plim(1) &
        .or. ru(:n)% plim(2) /= empty% plim(2) ) m = .false.
    endif

    if (present(zobs)) then
      where (ru(:n)% zlim(1) /= empty% zlim(1) .and. &
             ru(:n)% zlim(2) /= empty% zlim(2) .and. &
             ru(:n)% zlim(2) < ru(:n)% zlim(1)       )
        where ((ru(:n)% zlim(1) >         zobs   ).and. &
               (ru(:n)% zlim(2) <         zobs   )      ) m = .false.
      elsewhere
        where ((ru(:n)% zlim(1) /= empty% zlim(1) .and. &
                ru(:n)% zlim(1) >         zobs   ).or.  &
               (ru(:n)% zlim(2) /= empty% zlim(2) .and. &
                ru(:n)% zlim(2) <         zobs   )      ) m = .false.
      end where
    else
      where (ru(:n)% zlim(1) /= empty% zlim(1) &
        .or. ru(:n)% zlim(2) /= empty% zlim(2) ) m = .false.
    endif

    if (present(zsurf)) then
      where ((ru(:n)% mdzlim(1) /= empty% mdzlim(1) .and. &
              ru(:n)% mdzlim(1) >         zsurf   ) .or.  &
             (ru(:n)% mdzlim(2) /= empty% mdzlim(2) .and. &
              ru(:n)% mdzlim(2) <         zsurf   )      ) m = .false.
    else
      where (ru(:n)% mdzlim(1) /= empty% mdzlim(1) &
        .or. ru(:n)% mdzlim(2) /= empty% mdzlim(2) ) m = .false.
    endif

    if (present(sol_zenith)) then
      where ((ru(:n)% sol_zenith(1) /= empty% sol_zenith(1) .and. &
              ru(:n)% sol_zenith(1) >         sol_zenith   ) .or. &
             (ru(:n)% sol_zenith(2) /= empty% sol_zenith(2) .and. &
              ru(:n)% sol_zenith(2) <         sol_zenith   )      ) m = .false.
    else
      where (ru(:n)% sol_zenith(1) /= empty% sol_zenith(1) &
        .or. ru(:n)% sol_zenith(2) /= empty% sol_zenith(2) ) m = .false.
    endif

    if (present(sat_zenith)) then
      where ((ru(:n)% sat_zenith(1) /= empty% sat_zenith(1) .and. &
              ru(:n)% sat_zenith(1) >         sat_zenith   ) .or. &
             (ru(:n)% sat_zenith(2) /= empty% sat_zenith(2) .and. &
              ru(:n)% sat_zenith(2) <         sat_zenith   )      ) m = .false.
    else
      where (ru(:n)% sat_zenith(1) /= empty% sat_zenith(1) &
        .or. ru(:n)% sat_zenith(2) /= empty% sat_zenith(2) ) m = .false.
    endif

    if (present(psurf).and.present(pobs)) then
      where (empty% bsurf /= ru(:n)% bsurf .and. &
             pobs-psurf    > ru(:n)% bsurf) m = .false.
      where (empty% asurf /= ru(:n)% asurf .and. &
             psurf-pobs    > ru(:n)% asurf) m = .false.
    else
      where (ru(:n)% asurf /= empty% asurf &
        .or. ru(:n)% bsurf /= empty% bsurf) m = .false.
    endif

    if (present(phase)) then
      where ((ru(:n)% phase(1) /= empty% phase(1) .and. &
              ru(:n)% phase(1) >         phase   ) .or. &
             (ru(:n)% phase(2) /= empty% phase(2) .and. &
              ru(:n)% phase(2) <         phase   )      ) m = .false.
    else
      where (ru(:n)% phase(1) /= empty% phase(1) &
        .or. ru(:n)% phase(2) /= empty% phase(2) ) m = .false.
    endif

    if (present(instr)) then
      if (instr /= iud) then
        do i=1,n
          if (m(i) .and. ru(i)% instrs(1) /= iud) then
            if (.not.any(ru(i)% instrs(1:ninstr) == instr)) m(i) = .false.
          end if
        end do
      end if
    end if

    if (present(instrs)) then
      do i=1,n
        do j=1,size (ru(i)% instrs_miss)
          if (ru(i)% instrs_miss(j) /= iud) then
            if (any(instrs(:) == ru(i)% instrs_miss(j))) m(i) = .false.
          end if
        end do
      end do
    end if

    if (present(levs)) then
      do i=1,n
        ! The rule is to be applied if some level in levs_miss is missing
        if (any(ru(i)% levs_miss(:)%value > 0 .and. ru(i)% levs_miss(:)%instr > 0)) then
          l_avail = .true.
          lev_loop: do j=1,nlev
            if (ru(i)% levs_miss(j)%value > 0 .and. ru(i)% levs_miss(j)%instr > 0) then
              do k = 1, size(levs, 1)
                if (levs(k) == ru(i)% levs_miss(j)) then
                  ! level is available
                  cycle lev_loop
                end if
              end do
              ! level is missing
              l_avail = .false.
              exit lev_loop
            end if
          end do lev_loop
          if (l_avail) m(i) = .false. ! All levels in levs_miss available -> do not apply rule
        end if
      end do
    end if

    if (present(fr_land)) then
      where ((ru(:n)% fr_land /= empty% fr_land .and. &
              ru(:n)% fr_land >         fr_land) .or. &
             (ru(:n)% fr_sea  /= empty% fr_sea  .and. &
              ru(:n)% fr_sea  >   1._wp-fr_land)      ) m = .false.
    else
      where (ru(:n)% fr_land /= empty% fr_land &
        .or. ru(:n)% fr_sea  /= empty% fr_sea  ) m = .false.
    endif

    if (present(tsurf)) then
      where ((ru(:n)% tsurf(1) /= empty% tsurf(1) .and. &
              ru(:n)% tsurf(1) >  -t0c + tsurf   ) .or. &
             (ru(:n)% tsurf(2) /= empty% tsurf(2) .and. &
              ru(:n)% tsurf(2) <  -t0c + tsurf   )      ) m = .false.
    else
      where (ru(:n)% tsurf(1) /= empty% tsurf(1) &
        .or. ru(:n)% tsurf(2) /= empty% tsurf(2) ) m = .false.
    endif

    !--------------------------------
    ! rule 1 should never be disabled
    !--------------------------------
    if (n_ru == 0  ) call finish('get_rule','n_ru == 0')
    if (.not.(m(1))) call finish('get_rule','rule 1 is disabled')
    m(1) = .true.

    !------------------
    ! return parameters
    !------------------
    if (present(use)) then
      use = iud
      if (present (rule_num)) rule_num = 0
      do i=n,1,-1
        if (m(i) .and. ru(i)% use  /= iud) then
          use  = ru(i)% use
          if (present (rule_num)) rule_num = i
          exit
        end if
      end do
    end if

    if (present(verb)) then
      verb = iud
      do i=n,1,-1
        if (m(i) .and. ru(i)% verb /= iud) then
          verb = ru(i)% verb
          exit
        end if
      end do
    end if

    if (present(msl)) then
      msl = iud
      do i=n,1,-1
        if (m(i) .and. ru(i)% msl  /= iud) then
          msl  = ru(i)% msl
          exit
        end if
      end do
    end if

    if (present(t)) then
      do i=n,1,-1
        if (.not.m(i))         cycle
        call set_set (ru(i)% t, t)
      end do
    end if

!print *,'get_rule',m(:n)
    if (present(gp)) then
      do i=n,1,-1
!print *,'get_rule',i,m(i),ru(i)% gp, gp
        if (.not.m(i))         cycle
        call set_set (ru(i)% gp, gp)
!print *,'get_rule',i,m(i),ru(i)% gp, gp
      end do
    end if

    if (present(uv)) then
      do i=n,1,-1
        if (.not.m(i))         cycle
        call set_set (ru(i)% uv, uv)
      end do
    end if

    if (present(q)) then
      do i=n,1,-1
        if (.not.m(i))         cycle
        call set_set (ru(i)% q, q)
      end do
    end if

    if (present(p)) then
      do i=n,1,-1
        if (.not.m(i))         cycle
        call set_set (ru(i)% p, p)
      end do
    end if

    if (present(o)) then
      do i=n,1,-1
        if (.not.m(i))         cycle
        call set_set (ru(i)% o, o)
      end do
    end if

    if (present(c)) then
      do i=n,1,-1
        if (.not.m(i))         cycle
        call set_set_v (ru(i)% c, c)
      end do
    end if

  contains

    subroutine set_set (src, dst)
    type(t_set), intent(in)    :: src
    type(t_set), intent(inout) :: dst

    if   (src% use      /=iud .and. dst% use      == iud) dst% use      = src% use
    if   (src% frm_vq   /=iud .and. dst% frm_vq   == iud) dst% frm_vq   = src% frm_vq
    where(src% sgm_oc   /=rud .and. dst% sgm_oc   == rud) dst% sgm_oc   = src% sgm_oc
    where(src% sgm_fg   /=rud .and. dst% sgm_fg   == rud) dst% sgm_fg   = src% sgm_fg
    where(src% bnd_fg   /=rud .and. dst% bnd_fg   == rud) dst% bnd_fg   = src% bnd_fg
    if   (src% bnd_obs  /=rud .and. dst% bnd_obs  == rud) dst% bnd_obs  = src% bnd_obs
    if   (src% sgm_vq   /=rud .and. dst% sgm_vq   == rud) dst% sgm_vq   = src% sgm_vq
    if   (src% m_rej    /=iu2 .and. dst% m_rej    == iu2) dst% m_rej    = src% m_rej
    if   (src% prc      /=iu1 .and. dst% prc      == iu1) dst% prc      = src% prc
    if   (src% v_loc    /=0.  .and. dst% v_loc    == 0. ) dst% v_loc    = src% v_loc
    if   (src% h_loc    /=0.  .and. dst% h_loc    == 0. ) dst% h_loc    = src% h_loc
    if   (src% ekf_pass /=0   .and. dst% ekf_pass == 0  ) dst% ekf_pass = src% ekf_pass
    end subroutine set_set

    subroutine set_set_v (src, dst)
    type(t_set), intent(in)    :: src(nc)
    type(t_set), intent(inout) :: dst(nc)

    where (src% use   /= iud .and. dst% use   == iud) dst% use    = src% use
    where (src% frm_vq/= iud .and. dst% frm_vq== iud) dst% frm_vq = src% frm_vq
    where (src% sgm_vq/= rud .and. dst% sgm_vq== rud) dst% sgm_vq = src% sgm_vq
    where (src% m_rej /= iu2 .and. dst% m_rej == iu2) dst% m_rej  = src% m_rej
    where (src% prc   /= iu1 .and. dst% prc   == iu1) dst% prc    = src% prc

    where(src%sgm_oc(1)/=rud.and.dst%sgm_oc(1)==rud)dst%sgm_oc(1)=src%sgm_oc(1)
    where(src%sgm_oc(2)/=rud.and.dst%sgm_oc(2)==rud)dst%sgm_oc(2)=src%sgm_oc(2)
    where(src%sgm_oc(3)/=rud.and.dst%sgm_oc(3)==rud)dst%sgm_oc(3)=src%sgm_oc(3)
    where(src%sgm_oc(4)/=rud.and.dst%sgm_oc(4)==rud)dst%sgm_oc(4)=src%sgm_oc(4)

    where(src%sgm_fg(1)/=rud.and.dst%sgm_fg(1)==rud)dst%sgm_fg(1)=src%sgm_fg(1)
    where(src%sgm_fg(2)/=rud.and.dst%sgm_fg(2)==rud)dst%sgm_fg(2)=src%sgm_fg(2)
    where(src%sgm_fg(3)/=rud.and.dst%sgm_fg(3)==rud)dst%sgm_fg(3)=src%sgm_fg(3)

    where(src%bnd_fg(1)/=rud.and.dst%bnd_fg(1)==rud)dst%bnd_fg(1)=src%bnd_fg(1)
    where(src%bnd_fg(2)/=rud.and.dst%bnd_fg(2)==rud)dst%bnd_fg(2)=src%bnd_fg(2)

    where(src%bnd_obs  /=rud.and.dst%bnd_obs  ==rud)dst%bnd_obs  =src%bnd_obs

    where(src%v_loc    /=0. .and.dst%v_loc    ==0. )dst%v_loc    =src%v_loc
    where(src%h_loc    /=0. .and.dst%h_loc    ==0. )dst%h_loc    =src%h_loc
    where(src%ekf_pass /=0  .and.dst%ekf_pass ==0  )dst%ekf_pass =src%ekf_pass
    end subroutine set_set_v

  end subroutine get_rule
!------------------------------------------------------------------------------
  subroutine read_nml_rules
  !-----------------------------
  ! read namelist groups /rules/
  ! printout of settings
  !-----------------------------

    !--------------------------------------------
    ! namelist variables declaration
    ! corresponding to components of type t_rules
    !--------------------------------------------
    character(len=64) :: comment
    integer           :: type
    integer           :: obstype
    integer           :: codetype
    integer           :: bf_type
    integer           :: bf_subt
    integer           :: center
    integer           :: db_kz (nkz)
    character(len=10) :: stname
    character(len=8)  :: satids(nsat)   ! satellite ids (mnemonics)
    character(len=12) :: varnos(nv)

    real(wp)          :: lat        (2)
    real(wp)          :: lon        (2)
    logical           :: xlonlat        ! exclude area
    real(wp)          :: plim       (2)
    real(wp)          :: zlim       (2)
    real(wp)          :: mdzlim     (2)
    real(wp)          :: sol_zenith (2)
    real(wp)          :: sat_zenith (2)
    real(wp)          :: bsurf
    real(wp)          :: asurf
    real(wp)          :: fr_land
    real(wp)          :: fr_sea
    real(wp)          :: tsurf      (2)
    integer           :: phase      (2)
    integer           :: instrs     (ninstr)
    integer           :: instrs_miss(ninstm)
    type(t_ilev)      :: levs_miss  (nlev)
    integer           :: hours      (8) ! analysis hours

    character(len=8)  :: use
    integer           :: verb
    integer           :: msl
    type(t_set)       :: t
    type(t_set)       :: gp
    type(t_set)       :: uv
    type(t_set)       :: q
    type(t_set)       :: p
    type(t_set)       :: o
    type(t_set)       :: c (nc)

    namelist /RULES/ comment,                                                  &
                     type, obstype, codetype, bf_type, bf_subt, db_kz, varnos, &
                     satids, stname, lat, lon, xlonlat, plim, zlim, mdzlim,    &
                     sol_zenith, sat_zenith, bsurf, asurf, fr_land, fr_sea,    &
                     tsurf, phase, instrs, instrs_miss, levs_miss, hours, use, &
                     center, verb, msl, t, gp, uv, q, p, o, c
    !----------------------
    ! other local variables
    !----------------------
    logical :: first
    integer :: ierr
    integer :: n
    integer :: l
#if defined(__ibm__)
    integer :: ios
#endif
    !---------------------------------
    ! repeatedly read namelist /RULES/
    !---------------------------------
    call init_fdbk_tables ()
    first = .true.
    do
      !----------------------------------------
      ! preset namelist variables with defaults
      !----------------------------------------
      comment     = 'added by namelist /RULES/'
      type        = empty% modtype
      obstype     = empty% obstype
      codetype    = empty% codetype
      bf_type     = empty% bf_type
      bf_subt     = empty% bf_subt
      center      = empty% center
      db_kz       = empty% db_kz
      stname      = empty% stname
      lat         = empty% lat
      lon         = empty% lon
      plim        = empty% plim
      zlim        = empty% zlim
      mdzlim      = empty% mdzlim
      sol_zenith  = empty% sol_zenith
      sat_zenith  = empty% sat_zenith
      bsurf       = empty% bsurf
      asurf       = empty% asurf
      fr_land     = empty% fr_land
      fr_sea      = empty% fr_sea
      tsurf       = empty% tsurf
      phase       = empty% phase
      instrs      = empty% instrs
      instrs_miss = empty% instrs_miss
      levs_miss   = empty% levs_miss
      xlonlat     = empty% xlonlat
      satids      = ''
      varnos      = ''
      hours       = -1
      use         = ''
      verb        = empty% verb
      msl         = empty% msl
      t           = empty% t
      gp          = empty% gp
      uv          = empty% uv
      q           = empty% q
      p           = empty% p
      o           = empty% o
      c           = empty% c
      !--------------
      ! read namelist
      !--------------
      if (dace% lpio) then
        call position_nml ('RULES' ,lrewind=first ,status=ierr)
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=RULES, iostat=ios)
          if (ios/=0) call finish ('nml_run_flags','ERROR in namelist /RULES/')
#else
          read (nnml ,nml=RULES)
#endif
        end select
      endif
      first = .false.
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit
      !---------------------
      ! checks, unit changes
      !---------------------
      if    (bsurf /= empty% bsurf) bsurf = bsurf * 100._wp ! hPa -> Pa
      if    (asurf /= empty% asurf) asurf = asurf * 100._wp ! hPa -> Pa
      where (plim  /= empty% plim ) plim  = plim  * 100._wp
!      !----------------------------------
!      ! restrict ..%sgm_vq to values < 30
!      ! (prevent log(0) in VQC-routines)
!      !----------------------------------
!      t % sgm_vq = min (t% sgm_vq, 30._wp)
!      gp% sgm_vq = min (t% sgm_vq, 30._wp)
!      uv% sgm_vq = min (t% sgm_vq, 30._wp)
!      q % sgm_vq = min (t% sgm_vq, 30._wp)
!      p % sgm_vq = min (t% sgm_vq, 30._wp)
!      o % sgm_vq = min (t% sgm_vq, 30._wp)
!      c % sgm_vq = min (t% sgm_vq, 30._wp)
      !--------------------------------------
      ! store rule in structure and broadcast
      !--------------------------------------
      n_ru = n_ru + 1
      if (n_ru > size (ru)) &
        call finish ('read_nml_rules','max number of rules exhausted')
      n = n_ru
      if (dace% lpio) then
        if (lon(2) < lon(1)) call finish ('read_nml_rules','lon(2) < lon(1)')
        if (lat(2) < lat(1)) call finish ('read_nml_rules','lat(2) < lat(1)')
        ru(n)% comment  = comment
        ru(n)% modtype  = type
        ru(n)% obstype  = obstype
        ru(n)% codetype = codetype
        ru(n)% bf_type  = bf_type
        ru(n)% bf_subt  = bf_subt
        ru(n)% center   = center
        ru(n)% db_kz    = db_kz
        ru(n)% stname   = stname
        ru(n)% l_stname = stname /= ""
        ru(n)% satid    = satid (satids)
        if (any (ru(n)% satid == 0)) then
           l = minloc (abs (ru(n)% satid), dim=1)
           call finish ('read_nml_rules','invalid satellite name: '//satids(l))
        end if
#ifdef _CRAYFTN
!DIR$ NEXTSCALAR                        ! work around bug in crayftn <= 8.2.2
#endif
        do l = nsat, 1, -1
           if (ru(n)% satid(l) == -1) ru(n)% satid(l) = iud
        end do
        do l = 1, nv
           if (varnos(l) == '') exit
           ru(n)% varno(l) = value_name (tab_varno, varnos(l))
           if (ru(n)% varno(l) == -999) then
             call finish('read_nml_rules', 'Invalid varno: '//varnos(l))
           end if
        end do
        ru(n)% lat         = lat
        ru(n)% lon         = lon
        ru(n)% xlonlat     = xlonlat
        ru(n)% rlon        = .not. xlonlat .and. all (lon/= empty% lon)
        ru(n)% rlat        = .not. xlonlat .and. all (lat/= empty% lat)
        ru(n)% lhour       = .true.
        if (any (hours >= 0)) ru(n)% lhour = any (hours  == ihh (ana_time))
        ru(n)% plim        = plim
        ru(n)% zlim        = zlim
        ru(n)% mdzlim      = mdzlim
        ru(n)% sol_zenith  = sol_zenith
        ru(n)% sat_zenith  = sat_zenith
        ru(n)% bsurf       = bsurf
        ru(n)% asurf       = asurf
        ru(n)% fr_land     = fr_land
        ru(n)% fr_sea      = fr_sea
        ru(n)% tsurf       = tsurf
        ru(n)% phase       = phase
        ru(n)% instrs      = instrs
        ru(n)% instrs_miss = instrs_miss
        ru(n)% levs_miss   = levs_miss
        ru(n)% verb        = verb
        ru(n)% msl         = msl
        ru(n)% t           = t
        ru(n)% gp          = gp
        ru(n)% uv          = uv
        ru(n)% q           = q
        ru(n)% p           = p
        ru(n)% o           = o
        ru(n)% c           = c
        if (use /= '') then
          ru(n)% use   = stat_key (use)
          if (ru(n)% use < 1) call finish ('read_nml_rules, namelist /RULES/',&
            'invalid value for parameter "use" : '//use)
        else
          ru(n)% use   = empty% use
        endif
      endif
      l = size (transfer (empty ,(/' '/)))
      call p_bcast_derivedtype (ru(n), l, (dace% pio), (dace% comm))
    end do
    !---------------------
    ! printout of settings
    !---------------------
    call print_rules

  end subroutine read_nml_rules
!==============================================================================
  elemental function ilev_eq (ilev1, ilev2) result (y)
    logical                    :: y
    type (t_ilev) ,intent(in) :: ilev1
    type (t_ilev) ,intent(in) :: ilev2
    y = ilev1% value == ilev2% value .and. ilev1% instr == ilev2% instr
  end function ilev_eq
!------------------------------------------------------------------------------
end module mo_obs_rules
