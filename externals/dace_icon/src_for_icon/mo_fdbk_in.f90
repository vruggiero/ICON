!
!+ Read feedback files (new format, common with COSMO and LETKF)
!
MODULE mo_fdbk_in
!
! Description:
!   Read feedback files (new format, common with COSMO)
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_4         2009/03/26 Andreas Rhodin
!  New module to read feedback files in 3DVAR/LETKF
! V1_5         2009/05/25 Andreas Rhodin
!  finish routine to read feedback files
! V1_7         2009/08/24 Harald Anlauf
!  fix bookkeeping of observations, update observation statistics
! V1_8         2009/12/09 Andreas Rhodin
!  read VQC weight from 'veri_data'
! V1_9         2010/04/20 Andreas Rhodin
!  changed spot% phase from i1 to i2 to be able to hold GPSRO PCD flags
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  new namelist /FOF_INPUT/
! V1_19        2012-04-16 Andreas Rhodin
!  changes for COSMO fof files and RADAR operator
! V1_20        2012-06-18 Andreas Rhodin
!  correctly merge status flags from feedback file and 3dvar/LETKF
! V1_22        2013-02-13 Andreas Rhodin
!  changes for RADAR observation operator
! V1_27        2013-11-08 Andreas Rhodin
!  bugfix for aircraft input from feedback file: correctly set station pressure
! V1_28        2014/02/26 Andreas Rhodin
!  implement 'max_proc' for input from fof-files
! V1_29        2014/04/02 Andreas Rhodin
!  obs error modification for RREFL implemented
!  fix behaviour of LETKF namelist flag rm_fg_check
! V1_31        2014-08-21 Andreas Rhodin
!  read_feedback: add check for corrupt fof-files
!  diagnose empty feedback-(fof-)files and skip reading them
! V1_35        2014-11-07 Andreas Rhodin
!  /FOF_INPUT/, new: percent_fg_check (% of hits required for rejection)
! V1_37        2014-12-23 Harald Anlauf
!  Change ijdp, index_x from i2 to int
! V1_42        2015-06-08 Harald Anlauf
!  read/restore relevant data for GPSRO from feedback file
! V1_44        2015-09-30 Harald Anlauf
!  move type t_bufr_inv to mo_t_obs to simplify dependencies
! V1_48        2016-10-06 Andreas Rhodin
!  read_feedback: set body% ps (station height) correctly (single level data)
! V1_50        2017-01-09 Andreas Rhodin
!  new parameter 'check_member' in namelist /fof_input/
! V1_51        2017-02-24 Andreas Rhodin
!  read_feedback: set 'satid' for observation type SOIL (for COMET)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

!-------------
! modules used
!-------------
use mo_kind,        only: wp, sp, i2        ! kind parameters
use mo_exception,   only: finish            ! abort routine
use mo_dace_string, only: char5             ! convert integer -> char(len=5)
use mo_mpi_dace,    only: dace,            &! MPI group info
                          p_bcast,         &! MPI broadcast routine
                          p_or, p_max       ! MPI or/max routine
use utilities,      only: uvrot2uv          ! convert u,v from the rotated system
use mo_run_params,  only: path_file         ! concatenate path/filename
use mo_namelist,    only: position_nml,    &! routine to position nml group
                          nnml,            &! namelist fortran unit number
                          POSITIONED        ! position_nml: OK return flag
use mo_physics,     only: gacc              ! gravity acceleration
use mo_time,        only: t_time,          &! time derived type
                          init_time,       &! set time derived type variable
                          operator(+)       ! add times
use mo_t_obs,       only: source,          &! table of observation files
                          n_source,        &!    no of Report source files to read
                          bufr_inv,        &! observation file inventory
                          FT_FEEDBACK,     &! flag for feedback file
                          FT_SATPP,        &! flag for satpp file
                          FT_EMPTY,        &! flag for empty file
                          new_spot,        &! allocate mem for observation header
                          new_obs,         &! reserve memory for new observations
                          t_obs,           &! observation data type
                          t_spot,          &! component of t_obs
                          t_head,          &! component of t_spot
                          set_xuv,         &! set components of t_spot
                          frac_mdlsfc,     &! land/ice/snow fraction from flag
                          invalid           ! default for invalid value
use mo_obs_tables,  only: idb_dbk,         &! index in table rept_stat
                          check_report_0,  &! init. flags, standard checks
                          check_report_1,  &! standard checks
                          rept_char,       &! report type characteristics table
                          rept_use          ! observation usage flag table
use mo_obstypes,    only: t_obsid,         &! derived type for obstype,etc.
                          obstype_code      ! derive obstype,etc. from codetype
use mo_t_use,       only: reverse_code,    &! convert code back (feedback to 3dvar)
                          reverse_bits,    &! convert bit-field back
                          CHK_CONSIST,     &! consistency check flag
                          checks => chk,   &! table of 'checks'
                          stats             ! table of 'states'
use mo_t_datum,     only: t_datum,         &! data type for a single entity
                          inv_datum         ! invalid value settings
use mo_fdbk,        only: t_fdbk,          &! feedback file data type
                          setup_fdbk,      &! set up feedback-file data structure
                          open_fdbk_read  ,&! open feedback file for read access
                          read_meta       ,&!
                          get_veri_index  ,&! return indices of verification runs
                          get_veri        ,&! read verification entry
                          print_fdbk,      &! print feedback file
                          close_fdbk,      &! close feedback file
                          cleanup_fdbk      ! deallocate components
use mo_t_table,     only: name_value        ! get name from value
use mo_fdbk_tables, only: n_ot,            &! number of observation types
                          VE_VQC_WEIGHT,   &! VQC weight index
                          VE_DETERM,       &! deterministic run flag
                          VT_FIRSTGUESS,   &! first guess flag
                          VN_P,            &! pressure variable indicator
                          VN_PS,           &! surface pressure  indicator
                          VN_U, VN_V,      &! wind variable indicator
                          VN_U10M,VN_V10M, &! wind variable indicator
                          VN_T,  VN_T2M,   &! temperature variable indicator
                          VN_RH, VN_RH2M,  &! rel. hum (2m) indicator
                          VN_RADVEL,       &! radial velocity
                          VN_RREFL,        &! radar reflectivity
                          FL_DATASET,      &! check flag value:  dataset
                          FL_FG,           &!                    first guess
                          FL_NONE,         &!                    none
                          ST_ACTIVE,       &! status flag value: active
                          ST_PASSIVE,      &!                    passive
                          ST_REJECTED,     &!                    rejected
                          ST_PAS_REJ,      &!                    passive rej.
                          ST_OBS_ONLY,     &!                 no model equiv.
               flags_t => flags,           &! check/flags table
              status_t => status            ! status      table
use mo_netcdf_param,only: NF_NOERR        ,&!
                          NF_FILL_REAL      !
use mo_t_netcdf,    only: ncid,            &! NetCDF file id
                          stanc,           &! NetCDF start  parameter
                          counc,           &! NetCDF count  parameter
                          strnc,           &! NetCDF stride parameter
                          get_var           ! read variable
use mo_fdbk_tables, only: OT_TEMP,         &! TEMP              obstype value
                          OT_SATOB,        &! AMV               obstype value
                          OT_AIREP,        &! Aircraft          obstype value
                          OT_GPSRO,        &! radio occultation obstype value
                          OT_GPSGB,        &! GNSS ground based obstype value
                          OT_RAD,          &! radiances         obstype value
                          OT_SOIL,         &! Soil              obstype value
                          OT_RADAR          ! Radar             obstype value
use mo_wmo_tables,  only: WMO0_ECMWF,      &! generating center: ECMWF
                          WMO0_DWD          !                    DWD
use mo_rad,         only: rad_set, n_set    ! RADIANCES options
use mo_tovs,        only: read_fdbk_rad     ! read RADIANCES specific table
use mo_radar,       only: read_fdbk_radar   ! read RADAR     specific table
use mo_occ,         only: read_fdbk_occ     ! restore GPSRO  specific tables
use mo_std,         only: read_fdbk_gpsgb   ! restore GPSGB  specific tables
use mo_temp,        only: read_fdbk_temp    ! fix     TEMP   level count
implicit none

!----------------
! public entities
!----------------
private
public :: read_feedback       ! read the feedback file
public :: read_nml_fof_input  ! read namelist /FOF_INPUT/
public :: percent_fg_check    ! percentage of fg checks required for rejection
public :: check_member        ! flag to check for correct member in fof-file

!---------------------
! namelist /FOF_INPUT/
!---------------------
logical  :: flag_dataset    = .false. ! set FL_DATASET for rejected data
logical  :: rm_fg_check     = .false. ! remove first guess check
logical  :: force_active    = .false. ! force PASSIVE observations to ACTIVE
logical  :: check_member    = .false. ! check for correct member in fof-file
logical  :: replace_o_fg    = .false. ! replace observations by first guess
real(sp) :: percent_fg_check= 0.      ! percentage of fg checks for rejection
real(sp) :: ps_obs_error    = 0.      ! <0: factor on, >0: value of PS obs. error
real(sp) :: uv_obs_error    = 0.      ! <0: factor on, >0: value of wind obs.error
real(sp) :: uv10m_obs_error = 0.      ! <0: factor on, >0: value of wind obs.error
real(sp) :: radvel_obs_error= 0.      ! <0: factor on, >0: value of radial wind oe
real(sp) :: refl_obs_error  = 0.      ! <0: factor on, >0: value of radar reflectivity oe
real(sp) :: t_obs_error     = 0.      ! <0: factor on, >0: value of temp.obs.error
real(sp) :: t2m_obs_error   = 0.      ! <0: factor on, >0: value of temp.obs.error
real(sp) :: rh_obs_error    = 0.      ! <0: factor on, >0: value of hum. obs.error
real(sp) :: rh2m_obs_error  = 0.      ! <0: factor on, >0: value of hum. obs.error


namelist /fof_input/ flag_dataset, ps_obs_error, uv_obs_error, t_obs_error,   &
                     rh_obs_error, rh2m_obs_error, rm_fg_check, force_active, &
                     radvel_obs_error, refl_obs_error, percent_fg_check,      &
                     check_member, replace_o_fg, t2m_obs_error, uv10m_obs_error

!==============================================================================
contains
!==============================================================================
  subroutine read_feedback (pass, obs, qc, bg, id, chk, src, spec)
  integer     ,intent(in)              :: pass
  type(t_obs) ,intent(inout) ,optional :: obs(:)
  logical     ,intent(in)    ,optional :: qc   ! read status flags
  logical     ,intent(in)    ,optional :: bg   ! read background values
  logical     ,intent(in)    ,optional :: id   ! read observation id
  logical     ,intent(in)    ,optional :: chk  ! apply basic checks, rules
  integer     ,intent(in)    ,optional :: src  ! source for id = .false.
  logical     ,intent(in)    ,optional :: spec ! read specific data
  !-------------------------------------------------------
  ! read observations from feedback file
  ! pass = 1 : scan the file (called from set_input_boxes)
  ! pass = 2 : read the file, store in 'obs'
  !-------------------------------------------------------
    target :: obs

    !----------------
    ! local variables
    !----------------
    logical               :: lexicdf   ! existence of feedback input file
    integer               :: ifile     ! file index
!   integer               :: entry     ! position in source file (subset)
    character(len=192)    :: filename  ! file name
    type (t_fdbk)         :: fb        ! feedback file content
    type (t_time)         :: ref_time  ! reference date/time
    type (t_time)         :: minutes   ! observ.time - reference time
    logical               :: lqc
    logical               :: lbg
    logical               :: lid
    logical               :: lchk
    logical               :: lspec
    logical               :: luvbug
!   integer               :: obstype   ! observation type
    integer               :: ib        ! observation box index
    type(t_obs) ,pointer  :: bi

    integer  ,allocatable :: itmp    (:)
    integer  ,allocatable :: i_body  (:)
    integer  ,allocatable :: l_body  (:)
    integer  ,allocatable :: i_spec  (:)
    integer  ,allocatable :: l_spec  (:)
    integer  ,allocatable :: state(:), check(:), flags(:)
    real(wp) ,allocatable :: rtmp    (:)
    integer               :: ot
    integer               :: pe
    integer               :: n         ! number of header entries
    integer               :: is, ie    ! header entry index range
    integer               :: nb        ! number of body entries
    integer               :: isb, ieb  ! body entry index range
    type(t_spot) ,pointer :: s (:)
    type(t_head) ,pointer :: hd(:)
    type(t_datum),pointer :: bd(:)
    integer               :: ierr
    integer               :: i
    integer               :: j1, jn
    integer               :: i0, in, ni, j
    integer               :: idx (1)   ! index in verification data
    integer               :: nidx      ! number of indices
    type (t_obsid)        :: obsinfo
    real(wp)              :: u, v
    ! radiance specific stuff:
    real(sp), allocatable :: l2c(:), emiss(:)
    integer,  allocatable :: nwc_flag(:), tovs_flag(:)

    !---------------------------
    ! process optional arguments
    !---------------------------
    lqc   = .true.;  if (present (qc))   lqc   = qc
    lbg   = .true.;  if (present (bg))   lbg   = bg
    lid   = .true.;  if (present (id))   lid   = id
    lchk  = .false.; if (present (chk))  lchk  = chk
    lspec = .false.; if (present (spec)) lspec = spec
    !------------------------
    ! loop over files to read
    !------------------------
    do ifile = 1, n_source
!     entry = 0

      !----------------------------
      ! skip if not a feedback-file
      !----------------------------
      if (source(ifile)% filetype /= FT_FEEDBACK) cycle
      !---------
      ! chose PE
      !---------
      select case (pass)
      case (1)
        pe = dace% pio
      case (2)
        pe = source(ifile)% pe
      case default
        call finish('read_feedback','invalid pass')
      end select
      if (dace% pe == pe) then

        !------------------------------
        ! derive full pathname/filename
        !------------------------------
        filename = path_file (source(ifile)% path, source(ifile)% file)

        !---------------------------------
        ! check existence of feedback file
        !---------------------------------
        inquire (file=trim(filename), exist=lexicdf)
        if (.not. lexicdf) &
          call finish('read_feedback',trim(filename) // ' not existent')

        !-------------------
        ! open feedback file
        !-------------------
        call setup_fdbk     (fb% nc)
        call open_fdbk_read (fb, filename)
        if (fb% nc% error /= NF_NOERR) &
          call finish('read_feedback','cannot open '//trim(filename))

        !------------------------------------------------
        ! read global attributes + verification meta data
        !------------------------------------------------
        call read_meta      (fb, filename)
        if (pass==1) call print_fdbk (fb)
        if (fb% nc% error /= NF_NOERR) &
          call finish('read_feedback','read_global_attributes '//trim(filename))
        if (fb% reftime > 2400) then
           write(0,*) 'read_feedback: WARNING: verification_ref_time', &
                fb% reftime, '> 2400 found in file ', trim (filename)
        end if
        call init_time (ref_time, yyyymmdd = fb% refdate,    &
                                  hhmmss   = fb% reftime *100)
        !-------------------------------------
        ! check for bug in versions up to 1.01
        !-------------------------------------
        luvbug = ( fb% source      (1:5) =='COSMO' .and. &
                   fb% institution (1:5) /='CNMCA' .and. &
                   fb% institution (1:5) /='COMET' .and. &
                  (fb% version           ==' 1.01' .or.  &
                   fb% version           ==' 1.00').and. &
                  (fb% pole(1)           /=  0._sp .or.  &
                   fb% pole(2)           /=  0._sp)      )
        !---------------------
        ! pass = 1 : scan file
        !---------------------
        if (pass == 1) then
          !--------------------------------
          ! scan file for observation types
          !--------------------------------
          allocate (itmp (fb% n_hdr))
          ncid  = fb%nc% ncid
          stanc = 1
          counc = fb% n_hdr
          strnc = 1
          call get_var (itmp, 'obstype')
          !------------------------------
          ! store file inventory (pass 1)
          !------------------------------
          if (ifile>1) bufr_inv%subseto(ifile) = bufr_inv%subseto(ifile-1)
          do ot = 1, n_ot
            n = count (itmp == ot)
            if (n>0) then
              bufr_inv(ot)% file   (ifile) = .true.
              bufr_inv(ot)% nrec           = bufr_inv(ot)% nrec           + n
              bufr_inv(ot)% nsubset        = bufr_inv(ot)% nsubset        + n
              bufr_inv(ot)% subseto(ifile) = bufr_inv(ot)% subseto(ifile) + n
            endif
            if (ot == OT_RAD) then
              ! SAT_PP files (containing OT_RAD data) are scanned before in
              ! scan_satpp_feedback@mo_tovs, where the offset of the previous
              ! FT_FEEDBACK files with OT_RAD was not known. Thus, we have to add
              ! this offset here.
              do i = ifile+1, n_source
                if (source(i)% filetype == FT_SATPP) &
                     bufr_inv(ot)% subseto(i) = bufr_inv(ot)% subseto(i) + n
              end do
            end if
          end do
          source(ifile)% entries = fb% n_hdr
          if (fb% n_hdr > 0) then
             if (all (itmp == itmp(1))) then
                source(ifile)% obstype = itmp(1)
             end if
          else
            source(ifile)% filetype = FT_EMPTY
          end if
          deallocate (itmp)
          !---------------------------------------------------------------
          ! Get resolution, horizontal domain size from global atttributes
          !---------------------------------------------------------------
          if (fb% n_hdr > 0) then
!            source(ifile)% model      = fb% model
             source(ifile)% resolution = fb% resolution
             source(ifile)% domain     = fb% domain(1:2) ! nx,ny only
          end if
        !------------------------------------------------
        ! pass = 2 : read the file, store in t_obs/t_spot
        !------------------------------------------------
        else
          !--------------------------------
          ! determine observation box index
          !--------------------------------
          if (size(obs)==1) then
            ib = 1
          else
            ib = ifile
          endif
          bi => obs(ib)
          !-------------------------
          ! allocate header (t_spot)
          !-------------------------
          n  = fb% n_hdr
          if (source(ifile)% obstype > 0) &
            n = min (n, rept_use(source(ifile)% obstype)% max_proc)
          is = bi% n_spot + 1
          ie = bi% n_spot + n
          allocate (i_body (n))
          allocate (l_body (n))
          allocate (i_spec (n))
          allocate (l_spec (n))
          allocate (itmp   (n))
          allocate (check  (n))
          allocate (state  (n))
          allocate (flags  (n))
          call new_spot (bi , n, set_id=.true.)
          s  => bi% spot(is:ie)
          hd => s% hd
          !---------------------------------
          ! set mo_t_netcdf module variables
          !---------------------------------
          ncid  = fb%nc% ncid
          stanc = 1
          counc = n
          strnc = 1
          !---------------------------
          ! set indices header -> body
          !---------------------------
          call get_var   (    i_body,       'i_body'           ,fill=-1)
          call get_var   (    l_body,       'l_body'           ,fill=-1)
          call get_var   (    i_spec,       'i_spec' ,ierr=ierr,fill= 0)
          call get_var   (    l_spec,       'l_spec' ,ierr=ierr,fill= 0)
          s% o% n = l_body
          s% o% i = i_body - i_body(1) + bi% n_obs
          s% s% n = l_spec
          s% s% i = i_spec
          call get_var   (s% col% nlev,     'n_level'      )
          !------------------
          ! consistency check
          !------------------
          if (any (i_body < 0 .or.                            &
                   l_body < 0     ) .or.                      &
              any (i_body(2:n)-i_body(1:n-1) /= l_body(1:n-1))) then
            call finish('read_feedback','file is corrupt: '//trim(filename))
          endif
          !-----------------------
          ! read mandatory entries
          !-----------------------
          call get_var   (hd% buf_type,     'data_category')
          call get_var   (hd% buf_subtype,  'sub_category' )
          call get_var   (hd% center,       'center'       )
          call get_var   (hd% subcenter,    'sub_center'   )
          call get_var   (hd% obstype,      'obstype'      )
          call get_var   (hd% codetype,     'codetype'     )
          call get_var   (hd% dbkz,         'dbkz'     ,ierr=ierr,fill=-1)
          call get_var   (s % ident,        'ident'        )
          where (hd% obstype == OT_GPSRO .or. &
                 hd% obstype == OT_RAD   .or. &
                 hd% obstype == OT_SATOB .or. &
                 hd% obstype == OT_SOIL       ) hd% satid = s% ident
          call get_var   (s % statid,       'statid'       )
          call get_var   (s%col%c% dlat,    'lat'          )
          call get_var   (s%col%c% dlon,    'lon'          )
          call get_var   (s % z,            'z_station'    )
          call get_var   (s % corme,        'sta_corr'     )
          call get_var   (s % stlsf,        'surftype' ,ierr=ierr,fill=-1)
          call get_var   (s % soiltype,     'soiltype' ,ierr=ierr,fill=-1)
          call get_var   (    itmp,         'time'     ,fill=-huge(itmp))
          do i = 1,n
            if (itmp(i)==-huge(itmp)) &
              call finish('read_feedback','fillvalue given for "time"')
            call init_time(minutes,  mi = itmp(i))
            s(i)% actual_time = minutes + ref_time
          end do
          call get_var   (    itmp,         'time_nomi' ,fill=-huge(itmp))
          do i = 1,n
            if (itmp(i)==-huge(itmp)) &
              call finish('read_feedback','fillvalue given for "time_nomi"')
            call init_time(minutes,  mi = itmp(i))
            hd(i)% time = minutes + ref_time
          end do
          call get_var   (    itmp,         'time_dbase',fill=-huge(itmp))
          do i = 1,n
            if (itmp(i) /= -huge(itmp)) then
              call init_time(minutes,  mi = itmp(i))
              hd(i)% db_time = minutes + ref_time
            end if
          end do
          call get_var   (s% sozen,         'sun_zenit' ,fill=invalid)
          !------------------------------------------
          ! processing flags and quality control (qc)
          !------------------------------------------
          if (lqc) then
            call get_var      (state, 'r_state')
            call get_var      (flags, 'r_flags')
            call get_var      (check, 'r_check')
            call modify_fof_flags (state, check, flags, .true.)
            call reverse_code (state,    stats )
            s% use% state    = state
            call reverse_bits (flags,    checks, unknown = CHK_CONSIST)
            s% use% flags    = flags
            call reverse_code (check,    checks, unknown = CHK_CONSIST)
            s% use% check    = check
          endif
          !-------------------------------------
          ! model grid indices & parameters (bg)
          !-------------------------------------
          if(lbg) then
            call get_var (s % col%h%ijdp(1),'index_x')
            call get_var (s % col%h%ijdp(2),'index_y')
!           call get_var (s % col%h%ijdp(3),'index_d',ierr=ierr,fill=1_i2)
            call get_var (s % col%h%ijdp(3),'index_d',ierr=ierr,fill=1   )
            call get_var (s % gp_bg,        'z_modsurf')
            call get_var (s % ssd_bg,       'sso_stdh', ierr=ierr,fill=invalid)
            s % gp_bg =   s % gp_bg * gacc
            call get_var (    itmp,         'mdlsfc' ,ierr=ierr,fill=0)
            if (ierr == NF_NOERR) then
              s % mdlsfc = itmp
              call frac_mdlsfc (s, itmp)
            else
              call get_var (s % sl_bg,      'mdlsf'  ,ierr=ierr,fill=invalid)
            endif
            allocate (rtmp(fb% n_hdr))
            call get_var (rtmp,             'fr_land',ierr=ierr,fill=invalid)
            if (ierr == NF_NOERR) then
              where (rtmp /= invalid) s% sl_bg = rtmp
            end if
            deallocate (rtmp)
          endif
          !-----------------
          ! optional entries
          !-----------------
          call get_var (s% sttyp,    'instype'     ,ierr=ierr,fill=-1   ) ! TEMP RAD
          call get_var (s% stret,    'retrtype'    ,ierr=ierr,fill=-1   ) ! SATOB
          call get_var (   itmp,     'rad_corr'    ,ierr=ierr,fill=-1   ) ! TEMP
          where(itmp/=-1) s% stret = itmp
          if (any (s% hd% obstype == OT_GPSRO)) then
          call get_var (s% tracking, 'sat_class'   ,ierr=ierr,fill=-1   ) ! GPSRO
          call get_var (s% sender_id,'prn'         ,ierr=ierr,fill=-1   ) ! GPSRO
          else
          call get_var (s% tracking, 'tracking'    ,ierr=ierr,fill=-1   ) ! TEMP PILOT
          end if
          call get_var (s% meas_type,'meas_type'   ,ierr=ierr,fill=-1   ) ! TEMP PILOT
          call get_var (s% phase,    'phase'       ,ierr=ierr,fill=-1_i2) ! AIREP RAD
!         where(itmp/=-1) s% phase = itmp
          call get_var (s% stclf,    'flg_cld'     ,ierr=ierr,fill=-1   ) ! RAD
!         call get_var (   itmp,     'varno_back'  ,ierr=ierr,fill=-1   ) ! RADAR
          call get_var (s% params(1),'vnyquist'    ,ierr=ierr,fill=0._sp) ! RADAR
          call get_var (s% spec_rfl, 'spec_r_flags',ierr=ierr,fill= 0   ) ! RADAR
          call get_var (s% stzen,    'sat_zenit'   ,ierr=ierr,fill=invalid) ! RAD
          call get_var (s% stazi,    'sat_azimuth' ,ierr=ierr,fill=invalid) ! RAD
          call get_var (s% soazi,    'sun_azimuth' ,ierr=ierr,fill=invalid) ! RAD
          call get_var (s% center_id,'center_id'   ,ierr=ierr,fill=-1)      ! GPSGB
          !----------------------------
          ! 3dvar specific entries (id)
          !----------------------------
          if(lid) then
            call get_var (hd% id,           'obs_id' ,ierr=ierr,fill=0)
            call get_var (hd% source,       'source' ,ierr=ierr,fill=0)
            call get_var (hd% record,       'record' ,ierr=ierr,fill=0)
            call get_var (hd% subset,       'subset' ,ierr=ierr,fill=0)
          else
            hd% id     = 0
            hd% source = ifile
            hd% record = [(i,i=1,n)]
            hd% subset = 0
            if (present(src)) hd% source = src
          endif
          !--------------
          ! allocate body
          !--------------
          nb  = fb% n_body
          if (n < fb% n_hdr) nb = i_body(n) - i_body(1) + l_body(n)
          isb = bi% n_obs + 1
          ieb = bi% n_obs + nb
          call new_obs (bi , nb)
          bi% n_obs = ieb
          bd => bi %body (isb:ieb)
          deallocate (itmp,     state,     check,     flags    )
          allocate   (itmp(nb), state(nb), check(nb), flags(nb))
          !--------------------------
          ! remember position in file
          !--------------------------
          i0 = 0
          do i = 1, n
            hd(i)% mon_file = ifile
            hd(i)% mon_rec  = i
            ni = s(i)% o% n
            in = i0 + ni
            bd (i0+1:in)% mon_pos = (/(j,j=1,ni)/)
            i0 = in
          end do
          !---------------------------------
          ! set mo_t_netcdf module variables
          !---------------------------------
          ncid  = fb%nc% ncid
          stanc = 1
          counc = nb
          strnc = 1
          !-----------------------
          ! read mandatory entries
          !-----------------------
          call get_var   (bi %   varno(isb:ieb),'varno'    )
          call get_var   (bd %   o,             'obs'      )
          call get_var   (bd %   bc,            'bcor'     )
          call get_var   (bi %   olev (isb:ieb),'level'    )
          call get_var   (bd %   plev,          'plevel'   ,fill=-1._sp)
          call get_var   (bd %   lev_typ, 'level_typ')
          call get_var   (bd %   lev_sig, 'level_sig')
          call get_var       (   state,   'state'    )
          call get_var       (   flags,   'flags'    )
          call get_var       (   check,   'check'    )
          !------------------------------------------------
          ! convert to 3dvar convention, set plev from olev
          !------------------------------------------------
          call modify_fof_flags (state, check, flags, .false.)
          call reverse_code  (   state,    stats     )
          bd% use% state     =   state
          call reverse_bits  (   flags,    checks    ,unknown = CHK_CONSIST)
          bd% use% flags     =   flags
          call reverse_code  (   check,    checks    ,unknown = CHK_CONSIST)
          bd% use% check     =   check
          where (bd% lev_typ == VN_P .and. bd %plev < 0._sp) &
            bd %plev = bi %olev (isb:ieb)
          !----------------------
          ! read optional entries
          !----------------------
          call get_var   (itmp, 'qual'      ,ierr=ierr,fill=-1)
          if (ierr /= NF_NOERR) &
            call get_var (itmp, 'pcc'       ,ierr=ierr,fill=-1)
          bd % pcc = itmp
          call get_var   (itmp, 'spec_index',ierr=ierr,fill=-1)
          bd % spec_index = itmp
          call get_var (bd% lat,        'dlat'      ,fill=inv_datum% lat)
          call get_var (bd% lon,        'dlon'      ,fill=inv_datum% lon)
          call get_var (bd% obs_par(1), 'azimuth'   ,fill=inv_datum% obs_par(1))
          call get_var (bd% plev_width, 'plev_width',fill=inv_datum% plev_width)
          call get_var (bd% set% v_loc, 'v_loc'     ,fill=inv_datum% set% v_loc)
          call get_var (bd% set% h_loc, 'h_loc'     ,fill=inv_datum% set% h_loc)
          call get_var (bd% ac,         'accuracy'  ,fill=inv_datum% ac)
          !----------------------------
          ! read 3dvar specific entries
          !----------------------------
          call get_var (bd % eo,  'e_o' ,ierr=ierr,fill=0._sp)
          call get_var (bd % wqc, 'w_qc',ierr=ierr,fill=0._sp)
          !-------------------------------------
          ! read observation type specific table
          !-------------------------------------
          if (any (hd% obstype == OT_RAD)) then
            if (lspec) then
              call read_fdbk_rad('read',obs(ib:ib), i_s=is, i_e=ie, nh=fb%n_hdr, nb=nb)
            end if
          endif
          if (any (hd% obstype == OT_RADAR)) then
            call read_fdbk_radar (bi, fb% n_radar)
          endif
          if (any (hd% obstype == OT_GPSRO)) then
            call read_fdbk_occ (bi)
          endif
          if (any (hd% obstype == OT_GPSGB)) then
            call read_fdbk_gpsgb (bi)
          endif
          if (any (hd% obstype == OT_TEMP)) then
            call read_fdbk_temp (bi)
          endif
          !------------------------
          ! modify fof-file content
          !------------------------
          call modify_obs_error (bd% eo, bi% varno(isb:ieb))
          !-----------------------
          ! read verification data
          !-----------------------
          call get_veri_index (idx, nidx, fb, ens_member=VE_VQC_WEIGHT)
          if (idx(1)>0) bd% wqc = get_veri (fb, idx(1))
          where (bd% wqc == NF_FILL_REAL) bd% wqc = 0._sp
          if (replace_o_fg) then
            call get_veri_index (idx, nidx, fb, ens_member=VE_DETERM,   &
                                                  run_type=VT_FIRSTGUESS)
            if (idx(1)<=0) call finish('read_feedback','no fg present')
            bd% o  = get_veri (fb, idx(1))
            bd% bc = 0._sp
          endif
          !--------------------------
          ! read 3dvar internal flags
          !--------------------------
          hd% modtype = rept_char(hd% obstype)% mod
          call set_xuv (s% col)
!      id           =
          do i = 1, n
            !------------------------------
            ! set 3dvar internal parameters
            !------------------------------
            if (hd(i)% dbkz < 0) then
              obsinfo   = obstype_code (hd(i)% obstype,  &
                                        hd(i)% codetype, &
                                 centre=WMO0_ECMWF       )
              if (obsinfo% dbkz < 0) &
                obsinfo = obstype_code (hd(i)% obstype,  &
                                        hd(i)% codetype, &
                                 centre=WMO0_DWD         )
              hd(i)% dbkz = obsinfo% dbkz
            endif
            if (hd(i)% dbkz >= 0) then
              hd(i)% idbk = idb_dbk (hd(i)% dbkz, hd(i)% obstype)
            endif
            !----------------------------------------------
            ! station pressure (for thinning, height check)
            ! single level data only
            ! report quality information for SATOB
            !----------------------------------------------
            j1 = s(i)% o% i + 1
            jn = s(i)% o% i + s(i)% o% n
            if (all (bi% olev (j1:jn) == bi% olev (j1))) then
              if (bi% body (j1)% lev_typ == VN_P) then
                s(i)% ps = bi% olev (j1)
              else if     (bi% body (j1)% plev > 0._sp) then
                s(i)% ps = bi% body (j1)% plev
              endif
            endif
            select case (hd(i)% obstype)
            case (OT_SATOB)
              s(i)% pcc = bi% body (j1)% pcc
            end select
            !------------------------------------------------
            ! wind components: revert bug up to version 1.01
            ! rotate wind back from rotated model coordinates
            !------------------------------------------------
            if (luvbug) then
              do j = j1, jn-1
                select case (bi% varno(j))
                case (VN_U, VN_U10M)
                  select case (bi% varno(j+1))
                  case (VN_V, VN_V10M)
                    call uvrot2uv (real (bi% body(j)% o, wp), real (bi% body(j+1)% o, wp), &
                                   s(i)%col%c% dlat,          s(i)%col%c% dlon,            &
                                   real (fb% pole(1), wp),    real (fb% pole(2), wp),      &
                                   u, v                                                    )
                    bi% body(j)  % o = u
                    bi% body(j+1)% o = v
                  end select
                end select
              end do
            endif
          end do
          if (lchk) then
             !-----------------------
             ! Apply checks and rules
             !-----------------------
             do i = 1, n
                call check_report_0 (s(i)% use, hd(i), 1, keep=.true.)
                call check_report_1 (s(i))
             end do
          end if
          !---------
          ! clean up
          !---------
          deallocate (i_body)
          deallocate (l_body)
          deallocate (i_spec)
          deallocate (l_spec)
          deallocate (itmp, state, check, flags)
        endif

        !-----------
        ! close file
        !-----------
        call close_fdbk   (fb)
        call cleanup_fdbk (fb)

      endif
    enddo

    if (pass == 1) then
      do ifile = 1, n_source
         !----------------------------
         ! skip if not a feedback-file
         !----------------------------
         if (source(ifile)% filetype /= FT_FEEDBACK) cycle
         !--------------------------
         ! PE where file was scanned
         !--------------------------
         pe = dace% pio
         !--------------------------------------------------------------
         ! Update components of array 'source' not known at add_source()
         !--------------------------------------------------------------
         call p_bcast (source(ifile)% obstype,    pe)
!        call p_bcast (source(ifile)% model,      pe)
         call p_bcast (source(ifile)% resolution, pe)
         call p_bcast (source(ifile)% domain,     pe)
      enddo
    end if

  end subroutine read_feedback
!==============================================================================
  subroutine read_nml_fof_input
  !--------------------------
  ! read namelist /FOF_INPUT/
  !--------------------------
    integer :: ierr
    !-------------
    ! set defaults
    !-------------
    flag_dataset    = .false.
    rm_fg_check     = .false.
    force_active    = .false.
    check_member    = .false.
    replace_o_fg    = .false.
    percent_fg_check= 0.      ! % of fg checks required for rejection
    ps_obs_error    = 0.      ! <0 factor;>0 value of PS obs error
    uv_obs_error    = 0.      ! <0: factor on, >0: value of wind obs.error
    uv10m_obs_error = 0.      ! <0: factor on, >0: value of wind obs.error
    radvel_obs_error= 0.      ! <0: factor on, >0: value of radial wind oe
    refl_obs_error  = 0.      ! <0: factor on, >0: value of radar reflectivity oe
    t_obs_error     = 0.      ! <0: factor on, >0: value of temp.obs.error
    t2m_obs_error   = 0.      ! <0: factor on, >0: value of temp.obs.error
    rh_obs_error    = 0.      ! <0: factor on, >0: value of hum. obs.error
    rh2m_obs_error  = 0.      ! <0: factor on, >0: value of hum. obs.error
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('FOF_INPUT', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=FOF_INPUT, iostat=ierr)
        if (ierr/=0) call finish ('read_psas_nml',               &
                                  'ERROR in namelist /FOF_INPUT/')
#else
        read (nnml ,nml=FOF_INPUT)
#endif
      end select
      !---------
      ! printout
      !---------
      write (6,'(a     )') repeat('-',79)
      write (6,'(      )')
      write (6,'(a     )') '  namelist /FOF_INPUT/'
      write (6,'(      )')
      write (6,'(a,l6  )') '  flag_dataset     = ',flag_dataset
      write (6,'(a,l6  )') '  rm_fg_check      = ',rm_fg_check
      write (6,'(a,l6  )') '  force_active     = ',force_active
      write (6,'(a,l6  )') '  check_member     = ',check_member
      write (6,'(a,l6  )') '  replace_o_fg     = ',replace_o_fg
      write (6,'(a,f6.2)') '  percent_fg_check = ',percent_fg_check
      write (6,'(a,f6.2)') '  ps_obs_error     = ',ps_obs_error
      write (6,'(a,f6.2)') '  uv_obs_error     = ',uv_obs_error
      write (6,'(a,f6.2)') '   t_obs_error     = ',t_obs_error
      write (6,'(a,f6.2)') '  rh_obs_error     = ',rh_obs_error
      write (6,'(a,f6.2)') ' uv10m_obs_error   = ',uv10m_obs_error
      write (6,'(a,f6.2)') '   t2m_obs_error   = ',t2m_obs_error
      write (6,'(a,f6.2)') '  rh2m_obs_error   = ',rh2m_obs_error
      write (6,'(a,f6.2)') '  radvel_obs_error = ',radvel_obs_error
      write (6,'(a,f6.2)') '  refl_obs_error   = ',refl_obs_error
    endif
    !---------------------------
    ! broadcast namelist entries
    !---------------------------
    call p_bcast (flag_dataset     ,dace% pio)
    call p_bcast (rm_fg_check      ,dace% pio)
    call p_bcast (force_active     ,dace% pio)
    call p_bcast (check_member     ,dace% pio)
    call p_bcast (replace_o_fg     ,dace% pio)
    call p_bcast (percent_fg_check ,dace% pio)
    call p_bcast (ps_obs_error     ,dace% pio)
    call p_bcast (uv_obs_error     ,dace% pio)
    call p_bcast (radvel_obs_error ,dace% pio)
    call p_bcast (refl_obs_error   ,dace% pio)
    call p_bcast ( t_obs_error     ,dace% pio)
    call p_bcast (rh_obs_error     ,dace% pio)
    call p_bcast (uv10m_obs_error  ,dace% pio)
    call p_bcast ( t2m_obs_error   ,dace% pio)
    call p_bcast (rh2m_obs_error   ,dace% pio)
  end subroutine read_nml_fof_input
!------------------------------------------------------------------------------
  subroutine modify_fof_flags (state, check, flags, report)
  integer, intent(inout) :: state (:) ! status   variable to modify
  integer, intent(inout) :: check (:) ! check    variable to modify
  integer, intent(inout) :: flags (:) ! bit flag variable to modify
  logical, intent(in)    :: report    ! apply to report or single observations
  !------------------------------------------------------
  ! optionally modify flags according to namelist setting
  !------------------------------------------------------
    integer :: i, j
    !------------------------------------------
    ! reset the effect of the first guess check
    !------------------------------------------
    if (rm_fg_check .and. .not. report) then
      !-----------------------
      ! loop over observations
      !-----------------------
      do i=1,size(state)
        !-------------------------------------
        ! action required if FL_FG flag is set
        !-------------------------------------
        if (btest (flags(i), FL_FG)) then
          flags(i) = ibclr (flags(i), FL_FG)
          if (flags(i) == 0) then
            !----------------------
            ! only FG flag was set:
            !   set to ACTIVE
            !----------------------
            select case (state(i))
            case (ST_REJECTED)
              check(i) = FL_NONE
              state(i) = ST_ACTIVE
              case default
                call finish('rm_fg_check',                    &
                             name_value (status_t, state(i))//&
                        ' '//name_value ( flags_t, check(i))//&
                        ' '//char5      (          flags(i))  )
            end select
          else
            !--------------------------
            ! some other flags were set
            !--------------------------
            select case (state(i))
            case (ST_REJECTED)
              !---------------------
              ! REJECTED -> REJECTED
              !---------------------
            case (ST_PAS_REJ)
              !-------------------
              ! PAS_REJ -> PASSIVE
              !-------------------
              state(i) = ST_PASSIVE
            case default
              call finish('rm_fg_check',                    &
                           name_value (status_t, state(i))//&
                      ' '//name_value ( flags_t, check(i))//&
                      ' '//char5      (          flags(i))  )
            end select
            !------------------------
            ! set appropriate 'check'
            !------------------------
            check(i) = FL_NONE
            do j = 0, 31
              if (btest (flags(i), j)) check(i) = j
            end do
          endif
        endif
      end do
    endif
    !-------------------------------------
    ! force passive observations to active
    !-------------------------------------
    if (force_active) then
      !-----------------------
      ! loop over observations
      !-----------------------
      do i=1,size(state)
        if (state(i) == ST_PASSIVE) then
          state(i) = ST_ACTIVE
          check(i) = FL_NONE
        endif
      end do
    endif
    !-------------------------------------------------------------
    ! set all checks to 'FL_DATASET'
    ! for debugging, so that any degradation of the state variable
    ! present in the input file may be traced back
    !-------------------------------------------------------------
    if (flag_dataset) then
      where (state > ST_ACTIVE) check = FL_DATASET
    endif
  end subroutine modify_fof_flags
!------------------------------------------------------------------------------
  subroutine modify_obs_error (eo, varno)
  real(sp) ,intent(inout) :: eo   (:) !  observation error
  integer  ,intent(in)    :: varno(:) ! variable number
  !-----------------------------------------
  ! modify observational error from fof-file
  !-----------------------------------------
    if (ps_obs_error < 0.) then
      where (varno == VN_PS) eo = eo * abs (ps_obs_error)
    endif
    if (ps_obs_error > 0.) then
      where (varno == VN_PS) eo =           ps_obs_error
    endif

    if (uv_obs_error < 0.) then
      where (varno == VN_U .or. varno == VN_V) &
                             eo = eo * abs (uv_obs_error)
    endif
    if (uv_obs_error > 0.) then
      where (varno == VN_U .or. varno == VN_V) &
                             eo =           uv_obs_error
    endif

    if (uv10m_obs_error < 0.) then
      where (varno == VN_U10M .or. varno == VN_V10M) &
                             eo = eo * abs (uv10m_obs_error)
    endif
    if (uv10m_obs_error > 0.) then
      where (varno == VN_U10M .or. varno == VN_V10M) &
                             eo =           uv10m_obs_error
    endif

    if (radvel_obs_error < 0.) then
      where (varno == VN_RADVEL) &
                             eo = eo * abs (radvel_obs_error)
    endif
    if (radvel_obs_error > 0.) then
      where (varno == VN_RADVEL) &
                             eo =           radvel_obs_error
    endif

    if (refl_obs_error < 0.) then
      where (varno == VN_RREFL) &
                             eo = eo * abs (refl_obs_error)
    endif
    if (refl_obs_error > 0.) then
      where (varno == VN_RREFL) &
                             eo =           refl_obs_error
    endif

    if ( t_obs_error < 0.) then
      where (varno == VN_T) eo = eo * abs ( t_obs_error)
    endif
    if ( t_obs_error > 0.) then
      where (varno == VN_T) eo =            t_obs_error
    endif

    if ( t2m_obs_error < 0.) then
      where (varno == VN_T2M) eo = eo * abs (t2m_obs_error)
    endif
    if ( t2m_obs_error > 0.) then
      where (varno == VN_T2M) eo =           t2m_obs_error
    endif

    if ( rh_obs_error < 0.) then
      where (varno == VN_RH) eo = eo * abs ( rh_obs_error)
    endif
    if ( rh_obs_error > 0.) then
      where (varno == VN_RH) eo =            rh_obs_error
    endif

    if ( rh2m_obs_error < 0.) then
      where (varno == VN_RH2m) eo = eo * abs ( rh2m_obs_error)
    endif
    if ( rh2m_obs_error > 0.) then
      where (varno == VN_RH2m) eo =            rh2m_obs_error
    endif

  end subroutine modify_obs_error
!==============================================================================
end module mo_fdbk_in
