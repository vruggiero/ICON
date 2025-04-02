!
!+ Write feedback files (new format, common with COSMO and LETKF)
!
MODULE mo_fdbk_3dvar
!
! Description:
!   Write feedback files (new format, common with COSMO)
!   3DVAR/LETKF specific part
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  changes for lat/lon grids
! V1_5         2009/05/25 Andreas Rhodin
!  update feedback file interface
! V1_8         2009/12/09 Andreas Rhodin
!  Changes for verification: merge analysis data into monitoring file
! V1_9         2010/04/20 Andreas Rhodin
!  dont write dismissed reports to feedback file
!  bugfix in write_fdbk_3dv (array bound violation when writing cof-file)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  write correct veri_ref_date, veri_forecast_time as clarified in the doc.
!  write_fdbk_3dv: only write observations with status > DISMISSED
!  add_veri_3dv: bugfix (valid for 'feedback=8')
! V1_13        2011/11/01 Andreas Rhodin
!  ensemble output
! V1_15        2011/12/06 Andreas Rhodin
!  write (optional) body variable 'plevel' to feedback file
! V1_17        2011/12/21 Andreas Rhodin
!  for feedback=7: merge psas-analysis into mon-files
!  remove temporary fix for inconsistent labeling of channels
! V1_18        2012/01/06 Andreas Rhodin
!  bugfix (bcast array of different size on src and dest)
! V1_19        2012-04-16 Harald Anlauf
!  convert to p_gather(v)
!  revised IASI diagnostics in feedback file
!  write lin_ana (Ya linear operator on analysis)
! V1_20        2012-06-18 Andreas Rhodin
!  correctly write 'level_typ' to feedback file for KENDA
!  write_fdbk_3dv: properly round obs.time, heights to nearest integer
! V1_22        2013-02-13 Andreas Rhodin
!  changes for RADAR observation operator
! V1_23        2013-03-26 Andreas Rhodin
!  write localisation length scales (v_loc/h_loc) to LETKF feedback files
! V1_24        2013/04/26 Andreas Rhodin
!  fix bug (optional parameter not present, crash in routine V1_23 20130425)
! V1_26        2013/06/27 Harald Anlauf
!  Determine effective resolution for ICON grid written to feedback files.
!  Fallback to 'model' from nml /run/ if actual parameter is not present
! V1_27        2013-11-08 Andreas Rhodin
!  write variational bias correction
!  correctly write optional parameters for PSAS, PSAS+LETKF, VARENKF
! V1_31        2014-08-21 Andreas Rhodin
!  changes for MEC
! V1_35        2014-11-07 Andreas Rhodin
!  changes for MEC and GPSGB
! V1_42        2015-06-08 Andreas Rhodin
!  write variables 'dlat',dlon','azimuth' to feedback files (for GPSRO/TEMP)
! V1_43        2015-08-19 Andreas Rhodin
!  write_fdbk_3dv: abort if feedback file cannot be created
! V1_45        2015-12-15 Andreas Rhodin
!  replace pl_width (integer) with plev_width (real)
! V1_48        2016-10-06 Robin Faulwetter
!   Option for splitting feedback files according to satellite ID and grid ID.
! V1_50        2017-01-09 Robin Faulwetter
!  Minor bugfix for splitted feedback file writing.
! V1_51        2017-02-24 Andreas Rhodin
!  changed comment
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Author:
! Andreas Rhodin  DWD  2008  original code
! Harald Anlauf   DWD  2011
!==============================================================================

!-------------
! Modules used
!-------------
use mo_kind,          only: sp, wp, i2, i4    ! precision kind parameters
use mo_exception,     only: finish            ! return in case of errors
use mo_version,       only: var3d_version     ! fctn. returning version number
use mo_dace_string,   only: char4             ! Conversion: INT -> CHAR(LEN=4)
use mo_physics,       only: gacc              ! gravity acceleration
use mo_run_params,    only: path_file,       &! concatenate path/filename
                            aux,             &! path
                            run_time,        &! run time
                            ana_time,        &! analysis reference time
                            fc_ref_time,     &! forecast start time
                            run_type,        &! haupt=0,(vor=1),ass=2,test=3
                            nex,             &! experiment id
                            fdbk_addvar,     &! optional fdbk output parameters
                  rmodel => model             ! 'GME' or 'COSMO' from nml /run/
!                rmethod => method            ! 'PSAS', 'LETKF', ...
use mo_wmo_tables,    only: WMO6_LATLON,     &! latitude/Longitude         grid
                            WMO6_GAUSSIAN,   &! Gaussian                   grid
                            WMO6_ROTLL,      &! rotated latitude/longitude grid
                            DWD6_ICOSAHEDRON,&! icosahedral triangular     grid
                            DWD6_ICON         ! ICON (triangular)          grid
use mo_time,          only: t_time,          &! derived type for date+time
                            iyyyymmdd,       &! extract date (as integer)
                            ihhmm,           &! extract time (as integer)
                            cyyyymmddhhmm,   &! extract time (as character)
                            minutes,         &! convert to minutes
                            operator(-),     &! calculate time difference
                            operator(/=),    &! compare times
                            zero_time         ! invalid value
use mo_t_obs,         only: t_obs,           &! observation data type
                            t_spot,          &! header data type
                            add_source,      &! add source-file to list
!                           source,          &! list of source-files
                            n_source,        &! number of source files
                            FT_FEEDBACK,     &! flag for feedback file
                            FT_MISSING,      &! flag for missing  file
                            mstep,           &! # of feedback file types (cof, mon, ekf, ver)
                            pref_step,       &! prefix  for feedback files
                            comm_step,       &! comment for feedback files
                            fdbk_split        ! split feedback files
use mo_wigos,         only: wsi_to_text       ! Convert t_wsi to text format
use mo_obs_set,       only: t_obs_set         ! observation data type
use mo_t_datum,       only: t_datum           ! body   data type
use mo_t_use,         only: chk,             &! codes for 'checks'
                            stats,           &! codes for 'states'
                            STAT_DISMISS      ! flag for states to dismiss
use mo_radar,         only: t_radar,         &! radar obstype specific table entry
                            radar_tables      ! prepare the radar specific tables
use mo_atm_grid,      only: t_grid            ! grid derived type
use mo_dec_matrix,    only: t_vector,        &! type definition: real vector
                            t_bvector,       &!               boolean vector
                            construct,       &! derived type constructor routine
                            destruct,        &! derived type  destructor routine
                            delete_storage    ! delete function return args.
use mo_allocate,      only: enter_function,  &! called at start of routines
                            leave_function    ! called at end of routines
use mo_mpi_dace,      only: dace,            &! MPI group info
                            p_bcast,         &! generic MPI broadcast routine
                            p_gather,        &! generic MPI_gather    routine
                            p_gatherv,       &! generic MPI_gatherv   routine
                            p_max,           &! generic MPI maximum   routine
                            p_or,            &! generic MPI or        routine
                            p_ior,           &! generic MPI ior       routine
                            p_allgather       ! generic MPI_allgather routine
use mo_fdbk,          only: t_fdbk,          &! feedback file data type
                            setup_fdbk,      &!
                            create_fdbk,     &!
                            open_fdbk_write, &!
                            read_meta,       &! read  attributes + meta data
                            close_fdbk,      &!
                            cleanup_fdbk,    &!
                            add_verification,&!
                            get_veri_index,  &!
                            get_varid,       &!
                            check_all_fdbk_addvar,&!
                            write_global_attributes
use mo_t_table,       only: name_value        ! find name of table entry
use mo_fdbk_tables,   only: init_fdbk_tables,&! set up tables
                            obstype,         &! observation type table
                            VT_FIRSTGUESS, VT_PREL_ANA, VT_ANALYSIS,   &
                            VT_LIN_ANA,                                &
                            VE_DETERM, VE_BG_ERROR, VE_VQC_WEIGHT,     &
                            VE_ENS_MEAN, VE_ENS_SPREAD, VE_TALAGRAND,  &
                            VE_VQC_WEIGHT, VE_ENS_MEAN_OBS, VE_BIASCOR,&
                            VN_NUM, VN_IMPPAR, VN_P,                   &
                            ST_ACCEPTED, ST_ACTIVE,                    &
                            OT_RAD, OT_GPSRO,&! observation type ids.
                            OT_RADAR,        &! volume radar obstype id.
                            n_ot,            &! number of observation types
                            c_optv, n_optv    ! TOVS specific optional fdbk-file vars
use mo_t_netcdf_file, only: t_netcdf_var,    &! derived type definition
                            PLEN              ! max. length of filename/path
use mo_t_tovs,        only: t_tovs,          &! type to store radiance specific stuff
                            t_tovs_instr,    &! info on instruments in t_tovs
                            load, destruct    ! load/destuct t_tovs
use mo_rad,           only: t_rad_set,       &!
                            rad_set,         &!
                            n_set,           &!
                            set_indx,        &!
                            OPTV_L2C,        &!
                            OPTV_NWC_FLG,    &!
                            OPTV_ORB_PH ,    &!
                            OPTV_INS_TMP,    &!
                            OPTV_CLD_FRC,    &!
                            OPTV_CLD_FLG
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,         only: NF90_INT,        &!
                            NF90_SHORT,      &!
                            NF90_BYTE,       &!
                            NF90_REAL,       &!
                            NF90_CHAR,       &!
                            NF90_FILL_REAL,  &!
                            nf90_enddef,     &!
                            nf90_inq_varid,  &!
                            nf90_get_var,    &!
                            nf90_put_var,    &!
                            nf90_strerror,   &!
                            NF90_NOERR        ! status return value: no error


implicit none

!----------------
! Public entities
!----------------
private
public :: t_fdbk_3dv      ! 3dvar specific feedback file data type
public :: destruct        ! t_fdbk_3dv destructor routine
public :: write_fdbk_3dv  ! write 3dvar feedback file
public :: add_veri_3dv    ! add verification data

!-----------
! interfaces
!-----------
interface destruct
  module procedure destruct_fdbk_3dv
end interface destruct

!================================================
! feedback file meta data derived type definition
!================================================

  !--------------------------------------------------------------------
  ! Feedback file information is split over several files. File IDs are
  ! related to observation types by means of derived type t_obst:
  !--------------------------------------------------------------------
  type t_obst
    integer :: obstype  ! observation type (TEMP, SYNOP, ...)
    integer :: satid    ! satellite id
    integer :: grid     ! instrument grid
    integer :: file     ! file index
  end type t_obst

  !----------------------------------------------------------------------
  ! Feedback file information is split over several files. File names and
  ! dimension sizes are stored in derived type t_file.
  ! The offset in the header and body arrays are stored for each processor
  ! element.
  !----------------------------------------------------------------------
  type t_file
    character(len=14) :: name        = '' ! file name
    character(PLEN)   :: path        = '' ! path name + file name
    integer           :: n_hdr       =  0 ! size of header     (dimension n_hdr)
    integer           :: n_body      =  0 ! size of body      (dimension n_body)
    integer           :: n_spec(n_ot)=  0 ! size of observation specific table
    integer           :: n_veri      =  0 ! number of verification runs (n_veri)
    integer           :: o_hdr       =  0 ! offset in header (for data from ..
    integer           :: o_body      =  0 ! offset in body   ( .. this PE)
    integer           :: o_spec(n_ot)=  0 ! offset in observation specific table
  end type t_file

  !-----------------------------------------------------------------------
  ! In the 3dvar observations are spread over processor elements. They are
  ! stored in 'boxes' (with a typical size of 500 observations) and
  ! 'spots' (corresponding to one report or field of view). For each
  ! spot/box the corresponding output file and the position therein is
  ! given in derived type t_index. This information is held locally for
  ! each processor element.
  !-----------------------------------------------------------------------
  type t_index
    integer :: file   ! file index
    integer :: i_hdr  ! report header number
    integer :: o_body ! report body   offset
    integer :: n_body ! report body   length
    integer :: o_spec ! report table  offset
    integer :: n_spec ! report table  length
    integer :: ib     ! 3dvar  box number
    integer :: is     ! 3dvar spot number
  end type t_index

  type t_fdbk_3dv
    private
    type (t_file)  ,pointer :: file      (:)   => NULL() ! table of files
    type (t_obst)  ,pointer :: file_obst (:)   => NULL() ! map obstype -> file
    type (t_file)  ,pointer :: pelist    (:,:) => NULL() ! fileentries
    type (t_index) ,pointer :: index     (:)   => NULL() ! indices
    type (t_fdbk)           :: fb                        ! 3dvar independent ..
  end type t_fdbk_3dv                                    ! .. part

contains
!------------------------------------------------------------------------------
  subroutine destruct_fdbk_3dv (fdbk)
  type(t_fdbk_3dv) ,intent(inout) :: fdbk
  !------------------------------
  ! t_fdbk_3dv destructor routine
  !------------------------------
    if (associated (fdbk% file))      deallocate (fdbk% file)
    if (associated (fdbk% file_obst)) deallocate (fdbk% file_obst)
    if (associated (fdbk% pelist))    deallocate (fdbk% pelist)
    if (associated (fdbk% index))     deallocate (fdbk% index)
    call close_fdbk   (fdbk% fb)
    call cleanup_fdbk (fdbk% fb)
  end subroutine destruct_fdbk_3dv
!------------------------------------------------------------------------------
  subroutine setup_filenames (fb3d, basename, obs, step)
  type (t_fdbk_3dv) ,intent(inout) :: fb3d
  character(len=*)  ,intent(in)    :: basename
  type(t_obs_set)   ,intent(in)    :: obs      ! observational data
  integer           ,intent(in)    :: step
  !---------------------------------------------------------------
  ! set up table with file basenames
  ! and corresponding 3dvar observation types
  ! either:                      1   observation type  per    file
  ! or (if basename is present): all observation types into 1 file
  !---------------------------------------------------------------
    type (t_file)  ,pointer     :: file_dum      (:)   => NULL()
    type (t_obst)  ,pointer     :: file_obst_dum (:)   => NULL()
    integer        ,parameter   :: nsplit         = 100
    integer                     :: satid  (nsplit)= -1
    integer                     :: grid   (nsplit)= -1
    integer                     :: satid_r(dace% npe)
    integer                     :: grid_r (dace% npe)
    integer                     :: satid_, grid_
    integer                     :: i, j, k
    integer                     :: io  ! 3dvar observation type
    integer                     :: i_f ! file index
    integer                     :: n   ! number of files for obstype
    integer                     :: nm  ! number of files for obstype
    integer                     :: mf  ! max. number of files
    logical                     :: l_sat, l_grid, l_new

    type(t_obs)    ,pointer     :: o
    integer                     :: ib, is

    allocate (fb3d% file_obst (obstype% n))
    if (basename /= '') then
      !----------------------------------
      ! all observation types into 1 file
      !----------------------------------
      allocate (fb3d% file (1))
      fb3d%           file (1)% name = basename
      do io = 1, obstype% n
        fb3d% file_obst(io)% obstype = io
        fb3d% file_obst(io)% satid   = -1
        fb3d% file_obst(io)% grid    = -1
        fb3d% file_obst(io)% file    = 1
      end do
    else
      !----------------------------
      ! 1 observation type per file
      !----------------------------
      allocate (fb3d% file (obstype% n))
      i_f = 0
      do io = 1, obstype% n
        l_sat  = iand(fdbk_split(step,io), 1) /= 0
        l_grid = iand(fdbk_split(step,io), 2) /= 0
        n = 1
        satid(:) = -1
        grid(:)  = -1
        if (l_sat .or. l_grid) then
          !---------------------------------------
          ! 1 satellite and/or instr.grid per file
          !---------------------------------------
          do ib = 1, size(obs% o)                           ! loop over 'boxes'
            if (obs% o(ib)% pe /= dace% pe) cycle           ! handled on this PE only
            o => obs%  o(ib)
            do is = 1, o% n_spot                         ! loop over reports
              if (o% spot(is)% hd% obstype /= io) cycle  ! observation type
              l_new = .false.
              if (l_sat) then
                satid_ = o% spot(is)% hd% satid
                l_new  = all(satid(1:n) /= satid_)
              end if
              if (l_grid) then
                grid_ = o% spot(is)% hd% grid_id
                l_new = all(grid(1:n) /= grid_) .or. l_new
              end if
              if (l_new) then
                if (n>nsplit) call finish('setup_filenames', 'too many files for obstype')
                if (l_sat)  satid(n) = satid_
                if (l_grid) grid(n)  = grid_
                n = n + 1
              end if
            end do
          end do
          ! Communicate satids and grids to all other pes
          nm = p_max(n)
          do i = 1, nm
            call p_allgather(satid(i), satid_r(:))
            call p_allgather(grid(i),  grid_r(:))
            do j = 1, dace% npe
              if (satid_r(j) >= 0 .or. grid_r(j) >= 0) then
                if (.not.any(satid(1:n)==satid_r(j) .and. grid(1:n)==grid_r(j))) then
                  if (n>nsplit) call finish('setup_filenames', 'too many files for obstype')
                  satid(n) = satid_r(j)
                  grid (n) = grid_r (j)
                  n = n + 1
                end if
              end if
            end do
          end do
          ! Sort satid/grid combinations
          n = n - 1
          do i = n-1, 1, -1
            do j = 1, i
              if (satid(j+1) < satid(j) .or. satid(j+1) == satid(j) .and. grid(j+1) < grid(j)) then
                k          = satid(j+1)
                satid(j+1) = satid(j)
                satid(j)   = k
                k          = grid (j+1)
                grid (j+1) = grid (j)
                grid (j)   = k
              end if
            end do
          end do
        end if
        mf = size(fb3d% file)
        if (i_f + n > mf) then
          mf = (i_f + n) + (obstype% n - io)
          allocate(file_dum(mf), file_obst_dum(mf))
          file_dum     (1:i_f) = fb3d% file     (1:i_f)
          file_obst_dum(1:i_f) = fb3d% file_obst(1:i_f)
          deallocate(fb3d% file, fb3d% file_obst)
          fb3d% file_obst => file_obst_dum
          fb3d% file      => file_dum
        end if
        fb3d% file_obst(i_f+1:i_f+n)% obstype = io
        fb3d% file_obst(i_f+1:i_f+n)% satid   = satid(1:n)
        fb3d% file_obst(i_f+1:i_f+n)% grid    = grid (1:n)
        fb3d% file_obst(i_f+1:i_f+n)% file    = (/ (i, i=i_f+1,i_f+n) /)
        fb3d% file     (i_f+1:i_f+n)% name    = name_value (obstype, io)
        do i = 1, n
          if (satid(i) >= 0) then
            k = len_trim(fb3d% file(i_f+i)% name)
            write(fb3d% file(i_f+i)% name(k+1:k+3), '(I3.3)') satid(i)
          end if
          if (grid(i) >= 0) then
            k = len_trim(fb3d% file(i_f+i)% name)
            write(fb3d% file(i_f+i)% name(k+1:k+3), '(I3.3)') grid(i)
          end if
        end do
        i_f = i_f + n
      end do
    endif

  end subroutine setup_filenames
!==============================================================================
  subroutine write_fdbk_3dv (fdbk, obs, reftime, grid, step, name, prefix, &
     method, e_o, w_qc, fg, ana, ana_psas, varbc, e_fg, s_fg, s_ana, m_fg, &
     m_ana, mo_fg, mo_ana, tg_fg, tg_ana, en_fg, en_ana, lin_ana, comment  )
  type(t_fdbk_3dv) ,intent(out), target  :: fdbk     ! feedback file data
  type(t_obs_set)  ,intent(inout)        :: obs      ! observational data
  type(t_time)     ,intent(in)           :: reftime  ! reference time
  type(t_grid)     ,intent(in)           :: grid     ! model grid info
  integer          ,intent(in)           :: step     ! 1:mon 2:ana 3:enkf 4:ver
  character(len=*) ,intent(in)           :: name     ! file basename
  character(len=*) ,intent(in) ,optional :: prefix   ! file prefix
  character(len=*) ,intent(in) ,optional :: method   ! 'PSAS','LETKF'
  type(t_vector)   ,intent(in) ,optional :: e_o      ! observational error
  type(t_vector)   ,intent(in) ,optional :: w_qc     ! VQC weight
  type(t_vector)   ,intent(in) ,optional :: fg       ! first guess
  type(t_vector)   ,intent(in) ,optional :: ana      ! analysis
  type(t_vector)   ,intent(in) ,optional :: ana_psas ! analysis in psas-space
  type(t_vector)   ,intent(in) ,optional :: varbc    ! variational bias corr.
  type(t_vector)   ,intent(in) ,optional :: e_fg     ! first guess error
  type(t_vector)   ,intent(in) ,optional :: s_fg     ! spread first guess
  type(t_vector)   ,intent(in) ,optional :: s_ana    ! spread analysis
  type(t_vector)   ,intent(in) ,optional :: m_fg     ! first guess mean
  type(t_vector)   ,intent(in) ,optional :: m_ana    ! analysis    mean
  type(t_vector)   ,intent(in) ,optional :: mo_fg    ! fg  mean in obs.space
  type(t_vector)   ,intent(in) ,optional :: mo_ana   ! an. mean in obs.space
  type(t_vector)   ,intent(in) ,optional :: tg_fg    ! Talagrand position fg
  type(t_vector)   ,intent(in) ,optional :: tg_ana   ! Talagrand position ana
  type(t_vector)   ,intent(in) ,optional :: en_fg (:)! first guess ensemble
  type(t_vector)   ,intent(in) ,optional :: en_ana(:)! analysis ensemble
  type(t_vector)   ,pointer    ,optional :: lin_ana(:)! Kalman gain
  character(len=*) ,intent(in) ,optional :: comment  ! step comment
  !==========================================================================
  ! write 3dvar feedback file
  !==========================================================================
    !------------------------
    ! loop bounds and indices
    !------------------------
    integer            :: nb      ! number of 'boxes'
    integer            :: nh      ! number of reports (header entries) on the PE
    integer            :: nf      ! number of files to write
    integer            :: ib      ! 'box' index variable
    integer            :: is      ! report index (within a box)
    integer            :: ih      ! report index (on a PE)
    integer            :: ot      ! observation type
    integer            :: satid   ! satellite id
    integer            :: grid_id ! instrument grid
    integer            :: if      ! file index variable
    integer            :: ip      ! processor index variable
    integer            :: iv      ! NetCDF variable index
    integer            :: i,j,ii  ! index variable
    integer            :: l       ! level index
    integer            :: k       ! ensemble index
    integer            :: ihbv    ! flag: write to header (1) or body (2)
    integer            :: i0, il  ! indices in 'body'
    integer            :: clen    ! character string length
    !--------------------------
    ! sizes of temporary arrays
    !--------------------------
    integer            :: mh, mb, lh, lb, lbh, ms, ls, lo(n_ot)
    !------------------
    ! NetCDF processing
    !------------------
    integer            :: ir      ! NetCDF return (status) variable
    character(PLEN)    :: file    ! file name
    character(PLEN)    :: path    ! full path name
    character(len=8)   :: lprefix ! file prefix
    integer            :: start   ! start index in file
    integer            :: cnt     ! number of words to write
    integer            :: ftype   ! file type (NetCDF or missing)
    integer            :: ltr     ! length of filename without trailing numbers
    !-------------------------------------------------
    ! to be written to NetCDF file (global attributes)
    !-------------------------------------------------
    real(sp)           :: dx, dy  ! model resolution (degree)
    character(len=22)  :: lcommnt ! comment (for history)
    integer            :: nx, ny  ! number of model gridpoints
    real(sp)           :: pole       (2) ! location of pole    (degree)
    real(sp)           :: lower_left (2) ! lower left  (lat,lon degree)
    real(sp)           :: upper_right(2) ! upper right (lat,lon degree)
    integer                       :: i1, in      ! indices to pack body
    type(t_index)    ,allocatable :: ix     (:)  ! local index array
    character(len=64)             :: opt
    character(len=8)              :: lmodel
    character(len=11)             :: lmethod     ! '3DVAR' or blank
    !------------------------------------
    ! temporary arrays (for send/receive)
    !------------------------------------
    real(sp)         ,allocatable :: b_real (:)  ! real numbers
    integer(i4)      ,allocatable :: b_int  (:)  ! integer numbers
    character(len=10),allocatable :: b_c10  (:)  ! character strings
    character(len=32),allocatable :: b_c32  (:)  ! character strings
    type(t_spot)     ,allocatable :: spt    (:)  ! header entries
    type(t_datum)    ,allocatable :: bdy    (:)  ! body entries
    type(t_radar)    ,pointer     :: radar  (:)  ! specific table entries
    real(sp)                      :: llev        ! level to compare
    logical                       :: lwsi(n_ot)  ! obstype uses WIGOS station id
    !-------------------------------
    ! temporary arrays (for gatherv)
    !-------------------------------
    real(sp)         ,allocatable :: r_real (:)             ! real numbers
    integer(i4)      ,allocatable :: r_int  (:)             ! integer numbers
    character(len=10),allocatable :: r_c10  (:)             ! character strings
    character(len=32),allocatable :: r_c32  (:)             ! character strings
    integer                       :: counts (0:dace% npe-1) ! Receive counts
    !------------------------
    ! default values in 3dvar
    !------------------------
    type(t_datum) :: default_body
    type(t_spot)  :: default_head
    !-------------------------------------------------------
    ! mask for observations to be written (depends on state)
    !-------------------------------------------------------
    type(t_bvector)      :: mask
    logical ,pointer     :: m (:)
    !-----------------------------------
    ! pointer to current observation box
    !-----------------------------------
    type(t_obs) ,pointer :: o
    !----------------
    ! pointer to obst
    !----------------
    type(t_obst) ,pointer :: obst
    !------------------
    ! radiance specific
    !------------------
    type(t_tovs)                 :: ttovs
    type(t_tovs_instr)           :: ti
    type(t_rad_set), pointer     :: rs => null()
    real(sp),        allocatable :: l2c_t(:), l2c_t_(:), sinfl_t(:) ! Required for load_tovs
    real(wp),        allocatable :: emiss_t(:)            ! ...
    integer,         allocatable :: tflag_t(:),ci(:)      ! ...
    real(sp),        allocatable :: l2c(:),sinfl(:),emiss(:),ins_tmp(:),orb_ph(:) ! Store stuff from t_tovs
    real(sp),        allocatable :: cld_frc(:) ! ...
    integer(i4),     allocatable :: tflag(:), nwc_flg(:), scanl(:), cld_flg(:)
    integer                      :: mx_body, iset
    logical                      :: l_tovs
    integer,         allocatable :: rad_opt_vars(:)

    character(len=300)    :: msg = ''

    call enter_function

    !=================================
    ! prepare file names, indices etc.
    !=================================
    lmodel  = rmodel
    lmethod = '3DVAR'
    if (present (method)) then
      select case (method)
      case ('LETKF')
         lmethod = 'LETKF'
      case ('PSAS+LETKF','ENVAR+LETKF')
         lmethod = '3DVAR LETKF'
      end select
    endif

    !--------------------------------------------------
    ! monitoring or analysis step ?
    ! set file prefix and comment for history attribute
    !--------------------------------------------------
    if (step >= 1 .and. step <= mstep) then
      lprefix = trim(pref_step(step))
      lcommnt = trim(comm_step(step)); if (present (comment)) lcommnt = comment
      if (step==4) n_source = 0
    else
      call finish('write_fdbk_3dv','invalid step')
    end if
    if (present (prefix)) lprefix = prefix(1:min(len(lprefix),len_trim(prefix)))

    !------------------------------------
    ! set up the file-name table
    ! (relate observation types to files)
    !------------------------------------
    call init_fdbk_tables
    call setup_filenames (fdbk, name, obs, step)

    !-------------------------------------------------------------------
    ! count number of report and body entries, allocate array components
    !-------------------------------------------------------------------
    nb = size (obs% o)                                  ! number of 'boxes'
    nf = size (fdbk% file)                              ! number of files
    nh = sum (obs%o% n_spot, mask=(obs%o% pe==dace% pe))! no.observations/PE
    allocate (fdbk% pelist (nf, 0:dace% npe-1))
    allocate (fdbk% index  (nh))
    call construct (mask, obs% oi)
    allocate(rad_opt_vars(nf)) ; rad_opt_vars = 0
    !---------------------------------------------------------------
    ! relate report indices to the positions in file header and body
    !---------------------------------------------------------------
    lwsi = .false.
    ih = 0
    do ib = 1, nb                            ! loop over 'boxes'
      if (obs% o(ib)% pe /= dace% pe) cycle  ! handled on this PE only
      o => obs%  o(ib)
      if (o% n_spot == 0) cycle
      m => mask% s(ib)% x                    ! check for state
      m =  o% body(1:o% n_obs)% use% state > STAT_DISMISS
      do is = 1, o% n_spot                   ! loop over reports
        ot      = o% spot(is)% hd% obstype   ! observation type
        satid   = o% spot(is)% hd% satid     ! satellite ID
        grid_id = o% spot(is)% hd% grid_id   ! instr. grid

        if = -1
        do i = 1, size(fdbk% file_obst)
          obst => fdbk% file_obst(i)
          if (ot      /= obst% obstype                     ) cycle
          if (satid   /= obst% satid .and. obst% satid >= 0) cycle
          if (grid_id /= obst% grid  .and. obst% grid  >= 0) cycle
          if = i
        enddo
        !if = fdbk% file_obst(ot)% file       ! file index

        if (o% spot(is)% use% state <= STAT_DISMISS) then
           if = 0
        else
           if (if <= 0) then
              write(msg,'("Did not find file for spot ",I9," satid=",I3.3," grid=",I3.3)') &
                   o% spot(is)% hd% id, satid, grid_id
              call finish('write_fdbk_3dvar', trim(msg))
           end if
        end if
        if (if > 0) then                     ! valid file indices only
          fdbk% file(if)%  n_hdr      =  fdbk% file(if)% n_hdr + 1
          ih                          =  ih                    + 1
          i0                          =  o% spot(is)% o% i
          il                          =  o% spot(is)% o% n
          m                           => mask% s(ib)% x (i0+1:i0+il)
          fdbk% index(ih)% file       =  if
          fdbk% index(ih)% i_hdr      =  ih
          fdbk% index(ih)% o_body     =  fdbk% file(if)% n_body
          fdbk% index(ih)% n_body     =  count (m)
          fdbk% index(ih)% o_spec     =  fdbk% file(if)% n_spec(ot)
          fdbk% index(ih)% n_spec     =  o% spot(is)% s% n
          fdbk% index(ih)% ib         =  ib
          fdbk% index(ih)% is         =  is
          fdbk% file(if)%  n_body     =  fdbk% file (if)% n_body &
                                      +  fdbk% index(ih)% n_body
          fdbk% file(if)%  n_spec(ot) =  fdbk% file(if)%  n_spec(ot) &
                                      +  fdbk% index(ih)% n_spec
          !---------------------------------
          ! keep location in monitoring file
          !---------------------------------
          if (step == 1 .or. step == 4) then
            o% spot(is)% hd% mon_file = if
            o% spot(is)% hd% mon_rec  = fdbk% file(if)% n_hdr
            j = 0
            do i = 1, il
              if (m(i)) then
                j = j + 1
                o% body (i0+i)% mon_pos = j
              else
                o% body (i0+i)% mon_pos = -1
              endif
            end do
          endif
          !---------------------------------------------------
          ! Check if at least one station with WSI for obstype
          !---------------------------------------------------
          if (o% spot(is)% wsi% valid) lwsi(ot) = .true.
          !------------------------------------------------------------------------
          ! Check whether special radiance specific variables are required for file
          !------------------------------------------------------------------------
          if (ot == OT_RAD) then
            iset = set_indx(rad_set, satid=satid, grid=grid_id)
            if (iset > 0 .and. iset <= n_set) then
              rad_opt_vars(if) = ior(rad_opt_vars(if), rad_set(iset)%gopts%opt_vars)
            else
              ! call finish('write_fdbk_3dv','failed to find rad_set')
            end if
          end if
        else
          if (step == 1 .or. step == 4) o% spot(is)% hd% mon_file = 0
        endif
      end do
    end do
    fdbk% index(ih+1:)% file = 0
    !--------------------------------------------------------
    ! broadcast the header and body size required for this PE
    !--------------------------------------------------------
    fdbk% pelist (:,dace% pe) = fdbk% file
    do ip = 0, dace% npe-1
      call p_bcast_t_file (fdbk% pelist(:,ip), ip)
    end do
    !----------------------------------------
    ! sum of header and body sizes on all PEs
    !----------------------------------------
    do if = 1, nf
      fdbk% file(if)%   n_hdr      = sum (fdbk% pelist(if,:)% n_hdr)
      fdbk% file(if)%   n_body     = sum (fdbk% pelist(if,:)% n_body)
      do ot = 1,n_ot
        fdbk% file(if)% n_spec(ot) = sum (fdbk% pelist(if,:)% n_spec(ot))
      end do
      do ip = 1, dace% npe-1
        fdbk% pelist(:,ip)%   o_hdr      = fdbk% pelist(:,ip-1)% o_hdr  &
                                         + fdbk% pelist(:,ip-1)% n_hdr
        fdbk% pelist(:,ip)%   o_body     = fdbk% pelist(:,ip-1)% o_body &
                                         + fdbk% pelist(:,ip-1)% n_body
        do ot = 1,n_ot
          fdbk% pelist(:,ip)% o_spec(ot) = fdbk% pelist(:,ip-1)% o_spec(ot) &
                                         + fdbk% pelist(:,ip-1)% n_spec(ot)
        end do
      end do
    end do
    !------------------------------------
    ! correct location in monitoring file
    !------------------------------------
    if (step == 1 .or. step == 4) then
      do ib = 1, nb
        if (obs%o (ib)% pe /= dace% pe) cycle  ! handled on this PE only
        o => obs% o(ib)
        do is = 1, o% n_spot
          if = o% spot(is)% hd% mon_file
          if (if > 0) then
            o% spot(is)% hd% mon_rec = o% spot(is)% hd% mon_rec       &
                                     + fdbk% pelist(if,dace% pe)% o_hdr
            !-----------------------------------------
            ! for verification files set source/record
            !-----------------------------------------
            if (step == 4) then
              o% spot(is)% hd% source = o% spot(is)% hd% mon_file
              o% spot(is)% hd% record = o% spot(is)% hd% mon_rec
            endif

          endif
        end do
      end do
    endif
    !------------------------------------
    ! broadcast info on special variables
    !------------------------------------
    rad_opt_vars = p_ior(rad_opt_vars)
    !-------------------------------------------
    ! allocate buffers with the appropriate size
    !-------------------------------------------
    lh = maxval (fdbk% pelist(:,dace% pe)% n_hdr)
    mh = maxval (fdbk% pelist(:,:       )% n_hdr)
    lb = maxval (fdbk% pelist(:,dace% pe)% n_body)
    mb = maxval (fdbk% pelist(:,:       )% n_body)
    ls = 0
    ms = 0
    do ot = 1, n_ot
      ls = max(ls, maxval (fdbk% pelist(:,dace% pe)% n_spec(ot)))
      ms = max(ms, maxval (fdbk% pelist(:,:       )% n_spec(ot)))
    end do
    lwsi = p_or (lwsi)
    if (dace% lpio) then
      allocate (b_int  (max(mb,mh,ms)))
      allocate (b_real (max(mb,mh,ms)))
      allocate (b_c10  (max(mb,mh,ms)))
      if (any (lwsi))                 &
        allocate (b_c32(max(mb,mh,ms)))
    else
      allocate (b_int  (max(lb,lh,ls)))
      allocate (b_real (max(lb,lh,ls)))
      allocate (b_c10  (max(lb,lh,ls)))
      if (any (lwsi))                 &
        allocate (b_c32(max(lb,lh,ls)))
    endif
    allocate (ix  (lh))
    allocate (spt (lh))
    allocate (bdy (lb))

    !=============================
    ! now actually write the files
    !=============================

    !----------------
    ! loop over files
    !----------------
    do if = 1, nf
      file = trim(lprefix)//trim(fdbk% file(if)% name)//'.nc'
      !-----------------------------------
      ! bookkeeping for verification files
      !-----------------------------------
      if (step == 4) then
        ftype                                  = FT_FEEDBACK
        if (fdbk% file(if)% n_body == 0) ftype = FT_MISSING
        call add_source (aux, file, ftype, complete=.true.)
      endif
      !-------------------------------
      ! loop over non-empty files only
      !-------------------------------
      if (fdbk% file(if)% n_hdr > 0 .and. fdbk% file(if)% n_body > 0) then
        !-----------------------
        ! create the NetCDF file
        !-----------------------
        if (dace% lpio) then
          !-----------------------------------
          ! set up derived type, set file name
          !-----------------------------------
          call setup_fdbk (fdbk% fb% nc)
          path = path_file (aux, file)
          fdbk% file(if)% path = path
          !---------------------------------
          ! derive model resolution and size
          !---------------------------------
          pole        = 0._sp
          lower_left  = 0._sp
          upper_right = 0._sp
          select case (grid% gridtype)
          case (WMO6_LATLON, WMO6_ROTLL, WMO6_GAUSSIAN)
            nx          = grid% nx
            ny          = grid% ny
            dx          = grid% di
            dy          = grid% dj
            lower_left  = [grid% dlat(1 ), grid% dlon(1 )]
            upper_right = [grid% dlat(ny), grid% dlon(nx)]
            if (grid% gridtype == WMO6_ROTLL)            &
            pole        = [grid% dlatr   , grid% dlonr   ]
          case (DWD6_ICOSAHEDRON)
            nx = grid% ni
            ny = nx
            dx = 63.486_sp/nx
            dy = dx
          case (DWD6_ICON)
            nx = grid% ni
            ny = nx
            dx = 45.416_sp/nx
            dy = dx
          case default
            write(0,*) "gridtype =", grid% gridtype
            call finish('write_fdbk_3dv','unknown gridtype')
          end select
          !--------------------------------------------------------------
          ! create file (global attributes, dimensions, define variables)
          !--------------------------------------------------------------
          ! check if additional output variables are allowed
          call check_all_fdbk_addvar(fdbk_addvar)
          opt=" "//fdbk_addvar

          if (lwsi(fdbk% file_obst(if)% obstype)) then
             opt = " WSI" // trim (opt)
          end if

          if (fdbk% file_obst(if)% obstype == OT_RAD) then
            do i = 0, n_optv-1
              if (btest(rad_opt_vars(if), i)) opt = trim(opt)//' '//c_optv(i)
            end do
          end if
          do ltr = len_trim(fdbk% file(if)% name), 1, -1
            i = iachar(fdbk% file(if)% name(ltr:ltr))
            if (i < 48 .or. i > 57) exit   ! not a digit (0-9)
          end do
          call create_fdbk  (                         &!
                     fdbk% fb,                        &! feedbackfile meta data
                     path,                            &! path
                     trim(lmodel),                    &! model
                     var3d_version(),                 &! model version
                     'German Weather Service',        &! institution
                     fdbk% file(if)% n_hdr,           &! d_hdr
                     fdbk% file(if)% n_body,          &! d_body
                     iyyyymmdd(reftime),              &! reference date
                     ihhmm    (reftime),              &! reference time
                     0,0,                             &! start,stop of verific.
                     (/dx,dy/),                       &! resolution
                     (/nx,ny,grid%nz/),               &! domain
                     lcommnt,                         &! comment  for history
                     cyyyymmddhhmm(reftime),          &! time     for history
            runtime= cyyyymmddhhmm(run_time),         &! run time for history
               pole= pole,                            &!
         lower_left= lower_left,                      &!
        upper_right= upper_right,                     &!
            n_radar= fdbk% file(if)% n_spec(OT_RADAR),&!
                opt= trim(lmodel)//' '                &! flag
                     //trim(lmethod)//' '             &!  for
                     //fdbk% file(if)% name(1:ltr)    &!   optional
                     //opt)                            !   variables
          if (fdbk% fb%nc% error /= 0) &
            call finish('write_fdbk_3dv','cannot create: '//trim(path))
          !-------------------------------------------
          ! re-write size of header/body actually used
          !-------------------------------------------
          fdbk% fb% n_hdr   = fdbk% file(if)% n_hdr    ! 3dvar always uses
          fdbk% fb% n_body  = fdbk% file(if)% n_body   ! all allocated memory
          fdbk% fb% n_radar = fdbk% file(if)% n_spec(OT_RADAR)
          call write_global_attributes (fdbk% fb)
          ir = nf90_enddef (fdbk% fb% nc% ncid)
        endif
        !--------------------------
        ! prepare local index array
        !--------------------------
        lh = fdbk% pelist(if, dace% pe)% n_hdr
        lb = fdbk% pelist(if, dace% pe)% n_body
        lo = fdbk% pelist(if, dace% pe)% n_spec (:)
        if  (lh>0) &                  !+++ fix for pack with count(mask)==0 +++!
        ix  (1:lh) = pack (fdbk% index, mask = fdbk% index% file == if)
        !------------------
        ! pack header, body
        !------------------
        do i =  1, lh
          ib =  ix(i)% ib
          is =  ix(i)% is
          i1 =  ix(i)% o_body + 1
          in =  ix(i)% n_body - 1 + i1
          o  => obs% o(ib)
          i0 =  o% spot(is)% o% i + 1
          il =  o% spot(is)% o% n - 1 + i0
          spt (i) = o% spot(is)
          if (i1<=in) &               !+++ fix for pack with count(mask)==0 +++!
          bdy (i1:in) = pack (o% body (i0:il), mask% s(ib)% x(i0:il))
        end do
        !---------------------------------
        ! pack observation specific tables
        !---------------------------------
        call radar_tables (radar, spt(1:lh), obs% o, ix(1:lh)% o_spec, ix(1:lh)% ib, &
                           fdbk% pelist(if, dace% pe)% n_spec(OT_RADAR), mask)
        ! Unpack stuff stored in t_tovs for radiances
        l_tovs = any(spt(1:lh)% hd% obstype == OT_RAD .and. spt(1:lh)%p%n > 0)
        l_tovs = p_or (l_tovs)
        if (l_tovs) then
          mx_body = maxval(ix(1:lh)% n_body)
          if (lh > 0) allocate(l2c_t(mx_body), l2c_t_(mx_body), emiss_t(mx_body),&
                               tflag_t(mx_body), ci(mx_body), sinfl_t(mx_body))
          allocate(emiss(lb), tflag(lb), sinfl(lb), scanl(lh))
          if (btest(rad_opt_vars(if),OPTV_L2C    )) allocate(l2c    (lb))
          if (btest(rad_opt_vars(if),OPTV_NWC_FLG)) allocate(nwc_flg(lh))
          if (btest(rad_opt_vars(if),OPTV_ORB_PH )) allocate(orb_ph (lh))
          if (btest(rad_opt_vars(if),OPTV_INS_TMP)) allocate(ins_tmp(lh))
          if (btest(rad_opt_vars(if),OPTV_CLD_FRC)) allocate(cld_frc(lh))
          if (btest(rad_opt_vars(if),OPTV_CLD_FLG)) allocate(cld_flg(lh))
          do i =  1, lh
            i1 =  ix(i)% o_body + 1
            in =  ix(i)% n_body - 1 + i1
            if (spt(i)%hd%obstype == OT_RAD .and. spt(i)%p%n > 0) then
              ib =  ix(i)% ib
              is =  ix(i)% is
              o  => obs% o(ib)
              i0 =  spt(i)% o% i + 1
              il =  spt(i)% o% i + spt(i)% o% n
              if (n_set > 0) then
                call load(obs% o(ib), spt(i), tovs=ttovs, rs=rs, ti=ti, tovs_io=0, &
                     l2c=l2c_t, emis=emiss_t, sinfl=sinfl_t, flag=tflag_t, ci=ci)
              else ! (e.g. pure LETKF)
                call load(obs% o(ib), spt(i), tovs=ttovs, tovs_io=0, &
                     l2c=l2c_t, emis=emiss_t, sinfl=sinfl_t, flag=tflag_t, ci=ci)
                rs => null()
              end if
              ! ci is not required here, however the ti optional argument will cause the loading
              ! of ci. It is more efficient to supply ci as an external array than to
              ! allocate/deallocate t_tovs%ci for each individual spot.
              if (ttovs%nchan /= spt(i)%o%n) call finish('write_fdbk_3dv',&
                   'Inconcsistent t_tovs%nchan and spt%o%n.')
              ! if (count(mask%s(ib)%x(i0:il)) /= ix(i)%n_body) call finish('write_fdbk_3dv',&
              !      'Inconsistent mask and n_body.')
              j = 0
              if (btest(rad_opt_vars(if),OPTV_L2C)) then
                if (associated(ttovs%l2c) .and. associated(rs)) then
                  do ii = 1, ti%n_instr
                    if (rs% iopts(ti%ii(ii))% l2c_type > 0) then
                      l2c_t_(ti%o_ch_i(ii)+1:ti%o_ch_i(ii)+ti%n_ch_i(ii)) = &
                           l2c_t(j+1:j+ti%n_ch_i(ii))
                      j = j + ti%n_ch_i(ii)
                    else
                      l2c_t_(ti%o_ch_i(ii)+1:ti%o_ch_i(ii)+ti%n_ch_i(ii)) = -huge(b_real(1))
                    end if
                  end do
                else
                  l2c_t_(1:ttovs%nchan) = -huge(b_real(1))
                end if
                l2c(i1:in) = pack(l2c_t_ (1:ttovs%nchan),mask=mask%s(ib)%x(i0:il))
              end if
              if (btest(rad_opt_vars(if),OPTV_NWC_FLG)) nwc_flg(i) = int(ttovs%nwc_flg,kind=i4)
              if (btest(rad_opt_vars(if),OPTV_ORB_PH )) orb_ph (i) =     ttovs%orb_phase
              if (btest(rad_opt_vars(if),OPTV_INS_TMP)) ins_tmp(i) =     ttovs%instr_temp
              if (btest(rad_opt_vars(if),OPTV_CLD_FRC)) cld_frc(i) =     ttovs%cloud_imag
              if (btest(rad_opt_vars(if),OPTV_CLD_FLG)) cld_flg(i) =     ttovs%cld_flg
              scanl  (i)   = int(ttovs%scanl,  kind=i4)
              if (associated(ttovs%emis)) then
                emiss(i1:in) = pack(real(emiss_t(1:ttovs%nchan),kind=sp),mask=mask%s(ib)%x(i0:il))
              else
                emiss(i1:in) = -huge(b_real(1))
              end if
              if (associated(ttovs%sinfl)) then
                sinfl(i1:in) = pack(sinfl_t(1:ttovs%nchan),mask=mask%s(ib)%x(i0:il))
              else
                sinfl(i1:in) = -huge(b_real(1))
              end if
              if (associated(ttovs%flag)) then
                tflag(i1:in) = pack(int(tflag_t(1:ttovs%nchan),kind=i4),mask=mask%s(ib)%x(i0:il))
              else
                tflag(i1:in) = 0
              end if
            else
              scanl    (i) = -huge(b_int(1))
              emiss(i1:in) = -huge(b_real(1))
              sinfl(i1:in) = -huge(b_real(1))
              tflag(i1:in) = 0
              if (btest(rad_opt_vars(if),OPTV_L2C    )) l2c  (i1:in) = -huge(b_real(1))
              if (btest(rad_opt_vars(if),OPTV_NWC_FLG)) nwc_flg  (i) = 0
              if (btest(rad_opt_vars(if),OPTV_ORB_PH )) orb_ph   (i) = -huge(b_real(1))
              if (btest(rad_opt_vars(if),OPTV_INS_TMP)) ins_tmp  (i) = -huge(b_real(1))
              if (btest(rad_opt_vars(if),OPTV_CLD_FRC)) cld_frc  (i) = -huge(b_real(1))
              if (btest(rad_opt_vars(if),OPTV_CLD_FLG)) cld_flg  (i) = 0
            end if
          end do
        end if
        !--------------------
        ! loop over variables
        !--------------------
        call p_bcast (fdbk% fb% nc% nvar, dace% pio)
        if (.not.dace% lpio) allocate (fdbk% fb% nc% vars(fdbk% fb% nc% nvar))
        call p_bcast_t_var (fdbk% fb% nc% vars(1:fdbk% fb% nc% nvar), dace% pio)

        do iv = 1, fdbk% fb% nc% nvar
          if (.not. fdbk% fb% nc% vars(iv)% opt_used) cycle
          !-------------------------------
          ! fill variable into send-buffer
          !-------------------------------
          ihbv = 0
          clen = 0
          select case (fdbk% fb% nc% vars(iv)% name)
          !------------------
          ! integer variables
          !------------------
          case ('i_body')
            b_int(1:lh) = ix(1:lh)% o_body + fdbk% pelist(if, dace% pe)% o_body + 1
            ihbv = -1
          case ('l_body')
            b_int(1:lh) = ix(1:lh)% n_body
            ihbv = -1
          case ('n_level')
            do i = 1, lh
              if (spt(i)% col% nlev == 1) then
                b_int(i) = 1
              else if (spt(i)% hd% obstype == OT_GPSRO) then
                b_int(i) = ix(i)% n_body
              else
                i1 =  ix(i)% o_body + 1
                in =  ix(i)% n_body - 1 + i1
                llev = huge(1._sp)
                l    = 0
                do j = i1, in
                  if (bdy(j)% plev /= llev) then
                    l    = l + 1
                    llev = bdy(j)% plev
                  endif
                end do
                b_int(i) = l
              endif
            end do
            ihbv = -1
          case ('data_category')
            b_int(1:lh) = spt(1:lh)% hd% buf_type
            ihbv = -1
          case ('sub_category')
            b_int(1:lh) = spt(1:lh)% hd% buf_subtype
            ihbv = -1
          case ('obstype')
            b_int(1:lh) = spt(1:lh)% hd% obstype
            ihbv = -1
          case ('codetype')
            b_int(1:lh) = spt(1:lh)% hd% codetype
            ihbv = -1
          case ('ident')
            b_int(1:lh) = spt(1:lh)% ident
            ihbv = -1
          case ('statid')
            b_c10(1:lh) = spt(1:lh)% statid
            ihbv = -1
            clen = 10
          case ('wsi')
            do i = 1, lh
               if (spt(i)% wsi% valid) then
                  call wsi_to_text (spt(i)% wsi, b_c32(i))
               else
                  b_c32(i) = ""
               end if
            end do
            ihbv = -1
            clen = 32
          case ('lat')
            b_real(1:lh) = spt(1:lh)% col% c% dlat
            ihbv = -1
          case ('lon')
            b_real(1:lh) = spt(1:lh)% col% c% dlon
            ihbv = -1
          case ('time')
            b_int(1:lh) = nint (minutes (spt(1:lh)% actual_time - reftime))
            ihbv = -1
          case('time_nomi')
            b_int(1:lh) = nint (minutes (spt(1:lh)% hd%    time - reftime))
            ihbv = -1
          case('time_dbase')
            where (spt(1:lh)% hd% db_time /= zero_time)
              b_int(1:lh) = nint (minutes (spt(1:lh)% hd% db_time - reftime))
            elsewhere
              b_int(1:lh) = -huge(b_int(1))
            endwhere
            ! Catch case where newer test data are used from the database
            where (b_int(1:lh) > huge(0_i2)) b_int(1:lh) = huge(0_i2)
            ihbv = -1
          case ('z_station')
            b_int(1:lh) = nint (spt(1:lh)% z)
            ihbv = -1
          case ('z_modsurf')
            b_int(1:lh) = nint (spt(1:lh)% gp_bg / gacc)
            ihbv = -1
          case ('r_state')
            b_int(1:lh) = spt(1:lh)% use% state
            call change_code (b_int(1:lh), stats% code)
            ihbv = -1
          case ('r_flags')
            b_int(1:lh) = spt(1:lh)% use% flags
            call change_bits (b_int(1:lh), chk% code)
            ihbv = -1
          case ('r_check')
            b_int(1:lh) = spt(1:lh)% use% check
            call change_code (b_int(1:lh), chk% code)
            ihbv = -1
          case ('sta_corr')
            b_int(1:lh) = min(spt(1:lh)% corme, 1)
            ihbv = -1
          case ('obs_id')
            b_int(1:lh) = spt(1:lh)% hd% id
            ihbv = -1
          case ('source')
            b_int(1:lh) = spt(1:lh)% hd% source
            ihbv = -1
          case ('record')
            b_int(1:lh) = spt(1:lh)% hd% record
            ihbv = -1
          case ('scanline')
            b_int(1:lh) = scanl(1:lh)
            ihbv = -1
          case ('subset')
            b_int(1:lh) = spt(1:lh)% hd% subset
            ihbv = -1
          case ('dbkz')
            b_int(1:lh) = spt(1:lh)% hd% dbkz
            ihbv = -1
          case ('index_x')
            b_int(1:lh) = spt(1:lh)% col% h% ijdp(1)
            ihbv = -1
          case ('index_y')
            b_int(1:lh) = spt(1:lh)% col% h% ijdp(2)
            ihbv = -1
          case ('index_d')
            b_int(1:lh) = spt(1:lh)% col% h% ijdp(3)
            ihbv = -1
          case ('mdlsfc')
            b_int(1:lh) = spt(1:lh)% mdlsfc
            ihbv = -1
          case ('vnyquist')
            where (spt(1:lh)% hd% obstype == OT_RADAR)
              b_real(1:lh) = spt(1:lh)% params(1)
            elsewhere
              b_real(1:lh) = -huge (b_real(1))
            endwhere
            ihbv = -1
          case ('fr_land')
            where (spt(1:lh)% sl_bg >= 0._wp)
              b_real(1:lh) = spt(1:lh)% sl_bg
            elsewhere
              b_real(1:lh) = -huge (b_real(1))
            end where
            ihbv = -1
          case ('sso_stdh')
            where (spt(1:lh)% ssd_bg >= 0._wp)
              b_real(1:lh) = spt(1:lh)% ssd_bg
            elsewhere
              b_real(1:lh) = -huge (b_real(1))
            end where
            ihbv = -1
          case ('nwc_flag')
            b_int(1:lh) = nwc_flg(1:lh)
            ihbv = -1
          case ('orbit_phase')
            b_real(1:lh) = orb_ph(1:lh)
            ihbv = -1
          case ('instr_temp')
            b_real(1:lh) = ins_tmp(1:lh)
            ihbv = -1
          case ('cloud_frac')
            b_real(1:lh) = cld_frc(1:lh)
            ihbv = -1
          case ('cloud_flag')
            b_int(1:lh) = cld_flg(1:lh)
            ihbv = -1
          case ('varno')
            do i =  1, lh
              ib =  ix(i)% ib
              is =  ix(i)% is
              i1 =  ix(i)% o_body + 1
              in =  ix(i)% n_body - 1 + i1
              o  => obs%  o(ib)
              i0 =  o% spot(is)% o% i + 1
              il =  o% spot(is)% o% n - 1 + i0
              m  => mask% s(ib)% x (i0:il)
              b_int(i1:in) = pack (          o% varno (i0:il) , m)
!             !++++++++++++++++++++++++++++++++++++++++++++++++++
!             ! temporarily fix inconsistent labeling of channels
!             !++++++++++++++++++++++++++++++++++++++++++++++++++
!             select case (o% spot(is)% hd% obstype)
!             case (OT_RAD)
!               b_int(i1:in) = VN_RAWBT
!             end select
!             !++++++++++++++++++++++++++++++++++++++++++++++++++
            end do
            ihbv = -2
          case ('obs')
            b_real(1:lb) = bdy(1:lb)% o
            ihbv = -2
          case ('bcor')
            b_real(1:lb) = bdy(1:lb)% bc
            ihbv = -2
          case ('plevel')
            b_real(1:lb) = bdy(1:lb)% plev
            where (b_real(1:lb) < 0.) &
              b_real(1:lb) = fdbk% fb% nc% vars(iv)% rinvalid
            ihbv = -2
          case ('level')
            do i = 1, lh
              ib = ix(i)% ib
              is = ix(i)% is
              i1 = ix(i)% o_body + 1
              in = ix(i)% n_body - 1 + i1
              o  => obs% o(ib)
              i0 =  o% spot(is)% o% i + 1
              il =  o% spot(is)% o% n - 1 + i0
              m  => mask% s(ib)% x (i0:il)
              b_real(i1:in) = pack (o% olev (i0:il), m)
!             !++++++++++++++++++++++++++++++++++++++++++++++++++
!             ! temporarily fix inconsistent labeling of channels
!             !++++++++++++++++++++++++++++++++++++++++++++++++++
!             select case (o% spot(is)% hd% obstype)
!             case (OT_RAD)
!!              !+++ does not work on NEC
!!              where (.not.hss_instr(int(bdy(i1:in)% lev_sig), wmo=.true.)) &! not hyperspectral sounder
!!                  b_real(i1:in) = mod (b_real(i1:in),100._sp)
!               !+++ works on NEC
!               do j=i1,in
!                 if (.not.hss_instr(int(bdy(j)% lev_sig), wmo=.true.)) &! not hyperspectral sounder
!                     b_real(j) = mod (b_real(j),100._sp)
!               end do
!             end select
!             !++++++++++++++++++++++++++++++++++++++++++++++++++
            end do
            ihbv = -2

          case ('level_typ')
            !--------------------------
            ! level_typ: 3dvar settings
            !--------------------------
            do i = 1, lh
              ib = ix(i)% ib
              is = ix(i)% is
              i1 = ix(i)% o_body + 1
              in = ix(i)% n_body - 1 + i1
              select case (obs% o(ib)% spot(is)% hd% obstype)
              case (OT_RAD)
                b_int(i1:in) = VN_NUM
              case (OT_GPSRO)
                b_int(i1:in) = VN_IMPPAR
              case default
                b_int(i1:in) = VN_P
              end select
            end do
            !------------------------------------------
            ! level_typ: taken from fof files for KENDA
            !------------------------------------------
            where (bdy(1:lb)% lev_typ >= 0) &
              b_int(1:lb) = bdy(1:lb)% lev_typ
            ihbv = -2
          case ('level_sig')
            b_int(1:lb) = bdy(1:lb)% lev_sig
            where (b_int(1:lb) < 0)             &
                   b_int(1:lb) = -huge(b_int(1))
            ihbv = -2
          case ('state')
            b_int(1:lb) = bdy(1:lb)% use% state
            call change_code (b_int(1:lb), stats% code)
            ihbv = -2
          case ('flags')
            b_int(1:lb) = bdy(1:lb)% use% flags
            call change_bits (b_int(1:lb), chk% code)
            ihbv = -2
          case ('check')
            b_int(1:lb) = bdy(1:lb)% use% check
            call change_code (b_int(1:lb), chk% code)
            ihbv = -2
          case ('qual')                  ! SATOB, GPSRO: per cent confidence
            b_int(1:lb) = bdy(1:lb)% pcc ! AIREP       : roll angle
            where (b_int(1:lb) == default_body% pcc) &
                   b_int(1:lb) =  -huge(b_int(1))
            ihbv = -2
          case ('l2c')                   ! IASI level assignment
            do i = 1, lh
              select case (spt(i)% hd% obstype)
              case (OT_RAD)
                i1 = ix(i)% o_body + 1
                in = ix(i)% n_body - 1 + i1
                b_real(i1:in) = l2c(i1:in)
              case (OT_RADAR)
!                 i1 = ix(i)% o_body + 1
!                 in = ix(i)% n_body - 1 + i1
!                 b_real(i1:in) = bdy(i1:in)% l2c
!                 where (b_real(i1:in) == default_body% l2c) &
!                        b_real(i1:in) =  -huge(b_real(1))
                ii = ix(i)% o_body + 1
                do j = 1, obs% o(ib)% spot(is)% s% n
                  do k = 1, radar(j)% radar_nbody
                    if (bdy(ii)% spec_index /= default_body% spec_index) then
                      b_real(ii) = radar(j)% radar_range_start + radar(j)% radar_drange * &
                                   (bdy(ii)% spec_index - 1)
                    else
                      b_real(ii) = -huge(b_real(1))
                    end if
                  end do
                end do
              case default
                b_real(i1:in) = -huge(b_real(1))
              end select
            end do
            ihbv = -2
          case ('emissiv')                   ! emissivities
            b_real(1:lb) = emiss(1:lb)
            ihbv = -2
          case ('surf_infl') ! surface influence
            b_real(1:lb) = sinfl(1:lb)
            ihbv = -2
          case ('tovs_flag')
            b_int(1:lb) = tflag(1:lb)
            ihbv = -2
          case ('plev_width')            ! IASI pressure level width
            b_real(1:lb) = bdy(1:lb)% plev_width
            where (b_real(1:lb) == default_body% plev_width) &
                   b_real(1:lb) =  -huge(b_real(1))
            ihbv = -2
          case ('accuracy')
            b_real(1:lb) = bdy(1:lb)% ac ! PILOT (wind profiler)
            where (b_real(1:lb) == default_body% ac) &
                   b_real(1:lb) =  -huge(b_real(1))
            ihbv = -2
          case ('dlat')
            b_real(1:lb) = bdy(1:lb)% lat ! TEMP/GPSRO
            where (b_real(1:lb) == default_body% lat) &
                   b_real(1:lb) =  -huge(b_real(1))
            ihbv = -2
          case ('dlon')
            b_real(1:lb) = bdy(1:lb)% lon ! TEMP/GPSRO
            where (b_real(1:lb) == default_body% lon) &
                   b_real(1:lb) =  -huge(b_real(1))
            ihbv = -2
          case ('azimuth')
            b_real(1:lb) = bdy(1:lb)% obs_par(1) ! GPSGB/GPSRO azimuth
            where (b_real(1:lb) == default_body% obs_par(1)) &
                   b_real(1:lb) =  -huge(b_real(1))
            ihbv = -2
          case ('obs_par_1')
            b_real(1:lb) = bdy(1:lb)% obs_par(1) ! RAD surface influence
            where (b_real(1:lb) == default_body% obs_par(1)) &
                   b_real(1:lb) =  -huge(b_real(1))
            ihbv = -2
          case ('obs_par_2')
            b_real(1:lb) = bdy(1:lb)% obs_par(2) ! currently not used
            where (b_real(1:lb) == default_body% obs_par(2)) &
                   b_real(1:lb) =  -huge(b_real(1))
            ihbv = -2
          case ('v_loc')
            b_real(1:lb) = bdy(1:lb)% set% v_loc ! vert.loc.radius
            where (b_real(1:lb) == default_body% set% v_loc) &
                   b_real(1:lb) =  -huge(b_real(1))
            ihbv = -2
          case ('h_loc')
            b_real(1:lb) = bdy(1:lb)% set% h_loc ! hor.loc.radius
            where (b_real(1:lb) == default_body% set% h_loc) &
                   b_real(1:lb) =  -huge(b_real(1))
            ihbv = -2
          case ('e_o')
            if (present(e_o)) then
              !-------------------------------------
              ! possibly modified obs.error in 3dvar
              !-------------------------------------
              do i =  1, lh
                ib =  ix(i)% ib
                is =  ix(i)% is
                i1 =  ix(i)% o_body + 1
                in =  ix(i)% n_body - 1 + i1
                o  => obs% o(ib)
                i0 =  o% spot(is)% o% i + 1
                il =  o% spot(is)% o% n - 1 + i0
                m  => mask% s(ib)% x (i0:il)
                b_real(i1:in) = pack (e_o% s(ib)% x(i0:il), m)
              end do
              ihbv = -2
            else
              !-------------------
              ! obs.error in LETKF
              !-------------------
              b_real(1:lb) = bdy(1:lb)% eo
              where (b_real(1:lb) == default_body% eo) &
                     b_real(1:lb) =  -huge(b_real(1))
              ihbv = -2
            endif
          case ('w_qc')
            if (.not. present(w_qc)) cycle
            b_real(1:lb) = 0._sp
            if (associated(w_qc% s)) then
              do i =  1, lh
                ib =  ix(i)% ib
                is =  ix(i)% is
                i1 =  ix(i)% o_body + 1
                in =  ix(i)% n_body - 1 + i1
                o  => obs% o(ib)
                i0 =  o% spot(is)% o% i + 1
                il =  o% spot(is)% o% n - 1 + i0
                m  => mask% s(ib)% x (i0:il)
                b_real(i1:in) = pack (w_qc% s(ib)% x(i0:il), m)
              end do
            end if
            ihbv = -2
          case ('spec_index')
            b_int(1:lb) = bdy(1:lb)% spec_index  ! obstype specific table entry
            where (b_int(1:lb) == default_body% spec_index) &
                   b_int(1:lb) =  -huge(b_int(1))
            ihbv = -2
          case ('center')
            b_int(1:lh) = spt(1:lh)% hd% center
            ihbv = -1
          case ('sub_center')
            b_int(1:lh) = spt(1:lh)% hd% subcenter
            ihbv = -1
          case ('surftype')                    ! surface type, 1dvar convention
            b_int(1:lh) = spt(1:lh)% stlsf     ! +++ CHECK MEANING OF AMV FLAGS
            ihbv = -1
          case ('soiltype')
            b_int(1:lh) = spt(1:lh)% soiltype  ! soiltype
            ihbv = -1
          case ('phase','fov')                 ! AIREP:phase of aircraft flight
                                               ! RAD: fov (Field of View index)
                                               ! GPSRO: PCD flags
            b_int(1:lh) = spt(1:lh)% phase
            ihbv = -1
          case ('instype')                     ! TEMP: station or sat.instr. type, Table C-2, C-8
            b_int(1:lh) = spt(1:lh)% sttyp     ! RAD: interpolated/aggregated to instrument
            where (b_int(1:lh) == default_head% sttyp) &
                   b_int(1:lh) =  -huge(b_int(1))
            ihbv = -1
          case ('rad_corr')                    ! TEMP
            b_int(1:lh) = spt(1:lh)% stret     ! solar and infrared radiation correction
            where (b_int(1:lh) == default_head% stret) &
                   b_int(1:lh) =  -huge(b_int(1))
            ihbv = -1
          case ('sat_class')                   ! GPSRO
            b_int(1:lh) = spt(1:lh)% tracking  ! (GNSS) satellite class
            where (b_int(1:lh) == default_head% tracking) &
                   b_int(1:lh) =  -huge(b_int(1))
            ihbv = -1
          case ('prn')                         ! GPSRO
            b_int(1:lh) = spt(1:lh)% sender_id ! (GNSS) transmitter ID
            where (b_int(1:lh) == default_head% sender_id) &
                   b_int(1:lh) =  -huge(b_int(1))
            ihbv = -1
          case ('tracking')                    ! TEMP PILOT
            b_int(1:lh) = spt(1:lh)% tracking  ! tracking technique,status of system
            where (b_int(1:lh) == default_head% tracking) &
                   b_int(1:lh) =  -huge(b_int(1))
            ihbv = -1
          case ('meas_type')                   ! TEMP PILOT
            b_int(1:lh) = spt(1:lh)% meas_type ! solar and infrared radiation correction
            where (b_int(1:lh) == default_head% meas_type) &
                   b_int(1:lh) =  -huge(b_int(1))
            ihbv = -1
          case ('flg_cld')                     ! RAD
            b_int(1:lh) = spt(1:lh)% stclf     ! station cloud flag
            ihbv = -1
          case ('sat_zenit')                   ! RAD
            b_real(1:lh) = spt(1:lh)% stzen    ! satellite zenith angle
            ihbv = -1
          case ('sat_azimuth')                 ! RAD
            b_real(1:lh) = spt(1:lh)% stazi    ! satellite azimuth angle
            ihbv = -1
          case ('sun_zenit')                   ! RAD
            b_real(1:lh) = spt(1:lh)% sozen    ! solar zenith angle
            ihbv = -1
          case ('sun_azimuth')                 ! RAD
            b_real(1:lh) = spt(1:lh)% soazi    ! solar azimuth angle
            ihbv = -1
          case ('retrtype')                    ! SATOB, GPSGB
            b_int(1:lh) = spt(1:lh)% stret     ! satellite derived wind computation method
            ihbv = -1
          case ('center_id')                   ! GPSGB
            b_int(1:lh) = spt(1:lh)% center_id ! ID of processing center
            ihbv = -1
          !-------------------------------------------
          ! verification data not handled in this loop
          !-------------------------------------------
          case ('veri_data','veri_model','veri_run_type','veri_initial_date',&
                'veri_forecast_time','veri_resolution','veri_domain_size',   &
                'veri_description','veri_ens_member','veri_exp_id','veri_run_class')
            cycle
          !---------------------
          ! radar specific table
          !---------------------
          case ('radar_azimuth')
            b_real(1:lo(OT_RADAR)) = radar% radar_azimuth
            ihbv = OT_RADAR
          case ('radar_elevation')
            b_real(1:lo(OT_RADAR)) = radar% radar_elevation
            ihbv = OT_RADAR
          case ('radar_nrange')
            b_int(1:lo(OT_RADAR))  = radar% radar_nrange
            ihbv = OT_RADAR
          case ('radar_range_start')
            b_real(1:lo(OT_RADAR)) = radar% radar_range_start
            ihbv = OT_RADAR
          case ('radar_drange')
            b_real(1:lo(OT_RADAR)) = radar% radar_drange
            ihbv = OT_RADAR
          case ('radar_nbody')
            b_int(1:lo(OT_RADAR))  = radar% radar_nbody
            ihbv = OT_RADAR
          !------------------------------------
          ! abort for variables not implemented
          !------------------------------------
          case default

!           call finish('write_fdbk_3dv',&
!           'variable(NF_INT) not implemented: '//fdbk% fb% nc% vars(iv)% name)

            cycle
          end select
          !---------------------------------------------------------------
          ! check for correct 'ihbv', set number of words for send/receive
          ! Assert consistency of offsets and counts
          !---------------------------------------------------------------
          select case (ihbv)
          case (1:n_ot)
            lbh = lo(ihbv)
            if (.not. all ( (fdbk% pelist(if, 0:dace% npe-2)% o_spec(ihbv) + &
                             fdbk% pelist(if, 0:dace% npe-2)% n_spec(ihbv) ) &
                          == fdbk% pelist(if, 1:dace% npe-1)% o_spec(ihbv) ) ) then
               write(0,*) 'pe,offsets=', dace% pe, fdbk% pelist(if, :)% o_spec(ihbv)
               write(0,*) 'pe,counts =', dace% pe, fdbk% pelist(if, :)% n_spec(ihbv)
               call finish ('write_fdbk_3dv(spec)',&
                            'offsets and counts not consistent')
            end if
            start  = fdbk% pelist(if, 0)% o_spec(ihbv) + 1
            counts = fdbk% pelist(if, 0:dace% npe-1)% n_spec(ihbv)
          case (-1)
            lbh = lh
            if (.not. all ( (fdbk% pelist(if, 0:dace% npe-2)% o_hdr + &
                             fdbk% pelist(if, 0:dace% npe-2)% n_hdr ) &
                          == fdbk% pelist(if, 1:dace% npe-1)% o_hdr ) ) then
               write(0,*) 'pe,offsets=', dace% pe, fdbk% pelist(if, :)% o_hdr
               write(0,*) 'pe,counts =', dace% pe, fdbk% pelist(if, :)% n_hdr
               call finish ('write_fdbk_3dv(header)',&
                            'offsets and counts not consistent')
            end if
            start  = fdbk% pelist(if, 0)% o_hdr + 1
            counts = fdbk% pelist(if, 0:dace% npe-1)% n_hdr
          case (-2)
            lbh = lb
            if (.not. all ( (fdbk% pelist(if, 0:dace% npe-2)% o_body + &
                             fdbk% pelist(if, 0:dace% npe-2)% n_body ) &
                          == fdbk% pelist(if, 1:dace% npe-1)% o_body ) ) then
               write(0,*) 'pe,offsets=', dace% pe, fdbk% pelist(if, :)% o_body
               write(0,*) 'pe,counts =', dace% pe, fdbk% pelist(if, :)% n_body
               call finish ('write_fdbk_3dv(body)',&
                            'offsets and counts not consistent')
            end if
            start  = fdbk% pelist(if, 0)% o_body + 1
            counts = fdbk% pelist(if, 0:dace% npe-1)% n_body
          case (-3)
            cycle
          case default
            call finish('write_fdbk_3dv','ihbv not in 1..3')
            cycle
          end select

!!! New version:
          !--------------------------
          ! gather data/write at once
          !--------------------------
          cnt = sum (counts)
          if (cnt > 0) then
             if (.not.dace% lpio) cnt = 0   ! Receive buffer on non-I/O PEs
             select case (fdbk% fb% nc% vars(iv)% xtype)
             case (NF90_INT, NF90_SHORT, NF90_BYTE)
                where (b_int(1:lbh) == -huge(b_int(1))) &
                       b_int(1:lbh) =  fdbk% fb% nc% vars(iv)% invalid
                allocate (r_int(cnt))
                call p_gatherv (sendbuf=b_int(1:lbh), recvbuf=r_int, &
                                root=dace% pio, recvcounts=counts    )
                if (dace% lpio) then
                   ir = nf90_put_var (fdbk% fb% nc% ncid,           &
                                      fdbk% fb% nc% vars(iv)% varid,&
                                      r_int(1:cnt),                 &
                                      [start],                      &
                                      [cnt]                         )
                   if (ir/=0) then
                      write(0,*) fdbk% fb% nc% vars(iv)% name,'=',&
                           minval(r_int(1:cnt)),maxval(r_int(1:cnt)),&
                           ', file =',trim(path)
                      call finish ('write_fdbk_3dv:nf90_put_var::int',&
                                   fdbk% fb% nc% vars(iv)% name//':'//&
                                   nf90_strerror(ir))
                   endif
                endif
                deallocate (r_int)
             case (NF90_REAL)
                where (b_real(1:lbh) <  -huge(b_real(1))/2) &
                       b_real(1:lbh) =  fdbk% fb% nc% vars(iv)% rinvalid
                allocate (r_real(cnt))
                call p_gatherv (sendbuf=b_real(1:lbh), recvbuf=r_real, &
                                root=dace% pio, recvcounts=counts      )
                if (dace% lpio) then
                   ir = nf90_put_var (fdbk% fb% nc% ncid,           &
                                      fdbk% fb% nc% vars(iv)% varid,&
                                      r_real(1:cnt),                &
                                      [start],                      &
                                      [cnt]                         )
                   if (ir/=0) then
                      write(0,*) fdbk% fb% nc% vars(iv)% name,'=',&
                           minval(r_real(1:cnt)),maxval(r_real(1:cnt)),&
                           ', file =',trim(path)
                      call finish ('write_fdbk_3dv:nf90_put_var::real',&
                                   fdbk% fb% nc% vars(iv)% name//':'// &
                                   nf90_strerror(ir))
                   endif
                endif
                deallocate (r_real)
             case (NF90_CHAR)
               if (.not.dace% lpio) cnt = 1    ! Needed for p_gatherv!

               ! Handle predefined (fixed) character string lengths:
               select case (clen)
               case (10)
                allocate (r_c10(cnt))
                call p_gatherv (sendbuf=b_c10(1:lbh), recvbuf=r_c10, &
                                root=dace% pio, recvcounts=counts    )
                if (dace% lpio) then
                   ir = nf90_put_var (fdbk% fb% nc% ncid,           &
                                      fdbk% fb% nc% vars(iv)% varid,&
                                      r_c10(1:cnt),                 &
                                      (/1,start/),                  &
                                      (/10,cnt/)                    )
                   if (ir/=0) then
                     call finish ('write_fdbk_3dv:nf90_put_var::text',&
                                  fdbk% fb% nc% vars(iv)% name//':'// &
                                  nf90_strerror(ir))
                   endif
                endif
                deallocate (r_c10)

               case (32)
                allocate (r_c32(cnt))
                call p_gatherv (sendbuf=b_c32(1:lbh), recvbuf=r_c32, &
                                root=dace% pio, recvcounts=counts    )
                if (dace% lpio) then
                   ir = nf90_put_var (fdbk% fb% nc% ncid,           &
                                      fdbk% fb% nc% vars(iv)% varid,&
                                      r_c32(1:cnt),                 &
                                      (/1,start/),                  &
                                      (/32,cnt/)                    )
                   if (ir/=0) then
                     call finish ('write_fdbk_3dv:nf90_put_var::text', &
                                  fdbk% fb% nc% vars(iv)% name//':'//  &
                                  nf90_strerror(ir))
                   endif
                endif
                deallocate (r_c32)

               case default
                  call finish ('write_fdbk_3dv:nf90_put_var::text', &
                               fdbk% fb% nc% vars(iv)% name //      &
                               ': bad character length'             )
               end select

             end select
          endif
!!! Old version: (kept for debugging)
!         !-----
!         ! send
!         !-----
!         if (.not.dace% lpio) then
!           if (lbh > 0) then
!             select case (fdbk% fb% nc% vars(iv)% xtype)
!             case (NF_INT, NF_SHORT, NF_BYTE)
!               call p_send (b_int (1:lbh), dace% pio, 1)
!             case (NF_REAL)
!               call p_send (b_real(1:lbh), dace% pio, 1)
!             case (NF_CHAR)
!               call p_send (b_c10 (1:lbh), dace% pio, 1)
!             end select
!           endif
!         !------------------
!         ! receive and write
!         !------------------
!         else
!           do ip = 0, dace% npe-1
!             select case (ihbv)
!             case (1)
!               start = fdbk% pelist(if, ip)% o_hdr + 1
!               cnt   = fdbk% pelist(if, ip)% n_hdr
!             case (2)
!               start = fdbk% pelist(if, ip)% o_body + 1
!               cnt   = fdbk% pelist(if, ip)% n_body
!             end select
!
!             if (cnt > 0) then
!               select case (fdbk% fb% nc% vars(iv)% xtype)
!               case (NF_INT, NF_SHORT, NF_BYTE)
!                 if (ip /= dace% pio) then
!                   call p_recv (b_int(1:cnt), ip, 1)
!                 endif
!                 where (b_int(1:cnt) == -huge(b_int(1))) &
!                        b_int(1:cnt) =  fdbk% fb% nc% vars(iv)% invalid
!                 ir = nf_put_vara_int (fdbk% fb% nc% ncid,           &
!                                       fdbk% fb% nc% vars(iv)% varid,&
!                                       start,                        &
!                                       cnt,                          &
!                                       b_int(1:cnt)                  )
!                 if (ir/=0) then
!                   write(0,*) fdbk% fb% nc% vars(iv)% name,'=',&
!                     minval(b_int(1:cnt)),maxval(b_int(1:cnt)),&
!                     ', file =',trim(path)
!                   call finish ('write_fdbk_3dv:nf_put_vara_int',            &
!                                          fdbk% fb% nc% vars(iv)% name//':'//&
!                                          nf_strerror(ir))
!                 endif
!               case (NF_REAL)
!                 if (ip /= dace% pio) then
!                   call p_recv (b_real(1:cnt), ip, 1)
!                 endif
!                 where (b_real(1:cnt) <  -huge(b_real(1))/2) &
!                        b_real(1:cnt) =  fdbk% fb% nc% vars(iv)% rinvalid
!                 ir = nf_put_vara_real (fdbk% fb% nc% ncid,           &
!                                        fdbk% fb% nc% vars(iv)% varid,&
!                                        start,                        &
!                                        cnt,                          &
!                                        b_real(1:cnt)                 )
!                 if (ir/=0) &
!                   call finish ('write_fdbk_3dv:nf_put_vara_real:', &
!                                 fdbk% fb% nc% vars(iv)% name//':'//&
!                                 nf_strerror(ir))
!               case (NF_CHAR)
!                 if (ip /= dace% pio) then
!                   call p_recv (b_c10(1:cnt), ip, 1)
!                 endif
!                 ir = nf_put_vara_text (fdbk% fb% nc% ncid,           &
!                                        fdbk% fb% nc% vars(iv)% varid,&
!                                        (/1,start/),                  &
!                                        (/10,cnt/),                   &
!                                        b_c10(1:cnt)                  )
!                 if (ir/=0) &
!                   call finish ('write_fdbk_3dv:nf_put_vara_text:', &
!                                 fdbk% fb% nc% vars(iv)% name//':'//&
!                                 nf_strerror(ir))
!               end select
!             endif
!           end do
!         endif

        end do  ! Loop over iv
        if (associated (radar)) deallocate (radar)
        if (l_tovs) then
          if (lh > 0) deallocate(l2c_t, l2c_t_, emiss_t, sinfl_t, tflag_t, ci)
          deallocate(emiss, tflag, sinfl, scanl)
          if (allocated(l2c    )) deallocate(l2c    )
          if (allocated(nwc_flg)) deallocate(nwc_flg)
          if (allocated(orb_ph )) deallocate(orb_ph )
          if (allocated(ins_tmp)) deallocate(ins_tmp)
          if (allocated(cld_frc)) deallocate(cld_frc)
        end if
        !------------------------
        ! write verification data
        !------------------------
        if (present (fg))       call add_veri (fg,       VT_FIRSTGUESS, VE_DETERM,      lmodel, 'first guess'                  )
        if (present (e_fg))     call add_veri (e_fg,     VT_FIRSTGUESS, VE_BG_ERROR,    lmodel, 'first guess error'            )
        if (present (ana_psas)) call add_veri (ana_psas, VT_PREL_ANA,   VE_DETERM,      lmodel, 'analysis in observation space')
        if (present (ana))      call add_veri (ana,      VT_ANALYSIS,   VE_DETERM,      lmodel, '3dvar analysis'               )
        if (present (s_fg))     call add_veri (s_fg,     VT_FIRSTGUESS, VE_ENS_SPREAD,  lmodel, 'first guess ensemble spread'  )
        if (present (s_ana))    call add_veri (s_ana,    VT_ANALYSIS,   VE_ENS_SPREAD,  lmodel, 'analysis ensemble spread'     )
        if (present (m_fg))     call add_veri (m_fg,     VT_FIRSTGUESS, VE_ENS_MEAN,    lmodel, 'first guess ensemble mean'    )
        if (present (m_ana))    call add_veri (m_ana,    VT_ANALYSIS,   VE_ENS_MEAN,    lmodel, 'analysis ensemble mean'       )
        if (present (mo_fg))    call add_veri (mo_fg,    VT_FIRSTGUESS, VE_ENS_MEAN_OBS,lmodel, &
                                                                               'first guess ensemble mean in observation space')
        if (present (mo_ana))   call add_veri (mo_ana,   VT_ANALYSIS,   VE_ENS_MEAN_OBS,lmodel, &
                                                                               'analysis ensemble mean in observation space'   )
        if (present (tg_fg))    call add_veri (tg_fg,    VT_FIRSTGUESS, VE_TALAGRAND,   lmodel, 'first guess Talagrand index'  )
        if (present (tg_ana))   call add_veri (tg_ana,   VT_ANALYSIS,   VE_TALAGRAND,   lmodel, 'analysis Talagrand index'     )
        if (present (w_qc))     call add_veri (w_qc,     VT_ANALYSIS,   VE_VQC_WEIGHT,  lmodel, 'VQC weight'                   )
        if (present (varbc))    call add_veri (varbc,    VT_PREL_ANA,   VE_BIASCOR,     lmodel, 'variational bias correction'  )
        if (present (en_fg)) then
          do k = 1, size(en_fg)
                                call add_veri (en_fg(k), VT_FIRSTGUESS, k, lmodel, 'first guess ensemble member '//char4(k))
          end do
        endif
        if (present (en_ana)) then
          do k = 1, size(en_ana)
                                call add_veri (en_ana(k),VT_ANALYSIS,   k, lmodel, 'analysis ensemble member '//char4(k))
          end do
        endif

        if (present    (lin_ana)) then
        if (associated (lin_ana)) then
          do k = 1, size(lin_ana)
                                call add_veri (lin_ana(k),VT_LIN_ANA,k, lmodel, 'Ya linear operator on analysis '//char4(k))
          end do
        endif
        endif


!       if (fdbk% file(if)% name == 'RAD') then
!                               call add_veri (fg_1d,    VT_FIRSTGUESS, VE_DETERM,   '1DVAR', '1dvar first guess')
!                               call add_veri (an_1d,    VT_ANALYSIS,   VE_DETERM,   '1DVAR', '1dvar analysis')
!       endif

        !-------------------------------------
        ! close the file, deallocate meta data
        !-------------------------------------
        if (dace% lpio) &
          call close_fdbk (fdbk% fb)
      endif
      call cleanup_fdbk (fdbk% fb)
    end do
    !------------------
    ! deallocate buffer
    !------------------
    deallocate (b_int, b_real, b_c10, ix, spt, bdy)
    if (allocated (b_c32)) deallocate (b_c32)

    if (present(e_o))      call delete_storage (e_o)
    if (present(w_qc))     call delete_storage (w_qc)
    if (present(fg))       call delete_storage (fg)
    if (present(e_fg))     call delete_storage (e_fg)
    if (present(ana))      call delete_storage (ana)
    if (present(ana_psas)) call delete_storage (ana_psas)
    if (present(varbc))    call delete_storage (varbc)
    if (present(s_fg))     call delete_storage (s_fg)
    if (present(s_ana))    call delete_storage (s_ana)
    if (present(m_fg))     call delete_storage (m_fg)
    if (present(m_ana))    call delete_storage (m_ana)
    if (present(mo_fg))    call delete_storage (mo_fg)
    if (present(mo_ana))   call delete_storage (mo_ana)
    if (present(tg_fg))    call delete_storage (tg_fg)
    if (present(tg_ana))   call delete_storage (tg_ana)
    call destruct (mask)
    call leave_function
contains
!------------------------------------------------------------------------------
  subroutine add_veri (x, runtype, member, model, description)
  type(t_vector)   ,intent(in) :: x            ! values to write
  integer          ,intent(in) :: runtype      ! run type flag
  integer          ,intent(in) :: member       ! member flag
  character(len=*) ,intent(in) :: model        ! model
  character(len=*) ,intent(in) :: description  ! description

    integer      :: iv  ! index in extendable dimension 'n_veri'
    integer      :: i,ib,is,i1,in,ir
    integer      :: start, count
    integer      :: varid
    integer      :: ii1,iin
    type(t_time) :: time
    real(sp), allocatable :: recvbuf(:) ! Receive buffer for MPI gatherv
    type(t_obs) ,pointer  :: o          ! pointer to current observation box
    logical     ,pointer  :: m(:)

    !---------------------------------------------
    ! return if data components are not associated
    !---------------------------------------------
    if (.not. associated (x% s)) return
    !----------------
    ! write meta data
    !----------------
    if (dace% lpio) then
      select case (runtype)
      case (VT_FIRSTGUESS)
        time = fc_ref_time
      case default
        time =    ana_time
      end select
      select case (grid% gridtype)
      case default
        nx = grid% nx
        ny = grid% ny
        dx = grid% di
        dy = grid% dj
      case (DWD6_ICOSAHEDRON)
        nx = grid% ni
        ny = nx
        dx = 63.486_sp/nx
        dy = dx
      case (DWD6_ICON)
        nx = grid% ni
        ny = nx
        dx = 45.416_sp/nx
        dy = dx
      end select
      call add_verification (        &
        fdbk% fb,                    &! feedback file meta data
        model,                       &! model
        runtype,                     &! run type
        run_type,                    &! run class
        cyyyymmddhhmm(time),         &! initial date
        ihhmm(ana_time-time),        &! forecast time (hhhmm)
        (/dx, dy/),                  &! resolution
        (/nx, ny, grid%nz/),         &! domain
        description,                 &! description
        member,                      &! ensemble member flag
        nex,                         &! experiment id
        varid)                        ! var.id. of 'veri_data'
      iv = fdbk% fb% n_veri
    endif
    !-----------
    ! write data
    !-----------
    if (associated(x% s)) then
      !----------
      ! pack body
      !----------
      do i = 1, lh
      ib = ix(i)% ib
        o   => obs% o(ib)
        m   => mask% s(ib)% x
        is  =  ix(i)% is
        i1  =  ix(i)% o_body + 1
        in  =  ix(i)% n_body - 1 + i1
        ii1 =  o% spot(is)% o% i + 1
        iin =  o% spot(is)% o% i + o% spot(is)% o% n
        b_real(i1:in) = pack (x% s(ib)% x(ii1:iin), m(ii1:iin))
      end do

!!! Old version:
!     !------------------------
!     ! send/receive/write data
!     !------------------------
!     if (.not.dace% lpio) then
!       if (lb > 0) call p_send (b_real(1:lb), dace% pio, 1)
!     else
!       do ip = 0, dace% npe-1
!         start = fdbk% pelist(if, ip)% o_body + 1
!         count = fdbk% pelist(if, ip)% n_body
!         if (count > 0) then
!           if (ip /= dace% pio) call p_recv (b_real(1:count), ip, 1)
!           ir = nf_put_vara_real (fdbk% fb% nc% ncid,           &
!                                  varid,                        &
!                                  (/start,iv/),                 &
!                                  (/count, 1/),                 &
!                                  b_real(1:count)               )
!           if (ir/=0) call finish ('add_veri:nf_put_vara_real:',&
!                                  description//':'//            &
!                                  nf_strerror(ir))
!         endif
!       end do
!     endif
!!! New version:
      ! Assert consistency of offsets and counts, otherwise we need a fallback
      if (.not. all ( (fdbk% pelist(if, 0:dace% npe-2)% o_body + &
                       fdbk% pelist(if, 0:dace% npe-2)% n_body ) &
                    == fdbk% pelist(if, 1:dace% npe-1)% o_body ) ) then
         write(0,*) 'pe,offsets=', dace% pe, fdbk% pelist(if, :)% o_body
         write(0,*) 'pe,counts =', dace% pe, fdbk% pelist(if, :)% n_body
         call finish ('write_fdbk_3dv:add_veri',&
                      'offsets and counts not consistent')
      end if
      !--------------------------
      ! gather data/write at once
      !--------------------------
      start =      fdbk% pelist(if, 0)% o_body + 1
      count = sum (fdbk% pelist(if, 0:dace% npe-1)% n_body)
      if (count > 0) then
         if (.not.dace% lpio) count = 0
         allocate (recvbuf(count))
         call p_gatherv (sendbuf=b_real(1:lb), recvbuf=recvbuf, root=dace% pio, &
                         recvcounts=fdbk% pelist(if, 0:dace% npe-1)% n_body      )
         if (dace% lpio) then
            ir = nf90_put_var (fdbk% fb% nc% ncid, &
                               varid,              &
                               recvbuf(1:count),   &
                               (/start,iv/),       &
                               (/count, 1/)        )
            if (ir/=0) call finish ('add_veri:nf90_put_var::real',&
                                     description//':'//           &
                                     nf90_strerror(ir))
         end if
         deallocate (recvbuf)
      end if

    endif

  end subroutine add_veri
!------------------------------------------------------------------------------
end subroutine write_fdbk_3dv
!==============================================================================
  subroutine change_bits (bits, table)
  integer ,intent(inout) :: bits  (:)
  integer ,intent(in)    :: table (:)
    integer :: i,n,tmp
    do n=1,size(bits)
      tmp     = bits(n)
      bits(n) = 0
      do i=1,size(table)
        if (btest(tmp,i)) then
          if ((table(i)>=0) .and. (table(i)<=31)) then
            bits(n) = ibset(bits(n),table(i))
          endif
        endif
      end do
    end do
  end subroutine change_bits
!-----------------------------------------------------------------------------
  subroutine change_code (code, table)
  integer ,intent(inout) :: code  (:)
  integer ,intent(in)    :: table (:)
    integer :: n
    do n=1,size(code)
      code(n) = table(code(n))
    end do
  end subroutine change_code
!==============================================================================
  subroutine add_veri_3dv  (fdbk, obs, grid, w_qc, ana, prel_ana, varbc)
  !-----------------------------------------------------
  ! extend the feedback file: write forecast or analysis
  !-----------------------------------------------------
  type(t_fdbk_3dv),intent(inout)        :: fdbk     ! feedback file data
  type(t_obs)     ,intent(inout)        :: obs (:)  ! observational data
  type(t_grid)    ,intent(in)           :: grid     ! model grid info
  type(t_vector)  ,intent(in) ,optional :: w_qc     ! VQC weight
  type(t_vector)  ,intent(in) ,optional :: ana      ! analysis
  type(t_vector)  ,intent(in) ,optional :: prel_ana ! preliminary psas analysis
  type(t_vector)  ,intent(in) ,optional :: varbc    ! variational bias corr.

    integer                    :: if     ! file index
    integer                    :: nf     ! number of files
    integer                    :: lh     ! header index
    integer                    :: lb     ! body index
    integer                    :: lhpe (0:dace% npe-1, size (fdbk% file))
    integer                    :: lbpe (0:dace% npe-1, size (fdbk% file))
    integer                    :: mh     ! max no. header entries / pe
    integer                    :: mb     ! max no. body   entries / pe
    integer                    :: mr     ! max no. body   entries
    integer                    :: n_hdr  !     no. header entries
    integer                    :: n_body !     no. header entries
    integer                    :: ir     ! NetCDF return argument
    integer                    :: varid  ! NetCDF variable id
    real(wp)      ,allocatable :: b_s    (:) ! real numbers
    real(sp)      ,allocatable :: b_out  (:) ! output buffer
    integer       ,allocatable :: i_body (:) ! pointer hdr -> body
    integer       ,allocatable :: state  (:) ! status of observation
    integer       ,allocatable :: ih_s   (:) ! report index    (send buffer)
    integer       ,allocatable :: ib_s   (:) ! body index      (send buffer)
    integer       ,allocatable :: nb_s   (:) ! no.body entries (send buffer)
    !-------------------------------------------------
    ! to be written to NetCDF file (global attributes)
    !-------------------------------------------------
    real(sp)           :: dx, dy  ! model resolution (degree)
    integer            :: nx, ny  ! number of model gridpoints

    integer :: ib    ! observation box index
    integer :: is    ! record index
    integer :: i0,in ! observation index range

    call enter_function

    !---------------------------------
    ! derive model resolution and size
    !---------------------------------
    select case (grid% gridtype)
    case (WMO6_LATLON,WMO6_ROTLL,WMO6_GAUSSIAN)
      nx = grid% nx
      ny = grid% ny
      dx = grid% di
      dy = grid% dj
    case (DWD6_ICOSAHEDRON)
      nx = grid% ni
      ny = nx
      dx = 63.486_sp/nx
      dy = dx
    case (DWD6_ICON)
      nx = grid% ni
      ny = nx
      dx = 45.416_sp/nx
      dy = dx
    case default
      write(0,*) "gridtype =", grid% gridtype
      call finish('write_fdbk_3dv','unknown gridtype')
    end select

    !--------------------
    ! count array lengths
    !--------------------
    nf = size (fdbk% file)
    lhpe = 0
    lbpe = 0
    do if = 1, nf
      n_hdr  = fdbk% file(if)% n_hdr
      n_body = fdbk% file(if)% n_body
      if (n_hdr > 0 .and. n_body > 0) then
        lh = 0
        lb = 0
        do ib = 1, size(obs)
          if (obs(ib)% pe == dace% pe) then
            do is = 1, obs(ib)% n_spot
              if (obs(ib)% spot(is)% hd% mon_file == if) then
                in = obs(ib)% spot(is)% o% n
                lb = lb + obs(ib)% spot(is)% o% n
                lh = lh + 1
              endif
            end do
          endif
          lhpe(dace% pe, if) = lh
          lbpe(dace% pe, if) = lb
        end do
      endif
    enddo
    lhpe = p_max (lhpe)
    lbpe = p_max (lbpe)

    !----------------
    ! allocate arrays
    !----------------
    mh = maxval (lhpe (dace% pe, :))
    mb = maxval (lbpe (dace% pe, :))
    mr = maxval (fdbk% file% n_body)
    allocate (ih_s   (mh))
    allocate (nb_s   (mh))
    allocate (ib_s   (mb))
    allocate (b_s    (mb))
    if (dace% lpio) then
      mh = maxval (lhpe)
      mb = maxval (lbpe)
      allocate (b_out  (mr))
      allocate (i_body (mr))
    endif

    !--------------------------
    ! loop over non-empty files
    !--------------------------
    nf = size (fdbk% file)
    do if = 1, nf
      n_hdr  = fdbk% file(if)% n_hdr
      n_body = fdbk% file(if)% n_body
      if (n_hdr > 0 .and. n_body > 0) then
        !---------------------
        ! open the NetCDF file
        !---------------------
        if (dace% lpio) then
          call open_fdbk_write (fdbk% fb, fdbk% file(if)% path)
          call read_meta (fdbk% fb)
          ir = nf90_inq_varid  (fdbk% fb% nc% ncid,'i_body',varid)
          if (ir/=0) call finish ('add_veri_3dv:nf90_inq_varid', &
                                  'i_body:'//nf90_strerror(ir)   )
          ir = nf90_get_var    (fdbk% fb% nc% ncid, &
                                varid,              &
                                i_body(1:n_hdr)     )
          if (ir/=0) call finish ('add_veri_3dv:nf90_get_var::int',&
                                  'i_body:'//nf90_strerror(ir)     )
        endif

        !----------------
        ! set index array
        !----------------
        lh = 0
        lb = 0
        do ib = 1, size(obs)
          if (obs(ib)% pe == dace% pe) then
            do is = 1, obs(ib)% n_spot
              if (obs(ib)% spot(is)% hd% mon_file == if) then
                lh                = lh + 1
                i0                = obs(ib)% spot(is)% o% i
                in                = obs(ib)% spot(is)% o% n
                ih_s (lh)         = obs(ib)% spot(is)% hd% mon_rec
                nb_s (lh)         = in
                ib_s (lb+1:lb+in) = obs(ib)% body(i0+1:i0+in)% mon_pos
                lb                = lb + in
              endif
            end do
          endif
        end do

        !------------------------
        ! write verification data
        !------------------------
        if (present (ana))                                                    &
          call add_veri(ana,     VT_ANALYSIS,VE_DETERM,    rmodel,            &
                                                              '3dvar analysis')
        if (present (prel_ana))                                               &
          call add_veri(prel_ana,VT_PREL_ANA,VE_DETERM,    rmodel,            &
                                               'analysis in observation space')
        if (present (w_qc))                                                   &
          call add_veri(w_qc,    VT_ANALYSIS,VE_VQC_WEIGHT,rmodel,'VQC weight')
        if (present (varbc))                                                  &
          call add_veri(varbc,   VT_PREL_ANA,VE_BIASCOR,   rmodel,            &
                                                 'Variational Bias Correction')
        !----------------------
        ! close the NetCDF file
        !----------------------
        if (dace% lpio) then
          call close_fdbk   (fdbk% fb)
          call cleanup_fdbk (fdbk% fb)
        endif
      endif
    end do
    if (present (w_qc))     call delete_storage (w_qc)
    if (present (ana))      call delete_storage (ana)
    if (present (prel_ana)) call delete_storage (prel_ana)
    if (present (varbc))    call delete_storage (varbc)
    call leave_function
  contains
!------------------------------------------------------------------------------
    subroutine add_veri (x, runtype, member, model, description)
    type(t_vector)   ,intent(in) :: x            ! values to write
    integer          ,intent(in) :: runtype      ! run type flag
    integer          ,intent(in) :: member       ! member flag
    character(len=*) ,intent(in) :: model        ! model
    character(len=*) ,intent(in) :: description  ! description

      integer :: iv(1) ! index in extendable dimension 'n_veri'
      integer :: nidx  ! number of matches in verification entries
      integer :: pe    ! processor index
      integer :: lh_r  ! number of headers (received value)
      integer :: lb_r  ! number of body entries  (received)
      integer :: ih    ! header index
      integer :: iho   ! header index offset
      integer :: ihg   ! header index (global)
      integer :: ibd   ! body index
      integer :: ibo   ! body index offset
      integer :: ibdl  ! body index (local)
      integer :: ibdg  ! body index (global)
      integer :: vd_id ! NetCDF varid for 'veri_data'
      integer :: st_id ! NetCDF varid for 'state'
      integer :: i0,in ! observation index range
      integer :: lh    ! header index
      integer :: lb    ! body index
      integer :: ib, is

      integer :: nh (0:dace% npe-1) ! number of headers
      integer :: nb (0:dace% npe-1) ! number of body entries
      integer :: nhtot
      integer :: nbtot
      integer  ,allocatable :: ih_r(:)   ! report index    (recv buffer)
      integer  ,allocatable :: ib_r(:)   ! body index      (recv buffer)
      integer  ,allocatable :: nb_r(:)   ! no.body entries (recv buffer)
      real(wp) ,allocatable ::  b_r(:)   ! real numbers    (recv buffer)

      !----------------------------------------
      ! check if data components are associated
      !----------------------------------------
      if (associated(x% s)) then
        !----------------
        ! write meta data
        !----------------
        if (dace% lpio) then

          call read_meta (fdbk% fb)
          call get_veri_index (iv, nidx, fdbk% fb, &
                       model = model,              &
                    run_type = runtype,            &
                  ens_member = member              )

          select case (nidx)
          case default
            call finish('add_veri','cannot find unique verification index')
          case (1)
            vd_id = get_varid (fdbk% fb, 'veri_data' )
          case (0)
            call add_verification (        &
              fdbk% fb,                    &! feedback file meta data
              model,                       &! model
              runtype,                     &! run type
              run_type,                    &! run class
              cyyyymmddhhmm(ana_time),     &! initial date
              0,                           &! forecast time (hhhmm)
              (/dx,dy/),                   &! resolution
              (/grid%ni,grid%ni,grid%nz/), &! domain
              description,                 &! description
              member,                      &! ensemble member flag
              nex,                         &! experiment id
              vd_id)                        ! var.id. of 'veri_data'
            iv = fdbk% fb% n_veri
          end select
        endif

        !----------
        ! pack body
        !----------
        lh = 0
        lb = 0
        do ib = 1, size(obs)
          if (obs(ib)% pe == dace% pe) then
            do is = 1, obs(ib)% n_spot
              if (obs(ib)% spot(is)% hd% mon_file == if) then
                lh = lh + 1
                i0 = obs(ib)% spot(is)% o% i
                in = obs(ib)% spot(is)% o% n
                b_s (lb+1:lb+in) = x% s(ib)% x(i0+1:i0+in)
                lb = lb + in
              endif
            end do
          endif
        end do

        if (dace% lpio) then
          !----------
          ! read data
          !----------
          if (nidx == 1) then
            ir = nf90_get_var (fdbk% fb% nc% ncid, &
                               vd_id,              &
                               b_out(1:n_body),    &
                               (/1,iv(1)/),        &
                               (/n_body, 1/)       )
            if (ir/=0) call finish ('add_veri:nf90_get_var::real',&
                                     description//':'//           &
                                     nf90_strerror(ir))
          else
            b_out(1:n_body) = NF90_FILL_REAL
          endif
        end if
        !--------------------------------
        ! gather updates on I/O processor
        !--------------------------------
        nh = 0
        nb = 0
        call p_gather (lh, nh, root=dace% pio)
        call p_gather (lb, nb, root=dace% pio)
        nhtot = sum (nh)
        nbtot = sum (nb)
        allocate (ih_r(nhtot))
        allocate (nb_r(nhtot))
        allocate (ib_r(nbtot))
        allocate ( b_r(nbtot))
        call p_gatherv (ih_s(1:lh),ih_r,root=dace% pio) ! hd% mon_rec
        call p_gatherv (nb_s(1:lh),nb_r,root=dace% pio) ! o% n
        call p_gatherv (ib_s(1:lb),ib_r,root=dace% pio) ! body% mon_pos
        call p_gatherv ( b_s(1:lb), b_r,root=dace% pio)
        if (dace% lpio) then
          ibo = 0
          iho = 0
          do pe = 0, dace% npe-1
            lh_r = nh(pe)
            lb_r = nb(pe)
            !-------------
            ! replace data
            !-------------
            do ih = 1, lh_r
              ihg = iho + ih
              do ibd = 1, nb_r(ihg)
                ibdl = ibo + ibd
                if   ( ib_r (ibdl) < 0) cycle
                ibdg = ib_r (ibdl) + i_body (ih_r(ihg)) - 1
                b_out (ibdg) = b_r (ibdl)
              end do
              ibo    = ibo + nb_r(ihg)
            end do
            iho = iho + lh_r
          end do
          !----------------
          ! write back data
          !----------------
          ir = nf90_put_var (fdbk% fb% nc% ncid, &
                             vd_id,              &
                             b_out(1:n_body),    &
                             (/1,iv(1)/),        &
                             (/n_body, 1/)       )
          if (ir/=0) call finish ('add_veri:nf90_put_var::real',&
                                   description//':'//           &
                                   nf90_strerror(ir))
        end if
        deallocate (ih_r)
        deallocate (nb_r)
        deallocate (ib_r)
        deallocate ( b_r)

        if (dace% lpio) then
          !-------------------------------------------
          ! cross check status and set 'accepted' flag
          !-------------------------------------------
          if (member == VE_VQC_WEIGHT) then
            allocate (state (n_body))
            ir = nf90_inq_varid (fdbk% fb% nc% ncid,'state',st_id    )
            if (ir/=0) call finish ('add_veri_3dv:nf90_inq_varid',   &
                                    'state:'//nf90_strerror(ir)      )
            ir = nf90_get_var   (fdbk% fb% nc% ncid, st_id, state    )
            if (ir/=0) call finish ('add_veri_3dv:nf90_get_var::int',&
                                    'state:'//nf90_strerror(ir)      )
            do ibd = 1, n_body
              if (b_out(ibd) /= NF90_FILL_REAL) then
                if (state(ibd) /= ST_ACTIVE) then
                  write(6,*) 'add_veri: mismatch: ibd, state, weight =',&
                              ibd, state(ibd), b_out(ibd)
                  call finish('add_veri','mismatch')
                elseif (b_out(ibd) > 0.5_wp ) then
                  state (ibd) = ST_ACCEPTED
                endif
              else
                if (state(ibd) == ST_ACTIVE) then
                  write(6,*) 'add_veri: mismatch: ibd, state, weight =',&
                              ibd, state(ibd), b_out(ibd)
                  call finish('add_veri','mismatch')
                endif
              endif
            end do
            ir = nf90_put_var (fdbk% fb% nc% ncid, st_id, state      )
            if (ir/=0) call finish ('add_veri_3dv:nf90_put_var::int',&
                                    'state:'//nf90_strerror(ir)      )
            deallocate (state)
          endif ! member == VE_VQC_WEIGHT
        endif   ! dace% lpio
      endif
    end subroutine add_veri
!------------------------------------------------------------------------------
  end subroutine add_veri_3dv
!==============================================================================
! define MPI broadcast routines p_bcast_t_file and p_bcast_t_var
! for the derived types t_file and t_netcdf_var
!==============================================================================
#define VECTOR
#undef  DERIVED
#define DERIVED type(t_file),DIMENSION(:)
#undef  p_bcast_DERIVED
#define p_bcast_DERIVED p_bcast_t_file
#undef MPI_TYPE
#include "p_bcast.incf"
!------------------------------------------------------------------------------
#define VECTOR
#undef  DERIVED
#define DERIVED type(t_netcdf_var),DIMENSION(:)
#undef  p_bcast_DERIVED
#define p_bcast_DERIVED p_bcast_t_var
#undef MPI_TYPE
#include "p_bcast.incf"
!==============================================================================
end module mo_fdbk_3dvar
