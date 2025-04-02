!
!+ module to calculate thresholds for various cloud detection tests
!
MODULE mo_cloud_params
!
! Description:
! Module to derive thresholds for various flavors of cloud detection tests.
! Read parameter files for some special cloud detection tests and calculate
! appropriate thresholds for these tests. For more simple tests with
! type(t_range_fparse) thresholds, the get_threshold routine (provided by this
! module) should be used to evaluate the threshold.
!
! Heritage: mo_amsub.f90
!
! Current Maintainer: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Mashrab Kuvatov
!  Added two more cloud filtering methods
! V1_13        2011/11/01 Mashrab Kuvatov
!  changes for general # of scanlines
! V1_48        2016-10-06 Robin Faulwetter
!  add print statement
! V1_50        2017-01-09 Robin Faulwetter
!  Restructured/unified cloud detection for radiances.
! V2_22        2024-07-23 Robin Faulwetter
!  Renamed from mo_amsub to mo_cloud_params, added routines.
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin   DWD  2008 original code
! Mashrab Kuvatov  DWD  2008 original code
! Robin Faulwetter DWD  2016- modifications
!
!==============================================================================


  !=============
  ! Modules used
  !=============

  use mo_kind,              only: wp
  use mo_exception,         only: finish
  use mo_mpi_dace,          only: dace,               &
                                  p_bcast
  use mo_dace_string,       only: split2
  use mo_fortran_units,     only: get_unit_number,    &
                                  return_unit_number
  use mo_run_params,        only: data,               &
                                  path_file
  use mo_t_obs,             only: t_spot
  use mo_rad,               only: sat_poly2fparse
  use mo_range_fparse,      only: t_range_fparse,     &
                                  evaluate,           &
                                  destruct,           &
                                  init,               &
                                  p_bcast
  use mo_tovs_prof,         only: get_tovs_var,       &
                                  c_var
  use mo_physics,           only: d2r


  implicit none

  !================
  ! public entities
  !================
  private
  public :: inv_thr
  public :: get_threshold
  ! Scattering indices thresholds
  public :: param_file_scatt      ! File with thresholds for MHS-like scattering test
  public :: param_file_si_amsua   ! File with thresohlds for AMSUA-like scattering test
  public :: read_thresh_params    ! Read file with threshold parameters
  ! AMSUB/MHS thresholds
  public :: get_threshold_buehler
  public :: param_amsub       ! AMSUB cloud detection parameters
  public :: t_param_amsub     ! derived type
  public :: read_param_amsub  ! read AMSUB cloud detection parameters
  public :: p_bcast
  public :: dlat_thresh_def
  public :: dsz_thresh_def
  public :: buehler_version

  !============================================
  ! Module parameters, variables and data types
  !===========================================

  real(wp), parameter :: inv_thr = -huge(0._wp)

  !---------------------------
  ! Scattering index threshold
  !---------------------------
  character(len=300), save      :: param_file_scatt    = ''
  character(len=300), save      :: param_file_si_amsua = ''
  integer,            parameter :: mx_par = 30

  type t_thresh_par
    real(wp)                    :: lat   = -99._wp
    integer                     :: ndata = 0
    integer                     :: nfit  = 0
    type(t_range_fparse)        :: p
  end type t_thresh_par

  type t_thresh_set
    integer                     :: instr = -1
    integer                     :: satid = -1
    integer                     :: nlat  =  0
    integer                     :: npar  =  0
    type(t_thresh_par), pointer :: p(:)  => null()
  end type t_thresh_set

  type t_thresh
    integer                     :: file_version = -99
    integer                     :: nset         =   0
    type(t_thresh_set)          :: set(mx_par)
  end type t_thresh

  type(t_thresh), save, target  :: sc_par
  type(t_thresh), save, target  :: sia_par

  ! ------------------------------------
  ! AMSUB/MHS cloud detection thresholds
  ! ------------------------------------
  ! parameters for Buehler 2007 scheme
  integer,                 parameter :: nlat               = 36 ! number of latitude dependent coefficients
  integer,                 parameter :: nscan              = 45 ! number of scan dependant coefficients
  ! required for determining best parameter set
  real(wp),                parameter :: lat_scale          = 90._wp  ! scaling for diffs to nominal parameter latitudes
  real(wp),                parameter :: sz_scale           = 1._wp/cos(65._wp*d2r) - &  ! scaling for diffs to nominal parameter
                                                             1._wp/cos( 0._wp*d2r)      ! zenith angles. Diff between path through
                                                                                        ! atm. at a large zenith angle and nadir

  ! Modified Buehler scheme
  type t_buehler_par
    integer  :: n        = 0      ! number of data used to calc. parameters
    real(wp) :: lat      = 0._wp  ! latitude
    real(wp) :: stzen    = 0._wp  ! sat. zenith angle
    real(wp) :: bt_h_min = 0._wp  ! minimum accepted value for the high channel
    real(wp) :: a        = 0._wp  ! constant of BT_low-BT_high threshold
    real(wp) :: b        = 0._wp  ! slope of BT_low-BT_high threshold
                                  ! y = A + BT_high * B is a threshold on BT_low - BT_high
  end type t_buehler_par

  type t_buehler_mod
    integer                      :: instr       = -1
    integer                      :: satid       = -1
    integer                      :: nlat        =  0
    integer                      :: nfov        =  0
    integer                      :: npar        =  0
    type(t_buehler_par), pointer :: p(:)        => null()
    real(kind=wp)                :: dlat_thresh =  0._wp  ! maximum accepted scaled latitude difference
    real(kind=wp)                :: dsz_thresh  =  0._wp  ! maximum accepted scaled 1/cos(stzen) difference
    integer                      :: sz_mode     =  0      ! type sat.zen. angles in parameter file
                                                          ! 0: abs(sat.zen.)
                                                          ! 1:     sat.zen.
  end type t_buehler_mod

  real(kind=wp),       save      :: dlat_thresh_def =  0.1_wp  ! default value for maximum accepted scaled latitude difference
  real(kind=wp),       save      :: dsz_thresh_def  =  0.1_wp  ! default value for maximum accepted scaled 1/cos(stzen) difference
  real(kind=wp),       save      :: max_slope       =  0.0_wp  ! maximum accepted slope (parameter "B")
  integer,             save      :: buehler_version = 1        ! Buehler cloud check compatibility version

  ! ---------------------------
  ! type holding all parameters
  ! ---------------------------
  type t_param_amsub
    character(len=80)   :: check      = ''               ! FNCT_*
    real(wp)            :: ec_bnd (2) = (/-5._wp,5._wp/) ! bounds for ECMWF cloud check
    real(wp)            :: delta_ch20_ch18 = 0._wp       ! threshold of AMSU-B Ch(20) - Ch(18)
    integer             :: file_version = -99
    real(wp)            :: peak_ch_18 (nscan, nlat)      ! thresholds for old Buehler scheme
    integer             :: n_bm       = 0
    type(t_buehler_mod) :: bm(mx_par)                    ! thresholds for nodified Buehler scheme
  end type t_param_amsub

  type (t_param_amsub) ,save, target :: param_amsub


  real(wp) ,parameter     :: amsub_zen (nscan) =  & ! AMSU-B zenith angles
     (/ 0.62,  1.87,  3.11,  4.36,  5.60,  6.85,  8.09,  9.34, 10.59, 11.84, &
       13.09, 14.34, 15.59, 16.85, 18.11, 19.37, 20.63, 21.89, 23.16, 24.43, &
       25.71, 26.98, 28.26, 29.55, 30.84, 32.13, 33.43, 34.74, 36.05, 37.36, &
       38.69, 40.02, 41.36, 42.71, 44.07, 45.44, 46.83, 48.22, 49.64, 51.06, &
       52.51, 53.98, 55.47, 56.99, 58.54 /)

  interface p_bcast
    module procedure bcast_thresh
    module procedure bcast_thresh_par
    module procedure bcast_param_amsub
    module procedure bcast_buehler_par
  end interface

  interface destruct
    module procedure destruct_thresh_set
    module procedure destruct_thresh_par
    module procedure destruct_t_buehler_mod
  end interface

  interface construct
    module procedure construct_thresh_set
    module procedure construct_thresh_par
    module procedure construct_t_buehler_mod
  end interface

  interface get_threshold
    module procedure get_threshold_val
    module procedure get_threshold_par
  end interface get_threshold


#if (defined (__GFORTRAN__) && (__GNUC__ >= 10)) || defined (NAGFOR)
  interface
     subroutine p_bcast_derivedtype (buffer, count, source, comm)
       type(*) ,intent(inout)     :: buffer       ! variable to bcast
       integer ,intent(in)        :: count        ! len(byte) of variable
       integer ,intent(in)        :: source       ! source processor index
       integer ,intent(in)        :: comm         ! communicator
     end subroutine p_bcast_derivedtype
  end interface
#endif

contains


  subroutine read_thresh_params(stat)
    integer,           intent(out) :: stat

    type(t_thresh_set), pointer :: ts    => null()
    integer                     :: nline, satid

#define ERR(text) call finish('read_thresh_params',text)

    stat = 0
    call read_par_file(param_file_scatt,    sc_par)
    call read_par_file(param_file_si_amsua, sia_par)

  contains

    subroutine read_par_file(file, par)
      character(len=*), intent(inout)        :: file
      type(t_thresh),    intent(out),  target :: par

      type(t_thresh_par), pointer   :: p_p   => null()
      integer,           parameter :: n_str = 12
      character(len=300)           :: line, pstr
      character(len=300)           :: instr_str
      character(len=30)            :: str(n_str)
      integer                      :: instr, iu
      integer                      :: n, nfit, ndata, i
      integer                      :: version
      real(wp)                     :: lat

      if (file == '') RETURN
      if (dace%lpio) then
        stat = 13
        par%file_version = -99
        version          = -99
        file = path_file (data, file)
        iu = get_unit_number()
        write(*,*)
        write(*,*) 'Read scattering index test parameter file "'//trim(file)//'"'
        par%nset = 0
        open  (iu, file=trim(file), iostat=stat, action='read')
        do while (stat==0)
          read(iu,'(A)',iostat=stat) line
          if (is_iostat_end(stat)) exit
          if (stat/=0) ERR('Failed to read line')
          line = trim(adjustl(line))
          if (line == '') cycle
          if (line(1:1) == '#') cycle
          call split2(line,array=str,n=n,status=stat)
          if (stat/=0) ERR('Failed to split line: '//trim(line))
          if (n==0) cycle
          if (n >= 2) then
            select case(str(1))
            case('version')
              if (str(2) == '1') then
                version = 1
                stat = 0
              else
                ERR('file version '//trim(str(2))//' not supported.')
              end if
              cycle
            case('instr')
              if (par%nset >= 1) call sc_finish
              satid = -1
              read(str(2),*, iostat=stat) instr
              if (stat/=0) ERR('invalid instr '//trim(str(2)))
              if (n >= 4) then
                if (str(3) == 'sat') then
                  read(str(4),*, iostat=stat) satid
                  if (stat/=0) ERR('invalid satellite '//trim(str(4)))
                end if
              end if
              par%nset = par%nset + 1
              ts => par%set(par%nset)
              ts%instr       = instr
              ts%satid       = satid
              instr_str  = trim(line)
              nline      = 0
              cycle
            case('nlat','nlats')
              if (ts%nlat > 0) ERR('invalid (most like double) line "'//trim(line)//'"')
              read(str(2),*, iostat=stat) ts%nlat
              if (stat/=0) ERR('invalid nlat '//trim(str(2)))
              cycle
            end select
          end if
          if (n >= 4) then
            if (.not.associated(ts)) ERR('no instr (sat) before parameter line.')
            if (.not.associated(ts%p)) then
              if (ts%nlat <= 0)  ERR('no nlat for '//trim(instr_str))
              allocate(ts%p(ts%nlat))
              ts%npar = 0
            end if
            read(str(1),*, iostat=stat) lat
            if (stat/=0) ERR('invalid line "'//trim(line)//'" (lat)')
            read(str(2),*, iostat=stat) ndata
            if (stat/=0) ERR('invalid line "'//trim(line)//'" (ndata)')
            read(str(3),*, iostat=stat) nfit
            if (stat/=0) ERR('invalid line "'//trim(line)//'" (nfit)')
            if (n > 0 .and. nfit > 0) then
              ts%npar = ts%npar + 1
              if (ts%npar > size(ts%p)) ERR('parameter array too small')
              p_p => ts%p(ts%npar)
              p_p%lat   = lat
              p_p%ndata = ndata
              p_p%nfit  = nfit
              pstr     = ''
              do i = 4, n
                pstr = trim(pstr)//' '//trim(str(i))
              end do
              if (version == 1) pstr = trim(sat_poly2fparse(pstr))
              call init(pstr, c_var, p_p%p, stat)
              if (stat/=0) ERR('failed to initialize function : "'//trim(pstr)//'"')
            end if
            nline = nline + 1
          end if
        end do
        close (iu)
        call return_unit_number(iu)
        if (is_iostat_end(stat)) then
          par% file_version = version
          stat = 0
        end if
        if (par%nset >= 1) call sc_finish
        write(*,'(1x,"Read ",I2," parameter datasets.")') par%nset
      end if
      call p_bcast(par, dace%pio)
    end subroutine read_par_file

    subroutine sc_finish
      character(len=300) :: msg
      if (.not.associated(ts)) return
      if (satid > 0) then
        write(msg,'("instr ",I3,"  satid ",I3," :")') ts%instr, ts%satid
      else
        write(msg,'("instr ",I3," :")') ts%instr
      end if
      write(*,'(1x,A)') trim(msg)
      write(*,'(3x,I5," lines, ",I4," valid parameters ")') nline, ts%npar
      write(*,'(3x,I5," latitudes")') ts%nlat
    end subroutine sc_finish


  end subroutine read_thresh_params


  subroutine get_thresh_par(instr, satid, lat, p, thr, typ)
    integer,          intent(in)                   :: instr
    integer,          intent(in)                   :: satid
    real(kind=wp),    intent(in)                   :: lat
    type(t_range_fparse), intent(out)              :: p      ! Resulting function
    type(t_thresh),   intent(in), optional, target :: thr
    character(len=10),intent(in), optional         :: typ

    type(t_thresh_set), pointer :: ts => null()
    type(t_thresh),     pointer :: t  => null()
    integer                     :: iset, iset1, iset2, ip, i
    real(kind=wp)               :: md, d

    t => null()
    if (present(thr)) then
      t => thr
    elseif (present(typ)) then
      select case(typ)
      case('sc','scat','scatt')
        t => sc_par
      case('sia','si_amsua')
        t => sia_par
      case default
        call finish('get_thresh_par','typ="'//trim(typ)//'" not known')
      end select
    else
      call finish('get_thresh_par','Either typ or thr argument is required')
    end if

    if (t%file_version /= 1) RETURN
    iset1 = -1
    iset2 = -1
    do iset = 1, t%nset
      if (t%set(iset)%satid == satid .and. t%set(iset)%instr == instr) then
        iset1 = iset
      elseif (t%set(iset)%satid < 0 .and. t%set(iset)%instr == instr) then
        iset2 = iset
      end if
    end do
    if (iset1 > 0) then
      iset = iset1
    elseif (iset2 > 0) then
      iset = iset2  ! Backup set
    else
      iset = -1
    end if
    if (iset <= 0) RETURN
    ts => t%set(iset)
    md = huge(md)
    ip = -1
    do i = 1, ts%npar
      if (ts%p(i)%ndata > 20) then
        d = abs(lat - ts%p(i)%lat)
        if (d < md) then
          md = d
          ip = i
        end if
      end if
    end do
    if (ip > 0) then
      p = ts%p(ip)%p
    end if

  end subroutine get_thresh_par


  function get_threshold_val(thresh, s) result(x)
    real(wp)                            :: x
    type(t_range_fparse), intent(inout) :: thresh
    type(t_spot),         intent(inout) :: s
    real(wp) :: v(thresh%nvar)
    integer  :: ierr
    call get_tovs_var(thresh%vnames(1:thresh%nvar), v(:), ierr, s)
    if (ierr == 0) then
      x = evaluate(v(:),thresh)
    else
      x = inv_thr
    end if
  end function get_threshold_val


  function get_threshold_par(typ, instr, satid, s) result(x)
    real(wp)                        :: x
    character(len=*), intent(in)    :: typ
    integer,          intent(in)    :: instr
    integer,          intent(in)    :: satid
    type(t_spot),     intent(inout) :: s

    type(t_range_fparse) :: thresh_

    call get_thresh_par(instr, satid, s%col%c%dlat, thresh_, typ=typ)
    if (thresh_%func /= '') then
      x = get_threshold_val(thresh_, s)
    else
      x = inv_thr
    end if

  end function get_threshold_par


  subroutine read_param_amsub (file)
    ! Description:
    ! Reads the AMSU-B Channel 18 (183.3 +/- 1 GHz) clear-sky thresholds file.
    character(len=*) ,intent(in)   :: file  ! file name
    type(t_buehler_mod), pointer   :: p_bm => null()
    type(t_buehler_par), pointer   :: p_p  => null()
    character(len=120)             :: msg  = ''
    integer,             parameter :: n_str = 12
    character(len=300)             :: file_
    character(len=300)             :: line
    character(len=300)             :: instr_str
    character(len=30)              :: str(n_str)
    integer                        :: n, stat, instr, satid, fov, nline, iunit
    real(wp)                       :: lat, stzen, bth, a, b
    !--------------------------------------
    ! read AMSUB cloud detection parameters
    !--------------------------------------
#define ERR(text) call finish('read_param_amsub',text)

    if (file == '') RETURN

    if (dace%lpio) then
      file_ = path_file (data, file)
      iunit = get_unit_number()
      param_amsub% file_version = -99
      write(*,*)
      write(*,*) 'Read cloud detection parameter file "'//trim(file_)//'" for mhs-like instruments.'
      param_amsub%n_bm = 0
      open  (iunit, file=file_, iostat=stat, action='read')
      if (stat /= 0) ERR('Failed to read amsub cloud det. parameter file "'//trim(file_)//'"')
      read  (iunit,*, iostat=stat) param_amsub% peak_ch_18
      if (stat /= 0) then
        close (iunit)
        open  (iunit, file=file_, iostat=stat, action='read')
        do while(stat==0)
          read(iunit,'(A)',iostat=stat) line
          if (is_iostat_end(stat)) exit
          if (stat/=0) ERR('Failed to read line')
          line = trim(adjustl(line))
          if (line == '') cycle
          if (line(1:1) == '#') cycle
          call split2(line,array=str,n=n,status=stat)
          if (stat/=0) ERR('Failed to split line: '//trim(line))
          if (n==0) cycle
          if (n >= 2) then
            select case(str(1))
            case('version')
              if (str(2) == '1') then
                param_amsub% file_version = 1
                write(*,*) 'File format version 1'
              else
                ERR('file version '//trim(str(2))//' not supported.')
              end if
            case('instr')
              if (param_amsub%n_bm >= 1) call bm_finish
              satid = -1
              read(str(2),*, iostat=stat) instr
              if (stat/=0) ERR('invalid instr '//trim(str(2)))
              if (n >= 4) then
                if (str(3) == 'sat') then
                  read(str(4),*, iostat=stat) satid
                  if (stat/=0) ERR('invalid satellite '//trim(str(4)))
                end if
              end if
              param_amsub%n_bm = param_amsub%n_bm + 1
              p_bm => param_amsub%bm(param_amsub%n_bm)
              p_bm%instr       = instr
              p_bm%satid       = satid
              p_bm%dlat_thresh = dlat_thresh_def
              p_bm%dsz_thresh  = dsz_thresh_def
              instr_str  = trim(line)
              nline      = 0
            case('nlat','nlats')
              if (p_bm%nlat > 0) ERR('invalid (most like double) line "'//trim(line)//'"')
              read(str(2),*, iostat=stat) p_bm%nlat
              if (stat/=0) ERR('invalid nlat '//trim(str(2)))
            case('nfov','nfovs')
              if (p_bm%nfov > 0) ERR('invalid (most like double) line "'//trim(line)//'"')
              read(str(2),*, iostat=stat) p_bm%nfov
              if (stat/=0) ERR('invalid nfov '//trim(str(2)))
            case('dlat_thresh')
              read(str(2),*, iostat=stat) p_bm%dlat_thresh
              if (stat/=0) ERR('invalid dlat_thresh '//trim(str(2)))
            case('dsz_thresh')
              read(str(2),*, iostat=stat) p_bm%dsz_thresh
              if (stat/=0) ERR('invalid dsz_thresh '//trim(str(2)))
            end select
          end if
          if (n >= 7) then
            if (.not.associated(p_bm)) ERR('no instr (sat) before parameter line.')
            if (.not.associated(p_bm%p)) then
              if (p_bm%nlat <= 0)  ERR('no nlat for '//trim(instr_str))
              if (p_bm%nfov <= 0)  ERR('no nfov for '//trim(instr_str))
              allocate(p_bm%p   (p_bm%nlat * p_bm%nfov))
              p_bm%npar = 0
            end if
            read(line,*, iostat=stat) lat, fov, stzen, n, bth, a, b
            if (stat/=0) ERR('invalid line "'//trim(line)//'"')
            if (abs(stzen) < 90._wp) then
              p_bm%npar = p_bm%npar + 1
              if (p_bm%npar > size(p_bm%p)) ERR('parameter array too small')
              p_p => p_bm%p(p_bm%npar)
              p_p%lat      = lat
              p_p%stzen    = stzen
              p_p%n        = n
              p_p%bt_h_min = bth
              p_p%a        = a
              p_p%b        = b
            end if
            nline = nline + 1
          end if
        end do
        if (param_amsub%n_bm >= 1) call bm_finish
        write(*,'(1x,"Read ",I2," parameter datasets.")') param_amsub%n_bm
        close (iunit)
      else
        write(*,*) 'File format version 0'
        param_amsub% file_version = 0
      end if

      if (param_amsub% file_version < 0 .or. param_amsub% file_version > 1) &
           ERR('failed to read file (unknown version)')
    end if

  contains

    subroutine bm_finish
      integer       :: ncorr, ipar
      real(kind=wp) :: y

      if (.not.associated(p_bm)) return
      if (satid > 0) then
        write(msg,'("instr ",I3,"  satid ",I3," :")') p_bm%instr, p_bm%satid
      else
        write(msg,'("instr ",I3," :")') p_bm%instr
      end if
      write(*,'(1x,A)') trim(msg)
      write(*,'(3x,I4," lines, ",I4," valid parameters ")') nline, p_bm%npar
      write(*,'(3x,I4," latitudes, ",I3," FOVs ")') p_bm%nlat, p_bm%nfov
      ! Correct unphysical parameters
      ncorr = 0
      do ipar = 1, p_bm%npar
        p_p => p_bm%p(ipar)
        if (p_p%b > max_slope) then
          y = p_p%a + p_p%bt_h_min * p_p%b
          p_p%a = y - p_p%bt_h_min * max_slope
          p_p%b = max_slope
          ncorr = ncorr + 1
        end if
      end do
      if (ncorr > 0) write(*,'(3x,"unphysical (corrected) parameters: ",F6.2,"%")') &
           (100.*ncorr)/p_bm%npar
      ! Inspect sat.zen. angles
      p_bm%sz_mode = 0
      if (p_bm%npar > 0) then
        if (.not.all(p_bm%p(1:p_bm%npar)%stzen >= 0._wp)) p_bm%sz_mode = 1
      end if
      msg = 'Sat. zenith angles: '
      select case(p_bm%sz_mode)
      case(0)
        msg = trim(msg)//' positive'
      case(1)
        msg = trim(msg)//' positive and negative'
      end select
      write(*,'(3x,A)') trim(msg)
      ! Possibly we have to adapt the maximum acceptable diffs. for lat and stzen:
      if (p_bm%nlat > 0) then
        p_bm%dlat_thresh = max(180._wp/(p_bm%nlat*lat_scale), p_bm%dlat_thresh)
      end if
      if (p_bm%nfov > 0 .and. p_bm%npar > 0) then
        a = min(maxval(abs(p_bm%p(1:p_bm%npar)%stzen)), 87.5_wp)
        b = a * (1._wp - 1._wp/(2*p_bm%nfov))
        p_bm%dsz_thresh = max(1._wp/cos(a*d2r) - 1._wp/cos(b*d2r), p_bm%dsz_thresh)
      end if
      write(*,'(3x,"maximum accepted scaled differences: dlat_thresh=",F7.4,"  dsz_thresh=",F7.4)') &
           p_bm%dlat_thresh, p_bm%dsz_thresh
    end subroutine bm_finish

  end subroutine read_param_amsub


  subroutine get_threshold_buehler(stzen, lat, obs_high, satid, instr, bt_h_min, hl_max, ierr, ld)
    real(wp), intent(in)  :: stzen
    real(wp), intent(in)  :: lat
    real(wp), intent(in)  :: obs_high
    integer,  intent(in)  :: satid
    integer,  intent(in)  :: instr
    real(wp), intent(out) :: bt_h_min
    real(wp), intent(out) :: hl_max
    integer,  intent(out) :: ierr
    logical,  intent(in), optional :: ld

    real(wp),            parameter :: lat_band = 5        ! size of latitude band [degree]
    real(wp),            parameter :: lat_low_bound = -90 ! lowest latitude bound [degree]
    type(t_buehler_mod), pointer   :: p_bm => null()
    type(t_buehler_par), pointer   :: p_p  => null()
    integer                        :: fov_indx, lat_indx
    integer                        :: ibm1, ibm2, ipar
    integer                        :: i
    logical                        :: ld_

    if (present(ld)) then ; ld_=ld ; else ; ld_=.false. ; endif

    bt_h_min = inv_thr
    hl_max   = inv_thr
    ierr = 0

    select case(param_amsub% file_version)
    case(0)
      fov_indx = minloc(abs(stzen - amsub_zen), 1)
      lat_indx = int((lat - lat_low_bound + lat_band) / lat_band)
      if (lat_indx == 0)   lat_indx = 1
      if (lat_indx > nlat) lat_indx = nlat
      bt_h_min = param_amsub% peak_ch_18(fov_indx, lat_indx)
      hl_max   = param_amsub% delta_ch20_ch18
    case (1)
      ! Search parameters
      ibm1 = -1
      ibm2 = -1
      do i = 1, param_amsub% n_bm
        p_bm => param_amsub% bm(i)
        if (p_bm%instr == instr .and. p_bm%satid == satid) ibm1 = i
        if (p_bm%instr == instr                          ) ibm2 = i
      end do
      !if (ld_) print*,'get_threshold imb',ibm1,ibm2
      if (max(ibm1, ibm2) < 0) then
        ierr = 1
        return
      end if
      ipar = -1
      if (ibm1 > 0) then
        p_bm => param_amsub%bm(ibm1)
        call find_par(ipar)
      end if
      !if (ld_) print*,'get_threshold ipar',ipar
      if (ipar <= 0) then
        p_bm => param_amsub%bm(ibm2)
        call find_par(ipar)
      end if
      !if (ld_) print*,'get_threshold ipar',ipar
      if (ipar <= 0) then
        ierr = 1
        return
      end if
      p_p  => p_bm%p(ipar)

      bt_h_min = p_p%bt_h_min
      hl_max   = p_p%a + p_p%b*obs_high
      !if (ld_) print*,'get_threshold',bt_h_min,hl_max
    case default
      ierr = 2
    end select

    if (buehler_version >= 1) then
      if (bt_h_min <= 0._wp) bt_h_min = inv_thr
    end if

  contains

    subroutine find_par(ipar)
      integer,       intent(out) :: ipar

      real(kind=wp) :: dlat, dsz, cost, cost_min
      integer       :: i, np

      ipar     = -1
      np       = p_bm%npar
      cost_min = huge(cost_min)
      do i = 1, np
        p_p => p_bm%p(i)
        dlat = abs(lat-p_p%lat) / lat_scale
        if (dlat > p_bm%dlat_thresh) cycle
        select case(p_bm%sz_mode)
        case(0)
          dsz  = abs(1._wp/cos(stzen*d2r) - 1._wp/cos(p_p%stzen*d2r))
        case(1)
          if (sign(1._wp,stzen) /= sign(1._wp,p_p%stzen)) then
            dsz  = abs(1._wp/cos(stzen*d2r) + 1._wp/cos(p_p%stzen*d2r) - 2._wp)
          else
            dsz  = abs(1._wp/cos(stzen*d2r) - 1._wp/cos(p_p%stzen*d2r))
          end if
        end select
        dsz = dsz  / sz_scale
        if (dsz > p_bm%dsz_thresh) cycle
        cost = dlat**2 * dsz**2
        if (cost < cost_min) then
          ipar = i
          cost_min = cost
        end if
      end do

    end subroutine find_par

  end subroutine get_threshold_buehler


  subroutine bcast_param_amsub (s_p, source, comm)
    type (t_param_amsub) ,intent(inout) :: s_p
    integer              ,intent(in)    :: source
    integer    ,optional ,intent(in)    :: comm

    integer :: i, j
    integer :: lcom
    integer :: count = 0

    lcom = dace% comm ;if(present(comm)) lcom = comm

    if (count == 0) count = size(transfer(s_p,(/' '/)))

    if (dace% pe /= source) call destruct(s_p%bm)
    call p_bcast_derivedtype (s_p, count, source, lcom)
    if (dace% pe /= source) then
      do i = 1, s_p%n_bm
        if (associated(s_p%bm(i)%p)) allocate(s_p%bm(i)%p(s_p%bm(i)%npar))
      end do
    end if

    do i = 1, s_p%n_bm
      if (associated(s_p%bm(i)%p)) then
        do j = 1, s_p%bm(i)%npar
          call p_bcast(s_p%bm(i)%p(j), source, lcom)
        end do
      end if
    end do

  end subroutine bcast_param_amsub


  subroutine bcast_thresh (t, source, comm)
    type (t_thresh),   intent(inout) :: t
    integer,           intent(in)    :: source
    integer, optional, intent(in)    :: comm

    integer :: i, j
    integer :: lcom
    integer :: count = 0

    lcom = dace% comm ;if(present(comm)) lcom = comm

    if (count == 0) count = size(transfer(t,(/' '/)))

    if (dace% pe /= source) call destruct(t%set)
    call p_bcast_derivedtype(t, count, source, lcom)
    if (dace% pe /= source) then
      do i = 1, t%nset
        if (associated(t%set(i)%p)) allocate(t%set(i)%p(t%set(i)%npar))
      end do
    end if
    do i = 1, t%nset
      if (associated(t%set(i)%p)) then
        do j = 1, t%set(i)%npar
          call p_bcast(t%set(i)%p(j), source, lcom)
        end do
      end if
    end do

  end subroutine bcast_thresh

  subroutine bcast_thresh_par (tp, source, comm)
    type (t_thresh_par), intent(inout) :: tp
    integer,             intent(in)    :: source
    integer,  optional,  intent(in)    :: comm
    call p_bcast(tp%lat  , source, comm)
    call p_bcast(tp%ndata, source, comm)
    call p_bcast(tp%nfit , source, comm)
    call p_bcast(tp%p    , source, comm)
  end subroutine bcast_thresh_par

  elemental subroutine construct_thresh_set(s)
    type(t_thresh_set), intent(out) :: s
  end subroutine construct_thresh_set

  elemental subroutine destruct_thresh_set(s)
    type(t_thresh_set), intent(inout) :: s
    if (associated(s%p)) deallocate(s%p)
    call construct(s)
  end subroutine destruct_thresh_set

  elemental subroutine construct_thresh_par(p)
    type(t_thresh_par), intent(out) :: p
  end subroutine construct_thresh_par

  elemental subroutine destruct_thresh_par(p)
    type(t_thresh_par), intent(inout) :: p
    call destruct(p%p)
    call construct(p)
  end subroutine destruct_thresh_par


  elemental subroutine construct_t_buehler_mod(x)
    type(t_buehler_mod), intent(out) :: x
  end subroutine construct_t_buehler_mod

  elemental subroutine destruct_t_buehler_mod(x)
    type(t_buehler_mod), intent(inout) :: x
    if (associated(x%p)) deallocate(x%p)
    call construct(x)
  end subroutine destruct_t_buehler_mod

!
! specific MPI-bcast for derived type t_param_amsub
!
#undef  VECTOR
#undef  DERIVED
#define DERIVED type(t_buehler_par)
#define p_bcast_DERIVED bcast_buehler_par
#undef  MPI_TYPE
#include "p_bcast.incf"
!==============================================================================
end module mo_cloud_params
