!
!+ global arrays to read observation header from NetCDF
!
MODULE mo_head_netcdf
!
! Description:
!   Declaration of global arrays to hold temporary data for
!   reading observations header data from NetCDF.
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
! V1_8         2009/12/09 Harald Anlauf
!  Add "implicit none"
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  Move dimids, dimids_max to mo_head_netcdf
! V1_31        2014-08-21 Andreas Rhodin
!  allow up to 30 dimensions for (BUFR2)NetCDF file input
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Gerhard Paul   DWD  2008  original code
! Andreas Rhodin DWD  2008
! Harald Anlauf  DWD  2010
!------------------------------------------------------------------------------

!=============
! Modules used
!=============
  use mo_kind,          only: sp, dp, i8     ! kind parameters
  use mo_mpi_dace,      only: p_bcast        ! generic MPI-bcast routine
  use mo_exception,     only: finish         ! abort in case of error
  use mo_time,          only: t_time,       &! time data type
                              init_time      ! initialise time data type
  use mo_wigos,         only: t_wsi          ! WIGOS id data type
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,         only: nf90_Inquire_Variable, &!
                            nf90_inq_varid,        &!
                            nf90_get_var,          &!
                            nf90_get_att,          &!
                            nf90_strerror,         &!
                            NF90_FLOAT,            &!
                            NF90_DOUBLE,           &!
                            NF90_INT,              &!
                            NF90_NOERR              !
  implicit none

!================
! public entities
!================
  private
  public :: ncid       ! NetCDF file id
  public :: dimids_max ! max. number of NetCDF dimensiond
  public :: imissing   ! NetCDF _FillValue for integers
  public :: rmissing   ! NetCDF _FillValue for reals
  public :: ymissing   ! NetCDF _FillValue for character
  public :: iascii_ymissing ! number in ascii sorting sequence for ymissing value
  public :: s2ikz      ! DWD-internal classifier
  public :: edition    ! BUFR edition
  public :: s1date     ! synoptic date (section 1)
  public :: s1time     ! synoptic time (section 1)
  public :: s1cat      ! data category
  public :: s1cats     ! data sub-category
  public :: s1catls    ! local data sub-category
  public :: s1cent     ! data centre
  public :: s1cents    ! data subcentre
  public :: stime      ! header observation time (section1)
  public :: db_time    ! data bank time
  public :: s1updat    ! update sequence number
  public :: s2dcdate   ! data base decoding date
  public :: s2dctime   ! data base decoding time
  public :: mlah       ! latitude
  public :: mloh       ! longitude
  public :: obs_time   ! valid time of observation
  public :: mjjj       !   jear
  public :: mmm        !   month
  public :: myy        !   day
  public :: mgg        !   hour
  public :: ngg        !   minute
  public :: msec       !   second
  public :: istidn     ! station id as integer
  public :: ystidn     ! station id as character string
  public :: ilstidn    ! length of station identifier in characters
  public :: ystid      ! any type of station identifier as array
  public :: mii        ! WMO block number
  public :: niii       ! WMO station number
  public :: istidb     ! WMO buoy/platform identifier
  public :: istids     ! WMO satellite identifier
  public :: msec_r, bdy_date, bdy_time, &
            lvarid     ! needed in mo_obs_netcdf
  public :: lwsi       ! WIGOS station identifier provided
  public :: wsi        ! WIGOS station id stored representation
  public :: wsihash    ! Station id hash

  public :: get_int    ! read integer variable from NetCDF file
  public :: get_real   ! read real    variable from NetCDF file

!==================
! Module variables:
!==================
  !----------------------------
  ! general info on NetCDF file
  !----------------------------
  integer , PARAMETER       :: dimids_max=30         ! max number of dimension ids
  integer                   :: ncid                  ! Netcdf file id
  !---------------------------------------------------------------
  ! variables for exchange of Header information from NetCDF files
  !---------------------------------------------------------------
  integer ,allocatable      :: edition (:)  ! BUFR edition
  !=================================
  ! Header information section 1 / 2
  !=================================
  !-------------------------------------------------
  !    expansion  in NetCDF file BUFR- data section1
  !    ftn-fields                for   data section1  (date time category centre)
  !-------------------------------------------------
  integer ,allocatable      :: s1date  (:)  ! nominal (synoptic) date
  integer ,allocatable      :: s1time  (:)  ! nominal (synoptic) time
  integer ,allocatable      :: s1cat   (:)  ! data category
  integer ,allocatable      :: s1cats  (:)  ! data sub category
  integer ,allocatable      :: s1catls (:)  ! local data sub category
  integer ,allocatable      :: s1cent  (:)  ! data centre
  integer ,allocatable      :: s1cents (:)  ! data sub centre
  integer ,allocatable      :: s1updat (:)  ! update sequence no.
  !-------------------------------------------------
  !    expansion  in NetCDF file BUFR- data section2
  !    ftn-fields                for   data section2 (kennzahl decoding date time)
  !-------------------------------------------------
  integer ,allocatable      :: s2ikz   (:)  ! DWD-internal classifier
  integer ,allocatable      :: s2dcdate(:)  ! decoding date
  integer ,allocatable      :: s2dctime(:)  ! decoding time
  type (t_time),allocatable :: stime   (:)  ! header observation time (section1)
  type (t_time),allocatable :: db_time (:)  ! data bank time
  !===================
  ! body   information
  !===================
  !-------------------------------------------------
  !    expansion  in NetCDF file BUFR- data section4
  !    ftn-fields                for   data section4 (space time)
  ! for exchange
  !-------------------------------------------------
  real    ,allocatable      :: mlah    (:)  ! latitude
  real    ,allocatable      :: mloh    (:)  ! longitude
  type (t_time),allocatable :: obs_time(:)  ! body observation time
  !--------------
  ! for expansion
  !--------------
  integer ,allocatable      :: mjjj    (:)   ! year
  integer ,allocatable      :: mmm     (:)   ! month
  integer ,allocatable      :: myy     (:)   ! day
  integer ,allocatable      :: mgg     (:)   ! hour
  integer ,allocatable      :: ngg     (:)   ! minute
  integer ,allocatable      :: msec    (:)   ! second ; most observation types
  real    ,allocatable      :: msec_r  (:)   ! second ; radiooccultations in SSsss
  !---------------------
  !  combined  date time
  !---------------------
  integer ,allocatable      :: bdy_date(:) ! observation date: CCYYMMDD
  integer ,allocatable      :: bdy_time(:) ! observation time: HHMMSSsss; sec.123
  !
  !------------------------------
  ! station identifier (combined)
  ! for exchange
  !------------------------------
  integer ,allocatable      :: istidn (:)   ! WMO numeric station number combined
  integer , PARAMETER       :: ilstidn=10   ! length of station identifier in characters
  character (LEN=ilstidn),allocatable                                          &
                            :: ystidn(:)    ! any type of station identifier as variable
  !--------------
  ! for expansion
  !--------------
  integer ,allocatable      :: mii    (:)   ! WMO block number
  integer ,allocatable      :: niii   (:)   ! WMO station number
  integer ,allocatable      :: istidb (:)   ! WMO buoy/platform identifier
  integer ,allocatable      :: istids (:)   ! WMO satellite identifier
! character ,allocatable    :: ydddd(:,:)   ! WMO temp/pilot mobile station identifier
  character(LEN=1) ,allocatable ::                                             &
                               ystid(:,:)   ! any type of station identifier as array
  logical,allocatable       :: lvarid (:)   ! variable for station ID exists in NetCDF file

  !-------------------------
  ! WIGOS station identifier
  !-------------------------
  logical      ,allocatable :: lwsi   (:)   ! WIGOS station identifier valid
  type(t_wsi)  ,allocatable :: wsi    (:)   ! WIGOS station id stored representation
  integer(i8)  ,allocatable :: wsihash(:)   ! Station id hash

  !=======================================
  ! NetCDF variable IDs for header section
  !=======================================
  !-------------------------------------------------
  ! variable ID's in NetCDF file for
  ! no expansion  in NetCDF file BUFR- data section0 (total length in Bytes)
  !    expansion  in NetCDF file BUFR- data section1
  ! variable ID's in NetCDF file in    data section1  (date time category centre)
  !-------------------------------------------------
! integer              :: varid_s1date   ! nominal (synoptic) date
! integer              :: varid_s1time   ! nominal (synoptic) time
! integer              :: varid_s1cat    ! data category
! integer              :: varid_s1cats   ! data sub category  (undefined in bufr edition 3)
! integer              :: varid_s1catls  ! local data sub category
! integer              :: varid_s1cent   ! data centre
! integer              :: varid_s1cents  ! data sub centre
! integer              :: varid_s1updat  ! update sequence no.
  !-------------------------------------------------
  !    expansion  in NetCDF file BUFR- data section2
  ! variable ID's in NetCDF file  in   data section2 (kennzahl decoding date time)
  !-------------------------------------------------
! integer              :: varid_s2ikz    ! DWD-internal classifier
! integer              :: varid_s2dcdate ! decoding date
! integer              :: varid_s2dctime ! decoding time
  !-------------------------------------------------
  ! variable ID's in NetCDF file within data section3 (number of subsets,
  !                                                    still not defined in netCDF)
  ! variable ID's in NetCDF file within data section4
  !-------------------------------------------------
! integer              :: varid_mii      ! NetCDF variable  id for  WMO-Block-Nummer
! integer              :: varid_niii     ! NetCDF variable  id for  WMO-Stationsnummer
! integer              :: varid_mabnn    ! NetCDF variable  id for  WMO buoy/platform identifier
! integer              :: varid_mi1i2    ! NetCDF variable  id for  WMO satellite identifier
! integer              :: varid_ystid    ! NetCDF variable  id for  any type of character station id
! integer              :: varid_mlah     ! NetCDF variable  id for  Latitude
! integer              :: varid_mloh     ! NetCDF variable  id for  Longitude
! integer              :: varid_mjjj     ! NetCDF variable  id for  year
! integer              :: varid_mmm      ! NetCDF variable  id for  month
! integer              :: varid_myy      ! NetCDF variable  id for  day
! integer              :: varid_mgg      ! NetCDF variable  id for  hour
! integer              :: varid_ngg      ! NetCDF variable  id for  minute
! integer              :: varid_msec     ! NetCDF variable  id for  second
  !=======================================================
  ! missing values in NetCDF file for different data types
  !=======================================================
  integer              :: imissing       ! NetCDF _FillValue for integer
  real                 :: rmissing       ! NetCDF _FillValue for reals
  character (LEN=1)    :: ymissing       ! NetCDF _FillValue for character
  integer              :: iascii_ymissing! number in ascii sorting sequence for ymissing value
!------------------------------------------------------------------------------
  !===========
  ! Interfaces
  !===========
  interface get_real
    module procedure get_real     ! Read real variable from NetCDF file (1d)
    module procedure get_real_2d  ! Read real variable from NetCDF file (2d)
  end interface get_real

  interface get_int
     module procedure get_int     ! Read integer variable from NetCDF file (1d)
     module procedure get_int_2d  ! Read integer variable from NetCDF file (1d)
  end interface get_int
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  subroutine get_int (x, name, miss, count, start)
    integer          ,intent(out)          :: x(:)
    character(len=*) ,intent(in)           :: name
    integer          ,intent(in)           :: miss
    integer          ,intent(in) ,optional :: count(:)
    integer          ,intent(in) ,optional :: start(:)
    !--------------------------------
    ! read variable, store to integer
    !--------------------------------
      integer :: status
      integer :: varid
      !--------------
      ! read variable
      !--------------
      status = nf90_inq_varid (ncid, name,  varid)
      if (status /= nf90_noerr) &
        call finish('netcdf:get_int', 'cannot read '//trim(name)//' : '//&
                    trim(nf90_strerror(status))                          )
      status = nf90_get_var (ncid, varid, x, count=count, start=start)
      if (status /= nf90_noerr) &
        call finish('netcdf:get_int', 'cannot read '//trim(name)//' : '//&
                    trim(nf90_strerror(status))                          )
      !-------------------------
      ! check for missing values
      !-------------------------
      where (x == imissing) x = miss
  end subroutine get_int

!------------------------------------------------------------------------------

  subroutine get_int_2d (x, name, miss, count, start)
    integer          ,intent(out)          :: x(:,:)
    character(len=*) ,intent(in)           :: name
    integer          ,intent(in)           :: miss
    integer          ,intent(in) ,optional :: count(:)
    integer          ,intent(in) ,optional :: start(:)
    !--------------------------------
    ! read variable, store to integer
    !--------------------------------
      integer :: status
      integer :: varid
      !--------------
      ! read variable
      !--------------
      status = nf90_inq_varid (ncid, name,  varid)
      if (status /= nf90_noerr) &
        call finish('netcdf:get_int', 'cannot read '//trim(name)//' : '//&
                    trim(nf90_strerror(status))                          )
      status = nf90_get_var (ncid, varid, x, count=count, start=start)
      if (status /= nf90_noerr) &
        call finish('netcdf:get_int', 'cannot read '//trim(name)//' : '//&
                    trim(nf90_strerror(status))                          )
      !-------------------------
      ! check for missing values
      !-------------------------
      where (x == imissing) x = miss
  end subroutine get_int_2d

!------------------------------------------------------------------------------

  subroutine get_real (x, name, miss, count, start)
    real(sp)         ,intent(out)          :: x(:)
    character(len=*) ,intent(in)           :: name
    real(sp)         ,intent(in)           :: miss
    integer          ,intent(in) ,optional :: count(:)
    integer          ,intent(in) ,optional :: start(:)
    !-----------------------------
    ! read variable, store to real
    !-----------------------------
      integer              :: status
      integer              :: varid
      integer              :: xtype ! NetCDF variable type
      integer ,allocatable :: ix(:) ! temporary
      real(dp),allocatable :: dx(:) ! temporary
      real(dp)             :: dmiss
      integer              :: n
      !-----------------------------------
      ! check for data type in NetCDF file
      !-----------------------------------
      status = nf90_inq_varid (ncid, name, varid)
      if (status /= nf90_noerr) &
        call finish('netcdf:get_real', 'cannot read '//trim(name)//' : '//&
                    trim(nf90_strerror(status))                           )
      status = nf90_inquire_variable (ncid, varid, xtype=xtype)
      !----------------------------------------
      ! read variable, check for missing values
      !----------------------------------------
      n = size (x)
      select case (xtype)
      case (NF90_FLOAT)
        status = nf90_get_var (ncid, varid, x, count=count, start=start)
        if (status /= nf90_noerr) &
          call finish('netcdf:get_real', 'cannot read '//trim(name)//' : '//&
                      trim(nf90_strerror(status))                           )
        where (x == rmissing) x = miss
      case (NF90_DOUBLE)
        allocate (dx (n))
        status = nf90_get_var (ncid, varid, dx, count=count, start=start)
        if (status /= nf90_noerr) &
          call finish('netcdf:get_real', 'cannot read '//trim(name)//' : '//&
                      trim(nf90_strerror(status))                           )
        status = nf90_get_att (ncid, varid , '_FillValue', dmiss)
        if (status /= nf90_noerr) &
          call finish('netcdf:get_real', 'fillvalue?? '//trim(name)//' : '//&
                      trim(nf90_strerror(status))                           )
        where (dx == dmiss)
          x = miss
        elsewhere
          x = dx
        endwhere
      case (NF90_INT)
        allocate (ix (n))
        status = nf90_get_var (ncid, varid, ix, count=count, start=start)
        if (status /= nf90_noerr) &
          call finish('netcdf:get_real', 'cannot read '//trim(name)//' : '//&
                      trim(nf90_strerror(status))                           )
        where (ix == imissing)
          x = miss
        elsewhere
          x = ix
        endwhere
      case default
        call finish('netcdf:get_real', 'invalid xtype for '//trim(name))
      end select
  end subroutine get_real

!------------------------------------------------------------------------------
  subroutine get_real_2d (x, name, miss, count, start)
    real(sp)         ,intent(out)          :: x(:,:)
    character(len=*) ,intent(in)           :: name
    real(sp)         ,intent(in)           :: miss
    integer          ,intent(in) ,optional :: count(:)
    integer          ,intent(in) ,optional :: start(:)
    !-----------------------------
    ! read variable, store to real
    !-----------------------------
    integer              :: status
    integer              :: varid
    integer              :: xtype   ! NetCDF variable type
    integer              :: n(2)    ! Shape of variable
    integer ,allocatable :: ix(:,:) ! temporary
    real(dp),allocatable :: dx(:,:) ! temporary
    real(dp)             :: dmiss
    !-----------------------------------
    ! check for data type in NetCDF file
    !-----------------------------------
    status = nf90_inq_varid (ncid, name, varid)
    if (status /= nf90_noerr) goto 999
    status = nf90_inquire_variable (ncid, varid, xtype=xtype)
    !----------------------------------------
    ! read variable, check for missing values
    !----------------------------------------
    n = shape (x)
    select case (xtype)
    case (NF90_FLOAT)
       status = nf90_get_var (ncid, varid, x,  start=start, count=count)
       if (status /= nf90_noerr) goto 999
       where (x == rmissing) x = miss
    case (NF90_INT)
       allocate (ix (n(1),n(2)))
       status = nf90_get_var (ncid, varid, ix, start=start, count=count)
       if (status /= nf90_noerr) goto 999
       where (ix == imissing)
          x = miss
       elsewhere
          x = ix
       endwhere
    case (NF90_DOUBLE)
       allocate (dx (n(1),n(2)))
       status = nf90_get_var (ncid, varid, dx, start=start, count=count)
       if (status /= nf90_noerr) goto 999
       status = nf90_get_att (ncid, varid , '_FillValue', dmiss)
       if (status /= nf90_noerr) &
            call finish ('get_real_2d', 'fillvalue?? '//trim(name)//' : '//&
                         trim (nf90_strerror(status))                      )
       where (dx == dmiss)
          x = miss
       elsewhere
          x = dx
       endwhere
    case default
       call finish ('get_real_2d', 'invalid xtype for '//trim(name))
    end select
    return
999 continue
    if (status /= nf90_noerr) &
         call finish ('get_real_2d', 'cannot read '//trim(name)//' : '//&
                      trim (nf90_strerror(status))                      )
  end subroutine get_real_2d
!------------------------------------------------------------------------------
end module mo_head_netcdf
