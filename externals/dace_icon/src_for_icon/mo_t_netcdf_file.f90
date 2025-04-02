!
!+ 3DVAR/COSMO derived data type which mirrors meta data struct. of NetCDF files
!
!-------------------------------------------------------------------------------
!
MODULE mo_t_netcdf_file
!
!-------------------------------------------------------------------------------
! Description:
!   This module defines the derived type 't_netcdf_file' which
!   mirrors some of the meta data of a NetCDF-file. The subroutines
!   'add_dim' and 'add_var' add a NetCDF-dimension or a
!   NetCDF-variable to a variable of this type. Subroutine
!   'create_netcdf_file' creates (writes) an empty (dimension and
!   variable definitions) NetCDF file based on the information in the
!   variable. Subroutines to write the data (contents of the variables)
!   are not yet implemented.) Subroutine 'close_netcdf_file'
!   closes the NetCDF file. Subroutine 'destruct_netcdf_file'
!   finally deallocates pointer components of a variable of type
!   't_netcdf_file' when it is not used any more.
!   This module is jointly used by the COSMO model and 3DVAR program packages
!
! Current Code Owners:
!    For DWD 3DVAR:                        For COSMO:
!    DWD, Harald Anlauf                    DWD, Christoph Schraff
!    phone: +49 69 8062 4941               phone: +49 69 8062 2725
!    fax:   +49 69 8062 3721               fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de           email: christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  New routines: open_netcdf_file_read/write
! V1_5         2009/05/25 Andreas Rhodin
!  handle FillValues
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Christoph Schraff
!  modifications required for COSMO to write feedobs files
! V1_17        2011/12/21 Andreas Rhodin
!  extend error message in case of file creation failure
! V1_19        2012-04-16 Andreas Rhodin
!  create netCDF files with option NF_64BIT_OFFSET
!  increase size of optionals array (ODIM) to 8 (required for plevel)
! V1_20        2012-06-18 Andreas Rhodin
!  force zero sized dimensions to 1 (0 interpreted as unlimited dimension)
! V1_31        2014-08-21 Andreas Rhodin
!  increase PLEN (len of file path/name)
! V1_37        2014-12-23 Andreas Rhodin
!  homogenize interface with COSMO V5.1
! V1_51        2017-02-24 Andreas Rhodin
!  update shared files with COSMO 5.04d
!
! CAUTION: This module is used commonly by the 3DVAR and COSMO main programs.!!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2007  original code to be used in 3DVAR/LETKF/COSMO
!--------------------------------------------------------------------------

!=============
! Modules used
!=============
use kind_parameters, only: dp       ! real precision kind parameters
use mo_exception,    only: finish   ! abort in case of error
use mo_t_table,      only: t_table  ! derived types for tables
use mo_netcdf_param                 ! include 'netcdf.inc'

implicit none

!================
! public entities
!================
private
!-------------------------
! derived type definitions
!-------------------------
public :: t_netcdf_file          ! derived type, structure of NetCDF files
public :: t_netcdf_var           !   component of t_netcdf_file (variables)
public :: t_netcdf_dim           !   component of t_netcdf_file (dimensions)
public :: p_netcdf_dim           !   derived type, pointer to t_netcdf_dim
!------------
! subroutines
!------------
public :: add_dim                ! add a variable         to t_netcdf_file
public :: add_var                ! add a dimension        to t_netcdf_file
public :: create_netcdf_file     ! create a netcdf_file from t_netcdf_file
public :: close_netcdf_file      ! close  a netcdf_file
public :: destruct_netcdf_file   ! deallocate components  of t_netcdf_file
public :: open_netcdf_file_read  ! open file for read       access
public :: open_netcdf_file_write ! open file for read/write access
!----------
! constants
!----------
public :: MDIM                   ! max. rank of NetCDF variables handled
public :: NLEN                   ! max. length of NetCDF variable names
public :: PLEN                   ! max. length of NetCDF file path name
public :: CLOSED, CREATE,       &! file status flag values
          OPEN_R, OPEN_W
!==========
! constants
!==========
!-------------------------------------------------------------------
! size of components in t_netcdf_file (may be increased if required)
!-------------------------------------------------------------------
integer, parameter :: NLEN =  32     ! len  of variables name
integer, parameter :: LLEN =  64     ! len  of longname
integer, parameter :: ULEN =  16     ! len  of units attribute
integer, parameter :: OLEN =   8     ! len  of optionals string
integer, parameter :: ODIM =   8     ! size of optionals array
integer, parameter :: PLEN = 256     ! len  of file path/name
integer, parameter :: MDIM =   4     ! max. number of dimensions
!-----------------------------------------------
! values of t_netcdf_file% status  (file status)
!-----------------------------------------------
integer ,parameter :: UNDEFINED = 0  ! undefined (not yet used ?)
integer ,parameter :: CLOSED    = 1  ! file has been closed
integer ,parameter :: CREATE    = 2  ! file has been created
integer ,parameter :: OPEN_R    = 3  ! file has been opened for reading
integer ,parameter :: OPEN_W    = 4  ! file has been opened for writing


!=========================
! derived type definitions
!=========================
!----------
! container
!----------
type t_netcdf_file
  character(len=PLEN)          :: path     =  ''        ! file path/name
  integer                      :: status   =  UNDEFINED ! file status
  integer                      :: error    =  NF_NOERR  ! error return value
  integer                      :: ndim     =  0         ! number of dimensions
  integer                      :: nvar     =  0         ! number of variables
  integer                      :: ncid     =  0         ! NetCDF file id
  logical                      :: defmode  =  .false.   ! file in define mode?
  type (t_netcdf_dim) ,pointer :: dims (:) => NULL()    ! dimensions
  type (t_netcdf_var) ,pointer :: vars (:) => NULL()    ! variables
end type t_netcdf_file
!-----------
! dimensions
!-----------
type t_netcdf_dim
  character(len=NLEN)    :: name      = ''      ! name of dimension
  integer                :: len       = -1      ! length of dimension
  logical                :: unlimited = .false. ! unlimited dimension ?
  integer                :: dimid     = -1      ! NetCDF dimension id
  integer                :: pos       =  0      ! index of this entry in array
end type t_netcdf_dim
!--------------------------------
! pointer component to dimensions
!--------------------------------
type p_netcdf_dim
  type(t_netcdf_dim) ,pointer :: dim => NULL()
end type p_netcdf_dim
!----------
! variables
!----------
type t_netcdf_var
  character(len=NLEN)    :: name     =  ''         ! name of the variable
  character(len=LLEN)    :: longname =  ''         ! CF convention
  character(len=ULEN)    :: units    =  ''         ! CF convention
  type(t_table) ,pointer :: table    => NULL()     ! table of valid values
  integer                :: invalid  = -huge(1)    ! invalid value for ..
  real(dp)               :: rinvalid = -huge(1._dp)! .. this variable
  integer                :: nvdims   =  0          ! number of dimensions
  type(p_netcdf_dim)     :: p (MDIM)               ! pointer to dimensions
  character(len=OLEN)    :: optionals(ODIM) =  ''  ! mark optional variables
  logical                :: opt_used = .false.     ! indicate used opt.variable
  integer                :: xtype    = -1          ! NetCDF data type
  integer                :: varid    = -1          ! NetCDF variable id
end type t_netcdf_var

contains
!==============================================================================
  subroutine add_dim (file, name, len, pos, unlimited)
  !---------------------------------------------
  ! add a dimension to the NetCDF file meta data
  !---------------------------------------------
  type (t_netcdf_file),intent(inout)        :: file      ! NetCDF meta data
  character(len=*)    ,intent(in)           :: name      ! name of dimension
  integer             ,intent(in)           :: len       ! length of dimension
  integer             ,intent(out)          :: pos       ! index in file%dims
  logical             ,intent(in) ,optional :: unlimited ! unlimited dimension

    logical :: ul

    ul = .false.; if(present(unlimited)) ul = unlimited
    file% ndim = file% ndim + 1
    pos        = file% ndim
    if (pos > size (file% dims)) then
       call finish ("add_dim","table overflow: cannot add "//trim (name))
    end if
    file% dims(pos)% name      =  name
    file% dims(pos)% len       =  len
    file% dims(pos)% pos       =  pos
    file% dims(pos)% unlimited =  ul
  end subroutine add_dim

!------------------------------------------------------------------------------

  subroutine add_var (file, name, xtype, pos, longname, units, table, opt, &
                      invalid, rinvalid)
  !--------------------------------------------
  ! add a variable to the NetCDF file meta data
  !--------------------------------------------
  type (t_netcdf_file),intent(inout)        :: file     ! NetCDF meta data
  character(len=*)    ,intent(in)           :: name     ! name of the variable
  integer             ,intent(in)           :: xtype    ! NetCDF data type
  integer             ,intent(in)           :: pos (:)  ! pointer to dimensions
  character(len=*)    ,intent(in)           :: longname ! CF convention
  character(len=*)    ,intent(in) ,optional :: units    ! CF convention
  type(t_table)       ,pointer    ,optional :: table    ! table to valid values
  character(len=*)    ,intent(in) ,optional :: opt      ! optional flags
  integer             ,intent(in) ,optional :: invalid  ! fillvalue for ints
  integer             ,intent(in) ,optional :: rinvalid ! fillvalue for reals

    integer :: n, i, len
    type (t_netcdf_var) ,pointer :: var

    file% nvar    =  file% nvar + 1
    n             =  file% nvar
    if (n > size (file% vars)) then
       call finish ("add_var","table overflow: cannot add "//trim (name))
    end if
    var           => file% vars(n)
    var% name     =  name
    var% xtype    =  xtype
    var% longname =  longname
    var% nvdims   =  size (pos)
    do i = 1, var% nvdims
      var% p(i)% dim => file% dims(pos(i))
    end do
    select case (xtype)
    case (NF_BYTE)
      var%  invalid = NF_FILL_BYTE
      var% rinvalid = NF_FILL_BYTE
    case (NF_CHAR)
      var%  invalid = NF_FILL_CHAR
      var% rinvalid = NF_FILL_CHAR
    case (NF_SHORT)
      var%  invalid = NF_FILL_SHORT
      var% rinvalid = NF_FILL_SHORT
    case (NF_INT)
      var%  invalid = NF_FILL_INT
      var% rinvalid = NF_FILL_INT
    case (NF_FLOAT)
      var% rinvalid = NF_FILL_FLOAT
    case (NF_DOUBLE)
      var% rinvalid = NF_FILL_DOUBLE
    end select
    if (present(units   )) var% units    =  units
    if (present(table   )) var% table    => table
    if (present(opt     )) call split (var% optionals, opt, len)
    if (present(invalid )) var% invalid  =  invalid
    if (present(rinvalid)) var% rinvalid =  rinvalid
  end subroutine add_var

!------------------------------------------------------------------------------

  subroutine create_netcdf_file (file, cmode, opt)
  type (t_netcdf_file)       ,intent(inout) :: file  ! NetCDF meta data
  integer          ,optional ,intent(in)    :: cmode ! NetCDF creation mode
  character(len=*) ,optional ,intent(in)    :: opt   ! optional parameter flag
  !---------------------------------------------------------------------
  ! Create a NetCDF file from its meta data:
  !   1) create the file
  !   2) define the dimensions
  !   3) define the variables
  !   4) define the variable attributes
  !      (CF conventions: longname, units)
  !
  ! global attributes  must be defined separately
  ! variables contents must be written separately
  !
  ! variables are only defined if the optional parameter flags match the
  ! respective flags set by add_var.
  !---------------------------------------------------------------------

    integer                :: cmod
    integer                :: i, j, len
    integer                :: nvdims, vdims (MDIM)
    character(len=OLEN)    :: optionals (ODIM)

    if (present (opt)) call split (optionals, opt, len)
    !----------------
    ! create the file
    !----------------
    file% status = UNDEFINED
    cmod = NF_CLOBBER; if(present(cmode)) cmod = cmode
    cmod = ior(cmod,NF_64BIT_OFFSET)
    file% error = nf_create (file% path, cmod, file% ncid)
    if (file% error /= NF_NOERR) then
      write(0,*) 'create_netcdf_file: nf_create:',trim(file% path)
      write(0,*) nf_strerror  (file% error)
      return
    endif
    file% defmode = .true.
    !----------------------
    ! define the dimensions
    !----------------------
    do i = 1, file% ndim
      if (file% dims(i)% unlimited) then
        file% dims(i)% len = 0
      else
!       if (file% dims(i)% len == 0) cycle
        if (file% dims(i)% len == 0) file% dims(i)% len = 1
      endif
      file% error = nf_def_dim (file% ncid,          &
                                file% dims(i)% name, &
                                file% dims(i)% len,  &
                                file% dims(i)% dimid )
      if (file% error /= NF_NOERR) then
        write(0,*) 'create_netcdf_file: nf_def_dim ',file% dims(i)% name
        write(0,*) nf_strerror  (file% error)
        return
      endif
    end do
    !---------------------
    ! define the variables
    !---------------------
    do i = 1, file% nvar
      file% vars(i)% opt_used = select (file% vars(i)% optionals, opt)
      if (.NOT. file% vars(i)% opt_used) cycle
      nvdims = file% vars(i)% nvdims
      do j = 1, nvdims
        vdims (j) = file% vars(i)% p(j)% dim% dimid
      end do
      file% error = nf_def_var (file% ncid           ,&
                                file% vars(i)% name  ,&
                                file% vars(i)% xtype ,&
                                file% vars(i)% nvdims,&
                                vdims (1:nvdims)     ,&
                                file% vars(i)% varid  )
      if (file% error /= NF_NOERR) then
        write(0,*) 'ERROR --> create_netcdf_file: nf_def_var ',file% vars(i)% name
        write(0,*) nf_strerror  (file% error)
        return
      endif
    end do
    !--------------------------
    ! write variable attributes
    !--------------------------
    call write_local_attributes

    file% status = CREATE

!------------------------------------------------------------------------------
  contains

    subroutine write_local_attributes

      integer :: i, l

      do i = 1, file% nvar
        if (file% vars(i)% varid /= -1) then
          if (file% vars(i)% units /= '') then
            l = len_trim (file% vars(i)% units)
            file% error = nf_put_att_text (file% ncid,            &
                                           file% vars(i)% varid,  &
                                           'units',               &
                                           l,                     &
                                           file% vars(i)% units   )
            if (file% error /= NF_NOERR) then
              write(0,*) 'ERROR --> create_netcdf_file: nf_put_att_text ',&
                   trim(file% vars(i)% name)//' units='//trim(file% vars(i)% units)
              write(0,*) nf_strerror  (file% error)
              return
            endif
          endif
          if (file% vars(i)% longname /= '') then
            l = len_trim (file% vars(i)% longname)
            file% error = nf_put_att_text (file% ncid,            &
                                           file% vars(i)% varid,  &
                                           'longname',            &
                                           l,                     &
                                           file% vars(i)% longname)
            if (file% error /= NF_NOERR) then
              write(0,*) 'ERROR --> create_netcdf_file: nf_put_att_text ',&
                   trim(file% vars(i)% name)//' longname='//trim(file% vars(i)% longname)
              write(0,*) nf_strerror  (file% error)
              return
            endif
          endif
          select case (file% vars(i)% xtype)
          case default
            file% error = nf_put_att_int    (file% ncid,            &
                                             file% vars(i)% varid,  &
                                             '_FillValue',          &
                                             file% vars(i)% xtype,  &
                                             1,                     &
                                             file% vars(i)% invalid )
            if (file% error /= NF_NOERR) then
              write(0,*) 'ERROR --> create_netcdf_file: nf_put_att_int ',&
                   trim(file% vars(i)% name)//' _FillValue=',file% vars(i)% invalid
              write(0,*) nf_strerror  (file% error)
              return
            endif
          case (NF_FLOAT, NF_DOUBLE)
            file% error = nf_put_att_double (file% ncid,            &
                                             file% vars(i)% varid,  &
                                             '_FillValue',          &
                                             file% vars(i)% xtype,  &
                                             1,                     &
                                             file% vars(i)% rinvalid)
            if (file% error /= NF_NOERR) then
              write(0,*) 'ERROR --> create_netcdf_file: nf_put_att_int ',&
                   trim(file% vars(i)% name)//' _FillValue=',file% vars(i)% rinvalid
              write(0,*) nf_strerror  (file% error)
              return
            endif
          case(NF_CHAR)
          end select
        endif
      end do
    end subroutine write_local_attributes

!------------------------------------------------------------------------------

    logical function select (varopts, opt)
    character(len=OLEN)        :: varopts (ODIM)
    character(len=*) ,optional :: opt
    !----------------------------------------
    ! check if optional variables are written
    !----------------------------------------
      integer :: i
      select = .true.
      if (present (opt)) then
        !-------------------------------------------------
        ! the variable is selected if it is not 'optional'
        ! or no keywords are given
        !-------------------------------------------------
        if (all (varopts == ''   )) return
        if (any (varopts == 'all')) return
        if (all (      optionals == ''   )) return
        if (any (      optionals == 'all')) return
        !------------------------------------------------------
        ! for optional variables a matching keyword is required
        !------------------------------------------------------
        select = .false.
        do i=1, len
          if (    optionals(i) == ''      ) cycle
          if (any(optionals(i) == varopts)) then
            select = .true.
            exit
          end if
        end do
      end if
    end function select

  end subroutine create_netcdf_file

!------------------------------------------------------------------------------

  subroutine open_netcdf_file_read (file)
  type (t_netcdf_file) ,intent(inout) :: file  ! NetCDF meta data
  !-------------------------------------
  ! open the netCDF file for read access
  !-------------------------------------
    file% status = UNDEFINED
    if (file% error /= NF_NOERR) return
    file% error = nf_open (file% path, NF_NOWRITE, file% ncid)
    if (file% error /= NF_NOERR) return
    file% status = OPEN_R
  end subroutine open_netcdf_file_read

!------------------------------------------------------------------------------

  subroutine open_netcdf_file_write (file)
  type (t_netcdf_file) ,intent(inout) :: file  ! NetCDF meta data
  !-------------------------------------------
  ! open the netCDF file for read/write access
  !-------------------------------------------
    file% status = UNDEFINED
    if (file% error /= NF_NOERR) return
    file% error = nf_open (file% path, NF_WRITE, file% ncid)
    if (file% error /= NF_NOERR) return
    file% status = OPEN_W
  end subroutine open_netcdf_file_write

!------------------------------------------------------------------------------

  subroutine close_netcdf_file (file)
  type (t_netcdf_file) ,intent(inout) :: file
  !--------------------
  ! close a NetCDF file
  !--------------------
    file% status = UNDEFINED
    if (file% error /= NF_NOERR) return
    file% error  = nf_close(file% ncid)
    if (file% error /= NF_NOERR) return
    file% status = CLOSED
  end subroutine close_netcdf_file

!------------------------------------------------------------------------------

  subroutine destruct_netcdf_file (file)
  type (t_netcdf_file) ,intent(inout) :: file
  !---------------------------------------------------------
  ! clean up the NetCDF file meta data derived type variable
  ! deallocate components
  !---------------------------------------------------------
    type (t_netcdf_file) :: empty
    if (associated (file% dims)) deallocate (file% dims)
    if (associated (file% vars)) deallocate (file% vars)
    file = empty
  end subroutine destruct_netcdf_file

!------------------------------------------------------------------------------

  pure subroutine split (array, string, len)
  character (len=*) ,intent(out) :: array (:)
  character (len=*) ,intent(in)  :: string
  integer           ,intent(out) :: len
  !------------------------------------------------------------------------
  ! Put each word (separated by blanks) of STRING into an element of ARRAY.
  ! In LEN return the number of words returned or -1 in case of error
  ! (dimension of ARRAY too small).
  !------------------------------------------------------------------------

    integer :: l, i, n, i1, i2, i3

    l = len_trim (string)
    n = size     (array)

    array = ''

    i1  = 0                             ! start position in STRING
    i   = 0                             ! counter for number of words
    len = -1                            ! (error) return argument
    do
      i2 = index(string(i1+1:),' ')     ! position of blank
      if (i2 /= 1)  then                ! skip leading blank
        i3 = l                          ! i3 = end of string
        if (i2 > 1) i3 = i1+i2-1        !      or last char before blank
        if (i3-i1>0) then               ! word bounded by i1..i3
          i = i + 1                     ! next word
          if (i > n) return             ! error return, array too small
          array (i) = string (i1+1:i3)  ! store word
        endif
        if (i2 == 0) exit               ! exit loop, no further blank
      endif
      i1 = i1 + i2                      ! skip word
    end do
    len = i                             ! return number of words

  end subroutine split
!------------------------------------------------------------------------------

end module mo_t_netcdf_file
