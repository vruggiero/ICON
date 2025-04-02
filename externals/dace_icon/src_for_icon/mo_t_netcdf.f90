!
!+ Subroutines and global variables to handle NetCDF I/O (uses f90 interface)
!
MODULE mo_t_netcdf
!
! Description:
!   Subroutines and global variables to handle NetCDF I/O.
!   This module uses the Fortran90 interface to NetCDF.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  Cope with fillvalues
! V1_8         2009/12/09 Andreas Rhodin
!  new specific subroutine get_attr_real
!  new optional arguments: subroutine get_dim (ierr), get_attr_char (lena)
! V1_9         2010/04/20 Andreas Rhodin
!  get_var_int_1: correctly handle fillvalues
!  new specific subroutines: get_var_real_3, get_var_int_3
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  revised handling for fillvalues
! V1_20        2012-06-18 Harald Anlauf
!  mo_t_netcdf: USE netcdf, only; fix comments
! V1_23        2013-03-26 Robin Faulwetter
!  fix error message
! V1_31        2014-08-21 Andreas Rhodin
!  additional optional parametrs to get_var (real(:,:),...)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2002-2008
!==============================================================================

!#ifdef __SX__
!-----------------------------------------------------
! integer(i1) is not supported (in this case i1 == i2)
!-----------------------------------------------------
!#undef HAVE_I1
!#else
!-------------------------------------------------
! integer(i1) is supported (in this case i1 /= i2)
!-------------------------------------------------
#define HAVE_I1
!#endif

!=============
! Modules used
!=============
use mo_kind,      only: wp, sp, dp, i1,i2 ! kind parameters
use mo_exception, only: finish            ! abort routine
use mo_ascii,     only: NUL               ! achar ( 00)
!---------------------
! netCDF f90 interface
!---------------------
use netcdf,       only: nf90_Inquire_Attribute,&!
                        nf90_Inquire_Dimension,&!
                        nf90_inq_dimid,        &!
                        nf90_inq_varid,        &!
                        nf90_def_var,          &!
                        nf90_get_var,          &!
                        nf90_get_att,          &!
                        nf90_put_att,          &!
                        nf90_strerror,         &!
                        NF90_CHAR,             &!
                        NF90_INT,              &!
                        NF90_FLOAT,            &!
                        NF90_FILL_BYTE,        &!
                        NF90_FILL_INT,         &!
                        NF90_FILL_INT2,        &!
                        NF90_FILL_FLOAT,       &!
                        NF90_GLOBAL,           &!
                        NF90_MAX_NAME,         &!
                        NF90_NOERR              !
implicit none

!================
! Public entities
!================
private
!-----------------
! Module variables
!-----------------
public :: ncid        ! NetCDF file id
public :: vname       ! variable name or other comment
public :: xtype       ! type  of variable to be defined
public :: dimid       ! dimid of variable to be defined
public :: stanc       ! NetCDF start  parameter
public :: counc       ! NetCDF count  parameter
public :: strnc       ! NetCDF stride parameter
public :: rname       ! routine name + netcdf routine called
public :: dimid2      ! dimension of variable(2) to be defined
public :: dimid3      ! dimension of variable(3) to be defined
!------------------
! Module procedures
!------------------
public :: chk          ! checks error status after netcdf call
public :: def_var      ! calls nf90_def_var
public :: def_var_2    ! as def_var but with dimension dimid2(2)
public :: def_var_3    ! as def_var but with dimension dimid3(3)
public :: def_str      ! define a character string variable
public :: get_var      ! get variable
public :: get_attr     ! get attribute
public :: get_dim      ! get dimension from Netcdf file
!=================
! Module variables
!=================
  integer            :: ncid        ! NetCDF file id
  character(len=64)  :: vname =''   ! variable name ot other comment
  integer            :: xtype       ! type  of variable to be defined
  integer            :: dimid       ! dimid of variable to be defined
  integer            :: stanc(3)    ! NetCDF start  parameter
  integer            :: counc(3)    ! NetCDF count  parameter
  integer            :: strnc(3)    ! NetCDF stride parameter
  character(len=32)  :: rname       ! routine name + netcdf routine called
  integer            :: dimid2(2)   ! dimension of variable(2) to be defined
  integer            :: dimid3(3)   ! dimension of variable(3) to be defined

!===========
! Interfaces
!===========
interface get_var
  module procedure get_var_int_0   ! get scalar integer variable
  module procedure get_var_int_1   ! get 1-d    integer variable
  module procedure get_var_int_2   ! get 2-d    integer variable
#ifdef HAVE_I1
  module procedure get_var_byte_1  ! get 1-d    integer variable (byte size)
#endif
  module procedure get_var_int2_1  ! get 1-d    integer variable (2 byte)
  module procedure get_var_int_3   ! get 3-d    integer variable
  module procedure get_var_real_0  ! get scalar real    variable
  module procedure get_var_real4_1 ! get 1-d    real    variable
  module procedure get_var_real4_2 ! get 2-d    real    variable
  module procedure get_var_real8_1 ! get 1-d    real    variable
  module procedure get_var_real8_2 ! get 2-d    real    variable
  module procedure get_var_real8_3 ! get 3-d    real    variable
! module procedure get_var_real_3  ! get 3-d    real    variable
  module procedure get_var_char_0  ! get character      variable
  module procedure get_var_char_1  ! get character      variable
end interface get_var

interface get_attr
  !--------------
  ! get attribute
  !--------------
  module procedure get_attr_char   ! for character string
  module procedure get_attr_intp   ! for integer pointer array
  module procedure get_attr_int    ! for integer array
  module procedure get_attr_real   ! for real    scalar
end interface get_attr

contains
!==============================================================================
! Usefull netcdf subroutines
!==============================================================================
  subroutine chk (status, ierr)
    integer, intent(in)            :: status
    integer, intent(out) ,optional :: ierr
    !--------------------------------------------------------
    ! checks error status after each netcdf, prints out text
    ! message each time an error code is returned and aborts.
    !--------------------------------------------------------
    if (present(ierr)) then
      ierr = status
    else
      if(status /= NF90_NOERR) then
        call finish(trim(rname),trim(nf90_strerror(status))//' : '//trim(vname))
      end if
    endif
    vname = trim(vname)//'?'
  end subroutine chk
!----------------------------------------------------------------------------
  subroutine def_var (name, varid, units, long_name)
    character(len=*) ,intent(in)  :: name, units, long_name
    integer          ,intent(out) :: varid
    !--------------------------------------------------------------
    ! calls nf90_def_var with actual parameters: name, varid
    !               and module variables: ncid, xtype, dimid
    !--------------------------------------------------------------
    vname = name
    call chk(nf90_def_var (ncid, name, xtype, dimid, varid))
    if(units    /='')call chk(nf90_put_att(ncid,varid,'units',    units))
    if(long_name/='')call chk(nf90_put_att(ncid,varid,'long_name',long_name))
    end subroutine def_var
!----------------------------------------------------------------------------
  subroutine def_var_2 (name, varid, units, long_name)
    character(len=*) ,intent(in)  :: name, units,long_name
    integer          ,intent(out) :: varid
    !-----------------------------------------------------------------
    ! as def_var but with dimension dimid2(2) and missing values
    !-----------------------------------------------------------------
    vname = name
    call chk(nf90_def_var (ncid, name, xtype, dimid2, varid))
    if(units    /='')call chk(nf90_put_att(ncid,varid,'units',    units))
    if(long_name/='')call chk(nf90_put_att(ncid,varid,'long_name',long_name))
    if(xtype == NF90_INT) &
       call chk(nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_INT))
    if(xtype == NF90_FLOAT) &
       call chk(nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_FLOAT))
  end subroutine def_var_2
!----------------------------------------------------------------------------
  subroutine def_var_3 (name, varid, units, long_name)
    character(len=*) ,intent(in)  :: name, units,long_name
    integer          ,intent(out) :: varid
    !--------------------------------------------------------------
    ! as def_var but with dimension dimid3(3) and missing values
    !--------------------------------------------------------------
    vname = name
    call chk(nf90_def_var (ncid, name, xtype, dimid3, varid))
    if(units    /='')call chk(nf90_put_att(ncid,varid,'units',    units))
    if(long_name/='')call chk(nf90_put_att(ncid,varid,'long_name',long_name))
    if(xtype == NF90_INT) &
       call chk(nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_INT))
    if(xtype == NF90_FLOAT) &
       call chk(nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_FLOAT))
  end subroutine def_var_3
!----------------------------------------------------------------------------
  subroutine def_str (name, varid, dimstr,long_name)
    character(len=*) ,intent(in)  :: name, long_name
    integer          ,intent(in)  :: dimstr
    integer          ,intent(out) :: varid
    !--------------------------------------------------------------
    ! define a character string variable.
    ! calls nf90_def_var with actual parameters: name, varid, dimstr
    !               and module variables: ncid, xtype, dimid
    !--------------------------------------------------------------
    vname = name
    call chk(nf90_def_var (ncid, name, NF90_CHAR, (/dimstr,dimid/), varid))
    if(long_name/='')call chk(nf90_put_att(ncid,varid,'long_name',long_name))
  end subroutine def_str
!----------------------------------------------------------------------------
  subroutine get_dim (dim ,name, ierr)
  !-------------------------------
  ! get dimension from Netcdf file
  !-------------------------------
  integer           ,intent(out)           :: dim
  character(len=*)  ,intent(in)            :: name
  integer           ,intent(out) ,optional :: ierr
    integer                        :: vdimid
    character(len = NF90_MAX_NAME) :: voname
    rname = 'get_dim: nf90_inq_dimid'
    vname = name
    call chk (nf90_inq_dimid  (ncid , name , vdimid) ,ierr)    ! name input
    if(present(ierr)) then
      if (ierr/= NF90_NOERR) return
    endif
    rname = 'get_dim: nf90_inquire_dimension'
    call chk (nf90_inquire_dimension (ncid, vdimid, voname, dim))! vname output
  end subroutine get_dim
!----------------------------------------------------------------------------
  subroutine get_var_int_0 (int ,name)
    !---------------------
    ! get 0-D int variable
    !---------------------
    integer          ,intent(out) :: int
    character(len=*) ,intent(in)  :: name
    integer :: varid
    vname = name
    rname = 'get_var_int_0'
    call chk (nf90_inq_varid (ncid, name, varid))
    call chk (nf90_get_var (ncid, varid, int , stanc ))
  end subroutine get_var_int_0
!----------------------------------------------------------------------------
  subroutine get_var_int_1 (int ,name, ierr, fill)
  !---------------------
  ! get 1-D int variable
  !---------------------
  integer          ,intent(out)           :: int (:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
  integer          ,intent(in)  ,optional :: fill
    integer :: varid
    integer :: lc(3)
    integer :: fillvalue
    vname = name
    rname = 'get_var_int_1'
    varid = -1
    call chk (nf90_inq_varid (ncid, name, varid),ierr)
    if(present(ierr)) then
      if (ierr/= NF90_NOERR) then
        if (present(fill)) then
          int = fill
        else
          int = NF90_FILL_INT
        endif
        return
      endif
    endif
    lc = 1; lc(1) = size(int); if(counc(1)/=0) lc = counc
    call chk (nf90_get_var (ncid, varid, int, stanc, lc, strnc))
    if (present(fill)) then
      call chk (nf90_get_att (ncid, varid, '_FillValue', fillvalue))
      where (int == fillvalue) int = fill
    endif
  end subroutine get_var_int_1
!----------------------------------------------------------------------------
#ifdef HAVE_I1
  subroutine get_var_byte_1 (int ,name, ierr, fill)
  !---------------------
  ! get 1-D int variable
  !---------------------
  integer(i1)      ,intent(out)           :: int (:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
  integer(i1)      ,intent(in)  ,optional :: fill
    integer :: varid
    integer :: lc(3)
    vname = name
    rname = 'get_var_byte_1'
    varid = -1
    lc = 1; lc(1) = size(int); if(counc(1)/=0) lc = counc
    call chk (nf90_inq_varid (ncid, name, varid),ierr)
    if(present(ierr)) then
      if (ierr/= NF90_NOERR) then
        if (present(fill)) then
          int = fill
        else
          int = NF90_FILL_BYTE
        endif
        return
      endif
    endif
    call chk (nf90_get_var (ncid, varid, int, stanc, lc, strnc ))
  end subroutine get_var_byte_1
#endif
!----------------------------------------------------------------------------
  subroutine get_var_int2_1 (int ,name, ierr, fill)
  !---------------------
  ! get 1-D int variable
  !---------------------
  integer(i2)      ,intent(out)           :: int (:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
  integer(i2)      ,intent(in)  ,optional :: fill
    integer :: varid
    integer :: lc(3)
    vname = name
    rname = 'get_var_int2_1'
    varid = -1
    lc = 1; lc(1) = size(int); if(counc(1)/=0) lc = counc
    call chk (nf90_inq_varid (ncid, name, varid),ierr)
    if(present(ierr)) then
      if (ierr/= NF90_NOERR) then
        if (present(fill)) then
          int = fill
        else
          int = NF90_FILL_INT2
        endif
        return
      endif
    endif
    call chk (nf90_get_var (ncid, varid, int, stanc, lc, strnc ))
  end subroutine get_var_int2_1
!----------------------------------------------------------------------------
  subroutine get_var_real_0 (x ,name)
    !---------------------
    ! get 0-D real variable
    !---------------------
    real(wp)         ,intent(out) :: x
    character(len=*) ,intent(in)  :: name
    integer :: varid
    vname = name
    rname = 'get_var_real_0'
    call chk (nf90_inq_varid (ncid, name, varid))
    call chk (nf90_get_var (ncid, varid,  x, stanc ))
  end subroutine get_var_real_0
!----------------------------------------------------------------------------
  subroutine get_var_real8_1 (x ,name, ierr, fill)
  !----------------------
  ! get 1-D real variable
  !----------------------
  real(dp)         ,intent(out)           :: x (:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
  real(dp)         ,intent(in)  ,optional :: fill
    integer  :: varid
    integer  :: lc(3)
    integer  :: ie
    real(dp) :: fillvalue
    vname = name
    rname = 'get_var_real8_1'
    varid = -1
    if (present(ierr)) ierr = 0
    call chk (nf90_inq_varid (ncid, name, varid),ie)
    if (ie /= NF90_NOERR) then
      if (present(ierr)) ierr = ie
      if (present(fill)) then
                         x    = fill
      else
                         x    = NF90_FILL_FLOAT
      end if
      if (present(ierr).or.present(fill)) return
      call finish('get_var_real8_1','cannot read '//trim(name))
    endif
    lc = 1; lc(1) = size(x); if(counc(1)/=0) lc = counc
    x = -huge (x)          !+++ workaround for NAGf95 compiler bug (-nan)
    call chk (nf90_get_var (ncid, varid, x, stanc, lc, strnc))
    if (present(fill)) then
      call chk (nf90_get_att (ncid, varid, '_FillValue', fillvalue))
      where (x == fillvalue) x = fill
    endif
  end subroutine get_var_real8_1
!----------------------------------------------------------------------------
  subroutine get_var_real4_1 (x ,name, ierr, fill)
  !----------------------
  ! get 1-D real variable
  !----------------------
  real(sp)         ,intent(out)           :: x (:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
  real(sp)         ,intent(in)  ,optional :: fill
    integer  :: varid
    integer  :: lc(3)
    integer  :: ie
    real(sp) :: fillvalue
    vname = name
    rname = 'get_var_real4_1'
    varid = -1
    if (present(ierr)) ierr = 0
    call chk (nf90_inq_varid (ncid, name, varid),ie)
    if (ie /= NF90_NOERR) then
      if (present(ierr)) ierr = ie
      if (present(fill)) x    = fill
      if (present(ierr).or.present(fill)) return
      call finish('get_var_real4_1','cannot read '//trim(name))
    endif
    lc = 1; lc(1) = size(x); if(counc(1)/=0) lc = counc
    call chk (nf90_get_var (ncid, varid, x, stanc, lc, strnc))
    if (present(fill)) then
      call chk (nf90_get_att (ncid, varid, '_FillValue', fillvalue))
      where (x == fillvalue) x = fill
    endif

  end subroutine get_var_real4_1

!----------------------------------------------------------------------------
  subroutine get_var_real8_2 (x ,name, ierr, fill)
  !----------------------
  ! get 2-D real variable
  !----------------------
  real(dp)         ,intent(out)           :: x (:,:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
  real(dp)         ,intent(in)  ,optional :: fill
    integer  :: varid
    integer  :: lc(3)
    integer  :: ie
    real(dp) :: fillvalue
    vname = name
    rname = 'get_var_real8_2'
    varid = -1
    if (present(ierr)) ierr = 0
    call chk (nf90_inq_varid (ncid, name, varid),ie)
    if (ie /= NF90_NOERR) then
      if (present(ierr)) ierr = ie
      if (present(fill)) then
                         x    = fill
      else
                         x    = NF90_FILL_FLOAT
      end if
      if (present(ierr).or.present(fill)) return
      call finish('get_var_real8_2','cannot read '//trim(name))
    endif
    lc = 1; lc(1:2) = shape(x); if(counc(1)/=0) lc = counc
    x = -huge (x)          !+++ workaround for NAGf95 compiler bug (-nan)
    call chk (nf90_get_var (ncid, varid, x, stanc, lc, strnc))
    if (present(fill)) then
      call chk (nf90_get_att (ncid, varid, '_FillValue', fillvalue))
      where (x == fillvalue) x = fill
    endif
  end subroutine get_var_real8_2
!----------------------------------------------------------------------------
  subroutine get_var_real4_2 (x ,name, ierr, fill)
  !----------------------
  ! get 2-D real variable
  !----------------------
  real(sp)         ,intent(out)           :: x (:,:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
  real(sp)         ,intent(in)  ,optional :: fill
    integer  :: varid
    integer  :: lc(3)
    integer  :: ie
    real(sp) :: fillvalue
    vname = name
    rname = 'get_var_real4_2'
    varid = -1
    if (present(ierr)) ierr = 0
    call chk (nf90_inq_varid (ncid, name, varid),ie)
    if (ie /= NF90_NOERR) then
      if (present(ierr)) ierr = ie
      if (present(fill)) x    = fill
      if (present(ierr).or.present(fill)) return
      call finish('get_var_real4_2','cannot read '//trim(name))
    endif
    lc = 1; lc(1:2) = shape(x); if(counc(1)/=0) lc = counc
    call chk (nf90_get_var (ncid, varid, x, stanc, lc, strnc))
    if (present(fill)) then
      call chk (nf90_get_att (ncid, varid, '_FillValue', fillvalue))
      where (x == fillvalue) x = fill
    endif

  end subroutine get_var_real4_2
!----------------------------------------------------------------------------
  subroutine get_var_real8_3 (x ,name, ierr, fill)
  !----------------------
  ! get 3-D real variable
  !----------------------
  real(dp)         ,intent(out)           :: x (:,:,:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
  real(dp)         ,intent(in)  ,optional :: fill
    integer  :: varid
    integer  :: lc(3)
    integer  :: ie
    real(dp) :: fillvalue
    vname = name
    rname = 'get_var_real8_3'
    varid = -1
    if (present(ierr)) ierr = 0
    call chk (nf90_inq_varid (ncid, name, varid),ie)
    if (ie /= NF90_NOERR) then
      if (present(ierr)) ierr = ie
      if (present(fill)) then
                         x    = fill
      else
                         x    = NF90_FILL_FLOAT
      end if
      if (present(ierr).or.present(fill)) return
      call finish('get_var_real8_3','cannot read '//trim(name))
    endif
    lc = shape(x); if(counc(1)/=0) lc = counc
    call chk (nf90_get_var (ncid, varid, x, stanc, lc, strnc))
    if (present(fill)) then
      call chk (nf90_get_att (ncid, varid, '_FillValue', fillvalue))
      where (x == fillvalue) x = fill
    endif
  end subroutine get_var_real8_3
!----------------------------------------------------------------------------
  subroutine get_var_int_2 (x ,name, ierr, fill)
  !-------------------------
  ! get 2-D integer variable
  !-------------------------
  integer          ,intent(out)           :: x (:,:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
  integer          ,intent(in)  ,optional :: fill
    integer :: varid
    integer :: lc(3)
    integer :: fillvalue
    vname = name
    rname = 'get_var_int_2'
    varid = -1
    call chk (nf90_inq_varid (ncid, name, varid), ierr)
    if(present(ierr)) then
      if (ierr/= NF90_NOERR) then
        if (present(fill)) then
          x = fill
        else
          x = NF90_FILL_INT
        endif
        return
      endif
    endif
    lc = 1; lc(1:2) = shape(x); if(counc(1)/=0) lc = counc
    call chk (nf90_get_var (ncid, varid, x, stanc, lc, strnc))
    if (present(fill)) then
      call chk (nf90_get_att (ncid, varid, '_FillValue', fillvalue))
      where (x == fillvalue) x = fill
    endif
  end subroutine get_var_int_2
!----------------------------------------------------------------------------
! subroutine get_var_real_3 (x, name, ierr)
! !----------------------
! ! get 3-D real variable
! !----------------------
! real(wp)         ,intent(out)           :: x (:,:,:)
! character(len=*) ,intent(in)            :: name
! integer          ,intent(out) ,optional :: ierr
!   integer :: varid
!   integer :: lc(3)
!   vname = name
!   rname = 'get_var_real_3'
!   varid = -1
!   call chk (nf90_inq_varid (ncid, name, varid), ierr)
!   if (present(ierr)) then
!     if (ierr /= NF90_NOERR) then
!       x = NF90_FILL_FLOAT
!       return
!     endif
!   endif
!   lc = shape(x); if(counc(1)/=0) lc = counc
!   call chk (nf90_get_var (ncid, varid, x, stanc, lc, strnc))
! end subroutine get_var_real_3
!----------------------------------------------------------------------------
  subroutine get_var_int_3 (x, name, ierr)
  !-------------------------
  ! get 3-D integer variable
  !-------------------------
  integer          ,intent(out)           :: x (:,:,:)
  character(len=*) ,intent(in)            :: name
  integer          ,intent(out) ,optional :: ierr
    integer :: varid
    integer :: lc(3)
    vname = name
    rname = 'get_var_int_3'
    varid = -1
    call chk (nf90_inq_varid (ncid, name, varid), ierr)
    if (present(ierr)) then
      if (ierr /= NF90_NOERR) then
        x = NF90_FILL_INT
        return
      endif
    endif
    lc = shape(x); if(counc(1)/=0) lc = counc
    call chk (nf90_get_var (ncid, varid, x, stanc, lc, strnc))
  end subroutine get_var_int_3
!----------------------------------------------------------------------------
  subroutine get_var_char_0 (x ,name)
  !---------------------------
  ! get character variable
  !---------------------------
  character(len=*) ,intent(out) :: x
  character(len=*) ,intent(in)  :: name
    integer :: varid
    integer :: lc(3)
    integer :: i
    vname = name
    rname = 'get_var_char_0'
    lc=1; lc(1)=len(x)
    call chk (nf90_inq_varid (ncid, name, varid))
    call chk (nf90_get_var(ncid,varid,x,(/1,stanc/), lc, (/1,strnc/)))
    !--------------------------------------
    ! convert spurious NUL-character to ' '
    !--------------------------------------
    do i=1,len(x)
      if (x(i:i) == NUL) x(i:i) = ' '
    end do
  end subroutine get_var_char_0
!----------------------------------------------------------------------------
  subroutine get_var_char_1 (x ,name)
  !---------------------------
  ! get character variable
  !---------------------------
  character(len=*) ,intent(out) :: x (:)
  character(len=*) ,intent(in)  :: name
    integer               :: varid
    integer               :: lc(2)
    integer               :: i, j
    character(len=len(x)) :: z (size(x))
    vname = name
    rname = 'get_var_char_1'
    z     = ''

    lc(1)=len(x)
    lc(2)= size(x); if(counc(1)/=0) lc(2) = counc(1)
    if (lc(2) > size(x)) call finish('get_var_char_1','counc(1) > size(x)')

    call chk (nf90_inq_varid (ncid, name, varid))
    call chk (nf90_get_var(ncid,varid,z,(/1,stanc(1)/), lc, (/1,strnc(1)/)))

    !--------------------------------------
    ! convert spurious NUL-character to ' '
    !--------------------------------------
    do j = 1, size (z)
      do i=1,len(z)
        if (z(j)(i:i) == NUL) z(j)(i:i) = ' '
      end do
    end do

    x = z

  end subroutine get_var_char_1
!------------------------------------------------------------------------------
  subroutine get_attr_char (attr, varname, name, lena)
  character(len=*)  ,intent(out) :: attr    ! attribute
  character(len=*)  ,intent(in)  :: varname ! name of variable, '' for global
  character(len=*)  ,intent(in)  :: name    ! name of attribute
  integer ,optional ,intent(out) :: lena    ! len of attribute in NetCDF file
  !-------------------------------------------------------------------------
  ! get attribute (character string)
  ! checks for sufficient len of actual parameter 'attr'
  ! aborts if len is insufficient and actual parameter 'lena' is not present
  !-------------------------------------------------------------------------
    integer             :: varid
    integer             :: ln    ! len in NetCDF file
    integer             :: la    ! len of actual argument
    character(len=1024) :: temp  ! temporary
    !----------
    ! get varid
    !----------
    varid = NF90_GLOBAL
    if (varname /= '') then
      vname = varname
      rname = 'nf90_inq_varid'
      call chk (nf90_inq_varid (ncid, varname, varid))
    endif
    !----------------------------------------
    ! get attribute, check for sufficient len
    !----------------------------------------
    call chk (nf90_inquire_attribute(ncid, varid, name, len=ln))
    vname = trim(varname)//':'//name
    rname = 'nf90_get_att'
    attr  = ''
    la = len(attr)
    if (present(lena)) then
      lena = ln
      if (ln <= la) then
        call chk (nf90_get_att(ncid, varid, name, attr))
      elseif (ln > len(temp)) then
        call finish ('get_attr_char','ln > len(temp) : '//trim(vname))
      else
        call chk (nf90_get_att(ncid, varid, name, temp))
        attr = temp(1:la)
      endif
    else
      if (ln > la) call finish ('get_attr_char','ln > la : '//trim(vname))
      call chk (nf90_get_att(ncid, varid, name, attr))
    endif
    attr (ln+1:la) = ''
  end subroutine get_attr_char
!------------------------------------------------------------------------------
  subroutine get_attr_intp (attr, varname, name)
  integer          ,pointer     :: attr(:) ! attribute
  character(len=*) ,intent(in)  :: varname ! name of variable, '' for global
  character(len=*) ,intent(in)  :: name    ! name of attribute

    integer :: varid, attlen
    varid = NF90_GLOBAL
    if (varname /= '') then
      vname = varname
      rname = 'nf90_inq_varid'
      call chk (nf90_inq_varid (ncid, varname, varid))
    endif
    vname = trim(varname)//':'//name
    rname = 'nf90_Inquire_Attribute'
    call chk (nf90_Inquire_Attribute(ncid, varid, name, len=attlen))
    if (associated(attr)) deallocate (attr)
    allocate (attr(attlen))
    rname = 'nf90_get_att'
    call chk (nf90_get_att(ncid, varid, name, attr))
  end subroutine get_attr_intp
!------------------------------------------------------------------------------
  subroutine get_attr_int (attr, varname, name)
  integer          ,intent(out) :: attr    ! attribute
  character(len=*) ,intent(in)  :: varname ! name of variable, '' for global
  character(len=*) ,intent(in)  :: name    ! name of attribute

    integer :: varid, attlen
    varid = NF90_GLOBAL
    if (varname /= '') then
      vname = varname
      rname = 'nf90_inq_varid'
      call chk (nf90_inq_varid (ncid, varname, varid))
    endif
    vname = trim(varname)//':'//name
    rname = 'nf90_Inquire_Attribute'
    call chk (nf90_Inquire_Attribute(ncid, varid, name, len=attlen))
    if (attlen /= 1) &
      call finish('get_attr_int','attlen /= 1 :'//trim(varname))
    rname = 'nf90_get_att'
    call chk (nf90_get_att(ncid, varid, name, attr))
  end subroutine get_attr_int
!------------------------------------------------------------------------------
  subroutine get_attr_real (attr, varname, name)
  real(wp)         ,intent(out) :: attr    ! attribute
  character(len=*) ,intent(in)  :: varname ! name of variable, '' for global
  character(len=*) ,intent(in)  :: name    ! name of attribute

    integer :: varid, attlen
    varid = NF90_GLOBAL
    if (varname /= '') then
      vname = varname
      rname = 'nf90_inq_varid'
      call chk (nf90_inq_varid (ncid, varname, varid))
    endif
    vname = trim(varname)//':'//name
    rname = 'nf90_Inquire_Attribute'
    call chk (nf90_Inquire_Attribute(ncid, varid, name, len=attlen))
    if (attlen /= 1) &
      call finish('get_attr_real','attlen /= 1 :'//trim(varname))
    rname = 'nf90_get_att'
    call chk (nf90_get_att(ncid, varid, name, attr))
  end subroutine get_attr_real
!------------------------------------------------------------------------------

end module mo_t_netcdf
