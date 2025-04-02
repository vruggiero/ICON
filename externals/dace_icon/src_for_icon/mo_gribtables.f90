!
!+ GRIB table definitions
!
MODULE mo_gribtables
!
! Description:
!   GRIB table definition.
!   includes 'gme_gribtab.incf'
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
!  Add 'center' to gribtable entries
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_19        2012-04-16 Andreas Rhodin
!  new routines dwd_to_v3d, v3d_to_dwd: convert dwd<->3dvar convention
! V1_22        2013-02-13 Andreas Rhodin
!  extend v3d_to_dwd, dwd_to_v3d (work on character strings or arrays)
! V1_27        2013-11-08 Andreas Rhodin
!  add GRIB1 codes for ECHAM
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Author:
! Andreas Rhodin  DWD  2003  original source
!==============================================================================
  !-------------
  ! Modules used
  !-------------
  use mo_kind,        only: wp          ! real kind parameter
  use mo_exception,   only: finish      ! exit on error condition
  use mo_wmo_tables,  only: WMO0_DWD  ,&! center codes
                            WMO0_ECMWF,&! ..
                            WMO0_MPIFM  ! ..
  use mo_dace_string, only: split,     &! string         -> array of words
                            concat      ! array of words -> string
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: ar_des            ! GRIB table entry data type
  public :: dwd_codes         ! DWD GRIB table
  public :: search            ! search in GRIB table
  public :: setup_gribtables  ! set up GRIB tables
  public :: delete_gribtables ! deallocate GRIB tables
  public :: dwd_to_v3d        ! change product names to 3dvar convention
  public :: v3d_to_dwd        ! change product names to DWD  convention

  !------------------------------
  ! 'ar_des' data type definition
  !------------------------------
  include 'gme_pr_gribvar.incf'

  !----------------------
  ! GRIB table definition
  !----------------------
  type (ar_des) ,allocatable :: dwd_codes (:)

  !-------------------------
  ! 'empty' GRIB table entry
  !-------------------------
  type (ar_des) ,parameter :: empty = &
  ar_des('', 0,0,0,0,0,0, 1.0_wp,0.0_wp,'', 0, 0,'','','')

  !-----------
  ! interfaces
  !-----------
  interface v3d_to_dwd            ! change product names to DWD  convention
    module procedure v3d_to_dwd   ! works on character string (list)
    module procedure v3d_to_dwd_1 ! works on character string array
  end interface v3d_to_dwd

  interface dwd_to_v3d            ! change product names to 3dvar convention
    module procedure dwd_to_v3d   ! works on character string (list)
    module procedure dwd_to_v3d_1 ! works on character string array
  end interface dwd_to_v3d

contains
  !-------------------------------------------------------------------
  ! Routine to search in GRIB table
  ! The first entry matching the optional parameter values is returned
  !-------------------------------------------------------------------
  function search (code, table, levtyp, center, name, iname) result (entry)
  type (ar_des)                         :: entry  ! returned table entry
  integer         ,intent(in) ,optional :: code   ! code  number to search
  integer         ,intent(in) ,optional :: table  ! table number to search
  integer         ,intent(in) ,optional :: levtyp ! level type   to search
  integer         ,intent(in) ,optional :: center ! center code
  character(len=*),intent(in) ,optional :: name   ! postprocessing name
  character(len=*),intent(in) ,optional :: iname  ! name used internally
    integer :: i
    if (.not.allocated(dwd_codes)) call setup_gribtables
    do i=1,size (dwd_codes)
      if   (present(code)) then
        if (code   /= dwd_codes(i)% ee)     cycle
      endif
      if   (present(table)) then
        if (table  /= dwd_codes(i)% tabtyp) cycle
      endif
      if   (present(levtyp)) then
        if (levtyp /= dwd_codes(i)% levtyp) cycle
      endif
      if   (present(center)) then
        if (center /= dwd_codes(i)% center) cycle
      endif
      if   (present(name)) then
        if (name   /= dwd_codes(i)% name)   cycle
      endif
      if   (present(iname)) then
        if (iname  /= dwd_codes(i)% iname)  cycle
      endif
      entry = dwd_codes(i)
      return
    end do
    entry = empty
  end function search
!-------------------------------------------------------------------------------
  subroutine setup_gribtables

    include 'gme_gribtab.incf'

  end subroutine setup_gribtables
!-------------------------------------------------------------------------------
  subroutine delete_gribtables

    if(allocated(dwd_codes)) deallocate (dwd_codes)

  end subroutine delete_gribtables
!-------------------------------------------------------------------------------
  subroutine dwd_to_v3d (params)
  character(len=*) ,intent(inout) :: params
  !----------------------------------------------------------
  ! change name of GRIB products from DWD to 3dvar convention
  !----------------------------------------------------------
    integer ,parameter :: nm = 64  ! max. number of entries in 'params'
    character(len=16)  :: pars(nm) ! names of variables to derive
    integer            :: n        ! return code

    call split (pars, params, n)
    if (n<0) call finish('dwd_to_v3d','increase array size "nm"')
    call dwd_to_v3d_1   (pars(1:n))
    call concat (params, pars(1:n))

  end subroutine dwd_to_v3d
!-------------------------------------------------------------------------------
  subroutine dwd_to_v3d_1 (pars)
  character(len=*) ,intent(inout) :: pars (:)
  !----------------------------------------------------------
  ! change name of GRIB products from DWD to 3dvar convention
  !----------------------------------------------------------
    integer            :: m        ! loop index (parameters)
    type (ar_des)      :: entry    ! GRIB table entry

    do m = 1, size (pars)
      entry = search (name=pars(m))
      if (entry% iname /= '') pars(m) = entry% iname
    end do

  end subroutine dwd_to_v3d_1
!-------------------------------------------------------------------------------
  subroutine v3d_to_dwd (params)
  character(len=*) ,intent(inout) :: params
  !----------------------------------------------------------
  ! change name of GRIB products from 3dvar to DWD convention
  !----------------------------------------------------------
    integer ,parameter :: nm = 64  ! max. number of entries in 'params'
    character(len=16)  :: pars(nm) ! names of variables to derive
    integer            :: n        ! return code

    call split (pars, params, n)
    if (n<0) call finish('v3d_to_dwd','increase array size "nm"')
    call v3d_to_dwd_1   (pars(1:n))
    call concat (params, pars(1:n))

  end subroutine v3d_to_dwd
!-------------------------------------------------------------------------------
  subroutine v3d_to_dwd_1 (pars)
  character(len=*) ,intent(inout) :: pars (:)
  !----------------------------------------------------------
  ! change name of GRIB products from 3dvar to DWD convention
  !----------------------------------------------------------
    integer            :: m        ! loop index (parameters)
    type (ar_des)      :: entry    ! GRIB table entry

    do m = 1, size (pars)
      entry = search (iname=pars(m))
      if (entry% iname /= '') pars(m) = entry% name
    end do

  end subroutine v3d_to_dwd_1
!------------------------------------------------------------------------------
end module mo_gribtables
