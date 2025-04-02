!
!+ derived type to hold grid-point, Fourier, spectral, icosahedral fields
!
MODULE mo_memory
!
! Description:
!   Definition of data type 't_m' to hold 4D fields in
!   grid-point, Fourier, spectral, or icosahedral representation.
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
!  Changed printout
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_19        2012-04-16 Andreas Rhodin
!  loop collapsing (minval,maxval) for print(atmospheric state)
! V1_23        2013-03-26 Harald Anlauf
!  Handle GRIB2 triples for fields, setting of grib1/grib2 codes from GRIB_API
! V1_26        2013/06/27 Andreas Rhodin
!  print_atm_state: print mean in addition to minimum/maximum
! V1_29        2014/04/02 Andreas Rhodin
!  new specific function: integer * atmospheric_state
! V1_31        2014-08-21 Harald Anlauf
!  print_m: indicate level(s) with min. or max. value if verbose=true
! V1_42        2015-06-08 Andreas Rhodin
!  implement max(real, t_atm)
! V1_43        2015-08-19 Harald Anlauf
!  Improve GRIB encoding of single-level/surface fields for flake
! V1_44        2015-09-30 Andreas Rhodin
!  print_m: account for special values huge and -huge
! V1_45        2015-12-15 Andreas Rhodin
!  print(t_atm): saveguard for values of +-huge
! V1_47        2016-06-06 Harald Anlauf
!  Enable GRIB encoding of selected analysis fields with 24 bits.
!  Allow up to 5 time range entries.
! V1_51        2017-02-24 Andreas Rhodin
!  new subroutines p_min_atm,p_max_atm
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD  2001-2007
!============================================================
#include "tr15581.incf"
  !-------------
  ! modules used
  !-------------
  use mo_kind,      only: wp, i8           ! working precision kind parameter
  use mo_exception, only: finish           ! abort routine
  use mo_transform, only: lgtd, lgti,     &! direct and inverse legendre trans.
                          lgtd_ad,lgti_ad,&! adjoint legendre transforms
                          nsp_nn           ! get number of spectral coeffs.
! use mo_grads,     only: t_ctl,          &! data type to hold .ctl file
!                         write_var        ! write field
  use mo_mpi_dace,  only: dace,           &! dace communication info
                          p_sum,          &! sum over processor elements
                          p_min,          &! min over processor elements
                          p_max            ! max over processor elements
  implicit none
!------------------------------------------------------------------------------
  !----------------
  ! public entities
  !----------------
  private
  public :: t_m              ! table to hold 3d fields
  public :: t_mi             ! metadata component of t_m

  public :: construct        ! construct a table
  public :: destruct         ! destruct a table
  public :: new_entry        ! establish a new entry
  public :: del_entry        ! delete entry
  public :: index            ! return the index of 3d field from its name
  public :: pt4,pt3,pt2,pt1  ! return a pointer to 3d field from its name

  public :: update           ! update field bounds and allocation status
  public :: allocate         ! allocate fields in table
  public :: deallocate       ! deallocate fields in table

  public :: assignment (=)   ! assign real number to table
  public :: assign_m_reals   ! PRELIMINARY
  public :: plus             ! addition
  public :: minus            ! subtraction
  public :: times            ! multiplication
  public :: over             ! division
  public :: power            ! exponentiation
  public :: sqrt_m           ! square root
  public :: max_m            ! max value
  public :: sum              ! sum of field elements
  public :: operator (==)    ! equal t_mis
  public :: operator (/=)    ! not equal t_mis
  public :: lgtd, lgti       ! direct and inverse legendre transform
  public :: lgtd_ad, lgti_ad ! adjoint routines
  public :: print            ! print information on table entries
  public :: dump             ! dump table to binary file
  public :: retrieve         ! retrieve table from binary file
  public :: nwords           ! number of real array elements allocated
  public :: cut              ! cutout region
  public :: p_min_m          ! sum over processor elements
  public :: p_max_m          ! sum over processor elements
  public :: p_sum_m          ! sum over processor elements
!==============================================================================
  integer, parameter :: mt = 5 ! max. number of time range entries
  !---------------
  ! meta data type
  !---------------
  type t_mi
    integer           :: lb(4)  = (/1,1,1,1/) ! lower bound of grid
    integer           :: ub(4)  = (/0,1,1,1/) ! upper bound of grid
    integer           :: nn     = 0           ! triangular truncation
    integer           :: nsp    = 0           ! number of spectral coefficients
    character(len=16) :: name   = '        '  ! name
    integer           :: code   = 255         ! grib1 code  number
    integer           :: table  = 255         ! grib1 table number
    integer           :: parID  = -1          ! grib_api parameter ID
    integer           :: dis    = 255         ! grib2 discipline
    integer           :: cat    = 255         ! grib2 parameter category
    integer           :: num    = 255         ! grib2 parameter number
    integer           :: ctyp   = -1          ! grib2 constituentType
    integer           :: ffs    = 255         ! grib2 typeOfFirstFixedSurface
    integer           :: sfs    = 255         ! grib2 typeOfSecondFixedSurface
    logical           :: alloc  = .false.     ! true if allocated
    character(len=2)  :: rep    = 'gg'        ! representation: 'gg','ff','sh'
    logical           :: ref    = .false.     ! flag for reference value
    integer           :: tr(mt) =  0          ! time ranges (min)
    integer           :: range  = -1          ! time range indicator
    integer           :: bits   = -1          ! bits used for GRIB encoding
    logical           :: dummy  = .false.     ! dummy field, not eligible for I/O
  end type t_mi
!------------------------------------------------------------------------------
  !----------------------------
  ! data type to hold 3D fields
  !----------------------------
  type t_m
    type (t_mi)           :: i                         ! meta data
    real(wp) ,ALLOCATABLE :: ptr (:,:,:,:) DEFAULTNULL ! field
  end type t_m
!==============================================================================
  !-----------
  ! interfaces
  !-----------
  interface construct
    module procedure construct_m
  end interface construct
!------------------------------------------------------------------------------
  interface destruct
    module procedure destruct_m
#ifndef M3ELEMENTAL
    module procedure destruct_ms
#endif
  end interface destruct
!------------------------------------------------------------------------------
  interface update
    module procedure update_m
#ifndef M3ELEMENTAL
    module procedure update_ms
#endif
  end interface update
!------------------------------------------------------------------------------
  interface index
    module procedure index_m
  end interface index
!------------------------------------------------------------------------------
  interface pt4
    module procedure pt4_m
  end interface pt4
!------------------------------------------------------------------------------
  interface pt3
    module procedure pt3_m
  end interface pt3
!------------------------------------------------------------------------------
  interface pt2
    module procedure pt2_m
  end interface pt2
!------------------------------------------------------------------------------
  interface pt1
    module procedure pt1_m
  end interface pt1
!------------------------------------------------------------------------------
  interface allocate
    module  procedure allocate_m
    module  procedure allocate_m_ptr4
    module  procedure allocate_m_ptr3
    module  procedure allocate_m_ptr2
  end interface allocate
!------------------------------------------------------------------------------
  interface deallocate
    module  procedure deallocate_m
  end interface deallocate
!------------------------------------------------------------------------------
  interface new_entry
    module  procedure new_entry_m
  end interface new_entry
!------------------------------------------------------------------------------
  interface del_entry
    module  procedure del_entry_m
  end interface del_entry
!------------------------------------------------------------------------------
  interface assignment (=)
    module  procedure assign_m_real
#ifndef TR15581
    module  procedure assign_m
#endif
#ifndef M3ELEMENTAL
    module  procedure assign_m_reals
#ifndef TR15581
    module  procedure assign_ms
#endif
#endif
  end interface ! assignment (=)
!------------------------------------------------------------------------------
  interface plus
    module  procedure sub_m_plus_m
#ifndef M3ELEMENTAL
    module  procedure sub_ms_plus_ms
#endif
  end interface plus
!------------------------------------------------------------------------------
  interface minus
    module  procedure sub_m_minus_m
#ifndef M3ELEMENTAL
    module  procedure sub_ms_minus_ms
#endif
  end interface minus
!------------------------------------------------------------------------------
  interface times
    module  procedure sub_m_times_m
    module  procedure sub_m_times_r
    module  procedure sub_r_times_m
    module  procedure sub_i_times_m
#ifndef M3ELEMENTAL
    module  procedure sub_ms_times_ms
    module  procedure sub_ms_times_r
    module  procedure sub_r_times_ms
    module  procedure sub_i_times_ms
#endif
  end interface times
!------------------------------------------------------------------------------
  interface over
    module  procedure sub_m_over_m
#ifndef M3ELEMENTAL
    module  procedure sub_ms_over_ms
#endif
  end interface
!------------------------------------------------------------------------------
  interface power
    module  procedure sub_m_power_i
#ifndef M3ELEMENTAL
    module  procedure sub_ms_power_i
#endif
  end interface power
!------------------------------------------------------------------------------
  interface sqrt_m
    module  procedure sqrt_m1
#ifndef M3ELEMENTAL
    module  procedure sqrt_ms
#endif
  end interface sqrt_m
!------------------------------------------------------------------------------
  interface max_m
    module  procedure sub_max_r_m
    module  procedure sub_max_r_ms
  end interface max_m
!------------------------------------------------------------------------------
  interface sum
    module  procedure sum_m
    module  procedure sum_ms
  end interface sum
!------------------------------------------------------------------------------
  interface operator (==)
    module procedure equal_m_infos
  end interface ! operator (==)
!------------------------------------------------------------------------------
  interface operator (/=)
    module procedure nequal_m_infos
  end interface
!------------------------------------------------------------------------------
  interface lgtd
    module procedure lgt_m
  end interface lgtd
!------------------------------------------------------------------------------
  interface lgti
    module procedure lgti_m
  end interface lgti
!------------------------------------------------------------------------------
  interface lgtd_ad
    module procedure lgtad_m
  end interface lgtd_ad
!------------------------------------------------------------------------------
  interface lgti_ad
    module procedure lgtiad_m
  end interface lgti_ad
!------------------------------------------------------------------------------
  interface print
    module procedure print_m
    module procedure print_ms
  end interface print
!------------------------------------------------------------------------------
  interface dump
    module procedure dump_m
  end interface dump
!------------------------------------------------------------------------------
  interface retrieve
    module procedure retrieve_m
  end interface retrieve
!------------------------------------------------------------------------------
  interface nwords
    module procedure nwords_m
    module procedure nwords_ms
  end interface nwords
!------------------------------------------------------------------------------
  interface cut
    module procedure cut_m
    module procedure cut_ms
  end interface cut
!------------------------------------------------------------------------------
  interface p_min_m
    module procedure p_min_m
    module procedure p_min_ms
  end interface p_min_m
!------------------------------------------------------------------------------
  interface p_max_m
    module procedure p_max_m
    module procedure p_max_ms
  end interface p_max_m
!------------------------------------------------------------------------------
  interface p_sum_m
    module procedure p_sum_m
    module procedure p_sum_ms
  end interface p_sum_m
!==============================================================================
contains
!==============================================================================
  ELEMENTAL subroutine construct_m (m3)
  !-------------------------
  ! construct table entry(s)
  !-------------------------
  type (t_m) ,intent (out) :: m3
  end subroutine construct_m
!------------------------------------------------------------------------------
  ELEMENTAL subroutine destruct_m (m3)
  !------------------------
  ! destruct table entry(s)
  !------------------------
  type (t_m) ,intent (inout) :: m3
    type (t_m) :: empty
    if (ALLOCATED (m3% ptr)) deallocate (m3% ptr)
    m3 = empty
  end subroutine destruct_m
!------------------------------------------------------------------------------
  function index_m (m3, name)
  !----------------------------------
  ! return the index of a table entry
  !----------------------------------
  type (t_m)        ,intent(in) :: m3 (:)
  character (len=8) ,intent(in) :: name
  integer                       :: index_m
    integer :: i
    index_m = 0
    do i=1,size (m3)
      if (name==m3(i)%i%name) then
        index_m = i
        exit
      endif
    end do
  end function index_m
!------------------------------------------------------------------------------
  function pt4_m (m4,name)
  !-----------------------------------------
  ! return a rank 3 pointer to a table entry
  !-----------------------------------------
  type (t_m)        TARGET ,intent(in) :: m4 (:)
  character (len=8)        ,intent(in) :: name
  real(wp)                 ,pointer    :: pt4_m (:,:,:,:)
    pt4_m => m4(index(m4,name))%ptr
  end function pt4_m
!------------------------------------------------------------------------------
  function pt3_m (m3,name)
  !-----------------------------------------
  ! return a rank 3 pointer to a table entry
  !-----------------------------------------
  type (t_m)        TARGET ,intent(in)  :: m3 (:)
  character (len=8)         ,intent(in) :: name
  real(wp)                  ,pointer    :: pt3_m (:,:,:)
    pt3_m => m3(index(m3,name))%ptr(:,:,:,1)
  end function pt3_m
!------------------------------------------------------------------------------
  function pt2_m (m2,name)
  !-----------------------------------------
  ! return a rank 2 pointer to a table entry
  !-----------------------------------------
  type (t_m)        TARGET ,intent(in) :: m2 (:)
  character (len=8)        ,intent(in) :: name
  real(wp)                 ,pointer    :: pt2_m (:,:)
    pt2_m => m2(index(m2,name))%ptr(:,:,1,1)
  end function pt2_m
!------------------------------------------------------------------------------
  function pt1_m (m1,name)
  !-----------------------------------------
  ! return a rank 1 pointer to a table entry
  !-----------------------------------------
  type (t_m)        TARGET ,intent(in) :: m1 (:)
  character (len=8)        ,intent(in) :: name
  real(wp)                 ,pointer    :: pt1_m (:)
    pt1_m => m1(index(m1,name))%ptr(:,1,1,1)
  end function pt1_m
!------------------------------------------------------------------------------
  ELEMENTAL subroutine allocate_m (m3, nlev)
  !-----------------------
  ! allocate a table entry
  !-----------------------
  type (t_m) ,intent(inout)           :: m3
  integer    ,intent(in)    ,optional :: nlev
    if (present(nlev)) then
      m3%i%lb(3) = 1
      m3%i%ub(3) = nlev
      if (ALLOCATED(m3% ptr)) deallocate (m3% ptr)
    endif
    select case (m3%i% rep)
    case ('sh')
      if (.not.ALLOCATED(m3% ptr)) &
        allocate (m3% ptr(m3%i%nsp,1,m3%i%lb(3):m3%i%ub(3),1))
    case default
      if (.not.ALLOCATED(m3% ptr)) &
        allocate (m3% ptr(m3%i%lb(1):m3%i%ub(1),&
                          m3%i%lb(2):m3%i%ub(2),&
                          m3%i%lb(3):m3%i%ub(3),&
                          m3%i%lb(4):m3%i%ub(4)))
    end select
    m3%i% alloc = .true.
  end subroutine allocate_m
!------------------------------------------------------------------------------
  pure subroutine allocate_m_ptr4 (m3, ptr, nlev)
  !-----------------------------------------------------
  ! allocate a table entry and return a 4D pointer to it
  !-----------------------------------------------------
  type (t_m) TARGET ,intent(inout) :: m3
  real(wp)          ,pointer       :: ptr (:,:,:,:)
  integer ,optional ,intent(in)    :: nlev
    if (present(nlev)) then
      m3%i%lb(3) = 1
      m3%i%ub(3) = nlev
      if (ALLOCATED(m3% ptr)) deallocate (m3% ptr)
    endif
    select case (m3%i% rep)
    case ('sh')
      if (.not.ALLOCATED(m3% ptr)) &
        allocate (m3% ptr(m3%i%nsp,1,m3%i%lb(3):m3%i%ub(3),1))
    case default
      if (.not.ALLOCATED(m3% ptr)) &
        allocate (m3% ptr(m3%i%lb(1):m3%i%ub(1),&
                          m3%i%lb(2):m3%i%ub(2),&
                          m3%i%lb(3):m3%i%ub(3),&
                          m3%i%lb(4):m3%i%ub(4)))
    end select
    m3%i% alloc = .true.
    ptr => m3% ptr
  end subroutine allocate_m_ptr4
!------------------------------------------------------------------------------
  subroutine allocate_m_ptr3 (m3, ptr)
  !--------------------------------------------------
  ! allocate a table entry and return a pointer to it
  !--------------------------------------------------
  type (t_m) TARGET ,intent(inout) :: m3
  real(wp)          ,pointer       :: ptr (:,:,:)
    select case (m3%i% rep)
    case ('sh')
      if (.not.ALLOCATED(m3% ptr)) &
        allocate (m3% ptr(m3%i%nsp,1,m3%i%lb(3):m3%i%ub(3),1))
    case default
      if (.not.ALLOCATED(m3% ptr)) &
        allocate (m3% ptr(m3%i%lb(1):m3%i%ub(1),&
                          m3%i%lb(2):m3%i%ub(2),&
                          m3%i%lb(3):m3%i%ub(3),&
                          m3%i%lb(4):m3%i%ub(4)))
    end select
    m3%i% alloc = .true.
    ptr => m3% ptr (:,:,:,1)
  end subroutine allocate_m_ptr3
!------------------------------------------------------------------------------
  subroutine allocate_m_ptr2 (m3, ptr)
  !-----------------------------------------------------
  ! allocate a table entry and return a 2D pointer to it
  !-----------------------------------------------------
  type (t_m) TARGET ,intent(inout) :: m3
  real(wp)          ,pointer       :: ptr (:,:)
    select case (m3%i% rep)
    case ('sh')
      if (.not.ALLOCATED(m3% ptr)) &
        allocate (m3% ptr(m3%i%nsp,1,m3%i%lb(3):m3%i%ub(3),1))
    case default
      if (.not.ALLOCATED(m3% ptr)) &
        allocate (m3% ptr(m3%i%lb(1):m3%i%ub(1),&
                          m3%i%lb(2):m3%i%ub(2),&
                          m3%i%lb(3):m3%i%ub(3),&
                          m3%i%lb(4):m3%i%ub(4)))
    end select
    m3%i% alloc = .true.
    ptr => m3% ptr(:,:,1,1)
  end subroutine allocate_m_ptr2
!------------------------------------------------------------------------------
  ELEMENTAL pure subroutine deallocate_m (m3)
  !-------------------------
  ! deallocate a table entry
  !-------------------------
  type (t_m) ,intent(inout) :: m3
    if (ALLOCATED(m3% ptr)) deallocate (m3% ptr)
    m3%i% alloc = .false.
  end subroutine deallocate_m
!------------------------------------------------------------------------------
  ELEMENTAL subroutine update_m (m3)
  type (t_m) ,intent(inout) :: m3
  !----------------------------------------
  ! keep array bounds and allocation status
  ! consistent with metainformation
  !----------------------------------------
    if (m3%i% alloc) then
      !--------------------------
      ! field should be allocated
      !--------------------------
      select case (m3%i% rep)
      case ('sh')
        !------------------------
        ! spectral representation
        !------------------------
        if (ALLOCATED(m3% ptr)) then
          if (size(m3%ptr,1)/=m3%i%nsp.or.size(m3%ptr,2)/=1 &
              .or. m3%i%lb(3)/=lbound(m3%ptr,3)             &
              .or. m3%i%ub(3)/=ubound(m3%ptr,3))            &
                 deallocate (m3% ptr)
        endif
        if (.not.ALLOCATED(m3% ptr)) &
          allocate (m3% ptr(m3%i%nsp,1,m3%i%lb(3):m3%i%ub(3),1))
      case default
        !--------------------------
        ! grid-point representation
        !--------------------------
        if (ALLOCATED(m3% ptr)) then
          if (any (m3%i%lb/=lbound(m3%ptr).or.m3%i%ub/=ubound(m3%ptr))) &
            deallocate (m3% ptr)
        endif
        if (.not.ALLOCATED(m3% ptr)) then
          allocate (m3% ptr(m3%i%lb(1):m3%i%ub(1),&
                            m3%i%lb(2):m3%i%ub(2),&
                            m3%i%lb(3):m3%i%ub(3),&
                            m3%i%lb(4):m3%i%ub(4)))
        endif
      end select
    else
      !----------------------------
      ! field should be deallocated
      !----------------------------
      if (ALLOCATED(m3% ptr)) then
        deallocate (m3% ptr)
      endif
    endif
!print *,'update_m:',m3%i% name, m3%i% rep, m3%i% alloc, associated(m3% ptr)
  end subroutine update_m
!------------------------------------------------------------------------------
  subroutine new_entry_m (m3, name, lb, ub, nn, rep, ind, ptr)
  !------------------------------------------------------
  ! establish a new table entry, specify name,
  ! optionally specify bounds, truncation, representation
  ! or obtain index or pointer
  !------------------------------------------------------
  type (t_m)        ,intent(inout) TARGET    :: m3 (:)
  character (len=8) ,intent(in)              :: name
  integer           ,intent(in)    ,optional :: lb(:)
  integer           ,intent(in)    ,optional :: ub(:)
  integer           ,intent(in)    ,optional :: nn
  character(len=2)  ,intent(in)    ,optional :: rep
  integer           ,intent(out)   ,optional :: ind
  real(wp)          ,pointer       ,optional :: ptr (:,:,:)
    integer :: i
    i = index (m3, name)
    if (i/=0) call finish('new_entry_m','name is already used: '//name)
    i = index (m3, '        ')
    if (i==0) call finish('new_entry_m','table is full')
    m3(i)%i% name = name
    if(present(lb))  m3(i)%i% lb (1:size(lb)) = lb
    if(present(ub))  m3(i)%i% ub (1:size(ub)) = ub
    if(present(nn))  m3(i)%i% nn  = nn
    if(present(rep)) m3(i)%i% rep = rep
    m3(i)%i% nsp = nsp_nn (m3(i)%i% nn)
    call allocate (m3(i))
    if(present(ind)) ind = i
    if(present(ptr)) ptr => m3(i)% ptr (:,:,:,1)
  end subroutine new_entry_m
!==============================================================================
  !-----------
  ! Assignment
  !-----------
!------------------------------------------------------------------------------
  subroutine del_entry_m (m3, name)
  !---------------------
  ! delete a table entry
  !---------------------
  type (t_m)       ,intent(inout) :: m3 (:)
  character (len=8) ,intent(in)    :: name
    call destruct (m3 (index(m3,name)))
  end subroutine del_entry_m
!------------------------------------------------------------------------------
#ifndef TR15581
  ELEMENTAL subroutine assign_m (y, x)
  !------------------------------------------------
  ! overwrite intrinsic assignment routine: m3 = m3
  !------------------------------------------------
  type (t_m) ,intent(inout) :: y
  type (t_m) ,intent(in)    :: x
    y% i   = x% i
    call update (y)
    if (y%i%alloc) y% ptr = x% ptr
  end subroutine assign_m
#endif
!------------------------------------------------------------------------------
  ELEMENTAL subroutine assign_m_real (y, x)
  !-------------------------------
  ! assign a real to table entries
  !-------------------------------
  type (t_m) ,intent(inout) :: y
  real(wp)   ,intent(in)    :: x
    if (y%i% alloc) y% ptr = x
  end subroutine assign_m_real
!------------------------------------------------------------------------------
!
! NAG f95 beta (112) bug(?):
!
!Implicit pointer assignment from X1 illegal in elemental procedure M3_PLUS_M
!Implicit pointer assignment from X1 illegal in elemental procedure M3_MINUS_M
!Implicit pointer assignment from X1 illegal in elemental procedure M3_TIMES_R
!Implicit pointer assignment from X2 illegal in elemental procedure R_TIMES_M
!Implicit pointer assignment from X1 illegal in elemental procedure M3_POWER_I
!
#ifdef TR15581
#  define ASSIGN(Y,X) Y = X
#else
#  define ASSIGN(Y,X) call assign_m (Y, X)
#endif
!==============================================================================
  subroutine p_min_m (x)
  type (t_m) ,intent(inout) :: x
    if (x%i% alloc .and..not. x%i% ref) then
      x% ptr = p_min (x% ptr)
    endif
  end subroutine p_min_m
!------------------------------------------------------------------------------
  subroutine p_min_ms (x)
  type (t_m) ,intent(inout) :: x(:)
    integer :: l
    do l=1,size(x)
      if (x(l)%i% alloc .and..not. x(l)%i% ref) then
        x(l)% ptr = p_min (x(l)% ptr)
      endif
    end do
  end subroutine p_min_ms
!==============================================================================
  subroutine p_max_m (x)
  type (t_m) ,intent(inout) :: x
    if (x%i% alloc .and..not. x%i% ref) then
      x% ptr = p_max (x% ptr)
    endif
  end subroutine p_max_m
!------------------------------------------------------------------------------
  subroutine p_max_ms (x)
  type (t_m) ,intent(inout) :: x(:)
    integer :: l
    do l=1,size(x)
      if (x(l)%i% alloc .and..not. x(l)%i% ref) then
        x(l)% ptr = p_max (x(l)% ptr)
      endif
    end do
  end subroutine p_max_ms
!==============================================================================
  subroutine p_sum_m (x)
  type (t_m) ,intent(inout) :: x
    if (x%i% alloc .and..not. x%i% ref) then
      x% ptr = p_sum (x% ptr)
    endif
  end subroutine p_sum_m
!------------------------------------------------------------------------------
  subroutine p_sum_ms (x)
  type (t_m) ,intent(inout) :: x(:)
    integer :: l
    do l=1,size(x)
!     workaround bug in IBM xlf 7.1.1.4 compiler
!     if (x(l)%i% alloc .and..not. x(l)%i% ref) x(l)% ptr = p_sum (x(l)% ptr)
      if (x(l)%i% alloc .and..not. x(l)%i% ref) then
        x(l)% ptr = p_sum (x(l)% ptr)
      endif
    end do
  end subroutine p_sum_ms
!==============================================================================
  !----------------------
  ! arithmetic operations
  !----------------------
!
! There is no test on corresponding representation
!
!------------------------------------------------------------------------------
  ELEMENTAL subroutine sub_m_plus_m (y, x1, x2)
  type (t_m) ,intent(inout) :: y
  type (t_m) ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2
    ASSIGN(y, x1)
    if (y%i% alloc .and. x2%i% alloc .and..not.y%i% ref) &
      y% ptr = y% ptr + x2% ptr
  end subroutine sub_m_plus_m
!------------------------------------------------------------------------------
  ELEMENTAL subroutine sub_m_minus_m (y, x1, x2)
  type (t_m) ,intent(inout) :: y
  type (t_m) ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2
    ASSIGN(y, x1)
    if (y%i% alloc .and. x2%i% alloc .and..not.y%i% ref) &
      y% ptr = y% ptr - x2% ptr
  end subroutine sub_m_minus_m
!------------------------------------------------------------------------------
  ELEMENTAL subroutine sub_m_times_m (y, x1, x2)
  type (t_m) ,intent(inout) :: y
  type (t_m) ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2
    integer :: i,j,l
    if (x1% i% rep == 'vc') then
      !------------------------------------
      ! vertical correlation matrix * field
      !------------------------------------
      y%i        = x2%i
      y%i% alloc = x1%i% alloc .and. x2%i% alloc
      call update (y)
      if (y%i% alloc) then
        select case (x2% i% rep)
        case ('sh')
          do i=1,size(x2% ptr,1)
            y% ptr(i,1,:,1) = matmul (x1% ptr(1,:,:,1), x2% ptr(i,1,:,1))
          end do
        case ('gg')
          do l=1,size(x2% ptr,4)
            do j=1,size(x2% ptr,2)
              do i=1,size(x2% ptr,1)
                y% ptr(i,j,:,l) = matmul (x1% ptr(1,:,:,1), x2% ptr(i,j,:,l))
              end do
            end do
          end do
        end select
      endif
    else
      !--------------
      ! field * field
      !--------------
      if (x1%i% ref) then
        y = x1
      else
        y%i        = x1%i
        y%i% alloc = x1%i% alloc .and. x2%i% alloc
        call update (y)
        if (y%i% alloc) y% ptr = x1% ptr * x2% ptr
      endif
    endif
  end subroutine sub_m_times_m
!------------------------------------------------------------------------------
  ELEMENTAL subroutine sub_m_times_r (y, x1, x2)
  type (t_m) ,intent(inout) :: y
  type (t_m) ,intent(in)    :: x1
  real(wp)   ,intent(in)    :: x2
    ASSIGN(y, x1)
    if (y%i% alloc .and..not.y%i% ref) y% ptr = y% ptr * x2
  end subroutine sub_m_times_r
!------------------------------------------------------------------------------
  ELEMENTAL subroutine sub_r_times_m (y, x1, x2)
  type (t_m) ,intent(inout) :: y
  real(wp)   ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2
    ASSIGN(y, x2)
    if (y%i% alloc .and..not.y%i% ref) y% ptr = y% ptr * x1
  end subroutine sub_r_times_m
!------------------------------------------------------------------------------
  ELEMENTAL subroutine sub_i_times_m (y, x1, x2)
  type (t_m) ,intent(inout) :: y
  integer    ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2
    ASSIGN(y, x2)
    if (y%i% alloc .and..not.y%i% ref) y% ptr = y% ptr * x1
  end subroutine sub_i_times_m
!------------------------------------------------------------------------------
  ELEMENTAL subroutine sub_m_over_m (y, x1, x2)
  type (t_m) ,intent(inout) :: y
  type (t_m) ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2
    y%i        = x1%i
    y%i% alloc = x1%i% alloc .and. x2%i% alloc
    call update (y)
    if (y%i% alloc .and..not.y%i% ref) y% ptr = x1% ptr / x2% ptr
  end subroutine sub_m_over_m
!------------------------------------------------------------------------------
  ELEMENTAL subroutine sub_m_power_i (y, x1, i)
  type (t_m) ,intent(inout) :: y
  type (t_m) ,intent(in)    :: x1
  integer    ,intent(in)    :: i
    ASSIGN(y, x1)
    if (y%i% alloc .and..not.y%i% ref) y% ptr = y% ptr ** i
  end subroutine sub_m_power_i
!------------------------------------------------------------------------------
  ELEMENTAL subroutine sqrt_m1 (y, x)
  type (t_m) ,intent(inout) :: y
  type (t_m) ,intent(in)    :: x
    ASSIGN(y, x)
    if (y%i% alloc .and..not.y%i% ref) y% ptr = sqrt (y% ptr)
  end subroutine sqrt_m1
!------------------------------------------------------------------------------
  function sum_m (x) result (y)
  type (t_m) ,intent(in) :: x
  real(wp)               :: y
    y = 0._wp
    if (x%i% alloc) y = sum (x% ptr)
  end function sum_m
!------------------------------------------------------------------------------
  function sum_ms (x) result (y)
  type (t_m) ,intent(in) :: x (:)
  real(wp)               :: y
    integer :: i
    y = 0._wp
    do i=1,size(x)
      y = y + sum (x(i))
    end do
  end function sum_ms
!==============================================================================
  subroutine cut_m (y, x, lb, ub)
  !---------------
  ! Cutout regions
  !---------------
  type (t_m) ,intent(inout) :: y
  type (t_m) ,intent(in)    :: x
  integer    ,intent(in)    :: lb(2)
  integer    ,intent(in)    :: ub(2)
    y% i          = x% i
    y% i% lb(1:2) = lb
    y% i% ub(1:2) = ub
    call update (y)
    if (y%i%alloc) y% ptr = x% ptr (lb(1):ub(1),lb(2):ub(2),:,:)
  end subroutine cut_m
!------------------------------------------------------------------------------
  subroutine cut_ms (y, x, lb, ub)
  !---------------
  ! Cutout regions
  !---------------
  type (t_m) ,intent(inout) :: y(:)
  type (t_m) ,intent(in)    :: x(:)
  integer    ,intent(in)    :: lb(2)
  integer    ,intent(in)    :: ub(2)
    integer :: i
    do i=1,size(y)
      call cut (y(i), x(i), lb, ub)
    end do
  end subroutine cut_ms
!==============================================================================
  subroutine lgt_m (m3,nn,norm)
  !---------------------------------------
  ! direct legendre transformation routine
  !---------------------------------------
  type (t_m) ,intent(inout)           :: m3
  integer    ,intent(in)    ,optional :: nn
  integer    ,intent(in)    ,optional :: norm
    real(wp), allocatable :: sh (:,:)
    if (m3%i%rep == 'gg') then
      if (present(nn)) then
        m3%i% nn  = nn
      endif
      if (m3%i%alloc) then
        m3%i% nsp = nsp_nn (m3%i% nn)
        allocate (sh (m3%i% nsp,size(m3% ptr,3)))
        call lgtd (sh, m3% ptr(:,:,:,1), m3%i% nn, norm=norm)
        m3%i% rep = 'sh'
        call update (m3)
        m3% ptr(:,1,:,1) = sh
      endif
    endif
  end subroutine lgt_m
!------------------------------------------------------------------------------
  subroutine lgti_m (m3,norm)
  !----------------------------------------
  ! inverse legendre transformation routine
  !----------------------------------------
  type (t_m) ,intent(inout)           :: m3
  integer    ,intent(in)    ,optional :: norm
    real(wp), allocatable :: gg (:,:,:)
    integer :: n (3)
    if (m3%i%rep == 'sh' .and. m3%i%alloc) then
      n = m3%i% ub(1:3) - m3%i% lb(1:3) + 1
      allocate (gg (n(1),n(2),n(3)))
      call lgti (gg, m3% ptr(:,1,:,1), m3%i% nn, norm=norm)
      m3%i% rep = 'gg'
      call update (m3)
      m3% ptr(:,:,:,1) = gg
    endif
  end subroutine lgti_m
!------------------------------------------------------------------------------
  subroutine lgtiad_m (m3,nn,norm)
  !------------------------------------------------
  ! adjoint inverse legendre transformation routine
  ! (gg -> sh)
  !------------------------------------------------
  type (t_m) ,intent(inout)           :: m3
  integer    ,intent(in)    ,optional :: nn
  integer    ,intent(in)    ,optional :: norm
    real(wp), allocatable :: sh (:,:)
    if (m3%i%rep == 'gg') then
      if (present(nn)) then
        m3%i% nn  = nn
        m3%i% nsp = nsp_nn (m3%i% nn)
      endif
      if (m3%i%alloc) then
        allocate (sh (m3%i% nsp,size(m3% ptr,3)))
        call lgti_ad (m3% ptr(:,:,:,1), sh, m3%i% nn, norm=norm)
        m3%i% rep = 'sh'
        call update (m3)
        m3% ptr(:,1,:,1) = sh
      endif
    endif
  end subroutine lgtiad_m
!------------------------------------------------------------------------------
  subroutine lgtad_m (m3,norm)
  !-----------------------------------------------
  ! adjoint direct legendre transformation routine
  ! (sh -> gp)
  !-----------------------------------------------
  type (t_m) ,intent(inout)           :: m3
  integer    ,intent(in)    ,optional :: norm
    real(wp), allocatable :: gg (:,:,:)
    integer :: n (3)
    if (m3%i%rep == 'sh' .and. m3%i%alloc) then
      n = m3%i% ub(1:3) - m3%i% lb(1:3) + 1
      allocate (gg (n(1),n(2),n(3)))
      call lgtd_ad (m3% ptr(:,1,:,1), gg, m3%i% nn, norm=norm)
      m3%i% rep = 'gg'
      call update (m3)
      m3% ptr(:,:,:,1) = gg
    endif
  end subroutine lgtad_m
!==============================================================================
  subroutine print_ms (m3, iunit, intend, verbose)
  type (t_m)       ,intent(in)           :: m3 (:) ! table to print
  integer          ,intent(in) ,optional :: iunit  ! unit   (default=6)
  character(len=*) ,intent(in) ,optional :: intend ! intend string ('')
  logical          ,intent(in) ,optional :: verbose
    integer i
    do i=1,size(m3)
      call print (m3(i), iunit, intend, verbose)
    end do
  end subroutine print_ms
!------------------------------------------------------------------------------
  subroutine print_m (m3, iunit, intend, verbose)
  type (t_m)       ,intent(in)           :: m3     ! table entry to print
  integer          ,intent(in) ,optional :: iunit  ! unit   (default=6)
  character(len=*) ,intent(in) ,optional :: intend ! intend string ('')
  logical          ,intent(in) ,optional :: verbose
    integer           :: k, iu, n, l, m
    character(len=32) :: c
    logical           :: v
    real(wp)          :: mi, ma, sm, mm
    character(len=2)  :: cm
    !--------------------
    ! optional parameters
    !--------------------
    n=0
    c=''
    if (present(intend)) then
      c = intend
      n = len_trim(c)
    endif
    iu = 6; if (present(iunit)) iu = iunit
    v = .false.; if (present(verbose)) v = verbose
    if (m3% i% alloc) then
      if (.not. ALLOCATED(m3% ptr)) then
        write (iu,"(a,a,a,a)") c(:n+1), m3%i% name ,m3%i%rep,&
          'ERROR: is not associated !!!'
      else
        m  = size(m3% ptr)
        call mi_ma (mi, ma, sm, size(m3% ptr), m3% ptr, m)
        mi = p_min (mi)
        ma = p_max (ma)
        m  = p_sum (m )
        sm = p_sum (sm)
        if (m > 0) sm = sm / m
        if (dace% lpio) then
          write (iu,"(a,a,a,'(',a,')',3g11.3,5i4,i3,':',i3)") c(:n+1),      &
                  m3%i% name ,m3%i%rep,': min,max=',mi,sm,ma,               &
                  m3%i% code, m3%i% table, m3%i% dis, m3%i% cat, m3%i% num, &
                  lbound(m3% ptr,3),ubound(m3% ptr,3)
        endif
        if(v .and. size(m3% ptr,3)>1) then
          mm = max (abs (mi), abs (ma))
          if (mm == 0._wp) mm = -HUGE(0._wp)
          do k   = lbound(m3% ptr,3),ubound(m3% ptr,3)
            mi =  huge(mi)
            ma = -huge(ma)
            sm = 0._wp
            m  = size(m3% ptr(:,:,k,:))
!           do l = lbound(m3% ptr,4),ubound(m3% ptr,4)
!!CDIR COLLAPSE
!              mi = min (mi, minval(m3% ptr(:,:,k,l)))
!!CDIR COLLAPSE
!              ma = max (ma, maxval(m3% ptr(:,:,k,l)))
!!CDIR COLLAPSE
!              sm = sm +     sum   (m3% ptr(:,:,k,l))
!           end do
            call mi_ma (mi, ma, sm, size(m3% ptr(:,:,k,:)), m3% ptr(:,:,k,:), m)
            mi = p_min (mi)
            ma = p_max (ma)
            m  = p_sum (m )
            sm = p_sum (sm) / m
            if      (abs (ma) == mm) then
               cm = " >"
            else if (abs (mi) == mm) then
               cm = " <"
            else
               cm = "  "
            end if
            if (dace% lpio)                       &
              write (iu,"(a,a,a,i3,4x,3g11.3,a)") &
              c(:n+1), m3%i% name ,'  k=',k,mi,sm,ma, trim (cm)
          end do
        endif
      endif
!   else
!     if (v .and. dace% lpio)                  &
!       write(iu,"(a,a,a,'(',a,')')")c(:n+1),  &
!       m3%i% name,m3%i%rep,': not ALLOCATED!'
    endif
  contains

    subroutine mi_ma (mi, ma, sm, m, ptr, n)
    !--------------------------------------------------------------
    ! derive minval, maxval; collapse loop indices in ptr (for SX9)
    !--------------------------------------------------------------
    real(wp) ,intent(out) :: mi
    real(wp) ,intent(out) :: ma
    real(wp) ,intent(out) :: sm
    integer  ,intent(in)  :: m
    real(wp) ,intent(in)  :: ptr(m)
    integer  ,intent(out) :: n

      logical :: valid (m)

      mi    = minval (ptr)
      ma    = maxval (ptr)
      valid = (abs (ptr) < huge(1._wp))
      sm    = sum    (ptr, mask=valid)
      n     = count  (valid)
    end subroutine mi_ma

  end subroutine print_m
!------------------------------------------------------------------------------
  subroutine dump_m (m3, iunit)
  type (t_m) ,intent(in)  :: m3(:)  ! table to dump
  integer    ,intent(in)  :: iunit  ! unit
    integer           :: i
    character(len=11) :: formatted
    inquire(iunit,form=formatted)
    if (formatted=='FORMATTED') then
      write (iunit,*) m3% i
      do i=1,size(m3)
        if (m3(i)%i% alloc) write (iunit,*) m3(i)% ptr
      end do
    else
      write (iunit) m3% i
      do i=1,size(m3)
        if (m3(i)%i% alloc) write (iunit) m3(i)% ptr
      end do
    endif
  end subroutine dump_m
!------------------------------------------------------------------------------
  subroutine retrieve_m (m3, iunit)
  type (t_m)   ,intent(out) :: m3(:)  ! table to retrieve
  integer      ,intent(in)  :: iunit  ! unit
    integer           :: i
    character(len=11) :: formatted
    inquire(iunit,form=formatted)
    if (formatted=='FORMATTED') then
      read (iunit,*) m3% i
      call update (m3)
      do i=1,size(m3)
        if (m3(i)%i% alloc) read (iunit,*) m3(i)% ptr
      end do
    else
      read (iunit) m3% i
      call update (m3)
      do i=1,size(m3)
        if (m3(i)%i% alloc) read (iunit) m3(i)% ptr
      end do
    endif
  end subroutine retrieve_m
!------------------------------------------------------------------------------
  elemental function equal_m_infos (i1, i2) result (equal)
  type (t_mi) ,intent(in) :: i1, i2
  logical                 :: equal
    equal = ( all(i1% lb   == i2% lb)  .and. all(i1% ub  ==  i2% ub) .and. &
                  i1% name == i2% name .and. i1% alloc .eqv. i2% alloc)
  end function equal_m_infos
!------------------------------------------------------------------------------
  elemental function nequal_m_infos (i1, i2) result (nequal)
  type (t_mi) ,intent(in) :: i1, i2
  logical                 :: nequal
    nequal = ( any(i1% lb   /= i2% lb)  .or. any(i1% ub  /=  i2% ub) .or. &
                   i1% name /= i2% name .or. i1% alloc .neqv. i2% alloc)
  end function nequal_m_infos
!------------------------------------------------------------------------------
  pure function nwords_m (m3) result (n)
  type (t_m) ,intent(in) :: m3
  integer(i8)            :: n
    n = 0_i8
    if (m3%i%alloc) n = size (m3%ptr)
  end function nwords_m
!------------------------------------------------------------------------------
  pure function nwords_ms (m3s) result (n)
  type (t_m) ,intent(in) :: m3s(:)
  integer(i8)            :: n
    integer :: i
    n = 0_i8
    do i=1,size(m3s)
      if (m3s(i)%i%alloc) n = n + size (m3s(i)%ptr)
    end do
  end function nwords_ms
!------------------------------------------------------------------------------
  subroutine sub_max_r_ms (y, x1, x2)
  type (t_m) ,intent(inout) :: y (:)
  real(wp)   ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2(:)
    integer :: i
    do i=1,size(y)
      call sub_max_r_m (y(i),x1,x2(i))
    end do
  end subroutine sub_max_r_ms
!------------------------------------------------------------------------------
  subroutine sub_max_r_m (y, x1, x2)
  type (t_m) ,intent(inout) :: y
  real(wp)   ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2
    ASSIGN(y, x2)
    if (y%i% alloc .and..not.y%i% ref) y% ptr = max (y% ptr, x1)
  end subroutine sub_max_r_m
!------------------------------------------------------------------------------
#ifndef M3ELEMENTAL
  !-------------------------------------------------------------
  ! the following routines are required only if elemental doesnt
  ! work on type t_m
  !-------------------------------------------------------------
  subroutine update_ms (m3)
  type (t_m) ,intent(inout) :: m3(:)
    integer :: i
    do i=1,size(m3)
      call update_m (m3(i))
    end do
  end subroutine update_ms
!------------------------------------------------------------------------------
  subroutine destruct_ms (m3)
  !------------------------
  ! destruct table entry(s)
  !------------------------
  type (t_m) ,intent (inout) :: m3(:)
    integer :: i
    do i=1,size(m3)
      call destruct_m (m3(i))
    end do
  end subroutine destruct_ms
!------------------------------------------------------------------------------
  subroutine sub_ms_plus_ms (y, x1, x2)
  type (t_m) ,intent(inout) :: y (:)
  type (t_m) ,intent(in)    :: x1(:)
  type (t_m) ,intent(in)    :: x2(:)
    integer :: i
    do i=1,size(y)
      call sub_m_plus_m(y(i),x1(i),x2(i))
    end do
  end subroutine sub_ms_plus_ms
!------------------------------------------------------------------------------
  subroutine sub_ms_minus_ms (y, x1, x2)
  type (t_m) ,intent(inout) :: y (:)
  type (t_m) ,intent(in)    :: x1(:)
  type (t_m) ,intent(in)    :: x2(:)
    integer :: i
    do i=1,size(y)
      call sub_m_minus_m(y(i),x1(i),x2(i))
    end do
  end subroutine sub_ms_minus_ms
!------------------------------------------------------------------------------
  subroutine sub_ms_times_ms (y, x1, x2)
  type (t_m) ,intent(inout) :: y (:)
  type (t_m) ,intent(in)    :: x1(:)
  type (t_m) ,intent(in)    :: x2(:)
    integer :: i
    do i=1,size(y)
      call sub_m_times_m(y(i),x1(i),x2(i))
    end do
  end subroutine sub_ms_times_ms
!------------------------------------------------------------------------------
  subroutine sub_ms_times_r (y, x1, x2)
  type (t_m) ,intent(inout) :: y (:)
  type (t_m) ,intent(in)    :: x1(:)
  real(wp)   ,intent(in)    :: x2
    integer :: i
    do i=1,size(y)
      call sub_m_times_r (y(i),x1(i),x2)
    end do
  end subroutine sub_ms_times_r
!------------------------------------------------------------------------------
  subroutine sub_r_times_ms (y, x1, x2)
  type (t_m) ,intent(inout) :: y (:)
  real(wp)   ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2(:)
    integer :: i
    do i=1,size(y)
      call sub_r_times_m (y(i),x1,x2(i))
    end do
  end subroutine sub_r_times_ms
!------------------------------------------------------------------------------
  subroutine sub_i_times_ms (y, x1, x2)
  type (t_m) ,intent(inout) :: y (:)
  integer    ,intent(in)    :: x1
  type (t_m) ,intent(in)    :: x2(:)
    integer :: i
    do i=1,size(y)
      call sub_i_times_m (y(i),x1,x2(i))
    end do
  end subroutine sub_i_times_ms
!------------------------------------------------------------------------------
  subroutine sub_ms_over_ms (y, x1, x2)
  type (t_m) ,intent(inout) :: y (:)
  type (t_m) ,intent(in)    :: x1(:)
  type (t_m) ,intent(in)    :: x2(:)
    integer :: i
    do i=1,size(y)
      call sub_m_over_m(y(i),x1(i),x2(i))
    end do
  end subroutine sub_ms_over_ms
!------------------------------------------------------------------------------
  subroutine sub_ms_power_i (y, x1, i)
  type (t_m) ,intent(inout) :: y (:)
  type (t_m) ,intent(in)    :: x1(:)
  integer     ,intent(in)    :: i
    integer :: j
    do j=1,size(y)
      call sub_m_power_i (y(j), x1(j), i)
    end do
  end subroutine sub_ms_power_i
!------------------------------------------------------------------------------
  subroutine sqrt_ms (y, x)
  type (t_m) ,intent(inout) :: y (:)
  type (t_m) ,intent(in)    :: x (:)
    integer :: j
    do j=1,size(y)
      call sqrt_m1 (y(j), x(j))
    end do
  end subroutine sqrt_ms
!------------------------------------------------------------------------------
#ifndef TR15581
  subroutine assign_ms (y, x)
  !------------------------------------------------
  ! overwrite intrinsic assignment routine: m3 = m3
  !------------------------------------------------
  type (t_m) ,intent(inout) :: y(:)
  type (t_m) ,intent(in)    :: x(:)
    integer :: i
    do i=1,size(y)
      y(i) = x(i)
    end do
  end subroutine assign_ms
#endif
#endif
!------------------------------------------------------------------------------
  subroutine assign_m_reals (y, x)
  !-------------------------------
  ! assign a real to table entries
  !-------------------------------
  type (t_m) ,intent(inout) :: y(:)
  real(wp)    ,intent(in)    :: x
    integer :: i
    do i=1,size(y)
      call assign_m_real (y(i), x)
    end do
  end subroutine assign_m_reals
!------------------------------------------------------------------------------
end module mo_memory
!==============================================================================
#ifndef TR15581
pure subroutine delete_m (m,n)
use mo_memory, only: t_m
implicit none
integer  , intent(in)    :: n
type(t_m), intent(inout) :: m(n)
  integer :: i
  do i=1,n
    if (ALLOCATED (m(i)% ptr)) then
      deallocate (m(i)% ptr)
    endif
  end do
end subroutine delete_m
#endif
!==============================================================================
