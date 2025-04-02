!
!+ DWD modifications to the ropp_io library (for GNSS radio occultations)
!
! $Id$
!
MODULE ropp_io_dwd
!
! Description:
!   DWD modifications to the ropp_io library.
!   Handle additional components in ROprof% L1btype:
!     r_gns: coordinates of GPS satellite
!     v_gns: velocity    of GPS satellite
!     r_leo: coordinates of LEO satellite
!     v_leo: velocity    of LEO satellite
!   Activated only if compiler directive -DROPP is set.
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
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2006-2007  original code
!------------------------------------------------------------------------------
#ifndef ROPP

! empty module

#else
  !========================================
  ! use entities of original module ropp_io
  !========================================
  use ropp_io        ,only :                 &!
       ropp_io_init_orig  => ropp_io_init,   &! initialize data type ROprof
       ropp_io_free_orig  => ropp_io_free,   &! free components of   ROprof
                             ROprof,         &! ROPP RO-profile data type
                             DT7type,        &! ROPP time data type
                             L1btype,        &! component of ", level 1b data
       ropp_io_write_orig => ropp_io_write,  &! write data type ROprof
                             ropp_io_addvar, &! add user defined variables
                             ropp_io_occid,  &! derive occultation id
                             wp               ! working precision kind
  use mo_time        ,only : t_time,         &! 3dvar time data type
                             iyyyy,          &! years   from t_time
                             imm,            &! months  from t_time
                             idd,            &! days    from t_time
                             ihh,            &! hours   from t_time
                             imi,            &! minutes from t_time
                             iss              ! seconds from t_time
  use mo_run_params   ,only: data,           &! pathto ROPIC_Impact_Heights.dat
                             path_file        ! concatenate path+file name
  use mo_exception    ,only: finish           ! abort in case of error
  use mo_fortran_units,only: get_unit_number,&! get free unit number
                             return_unit_number
  implicit none

  !================
  ! Public entities
  !================
  private
  public  :: ROprof         ! data type ROprof, modified in ROPP library
  public  :: ropp_io_init   ! modified routine for  DWD extensions
  public  :: ropp_io_free   ! modified routine for  DWD extensions
  public  :: ropp_io_write  ! original routine from ROPP library
  public  :: ropp_io_addvar ! original routine from ROPP library
  public  :: ropp_io_occid  ! derive occultation id
  public  :: assignment(=)  ! ROPP time = 3dvar time
  !-------------------------------------------------
  ! parameters for the ROPIC intercomparison project
  !-------------------------------------------------
  public  ::  n_ropic_1b    ! number    of levels for level 1b data
  public  ::  n_ropic_2     ! number    of levels for level 2  data
  public  :: dz_ropic_2     ! increment of levels for level 2  data
  public  :: impact_heights ! read ROPIC impact heights

  !=====================================
  ! modified interface for ropp_io_init:
  ! use DWD extensions for level 1b data
  !=====================================
  interface ropp_io_init

      subroutine ropp_io_init_l1atype(var, n)
        use ropp_io_types
        type(L1atype), intent(inout) :: var
        integer                      :: n
      end subroutine ropp_io_init_l1atype

     subroutine ropp_io_init_l2atype(var, n)
       use ropp_io_types
       type(L2atype), intent(inout) :: var
       integer                      :: n
     end subroutine ropp_io_init_l2atype

      subroutine ropp_io_init_l2btype(var, n)
        use ropp_io_types
        type(L2btype), intent(inout) :: var
        integer                     :: n
      end subroutine ropp_io_init_l2btype

!     subroutine ropp_io_init_l2ctype(var, n)
!       use ropp_io_types
!       type(L2ctype), intent(inout) :: var
!       integer                      :: n
!     end subroutine ropp_io_init_l2ctype

!     subroutine ropp_io_init_l2dtype(var, n)
!       use ropp_io_types
!       type(L2dtype), intent(inout) :: var
!       integer                      :: n
!     end subroutine ropp_io_init_l2dtype

!     subroutine ropp_io_init_vlist(var)
!       use ropp_io_types
!       type(Vlisttype), intent(inout) :: var
!     end subroutine ropp_io_init_vlist

    module procedure ropp_io_init_dwd

    module procedure ropp_io_init_l1btype_dwd

  end interface ropp_io_init

  interface ropp_io_free
    module procedure ropp_io_free_dwd
  end interface ropp_io_free

  interface assignment(=)
    module procedure dt7_time   ! DT7type = t_time
  end interface assignment(=)

  !=================================================
  ! parameters for the ROPIC intercomparison project
  !=================================================
  integer  ::  n_ropic_1b = 247     ! number    of levels for level 1b data
  integer  ::  n_ropic_2  = 226     ! number    of levels for level 2  data
  real(wp) :: dz_ropic_2  = 200._wp ! increment of levels for level 2  data (m)

!==============================================================================
contains
!==============================================================================
  subroutine ropp_io_init_dwd (ro_data, n_lev_1a, &
                                        n_lev_1b, &
                                        n_lev_2a, &
                                        n_lev_2b, &
                                        n_lev_2c, &
                                        n_lev_2d  )
  !----------------------------------------------------------------------------
  ! allocate components of derived data type ROprof
  ! DWD: additional components for level 1b
  !      (satellte coordinates and velocities)
  !----------------------------------------------------------------------------
  type(ROprof) ,intent(inout) :: ro_data  ! GPS RO DATA (ROPP data type)
  integer      ,intent(in)    :: n_lev_1a ! size of profile for level 1a data
  integer      ,intent(in)    :: n_lev_1b ! size of profile for level 1b data
  integer      ,intent(in)    :: n_lev_2a ! size of profile for level 2a data
  integer      ,intent(in)    :: n_lev_2b ! size of profile for level 2b data
  integer      ,intent(in)    :: n_lev_2c ! size of profile for level 2c data
  integer      ,intent(in)    :: n_lev_2d ! size of profile for level 2d data

    !----------------------
    ! call original routine
    !----------------------
    call ropp_io_init_orig (ro_data,              &
      n_lev_1a,             &
      n_lev_1b,             &
      n_lev_2a,             &
      n_lev_2b,             &
      n_lev_2c,             &
      n_lev_2d)

    !------------------------
    ! allocate DWD extensions
    !------------------------
    call alloc_l1b_extensions (ro_data% Lev1b)

  end subroutine ropp_io_init_dwd
!------------------------------------------------------------------------------
  subroutine ropp_io_init_l1btype_dwd(var, n)
  type(L1btype), intent(inout) :: var
  integer      , intent(in)    :: n
    !----------------------
    ! call original routine
    !----------------------
    call ropp_io_init_orig (var, n)
    !------------------------
    ! allocate DWD extensions
    !------------------------
    call alloc_l1b_extensions (var)

  end subroutine ropp_io_init_l1btype_dwd
!------------------------------------------------------------------------------
  subroutine alloc_l1b_extensions (l1b)
  type(L1btype), intent(inout) :: l1b
  !------------------------
  ! allocate DWD extensions
  !------------------------
    allocate (l1b% r_gns (l1b% npoints, 3))
    allocate (l1b% v_gns (l1b% npoints, 3))
    allocate (l1b% r_leo (l1b% npoints, 3))
    allocate (l1b% v_leo (l1b% npoints, 3))

    l1b% r_gns = -999.9_wp
    l1b% v_gns = -999.9_wp
    l1b% r_leo = -999.9_wp
    l1b% v_leo = -999.9_wp

  end subroutine alloc_l1b_extensions
!------------------------------------------------------------------------------
  subroutine dealloc_l1b_extensions (l1b)
  type(L1btype), intent(inout) :: l1b
  !--------------------------
  ! deallocate DWD extensions
  !--------------------------
    if (associated(l1b% r_gns)) deallocate (l1b% r_gns)
    if (associated(l1b% v_gns)) deallocate (l1b% v_gns)
    if (associated(l1b% r_leo)) deallocate (l1b% r_leo)
    if (associated(l1b% v_leo)) deallocate (l1b% v_leo)

  end subroutine dealloc_l1b_extensions
!------------------------------------------------------------------------------
  subroutine ropp_io_free_dwd (ROdata)
  type(ROprof), intent(inout) :: ROdata
    call ropp_io_free_orig      (ROdata)
    call dealloc_l1b_extensions (ROdata% Lev1b)
  end subroutine ropp_io_free_dwd
!------------------------------------------------------------------------------
  subroutine ropp_io_write (ROdata, file, path, type, append, rec, ierr)
  !-----------------------------------------------------------
  ! Writing a single profile of ROprof data (netCDF and ASCII)
  !-----------------------------------------------------------
  type(ROprof),       intent(inout) :: ROdata
  character(len = *), optional      :: file
  character(len = *), optional      :: path
  character(len = *), optional      :: type
  logical,            optional      :: append
  integer,            optional      :: rec
  integer,            optional      :: ierr
    !----------------------
    ! handle DWD extensions
    !----------------------
    if (associated(rodata% Lev1b% r_gns))      &
      call ropp_io_addvar(                     &
         rodata,                               &
         'r_gns_1b',                           &
         'GNSS transmitter position level 1b', &
         'meters',                             &
         (/-43000000._wp, 43000000._wp/),      &
         rodata% Lev1b% r_gns)

    if (associated(rodata% Lev1b% r_leo))      &
      call ropp_io_addvar(                     &
         rodata,                               &
         'r_leo_1b',                           &
         'LEO transmitter position level 1b',  &
         'meters',                             &
         (/ -10000000._wp, 10000000._wp/),     &
         rodata% Lev1b% r_leo)

    if (associated(rodata% Lev1b% v_gns))      &
      call ropp_io_addvar(                     &
         rodata,                               &
         'v_gns_1b',                           &
         'GNSS transmitter velocity level 1b', &
         'meters',                             &
         (/ -10000._wp, 10000._wp/),           &
         rodata% Lev1b% v_gns)

    if (associated(rodata% Lev1b% v_leo))      &
      call ropp_io_addvar(                     &
         rodata,                               &
         'v_leo_1b',                           &
         'LEO transmitter velocity level 1b',  &
         'meters',                             &
         (/ -10000._wp, 10000._wp/),           &
         rodata% Lev1b% v_leo)

    !----------------------
    ! call original routine
    !----------------------
    call ropp_io_write_orig (ROdata, file, path, type, append, rec, ierr)
  end subroutine ropp_io_write
!------------------------------------------------------------------------------
  subroutine dt7_time (dt7, time) ! assignment: DT7type = t_time
  type (DT7type) ,intent(inout) :: dt7  ! ROPP  time data type
  type (t_time)  ,intent(in)    :: time ! 3dvar time data type
    dt7% year   = iyyyy (time)
    dt7% month  = imm   (time)
    dt7% day    = idd   (time)
    dt7% hour   = ihh   (time)
    dt7% minute = imi   (time)
    dt7% second = iss   (time)
    dt7% msec   = 0
  end subroutine dt7_time
!==============================================================================
  subroutine impact_heights (h, file)
  real(wp)         ,pointer               :: h (:)  ! returned impact heights
  character(len=*) ,intent(in)  ,optional :: file   ! file name
  !--------------------------
  ! read ROPIC impact heights
  !--------------------------
    integer            :: iunit, i, stat, ni, id
    character(len=256) :: fil
    real(wp)           :: rd

    iunit=get_unit_number()
    if (present (file)) then
      fil = file
    else
      fil = path_file (data, 'ROPIC_Impact_Heights.dat')
    endif
      open(iunit,file=fil,action='read',status='old')
      i=1
    do
      read (iunit,*,iostat=stat) id, rd
      if (stat<0) then
        exit
      else if (stat>0) then
        call finish('impact_heights','ERROR while reading '//trim(fil))
      endif
      i=i+1
    enddo
    ni = i-1

    allocate(h(ni))
    rewind(iunit)
    do i=1,ni
      read (iunit,*,iostat=stat) id, h(i)
    enddo
    close(iunit)
    call return_unit_number(iunit)
  end subroutine impact_heights
!==============================================================================
#endif
end module ropp_io_dwd
