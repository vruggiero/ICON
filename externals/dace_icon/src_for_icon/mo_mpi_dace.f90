!
!+ Fortran 90 interface to MPI library
!
MODULE mo_mpi_dace
!
! Description:
!   Fortran 90 interface to MPI library.
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
! V1_4         2009/03/26 Harald Anlauf
!  Implement p_gatherv,p_scatterv; Instrumentation for checksumming
! V1_5         2009/05/25 Harald Anlauf
!  new routines: p_min_i2_1d, p_min_i1_1d, p_gather_real
!  improve MPI checksumming
! V1_6         2009/06/10 Harald Anlauf
!  improve MPI checksumming
! V1_8         2009/12/09 Andreas Rhodin
!  new subroutine p_max_i2d
! V1_9         2010/04/20 Harald Anlauf
!  set public: MPI_ANY_SOURCE
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  p_bcast: implement p_bcast_real_0d_sp
! V1_13        2011/11/01 Andreas Rhodin
!  new generic and specific routines
! V1_15        2011/12/06 Harald Anlauf
!  Remove workarounds for NEC SX, AIX; new workaround for Vampir on MPI/SX
! V1_17        2011/12/21 Harald Anlauf
!  check_type_dp: initialize array
! V1_19        2012-04-16 Harald Anlauf
!  add_veri: use p_gatherv instead of p_send/p_recv
!  p_start: fix off-by-one error in selection of I/O processor
! V1_20        2012-06-18 Andreas Rhodin
!  new specific routine p_gather_real_sp
! V1_27        2013-11-08 Andreas Rhodin
!  new specific function p_max_1s: p_max (real(sp))
! V1_28        2014/02/26 Andreas Rhodin
!  implement p_allgatherv_int
! V1_29        2014/04/02 Andreas Rhodin
!  p_scatterv_real_1d: disable consistency check for the receiver
! V1_31        2014-08-21 Robin Faulwetter
!  new gatherv routines for new write_rttov_prof routine;
!  new specific p_send/p_recv
! V1_35        2014-11-07 Andreas Rhodin
!  new specific routines: p_min_int_1d, p_sum_i2_1d, p_sum_i1_1d, p_or_i4_1d
! V1_42        2015-06-08 Harald Anlauf
!  mo_mpi_dace::p_start: properly initialize hybrid MPI-OpenMP mode
! V1_45        2015-12-15 Harald Anlauf
!  Add non-blocking p_i{send,recv}_derivedtype; add generic p_ibcast
! V1_48        2016-10-06 Harald Anlauf
!  implement p_irecv_real_2d
! V1_51        2017-02-24 Andreas Rhodin
!  new subroutines p_min_atm,p_max_atm
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Luis Kornblueh  MPIfM     ????       original code
! Andreas Rhodin  MPIfM/DWD 2001-2008
! Harald Anlauf   DWD       2007-2009
!------------------------------------------------------------------------------
  USE mo_kind
  USE mo_fortran_units, ONLY: nerr ! standard error fortran unit
#ifdef NOMPI
  USE mo_system,        ONLY: abort
#endif
#ifdef _OPENMP
  USE omp_lib,          only: omp_get_max_threads,    &!
                              omp_get_max_active_levels
#endif
#ifndef NOMPI
#ifdef  HAVE_MPI_MOD    /* prefer MPI module over mpif.h */
  USE mpi
#endif
#endif

#define MO_MPI_SOURCE             ! Avoid recursive USE when including p_*.incf

!#if defined (__SX__)
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

  IMPLICIT NONE

  PRIVATE                   ! declarations are private by default
  !------------------------------------------------------------
  ! subroutines to initialise and shut down the MPI environment
  !------------------------------------------------------------
  PUBLIC :: p_start         ! initialse the parallel environment.
  PUBLIC :: p_stop          ! finish MPI and clean up all PEs
  PUBLIC :: p_abort         ! abort the program (with error code /= 0)
  PUBLIC :: set_dace_comm   ! set the DACE communicator and defaults
  !----------------------------------------------
  ! generic point to point communication routines
  !----------------------------------------------
  PUBLIC :: p_send          ! generic MPI_send      routine
  PUBLIC :: p_recv          ! generic MPI_recv      routine
  PUBLIC :: p_sendrecv      ! generic MPI_sendrecv  routine
  PUBLIC :: p_isend         ! generic MPI non-blocking send
  PUBLIC :: p_irecv         ! generic MPI non-blocking receive
  !--------------------------------
  ! requests for non blocking calls
  !--------------------------------
  PUBLIC :: p_wait          ! waits for an MPI send or receive to complete
  PUBLIC :: p_waitall       ! waits for all given MPI requests to complete
  PUBLIC :: p_probe         ! check incoming messages w/o actually receiving
  PUBLIC :: p_request       ! Array of handles of non-blocking requests
  PUBLIC :: p_irequest      ! Index to next request handle
  PUBLIC :: p_mrequest      ! Maximum of outstanding non-blocking requests
  !------------------------------------------
  ! generic collective communication routines
  !------------------------------------------
  PUBLIC :: p_barrier       !         MPI barrier   routine
  PUBLIC :: p_bcast         ! generic MPI broadcast routine
  PUBLIC :: p_bcast_ptr     ! generic MPI broadcast routine (pointer + content)
  PUBLIC :: p_ibcast        ! generic MPI broadcast routine (non-blocking)
  PUBLIC :: p_alltoall      ! generic MPI_alltoall  routine
  PUBLIC :: p_gather        ! generic MPI_gather    routine
  PUBLIC :: p_allgather     ! generic MPI_allgather routine
  PUBLIC :: p_scatterv      ! generic MPI_scatterv  routine
  PUBLIC :: p_gatherv       ! generic MPI_gatherv   routine
  PUBLIC :: alltoallv_args  ! set up/check arguments for MPI alltoallv
  PUBLIC :: allgatherv_args ! set up/check arguments for MPI allgatherv
  PUBLIC :: gatherv_args    ! set up/check arguments for MPI gatherv
  !------------------------------------------------------
  ! generic MPI reduction routines (minimum, maximum, ..)
  !------------------------------------------------------
  PUBLIC :: p_max           ! calculate maximum value of variable on all PEs
  PUBLIC :: p_min           ! calculate minimum value of variable on all PEs
  PUBLIC :: p_sum           ! calculate sum of values of variable on all PEs
  PUBLIC :: p_or            ! calculate logical .or.  of variable on all PEs
  PUBLIC :: p_and           ! calculate logical .and. of variable on all PEs
  PUBLIC :: p_ior           ! calculate bitwise .or.  of variable on all PEs
  !------------------------------
  ! processor element identifiers !+++ disabled, use dace%.. instead +++
  !------------------------------
! PUBLIC :: p_pe            ! Id of THIS pe (starting with 0)
! PUBLIC :: p_io            ! Id of processor which should do I/O
! PUBLIC :: p_nprocs        ! number of available processor elements
! PUBLIC :: npes
  !--------------
  ! logical flags
  !--------------
  PUBLIC :: p_parallel      ! true for parallel environment
  PUBLIC :: p_parallel_io   ! true if this PE is p_io
  !---------------------------------------
  ! pre-defined MPI communicator group Ids !+++ disabled use dace%.. world%.. self%..
  !---------------------------------------
! PUBLIC :: MPI_COMM_WORLD    ! group consisting of all PEs
! PUBLIC :: MPI_COMM_SELF     ! group consisting of this PE only
! PUBLIC :: MPI_COMM_NULL     ! group consisting of no PE
! PUBLIC :: p_communicator_a  !
! PUBLIC :: p_communicator_b
! PUBLIC :: p_communicator_d
  !------------------------------------
  ! pre-defined MPI communicator groups
  !------------------------------------
  PUBLIC :: world   ! MPI_COMM_WORLD group info
  PUBLIC :: self    ! MPI_COMM_SELF  group info
  PUBLIC :: dace    ! DACE           group info
  PUBLIC :: d_npe   ! copy of dace% npe
  PUBLIC :: d_comm  ! copy of dace% comm
  PUBLIC :: d_pe    ! copy of dace% pe
  !--------------------
  ! MPI type parameters
  !--------------------
  PUBLIC :: p_real          ! default real
  PUBLIC :: p_real_wp       ! real(wp)
  PUBLIC :: p_real_dp       ! real(dp)
  PUBLIC :: p_type          ! return MPI_type for real(sp,dp), integer(i4,i8)
#ifndef NOMPI               /* Crucial for inlining: */
  PUBLIC :: p_real_sp       ! real(sp)
  PUBLIC :: p_int           ! default integer
  PUBLIC :: p_int_i2
  PUBLIC :: p_int_i4
  PUBLIC :: p_int_i8
  PUBLIC :: p_bool
  PUBLIC :: p_char
#endif
  ! native types
  PUBLIC :: MPI_BYTE        !
  PUBLIC :: MPI_INTEGER     ! default integer
  PUBLIC :: MPI_INTEGER2    ! integer(*2)
  PUBLIC :: MPI_CHARACTER
  PUBLIC :: MPI_DOUBLE_PRECISION
  !--------------------
  ! other MPI constants
  !--------------------
  PUBLIC :: MPI_ANY_SOURCE
  PUBLIC :: MPI_ANY_TAG
  PUBLIC :: MPI_STATUS_SIZE
  PUBLIC :: MPI_SUCCESS
  PUBLIC :: MPI_UNDEFINED
  PUBLIC :: nbcast
  !--------------------------
  ! diagnostics and debugging
  !--------------------------
  PUBLIC :: crc      ! Checksum for a bitstream (integers, logicals, reals)
#ifdef DEBUG
  public :: mype     ! PE number of this task
#endif
  !------------------------------------------------------
  ! Enable checksums for selected communication routines:
  !   p_{send,recv,bcast}_derivedtype, p_bcast_DERIVED
!#define MPI_CHECKSUM
  !------------------------------------------------------
#ifdef MPI_CHECKSUM
  public :: dump_buffer_int
#endif

  ! old fashioned method (MPI-1)

#ifndef NOMPI

#ifndef HAVE_MPI_MOD    /* beware: prefer MPI module over mpif.h */
INCLUDE 'mpif.h'
#endif

#else

  ! define some dummy parameters used by derived-type send/recv routines

  INTEGER, PARAMETER :: MPI_ANY_SOURCE       = 0
  INTEGER, PARAMETER :: MPI_ANY_TAG          = 0
  INTEGER, PARAMETER :: MPI_COMM_WORLD       = 0
  INTEGER, PARAMETER :: MPI_COMM_SELF        = 0
  INTEGER, PARAMETER :: MPI_COMM_NULL        = 0
  INTEGER, PARAMETER :: MPI_BYTE             = 0
  INTEGER, PARAMETER :: MPI_STATUS_SIZE      = 0
  INTEGER, PARAMETER :: MPI_SUCCESS          = 0
  INTEGER, PARAMETER :: MPI_INTEGER          = 0
  INTEGER, PARAMETER :: MPI_INTEGER2         = 0
  INTEGER, PARAMETER :: MPI_DOUBLE_PRECISION = 0
  INTEGER, PARAMETER :: MPI_CHARACTER        = 0
  INTEGER, PARAMETER :: MPI_UNDEFINED        = 0
  INTEGER, PARAMETER :: p_real               = 0
  INTEGER, PARAMETER :: p_real_wp            = 0
  INTEGER, PARAMETER :: p_real_dp            = 0
#endif

  ! IBM RS/6000 SP MPI has not defined the optional MPI_INTEGER8.
  ! It is redefined to MPI_INTEGER.  (This is now obsolete?)

! #ifdef _AIX
!#ifdef AIX
!#define MPI_INTEGER8 MPI_INTEGER
!#endif

  ! general run time information

#ifndef NOMPI
  INTEGER :: version, subversion   ! MPI version
#endif

  ! MPI call inherent variables

#ifndef NOMPI
  INTEGER :: p_error                     ! MPI error number

  INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV
#endif

  ! public parallel run information (defaults provided in case mpi is not used)

  LOGICAL :: p_parallel    = .FALSE. ! flag for parallel run
  LOGICAL :: p_parallel_io = .TRUE.  ! IO is done on this PE
  INTEGER :: p_pe          = 0       ! processor index of this task
  INTEGER :: p_io          = 0       ! processor index of PE handling IO
  INTEGER :: p_nprocs      = 1       ! number of available PEs (processors)

  ! communicator sets

! INTEGER :: p_communicator_a ! for Set A
! INTEGER :: p_communicator_b ! for Set B
! INTEGER :: p_communicator_d ! for debug node

  ! non blocking calls
  ! the first p_irequest-1 values of p_request are in use

  INTEGER, ALLOCATABLE :: p_request(:)
  INTEGER :: p_irequest = 0   ! Index to next request handle
  INTEGER :: p_mrequest = 0   ! Maximum of outstanding non-blocking requests

  ! module intrinsic names

  INTEGER :: mype                  ! this is the PE number of this task
#ifndef NOMPI
  INTEGER :: iope                  ! PE able to do IO
#endif
  INTEGER :: npes = 1              ! number of available PEs

  INTEGER :: nbcast = 0            ! counter for broadcasts for debugging

  ! MPI transfer types

#ifndef NOMPI
  INTEGER :: p_real_sp  = MPI_REAL
  INTEGER :: p_real_dp  = MPI_DOUBLE_PRECISION
  INTEGER :: p_real_wp  = MPI_DOUBLE_PRECISION  ! DACE default
  INTEGER :: p_int_i2   = MPI_INTEGER2
  INTEGER :: p_int_i4   = MPI_INTEGER4
  INTEGER :: p_int_i8   = MPI_INTEGER8

  ! native types

  INTEGER :: p_int      = MPI_INTEGER  ! may depend on compiler options
  INTEGER :: p_real     = MPI_REAL     ! may depend on compiler options

  INTEGER :: p_bool     = MPI_LOGICAL
  INTEGER :: p_char     = MPI_CHARACTER

  ! temporary types for probing

  INTEGER :: p_i4, p_i8
  INTEGER :: p_sp, p_dp
#endif

  ! for checking out integer and real variables separately.
  ! KIND values usually overlap and give the byte size or are defined
  ! as sequence separate for both groups.

  INTEGER, PARAMETER :: real_type    = 1
  INTEGER, PARAMETER :: integer_type = 2

  ! define generic interfaces to allow proper compiling
  ! with picky compilers like NAG f95 for clean argument checking and
  ! shortening the call sequence.

  INTERFACE p_send
  !-------------------------------------------------
  ! generic interface for specific MPI send routines
  !-------------------------------------------------
    MODULE PROCEDURE p_send_real       ! send real(dp)
    MODULE PROCEDURE p_send_real_1d    ! send real(dp)        1d array
    MODULE PROCEDURE p_send_real_2d    ! send real(dp)        2d array
    MODULE PROCEDURE p_send_real_3d    ! send real(dp)        3d array
    MODULE PROCEDURE p_send_real_4d    ! send real(dp)        4d array
    MODULE PROCEDURE p_send_real_5d    ! send real(dp)        5d array
    MODULE PROCEDURE p_send_real_1d_sp ! send real(sp)        1d array
    MODULE PROCEDURE p_send_real_2d_sp ! send real(sp)        2d array
    MODULE PROCEDURE p_send_int        ! send default integer
    MODULE PROCEDURE p_send_int_1d     ! send integer         1d array
    MODULE PROCEDURE p_send_int_2d     ! send integer         2d array
    MODULE PROCEDURE p_send_int_3d     ! send integer         3d array
    MODULE PROCEDURE p_send_int_4d     ! send integer         4d array
    MODULE PROCEDURE p_send_int_1d_i2  ! send integer(i2)
#ifdef HAVE_I1
    MODULE PROCEDURE p_send_int_1d_i1  ! send integer(i1)
#endif
    MODULE PROCEDURE p_send_bool       ! send boolean
    MODULE PROCEDURE p_send_bool_1d    ! send boolean         1d array
    MODULE PROCEDURE p_send_bool_2d    ! send boolean         2d array
    MODULE PROCEDURE p_send_bool_3d    ! send boolean         3d array
    MODULE PROCEDURE p_send_bool_4d    ! send boolean         4d array
    MODULE PROCEDURE p_send_char       ! send character
    MODULE PROCEDURE p_send_char_1d    ! send character          array
  END INTERFACE

  INTERFACE p_isend
  !--------------------------------------------------
  ! generic interface for unblocked MPI send routines
  !--------------------------------------------------
    MODULE PROCEDURE p_isend_real      ! send real(dp)
    MODULE PROCEDURE p_isend_real_1d   ! send real(dp)        1d array
    MODULE PROCEDURE p_isend_real_2d   ! send real(dp)        2d array
    MODULE PROCEDURE p_isend_real_3d   ! send real(dp)        3d array
    MODULE PROCEDURE p_isend_real_4d   ! send real(dp)        4d array
    MODULE PROCEDURE p_isend_real_5d   ! send real(dp)        5d array
    MODULE PROCEDURE p_isend_int       ! send default integer
    MODULE PROCEDURE p_isend_int_1d    ! send integer         1d array
    MODULE PROCEDURE p_isend_int_2d    ! send integer         2d array
    MODULE PROCEDURE p_isend_int_3d    ! send integer         3d array
    MODULE PROCEDURE p_isend_int_4d    ! send integer         4d array
    MODULE PROCEDURE p_isend_bool      ! send boolean
    MODULE PROCEDURE p_isend_bool_1d   ! send boolean         1d array
    MODULE PROCEDURE p_isend_bool_2d   ! send boolean         2d array
    MODULE PROCEDURE p_isend_bool_3d   ! send boolean         3d array
    MODULE PROCEDURE p_isend_bool_4d   ! send boolean         4d array
    MODULE PROCEDURE p_isend_char      ! send character
    MODULE PROCEDURE p_isend_char_1d   ! send character          array
  END INTERFACE

  INTERFACE p_recv
  !----------------------------------------------------
  ! generic interface for specific MPI receive routines
  !----------------------------------------------------
    MODULE PROCEDURE p_recv_real       ! receive real(dp)
    MODULE PROCEDURE p_recv_real_1d    ! receive real(dp)        1d array
    MODULE PROCEDURE p_recv_real_2d    ! receive real(dp)        2d array
    MODULE PROCEDURE p_recv_real_3d    ! receive real(dp)        3d array
    MODULE PROCEDURE p_recv_real_4d    ! receive real(dp)        4d array
    MODULE PROCEDURE p_recv_real_5d    ! receive real(dp)        5d array
    MODULE PROCEDURE p_recv_real_1d_sp ! receive real(sp)        1d array
    MODULE PROCEDURE p_recv_real_2d_sp ! receive real(sp)        2d array
    MODULE PROCEDURE p_recv_int        ! receive default integer
    MODULE PROCEDURE p_recv_int_1d     ! receive integer         1d array
    MODULE PROCEDURE p_recv_int_2d     ! receive integer         2d array
    MODULE PROCEDURE p_recv_int_3d     ! receive integer         3d array
    MODULE PROCEDURE p_recv_int_4d     ! receive integer         4d array
    MODULE PROCEDURE p_recv_int_1d_i2  ! receive integer(i2)
#ifdef HAVE_I1
    MODULE PROCEDURE p_recv_int_1d_i1  ! receive integer(i1)
#endif
    MODULE PROCEDURE p_recv_bool       ! receive boolean
    MODULE PROCEDURE p_recv_bool_1d    ! receive boolean         1d array
    MODULE PROCEDURE p_recv_bool_2d    ! receive boolean         2d array
    MODULE PROCEDURE p_recv_bool_3d    ! receive boolean         3d array
    MODULE PROCEDURE p_recv_bool_4d    ! receive boolean         4d array
    MODULE PROCEDURE p_recv_char       ! receive character
    MODULE PROCEDURE p_recv_char_1d    ! receive character          array
  END INTERFACE

  INTERFACE p_irecv
  !-----------------------------------------------------
  ! generic interface for unblocked MPI receive routines
  !-----------------------------------------------------
    MODULE PROCEDURE p_irecv_real      ! receive real(dp)
    MODULE PROCEDURE p_irecv_real_1d   ! receive real(dp)        1d array
    MODULE PROCEDURE p_irecv_real_2d   ! receive real(dp)        2d array
!   MODULE PROCEDURE p_irecv_real_3d
!   MODULE PROCEDURE p_irecv_real_4d
    MODULE PROCEDURE p_irecv_int       ! receive default integer
  END INTERFACE

  INTERFACE p_sendrecv
  !---------------------------------------------------------
  ! generic interface for specific MPI send/receive routines
  !---------------------------------------------------------
     MODULE PROCEDURE p_sendrecv_real_2d ! send/receive real(dp) 2d array
     MODULE PROCEDURE p_sendrecv_real_3d ! send/receive real(dp) 3d array
     MODULE PROCEDURE p_sendrecv_real_4d ! send/receive real(dp) 4d array
  END INTERFACE

  INTERFACE p_bcast
  !------------------------------
  ! generic MPI broadcast routine
  !------------------------------
    MODULE PROCEDURE p_bcast_int_i4      ! broadcast integer(i4)
    MODULE PROCEDURE p_bcast_int_1d      ! broadcast integer     1d array
    MODULE PROCEDURE p_bcast_int_2d      ! broadcast integer     2d array
    MODULE PROCEDURE p_bcast_int_3d      ! broadcast integer     3d array
    MODULE PROCEDURE p_bcast_int_4d      ! broadcast integer     4d array
    MODULE PROCEDURE p_bcast_int_i8      ! broadcast integer(i8)
    MODULE PROCEDURE p_bcast_int_1d_i2   ! broadcast integer(i2) 1d array
#ifdef HAVE_I1
    MODULE PROCEDURE p_bcast_int_1d_i1   ! broadcast integer(i1) 1d array
    MODULE PROCEDURE p_bcast_int_2d_i1   ! broadcast integer(i1) 2d array
#endif
    MODULE PROCEDURE p_bcast_real        ! broadcast real(dp)
    MODULE PROCEDURE p_bcast_real_1d     ! broadcast real(dp)    1d array
    MODULE PROCEDURE p_bcast_real_2d     ! broadcast real(dp)    2d array
    MODULE PROCEDURE p_bcast_real_3d     ! broadcast real(dp)    3d array
    MODULE PROCEDURE p_bcast_real_4d     ! broadcast real(dp)    4d array
    MODULE PROCEDURE p_bcast_real_0d_sp  ! broadcast real(sp)
    MODULE PROCEDURE p_bcast_real_1d_sp  ! broadcast real(sp)    1d array
    MODULE PROCEDURE p_bcast_real_2d_sp  ! broadcast real(sp)    2d array
    MODULE PROCEDURE p_bcast_real_3d_sp  ! broadcast real(sp)    3d array
    MODULE PROCEDURE p_bcast_bool        ! broadcast boolean
    MODULE PROCEDURE p_bcast_bool_1d     ! broadcast boolean     1d array
    MODULE PROCEDURE p_bcast_bool_2d     ! broadcast boolean     2d array
    MODULE PROCEDURE p_bcast_bool_3d     ! broadcast boolean     3d array
    MODULE PROCEDURE p_bcast_bool_4d     ! broadcast boolean     4d array
    MODULE PROCEDURE p_bcast_char        ! broadcast character
    MODULE PROCEDURE p_bcast_char_1d     ! broadcast character   1d array
  END INTERFACE

  INTERFACE p_bcast_ptr
    !--------------------------------------------------------
    ! generic MPI broadcast routine (pointer + content)
    ! allocate the receive pointers with the appropriate size
    !--------------------------------------------------------
    MODULE PROCEDURE p_bcast_ptr_int_1d  ! broadcast integer     1d pointer
    MODULE PROCEDURE p_bcast_ptr_int_2d  ! broadcast integer     2d pointer
#ifdef HAVE_I1
    MODULE PROCEDURE p_bcast_ptr_int1_1d ! broadcast integer(i1) 1d pointer
    MODULE PROCEDURE p_bcast_ptr_int1_2d ! broadcast integer(i1) 2d pointer
#endif
    MODULE PROCEDURE p_bcast_ptr_real_1d ! broadcast real(dp)    1d pointer
    MODULE PROCEDURE p_bcast_ptr_real_2d ! broadcast real(dp)    2d pointer
    MODULE PROCEDURE p_bcast_ptr_sp_1d   ! broadcast real(sp)    1d pointer
    MODULE PROCEDURE p_bcast_ptr_sp_2d   ! broadcast real(sp)    2d pointer
  END INTERFACE p_bcast_ptr

  INTERFACE p_ibcast
  !-------------------------------------------
  ! generic non-blocking MPI broadcast routine
  !-------------------------------------------
    MODULE PROCEDURE p_ibcast_int_1d     ! broadcast integer     1d array
    MODULE PROCEDURE p_ibcast_int_2d     ! broadcast integer     2d array
    MODULE PROCEDURE p_ibcast_int_3d     ! broadcast integer     3d array
    MODULE PROCEDURE p_ibcast_int_4d     ! broadcast integer     4d array
    MODULE PROCEDURE p_ibcast_real_1d    ! broadcast real(dp)    1d array
    MODULE PROCEDURE p_ibcast_real_2d    ! broadcast real(dp)    2d array
    MODULE PROCEDURE p_ibcast_real_3d    ! broadcast real(dp)    3d array
    MODULE PROCEDURE p_ibcast_real_4d    ! broadcast real(dp)    4d array
  END INTERFACE p_ibcast

  INTERFACE p_probe
  !--------------------------------------------------------
  ! check incoming messages without actually receiving them
  !--------------------------------------------------------
    MODULE PROCEDURE p_probe_real
    MODULE PROCEDURE p_probe_int
    MODULE PROCEDURE p_probe_bool
    MODULE PROCEDURE p_probe_char
  END INTERFACE

  INTERFACE p_max
  !-----------------------------------------------
  ! generic MPI reduction routine:
  ! calculate maximum value of variable on all PEs
  !-----------------------------------------------
    MODULE PROCEDURE p_max_i0d     ! maximum of    integer
    MODULE PROCEDURE p_max_i1d     ! maximum of 1d integer array
    MODULE PROCEDURE p_max_i2d     ! maximum of 2d integer array
    MODULE PROCEDURE p_max_int_i8  ! maximum of    integer(i8)
    MODULE PROCEDURE p_max_0d      ! maximum of    real
    MODULE PROCEDURE p_max_1d      ! maximum of 1d real    array
    MODULE PROCEDURE p_max_1s      ! maximum of 1d real(sp)array
    MODULE PROCEDURE p_max_2d      ! maximum of 2d real    array
    MODULE PROCEDURE p_max_3d      ! maximum of 3d real    array
    MODULE PROCEDURE p_max_4d      ! maximum of 4d real    array
  END INTERFACE

  INTERFACE p_min
  !-----------------------------------------------
  ! generic MPI reduction routine:
  ! calculate minimum value of variable on all PEs
  !-----------------------------------------------
    MODULE PROCEDURE p_min_0d      ! minimum of    real
    MODULE PROCEDURE p_min_1d      ! minimum of 1d real        array
    MODULE PROCEDURE p_min_2d      ! minimum of 2d real        array
    MODULE PROCEDURE p_min_3d      ! minimum of 3d real        array
    MODULE PROCEDURE p_min_4d      ! minimum of 4d real        array
    MODULE PROCEDURE p_min_int_0d  ! minimum of    integer
    MODULE PROCEDURE p_min_int_1d  ! minimum of    integer(:)
    MODULE PROCEDURE p_min_int_i8  ! minimum of    integer(i8)
    MODULE PROCEDURE p_min_i2_1d   ! minimum of 1d integer(i2) array
#ifdef HAVE_I1
    MODULE PROCEDURE p_min_i1_1d   ! minimum of 1d integer(i1) array
#endif
  END INTERFACE

  INTERFACE p_sum
  !-------------------------------------
  ! generic MPI reduction routine:
  ! calculate sum of variable on all PEs
  !-------------------------------------
    MODULE PROCEDURE p_sum_0d      ! sum of    real
    MODULE PROCEDURE p_sum_1d      ! sum of 1d real    array
    MODULE PROCEDURE p_sum_2d      ! sum of 2d real    array
    MODULE PROCEDURE p_sum_3d      ! sum of 3d real    array
    MODULE PROCEDURE p_sum_4d      ! sum of 4d real    array
    MODULE PROCEDURE p_sum_int_0d  ! sum of    integer
    MODULE PROCEDURE p_sum_int_1d  ! sum of 1d integer array
    MODULE PROCEDURE p_sum_int_2d  ! sum of 2d integer array
    MODULE PROCEDURE p_sum_int_3d  ! sum of 3d integer array
    MODULE PROCEDURE p_sum_int_4d  ! sum of 4d integer array
    MODULE PROCEDURE p_sum_int_i8  ! sum of    integer(i8)
    MODULE PROCEDURE p_sum_i8_1d   ! sum of 1d integer(i8)
    MODULE PROCEDURE p_sum_i2_1d   ! sum of 1d integer(i2) array
#ifdef HAVE_I1
    MODULE PROCEDURE p_sum_i1_1d   ! sum of 1d integer(i1) array
#endif
  END INTERFACE

  INTERFACE p_or
  !----------------------------------------------
  ! generic MPI reduction routine:
  ! calculate logical .or. of variable on all PEs
  !----------------------------------------------
    MODULE PROCEDURE p_or_0d       ! logical .or.
    MODULE PROCEDURE p_or_1d       ! logical .or. of 1d boolean array
    MODULE PROCEDURE p_or_i4_1d    ! minimum of 1d integer(i2) array
  END INTERFACE

  INTERFACE p_and
  !-----------------------------------------------
  ! generic MPI reduction routine:
  ! calculate logical .and. of variable on all PEs
  !-----------------------------------------------
    MODULE PROCEDURE p_and_0d      ! logical .and.
    MODULE PROCEDURE p_and_1d      ! logical .and. of 1d boolean array
  END INTERFACE

  INTERFACE p_ior
  !----------------------------------------------
  ! generic MPI reduction routine:
  ! calculate bitwise .or. of variable on all PEs
  !----------------------------------------------
    MODULE PROCEDURE p_ior_0d_i4       ! bitwise .or.
    MODULE PROCEDURE p_ior_0d_i8       ! bitwise .or.
    MODULE PROCEDURE p_ior_1d_i4       ! bitwise .or.
    MODULE PROCEDURE p_ior_1d_i8       ! bitwise .or.
  END INTERFACE

  INTERFACE p_gather
  !---------------------------------
  ! one PE gathers data from all PEs
  !---------------------------------
    MODULE PROCEDURE p_gather_bool    ! gather boolean  data
    MODULE PROCEDURE p_gather_int     ! gather integer  data
    MODULE PROCEDURE p_gather_real_sp ! gather real(sp) data
    MODULE PROCEDURE p_gather_real    ! gather real(dp) data
  END INTERFACE

  INTERFACE p_allgather
  !---------------------------------
  ! all PEs gather data from all PEs
  !---------------------------------
    MODULE PROCEDURE p_allgather_int   ! gather integer data
    MODULE PROCEDURE p_allgatherv_int  ! gather integer data
  END INTERFACE

  INTERFACE p_alltoall
  !-----------------------------------------
  ! generic MPI_alltoall  routine
  ! all PEs send and receive to/from all PEs
  !-----------------------------------------
    MODULE PROCEDURE p_alltoall_int   ! send/receive integers
    MODULE PROCEDURE p_alltoall_real  ! send/receive reals
    MODULE PROCEDURE p_alltoall_bool  ! send/receive booleans
  END INTERFACE

  INTERFACE p_scatterv
  !---------------------------------------
  ! one PE scatters vector data to all PEs
  !---------------------------------------
    MODULE PROCEDURE p_scatterv_real_1d  ! scatter real vectors
    MODULE PROCEDURE p_scatterv_int_1d   ! scatter int vectors
  END INTERFACE

  INTERFACE p_gatherv
  !----------------------------------------
  ! one PE gathers vector data from all PEs
  !----------------------------------------
    MODULE PROCEDURE p_gatherv_int_1d_i4  ! gather integer(i4) vectors
    MODULE PROCEDURE p_gatherv_int_2d_i4  ! gather integer(i4) vectors
    MODULE PROCEDURE p_gatherv_real_1d    ! gather real(dp)    vectors
    MODULE PROCEDURE p_gatherv_real_1d_sp ! gather real(sp)    vectors
    MODULE PROCEDURE p_gatherv_real_2d    ! gather real(dp)    vectors
    MODULE PROCEDURE p_gatherv_real_3d    ! gather real(dp)    vectors
    MODULE PROCEDURE p_gatherv_char_1d    ! gather character   vectors
    MODULE PROCEDURE p_gatherv_bool_2d    ! gather logical     vectors
  END INTERFACE

  INTERFACE crc
  !-----------------------------------------------------
  ! Checksum for a bitstream (integers, logicals, reals)
  !-----------------------------------------------------
    MODULE PROCEDURE crc32_int      ! Checksum for array of integers
    MODULE PROCEDURE crc32_int_i2   ! Checksum for array of short integers
    MODULE PROCEDURE crc32_real     ! Checksum for array of reals
    MODULE PROCEDURE crc32_real_sp  ! Checksum for array of reals
    MODULE PROCEDURE crc32_logical  ! Checksum for array of logicals
    MODULE PROCEDURE crc32_char     ! Checksum for array of characters
  END INTERFACE

  INTERFACE p_type
  !------------------------------------------------
  ! return MPI_type for real(sp,dp), integer(i4,i8)
  !------------------------------------------------
    MODULE PROCEDURE p_type_i4
    MODULE PROCEDURE p_type_i8
    MODULE PROCEDURE p_type_sp
    MODULE PROCEDURE p_type_dp
  END INTERFACE

  !----------------------------------------
  ! derived type to hold commumication info
  ! defaults apply to single processor run
  !----------------------------------------
  type t_comm
    integer :: comm = 0      ! communicator group Id
    integer :: pe   = 0      ! Id of THIS pe (starting with 0)
    integer :: pio  = 0      ! Id of processor which should do I/O
    logical :: lpio = .true. ! .true. if this PE should do I/O
    integer :: npe  = 1      ! number of available processor elements
  end type t_comm

  !--------------------------------------------------
  ! derived type variables to hold communication info
  !--------------------------------------------------
  type (t_comm) :: world      ! MPI_COMM_WORLD
  type (t_comm) :: self       ! MPI_COMM_SELF
  type (t_comm) :: dace       ! communicator for DACE, usually == world
  integer       :: d_comm = 0 ! copy of dace% comm
  integer       :: d_npe  = 1 ! copy of dace% npe
  integer       :: d_pe   = 0 ! copy of dace% pe
  !--------------------------------------------------------
  ! prevent communication info from accidental modification
  !--------------------------------------------------------
  protected     :: world, self, dace, d_comm, d_npe, d_pe
  save          :: world, self, dace, d_comm, d_npe, d_pe

CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE p_start
  !----------------------------------------------------------------------
  ! Initialise the parallel environment.
  ! Should be called at the beginning of the program.
  !
  ! Variables are provided for determing I/O size in bytes of the defined
  ! KIND types for assigning the right MPI data types with the used kinds
  !----------------------------------------------------------------------

#ifndef NOMPI
    INTEGER :: io_size, integer_io_size, integer_byte_size

    INTEGER      :: iig = 0
    INTEGER (i4) :: ii4 = 0_i4
    INTEGER (i8) :: ii8 = 0_i8

    REAL         :: rrg = 0.0
    REAL (sp)    :: rsp = 0.0_sp
    REAL (dp)    :: rdp = 0.0_dp

    INTEGER      :: p_ig, p_rg

    ! temporary array to distibute the determined MPI types

    INTEGER :: p_send(8)

    ! variables used for determing the I/O PE

    integer :: i
    LOGICAL :: liope, linit
    INTEGER, ALLOCATABLE :: iope_table(:)
    LOGICAL, ALLOCATABLE :: liope_table(:)
#endif

    CHARACTER (132) :: io_pe_message
!$  integer :: provided

    ! temporaries to create DACE group

#ifndef NOMPI
    integer :: group_world, group_dace
#endif

    ! Executable statements:

    nbcast = 0

    io_pe_message(:) = ' '

    ! start MPI

#ifndef NOMPI
    CALL MPI_INITIALIZED (linit, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_INITIALIZED failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF

    IF (.NOT.linit) THEN
#ifdef _OPENMP
      ! We guarantee that only the master thread will call MPI routines
      CALL MPI_INIT_THREAD (MPI_THREAD_FUNNELED, provided, p_error)
#else
      CALL MPI_INIT (p_error)
#endif

      IF (p_error /= MPI_SUCCESS) THEN
         WRITE (nerr,'(a)') ' MPI_INIT failed.'
         WRITE (nerr,'(a,i4)') ' Error =  ', p_error
         STOP
      END IF

#ifdef _OPENMP
      ! Check if MPI_INIT_THREAD returned at least MPI_THREAD_FUNNELED in "provided"
      IF (provided < MPI_THREAD_FUNNELED) THEN
         WRITE (nerr,'(a)') &
              ' MPI_INIT_THREAD did not return desired level of thread support:'
         WRITE (nerr,'(a,i0)') " provided: ", provided
         WRITE (nerr,'(a,i0)') " required: ", MPI_THREAD_FUNNELED
         ! CALL MPI_Finalize(p_error)
         ! STOP
      END IF
#endif

    ENDIF
#endif

    ! get local PE identification

#ifndef NOMPI
    CALL MPI_COMM_RANK (MPI_COMM_WORLD, mype, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_RANK failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ELSE
#ifdef DEBUG
       WRITE (nerr,'(a,i4,a)') ' PE ', mype, ' started.'
#endif
    END IF
#else
    mype = 0
#endif

    ! get number of available PEs

#ifndef NOMPI
    CALL MPI_COMM_SIZE (MPI_COMM_WORLD, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error

       CALL p_abort
    END IF
#else
    npes = 1
#endif

#ifndef NOMPI
    ! for non blocking calls

    p_mrequest = 100*npes
    ALLOCATE (p_request(p_mrequest))
    p_irequest = 1

    ! look for a dedicated IO PE

    CALL MPI_ATTR_GET (MPI_COMM_WORLD, MPI_IO, iope, liope, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_ATTR_GET failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i0,2l2,1x,i0)') ' PE: ', mype,      &
         ' MPI_ATTR_GET returns iope, liope, any_source = ', &
         iope, liope, iope == MPI_ANY_SOURCE, MPI_ANY_SOURCE
#endif

    IF (iope == MPI_ANY_SOURCE) THEN

       ! all nodes can do IO

       IF (mype == 0) THEN
          WRITE (io_pe_message,'(a)') &
               '  All nodes can do I/O, selected PE 0 for I/O handling.'
       END IF
       p_io = 0

    ELSE

       ALLOCATE (iope_table(0:npes-1))

       IF (liope) THEN
          iope_table(mype) = iope
          CALL MPI_GATHER (iope      , 1, MPI_INTEGER, &
                           iope_table, 1, MPI_INTEGER, &
                           0, MPI_COMM_WORLD, p_error)

          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GATHER failed.'
             WRITE (nerr,'(a,i4)') ' Error =  ', p_error
             CALL p_abort
          END IF

          ! Now select the first given iope from table as IO PE.
          p_io = iope_table(0)
          CALL MPI_BCAST (p_io, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, p_error)
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_BCAST failed.'
             WRITE (nerr,'(a,i4)') ' Error =  ', p_error
             CALL p_abort
          END IF

          ! Workaround for a Vampir problem with MPI/SX
          if (p_io < 0) then
             ALLOCATE (liope_table(0:npes-1))
             liope_table       = .true.
             liope_table(mype) = liope
             CALL MPI_GATHER (liope      , 1, MPI_LOGICAL, &
                              liope_table, 1, MPI_LOGICAL, &
                              0, MPI_COMM_WORLD, p_error)

             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GATHER failed.'
                WRITE (nerr,'(a,i4)') ' Error =  ', p_error
                CALL p_abort
             END IF
             do i = 0, npes-1
                if (liope_table(i)) then
                   p_io = i
                   exit
                end if
             end do
             CALL MPI_BCAST (p_io, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, p_error)
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_BCAST failed.'
                WRITE (nerr,'(a,i4)') ' Error =  ', p_error
                CALL p_abort
             END IF
          end if

          IF (mype == 0) THEN
             WRITE (io_pe_message,'(a,i3,a)') &
                  '  Selected PE ', p_io, ' for I/O handling.'
          END IF

       ELSE
          ! if no dedicated IO PE is given, use PE 0
          p_io = 0

          WRITE (io_pe_message,'(a)') &
               '  No dedicated I/O PE, selected PE 0 for I/O handling.'
       END IF

       DEALLOCATE (iope_table)

    END IF
#else
    p_io = 0
#endif

    ! Information ...

    IF (mype == 0) THEN
       WRITE (nerr,'(/,a)') ' MPI interface runtime information:'
    END IF

     IF (npes < 2) THEN
       p_parallel = .FALSE.
       p_parallel_io = .TRUE.   ! can always do I/O
       IF (mype == 0) THEN
          WRITE (nerr,'(/)')
          WRITE (nerr,'(a)') '  Single processor run.'
       END IF
       p_pe = 0
       p_nprocs = 1
    ELSE
       p_parallel = .TRUE.
       IF (mype == p_io) THEN
          p_parallel_io = .TRUE.
       ELSE
          p_parallel_io = .FALSE.
       END IF
       IF (mype == 0) THEN
          WRITE (nerr,'(/)')
          WRITE (nerr,'(a,i4,a)') '  Run on ', npes, ' processors.'
       END IF
       p_pe = mype
       p_nprocs = npes
    END IF

#ifdef _OPENMP
    if (mype == 0) then
       if (p_parallel) then
          WRITE (nerr,'(a)') '  Running in hybrid MPI-OpenMP mode.'
       else
          WRITE (nerr,'(a)') '  Running in OpenMP mode.'
       end if
       WRITE (nerr,'(a,i0)')  &
            '  Maximum number of threads: ', omp_get_max_threads()
    end if
    if (omp_get_max_active_levels() /= 1) then
       if (mype == 0) then
          WRITE (nerr,'(a)')  '  Unsupported: OMP_MAX_ACTIVE_LEVELS /= 1'
       end if
    end if
#endif

    ! inform on I/O PE situation

    IF (mype == 0) THEN
       WRITE (nerr, '(a)') io_pe_message(1:LEN_TRIM(io_pe_message))
    END IF

#ifndef NOMPI

    IF (p_parallel) THEN

       ! let's check the available MPI version

       CALL MPI_GET_VERSION (version, subversion, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_GET_VERSION failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (mype == 0) THEN
          WRITE (nerr,'(a,i1,a1,i1)') &
               '  Used MPI version: ', version, '.', subversion
       END IF

       ! due to a possible circular dependency with mo_machine and other
       ! modules, we determine here locally the I/O size of the different
       ! kind types (assume 8 bit/byte. This is than used for determing
       ! the right MPI send/receive type parameters.

       ! first get the native INTEGER size

       integer_byte_size = BIT_SIZE(iig)/8

       ! and inquire for the I/O size (is independent of byte or word
       ! values ...)

       INQUIRE (iolength=io_size) iig
       integer_io_size = io_size
       p_ig = io_size/integer_io_size*integer_byte_size

       ! and the native REAL size

       INQUIRE (iolength=io_size) rrg
       p_rg = io_size/integer_io_size*integer_byte_size

       ! find now the size of usual 4 byte and 8 byte INTEGER
       ! (might be 8 byte both, or only the 4 byte available ...

       INQUIRE (iolength=io_size) ii4
       p_i4 = io_size/integer_io_size*integer_byte_size
       INQUIRE (iolength=io_size) ii8
       p_i8 = io_size/integer_io_size*integer_byte_size

       ! find now the size of usual 4 byte and 8 byte REAL
       ! (might be 8 byte both)

       INQUIRE (iolength=io_size) rsp
       p_sp = io_size/integer_io_size*integer_byte_size
       INQUIRE (iolength=io_size) rdp
       p_dp = io_size/integer_io_size*integer_byte_size

       ! testing this variables

       p_int_i4  = p_type_kind (i4, integer_type)
       p_int_i8  = p_type_kind (i8, integer_type)
       p_real_sp = p_type_kind (sp, real_type)
       p_real_dp = p_type_kind (dp, real_type)
       p_real_wp = p_type_kind (wp, real_type)

       IF (mype == 0) THEN

          IF (p_ig == p_i4) THEN
             p_int = p_int_i4
          ELSE IF (p_ig == p_i8) THEN
             p_int = p_int_i8
          END IF

          IF (p_rg == p_sp) THEN
             p_real = p_real_sp
          ELSE IF (p_rg == p_dp) THEN
             p_real = p_real_dp
          END IF

       END IF

       p_send(1) = p_int
       p_send(2) = p_int_i4
       p_send(3) = p_int_i8
       p_send(4) = 0         ! not used
       p_send(5) = p_real
       p_send(6) = p_real_sp
       p_send(7) = p_real_dp
       p_send(8) = p_real_wp

       CALL MPI_BCAST (p_send, 8, MPI_INTEGER, 0, MPI_COMM_WORLD, p_error)

#ifdef DEBUG
       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_BCAST for send/receive types failed.'
          WRITE (nerr,'(a,i4)') ' Error = ', p_error
       END IF
#endif

       p_int      = p_send(1)
       p_int_i4   = p_send(2)
       p_int_i8   = p_send(3)
       p_real     = p_send(5)
       p_real_sp  = p_send(6)
       p_real_dp  = p_send(7)
       p_real_wp  = p_send(8)

       ! set logical and character types to native types

       p_bool = MPI_LOGICAL
       p_char = MPI_CHARACTER

       IF (mype == 0) THEN
          WRITE (nerr,'(/)')
          IF (p_real_sp == MPI_REAL) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_REAL for KIND sp'
          ELSE IF (p_real_sp == MPI_DOUBLE_PRECISION) THEN
             WRITE (nerr,'(a)') &
                  ' Selected type: MPI_DOUBLE_PRECISION for KIND sp'
          END IF

          IF (p_real_dp == MPI_DOUBLE_PRECISION) THEN
             WRITE (nerr,'(a)') &
                  ' Selected type: MPI_DOUBLE_PRECISION for KIND dp'
          END IF

          IF (p_int_i4 == MPI_INTEGER4) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_INTEGER4 for KIND i4'
          ELSE IF (p_int_i4 == MPI_INTEGER8) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_INTEGER8 for KIND i4'
          END IF

          IF (p_int_i8 == MPI_INTEGER8) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_INTEGER8 for KIND i8'
          END IF
       END IF
    END IF
#endif

    IF (mype == 0) WRITE (nerr,'(/)')

#ifdef DEBUG
    IF (mype == 0) then
      WRITE (nerr,'(a)')    ' Transfer types:'
      WRITE (nerr,'(a,i0)') '  INTEGER generic: ', p_int
      WRITE (nerr,'(a,i0)') '  INTEGER 4 byte : ', p_int_i4
      WRITE (nerr,'(a,i0)') '  INTEGER 8 byte : ', p_int_i8
      WRITE (nerr,'(a,i0)') '  REAL generic   : ', p_real
      WRITE (nerr,'(a,i0)') '  REAL single    : ', p_real_sp
      WRITE (nerr,'(a,i0)') '  REAL double    : ', p_real_dp
      WRITE (nerr,'(/)')
    end IF
#endif

  !--------------------------------------------------------------
  ! set communication group derived type variable: MPI_COMM_WORLD
  !--------------------------------------------------------------
  world% comm = MPI_COMM_WORLD
  world% pe   = p_pe
  world% pio  = p_io
  world% lpio = p_pe == p_io
  world% npe  = p_nprocs
  !-------------------------------------------------------------
  ! set communication group derived type variable: MPI_COMM_SELF
  !-------------------------------------------------------------
  self % comm = MPI_COMM_SELF
  self % pe   = 0
  self % pio  = 0
  self % lpio = .true.
  self % npe  = 1
  !----------------------------------------------------
  ! set communication group derived type variable: DACE
  ! generally same as MPI_COMM_WORLD
  !----------------------------------------------------
  dace = world
  !------------------------------------
  ! for tests: exclude last PE for DACE
  !------------------------------------
#ifndef NOMPI
  if (.FALSE.) then
    dace % comm = -1
    dace % pe   = -1
    call MPI_COMM_GROUP  (MPI_COMM_WORLD, group_world, p_error)
    call MPI_GROUP_EXCL  (group_world, 1, [p_nprocs-1], group_dace, p_error)  ! last PE not member
    call MPI_COMM_CREATE (MPI_COMM_WORLD, group_dace, dace% comm, p_error)
    call set_dace_comm (dace% comm, 0)

!   if (dace% pe == dace% pio) then
!     do i=1,p_nprocs
!       print *,p_pe,'###',world,' - ',self,' - ',dace,' - ',group_dace, MPI_COMM_NULL
!     end do
!   endif
  endif
#endif
  !--------------------------------------------
  ! make copies for include-files used in COSMO
  !--------------------------------------------
  d_npe       = dace% npe
  d_comm      = dace% comm
  d_pe        = dace% pe

  END SUBROUTINE p_start
!------------------------------------------------------------------------------
  subroutine set_dace_comm (comm, p_io)
  integer ,intent(in) :: comm   ! communicator to use
  integer ,intent(in) :: p_io   ! processor to use for I/O
  !---------------------------------------
  ! set the DACE communicator and defaults
  !---------------------------------------
    integer :: ierr
    if  (comm /= MPI_COMM_NULL) then
#ifndef NOMPI
      call MPI_COMM_SIZE (comm, dace% npe, ierr)
      call MPI_COMM_RANK (comm, dace% pe,  ierr)
      dace % comm = comm
      dace % pio  = p_io
      dace % lpio = dace% pe == dace% pio
#else
      dace = t_comm ()  ! Default communicator for single-processor run
#endif
    else
      dace % comm = MPI_COMM_NULL
      dace % pe   = -1
      dace % pio  =  0
      dace % lpio = .false.
      dace % npe  = 0
    endif
    !--------------------------------------------
    ! make copies for include-files used in COSMO
    !--------------------------------------------
    d_npe       = dace% npe
    d_comm      = dace% comm
    d_pe        = dace% pe
    !-------------------------------------------------------
    ! initialize miscellaneous variables used by this module
    !-------------------------------------------------------
    p_pe        = dace% pe              ! Legacy
    mype        = dace% pe
    p_parallel  = dace% npe > 1
    if (p_mrequest == 0) then
       p_mrequest = 100 * dace% npe
       ALLOCATE (p_request(p_mrequest))
       p_irequest = 1
    end if
  end subroutine set_dace_comm
!------------------------------------------------------------------------------
  SUBROUTINE p_stop
  !------------------------------------------------
  ! finish MPI and clean up all PEs.
  ! should be called before the end of the program.
  !------------------------------------------------

#ifndef NOMPI

    CALL p_barrier (MPI_COMM_WORLD) ! to prevent abort due to unfinished communication

    CALL MPI_FINALIZE (p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_FINALIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
    p_parallel = .FALSE.
    DEALLOCATE (p_request)
#endif
!   STOP
  END SUBROUTINE p_stop
!------------------------------------------------------------------------------
  SUBROUTINE p_abort
  !-------------------------------------------------------------------
  ! this routine should be used instead of abort, util_abort() or STOP
  ! in all routines for proper clean up of all PEs
  !-------------------------------------------------------------------

!    EXTERNAL util_exit

#ifndef NOMPI
    CALL MPI_ABORT (MPI_COMM_WORLD, 1, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ABORT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
#else
!   CALL util_exit(1)
    CALL abort('p_abort called')
#endif

  END SUBROUTINE p_abort
!------------------------------------------------------------------------------

  FUNCTION p_type_kind (kind_type, var_type) RESULT (p_message_type)

    USE mo_kind

    INTEGER              :: p_message_type
    INTEGER, INTENT(in)  :: kind_type, var_type

    IF (var_type == integer_type) THEN
       IF (kind_type == i8) THEN
          p_message_type = check_type_i8 ()
       ELSE IF (kind_type == i4) THEN
          p_message_type = check_type_i4 ()
       END IF
    ELSE IF (var_type == real_type) THEN
       IF (kind_type == dp) THEN
          p_message_type = check_type_dp ()
       ELSE IF (kind_type == sp) THEN
          p_message_type = check_type_sp ()
       END IF
    END IF

  END FUNCTION p_type_kind

  FUNCTION p_type_i4 (int)
  !--------------------------------
  ! return MPI_type for integer(i4)
  !--------------------------------
  INTEGER(i4), INTENT(in) :: int
  INTEGER                 :: p_type_i4
    p_type_i4 = p_type_kind (i4, integer_type)
  END FUNCTION p_type_i4

  FUNCTION p_type_i8 (int)
  !--------------------------------
  ! return MPI_type for integer(i8)
  !--------------------------------
  INTEGER(i8), INTENT(in) :: int
  INTEGER                 :: p_type_i8
    p_type_i8 = p_type_kind (i8, integer_type)
  END FUNCTION p_type_i8

  FUNCTION p_type_sp (r)
  !-----------------------------
  ! return MPI_type for real(sp)
  !-----------------------------
  REAL(sp), INTENT(in) :: r
  INTEGER              :: p_type_sp
    p_type_sp = p_type_kind (sp, real_type)
  END FUNCTION p_type_sp

  FUNCTION p_type_dp (r)
  !-----------------------------
  ! return MPI_type for real(dp)
  !-----------------------------
  REAL(dp), INTENT(in) :: r
  INTEGER              :: p_type_dp
    p_type_dp = p_type_kind (dp, real_type)
  END FUNCTION p_type_dp

  FUNCTION check_type_i4 () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    INTEGER :: datatype

    INTEGER (i4) :: buf_int_i4(1,8)
    INTEGER (i4) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN

       a = HUGE(a)
       buf_int_i4(1,1) = a

       CALL MPI_SEND (buf_int_i4(1,1), 1, MPI_INTEGER4, 1, 1, &
            MPI_COMM_WORLD, p_error)
       CALL MPI_RECV (buf_int_i4(1,5), 1, MPI_INTEGER4, 1, 2, &
            MPI_COMM_WORLD, p_status, p_error)
       b = buf_int_i4(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (MPI_BYTE, 1, MPI_INTEGER, 1, 0, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_SEND (buf_int_i4(1,1), p_i4, MPI_BYTE, 1, 3, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_RECV (buf_int_i4(1,5), p_i4, MPI_BYTE, 1, 4, &
               MPI_COMM_WORLD, p_status, p_error)
          c =  buf_int_i4(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') &
                  ' Warning: MPI_INTEGER4 and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_INTEGER4 anyway'
          END IF

          p_type = MPI_INTEGER4

       ELSE
          CALL MPI_SEND (MPI_INTEGER8, 1, MPI_INTEGER, 1, 0, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_SEND (buf_int_i4(1,1), 1, MPI_INTEGER8, 1, 5, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_RECV (buf_int_i4(1,5), 1, MPI_INTEGER8, 1, 6, &
               MPI_COMM_WORLD, p_status, p_error)
          b = buf_int_i4(1,5)

          IF (a == b) THEN
             CALL MPI_SEND (buf_int_i4(1,1), p_i4, MPI_BYTE, 1, 7, &
                  MPI_COMM_WORLD, p_error)
             CALL MPI_RECV (buf_int_i4(1,5), p_i4, MPI_BYTE, 1, 8, &
                  MPI_COMM_WORLD, p_status, p_error)
             c =  buf_int_i4(1,5)

             IF (a /= c) THEN
                WRITE (nerr,'(a)') &
                     ' Warning: MPI_INTEGER8 and MPI_BYTE not equivalent'
                WRITE (nerr,'(a)') &
                     ' Using MPI_INTEGER8 anyway'
             END IF

             p_type = MPI_INTEGER8

          END IF
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_int_i4(1,7), 1, MPI_INTEGER4, 0, 1, &
            MPI_COMM_WORLD, p_status, p_error)
       buf_int_i4(1,3) = buf_int_i4(1,7)
       CALL MPI_SEND (buf_int_i4(1,3), 1, MPI_INTEGER4, 0, 2, &
            MPI_COMM_WORLD, p_error)

       CALL MPI_RECV (datatype, 1, MPI_INTEGER, 0, 0, &
            MPI_COMM_WORLD, p_status, p_error)

       IF (datatype == MPI_BYTE) THEN
          CALL MPI_RECV (buf_int_i4(1,7), p_i4, MPI_BYTE, 0, 3, &
               MPI_COMM_WORLD, p_status, p_error)
          buf_int_i4(1,3) = buf_int_i4(1,7)
          CALL MPI_SEND (buf_int_i4(1,3), p_i4, MPI_BYTE, 0, 4, &
               MPI_COMM_WORLD, p_error)
       ELSE IF (datatype == MPI_INTEGER8) THEN
          CALL MPI_RECV (buf_int_i4(1,7), 1, MPI_INTEGER8, 0, 5, &
               MPI_COMM_WORLD, p_status, p_error)
          buf_int_i4(1,3) = buf_int_i4(1,7)
          CALL MPI_SEND (buf_int_i4(1,3), 1, MPI_INTEGER8, 0, 6, &
               MPI_COMM_WORLD, p_error)
       END IF

       IF (datatype == MPI_INTEGER8) THEN
          CALL MPI_RECV (buf_int_i4(1,7), p_i4, MPI_BYTE, 0, 7, &
               MPI_COMM_WORLD, p_status, p_error)
          buf_int_i4(1,3) = buf_int_i4(1,7)
          CALL MPI_SEND (buf_int_i4(1,3), p_i4, MPI_BYTE, 0, 8, &
               MPI_COMM_WORLD, p_error)
       END IF
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_i4

  FUNCTION check_type_i8 () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    INTEGER (i8) :: buf_int_i8(1,8)
    INTEGER (i8) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN

       a = HUGE(a)
       buf_int_i8(1,1) = a

       CALL MPI_SEND (buf_int_i8(1,1), 1, MPI_INTEGER8, 1, 1, &
            MPI_COMM_WORLD, p_error)
       CALL MPI_RECV (buf_int_i8(1,5), 1, MPI_INTEGER8, 1, 2, &
            MPI_COMM_WORLD, p_status, p_error)
       b = buf_int_i8(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (buf_int_i8(1,1), p_i8, MPI_BYTE, 1, 3, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_RECV (buf_int_i8(1,5), p_i8, MPI_BYTE, 1, 4, &
               MPI_COMM_WORLD, p_status, p_error)
          c =  buf_int_i8(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') &
                  ' Warning: MPI_INTEGER8 and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_INTEGER8 anyway'
          END IF

          p_type = MPI_INTEGER8

       ELSE
          WRITE (nerr,'(a)') ' MPI_INTEGER8 not available.'
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_int_i8(1,7), 1, MPI_INTEGER8, 0, 1, &
            MPI_COMM_WORLD, p_status, p_error)
       buf_int_i8(1,3) = buf_int_i8(1,7)
       CALL MPI_SEND (buf_int_i8(1,3), 1, MPI_INTEGER8, 0, 2, &
            MPI_COMM_WORLD, p_error)

       CALL MPI_RECV (buf_int_i8(1,7), p_i8, MPI_BYTE, 0, 3, &
            MPI_COMM_WORLD, p_status, p_error)
       buf_int_i8(1,3) = buf_int_i8(1,7)
       CALL MPI_SEND (buf_int_i8(1,3), p_i8, MPI_BYTE, 0, 4, &
            MPI_COMM_WORLD, p_error)
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_i8

  FUNCTION check_type_sp () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    INTEGER :: datatype

    REAL (sp) :: buf_real_sp(1,8)
    REAL (sp) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN

       a = HUGE(a)
       buf_real_sp(1,1) = a

       CALL MPI_SEND (buf_real_sp(1,1), 1, MPI_REAL, 1, 1, &
            MPI_COMM_WORLD, p_error)
       CALL MPI_RECV (buf_real_sp(1,5), 1, MPI_REAL, 1, 2, &
            MPI_COMM_WORLD, p_status, p_error)
       b = buf_real_sp(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (MPI_BYTE, 1, MPI_INTEGER, 1, 0, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_SEND (buf_real_sp(1,1), p_sp, MPI_BYTE, 1, 3, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_RECV (buf_real_sp(1,5), p_sp, MPI_BYTE, 1, 4, &
               MPI_COMM_WORLD, p_status, p_error)
          c =  buf_real_sp(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') ' Warning: MPI_REAL and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_REAL anyway'
          END IF

          p_type = MPI_REAL

       ELSE
          CALL MPI_SEND (MPI_DOUBLE_PRECISION, 1, MPI_INTEGER, 1, 0, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_SEND (buf_real_sp(1,1), 1, MPI_DOUBLE_PRECISION, 1, 5, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_RECV (buf_real_sp(1,5), 1, MPI_DOUBLE_PRECISION, 1, 6, &
               MPI_COMM_WORLD, p_status, p_error)
          b = buf_real_sp(1,5)

          IF (a == b) THEN
             CALL MPI_SEND (buf_real_sp(1,1), p_sp, MPI_BYTE, 1, 7, &
                  MPI_COMM_WORLD, p_error)
             CALL MPI_RECV (buf_real_sp(1,5), p_sp, MPI_BYTE, 1, 8, &
                  MPI_COMM_WORLD, p_status, p_error)
             c =  buf_real_sp(1,5)

             IF (a /= c) THEN
                WRITE (nerr,'(a,a)') &
                     ' Warning: MPI_DOUBLE_PRECISION and MPI_BYTE ', &
                     'not equivalent'
                WRITE (nerr,'(a)') &
                     ' Using MPI_DOUBLE_PRECISION anyway'
             END IF

             p_type = MPI_DOUBLE_PRECISION

          END IF
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_real_sp(1,7), 1, MPI_REAL, 0, 1, &
            MPI_COMM_WORLD, p_status, p_error)
       buf_real_sp(1,3) = buf_real_sp(1,7)
       CALL MPI_SEND (buf_real_sp(1,3), 1, MPI_REAL, 0, 2, &
            MPI_COMM_WORLD, p_error)

       CALL MPI_RECV (datatype, 1, MPI_INTEGER, 0, 0, &
            MPI_COMM_WORLD, p_status, p_error)

       IF (datatype == MPI_BYTE) THEN
          CALL MPI_RECV (buf_real_sp(1,7), p_sp, MPI_BYTE, 0, 3, &
               MPI_COMM_WORLD, p_status, p_error)
          buf_real_sp(1,3) = buf_real_sp(1,7)
          CALL MPI_SEND (buf_real_sp(1,3), p_sp, MPI_BYTE, 0, 4, &
               MPI_COMM_WORLD, p_error)
       ELSE IF (datatype == MPI_DOUBLE_PRECISION) THEN
          CALL MPI_RECV (buf_real_sp(1,7), 1, MPI_DOUBLE_PRECISION, 0, 5, &
               MPI_COMM_WORLD, p_status, p_error)
          buf_real_sp(1,3) = buf_real_sp(1,7)
          CALL MPI_SEND (buf_real_sp(1,3), 1, MPI_DOUBLE_PRECISION, 0, 6, &
               MPI_COMM_WORLD, p_error)
       END IF

       IF (datatype == MPI_DOUBLE_PRECISION) THEN
          CALL MPI_RECV (buf_real_sp(1,7), p_sp, MPI_BYTE, 0, 7, &
               MPI_COMM_WORLD, p_status, p_error)
          buf_real_sp(1,3) = buf_real_sp(1,7)
          CALL MPI_SEND (buf_real_sp(1,3), p_sp, MPI_BYTE, 0, 8, &
               MPI_COMM_WORLD, p_error)
       END IF
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_sp

  FUNCTION check_type_dp () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    REAL (dp) :: buf_real_dp(1,8)
    REAL (dp) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    buf_real_dp = 0
    IF (mype == 0) THEN

       a = HUGE(a)
       buf_real_dp(1,1) = a

       CALL MPI_SEND (buf_real_dp(1,1), 1, MPI_DOUBLE_PRECISION, 1, 1, &
            MPI_COMM_WORLD, p_error)
       CALL MPI_RECV (buf_real_dp(1,5), 1, MPI_DOUBLE_PRECISION, 1, 2, &
            MPI_COMM_WORLD, p_status, p_error)
       b = buf_real_dp(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (buf_real_dp(1,1), p_dp, MPI_BYTE, 1, 3, &
               MPI_COMM_WORLD, p_error)
          CALL MPI_RECV (buf_real_dp(1,5), p_dp, MPI_BYTE, 1, 4, &
               MPI_COMM_WORLD, p_status, p_error)
          c =  buf_real_dp(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') &
                  ' Warning: MPI_DOUBLE_PRECISION and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_DOUBLE_PRECISION anyway'
          END IF

          p_type = MPI_DOUBLE_PRECISION

       ELSE
          WRITE (nerr,'(a)') ' MPI_DOUBLE_PRECISION not available.'
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_real_dp(1,7), 1, MPI_DOUBLE_PRECISION, 0, 1, &
            MPI_COMM_WORLD, p_status, p_error)
       buf_real_dp(1,3) = buf_real_dp(1,7)
       CALL MPI_SEND (buf_real_dp(1,3), 1, MPI_DOUBLE_PRECISION, 0, 2, &
            MPI_COMM_WORLD, p_error)

       CALL MPI_RECV (buf_real_dp(1,7), p_dp, MPI_BYTE, 0, 3, &
            MPI_COMM_WORLD, p_status, p_error)
       buf_real_dp(1,3) = buf_real_dp(1,7)
       CALL MPI_SEND (buf_real_dp(1,3), p_dp, MPI_BYTE, 0, 4, &
            MPI_COMM_WORLD, p_error)
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_dp

!------------------------------------------------------------------------------
!  SUBROUTINE p_set_communicator (nproca, nprocb, mapmesh, debug_parallel)
!  !-----------------------------------------------------------
!  ! set up the communicator for atmospheric field partitioning
!  !-----------------------------------------------------------
!  INTEGER, INTENT(in) :: nproca
!  INTEGER, INTENT(in) :: nprocb
!  INTEGER, INTENT(in) :: mapmesh (0:,0:)
!  INTEGER, INTENT(in) :: debug_parallel
!
!#ifndef NOMPI
!    INTEGER :: all_debug_pes(SIZE(mapmesh))
!
!    INTEGER :: group_world, group_a, group_b, group_d
!    INTEGER :: p_communicator_tmp
!
!    INTEGER :: n, members
!
!    ! first set global group
!
!    CALL MPI_COMM_GROUP (MPI_COMM_WORLD, group_world, p_error)
!
!    IF (p_error /= MPI_SUCCESS) THEN
!       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_GROUP failed.'
!       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!       CALL p_abort
!    END IF
!
!    ! communicator is MPI_COMM_WORLD
!
!    IF (debug_parallel >= 0 ) THEN
!
!       CALL MPI_GROUP_INCL (group_world, 1, [0], group_d, p_error)
!
!       IF (p_error /= MPI_SUCCESS) THEN
!          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
!          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!          CALL p_abort
!       END IF
!
!       CALL MPI_COMM_CREATE (MPI_COMM_WORLD, group_d, p_communicator_tmp, &
!            p_error)
!
!       IF (p_error /= MPI_SUCCESS) THEN
!          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
!          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!          CALL p_abort
!       END IF
!
!       IF (mype == 0) p_communicator_d = p_communicator_tmp
!
!       DO n = 1, SIZE(mapmesh)
!          all_debug_pes(n) = n
!       END DO
!
!       CALL MPI_GROUP_INCL (group_world, SIZE(mapmesh), all_debug_pes, &
!            group_d, p_error)
!
!       IF (p_error /= MPI_SUCCESS) THEN
!          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
!          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!          CALL p_abort
!       END IF
!
!       CALL MPI_COMM_CREATE (MPI_COMM_WORLD, group_d, p_communicator_tmp, &
!            p_error)
!
!       IF (p_error /= MPI_SUCCESS) THEN
!          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
!          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!          CALL p_abort
!       END IF
!
!       IF (mype /= 0) p_communicator_d = p_communicator_tmp
!
!    ELSE
!       p_communicator_d = MPI_COMM_WORLD
!    END IF
!
!    DO n = 0, nproca-1
!       members = nprocb
!       CALL MPI_GROUP_INCL (group_world, members, mapmesh(:,n), group_a, &
!            p_error)
!
!       IF (p_error /= MPI_SUCCESS) THEN
!          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
!          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!          CALL p_abort
!       END IF
!
!       CALL MPI_COMM_CREATE (MPI_COMM_WORLD, group_a, p_communicator_tmp, &
!            p_error)
!
!       IF (p_error /= MPI_SUCCESS) THEN
!          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
!          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!          CALL p_abort
!       END IF
!       IF(p_communicator_tmp/=MPI_COMM_NULL) &
!         p_communicator_a = p_communicator_tmp
!
!    END DO
!
!    ! create groups for set Bs
!
!    DO n = 0, nprocb-1
!       members = nproca
!       CALL MPI_GROUP_INCL (group_world, members, mapmesh(n,:), group_b, &
!            p_error)
!
!       IF (p_error /= MPI_SUCCESS) THEN
!          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
!          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!          CALL p_abort
!       END IF
!
!       CALL MPI_COMM_CREATE (MPI_COMM_WORLD, group_b, p_communicator_tmp, &
!            p_error)
!
!       IF (p_error /= MPI_SUCCESS) THEN
!          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
!          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!          CALL p_abort
!       END IF
!       IF(p_communicator_tmp/=MPI_COMM_NULL) &
!         p_communicator_b = p_communicator_tmp
!
!    END DO
!
!    CALL MPI_BARRIER (MPI_COMM_WORLD, p_error)
!
!    IF (p_error /= MPI_SUCCESS) THEN
!       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_BARRIER failed.'
!       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!       CALL p_abort
!    END IF
!
!    IF (debug_parallel >= 0 .AND. mype == 0) THEN
!      p_communicator_a = p_communicator_d
!      p_communicator_b = p_communicator_d
!    ENDIF
!
!    WRITE (nerr,'(a,i4,a,3i8)') &
!         'p_set_communicator on PE ', mype, ': ', &
!         p_communicator_d, &
!         p_communicator_a, &
!         p_communicator_b
!#endif
!  END SUBROUTINE p_set_communicator
!
!=========================================================================

  ! send implementation

  SUBROUTINE p_send_real (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real

  SUBROUTINE p_send_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_1d

  SUBROUTINE p_send_real_1d_sp (buffer, p_destination, p_tag, p_count, comm)

    REAL (sp), INTENT(in) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, P_REAL_SP, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), P_REAL_SP, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_1d_sp

  SUBROUTINE p_send_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_2d

  SUBROUTINE p_send_real_2d_sp (buffer, p_destination, p_tag, p_count, comm)

    REAL (sp), INTENT(in) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, P_REAL_SP, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), P_REAL_SP, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_2d_sp

  SUBROUTINE p_send_real_3d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_3d

  SUBROUTINE p_send_real_4d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_4d

  SUBROUTINE p_send_real_5d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_5d

  SUBROUTINE p_send_int (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int

  SUBROUTINE p_send_int_1d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm, count
#ifdef MPI_CHECKSUM
    integer :: checksum, chk
#endif

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       count = p_count
    ELSE
       count = size (buffer)
    END IF
    CALL MPI_SEND (buffer, count, p_int, p_destination, p_tag, &
                   p_comm, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#ifdef MPI_CHECKSUM
    checksum = crc (buffer(1:count))
    call p_send (checksum, p_destination, p_tag=32)
    call p_recv (chk     , p_destination, p_tag=33)
    if (checksum /= chk) then
       write (nerr,'()')
       write (nerr,'(a,i12,i3,a,i12,i3,1x,i0)') &
            "p_send_int_1d:BAD CHECKSUM:", chk, p_pe, &
            "  dest:", checksum, p_destination, count
       call dump_buffer_int (buffer(1:count), unit=2000+p_pe, &
                             title="p_send_int_1d: sent data")
       call p_send (chk,   p_destination, p_tag=34)  ! Synchronize before crash
       CALL p_abort
    end if
#endif /* MPI_CHECKSUM */
#endif /* NOMPI */

  END SUBROUTINE p_send_int_1d

!------------------------------------------------------------------------------
#ifdef MPI_CHECKSUM
  subroutine dump_buffer_int (v, unit, title)
    USE mo_system, ONLY: flush
    integer,          intent(in)           :: v(:)
    integer,          intent(in), optional :: unit
    character(len=*), intent(in), optional :: title
    integer :: n, i, j, lu
    lu = 0; if (present (unit)) lu = unit
    write (0,'()')
    write (0,'(a,i5)') "Dumping buffer to fortran unit", lu
    write (0,'()')
    if (lu>0 .and. present (title)) then
       write (lu,*) title
    end if
    n = size (v)
    write (lu,*) "pe=", p_pe, ": Buffer size (words):", n
    write (lu,'()')
    write (lu,*) "pe=", p_pe, ": Data:", v(:min (n,4)), "...", v(max (5,n-3):)
    write (lu,'()')
    if (lu == 0) return       ! Do not flood stderr...
    write (lu,'(A)') "Nontrivial contents:"
    do i = 1, n, 8
       j = min (i+7,n)
       if (any (v(i:j) /= 0)) write (lu,*) i,"..",j,":", v(i:j)
    end do
    call flush (lu)
#ifdef __SX__
    call SLEEP (5)
#endif
  end subroutine dump_buffer_int
#endif
!------------------------------------------------------------------------------

  SUBROUTINE p_send_int_1d_i2 (buffer, p_destination, p_tag, p_count, comm)

    INTEGER(i2), INTENT(in) :: buffer(:)
    INTEGER,     INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int_i2, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int_i2, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_1d_i2

!------------------------------------------------------------------------------

#ifdef HAVE_I1
  SUBROUTINE p_send_int_1d_i1 (buffer, p_destination, p_tag, p_count, comm)

    INTEGER(i1), INTENT(in) :: buffer(:)
    INTEGER,     INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, MPI_BYTE, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), MPI_BYTE, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_1d_i1
#endif

!------------------------------------------------------------------------------

  SUBROUTINE p_send_int_2d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_2d

  SUBROUTINE p_send_int_3d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_3d

  SUBROUTINE p_send_int_4d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_4d


  SUBROUTINE p_send_bool (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool

  SUBROUTINE p_send_bool_1d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_1d

  SUBROUTINE p_send_bool_2d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_2d

  SUBROUTINE p_send_bool_3d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_3d

  SUBROUTINE p_send_bool_4d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_4d

  SUBROUTINE p_send_char (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER (*), INTENT(in) :: buffer
    INTEGER,       INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, LEN(buffer), p_char, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_char

  SUBROUTINE p_send_char_1d (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER (*), INTENT(in) :: buffer (:)
    INTEGER,       INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm, l

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       l = LEN(buffer) * p_count
    ELSE
       l = LEN(buffer) * SIZE(buffer)
    END IF
    CALL MPI_SEND (buffer, l, p_char, p_destination, p_tag, p_comm, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_char_1d

! non-blocking sends

  SUBROUTINE p_isend_real (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in)           :: buffer
    INTEGER,   INTENT(in)           :: p_destination, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, 1, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real

  SUBROUTINE p_isend_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in)           :: buffer(:)
    INTEGER,   INTENT(in)           :: p_destination, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_1d

  SUBROUTINE p_isend_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in)           :: buffer(:,:)
    INTEGER,   INTENT(in)           :: p_destination, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_2d

  SUBROUTINE p_isend_real_3d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in)           :: buffer(:,:,:)
    INTEGER,   INTENT(in)           :: p_destination, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_3d

  SUBROUTINE p_isend_real_4d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in)           :: buffer(:,:,:,:)
    INTEGER,   INTENT(in)           :: p_destination, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_4d

  SUBROUTINE p_isend_real_5d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in)           :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)           :: p_destination, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), P_REAL_DP, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_5d

  SUBROUTINE p_isend_int (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in)           :: buffer
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int

  SUBROUTINE p_isend_int_1d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in)           :: buffer(:)
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_1d

  SUBROUTINE p_isend_int_2d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in)           :: buffer(:,:)
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_2d

  SUBROUTINE p_isend_int_3d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in)           :: buffer(:,:,:)
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_3d

  SUBROUTINE p_isend_int_4d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in)           :: buffer(:,:,:,:)
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_4d


  SUBROUTINE p_isend_bool (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in)           :: buffer
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool

  SUBROUTINE p_isend_bool_1d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in)           :: buffer(:)
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_1d

  SUBROUTINE p_isend_bool_2d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in)           :: buffer(:,:)
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_2d

  SUBROUTINE p_isend_bool_3d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in)           :: buffer(:,:,:)
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_3d

  SUBROUTINE p_isend_bool_4d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in)           :: buffer(:,:,:,:)
    INTEGER, INTENT(in)           :: p_destination, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    ELSE
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
       p_irequest = p_irequest + 1
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_4d

  SUBROUTINE p_isend_char (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER (*), INTENT(in)           :: buffer
    INTEGER,       INTENT(in)           :: p_destination, p_tag
    INTEGER,       INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm, l

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       l = p_count
    ELSE
       l = LEN (buffer)
    END IF
    CALL MPI_ISEND (buffer, l, p_char, p_destination, p_tag, &
                    p_comm, p_request(p_irequest), p_error   )
    p_irequest = p_irequest + 1

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_char

  SUBROUTINE p_isend_char_1d (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER (*), INTENT(in)           :: buffer (:)
    INTEGER,       INTENT(in)           :: p_destination, p_tag
    INTEGER,       INTENT(in), OPTIONAL :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm, l

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       l = LEN (buffer) * p_count
    ELSE
       l = LEN (buffer) * SIZE (buffer)
    END IF
    CALL MPI_ISEND (buffer, l, p_char, p_destination, p_tag, &
                    p_comm, p_request(p_irequest), p_error   )
    p_irequest = p_irequest + 1

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_char_1d

  ! recv implementation

  SUBROUTINE p_recv_real (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out)          :: buffer
    INTEGER,   INTENT(in)           :: p_source, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real

  SUBROUTINE p_recv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_1d

  SUBROUTINE p_recv_real_1d_sp (buffer, p_source, p_tag, p_count, comm)

    REAL (sp), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, P_REAL_SP, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), P_REAL_SP, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_1d_sp

  SUBROUTINE p_recv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_2d

  SUBROUTINE p_recv_real_2d_sp (buffer, p_source, p_tag, p_count, comm)

    REAL (sp), INTENT(out) :: buffer(:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, P_REAL_SP, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), P_REAL_SP, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_2d_sp

  SUBROUTINE p_recv_real_3d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_3d

  SUBROUTINE p_recv_real_4d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_4d

  SUBROUTINE p_recv_real_5d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), P_REAL_DP, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_5d

  SUBROUTINE p_recv_int (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int

  SUBROUTINE p_recv_int_1d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm, count
#ifdef MPI_CHECKSUM
    integer :: chk, checksum
#endif

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       count = p_count
    ELSE
       count = size (buffer)
    END IF
    CALL MPI_RECV (buffer, count, p_int, p_source, p_tag, &
                   p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#ifdef MPI_CHECKSUM
    chk = crc (buffer(1:count))
    call p_recv (checksum, p_source, p_tag=32)
    call p_send (chk     , p_source, p_tag=33)
    if (checksum /= chk) then
       write (nerr,'()')
       write (nerr,'(a,i12,i3,a,i12,i3,1x,i0)') &
            "p_recv_int_1d:BAD CHECKSUM:", chk, p_pe, &
            "  source:", checksum, p_source, count
       call dump_buffer_int (buffer(1:count), unit=3000+p_pe, &
                             title="p_recv_int_1d: received data")
       call p_recv (chk,   p_source, p_tag=34)  ! Synchronize before crash
       CALL p_abort
    end if
#endif /* MPI_CHECKSUM */
#endif /* NOMPI */

  END SUBROUTINE p_recv_int_1d

!------------------------------------------------------------------------------

  SUBROUTINE p_recv_int_1d_i2 (buffer, p_source, p_tag, p_count, comm)

    INTEGER(i2), INTENT(out) :: buffer(:)
    INTEGER,     INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int_i2, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int_i2, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_1d_i2

!------------------------------------------------------------------------------

#ifdef HAVE_I1
  SUBROUTINE p_recv_int_1d_i1 (buffer, p_source, p_tag, p_count, comm)

    INTEGER(i1), INTENT(out) :: buffer(:)
    INTEGER,     INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, MPI_BYTE, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), MPI_BYTE, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_1d_i1
#endif

!------------------------------------------------------------------------------

  SUBROUTINE p_recv_int_2d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_2d

  SUBROUTINE p_recv_int_3d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_3d

  SUBROUTINE p_recv_int_4d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_4d


  SUBROUTINE p_recv_bool (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool

  SUBROUTINE p_recv_bool_1d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_1d

  SUBROUTINE p_recv_bool_2d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_2d

  SUBROUTINE p_recv_bool_3d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_3d

  SUBROUTINE p_recv_bool_4d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_4d

  SUBROUTINE p_recv_char (buffer, p_source, p_tag, p_count, comm)

    CHARACTER (*), INTENT(out) :: buffer
    INTEGER,       INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, LEN(buffer), p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_char

  SUBROUTINE p_recv_char_1d (buffer, p_source, p_tag, p_count, comm)

    CHARACTER (*), INTENT(out) :: buffer (:)
    INTEGER,       INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm, l

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       l = LEN(buffer) * p_count
    ELSE
       l = LEN(buffer) * SIZE(buffer)
    END IF
    CALL MPI_RECV (buffer,l,p_char,p_source,p_tag,p_comm,p_status,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_char_1d

  ! non-blocking receives

  SUBROUTINE p_irecv_real (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout)        :: buffer
    INTEGER,   INTENT(in)           :: p_source, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm, l

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       l = p_count
    ELSE
       l = 1
    END IF

    CALL MPI_IRECV (buffer, l, p_real_dp, p_source, p_tag, &
                    p_comm, p_request(p_irequest), p_error)
    p_irequest = p_irequest + 1

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real

  SUBROUTINE p_irecv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout)        :: buffer(:)
    INTEGER,   INTENT(in)           :: p_source, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm, l

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       l = p_count
    ELSE
       l = SIZE(buffer)
    END IF

    CALL MPI_IRECV (buffer, l, p_real_dp, p_source, p_tag, &
                    p_comm, p_request(p_irequest), p_error)
    p_irequest = p_irequest + 1

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_1d

  SUBROUTINE p_irecv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout)        :: buffer(:,:)
    INTEGER,   INTENT(in)           :: p_source, p_tag
    INTEGER,   INTENT(in), OPTIONAL :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm, l

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       l = p_count
    ELSE
       l = SIZE(buffer)
    END IF

    CALL MPI_IRECV (buffer, l, p_real_dp, p_source, p_tag, &
                    p_comm, p_request(p_irequest), p_error)
    p_irequest = p_irequest + 1

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_2d

  SUBROUTINE p_irecv_int (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout)        :: buffer
    INTEGER, INTENT(in)           :: p_source, p_tag
    INTEGER, INTENT(in), OPTIONAL :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm, l

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       l = p_count
    ELSE
       l = 1
    END IF

    CALL MPI_IRECV (buffer, l, p_int, p_source, p_tag,    &
                    p_comm, p_request(p_irequest), p_error)
    p_irequest = p_irequest + 1

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int

  ! sendrecv implementation

  SUBROUTINE p_sendrecv_real_2d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), P_REAL_DP, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), P_REAL_DP, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_2d

  SUBROUTINE p_sendrecv_real_3d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), P_REAL_DP, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), P_REAL_DP, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_3d

  SUBROUTINE p_sendrecv_real_4d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), P_REAL_DP, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), P_REAL_DP, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_4d

!------------------------------------------------------------------------------
! broadcast character string arrays
!----------------------------------
  SUBROUTINE p_bcast_char_1d (buffer, p_source, comm)
  CHARACTER (*),     INTENT(inout) :: buffer(:)
  INTEGER,           INTENT(in)    :: p_source
  INTEGER, OPTIONAL, INTENT(in)    :: comm
    INTEGER :: i
    DO i=1,SIZE(buffer)
      CALL p_bcast (buffer(i), p_source, comm)
    END DO
  END SUBROUTINE p_bcast_char_1d

!------------------------------------------------------------------------------
! broadcast pointers (including allocation status)
!-------------------------------------------------
  SUBROUTINE p_bcast_ptr_int_1d (buffer, source, comm)
  INTEGER, POINTER :: buffer(:)
  INTEGER,           INTENT(in)    :: source
  INTEGER, OPTIONAL, INTENT(in)    :: comm
    INTEGER :: lu(3)
    IF (p_pe == source) THEN
      lu = (/0,-1,0/)
      IF (ASSOCIATED (buffer)) THEN
        lu (3)   = 1
        lu (1:1) = LBOUND (buffer)
        lu (2:2) = UBOUND (buffer)
      ENDIF
    ENDIF
    CALL p_bcast (lu, source, comm)
    IF (p_pe /= source) THEN
      IF (ASSOCIATED (buffer)) DEALLOCATE (buffer)
      IF (lu(3)==1) ALLOCATE (buffer (lu(1):lu(2)))
      IF(ASSOCIATED(buffer)) buffer=0               ! IBM xlf bug? workaround
    ENDIF
    IF (lu(2) >= lu(1)) CALL p_bcast (buffer, source, comm)
  END SUBROUTINE p_bcast_ptr_int_1d

  SUBROUTINE p_bcast_ptr_int_2d (buffer, source, comm)
  INTEGER, POINTER :: buffer(:,:)
  INTEGER,           INTENT(in)    :: source
  INTEGER, OPTIONAL, INTENT(in)    :: comm
    INTEGER :: lu(5)
    IF (p_pe == source) THEN
      lu = (/0,0,-1,-1,0/)
      IF (ASSOCIATED (buffer)) THEN
        lu (5)   = 1
        lu (1:2) = LBOUND (buffer)
        lu (3:4) = UBOUND (buffer)
      ENDIF
    ENDIF
    CALL p_bcast (lu, source, comm)
    IF (p_pe /= source) THEN
      IF (ASSOCIATED (buffer)) DEALLOCATE (buffer)
      IF (lu(5)==1) ALLOCATE (buffer (lu(1):lu(3),lu(2):lu(4)))
      IF(ASSOCIATED(buffer)) buffer=0               ! IBM xlf bug? workaround
    ENDIF
    IF (lu(3) >= lu(1) .and. lu(4) >= lu(2)) CALL p_bcast (buffer, source, comm)
  END SUBROUTINE p_bcast_ptr_int_2d

#ifdef HAVE_I1

  SUBROUTINE p_bcast_ptr_int1_1d (buffer, source, comm)
  INTEGER(i1),       POINTER    :: buffer(:)
  INTEGER,           INTENT(in) :: source
  INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER :: lu(3)
    IF (p_pe == source) THEN
      lu = (/0,-1,0/)
      IF (ASSOCIATED (buffer)) THEN
        lu (3)   = 1
        lu (1:1) = LBOUND (buffer)
        lu (2:2) = UBOUND (buffer)
      ENDIF
    ENDIF
    CALL p_bcast (lu, source, comm)
    IF (p_pe /= source) THEN
      IF (ASSOCIATED (buffer)) DEALLOCATE (buffer)
      IF (lu(3)==1) ALLOCATE (buffer (lu(1):lu(2)))
      IF(ASSOCIATED(buffer)) buffer=0               ! IBM xlf bug? workaround
    ENDIF
    IF (lu(2) >= lu(1)) CALL p_bcast (buffer, source, comm)
  END SUBROUTINE p_bcast_ptr_int1_1d

  SUBROUTINE p_bcast_ptr_int1_2d (buffer, source, comm)
  INTEGER(i1),       POINTER    :: buffer(:,:)
  INTEGER,           INTENT(in) :: source
  INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER :: lu(5)
    IF (p_pe == source) THEN
      lu = (/0,0,-1,-1,0/)
      IF (ASSOCIATED (buffer)) THEN
        lu (5)   = 1
        lu (1:2) = LBOUND (buffer)
        lu (3:4) = UBOUND (buffer)
      ENDIF
    ENDIF
    CALL p_bcast (lu, source, comm)
    IF (p_pe /= source) THEN
      IF (ASSOCIATED (buffer)) DEALLOCATE (buffer)
      IF (lu(5)==1) ALLOCATE (buffer (lu(1):lu(3),lu(2):lu(4)))
      IF(ASSOCIATED(buffer)) buffer=0               ! IBM xlf bug? workaround
    ENDIF
    IF (lu(3) >= lu(1) .and. lu(4) >= lu(2)) CALL p_bcast (buffer, source, comm)
  END SUBROUTINE p_bcast_ptr_int1_2d

#endif

  SUBROUTINE p_bcast_ptr_real_1d (buffer, source, comm)
  REAL(dp),POINTER :: buffer(:)
  INTEGER,           INTENT(in)    :: source
  INTEGER, OPTIONAL, INTENT(in)    :: comm
    INTEGER :: lu(3)
    IF (p_pe == source) THEN
      lu = (/0,-1,0/)
      IF (ASSOCIATED (buffer)) THEN
        lu (3)   = 1
        lu (1:1) = LBOUND (buffer)
        lu (2:2) = UBOUND (buffer)
      ENDIF
    ENDIF
    CALL p_bcast (lu, source, comm)
    IF (p_pe /= source) THEN
      IF (ASSOCIATED (buffer)) DEALLOCATE (buffer)
      IF (lu(3)==1) ALLOCATE (buffer (lu(1):lu(2)))
      IF(ASSOCIATED(buffer)) buffer=0               ! IBM xlf bug? workaround
    ENDIF
    IF (lu(2) >= lu(1)) CALL p_bcast (buffer, source, comm)
  END SUBROUTINE p_bcast_ptr_real_1d

  SUBROUTINE p_bcast_ptr_real_2d (buffer, source, comm)
  REAL(dp), POINTER :: buffer(:,:)
  INTEGER,           INTENT(in)    :: source
  INTEGER, OPTIONAL, INTENT(in)    :: comm
    INTEGER :: lu(5)
    IF (p_pe == source) THEN
      lu = (/0,0,-1,-1,0/)
      IF (ASSOCIATED (buffer)) THEN
        lu (5)   = 1
        lu (1:2) = LBOUND (buffer)
        lu (3:4) = UBOUND (buffer)
      ENDIF
    ENDIF
    CALL p_bcast (lu, source, comm)
    IF (p_pe /= source) THEN
      IF (ASSOCIATED (buffer)) DEALLOCATE (buffer)
      IF (lu(5)==1) ALLOCATE (buffer (lu(1):lu(3),lu(2):lu(4)))
      IF(ASSOCIATED(buffer)) buffer=0               ! IBM xlf bug? workaround
    ENDIF
    IF (lu(3) >= lu(1) .and. lu(4) >= lu(2)) CALL p_bcast (buffer, source, comm)
  END SUBROUTINE p_bcast_ptr_real_2d

  SUBROUTINE p_bcast_ptr_sp_1d (buffer, source, comm)
  REAL(sp),POINTER :: buffer(:)
  INTEGER,           INTENT(in)    :: source
  INTEGER, OPTIONAL, INTENT(in)    :: comm
    INTEGER :: lu(3)
    IF (p_pe == source) THEN
      lu = (/0,-1,0/)
      IF (ASSOCIATED (buffer)) THEN
        lu (3)   = 1
        lu (1:1) = LBOUND (buffer)
        lu (2:2) = UBOUND (buffer)
      ENDIF
    ENDIF
    CALL p_bcast (lu, source, comm)
    IF (p_pe /= source) THEN
      IF (ASSOCIATED (buffer)) DEALLOCATE (buffer)
      IF (lu(3)==1) ALLOCATE (buffer (lu(1):lu(2)))
      IF(ASSOCIATED(buffer)) buffer=0               ! IBM xlf bug? workaround
    ENDIF
    IF (lu(2) >= lu(1)) CALL p_bcast (buffer, source, comm)
  END SUBROUTINE p_bcast_ptr_sp_1d

  SUBROUTINE p_bcast_ptr_sp_2d (buffer, source, comm)
  REAL(sp), POINTER :: buffer(:,:)
  INTEGER,           INTENT(in)    :: source
  INTEGER, OPTIONAL, INTENT(in)    :: comm
    INTEGER :: lu(5)
    IF (p_pe == source) THEN
      lu = (/0,0,-1,-1,0/)
      IF (ASSOCIATED (buffer)) THEN
        lu (5)   = 1
        lu (1:2) = LBOUND (buffer)
        lu (3:4) = UBOUND (buffer)
      ENDIF
    ENDIF
    CALL p_bcast (lu, source, comm)
    IF (p_pe /= source) THEN
      IF (ASSOCIATED (buffer)) DEALLOCATE (buffer)
      IF (lu(5)==1) ALLOCATE (buffer (lu(1):lu(3),lu(2):lu(4)))
      IF(ASSOCIATED(buffer)) buffer=0               ! IBM xlf bug? workaround
    ENDIF
    IF (lu(3) >= lu(1) .and. lu(4) >= lu(2)) CALL p_bcast (buffer, source, comm)
  END SUBROUTINE p_bcast_ptr_sp_2d


  ! probe implementation

  SUBROUTINE p_probe_real (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    REAL (dp), INTENT(in)  :: buffer
    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i = 1, p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, P_REAL_DP, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO

#endif

  END SUBROUTINE p_probe_real

  SUBROUTINE p_probe_int (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    INTEGER,   INTENT(in)  :: buffer
    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i=1,p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_int, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO


#endif

  END SUBROUTINE p_probe_int

  SUBROUTINE p_probe_bool (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    LOGICAL,   INTENT(in)  :: buffer
    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i=1,p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_bool, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO


#endif

  END SUBROUTINE p_probe_bool

  SUBROUTINE p_probe_char (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    CHARACTER (*),     INTENT(in)  :: buffer
    INTEGER,           INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,           INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i=1,p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_char, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO


#endif

  END SUBROUTINE p_probe_char

!------------------------------------------------------------------------------

  SUBROUTINE p_wait ()
  !---------------------------------------------
  ! waits for an MPI send or receive to complete
  !---------------------------------------------
#ifndef NOMPI
    INTEGER :: i
    DO i = 1, p_irequest-1
       CALL MPI_WAIT (p_request(i), p_status, p_error)
    END DO
    p_irequest = 1
#endif
  END SUBROUTINE p_wait

!------------------------------------------------------------------------------

  SUBROUTINE p_waitall ()
  !---------------------------------------------
  ! waits for all given MPI requests to complete
  !---------------------------------------------
#ifndef NOMPI
    integer :: p_statuses(MPI_STATUS_SIZE,p_irequest-1)

    if (p_irequest > 1) then
       CALL MPI_WAITALL (p_irequest-1, p_request, p_statuses, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_WAITALL failed.'
          WRITE (nerr,'(a,i0)') ' Error =  ', p_error
          CALL p_abort
       ENDIF
    end if
    p_irequest = 1
#endif
  END SUBROUTINE p_waitall

!------------------------------------------------------------------------------

  SUBROUTINE p_barrier (comm)
  !--------------------------------------------------
  ! MPI barrier routine:
  ! wait until all processors have reached this point
  !--------------------------------------------------
  INTEGER ,INTENT(IN) ,OPTIONAL :: comm  ! communicator
#ifndef NOMPI
    INTEGER :: com
    IF(.NOT. p_parallel) RETURN
    com = d_comm; IF(PRESENT(comm)) com = comm
    CALL MPI_BARRIER (com, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BARRIER on ', mype, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif
  END SUBROUTINE p_barrier

!------------------------------------------------------------------------------

  FUNCTION p_or_0d (zfield, comm) RESULT (p_or)

    LOGICAL                       :: p_or
    LOGICAL,           INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_or, 1, p_bool, &
            MPI_LOR, p_comm, p_error)
    ELSE
       p_or = zfield
    END IF
#else
    p_or = zfield
#endif

  END FUNCTION p_or_0d

!------------------------------------------------------------------------------

  FUNCTION p_or_i4_1d (zfield, comm) RESULT (p_or)

    integer(i4),       INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    integer(i4)                   :: p_or (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_or, SIZE(zfield), p_int_i4, &
            MPI_BOR, p_comm, p_error)
    ELSE
       p_or = zfield
    END IF
#else
    p_or = zfield
#endif

  END FUNCTION p_or_i4_1d

!------------------------------------------------------------------------------

  FUNCTION p_or_1d (zfield, comm) RESULT (p_or)

    LOGICAL,           INTENT(in) :: zfield (:)
    LOGICAL                       :: p_or (size(zfield))
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_or, size(zfield), p_bool, &
            MPI_LOR, p_comm, p_error)
    ELSE
       p_or = zfield
    END IF
#else
    p_or = zfield
#endif

  END FUNCTION p_or_1d

!------------------------------------------------------------------------------

  FUNCTION p_ior_0d_i8 (zfield, comm) RESULT (p_ior)

    integer(i8)                   :: p_ior
    INTEGER(i8),       INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_ior, 1, p_int_i8, &
            MPI_BOR, p_comm, p_error)
    ELSE
       p_ior = zfield
    END IF
#else
    p_ior = zfield
#endif

  END FUNCTION p_ior_0d_i8


  FUNCTION p_ior_0d_i4 (zfield, comm) RESULT (p_ior)

    integer(i4)                   :: p_ior
    INTEGER(i4),       INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_ior, 1, p_int_i4, &
            MPI_BOR, p_comm, p_error)
    ELSE
       p_ior = zfield
    END IF
#else
    p_ior = zfield
#endif

  END FUNCTION p_ior_0d_i4
  
  FUNCTION p_ior_1d_i8 (zfield, comm) RESULT (p_ior)
    INTEGER(i8),       INTENT(in) :: zfield(:)
    integer(i8)                   :: p_ior(size(zfield))
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_ior, size(zfield), p_int_i8, &
            MPI_BOR, p_comm, p_error)
    ELSE
       p_ior = zfield
    END IF
#else
    p_ior = zfield
#endif

  END FUNCTION p_ior_1d_i8
  
  FUNCTION p_ior_1d_i4 (zfield, comm) RESULT (p_ior)
    INTEGER(i4),       INTENT(in) :: zfield(:)
    integer(i4)                   :: p_ior(size(zfield))
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_ior, size(zfield), p_int_i4, &
            MPI_BOR, p_comm, p_error)
    ELSE
       p_ior = zfield
    END IF
#else
    p_ior = zfield
#endif

  END FUNCTION p_ior_1d_i4 

!------------------------------------------------------------------------------

  FUNCTION p_and_0d (zfield, comm) RESULT (p_and)

    LOGICAL                       :: p_and
    LOGICAL,           INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_and, 1, p_bool, &
            MPI_LAND, p_comm, p_error)
    ELSE
       p_and = zfield
    END IF
#else
    p_and = zfield
#endif

  END FUNCTION p_and_0d

  FUNCTION p_and_1d (zfield, comm) RESULT (p_and)

    LOGICAL,           INTENT(in) :: zfield (:)
    LOGICAL                       :: p_and (size(zfield))
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_and, size(zfield), p_bool, &
            MPI_LAND, p_comm, p_error)
    ELSE
       p_and = zfield
    END IF
#else
    p_and = zfield
#endif

  END FUNCTION p_and_1d
!==============================================================================
  FUNCTION p_max_i0d (zfield, comm) RESULT (p_max)

    INTEGER                          :: p_max
    INTEGER,              INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_max, 1, MPI_INTEGER, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_i0d
!------------------------------------------------------------------------------
  FUNCTION p_max_int_i8 (zfield, comm) RESULT (p_max)

    INTEGER(i8)                   :: p_max
    INTEGER(i8),       INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_max, 1, p_int_i8, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_int_i8
!------------------------------------------------------------------------------
  FUNCTION p_max_0d (zfield, comm) RESULT (p_max)

    REAL(dp)                          :: p_max
    REAL(dp),              INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_max, 1, P_REAL_DP, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_0d
!------------------------------------------------------------------------------
  FUNCTION p_max_1d (zfield, comm) RESULT (p_max)

    REAL(dp),              INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                          :: p_max (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), P_REAL_DP, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_1d
!------------------------------------------------------------------------------
  FUNCTION p_max_1s (zfield, comm) RESULT (p_max)

    REAL(sp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(sp)                      :: p_max (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), P_REAL_SP, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_1s
!------------------------------------------------------------------------------
  FUNCTION p_max_2d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), P_REAL_DP, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_2d
!------------------------------------------------------------------------------
  FUNCTION p_max_i1d (zfield, comm) RESULT (p_max)

    INTEGER,           INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_max (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), MPI_INTEGER, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_i1d
!------------------------------------------------------------------------------
  FUNCTION p_max_i2d (zfield, comm) RESULT (p_max)

    INTEGER,           INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_max (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), MPI_INTEGER, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_i2d
!------------------------------------------------------------------------------
  FUNCTION p_max_3d (zfield, comm) RESULT (p_max)

    REAL(dp),              INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                          :: p_max (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), P_REAL_DP, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_3d

  FUNCTION p_max_4d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3),SIZE(zfield,4))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN

       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), P_REAL_DP, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_4d


  FUNCTION p_sum_0d (zfield, comm) RESULT (p_sum)

    REAL(dp)                          :: p_sum
    REAL(dp),              INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, P_REAL_DP, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif
  END FUNCTION p_sum_0d

  FUNCTION p_sum_int_i8 (zfield, comm) RESULT (p_sum)

    INTEGER(i8)                   :: p_sum
    INTEGER(i8),       INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, p_int_i8, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_int_i8

  FUNCTION p_sum_int_0d (zfield, comm) RESULT (p_sum)

    INTEGER                           :: p_sum
    INTEGER,               INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, P_INT, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_int_0d

  FUNCTION p_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),              INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                          :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), P_REAL_DP, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_1d

  FUNCTION p_sum_int_1d (zfield, comm) RESULT (p_sum)

    INTEGER,           INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), P_INT, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_int_1d

  FUNCTION p_sum_int_2d (zfield, comm) RESULT (p_sum)

    INTEGER,           INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), P_INT, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_int_2d

  FUNCTION p_sum_2d (zfield, comm) RESULT (p_sum)

    REAL(dp),              INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                          :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), P_REAL_DP, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_2d

  FUNCTION p_sum_3d (zfield, comm) RESULT (p_sum)

    REAL(dp),              INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                          :: p_sum (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), P_REAL_DP, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_3d

  FUNCTION p_sum_int_3d (zfield, comm) RESULT (p_sum)

    INTEGER,               INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                           :: p_sum (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), P_INT, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_int_3d

  FUNCTION p_sum_int_4d (zfield, comm) RESULT (p_sum)

    INTEGER,           INTENT(in) :: zfield(:,:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_sum (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3),SIZE(zfield,4))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), P_INT, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_int_4d

  FUNCTION p_sum_4d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3),SIZE(zfield,4))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN

       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), P_REAL_DP, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_4d

  FUNCTION p_sum_i8_1d (zfield, comm) RESULT (p_sum)

    INTEGER(i8),       INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER(i8)                   :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_int_i8, &
                           MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_i8_1d


!------------------------------------------------------------------------------

#ifdef HAVE_I1
  FUNCTION p_sum_i1_1d (zfield, comm) RESULT (p_sum)

    integer(i1),       INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    integer(i1)                   :: p_sum (SIZE(zfield))

#ifndef NOMPI

    integer(i2) :: tmp (size(zfield))
    IF(p_parallel)THEN
      tmp   = zfield
      p_sum = p_sum_i2_1d (tmp, comm)
    ELSE
      p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_i1_1d
#endif
!------------------------------------------------------------------------------

  FUNCTION p_sum_i2_1d (zfield, comm) RESULT (p_sum)

    integer(i2),       INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    integer(i2)                   :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_int_i2, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_i2_1d

!------------------------------------------------------------------------------

  FUNCTION p_min_int_0d (zfield, comm) RESULT (p_min)

    INTEGER,           INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_min
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_min, 1, P_INT, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_int_0d
!------------------------------------------------------------------------------
  FUNCTION p_min_int_1d (zfield, comm) RESULT (p_min)

    INTEGER,           INTENT(in) :: zfield (:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_min  (SIZE(zfield))
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), P_INT, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_int_1d
!------------------------------------------------------------------------------
  FUNCTION p_min_int_i8 (zfield, comm) RESULT (p_min)

    INTEGER(i8)                   :: p_min
    INTEGER(i8),       INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_min, 1, p_int_i8, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_int_i8
!------------------------------------------------------------------------------
#ifdef HAVE_I1
  FUNCTION p_min_i1_1d (zfield, comm) RESULT (p_min)

    integer(i1),       INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    integer(i1)                   :: p_min (SIZE(zfield))

#ifndef NOMPI

!***   An error occurred in MPI_Allreduce:
!***   the reduction operation MPI_MIN is not defined on the MPI_BYTE datatype
!thus: use integer(i2)

    integer(i2) :: tmp (size(zfield))
    IF(p_parallel)THEN
      tmp   = zfield
      p_min = p_min_i2_1d (tmp, comm)
    ELSE
      p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_i1_1d
#endif
!------------------------------------------------------------------------------
  FUNCTION p_min_i2_1d (zfield, comm) RESULT (p_min)

    integer(i2),       INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    integer(i2)                   :: p_min (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_int_i2, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_i2_1d
!------------------------------------------------------------------------------
  FUNCTION p_min_0d (zfield, comm) RESULT (p_min)

    REAL(dp),              INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                          :: p_min
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_min, 1, P_REAL_DP, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_0d

  FUNCTION p_min_1d (zfield, comm) RESULT (p_min)

    REAL(dp),              INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                          :: p_min (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), P_REAL_DP, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_1d

  FUNCTION p_min_2d (zfield, comm) RESULT (p_min)

    REAL(dp),              INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                          :: p_min (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), P_REAL_DP, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_2d

  FUNCTION p_min_3d (zfield, comm) RESULT (p_min)

    REAL(dp),              INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                          :: p_min (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), P_REAL_DP, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_3d

  FUNCTION p_min_4d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3),SIZE(zfield,4))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = d_comm
    ENDIF

!    IF (P_REAL_DP /= MPI_DATATYPE_NULL) THEN
IF(p_parallel)THEN

       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), P_REAL_DP, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_4d

!------------------------------------------------------------------------------
  SUBROUTINE p_gather_bool (sendbuf, recvbuf, root, comm)
  LOGICAL ,INTENT(in)           :: sendbuf     ! send buffer
  LOGICAL ,INTENT(out)          :: recvbuf (:) ! receive buffer
  INTEGER ,INTENT(in)           :: root        ! rank of receiving process
  INTEGER ,INTENT(in) ,OPTIONAL :: comm        ! communicator handle

    INTEGER :: lcom, npes

#ifndef NOMPI

    lcom = d_comm; IF(PRESENT(comm)) lcom = comm

    CALL MPI_COMM_SIZE (lcom, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

IF(p_parallel)THEN
    IF (size(recvbuf) /= npes) THEN
       WRITE (nerr,'(a)') ' p_gather_bool: size(recvbuf) /= size(communicator)'
       CALL p_abort
    ENDIF

    CALL MPI_GATHER (sendbuf, 1, MPI_LOGICAL, &
                     recvbuf, 1, MPI_LOGICAL, &
                     root, lcom, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_GATHER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF
ELSE
    recvbuf = sendbuf
END IF

#else
!   CALL util_exit(1)
    IF (size(recvbuf)/=1) &
      CALL abort('p_gather_bool: size(recvbuf)/=1 in single processor mode')
    recvbuf = sendbuf
#endif

  END SUBROUTINE p_gather_bool
!------------------------------------------------------------------------------
  SUBROUTINE p_gather_int (sendbuf, recvbuf, root, comm)
  INTEGER ,INTENT(in)           :: sendbuf     ! send buffer
  INTEGER ,INTENT(out)          :: recvbuf (:) ! receive buffer
  INTEGER ,INTENT(in)           :: root        ! rank of receiving process
  INTEGER ,INTENT(in) ,OPTIONAL :: comm        ! communicator handle

    INTEGER :: lcom, npes

#ifndef NOMPI

    lcom = d_comm; IF(PRESENT(comm)) lcom = comm

    CALL MPI_COMM_SIZE (lcom, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

    IF (size(recvbuf) /= npes) THEN
       WRITE (nerr,'(a)') ' p_gather_int: size(recvbuf) /= size(communicator)'
       CALL p_abort
    ENDIF

IF(p_parallel)THEN
    CALL MPI_GATHER (sendbuf, 1, MPI_INTEGER, &
                     recvbuf, 1, MPI_INTEGER, &
                     root, lcom, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_GATHER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF
ELSE
    recvbuf = sendbuf
END IF

#else
!   CALL util_exit(1)
    IF (size(recvbuf)/=1) &
      CALL abort('p_gather_int: size(recvbuf)/=1 in single processor mode')
    recvbuf = sendbuf
#endif

  END SUBROUTINE p_gather_int
!------------------------------------------------------------------------------
  SUBROUTINE p_gather_real (sendbuf, recvbuf, root, comm)
  REAL(dp) ,INTENT(in)           :: sendbuf     ! send buffer
  REAL(dp) ,INTENT(out)          :: recvbuf (:) ! receive buffer
  INTEGER  ,INTENT(in)           :: root        ! rank of receiving process
  INTEGER  ,INTENT(in) ,OPTIONAL :: comm        ! communicator handle

    INTEGER :: lcom, npes

#ifndef NOMPI

    lcom = d_comm; IF(PRESENT(comm)) lcom = comm

    CALL MPI_COMM_SIZE (lcom, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

    IF (size(recvbuf) /= npes) THEN
       WRITE (nerr,'(a)') ' p_gather_real: size(recvbuf) /= size(communicator)'
       CALL p_abort
    ENDIF

IF(p_parallel)THEN
    CALL MPI_GATHER (sendbuf, 1, P_REAL_DP, &
                     recvbuf, 1, P_REAL_DP, &
                     root, lcom, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_GATHER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF
ELSE
    recvbuf = sendbuf
END IF

#else
!   CALL util_exit(1)
    IF (size(recvbuf)/=1) &
      CALL abort('p_gather_real: size(recvbuf)/=1 in single processor mode')
    recvbuf = sendbuf
#endif

  END SUBROUTINE p_gather_real
!------------------------------------------------------------------------------
  SUBROUTINE p_gather_real_sp (sendbuf, recvbuf, root, comm)
  REAL(sp) ,INTENT(in)           :: sendbuf     ! send buffer
  REAL(sp) ,INTENT(out)          :: recvbuf (:) ! receive buffer
  INTEGER  ,INTENT(in)           :: root        ! rank of receiving process
  INTEGER  ,INTENT(in) ,OPTIONAL :: comm        ! communicator handle

    INTEGER :: lcom, npes

#ifndef NOMPI

    lcom = d_comm; IF(PRESENT(comm)) lcom = comm

    CALL MPI_COMM_SIZE (lcom, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

    IF (size(recvbuf) /= npes) THEN
       WRITE (nerr,'(a)') ' p_gather_real_sp: size(recvbuf) /= size(communicator)'
       CALL p_abort
    ENDIF

IF(p_parallel)THEN
    CALL MPI_GATHER (sendbuf, 1, P_REAL_SP, &
                     recvbuf, 1, P_REAL_SP, &
                     root, lcom, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_GATHER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF
ELSE
    recvbuf = sendbuf
END IF

#else
!   CALL util_exit(1)
    IF (size(recvbuf)/=1) &
      CALL abort('p_gather_real_sp: size(recvbuf)/=1 in single processor mode')
    recvbuf = sendbuf
#endif

  END SUBROUTINE p_gather_real_sp
!------------------------------------------------------------------------------
  SUBROUTINE p_allgather_int (sendbuf, recvbuf, comm)
  INTEGER ,INTENT(in)           :: sendbuf     ! send buffer
  INTEGER ,INTENT(out)          :: recvbuf (:) ! receive buffer
  INTEGER ,INTENT(in) ,OPTIONAL :: comm        ! communicator handle

    INTEGER :: lcom, npes

#ifndef NOMPI

    lcom = d_comm; IF(PRESENT(comm)) lcom = comm

    CALL MPI_COMM_SIZE (lcom, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

    IF (size(recvbuf) /= npes) THEN
       WRITE(nerr,'(a)')' p_allgather_int: size(recvbuf) /= size(communicator)'
       CALL p_abort
    ENDIF

    CALL MPI_ALLGATHER (sendbuf, 1, MPI_INTEGER, &
                        recvbuf, 1, MPI_INTEGER, &
                        lcom, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ALLGATHER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

#else
!   CALL util_exit(1)
    IF (size(recvbuf)/=1) &
      CALL abort('p_allgather_int: size(recvbuf)/=1 in single processor mode')
    recvbuf = sendbuf
#endif

  END SUBROUTINE p_allgather_int
!------------------------------------------------------------------------------
  SUBROUTINE p_scatterv_real_1d (sendbuf, sendcounts, recvbuf, root, comm)
  real(dp) ,INTENT(in)           :: sendbuf    (:)  ! send    buffer
  INTEGER  ,INTENT(in)           :: sendcounts(0:)  ! no. of elements to send
  real(dp) ,INTENT(out)          :: recvbuf    (:)  ! receive buffer
  INTEGER  ,INTENT(in)           :: root            ! rank of receiving process
  INTEGER  ,INTENT(in) ,OPTIONAL :: comm            ! communicator handle

    INTEGER :: lcom, npes
    INTEGER :: recvcnt                  ! Receive count for MPI_SCATTERV

#ifndef NOMPI
    integer              :: pe
    integer, allocatable :: displs(:)   ! Displacements

    lcom = d_comm; IF(PRESENT(comm)) lcom = comm

    if (.NOT. p_parallel) then
      if (size(sendbuf) /= size(recvbuf)) then
        WRITE (nerr,'(a)') &
        ' p_scatterv_real_1d: npes==1 .and. size(sendbuf)/=size(recvbuf)'
        CALL p_abort
      endif
      recvbuf = sendbuf
      return
    end if

    CALL MPI_COMM_SIZE (lcom, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

    if (npes==1) then
      if (size(sendbuf) /= size(recvbuf)) then
        WRITE (nerr,'(a)') &
        ' p_scatterv_real_1d: npes==1 .and. size(sendbuf)/=size(recvbuf)'
        CALL p_abort
      endif
      recvbuf = sendbuf
      return
    endif

    allocate (displs(0:npes))
    displs(0) = 0
    do pe = 1, npes
       displs(pe) = displs(pe-1) + sendcounts(pe-1)
    end do

!if (p_pe == root) then
!   write(0,*) "p_scatterv_real_1d: npes         =", npes
!   write(0,*) "p_scatterv_real_1d: sendcounts   =", sendcounts
!   write(0,*) "p_scatterv_real_1d: displs       =", displs
!   write(0,*) "p_scatterv_real_1d: size(sendbuf)=", size(sendbuf)
!   write(0,*) "p_scatterv_real_1d: size(recvbuf)=", size(recvbuf)
!end if
    !----------------------------------------
    ! Validate size of sendbuf for the sender
    !----------------------------------------
    if (p_pe == root) then
       IF (size(sendbuf) /= displs(npes)) THEN
          WRITE (nerr,'(a,i6,a,2i11)') ' p_scatterv_real_1d: pe=', p_pe, &
               'size(sendbuf) /= sum(sendcounts)', size (recvbuf), displs(npes)
          CALL p_abort
       ENDIF
    end if
    recvcnt = size (recvbuf)
!   recvcnt = sendcounts(p_pe)
!   IF (size (recvbuf) /= recvcnt) THEN
!      WRITE (nerr,'(a,i6,a,2i11)') ' p_scatterv_real_1d: pe=', p_pe, &
!           'size(recvbuf) /= recvcount', size (recvbuf), recvcnt
!      CALL p_abort
!   ENDIF

    CALL MPI_SCATTERV (sendbuf, sendcounts, displs, P_REAL_DP, &
                       recvbuf, recvcnt,            P_REAL_DP, &
                       root, lcom, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_SCATTERV failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

#else
    recvcnt = sendcounts(p_pe)
    IF (size (recvbuf) /= size (sendbuf) .or. size (sendbuf) /= recvcnt) THEN
       write (nerr,'(a,i11)') 'size(sendbuf) =', size (sendbuf)
       write (nerr,'(a,i11)') 'size(recvbuf) =', size (recvbuf)
       write (nerr,'(a,i11)') 'sendcounts(0) =', sendcounts(0)
       CALL abort('p_scatterv_real_1d: in single processor mode: bad args')
    END IF
    recvbuf(1:recvcnt) = sendbuf(1:recvcnt)
#endif

  END SUBROUTINE p_scatterv_real_1d
!------------------------------------------------------------------------------
  SUBROUTINE p_scatterv_int_1d (sendbuf, sendcounts, recvbuf, root, comm)
  INTEGER  ,INTENT(in)           :: sendbuf    (:)  ! send    buffer
  INTEGER  ,INTENT(in)           :: sendcounts(0:)  ! no. of elements to send
  INTEGER  ,INTENT(out)          :: recvbuf    (:)  ! receive buffer
  INTEGER  ,INTENT(in)           :: root            ! rank of receiving process
  INTEGER  ,INTENT(in) ,OPTIONAL :: comm            ! communicator handle

    INTEGER :: lcom, npes
    INTEGER :: recvcnt                  ! Receive count for MPI_SCATTERV

#ifndef NOMPI
    integer              :: pe
    integer, allocatable :: displs(:)   ! Displacements

    lcom = d_comm; IF(PRESENT(comm)) lcom = comm

    if (.NOT. p_parallel) then
      if (size(sendbuf) /= size(recvbuf)) then
        WRITE (nerr,'(a)') &
        ' p_scatterv_int_1d: npes==1 .and. size(sendbuf)/=size(recvbuf)'
        CALL p_abort
      endif
      recvbuf = sendbuf
      return
    end if

    CALL MPI_COMM_SIZE (lcom, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

    if (npes==1) then
      if (size(sendbuf) /= size(recvbuf)) then
        WRITE (nerr,'(a)') &
        ' p_scatterv_int_1d: npes==1 .and. size(sendbuf)/=size(recvbuf)'
        CALL p_abort
      endif
      recvbuf = sendbuf
      return
    endif

    allocate (displs(0:npes))
    displs(0) = 0
    do pe = 1, npes
       displs(pe) = displs(pe-1) + sendcounts(pe-1)
    end do

!if (p_pe == root) then
!   write(0,*) "p_scatterv_int_1d: npes         =", npes
!   write(0,*) "p_scatterv_int_1d: sendcounts   =", sendcounts
!   write(0,*) "p_scatterv_int_1d: displs       =", displs
!   write(0,*) "p_scatterv_int_1d: size(sendbuf)=", size(sendbuf)
!   write(0,*) "p_scatterv_int_1d: size(recvbuf)=", size(recvbuf)
!end if
    !----------------------------------------
    ! Validate size of sendbuf for the sender
    !----------------------------------------
    if (p_pe == root) then
       IF (size(sendbuf) /= displs(npes)) THEN
          WRITE (nerr,'(a,i6,a,2i11)') ' p_scatterv_int_1d: pe=', p_pe, &
               'size(sendbuf) /= sum(sendcounts)', size (recvbuf), displs(npes)
          CALL p_abort
       ENDIF
    end if
    recvcnt = size (recvbuf)
!   recvcnt = sendcounts(p_pe)
!   IF (size (recvbuf) /= recvcnt) THEN
!      WRITE (nerr,'(a,i6,a,2i11)') ' p_scatterv_int_1d: pe=', p_pe, &
!           'size(recvbuf) /= recvcount', size (recvbuf), recvcnt
!      CALL p_abort
!   ENDIF

    CALL MPI_SCATTERV (sendbuf, sendcounts, displs, P_INT, &
                       recvbuf, recvcnt,            P_INT, &
                       root, lcom, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_SCATTERV failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ENDIF

#else
    recvcnt = sendcounts(p_pe)
    IF (size (recvbuf) /= size (sendbuf) .or. size (sendbuf) /= recvcnt) THEN
       write (nerr,'(a,i11)') 'size(sendbuf) =', size (sendbuf)
       write (nerr,'(a,i11)') 'size(recvbuf) =', size (recvbuf)
       write (nerr,'(a,i11)') 'sendcounts(0) =', sendcounts(0)
       CALL abort('p_scatterv_int_1d: in single processor mode: bad args')
    END IF
    recvbuf(1:recvcnt) = sendbuf(1:recvcnt)
#endif

  END SUBROUTINE p_scatterv_int_1d

!==============================================================================
! specific p_bcast routines
!--------------------------
#undef  VECTOR
#define DERIVED INTEGER(i4)
#define p_bcast_DERIVED p_bcast_int_i4
#define MPI_TYPE p_int_i4
#include "p_bcast.incf"
#undef  DERIVED
#undef  p_bcast_DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define DERIVED INTEGER(i8)
#define p_bcast_DERIVED p_bcast_int_i8
#define MPI_TYPE p_int_i8
#include "p_bcast.incf"
#undef  DERIVED
#undef  p_bcast_DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define VECTOR
#define p_bcast_DERIVED p_bcast_int_1d
#define DERIVED INTEGER,DIMENSION(:)
#define MPI_TYPE p_int
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_int_2d
#define DERIVED INTEGER,DIMENSION(:,:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_int_3d
#define DERIVED INTEGER,DIMENSION(:,:,:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_int_4d
#define DERIVED INTEGER,DIMENSION(:,:,:,:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
#undef  MPI_TYPE
!==============================================================================
#define p_bcast_DERIVED p_bcast_int_1d_i2
#define DERIVED INTEGER(i2),DIMENSION(:)
#define MPI_TYPE p_int_i2
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
#undef  MPI_TYPE
!==============================================================================
#ifdef HAVE_I1
#define p_bcast_DERIVED p_bcast_int_1d_i1
#define DERIVED INTEGER(i1),DIMENSION(:)
#define MPI_TYPE MPI_BYTE
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_int_2d_i1
#define DERIVED INTEGER(i1),DIMENSION(:,:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
#undef  MPI_TYPE
#endif
!==============================================================================
#undef  VECTOR
#define MPI_TYPE p_real_wp
#define DERIVED REAL(wp)
#define p_bcast_DERIVED p_bcast_real
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define VECTOR
#define p_bcast_DERIVED p_bcast_real_3d
#define DERIVED REAL(wp),DIMENSION(:,:,:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_real_4d
#define DERIVED REAL(wp),DIMENSION(:,:,:,:)
#include "p_bcast.incf"
#undef  MPI_TYPE
#undef  p_bcast_DERIVED
#undef  DERIVED
!==============================================================================
#define MPI_TYPE p_real_dp
#define p_bcast_DERIVED p_bcast_real_1d
#define DERIVED REAL(dp),DIMENSION(:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_real_2d
#define DERIVED REAL(dp),DIMENSION(:,:)
#include "p_bcast.incf"
#undef  MPI_TYPE
#undef  p_bcast_DERIVED
#undef  DERIVED
!==============================================================================
#undef  VECTOR
#define MPI_TYPE p_real_sp
#define DERIVED REAL(sp)
#define p_bcast_DERIVED p_bcast_real_0d_sp
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define VECTOR
#define p_bcast_DERIVED p_bcast_real_1d_sp
#define DERIVED REAL(sp),DIMENSION(:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_real_2d_sp
#define DERIVED REAL(sp),DIMENSION(:,:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_real_3d_sp
#define DERIVED REAL(sp),DIMENSION(:,:,:)
#include "p_bcast.incf"
#undef  DERIVED
#undef  p_bcast_DERIVED
#undef  MPI_TYPE
!==============================================================================
#undef  VECTOR
#define DERIVED LOGICAL
#define p_bcast_DERIVED p_bcast_bool
#define MPI_TYPE p_bool
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define VECTOR
#define p_bcast_DERIVED p_bcast_bool_1d
#define DERIVED LOGICAL,DIMENSION(:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_bool_2d
#define DERIVED LOGICAL,DIMENSION(:,:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_bool_3d
#define DERIVED LOGICAL,DIMENSION(:,:,:)
#include "p_bcast.incf"
#undef  p_bcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_bcast_DERIVED p_bcast_bool_4d
#define DERIVED LOGICAL,DIMENSION(:,:,:,:)
#include "p_bcast.incf"
#undef  DERIVED
#undef  p_bcast_DERIVED
#undef  MPI_TYPE
!==============================================================================
#undef  VECTOR
#define DERIVED CHARACTER(len=*)
#define p_bcast_DERIVED p_bcast_char
#include "p_bcast.incf"
#undef  DERIVED
#undef  p_bcast_DERIVED
#undef  MPI_TYPE
!==============================================================================
! specific p_ibcast routines
!---------------------------
#define VECTOR
#define p_ibcast_DERIVED p_ibcast_int_1d
#define DERIVED INTEGER,DIMENSION(:)
#define MPI_TYPE p_int
#include "p_ibcast.incf"
#undef  p_ibcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_ibcast_DERIVED p_ibcast_int_2d
#define DERIVED INTEGER,DIMENSION(:,:)
#include "p_ibcast.incf"
#undef  p_ibcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_ibcast_DERIVED p_ibcast_int_3d
#define DERIVED INTEGER,DIMENSION(:,:,:)
#include "p_ibcast.incf"
#undef  p_ibcast_DERIVED
#undef  DERIVED
!------------------------------------------------------------------------------
#define p_ibcast_DERIVED p_ibcast_int_4d
#define DERIVED INTEGER,DIMENSION(:,:,:,:)
#include "p_ibcast.incf"
#undef  p_ibcast_DERIVED
#undef  DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define VECTOR
#define DERIVED REAL(wp),DIMENSION(:)
#define p_ibcast_DERIVED p_ibcast_real_1d
#define MPI_TYPE p_real_wp
#include "p_ibcast.incf"
#undef  DERIVED
#undef  p_ibcast_DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define DERIVED REAL(wp),DIMENSION(:,:)
#define p_ibcast_DERIVED p_ibcast_real_2d
#define MPI_TYPE p_real_wp
#include "p_ibcast.incf"
#undef  DERIVED
#undef  p_ibcast_DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define DERIVED REAL(wp),DIMENSION(:,:,:)
#define p_ibcast_DERIVED p_ibcast_real_3d
#define MPI_TYPE p_real_wp
#include "p_ibcast.incf"
#undef  DERIVED
#undef  p_ibcast_DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define DERIVED REAL(wp),DIMENSION(:,:,:,:)
#define p_ibcast_DERIVED p_ibcast_real_4d
#define MPI_TYPE p_real_wp
#include "p_ibcast.incf"
#undef  DERIVED
#undef  p_ibcast_DERIVED
#undef  MPI_TYPE
!==============================================================================
! specific p_alltoall routines
!-----------------------------
#undef  VECTOR
#define DERIVED INTEGER
#define p_alltoall_DERIVED p_alltoall_int
#define MPI_TYPE p_int
#include "p_alltoall_derived.incf"
#undef  DERIVED
#undef  p_alltoall_DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define DERIVED REAL(WP)
#define p_alltoall_DERIVED p_alltoall_real
#define MPI_TYPE p_real_wp
#include "p_alltoall_derived.incf"
#undef  DERIVED
#undef  p_alltoall_DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define DERIVED LOGICAL
#define p_alltoall_DERIVED p_alltoall_bool
#define MPI_TYPE p_bool
#include "p_alltoall_derived.incf"
#undef  DERIVED
#undef  p_alltoall_DERIVED
#undef  MPI_TYPE
!==============================================================================
! specific p_gather(v) routines
!------------------------------
#define DERIVED  REAL(WP)
#define MPI_TYPE p_real_wp
#define p_gather_DERIVED p_gatherv_real_1d
#include "p_gather_derived.incf"
#undef  p_gather_DERIVED
#undef  DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define DERIVED  REAL(SP)
#define MPI_TYPE p_real_sp
#define p_gather_DERIVED p_gatherv_real_1d_sp
#include "p_gather_derived.incf"
#undef  p_gather_DERIVED
#undef  DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define DERIVED  REAL(WP)
#define MPI_TYPE p_real_wp
#define DIMS (:,:)
#define p_gather_DERIVED p_gatherv_real_2d
#include "p_gather_derived.incf"
#undef  p_gather_DERIVED
#undef  DIMS
#undef  DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define DERIVED  REAL(WP)
#define MPI_TYPE p_real_wp
#define DIMS (:,:,:)
#define p_gather_DERIVED p_gatherv_real_3d
#include "p_gather_derived.incf"
#undef  p_gather_DERIVED
#undef  DIMS
#undef  DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define MPI_TYPE p_int_i4
#define DERIVED  INTEGER(i4)
#define p_gather_DERIVED p_gatherv_int_1d_i4
#include "p_gather_derived.incf"
#undef  p_gather_DERIVED
#undef  DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define MPI_TYPE p_int_i4
#define DERIVED  INTEGER(i4)
#define DIMS (:,:)
#define p_gather_DERIVED p_gatherv_int_2d_i4
#include "p_gather_derived.incf"
#undef  p_gather_DERIVED
#undef  DIMS
#undef  DERIVED
#undef  MPI_TYPE
!------------------------------------------------------------------------------
#define VECTOR   /* Use as hack for character(*) */
#define DERIVED  CHARACTER(len=*)
#define p_gather_DERIVED p_gatherv_char_1d
#include "p_gather_derived.incf"
#undef  p_gather_DERIVED
#undef  DERIVED
#undef  MPI_TYPE
#undef  VECTOR
!------------------------------------------------------------------------------
#define DERIVED  LOGICAL
#define MPI_TYPE p_bool
#define DIMS (:,:)
#define p_gather_DERIVED p_gatherv_bool_2d
#include "p_gather_derived.incf"
#undef  p_gather_DERIVED
#undef  DIMS
#undef  DERIVED
#undef  MPI_TYPE
!=============================
! specific p_alltoall routines
!-----------------------------
#define VECTOR
#define DERIVED INTEGER
#define p_allgather_DERIVED p_allgatherv_int
#define MPI_TYPE MPI_INTEGER
#include "p_allgather.incf"
#undef  p_allgather_DERIVED
#undef  DERIVED
#undef  MPI_TYPE
#undef  VECTOR
!==============================================================================

!==============================================================================
  SUBROUTINE alltoallv_args (bufsize, sbuflen, rbuflen, groupsize, comm,     &
                             sendcounts, recvcounts,                         &
                             p_sendcounts, p_sdispls, p_recvcounts, p_rdispls)
  !---------------------------------------------
  ! prepare arguments for MPI_ALLTOALLV routines
  !---------------------------------------------
  INTEGER ,INTENT(in)                   :: bufsize
  INTEGER ,INTENT(in)                   :: sbuflen
  INTEGER ,INTENT(in)                   :: rbuflen
  INTEGER ,INTENT(in)                   :: groupsize
  INTEGER ,INTENT(in)                   :: comm
  INTEGER ,INTENT(in)                   :: sendcounts   (:)
  INTEGER ,INTENT(in) ,OPTIONAL         :: recvcounts   (:)
  INTEGER ,POINTER                      :: p_sendcounts (:)
  INTEGER ,POINTER                      :: p_sdispls    (:)
  INTEGER ,POINTER                      :: p_recvcounts (:)
  INTEGER ,POINTER                      :: p_rdispls    (:)

#ifdef NOMPI
    call abort ('alltoallv_args: should never be called ifdef NOMPI')
#else

    INTEGER :: p_error, i

    IF (SUM (sendcounts) /= sbuflen) THEN
      WRITE (nerr,*) p_pe,' alltoallv_args:  sum (sendcounts) /= sbuflen',&
                                             sum (sendcounts) ,  sbuflen
      CALL p_abort
    ENDIF

    ALLOCATE (p_sendcounts (groupsize))
    ALLOCATE (p_recvcounts (groupsize))
    ALLOCATE (p_sdispls    (groupsize))
    ALLOCATE (p_rdispls    (groupsize))

    p_sendcounts = bufsize * sendcounts
    IF (PRESENT (recvcounts)) THEN
      p_recvcounts = recvcounts
    ELSE
      CALL MPI_ALLTOALL (  sendcounts, 1, p_int, &
                         p_recvcounts, 1, p_int, &
                         comm, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_ALLTOALL failed in alltoallv_args.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      ENDIF
    ENDIF
    IF (SUM (p_recvcounts) /= rbuflen) THEN
      WRITE (nerr,'(a)') ' alltoallv_args: sum (recvcounts) /= rbuflen'

!+++++++++++++++++++++++++++++++++
!write some additional diagnostics
!+++++++++++++++++++++++++++++++++
WRITE (nerr,*) '    p_pe               =',p_pe

WRITE (nerr,*) '    bufsize            =',bufsize
WRITE (nerr,*) '    sbuflen            =',sbuflen
WRITE (nerr,*) '    rbuflen            =',rbuflen
WRITE (nerr,*) '    groupsize          =',groupsize
WRITE (nerr,*) '    comm               =',comm
WRITE (nerr,*) '    size(sendcounts)   =',size(sendcounts)
WRITE (nerr,*) '    sendcounts         =',sendcounts

WRITE (nerr,*) '    p_sendcounts       =',p_sendcounts
WRITE (nerr,*) '    p_recvcounts       =',p_recvcounts

WRITE (nerr,*) '    PRESENT(recvcounts)=',PRESENT(recvcounts)
if(PRESENT(recvcounts)) then
WRITE (nerr,*) '    size(recvcounts)   =',size(recvcounts)
WRITE (nerr,*) '    recvcounts         =',recvcounts
endif

      CALL p_abort
    ENDIF
    p_recvcounts = bufsize * p_recvcounts

    p_sdispls (1) = 0
    p_rdispls (1) = 0
    DO i = 2, groupsize
      p_sdispls (i) = p_sdispls (i-1) + p_sendcounts (i-1)
      p_rdispls (i) = p_rdispls (i-1) + p_recvcounts (i-1)
    END DO

#endif

  END SUBROUTINE alltoallv_args
!------------------------------------------------------------------------------
  SUBROUTINE allgatherv_args (bufsize, sbuflen, rbuflen, groupsize, comm,     &
                              recvcounts, p_sendcount, p_recvcounts, p_rdispls)
  !----------------------------------------------
  ! prepare arguments for MPI_ALLGATHERV routines
  !----------------------------------------------
  INTEGER ,INTENT(in)                   :: bufsize
  INTEGER ,INTENT(in)                   :: sbuflen
  INTEGER ,INTENT(in)                   :: rbuflen
  INTEGER ,INTENT(in)                   :: groupsize
  INTEGER ,INTENT(in)                   :: comm
  INTEGER ,INTENT(in) ,OPTIONAL         :: recvcounts   (:)
  INTEGER ,INTENT(out)                  :: p_sendcount
  INTEGER ,POINTER                      :: p_recvcounts (:)
  INTEGER ,POINTER                      :: p_rdispls    (:)

#ifdef NOMPI
    call abort ('allgatherv_args: should never be called ifdef NOMPI')
#else

    INTEGER :: p_error, i

    ALLOCATE (p_recvcounts (groupsize))
    ALLOCATE (p_rdispls    (groupsize))

    p_sendcount = bufsize * sbuflen
    IF (PRESENT (recvcounts)) THEN
      p_recvcounts = recvcounts
    ELSE
      CALL MPI_ALLGATHER (sbuflen,      1, p_int, &
                          p_recvcounts, 1, p_int, &
                          comm, p_error)
    ENDIF
    IF (SUM (p_recvcounts) /= rbuflen) THEN
      WRITE (nerr,'(a)') ' allgatherv_args: sum (recvcounts) /= rbuflen'
      CALL p_abort
    ENDIF
    p_recvcounts = bufsize * p_recvcounts

    p_rdispls (1) = 0
    DO i = 2, groupsize
      p_rdispls (i) = p_rdispls (i-1) + p_recvcounts (i-1)
    END DO

#endif

  END SUBROUTINE allgatherv_args
!------------------------------------------------------------------------------
  SUBROUTINE gatherv_args (bufsize, sbuflen, rbuflen, groupsize, root, comm, &
                           recvcounts, p_sendcount, p_recvcounts, p_rdispls)
  !-------------------------------------------++++++++++++++++
  ! prepare arguments for MPI_GATHERV routines NOT YET TESTED
  !-------------------------------------------++++++++++++++++
  INTEGER ,INTENT(in)                   :: bufsize
  INTEGER ,INTENT(in)                   :: sbuflen
  INTEGER ,INTENT(in)                   :: rbuflen
  INTEGER ,INTENT(in)                   :: groupsize
  INTEGER ,INTENT(in)                   :: root
  INTEGER ,INTENT(in)                   :: comm
  INTEGER ,INTENT(in) ,OPTIONAL         :: recvcounts   (:)
  INTEGER ,INTENT(out)                  :: p_sendcount
  INTEGER ,POINTER                      :: p_recvcounts (:)
  INTEGER ,POINTER                      :: p_rdispls    (:)

#ifdef NOMPI
    call abort ('gatherv_args: should never be called ifdef NOMPI')
#else

    INTEGER :: p_error, i

    ALLOCATE (p_recvcounts (groupsize))
    ALLOCATE (p_rdispls    (groupsize))

    p_sendcount = bufsize * sbuflen
    IF (PRESENT (recvcounts)) THEN
      p_recvcounts = recvcounts
    ELSE
      CALL MPI_GATHER (sbuflen,      1, p_int, &
                       p_recvcounts, 1, p_int, &
                       root, comm, p_error)

!MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE,
!           ROOT, COMM, IERROR)

    ENDIF
    if (p_pe == root) then
       IF (SUM (p_recvcounts) /= rbuflen) THEN
!          WRITE (nerr,'(a)') ' gatherv_args: sum (recvcounts) /= rbuflen'
          WRITE (nerr, *) ' gatherv_args: sum (recvcounts) /= rbuflen', rbuflen,sum (p_recvcounts), p_recvcounts
          CALL p_abort
       ENDIF
       p_recvcounts = bufsize * p_recvcounts

       p_rdispls (1) = 0
       DO i = 2, groupsize
          p_rdispls (i) = p_rdispls (i-1) + p_recvcounts (i-1)
       END DO
    end if

#endif

  END SUBROUTINE gatherv_args
!------------------------------------------------------------------------------
!==============================================================================
  pure function crc32_int (v) result (crc)
    integer, intent(in) :: v(:)
    integer             :: crc
    !-----------------------------------------------------------
    ! 32-bit checksum (CRC32) of a bitstream (array of integers)
    !-----------------------------------------------------------
!#if defined (__SX__) || defined (NAG)
!   integer, parameter :: CRC32_POLY =      Z'edb88320'    ! Non-std.
!#else
    integer, parameter :: CRC32_POLY = int (Z'edb88320')   ! F2003
!#endif

    integer :: i, k, poly

    crc = 0
    do i = 1, size (v)
       crc  = ieor (v(i), crc)
       do k = 1, bit_size (v)
          if (btest (crc, 0)) then
             poly = CRC32_POLY
          else
             poly = 0
          end if
          crc = ieor (ishft (crc, -1), poly)
       end do
    end do
  end function crc32_int
!------------------------------------------------------------------------------
  pure function crc32_int_i2 (v, len) result (crc)
    integer(i2), intent(in)           :: v(:)
    integer,     intent(in), optional :: len  ! Max. # of integers to process
    integer                           :: crc
    !-----------------------------------------------------------
    ! 32-bit checksum (CRC32) of a bitstream (array of integers)
    !-----------------------------------------------------------
!#if defined (__SX__) || defined (NAG)
!   integer, parameter :: CRC32_POLY =      Z'edb88320'    ! Non-std.
!#else
    integer, parameter :: CRC32_POLY = int (Z'edb88320')   ! F2003
!#endif

    integer :: i, k, poly, n

    n = size (v)
    if (present (len)) n = min (n, len)
    crc = 0
    do i = 1, n
       crc  = ieor (int(v(i)), crc)
       do k = 1, bit_size (v)
          if (btest (crc, 0)) then
             poly = CRC32_POLY
          else
             poly = 0
          end if
          crc = ieor (ishft (crc, -1), poly)
       end do
    end do
  end function crc32_int_i2
!------------------------------------------------------------------------------
  pure function crc32_char (v, len) result (crc)
    character, intent(in)           :: v(:)
    integer,   intent(in), optional :: len  ! Max. # of characters to process
    integer                         :: crc
    !---------------------------------------------------------
    ! 32-bit checksum (CRC32) of a bitstream (character array)
    !---------------------------------------------------------
!#if defined (__SX__) || defined (NAG)
!   integer, parameter :: CRC32_POLY =      Z'edb88320'    ! Non-std.
!#else
    integer, parameter :: CRC32_POLY = int (Z'edb88320')   ! F2003
!#endif

    integer :: i, k, poly, c, n

    n = size (v)
    if (present (len)) n = min (n, len)
    crc = 0
    do i = 1, n
       c = ichar (v(i))
#if defined(__ibm__)
       !-------------------------------------------------------
       ! On IBM the ichar intrinsic may return negative numbers
       !-------------------------------------------------------
       if (c < 0) c = c + 256
#endif
       crc  = ieor (crc, c)
       do k = 1, 8
          if (btest (crc, 0)) then
             poly = CRC32_POLY
          else
             poly = 0
          end if
          crc = ieor (ishft (crc, -1), poly)
       end do
    end do
  end function crc32_char
!------------------------------------------------------------------------------
  pure function crc32_real (y)
    real(wp), intent(in) :: y(:)
    integer              :: crc32_real
    !---------------------------------------------
    ! 32-bit checksum (CRC32) of an array of reals
    !---------------------------------------------
    integer, parameter   :: b32(1) = 0

    crc32_real = crc (transfer (y, b32))
  end function crc32_real
!------------------------------------------------------------------------------
  pure function crc32_real_sp (y)
    real(sp), intent(in) :: y(:)
    integer              :: crc32_real_sp
    !---------------------------------------------
    ! 32-bit checksum (CRC32) of an array of reals
    !---------------------------------------------
    integer, parameter   :: b32(1) = 0

    crc32_real_sp = crc (transfer (y, b32))
  end function crc32_real_sp
!------------------------------------------------------------------------------
  pure function crc32_logical (y)
    logical, intent(in) :: y(:)
    integer             :: crc32_logical
    !---------------------------------------------
    ! 32-bit checksum (CRC32) of an array of reals
    !---------------------------------------------
    integer, parameter   :: b32(1) = 0

    crc32_logical = crc (transfer (y, b32))
  end function crc32_logical
!==============================================================================
END MODULE mo_mpi_dace
!==============================================================================
!
! The following routines can be called as external routines to perform
! send/receive/broadcast operations on derived data types. Suitable
! interfaces must be defined in the calling routines.
!
!------------------------------------------------------------------------------
!
! MPI send routine, recommended interface:
!
! INTERFACE p_send
!   SUBROUTINE p_send_derivedtype (buffer, count, dest, tag, comm)
!   TYPE (MYTYPE) ,INTENT(in)           :: buffer    ! variable to send
!   INTEGER       ,INTENT(in)           :: count     ! len(byte) of variable
!   INTEGER       ,INTENT(in)           :: dest      ! destination processor
!   INTEGER       ,INTENT(in)           :: tag       ! tag
!   INTEGER       ,INTENT(in)           :: comm      ! communicator
!   END SUBROUTINE p_send_derivedtype
! END interface p_send
!
  SUBROUTINE p_send_derivedtype (buffer, count, dest, tag, comm)

  USE mo_mpi_dace,      ONLY: MPI_SUCCESS, MPI_BYTE, dace, &
                              p_abort
#ifdef MPI_CHECKSUM
  USE mo_mpi_dace,      ONLY: crc, p_send, p_recv, dump_buffer_int
#endif
  USE mo_fortran_units, ONLY: nerr
  IMPLICIT NONE

  INTEGER, INTENT(in)           :: buffer(*) ! variable to send
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: dest      ! destination processor index
  INTEGER, INTENT(in)           :: tag       ! tag
  INTEGER, INTENT(in)           :: comm      ! communicator

#ifndef NOMPI
    INTEGER :: ierror
#ifdef MPI_CHECKSUM
    integer :: checksum, chk
#endif
    CALL MPI_SEND (buffer, count, MPI_BYTE ,dest ,tag ,comm ,ierror)
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', dace% pe, &
            ' to ', dest, ' for tag ', tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierror
       CALL p_abort
    END IF
#ifdef MPI_CHECKSUM
    if (mod (count, 4) /= 0) then
       WRITE (nerr,'(a,i8)') &
            "p_send_derivedtype:WARNING: count not multiple of 4:", count
    end if
    checksum = crc (buffer(1:count/4))
    call p_send (checksum, dest, p_tag=32)
    call p_recv (chk     , dest, p_tag=33)
    if (checksum /= chk) then
       write (nerr,'()')
       write (nerr,'(a,i12,i3,a,i12,i3,1x,i0)') &
            "p_send_derivedtype:BAD CHECKSUM:", chk, dace% pe, &
            "  dest:", checksum, dest, count
       call dump_buffer_int (buffer(1:count/4), unit=2000+dace% pe, &
                             title="p_send_derivedtype: sent data")
       call p_send (chk,   dest, p_tag=34)  ! Synchronize before crash
       CALL p_abort
    end if
#endif /* MPI_CHECKSUM */
#endif
  END SUBROUTINE p_send_derivedtype
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 2 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_send_derivedtype2 (buffer, count, dest, tag, comm)
  IMPLICIT NONE
  INTEGER, INTENT(in)           :: buffer(*) ! variable to send
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: dest      ! destination processor index
  INTEGER, INTENT(in)           :: tag       ! tag
  INTEGER, INTENT(in)           :: comm      ! communicator
  INTERFACE
    SUBROUTINE p_send_derivedtype (buffer, count, dest, tag, comm)
    INTEGER, INTENT(in)           :: buffer(*) ! variable to send
    INTEGER, INTENT(in)           :: count     ! len(byte) of variable
    INTEGER, INTENT(in)           :: dest      ! destination processor index
    INTEGER, INTENT(in)           :: tag       ! tag
    INTEGER, INTENT(in)           :: comm      ! communicator
    END SUBROUTINE p_send_derivedtype
  END INTERFACE
  CALL p_send_derivedtype (buffer, count, dest, tag, comm)
  END SUBROUTINE p_send_derivedtype2
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 3 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_send_derivedtype3 (buffer, count, dest, tag, comm)
  IMPLICIT NONE
  INTEGER, INTENT(in)           :: buffer(*) ! variable to send
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: dest      ! destination processor index
  INTEGER, INTENT(in)           :: tag       ! tag
  INTEGER, INTENT(in)           :: comm      ! communicator
  INTERFACE
    SUBROUTINE p_send_derivedtype (buffer, count, dest, tag, comm)
    INTEGER, INTENT(in)           :: buffer(*) ! variable to send
    INTEGER, INTENT(in)           :: count     ! len(byte) of variable
    INTEGER, INTENT(in)           :: dest      ! destination processor index
    INTEGER, INTENT(in)           :: tag       ! tag
    INTEGER, INTENT(in)           :: comm      ! communicator
    END SUBROUTINE p_send_derivedtype
  END INTERFACE
  CALL p_send_derivedtype (buffer, count, dest, tag, comm)
  END SUBROUTINE p_send_derivedtype3
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 4 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_send_derivedtype4 (buffer, count, dest, tag, comm)
  IMPLICIT NONE
  INTEGER, INTENT(in)           :: buffer(*) ! variable to send
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: dest      ! destination processor index
  INTEGER, INTENT(in)           :: tag       ! tag
  INTEGER, INTENT(in)           :: comm      ! communicator
  INTERFACE
    SUBROUTINE p_send_derivedtype (buffer, count, dest, tag, comm)
    INTEGER, INTENT(in)           :: buffer(*) ! variable to send
    INTEGER, INTENT(in)           :: count     ! len(byte) of variable
    INTEGER, INTENT(in)           :: dest      ! destination processor index
    INTEGER, INTENT(in)           :: tag       ! tag
    INTEGER, INTENT(in)           :: comm      ! communicator
    END SUBROUTINE p_send_derivedtype
  END INTERFACE
  CALL p_send_derivedtype (buffer, count, dest, tag, comm)
  END SUBROUTINE p_send_derivedtype4
!------------------------------------------------------------------------------
  !---------------------------------------------
  ! Non-blocking versions of p_send_derivedtype*
  !---------------------------------------------
  SUBROUTINE p_isend_derivedtype (buffer, count, dest, tag, comm)

  USE mo_mpi_dace,      ONLY: MPI_SUCCESS, MPI_BYTE, dace, &
                              p_abort, p_request, p_irequest
#ifdef MPI_CHECKSUM
  USE mo_mpi_dace,      ONLY: crc, p_send, p_recv, dump_buffer_int
#endif
  USE mo_fortran_units, ONLY: nerr
#ifndef NOMPI
! USE mpi
#endif
  IMPLICIT NONE

  INTEGER, INTENT(in)           :: buffer(*) ! variable to send
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: dest      ! destination processor index
  INTEGER, INTENT(in)           :: tag       ! tag
  INTEGER, INTENT(in)           :: comm      ! communicator

#ifndef NOMPI
    INTEGER :: ierror
#ifdef MPI_CHECKSUM
    integer :: checksum, chk
#endif
    CALL MPI_ISEND (buffer, count, MPI_BYTE, dest, tag, &
                    comm, p_request(p_irequest), ierror)
    p_irequest = p_irequest + 1
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', dace% pe, &
            ' to ', dest, ' for tag ', tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierror
       CALL p_abort
    END IF
#if 0 /* #ifdef MPI_CHECKSUM */
    ! Careful: code below has not been extended to non-blocking communication
    if (mod (count, 4) /= 0) then
       WRITE (nerr,'(a,i8)') &
            "p_isend_derivedtype:WARNING: count not multiple of 4:", count
    end if
    checksum = crc (buffer(1:count/4))
    call p_send (checksum, dest, p_tag=32)
    call p_recv (chk     , dest, p_tag=33)
    if (checksum /= chk) then
       write (nerr,'()')
       write (nerr,'(a,i12,i3,a,i12,i3,1x,i0)') &
            "p_isend_derivedtype:BAD CHECKSUM:", chk, dace% pe, &
            "  dest:", checksum, dest, count
       call dump_buffer_int (buffer(1:count/4), unit=2000+dace% pe, &
                             title="p_isend_derivedtype: sent data")
       call p_send (chk,   dest, p_tag=34)  ! Synchronize before crash
       CALL p_abort
    end if
#endif /* MPI_CHECKSUM */
#endif
  END SUBROUTINE p_isend_derivedtype
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 2 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_isend_derivedtype2 (buffer, count, dest, tag, comm)
    IMPLICIT NONE
    INTEGER, INTENT(in)           :: buffer(*) ! variable to send
    INTEGER, INTENT(in)           :: count     ! len(byte) of variable
    INTEGER, INTENT(in)           :: dest      ! destination processor index
    INTEGER, INTENT(in)           :: tag       ! tag
    INTEGER, INTENT(in)           :: comm      ! communicator
    INTERFACE
       SUBROUTINE p_isend_derivedtype (buffer, count, dest, tag, comm)
         INTEGER, INTENT(in)      :: buffer(*) ! variable to send
         INTEGER, INTENT(in)      :: count     ! len(byte) of variable
         INTEGER, INTENT(in)      :: dest      ! destination processor index
         INTEGER, INTENT(in)      :: tag       ! tag
         INTEGER, INTENT(in)      :: comm      ! communicator
       END SUBROUTINE p_isend_derivedtype
    END INTERFACE
    CALL p_isend_derivedtype (buffer, count, dest, tag, comm)
  END SUBROUTINE p_isend_derivedtype2
!------------------------------------------------------------------------------
!
! MPI recv routine, recommended interface:
!
! INTERFACE p_recv
!   SUBROUTINE p_recv_derivedtype (buffer, count, dest, tag, comm)
!   TYPE(MYTYPE) ,INTENT(out)          :: buffer    ! variable to receive
!   INTEGER      ,INTENT(in)           :: count     ! len(byte) of variable
!   INTEGER      ,INTENT(in)           :: dest      ! destination processor
!   INTEGER      ,INTENT(in)           :: tag       ! tag
!   INTEGER      ,INTENT(in)           :: comm      ! communicator
!   END SUBROUTINE p_recv_derivedtype
! END interface p_recv
!
  SUBROUTINE p_recv_derivedtype (buffer, count, source, tag, comm)

  USE mo_mpi_dace,      ONLY: MPI_SUCCESS, MPI_BYTE, dace, &
                              p_abort, MPI_STATUS_SIZE
#ifdef MPI_CHECKSUM
  USE mo_mpi_dace,      ONLY: crc, p_recv, p_send, dump_buffer_int
#endif
  USE mo_fortran_units, ONLY: nerr
  IMPLICIT NONE

    INTEGER, INTENT(out) :: buffer(*)
    INTEGER, INTENT(in)           :: count     ! len(byte) of variable
    INTEGER, INTENT(in)           :: source    ! destination processor index
    INTEGER, INTENT(in)           :: tag       ! tag
    INTEGER, INTENT(in)           :: comm      ! communicator

#ifndef NOMPI
    INTEGER :: ierror
    INTEGER :: status(MPI_STATUS_SIZE)
#ifdef MPI_CHECKSUM
    integer :: chk, checksum
#endif
    !--------
    ! receive
    !--------
    CALL MPI_RECV (buffer,count,MPI_BYTE,source,tag,comm,status,ierror)
    !---------------
    ! error handling
    !---------------
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', dace% pe, &
            ' from ', source, ' for tag ', tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierror
       CALL p_abort
    END IF
#ifdef MPI_CHECKSUM
    if (mod (count, 4) /= 0) then
       WRITE (nerr,'(a,i8)') &
            "p_recv_derivedtype:WARNING: count not multiple of 4:", count
    end if
    chk = crc (buffer(1:count/4))
    call p_recv (checksum, source, p_tag=32)
    call p_send (chk     , source, p_tag=33)
    if (checksum /= chk) then
       write (nerr,'()')
       write (nerr,'(a,i12,i3,a,i12,i3,1x,i0)') &
            "p_recv_derivedtype:BAD CHECKSUM:", chk, dace% pe, &
            "  source:", checksum, source, count
       call dump_buffer_int (buffer(1:count/4), unit=3000+dace% pe, &
                             title="p_recv_derivedtype: received data")
       call p_recv (chk,   source, p_tag=34)  ! Synchronize before crash
       CALL p_abort
    end if
#endif
#endif
  END SUBROUTINE p_recv_derivedtype
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 2 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_recv_derivedtype2 (buffer, count, source, tag, comm)
  INTEGER, INTENT(out)          :: buffer(*)
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: source    ! destination processor index
  INTEGER, INTENT(in)           :: tag       ! tag
  INTEGER, INTENT(in)            :: comm      ! communicator
    INTERFACE
      SUBROUTINE p_recv_derivedtype (buffer, count, source, tag, comm)
      INTEGER, INTENT(out)          :: buffer(*)
      INTEGER, INTENT(in)           :: count     ! len(byte) of variable
      INTEGER, INTENT(in)           :: source    ! destination processor index
      INTEGER, INTENT(in)           :: tag       ! tag
      INTEGER, INTENT(in)           :: comm      ! communicator
      END SUBROUTINE p_recv_derivedtype
    END INTERFACE
    CALL p_recv_derivedtype (buffer, count, source, tag, comm)
  END SUBROUTINE p_recv_derivedtype2
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 3 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_recv_derivedtype3 (buffer, count, source, tag, comm)
  INTEGER, INTENT(out)          :: buffer(*)
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: source    ! destination processor index
  INTEGER, INTENT(in)           :: tag       ! tag
  INTEGER, INTENT(in)           :: comm      ! communicator
    INTERFACE
      SUBROUTINE p_recv_derivedtype (buffer, count, source, tag, comm)
      INTEGER, INTENT(out)          :: buffer(*)
      INTEGER, INTENT(in)           :: count     ! len(byte) of variable
      INTEGER, INTENT(in)           :: source    ! destination processor index
      INTEGER, INTENT(in)           :: tag       ! tag
      INTEGER, INTENT(in)           :: comm      ! communicator
      END SUBROUTINE p_recv_derivedtype
    END INTERFACE
    CALL p_recv_derivedtype (buffer, count, source, tag, comm)
  END SUBROUTINE p_recv_derivedtype3
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 4 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_recv_derivedtype4 (buffer, count, source, tag, comm)
  INTEGER, INTENT(out)          :: buffer(*)
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: source    ! destination processor index
  INTEGER, INTENT(in)           :: tag       ! tag
  INTEGER, INTENT(in)           :: comm      ! communicator
    INTERFACE
      SUBROUTINE p_recv_derivedtype (buffer, count, source, tag, comm)
      INTEGER, INTENT(out)          :: buffer(*)
      INTEGER, INTENT(in)           :: count     ! len(byte) of variable
      INTEGER, INTENT(in)           :: source    ! destination processor index
      INTEGER, INTENT(in)           :: tag       ! tag
      INTEGER, INTENT(in)           :: comm      ! communicator
      END SUBROUTINE p_recv_derivedtype
    END INTERFACE
    CALL p_recv_derivedtype (buffer, count, source, tag, comm)
  END SUBROUTINE p_recv_derivedtype4
!------------------------------------------------------------------------------
  !---------------------------------------------
  ! Non-blocking versions of p_recv_derivedtype*
  !---------------------------------------------
  SUBROUTINE p_irecv_derivedtype (buffer, count, source, tag, comm)

  USE mo_mpi_dace,      ONLY: MPI_SUCCESS, MPI_BYTE, dace, &
                              p_abort, p_request, p_irequest
#ifdef MPI_CHECKSUM
  USE mo_mpi_dace,      ONLY: crc, p_recv, p_send, dump_buffer_int
#endif
  USE mo_fortran_units, ONLY: nerr
#ifndef NOMPI
! USE mpi
#endif
  IMPLICIT NONE

    INTEGER, INTENT(inout) :: buffer(*)
    INTEGER, INTENT(in)    :: count     ! len(byte) of variable
    INTEGER, INTENT(in)    :: source    ! destination processor index
    INTEGER, INTENT(in)    :: tag       ! tag
    INTEGER, INTENT(in)    :: comm      ! communicator

#ifndef NOMPI
    INTEGER :: ierror
#ifdef MPI_CHECKSUM
    integer :: chk, checksum
#endif
    !--------
    ! receive
    !--------
    CALL MPI_IRECV (buffer, count, MPI_BYTE, source, tag, &
                    comm, p_request(p_irequest), ierror)
    p_irequest = p_irequest + 1
    !---------------
    ! error handling
    !---------------
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', dace% pe, &
            ' from ', source, ' for tag ', tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierror
       CALL p_abort
    END IF
#if 0 /* #ifdef MPI_CHECKSUM */
    ! Careful: code below has not been extended to non-blocking communication
    if (mod (count, 4) /= 0) then
       WRITE (nerr,'(a,i8)') &
            "p_irecv_derivedtype:WARNING: count not multiple of 4:", count
    end if
    chk = crc (buffer(1:count/4))
    call p_recv (checksum, source, p_tag=32)
    call p_send (chk     , source, p_tag=33)
    if (checksum /= chk) then
       write (nerr,'()')
       write (nerr,'(a,i12,i3,a,i12,i3,1x,i0)') &
            "p_irecv_derivedtype:BAD CHECKSUM:", chk, dace% pe, &
            "  source:", checksum, source, count
       call dump_buffer_int (buffer(1:count/4), unit=3000+dace% pe, &
                             title="p_irecv_derivedtype: received data")
       call p_recv (chk,   source, p_tag=34)  ! Synchronize before crash
       CALL p_abort
    end if
#endif
#endif
  END SUBROUTINE p_irecv_derivedtype
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 2 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_irecv_derivedtype2 (buffer, count, source, tag, comm)
    INTEGER, INTENT(inout)        :: buffer(*)
    INTEGER, INTENT(in)           :: count     ! len(byte) of variable
    INTEGER, INTENT(in)           :: source    ! destination processor index
    INTEGER, INTENT(in)           :: tag       ! tag
    INTEGER, INTENT(in)           :: comm      ! communicator
    INTERFACE
       SUBROUTINE p_irecv_derivedtype (buffer, count, source, tag, comm)
         INTEGER, INTENT(inout)   :: buffer(*)
         INTEGER, INTENT(in)      :: count     ! len(byte) of variable
         INTEGER, INTENT(in)      :: source    ! destination processor index
         INTEGER, INTENT(in)      :: tag       ! tag
         INTEGER, INTENT(in)      :: comm      ! communicator
       END SUBROUTINE p_irecv_derivedtype
    END INTERFACE
    CALL p_irecv_derivedtype (buffer, count, source, tag, comm)
  END SUBROUTINE p_irecv_derivedtype2
!------------------------------------------------------------------------------
!
! MPI bcast routine, recommended interface:
!
! INTERFACE p_bcast
!   SUBROUTINE p_bcast_derivedtype (buffer, count, source, comm)
!   TYPE(MYTYPE) ,INTENT(inout)        :: buffer    ! variable to bcast
!   INTEGER      ,INTENT(in)           :: count     ! len(byte) of variable
!   INTEGER      ,INTENT(in)           :: source    ! source processor index
!   INTEGER      ,INTENT(in)           :: comm      ! communicator
!   END SUBROUTINE p_bcast_derivedtype
! END interface p_bcast
!
  SUBROUTINE p_bcast_derivedtype (buffer, count, source, comm)

#ifndef NOMPI
    USE mo_mpi_dace, ONLY: MPI_BYTE, MPI_SUCCESS, MPI_STATUS_SIZE, &
                           dace, p_abort
#ifdef DEBUG
    USE mo_mpi_dace, ONLY: nbcast
#endif
#ifdef MPI_CHECKSUM
    USE mo_mpi_dace, ONLY: crc, p_bcast, p_gather, p_barrier
    USE mo_kind,     ONLY: i2
#endif
    USE mo_fortran_units, ONLY: nerr
#endif
    IMPLICIT NONE

    INTEGER, INTENT(inout)        :: buffer(*)
    INTEGER, INTENT(in)           :: count     ! len(byte) of variable
    INTEGER, INTENT(in)           :: source    ! source processor index
    INTEGER, INTENT(in)           :: comm      ! communicator

#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: ierror
!   INTEGER :: status(MPI_STATUS_SIZE)
#ifdef MPI_CHECKSUM
    integer :: chk, checksum
    logical :: ok, bad(0:dace% npe-1)
    character, parameter :: c1(1) = "*"
    character :: charbuf(count)                ! Work around sxf90 bug
#endif

!#define MPISX_BCAST_WORKAROUND
#ifdef  MPISX_BCAST_WORKAROUND
    !-------------------------------------------
    ! Tentative workaround for MPI/SX bcast bug:
    ! Maximum message size for bcasts on MPI/SX:
    ! (must be multiple of 4)
    !-------------------------------------------
    integer, parameter :: maxmsg = 65536    ! appears to work
!   integer, parameter :: maxmsg = 1000000  ! fails
    integer            :: i, off, msgsize
#endif
    !----------------------------
    ! process optional parameters
    !----------------------------
    p_comm = comm
    !----------
    ! broadcast
    !----------
#ifdef DEBUG
    nbcast = nbcast+1
#endif
    IF (dace% npe == 1) RETURN

#ifdef MPI_CHECKSUM
    ! Write diagnostics for "large" bcasts
    if (dace% pe == source .and. count > 8000) then
       write (0,*) "p_bcast_derivedtype:INFO:pe, message size =", dace% pe, count
    end if
#endif
#ifdef MPISX_BCAST_WORKAROUND
    !----------------------------------------------
    ! Tentative workaround for MPI/SX bcast bug:
    ! Split large bcasts into multiple smaller ones
    !----------------------------------------------
    do i = 1, count, maxmsg
       off     = (i-1)/4 + 1
       msgsize = min (count-i+1, maxmsg)
       CALL MPI_BCAST (buffer(off), msgsize, MPI_BYTE, source, p_comm, ierror)
       IF (ierror /= MPI_SUCCESS) EXIT
    end do
#else
    !-------------------------------------------
    ! Original code, sending full buffer at once
    !-------------------------------------------
    CALL MPI_BCAST (buffer, count, MPI_BYTE, source, p_comm, ierror)
#endif
    !---------------
    ! error handling
    !---------------
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST on ', dace% pe, &
            ' from ', source, ' failed in p_bcast_derivedtype.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierror
       CALL p_abort
    END IF
#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', source, &
            ' with broadcast number ', nbcast, ' successful.'
#endif

#ifdef MPI_CHECKSUM
    if (count == 0) return
    if (dace% pe == source) then
       if (mod (count, 4) /= 0) then
          WRITE (nerr,'(a,i8)') &
               "p_bcast_derivedtype:WARNING: count not multiple of 4:", count
       end if
    end if
    select case (mod (count,4))
    case (0)
       checksum = crc (buffer(1:count/4))
    case (2)
       checksum = crc (transfer (buffer(1:count/4+1), (/0_i2/), size=count/2))
    case default
       ! The following line produces a segfault with sxf90 rev.381:
!!!    checksum = crc (transfer (buffer(1:count/4+1), c1, size=count))
       ! Workaround for sxf90 rev.381 bug:
       charbuf  = transfer (buffer(1:count/4+1), c1, size=count)
       checksum = crc (charbuf)
    end select
    chk = checksum
    call p_bcast (checksum, source)
    ok = (checksum == chk)      ! Compare checksums of sender and receiver
    bad = .false.
    call p_gather (.not. ok, bad, root=source)
    if (.not. ok) then
       WRITE (nerr,'(a,i12,i3,a,i12,i3,1x,i0)') &
            "p_bcast_derivedtype:BAD CHECKSUM:", chk, dace% pe, &
            "  source:", checksum, source, count
!      call print_head_tail (buffer(1:count/4))
       WRITE (nerr,'(a,i5)') "Dumping receive buffer to unit", 1000+dace% pe
       call print_head_tail (buffer(1:count/4), unit=1000+dace% pe)
       WRITE (nerr,'()')
    end if
    if (dace% pe == source) then
       if (any (bad)) then
          ok = .false.
          WRITE (nerr,'(a,i5)') "Dumping send buffer to unit", 1000+dace% pe
          call print_head_tail (buffer(1:count/4), unit=1000+dace% pe)
          WRITE (nerr,'(a)') "One or more bad checksums, aborting..."
          WRITE (nerr,'()')
       end if
    end if
    call p_barrier ()           ! Synchronize before crashing
    call p_bcast (ok, source)
    if (.not. ok) CALL p_abort ()
  contains
    subroutine print_head_tail (v, unit)
      integer, intent(in)           :: v(:)
      integer, intent(in), optional :: unit
      integer :: n, i, j, lu
      n = size (v)
      lu = 0; if (present (unit)) lu = unit
      if (lu>0) then
         if (dace% pe == source) then
            write (lu,*) "p_bcast_derivedtype: source data"
         else
            write (lu,*) "p_bcast_derivedtype: receiver data"
         end if
      end if
      WRITE (lu,*) dace% pe, ":", v(:min (n,4)), "...", v(max (5,n-3):)
      WRITE (lu,'()')
      WRITE (lu,'(A)') "Nontrivial contents:"
      do i = 1, n, 8
         j = min (i+7,n)
         if (any (v(i:j) /= 0)) write (lu,*) i,"..",j,":", v(i:j)
      end do
    end subroutine print_head_tail
#endif /* MPI_CHECKSUM */
#endif /* !NOMPI */
  END SUBROUTINE p_bcast_derivedtype
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 2 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_bcast_derivedtype2 (buffer, count, source, comm)
  INTEGER, INTENT(inout)        :: buffer(*)
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: source    ! source processor index
  INTEGER, INTENT(in)           :: comm      ! communicator
    INTERFACE
      SUBROUTINE p_bcast_derivedtype (buffer, count, source, comm)
      INTEGER, INTENT(inout)        :: buffer(*)
      INTEGER, INTENT(in)           :: count     ! len(byte) of variable
      INTEGER, INTENT(in)           :: source    ! source processor index
      INTEGER, INTENT(in)           :: comm      ! communicator
      END SUBROUTINE p_bcast_derivedtype
    END INTERFACE
    CALL p_bcast_derivedtype (buffer, count, source, comm)
  END SUBROUTINE p_bcast_derivedtype2
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 3 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_bcast_derivedtype3 (buffer, count, source, comm)
  INTEGER, INTENT(inout)        :: buffer(*)
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: source    ! source processor index
  INTEGER, INTENT(in)           :: comm      ! communicator
    INTERFACE
      SUBROUTINE p_bcast_derivedtype (buffer, count, source, comm)
      INTEGER, INTENT(inout)        :: buffer(*)
      INTEGER, INTENT(in)           :: count     ! len(byte) of variable
      INTEGER, INTENT(in)           :: source    ! source processor index
      INTEGER, INTENT(in)           :: comm      ! communicator
      END SUBROUTINE p_bcast_derivedtype
    END INTERFACE
    CALL p_bcast_derivedtype (buffer, count, source, comm)
  END SUBROUTINE p_bcast_derivedtype3
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 4 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_bcast_derivedtype4 (buffer, count, source, comm)
  INTEGER, INTENT(inout)        :: buffer(*)
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: source    ! source processor index
  INTEGER, INTENT(in)           :: comm      ! communicator
    INTERFACE
      SUBROUTINE p_bcast_derivedtype (buffer, count, source, comm)
      INTEGER, INTENT(inout)        :: buffer(*)
      INTEGER, INTENT(in)           :: count     ! len(byte) of variable
      INTEGER, INTENT(in)           :: source    ! source processor index
      INTEGER, INTENT(in)           :: comm      ! communicator
      END SUBROUTINE p_bcast_derivedtype
    END INTERFACE
    CALL p_bcast_derivedtype (buffer, count, source, comm)
  END SUBROUTINE p_bcast_derivedtype4
  !------------------------------------------------------------------------
  ! define the same routine with another name for the case that 5 different
  ! interfaces are required in a programming unit
  !------------------------------------------------------------------------
  SUBROUTINE p_bcast_derivedtype5 (buffer, count, source, comm)
  INTEGER, INTENT(inout)        :: buffer(*)
  INTEGER, INTENT(in)           :: count     ! len(byte) of variable
  INTEGER, INTENT(in)           :: source    ! source processor index
  INTEGER, INTENT(in)           :: comm      ! communicator
    INTERFACE
      SUBROUTINE p_bcast_derivedtype (buffer, count, source, comm)
      INTEGER, INTENT(inout)        :: buffer(*)
      INTEGER, INTENT(in)           :: count     ! len(byte) of variable
      INTEGER, INTENT(in)           :: source    ! source processor index
      INTEGER, INTENT(in)           :: comm      ! communicator
      END SUBROUTINE p_bcast_derivedtype
    END INTERFACE
    CALL p_bcast_derivedtype (buffer, count, source, comm)
  END SUBROUTINE p_bcast_derivedtype5
!==============================================================================
