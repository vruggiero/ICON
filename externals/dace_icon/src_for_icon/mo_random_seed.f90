!
!+ Generate deterministic seed for the random number generator
!
MODULE mo_random_seed
!
! Description:
!   Generate deterministic seed for the random number generator
!   derived from:
!                                                       default no. bits used
!   1) step (analysis, forecast)                        2^ 1 =    2
!   2) explicit seed (different experiments)            2^ 2 =    4
!   3) analysis date and time (hour)                    2^ 8 =  256 (10.6 days)
!   4) ensemble member number                           2^ 7 =  128
!   5) unique stream number within the parallel program 2^12 = 4096
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
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_31        2014-08-21 Andreas Rhodin
!  construct_stream: new optional parameters 'uniq', 'hour'
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2007  original source
!------------------------------------------------------------------------------

!-------------
! Modules used
!-------------
use mo_run_params, only: ana_time,       &! analysis time
                         expseed          ! non-default random number seed
use mo_time,       only: hours            ! convert date+time to hours
use mo_random,     only: random_state_t, &! random generator state derived type
                         construct        ! New random generator from seed
use mo_exception,  only: finish           ! abort in case of error
use mo_mpi_dace,   only: p_min, p_max     ! generic MPI min, max routines
use mo_mpi_dace,   only: dace
implicit none

!----------------
! Public entities
!----------------
private
public :: construct_stream   ! New random generator, derive seed
public :: init_streams       ! Initialize new sequence of streams

!---------------
! some constants
!---------------
integer ,parameter :: BITS_STEP   =  1 ! bits for analysis/forecast step
integer ,parameter :: BITS_EXP    =  2 ! bits for different experiments
integer ,parameter :: BITS_DATE   =  8 ! bits for date and time
integer ,parameter :: BITS_ENSMEM =  7 ! bits for ensemble member
integer ,parameter :: BITS_UNIQUE = 12 ! bits for unique stream

integer ,parameter :: SHFT_STEP   =  0
integer ,parameter :: SHFT_EXP    =  SHFT_STEP   + BITS_STEP
integer ,parameter :: SHFT_DATE   =  SHFT_EXP    + BITS_EXP
integer ,parameter :: SHFT_ENSMEM =  SHFT_DATE   + BITS_DATE
integer ,parameter :: SHFT_UNIQUE =  SHFT_ENSMEM + BITS_ENSMEM

!-----------------
! module variables
!-----------------
integer ,save :: unique = 0 ! unique stream number within the parallel program
integer ,save :: step   = 0 ! e.g. 0=analysis; 1=forecast

contains
!------------------------------------------------------------------------------
  subroutine construct_stream (state, ensmem, uniq, use, hour, comm)
  type(random_state_t),intent(out) ,optional :: state  ! random generator state
  integer             ,intent(in)  ,optional :: ensmem ! ensemble member number
  integer             ,intent(in)  ,optional :: uniq   ! unique value to use
  logical             ,intent(in)  ,optional :: use    ! use on this PE
  integer             ,intent(in)  ,optional :: hour   ! use instead of ana_time
  integer             ,intent(in)  ,optional :: comm   ! communicator to use
  !------------------------------------------------------
  ! Generate a seed for the random number generator
  !   derived from:
  !   1) step (analysis, forecast)
  !   2) explicit seed (different experiments)
  !   3) analysis date and time (hour)
  !   4) ensemble member number
  !   5) unique stream number within the parallel program
  !------------------------------------------------------

    !----------------
    ! local variables
    !----------------
    integer :: seed       ! unique derived seed
    logical :: luse       ! random stream used on this processor element
    integer :: seed_time  ! seed derived from analysis date & time
    integer :: seed_mem   ! seed derived from ensemble member number
    integer :: luniq      ! unique value to use here
    integer :: lcomm      ! local copy of communicator

    !-------------------
    ! process parameters
    !-------------------
    luniq = -1; if(present(uniq)) luniq = uniq
    select case (luniq)
    case (:-1)
      unique  = unique + 1
      luniq   = unique
    case (0)
      luniq   = unique
    end select
    if (present (hour)) then
      seed_time = mod (hour,                    2 ** BITS_DATE)
    else
      seed_time = mod (nint (hours (ana_time)), 2 ** BITS_DATE)
    endif
    luse      = present (state) ;if (present(use))    luse     = luse .and. use
    seed_mem  = 0               ;if (present(ensmem)) seed_mem = ensmem
    lcomm     = dace% comm      ;if (present(comm))   lcomm    = comm

    !---------------------------
    ! check seed for valid range
    !---------------------------
    if (seed_mem < 0 .or. seed_mem >= 2**BITS_ENSMEM)       &
      call finish('construct_stream','seed_mem out of range')
    if (luniq >= 2**BITS_UNIQUE)                              &
      call finish('construct_stream','unique seed out of range')
    if (expseed < 0 .or. expseed >= 2**BITS_EXP)           &
      call finish('construct_stream','expseed out of range')

    !-------------------
    ! derive unique seed
    !-------------------
    seed = 0
    seed = ieor (seed, step      * 2 ** SHFT_STEP)
    seed = ieor (seed, expseed   * 2 ** SHFT_EXP)
    seed = ieor (seed, seed_time * 2 ** SHFT_DATE)
    seed = ieor (seed, seed_mem  * 2 ** SHFT_ENSMEM)
    seed = ieor (seed, luniq     * 2 ** SHFT_UNIQUE)

    !----------------------------------------------------
    ! check for same seed on different processor elements
    !----------------------------------------------------
    if (p_min(seed,lcomm) /= p_max(seed,lcomm))             &
         call finish ('construct_stream','inconsistent seed')

    !------------------------------------
    ! finally construct the random stream
    !------------------------------------
    if (luse) then
      call construct (state, seed)
    endif

  end subroutine construct_stream
!------------------------------------------------------------------------------
  subroutine init_streams (stepnum, uniq)
    integer ,intent(in) ,optional :: stepnum  ! step or runtype
    integer ,intent(in) ,optional :: uniq     ! initial unique value
    !-----------------------------------
    ! Initialize new sequence of streams
    ! Increment step if stepnum not given.
    !-----------------------------------
    step   = step + 1
    unique = 0
    if (present (stepnum)) step   = stepnum
    if (present (uniq   )) unique = uniq
  end subroutine init_streams
!------------------------------------------------------------------------------
end module mo_random_seed
