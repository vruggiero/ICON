!
!+ Replacement for COSMO module 'environment'
!
MODULE environment
!
! Description:
!   Replacement for COSMO module environment
!   Defines:
!    - subroutine model_abort
!    - subroutine get_free_unit
!    - subroutine comm_barrier
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_9         2010/04/20 Andreas Rhodin
!  Replacement for COSMO module 'environment'
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_20        2012-06-18 Andreas Rhodin
!  some cleanup
! V1_26        2013/06/27 Andreas Rhodin
!  some more public entities for COSMO complient modules
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!------------------------------------------------------------------------------
  !-----------------------------
  ! se from 3D-Var/LETKF modules
  !-----------------------------
  use mo_exception,     only: finish          ! abort in case of error
  use mo_fortran_units, only: get_unit_number ! get Fortran unit number
  use mo_mpi_dace,      only: p_barrier,   &  ! MPI barrier routine
!                    nproc => p_nprocs,    &  ! # of MPI processor elements
                              MPI_BYTE,    &  ! MPI byte id
                              MPI_INTEGER, &  ! MPI integer id
                              MPI_SUCCESS     ! MPI return code
!                             MPI_COMM_WORLD  ! MPI communicator
  implicit none

  !---------------------------------
  ! Public entities, COSMO complient
  !---------------------------------
  private
  public :: model_abort    ! program abort in case of error
  public :: get_free_unit  ! get a free fortran unit
  public :: comm_barrier   ! MPI barrier routine
! public :: nproc          ! number of MPI processor elements
  public :: MPI_BYTE       ! MPI byte id
  public :: MPI_INTEGER    ! MPI integer id
  public :: MPI_SUCCESS    ! MPI return code
! public :: MPI_COMM_WORLD ! MPI communicator

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  subroutine model_abort (my_id, ierrorcode, yerrorstring, yroutine, implerrorcode)
  integer           ,intent(in) :: my_id         ! id of this processor
  integer           ,intent(in) :: ierrorcode    ! self-defined integer code of the error detected
  character(len=*)  ,intent(in) :: yerrorstring  ! self-defined error message
  character(len=*)  ,intent(in) :: yroutine      ! calling routine
  integer ,optional ,intent(in) :: implerrorcode ! error-code of the message passing library

    call finish (yroutine, yerrorstring)

  end subroutine model_abort

!------------------------------------------------------------------------------

  subroutine get_free_unit (iunit)
  integer ,intent(out) :: iunit ! is set to the next free unit number

    iunit = get_unit_number()

  end subroutine get_free_unit

!------------------------------------------------------------------------------

  subroutine comm_barrier (icomm, ierror, yerrmsg)
  integer          ,intent (in)  :: icomm    ! communicator to be used
  integer          ,intent (out) :: ierror   ! error-status variable
  character(len=*) ,intent (out) :: yerrmsg  ! for mpi error message
    call p_barrier (icomm)
    yerrmsg = ''            ! set intent out variables
    ierror  = 0
  end subroutine comm_barrier

!------------------------------------------------------------------------------

end module environment
