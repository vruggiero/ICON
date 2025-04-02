!
!+ Communication routines for type t_obs
!
MODULE mo_obs_sndrcv
!
! Description:
! This module defines communication routines for type t_obs.
! To be merged with mo_t_obs if F2003 IMPORT statement becomes available.
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
!  remove obsolete component T_obs% obs
! V1_5         2009/05/25 Harald Anlauf
!  scatter_mcol: new, use MPI_scatterv instead of send/recv
! V1_7         2009/08/24 Harald Anlauf
!  scatter_mcol: allocate dummy array (required by INTEL compiler)
! V1_8         2009/12/09 Andreas Rhodin
!  replace sort routine 'index' by simpler algorithm
! V1_9         2010/04/20 Andreas Rhodin
!  dont send/receive/broadcast components of type t_spot
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  technical changes for revised minimisation routine
! V1_20        2012-06-18 Andreas Rhodin
!  cleanup: remove unused variables
! V1_27        2013-11-08 Andreas Rhodin
!  for VarBC: handle bias correction predictors
! V1_28        2014/02/26 Andreas Rhodin
!  handle sink variable component in t_obs send/recv/bcast/alltoall routines
!  changed interface of subroutine scatter_mcol (now independent from t_obs)
! V1_29        2014/04/02 Andreas Rhodin
!  precautions to transfer first guess from monitoring to analysis pass
! V1_42        2015-06-08 Andreas Rhodin
!  preparations for time interpolation for COSMO MEC
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2005       extracted from module mo_obs
! Andreas Rhodin  DWD  2006-2008  extensions
!==============================================================================
!
! Diagnostics for vectorization
!
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!==============================================================================
  !=============
  ! modules used
  !=============
  use mo_exception,  only: finish            ! abort in case of error
  !-----------------------------------------------
  ! observation data type and predefined constants
  !-----------------------------------------------
  use mo_t_obs,      only: t_obs,           &! observation data type
                           t_spot,          &!   component of t_obs
                           t_mcol,          &!   component of t_obs
                           t_mcols,         &!   component of t_obs
                           t_datum,         &!   component of t_obs
                           t_box,           &!   component of t_obs
                           !-------------------------
                           ! byte sizes of data types
                           !-------------------------
                           obs_bytes, spot_bytes, body_bytes,&
                           set_byte_size     ! routine to setup sizes
  use mo_sink,       only: t_sink            !   component of t_obs
  !-----------------------
  ! communication routines
  !-----------------------
  use mo_mpi_dace,   only: dace,            &! MPI group info
                           p_or,            &! parallel or
                           p_recv,          &! generic MPI receive   routine
                           p_send,          &! generic MPI send      routine
                           p_bcast,         &! generic MPI bcast     routine
                           p_alltoall        ! generic MPI alltoall  routine
  implicit none
!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: p_send         ! (mpi-) send    variable(s) of data type obs
  public :: p_recv         ! (mpi-) receive variable(s) of data type obs
  public :: p_bcast        ! (mpi-) receive variable(s) of data type obs
  public :: p_alltoall     ! (mpi-) alltoall
  public :: send_mcol
  public :: recv_mcol
  public :: scatter_mcol   ! Scatter array of type t_mcol.
  public :: alltoall_mcol  ! Scatter array of type t_mcol.
  public :: p_mcol         ! Array of pointers to model column descriptor
!------------------------------------------------------------------------------
  !======================
  ! data type definitions
  !======================
  type p_mcol
     type(t_mcol), pointer :: p(:)  ! pointer to model column descriptor
  end type p_mcol
!------------------------------------------------------------------------------
  !===========
  ! Interfaces
  !===========
!------------------------------------------------------------------------------
  interface p_send
    module procedure send_obs
    module procedure send_obss
    module procedure send_sink
  end interface p_send
!------------------------------------------------------------------------------
  interface p_recv
    module procedure recv_obs
    module procedure recv_obss
    module procedure recv_sink
  end interface p_recv
!------------------------------------------------------------------------------
  interface p_bcast
    module procedure bcast_obs
    module procedure bcast_obss
    module procedure bcast_sink
  end interface p_bcast
!------------------------------------------------------------------------------
  interface p_alltoall
    module procedure alltoall_datum
    module procedure alltoall_spot
    module procedure alltoall_obs
    module procedure alltoall_sink
  end interface p_alltoall
!==============================================================================
contains
!==============================================================================
  subroutine send_mcol (mc, dest, tslot)
  type(t_mcols)      ,intent(in) :: mc
  integer            ,intent(in) :: dest
  integer            ,intent(in) :: tslot     ! time slot
    !---------------------------------------------
    ! include interfaces for external send routine
    !---------------------------------------------
    interface p_send
      subroutine p_send_derivedtype4 (buffer, count, dest, tag, comm)
      use mo_t_obs, only: t_mcol
      type(t_mcol) ,intent(in)           :: buffer(*) ! variable to send
      integer      ,intent(in)           :: count     ! len(byte) of variable
      integer      ,intent(in)           :: dest      ! destination processor
      integer      ,intent(in)           :: tag       ! tag
      integer      ,intent(in)           :: comm      ! communicator
      end subroutine p_send_derivedtype4
    end interface p_send

    integer                    :: cnt
    integer                    :: n
    integer                    :: ncol
    type (t_mcol) ,pointer     :: mcol (:)
    logical       ,allocatable :: mask (:)

    ncol = mc% n
    allocate (mask(ncol))
    mcol => mc% c(1:ncol)
    mask = (mcol% ijdtp(4) == tslot .and. mcol% ijdtp(5) == dest)
    n = count (mask)
    call p_send (n, dest, p_tag=1)
    if (n>0) then
      cnt = n * size(transfer(mcol(1),(/' '/)))
      call p_send (pack(mcol,mask), cnt, dest, 1, dace% comm)
    endif
    deallocate (mask)
  end subroutine send_mcol
!------------------------------------------------------------------------------
  subroutine recv_mcol (mcol, src)
  type(t_mcol)      ,pointer     :: mcol (:)
  integer           ,intent(in)  :: src
    !---------------------------------------------
    ! include interfaces for external recv routine
    !---------------------------------------------
    interface p_recv
      subroutine p_recv_derivedtype4 (buffer, count, src, tag, comm)
      use mo_t_obs, only: t_mcol
      type(t_mcol) ,intent(out)          :: buffer (*) ! variable to recv
      integer      ,intent(in)           :: count      ! len(byte) of variable
      integer      ,intent(in)           :: src        ! source processor el.
      integer      ,intent(in)           :: tag        ! tag
      integer      ,intent(in)           :: comm       ! communicator
      end subroutine p_recv_derivedtype4
    end interface p_recv

    integer :: n
    integer :: cnt

    call p_recv (n, src, p_tag=1)
    if (associated (mcol)) deallocate (mcol)
    allocate (mcol(n))
    if (n>0) then
      cnt = n * size(transfer(mcol(1),(/' '/)))
      call p_recv (mcol, cnt, src, 1, dace% comm)
    endif
  end subroutine recv_mcol

!==============================================================================
! Scatter array of type t_mcol.
! t_mcol is used by the interpolation operator to describe the model columns
! (i.e. indices i,j,diamond and pe) required for interpolation to the location
! of an observation (or destination grid). This routine distributes this
! information from the observation pe to the model grid pes so that the
! required columns may be collected.
!
! The implementation using MPI_scatterv is faster,
! but we provide a fallback implementation using explicit send/recv for
! easier debugging of problems with "real-existing" MPI implementations.
!------------------------------------------------------------------------------
#define USE_MPI_SCATTER
#if defined (USE_MPI_SCATTER)
  !----------------------------------
  ! Implementation using MPI_scatterv
  !----------------------------------
  subroutine scatter_mcol (mc, mcol, tslot)
  type(t_mcols) ,intent(in) :: mc      ! info at observation pe
  type(t_mcol)  ,pointer    :: mcol(:) ! info scattered to model pes
  integer       ,intent(in) :: tslot   ! time slot
    !----------------
    ! local variables
    !----------------
    integer                   :: ncol, pe
    integer                   :: ni(0:dace% npe-1)
    integer                   :: in(0:dace% npe-1)
    integer                   :: i
    type(t_mcol) ,allocatable :: mcols(:)

    !--------------------------------
    ! Check if box handled by this PE
    !--------------------------------
    if (mc% pe == dace% pe) then

       !---------------------------
       ! Order cols by receiving PE
       !  keep order within same PE
       !---------------------------
       do pe = 0, dace% npe-1
         ni(pe) = count (mc% c(1:mc% n)% ijdtp(5) == pe   &
                   .and. mc% c(1:mc% n)% ijdtp(4) == tslot)
       end do
       in (0) = 0
       do pe = 1, dace% npe-1
         in (pe) = in (pe-1) + ni (pe-1)
       end do

       ncol = sum (ni)
       allocate (mcols(ncol))

       do i = 1, mc% n
         if  (mc% c(i)% ijdtp(4) /= tslot) cycle
         pe = mc% c(i)% ijdtp(5)
         in (pe) = in(pe) + 1
         mcols (in (pe)) = mc% c(i)
       end do

       if (any (mcols(1:ncol-1)% ijdtp(5) > mcols(2:ncol)% ijdtp(5))) then
          write(0,*) 'ncol       :',ncol
          write(0,*) 'mc% pe     :',mc% c(1:ncol)% ijdtp(5)
          write(0,*) 'mcols% pe  :',mcols(1:ncol)% ijdtp(5)
          write(0,*) 'sorted     :',(mcols(1:ncol-1)% ijdtp(5) &
                                  <= mcols(2:ncol  )% ijdtp(5) )
          call finish ("scatter_mcol","INTERNAL ERROR: bad pe ordering!")
       end if
       ! Display communication pattern for potential optimizations
       !write(0,*) "scatter_mcol: pe, ni =", dace% pe, ni
    else
       allocate (mcols(0))              ! Needed on other PEs
    end if
    call p_bcast (ni, mc% pe)           ! Might use p_scatter instead

    if (associated (mcol)) deallocate (mcol)
    allocate (mcol(ni(dace% pe)))
    call p_scatter_mcol (mcols, mcol, sendcounts=ni, root=mc% pe)

  end subroutine scatter_mcol
!==============================================================================
#define DERIVED  type(t_mcol)
#undef  MPI_TYPE
#define p_scatter_DERIVED p_scatter_mcol
#include "p_scatter_derived.incf"
#undef  DERIVED
#undef  MPI_TYPE
!==============================================================================
#else  /* !USE_MPI_SCATTER */
  !----------------------------------------
  ! Implementation using explicit send/recv
  !----------------------------------------
  subroutine scatter_mcol (mc, mcol, tslot)
  type(t_mcols) ,intent(in) :: mc      ! info at observation pe
  type(t_mcol)  ,pointer    :: mcol(:) ! info scattered to model pes
  integer       ,intent(in) :: tslot   ! time slot

    !----------------
    ! local variables
    !----------------
    integer              :: j, n, ncol
    logical, allocatable :: mask(:)

    do j=0,dace% npe-1
       !----------------------------------
       ! send if box is handled by this PE
       !----------------------------------
       if (mc% pe == dace% pe) then
          if (j /= dace% pe) then
             call send_mcol (mc, j, tslot)
          else if (j == dace% pe) then
             !---------------------------
             ! copy if sender == receiver
             !---------------------------
             ncol = mc% n
             allocate (mask(ncol))
             mask = (mc% c(1:ncol)% ijdtp(5) == dace% pe &
               .and. mc% c(1:ncol)% ijdtp(4) == tslot)
             n = count (mask)
             allocate (mcol(n))
             mcol = pack (mc% c(1:ncol), mask)
             deallocate (mask)
          endif
       endif
       !-------------------------------------
       ! receive if box is handled by partner
       !-------------------------------------
       if (j == mc% pe .and. j /= dace% pe) then
          call recv_mcol (mcol, j)
       endif
    end do
  end subroutine scatter_mcol
#endif /* !USE_MPI_SCATTER */
!==============================================================================
  subroutine alltoall_mcol (mc, mcols, tslot)
    type(t_mcols) ,intent(in)  :: mc   (:)  ! model column descriptor sets/boxes
    type(p_mcol)  ,intent(out) :: mcols(:)  ! info scattered to model pes
    integer       ,intent(in)  :: tslot     ! time slot
    !----------------
    ! local variables
    !----------------
    integer                   :: ib

#ifdef ALLTOALL_MCOL_ORIGINAL

    !----------------
    ! Loop over boxes
    !----------------
    do ib = 1, size (mc)
      nullify (mcols(ib)%p)
FTRACE_BEGIN("alltoall_mcol:scatter_mcol")
      call scatter_mcol (mc(ib), mcols(ib)% p, tslot)
FTRACE_END  ("alltoall_mcol:scatter_mcol")
    end do

#else

    integer                   :: i, k                     ! Loop indices
    integer                   :: ncol, pe
    integer                   :: scnt(0:dace% npe-1)      ! Sendcounts
    integer                   :: rcnt(0:dace% npe-1)      ! Recvcounts
    integer                   :: soff(0:dace% npe-1)      ! Send displacements
!   integer                   :: roff(0:dace% npe-1)      ! Recv displacements
    integer                   :: idx (0:dace% npe-1)      ! Aux. index
    integer,      allocatable :: ixb (:)                  ! indices in boxes
    integer,      allocatable :: ncb (:)                  ! # cols in box
    integer,      allocatable :: map (:,:)                ! Send/recv pattern
    integer,      allocatable :: sbuf(:), rbuf(:)         ! temp. buffers
    type(t_mcol), allocatable :: scols(:)                 ! Send buffer
    type(t_mcol), allocatable :: rcols(:)                 ! Recv buffer

    if (size (mc) == 0) return

FTRACE_BEGIN("alltoall_mcol:prep")
    !----------------------------
    ! Loop over boxes: sendcounts
    !----------------------------
    allocate (map (size (mc),0:dace% npe-1))
    map  = 0
    scnt = 0
    do ib = 1, size (mc)
       !--------------------------------
       ! Check if box handled by this PE
       !--------------------------------
       if (mc(ib)% pe /= dace% pe) cycle
       do pe = 0, dace% npe-1
          map(ib,pe) = count (mc(ib)% c(1:mc(ib)% n)% ijdtp(5) == pe   &
                        .and. mc(ib)% c(1:mc(ib)% n)% ijdtp(4) == tslot)
          scnt  (pe) = scnt(pe) + map(ib,pe)
       end do ! pe
    end do ! ib

    soff(0) = 0
    do pe = 1, dace% npe-1
       soff(pe) = soff(pe-1) + scnt(pe-1)
    end do ! pe

    !--------------------
    ! Prepare send buffer
    !--------------------
    ncol = sum (scnt)
    allocate (scols(ncol))

    !---------------------------
    ! Pack cols for all boxes:
    ! order cols by receiving PE
    !  keep order within same PE
    !---------------------------
    idx(:) = soff(:)
    do ib = 1, size (mc)
       if (mc(ib)% pe /= dace% pe) cycle
       do i = 1, mc(ib)% n
          if  (mc(ib)% c(i)% ijdtp(4) /= tslot) cycle
          pe = mc(ib)% c(i)% ijdtp(5)
          idx(pe) = idx(pe) + 1
          scols(idx(pe)) = mc(ib)% c(i)
       end do
    end do ! ib

    if (any (idx /= soff+scnt)) &
         call finish ("alltoall_mcols","idx /= soff+scnt")

    if (any (scols(1:ncol-1)% ijdtp(5) > scols(2:ncol)% ijdtp(5))) then
       write(0,*) 'ncol       :',ncol
       write(0,*) 'mc% pe     :',mc(ib)% c(1:ncol)% ijdtp(5)
       write(0,*) 'mcols% pe  :',scols(1:ncol)% ijdtp(5)
       write(0,*) 'sorted     :',(scols(1:ncol-1)% ijdtp(5) &
                               <= scols(2:ncol  )% ijdtp(5) )
       call finish ("alltoall_mcols","INTERNAL ERROR: bad pe ordering!")
    end if

FTRACE_END  ("alltoall_mcol:prep")
FTRACE_BEGIN("alltoall_mcol:alltoall")
    !---------------------------------------------------
    ! Communicate mapping for unpacking on receiving pes
    !---------------------------------------------------
    allocate (sbuf(size (mc) * dace% npe))
    allocate (rbuf(size (mc) * dace% npe))
    sbuf = reshape (map,  shape (sbuf))
    call p_alltoall (sbuf, rbuf)
    map  = reshape (rbuf, shape (map))
    deallocate (sbuf, rbuf)

    call p_alltoall (scnt, rcnt)

    ncol = sum (rcnt)
    allocate (rcols(ncol))
    call p_alltoall_mcol (scols, rcols, sendcounts=scnt, recvcounts=rcnt)
    deallocate (scols)

FTRACE_END  ("alltoall_mcol:alltoall")
FTRACE_BEGIN("alltoall_mcol:unpack")
    !--------------------------
    ! Prepare unpacking of cols
    !--------------------------
    allocate (ncb(size (mc)))
    ncb(:) = sum (map, dim=2)
    do ib = 1, size (mc)
       allocate (mcols(ib)% p(ncb(ib)))
    end do ! ib

    !------------
    ! Unpack cols
    !------------
    allocate (ixb (size (mc)))
    ixb(:) = 0
    k = 0
    do pe = 0, dace% npe-1
       do ib = 1, size (mc)
          do i = 1, map(ib,pe)
             k = k + 1
             ixb(ib) = ixb(ib) + 1
             mcols(ib)% p(ixb(ib)) = rcols(k)
          end do
       end do ! ib
    end do ! pe

    if (any (ixb /= ncb)) then
       write(0,*) 'ncb :',ncb
       write(0,*) 'ixb :',ixb
       call finish ("alltoall_mcols","INTERNAL ERROR: ixb /= ncb!")
    end if

    deallocate (rcols)

FTRACE_END  ("alltoall_mcol:unpack")

#endif
  end subroutine alltoall_mcol
!==============================================================================
  subroutine send_obs (obs, dest)
  type(t_obs)       ,intent(in) :: obs
  integer           ,intent(in) :: dest
    !---------------------------------------------
    ! include interfaces for external send routine
    !---------------------------------------------
    interface p_send
      subroutine p_send_derivedtype (buffer, count, dest, tag, comm)
      use mo_t_obs, only: t_obs
      type(t_obs) ,intent(in)           :: buffer    ! variable to send
      integer     ,intent(in)           :: count     ! len(byte) of variable
      integer     ,intent(in)           :: dest      ! destination processor
      integer     ,intent(in)           :: tag       ! tag
      integer     ,intent(in)           :: comm      ! communicator
      end subroutine p_send_derivedtype
    end interface p_send

    interface p_send
      subroutine p_send_derivedtype2 (buffer, count, dest, tag, comm)
      use mo_t_obs, only: t_spot
      type(t_spot),intent(in)           :: buffer(*) ! variable to send
      integer     ,intent(in)           :: count     ! len(byte) of variable
      integer     ,intent(in)           :: dest      ! destination processor
      integer     ,intent(in)           :: tag       ! tag
      integer     ,intent(in)           :: comm      ! communicator
      end subroutine p_send_derivedtype2
    end interface p_send

    interface p_send
      subroutine p_send_derivedtype3 (buffer, count, dest, tag, comm)
      use mo_t_obs, only: t_datum
      type(t_datum),intent(in)           :: buffer(*) ! variable to send
      integer      ,intent(in)           :: count     ! len(byte) of variable
      integer      ,intent(in)           :: dest      ! destination processor
      integer      ,intent(in)           :: tag       ! tag
      integer      ,intent(in)           :: comm      ! communicator
      end subroutine p_send_derivedtype3
    end interface p_send
    !---------------------------------------------------------
    ! return if destination and source processors are the same
    !---------------------------------------------------------
    if (dest == dace% pe) return
    !--------------------------------------------------
    ! determine byte size of derived data types to send
    !--------------------------------------------------
    if (obs_bytes == 0) call set_byte_size
    !-----
    ! send
    !-----
    call p_send (obs       ,obs_bytes              ,dest, 1, dace% comm)
    call p_send (obs% spot ,spot_bytes*obs% n_spot ,dest, 2, dace% comm)
    if (associated (obs% body  )) &
       call p_send (obs% body, body_bytes*obs% n_obs,dest,3, dace% comm)
    if (associated (obs% varno )) &
       call p_send (obs% varno (:obs%n_obs) ,dest,p_tag=4)
    if (associated (obs% olev  )) &
       call p_send (obs% olev  (:obs%n_obs) ,dest,p_tag=6)
    if (associated (obs% bger  )) &
       call p_send (obs% bger  (:obs%n_obs) ,dest,p_tag=10)
    if (associated (obs% t_int )) &
       call p_send (obs% t_int (:obs%n_int) ,dest,p_tag=11)
    if (associated (obs% lev   )) &
       call p_send (obs% lev   (:obs%n_int) ,dest,p_tag=12)
    if (associated (obs% bgeri )) &
       call p_send (obs% bgeri (:obs%n_int) ,dest,p_tag=13)
    if (associated (obs% bgi   )) &
       call p_send (obs% bgi   (:obs%n_int) ,dest,p_tag=14)
    if (associated (obs% par   )) &
       call p_send (obs% par   (:obs%n_par) ,dest,p_tag=16)
    if (associated (obs% s_vqc )) &
       call p_send (obs% s_vqc (:obs%n_obs) ,dest,p_tag=17)
    if (associated (obs% f_vqc )) &
       call p_send (obs% f_vqc (:obs%n_obs) ,dest,p_tag=18)
    if (associated (obs% bcpred)) &
       call p_send (obs% bcpred(:obs%n_bcp) ,dest,p_tag=19)
    if (associated (obs% sink  )) &
       call p_send (obs% sink  (:obs%n_sink),dest  ,tag=20)
  end subroutine send_obs
!------------------------------------------------------------------------------
  subroutine send_obss (obs, dest)
  type(t_obs)       ,intent(in) :: obs  (:)
  integer           ,intent(in) :: dest (:)
    integer :: i
    do i=1,size(obs)
      call send_obs (obs(i), dest(i))
    end do
  end subroutine send_obss
!==============================================================================
  subroutine recv_obs (obs, src)
  type(t_obs)       ,intent(inout) :: obs
  integer           ,intent(in)    :: src
    integer :: i
    !---------------------------------------------
    ! include interfaces for external recv routine
    !---------------------------------------------
    interface p_recv
      subroutine p_recv_derivedtype (buffer, count, dest, tag, comm)
      use mo_t_obs, only: t_obs
      type(t_obs) ,intent(out)          :: buffer    ! variable to recv
      integer     ,intent(in)           :: count     ! len(byte) of variable
      integer     ,intent(in)           :: dest      ! destination processor
      integer     ,intent(in)           :: tag       ! tag
      integer     ,intent(in)           :: comm      ! communicator
      end subroutine p_recv_derivedtype
    end interface p_recv

    interface p_recv
      subroutine p_recv_derivedtype2 (buffer, count, dest, tag, comm)
      use mo_t_obs, only: t_spot
      type(t_spot),intent(inout)        :: buffer(*) ! variable to recv
      integer     ,intent(in)           :: count     ! len(byte) of variable
      integer     ,intent(in)           :: dest      ! destination processor
      integer     ,intent(in)           :: tag       ! tag
      integer     ,intent(in)           :: comm      ! communicator
      end subroutine p_recv_derivedtype2
    end interface p_recv

    interface p_recv
      subroutine p_recv_derivedtype3 (buffer, count, dest, tag, comm)
      use mo_t_obs, only: t_datum
      type(t_datum),intent(inout)        :: buffer(*) ! variable to recv
      integer      ,intent(in)           :: count     ! len(byte) of variable
      integer      ,intent(in)           :: dest      ! destination processor
      integer      ,intent(in)           :: tag       ! tag
      integer      ,intent(in)           :: comm      ! communicator
      end subroutine p_recv_derivedtype3
    end interface p_recv

    if (src == dace% pe) return
    !--------------------------------------------------
    ! determine byte size of derived data types to recv
    !--------------------------------------------------
    if (obs_bytes == 0) call set_byte_size
    !------------------
    ! receive container
    !------------------
    call p_recv (obs ,obs_bytes ,src , 1, dace% comm)
    call alloc_recvd_obs (obs)

    !----------------------------------------------
    ! receive spots, no pointer components so far
    !----------------------------------------------
    call p_recv (obs% spot,    spot_bytes*obs% n_spot ,src, 2, dace% comm)
    do i=1,size(obs% spot)
      nullify (obs% spot(i)% imcol)
      nullify (obs% spot(i)% pe_dest)
    end do

    !-----------------------------
    ! receive remaining components
    !-----------------------------
    if (associated (obs% body  )) &
       call p_recv (obs% body, body_bytes*obs% n_obs  ,src, 3, dace% comm)
    if (associated (obs% varno )) &
       call p_recv (obs% varno  ,src ,p_tag= 4 )
    if (associated (obs% olev  )) &
       call p_recv (obs% olev   ,src ,p_tag= 6 )
    if (associated (obs% bger  )) &
       call p_recv (obs% bger   ,src ,p_tag=10 )
    if (associated (obs% t_int )) &
       call p_recv (obs% t_int  ,src ,p_tag=11 )
    if (associated (obs% lev   )) &
       call p_recv (obs% lev    ,src ,p_tag=12 )
    if (associated (obs% bgeri )) &
       call p_recv (obs% bgeri  ,src ,p_tag=13 )
    if (associated (obs% bgi   )) &
       call p_recv (obs% bgi    ,src ,p_tag=14 )
    if (associated (obs% par   )) &
       call p_recv (obs% par    ,src ,p_tag=16 )
    if (associated (obs% s_vqc )) &
       call p_recv (obs% s_vqc  ,src ,p_tag=17 )
    if (associated (obs% f_vqc )) &
       call p_recv (obs% f_vqc  ,src ,p_tag=18 )
    if (associated (obs% bcpred)) &
       call p_recv (obs% bcpred ,src ,p_tag=19 )
    if (associated (obs% sink))   &
       call p_recv (obs% sink   ,src   ,tag=20 )

    nullify (obs% levs)
    obs% n_lev         = 0
    obs% spot(:)% l% i = 0
    obs% spot(:)% l% n = 0

  end subroutine recv_obs
!------------------------------------------------------------------------------
  subroutine alloc_recvd_obs (obs)
  type (t_obs) ,intent(inout) :: obs
  !---------------------------------------------------
  ! allocate components of received data of type t_obs
  !---------------------------------------------------
                                  allocate (obs% spot  ( obs% n_spot))
    if (associated (obs% body  )) allocate (obs% body   (obs% n_obs ))
    if (associated (obs% varno )) allocate (obs% varno  (obs% n_obs ))
    if (associated (obs% olev  )) allocate (obs% olev   (obs% n_obs ))
    if (associated (obs% bger  )) allocate (obs% bger   (obs% n_obs ))
    if (associated (obs% t_int )) allocate (obs% t_int  (obs% n_int ))
    if (associated (obs% lev   )) allocate (obs% lev    (obs% n_int ))
    if (associated (obs% bgeri )) allocate (obs% bgeri  (obs% n_int ))
    if (associated (obs% bgi   )) allocate (obs% bgi    (obs% n_int ))
    if (associated (obs% par   )) allocate (obs% par    (obs% n_par ))
    if (associated (obs% s_vqc )) allocate (obs% s_vqc  (obs% n_obs ))
    if (associated (obs% f_vqc )) allocate (obs% f_vqc  (obs% n_obs ))
    if (associated (obs% bcpred)) allocate (obs% bcpred (obs% n_bcp ))
    if (associated (obs% sink  )) allocate (obs% sink   (obs% n_sink))
    nullify (obs% mc% c)
  end subroutine alloc_recvd_obs
!------------------------------------------------------------------------------
  subroutine recv_obss (obs, src)
  type(t_obs)       ,intent(inout) :: obs (:)
  integer           ,intent(in)    :: src (:)
    integer :: i
    do i=1,size(obs)
      call recv_obs (obs(i), src(i))
    end do
  end subroutine recv_obss
!==============================================================================
  subroutine bcast_obs (obs, src)
  type(t_obs)       ,intent(inout) :: obs
  integer           ,intent(in)    :: src
    integer :: i
    !----------------------------------------------
    ! include interfaces for external bcast routine
    !----------------------------------------------
    interface p_bcast
      subroutine p_bcast_derivedtype (buffer, count, source, comm)
      use mo_t_obs, only: t_obs
      type(t_obs)  ,intent(inout)        :: buffer    ! variable to bcast
      integer      ,intent(in)           :: count     ! len(byte) of variable
      integer      ,intent(in)           :: source    ! source processor index
      integer      ,intent(in)           :: comm      ! communicator
      end subroutine p_bcast_derivedtype
    end interface p_bcast

    interface p_bcast
      subroutine p_bcast_derivedtype2 (buffer, count, source, comm)
      use mo_t_obs, only: t_spot
      type(t_spot) ,intent(inout)        :: buffer(*) ! variable to bcast
      integer      ,intent(in)           :: count     ! len(byte) of variable
      integer      ,intent(in)           :: source    ! source processor index
      integer      ,intent(in)           :: comm      ! communicator
      end subroutine p_bcast_derivedtype2
    end interface p_bcast

    interface p_bcast
      subroutine p_bcast_derivedtype3 (buffer, count, source, comm)
      use mo_t_obs, only: t_datum
      type(t_datum),intent(inout)        :: buffer(*) ! variable to bcast
      integer      ,intent(in)           :: count     ! len(byte) of variable
      integer      ,intent(in)           :: source    ! source processor index
      integer      ,intent(in)           :: comm      ! communicator
      end subroutine p_bcast_derivedtype3
    end interface p_bcast

    !--------------------------------------------------
    ! determine byte size of derived data types to recv
    !--------------------------------------------------
    if (obs_bytes == 0) call set_byte_size
    !--------------------
    ! broadcast container
    !--------------------
    call p_bcast (obs ,obs_bytes ,src , dace% comm)

    if (src /= dace% pe) call alloc_recvd_obs (obs)

    !----------------------------------------------
    ! broadcast spots, no pointer components so far
    !----------------------------------------------
    call p_bcast (obs% spot,    spot_bytes*obs% n_spot ,src, dace% comm)
    if (src /= dace% pe) then
      do i=1,size(obs% spot)
        nullify (obs% spot(i)% imcol)
        nullify (obs% spot(i)% pe_dest)
      end do
    endif

    !-------------------------------
    ! broadcast remaining components
    !-------------------------------
    if (associated  (obs% body ))   &
       call p_bcast (obs% body, body_bytes*obs% n_obs  ,src, dace% comm)
    if (associated  (obs% varno))   &
       call p_bcast (obs% varno ,src)
    if (associated  (obs% olev ))   &
       call p_bcast (obs% olev  ,src)
    if (associated  (obs% bger ))   &
       call p_bcast (obs% bger  ,src)
    if (associated  (obs% t_int))   &
       call p_bcast (obs% t_int ,src)
    if (associated  (obs% lev  ))   &
       call p_bcast (obs% lev   ,src)
    if (associated  (obs% bgeri))   &
       call p_bcast (obs% bgeri ,src)
    if (associated  (obs% bgi  ))   &
       call p_bcast (obs% bgi   ,src)
    if (associated  (obs% par  ))   &
       call p_bcast (obs% par   ,src)
    if (associated  (obs% s_vqc))   &
       call p_bcast (obs% s_vqc ,src)
    if (associated  (obs% f_vqc))   &
       call p_bcast (obs% f_vqc ,src)
    if (associated  (obs% bcpred))  &
       call p_bcast (obs% bcpred,src)
    if (associated  (obs% sink  ))  &
       call p_bcast (obs% sink  ,src)

    if (src == dace% pe) return
    nullify (obs% levs)
    obs% n_lev         = 0
    obs% spot(:)% l% i = 0
    obs% spot(:)% l% n = 0

  end subroutine bcast_obs
!------------------------------------------------------------------------------
  subroutine bcast_obss (obs, src)
  type(t_obs)       ,intent(inout) :: obs (:)
  integer           ,intent(in)    :: src (:)
    integer :: i
    do i=1,size(obs)
      call bcast_obs (obs(i), src(i))
    end do
  end subroutine bcast_obss
!==============================================================================
  subroutine alltoall_obs (sendbuf, recvbuf, ibx)
  type(t_obs)       ,INTENT(inout) :: sendbuf     ! send buffer
  type(t_obs)       ,INTENT(out)   :: recvbuf     ! receive buffer
  type(t_box)       ,INTENT(in)    :: ibx (:)     ! where to send each report
  !---------------------------------------------------------
  ! MPI_ALLTOALL for observation data
  !
  ! Used after the first guess step to redistribute
  ! observations for the analysis step.
  !
  ! The destination PE for each report is given in ibx% pe .
  ! Observations on each sending PE must be sorted
  ! according to the destination PE previously.
  !---------------------------------------------------------
    integer          :: ib                 ! index in send-count array
    integer          :: pe1                ! index: processor element + 1
    integer          :: i                  ! index
    integer          :: sndcnt (dace% npe) ! MPI send-count array
    integer          :: rcvcnt (dace% npe) ! MPI receive-count array
    integer          :: rsize, ssize       ! sizes of send and receive buffers
    logical          :: asso               ! association status of components
    logical          :: error              ! error flag

    !-------------------------------------------
    ! set up sendcount array, consistency checks
    !-------------------------------------------
    if (size(ibx) /= sendbuf% n_spot) &
      call finish ('alltoall_obs','size(ibx) /= n_spot')
    if (maxval(ibx% pe) > dace% npe -1 .or. &
        minval(ibx% pe) < 0           )    &
      call finish ('alltoall_obs','ibx% pe out of range')
    sndcnt = 0
    ib     = 1
    error  = .false.
    do i = 1, sendbuf% n_spot
      pe1 = ibx(i)% pe + 1
      if (pe1 > ib) ib = pe1
      if (pe1 < ib) error = .true.
      sndcnt(ib) = sndcnt(ib) + 1
    end do
    if (error) call finish ('alltoall_obs','ibx% pe not in consecutive order')
    recvbuf% n_time = sendbuf% n_time
    !-------------------------
    ! send/recv component spot
    !-------------------------
    call p_alltoall (sndcnt, rcvcnt)
    rsize = sum(rcvcnt)
    allocate (recvbuf% spot (rsize))
    call p_alltoall (sendbuf% spot, recvbuf% spot,       &
                     sendcounts=sndcnt, recvcounts=rcvcnt)
    recvbuf% n_spot = rsize

    !------------------------------------------
    ! send/recv components in observation space
    !------------------------------------------
    sndcnt = 0
    do i = 1, size(ibx)
      ib = ibx(i)% pe + 1
      sndcnt (ib) = sndcnt (ib) + sendbuf% spot(i)% o% n
    end do
    ssize = sum(sndcnt)
    rsize = sum(recvbuf% spot(:)% o% n)
    recvbuf% n_obs = rsize

    asso = p_or (associated(sendbuf% body) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% body))) allocate(sendbuf% body(0))
      allocate(recvbuf% body(rsize))
      call p_alltoall (sendbuf% body, recvbuf% body, sendcounts=sndcnt)
    endif

    asso = p_or (associated(sendbuf% varno) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% varno))) allocate(sendbuf% varno(0))
      allocate(recvbuf% varno(rsize))
      call p_alltoall (sendbuf% varno, recvbuf% varno, sendcounts=sndcnt)
    endif

    asso = p_or (associated(sendbuf% olev) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% olev))) allocate(sendbuf% olev(0))
      allocate(recvbuf% olev(rsize))
      call p_alltoall (sendbuf% olev, recvbuf% olev, sendcounts=sndcnt)
    endif

    asso = p_or (associated(sendbuf% bger) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% bger))) allocate(sendbuf% bger(0))
      allocate(recvbuf% bger(rsize))
      call p_alltoall (sendbuf% bger, recvbuf% bger, sendcounts=sndcnt)
    endif

    asso = p_or (associated(sendbuf% s_vqc) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% s_vqc))) allocate(sendbuf% s_vqc(0))
      allocate(recvbuf% s_vqc(rsize))
      call p_alltoall (sendbuf% s_vqc, recvbuf% s_vqc, sendcounts=sndcnt)
    endif

    asso = p_or (associated(sendbuf% f_vqc) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% f_vqc))) allocate(sendbuf% f_vqc(0))
      allocate(recvbuf% f_vqc(rsize))
      call p_alltoall (sendbuf% f_vqc, recvbuf% f_vqc, sendcounts=sndcnt)
    endif

    !--------------------------------------------
    ! send/recv components in interpolation space
    !--------------------------------------------
    sndcnt = 0
    do i = 1, size(ibx)
      ib = ibx(i)% pe + 1
      sndcnt (ib) = sndcnt (ib) + sendbuf% spot(i)% i% n
    end do
    ssize = sum(sndcnt)
    rsize = sum(recvbuf% spot(:)% i% n)
    recvbuf% n_int = rsize

    asso = p_or (associated(sendbuf% t_int) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% t_int))) allocate(sendbuf% t_int(0))
      allocate(recvbuf% t_int(rsize))
      call p_alltoall (sendbuf% t_int, recvbuf% t_int, sendcounts=sndcnt)
    endif

    asso = p_or (associated(sendbuf% lev) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% lev))) allocate(sendbuf% lev(0))
      allocate(recvbuf% lev(rsize))
      call p_alltoall (sendbuf% lev, recvbuf% lev, sendcounts=sndcnt)
    endif

    asso = p_or (associated(sendbuf% bgeri) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% bgeri))) allocate(sendbuf% bgeri(0))
      allocate(recvbuf% bgeri(rsize))
      call p_alltoall (sendbuf% bgeri, recvbuf% bgeri, sendcounts=sndcnt)
    endif

    asso = p_or (associated(sendbuf% bgi  ) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% bgi  ))) allocate(sendbuf% bgi  (0))
      allocate(recvbuf% bgi  (rsize))
      call p_alltoall (sendbuf% bgi  , recvbuf% bgi  , sendcounts=sndcnt)
    endif

    !-------------------------------------
    ! send/recv bias correction predictors
    !-------------------------------------
    sndcnt = 0
    do i = 1, size(ibx)
      ib = ibx(i)% pe + 1
      sndcnt (ib) = sndcnt (ib) + sendbuf% spot(i)% bcp% n
    end do
    ssize = sum(sndcnt)
    rsize = sum(recvbuf% spot(:)% bcp% n)
    recvbuf% n_bcp = rsize

    asso = p_or (associated(sendbuf% bcpred) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% bcpred))) allocate(sendbuf% bcpred(0))
      allocate(recvbuf% bcpred(rsize))
      call p_alltoall (sendbuf% bcpred, recvbuf% bcpred, sendcounts=sndcnt)
    endif

    !----------------------------------------
    ! send/recv components in parameter space
    !----------------------------------------
    sndcnt = 0
    do i = 1, size(ibx)
      ib = ibx(i)% pe + 1
      sndcnt (ib) = sndcnt (ib) + sendbuf% spot(i)% p% n
    end do
    ssize = sum(sndcnt)
    rsize = sum(recvbuf% spot(:)% p% n)
    recvbuf% n_par = rsize

    asso = p_or (associated(sendbuf% par) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% par))) allocate(sendbuf% par(0))
      allocate(recvbuf% par(rsize))
      call p_alltoall (sendbuf% par, recvbuf% par, sendcounts=sndcnt)
    endif

    !----------------------------------
    ! send/recv sink variable meta data
    !----------------------------------
    sndcnt = 0
    do i = 1, size(ibx)
      ib = ibx(i)% pe + 1
      sndcnt (ib) = sndcnt (ib) + sendbuf% spot(i)% d% n
    end do
    ssize = sum(sndcnt)
    rsize = sum(recvbuf% spot(:)% d% n)
    recvbuf% n_sink = rsize

    asso = p_or (associated(sendbuf% sink) .and. ssize > 0)
    if (asso) then
      if (.not.(associated(sendbuf% sink))) allocate(sendbuf% sink(0))
      allocate(recvbuf% sink(rsize))
      call p_alltoall (sendbuf% sink, recvbuf% sink, sendcounts=sndcnt)
    endif

    nullify (recvbuf% levs)
    recvbuf% n_lev = 0

    !------------------------------
    ! fix indices in component spot
    !------------------------------
    if (recvbuf% n_spot > 0) then
      recvbuf% spot(1)% o%   i =  0
      recvbuf% spot(1)% i%   i =  0
      recvbuf% spot(1)% p%   i =  0
      recvbuf% spot(1)% d%   i =  0
!     recvbuf% spot(1)% l%   i =  0
      recvbuf% spot(1)% bcp% i =  0
    endif
    do i = 1, recvbuf% n_spot-1
      recvbuf% spot(i+1)% o%  i = recvbuf% spot(i)% o%   i &
                                + recvbuf% spot(i)% o%   n
      recvbuf% spot(i+1)% i%  i = recvbuf% spot(i)% i%   i &
                                + recvbuf% spot(i)% i%   n
      recvbuf% spot(i+1)% p%  i = recvbuf% spot(i)% p%   i &
                                + recvbuf% spot(i)% p%   n
      recvbuf% spot(i+1)% d%  i = recvbuf% spot(i)% d%   i &
                                + recvbuf% spot(i)% d%   n
!     recvbuf% spot(i+1)% l%  i = recvbuf% spot(i)% l%   i &
!                               + recvbuf% spot(i)% l%   n
      recvbuf% spot(i+1)% bcp%i = recvbuf% spot(i)% bcp% i &
                                + recvbuf% spot(i)% bcp% n
    end do

    recvbuf% spot(:)% l% i = 0
    recvbuf% spot(:)% l% n = 0

  end subroutine alltoall_obs
!==============================================================================
#define DERIVED type(t_obs)
#undef MPI_TYPE
#define p_alltoall_DERIVED p_alltoall_obs
#include "p_alltoall_derived.incf"
!------------------------------------------------------------------------------
#undef  DERIVED
#define DERIVED type(t_datum)
#undef MPI_TYPE
#undef  p_alltoall_DERIVED
#define p_alltoall_DERIVED alltoall_datum
#include "p_alltoall_derived.incf"
!------------------------------------------------------------------------------
#undef  DERIVED
#define DERIVED type(t_spot)
#undef MPI_TYPE
#undef  p_alltoall_DERIVED
#define p_alltoall_DERIVED alltoall_spot
#include "p_alltoall_derived.incf"
!------------------------------------------------------------------------------
#undef  DERIVED
#define DERIVED type(t_sink)
#undef MPI_TYPE
#undef  p_alltoall_DERIVED
#define p_alltoall_DERIVED alltoall_sink
#include "p_alltoall_derived.incf"
!==============================================================================
#undef  DERIVED
#define DERIVED type(t_mcol)
#undef MPI_TYPE
#undef  p_alltoall_DERIVED
#define p_alltoall_DERIVED p_alltoall_mcol
#include "p_alltoall_derived.incf"
!==============================================================================
! send/recv/bcast t_sink
!-----------------------
#undef  DERIVED
#define DERIVED type(t_sink),dimension(:)
#undef  MPI_TYPE
#define VECTOR
#define p_send_DERIVED send_sink
#include "p_send.incf"
!------------------------------------------------------------------------------
#undef  DERIVED
#define DERIVED type(t_sink),dimension(:)
#undef  MPI_TYPE
#define VECTOR
#define p_recv_DERIVED recv_sink
#include "p_recv.incf"
!------------------------------------------------------------------------------
#undef  DERIVED
#define DERIVED type(t_sink),dimension(:)
#undef  MPI_TYPE
#define VECTOR
#define p_bcast_DERIVED bcast_sink
#include "p_bcast.incf"
!==============================================================================
end module mo_obs_sndrcv
