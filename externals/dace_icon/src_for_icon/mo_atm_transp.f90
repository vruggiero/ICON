!
!+ gather/scatter routines for atmospheric fields
!
MODULE mo_atm_transp
!
! Description:
!   transpositions and gather/scatter routines for atmospheric fields
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
!  Optimisations for NEC-SX
! V1_5         2009/05/25 Harald Anlauf
!  minor cleanups in gather_level
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  workaround for NEC compiler (specific interface for scatter_multi)
! V1_19        2012-04-16 Andreas Rhodin
!  Optimize communication for LETKF parallel I/O
! V1_20        2012-06-18 Andreas Rhodin
!  minor cleanup
! V1_45        2015-12-15 Harald Anlauf
!  alltoallv_multi: parallelize loop
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM  2002  original code
!------------------------------------------------------------------------------
! Platform-dependent settings, may depend on quality of MPI implementation.
! Choose between MPI_SCATTER(V)/MPI_GATHER(V) and explicit MPI_SEND/MPI_RECV:
! (MPI_SCATTER/MPI_GATHER has better performance on NEC SX).
#define USE_MPI_SCATTER
#define USE_MPI_GATHER
!------------------------------------------------------------------------------
#if defined (_FTRACE) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!------------------------------------------------------------------------------
#include "tr15581.incf"
  !-------------
  ! Modules used
  !-------------
  use mo_kind,       only : wp
  use mo_exception,  only : finish
  use mo_atm_decomp, only : t_atm_dec
  use mo_memory,     only : t_m                 ! table of atmospheric fields
!------------------------------------------------------------------------------
  use mo_mpi_dace,   only : dace,              &! dace communication info
                            p_alltoall          ! generic alltoall(v)
#if  defined (USE_MPI_SCATTER)
  use mo_mpi_dace,   only : p_scatterv          ! generic mpi_scatterv
#endif
#if  defined (USE_MPI_GATHER)
  use mo_mpi_dace,   only : p_gatherv           ! generic mpi_gatherv
#endif
#if !defined (USE_MPI_SCATTER) || !defined (USE_MPI_GATHER)
  use mo_mpi_dace,   only : p_send, p_recv      ! generic mpi_send/mpi_recv
#endif
!------------------------------------------------------------------------------
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: scatter_level
  public :: scatter_multi     ! +++ provide specific routine for NEC compiler
  public :: gather_level      ! +++ to be replaced by alltoallv_multi_w
  public :: gather_multi      ! +++ for preliminary tests
  public :: alltoallv_multi   ! scatter multi level ensemble (read)
  public :: alltoallv_multi_w ! gather  multi level ensemble (to write)

! interface scatter_level  !+++ provide specific routines for NEC compiler
!    module procedure scatter_level     ! Scatter single-level field
!    module procedure scatter_multi     ! Scatter multi -level field
! end interface

  interface gather_level
     module procedure gather_level      ! Gather single-level field
     module procedure gather_multi      ! Gather multi -level field
  end interface

contains
!==============================================================================
#ifdef  USE_MPI_SCATTER
!------------------------------------------------------------------------------
  subroutine scatter_level (xg, x, dc, source)
  real(wp)        ,pointer     :: xg (:,:,:) ! Global fields
  real(wp)        ,intent(out) ::  x (:,:,:) ! Local  fields
  type(t_atm_dec) ,intent(in)  :: dc         ! Domain decomposition information
  integer         ,intent(in)  :: source     ! Processor holding global fields

    integer               :: i, j       ! loop indices
    integer               :: pe         ! receiving pe
    integer               :: nd         ! # of diamonds
    integer,  pointer     :: counts(:)  ! sendcounts
    integer,  pointer     :: displs(:)  ! displacements
    real(wp), allocatable :: sendbuf(:) ! temporary send buffer
    real(wp), allocatable :: recvbuf(:) ! temporary recv buffer

    call setup_sendrecv (dc, counts, displs)
    nd = size (x, dim=3)
    counts(:) = counts(:) * nd
    displs(:) = displs(:) * nd

    if (counts(dace% pe) /= size (x)) then
       write(0,*) dace% pe, ":", counts(dace% pe), size (x)
       call finish ("scatter_level","counts(pe) /= size (x)")
    end if
    if (dace% pe == source) then
       if (sum (counts) /= size (xg)) then
          write(0,*) dace% pe, ":", sum (counts), size (xg)
          call finish ("scatter_level","sum (counts) /= size (xg)")
       end if
       allocate (sendbuf(sum (counts))) ! Allocate full buffer on sender
       !-----------------
       ! Fill send buffer
       !-----------------
       do j = 1, dc% nproc2
          do i = 1, dc% nproc1
             pe = dc% pe_12 (i,j)
             if (counts(pe) == 0) cycle ! Workaround NEC compiler bug
             sendbuf(1+displs(pe):counts(pe)+displs(pe)) =   &
                  reshape (xg (dc%ilim1(i-1):dc%ilim1(i)-1,  &
                               dc%ilim2(j-1):dc%ilim2(j)-1,  &
                                            :             ), &
                           shape = (/ counts(pe) /) )
          end do
       end do
    else
       allocate (sendbuf(0))            ! Allocate dummy buffer on receiver
    end if

    allocate (recvbuf(counts(dace% pe)))
    call p_scatterv (sendbuf, counts, recvbuf, root=source)
    x(:,:,:) = reshape (recvbuf, shape = shape (x))
    deallocate (counts, displs)

  end subroutine scatter_level
!------------------------------------------------------------------------------
  subroutine scatter_multi (xg, x, dc, source)
  real(wp)        ,pointer     :: xg (:,:,:,:) ! Global fields
  real(wp)        ,intent(out) ::  x (:,:,:,:) ! Local  fields
  type(t_atm_dec) ,intent(in)  :: dc           ! Domain decomposition inform.
  integer         ,intent(in)  :: source       ! Processor holding global flds.

    integer               :: i, j       ! loop indices
    integer               :: pe         ! receiving pe
    integer               :: nz         ! # of levels
    integer               :: nd         ! # of diamonds
    integer,  pointer     :: counts(:)  ! sendcounts
    integer,  pointer     :: displs(:)  ! displacements
    real(wp), allocatable :: sendbuf(:) ! temporary send buffer
    real(wp), allocatable :: recvbuf(:) ! temporary recv buffer
#ifdef __NEC__
    integer :: i_, j_, k_, l_, n3, n4, ind
#endif

    call setup_sendrecv (dc, counts, displs)
    nz = size (x, dim=3)
    nd = size (x, dim=4)
    counts(:) = counts(:) * (nz * nd)
    displs(:) = displs(:) * (nz * nd)

    if (counts(dace% pe) /= size (x)) then
       write(0,*) dace% pe, ":", counts(dace% pe), size (x)
       call finish ("scatter_multi","counts(pe) /= size (x)")
    end if
    if (dace% pe == source) then
       if (sum (counts) /= size (xg)) then
          write(0,*) dace% pe, ":", sum (counts), size (xg)
          call finish ("scatter_multi","sum (counts) /= size (xg)")
       end if
       allocate (sendbuf(sum (counts))) ! Allocate full buffer on sender
       !-----------------
       ! Fill send buffer
       !-----------------
#ifdef __NEC__
       n3 = size (xg, dim=3)
       n4 = size (xg, dim=4)
#endif
       do j = 1, dc% nproc2
          do i = 1, dc% nproc1
             pe = dc% pe_12 (i,j)

#ifndef __NEC__ /* original code */

             sendbuf(1+displs(pe):counts(pe)+displs(pe)) =   &
                  reshape (xg (dc%ilim1(i-1):dc%ilim1(i)-1,  &
                               dc%ilim2(j-1):dc%ilim2(j)-1,  &
                                            :             ,  &
                                            :             ), &
                           shape = (/ counts(pe) /) )

#else           /* workaround for slow implemenation of reshape() by NEC */

             ind = 1+displs(pe)
             do l_ = 1, n4
                do k_ = 1, n3
                   do j_ = dc%ilim2(j-1), dc%ilim2(j)-1
                      do i_ = dc%ilim1(i-1), dc%ilim1(i)-1
                         sendbuf(ind+i_-dc%ilim1(i-1)) = xg (i_,  j_, k_, l_ )
                      end do
                      ind = ind + dc%ilim1(i)-dc%ilim1(i-1)
                   end do
                end do
             end do
#endif
          end do
       end do
    else
       allocate (sendbuf(0))            ! Allocate dummy buffer on receiver
    end if

    allocate (recvbuf(counts(dace% pe)))
    call p_scatterv (sendbuf, counts, recvbuf, root=source)
    x(:,:,:,:) = reshape (recvbuf, shape = shape (x))
    deallocate (counts, displs)

  end subroutine scatter_multi
!------------------------------------------------------------------------------
#else
!------------------------------------------------------------------------------
  subroutine scatter_level (xg, x, dc, source)
  real(wp)        ,pointer     :: xg (:,:,:) ! Global fields
  real(wp)        ,intent(out) ::  x (:,:,:) ! Local  fields
  type(t_atm_dec) ,intent(in)  :: dc         ! Domain decomposition information
  integer         ,intent(in)  :: source     ! Processor holding global fields

    integer :: i,j,pe

    do i = 1, dc% nproc1
      do j = 1, dc% nproc2
        pe = dc% pe_12 (i,j)
        if (pe==source) then
          if (dace% pe == source) &
            x = xg (dc%ilim1(i-1):dc%ilim1(i)-1, &
                    dc%ilim2(j-1):dc%ilim2(j)-1, &
                                 :             )
        else
          if (dace% pe == source) then
            call p_send (xg (dc%ilim1(i-1):dc%ilim1(i)-1, &
                             dc%ilim2(j-1):dc%ilim2(j)-1, &
                                          :             ),&
                         pe, p_tag=1)
          else if (pe == dace% pe) then
            call p_recv (x, source, p_tag=1)
          endif
        endif
      end do
    end do

  end subroutine scatter_level
!------------------------------------------------------------------------------
#endif
!------------------------------------------------------------------------------
  subroutine alltoallv_multi (xg, m, dc, sources, lb3, ub3)
  real(wp)        ,pointer       :: xg (:,:,:,:) ! global fields
  type(t_m)       ,intent(inout) :: m(:)         ! local fields
  type(t_atm_dec) ,intent(in)    :: dc           ! domain decomposition info
  integer         ,intent(in)    :: sources (:)  ! PEs holding global fields
  integer,optional,intent(in)    :: lb3, ub3     ! lower/upper level bounds
  !---------------------------------------------------------------------
  ! Distribute multiple atmospheric fields (ensemble members xg) read on
  ! different PEs (sources) to all PEs (m) according to the domain
  ! decomposition info (dc). Use alltoallv MPI communication
  !---------------------------------------------------------------------
    integer               :: i, j, k    ! loop indices
    integer               :: ns         ! number of senders
    integer               :: pe         ! receiving pe
    integer               :: nz         ! # of levels
    integer               :: nd         ! # of diamonds
    integer               :: nr         ! # words received from each sender
    integer               :: l3, u3     ! lower/upper level bounds
    integer,  pointer     :: counts(:)  ! sendcounts
    integer,  pointer     :: displs(:)  ! displacements
    integer,  pointer     :: rcounts(:) ! receive counts
    real(wp), allocatable :: sendbuf(:) ! temporary send buffer
    real(wp), allocatable :: recvbuf(:) ! temporary recv buffer

    call setup_sendrecv (dc, counts, displs, rcounts, sources)
    ns = size   (sources)
    nd = size   (m(1)% ptr, dim=4)
    l3 = lbound (m(1)% ptr, dim=3); if (present (lb3)) l3 = max (l3, lb3)
    u3 = ubound (m(1)% ptr, dim=3); if (present (ub3)) u3 = min (u3, ub3)
    nz = max    (u3 - l3 + 1, 0)
    counts (:) = counts (:) * (nz * nd)
    displs (:) = displs (:) * (nz * nd)
    rcounts(:) = rcounts(:) * (nz * nd)
    nr = counts (dace% pe)

    if (nr /= size (m(1)% ptr(:,:,l3:u3,:))) then
       write(0,*) dace% pe, ":", counts(dace% pe), size (m(1)% ptr(:,:,l3:u3,:))
       call finish ("alltoallv_multi","counts(pe) /= size (x)")
    end if
    if (any(dace% pe == sources)) then
       if (sum (counts) /= size (xg)) then
          write(0,*) dace% pe, ":", sum (counts), size (xg)
          call finish ("alltoallv_multi","sum (counts) /= size (xg)")
       end if
       allocate (sendbuf(sum (counts))) ! Allocate full buffer on senders
       !-----------------
       ! Fill send buffer
       !-----------------
FTRACE_BEGIN("alltoallv_multi:fillbuf")
       do j = 1, dc% nproc2
          do i = 1, dc% nproc1
             pe = dc% pe_12 (i,j)
             if (counts(pe) == 0) cycle ! Workaround NEC compiler bug
#ifdef __NEC__ /* Hack to work around abysmal performance of RESHAPE implementation by NEC. */
#if 1
             ! Work around NEC nfort-3.1.0 issues with memory use
             block
               integer :: ii, jj, kk, dd ! auxiliary loop indices
               integer :: nx, off        ! sendbuf increment, offset
               nx  = dc%ilim1(i) - dc%ilim1(i-1)
               off = displs(pe)
               do          dd = lbound (xg,4), ubound (xg,4)
                  do       kk = lbound (xg,3), ubound (xg,3)
                     do    jj = dc%ilim2(j-1), dc%ilim2(j)-1
                        do ii =             1, nx
                           sendbuf(off+ii) = xg(dc%ilim1(i-1)-1+ii,jj,kk,dd)
                        end do
                        off = off + nx
                     end do
                  end do
               end do
               if (off /= counts(pe)+displs(pe)) &
                    call finish ("alltoallv_multi","off/=counts(pe)+displs(pe)")
             end block
#else
             sendbuf(1+displs(pe):counts(pe)+displs(pe)) =   &
                  transfer(xg (dc%ilim1(i-1):dc%ilim1(i)-1,  &
                               dc%ilim2(j-1):dc%ilim2(j)-1,  &
                                            :             ,  &
                                            :             ), &
                           size = counts(pe), mold = sendbuf )
#endif
#else
             sendbuf(1+displs(pe):counts(pe)+displs(pe)) =   &
                  reshape (xg (dc%ilim1(i-1):dc%ilim1(i)-1,  &
                               dc%ilim2(j-1):dc%ilim2(j)-1,  &
                                            :             ,  &
                                            :             ), &
                           shape = (/ counts(pe) /) )
#endif
          end do
       end do
FTRACE_END  ("alltoallv_multi:fillbuf")
    else
       counts(:) = 0
       displs(:) = 0
       allocate (sendbuf(0))            ! Allocate dummy buffer on receiver
    end if

    allocate (recvbuf(ns*nr))
    call p_alltoall (sendbuf, recvbuf, sendcounts=counts, recvcounts=rcounts)

    !-------------------------------
    ! fetch data from receive buffer
    !-------------------------------
!$omp parallel do private(i,k) schedule(static)
    do i = 1, ns
      k = (i-1)*nr
      m(i)% ptr(:,:,l3:u3,:) = reshape (recvbuf(k+1:k+nr),                  &
                                        shape=shape (m(i)% ptr(:,:,l3:u3,:)))
    end do
!$omp end parallel do

    deallocate (counts, displs, rcounts)

  end subroutine alltoallv_multi
!------------------------------------------------------------------------------
  subroutine alltoallv_multi_w (xg, m, dc, dests)
  real(wp)        _POINTER       :: xg (:,:,:,:) ! global fields (recv buffer)
  type(t_m)       ,intent(in)    :: m(:)         ! local fields
  type(t_atm_dec) ,intent(in)    :: dc           ! domain decomposition info
  integer         ,intent(in)    :: dests (:)    ! PEs holding global fields
  !--------------------------------------------------------------------------
  ! Gather multiple atmospheric fields (ensemble members xg) to be written on
  ! different PEs (dests) from all PEs (m) according to the domain
  ! decomposition info (dc). Use alltoallv MPI communication
  !--------------------------------------------------------------------------

    integer               :: i, j, k    ! loop indices
    integer               :: ns         ! number of senders
    integer               :: ms         ! number of senders (excluding 'holes')
    integer               :: pe         ! receiving pe
    integer               :: nz         ! # of levels to transfer
    integer               :: nzg, nzl   ! # of levels in xg, m%ptr
    integer               :: nd         ! # of diamonds
    integer               :: nr         ! # words received from each sender
    integer               :: kl, ku     ! lower and upper level index
    integer,  pointer     :: counts (:) ! receivecounts
    integer,  pointer     :: displs (:) ! displacements
    integer,  pointer     :: scounts(:) ! sendcounts
    real(wp), allocatable :: sendbuf(:) ! temporary send buffer
    real(wp), allocatable :: recvbuf(:) ! temporary recv buffer

    call setup_sendrecv (dc, counts, displs, scounts, dests)
    ns  = size   (dests)
    nzg = size   (xg       , dim=3)     ! Global fields (recv buffer #levels)
    nzl = size   (m(1)% ptr, dim=3)     ! Local  fields
    nd  = size   (m(1)% ptr, dim=4)
    nz  = min    (nzl, nzg)
    counts (:) = counts (:) * (nz * nd)
    displs (:) = displs (:) * (nz * nd)
    scounts(:) = scounts(:) * (nz * nd)
    nr = counts (dace% pe)

    kl = lbound (m(1)% ptr, dim=3)
    ku = ubound (m(1)% ptr, dim=3)
    if (nzl /= nzg) then
       kl = lbound (xg, dim=3)
       ku = ubound (xg, dim=3)
    end if
    if (nr /= size (m(1)% ptr(:,:,kl:ku,:))) then
       write(0,*) dace% pe, ":", counts(dace% pe), size (m(1)% ptr(:,:,kl:ku,:))
       call finish ("alltoallv_multi_w","counts(pe) /= size (x)")
    end if

    !-----------------
    ! fill send buffer
    !-----------------
    ms = count (dests >= 0)
    allocate (sendbuf(ms*nr))
    k       = 0
FTRACE_BEGIN("alltoallv_multi_w:fillbuf")
    do i = 1, ns
       if (dests(i) >= 0) then
#ifdef __NEC__ /* Hack to work around abysmal performance of RESHAPE implementation by NEC. */
          sendbuf(k+1:k+nr) = transfer(m(i)% ptr(:,:,kl:ku,:), size=nr, mold=sendbuf)
#else
          sendbuf(k+1:k+nr) = reshape (m(i)% ptr(:,:,kl:ku,:), [nr])
#endif
          k = k + nr
       end if
    end do
FTRACE_END  ("alltoallv_multi_w:fillbuf")

    !------------------------
    ! allocate receive buffer
    !------------------------
    if (any(dace% pe == dests)) then
       if (sum (counts) /= size (xg)) then
          write(0,*) dace% pe, ":", sum (counts), size (xg)
          call finish ("alltoallv_multi_w","sum (counts) /= size (xg)")
       end if
       allocate (recvbuf(sum (counts))) ! Allocate full buffer on receiver
    else
       counts(:) = 0
       displs(:) = 0
       allocate (recvbuf(0))            ! Allocate dummy buffer on receiver
    end if

    !--------------
    ! MPI alltoallv
    !--------------
    call p_alltoall (sendbuf, recvbuf, sendcounts=scounts, recvcounts=counts)

    !-------------------------------
    ! fetch data from receive buffer
    !-------------------------------
    if (any(dace% pe == dests)) then
       do j = 1, dc% nproc2
          do i = 1, dc% nproc1
             pe = dc% pe_12 (i,j)
             xg (dc%ilim1(i-1):dc%ilim1(i)-1,                       &
                 dc%ilim2(j-1):dc%ilim2(j)-1,                       &
                              :             ,                       &
                              :             ) =                     &
               reshape (recvbuf(1+displs(pe):counts(pe)+displs(pe)),&
                        shape = (/dc%ilim1(i)-dc%ilim1(i-1),        &
                                  dc%ilim2(j)-dc%ilim2(j-1),        &
                                  size(xg,3),size(xg,4)    /)       )
          end do
       end do
    end if

    deallocate (counts, displs, scounts)

  end subroutine alltoallv_multi_w
!==============================================================================
#ifdef  USE_MPI_GATHER
!------------------------------------------------------------------------------
  subroutine gather_level (xg, x, dc, dest)
  real(wp)        ,pointer     :: xg (:,:,:) ! Global fields
  real(wp)        ,intent(in)  ::  x (:,:,:) ! Local  fields
  type(t_atm_dec) ,intent(in)  :: dc         ! Domain decomposition information
  integer         ,intent(in)  :: dest       ! Processor holding global fields

    integer               :: i, j       ! loop indices
    integer               :: pe         ! sending pe
    integer               :: nd         ! # of diamonds
    integer,  pointer     :: counts(:)  ! recvcounts
    integer,  pointer     :: displs(:)  ! displacements
    real(wp), allocatable :: sendbuf(:) ! temporary send buffer
    real(wp), allocatable :: recvbuf(:) ! temporary recv buffer

    call setup_sendrecv (dc, counts, displs)
    nd = size (x, dim=3)
    counts(:) = counts(:) * nd
    displs(:) = displs(:) * nd

    if (counts(dace% pe) /= size (x)) then
       write(0,*) dace% pe, ":", counts(dace% pe), size (x)
       call finish ("gather_level","counts(pe) /= size (x)")
    end if
    if (dace% pe == dest) then
       if (sum (counts) /= size (xg)) then
          write(0,*) dace% pe, ":", sum (counts), size (xg)
          call finish ("gather_level","sum (counts) /= size (xg)")
       end if
    end if

FTRACE_BEGIN("gather_level:gatherv")

    if (dace% pe == dest) then
       allocate (recvbuf(sum (counts))) ! Allocate full buffer on receiver
    else
       allocate (recvbuf(0))            ! Allocate dummy buffers on sender
    end if
    !-----------------
    ! Fill send buffer
    !-----------------
    allocate (sendbuf(counts(dace% pe)))
    sendbuf(:) = reshape (x, shape (sendbuf))
    call p_gatherv (sendbuf, recvbuf, root=dest, recvcounts=counts)
    deallocate (sendbuf)
    if (dace% pe == dest) then
       !--------------------------------
       ! Gather data from receive buffer
       !--------------------------------
       do j = 1, dc% nproc2
          do i = 1, dc% nproc1
             pe = dc% pe_12 (i,j)
             xg (dc%ilim1(i-1):dc%ilim1(i)-1,   &
                 dc%ilim2(j-1):dc%ilim2(j)-1,   &
                              :             ) = &
                  reshape (recvbuf(1+displs(pe):counts(pe)+displs(pe)), &
                           shape = (/ dc%ilim1(i)-dc%ilim1(i-1),        &
                                      dc%ilim2(j)-dc%ilim2(j-1), nd /)  )
          end do
       end do
    end if
    deallocate (recvbuf, counts, displs)

FTRACE_END  ("gather_level:gatherv")

  end subroutine gather_level
!------------------------------------------------------------------------------
  subroutine gather_multi (xg, x, dc, dest)
  real(wp)        _POINTER     :: xg (:,:,:,:) ! Global fields
  real(wp)        ,intent(in)  ::  x (:,:,:,:) ! Local  fields
  type(t_atm_dec) ,intent(in)  :: dc           ! Domain decomposition inform.
  integer         ,intent(in)  :: dest         ! Processor holding global flds.

    integer               :: i, j       ! loop indices
    integer               :: pe         ! sending pe
    integer               :: nz         ! # of levels
    integer               :: nd         ! # of diamonds
    integer,  pointer     :: counts(:)  ! recvcounts
    integer,  pointer     :: displs(:)  ! displacements
    real(wp), allocatable :: sendbuf(:) ! temporary send buffer
    real(wp), allocatable :: recvbuf(:) ! temporary recv buffer

    call setup_sendrecv (dc, counts, displs)
    nz = size (x, dim=3)
    nd = size (x, dim=4)
    counts(:) = counts(:) * (nz * nd)
    displs(:) = displs(:) * (nz * nd)

    if (counts(dace% pe) /= size (x)) then
       write(0,*) dace% pe, ":", counts(dace% pe), size (x)
       call finish ("gather_multi","counts(pe) /= size (x)")
    end if
    if (dace% pe == dest) then
       if (sum (counts) /= size (xg)) then
          write(0,*) dace% pe, ":", sum (counts), size (xg)
          call finish ("gather_multi","sum (counts) /= size (xg)")
       end if
    end if

FTRACE_BEGIN("gather_multi:gatherv")

    if (dace% pe == dest) then
       allocate (recvbuf(sum (counts))) ! Allocate full buffer on receiver
    else
       allocate (recvbuf(0))            ! Allocate dummy buffers on sender
    end if
    !-----------------
    ! Fill send buffer
    !-----------------
    allocate (sendbuf(counts(dace% pe)))
    sendbuf(:) = reshape (x, shape (sendbuf))
    call p_gatherv (sendbuf, recvbuf, root=dest, recvcounts=counts)
    deallocate (sendbuf)
    if (dace% pe == dest) then
       !--------------------------------
       ! Gather data from receive buffer
       !--------------------------------
       do j = 1, dc% nproc2
          do i = 1, dc% nproc1
             pe = dc% pe_12 (i,j)
             xg (dc%ilim1(i-1):dc%ilim1(i)-1,   &
                 dc%ilim2(j-1):dc%ilim2(j)-1,   &
                              :             ,   &
                              :             ) = &
                  reshape (recvbuf(1+displs(pe):counts(pe)+displs(pe)),     &
                           shape = (/ dc%ilim1(i)-dc%ilim1(i-1),            &
                                      dc%ilim2(j)-dc%ilim2(j-1), nz, nd /)  )
          end do
       end do
    end if
    deallocate (recvbuf, counts, displs)

FTRACE_END  ("gather_multi:gatherv")

  end subroutine gather_multi
!------------------------------------------------------------------------------
#else
!------------------------------------------------------------------------------
  subroutine gather_level (xg, x, dc, dest)
  real(wp)        ,pointer     :: xg (:,:,:) ! Global fields
  real(wp)        ,intent(in)  ::  x (:,:,:) ! Local  fields
  type(t_atm_dec) ,intent(in)  :: dc         ! Domain decomposition information
  integer         ,intent(in)  :: dest       ! Processor holding global fields

    integer :: i,j,pe

    do i = 1, dc% nproc1
      do j = 1, dc% nproc2
        pe = dc% pe_12 (i,j)
        if (pe==dest) then
          if (dace% pe == dest) &
            xg (dc%ilim1(i-1):dc%ilim1(i)-1, &
                dc%ilim2(j-1):dc%ilim2(j)-1, &
                             :             ) = x
        else
          if (dace% pe == dest) then
            call p_recv (xg (dc%ilim1(i-1):dc%ilim1(i)-1, &
                             dc%ilim2(j-1):dc%ilim2(j)-1, &
                                          :             ),&
                         pe, p_tag=1)
          else if (pe == dace% pe) then
            call p_send (x, dest, p_tag=1)
          endif
        endif
      end do
    end do

  end subroutine gather_level
!------------------------------------------------------------------------------
#endif
!==============================================================================
  subroutine setup_sendrecv (dc, counts, displs, count2, pes)
    type(t_atm_dec)   ,intent(in) :: dc        ! Domain decomposition info
    integer           ,pointer    :: counts(:) ! (out) sendcounts/recvcounts
    integer           ,pointer    :: displs(:) ! (out) displacements
    integer ,optional ,pointer    :: count2(:) ! (out) alltoall recv/send cnts
    integer ,optional ,intent(in) :: pes   (:) ! I/O processors
    !-----------------------------------------------
    ! Set up sendcounts/recvcounts and displacements
    ! for the given horizontal domain decomposition
    !-----------------------------------------------
    integer :: i, j, pe
!   logical :: warned = .false.         ! warn only once
    logical :: warned = .true.          ! do not warn

    allocate (counts(0:dace% npe-1))
    allocate (displs(0:dace% npe  ))
    counts(:) = 0
    !-----------------------------
    ! Derive sendcounts/recvcounts
    !-----------------------------
    do j = 1, dc% nproc2
      do i = 1, dc% nproc1
        pe = dc% pe_12(i,j)
        counts(pe) = (dc% ilim1(i) - dc% ilim1(i-1)) &
                   * (dc% ilim2(j) - dc% ilim2(j-1))
      end do
    end do
    !---------------------
    ! Derive displacements
    !---------------------
    displs(0) = 0
    do pe = 1, dace% npe
      displs(pe) = displs(pe-1) + counts(pe-1)
    end do

    if (any (counts == 0) .and. dace% lpio .and. .not. warned) then
      write (0,*)
      write (0,*) "setup_sendrecv: WARNING: some sendcounts/recvcounts are 0!"
      write (0,*) counts(:)
      write (0,*) "setup_sendrecv: INFO: further warnings suppressed"
      write (0,*)
      warned = .true.
    end if

    !-------------------------------------------
    ! derive counts on remote side for alltoallv
    !-------------------------------------------
    if (present (count2)) then
      allocate  (count2(0:dace% npe-1))
      count2 = 0
      do i = 1, size(pes)
         if (pes(i) >= 0) count2 (pes(i)) = counts (dace% pe)
      end do
    endif

  end subroutine setup_sendrecv
!==============================================================================
end module mo_atm_transp
