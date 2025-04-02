!
!+ wavelet liftung scheme implementation
!
! $Id$
!
MODULE mo_lift_1d
!
! Description: Simple wavelet liftung scheme implementation.  The
!   number of coefficients in wavelet representation is 2x the number
!   of coefficients in physical space (frame).
!   Approch (analysis):
!   - Interpolate to coarser grid
!   - Re-interpolate fo fine grid
!   - optionally smooth the re-interpolated field
!   - keep the difference between orininal and re-interpolated (smoothed) field
!   - iterate the procedure on coarse grid
!   Approach (synthesis):
!   - Interpolate fo finer grid
!   - optionally smooth the interpolated field
!   - add coefficients at finer grid and interpolated (smoothed) field
!   - iterate the procedure on fine grid
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_42        2015-06-08 Andreas Rhodin
!  Simple wavelet liftung scheme implementation
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

  !-------------
  ! Modules used
  !-------------
  use mo_kind,       only: wp        ! working precision kind parameter
  use mo_exception,  only: finish    ! abort in case of error
  use mo_1dmra,      only: TR_SYN,  &! Synthesis  (w -> x)
                           TR_ANA,  &! Analysis   (x -> w)
!                          TR_DUAL, &! Dual       (adjoint analysis)
                           TR_ADJ    ! Adjoint    (adjoint synthesis)
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: lift_1d  ! 1-d lifting scheme implementation
  public :: smooth_d ! factor for smoothing in analysis step
  public :: smooth_u ! flag to smooth in synthesis step
  public :: cyclic   ! cyclic coordinate

  !-----------------
  ! module variables
  !-----------------
  real(wp) :: smooth_d = 1._wp
  logical  :: smooth_u = .false.
  logical  :: cyclic   = .false.

!==============================================================================
contains
!==============================================================================

  subroutine lift_1d (x, w, trans)
  real(wp) ,intent(inout)  :: x(0:)   ! array to transform
  real(wp) ,intent(inout)  :: w(0:)   ! wavelet transform
  integer  ,intent(in)     :: trans   ! transformation
  !------------------------------------------------------------------
  ! 1-d lifting scheme implementation.  The number of coefficients in
  ! wavelet representation is 2x the number of coefficients in
  ! physical space (frame). In dependence on the value of 'trans'
  ! either the analysis, synthesis or adjoint synthesis step is
  ! performed.
  !
  ! w(0:n-2) holds the wavelet coefficients
  ! x(0:n-1) holds the original field
  ! w(0:n-1) is not used
  ! n = size(w); m=n/2
  ! m is a power of 2
  !------------------------------------------------------------------

    integer :: m
    m = size(x)

    select case (trans)
    case (TR_ANA)
      !---------
      ! analysis
      !---------
      w (0:m-1) = x
      w (m:   ) = 0._wp
      call up (w)
    case (TR_SYN)
      !----------
      ! synthesis
      !----------
      call down      (x, w)
    case (TR_ADJ)
      !--------
      ! adjoint
      !--------
      call adj  (x, w)
    case default
      write (6,*)  'lift_1d:  invalid transform =',trans
      call finish ('lift_1d','invalid transform')
    end select
  end subroutine lift_1d
!------------------------------------------------------------------------------
  recursive subroutine up (w)
  real(wp) ,intent(inout) :: w(0:)
  !-----------------------------------------------
  ! analysis step of the 1d wavelet lifting scheme
  ! on exit : w(0:n-2) holds the wavelet coefficients
  ! on entry: w(0:m-1) holds the original field
  !           w(0:n-1) is not used
  !           n = size(w); m=n/2
  !           m is a power of 2
  !-----------------------------------------------
    integer  :: n                  ! size of w
    integer  :: m                  ! site of x
    integer  :: l                  ! size of next upper level
    integer  :: i,j,k
    real(wp) :: v (0:size(w)/4-1)  ! average
    real(wp) :: y (0:size(w)/2-1)  ! re-interpolated field
    real(wp) :: z (0:size(w)/2-1)  ! (optionally) smoothed field

    n = size (w)
    m = n / 2
    l = m / 2

    !--------------------
    ! return on top level
    !--------------------
    if (m == 0) call finish('lift_1d:up','m == 0')
    if (m == 1 .or. m == 3) then
      w(m:) = 0._wp
      return
    endif

    !--------
    ! average
    !--------
    do i = 0, l-1
      k = i * 2
      v (i) = (w(k) + w(k+1)) / 2._wp
    end do

    !-------
    ! smooth
    !-------
    if (l>1 .and. smooth_u) then
      do i = 1, l-2
        j = i + m
        w (j)  = 0.25_wp * v(i-1) + 0.5_wp * v(i)   + 0.25 * v(i+1)
      end do
      if (cyclic) then
        w(m)     = 0.25_wp * v(l-1) + 0.5_wp * v( 0 ) + 0.25 * v(1)
        w(m+l-1) = 0.25_wp * v(l-2) + 0.5_wp * v(l-1) + 0.25 * v(0)
      else
        w(m)     = v( 0 )
        w(m+l-1) = v(l-1)
      endif
    else
      w(m:m+l-1) = v
    endif

    !------------
    ! interpolate
    !------------
    do i = 0, l-2
      j = i + m
      k = i * 2
      y(k+1) = 0.75_wp * w(j) + 0.25_wp * w(j+1)
      y(k+2) = 0.25_wp * w(j) + 0.75_wp * w(j+1)
    end do
    if (cyclic) then
      y( 0 ) = 0.75_wp * w(m) + 0.25_wp * w(m+l-1)
      y(m-1) = 0.25_wp * w(m) + 0.75_wp * w(m+l-1)
    else
      if (l > 1) then
        y( 0 ) = 1.25_wp * w(m)     - 0.25 * w(m+1)
        y(m-1) = 1.25_wp * w(m+l-1) - 0.25 * w(m+l-2)
      else
        y( 0 ) = w(m)
        y(m-1) = w(m+l-1)
      endif
    endif

    !-------
    ! smooth
    !-------
    if (smooth_d > 0._wp) then
      do k = 1, m-2
        z(k)   = 0.25_wp * y(k-1) - 0.5_wp * y(k)   + 0.25_wp * y(k+1)
      end do
      if (cyclic) then
        z(0)   = 0.25_wp * y(m-1) - 0.5_wp * y( 0 ) + 0.25_wp * y(1)
        z(m-1) = 0.25_wp * y(m-2) - 0.5_wp * y(m-1) + 0.25_wp * y(0)
      else
        z(0)   = 0._wp
        z(m-1) = 0._wp
      endif
      z = y + smooth_d * z
    else
      z = y
    endif

    !------------------------------------------------
    ! keep deviation: original - re-constructed field
    !------------------------------------------------
    do k = 0, m-1
      w(k) = w(k) - z(k)
    end do
    !--------
    ! iterate
    !--------
    call up (w(m:))

  end subroutine up

!------------------------------------------------------------------------------
  recursive subroutine down (x, w)
  real(wp) ,intent(out) :: x(0:)  ! physical representation
  real(wp) ,intent(in ) :: w(0:)  ! wavelet coefficients
  !--------------------------------------------------
  ! synthesis step of the 1d wavelet lifting scheme
  ! w(0:n-2) holds the wavelet coefficients
  ! x(0:m-1) holds the field in physical representation
  ! w(0:n-1) is not used
  !           n = size(w); m=n/2
  !           m is a power of 2
  !--------------------------------------------------

    integer  :: n
    integer  :: m
    integer  :: l
    integer  :: i,j,k

    real(wp) :: y (0:size(x)-1)

    real(wp) :: z (0:size(x)/2-1)

    n = size (w)
    m = n / 2
    l = m / 2

    !--------------------
    ! return on top level
    !--------------------
    if (m == 0) call finish('lift_1d:up','m == 0')
    if (m == 1 .or. m == 3) then
      x = w(:m-1)
      return
    endif

    !--------
    ! iterate
    !--------
    call down (z, w(m:))

    !------------
    ! interpolate
    !------------
    do i = 0, l-2
      k = i * 2
      y(k+1) = 0.75_wp * z(i) + 0.25_wp * z(i+1)
      y(k+2) = 0.25_wp * z(i) + 0.75_wp * z(i+1)
    end do
    if (cyclic) then
      y( 0 ) = 0.75_wp * z(0) + 0.25_wp * z(l-1)
      y(m-1) = 0.25_wp * z(0) + 0.75_wp * z(l-1)
    else
      if (l > 1) then
        y( 0 ) = 1.25_wp * z(0)   - 0.25_wp * z(1)
        y(m-1) = 1.25_wp * z(l-1) - 0.25_wp * z(l-2)
      else
        y ( 0 ) = z (0)
        y (m-1) = z (l-1)
      endif
    endif

    !-------
    ! smooth
    !-------
    if (smooth_d > 0._wp) then
      do k = 1, m-2
        x(k)   = 0.25_wp * y(k-1) - 0.5_wp * y(k)   + 0.25_wp * y(k+1)
      end do
      if (cyclic) then
        x(0)   = 0.25_wp * y(m-1) - 0.5_wp * y( 0 ) + 0.25_wp * y(1)
        x(m-1) = 0.25_wp * y(m-2) - 0.5_wp * y(m-1) + 0.25_wp * y(0)
      else
        x(0)   = 0._wp
        x(m-1) = 0._wp
      endif
      x = y + smooth_d * x
    else
      x = y
    endif

    !--------------
    ! add deviation
    !--------------
    do k = 0, m-1
      x(k) = x(k) + w(k)
    end do

  end subroutine down

!------------------------------------------------------------------------------
  recursive subroutine adj (x, w)
  real(wp) ,intent(in ) :: x(0:)
  real(wp) ,intent(out) :: w(0:)
  !--------------------------------------------------------
  ! adjoint synthesis step of the 1d wavelet lifting scheme
  ! w(0:n-2) holds the wavelet coefficients
  ! x(0:m-1) holds the field in physical representation
  ! w(0:n-1) is not used
  !           n = size(w); m=n/2
  !           m is a power of 2
  !--------------------------------------------------------

    integer  :: n
    integer  :: m
    integer  :: l
    integer  :: i,j,k

    real(wp) :: y (0:size(x)-1)

    real(wp) :: z (0:size(x)/2-1)

    n = size (w)
    m = n / 2
    l = m / 2

    !----------
    ! top level
    !----------
    if (m == 0) call finish('lift_1d:up','m == 0')
    if (m == 1 .or. m == 3) then
      w(:m-1) = x
      w(m:)   = 0._wp
      return
    endif

    !--------------
    ! add deviation
    !--------------
    w(0:m-1) = x(0:m-1)

    !------------------
    ! adjoint smoothing
    !------------------
    if (smooth_d > 0._wp) then
      y = 0._wp
      do k = 1, m-2
        y(k)   = y(k)   - 0.5_wp  * x(k)
        y(k-1) = y(k-1) + 0.25_wp * x(k)
        y(k+1) = y(k+1) + 0.25_wp * x(k)
      end do
      if (cyclic) then
        y( 1 ) = y( 1 ) + 0.25_wp * x(0)
        y( 0 ) = y( 0 ) - 0.5_wp  * x(0) + 0.25_wp * x(m-1)
        y(m-1) = y(m-1) + 0.25_wp * x(0) - 0.5_wp  * x(m-1)
        y(m-2) = y(m-2) + 0.25_wp * x(m-1)
      endif
      y = x + smooth_d * y
    else
      y = x
    endif

    !----------------------
    ! adjoint interpolation
    !----------------------
    z = 0._wp
    do i = 0, l-2
      k = i * 2
      z( i ) = z( i ) + 0.75_wp * y(k+1) + 0.25_wp * y(k+2)
      z(i+1) = z(i+1) + 0.25_wp * y(k+1) + 0.75_wp * y(k+2)
    end do
    if (cyclic) then
      z( 0 ) = z( 0 ) + 0.75_wp * y( 0 ) + 0.25_wp * y(m-1)
      z(l-1) = z(l-1) + 0.25_wp * y( 0 ) + 0.75_wp * y(m-1)
    else
      if (l > 1) then
        z( 1 ) = z( 1 ) - 0.25_wp * y( 0 )
        z( 0 ) = z( 0 ) + 1.25_wp * y( 0 )
        z(l-1) = z(l-1) + 1.25_wp * y(m-1)
        z(l-2) = z(l-2) - 0.25_wp * y(m-1)
      else
        z ( 0 ) = z ( 0 ) + y ( 0 )
        z (l-1) = z (l-1) + y (m-1)
      endif
    endif

    !--------
    ! iterate
    !--------
    call adj (z, w(m:))

  end subroutine adj

!==============================================================================

!!$  subroutine loc_lift (w, l, r)
!!$  !----------------------------------------------
!!$  ! localise matrix in wavelet representation
!!$  !   l = 0 : no localisation
!!$  !       1 : diagonal localisation
!!$  !       2 : zero off-diagonal blocks
!!$  !   l = 3 : localisation on off diagonal blocks
!!$  !----------------------------------------------
!!$  real(wp) ,intent(inout) :: w(0:,0:)
!!$  integer  ,intent(in)    :: l
!!$  real(wp) ,intent(in)    :: r
!!$
!!$    integer  :: n,m1,l1,m2,l2
!!$    real(wp) :: rl
!!$
!!$    if (l==0) return
!!$
!!$    rl = sqrt(10._wp/3._wp) * r
!!$
!!$    n = size(w,1)
!!$    l1 = 0
!!$    m1 = n / 2
!!$    do
!!$      l2 = 0
!!$      m2 = n / 2
!!$      do
!!$        !--------------
!!$        !localise block
!!$        !--------------
!!$        call loc_block (w(l1:l1+m1-1,l2:l2+m2-1), l, rl)
!!$
!!$        if (m2 == 1 .or. m2 == 3) exit
!!$        l2 = l2 + m2
!!$        m2 = m2 / 2
!!$      end do
!!$
!!$      if (m1 == 1 .or. m1 == 3) exit
!!$      l1 = l1 + m1
!!$      m1 = m1 / 2
!!$    end do
!!$
!!$  end subroutine loc_lift
!!$!---------------------------------------------------------------------------
!!$
!!$  subroutine loc_block (w, l, rl)
!!$  real(wp) ,intent(inout) :: w(0:,0:)
!!$  integer  ,intent(in)    :: l
!!$  real(wp) ,intent(in)    :: rl
!!$
!!$    integer  :: n1,n2,i,j,k,m
!!$    real(wp) :: r
!!$
!!$    n1 = size (w,1)
!!$    n2 = size (w,2)
!!$
!!$    select case (l)
!!$    !---------
!!$    ! diagonal
!!$    !---------
!!$    case (1)
!!$      if (n1==n2) then
!!$        do i = 0,n1-1
!!$        do j = 0,n2-1
!!$          if (i/=j) w(i,j) = 0._wp
!!$        end do
!!$        end do
!!$      else
!!$        w = 0._wp
!!$      endif
!!$    !---------
!!$    ! localize
!!$    !---------
!!$    case (2,3)
!!$      !---------------
!!$      ! diagonal block
!!$      !---------------
!!$      if (n1==n2) then
!!$        do i = 0,n1-1
!!$        do j = 0,n2-1
!!$          w(i,j) = w(i,j) * gaspari_cohn (real(abs(i-j),wp), rl)
!!$        end do
!!$        end do
!!$      !------------------
!!$      ! offdiagonal block
!!$      !------------------
!!$      else
!!$        select case (l)
!!$        !-------------------------
!!$        ! not on offdiagonal block
!!$        !-------------------------
!!$        case (2)
!!$          w = 0._wp
!!$        !---------------------
!!$        ! on offdiagonal block
!!$        !---------------------
!!$        case (3)
!!$          if (n1 > n2) then
!!$            m = n1 / n2
!!$            r = rl / m
!!$            do i = 0,n1-1
!!$              k = i / m
!!$              do j = 0,n2-1
!!$                w(i,j) = w(i,j) * gaspari_cohn (real(abs(k-j),wp), r)
!!$              end do
!!$            end do
!!$          else
!!$            m = n2 / n1
!!$            r = rl / m
!!$            do j = 0,n2-1
!!$              k = j / m
!!$              do i = 0,n1-1
!!$                w(i,j) = w(i,j) * gaspari_cohn (real(abs(i-k),wp), r)
!!$              end do
!!$            end do
!!$          endif
!!$        end select
!!$      endif
!!$    end select
!!$
!!$  end subroutine loc_block

!==============================================================================
end module mo_lift_1d
