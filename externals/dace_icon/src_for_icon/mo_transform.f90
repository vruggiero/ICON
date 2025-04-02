!
!+ FFT and Legendre-transforms and adjoint operations
!
! $Id$
!
MODULE mo_transform
!
! Description:
!   Module 'mo_transform' holds generic routines to perform FFT or
!   Legendre-transformations on 2D or 3D fields and the adjoint
!   operations.
!
!   Normalization may be chosen as follows:
!     0: as defined in ECHAM (default)
!     1: strictly orthonormal transformation
!
!   Adjoints:
!
!     The adjoint routines overwrite the intent(out) argument !
!
!     Note that for the Legendre-transformation the scalar product
!     in gridpoint space is:
!
!       (A,B) = Sum_i=1,m,j=1,n ( A(i,j) * B(i,j) * gw(j) / m)
!
!     For the Fourier-transformation the scalar product is simply:
!
!       (A,B) = Sum_i=1,m,j=1,n ( A(i,j) * B(i,j) )
!
!     In Fourier and spectral space the scalar product is simply:
!
!       (A,B) = Sum_coefficients ( A(c) * B(c) )
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
! V1_4         2009/03/26 Harald Anlauf
!  fftc: further optimizations for NEC-SX
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Author:
! Andreas Rhodin  MPIfM  1999
!----------------------------------------------------------------
  !-------------
  ! modules used
  !-------------
  use mo_fft,      only: trig,     &! array used by 'fft991'.
                         nfax,     &! array used by 'fft991'.
                         nfft,     &! maximum length of calls to 'fft991'.
                         inifft,   &! initialization routine
                         fft991cy   ! multiple fast real periodic transform
  use mo_legendre, only: inileg   ,&! set up polynominals for Legendre transf.
                         gaw=>gauaw ! Calculate Gaussian abscissas and weights
  use mo_kind,     only: wp         ! working precision
  implicit none
  !----------------
  ! public entities
  !----------------
  private
  public :: fft             ! generic fft routine
  public :: fft_ad          ! adjoint fft routine
  public :: lgtd            ! direct Legendre transform
  public :: lgti            ! inverse Legendre transform
  public :: lgtd_ad         ! adjoint  Legendre transform
  public :: lgti_ad         ! adjoint inverse Legendre transform
  public :: nsp_nn          ! derive number of spectral coeff. from truncation
  public :: gauaw           ! Calculate Gaussian abscissas and weights
  public :: spread_m0_n     ! spread m=0 coefficients over n
  public :: spread_m0_l     ! spread m=0 coefficients over l
  public :: NORM_ECHAM      ! Flag: ECHAM4 normalization for transforms
  public :: NORM_ORTHO      ! Flag: Strictly orthonormal transforms
  public :: NORM_ECHAM_AD   ! Flag: ECHAM4 normalization for adjoint transforms
  public :: cnorm2
!==============================================================================
  !-----------
  ! interfaces
  !-----------
  interface fft
    module procedure fftc       ! real <-> complex
    module procedure fftc1      ! real <-> complex (1d)
    module procedure fftc3      ! real <-> complex
    module procedure fft2
    module procedure fft3
  end interface fft
!------------------------------------------------------------------------------
  interface fft_ad
    module procedure fftc_adj   ! real <-> complex
    module procedure fft2_adj
    module procedure fft3_adj
  end interface fft_ad
!------------------------------------------------------------------------------
  interface lgtd
    module procedure lgt1
    module procedure lgt2
    module procedure lgt3
  end interface lgtd
!------------------------------------------------------------------------------
  interface lgti
    module procedure lgti1
    module procedure lgti2
    module procedure lgti3
  end interface lgti
!------------------------------------------------------------------------------
  interface lgtd_ad
    module procedure lgt2_adj
    module procedure lgt3_adj
  end interface lgtd_ad
!------------------------------------------------------------------------------
  interface lgti_ad
    module procedure lgti2_adj
    module procedure lgti3_adj
  end interface lgti_ad
!------------------------------------------------------------------------------
  integer, parameter :: NORM_ECHAM    = 0     ! Normalization flags for
  integer, parameter :: NORM_ORTHO    = 1     ! transformation routines
  integer, parameter :: NORM_ECHAM_AD = 2     ! (do not change!)
  real(wp)           :: cnorm2        = 1._wp
!==============================================================================
contains
!==============================================================================
subroutine fft2 (x, isign, l1, l2, norm)
real(wp) ,intent(inout)           :: x(:,:) ! array to transform
integer  ,intent(in)              :: isign  ! +1: spectral -> gridpoint, -1: <-
logical  ,intent(in)    ,optional :: l1     ! (def=.true.)  transform 1st index
logical  ,intent(in)    ,optional :: l2     ! (def=.false.) transform 2nd index
integer  ,intent(in)    ,optional :: norm   ! normalization flag
!------------------------------------------------------------------------------
! Description:
!
! Calls fortran-versions of fft s.
!
! Method:
!
! Subroutine 'fft991cy' - multiple fast real periodic transform
!
! Real transform of length n performed by removing redundant
! operations from complex transform of length n.
!
! x       is the array containing input & output data.
! isign = +1 for transform from spectral to gridpoint
!       = -1 for transform from gridpoint to spectral.
! l1      set to .true. to perform fft along 1st index of x (default=.true.)
! l2      set to .true. to perform fft along 2nd index of x (default=.false.)
! norm    0(Default): ECHAM4 norm, 1: strictly orthonormal transform
!
! ordering of coefficients:
! a(0),a(1),b(1),a(2),b(2),.,a(n/2)   (where b(0)=b(n/2)=0).
!
! ordering of data:
! x(0),x(1),x(2),.,x(n-1).
!
! Vectorization is achieved on cray by doing the transforms
! in parallel.
!
! n must be composed of factors 2,3 & 5 but does not have to be even.
!
! For norm==1 This fft routine is strictly orthonormal. The definition of
! coefficients differs from that of fft991cy by a factor of SQRT(N) for all
! coefficiens and an additionlal factor of SQRT(2) for a(0) and a(n/2) .
! For norm==0 normalization is the same as in fft991cy .
!
! definition of transforms in fft991cy:
!
!  isign=+1: x(j)=sum(k=0,.,n-1)(c(k)*exp(2*i*j*k*pi/n))
!      where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
!
!  isign=-1: a(k)= (1/n)*sum(j=0,.,n-1)(x(j)*cos(2*j*k*pi/n))
!            b(k)=-(1/n)*sum(j=0,.,n-1)(x(j)*sin(2*j*k*pi/n))
!
! definition of transforms here (with norm=1):
!
!  isign=+1: x(j)=(1/SQRT(N))*sum(k=0,.,n-1)(c(k)*exp(2*i*j*k*pi/n))
!      where c(k)=SQRT(2)*(a(k)+i*b(k)) and c(n-k)=SQRT(2)*(a(k)-i*b(k))
!      but   c(0)=a(0) and c(n/2)=a(n/2)
!
!  isign=-1: a(k)=  (1/SQRT(n))*SQRT(2)*sum(j=0,.,n-1)(x(j)*cos(2*j*k*pi/n))
!            b(k)= -(1/SQRT(n))*SQRT(2)*sum(j=0,.,n-1)(x(j)*sin(2*j*k*pi/n))
!      but   a(0)=  (1/SQRT(n))        *sum(j=0,.,n-1)(x(j))
!            a(n/2)=(1/SQRT(n))        *sum(j=0,.,n-1)(x(j)*cos(  j  *pi))
!------------------------------------------------------------------------------
  !-----------------------------------
  ! local variables passed to fft991cy
  !-----------------------------------
  real(wp) ,allocatable :: a    (:) ! array containing input & output data.
  real(wp) ,allocatable :: work (:) ! area of size (n+1)*min(lot,nfft).
  integer               :: inc      ! increment within each data 'vector'
  integer               :: jump     ! incr. betw. start of each data vector.
  integer               :: n        ! length of the data vectors.
  integer               :: lot      ! is the number of data vectors.
  integer               :: i,i0
  real(wp)              :: rnorm1, rnorm2
  !--------------------
  ! optional parameters
  !--------------------
  logical :: ll1, ll2 ! transform 1st, 2nd index
  integer :: lno      ! transformation flag
  ll1 = .true.     ;if (present(l1))   ll1 = l1
  ll2 = .false.    ;if (present(l2))   ll2 = l2
  lno = NORM_ECHAM ;if (present(norm)) lno = norm
  !---------------
  ! transformation
  !---------------
  if(ll1) then
    !------------------
    ! fft along index 1
    !------------------
    n    = size(x,1)
    lot  = size(x,2)
    inc  = 1
    jump = n+2
    allocate (a ((n+2)*lot))
    allocate (work ((n+1)*min(lot,nfft)))
    if (nfax(10)/=n) call inifft(n)
    rnorm1 = sqrt(real(n,wp)) ** lno
    rnorm2 = sqrt(2._wp)      ** lno
    select case (isign)
    case (-1)
      !----------------------
      ! gridpoint to spectral
      !----------------------
!NEC$ ivdep
      do i=1,lot
        i0 = (i-1)*(n+2)
        a(i0+1:i0+n)     = x(:,i) * rnorm1
        a(i0+n+1:i0+n+2) = 0._wp
      end do
      call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
!NEC$ ivdep
      do i=1,lot
        i0 = (i-1)*(n+2)
        x(1    ,i) = a(i0+1)
        x(2:n-1,i) = a(i0+3:i0+n) * rnorm2
        x(n    ,i) = a(i0+n+1)
      end do
    case (+1)
      !----------------------
      ! spectral to gridpoint
      !----------------------
!NEC$ ivdep
      do i=1,lot
        i0 = (i-1)*(n+2)
        a(i0+1)        = x(1    ,i)
        a(i0+2)        = 0._wp
        a(i0+3:i0+n)   = x(2:n-1,i) / rnorm2
        a(i0+n+1)      = x(n    ,i)
        a(i0+n+2)      = 0._wp
      end do
      call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
!NEC$ ivdep
      do i=1,lot
        i0 = (i-1)*(n+2)
        x(:,i) = a(i0+1:i0+n) / rnorm1
      end do
    end select
    deallocate (a)
    deallocate (work)
  endif
  if(ll2) then
    !------------------
    ! fft along index 2
    !------------------
    n    = size(x,2)
    lot  = size(x,1)
    inc  = 1
    jump = n+2
    allocate (a ((n+2)*lot))
    allocate (work ((n+1)*min(lot,nfft)))
    if (nfax(10)/=n) call inifft(n)
    rnorm1 = sqrt(real(n,wp)) ** lno
    rnorm2 = sqrt(2._wp)      ** lno
    select case (isign)
    case (-1)
      !----------------------
      ! gridpoint to spectral
      !----------------------
!NEC$ ivdep
      do i=1,lot
        i0 = (i-1)*(n+2)
        a(i0+1:i0+n)     = x(i,:) * rnorm1
        a(i0+n+1:i0+n+2) = 0._wp
      end do
      call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
      do i=1,lot
        i0 = (i-1)*(n+2)
        x(i,1    ) = a(i0+1)
        x(i,2:n-1) = a(i0+3:i0+n) * rnorm2
        x(i,n    ) = a(i0+n+1)
      end do
    case (+1)
      !----------------------
      ! spectral to gridpoint
      !----------------------
!NEC$ ivdep
      do i=1,lot
        i0 = (i-1)*(n+2)
        a(i0+1)      = x(i,1)
        a(i0+2)      = 0._wp
        a(i0+3:i0+n) = x(i,2:n-1) / rnorm2
        a(i0+n+1)    = x(i,n)
        a(i0+n+2)    = 0._wp
      end do
      call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
!NEC$ ivdep
      do i=1,lot
        i0 = (i-1)*(n+2)
        x(i,:) = a(i0+1:i0+n) / rnorm1
      end do
    end select
    deallocate (a)
    deallocate (work)
  endif
end subroutine fft2
!------------------------------------------------------------------------------
subroutine fftc3 (x, c, isign, norm)
real(wp)    ,intent(inout)        :: x(:,:,:) ! array to transform
complex(wp) ,intent(inout)        :: c(:,:,:) ! array to transform
integer     ,intent(in)           :: isign    ! +1: spectral->gridpoint, -1: <-
integer     ,intent(in) ,optional :: norm     ! normalization flag
  integer :: i
  do i=1,size(x,3)
    call fftc (x(:,:,i), c(:,:,i), isign, norm)
  end do
end subroutine fftc3
!------------------------------------------------------------------------------
subroutine fftc1 (x, c, isign, norm)
  real(wp)    ,intent(inout)        :: x(:)  ! vector to transform
  complex(wp) ,intent(inout)        :: c(:)  ! vector to transform
  integer     ,intent(in)           :: isign ! +1: spectral->gridpoint, -1: <-
  integer     ,intent(in) ,optional :: norm  ! normalization flag
  real(wp)    :: x2(size (x),1)
  complex(wp) :: c2(size (c),1)
  x2(:,1) = x(:)
  c2(:,1) = c(:)
  call fftc (x2(:,:), c2(:,:), isign, norm)
  x(:) = x2(:,1)
  c(:) = c2(:,1)
end subroutine fftc1
!------------------------------------------------------------------------------
subroutine fftc (x, c, isign, norm)
real(wp)    ,intent(inout)        :: x(:,:) ! array to transform
complex(wp) ,intent(inout)        :: c(:,:) ! array to transform
integer     ,intent(in)           :: isign  ! +1: spectral -> gridpoint, -1: <-
integer     ,intent(in) ,optional :: norm   ! normalization flag
!------------------------------------------------------------------------------
! Description:
!
! Same as fft2, but Fourier coefficients are stored in complex array C .
!
! Method:
!
! Subroutine 'fft991cy' - multiple fast real periodic transform
!
! Real transform of length n performed by removing redundant
! operations from complex transform of length n.
! perform fft along 1st index of x.
!
! x       is the array containing input or output data (gridpoint real)
! c       is the array containing input or output data (Fourier coefficients)
! isign = +1 for transform from spectral to gridpoint
!       = -1 for transform from gridpoint to spectral.
! norm    0(Default): ECHAM4 norm, 1: strictly orthonormal transform
!
! ordering of coefficients c=a+ib:
! a(0),b(0),a(1),b(1),a(2),b(2),.,a(n/2),b(n/2)
! where b(0)=b(n/2)=0; (n/2+1) complex locations required.
!
! ordering of data:
! x(0),x(1),x(2),.,x(n-1).
!
! Vectorization is achieved on cray by doing the transforms
! in parallel.
!
! n must be composed of factors 2,3 & 5 and must be even.
!------------------------------------------------------------------------------
  !-----------------------------------
  ! local variables passed to fft991cy
  !-----------------------------------
  real(wp) ,allocatable :: a    (:) ! array containing input & output data.
  real(wp) ,allocatable :: work (:) ! area of size (n+1)*min(lot,nfft).
  integer               :: inc      ! increment within each data 'vector'
  integer               :: jump     ! incr. betw. start of each data vector.
  integer               :: n        ! length of the data vectors.
  integer               :: lot      ! is the number of data vectors.
  integer               :: i,i0
  real(wp)              :: rnorm1
  real(wp)              :: rnorm2
  logical               :: even
  integer               :: n1
  !--------------------
  ! optional parameters
  !--------------------
  integer :: lno      ! transformation flag
  lno = NORM_ECHAM ;if (present(norm)) lno = norm
  !------------------
  ! fft along index 1
  !------------------
  n    = size(x,1)
  lot  = size(x,2)
  inc  = 1
  jump = n+2
  allocate (a ((n+2)*lot))
  allocate (work ((n+1)*min(lot,nfft)))
  if (nfax(10)/=n) call inifft(n)
  rnorm1 = sqrt(real(n,wp)) ** lno
  rnorm2 = sqrt(cnorm2)     ** lno
  even   = mod(n,2)==0
  n1     = n+1
  if(even) n1 = n
  select case (isign)
  case (-1)
    !----------------------
    ! gridpoint to spectral
    !----------------------
!NEC$ ivdep
    do i=1,lot
      i0 = (i-1)*(n+2)
      a(i0+1:i0+n)     = x(:,i) * rnorm1
      a(i0+n+1:i0+n+2) = 0._wp
    end do
    call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
!NEC$ ivdep
    do i=1,lot
      i0 = (i-1)*(n+2)
      c(1     ,i) = a(i0+1)
#if defined (__SX__)
      c(2:n1/2,i) = cmplx    (a(i0+3:i0+n1-1:2) * rnorm2, &
                              a(i0+4:i0+n1  :2) * rnorm2, kind=wp)
#else
      c(2:n1/2,i) = transfer (a(i0+3:i0+n1) * rnorm2, c(1:1,1))
#endif
      if (even) c(n/2+1,i) =  a(i0+n+1)
    end do
  case (+1)
    !----------------------
    ! spectral to gridpoint
    !----------------------
!NEC$ ivdep
    do i=1,lot
      i0 = (i-1)*(n+2)
      a(i0+1)             = real     (c(1       ,i) ,wp)
      a(i0+2)             = 0._wp
#if defined (__SX__)
      a(i0+3:i0+n1-1:2)   = real  (c(2:n1/2  ,i), kind=wp) / rnorm2
      a(i0+4:i0+n1  :2)   = aimag (c(2:n1/2  ,i)         ) / rnorm2
#else
      a(i0+3:i0+n1)       = transfer (c(2:n1/2  ,i) / rnorm2, a(1:1))
#endif
      if (even) a(i0+n+1) = real     (c(  n /2+1,i) ,wp)
      a(i0+n+2)           = 0._wp
    end do
    call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
!NEC$ ivdep
    do i=1,lot
      i0 = (i-1)*(n+2)
      x(:,i) = a(i0+1:i0+n) / rnorm1
    end do
  end select
  deallocate (a)
  deallocate (work)
end subroutine fftc
!==============================================================================
subroutine fft3 (x, isign, l1, l2, norm)
real(wp) ,intent(inout)           :: x(:,:,:)! array to transform
integer  ,intent(in)              :: isign   !+1: spectral -> gridpoint, -1: <-
logical  ,intent(in)    ,optional :: l1      !(def=.true.)  transform 1st index
logical  ,intent(in)    ,optional :: l2      !(def=.false.) transform 2nd index
integer  ,intent(in)    ,optional :: norm    ! normalization flag
!------------------------------------------------------------------------------
! Description:
!
! Calls fortran-versions of fft s.
!
! Method:
!
! Subroutine 'fft991cy' - multiple fast real periodic transform
!
! Real transform of length n performed by removing redundant
! operations from complex transform of length n.
!
! x       is the array containing input & output data.
! isign = +1 for transform from spectral to gridpoint
!       = -1 for transform from gridpoint to spectral.
! l1      set to .true. to perform fft along 1st index of x (default=.true.)
! l2      set to .true. to perform fft along 2nd index of x (default=.false.)
! norm    0(Default): ECHAM4 norm, 1: strictly orthonormal transform
!
! ordering of coefficients:
! a(0),a(1),b(1),a(2),b(2),.,a(n/2)   (where b(0)=b(n/2)=0).
!
! ordering of data:
! x(0),x(1),x(2),.,x(n-1).
!
! Vectorization is achieved on cray by doing the transforms
! in parallel.
!
! n must be composed of factors 2,3 & 5 but does not have to be even.
!
! This fft routine is strictly orthonormal. Thus the definition of coefficients
! differs from that of fft991cy by a factor of SQRT(N) for all coefficiens
! and an additionlal factor of SQRT(2) for a(0) and a(n/2)
!
! definition of transforms in fft991cy:
!
!  isign=+1: x(j)=sum(k=0,.,n-1)(c(k)*exp(2*i*j*k*pi/n))
!      where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
!
!  isign=-1: a(k)= (1/n)*sum(j=0,.,n-1)(x(j)*cos(2*j*k*pi/n))
!            b(k)=-(1/n)*sum(j=0,.,n-1)(x(j)*sin(2*j*k*pi/n))
!
! definition of transforms here (with norm=1):
!
!  isign=+1: x(j)=(1/SQRT(N))*sum(k=0,.,n-1)(c(k)*exp(2*i*j*k*pi/n))
!      where c(k)=SQRT(2)*(a(k)+i*b(k)) and c(n-k)=SQRT(2)*(a(k)-i*b(k))
!      but   c(0)=a(0) and c(n/2)=a(n/2)
!
!  isign=-1: a(k)=  (1/SQRT(n))*SQRT(2)*sum(j=0,.,n-1)(x(j)*cos(2*j*k*pi/n))
!            b(k)= -(1/SQRT(n))*SQRT(2)*sum(j=0,.,n-1)(x(j)*sin(2*j*k*pi/n))
!      but   a(0)=  (1/SQRT(n))        *sum(j=0,.,n-1)(x(j))
!            a(n/2)=(1/SQRT(n))        *sum(j=0,.,n-1)(x(j)*cos(  j  *pi))
!------------------------------------------------------------------------------
  !-----------------------------------
  ! local variables passed to fft991cy
  !-----------------------------------
  real(wp) ,allocatable :: a    (:) ! array containing input & output data.
  real(wp) ,allocatable :: work (:) ! area of size (n+1)*min(lot,nfft).
  integer               :: inc      ! increment within each data 'vector'
  integer               :: jump     ! incr. betw. start of each data vector.
  integer               :: n        ! length of the data vectors.
  integer               :: lot      ! is the number of data vectors.
  integer               :: n2, n3
  integer               :: i,i2,i3,i0
  real(wp)              :: rnorm1, rnorm2
  !--------------------
  ! optional parameters
  !--------------------
  logical :: ll1, ll2 ! transform 1st, 2nd index
  integer :: lno      ! normalization flag
  ll1 = .true.     ;if (present(l1)) ll1 = l1
  ll2 = .false.    ;if (present(l2)) ll2 = l2
  lno = NORM_ECHAM ;if (present(norm)) lno = norm
  !---------------
  ! transformation
  !---------------
  if(ll1) then
    !------------------
    ! fft along index 1
    !------------------
    n    = size(x,1)
    n2   = size(x,2)
    n3   = size(x,3)
    lot  = n2*n3
    inc  = 1
    jump = n+2
    allocate (a ((n+2)*lot))
    allocate (work ((n+1)*min(lot,nfft)))
    if (nfax(10)/=n) call inifft(n)
    rnorm1 = sqrt(real(n,wp)) ** lno
    rnorm2 = sqrt(2._wp)      ** lno
    select case (isign)
    case (-1)
      !----------------------
      ! gridpoint to spectral
      !----------------------
      i = 0
      do i3=1,n3
#if defined (__SX__)
!NEC$ ivdep
#endif
        do i2=1,n2
#if defined (__SX__)
          i = (i3-1)*n2 + i2
#else
          i=i+1
#endif
          i0 = (i-1)*(n+2)
          a(i0+1:i0+n)     = x(:,i2,i3) * rnorm1
          a(i0+n+1:i0+n+2) = 0._wp
        end do
      end do
      call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
      i = 0
      do i3=1,n3
#if defined (__SX__)
!NEC$ ivdep
#endif
        do i2=1,n2
#if defined (__SX__)
          i = (i3-1)*n2 + i2
#else
          i=i+1
#endif
          i0 = (i-1)*(n+2)
          x(1    ,i2,i3) = a(i0+1)
          x(2:n-1,i2,i3) = a(i0+3:i0+n) * rnorm2
          x(n    ,i2,i3) = a(i0+n+1)
        end do
      end do
    case (+1)
      !----------------------
      ! spectral to gridpoint
      !----------------------
      i = 0
      do i3=1,n3
#if defined (__SX__)
!NEC$ ivdep
#endif
        do i2=1,n2
#if defined (__SX__)
          i = (i3-1)*n2 + i2
#else
          i=i+1
#endif
          i0 = (i-1)*(n+2)
          a(i0+1)        = x(1    ,i2,i3)
          a(i0+2)        = 0._wp
          a(i0+3:i0+n)   = x(2:n-1,i2,i3) / rnorm2
          a(i0+n+1)      = x(n    ,i2,i3)
          a(i0+n+2)      = 0._wp
        end do
      end do
      call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
      !----------------------
      i = 0
      do i3=1,n3
#if defined (__SX__)
!NEC$ ivdep
#endif
        do i2=1,n2
#if defined (__SX__)
          i = (i3-1)*n2 + i2
#else
          i=i+1
#endif
          i0 = (i-1)*(n+2)
          x(:,i2,i3) = a(i0+1:i0+n) / rnorm1
        end do
      end do
    end select
    deallocate (a)
    deallocate (work)
  endif
  if(ll2) then
    !------------------
    ! fft along index 2
    !------------------
    n    = size(x,2)
    n2   = size(x,1)
    n3   = size(x,3)
    lot  = n2*n3
    inc  = 1
    jump = n+2
    allocate (a ((n+2)*lot))
    allocate (work ((n+1)*min(lot,nfft)))
    if (nfax(10)/=n) call inifft(n)
    rnorm1 = sqrt(real(n,wp)) ** lno
    rnorm2 = sqrt(2._wp)      ** lno
    select case (isign)
    case (-1)
      !----------------------
      ! gridpoint to spectral
      !----------------------
      i = 0
      do i3=1,n3
#if defined (__SX__)
!NEC$ ivdep
#endif
        do i2=1,n2
#if defined (__SX__)
          i = (i3-1)*n2 + i2
#else
          i=i+1
#endif
          i0 = (i-1)*(n+2)
          a(i0+1:i0+n)     = x(i2,:,i3) * rnorm1
          a(i0+n+1:i0+n+2) = 0._wp
        end do
      end do
      call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
      i = 0
      do i3=1,n3
#if defined (__SX__)
!NEC$ ivdep
#endif
        do i2=1,n2
#if defined (__SX__)
          i = (i3-1)*n2 + i2
#else
          i=i+1
#endif
          i0 = (i-1)*(n+2)
          x(i2,1    ,i3) = a(i0+1)
          x(i2,2:n-1,i3) = a(i0+3:i0+n) * rnorm2
          x(i2,n    ,i3) = a(i0+n+1)
        end do
      end do
    case (+1)
      !----------------------
      ! spectral to gridpoint
      !----------------------
      i = 0
      do i3=1,n3
#if defined (__SX__)
!NEC$ ivdep
#endif
        do i2=1,n2
#if defined (__SX__)
          i = (i3-1)*n2 + i2
#else
          i=i+1
#endif
          i0 = (i-1)*(n+2)
          a(i0+1)      = x(i2,1    ,i3)
          a(i0+2)      = 0._wp
          a(i0+3:i0+n) = x(i2,2:n-1,i3) / rnorm2
          a(i0+n+1)    = x(i2,n    ,i3)
          a(i0+n+2)    = 0._wp
        end do
      end do
      call fft991cy (a,work,trig,nfax,inc,jump,n,lot,isign)
      i = 0
      do i3=1,n3
#if defined (__SX__)
!NEC$ ivdep
#endif
        do i2=1,n2
#if defined (__SX__)
          i = (i3-1)*n2 + i2
#else
          i=i+1
#endif
          i0 = (i-1)*(n+2)
          x(i2,:,i3) = a(i0+1:i0+n) / rnorm1
        end do
      end do
    end select
    deallocate (a)
    deallocate (work)
  endif
end subroutine fft3
!==============================================================================
  subroutine spread_m0_n (sh, f0, nn)
  real (wp) ,intent(out) :: sh (:)
  real (wp) ,intent(in)  :: f0 (:)
  integer   ,intent(in)  :: nn
    integer :: jl, jm, j1, j2
    j2 = 0
    do jl  = 1,2*nn+1
      jm = jl/2+1
      j1=j2+1
      j2=j2+nn+2-jm
      sh(j1  :j2) = f0 (jm:nn+1)
    end do
  end subroutine spread_m0_n
!------------------------------------------------------------------------------
  subroutine spread_m0_l (sh, f0, nn)
  real (wp) ,intent(out) :: sh (:)
  real (wp) ,intent(in)  :: f0 (:)
  integer   ,intent(in)  :: nn
    integer :: jl, jm, j1, j2, jo
    j2 = 0
    jo = 0
    sh = 0._wp
    do jl  = 1,2*nn! +1
!     if(jl<=ni-1) then
      jm = jl/2+1
      j1=j2+1
      j2=j2+nn+2-jm
!print *,'spread_m0_l: jl,jm,j1,j2=',jl,jm,j1,j2
      if(jo/=jm) sh(j1  :j2) = f0 (1:nn+2-jm)
      jo = jm
!     end if
    end do
  end subroutine spread_m0_l
!==============================================================================
  subroutine lgt2 (sh, gp, nn, norm, adj)
  !--------------------------------------------
  ! Legendre transformation routine (2D-fields)
  !--------------------------------------------
  use mo_legendre, only : nm, nmp, nnp,&! Legendre polynominal index descript.
                          pnmd,        &! Modified coef. for direct transform
                          legmod,      &! prepare polynominals (one latitude)
                          gw            ! gaussian weights
  !----------
  ! arguments
  !----------
  real(wp) ,intent(out)          :: sh(:)  ! spherical harmonics field
  real(wp) ,intent(in)           :: gp(:,:)! grid point field
  integer  ,intent(in)           :: nn     ! max meridional wave number for m=0
  integer  ,intent(in) ,optional :: norm   ! normalization flag
  logical  ,intent(in) ,optional :: adj    ! perform adjoint of inv. operation
    !----------------
    ! local variables
    !----------------
    real(wp) :: fsym(size(gp,1))          ! symmetric part of field
    real(wp) :: fasy(size(gp,1))          ! antisymmetric part of field
    integer  :: ni, nj                    ! 'gp' field shape
    real(wp) :: a (size(gp,1),size(gp,2)) ! Fourier transformed field 'gp'
    real(wp) :: rnorm, fac
    integer  :: inorth, isouth, jl, jm, i1,i2, j1,j2, lno
    logical  :: ad
    ni  = size (gp,1)
    nj  = size (gp,2)
    lno = NORM_ECHAM; if (present(norm)) lno = norm
    ad  = .false.   ; if (present(adj))  ad  = adj
    if (ni==1) then
      !----------------------
      ! only m=0 coefficients
      !----------------------
      sh = 0._wp
      call lgt1 (sh, gp(1,:), nn, adj)
    else
      !-------------------------------------
      ! Initialize Legendre transform
      ! up to now only triangular truncation
      !-------------------------------------
      call inileg (im=nn, in=nn, ik=nn, igl=nj)
      rnorm = sqrt(2._wp) ** lno
      !----
      ! fft
      !----
      a = gp
      call fft (a,-1)
      !--------------
      ! latitude loop
      !--------------
      sh  = 0._wp
      fac = 0.5_wp
      do inorth = 1,nj/2
        !------------------------------
        ! prepare legendre polynominals
        !------------------------------
        call legmod (inorth)
        !----------------------------------------
        ! derive symmetric and antisymmetric part
        !----------------------------------------
        isouth = nj+1-inorth
        if (ad) fac = 0.5_wp * ni / gw(inorth)
        fsym = fac * (a(:,inorth) + a(:,isouth))
        fasy = fac * (a(:,inorth) - a(:,isouth))
        !-------------------------------------
        ! longitude (fourier coefficient) loop
        !-------------------------------------
        j2=0
        do jl  = 1,2*nm+1
          if(jl<=ni-1) then
            jm = jl/2+1
            i1=nmp(jm)+1
            i2=nmp(jm)+nnp(jm)
            j1=j2+1
            j2=j2+nnp(jm)
            sh(j1  :j2:2) = sh(j1  :j2:2) + fsym(jl) * pnmd(i1  :i2:2)
            sh(j1+1:j2:2) = sh(j1+1:j2:2) + fasy(jl) * pnmd(i1+1:i2:2)
          end if
        end do
      end do
      sh(nmp(1)+nnp(1)+1:) = sh(nmp(1)+nnp(1)+1:) * rnorm
    endif
  end subroutine lgt2
!------------------------------------------------------------------------------
  subroutine lgt1 (sh, gp, nn, adj)
  !--------------------------------------------
  ! Legendre transformation routine (2D-fields)
  !--------------------------------------------
  use mo_legendre, only : nnp,         &! Legendre polynominal index descript.
                          pnmd,        &! Modified coef. for direct transform
                          legmod,      &! prepare polynominals (one latitude)
                          gw            ! gaussian weights
  !----------
  ! arguments
  !----------
  real(wp) ,intent(out)          :: sh(:)  ! spherical harmonics field
  real(wp) ,intent(in)           :: gp(:)  ! grid point field
  integer  ,intent(in)           :: nn     ! max meridional wave number for m=0
  logical  ,intent(in) ,optional :: adj    ! perform adjoint of inv. operation
    !----------------
    ! local variables
    !----------------
    real(wp) :: fsym                      ! symmetric part of field
    real(wp) :: fasy                      ! antisymmetric part of field
    integer  :: nj                        ! 'gp' field shape
    real(wp) :: fac
    integer  :: inorth, isouth, i2
    logical  :: ad
    nj  = size (gp)
    ad  = .false.; if (present(adj))  ad  = adj
    !-------------------------------------
    ! Initialize Legendre transform
    ! up to now only triangular truncation
    !-------------------------------------
    call inileg (im=nn, in=nn, ik=nn, igl=nj)
    !--------------
    ! latitude loop
    !--------------
    sh  = 0._wp
    fac = 0.5_wp
    do inorth = 1,nj/2
      !------------------------------
      ! prepare legendre polynominals
      !------------------------------
      call legmod (inorth)
      !----------------------------------------
      ! derive symmetric and antisymmetric part
      !----------------------------------------
      isouth = nj+1-inorth
      if (ad) fac = 0.5_wp / gw(inorth)
      fsym = fac * (gp(inorth) + gp(isouth))
      fasy = fac * (gp(inorth) - gp(isouth))
      !-------------------------------------
      ! longitude (fourier coefficient) loop
      !-------------------------------------
      i2=nnp(1)
      sh(1:i2:2) = sh(1:i2:2) + fsym * pnmd(1:i2:2)
      sh(2:i2:2) = sh(2:i2:2) + fasy * pnmd(2:i2:2)
    end do
  end subroutine lgt1
!==============================================================================
  subroutine lgti2 (gp, sh, nn, norm, adj)
  !----------------------------------------------------
  ! Inverse Legendre transformation routine (2D-fields)
  !----------------------------------------------------
  use mo_legendre, only : nm, nmp, nnp,&! Legendre polynominal index descript.
                          pnmi,        &! Modified coef. for inverse transform
                          leginv,      &! prepare polynominals (one latitude)
                          gw            ! gaussian weights
  !----------
  ! arguments
  !----------
  real(wp) ,intent(out)          :: gp(:,:)! grid point field
  real(wp) ,intent(in)           :: sh(:)  ! spherical harmonics field
  integer  ,intent(in)           :: nn     ! max meridional wave number for m=0
  integer  ,intent(in) ,optional :: norm   ! normalization flag
  logical  ,intent(in) ,optional :: adj    ! perform adjoint of dir. operation
    !----------------
    ! local variables
    !----------------
    real(wp) :: fsym(size(gp,1))          ! symmetric part of field
    real(wp) :: fasy(size(gp,1))          ! antisymmetric part of field
    integer  :: ni, nj                    ! 'gp' field shape
    integer  :: inorth, isouth, jl, jm, i1,i2, j1,j2, lno
    real(wp) :: rnorm
    logical  :: ad
    ni  = size (gp,1)
    nj  = size (gp,2)
    lno = NORM_ECHAM; if (present(norm)) lno = norm
    ad  = .false.   ; if (present(adj))  ad  = adj
    if (ni==1) then
      call lgti1 (gp(1,:), sh, nn, adj)
    else
      !-------------------------------------
      ! Initialize Legendre transform
      ! up to now only triangular truncation
      !-------------------------------------
      call inileg (im=nn, in=nn, ik=nn, igl=nj)
      rnorm = sqrt(2._wp) ** lno
      !--------------
      ! latitude loop
      !--------------
      do inorth = 1,nj/2
        !------------------------------
        ! prepare legendre polynominals
        !------------------------------
        call leginv (inorth)
        !-------------------------------------
        ! longitude (fourier coefficient) loop
        !-------------------------------------
        j2=0
        fsym = 0._wp; fasy = 0._wp
        do jl  = 1,2*nm+1
          if(jl<=ni-1) then
            jm = jl/2+1
            i1=nmp(jm)+1
            i2=nmp(jm)+nnp(jm)
            j1=j2+1
            j2=j2+nnp(jm)
            fsym(jl) = sum (pnmi(i1  :i2:2) * sh(j1  :j2:2))
            fasy(jl) = sum (pnmi(i1+1:i2:2) * sh(j1+1:j2:2))
          endif
        end do
        !----------------------------------------
        ! derive symmetric and antisymmetric part
        !----------------------------------------
        isouth = nj+1-inorth
        if (ad) then
          gp(:,inorth) = gw(inorth)/ni * (fsym + fasy)
          gp(:,isouth) = gw(inorth)/ni * (fsym - fasy)
        else
          gp(:,inorth) = fsym + fasy
          gp(:,isouth) = fsym - fasy
        endif
      end do
      !----
      ! fft
      !----
      gp(2:,:) = gp(2:,:) / rnorm
      call fft(gp,1)
    endif
  end subroutine lgti2
!==============================================================================
  subroutine lgti1 (gp, sh, nn, adj)
  !----------------------------------------------------
  ! Inverse Legendre transformation routine (2D-fields)
  !----------------------------------------------------
  use mo_legendre, only : nmp, nnp,    &! Legendre polynominal index descript.
                          pnmi,        &! Modified coef. for inverse transform
                          leginv,      &! prepare polynominals (one latitude)
                          gw            ! gaussian weights
  !----------
  ! arguments
  !----------
  real(wp) ,intent(out)          :: gp(:)  ! grid point field
  real(wp) ,intent(in)           :: sh(:)  ! spherical harmonics field
  integer  ,intent(in)           :: nn     ! max meridional wave number for m=0
  logical  ,intent(in) ,optional :: adj    ! perform adjoint of dir. operation
    !----------------
    ! local variables
    !----------------
    real(wp) :: fsym                      ! symmetric part of field
    real(wp) :: fasy                      ! antisymmetric part of field
    integer  :: nj                        ! 'gp' field shape
    integer  :: inorth, isouth, i1, i2, j2
    logical  :: ad
    nj  = size (gp,1)
    ad  = .false.; if (present(adj))  ad  = adj
    !-------------------------------------
    ! Initialize Legendre transform
    ! up to now only triangular truncation
    !-------------------------------------
    call inileg (im=nn, in=nn, ik=nn, igl=nj)
    !--------------
    ! latitude loop
    !--------------
    do inorth = 1,nj/2
      !------------------------------
      ! prepare legendre polynominals
      !------------------------------
      call leginv (inorth)
      !-------------------------------------
      ! longitude (fourier coefficient) loop
      !-------------------------------------
      i1=nmp(1)+1
      i2=nmp(1)+nnp(1)
      j2=nnp(1)
      fsym = sum (pnmi(i1  :i2:2) * sh(1  :j2:2))
      fasy = sum (pnmi(i1+1:i2:2) * sh(1+1:j2:2))
      !----------------------------------------
      ! derive symmetric and antisymmetric part
      !----------------------------------------
      isouth = nj+1-inorth
      if (ad) then
        gp(inorth) = gw(inorth) * (fsym + fasy)
        gp(isouth) = gw(inorth) * (fsym - fasy)
      else
        gp(inorth) = fsym + fasy
        gp(isouth) = fsym - fasy
      endif
    end do
  end subroutine lgti1
!==============================================================================
  subroutine lgt3 (sh, gp, nn, norm, adj)
  !---------------------------------------------
  ! Legendre transformation routine (3-D fields)
  !---------------------------------------------
  real(wp) ,intent(out)          :: sh(:,:)  ! spherical harmonics field
  real(wp) ,intent(in)           :: gp(:,:,:)! grid point field
  integer  ,intent(in)           :: nn       ! max meridional wave number m=0
  integer  ,intent(in) ,optional :: norm     ! normalization flag
  logical  ,intent(in) ,optional :: adj      ! perform adjoint of inv. op.
    integer :: k
    do k=1,size(gp,3)
      call lgtd (sh(:,k), gp(:,:,k), nn, norm, adj)
    end do
  end subroutine lgt3
!==============================================================================
  subroutine lgti3 (gp, sh, nn, norm, adj)
  !----------------------------------------------------
  ! Inverse Legendre transformation routine (3D-fields)
  !----------------------------------------------------
  real(wp) ,intent(out)          :: gp(:,:,:)! grid point field
  real(wp) ,intent(in)           :: sh(:,:)  ! spherical harmonics field
  integer  ,intent(in)           :: nn       ! max meridional wave number (m=0)
  integer  ,intent(in) ,optional :: norm     ! normalization flag
  logical  ,intent(in) ,optional :: adj      ! perform adjoint of dir. op.
    integer :: k
    do k=1,size(gp,3)
      call lgti (gp(:,:,k), sh(:,k), nn, norm, adj)
    end do
  end subroutine lgti3
!==============================================================================
  subroutine fft2_adj (x, isign, l1, l2, norm)
  !----------------------------------------------------------------------------
  ! Adjoint routine to fft2.
  !
  ! Method:
  ! Fft is a orthogonal operation. Thus the adjoint routine is the
  ! inverse routine (despite some normalization factors).
  !----------------------------------------------------------------------------
  real(wp) ,intent(inout)         :: x(:,:) ! array to transform
  integer  ,intent(in)            :: isign  ! +1: spectral -> gridpoint, -1: <-
  logical  ,intent(in)  ,optional :: l1     ! (def=.true.)  transform 1st index
  logical  ,intent(in)  ,optional :: l2     ! (def=.false.) transform 2nd index
  integer  ,intent(in)  ,optional :: norm   ! normalization flag
    integer :: lno
    lno = NORM_ECHAM_AD; if(present(norm)) lno=lno-norm
    call fft(x, -isign, l1=l1, l2=l2, norm=lno)
  end subroutine fft2_adj
!------------------------------------------------------------------------------
  subroutine fftc_adj (x, c, isign, norm)
  !----------------------------------------------------------------------------
  ! Adjoint routine to fftc.
  !
  ! Method:
  ! Fft is a orthogonal operation. Thus the adjoint routine is the
  ! inverse routine (despite some normalization factors).
  !----------------------------------------------------------------------------
  real(wp)    ,intent(inout)        :: x(:,:) ! array to transform
  complex(wp) ,intent(inout)        :: c(:,:) ! array to transform
  integer     ,intent(in)           :: isign  ! +1:spectral -> gridpoint, -1:<-
  integer     ,intent(in) ,optional :: norm   ! normalization flag
    integer :: lno
    lno = NORM_ECHAM_AD; if(present(norm)) lno=lno-norm
    call fftc(x, c, -isign, norm=lno)
  end subroutine fftc_adj
!==============================================================================
  subroutine fft3_adj (x, isign, l1, l2, norm)
  !----------------------------------------------------------------------------
  ! Adjoint routine to fft.
  !
  ! Method:
  ! Fft is a orthogonal operation. Thus the adjoint routine is the
  ! inverse routine (despite some normalization factors).
  !----------------------------------------------------------------------------
  real(wp) ,intent(inout)         :: x(:,:,:)! array to transform
  integer  ,intent(in)            :: isign   !+1: spectral -> gridpoint, -1: <-
  logical  ,intent(in)  ,optional :: l1      !(def=.true.)  transform 1st index
  logical  ,intent(in)  ,optional :: l2      !(def=.false.) transform 2nd index
  integer  ,intent(in)  ,optional :: norm    ! normalization flag
    integer :: lno
    lno = NORM_ECHAM_AD; if(present(norm)) lno=lno-norm
    call fft(x, -isign, l1=l1, l2=l2, norm=lno)
  end subroutine fft3_adj
!==============================================================================
  subroutine lgt2_adj (sh, gp, nn, norm)
  !----------------------------------------------------------------------------
  ! Adjoint Legendre transformation routine (2D-field)
  !
  ! Method:
  ! Lgt is a orthogonal operation. Thus the adjoint routine is the
  ! inverse routine lgti (despite some normalization factors).
  !----------------------------------------------------------------------------
  real(wp) ,intent(in)           :: sh(:)  ! spherical harmonics field
  real(wp) ,intent(out)          :: gp(:,:)! grid point field
  integer  ,intent(in)           :: nn     ! max meridional wave number for m=0
  integer  ,intent(in) ,optional :: norm   ! normalization flag
    integer :: lno
    lno = NORM_ECHAM_AD; if(present(norm)) lno=lno-norm
    call lgti (gp, sh, nn, lno, .true.)
  end subroutine lgt2_adj
!==============================================================================
  subroutine lgti2_adj (gp, sh, nn, norm)
  !----------------------------------------------------------------------------
  ! Adjoint inverse Legendre transformation routine (2D-field)
  !
  ! Method:
  ! Lgti is a orthogonal operation. Thus the adjoint routine is the
  ! inverse routine lgt (despite some normalization factors).
  !----------------------------------------------------------------------------
  real(wp) ,intent(in)           :: gp(:,:)! grid point field
  real(wp) ,intent(out)          :: sh(:)  ! spherical harmonics field
  integer  ,intent(in)           :: nn     ! max meridional wave number for m=0
  integer  ,intent(in) ,optional :: norm   ! normalization flag
    integer :: lno
    lno = NORM_ECHAM_AD; if(present(norm)) lno=lno-norm
    call lgtd (sh, gp, nn, lno, .true.)
  end subroutine lgti2_adj
!==============================================================================
  subroutine lgt3_adj (sh, gp, nn, norm)
  !----------------------------------------------------------------------------
  ! Adjoint Legendre transformation routine (3D-field)
  !----------------------------------------------------------------------------
  real(wp) ,intent(in)           :: sh(:,:)  ! spherical harmonics field
  real(wp) ,intent(out)          :: gp(:,:,:)! grid point field
  integer  ,intent(in)           :: nn       ! max meridional wave number (m=0)
  integer  ,intent(in) ,optional :: norm     ! normalization flag
    integer :: lno
    lno = NORM_ECHAM_AD; if(present(norm)) lno=lno-norm
    call lgti (gp, sh, nn, lno, .true.)
  end subroutine lgt3_adj
!==============================================================================
  subroutine lgti3_adj (gp, sh, nn, norm)
  !----------------------------------------------------------------------------
  ! Adjoint inverse Legendre transformation routine (3D-field)
  !----------------------------------------------------------------------------
  real(wp) ,intent(in)           :: gp(:,:,:)! grid point field
  real(wp) ,intent(out)          :: sh(:,:)  ! spherical harmonics field
  integer  ,intent(in)           :: nn       ! max meridional wave number (m=0)
  integer  ,intent(in) ,optional :: norm     ! normalization flag
    integer :: lno
    lno = NORM_ECHAM_AD; if(present(norm)) lno=lno-norm
    call lgtd (sh, gp, nn, lno, .true.)
  end subroutine lgti3_adj
!==============================================================================
  function nsp_nn (nn)
  !-------------------------------------------------------
  ! derive number of spectral coefficients from truncation
  ! (only triangular truncation jet)
  !-------------------------------------------------------
  integer, intent(in) :: nn
  integer             :: nsp_nn
    nsp_nn = (nn+1)**2
  end function nsp_nn
!------------------------------------------------------------------------------
  subroutine gauaw (ga, gw)
  real(wp), intent(out), optional :: ga (:)
  real(wp), intent(out), optional :: gw (:)
    integer               :: ngl
    real(wp) ,allocatable :: zga (:)
    real(wp) ,allocatable :: zgw (:)
    if (present(ga)) then
      ngl = size(ga)
    else if (present(gw)) then
      ngl = size(gw)
    else
      return
    endif
    allocate (zga (ngl), zgw(ngl))
    call gaw (zga, zgw, ngl)
    if (present(ga)) ga = zga
    if (present(gw)) gw = zgw
  end subroutine gauaw
!==============================================================================
end module mo_transform
