!
!+ Provides numerical algorithms (matrix inversion, SVD, index search, etc)
!
MODULE mo_algorithms
!
! Description:
!   Provides numerical algorithms (matrix inversion, SVD, index search, etc)
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
! V1_8         2009/12/09 Andreas Rhodin
!  splint_vec: remove if clause (x(i,k) <= xa(i,kdimin)) (same as splint)
! V1_9         2010/04/20 Harald Anlauf
!  Fallback implementation of Gaussian error function and companions
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_15        2011/12/06 Harald Anlauf
!  splint_vec: remove optional dummy argument 'sequential'; some cleanups
! V1_27        2013-11-08 Harald Anlauf
!  Workaround for crayftn optimization problem
! V1_35        2014-11-07 Robin Faulwetter
!  new subroutine fit_polynomial
! V1_42        2015-06-08 Harald Anlauf
!  OpenMP parallelization of vintp_spline, splint_vec
!  init_splinex_vec: use OpenMP workshare
! V1_45        2015-12-15 Harald Anlauf
!  minor optimization for OpenMP; Add 'index' function for integer(i8) arrays
! V1_46        2016-02-05 Michael Bender
!  some routines moved from the STD operator to mo_physics and mo_algorithms
! V1_51        2017-02-24 Harald Anlauf
!  init_splinex_vec: optimize OpenMP parallel region/workshare
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  1999       Routines gathered from different sources
! Andreas Rhodin  DWD  2000-2008
! Harald Anlauf   DWD  2007-2008
!=============================================================================
  !-------------
  ! Modules used
  !-------------
  use mo_kind      , only: wp, sp, dp, i8
  use mo_exception , only: finish
  use slatec_module, only: rs
  use slatec_module, only: sort
  use mo_random    , only: random_number,&! D.E.Knuth's Random number generator
                           random_static,&! Random number gener. without seed
                           random_state_t ! Derived type:random generator state
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: simplex        ! simplex minimisation
  public :: jacobi         ! eigenval.+vectors(sym.matrix)
  public :: eigen          ! eigenval.+vectors(sym.matrix)
  public :: index          ! index an array in asc. order
  public :: xgesvd         ! singular value decomposition (lapack)
  public :: random_gauss   ! returns normally distributed deviate
  public :: hunt           ! finds a location in a table of monotonic numbers.
  public :: init_splinex   ! set coefficients (2nd deriv.) for spline interpol.
  public :: splint         ! Spline interpol.with optimization for subseq.calls
  public :: splint_vec     ! Spline interpol. (vectorizing version)
  public :: init_spline2   ! 2nd derivatives for ...
  public :: vintp_spline   ! ... spline interpolation (atm comp.)
  public :: polint         ! Polynomial interpolation (Neville's algorithm)
  public :: fit_polynomial ! fit polinominal to some function
  public :: erf            ! returns the error function erf(x)
  public :: gammln         ! returns the value ln(Gamma(xx)) for xx > 0
  public :: init_splinex_vec ! Vectorized version of init_splinex
  public :: integpolycube    ! numerical integration
  public :: integpolycube_tl ! -> tangent-linear function
  public :: integpolycube_ad ! -> adjoint routine
  public :: binsearch        ! binary search in sorted list
  public :: runsort          ! sort array
  public :: mergearr         ! merge two sorted arrays
  !------------------------------------------------------------------
  ! Fallback implementation of Gaussian error function and companions
  !------------------------------------------------------------------
  public :: derf         ! (Double precision) Error function
  public :: derfc        ! (Double precision) Complementary error function
  public :: derfcx       ! (Double precision) derfc, scaled version

!-----------------------------------
! provide scalar and vector versions
!-----------------------------------
  interface random_gauss
    module procedure random_gauss0
    module procedure random_gauss1
    module procedure random_gauss2
    module procedure random_gauss3
  end interface
!------------------------------------------------------------------------------
  interface hunt
    module procedure hunt1
    module procedure huntn
  end interface
!------------------------------------------------------------------------------
  interface index
    module procedure indexsp, indexdp, indexi, indexi8, indexc
  end interface
!------------------------------------------------------------------------------
  interface init_splinex
    module procedure init_splinex
  end interface init_splinex
!------------------------------------------------------------------------------
  interface xgesvd
    module procedure xgesvd_r   ! Real    SVD wrapper
    module procedure xgesvd_c   ! Complex SVD wrapper
  end interface
!------------------------------------------------------------------------------
  interface binsearch
    module procedure binsearchi   ! Integer search
    module procedure binsearchi8  ! Integer search
    module procedure binsearchwp  ! Real search
  end interface binsearch
!------------------------------------------------------------------------------
  interface runsort
     module procedure runsorti   ! Sort integer array
     module procedure runsorti8  ! Sort integer array
     module procedure runsortc   ! Sort character array
  end interface runsort
!------------------------------------------------------------------------------
  interface mergearr
     module procedure mergearri   ! Merge two sorted integer arrays
     module procedure mergearri8  ! Merge two sorted integer arrays
     module procedure mergearrc   ! Merge two sorted character arrays
  end interface mergearr
!------------------------------------------------------------------------------

!--------------------------------------------
! saved variables for subroutine random_gauss
!--------------------------------------------
  logical  ,save :: valid = .false.
  real(wp) ,save :: gaus2

!==============================================================================
contains
!==============================================================================
  subroutine huntn(xx,x,jlo)
  real(wp) ,intent(in)  :: xx(:), x(:)
  integer  ,intent(out) :: jlo(size(x))
    integer :: i, j
    j = 1
    do i = 1, size(x)
      call hunt1(xx(:), x(i), j)
      jlo(i) = j
    end do
  end subroutine huntn
!------------------------------------------------------------------------------
  subroutine hunt1(xx,x,jlo)
  real(wp) ,intent(in)    :: xx(:), x
  integer  ,intent(inout) :: jlo
  !---------------------------------------------------------------------------
  ! Given an array XX(1:N), and given a value X, returns a value JLO such
  ! that x is between XX(JLO) and XX(JLO+1). XX(1:N) must be monotonic, either
  ! increasing or decreasing. JLO=0 or JLO=N is returned to indicate that x is
  !  out of range. JLO on input is taken as the initial guess for JLO on
  ! output.
  !---------------------------------------------------------------------------
    integer :: inc,jhi,jm,n
    logical :: ascnd
      n     = size(xx)
      ascnd = xx(n) >= xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
      else
        inc=1
        if(x.ge.xx(jlo).eqv.ascnd)then
          do
            jhi=jlo+inc
            if(jhi.gt.n)then
              jhi=n+1
            else if(x.ge.xx(jhi).eqv.ascnd)then
              jlo=jhi
              inc=inc+inc
              cycle
            endif
            exit
          end do
        else
          jhi=jlo
          do
            jlo=jhi-inc
            if(jlo.lt.1)then
              jlo=0
            else if(x.lt.xx(jlo).eqv.ascnd)then
              jhi=jlo
              inc=inc+inc
              cycle
            endif
            exit
          end do
        endif
      endif
      do
        if(jhi-jlo.eq.1)then
          if(x.eq.xx(n))jlo=n-1
          if(x.eq.xx(1))jlo=1
          return
        endif
        jm=(jhi+jlo)/2
        if(x.ge.xx(jm).eqv.ascnd)then
          jlo=jm
        else
          jhi=jm
        endif
      end do
  end subroutine hunt1
!==============================================================================
  function indexi(x) result (ix)
  integer,  intent(in)  :: x  (:)
  integer               :: ix (size(x))
  !--------------------------------------------------------------
  ! Indexes an array X, i.e. outputs the array INDEX such that
  ! X(INDEX(J)) is in ascending order for J=1,2,..
  !--------------------------------------------------------------
    integer :: ier
    call sort (x, ix, kflag=1, ier=ier)
    if (ier/=0) call finish ('indexi','ier/=0')
  end function indexi
!------------------------------------------------------------------------------
  function indexi8(x) result (ix)
  integer(i8),  intent(in)  :: x  (:)
  integer                   :: ix (size(x))
  !--------------------------------------------------------------
  ! Indexes an array X, i.e. outputs the array INDEX such that
  ! X(INDEX(J)) is in ascending order for J=1,2,..
  !--------------------------------------------------------------
    integer :: ier
    call sort (x, ix, kflag=1, ier=ier)
    if (ier/=0) call finish ('indexi','ier/=0')
  end function indexi8
!------------------------------------------------------------------------------
  function indexc(x) result (ix)
  character(len=*), intent(in) :: x  (:)
  integer                      :: ix (size(x))
  !--------------------------------------------------------------
  ! Indexes an array X, i.e. outputs the array INDEX such that
  ! X(INDEX(J)) is in ascending order for J=1,2,..
  !--------------------------------------------------------------
    integer :: ier
    call sort (x, ix, kflag=1, ier=ier)
    if (ier/=0) call finish ('indexc','ier/=0')
  end function indexc
!------------------------------------------------------------------------------
  function indexdp(x) result (ix)
  real(dp), intent(in)  :: x  (:)
  integer               :: ix (size(x))
  !--------------------------------------------------------------
  ! Indexes an array X, i.e. outputs the array INDEX such that
  ! X(INDEX(J)) is in ascending order for J=1,2,..
  !--------------------------------------------------------------
    integer :: ier
    call sort (x, ix, kflag=1, ier=ier)
    if (ier/=0) call finish ('indexdp','ier/=0')
  end function indexdp
!------------------------------------------------------------------------------
  function indexsp(x) result (ix)
  real(sp), intent(in)  :: x  (:)
  integer               :: ix (size(x))
  !--------------------------------------------------------------
  ! Indexes an array X, i.e. outputs the array INDEX such that
  ! X(INDEX(J)) is in ascending order for J=1,2,..
  !--------------------------------------------------------------
    integer :: ier
    call sort (x, ix, kflag=1, ier=ier)
    if (ier/=0) call finish ('indexsp','ier/=0')
  end function indexsp
!==============================================================================
  subroutine xgesvd_r (a,w,v)     ! Real LAPACK SVD driver
  real(wp), intent(inout)         :: a(:,:)
  real(wp), intent(out)           :: w(:)
  real(wp), intent(out)           :: v(:,:)
  !--------------------------------------------------------------------------
  ! Given a matrix A, with dimensions M by N, this routine computes its
  ! singular value decomposition, A=U.W.Vt. The Matrix U replaces A on
  ! output. The diagonal matrix of singular values W is output as a vector
  ! W. The matrix V (not the transpose Vt) is output as V.
  !--------------------------------------------------------------------------
    real(wp) ,allocatable :: work  (:)                ! work array
    real(wp)              :: work1 (1)                ! dummy work array
    real(wp)              :: vt (size(v,2),size(v,1)) ! v transposed
    real(wp)              :: u  (1)                   ! dummy, U is in A(out)
    integer               :: m,n,mn,info,i            ! dimensions, error flag
    integer               :: lwork                    ! work array size
    integer               :: ier                      ! error flag
    integer ,save         :: msave = -1, nsave = -1, lwsave
#if defined (__SX__)
    real(wp) ,allocatable :: a_(:,:)                  ! temporary copy of a
    integer               :: lda_                     ! leading dim. of a_
#endif
    !------------------------------------
    ! determine work array size, allocate
    !------------------------------------
    m  = size (a,1)
    n  = size (a,2)
    mn = min  (m,n)
    if (mn == 0) return
#if defined (__SX__)
    !--------------------------------------------------------
    ! Tune temporary array to reduce bank conflicts on NEC SX
    !--------------------------------------------------------
    lda_ = (m/2)*2 + 1
    allocate (a_(lda_,n))
#endif
    if (m == msave .and. n == nsave) then
      lwork = lwsave
    else
      lwork = -1
#if defined (__SX__)
      call dgesvd ('O','S',m,n,a_,lda_,&
                                   w,u,1,vt,size(vt,1),work1,lwork,info)
#else
      call dgesvd ('O','S',m,n,a,m,w,u,1,vt,size(vt,1),work1,lwork,info)
#endif
      lwork  = work1(1)
      nsave  = n
      msave  = m
      lwsave = lwork
    endif
    allocate (work(lwork), stat=ier)                ! optimum size
    if (ier /= 0) then
      lwork = max(3*mn+max(m,n),5*mn-4)
      allocate (work(lwork))                        ! minimum size
    endif
    !--------------------
    ! call LAPACK routine
    !--------------------
#if defined (__SX__)
    a_(1:m,:) = a(1:m,:)
    call dgesvd ('O','S',m,n,a_,lda_,&
                                 w,u,1,vt,size(vt,1),work, lwork,info)
    a(1:m,:mn) = a_(1:m,:mn)                        ! extract columns of U
#else
    call dgesvd ('O','S',m,n,a,m,w,u,1,vt,size(vt,1),work, lwork,info)
#endif
    deallocate (work)
    !------------
    ! error check
    !------------
    select case (info)
    case (0)  ! sucessful exit
    case (:-1)
      call finish ('xgesvd_r','invalid argument')
    case (1:)
      call finish ('xgesvd_r','dbdsqr did not converge')
    end select
    !-----------------
    ! return v, not vt
    !-----------------
     do i=1,size(v,1)
       v(i,:mn) = vt(:mn,i)
     end do
  end subroutine xgesvd_r
!==============================================================================
  subroutine xgesvd_c (a,w,v)           ! Complex LAPACK SVD driver (zgesvd)
    complex(wp), intent(inout)  :: a(:,:)
    real   (wp), intent(out)    :: w(:)
    complex(wp), intent(out)    :: v(:,:)
    !--------------------------------------------------------------------------
    ! Given a complex matrix A, with dimensions M by N, this routine computes
    ! its singular value decomposition, A=U.W.V^H.  The Matrix U replaces A
    ! on output.  The diagonal matrix of singular values W is output as vector
    ! W.  The matrix V (not the Hermitian conjugate V^H) is output as V.
    !--------------------------------------------------------------------------
    complex(wp) ,allocatable :: work (:)                ! work array
    complex(wp)              :: work1(1)                ! dummy work array
    real(wp)    ,allocatable :: rwork(:)                ! work array
    complex(wp)              :: vt(size(v,2),size(v,1)) ! v transposed
    real(wp)                 :: u    (1)                ! dummy, U is in A(out)
    integer                  :: m,n,mn,info,i           ! dims., error flag
    integer                  :: lwork                   ! work array size
    integer                  :: ier                     ! error flag
    integer ,save            :: msave = -1, nsave = -1, lwsave
#if defined (__SX__)
    complex(wp) ,allocatable :: a_(:,:)                 ! temporary copy of a
    integer                  :: lda_                    ! leading dim. of a_
#endif
    !------------------------------------
    ! Determine work array size, allocate
    !------------------------------------
    m  = size (a,1)
    n  = size (a,2)
    mn = min  (m,n)
    if (mn == 0) return
    allocate (rwork(5*mn))
#if defined (__SX__)
    !--------------------------------------------------------
    ! Tune temporary array to reduce bank conflicts on NEC SX
    !--------------------------------------------------------
    lda_ = (m/2)*2 + 1
    allocate (a_(lda_,n))
#endif
    if (m == msave .and. n == nsave) then
       lwork = lwsave
    else
       lwork = -1
#if defined (__SX__)
       call zgesvd ('O','S',m,n,a_,lda_,&
                                    w,u,1,vt,size(vt,1),work1,lwork,rwork,info)
#else
       call zgesvd ('O','S',m,n,a,m,w,u,1,vt,size(vt,1),work1,lwork,rwork,info)
#endif
       lwork  = work1(1)
       nsave  = n
       msave  = m
       lwsave = lwork
    endif
    allocate (work(lwork), stat=ier)                ! optimum size
    if (ier /= 0) then
       lwork = 2*mn + max (m,n)
       allocate (work(lwork))                       ! minimum size
    endif
    !--------------------
    ! Call LAPACK routine
    !--------------------
#if defined (__SX__)
    a_(1:m,:) = a(1:m,:)
    call zgesvd ('O','S',m,n,a_,lda_,&
                                 w,u,1,vt,size(vt,1),work,lwork,rwork,info)
    a(1:m,:mn) = a_(1:m,:mn)                        ! extract columns of U
#else
    call zgesvd ('O','S',m,n,a,m,w,u,1,vt,size(vt,1),work,lwork,rwork,info)
#endif
    deallocate (work)
    deallocate (rwork)
    !------------
    ! Error check
    !------------
    select case (info)
    case (0)  ! Sucessful exit
    case (:-1)
       call finish ('xgesvd_c','invalid argument')
    case (1:)
       call finish ('xgesvd_c','zbdsqr did not converge')
    end select
    !------------------
    ! Return v, not v^H
    !------------------
    do i=1, size (v,1)
       v(i,:mn) = conjg (vt(:mn,i))
    end do
  end subroutine xgesvd_c
!==============================================================================
  subroutine eigen (a,d,v)
  real(wp),          intent(in)    :: a(:,:)
  real(wp),          intent(out)   :: d(size(a,1))
  real(wp),          intent(out)   :: v(size(a,1),size(a,1))
  !------------------------------------------------------------------------
  ! Computes all eigenvalues and eigenvectors of a real symmetric matrix A.
  ! D returnes the eigenvalues of A.
  ! The columns of V contain on output the normalized eigenvectors of A.
  !-------------------------------------------------------------------------
    call rs(a,d,v)
  end subroutine eigen
!------------------------------------------------------------------------------
  subroutine jacobi (a,d,v)
  real(wp),          intent(in)    :: a(:,:)
  real(wp),          intent(out)   :: d(size(a,1))
  real(wp),          intent(out)   :: v(size(a,1),size(a,1))
  !------------------------------------------------------------------------
  ! Computes all eigenvalues and eigenvectors of a real symmetric matrix A.
  ! D returnes the eigenvalues of A.
  ! The columns of V contain on output the normalized eigenvectors of A.
  !-------------------------------------------------------------------------
    call rs(a,d,v)
  end subroutine jacobi
!==============================================================================
  subroutine random_gauss0 (gauss, seed)
  real(wp)             ,intent(out)             :: gauss
  type(random_state_t) ,intent(inout) ,optional :: seed
  !-----------------------------------------------------------------------
  ! Returns a normal distributed deviate with zero mean and unit variance.
  !-----------------------------------------------------------------------
    real(wp)       :: f, rr, r1, r2
    if (valid) then
      gauss = gaus2
      valid = .false.
    else
      do
        if (present(seed)) then
          call random_number (r1, seed)
        else
          call random_static (r1)
        endif
        r1 = 2._wp * r1 - 1._wp
        if (present(seed)) then
          call random_number (r2, seed)
        else
          call random_static (r2)
        endif
        r2 = 2._wp * r2 - 1._wp
        rr = r1**2 + r2**2
        if (rr >= 1._wp .or. rr == 0._wp) cycle
        exit
      end do
      f     = sqrt ( -2._wp * log(rr) / rr )
      gaus2 = r1 * f
      gauss = r2 * f
      valid = .true.
    endif
  end subroutine random_gauss0
!------------------------------------------------------------------------------
  subroutine random_gauss1 (gauss, seed)
  real(wp)             ,intent(out)             :: gauss(:)
  type(random_state_t) ,intent(inout) ,optional :: seed
  !-----------------------------------------------------------------------
  ! Returns a normal distributed deviate with zero mean and unit variance.
  !-----------------------------------------------------------------------
    integer :: i
    do i = 1, size (gauss)
      call random_gauss0 (gauss(i), seed)
    end do
  end subroutine random_gauss1
!------------------------------------------------------------------------------
  subroutine random_gauss2 (gauss, seed)
  real(wp)             ,intent(out)             :: gauss(:,:)
  type(random_state_t) ,intent(inout) ,optional :: seed
  !-----------------------------------------------------------------------
  ! Returns a normal distributed deviate with zero mean and unit variance.
  !-----------------------------------------------------------------------
    integer :: i
    do i = 1, size (gauss,2)
      call random_gauss1 (gauss(:,i), seed)
    end do
  end subroutine random_gauss2
!------------------------------------------------------------------------------
  subroutine random_gauss3 (gauss, seed)
  real(wp)             ,intent(out)             :: gauss(:,:,:)
  type(random_state_t) ,intent(inout) ,optional :: seed
  !-----------------------------------------------------------------------
  ! Returns a normal distributed deviate with zero mean and unit variance.
  !-----------------------------------------------------------------------
    integer :: i
    do i = 1, size (gauss,3)
      call random_gauss2 (gauss(:,:,i), seed)
    end do
  end subroutine random_gauss3
!==============================================================================
  subroutine init_splinex (x,y,y2,yp1,ypn)
  real(wp) ,intent(in)           :: x (:) ! abscissa
  real(wp) ,intent(in)           :: y (:) ! function to interpolate
  real(wp) ,intent(out)          :: y2(:) ! second derivative of function
  real(wp) ,intent(in) ,optional :: yp1   ! first derivative at x(1)
  real(wp) ,intent(in) ,optional :: ypn   ! first derivative at x(n)
  !-----------------------------------------------------------------------
  ! Given arrays X(1:N), Y(1:N) containing a tabulated function, i.e.,
  ! y_i=f(x_i), with x_1..x_n being monotonous, and given values YP1 and
  ! YPN for the first derivative of the interpolating function at points 1
  ! and n, respectively, this routine returns an array Y2(1:N) of length
  ! N=size(X) which contains the second derivatives of the interpolating
  ! function at the tabulated points x_i. If YP1 and/or YPN are equal to
  ! 1.e30 or larger or not present, the routine is signaled to set the
  ! corresponding boundary condition for a natural spline, with zero second
  ! derivative on that boundary.
  !-----------------------------------------------------------------------

    integer i,k,n
    real(wp) p, qn, sig, un, u(size(x))
    logical nat1, natn

    n = size(x)
    nat1 = .true.; natn = .true.
    if (present(yp1)) then
      if (yp1 < 1.e30_wp) nat1 = .false.
    endif
    if (present(ypn)) then
      if (ypn < 1.e30_wp) natn = .false.
    endif

    if (nat1) then
      y2(1)=0._wp
      u(1)=0._wp
    else
      y2(1)=-0.5_wp
      u(1)=(3._wp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif

    do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2._wp
      y2(i)=(sig-1._wp)/p
      u(i)=(6._wp*((y(i+1)-y(i))/(x(i+1)-x(i))                  &
                  -(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1)) &
                  -sig*u(i-1))/p
    end do

    if (natn) then
      qn=0._wp
      un=0._wp
    else
      qn=0.5_wp
      un=(3._wp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif

    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1._wp)
    do k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
    end do
  end subroutine init_splinex
!------------------------------------------------------------------------------
  subroutine init_splinex_vec (x, y, y2, np, n, yp1, ypn)
    integer  ,intent(in)           :: np        ! Number of points (columns)
    integer  ,intent(in)           :: n         ! No. of function values/point
    real(wp) ,intent(in)           :: x  (np,n) ! abscissa
    real(wp) ,intent(in)           :: y  (np,n) ! function to interpolate
    real(wp) ,intent(out)          :: y2 (np,n) ! second derivative of function
    real(wp) ,intent(in) ,optional :: yp1(np)   ! first derivative at x(:,1)
    real(wp) ,intent(in) ,optional :: ypn(np)   ! first derivative at x(:,n)
    !------------------------------------
    ! Vectorized version of init_splinex.
    !------------------------------------
    integer  :: i, j, k
#ifdef _OPENMP                          /* prefer large arrays on heap */
    real(wp), allocatable :: u(:,:)
#else
    real(wp) :: u(np,n-1)
#endif
    real(wp) :: un(np), p(np), sig(np), qn(np)
    logical  :: nat1(np), natn(np)

#ifdef _OPENMP
    allocate (u(np,n-1))
#endif

    nat1 = .true.; natn = .true.
    if (present (yp1)) then
       nat1 = (yp1 >= 1.e30_wp)
    endif
    if (present (ypn)) then
       natn = (ypn >= 1.e30_wp)
    endif

!---------------------------------------------------------------------
! Work around optimization bug with missing optional arguments (sxf90,crayftn)
!NEC$ nomove
!NEC$ NOMOVE
#ifdef _CRAYFTN
!DIR$ NEXTSCALAR
#endif
!---------------------------------------------------------------------
    do j = 1, np
       if (nat1(j)) then
          y2(j,1) = 0._wp
          u (j,1) = 0._wp
       else
          y2(j,1) = -0.5_wp
          u (j,1) = (3._wp/(x(j,2)-x(j,1))) * &
                    ((y(j,2)-y(j,1))/(x(j,2)-x(j,1))-yp1(j))
       endif
    end do

!$omp parallel private(sig,p)
!NEC$ novector
    do i = 2,n-1
!$omp workshare
       sig(:) = (x(:,i)-x(:,i-1))/(x(:,i+1)-x(:,i-1))
       p(:)   = sig(:)*y2(:,i-1)+2._wp
       y2(:,i)=(sig(:)-1._wp)/p(:)
       u (:,i)=(6._wp*((y(:,i+1)-y(:,i))/(x(:,i+1)-x(:,i))    &
                      -(y(:,i)-y(:,i-1))/(x(:,i)-x(:,i-1))) / &
                       (x(:,i+1)-x(:,i-1))                    &
                - sig(:)*u(:,i-1))/p(:)
!$omp end workshare nowait
    end do
!$omp end parallel

!---------------------------------------------------------------------
! Work around optimization bug with missing optional arguments on SX-6
!NEC$ nomove
!NEC$ NOMOVE
!---------------------------------------------------------------------
    do j = 1, np
       if (natn(j)) then
          qn(j) = 0._wp
          un(j) = 0._wp
       else
          qn(j) = 0.5_wp
          un(j) = (3._wp/(x(j,n)-x(j,n-1))) * &
                  (ypn(j)-(y(j,n)-y(j,n-1))/(x(j,n)-x(j,n-1)))
       endif
    end do

    y2(:,n) = (un(:)-qn(:)*u(:,n-1)) / (qn(:)*y2(:,n-1)+1._wp)
!NEC$ novector
    do k = n-1,1,-1
       y2(:,k) = y2(:,k)*y2(:,k+1) + u(:,k)
    end do
  end subroutine init_splinex_vec
!------------------------------------------------------------------------------
  subroutine splint_vec (xa, ya, y2a, x, y, idim, kdimin, kdimout)
    !--------------------------------------------------------------------
    ! Spline interpolation (adapted from Numerical Recipes, sec.3.3):
    ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate
    ! a function (with the xai's in order), and given the array y2a(1:n),
    ! which is the output from init_splinex above, and given a value
    ! of x, this routine returns a cubic-spline interpolated value y.
    !
    ! Vectorizing version.
    !
    ! Caveat: abscissa values must be monotonically increasing!
    !--------------------------------------------------------------------
    integer,  intent(in)                           :: idim, kdimout, kdimin
    real(wp), intent(in),  dimension(idim,kdimin)  :: xa, ya, y2a
    real(wp), intent(in),  dimension(idim,kdimout) :: x
    real(wp), intent(out), dimension(idim,kdimout) :: y          ! Interp.value
    !----------------
    ! Local variables
    !----------------
    integer  :: k, ierr, i, k_, kk, kkh
    integer  :: khi (idim,kdimout), klo (idim,kdimout)
    real(wp) :: a, b, h

    ierr = 0
    !--------------------------------------------------------------------
    ! We will find the right place in the table by means of bisection.
    ! This is the only option for a fully vectorizing search.
    !--------------------------------------------------------------------
    klo = 1
    khi = kdimin
    kkh = int(log (real(kdimin))/log(2.)) + 1
!$omp parallel private(k,kk,i,k_,h,a,b) shared(ierr)
!$omp do schedule(static)
    do k = 1,kdimout
       do kk = 1, kkh
          do i = 1, idim
             if (khi(i,k)-klo(i,k) > 1) then
                k_ = (khi(i,k)+klo(i,k))/2
                if (xa(i,k_) > x(i,k)) then
                   khi(i,k)=k_
                else
                   klo(i,k)=k_
                endif
             endif
          enddo
       end do
    end do
!$omp end do nowait
    !----------------------------------------------
    ! klo and khi now bracket the input value of x.
    !----------------------------------------------
!$omp do schedule(static)
    do k = 1,kdimout
       do i   = 1, idim
          h = xa(i,khi(i,k))-xa(i,klo(i,k)) ! The xa's must be distinct.
          if (h == 0) then
            ierr = 1    ! Set error flag
            h    = 1    ! Dummy value to aid loop optimization
          endif
!         else
            !----------------------------------
            ! Evaluate cubic spline polynomial.
            !----------------------------------
            a = (xa(i,khi(i,k))-x(i,k))/h
            b = 1._wp - a       ! == (x(i,k)-xa(i,klo(i,k)))/h
            y(i,k) = a*ya(i,klo(i,k)) + b*ya(i,khi(i,k)) &
                   + (a*(a**2-1)*y2a(i,klo(i,k))         &
                   +  b*(b**2-1)*y2a(i,khi(i,k))) * h**2/6
!         endif
       enddo
    enddo
!$omp end do
!$omp end parallel
    if ( ierr == 1 ) call finish ("splint_vec","bad xa input")
  end subroutine splint_vec
!==============================================================================
  subroutine splint (xa, ya, y2a, x, y, sequential)
    !--------------------------------------------------------------------
    ! Spline interpolation (adapted from Numerical Recipes, sec.3.3):
    ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate
    ! a function (with the xai's in order), and given the array y2a(1:n),
    ! which is the output from init_splinex above, and given a value
    ! of x, this routine returns a cubic-spline interpolated value y.
    !
    ! For sequential calls, the subroutine hunt may be used to speed up
    ! the table search for subsequent invocations of splint.
    ! (This option is recommended only for reasonably large arrays
    ! where hunt can beat the simple binary search (size > 100?)).
    !--------------------------------------------------------------------
    real(wp), intent(in)           :: xa(:), ya(:), y2a(:)
    real(wp), intent(in)           :: x             ! Abscissa
    real(wp), intent(out)          :: y             ! Interpol. function value
    logical,  intent(in), optional :: sequential    ! Assume sequential call
    !----------------
    ! Local variables
    !----------------
    integer  :: n, k, khi, klo
    integer  :: klo_ = 1                ! Saved position from last search
!$omp threadprivate(klo_)
    real(wp) :: a, b, h
    logical  :: seq

    seq = .false.
    if (present (sequential)) seq = sequential
    n   = size (xa)
    if (seq .and. klo_ < n) then
       !--------------------------------------------------------------------
       ! For long arrays xa, reusing the information from the previous call
       ! may speed up locating the proper interval in the table.
       !--------------------------------------------------------------------
       klo = klo_
       call hunt (xa, x, klo)
       khi = klo + 1
    else
       !--------------------------------------------------------------------
       ! We will find the right place in the table by means of bisection.
       ! This is optimal if sequential calls to this routine are at random
       ! values of x.  If sequential calls are in order, and closely spaced,
       ! one would do better to store previous values of klo and khi and
       ! test if they remain appropriate on the next call.  (See above)
       !--------------------------------------------------------------------
       klo = 1
       khi = n
       do while (khi-klo > 1)
          k = (khi+klo)/2
          if (xa(k) > x) then
             khi=k
          else
             klo=k
          endif
       end do
    end if
    if (seq) klo_ = klo
    !----------------------------------------------
    ! klo and khi now bracket the input value of x.
    !----------------------------------------------
    h = xa(khi)-xa(klo)         ! The xa's must be distinct.
    if (h == 0) then
       write (0,*) khi, xa(khi), klo, xa(klo)
       call finish ("splint","bad xa input")
    end if
    !----------------------------------
    ! Evaluate cubic spline polynomial.
    !----------------------------------
    a = (xa(khi)-x)/h
    b = 1._wp - a       ! == (x-xa(klo))/h
    y = a*ya(klo) + b*ya(khi) &
      + (a*(a**2-1)*y2a(klo) + b*(b**2-1)*y2a(khi)) * h**2/6
  end subroutine splint

!==============================================================================
  subroutine init_spline2 (f, x, d2, df)
  !----------------------------------
  ! derives second derivatives (ln p)
  !----------------------------------
  real(wp) ,intent(in)           :: x  (:,:,:,:)
  real(wp) ,intent(in)           :: f  (:,:,:,:)
  real(wp) ,intent(out)          :: d2 (:,:,:,:)
  real(wp) ,intent(in) ,optional :: df (:,:,  :)  ! Surface derivative (ln p)

    integer :: d
!   integer :: i, j
    integer :: k

    k = size(f,3) - size(x ,3) + 1

    if (size(f,3) <size(x ,3)) call finish ('init_spline2','dims(f) <dims(x)')
    if (size(x,3)/=size(d2,3)) call finish ('init_spline2','dims(x)=/dims(d2)')

    do     d = lbound(f,4), ubound(f,4)
#if 1
       ! Call vectorized version of init_splinex:
       ! (similar to better performance on almost all platforms)
       if (present (df)) then
         call init_splinex_vec (x(:,:,:,d), f(:,:,k:,d), d2(:,:,:,d), &
                                size(x,1)*size(x,2), size(x,3), ypn=df(:,:,d))
       else
         call init_splinex_vec (x(:,:,:,d), f(:,:,k:,d), d2(:,:,:,d), &
                                size(x,1)*size(x,2), size(x,3))
       endif
#else
       ! Call old (scalar) version of init_splinex:
      do   j = lbound(f,2), ubound(f,2)
       do  i = lbound(f,1), ubound(f,1)
        if (present (df)) then
         call init_splinex(x(i,j,:,d), f(i,j,k:,d), d2(i,j,:,d), ypn=df(i,j,d))
        else
         call init_splinex(x(i,j,:,d), f(i,j,k:,d), d2(i,j,:,d))
        endif
       end do
      end do
#endif
    end do

  end subroutine init_spline2
!------------------------------------------------------------------------------
  subroutine vintp_spline (x, f, d2, x_int, f_int, df, a, lgeo, lcubic)
  real(wp),intent(in)          :: x    (:,:,:,:)
  real(wp),intent(in)          :: f    (:,:,:,:)
  real(wp),intent(in)          :: d2   (:,:,:,:)
  real(wp),intent(in)          :: x_int(:,:,:,:)
  real(wp),intent(out)         :: f_int(:,:,:,:)
  real(wp),intent(in),optional :: df   (:,:  ,:) ! Surface derivative (ln p)
  real(wp),intent(in),optional :: a    (:,:  ,:) ! R*gamma/g for geopotential
  logical ,intent(in),optional :: lgeo
  logical ,intent(in),optional :: lcubic         ! Cubic interpol., not spline
  !--------------------------------
  ! vertical (spline) interpolation
  !--------------------------------
    integer :: i,j,k,d, ke
    logical :: lspline, seq
#ifdef __NEC__
    integer :: is, ie, js ,je
    integer :: idim, kdimout, kdimin
#else
    integer :: lb(4), ub(4)             ! Work around Crayftn 8.3.3 OpenMP bug
#endif

    lspline = .true.; if (present (lcubic)) lspline = .not. lcubic

    ke      = ubound(x,3)
#ifdef __NEC__
    kdimin  = size (x    ,3)
    kdimout = size (x_int,3)
    idim    = size (x_int,1) * size (x_int,2)
    if (idim /= size(x,1)*size(x,2)) &
      call finish('vintp_spline','size (x_int) /= size (x)')

    do     d   = lbound(f    ,4), ubound(f    ,4)

          seq = .false.
               if (lspline) then
                  !---------------------
                  ! Spline interpolation
                  !---------------------
                  call splint_vec       &
                    (x    ( :,:,:,d ),  & ! <-- Argument grid
                     f    ( :,:,:,d ),  & ! <-- Gridded function
                     d2   ( :,:,:,d ),  & ! <-- 2nd derivative of spline
                     x_int( :,:,:,d ),  & ! <-- Interpolation point
                     f_int( :,:,:,d ),  & ! --> Interpolated function value
                     idim,              &
                     kdimin,            &
                     kdimout            )
               else
                  !--------------------
                  ! Cubic interpolation
                  !--------------------
!                 stop "cubic_intp_vec not yet implemented"
                  call finish('vintp_spline',&
                              "cubic_intp_vec not yet implemented")
!                 call cubic_intp_vec &
!                   (x    (i,j,:,d),  & ! <-- Argument grid
!                    f    (i,j,:,d),  & ! <-- Gridded function
!                    x_int(i,j,k,d),  & ! <-- Interpolation point
!                    f_int(i,j,k,d)   & ! --> Interpolated function value
!                   )
               end if
               !---------------------------------------------
               ! Extrapolation using lower boundary condition
               !---------------------------------------------
      if (present(df)) then
      is = lbound(f    ,1)
      ie = ubound(f    ,1)
      js = lbound(f    ,2)
      je = ubound(f    ,2)
            if (lgeo) then
      do k = lbound(x_int,3), ubound(x_int,3)
        do j  = js, je
          do i  = is, ie
                  !-------------------------------------------------------
                  ! Extrapolate using a uniform lapse rate atmosphere with
                  !   z(p) = z(p_0) + dz/dln(p)_{p=p_0}*[((p/p_0)^a-1)/a]
                  !-------------------------------------------------------
               if (x_int(i,j,k,d) > x (i,j,ke,d)) then
                  f_int (i,j,k,d) = f(i,j,ke,d)                           &
                     + (exp(a(i,j,d)*(x_int(i,j,k,d)-x(i,j,ke,d)))-1._wp) &
                     * df (i,j,d) / a (i,j,d)
               endif
          end do
        end do
      end do
            else
      do k = lbound(x_int,3), ubound(x_int,3)
        do j = js, je
          do i = is, ie
               if (x_int(i,j,k,d) > x (i,j,ke,d)) then
                  f_int (i,j,k,d) = f(i,j,ke,d) &
                                  + (x_int(i,j,k,d) - x(i,j,ke,d)) * df (i,j,d)
               endif
          end do
        end do
      end do
            endif
      endif
    end do
#else
! Work around OpenMP bug in Crayftn 8.3.3
    lb    = lbound (f)
    ub    = ubound (f)
    lb(3) = lbound (x_int,3)
    ub(3) = ubound (x_int,3)
!$omp parallel do private(d,i,j,k,seq) collapse(3) schedule(static)
    do     d   = lb(4), ub(4)
      do   j   = lb(2), ub(2)
        do i   = lb(1), ub(1)
          seq = .false.
          do k = lb(3), ub(3)
            if (x_int(i,j,k,d) <= x (i,j,ke,d)) then
               if (lspline) then
                  !---------------------
                  ! Spline interpolation
                  !---------------------
                  call splint         &
                    (x    (i,j,:,d),  & ! <-- Argument grid
                     f    (i,j,:,d),  & ! <-- Gridded function
                     d2   (i,j,:,d),  & ! <-- 2nd derivative of spline
                     x_int(i,j,k,d),  & ! <-- Interpolation point
                     f_int(i,j,k,d)   & ! --> Interpolated function value
                   !,sequential=seq   & !     Subsequent call in a sequence
                                        !     (no speed gain for nz<100)
                    )
!                 seq = .true.
               else
                  !--------------------
                  ! Cubic interpolation
                  !--------------------
                  call cubic_intp     &
                    (x    (i,j,:,d),  & ! <-- Argument grid
                     f    (i,j,:,d),  & ! <-- Gridded function
                     x_int(i,j,k,d),  & ! <-- Interpolation point
                     f_int(i,j,k,d)   & ! --> Interpolated function value
                    )
               end if
            else
              if (present(df)) then
               !---------------------------------------------
               ! Extrapolation using lower boundary condition
               !---------------------------------------------
               if (lgeo) then
                  !-------------------------------------------------------
                  ! Extrapolate using a uniform lapse rate atmosphere with
                  !   z(p) = z(p_0) + dz/dln(p)_{p=p_0}*[((p/p_0)^a-1)/a]
                  !-------------------------------------------------------
                  f_int (i,j,k,d) = f(i,j,ke,d)                           &
                     + (exp(a(i,j,d)*(x_int(i,j,k,d)-x(i,j,ke,d)))-1._wp) &
                     * df (i,j,d) / a (i,j,d)
               else
                  f_int (i,j,k,d) = f(i,j,ke,d) &
                                  + (x_int(i,j,k,d) - x(i,j,ke,d)) * df (i,j,d)
               endif
              else
                call finish('vintp_spline','df is not present')
              endif
            endif
          end do
        end do
      end do
    end do
!$omp end parallel do
#endif

  end subroutine vintp_spline
!==============================================================================
  subroutine cubic_intp (x, f, x_int, f_int)
    !---------------------------------------------------------
    ! Cubic interpolation using the (4-point) Lagrange formula
    !---------------------------------------------------------
    real(wp), intent(in)  :: x(:), f(:), x_int
    real(wp), intent(out) :: f_int
    !----------------
    ! Local variables
    !----------------
    integer :: n, jlo
    n = size (x)
    if (n < 4) then
       call finish ("cubic_intp","cannot handle less than 4 points")
    end if
    !--------------------------------------
    ! Locate interval that x_int belongs to
    !--------------------------------------
    jlo = n/2
    call hunt (x, x_int, jlo)
    if (jlo <= 0 .or. jlo >= n) then
       write (0,*) "jlo, x_int, x(1), x(N) =", jlo, x_int, x(1), x(n)
       call finish ("cubic_intp", "abscissa not in table")
    else if (jlo == 1) then
       ! Leftmost interval
!      call lagrange4 (x(1:4), f(1:4), x_int, f_int)
       call lagrange3 (x(1:3), f(1:3), x_int, f_int)      ! Quadratic interpol.
    else if (jlo == n-1) then
       ! Rightmost interval
!      call lagrange4 (x(n-3:n), f(n-3:n), x_int, f_int)
       call lagrange3 (x(n-2:n), f(n-2:n), x_int, f_int)  ! Quadratic interpol.
    else
       ! Interval [jlo:jlo+1]
       call lagrange4 (x(jlo-1:jlo+2), f(jlo-1:jlo+2), x_int, f_int)
    end if
  contains
    subroutine lagrange4 (xa, ya, x, y)
      !---------------------------------------------------------
      ! Cubic interpolation using the (4-point) Lagrange formula
      !---------------------------------------------------------
      real(wp), intent(in)  :: xa(4), ya(4), x
      real(wp), intent(out) :: y
      !----------------
      ! Local variables
      !----------------
      real(wp) :: w(4), u(4)

      u(1:4) = (x - xa(1:4))
      w(1) = u(2)*u(3)*u(4) / ((xa(1)-xa(2)) * (xa(1)-xa(3)) * (xa(1)-xa(4)))
      w(2) = u(1)*u(3)*u(4) / ((xa(2)-xa(1)) * (xa(2)-xa(3)) * (xa(2)-xa(4)))
      w(3) = u(1)*u(2)*u(4) / ((xa(3)-xa(1)) * (xa(3)-xa(2)) * (xa(3)-xa(4)))
      w(4) = u(1)*u(2)*u(3) / ((xa(4)-xa(1)) * (xa(4)-xa(2)) * (xa(4)-xa(3)))
      y = dot_product (w, ya)
    end subroutine lagrange4
    !--
    subroutine lagrange3 (xa, ya, x, y)
      !-----------------------------------------------------------
      ! Quadratic interpolation using the 3-point Lagrange formula
      !-----------------------------------------------------------
      real(wp), intent(in)  :: xa(3), ya(3), x
      real(wp), intent(out) :: y
      !----------------
      ! Local variables
      !----------------
      real(wp) :: w(3), u(3)

      u(1:3) = (x - xa(1:3))
      w(1) = u(2)*u(3) / ((xa(1)-xa(2)) * (xa(1)-xa(3)))
      w(2) = u(1)*u(3) / ((xa(2)-xa(1)) * (xa(2)-xa(3)))
      w(3) = u(1)*u(2) / ((xa(3)-xa(1)) * (xa(3)-xa(2)))
      y = dot_product (w, ya)
    end subroutine lagrange3
  end subroutine cubic_intp
!==============================================================================
  subroutine fit_polynomial (x, y, n, stat, wgts, yp, p, rms)
  !---------------------------------
  ! fit polinominal to some function
  !---------------------------------
    real(wp), intent(in)                    :: x(:)    ! x values
    real(wp), intent(in)                    :: y(:)    ! y values
    integer,  intent(in)                    :: n       ! order of polynomial
    integer,  intent(out)                   :: stat    ! exit status
    real(wp), intent(in),  optional, target :: wgts(:) ! y values of fitted polynomial
    real(wp), intent(out), optional         :: yp(:)   ! y values of fitted polynomial
    real(wp), intent(out), optional         :: p(:)    ! fitted parameters of polynomial
    real(wp), intent(out), optional         :: rms     ! rms of fit

    integer               :: npar, i, j, k, m
    integer,  allocatable :: ipiv(:)
    real(wp), allocatable :: mat(:)
    real(wp), allocatable :: b(:)
    real(wp), allocatable :: yf(:)
    real(wp), pointer     :: w(:) => null()
    real(wp)              :: wsum

#define ERR_EXIT(text) call error_local(text) ; return

    stat = 0
    npar = n + 1
    npar = min(npar, size(x))

    if (size(x) .ne. size(y)) THEN
      ERR_EXIT('sizes of x and y are inconsistent')
    end if
    if (present(wgts)) then
      if (size(x) .ne. size(wgts)) THEN
        ERR_EXIT('sizes of x and wgts are inconsistent')
      end if
      w => wgts
    else
      allocate(w(size(x)))
      w = 1._wp
    end if
    m = size(x)

    ! Allocate arrays
    allocate(mat(npar*(npar+1)/2), b(npar), ipiv(npar), stat=stat)
    IF (stat /= 0) THEN
      ERR_EXIT('allocation failed')
    END IF

    ! Calculate sum(x_i^j * x_i^k) and sum(y_i * x_i^k)
!     forall (k=1:npar)
!       forall (j=1:k)
!         mat(j+(k-1)*k/2) = dot_product(x(:)**(k-1), x(:)**(j-1))
!       end forall
!       b(k) = dot_product(x(:)**(k-1), y(:))
!     end forall
    mat = 0._wp
    b   = 0._wp
    do i = 1, m
      do k=1,npar
        do j=1,k
          mat(j+(k-1)*k/2) = mat(j+(k-1)*k/2) + x(i)**(k-1) * x(i)**(j-1) * w(i)
        end do
        b(k) = b(k) + x(i)**(k-1) * y(i) * w(i)
      end do
    end do

    call DSPTRF('U', npar, mat, ipiv, stat)
    if (stat /= 0) then
      ERR_EXIT('DSPTRF failed')
    end if
    call DSPTRS('U', npar, 1, mat, ipiv, b, npar, stat)
    if (stat /= 0) then
      ERR_EXIT('DSPTRS failed')
    end if

    if (present(p)) then
      if (size(p) < npar) then
        ERR_EXIT('size of array p too small')
      end if
      p(1:npar) = b(:)
    end if
    if (present(yp) .or. present(rms)) then
      allocate(yf(size(x)))
      yf(:) = b(1)
      do k = 1, n
        yf(:) = yf(:) + b(k+1) * x(:)**k
      end do
      if (present(yp)) then
        j = min(size(yp), size(yf))
        yp(1:j) = yf(1:j)
      end if
      if (present(rms)) then
        rms  = 0._wp
        wsum = 0._wp
        do i = 1, m
          rms  = rms  + w(i) * (y(i) - yf(i))**2
          wsum = wsum + w(i)
        end do
        if (wsum > 0._wp) rms = rms / wsum
        if (rms  > 0._wp) rms = sqrt(rms)
      end if
    end if

    if (.not.present(wgts)) deallocate(w)

  contains

    subroutine error_local(text)
      character(len=*), intent(in) :: text
      write(0,*) '*** ERROR in fit_polynomial: '//trim(text)
      stat=13
      if (.not.present(wgts) .and. associated(w)) deallocate(w)
    end subroutine error_local

  end subroutine fit_polynomial
!==============================================================================
  subroutine polint (xa, ya, x, y, dy)
    !---------------------------------------------------------------
    ! Polynomial interpolation using Neville's algorithm.
    ! Given arrays xa(1:N) and ya(1:N), and given a value x,
    ! this routine returns a value y, and an error estimate dy.
    ! If P(x) is the polynomial of degree N-1 such that
    ! P(xa_i)=ya_i, i=1,...,N, then the returned value y=P(x).
    ! (Taken from Numerical Recipes, sec.3.1, and adapted to F95)
    !---------------------------------------------------------------
    real(wp), intent(in)  :: xa(:), ya(:), x
    real(wp), intent(out) :: y, dy
    !----------------
    ! Local variables
    !----------------
    integer  :: n, i, m, ns
    real(wp) :: den, dif, dift, ho, hp, w
    real(wp), dimension(size (ya)) :: c, d
    n = size (ya)
    ns=1
    dif=abs (x-xa(1))
    do i=1,n            ! Find the index ns of the closest table entry,
       dift=abs (x-xa(i))
       if (dift < dif) then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)       ! and initialize the tableau of c's and d's.
       d(i)=ya(i)
    enddo
    y=ya(ns)            ! Initial approximation to y.
    ns=ns-1
    do m=1,n-1          ! For each column of the tableau,
       do i=1,n-m       ! loop over the current c's and d's and update them.
          ho=xa(i)  -x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if (den == 0) then
             ! This error can occur only if two input xa's are
             ! (to within roundoff) identical.
             write(0,*) 'failure in polint'
             stop
          end if
          den=w/den
          d(i)=hp*den   ! Here the c's and d's are updated.
          c(i)=ho*den
       enddo
       ! After each column in the tableau is completed, we decide which
       ! correction, c or d, we want to add to our accumulating value
       ! of y, i.e., which path to take through the tableau--forking
       ! up or down. We do this in such a way as to take the most
       ! "straight line" route through the tableau to its apex, updating
       ! ns accordingly to keep track of where we are. This route keeps
       ! the partial approximations centered (insofar as possible) on
       ! the target x. The last dy added is thus the error indication.
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    enddo
  end subroutine polint
!==============================================================================
  function erf (x)
  !----------------------------------
  ! Returns the error function ERF(X)
  !----------------------------------
  real(wp) ,intent(in) :: x
  real(wp)             :: erf
    if (x<0.0_wp) then
      erf = -gammp (0.5_wp,x**2)
    else
      erf =  gammp (0.5_wp,x**2)
    endif
  end function erf
!------------------------------------------------------------------------------
  function gammp (a, x)
  real(wp) ,intent(in) :: a
  real(wp) ,intent(in) :: x
  real(wp)             :: gammp
  !----------------------------------------------
  ! Returns the incomplete gamma function P(a,x).
  !----------------------------------------------
    real(wp) :: gammcf ,gamser ,gln
    if ((x < 0.0_wp) .or. (a <= 0.0_wp)) call finish ('gammp','bad arguments')
    if (x < a + 1.0_wp) then
      call gser (gamser,a,x,gln) ! Use the series representation
      gammp = gamser
    else
      call gcf (gammcf,a,x,gln)  ! Use the continued fraction representation
      gammp = 1._wp - gammcf     ! and take its complement.
    endif
  end function gammp
!------------------------------------------------------------------------------
  subroutine gser (gamser,a,x,gln)
  real(wp) ,intent(out) :: gamser
  real(wp) ,intent(in)  :: a
  real(wp) ,intent(in)  :: x
  real(wp) ,intent(out) :: gln
  !---------------------------------------------------------------------
  ! Returns the incomplete gamma function P(a,x) evaluated by its series
  ! representation as GAMSER. Also returns ln(Gamma(a)) as GLN.
  !---------------------------------------------------------------------
    integer  ,parameter :: itmax = 100
    real(wp) ,parameter :: eps   = 3.e-7_wp
    integer             :: n
    real(wp)            :: ap, del, sum
    gln = gammln(a)
    if (x <= 0.0_wp) then
      if (x < 0.0_wp) call finish ('gser','x < 0')
      gamser = 0.0_wp
      return
    endif
    ap  = a
    sum = 1._wp / a
    del = sum
    do n=1,itmax
      ap  = ap + 1._wp
      del = del * x / ap
      sum = sum + del
      if (abs(del) < abs(sum)*eps) goto 10
    enddo
    write (0,*)
    write (0,*) 'gser:  a, x, itmax =',a,x,itmax
    write (0,*) 'gser:  A too large, ITMAX too small'
    call finish('gser','A too large, ITMAX too small')
10  gamser = sum * exp (-x + a * log(x) - gln)
  end subroutine gser
!------------------------------------------------------------------------------
  subroutine gcf(gammcf,a,x,gln)
  real(wp) ,intent(out) :: gammcf
  real(wp) ,intent(in)  :: a
  real(wp) ,intent(in)  :: x
  real(wp) ,intent(out) :: gln
  !-----------------------------------------------------------------------
  ! Returns the incomplete gamma function Q(a,x) evaluated by its
  ! continued fraction representation as GAMMCF. Also returns ln(Gamma(a))
  ! as GLN.
  ! Parameters: ITMAX is the maximum allowed number of iterations. EPS is
  ! the relative accuracy; FPMIN is a number near the smallest
  ! representable floating-point number.
  !-----------------------------------------------------------------------
    integer  ,parameter :: itmax = 100
    real(wp) ,parameter :: eps   = 3.e-7_wp
    real(wp) ,parameter :: fpmin = 1.e-30_wp
    integer             :: i
    real(wp)            :: an,b,c,d,del,h
    gln = gammln (a)
    b   = x + 1.0_wp - a
    c   = 1.0_wp / fpmin
    d   = 1.0_wp / b
    h   = d
    do i=1,itmax
      an  = -i * (i-a)
      b   = b + 2.0_wp
      d   = an * d + b
      if (abs(d) < fpmin) d = fpmin
      c   = b + an / c
      if (abs(c) < fpmin) c = fpmin
      d   = 1.0_wp / d
      del = d * c
      h   = h * del
      if (abs(del-1.0_wp) < eps) goto 20
    enddo
    write (0,*)
    write (0,*)  'gcf:  a, x, itmax =',a,x,itmax
    write (0,*)  'gcf:  A too large, ITMAX too small'
    call finish ('gcf','A too large, ITMAX too small')
20  gammcf = exp (-x + a * log(x) - gln) * h
  end subroutine gcf
!------------------------------------------------------------------------------
  elemental function gammln (xx)
  real(wp) ,intent(in) :: xx
  real(wp)             :: gammln
  !--------------------------------------------
  ! Returns the value ln(Gamma(xx)) for xx > 0.
  !--------------------------------------------
    !--------------------------------------------------------------------
    ! Internal arithmetic will be done in double precision, a nicety that
    ! you can omit if five-figure accuracy is good enough.
    !--------------------------------------------------------------------
    integer             :: j
    real(dp)            :: ser, tmp, x, y
    real(dp) ,parameter :: stp = 2.5066282746310005_dp
    real(dp) ,parameter :: cof(6) =                                           &
      (/ 76.18009172947146_dp,  -86.50532032941677_dp,  24.01409824083091_dp, &
        -1.231739572450155_dp,.1208650973866179e-2_dp, -.5395239384953e-5_dp /)
    x   = xx
    y   = x
    tmp = x + 5.5_dp
    tmp = (x + 0.5_dp) * log (tmp) - tmp
    ser = 1.000000000190015_dp
    do j=1,6
      y   = y + 1._dp
      ser = ser + cof(j) / y
    enddo
    gammln = tmp + log (stp * ser / x)
  end function gammln
!==============================================================================
  DOUBLE PRECISION ELEMENTAL FUNCTION DERF (X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for erf(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, January 8, 1985
!
!--------------------------------------------------------------------
    DOUBLE PRECISION, INTENT(IN)  :: X
!------------------------------------------------------------------
    DOUBLE PRECISION RESULT
!------------------------------------------------------------------
    CALL CALERF (X, RESULT, JINT=0)
    DERF = RESULT
  END FUNCTION DERF
!--------------------------------------------------------------------
  DOUBLE PRECISION ELEMENTAL FUNCTION DERFC (X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for erfc(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, January 8, 1985
!
!--------------------------------------------------------------------
    DOUBLE PRECISION, INTENT(IN)  :: X
!------------------------------------------------------------------
    DOUBLE PRECISION RESULT
!------------------------------------------------------------------
    CALL CALERF (X, RESULT, JINT=1)
    DERFC = RESULT
  END FUNCTION DERFC
!------------------------------------------------------------------
  DOUBLE PRECISION ELEMENTAL FUNCTION DERFCX (X)
!------------------------------------------------------------------
!
! This subprogram computes approximate values for exp(x*x) * erfc(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, March 30, 1987
!
!------------------------------------------------------------------
    DOUBLE PRECISION, INTENT(IN)  :: X
!------------------------------------------------------------------
    DOUBLE PRECISION RESULT
!------------------------------------------------------------------
    CALL CALERF (X, RESULT, JINT=2)
    DERFCX = RESULT
  END FUNCTION DERFCX
!=======================================================================
  ELEMENTAL SUBROUTINE CALERF (ARG, RESULT, JINT)
!------------------------------------------------------------------
!
! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!   for a real argument  x.  It contains three FUNCTION type
!   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
!   and one SUBROUTINE type subprogram, CALERF.  The calling
!   statements for the primary entries are:
!
!                   Y=ERF(X)     (or   Y=DERF(X)),
!
!                   Y=ERFC(X)    (or   Y=DERFC(X)),
!   and
!                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
!
!   The routine  CALERF  is intended for internal packet use only,
!   all computations within the packet being concentrated in this
!   routine.  The function subprograms invoke  CALERF  with the
!   statement
!
!          CALL CALERF(ARG,RESULT,JINT)
!
!   where the parameter usage is as follows
!
!      Function                     Parameters for CALERF
!       call              ARG                  Result          JINT
!
!     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
!     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
!     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
!
!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XMIN   = the smallest positive floating-point number.
!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!            the negative of the solution to the equation
!            2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!            A conservative value is the largest machine number X
!            such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!            W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!            machine precision.  A conservative value is
!            1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!            of XINF and 1/[sqrt(pi)*XMIN].
!
!   Approximate values for some important machines are:
!
!                          XMIN       XINF        XNEG     XSMALL
!
!  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
!
!
!                          XBIG       XHUGE       XMAX
!
!  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ERFC = 0      for  ARG .GE. XBIG;
!
!                       ERFCX = XINF  for  ARG .LT. XNEG;
!      and
!                       ERFCX = 0     for  ARG .GE. XMAX.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, EXP
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: March 19, 1990
!
!------------------------------------------------------------------
    DOUBLE PRECISION, INTENT(IN)  :: ARG
    INTEGER,          INTENT(IN)  :: JINT
    DOUBLE PRECISION, INTENT(OUT) :: RESULT
!------------------------------------------------------------------
    INTEGER I
    DOUBLE PRECISION :: DEL,X,XDEN,XNUM,Y,YSQ
!------------------------------------------------------------------
!  Mathematical constants
!------------------------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: ZERO   =  0.0D0
    DOUBLE PRECISION, PARAMETER :: HALF   =  0.5D0
    DOUBLE PRECISION, PARAMETER :: ONE    =  1.0D0
    DOUBLE PRECISION, PARAMETER :: TWO    =  2.0D0
    DOUBLE PRECISION, PARAMETER :: FOUR   =  4.0D0
    DOUBLE PRECISION, PARAMETER :: SIXTEN = 16.0D0
    DOUBLE PRECISION, PARAMETER :: SQRPI  = 5.6418958354775628695D-1
    DOUBLE PRECISION, PARAMETER :: THRESH = 0.46875D0
!------------------------------------------------------------------
!  Machine-dependent constants
!------------------------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: &
         XINF = 1.79D308, XNEG = -26.628D0, XSMALL = 1.11D-16, &
         XBIG = 26.543D0, XHUGE = 6.71D7,   XMAX = 2.53D307
!------------------------------------------------------------------
!  Coefficients for approximation to  erf  in first interval
!------------------------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: &
         A(5)=(/3.16112374387056560D00,1.13864154151050156D02, &
                3.77485237685302021D02,3.20937758913846947D03, &
                1.85777706184603153D-1/)
    DOUBLE PRECISION, PARAMETER :: &
         B(4)=(/2.36012909523441209D01,2.44024637934444173D02, &
                1.28261652607737228D03,2.84423683343917062D03/)
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in second interval
!------------------------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: &
         C(9)=(/5.64188496988670089D-1,8.88314979438837594D0,  &
                6.61191906371416295D01,2.98635138197400131D02, &
                8.81952221241769090D02,1.71204761263407058D03, &
                2.05107837782607147D03,1.23033935479799725D03, &
                2.15311535474403846D-8/)
    DOUBLE PRECISION, PARAMETER :: &
         D(8)=(/1.57449261107098347D01,1.17693950891312499D02, &
                5.37181101862009858D02,1.62138957456669019D03, &
                3.29079923573345963D03,4.36261909014324716D03, &
                3.43936767414372164D03,1.23033935480374942D03/)
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in third interval
!------------------------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: &
         P(6)=(/3.05326634961232344D-1,3.60344899949804439D-1, &
                1.25781726111229246D-1,1.60837851487422766D-2, &
                6.58749161529837803D-4,1.63153871373020978D-2/)
    DOUBLE PRECISION, PARAMETER :: &
         Q(5)=(/2.56852019228982242D00,1.87295284992346047D00, &
                5.27905102951428412D-1,6.05183413124413191D-2, &
                2.33520497626869185D-3/)
!------------------------------------------------------------------
    X = ARG
    Y = ABS(X)
    IF (Y .LE. THRESH) THEN
!------------------------------------------------------------------
!  Evaluate  erf  for  |X| <= 0.46875
!------------------------------------------------------------------
       YSQ = ZERO
       IF (Y .GT. XSMALL) YSQ = Y * Y
       XNUM = A(5)*YSQ
       XDEN = YSQ
       DO I = 1, 3
          XNUM = (XNUM + A(I)) * YSQ
          XDEN = (XDEN + B(I)) * YSQ
       END DO
       RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
       IF (JINT .NE. 0) RESULT = ONE - RESULT
       IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
       GO TO 800
!------------------------------------------------------------------
!  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
!------------------------------------------------------------------
    ELSE IF (Y .LE. FOUR) THEN
       XNUM = C(9)*Y
       XDEN = Y
       DO I = 1, 7
          XNUM = (XNUM + C(I)) * Y
          XDEN = (XDEN + D(I)) * Y
       END DO
       RESULT = (XNUM + C(8)) / (XDEN + D(8))
       IF (JINT .NE. 2) THEN
          YSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-YSQ)*(Y+YSQ)
          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
       END IF
!------------------------------------------------------------------
!  Evaluate  erfc  for |X| > 4.0
!------------------------------------------------------------------
    ELSE
       RESULT = ZERO
       IF (Y .GE. XBIG) THEN
          IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 300
          IF (Y .GE. XHUGE) THEN
             RESULT = SQRPI / Y
             GO TO 300
          END IF
       END IF
       YSQ = ONE / (Y * Y)
       XNUM = P(6)*YSQ
       XDEN = YSQ
       DO I = 1, 4
          XNUM = (XNUM + P(I)) * YSQ
          XDEN = (XDEN + Q(I)) * YSQ
       END DO
       RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
       RESULT = (SQRPI -  RESULT) / Y
       IF (JINT .NE. 2) THEN
          YSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-YSQ)*(Y+YSQ)
          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
       END IF
    END IF
!------------------------------------------------------------------
!  Fix up for negative argument, erf, etc.
!------------------------------------------------------------------
300 IF (JINT .EQ. 0) THEN
       RESULT = (HALF - RESULT) + HALF
       IF (X .LT. ZERO) RESULT = -RESULT
    ELSE IF (JINT .EQ. 1) THEN
       IF (X .LT. ZERO) RESULT = TWO - RESULT
    ELSE
       IF (X .LT. ZERO) THEN
          IF (X .LT. XNEG) THEN
             RESULT = XINF
          ELSE
             YSQ = AINT(X*SIXTEN)/SIXTEN
             DEL = (X-YSQ)*(X+YSQ)
             Y = EXP(YSQ*YSQ) * EXP(DEL)
             RESULT = (Y+Y) - RESULT
          END IF
       END IF
    END IF
800 RETURN
!---------- Last card of CALERF ----------
  END SUBROUTINE CALERF
!------------------------------------------------------------------------------
! Numerical integration of unequally spaced data
! A cubic interpolation polynomial with 4 supporting points is used to
! estimate the integral:
!
! P. E. Gill and G. F. Miller
! An algorithm for the integration of unequally spaced data
! The Computer Journal, Vol. 15, No. 1, p. 80-83, 1972
! http://comjnl.oxfordjournals.org/content/15/1/80.abstract
!
! The integral is always computed over the whole interval covered by the
! data. To obtain integrals over subintervals [ia, ib] call the function
! with the corresponding subinterval:
! integ = integpolycube (fx(ia:ib))
!
! Warning:
! The interpolation polynomial does work well only if the spacing of
! the data is is not too inhomogeneous. Single intervals with much larger
! spacings (> 100 times) can lead to very wrong results.
!
! Input
! fx(i,j) - i=1, n, index of nodes, n - number of array elements
!           j=1 - x, supporting points
!           j=2 - y, function y = f(x)
!
! Output
! integpolycube (return value) - integral
!
!---------------------------------------------------------------------
function integpolycube ( fx )

implicit none

! List of calling arguments:
real (wp)                              :: integpolycube
real (wp), dimension(:,:), intent(in)  :: fx

!  List of local variables:
real (wp) :: sum, a, b
real (wp) :: Dx0, Dx1, Dx2, DDx1, DDx2, DDx3
real (wp) :: Df0, Df1, Df2, DDf1, DDf2, DDf3
integer              :: i, j, k,  n
!---------------------------------------------------------------------

j = lbound(fx,1)
k = ubound(fx,1)
n = k - j + 1

if (n .lt. 4) then
   ! At least 4 points are required to compute the integral
   write(*,*) 'integpolycube> More than three points are required, n = ', n
   integpolycube = 0.0_wp
else

   sum = 0.0_wp
   do i=j, k-3

      if (i .eq. j) then
         ! Compute first set of differences and the integral over the
         ! first interval:

         ! Delta x ( x => fx(i,1) )
         Dx0 = fx(i+1,1) - fx(i,1)
         Dx1 = fx(i+2,1) - fx(i+1,1)
         Dx2 = fx(i+3,1) - fx(i+2,1)
         DDx1 = Dx0 + Dx1
         DDx2 = Dx1 + Dx2
         DDx3 = DDx2 + Dx0
         ! Delta y / Delta x  ( y => fx(i,2) )
         Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
         Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
         Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
         DDf1 = ( Df1 - Df0 ) / DDx1
         DDf2 = ( Df2 - Df1 ) / DDx2
         DDf3 = ( DDf2 - DDf1 ) / DDx3

         ! Integral
         a = (Dx0 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         b = DDf1 / 6.0_wp
         sum = Dx0 * ( fx(i,2) + 0.50_wp*Dx0*Df0 +    &
                       (a - b) * Dx0**2           )

      end if

      ! Compute integral using 4 point formula

      ! Delta x ( x => fx(i,1) )
      !Dx0 = fx(i+1,1) - fx(i,1)
      !Dx1 = fx(i+2,1) - fx(i+1,1)
      Dx2 = fx(i+3,1) - fx(i+2,1)
      !DDx1 = Dx0 + Dx1
      DDx2 = Dx1 + Dx2
      DDx3 = DDx2 + Dx0
      ! Delta y / Delta x  ( y => fx(i,2) )
      !Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
      !Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
      Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
      !DDf1 = ( Df1 - Df0 ) / DDx1
      DDf2 = ( Df2 - Df1 ) / DDx2
      DDf3 = ( DDf2 - DDf1 ) / DDx3

      ! Integral
      a = 0.50_wp * ( DDf1 + DDf2 )
      b = 0.50_wp * ( Dx0 - Dx2 ) * DDf3
      sum = sum + Dx1 * ( 0.50_wp * ( fx(i+1,2) + fx(i+2,2) ) -       &
                          Dx1**2 * (a + b) / 6.0_wp )

      if (i .eq. k-3) then
         ! Compute  integral of the last interval using
         ! the same differences as above:

         ! Integral
         a = (Dx2 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         b = DDf2 / 6.0_wp
         sum = sum + Dx2 * ( fx(n,2) - 0.50_wp*Dx2*Df2 -    &
                             (a + b) * Dx2**2           )

         exit

      end if

      ! Reuse differences from last loop:
      Dx0  = Dx1
      Dx1  = Dx2
      DDx1 = DDx2
      Df0  = Df1
      Df1  = Df2
      DDf1 = DDf2

   end do

   integpolycube = sum

end if

end function integpolycube
!------------------------------------------------------------------------------
! Tangent-linear of integpolycube
! Numerical integration of unequally spaced data
! A cubic interpolation polynomial with 4 supporting points is used to
! estimate the integral:
!
! P. E. Gill and G. F. Miller
! An algorithm for the integration of unequally spaced data
! The Computer Journal, Vol. 15, No. 1, p. 80-83, 1972
! http://comjnl.oxfordjournals.org/content/15/1/80.abstract
!
! The integral is always computed over the whole interval covered by the
! data. To obtain integrals over subintervals [ia, ib] call the function
! with the corresponding subinterval:
! integ = integpolycube (fx(ia:ib))
!
! Input
! fx(i,j) - i=1, n, index of nodes, n - number of array elements
!           j=1 - x, supporting points
!           j=2 - y, function y = f(x)
!
! fx_tl   - increments of fx
!
! Output
! integpolycube_tl (return value) - derivative of integral
!
!---------------------------------------------------------------------
function integpolycube_tl ( fx, fx_tl )

implicit none

! List of calling arguments:
real (wp)                              :: integpolycube_tl
real (wp), dimension (:,:), intent(in) :: fx, fx_tl

!  List of local variables:
real (wp) :: sum, a, b
real (wp) :: Dx0, Dx1, Dx2, DDx1, DDx2, DDx3
real (wp) :: Df0, Df1, Df2, DDf1, DDf2, DDf3
real (wp) :: tsum, ta, tb
real (wp) :: tDx0, tDx1, tDx2, tDDx1, tDDx2, tDDx3
real (wp) :: tDf0, tDf1, tDf2, tDDf1, tDDf2, tDDf3
integer   :: i, j, k,  n
!---------------------------------------------------------------------

j = lbound(fx,1)
k = ubound(fx,1)
n = k - j + 1

if (n .lt. 4) then
   ! At least 4 points are required to compute the integral
   write(*,*) 'integpolycube> More than three points are required, n = ', n
   integpolycube_tl = 0.0_wp
else

   sum = 0.0_wp
   do i=j, k-3

      if (i .eq. j) then
         ! Compute first set of differences and the integral over the
         ! first interval:

         ! Delta x ( x => fx(i,1) )
         Dx0 = fx(i+1,1) - fx(i,1)
         tDx0 =  fx_tl(i+1,1) - fx_tl(i,1)

         Dx1 = fx(i+2,1) - fx(i+1,1)
         tDx1 = fx_tl(i+2,1) - fx_tl(i+1,1)

         Dx2 = fx(i+3,1) - fx(i+2,1)
         tDx2 = fx_tl(i+3,1) - fx_tl(i+2,1)

         DDx1 = Dx0 + Dx1
         tDDx1 = tDx0 + tDx1

         DDx2 = Dx1 + Dx2
         tDDx2 = tDx1 + tDx2

         DDx3 = DDx2 + Dx0
         tDDx3 = tDDx2 + tDx0

         ! Delta y / Delta x  ( y => fx(i,2) )
         Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
         tDf0 = (fx_tl(i+1,2) - fx_tl(i,2))/Dx0 - &
                tDx0 * ( fx(i+1,2) - fx(i,2)   ) / Dx0**2

         Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
         tDf1 = ( fx_tl(i+2,2) - fx_tl(i+1,2) ) / Dx1 -    &
                tDx1 * ( fx(i+2,2) - fx(i+1,2) ) / Dx1**2

         Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
         tDf2 = ( fx_tl(i+3,2) - fx_tl(i+2,2) ) / Dx2 - &
                tDx2 * ( fx(i+3,2) - fx(i+2,2) ) / Dx2**2

         DDf1 = ( Df1 - Df0 ) / DDx1
         tDDf1 = ( tDf1 - tDf0 ) / DDx1 -          &
                  tDDx1 * ( Df1 - Df0 ) / DDx1**2

         DDf2 = ( Df2 - Df1 ) / DDx2
         tDDf2 = ( tDf2 - tDf1 ) / DDx2 -          &
                 tDDx2 * ( Df2 - Df1 ) / DDx2**2

         DDf3 = ( DDf2 - DDf1 ) / DDx3
         tDDf3 = ( tDDf2 - tDDf1 ) / DDx3 -           &
                 tDDx3 * ( DDf2 - DDf1 ) / DDx3**2

         ! Integral
         a = (Dx0 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         ta = (tDx0 * DDf3 + 2.0_wp*tDx1*DDf3 +    &
              (Dx0 + 2.0_wp*Dx1)*tDDf3) / 12.0_wp

         b = DDf1 / 6.0_wp
         tb = tDDf1 / 6.0_wp

         sum = Dx0 * ( fx(i,2) + 0.50_wp*Dx0*Df0 +    &
                       (a - b) * Dx0**2           )
         tsum = tDx0 * (fx(i,2) + Dx0*Df0 + 3.0_wp * (a - b) * Dx0**2) +  &
                fx_tl(i,2)*Dx0 + 0.50_wp*tDf0*Dx0**2 + (ta-tb)*Dx0**3

      end if

      ! Compute integral using 4 point formula

      ! Delta x ( x => fx(i,1) )
      !Dx0 = fx(i+1,1) - fx(i,1)
      !Dx1 = fx(i+2,1) - fx(i+1,1)
      Dx2 = fx(i+3,1) - fx(i+2,1)
      tDx2 = fx_tl(i+3,1) - fx_tl(i+2,1)

      !DDx1 = Dx0 + Dx1
      DDx2 = Dx1 + Dx2
      tDDx2 = tDx1 + tDx2

      DDx3 = DDx2 + Dx0
      tDDx3 = tDDx2 + tDx0

      ! Delta y / Delta x  ( y => fx(i,2) )
      !Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
      !Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
      Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
      tDf2 = ( fx_tl(i+3,2) - fx_tl(i+2,2) ) / Dx2 - &
               tDx2 * ( fx(i+3,2) - fx(i+2,2) ) / Dx2**2

      !DDf1 = ( Df1 - Df0 ) / DDx1
      DDf2 = ( Df2 - Df1 ) / DDx2
      tDDf2 = ( tDf2 - tDf1 ) / DDx2 -          &
                tDDx2 * ( Df2 - Df1 ) / DDx2**2

      DDf3 = ( DDf2 - DDf1 ) / DDx3
      tDDf3 = ( tDDf2 - tDDf1 ) / DDx3 -            &
                tDDx3 * ( DDf2 - DDf1 ) / DDx3**2

      ! Integral
      a = 0.50_wp * ( DDf1 + DDf2 )
      ta = 0.50_wp * ( tDDf1 + tDDf2 )

      b = 0.50_wp * ( Dx0 - Dx2 ) * DDf3
      tb = 0.50_wp * ( ( tDx0 - tDx2 ) * DDf3 +               &
                       tDDf3 * ( Dx0 - Dx2 )   )

      sum = sum + Dx1 * ( 0.50_wp * ( fx(i+1,2) + fx(i+2,2) ) -       &
                          Dx1**2 * (a + b) / 6.0_wp )
      tsum = tsum + tDx1 * ( 0.50_wp * (( fx(i+1,2) + fx(i+2,2) ) -       &
                             Dx1**2 * (a + b)) ) +                       &
                    0.50_wp*Dx1*(fx_tl(i+1,2)+fx_tl(i+2,2))       -       &
                    (tb + ta) * Dx1**3 / 6.0_wp

      if (i .eq. k-3) then
         ! Compute  integral of the last interval using
         ! the same differences as above:

         ! Integral
         a = (Dx2 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         ta = ( (tDx2 + 2.0_wp*tDx1) * DDf3 +                &
                (Dx2 + 2.0_wp*Dx1) * tDDf3   )  / 12.0_wp
         b = DDf2 / 6.0_wp
         tb = tDDf2 / 6.0_wp
         sum = sum + Dx2 * ( fx(n,2) - 0.50_wp*Dx2*Df2 -    &
                             (a + b) * Dx2**2           )
         tsum = tsum + Dx2*fx_tl(n,2) -  0.50_wp*tDf2*Dx2**2 -  &
                (ta + tb) * Dx2**3 +                           &
                ( fx(n,2) - Dx2*Df2 - 3.0_wp*(a + b) * Dx2**2 ) * tDx2

         exit

      end if

      ! Reuse differences from last loop:
      Dx0  = Dx1
      tDx0  = tDx1

      Dx1  = Dx2
      tDx1  = tDx2

      DDx1 = DDx2
      tDDx1 = tDDx2

      Df0  = Df1
      tDf0  = tDf1

      Df1  = Df2
      tDf1  = tDf2

      DDf1 = DDf2
      tDDf1 = tDDf2

   end do

   integpolycube_tl = tsum

   if (integpolycube_tl /= integpolycube_tl) then
      write(*,*) 'integpolycube_tl> integpolycube_tl = NaN', char(10),   &
           'max fx_tl = ', maxval(fx_tl)
   end if

end if

end function integpolycube_tl
!------------------------------------------------------------------------------
! Adjoint of integpolycube
! Numerical integration of unequally spaced data
! A cubic interpolation polynomial with 4 supporting points is used to
! estimate the integral:
!
! P. E. Gill and G. F. Miller
! An algorithm for the integration of unequally spaced data
! The Computer Journal, Vol. 15, No. 1, p. 80-83, 1972
! http://comjnl.oxfordjournals.org/content/15/1/80.abstract
!
! The integral is always computed over the whole interval covered by the
! data. To obtain integrals over subintervals [ia, ib] call the function
! with the corresponding subinterval:
! integ = integpolycube (fx(ia:ib))
!
! Input
! fx(i,j) - i=1, n, index of nodes, n - number of array elements
!           j=1 - x, supporting points
!           j=2 - y, function y = f(x)
! integpolycubeAdj - variation of the integral
!
! Output
! integpolycube_ad (return value) - derivatives
!
!---------------------------------------------------------------------
subroutine integpolycube_ad ( fx, fx_ad, integpolycubeadj )

implicit none

! List of calling arguments:
real (wp), intent(inout)                :: integpolycubeadj
real (wp), dimension (:,:), intent(in)  :: fx
real (wp), dimension (:,:), intent(out) :: fx_ad

!  List of local variables:
real (wp) :: a, b
real (wp) :: Dx0, Dx1, Dx2, DDx1, DDx2, DDx3
real (wp) :: Df0, Df1, Df2, DDf1, DDf2, DDf3
real (wp) :: asum, aa, ab
real (wp) :: aDx0, aDx1, aDx2, aDDx1, aDDx2, aDDx3
real (wp) :: aDf0, aDf1, aDf2, aDDf1, aDDf2, aDDf3

integer              :: i, j, k,  n
!---------------------------------------------------------------------

j = lbound(fx,1)
k = ubound(fx,1)
n = k - j + 1

if (n .lt. 4) then
   ! At least 4 points are required to compute the integral
   write(*,*) 'integpolycube> More than three points are required, n = ', n
   fx_ad = 0.0_wp
else

   ! Adjoint code
   fx_ad = 0.0_wp
   ! Delta x ( x => fx(i,1) )
   aDx0 = 0.0_wp
   aDx1 = 0.0_wp
   aDx2 = 0.0_wp
   aDDx1 = 0.0_wp
   aDDx2 = 0.0_wp
   aDDx3 = 0.0_wp
   ! Delta y / Delta x  ( y => fx(i,2) )
   aDf0 = 0.0_wp
   aDf1 = 0.0_wp
   aDf2 = 0.0_wp
   aDDf1 = 0.0_wp
   aDDf2 = 0.0_wp
   aDDf3 = 0.0_wp

   ! integpolycube = sum
   asum = integpolycubeAdj

   do i=k-3, j, -1

      ! Delta x ( x => fx(i,1) )
      Dx0 = fx(i+1,1) - fx(i,1)
      Dx1 = fx(i+2,1) - fx(i+1,1)
      Dx2 = fx(i+3,1) - fx(i+2,1)
      DDx1 = Dx0 + Dx1
      DDx2 = Dx1 + Dx2
      DDx3 = DDx2 + Dx0
      ! Delta y / Delta x  ( y => fx(i,2) )
      Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
      Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
      Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
      DDf1 = ( Df1 - Df0 ) / DDx1
      DDf2 = ( Df2 - Df1 ) / DDx2
      DDf3 = ( DDf2 - DDf1 ) / DDx3

      if (i .eq. k-3) then
         ! Compute  integral of the last interval using
         ! the same differences as above:

         ! Integral
         a = (Dx2 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         b = DDf2 / 6.0_wp

         ! sum = sum + Dx2 * ( fx(n,2) - 0.50D0*Dx2*Df2 -    &
         !                    (a + b) * Dx2**2           )
         aDx2 = aDx2 + ( fx(n,2)-Dx2*Df2-3.0_wp * (a+b)*Dx2**2 ) * asum
         fx_ad(n,2) = fx_ad(n,2) + Dx2 * asum
         aDf2 = aDf2 - asum * 0.50_wp * Dx2**2
         aa = -asum * Dx2**3
         ab = -asum * Dx2**3

         ! b = DDf2 / 6.0D0
         aDDf2 = aDDf2 + ab / 6.0_wp

         ! a = (Dx2 + 2.0D0*Dx1) * DDf3 / 12.0D0
         aDx2 =  aDx2 + (DDf3 / 12.0_wp) * aa
         aDx1 = aDx1 + (DDf3 / 6.0_wp) * aa
         aDDf3 = aDDf3 + ((Dx2 + 2.0_wp*Dx1)/12.0_wp) * aa

      end if

      ! Compute integral using 4 point formula
      a = 0.50_wp * ( DDf1 + DDf2 )
      b = 0.50_wp * ( Dx0 - Dx2 ) * DDf3

      ! Integral
      ! sum = sum + Dx1 * ( 0.50D0 * ( fx(i+1,2) + fx(i+2,2) ) -       &
      !                     Dx1**2 * (a + b) / 6.0D0 )
      aDx1 = aDx1 + 0.50_wp*( (fx(i+1,2)+fx(i+2,2)) - Dx1**2 * (a + b) ) * asum
      fx_ad(i+1,2) = fx_ad(i+1,2) + 0.50_wp*Dx1 * asum
      fx_ad(i+2,2) = fx_ad(i+2,2) + 0.50_wp*Dx1 * asum
      aa = -asum * Dx1**3  / 6.0_wp
      ab = -asum * Dx1**3  / 6.0_wp

      ! b = 0.50D0 * ( Dx0 - Dx2 ) * DDf3
      aDx0 = aDx0 + 0.50_wp * DDf3 * ab
      aDx2 = aDx2 - 0.50_wp * DDf3 * ab
      aDDf3 = aDDf3 + 0.50_wp * ( Dx0 - Dx2 ) * ab

      ! a = 0.50D0 * ( DDf1 + DDf2 )
      aDDf1 = aDDf1 + 0.50_wp * aa
      aDDf2 = aDDf2 + 0.50_wp * aa

      ! DDf3 = ( DDf2 - DDf1 ) / DDx3
      aDDf2 = aDDf2 + aDDf3 / DDx3
      aDDf1 = aDDf1 - aDDf3 / DDx3
      aDDx3 = aDDx3 - aDDf3 * ( DDf2 - DDf1 ) / DDx3**2

      ! DDf2 = ( Df2 - Df1 ) / DDx2
      aDf2 = aDf2 + aDDf2 / DDx2
      aDf1 = aDf1 - aDDf2 / DDx2
      aDDx2 = aDDx2 - aDDf2*( Df2 - Df1 ) / DDx2**2

      ! DDf1 = ( Df1 - Df0 ) / DDx1
      aDf1 = aDf1 + aDDf1 / DDx1
      aDf0 = aDf0 - aDDf1 / DDx1
      aDDx1 = aDDx1 - aDDf1*( Df1 - Df0 ) / DDx1**2

      ! Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
      fx_ad(i+3,2) = fx_ad(i+3,2) + aDf2 / Dx2
      fx_ad(i+2,2) = fx_ad(i+2,2) - aDf2 / Dx2
      aDx2 = aDx2 - aDf2*( fx(i+3,2) - fx(i+2,2) ) / Dx2**2

      ! Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
      fx_ad(i+2,2) = fx_ad(i+2,2) + aDf1 / Dx1
      fx_ad(i+1,2) = fx_ad(i+1,2) - aDf1 / Dx1
      aDx1 = aDx1 - aDf1*( fx(i+2,2) - fx(i+1,2) ) / Dx1**2

      ! Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
      fx_ad(i+1,2) = fx_ad(i+1,2) + aDf0 / Dx0
      fx_ad(i,2) = fx_ad(i,2) - aDf0 / Dx0
      aDx0 = aDx0 - aDf0*( fx(i+1,2) - fx(i,2)   ) / Dx0**2

      ! DDx3 = DDx2 + Dx0
      aDDx2 = aDDx2 + aDDx3
      aDx0 = aDx0 + aDDx3

      ! DDx2 = Dx1 + Dx2
      aDx1 = aDx1 + aDDx2
      aDx2 = aDx2  + aDDx2

      ! DDx1 = Dx0 + Dx1
      aDx0 = aDx0 + aDDx1
      aDx1 = aDx1 + aDDx1

      ! Dx2 = fx(i+3,1) - fx(i+2,1)
      fx_ad(i+3,1) = fx_ad(i+3,1) + aDx2
      fx_ad(i+2,1) = fx_ad(i+2,1) - aDx2

      ! Dx1 = fx(i+2,1) - fx(i+1,1)
      fx_ad(i+2,1) = fx_ad(i+2,1) + aDx1
      fx_ad(i+1,1) = fx_ad(i+1,1) - aDx1

      ! Dx0 = fx(i+1,1) - fx(i,1)
      fx_ad(i+1,1) = fx_ad(i+1,1) + aDx0
      fx_ad(i,1) = fx_ad(i,1) - aDx0

      if (i .eq. j) then
         ! Compute first set of differences and the integral over the
         ! first interval:

         ! Delta x ( x => fx(i,1) )
         aDx0 = 0.0_wp
         aDx1 = 0.0_wp
         aDx2 = 0.0_wp
         aDDx1 = 0.0_wp
         aDDx2 = 0.0_wp
         aDDx3 = 0.0_wp
         ! Delta y / Delta x  ( y => fx(i,2) )
         aDf0 = 0.0_wp
         aDf1 = 0.0_wp
         aDf2 = 0.0_wp
         aDDf1 = 0.0_wp
         aDDf2 = 0.0_wp
         aDDf3 = 0.0_wp

         a = (Dx0 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         b = DDf1 / 6.0_wp

         ! sum = Dx0 * ( fx(i,2) + 0.50D0*Dx0*Df0 +    &
         !               (a - b) * Dx0**2           )
         aDx0 = aDx0 + asum * ( fx(i,2) + Dx0*Df0 +       &
                                3.0D0*(a - b) * Dx0**2 )
         fx_ad(i,2) = fx_ad(i,2) + Dx0 * asum
         aDf0 = aDf0 + 0.50_wp * asum * Dx0**2
         aa = asum * Dx0**3
         ab = -asum * Dx0**3

         ! b = DDf1 / 6.0D0
         aDDf1 = aDDf1 + ab / 6.0_wp

         ! a = (Dx0 + 2.0D0*Dx1) * DDf3 / 12.0D0
         aDx0 = aDx0 + aa * DDf3 / 12.0_wp
         aDx1 = aDx1 + aa * DDf3 / 6.0_wp
         aDDf3 = aDDf3 + aa * (Dx0 + 2.0_wp*Dx1) / 12.0_wp

         ! DDf3 = ( DDf2 - DDf1 ) / DDx3
         aDDf2 = aDDf2 + aDDf3 / DDx3
         aDDf1 = aDDf1 - aDDf3 / DDx3
         aDDx3 = aDDx3 - aDDf3 * ( DDf2 - DDf1 ) / DDx3**2

         ! DDf2 = ( Df2 - Df1 ) / DDx2
         aDf2 = aDf2 + aDDf2 / DDx2
         aDf1 = aDf1 - aDDf2 / DDx2
         aDDx2 = aDDx2 - aDDf2*( Df2 - Df1 ) / DDx2**2

         ! DDf1 = ( Df1 - Df0 ) / DDx1
         aDf1 = aDf1 + aDDf1 / DDx1
         aDf0 = aDf0 - aDDf1 / DDx1
         aDDx1 = aDDx1 - aDDf1*( Df1 - Df0 ) / DDx1**2

         ! Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
         fx_ad(i+3,2) = fx_ad(i+3,2) + aDf2 / Dx2
         fx_ad(i+2,2) = fx_ad(i+2,2) - aDf2 / Dx2
         aDx2 = aDx2 - aDf2*( fx(i+3,2) - fx(i+2,2) ) / Dx2**2

         ! Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
         fx_ad(i+2,2) = fx_ad(i+2,2) + aDf1 / Dx1
         fx_ad(i+1,2) = fx_ad(i+1,2) - aDf1 / Dx1
         aDx1 = aDx1 - aDf1*( fx(i+2,2) - fx(i+1,2) ) / Dx1**2

         ! Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
         fx_ad(i+1,2) = fx_ad(i+1,2) + aDf0 / Dx0
         fx_ad(i,2) = fx_ad(i,2) - aDf0 / Dx0
         aDx0 = aDx0 - aDf0*( fx(i+1,2) - fx(i,2)   ) / Dx0**2

         ! DDx3 = DDx2 + Dx0
         aDDx2 = aDDx2 + aDDx3
         aDx0 = aDx0 + aDDx3

         ! DDx2 = Dx1 + Dx2
         aDx1 = aDx1 + aDDx2
         aDx2 = aDx2  + aDDx2

         ! DDx1 = Dx0 + Dx1
         aDx0 = aDx0 + aDDx1
         aDx1 = aDx1 + aDDx1

         ! Dx2 = fx(i+3,1) - fx(i+2,1)
         fx_ad(i+3,1) = fx_ad(i+3,1) + aDx2
         fx_ad(i+2,1) = fx_ad(i+2,1) - aDx2

         ! Dx1 = fx(i+2,1) - fx(i+1,1)
         fx_ad(i+2,1) = fx_ad(i+2,1) + aDx1
         fx_ad(i+1,1) = fx_ad(i+1,1) - aDx1

         ! Dx0 = fx(i+1,1) - fx(i,1)
         fx_ad(i+1,1) = fx_ad(i+1,1) + aDx0
         fx_ad(i,1) = fx_ad(i,1) - aDx0

      end if

      ! Delta x ( x => fx(i,1) )
      aDx0 = 0.0_wp
      aDx1 = 0.0_wp
      aDx2 = 0.0_wp
      aDDx1 = 0.0_wp
      aDDx2 = 0.0_wp
      aDDx3 = 0.0_wp
      ! Delta y / Delta x  ( y => fx(i,2) )
      aDf0 = 0.0_wp
      aDf1 = 0.0_wp
      aDf2 = 0.0_wp
      aDDf1 = 0.0_wp
      aDDf2 = 0.0_wp
      aDDf3 = 0.0_wp

   end do

end if

integpolycubeadj = 0.0_wp

end subroutine integpolycube_ad

!==============================================================================

subroutine simplex (x, dx, y, func, ftol, itmax, iter)
!------------------------------------------------------
! Multidimensional minimization of the function FUNC(X)
!------------------------------------------------------
real(wp) ,intent(inout) :: x (:) ! set of parameters
real(wp) ,intent(in)    :: dx(:) ! initioal variation of set of parameters
real(wp) ,intent(out)   :: y     ! minimum value
real(wp) ,intent(in)    :: ftol  ! tolerance
integer  ,intent(in)    :: itmax ! max. number of iterations allowed
integer  ,intent(out)   :: iter  ! number of iterations used
interface
  function func (x)
    use mo_kind, only: wp
    real(wp) ,intent(in) :: x (:)
    real(wp)             :: func
  end function func
end interface

  real(wp) :: z (size(x)+1)
  real(wp) :: p (size(x)+1,size(x))
  integer  :: i, n

  n = size(x)

  do i = 1, n
    p(i,:) = x
    p(i,i) = x(i) + dx(i)
    z(i)   = func (p(i,:))
  end do
  i = n + 1
  p(i,:) = x
  z(i)   = func (p(i,:))

  call amoeba (p, z, func, ftol, itmax, iter)

  i = sum (minloc (z))
  x = p(i,:)
  y = z(i)

end subroutine simplex

!------------------------------------------------------------------------------

subroutine amoeba (p, y, func, ftol, itmax, iter)
!-------------------------------------------------------------------
! (adapted from Numerical Recipes)
! Multidimensional minimization of the function FUNC(X) where X is
! an NDIM-dimensional vector, by the downhill simplex method of
! Nelder and Mead. Input is a matrix P whose NDIM+1 rows are NDIM-
! dimensional vectors which are the vertices of the starting simplex
! (Logical dimensions of P are P(NDIM+1,NDIM). Also input is the
! vector Y of length NDIM +1, whose components must be pre-initialized
! to the values of FUNC evaluated at the NDIM+1 vertices (rows) of P;
! and FTOL the fractional convergence tolerance to be achieved in the
! function value. On output, P and Y will have been reset to NDIM+1
! new points all within FTOL of a minimum function value, and ITER
! gives the number of iterations taken.
!--------------------------------------------------------------------
real(wp) ,intent(inout) :: p(:,:) ! (ndim+1,ndim) vertices of the simplex
real(wp) ,intent(inout) :: y(:)   ! (ndim+1)      values of func at vertices
real(wp) ,intent(in)    :: ftol   ! fractional convergence tolerance
integer  ,intent(in)    :: itmax  ! maximum allowed number of iterations
integer  ,intent(out)   :: iter   ! number of iterations taken
interface
  function func (x)
    use mo_kind, only: wp
    real(wp) ,intent(in) :: x (:)
    real(wp)             :: func
  end function func
end interface

  !--------------------------------------------------------------------
  ! parameters alpha,beta,gamma define the expansions and contractions.
  !--------------------------------------------------------------------
  real(wp) ,parameter :: alpha = 1.0_wp
  real(wp) ,parameter :: beta  = 0.5_wp
  real(wp) ,parameter :: gamma = 2.0_wp

  integer  :: i, ihi, ilo, inhi, j, mpts, ndim
  real(wp) :: ypr, rtol, yprr
  real(wp) :: pr   (size(p,2))
  real(wp) :: prr  (size(p,2))
  real(wp) :: pbar (size(p,2))

  ndim = size (p,2)
  mpts = ndim + 1
  iter = 0
  do
    ilo = 1
    if (y(1) > y(2)) then
      ihi  = 1
      inhi = 2
    else
      ihi  = 2
      inhi = 1
    endif
    do i=1, mpts
      if (y(i) < y(ilo)) ilo = i
      if (y(i) > y(ihi)) then
        inhi = ihi
        ihi  = i
      else if (y(i) > y(inhi)) then
        if (i /= ihi) inhi = i
      endif
    enddo
    !------------------------------------------------------------------
    ! Compute the fractional range from highest to lowest and return if
    ! satisfactory.
    !------------------------------------------------------------------
    rtol = 2._wp * abs (y(ihi) - y(ilo)) / (abs (y(ihi)) + abs (y(ilo)))
    if (rtol <  ftol)  return
    if (iter == itmax) call finish ('amoeba','exceeding maximum iterations.')
    iter    = iter + 1
    pbar(:) = 0._wp
    do i = 1, mpts
      if(i /= ihi) pbar(:) = pbar(:) + p(i,:)
    enddo
    pbar(:) = pbar(:) / ndim
    pr  (:) = (1._wp + alpha) * pbar(:) - alpha * p(ihi,:)
    ypr = func (pr)
    if (ypr <= y(ilo)) then
      prr(:) = gamma * pr(:) + (1._wp - gamma) * pbar(:)
      yprr   = func(prr)
      if (yprr < y(ilo)) then
        p(ihi,:) = prr(:)
        y(ihi)   = yprr
      else
        p(ihi,:) = pr(:)
        y(ihi)   = ypr
      endif
    else if (ypr >= y(inhi)) then
      if (ypr < y(ihi)) then
        p(ihi,:) = pr(:)
        y(ihi)   = ypr
      endif
      prr(:) = beta * p(ihi,:) + (1._wp - beta) * pbar(:)
      yprr   = func(prr)
      if (yprr < y(ihi)) then
        p(ihi,:) = prr(:)
        y(ihi)   = yprr
      else
        do i = 1, mpts
          if (i /= ilo) then
            pr(:)   = 0.5_wp * (p(i,:) + p(ilo,:))
            p (i,:) = pr(:)
            y(i) = func(pr)
          endif
        enddo
      endif
    else
      p(ihi,:) = pr(:)
      y(ihi)   = ypr
    endif
  enddo
end subroutine amoeba
!==============================================================================
! Binary search in a sorted list - integer
!
! Input:
! a     - sorted array
! value - value to be found in the array
!
! Return value:
! binsearchi - index of array element that contains the value
!            - 0 if value not in array
!
! The function is based on the bisection algorithm described in the
! Numerical Recipes. The value should be found in about
! log_2 n steps.
!---------------------------------------------------------------------------
function binsearchi (a, value)
  integer                     :: binsearchi
  integer, intent(in), target :: a(:)
  integer, intent(in)         :: value
  integer, pointer            :: p(:)
  integer                     :: mid, offset

  p => a
  binsearchi = 0
  offset = 0
  do while (size(p) > 0)
     mid = size(p)/2 + 1
     if (p(mid) > value) then
        p => p(:mid-1)
     else if (p(mid) < value) then
        offset = offset + mid
        p => p(mid+1:)
     else
        ! value has been found
        binsearchi = offset + mid
        return
     end if
  end do

end function binsearchi
!---------------------------------------------------------------------------
! Binary search in a sorted list - integer
!
! Input:
! a     - sorted array
! value - value to be found in the array
!
! Return value:
! binsearchi - index of array element that contains the value
!            - 0 if value not in array
!
! The function is based on the bisection algorithm described in the
! Numerical Recipes. The value should be found in about
! log_2 n steps.
!---------------------------------------------------------------------------
function binsearchi8 (a, value)
  integer                          :: binsearchi8
  integer (i8), intent(in), target :: a(:)
  integer (i8), intent(in)         :: value
  integer (i8), pointer            :: p(:)
  integer                          :: mid, offset

  p => a
  binsearchi8 = 0
  offset = 0
  do while (size(p) > 0)
     mid = size(p)/2 + 1
     if (p(mid) > value) then
        p => p(:mid-1)
     else if (p(mid) < value) then
        offset = offset + mid
        p => p(mid+1:)
     else
        ! value has been found
        binsearchi8 = offset + mid
        return
     end if
  end do

end function binsearchi8
!---------------------------------------------------------------------------
! Binary search in a sorted list - real
!
! Input:
! a     - sorted array
! value - value to be found in the array
!
! Return value:
! binsearchi - index of array element that contains the value
!            - 0 if value not in array
!
! The function is based on the bisection algorithm described in the
! Numerical Recipes. The value should be found in about
! log_2 n steps.
!---------------------------------------------------------------------------
function binsearchwp (a, value)
  integer                       :: binsearchwp
  real (wp), intent(in), target :: a(:)
  real (wp), intent(in)         :: value
  real (wp), pointer            :: p(:)
  integer                       :: mid, offset

  p => a
  binsearchwp = 0
  offset = 0
  do while (size(p) > 0)
     mid = size(p)/2 + 1
     if (p(mid) > value) then
        p => p(:mid-1)
     else if (p(mid) < value) then
        offset = offset + mid
        p => p(mid+1:)
     else
        ! value has been found
        binsearchwp = offset + mid
        return
     end if
  end do

end function binsearchwp
!==============================================================================
! Sorting using the merge sort algorithm
! The array is replaced by the sorted array
!
! runsorti   - driver routine that allocates a temporary array
!              and calls the actual merge sort routine.
! mergesorti - The merge sort routine
! mergearri  - Merge two sorted arrays into a new sorted array
!
! runsorti(a) - Input : Unsorted array a
!               Output: Sorted array a
!---------------------------------------------------------------------------
subroutine runsorti(a)
  integer, intent(inout) :: a(:)
  integer, allocatable   :: a2(:)

  allocate( a2((size(a) + 1) / 2) )
  call mergesorti(a, a2)

end subroutine runsorti

recursive subroutine mergesorti(a, a2)
  integer, intent(inout) :: a(:)
  integer, intent(inout) :: a2(:)
  integer                :: half

  half = (size(a) + 1) / 2
  if (size(a) < 2) then
     continue
  else if (size(a) == 2) then
     if (a(1) > a(2)) then
        call swap(a(1), a(2))
     end if
  else
     call mergesorti(a( : half), a2)
     call mergesorti(a(half + 1 :), a2)
     if (a(half) > a(half + 1)) then
        a2(1 : half) = a(1 : half)
        call mergearri(a2(1 : half), a(half + 1:), a)
     endif
  end if

  contains

    subroutine swap(x, y)
      integer, intent(inout) :: x, y
      integer :: tmp
      tmp = x
      x = y
      y = tmp
    end subroutine swap

end subroutine mergesorti

subroutine mergearri(a, b, c)

  ! The targe attribute is necessary, because A .or. B might overlap with C.
  integer, target, intent(in) :: A(:), B(:)
  integer, target, intent(inout) :: C(:)
  integer :: i, j, k

  if (size(a) + size(b) > size(c)) then
     call finish ('mergearri','cannot merge, array c too small')
  end if

  i = 1; j = 1
  do k = 1, size(c)
     if (i <= size(a) .and. j <= size(b)) then
        if (a(i) <= b(j)) then
           c(k) = a(i)
           i = i + 1
        else
           c(k) = b(j)
           j = j + 1
        end if
     else if (i <= size(a)) then
        c(k) = a(i)
        i = i + 1
     else if (j <= size(b)) then
        c(k) = b(j)
        j = j + 1
     end if
  end do

end subroutine mergearri
!---------------------------------------------------------------------------
! Sorting using the merge sort algorithm
! The array is replaced by the sorted array
!
! runsorti8   - driver routine that allocates a temporary array
!               and calls the actual merge sort routine.
! mergesorti8 - The merge sort routine
! mergearri8  - Merge two sorted arrays into a new sorted array
!
! runsorti8(a) - Input : Unsorted array a
!                Output: Sorted array a
!---------------------------------------------------------------------------
subroutine runsorti8(a)
  integer (i8), intent(inout) :: a(:)
  integer (i8), allocatable   :: a2(:)

  allocate( a2((size(a) + 1) / 2) )
  call mergesorti8(a, a2)

end subroutine runsorti8

recursive subroutine mergesorti8(a, a2)
  integer (i8), intent(inout) :: a(:)
  integer (i8), intent(inout) :: a2(:)
  integer                     :: half

  half = (size(a) + 1) / 2
  if (size(a) < 2) then
     continue
  else if (size(a) == 2) then
     if (a(1) > a(2)) then
        call swap(a(1), a(2))
     end if
  else
     call mergesorti8(a( : half), a2)
     call mergesorti8(a(half + 1 :), a2)
     if (a(half) > a(half + 1)) then
        a2(1 : half) = a(1 : half)
        call mergearri8(a2(1 : half), a(half + 1:), a)
     endif
  end if

  contains

    subroutine swap(x, y)
      integer (i8), intent(inout) :: x, y
      integer (i8):: tmp
      tmp = x
      x = y
      y = tmp
    end subroutine swap

end subroutine mergesorti8

subroutine mergearri8(a, b, c)

  ! The targe attribute is necessary, because A .or. B might overlap with C.
  integer (i8), target, intent(in) :: A(:), B(:)
  integer (i8), target, intent(inout) :: C(:)
  integer :: i, j, k

  if (size(a) + size(b) > size(c)) then
     call finish ('mergearri','cannot merge, array c too small')
  end if

  i = 1; j = 1
  do k = 1, size(c)
     if (i <= size(a) .and. j <= size(b)) then
        if (a(i) <= b(j)) then
           c(k) = a(i)
           i = i + 1
        else
           c(k) = b(j)
           j = j + 1
        end if
     else if (i <= size(a)) then
        c(k) = a(i)
        i = i + 1
     else if (j <= size(b)) then
        c(k) = b(j)
        j = j + 1
     end if
  end do

end subroutine mergearri8


subroutine runsortc(a)
  character (len=*), intent(inout)    :: a(:)
  character (len=len(a)), allocatable :: a2(:)

  allocate( a2((size(a) + 1) / 2) )
  call mergesortc(a, a2)

end subroutine runsortc

recursive subroutine mergesortc(a, a2)
  character (len=*), intent(inout) :: a(:)
  character (len=*), intent(inout) :: a2(:)
  integer                     :: half

  half = (size(a) + 1) / 2
  if (size(a) < 2) then
     continue
  else if (size(a) == 2) then
     if (a(1) > a(2)) then
        call swap(a(1), a(2))
     end if
  else
     call mergesortc(a( : half), a2)
     call mergesortc(a(half + 1 :), a2)
     if (a(half) > a(half + 1)) then
        a2(1 : half) = a(1 : half)
        call mergearrc(a2(1 : half), a(half + 1:), a)
     endif
  end if

contains

    subroutine swap(x, y)
      character (len=*), intent(inout) :: x, y
      character (len=len(x)) :: tmp
      tmp = x
      x = y
      y = tmp
    end subroutine swap

end subroutine mergesortc

subroutine mergearrc(a, b, c)

  ! The targe attribute is necessary, because A .or. B might overlap with C.
  character (len=*), target, intent(in) :: A(:), B(:)
  character (len=*), target, intent(inout) :: C(:)
  integer :: i, j, k

  if (size(a) + size(b) > size(c)) then
     !call finish ('mergearri','cannot merge, array c too small')
  end if

  i = 1; j = 1
  do k = 1, size(c)
     if (i <= size(a) .and. j <= size(b)) then
        if (a(i) <= b(j)) then
           c(k) = a(i)
           i = i + 1
        else
           c(k) = b(j)
           j = j + 1
        end if
     else if (i <= size(a)) then
        c(k) = a(i)
        i = i + 1
     else if (j <= size(b)) then
        c(k) = b(j)
        j = j + 1
     end if
  end do

end subroutine mergearrc

!==============================================================================
end module mo_algorithms
