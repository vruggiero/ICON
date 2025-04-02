! Source module for the radar forward operator EMVORADO
!
! ---------------------------------------------------------------
! Copyright (C) 2005-2024, DWD, KIT
! Contact information: ulrich.blahak (at) dwd.de 
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE radar_gamma_functions_vec

!------------------------------------------------------------------------------
!
! Description: Utilitiy functions for the for the computation of the complete and incomplete gamma
!              functions. These are used in the radar forward operator EMVORADO.
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!
!------------------------------------------------------------------------------

  USE radar_kind, ONLY : dp

  ! Modules used:


  IMPLICIT NONE

  PUBLIC

  ! Type to hold the lookup table for the incomplete gamma functions.
  ! The table is divided into a low resolution part, which spans the
  ! whole range of x-values up to the 99.5 % x-value, and a high resolution part for the
  ! smallest 1 % of these x-values, where the incomplete gamma function may increase
  ! very rapidly and nonlinearily, depending on paramter a.
  ! For some applications (e.g., Newtons Method in subroutine
  ! graupel_hail_conv_wetgrowth_Dg_gamlook() ), this rapid change requires a much higher
  ! accuracy of the table lookup as compared to be achievable with the low resolution table.

  TYPE gamlookuptable
    ! Number of bins in the tables:
    INTEGER                              :: n        ! Internal number of bins (low res part)
    INTEGER                              :: nhr      ! Internal number of bins (high res part)
    REAL(KIND=dp)                        :: a        ! a-parameter
    REAL(KIND=dp), DIMENSION(:), POINTER :: x        ! vector of x-parameters (limit of integration) -
                                                     ! always starts at 0 and has equidistant dx (low resolution part)
    REAL(KIND=dp), DIMENSION(:), POINTER :: xhr      ! vector of x-parameters (limit of integration) -
                                                     ! always starts at 0 and has equidistant dxhr (high resolution part)
    REAL(KIND=dp)                        :: dx       ! dx   (low resolution part)
    REAL(KIND=dp)                        :: dxhr     ! dxhr (high resolution part)
    REAL(KIND=dp)                        :: odx      ! one over dx
    REAL(KIND=dp)                        :: odxhr    ! one over dxhr
    REAL(KIND=dp), DIMENSION(:), POINTER :: igf      ! value of the inc. gamma function at (a,x) (low res)
    REAL(KIND=dp), DIMENSION(:), POINTER :: igfhr    ! value of the inc. gamma function at (a,x) (high res)
  END TYPE gamlookuptable

  PRIVATE :: gamma_help_cf, gamma_help_ser

CONTAINS

  REAL(KIND=dp) FUNCTION gammln(x)
    
    !*******************************************************************************
    !
    ! Log(gamma function)
    !
    ! Gamma function approximation from Lanczos (1964)
    !
    ! Lanczos, C., 1964: A Precision Approximation of the Gamma Function,
    ! Journal of the Society for Industrial and Applied Mathematics Series B
    ! Numerical Analysis, vol. 1, pp. 86-96, doi: 10.1137/0701008.
    !
    ! Coefficients are for g=5 and N=6.
    !
    ! See: https://en.wikipedia.org/wiki/Lanczos_approximation
    !
    !*******************************************************************************

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(IN) :: x

    REAL(KIND=dp) :: tmp, p

    REAL(KIND=dp), PARAMETER :: c0 =  1.000000000190015_dp
    REAL(KIND=dp), PARAMETER :: c1 = 76.18009172947146_dp
    REAL(KIND=dp), PARAMETER :: c2 = -86.50532032941677_dp
    REAL(KIND=dp), PARAMETER :: c3 = 24.01409824083091_dp
    REAL(KIND=dp), PARAMETER :: c4 = -1.231739572450155_dp
    REAL(KIND=dp), PARAMETER :: c5 = 0.1208650973866179e-2_dp
    REAL(KIND=dp), PARAMETER :: c6 = -0.5395239384953e-5_dp
    REAL(KIND=dp), PARAMETER :: squareroot_two_pi = SQRT(2.0_dp*4.0_dp*atan(1.0_dp))

    tmp = x + 4.5_dp;
    p = squareroot_two_pi * (c0 + c1/x + c2/(x+1.0_dp) + c3/(x+2.0_dp) + c4/(x+3.0_dp) + c5/(x+4.0_dp) + c6/(x+5.0_dp))
    gammln = (x-0.5_dp) * LOG(tmp) - tmp + LOG(p)

  END FUNCTION gammln

  REAL(KIND=dp) FUNCTION gfct(x)
    
    !*******************************************************************************
    !
    ! Gamma function approximation from Lanczos (1964)
    !
    ! Lanczos, C., 1964: A Precision Approximation of the Gamma Function,
    ! Journal of the Society for Industrial and Applied Mathematics Series B
    ! Numerical Analysis, vol. 1, pp. 86-96, doi: 10.1137/0701008.
    !
    ! Coefficients are for g=5 and N=6.
    !
    ! See: https://en.wikipedia.org/wiki/Lanczos_approximation
    !
    !*******************************************************************************

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(IN) :: x

    REAL(KIND=dp) :: tmp, p

    REAL(KIND=dp), PARAMETER :: c0 =  1.000000000190015_dp
    REAL(KIND=dp), PARAMETER :: c1 = 76.18009172947146_dp
    REAL(KIND=dp), PARAMETER :: c2 = -86.50532032941677_dp
    REAL(KIND=dp), PARAMETER :: c3 = 24.01409824083091_dp
    REAL(KIND=dp), PARAMETER :: c4 = -1.231739572450155_dp
    REAL(KIND=dp), PARAMETER :: c5 = 0.1208650973866179e-2_dp
    REAL(KIND=dp), PARAMETER :: c6 = -0.5395239384953e-5_dp
    REAL(KIND=dp), PARAMETER :: squareroot_two_pi = SQRT(2.0_dp*4.0_dp*atan(1.0_dp))

    tmp = x + 4.5_dp;
    p = squareroot_two_pi * (c0 + c1/x + c2/(x+1.0_dp) + c3/(x+2.0_dp) + c4/(x+3.0_dp) + c5/(x+4.0_dp) + c6/(x+5.0_dp))
    gfct = p * EXP( (x-0.5_dp) * LOG(tmp) - tmp )

  END FUNCTION gfct

  FUNCTION gfct_frac_approx(a,b) RESULT(x)

    !*******************************************************************************
    !
    ! This function approximates the following special fraction y of gamma-functions:
    !
    ! y = gamma(a+b) / gamma(a)
    !
    !      for  a >= 1, b >= 0
    !
    ! The approximation error is less than +-1.5 % over the whole range.
    !
    ! The approximation exploits the rule gamma(x+1) = x*gamma(x)
    ! in both nominator and denominator recursively to arrive
    ! at an x' between 1 and 2. The remaining gamma-function of this x'
    ! is approximated by a simple parabola.
    !
    ! The lower a and b, the less terms have to be computed. Therefore
    ! efficiency decreases with increasing a and b. However, if a <= 10
    ! and b <= 2, the performance increases by a factor of more than 10 compared
    ! to using gfct() from numerical recipes for the nominator and denominator.
    !
    !*******************************************************************************

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: a, b
    REAL(KIND=dp)             :: x
    INTEGER                   :: ii, jj, i, j
    REAL(KIND=dp), PARAMETER  :: c = -0.113773074547242d0    ! gamma(1.5) - 1
    REAL(KIND=dp), PARAMETER  :: c1 = 1.9101845963779360d0   ! 1 - 8*c
    REAL(KIND=dp), PARAMETER  :: c2 = -1.3652768945669040d0  ! 12*c
    REAL(KIND=dp), PARAMETER  :: c3 = 0.45509229818896801d0  ! -4*c

    REAL(KIND=dp)             :: p, xx
    REAL(KIND=dp)             :: nom, denom

    ! statement function: parabola interpolation of the gamma function using
    !                     the nodes xx = 1, 1.5, 2, usable for the range xx in [1,2]:
! original Newton's form by "divided differences":
!    p(xx) = 1.0d0 + (xx-1.0d0)*(2d0*c - (xx-1.5d0)*4d0*c)
! after reorganization:
    p(xx) = c1 + c2*xx + c3*xx*xx

    ii = INT(a+b) - 1
    jj = INT(a)  - 1

    nom = 1.0d0
    DO i=1, ii
      nom = nom * (a+b-i)
    END DO
    denom = 1.0d0
    DO j=1, jj
      denom = denom * (a-j)
    END DO

    x = nom*p(a+b-ii) / (denom*p(a-jj))

    RETURN
  END FUNCTION gfct_frac_approx

  FUNCTION gfct_approx(a) RESULT(x)
    
    !*******************************************************************************
    !
    ! This function approximates the gamma-function
    !
    ! x = gamma(a)
    !
    !      for  a >= 1
    !
    ! The approximation error is less than + 1.0 % - 0.5 % over the whole range.
    !
    ! The approximation exploits the rule gamma(x+1) = x*gamma(x)
    ! recursively to arrive at an x' between 1 and 2.
    ! The remaining gamma-function of this x'
    ! is approximated by a simple parabola.
    !
    ! The lower a, the less terms have to be computed. Therefore
    ! efficiency decreases with increasing a. However, if a <= 12 or so,
    ! the performance increases by a factor of more than 10 compared
    ! to using gfct() from numerical recipes.
    !
    !*******************************************************************************
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: a
    REAL(KIND=dp)             :: x
    INTEGER                      :: ii, i
    REAL(KIND=dp), PARAMETER  :: c  = -0.113773074547242d0   ! gamma(1.5) - 1
    REAL(KIND=dp), PARAMETER  :: c1 = 1.9101845963779360d0   ! 1 - 8*c
    REAL(KIND=dp), PARAMETER  :: c2 = -1.3652768945669040d0  ! 12*c
    REAL(KIND=dp), PARAMETER  :: c3 = 0.45509229818896801d0  ! -4*c

    REAL(KIND=dp)             :: p, xx
    REAL(KIND=dp)             :: nom

    ! statement function: parabola interpolation of the gamma function using
    !                     the nodes xx = 1, 1.5, 2, usable for the range xx in [1,2]:
! original Newton's form by "divided differences":
!    p(xx) = 1.0d0 + (xx-1.0d0)*(2d0*c - (xx-1.5d0)*4d0*c)
! after reorganization:
    p(xx) = c1 + c2*xx + c3*xx*xx

    ii = INT(a) - 1

    nom = 1.0d0
    DO i=1, ii
      nom = nom * (a-i)
    END DO

    x = nom*p(a-ii)

    RETURN
  END FUNCTION gfct_approx

  !*******************************************************************************
  !                                                                              *
  !       Incomplete gamma function                                              *
  !                                                                              *
  !*******************************************************************************

  !*******************************************************************************
  ! 1) some helper functions

  SUBROUTINE gamma_help_cf(gammcf,a,x,gln)

    REAL(KIND=dp), INTENT(in)  :: a, x
    REAL(KIND=dp), INTENT(out) :: gammcf, gln

    INTEGER,       PARAMETER   :: maxiter = 100
    REAL(KIND=dp), PARAMETER   :: eps = 3.d-7, fpmin = 1.d-30
    INTEGER                    :: i
    REAL(KIND=dp)              :: an, b, c, d, del, h

    gln = gammln(a)
    b   = x + 1.0_dp - a
    c   = 1.0_dp / fpmin
    d   = 1.0_dp / b
    h   = d
    DO i = 1, maxiter
      an = -i*(i-a)
      b  = b + 2.0_dp
      d  = an*d + b
      IF (ABS(d) < fpmin) d = fpmin
      c  = b + an/c
      IF (ABS(c) < fpmin) c = fpmin
      d  = 1.0_dp / d
      del= d * c
      h  = h * del
      IF (ABS(del-1.0_dp) < EPS) EXIT
    END DO

    IF (ABS(del-1.0_dp) >= eps) THEN
      WRITE (*,*) 'ERROR  radar_gamma_functions_vec.f90 in GAMMA_HELP_CF: a too large, maxiter too small'
      gammcf = 0.0_dp
    ELSE
      gammcf = EXP(-x + a*LOG(x) - gln) * h
    END IF

  END SUBROUTINE gamma_help_cf

  SUBROUTINE gamma_help_ser(gamser,a,x,gln)

    REAL(KIND=dp), INTENT(in)  :: a, x
    REAL(KIND=dp), INTENT(out) :: gamser, gln

    INTEGER,       PARAMETER   :: maxiter = 100
    REAL(KIND=dp), PARAMETER   :: eps = 3.d-7
    INTEGER                    :: n
    REAL(KIND=dp)              :: ap,del,sum

    gln = gammln(a)
    IF (x <= 0.0_dp) THEN
      IF (x < 0.0_dp) THEN
        WRITE (*,*) 'ERROR radar_gamma_functions_vec.f90 in GAMMA_HELP_SER: x < 0'
      END IF
      
      gamser = 0.0_dp
    
    ELSE

      ap  = a
      sum = 1.0_dp/a
      del = sum
      DO n = 1, maxiter
        ap  = ap + 1.0_dp
        del = del*x/ap
        sum = sum + del
        IF (ABS(del) < ABS(sum)*eps) EXIT
      END DO

      IF (ABS(del) >= ABS(sum)*eps) THEN
        WRITE (*,*) 'ERROR radar_gamma_functions_vec.f90 in GAMMA_HELP_SER: a too large, maxiter too small'
        gamser = 0.0_dp
      ELSE
        gamser = sum * EXP(-x + a*LOG(x) - gln)
      END IF

    END IF

  END SUBROUTINE gamma_help_ser

  REAL(KIND=dp) FUNCTION gamma_p(a,x,gln)

    REAL(KIND=dp), INTENT(in)  :: a, x
    REAL(KIND=dp), INTENT(out) :: gln
    REAL(KIND=dp)              :: gammcf, gamser

    IF (x < 0.0_dp .OR. a <= 0.0_dp) THEN
      WRITE(*,*) 'ERROR radar_gamma_functions_vec.f90 in GAMMA_P: bad arguments'

      gamma_p = 0.0d0

    ELSE
    
      IF (x < a+1.0_dp) THEN
        CALL gamma_help_ser(gamser,a,x,gln)
        gamma_p = gamser
      ELSE
        CALL gamma_help_cf(gammcf,a,x,gln)
        gamma_p = 1.0_dp - gammcf
      ENDIF
      
    END IF

  END FUNCTION gamma_p

  REAL(KIND=dp) FUNCTION gamma_q(a,x,gln)

    REAL(KIND=dp), INTENT(in)  :: a, x
    REAL(KIND=dp), INTENT(out) :: gln
    REAL(KIND=dp)              :: gammcf, gamser

    IF (x < 0.0_dp .OR. a <= 0.0_dp) THEN
      WRITE(*,*) 'ERROR radar_gamma_functions_vec.f90 in GAMMA_Q: bad arguments'

      gamma_q = 0.0_dp

    ELSE

      IF (x < a+1.0_dp) THEN
        CALL gamma_help_ser(gamser,a,x,gln)
        gamma_q = 1.0_dp - gamser
      ELSE
        CALL gamma_help_cf(gammcf,a,x,gln)
        gamma_q = gammcf
      ENDIF
    
    END IF

  END FUNCTION gamma_q

  ! End of helper functions
  !*******************************************************************************

  !*******************************************************************************
  !
  ! Incomplete gamma-function, upper integral
  !
  !              int(x)(oo) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  REAL(KIND=dp) FUNCTION incgfct_upper(a,x)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: a, x
    REAL(KIND=dp) :: gam, gln

    gam = gamma_q(a,x,gln)
    incgfct_upper = EXP(gln) * gam

  END FUNCTION incgfct_upper

  !*******************************************************************************
  !
  ! Incomplete gamma-function, lower integral
  !
  !              int(0)(x) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  REAL(KIND=dp) FUNCTION incgfct_lower(a,x)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: a, x
    REAL(KIND=dp) :: gam, gln

    gam = gamma_p(a,x,gln)
    incgfct_lower = EXP(gln) * gam

  END FUNCTION incgfct_lower

  !*******************************************************************************
  !
  ! Incomplete gamma-function, integral
  !
  !              int(x1)(x2) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  REAL(KIND=dp) FUNCTION incgfct(a,x1,x2)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: a, x1, x2

    incgfct = incgfct_lower(a,x2) - incgfct_lower(a,x1)

  END FUNCTION incgfct

  !*******************************************************************************
  !
  ! Create Lookup-table vectors for the lower incomplete gamma function,
  !              int(0)(x) exp(-t) t^(a-1) dt
  ! as function of x at constant a.
  ! The table runs from x=0 to the 99.5 % - value of the normalized
  ! incomplete gamma function. This 99.5 % - value has been fitted
  ! with high accuracy as function of a in the range a in [0;20], but can
  ! safely be applied also to higher values of a. (Fit created with the
  ! matlab-program "gamma_unvoll_lower_lookup.m" by Ulrich Blahak, 2008/11/13).
  !
  ! The last value in the table corresponds to x = infinity, so that
  ! during the reconstruction of incgfct-values from the table,
  ! the x-value can safely be truncated at the maximum table x-value.
  !
  !*******************************************************************************

  SUBROUTINE incgfct_lower_lookupcreate(a,ltable,nl,nlhr)
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: a  ! value of a
    TYPE(gamlookuptable), intent(inout) :: ltable
    INTEGER, INTENT(in) :: nl, nlhr

    INTEGER :: i, err

    REAL(KIND=dp), PARAMETER ::   &
         c1 =  36.629433904824623d0, &
         c2 = -0.119475603955226d0,  &
         c3 =  0.339332937820052d0,  &
         c4 =  1.156369000458310d0

    ! Store parameters in the structure ltable:
    ltable%a = a
    ltable%n = nl
    ltable%nhr = nlhr

    ! Allocate Memory for the table vectors:
    NULLIFY(ltable%x)
    NULLIFY(ltable%xhr)
    NULLIFY(ltable%igf)
    NULLIFY(ltable%igfhr)

    ALLOCATE(ltable%x(nl), STAT=err)
    IF (err /= 0) THEN
      WRITE (*,*) 'INCGFCT_LOWER_LOOKUPCREATE: Allocation error x'
      STOP
    END IF
    ALLOCATE(ltable%xhr(nlhr), STAT=err)
    IF (err /= 0) THEN
      WRITE (*,*) 'INCGFCT_LOWER_LOOKUPCREATE: Allocation error xhr'
      STOP
    END IF
    ALLOCATE(ltable%igf(nl), STAT=err)
    IF (err /= 0) THEN
      WRITE (*,*) 'INCGFCT_LOWER_LOOKUPCREATE: Allocation error igf'
      STOP
    END IF
    ALLOCATE(ltable%igfhr(nlhr), STAT=err)
    IF (err /= 0) THEN
      WRITE (*,*) 'INCGFCT_LOWER_LOOKUPCREATE: Allocation error igfhr'
      STOP
    END IF

    !==================================================================
    ! low resolution part of the table:
    !==================================================================

    ! maximum x-value of the lookup table (99.5-%-value):
    ltable%x(ltable%n-1) = c1 * ( 1.0d0 - EXP(c2*a**c3) ) + c4*a

    ! create lookup table vectors:
    ltable%dx = ltable%x(ltable%n-1) / (ltable%n-2.0d0)
    ltable%odx = 1.0d0 / ltable%dx
!!! This loop does not vectorize because of incgfct_lower():
    DO i = 1, ltable%n - 1
      ltable%x(i) = (i-1) * ltable%dx
      ltable%igf(i) = incgfct_lower(a,ltable%x(i))
    END DO

    ! The last value is for x = infinity:
    ltable%x(ltable%n) = (ltable%n-1) * ltable%dx
    ltable%igf(ltable%n) = gfct(a)

    !==================================================================
    ! high resolution part of the table (lowest 1 % of the X-values):
    !==================================================================

    ! create lookup table vectors:
    ltable%dxhr = ltable%x(NINT(0.01*(ltable%n-1))) / (ltable%nhr-1.0d0)
    ltable%odxhr = 1.0d0 / ltable%dxhr
!!! This loop does not vectorize because of incgfct_lower():
    DO i = 1, ltable%nhr
      ltable%xhr(i) = (i-1) * ltable%dxhr
      ltable%igfhr(i) = incgfct_lower(a,ltable%xhr(i))
    END DO

    RETURN
  END SUBROUTINE incgfct_lower_lookupcreate

  !*******************************************************************************
  !
  ! Retrieve values from a lookup table of the lower incomplete gamma function,
  ! as function of x at a constant a, for which the lookup table has been
  ! created. Before the below function can be applied, the table itself
  ! has to be created by function incgfct_lower_lookupcreate(a,ltable,nl,nlhr) above !!!
  !
  ! The last value in the table has to correspond to x = infinity, so that
  ! during the reconstruction of incgfct-values from the table,
  ! the x-value can safely be truncated at the maximum table x-value:
  !
  ! ltable%igf( ltable%x(ltable%n),...) = gfct(a)
  !
  ! Profiling with ifort on a Linux-PC shows, that table lookup for the
  ! incompl. gamma-Funktion is faster by a factor of about 15 compared
  ! to the original function without optimization (-O0). Using optimization
  ! could change this ratio (we encoutered up to 300 depending on function inlining).
  !
  ! Concerning the accuracy, comparisons show that the results of table lookup
  ! are accurate to within better than 0.1 % or even much less, except for
  ! very small values of X, for which the absolute values are however very
  ! close to 0. For X -> infinity (X > 99.5 % - value), accuracy may be
  ! somewhat reduced up to about 0.5 % ,
  ! because the table is truncated at the 99.5 % value (second-last value)
  ! and the last value is set to the ordinary gamma function.
  !
  !*******************************************************************************

  REAL(KIND=dp) FUNCTION incgfct_lower_lookup(x, ltable)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: x  ! value of x for table lookup
    TYPE(gamlookuptable), INTENT(in) :: ltable

    INTEGER :: iu, io
    REAL(KIND=dp) :: xt

    ! Trunkcate x to the range of the table:
    xt = MAX(MIN(x, ltable%x(ltable%n)), 0.0d0)

    ! calculate indices of the neighbouring regular x-values
    ! in the table:
    iu = MIN(FLOOR(xt * ltable%odx) + 1, ltable%n-1)
    io = iu + 1

    ! interpolate linearily:
    incgfct_lower_lookup = ltable%igf(iu) + &
         (ltable%igf(io) - ltable%igf(iu)) * ltable%odx * (xt-ltable%x(iu))

    RETURN
  END FUNCTION incgfct_lower_lookup


  !*******************************************************************************
  !
  ! Retrieve values of the upper incomplete gamma function
  ! from a lookup table of the lower incomplete gamma function,
  ! as function of x at a constant a, for which the lookup table has been
  ! created. Before the below function can be applied, the table itself
  ! has to be created by function incgfct_lower_lookupcreate(a,ltable,nl,nlhr) above !!!
  !
  ! The last value in the table has to correspond to x = infinity
  ! (the ordinary gamma function of a), so that
  ! during the reconstruction of incgfct-values from the table,
  ! the x-value can safely be truncated at the maximum table x-value:
  !
  ! ltable%igf( ltable%x(ltable%n),...) = gfct(a)
  !
  !*******************************************************************************

  REAL(KIND=dp) FUNCTION incgfct_upper_lookup(x, ltable)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: x  ! value of x for table lookup
    TYPE(gamlookuptable), INTENT(in) :: ltable

    INTEGER :: iu, io
    REAL(KIND=dp) :: xt

    ! Trunkcate x to the range of the table:
    xt = MAX(MIN(x, ltable%x(ltable%n)), 0.0d0)

    ! calculate indices of the neighbouring regular x-values
    ! in the table:
    iu = MIN(FLOOR(xt * ltable%odx) + 1, ltable%n-1)
    io = iu + 1

    ! interpolate lower inc. gamma function linearily and subtract from
    ! the ordinary gamma function to get the upper
    ! incomplete gamma function:
    incgfct_upper_lookup = ltable%igf(ltable%n) - ltable%igf(iu) -  &
         (ltable%igf(io) - ltable%igf(iu)) * ltable%odx * (xt-ltable%x(iu))

    ! Due to truncation errors (differences of 2 almost identical numbers) may it happen
    ! that incgfct_upper_lookup(x, ltable) gets negative, but actually should be positive.
    ! This will happens in case of very large numbers of x.
    ! Limitate incgfct_upper_lookup to positive values!

    incgfct_upper_lookup = MAX(incgfct_upper_lookup, 0.0d0)

    RETURN
  END FUNCTION incgfct_upper_lookup

END MODULE radar_gamma_functions_vec
