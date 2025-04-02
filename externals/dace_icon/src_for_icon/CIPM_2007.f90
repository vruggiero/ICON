!
!+ Density of moist air (CIPM-2007)
!
! $Id$
!
MODULE CIPM_2007
!
! Description:
!   This module contains expressions for the calculation of the density
!   of moist air, and auxiliary functions (including their derivatives)
!   for the computation of the geopotential of non-ideal gases.
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_43        2015-08-19 Harald Anlauf
!  new module
!
! Code Description:
! Language: Fortran.
! Software Standards:
!
! Reference:
!   A. Picard, R.S. Davis, M. Gläser, K. Fujii,
!   Revised formula for the density of moist air (CIPM-2007),
!   Metrologia 45 (2008) 149-155.
!
!==============================================================================
  !==============
  ! Modules used:
  !==============
  use Defaults, only: double
  !----------------------------------------------------------
  implicit none
  !----------------------------------------------------------
  !================
  ! Public entities
  !================
  private

  public :: Z           ! Compressibility factor Z
  public :: Z_adj       ! Compressibility factor Z and (P,T,Q) derivatives
  public :: Zeta        ! Compressibility integral
  public :: Zeta_adj    ! Compressibility integral and (P,T,Q) derivatives

  !---------------------------
  ! Public physical constants:
  !---------------------------
  real(double), parameter, public :: &
       Md_cipm = 28.96546e-3_double, & ! [kg/mol] dry air mass @ x(CO2)=400ppm
       Mv_cipm = 18.01528e-3_double    ! [kg/mol] water vapor mass (IUPAC)

  real(double), parameter, public :: &
       R_2006  =  8.314472_double      ! [J/mol/K] CODATA 2006

  real(double), parameter, public :: &
       Rd_cipm = R_2006 / Md_cipm,   & ! [J/(kg*K)] Dry air gas constant
       Rv_cipm = R_2006 / Mv_cipm      ! [J/(kg*K)] Water vapor gas constant

  !====================
  ! Private parameters:
  !====================

  real(double), parameter :: t0c = 273.15_double

  !-------------------
  ! Derived constants:
  !-------------------
  real(double), parameter :: &      ! Choice
       Rd  = Rd_cipm,        &
       Rnu = Rv_cipm

  real(double), parameter :: &      ! T = Tv/(1 + Eps*q)
       Eps = Rnu/Rd - 1._double     ! Eps ~ 0.607827

  real(double), parameter :: &      ! xv = q/(aq + bq*q)
       aq  = Rd / Rnu,       &      ! aq ~ 0.621957
       bq  = 1._double - aq         ! bq ~ 0.378043

  !-----------------------------
  ! Compressibility coefficients
  ! CIPM-2007, (A1.4)
  !-----------------------------
  real(double), parameter ::    &
       a0 =  1.58123e-6_double, &
       a1 = -2.9331e-8_double,  &
       a2 =  1.1043e-10_double, &
       b0 =  5.707e-6_double,   &
       b1 = -2.051e-8_double,   &
       c0 =  1.9898e-4_double,  &
       c1 = -2.376e-6_double,   &
       d0 =  1.83e-11_double,   &
       e0 = -0.765e-8_double

  !============================================================================
contains
  !============================================================================
  elemental function Z (P, T, Q)
    !----------------------------------------------------------
    ! Calculation of the compressibility factor Z of moist air.
    !----------------------------------------------------------
    real(double), intent(in)  :: P    ! Pressure               [Pa]
    real(double), intent(in)  :: T    ! Temperature            [K]
    real(double), intent(in)  :: Q    ! Specific humidity      [kg/kg]
    real(double)              :: Z    ! Compressibility factor [dimensionless]

    !----------------
    ! Local variables
    !----------------
    real(double) :: tc            ! Temperature [°C]
    real(double) :: pt            ! P/T
    real(double) :: xv            ! Water vapor mixing ratio

    tc   =   T  - t0c
    pt   =   P  / T
    xv   =   Q  / (aq + bq*Q)

    !------------------
    ! CIPM-2007, (A1.4)
    !------------------
    Z   = 1._double &
          - pt    * (a0 + (a1+a2*tc)*tc + (b0+b1*tc)*xv + (c0+c1*tc)*xv**2) &
          + pt**2 * (d0 + e0*xv**2)

  end function Z
  !============================================================================
  elemental function Zeta (P, T, Q)
    !----------------------------------------------------
    ! Calculation of the compressibility integral
    !   Zeta(P,T,Q) = int_0^P (dp/p) (Z(p,T,Q)-1)
    !----------------------------------------------------
    real(double), intent(in)  :: P       ! Pressure                 [Pa]
    real(double), intent(in)  :: T       ! Temperature              [K]
    real(double), intent(in)  :: Q       ! Specific humidity        [kg/kg]
    real(double)              :: Zeta    ! Compressibility integral [dim.less]

    !----------------
    ! Local variables
    !----------------
    real(double) :: tc            ! Temperature [°C]
    real(double) :: pt            ! P/T
    real(double) :: xv            ! Water vapor mixing ratio

    tc   =   T  - t0c
    pt   =   P  / T
    xv   =   Q  / (aq + bq*Q)

    Zeta   = - pt    * (a0 + (a1+a2*tc)*tc + (b0+b1*tc)*xv + (c0+c1*tc)*xv**2) &
             + pt**2 * (d0 + e0*xv**2) * 0.5_double

  end function Zeta
  !============================================================================
  elemental subroutine Z_adj (P, T, Q, Z, Z_P, Z_T, Z_Q)
    !----------------------------------------------------------
    ! Calculation of the compressibility factor Z of moist air.
    ! Adjoint version.
    !----------------------------------------------------------
    real(double), intent(in)  :: P    ! Pressure               [Pa]
    real(double), intent(in)  :: T    ! Temperature            [K]
    real(double), intent(in)  :: Q    ! Specific humidity      [kg/kg]
    real(double), intent(out) :: Z    ! Compressibility factor [dimensionless]
    real(double), intent(out) :: Z_P  ! d(Z)/d(P)              [1/Pa]
    real(double), intent(out) :: Z_T  ! d(Z)/d(T)              [1/K]
    real(double), intent(out) :: Z_Q  ! d(Z)/d(Q)              [kg/kg]

    !----------------
    ! Local variables
    !----------------
    real(double) :: tc            ! Temperature [°C]
   !real(double) :: tc_t = 1      ! d(tc)/T
    real(double) :: pt            ! P/T
    real(double) :: pt_T          ! d(p/T)/d(T)
    real(double) :: pt_P          ! d(p/T)/d(P)
    real(double) :: xv            ! Water vapor mixing ratio
    real(double) :: xv_Q          ! d(xv)/d(Q)

    tc   =   T  - t0c
    pt   =   P  / T
    pt_T = - pt / T
    pt_P =   1  / T
    xv   =   Q  / (aq + bq*Q)
    xv_Q =   aq / (aq + bq*Q)**2

    !------------------
    ! CIPM-2007, (A1.4)
    !------------------
    Z   = 1._double &
          - pt    * (a0 + (a1+a2*tc)*tc + (b0+b1*tc)*xv + (c0+c1*tc)*xv**2) &
          + pt**2 * (d0 + e0*xv**2)

    Z_T = - pt_T  * (a0 + (a1+a2*tc)*tc + (b0+b1*tc)*xv + (c0+c1*tc)*xv**2) &
          - pt    * (      a1+ 2*a2 *tc +     b1    *xv +     c1    *xv**2) &
          + 2*pt  * (d0 + e0*xv**2) * pt_T

    Z_P = - pt_P  * (a0 + (a1+a2*tc)*tc + (b0+b1*tc)*xv + (c0+c1*tc)*xv**2) &
          + 2*pt  * (d0 + e0*xv**2) * pt_P

    Z_Q = (- pt    *                     ((b0+b1*tc) +  2*(c0+c1*tc)*xv)    &
           + pt**2 *    2*e0*xv   ) * xv_Q

  end subroutine Z_adj
  !============================================================================
  elemental subroutine Zeta_adj (P, T, Q, Zeta, Zeta_P, Zeta_T, Zeta_Q)
    !----------------------------------------------------
    ! Calculation of the compressibility integral
    !   Zeta(P,T,Q) = int_0^P (dp/p) (Z(p,T,Q)-1)
    ! and its derivatives w.r.t. P, T, and Q.
    !----------------------------------------------------
    real(double), intent(in)  :: P       ! Pressure                 [Pa]
    real(double), intent(in)  :: T       ! Temperature              [K]
    real(double), intent(in)  :: Q       ! Specific humidity        [kg/kg]
    real(double), intent(out) :: Zeta    ! Compressibility integral [dim.less]
    real(double), intent(out) :: Zeta_P  ! d(Zeta)/d(P)             [1/Pa]
    real(double), intent(out) :: Zeta_T  ! d(Zeta)/d(T)             [1/K]
    real(double), intent(out) :: Zeta_Q  ! d(Zeta)/d(Q)             [kg/kg]

    !----------------
    ! Local variables
    !----------------
    real(double) :: tc            ! Temperature [°C]
   !real(double) :: tc_t = 1      ! d(tc)/T
    real(double) :: pt            ! P/T
    real(double) :: pt_T          ! d(p/T)/d(T)
    real(double) :: pt_P          ! d(p/T)/d(P)
    real(double) :: xv            ! Water vapor mixing ratio
    real(double) :: xv_Q          ! d(xv)/d(Q)

    tc   =   T  - t0c
    pt   =   P  / T
    pt_T = - pt / T
    pt_P =   1  / T
    xv   =   Q  / (aq + bq*Q)
    xv_Q =   aq / (aq + bq*Q)**2

    !-----------
    ! C.f. Z_adj
    !-----------
    Zeta   = - pt    * (a0 + (a1+a2*tc)*tc + (b0+b1*tc)*xv + (c0+c1*tc)*xv**2) &
             + pt**2 * (d0 + e0*xv**2) * 0.5_double

    Zeta_T = - pt_T  * (a0 + (a1+a2*tc)*tc + (b0+b1*tc)*xv + (c0+c1*tc)*xv**2) &
             - pt    * (      a1+ 2*a2 *tc +     b1    *xv +     c1    *xv**2) &
             + pt    * (d0 + e0*xv**2) * pt_T

    Zeta_P = - pt_P  * (a0 + (a1+a2*tc)*tc + (b0+b1*tc)*xv + (c0+c1*tc)*xv**2) &
             + pt    * (d0 + e0*xv**2) * pt_P

    Zeta_Q = (-pt    *                      ((b0+b1*tc) +  2*(c0+c1*tc)*xv)    &
             + pt**2 *       e0*xv   ) * xv_Q

  end subroutine Zeta_adj
  !============================================================================
end MODULE CIPM_2007
