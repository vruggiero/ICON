!
!+ Conversion routines for scalar quantities (mostly moisture)
!
!==============================================================================
!> Conversion routines for scalar quantities (mostly moisture)
!>
!>  This module provides:
!>  - Mathematical constants
!>  - Physical constants                         (originally taken from DM)
!>  - conversion routines for moist physics      (elemental functions)
!>  - adjoint conversion routines                (elemental subroutines)
!>
!>  Some subroutines or functions are provided in different flavours to
!>  allow vectorisation (on NEC SX).
!>
!> \todo move constants to module mo_constants
!> \todo move routines which depend on a model grid to another module
!> \todo do not use SATWAT preprocessor definition any more
!>
MODULE mo_physics
!
! Description:
!   This module provides:
!   - Mathematical constants
!   - Physical constants
!   - conversion routines for moist physics      (elemental functions)
!     currently DM (Deutschlandmodell) constants
!   - adjoint conversion routines                (elemental subroutines)
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
! V1_5         2009/05/25 Harald Anlauf
!  tq_tvrh_vec: vectorized variant of tq_tvrh
! V1_7         2009/08/24 Andreas Rhodin
!  add directives for doxygen;
!  change c_p from 1004 to 1005 for consistence with COSMO
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  q_rhw, rhw_q_adj, tq_tvrh: (optionally) use relative humidity over WATER
! V1_13        2011/11/01 Detlef Pingel
!  additional constants for IASI
! V1_22        2013-02-13 Harald Anlauf
!  changes for ICON: set_t, transformations for virtual potential temperature
!  changed value of c_p, set_p for compliance with ICON
! V1_28        2014/02/26 Harald Anlauf
!  Dry air gas constant: GME /= ICON, use #ifdef __ICON__
! V1_35        2014-11-07 Harald Anlauf
!  use proper constants introduced in rev.10717;
!  new q_rhi: specific humidity from relative humidity over ice
! V1_45        2015-12-15 Andreas Rhodin
!  tq_rh: choose convergence criterium less tight
! V1_46        2016-02-05 Michael Bender
!  some routines moved from the STD operator to mo_physics and mo_algorithms
! V1_47        2016-06-06 Harald Anlauf
!  Fix comment to agree with actual code
! V1_48        2016-10-06 Robin Faulwetter
!  conversion routines for ppmv  (for RTTOV12)
! V1_50        2017-01-09 Harald Anlauf
!  OpenMP parallelization of tq_tvrh_vec
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM 2000  original code
! Detlef  Pingel  DWD   2004  bug fixes, new subroutine tq_tvrh
! Andreas Rhodin  DWD   2004  split module into mo_physics and mo_cntrlvar
! Harald Anlauf   DWD   2008  fixes for more consistent cross-platform results
!==============================================================================
!!
!! define  SATWAT to use saturation water vapour over water
!!
#define SATWAT

  !=============
  ! modules used
  !=============
  use mo_kind,        only: wp       ! real working precision kind parameter
  use data_constants, only: pi,     &! 3.1415...
                            r_v,    &! gas constant for water vapor
                            r_earth  ! mean radius of the earth (m)
  use mo_dace_string, only: tolower
  implicit none
  !================
  ! public entities
  !================
  private
  !----------
  ! constants
  !----------
  public :: pi        !          3.1415...
  public :: e         !          2.718281828459045
  public :: d2r, r2d  !          conversion factors: degree <-> radians
  public :: d2x, x2d  !          conversion factors: degree <-> meter
  public :: R         !          gas constant of dry air [J/(kg*K)]
  public :: Rd        !          gas constant of water vapour
  public :: RdRd      !          R/Rd
  public :: gacc      !          gravity acceleration
  public :: rearth    !          earth radius
  public :: t0c       !          T [K] for T=0 degree Celsius
  public :: c_p       !          specific heat           [J/(kg*K)]
  public :: lapse_cl  !          Climatological lapse rate (-dT/dz) [K/m]
  public :: RDDRM1    !          for NEC SX6 workaround
  public :: B1        !>
  public :: B2W       ! >
  public :: B2E       !  >       constants used in magnus formula
  public :: B2I       !  >
  public :: B3        !  >
  public :: B4W       ! >
  public :: B4E       !>
  public :: B4I       !>
  public :: EMRDRD    !          1. - RD/RD
  public :: p0ref     !          Reference pressure (Exner function,pot.temp.)
  public :: gas_id_o3 !
  public :: gas_id_co2!
  public :: gas_id_n2o!
  public :: gas_id_co !
  public :: gas_id_ch4!
  public :: gas_id_so2!
  public :: ngases    !
  public :: RTRG
  public :: get_gas

  !------------------------------
  ! conversion routines, adjoints
  !------------------------------
  public :: es_t           !         saturation vapour pressure <-temperature
  public :: esw_t          !         saturation vapour pressure over water
  public :: esi_t          !         saturation vapour pressure over ice
  public :: desw_dt        !         derivative of esw_t
  public :: esw_t_hardy    !         like esw_t, but using Hardy (1998)
  public :: e_q            !         water vap.press.  <- specific humidity
  public :: t_es           !         temperature       <- sat. water vapour press.
  public :: t_esw          !         temperature       <- sat. pressure over water
  public :: q_e            !         specific humidity <- water vapour pressure
  public :: p_ps           !         pressure          <- surface pressure
  public :: p_ps_adj       !adjoint: pressure          <- surface pressure
  public :: rh_q           !         relative humidity <- specific humidity
  public :: rh_q_adj       !adjoint: relative humidity <- specific humidity
  public :: rhw_q          !         rel.hum. (water)  <- specific humidity
  public :: rhw_q_adj      !adjoint: rel.hum. (water)  <- specific humidity
  public :: rhw_q_hardy    !         rel.hum. (water)  <- specific humidity
  public :: rhw_m          !         rel.hum. (water)  <- mixing ratio
  public :: rhw_m_hardy    !         rel.hum. (water)  <- mixing ratio
  public :: q_rh           !         specific humidity <- relative humidity
  public :: q_rh_adj       !adjoint: specific humidity <- relative humidity
  public :: q_rhw          !         specific humidity <- relative humidity over water
  public :: q_rhi          !         specific humidity <- relative humidity over ice
  public :: ppmv_dry_q     !         ppmv (dry air)    <- specific humidity
  public :: dppmv_dry_dq   !deriv.:  ppmv (dry air)    <- specific humidity
  public :: q_ppmv_dry     !         ppmv (dry air)    -> specific humidity
  public :: ppmv_moist_q   !         ppmv (moist air)  <- specific humidity
  public :: dppmv_moist_dq !deriv.:  ppmv (moist air)  <- specific humidity
  public :: q_ppmv_moist   !         ppmv (moist air)  -> specific humidity
  public :: ppmv_dry_trg   !         ppmv (dry air)    <- spec. conc. of tracegas
  public :: ppmv_moist_trg !         ppmv (moist air)  <- spec. conc. of tracegas
  public :: q_ppmv_dry_trg !         ppmv (dry air)    -> spec. conc. of tracegas
  public :: tv_t_q         !         virtual temperat. <- temp., spec.hum., liq/ice
  public :: tv_t_q_adj     !adjoint: virtual temperat. <- temp., spec.hum., liq/ice
  public :: tv_t_rh        !         virtual temperat. <- temp.,  rel.hum., liq/ice
  public :: tv_t_rh_adj    !adjoint: virtual temperat. <- temp.,  rel.hum., liq/ice
  public :: t_tv_q         !         temperature       <- virt.temp., spec.hum.
  public :: t_tv_q_adj     !adjoint: temperature       <- virt.temp., spec.hum.
  public :: tq_tvrh        !         temp., spec.hum.  <- virt.temp., rel.hum.
  public :: tq_tvrh_vec    !         vectorized variant of tq_tvrh (for SX-9)
  public :: tq_tvrh_vec2   !         vectorized variant of tq_tvrh (for SX Aurora)
  public :: tq_tvrh_noadj, tq_tvrh_noxfail           !+++ work around NEC SX bug
  public :: t_tq_tvrh      !         debug info derived type
  ! public :: uv_fd          !         u,v from wind speed and direction
  public :: fd_uv          !         wind speed and direction from u,v
  public :: exner          !         Exner function
  public :: tv_thetav_p    !         virtual temperat. <- virt.pot.temp., pressure
  public :: thetav_tv_p    !         virtual pot.temp. <- virt.    temp., pressure
  public :: p_rho_thetav   !       pressure          <- density, virt.pot.temp.
  public :: geopot_geoid   !       geopotential height -> height above geoid
  public :: geoid_geopot   !       height above geoid -> geopotential height
  !------------------------------------------------------
  ! make specific routines public to be traced by doxygen
  !------------------------------------------------------
  public :: ph_pf_ps_1
  public :: ph_pf_ps_3

  public :: R_RTTOV

!==============================================================================

  !--------------------------------------------------
  ! physical constants (originally taken from the DM)
  !--------------------------------------------------
  !
  !>  general gas constant [J/(mol*K)]
  !
  real(wp), parameter :: RGEN   = 8.31446261815324_wp
  !
  !>  gas constant of dry air [J/(kg*K)]
  !
#ifdef __ICON__
  real(wp), parameter :: R      = 287.04_wp     ! COSMO/ICON
#else
  real(wp), parameter :: R      = 287.05_wp     ! HRM/GME
#endif
  real(wp), parameter :: R_RTTOV = RGEN/0.0289644_wp ! Only used for tracegas unit conversion
  !
  !>  gas constant of water vapor [J/(kg*K)]
  !
  real(wp), parameter :: RD     = r_v           ! 461.51 as in COSMO
  !
  !>  R/RD
  !
  real(wp), parameter :: RDRD   = R/RD
  !
  !>  RD/R - 1.
  !
  real(wp), parameter :: RDDRM1 = RD/R - 1._wp
  !
  !>  1. - RD/RD
  !
  real(wp), parameter :: EMRDRD = 1._wp - RDRD
  !
  !>  constant used in Magnus formula [Pa]
  !
  real(wp), parameter :: B1     = 610.78_wp
  !
  !>  constant used in Magnus formula [K]
  !
  real(wp), parameter :: B3     = 273.16_wp
  !
  !>  constant used in Magnus formula (water)
  !
  real(wp), parameter :: B2W    =  17.2693882_wp
  !
  !>  constant used in Magnus formula (water) [K]
  !
  real(wp), parameter :: B4W    =  35.86_wp
  !
  !>  constant used in Magnus formula (ice)
  !
  real(wp), parameter :: B2I    =  21.8745584_wp
  !
  !>  constant used in Magnus formula (ice)
  !
  real(wp), parameter :: B4I    =   7.66_wp

#ifdef SATWAT
  !>  constant used in Magnus formula (ice) forced to water!!
  !
  real(wp), parameter :: B2E    = B2W
  !
  !>  constant used in Magnus formula (ice) forced to water!! [K]
  !
  real(wp), parameter :: B4E    = B4W
#else
  !>  constant used in Magnus formula (ice)
  !
  real(wp), parameter :: B2E    = B2I
  !
  !>  constant used in Magnus formula (ice)
  !
  real(wp), parameter :: B4E    = B4I
#endif
  !
  !>  radius of the earth
  !
! real(wp),parameter :: rearth  = 6371000.0_wp ! radius of the earth  (ECHAM)
! real(wp),parameter :: rearth  = 6371229.0_wp ! radius of the earth  (GME)
  real(wp),parameter :: rearth  = r_earth      ! 6371229.             (COSMO)
  !
  !>  gravity acceleration (ECHAM,GME)
  !
  real(wp),parameter :: gacc    = 9.80665_wp
  !
  !>  tmelt in GME
  !
  real(wp),parameter :: t0c     = 273.15_wp
  !
  !>  specific heat [J/(kg*K)]
  !
  real(wp),parameter :: c_p     = R * 3.5_wp   ! c_p/R = 7/2 (same as ICON)
!!real(wp),parameter :: c_p     = 1005._wp     ! old value taken from COSMO
  !
  !>  Climatological lapse rate [K/m]
  !
  real(wp),parameter :: lapse_cl= 6.5e-3_wp
  !
  !>  Reference pressure of Exner function (c.f. potential temperature) [Pa]
  !
  real(wp),parameter :: p0ref   = 100000._wp
  !
  !>  Inverse of reference pressure of Exner function
  !
  real(wp),parameter :: p0_inv  = 1._wp / p0ref
  !------------------
  ! general constants
  !------------------
  !
  !>  e
  !
  real(wp),parameter :: e       = 2.718281828459045235_wp
  !
  !>  pi / 180           (factor degree  -> radians)
  !
  real(wp),parameter :: d2r     = pi     /180._wp ! factor degree  -> radians
  !
  !>  180 / pi           (factor radians -> degree)
  !
  real(wp),parameter :: r2d     = 180._wp/pi      ! factor radians -> degree
  !
  !>  pi / 180 * rearth  (factor degree  -> meter)
  !
  real(wp),parameter :: d2x     = d2r    * rearth ! factor degree  -> meter
  !
  !>  180 / pi / rearth  (factor meter   -> degree)
  !
  real(wp),parameter :: x2d     = r2d    / rearth ! factor meter   -> degree

  !> Trace gas IDs
  integer, parameter :: gas_id_o3  = 1
  integer, parameter :: gas_id_co2 = 2
  integer, parameter :: gas_id_n2o = 3
  integer, parameter :: gas_id_co  = 4
  integer, parameter :: gas_id_ch4 = 5
  integer, parameter :: gas_id_so2 = 6
  integer, parameter :: ngases     = 6
  !> Trace gas names
  character(len=3), parameter :: gas_names(ngases) = (/ 'o3 ', &
                                                        'co2', &
                                                        'n2o', &
                                                        'co ', &
                                                        'ch4', &
                                                        'so2' /)
  !> Trace gas gas-constants (taken from rttov_const (V13.1))
  real(kind=wp), parameter :: RTRG(ngases) = (/ RGEN/0.04799820_wp, &
                                                RGEN/0.04400950_wp, &
                                                RGEN/0.04401280_wp, &
                                                RGEN/0.02801010_wp, &
                                                RGEN/0.01604246_wp, &
                                                RGEN/0.06406400_wp /)

!==============================================================================
  !-----------
  ! interfaces
  !-----------
  !
  !>  derive p from ps for hybrid vertical coordinated (GME)
  !
  interface p_ps                ! derive p from ps
    module procedure ph_pf_ps_1 ! 1-d (column) version
    module procedure ph_pf_ps_3 ! 3-d (volume) version
  end interface
  !
  !> adjoint to p_ps
  !
  interface p_ps_adj               ! adjoint to p_ps
    module procedure ph_pf_ps_adj3 ! 3-d (volume) version
  end interface p_ps_adj
!==============================================================================
  !------------------------
  ! terived type definition
  !------------------------
  !
  !>  maximum number of iterations used in subroutine tq_tvrh
  !
  integer, parameter :: maxiter_tq_tvrh = 30
  !
  !>  return argument from tq_tvrh (for debugging purposes)
  !
  type t_tq_tvrh
    real(wp) :: tmin, tmax, qmin, qmax
    real(wp) :: t  (maxiter_tq_tvrh)
    real(wp) :: q  (maxiter_tq_tvrh)
    real(wp) :: tv (maxiter_tq_tvrh)
    real(wp) :: rh (maxiter_tq_tvrh)
  end type t_tq_tvrh
!==============================================================================
contains
!==============================================================================
  !---------------------------------------------------------------
  !>  Magnus formula for water
  !> \param[in] t temperature                                 [K]
  !> \return      saturation water vapour pressure over water [Pa]
  !---------------------------------------------------------------
  elemental function esw_t (t)  result (e)
  real(wp), intent (in) :: t ! temperature
  real(wp)              :: e ! saturation water vapour pressure over water
    e = b1 * exp (b2w * (t - b3) / (t - b4w))
  end function esw_t
!------------------------------------------------------------------------------
  elemental function desw_dt (t) result (de)
  !------------------------------------------------------------
  ! derivative of esw_t
  ! (of saturation vapour pressure with respect to temperature)
  !------------------------------------------------------------
  real(wp), intent (in)  :: t     ! temperature
  real(wp)               :: de    ! derivative

    real(wp) :: tb3, tb4, ex, e

!   e = b1 * exp (b2w * (t - b3) / (t - b4w)) ! original formular

    tb3 = t - b3
    tb4 = t - b4w
    ex  = b2w * tb3 / tb4
    e   = b1 * exp (ex)
    de  = e * b2w * (1._wp/tb4 - tb3/tb4**2)

  end function desw_dt
!------------------------------------------------------------------------------
  !-------------------------------------------------------------
  !>  Magnus formula for ice
  !> \param[in] t temperature                               [K]
  !> \return      saturation water vapour pressure over ice [Pa]
  !-------------------------------------------------------------
  elemental function esi_t (t)  result (e)
  real(wp) ,intent (in) :: t ! temperature
  real(wp)              :: e ! saturation water vapour pressure over ice
     e = b1 * exp (b2i * (t - b3) / (t - b4i))
  end function esi_t
!==============================================================================
  !----------------------------------------------------------------------
  !> Magnus formula over ice or water
  !> \param[in] t temperature                                        [K]
  !> \return      saturation water vapour pressure over ice or water [Pa]
  !----------------------------------------------------------------------
  elemental function es_t (t)  result (e)
  real(wp) ,intent (in) :: t ! temperature
  real(wp)              :: e ! sat. water vap. pressure
    if (t >= b3) then
      e = b1 * exp (b2w * (t - b3) / (t - b4w))
    else
      e = b1 * exp (b2e * (t - b3) / (t - b4e))
    endif
  end function es_t
!==============================================================================
  !----------------------------------------------------
  !> inverse Magnus formula
  !> \param[in] e saturation water vapour pressure [Pa]
  !> \return      temperature [K], (0 for e <= 0)  [K]
  !----------------------------------------------------
  elemental function t_es (e)  result (te)
  real(wp) ,intent (in) :: e  ! sat. water vapour pressure
  real(wp)              :: te ! temperature
    real(wp) :: aux
    if (e > 0._wp) then
      aux = log (e / b1)
      te = max ( (b2w * b3 - b4w * aux) / (b2w - aux), &
                 (b2e * b3 - b4e * aux) / (b2e - aux))
    else
      te = 0._wp
    endif
  end function t_es
!==============================================================================
  !----------------------------------------------------
  !> inverse Magnus formula (over water)
  !> \param[in] e saturation water vapour pressure [Pa]
  !> \return      temperature                      [K]
  !----------------------------------------------------
  elemental function t_esw (e)  result (t)
  real(wp) ,intent (in) :: e ! saturation water vapour pressure over water
  real(wp)              :: t ! temperature
    real(wp) :: aux
    aux = log (e / b1)
    T = (B2W * B3 - B4W * aux) / (B2W - aux)
  end function t_esw
!==============================================================================
  !--------------------------------------------------------
  !>  Hardy (1998) ITS-90 formulation for vapor pressure
  !> \param[in] t temperature                          [K]
  !> \return      saturation vapor pressure over water [Pa]
  !--------------------------------------------------------
  elemental function esw_t_hardy (t) result (e)
    real(wp), intent (in) :: t  ! temperature
    real(wp)              :: e  ! saturation vapor pressure over water

    ! Wexler's formula using Hardy's coefficients:
    !          6
    ! ln e_s = ∑ g_i T**(i−2) + g_7 ln T
    !         i=0
    real(wp), parameter :: g0 = -2.8365744e3_wp, &
                           g1 = -6.028076559e3_wp, &
                           g2 =  1.954263612e1_wp, &
                           g3 = -2.737830188e-2_wp, &
                           g4 =  1.6261698e-5_wp, &
                           g5 =  7.0229056e-10_wp, &
                           g6 = -1.8680009e-13_wp, &
                           g7 =  2.7150305_wp
    real(wp) :: tinv

    tinv = 1 / t
    e = exp ((g0 * tinv + g1) * tinv + g2 +                        &
             (((g6 * t + g5) * t + g4) * t + g3) * t + g7 * log (t))
  end function esw_t_hardy
!==============================================================================
  !-----------------------------------------------
  !>  specific humidity from water vapour pressure
  !> \param[in] e water vapour pressure [Pa]
  !> \param[in] p air pressure          [Pa]
  !> \return      specific humidity     [kg/kg]
  !-----------------------------------------------
  elemental function q_e (e, p)   result (q)
  real(wp) ,intent (in) :: e ! water vapour pressure
  real(wp) ,intent (in) :: p ! air pressure
  real(wp)              :: q ! specific humidity
    q = rdrd * e / (p - emrdrd * e)
  end function q_e
!==============================================================================
  !-----------------------------------------------
  !>  water vapour pressure from specific humidity
  !> \param[in] q specific humidity      [kg/kg]
  !> \param[in] p air pressure           [Pa]
  !> \return      water vapour pressure  [Pa]
  !-----------------------------------------------
  elemental function e_q (q, p)   result (e)
  real(wp), intent (in) :: q ! specific humidity
  real(wp), intent (in) :: p ! air pressure
  real(wp)              :: e ! water vapour pressure
    e = q * p / (rdrd + emrdrd * q)
  end function e_q
!==============================================================================
  !-----------------------------------------------------------
  !> specific humidity from relative humidity, adjoint routine
  !
  !> adjoints are added to arguments rh_a, t_a, p_a.
  !> q_a is set to zero on return
  !
  !> \param [in,out] q_a  specific humidity, adjoint variable
  !> \param [in,out] rh_a relative humidity, adjoint variable
  !> \param [in,out] t_a  temperature,       adjoint variable
  !> \param [in,out] p_a  pressure,          adjoint variable
  !> \param [in]     rh   relative humidity  [0..1]
  !> \param [in]     t    temperature        [K]
  !> \param [in]     p    pressure           [Pa]
  !-----------------------------------------------------------
  elemental subroutine q_rh_adj (q_a, rh_a, t_a, p_a, rh, t, p)
  real(wp) ,intent (inout) :: q_a        ! specific humidity, adjoint variable
  real(wp) ,intent (inout) :: rh_a       ! relative humidity, adjoint variable
  real(wp) ,intent (inout) :: t_a        ! temperature,       adjoint variable
  real(wp) ,intent (inout) :: p_a        ! pressure,          adjoint variable
  real(wp) ,intent (in)    :: rh         ! relative humidity
  real(wp) ,intent (in)    :: t          ! temperature
  real(wp) ,intent (in)    :: p          ! pressure
    real(wp) :: e, e_a, es, c1, c2, c3
!
! nonlinear code:
!
! E = RH * b1 * exp (b2w * (T - b3) / (T - b4w))
! Q = rdrd * E / (P - emrdrd * E)
!
! tangentlinear code (capital letters=nonl.variables, small l.=tl.variables):
!
! ES = b1 * exp (b2w * (T - b3) / (T - b4w))
! if (ES > P) ES = P                          ! saveguard for low p
! E  = RH * ES
! C2 = 1/(T-b4w)
! e  = rh * ES
! if (ES < P)
!   e = e + t  * E * b2w * C2 * (1 - (T-b3) * C2)
! else
!   e = e + RH * p
! endif
! C3 = 1/(P-emrdrd * E)
! q' = e' * rdrd * C3 * (1+E*emrdrd*C3)
!    - p' * rdrd * E * C3**2
!
! adjoint code:
!
    if (t >= b3) then
      es = min (b1 * exp (b2w * (t - b3) / (t - b4w)) ,p)
      c2 = 1._wp/ (t-b4w)
      c1 = b2w * c2
    else
      es = min (b1 * exp (b2e * (t - b3) / (t - b4e)) ,p)
      c2 = 1._wp/ (t-b4e)
      c1 = b2e * c2
    endif
    e    = rh * es
    c3   = 1._wp/(p-emrdrd*e)
    p_a  = p_a  - q_a * rdrd *e * c3**2
    e_a  =        q_a * rdrd * c3 * ( 1._wp + e * c3 * emrdrd)
    rh_a = rh_a + e_a * es
    if (es >= p) then
      p_a = p_a + e_a * rh
    else
      t_a = t_a + e_a * e * c1 * (1._wp - (t-b3)*c2)
    endif
    q_a  = 0._wp
  end subroutine q_rh_adj
!______________________________________________________________________________
  !-------------------------------------------
  !>  specific humidity from relative humidity
  !> \param [in] rh relative humidity  [0..1]
  !> \param [in] t  temperature        [K]
  !> \param [in] p  pressure           [Pa]
  !> \return        specific humidity  [kg/kg]
  !-------------------------------------------
  elemental function q_rh (rh, t, p) result (q)
  real(wp) ,intent (in) :: rh         ! relative humidity
  real(wp) ,intent (in) :: t          ! temperature
  real(wp) ,intent (in) :: p          ! pressure
  real(wp)              :: q          ! specific humidity
    real(wp) :: e ! water vapor pressure
!    q = rdrd * e / (p - emrdrd * e)
!    e = rH * es
!!   e = q * p / (rdrd + emrdrd * q)
!    if (t >= b3) then
!      es = b1 * exp (b2w * (t - b3) / (t - b4w))
!    else
!      es = b1 * exp (b2e * (t - b3) / (t - b4e))
!    endif
!    es = min (es, p)    ! saveguard for low p
    if (t >= b3) then
      e = rh * min (b1 * exp (b2w * (t - b3) / (t - b4w)) ,p)
    else
      e = rh * min (b1 * exp (b2e * (t - b3) / (t - b4e)) ,p)
    endif
    q = rdrd * e / (p - emrdrd * e)
  end function q_rh
!______________________________________________________________________________
  elemental function q_rhw (rh, t, p) result (q)
  !----------------------------------------------------
  ! specific humidity from relative humidity over water
  !----------------------------------------------------
  real(wp) ,intent (in) :: rh         ! relative humidity
  real(wp) ,intent (in) :: t          ! temperature
  real(wp) ,intent (in) :: p          ! pressure
  real(wp)              :: q          ! relative humidity
    real(wp) :: e ! water vapor pressure
    e = rh * min (b1 * exp (b2w * (t - b3) / (t - b4w)) ,p)
    q = rdrd * e / (p - emrdrd * e)
  end function q_rhw
!______________________________________________________________________________
  elemental function q_rhi (rh, t, p) result (q)
  !--------------------------------------------------
  ! specific humidity from relative humidity over ice
  !--------------------------------------------------
  real(wp) ,intent (in) :: rh         ! relative humidity
  real(wp) ,intent (in) :: t          ! temperature
  real(wp) ,intent (in) :: p          ! pressure
  real(wp)              :: q          ! relative humidity
    real(wp) :: e ! water vapor pressure
    e = rh * min (b1 * exp (b2i * (t - b3) / (t - b4i)) ,p)
    q = rdrd * e / (p - emrdrd * e)
  end function q_rhi
!==============================================================================
  !-------------------------------------------
  !>  relative humidity from specific humidity
  !> \param [in] q  specific humidity  [kg/kg]
  !> \param [in] t  temperature        [K]
  !> \param [in] p  pressure           [Pa]
  !> \return        relative humidity  [0..1]
  !-------------------------------------------
  elemental function rh_q (q, t, p) result (rh)
  real(wp) ,intent (in) :: q          ! specific humidity
  real(wp) ,intent (in) :: t          ! temperature
  real(wp) ,intent (in) :: p          ! pressure
  real(wp)              :: rh         ! relative humidity
!    e = q * p / (rdrd + emrdrd * q)
!    if (t >= b3) then
!      es = b1 * exp (b2w * (t - b3) / (t - b4w))
!    else
!      es = b1 * exp (b2e * (t - b3) / (t - b4e))
!    endif
!    es = min (es, p)    ! saveguard for low p
    if (t >= b3) then
      rh = (q*p / (rdrd + emrdrd*q)) / min(b1 * exp (b2w * (t-b3) / (t-b4w)),p)
    else
      rh = (q*p / (rdrd + emrdrd*q)) / min(b1 * exp (b2e * (t-b3) / (t-b4e)),p)
    endif
  end function rh_q
!------------------------------------------------------------------------------
  !-----------------------------------------------------
  !> relative humidity over water from specific humidity
  !> \param [in] q  specific humidity  [kg/kg]
  !> \param [in] t  temperature        [K]
  !> \param [in] p  pressure           [Pa]
  !> \return        relative humidity  [0..1]
  !-----------------------------------------------------
  elemental function rhw_q (q, t, p) result (rh)
  real(wp) ,intent (in) :: q          ! specific humidity
  real(wp) ,intent (in) :: t          ! temperature
  real(wp) ,intent (in) :: p          ! pressure
  real(wp)              :: rh         ! relative humidity
    rh = (q*p / (rdrd + emrdrd*q)) / min(b1 * exp (b2w * (t-b3) / (t-b4w)) ,p)
  end function rhw_q
!------------------------------------------------------------------------------
  !------------------------------------------------
  !> relative humidity over water from mixing ratio
  !> \param [in] m mixing ratio        [kg/kg]
  !> \param [in] t  temperature        [K]
  !> \param [in] p  pressure           [Pa]
  !> \return        relative humidity  [0..1]
  !------------------------------------------------
  elemental function rhw_m (m, t, p) result (rh)
  real(wp) ,intent (in) :: m          ! mixing ratio
  real(wp) ,intent (in) :: t          ! temperature
  real(wp) ,intent (in) :: p          ! pressure
  real(wp)              :: rh         ! relative humidity
    real(wp) :: q ! specific humidity
    q = m / (1._wp + m)
    rh = (q*p / (rdrd + emrdrd*q)) / min(b1 * exp (b2w * (t-b3) / (t-b4w)) ,p)
  end function rhw_m
!------------------------------------------------------------------------------
  !-------------------------------------------------------------
  !>  relative humidity from specific humidity, adjoint routine
  !
  !> \param [in,out] rh_a relative humidity, adjoint variable
  !> \param [in,out] q_a  specific humidity, adjoint variable
  !> \param [in,out] t_a  temperature,       adjoint variable
  !> \param [in,out] p_a  pressure,          adjoint variable
  !> \param [in]     q    specific humidity  [kg/kg]
  !> \param [in]     t    temperature        [K]
  !> \param [in]     p    pressure           [Pa]
  !>
  !> adjoints are added to arguments q_a, t_a, p_a.
  !> rh_a is set to zero on return.
  !-------------------------------------------------------------
  elemental subroutine rh_q_adj (rh_a, q_a, t_a, p_a, q, t, p)
  real(wp) ,intent (inout) :: rh_a       ! relative humidity, adjoint variable
  real(wp) ,intent (inout) :: q_a        ! specific humidity, adjoint variable
  real(wp) ,intent (inout) :: t_a        ! temperature,       adjoint variable
  real(wp) ,intent (inout) :: p_a        ! pressure,          adjoint variable
  real(wp) ,intent (in)    :: q          ! specific humidity
  real(wp) ,intent (in)    :: t          ! temperature
  real(wp) ,intent (in)    :: p          ! pressure
    real(wp) :: e, d, e_a, es_a, es, c1, c2
!
! nonlinear code:
! ES = b1 * exp (b2x * (T - b3) / (T - b4x))
! ES = min (ES, P)                                 ! saveguard for low p
! D  = rdrd + emrdrd * Q
! E  = Q*P / D
! RH = E / ES
!
! tangentlinear code
!
! es = ES * ( t * [b2x/(T-b4x) - b2x*(T-b3)/(T-b4x)**2]
! if (ES > P) ES = P
! d  = emrdrd * q
! e  = q * P/D + p * Q/D - d * Q*P/D**2
! rh = e/ES - es * E/ES**2
!
! adjoint code
!
   if (t >= b3) then
      es = min( b1 * exp (b2w * (t - b3) / (t - b4w)) ,p)
      c2 = 1._wp/ (t-b4w)
      c1 = b2w * c2
    else
      es = min( b1 * exp (b2e * (t - b3) / (t - b4e)) ,p)
      c2 = 1._wp/ (t-b4e)
      c1 = b2e * c2
    endif
    d    =   rdrd + emrdrd * q
    e    =   q*p / d
    e_a  =        rh_a / es
    es_a =      - rh_a * e/es**2
    p_a  =  p_a + e_a * Q/D
    q_a  =  q_a + e_a * (P/D - Q*P/D**2 * EMRDRD)
    if (es >= p) then
      p_a  =  p_a + es_a
    else
      t_a  =  t_a - e_a * e * c1 * (1._wp - (t-b3)*c2)
    endif
    rh_a =  0._wp

  end subroutine rh_q_adj
!------------------------------------------------------------------------------
  elemental subroutine rhw_q_adj (rh_a, q_a, t_a, p_a, q, t, p)
  real(wp) ,intent (inout) :: rh_a       ! relative humidity, adjoint variable
  real(wp) ,intent (inout) :: q_a        ! specific humidity, adjoint variable
  real(wp) ,intent (inout) :: t_a        ! temperature,       adjoint variable
  real(wp) ,intent (inout) :: p_a        ! pressure,          adjoint variable
  real(wp) ,intent (in)    :: q          ! specific humidity
  real(wp) ,intent (in)    :: t          ! temperature
  real(wp) ,intent (in)    :: p          ! pressure
    real(wp) :: e, d, e_a, es_a, es, c1, c2
!
! nonlinear code:
! ES = b1 * exp (b2x * (T - b3) / (T - b4x))
! ES = min (ES, P)                                 ! saveguard for low p
! D  = rdrd + emrdrd * Q
! E  = Q*P / D
! RH = E / ES
!
! tangentlinear code
!
! es = ES * ( t * [b2x/(T-b4x) - b2x*(T-b3)/(T-b4x)**2]
! if (ES > P) ES = P
! d  = emrdrd * q
! e  = q * P/D + p * Q/D - d * Q*P/D**2
! rh = e/ES - es * E/ES**2
!
! adjoint code
!
    es   = min( b1 * exp (b2w * (t - b3) / (t - b4w)) ,p)
    c2   = 1._wp/ (t-b4w)
    c1   = b2w * c2

    d    =   rdrd + emrdrd * q
    e    =   q*p / d
    e_a  =        rh_a / es
    es_a =      - rh_a * e/es**2
    p_a  =  p_a + e_a * Q/D
    q_a  =  q_a + e_a * (P/D - Q*P/D**2 * EMRDRD)
    if (es >= p) then
      p_a  =  p_a + es_a
    else
      t_a  =  t_a - e_a * e * c1 * (1._wp - (t-b3)*c2)
    endif
    rh_a =  0._wp

  end subroutine rhw_q_adj
!==============================================================================
  !-----------------------------------------------------
  !> relative humidity over water from specific humidity using Hardy (1998)
  !> \param [in] q  specific humidity  [kg/kg]
  !> \param [in] t  temperature        [K]
  !> \param [in] p  pressure           [Pa]
  !> \return        relative humidity  [0..1]
  !-----------------------------------------------------
  elemental function rhw_q_hardy (q, t, p) result (rh)
    real(wp) ,intent (in) :: q          ! specific humidity
    real(wp) ,intent (in) :: t          ! temperature
    real(wp) ,intent (in) :: p          ! pressure
    real(wp)              :: rh         ! relative humidity
    rh = q*p / ((rdrd + emrdrd*q) * min(esw_t_hardy(t), p))
  end function rhw_q_hardy
!------------------------------------------------------------------------------
  !------------------------------------------------
  !> relative humidity over water from mixing ratio using Hardy (1998)
  !> \param [in] m mixing ratio        [kg/kg]
  !> \param [in] t  temperature        [K]
  !> \param [in] p  pressure           [Pa]
  !> \return        relative humidity  [0..1]
  !------------------------------------------------
  elemental function rhw_m_hardy (m, t, p) result (rh)
    real(wp) ,intent (in) :: m          ! mixing ratio
    real(wp) ,intent (in) :: t          ! temperature
    real(wp) ,intent (in) :: p          ! pressure
    real(wp)              :: rh         ! relative humidity
    real(wp) :: q ! specific humidity
    q = m / (1._wp + m)
    rh = q*p / ((rdrd + emrdrd*q) * min(esw_t_hardy(t), p))
  end function rhw_m_hardy
!==============================================================================
  !-------------------------------------------
  !>  specific humidity to ppmv over dry air
  !> \param [in] q  specific humidity  [kg/kg]
  !> \return        ppmv over dry air  [0..1]
  !-------------------------------------------
  elemental function ppmv_dry_q (q) result (ppmv)
  real(wp) ,intent (in) :: q          ! specific humidity
  real(wp)              :: ppmv       ! ppmv over dry air
    ppmv = (1.d6 * q) / (RDRD * (1 - q))
  end function ppmv_dry_q
!------------------------------------------------------------------------------
  !----------------------------------------------------------
  !>  derivative d(ppmv over dry air)/d(specific humidity)
  !> \param [in] q  specific humidity  [kg/kg]
  !> \return        d(ppmv_dry)/dq
  !----------------------------------------------------------
  elemental function dppmv_dry_dq (q) result (dppmv_dq)
  real(wp) ,intent (in) :: q          ! specific humidity
  real(wp)              :: dppmv_dq   ! d(ppmv_dry)/dq
    dppmv_dq = 1.d6 / (RDRD * (1 - q)**2)
  end function dppmv_dry_dq
!------------------------------------------------------------------------------
  !-------------------------------------------
  !>  ppmv over dry air to specific humidity
  !> \param [in] q  ppmv over dry air  [0..1]
  !> \return        specific humidity  [kg/kg]
  !-------------------------------------------
  elemental function q_ppmv_dry (ppmv) result (q)
  real(wp) ,intent (in) :: ppmv       ! ppmv over dry air
  real(wp)              :: q          ! specific humidity
    q = (RDRD * ppmv) / (1.d6 + RDRD * ppmv)
  end function q_ppmv_dry
!------------------------------------------------------------------------------
  !---------------------------------------------
  !>  specific humidity to ppmv over moist air
  !> \param [in] q  specific humidity    [kg/kg]
  !> \return        ppmv over moist air  [0..1]
  !---------------------------------------------
  elemental function ppmv_moist_q (q) result (ppmv)
  real(wp) ,intent (in) :: q          ! specific humidity
  real(wp)              :: ppmv       ! ppmv over moist air
    ppmv = (1.d6 * q) / (RDRD * (1 - q) + q)
  end function ppmv_moist_q
!------------------------------------------------------------------------------
  !----------------------------------------------------------
  !>  derivative d(ppmv over moist air)/d(specific humidity)
  !> \param [in] q  specific humidity  [kg/kg]
  !> \return        d(ppmv_moist)/dq
  !----------------------------------------------------------
  elemental function dppmv_moist_dq (q) result (dppmv_dq)
  real(wp) ,intent (in) :: q          ! specific humidity
  real(wp)              :: dppmv_dq   ! d(ppmv_moist)/dq
    dppmv_dq = 1.d6 / (RDRD * (1 - q + q/RDRD)**2)
  end function dppmv_moist_dq
!------------------------------------------------------------------------------
  !-------------------------------------------
  !>  ppmv over moist air to specific humidity
  !> \param [in] q  ppmv over moist air  [0..1]
  !> \return        specific humidity  [kg/kg]
  !-------------------------------------------
  elemental function q_ppmv_moist (ppmv) result (q)
  real(wp) ,intent (in) :: ppmv       ! ppmv over moist air
  real(wp)              :: q          ! specific humidity
    q = (RDRD * ppmv) / (1.d6 + (RDRD-1.) * ppmv)
  end function q_ppmv_moist
!------------------------------------------------------------------------------
  !-------------------------------------------
  !>  specific conc. of tracegas to ppmv over dry air
  !> \param [in] trg specific tracegas conc.  [kg/kg]
  !> \return        ppmv over dry air  [0..1]
  !-------------------------------------------
  elemental function q_ppmv_dry_trg (trg,q,gas_id) result (qtrg)
    real(wp), intent(in) :: trg        ! ppmv of trace das
    real(wp), intent(in) :: q          ! specific humidity
    integer,  intent(in) :: gas_id     ! gas_id_* from rttov_const
    real(wp)             :: qtrg       ! specific conc. of trace gas
    if (gas_id > 0 .and. gas_id <= ngases) then
      qtrg = trg * (1 - q) * R_RTTOV/(RTRG(gas_id)*1.d6)
    else
      qtrg = huge(qtrg)
    end if
  end function q_ppmv_dry_trg
!------------------------------------------------------------------------------
  !-------------------------------------------
  !>  specific conc. of tracegas to ppmv over dry air
  !> \param [in] trg specific tracegas conc.  [kg/kg]
  !> \return        ppmv over dry air  [0..1]
  !-------------------------------------------
  elemental function ppmv_dry_trg (trg,q,gas_id) result (ppmv)
    real(wp), intent(in) :: trg        ! specific conc. of trace gas
    real(wp), intent(in) :: q          ! specific humidity
    integer,  intent(in) :: gas_id     ! gas_id_* from rttov_const
    real(wp)             :: ppmv       ! ppmv over dry air
    if (gas_id > 0 .and. gas_id <= ngases) then
      ppmv = (1.d6 * trg) / (1 - q) * RTRG(gas_id)/R_RTTOV
    else
      ppmv = huge(ppmv)
    end if
  end function ppmv_dry_trg
!------------------------------------------------------------------------------
  !-------------------------------------------
  !>  specific conc. of tracegas to ppmv over moist air
  !> \param [in] trg specific tracegas conc.  [kg/kg]
  !> \return        ppmv over moist air  [0..1]
  !-------------------------------------------
  elemental function ppmv_moist_trg (trg,q,gas_id) result (ppmv)
    real(wp), intent(in) :: trg        ! specific conc. of trace gas
    real(wp), intent(in) :: q          ! specific humidity
    integer,  intent(in) :: gas_id     ! gas_id_* from rttov_const
    real(wp)             :: ppmv       ! ppmv over moist air
    if (gas_id > 0 .and. gas_id <= ngases) then
      ppmv = (1.d6 * trg * RTRG(gas_id)) / (R_RTTOV - q * (R_RTTOV + RD))
    else
      ppmv = huge(ppmv)
    end if
  end function ppmv_moist_trg
!------------------------------------------------------------------------------

!==============================================================================
  !-----------------------------------------------------------------------
  !>  derive p from ps for hybrid vertical coordinated (GME), 1-d (column)
  !>
  !> \param [out] ph  pressure at half levels            [Pa]
  !> \param [out] pf  (optional) pressure at full levels [Pa]
  !> \param [in]  a   sigma coordinate coefficients
  !> \param [in]  b   sigma coordinate coefficients
  !> \param [in]  ps  surface pressure                   [Pa]
  !-----------------------------------------------------------------------
  subroutine ph_pf_ps_1 (ph, pf, a, b, ps)
  real(wp) ,intent (out)           :: ph (:) ! pressure at half levels
  real(wp) ,intent (out) ,optional :: pf (:) ! pressure at full levels
  real(wp) ,intent (in)            :: a  (:) ! sigma coordinate coeff.
  real(wp) ,intent (in)            :: b  (:) ! sigma coordinate coeff.
  real(wp) ,intent (in)            :: ps     ! surface pressure
    integer :: k, kep1
    kep1 = size (ph)
    do k = 1,kep1
      ph(k) = a(k) + b(k) * ps
    end do
    if(present(pf)) pf = 0.5_wp * (ph(:kep1-1) + ph(2:))
  end subroutine ph_pf_ps_1
!------------------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !>  derive p from ps for hybrid vertical coordinated (GME), 3-d (volume)
  !> \param [out] ph  pressure at half levels            [Pa]
  !> \param [out] pf  (optional) pressure at full levels [Pa]
  !> \param [in]  a   sigma coordinate coefficients
  !> \param [in]  b   sigma coordinate coefficients
  !> \param [in]  ps  surface pressure                   [Pa]
  !-----------------------------------------------------------------------
  subroutine ph_pf_ps_3 (ph, pf, a, b, ps)
  real(wp) ,intent (out)           :: ph (:,:,:) ! pressure at half levels
  real(wp) ,intent (out) ,optional :: pf (:,:,:) ! pressure at full levels
  real(wp) ,intent (in)            :: a  (:)     ! sigma coordinate coeff.
  real(wp) ,intent (in)            :: b  (:)     ! sigma coordinate coeff.
  real(wp) ,intent (in)            :: ps (:,:)   ! surface pressure
    integer :: k, kep1
    kep1 = size (ph,3)
    do k = 1,kep1
      ph(:,:,k) = a(k) + b(k) * ps
    end do
    if(present(pf)) pf = 0.5_wp * (ph(:,:,:kep1-1) + ph(:,:,2:))
  end subroutine ph_pf_ps_3
!------------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>  derive p from ps for hybrid vertical coordinated (GME), 3-d, adjoint
  !> \param [in,out]  ph_a  pressure at half levels,            adjoint variable
  !> \param [in,out]  pf_a  (optional) pressure at full levels, adjoint variable
  !> \param [in,out]  ps_a  surface pressure,                   adjoint variable
  !> \param [in]  b   sigma coordinate coefficients
  !-------------------------------------------------------------------------
  subroutine ph_pf_ps_adj3 (ph_a, pf_a, ps_a, b)
  real(wp) ,intent (inout)           :: ph_a (:,:,:) ! pressure at half levels
  real(wp) ,intent (inout) ,optional :: pf_a (:,:,:) ! pressure at full levels
  real(wp) ,intent (inout)           :: ps_a (:,:)   ! surface pressure
  real(wp) ,intent (in)              :: b  (:)       ! sigma coordinate coeff.
    integer  :: k, kep1
!---------------
! nonlinear code
!---------------
!   kep1 = size (ph,3)
!   do k = 1,kep1
!     ph(:,:,k) = a(k) + b(k) * ps
!   end do
!   if(present(pf)) pf = 0.5_wp * (ph(:,:,:kep1-1) + ph(:,:,2:))
!--------------------
! tangent linear code
!--------------------
!   kep1 = size (ph,3)
!   do k = 1,kep1
!     ph_t(:,:,k) = b(k) * ps_t
!   end do
!   if(present(pf_t)) pf_t = 0.5_wp * (ph_t(:,:,:kep1-1) + ph_t(:,:,2:))
!-------------
! adjoint code
!-------------
    kep1 = size (ph_a,3)
    if(present(pf_a)) then
      ph_a(:,:,:kep1-1) = ph_a(:,:,:kep1-1) + 0.5_wp * pf_a
      ph_a(:,:,2:)      = ph_a(:,:,2:)      + 0.5_wp * pf_a
      pf_a              = 0._wp
    endif
    do k = 1,kep1
      ps_a        = ps_a   + ph_a(:,:,k) * b(k)
      ph_a(:,:,k) = 0._wp
    end do
  end subroutine ph_pf_ps_adj3
!==============================================================================
  !------------------------------------------------------
  !>  derive t from tv, q
  !> \param [in]  tv virtual temperature        [K]
  !> \param [in]  q  specific humidity          [kg/kg]
  !> \param [in]  xl cloud liquid water content [kg/m^3]
  !> \param [in]  xi cloud ice content          [kg/m^3]
  !> \return         temperature                [K]
  !-----------------------------------------------------
  elemental function t_tv_q (tv, q, xl, xi) result (t)
    real(wp)                        :: t
    real(wp) ,intent (in)           :: tv
    real(wp) ,intent (in)           :: q
    real(wp) ,intent (in) ,optional :: xl
    real(wp) ,intent (in) ,optional :: xi
    if (present(xl) .and. present(xi)) then
      t = tv / (1.0_wp + RDDRM1 * q - (xl + xi))
    else if (present(xl)) then
      t = tv / (1.0_wp + RDDRM1 * q -  xl)
    else
      t = tv / (1.0_wp + RDDRM1 * q)
    endif
  end function t_tv_q
!------------------------------------------------------------------------------
  !------------------------------------------------------------
  !>  derive t from tv, q; adjoint routine
  !> \param [in,out] t_a  temperature         ;adjoint variable
  !> \param [in,out] tv_a virtual temperature ;adjoint variable
  !> \param [in,out] q_a  specific humidity   ;adjoint variable
  !> \param [in,out] xl_a cloud liquid water  ;adjoint variable
  !> \param [in,out] xi_a cloud ice           ;adjoint variable
  !> \param [in]     tv   virtual temperature [K]
  !> \param [in]     q    specific humidity   [kg/kg]
  !> \param [in]     xl   cloud liquid water  [kg/m^3]
  !> \param [in]     xi   cloud ice           [kg/m^3]
  !------------------------------------------------------------
  elemental subroutine  t_tv_q_adj (t_a, tv_a, q_a, xl_a, xi_a, tv, q, xl, xi)
    real(wp) ,intent (inout)           :: t_a
    real(wp) ,intent (inout)           :: tv_a
    real(wp) ,intent (inout)           :: q_a
    real(wp) ,intent (inout) ,optional :: xl_a
    real(wp) ,intent (inout) ,optional :: xi_a
    real(wp) ,intent (in)              :: tv
    real(wp) ,intent (in)              :: q
    real(wp) ,intent (in)    ,optional :: xl
    real(wp) ,intent (in)    ,optional :: xi
!------------------------------
! nonlinear code
!------------------------------
!   t = tv / (1.+RDDRM1* q - (xl +xi))
!----------------------------
! adjoint code:
!----------------------------
  real(wp) :: x, y

  x = 1.0_wp + RDDRM1 * q
  if (present(xl)) x = x - xl
  if (present(xi)) x = x - xi
  y = tv / x**2

  tv_a                    = tv_a + t_a          / x
  q_a                     = q_a  - t_a * RDDRM1 * y
  if (present(xl_a)) xl_a = xl_a +                y
  if (present(xi_a)) xi_a = xi_a +                y
  t_a                     = 0._wp

  end subroutine  t_tv_q_adj
!==============================================================================
  !-------------------------------------------------------------------
  !> Calculate t, q from tv, rh
  !>
  !> VERSION WITHOUT ADJOINT:
  !> NEC SX cannot handle optional intent(out) in elemental subroutine
  !>
  !> \param [out] t  temperature          [K]
  !> \param [out] q  specific humidity    [kg/kg]
  !> \param [in]  tv virtual temperature  [K]
  !> \param [in]  rh relative humidity    [0..1]
  !> \param [in]  p  pressure             [Pa]
  !> \param [in]  xl liquid water content [kg/m^3]
  !> \param [in]  xi ice content          [kg/m^3]
  !> \param [out] i_fail >0: number of iterations; <0: no convergence
  !> \param [out] x_fail debug information
  !-------------------------------------------------------------------
  elemental subroutine tq_tvrh_noadj (t, q, tv, rh, p, xl, xi, i_fail, x_fail)
  real(wp) ,intent (out)            :: t      ! temperature
  real(wp) ,intent (out)            :: q      ! spec. humidity
  real(wp) ,intent (in)             :: tv     ! virt. temperature
  real(wp) ,intent (in)             :: rh     ! rel. humidity
  real(wp) ,intent (in)             :: p      ! pressure
  real(wp) ,intent (in)  ,optional  :: xl     ! liquid water content
  real(wp) ,intent (in)  ,optional  :: xi     ! ice          content
  integer  ,intent (out)            :: i_fail ! >0: number of iterations
                                              ! <0: no convergence
   type(t_tq_tvrh)                   &
            ,intent (out)           :: x_fail ! debug information
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! VERSION WITHOUT ADJOINT
  ! NEC SX cannot handle optional intent(out) in elemental subroutine
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    call tq_tvrh (t, q, tv, rh, p, xl, xi, i_fail=i_fail, x_fail=x_fail)
  end subroutine tq_tvrh_noadj
!------------------------------------------------------------------------------
  !-------------------------------------------------------------------
  !> Calculate t, q from tv, rh
  !>
  !> VERSION WITHOUT ADJOINT and x_fail:
  !> NEC SX cannot handle optional intent(out) in elemental subroutine
  !>
  !> \param [out] t  temperature          [K]
  !> \param [out] q  specific humidity    [kg/kg]
  !> \param [in]  tv virtual temperature  [K]
  !> \param [in]  rh relative humidity    [0..1]
  !> \param [in]  p  pressure             [Pa]
  !> \param [in]  xl liquid water content [kg/m^3]
  !> \param [in]  xi ice content          [kg/m^3]
  !> \param [out] i_fail >0: number of iterations; <0: no convergence
  !-------------------------------------------------------------------
  elemental subroutine tq_tvrh_noxfail (t, q, tv, rh, p, xl, xi, i_fail)
  real(wp) ,intent (out)            :: t      ! temperature
  real(wp) ,intent (out)            :: q      ! spec. humidity
  real(wp) ,intent (in)             :: tv     ! virt. temperature
  real(wp) ,intent (in)             :: rh     ! rel. humidity
  real(wp) ,intent (in)             :: p      ! pressure
  real(wp) ,intent (in)  ,optional  :: xl     ! liquid water content
  real(wp) ,intent (in)  ,optional  :: xi     ! ice    water content
  integer  ,intent (out)            :: i_fail ! >0: number of iterations
                                              ! <0: no convergence
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! VERSION WITHOUT ADJOINT and x_fail
  ! NEC SX cannot handle optional intent(out) in elemental subroutine
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    call tq_tvrh (t, q, tv, rh, p, xl, xi, i_fail=i_fail)
  end subroutine tq_tvrh_noxfail
!------------------------------------------------------------------------------
  !---------------------------------------------------------------------
  !> Calculate t, q from tv, rh; with partial derivatives
  !>
  !> Transform virtual temperature (tv) and relative humidity (rh)
  !>                to temperature (t)  and specific humidity (q).
  !>
  !> Provide the Jacobi-matrix elements as well.
  !>
  !> The result is derived iteratively from the inverse relationships by
  !> the Newton method.
  !>
  !> \param [out] t      temperature          [K]
  !> \param [out] q      specific humidity    [kg/kg]
  !> \param [in]  tv     virtual temperature  [K]
  !> \param [in]  rh     relative humidity    [0..1]
  !> \param [in]  p      pressure             [Pa]
  !> \param [in]  xl     liquid water content [kg/m^3]
  !> \param [in]  xi     ice content          [kg/m^3]
  !> \param [out] i_fail >0: number of iterations; <0: no convergence
  !> \param [out] dt_tv  partial derivative: d t / d tv
  !> \param [out] dt_rh  partial derivative: d t / d rh
  !> \param [out] dq_tv  partial derivative: d q / d tv
  !> \param [out] dq_rh  partial derivative: d q / d rh
  !> \param [out] x_fail debug information
  !-------------------------------------------------------------------
  elemental subroutine tq_tvrh (t, q, tv, rh, p, xl, xi, lw, i_fail, &
                                dt_tv, dt_rh, dq_tv, dq_rh,          &
                                x_fail                               )
  real(wp) ,intent (out)            :: t      ! temperature
  real(wp) ,intent (out)            :: q      ! spec. humidity
  real(wp) ,intent (in)             :: tv     ! virt. temperature
  real(wp) ,intent (in)             :: rh     ! rel. humidity
  real(wp) ,intent (in)             :: p      ! pressure
  real(wp) ,intent (in)  ,optional  :: xl     ! liquid water content
  real(wp) ,intent (in)  ,optional  :: xi     ! ice    water content
  logical  ,intent (in)  ,optional  :: lw     ! rel. humidity over water
  integer  ,intent (out)            :: i_fail ! >0: number of iterations
                                              ! <0: no convergence
  real(wp) ,intent (out) ,optional  :: dt_tv  ! part. deriv. dt/dtv
  real(wp) ,intent (out) ,optional  :: dt_rh  ! part. deriv. dt/drh
  real(wp) ,intent (out) ,optional  :: dq_tv  ! part. deriv. dq/dtv
  real(wp) ,intent (out) ,optional  :: dq_rh  ! part. deriv. dq/drh
  type(t_tq_tvrh)                   &
           ,intent (out) ,optional  :: x_fail ! debug information

    real(wp) :: tv_q   ! part. deriv. dtv/dq
    real(wp) :: tv_t   ! part. deriv. dtv/dt
    real(wp) :: rh_t   ! part. deriv. drh/dt
    real(wp) :: rh_q_  ! part. deriv. dtrh/dq
    real(wp) :: det    ! determinant of Jacobi-matrix
    real(wp) :: qmin,qmax,tmin,tmax
    real(wp) :: p_a
    real(wp) :: tv_a
    real(wp) :: rh_a
    real(wp) :: tv_n
    real(wp) :: rh_n
    integer  :: i
    real(wp) :: eps_tv, eps_rh
    real(wp) :: del_tv
    real(wp) :: del_rh
    real(wp) :: x
    logical  :: llw

    llw = .false.; if (present(lw)) llw = lw
    !-------------------------------------------------
    ! accuracy demanded for the recalculated arguments
    !-------------------------------------------------
    eps_tv = 100._wp * spacing (tv)
    eps_rh = 100._wp * spacing (rh)
    !------------
    ! first guess
    !------------
    t = tv
    if (llw) then
      q = q_rhw (rh, t, p)
    else
      q = q_rh  (rh, t, p)
    endif
    !------------------------
    ! physical bounds on q, t
    !------------------------
    x = 0._wp
    if (present(xl)) x = x + xl
    if (present(xi)) x = x + xi
    qmin = 0._wp
    qmax = 1._wp
    tmin = tv * RdRd
    tmax = tv / (1._wp - x)
    if (present(x_fail)) then
      x_fail% qmax = qmax
      x_fail% qmin = qmin
      x_fail% tmax = tmax
      x_fail% tmin = tmin
    endif

    !--------------------------------
    ! iterations of the Newton method
    !--------------------------------
    do i = 1, maxiter_tq_tvrh

      q = max (q,qmin)
      q = min (q,qmax)
      t = max (t,tmin)
      t = min (t,tmax)

      tv_n = tv_t_q (t,  q, xl=xl, xi=xi)
      rh_n = rh_q   (q,  t, p)

      del_tv = tv_n - tv
      del_rh = rh_n - rh

      tv_a = 1._wp
      tv_t = 0._wp
      tv_q = 0._wp
      call tv_t_q_adj (tv_a, tv_t, tv_q, t=t, q=q, xl=xl, xi=xi)

      rh_a  = 1._wp
      rh_q_ = 0._wp
      rh_t  = 0._wp
      p_a   = 0._wp

      if (llw) then
        call rhw_q_adj (rh_a, rh_q_, rh_t, p_a, q, t, p)
      else
        call rh_q_adj  (rh_a, rh_q_, rh_t, p_a, q, t, p)
      endif

      call rh_q_adj (rh_a, rh_q_, rh_t, p_a, q, t, p)

      if (present (x_fail)) then
        x_fail % t (i) = t
        x_fail % q (i) = q
        x_fail % tv(i) = tv_n
        x_fail % rh(i) = rh_n
      endif

      det  = tv_t * rh_q_ - rh_t * tv_q
      t = t - (   rh_q_ * del_tv - tv_q  * del_rh) / det
      q = q - ( - rh_t  * del_tv + tv_t  * del_rh) / det
      if (q < 0._wp ) q = 0._wp

      i_fail = i
      if ((abs(del_tv) <= eps_tv).and.(abs(del_rh) <= eps_rh))   exit
      i_fail = -i

    enddo
    !-------------------------------------------------------------
    ! Derive Jacobi-matrix elements from the inverse Jacobi-matrix
    !-------------------------------------------------------------
    if (present(dt_tv)) dt_tv  =   rh_q_ / det
    if (present(dt_rh)) dt_rh  = - tv_q  / det
    if (present(dq_tv)) dq_tv  = - rh_t  / det
    if (present(dq_rh)) dq_rh  =   tv_t  / det

  end subroutine tq_tvrh
!------------------------------------------------------------------------------
  ! Vectorized versions of tq_tvrh:
  ! - for SX-9, using indirect addressing:    tq_tvrh_vec
  ! Vectorized versions of tq_tvrh:
  ! - for SX Aurora, using masked operations: tq_tvrh_vec2
  !-------------------------------------------------------------------
  !> Calculate t, q from tv, rh; with partial derivatives
  !>
  !> Vectorization-adapted version of tq_tvrh
  !>
  !> \param [in]  n      number of gridpoints
  !> \param [out] t      temperature          [K]
  !> \param [out] q      specific humidity    [kg/kg]
  !> \param [in]  tv     virtual temperature  [K]
  !> \param [in]  rh     relative humidity    [0..1]
  !> \param [in]  p      pressure             [Pa]
  !> \param [in]  xl     liquid water content [kg/m^3]
  !> \param [out] i_fail >0: number of iterations; <0: no convergence
  !> \param [out] dt_tv  partial derivative: d t / d tv
  !> \param [out] dt_rh  partial derivative: d t / d rh
  !> \param [out] dq_tv  partial derivative: d q / d tv
  !> \param [out] dq_rh  partial derivative: d q / d rh
  !> \param [out] x_fail debug information
  !-------------------------------------------------------------------
  subroutine tq_tvrh_vec (n, t, q, tv, rh, p, xl, i_fail, &
                             dt_tv, dt_rh, dq_tv, dq_rh,  &
                             x_fail                       )
  integer  ,intent(in)             :: n         ! Number of gridpoints
  real(wp) ,intent(out)            :: t     (n) ! temperature
  real(wp) ,intent(out)            :: q     (n) ! spec. humidity
  real(wp) ,intent(in)             :: tv    (n) ! virt. temperature
  real(wp) ,intent(in)             :: rh    (n) ! rel. humidity
  real(wp) ,intent(in)             :: p     (n) ! pressure
  real(wp) ,intent(in)  ,optional  :: xl    (n) ! liquid+ice water content
  integer  ,intent(out)            :: i_fail(n) ! >0: number of iterations
                                                ! <0: no convergence
  real(wp) ,intent(out) ,optional  :: dt_tv (n) ! part. deriv. dt/dtv
  real(wp) ,intent(out) ,optional  :: dt_rh (n) ! part. deriv. dt/drh
  real(wp) ,intent(out) ,optional  :: dq_tv (n) ! part. deriv. dq/dtv
  real(wp) ,intent(out) ,optional  :: dq_rh (n) ! part. deriv. dq/drh
  type(t_tq_tvrh)                  &
           ,intent(out) ,optional  :: x_fail(n) ! debug information

    real(wp) :: tv_q (n)        ! part. deriv. dtv/dq
    real(wp) :: tv_t (n)        ! part. deriv. dtv/dt
    real(wp) :: rh_t (n)        ! part. deriv. drh/dt
    real(wp) :: rh_q_(n)        ! part. deriv. dtrh/dq
    real(wp) :: det  (n)        ! determinant of Jacobi-matrix
    real(wp), dimension(n) ::  &! Work arrays:
                tmin, tmax,    &
                tv_n, rh_n,    &
                eps_tv, eps_rh
    real(wp) :: qmin, qmax, p_a, tv_a, rh_a, del_tv, del_rh
    logical  :: mask (n)        ! Mask for active gridpoints
    integer  :: idx  (n)        ! Index array of active gridpoints
    integer  :: i               ! Iteration
    integer  :: j, k            ! Loop index
    integer  :: nn              ! Count of active gridpoints
    real(wp) :: x, y            ! Workaround for sxf90 problems

    !-------------------------------------------------
    ! accuracy demanded for the recalculated arguments
    !-------------------------------------------------
    ! Work around a problem with sxf90 rev.360 (avoid fusion of assignments)
!NEC$ novector
    eps_tv = 100._wp * spacing (tv)
!NEC$ novector
    eps_rh = 100._wp * spacing (rh)
    !------------
    ! first guess
    !------------
    t = tv
    q = q_rh (rh, t, p)
    !------------------------
    ! physical bounds on q, t
    !------------------------
    qmin = 0._wp
    qmax = 1._wp
    tmin = tv * RdRd
    if (present (xl)) then
       tmax = tv / (1._wp - xl)
    else
       tmax = tv
    end if

    if (present (x_fail)) then
       x_fail% qmax = qmax
       x_fail% qmin = qmin
       x_fail% tmax = tmax
       x_fail% tmin = tmin
    endif

    mask   = .true.                     ! Set all gridpoints as active,
    i_fail = -maxiter_tq_tvrh           ! not yet converged
    !--------------------------------
    ! iterations of the Newton method
    !--------------------------------
newton: &
    do i = 1, maxiter_tq_tvrh
      nn = 0
      do k = 1, n
         if (mask(k)) then
            nn = nn+1
            idx(nn) = k
         end if
      end do
      if (nn == 0) exit newton

!$omp parallel
      if (present (xl)) then
!NEC$ nomove
!NEC$ ivdep
!$omp do private (j,k,tv_a,rh_a,p_a) schedule(static)
         do j = 1, nn
            k = idx(j)
            q(k) = max (q(k),qmin)
            q(k) = min (q(k),qmax)
            t(k) = max (t(k),tmin(k))
            t(k) = min (t(k),tmax(k))

            tv_n(k) = tv_t_q (t(k), q(k), xl=xl(k))
            rh_n(k) = rh_q   (q(k), t(k), p(k))

            tv_a    = 1._wp
            tv_t(k) = 0._wp
            tv_q(k) = 0._wp
            call tv_t_q_adj (tv_a, tv_t(k), tv_q(k), t=t(k), q=q(k), xl=xl(k))

            rh_a     = 1._wp
            rh_q_(k) = 0._wp
            rh_t (k) = 0._wp
            p_a      = 0._wp
            call rh_q_adj (rh_a, rh_q_(k), rh_t(k), p_a, q(k), t(k), p(k))
         end do
!$omp end do nowait
      else
!NEC$ ivdep
!$omp do private (j,k,tv_a,rh_a,p_a,x,y) schedule(static)
         do j = 1, nn
            k = idx(j)
            q(k) = max (q(k),qmin)
            q(k) = min (q(k),qmax)
            t(k) = max (t(k),tmin(k))
            t(k) = min (t(k),tmax(k))

#if 0
            tv_n(k) = tv_t_q (t(k), q(k))
#else
            ! Manually inlined for better optimization:
            tv_n(k) = t(k) * (1.0_wp + RDDRM1 * q(k))
#endif
            rh_n(k) = rh_q   (q(k), t(k), p(k))

            tv_a    = 1._wp
            tv_t(k) = 0._wp
            tv_q(k) = 0._wp
#if 0
            call tv_t_q_adj (tv_a, tv_t(k), tv_q(k), t=t(k), q=q(k))
#else
            ! Manually inlined for better optimization:
            x       = 1.0_wp  + RDDRM1 * q(k)
            y       =             t(k) * tv_a
            tv_t(k) = tv_t(k) + x      * tv_a
            tv_q(k) = tv_q(k) + RDDRM1 * y
            tv_a    = 0._wp
#endif

            rh_a     = 1._wp
            rh_q_(k) = 0._wp
            rh_t (k) = 0._wp
            p_a      = 0._wp
            call rh_q_adj (rh_a, rh_q_(k), rh_t(k), p_a, q(k), t(k), p(k))
         end do
!$omp end do nowait
      end if

      if (present (x_fail)) then
!NEC$ ivdep
!$omp do private (j,k) schedule(static)
         do j = 1, nn
            k = idx(j)
            x_fail(k) % t (i) = t(k)
            x_fail(k) % q (i) = q(k)
            x_fail(k) % tv(i) = tv_n(k)
            x_fail(k) % rh(i) = rh_n(k)
         end do
!$omp end do nowait
      endif

!NEC$ ivdep
!$omp do private (j,k,del_tv,del_rh) schedule(static)
      do j = 1, nn
         k = idx(j)
         del_tv = tv_n(k) - tv(k)
         del_rh = rh_n(k) - rh(k)

         det(k) = tv_t(k) * rh_q_(k) - rh_t(k) * tv_q(k)
         t  (k) = t(k) - (rh_q_(k) * del_tv - tv_q(k) * del_rh) / det(k)
         q  (k) = q(k) + (rh_t (k) * del_tv - tv_t(k) * del_rh) / det(k)
         if (q(k) < 0._wp) q(k) = 0._wp

         if ( (abs(del_tv) <= eps_tv(k)) .and. &
              (abs(del_rh) <= eps_rh(k))       ) then
            ! We're converged; remove this grid point from active set
            mask(k)   = .false.
            i_fail(k) = i
         end if
      end do
!$omp end do
!$omp end parallel

    enddo newton
    !-------------------------------------------------------------
    ! Derive Jacobi-matrix elements from the inverse Jacobi-matrix
    !-------------------------------------------------------------
    if (present(dt_tv)) dt_tv  =   rh_q_ / det
    if (present(dt_rh)) dt_rh  = - tv_q  / det
    if (present(dq_tv)) dq_tv  = - rh_t  / det
    if (present(dq_rh)) dq_rh  =   tv_t  / det

  end subroutine tq_tvrh_vec

  subroutine tq_tvrh_vec2 (n, t, q, tv, rh, p, xl, i_fail, &
                              dt_tv, dt_rh, dq_tv, dq_rh,  &
                              x_fail                       )
  integer  ,intent(in)             :: n         ! Number of gridpoints
  real(wp) ,intent(out)            :: t     (n) ! temperature
  real(wp) ,intent(out)            :: q     (n) ! spec. humidity
  real(wp) ,intent(in)             :: tv    (n) ! virt. temperature
  real(wp) ,intent(in)             :: rh    (n) ! rel. humidity
  real(wp) ,intent(in)             :: p     (n) ! pressure
  real(wp) ,intent(in)  ,optional  :: xl    (n) ! liquid+ice water content
  integer  ,intent(out)            :: i_fail(n) ! >0: number of iterations
                                                ! <0: no convergence
  real(wp) ,intent(out) ,optional  :: dt_tv (n) ! part. deriv. dt/dtv
  real(wp) ,intent(out) ,optional  :: dt_rh (n) ! part. deriv. dt/drh
  real(wp) ,intent(out) ,optional  :: dq_tv (n) ! part. deriv. dq/dtv
  real(wp) ,intent(out) ,optional  :: dq_rh (n) ! part. deriv. dq/drh
  type(t_tq_tvrh)                  &
           ,intent(out) ,optional  :: x_fail(n) ! debug information

    real(wp) :: tv_q (n)        ! part. deriv. dtv/dq
    real(wp) :: tv_t (n)        ! part. deriv. dtv/dt
    real(wp) :: rh_t (n)        ! part. deriv. drh/dt
    real(wp) :: rh_q_(n)        ! part. deriv. dtrh/dq
    real(wp) :: det  (n)        ! determinant of Jacobi-matrix
    real(wp), dimension(n) ::  &! Work arrays:
                tmin, tmax,    &
                tv_n, rh_n,    &
                eps_tv, eps_rh
    real(wp) :: qmin, qmax, p_a, tv_a, rh_a, del_tv, del_rh
    logical  :: mask (n)        ! Mask for active gridpoints
    integer  :: i               ! Iteration
    integer  :: j               ! Loop index
    real(wp) :: x, y            ! Workaround for sxf90 problems

    !-------------------------------------------------
    ! accuracy demanded for the recalculated arguments
    !-------------------------------------------------
    ! Work around a problem with sxf90 rev.360 (avoid fusion of assignments)
!NEC$ novector
    eps_tv = 100._wp * spacing (tv)
!NEC$ novector
    eps_rh = 100._wp * spacing (rh)
    !------------
    ! first guess
    !------------
    t = tv
    q = q_rh (rh, t, p)
    !------------------------
    ! physical bounds on q, t
    !------------------------
    qmin = 0._wp
    qmax = 1._wp
    tmin = tv * RdRd
    if (present (xl)) then
       tmax = tv / (1._wp - xl)
    else
       tmax = tv
    end if

    if (present (x_fail)) then
       x_fail% qmax = qmax
       x_fail% qmin = qmin
       x_fail% tmax = tmax
       x_fail% tmin = tmin
    endif

    mask   = .true.                     ! Set all gridpoints as active,
    i_fail = -maxiter_tq_tvrh           ! not yet converged
    !--------------------------------
    ! iterations of the Newton method
    !--------------------------------
newton: &
    do i = 1, maxiter_tq_tvrh

!$omp parallel
      if (present (xl)) then
!NEC$ nomove
!NEC$ ivdep
!$omp do private (j,tv_a,rh_a,p_a) schedule(static)
         do j = 1,n
            if (mask(j)) then
              q(j) = max (q(j),qmin)
              q(j) = min (q(j),qmax)
              t(j) = max (t(j),tmin(j))
              t(j) = min (t(j),tmax(j))

              tv_n(j) = tv_t_q (t(j), q(j), xl=xl(j))
              rh_n(j) = rh_q   (q(j), t(j), p(j))

              tv_a    = 1._wp
              tv_t(j) = 0._wp
              tv_q(j) = 0._wp
              call tv_t_q_adj (tv_a, tv_t(j), tv_q(j), t=t(j), q=q(j), xl=xl(j))

              rh_a     = 1._wp
              rh_q_(j) = 0._wp
              rh_t (j) = 0._wp
              p_a      = 0._wp
              call rh_q_adj (rh_a, rh_q_(j), rh_t(j), p_a, q(j), t(j), p(j))
            end if
         end do
!$omp end do nowait
      else
!NEC$ ivdep
!$omp do private (j,tv_a,rh_a,p_a,x,y) schedule(static)
         do j = 1,n
            if (mask(j)) then
              q(j) = max (q(j),qmin)
              q(j) = min (q(j),qmax)
              t(j) = max (t(j),tmin(j))
              t(j) = min (t(j),tmax(j))

#if 0
              tv_n(j) = tv_t_q (t(j), q(j))
#else
              ! Manually inlined for better optimization:
              tv_n(j) = t(j) * (1.0_wp + RDDRM1 * q(j))
#endif
              rh_n(j) = rh_q   (q(j), t(j), p(j))

              tv_a    = 1._wp
              tv_t(j) = 0._wp
              tv_q(j) = 0._wp
#if 0
              call tv_t_q_adj (tv_a, tv_t(j), tv_q(j), t=t(j), q=q(j))
#else
              ! Manually inlined for better optimization:
              x       = 1.0_wp  + RDDRM1 * q(j)
              y       =             t(j) * tv_a
              tv_t(j) = tv_t(j) + x      * tv_a
              tv_q(j) = tv_q(j) + RDDRM1 * y
              tv_a    = 0._wp
#endif

              rh_a     = 1._wp
              rh_q_(j) = 0._wp
              rh_t (j) = 0._wp
              p_a      = 0._wp
              call rh_q_adj (rh_a, rh_q_(j), rh_t(j), p_a, q(j), t(j), p(j))
            end if
         end do
!$omp end do nowait
      end if

      if (present (x_fail)) then
!NEC$ ivdep
!$omp do private (j) schedule(static)
         do j = 1,n
            if (mask(j)) then
              x_fail(j) % t (i) = t(j)
              x_fail(j) % q (i) = q(j)
              x_fail(j) % tv(i) = tv_n(j)
              x_fail(j) % rh(i) = rh_n(j)
            end if
         end do
!$omp end do nowait
      endif

!NEC$ ivdep
!$omp do private (j,del_tv,del_rh) schedule(static)
      do j = 1,n
         if (mask(j)) then
           del_tv = tv_n(j) - tv(j)
           del_rh = rh_n(j) - rh(j)

           det(j) = tv_t(j) * rh_q_(j) - rh_t(j) * tv_q(j)
           t  (j) = t(j) - (rh_q_(j) * del_tv - tv_q(j) * del_rh) / det(j)
           q  (j) = q(j) + (rh_t (j) * del_tv - tv_t(j) * del_rh) / det(j)
           if (q(j) < 0._wp) q(j) = 0._wp

           if ( (abs(del_tv) <= eps_tv(j)) .and. &
                (abs(del_rh) <= eps_rh(j))       ) then
             ! We're converged; remove this grid point from active set
             mask(j)   = .false.
             i_fail(j) = i
           end if
         end if
      end do

!$omp end do
!$omp end parallel

      if (count(mask) == 0) exit

    enddo newton
    !-------------------------------------------------------------
    ! Derive Jacobi-matrix elements from the inverse Jacobi-matrix
    !-------------------------------------------------------------
    if (present(dt_tv)) dt_tv  =   rh_q_ / det
    if (present(dt_rh)) dt_rh  = - tv_q  / det
    if (present(dq_tv)) dq_tv  = - rh_t  / det
    if (present(dq_rh)) dq_rh  =   tv_t  / det

  end subroutine tq_tvrh_vec2

!==============================================================================
  !-----------------------------------------------
  !> derive tv from t, q
  !> \param  [in] t   temperature         [K]
  !> \param  [in] q   specific humidity   [kg/kg]
  !> \param  [in] xl  cloud water content [kg/m^3]
  !> \param  [in] xi  cloud ice content   [kg/m^3]
  !> \return          virtual temperature [K]
  !-----------------------------------------------
  elemental function tv_t_q (t, q, xl, xi) result (tv)
    real(wp)                        :: tv
    real(wp) ,intent (in)           :: t
    real(wp) ,intent (in)           :: q
    real(wp) ,intent (in) ,optional :: xl
    real(wp) ,intent (in) ,optional :: xi
    if (present(xl) .and. present(xi)) then
      tv = t * (1.0_wp + RDDRM1 * q - (xl + xi))
    else if (present(xl)) then
      tv = t * (1.0_wp + RDDRM1 * q -  xl)
    else
      tv = t * (1.0_wp + RDDRM1 * q)
    endif
  end function tv_t_q
!------------------------------------------------------------------------------
  !-------------------------------------------------------------
  !> derive tv from t, q; adjoint routine
  !> \param [in,out] tv_a  virtual temperature, adjoint variable
  !> \param [in,out] t_a   temperature,         adjoint variable
  !> \param [in,out] q_a   specific humidity,   adjoint variable
  !> \param [in,out] xl_a  cloud water content, adjoint variable
  !> \param [in,out] xi_a  cloud ice content,   adjoint variable
  !> \param [in]     t     temperature          [K]
  !> \param [in]     q     specific humidity    [kg/kg]
  !> \param [in]     xl    cloud water content  [kg/m^3]
  !> \param [in]     xi    cloud ice content    [kg/m^3]
  !-------------------------------------------------------------
  elemental subroutine tv_t_q_adj (tv_a, t_a, q_a, xl_a, xi_a, t, q, xl, xi)
    real(wp) ,intent (inout)           :: tv_a
    real(wp) ,intent (inout)           :: t_a
    real(wp) ,intent (inout)           :: q_a
    real(wp) ,intent (inout) ,optional :: xl_a
    real(wp) ,intent (inout) ,optional :: xi_a
    real(wp) ,intent (in)              :: t
    real(wp) ,intent (in)              :: q
    real(wp) ,intent (in)    ,optional :: xl
    real(wp) ,intent (in)    ,optional :: xi
!---------------
! nonlinear code
!---------------
! tv = t * (1.+RDDRM1* q - (xl +xi))
!--------------------
! tangent linear code
!--------------------
! x = 1.+RDDRM1* q - (xl +xi)
! tv_a = x * t_a + t * ( RDDRM1* q_a - (xl_a +xi_a))
!-------------
! adjoint code
!-------------
  real(wp) :: x, y

  x = 1.0_wp + RDDRM1 * q
  y = t * tv_a
  if (present(xl)) x = x - xl
  if (present(xi)) x = x - xi

  t_a                     = t_a  + x      * tv_a
  q_a                     = q_a  + RDDRM1 * y
  if (present(xl_a)) xl_a = xl_a -          y
  if (present(xi_a)) xi_a = xi_a -          y
  tv_a                    = 0._wp

  end subroutine tv_t_q_adj
!==============================================================================
  !---------------------------------------------
  !> derive tv from t, rh
  !> \param [in] t  temperature         [K]
  !> \param [in] rh relative humidity   [0..1]
  !> \param [in] p  pressure            [Pa]
  !> \param [in] xl cloud water content [kg/m^3]
  !> \param [in] xi cloud ice content   [kg/m^3]
  !> \return        virtual temperature [K]
  !---------------------------------------------
  elemental function tv_t_rh (t, rh, p, xl, xi) result (tv)
    real(wp)                        :: tv
    real(wp) ,intent (in)           :: t
    real(wp) ,intent (in)           :: rh
    real(wp) ,intent (in)           :: p
    real(wp) ,intent (in) ,optional :: xl
    real(wp) ,intent (in) ,optional :: xi

      real(wp) :: q
      q  = q_rh   (rh, t, p)
      tv = tv_t_q (t, q, xl, xi)

  end function tv_t_rh
!------------------------------------------------------------------------------
  !-------------------------------------------------------------
  !> derive tv from t, rh; adjoint routine
  !> \param [in,out]  tv_a virtualtemperature,  adjoint variable
  !> \param [in,out]  t_a  temperature,         adjoint variable
  !> \param [in,out]  rh_a relative humidity,   adjoint variable
  !> \param [in,out]  p_a  pressure,            adjoint variable
  !> \param [in,out]  xl_a cloud water content, adjoint variable
  !> \param [in,out]  xi_a cloud ice content,   adjoint variable
  !> \param [in]      t    temperature          [K]
  !> \param [in]      rh   relative humidity    [0..1]
  !> \param [in]      p    pressure             [Pa]
  !> \param [in]      xl   cloud water content  [kg/m^3]
  !> \param [in]      xi   cloud ice content    [kg/m^3]
  !-------------------------------------------------------------

  elemental subroutine tv_t_rh_adj (tv_a, t_a, rh_a, p_a, xl_a, xi_a, &
                                          t,   rh,   p,   xl,   xi)
    real(wp) ,intent (inout)           :: tv_a
    real(wp) ,intent (inout)           :: t_a
    real(wp) ,intent (inout)           :: rh_a
    real(wp) ,intent (inout)           :: p_a
    real(wp) ,intent (inout) ,optional :: xl_a
    real(wp) ,intent (inout) ,optional :: xi_a
    real(wp) ,intent (in)              :: t
    real(wp) ,intent (in)              :: rh
    real(wp) ,intent (in)              :: p
    real(wp) ,intent (in)    ,optional :: xl
    real(wp) ,intent (in)    ,optional :: xi

    real(wp) :: q, q_a

    q  = q_rh   (rh, t, p)
    call tv_t_q_adj (tv_a, t_a, q_a, xl_a, xi_a, t, q, xl, xi)
    call q_rh_adj   (q_a, rh_a, t_a, p_a, rh, t, p)

  end subroutine tv_t_rh_adj
!==============================================================================
!  elemental subroutine uv_fd (u, v, f, d, du, dv, df, dd)
!  !----------------------------------------
!  ! convert wind speed and direction to u,v
!  !----------------------------------------
!  real(wp) ,intent(out)           :: u  ! wind component u       [m/s]
!  real(wp) ,intent(out)           :: v  ! wind component v       [m/s]
!  real(wp) ,intent(in)            :: f  ! wind speed             [m/s]
!  real(wp) ,intent(in)            :: d  ! wind direction         [degree]
!  real(wp) ,intent(out) ,optional :: du ! wind component u error [m/s]
!  real(wp) ,intent(out) ,optional :: dv ! wind component v error [m/s]
!  real(wp) ,intent(in)  ,optional :: df ! wind speed       error [m/s]
!  real(wp) ,intent(in)  ,optional :: dd ! wind direction   error [degree]
!    real(wp) :: eu, ev, r
!    !---------------
!    ! nonlinear code
!    !---------------
!    r  =   d2r * d
!    eu =   cos (r)
!    ev = - sin (r)
!    u  = eu * f
!    v  = ev * f
!    !--------------------
!    ! tangent linear code
!    !--------------------
!    if (present(du)) du =   f * d2r * dd * ev + eu * df
!    if (present(dv)) dv = - f * d2r * dd * eu + ev * df
!  end subroutine uv_fd
!------------------------------------------------------------------------------
  !----------------------------------------------------
  !> convert u,v to     wind speed and direction
  !>
  !> optionally return partial derivatives.
  !> \param [out] f     wind speed             [m/s]
  !> \param [out] d     wind direction         [degree]
  !> \param [in]  u     wind component u       [m/s]
  !> \param [in]  v     wind component v       [m/s]
  !> \param [out] df_du partial derivative: d f / d u
  !> \param [out] df_dv partial derivative: d f / d v
  !> \param [out] dd_du partial derivative: d d / d u
  !> \param [out] dd_dv partial derivative: d d / d v
  !-----------------------------------------------------
  elemental subroutine fd_uv (f, d, u, v, df_du, df_dv, dd_du, dd_dv)
  real(wp) ,intent(out)           :: f     ! wind speed             [m/s]
  real(wp) ,intent(out)           :: d     ! wind direction         [degree]
  real(wp) ,intent(in)            :: u     ! wind component u       [m/s]
  real(wp) ,intent(in)            :: v     ! wind component v       [m/s]
  real(wp) ,intent(out) ,optional :: df_du ! df / du          [degree / m/s]
  real(wp) ,intent(out) ,optional :: df_dv ! df / dv          [degree / m/s]
  real(wp) ,intent(out) ,optional :: dd_du ! dd / du
  real(wp) ,intent(out) ,optional :: dd_dv ! dd / dv
    real(wp) :: s
    !---------------
    ! nonlinear code
    !---------------
    s = u*u+v*v
    f = sqrt (s)
    d = r2d * atan2 (-u ,-v )
    if(d<0._wp)d=d+360._wp
    !--------------------
    ! tangent linear code
    !--------------------
    if (s == 0._wp) then
      d = 0._wp                         ! undefined
      if (present(df_du)) df_du = 0._wp ! undefined
      if (present(df_dv)) df_dv = 0._wp ! undefined
      if (present(dd_du)) dd_du = 0._wp ! undefined
      if (present(dd_dv)) dd_dv = 0._wp ! undefined
    else
      if (present(df_du)) df_du =   u/f
      if (present(df_dv)) df_dv =   v/f
      if (present(dd_du)) dd_du =   v/s * r2d
      if (present(dd_dv)) dd_dv = - u/s * r2d
    endif

!  subroutine lin_atan2_real_type (x1, x2, y)
!  type (real_type) :: x1, x2
!  type (real_type) :: y
!    y%p%value =                                                           &
!      (x1%p%value * x2%p%nonl%p%value - x2%p%value * x1%p%nonl%p%value) / &
!         -du              -v          -   -dv            -u
!          du               v          -    dv             u
!      (x1%p% nonl%p% value ** 2 + x2%p% nonl%p% value ** 2)
!  end subroutine lin_atan2_real_type

  end subroutine fd_uv
!==============================================================================
  !---------------------------------------------
  !> derive Exner function from p
  !> \param [in] p  pressure                [Pa]
  !> \return        dimensionless pressure
  !---------------------------------------------
  elemental function exner (p)
    real(wp)                        :: exner
    real(wp) ,intent (in)           :: p

    exner = exp ( (R/c_p) * log (p * p0_inv) )
  end function exner
!==============================================================================
  !---------------------------------------------------------------------
  !> derive virtual temperature from virtual potential temperature and p
  !> \param [in] thetav  virtual potential temperature  [K]
  !> \param [in] p       pressure                      [Pa]
  !> \return             virtual temperature            [K]
  !---------------------------------------------------------------------
  elemental function tv_thetav_p (thetav, p) result (tv)
    real(wp)                        :: tv
    real(wp) ,intent (in)           :: thetav
    real(wp) ,intent (in)           :: p

    ! T_v = Theta_v * exner (p)
    tv = thetav * exp ((R/c_p) * log (p * p0_inv) )
  end function tv_thetav_p
!==============================================================================
  !---------------------------------------------------------------------
  !> derive virtual potential temperature from virtual temperature and p
  !> \param [in] tv      virtual temperature            [K]
  !> \param [in] p       pressure                      [Pa]
  !> \return             virtual potential temperature  [K]
  !---------------------------------------------------------------------
  elemental function thetav_tv_p (tv, p) result (thetav)
    real(wp)                        :: thetav
    real(wp) ,intent (in)           :: tv
    real(wp) ,intent (in)           :: p

    ! Theta_v = T_v / exner (p)
    thetav = tv * exp ((- R/c_p) * log (p * p0_inv) )
  end function thetav_tv_p
!==============================================================================
  !---------------------------------------------------------------------
  !> derive pressure from density and virtual potential temperature
  !> \param [in] rho     density                   [kg/m^3]
  !> \param [in] thetav  virtual potential temperature  [K]
  !> \return             pressure                      [Pa]
  !---------------------------------------------------------------------
  elemental function p_rho_thetav (rho, thetav) result (p)
    real(wp)                        :: p
    real(wp) ,intent (in)           :: rho
    real(wp) ,intent (in)           :: thetav

    real(wp), parameter :: a = c_p / (c_p - R)

    p = p0ref * exp (a * log ((R*p0_inv) * rho * thetav) )
  end function p_rho_thetav
!==============================================================================
  !---------------------------------------------------------------------
  !> Convert geopotential or geopotential height in height above geoid
  !>
  !> <b> AltGeoid = geopot_geoid( lat, hgeop ) </b>
  !>
  !> The geopotential is the potential energy due to gravity and therefore
  !> related to the geoid or height above sea level.
  !>
  !> Geopotential: \f$ \Phi = \int_0^h g \, dh  \f$, g = g(h) - gravity \n
  !> Geopotential height: \f$ h_{g} = \frac{1}{g_0} \int_0^h g \, dh \f$,
  !>                      \f$g_0\f$ - standard gravity, \f$g_0 = 9.80665\f$
  !>
  !> Inverse function: ::geoid_geopot
  !>
  !> @param[in] lat
  !>            latitude in radian \n
  !> @param[in] geop
  !>            geopotential or geopotential height in meter,
  !>            depending on "geopot": \n
  !>            geopot = .false. or missing: geop is geopotential height
  !>                                                 in meter \n
  !>            geopot = .true.            : geop is geopotential
  !> @param[in] geopot
  !>            logical giving the quantity in "geopot_geoid" \n
  !>            optional
  !> @return geopot_geoid - height above geoid in meter
  !---------------------------------------------------------------------
  function geopot_geoid(lat, geop, geopot)

    implicit none

    ! List of calling arguments:
    real(wp)                       :: geopot_geoid
    real(wp), intent(in)           :: geop
    real(wp), intent(in)           :: lat
    logical, optional,  intent(in) :: geopot

    ! List of local variables:
    real(wp), parameter  :: g0 = 9.80665_wp
    logical              :: IsPot

    real(wp) :: rE
    real(wp) :: gG
    real(wp) :: SinLat2
    !---------------------------------------------------------------------

    ! decide if "geop" is geopotential or geopotential height
    if (present(geopot)) then
       IsPot = geopot
    else
       ! assume "geop" is geopotential height
       IsPot = .false.
    end if

    SinLat2 = sin(lat)**2

    ! ellipsoidal radius of Earth, depending on latitude, radius in m
    rE = 6378137.0_wp / (1.006803_wp - 0.006706_wp * SinLat2)

    ! local gravity on ellipsoid, depending on latitude
    gG = 9.7803267714_wp *                            &
         (1.0_wp + 0.00193185138639_wp * SinLat2) /   &
         sqrt( 1.0_wp - 0.00669437999013_wp * SinLat2)

    if (IsPot) then
       ! geopotential => height above geoid in m
       geopot_geoid = (rE * geop) / (gG * rE - geop)
    else
       ! geopotential height in m => height above geoid in m
       geopot_geoid = (g0 * rE * geop) / (gG * rE - g0 * geop)
    end if

  end function geopot_geoid
!==============================================================================
  !---------------------------------------------------------------------
  !> Convert height above geoid to geopotential or geopotential height in m.
  !>
  !> <b> Hgp = geoid_geopot(lat, Hgeo) </b>
  !>
  !> The geopotential is the potential energy due to gravity and therefore
  !> related to the geoid or height above sea level.
  !>
  !> Geopotential: \f$ \Phi = \int_0^h g \, dh  \f$, g = g(h) - gravity \n
  !> Geopotential height: \f$ h_{g} = \frac{1}{g_0} \int_0^h g \, dh \f$,
  !>                      \f$g_0\f$ - standard gravity, \f$g_0 = 9.80665\f$
  !>
  !> Inverse function: ::geopot_geoid
  !>
  !> \param[in] lat
  !>            latitude in radian \n
  !> \param[in] Hgeo
  !>            height above geoid in meter
  !> \param[in] geopot
  !>            logical defining the output quantity "geoid_geopot" \n
  !>            optional
  !> \return    geoid_geopot -
  !>            geopotential if geopot = .true. \n
  !>            geopotential height in m if geopot = .false. or missing
  !---------------------------------------------------------------------
  function geoid_geopot(lat, Hgeo, geopot)

    implicit none

    ! List of calling arguments:
    real(wp)                       :: geoid_geopot
    real(wp), intent(in)           :: Hgeo
    real(wp), intent(in)           :: lat
    logical, optional,  intent(in) :: geopot

    ! List of local variables:
    real(wp), parameter  :: g0 = 9.80665_wp
    logical              :: IsPot

    real(wp) :: rE
    real(wp) :: gG
    real(wp) :: SinLat2
    !---------------------------------------------------------------------

    ! decide if "geop" is geopotential or geopotential height
    if (present(geopot)) then
       IsPot = geopot
    else
       ! assume "geop" is geopotential height
       IsPot = .false.
    end if

    SinLat2 = sin(lat)**2

    ! ellipsoidal radius of Earth, depending on latitude, radius in m
    rE = 6378137.0_wp / (1.006803_wp - 0.006706_wp * SinLat2)

    ! local gravity on ellipsoid, depending on latitude
    gG = 9.7803267714_wp *                            &
         (1.0_wp + 0.00193185138639_wp * SinLat2) /   &
         sqrt( 1.0_wp - 0.00669437999013_wp * SinLat2)

    if (IsPot) then
       ! geopotential => height above geoid in m
       geoid_geopot =  (rE * Hgeo * gG) / (rE + Hgeo)
    else
       ! geopotential height in m => height above geoid in m
       geoid_geopot = (rE * Hgeo * gG) / (g0*(rE + Hgeo))
    end if

  end function geoid_geopot
!==============================================================================
  !---------------------------------------------------------------------
  !> TODO
  !>
  subroutine get_gas(str, name, id)
    character(len=*), intent(in)            :: str
    character(len=*), intent(out), optional :: name
    integer,          intent(out), optional :: id
    character(len=len(str)) :: aux
    integer                 :: i
    aux = trim(adjustl(tolower(str)))
    select case(aux)
    case('ozone')
      aux = 'o3'
    end select
    if (present(name)) name = ''
    if (present(id  )) id   = -1
    do i = 1, ngases
      if (trim(aux) == trim(gas_names(i))) then
        if (present(name)) name = gas_names(i)
        if (present(id  )) id   = i
      end if
    end do
  end subroutine get_gas


!==============================================================================
end module mo_physics
