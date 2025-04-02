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

#ifdef _DACE_
#define __COSMO__
#endif

MODULE radar_mie_iface_cosmo_utils

!------------------------------------------------------------------------------
!
! Description:
!      Utilities for the interface functions for the COSMO/ICON microphysics
!      to the EMVORADO libraries for polarimetric radar moments computations.
!
!      Applicable to microwave radiation, wavelength > 1 mm
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_kind , ONLY : dp, wp

#if (defined AUXOUT_OFFLINE && defined __COSMO__)
  USE radar_data , ONLY :   &
       ydate_ini_mod, cmaxlen
  USE radar_data_namelist, only: &
       ydirradarout, lmodfield_output, loutnwp, &
       loutdbz, loutpolstd, loutpolall, lextdbz
  USE radar_data_mie, ONLY : &
       Tmax_i_modelgrid, Tmax_s_modelgrid, Tmax_g_modelgrid, Tmax_h_modelgrid
  USE radar_interface, ONLY : get_datetime_act, &
       t, rho, qv, qc, qi, qr, qs, qg, &
       qh, qnc, qni, qnr, qns, qng, qnh, qgl, qhl, qnc_s
#endif

  USE radar_data_mie, ONLY : &
       inv_rhow, inv_rhow2, &
       inv_rhoi, inv_rhoi2, &
       particle, rain, cloud, snow, ice, graupel, hail, rain_coeffs,  &
       ray_const, vt_const, &
       pi => pi_dp, pi6 => pi6_dp, inv_pi6 => inv_pi6_dp, inv_pi6sq => inv_pi6sq_dp, &
       itype_gscp_fwo

  USE radar_mielib_vec, ONLY : &
       m_complex_water_ray, &
       m_complex_water_ray_vec, m_complex_ice_maetzler_vec

  USE radar_gamma_functions_vec, ONLY: gfct

  USE radar_mie_utils, ONLY : &
       gamfac_imom_DMGD, mu_d_relation_seifert, &
       particle_assign

!===============================================================================
!===============================================================================

  IMPLICIT NONE

!===============================================================================
!==============================================================================

  PUBLIC

!==============================================================================

CONTAINS

!===============================================================================
!===============================================================================

  !==============================================================================
  !==============================================================================
  !
  ! UTILITIES FOR GRID POINT REFLECTIVITY/EXTINCTION CALCULATION FOR COSMO/ICON
  !
  !=======================================================================
  !=======================================================================

  ! Initialization of prefactors for reflectivity calculations in the global struct "ray_const":
  SUBROUTINE init_radar_rayleigh_consts(ldo_qc, ldo_qr, ldo_qi, ldo_qs, ldo_qg, ldo_qh, &
       lambda_radar)

    IMPLICIT NONE

    ! .. Input/output variables:
    LOGICAL, INTENT(in)        :: ldo_qc, ldo_qr, ldo_qi, ldo_qs, ldo_qg, ldo_qh
    REAL(KIND=dp), INTENT(in)  :: lambda_radar

    ! .. Local variables:
    INTEGER, SAVE :: firstcall   = 0, &
                     firstcall_c = 0, &
                     firstcall_r = 0, &
                     firstcall_i = 0, &
                     firstcall_s = 0, &
                     firstcall_g = 0, &
                     firstcall_h = 0, &
                     gscps = 0

    REAL(KIND=dp), SAVE :: lambda_radars=0d0
    TYPE(particle), SAVE   :: cloud_speicher
    TYPE(particle), SAVE   :: rain_speicher
    TYPE(particle), SAVE   :: ice_speicher
    TYPE(particle), SAVE   :: snow_speicher
    TYPE(particle), SAVE   :: graupel_speicher
    TYPE(particle), SAVE   :: hail_speicher

    LOGICAL :: kw0_has_changed

    REAL(KIND=dp) :: n0_nu

    ! 0) Value for K_w_0 depending on radar wavelength:
    !
    IF (firstcall /= 1 .OR. lambda_radar /= lambda_radars) THEN
      ! Common value of Kw0**2, mainly applied in radar software
      ray_const%m_w_0 = m_complex_water_ray(lambda_radar, 0d0)
      ray_const%K_w_0 = ABS( (ray_const%m_w_0*ray_const%m_w_0-1d0)/&
                             (ray_const%m_w_0*ray_const%m_w_0+2d0) )
      ray_const%K_w_0 = ray_const%K_w_0*ray_const%K_w_0

      ray_const%pi6sq_K0_fac = 1d18 * inv_pi6sq / ray_const%K_w_0
      ray_const%C_fac = pi**5 / lambda_radar**4
      ray_const%Z_fac = 1d18 / (ray_const%C_fac * ray_const%K_w_0)

      lambda_radars = lambda_radar
      firstcall = 1
      kw0_has_changed = .TRUE.
      ! not dependent on cloud mgd params, hence
      ! - no need to recalc them unless kw0 has changed
      ! - no problem if cloud variable is not defined
      ! therefore do here instead below and don't put behind a ldo_qc check
      ! (also, we can easily affor to do the calc even if not needed - an
      ! if-check might cost more...)
      ray_const%cloud_extprefac = 6d0*pi*inv_rhow / lambda_radar
    ELSE
      kw0_has_changed = .FALSE.
    END IF

    ! 1) Constant (pre)factors for Rayleigh-Oguchi-Approx. of cloud droplets (common for 1&2mom):
    !
    IF (ldo_qc) THEN
      IF (firstcall_c /= 1 .OR. kw0_has_changed .OR. &
          cloud_speicher%b_geo /= cloud%b_geo .OR. &
          cloud_speicher%mu /= cloud%mu .OR. cloud_speicher%nu /= cloud%nu) THEN

        ray_const%cloud_Zprefac = gamfac_imom_DMGD(cloud,2d0) * ray_const%pi6sq_K0_fac
        ray_const%cloud_Zprefac_rho = ray_const%cloud_Zprefac * inv_rhow2

        CALL particle_assign(cloud_speicher,cloud)
        firstcall_c = 1
      END IF
    END IF

    ! 2) Constant (pre)factors for Rayleigh-Oguchi-Approx. of rain drops (separate for 1&2mom):
    !
    IF (ldo_qr) THEN
      IF (firstcall_r /= 1 .OR. kw0_has_changed .OR. &
          rain_speicher%n0_const /= rain%n0_const .OR. &
          rain_speicher%b_geo /= rain%b_geo .OR. rain_speicher%a_geo /= rain%a_geo .OR. &
          rain_speicher%mu /= rain%mu .OR. rain_speicher%nu /= rain%nu) THEN

        ray_const%rain2mom_Zprefac_rho = gamfac_imom_DMGD(rain,2d0) * &
                                         ray_const%pi6sq_K0_fac * inv_rhow2

        n0_nu = rain%n0_const/rain%nu
        ray_const%rain1mom_expo = (rain%mu+2d0*rain%b_geo+1d0) / (rain%mu+rain%b_geo+1d0)
        ! n0 of rain is constant, ie we can include its contribution directly into the Z-prefactor
        ray_const%rain1mom_Zprefac_n0_rho = gfct((rain%mu+2d0*rain%b_geo+1d0)/rain%nu) * &
                                            (gfct((rain%mu+rain%b_geo+1d0)/rain%nu)* &
                                              rain%a_geo*n0_nu)**(-ray_const%rain1mom_expo) * &
                                             rain%a_geo*rain%a_geo*n0_nu * &
                                             ray_const%pi6sq_K0_fac * inv_rhow2

        CALL particle_assign(rain_speicher,rain)
        firstcall_r = 1
      END IF
    END IF

    ! 3) Constant (pre)factors for Rayleigh-Oguchi-Approx. of cloud ice (1&2mom):
    !
    IF (ldo_qi) THEN
      IF (firstcall_i /= 1 .OR. kw0_has_changed .OR. &
          ice_speicher%b_geo /= ice%b_geo .OR. &
          ice_speicher%mu /= ice%mu .OR. ice_speicher%nu /= ice%nu) THEN

        ! 1mom: monodisperse size distribution, therefore gamfac_imom_DMGD(ice,2d0) = 1.0
        ray_const%ice1mom_Zprefac     = ray_const%pi6sq_K0_fac
        ray_const%ice1mom_Zprefac_rho = ray_const%ice1mom_Zprefac * inv_rhoi2

        ray_const%ice2mom_Zprefac     = gamfac_imom_DMGD(ice,2d0) * ray_const%pi6sq_K0_fac
        ray_const%ice2mom_Zprefac_rho = ray_const%ice2mom_Zprefac * inv_rhoi2
        ray_const%ice_Davfac          = gamfac_imom_DMGD(ice,(1d0/ice%b_geo))

        CALL particle_assign(ice_speicher,ice)
        firstcall_i = 1
      END IF
    END IF

    ! 4) Constant (pre)factors for Rayleigh-Oguchi-Approx. of snow (separate for 1&2mom):
    !
    IF (ldo_qs) THEN
      IF (firstcall_s /= 1 .OR. kw0_has_changed .OR. &
           snow_speicher%b_geo /= snow%b_geo .OR. snow_speicher%a_geo /= snow%a_geo .OR. &
           snow_speicher%mu /= snow%mu .OR. snow_speicher%nu /= snow%nu) THEN

        ray_const%snow2mom_Zprefac = gamfac_imom_DMGD(snow,2d0) * ray_const%pi6sq_K0_fac
        ray_const%snow2mom_Zprefac_rho = ray_const%snow2mom_Zprefac * inv_rhoi2
        ray_const%snow_Davfac = gamfac_imom_DMGD(snow,(1d0/snow%b_geo))

        n0_nu = 1d0/snow%nu
        ray_const%snow1mom_expo = (snow%mu+2d0*snow%b_geo+1d0) / (snow%mu+snow%b_geo+1d0)
        ray_const%snow_gam_mub1nu = gfct((snow%mu+snow%b_geo+1d0)/snow%nu)
        ray_const%snow1mom_Zprefac = gfct((snow%mu+2d0*snow%b_geo+1d0)/snow%nu) * &
                                     (ray_const%snow_gam_mub1nu*snow%a_geo*n0_nu)** &
                                       (-ray_const%snow1mom_expo) * &
                                      snow%a_geo*snow%a_geo*n0_nu * &
                                      ray_const%pi6sq_K0_fac
        ray_const%snow1mom_Zprefac_rho = ray_const%snow1mom_Zprefac * inv_rhoi2

        CALL particle_assign(snow_speicher,snow)
        firstcall_s = 1
      END IF
    END IF

    ! 5) Constant (pre)factors for Rayleigh-Oguchi-Approx. of graupel (separate for 1&2mom):
    !
    IF (ldo_qg) THEN
      IF (firstcall_g /= 1 .OR. kw0_has_changed .OR. &
          graupel_speicher%n0_const /= graupel%n0_const .OR. &
          graupel_speicher%b_geo /= graupel%b_geo .OR. graupel_speicher%a_geo /= graupel%a_geo .OR. &
          graupel_speicher%mu /= graupel%mu .OR. graupel_speicher%nu /= graupel%nu) THEN

        ray_const%graupel2mom_Zprefac = gamfac_imom_DMGD(graupel,2d0) * ray_const%pi6sq_K0_fac
        ray_const%graupel2mom_Zprefac_rho = ray_const%graupel2mom_Zprefac * inv_rhoi2
        ray_const%graupel_Davfac = gamfac_imom_DMGD(graupel,(1d0/graupel%b_geo))

        n0_nu = graupel%n0_const/graupel%nu
        ray_const%graupel_gam_mub1nu = gfct((graupel%mu+graupel%b_geo+1d0)/graupel%nu)
        ray_const%graupel1mom_expo = (graupel%mu+2d0*graupel%b_geo+1d0) / (graupel%mu+graupel%b_geo+1d0)
        ! n0 of graupel is constant, ie we can include its contribution directly into the Z-prefactor
        ray_const%graupel1mom_Zprefac_n0 = gfct((graupel%mu+2d0*graupel%b_geo+1d0)/graupel%nu) * &
                                           (ray_const%graupel_gam_mub1nu*graupel%a_geo*n0_nu)** &
                                             (-ray_const%graupel1mom_expo) * &
                                           graupel%a_geo*graupel%a_geo*n0_nu * &
                                           ray_const%pi6sq_K0_fac
        ray_const%graupel1mom_Zprefac_n0_rho = ray_const%graupel1mom_Zprefac_n0 * inv_rhoi2

        CALL particle_assign(graupel_speicher,graupel)
        firstcall_g = 1
      END IF
    END IF

    ! 6) Constant (pre)factors for Rayleigh-Oguchi-Approx. of hail (2-mom. scheme):
    !
    IF (ldo_qh) THEN
      IF (firstcall_h /= 1 .OR. kw0_has_changed .OR. &
          hail_speicher%b_geo /= hail%b_geo .OR. &
          hail_speicher%mu /= hail%mu .OR. hail_speicher%nu /= hail%nu) THEN

        ray_const%hail2mom_Zprefac = gamfac_imom_DMGD(hail,2d0) * ray_const%pi6sq_K0_fac
        ray_const%hail2mom_Zprefac_rho = ray_const%hail2mom_Zprefac * inv_rhoi2
        ray_const%hail_Davfac = gamfac_imom_DMGD(hail,(1d0/hail%b_geo))

        CALL particle_assign(hail_speicher,hail)
        firstcall_h = 1
      END IF
    END IF

  END SUBROUTINE init_radar_rayleigh_consts


  ! Initialization of prefactors for mean/weighted fallspeed calculations in the global struct vt_const:
  SUBROUTINE init_radar_vt_oguchi (ldo_qc, ldo_qr, ldo_qi, ldo_qs, ldo_qg, ldo_qh, &
       lambda_radar, lwdbz)

    IMPLICIT NONE

    ! .. Input/output variables:
    LOGICAL, INTENT(in)        :: ldo_qc, ldo_qr, ldo_qi, ldo_qs, ldo_qg, ldo_qh, lwdbz
    REAL(KIND=dp), INTENT(in)  :: lambda_radar

    ! .. Local variables:
    INTEGER, SAVE :: firstcall   = 0, &
                     firstcall_c = 0, &
                     firstcall_r = 0, &
                     firstcall_i = 0, &
                     firstcall_s = 0, &
                     firstcall_g = 0, &
                     firstcall_h = 0

    REAL(KIND=dp), SAVE :: lambda_radars=0d0
    LOGICAL,       SAVE :: lwdbzs
    REAL(KIND=dp)       :: param_mu, param_bg, param_v, param_bv, param_vl, param_bvl, &
                           fac_gam_bg, n0_nu, Kfac, rhofac
    COMPLEX(KIND=dp)    :: m_w_0

    TYPE(particle), SAVE :: cloud_speicher
    TYPE(particle), SAVE :: rain_speicher
    TYPE(particle), SAVE :: ice_speicher, rain_speicher_i
    TYPE(particle), SAVE :: snow_speicher, rain_speicher_s
    TYPE(particle), SAVE :: graupel_speicher, rain_speicher_g
    TYPE(particle), SAVE :: hail_speicher, rain_speicher_h

    LOGICAL :: kw0_has_changed, lwdbz_has_changed

    ! Blending exponent for the fall velocity of partially melted particles (all types):
    !  vt_melt = vt_dry*(1-fm)^am + vt_drop*fm^am
    !        fm      = degree of melting of the actual particle
    !        vt_dry  = fallspeed of a dry particle of the same mass
    !        vt_drop = fallspeed of a drop of the same mass
    !
    vt_const%am = 2d0    ! 2.0 is recommended by Axel
!    vt_const%am = 3.5d0  ! 3.5 is in the appendix of Zeng et al. (2016)

    ! 0) Value for K_w_0 depending on radar wavelength:
    !
    IF (firstcall /= 1 .OR. lambda_radar /= lambda_radars) THEN

      m_w_0 = m_complex_water_ray(lambda_radar, 0d0)
      vt_const%K_w_0 = ABS( (m_w_0*m_w_0-1d0) / (m_w_0*m_w_0+2d0) )
      vt_const%K_w_0 = vt_const%K_w_0*vt_const%K_w_0

      lambda_radars = lambda_radar
      lwdbzs = lwdbz
      firstcall = 1
      kw0_has_changed = .TRUE.
    ELSE
      kw0_has_changed = .FALSE.
    END IF

    IF (lwdbz) THEN
      Kfac  = 1d18 * inv_pi6sq / vt_const%K_w_0
      param_vl = 2d0+rain%b_vel
    ELSE
      Kfac  = 1d0
      rhofac   = 1d0
      param_vl = rain%b_vel
    END IF

    IF (lwdbzs .NEQV. lwdbz) THEN
      lwdbz_has_changed = .TRUE.
    ELSE
      lwdbz_has_changed = .FALSE.
    END IF
    
    ! 1) Constant (pre)factors for Rayleigh-Oguchi-Approx. of cloud droplets:
    !
    IF (ldo_qc) THEN
      IF (firstcall_c /= 1 .OR. kw0_has_changed .OR. lwdbz_has_changed .OR. &
          cloud_speicher%b_geo /= cloud%b_geo .OR. &
          cloud_speicher%a_vel /= cloud%a_vel .OR. cloud_speicher%b_vel /= cloud%b_vel .OR. &
          cloud_speicher%mu /= cloud%mu .OR. cloud_speicher%nu /= cloud%nu) THEN

        IF (lwdbz) THEN
          rhofac  = inv_rhow2
          param_v = 2d0+cloud%b_vel
        ELSE
          param_v = cloud%b_vel
        END IF

        vt_const%cloud_vtexpo   = param_v - 1d0
        vt_const%cloud_vtprefac = gamfac_imom_DMGD(cloud,param_v) * &
                                  cloud%a_vel * rhofac * Kfac

        CALL particle_assign(cloud_speicher,cloud)
        firstcall_c = 1
      END IF
    END IF

    ! 2) Constant (pre)factors for Rayleigh-Oguchi-Approx. of rain drops:
    !
    IF (ldo_qr) THEN
      IF (firstcall_r /= 1 .OR. kw0_has_changed .OR. lwdbz_has_changed .OR. &
          rain_speicher%a_geo /= rain%a_geo .OR. rain_speicher%b_geo /= rain%b_geo .OR. &
          rain_speicher%a_vel /= rain%a_vel .OR. rain_speicher%b_vel /= rain%b_vel .OR. &
          rain_speicher%mu /= rain%mu .OR. rain_speicher%nu /= rain%nu) THEN

        IF (lwdbz) THEN
          rhofac  = inv_rhow2
        END IF

        IF (itype_gscp_fwo < 200) THEN
          param_bv = (rain%mu+rain%b_geo*param_vl+1d0)/rain%nu
          param_bg = (rain%mu+rain%b_geo+1d0)/rain%nu
          param_mu = (rain%mu+1d0)/rain%nu

          n0_nu = rain%n0_const/rain%nu
          fac_gam_bg = rain%a_geo * n0_nu * gfct(param_bg)

          vt_const%rain_vtexpo   = param_bv / param_bg
          vt_const%rain_nexpo    = param_mu / param_bg

          vt_const%rain_vtprefac = &
                  rain%a_vel*rain%a_geo**param_vl * n0_nu * &
                  gfct(param_bv) * fac_gam_bg**(-vt_const%rain_vtexpo) * Kfac * rhofac
          vt_const%rain_nprefac  = &
                  n0_nu * gfct(param_mu) * fac_gam_bg**(-vt_const%rain_nexpo)
        ELSE
          vt_const%rain_vtexpo   = param_vl - 1d0
          vt_const%rain_vtprefac = gamfac_imom_DMGD(rain,param_vl) * &
                                   rain%a_vel * rhofac * Kfac
        END IF

        CALL particle_assign(rain_speicher,rain)
        firstcall_r = 1
      END IF
    END IF

    ! 3) Constant (pre)factors for Rayleigh-Oguchi-Approx. of cloud ice:
    !
    IF (ldo_qi) THEN
      IF (firstcall_i /= 1 .OR. kw0_has_changed .OR. lwdbz_has_changed .OR. &
          ice_speicher%b_geo /= ice%b_geo .OR. &
          ice_speicher%a_vel /= ice%a_vel .OR. ice_speicher%b_vel /= ice%b_vel .OR. &
          rain_speicher_i%a_vel /= rain%a_vel .OR. rain_speicher_i%b_vel /= rain%b_vel .OR. &
          ice_speicher%nu /= ice%nu .OR. ice_speicher%mu /= ice%mu) THEN

        vt_const%ice_Davfac = gamfac_imom_DMGD(ice,(1d0/ice%b_geo))
        IF (lwdbz) THEN
          rhofac  = inv_rhoi2
          param_v = 2d0+ice%b_vel
        ELSE
          param_v = ice%b_vel
        END IF

        vt_const%ice_vtexpo    = param_v  - 1d0
        vt_const%ice_vtlexpo   = param_vl - 1d0

        IF (itype_gscp_fwo < 200) THEN
          ! Monodisperse size distribution assumed, therfore gamfac_imom_DMGD(ice,...) = 1.0
          vt_const%ice_vtlprefac = rain%a_vel * Kfac
          vt_const%ice_vtfprefac = ice%a_vel  * Kfac
        ELSE
          vt_const%ice_vtlprefac = gamfac_imom_DMGD (ice,param_vl) * &
                                   rain%a_vel * Kfac
          vt_const%ice_vtfprefac = gamfac_imom_DMGD (ice,param_v) * &
                                   ice%a_vel  * Kfac
        END IF
        vt_const%ice_vtprefac  = vt_const%ice_vtfprefac * rhofac

        CALL particle_assign(rain_speicher_i,rain)
        CALL particle_assign(ice_speicher,ice)
        firstcall_i = 1
      END IF
    END IF

    ! 4) Constant (pre)factors for Rayleigh-Oguchi-Approx. of snow:
    !
    IF (ldo_qs) THEN
      IF (firstcall_s /= 1 .OR. kw0_has_changed .OR. lwdbz_has_changed .OR. &
          snow_speicher%a_geo /= snow%a_geo .OR. snow_speicher%b_geo /= snow%b_geo .OR. &
          snow_speicher%a_vel /= snow%a_vel .OR. snow_speicher%b_vel /= snow%b_vel .OR. &
          rain_speicher_s%a_vel /= rain%a_vel .OR. rain_speicher_s%b_vel /= rain%b_vel .OR. &
          snow_speicher%mu /= snow%mu .OR. snow_speicher%nu /= snow%nu) THEN

        vt_const%snow_Davfac = gamfac_imom_DMGD(snow,(1d0/snow%b_geo))
        vt_const%snow_gam_mub1nu = gfct((snow%mu+snow%b_geo+1d0)/snow%nu)
        IF (lwdbz) THEN
          rhofac   = inv_rhoi2
          param_v  = 2d0+snow%b_vel
        ELSE
          param_v  = snow%b_vel
        END IF

        IF (itype_gscp_fwo < 200) THEN
          param_bv  = (snow%mu+snow%b_geo*param_v+1d0)/snow%nu
          param_bvl = (snow%mu+snow%b_geo*param_vl+1d0)/snow%nu
          param_bg  = (snow%mu+snow%b_geo+1d0)/snow%nu
          param_mu  = (snow%mu+1d0)/snow%nu

          n0_nu = 1d0/snow%nu
          fac_gam_bg = snow%a_geo * n0_nu * gfct(param_bg)

          vt_const%snow_vtexpo    = param_bv / param_bg
          vt_const%snow_vtlexpo   = param_bvl / param_bg
          vt_const%snow_nexpo     = param_mu / param_bg

          vt_const%snow_vtlprefac = &
                  rain%a_vel*snow%a_geo**param_vl * n0_nu * &
                  gfct(param_bvl) * fac_gam_bg**(-vt_const%snow_vtlexpo) * Kfac
          vt_const%snow_vtfprefac = &
                  snow%a_vel*snow%a_geo**param_v * n0_nu * &
                  gfct(param_bv) * fac_gam_bg**(-vt_const%snow_vtexpo) * Kfac
          vt_const%snow_vtprefac  = vt_const%snow_vtfprefac * rhofac

          vt_const%snow_nprefac   = &
                  n0_nu * gfct(param_mu) * fac_gam_bg**(-vt_const%snow_nexpo)
        ELSE
          vt_const%snow_vtexpo    = param_v - 1d0
          vt_const%snow_vtlexpo   = param_vl - 1d0

          vt_const%snow_vtlprefac = gamfac_imom_DMGD (snow,param_vl) * &
                                    rain%a_vel * Kfac
          vt_const%snow_vtfprefac = gamfac_imom_DMGD (snow,param_v) * &
                                    snow%a_vel * Kfac
          vt_const%snow_vtprefac  = vt_const%snow_vtfprefac * rhofac
        END IF

        CALL particle_assign(rain_speicher_s,rain)
        CALL particle_assign(snow_speicher,snow)
        firstcall_s = 1
      END IF
    END IF

    ! 5) Constant (pre)factors for Rayleigh-Oguchi-Approx. of graupel:
    !
    IF (ldo_qg) THEN
      IF (firstcall_g /= 1 .OR. kw0_has_changed .OR. lwdbz_has_changed .OR. &
          graupel_speicher%a_geo /= graupel%a_geo .OR. graupel_speicher%b_geo /= graupel%b_geo .OR. &
          graupel_speicher%a_vel /= graupel%a_vel .OR. graupel_speicher%b_vel /= graupel%b_vel .OR. &
          rain_speicher_g%a_vel /= rain%a_vel .OR. rain_speicher_g%b_vel /= rain%b_vel .OR. &
          graupel_speicher%mu /= graupel%mu .OR. graupel_speicher%nu /= graupel%nu) THEN

        vt_const%graupel_Davfac = gamfac_imom_DMGD(graupel,(1d0/graupel%b_geo))
        vt_const%graupel_gam_mub1nu = gfct((graupel%mu+graupel%b_geo+1d0)/graupel%nu)
        IF (lwdbz) THEN
          rhofac   = inv_rhoi2
          param_v  = 2d0+graupel%b_vel
        ELSE
          param_v  = graupel%b_vel
        END IF

        IF (itype_gscp_fwo < 200) THEN
          param_bv  = (graupel%mu+graupel%b_geo*param_v+1d0)/graupel%nu
          param_bvl = (graupel%mu+graupel%b_geo*param_vl+1d0)/graupel%nu
          param_bg  = (graupel%mu+graupel%b_geo+1d0)/graupel%nu
          param_mu  = (graupel%mu+1d0)/graupel%nu

          n0_nu = graupel%n0_const/graupel%nu
          fac_gam_bg = graupel%a_geo * n0_nu * gfct(param_bg)

          vt_const%graupel_vtexpo  = param_bv / param_bg
          vt_const%graupel_vtlexpo = param_bvl / param_bg
          vt_const%graupel_nexpo   = param_mu / param_bg

          vt_const%graupel_vtlprefac = &
                  rain%a_vel*graupel%a_geo**param_vl * n0_nu * &
                  gfct(param_bvl) * fac_gam_bg**(-vt_const%graupel_vtlexpo) * Kfac
          vt_const%graupel_vtfprefac = &
                  graupel%a_vel*graupel%a_geo**param_v * n0_nu * &
                  gfct(param_bv) * fac_gam_bg**(-vt_const%graupel_vtexpo) * Kfac
          vt_const%graupel_vtprefac  = vt_const%graupel_vtfprefac * rhofac

          vt_const%graupel_nprefac   = &
                  n0_nu * gfct(param_mu) * fac_gam_bg**(-vt_const%graupel_nexpo)
        ELSE
          vt_const%graupel_vtexpo    = param_v - 1d0
          vt_const%graupel_vtlexpo   = param_vl - 1d0

          vt_const%graupel_vtlprefac = gamfac_imom_DMGD (graupel,param_vl) * &
                                       rain%a_vel * Kfac
          vt_const%graupel_vtfprefac = gamfac_imom_DMGD (graupel,param_v) * &
                                       graupel%a_vel * Kfac
          vt_const%graupel_vtprefac  = vt_const%graupel_vtfprefac * rhofac
        END IF

        CALL particle_assign(rain_speicher_g,rain)
        CALL particle_assign(graupel_speicher,graupel)
        firstcall_g = 1
      END IF
    END IF

    ! 6) Constant (pre)factors for Rayleigh-Oguchi-Approx. of hail:
    !
    IF (ldo_qh) THEN
      IF (firstcall_h /= 1 .OR. kw0_has_changed .OR. lwdbz_has_changed .OR. &
           hail_speicher%b_geo /= hail%b_geo .OR. &
           hail_speicher%a_vel /= hail%a_vel .OR. hail_speicher%b_vel /= hail%b_vel .OR. &
           rain_speicher_h%a_vel /= rain%a_vel .OR. rain_speicher_h%b_vel /= rain%b_vel .OR. &
           hail_speicher%mu /= hail%mu .OR. hail_speicher%nu /= hail%nu) THEN

        vt_const%hail_Davfac = gamfac_imom_DMGD(hail,(1d0/hail%b_geo))
        IF (lwdbz) THEN
          rhofac   = inv_rhoi2
          param_v  = 2d0+hail%b_vel
        ELSE
          param_v  = hail%b_vel
        END IF

        IF (itype_gscp_fwo >= 200) THEN
          vt_const%hail_vtexpo    = param_v - 1d0
          vt_const%hail_vtlexpo   = param_vl - 1d0

          vt_const%hail_vtlprefac = gamfac_imom_DMGD (hail,param_vl) * &
                                    rain%a_vel * Kfac
          vt_const%hail_vtfprefac = gamfac_imom_DMGD (hail,param_v) * &
                                    hail%a_vel * Kfac
          vt_const%hail_vtprefac  = vt_const%hail_vtfprefac * rhofac
        END IF

        CALL particle_assign(rain_speicher_h,rain)
        CALL particle_assign(hail_speicher,hail)
        firstcall_h = 1
      END IF
    END IF

  END SUBROUTINE init_radar_vt_oguchi


  FUNCTION vtradar_normalize(&
      vt_local, z_local, rho, rho_0, ni, nj, nk, ilow, iup, jlow, jup, ku, ko) &
      RESULT(vt_radar)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ni, nj, nk, ilow, iup, jlow, jup, ku, ko
    REAL(KIND=dp), INTENT(IN) :: rho_0
    REAL(KIND=wp), DIMENSION(ni,nj,nk), INTENT(IN) :: rho
    REAL(KIND=dp), DIMENSION(ni,nj,nk), INTENT(IN) :: vt_local, z_local

    REAL(KIND=dp), DIMENSION(ni,nj,nk) :: vt_radar

    INTEGER :: i, j, k

!$OMP PARALLEL DO PRIVATE(i,j,k)
    DO k= ku, ko
      DO j = jlow, jup
        DO i = ilow, iup

          IF (z_local(i,j,k) > 0d0) THEN
            vt_radar(i,j,k) = SQRT(rho_0/rho(i,j,k)) * vt_local(i,j,k)/z_local(i,j,k)
          ELSE
            vt_radar(i,j,k) = 0d0
          END IF

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END FUNCTION vtradar_normalize

  !===========================================================================================================
  !===========================================================================================================
  !
  !  End of interface routines.
  !
  !===========================================================================================================
  !===========================================================================================================

#if (defined AUXOUT_OFFLINE && defined __COSMO__)
  SUBROUTINE write_modelfield_output(station_id, itype_refl, &
                                     zh_radar, ah_radar, &
                                     zv_radar, rrhv_radar, irhv_radar, &
                                     kdp_radar, adp_radar, zvh_radar)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: station_id, itype_refl
    REAL(kind=dp), OPTIONAL, INTENT(IN) :: zh_radar(:,:,:), ah_radar(:,:,:), &
                                           zv_radar(:,:,:), &
                                           rrhv_radar(:,:,:), irhv_radar(:,:,:), &
                                           kdp_radar(:,:,:), adp_radar(:,:,:), &
                                           zvh_radar(:,:,:)

    INTEGER :: fehler, ni, nj
    CHARACTER (LEN=32) :: field_desc
    REAL(kind=wp), ALLOCATABLE :: Tmax(:,:,:)

    WRITE(field_desc,'(2A,I6.6)') 'modgrid','-id-',station_id

    IF (loutdbz) THEN
      IF (PRESENT(zh_radar)) &
        CALL write_modelfield_param(&
             TRIM(ydirradarout), 'zh_radar', TRIM(field_desc), zh_radar, &
             '# (Horizontal) Reflectivity [mm^6 m^−3]', .FALSE., fehler)
      IF (lextdbz .AND. PRESENT(ah_radar)) &
        CALL write_modelfield_param(&
             TRIM(ydirradarout), 'ah_radar', TRIM(field_desc), ah_radar, &
             '# (Horizontal) Attenuation coefficient [1/km (?)]', .FALSE., fehler)

      ! Polarimetric parameters only make sense for TMatrix (or: non-sphere permitting)
      IF ((loutpolstd .OR. loutpolall) .AND. itype_refl > 4) THEN
        IF (PRESENT(zv_radar)) &
          CALL write_modelfield_param(&
               TRIM(ydirradarout), 'zv_radar', TRIM(field_desc), zv_radar, &
               '# Vertical reflectivity [mm^6 m^−3]', .FALSE., fehler)
        IF (PRESENT(rrhv_radar)) &
          CALL write_modelfield_param(&
               TRIM(ydirradarout), 'rrhv_radar', TRIM(field_desc), rrhv_radar, &
               '# Real co-polar correlation coefficient nominator [mm^6 m^−3]', .FALSE., fehler)
        IF (PRESENT(irhv_radar)) &
          CALL write_modelfield_param(&
               TRIM(ydirradarout), 'irhv_radar', TRIM(field_desc), irhv_radar, &
               '# Imag co-polar correlation coefficient nominator [mm^6 m^−3]', .FALSE., fehler)
        IF (PRESENT(kdp_radar)) &
          CALL write_modelfield_param(&
               TRIM(ydirradarout), 'kdp_radar', TRIM(field_desc), kdp_radar, &
               '# Specific differential phase [deg/km (?)]', .FALSE., fehler)
        IF (lextdbz .AND. PRESENT(adp_radar)) &
          CALL write_modelfield_param(&
               TRIM(ydirradarout), 'adp_radar', TRIM(field_desc), adp_radar, &
               '# Differential attenuation coefficient [1/km (?)]', .FALSE., fehler)
        IF (loutpolall .AND. PRESENT(zvh_radar)) &
          CALL write_modelfield_param(&
               TRIM(ydirradarout), 'zvh_radar', TRIM(field_desc), zvh_radar, &
               '# Horizontal-turned-vertical reflectivity (for LDR) [mm^6 m^−3]', .FALSE., fehler)
      END IF
    END IF

    IF (loutnwp) THEN
      WRITE(field_desc,'(A)') 'modgrid'

      CALL write_modelfield_param(TRIM(ydirradarout), 'rho', TRIM(field_desc), rho, &
                                  '# density [kg/m3]', .FALSE., fehler)
      CALL write_modelfield_param(TRIM(ydirradarout), 't', TRIM(field_desc), t, &
                                  '# temperature [K]', .FALSE., fehler)
      CALL write_modelfield_param(TRIM(ydirradarout), 'qc', TRIM(field_desc), qc, &
                                  '# hydromet mass concentration [kg/kg]', .FALSE., fehler)
      CALL write_modelfield_param(TRIM(ydirradarout), 'qr', TRIM(field_desc), qr, &
                                  '# hydromet mass concentration [kg/kg]', .FALSE., fehler)

      ni = SIZE(Tmax_i_modelgrid, dim=1)
      nj = SIZE(Tmax_i_modelgrid, dim=2)
      ALLOCATE(Tmax(ni,nj,1))
      
      IF (itype_gscp_fwo >= 130) THEN
        CALL write_modelfield_param(TRIM(ydirradarout), 'qi', TRIM(field_desc), qi, &
                                    '# hydromet mass concentration [kg/kg]', .FALSE., fehler)
        Tmax(:,:,1) = Tmax_i_modelgrid
        CALL write_modelfield_param(TRIM(ydirradarout), 'Tmax_i', TRIM(field_desc), Tmax, &
                                    '# hydromet max melting temparature [K]', .FALSE., fehler)
      END IF
      IF (itype_gscp_fwo >= 140) THEN
        CALL write_modelfield_param(TRIM(ydirradarout), 'qs', TRIM(field_desc), qs, &
                                    '# hydromet mass concentration [kg/kg]', .FALSE., fehler)
        Tmax(:,:,1) = Tmax_s_modelgrid
        CALL write_modelfield_param(TRIM(ydirradarout), 'Tmax_s', TRIM(field_desc), Tmax, &
                                    '# hydromet max melting temparature [K]', .FALSE., fehler)
      END IF
      IF (itype_gscp_fwo >= 150) THEN
        CALL write_modelfield_param(TRIM(ydirradarout), 'qg', TRIM(field_desc), qg, &
                                    '# hydromet mass concentration [kg/kg]', .FALSE., fehler)
        Tmax(:,:,1) = Tmax_g_modelgrid
        CALL write_modelfield_param(TRIM(ydirradarout), 'Tmax_g', TRIM(field_desc), Tmax, &
                                    '# hydromet max melting temparature [K]', .FALSE., fehler)
      END IF

      IF (itype_gscp_fwo >= 200) THEN
        CALL write_modelfield_param(TRIM(ydirradarout), 'qnc', TRIM(field_desc), qnc, &
                                    '# hydromet number concentration [#/kg]', .FALSE., fehler)
        CALL write_modelfield_param(TRIM(ydirradarout), 'qni', TRIM(field_desc), qni, &
                                    '# hydromet number concentration [#/kg]', .FALSE., fehler)
        CALL write_modelfield_param(TRIM(ydirradarout), 'qnr', TRIM(field_desc), qnr, &
                                    '# hydromet number concentration [#/kg]', .FALSE., fehler)
        CALL write_modelfield_param(TRIM(ydirradarout), 'qns', TRIM(field_desc), qns, &
                                    '# hydromet number concentration [#/kg]', .FALSE., fehler)
        CALL write_modelfield_param(TRIM(ydirradarout), 'qng', TRIM(field_desc), qng, &
                                    '# hydromet number concentration [#/kg]', .FALSE., fehler)
      END IF

      IF (itype_gscp_fwo >= 260) THEN
        CALL write_modelfield_param(TRIM(ydirradarout), 'qh', TRIM(field_desc), qh, &
                                    '# hydromet mass concentration [kg/kg]', .FALSE., fehler)
        CALL write_modelfield_param(TRIM(ydirradarout), 'qnh', TRIM(field_desc), qnh, &
                                    '# hydromet number concentration [#/kg]', .FALSE., fehler)
        Tmax(:,:,1) = Tmax_h_modelgrid
        CALL write_modelfield_param(TRIM(ydirradarout), 'Tmax_h', TRIM(field_desc), Tmax, &
                                    '# hydromet max melting temparature [K]', .FALSE., fehler)
      END IF
      DEALLOCATE(Tmax)
    END IF

   RETURN
END SUBROUTINE write_modelfield_output

SUBROUTINE write_modelfield_param(outdir, cvar, levtyp, feld, header, gzipflag, fehler)

  !-------------------------------------------------------------------------------
  !
  ! Description:
  !   This module procedure of module "write_modelfield_param" writes out
  !   parameters given on the model grid.
  !   It is essentially a copy of schreibe_feld from the offline operator's main
  !   file (radvop-offline's src/cosmo_refl_offline.f90)
  !
  ! ...
  !
  !===============================================================================

  !===============================================================================
  !
  ! Local USE statements:
  !
  !===============================================================================

  USE data_parameters, ONLY :   &
       wp ! kind-type parameter for "normal" integer variables
  USE data_modelconfig, ONLY :   &
       ie,           & ! number of grid points in zonal direction
       je,           & ! number of grid points in meridional direction
       !ke,           & ! number of grid points in vertical direction
       ie_tot,       & ! total number of grid points in meridional direction
       je_tot,       & ! total number of grid points in meridional direction
       ke_tot !,       & ! total number of grid points in vertical direction
  USE data_parallel, ONLY :  &
       num_compute,     & ! number of compute PEs
       my_world_id !,     & ! rank of this subdomain in the global communicator
  USE parallel_utilities, ONLY : gather_field
  USE radar_utilities, ONLY : tolower

    IMPLICIT NONE

    CHARACTER(*), INTENT(in) :: cvar, levtyp, outdir
    CHARACTER(len=*), INTENT(in) :: header
    REAL(kind=wp), INTENT(in) :: feld(:,:,:)
    LOGICAL, INTENT(in) :: gzipflag
    INTEGER, INTENT(out) :: fehler

    INTEGER :: ke_loc
    INTEGER :: iunit, i, j, k, error
    CHARACTER(len=20) :: cie, cje, cke
    CHARACTER(len=cmaxlen) :: tmpdir
    CHARACTER(len=(cmaxlen+100)) :: dateiname
    CHARACTER(len=14) :: model_starttime, model_validtime
    REAL(kind=wp), ALLOCATABLE :: feld_tot(:,:)

  !------------ End of header ----------------------------------------------------


  !-------------------------------------------------------------------------------
  ! Begin Subroutine write_modelfield_param
  !-------------------------------------------------------------------------------

    fehler = 0

    ke_loc = SIZE(feld, 3)

    IF (my_world_id == 0) THEN

      iunit = 30

      tmpdir(:) = ' '
      IF (LEN_TRIM(outdir) > 0) THEN
        tmpdir = TRIM(ADJUSTL(outdir))//'/'
      END IF

      model_starttime = ydate_ini_mod
      model_validtime = get_datetime_act( l_round_to_minute=.TRUE. )
      dateiname = REPEAT(' ', len(dateiname))
      ! only use model times down to minutes
      dateiname = TRIM(tmpdir)//TRIM(tolower(cvar))//'_'//&
                  TRIM(model_starttime(1:12))//'_'//TRIM(model_validtime(1:12))//&
                  '_'//TRIM(levtyp)//'.dat'
      OPEN(iunit, file=TRIM(ADJUSTL(dateiname)), status='replace', iostat=error)
      IF (error /= 0) THEN
        fehler = 1
        WRITE (*,*) 'Error opening ', TRIM(ADJUSTL(dateiname)), ' !'
        RETURN
      END IF

      !WRITE(iunit,'(3A,3(I3,A))') '# size of ',cvar,' is (',&
      !      SIZE(feld,1),',',SIZE(feld,2),',',SIZE(feld,3),')'
      !WRITE(iunit,'(A,3(I3,A))') '# gathering to extent of (',ie_tot,',',je_tot,',',ke_loc,')'

      ! JM201125:
      ! This check is nonsense. Where did I get that from? feld is the ungathered subdomain field!
      ! We should probably check something, but I don't know what exactly...
      ! Or maybe we don't need and gather_field is taking care?
      !
      !IF ( (SIZE(feld, 1) /= ie_tot) .OR. (SIZE(feld, 2) /= je_tot) ) THEN
      !  OPEN(iunit, file=TRIM(ADJUSTL(dateiname)), status='replace', iostat=error)
      !  WRITE (iunit,*) 'Field has wrong extent. Aborting write.'
      !  CLOSE(iunit)
      !  RETURN
      !END IF

      cie = REPEAT(' ', LEN(cie))
      cje = REPEAT(' ', LEN(cje))
      cke = REPEAT(' ', LEN(cke))
      WRITE (cie, '(i20)') ie_tot
      WRITE (cje, '(i20)') je_tot
      WRITE (cke, '(i20)') ke_loc
      WRITE (iunit,'(a)') header
      WRITE (iunit,'(a,1x,a,1x,a)') TRIM(ADJUSTL(cie)), TRIM(ADJUSTL(cje)), TRIM(ADJUSTL(cke))
      CLOSE(iunit)
      OPEN(iunit, file=TRIM(ADJUSTL(dateiname)), status='old', position='append', iostat=error)
    ENDIF

    ALLOCATE(feld_tot(ie_tot, je_tot))
    feld_tot = -9999.99_wp

    DO k=1,ke_loc

      IF (num_compute > 1) THEN
        CALL gather_field(feld(:,:,k), ie, je, &
             feld_tot, ie_tot, je_tot, 0, error)
      ELSE
        feld_tot = feld(:,:,k)
      END IF

      IF (my_world_id == 0) THEN
        DO j=1,je_tot
          DO i=1,ie_tot
            IF (ABS(feld_tot(i,j)) < 1.0d-30) THEN
              WRITE (iunit, '(i1)') 0
            ELSE
              WRITE (iunit, '(es13.5)') feld_tot(i,j)
            END IF
          END DO
        END DO
      END IF

    END DO

    IF (my_world_id == 0) THEN

      CLOSE(iunit)

      IF (gzipflag) THEN
        CALL execute_command_line('gzip -f '//TRIM(ADJUSTL(dateiname)))
      END IF

    END IF

    IF (ALLOCATED(feld_tot)) DEALLOCATE (feld_tot)

   RETURN

  !-------------------------------------------------------------------------------
  ! End of module procedure write_modelfield_param
  !-------------------------------------------------------------------------------

!CONTAINS
!
!  FUNCTION tolower (c) RESULT (cl)
!
!    IMPLICIT NONE
!
!    CHARACTER(len=*) :: c
!    CHARACTER(len=LEN(c)) :: cl
!    INTEGER :: i, lacode, uacode, uzcode, icode
!
!    lacode=IACHAR('a')
!    uacode=IACHAR('A')
!    uzcode=IACHAR('Z')
!
!    cl = c
!    DO i = 1, LEN_TRIM(c)
!      icode = IACHAR(c(i:i))
!      IF (icode >= uacode .AND. icode <= uzcode) THEN
!        cl(i:i) = ACHAR(icode - uacode + lacode)
!      END IF
!    END DO
!
!  END FUNCTION tolower

END SUBROUTINE write_modelfield_param
#endif

END MODULE radar_mie_iface_cosmo_utils
