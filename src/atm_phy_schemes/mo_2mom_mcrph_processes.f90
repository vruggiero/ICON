! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!NEC$ options "-finline-max-depth=3 -finline-max-function-size=2000"

!
! Two-moment bulk microphysics after Seifert, Beheng and Blahak
!
! Description:
! Provides various modules and subroutines for two-moment bulk microphysics
!

MODULE mo_2mom_mcrph_processes

!===============================================================================!
! Re-write of sedimentation schemes 03/2019 by UB:
! - Technical re-write of sedi_icon_core() overtaken from COSMO src_twomom_sb.f90:
!   - sedi_icon_core() vectorized version, in principle reproducible on the CRAY
!     but only used #if defined (__SX__) || defined (__NEC_VH__) || defined (__NECSX__)
!   - scalar version sedi_icon_core() for scalar architectures
!   - sedi_icon_core_lwf() including liquid water fraction
! - New internal switches "lboxtracking=.true.", activating the new explicit
!   boxtracking sedimentation method from COSMO src_twomom_f90
!   (http://www.cosmo-model.org/content/model/documentation/core/docu_sedi_twomom.pdf)
!   instead of sedi_icon_core() et al.:
!   - sedi_icon_box_core() (vectorized version #if defined (__SX__) || defined (__NEC_VH__) || defined (__NECSX__))
!   - sedi_icon_box_core_lwf() (vectorized version #if defined (__SX__) || defined (__NEC_VH__) || defined (__NECSX__))
!===============================================================================!
! OpenACC compiler error workarounds (04/2023 by MJ):
! IPSF: Several intermediate pointers have been introduced to circumvent
!       segmentation faults with nvhpc 22.7. The affected lines are marked by
!       the following abbreviation:
!       ! ACCWA (nvhpc 22.7, IPSF, see above)
!       Without these pointer, the compiler or Nvidia runtime is otherwise
!       unable to find the derived type on the accelerator device.
!       This workaround also requires additional WAIT clauses.
!       04/2024: This bug also affects nvhpc 23.3
!===============================================================================!

  USE mo_kind,               ONLY: sp, wp
  USE mo_exception,          ONLY: finish, message, txt => message_text
  USE mo_math_constants,     ONLY: pi, pi4 => pi_4
  USE mo_physical_constants, ONLY: &
       & R_l   => rd,     & ! gas constant of dry air (luft)
       & R_d   => rv,     & ! gas constant of water vapor (dampf)
       & cp    => cpd,    & ! specific heat capacity of air at constant pressure
       & c_w   => clw,    & ! specific heat capacity of water
       & L_wd  => alv,    & ! specific heat of vaporization (wd: wasser->dampf)
       & L_ed  => als,    & ! specific heat of sublimation (ed: eis->dampf)
       & L_ew  => alf,    & ! specific heat of fusion (ew: eis->wasser)
       & T_3   => tmelt,  & ! melting temperature of ice
       & rho_w => rhoh2o, & ! density of liquid water
       & rho_ice => rhoice,&! density of pure ice
       & nu_l  => con_m,  & ! kinematic viscosity of air
       & D_v   => dv0,    & ! diffusivity of water vapor in air at 0 C
       & K_t   => con0_h, & ! heat conductivity of air
       & N_avo => avo,    & ! Avogadro number [1/mol]
       & k_b   => ak,     & ! Boltzmann constant [J/K]
       & grav               ! acceleration due to Earth's gravity
  USE mo_thdyn_functions, ONLY:     &
       & e_ws  => sat_pres_water,  & ! saturation pressure over liquid water
       & e_es  => sat_pres_ice       ! saturation pressure over ice
  USE mo_2mom_mcrph_types, ONLY: &
       & particle, particle_frozen, particle_lwf, atmosphere, &
       & particle_sphere, particle_rain_coeffs, particle_cloud_coeffs, aerosol_ccn, &
       & particle_ice_coeffs, particle_snow_coeffs, particle_graupel_coeffs, &
       & particle_coeffs, collection_coeffs, rain_riming_coeffs, dep_imm_coeffs, &
       & coll_coeffs_ir_pm, lookupt_1D, lookupt_4D
  USE mo_2mom_mcrph_setup, ONLY: &
       & particle_mass, particle_meanmass,  particle_diameter, particle_normdiameter, &
       & particle_velocity, particle_lwf_idx, &
       & particle_assign, particle_frozen_assign, particle_lwf_assign, &
       & coll_delta_11, coll_delta_12, coll_theta_11, coll_theta_12,   &
       & rain_mue_dm_relation, moment_gamma, n_f, n_sc
  
  USE mo_2mom_mcrph_config,         ONLY: t_cfg_2mom
  USE mo_2mom_mcrph_config_default, ONLY: cfg_2mom_default
  USE mo_2mom_mcrph_util, ONLY: &
       & rat2do3,                    &  ! rational function for lwf-melting scheme
       & dyn_visc_sutherland,        &  ! used in lwf melting scheme
       & Dv_Rasmussen,               &  ! used in lwf melting scheme
       & ka_Rasmussen,               &  ! used in lwf melting scheme
       & lh_evap_RH87,               &  ! used in lwf melting scheme
       & lh_melt_RH87,               &  ! used in lwf melting scheme
       & gamlookuptable,             &  ! For look-up table of incomplete Gamma function
       & incgfct_lower_lookup,       &  ! interpolation in table, lower incomplete Gamma function
       & incgfct_upper_lookup,       &  ! interpolation in talbe, upper incomplete Gamma function
       & dmin_wg_gr_ltab_equi,       &  ! For look-up table of wet growth diameter
       & dmin_wetgrowth_fun,         &  ! For 4d functional fit of wet growth diameter
       & luse_dmin_wetgrowth_table,  &
       & set_qnc,                    &
       & set_qni,                    &
       & set_qnr,                    &
       & set_qns,                    &
       & set_qng,                    &
       & set_qnh_expPSD_N0const,     &
       & estick_ltab_equi,           &
       & otab, tab, get_otab, equi_table

  USE mo_fortran_tools, ONLY: set_acc_host_or_device, assert_acc_device_only, init

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_processes'

  ! switches for ice scheme, ice nucleation, drop activation and autoconversion
  INTEGER  :: ice_typ, nuc_i_typ, nuc_c_typ, auto_typ

  ! Physical parameters and coefficients which occur only in the two-moment scheme

  ! .. lower limit of tke used in turbulent collision enhancement parameterization:
  REAL(wp), PARAMETER :: tke_min = 0.01_wp ** 2 ! [m^2 s^-2]
  
  ! .. some physical parameters not found in ICON
  REAL(wp), PARAMETER :: T_f     = 233.0_wp     !..below this temperature there is no liquid water

  ! .. for old saturation pressure relations (keep this for some time for testing)
  REAL(wp), PARAMETER :: A_e  = 2.18745584e1_wp !..Konst. Saettigungsdamppfdruck - Eis
  REAL(wp), PARAMETER :: A_w  = 1.72693882e1_wp !..Konst. Saettigungsdamppfdruck - Wasser
  REAL(wp), PARAMETER :: B_e  = 7.66000000e0_wp !..Konst. Saettigungsdamppfdruck - Eis
  REAL(wp), PARAMETER :: B_w  = 3.58600000e1_wp !..Konst. Saettigungsdamppfdruck - Wasser
  REAL(wp), PARAMETER :: e_3  = 6.10780000e2_wp !..Saettigungsdamppfdruck bei T = T_3

  ! .. Hallet-Mossop ice multiplication
  REAL(wp), PARAMETER ::           &
       &    C_mult     = 3.5e8_wp, &    !..Koeff. fuer Splintering
       &    T_mult_min = 265.0_wp, &    !..Minimale Temp. Splintering
       &    T_mult_max = 270.0_wp, &    !..Maximale Temp. Splintering
       &    T_mult_opt = 268.0_wp       !..Optimale Temp. Splintering

  ! .. Phillips et al. ice nucleation scheme, see ice_nucleation_homhet() for more details
  REAL(wp) ::                         &
       &    na_dust    = 160.e4_wp,   & ! initial number density of dust [1/m], Phillips08 value 162e3 (never used, reset later)
       &    na_soot    =  25.e6_wp,   & ! initial number density of soot [1/m], Phillips08 value 15e6 (never used, reset later)
       &    na_orga    =  30.e6_wp,   & ! initial number density of organics [1/m3], Phillips08 value 177e6 (never used, reset later)
       &    ni_het_max = 100.0e3_wp,  & ! max number of IN between 1-10 per liter, i.e. 1d3-10d3
       &    ni_hom_max = 5000.0e3_wp    ! number of liquid aerosols between 100-5000 per liter

  INTEGER, PARAMETER ::               & ! Look-up table for Phillips et al. nucleation
       &    ttmax  = 30,              & ! sets limit for temperature in look-up table
       &    ssmax  = 60,              & ! sets limit for ice supersaturation in look-up table
       &    ttstep = 2,               & ! increment for temperature in look-up table
       &    ssstep = 1                  ! increment for ice supersaturation in look-up table

  REAL(sp), DIMENSION(0:100,0:100)  :: &
       &    afrac_dust, &  ! look-up table of activated fraction of dust particles acting as ice nuclei
       &    afrac_soot, &  ! ... of soot particles
       &    afrac_orga     ! ... of organic material
  !$ACC DECLARE COPYIN(afrac_dust, afrac_soot, afrac_orga)

  INCLUDE 'phillips_nucleation_2010.incf'

  ! .. LWF melting scheme: coefficients for rational approximation functions
  REAL(wp) ::  &           
       &      adstarh(10), bdstarh(9), amelth(10), bmelth(9), aviwch(10), bviwch(9),  &
       &      avlwch(10),  bvlwch(9),  avnumh(10), bvnumh(9), convqh(4),  convnh(4),  &   
       &      adstarg(10), bdstarg(9), ameltg(10), bmeltg(9), aviwcg(10), bviwcg(9),  &
       &      avlwcg(10),  bvlwcg(9),  avnumg(10), bvnumg(9), convqg(4) , convng(4)   

  !..Include file with coefficients for lwf melting scheme
  INCLUDE 'hailcoeffs.incf'
  INCLUDE 'grplcoeffs.incf'

  ! Various parameters for collision and conversion rates
  REAL(wp), PARAMETER ::             &
       &    ecoll_min = 0.01_wp,     & ! min. eff. for graupel_cloud, ice_cloud and snow_cloud
       &    q_crit_ii = 1.000e-6_wp, & ! q-threshold for ice_selfcollection
       &    D_crit_ii = 5.0e-6_wp,   & ! D-threshold for ice_selfcollection  
       &    q_crit_r  = 1.000e-5_wp, & ! q-threshold for ice_rain_riming and snow_rain_riming
       &    D_crit_r  = 100.0e-6_wp, & ! D-threshold for ice_rain_riming and snow_rain_riming
       &    q_crit_fr = 1.000e-6_wp, & ! q-threshold for rain_freeze
       &    q_crit_c  = 1.000e-6_wp, & ! q-threshold for cloud water
       &    q_crit    = 1.000e-9_wp, & ! q-threshold elsewhere 1e-7 kg/m3 = 1e-4 g/m3 = 0.1 mg/m3
       &    D_conv_sg = 200.0e-6_wp, & ! D-threshold for conversion of snow to graupel
       &    D_conv_ig = 200.0e-6_wp, & ! D-threshold for conversion of ice to graupel 
       &    x_conv    = 0.100e-9_wp, & ! minimum mass of conversion due to riming
       &    D_crit_c  = 10.00e-6_wp, & ! D-threshold for cloud drop collection efficiency
       &    D_coll_c  = 40.00e-6_wp    ! upper bound for diameter in collision efficiency

  REAL(wp), PARAMETER ::           &
       &    T_nuc     = 268.15_wp, & ! lower temperature threshold for ice nucleation, -5 C
       &    T_freeze  = 273.15_wp    ! lower temperature threshold for raindrop freezing

  ! Parameter for evaporation of rain, determines change of n_rain during evaporation
  REAL(wp) :: rain_gfak   ! this is set in init_twomoment

  ! debug switches
  LOGICAL, PARAMETER     :: isdebug = .false.   ! use only when really desperate
  LOGICAL, PARAMETER     :: isprint = .true.    ! print-out initialization values
  
  ! some cloud microphysical switches
  LOGICAL, PARAMETER     :: ice_multiplication = .TRUE.  ! default is .true.
  LOGICAL, PARAMETER     :: enhanced_melting   = .TRUE.  ! default is .true.
  LOGICAL, PARAMETER     :: classic_melting_in_lwf_scheme = .False.

  REAL(wp), PARAMETER    :: pi6 = pi/6.0_wp, pi8 = pi/8.0_wp ! more pieces of pi

  TYPE(t_cfg_2mom) :: cfg_params !.. Container to hold some config params for the actual 2-mom call
  
  ! Parameters
  PUBLIC :: q_crit
  PUBLIC :: cfg_2mom_default, cfg_params
  ! Switches
  PUBLIC :: ice_typ, nuc_i_typ, nuc_c_typ, auto_typ
  PUBLIC :: isdebug, isprint
  ! Functions
  PUBLIC :: particle_assign, particle_frozen_assign, particle_lwf_assign
  ! Process Routines
  PUBLIC :: sedi_vel_rain, sedi_vel_sphere, sedi_vel_lwf, init_2mom_sedi_vel
  PUBLIC :: autoconversionSB, accretionSB, rain_selfcollectionSB
  PUBLIC :: autoconversionKB, accretionKB
  PUBLIC :: autoconversionKK, accretionKK
  PUBLIC :: rain_evaporation, evaporation
  PUBLIC :: cloud_freeze
  PUBLIC :: ice_nucleation_homhet
  PUBLIC :: vapor_dep_relaxation
  PUBLIC :: rain_freeze_gamlook
  PUBLIC :: ice_selfcollection
  PUBLIC :: snow_selfcollection
  PUBLIC :: snow_melting, ice_melting, graupel_melting, hail_melting_simple
  PUBLIC :: particle_melting_lwf,prepare_melting_lwf
  PUBLIC :: particle_particle_collection
  PUBLIC :: graupel_selfcollection
  PUBLIC :: particle_cloud_riming, particle_rain_riming
  PUBLIC :: graupel_hail_conv_wet_gamlook
  PUBLIC :: ice_riming, snow_riming
  PUBLIC :: ccn_activation_sk, ccn_activation_hdcp2, ccn_activation_sk_4d
  PUBLIC :: sedi_icon_rain, sedi_icon_sphere, sedi_icon_sphere_lwf
  PUBLIC :: set_default_n 

CONTAINS
  
  !********************************************************************************
  ! bulk sedimentation velocities
  !********************************************************************************
  
  SUBROUTINE sedi_vel_rain(this_in,thisCoeffs,q,x,rhocorr,vn,vq,its,ite,qc,lacc)
    CLASS(particle), INTENT(in), TARGET :: this_in
    CLASS(particle), POINTER :: this ! ACCWA (nvhpc 22.7, IPSF, see above)
    TYPE(particle_rain_coeffs), INTENT(in) :: thisCoeffs
    INTEGER,  INTENT(in)  :: its,ite
    REAL(wp), INTENT(in)  :: q(:),x(:), rhocorr(:)
    REAL(wp), INTENT(in), OPTIONAL  :: qc(:)
    REAL(wp), INTENT(inout) :: vn(:), vq(:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lacc
    
    INTEGER  :: i
    REAL(wp) :: D_m,mue,D_p
    LOGICAL :: lzacc
    
    CALL set_acc_host_or_device(lzacc, lacc)

    this => this_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) FIRSTPRIVATE(its, ite) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(D_m, mue, D_p)
    DO i=its,ite
      IF (q(i).GT.q_crit) THEN
        D_m = particle_diameter(this, x(i))
        
        IF (cfg_params%luse_mu_Dm_rain) THEN
          IF (PRESENT(qc)) THEN
            IF (qc(i) >= q_crit) THEN
              mue = (this%nu+1.0_wp)/this%b_geo - 1.0_wp
            ELSE
              mue = rain_mue_dm_relation(thisCoeffs, D_m)
            END IF
          ELSE
            mue = rain_mue_dm_relation(thisCoeffs, D_m)
          END IF
        ELSE
          mue = (this%nu+1.0_wp)/this%b_geo - 1.0_wp
        END IF
        
        D_p = D_m * EXP((-1./3.)*LOG((mue+3.)*(mue+2.)*(mue+1.)))
        vn(i) = thisCoeffs%alfa - thisCoeffs%beta * EXP(-(mue+1.)*LOG(1.0 + thisCoeffs%gama*D_p))
        vq(i) = thisCoeffs%alfa - thisCoeffs%beta * EXP(-(mue+4.)*LOG(1.0 + thisCoeffs%gama*D_p))
        vq(i) = MAX(vq(i),  this%vsedi_min)
        vn(i) = MAX(vn(i),  this%vsedi_min)
        vn(i) = vn(i) * rhocorr(i)
        vq(i) = vq(i) * rhocorr(i)
      ELSE
        vn(i) = 0.0_wp
        vq(i) = 0.0_wp
      END IF
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE sedi_vel_rain

  ! bulk sedimentation velocities
  SUBROUTINE sedi_vel_sphere(this,thisCoeffs,q,x,rhocorr,vn,vq,its,ite)
    CLASS(particle), INTENT(in)        :: this
    CLASS(particle_sphere), INTENT(in) :: thisCoeffs
    INTEGER,  INTENT(in)  :: its,ite
    REAL(wp), INTENT(in)  :: q(:),x(:),rhocorr(:)
    REAL(wp), INTENT(out) :: vn(:), vq(:)

    INTEGER  :: i
    REAL(wp) :: lam,v_n,v_q

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE(lam, v_n, v_q)
    DO i=its,ite
      IF (q(i).GT.q_crit) THEN
        lam = EXP(this%b_vel* LOG(thisCoeffs%coeff_lambda*x(i)))
        v_n = thisCoeffs%coeff_alfa_n * lam
        v_q = thisCoeffs%coeff_alfa_q * lam
        v_n = MAX(v_n,this%vsedi_min)
        v_q = MAX(v_q,this%vsedi_min)
        v_n = MIN(v_n,this%vsedi_max)
        v_q = MIN(v_q,this%vsedi_max)
        vn(i) = v_n * rhocorr(i)
        vq(i) = v_q * rhocorr(i)
      ELSE
        vn(i) = 0.0_wp
        vq(i) = 0.0_wp
      END IF
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE sedi_vel_sphere


  ! bulk sedimentation velocities
  SUBROUTINE sedi_vel_lwf(this,thisCoeffs,q,ql,x,rhocorr,vn,vq,vl,its,ite)
    CLASS(particle_lwf), INTENT(in)    :: this
    CLASS(particle_sphere), INTENT(in) :: thisCoeffs
    INTEGER,  INTENT(in)  :: its,ite
    REAL(wp), INTENT(in)  :: q(:),x(:),ql(:),rhocorr(:)
    REAL(wp), INTENT(out) :: vn(:),vq(:),vl(:)
    
    REAL(wp), DIMENSION(10) :: avq, avl, avn
    REAL(wp), DIMENSION(9)  :: bvq, bvl, bvn
    REAL(wp), PARAMETER     :: eps = 1e-20_wp
    
    INTEGER  :: i
    REAL(wp) :: v_n,v_q,v_l,D_p,D_n,lwf
    
    IF (this%name .EQ. 'hail_vivek') THEN
      avq = aviwch
      bvq = bviwch
      avl = avlwch
      bvl = bvlwch
      avn = avnumh
      bvn = bvnumh
    ELSEIF (this%name .EQ. 'graupel_vivek') THEN
      avq = aviwcg
      bvq = bviwcg
      avl = avlwcg
      bvl = bvlwcg
      avn = avnumg
      bvn = bvnumg
    ELSE
      CALL finish(TRIM(routine),'Error: unknown particle name in sedi_vel_lwf')      
    END IF
    
    DO i=its,ite
      IF (q(i).GT.q_crit) THEN
        lwf = ql(i)/(q(i)+eps)
        D_p = particle_diameter(this,x(i))
        D_n = particle_normdiameter(this,D_p)
        v_n = rat2do3(D_n,lwf,avn,bvn)
        v_q = rat2do3(D_n,lwf,avq,bvq)
        v_l = rat2do3(D_n,lwf,avl,bvl)
        v_n = MAX(v_n,this%vsedi_min)
        v_q = MAX(v_q,this%vsedi_min)
        v_l = MAX(v_l,this%vsedi_min)
        v_n = MIN(v_n,this%vsedi_max)
        v_q = MIN(v_q,this%vsedi_max)
        v_l = MIN(v_l,this%vsedi_max)
        vn(i) = v_n * rhocorr(i)
        vq(i) = v_q * rhocorr(i)
        vl(i) = v_l * rhocorr(i)
      ELSE
        vn(i) = 0.0_wp
        vq(i) = 0.0_wp
      END IF
    END DO
    
  END SUBROUTINE sedi_vel_lwf
  
  ! initialize coefficients for bulk sedimentation velocity
  SUBROUTINE init_2mom_sedi_vel(this,thisCoeffs)
    CLASS(particle), INTENT(in) :: this
    CLASS(particle_sphere), INTENT(out) :: thisCoeffs
    
    CHARACTER(len=*), PARAMETER :: sroutine = 'init_2mom_sedi_vel'
    
    thisCoeffs%coeff_alfa_n = this%a_vel * GAMMA((this%nu+this%b_vel+1.0)/this%mu) / GAMMA((this%nu+1.0)/this%mu)
    thisCoeffs%coeff_alfa_q = this%a_vel * GAMMA((this%nu+this%b_vel+2.0)/this%mu) / GAMMA((this%nu+2.0)/this%mu)
    thisCoeffs%coeff_lambda = GAMMA((this%nu+1.0)/this%mu)/GAMMA((this%nu+2.0)/this%mu)
    
    IF (isprint) THEN
      WRITE (txt,'(2A)') "    name  = ",this%name ; CALL message(sroutine,TRIM(txt))
      WRITE (txt,'(A,D14.7)') "    c_lam = ",thisCoeffs%coeff_lambda ; CALL message(sroutine,TRIM(txt))
      WRITE (txt,'(A,D14.7)') "    alf_n = ",thisCoeffs%coeff_alfa_n ; CALL message(sroutine,TRIM(txt))
      WRITE (txt,'(A,D14.7)') "    alf_q = ",thisCoeffs%coeff_alfa_q ; CALL message(sroutine,TRIM(txt))
    END IF
  END SUBROUTINE init_2mom_sedi_vel

  ! currently not used
  FUNCTION D_average_factor (parti)
    ! UB: Faktor zur Berechnung des mittleren Durchmessers von verallg. gammaverteilten Hydrometeoren:
    !     gueltig fuer D = a_geo * x^b_geo
    !     Berechnung des mittleren Durchmessers: D_average = parti%b_geo * D_average_factor * (q/qn)**parti%b_geo
    REAL(wp) :: D_average_factor
    CLASS(particle), INTENT(in) :: parti

    D_average_factor = &
         ( GAMMA( (parti%b_geo+parti%nu+1.0_wp)/parti%mu ) / &
           GAMMA( (parti%nu+1.0_wp)/parti%mu ) ) * &
         ( GAMMA( (parti%nu+1.0_wp)/parti%mu ) / GAMMA( (parti%nu+2.0_wp)/parti%mu ) ) ** parti%b_geo
  END FUNCTION D_average_factor

  !*******************************************************************************
  ! Fundamental physical relations
  !*******************************************************************************

  ! Molecular diffusivity of water vapor
  ELEMENTAL FUNCTION diffusivity(T,p) result(D_v)

    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(IN) :: T,p
    REAL(wp) :: D_v
    ! This is D_v = 8.7602e-5_wp * T_a**(1.81_wp) / p_a
    D_v = 8.7602e-5_wp * EXP(1.81_wp*LOG(T)) / p
    RETURN
  END FUNCTION diffusivity

  !*******************************************************************************
  ! saturation pressure over ice and liquid water                                *
  !*******************************************************************************

  ! ELEMENTAL REAL(wp) FUNCTION e_es (ta)
  !   !$ACC ROUTINE SEQ
  !   REAL(wp), INTENT(IN) :: ta
  !   e_es  = sat_pres_ice(ta)
  ! END FUNCTION e_es

  !  ELEMENTAL REAL(wp) FUNCTION e_ws (ta)
  !    REAL(wp), INTENT (IN) :: ta
  !    e_ws  = e_3 * EXP (A_w * (ta - T_3) / (ta - B_w))
  !  END FUNCTION e_ws_old

  ! FUNCTION e_ws_vec (ta,idim,jdim)
  !   INTEGER :: idim, jdim
  !   REAL(wp)               :: e_ws_vec(idim,jdim)
  !   REAL(wp), INTENT (IN)  :: ta(idim,jdim)
  !   e_ws_vec  = e_3 * EXP (A_w * (ta - T_3) / (ta - B_w))
  ! END FUNCTION e_ws_vec

  ! FUNCTION e_es_vec (ta,idim,jdim)
  !   INTEGER :: idim, jdim
  !   REAL(wp)               :: e_es_vec(idim,jdim)
  !   REAL(wp), INTENT (IN)  :: ta(idim,jdim)
  !   e_es_vec  = e_3 * EXP (A_e * (ta - T_3) / (ta - B_e))
  ! END FUNCTION e_es_vec

  SUBROUTINE autoconversionSB(ik_slice,dt,atmo,cloud_coeffs,cloud_in,rain)
    !*******************************************************************************
    ! Autoconversion of Seifert and Beheng (2001, Atmos. Res.)                     *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    CLASS(atmosphere), INTENT(in)   :: atmo
    TYPE(particle_cloud_coeffs), INTENT(in) :: cloud_coeffs
    CLASS(particle), INTENT(inout), TARGET :: cloud_in
    CLASS(particle), POINTER :: cloud ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle), INTENT(inout) :: rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER          :: i,k
    REAL(wp)         :: q_c, q_r, n_c, x_c, tau, phi, x_s_i, au, sc

    REAL(wp), PARAMETER :: k_1  = 6.00e+2_wp   !..Parameter for Phi
    REAL(wp), PARAMETER :: k_2  = 0.68e+0_wp   !..Parameter fof Phi
    REAL(wp), PARAMETER :: eps  = 1.00e-25_wp

   ! Onishi kernel (of 29 July 2015)                                                                                                          
    REAL(wp), PARAMETER :: kc1_a1 = 3.985e-03_wp
    REAL(wp), PARAMETER :: kc1_a2 = 6.210e-03_wp
    REAL(wp), PARAMETER :: kc1_a3 = 1.331e+00_wp
    REAL(wp), PARAMETER :: kc1_b1 = 1.381e+01_wp
    REAL(wp), PARAMETER :: kc1_b2 = 9.980e+00_wp
    REAL(wp), PARAMETER :: kc1_b3 = 5.018e-01_wp
    REAL(wp), PARAMETER :: kc1_c1 = 6.325e+00_wp
    REAL(wp), PARAMETER :: kc1_c2 = -9.238e-01_wp
    REAL(wp), PARAMETER :: kc1_c3 = -1.528e-01_wp
    REAL(wp), PARAMETER :: kc1_bet = 2.026e-03_wp

    REAL(wp)  :: kc_alf,kc_rad,kc_sig,kc_bet,prey,Re,tke,diss,k_turb,nu_c,D_c

    cloud => cloud_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    nu_c = cloud%nu

    kc_alf = ( kc1_a1 + kc1_a2 * nu_c )/ ( 1.0_wp + kc1_a3 * nu_c )
    kc_rad = ( kc1_b1 + kc1_b2 * nu_c )/ ( 1.0_wp + kc1_b3 * nu_c )
    kc_sig = ( kc1_c1 + kc1_c2 * nu_c )/ ( 1.0_wp + kc1_c3 * nu_c )
    kc_bet = kc1_bet
    prey   = -0.125_wp

    IF (isdebug) CALL message(routine, "autoconversionSB")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    x_s_i = 1.0_wp / cloud%x_max

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(q_c, n_c, q_r, x_c, au, tau, phi, sc, D_c, tke, diss, Re, k_turb)
    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          IF (q_c > q_crit) THEN

            n_c = cloud%n(i,k)
            q_r = rain%q(i,k)
            x_c = particle_meanmass(cloud, q_c,n_c)
   

            IF (cfg_params%lturb_enhc) THEN
               D_c = particle_diameter(cloud, x_c)
               tke = MAX(0.5_wp*(atmo%tke(i,k) + atmo%tke(i,k+1)), tke_min) ! from half levels to full levels
               ! Dissiplation rate according to Mellor-Yamada
               diss = MIN(3000.0_wp, EXP(1.5_wp*LOG(2.0_wp*tke))/(16.6_wp*cfg_params%turb_len) * 1.0E4_wp)
               Re = MAX(10000.0_wp * EXP ( 1.0_wp/6.0_wp*LOG(diss/100.0_wp) ), 2000.0_wp)
               k_turb = 1.0_wp + diss* EXP(prey*LOG(Re)) * (kc_bet + kc_alf * &
                        EXP( -1.0_wp* ((((D_c/2.0_wp)*1.e6_wp-kc_rad)/kc_sig)**2) ))
            ELSE
               k_turb = 1.0_wp
            END IF
            
            au  = cloud_coeffs%k_au * q_c**2 * x_c**2 * dt &
                 * cloud%rho_v(i,k)
            tau = MIN(MAX(1.0_wp-q_c/(q_c+q_r+eps),eps),0.9_wp)
            phi = k_1 * tau**k_2 * (1.0_wp - tau**k_2)**3
            au  = au * (1.0_wp + phi/(1.0_wp - tau)**2) * k_turb

            au  = MAX(MIN(q_c,au),0.0_wp)

            sc  = cloud_coeffs%k_sc * q_c**2 * dt * cloud%rho_v(i,k)

            rain%n(i,k)  = rain%n(i,k)  + au * x_s_i
            rain%q(i,k)  = rain%q(i,k)  + au
            cloud%n(i,k) = cloud%n(i,k) - MIN(n_c,sc)
            cloud%q(i,k) = cloud%q(i,k) - au

          ENDIF
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE autoconversionSB

  SUBROUTINE accretionSB(ik_slice, dt, atmo, cloud_in, rain_in)
    !*******************************************************************************
    ! Accretion of Seifert and Beheng (2001, Atmos. Res.)                          *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    CLASS(atmosphere), INTENT(in)   :: atmo
    CLASS(particle), INTENT(inout), TARGET :: cloud_in, rain_in
    CLASS(particle), POINTER :: cloud, rain ! ACCWA (nvhpc 22.7, IPSF, see above)
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER     :: i, k
    REAL(wp)    :: ac
    REAL(wp)    :: q_c, q_r, tau, phi, n_c, x_c,diss, tke, n_r, x_r, k_turb

    REAL(wp), PARAMETER :: k_r = 5.78_wp       ! kernel
    REAL(wp), PARAMETER :: k_1 = 5.00e-04_wp   ! Phi function
    REAL(wp), PARAMETER :: eps = 1.00e-25_wp

    IF (isdebug) CALL message(routine, "accretionSB")

    rain => rain_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    cloud => cloud_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(n_c, q_c, q_r, tau, phi, ac, x_c, diss, tke, n_r, x_r, k_turb)
    DO k = kstart,kend
       DO i = istart,iend

          n_c = cloud%n(i,k)
          q_c = cloud%q(i,k)
          q_r = rain%q(i,k)
          n_r = rain%n(i,k)
          ! Calculate dissipation rate assuming tur_len =300 (namelist) and B1=16.6 (same as in Mello-Yamada and in code)
          ! Approximation only valid when z >> tur_len and no stability corrections (which are small for large tke)
          IF (q_c > 0.0_wp.AND.q_r > 0.0_wp) THEN

             !..accretion rate of SB2001
             tau = MIN(MAX(1.0_wp-q_c/(q_c+q_r+eps),eps),1.0_wp)
             phi = (tau/(tau+k_1))**4
             x_r = particle_meanmass(rain, q_r,n_r)
             IF (cfg_params%lturb_enhc) THEN
               tke = MAX(0.5_wp*(atmo%tke(i,k) + atmo%tke(i,k+1)), tke_min) ! from half levels to full levels
               diss = MIN(3000.0_wp ,EXP(1.5_wp*LOG(2.0_wp*tke))/(16.6_wp*cfg_params%turb_len) * 1.0E4_wp)  ! Dissipation rate in cm-2 
               k_turb = 1.0_wp + 0.8e-3_wp * diss * (rain%x_min/x_r)**(2.0_wp/3.0_wp) ! Onishi-Seifert Kernel
             ELSE
               k_turb =1.0
             END IF
             ac  = k_r *  q_c * q_r * phi * dt * k_turb

             ac = MIN(q_c,ac)

             x_c = particle_meanmass(cloud, q_c,n_c)

             rain%q(i,k)  = rain%q(i,k)  + ac
             cloud%q(i,k) = cloud%q(i,k) - ac
             cloud%n(i,k) = cloud%n(i,k) - MIN(n_c,ac/x_c)
          ENDIF
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE accretionSB

  SUBROUTINE rain_selfcollectionSB(ik_slice, dt, atmo, rain_in)
    !*******************************************************************************
    ! Selfcollection of Seifert and Beheng (2001, Atmos. Res.)                     *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(inout), TARGET :: rain_in
    CLASS(particle), POINTER:: rain ! ACCWA (nvhpc 22.7, IPSF, see above)
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER        :: i, k
    REAL(wp)       :: sc, br
    REAL(wp)       :: q_r, n_r, x_r, d_r
    REAL(wp)       :: tke, diss, k_turb

    !..Parameters based on Seifert (2008, JAS)
    REAL(wp), PARAMETER :: D_br = 1.10e-3_wp
    REAL(wp), PARAMETER :: k_rr = 4.33e+0_wp
    REAL(wp), PARAMETER :: k_br = 1.00e+3_wp

    IF (isdebug) CALL message(routine, "rain_selfcollectionSB")

    rain => rain_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(n_r, q_r, x_r, D_r, sc, br, tke, diss, k_turb)
    DO k = kstart,kend
       DO i = istart,iend

          n_r = rain%n(i,k)
          q_r = rain%q(i,k)
          IF (q_r > 0.0_wp) THEN
            x_r = particle_meanmass(rain, q_r,n_r)
            D_r = particle_diameter(rain, x_r)
            IF (cfg_params%lturb_enhc) THEN
              tke = MAX(0.5_wp*(atmo%tke(i,k) + atmo%tke(i,k+1)), tke_min) ! from half levels to full levels
              diss = MIN(3000.0_wp ,EXP(1.5_wp*LOG(2.0_wp*tke))/(16.6_wp*cfg_params%turb_len) * 1.0E4_wp)  ! Dissipation rate in cm-2 
              k_turb = 1.0_wp + 0.8e-3_wp * diss * (rain%x_min/x_r)**(2.0_wp/3.0_wp)
            ELSE
              k_turb = 1.0_wp
            END IF
            !..Selfcollection as in SB2001
            sc = k_rr *  n_r * q_r * rain%rho_v(i,k) * dt * k_turb 

            !..Breakup as in Seifert (2008, JAS), Eq. (A13)
            br = 0.0_wp
            IF (D_r.GT.0.30e-3_wp) THEN
              br = MIN(k_br * (D_r - D_br) + 1.0_wp, 1.05_wp)  * sc ! Limit of rain breakup from Saleeby et al. JAS 2022
            ENDIF

            rain%n(i,k) = n_r - MIN(n_r,sc-br)
         ENDIF
      END DO
   END DO
   !$ACC END PARALLEL
   !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE rain_selfcollectionSB

  SUBROUTINE autoconversionKB(ik_slice, dt, cloud, rain)
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    CLASS(particle), INTENT(inout) :: cloud, rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER          :: i,k
    REAL(wp)         :: q_c, x_c, n_c, x_s_i, au
    REAL(wp), PARAMETER :: nu_c = 9.59_wp
    REAL(wp), PARAMETER :: k_a  = 6.0e+25_wp * nu_c**(-1.7_wp)

#ifdef _OPENACC
    CALL finish("autoconversionKB", "Routine has not been ported to openACC yet.")
#endif

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    x_s_i = 1.0_wp / cloud%x_max

    !..Parameterization of Beheng (1994)
    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          n_c = cloud%n(i,k)
          IF (q_c > q_crit) THEN

             x_c = particle_meanmass(cloud, q_c,n_c)

             !..Berechnung der Autokonversionsrate nach Beheng (1994)
             au = k_a * (x_c*1e3_wp )**(3.3_wp) * (q_c*1e-3)**(1.4_wp) * dt * 1e3_wp
             au = MIN(q_c,au)

             rain%n(i,k)  = rain%n(i,k)  + au * x_s_i
             rain%q(i,k)  = rain%q(i,k)  + au
             cloud%n(i,k) = cloud%n(i,k) - au * x_s_i * 2.0_wp
             cloud%q(i,k) = cloud%q(i,k) - au

          ENDIF
       END DO
    END DO

  END SUBROUTINE autoconversionKB

  SUBROUTINE accretionKB(ik_slice, dt, cloud, rain)
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    CLASS(particle), INTENT(inout)   :: cloud, rain

    !..Parameter of Beheng (1994)
    REAL(wp), PARAMETER :: k_r = 6.00d+00
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER      :: i,k
    REAL(wp)     :: ac
    REAL(wp)     :: q_c, q_r

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          q_r = rain%q(i,k)

          IF (q_c > q_crit .and. q_r > q_crit) THEN

             ac = k_r *  q_c * q_r * dt
             ac = MIN(q_c,ac)

             rain%q(i,k)  = rain%q(i,k)  + ac
             cloud%q(i,k) = cloud%q(i,k) - ac

          ENDIF
       END DO
    END DO

  END SUBROUTINE accretionKB

  SUBROUTINE autoconversionKK(ik_slice, dt, cloud,rain)

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    CLASS(particle), INTENT(inout) :: cloud, rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_c, x_s_i, au, n_c
    REAL(wp), PARAMETER :: k_a  = 1350.0_wp

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    x_s_i  = 1.0_wp / cloud%x_max

    !..Parameterization of Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          IF (q_c > q_crit) THEN

             n_c = cloud%n(i,k) * 1e6_wp  ! in 1/cm3

             !..autoconversion rate of KK2000
             au = k_a * q_c**(2.47_wp) * n_c**(-1.79_wp) * dt
             au = MIN(q_c,au)

             rain%n(i,k)  = rain%n(i,k)  + au * x_s_i
             rain%q(i,k)  = rain%q(i,k)  + au
             cloud%n(i,k) = cloud%n(i,k) - au * x_s_i * 2.0_wp
             cloud%q(i,k) = cloud%q(i,k) - au

          ENDIF
       END DO
    END DO

  END SUBROUTINE autoconversionKK

  SUBROUTINE accretionKK(ik_slice, dt, cloud, rain)
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    CLASS(particle), INTENT(inout) :: cloud, rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: ac, q_c, q_r
    REAL(wp), PARAMETER :: k_a = 6.70d+01

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !..Parameterization of Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          q_r = rain%q(i,k)

          IF (q_c > q_crit .AND. q_r > q_crit) THEN

             ac = k_a *  (q_c * q_r)**1.15 * dt * 1e3_wp
             ac = MIN(q_c,ac)

             rain%q(i,k)  = rain%q(i,k)  + ac
             cloud%q(i,k) = cloud%q(i,k) - ac
          ENDIF
       END DO
    END DO

  END SUBROUTINE accretionKK

  SUBROUTINE rain_evaporation(ik_slice, dt, rain_coeffs, rain_gfak, atmo, cloud, rain_in)
    !*******************************************************************************
    ! Evaporation of rain based on Seifert (2008, J. Atmos. Sci.)                  *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    
    REAL(wp), INTENT(in) :: rain_gfak   ! this is set in init_twomoment
    TYPE(particle_rain_coeffs), INTENT(in) :: rain_coeffs

    ! 2mom variables
    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle),  INTENT(in)    :: cloud
    CLASS(particle),  INTENT(inout), TARGET :: rain_in
    CLASS(particle), POINTER :: rain ! ACCWA (nvhpc 22.7, IPSF, see above)

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend

    ! local variables
    INTEGER             :: i,k
    REAL(wp)            :: T_a,p_a,e_sw,s_sw,g_d,eva_q,eva_n,eva_q_fak,vm
    REAL(wp)            :: q_r,n_r,x_r,e_d,f_v
    REAL(wp)            :: mue,d_m,gamma_eva,lam,d_vtp,gfak,mue6
    REAL(wp)            :: aa,bb,cc
    REAL(wp), PARAMETER :: D_br = 1.1e-3_wp

    REAL(wp), PARAMETER :: eva_q_fak_low        = 0.3_wp  ! \  Parameters of the
    REAL(wp), PARAMETER :: eva_q_fak_high       = 1.0_wp  !  Parameters of the
    REAL(wp), PARAMETER :: eva_q_fak_Dbr_minfak = 0.75_wp ! |  ramp-function eva_q_fak(D_m) for reduction
    REAL(wp), PARAMETER :: eva_q_fak_Dbr_maxfak = 0.9_wp  ! /  of evaporation of drizzle-like rain         

    LOGICAL, PARAMETER   :: reduce_evaporation = .false.

    rain => rain_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    aa = rain_coeffs%alfa
    bb = rain_coeffs%beta
    cc = rain_coeffs%gama

    IF (isdebug) CALL message(routine, "rain_evaporation")

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) FIRSTPRIVATE(kstart, kend, istart, iend, dt, rain_gfak)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC   PRIVATE(T_a, p_a, e_sw, s_sw, g_d, eva_q, eva_n, eva_q_fak, q_r, n_r, x_r, e_d, f_v) &
    !$ACC   PRIVATE(mue, D_m, gamma_eva, lam, D_vtp, gfak, vm, mue6)
    DO k = kstart,kend
       DO i = istart,iend

          q_r  = rain%q(i,k)
          n_r  = rain%n(i,k)
          p_a  = atmo%p(i,k)
          T_a  = atmo%T(i,k)

          e_d  = atmo%qv(i,k) * R_d * T_a
          e_sw = e_ws(T_a)
          s_sw = e_d / e_sw - 1.0_wp

          IF (s_sw < 0.0_wp .AND. q_r > 0.0_wp .AND. cloud%q(i,k) < q_crit) THEN

             D_vtp = diffusivity(T_a,p_a)

             ! note that 2*pi is correct, because c_r = 1/2 is assumed
             g_d = 2.0_wp*pi / ( L_wd**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_sw) )

             x_r = particle_meanmass(rain, q_r,n_r)
             D_m = particle_diameter(rain, x_r)

             ! Eq. (20) of Seifert (2008)
             IF (cfg_params%luse_mu_Dm_rain) THEN
               mue = rain_mue_dm_relation(rain_coeffs,D_m)
             ELSE
               mue = (rain%nu+1.0_wp)/rain%b_geo - 1.0_wp
             END IF
                
             ! Eq. (A8)
             lam = exp(1.0_wp/3.0_wp*log(pi6*rho_w*(mue+3.0_wp)*(mue+2.0_wp)*(mue+1.0_wp)/x_r))

             ! chebyshev approximation of Gamma(mue+5/2)/Gamma(mue+2)
             gfak =  0.1357940435E+01_wp &
                  &  + mue * ( +0.3033273220E+00_wp  &
                  &  + mue * ( -0.1299313363E-01_wp  &
                  &  + mue * ( +0.4002257774E-03_wp  &
                  &  - mue * 0.4856703981E-05_wp ) ) )
             
             ! Mean velocity using the mean diameter: (mue+1)/D
             ! This is more exact than the exponential (Chebychev) approximation for small droplets,
             ! and not too bad for large droplets
             vm = MAX(aa - bb * EXP (-cc*D_m ), 0.0_wp)
             ! Ventilation coeffcient using an averaged velocity
             f_v  = rain%a_ven + rain%b_ven * N_sc**n_f * gfak                     &
                  &            * SQRT(vm/nu_l * rain%rho_v(i,k) / lam)

             !! This has been the previous Chebychev Approximation (kept for reference):
             ! mm = mue+5.0/2.0

             ! Eq. (A7) rewritten with (A5) and (A9)
             !f_v  = rain%a_ven + rain%b_ven * N_sc**n_f * gfak                     &
             !     &       * SQRT(aa/nu_l * rain%rho_v(i,k) / lam)                  &
             !     &    * (1.0 - 1./2.  * (bb/aa)**1 * EXP(mm*LOG(lam/(1.*cc+lam))) &
             !     &           - 1./8.  * (bb/aa)**2 * EXP(mm*LOG(lam/(2.*cc+lam))) &
             !     &           - 1./16. * (bb/aa)**3 * EXP(mm*LOG(lam/(3.*cc+lam))) &
             !     &           - 5./127.* (bb/aa)**4 * EXP(mm*LOG(lam/(4.*cc+lam))) )

             IF (rain_gfak > 0) THEN
               ! Eq. (20) of Seifert (2008) for fixed parameters (needed for rain_gfak)
               ! The reason of fixing cmu0=6.0 is that the relation can provide strange values for other mue values
               IF (D_m <= 1.1e-3_wp) THEN
                 mue6 = 6.0_wp*TANH((4.0e3_wp*(D_m-1.1e-3_wp))**2) + 1.0_wp
               ELSE
                 mue6 = 30.0_wp*TANH((1.0e3_wp*(D_m-1.1e-3_wp))**2) + 1.0_wp
               ENDIF
                ! Eq. (23)
               gamma_eva = rain_gfak * (D_br/D_m) * EXP(-0.2_wp*mue6)
             ELSE
               gamma_eva = 1.0_wp
             END IF

             ! Eq (A5) with (A9) and Gamma(mue+2) from (A7)
             eva_q = g_d * n_r * (mue+1.0_wp) / lam * f_v * s_sw * dt

             ! UB: empirical correction factor to reduce evaporation of drizzle-like rain (D_m < D_br):
             IF ( reduce_evaporation ) THEN
               eva_q_fak = MIN( MAX( ( eva_q_fak_low + &
                    (eva_q_fak_high - eva_q_fak_low) / &
                    (eva_q_fak_Dbr_maxfak*D_br - eva_q_fak_Dbr_minfak*D_br) * &
                    (D_m - eva_q_fak_Dbr_minfak*D_br) ), eva_q_fak_low), eva_q_fak_high)
               eva_q = eva_q_fak * eva_q
             END IF

             eva_n = gamma_eva * eva_q / x_r

             eva_q = MAX(-eva_q,0.0_wp)
             eva_n = MAX(-eva_n,0.0_wp)

             eva_q = MIN(eva_q,q_r)
             eva_n = MIN(eva_n,n_r)

             atmo%qv(i,k) = atmo%qv(i,k) + eva_q
             rain%q(i,k)  = rain%q(i,k)  - eva_q
             rain%n(i,k)  = rain%n(i,k)  - eva_n
          END IF
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE rain_evaporation

  SUBROUTINE evaporation(ik_slice, dt, atmo, prtcl_in, coeffs)
    !*******************************************************************************
    ! Evaporation of melting snow/graupel/hail, see SB2006                                      *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(inout), TARGET :: prtcl_in
    CLASS(particle), POINTER :: prtcl ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle_coeffs), INTENT(in) :: coeffs

    LOGICAL, PARAMETER  :: reduce_melting = .true.

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a,e_sw,s_sw,g_d,eva_q,eva_n
    REAL(wp)            :: q,n,x,d,v,f_v,e_d

    prtcl => prtcl_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(q, T_a, e_d, e_sw, s_sw, g_d, n, x, D, v, f_v) &
    !$ACC   PRIVATE(eva_q, eva_n)
    DO k = kstart,kend
       DO i = istart,iend

          q = prtcl%q(i,k)
          T_a = atmo%T(i,k)

          IF (q > 0.0_wp .AND. T_a > T_3) THEN

            e_d  = atmo%qv(i,k) * R_d * T_a
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0_wp

            !.. Eq. (37) of SB2006, note that 4*pi is correct because c is used below
            g_d = 4.0_wp*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n = prtcl%n(i,k)
            x = particle_meanmass(prtcl, q,n)
            D = particle_diameter(prtcl, x)
            v = particle_velocity(prtcl, x) * prtcl%rho_v(i,k)

            f_v  = coeffs%a_f + coeffs%b_f * sqrt(v*D)

            eva_q = g_d * n * coeffs%c_i * d * f_v * s_sw * dt

            eva_q = MAX(-eva_q,0.0_wp) 

            !.. Complete evaporation of some of the melting frozen particles: parameterized in a way
            !   to conserve the mean mass, similar to the case "gamma_eva = 1.0" in rain_evaporation() above:
            IF (reduce_melting) THEN
              eva_n = eva_q / x
              eva_n = MIN(eva_n,n)
              prtcl%n(i,k) = prtcl%n(i,k) - eva_n
            END IF

            eva_q = MIN(eva_q,q)

            atmo%qv(i,k) = atmo%qv(i,k) + eva_q
            prtcl%q(i,k) = prtcl%q(i,k) - eva_q

          END IF
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE evaporation

  SUBROUTINE cloud_freeze(ik_slice, dt, cloud_coeffs, qnc_const, atmo, cloud_in, ice)
    !*******************************************************************************
    ! This is only the homogeneous freezing of liquid water droplets.              *
    ! Immersion freezing and homogeneous freezing of liquid aerosols are           *
    ! treated in the subroutine ice_nucleation_homhet()                            *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    REAL(wp), INTENT(in) :: qnc_const
    TYPE(particle_cloud_coeffs), INTENT(in) :: cloud_coeffs
    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(inout), TARGET :: cloud_in
    CLASS(particle), POINTER :: cloud ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle), INTENT(inout) :: ice
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER            :: i, k
    REAL(wp)           :: fr_q, fr_n, T_a, q_c, x_c, n_c, j_hom, T_c

    REAL(wp), PARAMETER :: log_10 = LOG(10.0_wp)

    cloud => cloud_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, T_c, q_c, n_c, fr_q, fr_n, x_c, j_hom)
    DO k = kstart,kend
      DO i = istart,iend

        T_a = atmo%T(i,k)
        IF (T_a < T_3) THEN

          T_c = T_a - T_3
          q_c = cloud%q(i,k)
          n_c = cloud%n(i,k)
          IF (q_c > 0.0_wp .and. T_c < -30.0_wp) THEN
            IF (T_c < -50.0_wp) THEN
              fr_q = q_c             !..instantaneous freezing
              fr_n = n_c             !..below -50 C
            ELSE
              x_c = particle_meanmass(cloud, q_c, n_c)

              !..Hom. freezing based on Jeffrey und Austin (1997), see also Cotton und Field (2001)
              !  (note that log in Cotton and Field is log10, not ln)
              IF (T_c > -30.0_wp) THEN
!                 j_hom = 1.0e6_wp/rho_w * 10**(-7.63-2.996*(T_c+30.0))           !..J in 1/(kg s)
                 j_hom = 1.0e6_wp/rho_w * EXP((-7.63_wp-2.996_wp*(T_c+30.0_wp))*log_10)
              ELSE
!                 j_hom = 1.0e6_wp/rho_w &
!                      &  * 10**(-243.4-14.75*T_c-0.307*T_c**2-0.00287*T_c**3-0.0000102*T_c**4)
                 j_hom = 1.0e6_wp/rho_w &
                      &  * EXP((-243.4_wp-14.75_wp*T_c-0.307_wp*T_c**2-0.00287_wp*T_c**3-0.0000102_wp*T_c**4)*log_10)
              ENDIF

              fr_n  = j_hom * q_c *  dt
              fr_q  = j_hom * q_c * x_c * dt *  cloud_coeffs%c_z
              fr_q  = MIN(fr_q,q_c)
              fr_n  = MIN(fr_n,n_c)
            END IF

            cloud%q(i,k) = cloud%q(i,k) - fr_q
            cloud%n(i,k) = cloud%n(i,k) - fr_n

            fr_n  = MAX(fr_n, fr_q/cloud%x_max)

            !..special treatment for constant drop number
            IF (nuc_c_typ .EQ. 0) THEN
              ! ... force upper bound in cloud_freeze'
              fr_n = MAX(MIN(fr_n, qnc_const-ice%n(i,k)), 0.0_wp)
            ENDIF

            ice%q(i,k)   = ice%q(i,k) + fr_q
            ice%n(i,k)   = ice%n(i,k) + fr_n
          ENDIF
        END IF
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE cloud_freeze

  SUBROUTINE ice_nucleation_homhet(ik_slice, use_prog_in, &
       atmo, cloud, ice_in, n_inact, n_inpot)
    !*******************************************************************************
    !                                                                              *
    ! Homogeneous and heterogeneous ice nucleation                                 *
    !                                                                              *
    ! Nucleation scheme is based on the papers:                                    *
    !                                                                              *
    ! "A parametrization of cirrus cloud formation: Homogenous                     *
    ! freezing of supercooled aerosols" by B. Kaercher and                         *
    ! U. Lohmann 2002 (KL02 hereafter)                                             *
    !                                                                              *
    ! "Physically based parameterization of cirrus cloud formation                 *
    ! for use in global atmospheric models" by B. Kaercher, J. Hendricks           *
    ! and U. Lohmann 2006 (KHL06 hereafter)                                        *
    !                                                                              *
    ! and Phillips et al. (2008) with extensions                                   *
    !                                                                              *
    ! implementation by Carmen Koehler and AS                                      *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    LOGICAL, INTENT(in) :: use_prog_in

    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(inout), TARGET :: ice_in
    CLASS(particle), POINTER :: ice ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle), INTENT(inout) :: cloud
    REAL(wp), DIMENSION(:,:) :: n_inact
    REAL(wp), DIMENSION(:,:), OPTIONAL :: n_inpot

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER              :: i,k,nuc_typ
    REAL(wp)             :: nuc_n, nuc_q
    REAL(wp)             :: T_a, p_a, ssi
    REAL(wp)             :: q_i,n_i,x_i,r_i

    ! switch for version of Phillips et al. scheme
    ! (but make sure you have the correct INCLUDE file)
    INTEGER, PARAMETER :: iphillips = 2010

    ! switch for Hande et al. ice nucleation, if .true. this turns off Phillips scheme
    LOGICAL  :: use_hdcp2_het

    ! some more constants needed for homogeneous nucleation scheme
    REAL(wp), PARAMETER ::            &
         r_0     = 0.25e-6_wp          , &    ! aerosol particle radius prior to freezing
         alpha_d = 0.5_wp              , &    ! deposition coefficient (KL02; Spichtinger & Gierens 2009)
         M_w     = 18.01528e-3_wp      , &    ! molecular mass of water [kg/mol]
         M_a     = 28.96e-3_wp         , &    ! molecular mass of air [kg/mol]
         ma_w    = M_w / N_avo         , &    ! mass of water molecule [kg]
         svol    = ma_w / rho_ice             ! specific volume of a water molecule in ice

    REAL(wp)  :: e_si
    REAL(wp)  :: ni_hom,ri_hom,mi_hom
    REAL(wp)  :: v_th,n_sat,flux,phi,cool,tau,delta,w_pre,scr
    REAL(wp)  :: ctau, acoeff(3),bcoeff(2), ri_dot
    REAL(wp)  :: kappa,sqrtkap,ren,R_imfc,R_im,R_ik,ri_0

    ! parameters for Hande et al. nucleation parameterization for HDCP2 simulations
    TYPE(dep_imm_coeffs), PARAMETER :: hdcp2_nuc_coeffs(5) &
         = (/ dep_imm_coeffs(   &  ! Spring of Table 1
         nin_imm = 1.5684e5_wp, &
         alf_imm = 0.2466_wp,   &
         bet_imm = 1.2293_wp,   &
         nin_dep = 1.7836e5_wp, &
         alf_dep = 0.0075_wp,   &
         bet_dep = 2.0341_wp),  &
         dep_imm_coeffs(        &  ! Summer
         nin_imm = 2.9694e4_wp, &
         alf_imm = 0.2813_wp,   &
         bet_imm = 1.1778_wp,   &
         nin_dep = 2.6543e4_wp, &
         alf_dep = 0.0020_wp,   &
         bet_dep = 2.5128_wp),  &
         dep_imm_coeffs(        &  ! Autumn
         nin_imm = 4.9920e4_wp, &
         alf_imm = 0.2622_wp,   &
         bet_imm = 1.2044_wp,   &
         nin_dep = 7.7167e4_wp, &
         alf_dep = 0.0406_wp,   &
         bet_dep = 1.4705_wp),  &
         dep_imm_coeffs(        &  ! Winter
         nin_imm = 1.0259e5_wp, &
         alf_imm = 0.2073_wp,   &
         bet_imm = 1.2873_wp,   &
         nin_dep = 1.1663e4_wp, &
         alf_dep = 0.0194_wp,   &
         bet_dep = 1.6943),     &
         dep_imm_coeffs(        &  ! Spring with 95th percentile scaling factor
         nin_imm = 1.5684e5_wp * 17.82_wp, &
         alf_imm = 0.2466_wp,   &
         bet_imm = 1.2293_wp,   &
         nin_dep = 1.7836e5_wp * 5.87_wp, &
         alf_dep = 0.0075_wp,   &
         bet_dep = 2.0341_wp) /)

    LOGICAL   :: use_homnuc

    LOGICAL  :: ndiag_mask(ik_slice(1):ik_slice(2), ik_slice(3):ik_slice(4))
    REAL(wp) :: nuc_n_a(ik_slice(1):ik_slice(2), ik_slice(3):ik_slice(4))

    !$ACC DATA CREATE(nuc_n_a, ndiag_mask, acoeff, bcoeff)

    ice => ice_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    nuc_typ = nuc_i_typ

    SELECT CASE (nuc_typ)
    CASE(0)
      ! Heterogeneous nucleation ONLY
      use_homnuc = .FALSE.
    CASE(1:9)
      ! Homog. and het. nucleation"
      use_homnuc = .TRUE.
    END SELECT

    IF (isdebug) THEN
      IF (.NOT.use_homnuc) THEN
        WRITE(txt,*) "ice_nucleation_homhet: Heterogeneous nucleation only"
      ELSE
        WRITE(txt,*) "ice_nucleation_homhet: Homogeneous and heterogeneous nucleation"
      END IF
      CALL message(routine,TRIM(txt))
    END IF

    ! switch for Hande et al. ice nucleation, if .true. this turns off Phillips scheme
    use_hdcp2_het = (nuc_typ.le.5)

    ! Heterogeneous nucleation using Hande et al. scheme
    IF (use_hdcp2_het) THEN
#ifdef _OPENACC
      CALL finish('mo_2mom_mcrph_processes:','ice_nucleation_het_hdcp2 not available on GPU')
#endif
      IF (nuc_typ < 1 .OR. nuc_typ > 5) THEN
        CALL finish(TRIM(routine), &
             & 'Error in two_moment_mcrph: Invalid value nuc_typ in case of &
             &use_hdcp2_het=.true.')
      END IF
      CALL ice_nucleation_het_hdcp2(ik_slice, atmo, ice, cloud, &
           use_prog_in, hdcp2_nuc_coeffs(nuc_typ), n_inact, ndiag_mask, nuc_n_a)
    ELSE
      ! Heterogeneous nucleation using Phillips et al. scheme
      IF (iphillips == 2010) THEN
        ! possible pre-defined choices
        IF (nuc_typ.EQ.6) THEN  ! with no organics and rather high soot, coming close to Meyers formula at -20 C
          na_dust  = 160.e4_wp    ! initial number density of dust [1/m3]
          na_soot  =  30.e6_wp    ! initial number density of soot [1/m3]
          na_orga  =   0.e0_wp    ! initial number density of organics [1/m3]
        ELSEIF (nuc_typ.EQ.7) THEN     ! with some organics and rather high soot,
          na_dust  = 160.e4_wp    !          coming close to Meyers formula at -20 C
          na_soot  =  25.e6_wp
          na_orga  =  30.e6_wp
        ELSE IF (nuc_typ.EQ.8) THEN     ! no organics, no soot, coming close to DeMott et al. 2010 at -20 C
          na_dust  =  70.e4_wp    ! i.e. roughly one order in magnitude lower than Meyers
          na_soot  =   0.e6_wp
          na_orga  =   0.e6_wp
        ELSE
          CALL finish(TRIM(routine),&
               & 'Error in two_moment_mcrph: Invalid value nuc_typ in case of use_hdcp2_het=.false.')
        END IF
      END IF
      CALL ice_nucleation_het_philips(ik_slice, atmo, ice, cloud, &
           use_prog_in, n_inact, ndiag_mask, nuc_n_a, n_inpot)
    END IF

    IF (use_prog_in) THEN
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k = kstart, kend
        DO i = istart, iend
          n_inpot(i,k) = MERGE(MAX(n_inpot(i,k) - nuc_n_a(i, k), 0.0_wp), &
               &               n_inpot(i, k), &
               &               ndiag_mask(i, k))
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    ! Homogeneous nucleation using KHL06 approach
    IF (use_homnuc) THEN
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC   PRIVATE(acoeff, bcoeff, p_a, T_a, e_si, ssi, scr, n_i, q_i, x_i, r_i) &
      !$ACC   PRIVATE(v_th, flux, n_sat, ri_dot, R_ik, w_pre, cool, ctau, tau, delta, phi) &
      !$ACC   PRIVATE(kappa, sqrtkap, ren, R_imfc, R_im, ni_hom, ri_0, ri_hom, mi_hom, nuc_n, nuc_q)
      DO k = kstart,kend
        DO i = istart,iend
          p_a  = atmo%p(i,k)
          T_a  = atmo%T(i,k)
          e_si = e_es(T_a)
          ssi  = atmo%qv(i,k) * R_d * T_a / e_si

          ! critical supersaturation for homogeneous nucleation
          scr  = 2.349 - T_a * (1.0_wp/ 259.00_wp)

          IF (ssi > scr .AND. T_a < 235.0 .AND. ice%n(i,k) < ni_hom_max ) THEN

            n_i = ice%n(i,k)
            q_i = ice%q(i,k)
            x_i = particle_meanmass(ice, q_i,n_i)
!            r_i = (x_i/(4./3.*pi*rho_ice))**(1./3.)
            r_i = EXP( (1.0_wp/3.0_wp)*LOG(x_i/(4.0_wp/3.0_wp*pi*rho_ice)) )

            v_th  = SQRT( 8.0_wp*k_b*T_a/(pi*ma_w) )
            flux  = alpha_d * v_th/4.
            n_sat = e_si / (k_b*T_a)

            ! coeffs of supersaturation equation
            acoeff(1) = (L_ed * grav) / (cp * R_d * T_a**2) - grav/(R_l * T_a)
            acoeff(2) = 1.0_wp/n_sat
            acoeff(3) = (L_ed**2 * M_w * ma_w)/(cp * p_a * T_a * M_a)

            ! coeffs of depositional growth equation
            bcoeff(1) = flux * svol * n_sat * (ssi - 1.0_wp)
            bcoeff(2) = flux / diffusivity(T_a,p_a)

            ! pre-existing ice crystals included as reduced updraft speed
            ri_dot = bcoeff(1) / (1.0_wp + bcoeff(2) * r_i)
            R_ik   = (4.0_wp * pi) / svol * n_i * r_i**2 * ri_dot
            w_pre  = (acoeff(2) + acoeff(3) * ssi)/(acoeff(1) * ssi) * R_ik  ! KHL06 Eq. 19
            w_pre  = MAX(w_pre,0.0_wp)

            IF (atmo%w(i,k) > w_pre) THEN   ! homogenous nucleation event

              ! timescales of freezing event (see KL02, RM05, KHL06)
              cool    = grav / cp * atmo%w(i,k)
              ctau    = T_a * ( 0.004_wp*T_a - 2.0_wp ) + 304.4_wp
              tau     = 1.0_wp / (ctau * cool)                       ! freezing timescale, eq. (5)
              delta   = (bcoeff(2) * r_0)                         ! dimless aerosol radius, eq.(4)
              phi     = acoeff(1)*ssi / ( acoeff(2) + acoeff(3)*ssi) * (atmo%w(i,k) - w_pre)

              ! monodisperse approximation following KHL06
              kappa   = 2.0_wp * bcoeff(1) * bcoeff(2) * tau / (1.0_wp+ delta)**2  ! kappa, Eq. 8 KHL06
              sqrtkap = SQRT(kappa)                                        ! root of kappa
              ren     = 3.0_wp * sqrtkap / ( 2.0_wp + SQRT(1.0_wp+9.0_wp*kappa/pi) )       ! analy. approx. of erfc by RM05
              R_imfc  = 4.0_wp * pi * bcoeff(1)/bcoeff(2)**2 / svol
              R_im    = R_imfc / (1.0_wp+ delta) * ( delta**2 - 1.0_wp &
                   & + (1.0_wp+0.5_wp*kappa*(1.0_wp+ delta)**2) * ren/sqrtkap)           ! RIM Eq. 6 KHL06

              ! number concentration and radius of ice particles
              ni_hom  = phi / R_im                                         ! ni Eq.9 KHL06
              ri_0    = 1.0_wp + 0.5_wp * sqrtkap * ren                           ! for Eq. 3 KHL06
              ri_hom  = (ri_0 * (1.0_wp + delta) - 1.0_wp ) / bcoeff(2)            ! Eq. 3 KHL06 * REN = Eq.23 KHL06
              mi_hom  = (4.0_wp/3.0_wp * pi * rho_ice) * ni_hom * ri_hom**3
              mi_hom  = MAX(mi_hom,ice%x_min)

              nuc_n = MAX(MIN(ni_hom, ni_hom_max), 0.0_wp)
              nuc_q = MIN(nuc_n * mi_hom, atmo%qv(i,k))

              ice%n(i,k) = ice%n(i,k) + nuc_n
              ice%q(i,k) = ice%q(i,k) + nuc_q
              atmo%qv(i,k)  = atmo%qv(i,k)  - nuc_q

            END IF
          END IF
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    END IF

    !$ACC WAIT
    !$ACC END DATA ! nuc_n_a, ndiag_mask, acoeff, bcoeff, n_inpot, atmo, ice

  END SUBROUTINE ice_nucleation_homhet

  SUBROUTINE ice_nucleation_het_philips(ik_slice, atmo, ice, cloud, &
       use_prog_in, n_inact, ndiag_mask, nuc_n_a, n_inpot)
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(in) :: ice, cloud
    LOGICAL, INTENT(in) :: use_prog_in
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: n_inact
    REAL(wp), INTENT(out) :: &
         nuc_n_a(ik_slice(1):ik_slice(2), ik_slice(3):ik_slice(4))
    LOGICAL, INTENT(out) :: &
         ndiag_mask(ik_slice(1):ik_slice(2), ik_slice(3):ik_slice(4))
    REAL(wp), DIMENSION(:,:), OPTIONAL :: n_inpot

    REAL(wp)             :: nuc_n, nuc_q
    REAL(wp)             :: T_a, ssi, e_si
    REAL(wp)             :: ndiag, ndiag_dust, ndiag_all
    REAL(wp), PARAMETER  :: eps  = 1.0e-20_wp
    REAL(wp) :: infrac(3)
    LOGICAL :: lwrite_n_inpot

    ! variables for interpolation in look-up table (real is good enough here)
    REAL      :: xt,xs,ssr
    INTEGER   :: ss,tt

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER :: i, k

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC   PRIVATE(T_a, e_si, ssi, xt, tt, infrac, xs, ss, ssr) &
    !$ACC   PRIVATE(ndiag, ndiag_dust, ndiag_all, nuc_n, nuc_q, lwrite_n_inpot)
    DO k = kstart,kend
!NEC$ ivdep
      DO i = istart,iend

        T_a  = atmo%T(i,k)
        e_si = e_es(T_a)
        ssi  = atmo%qv(i,k) * R_d * T_a / e_si

        IF (T_a < T_nuc .AND. T_a > 180.0_wp .AND. ssi > 1.0_wp  &
             & .AND. ( n_inact(i,k) < ni_het_max*cfg_params%in_fact ) ) THEN

          xt = (274.- REAL(atmo%T(i,k)))  / ttstep
          xt = MIN(xt,REAL(ttmax-1))
          tt = INT(xt)

          IF (cloud%q(i,k) > eps) THEN
            ! immersion freezing at water saturation
            ! Phillips scheme
            ! immersion freezing at water saturation
            infrac(1) = (AINT(xt) + 1.0_wp - xt) * afrac_dust(tt,99) &
                 &        + (xt - AINT(xt)) * afrac_dust(tt+1,99)
            infrac(2) = (AINT(xt) + 1.0_wp - xt) * afrac_soot(tt,99) &
                 &        + (xt - AINT(xt)) * afrac_soot(tt+1,99)
            infrac(3) = (AINT(xt) + 1.0_wp - xt) * afrac_orga(tt,99) &
                 &        + (xt-AINT(xt)) * afrac_orga(tt+1,99)
          ELSE
            ! deposition nucleation below water saturation
            ! calculate indices used for 2D look-up tables
            xs = 100. * REAL(ssi-1.0_wp) / ssstep
            xs = MIN(xs,REAL(ssmax-1))
            ss = MAX(1,INT(xs))
            ssr = MAX(1.0, AINT(xs))
            ! bi-linear interpolation in look-up tables
            infrac(1) =   (AINT(xt) + 1.0_wp - xt) * (ssr + 1.0_wp - xs) &
                 &        * afrac_dust(tt, ss) &
                 &      + (xt - AINT(xt)) * (ssr + 1.0_wp - xs) &
                 &        * afrac_dust(tt+1, ss) &
                 &      + (AINT(xt) + 1.0_wp - xt) * (xs - ssr) &
                 &        * afrac_dust(tt, ss+1) &
                 &      + (xt - AINT(xt)) * (xs - ssr) &
                 &        * afrac_dust(tt+1, ss+1)
            infrac(2) =   (AINT(xt) + 1.0_wp - xt) * (ssr + 1.0_wp - xs) &
                 &        * afrac_soot(tt, ss) &
                 &      + (xt - AINT(xt)) * (ssr + 1.0_wp - xs) &
                 &        * afrac_soot(tt+1, ss  ) &
                 &      + (AINT(xt) + 1.0_wp - xt) * (xs - ssr) &
                 &        * afrac_soot(tt, ss + 1) &
                 &      + (xt - AINT(xt)) * (xs - ssr) &
                 &        * afrac_soot(tt+1, ss+1)
            infrac(3) = (AINT(xt) + 1.0_wp - xt) * (ssr + 1.0_wp - xs) &
                 &        * afrac_orga(tt,ss) &
                 &      + (xt - AINT(xt)) * (ssr + 1.0_wp - xs) &
                 &        * afrac_orga(tt+1, ss) &
                 &      + (AINT(xt) + 1.0_wp - xt) * (xs - ssr) &
                 &        * afrac_orga(tt, ss+1) &
                 &      + (xt - AINT(xt)) * (xs - ssr) &
                 &        * afrac_orga(tt+1, ss+1)
          END IF
          
          ! sum up the three modes
          IF (use_prog_in) THEN
            ! n_inpot replaces na_dust, na_soot and na_orga are assumed to be constant
            ndiag  = n_inpot(i,k) * infrac(1) + na_soot * infrac(2) + na_orga * infrac(3)
            ndiag_dust = n_inpot(i,k) * infrac(1)
            ndiag_all = ndiag
          ELSE
            ! all aerosol species are diagnostic
            ndiag = na_dust * infrac(1) + na_soot * infrac(2) + na_orga * infrac(3)
            ndiag_dust = ndiag
            ndiag_all = ndiag
          END IF
          ndiag = MIN(ndiag,ni_het_max)

          nuc_n = MAX(ndiag*cfg_params%in_fact - n_inact(i,k),0.0_wp)
          nuc_q = MIN(nuc_n * ice%x_min, atmo%qv(i,k))
          nuc_n = nuc_q / ice%x_min

          ice%n(i,k)   = ice%n(i,k)   + nuc_n
          ice%q(i,k)   = ice%q(i,k)   + nuc_q
          atmo%qv(i,k) = atmo%qv(i,k) - nuc_q
          n_inact(i,k) = n_inact(i,k) + nuc_n

          lwrite_n_inpot = use_prog_in .AND. ndiag .GT. 1.0e-12_wp
          ndiag_mask(i, k) = lwrite_n_inpot

          IF (lwrite_n_inpot) THEN
            nuc_n = nuc_n * ndiag_dust / ndiag_all
          END IF
          nuc_n_a(i, k) = nuc_n

        ELSE
          nuc_n_a(i, k) = 0.0_wp
          ndiag_mask(i, k) = .FALSE.
        ENDIF

      END DO
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE ice_nucleation_het_philips

  SUBROUTINE ice_nucleation_het_hdcp2(ik_slice, atmo, ice, cloud, &
       use_prog_in, nuc_coeffs, n_inact, ndiag_mask, nuc_n_a)
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(in) :: ice, cloud
    LOGICAL, INTENT(in) :: use_prog_in
    TYPE(dep_imm_coeffs), INTENT(in) :: nuc_coeffs
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: n_inact
    REAL(wp), INTENT(out) :: &
         nuc_n_a(ik_slice(1):ik_slice(2), ik_slice(3):ik_slice(4))
    LOGICAL, INTENT(out) :: &
         ndiag_mask(ik_slice(1):ik_slice(2), ik_slice(3):ik_slice(4))


    ! parameters for deposition formula, Eq (2) of Hande et al.
    REAL(wp), PARAMETER :: a_dep =  2.7626_wp * 0.1_wp
    REAL(wp), PARAMETER :: b_dep =  6.2100_wp  ! 0.0621*100, because we use ssi instead of RHi
    REAL(wp), PARAMETER :: c_dep = -1.3107_wp
    REAL(wp), PARAMETER :: d_dep =  2.6789_wp * 0.1_wp

    REAL(wp)             :: nuc_n, nuc_q
    REAL(wp)             :: T_a, ssi, e_si
    REAL(wp)             :: ndiag
    REAL(wp), PARAMETER  :: eps  = 1.0e-20_wp

    LOGICAL :: lwrite_n_inpot

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER :: i, k

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
      DO i = istart,iend

        T_a  = atmo%T(i,k)
        e_si = e_es(T_a)
        ssi  = atmo%qv(i,k) * R_d * T_a / e_si

        IF (T_a < T_nuc .AND. T_a > 180.0_wp .AND. ssi > 1.0_wp  &
             & .AND. ( n_inact(i,k) < ni_het_max*cfg_params%in_fact ) ) THEN

          IF (cloud%q(i,k) > eps) THEN
            ! Hande et al. scheme, Eq. (1)
            T_a = MAX(T_a,237.1501_wp)
            IF (T_a.LT.261.15_wp) THEN
              ndiag = nuc_coeffs%nin_imm * EXP( - nuc_coeffs%alf_imm &
                   * EXP(nuc_coeffs%bet_imm*LOG(T_a - 237.15_wp)) )
            ELSE
              ndiag = 0.0_wp
            END IF
          ELSE
            ! Hande et al. scheme, Eq. (3) with (2) and (1)
            T_a = max(T_a,220.001_wp)
            IF (T_a.LT.253.0_wp) THEN
              ndiag = nuc_coeffs%nin_dep * EXP( - nuc_coeffs%alf_dep &
                   * EXP(nuc_coeffs%bet_dep*LOG(T_a - 220.0_wp)) )
              ndiag = ndiag * (a_dep * atan(b_dep*(ssi-1.0_wp)+c_dep) + d_dep)
            ELSE
              ndiag = 0.0_wp
            END IF
          END IF


          nuc_n = MAX(ndiag*cfg_params%in_fact - n_inact(i,k),0.0_wp)
          nuc_q = MIN(nuc_n * ice%x_min, atmo%qv(i,k))
          nuc_n = nuc_q / ice%x_min

          ice%n(i,k)   = ice%n(i,k)   + nuc_n
          ice%q(i,k)   = ice%q(i,k)   + nuc_q
          atmo%qv(i,k) = atmo%qv(i,k) - nuc_q
          n_inact(i,k) = n_inact(i,k) + nuc_n

          lwrite_n_inpot = use_prog_in .AND. ndiag .GT. 1.0e-12_wp
          ndiag_mask(i, k) = lwrite_n_inpot

          nuc_n_a(i, k) = nuc_n

        ELSE
          nuc_n_a(i, k) = 0.0_wp
          ndiag_mask(i, k) = .FALSE.
        ENDIF

      END DO
    END DO
  END SUBROUTINE ice_nucleation_het_hdcp2

  SUBROUTINE vapor_dep_relaxation(ik_slice, dt_local, &
       &               ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs, &
       &               atmo, ice_in, snow_in, graupel_in, hail_in, dep_rate_ice, dep_rate_snow)
    !*******************************************************************************
    ! Deposition and sublimation                                                   *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    TYPE(atmosphere)    :: atmo
    CLASS(particle), INTENT(INOUT), TARGET :: ice_in, snow_in, graupel_in, hail_in
    CLASS(particle), POINTER :: ice, snow, graupel, hail  ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle_sphere), INTENT(IN) :: ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs
    REAL(wp), INTENT(IN) :: dt_local
    REAL(wp), INTENT(INOUT), DIMENSION(:,:) :: dep_rate_ice, dep_rate_snow

    REAL(wp), DIMENSION(size(dep_rate_ice,1),size(dep_rate_ice,2)) :: &
                                    & s_si,g_i,dep_ice,dep_snow,dep_graupel,dep_hail

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: D_vtp
    REAL(wp)            :: zdt,qvsidiff,Xi_i,Xfac
    REAL(wp)            :: tau_i_i,tau_s_i,tau_g_i,tau_h_i
    REAL(wp), PARAMETER :: eps  = 1.0e-20_wp
    REAL(wp)            :: T_a
    REAL(wp)            :: e_si            !..saturation water pressure over ice
    REAL(wp)            :: e_d,p_a,dep_sum !,weight
    REAL(wp)            :: dep_ice_n,dep_snow_n,dep_graupel_n,dep_hail_n,x_i,x_s,x_g,x_h

    LOGICAL, PARAMETER  :: reduce_sublimation = .TRUE.
    REAL(wp), PARAMETER :: dep_n_fac = 0.5_wp  ! UB: if this new parameterization of n-reduction during sublimation
                                               !     really makes sense, move to a global constant or into the particle types
    
    IF (isdebug) CALL message(routine, "vapor_deposition_growth")

    ice => ice_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    snow => snow_in
    graupel => graupel_in
    hail => hail_in

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC DATA CREATE(s_si, g_i, dep_ice, dep_snow, dep_graupel, dep_hail)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(p_a, T_a, e_d, e_si, D_vtp)
    DO k = kstart,kend
       DO i = istart,iend
          p_a  = atmo%p(i,k)
          T_a  = atmo%T(i,k)
          IF (T_a < T_3) THEN
             e_d  = atmo%qv(i,k) * R_d * T_a
             e_si = e_es(T_a)
             s_si(i,k) = e_d / e_si - 1.0_wp    !..supersaturation over ice
             D_vtp = diffusivity(T_a,p_a)    !  D_v = 8.7602e-5 * T_a**(1.81) / p_a
             g_i(i,k) = 4.0_wp*pi / ( L_ed**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_si) )
          ELSE
             g_i(i,k)  = 0.0_wp
             s_si(i,k) = 0.0_wp
          ENDIF
       ENDDO
    ENDDO
    !$ACC END PARALLEL


    CALL vapor_deposition_generic(ik_slice, ice, ice_coeffs, g_i, s_si,dt_local, dep_ice)
    CALL vapor_deposition_generic(ik_slice, snow, snow_coeffs, g_i, s_si, dt_local, dep_snow)
    CALL vapor_deposition_generic(ik_slice, graupel, graupel_coeffs, g_i, s_si, dt_local, dep_graupel)
    CALL vapor_deposition_generic(ik_slice, hail, hail_coeffs, g_i, s_si, dt_local, dep_hail)

    zdt = 1.0/dt_local

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(qvsidiff, tau_i_i, tau_s_i, tau_g_i, tau_h_i, Xi_i, Xfac) &
    !$ACC   PRIVATE(dep_ice_n, dep_snow_n, dep_graupel_n, dep_hail_n, dep_sum, x_i, x_s, x_g, x_h, T_a)
    DO k = kstart,kend
       DO i = istart,iend

          T_a  = atmo%T(i,k)

          ! Deposition only below T_3, evaporation of melting particles at warmer T is treated elsewhere
          IF (T_a < T_3) THEN

             ! Depositional growth with relaxation time-scale approach based on:
             ! "A New Double-Moment Microphysics Parameterization for Application in Cloud and
             ! Climate Models. Part 1: Description" by H. Morrison, J.A.Curry, V.I. Khvorostyanov

             qvsidiff  = atmo%qv(i,k) - e_es(T_a)/(R_d*T_a)

             if (abs(qvsidiff).gt.eps) then

                ! deposition rates are already multiplied with dt_local, therefore divide them here
                tau_i_i  = zdt/qvsidiff*dep_ice(i,k)
                tau_s_i  = zdt/qvsidiff*dep_snow(i,k)
                tau_g_i  = zdt/qvsidiff*dep_graupel(i,k)
                tau_h_i  = zdt/qvsidiff*dep_hail(i,k)

                Xi_i = ( tau_i_i + tau_s_i + tau_g_i + tau_h_i )

                if (Xi_i.lt.eps) then
                   Xfac = 0.0_wp
                else
                   Xfac =  qvsidiff / Xi_i * (1.0_wp - EXP(- dt_local*Xi_i))
                end if

                dep_ice(i,k)     = Xfac * tau_i_i
                dep_snow(i,k)    = Xfac * tau_s_i
                dep_graupel(i,k) = Xfac * tau_g_i
                dep_hail(i,k)    = Xfac * tau_h_i

                ! this limiter should not be necessary
                IF (qvsidiff < 0.0_wp) THEN
                   dep_ice(i,k)     = MAX(dep_ice(i,k),    -ice%q(i,k))
                   dep_snow(i,k)    = MAX(dep_snow(i,k),   -snow%q(i,k))
                   dep_graupel(i,k) = MAX(dep_graupel(i,k),-graupel%q(i,k))
                   dep_hail(i,k)    = MAX(dep_hail(i,k),   -hail%q(i,k))
                END IF

                dep_sum = dep_ice(i,k) + dep_graupel(i,k) + dep_snow(i,k) + dep_hail(i,k)
                
                IF ( reduce_sublimation) THEN
                  x_i = particle_meanmass(ice    , ice%    q(i,k), ice%    n(i,k))
                  x_s = particle_meanmass(snow   , snow%   q(i,k), snow%   n(i,k))
                  x_g = particle_meanmass(graupel, graupel%q(i,k), graupel%n(i,k))
                  x_h = particle_meanmass(hail   , hail%   q(i,k), hail%   n(i,k))
                END IF

                ice%q(i,k)     = ice%q(i,k)     + dep_ice(i,k)
                snow%q(i,k)    = snow%q(i,k)    + dep_snow(i,k)
                graupel%q(i,k) = graupel%q(i,k) + dep_graupel(i,k)
                hail%q(i,k)    = hail%q(i,k)    + dep_hail(i,k)

                atmo%qv(i,k) = atmo%qv(i,k) - dep_sum

                ! .. If deposition rate is negative, parameterize the complete evaporation of some of the particles in a way
                !    that mean size is conserved times a tuning factor < 1:
                IF ( reduce_sublimation) THEN
                  dep_ice_n      = MIN(dep_ice(i,k),0.0_wp) / x_i
                  dep_snow_n     = MIN(dep_snow(i,k),0.0_wp) / x_s
                  dep_graupel_n  = MIN(dep_graupel(i,k),0.0_wp) / x_g
                  dep_hail_n     = MIN(dep_hail(i,k),0.0_wp) / x_h

                  ice%n(i,k)     = MAX(ice%n(i,k)     + dep_n_fac*dep_ice_n    , 0.0_wp)
                  snow%n(i,k)    = MAX(snow%n(i,k)    + dep_n_fac*dep_snow_n   , 0.0_wp)
                  graupel%n(i,k) = MAX(graupel%n(i,k) + dep_n_fac*dep_graupel_n, 0.0_wp)
                  hail%n(i,k)    = MAX(hail%n(i,k)    + dep_n_fac*dep_hail_n   , 0.0_wp)
                END IF
                
                dep_rate_ice(i,k)  = dep_rate_ice(i,k)  + dep_ice(i,k)
                dep_rate_snow(i,k) = dep_rate_snow(i,k) + dep_snow(i,k)

             END IF

          ENDIF
       ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE vapor_dep_relaxation

  SUBROUTINE vapor_deposition_generic(ik_slice, prtcl_in, coeffs, g_i, s_si, &
       dt, dep_q)
    INTEGER, INTENT(in) :: ik_slice(4)
    CLASS(particle), INTENT(in), TARGET :: prtcl_in
    CLASS(particle), POINTER :: prtcl ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle_coeffs), INTENT(in) :: coeffs
    REAL(wp), INTENT(in) :: g_i(:, :), s_si(:, :)
    REAL(wp), INTENT(in) :: dt
    REAL(wp), INTENT(out) :: dep_q(:, :)
    REAL(wp)            :: q,n,x,d,v,f_v
    INTEGER             :: i,k
    INTEGER :: istart, iend, kstart, kend

    prtcl => prtcl_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(n, q, x, D, v, f_v)
    DO k = kstart,kend
      DO i = istart,iend
        IF (prtcl%q(i,k) == 0.0_wp) THEN
          dep_q(i,k) = 0.0_wp
        ELSE
          n = prtcl%n(i,k)
          q = prtcl%q(i,k)

          x = particle_meanmass(prtcl,q,n)
          D = particle_diameter(prtcl,x)
          v = particle_velocity(prtcl,x) * prtcl%rho_v(i,k)

          f_v  = coeffs%a_f + coeffs%b_f * SQRT(D*v)
          f_v  = MAX(f_v,coeffs%a_f/prtcl%a_ven)

          dep_q(i,k) = g_i(i,k) * n * coeffs%c_i * d * f_v * s_si(i,k) * dt
        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE vapor_deposition_generic

  SUBROUTINE rain_freeze_gamlook(ik_slice, dt, rain_ltable1, rain_ltable2, rain_ltable3, &
                                 rain_nm1, rain_nm2, rain_nm3, rain_g1, rain_g2,         &
                                 rain_coeffs, atmo, rain_in, ice, snow, graupel, hail)
    !*******************************************************************************
    ! Freezing of raindrops                                                        *
    ! by Uli Blahak                                                                *
    !                                                                              *
    ! incomplete gamma functions are implemented as look-up tables                 *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    ! coefficients and tables
    TYPE(gamlookuptable), INTENT(in) :: rain_ltable1, rain_ltable2, rain_ltable3
    REAL(wp),INTENT(in)              :: rain_nm1, rain_nm2, rain_nm3, rain_g1, rain_g2
    TYPE(particle_rain_coeffs),INTENT(in) :: rain_coeffs
    ! prognostic variables
    TYPE(atmosphere), INTENT(inout)  :: atmo
    CLASS(particle), INTENT(inout), TARGET :: rain_in
    CLASS(particle), POINTER :: rain ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle), INTENT(inout) :: ice, snow, graupel, hail

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: fr_q,fr_n,T_a,q_r,x_r,n_r,j_het,               &
         &                 fr_q_i,fr_n_i,fr_q_g,fr_n_g,fr_q_h,fr_n_h,n_0, &
         &                 lam,xmax_ice,xmax_gr,fr_q_tmp,fr_n_tmp,        &
         &                 lam_rnm1, lam_rnm2, lam_rnm3

    REAL(wp), PARAMETER :: a_HET = 6.5d-1      ! Data of Barklie and Gokhale (PK S.350)
    REAL(wp), PARAMETER :: b_HET = 2.0d+2      !         Barklie and Gokhale (PK S.350)
    REAL(wp), PARAMETER :: eps = 1e-15_wp      ! for clipping
    LOGICAL,  PARAMETER :: lclipping = .true.

    IF (isdebug) CALL message(routine, "rain_freeze_gamlook")

    rain => rain_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    xmax_ice = ( (cfg_params%D_rainfrz_ig/rain%a_geo)**(1.0_wp/rain%b_geo) )**rain%mu
    xmax_gr  = ( (cfg_params%D_rainfrz_gh/rain%a_geo)**(1.0_wp/rain%b_geo) )**rain%mu

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    ! ACCWA (nvhpc 21.3): Explicit PRESENT statement required for rain_ltable*
    ! - Otherwise error: illegal address during kernel execution
    ! - Reason unknown as inlining of incgfct_lower_lookup solves issue
    !$ACC DATA PRESENT(rain_ltable1, rain_ltable2, rain_ltable3)
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, q_r, n_r, fr_q, fr_n, fr_n_i, fr_q_i, fr_n_g) &
    !$ACC   PRIVATE(fr_q_g, fr_n_h, fr_q_h, fr_n_tmp, fr_q_tmp) &
    !$ACC   PRIVATE(x_r, lam, lam_rnm1, lam_rnm2, lam_rnm3, n_0, j_het)
    DO k = kstart,kend
!NEC$ ivdep
       DO i = istart,iend

          T_a = atmo%T(i,k)
          q_r = rain%q(i,k)
          n_r = rain%n(i,k)

          IF (T_a < T_freeze) THEN
             IF (q_r <= q_crit_fr) THEN
                IF (T_a < T_f) THEN  ! instantaneous freezing below T_f or -40 C
                   fr_q = q_r
                   fr_n = n_r
                   fr_n_i= n_r
                   fr_q_i= q_r
                   fr_n_g= 0.0_wp
                   fr_q_g= 0.0_wp
                   fr_n_h= 0.0_wp
                   fr_q_h= 0.0_wp
                   fr_n_tmp = 1.0_wp
                   fr_q_tmp = 1.0_wp
                ELSE
                   fr_q = 0.0   ; fr_n = 0.0_wp
                   fr_n_i = 0.0 ; fr_q_i = 0.0_wp
                   fr_n_g = 0.0 ; fr_q_g = 0.0_wp
                   fr_n_h = 0.0 ; fr_q_h = 0.0_wp
                   fr_n_tmp = 0.0 ; fr_q_tmp = 0.0_wp
                END IF
             ELSE
                x_r = particle_meanmass(rain, q_r,n_r)
                n_r = q_r / x_r
                IF (T_a < T_f) THEN
                   ! Diesen Zweig koennte man auch weglassen. ist zudem zwar quantitativ richtig,
                   ! aber nicht konsistent zum
                   ! Grenzfall fuer komplettes Gefrieren der Rechnung im T_a >= T_f - Zweig weiter unten
                   fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
                   fr_n = n_r
                   ! Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                   ! oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                   ! bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                   ! und von xmax_gr bis unendlich (--> Hagel).
                   lam = EXP( LOG( rain_g1/rain_g2*x_r ) * (-rain%mu) )
                   lam_rnm1 = EXP(rain_nm1*LOG(lam))  ! lam**rain_nm1
                   lam_rnm2 = EXP(rain_nm2*LOG(lam))  ! lam**rain_nm2
                   n_0 = rain%mu * n_r * lam_rnm1 / rain_g1
                   fr_n_i = n_0/(rain%mu*lam_rnm1) * incgfct_lower_lookup(lam*xmax_ice,rain_ltable1)
                   fr_q_i = n_0/(rain%mu*lam_rnm2) * incgfct_lower_lookup(lam*xmax_ice,rain_ltable2)
                   fr_n_g = n_0/(rain%mu*lam_rnm1) * incgfct_lower_lookup(lam*xmax_gr, rain_ltable1)
                   fr_q_g = n_0/(rain%mu*lam_rnm2) * incgfct_lower_lookup(lam*xmax_gr, rain_ltable2)

                   fr_n_h = fr_n - fr_n_g
                   fr_q_h = fr_q - fr_q_g
                   fr_n_g = fr_n_g - fr_n_i
                   fr_q_g = fr_q_g - fr_q_i
                   fr_n_tmp = n_r/MAX(fr_n,n_r)
                   fr_q_tmp = q_r/MAX(fr_q,q_r)
                ELSE
                   !..heterogeneous freezing
                   j_het = MAX(b_HET * ( EXP( a_HET * (T_3 - T_a)) - 1.0_wp ),0.0_wp) / rho_w * dt
!!                   if (use_prog_in) j_het = MIN(j_het, n_inact(i,k)/q_r)

                   ! Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                   ! oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                   ! bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                   ! und von xmax_gr bis unendlich (--> Hagel).
                   IF (j_het >= 1.0e-20_wp) THEN
                      fr_n  = j_het * q_r
                      fr_q  = j_het * q_r * x_r * rain_coeffs%c_z

!                      lam = ( rain_g1 / rain_g2 * x_r)**(-rain%mu)
                      lam = EXP( LOG( rain_g1/rain_g2*x_r ) * (-rain%mu) )
                      lam_rnm1 = EXP(rain_nm1*LOG(lam))  ! lam**rain_nm1
                      lam_rnm2 = EXP(rain_nm2*LOG(lam))  ! lam**rain_nm2
                      lam_rnm3 = EXP(rain_nm3*LOG(lam))  ! lam**rain_nm3
                      n_0 = rain%mu * n_r * lam_rnm1 / rain_g1
                      fr_n_i = j_het * n_0/(rain%mu*lam_rnm2) * incgfct_lower_lookup(lam*xmax_ice,rain_ltable2)
                      fr_q_i = j_het * n_0/(rain%mu*lam_rnm3) * incgfct_lower_lookup(lam*xmax_ice,rain_ltable3)
                      fr_n_g = j_het * n_0/(rain%mu*lam_rnm2) * incgfct_lower_lookup(lam*xmax_gr, rain_ltable2)
                      fr_q_g = j_het * n_0/(rain%mu*lam_rnm3) * incgfct_lower_lookup(lam*xmax_gr, rain_ltable3)

                      fr_n_h = fr_n - fr_n_g
                      fr_q_h = fr_q - fr_q_g
                      fr_n_g = fr_n_g - fr_n_i
                      fr_q_g = fr_q_g - fr_q_i
                      fr_n_tmp = n_r/MAX(fr_n,n_r)
                      fr_q_tmp = q_r/MAX(fr_q,q_r)
                   ELSE
                      fr_n= 0.0
                      fr_q= 0.0
                      fr_n_i= 0.0
                      fr_q_i= 0.0
                      fr_n_g= 0.0
                      fr_q_g= 0.0
                      fr_n_h= 0.0
                      fr_q_h= 0.0
                      fr_n_tmp = 0.0
                      fr_q_tmp = 0.0
                   END IF

                END IF

                fr_n = fr_n * fr_n_tmp
                fr_q = fr_q * fr_q_tmp
                fr_n_i = fr_n_i * fr_n_tmp
                fr_n_g = fr_n_g * fr_n_tmp
                fr_n_h = fr_n_h * fr_n_tmp
                fr_q_i = fr_q_i * fr_q_tmp
                fr_q_g = fr_q_g * fr_q_tmp
                fr_q_h = fr_q_h * fr_q_tmp
             END IF

             rain%q(i,k) = rain%q(i,k) - fr_q
             rain%n(i,k) = n_r - fr_n
             !if (use_prog_in) then
             !   n_inact(i,k) = n_inact(i,k) + fr_n
             !end if

             ! mit Hagelklasse, gefrierender Regen wird Eis, Graupel oder Hagel
!!$ UB            snow%q(i,k) = snow%q(i,k)  + fr_q_i
!!$ UB            snow%n(i,k) = snow%n(i,k)  + fr_n_i   ! put this into snow
             ice%q(i,k) = ice%q(i,k)  + fr_q_i    ! ... or into ice? --> UB: original idea was to put it into ice
             ice%n(i,k) = ice%n(i,k)  + fr_n_i

             graupel%q(i,k) = graupel%q(i,k)  + fr_q_g
             graupel%n(i,k) = graupel%n(i,k)  + fr_n_g
             hail%q(i,k) = hail%q(i,k)  + fr_q_h
             hail%n(i,k) = hail%n(i,k)  + fr_n_h

             ! clipping of small negatives is necessary here
             if (lclipping) then
                IF (rain%q(i,k) < 0.0_wp .and. abs(rain%q(i,k)) < eps) rain%q(i,k) = 0.0_wp
                IF (rain%n(i,k) < 0.0_wp .and. abs(rain%n(i,k)) < eps) rain%n(i,k) = 0.0_wp
                IF (graupel%q(i,k) < 0.0_wp .and. abs(graupel%q(i,k)) < eps) graupel%q(i,k) = 0.0_wp
                IF (graupel%n(i,k) < 0.0_wp .and. abs(graupel%q(i,k)) < eps) graupel%n(i,k) = 0.0_wp
                IF (hail%q(i,k) < 0.0_wp .and. abs(hail%q(i,k)) < eps) hail%q(i,k) = 0.0_wp
                IF (hail%n(i,k) < 0.0_wp .and. abs(hail%n(i,k)) < eps) hail%n(i,k) = 0.0_wp
             end if

          END IF
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer
    !$ACC END DATA

  END SUBROUTINE rain_freeze_gamlook

  SUBROUTINE ice_selfcollection(ik_slice, dt, atmo, ice_in, snow, ice_coeffs, ltab_estick_ice)
    !*******************************************************************************
    ! selfcollection of ice crystals, see SB2006 or Seifert (2002)                 *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt    
    TYPE(particle_ice_coeffs), INTENT(in) :: ice_coeffs

    ! 2mom variables
    TYPE(atmosphere),       INTENT(inout) :: atmo
    CLASS(particle_frozen), INTENT(inout), TARGET :: ice_in
    CLASS(particle_frozen), POINTER :: ice ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle_frozen), INTENT(inout) :: snow
    TYPE(lookupt_1D), INTENT(in), TARGET :: ltab_estick_ice
    TYPE(lookupt_1D), POINTER :: p_ltab_estick_ice

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend

    ! local variables
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_i,n_i,x_i,d_i,v_i,e_coll,x_conv_ii
    REAL(wp)            :: self_n,self_q

    IF (isdebug) CALL message(routine, "ice_selfcollection")

    ice => ice_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    p_ltab_estick_ice => ltab_estick_ice
    
    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    x_conv_ii = (cfg_params%D_conv_ii/snow%a_geo)**(1./snow%b_geo)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, q_i, n_i, x_i, d_i, v_i, e_coll, self_n, self_q)
    DO k = kstart,kend
       DO i = istart,iend

          q_i = ice%q(i,k)
          n_i = ice%n(i,k)

          x_i = particle_meanmass(ice,q_i,n_i)
          D_i = particle_diameter(ice,x_i)

          IF ( n_i > 0.0_wp .AND. q_i > q_crit_ii .AND. D_i > D_crit_ii ) THEN

             T_a = atmo%T(i,k)

             !.. Sticking efficiency depending on temperature:
             e_coll = estick_ltab_equi(T_a, p_ltab_estick_ice)  ! equidistant lookup table
             
             v_i = ice%a_vel * x_i**ice%b_vel * ice%rho_v(i,k)

             self_n = pi4 * e_coll * ice_coeffs%sc_delta_n * n_i * n_i * D_i * D_i &
                  & * SQRT( ice_coeffs%sc_theta_n * v_i * v_i + 2.0_wp*ice%s_vel**2 ) * dt

             self_q = pi4 * e_coll * ice_coeffs%sc_delta_q * n_i * q_i * D_i * D_i &
                  & * SQRT( ice_coeffs%sc_theta_q * v_i * v_i + 2.0_wp*ice%s_vel**2 ) * dt

             self_q = MIN(self_q,q_i)
             self_n = MIN(MIN(self_n,self_q/x_conv_ii),n_i)

             ice%q(i,k)  = ice%q(i,k)  - self_q
             snow%q(i,k) = snow%q(i,k) + self_q

             ice%n(i,k)  = ice%n(i,k)  - self_n
             snow%n(i,k) = snow%n(i,k) + self_n / 2.0_wp

          ENDIF
       ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE ice_selfcollection

  SUBROUTINE snow_selfcollection(ik_slice, dt, atmo, snow_in, snow_coeffs, ltab_estick_snow)
    !*******************************************************************************
    ! Selfcollection of snow                                                       *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    
    TYPE(particle_snow_coeffs), INTENT(in) :: snow_coeffs
    TYPE(atmosphere), INTENT(inout)        :: atmo
    CLASS(particle_frozen), INTENT(inout), TARGET :: snow_in
    CLASS(particle_frozen), POINTER :: snow ! ACCWA (nvhpc 22.7, IPSF, see above)
    TYPE(lookupt_1D), INTENT(in), TARGET :: ltab_estick_snow
    TYPE(lookupt_1D), POINTER :: p_ltab_estick_snow

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,e_coll
    REAL(wp)            :: self_n

    IF (isdebug) CALL message(routine, "snow_selfcollection")

    snow => snow_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    p_ltab_estick_snow => ltab_estick_snow

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, q_s, n_s, x_s, d_s, v_s, e_coll, self_n)
    DO k = kstart,kend
      DO i = istart,iend

          q_s = snow%q(i,k)
          T_a = atmo%T(i,k)

          IF ( q_s > q_crit ) THEN

             !.. Sticking efficiency depending on temperature:
             e_coll = estick_ltab_equi(T_a, p_ltab_estick_snow)  ! equidistant lookup table

             n_s = snow%n(i,k)
             x_s = particle_meanmass(snow,q_s,n_s)
             D_s = particle_diameter(snow,x_s)
             v_s = particle_velocity(snow,x_s) * snow%rho_v(i,k)

             self_n = pi8 * e_coll * n_s * n_s * snow_coeffs%sc_delta_n * D_s * D_s * &
                  &          SQRT(  snow_coeffs%sc_theta_n * v_s * v_s                &
                  &               + 2.0_wp * snow%s_vel**2 ) * dt

             self_n = MIN(self_n,n_s)

             snow%n(i,k) = snow%n(i,k) - self_n

          ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE snow_selfcollection

  SUBROUTINE snow_melting(ik_slice, dt, snow_coeffs, atmo, snow_in, rain)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    CLASS(particle_sphere), INTENT(in) :: snow_coeffs
    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(inout), TARGET :: snow_in
    CLASS(particle), POINTER :: snow ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle), INTENT(inout) :: rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,T_a,e_a
    REAL(wp)            :: melt,melt_v,melt_h,melt_n,melt_q
    REAL(wp)            :: fh_q, fv_q

    IF (isdebug) CALL message(routine, "snow_melting")

    snow => snow_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, q_s, e_a, n_s, x_s, D_s, v_s) &
    !$ACC   PRIVATE(fv_q, fh_q, melt, melt_h, melt_v, melt_q, melt_n)
    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)
          q_s = snow%q(i,k)
          IF (T_a > T_3 .AND. q_s > 0.0_wp) THEN
!            e_a = e_ws(T_a)            ! saturation pressure IS WRONG HERE, need actual e_v
            e_a  = atmo%qv(i,k) * R_d * T_a
            n_s = snow%n(i,k)

            x_s = particle_meanmass(snow,q_s,n_s)
            D_s = particle_diameter(snow,x_s)
            v_s = particle_velocity(snow,x_s) * snow%rho_v(i,k)

            fv_q = snow_coeffs%a_f + snow_coeffs%b_f * SQRT(v_s*D_s)

            ! UB: Based on Rasmussen and Heymsfield (1987) the ratio fh_q / fv_q is approx. 1.05
            !     for a wide temperature- and pressure range:
            fh_q = 1.05 * fv_q

            melt   = 2.0*pi/L_ew * D_s * n_s * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)
            melt_q = cfg_params%melt_g_tune_fak * (melt_h * fh_q + melt_v * fv_q)

            ! UB: for melt_n we assume that x_s is constant during melting
            melt_n = MIN(MAX( (melt_q - q_s) / x_s + n_s, 0.0_wp), n_s)

            melt_q = MIN(q_s,MAX(melt_q,0.0_wp))
            melt_n = MIN(n_s,MAX(melt_n,0.0_wp))

            ! UB: snow melts instantaneously at 10 C
            IF (T_a - T_3 > 10.0_wp) THEN
              melt_q = q_s
              melt_n = n_s
            ENDIF

            snow%q(i,k) = snow%q(i,k) - melt_q
            rain%q(i,k) = rain%q(i,k) + melt_q

            snow%n(i,k) = snow%n(i,k) - melt_n
            rain%n(i,k) = rain%n(i,k) + melt_n

            snow%n(i,k) = MAX(snow%n(i,k), snow%q(i,k)/snow%x_max)

         ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE snow_melting

  SUBROUTINE particle_particle_collection(ik_slice, dt, atmo, ctype_in, ptype_in, coeffs, ltab_estick_parti)
    !*******************************************************************************
    !  Most simple particle-particle collection for ice particles, e.g.,           *
    !    graupel+ice  -> graupel                                                   *
    !    graupel+snow -> graupel                                                   *
    !    hail+ice  -> hail                                                         *
    !    hail+snow -> hail                                                         *
    !    snow+ice  -> snow                                                         *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt

    ! 2mom variables and coefficients
    TYPE(atmosphere), INTENT(inout)        :: atmo
    CLASS(particle_frozen), INTENT(inout), TARGET :: ctype_in, ptype_in
    CLASS(particle_frozen), POINTER :: ctype, ptype ! ACCWA (nvhpc 22.7, IPSF, see above)
    TYPE(collection_coeffs), INTENT(in)    :: coeffs
    TYPE(lookupt_1D), INTENT(in), TARGET :: ltab_estick_parti
    TYPE(lookupt_1D), POINTER :: p_ltab_estick_parti

    ! Locale Variablen
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_p,n_p,x_p,d_p,v_p
    REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
    REAL(wp)            :: coll_n,coll_q,e_coll

    IF (isdebug) CALL message(routine, "particle_particle_collection")

    ctype => ctype_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    ptype => ptype_in
    p_ltab_estick_parti => ltab_estick_parti

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC   PRIVATE(T_a, q_p, n_p, x_p, d_p, v_p, q_i, n_i, x_i, d_i, v_i, coll_n, coll_q, e_coll)
    DO k = kstart,kend
      DO i = istart,iend
        
        q_i = ctype%q(i,k)
        q_p = ptype%q(i,k)
        
        IF (q_i > q_crit .AND. q_p > q_crit) THEN
          T_a = atmo%T(i,k)
          
          !.. Sticking efficiency depending on temperature:
          e_coll = estick_ltab_equi(T_a, p_ltab_estick_parti)  ! equidistant lookup table
          
          n_i = ctype%n(i,k)
          n_p = ptype%n(i,k)
          
          x_p = particle_meanmass(ptype, q_p,n_p)
          d_p = particle_diameter(ptype, x_p)
          v_p = particle_velocity(ptype, x_p) * ptype%rho_v(i,k)
          
          x_i = particle_meanmass(ctype, q_i,n_i)
          d_i = particle_diameter(ctype, x_i)
          v_i = particle_velocity(ctype, x_i) * ctype%rho_v(i,k)
          
          coll_n = pi4 * n_p * n_i * e_coll * dt      &
               & *     ( coeffs%delta_n_aa * D_p**2   &
               &       + coeffs%delta_n_ab * D_p*D_i  &
               &       + coeffs%delta_n_bb * D_i**2)  &
               & * sqrt( coeffs%theta_n_aa * v_p**2   &
               &       - coeffs%theta_n_ab * v_p*v_i  &
               &       + coeffs%theta_n_bb * v_i**2   &
               &       + ctype%s_vel**2)
          
          coll_q = pi4 * n_p * q_i * e_coll * dt      &
               & *     ( coeffs%delta_q_aa * D_p**2   &
               &       + coeffs%delta_q_ab * D_p*D_i  &
               &       + coeffs%delta_q_bb * D_i**2)  &
               & * sqrt( coeffs%theta_q_aa * v_p**2   &
               &       - coeffs%theta_q_ab * v_p*v_i  &
               &       + coeffs%theta_q_bb * v_i**2   &
               &       + ctype%s_vel**2)

          coll_n = MIN(n_i,coll_n)
          coll_q = MIN(q_i,coll_q)

          ptype%q(i,k) = ptype%q(i,k) + coll_q
          ctype%q(i,k)  = ctype%q(i,k)  - coll_q
          ctype%n(i,k)  = ctype%n(i,k)  - coll_n

        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE particle_particle_collection

  SUBROUTINE graupel_selfcollection(ik_slice, dt, atmo, graupel_in, graupel_coeffs)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    TYPE(particle_graupel_coeffs), INTENT(in) :: graupel_coeffs
    TYPE(atmosphere), INTENT(inout)           :: atmo
    CLASS(particle), INTENT(inout), TARGET    :: graupel_in
    CLASS(particle), POINTER :: graupel ! ACCWA (nvhpc 22.7, IPSF, see above)

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: self_n

    IF (isdebug) CALL message(routine, "graupel_selfcollection")

    graupel => graupel_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(q_g, n_g, x_g, d_g, v_g, self_n)
    DO k = kstart,kend
       DO i = istart,iend

          q_g = graupel%q(i,k)
          IF ( q_g > q_crit ) THEN

             n_g = graupel%n(i,k)
             x_g = particle_meanmass(graupel, q_g,n_g)
             d_g = particle_diameter(graupel, x_g)
             v_g = particle_velocity(graupel, x_g) * graupel%rho_v(i,k)

             self_n = graupel_coeffs%sc_coll_n * n_g**2 * D_g**2 * v_g * dt

             ! sticking efficiency does only distinguish dry and wet based on T_3
             self_n = self_n * MERGE(cfg_params%ecoll_gg_wet, cfg_params%ecoll_gg, atmo%T(i,k) > cfg_params%Tcoll_gg_wet )
!             self_n = self_n * MERGE(cfg_params%ecoll_gg_wet, cfg_params%ecoll_gg, atmo%T(i,k) > 270.16 )
!             self_n = self_n * MERGE(cfg_params%ecoll_gg_wet, cfg_params%ecoll_gg, atmo%T(i,k) > T_3)
             self_n = MIN(self_n,n_g)

             graupel%n(i,k) = graupel%n(i,k) - self_n
          ENDIF
       ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE graupel_selfcollection

  SUBROUTINE ice_melting(ik_slice, atmo, ice_in, cloud, rain)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(inout), TARGET :: ice_in
    CLASS(particle), POINTER :: ice ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle), INTENT(inout) :: cloud, rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_i,x_i,n_i
    REAL(wp)            :: melt_q,melt_n

    IF (isdebug) CALL message(routine, "ice_melting")

    ice => ice_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(q_i, n_i, x_i, melt_q, melt_n)
    DO k = kstart,kend
       DO i = istart,iend

          q_i = ice%q(i,k)

          IF (atmo%T(i,k) > T_3 .AND. q_i > 0.0) THEN

            n_i = ice%n(i,k)
            x_i = particle_meanmass(ice, q_i,n_i)

            ! complete melting within this time step
            melt_q = q_i
            melt_n = n_i
            ice%q(i,k) = 0.0_wp
            ice%n(i,k) = 0.0_wp

            ! ice either melts into cloud droplets or rain depending on x_i
            IF (x_i > cloud%x_max) THEN
               rain%q(i,k)  = rain%q(i,k)  + melt_q
               rain%n(i,k)  = rain%n(i,k)  + melt_n
            ELSE
               cloud%q(i,k) = cloud%q(i,k) + melt_q
               cloud%n(i,k) = cloud%n(i,k) + melt_n
            ENDIF

          END IF
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE ice_melting

  SUBROUTINE particle_cloud_riming(ik_slice, dt, atmo, ptype_in, coeffs, cloud_in, rain_in, ice_in, &
       shed_coeffs, snow_in, ltabdminwgp, &
       shed_ltab_dpp_03, shed_ltab_dpc_02, shed_ltab_dcc_01, &
       shed_ltab_tpp_05, shed_ltab_tpp_03, shed_ltab_tpc_04)
    !*******************************************************************************
    ! Riming of graupel or hail with cloud droplets                                *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt

    TYPE(collection_coeffs), INTENT(in)   :: coeffs
    TYPE(atmosphere), INTENT(inout)       :: atmo
    CLASS(particle_frozen), INTENT(inout), TARGET :: ice_in
    CLASS(particle), INTENT(inout), TARGET :: cloud_in, rain_in
    CLASS(particle_frozen), INTENT(inout), TARGET :: ptype_in

    ! coefficients, add. hydrometeors, incomplete gamma functions
    !  and wet growht LUT for OPTIONAL droplet shedding:
    TYPE(coll_coeffs_ir_pm), INTENT(in), OPTIONAL :: shed_coeffs
    TYPE(gamlookuptable), INTENT(in), OPTIONAL :: &
       shed_ltab_dpp_03, shed_ltab_dpc_02, shed_ltab_dcc_01, &
       shed_ltab_tpp_05, shed_ltab_tpp_03, shed_ltab_tpc_04
    CLASS(particle_frozen), INTENT(in), TARGET, OPTIONAL :: snow_in
    TYPE(lookupt_4d), INTENT(in), TARGET, OPTIONAL :: ltabdminwgp

    ! UB: why particle_frozen and not particle like in particle rain riming?
    CLASS(particle_frozen), POINTER :: ptype       ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle_frozen), POINTER :: ice, snow   ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle)       , POINTER :: cloud, rain ! ACCWA (nvhpc 22.7, IPSF, see above)
    TYPE(lookupt_4D)      , POINTER :: p_ltabdminwgp

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend

    INTEGER             :: i,k
    REAL(wp)            :: T_a, p_a
    REAL(wp)            :: q_p,n_p,x_p,d_p,v_p
    REAL(wp)            :: q_c,n_c,x_c,d_c,v_c
    REAL(wp)            :: q_r,n_r,x_r,q_i,q_s,d_wg,d_shed_p,x_shed_p,lam_p,um_p
    REAL(wp)            :: rime_n,rime_q,e_coll_n
    REAL(wp)            :: melt_n,melt_q,e_coll_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp)            :: const1
    REAL(wp), PARAMETER :: &
         const0 = 1.0/(D_coll_c - D_crit_c),     &
         const2 = 1.0/(T_mult_opt - T_mult_min), &
         const3 = 1.0/(T_mult_opt - T_mult_max), &
         const4 = c_w / L_ew

    LOGICAL  :: shedding_enabled
    ! shedding parameterization:
    REAL(wp)            :: shed_q, shed_n, vchar, nenner
    REAL(wp)            :: delta_aa_var, delta_ab_var, delta_bb_var, theta_aa_var, theta_ab_var, theta_bb_var
    REAL(wp), PARAMETER :: T_shed = 263.15_wp
    REAL(wp), PARAMETER :: D_shedding = 1000.0e-6_wp     !..mittlerer Durchmesser Shedding
    REAL(wp), PARAMETER :: x_shed = pi/6.0_wp * rho_w * D_shedding**3

    IF (isdebug) CALL message(routine, "particle_cloud_riming")

    ! itype_shedding = 2 can be computed if the optional parameters are present. Check a few of them:
    shedding_enabled = PRESENT(shed_coeffs) .AND. PRESENT(ltabdminwgp) .AND. PRESENT(shed_ltab_dpp_03)
    
    cloud => cloud_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    rain  => rain_in  ! ACCWA (nvhpc 22.7, IPSF, see above)
    ptype => ptype_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    ice   => ice_in
    IF (shedding_enabled) THEN
      p_ltabdminwgp => ltabdminwgp
      snow => snow_in
    END IF

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    const1 = const0 * ptype%ecoll_c

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, p_a, q_c, q_p, n_c, n_p, x_p, D_p, x_c, D_c) &
    !$ACC   PRIVATE(q_r, n_r, x_r, q_i, q_s, d_wg, d_shed_p, x_shed_p, lam_p, um_p) &
    !$ACC   PRIVATE(shed_q, shed_n, vchar, nenner) &
    !$ACC   PRIVATE(v_p, v_c, e_coll_n, e_coll_q, rime_n, rime_q) &
    !$ACC   PRIVATE(mult_1, mult_2, mult_n, mult_q, melt_n, melt_q) &
    !$ACC   PRIVATE(delta_aa_var, delta_ab_var, delta_bb_var, theta_aa_var, theta_ab_var, theta_bb_var)
    DO k = kstart,kend
      DO i = istart,iend

        q_c = cloud%q(i,k)
        q_p = ptype%q(i,k)
        n_c = cloud%n(i,k)
        n_p = ptype%n(i,k)

        x_p = particle_meanmass(ptype,q_p,n_p)
        D_p = particle_diameter(ptype,x_p)
        x_c = particle_meanmass(cloud,q_c,n_c)
        D_c = particle_diameter(cloud,x_c)

        T_a = atmo%T(i,k)
        IF (q_c > q_crit_c .AND. q_p > ptype%q_crit_c &
             &             .AND. D_p > ptype%D_crit_c .AND. D_c > D_crit_c) THEN

          v_p = particle_velocity(ptype,x_p) * ptype%rho_v(i,k)
          v_c = particle_velocity(cloud,x_c) * cloud%rho_v(i,k)

          e_coll_n = MIN(ptype%ecoll_c, MAX(const1*(D_c - D_crit_c),ecoll_min))
          e_coll_q = e_coll_n

          rime_n = pi4 * e_coll_n * n_p * n_c * dt &
               & *     (coeffs%delta_n_aa * D_p**2 + coeffs%delta_n_ab * D_p*D_c + coeffs%delta_n_bb * D_c**2) &
               & * SQRT(coeffs%theta_n_aa * v_p**2 - coeffs%theta_n_ab * v_p*v_c + coeffs%theta_n_bb * v_c**2)

          rime_q = pi4 * e_coll_q * n_p * q_c * dt &
               & *     (coeffs%delta_q_aa * D_p**2 + coeffs%delta_q_ab * D_p*D_c + coeffs%delta_q_bb * D_c**2) &
               & * SQRT(coeffs%theta_q_aa * v_p**2 - coeffs%theta_q_ab * v_p*v_c + coeffs%theta_q_bb * v_c**2)

          rime_q = MIN(q_c,rime_q)
          rime_n = MIN(n_c,rime_n)

          ptype%q(i,k) = ptype%q(i,k) + rime_q
          cloud%q(i,k) = cloud%q(i,k) - rime_q
          cloud%n(i,k) = cloud%n(i,k) - rime_n

          ! ice multiplication based on Hallet and Mossop
          IF (T_a < T_3 .AND. ice_multiplication) THEN
            mult_1 = const2*(T_a - T_mult_min)
            mult_2 = const3*(T_a - T_mult_max)
            mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
            mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
            mult_n = C_mult * mult_1 * mult_2 * rime_q
            mult_q = mult_n * ice%x_min
            mult_q = MIN(rime_q,mult_q)
            mult_n = mult_q / ice%x_min

            ice%n(i,k)   = ice%n(i,k)   + mult_n
            ice%q(i,k)   = ice%q(i,k)   + mult_q
            ptype%q(i,k) = ptype%q(i,k) - mult_q
          ENDIF

          ! enhancement of melting 
          IF (T_a > T_3 .AND. enhanced_melting) THEN
            melt_q = const4 * (T_a - T_3) * rime_q
            melt_n = melt_q / x_p

            melt_q = MIN(ptype%q(i,k),melt_q)
            melt_n = MIN(ptype%n(i,k),melt_n)

            ptype%q(i,k) = ptype%q(i,k) - melt_q
            rain%q(i,k)  = rain%q(i,k)  + melt_q
            ptype%n(i,k) = ptype%n(i,k) - melt_n
            rain%n(i,k)  = rain%n(i,k)  + melt_n
          ENDIF
          
          ! Shedding:
          IF (shedding_enabled .OR. cfg_params%itype_shedding_gh == 1) THEN
            
            SELECT CASE (cfg_params%itype_shedding_gh)
            CASE (1)
              ! simple version from COSMO: Complete shedding for T > T_3, else shedding if particle
              !  mean mass diameter is larger than D_shed_gh:
              IF ( (D_p > cfg_params%D_shed_gh .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN

                q_r = rain%q(i,k)
                n_r = rain%n(i,k)
                x_r = particle_meanmass(rain,q_r,n_r)

                shed_q = MIN(ptype%q(i,k),rime_q)
                IF (T_a <= T_3) THEN
                  shed_n = shed_q / MIN(x_shed,x_p)
                ELSE
                  shed_n = shed_q / MAX(x_r,x_p)
                ENDIF
                
                ptype%q(i,k) = ptype%q(i,k) - shed_q
                rain %q(i,k) = rain %q(i,k) + shed_q
                rain %n(i,k) = rain %n(i,k) + shed_n
              ENDIF
              
            CASE (2)
              
              ! more physical solution involving the wet growth diameter and upper incomplete gamma function:
              IF (T_a > T_shed) THEN

                p_a = atmo%p(i,k)
                q_r = rain%q(i,k)
                n_r = rain%n(i,k)
                x_r = particle_meanmass(rain,q_r,n_r)

                q_i = ice%q(i,k)
                q_s = snow%q(i,k)

                d_wg = dmin_wg_gr_ltab_equi(p_a,T_a,q_c+q_r,q_i+q_s,p_ltabdminwgp)
                ! this function sets d_wg=0.0 for T > T_3

                ! Shedding occurs only for cloud drops collected from particles > D_shed, which have a wet surface.
                ! Wet surface can happen in conditions of wet growth or for T > T_3.
                ! For single particles, D_shed is theoretically at least the critical diameter
                !  of about 9 mm found by Rasmussen und Heymsfield (1987), below which the water coating on
                !  an ice clump is aerodynamically stable:
                d_shed_p = MAX(d_wg, cfg_params%D_shed_gh)
                ! We use a fixed namelist parameter cfg_params%D_shed_gh to mimick RH87 stability
                !  regardless of environmental temperature T_a,
                !  but with the possibility of tuning to compensate for the effect that
                !  a bulk scheme re-creates the large tail of the PSD in every time step,
                !  which is a kind of numerical diffusion in size space.

                x_shed_p = particle_mass(ptype, d_shed_p)
                lam_p = shed_coeffs%lamfakt_a * x_p**(-ptype%mu)
                um_p  = lam_p * x_shed_p**ptype%mu

                ! mass of shed drops, which have been collected by particles > D_shed_p:
                delta_aa_var = incgfct_upper_lookup( um_p, shed_ltab_dpp_03)
                delta_ab_var = incgfct_upper_lookup( um_p, shed_ltab_dpc_02)
                delta_bb_var = incgfct_upper_lookup( um_p, shed_ltab_dcc_01)
                theta_aa_var = incgfct_upper_lookup( um_p, shed_ltab_tpp_05)
                theta_ab_var = incgfct_upper_lookup( um_p, shed_ltab_tpc_04)
                theta_bb_var = 1.0_wp
                nenner = incgfct_upper_lookup( um_p, shed_ltab_tpp_03)
                ! The nenner should not be 0 for mathematical reasons. Limit it to the second last
                !  table node. The table is for the lower incgft, so we have to take
                !  the difference to the uppermost table value which is the ordinary gamma function:
                nenner = MAX(nenner, shed_ltab_tpp_03%igf(shed_ltab_tpp_03%n) - &
                                     shed_ltab_tpp_03%igf(shed_ltab_tpp_03%n-1)) 

                vchar = shed_coeffs%theta_aa(0,1) * theta_aa_var / nenner * v_p**2 - &
                        shed_coeffs%theta_ab(0,1) * theta_ab_var / nenner * v_p*v_c + &
                        shed_coeffs%theta_bb(0,1) * theta_bb_var * v_c**2

                shed_q = pi4 * e_coll_q * n_p * q_c * dt * &
                     ( shed_coeffs%delta_aa(0,1) * delta_aa_var * D_p**2 + &
                       shed_coeffs%delta_ab(0,1) * delta_ab_var * D_p*D_c + &
                       shed_coeffs%delta_bb(0,1) * delta_bb_var * D_c**2 ) * &
                     SQRT( MAX(vchar, 0.0d0) )

                shed_q = MIN(ptype%q(i,k),shed_q)
                shed_n = shed_q / x_shed
                
                ptype%q(i,k) = ptype%q(i,k) - shed_q
                rain %q(i,k) = rain %q(i,k) + shed_q
                rain %n(i,k) = rain %n(i,k) + shed_n

              END IF
  
            END SELECT
            
          END IF

        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE particle_cloud_riming

  SUBROUTINE particle_rain_riming(ik_slice, dt, atmo, ptype_in, &
       coeffs, rain_in, ice_in, &
       shed_coeffs, snow_in, cloud_in, ltabdminwgp, &
       shed_ltab_dpp_03, shed_ltab_dpr_02, shed_ltab_drr_01, &
       shed_ltab_tpp_05, shed_ltab_tpp_03, shed_ltab_tpr_04)
    !*******************************************************************************
    ! Riming of graupel or hail with rain drops including shedding                 *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    TYPE(collection_coeffs), INTENT(in) :: coeffs
    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(inout), TARGET :: ice_in
    CLASS(particle), INTENT(inout), TARGET :: rain_in, ptype_in
    
    ! coefficients, add. hydrometeors, incomplete gamma functions
    !  and wet growht LUT for OPTIONAL droplet shedding:
    TYPE(coll_coeffs_ir_pm), INTENT(in), OPTIONAL :: shed_coeffs
    TYPE(gamlookuptable), INTENT(in), OPTIONAL :: &
       shed_ltab_dpp_03, shed_ltab_dpr_02, shed_ltab_drr_01, &
       shed_ltab_tpp_05, shed_ltab_tpp_03, shed_ltab_tpr_04
    CLASS(particle), INTENT(in), TARGET, OPTIONAL :: cloud_in, snow_in
    TYPE(lookupt_4d), INTENT(in), TARGET, OPTIONAL :: ltabdminwgp

    CLASS(particle), POINTER :: rain, ptype, ice, cloud, snow ! ACCWA (nvhpc 22.7, IPSF, see above)
    TYPE(lookupt_4D), POINTER            :: p_ltabdminwgp
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a,p_a
    REAL(wp)            :: q_p,n_p,x_p,d_p,v_p
    REAL(wp)            :: q_r,n_r,x_r,d_r,v_r
    REAL(wp)            :: q_c,q_i,q_s,d_wg,d_shed_p,x_shed_p,lam_p,um_p
    REAL(wp)            :: rime_n,rime_q,melt_n,melt_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp), PARAMETER :: &
         const2 = 1/(T_mult_opt - T_mult_min), &
         const3 = 1/(T_mult_opt - T_mult_max), &
         const4 = c_w / L_ew

    LOGICAL  :: shedding_enabled
    ! shedding parameterization:
    REAL(wp)            :: shed_q, shed_n, vchar, nenner
    REAL(wp)            :: delta_aa_var, delta_ab_var, delta_bb_var, theta_aa_var, theta_ab_var, theta_bb_var
    REAL(wp), PARAMETER :: T_shed = 263.15_wp
    REAL(wp), PARAMETER :: D_shedding = 1000.0e-6_wp     !..mittlerer Radius Shedding
    REAL(wp), PARAMETER :: x_shed = pi/6.0_wp * rho_w * D_shedding**3

    IF (isdebug) CALL message(routine, "particle_rain_riming")

    ! itype_shedding = 2 can be computed if the optional parameters are present. Check a few of them:
    shedding_enabled = PRESENT(shed_coeffs) .AND. PRESENT(ltabdminwgp) .AND. PRESENT(shed_ltab_dpp_03)
    
    rain => rain_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    ptype => ptype_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    ice => ice_in
    IF (shedding_enabled) THEN
      p_ltabdminwgp => ltabdminwgp
      cloud => cloud_in
      snow => snow_in
    END IF

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, p_a, q_r, q_p, n_r, n_p, x_p, d_p, v_p, x_r, D_r) &
    !$ACC   PRIVATE(v_r, rime_n, rime_q, mult_1, mult_2, mult_n, mult_q, vchar, nenner) &
    !$ACC   PRIVATE(melt_q, melt_n, shed_q, shed_n, d_wg, d_shed_p, x_shed_p, q_c, q_i, q_s, lam_p, um_p) &
    !$ACC   PRIVATE(delta_aa_var, delta_ab_var, delta_bb_var, theta_aa_var, theta_ab_var, theta_bb_var)
    DO k = kstart,kend
      DO i = istart,iend

        T_a = atmo%T(i,k)
        q_r = rain%q(i,k)
        q_p = ptype%q(i,k)

        IF (q_r > q_crit .AND. q_p > q_crit) THEN

          n_r = rain%n(i,k)
          n_p = ptype%n(i,k)

          x_p = particle_meanmass(ptype, q_p,n_p)
          D_p = particle_diameter(ptype, x_p)
          v_p = particle_velocity(ptype, x_p) * ptype%rho_v(i,k)
          x_r = particle_meanmass(rain, q_r,n_r)
          D_r = particle_diameter(rain, x_r)
          v_r = particle_velocity(rain, x_r) * rain%rho_v(i,k)

          rime_n = pi4 * n_p * n_r * dt &
               & *     (coeffs%delta_n_aa * D_p**2 + coeffs%delta_n_ab * D_p*D_r + coeffs%delta_n_bb * D_r**2) &
               & * SQRT(coeffs%theta_n_aa * v_p**2 - coeffs%theta_n_ab * v_p*v_r + coeffs%theta_n_bb * v_r**2)

          rime_q = pi4 * n_p * q_r * dt &
               & *     (coeffs%delta_n_aa * D_p**2 + coeffs%delta_q_ab * D_p*D_r + coeffs%delta_q_bb * D_r**2) &
               & * SQRT(coeffs%theta_n_aa * v_p**2 - coeffs%theta_q_ab * v_p*v_r + coeffs%theta_q_bb * v_r**2)

          rime_q = MIN(q_r,rime_q)
          rime_n = MIN(n_r,rime_n)

          ptype%q(i,k) = ptype%q(i,k) + rime_q
          rain%q(i,k)    = rain%q(i,k)    - rime_q
          rain%n(i,k)    = rain%n(i,k)    - rime_n

          ! ice multiplication based on Hallet and Mossop
          IF (T_a < T_3 .AND. ice_multiplication) THEN
            mult_1 = (T_a - T_mult_min) * const2
            mult_2 = (T_a - T_mult_max) * const3
            mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
            mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
            mult_n = C_mult * mult_1 * mult_2 * rime_q
            mult_q = mult_n * ice%x_min
            mult_q = MIN(rime_q,mult_q)

            ice%n(i,k)   = ice%n(i,k)   + mult_n
            ice%q(i,k)   = ice%q(i,k)   + mult_q
            ptype%q(i,k) = ptype%q(i,k) - mult_q
          ENDIF

          ! enhancement of melting of ptype
          IF (T_a > T_3 .AND. enhanced_melting) THEN
            melt_q = const4 * (T_a - T_3) * rime_q
            melt_n = melt_q / x_p

            melt_q = MIN(ptype%q(i,k),melt_q)
            melt_n = MIN(ptype%n(i,k),melt_n)

            ptype%q(i,k) = ptype%q(i,k) - melt_q
            rain%q(i,k)  = rain%q(i,k)  + melt_q

            ptype%n(i,k) = ptype%n(i,k) - melt_n
            rain%n(i,k)  = rain%n(i,k)  + melt_n
          ENDIF

          ! Shedding:
          IF (shedding_enabled .OR. cfg_params%itype_shedding_gh == 1) THEN
            
            SELECT CASE (cfg_params%itype_shedding_gh)
            CASE (1)
              ! simple version from COSMO: Complete shedding for T > T_3, else shedding if particle
              !  mean mass diameter is larger than D_shed_gh:
              IF ( (D_p > cfg_params%D_shed_gh .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN

                shed_q = MIN(ptype%q(i,k),rime_q)
                IF (T_a <= T_3) THEN
                  shed_n = shed_q / MIN(x_shed,x_p)
                ELSE
                  shed_n = shed_q / MAX(x_r,x_p)
                ENDIF
                
                ptype%q(i,k) = ptype%q(i,k) - shed_q
                rain %q(i,k) = rain %q(i,k) + shed_q
                rain %n(i,k) = rain %n(i,k) + shed_n
              ENDIF
              
            CASE (2)
              
              ! more physical solution involving the wet growth diameter and upper incomplete gamma function:
              IF (T_a > T_shed) THEN

                p_a = atmo%p(i,k)
                q_c = cloud%q(i,k)
                q_i = ice%q(i,k)
                q_s = snow%q(i,k)

                d_wg = dmin_wg_gr_ltab_equi(p_a,T_a,q_c+q_r,q_i+q_s,p_ltabdminwgp)
                ! this function sets d_wg=0.0 for T > T_3

                ! Shedding occurs only for cloud drops collected from particles > D_shed, which have a wet surface.
                ! Wet surface can happen in conditions of wet growth or for T > T_3.
                ! For single particles, D_shed is theoretically at least the critical diameter
                !  of about 9 mm found by Rasmussen und Heymsfield (1987), below which the water coating on
                !  an ice clump is aerodynamically stable:
                d_shed_p = MAX(d_wg, cfg_params%D_shed_gh)
                ! We use a fixed namelist parameter cfg_params%D_shed_gh to mimick RH87 stability
                !  regardless of environmental temperature T_a,
                !  but with the possibility of tuning to compensate for the effect that
                !  a bulk scheme re-creates the large tail of the PSD in every time step,
                !  which is a kind of numerical diffusion in size space.

                x_shed_p = particle_mass(ptype, d_shed_p)
                lam_p = shed_coeffs%lamfakt_a * x_p**(-ptype%mu)
                um_p  = lam_p * x_shed_p**ptype%mu

                ! mass of shed drops, which have been collected by particles > D_shed_p:
                delta_aa_var = incgfct_upper_lookup( um_p, shed_ltab_dpp_03)
                delta_ab_var = incgfct_upper_lookup( um_p, shed_ltab_dpr_02)
                delta_bb_var = incgfct_upper_lookup( um_p, shed_ltab_drr_01)
                theta_aa_var = incgfct_upper_lookup( um_p, shed_ltab_tpp_05)
                theta_ab_var = incgfct_upper_lookup( um_p, shed_ltab_tpr_04)
                theta_bb_var = 1.0_wp
                nenner = incgfct_upper_lookup( um_p, shed_ltab_tpp_03)
                ! The nenner should not be 0 for mathematical reasons. Limit it to the second last
                !  table node. The table is for the lower incgft, so we have to take
                !  the difference to the uppermost table value which is the ordinary gamma function:
                nenner = MAX(nenner, shed_ltab_tpp_03%igf(shed_ltab_tpp_03%n) - &
                                     shed_ltab_tpp_03%igf(shed_ltab_tpp_03%n-1)) 

                vchar = shed_coeffs%theta_aa(0,1) * theta_aa_var / nenner * v_p**2 - &
                        shed_coeffs%theta_ab(0,1) * theta_ab_var / nenner * v_p*v_r + &
                        shed_coeffs%theta_bb(0,1) * theta_bb_var * v_r**2

                shed_q = pi4 * n_p * q_r * dt * &
                     ( shed_coeffs%delta_aa(0,1) * delta_aa_var * D_p**2 + &
                       shed_coeffs%delta_ab(0,1) * delta_ab_var * D_p*D_r + &
                       shed_coeffs%delta_bb(0,1) * delta_bb_var * D_r**2 ) * &
                     SQRT( MAX(vchar, 0.0d0) )

                shed_q = MIN(ptype%q(i,k),shed_q)
                shed_n = shed_q / x_shed
                
                ptype%q(i,k) = ptype%q(i,k) - shed_q
                rain %q(i,k) = rain %q(i,k) + shed_q
                rain %n(i,k) = rain %n(i,k) + shed_n

              END IF
  
            END SELECT
            
          END IF
          
        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE particle_rain_riming

  SUBROUTINE graupel_melting(ik_slice, dt, graupel_coeffs, atmo, graupel_in, rain)
    !*******************************************************************************
    ! Melting of graupel                                                           *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    CLASS(particle_sphere), INTENT(in) :: graupel_coeffs
    TYPE(atmosphere), INTENT(inout)    :: atmo
    CLASS(particle), INTENT(inout)     :: rain    
    CLASS(particle), INTENT(inout), TARGET :: graupel_in
    CLASS(particle), POINTER :: graupel ! ACCWA (nvhpc 22.7, IPSF, see above)

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g,T_a,e_a
    REAL(wp)            :: melt,melt_v,melt_h,melt_n,melt_q
    REAL(wp)            :: fh_q,fv_q

    IF (isdebug) CALL message(routine, "graupel_melting")

    graupel => graupel_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    ! ACCWA (nvhpc 21.3): Explicit PRESENT statement required for graupel:
    ! - Otherwise error during runtime: variable in data clause is partially present on device: name=descriptor
    ! - Reason unknown as inlining of particle_* solves issue
    !$ACC DATA PRESENT(graupel)
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, q_g, e_a, n_g, x_g, D_g, v_g) &
    !$ACC   PRIVATE(fv_q, fh_q, melt, melt_h, melt_v, melt_q, melt_n)
    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)
          q_g = graupel%q(i,k)

          IF (T_a > T_3 .AND. q_g > 0.0_wp) THEN
!             e_a = e_ws(T_a)
             e_a  = atmo%qv(i,k) * R_d * T_a  ! need actual e_v here, not satur. vapor pres.
             n_g = graupel%n(i,k)

             x_g = particle_meanmass(graupel, q_g,n_g)
             D_g = particle_diameter(graupel, x_g)
             v_g = particle_velocity(graupel, x_g) * graupel%rho_v(i,k)

             fv_q = graupel_coeffs%a_f + graupel_coeffs%b_f * SQRT(v_g*D_g)
             fh_q = 1.05_wp * fv_q

             melt   = 2.0_wp*pi / L_ew * D_g * n_g * dt

             melt_h = melt * K_T * (T_a - T_3)
             melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)

             melt_q = cfg_params%melt_g_tune_fak* (melt_h * fh_q + melt_v * fv_q)

             ! UB: assume that x_g is constant during melting
             melt_n = MIN(MAX( (melt_q - q_g) / x_g + n_g, 0.0_wp), n_g)

             melt_q = MIN(q_g,melt_q)
             melt_n = MIN(n_g,melt_n)
             melt_q = MAX(0.0_wp,melt_q)
             melt_n = MAX(0.0_wp,melt_n)

             graupel%q(i,k) = graupel%q(i,k) - melt_q
             rain%q(i,k)    = rain%q(i,k)    + melt_q

             graupel%n(i,k) = graupel%n(i,k) - melt_n
             rain%n(i,k)    = rain%n(i,k)    + melt_n

          ENDIF
       ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer
    !$ACC END DATA

  END SUBROUTINE graupel_melting

!!$  INTERFACE hail_melting
!!$     SUBROUTINE hail_melting_simple(ik_slice, dt, atmo, hail, rain)
!!$       INTEGER,  INTENT(in) :: ik_slice(4)
!!$       REAL(wp), INTENT(in) :: dt
!!$       TYPE(atmosphere), INTENT(inout) :: atmo
!!$       CLASS(particle), INTENT(inout)  :: rain
!!$       TYPE(particle), INTENT(inout)   :: hail
!!$     END SUBROUTINE hail_melting_simple
!!$
!!$     SUBROUTINE hail_melting_lwf(ik_slice, dt, atmo, hail, rain)
!!$       INTEGER, INTENT(in) :: ik_slice(4)
!!$       REAL(wp), INTENT(in) :: dt
!!$       TYPE(atmosphere), INTENT(inout)   :: atmo
!!$       CLASS(particle),  INTENT(inout)   :: rain
!!$       TYPE(particle_lwf), INTENT(inout) :: hail
!!$     END SUBROUTINE HAIL_MELTING_LWF
!!$  END INTERFACE hail_melting

  SUBROUTINE hail_melting_simple(ik_slice, dt, hail_coeffs, atmo, hail_in, rain)
    !*******************************************************************************
    ! Melting of hail                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    TYPE(particle_sphere), INTENT(in) :: hail_coeffs
    TYPE(atmosphere), INTENT(inout)   :: atmo
    CLASS(particle), INTENT(inout)    :: rain
    CLASS(particle), INTENT(inout), TARGET:: hail_in
    CLASS(particle), POINTER :: hail ! ACCWA (nvhpc 22.7, IPSF, see above)
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    ! local variables
    INTEGER             :: i,k
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h,T_a,e_a
    REAL(wp)            :: melt,melt_v,melt_h,melt_n,melt_q
    REAL(wp)            :: fh_q,fv_q

    IF (isdebug) CALL message(routine, "hail_melting")

    hail => hail_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)


    ! ACCWA (nvhpc 21.3): Explicit PRESENT statement required for hail:
    ! - Otherwise error during runtime: variable in data clause is partially present on device: name=descriptor
    ! - Reason unknown as inlining of particle_* solves issue
    !$ACC DATA PRESENT(hail)
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, q_h, e_a, n_h, x_h, D_h, v_h) &
    !$ACC   PRIVATE(fv_q, fh_q, melt, melt_h, melt_v, melt_q, melt_n)
    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)
          q_h = hail%q(i,k)

          IF (T_a > T_3 .AND. q_h > 0.0_wp) THEN
!            e_a = e_ws(T_a)
            e_a  = atmo%qv(i,k) * R_d * T_a
            n_h = hail%n(i,k)

            x_h = particle_meanmass(hail,q_h,n_h)
            D_h = particle_diameter(hail,x_h)
            v_h = particle_velocity(hail,x_h) * hail%rho_v(i,k)

            fv_q = hail_coeffs%a_f + hail_coeffs%b_f * sqrt(v_h*D_h)
            fh_q = 1.05_wp * fv_q                            ! UB: based on Rasmussen and Heymsfield

            melt   = 2.0_wp*pi / L_ew * D_h * n_h * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)

            melt_q = cfg_params%melt_h_tune_fak * (melt_h * fh_q + melt_v * fv_q)

            ! UB: assume that x_h is constant during melting
            melt_n = MIN(MAX( (melt_q - q_h) / x_h + n_h, 0.0_wp), n_h)

            melt_q = MIN(q_h,melt_q)
            melt_n = MIN(n_h,melt_n)
            melt_q = MAX(0.0_wp,melt_q)
            melt_n = MAX(0.0_wp,melt_n)

            hail%q(i,k) = hail%q(i,k) - melt_q
            rain%q(i,k) = rain%q(i,k) + melt_q

            hail%n(i,k) = hail%n(i,k) - melt_n
            rain%n(i,k) = rain%n(i,k) + melt_n

         ENDIF
      ENDDO
   ENDDO
   !$ACC END PARALLEL
   !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer
   !$ACC END DATA

  END SUBROUTINE hail_melting_simple

  SUBROUTINE particle_melting_lwf(ik_slice, dt, ptype, rain, gta)
    !*******************************************************************************
    ! Melting of hail for lwf scheme                                               *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt

    ! 2mom variables
    CLASS(particle),  INTENT(inout)      :: rain
    TYPE(particle_lwf), INTENT(inout)    :: ptype
    REAL(wp), INTENT(IN), DIMENSION(:,:) :: gta   ! Thermodynamic environment function

    ! local variables
    INTEGER   :: i,k
    REAL(wp)  :: q_p,n_p,x_p,d_p,d_n,qliq,qice,lwf,qrain,nrain
    REAL(wp)  :: x_shed,D_star
    REAL(wp)  :: qmlt,qconv,nconv,qfrac,nfrac,qshed
    REAL(wp)  :: D_lim,lwf_lim,cmlt1,cmlt2,max_qmlt,melt

    ! local parameters
    REAL(wp), PARAMETER :: D_shed = 1.e-3_wp

    ! local arrays
    REAL(wp), DIMENSION(10) :: adst, amlt
    REAL(wp), DIMENSION(9)  :: bdst, bmlt
    REAL(wp), DIMENSION(4)  :: cnvq, cnvn

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend

    IF (isdebug) CALL message(routine, "particle_melting_lwf")

#ifdef _OPENACC
    CALL finish("particle_melting_lwf", "Routine has not been ported to openACC yet.")
#endif

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    if (ptype%name .eq. 'hail_vivek') then
      adst = adstarh
      bdst = bdstarh
      amlt = amelth
      bmlt = bmelth
      cnvq = convqh
      cnvn = convnh
      D_lim = 8.e-3_wp
      lwf_lim = 0.3_wp
    elseif (ptype%name .eq. 'graupel_vivek') then
      adst = adstarg
      bdst = bdstarg
      amlt = ameltg
      bmlt = bmeltg
      cnvq = convqg
      cnvn = convng
      D_lim = 1.e-2_wp
      lwf_lim = 0.6_wp
    else
      cnvq = convqg  ! avoid compiler warning
      cnvn = convng
      D_lim   = 0.0_wp
      lwf_lim = 0.0_wp
      CALL finish(TRIM(routine),'Error unknown particle name in LWF melting')
    end if

    cmlt1 = ptype%lwf_cmelt1
    cmlt2 = ptype%lwf_cmelt2

    x_shed = (D_shed/rain%a_geo)**(1._wp/rain%b_geo)

    DO k = kstart,kend
      DO i = istart,iend
        
        n_p  = ptype%n(i,k)
        q_p  = ptype%q(i,k)
        qliq = ptype%l(i,k)
        qice = ptype%q(i,k) - ptype%l(i,k)

        max_qmlt = qice / dt

        qrain = rain%q(i,k)
        nrain = rain%n(i,k)

        qmlt  = 0.0_wp
        qconv = 0.0_wp
        nconv = 0.0_wp
        qshed = 0._wp

        !! NOTE: LWF has to be given an upper limit!
        !!       Something like 0.85, in which case all meltwater is transfered!
        
        if (qice > 0.0_wp .and. gta(i,k) > 0.0_wp) then
 
          x_p = particle_meanmass(ptype,q_p,n_p)
          D_p = particle_diameter(ptype,x_p)
          D_n = particle_normdiameter(ptype,D_p)
          lwf = particle_lwf_idx(ptype,i,k)

          D_star = 1.e-2_wp * rat2do3(D_n,lwf,adst,bdst)  ! 1e-2 for conversion from cm to m

          melt = rat2do3(D_n,lwf,amlt,bmlt)
          !! UB: exp(log()) crashes for very small values of qice!
          !! Either include a "security-eps" or revert to original power-function 
          !! qmlt = gta(i,k) * melt * cmlt1 * exp( cmlt2*log(qliq+qice) )
          qmlt = gta(i,k) * melt * cmlt1 * (qliq+qice)**cmlt2
          qmlt = max( min(qmlt,max_qmlt), 0.0_wp)

          if (qmlt.eq.max_qmlt) then
            qconv = qliq/dt + qmlt
            nconv = n_p/dt
            qshed = 0.0_wp
          else

            ! conversion of meltwater to rain
            if (classic_melting_in_lwf_scheme .OR. lwf > 0.85_wp) then
              qfrac = 1.0_wp
            else
              qfrac = cnvq(1)*lwf**cnvq(2) + cnvq(3)*lwf**cnvq(4)
              qfrac = max(qfrac,0.0_wp)
            end if
            nfrac = cnvn(1)*lwf**cnvn(2) + cnvn(3)*lwf**cnvn(4)
            nfrac = max(nfrac,0.0_wp)

            qconv = qfrac * (qmlt + qliq/dt)
            nconv = nfrac * n_p/dt

            ! shedding of meltwater
            if (D_p.gt.D_lim .and. lwf.gt.lwf_lim) then
              qshed = max(lwf-0.3_wp, 0.0_wp) * qconv
            else
              qshed = 0.0_wp
            end if
            
          end if

          ! Forward integration
          n_p  = max(n_p - nconv*dt, 0.0_wp)
          qliq = max(qliq + (qmlt - qconv)*dt, 0.0_wp)
          qice = max(qice - qmlt*dt, 0.0_wp)

          qrain = qrain + (qconv + qshed)*dt
          nrain = nrain + (nconv + qshed/x_shed)*dt

          ptype%n(i,k) = n_p
          ptype%q(i,k) = qliq + qice
          ptype%l(i,k) = qliq
             
          rain%q(i,k) = qrain
          rain%n(i,k) = nrain

        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE particle_melting_lwf

  SUBROUTINE prepare_melting_lwf(ik_slice, atmo, gta)
    !*******************************************************************************
    ! Melting of hail for lwf scheme                                               *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)

    ! 2mom variables
    TYPE(atmosphere), INTENT(in) :: atmo
    REAL(wp), INTENT(OUT), DIMENSION(:,:) :: gta

    ! local variables
    INTEGER  :: i,k
    REAL(wp) :: T_a,p_a,rho_a,e_d,e_sw,rh_a,delta_T,delta_q,eta,Dv,ka,Le,Lm,nu,Sc,kt,Pr

    REAL(wp) :: e_ws_T3_o_T3  ! this could be a parameter

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend

    IF (isdebug) CALL message(routine, "hail_melting_lwf")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    e_ws_T3_o_T3 = e_ws(T_3)/T_3

    DO k = kstart,kend
      DO i = istart,iend

        p_a   = atmo%p(i,k)
        T_a   = atmo%T(i,k)
        rho_a = atmo%rho(i,k)

        e_d   = atmo%qv(i,k) * R_d * T_a
        e_sw  = e_ws(T_a)
        rh_a  = e_d / e_sw
    
        delta_T = T_a - T_3
        delta_q = rh_a * e_sw/T_a - e_ws_T3_o_T3

        IF (delta_q < 0.0_wp) delta_q = 0.0_wp
        
        IF  (delta_T > 0.0_wp .OR. delta_q > 0.0_wp) THEN

          eta = dyn_visc_sutherland(T_a)  !..dynamic viscosity        
          Dv  = Dv_Rasmussen(T_a,p_a)     !..diffusivity of water vapour
          ka  = ka_Rasmussen(T_a)         !..conductivity of air
          Le  = lh_evap_RH87(T_a)         !..latent heat of evaporation
          Lm  = lh_melt_RH87(T_a)         !..latent heat of melting

          nu  = eta/rho_a                 !..kinematic viscosity   
          Sc  = (nu/Dv)**(1.0_wp/3.0_wp)    !..Schmidt number
          kt  = ka / (rho_a*cp)           !..diffusivity of air    
          Pr  = (nu/kt)**(1.0_wp/3.0_wp)    !..Prandtl number    

          !.thermodynamic, i.e. environmental, function for melting
          gta(i,k) = 2.0_wp*pi/Lm * ( Pr*ka*delta_T + Sc*Le*Dv*delta_q/R_d )

        ELSE
          gta(i,k) = 0.0_wp
        END IF

      ENDDO
    ENDDO

  END SUBROUTINE prepare_melting_lwf

  SUBROUTINE graupel_hail_conv_wet_gamlook(ik_slice, graupel_ltable1, graupel_ltable2,       &
       &                                   graupel_nm1, graupel_nm2, graupel_g1, graupel_g2, &
       &                                   ltabdminwgg, atmo, graupel_in, cloud, rain, ice, snow, hail)
    !*******************************************************************************
    !  Wet growth and conversion of graupel to hail                                *
    !  (uses look-up table for incomplete gamma functions)                         *
    !  by Uli Blahak                                                               *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)

    ! Parameters and tables
    TYPE(gamlookuptable),INTENT(in) :: graupel_ltable1, graupel_ltable2
    REAL(wp), INTENT(in)            :: graupel_nm1, graupel_nm2, graupel_g1, graupel_g2
    TYPE(lookupt_4d), INTENT(in), TARGET    :: ltabdminwgg
    TYPE(lookupt_4D), POINTER       :: p_ltabdminwgg

    ! 2mom variables
    TYPE(atmosphere), INTENT(inout)       :: atmo
    CLASS(particle),  INTENT(inout)       :: cloud, rain
    CLASS(particle_frozen), INTENT(inout) :: ice, snow, hail    
    CLASS(particle_frozen), INTENT(inout), TARGET :: graupel_in
    CLASS(particle_frozen), POINTER :: graupel ! ACCWA (nvhpc 22.7, IPSF, see above)

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend

    ! local variables
    INTEGER             :: i, k
    REAL(wp)            :: T_a, p_a, d_trenn, qw_a, qi_a, n_0, lam, xmin
    REAL(wp)            :: q_g, n_g, x_g, d_g
    REAL(wp)            :: q_c, q_r
    REAL(wp)            :: conv_n, conv_q

    IF (isdebug) CALL message(routine, "graupel_hail_conv_wet_gamlook")

    graupel => graupel_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    p_ltabdminwgg => ltabdminwgg
    
    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    ! ACCWA (nvhpc 21.3): Explicit PRESENT statement required for graupel_ltable1 and graupel_ltable2:
    ! - Otherwise error: illegal address during kernel execution
    ! - Reason unknown as inlining of incgfct_upper_lookup solves issue
    !$ACC DATA PRESENT(graupel_ltable1, graupel_ltable2)
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, p_a, d_trenn, qw_a, qi_a, n_0, lam, xmin) &
    !$ACC   PRIVATE(q_g, n_g, x_g, d_g, q_c, q_r, conv_n, conv_q)
    DO k = kstart,kend
!NEC$ ivdep
       DO i = istart,iend

          q_c = cloud%q(i,k)
          q_r = rain%q(i,k)
          q_g = graupel%q(i,k)
          n_g = graupel%n(i,k)

          x_g = particle_meanmass(graupel, q_g,n_g)
          d_g = particle_diameter(graupel, x_g)
          n_g = q_g / x_g  ! for consistency for limiters, n_g is used explicitly below

          T_a = atmo%T(i,k)
          p_a = atmo%p(i,k)

          !..supercooled liquid water in the cloud environment = sum of rain and cloud water
          qw_a = q_r + q_c

          IF (T_a < T_3 .AND. q_g > graupel%q_crit_c .AND. qw_a > 1e-3_wp) THEN

            !.. Umgebungsgehalt Eispartikel (vernachl. werden Graupel und Hagel wg. geringer Kollisionseff.)
            !.. koennte problematisch sein, weil in konvekt. Wolken viel mehr Graupel und Hagel enthalten ist!!!
            qi_a = ice%q(i,k) + snow%q(i,k)

            IF (luse_dmin_wetgrowth_table) THEN
              d_trenn = dmin_wg_gr_ltab_equi(p_a,T_a,qw_a,qi_a,p_ltabdminwgg)
            ELSE
#ifdef _OPENACC
              CALL finish(routine, 'dmin_wetgrowth_fun not available on GPU')
#endif
              d_trenn = dmin_wetgrowth_fun(p_a,T_a,qw_a,qi_a)
            END IF

            IF (d_trenn > 0.0_wp .AND. d_trenn < 10.0_wp * D_g) THEN

               xmin = exp(log(d_trenn/graupel%a_geo)*(1.0_wp/graupel%b_geo))

               lam  = exp(log(graupel_g2/(graupel_g1*x_g))*(graupel%mu))
               xmin = exp(log(xmin)*graupel%mu)
               n_0  = graupel%mu * n_g * exp(log(lam)*graupel_nm1) / graupel_g1

               conv_n = n_0 / (graupel%mu * exp(log(lam)*graupel_nm1)) * incgfct_upper_lookup(lam*xmin,graupel_ltable1)
               conv_q = n_0 / (graupel%mu * exp(log(lam)*graupel_nm2)) * incgfct_upper_lookup(lam*xmin,graupel_ltable2)

               conv_n = MIN(conv_n,n_g)
               conv_q = MIN(conv_q,q_g)

               graupel%q(i,k) = q_g - conv_q
               graupel%n(i,k) = n_g - conv_n

               hail%q(i,k) = hail%q(i,k) + conv_q
               hail%n(i,k) = hail%n(i,k) + conv_n

            END IF
          ENDIF
       ENDDO
    ENDDO
   !$ACC END PARALLEL
   !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer
   !$ACC END DATA

  END SUBROUTINE graupel_hail_conv_wet_gamlook

  SUBROUTINE ice_riming(ik_slice, dt, icr_coeffs, irr_coeffs, &
       &                atmo, ice_in, cloud, rain_in, graupel_in, dep_rate_ice)
    !*******************************************************************************
    !  Riming of ice with cloud droplet and rain drops. First the process rates    *
    !  are calculated in                                                           *
    !      snow_cloud_riming ()                                                    *
    !      snow_rain_riming ()                                                     *
    !  using those rates and the previously calculated and stored deposition       *
    !  rate the conversion of snow to graupel and rain is done.                    *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    ! parameters
    TYPE(collection_coeffs), INTENT(in)  :: icr_coeffs
    TYPE(rain_riming_coeffs), INTENT(in) :: irr_coeffs
    ! progn. variables
    TYPE(atmosphere), INTENT(inout)       :: atmo
    CLASS(particle), INTENT(inout)        :: cloud
    CLASS(particle_frozen), INTENT(inout), TARGET :: graupel_in
    CLASS(particle_frozen), POINTER :: graupel ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle), INTENT(inout), TARGET :: rain_in
    CLASS(particle), POINTER :: rain ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle_frozen), INTENT(inout), TARGET :: ice_in
    CLASS(particle_frozen), POINTER:: ice ! ACCWA (nvhpc 22.7, IPSF, see above)

    REAL(wp), INTENT (IN), DIMENSION(:,:) :: dep_rate_ice
    REAL(wp), DIMENSION(size(dep_rate_ice,1),size(dep_rate_ice,2)) ::       &
         &               rime_rate_qc, rime_rate_nc,                        &
         &               rime_rate_qi, rime_rate_qr, rime_rate_nr

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_i,n_i,x_i,d_i
    REAL(wp)            :: T_a,x_r,D_r
    REAL(wp)            :: x_coll, D_coll, rho_coll, D_i_coll, D_g_coll, rho_i_coll, rho_g_coll, rho_limit
    REAL(wp)            :: rime_n,rime_q,rime_qr,rime_qi
    REAL(wp)            :: conv_n,conv_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2,const5
    LOGICAL             :: grconvflag

    REAL(wp), PARAMETER :: &
         const3 = 1.0_wp/(T_mult_opt - T_mult_min), &
         const4 = 1.0_wp/(T_mult_opt - T_mult_max)

    IF (isdebug) CALL message(routine, "ice riming")

    ice => ice_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    rain => rain_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    graupel => graupel_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)


    !$ACC DATA CREATE(rime_rate_qc, rime_rate_nc, rime_rate_qi, rime_rate_qr, rime_rate_nr)
    CALL riming_cloud_core(ik_slice, ice, cloud, icr_coeffs, dt, &
         &                 rime_rate_qc, rime_rate_nc)
    CALL riming_rain_core(ik_slice, ice, rain, irr_coeffs, dt, &
         &                rime_rate_qi, rime_rate_qr, rime_rate_nr)

    !
    ! Complete ice-cloud and ice-rain riming

!!$ This changes the results: !!!    const5 = rho_w/rho_ice * cfg_params%alpha_spacefilling
    const5 = cfg_params%alpha_spacefilling * rho_w/rho_ice

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(q_i, n_i, x_i, d_i, T_a, D_r, x_r, rime_n, rime_q) &
    !$ACC   PRIVATE(rime_qr, rime_qi, conv_n, conv_q, mult_n, mult_q, mult_1, mult_2) &
    !$ACC   PRIVATE(x_coll, D_coll, rho_coll, D_i_coll, D_g_coll, rho_i_coll, rho_g_coll, rho_limit, grconvflag)
    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)

          IF (dep_rate_ice(i,k) > 0.0_wp &
               & .and. dep_rate_ice(i,k) .ge. rime_rate_qc(i,k)+rime_rate_qr(i,k)) THEN

            ! 1) Depositional growth is stronger than riming growth, therefore ice stays ice

            !.. ice_cloud_riming

            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(cloud%q(i,k),rime_q)
              rime_n = MIN(cloud%n(i,k),rime_n)

              ice%q(i,k)   = ice%q(i,k)  + rime_q
              cloud%q(i,k) = cloud%q(i,k) - rime_q
              cloud%n(i,k) = cloud%n(i,k) - rime_n

              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                ice%n(i,k) = ice%n(i,k)  + mult_n
              ENDIF
            END IF

            !.. ice_rain_riming
            IF (rime_rate_qr(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qr(i,k)
              rime_n = rime_rate_nr(i,k)
              rime_q = MIN(rain%q(i,k),rime_q)
              rime_n = MIN(rain%n(i,k),rime_n)

              ice%q(i,k)  = ice%q(i,k)  + rime_q
              rain%q(i,k) = rain%q(i,k) - rime_q
              rain%n(i,k) = rain%n(i,k) - rime_n

              !..ice multiplication
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                ice%n(i,k) = ice%n(i,k)  + mult_n
              ENDIF

            END IF

          ELSE
            
            !.. 2) Depositional growth negative or smaller than riming growth, therefore ice is
            !      allowed to convert to graupel and / or hail

            !.. ice_cloud_riming
           
            n_i = ice%n(i,k)
            q_i = ice%q(i,k)
            x_i = particle_meanmass(ice, q_i,n_i)
            d_i = particle_diameter(ice, x_i)
              
            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(cloud%q(i,k),rime_q)
              rime_n = MIN(cloud%n(i,k),rime_n)

              ice%q(i,k)   = ice%q(i,k)   + rime_q
              cloud%q(i,k) = cloud%q(i,k) - rime_q
              cloud%n(i,k) = cloud%n(i,k) - rime_n

              ! ice multiplication
              mult_q = 0.0_wp
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                ice%n(i,k) = ice%n(i,k)  + mult_n
              ENDIF

              ! conversion ice -> graupel (depends on alpha_spacefilling)
              IF (D_i > D_conv_ig .AND. T_a < cfg_params%Tmax_gr_rime) THEN
                 q_i = ice%q(i,k)
                 conv_q = (rime_q - mult_q) / ( const5 * (pi6*rho_ice*d_i**3/x_i - 1.0_wp) )
                 conv_q = MIN(q_i,conv_q)
                 x_i    = particle_meanmass(ice, q_i,n_i)
                 conv_n = conv_q / MAX(x_i,x_conv)
                 conv_n = MIN(ice%n(i,k),conv_n)
              ELSE
                 conv_q = 0.0_wp
                 conv_n = 0.0_wp
              ENDIF

              ice%q(i,k)     = ice%q(i,k)     - conv_q
              graupel%q(i,k) = graupel%q(i,k) + conv_q
              ice%n(i,k)     = ice%n(i,k)     - conv_n
              graupel%n(i,k) = graupel%n(i,k) + conv_n
            END IF

            !.. ice_rain_riming
            IF (rime_rate_qi(i,k) > 0.0_wp) THEN

              x_r = particle_meanmass(rain, rain%q(i,k),rain%n(i,k))
              D_r = particle_diameter(rain, x_r)
              
              rime_qi = rime_rate_qi(i,k)
              rime_qr = rime_rate_qr(i,k)
              rime_n  = rime_rate_nr(i,k)
              rime_n  = MIN(MIN(rain%n(i,k),ice%n(i,k)),rime_n)
              rime_qr = MIN(rain%q(i,k),rime_qr)
              rime_qi = MIN(ice%q(i,k),rime_qi)

              ice%n(i,k)  = ice%n(i,k)  - rime_n
              rain%n(i,k) = rain%n(i,k) - rime_n
              ice%q(i,k)  = ice%q(i,k)  - rime_qi
              rain%q(i,k) = rain%q(i,k) - rime_qr

              ! ice multiplication
              mult_q = 0.0
              mult_n = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)
              ENDIF

              IF (T_a >= T_3) THEN
                 ! shedding of rain at warm temperatures
                 ! i.e. undo time integration
                 ice%n(i,k)  = ice%n(i,k)  + rime_n
                 rain%n(i,k) = rain%n(i,k) + rime_qr / x_r
                 ice%q(i,k)  = ice%q(i,k)  + rime_qi
                 rain%q(i,k) = rain%q(i,k) + rime_qr
              ELSE
                 ! new ice particles from multiplication
                 ice%n(i,k) = ice%n(i,k) + mult_n
                 ice%q(i,k) = ice%q(i,k) + mult_q
                 
                 IF (cfg_params%llim_gr_prod_rain_riming) THEN
                   ! riming to graupel, if bulk density of the collided mean mass
                   ! particles is nearer to equivalent graupel- than to ice bulk density of the collided particle:
                   D_coll = MAX(D_i,D_r) + cfg_params%wgt_D_coll_limgrprod * MIN(D_i,D_r) ! collided diameter, weighted sum of both partners
                   x_coll = x_i + x_r  ! mass of collided particles (the representative mean mass particles)
                   rho_coll = x_coll / (pi6*D_coll**3)
                   D_i_coll = particle_diameter(ice, x_coll)     ! mean mass diameter of ice particle having mass x_coll
                   rho_i_coll = x_coll / (pi6*D_i_coll**3)
                   D_g_coll = particle_diameter(graupel, x_coll) ! mean mass diameter of graupel particle having mass x_coll
                   rho_g_coll = x_coll / (pi6*D_g_coll**3)
                   rho_limit = (1.0_wp-cfg_params%wgt_rho_coll_limgrprod) * rho_i_coll + &
                        cfg_params%wgt_rho_coll_limgrprod  * rho_g_coll

                   IF (rho_i_coll <  rho_g_coll) THEN
                     grconvflag = (rho_coll > rho_limit)
                   ELSE
                     grconvflag = (rho_coll < rho_limit)
                   END IF

                   ! if the average collided particle is larger than the upper mass limit for ice,
                   ! also convert to graupel:
                   IF (x_coll > ice%x_max) grconvflag = .TRUE.
                   
                 ELSE
                   grconvflag = .TRUE.
                 END IF

                 IF (T_a < cfg_params%Tmax_gr_rime .AND. grconvflag )THEN
                   graupel%n(i,k) = graupel%n(i,k) + rime_n
                   graupel%q(i,k) = graupel%q(i,k) + rime_qi + rime_qr - mult_q
                 ELSE
                   ! Ice + frozen liquid stays ice:
                   ice%n(i,k) = ice%n(i,k) + rime_n
                   ice%q(i,k) = ice%q(i,k) + rime_qi + rime_qr - mult_q
                 END IF
              END IF
            END IF
          END IF
       END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE ice_riming

  SUBROUTINE snow_riming(ik_slice, dt, scr_coeffs, srr_coeffs, &
       &                 atmo, snow_in, cloud, rain_in, ice, graupel_in, dep_rate_snow)
    !*******************************************************************************
    !  Riming of snow with cloud droplet and rain drops. First the process rates   *
    !  are calculated in                                                           *
    !      snow_cloud_riming ()                                                    *
    !      snow_rain_riming ()                                                     *
    !  using those rates and the previously calculated and stored deposition       *
    !  rate the conversion of snow to graupel and rain is done.                    *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    ! parameters
    TYPE(collection_coeffs), INTENT(in) :: scr_coeffs  ! snow cloud riming
    TYPE(rain_riming_coeffs),INTENT(in) :: srr_coeffs  ! snow rain riming
    ! 2mom variables
    TYPE(atmosphere), INTENT(inout)       :: atmo
    CLASS(particle), INTENT(inout)        :: cloud
    CLASS(particle_frozen), INTENT(inout) :: ice
    CLASS(particle_frozen), INTENT(inout), TARGET :: graupel_in
    CLASS(particle_frozen), POINTER :: graupel ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle), INTENT(inout), TARGET :: rain_in
    CLASS(particle), POINTER :: rain ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle_frozen), INTENT(inout), TARGET :: snow_in
    CLASS(particle_frozen), POINTER:: snow ! ACCWA (nvhpc 22.7, IPSF, see above)

    REAL(wp), INTENT(IN), DIMENSION(:,:)  :: dep_rate_snow
    REAL(wp), DIMENSION(size(dep_rate_snow,1),size(dep_rate_snow,2)) ::       &
         & rime_rate_qc, rime_rate_nc,                                        &
         & rime_rate_qs, rime_rate_qr, rime_rate_nr

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_s,n_s,x_s,d_s, x_r, D_r
    REAL(wp)            :: x_coll, D_coll, rho_coll, D_s_coll, D_g_coll, rho_s_coll, rho_g_coll, rho_limit
    REAL(wp)            :: rime_n,rime_q,rime_qr,rime_qs
    REAL(wp)            :: conv_n,conv_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2,const5
    LOGICAL             :: grconvflag

    REAL(wp), PARAMETER :: &
         const3 = 1.0/(T_mult_opt - T_mult_min), &
         const4 = 1.0/(T_mult_opt - T_mult_max)

    IF (isdebug) CALL message(routine, "snow_riming")

    graupel => graupel_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    snow => snow_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    rain => rain_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)


    !$ACC DATA CREATE(rime_rate_qc, rime_rate_nc, rime_rate_qs, rime_rate_qr, rime_rate_nr)
    CALL riming_cloud_core(ik_slice, snow, cloud, scr_coeffs, dt, &
         &                 rime_rate_qc, rime_rate_nc)
    CALL riming_rain_core(ik_slice, snow, rain, srr_coeffs, dt, &
         &                rime_rate_qs, rime_rate_qr, rime_rate_nr)

!!$ This changes the results: !!!    const5 = rho_w/rho_ice * cfg_params%alpha_spacefilling
    const5 = cfg_params%alpha_spacefilling * rho_w/rho_ice 

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(T_a, q_s, n_s, x_s, d_s, x_r, D_r, rime_n, rime_q) &
    !$ACC   PRIVATE(rime_qr, rime_qs, conv_n, conv_q, mult_n, mult_q, mult_1, mult_2) &
    !$ACC   PRIVATE(x_coll, D_coll, rho_coll, D_s_coll, D_g_coll, rho_s_coll, rho_g_coll, rho_limit, grconvflag)
    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)

          IF (dep_rate_snow(i,k) > 0.0_wp &
               & .AND. dep_rate_snow(i,k) .ge. rime_rate_qc(i,k) + rime_rate_qr(i,k)) THEN

            ! 1) Depositional growth is stronger than riming growth, therefore snow stays snow:

            !.. time integration of snow_cloud_riming

            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(cloud%q(i,k),rime_q)
              rime_n = MIN(cloud%n(i,k),rime_n)

              snow%q(i,k)  = snow%q(i,k)  + rime_q
              cloud%q(i,k) = cloud%q(i,k) - rime_q
              cloud%n(i,k) = cloud%n(i,k) - rime_n

              ! ice multiplication
              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)
                mult_n = mult_q / ice%x_min

                ice%n(i,k)  = ice%n(i,k)  + mult_n
                ice%q(i,k)  = ice%q(i,k)  + mult_q
                snow%q(i,k) = snow%q(i,k) - mult_q
              ENDIF

            END IF

            !.. time integration snow_rain_riming

            IF (rime_rate_qr(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qr(i,k)
              rime_n = rime_rate_nr(i,k)
              rime_q = MIN(rain%q(i,k),rime_q)
              rime_n = MIN(rain%n(i,k),rime_n)

              snow%q(i,k) = snow%q(i,k) + rime_q
              rain%q(i,k) = rain%q(i,k) - rime_q
              rain%n(i,k) = rain%n(i,k) - rime_n

              ! ice multiplication
              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)
                mult_n = mult_q / ice%x_min

                ice%n(i,k)  = ice%n(i,k)  + mult_n
                ice%q(i,k)  = ice%q(i,k)  + mult_q
                snow%q(i,k) = snow%q(i,k) - mult_q
              ENDIF
            END IF

          ELSE

            !.. 2) Depositional growth is negative or smaller than riming growth, therefore snow is
            !      allowed to convert to graupel and / or hail:

            !.. time integration of snow_cloud_riming

            n_s = snow%n(i,k)
            q_s = snow%q(i,k)
            x_s = particle_meanmass(snow, q_s,n_s)
            d_s = particle_diameter(snow, x_s)

            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(cloud%q(i,k),rime_q)
              rime_n = MIN(cloud%n(i,k),rime_n)

              snow%q(i,k)  = snow%q(i,k)  + rime_q
              cloud%q(i,k) = cloud%q(i,k) - rime_q
              cloud%n(i,k) = cloud%n(i,k) - rime_n

              ! ice multiplication
              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)
                mult_n = mult_q / ice%x_min

                ice%n(i,k)  = ice%n(i,k)  + mult_n
                ice%q(i,k)  = ice%q(i,k)  + mult_q
                snow%q(i,k) = snow%q(i,k) - mult_q
              ENDIF

              !.. conversion of snow to graupel, depends on alpha_spacefilling

              IF (D_s > D_conv_sg .AND. T_a < cfg_params%Tmax_gr_rime) THEN
                 q_s = snow%q(i,k)  
                 conv_q = (rime_q - mult_q) / ( const5 * (pi6*rho_ice*d_s**3/x_s - 1.0_wp) )
                 conv_q = MIN(q_s,conv_q)
                 x_s    = particle_meanmass(snow, q_s,n_s)
                 conv_n = conv_q / MAX(x_s,x_conv)
                 conv_n = MIN(snow%n(i,k),conv_n)
              ELSE
                 conv_q = 0.0_wp
                 conv_n = 0.0_wp
              ENDIF

              snow%q(i,k)    = snow%q(i,k)    - conv_q
              graupel%q(i,k) = graupel%q(i,k) + conv_q

              snow%n(i,k)    = snow%n(i,k)    - conv_n
              graupel%n(i,k) = graupel%n(i,k) + conv_n

            END IF

            !.. time integration of snow_rain_riming

            IF (rime_rate_qs(i,k) > 0.0_wp) THEN

              x_r = particle_meanmass(rain, rain%q(i,k),rain%n(i,k))
              D_r = particle_diameter(rain, x_r)

              rime_qs = rime_rate_qs(i,k)
              rime_qr = rime_rate_qr(i,k)
              rime_n  = rime_rate_nr(i,k)
              rime_qr = MIN(rain%q(i,k),rime_qr)
              rime_qs = MIN(snow%q(i,k),rime_qs)
              rime_n  = MIN(rain%n(i,k),rime_n)
              rime_n  = MIN(snow%n(i,k),rime_n)

              snow%n(i,k) = snow%n(i,k) - rime_n
              rain%n(i,k) = rain%n(i,k) - rime_n
              snow%q(i,k) = snow%q(i,k) - rime_qs
              rain%q(i,k) = rain%q(i,k) - rime_qr

              ! ice multiplication
              mult_q = 0.0_wp
              mult_n = 0.0_wp
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)
                mult_n = mult_q / ice%x_min
              ENDIF

              IF (T_a >= T_3) THEN
                 ! shedding of rain at warm temperatures
                 ! i.e. undo time integration
                 snow%n(i,k) = snow%n(i,k) + rime_n
                 rain%n(i,k) = rain%n(i,k) + rime_qr / x_r
                 snow%q(i,k) = snow%q(i,k) + rime_qs
                 rain%q(i,k) = rain%q(i,k) + rime_qr
              ELSE
                 ! new ice particles from multiplication
                 ice%n(i,k)  = ice%n(i,k)  + mult_n
                 ice%q(i,k)  = ice%q(i,k)  + mult_q

                 IF (cfg_params%llim_gr_prod_rain_riming) THEN
                   ! riming to graupel, if bulk density of the collided mean mass
                   ! particles is nearer to equivalent graupel- than to ice bulk density of the collided particle:
                   D_coll = MAX(D_s,D_r) + cfg_params%wgt_D_coll_limgrprod * MIN(D_s,D_r) ! collided diameter, weighted sum of both partners
                   x_coll = x_s + x_r  ! mass of collided particles (the representative mean mass particles)
                   rho_coll = x_coll / (pi6*D_coll**3)
                   D_s_coll = particle_diameter(snow, x_coll)     ! mean mass diameter of ice particle having mass x_coll
                   rho_s_coll = x_coll / (pi6*D_s_coll**3)
                   D_g_coll = particle_diameter(graupel, x_coll) ! mean mass diameter of graupel particle having mass x_coll
                   rho_g_coll = x_coll / (pi6*D_g_coll**3)
                   rho_limit = (1.0_wp-cfg_params%wgt_rho_coll_limgrprod) * rho_s_coll + &
                                       cfg_params%wgt_rho_coll_limgrprod  * rho_g_coll

                   IF (rho_s_coll <  rho_g_coll) THEN
                     grconvflag = (rho_coll > rho_limit)
                   ELSE
                     grconvflag = (rho_coll < rho_limit)
                   END IF

                 ELSE
                   grconvflag = .TRUE.
                 END IF

                 IF (T_a < cfg_params%Tmax_gr_rime .AND. grconvflag )THEN
                   graupel%n(i,k) = graupel%n(i,k) + rime_n
                   graupel%q(i,k) = graupel%q(i,k) + rime_qr + rime_qs - mult_q
                 ELSE
                   ! Snow + frozen liquid stays snow:
                   snow%n(i,k) = snow%n(i,k) + rime_n
                   snow%q(i,k) = snow%q(i,k) + rime_qr + rime_qs - mult_q
                 END IF
              END IF

           END IF
        END IF

      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

 END SUBROUTINE snow_riming

 SUBROUTINE riming_cloud_core(ik_slice, ptype_in, cloud_in, coeffs, dt, &
      &                     rime_rate_qb, rime_rate_nb)
   !*******************************************************************************
   !  Riming rate of ice or snow collecting cloud droplets                        *
   !  This is a process of the form a+b->c, but b is here always cloud and        *
   !  therefore some parameters are hardcoded for cloud water                     *
   !*******************************************************************************
   ! start and end indices for 2D slices
   ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
   INTEGER, INTENT(in) :: ik_slice(4)
   CLASS(particle_frozen), INTENT(in), TARGET:: ptype_in
   CLASS(particle_frozen), POINTER :: ptype ! ACCWA (nvhpc 22.7, IPSF, see above)
   CLASS(particle), INTENT(in), TARGET :: cloud_in
   CLASS(particle), POINTER :: cloud ! ACCWA (nvhpc 22.7, IPSF, see above)

   REAL(wp), INTENT(in)                :: dt 
   TYPE(collection_coeffs), INTENT(in) :: coeffs
   REAL(wp), INTENT(out)               :: rime_rate_qb(:, :), rime_rate_nb(:, :)
   
   ! start and end indices for 2D slices
   INTEGER :: istart, iend, kstart, kend
   INTEGER             :: i,k
   REAL(wp)            :: q_p,n_p,x_p,d_p,v_p
   REAL(wp)            :: q_c,n_c,x_c,d_c,v_c,e_coll,const1
   REAL(wp)            :: rime_n,rime_q
   REAL(wp), PARAMETER :: &
        &  const0 = 1.0/(D_coll_c - D_crit_c)

   IF (isdebug) CALL message(routine, "riming_cloud_core")
   
   ptype => ptype_in ! ACCWA (nvhpc 22.7, IPSF, see above)
   cloud => cloud_in ! ACCWA (nvhpc 22.7, IPSF, see above)

   istart = ik_slice(1)
   iend   = ik_slice(2)
   kstart = ik_slice(3)
   kend   = ik_slice(4)

   const1 = const0 * ptype%ecoll_c

   !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
   !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(q_p, n_p, x_p, d_p, v_p) &
   !$ACC   PRIVATE(q_c, n_c, x_c, d_c, v_c, e_coll, rime_n, rime_q)
   DO k = kstart,kend
     DO i = istart,iend
       
       n_p = ptype%n(i,k)
       q_p = ptype%q(i,k)
       n_c = cloud%n(i,k)
       q_c = cloud%q(i,k)
       
       x_p = particle_meanmass(ptype, q_p, n_p)
       d_p = particle_diameter(ptype, x_p)
       x_c = particle_meanmass(cloud, q_c, n_c)
       d_c = particle_diameter(cloud, x_c)
       
       IF (q_c > q_crit_c .AND. q_p > ptype%q_crit_c &
            &             .AND. d_p > ptype%D_crit_c .AND. D_c > D_crit_c) THEN
          
         v_c = particle_velocity(cloud,x_c) * cloud%rho_v(i,k)
         v_p = particle_velocity(ptype,x_p) * ptype%rho_v(i,k)
         
         e_coll = MIN(ptype%ecoll_c, MAX(const1*(d_c - D_crit_c), ecoll_min))
         
         rime_n = pi4 * e_coll * n_p * n_c * dt &
              & *     (  coeffs%delta_n_aa * d_p**2 &
              &        + coeffs%delta_n_ab * d_p*d_c &
              &        + coeffs%delta_n_bb * d_c**2) &
              & * SQRT(  coeffs%theta_n_aa * v_p**2 &
              &        - coeffs%theta_n_ab * v_p*v_c &
              &        + coeffs%theta_n_bb * v_c**2 &
              &        + ptype%s_vel**2)
         
         rime_q = pi4 * e_coll * n_p * q_c * dt &
              & *     (  coeffs%delta_q_aa * d_p**2 &
              &        + coeffs%delta_q_ab * d_p*D_c &
              &        + coeffs%delta_q_bb * d_c**2) &
              & * SQRT(  coeffs%theta_q_aa * v_p**2 &
              &        - coeffs%theta_q_ab * v_p*v_c &
              &        + coeffs%theta_q_bb * v_c**2 &
              &        + ptype%s_vel**2)
         
         rime_rate_qb(i,k) = rime_q
         rime_rate_nb(i,k) = rime_n
       ELSE
         rime_rate_qb(i,k) = 0.0_wp
         rime_rate_nb(i,k) = 0.0_wp
       ENDIF
     ENDDO
   ENDDO
   !$ACC END PARALLEL
   !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

 END SUBROUTINE riming_cloud_core

  SUBROUTINE riming_rain_core(ik_slice, ptype_in, rain_in, coeffs, dt, &
       &                      rime_rate_qa, rime_rate_qb, rime_rate_nb)
    !*******************************************************************************
    !  Riming rate of ice collecting rain drop, or rain collecting ice             *
    !  and riming of snow with raindrops.                                          *
    !  This is a process of the form a+b->c, but b is here always rain and         *
    !  therefore some parameters are hardcoded for rain                            *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER,  INTENT(in) :: ik_slice(4)
    REAL(wp), INTENT(in) :: dt
    
    ! 2mom variables
    CLASS(particle_frozen), INTENT(in), TARGET :: ptype_in
    CLASS(particle_frozen), POINTER :: ptype ! ACCWA (nvhpc 22.7, IPSF, see above)
    CLASS(particle), INTENT(in), TARGET :: rain_in
    CLASS(particle), POINTER :: rain ! ACCWA (nvhpc 22.7, IPSF, see above)

    TYPE(rain_riming_coeffs), INTENT(in) :: coeffs
    REAL(wp), INTENT(out)                :: rime_rate_qa(:,:), rime_rate_qb(:,:), &
         &                                  rime_rate_nb(:,:)

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_a,n_a,x_a,d_a,v_a
    REAL(wp)            :: q_r,n_r,x_r,d_r,v_r
    REAL(wp)            :: rime_n,rime_qi,rime_qr

    IF (isdebug) CALL message(routine, "ice_rain_riming")

    ptype => ptype_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    rain => rain_in ! ACCWA (nvhpc 22.7, IPSF, see above)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) FIRSTPRIVATE(kstart, kend, istart, iend, dt)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(q_a, n_a, x_a, d_a, v_a) &
    !$ACC   PRIVATE(q_r, n_r, x_r, d_r, v_r, rime_n, rime_qi, rime_qr)
    DO k = kstart,kend
      DO i = istart,iend

        q_r = rain%q(i,k)
        n_r = rain%n(i,k)
        q_a = ptype%q(i,k)
        n_a = ptype%n(i,k)

        x_a = particle_meanmass(ptype, q_a,n_a)
        d_a = particle_diameter(ptype, x_a)

        IF (q_r > q_crit .AND. q_a > q_crit_r .AND. d_a > D_crit_r) THEN

          x_r = particle_meanmass(rain,q_r,n_r)
          d_r = particle_diameter(rain,x_r)
          v_r = particle_velocity(rain,x_r) * rain%rho_v(i,k)
          v_a = particle_velocity(ptype,x_a) * ptype%rho_v(i,k)

          rime_n  = pi4 * n_a * n_r * dt &
               &  *     (  coeffs%delta_n_aa * d_a * d_a &
               &         + coeffs%delta_n_ab * d_a * d_r &
               &         + coeffs%delta_n_bb * d_r * d_r) &
               &  * SQRT(coeffs%theta_n_aa * v_a * v_a &
               &         - coeffs%theta_n_ab * v_a * v_r &
               &         + coeffs%theta_n_bb * v_r * v_r &
               &         + ptype%s_vel**2)

          rime_qr = pi4 * n_a * q_r * dt &
               &  *     (  coeffs%delta_n_aa * d_a * d_a &
               &         + coeffs%delta_q_ab * d_a * d_r &
               &         + coeffs%delta_q_bb * d_r * d_r) &
               &  * SQRT(  coeffs%theta_n_aa * v_a * v_a &
               &         - coeffs%theta_q_ab * v_a * v_r &
               &         + coeffs%theta_q_bb * v_r * v_r &
               &         + ptype%s_vel**2)

          rime_qi = pi4 * n_r * q_a * dt &
               &  *     (  coeffs%delta_q_aa * d_a * d_a &
               &         + coeffs%delta_q_ba * d_a * d_r &
               &         + coeffs%delta_n_bb * d_r * d_r) &
               &  * SQRT(  coeffs%theta_q_aa * v_a * v_a &
               &         - coeffs%theta_q_ba * v_a * v_r &
               &         + coeffs%theta_n_bb * v_r * v_r &
               &         + ptype%s_vel**2)

          rime_rate_nb(i,k) = rime_n
          rime_rate_qa(i,k) = rime_qi
          rime_rate_qb(i,k) = rime_qr
        ELSE
          rime_rate_nb(i,k) = 0.0_wp
          rime_rate_qa(i,k) = 0.0_wp
          rime_rate_qb(i,k) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! ACCWA (nvhpc 22.7): wait is required for intermediate pointer

  END SUBROUTINE riming_rain_core

 SUBROUTINE ccn_activation_sk(ik_slice, ccn_coeffs, atmo, cloud, n_cn)
   !*******************************************************************************
   !       Calculation of ccn activation                                          *
   !       using the look-up tables by Segal and Khain 2006 (JGR, vol.11)         *
   !       (implemented by Heike Noppel) NOT USED ANYMORE!                        *
   !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)

    ! parameters
    TYPE(aerosol_ccn),INTENT(in)    :: ccn_coeffs

    ! 2mom variables
    TYPE(atmosphere), INTENT(inout)     :: atmo
    CLASS(particle), INTENT(inout)      :: cloud
    REAL(wp), DIMENSION(:,:), OPTIONAL  :: n_cn

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend

    ! local variables
    INTEGER             :: i,k,nuc_typ
    REAL(wp)            :: n_c,q_c
    REAL(wp)            :: nuc_n, nuc_q, zf
    REAL(wp), PARAMETER :: eps = 1e-10_wp
    
    ! for activation tables
    INTEGER, PARAMETER :: n_ncn=8, n_r2=3, n_lsigs=5, n_wcb=4
    INTEGER            :: i_lsigs, i_R2
    REAL(wp)           :: Ncn, wcb, &
         tab_Ndrop(n_wcb,n_ncn), & ! number of cloud droplets in look_up-Table
         tab_Ndrop_i(n_wcb)        ! number of cloud droplets in look_up-Table,
                                   ! interpolated with respect to Ncn

    REAL(wp), PARAMETER :: &
          ! look-up-table for Ncn
         tab_Ncn(n_ncn) = (/   50.d06,  100.d06,  200.d06,  400.d06, &
         &                    800.d06, 1600.d06, 3200.d06, 6400.d06 /), &
         ! look-up_tbale for R2, in 10^(-6) m
         tab_R2(n_r2) = (/0.02d0, 0.03d0, 0.04d0/), &
         ! look-up-table for log(sigma_s)
         tab_lsigs(n_lsigs)  = (/0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0/), &
         ! look-up-table for w at cloud base
         tab_wcb(n_wcb) = (/0.5d0, 1.0d0, 2.5d0, 5.0d0/)

    nuc_typ = nuc_c_typ

    ! ATTENTION: At the moment only the values given above can be chosen for R2,
    ! and lsigs (see below).
    ! Only wcb and N_cn0 can be interpolated. For the others this possibility is still missing.
    ! For high N_cn0 and small values of wcb a kind of "saturation" for cloud droplet
    ! nucleation was assumed, because the original look-up tables don't give these values.
    ! This assumption might be wrong!!! (comment by Heike Noppel)

    IF(isdebug) THEN
       WRITE(txt,*) "cloud_activation_SK: nuc_typ = ",nuc_typ
       CALL message(routine, TRIM(txt))
    ENDIF

    ! this could be done in an IF(firstcall.eq.false) block and i-value stored in ccn_coeffs
    i_lsigs = 0
    i_R2   = 0
    DO k=1,n_lsigs
      IF (ccn_coeffs%lsigs==tab_lsigs(k)) THEN
         i_lsigs=k
         EXIT
      ENDIF
    END DO
    DO k=1,n_r2
      IF (ccn_coeffs%R2==tab_R2(k)) THEN
         i_R2=k
         EXIT
      ENDIF
    END DO
    IF (i_lsigs==0) THEN
       CALL finish(TRIM(routine),'Error in two_moment_mcrph: Invalid value for LSIGS in ccn_activation_sk')
    END IF
    IF (i_R2==0) THEN
       CALL finish(TRIM(routine),'Error in two_moment_mcrph: Invalid value for R2 in ccn_activation_sk')
    END IF

    CALL lookuptable(tab_Ndrop,i_lsigs,i_R2)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          nuc_q = 0.0d0
          nuc_n = 0.d0
          n_c   = cloud%n(i,k)
          q_c   = cloud%q(i,k)
          wcb   = atmo%w(i,k+1)  ! This is the velocity at cloud base, therefore k+1

          if (q_c > eps .and. wcb > 0.0_wp) then

             ! hard upper limit for number cnc that
             ! eliminates also unrealistic high value
             ! that would come from the dynamical core

             cloud%n(i,k) = MIN(cloud%n(i,k),ccn_coeffs%Ncn0)

             IF (PRESENT(n_cn)) THEN
               Ncn = n_cn(i,k) ! number of CN from prognostic variable
             ELSE
               zf = 0.5_wp*(atmo%zh(i,k)+atmo%zh(i,k+1))             
               IF(zf > ccn_coeffs%z0) THEN
                 Ncn = ccn_coeffs%Ncn0 * EXP((ccn_coeffs%z0 - zf)/ccn_coeffs%z1e)
               ELSE
                 Ncn = ccn_coeffs%Ncn0
               END IF
             END IF
             
             ! min value for vertical velocity (instead of Nmin of older code)
             wcb = MAX(wcb,ccn_coeffs%wcb_min)

             ! Interpolation of the look-up tables with respect to Ncn

             ! If Ncn is outside the range of the lookup table values, resulting
             ! NCCN are clipped to the margin values. For the case of these margin values
             ! being larger than Ncn, limit NCCN by Ncn:
             tab_Ndrop_i = MIN(ip_ndrop_ncn(tab_Ndrop,tab_Ncn,Ncn), Ncn)

             ! interpol. with respect to wcb = ip_ndrop_wcb
             nuc_n = MAX(ccn_coeffs%etas * ip_ndrop_wcb(tab_Ndrop_i,tab_wcb,wcb),ccn_coeffs%Nmin) - n_c

             nuc_n = MAX(nuc_n,0.0d0)

             nuc_q = MIN(nuc_n * cloud%x_min, atmo%qv(i,k))
             nuc_n = nuc_q / cloud%x_min

             cloud%n(i,k) = cloud%n(i,k) + nuc_n
             cloud%q(i,k) = cloud%q(i,k) + nuc_q
             atmo%qv(i,k) = atmo%qv(i,k) - nuc_q

             if (present(n_cn)) then
               n_cn(i,k)    = n_cn(i,k)    - MIN(Ncn,nuc_n)
             end if

          END IF

       END DO
    END DO

  CONTAINS

    SUBROUTINE lookuptable(tab_ndrop,i_lsigs,i_R2)
      REAL(wp), INTENT(out) :: tab_ndrop(n_wcb,n_ncn)
      INTEGER, INTENT(in) :: i_lsigs, i_R2
      INTEGER :: i
      REAL(wp), PARAMETER :: ndrop(n_ncn,5,n_wcb,3) &
           ! look up tables
           ! Ncn              50       100       200       400       800       1600      3200      6400
           ! table4a (R2=0.02mum, wcb=0.5m/s) (for Ncn=3200  and Ncn=6400 "extrapolated")
           = RESHAPE((/ &
           !ndrop1_11 =
           42.2d06,  70.2d06, 112.2d06, 173.1d06, 263.7d06, 397.5d06, 397.5d06, 397.5d06, &
           !ndrop1_12 =
           35.5d06,  60.1d06, 100.0d06, 163.9d06, 264.5d06, 418.4d06, 418.4d06, 418.4d06, &
           !ndrop1_13
           32.6d06,  56.3d06,  96.7d06, 163.9d06, 272.0d06, 438.5d06, 438.5d06, 438.5d06, &
           !ndrop1_14
           30.9d06,  54.4d06,  94.6d06, 162.4d06, 271.9d06, 433.5d06, 433.5d06, 433.5d06, &
           !ndrop1_15
           29.4d06,  51.9d06,  89.9d06, 150.6d06, 236.5d06, 364.4d06, 364.4d06, 364.4d06, &
           ! table4b (R2=0.02mum, wcb=1.0m/s) (for Ncn=50 "interpolted" and Ncn=6400 extrapolated)
           !ndrop1_21
           45.3d06,  91.5d06, 158.7d06, 264.4d06, 423.1d06, 672.5d06, 397.5d06, 397.5d06, &
           !ndrop1_22
           38.5d06,  77.1d06, 133.0d06, 224.9d06, 376.5d06, 615.7d06, 418.4d06, 418.4d06, &
           !ndrop1_23
           35.0d06,  70.0d06, 122.5d06, 212.0d06, 362.1d06, 605.3d06, 438.5d06, 438.5d06, &
           !ndrop1_24
           32.4d06,  65.8d06, 116.4d06, 204.0d06, 350.6d06, 584.4d06, 433.5d06, 433.5d06, &
           !ndrop1_25
           31.2d06,  62.3d06, 110.1d06, 191.3d06, 320.6d06, 501.3d06, 364.4d06, 364.4d06, &
           ! table4c (R2=0.02mum, wcb=2.5m/s) (for Ncn=50 and Ncn=100 "interpolated")
           !ndrop1_31
           50.3d06, 100.5d06, 201.1d06, 373.1d06, 664.7d06,1132.8d06,1876.8d06,2973.7d06, &
           !ndrop1_32
           44.1d06,  88.1d06, 176.2d06, 314.0d06, 546.9d06, 941.4d06,1579.2d06,2542.2d06, &
           !ndrop1_33
           39.7d06,  79.5d06, 158.9d06, 283.4d06, 498.9d06, 865.9d06,1462.6d06,2355.8d06, &
           !ndrop1_34
           37.0d06,  74.0d06, 148.0d06, 264.6d06, 468.3d06, 813.3d06,1371.3d06,2137.2d06, &
           !ndrop1_35
           34.7d06,  69.4d06, 138.8d06, 246.9d06, 432.9d06, 737.8d06,1176.7d06,1733.0d06, &
           ! table4d (R2=0.02mum, wcb=5.0m/s) (for Ncn=50,100,200 "interpolated")
           !ndrop1_41
           51.5d06, 103.1d06, 206.1d06, 412.2d06, 788.1d06,1453.1d06,2585.1d06,4382.5d06, &
           !ndrop1_42
           46.6d06,  93.2d06, 186.3d06, 372.6d06, 657.2d06,1202.8d06,2098.0d06,3556.9d06, &
           !ndrop1_43
           70.0d06,  70.0d06, 168.8d06, 337.6d06, 606.7d06,1078.5d06,1889.0d06,3206.9d06, &
           !ndrop1_44
           42.2d06,  84.4d06, 166.4d06, 312.7d06, 562.2d06,1000.3d06,1741.1d06,2910.1d06, &
           !ndrop1_45
           36.5d06,  72.9d06, 145.8d06, 291.6d06, 521.0d06, 961.1d06,1551.1d06,2444.6d06, &
           ! table5a (R2=0.03mum, wcb=0.5m/s)
           !ndrop2_11
           50.0d06,  95.8d06, 176.2d06, 321.6d06, 562.3d06, 835.5d06, 835.5d06, 835.5d06, &
           !ndrop2_12
           44.7d06,  81.4d06, 144.5d06, 251.5d06, 422.7d06, 677.8d06, 677.8d06, 677.8d06, &
           !ndrop2_13
           40.2d06,  72.8d06, 129.3d06, 225.9d06, 379.9d06, 606.5d06, 606.5d06, 606.5d06, &
           !ndrop2_14
           37.2d06,  67.1d06, 119.5d06, 206.7d06, 340.5d06, 549.4d06, 549.4d06, 549.4d06, &
           !ndrop2_15
           33.6d06,  59.0d06,  99.4d06, 150.3d06, 251.8d06, 466.0d06, 466.0d06, 466.0d06, &
           ! table5b (R2=0.03mum, wcb=1.0m/s) (Ncn=50 "interpolated", Ncn=6400 "extrapolated)
           !ndrop2_21
           50.7d06, 101.4d06, 197.6d06, 357.2d06, 686.6d06,1186.4d06,1892.2d06,1892.2d06, &
           !ndrop2_22
           46.6d06,  93.3d06, 172.2d06, 312.1d06, 550.7d06, 931.6d06,1476.6d06,1476.6d06, &
           !ndrop2_23
           42.2d06,  84.4d06, 154.0d06, 276.3d06, 485.6d06, 811.2d06,1271.7d06,1271.7d06, &
           !ndrop2_24
           39.0d06,  77.9d06, 141.2d06, 251.8d06, 436.7d06, 708.7d06,1117.7d06,1117.7d06, &
           !ndrop2_25
           35.0d06,  70.1d06, 123.9d06, 210.2d06, 329.9d06, 511.9d06, 933.4d06, 933.4d06, &
           ! table5c (R2=0.03mum, wcb=2.5m/s)
           !ndrop2_31
           51.5d06, 103.0d06, 205.9d06, 406.3d06, 796.4d06,1524.0d06,2781.4d06,4609.3d06, &
           !ndrop2_32
           49.6d06,  99.1d06, 198.2d06, 375.5d06, 698.3d06,1264.1d06,2202.8d06,3503.6d06, &
           !ndrop2_33
           45.8d06,  91.6d06, 183.2d06, 339.5d06, 618.9d06,1105.2d06,1881.8d06,2930.9d06, &
           !ndrop2_34
           42.3d06,  84.7d06, 169.3d06, 310.3d06, 559.5d06, 981.7d06,1611.6d06,2455.6d06, &
           !ndrop2_35
           38.2d06,  76.4d06, 152.8d06, 237.3d06, 473.3d06, 773.1d06,1167.9d06,1935.0d06, &
           ! table5d (R2=0.03mum, wcb=5.0m/s)
           !ndrop2_41
           51.9d06, 103.8d06, 207.6d06, 415.1d06, 819.6d06,1616.4d06,3148.2d06,5787.9d06, &
           !ndrop2_42
           50.7d06, 101.5d06, 203.0d06, 405.9d06, 777.0d06,1463.8d06,2682.6d06,4683.0d06, &
           !ndrop2_43
           47.4d06,  94.9d06, 189.7d06, 379.4d06, 708.7d06,1301.3d06,2334.3d06,3951.8d06, &
           !ndrop2_44
           44.0d06,  88.1d06, 176.2d06, 352.3d06, 647.8d06,1173.0d06,2049.7d06,3315.6d06, &
           !ndrop2_45
           39.7d06,  79.4d06, 158.8d06, 317.6d06, 569.5d06, 988.5d06,1615.6d06,2430.3d06, &
           ! table6a (R2=0.04mum, wcb=0.5m/s)
           !ndrop3_11
           50.6d06, 100.3d06, 196.5d06, 374.7d06, 677.3d06,1138.9d06,1138.9d06,1138.9d06, &
           !ndrop3_12
           48.4d06,  91.9d06, 170.6d06, 306.9d06, 529.2d06, 862.4d06, 862.4d06, 862.4d06, &
           !ndrop3_13
           44.4d06,  82.5d06, 150.3d06, 266.4d06, 448.0d06, 740.7d06, 740.7d06, 740.7d06, &
           !ndrop3_14
           40.9d06,  75.0d06, 134.7d06, 231.9d06, 382.1d06, 657.6d06, 657.6d06, 657.6d06, &
           !ndrop3_15
           34.7d06,  59.3d06,  93.5d06, 156.8d06, 301.9d06, 603.8d06, 603.8d06, 603.8d06, &
           ! table6b (R2=0.04mum, wcb=1.0m/s)
           !ndrop3_21
           50.9d06, 101.7d06, 201.8d06, 398.8d06, 773.7d06,1420.8d06,2411.8d06,2411.8d06, &
           !ndrop3_22
           49.4d06,  98.9d06, 189.7d06, 356.2d06, 649.5d06,1117.9d06,1805.2d06,1805.2d06, &
           !ndrop3_23
           45.6d06,  91.8d06, 171.5d06, 214.9d06, 559.0d06, 932.8d06,1501.6d06,1501.6d06, &
           !ndrop3_24
           42.4d06,  84.7d06, 155.8d06, 280.5d06, 481.9d06, 779.0d06,1321.9d06,1321.9d06, &
           !ndrop3_25
           36.1d06,  72.1d06, 124.4d06, 198.4d06, 319.1d06, 603.8d06,1207.6d06,1207.6d06, &
           ! table6c (R2=0.04mum, wcb=2.5m/s)
           !ndrop3_31
           51.4d06, 102.8d06, 205.7d06, 406.9d06, 807.6d06,1597.5d06,3072.2d06,5393.9d06, &
           !ndrop3_32
           50.8d06, 101.8d06, 203.6d06, 396.0d06, 760.4d06,1422.1d06,2517.4d06,4062.8d06, &
           !ndrop3_33
           48.2d06,  96.4d06, 193.8d06, 367.3d06, 684.0d06,1238.3d06,2087.3d06,3287.1d06, &
           !ndrop3_34
           45.2d06,  90.4d06, 180.8d06, 335.7d06, 611.2d06,1066.3d06,1713.4d06,2780.3d06, &
           !ndrop3_35
           38.9d06,  77.8d06, 155.5d06, 273.7d06, 455.2d06, 702.2d06,1230.7d06,2453.7d06, &
           ! table6d (R2=0.04mum, wcb=5.0m/s)
           !ndrop3_41
           53.1d06, 106.2d06, 212.3d06, 414.6d06, 818.3d06,1622.2d06,3216.8d06,6243.9d06, &
           !ndrop3_42
           51.6d06, 103.2d06, 206.3d06, 412.5d06, 805.3d06,1557.4d06,2940.4d06,5210.1d06, &
           !ndrop3_43
           49.6d06,  99.2d06, 198.4d06, 396.7d06, 755.5d06,1414.5d06,2565.3d06,4288.1d06, &
           !ndrop3_44
           46.5d06,  93.0d06, 186.0d06, 371.9d06, 692.9d06,1262.0d06,2188.3d06,3461.2d06, &
           !ndrop3_45
           39.9d06,  79.9d06, 159.7d06, 319.4d06, 561.7d06, 953.9d06,1493.9d06,2464.7d06 /), (/ n_ncn, 5, n_wcb, 3 /) )
      IF (i_lsigs < 1 .OR. i_lsigs > 5 .OR. i_r2 < 1 .OR. i_r2 > 3) THEN
        WRITE(0, '(a)') "!!!! wrong value for lsigs or R2 in cloud_nucleation_SK !!!!!"
      END IF
      DO i = 1, 4
        tab_ndrop(i, :) = ndrop(:, i_lsigs, i, i_r2)
      END DO
    END SUBROUTINE lookuptable

    !--------------

    FUNCTION ip_ndrop_ncn(tab_ndrop,tab_ncn,Ncn)

      ! Interpolation of the look-up table with respect to aerosol concentration Ncn
      REAL(wp) :: ip_ndrop_ncn(n_wcb)
      REAL(wp), INTENT(in) :: Ncn, tab_ndrop(n_wcb,n_ncn), tab_ncn(n_ncn)
      INTEGER            ::  ki
      LOGICAL :: found

      ! Interpolation of Ndrop from the values of Ncn in the look-up-table to given Ncn
      found = .FALSE.
      IF (Ncn <= tab_ncn(1)) THEN
        ip_ndrop_ncn = tab_ndrop(:,1)
        found = .TRUE.
      ELSE IF (Ncn >= tab_ncn(n_ncn)) then
        ip_ndrop_ncn = tab_ndrop(:,n_ncn)
        found = .TRUE.
      ELSE
        DO ki = 1,n_ncn-1
          IF (Ncn >= tab_ncn(ki) .AND.  Ncn <= tab_ncn(ki+1)) THEN
            ip_ndrop_ncn(:) = tab_ndrop(:,ki) + &
                 (tab_ndrop(:,ki+1)-tab_ndrop(:,ki)) / (tab_ncn(ki+1)-tab_ncn(ki)) * (Ncn-tab_ncn(ki))
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF (.NOT.found) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, ip_ndrop_ncn: lookup table interpolation failed! ')
      END IF

      RETURN

    END FUNCTION  ip_ndrop_ncn

    !------------------
    ! Interpolation of the interpolated look-up table with respect to w_cb
    FUNCTION ip_ndrop_wcb(tab_Ndrop_i,tab_wcb,wcb)
      REAL(wp) :: ip_ndrop_wcb
      REAL(wp), INTENT(in) :: tab_wcb(1:n_wcb), tab_Ndrop_i(1:n_wcb), wcb
      INTEGER            :: ki
      LOGICAL :: found

      !... Interpolation for Ndrop from the values of wcb in the look-up-tabel to detected wcb

      found = .FALSE.
      ip_ndrop_wcb = 0
      IF (wcb <= tab_wcb(1)) THEN
        ip_ndrop_wcb = tab_ndrop_i(1) / tab_wcb(1) * wcb
        found = .TRUE.
      ELSE IF (wcb  >= tab_wcb(n_wcb)) then
        ip_ndrop_wcb = tab_ndrop_i(n_wcb)
        found = .TRUE.
      ELSE
        DO ki = 1,n_wcb-1
          IF (wcb >= tab_wcb(ki) .AND.  wcb <= tab_wcb(ki+1)) THEN
            ip_ndrop_wcb = tab_ndrop_i(ki) + &
              (tab_ndrop_i(ki+1)-tab_ndrop_i(ki)) / (tab_wcb(ki+1)-tab_wcb(ki)) * (wcb-tab_wcb(ki))
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF (.NOT.found) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, ip_ndrop_wcb: lookup table interpolation failed! ')
      END IF

      RETURN

    END FUNCTION  ip_ndrop_wcb

  END SUBROUTINE ccn_activation_sk

 SUBROUTINE ccn_activation_hdcp2(ik_slice, atmo, cloud)
   !*******************************************************************************
   !       Calculation of ccn activation                                          *
   !       using the approach of Hande et al 2015                                 *
   !*******************************************************************************
    IMPLICIT NONE
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)

    TYPE(atmosphere), INTENT(inout) :: atmo
    CLASS(particle), INTENT(inout)  :: cloud

    ! Locale Variablen
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend

    INTEGER            :: i,k,nuc_typ
    REAL(wp)           :: n_c,q_c
    REAL(wp)           :: nuc_n,nuc_q
    REAL(wp)           :: wcb,pres
    REAL(wp)           :: acoeff,bcoeff,ccoeff,dcoeff
    REAL(wp), PARAMETER:: eps = 1e-20_wp

    ! Data from HDCP2_CCN_params.txt for 20130417
    REAL(wp), PARAMETER :: &
         a_ccn(4) = (/  183230691.161_wp, 0.10147358938_wp, &
         &             -0.2922395814_wp, 229189886.226_wp /), &
         b_ccn(4) = (/ 0.0001984051994_wp, 4.473190485e-05_wp, &
         &             0.0001843225275_wp, 0.0001986158191_wp /), &
         c_ccn(4) = (/ 16.2420263911_wp, 3.22011836758_wp, &
         &             13.8499423719_wp, 16.2461600644_wp /), &
         d_ccn(4) = (/ 287736034.13_wp, 0.6258809883_wp, &
         &             0.8907491812_wp, 360848977.55_wp /)

#ifdef _OPENACC
    CALL finish(routine, 'ccn_activation_hdcp2 not available on GPU')
#endif

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    nuc_typ = nuc_c_typ

    IF(isdebug) THEN
       WRITE(txt,*) "cloud_activation_hdcp2: nuc_typ = ",nuc_typ ; CALL message(routine,TRIM(txt))
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          nuc_q = 0.0d0
          nuc_n = 0.d0
          n_c   = cloud%n(i,k)
          q_c   = cloud%q(i,k)
          pres  = atmo%p(i,k)
          wcb   = atmo%w(i,k)

          if (q_c > eps .and. wcb > 0.0_wp) then

             ! Based on write-up of Luke Hande of 6 May 2015

             acoeff = a_ccn(1) * atan(b_ccn(1) * pres - c_ccn(1)) + d_ccn(1)
             bcoeff = a_ccn(2) * atan(b_ccn(2) * pres - c_ccn(2)) + d_ccn(2)
             ccoeff = a_ccn(3) * atan(b_ccn(3) * pres - c_ccn(3)) + d_ccn(3)
             dcoeff = a_ccn(4) * atan(b_ccn(4) * pres - c_ccn(4)) + d_ccn(4)

             nuc_n = acoeff * atan(bcoeff * log(wcb) + ccoeff) + dcoeff

             nuc_n = MAX(MAX(nuc_n,1.0e7_wp) - n_c,0.0_wp)

             nuc_q = MIN(nuc_n * cloud%x_min, atmo%qv(i,k))
             nuc_n = nuc_q / cloud%x_min

             cloud%n(i,k) = cloud%n(i,k) + nuc_n
             cloud%q(i,k) = cloud%q(i,k) + nuc_q
             atmo%qv(i,k) = atmo%qv(i,k) - nuc_q

          END IF

       END DO
    END DO

  END SUBROUTINE ccn_activation_hdcp2

  SUBROUTINE ccn_activation_sk_4d(ik_slice, ccn_coeffs, atmo, cloud, n_cn)
    !*******************************************************************************
    !       Calculation of cloud droplet nucleation                                *
    !       using the look-up tables by Segal and Khain 2006 (JGR, vol.11)         *
    !                                                                              *
    !       Difference to Heikes routine cloud_nucleation_SK()                     *
    !       Equidistant lookup table is used to enable better vectorization        *
    !       properties.                                                            *
    !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in), OPTIONAL :: ik_slice(4)

    ! parameters
    TYPE(aerosol_ccn),INTENT(in), OPTIONAL    :: ccn_coeffs

    ! 2mom variables
    TYPE(atmosphere), INTENT(inout), OPTIONAL :: atmo
    CLASS(particle),  INTENT(inout), OPTIONAL :: cloud
    REAL(wp), DIMENSION(:,:), OPTIONAL        :: n_cn

    ! local variables
    REAL(wp), PARAMETER :: nuc_eps = 1e-20_wp

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend

    ! grid sizes of the original table:
    INTEGER, PARAMETER   :: n_r2 = 3, n_lsigs = 5, n_ncn = 8 , n_wcb = 4

    ! desired grid sizes of the new equidistant table:
    INTEGER, PARAMETER   :: nr2  = 3, nlsigs  = 5, nncn  = 129, nwcb  = 11

    ! more local variables
    REAL(wp)             :: n_c, q_c
    REAL(wp)             :: nuc_n, nuc_q
    REAL(wp)             :: ncn, n_cn0, lsigs, nccn, r2, wcb, wcb_min
    REAL(wp)             :: r2_loc, lsigs_loc, ncn_loc, wcb_loc
    REAL(wp)             :: z0_nccn, z1e_nccn, zf, etas
    INTEGER              :: i, k, kp1_fl
    INTEGER              :: iu, ju, ku, lu
    REAL(wp)             :: hilf1(2,2,2,2), hilf2(2,2,2), hilf3(2,2), hilf4(2)

    LOGICAL, PARAMETER   :: lincloud_nuc = .TRUE.

    ! call from init_2mom_scheme_once without arguments for initialization of tables
    IF (.NOT.PRESENT(ik_slice)) THEN
      CALL get_otab(n_r2,n_lsigs,n_ncn,n_wcb)     ! original look-up-table from Segal and Khain      
      CALL equi_table(nr2,nlsigs,nncn,nwcb)   ! construct the new equidistant table tab:
      !$ACC ENTER DATA COPYIN(tab)
      !$ACC ENTER DATA COPYIN(tab%ltable, tab%x1, tab%x2, tab%x3, tab%x4)
      RETURN
    END IF

    IF(isdebug) THEN
       WRITE(txt,*) "cloud_activation_sk_4d'"
       CALL message(routine, TRIM(txt))
    ENDIF

    !..parameter for exponential decrease of N_ccn with height:
    !  1) up to this height (m) constant unchanged value:
    !  2)  height interval at which N_ccn decreses by factor 1/e above z0_nccn:

    z0_nccn  = ccn_coeffs%z0
    z1e_nccn = ccn_coeffs%z1e
    n_cn0    = ccn_coeffs%Ncn0
    etas     = ccn_coeffs%etas
    wcb_min  = ccn_coeffs%wcb_min

    !..values for aerosol properties
    r2    = ccn_coeffs%R2
    lsigs = ccn_coeffs%lsigs

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)
    
    !$ACC WAIT(1)
    !$ACC EXIT DATA DETACH(tab%ltable, tab%x1, tab%x2, tab%x3, tab%x4) FINALIZE
    !$ACC ENTER DATA ATTACH(tab%ltable, tab%x1, tab%x2, tab%x3, tab%x4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) CREATE(hilf1, hilf2, hilf3, hilf4)
    !$ACC LOOP SEQ
    DO k = kstart,kend
      kp1_fl = MIN(k+1,SIZE(atmo%rho,dim=2))
!NEC$ ivdep
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(n_c, q_c, nuc_n, nuc_q, ncn, nccn, wcb) &
      !$ACC   PRIVATE(r2_loc, lsigs_loc, ncn_loc, wcb_loc, zf, iu, ju, ku, lu) &
      !$ACC   PRIVATE(hilf1, hilf2, hilf3, hilf4)
      DO i = istart,iend

        ! hard upper limit for number conc that
        ! eliminates also unrealistic high value
        ! that would come from the dynamical core

        cloud%n(i,k) = MIN(cloud%n(i,k),n_cn0)

!!$ UB: determine wcb as in old COSMO version:
        ! determine vertical velocity for Segal&Khain nucleation parameterization:
        IF (lincloud_nuc) THEN
          ! ... incloud nucleation is allowed, look for height layers where qc (mass specific) increases with height and w is positive:
          !      (at the lowest model level, nucleation happens without the gradient check)
          IF ( cloud%q(i,k) > nuc_eps .AND. &
               ( k == kp1_fl .OR. cloud%q(i,k)/atmo%rho(i,k) > cloud%q(i,kp1_fl)/atmo%rho(i,kp1_fl) ) .AND. &
               atmo%w(i,k+1) > 0.0_wp ) THEN
            wcb = atmo%w(i,k+1)  ! take w of the lower cell face
          ELSE
            wcb = 0.0_wp ! set w for nucleation to 0.0, so that no new nucleation will take place below
          END IF
! We still miss the nucleation in fog situations during radiative cooling. This
! would require the inclusion of -cp/g*dT/dt|_diabatic in the effective
! nucleation velocity and allowing incloud nucleation everywhere, not only if qc
! increases with height.
        ELSE
          ! ... nucleation is allowed only in the model layer above cloud base:
          !      (the lowest model level always counts as cloud base if qc > nuc_eps and w > 0)
          IF ( (k == kp1_fl .OR. cloud%q(i,kp1_fl) <= nuc_eps) .AND. cloud%q(i,k) > nuc_eps .AND. atmo%w(i,k+1) > 0.0_wp) THEN
            wcb = atmo%w(i,k+1)    ! take w of the lower cell face
          ELSE
            wcb = 0.0_wp ! set w for nucleation to 0.0, so that no new nucleation will take place below
          END IF
        END IF

! previous formulation without the new lincloud_nuc mechanism:
!        IF (cloud%q(i,k) > eps .and. atmo%w(i,k) > 0.0_wp) THEN
! new formulation:
        IF ( wcb > 0.0_wp ) THEN

          nuc_q = 0.0_wp
          nuc_n = 0.0_wp
          n_c   = cloud%n(i,k)
          q_c   = cloud%q(i,k)
          wcb   = MAX(wcb, wcb_min)  ! enforce a minimal updraft for nucleation
 
          IF (PRESENT(n_cn)) THEN
            Ncn = n_cn(i,k) ! number of CN from prognostic variable
          ELSE
            zf = 0.5_wp*(atmo%zh(i,k)+atmo%zh(i,k+1))
            IF(zf > z0_nccn) THEN
              Ncn = n_cn0 * MIN(EXP((z0_nccn - zf)/z1e_nccn),1.0_wp)
            ELSE
              Ncn = n_cn0
            END IF
          END IF

          ! Interpolation of the look-up tables with respect to all 4 parameters:
          ! (clip values outside range to the marginal values)
          r2_loc    = MIN(MAX(r2,     tab%x1(1)), tab%x1(tab%n1))
          iu = MIN(FLOOR((r2_loc -    tab%x1(1)) * tab%odx1 ) + 1, tab%n1-1)
          lsigs_loc = MIN(MAX(lsigs,  tab%x2(1)), tab%x2(tab%n2))
          ju = MIN(FLOOR((lsigs_loc - tab%x2(1)) * tab%odx2 ) + 1, tab%n2-1)
          ncn_loc   = MIN(MAX(ncn,    tab%x3(1)), tab%x3(tab%n3))
          ku = MIN(FLOOR((ncn_loc -   tab%x3(1)) * tab%odx3 ) + 1, tab%n3-1)
          wcb_loc   = MIN(MAX(wcb,    tab%x4(1)), tab%x4(tab%n4))
          lu = MIN(FLOOR((wcb_loc -   tab%x4(1)) * tab%odx4 ) + 1, tab%n4-1)

          hilf1 = tab%ltable( iu:iu+1, ju:ju+1, ku:ku+1, lu:lu+1)
          hilf2 = hilf1(1,:,:,:) + (hilf1(2,:,:,:) - hilf1(1,:,:,:)) * tab%odx1 * ( r2_loc    - tab%x1(iu) )
          hilf3 = hilf2(1,:,:)   + (hilf2(2,:,:)   - hilf2(1,:,:)  ) * tab%odx2 * ( lsigs_loc - tab%x2(ju) )
          hilf4 = hilf3(1,:)     + (hilf3(2,:)     - hilf3(1,:)    ) * tab%odx3 * ( ncn_loc   - tab%x3(ku) )
          nccn  = hilf4(1)       + (hilf4(2)       - hilf4(1)      ) * tab%odx4 * ( wcb_loc   - tab%x4(lu) )

          ! If n_cn is outside the range of the lookup table values, resulting 
          ! NCCN are clipped to the margin values. For the case of these margin values
          ! beeing larger than n_cn (which happens sometimes, unfortunately), limit NCCN by n_cn:
          nccn = MIN(nccn, n_cn0)

          nuc_n = etas * nccn - n_c

          nuc_n = MAX(nuc_n,0.0d0)

          nuc_q = MIN(nuc_n * cloud%x_min,atmo%qv(i,k))
          nuc_n = nuc_q / cloud%x_min

          cloud%n(i,k) = cloud%n(i,k) + nuc_n
          cloud%q(i,k) = cloud%q(i,k) + nuc_q
          atmo%qv(i,k) = atmo%qv(i,k) - nuc_q

          IF (PRESENT(n_cn)) THEN
            n_cn(i,k) = n_cn(i,k) - MIN(Ncn,nuc_n)
          END IF

        END IF
      END DO
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE ccn_activation_sk_4d

  !*******************************************************************************
  ! Sedimentation subroutines for ICON
  !*******************************************************************************

  SUBROUTINE sedi_icon_rain (rain_in,rain_coeffs,qp,np,precrate,precrate3D,qc,rhocorr,adz,dt, &
      &                      its,ite,kts,kte,cmax,lacc)

    CLASS(particle), TARGET,INTENT(in)      :: rain_in
    TYPE(particle_rain_coeffs), INTENT(in)  :: rain_coeffs
    INTEGER,  INTENT(IN)                    :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: qp,np,precrate3D
    REAL(wp), DIMENSION(:,:), INTENT(IN)    :: adz,qc,rhocorr
    REAL(wp), DIMENSION(:),   INTENT(INOUT) :: precrate
    REAL(wp), INTENT(IN)                    :: dt
    REAL(wp), INTENT(INOUT), OPTIONAL       :: cmax
    CLASS(particle),POINTER                 :: rain ! ACCWA (nvhpc 22.7, IPSF, see above)

    INTEGER  :: i, k
    REAL(wp) :: x_p,D_m,D_p,mue,v_n,v_q

    REAL(wp), DIMENSION(its:ite,kts-1:kte) :: v_n_sedi,v_q_sedi

!!$ Activate the new explicit and more stable boxtracking sedimentation method:
!!$  (http://www.cosmo-model.org/content/model/documentation/core/docu_sedi_twomom.pdf)
    LOGICAL, PARAMETER :: lboxtracking = .true. ! lboxtracking = .false. is not supported on GPU

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    !-----------------------------------------------------------------------
    CALL assert_acc_device_only("sedi_icon_rain", lacc)

    !$ACC DATA CREATE(v_n_sedi, v_q_sedi)

    rain => rain_in ! ACCWA (nvhpc 22.7, IPSF, see above)
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO i = its,ite
      v_n_sedi(i,kts-1) = 0.0_wp
      v_q_sedi(i,kts-1) = 0.0_wp
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(x_p, D_m, mue, D_p, v_n, v_q)
    DO k = kts,kte
      DO i = its,ite

        IF (qp(i,k) > q_crit) THEN

          x_p = particle_meanmass(rain, qp(i,k) ,np(i,k))
          D_m = particle_diameter(rain, x_p)
          IF (cfg_params%luse_mu_Dm_rain) THEN
            IF (qc(i,k) >= q_crit) THEN
              mue = (rain%nu+1.0_wp)/rain%b_geo - 1.0_wp
            ELSE
              mue = rain_mue_dm_relation(rain_coeffs, D_m)
            END IF
          ELSE
            mue = (rain%nu+1.0_wp)/rain%b_geo - 1.0_wp
          END IF
          D_p = D_m * EXP((-1.0_wp/3.0_wp)*LOG((mue+3.0_wp)*(mue+2.0_wp)*(mue+1.0_wp)))

          v_n = rain_coeffs%alfa - rain_coeffs%beta * EXP(-(mue+1.)*LOG(1.0_wp + rain_coeffs%gama*D_p))
          v_q = rain_coeffs%alfa - rain_coeffs%beta * EXP(-(mue+4.)*LOG(1.0_wp + rain_coeffs%gama*D_p))
          v_q = MAX(v_q,  rain%vsedi_min)
          v_n = MAX(v_n,  rain%vsedi_min)
          v_n = v_n * rhocorr(i,k)
          v_q = v_q * rhocorr(i,k)

          v_n_sedi(i,k) = - v_n
          v_q_sedi(i,k) = - v_q
        ELSE
          v_n_sedi(i,k) = 0.0_wp
          v_q_sedi(i,k) = 0.0_wp
        ENDIF
      END DO
    END DO
    !$ACC END PARALLEL

    IF (lboxtracking) THEN
      CALL sedi_icon_box_core(v_n_sedi, v_q_sedi, adz, dt, its, ite, kts, kte, &
                              np, qp, precrate, precrate3D, cmax)
    ELSE
      CALL sedi_icon_core(v_n_sedi, v_q_sedi, adz, dt, its, ite, kts, kte, &
                          np, qp, precrate, precrate3D, cmax)
    END IF

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE sedi_icon_rain

  SUBROUTINE sedi_icon_sphere (ptype_in,pcoeffs,qp,np,precrate,precrate3D,rhocorr,adz,dt, &
      &                  its,ite,kts,kte,cmax,lacc)

    CLASS(particle),TARGET, INTENT(in)      :: ptype_in
    CLASS(particle_sphere), INTENT(in)      :: pcoeffs
    INTEGER, INTENT(IN)                     :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: qp,np,precrate3D
    REAL(wp), DIMENSION(:,:), INTENT(IN)    :: adz,rhocorr
    REAL(wp), DIMENSION(:),   INTENT(INOUT) :: precrate
    REAL(wp), INTENT(IN)                    :: dt
    REAL(wp), INTENT(INOUT), OPTIONAL       :: cmax
    CLASS(particle),POINTER                 :: ptype ! ACCWA (nvhpc 22.7, IPSF, see above)

    INTEGER  :: i, k
    REAL(wp) :: x_p,v_n,v_q,lam

    REAL(wp), DIMENSION(its:ite,kts-1:kte) :: v_n_sedi,v_q_sedi

!!$ Activate the new explicit and more stable boxtracking sedimentation method:
!!$  (http://www.cosmo-model.org/content/model/documentation/core/docu_sedi_twomom.pdf)
    LOGICAL, PARAMETER :: lboxtracking = .true. ! lboxtracking = .false. is not supported on GPU

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    !-----------------------------------------------------------------------
    CALL assert_acc_device_only("sedi_icon_sphere", lacc)

    !$ACC DATA CREATE(v_n_sedi, v_q_sedi)

    ptype => ptype_in ! ACCWA (nvhpc 22.7, IPSF, see above)


    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO i = its,ite
      v_n_sedi(i, kts-1) = 0.0_wp
      v_q_sedi(i, kts-1) = 0.0_wp
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(x_p, lam, v_n, v_q)
    DO k = kts,kte
      DO i = its,ite
        IF (qp(i,k) > q_crit) THEN

          x_p = particle_meanmass(ptype, qp(i,k),np(i,k))
          lam = EXP(ptype%b_vel* LOG(pcoeffs%coeff_lambda*x_p))

          v_n = pcoeffs%coeff_alfa_n * lam
          v_q = pcoeffs%coeff_alfa_q * lam
          v_n = MAX(v_n,ptype%vsedi_min)
          v_q = MAX(v_q,ptype%vsedi_min)
          v_n = MIN(v_n,ptype%vsedi_max)
          v_q = MIN(v_q,ptype%vsedi_max)
          v_n = v_n * rhocorr(i,k)
          v_q = v_q * rhocorr(i,k)

          v_n_sedi(i,k) = -v_n
          v_q_sedi(i,k) = -v_q
        ELSE
          v_n_sedi(i,k) = 0.0_wp
          v_q_sedi(i,k) = 0.0_wp
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    IF (lboxtracking) THEN
      CALL sedi_icon_box_core(v_n_sedi, v_q_sedi, adz, dt, its, ite, kts, kte, &
                              np, qp, precrate, precrate3D, cmax)
    ELSE
      CALL sedi_icon_core(v_n_sedi, v_q_sedi, adz, dt, its, ite, kts, kte, &
                          np, qp, precrate, precrate3D, cmax)
    END IF

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=kts,kte
      DO i = its,ite
        np(i,k) = MIN(np(i,k), qp(i,k)/ptype%x_min)
        np(i,k) = MAX(np(i,k), qp(i,k)/ptype%x_max)
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE sedi_icon_sphere

  SUBROUTINE sedi_icon_sphere_lwf (ptype,pcoeffs,qp,np,ql,precrate,precrate3D,rhocorr,adz,dt, &
      &                  its,ite,kts,kte,cmax) !

    TYPE(particle_lwf), INTENT(in)          :: ptype
    CLASS(particle_sphere), INTENT(in)      :: pcoeffs
    INTEGER, INTENT(IN)                     :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: qp,np,ql,precrate3D
    REAL(wp), DIMENSION(:,:), INTENT(IN)    :: adz,rhocorr
    REAL(wp), DIMENSION(:),   INTENT(INOUT) :: precrate
    REAL(wp), INTENT(IN)                    :: dt
    REAL(wp), INTENT(INOUT), OPTIONAL       :: cmax

    INTEGER  :: i, k
    REAL(wp) :: x_p,D_p,D_n,v_n,v_q,v_l,lam,lwf

    REAL(wp), DIMENSION(its:ite,kts-1:kte) :: v_n_sedi,v_q_sedi,v_ql_sedi

    REAL(wp), DIMENSION(10) :: avq, avl, avn
    REAL(wp), DIMENSION(9)  :: bvq, bvl, bvn

!!$ Activate the new explicit and more stable boxtracking sedimentation method:
!!$  (http://www.cosmo-model.org/content/model/documentation/core/docu_sedi_twomom.pdf)
    LOGICAL, PARAMETER :: lboxtracking = .true. ! lboxtracking = .false. is not supported on GPU

#ifdef _OPENACC
    CALL finish("sedi_icon_sphere_lwf", "Routine has not been ported to openACC yet.")
#endif

    if (ptype%name .eq. 'hail_vivek') then
      avq = aviwch
      bvq = bviwch
      avl = avlwch
      bvl = bvlwch
      avn = avnumh
      bvn = bvnumh
    elseif (ptype%name .eq. 'graupel_vivek') then
      avq = aviwcg
      bvq = bviwcg
      avl = avlwcg
      bvl = bvlwcg
      avn = avnumg
      bvn = bvnumg
    end if

    v_n_sedi(:, kts-1) = 0.0_wp
    v_q_sedi(:, kts-1) = 0.0_wp
    v_ql_sedi(:,kts-1) = 0.0_wp

    DO k = kts,kte
      DO i = its,ite

        IF (qp(i,k) > q_crit) THEN

          x_p = particle_meanmass(ptype,qp(i,k),np(i,k))
          lwf = particle_lwf_idx(ptype,i,k)

          IF (.false.) THEN
            lam = EXP(ptype%b_vel* LOG(pcoeffs%coeff_lambda*x_p))            
            v_n = pcoeffs%coeff_alfa_n * lam
            v_q = pcoeffs%coeff_alfa_q * lam
            v_l = pcoeffs%coeff_alfa_q * lam
          ELSE
            D_p = particle_diameter(ptype,x_p)
            D_n = particle_normdiameter(ptype,D_p)
            v_n = rat2do3(D_n,lwf,avn,bvn)
            v_q = rat2do3(D_n,lwf,avq,bvq)
            v_l = rat2do3(D_n,lwf,avl,bvl)
          END IF

          v_n = MAX(v_n,ptype%vsedi_min)
          v_q = MAX(v_q,ptype%vsedi_min)
          v_l = MAX(v_l,ptype%vsedi_min)
          v_n = MIN(v_n,ptype%vsedi_max)
          v_q = MIN(v_q,ptype%vsedi_max)
          v_l = MIN(v_l,ptype%vsedi_max)
          v_n = v_n * rhocorr(i,k)
          v_q = v_q * rhocorr(i,k)
          v_l = v_l * rhocorr(i,k)

          v_n_sedi(i,k)  = -v_n
          v_q_sedi(i,k)  = -v_q
          v_ql_sedi(i,k) = -v_l
        ELSE
          v_n_sedi(i,k) = 0.0_wp
          v_q_sedi(i,k) = 0.0_wp
          v_ql_sedi(i,k) = 0.0_wp
        END IF
      END DO
    END DO

    IF (lboxtracking) THEN
      CALL sedi_icon_box_core_lwf (v_n_sedi, v_q_sedi, v_ql_sedi, adz, dt, &
           &                       its, ite, kts, kte, np, qp, ql, precrate, precrate3D, cmax)
    ELSE

      CALL sedi_icon_core_lwf(v_n_sedi, v_q_sedi, v_ql_sedi, adz, dt, &
           &                  its, ite, kts, kte, np, qp, ql, precrate, precrate3D, cmax)

    END IF

    DO k=kts,kte
      DO i = its,ite
        np(i,k) = MIN(np(i,k), qp(i,k)/ptype%x_min)
        np(i,k) = MAX(np(i,k), qp(i,k)/ptype%x_max)
      END DO
    END DO

  END SUBROUTINE sedi_icon_sphere_lwf

  !*******************************************************************************
  !       Set to a default number concentration in places with qnx = 0 and qx !=0*
  !       (implemented by Alberto de Lozar)                                      *
  !*******************************************************************************
  SUBROUTINE set_default_n(ik_slice, cloud, ice, rain, snow, graupel, hail, n_cn)
    INTEGER, INTENT(in) :: ik_slice(4)
    CLASS(particle), INTENT(inout)      :: cloud
    CLASS(particle), INTENT(inout)      :: ice
    CLASS(particle), INTENT(inout)      :: rain
    CLASS(particle), INTENT(inout)      :: snow
    CLASS(particle), INTENT(inout)      :: graupel
    CLASS(particle), INTENT(inout)      :: hail
    REAL(wp), DIMENSION(:,:), OPTIONAL  :: n_cn
    LOGICAL                             :: n_cn_pres

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER :: i,k
    REAL(wp), PARAMETER :: eps = 1e-3_wp

    IF (PRESENT(n_cn)) THEN
      n_cn_pres = .TRUE.
    ELSE
      n_cn_pres = .FALSE.
    ENDIF

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO k = kstart,kend

      IF ( .NOT. n_cn_pres) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i = istart,iend
          IF ( cloud%q(i,k) > 0.0_wp .AND. cloud%n(i,k) < eps) THEN
            cloud%n(i,k) = set_qnc(cloud%q(i,k)) 
          END IF
        END DO
      END IF

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = istart,iend
        IF ( ice%q(i,k) > 0.0_wp .AND. ice%n(i,k) < eps) THEN
          ice%n(i,k) = set_qni(ice%q(i,k)) 
        END IF
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = istart,iend
        IF ( rain%q(i,k) > 0.0_wp .AND. rain%n(i,k) < eps) THEN
          rain%n(i,k) = set_qnr(rain%q(i,k)) 
        END IF
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = istart,iend
        IF ( snow%q(i,k) > 0.0_wp .AND. snow%n(i,k) < eps) THEN
          snow%n(i,k) = set_qns(snow%q(i,k)) 
        END IF
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = istart,iend
        IF ( graupel%q(i,k) > 0.0_wp .AND. graupel%n(i,k) < eps) THEN
          graupel%n(i,k) = set_qng(graupel%q(i,k)) 
        END IF
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = istart,iend
        IF ( hail%q(i,k) > 0.0_wp .AND. hail%n(i,k) < eps) THEN
          hail%n(i,k) = set_qnh_expPSD_N0const(hail%q(i,k),750.0_wp,1.0e6_wp) 
        END IF
      END DO

    END DO
    !$ACC END PARALLEL

  END SUBROUTINE set_default_n



  ! UB: This is according to the actual version from the COSMO code, a more efficient and
  !     simplified re-write of the above non-reproducible sedi_icon_core:

#if defined (__SX__) || defined (__NEC_VH__) || defined (__NECSX__)
  
  ! Vectorized version for NEC SX
  SUBROUTINE sedi_icon_core(v_n_sedi, v_q_sedi, adz, dt, its, ite, kts, kte, &
                            np, qp, precrate, precrate3D, cmax)

    REAL(wp), DIMENSION(:,:), INTENT(IN) :: adz
    REAL(wp), INTENT(in) :: dt
    INTEGER, INTENT(IN) :: its,ite,kts,kte
    REAL(wp), DIMENSION(its:ite,kts-1:kte), INTENT(in) :: v_n_sedi,v_q_sedi
    REAL(wp), INTENT(INOUT), OPTIONAL       :: cmax
    REAL(wp), DIMENSION(:), INTENT(INOUT)   :: precrate
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: np,qp,precrate3D

    REAL(wp), DIMENSION(its:ite, 0:1) :: q_fluss, n_fluss
    REAL(wp), DIMENSION(its:ite) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    REAL(wp), DIMENSION(its:ite,kts:kte) :: dz
    INTEGER :: i, k, kk, k_c, k_p
    REAL(wp) :: cmax_temp, odt, cmax_j

#ifdef _OPENACC
    CALL finish("mo_2mom_mcrph_processes", "sedi_icon_core is not supported on GPU")
#endif

    dz(its:ite,kts:kte) = 1.0_wp / adz(its:ite,kts:kte)

    odt = 1.0_wp / dt

    q_fluss(:, 1-IAND(kts, 1))  = 0.0_wp
    n_fluss(:, 1-IAND(kts, 1))  = 0.0_wp

    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF

    DO k = kts, kte

      k_c = IAND(k, 1)
      k_p = 1-IAND(k, 1)
      
      DO i = its,ite
        v_nv(i) = v_n_sedi(i,k)
        v_qv(i) = v_q_sedi(i,k)

        ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
        c_nv(i) = -v_nv(i) * adz(i,k) * dt
        c_qv(i) = -v_qv(i) * adz(i,k) * dt
      END DO
      cmax_j = MAXVAL( c_qv(its:ite) )
      cmax_temp = MAX(cmax_temp, cmax_j)

      kk = k
      DO i = its, ite
        s_nv(i) = np(i,k) * dz(i,k) * MIN(c_nv(i),1.0_wp)
      END DO
      DO
        IF (kk <= kts) EXIT
        IF (MAXVAL(c_nv) <= 1.0_wp) EXIT
        kk  = kk - 1
        DO i = its, ite
          IF (c_nv(i) > 1.0_wp) THEN
            c_nv(i) = (c_nv(i) - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
            s_nv(i) = s_nv(i) + np(i,kk) * dz(i,kk) * MIN(c_nv(i),1.0_wp)
          END IF
        END DO
      END DO

      kk = k
      DO i = its, ite
        s_qv(i) = qp(i,k) * dz(i,k) * MIN(c_qv(i),1.0_wp)
      END DO
      DO 
        IF (kk <= kts) EXIT
        IF (MAXVAL(c_qv) <= 1.0_wp) EXIT
        kk  = kk - 1
        DO i = its, ite
          IF (c_qv(i) > 1.0_wp) THEN
            c_qv(i) = (c_qv(i) - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
            s_qv(i) = s_qv(i) + qp(i,kk) * dz(i,kk) * MIN(c_qv(i),1.0_wp)
          END IF
        END DO
      END DO

      ! Flux-limiter to avoid negative values
      DO i = its,ite
        s_nv(i) = -s_nv(i) * odt
        s_qv(i) = -s_qv(i) * odt
        n_fluss(i,k_c) = MAX(s_nv(i),n_fluss(i,k_p)-np(i,k) * dz(i,k)*odt)
        q_fluss(i,k_c) = MAX(s_qv(i),q_fluss(i,k_p)-qp(i,k) * dz(i,k)*odt)
      END DO

      DO i = its,ite
        np(i,k) = np(i,k) + ( n_fluss(i,k_c) - n_fluss(i,k_p) )*adz(i,k)*dt
        qp(i,k) = qp(i,k) + ( q_fluss(i,k_c) - q_fluss(i,k_p) )*adz(i,k)*dt
        precrate3D(i,k) = - q_fluss(i,k_c) ! precipitation rate at lower level boundary        
      ENDDO

    END DO

    IF (PRESENT(cmax)) cmax = cmax_temp

    precrate(its:ite) = - q_fluss(its:ite,IAND(kte, 1)) ! precipitation rate at ground

  END SUBROUTINE sedi_icon_core

#else

  ! UB: scalar version for non-vector-architectures:
  SUBROUTINE sedi_icon_core(v_n_sedi, v_q_sedi, adz, dt, its, ite, kts, kte, &
                            np, qp, precrate, precrate3D, cmax)

    INTEGER, INTENT(IN) :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: adz
    REAL(wp), INTENT(in) :: dt
    REAL(wp), DIMENSION(its:ite,kts-1:kte), INTENT(in) :: v_n_sedi,v_q_sedi
    REAL(wp), INTENT(INOUT), OPTIONAL       :: cmax
    REAL(wp), DIMENSION(:), INTENT(INOUT)   :: precrate
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: np,qp,precrate3D

    REAL(wp), DIMENSION(its:ite, 0:1) :: q_fluss, n_fluss
    REAL(wp)                          :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    REAL(wp), DIMENSION(its:ite,kts:kte) :: dz
    INTEGER :: i, k, kk, k_c, k_p
    REAL(wp) :: cmax_temp, odt

#ifdef _OPENACC
    CALL finish("mo_2mom_mcrph_processes", "sedi_icon_core is not supported on GPU")
#endif

    dz(its:ite,kts:kte) = 1.0_wp / adz(its:ite,kts:kte)

    odt = 1.0_wp / dt

    q_fluss(:, 1-IAND(kts, 1))  = 0.0_wp
    n_fluss(:, 1-IAND(kts, 1))  = 0.0_wp

    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF

    DO k = kts, kte

      k_c = IAND(k, 1)
      k_p = 1-IAND(k, 1)

      DO i = its,ite

        v_nv = v_n_sedi(i,k)
        v_qv = v_q_sedi(i,k)

        ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
        c_nv = -v_nv * adz(i,k) * dt
        c_qv = -v_qv * adz(i,k) * dt
        cmax_temp = MAX(cmax_temp, c_qv)

        kk = k
        s_nv = np(i,kk) * dz(i,kk) * MIN(c_nv,1.0_wp);
        DO WHILE (c_nv > 1.0_wp .AND. kk > kts)
          kk = kk - 1
          c_nv = (c_nv - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
          s_nv = s_nv + np(i,kk) * dz(i,kk) * MIN(c_nv,1.0_wp)
        END DO
        s_nv = -s_nv * odt;

        kk = k
        s_qv = qp(i,kk) * dz(i,kk) * MIN(c_qv,1.0_wp);
        DO WHILE (c_qv > 1.0_wp .AND. kk > kts)
          kk = kk - 1
          c_qv = (c_qv - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
          s_qv = s_qv + qp(i,kk) * dz(i,kk) * MIN(c_qv,1.0_wp)
        END DO
        s_qv = -s_qv * odt;

        ! Flux-limiter to avoid negative values
        n_fluss(i,k_c) = MAX(s_nv,n_fluss(i,k_p)-np(i,k) * dz(i,k)*odt)
        q_fluss(i,k_c) = MAX(s_qv,q_fluss(i,k_p)-qp(i,k) * dz(i,k)*odt)

        np(i,k) = np(i,k) + ( n_fluss(i,k_c) - n_fluss(i,k_p) )*adz(i,k)*dt
        qp(i,k) = qp(i,k) + ( q_fluss(i,k_c) - q_fluss(i,k_p) )*adz(i,k)*dt

        precrate3D(i,k) = - q_fluss(i,k_c) ! precipitation rate at lower level boundary

      ENDDO
    END DO

    IF (PRESENT(cmax)) cmax = cmax_temp

    precrate(its:ite) = - q_fluss(its:ite,IAND(kte, 1)) ! precipitation rate at ground

  END SUBROUTINE sedi_icon_core
#endif

  ! UB: This is the new explicit and more stable boxtracking sedimentation method
  !  (http://www.cosmo-model.org/content/model/documentation/core/docu_sedi_twomom.pdf)

#if defined (__SX__) || defined (__NEC_VH__) || defined (__NECSX__)

  ! Vectorized version for the NEC:
  SUBROUTINE sedi_icon_box_core(v_n_sedi, v_q_sedi, adz, dt, its, ite, kts, kte, &
                                np, qp, precrate, precrate3D, cmax)

    INTEGER,  INTENT(IN)                               :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(IN)               :: adz
    REAL(wp), INTENT(in)                               :: dt
    REAL(wp), DIMENSION(its:ite,kts-1:kte), INTENT(in) :: v_n_sedi,v_q_sedi
    REAL(wp), INTENT(INOUT), OPTIONAL                  :: cmax
    REAL(wp), DIMENSION(:), INTENT(INOUT)              :: precrate
    REAL(wp), DIMENSION(:,:), INTENT(INOUT)            :: np,qp,precrate3D

    INTEGER                                :: i, k, kk, k_c, k_p
    REAL(wp), DIMENSION(its:ite, 0:1)      :: q_fluss, n_fluss
    REAL(wp), DIMENSION(its:ite)           :: v_nv, v_qv, c_qv, dz_loc
    REAL(wp), DIMENSION(its:ite,kts:kte)   :: s_nv, s_qv
    REAL(wp), DIMENSION(its:ite,kts:kte+1) :: dz
    REAL(wp)                               :: cmax_temp, cmax_j, odt

#ifdef _OPENACC
    CALL finish("mo_2mom_mcrph_processes", "The NEC version of sedi_icon_box_core is not supported on GPU")
#endif

    dz(its:ite,kts:kte) = 1.0_wp / adz(its:ite,kts:kte)
    ! .. dummy value for level kte+1, does not have an effect but is needed
    !    to prevent an array bound violation below:
    dz(its:ite,kte+1) = dz(its:ite,kte)

    odt = 1.0_wp / dt

    q_fluss(:, 1-IAND(kts, 1))  = 0.0_wp
    n_fluss(:, 1-IAND(kts, 1))  = 0.0_wp

    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF

    s_nv(:,:) = 0.0_wp
    s_qv(:,:) = 0.0_wp

    DO k = kts, kte
      
      DO i = its,ite

        v_nv(i) = v_n_sedi(i,k)
        v_qv(i) = v_q_sedi(i,k)

        ! Formulated under the assumption that v_nv, v_qv always negative:
        c_qv(i) = -v_qv(i) * adz(i,k) * dt
      END DO
      cmax_j = MAXVAL( c_qv(its:ite) )
      cmax_temp = MAX(cmax_temp, cmax_j)

      kk = 0              ! Loop index for the flux aggregation in the next boxes below the k'th box
      dz_loc(:) = 0.0_wp  ! Distance from the k'th lower cell face to the k+kk'th
                          !   lower cell face for downward processing starting from k
      DO
        IF ( k+kk > kte ) EXIT 
        IF ( ALL( dz_loc(:) >= -v_nv(:)*dt ) ) EXIT
        DO i = its, ite
          IF ( dz_loc(i) < -v_nv(i)*dt ) THEN
            s_nv(i,k+kk) = s_nv(i,k+kk) + np(i,k) * MIN(dz(i,k),-dz_loc(i)-v_nv(i)*dt)
            dz_loc(i) = dz_loc(i) + dz(i,k+kk+1)
          END IF
        END DO
        kk = kk + 1
      END DO

      ! .. The same for the time-averaged mass density flux:

      kk = 0
      dz_loc(:) = 0.0_wp
      DO
        IF ( k+kk > kte ) EXIT 
        IF ( ALL( dz_loc(:) >= -v_qv(:)*dt ) ) EXIT
        DO i = its, ite
          IF ( dz_loc(i) < -v_qv(i)*dt ) THEN
            s_qv(i,k+kk) = s_qv(i,k+kk) + qp(i,k) * MIN(dz(i,k),-dz_loc(i)-v_qv(i)*dt)
            dz_loc(i) = dz_loc(i) + dz(i,k+kk+1)
          END IF
        END DO
        kk = kk + 1
      END DO
      
    END DO

    ! .. Divide the time-aggregated flux by dt to get the time-averaged
    !    flux and give a negative sign because fluxes are directed downward:
    s_nv(:,:) = -s_nv(:,:) * odt
    s_qv(:,:) = -s_qv(:,:) * odt

    DO k = kts, kte

      k_c = IAND(k, 1)
      k_p = 1-IAND(k, 1)

      DO i = its,ite

        ! .. Flux-limiter to avoid negative values:
        n_fluss(i,k_c) = MAX(s_nv(i,k),n_fluss(i,k_p)-np(i,k) * dz(i,k)*odt)
        q_fluss(i,k_c) = MAX(s_qv(i,k),q_fluss(i,k_p)-qp(i,k) * dz(i,k)*odt)

        ! .. Update of nx and qx due to sedimenation flux divergences:
        np(i,k) = np(i,k) + ( n_fluss(i,k_c) - n_fluss(i,k_p) )*adz(i,k)*dt
        qp(i,k) = qp(i,k) + ( q_fluss(i,k_c) - q_fluss(i,k_p) )*adz(i,k)*dt

        precrate3D(i,k) = - q_fluss(i,k_c) ! precipitation rate at lower level boundary

      ENDDO
    END DO

    IF (PRESENT(cmax)) cmax = cmax_temp

    precrate(its:ite) = - q_fluss(its:ite,IAND(kte, 1)) ! precipitation rate at ground

  END SUBROUTINE sedi_icon_box_core

#else

  ! UB: scalar version for non-vector-architectures:
  SUBROUTINE sedi_icon_box_core(v_n_sedi, v_q_sedi, adz, dt, its, ite, kts, kte, &
                                np, qp, precrate, precrate3D, cmax)

    INTEGER, INTENT(IN)                          :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(IN)         :: adz
    REAL(wp), INTENT(in)                         :: dt
    REAL(wp), DIMENSION(its:ite,kts-1:kte), INTENT(in) :: v_n_sedi,v_q_sedi
    REAL(wp), INTENT(INOUT), OPTIONAL            :: cmax
    REAL(wp), DIMENSION(:), INTENT(INOUT)        :: precrate
    REAL(wp), DIMENSION(:,:), INTENT(INOUT)      :: np,qp,precrate3D

    INTEGER                                :: i, k, kk, k_c, k_p
    REAL(wp), DIMENSION(its:ite, 0:1)      :: q_fluss, n_fluss
    REAL(wp)                               :: v_nv, v_qv, c_qv
    REAL(wp), DIMENSION(its:ite,kts:kte)   :: s_nv, s_qv
    REAL(wp), DIMENSION(its:ite,kts:kte+1) :: dz
    REAL(wp)                               :: cmax_temp, odt, dz_loc

    !$ACC DATA CREATE(q_fluss, n_fluss, s_nv, s_qv, dz)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k = kts, kte
      DO i = its,ite
        dz(i,k) = 1.0_wp / adz(i,k)
        IF (k == kte) THEN
          ! .. dummy value for level kte+1, does not have an effect but is needed
          !    to prevent an array bound violation below:
          dz(i,k+1) = dz(i,k)
        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    odt = 1.0_wp / dt


    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO i = its,ite
      ! .. Upper boundary condition on fluxes:
      q_fluss(i, 1-IAND(kts, 1))  = 0.0_wp
      n_fluss(i, 1-IAND(kts, 1))  = 0.0_wp
    ENDDO
    !$ACC END PARALLEL

    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF

    CALL init(s_nv, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(s_qv, lacc=.TRUE., opt_acc_async=.TRUE.)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) REDUCTION(MAX: cmax_temp)
    !$ACC LOOP SEQ
    DO k = kts, kte
      !$ACC LOOP GANG VECTOR PRIVATE(v_nv, v_qv, kk, dz_loc) REDUCTION(MAX: cmax_temp)
      DO i = its,ite

        v_nv = v_n_sedi(i,k)
        v_qv = v_q_sedi(i,k)

        ! Formulated under the assumption that v_nv, v_qv always negative:
        c_qv = -v_qv * adz(i,k) * dt
        cmax_temp = MAX(cmax_temp, c_qv)

        kk = 0           ! Loop index for the flux aggregation in the next boxes below the k'th box
        dz_loc = 0.0_wp  ! Distance from the k'th lower cell face to the k+kk'th
                         !   lower cell face for downward processing starting from k
        DO WHILE ( dz_loc < -v_nv*dt .AND. k+kk <= kte)
          s_nv(i,k+kk) = s_nv(i,k+kk) + np(i,k) * MIN(dz(i,k),-dz_loc-v_nv*dt)
          kk = kk + 1
          dz_loc = dz_loc + dz(i,k+kk)
        END DO

        ! .. The same for the time-averaged mass density flux:

        kk = 0
        dz_loc = 0.0_wp
        DO WHILE ( dz_loc < -v_qv*dt .AND. k+kk <= kte)
          s_qv(i,k+kk) = s_qv(i,k+kk) + qp(i,k) * MIN(dz(i,k),-dz_loc-v_qv*dt)
          kk = kk + 1
          dz_loc = dz_loc + dz(i,k+kk)
        END DO

      END DO
    END DO
    !$ACC END PARALLEL

    ! .. Divide the time-aggregated flux by dt to get the time-averaged
    !    flux and give a negative sign because fluxes are directed downward:

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k = kts, kte
      DO i = its,ite
        s_nv(i,k) = -s_nv(i,k) * odt
        s_qv(i,k) = -s_qv(i,k) * odt
      END DO
    END DO
    !$ACC END PARALLEL


    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO k = kts, kte
      k_c = IAND(k, 1)
      k_p = 1-IAND(k, 1)
      !$ACC LOOP GANG VECTOR
      DO i = its,ite

        ! .. Flux-limiter to avoid negative values:
        n_fluss(i,k_c) = MAX(s_nv(i,k),n_fluss(i,k_p)-np(i,k) * dz(i,k)*odt)
        q_fluss(i,k_c) = MAX(s_qv(i,k),q_fluss(i,k_p)-qp(i,k) * dz(i,k)*odt)

        ! .. Update of nx and qx due to sedimenation flux divergences:
        np(i,k) = np(i,k) + ( n_fluss(i,k_c) - n_fluss(i,k_p) )*adz(i,k)*dt
        qp(i,k) = qp(i,k) + ( q_fluss(i,k_c) - q_fluss(i,k_p) )*adz(i,k)*dt

        precrate3D(i,k) = - q_fluss(i,k_c) ! precipitation rate at lower level boundary

      ENDDO
    END DO
    !$ACC END PARALLEL

    IF (PRESENT(cmax)) cmax = cmax_temp

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO i = its,ite
      precrate(i) = - q_fluss(i,IAND(kte, 1)) ! precipitation rate at ground
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE sedi_icon_box_core
#endif

#if defined (__SX__) || defined (__NEC_VH__) || defined (__NECSX__)

  ! Vectorized version for NEC:
  SUBROUTINE sedi_icon_core_lwf(v_n_sedi, v_q_sedi, v_ql_sedi, adz, dt, its, ite, kts, kte, &
                                np, qp, ql, precrate, precrate3D, cmax)

    INTEGER, INTENT(IN)                                :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(IN)               :: adz
    REAL(wp), INTENT(in)                               :: dt
    REAL(wp), DIMENSION(its:ite,kts-1:kte), INTENT(in) :: v_n_sedi,v_q_sedi,v_ql_sedi
    REAL(wp), INTENT(INOUT), OPTIONAL                  :: cmax
    REAL(wp), DIMENSION(:), INTENT(INOUT)              :: precrate
    REAL(wp), DIMENSION(:,:), INTENT(INOUT)            :: np,qp,ql,precrate3D

    INTEGER                              :: i, k, kk, k_c, k_p
    REAL(wp), DIMENSION(its:ite, 0:1)    :: q_fluss, n_fluss, ql_fluss
    REAL(wp), DIMENSION(its:ite)         :: v_nv, v_qv, v_ql, s_nv, s_qv, s_ql, c_nv, c_qv, c_ql
    REAL(wp), DIMENSION(its:ite,kts:kte) :: dz
    REAL(wp)                             :: cmax_temp, odt, cmax_j

    dz(its:ite,kts:kte) = 1.0_wp / adz(its:ite,kts:kte)

    odt = 1.0_wp / dt

    q_fluss(:, 1-IAND(kts, 1))  = 0.0_wp
    ql_fluss(:, 1-IAND(kts, 1))  = 0.0_wp
    n_fluss(:, 1-IAND(kts, 1))  = 0.0_wp

    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF

    DO k = kts, kte

      DO i = its,ite
        v_nv(i) = v_n_sedi(i,k)
        v_qv(i) = v_q_sedi(i,k)
        v_ql(i) = v_ql_sedi(i,k)

        c_nv(i) = -v_nv(i) * adz(i,k) * dt
        c_qv(i) = -v_qv(i) * adz(i,k) * dt
        c_ql(i) = -v_ql(i) * adz(i,k) * dt
      END DO
      cmax_j = MAXVAL( c_qv(its:ite) )
      cmax_temp = MAX(cmax_temp, cmax_j)

      kk = k
      DO i = its, ite
        s_nv(i) = np(i,k) * dz(i,k) * MIN(c_nv(i),1.0_wp)
      END DO
      DO
        IF (kk <= kts) EXIT
        IF (MAXVAL(c_nv) <= 1.0_wp) EXIT
        kk  = kk - 1
        DO i = its, ite
          IF (c_nv(i) > 1.0_wp) THEN
            c_nv(i) = (c_nv(i) - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
            s_nv(i) = s_nv(i) + np(i,kk) * dz(i,kk) * MIN(c_nv(i),1.0_wp)
          END IF
        END DO
      END DO

      kk = k
      DO i = its, ite
        s_qv(i) = qp(i,k) * dz(i,k) * MIN(c_qv(i),1.0_wp)
      END DO
      DO 
        IF (kk <= kts) EXIT
        IF (MAXVAL(c_qv) <= 1.0_wp) EXIT
        kk  = kk - 1
        DO i = its, ite
          IF (c_qv(i) > 1.0_wp) THEN
            c_qv(i) = (c_qv(i) - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
            s_qv(i) = s_qv(i) + qp(i,kk) * dz(i,kk) * MIN(c_qv(i),1.0_wp)
          END IF
        END DO
      END DO

      kk = k
      DO i = its, ite
        s_ql(i) = ql(i,k) * dz(i,k) * MIN(c_ql(i),1.0_wp)
      END DO
      DO 
        IF (kk <= kts) EXIT
        IF (MAXVAL(c_ql) <= 1.0_wp) EXIT
        kk  = kk - 1
        DO i = its, ite
          IF (c_ql(i) > 1.0_wp) THEN
            c_ql(i) = (c_ql(i) - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
            s_ql(i) = s_ql(i) + ql(i,kk) * dz(i,kk) * MIN(c_ql(i),1.0_wp)
          END IF
        END DO
      END DO

      ! Flux-limiter to avoid negative values
      k_c = IAND(k, 1)
      k_p = 1-IAND(k, 1)
      
      DO i = its,ite

        s_nv(i) = -s_nv(i) * odt
        s_qv(i) = -s_qv(i) * odt
        s_ql(i) = -s_ql(i) * odt
        n_fluss(i,k_c)  = MAX(s_nv(i),n_fluss (i,k_p)-np(i,k) * dz(i,k)*odt)
        q_fluss(i,k_c)  = MAX(s_qv(i),q_fluss (i,k_p)-qp(i,k) * dz(i,k)*odt)
        ql_fluss(i,k_c) = MAX(s_ql(i),ql_fluss(i,k_p)-ql(i,k) * dz(i,k)*odt)

      END DO

      DO i = its,ite
        
        np(i,k) = np(i,k) + ( n_fluss(i,k_c)  - n_fluss (i,k_p) )*adz(i,k)*dt
        qp(i,k) = qp(i,k) + ( q_fluss(i,k_c)  - q_fluss (i,k_p) )*adz(i,k)*dt
        ql(i,k) = ql(i,k) + ( ql_fluss(i,k_c) - ql_fluss(i,k_p) )*adz(i,k)*dt
        
        precrate3D(i,k) = - q_fluss(i,k_c) - ql_fluss(i,k_c) ! precipitation rate at lower level boundary

      ENDDO

    END DO

    IF (PRESENT(cmax)) cmax = cmax_temp

    precrate(its:ite) = - q_fluss(its:ite,IAND(kte, 1)) - ql_fluss(its:ite,IAND(kte, 1)) 

  END SUBROUTINE sedi_icon_core_lwf

#else

  ! UB: scalar version for non-vector-architectures:
  SUBROUTINE sedi_icon_core_lwf(v_n_sedi, v_q_sedi, v_ql_sedi, adz, dt, its, ite, kts, kte, &
                                np, qp, ql, precrate, precrate3D, cmax)

    INTEGER, INTENT(IN) :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: adz
    REAL(wp), INTENT(in) :: dt
    REAL(wp), DIMENSION(its:ite,kts-1:kte), INTENT(in) :: v_n_sedi,v_q_sedi,v_ql_sedi
    REAL(wp), INTENT(INOUT), OPTIONAL       :: cmax
    REAL(wp), DIMENSION(:), INTENT(INOUT)   :: precrate
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: np,qp,ql,precrate3D

    REAL(wp), DIMENSION(its:ite, 0:1) :: q_fluss, n_fluss, ql_fluss
    REAL(wp)                          :: v_nv, v_qv, v_ql, s_nv, s_qv, s_ql, c_nv, c_qv, c_ql
    REAL(wp), DIMENSION(its:ite,kts:kte) :: dz
    INTEGER :: i, k, kk, k_c, k_p
    REAL(wp) :: cmax_temp, odt

    dz(its:ite,kts:kte) = 1.0_wp / adz(its:ite,kts:kte)

    odt = 1.0_wp / dt

    q_fluss(:, 1-IAND(kts, 1))  = 0.0_wp
    ql_fluss(:, 1-IAND(kts, 1))  = 0.0_wp
    n_fluss(:, 1-IAND(kts, 1))  = 0.0_wp

    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF

    DO k = kts, kte
      DO i = its,ite

        v_nv = v_n_sedi(i,k)
        v_qv = v_q_sedi(i,k)
        v_ql = v_ql_sedi(i,k)

        c_nv = -v_nv * adz(i,k) * dt
        c_qv = -v_qv * adz(i,k) * dt
        c_ql = -v_ql * adz(i,k) * dt
        cmax_temp = MAX(cmax_temp, c_qv)

        kk = k
        s_nv = np(i,kk) * dz(i,kk) * MIN(c_nv,1.0_wp);
        DO WHILE (c_nv > 1.0_wp .AND. kk > kts)
          kk = kk - 1
          c_nv = (c_nv - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
          s_nv = s_nv + np(i,kk) * dz(i,kk) * MIN(c_nv,1.0_wp)
        END DO
        s_nv = -s_nv * odt;

        kk = k
        s_qv = qp(i,kk) * dz(i,kk) * MIN(c_qv,1.0_wp);
        DO WHILE (c_qv > 1.0_wp .AND. kk > kts)
          kk = kk - 1
          c_qv = (c_qv - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
          s_qv = s_qv + qp(i,kk) * dz(i,kk) * MIN(c_qv,1.0_wp)
        END DO
        s_qv = -s_qv * odt;

        kk = k
        s_ql = ql(i,kk) * dz(i,kk) * MIN(c_ql,1.0_wp);
        DO WHILE (c_ql > 1.0_wp .AND. kk > kts)
          kk = kk - 1
          c_ql = (c_ql - 1.0_wp) * adz(i,kk) * dz(i,kk+1)
          s_ql = s_ql + ql(i,kk) * dz(i,kk) * MIN(c_ql,1.0_wp)
        END DO
        s_ql = -s_ql * odt;

        ! Flux-limiter to avoid negative values
        k_c = IAND(k, 1)
        k_p = 1-IAND(k, 1)
        n_fluss(i,k_c)  = MAX(s_nv,n_fluss (i,k_p)-np(i,k) * dz(i,k)*odt)
        q_fluss(i,k_c)  = MAX(s_qv,q_fluss (i,k_p)-qp(i,k) * dz(i,k)*odt)
        ql_fluss(i,k_c) = MAX(s_ql,ql_fluss(i,k_p)-ql(i,k) * dz(i,k)*odt)

        np(i,k) = np(i,k) + ( n_fluss (i,k_c) - n_fluss (i,k_p) )*adz(i,k)*dt
        qp(i,k) = qp(i,k) + ( q_fluss (i,k_c) - q_fluss (i,k_p) )*adz(i,k)*dt
        ql(i,k) = ql(i,k) + ( ql_fluss(i,k_c) - ql_fluss(i,k_p) )*adz(i,k)*dt
        
        precrate3D(i,k) = - q_fluss(i,k_c) - ql_fluss(i,k_c) ! precipitation rate at lower level boundary

      ENDDO
    END DO

    IF (PRESENT(cmax)) cmax = cmax_temp

    precrate(its:ite) = - q_fluss(its:ite,IAND(kte, 1)) - ql_fluss(its:ite,IAND(kte, 1))

  END SUBROUTINE sedi_icon_core_lwf
#endif


  ! UB: This is the new explicit and more stable boxtracking sedimentation method
  !  (http://www.cosmo-model.org/content/model/documentation/core/docu_sedi_twomom.pdf)

#if defined (__SX__) || defined (__NEC_VH__) || defined (__NECSX__)

  ! Vectorized version for NEC:
  SUBROUTINE sedi_icon_box_core_lwf(v_n_sedi, v_q_sedi, v_ql_sedi, adz, dt, its, ite, kts, kte, &
                                    np, qp, ql, precrate, precrate3D, cmax)

    INTEGER,  INTENT(IN)                               :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(IN)               :: adz
    REAL(wp), INTENT(in)                               :: dt
    REAL(wp), DIMENSION(its:ite,kts-1:kte), INTENT(in) :: v_n_sedi,v_q_sedi,v_ql_sedi
    REAL(wp), INTENT(INOUT), OPTIONAL                  :: cmax
    REAL(wp), DIMENSION(:), INTENT(INOUT)              :: precrate
    REAL(wp), DIMENSION(:,:), INTENT(INOUT)            :: np,qp,ql,precrate3D

    INTEGER                                :: i, k, kk, k_c, k_p
    REAL(wp), DIMENSION(its:ite, 0:1)      :: q_fluss, n_fluss, ql_fluss
    REAL(wp), DIMENSION(its:ite)           :: v_nv, v_qv, v_ql, c_qv, dz_loc
    REAL(wp), DIMENSION(its:ite,kts:kte)   :: s_nv, s_qv, s_ql
    REAL(wp), DIMENSION(its:ite,kts:kte+1) :: dz
    REAL(wp)                               :: cmax_temp, cmax_j, odt

    dz(its:ite,kts:kte) = 1.0_wp / adz(its:ite,kts:kte)
    ! .. dummy value for level kte+1, does not have an effect but is needed
    !    to prevent an array bound violation below:
    dz(its:ite,kte+1) = dz(its:ite,kte)

    odt = 1.0_wp / dt

    ! .. Upper boundary condition on fluxes:
    q_fluss(:, 1-IAND(kts, 1))   = 0.0_wp
    ql_fluss(:, 1-IAND(kts, 1))  = 0.0_wp
    n_fluss(:, 1-IAND(kts, 1))   = 0.0_wp

    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF

    s_nv(:,:) = 0.0_wp
    s_qv(:,:) = 0.0_wp
    s_ql(:,:) = 0.0_wp

    DO k = kts, kte

      DO i = its,ite

        v_nv(i) = v_n_sedi(i,k)
        v_qv(i) = v_q_sedi(i,k)
        v_ql(i) = v_ql_sedi(i,k)

        ! Formulated under the assumption that v_nv, v_qv always negative:
        c_qv(i) = -v_qv(i) * adz(i,k) * dt
      END DO
      cmax_j = MAXVAL( c_qv(its:ite) )
      cmax_temp = MAX(cmax_temp, cmax_j)

      kk = 0              ! Loop index for the flux aggregation in the next boxes below the k'th box
      dz_loc(:) = 0.0_wp  ! Distance from the k'th lower cell face to the k+kk'th
                          !   lower cell face for downward processing starting from k
      DO
        IF ( k+kk > kte ) EXIT 
        IF ( ALL( dz_loc(:) >= -v_nv(:)*dt ) ) EXIT
        DO i = its, ite
          IF ( dz_loc(i) < -v_nv(i)*dt ) THEN
            s_nv(i,k+kk) = s_nv(i,k+kk) + np(i,k) * MIN(dz(i,k),-dz_loc(i)-v_nv(i)*dt)
            dz_loc(i) = dz_loc(i) + dz(i,k+kk+1)
          END IF
        END DO
        kk = kk + 1
      END DO

      ! .. The same for the mass density flux:
      kk = 0
      dz_loc(:) = 0.0_wp
      DO
        IF ( k+kk > kte ) EXIT 
        IF ( ALL( dz_loc(:) >= -v_qv(:)*dt ) ) EXIT
        DO i = its, ite
          IF ( dz_loc(i) < -v_qv(i)*dt ) THEN
            s_qv(i,k+kk) = s_qv(i,k+kk) + qp(i,k) * MIN(dz(i,k),-dz_loc(i)-v_qv(i)*dt)
            dz_loc(i) = dz_loc(i) + dz(i,k+kk+1)
          END IF
        END DO
        kk = kk + 1
      END DO
      
      ! .. The same for the liquid partial mass density flux:
      kk = 0
      dz_loc(:) = 0.0_wp
      DO
        IF ( k+kk > kte ) EXIT 
        IF ( ALL( dz_loc(:) >= -v_ql(:)*dt ) ) EXIT
        DO i = its, ite
          IF ( dz_loc(i) < -v_ql(i)*dt ) THEN
            s_ql(i,k+kk) = s_ql(i,k+kk) + ql(i,k) * MIN(dz(i,k),-dz_loc(i)-v_ql(i)*dt)
            dz_loc(i) = dz_loc(i) + dz(i,k+kk+1)
          END IF
        END DO
        kk = kk + 1
      END DO
      
    END DO

    ! .. Divide the time-aggregated flux by dt to get the time-averaged
    !    flux and give a negative sign because fluxes are directed downward:
    s_nv(:,:) = -s_nv(:,:) * odt
    s_qv(:,:) = -s_qv(:,:) * odt
    s_ql(:,:) = -s_ql(:,:) * odt

    DO k = kts, kte

      k_c = IAND(k, 1)
      k_p = 1-IAND(k, 1)

      DO i = its,ite

        ! .. Flux-limiter to avoid negative values:
        n_fluss(i,k_c)  = MAX(s_nv(i,k),n_fluss (i,k_p)-np(i,k) * dz(i,k)*odt)
        q_fluss(i,k_c)  = MAX(s_qv(i,k),q_fluss (i,k_p)-qp(i,k) * dz(i,k)*odt)
        ql_fluss(i,k_c) = MAX(s_ql(i,k),ql_fluss(i,k_p)-ql(i,k) * dz(i,k)*odt)

        ! .. Update of nx and qx due to sedimenation flux divergences:
        np(i,k) = np(i,k) + ( n_fluss (i,k_c) - n_fluss (i,k_p) )*adz(i,k)*dt
        qp(i,k) = qp(i,k) + ( q_fluss (i,k_c) - q_fluss (i,k_p) )*adz(i,k)*dt
        ql(i,k) = ql(i,k) + ( ql_fluss(i,k_c) - ql_fluss(i,k_p) )*adz(i,k)*dt

        precrate3D(i,k) = - q_fluss(i,k_c) - ql_fluss(i,k_c) ! precipitation rate at lower level boundary

      ENDDO
    END DO

    IF (PRESENT(cmax)) cmax = cmax_temp

    precrate(its:ite) = - q_fluss(its:ite,IAND(kte, 1)) - ql_fluss(its:ite,IAND(kte, 1))

  END SUBROUTINE sedi_icon_box_core_lwf

#else

  ! UB: scalar version for non-vector-architectures:
  SUBROUTINE sedi_icon_box_core_lwf(v_n_sedi, v_q_sedi, v_ql_sedi, adz, dt, its, ite, kts, kte, &
                                    np, qp, ql, precrate, precrate3D, cmax)

    INTEGER, INTENT(IN)                                :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(IN)               :: adz
    REAL(wp), INTENT(in)                               :: dt
    REAL(wp), DIMENSION(its:ite,kts-1:kte), INTENT(in) :: v_n_sedi,v_q_sedi,v_ql_sedi
    REAL(wp), INTENT(INOUT), OPTIONAL                  :: cmax
    REAL(wp), DIMENSION(:), INTENT(INOUT)              :: precrate
    REAL(wp), DIMENSION(:,:), INTENT(INOUT)            :: np,qp,ql,precrate3D

    INTEGER                                :: i, k, kk, k_c, k_p
    REAL(wp), DIMENSION(its:ite, 0:1)      :: q_fluss, n_fluss, ql_fluss
    REAL(wp)                               :: v_nv, v_qv, v_ql, c_qv
    REAL(wp), DIMENSION(its:ite,kts:kte)   :: s_nv, s_qv, s_ql
    REAL(wp), DIMENSION(its:ite,kts:kte+1) :: dz
    REAL(wp)                               :: cmax_temp, odt, dz_loc

    dz(its:ite,kts:kte) = 1.0_wp / adz(its:ite,kts:kte)

    odt = 1.0_wp / dt

    ! .. Upper boundary condition on fluxes:
    q_fluss(:, 1-IAND(kts, 1))   = 0.0_wp
    ql_fluss(:, 1-IAND(kts, 1))  = 0.0_wp
    n_fluss(:, 1-IAND(kts, 1))   = 0.0_wp

    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF

    s_nv(:,:) = 0.0_wp
    s_qv(:,:) = 0.0_wp
    s_ql(:,:) = 0.0_wp

    DO k = kts, kte
      DO i = its,ite

        v_nv = v_n_sedi(i,k)
        v_qv = v_q_sedi(i,k)
        v_ql = v_ql_sedi(i,k)

        c_qv = -v_qv * adz(i,k) * dt
        cmax_temp = MAX(cmax_temp, c_qv)

        kk = 0           ! Loop index for the flux aggregation in the next boxes below the k'th box
        dz_loc = 0.0_wp  ! Distance from the k'th lower cell face to the k+kk'th
                         !   lower cell face for downward processing starting from k
        DO WHILE ( dz_loc < -v_nv*dt .AND. k+kk <= kte)
          s_nv(i,k+kk) = s_nv(i,k+kk) + np(i,k) * MIN(dz(i,k),-dz_loc-v_nv*dt)
          kk = kk + 1
          dz_loc = dz_loc + dz(i,k+kk)
        END DO

        ! .. The same for the time-averaged mass density flux:

        kk = 0
        dz_loc = 0.0_wp
        DO WHILE ( dz_loc < -v_qv*dt .AND. k+kk <= kte)
          s_qv(i,k+kk) = s_qv(i,k+kk) + qp(i,k) * MIN(dz(i,k),-dz_loc-v_qv*dt)
          kk = kk + 1
          dz_loc = dz_loc + dz(i,k+kk)
        END DO

        kk = 0
        dz_loc = 0.0_wp
        DO WHILE ( dz_loc < -v_ql*dt .AND. k+kk <= kte)
          s_ql(i,k+kk) = s_ql(i,k+kk) + ql(i,k) * MIN(dz(i,k),-dz_loc-v_ql*dt)
          kk = kk + 1
          dz_loc = dz_loc + dz(i,k+kk)
        END DO

      END DO
    END DO

    ! .. Divide the time-aggregated flux by dt to get the time-averaged
    !    flux and give a negative sign because fluxes are directed downward:

    s_nv(:,:) = -s_nv(:,:) * odt
    s_qv(:,:) = -s_qv(:,:) * odt
    s_ql(:,:) = -s_ql(:,:) * odt


    DO k = kts, kte

      k_c = IAND(k, 1)
      k_p = 1-IAND(k, 1)

      DO i = its,ite

        ! .. Flux-limiter to avoid negative values:
        n_fluss(i,k_c)  = MAX(s_nv(i,k),n_fluss (i,k_p)-np(i,k) * dz(i,k)*odt)
        q_fluss(i,k_c)  = MAX(s_qv(i,k),q_fluss (i,k_p)-qp(i,k) * dz(i,k)*odt)
        ql_fluss(i,k_c) = MAX(s_ql(i,k),ql_fluss(i,k_p)-ql(i,k) * dz(i,k)*odt)

        ! .. Update of nx and qx due to sedimenation flux divergences:
        np(i,k) = np(i,k) + ( n_fluss (i,k_c) - n_fluss (i,k_p) )*adz(i,k)*dt
        qp(i,k) = qp(i,k) + ( q_fluss (i,k_c) - q_fluss (i,k_p) )*adz(i,k)*dt
        ql(i,k) = ql(i,k) + ( ql_fluss(i,k_c) - ql_fluss(i,k_p) )*adz(i,k)*dt

        precrate3D(i,k) = - q_fluss(i,k_c) - ql_fluss(i,k_c) ! precipitation rate at lower level boundary

      ENDDO
    END DO

    IF (PRESENT(cmax)) cmax = cmax_temp

    precrate(its:ite) = - q_fluss(its:ite,IAND(kte, 1)) - ql_fluss(its:ite,IAND(kte, 1))

  END SUBROUTINE sedi_icon_box_core_lwf
#endif

END MODULE mo_2mom_mcrph_processes
