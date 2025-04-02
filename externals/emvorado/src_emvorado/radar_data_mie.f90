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

#if defined TWOMOM_SB_NEW || defined TWOMOM_SB_OLD
#define TWOMOM_SB
#endif

MODULE radar_data_mie

!------------------------------------------------------------------------------
!
! Description:
!      Data declarations for the model specific interface(s) (at the moment only COSMO)
!      to the EMVORADO libraries for radar reflectivity calculations based on Mie-Theory
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_kind, ONLY : dp, wp
  
!!$ UB: this is unfortunately not possible other frameworks than COSMO:
!!$  USE data_constants, ONLY : T0_melt

!===============================================================================

  IMPLICIT NONE

!===============================================================================

  PUBLIC

  !----------------------------------------------------------------------------
  ! Version string for all lookup tables (rain, ice, snow, graupel, hail).
  ! Increase counter by 1, if something relevant has changed in the lookup table
  ! setup, such as:
  !  - Number of nodes and scaling of the table vectors
  !  - Particle models for melting stuff
  !  - Hydrometeor shape and orientation parametrizations
  !  - Adaptions/bugfixes in the basic Mie or Tmatrix codes
  !  - Format changes of the Mie table files
  !----------------------------------------------------------------------------

  CHARACTER(len=*), PARAMETER :: versionstring_lookup = 'version020'
  ! Lookuptable dimensions
  !                      #x-pts (2mom) #q-pts (1mom) #T-pts (dry) #T-pts (wet) #Tmelt-pts
  INTEGER, PARAMETER ::  nx_r = 40,    nq_r = 40,    nTa_r  = 10,               nmuD  = 40, &
                         nx_i = 30,    nq_i = 30,    nTa_id = 10, nTa_iw = 30 , nTm_i = 15, &
                         nx_s = 35,    nq_s = 35,    nTa_sd = 20, nTa_sw = 100, nTm_s = 50, &
                         nx_g = 40,    nq_g = 40,    nTa_gd = 20, nTa_gw = 100, nTm_g = 50, &
                         nx_h = 40,                  nTa_hd = 20, nTa_hw = 100, nTm_h = 50

  ! Size grid determining parameters
  REAL(dp), PARAMETER :: Dmin_r =  50.0d-6, &  ! lower limit of integration of reflectivity for: rain
                         Dmax_r =  10.0d-3, &  ! upper ------------------- " ------------------- --- " ---
                         Dmin_i =   5.0d-6, &  ! lower ------------------- " ------------------- cloud ice
                         Dmax_i =  10.0d-3, &  ! upper ------------------- " ------------------- --- " ---
                         Dmin_s =  50.0d-6, &  ! lower ------------------- " ------------------- snow
                         Dmax_s =  50.0d-3, &  ! upper ------------------- " ------------------- --- " ---
                         Dmin_g =  10.0d-6, &  ! lower ------------------- " ------------------- graupel
                         Dmax_g =  30.0d-3, &  ! upper ------------------- " ------------------- --- " ---
                         Dmin_h =  50.0d-6, &  ! lower ------------------- " ------------------- hail
                         Dmax_h = 100.0d-3     ! upper ------------------- " ------------------- --- " ---
  INTEGER, PARAMETER  :: n_stuetz = 250        ! Number size integration bins, equi-width distributed
                                               ! between Dmin and Dmax.
                                               !!! Has to be even (Simpson rule) !!!

! Local scalars:
!---------------

  REAL   (dp), PARAMETER :: quasi_zero = 1e-20_dp
  REAL   (dp), PARAMETER :: Deps = 2e-12_dp ! low bound on particle diameter to do e.g. Mie calcs on
  REAL   (dp), PARAMETER :: xeps = 1e-12_dp ! low bound on particle size parameter (alf,nue) to do e.g. Mie calcs on
  REAL   (dp), PARAMETER :: qeps = 1e-12_dp ! low bound on mass content


  REAL   (dp), PARAMETER :: third      = 1.0_dp/3.0_dp

  REAL   (dp), PARAMETER :: pi_dp      = 4.0_dp * ATAN(1.0_dp)
  REAL   (dp), PARAMETER :: pih_dp     = pi_dp*0.5_dp
  REAL   (dp), PARAMETER :: pi6_dp     = pi_dp/6.0_dp
  REAL   (dp), PARAMETER :: inv_pi6_dp = 6.0_dp/pi_dp
  REAL   (dp), PARAMETER :: inv_pi6sq_dp = inv_pi6_dp*inv_pi6_dp

  REAL   (dp), PARAMETER :: degrad_dp = pi_dp / 180.0_dp
  REAL   (dp), PARAMETER :: raddeg_dp = 180.0_dp / pi_dp

  COMPLEX(dp), PARAMETER :: nci_dp = (0.0_dp, 1.0_dp) ! negative complex i

  COMPLEX(dp), PARAMETER :: m_air       = (1.0_dp, 0.0_dp)
  REAL   (dp), PARAMETER :: T0C_fwo     = 273.15_dp ! 0C in K
  REAL   (dp), PARAMETER :: rho_0       = 1.225_dp  ! kg/m3
  REAL   (dp), PARAMETER :: rho_w_fwo   = 1000.0_dp ! kg/m3
  REAL   (dp), PARAMETER :: rho_ice_fwo = 900.0_dp  ! kg/m3
  REAL   (dp), PARAMETER :: inv_rhow    = 1.0_dp/rho_w_fwo
  REAL   (dp), PARAMETER :: inv_rhow2   = inv_rhow*inv_rhow
  REAL   (dp), PARAMETER :: inv_rhoi    = 1.0_dp/rho_ice_fwo
  REAL   (dp), PARAMETER :: inv_rhoi2   = inv_rhoi*inv_rhoi

  REAL   (dp), PARAMETER :: mw_Tmin = 0.0_dp   ! min temp[°C] for complex refindex calcs for water
  REAL   (dp), PARAMETER :: mw_Tmax = 30.0_dp  ! max temp[°C] for complex refindex calcs for water
  REAL   (dp), PARAMETER :: mi_Tmin = -80.0_dp ! min temp[°C] for complex refindex calcs for ice
  REAL   (dp), PARAMETER :: mi_Tmax = 0.0_dp   ! max temp[°C] for complex refindex calcs for ice

  ! .. Melting/freezing point temperature of ice/water [K]
  !    Used as Tmin in melting scheme.
!!$ This is unfortunately not possible in other frameworks than COSMO. Changed back to previous solution T0C_fwo:
!!$  REAL   (dp), PARAMETER :: Tmin_f = T0_melt
  REAL   (dp), PARAMETER :: Tmin_f = T0C_fwo

  ! .. Local domain bounds for computations (are set in the interface routines
  !    to the COSMO-model, e.g., radar_mie_2mom_vec(), ...):
  INTEGER :: ilow_modelgrid, iup_modelgrid, jlow_modelgrid, jup_modelgrid, klow_modelgrid, kup_modelgrid

  ! .. 2D fields for the temperature of actual melt start due to wet growth or T0 for graupel and hail:
  REAL(wp), ALLOCATABLE :: Tmin_g_modelgrid(:,:), Tmin_h_modelgrid(:,:)

  ! .. 2D fields for the temperature of "just-beeing-melted" for each hydrometeor type:
  REAL(wp), ALLOCATABLE :: Tmax_i_modelgrid(:,:), Tmax_s_modelgrid(:,:), &
                           Tmax_g_modelgrid(:,:), Tmax_h_modelgrid(:,:)

  ! .. Debug flag with default initialization:
  !    This parameter has to be set externally by the calling program,
  !    at the moment this is calc_dbz_vec(), and
  !    this subroutine itself receives the flag via its argument list:
  LOGICAL :: ldebug_dbz = .TRUE.

  ! .. Konfiguration switches, depending on the model's microphysics.
  !    Have to be mapped to the actual model namelist parameters by
  !    an interface routine. Currently, this is in get_model_config_for_radar()
  !    in radar_interface.f90:
  LOGICAL :: lgsp_fwo         ! Microphysics are on (.true.) or off (.false.)
  INTEGER :: itype_gscp_fwo   ! Type of microphysics parameterization
  LOGICAL :: luse_muD_relation_rain_fwo ! for 2mom-scheme: if if the mu-Dm-Relation of Seifert (2008) is applied outside the cloud cores
#ifdef __COSMO__
  INTEGER :: klv850_fwo       ! Index of the vertical level which is approx. at 850 hPa,
                              !  needed to define the level for DBZ_850 output.
#endif

  ! Derived type to hold the hydrometeor class-specific parameters:
  TYPE particle
    CHARACTER(20)    :: name  !..Bezeichnung der Partikelklasse
    REAL(dp) :: nu    !..Exp.-parameter der Verteil.
    REAL(dp) :: mu    !..Breiteparameter der Verteil.
    REAL(dp) :: x_max !..maximale Teilchenmasse
    REAL(dp) :: x_min !..minimale Teilchenmasse
    REAL(dp) :: a_geo !..Koeff. Geometrie
    REAL(dp) :: b_geo !..Koeff. Geometrie = 1/3
    REAL(dp) :: a_vel !..Koeff. Fallgesetz
    REAL(dp) :: b_vel !..Koeff. Fallgesetz
    REAL(dp) :: a_ven !..Koeff. Ventilationsparam.
    REAL(dp) :: b_ven !..Koeff. Ventilationsparam.
    REAL(dp) :: n0_const !..Const. n0 (optional; only relevant for (some) 1mom hydromets)
  END TYPE particle

  TYPE(particle), TARGET :: cloud, rain, ice, snow, graupel, hail

  ! Derived type to hold the MGD parameters:
  TYPE t_mgd_params
    ! for MGD with N(x) = n0 * x^mu * exp(-lam * x^nu)
    ! q and qn are total mass and number (concentration), respectively.
    REAL(dp) :: n0
    REAL(dp) :: mu
    REAL(dp) :: lam
    REAL(dp) :: nu
    REAL(dp) :: q
    REAL(dp) :: qn
  END TYPE t_mgd_params

  ! Derived type for the mu-D-relation in the 2-moment scheme:
  TYPE particle_rain_coeffs
    REAL(dp)   :: cmu0   !..Parameters for mu-D-relation of rain
    REAL(dp)   :: cmu1   !     max of left branch
    REAL(dp)   :: cmu2   !     max of right branch
    REAL(dp)   :: cmu3   !     min value of relation
    REAL(dp)   :: cmu4   !     location of min value = breakup equilibrium diameter
    INTEGER    :: cmu5   !     exponent
  END TYPE particle_rain_coeffs

  TYPE(particle_rain_coeffs) :: rain_coeffs

  ! Usage of structure entries:
  ! ---
  ! 'rayleigh' for any rayleigh calculations
  ! 'liquid' and 'frozen' in on-the-run Mie calculations for liquid and
  !    (dry&wet) frozen hydrometeor species
  !
  ! NOTE: Since cloud so far is calculated from Rayleigh approach,
  !       it uses always 'rayleigh' limit.
  TYPE qlimits
    REAL(dp) :: rayleigh
    REAL(dp) :: liquid
    REAL(dp) :: frozen
  END TYPE qlimits

  ! Type declaration for a general 4D equidistant lookup table. Used, e.g., for
  !  wet growth diameter lookup table
  TYPE lookupt_4D
    LOGICAL :: is_initialized = .FALSE.
    INTEGER :: iflag = -HUGE(1) ! general-purpose flag
    CHARACTER(len=40) :: name   ! general-purpose name
    INTEGER :: n1  ! number of grid points in x1-direction
    INTEGER :: n2  ! number of grid points in x2-direction
    INTEGER :: n3  ! number of grid points in x3-direction
    INTEGER :: n4  ! number of grid points in x4-direction
    REAL(dp), DIMENSION(:), POINTER :: x1 => NULL()  ! grid vector in x1-direction
    REAL(dp), DIMENSION(:), POINTER :: x2 => NULL()  ! grid vector in x1-direction
    REAL(dp), DIMENSION(:), POINTER :: x3 => NULL()  ! grid vector in x1-direction
    REAL(dp), DIMENSION(:), POINTER :: x4 => NULL()  ! grid vector in x1-direction
    REAL(dp)                     :: dx1          ! dx1   (grid distance w.r.t. x1)
    REAL(dp)                     :: dx2          ! dx2   (grid distance w.r.t. x2)
    REAL(dp)                     :: dx3          ! dx3   (grid distance w.r.t. x3)
    REAL(dp)                     :: dx4          ! dx4   (grid distance w.r.t. x4)
    REAL(dp)                     :: odx1         ! one over dx1
    REAL(dp)                     :: odx2         ! one over dx2
    REAL(dp)                     :: odx3         ! one over dx3
    REAL(dp)                     :: odx4         ! one over dx4
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: ltable => NULL()
  END TYPE lookupt_4D

  TYPE(qlimits), PARAMETER :: q_crit_radar = qlimits(1d-7, 1d-7, 1d-7)
  REAL(dp)     , PARAMETER :: n_crit_radar = 1d-3

  !===============================================================================
  !===============================================================================

  ! Flags for indicating the scaling of qi-values and table values:
  INTEGER, PARAMETER :: i_scal_log=1, i_scal_fscal=2, i_scal_lin=3, i_scal_dbz=4

  ! Flags for indicating the type of interpolation with respect to the scaled q_i:
  INTEGER, PARAMETER :: i_interp_lin=1, i_interp_cubic=2

  ! Zero-values and critical values to distinguish from zero-value for the different scalings of lookup table values:
  REAL(dp), PARAMETER :: zero_value_lut(4) = &
       [LOG10(quasi_zero), quasi_zero, 0.0_dp, 10.0_dp*LOG10(quasi_zero)]
  REAL(dp), PARAMETER :: scal_crit_lut(4) = &
       [LOG10(quasi_zero)+1.0_dp, 10.0_dp*quasi_zero, quasi_zero, 10.0_dp*(LOG10(quasi_zero)+1.0_dp)]
  
  ! Type for usage within dbzlookuptable to hold the pointer to a sub-block of memory mem for one parameter:
  TYPE t_tabparam
    REAL(dp), DIMENSION(:)    , POINTER   :: q_i => NULL()       ! value of the scaled q_i (q_i,T_a,T_m), points to qmem(:,i_<param>)
    REAL(dp), DIMENSION(:,:,:), POINTER   :: val => NULL()      ! value of the param at (q_i,T_a,T_m), points to mem(:,:,:,i_<param>)
    REAL(dp), DIMENSION(:,:,:), POINTER   :: dval => NULL()     ! value of the derivative dparam/dqi
    INTEGER,  POINTER                     :: qi_scal => NULL()  ! flag for desired scaling of qi for interpolation, points to flag_qi_scal(i_<param>)
    REAL(dp), POINTER                     :: f_qi => NULL(), if_qi => NULL()   ! exponent and its inverse for desired qi-scaling in interpolation (if qi_scal = i_scal_fscal)
    INTEGER,  POINTER                     :: val_scal => NULL() ! flag for type of scaling of val and dval, points to flag_mem_scal(i_<param>)
    REAL(dp), POINTER                     :: f_val => NULL(), if_val => NULL() ! exponent and its inverse, which have been used for scaling val and dval  (if val_scal = i_scal_fscal)
    INTEGER,  POINTER                     :: val_interp => NULL() ! flag for type of interpolation with respect to scaled qi
  END TYPE t_tabparam

  !     Type declaration for a lookup table calculating Z, ext, and polarimetric parameters
  !     as functions of (T_a,q_i,T_m) with q_i = q_r,q_s,q_g (rain,snow,graupel)
  TYPE t_dbzlookuptable
    ! config parameters for block of memory, do not change:
    INTEGER                                 :: i_zh=1, i_ah=2, i_zv=3, i_rrhv=4, i_irhv=5, i_kdp=6, i_adp=7, i_zvh=8
    INTEGER                                 :: nparams = 8
    ! other config parameters, partly with default init:
    LOGICAL                                 :: is_initialized = .FALSE.
    LOGICAL                                 :: luse_tmatrix = .TRUE.
    LOGICAL                                 :: ldo_nonsphere = .TRUE.
    INTEGER                                 :: itype_refl = 1
    INTEGER                                 :: itype_Dref_fmelt = -99
    INTEGER                                 :: magicnr = -1365623035 + HUGE(1) ! hash2 for text(1:3000) = ' '
    CHARACTER(len=20)                       :: cversion_lt  ! Version string of the table, to be defined by the developers
    CHARACTER(len=20)                       :: chydrotype   ! Hydrometeor type identifier
    CHARACTER(len=200)                      :: chydroconfig ! String representation of hydrometeor particle type parameters
    INTEGER                                 :: nTa      ! number of lookup bins of temperature
    INTEGER                                 :: nqi      ! number of lookup bins of q_i
    INTEGER                                 :: nTm      ! number of lookup bins of Tmax (needed for melting particles)
    REAL(dp)                                :: dTa      ! increment of temperature
    REAL(dp)                                :: dqi      ! increment of q_i (f-scaled)
    REAL(dp)                                :: dTm      ! increment of Tmax
    REAL(dp)                                :: idTa     ! 1/dTa
    REAL(dp)                                :: idqi     ! 1/dqi (f-scaled)
    REAL(dp)                                :: idTm     ! 1/dTm
    REAL(dp)                                :: Ta0      ! initial value of temperature
    REAL(dp)                                :: qi0      ! initial value of q_i (f-scaled)
    REAL(dp)                                :: Tm0      ! initial value of Tmax
    INTEGER                                 :: flag_qi_scal_eq ! flag for scaling type of equidistant qi table vector
    REAL(dp)                                :: f_eq     ! exponent of q_i scaling for equidistant qi table vector (if flag_qi_scal_eq = i_scal_fscal)
    REAL(dp)                                :: if_eq    ! 1/f_eq
    REAL(dp), DIMENSION(:), ALLOCATABLE     :: T_a      ! vector of temperature-values
    REAL(dp), DIMENSION(:), ALLOCATABLE     :: q_i      ! vector of equidistant q_i-values in the desired scaling;
                                                        !  flag flag_qi_scal_eq indicates the scaling: 1=log, 2=q_i^f_eq, 3=lin
    REAL(dp), DIMENSION(:), ALLOCATABLE     :: q_i_lin  ! vector of (linear) q_i-values
    REAL(dp), DIMENSION(:), ALLOCATABLE     :: T_m      ! vector of Tmax
!!$ UB: POINTER instead of ALLOCATABLE, so that pointers can point to them:
!!$     (Unfortunately no ALLOCATABLE, TARGET allowed in derived types!)
    REAL(dp), DIMENSION(:,:),     POINTER   :: qmem => NULL()   ! block of memory for the scaled qi table vectors: qmem(nqi,nparams)
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: mem => NULL()    ! block of memory for the table values: mem(nqi,nTa,nTm,nparams)
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: dmem => NULL()   ! block of memory for the table values of the derivative dmem/dqi
    INTEGER,  DIMENSION(:),       POINTER   :: flag_qi_scal => NULL()  ! for each param in mem (last index), store the qi scaling flag
    REAL(dp), DIMENSION(:),       POINTER   :: f_scal_qi => NULL()
    REAL(dp), DIMENSION(:),       POINTER   :: if_scal_qi => NULL()
    INTEGER,  DIMENSION(:),       POINTER   :: flag_mem_scal => NULL() ! for each param in mem (last index), store the parameter scaling flag
    REAL(dp), DIMENSION(:),       POINTER   :: f_scal_mem => NULL()
    REAL(dp), DIMENSION(:),       POINTER   :: if_scal_mem => NULL()
    INTEGER,  DIMENSION(:),       POINTER   :: flag_mem_interp_qi => NULL() ! for each param in mem (last index), store the interp method w.r.t. scaled qi
   
    TYPE(t_tabparam) :: zh   ! value of H-polarization radar reflectivity at (q_i,T_a,T_m)
    TYPE(t_tabparam) :: ah   ! value of H-polarization extinction at (q_i,T_a,T_m)
    TYPE(t_tabparam) :: zv   ! value of V-polarization radar reflectivity at (q_i,T_a,T_m)
    TYPE(t_tabparam) :: rrhv ! value of real part of ShhSvv* at (q_i,T_a,T_m)
    TYPE(t_tabparam) :: irhv ! value of imag part of ShhSvv* at (q_i,T_a,T_m)
    TYPE(t_tabparam) :: kdp  ! value of spec. diff. phase at (q_i,T_a,T_m)
    TYPE(t_tabparam) :: adp  ! value of specific diff. attenuation coeff. at (q_i,T_a,T_m)
    TYPE(t_tabparam) :: zvh  ! value of H-to-V-polarization radar reflectivity at (q_i,T_a,T_m)
  END TYPE t_dbzlookuptable

  TYPE t_tabledef
    INTEGER       :: nqi, nTa, nTm
    REAL(KIND=dp) :: qilow, qiup
  END TYPE t_tabledef

  !===============================================================================

  ! Maximum number of configs for polarimetric parameter computation, for which lookup tables can be created:
  INTEGER, PARAMETER :: nmax_lookup = 50

  ! Vectors of lookup table types for each hydrometeor type and each dBZ-configuration
  !  (a maximum of nmax_lookup configurations is supported):
  TYPE(t_dbzlookuptable), TARGET :: look_Z_rain(nmax_lookup)
  TYPE(t_dbzlookuptable), TARGET :: look_Z_rain_muD(nmax_lookup)
  TYPE(t_dbzlookuptable), TARGET :: look_Z_ice(nmax_lookup)
  TYPE(t_dbzlookuptable), TARGET :: look_Z_snow(nmax_lookup)
  TYPE(t_dbzlookuptable), TARGET :: look_Z_graupel(nmax_lookup)
  TYPE(t_dbzlookuptable), TARGET :: look_Z_hail(nmax_lookup)
  TYPE(t_dbzlookuptable), TARGET :: look_Z_meltice(nmax_lookup)
  TYPE(t_dbzlookuptable), TARGET :: look_Z_meltsnow(nmax_lookup)
  TYPE(t_dbzlookuptable), TARGET :: look_Z_meltgraupel(nmax_lookup)       ! table for actual Tmeltbegin_g from namelist
  TYPE(t_dbzlookuptable), TARGET :: look_Z_meltgraupel_0C(nmax_lookup) ! table for Tmeltbegin_g = 273.15 K
  TYPE(t_dbzlookuptable), TARGET :: look_Z_melthail(nmax_lookup)       ! table for actual Tmeltbegin_h from namelist
  TYPE(t_dbzlookuptable), TARGET :: look_Z_melthail_0C(nmax_lookup) ! table for Tmeltbegin_h = 273.15 K

  ! Lookup tables for wet growth diameters:
  TYPE(lookupt_4d), TARGET :: look_Dmin_wg_graupel
  TYPE(lookupt_4d), TARGET :: look_Dmin_wg_hail
  
  !===============================================================================
  !===============================================================================

  ! Type to hold constants and prefactors for the Rayleigh-Oguchi-Approximation:
  TYPE rayleigh_consts
    COMPLEX(dp) :: m_w_0
    REAL(dp)    :: K_w_0
    REAL(dp)    :: pi6sq_K0_fac
    REAL(dp)    :: C_fac
    REAL(dp)    :: Z_fac
    ! .. For cloud droplets:
    REAL(dp)    :: cloud_Zprefac
    REAL(dp)    :: cloud_Zprefac_rho
    REAL(dp)    :: cloud_extprefac
    ! .. For rain drops:
    REAL(dp)    :: rain2mom_Zprefac_rho
    REAL(dp)    :: rain1mom_Zprefac_n0_rho
    REAL(dp)    :: rain1mom_expo
    ! .. For cloud ice:
    REAL(dp)    :: ice1mom_Zprefac
    REAL(dp)    :: ice1mom_Zprefac_rho
    REAL(dp)    :: ice2mom_Zprefac
    REAL(dp)    :: ice2mom_Zprefac_rho
    REAL(dp)    :: ice_Davfac
    ! .. For snow:
    REAL(dp)    :: snow2mom_Zprefac
    REAL(dp)    :: snow2mom_Zprefac_rho
    REAL(dp)    :: snow1mom_Zprefac
    REAL(dp)    :: snow1mom_Zprefac_rho
    REAL(dp)    :: snow1mom_expo
    REAL(dp)    :: snow_Davfac
    REAL(dp)    :: snow_gam_mub1nu
    ! .. For graupel:
    REAL(dp)    :: graupel2mom_Zprefac
    REAL(dp)    :: graupel2mom_Zprefac_rho
    REAL(dp)    :: graupel1mom_Zprefac_n0
    REAL(dp)    :: graupel1mom_Zprefac_n0_rho
    REAL(dp)    :: graupel1mom_expo
    REAL(dp)    :: graupel_gam_mub1nu
    REAL(dp)    :: graupel_Davfac
    ! .. For hail:
    REAL(dp)    :: hail2mom_Zprefac
    REAL(dp)    :: hail2mom_Zprefac_rho
    REAL(dp)    :: hail_Davfac
  END TYPE rayleigh_consts

  TYPE(rayleigh_consts) :: ray_const

  ! Type to hold constants and prefactors for the radial winds with Rayleigh-Oguchi-Approximation:
  TYPE vt_oguchi_consts
    REAL(dp)    :: K_w_0
    REAL(dp)    :: am      ! blending exponent for the fall velocity of partially melted particles

    ! .. For cloud droplets:
    REAL(dp)    :: cloud_vtprefac
    REAL(dp)    :: cloud_vtexpo

    ! .. For rain drops:
    REAL(dp)    :: rain_vtexpo
    REAL(dp)    :: rain_nexpo
    REAL(dp)    :: rain_vtprefac
    REAL(dp)    :: rain_nprefac

    ! .. For cloud ice:
    REAL(dp)    :: ice_vtexpo
    REAL(dp)    :: ice_vtlexpo
    REAL(dp)    :: ice_vtprefac
    REAL(dp)    :: ice_vtfprefac
    REAL(dp)    :: ice_vtlprefac
    REAL(dp)    :: ice_Davfac

    ! .. For snow:
    REAL(dp)    :: snow_vtexpo
    REAL(dp)    :: snow_vtlexpo
    REAL(dp)    :: snow_nexpo
    REAL(dp)    :: snow_vtprefac
    REAL(dp)    :: snow_vtfprefac
    REAL(dp)    :: snow_vtlprefac
    REAL(dp)    :: snow_nprefac
    REAL(dp)    :: snow_Davfac
    REAL(dp)    :: snow_gam_mub1nu ! Factor gcft((mu+bgeo+1)/nu), needed diameter of mean mass for exp. distr.

    ! .. For graupel:
    REAL(dp)    :: graupel_vtexpo
    REAL(dp)    :: graupel_vtlexpo
    REAL(dp)    :: graupel_nexpo
    REAL(dp)    :: graupel_vtprefac
    REAL(dp)    :: graupel_vtfprefac
    REAL(dp)    :: graupel_vtlprefac
    REAL(dp)    :: graupel_nprefac
    REAL(dp)    :: graupel_Davfac
    REAL(dp)    :: graupel_gam_mub1nu ! Factor gcft((mu+bgeo+1)/nu), needed diameter of mean mass for exp. distr.

    ! .. For hail:
    REAL(dp)    :: hail_vtexpo
    REAL(dp)    :: hail_vtlexpo
    REAL(dp)    :: hail_vtprefac
    REAL(dp)    :: hail_vtfprefac
    REAL(dp)    :: hail_vtlprefac
    REAL(dp)    :: hail_Davfac
  END TYPE vt_oguchi_consts

  TYPE(vt_oguchi_consts) :: vt_const

  !===============================================================================
  !===============================================================================
  
END MODULE radar_data_mie
