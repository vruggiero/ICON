!NEC$ options "-finline-max-depth=4 -finline-max-function-size=500"

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

!!$ TODO Lookup tables:
!!$ - DONE: Update type t_dbzlookuptables: scaling expo f, 1 block of memory + pointers to this block,
!!$   integer indices for different parameters (i_dbz, i_ah, ...), initialization routine
!!$   to allocate space and to initialize the params with neutral values, space and pointers for derivative dP_x/dq_x
!!$ - Add computation of derivative via centered differences over q_i +- 0.05 dq_i
!!$ - Add cubic interpolation instead of linear, also add output of derivative from 3. order polynome
!!$ - Do cubic interpolation for all parameters, in log-log-space for zh, zv, ah, zvh,
!!$   in lin-lin-space for rrhv, irhv, kdp, adp.
!!$ 
!!$ 

MODULE radar_mie_specint

!------------------------------------------------------------------------------
!
! Description: Routines for computation of radar reflectivity and reflectivity-
!              weighted terminal fall speed as well as polarization parameters
!              at single model gridpoints for the interface functions for
!              the COSMO/ICON model based on Mie- or T-matrix theory
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!

  USE radar_kind , ONLY : dp

  USE radar_data , ONLY : miss_value, miss_value_rhv, cmaxlen, &
       lcompute_pe_fwo, icomm_cart_fwo, my_cart_id_fwo, num_compute_fwo

  USE radar_dbzcalc_params_type, ONLY : t_polMP

  USE radar_data_mie, ONLY : &
       rho_w_fwo, rho_ice_fwo, inv_rhow, inv_rhoi, T0C_fwo, Tmin_f, &
       pi_dp, &
       t_dbzlookuptable, t_tabparam, t_tabledef, &
       i_scal_log, i_scal_fscal, i_scal_lin, i_scal_dbz, &
       i_interp_lin, i_interp_cubic, &
       zero_value_lut, nmax_lookup, ldebug_dbz, &
       particle, rain, cloud, snow, ice, graupel, hail, &
       q_crit_radar, n_crit_radar, quasi_zero, Deps,&
       mw_Tmin, mw_Tmax, mi_Tmin, mi_Tmax, &
       t_mgd_params, n_stuetz

  USE radar_mielib_vec, ONLY: &
       m_air, &
       SPHERE_SCATTER_BH_VEC, &
       SPHEROID_SCATTER_TMAT_VEC, &
       MIE_DRYGR_VEC, &
       MIE_DRYSNOW_TWOSPH_VEC, &
       MIE_SOAK_WETGR_VEC, &
       MIE_SOAK_TWOSPH_WETGR_BH_VEC, &
       MIE_WATERSPH_WETGR_VEC, &
       MIE_WATERSPH_WETGR_BH_VEC, &
       MIE_MEAN_WETGR_VEC, &
       MIE_WETSNOW_TWOSPH_VEC, &
       MIE_SPONGY_WETHAIL_BH_VEC, &
       m_complex_water_ray, m_complex_ice_maetzler, &
       m_complex_water_ray_vec, m_complex_ice_maetzler_vec

  USE radar_utilities, ONLY : sub2ind2D, sub2ind3D, round_robin, &
                              tolower, toupper
  
  USE radar_mie_utils, ONLY : &
       nice_mono_1mom, calc_n0_snow, particle_assign, &
       simpson_weights, integweights_mgd, &
       mgd_1mom, mgd_2mom, &
       gamfac_imom_DMGD, D_of_X_average_1M_exp, &
       scale_val_lut, descale_val_lut, scale_dval_lut, &
       init_dbzlookuptable, init_scaling_tabparams

  USE radar_mie_meltdegree, ONLY : &
       degree_of_melting_xxx, degree_of_melting_xxx_Dref

#ifdef NETCDF
  USE netcdf, ONLY :  &
       NF90_NOERR, &
       NF90_CLOBBER, &
       NF90_NOWRITE, &
       nf90_create, &
       nf90_open, &
       nf90_def_dim, &
       nf90_inq_dimid, &
       nf90_inquire_dimension, &
       nf90_inq_varid, &
       nf90_close, &
       nf90_put_att, &
       nf90_strerror, &
       nf90_global, &
       nf90_enddef, &
       nf90_put_var, &
       nf90_def_var, &
       nf90_get_var, &
       nf90_get_att, &
       NF90_NETCDF4, &
       NF90_CLASSIC_MODEL, &
       NF90_REAL, &
       NF90_DOUBLE
#endif

  USE radar_parallel_utilities, ONLY : global_values_radar

!===============================================================================

#ifndef NOMPI
  USE mpi
#endif

!===============================================================================

  IMPLICIT NONE

!==============================================================================

! include statements

#ifdef NOMPI
  INCLUDE "nompi_mpif.h"
#endif

!===============================================================================

  PUBLIC

  ! These PRIVATE functions are not used any more but kept in the code
  !  for reference. There is a more modern vectorized and unified versions of these
  !  procedures in the meantime: zradar_triinterp_lookup_add_vec.
  PRIVATE :: zradar_biinterpolate_lookup, zradar_triinterpolate_lookup, &
             interp2d_lut, interp3d_lut
  
!===============================================================================
!===============================================================================

#ifdef NETCDF
  ! The format for lookup table files is netcdf:
  LOGICAL, PARAMETER :: llut_format_netcdf = .TRUE.
#else
  ! The format for lookup table files is fortran binary:
  LOGICAL, PARAMETER :: llut_format_netcdf = .FALSE.
#endif
  
!===============================================================================
!===============================================================================

CONTAINS

!===============================================================================
!===============================================================================


!==========================================================================================
!
! Calculation of reflectivity for the different hydrometeor categories and microphysics
!  schemes.
!
!==========================================================================================


  !----------------------------------------------------------------------------------------!
  !  Calculate radar reflectivity of rain/cloud droplets using full mie theory             !
  !  as a function of mass density and number density                                      !
  !  assuming a gamma mass distribution (in case of 2-moment-scheme of Seifert and Beheng)
  !  or an exponential size distribution (in case of LM 1-moment-scheme)
  !
  !  The integration of the reflectivity integral is done in diameter-space
  !  by means of Simpson's rule. Diameter space is chosen
  !  since the distribution of equally spaced sampling points is better
  !  suited for the integrand function compared to mass space, in that more sampling
  !  points are within the "small drop" end of the size distribution and less sampling
  !  points at the "large drop" end, compared to mass space.
  !
  !  The procedure contains a minor inconsistency:
  !  To deduce the parameters of the DSD with respect to diameter, it is assumed in accordance
  !  with the Microphysics parameterizations that L (and N in CASE of 2M) are representative
  !  for a DSD spanning over an infinite diameter range (to infinity).
  !  In contrast, when integrating over the spectrum, the integral is taken from Dmin to Dmax,
  !  neglecting the part to infinity.
  !  In fact, the latter is consistent with the physical argument, that no drops larger than
  !  a certain Dmax can exist for longer time because of spontaneous breakup, whereas it is
  !  felt that the assumption of an "infinite" diameter range for L (and N) in microphysics
  !  parameterizations is rather problematic. However, reflectivity is slightly underestimated,
  !  since the mass fraction above Dmax is not redistributed to smaller D, as it should be.
  !  This would be accomplished if DSD-Parameters N0 and lambda would be derived by partial
  !  moments (=incomplete gamma functions)  rather than "normal" gamma functions. However,
  !  that would require an iterative Procedure since the bounds in the incomplete gamma
  !  functions themselves depend on lambda; this is omitted here.
  !----------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------!
  ! Generalized version for both 1- and 2-moment schemes:
  ! (applying generalized (aka modified) gamma distribution)
  !----------------------------------------------------------------------------------------!

  SUBROUTINE zradar_rain_mie_vec(mgd,m_w,lambda,&
                                 luse_tmatrix,ldo_nonsphere,pMP,&
                                 Zfac,Dmin,Dmax,&
                                 zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                                 llookupgen_mode)

    IMPLICIT NONE

    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP
    REAL(KIND=dp), INTENT(in)      :: lambda, Zfac, Dmin, Dmax
    COMPLEX(kind=dp)               :: m_w
    TYPE(t_mgd_params), INTENT(in) :: mgd
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL            :: lmode_lgen

    INTEGER            :: i
    COMPLEX(kind=dp)   :: m_wv(n_stuetz+1), &
                          fa(n_stuetz+1), fb(n_stuetz+1), &
                          fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)      :: dD, &
                          D_w(n_stuetz+1), itw(n_stuetz+1), &
                          aspectratio(n_stuetz+1), sig

    INTEGER,          SAVE :: firstcall = 0, n_stuetz_last
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1), arain(7,n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last, Dmin_last, Dmax_last
    COMPLEX(kind=dp), SAVE :: m_w_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)
    !TYPE(t_polMP),      SAVE :: pMP_last

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF

    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      m_w_last    = CMPLX(-HUGE(1.0_dp),-HUGE(1.0_dp),kind=dp)

!      IF (luse_tmatrix) THEN
!        !CALL calc_orient_water(n_stuetz+1,arain)
!        CALL calc_orient_singlephase(n_stuetz+1,pMP%sig,arain)
!      ELSE
!        arain = 0d0
!      ENDIF

      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_w and, hence
    !           aspectratio, as SAVEd variables and recalc them only if Dmin or
    !           Dmax have changed (or n_stuetz, but that is not an input param
    !           but hardcoded herein, hence not changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_w(i) = Dmin + DBLE(i-1) * dD
    END DO

    m_wv = m_w

    IF (luse_tmatrix) THEN
      aspectratio = aspectratio_singlephase(D_w,n_stuetz+1,'rain',rain,pMP,sig)
      IF (.NOT. ldo_nonsphere) aspectratio = 1d0
      CALL calc_orient_singlephase(n_stuetz+1,sig,arain)
    ELSE
      aspectratio = 1d0
      arain = 0d0
    ENDIF

!    IF (firstcall /= 1) THEN
!      lambda_last   = -HUGE(1.0_dp)
!      m_w_last      = CMPLX(-HUGE(1.0_dp),-HUGE(1.0_dp),kind=dp)
!      Dmin_last     = -HUGE(1.0_dp)
!      Dmax_last     = -HUGE(1.0_dp)
!      n_stuetz_last = -HUGE(1)
!      !pMP_last%ARmodel = 'XXX'
!      !pMP_last%c0      = -HUGE(1.0_dp)
!      !pMP_last%c1      = -HUGE(1.0_dp)
!      !pMP_last%ARmin   = -HUGE(1.0_dp)
!      !pMP_last%sig     = -HUGE(1.0_dp)
!
!      firstcall = 1
!    ENDIF
!
!    IF (n_stuetz_last /= n_stuetz) THEN
!      simpson = simpson_weights(n_stuetz)
!    END IF
!
!    IF (n_stuetz_last /= n_stuetz .OR. &
!        Dmin_last /= Dmin .OR. Dmax_last /= Dmax) THEN
!      dD   = (Dmax-Dmin)/n_stuetz
!      DO i=1,n_stuetz+1
!        D_w(i) = Dmin + DBLE(i-1) * dD
!      END DO
!    END IF
!
!    IF (luse_tmatrix) THEN
!!      IF (n_stuetz_last /= n_stuetz .OR. &
!!          Dmin_last /= Dmin .OR. Dmax_last /= Dmax .OR. &
!!          pMP_last%ARmodel /= pMP%ARmodel .OR. &
!!          pMP_last%c0 /= pMP%c0 .OR. pMP_last%c1 /= pMP%c1 .OR.
!!          pMP_last%ARmin /= pMP%ARmin .OR. pMP_last%sig /= pMP%sig) THEN
!        ! aspectratio_singlephase adjusts pMP%sig in case of "full-MP model" like R11.
!        !   hence, needs to be called before calc_orient_* routines
!        aspectratio = aspectratio_singlephase(D_w,n_stuetz+1,'rain',rain,pMP)
!        CALL calc_orient_singlephase(n_stuetz+1,pMP%sig,arain)
!!      ENDIF
!      IF (.NOT. ldo_nonsphere) THEN
!        aspectratio = 1d0
!      ENDIF
!    ELSE
!      aspectratio = 1d0
!      arain = 0d0
!    ENDIF
!
!    n_stuetz_last = n_stuetz
!    Dmin_last = Dmin
!    Dmax_last = Dmax
!!    pMP_last = pMP
!
!    m_wv = m_w

    IF (.NOT.lmode_lgen .OR. lambda_last /= lambda .OR. m_w_last /= m_w) THEN
    
      IF (luse_tmatrix) THEN
        CALL SPHEROID_SCATTER_TMAT_VEC(MAX(D_w,Deps),m_wv,DBLE(m_air),aspectratio,lambda,&
                                       fa,fb,fa0,fb0,n_stuetz+1)
      ELSE
        CALL SPHERE_SCATTER_BH_VEC(MAX(D_w,Deps),m_wv,DBLE(m_air),lambda,&
                                   fa,fb,fa0,fb0,n_stuetz+1)
      ENDIF

      IF (lmode_lgen) THEN
        lambda_last = lambda
        m_w_last    = m_w
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    itw = integweights_mgd(mgd,D_w,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,arain,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_rain_mie_vec

  !----------------------------------------------------------------------------------------!
  !  Calculate radar reflectivity of dry graupel/hail using full mie theory                !
  !  as a function of mass density and number density                                      !
  !  assuming a gamma distribution (in case of 2-moment-scheme of Seifert and Beheng) or an
  !  exponential size distribution (in case of LM 1-moment-scheme)
  !
  !  The integration of the reflectivity integral is done in diameter-space
  !  by means of Simpson's rule. Diameter space is chosen
  !  since the distribution of equally spaced sampling points is better
  !  suited for the integrand function compared to mass space, in that more sampling
  !  points are within the "small drop" end of the size distribution and less sampling
  !  points at the "large drop" end.
  !
  !  The procedure contains a minor inconsistency:
  !  To deduce the parameters of the DSD with respect to diameter,
  !  it is assumed in accordance with the Microphysics parameterizations
  !  that L (and N in CASE of 2M) are representative
  !  for a DSD spanning over an infinite diameter range (to infinity).
  !  In contrast, when integrating over the spectrum, the integral is taken
  !  from Dmin to Dmax, neglecting the part to infinity.
  !  In fact, the latter is consistent with the physical argument, that no drops larger
  !  than a Dmax can exist for longer time because of spontaneous breakup, whereas it is
  !  felt that the assumption of an "infinite" diameter range for L (and N) in
  !  microphysics parameterizations is rather problematic.
  !
  !----------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------!
  ! Generalized version for both 1- and 2-moment schemes:
  ! (applying generalized (aka modified) gamma distribution)
  !----------------------------------------------------------------------------------------!

  ! Dry cloud ice:
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_ice_mie_vec(&
      mgd,m_i,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,&
      Zfac,parti,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring,matrixstring,inclusionstring,&
      llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)     :: parti
    TYPE(t_mgd_params), INTENT(in) :: mgd
    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP
    REAL(KIND=dp), INTENT(in)      :: lambda ,Zfac, Dmin, Dmax
    COMPLEX(kind=dp), INTENT(in)   :: m_i
    CHARACTER(len=*), INTENT(in)   :: mixingrulestring, matrixstring, inclusionstring
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL            :: lmode_lgen

    INTEGER            :: i
    COMPLEX(kind=dp)   :: m_iv(n_stuetz+1), &
                          fa(n_stuetz+1), fb(n_stuetz+1), &
                          fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)      :: dD, &
                          D_i(n_stuetz+1), itw(n_stuetz+1), &
                          aspectratio(n_stuetz+1), sig

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1), aice(7,n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last
    COMPLEX(kind=dp), SAVE :: m_i_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      m_i_last    = CMPLX(-HUGE(1.0_dp),-HUGE(1.0_dp),kind=dp)

!      IF (luse_tmatrix) THEN
!        fmelt=0.0d0
!        CALL calc_orient_singlephase(n_stuetz+1,sig,aice)
!      ELSE
!        aice = 0d0
!      ENDIF

      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_i and, hence
    !           aspectratio, as SAVEd variables and recalc them only if Dmin or
    !           Dmax have changed (or n_stuetz, but that is not an input param
    !           but hardcoded herein, hence not changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_i(i)   = Dmin + DBLE(i-1) * dD
    ENDDO

    m_iv = m_i

    IF (luse_tmatrix) THEN
      aspectratio = aspectratio_singlephase(D_i,n_stuetz+1,'ice',ice,pMP,sig)
      IF (.NOT. ldo_nonsphere) aspectratio = 1d0
      CALL calc_orient_singlephase(n_stuetz+1,sig,aice)
    ELSE
      aspectratio = 1d0
      aice = 0d0
    ENDIF

    IF (.NOT.lmode_lgen .OR. lambda_last /= lambda .OR. m_i_last /= m_i) THEN
    
      CALL MIE_DRYGR_VEC(&
           D_i,parti%a_geo,parti%b_geo,&
           m_iv,aspectratio,lambda,&
           fa,fb,fa0,fb0, &
           mixingrulestring,matrixstring,inclusionstring,&
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        m_i_last    = m_i
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    itw = integweights_mgd(mgd,D_i,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,aice,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_ice_mie_vec
  ! Dry snow: Two-sphere-model with radius ratio = r_inner/r_total
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_snow_mie_vec(&
      mgd,m_i,Rratio,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,&
      Zfac,parti,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring_core,matrixstring_core,inclusionstring_core,&
      mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)     :: parti
    TYPE(t_mgd_params), INTENT(in) :: mgd
    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP
    COMPLEX(kind=dp), INTENT(in)   :: m_i
    REAL(KIND=dp), INTENT(in)      :: lambda, Zfac, Rratio, Dmin, Dmax
    CHARACTER(len=*), INTENT(in)   :: mixingrulestring_core, matrixstring_core, inclusionstring_core, &
                                      mixingrulestring_shell, matrixstring_shell, inclusionstring_shell
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL            :: lmode_lgen

    INTEGER            :: i
    COMPLEX(kind=dp)   :: m_iv(n_stuetz+1), &
                          fa(n_stuetz+1), fb(n_stuetz+1), &
                          fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)      :: dD, &
                          D_s(n_stuetz+1), itw(n_stuetz+1), &
                          aspectratio(n_stuetz+1), sig

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1), asnow(7,n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last, Rratio_last
    COMPLEX(kind=dp), SAVE :: m_i_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      Rratio_last = -HUGE(1.0_dp)
      m_i_last    = CMPLX(-HUGE(1.0_dp),-HUGE(1.0_dp),kind=dp)

!      IF (luse_tmatrix) THEN
!        fmelt=0.0d0
!        CALL calc_orient_singlephase(n_stuetz+1,sig,asnow)
!      ELSE
!        asnow = 0d0
!      ENDIF

      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_s and, hence
    !           aspectratio, as SAVEd variables and recalc them only if Dmin or
    !           Dmax have changed (or n_stuetz, but that is not an input param
    !           but hardcoded herein, hence not changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_s(i)   = Dmin + DBLE(i-1) * dD
    ENDDO

    m_iv = m_i

    IF (luse_tmatrix) THEN
      aspectratio = aspectratio_singlephase(D_s,n_stuetz+1,'snow',snow,pMP,sig)
      IF (.NOT. ldo_nonsphere) aspectratio = 1d0
      CALL calc_orient_singlephase(n_stuetz+1,sig,asnow)
    ELSE
      aspectratio = 1d0
      asnow = 0d0
    ENDIF

    IF (.NOT.lmode_lgen .OR. lambda_last /= lambda .OR. m_i_last /= m_i .OR. &
                             Rratio_last /= Rratio ) THEN
    
      CALL MIE_DRYSNOW_TWOSPH_VEC(&
           D_s,parti%a_geo,parti%b_geo,&
           m_iv,Rratio,aspectratio,lambda,&
           fa,fb,fa0,fb0, &
           mixingrulestring_core,matrixstring_core,inclusionstring_core,&
           mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        Rratio_last = Rratio
        m_i_last    = m_i
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    itw = integweights_mgd(mgd,D_s,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,asnow,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_snow_mie_vec

  ! Dry graupel:
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_graupel_mie_vec(&
      mgd,m_i,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,&
      Zfac,parti,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring,matrixstring,inclusionstring,&
      llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)     :: parti
    TYPE(t_mgd_params), INTENT(in) :: mgd
    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP
    REAL(KIND=dp), INTENT(in)      :: lambda, Zfac, Dmin, Dmax
    COMPLEX(kind=dp), INTENT(in)   :: m_i
    CHARACTER(len=*), INTENT(in)   :: mixingrulestring, matrixstring, inclusionstring
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL            :: lmode_lgen

    INTEGER            :: i
    COMPLEX(kind=dp)   :: m_iv(n_stuetz+1), &
                          fa(n_stuetz+1), fb(n_stuetz+1), &
                          fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)      :: dD, &
                          D_g(n_stuetz+1), itw(n_stuetz+1), &
                          aspectratio(n_stuetz+1), sig

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1), agraupel(7,n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last
    COMPLEX(kind=dp), SAVE :: m_i_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      m_i_last    = CMPLX(-HUGE(1.0_dp),-HUGE(1.0_dp),kind=dp)

!      IF (luse_tmatrix) THEN
!        fmelt=0.0d0
!        CALL calc_orient_singlephase(n_stuetz+1,sig,agraupel)
!      ELSE
!        agraupel = 0d0
!      ENDIF

      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_g and, hence
    !           aspectratio, as SAVEd variables and recalc them only if Dmin or
    !           Dmax have changed (or n_stuetz, but that is not an input param
    !           but hardcoded herein, hence not changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_g(i)   = Dmin + DBLE(i-1) * dD
    ENDDO

    m_iv = m_i

    IF (luse_tmatrix) THEN
      aspectratio = aspectratio_singlephase(D_g,n_stuetz+1,'graupel',graupel,pMP,sig)
      IF (.NOT. ldo_nonsphere) aspectratio = 1d0
      CALL calc_orient_singlephase(n_stuetz+1,sig,agraupel)
    ELSE
      aspectratio = 1d0
      agraupel = 0d0
    ENDIF

    IF (.NOT.lmode_lgen .OR. lambda_last /= lambda .OR. m_i_last /= m_i) THEN
    
      CALL MIE_DRYGR_VEC(&
           D_g,parti%a_geo,parti%b_geo,&
           m_iv,aspectratio,lambda,&
           fa,fb,fa0,fb0, &
           mixingrulestring,matrixstring,inclusionstring,&
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        m_i_last    = m_i
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    itw = integweights_mgd(mgd,D_g,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,agraupel,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_graupel_mie_vec

  ! Dry hail:
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_hail_mie_vec(&
      mgd,m_i,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,&
      Zfac,parti,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring,matrixstring,inclusionstring,&
      llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)     :: parti
    TYPE(t_mgd_params), INTENT(in) :: mgd
    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP
    REAL(KIND=dp), INTENT(in)      :: lambda, Zfac, Dmin, Dmax
    COMPLEX(kind=dp), INTENT(in)   :: m_i
    CHARACTER(len=*), INTENT(in)   :: mixingrulestring, matrixstring, inclusionstring
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL            :: lmode_lgen

    INTEGER            :: i
    COMPLEX(kind=dp)   :: m_iv(n_stuetz+1), &
                          fa(n_stuetz+1), fb(n_stuetz+1), &
                          fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)      :: dD, &
                          D_g(n_stuetz+1), itw(n_stuetz+1), &
                          aspectratio(n_stuetz+1), sig

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1), ahail(7,n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last
    COMPLEX(kind=dp), SAVE :: m_i_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson     = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      m_i_last    = CMPLX(-HUGE(1.0_dp),-HUGE(1.0_dp),kind=dp)

!      IF (luse_tmatrix) THEN
!        fmelt=0.0d0
!        CALL calc_orient_singlephase(n_stuetz+1,sig,ahail)
!      ELSE
!        ahail = 0d0
!      ENDIF

      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_g and, hence
    !           aspectratio, as SAVEd variables and recalc them only if Dmin or
    !           Dmax have changed (or n_stuetz, but that is not an input param
    !           but hardcoded herein, hence not changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_g(i)   = Dmin + DBLE(i-1) * dD
    ENDDO

    m_iv = m_i

    IF (luse_tmatrix) THEN
      aspectratio = aspectratio_singlephase(D_g,n_stuetz+1,'hail',hail,pMP,sig)
      IF (.NOT. ldo_nonsphere) aspectratio = 1d0
      CALL calc_orient_singlephase(n_stuetz+1,sig,ahail)
    ELSE
      aspectratio = 1d0
      ahail = 0d0
    ENDIF

    IF (.NOT.lmode_lgen .OR. lambda_last /= lambda .OR. m_i_last /= m_i) THEN
    
      CALL MIE_DRYGR_VEC(&
           D_g,parti%a_geo,parti%b_geo,&
           m_iv,aspectratio,lambda,&
           fa,fb,fa0,fb0,&
           mixingrulestring,matrixstring,inclusionstring,&
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        m_i_last    = m_i
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    ! Size distribution integration weights
    itw = integweights_mgd(mgd,D_g,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,ahail,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_hail_mie_vec

  !----------------------------------------------------------------------------------------!
  !  Calculate radar reflectivity of wet graupel/hail using full mie theory                !
  !  as a function of mass density and number density                                      !
  !
  !----------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------!
  ! Generalized version for both 1- and 2-moment schemes:
  ! (applying generalized (aka modified) gamma distribution)
  !----------------------------------------------------------------------------------------!

  ! Melting cloud ice: soaked version
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_wetice_mie_vec(&
      mgd,T_a,m_i,m_w,&
      itype_Dref_fmelt,Tmeltbegin,meltdegTmin,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
      Zfac,parti,rain,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring,matrixstring,inclusionstring,&
      hoststring,hostmatrixstring,hostinclusionstring,Tmax,&
      llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)     :: parti, rain
    TYPE(t_mgd_params), INTENT(in) :: mgd
    INTEGER, INTENT(in)            :: itype_Dref_fmelt
    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP, pMPr
    COMPLEX(kind=dp), INTENT(in)   :: m_i, m_w
    REAL(KIND=dp), INTENT(in)      :: T_a, lambda, Zfac, Dmin, Dmax,&
                                      Tmeltbegin, meltdegTmin ,Tmax
    CHARACTER(len=*), INTENT(in)   :: mixingrulestring, matrixstring, inclusionstring,&
                                      hoststring, hostmatrixstring, hostinclusionstring
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL                  :: lmode_lgen

    ! FIXME UB: needs to be looked at! So far only slightly modified from snow fmelt:
    ! .. Parameters for the degree of melting parameterization with fixed Dref:
    !   ... for itype_Dref_fmelt = 1, the old method with Dref=fct(PSD):
    REAL(KIND=dp), PARAMETER :: Dref_max = 5e-4_dp,   Dexpo_type1 = 0.5_dp
    !   ... for itype_Dref_fmelt = 2, Dref and Dexpo are tuned to approximately reproduce the old method:
    REAL(KIND=dp), PARAMETER :: Dref_type2 = 5e-5_dp, Dexpo_type2 = 0.3_dp

    INTEGER                  :: i
    COMPLEX(kind=dp)         :: m_iv(n_stuetz+1), m_wv(n_stuetz+1), &
                                fa(n_stuetz+1), fb(n_stuetz+1), &
                                fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)            :: dD, meltratio_outside, &
                                D_r(n_stuetz+1), D_i(n_stuetz+1), &
                                itw(n_stuetz+1), fmelt(n_stuetz+1), &
                                aspectratio(n_stuetz+1), sig, sigr, &
                                arrain(n_stuetz+1), ardry(n_stuetz+1), awetice(7,n_stuetz+1)

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last, Tmax_last, T_a_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      Tmax_last   = -HUGE(1.0_dp)
      T_a_last    = -HUGE(1.0_dp)
      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_i as SAVEd variables
    !           and recalc them only if Dmin or Dmax have changed (or n_stuetz,
    !           but that is not an input param but hardcoded herein, hence not
    !           changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_i(i) = Dmin + DBLE(i-1) * dD
    END DO

    m_iv = m_i
    m_wv = m_w

    ! Parameterisation of degree of melting
    !
    SELECT CASE (itype_Dref_fmelt)
    CASE (1)
      ! Old method with Dref=fct(PSD mean size):
      CALL degree_of_melting_xxx(mgd,T_a,D_i,n_stuetz+1,Dexpo_type1,Dref_max,parti,&
                                 Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    CASE (2)
      ! New method with fixed Dref = 5e-4 m instead of Dref from PSD, to save computing time below:
      CALL degree_of_melting_xxx_Dref(T_a,D_i,n_stuetz+1,Dexpo_type2,Dref_type2,&
                                      Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    END SELECT

    IF (luse_tmatrix) THEN
      ardry = aspectratio_singlephase(D_i,n_stuetz+1,'ice',ice,pMP,sig)

      D_r = ((parti%a_geo * D_i**parti%b_geo)/rain%a_geo)**(1d0/rain%b_geo)
      arrain = aspectratio_singlephase(D_r,n_stuetz+1,'rain',rain,pMPr,sigr)

      IF (.NOT. ldo_nonsphere) THEN
        aspectratio = 1d0
        ardry = 1d0
      ELSE
        aspectratio = ardry + fmelt*(arrain - ardry)
      END IF

      CALL calc_orient_mixedphase(n_stuetz+1,fmelt,sig,sigr,awetice)
    ELSE
      aspectratio = 1d0
      ardry = 1d0
      awetice = 0d0
    ENDIF

    IF ( .NOT.lmode_lgen .OR. itype_Dref_fmelt /= 2 .OR. lambda_last /= lambda .OR. &
         T_a_last /= T_a .OR. Tmax_last /= Tmax ) THEN
    
      ! .. Compute the scattering amplitudes for all sizes in D_S-vector only if something relevant has changed:

      ! Ratio of melting mass on surface to total melted mass.
      ! Complementary part melts in the inner of the particle. This ration increases
      ! linearly from given value to 1.0 with increasing melting degree.
      meltratio_outside = 0.7d0

      ! Alternatives: MIE_WATERSPH_WETGR_VEC(), MIE_SOAK_TWOSPH_WETGR_VEC(),
      ! MIE_MEAN_WETGR_VEC(), MIE_WETSNOW_TWOSPH_VEC()
      CALL MIE_SOAK_WETGR_VEC(&
           D_i,parti%a_geo,parti%b_geo,fmelt,meltratio_outside,&
           m_wv,m_iv,aspectratio,ardry,lambda,&
           fa,fb,fa0,fb0,&
           mixingrulestring,matrixstring,inclusionstring,&
           hoststring,hostmatrixstring,hostinclusionstring,&
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        Tmax_last   = Tmax
        T_a_last    = T_a
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    itw = integweights_mgd(mgd,D_i,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,awetice,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_wetice_mie_vec

  ! Wet/melting snow
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_wetsnow_mie_vec(&
      mgd,T_a,m_i,m_w,Rratio,itype_Dref_fmelt,Tmeltbegin,meltdegTmin,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
      Zfac,parti,rain,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      mixingrulestring_core,matrixstring_core,inclusionstring_core,&
      hoststring_shell,hostmatrixstring_shell,hostinclusionstring_shell,&
      hoststring_core,hostmatrixstring_core,hostinclusionstring_core,Tmax,&
      llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)     :: parti, rain
    TYPE(t_mgd_params), INTENT(in) :: mgd
    INTEGER, INTENT(in)            :: itype_Dref_fmelt
    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP, pMPr
    COMPLEX(kind=dp), INTENT(in)   :: m_i, m_w
    REAL(KIND=dp), INTENT(in)      :: T_a, lambda, Rratio, Dmin, Dmax, Zfac, &
                                      Tmeltbegin, meltdegTmin,Tmax
    CHARACTER(len=*), INTENT(in)   :: mixingrulestring_core, matrixstring_core, inclusionstring_core,&
                                      mixingrulestring_shell, matrixstring_shell, inclusionstring_shell,&
                                      hoststring_shell, hostmatrixstring_shell, hostinclusionstring_shell,&
                                      hoststring_core, hostmatrixstring_core, hostinclusionstring_core
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL                  :: lmode_lgen

    ! .. Parameters for the degree of melting parameterization with fixed Dref:
    !   ... for itype_Dref_fmelt = 1, the old method with Dref=fct(PSD):
    REAL(KIND=dp), PARAMETER :: Dref_max = 6e-3_dp,   Dexpo_type1 = 0.5_dp
    !   ... for itype_Dref_fmelt = 2, Dref and Dexpo are tuned to approximately reproduce the old method:
    REAL(KIND=dp), PARAMETER :: Dref_type2 = 1e-6_dp, Dexpo_type2 = 0.1_dp

    INTEGER                  :: i
    COMPLEX(kind=dp)         :: m_iv(n_stuetz+1), m_wv(n_stuetz+1), &
                                fa(n_stuetz+1), fb(n_stuetz+1), &
                                fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)            :: dD, meltratio_outside, &
                                ardry(n_stuetz+1), arrain(n_stuetz+1), D_r(n_stuetz+1), &
                                D_s(n_stuetz+1), itw(n_stuetz+1), fmelt(n_stuetz+1), &
                                aspectratio(n_stuetz+1), sig, sigr, &
                                awetsnow(7,n_stuetz+1)

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last, Tmax_last, T_a_last, Rratio_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      Tmax_last   = -HUGE(1.0_dp)
      T_a_last    = -HUGE(1.0_dp)
      Rratio_last = -HUGE(1.0_dp)
      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_s as SAVEd variables
    !           and recalc them only if Dmin or Dmax have changed (or n_stuetz,
    !           but that is not an input param but hardcoded herein, hence not
    !           changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_s(i) = Dmin + DBLE(i-1) * dD
    END DO


    m_iv = m_i
    m_wv = m_w

    ! Parameterisation of degree of melting
    !
    SELECT CASE (itype_Dref_fmelt)
    CASE (1)
      ! Old method with Dref=fct(PSD mean size):
      CALL degree_of_melting_xxx(mgd,T_a,D_s,n_stuetz+1,Dexpo_type1,Dref_max,parti,&
                                 Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    CASE (2)
      ! New method with fixed Dref = 1e-6 m instead of Dref from PSD, to save computing time below:
      CALL degree_of_melting_xxx_Dref(T_a,D_s,n_stuetz+1,Dexpo_type2,Dref_type2,&
                                      Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    END SELECT
    
    IF (luse_tmatrix) THEN
      ardry = aspectratio_singlephase(D_s,n_stuetz+1,'snow',snow,pMP,sig)

      D_r = ((parti%a_geo * D_s**parti%b_geo)/rain%a_geo)**(1d0/rain%b_geo)
      arrain = aspectratio_singlephase(D_r,n_stuetz+1,'rain',rain,pMPr,sigr)

      IF (.NOT. ldo_nonsphere) THEN
        aspectratio = 1d0
        ardry = 1d0
      ELSE
        aspectratio = ardry + fmelt*(arrain - ardry)
      END IF

      CALL calc_orient_mixedphase(n_stuetz+1,fmelt,sig,sigr,awetsnow)
    ELSE
      aspectratio = 1d0
      ardry = 1d0
      awetsnow = 0d0
    END IF

    IF ( .NOT.lmode_lgen .OR. itype_Dref_fmelt /= 2 .OR. lambda_last /= lambda .OR. &
         T_a_last /= T_a .OR. Tmax_last /= Tmax .OR. Rratio_last /= Rratio ) THEN
    
      ! .. Compute the scattering amplitudes for all sizes in D_S-vector only if something relevant has changed:

      ! Ratio of melting mass on surface of core and shell to total melted mass of each.
      ! Identical for both core and shell.
      ! Complementary part melts in the inner of core and shell, respectively. This
      ! ration increases linearly from given value to 1.0 with increasing melting degree.
      meltratio_outside = 0.7d0

      ! Alternatives: MIE_WATERSPH_WETGR_VEC(), MIE_SOAK_TWOSPH_WETGR_VEC(),
      !               MIE_SOAK_WETGR_VEC(), MIE_MEAN_WETGR_VEC()
      CALL MIE_WETSNOW_TWOSPH_VEC(&
           D_s,parti%a_geo,parti%b_geo,fmelt,meltratio_outside,m_wv,m_iv,&
           aspectratio,ardry,lambda,Rratio,&
           fa,fb,fa0,fb0,&
           mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
           mixingrulestring_core,matrixstring_core,inclusionstring_core,&
           hoststring_shell,hostmatrixstring_shell,hostinclusionstring_shell,&
           hoststring_core,hostmatrixstring_core,hostinclusionstring_core,&
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        Tmax_last   = Tmax
        T_a_last    = T_a
        Rratio_last = Rratio
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF
    
    itw = integweights_mgd(mgd,D_s,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,awetsnow,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_wetsnow_mie_vec

  ! Melting graupel: soaked version
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_wetgr_mie_vec(&
      mgd,T_a,m_i,m_w,&
      itype_Dref_fmelt,Tmeltbegin,meltdegTmin,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
      Zfac,parti,rain,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring,matrixstring,inclusionstring,&
      hoststring,hostmatrixstring,hostinclusionstring,Tmax,&
      llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)     :: parti, rain
    TYPE(t_mgd_params), INTENT(in) :: mgd
    INTEGER, INTENT(in)            :: itype_Dref_fmelt
    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP, pMPr
    COMPLEX(kind=dp), INTENT(in)   :: m_i, m_w
    REAL(KIND=dp), INTENT(in)      :: T_a,lambda, Zfac, Dmin, Dmax,&
                                      Tmeltbegin, meltdegTmin,Tmax
    CHARACTER(len=*), INTENT(in)   :: mixingrulestring, matrixstring, inclusionstring,&
                                      hoststring, hostmatrixstring, hostinclusionstring
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL                  :: lmode_lgen

    ! .. Parameters for the degree of melting parameterization with fixed Dref:
    !   ... for itype_Dref_fmelt = 1, the old method with Dref=fct(PSD):
    REAL(KIND=dp), PARAMETER :: Dref_max = 6e-3_dp,   Dexpo_type1 = 0.6_dp
    !   ... for itype_Dref_fmelt = 2, Dref and Dexpo are tuned to approximately reproduce the old method:
    REAL(KIND=dp), PARAMETER :: Dref_type2 = 5e-4_dp, Dexpo_type2 = 0.5_dp

    INTEGER                  :: i
    COMPLEX(kind=dp)         :: m_iv(n_stuetz+1), m_wv(n_stuetz+1), &
                                fa(n_stuetz+1), fb(n_stuetz+1), &
                                fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)            :: dD, meltratio_outside, &
                                D_r(n_stuetz+1), D_g(n_stuetz+1), &
                                itw(n_stuetz+1), fmelt(n_stuetz+1), &
                                aspectratio(n_stuetz+1), sig, sigr, &
                                arrain(n_stuetz+1), ardry(n_stuetz+1), awetgr(7,n_stuetz+1)

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last, Tmax_last, T_a_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      Tmax_last   = -HUGE(1.0_dp)
      T_a_last    = -HUGE(1.0_dp)
      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_g as SAVEd variables
    !           and recalc them only if Dmin or Dmax have changed (or n_stuetz,
    !           but that is not an input param but hardcoded herein, hence not
    !           changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_g(i) = Dmin + DBLE(i-1) * dD
    END DO

    m_iv = m_i
    m_wv = m_w

    ! Parameterisation of degree of melting
    !
    SELECT CASE (itype_Dref_fmelt)
    CASE (1)
      ! Old method with Dref=fct(PSD mean size):
      CALL degree_of_melting_xxx(mgd,T_a,D_g,n_stuetz+1,Dexpo_type1,Dref_max,parti,&
                                 Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    CASE (2)
      ! New method with fixed Dref = 5e-4 m instead of Dref from PSD, to save computing time below:
      CALL degree_of_melting_xxx_Dref(T_a,D_g,n_stuetz+1,Dexpo_type2,Dref_type2,&
                                      Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    END SELECT

    IF (luse_tmatrix ) THEN
      ardry = aspectratio_singlephase(D_g,n_stuetz+1,'graupel',graupel,pMP,sig)

      D_r = ((parti%a_geo * D_g**parti%b_geo)/rain%a_geo)**(1d0/rain%b_geo)
      arrain = aspectratio_singlephase(D_r,n_stuetz+1,'rain',rain,pMPr,sigr)

      IF (.NOT. ldo_nonsphere) THEN
        aspectratio = 1d0
        ardry = 1d0
      ELSE
        ! JM200515: Would be nice to use a WHERE construct here, too. But I'm not
        !           sure how nested WHEREs (here fmelt-limits on a first level,
        !           D_g limits on a second (for fmelt<0.2).

        ! Or, we do two separate WHERE constructs instead...
        WHERE (fmelt < 0.2d0)
          aspectratio = ardry-5.0d0*(ardry - 0.8d0)*fmelt
        ELSEWHERE (fmelt < 0.8)
          aspectratio = 0.88d0 - 0.40d0*fmelt
        ELSEWHERE
          aspectratio = 2.8d0 - 4.0d0*arrain+5.0d0*(arrain-0.56d0)*fmelt
        END WHERE
      END IF

      CALL calc_orient_mixedphase(n_stuetz+1,fmelt,sig,sigr,awetgr)
    ELSE
      aspectratio = 1d0
      ardry = 1d0
      awetgr = 0d0
    ENDIF

    IF ( .NOT.lmode_lgen .OR. itype_Dref_fmelt /= 2 .OR. lambda_last /= lambda .OR. &
         T_a_last /= T_a .OR. Tmax_last /= Tmax ) THEN
    
      ! .. Compute the scattering amplitudes for all sizes in D_S-vector only if something relevant has changed:

      ! Ratio of melting mass on surface to total melted mass.
      ! Complementary part melts in the inner of the particle. This ration increases
      ! linearly from given value to 1.0 with increasing melting degree.
      meltratio_outside = 0.7d0

      ! Alternatives: MIE_WATERSPH_WETGR_VEC(), MIE_SOAK_TWOSPH_WETGR_VEC(),
      ! MIE_MEAN_WETGR_VEC(), MIE_WETSNOW_TWOSPH_VEC()
      CALL MIE_SOAK_WETGR_VEC(&
           D_g,parti%a_geo,parti%b_geo,fmelt,meltratio_outside,&
           m_wv,m_iv,aspectratio,ardry,lambda,&
           fa,fb,fa0,fb0,&
           mixingrulestring,matrixstring,inclusionstring,&
           hoststring,hostmatrixstring,hostinclusionstring,&
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        Tmax_last   = Tmax
        T_a_last    = T_a
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    itw = integweights_mgd(mgd,D_g,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,awetgr,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_wetgr_mie_vec

  ! Melting graupel: soaked two-sphere version
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_wetgr_twosph_mie_vec(&
      mgd,T_a,m_i,m_w,&
      itype_Dref_fmelt,Tmeltbegin,meltdegTmin,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
      Zfac,parti,rain,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring_iceair, matrixstring_iceair, inclusionstring_iceair,&
      mixingrulestring_icewater, matrixstring_icewater, inclusionstring_icewater,&
      Tmax,llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)     :: parti, rain
    TYPE(t_mgd_params), INTENT(in) :: mgd
    INTEGER, INTENT(in)            :: itype_Dref_fmelt
    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP, pMPr
    COMPLEX(kind=dp), INTENT(in)   :: m_i, m_w
    REAL(KIND=dp), INTENT(in)      :: T_a,lambda, Zfac, Dmin, Dmax,&
                                      Tmeltbegin, meltdegTmin, Tmax
    CHARACTER(len=*), INTENT(in)   :: mixingrulestring_iceair, matrixstring_iceair,&
                                      inclusionstring_iceair, &
                                      mixingrulestring_icewater, matrixstring_icewater,&
                                      inclusionstring_icewater
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL                  :: lmode_lgen

    ! .. Parameters for the degree of melting parameterization with fixed Dref:
    !   ... for itype_Dref_fmelt = 1, the old method with Dref=fct(PSD):
    REAL(KIND=dp), PARAMETER :: Dref_max = 6e-3_dp,   Dexpo_type1 = 0.6_dp
    !   ... for itype_Dref_fmelt = 2, Dref and Dexpo are tuned to approximately reproduce the old method:
    REAL(KIND=dp), PARAMETER :: Dref_type2 = 5e-4_dp, Dexpo_type2 = 0.5_dp

    INTEGER                  :: i
    COMPLEX(kind=dp)         :: m_iv(n_stuetz+1), m_wv(n_stuetz+1), &
                                fa(n_stuetz+1), fb(n_stuetz+1), &
                                fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)            :: dD, &
                                D_r(n_stuetz+1), D_g(n_stuetz+1), &
                                itw(n_stuetz+1),fmelt(n_stuetz+1), &
                                aspectratio(n_stuetz+1), sig, sigr, &
                                arrain(n_stuetz+1),ardry(n_stuetz+1), awetgr(7,n_stuetz+1)

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last, Tmax_last, T_a_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      Tmax_last   = -HUGE(1.0_dp)
      T_a_last    = -HUGE(1.0_dp)
      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_g as SAVEd variables
    !           and recalc them only if Dmin or Dmax have changed (or n_stuetz,
    !           but that is not an input param but hardcoded herein, hence not
    !           changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_g(i) = Dmin + DBLE(i-1) * dD
    END DO

    m_iv = m_i
    m_wv = m_w

    ! Parameterisation of degree of melting
    !
    SELECT CASE (itype_Dref_fmelt)
    CASE (1)
      ! Old method with Dref=fct(PSD mean size):
      CALL degree_of_melting_xxx(mgd,T_a,D_g,n_stuetz+1,Dexpo_type1,Dref_max,parti,&
                                 Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    CASE (2)
      ! New method with fixed Dref = 5e-4 m instead of Dref from PSD, to save computing time below:
      CALL degree_of_melting_xxx_Dref(T_a,D_g,n_stuetz+1,Dexpo_type2,Dref_type2,&
                                      Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    END SELECT

    IF (luse_tmatrix) THEN
      ardry = aspectratio_singlephase(D_g,n_stuetz+1,'graupel',graupel,pMP,sig)

      D_r = ((parti%a_geo * D_g**parti%b_geo)/rain%a_geo)**(1d0/rain%b_geo)
      arrain = aspectratio_singlephase(D_r,n_stuetz+1,'rain',rain,pMPr,sigr)

      IF (.NOT. ldo_nonsphere) THEN
        aspectratio = 1d0
        ardry = 1d0
      ELSE
        ! JM200515: Would be nice to use a WHERE construct here, too. But I'm not
        !           sure how nested WHEREs (here fmelt-limits on a first level,
        !           D_g limits on a second (for fmelt<0.2).

        ! Or, we do two separate WHERE constructs instead...
        WHERE (fmelt < 0.2d0)
          aspectratio = ardry-5.0d0*(ardry - 0.8d0)*fmelt
        ELSEWHERE (fmelt < 0.8)
          aspectratio = 0.88d0 - 0.40d0*fmelt
        ELSEWHERE
          aspectratio = 2.8d0 - 4.0d0*arrain+5.0d0*(arrain-0.56d0)*fmelt
        END WHERE
      END IF

      CALL calc_orient_mixedphase(n_stuetz+1,fmelt,sig,sigr,awetgr)
    ELSE
      aspectratio = 1d0
      ardry = 1d0
      awetgr = 0d0
    ENDIF

    IF ( .NOT.lmode_lgen .OR. itype_Dref_fmelt /= 2 .OR. lambda_last /= lambda .OR. &
         T_a_last /= T_a .OR. Tmax_last /= Tmax ) THEN
    
      ! .. Compute the scattering amplitudes for all sizes in D_S-vector only if something relevant has changed:

      ! Alternativen: MIE_WATERSPH_WETGR_VEC(), MIE_SOAK_WETGR_VEC(),
      !               MIE_MEAN_WETGR_VEC(), MIE_WETSNOW_TWOSPH_VEC()
      CALL MIE_SOAK_TWOSPH_WETGR_BH_VEC(&
           D_g,parti%a_geo,parti%b_geo,fmelt,&
           m_wv,m_iv,aspectratio,ardry,lambda,&
           fa,fb,fa0,fb0,&
           mixingrulestring_iceair, matrixstring_iceair, inclusionstring_iceair, &
           mixingrulestring_icewater, matrixstring_icewater, inclusionstring_icewater, &
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        Tmax_last   = Tmax
        T_a_last    = T_a
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    itw = integweights_mgd(mgd,D_g,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,awetgr,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_wetgr_twosph_mie_vec

  ! Melting graupel: watersphere version
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_wetgr_wsph_mie_vec(&
      mgd,T_a,m_i,m_w,&
      itype_Dref_fmelt,Tmeltbegin,meltdegTmin,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
      Zfac,parti,rain,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring_iceair,matrixstring_iceair,inclusionstring_iceair,Tmax,&
      llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)    :: parti, rain
    TYPE(t_mgd_params), INTENT(in)  :: mgd
    INTEGER, INTENT(in)           :: itype_Dref_fmelt
    LOGICAL, INTENT(in)           :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)     :: pMP, pMPr
    COMPLEX(kind=dp), INTENT(in)  :: m_i, m_w
    REAL(KIND=dp), INTENT(in)     :: T_a,lambda, Zfac, Dmin, Dmax,&
                                     Tmeltbegin, meltdegTmin,Tmax
    CHARACTER(len=*), INTENT(in)  :: mixingrulestring_iceair, matrixstring_iceair,&
                                     inclusionstring_iceair
    LOGICAL, INTENT(in), OPTIONAL :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)    :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL                  :: lmode_lgen

    ! .. Parameters for the degree of melting parameterization with fixed Dref:
    !   ... for itype_Dref_fmelt = 1, the old method with Dref=fct(PSD):
    REAL(KIND=dp), PARAMETER :: Dref_max = 6e-3_dp,   Dexpo_type1 = 0.6_dp
    !   ... for itype_Dref_fmelt = 2, Dref and Dexpo are tuned to approximately reproduce the old method:
    REAL(KIND=dp), PARAMETER :: Dref_type2 = 5e-4_dp, Dexpo_type2 = 0.5_dp

    INTEGER                  :: i
    COMPLEX(kind=dp)         :: m_iv(n_stuetz+1), m_wv(n_stuetz+1), &
                                fa(n_stuetz+1), fb(n_stuetz+1), &
                                fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)            :: dD, &
                                D_r(n_stuetz+1), D_g(n_stuetz+1), &
                                itw(n_stuetz+1), fmelt(n_stuetz+1), &
                                aspectratio(n_stuetz+1), sig, sigr, &
                                arrain(n_stuetz+1),ardry(n_stuetz+1),awetgr(7,n_stuetz+1)

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last, Tmax_last, T_a_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      Tmax_last   = -HUGE(1.0_dp)
      T_a_last    = -HUGE(1.0_dp)
      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_g as SAVEd variables
    !           and recalc them only if Dmin or Dmax have changed (or n_stuetz,
    !           but that is not an input param but hardcoded herein, hence not
    !           changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_g(i) = Dmin + DBLE(i-1) * dD
    END DO

    m_iv = m_i
    m_wv = m_w

    ! Parameterisation of degree of melting
    !
    SELECT CASE (itype_Dref_fmelt)
    CASE (1)
      ! Old method with Dref=fct(PSD mean size):
      CALL degree_of_melting_xxx(mgd,T_a,D_g,n_stuetz+1,Dexpo_type1,Dref_max,parti,&
                                 Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    CASE (2)
      ! New method with fixed Dref = 5e-4 m instead of Dref from PSD, to save computing time below:
      CALL degree_of_melting_xxx_Dref(T_a,D_g,n_stuetz+1,Dexpo_type2,Dref_type2,&
                                      Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    END SELECT

    IF (luse_tmatrix) THEN
      ardry = aspectratio_singlephase(D_g,n_stuetz+1,'graupel',graupel,pMP,sig)

      D_r = ((parti%a_geo * D_g**parti%b_geo)/rain%a_geo)**(1d0/rain%b_geo)
      arrain = aspectratio_singlephase(D_r,n_stuetz+1,'rain',rain,pMPr,sigr)

      IF (.NOT. ldo_nonsphere) THEN
        aspectratio = 1d0
        ardry = 1d0
      ELSE
        ! JM200515: Would be nice to use a WHERE construct here, too. But I'm not
        !           sure how nested WHEREs (here fmelt-limits on a first level,
        !           D_g limits on a second (for fmelt<0.2).

        ! Or, we do two separate WHERE constructs instead...
        WHERE (fmelt < 0.2d0)
          aspectratio = ardry-5.0d0*(ardry - 0.8d0)*fmelt
        ELSEWHERE (fmelt < 0.8)
          aspectratio = 0.88d0 - 0.40d0*fmelt
        ELSEWHERE
          aspectratio = 2.8d0 - 4.0d0*arrain+5.0d0*(arrain-0.56d0)*fmelt
        END WHERE
      END IF

      CALL calc_orient_mixedphase(n_stuetz+1,fmelt,sig,sigr,awetgr)
    ELSE
      aspectratio = 1d0
      ardry = 1d0
      awetgr = 0d0
    END IF

    IF ( .NOT.lmode_lgen .OR. itype_Dref_fmelt /= 2 .OR. lambda_last /= lambda .OR. &
         T_a_last /= T_a .OR. Tmax_last /= Tmax ) THEN
    
      ! .. Compute the scattering amplitudes for all sizes in D_S-vector only if something relevant has changed:

      ! Alternativen: MIE_SOAK_TWOSPH_WETGR_VEC(), MIE_SOAK_WETGR_VEC(),
      !               MIE_MEAN_WETGR_VEC(), MIE_WETSNOW_TWOSPH_VEC()
      CALL MIE_WATERSPH_WETGR_BH_VEC(&
           D_g,parti%a_geo,parti%b_geo,fmelt,&
           m_wv,m_iv,aspectratio,ardry,lambda,&
           fa,fb,fa0,fb0,&
           mixingrulestring_iceair, matrixstring_iceair, inclusionstring_iceair, &
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        Tmax_last   = Tmax
        T_a_last    = T_a
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    itw = integweights_mgd(mgd,D_g,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,awetgr,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_wetgr_wsph_mie_vec

  ! Melting hail:
  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_wethail_mie_vec(&
      mgd,T_a,m_i,m_w,&
      itype_Dref_fmelt,Tmeltbegin,meltdegTmin,lambda,&
      luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
      Zfac,parti,rain,Dmin,Dmax,&
      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
      mixingrulestring,matrixstring,inclusionstring,Tmax,&
      llookupgen_mode)

    IMPLICIT NONE

    TYPE(particle), INTENT(in)     :: parti, rain
    TYPE(t_mgd_params), INTENT(in) :: mgd
    INTEGER, INTENT(in)            :: itype_Dref_fmelt
    LOGICAL, INTENT(in)            :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)      :: pMP, pMPr
    COMPLEX(kind=dp), INTENT(in)   :: m_i, m_w
    REAL(KIND=dp), INTENT(in)      :: T_a,lambda, Zfac, Dmin, Dmax,&
                                      Tmeltbegin, meltdegTmin, Tmax
    CHARACTER(len=*), INTENT(in)   :: mixingrulestring,matrixstring,inclusionstring
    LOGICAL, INTENT(in), OPTIONAL  :: llookupgen_mode

    REAL(KIND=dp), INTENT(out)     :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    LOGICAL                  :: lmode_lgen

    ! .. Parameters for the degree of melting parameterization with fixed Dref:
    !   ... for itype_Dref_fmelt = 1, the old method with Dref=fct(PSD):
    REAL(KIND=dp), PARAMETER :: Dref_max = 6e-3_dp,   Dexpo_type1 = 0.8_dp
    !   ... for itype_Dref_fmelt = 2, Dref and Dexpo are tuned to approximately reproduce the old method:
    REAL(KIND=dp), PARAMETER :: Dref_type2 = 5e-4_dp, Dexpo_type2 = 0.7_dp

    INTEGER                  :: i
    COMPLEX(kind=dp)         :: m_iv(n_stuetz+1), m_wv(n_stuetz+1), &
                                fa(n_stuetz+1), fb(n_stuetz+1), &
                                fa0(n_stuetz+1), fb0(n_stuetz+1)
    REAL(KIND=dp)            :: dD, &
                                D_r(n_stuetz+1), D_g(n_stuetz+1), &
                                itw(n_stuetz+1), &
                                fwater_shell(n_stuetz+1), fmelt(n_stuetz+1), &
                                aspectratio(n_stuetz+1), sig, sigr, &
                                arrain(n_stuetz+1), ardry(n_stuetz+1), awethail(7,n_stuetz+1)

    INTEGER,          SAVE :: firstcall = 0
    REAL(KIND=dp),    SAVE :: simpson(n_stuetz+1)
    REAL(KIND=dp),    SAVE :: lambda_last, Tmax_last, T_a_last
    COMPLEX(kind=dp), SAVE :: fa_last(n_stuetz+1), fb_last(n_stuetz+1), &
                              fa0_last(n_stuetz+1), fb0_last(n_stuetz+1)

    IF (PRESENT(llookupgen_mode)) THEN
      lmode_lgen = llookupgen_mode
    ELSE
      lmode_lgen = .FALSE.
    END IF
    
    IF (firstcall /= 1) THEN
      simpson = simpson_weights(n_stuetz)
      lambda_last = -HUGE(1.0_dp)
      Tmax_last   = -HUGE(1.0_dp)
      T_a_last    = -HUGE(1.0_dp)
      firstcall = 1
    ENDIF

    ! FIXME
    ! JM200515: we should also be able to have the dD and D_g as SAVEd variables
    !           and recalc them only if Dmin or Dmax have changed (or n_stuetz,
    !           but that is not an input param but hardcoded herein, hence not
    !           changing)
    dD   = (Dmax-Dmin)/n_stuetz
    DO i=1,n_stuetz+1
      D_g(i) = Dmin + DBLE(i-1) * dD
    END DO

    m_iv = m_i
    m_wv = m_w

    ! Parameterisation of degree of melting
    !
    SELECT CASE (itype_Dref_fmelt)
    CASE (1)
      ! Old method with Dref=fct(PSD mean size):
      CALL degree_of_melting_xxx(mgd,T_a,D_g,n_stuetz+1,Dexpo_type1,Dref_max,parti,&
                                 Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    CASE (2)
      ! New method with fixed Dref = 5e-4 m instead of Dref from PSD, to save computing time below:
      CALL degree_of_melting_xxx_Dref(T_a,D_g,n_stuetz+1,Dexpo_type2,Dref_type2,&
                                      Tmeltbegin,meltdegTmin,Tmin_f,Tmax,fmelt)
    END SELECT

    IF (luse_tmatrix) THEN
      ardry = aspectratio_singlephase(D_g,n_stuetz+1,'hail',hail,pMP,sig)

      D_r = ((parti%a_geo * D_g**parti%b_geo)/rain%a_geo)**(1d0/rain%b_geo)
      arrain = aspectratio_singlephase(D_r,n_stuetz+1,'rain',rain,pMPr,sigr)

      IF (.NOT. ldo_nonsphere) THEN
        aspectratio = 1d0
        ardry = 1d0
      ELSE
        ! JM200515: Would be nice to use a WHERE construct here, too. But I'm not
        !           sure how nested WHEREs (here fmelt-limits on a first level,
        !           D_g limits on a second (for fmelt<0.2).

        WHERE (fmelt < 0.2d0)
          aspectratio = ardry-5.0d0*(ardry - 0.8d0)*fmelt
        ELSEWHERE (fmelt < 0.8)
          aspectratio = 0.88d0 - 0.40d0*fmelt
        ELSEWHERE
          aspectratio = 2.8d0 - 4.0d0*arrain+5.0d0*(arrain-0.56d0)*fmelt
        END WHERE
      END IF

      CALL calc_orient_mixedphase(n_stuetz+1,fmelt,sig,sigr,awethail)
    ELSE
      aspectratio = 1d0
      ardry = 1d0
      awethail = 0d0
    END IF

    IF ( .NOT.lmode_lgen .OR. itype_Dref_fmelt /= 2 .OR. lambda_last /= lambda .OR. &
         T_a_last /= T_a .OR. Tmax_last /= Tmax   ) THEN
    
      ! .. Compute the scattering amplitudes for all sizes in D_S-vector only if something relevant has changed:

      fwater_shell = 0.5d0

      ! Alternativen: MIE_SPONGY_WETHAIL(), MIE_WATERSPH_WETHAIL(), MIE_SOAK_WETGR()
      !
      ! MIE_SPONGY_WETHAIL assumes solid ice dry particle.
      ! That is, D_g here is (already) the Dveq of the particle, hence there is
      ! no need of ardry to derive volume/density/Dveq.
      CALL MIE_SPONGY_WETHAIL_BH_VEC(&
           D_g,parti%a_geo,parti%b_geo,fmelt,fwater_shell,&
           m_iv,m_wv,aspectratio,lambda,&
           fa,fb,fa0,fb0,&
           mixingrulestring,matrixstring,inclusionstring,&
           n_stuetz+1,luse_tmatrix)

      IF (lmode_lgen) THEN
        lambda_last = lambda
        Tmax_last   = Tmax
        T_a_last    = T_a
      
        ! .. Store the scattering amplitudes for subsequent computations with the same relevant input parameters:
        fa_last  = fa
        fb_last  = fb
        fa0_last = fa0
        fb0_last = fb0
      END IF
      
    ELSE

      ! .. Otherwise, overtake the save'd scattering amplitudes from the last computation:
      fa  = fa_last 
      fb  = fb_last 
      fa0 = fa0_last
      fb0 = fb0_last
      
    END IF

    ! Size distribution integration weights
    itw = integweights_mgd(mgd,D_g,dD,simpson,n_stuetz+1)

    CALL calc_radar_vars(lambda,Zfac,awethail,itw,fa,fb,fa0,fb0,n_stuetz+1,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh)

  END SUBROUTINE zradar_wethail_mie_vec


  !*******************************************************************************
  !
  ! Polarimetric radar parameter & oriented particle related routines
  !
  !*******************************************************************************

  FUNCTION aspectratio_singlephase(D,anz,hymet,parti,pMP,sig) & !,ARmin_safeTmat) &
      RESULT (AR)

    IMPLICIT NONE

    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D(anz)
    CHARACTER(len=*), INTENT(in) :: hymet
    TYPE(particle), INTENT(in)   :: parti
    TYPE(t_polMP), INTENT(in)    :: pMP

    REAL(KIND=dp), INTENT(out) :: sig

    REAL(KIND=dp) :: AR(anz)

    REAL(KIND=dp)            :: rho(anz)
    REAL(KIND=dp), PARAMETER :: ARmin_safeTmat = 0.2d0

    sig = pMP%sig

    SELECT CASE (toupper(TRIM(ADJUSTL(pMP%ARmodel))))
    CASE('POLY')
      AR = MAX(pMP%c0 + pMP%c1*D,pMP%ARmin)

    CASE('B02R')
      !IF (hymet .EQ. 'rain') THEN
        ! rain aspect ratio according to Brandes et al, JAM, 2002
        AR = MIN(MAX(0.9951d0 + 0.0251d0*(D*1d3) - 0.03644d0*(D*1d3)**2 +  &
                     0.005303d0*(D*1d3)**3 - 0.0002492d0*(D*1d3)**4, &
                   0.56d0), 1.0d0)
      !ELSE
      !  errormessage & abort

    ! JM: Use Andric13 with habit-specific rho-D-parametrization and derive
    ! AR consistently from that Andric13 rho-D relation
    ! HOWEVER, the choice is practically irrelevant with the low-AR threshold
    ! for TMatrix stability applied below
    CASE('A13IP')
      ! Andric13 plates
      rho = 769.8318d0 * D**(-0.014d0)
      ! adjust predicted rho if unphysical
      rho = MIN(rho,rho_ice_fwo)
      ! general AR relation for oblate(!!!) spheroids
      AR = 6d0/(pi_dp*rho) * (1d0/parti%a_geo)**(1d0/parti%b_geo) * &
           D**(1d0/parti%b_geo-3d0)
    CASE('A13ID')
      ! Andric13 dendrites
      rho =  43.4888d0 * D**(-0.377d0)
      ! adjust predicted rho if unphysical
      rho = MIN(rho,rho_ice_fwo)
      ! general AR relation for oblate(!!!) spheroids
      AR = 6d0/(pi_dp*rho) * (1d0/parti%a_geo)**(1d0/parti%b_geo) * &
           D**(1d0/parti%b_geo-3d0)

    CASE('M96ITP')
      ! Matrosov96 solid thick plates
      ! h = a2 * D**f
      ! AR = h / D = a2 * D**f / D = a2 * D**(f-1) with D [cm]
      AR =  0.009d0 * (D*1d2)**(0.377d0-1d0)

    CASE('R11')
      ! R11 is used as "full-MP" model, ie also modifying sig
      SELECT CASE (TRIM(ADJUSTL(hymet)))
      CASE('rain')
        ! rain aspect ratio according to Brandes et al, JAM, 2002 as suggested by R11
        AR = MIN(MAX(0.9951d0 + 0.0251d0*(D*1d3) - 0.03644d0*(D*1d3)**2 +  &
                     0.005303d0*(D*1d3)**3 - 0.0002492d0*(D*1d3)**4, &
                   0.56d0), 1.0d0)
        IF (sig .LT. 0d0) sig = 10d0

      CASE('ice')
        ! cloud ice AR from Matrosov96 as suggested by R11
        ! Assuming solid thick plates
        !   (rho-D relation of dendrites better fit ICON/COSMO, but both dendrites
        !   and hexagonal plates provide very small AR, ie would be set to
        !   TMat-safe ARmin_safeTmat(=0.2) everywhere but the ~3 smallest sizes)
        ! h = a2 * D**f
        ! AR = h / D = a2 * D**f / D = a2 * D**(f-1) with D [cm]
        AR =  0.009d0 * (D*1d2)**(0.377d0-1d0)
        IF (sig .LT. 0d0) sig = 10d0

      CASE('snow')
        ! that's our original implementation.
        !WHERE (D < 10.0d-3)
        !  AR = 1.0d0 - 0.02d0*(D*1d3)
        !ELSEWHERE
        !  AR = 0.8d0
        !ENDWHERE
        ! that's the actual R11 suggestion:
        AR = 0.8d0
        IF (sig .LT. 0d0) sig = 40d0

      CASE('graupel','hail')
        WHERE (D < 10.0d-3)
          AR = 1.0d0 - 0.02d0*(D*1d3)
        ELSEWHERE
          AR = 0.8d0
        ENDWHERE
        IF (sig .LT. 0d0) sig = 40d0

      END SELECT

    CASE('BPRO')
      ! R11 is used as "full-MP" model, ie also modifying sig
      SELECT CASE (TRIM(ADJUSTL(hymet)))
      CASE('rain')
        ! rain aspect ratio according to Brandes et al, JAM, 2002 as suggested by R11
        AR = MIN(MAX(0.9951d0 + 0.0251d0*(D*1d3) - 0.03644d0*(D*1d3)**2 +  &
                     0.005303d0*(D*1d3)**3 - 0.0002492d0*(D*1d3)**4, &
                   0.56d0), 1.0d0)
        IF (sig .LT. 0d0) sig = 10d0

      CASE('ice')
        ! Andric13 plates
        rho = 769.8318d0 * D**(-0.014d0)
        ! adjust predicted rho if unphysical
        rho = MIN(rho,rho_ice_fwo)
        ! general AR relation for oblate(!!!) spheroids
        AR = 6d0/(pi_dp*rho) * (1d0/parti%a_geo)**(1d0/parti%b_geo) * &
             D**(1d0/parti%b_geo-3d0)
        IF (sig .LT. 0d0) sig = 12d0

      CASE('snow')
        AR = MAX(0.7d0 - 0.01d0*(D*1e3), 0.5d0)
        IF (sig .LT. 0d0) sig = 40d0

      CASE('graupel','hail')
        AR = MAX(1.0d0 - 0.02d0*(D*1e3), 0.8d0)
        IF (sig .LT. 0d0) sig = 40d0
        
      END SELECT

    CASE default
      ! JM221103: replace this by a proper EMVORADO abort
      WRITE (*,*) 'ERROR in aspectratio_singlephase: '//&
                  'Invalid '//TRIM(ADJUSTL(hymet))//' ARmodel ', pMP%ARmodel
      STOP
    END SELECT

    IF (sig .LT. 0) THEN
      ! JM221103: replace this by a proper EMVORADO abort
      WRITE (*,*) 'ERROR in aspectratio_singlephase: '//&
                  'Missing setting of '//TRIM(ADJUSTL(hymet))//' orientation distribution width'
      STOP
    END IF

    ! global max
    WHERE (AR.GT.1d0)
      AR = 1d0
    ! global min
    ELSEWHERE (AR.LT.ARmin_safeTmat)
      AR = ARmin_safeTmat
    END WHERE
  END FUNCTION aspectratio_singlephase

  SUBROUTINE calc_orient_singlephase(anz,sigma,a)
    IMPLICIT NONE
    INTEGER,       INTENT(IN)  :: anz
    REAL(KIND=dp), INTENT(IN)  :: sigma
    REAL(KIND=dp), INTENT(OUT) :: a(7,anz)

    REAL(KIND=dp) :: sigma_eff, r

    sigma_eff = sigma/180.0d0*pi_dp  ! witdh of hydrometeor canting angle distr. [radians]
    r         = EXP(-2.0d0*sigma_eff**2)

    CALL calc_angmom(anz,r,a)

  END SUBROUTINE calc_orient_singlephase

  SUBROUTINE calc_orient_mixedphase(anz,vmelt,sigma_i,sigma_w,a)
    IMPLICIT NONE
    INTEGER,       INTENT(in)  :: anz
    REAL(KIND=dp), INTENT(IN)  :: vmelt(anz)       ! volume meltfraction
    REAL(KIND=dp), INTENT(in)  :: sigma_i, sigma_w ! width of canting angle distr. of (dry)frozen and liquid hydrometeors [degrees]
    REAL(KIND=dp), INTENT(OUT) :: a(7,anz)

    INTEGER       :: i
    REAL(KIND=dp) :: sigma_eff, r

    do i = 1, anz
       sigma_eff = (sigma_i + (sigma_w - sigma_i)*vmelt(i)) / 180d0*pi_dp
       r         = EXP(-2.0d0*sigma_eff**2)

      CALL calc_angmom(1,r,a(:,i:i))
    end do

    return
  END SUBROUTINE calc_orient_mixedphase

  SUBROUTINE calc_angmom(anz,r,a)
    IMPLICIT NONE
    INTEGER,       INTENT(IN)  :: anz
    REAL(KIND=dp), INTENT(IN)  :: r
    REAL(KIND=dp), INTENT(OUT) :: a(7,anz)

    a(1,:)= 0.25d0*(1.0+r)**2
    a(2,:)= 0.25d0*(1.0-r**2)
    a(3,:)= (0.375d0+0.5d0*r+0.125d0*r**4)**2
    a(4,:)= (0.375d0-0.5d0*r+0.125d0*r**4)* &
          (0.375d0+0.5d0*r+0.125d0*r**4)
    a(5,:)= 0.125d0*(0.375d0+0.5d0*r+0.125d0*r**4)*(1.0d0-r**4)
    a(6,:)= 0.0d0
    a(7,:)= 0.5d0*r*(1.0d0+r)
  END SUBROUTINE calc_angmom


  ! JCS - Calculate relevant radar variables
  ! JM/UB: Zfac has to be 1e18 * lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE calc_radar_vars(lambda,Zfac,a,NddD,fa,fb,fa0,fb0,anz,&
                             zh,ah,zv,rrhv,irhv,kdp,adp,zvh)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: lambda, Zfac
    INTEGER, INTENT(in) :: anz
    REAL(KIND=dp), DIMENSION(7,anz), INTENT(IN) :: a
    COMPLEX(KIND=dp), DIMENSION(anz), INTENT(IN) :: fa, fb, fa0, fb0
    REAL(KIND=dp), DIMENSION(anz), INTENT(IN) :: NddD

    REAL(KIND=dp), INTENT(OUT) :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh

    COMPLEX(KIND=dp) :: crhv
    REAL(KIND=dp) :: c1, c2, c3
    INTEGER :: i

    ! JM200513: Remove order-of-mag conversion factors from the f*'s, and hence
    !           equivalently from the c* (why all that back-and-forth converting
    !           in the first place???). Saves local storing of lambda and the
    !           f*'s and brings Z-results in agreement with singlepol reference
    !           without further factors onto c1.
    c1 = 4d0*pi_dp * Zfac
    c2 = 180d0*lambda/pi_dp
    c3 = 2d0*lambda

    zh = 0d0
    ah = 0d0
    zv = 0d0
    crhv = CMPLX(0d0, 0d0, kind=dp)
    kdp = 0d0
    adp = 0d0
    zvh = 0d0

    IF (ANY(fa == miss_value) .OR. ANY(fb == miss_value)) THEN
      zh = miss_value
      zv = miss_value
      rrhv = miss_value_rhv
      irhv = miss_value_rhv
      zvh = miss_value
    ELSE

!!$      DO i=1,anz
!!$        ! JM191021: multiplication with c1/2/3 constants moved out of loop.
!!$        !           should be (slightly) more computationally effective.
!!$        zh = zh + NddD(i) * &
!!$                  ( ABS(fb(i))**2 - &
!!$                    2d0*a(2,i)*REAL(CONJG(fb(i))*(fb(i)-fa(i))) + &
!!$                    a(4,i)*ABS(fb(i)-fa(i))**2 )
!!$        zv = zv + NddD(i) * &
!!$                  ( ABS(fb(i))**2 - &
!!$                    2d0*a(1,i)*REAL(CONJG(fb(i))*(fb(i)-fa(i))) + &
!!$                    a(3,i)*ABS(fb(i)-fa(i))**2 )
!!$        zvh = zvh + NddD(i) * &
!!$                    ( a(5,i) * (ABS(fb(i)-fa(i))**2) )
!!$        crhv = crhv + NddD(i) * &
!!$                     ( ABS(fb(i))**2 + &
!!$                       a(5,i)*ABS(fb(i)-fa(i))**2 - &
!!$                       a(1,i)*CONJG(fb(i))*(fb(i)-fa(i)) - &
!!$                       a(2,i)*fb(i)*CONJG(fb(i)-fa(i)) )
!!$      END DO
!!$      zh = zh*c1
!!$      zv = zv*c1
!!$      zvh = zvh*c1
!!$      crhv = crhv*c1
!!$
      
      zh = c1 * SUM( NddD(:) * &
                ( ABS(fb(:))**2 - &
                  2d0*a(2,:)*REAL(CONJG(fb(:))*(fb(:)-fa(:))) + &
                  a(4,:)*ABS(fb(:)-fa(:))**2 ) )
      zv = c1 * SUM( NddD(:) * &
                ( ABS(fb(:))**2 - &
                 2d0*a(1,:)*REAL(CONJG(fb(:))*(fb(:)-fa(:))) + &
                 a(3,:)*ABS(fb(:)-fa(:))**2 ) )
      zvh = c1 * SUM( NddD(:) * a(5,:) * ABS(fb(:)-fa(:))**2 )
      crhv = c1 * SUM( NddD(:) * &
                  ( ABS(fb(:))**2 + &
                   a(5,:)*ABS(fb(:)-fa(:))**2 - &
                   a(1,:)*CONJG(fb(:))*(fb(:)-fa(:)) - &
                   a(2,:)*fb(:)*CONJG(fb(:)-fa(:)) ) )
      ! JM191021:
      ! Don't do any non-additive conversions here. This needs to wait until the very end!
      !zdr = zh/zv !10*log10(zh(i)/zv(i))
      !ldr = 10.0d0*log10(ldr)
      !rhv = rhv/((zh*zv)**(0.5d0))

      rrhv = REAL(crhv)
      irhv = AIMAG(crhv)

    END IF

    IF (ANY(fa0 == miss_value) .OR. ANY(fb0 == miss_value)) THEN
      kdp = miss_value
      ah = miss_value
      adp = miss_value
    ELSE
!!$      DO i=1,anz
!!$        kdp = kdp + NddD(i) * &
!!$!                    ( a(7,i) * real(fb0(i)-fa0(i)) )
!!$                    ( a(7,i) * REAL(fa0(i)-fb0(i)) ) ! signs of fb0 and fa0 seem inverted with our
!!$                                                     ! scat.amp. setups.
!!$                                                     ! kdp of horizontal oblates should be positive.
!!$        ah =   ah + NddD(i) * &
!!$                   ( AIMAG(fb0(i)) - &
!!$                     a(2,i)*AIMAG(fb0(i)-fa0(i)) )
!!$        adp = adp + NddD(i) * &
!!$                    ( a(7,i) * AIMAG(fb0(i)-fa0(i)) )
!!$      END DO
!!$      kdp = kdp*c2
!!$      ah = ah*c3
!!$      adp = adp*c3

      kdp = c2 * SUM( NddD(:) * &
!                ( a(7,:) * real(fb0(:)-fa0(:)) ) )
                 ( a(7,:) * REAL(fa0(:)-fb0(:)) ) ) ! signs of fb0 and fa0 seem inverted with our
                                                     ! scat.amp. setups.
                                                     ! kdp of horizontal oblates should be positive.
      ah = c3 * SUM( NddD(:) * &
                ( AIMAG(fb0(:)) - &
                     a(2,:)*AIMAG(fb0(:)-fa0(:)) ) )
      adp = c3 * SUM( NddD(:) * &
                 ( a(7,:) * AIMAG(fb0(:)-fa0(:)) ) )

    END IF

  END SUBROUTINE calc_radar_vars


  !=========================================================================================
  ! begin lookup tables for radar reflectivity following Mie theory

  !*******************************************************************************
  !
  ! Create Lookup-table vectors for radar reflectivity and its extinction
  ! for dry hydrometeors using Mie theory and 1-moment-scheme as a function
  ! of temperature (T_a) and specific content (q_i) of rain water (q_i=q_r),
  ! dry snow (q_i=q_s) and dry graupel (q_i=q_g)
  !
  !*******************************************************************************

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_rain_1mom_lookupcreate (   &
       look_Z_rain,tableprops,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,&
       Zfac,Dmin,Dmax,&
       impipar_lookupgen,pe_start,pe_end,&
       linterp_mode_dualpol,&
       savepath_read,savepath_write,unitnr,hashtext)

    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_rain
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP
    INTEGER, INTENT(in)                   :: unitnr
    REAL(KIND=dp), INTENT(in)             :: lambda_radar, Zfac, Dmin, Dmax
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_w(tableprops%nTa+1)
    REAL(KIND=dp)       :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh
    REAL(KIND=dp)       :: zhu,ahu,zvu,rrhvu,irhvu,kdpu,adpu,zvhu,dqi
    INTEGER             :: i,j,k,kk,ierr,work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_rain_1mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile,tmppath_read,tmppath_write
    CHARACTER(len=10)   :: cmagicnr,cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist,print_debug

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_rain, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = 0, &                 ! table T_m will have 1 node
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = mw_Tmin + T0C_fwo, & ! K
         &                    Taup  = mw_Tmax + T0C_fwo, & ! K
         &                    Tmlow = T0C_fwo, &           ! for rain, this is a dummy
         &                    Tmup  = T0C_fwo, &           ! for rain, this is a dummy
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qr
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_rain%zh,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%ah,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%zv,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%rrhv, lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%irhv, lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%kdp,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%adp,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%zvh,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_rain%zh,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_rain%ah,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_rain%zv,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_rain%rrhv, lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%irhv, lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%kdp,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%adp,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%zvh,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF
    
    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_rain%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_1mom_rain_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_rain%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_rain%mem   = 0.0_dp   ! Entire memory blocks in one go
      look_Z_rain%dmem  = 0.0_dp   ! ... also for the derivatives w.r.t. qi
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for rain, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_w = m_complex_water_ray_vec(DBLE(lambda_radar), &
                                    MIN(MAX(look_Z_rain%T_a-T0C_fwo,mw_Tmin),mw_Tmax), &
                                    look_Z_rain%nTa)

      ! fill lookup table body
      DO j=1, look_Z_rain%nTa

        IF (impipar_lookupgen == 2) THEN
          work_pe = j
          work_pe = round_robin(work_pe-1, pe_start, pe_end)
        ELSE
          work_pe = my_cart_id_fwo
        END IF

        IF (my_cart_id_fwo == work_pe) THEN

          DO i=1, look_Z_rain%nqi

            IF (luse_tmatrix .AND. ldebug_dbz) THEN
              ! Computation takes long, so print some progress information
              CALL progress_information (pos_=[i,j], &
                   shape_=[look_Z_rain%nqi,look_Z_rain%nTa], &
                   stride_=[look_Z_rain%nqi,1], cident='Rain lookup creation')
            END IF
            
            ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
            ! .. Suitable value for dqi, which is smaller than the local grid spacing:
            IF (i > 1) THEN
              ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
              !  take the dqi for differentiation as a fraction of the "lower" diff (to lower neighbour):
              dqi = 0.01_dp * (look_Z_rain%q_i_lin(i)-look_Z_rain%q_i_lin(i-1))
            ELSE
              ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
              !  to avoid negative q_i_lin below:
              dqi = 0.01_dp * (look_Z_rain%q_i_lin(i+1)-look_Z_rain%q_i_lin(i))
              dqi = MIN(dqi, 0.05_dp * look_Z_rain%q_i_lin(i))
            END IF

            DO kk=1, 3

              SELECT CASE (kk)
              CASE (1)
                mgd = mgd_1mom(rain,look_Z_rain%q_i_lin(i),rain%n0_const)
              CASE (2)
                mgd = mgd_1mom(rain,look_Z_rain%q_i_lin(i)-dqi,rain%n0_const)
              CASE (3)
                mgd = mgd_1mom(rain,look_Z_rain%q_i_lin(i)+dqi,rain%n0_const)
              END SELECT
            
              CALL zradar_rain_mie_vec(mgd, m_w(j), lambda_radar,              &
                                       luse_tmatrix, ldo_nonsphere, pMP,       &
                                       Zfac, Dmin, Dmax,                       &
                                       zh, ah, zv, rrhv, irhv, kdp, adp, zvh , &
                                       llookupgen_mode=.TRUE. )

              
              SELECT CASE (kk)
              CASE (1)
                ! Store parameters in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_rain % zh   % val(i,j,1) = zh
                look_Z_rain % ah   % val(i,j,1) = ah
                look_Z_rain % zv   % val(i,j,1) = zv
                look_Z_rain % rrhv % val(i,j,1) = rrhv
                look_Z_rain % irhv % val(i,j,1) = irhv
                look_Z_rain % kdp  % val(i,j,1) = kdp
                look_Z_rain % adp  % val(i,j,1) = adp
                look_Z_rain % zvh  % val(i,j,1) = zvh
              CASE (2)
                zhu   = zh
                ahu   = ah
                zvu   = zv
                rrhvu = rrhv
                irhvu = irhv
                kdpu  = kdp
                adpu  = adp
                zvhu  = zvh
              CASE (3)
                ! Store derivatives in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_rain % zh   % dval(i,j,1) = (zh-zhu)     / dqi * 0.5_dp
                look_Z_rain % ah   % dval(i,j,1) = (ah-ahu)     / dqi * 0.5_dp
                look_Z_rain % zv   % dval(i,j,1) = (zv-zvu)     / dqi * 0.5_dp
                look_Z_rain % rrhv % dval(i,j,1) = (rrhv-rrhvu) / dqi * 0.5_dp
                look_Z_rain % irhv % dval(i,j,1) = (irhv-irhvu) / dqi * 0.5_dp
                look_Z_rain % kdp  % dval(i,j,1) = (kdp-kdpu)   / dqi * 0.5_dp
                look_Z_rain % adp  % dval(i,j,1) = (adp-adpu)   / dqi * 0.5_dp
                look_Z_rain % zvh  % dval(i,j,1) = (zvh-zvhu)   / dqi * 0.5_dp
              END SELECT

            END DO
            
          END DO

        END IF
        
      END DO

      ! JM201020: FIXME ( !!! this applies to ALL *_lookupcreate routines !!! )
      ! Should other polarimetric parameters, at least the rest of the Z-based
      ! ones (zvh, rrhv, irhv), be converted and stored in log-space, too?
      !
      ! JM201021 after discussion with UB:
      ! As also ext is stored (and interpolated) in log-space, it makes sense to
      ! do that for all other parameters, too. Shouldn't be a prob for Z-type
      ! parameters (zv, zvh, rrhv, irhv). With adp and kdp, an issue is posed by
      ! that they can be negative, ie no log possible. Possible solution or
      ! workaround can be to store (and interpolate) av instead of adp and the
      ! separate components Re(Shh0) and Re(Svv0) that form kdp. Makes yet an
      ! additional parameter to be handled, but on the other side allows for
      ! identical treatment of all parameters.

      look_Z_rain%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) look_Z_rain%magicnr,nTa,nqr,luse_tmatrix,ldo_nonsphere,&
      !     lambda_radar,Dmin,Dmax,mgd%n0,mgd%mu,' '//TRIM(look_Z_rain%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_rain, yzroutine, TRIM(hashtext))

    ELSE
      
      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: reading lookup table for rain, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_rain, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_rain%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_rain%nparams
      ! compute df(p(qi(qi_scaled)))/dqi_scaled = df/dp * dp/dqi * dqi/dqi_scaled = df/dp * dval * dqi/dqi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_rain%dmem(:,:,:,k) =  scale_dval_lut(look_Z_rain%mem(:,:,:,k), look_Z_rain%dmem(:,:,:,k), &
           look_Z_rain%qmem(:,k), &
           look_Z_rain%flag_qi_scal(k), look_Z_rain%f_scal_qi(k), &
           look_Z_rain%flag_mem_scal(k), look_Z_rain%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_rain%mem(:,:,:,k) =  scale_val_lut(look_Z_rain%mem(:,:,:,k), &
           look_Z_rain%flag_mem_scal(k), look_Z_rain%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_rain_1mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_ice_1mom_lookupcreate(&
       look_Z_ice,tableprops,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,&
       Zfac,Dmin,Dmax,ice,&
       Tmeltbegin,mixingrulestring,matrixstring,inclusionstring,&
       impipar_lookupgen,pe_start,pe_end,&
       linterp_mode_dualpol,&
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_ice
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    TYPE(particle), INTENT(in)            :: ice
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP
    INTEGER, INTENT(in)                   :: unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, lambda_radar, Zfac, Dmin, Dmax
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring,&
                                             savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1)
    REAL(KIND=dp)       :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh
    REAL(KIND=dp)       :: zhu,ahu,zvu,rrhvu,irhvu,kdpu,adpu,zvhu,dqi,n_i
    INTEGER             :: i,j,k,kk,ierr,work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_ice_1mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile,tmppath_read,tmppath_write
    CHARACTER(len=10)   :: cmagicnr,cscattheo
    CHARACTER(len=4)    :: saveext
    CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist,print_debug

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_ice, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = 0, &                 ! table T_m will have 1 node
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = mi_Tmin + T0C_fwo, & ! K
         &                    Taup  = Tmeltbegin, &        ! K
         &                    Tmlow = T0C_fwo, &           ! for dry ice, this is a dummy
         &                    Tmup  = T0C_fwo, &           ! for dry ice, this is a dummy
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qi
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_ice%zh,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%ah,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%zv,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%rrhv, lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%irhv, lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%kdp,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%adp,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%zvh,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_ice%zh,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_ice%ah,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_ice%zv,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_ice%rrhv, lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%irhv, lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%kdp,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%adp,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%zvh,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_ice%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_1mom_dryice_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_ice%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_ice%mem  = 0.0_dp
      look_Z_ice%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT

    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for dry ice, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_ice%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_ice%nTa)

      ! fill lookup table body
      DO j=1, look_Z_ice%nTa

        IF (impipar_lookupgen == 2) THEN
          work_pe = j
          work_pe = round_robin(work_pe-1, pe_start, pe_end)
        ELSE
          work_pe = my_cart_id_fwo
        END IF

        IF (my_cart_id_fwo == work_pe) THEN

          n_i = nice_mono_1mom(look_Z_ice%T_a(j))
              
          DO i=1, look_Z_ice%nqi
            
            IF (luse_tmatrix .AND. ldebug_dbz) THEN
              ! Computation takes long, so print some progress information
              CALL progress_information (pos_=[i,j], &
                   shape_=[look_Z_ice%nqi,look_Z_ice%nTa], &
                   stride_=[look_Z_ice%nqi,1], cident='Dry ice lookup creation')
            END IF
            
            ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
            ! .. Suitable value for dqi, which is smaller than the local grid spacing:
            IF (i > 1) THEN
              ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
              !  take the dqi for differentiation as a fraction of the "lower" diff (to lower neighbour):
              dqi = 0.01_dp * (look_Z_ice%q_i_lin(i)-look_Z_ice%q_i_lin(i-1))
            ELSE
              ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
              !  to avoid negative q_i_lin below:
              dqi = 0.01_dp * (look_Z_ice%q_i_lin(i+1)-look_Z_ice%q_i_lin(i))
              dqi = MIN(dqi, 0.05_dp * look_Z_ice%q_i_lin(i))
            END IF

            DO kk=1, 3

              ! Mimic a monodisperse size distribution around a mean mass of q_i/n_i by using 2-mom mgd method and
              !  choosing very large values for ice%mu and ice%nu (in init_1mom_types(), radar_interface.f90):
              SELECT CASE (kk)
              CASE (1)
                mgd = mgd_2mom(ice,look_Z_ice%q_i_lin(i),n_i)
              CASE (2)
                mgd = mgd_2mom(ice,look_Z_ice%q_i_lin(i)-dqi,n_i)
              CASE (3)
                mgd = mgd_2mom(ice,look_Z_ice%q_i_lin(i)+dqi,n_i)
              END SELECT

              CALL zradar_ice_mie_vec(&
                   mgd,m_i(j),lambda_radar,&
                   luse_tmatrix,ldo_nonsphere,pMP,&
                   Zfac,ice,Dmin,Dmax,&
                   zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                   mixingrulestring,matrixstring,inclusionstring,&
                   llookupgen_mode=.TRUE.)

              SELECT CASE (kk)
              CASE (1)
                ! Store parameters in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_ice % zh   % val(i,j,1) = zh
                look_Z_ice % ah   % val(i,j,1) = ah
                look_Z_ice % zv   % val(i,j,1) = zv
                look_Z_ice % rrhv % val(i,j,1) = rrhv
                look_Z_ice % irhv % val(i,j,1) = irhv
                look_Z_ice % kdp  % val(i,j,1) = kdp
                look_Z_ice % adp  % val(i,j,1) = adp
                look_Z_ice % zvh  % val(i,j,1) = zvh
              CASE (2)
                zhu   = zh
                ahu   = ah
                zvu   = zv
                rrhvu = rrhv
                irhvu = irhv
                kdpu  = kdp
                adpu  = adp
                zvhu  = zvh
              CASE (3)
                ! Store derivatives in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_ice % zh   % dval(i,j,1) = (zh-zhu)     / dqi * 0.5_dp
                look_Z_ice % ah   % dval(i,j,1) = (ah-ahu)     / dqi * 0.5_dp
                look_Z_ice % zv   % dval(i,j,1) = (zv-zvu)     / dqi * 0.5_dp
                look_Z_ice % rrhv % dval(i,j,1) = (rrhv-rrhvu) / dqi * 0.5_dp
                look_Z_ice % irhv % dval(i,j,1) = (irhv-irhvu) / dqi * 0.5_dp
                look_Z_ice % kdp  % dval(i,j,1) = (kdp-kdpu)   / dqi * 0.5_dp
                look_Z_ice % adp  % dval(i,j,1) = (adp-adpu)   / dqi * 0.5_dp
                look_Z_ice % zvh  % dval(i,j,1) = (zvh-zvhu)   / dqi * 0.5_dp
              END SELECT

            END DO
            
          END DO
          
        END IF
      END DO

      look_Z_ice%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_ice%magicnr,nTa,nqi,luse_tmatrix,ldo_nonsphere,&
      !     lambda_radar,Dmin,Dmax,&
      !     ice%a_geo,ice%b_geo,Tmeltbegin,&
      !     mixingrulestring,matrixstring,inclusionstring,' '//TRIM(look_Z_ice%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_ice, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_ice, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_ice%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_ice%nparams
      ! compute df(p(qi(qi_scaled)))/dqi_scaled = df/dp * dp/dqi * dqi/dqi_scaled = df/dp * dval * dqi/dqi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_ice%dmem(:,:,:,k) =  scale_dval_lut(look_Z_ice%mem(:,:,:,k), look_Z_ice%dmem(:,:,:,k), &
           look_Z_ice%qmem(:,k), &
           look_Z_ice%flag_qi_scal(k), look_Z_ice%f_scal_qi(k), &
           look_Z_ice%flag_mem_scal(k), look_Z_ice%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_ice%mem(:,:,:,k) =  scale_val_lut(look_Z_ice%mem(:,:,:,k), &
           look_Z_ice%flag_mem_scal(k), look_Z_ice%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_ice_1mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_snow_1mom_lookupcreate(&
       look_Z_snow,tableprops,isnow_n0temp,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,&
       Zfac,Dmin,Dmax,snow,Tmeltbegin,&
       mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
       mixingrulestring_core,matrixstring_core,inclusionstring_core,&
       impipar_lookupgen,pe_start,pe_end,&
       linterp_mode_dualpol,&
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_snow
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    TYPE(particle), INTENT(in)            :: snow
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP
    INTEGER, INTENT(in)                   :: isnow_n0temp, unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, lambda_radar, Zfac, Dmin, Dmax
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring_shell, matrixstring_shell, &
                                             inclusionstring_shell, &
                                             mixingrulestring_core, matrixstring_core, &
                                             inclusionstring_core, &
                                             savepath_read, savepath_write,hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1)
    REAL(KIND=dp)       :: n0_s
    REAL(KIND=dp)       :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh
    REAL(KIND=dp)       :: zhu,ahu,zvu,rrhvu,irhvu,kdpu,adpu,zvhu,dqi
    INTEGER             :: i,j,k,kk,ierr,work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_snow_1mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile,tmppath_read,tmppath_write
    CHARACTER(len=10)   :: cmagicnr,cscattheo
    CHARACTER(len=4)    :: saveext
    CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist,print_debug

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_snow, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = 0, &                 ! table T_m will have 1 node
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = mi_Tmin + T0C_fwo, & ! K
         &                    Taup  = Tmeltbegin, &        ! K
         &                    Tmlow = T0C_fwo, &           ! for dry snow, this is a dummy
         &                    Tmup  = T0C_fwo, &           ! for dry snow, this is a dummy
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qr
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_snow%zh,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%ah,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%zv,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%rrhv, lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%irhv, lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%kdp,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%adp,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%zvh,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_snow%zh,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_snow%ah,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_snow%zv,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_snow%rrhv, lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%irhv, lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%kdp,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%adp,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%zvh,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_snow%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_1mom_drysnow_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_snow%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_snow%mem = 0.0_dp
      look_Z_snow%dmem = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for dry snow, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_snow%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_snow%nTa)

      ! fill lookup table body
      DO j=1, look_Z_snow%nTa

        IF (impipar_lookupgen == 2) THEN
          work_pe = j
          work_pe = round_robin(work_pe-1, pe_start, pe_end)
        ELSE
          work_pe = my_cart_id_fwo
        END IF

        IF (my_cart_id_fwo == work_pe) THEN
          
          DO i=1, look_Z_snow%nqi
            
            IF (luse_tmatrix .AND. ldebug_dbz) THEN
              ! Computation takes long, so print some progress information
              CALL progress_information (pos_=[i,j], &
                   shape_=[look_Z_snow%nqi,look_Z_snow%nTa], &
                   stride_=[look_Z_snow%nqi,1], cident='Dry snow lookup creation')
            END IF
            
            ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
            ! .. Suitable value for dqi, which is smaller than the local grid spacing:
            IF (i > 1) THEN
              ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
              !  take the dqi for differentiation as a fraction of the "lower" diff (to lower neighbour):
              dqi = 0.01_dp * (look_Z_snow%q_i_lin(i)-look_Z_snow%q_i_lin(i-1))
            ELSE
              ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
              !  to avoid negative q_i_lin below:
              dqi = 0.01_dp * (look_Z_snow%q_i_lin(i+1)-look_Z_snow%q_i_lin(i))
              dqi = MIN(dqi, 0.05_dp * look_Z_snow%q_i_lin(i))
            END IF

            DO kk=1, 3

              SELECT CASE (kk)
              CASE (1)
                ! calculate temperature dependent value of n0_s
                CALL calc_n0_snow(isnow_n0temp,look_Z_snow%T_a(j),&
                                  look_Z_snow%q_i_lin(i),snow%a_geo,n0_s)
                mgd = mgd_1mom(snow,look_Z_snow%q_i_lin(i),n0_s)
              CASE (2)
                ! calculate temperature dependent value of n0_s
                CALL calc_n0_snow(isnow_n0temp,look_Z_snow%T_a(j),&
                                  look_Z_snow%q_i_lin(i)-dqi,snow%a_geo,n0_s)
                mgd = mgd_1mom(snow,look_Z_snow%q_i_lin(i)-dqi,n0_s)
              CASE (3)
                ! calculate temperature dependent value of n0_s
                CALL calc_n0_snow(isnow_n0temp,look_Z_snow%T_a(j),&
                                  look_Z_snow%q_i_lin(i)+dqi,snow%a_geo,n0_s)
                mgd = mgd_1mom(snow,look_Z_snow%q_i_lin(i)+dqi,n0_s)
              END SELECT


              CALL zradar_snow_mie_vec(mgd, m_i(j), 0.5d0, lambda_radar, &
                                       luse_tmatrix, ldo_nonsphere, pMP, &
                                       Zfac, snow, Dmin, Dmax, &
                                       zh, ah, zv, rrhv, irhv, kdp, adp, zvh, &
                                       mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                                       mixingrulestring_core, matrixstring_core, inclusionstring_core, &
                                       llookupgen_mode=.TRUE.)

              SELECT CASE (kk)
              CASE (1)
                ! Store parameters in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_snow % zh   % val(i,j,1) = zh
                look_Z_snow % ah   % val(i,j,1) = ah
                look_Z_snow % zv   % val(i,j,1) = zv
                look_Z_snow % rrhv % val(i,j,1) = rrhv
                look_Z_snow % irhv % val(i,j,1) = irhv
                look_Z_snow % kdp  % val(i,j,1) = kdp
                look_Z_snow % adp  % val(i,j,1) = adp
                look_Z_snow % zvh  % val(i,j,1) = zvh
              CASE (2)
                zhu   = zh
                ahu   = ah
                zvu   = zv
                rrhvu = rrhv
                irhvu = irhv
                kdpu  = kdp
                adpu  = adp
                zvhu  = zvh
              CASE (3)
                ! Store derivatives in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_snow % zh   % dval(i,j,1) = (zh-zhu)     / dqi * 0.5_dp
                look_Z_snow % ah   % dval(i,j,1) = (ah-ahu)     / dqi * 0.5_dp
                look_Z_snow % zv   % dval(i,j,1) = (zv-zvu)     / dqi * 0.5_dp
                look_Z_snow % rrhv % dval(i,j,1) = (rrhv-rrhvu) / dqi * 0.5_dp
                look_Z_snow % irhv % dval(i,j,1) = (irhv-irhvu) / dqi * 0.5_dp
                look_Z_snow % kdp  % dval(i,j,1) = (kdp-kdpu)   / dqi * 0.5_dp
                look_Z_snow % adp  % dval(i,j,1) = (adp-adpu)   / dqi * 0.5_dp
                look_Z_snow % zvh  % dval(i,j,1) = (zvh-zvhu)   / dqi * 0.5_dp
              END SELECT

            END DO
            
          END DO

        END IF
      END DO

      look_Z_snow%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_snow%magicnr,isnow_n0temp,nTa,nqs,luse_tmatrix,ldo_nonsphere,&
      !     lambda_radar,Dmin,Dmax,snow%a_geo,snow%b_geo,Tmeltbegin,&
      !     mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      !     mixingrulestring_core,matrixstring_core,inclusionstring_core,' '//TRIM(look_Z_snow%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_snow, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_snow, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_snow%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_snow%nparams
      ! compute df(p(qi(qi_scaled)))/dqi_scaled = df/dp * dp/dqi * dqi/dqi_scaled = df/dp * dval * dqi/dqi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_snow%dmem(:,:,:,k) =  scale_dval_lut(look_Z_snow%mem(:,:,:,k), look_Z_snow%dmem(:,:,:,k), &
           look_Z_snow%qmem(:,k), &
           look_Z_snow%flag_qi_scal(k), look_Z_snow%f_scal_qi(k), &
           look_Z_snow%flag_mem_scal(k), look_Z_snow%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_snow%mem(:,:,:,k) =  scale_val_lut(look_Z_snow%mem(:,:,:,k), &
           look_Z_snow%flag_mem_scal(k), look_Z_snow%f_scal_mem(k))
    END DO
    
  END SUBROUTINE zradar_snow_1mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_graupel_1mom_lookupcreate(&
       look_Z_graupel,tableprops,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,&
       Zfac,Dmin,Dmax,graupel,&
       Tmeltbegin,mixingrulestring,matrixstring,inclusionstring,&
       impipar_lookupgen,pe_start,pe_end,&
       linterp_mode_dualpol,&
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_graupel
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    TYPE(particle), INTENT(in)            :: graupel
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP
    INTEGER, INTENT(in)                   :: unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, lambda_radar, Zfac, Dmin, Dmax
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring, &
                                             savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1)
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dqi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_graupel_1mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_graupel, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = 0, &                 ! table T_m will have 1 node
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = mi_Tmin + T0C_fwo, & ! K
         &                    Taup  = Tmeltbegin, &        ! K
         &                    Tmlow = T0C_fwo, &           ! for dry graupel, this is a dummy
         &                    Tmup  = T0C_fwo, &           ! for dry graupel, this is a dummy
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qg
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zh,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%ah,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zv,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%rrhv, lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%irhv, lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%kdp,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%adp,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zvh,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zh,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%ah,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zv,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%rrhv, lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%irhv, lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%kdp,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%adp,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zvh,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_graupel%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_1mom_drygraupel_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_graupel%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_graupel%mem  = 0.0_dp
      look_Z_graupel%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT

    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for dry graupel, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_graupel%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_graupel%nTa)

      ! fill lookup table body
      DO j=1, look_Z_graupel%nTa

        IF (impipar_lookupgen == 2) THEN
          work_pe = j
          work_pe = round_robin(work_pe-1, pe_start, pe_end)
        ELSE
          work_pe = my_cart_id_fwo
        END IF

        IF (my_cart_id_fwo == work_pe) THEN

          DO i=1, look_Z_graupel%nqi
            
            IF (luse_tmatrix .AND. ldebug_dbz) THEN
              ! Computation takes long, so print some progress information
              CALL progress_information (pos_=[i,j], &
                   shape_=[look_Z_graupel%nqi,look_Z_graupel%nTa], &
                   stride_=[look_Z_graupel%nqi,1], cident='Dry graupel lookup creation')
            END IF
            
            ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
            ! .. Suitable value for dqi, which is smaller than the local grid spacing:
            IF (i > 1) THEN
              ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
              !  take the dqi for differentiation as a fraction of the "lower" diff (to lower neighbour):
              dqi = 0.01_dp * (look_Z_graupel%q_i_lin(i)-look_Z_graupel%q_i_lin(i-1))
            ELSE
              ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
              !  to avoid negative q_i_lin below:
              dqi = 0.01_dp * (look_Z_graupel%q_i_lin(i+1)-look_Z_graupel%q_i_lin(i))
              dqi = MIN(dqi, 0.05_dp * look_Z_graupel%q_i_lin(i))
            END IF

            DO kk=1, 3

              SELECT CASE (kk)
              CASE (1)
                mgd = mgd_1mom(graupel,look_Z_graupel%q_i_lin(i),graupel%n0_const)
              CASE (2)
                mgd = mgd_1mom(graupel,look_Z_graupel%q_i_lin(i)-dqi,graupel%n0_const)
              CASE (3)
                mgd = mgd_1mom(graupel,look_Z_graupel%q_i_lin(i)+dqi,graupel%n0_const)
              END SELECT

              CALL zradar_graupel_mie_vec(mgd,m_i(j),lambda_radar,&
                   luse_tmatrix,ldo_nonsphere,pMP,&
                   Zfac,graupel,Dmin,Dmax,&
                   zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                   mixingrulestring,matrixstring,inclusionstring,&
                   llookupgen_mode=.TRUE.)

              SELECT CASE (kk)
              CASE (1)
                ! Store parameters in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_graupel % zh   % val(i,j,1) = zh
                look_Z_graupel % ah   % val(i,j,1) = ah
                look_Z_graupel % zv   % val(i,j,1) = zv
                look_Z_graupel % rrhv % val(i,j,1) = rrhv
                look_Z_graupel % irhv % val(i,j,1) = irhv
                look_Z_graupel % kdp  % val(i,j,1) = kdp
                look_Z_graupel % adp  % val(i,j,1) = adp
                look_Z_graupel % zvh  % val(i,j,1) = zvh
              CASE (2)
                zhu   = zh
                ahu   = ah
                zvu   = zv
                rrhvu = rrhv
                irhvu = irhv
                kdpu  = kdp
                adpu  = adp
                zvhu  = zvh
              CASE (3)
                ! Store derivatives in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_graupel % zh   % dval(i,j,1) = (zh-zhu)     / dqi * 0.5_dp
                look_Z_graupel % ah   % dval(i,j,1) = (ah-ahu)     / dqi * 0.5_dp
                look_Z_graupel % zv   % dval(i,j,1) = (zv-zvu)     / dqi * 0.5_dp
                look_Z_graupel % rrhv % dval(i,j,1) = (rrhv-rrhvu) / dqi * 0.5_dp
                look_Z_graupel % irhv % dval(i,j,1) = (irhv-irhvu) / dqi * 0.5_dp
                look_Z_graupel % kdp  % dval(i,j,1) = (kdp-kdpu)   / dqi * 0.5_dp
                look_Z_graupel % adp  % dval(i,j,1) = (adp-adpu)   / dqi * 0.5_dp
                look_Z_graupel % zvh  % dval(i,j,1) = (zvh-zvhu)   / dqi * 0.5_dp
              END SELECT

            END DO
            
          END DO
          
        END IF
      END DO

      look_Z_graupel%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_graupel%magicnr,nTa,nqg,luse_tmatrix,ldo_nonsphere,&
      !     lambda_radar,Dmin,Dmax,&
      !     graupel%n0_const,graupel%a_geo,graupel%b_geo,Tmeltbegin,&
      !     mixingrulestring,matrixstring,inclusionstring,' '//TRIM(look_Z_graupel%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_graupel, yzroutine, TRIM(hashtext))


    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_graupel, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_graupel%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_graupel%nparams
      ! compute df(p(qi(qi_scaled)))/dqi_scaled = df/dp * dp/dqi * dqi/dqi_scaled = df/dp * dval * dqi/dqi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_graupel%dmem(:,:,:,k) =  scale_dval_lut(look_Z_graupel%mem(:,:,:,k), look_Z_graupel%dmem(:,:,:,k), &
           look_Z_graupel%qmem(:,k), &
           look_Z_graupel%flag_qi_scal(k), look_Z_graupel%f_scal_qi(k), &
           look_Z_graupel%flag_mem_scal(k), look_Z_graupel%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_graupel%mem(:,:,:,k) =  scale_val_lut(look_Z_graupel%mem(:,:,:,k), &
           look_Z_graupel%flag_mem_scal(k), look_Z_graupel%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_graupel_1mom_lookupcreate

  !*******************************************************************************
  !
  ! Create Lookup-table vectors for radar reflectivity and its extinction for
  ! melting hydrometeors using Mie theory and 1-moment-scheme as a function of
  ! temperature (T_a), specific content of snow (q_i=q_s) or graupel (q_i=q_g)
  ! and maximum temperature (Tmax) where all particles are assumed as totally melted.
  !
  !*******************************************************************************

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_meltice_1mom_lookupcreate(&
      look_Z_meltice,tableprops,&
      itype_Dref_fmelt,&
      Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,&
      lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
      Zfac,Dmin,Dmax,ice,&
      mixingrulestring,matrixstring,inclusionstring,&
      hoststring,hostmatrixstring,hostinclusionstring,&
      impipar_lookupgen,pe_start,pe_end,&
      linterp_mode_dualpol,&
      savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_meltice
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    TYPE(particle), INTENT(in)            :: ice
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP, pMPr
    INTEGER, INTENT(in)                   :: itype_Dref_fmelt, unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, meltdegTmin, Tmax_min, Tmax_max, &
                                             lambda_radar, Zfac, &
                                             Dmin, Dmax
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring, &
                                             hoststring, hostmatrixstring, hostinclusionstring, &
                                             savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1), m_w(tableprops%nTa+1)
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dqi, n_i
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_meltice_1mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_meltice, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = tableprops%nTm, &    ! table T_m will have nTm+1 nodes
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = Tmeltbegin,               &  ! K
         &                    Taup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    Tmlow = MAX(Tmeltbegin,Tmax_min), &  ! K
         &                    Tmup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qi
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zh,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%ah,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zv,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%rrhv, lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%irhv, lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%kdp,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%adp,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zvh,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zh,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%ah,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zv,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%rrhv, lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%irhv, lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%kdp,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%adp,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zvh,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_meltice%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_1mom_meltice_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_meltice%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_meltice%mem   = 0.0_dp
      look_Z_meltice%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for melting ice, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_w = m_complex_water_ray_vec(DBLE(lambda_radar), &
                                    MIN(MAX(look_Z_meltice%T_a-T0C_fwo,mw_Tmin),mw_Tmax), &
                                    look_Z_meltice%nTa)
      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_meltice%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_meltice%nTa)

      ! fill lookup table body
      DO k=1, look_Z_meltice%nTm
        DO j=1, look_Z_meltice%nTa

          IF (impipar_lookupgen == 2) THEN
            ! Parallelize lookup table nodes over nTa and nTm  (j and k) and feed the different q_i
            !  nodes to the same PE at constant T_a and T_m. Then, the Mie routine below will
            !  compute the scattering amplitudes for each particle size D_s only once and
            !  re-use it in the integration over size distribution instead of re-computing the same
            !  thing over and over.
            CALL sub2ind2D (j, k, look_Z_meltice%nTa, work_pe)
            work_pe = round_robin(work_pe-1, pe_start, pe_end)
          ELSE
            work_pe = my_cart_id_fwo
          END IF

          IF (my_cart_id_fwo == work_pe) THEN

            n_i = nice_mono_1mom(look_Z_meltice%T_a(j))

            DO i=1, look_Z_meltice%nqi

              IF (luse_tmatrix .AND. ldebug_dbz) THEN
                ! Computation takes long, so print some progress information
                CALL progress_information (pos_=[i,j,k], &
                     shape_=[look_Z_meltice%nqi,look_Z_meltice%nTa,look_Z_meltice%nTm], &
                     stride_=[look_Z_meltice%nqi,1,1], cident='Melting ice lookup creation')
              END IF

              ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
              ! .. Suitable value for dqi, which is smaller than the local grid spacing:
              IF (i > 1) THEN
                ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
                !  take the dqi for differentiation as a fraction of the "lower" diff (to lower neighbour):
                dqi = 0.01_dp * (look_Z_meltice%q_i_lin(i)-look_Z_meltice%q_i_lin(i-1))
              ELSE
                ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
                !  to avoid negative q_i_lin below:
                dqi = 0.01_dp * (look_Z_meltice%q_i_lin(i+1)-look_Z_meltice%q_i_lin(i))
                dqi = MIN(dqi, 0.05_dp * look_Z_meltice%q_i_lin(i))
              END IF

              DO kk=1, 3

                ! Mimic a monodisperse size distribution around a mean mass of q_i/n_i by using 2-mom mgd method and
                !  choosing very large values for ice%mu and ice%nu (in init_1mom_types(), radar_interface.f90):
                SELECT CASE (kk)
                CASE (1)
                  mgd = mgd_2mom(ice,look_Z_meltice%q_i_lin(i),n_i)
                CASE (2)
                  mgd = mgd_2mom(ice,look_Z_meltice%q_i_lin(i)-dqi,n_i)
                CASE (3)
                  mgd = mgd_2mom(ice,look_Z_meltice%q_i_lin(i)+dqi,n_i)
                END SELECT
              
                CALL zradar_wetice_mie_vec(mgd,look_Z_meltice%T_a(j),&
                                           m_i(j),m_w(j),&
                                           itype_Dref_fmelt,Tmeltbegin,meltdegTmin,&
                                           lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                                           Zfac,ice,rain,Dmin,Dmax,&
                                           zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                                           mixingrulestring,matrixstring,inclusionstring,&
                                           hoststring,hostmatrixstring,hostinclusionstring,&
                                           look_Z_meltice%T_m(k))

                SELECT CASE (kk)
                CASE (1)
                  ! Store parameters in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltice % zh   % val(i,j,k) = zh
                  look_Z_meltice % ah   % val(i,j,k) = ah
                  look_Z_meltice % zv   % val(i,j,k) = zv
                  look_Z_meltice % rrhv % val(i,j,k) = rrhv
                  look_Z_meltice % irhv % val(i,j,k) = irhv
                  look_Z_meltice % kdp  % val(i,j,k) = kdp
                  look_Z_meltice % adp  % val(i,j,k) = adp
                  look_Z_meltice % zvh  % val(i,j,k) = zvh
                CASE (2)
                  zhu   = zh
                  ahu   = ah
                  zvu   = zv
                  rrhvu = rrhv
                  irhvu = irhv
                  kdpu  = kdp
                  adpu  = adp
                  zvhu  = zvh
                CASE (3)
                  ! Store derivatives in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltice % zh   % dval(i,j,k) = (zh-zhu)     / dqi * 0.5_dp
                  look_Z_meltice % ah   % dval(i,j,k) = (ah-ahu)     / dqi * 0.5_dp
                  look_Z_meltice % zv   % dval(i,j,k) = (zv-zvu)     / dqi * 0.5_dp
                  look_Z_meltice % rrhv % dval(i,j,k) = (rrhv-rrhvu) / dqi * 0.5_dp
                  look_Z_meltice % irhv % dval(i,j,k) = (irhv-irhvu) / dqi * 0.5_dp
                  look_Z_meltice % kdp  % dval(i,j,k) = (kdp-kdpu)   / dqi * 0.5_dp
                  look_Z_meltice % adp  % dval(i,j,k) = (adp-adpu)   / dqi * 0.5_dp
                  look_Z_meltice % zvh  % dval(i,j,k) = (zvh-zvhu)   / dqi * 0.5_dp
                END SELECT
                
              END DO

            END DO
              
          END IF
          
        END DO
      END DO

      look_Z_meltice%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_meltice%magicnr,nTa,nqi,nTm,luse_tmatrix,ldo_nonsphere,&
      !     itype_Dref_fmelt,Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,lambda_radar,&
      !     Dmin,Dmax,ice%a_geo,ice%b_geo,&
      !     mixingrulestring,matrixstring,inclusionstring,&
      !     hoststring,hostmatrixstring,hostinclusionstring,' '//TRIM(look_Z_meltice%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_meltice, yzroutine, TRIM(hashtext))


    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_meltice, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_meltice%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_meltice%nparams
      ! compute df(p(qi(qi_scaled)))/dqi_scaled = df/dp * dp/dqi * dqi/dqi_scaled = df/dp * dval * dqi/dqi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_meltice%dmem(:,:,:,k) =  scale_dval_lut(look_Z_meltice%mem(:,:,:,k), look_Z_meltice%dmem(:,:,:,k), &
           look_Z_meltice%qmem(:,k), &
           look_Z_meltice%flag_qi_scal(k), look_Z_meltice%f_scal_qi(k), &
           look_Z_meltice%flag_mem_scal(k), look_Z_meltice%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_meltice%mem(:,:,:,k) =  scale_val_lut(look_Z_meltice%mem(:,:,:,k), &
           look_Z_meltice%flag_mem_scal(k), look_Z_meltice%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_meltice_1mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_meltsnow_1mom_lookupcreate(&
       look_Z_meltsnow,tableprops,isnow_n0temp,&
       itype_Dref_fmelt,&
       Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
       Zfac,Dmin,Dmax,snow,&
       mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
       mixingrulestring_core,matrixstring_core,inclusionstring_core,&
       hoststring_shell,hostmatrixstring_shell,hostinclusionstring_shell,&
       hoststring_core,hostmatrixstring_core,hostinclusionstring_core,&
       impipar_lookupgen,pe_start,pe_end,&
       linterp_mode_dualpol,&
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_meltsnow
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    TYPE(particle), INTENT(in)            :: snow
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP, pMPr
    INTEGER, INTENT(in)                   :: isnow_n0temp, itype_Dref_fmelt, unitnr
    REAL(KIND=dp),  INTENT(in)            :: Tmeltbegin, meltdegTmin, Tmax_min, Tmax_max, &
                                             lambda_radar, Zfac, &
                                             Dmin, Dmax
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                                             mixingrulestring_core, matrixstring_core, inclusionstring_core, &
                                             hoststring_shell, hostmatrixstring_shell, hostinclusionstring_shell, &
                                             hoststring_core, hostmatrixstring_core, hostinclusionstring_core, &
                                             savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1), m_w(tableprops%nTa+1)
    REAL(KIND=dp)       :: n0_s
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dqi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_meltsnow_1mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_meltsnow, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = tableprops%nTm, &    ! table T_m will have nTm+1 nodes
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = Tmeltbegin,               &  ! K
         &                    Taup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    Tmlow = MAX(Tmeltbegin,Tmax_min), &  ! K
         &                    Tmup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qr
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zh,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%ah,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zv,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%rrhv, lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%irhv, lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%kdp,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%adp,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zvh,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zh,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%ah,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zv,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%rrhv, lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%irhv, lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%kdp,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%adp,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zvh,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_meltsnow%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_1mom_meltsnow_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_meltsnow%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_meltsnow%mem   = 0.0_dp
      look_Z_meltsnow%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for melting snow, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_w = m_complex_water_ray_vec(DBLE(lambda_radar), &
                                    MIN(MAX(look_Z_meltsnow%T_a-T0C_fwo,mw_Tmin),mw_Tmax), &
                                    look_Z_meltsnow%nTa)
      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_meltsnow%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_meltsnow%nTa)

      ! fill lookup table body
      DO k=1, look_Z_meltsnow%nTm
        DO j=1, look_Z_meltsnow%nTa

          IF (impipar_lookupgen == 2) THEN
            ! Parallelize lookup table nodes over nTa and nTm  (j and k) and feed the different q_i
            !  nodes to the same PE at constant T_a and T_m. Then, the Mie routine below will
            !  compute the scattering amplitudes for each particle size D_s only once and
            !  re-use it in the integration over size distribution instead of re-computing the same
            !  thing over and over.
            CALL sub2ind2D (j, k, look_Z_meltsnow%nTa, work_pe)
            work_pe = round_robin(work_pe-1, pe_start, pe_end)
          ELSE
            work_pe = my_cart_id_fwo
          END IF

          IF (my_cart_id_fwo == work_pe) THEN

            DO i=1, look_Z_meltsnow%nqi

              IF (luse_tmatrix .AND. ldebug_dbz) THEN
                ! Computation takes long, so print some progress information
                CALL progress_information (pos_=[i,j,k], &
                     shape_=[look_Z_meltsnow%nqi,look_Z_meltsnow%nTa,look_Z_meltsnow%nTm], &
                     stride_=[look_Z_meltsnow%nqi,1,1], cident='Melting snow lookup creation')
              END IF
            
              ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
              ! .. Suitable value for dqi, which is smaller than the local grid spacing:
              IF (i > 1) THEN
                ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
                !  take the dqi for differentiation as a fraction of the "lower" diff (to lower neighbour):
                dqi = 0.01_dp * (look_Z_meltsnow%q_i_lin(i)-look_Z_meltsnow%q_i_lin(i-1))
              ELSE
                ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
                !  to avoid negative q_i_lin below:
                dqi = 0.01_dp * (look_Z_meltsnow%q_i_lin(i+1)-look_Z_meltsnow%q_i_lin(i))
                dqi = MIN(dqi, 0.05_dp * look_Z_meltsnow%q_i_lin(i))
              END IF


              DO kk=1, 3
                
                ! calculate temperature dependent value of n0_s
                SELECT CASE (kk)
                CASE (1)
                  CALL calc_n0_snow(isnow_n0temp, look_Z_meltsnow%T_a(j), &
                       look_Z_meltsnow%q_i_lin(i), snow%a_geo, n0_s)
                  mgd = mgd_1mom(snow,look_Z_meltsnow%q_i_lin(i),n0_s)
                CASE (2)
                  CALL calc_n0_snow(isnow_n0temp, look_Z_meltsnow%T_a(j), &
                       look_Z_meltsnow%q_i_lin(i)-dqi, snow%a_geo, n0_s)
                  mgd = mgd_1mom(snow,look_Z_meltsnow%q_i_lin(i)-dqi,n0_s)
                CASE (3)
                  CALL calc_n0_snow(isnow_n0temp, look_Z_meltsnow%T_a(j), &
                       look_Z_meltsnow%q_i_lin(i)+dqi, snow%a_geo, n0_s)
                  mgd = mgd_1mom(snow,look_Z_meltsnow%q_i_lin(i)+dqi,n0_s)
                END SELECT
              
                ! 0.5d0 = radienverh, Tmeltbegin = Tmin_f, meltdegTmin = 0d0
                CALL zradar_wetsnow_mie_vec(mgd, look_Z_meltsnow%T_a(j),&
                     m_i(j), m_w(j),&
                     0.5d0, itype_Dref_fmelt, Tmeltbegin, meltdegTmin, lambda_radar, &
                     luse_tmatrix, ldo_nonsphere, pMP, pMPr, &
                     Zfac, snow, rain, Dmin, Dmax, &
                     zh, ah, zv, rrhv, irhv, kdp, adp, zvh, &
                     mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                     mixingrulestring_core, matrixstring_core, inclusionstring_core, &
                     hoststring_shell, hostmatrixstring_shell, hostinclusionstring_shell, &
                     hoststring_core, hostmatrixstring_core, hostinclusionstring_core, &
                     look_Z_meltsnow%T_m(k),llookupgen_mode=.TRUE.)

                SELECT CASE (kk)
                CASE (1)
                  ! Store parameters in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltsnow % zh   % val(i,j,k) = zh
                  look_Z_meltsnow % ah   % val(i,j,k) = ah
                  look_Z_meltsnow % zv   % val(i,j,k) = zv
                  look_Z_meltsnow % rrhv % val(i,j,k) = rrhv
                  look_Z_meltsnow % irhv % val(i,j,k) = irhv
                  look_Z_meltsnow % kdp  % val(i,j,k) = kdp
                  look_Z_meltsnow % adp  % val(i,j,k) = adp
                  look_Z_meltsnow % zvh  % val(i,j,k) = zvh
                CASE (2)
                  zhu   = zh
                  ahu   = ah
                  zvu   = zv
                  rrhvu = rrhv
                  irhvu = irhv
                  kdpu  = kdp
                  adpu  = adp
                  zvhu  = zvh
                CASE (3)
                  ! Store derivatives in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltsnow % zh   % dval(i,j,k) = (zh-zhu)     / dqi * 0.5_dp
                  look_Z_meltsnow % ah   % dval(i,j,k) = (ah-ahu)     / dqi * 0.5_dp
                  look_Z_meltsnow % zv   % dval(i,j,k) = (zv-zvu)     / dqi * 0.5_dp
                  look_Z_meltsnow % rrhv % dval(i,j,k) = (rrhv-rrhvu) / dqi * 0.5_dp
                  look_Z_meltsnow % irhv % dval(i,j,k) = (irhv-irhvu) / dqi * 0.5_dp
                  look_Z_meltsnow % kdp  % dval(i,j,k) = (kdp-kdpu)   / dqi * 0.5_dp
                  look_Z_meltsnow % adp  % dval(i,j,k) = (adp-adpu)   / dqi * 0.5_dp
                  look_Z_meltsnow % zvh  % dval(i,j,k) = (zvh-zvhu)   / dqi * 0.5_dp
                END SELECT
                
              END DO

            END DO
              
          END IF
          
        END DO
      END DO

      look_Z_meltsnow%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_meltsnow%magicnr,isnow_n0temp,nTa,nqs,nTm,luse_tmatrix,ldo_nonsphere, &
      !     itype_Dref_fmelt,Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,&
      !     lambda_radar,Dmin,Dmax,snow%a_geo,snow%b_geo, &
      !     mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      !     mixingrulestring_core,matrixstring_core,inclusionstring_core,&
      !     hoststring_shell,hostmatrixstring_shell,hostinclusionstring_shell,&
      !     hoststring_core,hostmatrixstring_core,hostinclusionstring_core,&
      !     ' '//TRIM(look_Z_meltsnow%chydroconfig)
      
      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_meltsnow, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_meltsnow, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_meltsnow%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_meltsnow%nparams
      ! compute df(p(qi(qi_scaled)))/dqi_scaled = df/dp * dp/dqi * dqi/dqi_scaled = df/dp * dval * dqi/dqi_scaled,
      !  where df/dp may depend on val itself.
       look_Z_meltsnow%dmem(:,:,:,k) =  scale_dval_lut(look_Z_meltsnow%mem(:,:,:,k), look_Z_meltsnow%dmem(:,:,:,k), &
           look_Z_meltsnow%qmem(:,k), &
           look_Z_meltsnow%flag_qi_scal(k), look_Z_meltsnow%f_scal_qi(k), &
           look_Z_meltsnow%flag_mem_scal(k), look_Z_meltsnow%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_meltsnow%mem(:,:,:,k) =  scale_val_lut(look_Z_meltsnow%mem(:,:,:,k), &
           look_Z_meltsnow%flag_mem_scal(k), look_Z_meltsnow%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_meltsnow_1mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_meltgraupel_1mom_lookupcreate(&
       look_Z_meltgraupel,tableprops,igraupel_type,&
       itype_Dref_fmelt,&
       Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
       Zfac,Dmin,Dmax,graupel,&
       mixingrulestring,matrixstring,inclusionstring,&
       mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
       hoststring,hostmatrixstring,hostinclusionstring,&
       impipar_lookupgen,pe_start,pe_end,&
       linterp_mode_dualpol,&
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_meltgraupel
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    TYPE(particle), INTENT(in)            :: graupel
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP, pMPr
    INTEGER, INTENT(in)                   :: igraupel_type, itype_Dref_fmelt, unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, meltdegTmin, Tmax_min, Tmax_max, &
                                             lambda_radar, Zfac, &
                                             Dmin, Dmax
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring, &
                                             mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                                             hoststring, hostmatrixstring, hostinclusionstring, &
                                             savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1), m_w(tableprops%nTa+1)
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dqi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_meltgraupel_1mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_meltgraupel, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = tableprops%nTm, &    ! table T_m will have nTm+1 nodes
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = Tmeltbegin,               &  ! K
         &                    Taup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    Tmlow = MAX(Tmeltbegin,Tmax_min), &  ! K
         &                    Tmup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qi
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zh,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%ah,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zv,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%rrhv, lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%irhv, lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%kdp,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%adp,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zvh,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zh,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%ah,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zv,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%rrhv, lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%irhv, lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%kdp,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%adp,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zvh,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_meltgraupel%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_1mom_meltgraupel_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_meltgraupel%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_meltgraupel%mem   = 0.0_dp
      look_Z_meltgraupel%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for melting graupel, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_w = m_complex_water_ray_vec(DBLE(lambda_radar), &
                                    MIN(MAX(look_Z_meltgraupel%T_a-T0C_fwo,mw_Tmin),mw_Tmax), &
                                    look_Z_meltgraupel%nTa)
      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_meltgraupel%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_meltgraupel%nTa)

      ! fill lookup table body
      DO k=1, look_Z_meltgraupel%nTm
        DO j=1, look_Z_meltgraupel%nTa

          IF (impipar_lookupgen == 2) THEN
            ! Parallelize lookup table nodes over nTa and nTm  (j and k) and feed the different q_i
            !  nodes to the same PE at constant T_a and T_m. Then, the Mie routine below will
            !  compute the scattering amplitudes for each particle size D_s only once and
            !  re-use it in the integration over size distribution instead of re-computing the same
            !  thing over and over.
            CALL sub2ind2D (j, k, look_Z_meltgraupel%nTa, work_pe)
            work_pe = round_robin(work_pe-1, pe_start, pe_end)
          ELSE
            work_pe = my_cart_id_fwo
          END IF

          IF (my_cart_id_fwo == work_pe) THEN

            DO i=1, look_Z_meltgraupel%nqi

              IF (luse_tmatrix .AND. ldebug_dbz) THEN
                ! Computation takes long, so print some progress information
                CALL progress_information (pos_=[i,j,k], &
                     shape_=[look_Z_meltgraupel%nqi,look_Z_meltgraupel%nTa,look_Z_meltgraupel%nTm], &
                     stride_=[look_Z_meltgraupel%nqi,1,1], cident='Melting graupel lookup creation')
              END IF

              ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
              ! .. Suitable value for dqi, which is smaller than the local grid spacing:
              IF (i > 1) THEN
                ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
                !  take the dqi for differentiation as a fraction of the "lower" diff (to lower neighbour):
                dqi = 0.01_dp * (look_Z_meltgraupel%q_i_lin(i)-look_Z_meltgraupel%q_i_lin(i-1))
              ELSE
                ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
                !  to avoid negative q_i_lin below:
                dqi = 0.01_dp * (look_Z_meltgraupel%q_i_lin(i+1)-look_Z_meltgraupel%q_i_lin(i))
                dqi = MIN(dqi, 0.05_dp * look_Z_meltgraupel%q_i_lin(i))
              END IF

              DO kk=1, 3

                SELECT CASE (kk)
                CASE (1)
                  mgd = mgd_1mom(graupel,look_Z_meltgraupel%q_i_lin(i),graupel%n0_const)
                CASE (2)
                  mgd = mgd_1mom(graupel,look_Z_meltgraupel%q_i_lin(i)-dqi,graupel%n0_const)
                CASE (3)
                  mgd = mgd_1mom(graupel,look_Z_meltgraupel%q_i_lin(i)+dqi,graupel%n0_const)
                END SELECT
              
                ! Tmeltbegin = 263.16d0, meltdegTmin = 0.2d0
                SELECT CASE (igraupel_type)
                CASE(1)
                  CALL zradar_wetgr_mie_vec(mgd,look_Z_meltgraupel%T_a(j),&
                                            m_i(j),m_w(j),&
                                            itype_Dref_fmelt,Tmeltbegin,meltdegTmin,&
                                            lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                                            Zfac,graupel,rain,Dmin,Dmax,&
                                            zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                                            mixingrulestring,matrixstring,inclusionstring,&
                                            hoststring,hostmatrixstring,hostinclusionstring,&
                                            look_Z_meltgraupel%T_m(k))
                CASE(2)
                  CALL zradar_wetgr_twosph_mie_vec(mgd,look_Z_meltgraupel%T_a(j),&
                                            m_i(j),m_w(j),&
                                            itype_Dref_fmelt,Tmeltbegin,meltdegTmin,&
                                            lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                                            Zfac,graupel,rain,Dmin,Dmax,&
                                            zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                                            mixingrulestring,matrixstring,inclusionstring,&
                                            mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
                                            look_Z_meltgraupel%T_m(k))
                CASE(3)
                  CALL zradar_wetgr_wsph_mie_vec(mgd,look_Z_meltgraupel%T_a(j),&
                                            m_i(j),m_w(j),&
                                            itype_Dref_fmelt,Tmeltbegin,meltdegTmin,&
                                            lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                                            Zfac,graupel,rain,Dmin,Dmax,&
                                            zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                                            mixingrulestring,matrixstring,inclusionstring,&
                                            look_Z_meltgraupel%T_m(k))
                END SELECT

                SELECT CASE (kk)
                CASE (1)
                  ! Store parameters in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltgraupel % zh   % val(i,j,k) = zh
                  look_Z_meltgraupel % ah   % val(i,j,k) = ah
                  look_Z_meltgraupel % zv   % val(i,j,k) = zv
                  look_Z_meltgraupel % rrhv % val(i,j,k) = rrhv
                  look_Z_meltgraupel % irhv % val(i,j,k) = irhv
                  look_Z_meltgraupel % kdp  % val(i,j,k) = kdp
                  look_Z_meltgraupel % adp  % val(i,j,k) = adp
                  look_Z_meltgraupel % zvh  % val(i,j,k) = zvh
                CASE (2)
                  zhu   = zh
                  ahu   = ah
                  zvu   = zv
                  rrhvu = rrhv
                  irhvu = irhv
                  kdpu  = kdp
                  adpu  = adp
                  zvhu  = zvh
                CASE (3)
                  ! Store derivatives in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltgraupel % zh   % dval(i,j,k) = (zh-zhu)     / dqi * 0.5_dp
                  look_Z_meltgraupel % ah   % dval(i,j,k) = (ah-ahu)     / dqi * 0.5_dp
                  look_Z_meltgraupel % zv   % dval(i,j,k) = (zv-zvu)     / dqi * 0.5_dp
                  look_Z_meltgraupel % rrhv % dval(i,j,k) = (rrhv-rrhvu) / dqi * 0.5_dp
                  look_Z_meltgraupel % irhv % dval(i,j,k) = (irhv-irhvu) / dqi * 0.5_dp
                  look_Z_meltgraupel % kdp  % dval(i,j,k) = (kdp-kdpu)   / dqi * 0.5_dp
                  look_Z_meltgraupel % adp  % dval(i,j,k) = (adp-adpu)   / dqi * 0.5_dp
                  look_Z_meltgraupel % zvh  % dval(i,j,k) = (zvh-zvhu)   / dqi * 0.5_dp
                END SELECT
                
              END DO

            END DO
              
          END IF
          
        END DO
      END DO

      look_Z_meltgraupel%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_meltgraupel%magicnr,igraupel_type,nTa,nqg,nTm,luse_tmatrix,ldo_nonsphere,&
      !     itype_Dref_fmelt,Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,lambda_radar,&
      !     Dmin,Dmax,graupel%n0_const,graupel%a_geo,graupel%b_geo,&
      !     mixingrulestring,matrixstring,inclusionstring,&
      !     mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      !     hoststring,hostmatrixstring,hostinclusionstring,' '//TRIM(look_Z_meltgraupel%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_meltgraupel, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_meltgraupel, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_meltgraupel%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_meltgraupel%nparams
      ! compute df(p(qi(qi_scaled)))/dqi_scaled = df/dp * dp/dqi * dqi/dqi_scaled = df/dp * dval * dqi/dqi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_meltgraupel%dmem(:,:,:,k) =  scale_dval_lut(look_Z_meltgraupel%mem(:,:,:,k), look_Z_meltgraupel%dmem(:,:,:,k), &
           look_Z_meltgraupel%qmem(:,k), &
           look_Z_meltgraupel%flag_qi_scal(k), look_Z_meltgraupel%f_scal_qi(k), &
           look_Z_meltgraupel%flag_mem_scal(k), look_Z_meltgraupel%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_meltgraupel%mem(:,:,:,k) =  scale_val_lut(look_Z_meltgraupel%mem(:,:,:,k), &
           look_Z_meltgraupel%flag_mem_scal(k), look_Z_meltgraupel%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_meltgraupel_1mom_lookupcreate

#ifdef TWOMOM_SB

  !*******************************************************************************
  !
  ! Create Lookup-table vectors for radar reflectivity and its extinction
  ! for dry hydrometeors using Mie theory and 2-moment-scheme as a function
  ! of temperature (T_a) and specific content (q_i) of rain water (q_i=q_r),
  ! dry snow (q_i=q_s) and dry graupel (q_i=q_g)
  !
  !*******************************************************************************

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_rain_2mom_lookupcreate(&
       look_Z_rain,tableprops,do_muD,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,&
       Zfac,Dmin,Dmax,rain,&
       impipar_lookupgen,pe_start,pe_end,&
       linterp_mode_dualpol,&
       savepath_read,savepath_write,unitnr,hashtext)

    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_rain
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP
    INTEGER, INTENT(in)                   :: unitnr
    REAL(KIND=dp), INTENT(in)             :: lambda_radar, Zfac, Dmin, Dmax
    TYPE(particle), INTENT(in)            :: rain
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: savepath_read, savepath_write, hashtext
    LOGICAL, INTENT(in)                   :: do_muD

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_w(tableprops%nTa+1)
    REAL(KIND=dp)       :: q_r, nuD, x_r, x_ru, x_ro
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dxi
    INTEGER             :: i, j, k, kk, nk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_rain_2mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cmudstr, cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, is_for_mud_relation, print_debug

    ! .. hypothetical fixed rain number concentration for computing
    !     the lookup table in terms of x_r with a subroutine that
    !     takes as input q_r and n_r:
    REAL(KIND=dp)    :: n_r


    ! compute n_r in such a way that all q_r >= q_crit_radar and at the
    ! same time n_r >= n_crit_radar, so that reflectivity is guaranteed to
    ! be computed for all points of the table:
    n_r = MAX(1d-6 / tableprops%qilow, n_crit_radar)

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    cmudstr(:) = ' '
    IF (do_muD) THEN
      is_for_mud_relation = .TRUE.
      cmudstr = '-mud'
      CALL init_dbzlookuptable (lut   = look_Z_rain, &
           &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
           &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
           &                    ntm   = tableprops%nTm, &    ! table T_m will have nmuD+1 nodes
           &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
           &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
           &                    Talow = mw_Tmin + T0C_fwo, & ! K
           &                    Taup  = mw_Tmax + T0C_fwo, & ! K
           &                    Tmlow = 0d0, &               ! for rain, this is the range for mu in the mu(D_m)-relation
           &                    Tmup  = 31d0, &              ! for rain, this is the range for mu in the mu(D_m)-relation
           &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qr
           &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
           &                    ierr = ierr)
     
      nk = look_Z_rain%nTm
      nuD = 1d0
    ELSE
      is_for_mud_relation = .FALSE.
      cmudstr = '-mufixed'
      CALL init_dbzlookuptable (lut   = look_Z_rain, &
           &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
           &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
           &                    ntm   = 0, &                 ! table T_m will have 1 node
           &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
           &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
           &                    Talow = mw_Tmin + T0C_fwo, & ! K
           &                    Taup  = mw_Tmax + T0C_fwo, & ! K
           &                    Tmlow = 0d0, &               ! for rain with fixed mu, this is a dummy
           &                    Tmup  = 0d0, &               ! for rain with fixed mu, this is a dummy
           &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qr
           &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
           &                    ierr = ierr)

      nk = 1
    END IF

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_rain%zh,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%ah,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%zv,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%rrhv, lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%irhv, lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%kdp,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%adp,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%zvh,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_rain%zh,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_rain%ah,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_rain%zv,   lut=look_Z_rain, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_rain%rrhv, lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%irhv, lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%kdp,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%adp,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_rain%zvh,  lut=look_Z_rain, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_rain%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_2mom_rain'//TRIM(cmudstr)//'_'//&
               TRIM(ADJUSTL(cmagicnr))//'_'//TRIM(look_Z_rain%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_rain%mem   = 0.0_dp
      look_Z_rain%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for rain, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_w = m_complex_water_ray_vec(DBLE(lambda_radar), &
                                    MIN(MAX(look_Z_rain%T_a-T0C_fwo,mw_Tmin),mw_Tmax), &
                                    look_Z_rain%nTa)

      ! fill lookup table body
      DO k=1, nk
        DO j=1, look_Z_rain%nTa
          
          IF (impipar_lookupgen == 2) THEN
            CALL sub2ind2D (j, k, look_Z_rain%nTa, work_pe)
            work_pe = round_robin(work_pe-1, pe_start, pe_end)
          ELSE
            work_pe = my_cart_id_fwo
          END IF

          IF (my_cart_id_fwo == work_pe) THEN

            DO i=1, look_Z_rain%nqi

              IF (luse_tmatrix .AND. ldebug_dbz) THEN
                ! Computation takes long, so print some progress information
                CALL progress_information (pos_=[i,j,k], &
                     shape_=[look_Z_rain%nqi,look_Z_rain%nTa,nk], &
                     stride_=[look_Z_rain%nqi,1,1], cident='Rain lookup creation')
              END IF
            
              ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
              ! .. Suitable value for dxi, which is smaller than the local grid spacing:
              IF (i > 1) THEN
                ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
                !  take the dxi for differentiation as a fraction of the "lower" diff (to lower neighbour):
                x_r = look_Z_rain%q_i_lin(i)
                x_ru = look_Z_rain%q_i_lin(i-1)
                dxi = 0.01_dp * (x_r - x_ru)
              ELSE
                ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
                !  to avoid negative q_i_lin below:
                x_r = look_Z_rain%q_i_lin(i)
                x_ro = look_Z_rain%q_i_lin(i+1)
                dxi = 0.01_dp * (x_ro - x_r)
                dxi = MIN(dxi, 0.05_dp * x_r)
              END IF

              DO kk=1, 3
                
                SELECT CASE (kk)
                CASE (1)
                  q_r = look_Z_rain%q_i_lin(i) * n_r
                CASE (2)
                  q_r = (look_Z_rain%q_i_lin(i) - dxi) * n_r
                CASE (3)
                  q_r = (look_Z_rain%q_i_lin(i) + dxi) * n_r
                END SELECT
                
                IF (is_for_mud_relation) THEN
                  mgd = mgd_2mom(rain,q_r,n_r,look_Z_rain%T_m(k),nuD)
                ELSE
                  mgd = mgd_2mom(rain,q_r,n_r)
                END IF

                CALL zradar_rain_mie_vec(mgd,m_w(j),lambda_radar,&
                                         luse_tmatrix,ldo_nonsphere,pMP,&
                                         Zfac,Dmin,Dmax,&
                                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                                         llookupgen_mode=.TRUE.)

                
                SELECT CASE (kk)
                CASE (1)
                  ! Store parameters in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_rain % zh   % val(i,j,k) = zh   / n_r
                  look_Z_rain % ah   % val(i,j,k) = ah   / n_r
                  look_Z_rain % zv   % val(i,j,k) = zv   / n_r
                  look_Z_rain % rrhv % val(i,j,k) = rrhv / n_r
                  look_Z_rain % irhv % val(i,j,k) = irhv / n_r
                  look_Z_rain % kdp  % val(i,j,k) = kdp  / n_r
                  look_Z_rain % adp  % val(i,j,k) = adp  / n_r
                  look_Z_rain % zvh  % val(i,j,k) = zvh  / n_r
                CASE (2)
                  zhu   = zh
                  ahu   = ah
                  zvu   = zv
                  rrhvu = rrhv
                  irhvu = irhv
                  kdpu  = kdp
                  adpu  = adp
                  zvhu  = zvh
                CASE (3)
                  ! Store derivatives in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_rain % zh   % dval(i,j,k) = (zh-zhu)     / dxi * 0.5_dp / n_r
                  look_Z_rain % ah   % dval(i,j,k) = (ah-ahu)     / dxi * 0.5_dp / n_r
                  look_Z_rain % zv   % dval(i,j,k) = (zv-zvu)     / dxi * 0.5_dp / n_r
                  look_Z_rain % rrhv % dval(i,j,k) = (rrhv-rrhvu) / dxi * 0.5_dp / n_r
                  look_Z_rain % irhv % dval(i,j,k) = (irhv-irhvu) / dxi * 0.5_dp / n_r
                  look_Z_rain % kdp  % dval(i,j,k) = (kdp-kdpu)   / dxi * 0.5_dp / n_r
                  look_Z_rain % adp  % dval(i,j,k) = (adp-adpu)   / dxi * 0.5_dp / n_r
                  look_Z_rain % zvh  % dval(i,j,k) = (zvh-zvhu)   / dxi * 0.5_dp / n_r
                END SELECT
                
              END DO

            END DO

          END IF
          
        END DO
      END DO

      look_Z_rain%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !IF (is_for_mud_relation) THEN
      !  WRITE(extrastr,*) look_Z_rain%magicnr,nTa,nxr,luse_tmatrix,ldo_nonsphere,&
      !       nmuD,lambda_radar,Dmin,Dmax,' '//TRIM(look_Z_rain%chydroconfig)
      !ELSE
      !  WRITE(extrastr,*) look_Z_rain%magicnr,nTa,nxr,luse_tmatrix,ldo_nonsphere,&
      !       lambda_radar,Dmin,Dmax,' '//TRIM(look_Z_rain%chydroconfig)
      !END IF

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_rain, yzroutine, TRIM(hashtext))

    ELSE
      
      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: reading lookup table for rain, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_rain, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_rain%is_initialized = .TRUE.
      END IF

    END IF

    
    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_rain%nparams
      ! compute df(p(xi(xi_scaled)))/dxi_scaled = df/dp * dp/dxi * dxi/dxi_scaled = df/dp * dval * dxi/dxi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_rain%dmem(:,:,:,k) =  scale_dval_lut(look_Z_rain%mem(:,:,:,k), look_Z_rain%dmem(:,:,:,k), &
           look_Z_rain%qmem(:,k), &
           look_Z_rain%flag_qi_scal(k), look_Z_rain%f_scal_qi(k), &
           look_Z_rain%flag_mem_scal(k), look_Z_rain%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_rain%mem(:,:,:,k) =  scale_val_lut(look_Z_rain%mem(:,:,:,k), &
           look_Z_rain%flag_mem_scal(k), look_Z_rain%f_scal_mem(k))
    END DO

    
  END SUBROUTINE zradar_rain_2mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_ice_2mom_lookupcreate(&
       look_Z_ice,tableprops,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,&
       Zfac,Dmin,Dmax,ice,&
       Tmeltbegin,mixingrulestring,matrixstring,inclusionstring,&
       impipar_lookupgen,pe_start,pe_end,&
       linterp_mode_dualpol,&
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_ice
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP
    INTEGER, INTENT(in)                   :: unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, lambda_radar, Zfac, Dmin, Dmax
    CLASS(particle), INTENT(in)           :: ice
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1)
    REAL(KIND=dp)       :: q_i, x_i, x_iu, x_io
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dxi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_ice_2mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! .. hypothetical fixed snow number concentration for computing
    !     the lookup table in terms of x_g with a subroutine that
    !     takes as input q_i and n_i:
    REAL(KIND=dp)    :: n_i

    ! compute n_i in such a way that all q_i >= q_crit_radar and at the
    ! same time n_i >= n_crit_radar, so that reflectivity is guaranteed to
    ! be computed for all points of the table:
    !n_i = MAX(q_crit_radar%ice / ice%x_min, n_crit_radar%ice)
    n_i = MAX(1d-6 / tableprops%qilow, n_crit_radar)

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_ice, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = 0, &                 ! table T_m will have 1 node
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = mi_Tmin + T0C_fwo, & ! K
         &                    Taup  = Tmeltbegin, &        ! K
         &                    Tmlow = T0C_fwo, &           ! for dry ice, this is a dummy
         &                    Tmup  = T0C_fwo, &           ! for dry ice, this is a dummy
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qi
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_ice%zh,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%ah,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%zv,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%rrhv, lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%irhv, lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%kdp,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%adp,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%zvh,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_ice%zh,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_ice%ah,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_ice%zv,   lut=look_Z_ice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_ice%rrhv, lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%irhv, lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%kdp,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%adp,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_ice%zvh,  lut=look_Z_ice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_ice%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_2mom_dryice_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_ice%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_ice%mem   = 0.0_dp
      look_Z_ice%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for dry ice, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_ice%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_ice%nTa)

      ! fill lookup table body
      DO j=1, look_Z_ice%nTa

        IF (impipar_lookupgen == 2) THEN
          work_pe = j
          work_pe = round_robin(work_pe-1, pe_start, pe_end)
        ELSE
          work_pe = my_cart_id_fwo
        END IF

        IF (my_cart_id_fwo == work_pe) THEN
          
          DO i=1, look_Z_ice%nqi

            IF (luse_tmatrix .AND. ldebug_dbz) THEN
              ! Computation takes long, so print some progress information
              CALL progress_information (pos_=[i,j], &
                   shape_=[look_Z_ice%nqi,look_Z_ice%nTa], &
                   stride_=[look_Z_ice%nqi,1], cident='Dry ice lookup creation')
            END IF
            
            ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
            ! .. Suitable value for dxi, which is smaller than the local grid spacing:
            IF (i > 1) THEN
              ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
              !  take the dxi for differentiation as a fraction of the "lower" diff (to lower neighbour):
              x_i = look_Z_ice%q_i_lin(i)
              x_iu = look_Z_ice%q_i_lin(i-1)
              dxi = 0.01_dp * (x_i - x_iu)
            ELSE
              ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
              !  to avoid negative q_i_lin below:
              x_i = look_Z_ice%q_i_lin(i)
              x_io = look_Z_ice%q_i_lin(i+1)
              dxi = 0.01_dp * (x_io - x_i)
              dxi = MIN(dxi, 0.05_dp * x_i)
            END IF

            DO kk=1, 3
                
              SELECT CASE (kk)
              CASE (1)
                q_i = look_Z_ice%q_i_lin(i) * n_i
              CASE (2)
                q_i = (look_Z_ice%q_i_lin(i) - dxi) * n_i
              CASE (3)
                q_i = (look_Z_ice%q_i_lin(i) + dxi) * n_i
              END SELECT

              mgd = mgd_2mom(ice,q_i,n_i)
            
              CALL zradar_ice_mie_vec(&
                      mgd,m_i(j),lambda_radar,&
                      luse_tmatrix,ldo_nonsphere,pMP,&
                      Zfac,ice,Dmin,Dmax,&
                      zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                      mixingrulestring,matrixstring,inclusionstring,&
                      llookupgen_mode=.TRUE.)

              SELECT CASE (kk)
              CASE (1)
                ! Store parameters in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_ice % zh   % val(i,j,1) = zh   / n_i
                look_Z_ice % ah   % val(i,j,1) = ah   / n_i
                look_Z_ice % zv   % val(i,j,1) = zv   / n_i
                look_Z_ice % rrhv % val(i,j,1) = rrhv / n_i
                look_Z_ice % irhv % val(i,j,1) = irhv / n_i
                look_Z_ice % kdp  % val(i,j,1) = kdp  / n_i
                look_Z_ice % adp  % val(i,j,1) = adp  / n_i
                look_Z_ice % zvh  % val(i,j,1) = zvh  / n_i
              CASE (2)
                zhu   = zh
                ahu   = ah
                zvu   = zv
                rrhvu = rrhv
                irhvu = irhv
                kdpu  = kdp
                adpu  = adp
                zvhu  = zvh
              CASE (3)
                ! Store derivatives in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_ice % zh   % dval(i,j,1) = (zh-zhu)     / dxi * 0.5_dp / n_i
                look_Z_ice % ah   % dval(i,j,1) = (ah-ahu)     / dxi * 0.5_dp / n_i
                look_Z_ice % zv   % dval(i,j,1) = (zv-zvu)     / dxi * 0.5_dp / n_i
                look_Z_ice % rrhv % dval(i,j,1) = (rrhv-rrhvu) / dxi * 0.5_dp / n_i
                look_Z_ice % irhv % dval(i,j,1) = (irhv-irhvu) / dxi * 0.5_dp / n_i
                look_Z_ice % kdp  % dval(i,j,1) = (kdp-kdpu)   / dxi * 0.5_dp / n_i
                look_Z_ice % adp  % dval(i,j,1) = (adp-adpu)   / dxi * 0.5_dp / n_i
                look_Z_ice % zvh  % dval(i,j,1) = (zvh-zvhu)   / dxi * 0.5_dp / n_i
              END SELECT
                
            END DO

          END DO
          
        END IF
        
      END DO

      look_Z_ice%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,'(I12,1X,2(I4,1X),2(L1,1X),4(ES14.6E2,1X),3(A12,1X),A)') &
      !     look_Z_ice%magicnr,&
      !     nTa,nxi,luse_tmatrix,ldo_nonsphere,&
      !     Tmeltbegin,lambda_radar,Dmin,Dmax,&
      !     mixingrulestring,matrixstring,inclusionstring,&
      !     TRIM(look_Z_ice%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_ice, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_ice, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_ice%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_ice%nparams
      ! compute df(p(xi(xi_scaled)))/dxi_scaled = df/dp * dp/dxi * dxi/dxi_scaled = df/dp * dval * dxi/dxi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_ice%dmem(:,:,:,k) =  scale_dval_lut(look_Z_ice%mem(:,:,:,k), look_Z_ice%dmem(:,:,:,k), &
           look_Z_ice%qmem(:,k), &
           look_Z_ice%flag_qi_scal(k), look_Z_ice%f_scal_qi(k), &
           look_Z_ice%flag_mem_scal(k), look_Z_ice%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_ice%mem(:,:,:,k) =  scale_val_lut(look_Z_ice%mem(:,:,:,k), &
           look_Z_ice%flag_mem_scal(k), look_Z_ice%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_ice_2mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_snow_2mom_lookupcreate(&
       look_Z_snow,tableprops,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,&
       Zfac,Dmin,Dmax,snow,Tmeltbegin,&
       mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
       mixingrulestring_core,matrixstring_core,inclusionstring_core,&
       impipar_lookupgen, pe_start, pe_end,    &
       linterp_mode_dualpol,                   &
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_snow
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP
    INTEGER, INTENT(in)                   :: unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, lambda_radar, Zfac, Dmin, Dmax
    CLASS(particle), INTENT(in)           :: snow
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                                             mixingrulestring_core, matrixstring_core, inclusionstring_core
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1)
    REAL(KIND=dp)       :: q_s, x_s, x_su, x_so
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dxi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_snow_2mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! .. hypothetical fixed snow number concentration for computing
    !     the lookup table in terms of x_s with a subroutine that
    !     takes as input q_s and n_s:
    REAL(KIND=dp)    :: n_s

    ! compute n_s in such a way that all q_s >= q_crit_radar and at the
    ! same time n_r >= n_crit_radar, so that reflectivity is guaranteed to
    ! be computed for all points of the table:
    !n_s = MAX(q_crit_radar%snow / snow%x_min, n_crit_radar%snow)
    n_s = MAX(1d-6 / tableprops%qilow, n_crit_radar)

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_snow, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = 0, &                 ! table T_m will have 1 node
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = mi_Tmin + T0C_fwo, & ! K
         &                    Taup  = Tmeltbegin, &        ! K
         &                    Tmlow = T0C_fwo, &           ! for dry snow, this is a dummy
         &                    Tmup  = T0C_fwo, &           ! for dry snow, this is a dummy
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qr
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_snow%zh,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%ah,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%zv,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%rrhv, lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%irhv, lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%kdp,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%adp,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%zvh,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_snow%zh,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_snow%ah,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_snow%zv,   lut=look_Z_snow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_snow%rrhv, lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%irhv, lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%kdp,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%adp,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_snow%zvh,  lut=look_Z_snow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_snow%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_2mom_drysnow_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_snow%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_snow%mem   = 0.0_dp
      look_Z_snow%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for dry snow, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_snow%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_snow%nTa)

      ! fill lookup table body
      DO j=1, look_Z_snow%nTa

        IF (impipar_lookupgen == 2) THEN
          work_pe = j
          work_pe = round_robin(work_pe-1, pe_start, pe_end)
        ELSE
          work_pe = my_cart_id_fwo
        END IF

        IF (my_cart_id_fwo == work_pe) THEN

          DO i=1, look_Z_snow%nqi

            IF (luse_tmatrix .AND. ldebug_dbz) THEN
              ! Computation takes long, so print some progress information
              CALL progress_information (pos_=[i,j], &
                   shape_=[look_Z_snow%nqi,look_Z_snow%nTa], &
                   stride_=[look_Z_snow%nqi,1], cident='Dry snow lookup creation')
            END IF

            ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
            ! .. Suitable value for dxi, which is smaller than the local grid spacing:
            IF (i > 1) THEN
              ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
              !  take the dxi for differentiation as a fraction of the "lower" diff (to lower neighbour):
              x_s = look_Z_snow%q_i_lin(i)
              x_su = look_Z_snow%q_i_lin(i-1)
              dxi = 0.01_dp * (x_s - x_su)
            ELSE
              ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
              !  to avoid negative q_i_lin below:
              x_s = look_Z_snow%q_i_lin(i)
              x_so = look_Z_snow%q_i_lin(i+1)
              dxi = 0.01_dp * (x_so - x_s)
              dxi = MIN(dxi, 0.05_dp * x_s)
            END IF

            DO kk=1, 3
                
              SELECT CASE (kk)
              CASE (1)
                q_s = look_Z_snow%q_i_lin(i) * n_s
              CASE (2)
                q_s = (look_Z_snow%q_i_lin(i) - dxi) * n_s
              CASE (3)
                q_s = (look_Z_snow%q_i_lin(i) + dxi) * n_s
              END SELECT

              mgd = mgd_2mom(snow,q_s,n_s)

              CALL zradar_snow_mie_vec(mgd, m_i(j), 0.5d0, lambda_radar, &
                                       luse_tmatrix, ldo_nonsphere, pMP, &
                                       Zfac, snow, Dmin, Dmax, &
                                       zh, ah, zv, rrhv, irhv, kdp, adp, zvh, &
                                       mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                                       mixingrulestring_core, matrixstring_core, inclusionstring_core,&
                                       llookupgen_mode=.TRUE.)

              SELECT CASE (kk)
              CASE (1)
                ! Store parameters in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_snow % zh   % val(i,j,1) = zh   / n_s
                look_Z_snow % ah   % val(i,j,1) = ah   / n_s
                look_Z_snow % zv   % val(i,j,1) = zv   / n_s
                look_Z_snow % rrhv % val(i,j,1) = rrhv / n_s
                look_Z_snow % irhv % val(i,j,1) = irhv / n_s
                look_Z_snow % kdp  % val(i,j,1) = kdp  / n_s
                look_Z_snow % adp  % val(i,j,1) = adp  / n_s
                look_Z_snow % zvh  % val(i,j,1) = zvh  / n_s
              CASE (2)
                zhu   = zh
                ahu   = ah
                zvu   = zv
                rrhvu = rrhv
                irhvu = irhv
                kdpu  = kdp
                adpu  = adp
                zvhu  = zvh
              CASE (3)
                ! Store derivatives in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_snow % zh   % dval(i,j,1) = (zh-zhu)     / dxi * 0.5_dp / n_s
                look_Z_snow % ah   % dval(i,j,1) = (ah-ahu)     / dxi * 0.5_dp / n_s
                look_Z_snow % zv   % dval(i,j,1) = (zv-zvu)     / dxi * 0.5_dp / n_s
                look_Z_snow % rrhv % dval(i,j,1) = (rrhv-rrhvu) / dxi * 0.5_dp / n_s
                look_Z_snow % irhv % dval(i,j,1) = (irhv-irhvu) / dxi * 0.5_dp / n_s
                look_Z_snow % kdp  % dval(i,j,1) = (kdp-kdpu)   / dxi * 0.5_dp / n_s
                look_Z_snow % adp  % dval(i,j,1) = (adp-adpu)   / dxi * 0.5_dp / n_s
                look_Z_snow % zvh  % dval(i,j,1) = (zvh-zvhu)   / dxi * 0.5_dp / n_s
              END SELECT
                
            END DO

          END DO
          
        END IF
        
      END DO

      look_Z_snow%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_snow%magicnr,nTa,nxs,luse_tmatrix,ldo_nonsphere,&
      !     Tmeltbegin,lambda_radar,Dmin,Dmax, &
      !     mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      !     mixingrulestring_core,matrixstring_core,inclusionstring_core,' '//TRIM(look_Z_snow%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_snow, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_snow, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_snow%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_snow%nparams
      ! compute df(p(xi(xi_scaled)))/dxi_scaled = df/dp * dp/dxi * dxi/dxi_scaled = df/dp * dval * dxi/dxi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_snow%dmem(:,:,:,k) =  scale_dval_lut(look_Z_snow%mem(:,:,:,k), look_Z_snow%dmem(:,:,:,k), &
           look_Z_snow%qmem(:,k), &
           look_Z_snow%flag_qi_scal(k), look_Z_snow%f_scal_qi(k), &
           look_Z_snow%flag_mem_scal(k), look_Z_snow%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_snow%mem(:,:,:,k) =  scale_val_lut(look_Z_snow%mem(:,:,:,k), &
           look_Z_snow%flag_mem_scal(k), look_Z_snow%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_snow_2mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_graupel_2mom_lookupcreate(&
       look_Z_graupel,tableprops,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,&
       Zfac,Dmin,Dmax,graupel,&
       Tmeltbegin,mixingrulestring,matrixstring,inclusionstring,&
       impipar_lookupgen, pe_start, pe_end,    &
       linterp_mode_dualpol,                   &
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_graupel
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP
    INTEGER, INTENT(in)                   :: unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, lambda_radar, Zfac, Dmin, Dmax
    CLASS(particle), INTENT(in)           :: graupel
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1)
    REAL(KIND=dp)       :: q_g, x_g, x_gu, x_go
    REAL(KIND=dp)       :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dxi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_graupel_2mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! .. hypothetical fixed snow number concentration for computing
    !     the lookup table in terms of x_g with a subroutine that
    !     takes as input q_g and n_g:
    REAL(KIND=dp)    :: n_g

    ! compute n_g in such a way that all q_g >= q_crit_radar and at the
    ! same time n_g >= n_crit_radar, so that reflectivity is guaranteed to
    ! be computed for all points of the table:
    !n_g = MAX(q_crit_radar%graupel / graupel%x_min, n_crit_radar%graupel)
    n_g = MAX(1d-6 / tableprops%qilow, n_crit_radar)

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_graupel, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = 0, &                 ! table T_m will have 1 node
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = mi_Tmin + T0C_fwo, & ! K
         &                    Taup  = Tmeltbegin, &        ! K
         &                    Tmlow = T0C_fwo, &           ! for dry graupel, this is a dummy
         &                    Tmup  = T0C_fwo, &           ! for dry graupel, this is a dummy
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qr
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zh,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%ah,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zv,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%rrhv, lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%irhv, lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%kdp,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%adp,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zvh,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zh,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%ah,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zv,   lut=look_Z_graupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%rrhv, lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%irhv, lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%kdp,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%adp,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_graupel%zvh,  lut=look_Z_graupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_graupel%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_2mom_drygraupel_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_graupel%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_graupel%mem   = 0.0_dp
      look_Z_graupel%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for dry graupel, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_graupel%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_graupel%nTa)

      ! fill lookup table body
      DO j=1, look_Z_graupel%nTa

        IF (impipar_lookupgen == 2) THEN
          work_pe = j
          work_pe = round_robin(work_pe-1, pe_start, pe_end)
        ELSE
          work_pe = my_cart_id_fwo
        END IF

        IF (my_cart_id_fwo == work_pe) THEN
          
          DO i=1, look_Z_graupel%nqi

            IF (luse_tmatrix .AND. ldebug_dbz) THEN
              ! Computation takes long, so print some progress information
              CALL progress_information (pos_=[i,j], &
                   shape_=[look_Z_graupel%nqi,look_Z_graupel%nTa], &
                   stride_=[look_Z_graupel%nqi,1], cident='Dry graupel lookup creation')
            END IF
            
            ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
            ! .. Suitable value for dxi, which is smaller than the local grid spacing:
            IF (i > 1) THEN
              ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
              !  take the dxi for differentiation as a fraction of the "lower" diff (to lower neighbour):
              x_g = look_Z_graupel%q_i_lin(i)
              x_gu = look_Z_graupel%q_i_lin(i-1)
              dxi = 0.01_dp * (x_g - x_gu)
            ELSE
              ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
              !  to avoid negative q_i_lin below:
              x_g = look_Z_graupel%q_i_lin(i)
              x_go = look_Z_graupel%q_i_lin(i+1)
              dxi = 0.01_dp * (x_go - x_g)
              dxi = MIN(dxi, 0.05_dp * x_g)
            END IF

            DO kk=1, 3
                
              SELECT CASE (kk)
              CASE (1)
                q_g = look_Z_graupel%q_i_lin(i) * n_g
              CASE (2)
                q_g = (look_Z_graupel%q_i_lin(i) - dxi) * n_g
              CASE (3)
                q_g = (look_Z_graupel%q_i_lin(i) + dxi) * n_g
              END SELECT

              mgd = mgd_2mom(graupel,q_g,n_g)
            
              CALL zradar_graupel_mie_vec(mgd,m_i(j),lambda_radar,&
                                          luse_tmatrix,ldo_nonsphere,pMP,&
                                          Zfac,graupel,Dmin,Dmax,&
                                          zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                                          mixingrulestring,matrixstring,inclusionstring,&
                                          llookupgen_mode=.TRUE.)

              SELECT CASE (kk)
              CASE (1)
                ! Store parameters in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_graupel % zh   % val(i,j,1) = zh   / n_g
                look_Z_graupel % ah   % val(i,j,1) = ah   / n_g
                look_Z_graupel % zv   % val(i,j,1) = zv   / n_g
                look_Z_graupel % rrhv % val(i,j,1) = rrhv / n_g
                look_Z_graupel % irhv % val(i,j,1) = irhv / n_g
                look_Z_graupel % kdp  % val(i,j,1) = kdp  / n_g
                look_Z_graupel % adp  % val(i,j,1) = adp  / n_g
                look_Z_graupel % zvh  % val(i,j,1) = zvh  / n_g
              CASE (2)
                zhu   = zh
                ahu   = ah
                zvu   = zv
                rrhvu = rrhv
                irhvu = irhv
                kdpu  = kdp
                adpu  = adp
                zvhu  = zvh
              CASE (3)
                ! Store derivatives in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_graupel % zh   % dval(i,j,1) = (zh-zhu)     / dxi * 0.5_dp / n_g
                look_Z_graupel % ah   % dval(i,j,1) = (ah-ahu)     / dxi * 0.5_dp / n_g
                look_Z_graupel % zv   % dval(i,j,1) = (zv-zvu)     / dxi * 0.5_dp / n_g
                look_Z_graupel % rrhv % dval(i,j,1) = (rrhv-rrhvu) / dxi * 0.5_dp / n_g
                look_Z_graupel % irhv % dval(i,j,1) = (irhv-irhvu) / dxi * 0.5_dp / n_g
                look_Z_graupel % kdp  % dval(i,j,1) = (kdp-kdpu)   / dxi * 0.5_dp / n_g
                look_Z_graupel % adp  % dval(i,j,1) = (adp-adpu)   / dxi * 0.5_dp / n_g
                look_Z_graupel % zvh  % dval(i,j,1) = (zvh-zvhu)   / dxi * 0.5_dp / n_g
              END SELECT
                
            END DO

          END DO
          
        END IF
        
      END DO

      look_Z_graupel%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_graupel%magicnr,nTa,nxg,luse_tmatrix,ldo_nonsphere,&
      !     Tmeltbegin,lambda_radar,Dmin,Dmax,&
      !     mixingrulestring,matrixstring,inclusionstring,' '//TRIM(look_Z_graupel%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_graupel, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_graupel, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_graupel%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_graupel%nparams
      ! compute df(p(xi(xi_scaled)))/dxi_scaled = df/dp * dp/dxi * dxi/dxi_scaled = df/dp * dval * dxi/dxi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_graupel%dmem(:,:,:,k) =  scale_dval_lut(look_Z_graupel%mem(:,:,:,k), look_Z_graupel%dmem(:,:,:,k), &
           look_Z_graupel%qmem(:,k), &
           look_Z_graupel%flag_qi_scal(k), look_Z_graupel%f_scal_qi(k), &
           look_Z_graupel%flag_mem_scal(k), look_Z_graupel%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_graupel%mem(:,:,:,k) =  scale_val_lut(look_Z_graupel%mem(:,:,:,k), &
           look_Z_graupel%flag_mem_scal(k), look_Z_graupel%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_graupel_2mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_hail_2mom_lookupcreate(&
       look_Z_hail,tableprops,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,&
       Zfac,Dmin,Dmax,hail,&
       Tmeltbegin,mixingrulestring,matrixstring,inclusionstring,&
       impipar_lookupgen, pe_start, pe_end,    &
       linterp_mode_dualpol,                   &
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_hail
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP
    INTEGER, INTENT(in)                   :: unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, lambda_radar, Zfac, Dmin, Dmax
    CLASS(particle), INTENT(in)           :: hail
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1)
    REAL(KIND=dp)       :: q_h, x_h, x_hu, x_ho
    REAL(KIND=dp)       :: zh, ah, zv, rrhv, irhv, kdp, adp, zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dxi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_hail_2mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! .. hypothetical fixed snow number concentration for computing
    !     the lookup table in terms of x_h with a subroutine that
    !     takes as input q_h and n_h:
    REAL(KIND=dp)    :: n_h

    ! compute n_h in such a way that all q_h >= q_crit_radar and at the
    ! same time n_h >= n_crit_radar, so that reflectivity is guaranteed to
    ! be computed for all points of the table:
    !n_h = MAX(q_crit_radar%hail / hail%x_min, n_crit_radar%hail)
    n_h = MAX(1d-6 / tableprops%qilow, n_crit_radar)

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_hail, &
         &                    nqi   = tableprops%nqi, &    ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &    ! table T_a will have nTa+1 nodes
         &                    ntm   = 0, &                 ! table T_m will have 1 node
         &                    qilow = tableprops%qilow, &  ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &   ! kg, linear scaling, not log10!
         &                    Talow = mi_Tmin + T0C_fwo, & ! K
         &                    Taup  = Tmeltbegin, &        ! K
         &                    Tmlow = T0C_fwo, &           ! for dry hail, this is a dummy
         &                    Tmup  = T0C_fwo, &           ! for dry hail, this is a dummy
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qr
         &                    f_eq  = 1.0_dp, &            ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_hail%zh,   lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%ah,   lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%zv,   lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%rrhv, lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%irhv, lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%kdp,  lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%adp,  lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%zvh,  lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_hail%zh,   lut=look_Z_hail, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_hail%ah,   lut=look_Z_hail, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_hail%zv,   lut=look_Z_hail, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_hail%rrhv, lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%irhv, lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%kdp,  lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%adp,  lut=look_Z_hail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_hail%zvh,  lut=look_Z_hail, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_hail%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_2mom_dryhail_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_hail%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_hail%mem   = 0.0_dp
      look_Z_hail%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for dry hail, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_hail%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_hail%nTa)

      ! fill lookup table body
      DO j=1, look_Z_hail%nTa

        IF (impipar_lookupgen == 2) THEN
          work_pe = j
          work_pe = round_robin(work_pe-1, pe_start, pe_end)
        ELSE
          work_pe = my_cart_id_fwo
        END IF

        IF (my_cart_id_fwo == work_pe) THEN

          DO i=1, look_Z_hail%nqi
            
            IF (luse_tmatrix .AND. ldebug_dbz) THEN
              ! Computation takes long, so print some progress information
              CALL progress_information (pos_=[i,j], &
                   shape_=[look_Z_hail%nqi,look_Z_hail%nTa], &
                   stride_=[look_Z_hail%nqi,1], cident='Dry hail lookup creation')
            END IF

            ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
            ! .. Suitable value for dxi, which is smaller than the local grid spacing:
            IF (i > 1) THEN
              ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
              !  take the dxi for differentiation as a fraction of the "lower" diff (to lower neighbour):
              x_h = look_Z_hail%q_i_lin(i)
              x_hu = look_Z_hail%q_i_lin(i-1)
              dxi = 0.01_dp * (x_h - x_hu)
            ELSE
              ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
              !  to avoid negative q_i_lin below:
              x_h = look_Z_hail%q_i_lin(i)
              x_ho = look_Z_hail%q_i_lin(i+1)
              dxi = 0.01_dp * (x_ho - x_h)
              dxi = MIN(dxi, 0.05_dp * x_h)
            END IF

            DO kk=1, 3
                
              SELECT CASE (kk)
              CASE (1)
                q_h = look_Z_hail%q_i_lin(i) * n_h
              CASE (2)
                q_h = (look_Z_hail%q_i_lin(i) - dxi) * n_h
              CASE (3)
                q_h = (look_Z_hail%q_i_lin(i) + dxi) * n_h
              END SELECT
              
              mgd = mgd_2mom(hail,q_h,n_h)

              CALL  zradar_hail_mie_vec(mgd,m_i(j),lambda_radar,&
                   luse_tmatrix,ldo_nonsphere,pMP,&
                   Zfac,hail,Dmin,Dmax,&
                   zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                   mixingrulestring,matrixstring,inclusionstring,&
                   llookupgen_mode=.TRUE.)

              SELECT CASE (kk)
              CASE (1)              
                ! Store parameters in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_hail % zh   % val(i,j,1) = zh   / n_h
                look_Z_hail % ah   % val(i,j,1) = ah   / n_h
                look_Z_hail % zv   % val(i,j,1) = zv   / n_h
                look_Z_hail % rrhv % val(i,j,1) = rrhv / n_h
                look_Z_hail % irhv % val(i,j,1) = irhv / n_h
                look_Z_hail % kdp  % val(i,j,1) = kdp  / n_h
                look_Z_hail % adp  % val(i,j,1) = adp  / n_h
                look_Z_hail % zvh  % val(i,j,1) = zvh  / n_h
              CASE (2)
                zhu   = zh
                ahu   = ah
                zvu   = zv
                rrhvu = rrhv
                irhvu = irhv
                kdpu  = kdp
                adpu  = adp
                zvhu  = zvh
              CASE (3)
                ! Store derivatives in original linear scaling in the tables. Transform to
                ! desired scaling after table writing/reading:
                look_Z_hail % zh   % dval(i,j,1) = (zh-zhu)     / dxi * 0.5_dp / n_h
                look_Z_hail % ah   % dval(i,j,1) = (ah-ahu)     / dxi * 0.5_dp / n_h
                look_Z_hail % zv   % dval(i,j,1) = (zv-zvu)     / dxi * 0.5_dp / n_h
                look_Z_hail % rrhv % dval(i,j,1) = (rrhv-rrhvu) / dxi * 0.5_dp / n_h
                look_Z_hail % irhv % dval(i,j,1) = (irhv-irhvu) / dxi * 0.5_dp / n_h
                look_Z_hail % kdp  % dval(i,j,1) = (kdp-kdpu)   / dxi * 0.5_dp / n_h
                look_Z_hail % adp  % dval(i,j,1) = (adp-adpu)   / dxi * 0.5_dp / n_h
                look_Z_hail % zvh  % dval(i,j,1) = (zvh-zvhu)   / dxi * 0.5_dp / n_h
              END SELECT
                
            END DO

          END DO
          
        END IF
        
      END DO

      look_Z_hail%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_hail%magicnr,nTa,nxh,luse_tmatrix,ldo_nonsphere,&
      !     Tmeltbegin,lambda_radar,Dmin,Dmax,&
      !     mixingrulestring,matrixstring,inclusionstring,' '//TRIM(look_Z_hail%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_hail, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_hail, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_hail%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_hail%nparams
      ! compute df(p(xi(xi_scaled)))/dxi_scaled = df/dp * dp/dxi * dxi/dxi_scaled = df/dp * dval * dxi/dxi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_hail%dmem(:,:,:,k) =  scale_dval_lut(look_Z_hail%mem(:,:,:,k), look_Z_hail%dmem(:,:,:,k), &
           look_Z_hail%qmem(:,k), &
           look_Z_hail%flag_qi_scal(k), look_Z_hail%f_scal_qi(k), &
           look_Z_hail%flag_mem_scal(k), look_Z_hail%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_hail%mem(:,:,:,k) =  scale_val_lut(look_Z_hail%mem(:,:,:,k), &
           look_Z_hail%flag_mem_scal(k), look_Z_hail%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_hail_2mom_lookupcreate


  !*******************************************************************************
  !
  ! Create Lookup-table vectors for radar reflectivity and its extinction for
  ! melting hydrometeors using Mie theory and 1-moment-scheme as a function of
  ! temperature (T_a), specific content of snow (q_i=q_s), graupel (q_i=q_g) or
  ! hail (q_i=q_h), and maximum temperature (Tmax) where all particles are
  ! assumed as totally melted.
  !
  !*******************************************************************************

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_meltice_2mom_lookupcreate(&
       look_Z_meltice,tableprops,&
       itype_Dref_fmelt,&
       Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
       Zfac,Dmin,Dmax,ice,&
       mixingrulestring,matrixstring,inclusionstring,&
       hoststring,hostmatrixstring,hostinclusionstring,&
       impipar_lookupgen, pe_start, pe_end,    &
       linterp_mode_dualpol,                   &
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_meltice
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP, pMPr
    INTEGER, INTENT(in)                   :: itype_Dref_fmelt, unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, meltdegTmin, Tmax_min, Tmax_max, &
                                             lambda_radar, Zfac, Dmin, Dmax
    CLASS(particle), INTENT(in)           :: ice
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring, &
                                             hoststring, hostmatrixstring, hostinclusionstring, &
                                             savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1), m_w(tableprops%nTa+1)
    REAL(KIND=dp)       :: q_i, x_i, x_iu, x_io
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dxi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_meltice_2mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! .. hypothetical fixed snow number concentration for computing
    !     the lookup table in terms of x_i with a subroutine that
    !     takes as input q_i and n_i:
    REAL(KIND=dp)    :: n_i

    ! compute n_i in such a way that all q_i >= q_crit_radar and at the
    ! same time n_i >= n_crit_radar, so that reflectivity is guaranteed to
    ! be computed for all points of the table:
    !n_i = MAX(q_crit_radar%ice / ice%x_min, n_crit_radar%ice)
    n_i = MAX(1d-6 / tableprops%qilow, n_crit_radar)

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_meltice, &
         &                    nqi   = tableprops%nqi, &   ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &   ! table T_a will have nTa+1 nodes
         &                    ntm   = tableprops%nTm, &   ! table T_m will have nTm+1 nodes
         &                    qilow = tableprops%qilow, & ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &  ! kg, linear scaling, not log10!
         &                    Talow = Tmeltbegin,               &  ! K
         &                    Taup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    Tmlow = MAX(Tmeltbegin,Tmax_min), &  ! K
         &                    Tmup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qi
         &                    f_eq  = 1.0_dp, &           ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zh,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%ah,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zv,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%rrhv, lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%irhv, lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%kdp,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%adp,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zvh,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zh,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%ah,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zv,   lut=look_Z_meltice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%rrhv, lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%irhv, lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%kdp,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%adp,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltice%zvh,  lut=look_Z_meltice, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_meltice%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_2mom_meltice_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_meltice%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_meltice%mem   = 0.0_dp
      look_Z_meltice%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for melting ice, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_w = m_complex_water_ray_vec(DBLE(lambda_radar), &
                                    MIN(MAX(look_Z_meltice%T_a-T0C_fwo,mw_Tmin),mw_Tmax), &
                                    look_Z_meltice%nTa)
      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_meltice%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_meltice%nTa)

      ! fill lookup table body
      DO k=1, look_Z_meltice%nTm
        DO j=1, look_Z_meltice%nTa

          IF (impipar_lookupgen == 2) THEN
            ! Parallelize lookup table nodes over nTa and nTm  (j and k) and feed the different q_i
            !  nodes to the same PE at constant T_a and T_m. Then, the Mie routine below will
            !  compute the scattering amplitudes for each particle size D_s only once and
            !  re-use it in the integration over size distribution instead of re-computing the same
            !  thing over and over.
            CALL sub2ind2D (j, k, look_Z_meltice%nTa, work_pe)
            work_pe = round_robin(work_pe-1, pe_start, pe_end)
          ELSE
            work_pe = my_cart_id_fwo
          END IF

          IF (my_cart_id_fwo == work_pe) THEN

            DO i=1, look_Z_meltice%nqi

              IF (luse_tmatrix .AND. ldebug_dbz) THEN
                ! Computation takes long, so print some progress information
                CALL progress_information (pos_=[i,j,k], &
                     shape_=[look_Z_meltice%nqi,look_Z_meltice%nTa,look_Z_meltice%nTm], &
                     stride_=[look_Z_meltice%nqi,1,1], cident='Melting ice lookup creation')
              END IF

              ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
              ! .. Suitable value for dxi, which is smaller than the local grid spacing:
              IF (i > 1) THEN
                ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
                !  take the dxi for differentiation as a fraction of the "lower" diff (to lower neighbour):
                x_i = look_Z_meltice%q_i_lin(i)
                x_iu = look_Z_meltice%q_i_lin(i-1)
                dxi = 0.01_dp * (x_i - x_iu)
              ELSE
                ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
                !  to avoid negative q_i_lin below:
                x_i = look_Z_meltice%q_i_lin(i)
                x_io = look_Z_meltice%q_i_lin(i+1)
                dxi = 0.01_dp * (x_io - x_i)
                dxi = MIN(dxi, 0.05_dp * x_i)
              END IF

              DO kk=1, 3
                
                SELECT CASE (kk)
                CASE (1)
                  q_i = look_Z_meltice%q_i_lin(i) * n_i
                CASE (2)
                  q_i = (look_Z_meltice%q_i_lin(i) - dxi) * n_i
                CASE (3)
                  q_i = (look_Z_meltice%q_i_lin(i) + dxi) * n_i
                END SELECT

                mgd = mgd_2mom(ice,q_i,n_i)

                CALL zradar_wetice_mie_vec(mgd,look_Z_meltice%T_a(j),&
                     m_i(j),m_w(j),&
                     itype_Dref_fmelt,Tmeltbegin,meltdegTmin,&
                     lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                     Zfac,ice,rain,Dmin,Dmax,&
                     zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                     mixingrulestring,matrixstring,inclusionstring,&
                     hoststring,hostmatrixstring,hostinclusionstring,&
                     look_Z_meltice%T_m(k),llookupgen_mode=.TRUE.)
                
                SELECT CASE (kk)
                CASE (1)
                  ! Store parameters in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltice % zh   % val(i,j,k) = zh   / n_i
                  look_Z_meltice % ah   % val(i,j,k) = ah   / n_i
                  look_Z_meltice % zv   % val(i,j,k) = zv   / n_i
                  look_Z_meltice % rrhv % val(i,j,k) = rrhv / n_i
                  look_Z_meltice % irhv % val(i,j,k) = irhv / n_i
                  look_Z_meltice % kdp  % val(i,j,k) = kdp  / n_i
                  look_Z_meltice % adp  % val(i,j,k) = adp  / n_i
                  look_Z_meltice % zvh  % val(i,j,k) = zvh  / n_i
                CASE (2)
                  zhu   = zh
                  ahu   = ah
                  zvu   = zv
                  rrhvu = rrhv
                  irhvu = irhv
                  kdpu  = kdp
                  adpu  = adp
                  zvhu  = zvh
                CASE (3)
                  ! Store derivatives in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltice % zh   % dval(i,j,k) = (zh-zhu)     / dxi * 0.5_dp / n_i
                  look_Z_meltice % ah   % dval(i,j,k) = (ah-ahu)     / dxi * 0.5_dp / n_i
                  look_Z_meltice % zv   % dval(i,j,k) = (zv-zvu)     / dxi * 0.5_dp / n_i
                  look_Z_meltice % rrhv % dval(i,j,k) = (rrhv-rrhvu) / dxi * 0.5_dp / n_i
                  look_Z_meltice % irhv % dval(i,j,k) = (irhv-irhvu) / dxi * 0.5_dp / n_i
                  look_Z_meltice % kdp  % dval(i,j,k) = (kdp-kdpu)   / dxi * 0.5_dp / n_i
                  look_Z_meltice % adp  % dval(i,j,k) = (adp-adpu)   / dxi * 0.5_dp / n_i
                  look_Z_meltice % zvh  % dval(i,j,k) = (zvh-zvhu)   / dxi * 0.5_dp / n_i
                END SELECT
                
              END DO

            END DO
            
          END IF

        END DO
      END DO

      look_Z_meltice%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,'(I12,1X,3(I4,1X),2(L1,1X),I4,1X,7(ES14.6E2,1X),6(A12,1X),A)') &
      !     look_Z_meltice%magicnr,&
      !     nTa,nxi,nTm,luse_tmatrix,ldo_nonsphere,&
      !     itype_Dref_fmelt,&
      !     Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,lambda_radar,Dmin,Dmax,&
      !     mixingrulestring,matrixstring,inclusionstring,&
      !     hoststring,hostmatrixstring,hostinclusionstring,&
      !     TRIM(look_Z_meltice%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_meltice, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_meltice, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_meltice%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_meltice%nparams
      ! compute df(p(xi(xi_scaled)))/dxi_scaled = df/dp * dp/dxi * dxi/dxi_scaled = df/dp * dval * dxi/dxi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_meltice%dmem(:,:,:,k) =  scale_dval_lut(look_Z_meltice%mem(:,:,:,k), look_Z_meltice%dmem(:,:,:,k), &
           look_Z_meltice%qmem(:,k), &
           look_Z_meltice%flag_qi_scal(k), look_Z_meltice%f_scal_qi(k), &
           look_Z_meltice%flag_mem_scal(k), look_Z_meltice%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_meltice%mem(:,:,:,k) =  scale_val_lut(look_Z_meltice%mem(:,:,:,k), &
           look_Z_meltice%flag_mem_scal(k), look_Z_meltice%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_meltice_2mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_meltsnow_2mom_lookupcreate(&
       look_Z_meltsnow,tableprops,&
       itype_Dref_fmelt,&
       Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
       Zfac,Dmin,Dmax,snow,&
       mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
       mixingrulestring_core,matrixstring_core,inclusionstring_core,&
       hoststring_shell,hostmatrixstring_shell,hostinclusionstring_shell,&
       hoststring_core,hostmatrixstring_core,hostinclusionstring_core,&
       impipar_lookupgen, pe_start, pe_end,    &
       linterp_mode_dualpol,                   &
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_meltsnow
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP, pMPr
    INTEGER, INTENT(in)                   :: itype_Dref_fmelt, unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, meltdegTmin, Tmax_min, Tmax_max, &
                                             lambda_radar, Zfac, Dmin, Dmax
    CLASS(particle), INTENT(in)           :: snow
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                                             mixingrulestring_core, matrixstring_core, inclusionstring_core, &
                                               hoststring_shell, hostmatrixstring_shell, hostinclusionstring_shell, &
                                           hoststring_core, hostmatrixstring_core, hostinclusionstring_core
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1), m_w(tableprops%nTa+1)
    REAL(KIND=dp)       :: q_s, x_s, x_su, x_so
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dxi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_meltsnow_2mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! .. hypothetical fixed snow number concentration for computing
    !     the lookup table in terms of x_s with a subroutine that
    !     takes as input q_s and n_s:
    REAL(KIND=dp)    :: n_s

    ! compute n_s in such a way that all q_s >= q_crit_radar and at the
    ! same time n_r >= n_crit_radar, so that reflectivity is guaranteed to
    ! be computed for all points of the table:
    !n_s = MAX(q_crit_radar%snow / snow%x_min, n_crit_radar%snow)
    n_s = MAX(1d-6 / tableprops%qilow, n_crit_radar)

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_meltsnow, &
         &                    nqi   = tableprops%nqi, &   ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &   ! table T_a will have nTa+1 nodes
         &                    ntm   = tableprops%nTm, &   ! table T_m will have nTm+1 nodes
         &                    qilow = tableprops%qilow, & ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &  ! kg, linear scaling, not log10!
         &                    Talow = Tmeltbegin,               &  ! K
         &                    Taup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    Tmlow = MAX(Tmeltbegin,Tmax_min), &  ! K
         &                    Tmup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qi
         &                    f_eq  = 1.0_dp, &           ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zh,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%ah,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zv,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%rrhv, lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%irhv, lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%kdp,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%adp,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zvh,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zh,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%ah,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zv,   lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%rrhv, lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%irhv, lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%kdp,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%adp,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltsnow%zvh,  lut=look_Z_meltsnow, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_meltsnow%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_2mom_meltsnow_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_meltsnow%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_meltsnow%mem   = 0.0_dp
      look_Z_meltsnow%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for melting snow, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_w = m_complex_water_ray_vec(DBLE(lambda_radar), &
                                    MIN(MAX(look_Z_meltsnow%T_a-T0C_fwo,mw_Tmin),mw_Tmax), &
                                    look_Z_meltsnow%nTa)
      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_meltsnow%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_meltsnow%nTa)

      ! fill lookup table body
      DO k=1, look_Z_meltsnow%nTm
        DO j=1, look_Z_meltsnow%nTa

          IF (impipar_lookupgen == 2) THEN
            ! Parallelize lookup table nodes over nTa and nTm  (j and k) and feed the different q_i
            !  nodes to the same PE at constant T_a and T_m. Then, the Mie routine below will
            !  compute the scattering amplitudes for each particle size D_s only once and
            !  re-use it in the integration over size distribution instead of re-computing the same
            !  thing over and over.
            CALL sub2ind2D (j, k, look_Z_meltsnow%nTa, work_pe)
            work_pe = round_robin(work_pe-1, pe_start, pe_end)
          ELSE
            work_pe = my_cart_id_fwo
          END IF

          IF (my_cart_id_fwo == work_pe) THEN
            
            DO i=1, look_Z_meltsnow%nqi

              IF (luse_tmatrix .AND. ldebug_dbz) THEN

                ! Computation takes long, so print some progress information
                CALL progress_information (pos_=[i,j,k], &
                     shape_=[look_Z_meltsnow%nqi,look_Z_meltsnow%nTa,look_Z_meltsnow%nTm], &
                     stride_=[look_Z_meltsnow%nqi,1,1], cident='Melting snow lookup creation')
              END IF
            
              ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
              ! .. Suitable value for dxi, which is smaller than the local grid spacing:
              IF (i > 1) THEN
                ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
                !  take the dxi for differentiation as a fraction of the "lower" diff (to lower neighbour):
                x_s = look_Z_meltsnow%q_i_lin(i)
                x_su = look_Z_meltsnow%q_i_lin(i-1)
                dxi = 0.01_dp * (x_s - x_su)
              ELSE
                ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
                !  to avoid negative q_i_lin below:
                x_s = look_Z_meltsnow%q_i_lin(i)
                x_so = look_Z_meltsnow%q_i_lin(i+1)
                dxi = 0.01_dp * (x_so - x_s)
                dxi = MIN(dxi, 0.05_dp * x_s)
              END IF

              DO kk=1, 3
                
                SELECT CASE (kk)
                CASE (1)
                  q_s = look_Z_meltsnow%q_i_lin(i) * n_s
                CASE (2)
                  q_s = (look_Z_meltsnow%q_i_lin(i) - dxi) * n_s
                CASE (3)
                  q_s = (look_Z_meltsnow%q_i_lin(i) + dxi) * n_s
                END SELECT

                mgd = mgd_2mom(snow,q_s,n_s)

                ! 0.5d0 = radienverh
                CALL zradar_wetsnow_mie_vec(&
                     mgd,look_Z_meltsnow%T_a(j),m_i(j),m_w(j),&
                     0.5d0,itype_Dref_fmelt,Tmeltbegin, meltdegTmin,&
                     lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                     Zfac,snow,rain,Dmin,Dmax,&
                     zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                     mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
                     mixingrulestring_core,matrixstring_core,inclusionstring_core,&
                     hoststring_shell,hostmatrixstring_shell,hostinclusionstring_shell,&
                     hoststring_core,hostmatrixstring_core,hostinclusionstring_core,&
                     look_Z_meltsnow%T_m(k),llookupgen_mode=.TRUE.)

                SELECT CASE (kk)
                CASE (1)              
                  ! Store parameters in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltsnow % zh   % val(i,j,k) = zh   / n_s
                  look_Z_meltsnow % ah   % val(i,j,k) = ah   / n_s
                  look_Z_meltsnow % zv   % val(i,j,k) = zv   / n_s
                  look_Z_meltsnow % rrhv % val(i,j,k) = rrhv / n_s
                  look_Z_meltsnow % irhv % val(i,j,k) = irhv / n_s
                  look_Z_meltsnow % kdp  % val(i,j,k) = kdp  / n_s
                  look_Z_meltsnow % adp  % val(i,j,k) = adp  / n_s
                  look_Z_meltsnow % zvh  % val(i,j,k) = zvh  / n_s
                CASE (2)
                  zhu   = zh
                  ahu   = ah
                  zvu   = zv
                  rrhvu = rrhv
                  irhvu = irhv
                  kdpu  = kdp
                  adpu  = adp
                  zvhu  = zvh
                CASE (3)
                  ! Store derivatives in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltsnow % zh   % dval(i,j,k) = (zh-zhu)     / dxi * 0.5_dp / n_s
                  look_Z_meltsnow % ah   % dval(i,j,k) = (ah-ahu)     / dxi * 0.5_dp / n_s
                  look_Z_meltsnow % zv   % dval(i,j,k) = (zv-zvu)     / dxi * 0.5_dp / n_s
                  look_Z_meltsnow % rrhv % dval(i,j,k) = (rrhv-rrhvu) / dxi * 0.5_dp / n_s
                  look_Z_meltsnow % irhv % dval(i,j,k) = (irhv-irhvu) / dxi * 0.5_dp / n_s
                  look_Z_meltsnow % kdp  % dval(i,j,k) = (kdp-kdpu)   / dxi * 0.5_dp / n_s
                  look_Z_meltsnow % adp  % dval(i,j,k) = (adp-adpu)   / dxi * 0.5_dp / n_s
                  look_Z_meltsnow % zvh  % dval(i,j,k) = (zvh-zvhu)   / dxi * 0.5_dp / n_s
                END SELECT
                
              END DO

           END DO
            
          END IF
          
        END DO
      END DO

      look_Z_meltsnow%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_meltsnow%magicnr,nTa,nxs,nTm,luse_tmatrix,ldo_nonsphere, &
      !     itype_Dref_fmelt,Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,lambda_radar,Dmin,Dmax, &
      !     mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      !     mixingrulestring_core,matrixstring_core,inclusionstring_core,&
      !     hoststring_shell,hostmatrixstring_shell,hostinclusionstring_shell,&
      !     hoststring_core,hostmatrixstring_core,hostinclusionstring_core,&
      !     ' '//TRIM(look_Z_meltsnow%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_meltsnow, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_meltsnow, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_meltsnow%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_meltsnow%nparams
      ! compute df(p(xi(xi_scaled)))/dxi_scaled = df/dp * dp/dxi * dxi/dxi_scaled = df/dp * dval * dxi/dxi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_meltsnow%dmem(:,:,:,k) =  scale_dval_lut(look_Z_meltsnow%mem(:,:,:,k), look_Z_meltsnow%dmem(:,:,:,k), &
           look_Z_meltsnow%qmem(:,k), &
           look_Z_meltsnow%flag_qi_scal(k), look_Z_meltsnow%f_scal_qi(k), &
           look_Z_meltsnow%flag_mem_scal(k), look_Z_meltsnow%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_meltsnow%mem(:,:,:,k) =  scale_val_lut(look_Z_meltsnow%mem(:,:,:,k), &
           look_Z_meltsnow%flag_mem_scal(k), look_Z_meltsnow%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_meltsnow_2mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_meltgraupel_2mom_lookupcreate(&
       look_Z_meltgraupel,tableprops,&
       igraupel_type,&
       itype_Dref_fmelt,&
       Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
       Zfac,Dmin,Dmax,graupel,&
       mixingrulestring,matrixstring,inclusionstring,&
       mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
       hoststring,hostmatrixstring,hostinclusionstring,&
       impipar_lookupgen, pe_start, pe_end,    &
       linterp_mode_dualpol,                   &
       savepath_read,savepath_write,unitnr,hashtext)
    
    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_meltgraupel
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP, pMPr
    INTEGER, INTENT(in)                   :: igraupel_type, itype_Dref_fmelt, unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, meltdegTmin, Tmax_min, Tmax_max, &
                                             lambda_radar, Zfac, Dmin, Dmax
    CLASS(particle), INTENT(in)           :: graupel
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring, &
                                             mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                                             hoststring, hostmatrixstring, hostinclusionstring, &
                                             savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1), m_w(tableprops%nTa+1)
    REAL(KIND=dp)       :: q_g, x_g, x_gu, x_go
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dxi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_meltgraupel_2mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! .. hypothetical fixed snow number concentration for computing
    !     the lookup table in terms of x_g with a subroutine that
    !     takes as input q_g and n_g:
    REAL(KIND=dp)    :: n_g

    ! compute n_g in such a way that all q_g >= q_crit_radar and at the
    ! same time n_g >= n_crit_radar, so that reflectivity is guaranteed to
    ! be computed for all points of the table:
    !n_g = MAX(q_crit_radar%graupel / graupel%x_min, n_crit_radar%graupel)
    n_g = MAX(1d-6 / tableprops%qilow, n_crit_radar)

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_meltgraupel, &
         &                    nqi   = tableprops%nqi, &   ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &   ! table T_a will have nTa+1 nodes
         &                    ntm   = tableprops%nTm, &   ! table T_m will have nTm+1 nodes
         &                    qilow = tableprops%qilow, & ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &  ! kg, linear scaling, not log10!
         &                    Talow = Tmeltbegin,               &  ! K
         &                    Taup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    Tmlow = MAX(Tmeltbegin,Tmax_min), &  ! K
         &                    Tmup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qg
         &                    f_eq  = 1.0_dp, &           ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zh,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%ah,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zv,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%rrhv, lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%irhv, lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%kdp,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%adp,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zvh,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zh,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%ah,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zv,   lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%rrhv, lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%irhv, lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%kdp,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%adp,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_meltgraupel%zvh,  lut=look_Z_meltgraupel, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_meltgraupel%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_2mom_meltgraupel_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_meltgraupel%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_meltgraupel%mem   = 0.0_dp
      look_Z_meltgraupel%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for melting graupel, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_w = m_complex_water_ray_vec(DBLE(lambda_radar), &
                                    MIN(MAX(look_Z_meltgraupel%T_a-T0C_fwo,mw_Tmin),mw_Tmax), &
                                    look_Z_meltgraupel%nTa)
      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_meltgraupel%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_meltgraupel%nTa)

      ! fill lookup table body
      DO k=1, look_Z_meltgraupel%nTm
        DO j=1, look_Z_meltgraupel%nTa

          IF (impipar_lookupgen == 2) THEN
            ! Parallelize lookup table nodes over nTa and nTm  (j and k) and feed the different q_i
            !  nodes to the same PE at constant T_a and T_m. Then, the Mie routine below will
            !  compute the scattering amplitudes for each particle size D_s only once and
            !  re-use it in the integration over size distribution instead of re-computing the same
            !  thing over and over.
            CALL sub2ind2D (j, k, look_Z_meltgraupel%nTa, work_pe)
            work_pe = round_robin(work_pe-1, pe_start, pe_end)
          ELSE
            work_pe = my_cart_id_fwo
          END IF

          IF (my_cart_id_fwo == work_pe) THEN

            DO i=1, look_Z_meltgraupel%nqi

              IF (luse_tmatrix .AND. ldebug_dbz) THEN
                ! Computation takes long, so print some progress information
                CALL progress_information (pos_=[i,j,k], &
                     shape_=[look_Z_meltgraupel%nqi,look_Z_meltgraupel%nTa,look_Z_meltgraupel%nTm], &
                     stride_=[look_Z_meltgraupel%nqi,1,1], cident='Melting graupel lookup creation')
              END IF

              ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
              ! .. Suitable value for dxi, which is smaller than the local grid spacing:
              IF (i > 1) THEN
                ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
                !  take the dxi for differentiation as a fraction of the "lower" diff (to lower neighbour):
                x_g = look_Z_meltgraupel%q_i_lin(i)
                x_gu = look_Z_meltgraupel%q_i_lin(i-1)
                dxi = 0.01_dp * (x_g - x_gu)
              ELSE
                ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
                !  to avoid negative q_i_lin below:
                x_g = look_Z_meltgraupel%q_i_lin(i)
                x_go = look_Z_meltgraupel%q_i_lin(i+1)
                dxi = 0.01_dp * (x_go - x_g)
                dxi = MIN(dxi, 0.05_dp * x_g)
              END IF

              DO kk=1, 3
                
                SELECT CASE (kk)
                CASE (1)
                  q_g = look_Z_meltgraupel%q_i_lin(i) * n_g
                CASE (2)
                  q_g = (look_Z_meltgraupel%q_i_lin(i) - dxi) * n_g
                CASE (3)
                  q_g = (look_Z_meltgraupel%q_i_lin(i) + dxi) * n_g
                END SELECT

                mgd = mgd_2mom(graupel,q_g,n_g)

                SELECT CASE (igraupel_type)
                CASE(1)
                  CALL zradar_wetgr_mie_vec(mgd,look_Z_meltgraupel%T_a(j),&
                       m_i(j),m_w(j),&
                       itype_Dref_fmelt,Tmeltbegin,meltdegTmin,&
                       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                       Zfac,graupel,rain,Dmin,Dmax,&
                       zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                       mixingrulestring,matrixstring,inclusionstring,&
                       hoststring,hostmatrixstring,hostinclusionstring,&
                       look_Z_meltgraupel%T_m(k),llookupgen_mode=.TRUE.)
                CASE(2)
                  CALL zradar_wetgr_twosph_mie_vec(mgd,look_Z_meltgraupel%T_a(j),&
                       m_i(j),m_w(j),&
                       itype_Dref_fmelt,Tmeltbegin,meltdegTmin,&
                       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                       Zfac,graupel,rain,Dmin,Dmax,&
                       zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                       mixingrulestring,matrixstring,inclusionstring,&
                       mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
                       look_Z_meltgraupel%T_m(k),llookupgen_mode=.TRUE.)
                CASE(3)
                  CALL zradar_wetgr_wsph_mie_vec(mgd,look_Z_meltgraupel%T_a(j),&
                       m_i(j),m_w(j),&
                       itype_Dref_fmelt,Tmeltbegin,meltdegTmin,&
                       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                       Zfac,graupel,rain,Dmin,Dmax,&
                       zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                       mixingrulestring,matrixstring,inclusionstring,&
                       look_Z_meltgraupel%T_m(k),llookupgen_mode=.TRUE.)
                END SELECT

                SELECT CASE (kk)
                CASE (1)
                  ! Store parameters in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltgraupel % zh   % val(i,j,k) = zh   / n_g
                  look_Z_meltgraupel % ah   % val(i,j,k) = ah   / n_g
                  look_Z_meltgraupel % zv   % val(i,j,k) = zv   / n_g
                  look_Z_meltgraupel % rrhv % val(i,j,k) = rrhv / n_g
                  look_Z_meltgraupel % irhv % val(i,j,k) = irhv / n_g
                  look_Z_meltgraupel % kdp  % val(i,j,k) = kdp  / n_g
                  look_Z_meltgraupel % adp  % val(i,j,k) = adp  / n_g
                  look_Z_meltgraupel % zvh  % val(i,j,k) = zvh  / n_g
                CASE (2)
                  zhu   = zh
                  ahu   = ah
                  zvu   = zv
                  rrhvu = rrhv
                  irhvu = irhv
                  kdpu  = kdp
                  adpu  = adp
                  zvhu  = zvh
                CASE (3)
                  ! Store derivatives in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_meltgraupel % zh   % dval(i,j,k) = (zh-zhu)     / dxi * 0.5_dp / n_g
                  look_Z_meltgraupel % ah   % dval(i,j,k) = (ah-ahu)     / dxi * 0.5_dp / n_g
                  look_Z_meltgraupel % zv   % dval(i,j,k) = (zv-zvu)     / dxi * 0.5_dp / n_g
                  look_Z_meltgraupel % rrhv % dval(i,j,k) = (rrhv-rrhvu) / dxi * 0.5_dp / n_g
                  look_Z_meltgraupel % irhv % dval(i,j,k) = (irhv-irhvu) / dxi * 0.5_dp / n_g
                  look_Z_meltgraupel % kdp  % dval(i,j,k) = (kdp-kdpu)   / dxi * 0.5_dp / n_g
                  look_Z_meltgraupel % adp  % dval(i,j,k) = (adp-adpu)   / dxi * 0.5_dp / n_g
                  look_Z_meltgraupel % zvh  % dval(i,j,k) = (zvh-zvhu)   / dxi * 0.5_dp / n_g
                END SELECT
                
              END DO

            END DO
            
          END IF

        END DO
      END DO

      look_Z_meltgraupel%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_meltgraupel%magicnr,igraupel_type,nTa,nxg,nTm,luse_tmatrix,ldo_nonsphere,&
      !     itype_Dref_fmelt,Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,lambda_radar,Dmin,Dmax,&
      !     mixingrulestring,matrixstring,inclusionstring,&
      !     mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      !     hoststring,hostmatrixstring,hostinclusionstring,' '//TRIM(look_Z_meltgraupel%chydroconfig)

      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_meltgraupel, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_meltgraupel, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_meltgraupel%is_initialized = .TRUE.
      END IF

    END IF


!!$    IF (my_cart_id_fwo == 0) THEN
!!$      DO k=1, look_Z_meltgraupel%nTm
!!$        DO j=1, look_Z_meltgraupel%nTa
!!$          DO i=1, look_Z_meltgraupel%nqi
!!$            WRITE (*,'(a,2es15.6)') 'ULImeltgraupel ', look_Z_meltgraupel%q_i_lin(i), &
!!$                 SQRT( (look_Z_meltgraupel%rrhv%val(i,j,k)**2+look_Z_meltgraupel%irhv%val(i,j,k)**2) / &
!!$                       (look_Z_meltgraupel%zh%val(i,j,k)*look_Z_meltgraupel%zv%val(i,j,k)) )
!!$          END DO
!!$        END DO
!!$      END DO
!!$    END IF


    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_meltgraupel%nparams
      ! compute df(p(xi(xi_scaled)))/dxi_scaled = df/dp * dp/dxi * dxi/dxi_scaled = df/dp * dval * dxi/dxi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_meltgraupel%dmem(:,:,:,k) =  scale_dval_lut(look_Z_meltgraupel%mem(:,:,:,k), look_Z_meltgraupel%dmem(:,:,:,k), &
           look_Z_meltgraupel%qmem(:,k), &
           look_Z_meltgraupel%flag_qi_scal(k), look_Z_meltgraupel%f_scal_qi(k), &
           look_Z_meltgraupel%flag_mem_scal(k), look_Z_meltgraupel%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_meltgraupel%mem(:,:,:,k) =  scale_val_lut(look_Z_meltgraupel%mem(:,:,:,k), &
           look_Z_meltgraupel%flag_mem_scal(k), look_Z_meltgraupel%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_meltgraupel_2mom_lookupcreate

  ! UB: Zfac has to be lambda^4 / (pi^5 * |K_w_0|^2)
  SUBROUTINE zradar_melthail_2mom_lookupcreate(&
       look_Z_melthail,tableprops,&
       itype_Dref_fmelt,&
       Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,&
       lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
       Zfac,Dmin,Dmax,hail,&
       mixingrulestring,matrixstring,inclusionstring,&
       impipar_lookupgen, pe_start, pe_end,    &
       linterp_mode_dualpol,                   &
       savepath_read,savepath_write,unitnr,hashtext)

    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_melthail
    TYPE(t_tabledef), INTENT(in)          :: tableprops
    LOGICAL, INTENT(in)                   :: luse_tmatrix, ldo_nonsphere
    TYPE(t_polMP), INTENT(in)             :: pMP, pMPr
    INTEGER, INTENT(in)                   :: itype_Dref_fmelt, unitnr
    REAL(KIND=dp), INTENT(in)             :: Tmeltbegin, meltdegTmin, Tmax_min, Tmax_max, &
                                             lambda_radar, Zfac, &
                                             Dmin, Dmax
    CLASS(particle), INTENT(in)           :: hail
    INTEGER, INTENT(in)                   :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                   :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)          :: mixingrulestring, matrixstring, inclusionstring, &
                                             savepath_read, savepath_write, hashtext

    TYPE(t_mgd_params)  :: mgd
    COMPLEX(kind=dp)    :: m_i(tableprops%nTa+1), m_w(tableprops%nTa+1)
    REAL(KIND=dp)       :: q_h, x_h, x_hu, x_ho
    REAL(KIND=dp)       :: zh,  ah,  zv,  rrhv,  irhv,  kdp,  adp,  zvh
    REAL(KIND=dp)       :: zhu, ahu, zvu, rrhvu, irhvu, kdpu, adpu, zvhu, dxi
    INTEGER             :: i, j, k, kk, ierr, work_pe
    CHARACTER(len=*), PARAMETER :: yzroutine = 'zradar_melthail_2mom_lookupcreate'
    CHARACTER(len=cmaxlen) :: errmsg
    CHARACTER(len=cmaxlen) :: savefile, tmppath_read, tmppath_write
    CHARACTER(len=10)   :: cmagicnr, cscattheo
    CHARACTER(len=4)    :: saveext
    !CHARACTER(len=3000) :: extrastr
    LOGICAL             :: fileexist, print_debug

    ! .. hypothetical fixed hail number concentration for computing
    !     the lookup table in terms of x_h with a subroutine that
    !     takes as input q_h and n_h:
    REAL(KIND=dp)    :: n_h

    ! compute n_h in such a way that all q_h >= q_crit_radar and at the
    ! same time n_h >= n_crit_radar, so that reflectivity is guaranteed to
    ! be computed for all points of the table:
    !n_h = MAX(q_crit_radar%hail / hail%x_min, n_crit_radar%hail)
    n_h = MAX(1d-6 / tableprops%qilow, n_crit_radar)

    ! NOTE: if you change something with the table vectors or table values below,
    !       please increase the version counter of the tables in radar_mie_iface_cosmo.f90!

    ! Allocate memory blocks, set up internal table pointers and define table nodes:
    CALL init_dbzlookuptable (lut   = look_Z_melthail, &
         &                    nqi   = tableprops%nqi, &   ! table q_i will have nqr+1 nodes
         &                    nta   = tableprops%nTa, &   ! table T_a will have nTa+1 nodes
         &                    ntm   = tableprops%nTm, &   ! table T_m will have nTm+1 nodes
         &                    qilow = tableprops%qilow, & ! kg, linear scaling, not log10!
         &                    qiup  = tableprops%qiup, &  ! kg, linear scaling, not log10!
         &                    Talow = Tmeltbegin,               &  ! K
         &                    Taup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    Tmlow = MAX(Tmeltbegin,Tmax_min), &  ! K
         &                    Tmup  = MAX(Tmeltbegin,Tmax_max), &  ! K
         &                    flag_qi_scal_eq = i_scal_log, &  ! logarithmic equidistant table for qi
         &                    f_eq  = 1.0_dp, &           ! for log scaling, this is a dummy
         &                    ierr = ierr)

    IF (ierr /= 0) THEN
      ! an error message has already been written to stdout by init_dbzlookuptable()
      STOP
    END IF
    
    ! Define interpolation scaling for each parameter w.r.t. qr (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    IF (linterp_mode_dualpol) THEN
      
      CALL init_scaling_tabparams (tpar=look_Z_melthail%zh,   lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%ah,   lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%zv,   lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%rrhv, lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%irhv, lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%kdp,  lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%adp,  lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%zvh,  lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)

    ELSE
    
      CALL init_scaling_tabparams (tpar=look_Z_melthail%zh,   lut=look_Z_melthail, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%ah,   lut=look_Z_melthail, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_log, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%zv,   lut=look_Z_melthail, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%rrhv, lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%irhv, lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%kdp,  lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%adp,  lut=look_Z_melthail, &
           flag_scal_qi=i_scal_lin, f_scal_qi=1.0_dp, flag_scal_val=i_scal_lin, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_cubic)
      CALL init_scaling_tabparams (tpar=look_Z_melthail%zvh,  lut=look_Z_melthail, &
           flag_scal_qi=i_scal_log, f_scal_qi=1.0_dp, flag_scal_val=i_scal_dbz, f_scal_val=1.0_dp, &
           flag_interp_val_qi=i_interp_lin)

    END IF

    ! Put together the filename for the binary table storage file,
    !   and, if it already exists, read table from file, otherwise
    !   create table and write it to file:
    cmagicnr(:) = ' '
    WRITE(cmagicnr, '(i0.10)') look_Z_melthail%magicnr
    cscattheo(:) = ' '
    IF (luse_tmatrix) THEN
      cscattheo = 'tmatrix'
    ELSE
      cscattheo = 'mie'
    END IF
    savefile(:) = ' '
    savefile = 'radar_'//TRIM(cscattheo)//'tab_2mom_melthail_'//TRIM(ADJUSTL(cmagicnr))//'_'//&
               TRIM(look_Z_melthail%cversion_lt)
    tmppath_read(:) = ' '
    IF (LEN_TRIM(savepath_read) > 0) THEN
      tmppath_read = TRIM(savepath_read)//'/'
    END IF
    tmppath_write(:) = ' '
    IF (LEN_TRIM(savepath_write) > 0) THEN
      tmppath_write = TRIM(savepath_write)//'/'
    END IF
    saveext(:) = ' '
    IF (llut_format_netcdf) THEN
      saveext = '.nc'
    ELSE
      saveext = '.bin'
    END IF
    INQUIRE(file=TRIM(tmppath_read)//TRIM(savefile)//TRIM(saveext), exist=fileexist)

    SELECT CASE (impipar_lookupgen)
    CASE (1)
      ! In this case, either the file exists and this routine is called by all PEs
      !  to read the table from the file. Or the file does not exist and this
      !  routine is called by only one PE to create the table.
      print_debug = ldebug_dbz
    CASE (2)
      ! In this case, this routine is called by all PEs to create the table in parallel
      !  with sparse computations on some table nodes of the entire table array.
      ! A collective mpi_reduce('SUM')  (hidden in global_values_radar('SUM') is used
      ! to gather all field elements on the output PE.
      ! Therefore we need to initialize the table values with 0.0:
      look_Z_melthail%mem   = 0.0_dp
      look_Z_melthail%dmem  = 0.0_dp
      !  ... and for safety a parallel synchronization of the check on file existence is done here:
      ierr = 0
      errmsg(:) = ' '
      CALL global_values_radar (fileexist, 'OR', icomm_cart_fwo, -1, errmsg, ierr)
      print_debug = ldebug_dbz .AND. pe_start <= my_cart_id_fwo .AND. my_cart_id_fwo <= pe_end
    END SELECT
    
    IF (.NOT. fileexist) THEN

      IF (print_debug) THEN
        WRITE (*,'(a,i0)') 'DEBUG: computing lookup table for melting hail, '//&
                    TRIM(ADJUSTL(cmagicnr))//' on proc ', my_cart_id_fwo
      END IF

      m_w = m_complex_water_ray_vec(DBLE(lambda_radar), &
                                    MIN(MAX(look_Z_melthail%T_a-T0C_fwo,mw_Tmin),mw_Tmax), &
                                    look_Z_melthail%nTa)
      m_i = m_complex_ice_maetzler_vec(DBLE(lambda_radar), &
                                       MIN(MAX(look_Z_melthail%T_a-T0C_fwo,mi_Tmin),mi_Tmax), &
                                       look_Z_melthail%nTa)

      ! fill lookup table body
      DO k=1, look_Z_melthail%nTm
        DO j=1, look_Z_melthail%nTa

          IF (impipar_lookupgen == 2) THEN
            CALL sub2ind2D (j, k, look_Z_melthail%nTa, work_pe)
            work_pe = round_robin(work_pe-1, pe_start, pe_end)
          ELSE
            work_pe = my_cart_id_fwo
          END IF

          IF (my_cart_id_fwo == work_pe) THEN
              
            DO i=1, look_Z_melthail%nqi
              
              IF (luse_tmatrix .AND. ldebug_dbz) THEN
                ! Computation takes long, so print some progress information
                CALL progress_information (pos_=[i,j,k], &
                     shape_=[look_Z_melthail%nqi,look_Z_melthail%nTa,look_Z_melthail%nTm], &
                     stride_=[look_Z_melthail%nqi,1,1], cident='Melting hail lookup creation')
              END IF
            
              ! Prepare computation of the derivatives w.r.t. q_i_lin (centered diff):
              ! .. Suitable value for dxi, which is smaller than the local grid spacing:
              IF (i > 1) THEN
                ! grid spacing of q_i_lin normally increases with increasing q_i_lin, so
                !  take the dxi for differentiation as a fraction of the "lower" diff (to lower neighbour):
                x_h = look_Z_melthail%q_i_lin(i)
                x_hu = look_Z_melthail%q_i_lin(i-1)
                dxi = 0.01_dp * (x_h - x_hu)
              ELSE
                ! lowest q_i_lin value. Estimate from "higher" diff, but limit to 0.1 lowest value
                !  to avoid negative q_i_lin below:
                x_h = look_Z_melthail%q_i_lin(i)
                x_ho = look_Z_melthail%q_i_lin(i+1)
                dxi = 0.01_dp * (x_ho - x_h)
                dxi = MIN(dxi, 0.05_dp * x_h)
              END IF

              DO kk=1, 3
                
                SELECT CASE (kk)
                CASE (1)
                  q_h = look_Z_melthail%q_i_lin(i) * n_h
                CASE (2)
                  q_h = (look_Z_melthail%q_i_lin(i) - dxi) * n_h
                CASE (3)
                  q_h = (look_Z_melthail%q_i_lin(i) + dxi) * n_h
                END SELECT

                mgd = mgd_2mom(hail,q_h,n_h)

                ! Spongy Hail --- Two-sphere-model with a spherical ice core and a film
                !                 of ice-water-mixture with 70% water
                CALL zradar_wethail_mie_vec(mgd,look_Z_melthail%T_a(j),&
                     m_i(j),m_w(j),&
                     itype_Dref_fmelt,Tmeltbegin,meltdegTmin,&
                     lambda_radar,luse_tmatrix,ldo_nonsphere,pMP,pMPr,&
                     Zfac,hail,rain,Dmin,Dmax,&
                     zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                     mixingrulestring,matrixstring,inclusionstring,&
                     look_Z_melthail%T_m(k),llookupgen_mode=.TRUE.)

                SELECT CASE (kk)
                CASE (1)
                  ! Store parameters in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_melthail % zh   % val(i,j,k) = zh   / n_h 
                  look_Z_melthail % ah   % val(i,j,k) = ah   / n_h 
                  look_Z_melthail % zv   % val(i,j,k) = zv   / n_h 
                  look_Z_melthail % rrhv % val(i,j,k) = rrhv / n_h 
                  look_Z_melthail % irhv % val(i,j,k) = irhv / n_h 
                  look_Z_melthail % kdp  % val(i,j,k) = kdp  / n_h 
                  look_Z_melthail % adp  % val(i,j,k) = adp  / n_h 
                  look_Z_melthail % zvh  % val(i,j,k) = zvh  / n_h 
                CASE (2)
                  zhu   = zh
                  ahu   = ah
                  zvu   = zv
                  rrhvu = rrhv
                  irhvu = irhv
                  kdpu  = kdp
                  adpu  = adp
                  zvhu  = zvh
                CASE (3)
                  ! Store derivatives in original linear scaling in the tables. Transform to
                  ! desired scaling after table writing/reading:
                  look_Z_melthail % zh   % dval(i,j,k) = (zh-zhu)     / dxi * 0.5_dp / n_h
                  look_Z_melthail % ah   % dval(i,j,k) = (ah-ahu)     / dxi * 0.5_dp / n_h
                  look_Z_melthail % zv   % dval(i,j,k) = (zv-zvu)     / dxi * 0.5_dp / n_h
                  look_Z_melthail % rrhv % dval(i,j,k) = (rrhv-rrhvu) / dxi * 0.5_dp / n_h
                  look_Z_melthail % irhv % dval(i,j,k) = (irhv-irhvu) / dxi * 0.5_dp / n_h
                  look_Z_melthail % kdp  % dval(i,j,k) = (kdp-kdpu)   / dxi * 0.5_dp / n_h
                  look_Z_melthail % adp  % dval(i,j,k) = (adp-adpu)   / dxi * 0.5_dp / n_h
                  look_Z_melthail % zvh  % dval(i,j,k) = (zvh-zvhu)   / dxi * 0.5_dp / n_h
                END SELECT
                
              END DO

            END DO
            
          END IF

        END DO
      END DO

      look_Z_melthail%is_initialized = .TRUE.

      ! JM220224:
      ! changed write_lookup_to_file's extrastr to contain the full hash-forming string
      ! however, before removing old setting entirely, check that the contents are
      ! represented in the hash. or in any other sufficiently explicit way in the lookup
      ! file. and decide whether it should be.
      !extrastr(:) = ' '
      !WRITE(extrastr,*) &
      !     look_Z_melthail%magicnr,nTa,nxh,nTm,luse_tmatrix,ldo_nonsphere,&
      !     itype_Dref_fmelt,Tmeltbegin,meltdegTmin,Tmax_min,Tmax_max,lambda_radar,Dmin,Dmax, &
      !     mixingrulestring,matrixstring,inclusionstring,' '//TRIM(look_Z_melthail%chydroconfig)
      
      CALL write_lookup_to_file(unitnr, impipar_lookupgen, pe_start, &
                                TRIM(tmppath_write)//TRIM(savefile), &
                                look_Z_melthail, yzroutine, TRIM(hashtext))

    ELSE

      CALL read_lookup_from_file(unitnr, TRIM(tmppath_read)//TRIM(savefile), &
                                 look_Z_melthail, yzroutine, pe_start, ierr)
      IF (ierr == 0) THEN
        look_Z_melthail%is_initialized = .TRUE.
      END IF

    END IF

    ! Convert parameters from linear space to the desired scaling for interpolation:
    DO k=1, look_Z_melthail%nparams
      ! compute df(p(xi(xi_scaled)))/dxi_scaled = df/dp * dp/dxi * dxi/dxi_scaled = df/dp * dval * dxi/dxi_scaled,
      !  where df/dp may depend on val itself.
      look_Z_melthail%dmem(:,:,:,k) =  scale_dval_lut(look_Z_melthail%mem(:,:,:,k), look_Z_melthail%dmem(:,:,:,k), &
           look_Z_melthail%qmem(:,k), &
           look_Z_melthail%flag_qi_scal(k), look_Z_melthail%f_scal_qi(k), &
           look_Z_melthail%flag_mem_scal(k), look_Z_melthail%f_scal_mem(k))
      ! compute f(p) = scaled value of mem:
      look_Z_melthail%mem(:,:,:,k) =  scale_val_lut(look_Z_melthail%mem(:,:,:,k), &
           look_Z_melthail%flag_mem_scal(k), look_Z_melthail%f_scal_mem(k))
    END DO

  END SUBROUTINE zradar_melthail_2mom_lookupcreate

#endif

  !===================================================================================================
  ! 
  ! Subroutine for printing progress information for dualpol computations, if the actual
  ! coordinates are given in pos_, the shape of the table is given in shape_ and the interval
  ! for printing the info string in each dimension is given in stride_.
  !
  ! INPUT:
  !   pos_     INTEGER(N)   : N-dimensional coordinate vector
  !   shape_   INTEGER(N)   : N-dimensional vector containing the size of the lookup tables in each direction
  !   stride_  INTEGER(N)   : N-dimensional vector containing the interval for printing in each direction
  !   cident   CHARACTER(*) : String to identify the caller, used as prefix of the info string
  !
  !   where N is an arbitrary number of dimensions (2, 3, 4 or even more dimensions possible).
  !
  ! OUTPUT:
  !   A character string is printed on standard output.
  ! 
  !===================================================================================================

  SUBROUTINE progress_information (pos_, shape_, stride_, cident)

    INTEGER,          INTENT(in) :: pos_(:)     ! actual position in array
    INTEGER,          INTENT(in) :: shape_(:)   ! shape of array
    INTEGER,          INTENT(in) :: stride_(:)  ! stride along each dimensions after which progress message should be printed
    CHARACTER(len=*), INTENT(in) :: cident      ! String to identify the caller, will be part of the progress message

    INTEGER            :: ndims, i
    CHARACTER(len=3)   :: cndims
    CHARACTER(len=100) :: posstr, shapestr

    IF ( ALL(MODULO(pos_(:), MIN(stride_(:),shape_(:))) == 0) .OR. ALL(pos_(:) == 1) .OR. ALL(pos_(:) == shape_(:))) THEN

      ndims = SIZE(shape_)
      WRITE (cndims, '(i3.3)') ndims

      ! string for "(i,j,k,...)"
      posstr(:) = ' '
      WRITE (posstr, '("(",'//TRIM(cndims)//'a,")")') &
           (ACHAR(IACHAR('i')+i)//',', i=0,ndims-2), ACHAR(IACHAR('i')+ndims-1)
    
      ! string for (ni,nj,nk,...)
      shapestr(:) = ' '
      WRITE (shapestr, '("(",'//TRIM(cndims)//'a,")")') &
           ('n'//ACHAR(IACHAR('i')+i)//',', i=0,ndims-2), 'n'//ACHAR(IACHAR('i')+ndims-1)
    
      WRITE (*,'(a,": Tmatrix-computations on proc ",i4," for table node ",a," ",'//TRIM(cndims)//'i5,' // &
           '"   from ",a," ",'//TRIM(cndims)//'i5," ...")') &
           '[PROGRESS] '//TRIM(cident), my_cart_id_fwo, TRIM(posstr), pos_(:), TRIM(shapestr), shape_(:)
      
    END IF

  END SUBROUTINE progress_information
    
  !===================================================================================================
  !===================================================================================================

  !*******************************************************************************
  !
  ! Routines for writing and reading of lookup tables
  !
  !*******************************************************************************

  SUBROUTINE write_lookup_to_file(&
       unitnr, impipar_lookupgen, pe_start, savefile, look_Z_type, yzroutine, extrastr)

    IMPLICIT NONE

    INTEGER, INTENT(in)                   :: unitnr, impipar_lookupgen, pe_start
    CHARACTER(len=*), INTENT(in)          :: savefile, yzroutine, extrastr
    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_type

    INTEGER :: ierr, k, status, cmode, output_pe

    INTEGER :: fid, nqi_dimid, nqilin_dimid, nTa_dimid, nTm_dimid, dimid_qi
    INTEGER :: varid_qi, varid_qi_lin, varid_Ta, varid_Tm
    INTEGER :: varid_zh, varid_ah
    INTEGER :: varid_dzh, varid_dah
    INTEGER :: varid_zv, varid_rrhv, varid_irhv, varid_kdp, varid_adp, varid_zvh
    INTEGER :: varid_dzv, varid_drrhv, varid_dirhv, varid_dkdp, varid_dadp, varid_dzvh
    LOGICAL :: fileexist

    CHARACTER(len=cmaxlen) :: errmsg

    CHARACTER(len=200) :: formatnqi, formatnTa, formatnTm, format2d, format3d, formattmp

    ! The extra string will be written as is to the debug ASCII table file below.
    ! Normally it should contain the dbz configuration parameters which were applied
    ! to build the table.

    IF (impipar_lookupgen == 2 .AND. num_compute_fwo > 1) THEN

      ! In this case, the table values are sparsely scattered over all compute PEs in icomm_cart_fwo
      !  and this subroutine is called by all PEs in icomm_cart_fwo.
      ! Gather these sparsely scattered data on all PEs (so that all PEs have a valid table)
      !  by a collective MPI_REDUCE('SUM') (hidden in global_values_radar('SUM'):

      IF (ASSOCIATED(look_Z_type%mem)) THEN
        DO k = 1, look_Z_type%nparams
          CALL global_values_radar(look_Z_type%mem(:,:,:,k), SHAPE(look_Z_type%mem(:,:,:,k)), 'SUM', &
               icomm_cart_fwo, -1, errmsg, ierr)
        END DO
      END IF
      IF (ASSOCIATED(look_Z_type%dmem)) THEN
        DO k = 1, look_Z_type%nparams
          CALL global_values_radar(look_Z_type%dmem(:,:,:,k), SHAPE(look_Z_type%dmem(:,:,:,k)), 'SUM', &
               icomm_cart_fwo, -1, errmsg, ierr)
        END DO
      END IF

      ! The tables are output to file by just one output PE:
      output_pe = pe_start
    ELSE
      ! In this case, the table values are already completely on my_cart_id_fwo and no other
      ! compute PE calls this subroutine. So we can just output the table without any
      ! further parallel communication:
      output_pe = my_cart_id_fwo
    END IF

    ! Write table to file, if the file does not yet exist:
    IF (my_cart_id_fwo == output_pe) THEN

      ! Inquire, if the table file still does not exist. Maybe another Ensemble member
      ! has produced the table inbetween in case of ...
      IF (llut_format_netcdf) THEN
        INQUIRE(file=TRIM(savefile)//'.nc', exist=fileexist)
      ELSE
        INQUIRE(file=TRIM(savefile)//'.bin', exist=fileexist)
      END IF

      IF (.NOT. fileexist) THEN

        IF (.NOT. llut_format_netcdf) THEN

          WRITE (*,'(a,i0)') 'INFO '//TRIM(yzroutine)//&
               ': Writing binary Mie- or T-matrix lookup table to file '//TRIM(savefile)//'.bin on proc ', my_cart_id_fwo

          OPEN(unitnr, file=TRIM(savefile)//'.bin', status='replace', form='unformatted', iostat=ierr)
          IF (ierr /= 0) THEN

            WRITE (*,*) 'ERROR '//TRIM(yzroutine)//': binary lookup-table file ' &
                 // TRIM(savefile)// '.bin could not be created!'
            STOP
           
          ELSE

            WRITE(unitnr) look_Z_type%cversion_lt
            WRITE(unitnr) look_Z_type%magicnr
            WRITE(unitnr) look_Z_type%nqi
            WRITE(unitnr) look_Z_type%nTa
            WRITE(unitnr) look_Z_type%nTm
            WRITE(unitnr) look_Z_type%nparams
            WRITE(unitnr) look_Z_type%i_zh
            WRITE(unitnr) look_Z_type%i_ah
            WRITE(unitnr) look_Z_type%i_zv
            WRITE(unitnr) look_Z_type%i_rrhv
            WRITE(unitnr) look_Z_type%i_irhv
            WRITE(unitnr) look_Z_type%i_kdp
            WRITE(unitnr) look_Z_type%i_adp
            WRITE(unitnr) look_Z_type%i_zvh
            WRITE(unitnr) look_Z_type%flag_qi_scal_eq
            WRITE(unitnr) look_Z_type%f_eq
            WRITE(unitnr) look_Z_type%if_eq

            WRITE(unitnr) look_Z_type%qi0
            WRITE(unitnr) look_Z_type%Ta0
            WRITE(unitnr) look_Z_type%Tm0

            WRITE(unitnr) look_Z_type%dqi
            WRITE(unitnr) look_Z_type%dTa
            WRITE(unitnr) look_Z_type%dTm

            WRITE(unitnr) look_Z_type%idqi
            WRITE(unitnr) look_Z_type%idTa
            WRITE(unitnr) look_Z_type%idTm

            WRITE(unitnr) look_Z_type%q_i
            WRITE(unitnr) look_Z_type%q_i_lin
            WRITE(unitnr) look_Z_type%T_a
            WRITE(unitnr) look_Z_type%T_m

            WRITE(unitnr) look_Z_type%mem
            WRITE(unitnr) look_Z_type%dmem

            CLOSE(unitnr)

            WRITE (*,'(a,i0)') 'INFO '//TRIM(yzroutine)//&
                 ': Finished writing binary lookup table to file '//TRIM(savefile)//'.bin on proc ', my_cart_id_fwo

          END IF

        END IF

        IF (ldebug_dbz) THEN

          !---------------------------------------------------------------------
          ! Debug output of the entire lookup table as an ASCII file
          !---------------------------------------------------------------------

          formatnqi(:) = ' '
          WRITE (formatnqi, '("(a,T10,",i3.3,"(es14.5),/)")') look_Z_type%nqi
          formatnTa(:) = ' '
          WRITE (formatnTa, '("(a,T10,",i3.3,"(es14.5),/)")') look_Z_type%nTa
          formatnTm(:) = ' '
          WRITE (formatnTm, '("(a,T10,",i3.3,"(es14.5),/)")') look_Z_type%nTm
          format2d(:) = ' '
          WRITE (format2d, '("(a,T10,",i3.3,"(T10,",i3.3,"(es14.5),/))")') look_Z_type%nTa, look_Z_type%nqi
          format3d(:) = ' '
          WRITE (format3d, '("(a,T10,",i3.3,"(",i3.3,"(T10,",i3.3,"(es14.5),/),/))")') &
               look_Z_type%nTm, look_Z_type%nTa, look_Z_type%nqi

          OPEN(unitnr, file=TRIM(savefile)//'.dat', status='replace', form='formatted', iostat=ierr)
          IF (ierr /= 0) THEN
            
            WRITE (*,*) 'ERROR '//TRIM(yzroutine)//&
                 ': ASCII lookup-table file for checking purposes '//&
                 TRIM(savefile)// '.dat could not be created!'
            
          ELSE

            WRITE(unitnr,'(a,T30,a)')  'Mie-lt-Version', look_Z_type%cversion_lt
            WRITE(unitnr,'(a,T30,i0)') 'Magic-Number', look_Z_type%magicnr
            WRITE(unitnr,'(a,T30,a)')  'Config-string', TRIM(extrastr)
            WRITE(unitnr,'(a,T30,i0)') 'itype_refl', look_Z_type%itype_refl
            WRITE(unitnr,'(a,T30,L1)') 'luse_tmatrix', look_Z_type%luse_tmatrix
            WRITE(unitnr,'(a,T30,L1)') 'ldo_nonsphere', look_Z_type%ldo_nonsphere
            WRITE(unitnr,'(a,T30,i0)') 'itype_Dref_fmelt', look_Z_type%itype_Dref_fmelt

            WRITE(unitnr,'(a,T30,i0)')   'flag_qi_scal_eq', look_Z_type%flag_qi_scal_eq
            WRITE(unitnr,'(a,T30,f0.5)') 'f_eq', look_Z_type%f_eq
            WRITE(unitnr,'(a,T30,f0.5)') 'if_eq', look_Z_type%if_eq

            WRITE(unitnr,'(a,T10,i0)') 'nqi', look_Z_type%nqi
            WRITE(unitnr,'(a,T10,i0)') 'nTa', look_Z_type%nTa
            WRITE(unitnr,'(a,T10,i0)') 'nTm', look_Z_type%nTm
            WRITE(unitnr,'(a,T10,i0)') 'nparams', look_Z_type%nparams

            WRITE(unitnr,'(a,T10,es14.5)') 'qi0', look_Z_type%qi0
            WRITE(unitnr,'(a,T10,es14.5)') 'Ta0', look_Z_type%Ta0
            WRITE(unitnr,'(a,T10,es14.5)') 'Tm0', look_Z_type%Tm0
          
            WRITE(unitnr,'(a,T10,es14.5)') 'dqi', look_Z_type%dqi
            WRITE(unitnr,'(a,T10,es14.5)') 'dTa', look_Z_type%dTa
            WRITE(unitnr,'(a,T10,es14.5)') 'dTm', look_Z_type%dTm

            WRITE(unitnr,'(a,T10,es14.5)') 'idqi', look_Z_type%idqi
            WRITE(unitnr,'(a,T10,es14.5)') 'idTa', look_Z_type%idTa
            WRITE(unitnr,'(a,T10,es14.5,/)') 'idTm', look_Z_type%idTm

            WRITE(unitnr,TRIM(formatnqi)) 'q_i', REAL(look_Z_type%q_i)
            WRITE(unitnr,TRIM(formatnqi)) 'q_i_lin', REAL(look_Z_type%q_i_lin)
            WRITE(unitnr,TRIM(formatnTa)) 'T_a', REAL(look_Z_type%T_a)
            WRITE(unitnr,TRIM(formatnTm)) 'T_m', REAL(look_Z_type%T_m)

            formattmp(:) = ' '
            IF (look_Z_type%ntm <= 1) THEN
              formattmp(:) = TRIM(format2d)
            ELSE
              formattmp(:) = TRIM(format3d)
            END IF
            IF (ASSOCIATED(look_Z_type%zh%val))   WRITE(unitnr,TRIM(formattmp)) 'zh',   REAL(look_Z_type%zh %val   )
            IF (ASSOCIATED(look_Z_type%zh%val))   WRITE(unitnr,TRIM(formattmp)) 'dzh',   REAL(look_Z_type%zh%dval  )
            IF (ASSOCIATED(look_Z_type%ah%val))   WRITE(unitnr,TRIM(formattmp)) 'ah',   REAL(look_Z_type%ah %val   )
            IF (ASSOCIATED(look_Z_type%ah%val))   WRITE(unitnr,TRIM(formattmp)) 'dah',   REAL(look_Z_type%ah%dval  )
            IF (ASSOCIATED(look_Z_type%zv%val))   WRITE(unitnr,TRIM(formattmp)) 'zv',   REAL(look_Z_type%zv %val   )
            IF (ASSOCIATED(look_Z_type%zv%val))   WRITE(unitnr,TRIM(formattmp)) 'dzv',   REAL(look_Z_type%zv%dval  )
            IF (ASSOCIATED(look_Z_type%rrhv%val)) WRITE(unitnr,TRIM(formattmp)) 'rrhv', REAL(look_Z_type%rrhv%val  )
            IF (ASSOCIATED(look_Z_type%rrhv%val)) WRITE(unitnr,TRIM(formattmp)) 'drrhv', REAL(look_Z_type%rrhv%dval)
            IF (ASSOCIATED(look_Z_type%irhv%val)) WRITE(unitnr,TRIM(formattmp)) 'irhv', REAL(look_Z_type%irhv%val  )
            IF (ASSOCIATED(look_Z_type%irhv%val)) WRITE(unitnr,TRIM(formattmp)) 'dirhv', REAL(look_Z_type%irhv%dval)
            IF (ASSOCIATED(look_Z_type%kdp%val))  WRITE(unitnr,TRIM(formattmp)) 'kdp',  REAL(look_Z_type%kdp %val  )
            IF (ASSOCIATED(look_Z_type%kdp%val))  WRITE(unitnr,TRIM(formattmp)) 'dkdp',  REAL(look_Z_type%kdp%dval )
            IF (ASSOCIATED(look_Z_type%adp%val))  WRITE(unitnr,TRIM(formattmp)) 'adp',  REAL(look_Z_type%adp %val  )
            IF (ASSOCIATED(look_Z_type%adp%val))  WRITE(unitnr,TRIM(formattmp)) 'dadp',  REAL(look_Z_type%adp%dval )
            IF (ASSOCIATED(look_Z_type%zvh%val))  WRITE(unitnr,TRIM(formattmp)) 'zvh',  REAL(look_Z_type%zvh %val  )
            IF (ASSOCIATED(look_Z_type%zvh%val))  WRITE(unitnr,TRIM(formattmp)) 'dzvh',  REAL(look_Z_type%zvh%dval )

            CLOSE(unitnr)
            
          END IF

        END IF


#ifdef NETCDF
        IF (ldebug_dbz .OR. llut_format_netcdf) THEN

          IF (llut_format_netcdf) THEN
            WRITE (*,'(a,i0)') 'INFO '//TRIM(yzroutine)//&
                 ': Writing netcdf Mie- or T-matrix lookup table to file '//TRIM(savefile)//'.nc on proc ', my_cart_id_fwo
          END IF

          !---------------------------------------------------------------------
          ! Output of the entire lookup table as a NetCDF file
          !---------------------------------------------------------------------

!!$ FIXME: after having checked equality to former code versions, adjust variable names and add units

          ! 1) Open and create Netcdf File
          ! ------------------------------

          cmode = NF90_CLOBBER
!!$ UB: We go back to netcdf-3, because it is faster on input. For this, the special netcdf-4 flags are commented out:
!!$        cmode = IOR(cmode, NF90_NETCDF4)
!!$        cmode = IOR(cmode, NF90_CLASSIC_MODEL)  ! no fancy netcdf4 features (new data types, compounds, multiple unlimited dims, etc.), therefore smaller files

          fid = -1
          status = check_nc( nf90_create (TRIM(savefile)//'.nc', cmode, fid), 'nf90_create '//TRIM(savefile)//'.nc' )

          IF (status /= nf90_noerr) THEN
            
            WRITE(*,'(a)') 'ERROR creating NetCDF file '//TRIM(savefile)//'.nc '// &
                 TRIM(nf90_strerror(status))
            fid = -1
            IF (llut_format_netcdf) STOP
            
          ELSE

            status = check_nc( nf90_def_dim(fid, 'q_i',  look_Z_type%nqi, nqi_dimid), 'nf90_def_dim q_i')
            status = check_nc( nf90_def_dim(fid, 'T_a',  look_Z_type%nTa, nTa_dimid), 'nf90_def_dim T_a')
            status = check_nc( nf90_def_dim(fid, 'T_m',  look_Z_type%nTm, nTm_dimid), 'nf90_def_dim T_m')

            status = check_nc( nf90_def_var(fid, 'q_i',    NF90_DOUBLE, (/nqi_dimid/), varid_qi), 'nf90_def_var q_i')
            status = check_nc( nf90_put_att(fid,varid_qi , "qi0", look_Z_type%qi0),  'nf90_put_att qi0 for q_i')
            status = check_nc( nf90_put_att(fid,varid_qi , "dqi", look_Z_type%dqi),  'nf90_put_att dqi for q_i')

            status = check_nc( nf90_put_att(fid,varid_qi , "idqi", look_Z_type%idqi),  'nf90_put_att idqi for q_i')
            status = check_nc( nf90_def_dim(fid, 'q_i_lin',  look_Z_type%nqi, nqilin_dimid), 'nf90_def_dim q_i')
            status = check_nc( nf90_def_var(fid, 'q_i_lin',  NF90_DOUBLE, (/nqilin_dimid/), varid_qi_lin), 'nf90_def_var q_i_lin')

            status = check_nc( nf90_def_var(fid, 'T_a',    NF90_DOUBLE, (/nTa_dimid/), varid_Ta), 'nf90_def_var T_a')
            status = check_nc( nf90_put_att(fid,varid_Ta , "Ta0", look_Z_type%Ta0),  'nf90_put_att Ta0 for T_a')
            status = check_nc( nf90_put_att(fid,varid_Ta , "dTa", look_Z_type%dTa),  'nf90_put_att dTa for T_a')
            status = check_nc( nf90_put_att(fid,varid_Ta , "idTa", look_Z_type%idTa),  'nf90_put_att idTa for T_a')

            status = check_nc( nf90_def_var(fid, 'T_m',    NF90_DOUBLE, (/nTm_dimid/), varid_Tm), 'nf90_def_var T_m')
            status = check_nc( nf90_put_att(fid,varid_Tm , "Tm0", look_Z_type%Tm0),  'nf90_put_att Tm0 for T_m')
            status = check_nc( nf90_put_att(fid,varid_Tm , "dTm", look_Z_type%dTm),  'nf90_put_att dTm for T_m')
            status = check_nc( nf90_put_att(fid,varid_Tm , "idTm", look_Z_type%idTm),  'nf90_put_att idTm for T_m')

            ! Which qi-axis scaling should be used for the variables? q_i (log) or q_i_lin?
            !        dimid_qi = nqi_dimid
            dimid_qi = nqilin_dimid

            status = check_nc( nf90_def_var(fid, 'zh',   NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_zh), &
                 'nf90_def_var zh')
            status = check_nc( nf90_def_var(fid, 'ah',   NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_ah), &
                 'nf90_def_var ah')
            status = check_nc( nf90_def_var(fid, 'zv',   NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_zv), &
                 'nf90_def_var zv')
            status = check_nc( nf90_def_var(fid, 'rrhv', NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_rrhv), &
                 'nf90_def_var rrhv')
            status = check_nc( nf90_def_var(fid, 'irhv', NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_irhv), &
                 'nf90_def_var irhvm')
            status = check_nc( nf90_def_var(fid, 'kdp',  NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_kdp), &
                 'nf90_def_var kdp')
            status = check_nc( nf90_def_var(fid, 'adp',  NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_adp), &
                 'nf90_def_var adp')
            status = check_nc( nf90_def_var(fid, 'zvh',  NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_zvh), &
                 'nf90_def_var zvh')

            status = check_nc( nf90_def_var(fid, 'dzh',   NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_dzh), &
                 'nf90_def_var dzh')
            status = check_nc( nf90_def_var(fid, 'dah',   NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_dah), &
                 'nf90_def_var dah')
            status = check_nc( nf90_def_var(fid, 'dzv',   NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_dzv), &
                 'nf90_def_var dzv')
            status = check_nc( nf90_def_var(fid, 'drrhv', NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_drrhv), &
                 'nf90_def_var drrhv')
            status = check_nc( nf90_def_var(fid, 'dirhv', NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_dirhv), &
                 'nf90_def_var dirhvm')
            status = check_nc( nf90_def_var(fid, 'dkdp',  NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_dkdp), &
                 'nf90_def_var dkdp')
            status = check_nc( nf90_def_var(fid, 'dadp',  NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_dadp), &
                 'nf90_def_var dadp')
            status = check_nc( nf90_def_var(fid, 'dzvh',  NF90_DOUBLE, (/dimid_qi, nTa_dimid, nTm_dimid/), varid_dzvh), &
                 'nf90_def_var dzvh')


            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'Table_version', TRIM(look_Z_type%cversion_lt)), &
                 'nf90_put_att Table_version')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'Magic_number', look_Z_type%magicnr), &
                 'nf90_put_att Magic_number')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'Setup_string', TRIM(extrastr)), &
                 'nf90_put_att Setup_string')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'itype_refl', look_Z_type%itype_refl), &
                 'nf90_put_att itype_refl')
            IF (look_Z_type%luse_tmatrix) THEN
              status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'luse_tmatrix', 'TRUE'), &
                   'nf90_put_att luse_tmatrix')
            ELSE
              status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'luse_tmatrix', 'FALSE'), &
                   'nf90_put_att luse_tmatrix')
            END IF
            IF (look_Z_type%ldo_nonsphere) THEN
              status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'ldo_nonsphere', 'TRUE'), &
                   'nf90_put_att ldo_nonsphere')
            ELSE
              status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'ldo_nonsphere', 'FALSE'), &
                   'nf90_put_att ldo_nonsphere')
            END IF
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'itype_Dref_fmelt', look_Z_type%itype_Dref_fmelt), &
                 'nf90_put_att itype_Dref_fmelt')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'nparams', look_Z_type%nparams), &
                 'nf90_put_att nparams')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'i_zh', look_Z_type%i_zh), &
                 'nf90_put_att i_zh')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'i_ah', look_Z_type%i_ah), &
                 'nf90_put_att i_ah')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'i_zv', look_Z_type%i_zv), &
                 'nf90_put_att i_zv')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'i_rrhv', look_Z_type%i_rrhv), &
                 'nf90_put_att i_rrhv')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'i_irhv', look_Z_type%i_irhv), &
                 'nf90_put_att i_irhv')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'i_kdp', look_Z_type%i_kdp), &
                 'nf90_put_att i_kdp')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'i_adp', look_Z_type%i_adp), &
                 'nf90_put_att i_adp')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'i_zvh', look_Z_type%i_zvh), &
                 'nf90_put_att i_zvh')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'flag_qi_scal_eq', look_Z_type%flag_qi_scal_eq), &
                 'nf90_put_att flag_qi_scal_eq')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'f_eq', look_Z_type%f_eq), &
                 'nf90_put_att f_eq')
            status = check_nc( nf90_put_att(fid, NF90_GLOBAL, 'if_eq', look_Z_type%if_eq), &
                 'nf90_put_att if_eq')

            status = check_nc( nf90_enddef(fid) , 'nf90_enddef')

            ! 2) Write data to Netcdf File
            ! ----------------------------

            status = check_nc( nf90_put_var(fid, varid_qi, look_Z_type%q_i, start=(/1/)), 'nf90_put_var q_i')
            status = check_nc( nf90_put_var(fid, varid_qi_lin, look_Z_type%q_i_lin, start=(/1/)), 'nf90_put_var q_i_lin')
            status = check_nc( nf90_put_var(fid, varid_Ta, look_Z_type%T_a, start=(/1/)), 'nf90_put_var T_a')
            status = check_nc( nf90_put_var(fid, varid_Tm, look_Z_type%T_m, start=(/1/)), 'nf90_put_var T_m')

            status = check_nc( nf90_put_var(fid, varid_zh, look_Z_type%zh%val, start=(/1,1,1/)), 'nf90_put_var zh')
            status = check_nc( nf90_put_var(fid, varid_ah, look_Z_type%ah%val, start=(/1,1,1/)), 'nf90_put_var ah')
            status = check_nc( nf90_put_var(fid, varid_zv, look_Z_type%zv%val, start=(/1,1,1/)), 'nf90_put_var zv')
            status = check_nc( nf90_put_var(fid, varid_rrhv, look_Z_type%rrhv%val, start=(/1,1,1/)), 'nf90_put_var rrhv')
            status = check_nc( nf90_put_var(fid, varid_irhv, look_Z_type%irhv%val, start=(/1,1,1/)), 'nf90_put_var irhv')
            status = check_nc( nf90_put_var(fid, varid_kdp, look_Z_type%kdp%val, start=(/1,1,1/)), 'nf90_put_var kdp')
            status = check_nc( nf90_put_var(fid, varid_adp, look_Z_type%adp%val, start=(/1,1,1/)), 'nf90_put_var adp')
            status = check_nc( nf90_put_var(fid, varid_zvh, look_Z_type%zvh%val, start=(/1,1,1/)), 'nf90_put_var zvh')

            status = check_nc( nf90_put_var(fid, varid_dzh, look_Z_type%zh%dval, start=(/1,1,1/)), 'nf90_put_var dzh')
            status = check_nc( nf90_put_var(fid, varid_dah, look_Z_type%ah%dval, start=(/1,1,1/)), 'nf90_put_var dah')
            status = check_nc( nf90_put_var(fid, varid_dzv, look_Z_type%zv%dval, start=(/1,1,1/)), 'nf90_put_var dzv')
            status = check_nc( nf90_put_var(fid, varid_drrhv, look_Z_type%rrhv%dval, start=(/1,1,1/)), 'nf90_put_var drrhv')
            status = check_nc( nf90_put_var(fid, varid_dirhv, look_Z_type%irhv%dval, start=(/1,1,1/)), 'nf90_put_var dirhv')
            status = check_nc( nf90_put_var(fid, varid_dkdp, look_Z_type%kdp%dval, start=(/1,1,1/)), 'nf90_put_var dkdp')
            status = check_nc( nf90_put_var(fid, varid_dadp, look_Z_type%adp%dval, start=(/1,1,1/)), 'nf90_put_var dadp')
            status = check_nc( nf90_put_var(fid, varid_dzvh, look_Z_type%zvh%dval, start=(/1,1,1/)), 'nf90_put_var dzvh')

            ! 3) Close Netcdf File
            ! --------------------

            status = check_nc( nf90_close(fid), 'nf90_close '//TRIM(savefile)//'.nc' )
            IF (status /= nf90_noerr) THEN
              WRITE (*,'(a)') 'Error in closing file '//TRIM(savefile)//'.nc '// &
                   TRIM(nf90_strerror(status))
            END IF

          END IF   ! status == nf90_noerr

        END IF
#endif

      ELSE

        IF (llut_format_netcdf) THEN

          WRITE (*,'(a)') 'INFO '//TRIM(yzroutine)// &
               ': Tried to write binary Mie- or T-matrix lookup table to file '//TRIM(savefile)//'.nc '// &
               'but file already exists! Is not overwritten!'

        ELSE

          WRITE (*,'(a)') 'INFO '//TRIM(yzroutine)// &
               ': Tried to write binary Mie- or T-matrix lookup table to file '//TRIM(savefile)//'.bin '// &
               'but file already exists! Is not overwritten!'

        END IF

      END IF   ! .not. fileexist

    END IF   ! my_cart_id_fwo == output_pe

  CONTAINS

    FUNCTION check_nc ( istat, rinfo ) RESULT (ostat)
      INTEGER, INTENT(in)          :: istat
      CHARACTER(len=*), INTENT(in) :: rinfo
      INTEGER                      :: ostat

      ! .. return the error status, so that the below error handling might be
      !     relaxed in the future (e.g., WARNING instead of abort_run())
      !     and the calling routine can decide what to do:
      ostat = istat
      IF (istat /= NF90_NOERR) THEN
        WRITE(*,*) &
             'ERROR write_lookup_to_file() in writing NetCDF: '//TRIM(rinfo)//': '//TRIM(NF90_strerror(istat))
        STOP
      END IF

    END FUNCTION check_nc

  END SUBROUTINE write_lookup_to_file

  !===================================================================================================

  SUBROUTINE read_lookup_from_file(&
      unitnr, savefile, look_Z_type, yzroutine, pe_for_info, ierr)

    IMPLICIT NONE

    INTEGER, INTENT(in)                   :: unitnr, pe_for_info
    CHARACTER(len=*), INTENT(in)          :: savefile, yzroutine
    TYPE(t_dbzlookuptable), INTENT(inout) :: look_Z_type
    INTEGER, INTENT(out) :: ierr

    INTEGER :: nparams_lut, i_zh_lut, i_ah_lut, i_zv_lut, i_rrhv_lut, i_irhv_lut, i_kdp_lut, i_adp_lut, i_zvh_lut

    INTEGER :: fid, status
    INTEGER :: nqi_dimid, nqilin_dimid, nTa_dimid, nTm_dimid
    INTEGER :: varid_qi, varid_qi_lin, varid_Ta, varid_Tm
    INTEGER :: varid_zh, varid_ah
    INTEGER :: varid_dzh, varid_dah
    INTEGER :: varid_zv, varid_rrhv, varid_irhv, varid_kdp, varid_adp, varid_zvh
    INTEGER :: varid_dzv, varid_drrhv, varid_dirhv, varid_dkdp, varid_dadp, varid_dzvh
    
    ierr = 0

    IF (.NOT. llut_format_netcdf) THEN

      IF (ldebug_dbz .OR. my_cart_id_fwo == pe_for_info) THEN
        WRITE (*,'(a,i0)') 'DEBUG '//TRIM(yzroutine)//&
             ': Reading binary Mie- or T-matrix lookup table from file '//TRIM(savefile)//'.bin on proc ', my_cart_id_fwo
      END IF

      OPEN(unitnr, file=TRIM(savefile)//'.bin', status='old', form='unformatted', iostat=ierr)
      IF (ierr /= 0) THEN
        WRITE (*,*) 'ERROR '//TRIM(yzroutine)//': binary lookup-table file ' &
             // TRIM(savefile)// '.bin could not be opened for reading!'
        RETURN
      END IF

      READ(unitnr) look_Z_type%cversion_lt
      READ(unitnr) look_Z_type%magicnr
      READ(unitnr) look_Z_type%nqi
      READ(unitnr) look_Z_type%nTa
      READ(unitnr) look_Z_type%nTm
      READ(unitnr) nparams_lut
      READ(unitnr) i_zh_lut
      READ(unitnr) i_ah_lut
      READ(unitnr) i_zv_lut
      READ(unitnr) i_rrhv_lut
      READ(unitnr) i_irhv_lut
      READ(unitnr) i_kdp_lut
      READ(unitnr) i_adp_lut
      READ(unitnr) i_zvh_lut
      READ(unitnr) look_Z_type%flag_qi_scal_eq
      READ(unitnr) look_Z_type%f_eq
      READ(unitnr) look_Z_type%if_eq

      IF (nparams_lut /= look_Z_type%nparams) THEN
        WRITE (*,'(a,i0,a,i0,a)') 'ERROR '//TRIM(yzroutine)//': binary lookup-table file ' &
             // TRIM(savefile)// '.bin contained wrong nparams=', nparams_lut, &
             ' instead of ',look_Z_type%nparams ,'!'
        ierr = 2
        RETURN
      END IF
      IF ( i_zh_lut /= look_Z_type%i_zh     .OR. i_ah_lut /= look_Z_type%i_ah     .OR. &
           i_zv_lut /= look_Z_type%i_zv     .OR. i_rrhv_lut /= look_Z_type%i_rrhv .OR. &
           i_irhv_lut /= look_Z_type%i_irhv .OR.  i_kdp_lut /= look_Z_type%i_kdp  .OR. &
           i_adp_lut /= look_Z_type%i_adp   .OR.  i_zvh_lut /= look_Z_type%i_zvh ) THEN
        WRITE (*,'(a)') 'ERROR '//TRIM(yzroutine)//': binary lookup-table file ' &
             // TRIM(savefile)// '.bin contained wrong memory block indices for parameters!'
        ierr = 3
        RETURN
      END IF
    
      READ(unitnr) look_Z_type%qi0
      READ(unitnr) look_Z_type%Ta0
      READ(unitnr) look_Z_type%Tm0
      
      READ(unitnr) look_Z_type%dqi
      READ(unitnr) look_Z_type%dTa
      READ(unitnr) look_Z_type%dTm

      READ(unitnr) look_Z_type%idqi
      READ(unitnr) look_Z_type%idTa
      READ(unitnr) look_Z_type%idTm

      READ(unitnr) look_Z_type%q_i
      READ(unitnr) look_Z_type%q_i_lin
      READ(unitnr) look_Z_type%T_a
      READ(unitnr) look_Z_type%T_m

      READ(unitnr) look_Z_type%mem
      READ(unitnr) look_Z_type%dmem

      CLOSE(unitnr)

    ELSE

      IF (ldebug_dbz .OR. my_cart_id_fwo == pe_for_info) THEN
        WRITE (*,'(a,i0)') 'INFO '//TRIM(yzroutine)//&
             ': Reading netcdf Mie- or T-matrix lookup table from file '//TRIM(savefile)//'.nc on proc ', my_cart_id_fwo
      END IF
      
      fid = -1
      status = check_nc( nf90_open (TRIM(savefile)//'.nc', NF90_NOWRITE, fid), 'nf90_open '//TRIM(savefile)//'.nc' )
      
      IF (status /= NF90_NOERR) THEN
        WRITE (*,*) 'ERROR '//TRIM(yzroutine)//': netcdf lookup-table file ' &
             // TRIM(savefile)// '.nc could not be opened for reading!'
        fid = -1
        ierr = 1
        RETURN
      END IF

      status = check_nc( nf90_inq_dimid         (fid, 'q_i_lin', nqilin_dimid), 'nf90_inq_dimid q_i_lin')
      status = check_nc( nf90_inquire_dimension (fid, nqilin_dimid, len = look_Z_type%nqi), 'nf90_inquire_dimension q_i_lin')
      status = check_nc( nf90_inq_dimid         (fid, 'q_i', nqi_dimid), 'nf90_inq_dimid q_i')
      status = check_nc( nf90_inquire_dimension (fid, nqi_dimid, len = look_Z_type%nqi), 'nf90_inquire_dimension q_i')
      status = check_nc( nf90_inq_dimid         (fid, 'T_a', nTa_dimid), 'nf90_inq_dimid T_a')
      status = check_nc( nf90_inquire_dimension (fid, nTa_dimid, len = look_Z_type%nTa), 'nf90_inquire_dimension T_a')
      status = check_nc( nf90_inq_dimid         (fid, 'T_m', nTm_dimid), 'nf90_inq_dimid nTm')
      status = check_nc( nf90_inquire_dimension (fid, nTm_dimid, len = look_Z_type%nTm), 'nf90_inquire_dimension T_m')

      look_Z_type%cversion_lt(:) = ' '
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'Table_version', look_Z_type%cversion_lt), &
           'nf90_get_att Table_version')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'Magic_number', look_Z_type%magicnr), &
           'nf90_get_att Magic_number')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'itype_refl', look_Z_type%itype_refl), &
           'nf90_get_att itype_refl')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'itype_Dref_fmelt', look_Z_type%itype_Dref_fmelt), &
           'nf90_get_att itype_Dref_fmelt')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'nparams', nparams_lut), &
           'nf90_get_att nparams_lut')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'i_zh', i_zh_lut), &
           'nf90_get_att i_zh_lut')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'i_ah', i_ah_lut), &
           'nf90_get_att i_ah_lut')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'i_zv', i_zv_lut), &
           'nf90_get_att i_zv_lut')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'i_rrhv', i_rrhv_lut), &
           'nf90_get_att i_rrhv_lut')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'i_irhv', i_irhv_lut), &
           'nf90_get_att i_irhv_lut')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'i_kdp', i_kdp_lut), &
           'nf90_get_att i_kdp_lut')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'i_adp', i_adp_lut), &
           'nf90_get_att i_adp_lut')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'i_zvh', i_zvh_lut), &
           'nf90_get_att i_zvh_lut')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'flag_qi_scal_eq', look_Z_type%flag_qi_scal_eq), &
           'nf90_get_att flag_qi_scal_eq')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'f_eq', look_Z_type%f_eq), &
           'nf90_get_att f_eq')
      status = check_nc( nf90_get_att(fid, NF90_GLOBAL, 'if_eq', look_Z_type%if_eq), &
           'nf90_get_att if_eq')

      IF (nparams_lut /= look_Z_type%nparams) THEN
        WRITE (*,'(a,i0,a,i0,a)') 'ERROR '//TRIM(yzroutine)//': netcdf lookup-table file ' &
             // TRIM(savefile)// '.nc  contained wrong nparams=', nparams_lut, &
             ' instead of ',look_Z_type%nparams ,'!'
        ierr = 2
        status = check_nc( nf90_close(fid), 'nf90_close '//TRIM(savefile)//'.nc' )
        RETURN
      END IF
      IF ( i_zh_lut /= look_Z_type%i_zh     .OR. i_ah_lut /= look_Z_type%i_ah     .OR. &
           i_zv_lut /= look_Z_type%i_zv     .OR. i_rrhv_lut /= look_Z_type%i_rrhv .OR. &
           i_irhv_lut /= look_Z_type%i_irhv .OR.  i_kdp_lut /= look_Z_type%i_kdp  .OR. &
           i_adp_lut /= look_Z_type%i_adp   .OR.  i_zvh_lut /= look_Z_type%i_zvh ) THEN
        WRITE (*,'(a)') 'ERROR '//TRIM(yzroutine)//': binary lookup-table file ' &
             // TRIM(savefile)// '.bin contained wrong memory block indices for parameters!'
        ierr = 3
        status = check_nc( nf90_close(fid), 'nf90_close '//TRIM(savefile)//'.nc' )
        RETURN
      END IF

      status = check_nc( nf90_inq_varid (fid, 'q_i'   , varid_qi), 'nf90_inq_varid q_i')
      status = check_nc( nf90_get_var   (fid, varid_qi, look_Z_type%q_i, (/1/)), 'nf90_get_var q_i')
      status = check_nc( nf90_get_att (fid, varid_qi, 'qi0', look_Z_type%qi0), 'nf90_get_att qi0 for q_i')
      status = check_nc( nf90_get_att (fid, varid_qi, 'dqi', look_Z_type%dqi), 'nf90_get_att dqi for q_i')
      status = check_nc( nf90_get_att (fid, varid_qi, 'idqi', look_Z_type%idqi), 'nf90_get_att idqi for q_i')

      status = check_nc( nf90_inq_varid (fid, 'q_i_lin'   , varid_qi_lin), 'nf90_inq_varid q_i_lin')
      status = check_nc( nf90_get_var   (fid, varid_qi_lin, look_Z_type%q_i_lin, (/1/)), 'nf90_get_var q_i_lin')

      status = check_nc( nf90_inq_varid (fid, 'T_a'   , varid_Ta), 'nf90_inq_varid T_a')
      status = check_nc( nf90_get_var   (fid, varid_Ta, look_Z_type%T_a, (/1/)), 'nf90_get_var T_a')
      status = check_nc( nf90_get_att (fid, varid_Ta, 'Ta0', look_Z_type%Ta0), 'nf90_get_att Ta0 for T_a')
      status = check_nc( nf90_get_att (fid, varid_Ta, 'dTa', look_Z_type%dTa), 'nf90_get_att dTa for T_a')
      status = check_nc( nf90_get_att (fid, varid_Ta, 'idTa', look_Z_type%idTa), 'nf90_get_att idTa for T_a')

      status = check_nc( nf90_inq_varid (fid, 'T_m'   , varid_Tm), 'nf90_inq_varid T_m')
      status = check_nc( nf90_get_var   (fid, varid_Tm, look_Z_type%T_m, (/1/)), 'nf90_get_var T_m')
      status = check_nc( nf90_get_att (fid, varid_Tm, 'Tm0', look_Z_type%Tm0), 'nf90_get_att Tm0 for T_m')
      status = check_nc( nf90_get_att (fid, varid_Tm, 'dTm', look_Z_type%dTm), 'nf90_get_att dTm for T_m')
      status = check_nc( nf90_get_att (fid, varid_Tm, 'idTm', look_Z_type%idTm), 'nf90_get_att idTm for T_m')


      status = check_nc( nf90_inq_varid (fid, 'zh'   , varid_zh), 'nf90_inq_varid zh')
      status = check_nc( nf90_get_var(fid, varid_zh, look_Z_type%zh%val, start=(/1,1,1/)), 'nf90_get_var zh')
      status = check_nc( nf90_inq_varid (fid, 'ah'   , varid_ah), 'nf90_inq_varid ah')
      status = check_nc( nf90_get_var(fid, varid_ah, look_Z_type%ah%val, start=(/1,1,1/)), 'nf90_get_var ah')
      status = check_nc( nf90_inq_varid (fid, 'zv'   , varid_zv), 'nf90_inq_varid zv')
      status = check_nc( nf90_get_var(fid, varid_zv, look_Z_type%zv%val, start=(/1,1,1/)), 'nf90_get_var zv')
      status = check_nc( nf90_inq_varid (fid, 'rrhv'   , varid_rrhv), 'nf90_inq_varid rrhv')
      status = check_nc( nf90_get_var(fid, varid_rrhv, look_Z_type%rrhv%val, start=(/1,1,1/)), 'nf90_get_var rrhv')
      status = check_nc( nf90_inq_varid (fid, 'irhv'   , varid_irhv), 'nf90_inq_varid irhv')
      status = check_nc( nf90_get_var(fid, varid_irhv, look_Z_type%irhv%val, start=(/1,1,1/)), 'nf90_get_var irhv')
      status = check_nc( nf90_inq_varid (fid, 'kdp'   , varid_kdp), 'nf90_inq_varid kdp')
      status = check_nc( nf90_get_var(fid, varid_kdp, look_Z_type%kdp%val, start=(/1,1,1/)), 'nf90_get_var kdp')
      status = check_nc( nf90_inq_varid (fid, 'adp'   , varid_adp), 'nf90_inq_varid adp')
      status = check_nc( nf90_get_var(fid, varid_adp, look_Z_type%adp%val, start=(/1,1,1/)), 'nf90_get_var adp')
      status = check_nc( nf90_inq_varid (fid, 'zvh'   , varid_zvh), 'nf90_inq_varid zvh')
      status = check_nc( nf90_get_var(fid, varid_zvh, look_Z_type%zvh%val, start=(/1,1,1/)), 'nf90_get_var zvh')

      status = check_nc( nf90_inq_varid (fid, 'dzh'   , varid_dzh), 'nf90_inq_varid dzh')
      status = check_nc( nf90_get_var(fid, varid_dzh, look_Z_type%zh%dval, start=(/1,1,1/)), 'nf90_get_var dzh')
      status = check_nc( nf90_inq_varid (fid, 'dah'   , varid_dah), 'nf90_inq_varid dah')
      status = check_nc( nf90_get_var(fid, varid_dah, look_Z_type%ah%dval, start=(/1,1,1/)), 'nf90_get_var dah')
      status = check_nc( nf90_inq_varid (fid, 'dzv'   , varid_dzv), 'nf90_inq_varid dzv')
      status = check_nc( nf90_get_var(fid, varid_dzv, look_Z_type%zv%dval, start=(/1,1,1/)), 'nf90_get_var dzv')
      status = check_nc( nf90_inq_varid (fid, 'drrhv' , varid_drrhv), 'nf90_inq_varid drrhv')
      status = check_nc( nf90_get_var(fid, varid_drrhv, look_Z_type%rrhv%dval, start=(/1,1,1/)), 'nf90_get_var drrhv')
      status = check_nc( nf90_inq_varid (fid, 'dirhv' , varid_dirhv), 'nf90_inq_varid dirhv')
      status = check_nc( nf90_get_var(fid, varid_dirhv, look_Z_type%irhv%dval, start=(/1,1,1/)), 'nf90_get_var dirhv')
      status = check_nc( nf90_inq_varid (fid, 'dkdp'  , varid_dkdp), 'nf90_inq_varid dkdp')
      status = check_nc( nf90_get_var(fid, varid_dkdp, look_Z_type%kdp%dval, start=(/1,1,1/)), 'nf90_get_var dkdp')
      status = check_nc( nf90_inq_varid (fid, 'dadp'  , varid_dadp), 'nf90_inq_varid dadp')
      status = check_nc( nf90_get_var(fid, varid_dadp, look_Z_type%adp%dval, start=(/1,1,1/)), 'nf90_get_var dadp')
      status = check_nc( nf90_inq_varid (fid, 'dzvh'  , varid_dzvh), 'nf90_inq_varid dzvh')
      status = check_nc( nf90_get_var(fid, varid_dzvh, look_Z_type%zvh%dval, start=(/1,1,1/)), 'nf90_get_var dzvh')

      status = check_nc( nf90_close(fid), 'nf90_close '//TRIM(savefile)//'.nc' )
      IF (status /= nf90_noerr) THEN
        WRITE (*,'(a)') 'Error in closing file '//TRIM(savefile)//'.nc '// &
             TRIM(nf90_strerror(status))
      END IF

    END IF
    
  CONTAINS

    FUNCTION check_nc ( istat, rinfo ) RESULT (ostat)
      INTEGER, INTENT(in)          :: istat
      CHARACTER(len=*), INTENT(in) :: rinfo
      INTEGER                      :: ostat

      ! .. return the error status, so that the below error handling might be
      !     relaxed in the future (e.g., WARNING instead of abort_run())
      !     and the calling routine can decide what to do:
      ostat = istat
      IF (istat /= NF90_NOERR) THEN
        WRITE(*,*) &
             'ERROR read_lookup_from_file() during reading NetCDF: '//TRIM(rinfo)//': '//TRIM(NF90_strerror(istat))
        STOP
      END IF

    END FUNCTION check_nc

  END SUBROUTINE read_lookup_from_file

  !===================================================================================================

  !*******************************************************************************
  !
  ! Vectorized unified version of bi-/triinterpolate_lookup().
  !
  ! Unified code for both, because since LUT version006 all tables are 3-D.
  ! The former 2-D tables (1-mom rain, dry snow/graupel/hail) now contain
  ! only one node in T_m direction, but are 3-D.
  !
  ! The original code of bi-/triinterpolate_lookup() is kept in the code
  !  below for reference, but is not used any more.
  !
  !*******************************************************************************

  SUBROUTINE zradar_triinterp_lookup_add_vec(&
       look_Z_type,q_i,T_a,T_m,flag_interp,&
       ilow,iup,&
       zh_radar_int,ah_radar_int,&
       zv_radar_int,rrhv_radar_int,irhv_radar_int,&
       kdp_radar_int,adp_radar_int,zvh_radar_int,&
       n_i, ext_tune_fac )  ! ,&
!!$ ! in preparation for cubic splines and adjoint operator
!!$       dzh_dqi, dah_dqi, dzv_dqi, drrhv_dqi, dirhv_dqi, dkdp_dqi, dadp_dqi, dzvh_dqi)

    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(in)       :: look_Z_type
    REAL(KIND=dp), INTENT(in), DIMENSION(:)  :: q_i, T_a, T_m  ! q_i in kg/m^3 resp. m; T_a, T_m in K
    LOGICAL, INTENT(in)  :: flag_interp(:)
    INTEGER,                      INTENT(IN) :: ilow, iup

    REAL(KIND=dp), INTENT(inout), DIMENSION(:), TARGET :: &
         &                                      zh_radar_int, ah_radar_int, &
         &                                      zv_radar_int, rrhv_radar_int, irhv_radar_int, &
         &                                      kdp_radar_int, adp_radar_int, zvh_radar_int

    REAL(KIND=dp), INTENT(in), DIMENSION(:), OPTIONAL  :: n_i  ! n_i in 1/m^3; if present, multiply increments by n_i
    REAL(KIND=dp), INTENT(in), OPTIONAL  :: ext_tune_fac  ! tuning factor for attenuation coefficient

!!$ ! in preparation for cubic splines and adjoint operator
!!$    REAL(KIND=dp), INTENT(in), DIMENSION(:), OPTIONAL  :: dzh_dqi, dah_dqi, &
!!$         &                                                dzv_dqi, drrhv_dqi, &
!!$         &                                                dirhv_dqi, dkdp_dqi, &
!!$         &                                                dadp_dqi, dzvh_dqi

#ifdef __SX__
    ! This dummy save variable inhibits inlining of this routine. Inlining
    ! leads to a crash on the NEC in the calling routine when compiled with OpenMP support.
    LOGICAL, SAVE :: firstcall = .TRUE.
#endif

    INTEGER       :: i, itype, ntypes
    INTEGER       :: iu, io

    INTEGER,       DIMENSION(SIZE(q_i,dim=1)) :: ju, jo, ku, ko
    REAL(KIND=dp), DIMENSION(SIZE(q_i,dim=1)) :: qi_eq, qi_loc, qi_scal_loc, Ta_du, Ta_do, Tm_du, Tm_do, val_int
!!$    REAL(KIND=dp), DIMENSION(SIZE(q_i,dim=1)) :: dval_int

    REAL(kind=dp) :: Ta_loc, Tm_loc,         &
         &           dqi_loc, idqi_loc, qi_du, qi_do, &
         &           idTa_loc,               &
         &           idTm_loc
    REAL(kind=dp) :: val_int_u, val_int_o, dval_int_u, dval_int_o
    REAL(kind=dp) :: t, c0, c1, c2, c3

    REAL(KIND=dp), POINTER     :: val(:,:,:), dval(:,:,:), qi_for_val(:)
    

    TYPE t_val_list
      REAL(kind=dp), POINTER :: val(:)   ! pointer to output parameter vector
      REAL(kind=dp), POINTER :: dval(:)  ! pointer to vector of derivative of output parameter
    END TYPE t_val_list
    TYPE(t_val_list) :: p_val (look_Z_type%nparams)
    
    IF (look_Z_type%luse_tmatrix) THEN
      ! loop over all polarimetric parameters:
      ntypes = look_Z_type%nparams
    ELSE
      ! only i_zh and i_ah:
      ntypes = 2
    END IF

    ! Associate pointers to the output fields, so that we can loop over parameters below:
    p_val(look_Z_type%i_zh)%val => zh_radar_int
    p_val(look_Z_type%i_ah)%val => ah_radar_int
    p_val(look_Z_type%i_zv)%val => zv_radar_int
    p_val(look_Z_type%i_rrhv)%val => rrhv_radar_int
    p_val(look_Z_type%i_irhv)%val => irhv_radar_int
    p_val(look_Z_type%i_kdp)%val => kdp_radar_int
    p_val(look_Z_type%i_adp)%val => adp_radar_int
    p_val(look_Z_type%i_zvh)%val => zvh_radar_int
    
!!$ ! in preparation for cubic splines and adjoint operator
!!$    IF (PRESENT(dzh_dqi)) THEN
!!$      ALLOCATE(dval_int(ni)
!!$      p_val(look_Z_type%i_zh)%dval => dzh_dqi
!!$      p_val(look_Z_type%i_ah)%dval => dah_dqi
!!$      p_val(look_Z_type%i_zv)%dval => dzv_dqi
!!$      p_val(look_Z_type%i_rrhv)%dval => drrhv_dqi
!!$      p_val(look_Z_type%i_irhv)%dval => dirhv_dqi
!!$      p_val(look_Z_type%i_kdp)%dval => dkdp_dqi
!!$      p_val(look_Z_type%i_adp)%dval => dadp_dqi
!!$      p_val(look_Z_type%i_zvh)%dval => dzvh_dqi
!!$    END IF
    
    DO i = ilow, iup

      ! Truncate input values of q_i to the range of the table:
      IF (q_i(i) >= look_Z_type%q_i_lin(1) .AND. flag_interp(i)) THEN
        qi_loc(i) = MAX(MIN(q_i(i), look_Z_type%q_i_lin(look_Z_type%nqi)), look_Z_type%q_i_lin(1))

        ! Truncate input values of T_a to the range of the table:
        Ta_loc = MAX(MIN(T_a(i), look_Z_type%T_a(look_Z_type%nTa)), look_Z_type%Ta0)
                
        ! Calculate indices of the neighbouring regular T_a-values in the table:
        ju(i) = MAX(MIN(FLOOR((Ta_loc - look_Z_type%Ta0) * look_Z_type%idTa)+1, look_Z_type%nTa-1),1)
        jo(i) = ju(i) + 1

        ! Calculate interpolation weights (relative differences to the lower and upper coordinates):
        idTa_loc = 1.0_dp / (look_Z_type%T_a(jo(i)) - look_Z_type%T_a(ju(i)))
        Ta_du(i) = (Ta_loc      - look_Z_type%T_a(ju(i)))*idTa_loc
        Ta_do(i) = (look_Z_type%T_a(jo(i)) - Ta_loc     )*idTa_loc
        
      ELSE
        qi_loc(i) = -HUGE(1.0_dp)
      END IF

    END DO

    
    IF (look_Z_type%nTm > 1) THEN
      DO i = ilow, iup
        IF (q_i(i) >= look_Z_type%q_i_lin(1) .AND. flag_interp(i)) THEN
          ! Truncate input values of T_m to the range of the table:
          Tm_loc = MAX(MIN(T_m(i), look_Z_type%T_m(look_Z_type%nTm)), look_Z_type%Tm0)
          
          ! Calculate indices of the neighbouring regular T_m-values in the table:
          ku(i) = MAX(MIN(FLOOR((Tm_loc - look_Z_type%Tm0) * look_Z_type%idTm)+1, &
               look_Z_type%nTm-1),1)
          ko(i) = ku(i) + 1

          ! Calculate interpolation weights (relative differences to the lower and upper coordinates):
          idTm_loc = 1.0_dp / (look_Z_type%T_m(ko(i)) - look_Z_type%T_m(ku(i)))
          Tm_du(i) = (Tm_loc  - look_Z_type%T_m(ku(i)))*idTm_loc
          Tm_do(i) = (look_Z_type%T_m(ko(i)) - Tm_loc )*idTm_loc
        END IF
      END DO
    ELSE
      ! Neutral values for neighbouring indices and interpolation weights:
      ku(:) = 1
      ko(:) = 1
      Tm_du(:) = 1.0_dp
      Tm_do(:) = 0.0_dp
    END IF

        
    ! Scale qi_loc (linear value) to the scaling of the equidistant grid:
    qi_eq(:) = scale_val_lut (qi_loc(:), look_Z_type%flag_qi_scal_eq, look_Z_type%f_eq)

    DO itype = 1, ntypes

      val        => look_Z_type%mem  (:,:,:,itype)
      dval       => look_Z_type%dmem (:,:,:,itype)
      qi_for_val => look_Z_type%qmem  (:,itype)
          
      ! scale qi_loc to the scaling of the qi table vector:
      qi_scal_loc(:) = scale_val_lut (qi_loc(:), look_Z_type%flag_qi_scal(itype), look_Z_type%f_scal_qi(itype))

      SELECT CASE (look_Z_type%flag_mem_interp_qi(itype))
      CASE (i_interp_lin)

        ! Trilinear interpolation of scaled parameter with respect to scaled q_i, T_a and T_m:
        !-------------------------------------------------------------------------------------
        
        DO i = ilow, iup

          IF (qi_loc(i) > 0.0_dp) THEN

            ! calculate indices of the neighbouring regular q_i-nodes in the table:
            !  (the spacing of the q_i-nodes is regular, but not necessarily in the
            !   same scaling than the q_i-values for the actual table interplation)
            iu = MAX(MIN(FLOOR((qi_eq(i) - look_Z_type%qi0) * look_Z_type%idqi)+1, look_Z_type%nqi-1),1)
            io = iu + 1

            ! calculate local internodal distance in the given scaling of the node values:
            idqi_loc = 1.0_dp / (qi_for_val(io) - qi_for_val(iu))
      
            ! calculate difference to the lower coordinate
            qi_du = (qi_scal_loc(i) - qi_for_val(iu))*idqi_loc

            ! calculate difference to the upper coordinate
            qi_do = (qi_for_val(io) - qi_scal_loc(i))*idqi_loc

            ! trilinear interpolation in one go of val(:,:,:):
            val_int(i) = ((val(iu,ju(i),ku(i))*qi_do + val(io,ju(i),ku(i))*qi_du)*Ta_do(i) + &
                 &        (val(iu,jo(i),ku(i))*qi_do + val(io,jo(i),ku(i))*qi_du)*Ta_du(i))*Tm_do(i) + &
                 &       ((val(iu,ju(i),ko(i))*qi_do + val(io,ju(i),ko(i))*qi_du)*Ta_do(i) + &
                 &        (val(iu,jo(i),ko(i))*qi_do + val(io,jo(i),ko(i))*qi_du)*Ta_du(i))*Tm_du(i)

!!$ print-output for control plots to inspect the fits:
!!$PRINT '(a,1x,i3,8(1x,es14.5))', 'ULI '//TRIM(look_Z_type%chydrotype), itype, qi_for_val(iu), val(iu,ju(i),ku(i)), 0.0, &
!!$                 qi_for_val(io), val(io,jo(i),ko(i)), 0.0, &
!!$                 qi_scal_loc(i), val_int(i)
            
!!$ ! in preparation for cubic splines and adjoint operator
!!$          IF (PRESENT(dzh_dqi)) THEN
!!$            dval_int(i) = ...   ! Actual computation of derivative is still missing! Use cubic spline, see below!
!!$          END IF

          ELSE
            
            val_int (i) = zero_value_lut(look_Z_type%flag_mem_scal(itype))
!!$ ! in preparation for cubic splines and adjoint operator
!!$            dval_int(i) = 0.0_dp

          END IF

        END DO

      CASE (i_interp_cubic)

        ! Bilinear interpolation of scaled parameter with respect to T_a and T_m,
        !  followed by cubic interpolation with respect to scaled q_i
        !-----------------------------------------------------------------------

        DO i = ilow, iup

          IF (qi_loc(i) > 0.0_dp) THEN

            ! calculate indices of the neighbouring regular q_i-nodes in the table:
            !  (the spacing of the q_i-nodes is regular, but not necessarily in the
            !   same scaling than the q_i-values for the actual table interplation)
            iu = MAX(MIN(FLOOR((qi_eq(i) - look_Z_type%qi0) * look_Z_type%idqi)+1, look_Z_type%nqi-1),1)
            io = iu + 1

            ! 1) bilinear interpolation w.r.t Ta and Tm for values and derivatives at iu and io:
            val_int_u = (val(iu,ju(i),ku(i))*Ta_do(i) + val(iu,jo(i),ku(i))*Ta_du(i))*Tm_do(i) + &
                 &      (val(iu,ju(i),ko(i))*Ta_do(i) + val(iu,jo(i),ko(i))*Ta_du(i))*Tm_du(i)
            
            val_int_o = (val(io,ju(i),ku(i))*Ta_do(i) + val(io,jo(i),ku(i))*Ta_du(i))*Tm_do(i) + &
                 &      (val(io,ju(i),ko(i))*Ta_do(i) + val(io,jo(i),ko(i))*Ta_du(i))*Tm_du(i)

!!$ Das ist bisher dp/dqi, nicht df(p)/g(q)!!!
            dval_int_u = (dval(iu,ju(i),ku(i))*Ta_do(i) + dval(iu,jo(i),ku(i))*Ta_du(i))*Tm_do(i) + &
                 &       (dval(iu,ju(i),ko(i))*Ta_do(i) + dval(iu,jo(i),ko(i))*Ta_du(i))*Tm_du(i)
            
            dval_int_o = (dval(io,ju(i),ku(i))*Ta_do(i) + dval(io,jo(i),ku(i))*Ta_du(i))*Tm_do(i) + &
                 &       (dval(io,ju(i),ko(i))*Ta_do(i) + dval(io,jo(i),ko(i))*Ta_du(i))*Tm_du(i)
            

            ! 2) Cubic Hermite spline in a segment between 2 points (x0, p0, dp0) and (x1, p1, dp1): 
            !
            ! f(x) = c3*t^3 + c2*t^2 + c1*t + c0 
            ! with:
            !      t = (x-x0) / d
            !      c0 = p0
            !      c1 = dp0*d
            !      c2 = -3*p0 + 3*p1 - 2*dp0*d - dp1*d
            !      c3 =  2*p0 - 2*p1 +   dp0*d + dp1*d
            !      d  = x1 - x0

            ! calculate local internodal distance in the given scaling of the node values:
            dqi_loc = qi_for_val(io) - qi_for_val(iu)
      
            ! calculate difference to the lower coordinate:
            t = (qi_scal_loc(i) - qi_for_val(iu)) / dqi_loc

            ! coefficients of the Hermite spline:
            c0 = val_int_u
            c1 = dval_int_u * dqi_loc
            c2 = -3.0_dp*val_int_u + 3.0_dp*val_int_o - 2.0_dp*dval_int_u*dqi_loc - dval_int_o*dqi_loc
            c3 =  2.0_dp*val_int_u - 2.0_dp*val_int_o +        dval_int_u*dqi_loc + dval_int_o*dqi_loc

            val_int(i) = c3*t**3 + c2*t**2 + c1*t + c0 

!!$ print-output for control plots to inspect the fits:
!!$PRINT '(a,1x,i3,8(1x,es14.5))', 'ULI '//TRIM(look_Z_type%chydrotype), itype, qi_for_val(iu), val_int_u, dval_int_u, &
!!$                 qi_for_val(io), val_int_o, dval_int_o, &
!!$                 qi_scal_loc(i), val_int(i)
            
!!$ ! in preparation for cubic splines and adjoint operator
!!$          IF (PRESENT(dzh_dqi)) THEN
!!$            dval_int(i) = 3.0_dp*c3*t**2 + 2.0_dp*c2*t + c1
!!$          END IF

          ELSE
            
            val_int (i) = zero_value_lut(look_Z_type%flag_mem_scal(itype))
!!$ ! in preparation for cubic splines and adjoint operator
!!$            dval_int(i) = 0.0_dp

          END IF

        END DO

      CASE default
        
      END SELECT
      
      ! Re-scale parameter to linear range:
      val_int(:) = descale_val_lut (val_int(:), look_Z_type%flag_mem_scal(itype), look_Z_type%if_scal_mem(itype))

!!$ ! in preparation for cubic splines and adjoint operator
!!$      dval_int(:) = ...  ! computation of re-scaled derivative still missing!

      ! Update INOUT parameters:

      IF (PRESENT(ext_tune_fac)) THEN
        IF (ANY(itype == [look_Z_type%i_ah, look_Z_type%i_adp])) THEN
          ! to handwavingly tune the simulated extinction to reproduce the behaviour
          ! of an imperfect, conservative attenuation correction in comparisons
          ! with attenuation corrected observations, we multiply the tuning factor
          ! ext_tune_fac with values between 0 and 1:
          val_int(:) = val_int(:) * ext_tune_fac
        END IF
      END IF
      
      IF (PRESENT(n_i)) THEN
        p_val(itype)%val(:) = p_val(itype)%val(:) + n_i(:) * val_int(:)
!!$ ! in preparation for cubic splines and adjoint operator
!!$        IF (PRESENT(dzh_dqi)) THEN
!!$          p_val(itype)%dval(:) = ...
!!$        END IF
      ELSE
        p_val(itype)%val(:) = p_val(itype)%val(:) + val_int(:)
!!$ ! in preparation for cubic splines and adjoint operator
!!$        IF (PRESENT(dzh_dqi)) THEN
!!$          p_val(itype)%dval(:) = ...
!!$        END IF
      END IF

    END DO  ! itype

  END SUBROUTINE zradar_triinterp_lookup_add_vec


  !*******************************************************************************
  !
  ! NOT USED ANY MORE BUT KEPT IN THE CODE FOR REFERENCE!
  !
  ! Retrieve values of radar reflectivity and its extinction
  ! for not melting hydrometeors (rain, dry snow, dry graupel,...)
  ! by bilinear interpolation from a lookup table as function of
  ! T_a and q_i(=q_r,q_s,q_g,q_h) for which the lookup table has been
  ! created from the output of the subroutine zradar_*_lookupcreate
  ! (* = rain, drysnow, drygraupel, dryhail)
  !
  ! look_Z_type is 3D, but interpolation w.r.t. T_m will not be done.
  ! The third index into the table for interpolation is chosen to be look_Z_type%nTm.
  !
  !*******************************************************************************

  SUBROUTINE zradar_biinterpolate_lookup(&
      look_Z_type,q_i,T_a,&
      zh_radar_int,ah_radar_int,&
      zv_radar_int,rrhv_radar_int,irhv_radar_int,&
      kdp_radar_int,adp_radar_int,zvh_radar_int)

    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(in) :: look_Z_type
    REAL(KIND=dp), INTENT(in)  :: T_a, q_i  ! original values in linear scaling

    REAL(KIND=dp), INTENT(out) :: zh_radar_int, ah_radar_int, &
                                  zv_radar_int, rrhv_radar_int, irhv_radar_int, &
                                  kdp_radar_int, adp_radar_int, zvh_radar_int

    INTEGER       :: iu, io, ju, jo
    REAL(KIND=dp) :: Ta_loc, qi_loc, idqi_loc, qi_eq


    IF (q_i >= look_Z_type%q_i_lin(1)) THEN

      ! Truncate q_i-values and T_a-values to the range of the table:
      qi_loc = MAX(MIN(q_i, look_Z_type%q_i_lin(look_Z_type%nqi)), look_Z_type%q_i_lin(1))
      Ta_loc = MAX(MIN(T_a, look_Z_type%T_a(look_Z_type%nTa)), look_Z_type%Ta0)

      ! Scale qi_loc (linear value) to the scaling of the equidistant grid:
      qi_eq = scale_val_lut (qi_loc, look_Z_type%flag_qi_scal_eq, look_Z_type%f_eq)
      
      ! calculate indices of the neighbouring regular q_i-nodes in the table:
      !  (the spacing of the q_i-nodes is regular, but not necessarily in the
      !   same scaling than the q_i-values for the actual table interplation)
      iu = MAX(MIN(FLOOR((qi_eq - look_Z_type%qi0) * look_Z_type%idqi)+1, &
                   look_Z_type%nqi-1),1)
      io = iu + 1

      ! calculate indices of the neighbouring regular T_a-values in the table:
      ju = MAX(MIN(FLOOR((Ta_loc - look_Z_type%Ta0) * look_Z_type%idTa)+1, &
                   look_Z_type%nTa-1),1)
      jo = ju + 1

      ! bilinear interpolation and re-scaling to linear values of ...

      ! ... a) reflectivity horiz. polar.
      CALL interp2d_lut(look_Z_type%zh, look_Z_type, qi_loc, Ta_loc, iu, io, ju, jo, look_Z_type%nTm, zh_radar_int)

      ! ... b) extinction coeff. horizontal polar.
      CALL interp2d_lut(look_Z_type%ah, look_Z_type, qi_loc, Ta_loc, iu, io, ju, jo, look_Z_type%nTm, ah_radar_int)

      IF (look_Z_type%luse_tmatrix) THEN

        ! ... c) reflectivity vert. polar.
        CALL interp2d_lut(look_Z_type%zv, look_Z_type, qi_loc, Ta_loc, iu, io, ju, jo, look_Z_type%nTm, zv_radar_int)

        ! ... d1) RHOHV real part
        CALL interp2d_lut(look_Z_type%rrhv, look_Z_type, qi_loc, Ta_loc, iu, io, ju, jo, look_Z_type%nTm, rrhv_radar_int)

        ! ... d2) RHOHV imag. part
        CALL interp2d_lut(look_Z_type%irhv, look_Z_type, qi_loc, Ta_loc, iu, io, ju, jo, look_Z_type%nTm, irhv_radar_int)

        ! ... e) KDP
        CALL interp2d_lut(look_Z_type%kdp, look_Z_type, qi_loc, Ta_loc, iu, io, ju, jo, look_Z_type%nTm, kdp_radar_int)

        ! ... f) differential attenuation
        CALL interp2d_lut(look_Z_type%adp, look_Z_type, qi_loc, Ta_loc, iu, io, ju, jo, look_Z_type%nTm, adp_radar_int)

        ! ... g) cross-polar reflectivity (hor-turned-vertical)
        CALL interp2d_lut(look_Z_type%zvh, look_Z_type, qi_loc, Ta_loc, iu, io, ju, jo, look_Z_type%nTm, zvh_radar_int)

      ELSE

        zv_radar_int   = zh_radar_int
        rrhv_radar_int = 0.0_dp
        irhv_radar_int = 0.0_dp
        kdp_radar_int  = 0.0_dp
        adp_radar_int  = 0.0_dp
        zvh_radar_int  = 0.0_dp

      END IF

    ELSE

      zh_radar_int   = 0.0_dp
      ah_radar_int   = 0.0_dp
      zv_radar_int   = 0.0_dp
      rrhv_radar_int = 0.0_dp
      irhv_radar_int = 0.0_dp
      kdp_radar_int  = 0.0_dp
      adp_radar_int  = 0.0_dp
      zvh_radar_int  = 0.0_dp

    END IF

  END SUBROUTINE zradar_biinterpolate_lookup

  SUBROUTINE interp2d_lut(tval, lut, qi_loc, Ta_loc, iu, io, ju, jo, k, val_int, dval_int)

    TYPE (t_tabparam),       INTENT(in)  :: tval
    TYPE (t_dbzlookuptable), INTENT(in)  :: lut   ! only needed for lut%T_a(:)
    REAL(kind=dp),           INTENT(in)  :: qi_loc, Ta_loc
    INTEGER,                 INTENT(in)  :: iu, io, ju, jo, k

    REAL(kind=dp),           INTENT(out) :: val_int   ! interpolated value
    REAL(kind=dp), OPTIONAL, INTENT(out) :: dval_int  ! partial derivative w.r.t. scaled qi (PREPARATION FOR ADJOINT OPERATOR)

    REAL(kind=dp)                        :: qi_scal_loc, &
         &                                  idqi_loc, qi_du, qi_do, &
         &                                  idTa_loc, Ta_du, Ta_do

    ! scale value of qi to the scaling of the qi table vector:
    qi_scal_loc = scale_val_lut (qi_loc, tval%qi_scal, tval%f_qi)
    
    ! calculate local internodal distance in the given scaling of the node values:
    idqi_loc = 1.0_dp / (tval%q_i(io) - tval%q_i(iu))
    idTa_loc = 1.0_dp / (lut %T_a(jo) - lut %T_a(ju))
      
    ! calculate difference to the lower coordinate
    qi_du = (qi_scal_loc - tval%q_i(iu))*idqi_loc
    Ta_du = (Ta_loc      - lut %T_a(ju))*idTa_loc

    ! calculate difference to the upper coordinate
    qi_do = (tval%q_i(io) - qi_scal_loc)*idqi_loc
    Ta_do = (lut %T_a(jo) - Ta_loc     )*idTa_loc

    ! bilinear interpolation in one go of tval%val(:,:,k):
    val_int = (tval%val(iu,ju,k)*qi_do + tval%val(io,ju,k)*qi_du)*Ta_do + &
         &    (tval%val(iu,jo,k)*qi_do + tval%val(io,jo,k)*qi_du)*Ta_du
    
    ! Re-scale parameter to linear range:
    val_int = descale_val_lut (val_int, tval%val_scal, tval%if_val)
    
    IF (PRESENT(dval_int)) THEN
      dval_int = 0.0_dp  ! ADD HERE COMPUTATION OF DERIVATIVE W.R.T. qi_scal_loc FOR ADJOINT OPERATOR
    END IF

    ! More accurate would be cubic interpolation w.r.t. qi. This will require tval%dval(:,:,:) derivatives
    ! and computation of a 3rd order interpolating polynome after bilinear interpolation w.r.t. T_a and T_m.
    ! This might no longer be invariant against interchanging the interpolation order, but I'm not sure.

    ! Cubic Hermite spline in a segment between 2 points (x0, p0, dp0) and (x1, p1, dp1): 
    !
    ! f(x) = c3*t^3 + c2*t^2 + c1*t + c0 
    ! mit:
    !      t = (x-x0) / d
    !      c0 = p0
    !      c1 = dp0*d
    !      c2 = -3*p0 + 3*p1 - 2*dp0*d - dp1*d
    !      c3 = 2*p1 - 2*p1 + dp0*d + dp1*d
    !      d  = x1 - x0

    
  END SUBROUTINE interp2d_lut

  !*******************************************************************************
  !
  ! NOT USED ANY MORE BUT KEPT IN THE CODE FOR REFERENCE!
  !
  ! Retrieve values of radar reflectivity and its extinction
  ! for melting hydrometeors (meltsnow, meltgraupel,...)
  ! by trilinear interpolation from a lookup table as function of
  ! T_a and q_i(=q_r,q_s,q_g,q_h) and T_m for which the lookup table has been
  ! created from the output of the subroutine zradar_*_lookupcreate
  ! (* = meltsnow, meltgraupel, melthail)
  !
  ! look_Z_type%nTm has to be larger than 1, otherwise a crash will happen!
  !
  !*******************************************************************************

  SUBROUTINE zradar_triinterpolate_lookup(&
      look_Z_type,q_i,T_a,T_m,&
      zh_radar_int,ah_radar_int,&
      zv_radar_int,rrhv_radar_int,irhv_radar_int,&
      kdp_radar_int,adp_radar_int,zvh_radar_int)

    IMPLICIT NONE

    TYPE(t_dbzlookuptable), INTENT(in) :: look_Z_type
    REAL(KIND=dp), INTENT(in)  :: T_a, q_i, T_m  ! original values in linear scaling

    REAL(KIND=dp), INTENT(out) :: zh_radar_int, ah_radar_int, &
                                  zv_radar_int, rrhv_radar_int, irhv_radar_int, &
                                  kdp_radar_int, adp_radar_int, zvh_radar_int

    INTEGER       :: iu, io, ju, jo, ku, ko
    REAL(KIND=dp) :: Ta_loc, qi_loc, Tm_loc, qi_eq

    IF (q_i >= look_Z_type%q_i_lin(1)) THEN

      ! Truncate given values of q_i, T_a and T_m to the range of the table:
      qi_loc = MAX(MIN(q_i, look_Z_type%q_i_lin(look_Z_type%nqi)), look_Z_type%q_i_lin(1))
      Ta_loc = MAX(MIN(T_a, look_Z_type%T_a(look_Z_type%nTa)), look_Z_type%Ta0)
      Tm_loc = MAX(MIN(T_m, look_Z_type%T_m(look_Z_type%nTm)), look_Z_type%Tm0)

      ! Scale qi_loc (linear value) to the scaling of the equidistant grid:
      qi_eq = scale_val_lut (qi_loc, look_Z_type%flag_qi_scal_eq, look_Z_type%f_eq)
      
      ! calculate indices of the neighbouring regular q_i-values in the table:
      !  (the spacing of the q_i-nodes is regular, but not necessarily in the
      !   same scaling than the q_i-values for the actual table interplation)
      iu = MAX(MIN(FLOOR((qi_eq - look_Z_type%qi0) * look_Z_type%idqi)+1, &
                   look_Z_type%nqi-1),1)
      io = iu + 1

      ! calculate indices of the neighbouring regular T_a-values in the table:
      ju = MAX(MIN(FLOOR((Ta_loc - look_Z_type%Ta0) * look_Z_type%idTa)+1, &
                   look_Z_type%nTa-1),1)
      jo = ju + 1

      ! calculate indices of the neighbouring regular T_m-values in the table:
      ku = MAX(MIN(FLOOR((Tm_loc - look_Z_type%Tm0) * look_Z_type%idTm)+1, &
                   look_Z_type%nTm-1),1)
      ko = ku + 1

      ! Interpolation of ...
      
      ! ... a) reflectivity horiz. polar.
      CALL interp3d_lut(look_Z_type%zh, look_Z_type, qi_loc, Ta_loc, Tm_loc, iu, io, ju, jo, ku, ko, zh_radar_int)

      ! ... b) extinction coeff. horizontal polar.
      CALL interp3d_lut(look_Z_type%ah, look_Z_type, qi_loc, Ta_loc, Tm_loc, iu, io, ju, jo, ku, ko, ah_radar_int)

      IF (look_Z_type%luse_tmatrix) THEN

        ! ... c) reflectivity vert. polar.
        CALL interp3d_lut(look_Z_type%zv, look_Z_type, qi_loc, Ta_loc, Tm_loc, iu, io, ju, jo, ku, ko, zv_radar_int)

        ! ... d1) RHOHV real part
        CALL interp3d_lut(look_Z_type%rrhv, look_Z_type, qi_loc, Ta_loc, Tm_loc, iu, io, ju, jo, ku, ko, rrhv_radar_int)
        
        ! ... d2) RHOHV imag. part
        CALL interp3d_lut(look_Z_type%irhv, look_Z_type, qi_loc, Ta_loc, Tm_loc, iu, io, ju, jo, ku, ko, irhv_radar_int)
        
        ! ... e) KDP
        CALL interp3d_lut(look_Z_type%kdp, look_Z_type, qi_loc, Ta_loc, Tm_loc, iu, io, ju, jo, ku, ko, kdp_radar_int)
        
        ! ... f) differential attenuation
        CALL interp3d_lut(look_Z_type%adp, look_Z_type, qi_loc, Ta_loc, Tm_loc, iu, io, ju, jo, ku, ko, adp_radar_int)
        
        ! ... g) cross-polar reflectivity (hor-turned-vertical)
        CALL interp3d_lut(look_Z_type%zvh, look_Z_type, qi_loc, Ta_loc, Tm_loc, iu, io, ju, jo, ku, ko, zvh_radar_int)

      ELSE

        zv_radar_int   = zh_radar_int
        rrhv_radar_int = 0.0_dp
        irhv_radar_int = 0.0_dp
        kdp_radar_int  = 0.0_dp
        adp_radar_int  = 0.0_dp
        zvh_radar_int  = 0.0_dp

      END IF

    ELSE

      zh_radar_int   = 0.0_dp
      ah_radar_int   = 0.0_dp
      zv_radar_int   = 0.0_dp
      rrhv_radar_int = 0.0_dp
      irhv_radar_int = 0.0_dp
      kdp_radar_int  = 0.0_dp
      adp_radar_int  = 0.0_dp
      zvh_radar_int  = 0.0_dp

    END IF

  END SUBROUTINE zradar_triinterpolate_lookup

  SUBROUTINE interp3d_lut(tval, lut, qi_loc, Ta_loc, Tm_loc, iu, io, ju, jo, ku, ko, val_int, dval_int)

    TYPE (t_tabparam),       INTENT(in)  :: tval
    TYPE (t_dbzlookuptable), INTENT(in)  :: lut   ! only needed for lut%T_a(:) and lut%T_m(:)
    REAL(kind=dp),           INTENT(in)  :: qi_loc, Ta_loc, Tm_loc
    INTEGER,                 INTENT(in)  :: iu, io, ju, jo, ku, ko

    REAL(kind=dp),           INTENT(out) :: val_int   ! interpolated value
    REAL(kind=dp), OPTIONAL, INTENT(out) :: dval_int  ! partial derivative w.r.t. scaled qi (PREPARATION FOR ADJOINT OPERATOR)

    REAL(kind=dp)                        :: qi_scal_loc, &
         &                                  idqi_loc, qi_du, qi_do, &
         &                                  idTa_loc, Ta_du, Ta_do, &
         &                                  idTm_loc, Tm_du, Tm_do

    ! scale value of qi to the scaling of the qi table vector:
    qi_scal_loc = scale_val_lut (qi_loc, tval%qi_scal, tval%f_qi)
    
    ! calculate local internodal distance in the given scaling of the node values:
    idqi_loc = 1.0_dp / (tval%q_i(io) - tval%q_i(iu))
    idTa_loc = 1.0_dp / (lut %T_a(jo) - lut %T_a(ju))
    idTm_loc = 1.0_dp / (lut %T_m(ko) - lut %T_m(ku))
      
    ! calculate difference to the lower coordinate
    qi_du = (qi_scal_loc - tval%q_i(iu))*idqi_loc
    Ta_du = (Ta_loc      - lut %T_a(ju))*idTa_loc
    Tm_du = (Tm_loc      - lut %T_m(ku))*idTm_loc

    ! calculate difference to the upper coordinate
    qi_do = (tval%q_i(io) - qi_scal_loc)*idqi_loc
    Ta_do = (lut %T_a(jo) - Ta_loc     )*idTa_loc
    Tm_do = (lut %T_m(ko) - Tm_loc     )*idTm_loc

    ! trilinear interpolation in one go of tval%val(:,:,:):
    val_int = ((tval%val(iu,ju,ku)*qi_do + tval%val(io,ju,ku)*qi_du)*Ta_do + &
         &     (tval%val(iu,jo,ku)*qi_do + tval%val(io,jo,ku)*qi_du)*Ta_du)*Tm_do + &
         &    ((tval%val(iu,ju,ko)*qi_do + tval%val(io,ju,ko)*qi_du)*Ta_do + &
         &     (tval%val(iu,jo,ko)*qi_do + tval%val(io,jo,ko)*qi_du)*Ta_du)*Tm_du
    
    ! Re-scale parameter to linear range:
    val_int = descale_val_lut (val_int, tval%val_scal, tval%if_val)
    
    IF (PRESENT(dval_int)) THEN
      dval_int = 0.0_dp  ! ADD HERE COMPUTATION OF DERIVATIVE W.R.T. qi_scal_loc FOR ADJOINT OPERATOR
    END IF

    ! More accurate would be cubic interpolation w.r.t. qi. This will require tval%dval(:,:,:) derivatives
    ! and computation of a 3rd order interpolating polynome after bilinear interpolation w.r.t. T_a and T_m.
    ! This might no longer be invariant against interchanging the interpolation order, but I'm not sure.

    ! Cubic Hermite spline in a segment between 2 points (x0, p0, dp0) and (x1, p1, dp1): 
    !
    ! f(x) = c3*t^3 + c2*t^2 + c1*t + c0 
    ! mit:
    !      t = (x-x0) / d
    !      c0 = p0
    !      c1 = dp0*d
    !      c2 = -3*p0 + 3*p1 - 2*dp0*d - dp1*d
    !      c3 = 2*p1 - 2*p1 + dp0*d + dp1*d
    !      d  = x1 - x0
    
  END SUBROUTINE interp3d_lut

  ! end lookup tables for radar reflectivity (subroutines for creating as well as interpolation)
  !=========================================================================================
  !=========================================================================================


  !----------------------------------------------------------------------------------------!
  !  Calculate effective radar reflectivity using Rayleigh approximation assuming a        !
  !  modified gamma distribution as a function of the gamma distribution parameters.       !
  !  Input: mgd    ( parameters n0 & lam, specifically )                                   !
  !         lamexp ( exponent to lam, f(mu,nu,b) )                                         !
  !         Zfac   ( prefactor, f(mu,nu,b,K,Kw0,rho) )                                     !
  !                                                                                        !
  !  Version works for any hydrometeor of any shape and any moment scheme (these           !
  !  dependencies are imprinted in the input parameters).                                  !
  !----------------------------------------------------------------------------------------!

  FUNCTION zradar_rayleigh_mgd(mgd,Zfac,lamexp) &
      RESULT (Z)

    IMPLICIT NONE

    REAL(KIND=dp),      INTENT(in) :: Zfac,lamexp
    TYPE(t_mgd_params), INTENT(in) :: mgd

    REAL(KIND=dp) :: Z

    ! Results from general Rayleigh-Ze equation:
    ! Ze = K^2/Kw0^2 * n0/nu * (6a/(Pi*rho_hydromedium))^2
    !                                   * gam((mu+1-2b)/nu) / lam^((mu+1-2b)/nu)
    !    = K^2/Kw0^2 * (6a/(Pi*rho_hydromedium))^2
    !                     / nu * gam((mu+1-2b)/nu) *   n0   / lam^((mu+1-2b)/nu)
    !    =               Zfac                      * mgd%n0 / mgd%lam^lamexp
    ! Note that Zfac includes (location-dependent) K^2 in addition to otherwise
    ! typically (location-)constant parameters.
    Z = Zfac * mgd%n0 / mgd%lam**lamexp

  END FUNCTION zradar_rayleigh_mgd

  !----------------------------------------------------------------------------------------!
  !  Calculate effective radar reflectivity using Rayleigh approximation assuming a        !
  !  modified gamma distribution
  !
  !  FIXME: fix description
  !                                                                                        !
  !  Version works for any hydrometeor of any shape and any moment scheme (these           !
  !  dependencies are imprinted in the input parameters), but requires x =^= L/N be known  !
  !  (hence most practical for 2mom, but also usable for 1mom with known x_mean, like      !
  !  cloud and ice.                                                                        !
  !----------------------------------------------------------------------------------------!

  FUNCTION zradar_rayleigh_L_x(L,x,K2divRho2,gamfac) &
      RESULT (Z)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: L,x,K2divRho2,gamfac
    REAL(KIND=dp)             :: Z

    ! FIXME:
    ! add explanation/derivation of the equation and the meanings of its parameters
    ! K2divRho2 or gamfac need to include the 1d18 factor
    Z = K2divRho2 * gamfac * L * x

  END FUNCTION zradar_rayleigh_L_x

  !----------------------------------------------------------------------------------------!
  !  Calculate effective radar reflectivity using Rayleigh approximation assuming a        !
  !  modified gamma distribution
  !
  !  FIXME: fix description
  !                                                                                        !
  !  Version works for any hydrometeor of any shape and any moment scheme (these           !
  !  dependencies are imprinted in the input parameters).                                  !
  !----------------------------------------------------------------------------------------!

  FUNCTION zradar_rayleigh_L_Lexp(L,Lexp,K2divRho2,gamfac) &
      RESULT (Z)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: L,Lexp,K2divRho2,gamfac
    REAL(KIND=dp)             :: Z

    ! FIXME:
    ! add explanation/derivation of the equation and the meanings of its parameters
    ! K2divRho2 or gamfac need to include the 1d18 factor
    !Z = K2divRho2 * gamfac * L**Lexp

    Z = 0d0
    ! avoid log(0)
    IF (L > quasi_zero) Z = K2divRho2 * gamfac * EXP(LOG(L)*Lexp)

  END FUNCTION zradar_rayleigh_L_Lexp


  !=================================================================================================
  !
  ! Routines for the calculation of the number- and reflectivity aggregated terminal velocity
  !  of hydrometeors, i.e., the integral
  !
  !  \int_0_{\infty} sigma_b(D) vt(D) N(D) dD
  !
  ! respectively
  !
  !  \int_0_{\infty} vt(D) N(D) dD
  !
  !=================================================================================================


  ! All single-phase Vt formulas, both Vt*Ze and Vt*n, for 1- as well as 2 mom
  ! (or unknown and known x) can be reduced to the functional for of
  !
  ! Vt*weight = fac * base**expo
  !
  ! Formulas for all three parameters vary depending on case. When defining
  !
  !     bvel'   = bvelx (= bvelx+2) and
  !     Krhofac = 1     (= K_m**2/K_w0**2 * (6/Pi)**2 / rho_m**2)
  ! for weight  = n     (= Ze), respectively,
  !
  ! (where bvelx, as well as later occuring avelx are - in contrast to all other
  ! parameters - defined in x-notation)
  !
  ! we get for known-x cases (2mom scheme and monodisperse 1mom PSD):
  ! base = x = q/n
  ! expo = bvel'-1
  ! fac  = avelx * q * gamfac(bvel') * Krhofac
  !
  ! and for unknown-x cases (1mom scheme):
  ! base = q
  ! expo = (mu+bvel'*bgeo+1) / (mu+bgeo+1)
  ! fac  = avelx * ageo**bvel' * n0/nu *
  !        [ageo * n0/nu * gam( (mu+bgeo+1)/nu )]**-expo *
  !        gam( (mu+bvel'*bgeo+1)/nu ) * Krhofac
  !
  ! Depending on which microphysical parameters are constant, expo, fac, or parts of
  ! fac can be precalculated and re-used with arbitrary q.
  !
  FUNCTION vtradar_pure_phase(base,expo,fac) &
      RESULT (Vt)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: base, expo, fac

    REAL(KIND=dp) :: Vt

    Vt = 0d0
    IF (base > quasi_zero) THEN
      Vt = fac * EXP( LOG(base) * expo )
    END IF

  END FUNCTION vtradar_pure_phase

  ! All mixed-phase Vt formulas, both Vt*Ze and Vt*n, for 1- as well as 2 mom
  ! (or unknown and known x) can be reduced to the functional for of
  !
  ! Vt*weight = globfac * (    fmfac  * fac_liq * base**expo_liq +
  !                         (1-fmfac) * fac_fro * base**expo_fro )
  !
  ! Formulas for all three parameters vary depending on case. When defining
  !
  !     bvel'   = bvelx (= bvelx+2) and
  !     Krhofac = 1     (= Keffrho**2/K_w0**2 * (6/Pi)**2)
  ! for weight  = n     (= Ze), respectively,
  !
  ! (where bvelx, as well as later occuring avelx are - in contrast to all other
  ! parameters - defined in x-notation)
  !
  ! we get for known-x cases (2mom scheme and monodisperse 1mom PSD):
  ! base     = x = q/n
  ! expo_liq = bvel_liq'-1
  ! expo_fro = bvel_fro'-1
  ! fac_liq  = avelx_liq * gamfac(fro,bvel_liq')
  ! fac_fro  = avelx_fro * gamfac(fro,bvel_fro')
  ! globfac  = q * Krhofac
  !
  ! and for unknown-x cases (1mom scheme):
  ! base = q
  ! expo_liq = (mu+bvel_liq'*bgeo+1) / (mu+bgeo+1)
  ! expo_fro = (mu+bvel_fro'*bgeo+1) / (mu+bgeo+1)
  ! fac_liq  = avelx_liq * ageo**bvel_liq' * n0/nu *
  !            [ageo * n0/nu * gam( (mu+bgeo+1)/nu )]**-expo_liq *
  !            gam( (mu+bvel_liq'*bgeo+1)/nu )
  ! fac_fro  = avelx_fro * ageo**bvel_fro' * n0/nu *
  !            [ageo * n0/nu * gam( (mu+bgeo+1)/nu )]**-expo_fro *
  !            gam( (mu+bvel_fro'*bgeo+1)/nu )
  ! globfac  = Krhofac
  !
  ! Depending on which microphysical parameters are constant, expo_liq/fro,
  ! fac_liq/fro, and globfac or parts thereof can be precalculated and re-used
  ! with arbitrary q.
  ! Considering that Keffrho is T-, hence location dependent, it can be
  ! efficient join constant parts of globfac and fac_liq/fro into precalculable
  ! fully constant fac_liq/fro'.
  !
  FUNCTION vtradar_mixed_phase(&
      base,expo_fro,fac_fro,expo_liq,fac_liq,fac_glo,fmd) &
      RESULT (Vt)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: base, expo_fro, fac_fro, expo_liq, fac_liq, &
                                 fac_glo, fmd

    REAL(KIND=dp) :: Vt

    REAL(KIND=dp) :: logbase

    Vt = 0d0
    IF (base > quasi_zero) THEN
      logbase = LOG(base)
      Vt = fac_glo * &
           ( fmd      * fac_liq * EXP( logbase * expo_liq ) + &
            (1d0-fmd) * fac_fro * EXP( logbase * expo_fro ) )
    END IF

  END FUNCTION vtradar_mixed_phase


  ! For 1mom- cases, where generally mu, nu, and n0 of the MGD are fixed/determined,
  ! with all parameters in D-notation:
  !
  ! n = [ageo * n0/nu * gam( (mu+bgeo+1)/nu )]**-((mu+1)/(mu+bgeo+1)) *
  !     n0/nu * gam( (mu+1)/nu ) * L**((mu+1)/(mu+bgeo+1))
  !
  ! which can be reduced to a funcional form of
  ! n = fac * base**expo with
  !  base = q
  !  expo = (mu+1) / (mu+bgeo+1)
  !  fac  = [ageo * n0/nu * gam( (mu+bgeo+1)/nu )]**-expo *  n0/nu * gam( (mu+1)/nu )
  !
  ! Depending on which microphysical parameters are constant, expo, fac, or parts of
  ! fac can be precalculated and re-used with arbitrary q.
  !
  FUNCTION nradar(base,expo,fac) &
      RESULT (n)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: base, expo, fac

    REAL(KIND=dp) :: n

    n = 0d0
    IF (base > quasi_zero) THEN
      n = fac * EXP( LOG(base) * expo )
    END IF

  END FUNCTION nradar


  FUNCTION K_rho_fac_oguchi(K_i,K_w,inv_rho_ice,inv_rho_w,fm) &
      RESULT (K2divRho2)

    IMPLICIT NONE

    COMPLEX(kind=dp) :: K_i, K_w
    REAL(KIND=dp) :: inv_rho_ice,inv_rho_w,fm

    REAL(KIND=dp) :: K2divRho2

    ! FIXME:
    ! K_i, K_w, inv_rho_ice, inv_rho_w are the same for all hydrometeors,
    ! only fm is different
    ! Hence, instead of passing the K and inv_rho, we could precalculate (once
    ! per model/radar grid point) and pass their (factor) contribution to each
    ! of the 3 additive components. Or does that get too messy?
    !
    ! FIXME:
    ! Try calculating this way: ( K_i*(1-fm)*inv_rho_ice + K_w*fm*inv_rho_w )**2
    ! If that works, we could precalulate and pass the two K-rho-factors only
    ! (note, however, that they are complex!)
    K2divRho2 = ( (ABS(K_i)*(1d0-fm)*inv_rho_ice)**2 + &
                  2d0*(DBLE(K_i)*DBLE(K_w)+AIMAG(K_i)*AIMAG(K_w))* &
                    fm*(1d0-fm)*inv_rho_ice*inv_rho_w + &
                  (ABS(K_w)*fm*inv_rho_w)**2 )

  END FUNCTION K_rho_fac_oguchi


END MODULE radar_mie_specint

