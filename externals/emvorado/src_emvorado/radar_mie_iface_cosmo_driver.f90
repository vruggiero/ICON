!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"

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

#if (defined TWOMOM_SB_NEW || defined TWOMOM_SB_OLD)
#define TWOMOM_SB
#endif

#ifdef _DACE_
#define __COSMO__
#endif

MODULE radar_mie_iface_cosmo_driver

!------------------------------------------------------------------------------
!
! Description:
!      Drivers for the interface functions for the COSMO/ICON microphysics
!      to the EMVORADO libraries for polarimetric radar moments computations.
!
!      Calculates polarization parameters of precipitation particles using
!      full Tmatrix theory. An option for simpler Mie theorie for reflectivity and
!      attenuation only is available. All methods include scattering by wet and melting particles.
!
!      Applicable to microwave radiation, wavelength > 1 mm
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_kind , ONLY :   &
       dp,          &  !           "--"           double precision
       sp,          &  !           "--"           single precision
       wp              !           "--"           model working precision
                       !      (comes directly from the model via USE assoc. in radar_kind.f90)

  USE radar_data , ONLY :   &
       miss_value,  &
       miss_value_rhv,  &
       ie_fwo,      &  ! number of model grid points in zonal direction on the PE
       je_fwo,      &  ! number of model grid points in meridional direction on the PE
       ke_fwo,      &  ! number of model grid points in vertical direction
       i_fwo_ini,       &  ! Timing flag for the initialization of the forward operator
       i_fwo_barrier,   &  ! Timing flag for barrier waiting in MPI-communications (measure for load imbalance)
       rho_w_model, rho_ice_model, K_w_model, K_ice_model, pi_model, &
       t0_melt_model, itype_gscp_model, &
       itlrad_dyn, itlrad_qx, &
       icomm_cart_fwo, my_cart_id_fwo, num_compute_fwo, lcompute_pe_fwo, &
       cnmlstlen

  USE radar_dbzcalc_params_type, ONLY : t_dbzcalc_params, dbz_namlst_d
  
  USE radar_data_namelist, only: &
       ltestpattern_hydrometeors, &
       itype_mpipar_lookupgen, pe_start_lookupgen, pe_end_lookupgen, &
       llookup_interp_mode_dualpol

#if (defined AUXOUT_OFFLINE && defined __COSMO__)
  USE radar_data , ONLY :   &
       ydate_ini_mod, cmaxlen
  USE radar_data_namelist, only: &
       ydirradarout, lmodfield_output, &
       loutdbz, loutpolstd, loutpolall, lextdbz
#endif

  USE radar_data_mie, ONLY : &
       ldebug_dbz,    &  ! This value is set in calc_dbz_vec()
       ilow_modelgrid, iup_modelgrid, jlow_modelgrid, jup_modelgrid, klow_modelgrid, kup_modelgrid, &
       Tmax_i_modelgrid, Tmax_s_modelgrid, Tmax_g_modelgrid, Tmax_h_modelgrid, &
       Tmin_g_modelgrid, Tmin_h_modelgrid, &
       lgsp_fwo, itype_gscp_fwo, luse_muD_relation_rain_fwo

#ifdef __COSMO__
  USE radar_data_mie, ONLY :  klv850_fwo
#endif


  USE radar_interface, ONLY : get_runtime_timings, abort_run, get_datetime_act, &
       get_loc_domain, get_model_hydrometeors, get_model_variables, &
       get_model_config_for_radar, ndoms_max_model,     &
       set_testpattern_hydrometeors_mg, &
       initialize_tmax_1mom_vec_par, finalize_tmax, &
       t, rho, qv, qc, qi, qr, qs, qg, &
       qh, qnc, qni, qnr, qns, qng, qnh, qgl, qhl, qnc_s, &
       lalloc_qi, lalloc_qs, lalloc_qg, lalloc_qh, &
       get_dbz3dlin_with_model_method_1mom

#ifdef TWOMOM_SB
  USE radar_interface, ONLY : get_dbz3dlin_with_model_method_2mom, &
       initialize_tmax_2mom_vec_par
#endif

  USE radar_utilities, ONLY : &
       get_free_funit

  USE radar_mie_utils, ONLY : hash_radar_lookup, get_hashbase

  USE radar_mie_iface_cosmo_1mom, ONLY : &
       radar_mie_1mom_vec, &
       radar_rayleigh_oguchi_1mom_vec, &
       vtradar_rayleigh_oguchi_1mom_vec

#ifdef TWOMOM_SB
  USE radar_mie_iface_cosmo_2mom, ONLY : &
       radar_mie_2mom_vec, &
       radar_rayleigh_oguchi_2mom_vec, &
       vtradar_rayleigh_oguchi_2mom_vec
#endif

#if (defined AUXOUT_OFFLINE && defined __COSMO__)  
  USE radar_mie_iface_cosmo_utils, ONLY : write_modelfield_output
#endif

!===============================================================================

#ifndef NOMPI
  USE mpi
#endif

!===============================================================================

  IMPLICIT NONE

!===============================================================================

#ifdef NOMPI
  INCLUDE "nompi_mpif.h"
#endif

!==============================================================================

  PUBLIC

!==============================================================================

CONTAINS

!===============================================================================
!===============================================================================

  !==============================================================================
  !==============================================================================
  !
  ! INTERFACES FOR GRID POINT REFLECTIVITY/EXTINCTION CALCULATION FOR COSMO/ICON
  !
  !==================================================================================
  !
  ! Calculation of radar reflectivity.
  !
  ! Calculation is only done when it is necessary, i.e.,
  ! only once during one timestep, even if DBZ, DBZ_850 and DZB_CMAX
  ! are specified as output variables. However, if more than one
  ! output namelist is present, it might be that the parameters
  ! specifying the type of DBZ calculation are different among
  ! those namelists, and in that case, calculation has to be done
  ! more than once during one timestep.
  !
  !==================================================================================

  ! Uses vectorized mie subroutines.

  ! Interface for model grid variables:
  SUBROUTINE calc_dbz_vec_modelgrid (itl_dyn, itl_qx, idom, namlist_in, &      ! mandatory input
       l_use_neigh_tmax_melt, ldebug,       &      ! mandatory input
       ydir_lookup_read,                    &      ! mandatory input
       ydir_lookup_write,                   &      ! mandatory input
       lnostore,                            &      ! optional input
       zh_radar, zh_radar_wp,               &
#ifdef __COSMO__
       zh_radar_850, zh_radar_cmax,         &      ! optional output for COSMO only and for namlist_in%itype_refl=4
       zh_radar_850_wp, zh_radar_cmax_wp,   &
#endif
       ah_radar, ah_radar_wp,               &      ! optional output, only for namlist_in%itype_refl=1 or >4
       zv_radar,                            &      ! optional output, only for namlist_in%itype_refl>4 and loutpol
       rrhv_radar, irhv_radar,              &      !                              "
       kdp_radar, adp_radar,                &      !                              "
       zvh_radar)                                  !                              "                    and loutpolall

    !------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,                 INTENT(IN)   :: itl_dyn, itl_qx, idom
    TYPE(t_dbzcalc_params),  INTENT(IN)   :: namlist_in
    LOGICAL, INTENT(in)                   :: l_use_neigh_tmax_melt, ldebug
    CHARACTER(len=*),        INTENT(IN)   :: ydir_lookup_read, ydir_lookup_write
    LOGICAL, INTENT(in), OPTIONAL         :: lnostore

    REAL(kind=dp), OPTIONAL, INTENT(OUT)  :: &
         zh_radar(:,:,:), ah_radar(:,:,:),   &
         zv_radar(:,:,:), rrhv_radar(:,:,:), irhv_radar(:,:,:), &
         kdp_radar(:,:,:), adp_radar(:,:,:), &
         zvh_radar(:,:,:)
    REAL(kind=wp), OPTIONAL, INTENT(OUT)  :: &
         zh_radar_wp(:,:,:), ah_radar_wp(:,:,:)

#ifdef __COSMO__
    REAL(kind=dp), OPTIONAL, INTENT(OUT)  :: &
         zh_radar_850(:,:), zh_radar_cmax(:,:)
    REAL(kind=wp), OPTIONAL, INTENT(OUT)  :: &
         zh_radar_850_wp(:,:), zh_radar_cmax_wp(:,:)
#endif

    ! Helper strings for comparison of actual struct namlist with
    ! that of the last call to calc_dbz_vec_modelgrid, namlist_dbzcalc_vec_store:
    CHARACTER(len=cnmlstlen)              :: cnamnew, cnamold

    ! How to diagnose the N0 of snow PSD in 1-mom cases:
    INTEGER,                 PARAMETER    :: isnow_n0temp = 2

    ! Default value for neigh_tmax_melt (computation neighbourhood for the Tmax
    !  in the degree-of-melting parameterization of radar_mie_meltdegree.f90)
#ifdef __ICON__
    REAL(kind=dp), PARAMETER              :: neigh_tmax_melt_d = miss_value ! m
#else
    REAL(kind=dp), PARAMETER              :: neigh_tmax_melt_d = 5000.0_dp ! m
#endif

    TYPE(t_dbzcalc_params)                :: namlist
    CHARACTER(len=20)                     :: mode

    ! Flags and variables to determine if radar reflectivity calculation is necessary in this
    ! timestep or if it has been already calculated and we can use a SAVE storage array:
    CHARACTER(len=14)                     :: datetime
    LOGICAL                               :: zlnew_timestep
    TYPE(t_dbzcalc_params)                :: namlist_tmp
    CHARACTER (LEN=300)                   :: errstring, dateiname
    LOGICAL                               :: zlcalc_dbz_vec, zl_same_on_all_pes(1), zlcheck(num_compute_fwo)
    INTEGER                               :: mpierr, funit, k
    LOGICAL                               :: fileexist
    REAL(kind=dp)                         :: neigh_tmax_melt
    ! JCS: New logical luse_tmatrix to determine whether to use T-matrix (instead
    ! of existing Mie calculations)
    LOGICAL                               :: luse_tmatrix, ldo_nonsphere, lnostore_loc

    TYPE dbzstore_t
      INTEGER                     :: ie  = -1
      INTEGER                     :: je  = -1
      INTEGER                     :: ke  = -1
      CHARACTER(len=14)           :: datetime = '99999999999999'
      INTEGER                     :: itl_dyn = -HUGE(1)
      INTEGER                     :: itl_qx  = -HUGE(1)
      CHARACTER(len=cnmlstlen)    :: cnamelist
      REAL (KIND=dp), ALLOCATABLE :: dbzvecstore(:,:,:)
      REAL (KIND=dp), ALLOCATABLE :: dbzextstore(:,:,:)
      REAL (KIND=dp), ALLOCATABLE :: zvvecstore(:,:,:)
      REAL (KIND=dp), ALLOCATABLE :: rrhvvecstore(:,:,:)
      REAL (KIND=dp), ALLOCATABLE :: irhvvecstore(:,:,:)
      REAL (KIND=dp), ALLOCATABLE :: kdpvecstore(:,:,:)
      REAL (KIND=dp), ALLOCATABLE :: adpvecstore(:,:,:)
      REAL (KIND=dp), ALLOCATABLE :: zvhvecstore(:,:,:)
    END TYPE dbzstore_t
    TYPE(dbzstore_t), SAVE                :: dbzstore(ndoms_max_model)

    !JM190902: This has been introduced by UB in singlepol in 1908.
    !          Do we need this for other than zh?
    REAL(kind=wp), ALLOCATABLE    :: zh_radar_wp_tmp(:,:,:)  ! only alloc. for  namlist_in%itype_refl=4

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'emvorado::calc_dbz_vec_modelgrid()'

    datetime = get_datetime_act ()

    ! .. Set global module switch ldebug_dbz:
    IF (ldebug) THEN
      WRITE (*,*)  TRIM(yzroutine)//' on proc ', my_cart_id_fwo
      ldebug_dbz = .TRUE.
    ELSE
      ldebug_dbz = .FALSE.
    END IF

    IF (PRESENT(lnostore)) THEN
      lnostore_loc = lnostore
    ELSE
      lnostore_loc = .FALSE.
    END IF

    CALL get_model_config_for_radar ( idom )  ! sets ie, je, ke, itype_gscp_model, rho_w_model, ...

    !.. Initializations for the beginning and each new timestep:
    IF (.NOT. ALLOCATED(dbzstore(idom)%dbzvecstore) .OR. dbzstore(idom)%datetime /= datetime &
           .OR. dbzstore(idom)%itl_dyn /= itl_dyn .OR. dbzstore(idom)%itl_qx /= itl_qx) THEN

      zlnew_timestep = .TRUE.

      IF (.NOT.ALLOCATED(dbzstore(idom)%dbzvecstore)) THEN
        dbzstore(idom)%ie = ie_fwo
        dbzstore(idom)%je = je_fwo
        dbzstore(idom)%ke = ke_fwo
        ALLOCATE(dbzstore(idom)%dbzvecstore(ie_fwo,je_fwo,ke_fwo))
        dbzstore(idom)%dbzvecstore = miss_value
        ALLOCATE(dbzstore(idom)%dbzextstore(ie_fwo,je_fwo,ke_fwo))
        dbzstore(idom)%dbzextstore = 0.0_dp

        ALLOCATE(dbzstore(idom)%zvvecstore(ie_fwo,je_fwo,ke_fwo))
        dbzstore(idom)%zvvecstore = miss_value
        ALLOCATE(dbzstore(idom)%rrhvvecstore(ie_fwo,je_fwo,ke_fwo))
        dbzstore(idom)%rrhvvecstore = miss_value_rhv
        ALLOCATE(dbzstore(idom)%irhvvecstore(ie_fwo,je_fwo,ke_fwo))
        dbzstore(idom)%irhvvecstore = miss_value_rhv
        ALLOCATE(dbzstore(idom)%kdpvecstore(ie_fwo,je_fwo,ke_fwo))
        dbzstore(idom)%kdpvecstore = 0.0_dp
        ALLOCATE(dbzstore(idom)%adpvecstore(ie_fwo,je_fwo,ke_fwo))
        dbzstore(idom)%adpvecstore = 0.0_dp
        ALLOCATE(dbzstore(idom)%zvhvecstore(ie_fwo,je_fwo,ke_fwo))
        dbzstore(idom)%zvhvecstore = miss_value
      END IF

      namlist_tmp = dbz_namlst_d
      ! change struct namlist_dbzcalc_vec_store a little bit,
      ! that it is not equal to that of the last timestep:
      namlist_tmp%ctype_wetsnow_mie = 'gagagagagaga'
      dbzstore(idom)%cnamelist(:) = ' '
      WRITE(dbzstore(idom)%cnamelist,*) namlist_tmp,  .FALSE.

      dbzstore(idom)%datetime = datetime
      dbzstore(idom)%itl_dyn  = itl_dyn
      dbzstore(idom)%itl_qx   = itl_qx

    ELSE

      zlnew_timestep = .FALSE.

    END IF

    IF ( dbzstore(idom)%ie /= ie_fwo .OR.  &
         dbzstore(idom)%je /= je_fwo .OR.  &
         dbzstore(idom)%ke /= ke_fwo ) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,i2,a,i6,a)') 'Not the same size for domain ',idom,' compared to previous call on PE ',my_cart_id_fwo,'!'
      CALL abort_run (my_cart_id_fwo, 20569, &
           'ERROR in '//TRIM(yzroutine)//': '//TRIM(ADJUSTL(errstring)), &
           'src_radar.f90, '//TRIM(yzroutine))
    END IF

    ! If the optional flag lnostore is set TRUE, then there should be no
    !  re-use of the reflectivity and extinction from the storage arrays on the next call to calc_dbz_vec_modelgrid().
    !  To ensure this, we simply set dbzstore(idom)%itl to a non-sense value, so
    !  that zlnew_timestep will be set TRUE on the next call and
    !  reflectivity calculation is triggered rather than using the storage arrays.
    IF (lnostore_loc) THEN
      dbzstore(idom)%itl_dyn = -HUGE(1)
      dbzstore(idom)%itl_qx  = -HUGE(1)
    END IF

    ! If a horizontal neighbourhood should be incorporated to
    !  search for the maximum temperature Tmax (= minimum height)
    !  at which melting hydrometeors exist in a vertical column,
    !  set it to the value of the parameter neigh_tmax_melt_d,
    !  otherwise set it to a value < 0.0 to disable the
    !  horizontal neighbourhood when searching for Tmax in each grid column.
    IF (l_use_neigh_tmax_melt) THEN
      neigh_tmax_melt = neigh_tmax_melt_d
    ELSE
      neigh_tmax_melt = miss_value
    END IF

    ! Write namelists to the helper strings, disregard station_id:
    cnamnew(:) = ' '
    namlist = namlist_in
    ! This is to avoid that different station IDs lead to a recalculation of dBZ.
    !  although the rest of the configuration is the same:
    namlist%station_id = dbz_namlst_d%station_id
    WRITE(cnamnew,*) namlist, l_use_neigh_tmax_melt
    cnamold = dbzstore(idom)%cnamelist

    ! store namlist of the actual call for the next call:
    dbzstore(idom)%cnamelist(:) = ' '
    WRITE(dbzstore(idom)%cnamelist,*) namlist,  l_use_neigh_tmax_melt

    IF ( .NOT. lgsp_fwo ) THEN

      IF (my_cart_id_fwo == 0 .AND. ldebug) WRITE (*,*) '  DBZ is not calculated because of lgsp_fwo = .false.'

    ELSE

      ! .. Criterion for re-calculation of DBZ or copying the values from storage array:
      zlcalc_dbz_vec = ( zlnew_timestep .OR. TRIM(cnamnew) /= TRIM(cnamold) )

      IF (ldebug) THEN

        ! .. Debug output of dbzparams-namelists (one file for each PE and domain):
        dateiname(:) = ' '
        WRITE(dateiname,'("calc_dbz_vec_err_pe-",i3.3,"_dom-",i2.2,".txt")') my_cart_id_fwo, idom

        CALL get_free_funit(funit)

        INQUIRE(file=TRIM(dateiname), exist=fileexist)
        IF (fileexist) THEN
          OPEN(funit , FILE=TRIM(dateiname), FORM='FORMATTED', STATUS='OLD', &
                       action='WRITE', position='APPEND')
        ELSE
          OPEN(funit , FILE=TRIM(dateiname), FORM='FORMATTED', STATUS='REPLACE', &
                       action='WRITE')
          fileexist = .TRUE.
        END IF

        WRITE (funit,'(a)') 'New: '//TRIM(cnamnew)
        WRITE (funit,'(a)') 'Old: '//TRIM(cnamold)

        CLOSE(funit)

        IF (num_compute_fwo > 1 .AND. lcompute_pe_fwo .AND. neigh_tmax_melt > 0.0_dp) THEN

          ! Check if we have the same settings of zlcalc_dbz_vec
          ! on all PROCESSORS. IF not, the PROGRAM might
          ! hang up, since there is an MPI-exchange of temperature in the reflectivity calculations
          ! for namlist%itype_refl /= 4 below:

          zl_same_on_all_pes(1) = zlcalc_dbz_vec

          CALL MPI_GATHER (zl_same_on_all_pes, 1, MPI_LOGICAL,      &
                           zlcheck,            1, MPI_LOGICAL,      &
                           0, icomm_cart_fwo, mpierr)

          IF (my_cart_id_fwo == 0) THEN

            IF ( .NOT.ALL( zlcheck(2:num_compute_fwo) .EQV. zlcheck(1) ) ) THEN

              errstring(:) = ' '
              WRITE (errstring,'(a,i2,a)') 'Not the same settings on all PEs for DBZ re-calculation on domain ',idom,'!'
              CALL abort_run (my_cart_id_fwo, 20573, &
                   'ERROR in '//TRIM(yzroutine)//': '//TRIM(ADJUSTL(errstring)), &
                   'src_radar.f90, '//TRIM(yzroutine))

            END IF

          END IF

        END IF

      END IF  ! ldebug

      IF ( zlcalc_dbz_vec ) THEN

        CALL get_model_variables (itl_dyn, itl_qx, idom)
        CALL get_model_hydrometeors (itl_qx, idom)

        IF (ltestpattern_hydrometeors) THEN
          CALL set_testpattern_hydrometeors_mg   ! mg = modelgrid
        END IF

        ! .. Set ilow_modelgrid, iup_modelgrid, jlow_modelgrid, jup_modelgrid, klow,_modelgrid kup_modelgrid
        !     for dbz-computations in this sub-domain including the halo PE boundaries:
        CALL get_loc_domain( idom, zlinc_boundlines = .TRUE. )

        IF (itype_gscp_fwo < 200) THEN
          ! Set the Tmax_XX_modelgrid parameter for the parameterization of the degree of melting
          !  of each hydrometeor category as 2D fields on the model grid:
          CALL initialize_tmax_1mom_vec_par(neigh_tmax_melt, namlist &
#if (defined AUXOUT_OFFLINE && defined __COSMO__)
               ,.true.&
#endif
                       )
        ELSE
#ifdef TWOMOM_SB
          CALL initialize_tmax_2mom_vec_par(neigh_tmax_melt, namlist &
#if (defined AUXOUT_OFFLINE && defined __COSMO__)
               ,.true.&
#endif
                       )
#endif
        END IF

        SELECT CASE (namlist%itype_refl)
        CASE(1,5,6)
          
          SELECT CASE (namlist%itype_refl)
          CASE (1)
            IF (my_cart_id_fwo == 0 .AND. ldebug) THEN
              WRITE (*,*) '  Calculating DBZ with Mie-Scattering functions, vectorized'
              IF (itype_gscp_fwo < 200 .AND. namlist%llookup_mie) THEN
                WRITE (*,*) '  *** Using lookup tables for Mie-Scattering functions *** '
              END IF
            END IF
            luse_tmatrix = .FALSE.
            ldo_nonsphere = .FALSE.
          CASE (5,6)
            luse_tmatrix = .TRUE.
            IF (namlist%itype_refl == 5) THEN
              ldo_nonsphere = .TRUE.
              mode(:) = ' '
            ELSE
              ldo_nonsphere = .FALSE.
              mode = ' in quasi-Mie mode'
            END IF
            IF (my_cart_id_fwo == 0 .AND. ldebug) THEN
              WRITE (*,*) '  Calculating polarimetry using T-matrix code'//TRIM(ADJUSTL(mode))//', vectorized'
              IF (itype_gscp_fwo < 200 .AND. namlist%llookup_mie) THEN
                WRITE (*,*) '  *** Using lookup tables for T-matrix code *** '
              END IF
            END IF
          END SELECT

          IF (itype_gscp_fwo < 200) THEN
            CALL radar_mie_1mom_vec( my_cart_id_fwo, namlist%lambda_radar, &
                 itype_gscp_fwo, namlist%itype_refl, luse_tmatrix, ldo_nonsphere, &
                 isnow_n0temp, namlist%igraupel_type, namlist%itype_Dref_fmelt, &
                 namlist%ctype_dryice_mie, namlist%ctype_wetice_mie, &
                 namlist%ctype_drysnow_mie, namlist%ctype_wetsnow_mie, &
                 namlist%ctype_drygraupel_mie, namlist%ctype_wetgraupel_mie, &
                 namlist%ldynamic_wetgrowth_gh, &
                 namlist%Tmeltbegin_i, namlist%meltdegTmin_i, namlist%Tmax_min_i, namlist%Tmax_max_i, &
                 namlist%Tmeltbegin_s, namlist%meltdegTmin_s, namlist%Tmax_min_s, namlist%Tmax_max_s, &
                 namlist%Tmeltbegin_g, namlist%meltdegTmin_g, namlist%Tmax_min_g, namlist%Tmax_max_g, &
                 namlist%polMP_r, namlist%polMP_i, namlist%polMP_s, namlist%polMP_g, &
                 rho=rho(:,:,:), t=t(:,:,:), &
                 qc=qc(:,:,:), qr=qr(:,:,:), &
                 qi=qi(:,:,:), qs=qs(:,:,:), &
                 qg=qg(:,:,:), &
                 Tmax_i=Tmax_i_modelgrid(:,:), Tmax_s=Tmax_s_modelgrid(:,:), &
                 Tmax_g=Tmax_g_modelgrid(:,:), &
                 Tmin_g=Tmin_g_modelgrid(:,:), &
                 ilow=ilow_modelgrid, iup=iup_modelgrid, &
                 jlow=jlow_modelgrid, jup=jup_modelgrid, &
                 klow=klow_modelgrid, kup=kup_modelgrid, &
                 lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, &
                 llookup=namlist%llookup_mie, &
                 impipar_lookupgen=itype_mpipar_lookupgen, &
                 pe_start=pe_start_lookupgen, pe_end=pe_end_lookupgen,    &
                 linterp_mode_dualpol=llookup_interp_mode_dualpol,          &
                 ydir_lookup_read=TRIM(ydir_lookup_read), ydir_lookup_write=TRIM(ydir_lookup_write), &
                 ext_tune_fac_pure = namlist%ext_tune_fac_pure, &
                 ext_tune_fac_melt = namlist%ext_tune_fac_melt, &
                 zh_radar = dbzstore(idom)%dbzvecstore(:,:,:), &
                 ah_radar = dbzstore(idom)%dbzextstore(:,:,:), &
                 zv_radar = dbzstore(idom)%zvvecstore(:,:,:), &
                 rrhv_radar = dbzstore(idom)%rrhvvecstore(:,:,:), &
                 irhv_radar = dbzstore(idom)%irhvvecstore(:,:,:), &
                 kdp_radar = dbzstore(idom)%kdpvecstore(:,:,:), &
                 adp_radar = dbzstore(idom)%adpvecstore(:,:,:), &
                 zvh_radar = dbzstore(idom)%zvhvecstore(:,:,:), &
                 lhydrom_choice_testing=namlist%lhydrom_choice_testing )
#ifdef TWOMOM_SB
          ELSE
            CALL radar_mie_2mom_vec( my_cart_id_fwo, namlist%lambda_radar, &
                 itype_gscp_fwo, namlist%itype_refl, luse_tmatrix, ldo_nonsphere, &
                 namlist%igraupel_type, namlist%itype_Dref_fmelt, &
                 namlist%ctype_dryice_mie, namlist%ctype_wetice_mie, &
                 namlist%ctype_drysnow_mie, namlist%ctype_wetsnow_mie, &
                 namlist%ctype_drygraupel_mie, namlist%ctype_wetgraupel_mie, &
                 namlist%ctype_dryhail_mie, namlist%ctype_wethail_mie, &
                 namlist%ldynamic_wetgrowth_gh, &
                 namlist%Tmeltbegin_i, namlist%meltdegTmin_i, namlist%Tmax_min_i, namlist%Tmax_max_i, &
                 namlist%Tmeltbegin_s, namlist%meltdegTmin_s, namlist%Tmax_min_s, namlist%Tmax_max_s, &
                 namlist%Tmeltbegin_g, namlist%meltdegTmin_g, namlist%Tmax_min_g, namlist%Tmax_max_g, &
                 namlist%Tmeltbegin_h, namlist%meltdegTmin_h, namlist%Tmax_min_h, namlist%Tmax_max_h, &
                 namlist%polMP_r, namlist%polMP_i, namlist%polMP_s, namlist%polMP_g, namlist%polMP_h, &
                 rho=rho(:,:,:), t=t(:,:,:), &
                 qc=qc(:,:,:), qr=qr(:,:,:), &
                 qi=qi(:,:,:), qs=qs(:,:,:), &
                 qg=qg(:,:,:), qh=qh(:,:,:), &
                 qnc=qnc(:,:,:), qnr=qnr(:,:,:), &
                 qni=qni(:,:,:), qns=qns(:,:,:), &
                 qng=qng(:,:,:), qnh=qnh(:,:,:), &
                 qgl=qgl(:,:,:), qhl=qhl(:,:,:), &
                 Tmax_i=Tmax_i_modelgrid(:,:), Tmax_s=Tmax_s_modelgrid(:,:), &
                 Tmax_g=Tmax_g_modelgrid(:,:), Tmax_h=Tmax_h_modelgrid(:,:), &
                 Tmin_g=Tmin_g_modelgrid(:,:), Tmin_h=Tmin_h_modelgrid(:,:), &
                 ilow=ilow_modelgrid, iup=iup_modelgrid, &
                 jlow=jlow_modelgrid, jup=jup_modelgrid, &
                 klow=klow_modelgrid, kup=kup_modelgrid, &
                 lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, &
                 lalloc_qh=lalloc_qh, &
                 llookup=namlist%llookup_mie, &
                 impipar_lookupgen=itype_mpipar_lookupgen, &
                 pe_start=pe_start_lookupgen, pe_end=pe_end_lookupgen,    &
                 linterp_mode_dualpol=llookup_interp_mode_dualpol,          &
                 luse_muD_relation_rain = luse_muD_relation_rain_fwo, &
                 ydir_lookup_read=TRIM(ydir_lookup_read), ydir_lookup_write=TRIM(ydir_lookup_write), &
                 ext_tune_fac_pure = namlist%ext_tune_fac_pure, &
                 ext_tune_fac_melt = namlist%ext_tune_fac_melt, &
                 zh_radar = dbzstore(idom)%dbzvecstore(:,:,:), &
                 ah_radar = dbzstore(idom)%dbzextstore(:,:,:), &
                 zv_radar = dbzstore(idom)%zvvecstore(:,:,:), &
                 rrhv_radar = dbzstore(idom)%rrhvvecstore(:,:,:), &
                 irhv_radar = dbzstore(idom)%irhvvecstore(:,:,:), &
                 kdp_radar = dbzstore(idom)%kdpvecstore(:,:,:), &
                 adp_radar = dbzstore(idom)%adpvecstore(:,:,:), &
                 zvh_radar = dbzstore(idom)%zvhvecstore(:,:,:), &
                 lhydrom_choice_testing=namlist%lhydrom_choice_testing )
#endif
          END IF
        CASE(3)
          IF (my_cart_id_fwo == 0 .AND. ldebug) WRITE (*,*) '  Calculating DBZ ', &
               'with Rayleigh approximation and Oguchi refractive index, vectorized'
          IF (itype_gscp_fwo < 200) THEN
            CALL radar_rayleigh_oguchi_1mom_vec( my_cart_id_fwo, namlist%lambda_radar, &
                 itype_gscp_fwo, isnow_n0temp, &
                 namlist%Tmeltbegin_i, namlist%meltdegTmin_i, &
                 namlist%Tmeltbegin_s, namlist%meltdegTmin_s, &
                 namlist%Tmeltbegin_g, namlist%meltdegTmin_g, &
                 rho=rho(:,:,:), t=t(:,:,:), &
                 qc=qc(:,:,:), qr=qr(:,:,:), qi=qi(:,:,:), qs=qs(:,:,:), qg=qg(:,:,:), &
                 Tmax_i=Tmax_i_modelgrid(:,:), Tmax_s=Tmax_s_modelgrid(:,:), &
                 Tmax_g=Tmax_g_modelgrid(:,:), &
                 Tmin_g=Tmin_g_modelgrid(:,:), &
                 ilow=ilow_modelgrid, iup=iup_modelgrid, jlow=jlow_modelgrid, jup=jup_modelgrid, &
                 klow=klow_modelgrid, kup=kup_modelgrid, &
                 lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, &
                 zh_radar = dbzstore(idom)%dbzvecstore(:,:,:), &
                 lhydrom_choice_testing=namlist%lhydrom_choice_testing )
#ifdef TWOMOM_SB
          ELSE
            CALL radar_rayleigh_oguchi_2mom_vec( my_cart_id_fwo, namlist%lambda_radar, &
                 namlist%Tmeltbegin_i, namlist%meltdegTmin_i, &
                 namlist%Tmeltbegin_s, namlist%meltdegTmin_s, &
                 namlist%Tmeltbegin_g, namlist%meltdegTmin_g, &
                 namlist%Tmeltbegin_h, namlist%meltdegTmin_h, &
                 rho=rho(:,:,:), t=t(:,:,:), &
                 qc=qc(:,:,:), qr=qr(:,:,:), &
                 qi=qi(:,:,:), qs=qs(:,:,:), &
                 qg=qg(:,:,:), qh=qh(:,:,:), &
                 qnc=qnc(:,:,:), qnr=qnr(:,:,:), &
                 qni=qni(:,:,:), qns=qns(:,:,:), &
                 qng=qng(:,:,:), qnh=qnh(:,:,:), &
                 qgl=qgl(:,:,:), qhl=qhl(:,:,:), &
                 Tmax_i=Tmax_i_modelgrid(:,:), Tmax_s=Tmax_s_modelgrid(:,:), &
                 Tmax_g=Tmax_g_modelgrid(:,:), Tmax_h=Tmax_h_modelgrid(:,:), &
                 Tmin_g=Tmin_g_modelgrid(:,:), Tmin_h=Tmin_h_modelgrid(:,:), &
                 ilow=ilow_modelgrid, iup=iup_modelgrid, &
                 jlow=jlow_modelgrid, jup=jup_modelgrid, &
                 klow=klow_modelgrid, kup=kup_modelgrid, &
                 lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, lalloc_qh=lalloc_qh, &
                 luse_muD_relation_rain = luse_muD_relation_rain_fwo, &                 
                 zh_radar = dbzstore(idom)%dbzvecstore(:,:,:), &
                 lhydrom_choice_testing=namlist%lhydrom_choice_testing )
#endif
          END IF
          dbzstore(idom)%dbzextstore(:,:,:) = 0.0_dp
        CASE(4)

          ALLOCATE (zh_radar_wp_tmp(ie_fwo, je_fwo, ke_fwo))

          IF (my_cart_id_fwo == 0 .AND. ldebug) WRITE (*,*) '  Calculating DBZ ', &
               'with Rayleigh approximation from the hosting model'

          IF (itype_gscp_fwo == 140) THEN
            CALL get_dbz3dlin_with_model_method_1mom (ni=ie_fwo, nj=je_fwo, nk=ke_fwo, &
                 itype_gscp_model_in=itype_gscp_model, &
                 rho_water=rho_w_model, rho_ice=rho_ice_model, &
                 K_water=K_w_model, K_ice=K_ice_model, t0melt=t0_melt_model, &
                 t_in=t(:,:,:), rho_in=rho(:,:,:), &
                 qc_in=qc(:,:,:), qr_in=qr(:,:,:), &
                 qi_in=qi(:,:,:), qs_in=qs(:,:,:), &
                 qnc_s=qnc_s(:,:), &
                 z_radar=zh_radar_wp_tmp)
          ELSEIF (itype_gscp_fwo == 150) THEN
            CALL get_dbz3dlin_with_model_method_1mom (ni=ie_fwo, nj=je_fwo, nk=ke_fwo, &
                 itype_gscp_model_in=itype_gscp_model, &
                 rho_water=rho_w_model, rho_ice=rho_ice_model, &
                 K_water=K_w_model, K_ice=K_ice_model, t0melt=t0_melt_model, &
                 t_in=t(:,:,:), rho_in=rho(:,:,:), &
                 qc_in=qc(:,:,:), qr_in=qr(:,:,:), &
                 qi_in=qi(:,:,:), qs_in=qs(:,:,:), &
                 qg_in=qg(:,:,:), &
                 qnc_s=qnc_s(:,:), &
                 z_radar=zh_radar_wp_tmp)
#ifdef TWOMOM_SB
          ELSEIF (itype_gscp_fwo >= 200) THEN
            CALL get_dbz3dlin_with_model_method_2mom (ni=ie_fwo, nj=je_fwo, nk=ke_fwo, &
                 itype_gscp_model_in=itype_gscp_model, &
                 rho_water=rho_w_model, rho_ice=rho_ice_model, &
                 K_water=K_w_model, K_ice=K_ice_model, t0melt=t0_melt_model, &
                 luse_muD_relation_rain = luse_muD_relation_rain_fwo, &
                 t_in=t(:,:,:), rho_in=rho(:,:,:), &
                 qc_in=qc(:,:,:), qr_in=qr(:,:,:), &
                 qi_in=qi(:,:,:), qs_in=qs(:,:,:), &
                 qg_in=qg(:,:,:), qh_in=qh(:,:,:), &
                 qnc_in=qnc(:,:,:), qnr_in=qnr(:,:,:), &
                 qni_in=qni(:,:,:), qns_in=qns(:,:,:), &
                 qng_in=qng(:,:,:), qnh_in=qnh(:,:,:), &
                 qgl_in=qgl(:,:,:), qhl_in=qhl(:,:,:), &
                 z_radar=zh_radar_wp_tmp)
#endif
          ELSE
            IF (my_cart_id_fwo == 0) THEN
              errstring(:) = ' '
              WRITE (errstring,'(a,i1,a,i4)') &
                'Error itype_gscp_model: For itype_refl=',namlist%itype_refl, &
                ' gscp scheme not implemented in EMVORADO: itype_gscp_model = ', &
                itype_gscp_model
                CALL abort_run (my_cart_id_fwo, 35620, errstring, yzroutine)
            END IF
          END IF

          IF (wp == sp) THEN
            dbzstore(idom)%dbzvecstore(:,:,:) = REAL(zh_radar_wp_tmp, kind=dp)
          ELSE
            dbzstore(idom)%dbzvecstore(:,:,:) = zh_radar_wp_tmp
          END IF
          dbzstore(idom)%dbzextstore(:,:,:) = 0.0_dp
          DEALLOCATE (zh_radar_wp_tmp)

        CASE default
          IF (my_cart_id_fwo == 0) THEN
            WRITE (*,*) '  SRC_OUTPUT: Unknown type of radar reflectivity calculation: ',namlist%itype_refl
          END IF
          dbzstore(idom)%dbzvecstore(:,:,:) = miss_value
          dbzstore(idom)%dbzextstore(:,:,:) = 0.0_dp
          dbzstore(idom)%zvvecstore = miss_value
          dbzstore(idom)%rrhvvecstore = miss_value_rhv
          dbzstore(idom)%irhvvecstore = miss_value_rhv
          dbzstore(idom)%kdpvecstore = 0.0_dp
          dbzstore(idom)%adpvecstore = 0.0_dp
          dbzstore(idom)%zvhvecstore = miss_value
        END SELECT

      ELSE    ! ... IF (zlcalc_dbz_vec)

        IF (my_cart_id_fwo == 0 .AND. ldebug) THEN

          SELECT CASE (namlist%itype_refl)
          CASE(1)
            WRITE (*,*) '  Taking previous DBZ with Mie-Scattering functions'
          CASE(2)
            WRITE (*,*) '  Taking previous DBZ with Rayleigh approximation functions'
          CASE(3)
            WRITE (*,*) '  Taking previous DBZ with Rayleigh/Oguchi approximation'
          CASE(4)
            WRITE (*,*) '  Taking previous DBZ with Rayleigh approximation from hosting model'
          CASE(5)
            WRITE (*,*) '  Taking previous DBZ with T-matrix method'
          CASE(6)
            WRITE (*,*) '  Taking previous DBZ with T-matrix in quasi-Mie mode method'
          CASE default
            WRITE (*,*) '  SRC_OUTPUT: Unknown type of radar reflectivity calculation: ',namlist%itype_refl
          END SELECT

        END IF

      END IF  ! ... IF (zlcalc_dbz_vec)

    END IF    ! ... IF ( .NOT.lgsp_fwo )

!$OMP PARALLEL PRIVATE(k)
    IF (PRESENT(zh_radar)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        zh_radar(:,:,k) = dbzstore(idom)%dbzvecstore(:,:,k)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(ah_radar)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        ah_radar(:,:,k) = dbzstore(idom)%dbzextstore(:,:,k)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(zh_radar_wp)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        zh_radar_wp(:,:,k) = REAL(dbzstore(idom)%dbzvecstore(:,:,k), kind=wp)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(ah_radar_wp)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        ah_radar_wp(:,:,k) = REAL(dbzstore(idom)%dbzextstore(:,:,k), kind=wp)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(zv_radar)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        zv_radar(:,:,k) = dbzstore(idom)%zvvecstore(:,:,k)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(rrhv_radar)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        rrhv_radar(:,:,k) = dbzstore(idom)%rrhvvecstore(:,:,k)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(irhv_radar)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        irhv_radar(:,:,k) = dbzstore(idom)%irhvvecstore(:,:,k)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(kdp_radar)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        kdp_radar(:,:,k) = dbzstore(idom)%kdpvecstore(:,:,k)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(adp_radar)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        adp_radar(:,:,k) = dbzstore(idom)%adpvecstore(:,:,k)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(zvh_radar)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        zvh_radar(:,:,k) = dbzstore(idom)%zvhvecstore(:,:,k)
      ENDDO
!$OMP END DO
    END IF
!$OMP END PARALLEL

#ifdef __COSMO__
    IF ( PRESENT(zh_radar_850) ) THEN
      zh_radar_850(:,:) = dbzstore(idom)%dbzvecstore(:,:,klv850_fwo)
    ENDIF

    IF ( PRESENT(zh_radar_850_wp) ) THEN
      zh_radar_850_wp(:,:) = REAL(dbzstore(idom)%dbzvecstore(:,:,klv850_fwo), kind=wp)
    ENDIF

    IF ( PRESENT(zh_radar_cmax) ) THEN
      zh_radar_cmax(:,:) = MAXVAL(dbzstore(idom)%dbzvecstore,3)
    ENDIF

    IF ( PRESENT(zh_radar_cmax_wp) ) THEN
      zh_radar_cmax_wp(:,:) = REAL(MAXVAL(dbzstore(idom)%dbzvecstore,3), kind=wp)
    ENDIF
#endif

#if (defined AUXOUT_OFFLINE && defined __COSMO__)
    IF (zlcalc_dbz_vec .AND. .NOT.lnostore_loc) THEN
      IF (lmodfield_output) THEN
        CALL write_modelfield_output(namlist_in%station_id, namlist%itype_refl, &
                                     zh_radar, ah_radar, &
                                     zv_radar, rrhv_radar, irhv_radar, &
                                     kdp_radar, adp_radar, zvh_radar)
      END IF
    END IF
#endif

    CALL finalize_tmax ()

    IF (ldebug) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_cart_id_fwo

    RETURN
  END SUBROUTINE calc_dbz_vec_modelgrid


  ! General interface:
  !  (NOTE:  CALL get_model_config_for_radar ( idom ) before this!)

  SUBROUTINE calc_dbz_vec_generic (namlist,          & ! mandatory input
       ldebug,                                       & ! mandatory input
       ydir_lookup_read,                             & ! mandatory input
       ydir_lookup_write,                            & ! mandatory input
       rho, t, qc, qr, qi, qs, qg, qh,               & ! mandatory input
       qnc, qnr, qni, qns, qng, qnh, qgl, qhl, qnc_s,& ! mandatory input
       Tmax_i, Tmax_s, Tmax_g, Tmax_h,               & ! mandatory input
       Tmin_g, Tmin_h,                               & ! mandatory input
       zh_radar, ah_radar,                           & ! optional output
       zv_radar, rrhv_radar, irhv_radar,             & ! optional output
       kdp_radar, adp_radar, zvh_radar)                ! optional output

    !------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(t_dbzcalc_params),  INTENT(IN)     :: namlist
    LOGICAL, INTENT(in)                     :: ldebug
    CHARACTER(len=*),        INTENT(IN)     :: ydir_lookup_read, ydir_lookup_write
    REAL(kind=wp), DIMENSION(:,:,:), INTENT(in) :: &
         rho, t, qc, qr

    REAL(kind=wp), DIMENSION(:,:,:), INTENT(in) :: &
         qnc, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, qgl, qhl

    REAL(kind=wp), DIMENSION(:,:), INTENT(inout) :: &
         Tmax_i, Tmax_s, Tmax_g, Tmax_h, Tmin_g, Tmin_h, qnc_s

    REAL(kind=dp), OPTIONAL, INTENT(OUT) :: &
         zh_radar(:,:,:), ah_radar(:,:,:), &
         zv_radar(:,:,:), rrhv_radar(:,:,:), irhv_radar(:,:,:), &
         kdp_radar(:,:,:), adp_radar(:,:,:), zvh_radar(:,:,:)


    ! How to diagnose the N0 of snow PSD in 1-mom cases:
    INTEGER,                 PARAMETER     :: isnow_n0temp = 2

    CHARACTER(len=20)                     :: mode
    
    ! Flags and variables to determine if radar reflectivity calculation is necessary in this
    ! timestep or if it has been already calculated and we can use a SAVE storage array:
    REAL (KIND=dp), ALLOCATABLE  :: dbzvecstore(:,:,:)
    REAL (KIND=dp), ALLOCATABLE  :: dbzextstore(:,:,:)
    REAL (KIND=dp), ALLOCATABLE  :: zvvecstore(:,:,:)
    REAL (KIND=dp), ALLOCATABLE  :: rrhvvecstore(:,:,:)
    REAL (KIND=dp), ALLOCATABLE  :: irhvvecstore(:,:,:)
    REAL (KIND=dp), ALLOCATABLE  :: kdpvecstore(:,:,:)
    REAL (KIND=dp), ALLOCATABLE  :: adpvecstore(:,:,:)
    REAL (KIND=dp), ALLOCATABLE  :: zvhvecstore(:,:,:)

    INTEGER                         ni, nj, nk
    LOGICAL                      :: lalloc_qi, lalloc_qs, lalloc_qg, lalloc_qh
    ! JCS: New logical luse_tmatrix to determine whether to use T-matrix (instead
    ! of existing Mie calculations)
    LOGICAL                      :: luse_tmatrix,ldo_nonsphere

    REAL(kind=wp), ALLOCATABLE   :: zh_radar_wp_tmp(:,:,:)  ! only alloc. for  namlist_in%itype_refl=4

    CHARACTER (LEN=300)          :: errstring
    CHARACTER (len=*), PARAMETER :: yzroutine = 'emvorado::calc_dbz_vec_generic()'

    ! .. Set global module switch ldebug_dbz:
    IF (ldebug) THEN
      WRITE (*,*)  TRIM(yzroutine)//' on proc ', my_cart_id_fwo
      ldebug_dbz = .TRUE.
    ELSE
      ldebug_dbz = .FALSE.
    END IF

    ni = SIZE(rho, dim=1)
    nj = SIZE(rho, dim=2)
    nk = SIZE(rho, dim=3)

    ALLOCATE(dbzvecstore(ni,nj,nk), dbzextstore(ni,nj,nk))
    ALLOCATE(zvvecstore(ni,nj,nk), &
             rrhvvecstore(ni,nj,nk), irhvvecstore(ni,nj,nk), &
             kdpvecstore(ni,nj,nk), adpvecstore(ni,nj,nk))
    ALLOCATE(zvhvecstore(ni,nj,nk))
    dbzvecstore  = miss_value
    dbzextstore  = 0.0_dp
    zvvecstore   = miss_value
    rrhvvecstore = miss_value_rhv
    irhvvecstore = miss_value_rhv
    kdpvecstore  = 0.0_dp
    adpvecstore  = 0.0_dp
    zvhvecstore  = miss_value

    IF ( .NOT.lgsp_fwo ) THEN

      IF (my_cart_id_fwo == 0 .AND. ldebug) WRITE (*,*) '  DBZ is not calculated because of lgsp_fwo = .false.'

    ELSE

      SELECT CASE (itype_gscp_fwo)
      CASE(120)
        lalloc_qi = .FALSE.
        lalloc_qs = .FALSE.
        lalloc_qg = .FALSE.
      CASE(130)
        lalloc_qi = .FALSE.
        lalloc_qs = .TRUE.
        lalloc_qg = .FALSE.
      CASE(140)
        lalloc_qi = .TRUE.
        lalloc_qs = .TRUE.
        lalloc_qg = .FALSE.
      CASE(150)
        lalloc_qi = .TRUE.
        lalloc_qs = .TRUE.
        lalloc_qg = .TRUE.
      CASE(250)
        lalloc_qh = .FALSE.
      CASE(260)
        lalloc_qh = .TRUE.
      CASE default
        IF (my_cart_id_fwo == 0) THEN
          WRITE (*,*) 'ERROR '//TRIM(yzroutine)//': Unknown itype_gscp_fwo = ', &
                      itype_gscp_fwo, ', not implemented!'
        END IF
      END SELECT

      SELECT CASE (namlist%itype_refl)
      CASE(1,5,6)

        SELECT CASE (namlist%itype_refl)
        CASE (1)
          IF (my_cart_id_fwo == 0 .AND. ldebug) THEN
            WRITE (*,*) '  Calculating DBZ with Mie-Scattering functions, vectorized'
            IF (itype_gscp_fwo < 200 .AND. namlist%llookup_mie) THEN
              WRITE (*,*) '  *** Using lookup tables for Mie-Scattering functions *** '
            END IF
          END IF
          luse_tmatrix = .FALSE.
          ldo_nonsphere = .FALSE.
        CASE (5,6)
          luse_tmatrix = .TRUE.
          IF (namlist%itype_refl == 5) THEN
            ldo_nonsphere = .TRUE.
            mode(:) = ' '
          ELSE
            ldo_nonsphere = .FALSE.
            mode = ' in quasi-Mie mode'
          END IF
          IF (my_cart_id_fwo == 0 .AND. ldebug) THEN
            WRITE (*,*) '  Calculating polarimetry using T-matrix code'//TRIM(ADJUSTL(mode))//', vectorized'
            IF (itype_gscp_fwo < 200 .AND. namlist%llookup_mie) THEN
              WRITE (*,*) '  *** Using lookup tables for T-matrix code *** '
            END IF
          END IF
        END SELECT

        IF (itype_gscp_fwo < 200) THEN
          CALL radar_mie_1mom_vec( my_cart_id_fwo, namlist%lambda_radar, &
               itype_gscp_fwo, namlist%itype_refl, luse_tmatrix, ldo_nonsphere, &
               isnow_n0temp, namlist%igraupel_type, namlist%itype_Dref_fmelt, &
               namlist%ctype_dryice_mie, namlist%ctype_wetice_mie, &
               namlist%ctype_drysnow_mie, namlist%ctype_wetsnow_mie, &
               namlist%ctype_drygraupel_mie, namlist%ctype_wetgraupel_mie, &
               namlist%ldynamic_wetgrowth_gh, &
               namlist%Tmeltbegin_i, namlist%meltdegTmin_i, namlist%Tmax_min_i, namlist%Tmax_max_i, &
               namlist%Tmeltbegin_s, namlist%meltdegTmin_s, namlist%Tmax_min_s, namlist%Tmax_max_s, &
               namlist%Tmeltbegin_g, namlist%meltdegTmin_g, namlist%Tmax_min_g, namlist%Tmax_max_g, &
               namlist%polMP_r, namlist%polMP_i, namlist%polMP_s, namlist%polMP_g, &
               rho=rho(:,:,:), t=t(:,:,:), &
               qc=qc(:,:,:), qr=qr(:,:,:), &
               qi=qi(:,:,:), qs=qs(:,:,:), &
               qg=qg(:,:,:), &
               Tmax_i=Tmax_i(:,:), Tmax_s=Tmax_s(:,:), &
               Tmax_g=Tmax_g(:,:), &
               Tmin_g=Tmin_g(:,:), &
               ilow=1, iup=ni, jlow=1, jup=nj, klow=1, kup=nk, &
               lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, &
               llookup=namlist%llookup_mie, &
               impipar_lookupgen=itype_mpipar_lookupgen, &
               pe_start=pe_start_lookupgen, pe_end=pe_end_lookupgen,    &
               linterp_mode_dualpol=llookup_interp_mode_dualpol,          &
               ydir_lookup_read=TRIM(ydir_lookup_read), ydir_lookup_write=TRIM(ydir_lookup_write), &
               ext_tune_fac_pure = namlist%ext_tune_fac_pure, &
               ext_tune_fac_melt = namlist%ext_tune_fac_melt, &
               zh_radar = dbzvecstore(:,:,:), &
               ah_radar = dbzextstore(:,:,:), &
               zv_radar = zvvecstore(:,:,:), &
               rrhv_radar = rrhvvecstore(:,:,:), &
               irhv_radar = irhvvecstore(:,:,:), &
               kdp_radar = kdpvecstore(:,:,:), &
               adp_radar = adpvecstore(:,:,:), &
               zvh_radar = zvhvecstore(:,:,:), &
               lhydrom_choice_testing=namlist%lhydrom_choice_testing )
#ifdef TWOMOM_SB
        ELSE
          CALL radar_mie_2mom_vec( my_cart_id_fwo, namlist%lambda_radar, &
               itype_gscp_fwo, namlist%itype_refl, luse_tmatrix, ldo_nonsphere, &
               namlist%igraupel_type, namlist%itype_Dref_fmelt, &
               namlist%ctype_dryice_mie, namlist%ctype_wetice_mie, &
               namlist%ctype_drysnow_mie, namlist%ctype_wetsnow_mie, &
               namlist%ctype_drygraupel_mie, namlist%ctype_wetgraupel_mie, &
               namlist%ctype_dryhail_mie, namlist%ctype_wethail_mie, &
               namlist%ldynamic_wetgrowth_gh, &
               namlist%Tmeltbegin_i, namlist%meltdegTmin_i, namlist%Tmax_min_i, namlist%Tmax_max_i, &
               namlist%Tmeltbegin_s, namlist%meltdegTmin_s, namlist%Tmax_min_s, namlist%Tmax_max_s, &
               namlist%Tmeltbegin_g, namlist%meltdegTmin_g, namlist%Tmax_min_g, namlist%Tmax_max_g, &
               namlist%Tmeltbegin_h, namlist%meltdegTmin_h, namlist%Tmax_min_h, namlist%Tmax_max_h, &
               namlist%polMP_r, namlist%polMP_i, namlist%polMP_s, namlist%polMP_g, namlist%polMP_h, &
               rho=rho(:,:,:), t=t(:,:,:), &
               qc=qc(:,:,:), qr=qr(:,:,:), &
               qi=qi(:,:,:), qs=qs(:,:,:), &
               qg=qg(:,:,:), qh=qh(:,:,:), &
               qnc=qnc(:,:,:), qnr=qnr(:,:,:), &
               qni=qni(:,:,:), qns=qns(:,:,:), &
               qng=qng(:,:,:), qnh=qnh(:,:,:), &
               qgl=qgl(:,:,:), qhl=qhl(:,:,:), &
               Tmax_i=Tmax_i(:,:), Tmax_s=Tmax_s(:,:), &
               Tmax_g=Tmax_g(:,:), Tmax_h=Tmax_h(:,:), &
               Tmin_g=Tmin_g(:,:), Tmin_h=Tmin_h(:,:), &
               ilow=1, iup=ni, jlow=1, jup=nj, klow=1, kup=nk, &
               lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, &
               lalloc_qh=lalloc_qh, &
               llookup=namlist%llookup_mie, &
               impipar_lookupgen=itype_mpipar_lookupgen, &
               pe_start=pe_start_lookupgen, pe_end=pe_end_lookupgen,    &
               linterp_mode_dualpol=llookup_interp_mode_dualpol,          &
               luse_muD_relation_rain = luse_muD_relation_rain_fwo, &
               ydir_lookup_read=TRIM(ydir_lookup_read), ydir_lookup_write=TRIM(ydir_lookup_write), &
               ext_tune_fac_pure = namlist%ext_tune_fac_pure, &
               ext_tune_fac_melt = namlist%ext_tune_fac_melt, &
               zh_radar = dbzvecstore(:,:,:), &
               ah_radar = dbzextstore(:,:,:), &
               zv_radar = zvvecstore(:,:,:), &
               rrhv_radar = rrhvvecstore(:,:,:), &
               irhv_radar = irhvvecstore(:,:,:), &
               kdp_radar = kdpvecstore(:,:,:), &
               adp_radar = adpvecstore(:,:,:), &
               zvh_radar = zvhvecstore(:,:,:), &
               lhydrom_choice_testing=namlist%lhydrom_choice_testing )
#endif
        END IF
      CASE(3)
        IF (my_cart_id_fwo == 0 .AND. ldebug) WRITE (*,*) '  Calculating DBZ ', &
             'with Rayleigh approximation and Oguchi refractive index, vectorized'

        IF (itype_gscp_fwo < 200) THEN
          CALL radar_rayleigh_oguchi_1mom_vec( my_cart_id_fwo, namlist%lambda_radar, &
               itype_gscp_fwo, isnow_n0temp, &
               namlist%Tmeltbegin_i, namlist%meltdegTmin_i, &
               namlist%Tmeltbegin_s, namlist%meltdegTmin_s, &
               namlist%Tmeltbegin_g, namlist%meltdegTmin_g, &
               rho=rho(:,:,:), t=t(:,:,:), &
               qc=qc(:,:,:), qr=qr(:,:,:), &
               qi=qi(:,:,:), qs=qs(:,:,:), &
               qg=qg(:,:,:), &
               Tmax_i=Tmax_i(:,:), Tmax_s=Tmax_s(:,:), &
               Tmax_g=Tmax_g(:,:), &
               Tmin_g=Tmin_g(:,:), &
               ilow=1, iup=ni, jlow=1, jup=nj, klow=1, kup=nk, &
               lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, &
               zh_radar = dbzvecstore(:,:,:), &
               lhydrom_choice_testing=namlist%lhydrom_choice_testing )
#ifdef TWOMOM_SB
        ELSE
          CALL radar_rayleigh_oguchi_2mom_vec( my_cart_id_fwo, namlist%lambda_radar, &
               namlist%Tmeltbegin_i, namlist%meltdegTmin_i, &
               namlist%Tmeltbegin_s, namlist%meltdegTmin_s, &
               namlist%Tmeltbegin_g, namlist%meltdegTmin_g, &
               namlist%Tmeltbegin_h, namlist%meltdegTmin_h, &
               rho=rho(:,:,:), t=t(:,:,:), &
               qc=qc(:,:,:), qr=qr(:,:,:), &
               qi=qi(:,:,:), qs=qs(:,:,:), &
               qg=qg(:,:,:), qh=qh(:,:,:), &
               qnc=qnc(:,:,:), qnr=qnr(:,:,:), &
               qni=qni(:,:,:), qns=qns(:,:,:), &
               qng=qng(:,:,:), qnh=qnh(:,:,:), &
               qgl=qgl(:,:,:), qhl=qhl(:,:,:), &
               Tmax_i=Tmax_i(:,:), Tmax_s=Tmax_s(:,:), &
               Tmax_g=Tmax_g(:,:), Tmax_h=Tmax_h(:,:), &
               Tmin_g=Tmin_g(:,:), Tmin_h=Tmin_h(:,:), &
               ilow=1, iup=ni, jlow=1, jup=nj, klow=1, kup=nk, &
               lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, lalloc_qh=lalloc_qh, &
               luse_muD_relation_rain = luse_muD_relation_rain_fwo, &
               zh_radar = dbzvecstore(:,:,:), &
               lhydrom_choice_testing=namlist%lhydrom_choice_testing )
#endif
        END IF
        dbzextstore(:,:,:) = 0.0_dp
      CASE(4)

        ALLOCATE (zh_radar_wp_tmp(ni, nj, nk))

        IF (my_cart_id_fwo == 0 .AND. ldebug) WRITE (*,*) '  Calculating DBZ ', &
             'with Rayleigh approximation from hosting model'

        IF (itype_gscp_fwo == 140) THEN
          CALL get_dbz3dlin_with_model_method_1mom (ni=ni, nj=nj, nk=nk, &
               itype_gscp_model_in=itype_gscp_model, &
               rho_water=rho_w_model, rho_ice=rho_ice_model, &
               K_water=K_w_model, K_ice=K_ice_model, t0melt=t0_melt_model, &
               t_in=t(:,:,:), rho_in=rho(:,:,:), &
               qc_in=qc(:,:,:), qr_in=qr(:,:,:), &
               qi_in=qi(:,:,:), qs_in=qs(:,:,:), &
               qnc_s=qnc_s(:,:), &
               z_radar=zh_radar_wp_tmp)
        ELSEIF (itype_gscp_fwo == 150) THEN
          CALL get_dbz3dlin_with_model_method_1mom (ni=ni, nj=nj, nk=nk, &
               itype_gscp_model_in=itype_gscp_model, &
               rho_water=rho_w_model, rho_ice=rho_ice_model, &
               K_water=K_w_model, K_ice=K_ice_model, t0melt=t0_melt_model, &
               t_in=t(:,:,:), rho_in=rho(:,:,:), &
               qc_in=qc(:,:,:), qr_in=qr(:,:,:), &
               qi_in=qi(:,:,:), qs_in=qs(:,:,:), &
               qg_in=qg(:,:,:), &
               qnc_s=qnc_s(:,:), &
               z_radar=zh_radar_wp_tmp)
#ifdef TWOMOM_SB
        ELSEIF (itype_gscp_fwo >= 200) THEN
          CALL get_dbz3dlin_with_model_method_2mom (ni=ni, nj=nj, nk=nk, &
               itype_gscp_model_in=itype_gscp_model, &
               rho_water=rho_w_model, rho_ice=rho_ice_model, &
               K_water=K_w_model, K_ice=K_ice_model, t0melt=t0_melt_model, &
               luse_muD_relation_rain = luse_muD_relation_rain_fwo, &
               t_in=t(:,:,:), rho_in=rho(:,:,:), &
               qc_in=qc(:,:,:), qr_in=qr(:,:,:), &
               qi_in=qi(:,:,:), qs_in=qs(:,:,:), &
               qg_in=qg(:,:,:), qh_in=qh(:,:,:), &
               qnc_in=qnc(:,:,:), qnr_in=qnr(:,:,:), &
               qni_in=qni(:,:,:), qns_in=qns(:,:,:), &
               qng_in=qng(:,:,:), qnh_in=qnh(:,:,:), &
               qgl_in=qgl(:,:,:), qhl_in=qhl(:,:,:), &
               z_radar=zh_radar_wp_tmp)
#endif
        ELSE
          IF (my_cart_id_fwo == 0) THEN
            errstring(:) = ' '
            WRITE (errstring,'(a,i1,a,i4)') &
              'Error itype_gscp_model: For itype_refl=',namlist%itype_refl, &
              ' gscp scheme not implemented in EMVORADO: itype_gscp_model = ', &
              itype_gscp_model
              CALL abort_run (my_cart_id_fwo, 35620, errstring, yzroutine)
          END IF
        END IF


        IF (wp == sp) THEN
          dbzvecstore(:,:,:) = REAL(zh_radar_wp_tmp, kind=dp)
        ELSE
          dbzvecstore(:,:,:) = zh_radar_wp_tmp
        END IF
        DEALLOCATE (zh_radar_wp_tmp)
        dbzextstore(:,:,:) = 0.0_dp

      CASE default
        IF (my_cart_id_fwo == 0) THEN
          WRITE (*,*) 'ERROR '//TRIM(yzroutine)//&
                      ': Unknown type of radar reflectivity calculation: ',namlist%itype_refl
        END IF
        dbzvecstore(:,:,:) = miss_value
        dbzextstore(:,:,:) = 0.0_dp
        zvvecstore(:,:,:) = miss_value
        rrhvvecstore(:,:,:) = miss_value_rhv
        irhvvecstore(:,:,:) = miss_value_rhv
        kdpvecstore(:,:,:) = 0.0_dp
        adpvecstore(:,:,:) = 0.0_dp
        zvhvecstore(:,:,:) = miss_value
      END SELECT

    END IF    ! ... IF ( .NOT.lgsp_fwo )

    IF (PRESENT(zh_radar)) THEN
      zh_radar(:,:,:) = dbzvecstore(:,:,:)
    END IF
    IF (PRESENT(ah_radar)) THEN
      ah_radar(:,:,:) = dbzextstore(:,:,:)
    END IF

    IF (PRESENT(zv_radar)) THEN
      zv_radar(:,:,:) = zvvecstore(:,:,:)
    END IF
    IF (PRESENT(rrhv_radar)) THEN
      rrhv_radar(:,:,:) = rrhvvecstore(:,:,:)
    END IF
    IF (PRESENT(irhv_radar)) THEN
      irhv_radar(:,:,:) = irhvvecstore(:,:,:)
    END IF
    IF (PRESENT(kdp_radar)) THEN
      kdp_radar(:,:,:) = kdpvecstore(:,:,:)
    END IF
    IF (PRESENT(adp_radar)) THEN
      adp_radar(:,:,:) = adpvecstore(:,:,:)
    END IF
    IF (PRESENT(zvh_radar)) THEN
      zvh_radar(:,:,:) = zvhvecstore(:,:,:)
    END IF

    DEALLOCATE (dbzvecstore, dbzextstore)
    DEALLOCATE (zvvecstore, rrhvvecstore, irhvvecstore, &
                kdpvecstore, adpvecstore)
    DEALLOCATE (zvhvecstore)

    IF (ldebug) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_cart_id_fwo

    RETURN
  END SUBROUTINE calc_dbz_vec_generic


  !==================================================================================
  !
  ! Routine to initialize the Mie-lookup tables.
  !
  ! Is called at model start from organize_data.f90, section 'start', and/or
  !   from organize_radar.f90, subroutine init_radar().
  !
  ! This is intended for performance increase in parallel runs, so that lookup
  !   tables for different reflectivity configurations can be computed at the
  !   same time by different processors, stored as files on disk and read in
  !   again by other processors, when those have to generate the lookup table
  !   for these configurations. The purpose is just to compute the lookup
  !   tables, not to compute any useful reflectivity fields.
  !
  ! It is simply done by calling calc_dbz_vec() for every different
  !   reflectivity configuration in the input TYPE(dbzcalc_params) vector
  !   "dbzm", which must have at least "nsets" elements.
  !   calc_dbz_vec() is called using its optional "lnostore=.true."
  !   and "l_use_neigh_tmax_melt=.false." switches. The former is necessary to
  !   inhibit, that in the next call the reflectivity is not calculated but
  !   taken from the internal storage field "dbzstore", and the latter is
  !   necessary to suppress a search in a certain horizontal
  !   neighbourhood for the parameter Tmax, which is the maximum temperature
  !   ( = minimum height) where melting hydrometeors exist in each vertical
  !   grid column.
  !
  !
  ! The house-keeping of the lookup tables (computation, file-write-and-read,
  !   internal storage) is fully done in module radar_mie_specint.
  !
  ! To distinguish the different configurations, hash values are calculated
  !   using "get_hashbase" and "hash_radar_lookup" from module radar_mie_utils.
  !
  !==================================================================================

  SUBROUTINE init_lookup_mie (nsets, dbzm, idom, ydir_lookup_read, ydir_lookup_write, ldebug, debugstring)

    IMPLICIT NONE

    INTEGER,                 INTENT(in) :: nsets   ! dimension of meta data structure vector dbz_meta
    TYPE(t_dbzcalc_params),  INTENT(in) :: dbzm(nsets)
    INTEGER,                 INTENT(in) :: idom    ! model domain identifier, needed in the call to calc_dbz_vec_modelgrid()
    LOGICAL, INTENT(in)                 :: ldebug
    CHARACTER (LEN=*), INTENT(in)       :: ydir_lookup_read, ydir_lookup_write, debugstring

    INTEGER                             :: i, j, magicnr, magicnr_old, magicnr_unique(nsets), &
                                           nsets_unique, work_pe, izerror, n_hydrometeors, &
                                           pe_start, pe_end
    TYPE(t_dbzcalc_params)              :: dbz_unique(nsets)
    CHARACTER (LEN=3000)                :: hashtext, cbaseconfig
    CHARACTER (LEN=*), PARAMETER        :: yzroutine = 'emvorado::init_lookup_mie()'
    LOGICAL                             :: ltestpattern_tmp

    IF (ldebug .OR. my_cart_id_fwo == 0) WRITE (*,*)  TRIM(yzroutine)//' from '// TRIM(debugstring) //&
                                                      ' on proc ', my_cart_id_fwo

    IF (num_compute_fwo > 1 .AND. &
        (ANY(dbzm%itype_refl == 1) .OR. ANY(dbzm%itype_refl > 4)) .AND. & !== 5)) .AND.
        ANY(dbzm%llookup_mie)) THEN

      ! 1) find and store the unique dbzm-sets:

      magicnr_old = -1365623035   ! hash2 for text(1:3000) = ' '
      nsets_unique = 0
      IF (ldebug .AND. my_cart_id_fwo == 0) &
        WRITE(*,'(a,i3,a)') '### Identifying unique dbzm-sets (out of ',nsets,' in total) ###'
      DO i=1, nsets

        IF ((dbzm(i)%itype_refl == 1 .OR. dbzm(i)%itype_refl > 4) .AND. & !== 5) .AND. &
            dbzm(i)%llookup_mie) THEN

          ! hash-value of the dbz-configuration to find all the different reflectivity configs
          ! (without taking into account any other differences in the microphysics configuration or the number of table nodes):
#ifdef TWOMOM_SB
          IF (itype_gscp_fwo  < 200) THEN
#endif
          cbaseconfig = get_hashbase ('versionXXX', dbzm(i)%lambda_radar, itype_gscp_fwo, &
                                      .TRUE., 2, &
                                      dbzm(i)%itype_refl, dbzm(i)%igraupel_type, dbzm(i)%itype_Dref_fmelt,&
                                      dbzm(i)%Tmeltbegin_i, dbzm(i)%meltdegTmin_i, dbzm(i)%Tmax_min_i, dbzm(i)%Tmax_max_i, &
                                      dbzm(i)%Tmeltbegin_s, dbzm(i)%meltdegTmin_s, dbzm(i)%Tmax_min_s, dbzm(i)%Tmax_max_s, &
                                      dbzm(i)%Tmeltbegin_g, dbzm(i)%meltdegTmin_g, dbzm(i)%Tmax_min_g, dbzm(i)%Tmax_max_g, &
                                      0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                      dbzm(i)%ctype_dryice_mie, dbzm(i)%ctype_wetice_mie, &
                                      dbzm(i)%ctype_drysnow_mie, dbzm(i)%ctype_wetsnow_mie, &
                                      dbzm(i)%ctype_drygraupel_mie, dbzm(i)%ctype_wetgraupel_mie, &
                                      '', '')
#ifdef TWOMOM_SB
          ELSE
          cbaseconfig = get_hashbase ('versionXXX', dbzm(i)%lambda_radar, itype_gscp_fwo, &
                                      .FALSE., 0, &
                                      dbzm(i)%itype_refl, dbzm(i)%igraupel_type, dbzm(i)%itype_Dref_fmelt,&
                                      dbzm(i)%Tmeltbegin_i, dbzm(i)%meltdegTmin_i, dbzm(i)%Tmax_min_i, dbzm(i)%Tmax_max_i, &
                                      dbzm(i)%Tmeltbegin_s, dbzm(i)%meltdegTmin_s, dbzm(i)%Tmax_min_s, dbzm(i)%Tmax_max_s, &
                                      dbzm(i)%Tmeltbegin_g, dbzm(i)%meltdegTmin_g, dbzm(i)%Tmax_min_g, dbzm(i)%Tmax_max_g, &
                                      dbzm(i)%Tmeltbegin_h, dbzm(i)%meltdegTmin_h, dbzm(i)%Tmax_min_h, dbzm(i)%Tmax_max_h, &
                                      dbzm(i)%ctype_dryice_mie, dbzm(i)%ctype_wetice_mie, &
                                      dbzm(i)%ctype_drysnow_mie, dbzm(i)%ctype_wetsnow_mie, &
                                      dbzm(i)%ctype_drygraupel_mie, dbzm(i)%ctype_wetgraupel_mie, &
                                      dbzm(i)%ctype_dryhail_mie, dbzm(i)%ctype_wethail_mie)
          END IF
#endif
          CALL hash_radar_lookup ( '', '', cbaseconfig, magicnr, hashtext)

          ! check if this hash is unique, and if yes, store it and the corresponding dbzm(i):
          IF (magicnr /= magicnr_old) THEN
            IF (nsets_unique > 0) THEN
              IF (.NOT.ANY(magicnr_unique(1:nsets_unique) == magicnr)) THEN
                nsets_unique = nsets_unique + 1
                dbz_unique(nsets_unique) = dbzm(i)
                magicnr_unique(nsets_unique) = magicnr
                IF (ldebug .AND. my_cart_id_fwo == 0) &
                  WRITE(*,'(a,i3,a,i10,a)') 'For set #',i,' magicnr=',magicnr,' ---> IS unique'
              ELSE
                IF (ldebug .AND. my_cart_id_fwo == 0) &
                  WRITE(*,'(a,i3,a,i10,a)') 'For set #',i,' magicnr=',magicnr,' ---> NOT unique'
              END IF
            ELSE
              nsets_unique = nsets_unique + 1
              dbz_unique(nsets_unique) = dbzm(i)
              magicnr_unique(nsets_unique) = magicnr
              IF (ldebug .AND. my_cart_id_fwo == 0) &
                WRITE(*,'(a,i3,a,i10,a)') 'For set #',i,' magicnr=',magicnr,' ---> IS unique'
            END IF
          ELSE
            IF (ldebug .AND. my_cart_id_fwo == 0) &
              WRITE(*,'(a,i3,a,i10,a)') 'For set #',i,' magicnr=',magicnr,' ---> NOT unique'
          END IF

          magicnr_old = magicnr

        END IF

      END DO

      ! 2) round robin calls to calc_dbz_vec() with the unique sets, at the same
      !    time deactivating the ltestpattern_hydrometeor-option.
      !    This well check if, for each unique config set, there already exists a lookup table for each
      !    hydrometeor type. If yes, nothing will happen. If no, the lookup tables will be created.
      !    Further parallelization over hydrometeor types or elements of the tables.

      pe_start = pe_start_lookupgen
      pe_end   = pe_end_lookupgen

      SELECT CASE (itype_mpipar_lookupgen)
        
      CASE (1)

        !.. Parallelization over unique dbz-sets and hydrometeors:
        
        IF (nsets_unique < 1) THEN
          n_hydrometeors = 0
        ELSE
          n_hydrometeors = SIZE(dbz_unique(1)%lhydrom_choice_testing)
        END IF
        DO i=1, nsets_unique
          DO j=1, n_hydrometeors
            work_pe = pe_start + MOD(j+(i-1)*n_hydrometeors-1, pe_end-pe_start+1)
            IF (my_cart_id_fwo == work_pe) THEN
              dbz_unique(i)%lhydrom_choice_testing(:) = .FALSE.
!#ifdef AUXOUT_OFFLINE
              dbz_unique(i)%lhydrom_choice_testing(j) = dbzm(i)%lhydrom_choice_testing(j)
!#else
!              dbz_unique(i)%lhydrom_choice_testing(j) = .TRUE.
!#endif
              IF (ldebug) WRITE (*,'(5(a,i0))') TRIM(yzroutine)//' from '// TRIM(debugstring) //&
                   ': generating Mie-lookup table for dbz config set ', i, ' (from ', nsets_unique, &
                   ') for hydrotype ', j, ' (from ', n_hydrometeors, ') on on proc ', my_cart_id_fwo
              ltestpattern_tmp = ltestpattern_hydrometeors
              ltestpattern_hydrometeors = .FALSE.
              CALL calc_dbz_vec_modelgrid(itlrad_dyn, itlrad_qx, idom, dbz_unique(i), .FALSE., ldebug, &
                   TRIM(ydir_lookup_read), TRIM(ydir_lookup_write), lnostore=.TRUE.)
              ltestpattern_hydrometeors = ltestpattern_tmp
            END IF
          END DO
        END DO
        
      CASE (2)
        
        !.. Parallelization over table elements:
        
        DO i=1, nsets_unique
          IF (lcompute_pe_fwo) THEN
!#ifndef AUXOUT_OFFLINE
!            dbz_unique(i)%lhydrom_choice_testing(:) = .TRUE.
!#endif
            dbz_unique(i)%lhydrom_choice_testing = dbzm(i)%lhydrom_choice_testing
            IF (ldebug .AND. my_cart_id_fwo == 0) WRITE (*,'(a,2(i0,a))') TRIM(yzroutine)//' from '// TRIM(debugstring) //&
                 ': generating Mie-lookup table for dbz config set ', i, ' (from ', nsets_unique, &
                 ') on all procs in parallel (itype_mpipar_lookupgen=2)'
            ltestpattern_tmp = ltestpattern_hydrometeors
            ltestpattern_hydrometeors = .FALSE.
            CALL calc_dbz_vec_modelgrid(itlrad_dyn, itlrad_qx, idom, dbz_unique(i), .FALSE., ldebug, &
                 TRIM(ydir_lookup_read), TRIM(ydir_lookup_write), lnostore=.TRUE.)
            ltestpattern_hydrometeors = ltestpattern_tmp
          END IF
        END DO
        
      END SELECT

      IF (ldebug) WRITE (*,*)  'Finished with '//TRIM(yzroutine)//' from '// TRIM(debugstring) //&
           ' on proc ', my_cart_id_fwo

      CALL get_runtime_timings (i_fwo_ini)

      ! MPI-Barrier for synchronisation, so that subsequent timing results
      !  are not biased. if the number of configs is not a multiple
      !  of the number of PEs, there will be load imbalance!
      IF (num_compute_fwo > 1 .AND. lcompute_pe_fwo) THEN
        CALL MPI_BARRIER(icomm_cart_fwo, izerror)
      END IF

      CALL get_runtime_timings (i_fwo_barrier)

    ELSE

      IF (ldebug) WRITE (*,*)  'Nothing done within '//TRIM(yzroutine)//' from '// TRIM(debugstring) //&
           ' on proc ', my_cart_id_fwo

    END IF

  END SUBROUTINE init_lookup_mie


  SUBROUTINE calc_fallspeed_vec_modelgrid (itl_dyn, itl_qx, idom, namlist_in, lwdbz, ldebug, vt_radar, vt_radar_wp)

!!$ if at a later time Mie-scattering will be considered,
!!$  add "l_use_neigh_tmax_melt, lnostore" to enable parallel
!!$  creation of lookup tables!

    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------

    IMPLICIT NONE


    INTEGER,                 INTENT(IN)    :: itl_dyn, itl_qx, idom
    TYPE(t_dbzcalc_params),  INTENT(IN)    :: namlist_in
    LOGICAL, INTENT(in)                    :: lwdbz    ! swith on/off reflectivity weighting for fallspeed
    LOGICAL, INTENT(in)                    :: ldebug

    REAL(kind=dp), OPTIONAL, INTENT(OUT)   :: vt_radar(:,:,:)
    REAL(kind=wp), OPTIONAL, INTENT(OUT)   :: vt_radar_wp(:,:,:)

    ! Helper strings for comparison of actual struct namlist with
    ! that of the last call to calc_dbz_vec, namlist_dbzcalc_vec_store:
    CHARACTER(len=cnmlstlen)               :: cnamnew, cnamold

    ! How to diagnose the N0 of snow PSD in 1-mom cases:
    INTEGER,                 PARAMETER     :: isnow_n0temp = 2

    LOGICAL, PARAMETER                     :: l_use_neigh_tmax_melt = .TRUE.

    ! Default value for neigh_tmax_melt (computation neighbourhood for the Tmax
    !  in the degree-of-melting parameterization of radar_mie_meltdegree.f90)
#ifdef __ICON__
    REAL(kind=dp), PARAMETER               :: neigh_tmax_melt_d = miss_value ! m
#else
    REAL(kind=dp), PARAMETER               :: neigh_tmax_melt_d = 5000.0_dp ! m
#endif

    TYPE(t_dbzcalc_params)                 :: namlist

    ! Flags and variables to determine if radar reflectivity calculation is necessary in this
    ! timestep or if it has been already calculated and we can use a SAVE storage array:
    CHARACTER(len=14)                      :: datetime
    LOGICAL                                :: zlnew_timestep
    TYPE(t_dbzcalc_params)                 :: namlist_tmp
    CHARACTER (LEN=*), PARAMETER           :: yzroutine = 'emvorado::calc_fallspeed_vec_modelgrid()'
    CHARACTER (LEN=300)                    :: errstring, dateiname
    LOGICAL                                :: zlcalc_vt_vec, zl_same_on_all_pes(1), zlcheck(num_compute_fwo)
    INTEGER                                :: k, mpierr, funit
    LOGICAL                                :: fileexist
    REAL(kind=dp)                          :: neigh_tmax_melt

    TYPE vtstore_t
      INTEGER                     :: ie  = -1
      INTEGER                     :: je  = -1
      INTEGER                     :: ke  = -1
      CHARACTER(len=14)           :: datetime = '99999999999999'
      INTEGER                     :: itl_dyn = -HUGE(1)
      INTEGER                     :: itl_qx  = -HUGE(1)
      CHARACTER(len=cnmlstlen)    :: cnamelist
!!$ The entire namelist does not have to be checked, as it is done here. Sufficient would be:
!!$  lambda_radar, tmeltbegin_s, meltdegtmin_s, tmeltbegin_g, meltdegtmin_g, tmeltbegin_h, meltdegtmin_h, lhydrom_choice_testing
      REAL (KIND=dp), ALLOCATABLE :: vt(:,:,:)
    END TYPE vtstore_t
    TYPE(vtstore_t), SAVE                  :: vtstore(ndoms_max_model)

    datetime = get_datetime_act ()

    CALL get_model_config_for_radar ( idom )   ! sets idom, ie_fwo, je_fwo, ke_fwo

    !.. Initializations for the beginning and each new timestep:
    IF (.NOT. ALLOCATED(vtstore(idom)%vt) .OR. vtstore(idom)%datetime /= datetime &
           .OR. vtstore(idom)%itl_dyn /= itl_dyn .OR. vtstore(idom)%itl_qx /= itl_qx) THEN

      zlnew_timestep = .TRUE.

      IF (.NOT.ALLOCATED(vtstore(idom)%vt)) THEN
        vtstore(idom)%ie = ie_fwo
        vtstore(idom)%je = je_fwo
        vtstore(idom)%ke = ke_fwo
        ALLOCATE(vtstore(idom)%vt(ie_fwo,je_fwo,ke_fwo))
        vtstore(idom)%vt = 0.0_dp
      END IF

      namlist_tmp = dbz_namlst_d
      ! change struct namlist_dbzcalc_vec_store a little bit,
      ! that it is not equal to that of the last timestep:
      namlist_tmp%ctype_wetsnow_mie = 'gagagagagaga'
      vtstore(idom)%cnamelist(:) = ' '
      WRITE(vtstore(idom)%cnamelist,*) namlist_tmp,  .FALSE.

      vtstore(idom)%datetime = datetime
      vtstore(idom)%itl_dyn  = itl_dyn
      vtstore(idom)%itl_qx   = itl_qx

    ELSE

      zlnew_timestep = .FALSE.

    END IF

    IF ( vtstore(idom)%ie /= ie_fwo .OR.  &
         vtstore(idom)%je /= je_fwo .OR.  &
         vtstore(idom)%ke /= ke_fwo ) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,i2,a,i6,a)') 'Not the same size for domain ',idom,' compared to previous call on PE ',my_cart_id_fwo,'!'
      CALL abort_run (my_cart_id_fwo, 20469, &
           'ERROR in '//TRIM(yzroutine)//': '//TRIM(ADJUSTL(errstring)), &
           'src_radar.f90, '//TRIM(yzroutine))
    END IF


    ! If a horizontal neighbourhood should be incorporated to
    !  search for the maximum temperature Tmax (= minimum height)
    !  at which melting hydrometeors exist in a vertical column,
    !  set it to the value of the parameter neigh_tmax_melt_d,
    !  otherwise set it to a value < 0.0 to disable the
    !  horizontal neighbourhood when searching for Tmax in each grid column.
    IF (l_use_neigh_tmax_melt) THEN
      neigh_tmax_melt = neigh_tmax_melt_d
    ELSE
      neigh_tmax_melt = miss_value
    END IF

    ! Write namelists to the helper strings, disregard station_id:
    cnamnew(:) = ' '
    namlist = namlist_in
    ! This is to avoid that different station IDs lead to a recalculation of vterm.
    !  although the rest of the configuration is the same:
    namlist%station_id = dbz_namlst_d%station_id
    WRITE(cnamnew,*) namlist, l_use_neigh_tmax_melt
    cnamold = vtstore(idom)%cnamelist

    ! store namlist of the actual call for the next call:
    vtstore(idom)%cnamelist(:) = ' '
    WRITE(vtstore(idom)%cnamelist,*) namlist,  l_use_neigh_tmax_melt

    IF ( .NOT.lgsp_fwo ) THEN

      IF (my_cart_id_fwo == 0) WRITE (*,*) '  VTERM is not calculated because of lgsp_fwo = .false.'

    ELSE

      ! .. Criterion for re-calculation of Vterm or copying the values from storage array:
      zlcalc_vt_vec = ( zlnew_timestep .OR. TRIM(cnamnew) /= TRIM(cnamold) )

      IF (ldebug) THEN

        ! .. Debug output of dbzparams-namelists (one file for each PE):
        dateiname(:) = ' '
        WRITE(dateiname,'("calc_vt_vec_err_pe-",i3.3,"_dom-",i2.2,".txt")') my_cart_id_fwo, idom

        CALL get_free_funit(funit)

        INQUIRE(file=TRIM(dateiname), exist=fileexist)
        IF (fileexist) THEN
          OPEN(funit , FILE=TRIM(dateiname), FORM='FORMATTED', STATUS='OLD', &
               action='WRITE', position='APPEND')
        ELSE
          OPEN(funit , FILE=TRIM(dateiname), FORM='FORMATTED', STATUS='REPLACE', &
               action='WRITE')
          fileexist = .TRUE.
        END IF

        WRITE (funit,'(a)') 'New: '//TRIM(cnamnew)
        WRITE (funit,'(a)') 'Old: '//TRIM(cnamold)

        CLOSE(funit)

        IF (num_compute_fwo > 1 .AND. lcompute_pe_fwo .AND. neigh_tmax_melt > 0.0_dp) THEN

          ! Check if we have the same settings of zlcalc_dbz_vec
          ! on all PROCESSORS. IF not, the PROGRAM might
          ! hang up, since there is an MPI-exchange of temperature in the reflectivity calculations
          ! for namlist%itype_refl <= 3 below:

          zl_same_on_all_pes(1) = zlcalc_vt_vec

          CALL MPI_GATHER (zl_same_on_all_pes, 1, MPI_LOGICAL,      &
                           zlcheck,            1, MPI_LOGICAL,      &
                           0, icomm_cart_fwo, mpierr)

          IF (my_cart_id_fwo == 0) THEN

            IF ( .NOT.ALL( zlcheck(2:num_compute_fwo) .EQV. zlcheck(1) ) ) THEN

              errstring(:) = ' '
              WRITE (errstring,'(a,i2,a)') 'Not the same settings on all PEs for DBZ re-calculation on domain ',idom,'!'
              CALL abort_run (my_cart_id_fwo, 20473, &
                   'ERROR in '//TRIM(yzroutine)//': '//TRIM(ADJUSTL(errstring)), &
                   'src_radar.f90, '//TRIM(yzroutine))

            END IF

          END IF

        END IF

      END IF  ! ldebug

      IF ( zlcalc_vt_vec ) THEN

        CALL get_model_variables (itl_dyn, itl_qx, idom)
        CALL get_model_hydrometeors (itl_qx, idom)

        IF (ltestpattern_hydrometeors) THEN
          CALL set_testpattern_hydrometeors_mg   ! mg = modelgrid
        END IF

        ! .. Set ilow_modelgrid, iup_modelgrid, jlow_modelgrid, jup_modelgrid, klow,_modelgrid kup_modelgrid
        !     for dbz-computations in this sub-domain including the halo PE boundaries:
        CALL get_loc_domain( idom, zlinc_boundlines = .TRUE. )

        IF (itype_gscp_fwo < 200) THEN
          ! Set the Tmax_XX_modelgrid parameter for the parameterization of the degree of melting
          !  of each hydrometeor category as 2D fields on the model grid:
          CALL initialize_tmax_1mom_vec_par(neigh_tmax_melt, namlist)
        ELSE
#ifdef TWOMOM_SB
          CALL initialize_tmax_2mom_vec_par(neigh_tmax_melt, namlist)
#endif
        END IF

        IF (my_cart_id_fwo == 0 .AND. ldebug) WRITE (*,*) '  Calculating VTERM ', &
             'with Rayleigh approximation and Oguchi refractive index, vectorized'

        IF (itype_gscp_fwo < 200) THEN

          CALL vtradar_rayleigh_oguchi_1mom_vec( my_cart_id_fwo, namlist%lambda_radar,  &
               itype_gscp_fwo, isnow_n0temp, &
               namlist%Tmeltbegin_i, namlist%meltdegTmin_i, &
               namlist%Tmeltbegin_s, namlist%meltdegTmin_s, &
               namlist%Tmeltbegin_g, namlist%meltdegTmin_g, &
               rho=rho(:,:,:), t=t(:,:,:), &
               qc=qc(:,:,:), qr=qr(:,:,:), qi=qi(:,:,:), &
               qs=qs(:,:,:), qg=qg(:,:,:), &
               Tmax_i=Tmax_i_modelgrid(:,:), Tmax_s=Tmax_s_modelgrid(:,:), &
               Tmax_g=Tmax_g_modelgrid(:,:), &
               Tmin_g=Tmin_g_modelgrid(:,:), &
               ilow=ilow_modelgrid, iup=iup_modelgrid, &
               jlow=jlow_modelgrid, jup=jup_modelgrid, &
               klow=klow_modelgrid, kup=kup_modelgrid, &
               lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, &
               lwdbz=lwdbz, vt_radar=vtstore(idom)%vt, &
               lhydrom_choice_testing=namlist%lhydrom_choice_testing )

#ifdef TWOMOM_SB
        ELSE

          CALL vtradar_rayleigh_oguchi_2mom_vec( my_cart_id_fwo, namlist%lambda_radar, &
               namlist%Tmeltbegin_i, namlist%meltdegTmin_i, &
               namlist%Tmeltbegin_s, namlist%meltdegTmin_s, &
               namlist%Tmeltbegin_g, namlist%meltdegTmin_g, &
               namlist%Tmeltbegin_h, namlist%meltdegTmin_h, &
               rho=rho(:,:,:), t=t(:,:,:), &
               qc=qc(:,:,:)  , qr=qr(:,:,:)  , qi=qi(:,:,:), &
               qs=qs(:,:,:)  , qg=qg(:,:,:)  , qh=qh(:,:,:), &
               qnc=qnc(:,:,:), qnr=qnr(:,:,:), qni=qni(:,:,:), &
               qns=qns(:,:,:), qng=qng(:,:,:), qnh=qnh(:,:,:), &
               qgl=qgl(:,:,:), qhl=qhl(:,:,:), &
               Tmax_i=Tmax_i_modelgrid(:,:), Tmax_s=Tmax_s_modelgrid(:,:), &
               Tmax_g=Tmax_g_modelgrid(:,:), Tmax_h=Tmax_h_modelgrid(:,:), &
               Tmin_g=Tmin_g_modelgrid(:,:), Tmin_h=Tmin_h_modelgrid(:,:), &
               ilow=ilow_modelgrid, iup=iup_modelgrid, &
               jlow=jlow_modelgrid, jup=jup_modelgrid, &
               klow=klow_modelgrid, kup=kup_modelgrid, &
               lalloc_qi=lalloc_qi, lalloc_qs=lalloc_qs, lalloc_qg=lalloc_qg, lalloc_qh=lalloc_qh, &
               luse_muD_relation_rain = luse_muD_relation_rain_fwo, &
               lwdbz=lwdbz, vt_radar=vtstore(idom)%vt, &
               lhydrom_choice_testing=namlist%lhydrom_choice_testing )
#endif

        END IF

        CALL finalize_tmax ()

      ELSE    ! ... IF (zlcalc_vt_vec)

        IF (my_cart_id_fwo == 0 .AND. ldebug) THEN

          WRITE (*,*) '  Taking previous VTERM with Rayleigh/Oguchi approximation'

        END IF

      END IF

    END IF

!$OMP PARALLEL PRIVATE(k)
    IF (PRESENT(vt_radar)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        vt_radar(:,:,k) = vtstore(idom)%vt(:,:,k)
      ENDDO
!$OMP END DO
    END IF

    IF (PRESENT(vt_radar_wp)) THEN
!$OMP DO
      DO k = 1, ke_fwo
        vt_radar_wp(:,:,k) = REAL(vtstore(idom)%vt(:,:,k), kind=wp)
      ENDDO
!$OMP END DO
    END IF
!$OMP END PARALLEL

    IF (ldebug) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_cart_id_fwo

    RETURN
  END SUBROUTINE calc_fallspeed_vec_modelgrid

END MODULE radar_mie_iface_cosmo_driver
