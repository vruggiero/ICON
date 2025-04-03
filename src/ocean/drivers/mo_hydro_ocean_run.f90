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

! Contains the main stepping routine the 3-dim hydrostatic ocean model.

!----------------------------
#include "icon_definitions.inc"
!----------------------------
MODULE mo_hydro_ocean_run
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_parallel_config,        ONLY: nproma
#ifdef _OPENACC
  USE mo_mpi,                      ONLY: i_am_accel_node, my_process_is_work
  USE openacc,                   ONLY : acc_is_present
#endif
  USE mo_impl_constants,         ONLY: max_char_length, success
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d
  USE mo_grid_config,            ONLY: n_dom
  USE mo_coupling_config,        ONLY: is_coupled_to_atmo, is_coupled_to_output
  USE mo_memory_log,             ONLY: memory_log_add
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    &  i_sea_ice, cfl_check, cfl_threshold, cfl_stop_on_violation,   &
    &  cfl_write, surface_module, &
    &  lswr_jerlov, lsediment_only, &
    &  Cartesian_Mixing, GMRedi_configuration, OceanReferenceDensity_inv, &
    &  atm_pressure_included_in_ocedyn, &
    &  vert_mix_type, vmix_pp, lcheck_salt_content, &
    &  use_age_tracer, & ! by_nils
    &  use_layers, & ! by_nils
    &  use_draftave_for_transport_h, &
    & vert_cor_type, use_tides, check_total_volume
  USE mo_ocean_nml,              ONLY: iforc_oce, Coupled_FluxFromAtmo
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: n_checkpoints, write_last_restart
  USE mo_run_config,             ONLY: dtime, ltimer, output_mode, debug_check_level
  USE mtime,                     ONLY: datetime, datetimeToString, deallocateDatetime
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total, timer_solve_ab,  &
    & timer_tracer_ab, timer_vert_veloc, timer_normal_veloc,     &
    & timer_upd_flx, timer_extra20, timers_level, &
    & timer_scalar_prod_veloc, timer_extra21, timer_extra22, timer_bgc_ini, &
    & timer_bgc_inv, timer_bgc_tot, timer_coupling
  USE mo_ocean_ab_timestepping,    ONLY: solve_free_surface_eq_ab, &
    &                                    calc_normal_velocity_ab,  &
    &                                    calc_vert_velocity,       &
    &                                    update_time_indices
  USE mo_ocean_types,              ONLY: t_hydro_ocean_state, &
    & t_operator_coeff, t_solvercoeff_singleprecision
  USE mo_ocean_tracer_transport_types, ONLY: t_tracer_collection, t_ocean_transport_state
  USE mo_ocean_math_operators,   ONLY: update_height_depdendent_variables, check_cfl_horizontal, check_cfl_vertical
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d, map_edges2edges_viacell_3d_const_z
  USE mo_ocean_tracer,           ONLY: advect_ocean_tracers
  USE mo_ocean_tracer_dev,       ONLY: advect_ocean_tracers_dev
  USE mo_ocean_nudging,          ONLY: nudge_ocean_tracers
  USE mo_restart,                ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_restart_nml_and_att,    ONLY: getAttributesForRestarting
  USE mo_key_value_store,        ONLY: t_key_value_store
  USE mo_ocean_surface_refactor, ONLY: update_ocean_surface_refactor
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_ice_fem_interface,      ONLY: ice_fem_init_vel_restart, ice_fem_update_vel_restart
  USE mo_sea_ice_types,          ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_ocean_physics,          ONLY: update_ho_params, update_ho_params_zstar
  USE mo_ocean_physics_types,    ONLY: t_ho_params
  USE mo_ocean_thermodyn,        ONLY: calc_potential_density, calculate_density, &
    & calculate_density_zstar
  USE mo_name_list_output,       ONLY: write_name_list_output
  USE mo_name_list_output_init,  ONLY: output_file
  USE mo_name_list_output_types, ONLY: t_output_file
  USE mo_ocean_diagnostics,      ONLY: calc_fast_oce_diagnostics, calc_psi
  USE mo_ocean_check_total_content, ONLY: check_total_salt_content, &
    & check_total_salt_content_zstar, check_accumulated_volume_difference
  USE mo_master_config,          ONLY: isRestart
  USE mo_master_control,         ONLY: get_my_process_name
  USE mo_time_config,            ONLY: time_config, t_time_config
  USE mo_util_dbg_prnt,          ONLY: dbg_print, debug_printValue
  USE mo_dbg_nml,                ONLY: idbg_mxmn
  USE mo_statistics
  USE mo_var_list
  USE mo_swr_absorption,         ONLY: jerlov_swr_absorption
  USE mo_ocean_hamocc_interface, ONLY: ocean_to_hamocc_interface
  USE mo_derived_variable_handling, ONLY: update_statistics
  USE mo_ocean_output
  USE mo_ocean_atmo_coupling,    ONLY: couple_ocean_toatmo_fluxes
  USE mo_hamocc_nml,             ONLY: l_cpl_co2
  USE mo_ocean_time_events,      ONLY: ocean_time_nextStep, isCheckpoint, isEndOfThisRun, newNullDatetime
  USE mo_ocean_ab_timestepping_mimetic, ONLY: clear_ocean_ab_timestepping_mimetic
  USE mo_ocean_ab_timestepping_zstar,  ONLY: calc_vert_velocity_bottomup_zstar, calc_normal_velocity_ab_zstar, &
    & solve_free_surface_eq_zstar, update_zstar_variables
  USE mo_ocean_tracer_zstar,     ONLY:advect_individual_tracers_zstar, advect_ocean_tracers_zstar
  USE mo_swr_absorption,         ONLY: subsurface_swr_absorption_zstar
  USE mo_ocean_tracer_dev,       ONLY: advect_ocean_tracers_GMRedi_zstar

  USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
  USE mo_ocean_layers,           ONLY: calc_layers ! by_nils
  USE mo_physical_constants,     ONLY: rho_ref, grav
  USE mo_ocean_pressure_bc_conditions,  ONLY: create_pressure_bc_conditions
  USE mo_ocean_state,            ONLY: transfer_ocean_state
  USE mo_ocean_ext_data,         ONLY: ext_data
  USE mo_ocean_state,            ONLY: v_base
  USE mo_ocean_nudging,          ONLY: ocean_nudge
  USE mo_fortran_tools,          ONLY: set_acc_host_or_device
  USE mo_ocean_age_tracer,       ONLY: calc_age_tracer
  USE mo_output_coupling,        ONLY: output_coupling

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: perform_ho_stepping
  PUBLIC  :: prepare_ho_stepping, end_ho_stepping
  PUBLIC  :: write_initial_ocean_timestep
  PUBLIC  :: update_time_g_n, update_time_indices

  CHARACTER(LEN=12)  :: str_module = 'HYDRO-ocerun'  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug
  !-------------------------------------------------------------------------


CONTAINS

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE prepare_ho_stepping(patch_3d, operators_coefficients, ocean_state, &
    & p_oce_sfc, p_as, sea_ice, ext_data, is_restart, solvercoeff_sp)

    TYPE(t_patch_3d ), INTENT(in)     :: patch_3d
    TYPE(t_ocean_surface)             :: p_oce_sfc
    TYPE(t_atmos_for_ocean),  INTENT(in)   :: p_as
    TYPE(t_operator_coeff), INTENT(INOUT) :: operators_coefficients
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE (t_sea_ice),   INTENT(inout) :: sea_ice
    TYPE(t_external_data), TARGET, INTENT(inout) :: ext_data
! !   TYPE (t_ho_params)                :: p_phys_param
    LOGICAL, INTENT(in)               :: is_restart
    TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp

    REAL(wp) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e)           !!
    LOGICAL :: lzacc

    lzacc = .FALSE.

    IF (is_restart .AND. is_coupled_to_atmo() ) THEN
        ! Initialize 10m Wind Speed from restart file when run in coupled mode
        p_as%fu10 = p_oce_sfc%Wind_Speed_10m
        p_as%pao = p_oce_sfc%sea_level_pressure
    ENDIF

    IF (is_restart .AND. is_coupled_to_atmo() .AND. l_cpl_co2 ) THEN
        ! Initialize CO" Mixing Ration from restart file when run in coupled mode with HAMOCC
        p_as%co2 = p_oce_sfc%CO2_Mixing_Ratio
    ENDIF
    IF (is_restart .AND. (i_ice_dyn == 1)) THEN
        ! Initialize u_ice, v_ice with sea_ice vals read from the restart file
        CALL ice_fem_init_vel_restart(patch_3d%p_patch_2D(1), sea_ice)
    END IF
!
!     IF (is_restart) THEN
!       ! Prepare ocean_state%p_prog, since it is needed by the sea ice model (e.g. wind stress computation)
!       IF ( i_sea_ice > 0 )         &
!       CALL update_height_depdendent_variables( patch_3d, ocean_state, ext_data, operators_coefficients, solvercoeff_sp)
!
!       CALL calc_scalar_product_veloc_3d( patch_3d,  &
!         & ocean_state%p_prog(nold(1))%vn,         &
!         & ocean_state%p_diag,                     &
!         & operators_coefficients)
!     ELSE
!     ENDIF
!
!     !    CALL update_diffusion_matrices( patch_3d,         &
!     !      & p_phys_param,                 &
!     !      & operators_coefficients%matrix_vert_diff_e,&
!     !      & operators_coefficients%matrix_vert_diff_c)
!
    CALL update_height_depdendent_variables( patch_3d, ocean_state, ext_data, operators_coefficients, solvercoeff_sp)

    !! Initialize stretch variable for zstar
    IF ( ( .NOT. isRestart()  ) .AND. ( vert_cor_type == 1 ) ) THEN
        ocean_state%p_prog(nold(1))%stretch_c = 1.0_wp
    ENDIF


    ! this is needed as initial condition or restart
!     CALL calc_scalar_product_veloc_3d( patch_3d,  &
!       & ocean_state(1)%p_prog(nold(1))%vn,         &
!       & ocean_state(1)%p_diag,                     &
!       & operators_coefficients)
!     CALL update_ho_params(patch_3d, ocean_state(1), p_as%fu10, sea_ice%concsum, p_phys_param, operators_coefficients)

    ! calculate the constant jerlov type subsurface sw absoption factor

    IF ( lswr_jerlov ) THEN

      CALL jerlov_swr_absorption(patch_3d, ocean_state)

    ENDIF

    lzacc = .TRUE.
    !$ACC UPDATE DEVICE(sea_ice%hi, sea_ice%conc, sea_ice%hs, sea_ice%zunderice, sea_ice%draftave) IF(lzacc)
    !$ACC UPDATE DEVICE(patch_3d%p_patch_1d(1)%prism_volume, ocean_state%p_diag%swsum, ocean_state%p_diag%swrab) IF(lzacc)
  END SUBROUTINE prepare_ho_stepping
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE end_ho_stepping()


  END SUBROUTINE end_ho_stepping
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Main stepping routine for call of hydrostatic ocean model
  !!
  !
!<Optimize:inUse>
  SUBROUTINE perform_ho_stepping( patch_3d, ocean_state, p_ext_data,    &
    & p_oce_sfc, p_phys_param,                                          &
    & p_as, p_atm_f, sea_ice, operators_coefficients,     &
    & solvercoeff_sp)

    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(t_ocean_surface)   ,  INTENT(inout)                         :: p_oce_sfc
    TYPE(t_ho_params)                       ,  INTENT(inout)         :: p_phys_param
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: p_atm_f
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_operator_coeff),   INTENT(inout), TARGET :: operators_coefficients
    TYPE(t_solvercoeff_singleprecision), INTENT(inout), TARGET :: solvercoeff_sp

    ! local variables
    INTEGER :: jstep, jg, return_status
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0 ! start counter for time loop
    REAL(wp) :: mean_height, old_mean_height
    REAL(wp) :: verticalMeanFlux(n_zlev+1)
    INTEGER :: level,ifiles,i,j
    REAL(wp) :: r

    REAL(wp) :: eta_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! Surface height after time step
    REAL(wp) :: stretch_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! stretch factor
    REAL(wp) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e)           !!

    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CLASS(t_RestartDescriptor), POINTER :: restartDescriptor
    CHARACTER(LEN = *), PARAMETER :: routine = 'mo_hydro_ocean_run:perform_ho_stepping'

    TYPE(datetime), POINTER             :: current_time     => NULL()
    LOGICAL :: lzacc

    lzacc = .FALSE.


#ifdef __COMPAD_DECLARATIONS__

#include "adify_oes_hydro_ocean_run_perform_ho_stepping_declarations.inc"

#endif /*  __COMPAD_DECLARATIONS__  */

    !------------------------------------------------------------------
    patch_2d      => patch_3d%p_patch_2d(1)

    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------
    IF (n_dom > 1 ) THEN
      CALL finish(routine, ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    patch_2d => patch_3d%p_patch_2d(jg)

    !------------------------------------------------------------------
    jstep0 = 0

    CALL getAttributesForRestarting(restartAttributes)
    ! get start counter for time loop from restart file:
    IF (restartAttributes%is_init) CALL restartAttributes%get("jstep", jstep0)
    IF (isRestart() .AND. mod(nold(jg),2) /=1 ) THEN
      ! swap the g_n and g_nm1
      CALL update_time_g_n(ocean_state(jg))
    ENDIF

    restartDescriptor => createRestartDescriptor(TRIM(get_my_process_name()) )

    ! local time var to be passed along, so the global is kept safe
    current_time => newNullDatetime()

    !------------------------------------------------------------------
    ! call the dynamical core: start the time loop
    !------------------------------------------------------------------
    CALL timer_start(timer_total)

    CALL transfer_ocean_state(patch_3d, operators_coefficients)

    !------------------------------------------------------------------
    !! Update stretch variable for zstar
    IF ( ( vert_cor_type == 1 ) ) THEN
#ifdef _OPENACC
        i_am_accel_node = my_process_is_work()    ! Activate GPUs
        lzacc = .TRUE.

        !$ACC DATA COPYIN(operators_coefficients%edge2cell_coeff_cc) &
        !$ACC   COPY(ocean_state(jg)%p_prog(nold(1))%eta_c, ocean_state(jg)%p_prog(nnew(1))%eta_c) &
        !$ACC   COPY(ocean_state(jg)%p_prog(nold(1))%stretch_c, stretch_e) IF(lzacc)
#endif
      CALL update_zstar_variables( patch_3d, ocean_state(jg), operators_coefficients, &
        & ocean_state(jg)%p_prog(nold(1))%eta_c, &
        & ocean_state(jg)%p_prog(nold(1))%stretch_c, stretch_e, lacc=lzacc)
#ifdef _OPENACC
        !$ACC END DATA

        lzacc = .FALSE.
        i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif
    ENDIF
    !------------------------------------------------------------------

    jstep = jstep0
    TIME_LOOP: DO

      IF(lsediment_only) THEN
        CALL sed_only_time_step()
      ELSE
        IF ( vert_cor_type == 0 ) THEN
          CALL ocean_time_step()
        ELSEIF ( vert_cor_type == 1 ) THEN
          CALL ocean_time_step_zstar()
        ENDIF
      END IF

      IF (isEndOfThisRun()) THEN
        ! leave time loop
        EXIT TIME_LOOP
      END IF

    ENDDO TIME_LOOP


    CALL clear_ocean_ab_timestepping_mimetic()

    CALL restartDescriptor%destruct()

    CALL deallocateDatetime(current_time)

    CALL timer_stop(timer_total)

  CONTAINS

    !-------------------------------------------------------------------------
    SUBROUTINE ocean_time_step()
        REAL(wp) :: total_salt, total_saltinseaice, total_saltinliquidwater
        INTEGER  :: blockNo, i
        LOGICAL  :: lzacc
        CHARACTER(LEN = *), PARAMETER :: routine = 'mo_hydro_ocean_run:ocean_time_step'

        lzacc = .FALSE.

        ! optional memory loggin
        CALL memory_log_add

        jstep = jstep + 1
        ! update model date and time mtime based
        current_time = ocean_time_nextStep()

        IF (idbg_mxmn >= 1 .OR. debug_check_level > 3) THEN

          CALL datetimeToString(current_time, datestring)
          WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
          CALL message (TRIM(routine), message_text)

        END IF

!        IF (lcheck_salt_content) CALL check_total_salt_content(100,ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2), patch_2d, &
!         ocean_state(jg)%p_prog(nold(1))%h(:,:), patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!         sea_ice, 0)

#ifdef _OPENACC
        i_am_accel_node = my_process_is_work()    ! Activate GPUs
        lzacc = .TRUE.
#endif

        start_detail_timer(timer_extra22,6)
        CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, &
                                                 solvercoeff_sp, lacc=lzacc)
        stop_detail_timer(timer_extra22,6)

        start_timer(timer_scalar_prod_veloc,2)
        CALL calc_scalar_product_veloc_3d( patch_3d,  &
          & ocean_state(jg)%p_prog(nold(1))%vn,       &
          & ocean_state(jg)%p_diag,                   &
          & operators_coefficients,                   &
          & lacc=lzacc)
        stop_timer(timer_scalar_prod_veloc,2)

        !In case of a time-varying forcing:
        ! update_surface_flux or update_ocean_surface has changed p_prog(nold(1))%h, SST and SSS
        start_timer(timer_upd_flx,3)
        IF (lcheck_salt_content) CALL check_total_salt_content(105,ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2), patch_2d, &
         ocean_state(jg)%p_prog(nold(1))%h(:,:), patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
         sea_ice, 0)

        CALL update_ocean_surface_refactor( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, p_oce_sfc, &
             & current_time, operators_coefficients, lacc = lzacc )

        IF (lcheck_salt_content) CALL check_total_salt_content(110,ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2), patch_2d, &
         ocean_state(jg)%p_prog(nold(1))%h(:,:), patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
         sea_ice, 0)

        stop_timer(timer_upd_flx,3)

        start_detail_timer(timer_extra22,4)
        CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, &
                                                 solvercoeff_sp, lacc=lzacc)
        stop_detail_timer(timer_extra22,4)

        !--------------------------------------------------------------------------
        ! calculate in situ density here, as it may be used fotr the tides load
        CALL calculate_density( patch_3d,                         &
          & ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,:),      &
          & ocean_state(jg)%p_diag%rho(:,:,:), lacc = lzacc)

        !--------------------------------------------------------------------------
        CALL create_pressure_bc_conditions(patch_3d,ocean_state(jg), p_as, sea_ice, current_time, lacc=lzacc)
        !------------------------------------------------------------------------

  !       IF (timers_level > 2) CALL timer_start(timer_scalar_prod_veloc)
  !       CALL calc_scalar_product_veloc_3d( patch_3d,  &
  !         & ocean_state(jg)%p_prog(nold(1))%vn,         &
  !         & ocean_state(jg)%p_diag,                     &
  !         & operators_coefficients)
  !       IF (timers_level > 2) CALL timer_stop(timer_scalar_prod_veloc)

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('on entry: h-old'           ,ocean_state(jg)%p_prog(nold(1))%h ,str_module,idt_src, &
          & patch_2d%cells%owned )
        CALL dbg_print('on entry: h-new'           ,ocean_state(jg)%p_prog(nnew(1))%h ,str_module,idt_src, &
          & patch_2d%cells%owned )
        CALL dbg_print('HydOce: ScaProdVel kin'    ,ocean_state(jg)%p_diag%kin        ,str_module,idt_src, &
          & patch_2d%cells%owned )
        CALL dbg_print('HydOce: ScaProdVel ptp_vn' ,ocean_state(jg)%p_diag%ptp_vn     ,str_module,idt_src, &
          & patch_2d%edges%owned )
        CALL dbg_print('HydOce: fu10'              ,p_as%fu10                         ,str_module,idt_src, &
          & in_subset=patch_2d%cells%owned)
        CALL dbg_print('HydOce: concsum'           ,sea_ice%concsum                   ,str_module,idt_src, &
          & in_subset=patch_2d%cells%owned)

        !---------------------------------------------------------------------
        !by_ogut: added p_oce_sfc
        CALL update_ho_params(patch_3d, ocean_state(jg), p_as%fu10, sea_ice%concsum, p_phys_param, operators_coefficients, &
                              p_atm_f, p_oce_sfc, lacc=lzacc)

        !------------------------------------------------------------------------
        IF (debug_check_level > 5) THEN
          CALL horizontal_mean(values=ocean_state(jg)%p_prog(nold(1))%h(:,:), weights=patch_2d%cells%area(:,:), &
            & in_subset=patch_2d%cells%owned, mean=old_mean_height)
        END IF
        !------------------------------------------------------------------------
        IF (lcheck_salt_content) CALL check_total_salt_content(125,ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2), patch_2d, &
         ocean_state(jg)%p_prog(nold(1))%h(:,:), patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
         sea_ice,0)
        ! solve for new free surface
        start_timer(timer_solve_ab,1)
        CALL solve_free_surface_eq_ab (patch_3d, ocean_state(jg), p_ext_data(jg), &
          & p_as, p_oce_sfc, p_phys_param, jstep, operators_coefficients, solvercoeff_sp, return_status, lacc = lzacc)!, p_int(jg))
        IF (return_status /= 0) THEN
         CALL output_ocean(              &
           & patch_3d=patch_3d,          &
           & ocean_state=ocean_state,    &
           & this_datetime=current_time, &
           & surface_fluxes=p_oce_sfc, &
           & sea_ice=sea_ice,            &
           & jstep=jstep, jstep0=jstep0, &
           & force_output=.true.,        &
           & lacc=lzacc)
          CALL finish(routine, 'solve_free_surface_eq_ab  returned error')
        ENDIF

        stop_timer(timer_solve_ab,1)

!         IF (lcheck_salt_content) CALL check_total_salt_content(126,ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2), patch_2d, &
!          ocean_state(jg)%p_prog(nnew(1))%h(:,:), patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!          sea_ice,0)

        ! fill diffusion coefficients
        DO i = 1, ocean_state(jg)%p_prog(nold(1))%tracer_collection%no_of_tracers
          !$ACC ENTER DATA COPYIN(ocean_state(jg)%p_prog(nold(1))%tracer_collection%tracer(i)%concentration) IF(lzacc)
          !$ACC ENTER DATA COPYIN(ocean_state(jg)%p_prog(nnew(1))%tracer_collection%tracer(i)%concentration) IF(lzacc)
        END DO

        !------------------------------------------------------------------------
        ! Step 4: calculate final normal velocity from predicted horizontal
        ! velocity vn_pred and updated surface height
        start_timer(timer_normal_veloc,4)
        CALL calc_normal_velocity_ab(patch_3d, ocean_state(jg),&
          & operators_coefficients, solvercoeff_sp,  p_ext_data(jg), p_phys_param, lacc = lzacc)
        stop_timer(timer_normal_veloc,4)

        !------------------------------------------------------------------------
        ! Step 5: calculate vertical velocity and mass_flx_e from continuity equation under
        ! incompressiblity condition in the non-shallow-water case
        IF ( iswm_oce /= 1 ) THEN
          start_timer(timer_vert_veloc,4)
          CALL calc_vert_velocity( patch_3d, ocean_state(jg),operators_coefficients, lacc = lzacc)
          stop_timer(timer_vert_veloc,4)
        ELSE

#ifdef _OPENACC
          IF (lzacc) CALL finish(routine, 'OpenACC version currently not tested/validated')
#endif

          CALL map_edges2edges_viacell_3d_const_z( patch_3d, ocean_state(jg)%p_diag%vn_time_weighted, operators_coefficients, &
              & ocean_state(jg)%p_diag%mass_flx_e)
        ENDIF

        !------------------------------------------------------------------------
        IF (check_total_volume) THEN
          CALL check_accumulated_volume_difference(1, patch_2d, &
            & ocean_state(jg)%p_prog(nnew(1))%h(:,:), ocean_state(jg)%p_prog(nold(1))%h(:,:))
        ENDIF

        !------------------------------------------------------------------------
        IF (idbg_mxmn >= 2 .OR. debug_check_level > 5) THEN
          CALL horizontal_mean(values=ocean_state(jg)%p_prog(nnew(1))%h(:,:), weights=patch_2d%cells%area(:,:), &
            & in_subset=patch_2d%cells%owned, mean=mean_height)
          CALL debug_printValue(description="Mean Height", val=mean_height, detail_level=2)
        ENDIF
        IF (debug_check_level > 5 .AND. idbg_mxmn >= 2) THEN
          ! check difference from old_mean_height
  !         CALL debug_printValue(description="Old/New Mean Height", val=old_mean_height, &
  !           & value1=mean_height, value2=(mean_height-old_mean_height) / old_mean_height, &
  !           & detail_level=2)
          CALL debug_printValue(description="Old/New Mean Height", &
            & val=old_mean_height, value1=mean_height, detail_level=2)
          ! check if vertical and horizontal fluxes add to 0
  !         ocean_state(jg)%p_diag%w
          CALL horizontal_mean(values=ocean_state(jg)%p_diag%w, weights=patch_2d%cells%area(:,:), &
            & in_subset=patch_2d%cells%owned, mean=verticalMeanFlux, start_level=2, end_level=n_zlev)

          DO level=2, n_zlev-1
            CALL debug_printValue(description="Mean vertical flux at", val=REAL(level,wp),  &
              & value1=verticalMeanFlux(level), detail_level=2)
          ENDDO
        END IF

!         IF (lcheck_salt_content) CALL check_total_salt_content(130,ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2), patch_2d, &
!          ocean_state(jg)%p_prog(nnew(1))%h(:,:), patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!          sea_ice,0)
        !------------------------------------------------------------------------

        ! fill diffusion coefficients
        DO i = 1, ocean_state(jg)%p_prog(nold(1))%tracer_collection%no_of_tracers
          ocean_state(jg)%p_prog(nold(1))%tracer_collection%tracer(i)%hor_diffusion_coeff => p_phys_param%TracerDiffusion_coeff(:,:,:,i)
          ocean_state(jg)%p_prog(nold(1))%tracer_collection%tracer(i)%ver_diffusion_coeff => p_phys_param%a_tracer_v(:,:,:,i)
        END DO

#ifdef _OPENACC
        IF (GMRedi_configuration /= Cartesian_Mixing ) THEN
           IF (lzacc) CALL finish(routine, 'OpenACC version with GMRedi currently not tested/validated')
        END IF
#endif

        CALL tracer_transport(patch_3d, ocean_state(jg), p_as, sea_ice, p_oce_sfc, &
          & p_phys_param, operators_coefficients, current_time, lacc=lzacc)

        IF (use_age_tracer) THEN
          !$ACC DATA COPYIN(p_as%pao) IF(lzacc)
          CALL calc_age_tracer(patch_3d, ocean_state(jg), sea_ice, lacc=lzacc)
          !$ACC END DATA
        ENDIF

        IF (lcheck_salt_content) CALL check_total_salt_content(140,ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,2), patch_2d, &
          ocean_state(jg)%p_prog(nnew(1))%h(:,:), patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
          sea_ice,0)
        !----------------------------------------------------------------------

      !------------------------------------------------------------------------
        ! Optional : nudge temperature and salinity
        IF (no_tracer>=1) THEN
          CALL nudge_ocean_tracers( patch_3d, ocean_state(jg), lacc=lzacc)
        ENDIF
        !------------------------------------------------------------------------
        ! perform accumulation for special variables
        start_detail_timer(timer_extra20,5)

        IF (no_tracer>=1) THEN

          CALL calc_potential_density( patch_3d,                            &
            & ocean_state(jg)%p_prog(nold(1))%tracer,                       &
            & ocean_state(jg)%p_diag%rhopot, lacc=lzacc )

          ! calculate diagnostic barotropic stream function
          CALL calc_psi (patch_3d, ocean_state(jg)%p_diag%u(:,:,:),         &
            & patch_3D%p_patch_1d(1)%prism_thick_c(:,:,:),                  &
            & ocean_state(jg)%p_diag%u_vint, current_time, lacc=lzacc)


          CALL dbg_print('calc_psi: u_vint' ,ocean_state(jg)%p_diag%u_vint, str_module, 3, in_subset=patch_2d%cells%owned)

        ! calculate diagnostic barotropic stream function with vn
    !  not yet mature
    !   CALL calc_psi_vn (patch_3d, ocean_state(jg)%p_prog(nold(1))%vn,   &
    !     & patch_3D%p_patch_1d(1)%prism_thick_e(:,:,:),                  &
    !     & operators_coefficients,                                       &
    !     & ocean_state(jg)%p_diag%u_vint, ocean_state(jg)%p_diag%v_vint, current_time)
    !   CALL dbg_print('calc_psi_vn: u_vint' ,ocean_state(jg)%p_diag%u_vint, str_module, 5, in_subset=patch_2d%cells%owned)
    !   CALL dbg_print('calc_psi_vn: v_vint' ,ocean_state(jg)%p_diag%v_vint, str_module, 5, in_subset=patch_2d%cells%owned)
        ENDIF

        DO i = 1, ocean_state(jg)%p_prog(nold(1))%tracer_collection%no_of_tracers
          !$ACC EXIT DATA COPYOUT(ocean_state(jg)%p_prog(nold(1))%tracer_collection%tracer(i)%concentration) IF(lzacc)
          !$ACC EXIT DATA COPYOUT(ocean_state(jg)%p_prog(nnew(1))%tracer_collection%tracer(i)%concentration) IF(lzacc)
        END DO

#ifdef _OPENACC
        !$ACC UPDATE SELF(patch_3D%p_patch_1d(1)%prism_center_dist_c, patch_3D%p_patch_1d(1)%prism_volume) &
        !$ACC   SELF(patch_3D%p_patch_1d(1)%depth_cellmiddle, patch_3D%p_patch_1d(1)%depth_cellinterface) &
        !$ACC   SELF(patch_3D%p_patch_1d(1)%prism_thick_c, patch_3D%p_patch_1d(1)%prism_thick_e) &
        !$ACC   SELF(patch_3D%p_patch_1d(1)%inv_prism_thick_c, patch_3D%p_patch_1d(1)%inv_prism_thick_e) &
        !$ACC   SELF(patch_3D%p_patch_1d(1)%inv_prism_center_dist_c, patch_3D%p_patch_1d(1)%inv_prism_center_dist_e) IF(lzacc)

        !$ACC UPDATE SELF(sea_ice%hi, sea_ice%conc, sea_ice%concsum, sea_ice%hs, sea_ice%vn_e, sea_ice%vol) &
        !$ACC   SELF(sea_ice%zunderice, sea_ice%draftave, sea_ice%totalsnowfall) IF(lzacc)

        !$ACC UPDATE SELF(ocean_state(jg)%p_prog(nnew(1))%h, ocean_state(jg)%p_prog(nnew(1))%h, ocean_state(jg)%p_prog(nnew(1))%vn) &
        !$ACC   SELF(ocean_state(jg)%p_prog(nnew(1))%tracer) IF(lzacc)

        !$ACC UPDATE SELF(ocean_state(jg)%p_diag%kin, ocean_state(jg)%p_diag%rho) &
        !$ACC   SELF(ocean_state(jg)%p_diag%vn_pred, ocean_state(jg)%p_diag%zgrad_rho) &
        !$ACC   SELF(ocean_state(jg)%p_diag%div_mass_flx_c, ocean_state(jg)%p_diag%mass_flx_e) &
        !$ACC   SELF(ocean_state(jg)%p_diag%vn_time_weighted, ocean_state(jg)%p_diag%u_vint) &
        !$ACC   SELF(ocean_state(jg)%p_diag%rhopot, ocean_state(jg)%p_diag%w, ocean_state(jg)%p_diag%swrab) IF(lzacc)

        !$ACC UPDATE SELF(p_oce_sfc%heatflux_total, p_oce_sfc%frshflux_volumetotal) &
        !$ACC   SELF(ocean_state(jg)%p_diag%delta_thetao, ocean_state(jg)%p_diag%delta_so) &
        !$ACC   SELF(ocean_state(jg)%p_diag%delta_snow, ocean_state(jg)%p_diag%delta_ice) IF(lzacc)

        !$ACC UPDATE SELF(p_oce_sfc%FrshFlux_Precipitation, p_oce_sfc%FrshFlux_Evaporation, p_oce_sfc%FrshFlux_Runoff) &
        !$ACC   SELF(p_oce_sfc%FrshFlux_VolumeIce) &
        !$ACC   SELF(p_oce_sfc%FrshFlux_TotalOcean, p_oce_sfc%FrshFlux_TotalIce) IF(lzacc)

        lzacc = .FALSE.
        i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif

        CALL fill_auxiliary_diagnostics(patch_3d, ocean_state(1))

        CALL calc_fast_oce_diagnostics( patch_2d, &
            & patch_3d, &
            & ocean_state(1), &
            & patch_3d%p_patch_1d(1)%dolic_c, &
            & patch_3d%p_patch_1d(1)%prism_thick_c, &
            & patch_3d%p_patch_1d(1)%zlev_m, &
            & ocean_state(jg)%p_diag, &
            & ocean_state(jg)%p_prog(nnew(1))%h, &
            & ocean_state(jg)%p_prog(nnew(1))%vn, &
            & ocean_state(jg)%p_prog(nnew(1))%tracer, &
            & p_atm_f, &
            & p_oce_sfc, &
            & sea_ice)

        stop_detail_timer(timer_extra20,5)

#ifdef _OPENACC
        i_am_accel_node = my_process_is_work()    ! Activate GPUs
        lzacc = .TRUE.

        !$ACC UPDATE SELF(p_as%fu10, p_as%pao) IF(lzacc)
#endif
        CALL update_statistics

        CALL output_ocean( patch_3d, ocean_state, &
          &                current_time,              &
          &                p_oce_sfc,             &
          &                sea_ice,                 &
          &                jstep, jstep0, lacc=lzacc)

        IF ( is_coupled_to_output() ) THEN
          IF (ltimer) CALL timer_start(timer_coupling)
          CALL output_coupling(patch_3d%wet_c)
          IF (ltimer) CALL timer_stop(timer_coupling)
        END IF

#ifdef _OPENACC
        lzacc = .FALSE.
        i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif

        ! send and receive coupling fluxes for ocean at the end of time stepping loop
        IF (iforc_oce == Coupled_FluxFromAtmo) THEN  !  14

          IF ( is_coupled_to_atmo() ) THEN
            IF (ltimer) CALL timer_start(timer_coupling)
            CALL couple_ocean_toatmo_fluxes( &
              patch_3D, ocean_state(jg), sea_ice, p_atm_f, p_as)
            IF (ltimer) CALL timer_stop(timer_coupling)
          END IF

          ! copy fluxes updated in coupling from p_atm_f into p_oce_sfc
          p_oce_sfc%FrshFlux_Precipitation = p_atm_f%FrshFlux_Precipitation
          p_oce_sfc%FrshFlux_Evaporation   = p_atm_f%FrshFlux_Evaporation
          p_oce_sfc%FrshFlux_SnowFall      = p_atm_f%FrshFlux_SnowFall
          p_oce_sfc%HeatFlux_Total         = p_atm_f%HeatFlux_Total
          p_oce_sfc%HeatFlux_ShortWave     = p_atm_f%HeatFlux_ShortWave
          p_oce_sfc%HeatFlux_Longwave      = p_atm_f%HeatFlux_Longwave
          p_oce_sfc%HeatFlux_Sensible      = p_atm_f%HeatFlux_Sensible
          p_oce_sfc%HeatFlux_Latent        = p_atm_f%HeatFlux_Latent
          p_oce_sfc%FrshFlux_Runoff        = p_atm_f%FrshFlux_Runoff
          IF ( l_cpl_co2 ) THEN
            p_oce_sfc%CO2_Mixing_Ratio     = p_as%co2
          ENDIF

        ENDIF

        ! copy atmospheric wind speed of coupling from p_as%fu10 into forcing to be written by restart
        p_oce_sfc%Wind_Speed_10m(:,:) = p_as%fu10(:,:)
        p_oce_sfc%sea_level_pressure(:,:) = p_as%pao(:,:)

!        IF (lcheck_salt_content) CALL check_total_salt_content(150,ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,2), patch_2d, &
!         ocean_state(jg)%p_prog(nnew(1))%h(:,:), patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!         sea_ice,0)

#ifdef _OPENACC
        i_am_accel_node = my_process_is_work()    ! Activate GPUs
        lzacc = .TRUE.

        !$ACC UPDATE DEVICE(p_oce_sfc%Wind_Speed_10m, p_oce_sfc%sea_level_pressure) IF(lzacc)

        lzacc = .FALSE.
        i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif

        start_detail_timer(timer_extra21,5)

        ! Shift time indices for the next loop
        ! this HAS to ge into the restart files, because the start with the following loop
        CALL update_time_indices(jg)
        !$ACC UPDATE DEVICE(nold, nnew)

        ! update intermediate timestepping variables for the tracers
        CALL update_time_g_n(ocean_state(jg))
!        IF (lcheck_salt_content) CALL check_total_salt_content(160,ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2), patch_2d, &
!         ocean_state(jg)%p_prog(nold(1))%h(:,:), patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!         sea_ice,0)

        ! by_nils: layer diagnostic
        IF (use_layers) THEN
          CALL calc_layers(patch_3d, ocean_state(jg), p_oce_sfc, operators_coefficients, p_phys_param)
        END IF

        ! check whether time has come for writing restart file
        IF (isCheckpoint()) THEN
          IF (.NOT. output_mode%l_none ) THEN
            !
            ! For multifile restart (restart_write_mode = "joint procs multifile")
            ! the domain flag must be set to .TRUE. in order to activate the domain,
            ! even though we have one currently in the ocean. Without this the
            ! processes won't write out their data into a the patch restart files.
            !
            patch_2d%ldom_active = .TRUE.
#ifdef _OPENACC
            i_am_accel_node = my_process_is_work()    ! Activate GPUs
            lzacc = .TRUE.
#endif
            !
            IF (i_ice_dyn == 1) CALL ice_fem_update_vel_restart(patch_2d, sea_ice, lacc=lzacc) ! write FEM vel to restart or checkpoint file
            CALL restartDescriptor%updatePatch(patch_2d, &
                                              &opt_nice_class=1, &
                                              &opt_ocean_zlevels=n_zlev, &
                                              &opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
                                              &opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
            CALL restartDescriptor%writeRestart(current_time, jstep)
#ifdef _OPENACC
            lzacc = .FALSE.
            i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif
          END IF
        END IF

        stop_detail_timer(timer_extra21,5)

        IF (isEndOfThisRun()) THEN
          ! leave time loop
          RETURN
        END IF

        ! check cfl criterion
        IF (cfl_check) THEN
          CALL check_cfl_horizontal(ocean_state(jg)%p_prog(nnew(1))%vn, &
            & patch_2d%edges%inv_dual_edge_length, &
            & dtime, &
            & patch_2d%edges%ALL, &
            & cfl_threshold, &
            & ocean_state(jg)%p_diag%cfl_horz, &
            & cfl_stop_on_violation,&
            & cfl_write)
          CALL check_cfl_vertical(ocean_state(jg)%p_diag%w, &
            & patch_3d%p_patch_1d(1)%prism_center_dist_c, &
            & dtime, &
            & patch_2d%cells%ALL,&
            & cfl_threshold, &
            & ocean_state(jg)%p_diag%cfl_vert, &
            & cfl_stop_on_violation,&
            & cfl_write)
        END IF

    END SUBROUTINE ocean_time_step

    !-------------------------------------------------------------------------
    SUBROUTINE ocean_time_step_zstar()
        INTEGER  :: blockNo, i
        LOGICAL  :: lzacc
        CHARACTER(LEN = *), PARAMETER :: routine = 'mo_hydro_ocean_run:ocean_time_step_zstar'

        lzacc = .FALSE.

        ! fill transport state
        ocean_state(jg)%transport_state%patch_3d    => patch_3d

        ! optional memory loggin
        CALL memory_log_add

        jstep = jstep + 1
        ! update model date and time mtime based
        current_time = ocean_time_nextStep()

        IF (idbg_mxmn >= 1 .OR. debug_check_level > 3) THEN

          CALL datetimeToString(current_time, datestring)
          WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
          CALL message (TRIM(routine), message_text)

        END IF

!        IF (lcheck_salt_content) CALL check_total_salt_content_zstar(110, &
!          & ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2), patch_2d, &
!          & ocean_state(jg)%p_prog(nold(1))%stretch_c(:,:), &
!          & patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:), sea_ice, p_oce_sfc)

#ifdef _OPENACC
        i_am_accel_node = my_process_is_work()    ! Activate GPUs
        lzacc = .TRUE.

        !$ACC DATA COPYIN(ocean_state(jg)%p_aux%bc_tides_potential, ocean_state(jg)%p_aux%bc_SAL_potential) &
        !$ACC   COPY(ocean_state(jg)%p_prog(nold(1))%eta_c, ocean_state(jg)%p_prog(nnew(1))%eta_c) &
        !$ACC   COPY(ocean_state(jg)%p_prog(nold(1))%stretch_c, ocean_state(jg)%p_prog(nold(1))%h) &
        !$ACC   COPY(stretch_e) IF(lzacc)

        !$ACC UPDATE DEVICE(p_oce_sfc%data_surfRelax_Salt, p_oce_sfc%data_surfRelax_Temp) &
        !$ACC   DEVICE(p_oce_sfc%HeatFlux_ShortWave, p_oce_sfc%HeatFlux_LongWave, p_oce_sfc%HeatFlux_Sensible) &
        !$ACC   DEVICE(p_oce_sfc%HeatFlux_Latent, p_oce_sfc%HeatFlux_Total, p_oce_sfc%FrshFlux_IceSalt) &
        !$ACC   DEVICE(p_oce_sfc%FrshFlux_VolumeIce, p_oce_sfc%FrshFlux_TotalIce, p_oce_sfc%FrshFlux_VolumeTotal) &
        !$ACC   DEVICE(p_oce_sfc%FrshFlux_TotalSalt, p_oce_sfc%top_dilution_coeff) &
        !$ACC   DEVICE(p_oce_sfc%TopBC_WindStress_u, p_oce_sfc%TopBC_WindStress_v, p_oce_sfc%TopBC_WindStress_cc) &
        !$ACC   DEVICE(p_oce_sfc%Wind_Speed_10m, p_oce_sfc%FrshFlux_Precipitation, p_oce_sfc%cellThicknessUnderIce) &
        !$ACC   DEVICE(p_oce_sfc%FrshFlux_Evaporation, p_oce_sfc%FrshFlux_Runoff, p_oce_sfc%FrshFlux_TotalOcean) &
        !$ACC   DEVICE(p_oce_sfc%FrshFlux_SnowFall, p_oce_sfc%TempFlux_Relax) &
        !$ACC   DEVICE(p_oce_sfc%HeatFlux_Relax, p_oce_sfc%SaltFlux_Relax, p_oce_sfc%FrshFlux_Relax) IF(lzacc)

        !$ACC UPDATE DEVICE(ocean_state(jg)%p_diag%rsdoabsorb) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%swsum, ocean_state(jg)%p_diag%swrab) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%heatabs, ocean_state(jg)%p_diag%heatflux_rainevaprunoff) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%delta_ice, ocean_state(jg)%p_diag%delta_snow) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%delta_thetao, ocean_state(jg)%p_diag%delta_so) &
        !$ACC   DEVICE(ocean_state(jg)%p_prog(nold(1))%tracer, ocean_state(jg)%p_prog(nold(1))%vn) IF(lzacc)

        !$ACC UPDATE DEVICE(p_as%pao) &
        !$ACC   DEVICE(p_as%topBoundCond_windStress_u, p_as%topBoundCond_windStress_v, p_as%u, p_as%v) &
        !$ACC   DEVICE(p_as%tafo, p_as%ftdew, p_as%fu10, p_as%fclou, p_as%fswr) &
        !$ACC   DEVICE(p_as%FrshFlux_Precipitation, p_as%FrshFlux_Runoff, p_as%data_surfRelax_Temp) IF(lzacc)

        !$ACC UPDATE DEVICE(p_phys_param%vmix_params%iwe_Tdis) IF(lzacc)

        !$ACC UPDATE DEVICE(sea_ice%draftave, sea_ice%concsum) IF(lzacc)
#endif

        start_timer(timer_scalar_prod_veloc,2)
        CALL calc_scalar_product_veloc_3d( patch_3d,  &
          & ocean_state(jg)%p_prog(nold(1))%vn,         &
          & ocean_state(jg)%p_diag,                     &
          & operators_coefficients,                     &
          & lacc=lzacc)
        stop_timer(timer_scalar_prod_veloc,2)

        !In case of a time-varying forcing:
        ! update_surface_flux or update_ocean_surface has changed p_prog(nold(1))%h, SST and SSS
        start_timer(timer_upd_flx,3)

        !! Updates velocity, tracer boundary condition
        !! Changes height based on ice etc
        !! Ice eqn sends back heat fluxes and volume fluxes
        !! Tracer relaxation and surface flux boundary conditions
        CALL update_ocean_surface_refactor( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, p_oce_sfc, &
            & current_time, operators_coefficients, ocean_state(jg)%p_prog(nold(1))%eta_c, &
            & ocean_state(jg)%p_prog(nold(1))%stretch_c, lacc=lzacc)

        stop_timer(timer_upd_flx,3)

        !------------------------------------------------------------------------
        !! Update stretch variables
        CALL update_zstar_variables( patch_3d, ocean_state(jg), operators_coefficients, &
          & ocean_state(jg)%p_prog(nold(1))%eta_c, &
          & ocean_state(jg)%p_prog(nold(1))%stretch_c, stretch_e, lacc=lzacc)

        !------------------------------------------------------------------------
        ! calculate in situ density here, as it may be used fotr the tides load
!        CALL calculate_density_zstar( patch_3d,                         &
!          & ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,:),      &
!          & ocean_state(jg)%p_prog(nold(1))%eta_c, &
!          & ocean_state(jg)%p_prog(nold(1))%stretch_c, &
!          & ocean_state(jg)%p_diag%rho(:,:,:), &
!          & lacc=lzacc )
!--------------------------------------------------------------------------
        !! zlev density is calculated on fixed depth and ignores any surface
        !! variation ; to simplify the implementation for zstar we do the same ; 
        !! calculate in situ density on fixed levels
        CALL calculate_density( patch_3d,                         &
          & ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,:),      &
          & ocean_state(jg)%p_diag%rho(:,:,:), lacc = lzacc)


        !--------------------------------------------------------------------------
        CALL create_pressure_bc_conditions(patch_3d,ocean_state(jg), p_as, sea_ice, current_time, lacc=lzacc)
        !------------------------------------------------------------------------


        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('on entry: h-old'           ,ocean_state(jg)%p_prog(nold(1))%h ,str_module,idt_src, &
          & patch_2d%cells%owned )
        CALL dbg_print('on entry: h-new'           ,ocean_state(jg)%p_prog(nnew(1))%h ,str_module,idt_src, &
          & patch_2d%cells%owned )
        CALL dbg_print('HydOce: ScaProdVel kin'    ,ocean_state(jg)%p_diag%kin        ,str_module,idt_src, &
          & patch_2d%cells%owned )
        CALL dbg_print('HydOce: ScaProdVel ptp_vn' ,ocean_state(jg)%p_diag%ptp_vn     ,str_module,idt_src, &
          & patch_2d%edges%owned )
        CALL dbg_print('HydOce: fu10'              ,p_as%fu10                         ,str_module,idt_src, &
          & in_subset=patch_2d%cells%owned)
        CALL dbg_print('HydOce: concsum'           ,sea_ice%concsum                   ,str_module,idt_src, &
          & in_subset=patch_2d%cells%owned)

        !---------------------------------------------------------------------
        !by_ogut: added p_oce_sfc
        CALL update_ho_params_zstar(patch_3d, ocean_state(jg), p_as%fu10, sea_ice%concsum, p_phys_param, operators_coefficients, &
                            & p_atm_f, p_oce_sfc, ocean_state(jg)%p_prog(nold(1))%eta_c, &
                            & ocean_state(jg)%p_prog(nold(1))%stretch_c, stretch_e, lacc=lzacc)

#ifdef _OPENACC
        !$ACC UPDATE SELF(ocean_state(jg)%p_diag%heatabs, ocean_state(jg)%p_diag%heatflux_rainevaprunoff) &
        !$ACC   SELF(ocean_state(jg)%p_diag%rsdoabsorb, ocean_state(jg)%p_diag%delta_ice) &
        !$ACC   SELF(ocean_state(jg)%p_diag%delta_snow, ocean_state(jg)%p_diag%delta_thetao) &
        !$ACC   SELF(ocean_state(jg)%p_diag%delta_so, ocean_state(jg)%p_diag%rho) &
        !$ACC   SELF(ocean_state(jg)%p_diag%p_vn, ocean_state(jg)%p_diag%ptp_vn, ocean_state(jg)%p_diag%kin) &
        !$ACC   SELF(ocean_state(jg)%p_diag%p_vn_dual, ocean_state(jg)%p_diag%u, ocean_state(jg)%p_diag%v) &
        !$ACC   SELF(ocean_state(jg)%p_diag%w_deriv, ocean_state(jg)%p_prog(nold(1))%tracer) &
        !$ACC   SELF(ocean_state(jg)%p_diag%Richardson_Number, ocean_state(jg)%p_diag%zgrad_rho) &
        !$ACC   SELF(ocean_state(jg)%p_aux%bc_total_top_potential) IF(lzacc)

        !$ACC UPDATE SELF(p_oce_sfc%HeatFlux_ShortWave, p_oce_sfc%HeatFlux_LongWave, p_oce_sfc%HeatFlux_Sensible) &
        !$ACC   SELF(p_oce_sfc%HeatFlux_Latent, p_oce_sfc%HeatFlux_Total, p_oce_sfc%FrshFlux_IceSalt) &
        !$ACC   SELF(p_oce_sfc%FrshFlux_VolumeIce, p_oce_sfc%FrshFlux_TotalIce, p_oce_sfc%FrshFlux_VolumeTotal) &
        !$ACC   SELF(p_oce_sfc%FrshFlux_TotalSalt, p_oce_sfc%top_dilution_coeff) &
        !$ACC   SELF(p_oce_sfc%TopBC_WindStress_u, p_oce_sfc%TopBC_WindStress_v, p_oce_sfc%TopBC_WindStress_cc) &
        !$ACC   SELF(p_oce_sfc%Wind_Speed_10m, p_oce_sfc%FrshFlux_Precipitation, p_oce_sfc%cellThicknessUnderIce) &
        !$ACC   SELF(p_oce_sfc%FrshFlux_Evaporation, p_oce_sfc%FrshFlux_Runoff, p_oce_sfc%FrshFlux_TotalOcean) &
        !$ACC   SELF(p_oce_sfc%FrshFlux_SnowFall, p_oce_sfc%TempFlux_Relax) &
        !$ACC   SELF(p_oce_sfc%HeatFlux_Relax, p_oce_sfc%SaltFlux_Relax, p_oce_sfc%FrshFlux_Relax) IF(lzacc)

        !$ACC UPDATE SELF(p_as%topBoundCond_windStress_u, p_as%topBoundCond_windStress_v, p_as%u, p_as%v) &
        !$ACC   SELF(p_as%tafo, p_as%ftdew, p_as%fu10, p_as%fclou, p_as%pao, p_as%fswr) &
        !$ACC   SELF(p_as%FrshFlux_Precipitation, p_as%FrshFlux_Runoff, p_as%data_surfRelax_Temp) IF(lzacc)

        !$ACC UPDATE SELF(sea_ice%hi, sea_ice%conc, sea_ice%hs, sea_ice%zunderice, sea_ice%u, sea_ice%v) &
        !$ACC   SELF(sea_ice%qtop, sea_ice%qbot, sea_ice%qbot_slow, sea_ice%vn_e, sea_ice%draftave) &
        !$ACC   SELF(sea_ice%totalsnowfall, sea_ice%zheatocei, sea_ice%heatocei, sea_ice%heatocew, sea_ice%concSum) &
        !$ACC   SELF(sea_ice%Tsurf, sea_ice%vol, sea_ice%vols, sea_ice%draft) &
        !$ACC   SELF(sea_ice%delhi, sea_ice%delhs, sea_ice%draftave_old, sea_ice%hiold, sea_ice%newice) IF(lzacc)

        !$ACC UPDATE SELF(p_phys_param%a_veloc_v, p_phys_param%a_tracer_v, p_phys_param%vmix_params%tke) &
        !$ACC   SELF(p_phys_param%vmix_params%u_stokes, p_phys_param%vmix_params%hlc) &
        !$ACC   SELF(p_phys_param%vmix_params%wlc, p_phys_param%vmix_params%tke_plc) &
        !$ACC   SELF(p_phys_param%vmix_params%vmix_dummy_1, p_phys_param%vmix_params%vmix_dummy_2) &
        !$ACC   SELF(p_phys_param%vmix_params%vmix_dummy_3, p_phys_param%vmix_params%tke_Tbpr) &
        !$ACC   SELF(p_phys_param%vmix_params%tke_Tspr, p_phys_param%vmix_params%tke_Tdif) &
        !$ACC   SELF(p_phys_param%vmix_params%tke_Tdis, p_phys_param%vmix_params%tke_Twin) &
        !$ACC   SELF(p_phys_param%vmix_params%tke_Tiwf, p_phys_param%vmix_params%tke_Tbck) &
        !$ACC   SELF(p_phys_param%vmix_params%tke_Ttot, p_phys_param%vmix_params%tke_Lmix) &
        !$ACC   SELF(p_phys_param%vmix_params%tke_Pr) IF(lzacc)

        !$ACC UPDATE SELF(p_atm_f%stress_x, p_atm_f%stress_y, p_atm_f%stress_xw, p_atm_f%stress_yw) &
        !$ACC   SELF(p_atm_f%albvisdir, p_atm_f%albvisdif, p_atm_f%albnirdir, p_atm_f%albnirdif) IF(lzacc)

        !$ACC END DATA

        lzacc = .FALSE.
        i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif

        !------------------------------------------------------------------------
        ! solve for new free surface
        start_timer(timer_solve_ab,1)
        CALL solve_free_surface_eq_zstar( patch_3d, ocean_state, p_ext_data,  &
          & p_oce_sfc , p_as, p_phys_param, operators_coefficients, solvercoeff_sp, &
          & jstep, ocean_state(jg)%p_prog(nold(1))%eta_c, ocean_state(jg)%p_prog(nold(1))%stretch_c, &
          & stretch_e, ocean_state(jg)%p_prog(nnew(1))%eta_c, ocean_state(jg)%p_prog(nnew(1))%stretch_c, &
          & lacc=lzacc)

        stop_timer(timer_solve_ab,1)

#ifdef _OPENACC
        i_am_accel_node = my_process_is_work()    ! Activate GPUs
        lzacc = .TRUE.

        !$ACC DATA COPYIN(ocean_state(jg)%p_prog(nnew(1))%eta_c, stretch_e) &
        !$ACC   COPYIN(ocean_state(jg)%p_prog(nold(1))%stretch_c, ocean_state(jg)%p_prog(nnew(1))%stretch_c) &
        !$ACC   COPYIN(ocean_nudge%data_3dimRelax_Temp, ocean_nudge%data_3dimRelax_Salt) &
        !$ACC   COPY(ocean_nudge%forc_3dimRelax_Temp, ocean_nudge%forc_3dimRelax_Salt) IF(lzacc)

        !$ACC UPDATE DEVICE(patch_3d%p_patch_2d(1)%edges%cell_idx, patch_3d%p_patch_2d(1)%edges%cell_blk) &
        !$ACC   DEVICE(patch_3d%p_patch_2d(1)%cells%edge_idx, patch_3d%p_patch_2d(1)%cells%edge_blk) &
        !$ACC   DEVICE(patch_3d%p_patch_2d(1)%edges%in_domain%start_block, patch_3d%p_patch_2d(1)%edges%in_domain%end_block) &
        !$ACC   DEVICE(patch_3d%p_patch_2d(1)%cells%in_domain%start_block, patch_3d%p_patch_2d(1)%cells%in_domain%end_block) &
        !$ACC   DEVICE(patch_3d%p_patch_1d(1)%dolic_c, patch_3d%p_patch_1d(1)%dolic_e) &
        !$ACC   DEVICE(patch_3d%p_patch_1d(1)%depth_CellInterface, patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c) &
        !$ACC   DEVICE(patch_3d%p_patch_1d(1)%prism_thick_e) IF(lzacc)

        !$ACC UPDATE DEVICE(ocean_state(jg)%p_prog(nold(1))%vn, ocean_state(jg)%p_diag%vn_pred) &
        !$ACC   DEVICE(ocean_state(jg)%p_aux%bc_bot_vn, ocean_state(jg)%p_aux%bc_top_vn) &
        !$ACC   DEVICE(ocean_state(jg)%p_aux%bc_top_u, ocean_state(jg)%p_aux%bc_top_v) &
        !$ACC   DEVICE(ocean_state(jg)%p_aux%g_n, ocean_state(jg)%p_aux%g_nimd) &
        !$ACC   DEVICE(ocean_state(jg)%p_aux%p_rhs_sfc_eq, ocean_state(jg)%p_diag%rho) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%press_grad, ocean_state(jg)%p_diag%press_hyd) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%w_deriv) IF(lzacc)

        !$ACC UPDATE DEVICE(operators_coefficients%grad_coeff, operators_coefficients%div_coeff) &
        !$ACC   DEVICE(operators_coefficients%edge2edge_viacell_coeff) IF(lzacc)
#endif

        !------------------------------------------------------------------------
        ! Step 4: calculate final normal velocity from predicted horizontal
        ! velocity vn_pred and updated surface height
        start_timer(timer_normal_veloc,4)
        CALL calc_normal_velocity_ab_zstar(patch_3d, ocean_state(jg), operators_coefficients, &
          & ocean_state(jg)%p_prog(nnew(1))%eta_c, lacc=lzacc)
        stop_timer(timer_normal_veloc,4)

        !------------------------------------------------------------------------
        ! Step 5: calculate vertical velocity and mass_flx_e from continuity equation under
        ! incompressiblity condition in the non-shallow-water case
        start_timer(timer_vert_veloc,4)
        CALL calc_vert_velocity_bottomup_zstar( patch_3d, ocean_state(jg),operators_coefficients, &
          & ocean_state(jg)%p_prog(nold(1))%stretch_c, stretch_e, lacc=lzacc)
        stop_timer(timer_vert_veloc,4)

        IF (idbg_mxmn >= 2 .OR. debug_check_level > 5) THEN
          CALL horizontal_mean(values=ocean_state(jg)%p_prog(nnew(1))%eta_c(:,:), &
            & weights=patch_2d%cells%area(:,:), &
            & in_subset=patch_2d%cells%owned, mean=mean_height)
          CALL debug_printValue(description="Mean Height", val=mean_height, detail_level=2)
        ENDIF
        IF (debug_check_level > 5 .AND. idbg_mxmn >= 2) THEN
          ! check difference from old_mean_height
          CALL debug_printValue(description="Old/New Mean Height", &
            & val=old_mean_height, value1=mean_height, detail_level=2)
          ! check if vertical and horizontal fluxes add to 0
          CALL horizontal_mean(values=ocean_state(jg)%p_diag%w, weights=patch_2d%cells%area(:,:), &
            & in_subset=patch_2d%cells%owned, mean=verticalMeanFlux, start_level=2, end_level=n_zlev)

          DO level=2, n_zlev-1
            CALL debug_printValue(description="Mean vertical flux at", val=REAL(level,wp),  &
              & value1=verticalMeanFlux(level), detail_level=2)
          ENDDO
        END IF

        !------------------------------------------------------------------------
        ! fill diffusion coefficients
        DO i = 1, ocean_state(jg)%p_prog(nold(1))%tracer_collection%no_of_tracers
          ocean_state(jg)%p_prog(nold(1))%tracer_collection%tracer(i)%hor_diffusion_coeff => p_phys_param%TracerDiffusion_coeff(:,:,:,i)
          ocean_state(jg)%p_prog(nold(1))%tracer_collection%tracer(i)%ver_diffusion_coeff => p_phys_param%a_tracer_v(:,:,:,i)
          !$ACC ENTER DATA COPYIN(ocean_state(jg)%p_prog(nold(1))%tracer_collection%tracer(i)%concentration) IF(lzacc)
          !$ACC ENTER DATA COPYIN(ocean_state(jg)%p_prog(nnew(1))%tracer_collection%tracer(i)%concentration) IF(lzacc)
        END DO

#ifdef _OPENACC
        IF (GMRedi_configuration /= Cartesian_Mixing ) THEN
           IF (lzacc) CALL finish(routine, 'OpenACC version with GMRedi currently not tested/validated')
        END IF
#endif

        CALL tracer_transport_zstar(patch_3d, ocean_state(jg), p_as, sea_ice, &
          & p_oce_sfc, p_phys_param, operators_coefficients, current_time, &
          & ocean_state(jg)%p_prog(nold(1))%stretch_c, stretch_e, ocean_state(jg)%p_prog(nnew(1))%stretch_c, &
          & lacc=lzacc)

        IF (use_age_tracer) THEN
          CALL calc_age_tracer(patch_3d, ocean_state(jg), sea_ice, lacc=lzacc)
        ENDIF

        !$ACC WAIT(1)
        DO i = 1, ocean_state(jg)%p_prog(nold(1))%tracer_collection%no_of_tracers
          !$ACC EXIT DATA COPYOUT(ocean_state(jg)%p_prog(nold(1))%tracer_collection%tracer(i)%concentration) IF(lzacc)
          !$ACC EXIT DATA COPYOUT(ocean_state(jg)%p_prog(nnew(1))%tracer_collection%tracer(i)%concentration) IF(lzacc)
        END DO
        !------------------------------------------------------------------------

!        IF (lcheck_salt_content) CALL check_total_salt_content_zstar(130, &
!          & ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,2), patch_2d, &
!          & ocean_state(jg)%p_prog(nnew(1))%stretch_c(:,:), &
!          & patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:), sea_ice, p_oce_sfc)

        !------------------------------------------------------------------------
        ! Optional : nudge temperature and salinity
        !! FIXME zstar: Not adapted to zstar
        IF (no_tracer>=1) THEN
          CALL nudge_ocean_tracers( patch_3d, ocean_state(jg), lacc=lzacc)
        ENDIF

        !------------------------------------------------------------------------
        ! perform accumulation for special variables
        ! FIXME zstar: Not adapted to zstar
        start_detail_timer(timer_extra20,5)
        IF (no_tracer>=1) THEN
          CALL calc_potential_density( patch_3d,                            &
            & ocean_state(jg)%p_prog(nold(1))%tracer,                       &
            & ocean_state(jg)%p_diag%rhopot, lacc=lzacc)

          ! calculate diagnostic barotropic stream function
          CALL calc_psi (patch_3d, ocean_state(jg)%p_diag%u(:,:,:),         &
            & patch_3D%p_patch_1d(1)%prism_thick_c(:,:,:),                  &
            & ocean_state(jg)%p_diag%u_vint, current_time, lacc=lzacc)

          CALL dbg_print('calc_psi: u_vint' ,ocean_state(jg)%p_diag%u_vint, str_module, 3, in_subset=patch_2d%cells%owned)
        ENDIF

#ifdef _OPENACC
        !$ACC UPDATE SELF(ocean_state(jg)%p_diag%rhopot, ocean_state(jg)%p_diag%u_vint) &
        !$ACC   SELF(ocean_state(jg)%p_prog(nnew(1))%tracer, ocean_state(jg)%p_prog(nnew(1))%vn) &
        !$ACC   SELF(ocean_state(jg)%p_diag%vn_time_weighted, ocean_state(jg)%p_diag%w) &
        !$ACC   SELF(ocean_state(jg)%p_diag%div_mass_flx_c, ocean_state(jg)%p_diag%mass_flx_e) IF(lzacc)

        !$ACC END DATA

        lzacc = .FALSE.
        i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif

        CALL calc_fast_oce_diagnostics( patch_2d, &
            & patch_3d, &
            & ocean_state(1), &
            & patch_3d%p_patch_1d(1)%dolic_c, &
            & patch_3d%p_patch_1d(1)%prism_thick_c, &
            & patch_3d%p_patch_1d(1)%zlev_m, &
            & ocean_state(jg)%p_diag, &
            & ocean_state(jg)%p_prog(nnew(1))%eta_c, &
            & ocean_state(jg)%p_prog(nnew(1))%vn, &
            & ocean_state(jg)%p_prog(nnew(1))%tracer, &
            & p_atm_f, &
            & p_oce_sfc, &
            & sea_ice, &
            & lacc=lzacc)

        stop_detail_timer(timer_extra20,5)

#ifdef _OPENACC
        i_am_accel_node = my_process_is_work()    ! Activate GPUs
        lzacc = .TRUE.

        !$ACC UPDATE SELF(ocean_state(jg)%p_diag%delta_ice, ocean_state(jg)%p_diag%delta_snow) &
        !$ACC   SELF(ocean_state(jg)%p_diag%delta_thetao, ocean_state(jg)%p_diag%delta_so) IF(lzacc)

        !$ACC UPDATE DEVICE(sea_ice%hi, sea_ice%conc, sea_ice%hs, sea_ice%zunderice) &
        !$ACC   DEVICE(sea_ice%u, sea_ice%v, sea_ice%qtop, sea_ice%qbot, sea_ice%qbot_slow) &
        !$ACC   DEVICE(sea_ice%vn_e, sea_ice%draftave, sea_ice%totalsnowfall) &
        !$ACC   DEVICE(sea_ice%zheatocei, sea_ice%heatocei, sea_ice%heatocew) IF(lzacc)

        !$ACC UPDATE DEVICE(sea_ice%concSum, sea_ice%Tsurf, sea_ice%vol, sea_ice%vols, sea_ice%draft) IF(lzacc)

        !$ACC UPDATE DEVICE(p_oce_sfc%FrshFlux_Precipitation, p_oce_sfc%FrshFlux_Evaporation) &
        !$ACC   DEVICE(p_oce_sfc%FrshFlux_Runoff, p_oce_sfc%FrshFlux_VolumeTotal, p_oce_sfc%FrshFlux_TotalOcean) &
        !$ACC   DEVICE(p_oce_sfc%HeatFlux_Total, p_oce_sfc%HeatFlux_ShortWave, p_oce_sfc%HeatFlux_Longwave) &
        !$ACC   DEVICE(p_oce_sfc%HeatFlux_Sensible, p_oce_sfc%HeatFlux_Latent, p_oce_sfc%Wind_Speed_10m) &
        !$ACC   DEVICE(p_oce_sfc%FrshFlux_VolumeIce, p_oce_sfc%sea_level_pressure) IF(lzacc)

        !$ACC UPDATE DEVICE(ocean_state(jg)%p_prog(nnew(1))%vn, ocean_state(jg)%p_prog(nnew(1))%tracer) IF(lzacc)

        !$ACC UPDATE DEVICE(ocean_state(jg)%p_diag%vn_pred, ocean_state(jg)%p_diag%vn_time_weighted) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%w, ocean_state(jg)%p_diag%rhopot, ocean_state(jg)%p_diag%rho) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%mass_flx_e, ocean_state(jg)%p_diag%rsdoabsorb) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%heatabs, ocean_state(jg)%p_diag%press_grad) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%press_hyd, ocean_state(jg)%p_diag%ptp_vn) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%kin, ocean_state(jg)%p_diag%vort) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%u, ocean_state(jg)%p_diag%v) &
        !$ACC   DEVICE(ocean_state(jg)%p_diag%veloc_adv_horz) IF(lzacc)

        !$ACC UPDATE DEVICE(ocean_state(jg)%p_aux%bc_top_u, ocean_state(jg)%p_aux%bc_top_v) &
        !$ACC   DEVICE(ocean_state(jg)%p_aux%bc_top_vn, ocean_state(jg)%p_aux%bc_bot_vn) &
        !$ACC   DEVICE(ocean_state(jg)%p_aux%g_n, ocean_state(jg)%p_aux%g_nimd) &
        !$ACC   DEVICE(ocean_state(jg)%p_aux%p_rhs_sfc_eq) IF(lzacc)

        !$ACC UPDATE DEVICE(p_atm_f%stress_x, p_atm_f%stress_y, p_atm_f%stress_xw, p_atm_f%stress_yw) IF(lzacc)

        !$ACC UPDATE DEVICE(p_atm_f%albvisdir, p_atm_f%albvisdif, p_atm_f%albnirdir, p_atm_f%albnirdif) IF(lzacc)

        !$ACC UPDATE DEVICE(p_phys_param%A_veloc_v, p_phys_param%a_tracer_v) IF(lzacc)
  
        !$ACC UPDATE DEVICE(p_phys_param%vmix_params%tke, p_phys_param%vmix_params%tke_Lmix) &
        !$ACC   DEVICE(p_phys_param%vmix_params%tke_Ttot, p_phys_param%vmix_params%tke_Tbck) &
        !$ACC   DEVICE(p_phys_param%vmix_params%tke_Tbpr, p_phys_param%vmix_params%tke_Tspr) &
        !$ACC   DEVICE(p_phys_param%vmix_params%tke_Tdif, p_phys_param%vmix_params%tke_Tdis) &
        !$ACC   DEVICE(p_phys_param%vmix_params%tke_Twin, p_phys_param%vmix_params%tke_Pr) &
        !$ACC   DEVICE(p_phys_param%vmix_params%vmix_dummy_1, p_phys_param%vmix_params%vmix_dummy_2) &
        !$ACC   DEVICE(p_phys_param%vmix_params%vmix_dummy_3) IF(lzacc)
#endif

        CALL update_statistics

        CALL output_ocean( patch_3d, ocean_state, &
          &                current_time,              &
          &                p_oce_sfc,             &
          &                sea_ice,                 &
          &                jstep, jstep0, lacc=lzacc)

        IF ( is_coupled_to_output() ) THEN
          IF (ltimer) CALL timer_start(timer_coupling)
          CALL output_coupling(patch_3d%wet_c)
          IF (ltimer) CALL timer_stop(timer_coupling)
        END IF

#ifdef _OPENACC
        lzacc = .FALSE.
        i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif

        ! send and receive coupling fluxes for ocean at the end of time stepping loop
        ! FIXME zstar: Does this make sense for zstar
        IF (iforc_oce == Coupled_FluxFromAtmo) THEN  !  14

          IF ( is_coupled_to_atmo() ) THEN
            IF (ltimer) CALL timer_start(timer_coupling)
            CALL couple_ocean_toatmo_fluxes( &
              patch_3D, ocean_state(jg), sea_ice, p_atm_f, p_as)
            IF (ltimer) CALL timer_stop(timer_coupling)
          END IF

          ! copy fluxes updated in coupling from p_atm_f into p_oce_sfc
          p_oce_sfc%FrshFlux_Precipitation = p_atm_f%FrshFlux_Precipitation
          p_oce_sfc%FrshFlux_Evaporation   = p_atm_f%FrshFlux_Evaporation
          p_oce_sfc%FrshFlux_SnowFall      = p_atm_f%FrshFlux_SnowFall
          p_oce_sfc%HeatFlux_Total         = p_atm_f%HeatFlux_Total
          p_oce_sfc%HeatFlux_ShortWave     = p_atm_f%HeatFlux_ShortWave
          p_oce_sfc%HeatFlux_Longwave      = p_atm_f%HeatFlux_Longwave
          p_oce_sfc%HeatFlux_Sensible      = p_atm_f%HeatFlux_Sensible
          p_oce_sfc%HeatFlux_Latent        = p_atm_f%HeatFlux_Latent
          p_oce_sfc%FrshFlux_Runoff        = p_atm_f%FrshFlux_Runoff

        ENDIF

        ! copy atmospheric wind speed of coupling from p_as%fu10 into forcing to be written by restart
        p_oce_sfc%Wind_Speed_10m(:,:) = p_as%fu10(:,:)
        p_oce_sfc%sea_level_pressure(:,:) = p_as%pao(:,:)


#ifdef _OPENACC
        i_am_accel_node = my_process_is_work()    ! Activate GPUs
        lzacc = .TRUE.

        !$ACC UPDATE DEVICE(p_oce_sfc%Wind_Speed_10m, p_oce_sfc%sea_level_pressure) IF(lzacc)

        lzacc = .FALSE.
        i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif

        start_detail_timer(timer_extra21,5)

        ! Shift time indices for the next loop
        ! this HAS to ge into the restart files, because the start with the following loop
        CALL update_time_indices(jg)
        !$ACC UPDATE DEVICE(nold, nnew)

        ! update intermediate timestepping variables for the tracers
        CALL update_time_g_n(ocean_state(jg))

        ! by_nils: layer diagnostic
        IF (use_layers) THEN
          CALL calc_layers(patch_3d, ocean_state(jg), p_oce_sfc, operators_coefficients, p_phys_param, ocean_state(jg)%p_prog(nold(1))%stretch_c, stretch_e)
        END IF

        ! check whether time has come for writing restart file
        IF (isCheckpoint()) THEN
          IF (.NOT. output_mode%l_none ) THEN
            !
            ! For multifile restart (restart_write_mode = "joint procs multifile")
            ! the domain flag must be set to .TRUE. in order to activate the domain,
            ! even though we have one currently in the ocean. Without this the
            ! processes won't write out their data into a the patch restart files.
            !
            patch_2d%ldom_active = .TRUE.
            !

#ifdef _OPENACC
            i_am_accel_node = my_process_is_work()    ! Activate GPUs
            lzacc = .TRUE.

            !$ACC UPDATE DEVICE(ocean_state(jg)%p_prog(nnew(1))%vn, ocean_state(jg)%p_prog(nnew(1))%tracer) IF(lzacc)
#endif
            IF (i_ice_dyn == 1) CALL ice_fem_update_vel_restart(patch_2d, sea_ice, lacc=lzacc) ! write FEM vel to restart or checkpoint file

            CALL restartDescriptor%updatePatch(patch_2d, &
                                              &opt_nice_class=1, &
                                              &opt_ocean_zlevels=n_zlev, &
                                              &opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
                                              &opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
            CALL restartDescriptor%writeRestart(current_time, jstep)
#ifdef _OPENACC
            lzacc = .FALSE.
            i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif
          END IF
        END IF

        stop_detail_timer(timer_extra21,5)

        IF (isEndOfThisRun()) THEN
          ! leave time loop
          RETURN
        END IF

        ! check cfl criterion
        IF (cfl_check) THEN
          CALL check_cfl_horizontal(ocean_state(jg)%p_prog(nnew(1))%vn, &
            & patch_2d%edges%inv_dual_edge_length, &
            & dtime, &
            & patch_2d%edges%ALL, &
            & cfl_threshold, &
            & ocean_state(jg)%p_diag%cfl_horz, &
            & cfl_stop_on_violation,&
            & cfl_write)
          CALL check_cfl_vertical(ocean_state(jg)%p_diag%w, &
            & patch_3d%p_patch_1d(1)%prism_center_dist_c, &
            & dtime, &
            & patch_2d%cells%ALL,&
            & cfl_threshold, &
            & ocean_state(jg)%p_diag%cfl_vert, &
            & cfl_stop_on_violation,&
            & cfl_write)
        END IF


    END SUBROUTINE ocean_time_step_zstar



    SUBROUTINE sed_only_time_step()

        CHARACTER(LEN = *), PARAMETER :: routine = 'mo_hydro_ocean_run:sed_only_time_step'

        ! fill transport state
        ocean_state(jg)%transport_state%patch_3d    => patch_3d

        ! optional memory loggin
        CALL memory_log_add

        jstep = jstep + 1
        ! update model date and time mtime based
        current_time = ocean_time_nextStep()

        IF (idbg_mxmn >= 1 .OR. debug_check_level > 3) THEN
          CALL datetimeToString(current_time, datestring)
          WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
          CALL message (TRIM(routine), message_text)
        END IF

        ! Add here the calls for Hamocc
        CALL ocean_to_hamocc_interface(ocean_state(jg), ocean_state(jg)%transport_state, &
          & p_oce_sfc, p_as, sea_ice, p_phys_param, operators_coefficients, current_time)

        CALL update_statistics

        CALL output_ocean( patch_3d, ocean_state, &
          &                current_time,              &
          &                p_oce_sfc,             &
          &                sea_ice,                 &
          &                jstep, jstep0, lacc=lzacc)

        ! check whether time has come for writing restart file
        IF (isCheckpoint()) THEN
          IF (.NOT. output_mode%l_none ) THEN
            !
            ! For multifile restart (restart_write_mode = "joint procs multifile")
            ! the domain flag must be set to .TRUE. in order to activate the domain,
            ! even though we have one currently in the ocean. Without this the
            ! processes won't write out their data into a the patch restart files.
            !
            patch_2d%ldom_active = .TRUE.
            !
            CALL restartDescriptor%updatePatch(patch_2d, &
                                              &opt_nice_class=1, &
                                              &opt_ocean_zlevels=n_zlev, &
                                              &opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
                                              &opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
            CALL restartDescriptor%writeRestart(current_time, jstep)
          END IF
        END IF

        IF (isEndOfThisRun()) THEN
          ! leave time loop
          RETURN
        END IF

    END SUBROUTINE sed_only_time_step

  END SUBROUTINE perform_ho_stepping
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE tracer_transport(patch_3d, ocean_state, p_as, sea_ice, p_oce_sfc, p_phys_param, &
      operators_coefficients, current_time, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: p_phys_param
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(datetime), POINTER, INTENT(in)              :: current_time
    LOGICAL, INTENT(in), OPTIONAL                    :: lacc

    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
    TYPE(t_ocean_transport_state), POINTER           :: transport_state
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: i, jg
    LOGICAL :: lzacc

    patch_2D              => patch_3D%p_patch_2D(1)
    all_cells             => patch_2D%cells%ALL
    old_tracer_collection => ocean_state%p_prog(nold(1))%tracer_collection
    new_tracer_collection => ocean_state%p_prog(nnew(1))%tracer_collection
    transport_state       => ocean_state%transport_state

    CALL set_acc_host_or_device(lzacc, lacc)

!     IF (no_tracer>=1) THEN
!       !calculate some information that is used for all tracers
!       IF (GMRedi_configuration==Cartesian_Mixing) THEN
!         CALL prepare_tracer_transport( patch_3d, &
!           & ocean_state, operators_coefficients)
!       ELSE
!         CALL prepare_tracer_transport_GMRedi( patch_3d, &
!           & ocean_state, p_phys_param,                  &
!           & operators_coefficients)
!       ENDIF
!     ENDIF

    !------------------------------------------------------------------------
    IF (no_tracer>=1) THEN

      ! fill transport_state
      transport_state%patch_3d    => patch_3d
      transport_state%w           => ocean_state%p_diag%w  ! w_time_weighted
      transport_state%mass_flux_e => ocean_state%p_diag%mass_flx_e
      transport_state%vn          => ocean_state%p_diag%vn_time_weighted

      IF (use_draftave_for_transport_h) THEN
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        transport_state%h_old     = ocean_state%p_prog(nold(1))%h - sea_ice%draftave
        transport_state%h_new     = ocean_state%p_prog(nnew(1))%h - sea_ice%draftave
        !$ACC END KERNELS
        !$ACC WAIT(1)
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        transport_state%h_old     = ocean_state%p_prog(nold(1))%h
        transport_state%h_new     = ocean_state%p_prog(nnew(1))%h
        !$ACC END KERNELS
        !$ACC WAIT(1)
      ENDIF

      ! fill boundary conditions
      old_tracer_collection%tracer(1)%top_bc => p_oce_sfc%TopBC_Temp_vdiff
      IF (no_tracer > 1) &
        old_tracer_collection%tracer(2)%top_bc => p_oce_sfc%TopBC_Salt_vdiff
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! transport tracers and diffuse them
      start_timer(timer_tracer_ab,1)

      IF (GMRedi_configuration == Cartesian_Mixing ) THEN
        CALL advect_ocean_tracers(old_tracer_collection, new_tracer_collection, &
          & transport_state, operators_coefficients, lacc=lzacc)
      ELSE
        CALL  advect_ocean_tracers_dev(old_tracer_collection, new_tracer_collection, &
          &  ocean_state, transport_state, p_phys_param, operators_coefficients)
      ENDIF

      stop_timer(timer_tracer_ab,1)
    !------------------------------------------------------------------------

      ! Call the biogeochemistry after transporting for GMRedi
      CALL ocean_to_hamocc_interface(ocean_state, transport_state, &
        & p_oce_sfc, p_as, sea_ice, p_phys_param, operators_coefficients, current_time)

    END IF

!     CALL dbg_print('Tr3:new adv', ocean_state%p_prog(nnew(1))%tracer(:,:,:,3),str_module,1, &
!       & patch_3d%p_patch_2d(1)%cells%owned )
!     CALL dbg_print('Tr20:new adv', ocean_state%p_prog(nnew(1))%tracer(:,:,:,20),str_module,1, &
!       & patch_3d%p_patch_2d(1)%cells%owned )

  END SUBROUTINE tracer_transport
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE tracer_transport_zstar(patch_3d, ocean_state, p_as, sea_ice, &
      & p_oce_sfc, p_phys_param, operators_coefficients, current_time, &
      & stretch_c, stretch_e, stretch_c_new, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: p_phys_param
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(datetime), POINTER, INTENT(in)              :: current_time
    REAL(wp), INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor
    REAL(wp), INTENT(IN) :: stretch_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL                    :: lacc

    TYPE(t_ocean_transport_state) , POINTER                   :: transport_state
    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
    TYPE(t_patch), POINTER                           :: patch_2D
    TYPE(t_subset_range), POINTER                    :: all_cells

    INTEGER :: i
    LOGICAL :: lzacc

    !------------------------------------------------------------------------
    !Tracer transport

    patch_2D              => patch_3D%p_patch_2D(1)
    all_cells             => patch_2D%cells%ALL
    old_tracer_collection => ocean_state%p_prog(nold(1))%tracer_collection
    new_tracer_collection => ocean_state%p_prog(nnew(1))%tracer_collection
    transport_state => ocean_state%transport_state

    CALL set_acc_host_or_device(lzacc, lacc)

    !------------------------------------------------------------------------
    IF (no_tracer>=1) THEN

      ! fill transport_state
      transport_state%patch_3d    => patch_3d
      transport_state%h_old       => ocean_state%p_prog(nold(1))%h
      transport_state%h_new       => ocean_state%p_prog(nnew(1))%h
      transport_state%w           => ocean_state%p_diag%w  ! w_time_weighted
      transport_state%mass_flux_e => ocean_state%p_diag%mass_flx_e
      transport_state%vn          => ocean_state%p_diag%vn_time_weighted


      ! fill boundary conditions
      old_tracer_collection%tracer(1)%top_bc => p_oce_sfc%TopBC_Temp_vdiff
      IF (no_tracer > 1) &
        old_tracer_collection%tracer(2)%top_bc => p_oce_sfc%TopBC_Salt_vdiff
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! FIXME zstar: cartesian diffusion not implemented
      ! transport tracers and diffuse them
      start_timer(timer_tracer_ab,1)

      IF (GMRedi_configuration==Cartesian_Mixing) THEN
        !! Note that zstar has no horizontal diffusion
        CALL advect_ocean_tracers_zstar(old_tracer_collection, new_tracer_collection, &
          & transport_state, operators_coefficients, stretch_e, stretch_c, stretch_c_new, lacc=lzacc)
      ELSE
        CALL  advect_ocean_tracers_GMRedi_zstar(old_tracer_collection, new_tracer_collection, &
          &  ocean_state, transport_state, p_phys_param, operators_coefficients, &
          &  stretch_c, stretch_e, stretch_c_new)
      ENDIF

      stop_timer(timer_tracer_ab,1)
      !------------------------------------------------------------------------
    END IF

    CALL ocean_to_hamocc_interface(ocean_state, transport_state, &
      & p_oce_sfc, p_as, sea_ice, p_phys_param, operators_coefficients, current_time)

  END SUBROUTINE tracer_transport_zstar
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE write_initial_ocean_timestep(patch_3d,ocean_state,p_oce_sfc,sea_ice, &
          & operators_coefficients, p_phys_param)

    TYPE(t_patch_3D), INTENT(IN) :: patch_3d
    TYPE(t_hydro_ocean_state), INTENT(INOUT)    :: ocean_state
    TYPE(t_ocean_surface) , INTENT(INOUT)       :: p_oce_sfc
    TYPE(t_sea_ice),          INTENT(INOUT)     :: sea_ice
    TYPE(t_operator_coeff),   INTENT(in)     :: operators_coefficients
    TYPE(t_ho_params), INTENT(IN), OPTIONAL     :: p_phys_param

    TYPE(t_patch), POINTER :: patch_2d

    patch_2d => patch_3d%p_patch_2d(1)

    ! in general nml output is writen based on the nnew status of the
    ! prognostics variables. Unfortunately, the initialization has to be written
    ! to the nold state. That's why the following manual copying is nec.
    ocean_state%p_prog(nnew(1))%h      = ocean_state%p_prog(nold(1))%h

    ocean_state%p_prog(nnew(1))%vn     = ocean_state%p_prog(nold(1))%vn

    CALL calc_scalar_product_veloc_3d( patch_3d,  ocean_state%p_prog(nnew(1))%vn,&
      & ocean_state%p_diag, operators_coefficients)
    ! CALL update_height_depdendent_variables( patch_3d, ocean_state, p_ext_data, operators_coefficients, solvercoeff_sp)
#ifdef _OPENACC
        i_am_accel_node = my_process_is_work()    ! Activate GPUs
#endif
    CALL update_statistics

#ifdef _OPENACC
    CALL write_name_list_output(jstep=0, lacc=i_am_accel_node)
#else
    CALL write_name_list_output(jstep=0)
#endif
#ifdef _OPENACC
        i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif

  END SUBROUTINE write_initial_ocean_timestep
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE fill_auxiliary_diagnostics(patch_3d, ocean_state)
    TYPE(t_patch_3D), INTENT(IN) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells, owned_cells
    INTEGER :: block,j,cell_idx,cells_startidx,cells_endidx, neigbor, vertex_idx, vertex_blk

    patch_2d => patch_3d%p_patch_2d(1)
    owned_cells => patch_2D%cells%owned

! !ICON_OMP_PARALLEL_DO PRIVATE(block,j,cell_idx,cells_startidx,cells_endidx, neigbor, vertex_idx, vertex_blk)
!     DO block = owned_cells%start_block, owned_cells%end_block
!       CALL get_index_range(owned_cells, block, cells_startidx, cells_endidx)
!       DO cell_idx = cells_startidx, cells_endidx
!         ocean_state%p_diag%vort_f_cells_50m(cell_idx,block) = 0.0_wp
!         DO neigbor=1,patch_2D%cells%max_connectivity
!           vertex_blk = patch_2d%cells%vertex_blk(cell_idx,block,neigbor)
!           vertex_idx = patch_2d%cells%vertex_idx(cell_idx,block,neigbor)
!           ocean_state%p_diag%vort_f_cells_50m(cell_idx,block) =  ocean_state%p_diag%vort_f_cells_50m(cell_idx,block) + &
!             & ocean_state%p_diag%vort(vertex_idx,15,vertex_blk)
!         ENDDO
!         ocean_state%p_diag%vort_f_cells_50m(cell_idx,block) = ocean_state%p_diag%vort_f_cells_50m(cell_idx,block) / &
!           & (REAL(patch_2D%cells%max_connectivity, wp) *  patch_2D%cells%f_c(cell_idx,block))
!
! #ifndef NAGFOR
!         IF (isnan(ocean_state%p_diag%vort_f_cells_50m(cell_idx,block))) &
!           & ocean_state%p_diag%vort_f_cells_50m(cell_idx,block) = 999.0_wp
! #endif
!
!       ENDDO
!     ENDDO
! !ICON_OMP_END_PARALLEL_DO
!
! !     ocean_state%p_diag%vort_50m(:,:) = ocean_state%p_diag%vort(:,15,:)
!     ocean_state%p_diag%T_50m(:,:)    = ocean_state%p_prog(nnew(1))%tracer_collection%tracer(1)%concentration(:,15,:)
!     ocean_state%p_diag%u_50m(:,:)    = ocean_state%p_diag%u(:,15,:)
!     ocean_state%p_diag%v_50m(:,:)    = ocean_state%p_diag%v(:,15,:)

  END SUBROUTINE fill_auxiliary_diagnostics
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE update_time_g_n(ocean_state)
    TYPE(t_hydro_ocean_state), INTENT(inout) :: ocean_state
    REAL(wp), POINTER ::  tmp(:,:,:)

    ! velocity
    ! just exchange the pointers
    ! ocean_state%p_aux%g_nm1 = ocean_state%p_aux%g_n
    ! ocean_state%p_aux%g_n   = 0.0_wp
    tmp => ocean_state%p_aux%g_n
    ocean_state%p_aux%g_n => ocean_state%p_aux%g_nm1
    ocean_state%p_aux%g_nm1 => tmp

  END SUBROUTINE update_time_g_n
END MODULE mo_hydro_ocean_run
