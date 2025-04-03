!> Initialization of the the phenology memory.
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_pheno_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,              ONLY: missval

  dsl4jsb_Use_processes PHENO_
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_memory(PHENO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: pheno_init

  TYPE t_pheno_init_vars
    REAL(wp), POINTER ::  &
      & lai_mon(:,:,:), &
      & fract_fpc_mon(:,:,:), &
      & fract_forest(:,:), &
      & fract_fpc_max(:,:)
  END TYPE t_pheno_init_vars

  TYPE(t_pheno_init_vars) :: pheno_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_init'

CONTAINS

  ! ======================================================================================================= !
  !> Run phenology init
  !>
  SUBROUTINE pheno_init(tile)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init'

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL pheno_read_init_vars(tile)
    END IF
    CALL pheno_init_bc(tile)
    CALL pheno_init_ic(tile)
    IF (tile%Is_last_process_tile(PHENO_)) THEN
      CALL pheno_finalize_init_vars()
    END IF

  END SUBROUTINE pheno_init

  SUBROUTINE pheno_finalize_init_vars

    DEALLOCATE( &
      & pheno_init_vars%lai_mon,       &
      & pheno_init_vars%fract_fpc_mon, &
      & pheno_init_vars%fract_forest,  &
      & pheno_init_vars%fract_fpc_max  &
      & )

  END SUBROUTINE pheno_finalize_init_vars

  SUBROUTINE pheno_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    dsl4jsb_Def_config(PHENO_)

    REAL(wp), POINTER ::  &
      & ptr_mon(:,:,:),   &
      & ptr_2D(:,:)
    TYPE(t_input_file) :: input_file

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ALLOCATE( &
      & pheno_init_vars%lai_mon      (nproma, nblks, 0:13), &
      & pheno_init_vars%fract_fpc_mon(nproma, nblks, 0:13), &
      & pheno_init_vars%fract_forest (nproma, nblks      ), &
      & pheno_init_vars%fract_fpc_max(nproma, nblks      )  &
      & )

    pheno_init_vars%lai_mon      (:,:,:) = missval
    pheno_init_vars%fract_fpc_mon(:,:,:) = missval
    pheno_init_vars%fract_forest (:,:)   = missval
    pheno_init_vars%fract_fpc_max(:,:)   = missval

    dsl4jsb_Get_config(PHENO_)

    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting phenology init vars from ' &
      &                          //TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(PHENO_)%bc_filename), model%grid_id)

    ! Leaf area index climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading lai_clim ...')
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='lai_clim')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    pheno_init_vars%lai_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    pheno_init_vars%lai_mon(:,:,0   ) = ptr_mon(:,:,12  )
    pheno_init_vars%lai_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)

    ! Foliage projected cover climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_fract ...')
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='veg_fract')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    pheno_init_vars%fract_fpc_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    pheno_init_vars%fract_fpc_mon(:,:,0   ) = ptr_mon(:,:,12  )
    pheno_init_vars%fract_fpc_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)

    ! Forest fraction
    ! @todo: only used for use_alb_veg_simple=.TRUE. and HYDRO_ init if l_organic=.TRUE.
    IF (debug_on()) CALL message(TRIM(routine), 'reading forest_fract ...')
    ptr_2D => input_file%Read_2d(      &
      & variable_name='forest_fract',  &
      & fill_array = pheno_init_vars%fract_forest)
    ptr_2D = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

    CALL input_file%Close()

    input_file = jsb_netcdf_open_input(TRIM(model%config%fract_filename), model%grid_id)     ! Temporary

    ! Maximum foliage projected cover
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_ratio_max ...')
    ptr_2D => input_file%Read_2d(      &
      & variable_name='veg_ratio_max', &
      & fill_array = pheno_init_vars%fract_fpc_max)
    ptr_2D = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

    CALL input_file%Close()

  END SUBROUTINE pheno_read_init_vars

  SUBROUTINE pheno_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    IF (debug_on()) CALL message(TRIM(routine), 'Setting boundary conditions of pheno memory for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    ! Leaf area index climatology
    IF (debug_on()) CALL message(TRIM(routine), 'setting lai_clim ...')

    dsl4jsb_memory(PHENO_)%lai_mon(:,:,:) = pheno_init_vars%lai_mon(:,:,:)
    !$ACC ENTER DATA ATTACH(dsl4jsb_memory(PHENO_)%lai_mon)
    !$ACC UPDATE DEVICE(dsl4jsb_memory(PHENO_)%lai_mon) ASYNC(1)

    ! Foliage projected cover climatology
    IF (debug_on()) CALL message(TRIM(routine), 'setting veg_fract ...')
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,:) = pheno_init_vars%fract_fpc_mon(:,:,:)
    !$ACC UPDATE DEVICE(dsl4jsb_memory(PHENO_)%fract_fpc_mon) ASYNC(1)

    ! Forest fraction
    ! TODO: only used for use_alb_veg_simple=.TRUE. and HYDRO_ init if l_organic=.TRUE.
    IF (debug_on()) CALL message(TRIM(routine), 'setting forest_fract ...')
    dsl4jsb_var2D_onDomain(PHENO_, fract_forest) = pheno_init_vars%fract_forest
    !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(PHENO_, fract_forest)) ASYNC(1)

    ! Maximum foliage projected cover
    IF (debug_on()) CALL message(TRIM(routine), 'setting veg_ratio_max ...')
    dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) = pheno_init_vars%fract_fpc_max
    !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max)) ASYNC(1)

    IF (tile%lcts(1)%lib_id /= 0)  THEN ! only if the present tile is a pft
      !JN 10.01.18: for now: set the maxLai for non forest types to the one from the Lctlib
      !             later one might think about allometric relation ships in other multi annual plants, too?
      IF ( (dsl4jsb_Config(PHENO_)%l_forestRegrowth) .AND. (.NOT. dsl4jsb_Lctlib_param(ForestFlag)) ) THEN
        dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom) = dsl4jsb_Lctlib_param(MaxLAI)
        !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom)) ASYNC(1)
      END IF

      !JN 10.01.18: with l_forestRegrowth the maxLAI can change over time for forest pfts
      IF ( (dsl4jsb_Config(PHENO_)%l_forestRegrowth) .AND. (dsl4jsb_Lctlib_param(ForestFlag)) ) THEN
        dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =   &
          & 1.0_wp - exp(-dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom) / dsl4jsb_Lctlib_param(clumpinessFactor))
      !@todo: do not rely on the lib_id but use a flag instead?
      ELSE
        dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =   &
          & 1.0_wp - exp(-dsl4jsb_Lctlib_param(MaxLAI) / dsl4jsb_Lctlib_param(clumpinessFactor))
      END IF
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction)) ASYNC(1)
    END IF

  END SUBROUTINE pheno_init_bc

  SUBROUTINE pheno_init_ic(tile)

    USE mo_pheno_process, ONLY: get_foliage_projected_cover
    USE mo_jsb_time,      ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Real2D_onDomain :: lai
    dsl4jsb_Real2D_onDomain :: fract_fpc
    dsl4jsb_Real2D_onDomain :: fract_fpc_max
    dsl4jsb_Real2D_onDomain :: MaxLai_allom
    REAL(wp), POINTER :: &
      & lai_mon(:,:,:), &
      & fract_fpc_mon(:,:,:)

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2
    INTEGER  :: nblks, ib, nc, ic

    REAL(wp) :: ClumpinessFactor_param, MaxLai_param

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init_ic'

    model => get_model(tile%owner_model_id)

    IF (tile%lcts(1)%lib_id /= 0) THEN
      ClumpinessFactor_param = dsl4jsb_Lctlib_param(ClumpinessFactor)
      MaxLai_param           = dsl4jsb_Lctlib_param(MaxLai)
    END IF

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of pheno memory (lai) for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    nblks = SIZE(dsl4jsb_var_ptr(PHENO_, lai), 2)
    nc    = SIZE(dsl4jsb_var_ptr(PHENO_, lai), 1)

    dsl4jsb_Get_var2D_onDomain(PHENO_, lai)
    dsl4jsb_Get_var2D_onDomain(PHENO_, fract_fpc)
    dsl4jsb_Get_var2D_onDomain(PHENO_, fract_fpc_max)
    dsl4jsb_Get_var2D_onDomain(PHENO_, maxLAI_allom)
    lai_mon       => dsl4jsb_memory(PHENO_)%lai_mon
    fract_fpc_mon => dsl4jsb_memory(PHENO_)%fract_fpc_mon

    ! Set initial leaf area index from climatology (read in pheno_init_bc)
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
    DO ib=1,nblks
      DO ic=1,nc
        lai(ic,ib) = wgt1 * lai_mon(ic,ib,nmw1) + wgt2 * lai_mon(ic,ib,nmw2)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    IF (model%config%l_compat401) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
      DO ib=1,nblks
        DO ic=1,nc
          fract_fpc(ic,ib) = get_foliage_projected_cover( &
            &                  fract_fpc_max(ic,ib), lai(ic,ib), 2._wp)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      IF (       TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'climatology' &
          & .OR. tile%lcts(1)%lib_id == 0) THEN         ! lib_id=0: general vegetation tile, no specific PFT
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
        DO ib=1,nblks
          DO ic=1,nc
            fract_fpc(ic,ib) = wgt1 * fract_fpc_mon(ic,ib,nmw1) + wgt2 * fract_fpc_mon(ic,ib,nmw2)
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        IF ( .NOT. dsl4jsb_Config(PHENO_)%l_forestRegrowth )  THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
          DO ib=1,nblks
            DO ic=1,nc
              lai(ic,ib) = MIN(lai(ic,ib), MaxLai_param)
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
          DO ib=1,nblks
            DO ic=1,nc
              lai(ic,ib) = MIN(lai(ic,ib), maxLAI_allom(ic,ib))
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
        DO ib=1,nblks
          DO ic=1,nc
            ! fract_fpc = get_foliage_projected_cover(                &
            !   & fract_fpc_max, lai, &
            !   & dsl4jsb_Lctlib_param(ClumpinessFactor))
            fract_fpc(ic,ib) = get_foliage_projected_cover(                &
              & fract_fpc_max(ic,ib), lai(ic,ib), ClumpinessFactor_param)
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END IF

  END SUBROUTINE pheno_init_ic

#endif
END MODULE mo_pheno_init
