!> Initialization of the the hd memory.
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
MODULE mo_hd_init
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_decomposition_tools, ONLY: get_local_index
  USE mo_impl_constants,      ONLY: SEA_BOUNDARY

  dsl4jsb_Use_processes HD_
  dsl4jsb_Use_config(HD_)
  dsl4jsb_Use_memory(HD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: hd_init_ic, hd_init_bc      ! initialization of initial and boundary conditions

  INTEGER, PARAMETER :: nneigh_max = 12 ! Maximum number of inflow neighbor cells

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_init'

CONTAINS

  SUBROUTINE hd_init_bc(tile)

    CLASS(t_jsb_tile_abstract),  INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER       :: model
    TYPE(t_jsb_grid),  POINTER       :: grid

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)

    INTEGER, POINTER ::         &
      & hd_mask_ptr(:,:),       &
      & cell_sea_land_mask(:,:)
    REAL(wp), POINTER ::         &
      & return_pointer   (:,:),  &   !< temporary pointer
      & cells_up_ptr (:,:,:)     !< temporary pointer

    TYPE(t_input_file)  :: input_file

    INTEGER :: n, blk, blks, blke, jc, jcs, jce, idx
    INTEGER :: nblks, nproma, nneigh

    CHARACTER(len=*),  PARAMETER     :: routine = modname//':hd_init_bc'

    model => get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    ! Get hd config
    dsl4jsb_Get_config(HD_)

    IF (.NOT. model%Is_process_enabled(HD_)) RETURN

    ! Get hd memory of the tile (only the root tile is used)
    dsl4jsb_Get_memory(HD_)

    SELECT CASE (TRIM(dsl4jsb_Config(HD_)%routing_scheme))
    CASE ('full')

      IF (debug_on()) CALL message(TRIM(routine), 'Initializing hd memory for tile '//TRIM(tile%name)//&
                                  ' from '//TRIM(dsl4jsb_Config(HD_)%bc_filename))

      ! Initialize boundary conditions
      !
      ! Number of base flow reservoirs
      dsl4jsb_var2D_onDomain(HD_, nres_baseflow) = 1._wp

      ! Maximum number of reservoirs
      ! If this is changed one has to change the definition of the three grids in mo_hd_config_class!!!
      dsl4jsb_memory(HD_)%nres_o_max = 1
      dsl4jsb_memory(HD_)%nres_b_max = 1
      dsl4jsb_memory(HD_)%nres_r_max = 5

      !CALL message(TRIM(routine), 'Initializing HD fractions for tile '//TRIM(tile%name)//&
      !  &                         ' from '//TRIM(dsl4jsb_Config(HD_)%fract_filneame))
      !input_file = jsb_netcdf_open_input('TRIM(dsl4jsb_Config(HD_)%fract_filneame', model%id)

      ! Initialize data from input file
      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HD_)%bc_filename), model%grid_id)

      ! HD model mask: -1 ocean, 0: ocean inflow cells, 1: land with outflow, 2: internal drainage
      hd_mask_ptr => input_file%Read_2d_int(        &
        & variable_name='MASK')
      IF (ASSOCIATED(hd_mask_ptr)) THEN
        dsl4jsb_var_ptr(HD_,hd_mask)(:,:) = REAL(hd_mask_ptr(:,:), wp)
        DEALLOCATE(hd_mask_ptr)
      END IF

      ! Retention constant of overland flow [day]
      return_pointer => input_file%Read_2d(         &
        & variable_name='ALF_K',                    &
        & fill_array = dsl4jsb_var_ptr(HD_,ret_overlflow))

      ! Number of overland flow reservoirs
      return_pointer => input_file%Read_2d(         &
        & variable_name='ALF_N',                    &
        & fill_array = dsl4jsb_var_ptr(HD_,nres_overlflow))

      ! Retention constant of river flow [day]
      return_pointer => input_file%Read_2d(         &
        & variable_name='ARF_K',                    &
        & fill_array = dsl4jsb_var_ptr(HD_,ret_riverflow))

      ! Number of river flow reservoirs
      return_pointer => input_file%Read_2d(         &
        & variable_name='ARF_N',                    &
        & fill_array = dsl4jsb_var_ptr(HD_,nres_riverflow))

      ! Retention constant for base flow [day]
      return_pointer => input_file%Read_2d(         &
        & variable_name='AGF_K',                    &
        & fill_array = dsl4jsb_var_ptr(HD_,ret_baseflow))

      ! Global index neighbouring grid cells upstream
      cells_up_ptr => input_file%Read_2d_extdim(  &
        & variable_name='CELLS_UP',                   &
        & start_extdim=1, end_extdim=nneigh_max,    &
        & extdim_name='nneigh')

      IF (dsl4jsb_config(HD_)%use_bifurcated_rivers) THEN
        ! Number of channel river splits into
        return_pointer => input_file%Read_2d(         &
          & variable_name='NSPLIT',                   &
          & fill_array = dsl4jsb_var_ptr(HD_,nsplit))
      ELSE
        dsl4jsb_var_ptr(HD_,nsplit)(:,:) = 1._wp
      END IF

      nproma = grid%nproma
      nblks  = grid%nblks

      nneigh = 0
      DO n = nneigh_max, 1, -1
        IF (ANY(cells_up_ptr(:,:,n) > 0)) THEN
          nneigh = n
          EXIT
        END IF
      END DO
      dsl4jsb_memory(HD_)%nneigh = nneigh

      ! find local idices of grid cells upstream
      ALLOCATE(dsl4jsb_memory(HD_)%nidx_upstream%ptr (nproma, nneigh, nblks))
      ALLOCATE(dsl4jsb_memory(HD_)%bidx_upstream%ptr (nproma, nneigh, nblks))
      !$ACC ENTER DATA &
      !$ACC   CREATE(dsl4jsb_memory(HD_)%nidx_upstream%ptr) &
      !$ACC   CREATE(dsl4jsb_memory(HD_)%bidx_upstream%ptr)
      dsl4jsb_memory(HD_)%nidx_upstream%ptr(:,:,:) = -1
      dsl4jsb_memory(HD_)%bidx_upstream%ptr(:,:,:) = -1

      blks = grid%get_blk_start()
      blke = grid%get_blk_end()
      DO blk = blks, blke
        jcs = grid%get_col_start(blk)
        jce = grid%get_col_end  (blk)
        DO n = 1, nneigh
          DO jc = jcs, jce
            IF (cells_up_ptr(jc,blk,n) > 0) THEN
              idx = get_local_index(grid%patch%cells%decomp_info%glb2loc_index, NINT(cells_up_ptr(jc,blk,n)))
              !print*, 'vg: p_pe: ', p_pe, ' n: ', n, ' blk: ', blk, ' FDIR: ', NINT(return_pointer(n,blk)), 'idx: ',idx
              IF (idx ==  0) THEN
                WRITE (message_text,*) "index of grid cell upstream (", NINT(cells_up_ptr(jc,blk,n)), ") not in global domain"
                CALL finish(routine, message_text)
              ELSE IF (idx == -1) THEN
                WRITE (message_text,*) "index of grid cell upstream (", NINT(cells_up_ptr(jc,blk,n)), ") not in local domain"
                CALL finish(routine, message_text)
              END IF
              dsl4jsb_memory(HD_)%bidx_upstream%ptr(jc,n,blk) = (idx-1)/nproma + 1
              dsl4jsb_memory(HD_)%nidx_upstream%ptr(jc,n,blk) = idx - (dsl4jsb_memory(HD_)%bidx_upstream%ptr(jc,n,blk)-1)*nproma
            END IF
          END DO
        END DO
      END DO

      !$ACC UPDATE &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HD_,hd_mask), dsl4jsb_var2D_onDomain(HD_,nsplit))               &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HD_,ret_overlflow), dsl4jsb_var2D_onDomain(HD_,nres_overlflow)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HD_,ret_riverflow), dsl4jsb_var2D_onDomain(HD_,nres_riverflow)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HD_,ret_baseflow), dsl4jsb_var2D_onDomain(HD_,nres_baseflow))   &
      !$ACC   DEVICE(dsl4jsb_memory(HD_)%bidx_upstream%ptr(:,:,:)) &
      !$ACC   DEVICE(dsl4jsb_memory(HD_)%nidx_upstream%ptr(:,:,:)) &
      !$ACC   ASYNC(1)

      !$ACC WAIT(1)

      CALL input_file%Close()
      IF (ASSOCIATED(cells_up_ptr)) DEALLOCATE(cells_up_ptr)

    CASE ('weighted_to_coast')

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HD_)%bc_filename), model%grid_id)

      ! coastal ocean cells
      cell_sea_land_mask => input_file%Read_2d_int( &
        & variable_name='cell_sea_land_mask')
      dsl4jsb_var2D_onDomain(HD_, coast_ocean) = &
        & MERGE(1._wp, 0._wp, cell_sea_land_mask(:,:) == SEA_BOUNDARY)
      dsl4jsb_var2D_onDomain(HD_, hd_mask) = &
        & MERGE(0._wp, 1._wp, cell_sea_land_mask(:,:) == SEA_BOUNDARY)
      DEALLOCATE(cell_sea_land_mask)

      CALL input_file%Close()

      !$ACC UPDATE &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HD_,hd_mask), dsl4jsb_var2D_onDomain(HD_,coast_ocean)) &
      !$ACC   ASYNC(1)

      !$ACC WAIT(1)

    CASE ('zero')

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HD_)%bc_filename), model%grid_id)

      ! coastal ocean cells
      cell_sea_land_mask => input_file%Read_2d_int( &
        & variable_name='cell_sea_land_mask')
      dsl4jsb_var2D_onDomain(HD_, hd_mask) = &
        & MERGE(0._wp, 1._wp, cell_sea_land_mask(:,:) == SEA_BOUNDARY)
      DEALLOCATE(cell_sea_land_mask)

      CALL input_file%Close()

      !$ACC UPDATE &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HD_,hd_mask)) &
      !$ACC   ASYNC(1)

      !$ACC WAIT(1)

    CASE ('none')

    END SELECT

    ! avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(return_pointer)) CONTINUE
    NULLIFY(return_pointer)

  END SUBROUTINE hd_init_bc

  SUBROUTINE hd_init_ic(tile)

    CLASS(t_jsb_tile_abstract),  INTENT(inout) :: tile

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)

    REAL(wp), POINTER :: &
      & tile_fract(:,:), &
      & ptr(:,:,:)            !< temporary pointer

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    INTEGER :: n, blk, blks, blke, jc, jcs, jce
    INTEGER :: nblks, nproma
    INTEGER :: idx

    CHARACTER(len=*), PARAMETER :: routine = modname//':hd_init_ic'

    model => get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    ! Get hd config
    dsl4jsb_Get_config(HD_)

    IF (.NOT. model%Is_process_enabled(HD_)) RETURN

    ! Get hd memory of the tile (only the root tile is used)
    dsl4jsb_Get_memory(HD_)

    IF (TRIM(dsl4jsb_Config(HD_)%routing_scheme) /= 'full') RETURN

    ! read initial conditions

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing hd memory from '//TRIM(dsl4jsb_Config(HD_)%ic_filename))

    ! Initialize boundary conditions
    !
    IF (dsl4jsb_config(HD_)%read_initial_reservoirs) THEN

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HD_)%ic_filename), model%grid_id)

      nproma = grid%nproma
      nblks  = grid%nblks
      blks = grid%get_blk_start()
      blke = grid%get_blk_end()

      ! Overlandflow reservoir
      ptr => input_file%Read_2d_extdim(               &
        & variable_name='FLFMEM',                     &
        & start_extdim=1,                             &
        & end_extdim= dsl4jsb_memory(HD_)%nres_o_max, &
        & extdim_name='oresnum')


      DO blk = blks, blke
         jcs = grid%get_col_start(blk)
         jce = grid%get_col_end  (blk)
         DO n = 1, dsl4jsb_memory(HD_)%nres_o_max
           DO jc = jcs, jce
             dsl4jsb_var_ptr(HD_,overlflow_res)(jc,n,blk) = ptr(jc,blk,n)
           END DO
         END DO
      END DO

      ! Baseflow reservoir
      ptr => input_file%Read_2d_extdim(               &
        & variable_name='FGMEM',                      &
        & start_extdim=1,                             &
        & end_extdim= dsl4jsb_memory(HD_)%nres_b_max, &
        & extdim_name='bresnum')

      DO blk = blks, blke
         jcs = grid%get_col_start(blk)
         jce = grid%get_col_end  (blk)
         DO n = 1, dsl4jsb_memory(HD_)%nres_b_max
           DO jc = jcs, jce
              dsl4jsb_var_ptr(HD_,baseflow_res)(jc,n,blk) = ptr(jc,blk,n)
            END DO
          END DO
      END DO

      ! Riverflow reservoir
      ptr => input_file%Read_2d_extdim(               &
        & variable_name='FRFMEM',                     &
        & start_extdim=1,                             &
        & end_extdim= dsl4jsb_memory(HD_)%nres_r_max, &
        & extdim_name='rresnum')

      DO blk = blks, blke
         jcs = grid%get_col_start(blk)
         jce = grid%get_col_end  (blk)
         DO n = 1, dsl4jsb_memory(HD_)%nres_r_max
           DO jc = jcs, jce
             dsl4jsb_var_ptr(HD_,riverflow_res)(jc,n,blk) = ptr(jc,blk,n)
           END DO
         END DO
      END DO

      CALL input_file%Close()

    ELSE

      ALLOCATE(tile_fract(grid%nproma, grid%nblks))
      CALL tile%Get_fraction(fract=tile_fract(:,:))
      dsl4jsb_var3D_onDomain(HD_,overlflow_res) = 1000._wp * SPREAD(tile_fract(:,:), DIM=2, NCOPIES=dsl4jsb_memory(HD_)%nres_o_max)
      dsl4jsb_var3D_onDomain(HD_,baseflow_res)  = 5000._wp * SPREAD(tile_fract(:,:), DIM=2, NCOPIES=dsl4jsb_memory(HD_)%nres_b_max)
      dsl4jsb_var3D_onDomain(HD_,riverflow_res) = 1000._wp * SPREAD(tile_fract(:,:), DIM=2, NCOPIES=dsl4jsb_memory(HD_)%nres_r_max)

    END IF

    !$ACC UPDATE &
    !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HD_,overlflow_res), dsl4jsb_var3D_onDomain(HD_,baseflow_res)) &
    !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HD_,riverflow_res)) &
    !$ACC   ASYNC(1)

    !$ACC WAIT(1)

  END SUBROUTINE hd_init_ic

#endif
END MODULE mo_hd_init
