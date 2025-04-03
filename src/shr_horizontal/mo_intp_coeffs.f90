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

! Contains the implementation of interpolation and reconstruction
! routines used by the shallow water model, including the RBF
! reconstruction routines.
  
#ifdef __xlC__
! @PROCESS nosmp
! @PROCESS NOOPTimize
! @PROCESS smp=noopt
@process noopt
#endif
!#ifdef __PGI
! !pgi$g opt=1
!#endif

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_intp_coeffs
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi2
  USE mo_exception,           ONLY: message, finish
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, min_rlcell_int, min_rledge_int
  USE mo_impl_constants_grf,  ONLY: grf_nudge_start_c, grf_nudge_start_e
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: lfeedback, grid_sphere_radius
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates, &
    &                               t_tangent_vectors
  USE mo_math_utilities,      ONLY: gc2cc, cc2gc, gnomonic_proj,               &
    & gvec2cvec, cvec2gvec,                      &
    & rotate_latlon, arc_length,                 &
    & plane_torus_closest_coordinates
  USE mo_dynamics_config,     ONLY: divavg_cntrwgt
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, sync_idx, global_max
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_interpol_config,     ONLY: nudge_zone_width, nudge_max_coeff, nudge_efold_width
  USE mo_grid_subset,         ONLY: get_index_range
  USE mo_lib_grid_geometry_info,  ONLY: planar_torus_geometry, sphere_geometry
  USE mo_grid_subset,         ONLY: get_index_range
  USE mo_fortran_tools,       ONLY: init

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_cellavg_wgt
  PUBLIC :: init_geo_factors
  PUBLIC :: complete_patchinfo
  PUBLIC :: init_tplane_e
  PUBLIC :: init_tplane_c
  PUBLIC :: init_nudgecoeffs
  PUBLIC :: tri_quadrature_pts


CONTAINS


  !-------------------------------------------------------------------------
  !! Computes the weighting coefficients for cell averaging with
  !! variable interpolation factors. Results are stored in ptr_patch%cells%avg_wgt
  !!
  SUBROUTINE init_cellavg_wgt( ptr_patch, ptr_int )
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    TYPE(t_int_state), INTENT(inout):: ptr_int

    INTEGER                     :: max_iter !max no. of iterations for forcing
                                            !mass conservation
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:init_cellavg_wgt'

    SELECT CASE(ptr_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)
      max_iter = 1
      CALL calculate_uniform_bilinear_cellavg_wgt( ptr_patch, ptr_int )
      CALL force_mass_conservation_to_cellavg_wgt( ptr_patch, ptr_int, max_iter)

    CASE (sphere_geometry)
      max_iter = 1000
      CALL calculate_bilinear_cellavg_wgt( ptr_patch, ptr_int )
      CALL force_mass_conservation_to_cellavg_wgt( ptr_patch, ptr_int, max_iter )

    CASE default
      CALL finish(method_name, "Undefined geometry type")

    END SELECT

  END SUBROUTINE init_cellavg_wgt
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes the weighting coefficients for cell averaging with
  !! variable interpolation factors. Results are stored in ptr_patch%cells%avg_wgt
  !!
  !! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
  !! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
  !! are not aliased into a checkerboard pattern between upward- and downward
  !! directed cells. The third condition is sum(w(i)) = 1., and the weight
  !! of the local point is 0.5 (see above).
  SUBROUTINE calculate_uniform_bilinear_cellavg_wgt( patch, interpolation_state )
    !  patch on which computation is performed
    TYPE(t_patch), TARGET, INTENT(inout) :: patch
    ! Interpolation state
    TYPE(t_int_state), INTENT(inout):: interpolation_state

    INTEGER :: cell_block, cell_index, start_index, end_index
    REAL(wp) :: local_weight, neigbor_weight

    !-----------------------------------------------------------------------
    ! Initial weighting factor of the local grid point
    local_weight = divavg_cntrwgt
    neigbor_weight = (1.0_wp - local_weight) / 3.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(cell_block, cell_index, start_index, end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_block = patch%cells%all%start_block, patch%cells%all%end_block
      CALL get_index_range(patch%cells%all, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        ! Simple for plane torus
        interpolation_state%c_bln_avg(cell_index,1,  cell_block) = local_weight
        interpolation_state%c_bln_avg(cell_index,2:4,cell_block) = neigbor_weight

      ENDDO !cell loop
  END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ! Note: no sync is required since the weights are calculated for all cells
  !       including halos
  END SUBROUTINE calculate_uniform_bilinear_cellavg_wgt
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes the weighting coefficients for cell averaging with
  !! variable interpolation factors. Results are stored in ptr_patch%cells%avg_wgt
  !!
  !! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
  !! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
  !! are not aliased into a checkerboard pattern between upward- and downward
  !! directed cells. The third condition is sum(w(i)) = 1., and the weight
  !! of the local point is 0.5 (see above).
  !!
  SUBROUTINE calculate_bilinear_cellavg_wgt( ptr_patch, ptr_int )
    !  patch on which computation is performed
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    ! Interpolation state
    TYPE(t_int_state), INTENT(inout):: ptr_int

    INTEGER :: jc, jb
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: ilc1, ibc1, ilc2, ibc2, ilc3, ibc3

    REAL(wp) :: xtemp,ytemp,wgt(3),xloc,yloc,x(3),y(3), &
      & pollat,pollon,wgt_loc

    REAL(wp) :: cell_area   ! area of triangular cell made up by 3 neighboring cell centers
    REAL(wp) :: mfac        ! 1 or -1, depending on the sign of yloc
                            ! takes care of correct sign of B-matrix
    !-----------------------------------------------------------------------

    ! Initial weighting factor of the local grid point
    wgt_loc = divavg_cntrwgt

    ! values for the blocking
    rl_start = 2
    rl_end = min_rlcell

    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

    ! Compute coefficients for bilinear interpolation
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,yloc,xloc,pollat,pollon,&
!$OMP            ilc1,ibc1,ilc2,ibc2,ilc3,ibc3,xtemp,ytemp,wgt,x,y,cell_area,mfac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx

        IF(.NOT.ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        yloc = ptr_patch%cells%center(jc,jb)%lat
        xloc = ptr_patch%cells%center(jc,jb)%lon

        ! Rotate local point into the equator for better accuracy of bilinear weights
        IF (yloc >= 0._wp) THEN
          pollat = yloc - pi2/4._wp
          mfac   = -1._wp
        ELSE
          pollat = yloc + pi2/4._wp
          mfac   = 1._wp
        ENDIF
        pollon = xloc

        CALL rotate_latlon( yloc, xloc, pollat, pollon )

        ! line and block indices of the neighbouring cells

        ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
        ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
        ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
        ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
        ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
        ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

        ! x and y are the zonal and meridional distances from the local
        ! cell point (ignoring the earth's radius, which drops out anyway)

        xtemp = ptr_patch%cells%center(ilc1,ibc1)%lon
        ytemp = ptr_patch%cells%center(ilc1,ibc1)%lat
        CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

        y(1)  = ytemp-yloc
        x(1)  = xtemp-xloc
        ! This is needed when the date line is crossed
        IF (x(1) >  3.5_wp) x(1) = x(1) - pi2
        IF (x(1) < -3.5_wp) x(1) = x(1) + pi2

        xtemp = ptr_patch%cells%center(ilc2,ibc2)%lon
        ytemp = ptr_patch%cells%center(ilc2,ibc2)%lat
        CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

        y(2)  = ytemp-yloc
        x(2)  = xtemp-xloc
        ! This is needed when the date line is crossed
        IF (x(2) >  3.5_wp) x(2) = x(2) - pi2
        IF (x(2) < -3.5_wp) x(2) = x(2) + pi2

        xtemp = ptr_patch%cells%center(ilc3,ibc3)%lon
        ytemp = ptr_patch%cells%center(ilc3,ibc3)%lat
        CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

        y(3)  = ytemp-yloc
        x(3)  = xtemp-xloc
        ! This is needed when the date line is crossed
        IF (x(3) >  3.5_wp) x(3) = x(3) - pi2
        IF (x(3) < -3.5_wp) x(3) = x(3) + pi2

        ! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
        ! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
        ! are not aliased into a checkerboard pattern between upward- and downward
        ! directed cells. The third condition is sum(w(i)) = 1., and the weight
        ! of the local point is 0.5 (see above). Analytical elimination yields...

        IF (ABS(x(2)-x(1)) > 1.e-11_wp .AND. ABS(y(3)-y(1)) > 1.e-11_wp ) THEN
          wgt(3) = 1._wp/( (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))/(x(2)-x(1)) ) * &
            & (1._wp-wgt_loc)*( -y(1) + x(1)*(y(2)-y(1))/(x(2)-x(1)) )
          wgt(2) = (-(1._wp-wgt_loc)*x(1) - wgt(3)*(x(3)-x(1)))/(x(2)-x(1))
          wgt(1) = 1._wp - wgt_loc - wgt(2) - wgt(3)
        ELSE
          wgt(2) = 1._wp/( (y(2)-y(1)) - (x(2)-x(1))*(y(3)-y(1))/(x(3)-x(1)) ) * &
            & (1._wp-wgt_loc)*( -y(1) + x(1)*(y(3)-y(1))/(x(3)-x(1)) )
          wgt(3) = (-(1._wp-wgt_loc)*x(1) - wgt(2)*(x(2)-x(1)))/(x(3)-x(1))
          wgt(1) = 1._wp - wgt_loc - wgt(2) - wgt(3)
        ENDIF

        ! Store results in ptr_patch%cells%avg_wgt
        ptr_int%c_bln_avg(jc,1,jb) = wgt_loc
        ptr_int%c_bln_avg(jc,2,jb) = wgt(1)
        ptr_int%c_bln_avg(jc,3,jb) = wgt(2)
        ptr_int%c_bln_avg(jc,4,jb) = wgt(3)


        ! B-matrix for cell based gradient (based on linear triangular finite element)
        !
        ! !!Attention!!: pollat IF-statement flips sign of (xi-xj), (yi-yj) terms
        ! That's why we need to multiply by -1 for yloc >= 0._wp
        !
        cell_area = 0.5_wp*grid_sphere_radius**2  &
          &       *((x(2)*y(3)-x(3)*y(2)) - (x(1)*y(3)-x(3)*y(1)) + (x(1)*y(2)-x(2)*y(1)))
        ! compute b-matrix
        ptr_int%gradc_bmat(jc,1,1,jb) = mfac*grid_sphere_radius*(y(2) - y(3))/(2._wp*cell_area)
        ptr_int%gradc_bmat(jc,1,2,jb) = mfac*grid_sphere_radius*(y(3) - y(1))/(2._wp*cell_area)
        ptr_int%gradc_bmat(jc,1,3,jb) = mfac*grid_sphere_radius*(y(1) - y(2))/(2._wp*cell_area)

        ptr_int%gradc_bmat(jc,2,1,jb) = mfac*grid_sphere_radius*(x(3) - x(2))/(2._wp*cell_area)
        ptr_int%gradc_bmat(jc,2,2,jb) = mfac*grid_sphere_radius*(x(1) - x(3))/(2._wp*cell_area)
        ptr_int%gradc_bmat(jc,2,3,jb) = mfac*grid_sphere_radius*(x(2) - x(1))/(2._wp*cell_area)

      ENDDO !cell loop

    END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%c_bln_avg)
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gradc_bmat(:,1,:,:))
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gradc_bmat(:,2,:,:))

  END SUBROUTINE calculate_bilinear_cellavg_wgt
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  ! The coefficients for bilinear interpolation are iteratively modified
  ! in order to obtain mass conservation.
  ! The criterion for conservation is that the three-point divergence
  ! calculated for any given grid point is used with a total factor of 1
  SUBROUTINE force_mass_conservation_to_cellavg_wgt( ptr_patch, ptr_int, niter )
    !  patch on which computation is performed
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    ! Interpolation state
    TYPE(t_int_state), INTENT(inout)     :: ptr_int
    ! max number of iterations
    INTEGER, INTENT(in)                  ::  niter

    INTEGER :: jc, je, jb
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: ilc1, ibc1, ilc2, ibc2, ilc3, ibc3, inb1, inb2, inb3, ie4, ie5
    INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3, ile4, ibe4
    INTEGER, DIMENSION(nproma) :: iie1, iie2, iie3, iie4
    REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2, z_nx3, z_nx4, z_nx5
    REAL(wp) :: checksum(nproma,ptr_patch%nblks_e)

    REAL(wp) :: relax_coeff
    INTEGER ::  iter

    REAL(wp) :: maxwgt_loc,minwgt_loc

    REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c)  :: wgt_loc_sum, resid

    INTEGER, DIMENSION(nproma,ptr_patch%nblks_c,3) :: inv_neighbor_id
    REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c,3) :: z_inv_neighbor_id

#ifdef DEBUG_COEFF
    REAL(wp) :: sum1
#endif
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:force_mass_conservation_to_cellavg_wgt'

    ! write(0,*) "force_mass_conservation_to_cellavg_wgt, processing ", ptr_patch%grid_filename
    !-----------------------------------------------------------------------
    ! The coefficients for bilinear interpolation are now iteratively modified
    ! in order to obtain mass conservation.
    ! The criterion for conservation is that the three-point divergence
    ! calculated for any given grid point is used with a total factor of 1
    ! Number of iterations for computation of bilinear weights
    i_nchdom   = MAX(1,ptr_patch%n_childdom)

    !niter = 1000 !now this value is passed as an arguement

    ! Relaxation coefficient for adaptation of local weight (empirically determined)
    relax_coeff = 0.46_wp

    ! Maximum/minimum  weighting factors of the local grid point
    maxwgt_loc = divavg_cntrwgt + 0.003_wp
    minwgt_loc = divavg_cntrwgt - 0.003_wp

    ! Initialization of the residuum  field
    resid(:,:) = 0._wp

    ! values for the blocking
    rl_start = 2
    rl_end = min_rlcell

    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! Compute inverse neighbor ID's
    ! The inverse neigbor ID of a neighbor cell (ilc1,ibc1) is the neighbor ID
    ! the local cell (jc,jb) has from the point of view of the neighbor cell

    inv_neighbor_id = 0
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,&
!$OMP ibc3) ICON_OMP_DEFAULT_SCHEDULE
!    DO jb = i_startblk, i_endblk
!
!      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
!        & i_startidx, i_endidx, rl_start, rl_end)
!
!      DO jc = i_startidx, i_endidx
     DO jb = ptr_patch%cells%in_domain%start_block, ptr_patch%cells%in_domain%end_block
       CALL get_index_range(ptr_patch%cells%in_domain, jb, i_startidx, i_endidx)
       DO jc = i_startidx, i_endidx

        IF(.NOT.ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        ! line and block indices of the neighbouring cells

        ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
        ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
        ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
        ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
        ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
        ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

        IF ( (ilc1>0) .AND. (ibc1>0) ) THEN
          IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,1) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc1,ibc1,1) == jb)) THEN
            inv_neighbor_id(jc,jb,1) = 1
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,2) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc1,ibc1,2) == jb)) THEN
            inv_neighbor_id(jc,jb,1) = 2
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,3) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc1,ibc1,3) == jb)) THEN
            inv_neighbor_id(jc,jb,1) = 3
          ELSE
            CALL finish(method_name, "Undefined inv_neighbor_id 1")
          ENDIF
        ENDIF
        IF ( (ilc2>0) .AND. (ibc2>0) ) THEN
          IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,1) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc2,ibc2,1) == jb)) THEN
            inv_neighbor_id(jc,jb,2)  = 1
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,2) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc2,ibc2,2) == jb)) THEN
            inv_neighbor_id(jc,jb,2)  = 2
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,3) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc2,ibc2,3) == jb)) THEN
            inv_neighbor_id(jc,jb,2)  = 3
          ELSE
            CALL finish(method_name, "Undefined inv_neighbor_id 2")
          ENDIF
        ENDIF
        IF ( (ilc3>0) .AND. (ibc3>0) ) THEN
          IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,1) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc3,ibc3,1) == jb)) THEN
            inv_neighbor_id(jc,jb,3)  = 1
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,2) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc3,ibc3,2) == jb)) THEN
            inv_neighbor_id(jc,jb,3)  = 2
          ELSE IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,3) == jc) .AND. &
            & (ptr_patch%cells%neighbor_blk(ilc3,ibc3,3) == jb)) THEN
            inv_neighbor_id(jc,jb,3)  = 3
          ELSE
            CALL finish(method_name, "Undefined inv_neighbor_id 3")
          ENDIF
        ENDIF

      ENDDO !cell loop

    END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    z_inv_neighbor_id = REAL(inv_neighbor_id,wp)
    CALL sync_patch_array(sync_c,ptr_patch,z_inv_neighbor_id(:,:,1))
    CALL sync_patch_array(sync_c,ptr_patch,z_inv_neighbor_id(:,:,2))
    CALL sync_patch_array(sync_c,ptr_patch,z_inv_neighbor_id(:,:,3))
    inv_neighbor_id = NINT(z_inv_neighbor_id)

    DO iter = 1, niter

      ! Compute sum of weighting coefficients with which
      ! each local divergence value is used
      ! Note: the summation needs to be split into 4 loops in order to
      ! allow for vectorization and parallelization
      wgt_loc_sum = 0._wp

      rl_start = 2
      rl_end = min_rlcell
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,ibc3,inb1,&
!$OMP inb2,inb3) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

          ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
          ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
          ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
          ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
          ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
          ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)
          inb1 = inv_neighbor_id(jc,jb,1) + 1
          inb2 = inv_neighbor_id(jc,jb,2) + 1
          inb3 = inv_neighbor_id(jc,jb,3) + 1

          wgt_loc_sum(jc,jb) = &
            & ptr_int%c_bln_avg(jc,1,jb)*ptr_patch%cells%area(jc,jb)          + &
            & ptr_int%c_bln_avg(ilc1,inb1,ibc1)*ptr_patch%cells%area(ilc1,ibc1)  + &
            & ptr_int%c_bln_avg(ilc2,inb2,ibc2)*ptr_patch%cells%area(ilc2,ibc2)  + &
            & ptr_int%c_bln_avg(ilc3,inb3,ibc3)*ptr_patch%cells%area(ilc3,ibc3)

        ENDDO !cell loop

      END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      rl_start = 3
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

          ! For mass conservation, wgt_loc_sum/area should be 1 for each cell
          ! The deviation therefrom is termed residuum here.

          resid(jc,jb) = wgt_loc_sum(jc,jb)/ptr_patch%cells%area(jc,jb)-1._wp

        ENDDO !cell loop

      END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      CALL sync_patch_array(sync_c,ptr_patch,resid)

      IF (iter < niter) THEN ! Apply iterative correction to weighting coefficients
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,&
!$OMP ilc3,ibc3) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx

            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

            ! line and block indices of the neighbouring cells

            ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
            ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
            ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
            ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
            ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
            ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

            ! Modify weighting coefficients

            ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - relax_coeff*resid(jc,jb)
            ptr_int%c_bln_avg(jc,2,jb) = ptr_int%c_bln_avg(jc,2,jb) - relax_coeff*resid(ilc1,ibc1)
            ptr_int%c_bln_avg(jc,3,jb) = ptr_int%c_bln_avg(jc,3,jb) - relax_coeff*resid(ilc2,ibc2)
            ptr_int%c_bln_avg(jc,4,jb) = ptr_int%c_bln_avg(jc,4,jb) - relax_coeff*resid(ilc3,ibc3)

            wgt_loc_sum(jc,jb) = SUM(ptr_int%c_bln_avg(jc,1:4,jb)) - 1._wp

            ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - 0.25_wp*wgt_loc_sum(jc,jb)
            ptr_int%c_bln_avg(jc,2,jb) = ptr_int%c_bln_avg(jc,2,jb) - 0.25_wp*wgt_loc_sum(jc,jb)
            ptr_int%c_bln_avg(jc,3,jb) = ptr_int%c_bln_avg(jc,3,jb) - 0.25_wp*wgt_loc_sum(jc,jb)
            ptr_int%c_bln_avg(jc,4,jb) = ptr_int%c_bln_avg(jc,4,jb) - 0.25_wp*wgt_loc_sum(jc,jb)

            ! To be safe: Avoid runaway of central weight
            ptr_int%c_bln_avg(jc,1,jb) = MAX(ptr_int%c_bln_avg(jc,1,jb),minwgt_loc)
            ptr_int%c_bln_avg(jc,1,jb) = MIN(ptr_int%c_bln_avg(jc,1,jb),maxwgt_loc)

          ENDDO !cell loop

        END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CALL sync_patch_array(sync_c,ptr_patch,ptr_int%c_bln_avg)

      ELSE ! In the last iteration, enforce the mass conservation condition
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx

            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

            ! Modify weighting coefficients

            ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - resid(jc,jb)

          ENDDO !cell loop

        END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CALL sync_patch_array(sync_c,ptr_patch,ptr_int%c_bln_avg)

        ! Compute coefficients needed to reconstruct averaged mass fluxes
        ! for approximately mass-consistent transport with divergence-averaging
        ! They can alternatively be used to average the velocity going into the divergence
        ! computation (without div averaging), yielding exact mass consistency but somewhat
        ! larger discretization errors for divergence

        rl_start = 4
        rl_end   = min_rledge
        i_startblk = ptr_patch%edges%start_blk(rl_start,1)
        i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,inb1,&
!$OMP inb2,inb3,ie4,ie5) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

            ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
            ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
            ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
            ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

            IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,1) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,1)) THEN

              inb1 = inv_neighbor_id(ilc1,ibc1,1)
              ie4  = MOD(inb1,  3)+1
              ie5  = MOD(inb1+1,3)+1

              ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb1+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,2,ibc1)/ptr_int%geofac_div(ilc2,inb1,ibc2)

              ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb1+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,3,ibc1)/ptr_int%geofac_div(ilc2,inb1,ibc2)

              ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,2,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,1,ibc1)

              ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,2,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,1,ibc1)


            ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,2) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,2)) THEN

              inb2 = inv_neighbor_id(ilc1,ibc1,2)
              ie4  = MOD(inb2  ,3)+1
              ie5  = MOD(inb2+1,3)+1

              ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb2+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,3,ibc1)/ptr_int%geofac_div(ilc2,inb2,ibc2)

              ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb2+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,1,ibc1)/ptr_int%geofac_div(ilc2,inb2,ibc2)

              ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,3,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,2,ibc1)

              ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,3,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,2,ibc1)


            ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,3) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,3)) THEN

              inb3 = inv_neighbor_id(ilc1,ibc1,3)
              ie4  = MOD(inb3  ,3)+1
              ie5  = MOD(inb3+1,3)+1

              ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb3+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,1,ibc1)/ptr_int%geofac_div(ilc2,inb3,ibc2)

              ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb3+1,ibc2) &
                & *ptr_int%geofac_div(ilc1,2,ibc1)/ptr_int%geofac_div(ilc2,inb3,ibc2)

              ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,4,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,3,ibc1)

              ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,4,ibc1) &
                & *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,3,ibc1)

            ENDIF

          ENDDO !edge loop

        END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%e_flx_avg)

        rl_start = 5
        i_startblk = ptr_patch%edges%start_blk(rl_start,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,inb1,inb2,inb3,ie4,ie5, &
!$OMP            ile1,ibe1,ile2,ibe2,ile3,ibe3,ile4,ibe4,iie1,iie2,iie3,iie4) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

            ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
            ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
            ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
            ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

            ile1 = ptr_patch%edges%quad_idx(je,jb,1)
            ibe1 = ptr_patch%edges%quad_blk(je,jb,1)
            ile2 = ptr_patch%edges%quad_idx(je,jb,2)
            ibe2 = ptr_patch%edges%quad_blk(je,jb,2)
            ile3 = ptr_patch%edges%quad_idx(je,jb,3)
            ibe3 = ptr_patch%edges%quad_blk(je,jb,3)
            ile4 = ptr_patch%edges%quad_idx(je,jb,4)
            ibe4 = ptr_patch%edges%quad_blk(je,jb,4)

            IF (ptr_patch%edges%cell_idx(ile1,ibe1,1) == ilc1 .AND. &
              & ptr_patch%edges%cell_blk(ile1,ibe1,1) == ibc1 ) THEN
              iie1(je) = 3
            ELSE IF (ptr_patch%edges%cell_idx(ile1,ibe1,2) == ilc1 .AND. &
              & ptr_patch%edges%cell_blk(ile1,ibe1,2) == ibc1 ) THEN
              iie1(je) = 5
            ENDIF
            IF (ptr_patch%edges%cell_idx(ile2,ibe2,1) == ilc1 .AND. &
              & ptr_patch%edges%cell_blk(ile2,ibe2,1) == ibc1 ) THEN
              iie2(je) = 2
            ELSE IF (ptr_patch%edges%cell_idx(ile2,ibe2,2) == ilc1 .AND. &
              & ptr_patch%edges%cell_blk(ile2,ibe2,2) == ibc1 ) THEN
              iie2(je) = 4
            ENDIF
            IF (ptr_patch%edges%cell_idx(ile3,ibe3,1) == ilc2 .AND. &
              & ptr_patch%edges%cell_blk(ile3,ibe3,1) == ibc2 ) THEN
              iie3(je) = 3
            ELSE IF (ptr_patch%edges%cell_idx(ile3,ibe3,2) == ilc2 .AND. &
              & ptr_patch%edges%cell_blk(ile3,ibe3,2) == ibc2 ) THEN
              iie3(je) = 5
            ENDIF
            IF (ptr_patch%edges%cell_idx(ile4,ibe4,1) == ilc2 .AND. &
              & ptr_patch%edges%cell_blk(ile4,ibe4,1) == ibc2 ) THEN
              iie4(je) = 2
            ELSE IF (ptr_patch%edges%cell_idx(ile4,ibe4,2) == ilc2 .AND. &
              & ptr_patch%edges%cell_blk(ile4,ibe4,2) == ibc2 ) THEN
              iie4(je) = 4
            ENDIF

            IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,1) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,1)) THEN

              inb1 = inv_neighbor_id(ilc1,ibc1,1)
              ie4  = MOD(inb1  ,3)+1
              ie5  = MOD(inb1+1,3)+1

              ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
                & ( ptr_int%geofac_div(ilc1,1,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
                & + ptr_int%geofac_div(ilc2,inb1,ibc2)*ptr_int%c_bln_avg(ilc1,2,ibc1)         &
                & - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,2,ibc1)     &
                & - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,3,ibc1) )   &
                & / ptr_int%geofac_div(ilc1,1,ibc1)   +                                       &
                & ( ptr_int%geofac_div(ilc2,inb1,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
                & + ptr_int%geofac_div(ilc1,1,ibc1)*ptr_int%c_bln_avg(ilc2,inb1+1,ibc2)       &
                & - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
                & - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
                & / ptr_int%geofac_div(ilc2,inb1,ibc2) )


            ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,2) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,2)) THEN

              inb2 = inv_neighbor_id(ilc1,ibc1,2)
              ie4  = MOD(inb2  ,3)+1
              ie5  = MOD(inb2+1,3)+1

              ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
                & ( ptr_int%geofac_div(ilc1,2,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
                & + ptr_int%geofac_div(ilc2,inb2,ibc2)*ptr_int%c_bln_avg(ilc1,3,ibc1)         &
                & - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,3,ibc1)     &
                & - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,1,ibc1) )   &
                & / ptr_int%geofac_div(ilc1,2,ibc1)   +                                       &
                & ( ptr_int%geofac_div(ilc2,inb2,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
                & + ptr_int%geofac_div(ilc1,2,ibc1)*ptr_int%c_bln_avg(ilc2,inb2+1,ibc2)       &
                & - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
                & - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
                & / ptr_int%geofac_div(ilc2,inb2,ibc2) )


            ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,3) .AND. &
              & jb == ptr_patch%cells%edge_blk(ilc1,ibc1,3)) THEN

              inb3 = inv_neighbor_id(ilc1,ibc1,3)
              ie4  = MOD(inb3  ,3)+1
              ie5  = MOD(inb3+1,3)+1

              ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
                & ( ptr_int%geofac_div(ilc1,3,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
                & + ptr_int%geofac_div(ilc2,inb3,ibc2)*ptr_int%c_bln_avg(ilc1,4,ibc1)         &
                & - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,1,ibc1)     &
                & - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,2,ibc1) )   &
                & / ptr_int%geofac_div(ilc1,3,ibc1)   +                                       &
                & ( ptr_int%geofac_div(ilc2,inb3,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
                & + ptr_int%geofac_div(ilc1,3,ibc1)*ptr_int%c_bln_avg(ilc2,inb3+1,ibc2)       &
                & - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
                & - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
                & / ptr_int%geofac_div(ilc2,inb3,ibc2) )

            ENDIF

          ENDDO !edge loop

        END DO !block loop
!$OMP END DO

        ! Finally, the weighting coefficients are scaled in order to
        ! yield the right result for a constant wind field

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ile1,ibe1,ile2,ibe2,ile3,ibe3,ile4,ibe4, &
!$OMP            z_nx1,z_nx2,z_nx3,z_nx4,z_nx5) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

            ile1 = ptr_patch%edges%quad_idx(je,jb,1)
            ibe1 = ptr_patch%edges%quad_blk(je,jb,1)
            ile2 = ptr_patch%edges%quad_idx(je,jb,2)
            ibe2 = ptr_patch%edges%quad_blk(je,jb,2)
            ile3 = ptr_patch%edges%quad_idx(je,jb,3)
            ibe3 = ptr_patch%edges%quad_blk(je,jb,3)
            ile4 = ptr_patch%edges%quad_idx(je,jb,4)
            ibe4 = ptr_patch%edges%quad_blk(je,jb,4)

            z_nx1(je,1:3) = ptr_patch%edges%primal_cart_normal(je,jb)%x(1:3)
            z_nx2(je,1:3) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(1:3)
            z_nx3(je,1:3) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(1:3)
            z_nx4(je,1:3) = ptr_patch%edges%primal_cart_normal(ile3,ibe3)%x(1:3)
            z_nx5(je,1:3) = ptr_patch%edges%primal_cart_normal(ile4,ibe4)%x(1:3)

            ! The sum of the coefficients - multiplied by the projection factors -
            ! is enforced to be 1 so that a constant vector field is processed correctly

            checksum(je,jb) = ptr_int%e_flx_avg(je,1,jb)                            &
              & + DOT_PRODUCT(z_nx1(je,1:3),z_nx2(je,1:3))*ptr_int%e_flx_avg(je,2,jb) &
              & + DOT_PRODUCT(z_nx1(je,1:3),z_nx3(je,1:3))*ptr_int%e_flx_avg(je,3,jb) &
              & + DOT_PRODUCT(z_nx1(je,1:3),z_nx4(je,1:3))*ptr_int%e_flx_avg(je,4,jb) &
              & + DOT_PRODUCT(z_nx1(je,1:3),z_nx5(je,1:3))*ptr_int%e_flx_avg(je,5,jb)

            ptr_int%e_flx_avg(je,1,jb) = ptr_int%e_flx_avg(je,1,jb)/checksum(je,jb)
            ptr_int%e_flx_avg(je,2,jb) = ptr_int%e_flx_avg(je,2,jb)/checksum(je,jb)
            ptr_int%e_flx_avg(je,3,jb) = ptr_int%e_flx_avg(je,3,jb)/checksum(je,jb)
            ptr_int%e_flx_avg(je,4,jb) = ptr_int%e_flx_avg(je,4,jb)/checksum(je,jb)
            ptr_int%e_flx_avg(je,5,jb) = ptr_int%e_flx_avg(je,5,jb)/checksum(je,jb)

          ENDDO !edge loop

        END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CALL sync_patch_array(sync_c,ptr_patch,ptr_int%c_bln_avg)
        CALL sync_patch_array(sync_e,ptr_patch,ptr_int%e_flx_avg)

      ENDIF ! end of last iteration
    ENDDO ! iteration loop

    ! Optional debug output for bilinear averaging coefficients
#ifdef DEBUG_COEFF

    rl_start = 2
    rl_end = min_rlcell
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

    sum1 = 0._wp
    wgt_loc_sum = 1._wp

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx

        sum1 = sum1 + resid(jc,jb)**2
        wgt_loc_sum(jc,jb) = SUM(ptr_int%c_bln_avg(jc,1:4,jb))

        WRITE(710+ptr_patch%id,'(2i5,5f12.6,e13.5)') jb,jc,ptr_int%c_bln_avg(jc,1:4,jb),&
          & wgt_loc_sum(jc,jb),resid(jc,jb)

      END DO
    END DO
    WRITE(710+ptr_patch%id,'(4e13.5)') MAXVAL(resid),SQRT(sum1/ptr_patch%n_patch_cells),&
      & MAXVAL(wgt_loc_sum)-1._wp,MINVAL(wgt_loc_sum)-1._wp
    CLOSE (710+ptr_patch%id)

    ! Debug output for mass flux averaging weights

    rl_start = 5
    rl_end = min_rledge
    i_startblk = ptr_patch%edges%start_blk(rl_start,1)
    i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = i_startidx, i_endidx

        WRITE(720+ptr_patch%id,'(2i5,6f12.6)') jb,je,ptr_int%e_flx_avg(je,1:5,jb),&
          & checksum(je,jb)

      END DO
    END DO
    CLOSE (720+ptr_patch%id)

#endif

  END SUBROUTINE force_mass_conservation_to_cellavg_wgt
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! Computes the coefficients for lateral boundary nudging needed for
  !! one-way nesting and the limited-area mode
  !! The nudging coefficients are defined via three namelist variables:
  !! nudge_max_coeff: Maximum relaxation coefficient in the cell row bordering to
  !! the boundary interpolation zone
  !! nudge_efold_width: e-folding width of exponential decay of coefficients
  !! (in units of grid cell rows)
  !! nudge_zone_width: Total width of nudging zone (in units of grid cell rows)
  !!
  SUBROUTINE init_nudgecoeffs( ptr_patch, ptr_int )
    !
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

    ! Interpolation state
    TYPE(t_int_state), INTENT(inout):: ptr_int
    !

    INTEGER :: jc, je, jb
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: max_rlval

    !-----------------------------------------------------------------------

    ! Check if required refin_ctrl information is available
    max_rlval = MAXVAL(ptr_patch%cells%refin_ctrl(:,:))
    max_rlval = NINT(global_max(REAL(max_rlval,wp)))

    ! write(0,*) nudge_zone_width, grf_nudge_start_c, max_rlval

    IF (max_rlval < nudge_zone_width+grf_nudge_start_c-1) THEN
      CALL finish('init_nudgecoeffs',&
        & 'bdy_indexing_depth in prepare_gridref must be at least nudge_zone_width+4')
    ENDIF

    i_nchdom   = MAX(1,ptr_patch%n_childdom)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! a) Nudging coefficients for cells
    i_startblk = ptr_patch%cells%start_blk(grf_nudge_start_c,1)
    i_endblk   = ptr_patch%cells%end_blk(min_rlcell_int,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, grf_nudge_start_c, min_rlcell_int)

      DO jc = i_startidx, i_endidx

        IF (ptr_patch%cells%refin_ctrl(jc,jb) > 0 .AND. &
          & ptr_patch%cells%refin_ctrl(jc,jb) <= nudge_zone_width+grf_nudge_start_c-1) THEN
          ptr_int%nudgecoeff_c(jc,jb) = &
            & nudge_max_coeff*EXP(-REAL(ptr_patch%cells%refin_ctrl(jc,jb)-grf_nudge_start_c,wp) / &
            & nudge_efold_width)
        ENDIF

      ENDDO !cell loop

    END DO !block loop
!$OMP END DO

    ! b) Nudging coefficients for edges
    i_startblk = ptr_patch%edges%start_blk(grf_nudge_start_e,1)
    i_endblk   = ptr_patch%edges%end_blk(min_rledge_int,i_nchdom)

    IF (ptr_patch%id > 1 .AND. lfeedback(ptr_patch%id)) THEN
      ! Use nudging coefficients optimized for velocity boundary diffusion
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, grf_nudge_start_e, min_rledge_int)

        DO je = i_startidx, i_endidx

          IF (ptr_patch%edges%refin_ctrl(je,jb) > 0 .AND. &
            & ptr_patch%edges%refin_ctrl(je,jb) <= grf_nudge_start_e+9) THEN
            ptr_int%nudgecoeff_e(je,jb) = nudge_max_coeff* &
              & EXP(-REAL(ptr_patch%edges%refin_ctrl(je,jb)-grf_nudge_start_e,wp) / 4._wp)
          ENDIF

        ENDDO !edge loop

      END DO !block loop
!$OMP END DO NOWAIT
    ELSE
      ! Use nudging coefficients from namelist
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, grf_nudge_start_e, min_rledge_int)

        DO je = i_startidx, i_endidx

          IF (ptr_patch%edges%refin_ctrl(je,jb) > 0 .AND. &
            & ptr_patch%edges%refin_ctrl(je,jb) <= 2*nudge_zone_width+grf_nudge_start_e-3) THEN
            ptr_int%nudgecoeff_e(je,jb) = &
              & nudge_max_coeff*EXP(-REAL(ptr_patch%edges%refin_ctrl(je,jb)-grf_nudge_start_e,wp) / &
              & (2._wp*nudge_efold_width))
          ENDIF

        ENDDO !edge loop

      END DO !block loop
!$OMP END DO NOWAIT
    ENDIF
!$OMP END PARALLEL

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%nudgecoeff_c)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_int%nudgecoeff_e)

  END SUBROUTINE init_nudgecoeffs
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Precomputes the geometrical factors used in the divergence, rotation.
  !!
  !! Precomputes the geometrical factors used in the divergence, rotation
  !! and nabla_2_scalar operators in order to improve computational efficiency
  !!
  SUBROUTINE init_geo_factors( ptr_patch, ptr_int )
    !
    IMPLICIT NONE
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

    ! Interpolation state
    TYPE(t_int_state), INTENT(inout):: ptr_int
    !

    INTEGER :: jc, jb, je, jv, je1
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: ile, ibe, ilc1, ibc1, ilc2, ibc2, ifac, ic, ilnc, ibnc, &
      & ile1, ibe1, ile2, ibe2, ile3, ibe3

    !-----------------------------------------------------------------------

    i_nchdom   = MAX(1,ptr_patch%n_childdom)


!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk,ifac)
    ! a) Geometrical factor for divergence
    rl_start = 1
    rl_end = min_rlcell

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jc,i_startidx,i_endidx,ile,ibe) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = 1, ptr_patch%geometry_info%cell_type
        DO jc = i_startidx, i_endidx

          IF (je > ptr_patch%cells%num_edges(jc,jb)) CYCLE ! relevant for hexagons
!             write(0,*) jc, jb, ":", ptr_patch%cells%num_edges(jc,jb)
!             STOP
!           ENDIF

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)

          ptr_int%geofac_div(jc,je,jb) = &
            & ptr_patch%edges%primal_edge_length(ile,ibe) * &
            & ptr_patch%cells%edge_orientation(jc,jb,je)  / &
            & ptr_patch%cells%area(jc,jb)

        ENDDO !cell loop
      ENDDO

    END DO !block loop
!$OMP END DO

    ! b) Geometrical factor for rotation
    rl_start = 2
    rl_end = min_rlvert

    ! Vorticity should have the right sign
    SELECT CASE (ptr_patch%geometry_info%cell_type)
    CASE (3)
      ifac = 1
    CASE (6)
      ifac = -1
    END SELECT
    ! values for the blocking
    i_startblk = ptr_patch%verts%start_blk(rl_start,1)
    i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jv,i_startidx,i_endidx,ile,ibe) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = 1, 9-ptr_patch%geometry_info%cell_type
        DO jv = i_startidx, i_endidx

          IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE

          IF (je > ptr_patch%verts%num_edges(jv,jb)) CYCLE

          ile = ptr_patch%verts%edge_idx(jv,jb,je)
          ibe = ptr_patch%verts%edge_blk(jv,jb,je)

          ptr_int%geofac_rot(jv,je,jb) =                &
            & ptr_patch%edges%dual_edge_length(ile,ibe) * &
            & ptr_patch%verts%edge_orientation(jv,jb,je)/ &
            & ptr_patch%verts%dual_area(jv,jb) * REAL(ifac,wp)

        ENDDO !vertex loop
      ENDDO

    END DO !block loop
!$OMP END DO

    ! c) Geometrical factor for nabla2_scalar
    rl_start = 2
    rl_end = min_rlcell

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilnc,ibnc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = 1, ptr_patch%geometry_info%cell_type
        DO jc = i_startidx, i_endidx

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)

          ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
          ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
          ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

          IF (jc == ilc1 .AND. jb == ibc1) THEN
            IF (ptr_patch%geometry_info%cell_type == 3) THEN
              ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) - &
                & ptr_int%geofac_div(jc,je,jb) /                            &
                & ptr_patch%edges%dual_edge_length(ile,ibe)
            ENDIF
          ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
            IF (ptr_patch%geometry_info%cell_type == 3) THEN
              ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) + &
                & ptr_int%geofac_div(jc,je,jb) /                            &
                & ptr_patch%edges%dual_edge_length(ile,ibe)
            ENDIF
          ENDIF
          DO ic = 1, ptr_patch%geometry_info%cell_type
            ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
            ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
            IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
              IF (ptr_patch%geometry_info%cell_type == 3) THEN
                ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) - &
                  & ptr_int%geofac_div(jc,je,jb) /                                  &
                  & ptr_patch%edges%dual_edge_length(ile,ibe)
              ENDIF
            ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
              IF (ptr_patch%geometry_info%cell_type == 3) THEN
                ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) + &
                  & ptr_int%geofac_div(jc,je,jb) /                                  &
                  & ptr_patch%edges%dual_edge_length(ile,ibe)
              ENDIF
            ENDIF
          ENDDO

          ! To ensure that dummy edges have a factor of 0:
          IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
            ptr_int%geofac_n2s(jc,je+1,jb) = 0._wp
          ENDIF

        ENDDO !cell loop
      ENDDO

    END DO !block loop
!$OMP END DO


    ! d) Geometrical factor for gradient of divergence (triangles only)

    ! sync does not work on patch 0 for some unknown reason. But we don't need this
    ! field on the radiation grid anyway, so let's just skip it
    IF (ptr_patch%geometry_info%cell_type == 3 .AND. ptr_patch%id >= 1) THEN

      rl_start = 2
      rl_end = min_rledge

      ! values for the blocking
      i_startblk = ptr_patch%edges%start_blk(rl_start,1)
      i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,je,je1,i_startidx,i_endidx,ile,ibe,ile1,ibe1,ile2,ibe2,ile3,ibe3,&
!$OMP            ilc1,ilc2,ibc1,ibc2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO je = i_startidx, i_endidx

          IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

          ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
          ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
          ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
          ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

          ! First consider edges of neighbor cell 2
          ile1 = ptr_patch%cells%edge_idx(ilc2,ibc2,1)
          ibe1 = ptr_patch%cells%edge_blk(ilc2,ibc2,1)
          ile2 = ptr_patch%cells%edge_idx(ilc2,ibc2,2)
          ibe2 = ptr_patch%cells%edge_blk(ilc2,ibc2,2)
          ile3 = ptr_patch%cells%edge_idx(ilc2,ibc2,3)
          ibe3 = ptr_patch%cells%edge_blk(ilc2,ibc2,3)

          IF (je == ile1 .AND. jb == ibe1) THEN
            ptr_int%geofac_grdiv(je,1,jb) = ptr_int%geofac_div(ilc2,1,ibc2)
          ELSE IF (je == ile2 .AND. jb == ibe2) THEN
            ptr_int%geofac_grdiv(je,1,jb) = ptr_int%geofac_div(ilc2,2,ibc2)
          ELSE IF (je == ile3 .AND. jb == ibe3) THEN
            ptr_int%geofac_grdiv(je,1,jb) = ptr_int%geofac_div(ilc2,3,ibc2)
          ENDIF

          ! Now consider edges of neighbor cell 1 and compute gradient
          ile1 = ptr_patch%cells%edge_idx(ilc1,ibc1,1)
          ibe1 = ptr_patch%cells%edge_blk(ilc1,ibc1,1)
          ile2 = ptr_patch%cells%edge_idx(ilc1,ibc1,2)
          ibe2 = ptr_patch%cells%edge_blk(ilc1,ibc1,2)
          ile3 = ptr_patch%cells%edge_idx(ilc1,ibc1,3)
          ibe3 = ptr_patch%cells%edge_blk(ilc1,ibc1,3)

          IF (je == ile1 .AND. jb == ibe1) THEN
            ptr_int%geofac_grdiv(je,1,jb) = (ptr_int%geofac_grdiv(je,1,jb) - &
              & ptr_int%geofac_div(ilc1,1,ibc1))*ptr_patch%edges%inv_dual_edge_length(je,jb)
          ELSE IF (je == ile2 .AND. jb == ibe2) THEN
            ptr_int%geofac_grdiv(je,1,jb) = (ptr_int%geofac_grdiv(je,1,jb) - &
              & ptr_int%geofac_div(ilc1,2,ibc1))*ptr_patch%edges%inv_dual_edge_length(je,jb)
          ELSE IF (je == ile3 .AND. jb == ibe3) THEN
            ptr_int%geofac_grdiv(je,1,jb) = (ptr_int%geofac_grdiv(je,1,jb) - &
              & ptr_int%geofac_div(ilc1,3,ibc1))*ptr_patch%edges%inv_dual_edge_length(je,jb)
          ENDIF
        ENDDO

        ! The quad edge indices are computed such that edges 1 and 2 border
        ! to neighbor cell 1 and edges 3 and 4 border to neighbor cell 2.
        ! Thus, splitting the following loop saves case discriminations
        DO je1 = 1, 2
          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

            ile = ptr_patch%edges%quad_idx(je,jb,je1)
            ibe = ptr_patch%edges%quad_blk(je,jb,je1)

            ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
            ibc1 = ptr_patch%edges%cell_blk(je,jb,1)

            ile1 = ptr_patch%cells%edge_idx(ilc1,ibc1,1)
            ibe1 = ptr_patch%cells%edge_blk(ilc1,ibc1,1)
            ile2 = ptr_patch%cells%edge_idx(ilc1,ibc1,2)
            ibe2 = ptr_patch%cells%edge_blk(ilc1,ibc1,2)
            ile3 = ptr_patch%cells%edge_idx(ilc1,ibc1,3)
            ibe3 = ptr_patch%cells%edge_blk(ilc1,ibc1,3)

            IF (ile == ile1 .AND. ibe == ibe1) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = - &
                & ptr_int%geofac_div(ilc1,1,ibc1)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ELSE IF (ile == ile2 .AND. ibe == ibe2) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = - &
                & ptr_int%geofac_div(ilc1,2,ibc1)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ELSE IF (ile == ile3 .AND. ibe == ibe3) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = - &
                & ptr_int%geofac_div(ilc1,3,ibc1)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ENDIF

          ENDDO !edge loop
        ENDDO

        DO je1 = 3, 4
          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

            ile = ptr_patch%edges%quad_idx(je,jb,je1)
            ibe = ptr_patch%edges%quad_blk(je,jb,je1)

            ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
            ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

            ile1 = ptr_patch%cells%edge_idx(ilc2,ibc2,1)
            ibe1 = ptr_patch%cells%edge_blk(ilc2,ibc2,1)
            ile2 = ptr_patch%cells%edge_idx(ilc2,ibc2,2)
            ibe2 = ptr_patch%cells%edge_blk(ilc2,ibc2,2)
            ile3 = ptr_patch%cells%edge_idx(ilc2,ibc2,3)
            ibe3 = ptr_patch%cells%edge_blk(ilc2,ibc2,3)

            IF (ile == ile1 .AND. ibe == ibe1) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = &
                & ptr_int%geofac_div(ilc2,1,ibc2)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ELSE IF (ile == ile2 .AND. ibe == ibe2) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = &
                & ptr_int%geofac_div(ilc2,2,ibc2)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ELSE IF (ile == ile3 .AND. ibe == ibe3) THEN
              ptr_int%geofac_grdiv(je,je1+1,jb) = &
                & ptr_int%geofac_div(ilc2,3,ibc2)*ptr_patch%edges%inv_dual_edge_length(je,jb)
            ENDIF

          ENDDO !edge loop
        ENDDO

      END DO !block loop
!$OMP END DO

    ENDIF


    ! e) Geometrical factor for Green-Gauss gradient
    rl_start = 2
    rl_end = min_rlcell

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilnc,ibnc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = 1, ptr_patch%geometry_info%cell_type
        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)

          ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
          ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
          ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

          IF (jc == ilc1 .AND. jb == ibc1) THEN
            ptr_int%geofac_grg(jc,1,jb,1) = ptr_int%geofac_grg(jc,1,jb,1) +      &
              & ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
              & ptr_int%c_lin_e(ile,1,ibe)
            ptr_int%geofac_grg(jc,1,jb,2) = ptr_int%geofac_grg(jc,1,jb,2) +      &
              & ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
              & ptr_int%c_lin_e(ile,1,ibe)
          ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
            ptr_int%geofac_grg(jc,1,jb,1) = ptr_int%geofac_grg(jc,1,jb,1) +      &
              & ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
              & ptr_int%c_lin_e(ile,2,ibe)
            ptr_int%geofac_grg(jc,1,jb,2) = ptr_int%geofac_grg(jc,1,jb,2) +      &
              & ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
              & ptr_int%c_lin_e(ile,2,ibe)
          ENDIF
          DO ic = 1, ptr_patch%geometry_info%cell_type
            ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
            ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
            IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
              ptr_int%geofac_grg(jc,ic+1,jb,1) = ptr_int%geofac_grg(jc,ic+1,jb,1)+ &
                & ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
                & ptr_int%c_lin_e(ile,1,ibe)
              ptr_int%geofac_grg(jc,ic+1,jb,2) = ptr_int%geofac_grg(jc,ic+1,jb,2)+ &
                & ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
                & ptr_int%c_lin_e(ile,1,ibe)
            ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
              ptr_int%geofac_grg(jc,ic+1,jb,1) = ptr_int%geofac_grg(jc,ic+1,jb,1)+ &
                & ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
                & ptr_int%c_lin_e(ile,2,ibe)
              ptr_int%geofac_grg(jc,ic+1,jb,2) = ptr_int%geofac_grg(jc,ic+1,jb,2)+ &
                & ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
                & ptr_int%c_lin_e(ile,2,ibe)
            ENDIF
          ENDDO

          ! To ensure that dummy edges have a factor of 0:
          IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
            ptr_int%geofac_grg(jc,je+1,jb,1:2) = 0._wp
          ENDIF

        ENDDO !cell loop
      ENDDO

    END DO !block loop
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%geofac_div)
    CALL sync_patch_array(sync_v,ptr_patch,ptr_int%geofac_rot)
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%geofac_n2s)

    IF (ptr_patch%geometry_info%cell_type == 3) THEN
      IF (ptr_patch%id >= 1) CALL sync_patch_array(sync_e,ptr_patch,ptr_int%geofac_grdiv)
    ENDIF

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%geofac_grg(:,:,:,1))
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%geofac_grg(:,:,:,2))

  END SUBROUTINE init_geo_factors
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! Completes the computation of geometric information for a given patch, such as
  !! - the control volume associated to an edge
  !! - the local orientation of the edge primal normal and dual normal at the
  !!   location of the cell centers and vertices.
  !! - the inverse primal and dual edge lengths
  !!
  !! edges% area_edge
  !!        inv_primal_edge_length
  !!        inv_dual_edge_length
  !!        vertex_idx/blk
  !!        inv_vert_vert_length
  !!        primal_normal_cell
  !!        dual_normal_cell
  !!        primal_normal_vert
  !!        dual_normal_vert
  !!
  !! Moreover, the local orientation of the edge primal normal at the location
  !! of the cell centers in a cell-based data structure is computed and stored
  !! in the interpolation state together with some distance fields.
  !!
  !! ptr_int% primal_normal_ec
  !!          edge_cell_length
  !!          cell_vert_dist
  SUBROUTINE complete_patchinfo( ptr_patch, ptr_int )
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch),     INTENT(inout) :: ptr_patch

    ! Interpolation state
    TYPE(t_int_state), INTENT(inout) :: ptr_int
    !

    INTEGER :: jb, je, jc
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

    INTEGER :: ilc1, ibc1, ilv1, ibv1, ilc2, ibc2, ilv2, ibv2, &
      & ilv3, ibv3, ilv4, ibv4, ile1, ibe1

    REAL(wp) :: z_nu, z_nv, z_lon, z_lat, z_nx1(3), z_nx2(3), z_norm

    TYPE(t_cartesian_coordinates) :: cc_cell, cc_ev3, cc_ev4, cc_v1, cc_v2, cc_v3, &
      & cc_dis1, cc_dis2, cc_dis3

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:complete_patchinfo'

    !-----------------------------------------------------------------------


!$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
    rl_start = 1
    rl_end = min_rledge

    ! values for the blocking
    i_startblk = ptr_patch%edges%start_block(rl_start)
    i_endblk   = ptr_patch%edges%end_block(rl_end)
    !
    ! The fields for the inverse primal and dual edge lengths are
    ! initialized here.
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je =  i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        ! compute the control volume associated to each edge.
        ! It is defined as the quadrilateral whose edges are the primal edge
        ! and the associated dual edge
        ptr_patch%edges%area_edge(je,jb) =  &
          &    ptr_patch%edges%primal_edge_length(je,jb)  &
          &  * ptr_patch%edges%dual_edge_length(je,jb)

        ! compute inverse primal edge length
        ! (dual follows below in the rl_start=2 section)
        ptr_patch%edges%inv_primal_edge_length(je,jb) = &
          & 1._wp/ptr_patch%edges%primal_edge_length(je,jb)

      ENDDO

    END DO !block loop
!$OMP END DO

    rl_start = 2
    rl_end = min_rledge

    ! Second step: computed projected orientation vectors and related information
    i_startblk = ptr_patch%edges%start_block(rl_start)
    ! i_endblk   = ptr_patch%edges%end_block(rl_end)
    i_endblk   = ptr_patch%nblks_e

    ! Initialization of lateral boundary points
    IF (ptr_patch%id > 1) THEN
      CALL init(ptr_patch%edges%inv_dual_edge_length(:,1:i_startblk), lacc=.FALSE.)
      CALL init(ptr_patch%edges%vertex_idx(:,1:i_startblk,3), lacc=.FALSE.)
      CALL init(ptr_patch%edges%vertex_idx(:,1:i_startblk,4), lacc=.FALSE.)
      CALL init(ptr_patch%edges%vertex_blk(:,1:i_startblk,3), lacc=.FALSE.)
      CALL init(ptr_patch%edges%vertex_blk(:,1:i_startblk,4), lacc=.FALSE.)
      CALL init(ptr_patch%edges%inv_vert_vert_length(:,1:i_startblk), lacc=.FALSE.)
      CALL init(ptr_patch%edges%primal_normal_cell(:,1:i_startblk,:)%v1, lacc=.FALSE.)
      CALL init(ptr_patch%edges%dual_normal_cell  (:,1:i_startblk,:)%v1, lacc=.FALSE.)
      CALL init(ptr_patch%edges%primal_normal_vert(:,1:i_startblk,:)%v1, lacc=.FALSE.)
      CALL init(ptr_patch%edges%dual_normal_vert  (:,1:i_startblk,:)%v1, lacc=.FALSE.)
      CALL init(ptr_patch%edges%primal_normal_cell(:,1:i_startblk,:)%v2, lacc=.FALSE.)
      CALL init(ptr_patch%edges%dual_normal_cell  (:,1:i_startblk,:)%v2, lacc=.FALSE.)
      CALL init(ptr_patch%edges%primal_normal_vert(:,1:i_startblk,:)%v2, lacc=.FALSE.)
      CALL init(ptr_patch%edges%dual_normal_vert  (:,1:i_startblk,:)%v2, lacc=.FALSE.)
!$OMP BARRIER
    ENDIF
    !
    ! loop through all patch edges
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,ilc1,ibc1,ilv1,ibv1,ilc2,ibc2,ilv2, &
!$OMP            ibv2,ilv3,ibv3,ilv4,ibv4,z_nu,z_nv,z_lon,z_lat,z_nx1,z_nx2,   &
!$OMP            cc_ev3,cc_ev4,z_norm) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je =  i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        ! compute inverse dual edge length (undefined for refin_ctrl=1)

        ptr_patch%edges%inv_dual_edge_length(je,jb) = &
          & 1._wp/ptr_patch%edges%dual_edge_length(je,jb)

        ! compute edge-vertex indices (and blocks) 3 and 4, which
        ! are the outer vertices of cells 1 and 2, respectively,
        ! and the inverse length bewtween vertices 3 and 4

        ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
        ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
        ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
        ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

        ilv1 = ptr_patch%edges%vertex_idx(je,jb,1)
        ibv1 = ptr_patch%edges%vertex_blk(je,jb,1)
        ilv2 = ptr_patch%edges%vertex_idx(je,jb,2)
        ibv2 = ptr_patch%edges%vertex_blk(je,jb,2)

        IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,1) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc1,ibc1,1) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,1)
          ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,1)

        ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,2)
          ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,2)

        ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,3)
          ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,3)

        ELSE
          CALL finish(method_name, "Unresolved edges%vertex(3)")
        ENDIF

        IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,1)
          ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,1)

        ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,2)
          ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,2)

        ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
          & ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
          & ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
          & (ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
          & ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
          & ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
          & ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

          ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,3)
          ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,3)

        ELSE
          CALL finish(method_name, "Unresolved edges%vertex(4)")
        ENDIF

        ilv3 = ptr_patch%edges%vertex_idx(je,jb,3)
        ibv3 = ptr_patch%edges%vertex_blk(je,jb,3)
        ilv4 = ptr_patch%edges%vertex_idx(je,jb,4)
        ibv4 = ptr_patch%edges%vertex_blk(je,jb,4)

        cc_ev3 = ptr_patch%verts%cartesian(ilv3,ibv3)
        cc_ev4 = ptr_patch%verts%cartesian(ilv4,ibv4)

        ! inverse length bewtween vertices 3 and 4
        IF (ptr_patch%geometry_info%cell_type == 3 ) THEN
          ptr_patch%edges%inv_vert_vert_length(je,jb) = 1._wp/&
            & (grid_sphere_radius*arc_length(cc_ev3,cc_ev4,ptr_patch%geometry_info))
        ENDIF

        ! next step: compute projected orientation vectors for cells and vertices
        ! bordering to each edge (incl. vertices 3 and 4 intorduced above)

        ! transform orientation vectors at local edge center to Cartesian space
        z_lon = ptr_patch%edges%center(je,jb)%lon
        z_lat = ptr_patch%edges%center(je,jb)%lat

        ! transform primal normal to cartesian vector z_nx1
        z_nx1(:)=ptr_patch%edges%primal_cart_normal(je,jb)%x(:)

        ! transform dual normal to cartesian vector z_nx2
        z_nx2(:)=ptr_patch%edges%dual_cart_normal(je,jb)%x(:)

        ! get location of cell 1

        z_lon = ptr_patch%cells%center(ilc1,ibc1)%lon
        z_lat = ptr_patch%cells%center(ilc1,ibc1)%lat

        ! compute local primal and dual normals at cell 1

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_cell(je,jb,1)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_cell(je,jb,1)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_cell(je,jb,1)%v2 = z_nv/z_norm

        ! get location of cell 2

        z_lon = ptr_patch%cells%center(ilc2,ibc2)%lon
        z_lat = ptr_patch%cells%center(ilc2,ibc2)%lat

        ! compute local primal and dual normals at cell 2

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_cell(je,jb,2)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_cell(je,jb,2)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_cell(je,jb,2)%v2 = z_nv/z_norm

        ! get location of vertex 1

        z_lon = ptr_patch%verts%vertex(ilv1,ibv1)%lon
        z_lat = ptr_patch%verts%vertex(ilv1,ibv1)%lat

        ! compute local primal and dual normals at vertex 1

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_vert(je,jb,1)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_vert(je,jb,1)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_vert(je,jb,1)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_vert(je,jb,1)%v2 = z_nv/z_norm

        ! get location of vertex 2

        z_lon = ptr_patch%verts%vertex(ilv2,ibv2)%lon
        z_lat = ptr_patch%verts%vertex(ilv2,ibv2)%lat

        ! compute local primal and dual normals at vertex 2

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_vert(je,jb,2)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_vert(je,jb,2)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_vert(je,jb,2)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_vert(je,jb,2)%v2 = z_nv/z_norm

        ! get location of vertex 3

        z_lon = ptr_patch%verts%vertex(ilv3,ibv3)%lon
        z_lat = ptr_patch%verts%vertex(ilv3,ibv3)%lat

        ! compute local primal and dual normals at vertex 3

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_vert(je,jb,3)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_vert(je,jb,3)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_vert(je,jb,3)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_vert(je,jb,3)%v2 = z_nv/z_norm

        ! get location of vertex 4

        z_lon = ptr_patch%verts%vertex(ilv4,ibv4)%lon
        z_lat = ptr_patch%verts%vertex(ilv4,ibv4)%lat

        ! compute local primal and dual normals at vertex 2

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%primal_normal_vert(je,jb,4)%v1 = z_nu/z_norm
        ptr_patch%edges%primal_normal_vert(je,jb,4)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv,ptr_patch%geometry_info)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        ptr_patch%edges%dual_normal_vert(je,jb,4)%v1 = z_nu/z_norm
        ptr_patch%edges%dual_normal_vert(je,jb,4)%v2 = z_nv/z_norm

      ENDDO

    END DO !block loop
!$OMP END DO NOWAIT

!$OMP END PARALLEL


    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%area_edge)

    ! primal_normal_cell must be sync'd before next loop,
    ! so do a sync for all above calculated quantities

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%inv_primal_edge_length)

    CALL sync_idx(sync_e,sync_v,ptr_patch,ptr_patch%edges%vertex_idx(:,:,3), &
      & ptr_patch%edges%vertex_blk(:,:,3))
    CALL sync_idx(sync_e,sync_v,ptr_patch,ptr_patch%edges%vertex_idx(:,:,4), &
      & ptr_patch%edges%vertex_blk(:,:,4))

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%inv_dual_edge_length)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%inv_vert_vert_length)

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,1)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,2)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,1)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,2)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,3)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,4)%v1)

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,1)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,2)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,1)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,2)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,3)%v1)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,4)%v1)

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,1)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_cell(:,:,2)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,1)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,2)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,3)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%primal_normal_vert(:,:,4)%v2)

    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,1)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_cell(:,:,2)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,1)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,2)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,3)%v2)
    CALL sync_patch_array(sync_e,ptr_patch,ptr_patch%edges%dual_normal_vert(:,:,4)%v2)


!$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    rl_start = 2
    rl_end = min_rlcell

    ! Final step: store primal_normal_cell also with respect to cell points
    ! in order to reduce indirect addressing during runtime
    i_startblk = ptr_patch%cells%start_block(rl_start)
    i_endblk   = ptr_patch%cells%end_block(rl_end)
    !
    ! loop through all patch cells
    !
!$OMP DO PRIVATE(jb,jc,je,i_startidx,i_endidx,ile1,ibe1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc =  i_startidx, i_endidx

        IF(.NOT.ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        DO je = 1, ptr_patch%cells%num_edges(jc,jb)

          ile1 = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe1 = ptr_patch%cells%edge_blk(jc,jb,je)


          IF ((ptr_patch%edges%cell_idx(ile1,ibe1,1) == jc) .AND. &
            & (ptr_patch%edges%cell_blk(ile1,ibe1,1) == jb)) THEN

            ptr_int%primal_normal_ec(jc,jb,je,1) = &
              & ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v1
            ptr_int%primal_normal_ec(jc,jb,je,2) = &
              & ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v2
            ptr_int%edge_cell_length(jc,jb,je) = &
              & ptr_patch%edges%edge_cell_length(ile1,ibe1,1)

          ELSE IF ((ptr_patch%edges%cell_idx(ile1,ibe1,2) == jc) .AND. &
            & (ptr_patch%edges%cell_blk(ile1,ibe1,2) == jb)) THEN

            ptr_int%primal_normal_ec(jc,jb,je,1) = &
              & ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v1
            ptr_int%primal_normal_ec(jc,jb,je,2) = &
              & ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v2
            ptr_int%edge_cell_length(jc,jb,je) = &
              & ptr_patch%edges%edge_cell_length(ile1,ibe1,2)

          ENDIF

        ENDDO

      ENDDO

    END DO !block loop
!$OMP END DO

    rl_start = 1
    rl_end = min_rlcell

    ! Compute cell-vertex distances for gradient limiter
    i_startblk = ptr_patch%cells%start_block(rl_start)
    i_endblk   = ptr_patch%cells%end_block(rl_end)
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,cc_cell,cc_v1,cc_v2,cc_v3,cc_dis1,cc_dis2,cc_dis3, &
!$OMP            z_lon,z_lat,z_nx1,z_nx2,z_norm,ilv1,ibv1,ilv2,ibv2,ilv3,ibv3) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc =  i_startidx, i_endidx

        cc_cell = ptr_patch%cells%cartesian_center(jc,jb)
        z_lon   = ptr_patch%cells%center(jc,jb)%lon
        z_lat   = ptr_patch%cells%center(jc,jb)%lat

        CALL gvec2cvec(1._wp,0._wp,z_lon,z_lat,z_nx1(1),z_nx1(2),z_nx1(3),ptr_patch%geometry_info)
        z_norm = SQRT( DOT_PRODUCT(z_nx1(1:3),z_nx1(1:3)) )
        z_nx1(1:3)  = 1._wp/z_norm * z_nx1(1:3)

        CALL gvec2cvec(0._wp,1._wp,z_lon,z_lat,z_nx2(1),z_nx2(2),z_nx2(3),ptr_patch%geometry_info)
        z_norm = SQRT( DOT_PRODUCT(z_nx2(1:3),z_nx2(1:3)) )
        z_nx2(1:3)  = 1._wp/z_norm * z_nx2(1:3)

        ilv1 = ptr_patch%cells%vertex_idx(jc,jb,1)
        ibv1 = ptr_patch%cells%vertex_blk(jc,jb,1)
        ilv2 = ptr_patch%cells%vertex_idx(jc,jb,2)
        ibv2 = ptr_patch%cells%vertex_blk(jc,jb,2)
        ilv3 = ptr_patch%cells%vertex_idx(jc,jb,3)
        ibv3 = ptr_patch%cells%vertex_blk(jc,jb,3)

        cc_v1  =  ptr_patch%verts%cartesian(ilv1,ibv1)
        cc_v2  =  ptr_patch%verts%cartesian(ilv2,ibv2)
        cc_v3  =  ptr_patch%verts%cartesian(ilv3,ibv3)

        !Find closest coordinate for flat torus case: no need as
        !these points belong to 1 triangle so they can not be cyclic
        !IF(ptr_patch%geometry_info%geometry_type==planar_torus_geometry)THEN
        !  cc_v1  = plane_torus_closest_coordinates(cc_cell%x(1:3),cc_v1%x(1:3),ptr_patch%geometry_info)
        !  cc_v2  = plane_torus_closest_coordinates(cc_cell%x(1:3),cc_v2%x(1:3),ptr_patch%geometry_info)
        !  cc_v3  = plane_torus_closest_coordinates(cc_cell%x(1:3),cc_v3%x(1:3),ptr_patch%geometry_info)
        !END IF

        cc_dis1%x(1:3) = cc_v1%x(1:3) - cc_cell%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis1%x(1:3),cc_dis1%x(1:3)) )
        cc_dis1%x(1:3) = cc_dis1%x(1:3)/z_norm

        cc_dis2%x(1:3) = cc_v2%x(1:3) - cc_cell%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis2%x(1:3),cc_dis2%x(1:3)) )
        cc_dis2%x(1:3) = cc_dis2%x(1:3)/z_norm

        cc_dis3%x(1:3) = cc_v3%x(1:3) - cc_cell%x(1:3)
        z_norm = SQRT( DOT_PRODUCT(cc_dis3%x(1:3),cc_dis3%x(1:3)) )
        cc_dis3%x(1:3) = cc_dis3%x(1:3)/z_norm

        ptr_int%cell_vert_dist(jc,1,1,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v1,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx1(1:3),cc_dis1%x(1:3))
        ptr_int%cell_vert_dist(jc,2,1,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v2,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx1(1:3),cc_dis2%x(1:3))
        ptr_int%cell_vert_dist(jc,3,1,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v3,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx1(1:3),cc_dis3%x(1:3))

        ptr_int%cell_vert_dist(jc,1,2,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v1,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx2(1:3),cc_dis1%x(1:3))
        ptr_int%cell_vert_dist(jc,2,2,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v2,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx2(1:3),cc_dis2%x(1:3))
        ptr_int%cell_vert_dist(jc,3,2,jb) = grid_sphere_radius* &
          & arc_length(cc_cell,cc_v3,ptr_patch%geometry_info)*DOT_PRODUCT(z_nx2(1:3),cc_dis3%x(1:3))

      ENDDO

    END DO !block loop
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    DO je = 1, ptr_patch%geometry_info%cell_type
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%primal_normal_ec(:,:,je,1))
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%primal_normal_ec(:,:,je,2))
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%edge_cell_length(:,:,je))
    ENDDO

  END SUBROUTINE complete_patchinfo
  !----------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Calls routines to calculate coefficients "p_int%pos_on_tplane_e" for backward tracjectory
  !! calculations depending on grid geometry.
  !!
  SUBROUTINE init_tplane_e( p_patch, p_int )
    TYPE(t_patch),     INTENT(in)    :: p_patch
    TYPE(t_int_state), INTENT(inout) :: p_int

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:init_tplane_e'

    !
    SELECT CASE(p_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)
      CALL calculate_planar_distance_at_edge( p_patch, p_int )

    CASE (sphere_geometry)
      CALL calculate_tangent_plane_at_edge( p_patch, p_int )

    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")

    END SELECT

  END SUBROUTINE init_tplane_e
  !-------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !! Initializes a tangential plane at each edge midpoint. Necessary for efficient
  !! calculation of backward trajectories and the corresponding 'flux area'.
  !!
  !! For FFSL-schemes like Miura it is necessary to calculate backward trajectories
  !! in order to determine an approximation to the flux area which is advected
  !! across each cell edge during the time step $\delta t$. In our case, this
  !! calculation is perfomed on a plane which is tangent to the edge midpoint.
  !! The coordinate axes point to the local normal and tangential direction.
  !!
  !! The position of points on this tangential plane (like the circumcenters
  !! of the neighbor cells and edge vertices) is precomputed using the
  !! gnomonic projection.
  !!
  !!
  !! Order of storage for p_int%pos_on_tplane_e(nproma,4,2,nblks_e)
  !! pos_on_tplane_e(:,1:2,:,:) :: neighboring cell centers
  !! pos_on_tplane_e(:,3:4,:,:) :: edge vertices
  !!
  SUBROUTINE calculate_tangent_plane_at_edge (p_patch, p_int)

    TYPE(t_patch),     INTENT(in)    :: p_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: p_int    !< interpolation state

    REAL(wp) ::                  &    !< origin of tangent plane in geographical coordinates
      &  origin_lon, origin_lat

    REAL(wp) ::                  &    !< point to be projected given in geographical coordinates
      &  src_lon, src_lat

    REAL(wp) ::                  &    !< cartesian coordinates of cell centers after projection
      &  xyloc_plane_c(2,2)

    REAL(wp) ::                  &    !< cartesian coordinates of edge vertices after projection
      &  xyloc_plane_ve(2,2)

    INTEGER :: ilc, ibc               !< line and block indices of neighbor cell centers
    INTEGER :: ilv, ibv               !< line and block indices of edge vertices
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: jb, je                 !< loop indices for block and edges
    INTEGER :: nc, nv                 !< loop indices for cells and and edge vertices
    INTEGER :: ne
    !-------------------------------------------------------------------------

    CALL message('mo_intp_coeffs:calculate_tangent_plane_at_edge', '')

    i_rcstartlev = 2

    ! start and end block
    i_startblk = p_patch%edges%start_block(i_rcstartlev)
    i_endblk   = p_patch%nblks_e


    !
    ! compute position of adjacent cell circumcenters and edge vertices
    ! on the tangent plane. The gnomonic projection is used. Note that we
    ! first project the points onto a local geographical (\lambda-\Phi)
    ! system. Then we rotate the coordinates into a new system with
    ! unit vectors normal and tangential to the edge.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,nc,nv,ilc,ibc,ilv,ibv,i_startidx,i_endidx,  &
!$OMP            origin_lon,origin_lat,src_lon,src_lat,            &
!$OMP            xyloc_plane_c,xyloc_plane_ve) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.p_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        ! get geographical coordinates of the edge midpoint which will be
        ! the origin of our cartesian coordinate system.
        origin_lon = p_patch%edges%center(je,jb)%lon
        origin_lat = p_patch%edges%center(je,jb)%lat

        !
        ! 1. project cell centers sharing edge je
        !
        DO nc = 1,2
          ! get line and block indices of neighbor cell
          ilc = p_patch%edges%cell_idx(je,jb,nc)
          ibc = p_patch%edges%cell_blk(je,jb,nc)

          ! get geographical coordinates of cell center
          src_lon = p_patch%cells%center(ilc,ibc)%lon
          src_lat = p_patch%cells%center(ilc,ibc)%lat

          ! project cell center into local \lambda-\Phi-system
          CALL gnomonic_proj( origin_lon, origin_lat, src_lon, src_lat, &! in
            &  xyloc_plane_c(nc,1), xyloc_plane_c(nc,2) )                ! out
        ENDDO

        !
        ! 2. project edge vertices
        !
        DO nv = 1,2
          ! get line and block indices of edge vertex
          ilv = p_patch%edges%vertex_idx(je,jb,nv)
          ibv = p_patch%edges%vertex_blk(je,jb,nv)

          ! get geographical coordinates of edge vertex
          src_lon = p_patch%verts%vertex(ilv,ibv)%lon
          src_lat = p_patch%verts%vertex(ilv,ibv)%lat

          ! projection edge vertex into local \lambda-\Phi-system
          CALL gnomonic_proj( origin_lon, origin_lat, src_lon, src_lat, &! in
            & xyloc_plane_ve(nv,1), xyloc_plane_ve(nv,2)   ) ! out
        END DO


        !
        ! 3. rotate vectors into a cartesian system where the coordinate axes
        !    point to the local normal and tangential direction.
        !

        ! centers
        !
        DO nc = 1,2
          p_int%pos_on_tplane_e(je,nc,1,jb) = grid_sphere_radius * (      &
            & xyloc_plane_c(nc,1)  * p_patch%edges%primal_normal(je,jb)%v1  &
            & + xyloc_plane_c(nc,2)  * p_patch%edges%primal_normal(je,jb)%v2 )

          p_int%pos_on_tplane_e(je,nc,2,jb) = grid_sphere_radius * (      &
            & xyloc_plane_c(nc,1)  * p_patch%edges%dual_normal(je,jb)%v1    &
            & + xyloc_plane_c(nc,2)  * p_patch%edges%dual_normal(je,jb)%v2 )
        ENDDO

        ! vertices
        !
        DO nv = 1,2
          p_int%pos_on_tplane_e(je,2+nv,1,jb) = grid_sphere_radius * (     &
            & xyloc_plane_ve(nv,1)  * p_patch%edges%primal_normal(je,jb)%v1 &
            & + xyloc_plane_ve(nv,2)  * p_patch%edges%primal_normal(je,jb)%v2 )

          p_int%pos_on_tplane_e(je,2+nv,2,jb) = grid_sphere_radius * (     &
            & xyloc_plane_ve(nv,1)  * p_patch%edges%dual_normal(je,jb)%v1   &
            & + xyloc_plane_ve(nv,2)  * p_patch%edges%dual_normal(je,jb)%v2 )
        END DO

      ENDDO ! edges
    ENDDO  ! blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO ne=1,SIZE(p_int%pos_on_tplane_e,2)
      CALL sync_patch_array(SYNC_E, p_patch, p_int%pos_on_tplane_e(:,ne,:,:))
    ENDDO

  END SUBROUTINE calculate_tangent_plane_at_edge


  !----------------------------------------------------------------------------
  !!
  SUBROUTINE calculate_planar_distance_at_edge (p_patch, p_int)

    TYPE(t_patch),     INTENT(in)    :: p_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: p_int    !< interpolation state

    !CC of points on the plane torus grid
    TYPE(t_cartesian_coordinates) :: cc_c(2), cc_ve(2), cc_origin

    !relative location of those points w.r.t cc_origin
    TYPE(t_cartesian_coordinates) :: cc_plane_c(2), cc_plane_ve(2)

    INTEGER :: ilc, ibc            !< line and block indices of neighbor cell centers
    INTEGER :: ilv, ibv            !< line and block indices of edge vertices
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: jb, je              !< loop indices for block and edges
    INTEGER :: nc, nv              !< loop index for cell center and edge vertices
    INTEGER :: ne
    !-------------------------------------------------------------------------

    CALL message('mo_intp_coeffs:calculate_planar_distance_at_edge', '')

    i_rcstartlev = 2

    ! start and end block
    i_startblk = p_patch%edges%start_block(i_rcstartlev)
    i_endblk   = p_patch%nblks_e

    !<< AD
    ! Modification for is_plane_torus for HDCP2: the planar distance between any two
    ! points is same as the Gnomonic projected distance
    ! AD >>

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,nc,nv,ilc,ibc,ilv,ibv,i_startidx,i_endidx,  &
!$OMP            cc_origin,cc_c,cc_plane_c,cc_ve,cc_plane_ve) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.p_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        ! the edge midpoint will be the origin of our local 2D
        ! cartesian system.
        cc_origin = p_patch%edges%cartesian_center(je,jb)

        !
        ! 1. CC of neighboring cell centers; z index is 0
        !
        DO nc = 1,2
          ! get line and block indices of neighbor cells
          ilc = p_patch%edges%cell_idx(je,jb,nc)
          ibc = p_patch%edges%cell_blk(je,jb,nc)

          ! get CC of the cell centers
          cc_c(nc) = p_patch%cells%cartesian_center(ilc,ibc)

          ! now calculate the separation vector between the edge and cell centers
          cc_c(nc) = plane_torus_closest_coordinates(cc_origin%x,cc_c(nc)%x,p_patch%geometry_info)
          cc_plane_c(nc)%x(:) = cc_c(nc)%x(:) - cc_origin%x(:)
        ENDDO


        !
        ! 2. Edge vertices
        !
        DO nv = 1,2
          ! get line and block indices of edge vertices
          ilv = p_patch%edges%vertex_idx(je,jb,nv)
          ibv = p_patch%edges%vertex_blk(je,jb,nv)

          ! get CC coordinates of edge vertices
          cc_ve(nv) = p_patch%verts%cartesian(ilv,ibv)

          !now calculate the separation vector between the edge and cell centers
          cc_ve(nv) =  plane_torus_closest_coordinates(cc_origin%x,cc_ve(nv)%x,p_patch%geometry_info)
          cc_plane_ve(nv)%x(:) = cc_ve(nv)%x(:) - cc_origin%x(:)
        END DO

        !
        ! 3. rotate these vectors into a new local cartesian system. In this rotated
        !    system the coordinate axes point into the local normal and tangential
        !    direction at each edge. All these vectors are 2D so no point using the
        !    z coordinate
        !

        ! centers
        !
        DO nc = 1,2
          p_int%pos_on_tplane_e(je,nc,1,jb) = SUM(      &
            &  cc_plane_c(nc)%x(1:2)  * p_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

          p_int%pos_on_tplane_e(je,nc,2,jb) = SUM(      &
            &  cc_plane_c(nc)%x(1:2)  * p_patch%edges%dual_cart_normal(je,jb)%x(1:2) )
        ENDDO

        ! vertices
        !
        DO nv = 1,2
          p_int%pos_on_tplane_e(je,2+nv,1,jb) = SUM(      &
            &  cc_plane_ve(nv)%x(1:2)  * p_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

          p_int%pos_on_tplane_e(je,2+nv,2,jb) = SUM(      &
            &  cc_plane_ve(nv)%x(1:2)  * p_patch%edges%dual_cart_normal(je,jb)%x(1:2) )
        END DO

      ENDDO ! edges
    ENDDO  ! blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO ne=1,SIZE(p_int%pos_on_tplane_e,2)
      CALL sync_patch_array(SYNC_E, p_patch, p_int%pos_on_tplane_e(:,ne,:,:))
    ENDDO


  END SUBROUTINE calculate_planar_distance_at_edge



  !----------------------------------------------------------------------------
  !>
  !! Calls routines to calculate coefficients "ptr_int%pos_on_tplane_c_edge".
  !! Bifurcation for calculations depending on grid geometry.
  !!
  SUBROUTINE init_tplane_c( ptr_patch, ptr_int )
    TYPE(t_patch),     INTENT(in)    :: ptr_patch
    TYPE(t_int_state), INTENT(inout) :: ptr_int

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:init_tplane_c'

    !
    SELECT CASE(ptr_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)
      CALL init_tplane_c_torus ( ptr_patch, ptr_int%pos_on_tplane_c_edge )
    CASE (sphere_geometry)
      CALL init_tplane_c_sphere ( ptr_patch, ptr_int%pos_on_tplane_c_edge )
    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT

  END SUBROUTINE init_tplane_c


  !----------------------------------------------------------------------------
  !!
  !! Initializes a plane tangent to the edge midpoint for every edge.
  !! It is used to compute overlap regions between flux areas and
  !! the model grid.
  !!
  !! Projection method: gnomonic
  !!
  !! For every edge, the following set of points is projected onto a
  !! plane, which has its origin at the edge midpoint:
  !! - vertices of the two grid cells which share the given edge
  !! - cell circumcenters of the four cells which are adjacent to these
  !!   two grid cells (butterfly cells).
  !!
  !! In a second step, these points are expressed in two different cartesian
  !! systems, with the origin at the cell circumcenter of the cells that share
  !! the given edge. Unit vectors point into local north and local east
  !! direction. This change of coordinate basis is done for convenience.
  !!
  !! Order of storage for pos_on_tplane_c_edge:
  !! pos_on_tplane_c_edge(nproma,nblks_e,ncells=2,npts=5)%lon/lat
  !!
  !! - cell ordering according to edge%cell_idx/blk
  !! - npts 1-3: cell vertices
  !!   - vertex one and two bound the given edge.
  !!     Ordering according to edge%vertex_idx/blk
  !! - npts 4-5: coordinates of neighboring cell centers
  !!   - those 2 neighbors that do not share the given edge,
  !!     but share vertex 1 and 2, respectively
  !! - Ordering of coordinates: lon=x, lat=y
  !!
  SUBROUTINE init_tplane_c_sphere (ptr_patch, pos_on_tplane_c_edge)
#ifdef __INTEL_COMPILER
!DIR$ OPTIMIZE:2
#endif
    TYPE(t_patch),                    INTENT(in)    :: &
      &  ptr_patch

    TYPE(t_geographical_coordinates), INTENT(inout) :: &
      &  pos_on_tplane_c_edge(:,:,:,:)  !< positions of vertices and butterfly neighbors
                                        !< on local plane tangent to the edge midpoint
                                        !< (nproma,nblks_e,2,5)
    ! local
    TYPE(t_geographical_coordinates) :: &
      &  origin                    !< origin of plane in geographical coordinates

    TYPE(t_geographical_coordinates) :: &
      &  gc_v, &                   !< geographical coords of edge-vertices
      &  gc_c, &                   !< geographical coords of adjacent cells
      &  gc_bf                     !< geographical coords of butterfly cell centers

    REAL(wp) ::               &
      &  xyloc_v_plane(4,2),  &    !< local 2D cartesian coords of edge-vertices
      &  xyloc_c_plane(2,2),  &    !< local 2D cartesian coords of adjacent cells
      &  xyloc_bf1_plane(2,2),&    !< local 2D cartesian coords of butterfly cell centers
      &  xyloc_bf2_plane(2,2)      !< local north and local east directions
                                   !< (nelements,(/x,y/))

    REAL(wp) ::                  & !< same, but for rotated system
      &  xyloc_v_plane_nt(4,2),  & !< normal and tangential directions
      &  xyloc_c_plane_nt(2,2),  & !< (nelements,(/x,y/))
      &  xyloc_bf1_plane_nt(2,2),&
      &  xyloc_bf2_plane_nt(2,2)

    REAL(wp) ::                 &  !< same, but for translated rotated system
      &  xyloc_v_trans1(4,2),   &  !< origin at adjacent cell center 1
      &  xyloc_v_trans2(4,2),   &  !< origin at adjacent cell center 2
      &  xyloc_bf1_trans1(2,2), &  !< origin at adjacent cell center 1
      &  xyloc_bf2_trans2(2,2)     !< origin at adjacent cell center 2
                                   !< (nelements,(/x,y/))

    TYPE(t_tangent_vectors) ::  &  !< primal and dual normal for adjacent cells
      &  pn_cell1, pn_cell2,    &
      &  dn_cell1, dn_cell2

    INTEGER :: jb, je ,nv, nc, nb              ! loop variables
    INTEGER :: ilv, ibv, ilc, ibc, ilbf, ibbf  ! line and block indices
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end

    !-------------------------------------------------------------------------

    CALL message('mo_interpolation:init_tplane_c_sphere','')

    rl_start = 3
    rl_end   = min_rledge

    i_startblk = ptr_patch%edges%start_block(rl_start)
    i_endblk   = ptr_patch%edges%end_block(rl_end)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,nv,nc,nb,i_startidx,i_endidx,origin,                   &
!$OMP            ilv,ibv,gc_v,xyloc_v_plane,ilc,ibc,gc_c,xyloc_c_plane,       &
!$OMP            ilbf,ibbf,gc_bf,xyloc_bf1_plane,xyloc_bf2_plane,             &
!$OMP            xyloc_v_plane_nt,xyloc_c_plane_nt,xyloc_bf1_plane_nt,        &
!$OMP            xyloc_bf2_plane_nt,xyloc_v_trans1,xyloc_v_trans2,            &
!$OMP            xyloc_bf1_trans1,xyloc_bf2_trans2,pn_cell1,pn_cell2,         &
!$OMP            dn_cell1,dn_cell2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        !
        ! 1. project edge-vertices, the centers of the adjacent cells
        !    and the butterfly cell centers onto a plane tangent
        !    to the edge midpoint.
        !

        ! get geographical coordinates of the edge midpoint, which will be
        ! the origin of our 2D cartesian coordinate system.
        origin%lon = ptr_patch%edges%center(je,jb)%lon
        origin%lat = ptr_patch%edges%center(je,jb)%lat

        ! get geographical coordinates of edge vertices
        ! and project them into the edge-based local system.
        ! (including the non-edge-aligned vertices of the adjacent cells)
        DO nv=1,4
          ilv = ptr_patch%edges%vertex_idx(je,jb,nv)
          ibv = ptr_patch%edges%vertex_blk(je,jb,nv)

          ! vertex coordinates
          gc_v%lon = ptr_patch%verts%vertex(ilv,ibv)%lon
          gc_v%lat = ptr_patch%verts%vertex(ilv,ibv)%lat

          ! project vertex into edge-based local system
          CALL gnomonic_proj( origin%lon, origin%lat,                  &! in
            &                 gc_v%lon, gc_v%lat,                      &! in
            &                 xyloc_v_plane(nv,1), xyloc_v_plane(nv,2) )! out
        ENDDO


        ! get geographical coordinates of adjacent cell centers
        ! and project them into the edge-based local system.
        DO nc=1,2
          ! get line and block indices of adjacent cells
          ilc = ptr_patch%edges%cell_idx(je,jb,nc)
          ibc = ptr_patch%edges%cell_blk(je,jb,nc)

          ! get geographical coordinates of cell center
          gc_c%lon = ptr_patch%cells%center(ilc,ibc)%lon
          gc_c%lat = ptr_patch%cells%center(ilc,ibc)%lat

          ! project cell center into edge-based local system
          CALL gnomonic_proj( origin%lon, origin%lat,                  &! in
            &                 gc_c%lon, gc_c%lat,                      &! in
            &                 xyloc_c_plane(nc,1), xyloc_c_plane(nc,2) )! out
        ENDDO


        !
        ! Projection of butterfly cell centers
        !
        ! project cell centers which are adjacent to the cells sharing edge je.
        ! (from butterfly_idx)

        ! neighbors of edge neighbor 1
        !
        DO nb=1,2
          ! get line and block indices of the neighbors of edge-neighbor 1
          ilbf = ptr_patch%edges%butterfly_idx(je,jb,1,nb)
          ibbf = ptr_patch%edges%butterfly_blk(je,jb,1,nb)

          ! get geographical coordinates of cell centers
          gc_bf%lon = ptr_patch%cells%center(ilbf,ibbf)%lon
          gc_bf%lat = ptr_patch%cells%center(ilbf,ibbf)%lat

          ! project butterfly cell centers into edge-based local system
          CALL gnomonic_proj( origin%lon, origin%lat,                      &! in
            &                 gc_bf%lon, gc_bf%lat,                        &! in
            &                 xyloc_bf1_plane(nb,1), xyloc_bf1_plane(nb,2) )! out
        ENDDO

        ! neighbors of edge neighbor 2
        !
        DO nb=1,2
          ! get line and block indices of the neighbors of edge-neighbor 2
          ilbf = ptr_patch%edges%butterfly_idx(je,jb,2,nb)
          ibbf = ptr_patch%edges%butterfly_blk(je,jb,2,nb)

          ! get geographical coordinates of cell centers
          gc_bf%lon = ptr_patch%cells%center(ilbf,ibbf)%lon
          gc_bf%lat = ptr_patch%cells%center(ilbf,ibbf)%lat

          ! project butterfly cell centers into edge-based local system
          CALL gnomonic_proj( origin%lon, origin%lat,                      &! in
            &                 gc_bf%lon, gc_bf%lat,                        &! in
            &                 xyloc_bf2_plane(nb,1), xyloc_bf2_plane(nb,2) )! out
        ENDDO



        !
        ! 2. rotate these vectors into a new local cartesian system, with
        !    the coordinate axes pointing into the direction normal and
        !    tangential to the edge.
        !

        ! vertices
        !
        DO nv=1,4
          xyloc_v_plane_nt(nv,1) =                                            &
            &   xyloc_v_plane(nv,1) * ptr_patch%edges%primal_normal(je,jb)%v1 &
            & + xyloc_v_plane(nv,2) * ptr_patch%edges%primal_normal(je,jb)%v2

          xyloc_v_plane_nt(nv,2) =                                            &
            &   xyloc_v_plane(nv,1) * ptr_patch%edges%dual_normal(je,jb)%v1   &
            & + xyloc_v_plane(nv,2) * ptr_patch%edges%dual_normal(je,jb)%v2
        END DO

        ! centers
        !
        DO nc=1,2
          xyloc_c_plane_nt(nc,1) =                                              &
            &   xyloc_c_plane(nc,1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
            & + xyloc_c_plane(nc,2)  * ptr_patch%edges%primal_normal(je,jb)%v2

          xyloc_c_plane_nt(nc,2) =                                              &
            &   xyloc_c_plane(nc,1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
            & + xyloc_c_plane(nc,2)  * ptr_patch%edges%dual_normal(je,jb)%v2
        ENDDO



        ! butterfly cells
        !
        DO nb=1,2
          xyloc_bf1_plane_nt(nb,1) = &
            &   xyloc_bf1_plane(nb,1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
            & + xyloc_bf1_plane(nb,2)  * ptr_patch%edges%primal_normal(je,jb)%v2

          xyloc_bf1_plane_nt(nb,2) = &
            &   xyloc_bf1_plane(nb,1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
            & + xyloc_bf1_plane(nb,2)  * ptr_patch%edges%dual_normal(je,jb)%v2

          xyloc_bf2_plane_nt(nb,1) = &
            &   xyloc_bf2_plane(nb,1)  * ptr_patch%edges%primal_normal(je,jb)%v1  &
            & + xyloc_bf2_plane(nb,2)  * ptr_patch%edges%primal_normal(je,jb)%v2

          xyloc_bf2_plane_nt(nb,2) = &
            &   xyloc_bf2_plane(nb,1)  * ptr_patch%edges%dual_normal(je,jb)%v1    &
            & + xyloc_bf2_plane(nb,2)  * ptr_patch%edges%dual_normal(je,jb)%v2
        ENDDO



        ! 3. Calculate distance vectors in a translated coordinate system.
        !    This is done twice. The origin is located once at the circumcenter
        !    of the adjacent cell 1 and once cell 2. The distance vectors
        !    point from the cell center to the vertices and butterfly cell centers.
        !
        ! vertices
        DO nv=1,4
          xyloc_v_trans1(nv,1) = xyloc_v_plane_nt(nv,1) - xyloc_c_plane_nt(1,1)
          xyloc_v_trans1(nv,2) = xyloc_v_plane_nt(nv,2) - xyloc_c_plane_nt(1,2)
          xyloc_v_trans2(nv,1) = xyloc_v_plane_nt(nv,1) - xyloc_c_plane_nt(2,1)
          xyloc_v_trans2(nv,2) = xyloc_v_plane_nt(nv,2) - xyloc_c_plane_nt(2,2)
        ENDDO

        ! butterfly cells
        DO nb=1,2
          xyloc_bf1_trans1(nb,1) = xyloc_bf1_plane_nt(nb,1) - xyloc_c_plane_nt(1,1)
          xyloc_bf1_trans1(nb,2) = xyloc_bf1_plane_nt(nb,2) - xyloc_c_plane_nt(1,2)

          xyloc_bf2_trans2(nb,1) = xyloc_bf2_plane_nt(nb,1) - xyloc_c_plane_nt(2,1)
          xyloc_bf2_trans2(nb,2) = xyloc_bf2_plane_nt(nb,2) - xyloc_c_plane_nt(2,2)
        ENDDO


        ! 4. Rotate points into coordinate system pointing into local north
        !    and local east direction. Store in edge-based data structure. This
        !    is done twice (for both adjacent cell centers).
        ! e_n= pn_cell1%v1 * e_\lambda + pn_cell1%v2 * e_\phi
        ! e_t= dn_cell1%v1 * e_\lambda + dn_cell1%v2 * e_\phi
        pn_cell1%v1 = ptr_patch%edges%primal_normal_cell(je,jb,1)%v1
        pn_cell1%v2 = ptr_patch%edges%primal_normal_cell(je,jb,1)%v2
        dn_cell1%v1 = ptr_patch%edges%dual_normal_cell(je,jb,1)%v1
        dn_cell1%v2 = ptr_patch%edges%dual_normal_cell(je,jb,1)%v2

        pn_cell2%v1 = ptr_patch%edges%primal_normal_cell(je,jb,2)%v1
        pn_cell2%v2 = ptr_patch%edges%primal_normal_cell(je,jb,2)%v2
        dn_cell2%v1 = ptr_patch%edges%dual_normal_cell(je,jb,2)%v1
        dn_cell2%v2 = ptr_patch%edges%dual_normal_cell(je,jb,2)%v2

        ! vertices
        !
        ! components in longitudinal direction (cell 1)
        pos_on_tplane_c_edge(je,jb,1,1:3)%lon =  grid_sphere_radius &
          & *( xyloc_v_trans1(1:3,1) * pn_cell1%v1                  &
          & +  xyloc_v_trans1(1:3,2) * dn_cell1%v1 )

        ! components in latitudinal direction (cell 1)
        pos_on_tplane_c_edge(je,jb,1,1:3)%lat =  grid_sphere_radius &
          & *( xyloc_v_trans1(1:3,1) * pn_cell1%v2                  &
          & +  xyloc_v_trans1(1:3,2) * dn_cell1%v2 )

        ! components in longitudinal direction (cell 2)
        pos_on_tplane_c_edge(je,jb,2,1:3)%lon =  grid_sphere_radius &
          & *( xyloc_v_trans2((/1,2,4/),1) * pn_cell2%v1            &
          & +  xyloc_v_trans2((/1,2,4/),2) * dn_cell2%v1 )

        ! components in latitudinal direction (cell 2)
        pos_on_tplane_c_edge(je,jb,2,1:3)%lat =  grid_sphere_radius &
          & *( xyloc_v_trans2((/1,2,4/),1) * pn_cell2%v2            &
          & +  xyloc_v_trans2((/1,2,4/),2) * dn_cell2%v2 )


        ! butterfly cells
        !
        ! components in longitudinal direction (cell 1)
        pos_on_tplane_c_edge(je,jb,1,4:5)%lon =  grid_sphere_radius * &
          & ( xyloc_bf1_trans1(1:2,1) * pn_cell1%v1                   &
          & + xyloc_bf1_trans1(1:2,2) * dn_cell1%v1 )

        ! components in latitudinal direction (cell 1)
        pos_on_tplane_c_edge(je,jb,1,4:5)%lat =  grid_sphere_radius * &
          & ( xyloc_bf1_trans1(1:2,1) * pn_cell1%v2                   &
          & + xyloc_bf1_trans1(1:2,2) * dn_cell1%v2 )

        ! components in longitudinal direction (cell 2)
        pos_on_tplane_c_edge(je,jb,2,4:5)%lon =  grid_sphere_radius * &
          & ( xyloc_bf2_trans2(1:2,1) * pn_cell2%v1                   &
          & + xyloc_bf2_trans2(1:2,2) * dn_cell2%v1 )

        ! components in latitudinal direction (cell 2)
        pos_on_tplane_c_edge(je,jb,2,4:5)%lat =  grid_sphere_radius * &
          & ( xyloc_bf2_trans2(1:2,1) * pn_cell2%v2                   &
          & + xyloc_bf2_trans2(1:2,2) * dn_cell2%v2 )
      ENDDO  ! je
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


    DO nv=1,5
      DO nc=1,2
        CALL sync_patch_array(sync_e, ptr_patch, pos_on_tplane_c_edge(:,:,nc,nv)%lon)
        CALL sync_patch_array(sync_e, ptr_patch, pos_on_tplane_c_edge(:,:,nc,nv)%lat)
      ENDDO
    ENDDO

  END SUBROUTINE init_tplane_c_sphere


  !----------------------------------------------------------------------------
  !!AD> This is torus version of "init_tplane_c_sphere"- everything remains same
  !!    just that for flat torus geometry there is no need of gnomonic projection.
  !!    Radial vector between two point is just coordinate difference.
  !!
  !!    Please note that the latest changes to init_tplane_c_sphere have not been
  !!    taken over to init_tplane_c_torus. The changes are related to the projection
  !!    of the butterfly cell centers, which are now projected in the same way
  !!    as all other points. Even though this modification has not been included
  !!    the torus version is still fully functionable.
  !!
  SUBROUTINE init_tplane_c_torus (ptr_patch, pos_on_tplane_c_edge)

    TYPE(t_patch),                    INTENT(in)    :: ptr_patch  !< patch

    TYPE(t_geographical_coordinates), INTENT(inout) :: &
      &  pos_on_tplane_c_edge(:,:,:,:)  !< positions of vertices and butterfly neighbors
                                        !< on local plane tangent to the edge midpoint
                                        !< (nproma,nblks_e,2,5)

    !CC of points in the stencil
    TYPE(t_cartesian_coordinates) :: cc_edge, cc_vert(4), cc_cell(2), &
                                     cc_bf1(2), cc_bf2(2)

    !Distance vectors
    TYPE(t_cartesian_coordinates) :: dist_edge_cell(2), dist_edge_vert(4), &
                                     dist_bf1_cell1(2), dist_bf2_cell2(2)

    REAL(wp) ::   &              !< same, but for rotated system (normal-tangential)
      & xyloc_plane_nt_v(4,2)

    REAL(wp) ::   &              !< same but for rotated system (normal-tangential)
      & xyloc_plane_nt_n1(2), xyloc_plane_nt_n2(2)

    REAL(wp) ::   &              !< coords of edge-vertices in translated system
      & xyloc_trans1_v(4,2), xyloc_trans2_v(4,2)

    REAL(wp) ::   &              !< primal and dual normals for neighboring cells
      & pn_cell1(2), pn_cell2(2), dn_cell1(2), dn_cell2(2)

    INTEGER :: ilv(4), ibv(4)
    INTEGER :: ilc1, ilc2, ibc1, ibc2
    INTEGER :: ilc_bf1(2), ilc_bf2(2), ibc_bf1(2), ibc_bf2(2)
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: nv, nc            !< loop indices for vertices and cells
    INTEGER :: jb, je            !< loop indices for block and edges

    !-------------------------------------------------------------------------

    CALL message('mo_interpolation:', 'init_tplane_c_torus')

    i_rcstartlev = 2

    ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_e

    !<< AD
    ! Modification for is_plane_torus for HDCP2: the planar distance between any two
    ! points is same as the Gnomonic projected distance
    ! AD >>

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,nv,i_startidx,i_endidx,cc_edge,ilv,ibv,          &
!$OMP            cc_vert,dist_edge_vert,ilc1,ilc2,ibc1,ibc2,cc_cell,    &
!$OMP            dist_edge_cell,xyloc_plane_nt_n1,xyloc_plane_nt_n2,    &
!$OMP            xyloc_plane_nt_v,xyloc_trans1_v,        &
!$OMP            xyloc_trans2_v,pn_cell1,pn_cell2,dn_cell1,dn_cell2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        !
        ! 1. Get distance vector of edge-vertices and the centers of the neighboring cells
        !    w.r.t to the edge midpoint.
        !

        ! get CC of edge midpoint
        cc_edge = ptr_patch%edges%cartesian_center(je,jb)

        ! get line and block indices of edge vertices (including the
        ! non-edge-aligned vertices of the neighboring cells
        ilv(1:4)=ptr_patch%edges%vertex_idx(je,jb,1:4)
        ibv(1:4)=ptr_patch%edges%vertex_blk(je,jb,1:4)

        ! get CC of edge vertices
        DO nv=1,4
          cc_vert(nv) = ptr_patch%verts%cartesian(ilv(nv),ibv(nv))

          !Get the actual location of the vertex w.r.t the edge cc
          cc_vert(nv) = plane_torus_closest_coordinates(cc_edge%x,cc_vert(nv)%x, &
                                                         ptr_patch%geometry_info)
          !Get the distance vector between cc_edge and cc_vets
          dist_edge_vert(nv)%x(:) = cc_vert(nv)%x(:) - cc_edge%x(:)
        ENDDO

        ! get line and block indices of neighbour cells
        ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
        ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
        ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
        ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

        ! get CC of the two cell centers
        cc_cell(1)   = ptr_patch%cells%cartesian_center(ilc1,ibc1)
        cc_cell(2)   = ptr_patch%cells%cartesian_center(ilc2,ibc2)

        !Get the actual location of the cell w.r.t the edge cc
        cc_cell(1) = plane_torus_closest_coordinates(cc_edge%x,cc_cell(1)%x, &
                                                     ptr_patch%geometry_info)
        cc_cell(2) = plane_torus_closest_coordinates(cc_edge%x,cc_cell(2)%x, &
                                                     ptr_patch%geometry_info)

        !Get the distance vector between cc_edge and cc_cell
        dist_edge_cell(1)%x(:) = cc_cell(1)%x(:) - cc_edge%x(:)
        dist_edge_cell(2)%x(:) = cc_cell(2)%x(:) - cc_edge%x(:)


        !
        ! 2. rotate these vectors into a new local cartesian system. In this rotated
        !    system the coordinate axes point into the local normal and tangential
        !    direction at each edge. All these vectors are 2D so no point using the
        !    z coordinate
        !

        ! centers
        !
        xyloc_plane_nt_n1(1) =  SUM( dist_edge_cell(1)%x(1:2) *  &
             ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

        xyloc_plane_nt_n1(2) =  SUM( dist_edge_cell(1)%x(1:2) *  &
             ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )

        xyloc_plane_nt_n2(1) =  SUM( dist_edge_cell(2)%x(1:2) *  &
             ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

        xyloc_plane_nt_n2(2) =  SUM( dist_edge_cell(2)%x(1:2) *  &
             ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )

        ! vertices
        !
        DO nv = 1,4
           xyloc_plane_nt_v(nv,1) =  SUM( dist_edge_vert(nv)%x(1:2) *  &
                ptr_patch%edges%primal_cart_normal(je,jb)%x(1:2) )

           xyloc_plane_nt_v(nv,2) =  SUM( dist_edge_vert(nv)%x(1:2) *  &
                ptr_patch%edges%dual_cart_normal(je,jb)%x(1:2) )
        END DO

        !
        !AD: From this point on the sphere and the torus version are same
        !    except for the multiplication with grid_sphere_radius

        ! 3. Calculate position of vertices in a translated coordinate system.
        !    This is done twice. The origin is located once at the circumcenter
        !    of the neighboring cell 1 and once cell 2. The distance vectors point
        !    from the cell center to the vertices.
        xyloc_trans1_v(1:4,1) = xyloc_plane_nt_v(1:4,1) - xyloc_plane_nt_n1(1)
        xyloc_trans1_v(1:4,2) = xyloc_plane_nt_v(1:4,2) - xyloc_plane_nt_n1(2)

        xyloc_trans2_v(1:4,1) = xyloc_plane_nt_v(1:4,1) - xyloc_plane_nt_n2(1)
        xyloc_trans2_v(1:4,2) = xyloc_plane_nt_v(1:4,2) - xyloc_plane_nt_n2(2)



        ! 4. Rotate points into coordinate system pointing into local north
        !    and local east direction. Store in edge-based data structure. This
        !    is done twice (for both neighboring cells).
        ! e_n= pn_cell1(1)*e_\lambda + pn_cell1(2)*e_\phi
        ! e_t= dn_cell1(1)*e_\lambda + dn_cell1(2)*e_\phi
        pn_cell1(1) = ptr_patch%edges%primal_normal_cell(je,jb,1)%v1
        pn_cell1(2) = ptr_patch%edges%primal_normal_cell(je,jb,1)%v2
        dn_cell1(1) = ptr_patch%edges%dual_normal_cell(je,jb,1)%v1
        dn_cell1(2) = ptr_patch%edges%dual_normal_cell(je,jb,1)%v2

        pn_cell2(1) = ptr_patch%edges%primal_normal_cell(je,jb,2)%v1
        pn_cell2(2) = ptr_patch%edges%primal_normal_cell(je,jb,2)%v2
        dn_cell2(1) = ptr_patch%edges%dual_normal_cell(je,jb,2)%v1
        dn_cell2(2) = ptr_patch%edges%dual_normal_cell(je,jb,2)%v2

        ! components in longitudinal direction (cell 1)
        pos_on_tplane_c_edge(je,jb,1,1:3)%lon =          &
          &  ( xyloc_trans1_v(1:3,1) * pn_cell1(1)       &
          & +  xyloc_trans1_v(1:3,2) * dn_cell1(1) )

        ! components in latitudinal direction (cell 1)
        pos_on_tplane_c_edge(je,jb,1,1:3)%lat =          &
          &  ( xyloc_trans1_v(1:3,1) * pn_cell1(2)       &
          & +  xyloc_trans1_v(1:3,2) * dn_cell1(2) )

        ! components in longitudinal direction (cell 2)
        pos_on_tplane_c_edge(je,jb,2,1:3)%lon =          &
          &  ( xyloc_trans2_v((/1,2,4/),1) * pn_cell2(1) &
          & +  xyloc_trans2_v((/1,2,4/),2) * dn_cell2(1) )

        ! components in latitudinal direction (cell 2)
        pos_on_tplane_c_edge(je,jb,2,1:3)%lat =          &
          &  ( xyloc_trans2_v((/1,2,4/),1) * pn_cell2(2) &
          & +  xyloc_trans2_v((/1,2,4/),2) * dn_cell2(2) )

      ENDDO  ! je
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL



    !
    ! Projection of butterfly cell centers
    !

    i_rcstartlev = 3

    ! start and end block
    i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_e

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ilc2,ibc1,ibc2,cc_cell,&
!$OMP            ilc_bf1,ibc_bf1,ilc_bf2,ibc_bf2,cc_bf1,    &
!$OMP            cc_bf2,dist_bf1_cell1,dist_bf2_cell2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO je = i_startidx, i_endidx

        IF(.NOT.ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE


        ! get line and block indices of edge-neighbours (cells)
        ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
        ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
        ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
        ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

        ! get CC of edge-neighbor 1
        cc_cell(1) = ptr_patch%cells%cartesian_center(ilc1,ibc1)

        ! get CC of edge-neighbor 2
        cc_cell(2) = ptr_patch%cells%cartesian_center(ilc2,ibc2)


        ! 5a. project cell centers adjacent to the cells sharing edge je.
        !    (from butterfly_idx)


        ! get line and block indices of the neighbors of edge-neighbor 1
        ilc_bf1(1:2) = ptr_patch%edges%butterfly_idx(je,jb,1,1:2)
        ibc_bf1(1:2) = ptr_patch%edges%butterfly_blk(je,jb,1,1:2)

        ! get line and block indices of the neighbors of edge-neighbor 2
        ilc_bf2(1:2) = ptr_patch%edges%butterfly_idx(je,jb,2,1:2)
        ibc_bf2(1:2) = ptr_patch%edges%butterfly_blk(je,jb,2,1:2)


        ! get CC cell centers (neighbor 1)
        cc_bf1(1) = ptr_patch%cells%cartesian_center(ilc_bf1(1),ibc_bf1(1))
        cc_bf1(2) = ptr_patch%cells%cartesian_center(ilc_bf1(2),ibc_bf1(2))

        ! get CC cell centers (neighbor 1)
        cc_bf2(1) = ptr_patch%cells%cartesian_center(ilc_bf2(1),ibc_bf2(1))
        cc_bf2(2) = ptr_patch%cells%cartesian_center(ilc_bf2(2),ibc_bf2(2))


        ! project cell centers into cell-based local system (edge-neighbor 1)
        cc_bf1(1) = plane_torus_closest_coordinates(cc_cell(1)%x, cc_bf1(1)%x, &
                                                   ptr_patch%geometry_info)
        dist_bf1_cell1(1)%x(:) = cc_bf1(1)%x(:) - cc_cell(1)%x(:)

        cc_bf1(2) = plane_torus_closest_coordinates(cc_cell(1)%x, cc_bf1(2)%x, &
                                                   ptr_patch%geometry_info)
        dist_bf1_cell1(2)%x(:) = cc_bf1(2)%x(:) - cc_cell(1)%x(:)


        ! project cell centers into cell-based local system (edge-neighbor 2)

        cc_bf2(1) = plane_torus_closest_coordinates(cc_cell(2)%x, cc_bf2(1)%x, &
                                                   ptr_patch%geometry_info)
        dist_bf2_cell2(1)%x(:) = cc_bf2(1)%x(:) - cc_cell(2)%x(:)

        cc_bf2(2) = plane_torus_closest_coordinates(cc_cell(2)%x, cc_bf2(2)%x, &
                                                   ptr_patch%geometry_info)
        dist_bf2_cell2(2)%x(:) = cc_bf2(2)%x(:) - cc_cell(2)%x(:)


        ! components in longitudinal direction (or X for torus) (cell 1)
        pos_on_tplane_c_edge(je,jb,1,4:5)%lon = dist_bf1_cell1(1:2)%x(1)

        ! components in latitudinal direction (or Y for torus) (cell 1)
        pos_on_tplane_c_edge(je,jb,1,4:5)%lat = dist_bf1_cell1(1:2)%x(2)


        ! components in longitudinal (or X for torus) direction (cell 2)
        pos_on_tplane_c_edge(je,jb,2,4:5)%lon =  dist_bf2_cell2(1:2)%x(1)

        ! components in latitudinal (or Y for torus) direction (cell 2)
        pos_on_tplane_c_edge(je,jb,2,4:5)%lat =  dist_bf2_cell2(1:2)%x(2)

      ENDDO  ! je
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    DO nv=1,5
      DO nc=1,2
        CALL sync_patch_array(sync_e, ptr_patch, pos_on_tplane_c_edge(:,:,nc,nv)%lon)
        CALL sync_patch_array(sync_e, ptr_patch, pos_on_tplane_c_edge(:,:,nc,nv)%lat)
      ENDDO
    ENDDO


  END SUBROUTINE init_tplane_c_torus

  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  !! Primal cell quadrature points and weights
  !!
  !! Computes quadrature points and weights for triangular grid cells.
  !! Quadrature points and weights are provided for accurately integrating
  !! linear, quadratic and cubic functions. This is necessary for initializing
  !! idealized testcases.
  !!
  !! @par Literature
  !! Numerical Methods in Engineering with Python, Jaan Kiusalaas (2005),
  !! 233-247
  SUBROUTINE tri_quadrature_pts (ptr_patch, ptr_int)

    TYPE(t_patch), INTENT(inout) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(inout) :: ptr_int  !< interpolation state

    REAL(wp) ::  alpha_l(3),        & !< area coordinates for quadrature up to
      & alpha_q(3,3), alpha_c(3,4)   !< fourth order
    !< (n_area_coords,n_pts))

    TYPE(t_cartesian_coordinates)    :: z_vert_cc(3) ! cell vertices in cartesian
    ! coordinates
    TYPE(t_cartesian_coordinates)    :: z_quad_cc  ! triangle quadrature point in cartesian
    ! coordinates
    TYPE(t_geographical_coordinates) :: z_quad_gg  ! triangle quadrature point in geographical
    ! coordinates

    INTEGER :: ilv, ibv               !< line and block indices of cell vertices
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: nq, nv                 !< loop index for quadrature points and
    !< cell vertices
    INTEGER :: jc, jb                 !< loop index for cells

    !-------------------------------------------------------------------------

    CALL message('mo_interpolation:tri_quadrature_pts', '')

    ! set area coordinates
    !
    ! linear
    alpha_l(1) = 1._wp/3._wp
    alpha_l(2) = 1._wp/3._wp
    alpha_l(3) = 1._wp/3._wp

    !
    ! quadratic
    !
    alpha_q(1,1) = 0.5_wp
    alpha_q(2,1) = 0._wp
    alpha_q(3,1) = 0.5_wp

    alpha_q(1,2) = 0.5_wp
    alpha_q(2,2) = 0.5_wp
    alpha_q(3,2) = 0._wp

    alpha_q(1,3) = 0._wp
    alpha_q(2,3) = 0.5_wp
    alpha_q(3,3) = 0.5_wp

    !
    ! cubic
    !
    alpha_c(1,1) = 1._wp/3._wp
    alpha_c(2,1) = 1._wp/3._wp
    alpha_c(3,1) = 1._wp/3._wp

    alpha_c(1,2) = 1._wp/5._wp
    alpha_c(2,2) = 1._wp/5._wp
    alpha_c(3,2) = 3._wp/5._wp

    alpha_c(1,3) = 3._wp/5._wp
    alpha_c(2,3) = 1._wp/5._wp
    alpha_c(3,3) = 1._wp/5._wp

    alpha_c(1,4) = 1._wp/5._wp
    alpha_c(2,4) = 3._wp/5._wp
    alpha_c(3,4) = 1._wp/5._wp

    ! note that the linear weighting factor is 1 (not stored)
    ptr_int%gquad%weights_tri_q(1:3) = 1._wp/3._wp
    ptr_int%gquad%weights_tri_c(1)   = -27._wp/48._wp
    ptr_int%gquad%weights_tri_c(2:4) = 25._wp/48._wp


    i_rcstartlev = 2

    ! start and end block
    i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_c

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nv,nq,i_startidx,i_endidx,ilv,ibv,z_vert_cc,z_quad_cc, &
!$OMP           z_quad_gg) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, i_rcstartlev)

      DO jc = i_startidx, i_endidx

        IF(.NOT.ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        ! loop over triangle vertices
!CDIR EXPAND=3
        DO nv=1,3
          ! get line and block indices of cell vertices
          ilv= ptr_patch%cells%vertex_idx(jc,jb,nv)
          ibv= ptr_patch%cells%vertex_blk(jc,jb,nv)

          ! Transform geographical coordinates to cartesian coordinates for vertices
          z_vert_cc(nv)=gc2cc(ptr_patch%verts%vertex(ilv,ibv))
        ENDDO

        !
        ! Linear
        !
        ! Compute quadrature point in cartesian coordinates (= triangle centroid)
        ! i.e. map area coordinates into cartesian coordinates
        z_quad_cc%x(1)= alpha_l(1)*z_vert_cc(1)%x(1)  &
          & + alpha_l(2)*z_vert_cc(2)%x(1)  &
          & + alpha_l(3)*z_vert_cc(3)%x(1)

        z_quad_cc%x(2)= alpha_l(1)*z_vert_cc(1)%x(2)  &
          & + alpha_l(2)*z_vert_cc(2)%x(2)  &
          & + alpha_l(3)*z_vert_cc(3)%x(2)

        z_quad_cc%x(3)= alpha_l(1)*z_vert_cc(1)%x(3)  &
          & + alpha_l(2)*z_vert_cc(2)%x(3)  &
          & + alpha_l(3)*z_vert_cc(3)%x(3)


        ! Transform back to geographical coordinates
        z_quad_gg = cc2gc(z_quad_cc)

        ! store
        ptr_int%gquad%qpts_tri_l(jc,jb)%lat = z_quad_gg%lat
        ptr_int%gquad%qpts_tri_l(jc,jb)%lon = z_quad_gg%lon


        !
        ! quadratic
        !
        ! Loop over quadrature points
!CDIR EXPAND=3
        DO nq=1,3
          ! map area coordinates into cartesian coordinates
          z_quad_cc%x(1)= alpha_q(1,nq)*z_vert_cc(1)%x(1)  &
            & + alpha_q(2,nq)*z_vert_cc(2)%x(1)  &
            & + alpha_q(3,nq)*z_vert_cc(3)%x(1)

          z_quad_cc%x(2)= alpha_q(1,nq)*z_vert_cc(1)%x(2)  &
            & + alpha_q(2,nq)*z_vert_cc(2)%x(2)  &
            & + alpha_q(3,nq)*z_vert_cc(3)%x(2)

          z_quad_cc%x(3)= alpha_q(1,nq)*z_vert_cc(1)%x(3)  &
            & + alpha_q(2,nq)*z_vert_cc(2)%x(3)  &
            & + alpha_q(3,nq)*z_vert_cc(3)%x(3)


          ! Transform back to geographical coordinates
          z_quad_gg = cc2gc(z_quad_cc)

          ! store
          ptr_int%gquad%qpts_tri_q(jc,jb,nq)%lat = z_quad_gg%lat
          ptr_int%gquad%qpts_tri_q(jc,jb,nq)%lon = z_quad_gg%lon
        ENDDO


        !
        ! cubic
        !
        ! Loop over quadrature points
!CDIR EXPAND=4
        DO nq=1,4
          ! map area coordinates into cartesian coordinates
          z_quad_cc%x(1)= alpha_c(1,nq)*z_vert_cc(1)%x(1)  &
            & + alpha_c(2,nq)*z_vert_cc(2)%x(1)  &
            & + alpha_c(3,nq)*z_vert_cc(3)%x(1)

          z_quad_cc%x(2)= alpha_c(1,nq)*z_vert_cc(1)%x(2)  &
            & + alpha_c(2,nq)*z_vert_cc(2)%x(2)  &
            & + alpha_c(3,nq)*z_vert_cc(3)%x(2)

          z_quad_cc%x(3)= alpha_c(1,nq)*z_vert_cc(1)%x(3)  &
            & + alpha_c(2,nq)*z_vert_cc(2)%x(3)  &
            & + alpha_c(3,nq)*z_vert_cc(3)%x(3)


          ! Transform back to geographical coordinates
          z_quad_gg = cc2gc(z_quad_cc)

          ! store
          ptr_int%gquad%qpts_tri_c(jc,jb,nq)%lat = z_quad_gg%lat
          ptr_int%gquad%qpts_tri_c(jc,jb,nq)%lon = z_quad_gg%lon
        ENDDO

      END DO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_l(:,:)%lat)
    CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_l(:,:)%lon)
    DO nq=1,3
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_q(:,:,nq)%lat)
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_q(:,:,nq)%lon)
    ENDDO
    DO nq=1,4
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_c(:,:,nq)%lat)
      CALL sync_patch_array(sync_c,ptr_patch,ptr_int%gquad%qpts_tri_c(:,:,nq)%lon)
    ENDDO

  END SUBROUTINE tri_quadrature_pts
  !-------------------------------------------------------------------------


END MODULE mo_intp_coeffs
