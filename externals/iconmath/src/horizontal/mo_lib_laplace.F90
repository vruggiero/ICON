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

!>
!!   Contains the implementation of the nabla mathematical operators.
!!
!!   Contains the implementation of the mathematical operators
!!   employed by the shallow water prototype.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_lib_laplace
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64
  USE mo_lib_interpolation_scalar, ONLY: edges2verts_scalar_lib, verts2edges_scalar_lib
  USE mo_lib_loopindices, ONLY: get_indices_c_lib, get_indices_e_lib
  USE mo_lib_divrot, ONLY: div_lib, rot_vertex_atmos_lib
  USE mo_lib_gradients, ONLY: grad_fd_norm_lib
  USE mo_fortran_tools, ONLY: copy, set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nabla2_vec_atmos_lib
  PUBLIC :: nabla2_scalar_lib, nabla2_scalar_avg_lib

!INTERFACE nabla2_vec_lib
!
!  MODULE PROCEDURE nabla2_vec_atmos_lib
!
!END INTERFACE

CONTAINS

!-------------------------------------------------------------------------
!
!>
!!  Computes  laplacian of a vector field.
!!
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P. Korn (2005).
!! Modified by Th.Heinze (2006-06-20):
!! - changed u_out(ie1,jn) to u_out(j,jn) according to hint of P.Korn
!! Modifications by P. Korn, MPI-M(2007-2)
!! -Switch from array arguments to pointers
!! Modified by P Ripodas (2007-02):
!! - include the system orientation factor in the vorticity term
!! Modified by Almut Gassmann (2007-04-20)
!! - abandon grid for the sake of patch
!!
  SUBROUTINE nabla2_vec_atmos_lib(vec_e, &
    &                             edge_cell_idx, edge_cell_blk, & ! require to calculate nabla2
    &                             edge_vertex_idx, edge_vertex_blk, & ! required to calculate nabla2
    &                             cell_edge_idx, cell_edge_blk, & ! required for div_lib
    &                             vert_edge_idx, vert_edge_blk, & ! required for rot_vertex_lib
    &                             tangent_orientation, inv_primal_edge_length, &
    &                             inv_dual_edge_length, geofac_div, geofac_rot, &
    &                             nabla2_vec_e, & ! main output vector
    &                             i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, & ! required for div_lib
    &                             i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v, & ! required for rot_vertex_lib
    &                             i_startblk, i_endblk, i_startidx_in, i_endidx_in, & ! this four are needed to call get_indices_e_lib
    &                             nlev, nblks_c, nblks_v, slev, elev, nproma, lacc)

    ! edge based variable of which laplacian is computed, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: vec_e(:, :, :)

    ! dimension: (nproma,nblks_e,2)
    INTEGER, TARGET, INTENT(IN) :: edge_cell_idx(:, :, :), edge_cell_blk(:, :, :)

    ! dimension: (nproma,nblks_e,4)
    INTEGER, TARGET, INTENT(IN) :: edge_vertex_idx(:, :, :), edge_vertex_blk(:, :, :)

    ! dimension: (nproma,nblks_c,3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_idx(:, :, :), cell_edge_blk(:, :, :)

    ! dimension: (nproma,nblks_v,6)
    INTEGER, TARGET, INTENT(IN) :: vert_edge_idx(:, :, :), vert_edge_blk(:, :, :)

    ! dimension: (nproma,nblks_e)
    REAL(wp), INTENT(IN) :: tangent_orientation(:, :)

    ! dimension: (nproma,nblks_e)
    REAL(wp), INTENT(IN) :: inv_primal_edge_length(:, :)

    ! dimension: (nproma,nblks_e)
    REAL(wp), INTENT(IN) :: inv_dual_edge_length(:, :)

    ! factor for divergence (nproma,cell_type,nblks_c)
    REAL(wp), INTENT(IN) :: geofac_div(:, :, :)

    ! factor for rotation (nproma,9-cell_type,nblks_v)
    REAL(wp), INTENT(IN) :: geofac_rot(:, :, :)

    ! edge based variable in which laplacian is stored, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(INOUT) :: nabla2_vec_e(:, :, :)

    ! required for div_lib
    INTEGER, INTENT(IN) :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c

    ! required for rot_vertex_lib
    INTEGER, INTENT(IN) :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v

    ! this four are needed to call get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

    INTEGER, INTENT(IN) :: nlev, nblks_c, nblks_v

    INTEGER, INTENT(IN) :: slev, elev

    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: je, jk, jb
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    REAL(wp) ::  &
      &  z_div_c(nproma, nlev, nblks_c),  &
      &  z_rot_v(nproma, nlev, nblks_v)

    INTEGER, DIMENSION(:, :, :), POINTER :: icidx, icblk, ividx, ivblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    icidx => edge_cell_idx
    icblk => edge_cell_blk
    ividx => edge_vertex_idx
    ivblk => edge_vertex_blk

!$ACC DATA CREATE(z_div_c, z_rot_v) PRESENT(vec_e) PRESENT(nabla2_vec_e) &
!$ACC   PRESENT(tangent_orientation, inv_primal_edge_length, inv_dual_edge_length, geofac_div, geofac_rot) &
!$ACC   PRESENT(icidx, icblk, ividx, ivblk) IF(lzacc)

! compute divergence of vector field
    CALL div_lib(vec_e, cell_edge_idx, cell_edge_blk, &
      &          geofac_div, z_div_c, &
      &          i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
      &          slev, elev, nproma)

!
!  loop through over all patch edges (and blocks)
!

! The special treatment of 2D fields is essential for efficiency on the NEC

    ! compute rotation of vector field
    CALL rot_vertex_atmos_lib(vec_e, vert_edge_idx, vert_edge_blk, &
      &                       geofac_rot, z_rot_v, &
      &                       i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v, &
      &                       slev, elev, nproma)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO je = i_startidx, i_endidx
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
!CDIR UNROLL=3
      !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO je = i_startidx, i_endidx
#endif

          nabla2_vec_e(je, jk, jb) =  &
            &   tangent_orientation(je, jb)*  &
            &   (z_rot_v(ividx(je, jb, 2), jk, ivblk(je, jb, 2))  &
            &   - z_rot_v(ividx(je, jb, 1), jk, ivblk(je, jb, 1)))  &
            &   *inv_primal_edge_length(je, jb)  &
            & + (z_div_c(icidx(je, jb, 2), jk, icblk(je, jb, 2))    &
            &   - z_div_c(icidx(je, jb, 1), jk, icblk(je, jb, 1)))  &
            &   *inv_dual_edge_length(je, jb)

        END DO
      END DO
      !$ACC END PARALLEL

    END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC WAIT(1)
!$ACC END DATA

  END SUBROUTINE nabla2_vec_atmos_lib
!-----------------------------------------------------------------------
!>
!!  Computes laplacian @f$\nabla ^2 @f$ of a scalar field.
!!
!! input:  lives on cells (mass points)
!! output: lives on cells
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P.Korn (2005).
!! Derived from nabla4_scalar by Th.Heinze (2006-07-24).
!! Modifications by P. Korn, MPI-M(2007-2)
!! -Switch from array arguments to pointers
!! Modifications by Almut Gassmann, MPI-M (2007-04-20)
!! -abandon grid for the sake of patch
!!
  SUBROUTINE nabla2_scalar_lib(psi_c, cell_neighbor_idx, cell_neighbor_blk, &
    &                          edge_cell_idx, edge_cell_blk, inv_dual_edge_length, & ! required for grad_fd_norm_lib
    &                          cell_edge_idx, cell_edge_blk, & ! required for div_lib
    &                          geofac_n2s, geofac_div, nabla2_psi_c, &
    &                          i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                          i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
    &                          nlev, slev, elev, nproma, nblks_e, cell_type, lacc)

    ! cells based variable of which biharmonic laplacian is computed, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(IN) :: psi_c(:, :, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c,3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_idx(:, :, :)

    ! block indices of triangles next to each cell, dim: (nproma,nblks_c,3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_blk(:, :, :)

    ! line indices of adjacent cells, dim: (nproma,nblks_e,2)
    INTEGER, TARGET, INTENT(IN) :: edge_cell_idx(:, :, :)

    ! block indices of adjacent cells, dim: (nproma,nblks_e,2)
    INTEGER, TARGET, INTENT(IN) :: edge_cell_blk(:, :, :)

    ! inverse of dual edge length, dim: (nproma,nblks_e)
    REAL(wp), INTENT(IN) :: inv_dual_edge_length(:, :)

    ! line indices of edges of triangles, dim: (nproma,nblks_c,3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_idx(:, :, :)

    ! block indices of edges of triangles, dim: (nproma,nblks_c,3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_blk(:, :, :)

    ! factor for nabla2-scalar (nproma,cell_type+1,nblks_c)
    REAL(wp), INTENT(IN) :: geofac_n2s(:, :, :)

    ! factor for divergence (nproma,cell_type,nblks_c)
    REAL(wp), INTENT(IN) :: geofac_div(:, :, :)

    ! cell based variable in which biharmonic laplacian is stored, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) :: nabla2_psi_c(:, :, :)


    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! start_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startblk_e

    ! end_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endblk_e

    ! start_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startidx_e

    ! end_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endidx_e

    ! number of vertical levels
    INTEGER, INTENT(IN) :: nlev

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! number of blocks for the edges: t_patch%nblks_e
    INTEGER, INTENT(IN) :: nblks_e

    ! the grid domain geometry parameter, cell_type = 3 or 6
    INTEGER, INTENT(IN) :: cell_type

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jb, jc, jk, i_startidx, i_endidx
    REAL(wp) :: z_grad_fd_norm_e(nproma, nlev, nblks_e)
    LOGICAL :: lzacc
    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk

!$ACC DATA CREATE(z_grad_fd_norm_e) PRESENT(psi_c) PRESENT(nabla2_psi_c) &
!$ACC   PRESENT(iidx, iblk, inv_dual_edge_length, geofac_n2s) IF(lzacc)

    SELECT CASE (cell_type)

    CASE (3) ! (cell_type == 3)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
 
        CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx)
 
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
 
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = slev, elev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
#endif
            !
            !  calculate div(grad) in one step
            !
            nabla2_psi_c(jc, jk, jb) =  &
              &    psi_c(jc, jk, jb)*geofac_n2s(jc, 1, jb) &
              &  + psi_c(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))*geofac_n2s(jc, 2, jb) &
              &  + psi_c(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))*geofac_n2s(jc, 3, jb) &
              &  + psi_c(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))*geofac_n2s(jc, 4, jb)
          END DO !cell loop
        END DO !vertical levels loop
        !$ACC END PARALLEL

      END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CASE (6) ! (cell_type == 6) THEN ! Use unoptimized version for the time being

    ! compute finite difference gradient in normal direction
      CALL grad_fd_norm_lib(psi_c, edge_cell_idx, edge_cell_blk, &
        &                   inv_dual_edge_length, z_grad_fd_norm_e, &
        &                   i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
        &                   slev, elev, nproma)
    
      ! compute divergence of resulting vector field
      CALL div_lib(z_grad_fd_norm_e, cell_edge_idx, cell_edge_blk, &
        &          geofac_div, nabla2_psi_c, &
        &          i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
        &          slev, elev, nproma)
  
    END SELECT
!$ACC WAIT(1)

!$ACC END DATA

  END SUBROUTINE nabla2_scalar_lib

!-------------------------------------------------------------------------
!

!>
!!  Computes Laplacian @f$\nabla ^2 @f$ of a scalar field, followed by weighted averaging.
!!
!!  Computes Laplacian @f$\nabla ^2 @f$ of a scalar field, followed by weighted averaging
!!  with the neighboring cells to increase computing efficiency.
!!  NOTE: This optimized routine works for triangular grids only.
!! input:  lives on cells (mass points)
!! output: lives on cells
!!
!! @par Revision History
!! Developed by Guenther Zaengl, DWD, 2009-05-19
!!
  SUBROUTINE nabla2_scalar_avg_lib(psi_c, cell_neighbor_idx, cell_neighbor_blk, &
    &                              geofac_n2s, avg_coeff, nabla2_psi_c, &
    &                              i_startblk_in, i_endblk_in, i_startidx_in, i_endidx_in, &
    &                              nblks_c, cell_type, patch_id, &
    &                              nlev, slev, elev, nproma, l_limited_area, lacc)

    ! cells based variable of which biharmonic laplacian is computed, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(IN) :: psi_c(:, :, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c,3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_idx(:, :, :)

    ! block indices of triangles next to each cell, dim: (nproma,nblks_c,3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_blk(:, :, :)

    ! factor for nabla2-scalar (nproma,cell_type+1,nblks_c)
    REAL(wp), INTENT(IN) :: geofac_n2s(:, :, :)

    ! averaging coefficients dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(IN) :: avg_coeff(:, :, :)

    ! cell based variable in which biharmonic laplacian is stored, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) :: nabla2_psi_c(:, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk_in(3)

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk_in(3)

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in(3)

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in(3)

    ! number of blocks for the cells: t_patch%nblks_c
    INTEGER, INTENT(IN) :: nblks_c

    ! the grid domain geometry parameter, cell_type = 3 or 6
    INTEGER, INTENT(IN) :: cell_type

    ! domain ID of current domain
    INTEGER, INTENT(IN) :: patch_id

    ! number of vertical levels
    INTEGER, INTENT(IN) :: nlev

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! limited area setup where forcing comes in from sides
    LOGICAL, INTENT(IN) :: l_limited_area

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jb, jc, jk, i_startblk, i_endblk, i_startidx, i_endidx
    LOGICAL :: lzacc

    REAL(wp), DIMENSION(nproma, nlev, nblks_c) :: aux_c

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk

! The special treatment of 2D fields is essential for efficiency on the NEC

!$ACC DATA CREATE(aux_c) PRESENT(avg_coeff, psi_c) PRESENT(nabla2_psi_c) &
!$ACC   PRESENT(iidx, iblk, geofac_n2s) IF(lzacc)

    SELECT CASE (cell_type)

    CASE (3) ! (cell_type == 3)

      IF (slev == elev) THEN
        jk = slev

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! values for the blocking
        i_startblk = i_startblk_in(1)
        i_endblk = i_endblk_in(1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c_lib(i_startidx_in(1), i_endidx_in(1), nproma, jb, i_startblk, i_endblk, &
                                 i_startidx, i_endidx)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx

            !
            !  calculate div(grad) in one step
            !
            aux_c(jc, jk, jb) =  &
              &    psi_c(jc, jk, jb)*geofac_n2s(jc, 1, jb) &
              &  + psi_c(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))*geofac_n2s(jc, 2, jb) &
              &  + psi_c(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))*geofac_n2s(jc, 3, jb) &
              &  + psi_c(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))*geofac_n2s(jc, 4, jb)

          END DO
          !$ACC END PARALLEL

        END DO
!$OMP END DO

        IF (l_limited_area .OR. patch_id > 1) THEN
          ! Fill nabla2_psi_c along the lateral boundaries of nests

          i_startblk = i_startblk_in(2)
          i_endblk = i_endblk_in(2)

          CALL copy(aux_c(:, jk, i_startblk:i_endblk), &
                    nabla2_psi_c(:, jk, i_startblk:i_endblk), lacc=lzacc)
!$OMP BARRIER
        END IF

!
! Now do averaging with weights given by avg_coeff

        ! values for the blocking
        i_startblk = i_startblk_in(3)
        i_endblk = i_endblk_in(3)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c_lib(i_startidx_in(3), i_endidx_in(3), nproma, jb, i_startblk, i_endblk, &
                                 i_startidx, i_endidx)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP VECTOR
          DO jc = i_startidx, i_endidx
            !
            !  calculate the weighted average
            !
            nabla2_psi_c(jc, jk, jb) =  &
              &    aux_c(jc, jk, jb)*avg_coeff(jc, 1, jb) &
              &  + aux_c(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))*avg_coeff(jc, 2, jb) &
              &  + aux_c(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))*avg_coeff(jc, 3, jb) &
              &  + aux_c(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))*avg_coeff(jc, 4, jb)

          END DO !cell loop
          !$ACC END PARALLEL

        END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ELSE

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! values for the blocking
        i_startblk = i_startblk_in(1)
        i_endblk = i_endblk_in(1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c_lib(i_startidx_in(1), i_endidx_in(1), nproma, jb, i_startblk, i_endblk, &
                                 i_startidx, i_endidx)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            !$ACC LOOP VECTOR
            DO jk = slev, elev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
          DO jk = slev, elev
            !$ACC LOOP VECTOR
            DO jc = i_startidx, i_endidx
#endif
              !
              !  calculate div(grad) in one step
              !
              aux_c(jc, jk, jb) =  &
                &    psi_c(jc, jk, jb)*geofac_n2s(jc, 1, jb) &
                &  + psi_c(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))*geofac_n2s(jc, 2, jb) &
                &  + psi_c(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))*geofac_n2s(jc, 3, jb) &
                &  + psi_c(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))*geofac_n2s(jc, 4, jb)
            END DO !cell loop
          END DO !vertical levels loop
          !$ACC END PARALLEL

        END DO
!$OMP END DO

        IF (l_limited_area .OR. patch_id > 1) THEN
          ! Fill nabla2_psi_c along the lateral boundaries of nests

          i_startblk = i_startblk_in(2)
          i_endblk = i_endblk_in(2)

          CALL copy(aux_c(:, :, i_startblk:i_endblk), &
                    nabla2_psi_c(:, :, i_startblk:i_endblk), lacc=lzacc)
!$OMP BARRIER
        END IF

!
! Now do averaging with weights given by avg_coeff

         ! values for the blocking
        i_startblk = i_startblk_in(3)
        i_endblk = i_endblk_in(3)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c_lib(i_startidx_in(3), i_endidx_in(3), nproma, jb, i_startblk, i_endblk, &
                                 i_startidx, i_endidx)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            !$ACC LOOP VECTOR
            DO jk = slev, elev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
          DO jk = slev, elev
            !$ACC LOOP VECTOR
            DO jc = i_startidx, i_endidx
#endif
              !
              !  calculate the weighted average
              !
              nabla2_psi_c(jc, jk, jb) =  &
                &    aux_c(jc, jk, jb)*avg_coeff(jc, 1, jb) &
                &  + aux_c(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))*avg_coeff(jc, 2, jb) &
                &  + aux_c(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))*avg_coeff(jc, 3, jb) &
                &  + aux_c(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))*avg_coeff(jc, 4, jb)

            END DO !cell loop
          END DO !vertical levels loop
          !$ACC END PARALLEL

        END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      END IF

    END SELECT
!$ACC WAIT(1)

!$ACC END DATA

  END SUBROUTINE nabla2_scalar_avg_lib

END MODULE mo_lib_laplace
