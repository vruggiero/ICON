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
!!   Contains the implementation of the mathematical grad operators.
!!
!!   Contains the implementation of the mathematical operators
!!   employed by the shallow water prototype.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_lib_gradients
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64, &
                                           sp => real32
  USE mo_lib_loopindices, ONLY: get_indices_c_lib, get_indices_e_lib
  USE mo_fortran_tools, ONLY: init, set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: grad_fd_norm_lib
  PUBLIC :: grad_fd_tang_lib
  PUBLIC :: grad_green_gauss_cell_lib
  PUBLIC :: grad_fe_cell_lib

#ifdef __MIXED_PRECISION
  INTEGER, PARAMETER :: vp = sp
#else
  INTEGER, PARAMETER :: vp = wp
#endif

  INTERFACE grad_green_gauss_cell_lib
    MODULE PROCEDURE grad_green_gauss_cell_adv_lib
    MODULE PROCEDURE grad_green_gauss_cell_dycore_lib
  END INTERFACE
!
  INTERFACE grad_fe_cell_lib
    MODULE PROCEDURE grad_fe_cell_adv_lib
    MODULE PROCEDURE grad_fe_cell_adv_2d_lib
    MODULE PROCEDURE grad_fe_cell_dycore_lib
  END INTERFACE

CONTAINS

!-------------------------------------------------------------------------
!

!-------------------------------------------------------------------------
!
!>
!!  Computes directional  derivative of a cell centered variable.
!!
!!  Computes directional  derivative of a cell centered variable
!!  with respect to direction normal to triangle edge.
!! input: lives on centres of triangles
!! output:  lives on edges (velocity points)
!!
  SUBROUTINE grad_fd_norm_lib(psi_c, edge_cell_idx, edge_cell_blk, &
    &                         inv_dual_edge_length, grad_norm_psi_e, &
    &                         i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                         slev, elev, nproma, lacc)

    ! cell based variable of which normal derivative is computed, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(in) :: psi_c(:, :, :)

    ! line indices of triangles next to each edge, dim: (nproma,nblks_e, 2)
    INTEGER, TARGET, INTENT(IN) :: edge_cell_idx(:, :, :)

    ! block indices of triangles next to each edge, dim: (nproma,nblks_e, 2)
    INTEGER, TARGET, INTENT(IN) :: edge_cell_blk(:, :, :)
!
    ! inverse of dual edge length, dim: (nproma, nblks_e)
    REAL(wp), INTENT(IN) :: inv_dual_edge_length(:, :)

    ! edge based variable in which normal derivative is stored, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(inout) :: grad_norm_psi_e(:, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: je, jk, jb
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!
!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => edge_cell_idx
    iblk => edge_cell_blk

!
!  loop through all patch edges (and blocks)
!

    !$ACC DATA PRESENT(psi_c, grad_norm_psi_e, inv_dual_edge_length, iidx, iblk) IF(lzacc)

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
          !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO je = i_startidx, i_endidx
#endif
          !
          ! compute the normal derivative
          ! by the finite difference approximation
          ! (see Bonaventura and Ringler MWR 2005)
          !
          grad_norm_psi_e(je, jk, jb) =  &
             &  (psi_c(iidx(je, jb, 2), jk, iblk(je, jb, 2)) - &
             &    psi_c(iidx(je, jb, 1), jk, iblk(je, jb, 1)))  &
             &  *inv_dual_edge_length(je, jb)

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE grad_fd_norm_lib

!-------------------------------------------------------------------------
!
! RESTRUCT: @Marco: please adjust calls to this routine to your needs.
!>
!! Computes directional derivative of a vertex centered variable with.
!!
!! Computes directional derivative of a vertex centered variable with
!! respect to direction tangent to triangle edge. Notice that the
!! tangential direction is defined by
!!   iorient*(vertex2 - vertex1)
!! input: lives on vertices of triangles
!! output: lives on edges (velocity points)
!!
  SUBROUTINE grad_fd_tang_lib(psi_v, edge_vertex_idx, edge_vertex_blk, &
    &                          primal_edge_length, tangent_orientation, grad_tang_psi_e, &
    &                          i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                          slev, elev, nproma, lacc)

    ! vertex based variable of which tangential derivative is computed, dim: (nproma,nlev,nblks_v)
    REAL(wp), INTENT(in) ::  psi_v(:, :, :)

    ! line indices of vertices next to each edge, dim: (nproma,nblks_e, 2)
    INTEGER, TARGET, INTENT(IN) :: edge_vertex_idx(:, :, :)

    ! block indices of vertices next to each edge, dim: (nproma,nblks_e, 2)
    INTEGER, TARGET, INTENT(IN) :: edge_vertex_blk(:, :, :)

    ! primal edge length, dim: (nproma, nblks_e)
    REAL(wp), INTENT(IN) :: primal_edge_length(:, :)

    ! tangent_orientation, dim: (nproma, nblks_e)
    REAL(wp), INTENT(IN) :: tangent_orientation(:, :)

    ! edge based variable in which tangential derivative is stored, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(inout) :: grad_tang_psi_e(:, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp) :: iorient

    INTEGER :: je, jk, jb
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(nproma) ::  &
      &  ilv1, ibv1, ilv2, ibv2
!
!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

!$ACC DATA PRESENT(psi_v, grad_tang_psi_e, edge_vertex_idx, edge_vertex_blk, tangent_orientation, primal_edge_length) CREATE(ilv1, ibv1, ilv2, ibv2) IF(lzacc)

!
! TODO: OpenMP
!
!
!  loop through all patch edges (and blocks)
!
    DO jb = i_startblk, i_endblk

      CALL get_indices_e_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO je = i_startidx, i_endidx
        !
        !  get the line and block indices of the vertices of edge je
        !
        ilv1(je) = edge_vertex_idx(je, jb, 1)
        ibv1(je) = edge_vertex_blk(je, jb, 1)
        ilv2(je) = edge_vertex_idx(je, jb, 2)
        ibv2(je) = edge_vertex_blk(je, jb, 2)
      END DO

      DO jk = slev, elev

        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(iorient)
        DO je = i_startidx, i_endidx
          !
          ! compute the tangential derivative
          ! by the finite difference approximation
          iorient = tangent_orientation(je, jb)
          grad_tang_psi_e(je, jk, jb) = iorient  &
            &  *(psi_v(ilv2(je), jk, ibv2(je)) - psi_v(ilv1(je), jk, ibv1(je)))  &
            &    /primal_edge_length(je, jb)
        END DO

      END DO
      !$ACC END PARALLEL

    END DO
!
! TODO: OpenMP
!

!$ACC WAIT(1)
!$ACC END DATA

  END SUBROUTINE grad_fd_tang_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons
!!
!!
  SUBROUTINE grad_fe_cell_adv_lib(p_cc, cell_neighbor_idx, cell_neighbor_blk, &
    &                              gradc_bmat, p_grad, &
    &                              i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                              slev, elev, nproma, patch_id, lacc)

    REAL(wp), INTENT(in) :: p_cc(:, :, :) !  cell centered variable

    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_idx(:, :, :) ! line indices of triangles next to each cell
    ! dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_blk(:, :, :) ! block indices of triangles next to each cell
    ! dim: (nproma,nblks_c, 3)
    REAL(wp), INTENT(in) ::   gradc_bmat(:, :, :, :) ! Bmatrix for cell centered shape function based
    ! gradient (nproma,2,3,nblks_c)

    REAL(vp), INTENT(inout) :: p_grad(:, :, :, :) ! cell based Green-Gauss reconstructed geographical gradient vector
    ! dim:(2,nproma,nlev,nblks_c)

    INTEGER, INTENT(IN) :: i_startblk ! start_block needed for get_indices_c_lib

    INTEGER, INTENT(IN) :: i_endblk ! end_block needed for get_indices_c_lib

    INTEGER, INTENT(IN) :: i_startidx_in ! start_index needed for get_indices_c_lib

    INTEGER, INTENT(IN) :: i_endidx_in ! end_index needed for get_indices_c_lib

    INTEGER, INTENT(IN) :: slev ! vertical start level

    INTEGER, INTENT(IN) :: elev ! vertical end level

    INTEGER, INTENT(IN) :: nproma ! inner loop length/vector length

    INTEGER, INTENT(IN) :: patch_id ! domain ID of current domain

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! if true, use OpenAcc

    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk

!
! 2. reconstruction of cell based geographical gradient
!

    !$ACC DATA PRESENT(p_cc, p_grad, gradc_bmat, iidx, iblk) IF(lzacc)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    IF (patch_id > 1) THEN
      ! Fill nest boundaries with zero to avoid trouble with MPI synchronization

#ifdef _OPENACC
      !$ACC KERNELS PRESENT(p_grad) ASYNC(1) IF(lzacc)
      p_grad(:, :, :, 1:i_startblk) = 0._wp
      !$ACC END KERNELS
#else
      CALL init(p_grad(:, :, :, 1:i_startblk), lacc=lzacc)
!$OMP BARRIER
#endif
    END IF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          ! We do not make use of the intrinsic function DOT_PRODUCT on purpose,
          ! since it is extremely slow on the SX9, when combined with indirect
          ! addressing.

          ! multiply cell-based input values with precomputed grid geometry factor

          ! zonal(u)-component of Green-Gauss gradient
          p_grad(1, jc, jk, jb) = &
            &    gradc_bmat(jc, 1, 1, jb)*p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))  &
            &  + gradc_bmat(jc, 1, 2, jb)*p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))  &
            &  + gradc_bmat(jc, 1, 3, jb)*p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

          ! meridional(v)-component of Green-Gauss gradient
          p_grad(2, jc, jk, jb) =  &
            &    gradc_bmat(jc, 2, 1, jb)*p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))  &
            &  + gradc_bmat(jc, 2, 2, jb)*p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))  &
            &  + gradc_bmat(jc, 2, 3, jb)*p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE grad_fe_cell_adv_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!! 2D version, i.e. for a single vertical level
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons
!!
!!
  SUBROUTINE grad_fe_cell_adv_2d_lib(p_cc, cell_neighbor_idx, cell_neighbor_blk, &
    &                                gradc_bmat, p_grad, &
    &                                i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                nproma, patch_id, lacc)

    ! cell centered variable, dim: (nproma,nblks_c)
    REAL(wp), INTENT(in) :: p_cc(:, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_idx(:, :, :)

    ! block indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_blk(:, :, :)

    ! Bmatrix for cell centered shape function based gradient (nproma,2,3,nblks_c)
    REAL(wp), INTENT(in) :: gradc_bmat(:, :, :, :)

    ! cell based Green-Gauss reconstructed geographical gradient vector, dim:(2,nproma,nblks_c)
    REAL(wp), INTENT(inout) :: p_grad(:, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! domain ID of current domain
    INTEGER, INTENT(IN) :: patch_id

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jc, jb
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk

!
! 2. reconstruction of cell based geographical gradient
!
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    IF (patch_id > 1) THEN
      ! Fill nest boundaries with zero to avoid trouble with MPI synchronization
      CALL init(p_grad(:, :, 1:i_startblk), lacc=lzacc)
!$OMP BARRIER
    END IF

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      DO jc = i_startidx, i_endidx

        ! We do not make use of the intrinsic function DOT_PRODUCT on purpose,
        ! since it is extremely slow on the SX9, when combined with indirect
        ! addressing.

        ! multiply cell-based input values with precomputed grid geometry factor

        ! zonal(u)-component of gradient
        p_grad(1, jc, jb) = &
          &    gradc_bmat(jc, 1, 1, jb)*p_cc(iidx(jc, jb, 1), iblk(jc, jb, 1))  &
          &  + gradc_bmat(jc, 1, 2, jb)*p_cc(iidx(jc, jb, 2), iblk(jc, jb, 2))  &
          &  + gradc_bmat(jc, 1, 3, jb)*p_cc(iidx(jc, jb, 3), iblk(jc, jb, 3))

        ! meridional(v)-component of gradient
        p_grad(2, jc, jb) =  &
          &    gradc_bmat(jc, 2, 1, jb)*p_cc(iidx(jc, jb, 1), iblk(jc, jb, 1))  &
          &  + gradc_bmat(jc, 2, 2, jb)*p_cc(iidx(jc, jb, 2), iblk(jc, jb, 2))  &
          &  + gradc_bmat(jc, 2, 3, jb)*p_cc(iidx(jc, jb, 3), iblk(jc, jb, 3))

      END DO ! end loop over cells

    END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE grad_fe_cell_adv_2d_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!! Special dycore version, which handles two fields at a time.
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons
!!
!!
  SUBROUTINE grad_fe_cell_dycore_lib(p_ccpr, cell_neighbor_idx, cell_neighbor_blk, &
    &                                gradc_bmat, p_grad, &
    &                                i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                slev, elev, nproma, lacc)

    ! perturbation fields passed from dycore, dim: (nproma,2,nlev,nblks_c)
    REAL(vp), INTENT(in) :: p_ccpr(:, :, :, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN)  :: cell_neighbor_idx(:, :, :)

    ! block indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN)  :: cell_neighbor_blk(:, :, :)

    ! Bmatrix for cell centered shape function based gradient (nproma,2,3,nblks_c)
    REAL(wp), TARGET, INTENT(in) :: gradc_bmat(:, :, :, :)

    ! cell based Green-Gauss reconstructed geographical gradient vector, dim:(2,nproma,nlev,nblks_c)
    REAL(vp), INTENT(inout) :: p_grad(:, :, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk
!
! 2. reconstruction of cell based geographical gradient
!

    !$ACC DATA PRESENT(p_ccpr, p_grad, gradc_bmat, iidx, iblk) IF(lzacc)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          ! We do not make use of the intrinsic function DOT_PRODUCT on purpose,
          ! since it is extremely slow on the SX9, when combined with indirect
          ! addressing.

          ! multiply cell-based input values with shape function derivatives

          ! zonal(u)-component of gradient, field 1
          p_grad(1, jc, jk, jb) = &
            &    gradc_bmat(jc, 1, 1, jb)*p_ccpr(1, iidx(jc, jb, 1), jk, iblk(jc, jb, 1))  &
            &  + gradc_bmat(jc, 1, 2, jb)*p_ccpr(1, iidx(jc, jb, 2), jk, iblk(jc, jb, 2))  &
            &  + gradc_bmat(jc, 1, 3, jb)*p_ccpr(1, iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

          ! meridional(v)-component of gradient, field 1
          p_grad(2, jc, jk, jb) =  &
            &    gradc_bmat(jc, 2, 1, jb)*p_ccpr(1, iidx(jc, jb, 1), jk, iblk(jc, jb, 1))  &
            &  + gradc_bmat(jc, 2, 2, jb)*p_ccpr(1, iidx(jc, jb, 2), jk, iblk(jc, jb, 2))  &
            &  + gradc_bmat(jc, 2, 3, jb)*p_ccpr(1, iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

          ! zonal(u)-component of gradient, field 2
          p_grad(3, jc, jk, jb) = &
            &    gradc_bmat(jc, 1, 1, jb)*p_ccpr(2, iidx(jc, jb, 1), jk, iblk(jc, jb, 1))  &
            &  + gradc_bmat(jc, 1, 2, jb)*p_ccpr(2, iidx(jc, jb, 2), jk, iblk(jc, jb, 2))  &
            &  + gradc_bmat(jc, 1, 3, jb)*p_ccpr(2, iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

          ! meridional(v)-component of gradient, field 2
          p_grad(4, jc, jk, jb) =  &
            &    gradc_bmat(jc, 2, 1, jb)*p_ccpr(2, iidx(jc, jb, 1), jk, iblk(jc, jb, 1))  &
            &  + gradc_bmat(jc, 2, 2, jb)*p_ccpr(2, iidx(jc, jb, 2), jk, iblk(jc, jb, 2))  &
            &  + gradc_bmat(jc, 2, 3, jb)*p_ccpr(2, iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE grad_fe_cell_dycore_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The Green-Gauss approach is used. See for example:
!! http://www.cfd-online.com/Wiki/Gradient_computation
!!
  SUBROUTINE grad_green_gauss_cell_adv_lib(p_cc, cell_neighbor_idx, cell_neighbor_blk, &
    &                                       geofac_grg, p_grad, &
    &                                       i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                       slev, elev, nproma, patch_id, lacc)

    ! cell centered variable, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(in) :: p_cc(:, :, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_idx(:, :, :)
    ! dim: (nproma,nblks_c, 3)
    ! block indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_blk(:, :, :)
    ! dim: (nproma,nblks_c, 3)
    ! factor for Green-Gauss gradient (nproma,4,nblks_c,2)
    REAL(wp), INTENT(in) ::  geofac_grg(:, :, :, :)

    ! cell based Green-Gauss reconstructed geographical gradient vector, dim:(2,nproma,nlev,nblks_c)
    REAL(vp), INTENT(inout) :: p_grad(:, :, :, :)
    ! dim:(2,nproma,nlev,nblks_c)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk ! end_block needed for get_indices_c_lib

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in ! start_index needed for get_indices_c_lib

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in ! end_index needed for get_indices_c_lib

    ! vertical start level
    INTEGER, INTENT(IN) :: slev ! vertical start level

    ! vertical end level
    INTEGER, INTENT(IN) :: elev ! vertical end level

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma ! inner loop length/vector length

    ! domain ID of current domain
    INTEGER, INTENT(IN) :: patch_id ! domain ID of current domain

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! if true, use OpenAcc

    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk

!
! 2. reconstruction of cell based geographical gradient
!
    !$ACC DATA PRESENT(p_cc, p_grad, geofac_grg, iidx, iblk) IF(lzacc)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    IF (patch_id > 1) THEN
      ! Fill nest boundaries with zero to avoid trouble with MPI synchronization
#ifdef _OPENACC
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      p_grad(:, :, :, 1:i_startblk) = 0._wp
      !$ACC END KERNELS
#else
      CALL init(p_grad(:, :, :, 1:i_startblk), lacc=lzacc)
!$OMP BARRIER
#endif
    END IF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          ! multiply cell-based input values with precomputed grid geometry factor

          ! zonal(u)-component of Green-Gauss gradient
          p_grad(1, jc, jk, jb) = geofac_grg(jc, 1, jb, 1)*p_cc(jc, jk, jb) + &
            geofac_grg(jc, 2, jb, 1)*p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            geofac_grg(jc, 3, jb, 1)*p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            geofac_grg(jc, 4, jb, 1)*p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

          ! meridional(v)-component of Green-Gauss gradient
          p_grad(2, jc, jk, jb) = geofac_grg(jc, 1, jb, 2)*p_cc(jc, jk, jb) + &
            geofac_grg(jc, 2, jb, 2)*p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            geofac_grg(jc, 3, jb, 2)*p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            geofac_grg(jc, 4, jb, 2)*p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE grad_green_gauss_cell_adv_lib

  SUBROUTINE grad_green_gauss_cell_dycore_lib(p_ccpr, cell_neighbor_idx, cell_neighbor_blk, &
    &                                         geofac_grg, p_grad, &
    &                                         i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                         slev, elev, nproma, lacc, acc_async)

    !
    !  cell centered I/O variables
    !
    ! perturbation fields passed from dycore, dim: (nproma,2,nlev,nblks_c)
    REAL(vp), INTENT(in) :: p_ccpr(:, :, :, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_idx(:, :, :)

    ! block indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_blk(:, :, :)

    ! factor for Green-Gauss gradient (nproma,4,nblks_c,2)
    REAL(wp), INTENT(in) :: geofac_grg(:, :, :, :)

    ! cell based Green-Gauss reconstructed geographical gradient vector, dim:(2,nproma,nlev,nblks_c)
    REAL(vp), INTENT(inout) :: p_grad(:, :, :, :)
    ! dim:(4,nproma,nlev,nblks_c)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! async OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: acc_async

    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk

    !
    ! 2. reconstruction of cell based geographical gradient
    !

    !$ACC DATA PRESENT(p_ccpr, p_grad, geofac_grg, iidx, iblk) IF(lzacc)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else

      !$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(8)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

#ifdef __SWAPDIM
          ! zonal(u)-component of Green-Gauss gradient, field 1
          p_grad(jc, jk, jb, 1) = &
            geofac_grg(jc, 1, jb, 1)*p_ccpr(jc, jk, jb, 1) + &
            geofac_grg(jc, 2, jb, 1)*p_ccpr(iidx(jc, jb, 1), jk, iblk(jc, jb, 1), 1) + &
            geofac_grg(jc, 3, jb, 1)*p_ccpr(iidx(jc, jb, 2), jk, iblk(jc, jb, 2), 1) + &
            geofac_grg(jc, 4, jb, 1)*p_ccpr(iidx(jc, jb, 3), jk, iblk(jc, jb, 3), 1)
          ! meridional(v)-component of Green-Gauss gradient, field 1
          p_grad(jc, jk, jb, 2) = &
            geofac_grg(jc, 1, jb, 2)*p_ccpr(jc, jk, jb, 1) + &
            geofac_grg(jc, 2, jb, 2)*p_ccpr(iidx(jc, jb, 1), jk, iblk(jc, jb, 1), 1) + &
            geofac_grg(jc, 3, jb, 2)*p_ccpr(iidx(jc, jb, 2), jk, iblk(jc, jb, 2), 1) + &
            geofac_grg(jc, 4, jb, 2)*p_ccpr(iidx(jc, jb, 3), jk, iblk(jc, jb, 3), 1)
          ! zonal(u)-component of Green-Gauss gradient, field 2
          p_grad(jc, jk, jb, 3) = &
            geofac_grg(jc, 1, jb, 1)*p_ccpr(jc, jk, jb, 2) + &
            geofac_grg(jc, 2, jb, 1)*p_ccpr(iidx(jc, jb, 1), jk, iblk(jc, jb, 1), 2) + &
            geofac_grg(jc, 3, jb, 1)*p_ccpr(iidx(jc, jb, 2), jk, iblk(jc, jb, 2), 2) + &
            geofac_grg(jc, 4, jb, 1)*p_ccpr(iidx(jc, jb, 3), jk, iblk(jc, jb, 3), 2)
          ! meridional(v)-component of Green-Gauss gradient, field 2
          p_grad(jc, jk, jb, 4) = &
            geofac_grg(jc, 1, jb, 2)*p_ccpr(jc, jk, jb, 2) + &
            geofac_grg(jc, 2, jb, 2)*p_ccpr(iidx(jc, jb, 1), jk, iblk(jc, jb, 1), 2) + &
            geofac_grg(jc, 3, jb, 2)*p_ccpr(iidx(jc, jb, 2), jk, iblk(jc, jb, 2), 2) + &
            geofac_grg(jc, 4, jb, 2)*p_ccpr(iidx(jc, jb, 3), jk, iblk(jc, jb, 3), 2)
#else
          ! zonal(u)-component of Green-Gauss gradient, field 1
          p_grad(1, jc, jk, jb) = geofac_grg(jc, 1, jb, 1)*p_ccpr(1, jc, jk, jb) + &
            geofac_grg(jc, 2, jb, 1)*p_ccpr(1, iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            geofac_grg(jc, 3, jb, 1)*p_ccpr(1, iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            geofac_grg(jc, 4, jb, 1)*p_ccpr(1, iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

          ! meridional(v)-component of Green-Gauss gradient, field 1
          p_grad(2, jc, jk, jb) = geofac_grg(jc, 1, jb, 2)*p_ccpr(1, jc, jk, jb) + &
            geofac_grg(jc, 2, jb, 2)*p_ccpr(1, iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            geofac_grg(jc, 3, jb, 2)*p_ccpr(1, iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            geofac_grg(jc, 4, jb, 2)*p_ccpr(1, iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

          ! zonal(u)-component of Green-Gauss gradient, field 2
          p_grad(3, jc, jk, jb) = geofac_grg(jc, 1, jb, 1)*p_ccpr(2, jc, jk, jb) + &
            geofac_grg(jc, 2, jb, 1)*p_ccpr(2, iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            geofac_grg(jc, 3, jb, 1)*p_ccpr(2, iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            geofac_grg(jc, 4, jb, 1)*p_ccpr(2, iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

          ! meridional(v)-component of Green-Gauss gradient, field 2
          p_grad(4, jc, jk, jb) = geofac_grg(jc, 1, jb, 2)*p_ccpr(2, jc, jk, jb) + &
            geofac_grg(jc, 2, jb, 2)*p_ccpr(2, iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            geofac_grg(jc, 3, jb, 2)*p_ccpr(2, iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            geofac_grg(jc, 4, jb, 2)*p_ccpr(2, iidx(jc, jb, 3), jk, iblk(jc, jb, 3))
#endif
        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    END DO ! end loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (PRESENT(acc_async)) THEN
      IF (.NOT. acc_async) THEN
        !$ACC WAIT
      END IF
    ELSE
      !$ACC WAIT
    END IF
    !$ACC END DATA
  END SUBROUTINE grad_green_gauss_cell_dycore_lib

END MODULE mo_lib_gradients
