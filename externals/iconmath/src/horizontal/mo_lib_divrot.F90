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
!!   Contains the implementation of the div,rot,recon mathematical operators.
!!
!!   Contains the implementation of the mathematical operators
!!   employed by the shallow water prototype.
!!
MODULE mo_lib_divrot
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
  USE mo_lib_loopindices, ONLY: get_indices_c_lib, get_indices_e_lib, get_indices_v_lib
  USE mo_fortran_tools, ONLY: init, set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: recon_lsq_cell_l_lib, recon_lsq_cell_l_svd_lib
  PUBLIC :: recon_lsq_cell_q_lib, recon_lsq_cell_q_svd_lib
  PUBLIC :: recon_lsq_cell_c_lib, recon_lsq_cell_c_svd_lib
  PUBLIC :: div_lib, div_avg_lib
  PUBLIC :: rot_vertex_lib, rot_vertex_ri_lib
  PUBLIC :: rot_vertex_atmos_lib

#ifdef __MIXED_PRECISION
  INTEGER, PARAMETER :: vp = sp
#else
  INTEGER, PARAMETER :: vp = wp
#endif

  INTERFACE rot_vertex_lib

    MODULE PROCEDURE rot_vertex_atmos_lib

  END INTERFACE
!
!
  INTERFACE div_lib

    MODULE PROCEDURE div3d_lib
    MODULE PROCEDURE div3d_2field_lib
    MODULE PROCEDURE div4d_lib
!
  END INTERFACE

CONTAINS

!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered linear
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! l    : linear reconstruction
!!
!! The least squares approach is used. Solves Rx = Q^T d.
!! R: upper triangular matrix (2 x 2)
!! Q: orthogonal matrix (3 x 2)
!! d: input vector (3 x 1)
!! x: solution vector (2 x 1)
!! works only on triangular grid yet
!!
  SUBROUTINE recon_lsq_cell_l_lib(p_cc, cell_neighbor_idx, cell_neighbor_blk, &
    &                             lsq_qtmat_c, lsq_rmat_rdiag_c, &
    &                             lsq_rmat_utri_c, lsq_moments, p_coeff, &
    &                             i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                             nlev, slev, elev, nproma, l_consv, lacc, acc_async)

    ! edge based cell centered variable of which divergence is computed
    REAL(wp), INTENT(IN) :: p_cc(:, :, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN)  :: cell_neighbor_idx(:, :, :)

    ! block indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN)  :: cell_neighbor_blk(:, :, :)

    ! transposed Q of QR-factorization for lsq reconstruction, dim: (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
    REAL(wp), INTENT(IN) :: lsq_qtmat_c(:, :, :, :)

    ! reciprocal diagonal elements of R-matrix resulting from QR-decomposition, dim: (nproma,lsq_dim_unk,nblks_c)
    REAL(wp), INTENT(IN) :: lsq_rmat_rdiag_c(:, :, :)

    ! upper triangular elements without diagonal elements of R-matrix (starting from the bottom right),
    ! dim: (nproma,(lsq_dim_unk^2-lsq_dim_unk)/2,nblks_c)
    REAL(wp), INTENT(IN) :: lsq_rmat_utri_c(:, :, :)

    ! Moments (x^ny^m)_{i} for control volume, dim: (nproma,nblks_c,lsq_dim_unk)
    REAL(wp), INTENT(IN) :: lsq_moments(:, :, :)

    ! cell based coefficients (geographical components)
    ! (constant and gradients in latitudinal and longitudinal direction)
    REAL(wp), INTENT(INOUT) :: p_coeff(:, :, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! total vertical level
    INTEGER, INTENT(IN) :: nlev

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, conservative reconstruction is used
    LOGICAL, INTENT(IN) :: l_consv

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! async OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: acc_async

#ifdef __LOOP_EXCHANGE
    REAL(wp)  ::   & !< weights * difference of scalars i j
      &  z_d(3, nproma, nlev)
#else
    REAL(wp)  :: z_d(3)
#endif
    REAL(wp)  ::   & !< matrix product of transposed Q matrix and d
      &  z_qt_times_d(2)

    INTEGER, POINTER ::   & !< Pointer to line and block indices of
      &  iidx(:, :, :), iblk(:, :, :) !< required stencil
    INTEGER :: jc, jk, jb !< index of cell, vertical level and block
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    ! pointers to line and block indices of stencil
    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk

    !
    ! 1. reconstruction of cell based gradient (geographical components)
    !
    !$ACC DATA PRESENT(p_coeff, p_cc, iidx, iblk, lsq_qtmat_c, lsq_rmat_rdiag_c, lsq_rmat_utri_c, lsq_moments) IF(lzacc)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

#ifdef __LOOP_EXCHANGE
      !$ACC DATA CREATE(z_d) IF(lzacc)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev

          ! note that the multiplication with lsq_weights_c(jc,js,jb) at
          ! runtime is now avoided. Instead, the multiplication with
          ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
          ! Q-matrix.
          z_d(1, jc, jk) = p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) - p_cc(jc, jk, jb)
          z_d(2, jc, jk) = p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) - p_cc(jc, jk, jb)
          z_d(3, jc, jk) = p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) - p_cc(jc, jk, jb)

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_qt_times_d)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! matrix multiplication Q^T d (partitioned into 2 dot products)
          z_qt_times_d(1) = lsq_qtmat_c(jc, 1, 1, jb)*z_d(1, jc, jk)  &
            &             + lsq_qtmat_c(jc, 1, 2, jb)*z_d(2, jc, jk)  &
            &             + lsq_qtmat_c(jc, 1, 3, jb)*z_d(3, jc, jk)
          z_qt_times_d(2) = lsq_qtmat_c(jc, 2, 1, jb)*z_d(1, jc, jk)  &
            &             + lsq_qtmat_c(jc, 2, 2, jb)*z_d(2, jc, jk)  &
            &             + lsq_qtmat_c(jc, 2, 3, jb)*z_d(3, jc, jk)

          ! Solve linear system by backward substitution
          ! Gradient in zonal and meridional direction
          !
          ! meridional
          p_coeff(3, jc, jk, jb) = lsq_rmat_rdiag_c(jc, 2, jb)*z_qt_times_d(2)

          ! zonal
          p_coeff(2, jc, jk, jb) = lsq_rmat_rdiag_c(jc, 1, jb)                  &
                                  & *(z_qt_times_d(1) - lsq_rmat_utri_c(jc, 1, jb)                &
                                     & *p_coeff(3, jc, jk, jb))

          ! constant
          p_coeff(1, jc, jk, jb) = p_cc(jc, jk, jb)

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL
      !$ACC WAIT
      !$ACC END DATA

#else
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_d, z_qt_times_d)
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          ! note that the multiplication with lsq_weights_c(jc,js,jb) at
          ! runtime is now avoided. Instead, the multiplication with
          ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
          ! Q-matrix.
          z_d(1) = p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) - p_cc(jc, jk, jb)
          z_d(2) = p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) - p_cc(jc, jk, jb)
          z_d(3) = p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) - p_cc(jc, jk, jb)

          ! matrix multiplication Q^T d (partitioned into 2 dot products)
          z_qt_times_d(1) = lsq_qtmat_c(jc, 1, 1, jb)*z_d(1)  &
            &             + lsq_qtmat_c(jc, 1, 2, jb)*z_d(2)  &
            &             + lsq_qtmat_c(jc, 1, 3, jb)*z_d(3)
          z_qt_times_d(2) = lsq_qtmat_c(jc, 2, 1, jb)*z_d(1)  &
            &             + lsq_qtmat_c(jc, 2, 2, jb)*z_d(2)  &
            &             + lsq_qtmat_c(jc, 2, 3, jb)*z_d(3)

          ! Solve linear system by backward substitution
          ! Gradient in zonal and meridional direction
          !
          ! meridional
          p_coeff(3, jc, jk, jb) = lsq_rmat_rdiag_c(jc, 2, jb)*z_qt_times_d(2)

          ! zonal
          p_coeff(2, jc, jk, jb) = lsq_rmat_rdiag_c(jc, 1, jb)                  &
                                  & *(z_qt_times_d(1) - lsq_rmat_utri_c(jc, 1, jb)                &
                                     & *p_coeff(3, jc, jk, jb))

          ! constant
          p_coeff(1, jc, jk, jb) = p_cc(jc, jk, jb)

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

#endif

      IF (l_consv) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
            ! constant
            p_coeff(1, jc, jk, jb) = p_coeff(1, jc, jk, jb)                                    &
              &                 - p_coeff(2, jc, jk, jb)*lsq_moments(jc, jb, 1) &
              &                 - p_coeff(3, jc, jk, jb)*lsq_moments(jc, jb, 2)

          END DO ! end loop over cells
        END DO ! end loop over vertical levels
        !$ACC END PARALLEL

      END IF

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

  END SUBROUTINE recon_lsq_cell_l_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered linear
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! l    : linear reconstruction
!!
!! The least squares approach is used. Solves Ax = b via Singular
!! Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape:
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (2 x 3)
!! b: input vector (3 x 1)
!! x: solution vector (2 x 1)
!! only works on triangular grid yet
!!
  SUBROUTINE recon_lsq_cell_l_svd_lib(p_cc, cell_neighbor_idx, cell_neighbor_blk, &
    &                                 lsq_pseudoinv, lsq_moments, p_coeff, &
    &                                 i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                 nlev, slev, elev, nproma, l_consv, lacc, acc_async)

    ! edge based cell centered variable of which divergence is computed
    REAL(wp), INTENT(IN)         :: p_cc(:, :, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN)  :: cell_neighbor_idx(:, :, :)

    ! block indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN)  :: cell_neighbor_blk(:, :, :)

    ! pseudo (or Moore-Penrose) inverse of lsq design matrix A, 
    ! dim: (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
    REAL(wp), INTENT(IN) :: lsq_pseudoinv(:, :, :, :)

    ! Moments (x^ny^m)_{i} for control volume, dim: (nproma,nblks_c,lsq_dim_unk)
    REAL(wp), INTENT(IN) :: lsq_moments(:, :, :)

    ! cell based coefficients (geographical components)
    ! (constant and gradients in latitudinal and longitudinal direction)
    REAL(wp), INTENT(INOUT) :: p_coeff(:, :, :, :)

    ! start_block needed for get_indices_c_lib 
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! total vertical level
    INTEGER, INTENT(IN) :: nlev 

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, conservative reconstruction is used
    LOGICAL, INTENT(IN) :: l_consv

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! async OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: acc_async

#ifdef __LOOP_EXCHANGE
    REAL(wp)  ::   & !< weights * difference of scalars i j
      &  z_b(3, nproma, nlev)
#else
    REAL(wp)  ::  z_b(3)
#endif

    INTEGER, POINTER ::   & !< Pointer to line and block indices of
      &  iidx(:, :, :), iblk(:, :, :) !< required stencil

    INTEGER :: jc, jk, jb !< index of cell, vertical level and block
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! pointers to line and block indices of stencil
    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

#ifdef __LOOP_EXCHANGE
      !$ACC DATA CREATE(z_b) IF(lzacc)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev

          ! note that the multiplication with lsq_weights_c(jc,js,jb) at
          ! runtime is now avoided. Instead, the weights have been shifted
          ! into the pseudoinverse.
          z_b(1, jc, jk) = p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) - p_cc(jc, jk, jb)
          z_b(2, jc, jk) = p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) - p_cc(jc, jk, jb)
          z_b(3, jc, jk) = p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) - p_cc(jc, jk, jb)

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

      !
      ! 2. compute cell based coefficients for linear reconstruction
      !    calculate matrix vector product PINV(A) * b
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! meridional
          p_coeff(3, jc, jk, jb) = lsq_pseudoinv(jc, 2, 1, jb)*z_b(1, jc, jk)  &
            &                 + lsq_pseudoinv(jc, 2, 2, jb)*z_b(2, jc, jk)  &
            &                 + lsq_pseudoinv(jc, 2, 3, jb)*z_b(3, jc, jk)

          ! zonal
          p_coeff(2, jc, jk, jb) = lsq_pseudoinv(jc, 1, 1, jb)*z_b(1, jc, jk)  &
            &                 + lsq_pseudoinv(jc, 1, 2, jb)*z_b(2, jc, jk)  &
            &                 + lsq_pseudoinv(jc, 1, 3, jb)*z_b(3, jc, jk)

          ! constant
          p_coeff(1, jc, jk, jb) = p_cc(jc, jk, jb)

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL
      !$ACC WAIT
      !$ACC END DATA

#else
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_b)
!$NEC outerloop_unroll(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! note that the multiplication with lsq_weights_c(jc,js,jb) at
          ! runtime is now avoided. Instead, the weights have been shifted
          ! into the pseudoinverse.
          z_b(1) = p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) - p_cc(jc, jk, jb)
          z_b(2) = p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) - p_cc(jc, jk, jb)
          z_b(3) = p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) - p_cc(jc, jk, jb)

          ! meridional
          p_coeff(3, jc, jk, jb) = lsq_pseudoinv(jc, 2, 1, jb)*z_b(1)  &
            &                 + lsq_pseudoinv(jc, 2, 2, jb)*z_b(2)  &
            &                 + lsq_pseudoinv(jc, 2, 3, jb)*z_b(3)

          ! zonal
          p_coeff(2, jc, jk, jb) = lsq_pseudoinv(jc, 1, 1, jb)*z_b(1)  &
            &                 + lsq_pseudoinv(jc, 1, 2, jb)*z_b(2)  &
            &                 + lsq_pseudoinv(jc, 1, 3, jb)*z_b(3)

          ! constant
          p_coeff(1, jc, jk, jb) = p_cc(jc, jk, jb)

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL
#endif

      IF (l_consv) THEN

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx

            ! In the case of a conservative reconstruction,
            ! the coefficient c0 is derived from the linear constraint
            !
            p_coeff(1, jc, jk, jb) = p_coeff(1, jc, jk, jb)                                    &
              &                 - p_coeff(2, jc, jk, jb)*lsq_moments(jc, jb, 1) &
              &                 - p_coeff(3, jc, jk, jb)*lsq_moments(jc, jb, 2)

          END DO ! end loop over cells
        END DO ! end loop over vertical levels
        !$ACC END PARALLEL

      END IF

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

  END SUBROUTINE recon_lsq_cell_l_svd_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered quadratic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! q    : quadratic reconstruction
!!
!! Computes the coefficients (derivatives) for a quadratic reconstruction,
!! using the the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!! Solves linear system Rx = Q^T d.
!! The matrices have the following size and shape:
!! R  : upper triangular matrix (5 x 5)
!! Q  : orthogonal matrix (9 x 5)
!! Q^T: transposed of Q (5 x 9)
!! d  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (5 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (0.5*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (0.5*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!!
!! works only on triangular grid yet
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
  SUBROUTINE recon_lsq_cell_q_lib(p_cc, lsq_idx_c, lsq_blk_c, &
    &                             lsq_rmat_rdiag_c, lsq_rmat_utri_c, &
    &                             lsq_moments, lsq_qtmat_c, p_coeff, &
    &                             i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                             nlev, slev, elev, nproma, patch_id, lsq_high_set_dim_c, l_limited_area, lacc)

    ! edge based cell centered variable of which divergence is computed
    REAL(wp), INTENT(IN) ::  p_cc(:, :, :)

    ! index array defining the stencil for lsq reconstruction (nproma,nblks_c,lsq_dim_c)
    INTEGER, TARGET, INTENT(IN)  :: lsq_idx_c(:, :, :)

    ! block index array defining the stencil for lsq reconstruction (nproma,nblks_c,lsq_dim_c)
    INTEGER, TARGET, INTENT(IN)  :: lsq_blk_c(:, :, :)

    ! reciprocal diagonal elements of R-matrix resulting from QR-decomposition (nproma,lsq_dim_unk,nblks_c)
    REAL(wp), TARGET, INTENT(IN) :: lsq_rmat_rdiag_c(:, :, :)

    ! upper triangular elements without diagonal elements of R-matrix (starting from the bottom right)
    ! (nproma,(lsq_dim_unk^2-lsq_dim_unk)/2,nblks_c)
    REAL(wp), TARGET, INTENT(IN) :: lsq_rmat_utri_c(:, :, :)

    ! Moments (x^ny^m)_{i} for control volume (nproma,nblks_c,lsq_dim_unk)
    REAL(wp), INTENT(IN) :: lsq_moments(:, :, :)

    ! transposed Q of QR-factorization for lsq reconstruction (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
    REAL(wp), INTENT(IN) :: lsq_qtmat_c(:, :, :, :)

    ! cell based coefficients (geographical components) phsically this vector contains gradients, second
    ! derivatives, one mixed derivative and a constant coefficient for zonal and meridional direction
    REAL(wp), INTENT(INOUT) :: p_coeff(:, :, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! total vertical level
    INTEGER, INTENT(IN) :: nlev

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! domain ID of current domain
    INTEGER, INTENT(IN)  :: patch_id

    ! parameter determining the size of the lsq stencil
    INTEGER, INTENT(IN)  :: lsq_high_set_dim_c

    ! limited area setup where forcing comes in from sides
    LOGICAL, INTENT(IN)  :: l_limited_area

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp)  ::           & !< difference of scalars i j
      &  z_d(lsq_high_set_dim_c, nproma, nlev)
    REAL(wp)  ::           & !< matrix-vector product of transposed
      &  z_qt_times_d(5) !< Q matrix and d

    REAL(wp), POINTER ::   & !< Pointer to reciprocal diagonal R-matrix-elements
      &  ptr_rrdiag(:, :, :)
    REAL(wp), POINTER ::   & !< Pointer to upper triangular R-matrix-elements
      &  ptr_rutri(:, :, :)

    INTEGER, POINTER  ::   & !< Pointer to line and block indices of
      &  iidx(:, :, :), iblk(:, :, :) !< required stencil
    INTEGER :: jc, jk, jb !< index of cell, vertical level and block
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! pointers to line and block indices of required stencil
    iidx => lsq_idx_c
    iblk => lsq_blk_c

    ! pointer to reciprocal diagonal R-elements
    ptr_rrdiag => lsq_rmat_rdiag_c(:, :, :)

    ! pointer to upper triangular R-elements
    ptr_rutri => lsq_rmat_utri_c(:, :, :)

    !$ACC DATA PRESENT(p_cc, p_coeff, lsq_moments, lsq_qtmat_c, iidx, iblk, ptr_rrdiag, ptr_rutri) &
    !$ACC   CREATE(z_d) IF(lzacc)
!$OMP PARALLEL

    IF (patch_id > 1 .OR. l_limited_area) THEN
      CALL init(p_coeff(:, :, 1:6, 1:i_startblk), lacc=lzacc)
!$OMP BARRIER
    END IF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !
      ! 1. compute right hand side of linear system
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          z_d(1, jc, jk) = p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) - p_cc(jc, jk, jb)
          z_d(2, jc, jk) = p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) - p_cc(jc, jk, jb)
          z_d(3, jc, jk) = p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) - p_cc(jc, jk, jb)
          z_d(4, jc, jk) = p_cc(iidx(jc, jb, 4), jk, iblk(jc, jb, 4)) - p_cc(jc, jk, jb)
          z_d(5, jc, jk) = p_cc(iidx(jc, jb, 5), jk, iblk(jc, jb, 5)) - p_cc(jc, jk, jb)
          z_d(6, jc, jk) = p_cc(iidx(jc, jb, 6), jk, iblk(jc, jb, 6)) - p_cc(jc, jk, jb)
          z_d(7, jc, jk) = p_cc(iidx(jc, jb, 7), jk, iblk(jc, jb, 7)) - p_cc(jc, jk, jb)
          z_d(8, jc, jk) = p_cc(iidx(jc, jb, 8), jk, iblk(jc, jb, 8)) - p_cc(jc, jk, jb)
          z_d(9, jc, jk) = p_cc(iidx(jc, jb, 9), jk, iblk(jc, jb, 9)) - p_cc(jc, jk, jb)

        END DO
      END DO
      !$ACC END PARALLEL

      !
      ! 2. compute cell based coefficients for quadratic reconstruction
      !

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG
      DO jk = slev, elev

        !$ACC LOOP VECTOR PRIVATE(z_qt_times_d)
        DO jc = i_startidx, i_endidx

          ! calculate matrix vector product Q^T d (transposed of Q times LHS)
          ! (intrinsic function matmul not applied, due to massive
          ! performance penalty on the NEC. Instead the intrinsic dot product
          ! function is applied
          z_qt_times_d(1) = DOT_PRODUCT(lsq_qtmat_c(jc, 1, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(2) = DOT_PRODUCT(lsq_qtmat_c(jc, 2, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(3) = DOT_PRODUCT(lsq_qtmat_c(jc, 3, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(4) = DOT_PRODUCT(lsq_qtmat_c(jc, 4, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(5) = DOT_PRODUCT(lsq_qtmat_c(jc, 5, 1:9, jb), z_d(1:9, jc, jk))

          !
          ! Solve linear system Rx=Q^T d by back substitution
          !
          p_coeff(6, jc, jk, jb) = ptr_rrdiag(jc, 5, jb)*z_qt_times_d(5)
          p_coeff(5, jc, jk, jb) = ptr_rrdiag(jc, 4, jb)                                       &
            &                 *(z_qt_times_d(4) - ptr_rutri(jc, 1, jb)*p_coeff(6, jc, jk, jb))
          p_coeff(4, jc, jk, jb) = ptr_rrdiag(jc, 3, jb)                                       &
            &                 *(z_qt_times_d(3) - ptr_rutri(jc, 2, jb)*p_coeff(5, jc, jk, jb)  &
            &                 - ptr_rutri(jc, 3, jb)*p_coeff(6, jc, jk, jb))
          p_coeff(3, jc, jk, jb) = ptr_rrdiag(jc, 2, jb)                                       &
            &                 *(z_qt_times_d(2) - ptr_rutri(jc, 4, jb)*p_coeff(4, jc, jk, jb)  &
            &                 - ptr_rutri(jc, 5, jb)*p_coeff(5, jc, jk, jb)                    &
            &                 - ptr_rutri(jc, 6, jb)*p_coeff(6, jc, jk, jb))
          p_coeff(2, jc, jk, jb) = ptr_rrdiag(jc, 1, jb)                                       &
            &                 *(z_qt_times_d(1) - ptr_rutri(jc, 7, jb)*p_coeff(3, jc, jk, jb)  &
            &                 - ptr_rutri(jc, 8, jb)*p_coeff(4, jc, jk, jb)                    &
            &                 - ptr_rutri(jc, 9, jb)*p_coeff(5, jc, jk, jb)                    &
            &                 - ptr_rutri(jc, 10, jb)*p_coeff(6, jc, jk, jb))

          p_coeff(1, jc, jk, jb) = p_cc(jc, jk, jb)                                            &
            &                 - p_coeff(2, jc, jk, jb)*lsq_moments(jc, jb, 1)                  &
            &                 - p_coeff(3, jc, jk, jb)*lsq_moments(jc, jb, 2)                  &
            &                 - p_coeff(4, jc, jk, jb)*lsq_moments(jc, jb, 3)                  &
            &                 - p_coeff(5, jc, jk, jb)*lsq_moments(jc, jb, 4)                  &
            &                 - p_coeff(6, jc, jk, jb)*lsq_moments(jc, jb, 5)   

        END DO ! end loop over cells

      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE recon_lsq_cell_q_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered quadratic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! q    : quadratic reconstruction
!!
!! Computes unknown coefficients (derivatives) of a quadratic polynomial,
!! using the least-squares method. The coefficients are provided at cell
!! centers in a local 2D cartesian system (tangential plane).
!!
!! Mathematically we solve Ax = b via Singular Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape (triangular grid) :
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (5 x 9)
!! b  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (5 x 1)
!!
!! Coefficients:
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (0.5*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (0.5*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
  SUBROUTINE recon_lsq_cell_q_svd_lib(p_cc, lsq_idx_c, lsq_blk_c, &
      &                                lsq_moments, lsq_pseudoinv, p_coeff, &
      &                                i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
      &                                nlev, slev, elev, nproma, patch_id, lsq_high_set_dim_c, l_limited_area, lacc)

    ! edge based cell centered variable of which divergence is computed
    REAL(wp), INTENT(IN) ::  p_cc(:, :, :)

    ! index array defining the stencil for lsq reconstruction (nproma,nblks_c,lsq_dim_c)
    INTEGER, TARGET, INTENT(IN)  :: lsq_idx_c(:, :, :)

    ! block index array defining the stencil for lsq reconstruction (nproma,nblks_c,lsq_dim_c)
    INTEGER, TARGET, INTENT(IN)  :: lsq_blk_c(:, :, :)

    ! Moments (x^ny^m)_{i} for control volume (nproma,nblks_c,lsq_dim_unk)
    REAL(wp), INTENT(IN) :: lsq_moments(:, :, :)

    ! pseudo (or Moore-Penrose) inverse of lsq design matrix A, (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
    REAL(wp), INTENT(IN) :: lsq_pseudoinv(:, :, :, :)

    ! cell based coefficients (geographical components) phsically this vector contains gradients, second
    ! derivatives, one mixed derivative and a constant coefficient for zonal and meridional direction
    REAL(wp), INTENT(INOUT) :: p_coeff(:, :, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! total vertical level
    INTEGER, INTENT(IN) :: nlev

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! domain ID of current domain
    INTEGER, INTENT(IN) :: patch_id

    ! parameter determining the size of the lsq stencil
    INTEGER, INTENT(IN) :: lsq_high_set_dim_c

    ! limited area setup where forcing comes in from sides
    LOGICAL, INTENT(IN) :: l_limited_area

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp)  ::           & !< difference of scalars i j
      &  z_b(lsq_high_set_dim_c, nproma, nlev)

    INTEGER, POINTER  ::   & !< Pointer to line and block indices of
      &  iidx(:, :, :), iblk(:, :, :) !< required stencil
    INTEGER :: jc, jk, jb !< index of cell, vertical level and block
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! pointers to line and block indices of required stencil
    iidx => lsq_idx_c
    iblk => lsq_blk_c

    !$ACC DATA PRESENT(p_cc, p_coeff, lsq_moments, lsq_pseudoinv, iidx, iblk) &
    !$ACC   CREATE(z_b) IF(lzacc)
!$OMP PARALLEL

    IF (patch_id > 1 .OR. l_limited_area) THEN
      CALL init(p_coeff(:, :, 1:6, 1:i_startblk), lacc=lzacc)
!$OMP BARRIER
    END IF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !
      ! 1. compute right hand side of linear system
      !

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          z_b(1, jc, jk) = p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) - p_cc(jc, jk, jb)
          z_b(2, jc, jk) = p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) - p_cc(jc, jk, jb)
          z_b(3, jc, jk) = p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) - p_cc(jc, jk, jb)
          z_b(4, jc, jk) = p_cc(iidx(jc, jb, 4), jk, iblk(jc, jb, 4)) - p_cc(jc, jk, jb)
          z_b(5, jc, jk) = p_cc(iidx(jc, jb, 5), jk, iblk(jc, jb, 5)) - p_cc(jc, jk, jb)
          z_b(6, jc, jk) = p_cc(iidx(jc, jb, 6), jk, iblk(jc, jb, 6)) - p_cc(jc, jk, jb)
          z_b(7, jc, jk) = p_cc(iidx(jc, jb, 7), jk, iblk(jc, jb, 7)) - p_cc(jc, jk, jb)
          z_b(8, jc, jk) = p_cc(iidx(jc, jb, 8), jk, iblk(jc, jb, 8)) - p_cc(jc, jk, jb)
          z_b(9, jc, jk) = p_cc(iidx(jc, jb, 9), jk, iblk(jc, jb, 9)) - p_cc(jc, jk, jb)

        END DO
      END DO
      !$ACC END PARALLEL

      !
      ! 2. compute cell based coefficients for quadratic reconstruction
      !    calculate matrix vector product PINV(A) * b
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG
      DO jk = slev, elev

        !$ACC LOOP VECTOR
!$NEC ivdep
        DO jc = i_startidx, i_endidx

          ! (intrinsic function matmul not applied, due to massive
          ! performance penalty on the NEC. Instead the intrinsic dot product
          ! function is used.
          p_coeff(6, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 5, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(5, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 4, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(4, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 3, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(3, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 2, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(2, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 1, 1:9, jb), &
            &                               z_b(1:9, jc, jk))

          ! At the end, the coefficient c0 is derived from the linear constraint
          !
          p_coeff(1, jc, jk, jb) = p_cc(jc, jk, jb) - DOT_PRODUCT(p_coeff(2:6, jc, jk, jb), &
            &                   lsq_moments(jc, jb, 1:5))

        END DO ! end loop over cells

      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE recon_lsq_cell_q_svd_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered cubic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes the coefficients (derivatives) for a cubic reconstruction,
!! using the the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!! Solves linear system Rx = Q^T d.
!! The matrices have the following size and shape:
!! R  : upper triangular matrix (9 x 9)
!! Q  : orthogonal matrix (9 x 9)
!! Q^T: transposed of Q (9 x 9)
!! d  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (9 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb, 1) : C0
!! p_coeff(jc,jk,jb, 2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb, 3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb, 4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb, 5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb, 6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb, 7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb, 8) : C7 (1/6*dddPhi_dddy)
!! p_coeff(jc,jk,jb, 9) : C8 (1/2*dddPhi_ddxdy)
!! p_coeff(jc,jk,jb,10) : C9 (1/2*dddPhi_dxddy)
!!
!! works only on triangular grid yet
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
  SUBROUTINE recon_lsq_cell_c_lib(p_cc, lsq_idx_c, lsq_blk_c, &
      &                        lsq_rmat_rdiag_c, lsq_rmat_utri_c, &
      &                        lsq_moments, lsq_qtmat_c, p_coeff, &
      &                        i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
      &                        nlev, slev, elev, nproma, patch_id, lsq_high_set_dim_c, l_limited_area, lacc)

    ! edge based cell centered variable of which divergence is computed
    REAL(wp), INTENT(IN) ::  p_cc(:, :, :)

    ! index array defining the stencil for lsq reconstruction, dim: (nproma,nblks_c,lsq_dim_c)
    INTEGER, TARGET, INTENT(IN)  :: lsq_idx_c(:, :, :)

    ! block array defining the stencil for lsq reconstruction, dim: (nproma,nblks_c,lsq_dim_c)
    INTEGER, TARGET, INTENT(IN)  :: lsq_blk_c(:, :, :)

    ! reciprocal diagonal elements of R-matrix resulting from QR-decomposition, 
    ! dim: (nproma,lsq_dim_unk,nblks_c)
    REAL(wp), TARGET, INTENT(IN) :: lsq_rmat_rdiag_c(:, :, :)

    ! upper triangular elements without diagonal elements of R-matrix (starting from the bottom right),
    ! dim: (nproma,(lsq_dim_unk^2-lsq_dim_unk)/2,nblks_c)
    REAL(wp), TARGET, INTENT(IN) :: lsq_rmat_utri_c(:, :, :)

    ! Moments (x^ny^m)_{i} for control volume, dim: (nproma,nblks_c,lsq_dim_unk)
    REAL(wp), INTENT(IN) :: lsq_moments(:, :, :)

    ! transposed Q of QR-factorization for lsq reconstruction,
    ! dim: (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
    REAL(wp), INTENT(IN) :: lsq_qtmat_c(:, :, :, :)

    ! cell based coefficients (geographical components) phsically this vector contains gradients, second
    ! derivatives, one mixed derivative and a constant coefficient for zonal and meridional direction
    REAL(wp), INTENT(INOUT) :: p_coeff(:, :, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! total vertical level
    INTEGER, INTENT(IN) :: nlev

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! domain ID of current domain
    INTEGER, INTENT(IN) :: patch_id

    ! parameter determining the size of the lsq stencil
    INTEGER, INTENT(IN) :: lsq_high_set_dim_c

    ! limited area setup where forcing comes in from sides
    LOGICAL, INTENT(IN) :: l_limited_area

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp)  ::           & !< difference of scalars i j
      &  z_d(lsq_high_set_dim_c, nproma, nlev)
    REAL(wp)  ::           & !< matrix-vector product of transposed
      &  z_qt_times_d(9) !< Q matrix and d

    REAL(wp), POINTER ::   & !< Pointer to reciprocal diagonal R-matrix-elements
      &  ptr_rrdiag(:, :, :)
    REAL(wp), POINTER ::   & !< Pointer to upper triangular R-matrix-elements
      &  ptr_rutri(:, :, :)

    INTEGER, POINTER  ::   & !< Pointer to line and block indices of
      &  iidx(:, :, :), iblk(:, :, :) !< required stencil
    INTEGER :: jc, jk, jb !< index of cell, vertical level and block
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! pointers to line and block indices of required stencil
    iidx => lsq_idx_c
    iblk => lsq_blk_c

    ! pointer to reciprocal diagonal R-elements
    ptr_rrdiag => lsq_rmat_rdiag_c(:, :, :)

    ! pointer to upper triangular R-elements
    ptr_rutri => lsq_rmat_utri_c(:, :, :)

    !$ACC DATA PRESENT(p_cc, p_coeff, lsq_moments, lsq_qtmat_c, iidx, iblk) &
    !$ACC   CREATE(z_d) IF(lzacc)
!$OMP PARALLEL

    IF (patch_id > 1 .OR. l_limited_area) THEN
      CALL init(p_coeff(:, :, 1:10, 1:i_startblk), lacc=lzacc)
!$OMP BARRIER
    END IF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !
      ! 1. compute right hand side of linear system
      !
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        !$ACC LOOP VECTOR
!NEC$ ivdep
        DO jc = i_startidx, i_endidx
#endif

          z_d(1, jc, jk) = p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) - p_cc(jc, jk, jb)
          z_d(2, jc, jk) = p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) - p_cc(jc, jk, jb)
          z_d(3, jc, jk) = p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) - p_cc(jc, jk, jb)
          z_d(4, jc, jk) = p_cc(iidx(jc, jb, 4), jk, iblk(jc, jb, 4)) - p_cc(jc, jk, jb)
          z_d(5, jc, jk) = p_cc(iidx(jc, jb, 5), jk, iblk(jc, jb, 5)) - p_cc(jc, jk, jb)
          z_d(6, jc, jk) = p_cc(iidx(jc, jb, 6), jk, iblk(jc, jb, 6)) - p_cc(jc, jk, jb)
          z_d(7, jc, jk) = p_cc(iidx(jc, jb, 7), jk, iblk(jc, jb, 7)) - p_cc(jc, jk, jb)
          z_d(8, jc, jk) = p_cc(iidx(jc, jb, 8), jk, iblk(jc, jb, 8)) - p_cc(jc, jk, jb)
          z_d(9, jc, jk) = p_cc(iidx(jc, jb, 9), jk, iblk(jc, jb, 9)) - p_cc(jc, jk, jb)

        END DO
      END DO
      !$ACC END PARALLEL

      !
      ! 2. compute cell based coefficients for quadratic reconstruction
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG
      DO jk = slev, elev

        !$ACC LOOP VECTOR PRIVATE(z_qt_times_d)
!NEC$ ivdep
        DO jc = i_startidx, i_endidx

          ! calculate matrix vector product Q^T d (transposed of Q times LHS)
          ! (intrinsic function matmul not applied, due to massive
          ! performance penalty on the NEC. Instead the intrinsic dot product
          ! function is applied
!TODO:  these should be nine scalars, since they should reside in registers
          z_qt_times_d(1) = DOT_PRODUCT(lsq_qtmat_c(jc, 1, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(2) = DOT_PRODUCT(lsq_qtmat_c(jc, 2, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(3) = DOT_PRODUCT(lsq_qtmat_c(jc, 3, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(4) = DOT_PRODUCT(lsq_qtmat_c(jc, 4, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(5) = DOT_PRODUCT(lsq_qtmat_c(jc, 5, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(6) = DOT_PRODUCT(lsq_qtmat_c(jc, 6, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(7) = DOT_PRODUCT(lsq_qtmat_c(jc, 7, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(8) = DOT_PRODUCT(lsq_qtmat_c(jc, 8, 1:9, jb), z_d(1:9, jc, jk))
          z_qt_times_d(9) = DOT_PRODUCT(lsq_qtmat_c(jc, 9, 1:9, jb), z_d(1:9, jc, jk))

          !
          ! Solve linear system Rx=Q^T d by back substitution
          !
          p_coeff(10, jc, jk, jb) = ptr_rrdiag(jc, 9, jb)*z_qt_times_d(9)
          p_coeff(9, jc, jk, jb) = ptr_rrdiag(jc, 8, jb)                                         &
            &                  *(z_qt_times_d(8) - ptr_rutri(jc, 1, jb)*p_coeff(10, jc, jk, jb))
          p_coeff(8, jc, jk, jb) = ptr_rrdiag(jc, 7, jb)                                         &
            &                  *(z_qt_times_d(7) - (ptr_rutri(jc, 2, jb)*p_coeff(9, jc, jk, jb) &
            &                  + ptr_rutri(jc, 3, jb)*p_coeff(10, jc, jk, jb)))
          p_coeff(7, jc, jk, jb) = ptr_rrdiag(jc, 6, jb)                                         &
            &                  *(z_qt_times_d(6) - (ptr_rutri(jc, 4, jb)*p_coeff(8, jc, jk, jb) &
            &                  + ptr_rutri(jc, 5, jb)*p_coeff(9, jc, jk, jb)                    &
            &                  + ptr_rutri(jc, 6, jb)*p_coeff(10, jc, jk, jb)))
          p_coeff(6, jc, jk, jb) = ptr_rrdiag(jc, 5, jb)                                         &
            &                  *(z_qt_times_d(5) - (ptr_rutri(jc, 7, jb)*p_coeff(7, jc, jk, jb) &
            &                  + ptr_rutri(jc, 8, jb)*p_coeff(8, jc, jk, jb)                    &
            &                  + ptr_rutri(jc, 9, jb)*p_coeff(9, jc, jk, jb)                    &
            &                  + ptr_rutri(jc, 10, jb)*p_coeff(10, jc, jk, jb)))
          p_coeff(5, jc, jk, jb) = ptr_rrdiag(jc, 4, jb)                                         &
            &                  *(z_qt_times_d(4) - (ptr_rutri(jc, 11, jb)*p_coeff(6, jc, jk, jb)&
            &                  + ptr_rutri(jc, 12, jb)*p_coeff(7, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 13, jb)*p_coeff(8, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 14, jb)*p_coeff(9, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 15, jb)*p_coeff(10, jc, jk, jb)))
          p_coeff(4, jc, jk, jb) = ptr_rrdiag(jc, 3, jb)                                         &
            &                  *(z_qt_times_d(3) - (ptr_rutri(jc, 16, jb)*p_coeff(5, jc, jk, jb)&
            &                  + ptr_rutri(jc, 17, jb)*p_coeff(6, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 18, jb)*p_coeff(7, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 19, jb)*p_coeff(8, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 20, jb)*p_coeff(9, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 21, jb)*p_coeff(10, jc, jk, jb)))
          p_coeff(3, jc, jk, jb) = ptr_rrdiag(jc, 2, jb)                                         &
            &                  *(z_qt_times_d(2) - (ptr_rutri(jc, 22, jb)*p_coeff(4, jc, jk, jb)&
            &                  + ptr_rutri(jc, 23, jb)*p_coeff(5, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 24, jb)*p_coeff(6, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 25, jb)*p_coeff(7, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 26, jb)*p_coeff(8, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 27, jb)*p_coeff(9, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 28, jb)*p_coeff(10, jc, jk, jb)))
          p_coeff(2, jc, jk, jb) = ptr_rrdiag(jc, 1, jb)                                         &
            &                  *(z_qt_times_d(1) - (ptr_rutri(jc, 29, jb)*p_coeff(3, jc, jk, jb)&
            &                  + ptr_rutri(jc, 30, jb)*p_coeff(4, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 31, jb)*p_coeff(5, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 32, jb)*p_coeff(6, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 33, jb)*p_coeff(7, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 34, jb)*p_coeff(8, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 35, jb)*p_coeff(9, jc, jk, jb)                   &
            &                  + ptr_rutri(jc, 36, jb)*p_coeff(10, jc, jk, jb)))

          p_coeff(1, jc, jk, jb) = p_cc(jc, jk, jb) - (                                          &
            &                    p_coeff(2, jc, jk, jb)*lsq_moments(jc, jb, 1)     &
            &                  + p_coeff(3, jc, jk, jb)*lsq_moments(jc, jb, 2)     &
            &                  + p_coeff(4, jc, jk, jb)*lsq_moments(jc, jb, 3)     &
            &                  + p_coeff(5, jc, jk, jb)*lsq_moments(jc, jb, 4)     &
            &                  + p_coeff(6, jc, jk, jb)*lsq_moments(jc, jb, 5)     &
            &                  + p_coeff(7, jc, jk, jb)*lsq_moments(jc, jb, 6)     &
            &                  + p_coeff(8, jc, jk, jb)*lsq_moments(jc, jb, 7)     &
            &                  + p_coeff(9, jc, jk, jb)*lsq_moments(jc, jb, 8)     &
            &                  + p_coeff(10, jc, jk, jb)*lsq_moments(jc, jb, 9))

        END DO ! end loop over cells

      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE recon_lsq_cell_c_lib

!--------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered cubic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes unknown coefficients (derivatives) of a cubic polynomial,
!! using the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!!
!! Mathematically we solve Ax = b via Singular Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape (triangular grid) :
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (9 x 9)
!! b  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (9 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb, 1) : C0
!! p_coeff(jc,jk,jb, 2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb, 3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb, 4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb, 5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb, 6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb, 7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb, 8) : C7 (1/6*dddPhi_dddy)
!! p_coeff(jc,jk,jb, 9) : C8 (1/2*dddPhi_ddxdy)
!! p_coeff(jc,jk,jb,10) : C9 (1/2*dddPhi_dxddy)
!!
!! works only on triangular grid yet
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
  SUBROUTINE recon_lsq_cell_c_svd_lib(p_cc, lsq_idx_c, lsq_blk_c, &
    &                                 lsq_moments, lsq_pseudoinv, p_coeff, &
    &                                 i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                 i_startblk_init, i_endblk_init, i_startidx_init, i_endidx_init, &
    &                                 nlev, slev, elev, nproma, patch_id, lsq_high_set_dim_c, l_limited_area, lacc)

    ! edge based cell centered variable of which divergence is computed
    REAL(wp), INTENT(IN) ::  p_cc(:, :, :)

    ! index array defining the stencil for lsq reconstruction, dim: (nproma,nblks_c,lsq_dim_c)
    INTEGER, TARGET, INTENT(IN)  :: lsq_idx_c(:, :, :)

    ! block array defining the stencil for lsq reconstruction, dim: (nproma,nblks_c,lsq_dim_c)
    INTEGER, TARGET, INTENT(IN)  :: lsq_blk_c(:, :, :)

    ! Moments (x^ny^m)_{i} for control volume, dim: (nproma,nblks_c,lsq_dim_unk)
    REAL(wp), INTENT(IN) :: lsq_moments(:, :, :)

    ! pseudo (or Moore-Penrose) inverse of lsq design matrix A, dim: (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
    REAL(wp), INTENT(IN) :: lsq_pseudoinv(:, :, :, :)

    ! cell based coefficients (geographical components) phsically this vector contains gradients, second
    ! derivatives, one mixed derivative and a constant coefficient for zonal and meridional direction
    REAL(wp), INTENT(INOUT) :: p_coeff(:, :, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk_init

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk_init

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_init

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_init

    ! total vertical level
    INTEGER, INTENT(IN) :: nlev

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! domain ID of current domain
    INTEGER, INTENT(IN)  :: patch_id

    ! parameter determining the size of the lsq stencil
    INTEGER, INTENT(IN)  :: lsq_high_set_dim_c

    ! limited area setup where forcing comes in from sides
    LOGICAL, INTENT(IN)  :: l_limited_area

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

#ifdef __LOOP_EXCHANGE
    REAL(wp)  ::           & !< difference of scalars i j
      &  z_b(lsq_high_set_dim_c, nproma, nlev)
#else
    REAL(wp)  ::           & !< difference of scalars i j
      &  z_b(9)
#endif

    INTEGER, POINTER  ::   & !< Pointer to line and block indices of
      &  iidx(:, :, :), iblk(:, :, :) !< required stencil
    INTEGER :: jc, jk, jb !< index of cell, vertical level and block
    INTEGER :: i_startidx, i_endidx
    INTEGER :: ji
    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! pointers to line and block indices of required stencil
    iidx => lsq_idx_c
    iblk => lsq_blk_c

    !$ACC DATA PRESENT(p_cc, p_coeff, lsq_moments, lsq_pseudoinv, iidx, iblk) IF(lzacc)
!$OMP PARALLEL

    IF (patch_id > 1 .OR. l_limited_area) THEN
      !CALL init(p_coeff(:,:,:,1:i_startblk))
      ! Only zero-init the lateral boundary points

!$OMP DO PRIVATE(jb,jc,jk,ji,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk_init, i_endblk_init

        CALL get_indices_c_lib(i_startidx_init, i_endidx_init, nproma, jb, i_startblk_init, i_endblk_init, &
                               i_startidx, i_endidx)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(lzacc)
!NEC$ forced_collapse
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
            DO ji = 1, 10
              p_coeff(ji, jc, jk, jb) = 0._wp
            END DO
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP BARRIER
!$OMP BARRIER
    END IF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !
      ! 1. compute right hand side of linear system
      !

#ifdef __LOOP_EXCHANGE
      !$ACC DATA CREATE(z_b) IF(lzacc)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev

          z_b(1, jc, jk) = p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) - p_cc(jc, jk, jb)
          z_b(2, jc, jk) = p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) - p_cc(jc, jk, jb)
          z_b(3, jc, jk) = p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) - p_cc(jc, jk, jb)
          z_b(4, jc, jk) = p_cc(iidx(jc, jb, 4), jk, iblk(jc, jb, 4)) - p_cc(jc, jk, jb)
          z_b(5, jc, jk) = p_cc(iidx(jc, jb, 5), jk, iblk(jc, jb, 5)) - p_cc(jc, jk, jb)
          z_b(6, jc, jk) = p_cc(iidx(jc, jb, 6), jk, iblk(jc, jb, 6)) - p_cc(jc, jk, jb)
          z_b(7, jc, jk) = p_cc(iidx(jc, jb, 7), jk, iblk(jc, jb, 7)) - p_cc(jc, jk, jb)
          z_b(8, jc, jk) = p_cc(iidx(jc, jb, 8), jk, iblk(jc, jb, 8)) - p_cc(jc, jk, jb)
          z_b(9, jc, jk) = p_cc(iidx(jc, jb, 9), jk, iblk(jc, jb, 9)) - p_cc(jc, jk, jb)

        END DO
      END DO
      !$ACC END PARALLEL

      !
      ! 2. compute cell based coefficients for cubic reconstruction
      !    calculate matrix vector product PINV(A) * b
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! (intrinsic function matmul not applied, due to massive
          ! performance penalty on the NEC. Instead the intrinsic dot product
          ! function is applied

          p_coeff(10, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 9, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(9, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 8, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(8, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 7, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(7, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 6, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(6, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 5, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(5, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 4, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(4, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 3, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(3, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 2, 1:9, jb), &
            &                               z_b(1:9, jc, jk))
          p_coeff(2, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 1, 1:9, jb), &
            &                               z_b(1:9, jc, jk))

          p_coeff(1, jc, jk, jb) = p_cc(jc, jk, jb) - DOT_PRODUCT(p_coeff(2:10, jc, jk, jb), &
            &                    lsq_moments(jc, jb, 1:9))

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL
      !$ACC WAIT
      !$ACC END DATA

#else
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_b)
      DO jk = slev, elev
!$NEC ivdep
        DO jc = i_startidx, i_endidx

          !
          ! 1. compute right hand side of linear system
          !
          z_b(1) = p_cc(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) - p_cc(jc, jk, jb)
          z_b(2) = p_cc(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) - p_cc(jc, jk, jb)
          z_b(3) = p_cc(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) - p_cc(jc, jk, jb)
          z_b(4) = p_cc(iidx(jc, jb, 4), jk, iblk(jc, jb, 4)) - p_cc(jc, jk, jb)
          z_b(5) = p_cc(iidx(jc, jb, 5), jk, iblk(jc, jb, 5)) - p_cc(jc, jk, jb)
          z_b(6) = p_cc(iidx(jc, jb, 6), jk, iblk(jc, jb, 6)) - p_cc(jc, jk, jb)
          z_b(7) = p_cc(iidx(jc, jb, 7), jk, iblk(jc, jb, 7)) - p_cc(jc, jk, jb)
          z_b(8) = p_cc(iidx(jc, jb, 8), jk, iblk(jc, jb, 8)) - p_cc(jc, jk, jb)
          z_b(9) = p_cc(iidx(jc, jb, 9), jk, iblk(jc, jb, 9)) - p_cc(jc, jk, jb)

          !
          ! 2. compute cell based coefficients for cubic reconstruction
          !    calculate matrix vector product PINV(A) * b
          !

          ! (intrinsic function matmul not applied, due to massive
          ! performance penalty on the NEC. Instead the intrinsic dot product
          ! function is applied

          p_coeff(10, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 9, 1:9, jb), &
            &                               z_b(1:9))
          p_coeff(9, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 8, 1:9, jb), &
            &                               z_b(1:9))
          p_coeff(8, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 7, 1:9, jb), &
            &                               z_b(1:9))
          p_coeff(7, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 6, 1:9, jb), &
            &                               z_b(1:9))
          p_coeff(6, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 5, 1:9, jb), &
            &                               z_b(1:9))
          p_coeff(5, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 4, 1:9, jb), &
            &                               z_b(1:9))
          p_coeff(4, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 3, 1:9, jb), &
            &                               z_b(1:9))
          p_coeff(3, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 2, 1:9, jb), &
            &                               z_b(1:9))
          p_coeff(2, jc, jk, jb) = DOT_PRODUCT(lsq_pseudoinv(jc, 1, 1:9, jb), &
            &                               z_b(1:9))

          p_coeff(1, jc, jk, jb) = p_cc(jc, jk, jb) - DOT_PRODUCT(p_coeff(2:10, jc, jk, jb), &
            &                    lsq_moments(jc, jb, 1:9))

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

#endif

    END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE recon_lsq_cell_c_svd_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field.
!!
!! Computes discrete divergence of a vector field
!! given by its components in the directions normal to triangle edges.
!! The midpoint rule is used for quadrature.
!! input:  lives on edges (velocity points)
!! output: lives on centers of triangles
!!
  SUBROUTINE div3d_lib(vec_e, cell_edge_idx, cell_edge_blk, &
    &                  geofac_div, div_vec_c, &
    &                  i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                  slev, elev, nproma, lacc)

    ! edge based variable of which divergence is computed, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: vec_e(:, :, :)

    ! line indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_idx(:, :, :)

    ! block indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_blk(:, :, :)

    !Interpolation state
    ! factor for divergence, dim: (nproma,cell_type,nblks_c)
    REAL(wp), INTENT(IN) :: geofac_div(:, :, :)

    ! cell based variable in which divergence is stored, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(inout) :: div_vec_c(:, :, :)

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

    iidx => cell_edge_idx
    iblk => cell_edge_blk

! loop through all patch cells (and blocks)
!

    !$ACC DATA PRESENT(vec_e, div_vec_c, geofac_div, iidx, iblk) IF(lzacc)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      ! original comment for divergence computation;
      ! everything that follows in this explanation has been combined into geofac_div

      ! compute the discrete divergence for cell jc by finite volume
      ! approximation (see Bonaventura and Ringler MWR 2005);
      ! multiplication of the normal vector component vec_e at the edges
      ! by the appropriate cell based edge_orientation is required to
      ! obtain the correct value for the application of Gauss theorem
      ! (which requires the scalar product of the vector field with the
      ! OUTWARD pointing unit vector with respect to cell jc; since the
      ! positive direction for the vector components is not necessarily
      ! the outward pointing one with respect to cell jc, a correction
      ! coefficient (equal to +-1) is necessary, given by
      ! ptr_patch%grid%cells%edge_orientation)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc, jk, jb) = &
            vec_e(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))*geofac_div(jc, 1, jb) + &
            vec_e(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))*geofac_div(jc, 2, jb) + &
            vec_e(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))*geofac_div(jc, 3, jb)

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC END DATA

  END SUBROUTINE div3d_lib

  SUBROUTINE div3d_2field_lib(vec_e, cell_edge_idx, cell_edge_blk, &
    &                         geofac_div, div_vec_c, in2, out2, &
    &                         i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                         slev, elev, nproma, lacc)

    ! edge based variable of which divergence is computed, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: vec_e(:, :, :)

    ! line indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_idx(:, :, :)

    ! block indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_blk(:, :, :)

    ! factor for divergence, dim: (nproma,cell_type,nblks_c)
    REAL(wp), INTENT(IN) :: geofac_div(:, :, :)

    ! cell based variable in which divergence is stored, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(inout) :: div_vec_c(:, :, :)

    ! second input field for more efficient processing in NH core, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: in2(:, :, :)

    ! second output field, dim: (nproma,nlev,nblks_c)
    REAL(wp), OPTIONAL, INTENT(inout) :: out2(:, :, :)

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
!
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_edge_idx
    iblk => cell_edge_blk

! loop through all patch cells (and blocks)
!
    !$ACC DATA PRESENT(vec_e, in2, div_vec_c, out2, geofac_div, iidx, iblk) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      ! original comment for divergence computation;
      ! everything that follows in this explanation has been combined into geofac_div

      ! compute the discrete divergence for cell jc by finite volume
      ! approximation (see Bonaventura and Ringler MWR 2005);
      ! multiplication of the normal vector component vec_e at the edges
      ! by the appropriate cell based edge_orientation is required to
      ! obtain the correct value for the application of Gauss theorem
      ! (which requires the scalar product of the vector field with the
      ! OUTWARD pointing unit vector with respect to cell jc; since the
      ! positive direction for the vector components is not necessarily
      ! the outward pointing one with respect to cell jc, a correction
      ! coefficient (equal to +-1) is necessary, given by
      ! ptr_patch%grid%cells%edge_orientation)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
          !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc, jk, jb) = &
            vec_e(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))*geofac_div(jc, 1, jb) + &
            vec_e(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))*geofac_div(jc, 2, jb) + &
            vec_e(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))*geofac_div(jc, 3, jb)

          out2(jc, jk, jb) = &
            in2(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))*geofac_div(jc, 1, jb) + &
            in2(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))*geofac_div(jc, 2, jb) + &
            in2(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))*geofac_div(jc, 3, jb)

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC END DATA

  END SUBROUTINE div3d_2field_lib

!-------------------------------------------------------------------------
!
!
!>
!! Special version of div that processes 4D fields in one step
!!
!! See standard routine (div3d) for further description
!!
  SUBROUTINE div4d_lib(cell_edge_idx, cell_edge_blk, &
    &                  geofac_div, f4din, f4dout, dim4d, &
    &                  i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                  slev, elev, nproma, lacc)

    ! line indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_idx(:, :, :)

    ! block indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_blk(:, :, :)

    ! factor for divergence, dim: (nproma,cell_type,nblks_c)
    REAL(wp), INTENT(IN) :: geofac_div(:, :, :)

    ! edge based 4D input field of which divergence is computed, dim: (nproma,nlev,nblks_e,dim4d)
    REAL(wp), INTENT(IN) :: f4din(:, :, :, :)

    ! cell based 4D output field in which divergence is stored, dim: (nproma,nlev,nblks_c,dim4d)
    REAL(vp), INTENT(inout) :: f4dout(:, :, :, :)

    ! last dimension of the input/output fields
    INTEGER, INTENT(IN) :: dim4d

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev(:)

    ! vertical end level
    INTEGER, INTENT(IN) :: elev(:)

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jc, jk, jb, ji
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_edge_idx
    iblk => cell_edge_blk

! loop through all patch cells (and blocks)
!

    !$ACC DATA PRESENT(f4din, f4dout, geofac_div, iidx, iblk) &
    !$ACC   COPY(slev, elev) IF(lzacc)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,ji) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
        DO ji = 1, dim4d
          !$ACC LOOP VECTOR
          DO jk = slev(ji), elev(ji)
#else
      !$ACC LOOP GANG
      DO ji = 1, dim4d
        DO jk = slev(ji), elev(ji)
          !$ACC LOOP VECTOR
          DO jc = i_startidx, i_endidx
#endif

            f4dout(jc, jk, jb, ji) = &
              f4din(iidx(jc, jb, 1), jk, iblk(jc, jb, 1), ji)*geofac_div(jc, 1, jb) + &
              f4din(iidx(jc, jb, 2), jk, iblk(jc, jb, 2), ji)*geofac_div(jc, 2, jb) + &
              f4din(iidx(jc, jb, 3), jk, iblk(jc, jb, 3), ji)*geofac_div(jc, 3, jb)

          END DO
        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC END DATA

!IF(ltimer) CALL timer_stop(timer_div)

  END SUBROUTINE div4d_lib

!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field.
!!
!! Computes discrete divergence of a vector field
!! given by its components in the directions normal to triangle edges,
!! followed by bilinear averaging to remove checkerboard noise
!! (Combines div_midpoint and cell_avg_varwgt to increase computing efficiency)
!!
  SUBROUTINE div_avg_lib(vec_e, cell_neighbor_idx, cell_neighbor_blk, cell_edge_idx, cell_edge_blk, &
   &                     geofac_div, avg_coeff, div_vec_c, opt_in2, opt_out2, &
   &                     i_startblk_in, i_endblk_in, i_startidx_in, i_endidx_in, &
   &                     nlev, nblks_c, patch_id, l_limited_area, slev, elev, nproma, l2fields, lacc)

    ! edge based variable of which divergence is computed, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) ::  vec_e(:, :, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN)  :: cell_neighbor_idx(:, :, :)

    ! block indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN)  :: cell_neighbor_blk(:, :, :)

    ! line indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_idx(:, :, :)

    ! block indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_blk(:, :, :)

    ! factor for quad-cell divergence  (nproma,4,nblks_e)
    REAL(wp), INTENT(IN) :: geofac_div(:, :, :)

    ! averaging coefficients dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(IN) :: avg_coeff(:, :, :)

    ! cell based variable in which divergence is stored, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(inout) :: div_vec_c(:, :, :)

    ! optional second input field for more efficient processing in NH core, dim: (nproma,nlev,nblks_e)
    REAL(wp), OPTIONAL, INTENT(IN) :: opt_in2(:, :, :)

    ! optional second output field, dim: (nproma,nlev,nblks_c)
    REAL(wp), OPTIONAL, INTENT(inout) :: opt_out2(:, :, :) ! optional second output field
    ! dim: (nproma,nlev,nblks_c)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk_in(3)

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk_in(3)

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in(3)

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in(3)

    ! total vertical level
    INTEGER, INTENT(IN) :: nlev

    ! number of blocks for the cells: t_patch%nblks_c
    INTEGER, INTENT(IN) :: nblks_c

    ! domain ID of current domain
    INTEGER, INTENT(IN) :: patch_id

    ! patch on which computation is performed
    LOGICAL, INTENT(IN) :: l_limited_area

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! whether the second field is present or not
    LOGICAL, INTENT(IN) :: l2fields

    ! if true, use OpenAcc
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jc, jk, jb
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    LOGICAL :: lzacc

    REAL(wp), DIMENSION(nproma, nlev, nblks_c) :: aux_c, aux_c2

    INTEGER, DIMENSION(:, :, :), POINTER :: inidx, inblk, ieidx, ieblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    inidx => cell_neighbor_idx
    inblk => cell_neighbor_blk
    ieidx => cell_edge_idx
    ieblk => cell_edge_blk

! First compute divergence
!
!$ACC DATA PRESENT(vec_e, avg_coeff, div_vec_c) CREATE(aux_c) &
!$ACC   PRESENT(geofac_div, ieidx, ieblk, inidx, inblk) IF(lzacc)
!$ACC DATA PRESENT(opt_in2, opt_out2) CREATE(aux_c2) IF(l2fields)
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk = i_startblk_in(1)
    i_endblk = i_endblk_in(1)

    IF (l2fields) THEN

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c_lib(i_startidx_in(1), i_endidx_in(1), nproma, jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)

#ifdef __LOOP_EXCHANGE
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jc = i_startidx, i_endidx
          DO jk = slev, elev
#else
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
#endif

            aux_c(jc, jk, jb) = &
              vec_e(ieidx(jc, jb, 1), jk, ieblk(jc, jb, 1))*geofac_div(jc, 1, jb) + &
              vec_e(ieidx(jc, jb, 2), jk, ieblk(jc, jb, 2))*geofac_div(jc, 2, jb) + &
              vec_e(ieidx(jc, jb, 3), jk, ieblk(jc, jb, 3))*geofac_div(jc, 3, jb)

            aux_c2(jc, jk, jb) = &
              opt_in2(ieidx(jc, jb, 1), jk, ieblk(jc, jb, 1))*geofac_div(jc, 1, jb) + &
              opt_in2(ieidx(jc, jb, 2), jk, ieblk(jc, jb, 2))*geofac_div(jc, 2, jb) + &
              opt_in2(ieidx(jc, jb, 3), jk, ieblk(jc, jb, 3))*geofac_div(jc, 3, jb)

          END DO
        END DO
        !$ACC END PARALLEL

      END DO
!$OMP END DO

    ELSE

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c_lib(i_startidx_in(1), i_endidx_in(1), nproma, jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jc = i_startidx, i_endidx
          DO jk = slev, elev
#else
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
#endif
            aux_c(jc, jk, jb) = &
              vec_e(ieidx(jc, jb, 1), jk, ieblk(jc, jb, 1))*geofac_div(jc, 1, jb) + &
              vec_e(ieidx(jc, jb, 2), jk, ieblk(jc, jb, 2))*geofac_div(jc, 2, jb) + &
              vec_e(ieidx(jc, jb, 3), jk, ieblk(jc, jb, 3))*geofac_div(jc, 3, jb)

          END DO
        END DO
        !$ACC END PARALLEL

      END DO
!$OMP END DO

    END IF

    IF (l_limited_area .OR. patch_id > 1) THEN
      ! Fill div_vec_c along the lateral boundaries of nests

      i_startblk = i_startblk_in(2)
      i_endblk = i_endblk_in(2)

      !$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk ! like copy(aux_c, div_vec_c)

        CALL get_indices_c_lib(i_startidx_in(2), i_endidx_in(2), nproma, jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
            div_vec_c(jc, jk, jb) = aux_c(jc, jk, jb)
          END DO
        END DO
        !$ACC END PARALLEL
      END DO
      !$OMP END DO

      IF (l2fields) THEN
        !$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk, i_endblk ! like copy(aux_c2, opt_out2)

          CALL get_indices_c_lib(i_startidx_in(2), i_endidx_in(2), nproma, jb, i_startblk, i_endblk, &
                                 i_startidx, i_endidx)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = slev, elev
            DO jc = i_startidx, i_endidx
              opt_out2(jc, jk, jb) = aux_c2(jc, jk, jb)
            END DO
          END DO
          !$ACC END PARALLEL
        END DO
        !$OMP END DO
      END IF
    END IF

!
! Now do averaging with weights given by avg_coeff

! values for the blocking
    i_startblk = i_startblk_in(3)
    i_endblk = i_endblk_in(3)

! loop through all patch cells (and blocks)
!

    IF (l2fields) THEN

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c_lib(i_startidx_in(3), i_endidx_in(3), nproma, jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jc = i_startidx, i_endidx
          DO jk = slev, elev
#else
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
#endif

            !  calculate the weighted average
            div_vec_c(jc, jk, jb) =  &
              &    aux_c(jc, jk, jb)*avg_coeff(jc, 1, jb) &
              &  + aux_c(inidx(jc, jb, 1), jk, inblk(jc, jb, 1))*avg_coeff(jc, 2, jb) &
              &  + aux_c(inidx(jc, jb, 2), jk, inblk(jc, jb, 2))*avg_coeff(jc, 3, jb) &
              &  + aux_c(inidx(jc, jb, 3), jk, inblk(jc, jb, 3))*avg_coeff(jc, 4, jb)

            opt_out2(jc, jk, jb) =  &
              &    aux_c2(jc, jk, jb)*avg_coeff(jc, 1, jb) &
              &  + aux_c2(inidx(jc, jb, 1), jk, inblk(jc, jb, 1))*avg_coeff(jc, 2, jb) &
              &  + aux_c2(inidx(jc, jb, 2), jk, inblk(jc, jb, 2))*avg_coeff(jc, 3, jb) &
              &  + aux_c2(inidx(jc, jb, 3), jk, inblk(jc, jb, 3))*avg_coeff(jc, 4, jb)

          END DO !cell loop
        END DO !vertical levels loop
        !$ACC END PARALLEL

      END DO !block loop

!$OMP END DO NOWAIT

    ELSE

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c_lib(i_startidx_in(3), i_endidx_in(3), nproma, jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jc = i_startidx, i_endidx
          DO jk = slev, elev
#else
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
#endif

            !  calculate the weighted average
            div_vec_c(jc, jk, jb) =  &
              &    aux_c(jc, jk, jb)*avg_coeff(jc, 1, jb) &
              &  + aux_c(inidx(jc, jb, 1), jk, inblk(jc, jb, 1))*avg_coeff(jc, 2, jb) &
              &  + aux_c(inidx(jc, jb, 2), jk, inblk(jc, jb, 2))*avg_coeff(jc, 3, jb) &
              &  + aux_c(inidx(jc, jb, 3), jk, inblk(jc, jb, 3))*avg_coeff(jc, 4, jb)

          END DO !cell loop
        END DO !vertical levels loop
        !$ACC END PARALLEL

      END DO !block loop

!$OMP END DO NOWAIT

    END IF

!$OMP END PARALLEL
!$ACC WAIT
!$ACC END DATA
!$ACC END DATA

  END SUBROUTINE div_avg_lib

!-------------------------------------------------------------------------
!
!>
!! Computes discrete rotation.
!!
!! Computes discrete rotation at
!! vertices of triangle cells (centers of dual hexagon cells)
!! from a vector field given by its components in the directions normal
!! to triangle edges.
!! input:  lives on edges (velocity points)
!! output: lives on dual of cells (vertices for triangular grid)
!!
  SUBROUTINE rot_vertex_atmos_lib(vec_e, vert_edge_idx, vert_edge_blk, &
    &                             geofac_rot, rot_vec, &
    &                             i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                             slev, elev, nproma, lacc)

    ! edge based variable of which rotation is computed, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: vec_e(:, :, :)

    ! line indices of edges around a vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_edge_idx(:, :, :)

    ! block indices of edges around a vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_edge_blk(:, :, :)

    ! factor for rotation (nproma,9-cell_type,nblks_v)
    REAL(wp), INTENT(IN) :: geofac_rot(:, :, :)

    ! vertex based variable in which rotation is stored, dim: (nproma,nlev,nblks_v)
    REAL(wp), INTENT(inout) :: rot_vec(:, :, :)

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
    
    INTEGER :: jv, jk, jb
    
    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc
    
    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)
    
    iidx => vert_edge_idx
    iblk => vert_edge_blk

!$ACC DATA PRESENT(vec_e, rot_vec, geofac_rot, iidx, iblk) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
    
      CALL get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)
    
      !
      ! compute the discrete rotation for vertex jv by
      ! finite volume approximation
      ! (see Bonaventura and Ringler MWR 2005);
      ! multiplication of the vector component vec_e by
      ! the appropriate dual cell based verts%edge_orientation
      ! is required to obtain the correct value for the
      ! application of Stokes theorem (which requires the scalar
      ! product of the vector field with the tangent unit vectors
      ! going around dual cell jv COUNTERCLOKWISE;
      ! since the positive direction for the vec_e components is
      ! not necessarily the one yielding counterclockwise rotation
      ! around dual cell jv, a correction coefficient (equal to +-1)
      ! is necessary, given by g%verts%edge_orientation
      !
    
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jv = i_startidx, i_endidx
        DO jk = slev, elev
#else
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jv = i_startidx, i_endidx
#endif
          !
          ! calculate rotation, i.e.
          ! add individual edge contributions to rotation
          ! (remark: for pentagon points the 6th weighting is 0)
          !
          
          rot_vec(jv, jk, jb) = &
            vec_e(iidx(jv, jb, 1), jk, iblk(jv, jb, 1))*geofac_rot(jv, 1, jb) + &
            vec_e(iidx(jv, jb, 2), jk, iblk(jv, jb, 2))*geofac_rot(jv, 2, jb) + &
            vec_e(iidx(jv, jb, 3), jk, iblk(jv, jb, 3))*geofac_rot(jv, 3, jb) + &
            vec_e(iidx(jv, jb, 4), jk, iblk(jv, jb, 4))*geofac_rot(jv, 4, jb) + &
            vec_e(iidx(jv, jb, 5), jk, iblk(jv, jb, 5))*geofac_rot(jv, 5, jb) + &
            vec_e(iidx(jv, jb, 6), jk, iblk(jv, jb, 6))*geofac_rot(jv, 6, jb)

        END DO
 
      END DO
      !$ACC END PARALLEL
 
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!$ACC END DATA

  END SUBROUTINE rot_vertex_atmos_lib

!>
!! Same as above routine, but expects reversed index order (vertical first)
!! of the output field if __LOOP_EXCHANGE is specified. In addition, the
!! output field (vorticity) has single precision if __MIXED_PRECISION is specified
!!
!!
  SUBROUTINE rot_vertex_ri_lib(vec_e, vert_edge_idx, vert_edge_blk, &
    &                          geofac_rot, rot_vec, &
    &                          i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                          slev, elev, nproma, lacc, acc_async)
!

!
    ! edge based variable of which rotation is computed, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: vec_e(:, :, :)

    ! line indices of edges around a vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_edge_idx(:, :, :)

    ! block indices of edges around a vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_edge_blk(:, :, :)

    ! factor for rotation (nproma,9-cell_type,nblks_v)
    REAL(wp), INTENT(IN) :: geofac_rot(:, :, :)

    ! vertex based variable in which rotation is stored, dim: (nproma,nlev,nblks_v)
    REAL(vp), INTENT(inout) :: rot_vec(:, :, :)

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

    INTEGER :: jv, jk, jb

    INTEGER :: i_startidx, i_endidx
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => vert_edge_idx
    iblk => vert_edge_blk

!  loop through over all patch vertices (and blocks)
!

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)

      ! calculate rotation, i.e.
      ! add individual edge contributions to rotation
      !
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jv = i_startidx, i_endidx
        DO jk = slev, elev
          rot_vec(jk, jv, jb) = &
#else
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jv = i_startidx, i_endidx
          rot_vec(jv, jk, jb) = &
#endif
            vec_e(iidx(jv, jb, 1), jk, iblk(jv, jb, 1))*geofac_rot(jv, 1, jb) + &
            vec_e(iidx(jv, jb, 2), jk, iblk(jv, jb, 2))*geofac_rot(jv, 2, jb) + &
            vec_e(iidx(jv, jb, 3), jk, iblk(jv, jb, 3))*geofac_rot(jv, 3, jb) + &
            vec_e(iidx(jv, jb, 4), jk, iblk(jv, jb, 4))*geofac_rot(jv, 4, jb) + &
            vec_e(iidx(jv, jb, 5), jk, iblk(jv, jb, 5))*geofac_rot(jv, 5, jb) + &
            vec_e(iidx(jv, jb, 6), jk, iblk(jv, jb, 6))*geofac_rot(jv, 6, jb)

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (PRESENT(acc_async)) THEN
      IF (.NOT. acc_async) THEN
        !$ACC WAIT
      END IF
    ELSE
      !$ACC WAIT
    END IF

  END SUBROUTINE rot_vertex_ri_lib

END MODULE mo_lib_divrot
