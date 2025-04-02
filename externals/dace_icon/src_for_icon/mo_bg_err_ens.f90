!
!+ Operations on the VarEnKF ensemble covariance matrix
!
MODULE mo_bg_err_ens
!
! Description:
!   Operations on the VarEnKF ensemble covariance matrix,
!   called from generic (NMC + EnKF) routines from module mo_bg_err_2d.
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_37        2014-12-23 Andreas Rhodin
!  New routines for Variational Ensemble Kalman Filter (VarEnKF)
! V1_42        2015-06-08 Andreas Rhodin
!  optimisation for EnKF/PSAS; temporal interpolation for COSMO MEC
! V1_46        2016-02-05 Andreas Rhodin
!  base decisions on new flag 'vct'; fix for passive SYNOP t2m observations
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!------------------------------------------------------------------------------
  !=============
  ! modules used
  !=============
  !-----------------
  ! general routines
  !-----------------
  use mo_kind,       only: wp, i8          ! precision kind parameter
  use mo_exception,  only: finish          ! abort on error condition
  use mo_cpu_time,   only: stop_time       ! determine cpu and wall time
  !--------------------------------------
  ! parallel matrix and vector operations
  !--------------------------------------
  use mo_dec_matrix, only: t_vector,      &! partitioned vector data type
                           assignment(=)   ! vector assignment
  !------------------
  ! atmospheric state
  !------------------
  use mo_t_col,      only: t_cols,        &! data type to hold columns
                           get_cols,      &! redistribute model columns
                           get_cols_adj,  &! adjoint routine
                           alloc_cols,    &! allocate components of t_cols
                           dealloc_cols,  &! deallocate t_cols components
                           assignment(=), &! assign t_cols
                           COL_T, COL_UV, &! identifier for temp., hor.wind,
                           COL_RH,        &!   rel.humidity
                           COL_TV2,       &!   virtual temperature
                           COL_GEOH,      &!   geopot.height (half levels)
                           COL_GEO, COL_P  !   geopot.height (full levels)
  use mo_atm_state,  only: t_atm,         &! atmospheric state derived type
                           assignment(=), &! t_atm = real
                           construct,     &! set up t_atm
                           destruct        ! deallocate t_atm components
  use mo_atm_grid,   only: VCT_P_HYB       ! GME/HRM Vertical coordinate
  use mo_wmo_tables, only: WMO3_ISOBARIC   ! level type
  use mo_physics,    only: gacc,          &! gravity acceleration
                           R               ! gas constant of dry air
  !-----------------------------------
  ! background error covariance module
  !-----------------------------------
  use mo_bg_err_io,  only: t_box,         &! derived type definition
                           t_intop,       &! interpolation operator data type
                           construct,     &! set up type t_intop
                           destruct,      &! deallocate t_intop components
                           set_vic_ps,    &! set vert.interp.coeffs. for ps
                           set_vic_atm,   &! set vert.interp.coeffs. for atm.
                           IL_H, IL_RH, IL_U,   IL_V, &! positions in xl(,:,)
                           IV_TV,IL_TV, IV_U,   IV_V, &!   in xv(,:,) (Ens-B)
                           IV_H, IV_RH, IV_PSI, IV_CHI !   in xv(,:,) (NMC_B)
  use mo_varenkf,    only: Benkf,         &! Ensemble background error Cov.M.
                           apply_varenkf, &! apply B_enkf to atmospheric state
                           v_lift_synth,  &! vertical wavelet synthesis
                           v_lift_adj,    &! adjoint vertical wavelet synthesis
                           w_ens_b         ! weight of ensemble B
  !-------------
  ! observations
  !-------------
  use mo_t_obs,      only: OBS_TV,        &! virtual temperature  id
                           OBS_T,         &! temperature          id
                           OBS_RH,        &! relative humidity    id
                           OBS_H,         &! geopotential height  id
                           OBS_HS,        &! surface geopotential id
                           OBS_U,         &! u-wind               id
                           OBS_V,         &! v-wind               id
                           OBS_DUM         ! dummy type
  implicit none
  !================
  ! public entities
  !================
  private
  public :: apply_Hi_ens       ! apply Io        (EnKF model -> obs. input)
  public :: apply_B_mi_ensb    ! apply Im B Io^t (obs. input -> det. model)
  public :: apply_B_ii_ensb    ! apply Io B Io^t (obs. input -> obs. input)
!==============================================================================
contains
!==============================================================================
  subroutine apply_B_ii_ensb (x, y)
  type (t_vector)  ,intent(in)    :: x     ! input
  type (t_vector)  ,intent(inout) :: y     ! output

    real(wp) ,allocatable  :: xv(:,:,:) !  'vertical mode' space
    type(t_atm)  ,target   :: z         ! gradient           on coarse grid
    type(t_atm)  ,target   :: a         ! analysis increment on coarse grid
    type(t_atm)  ,pointer  :: zz        ! gradient           on coarse grid
    type(t_atm)  ,pointer  :: aa        ! analysis increment on coarse grid
    type (t_cols)          :: cols      ! model columns
    integer                :: ic        ! column index
    integer(i8)            :: ids
    integer                :: l         ! index for optimised PSAS columns

    if (w_ens_b == 0._wp) return

    !---------------------
    ! allocate temporaries
    !---------------------
    l = Benkf% l_set
    allocate (xv (Benkf% grid% nz, Benkf% io% np_v, Benkf% io% nc_v))
    call construct (z, template = Benkf% mean(l), alloc ='u v rh tv geof')
    call construct (a, template = Benkf% mean(l), alloc ='u v rh tv geof')

    if (Benkf% n_vert == 0) then
      aa => a
      zz => z
    else
      allocate (aa)
      allocate (zz)
      call construct (zz, template = Benkf% mean(l), nz = Benkf% j_vert,&
                             alloc ='u v rh tv geof'                    )
      call construct (aa, template = Benkf% mean(l), nz = Benkf% j_vert,&
                             alloc ='u v rh tv geof'                    )
    endif

    !----------------------------------------------------------
    ! adjoint interpolation (coarse model grid <- observations)
    !----------------------------------------------------------
    call apply_IhIv_ensb_t (Benkf% io, xv, x)

    !----------------------
    ! adjoint transposition
    !----------------------
    ids = COL_UV+COL_RH+COL_TV2+COL_GEO
    if (Benkf% l_set == 1) then
      call alloc_cols (cols, ncol=Benkf% io% nc_v, ke=Benkf% grid% nz,ids=ids)
      do ic = 1, cols% ncol
        cols% col(ic)% u   = xv (:, IV_U , ic)
        cols% col(ic)% v   = xv (:, IV_V , ic)
        cols% col(ic)% rh  = xv (:, IV_RH, ic)
        cols% col(ic)% geo = xv (:, IV_H , ic)
        cols% col(ic)% tv  = xv (:, IV_TV, ic)
      end do
      z = 0._wp
      call get_cols_adj (Benkf% io% mc, z, cols, ids)
      call dealloc_cols (cols)
    else
      do ic = 1, Benkf% io% nc_v
        z% u    (ic,1,:,1) = xv (:, IV_U , ic)
        z% v    (ic,1,:,1) = xv (:, IV_V , ic)
        z% rh   (ic,1,:,1) = xv (:, IV_RH, ic)
        z% geof (ic,1,:,1) = xv (:, IV_H , ic)
        z% tv   (ic,1,:,1) = xv (:, IV_TV, ic)
      end do
    endif

    !-------------------------------------
    ! vertical wavelet transform (adjoint)
    !-------------------------------------
    if (Benkf% n_vert > 0) call v_lift_adj (zz, z)

    !------------------------------------------
    ! localised ensemble B on coarse model grid
    !------------------------------------------
    call apply_varenkf (Benkf, aa, zz, l, l)

    !---------------------------------------
    ! vertical wavelet transform (synthesis)
    !---------------------------------------
    if (Benkf% n_vert > 0) call v_lift_synth (aa, a)

    !--------------
    ! transposition
    !--------------
    if (Benkf% l_set == 1) then
      call get_cols (Benkf% io% mc, a, cols, ids)
      do ic = 1, cols% ncol
        xv (:, IV_U , ic) = cols% col(ic)% u
        xv (:, IV_V , ic) = cols% col(ic)% v
        xv (:, IV_RH, ic) = cols% col(ic)% rh
        xv (:, IV_H , ic) = cols% col(ic)% geo
        xv (:, IV_TV, ic) = cols% col(ic)% tv
      end do
      call dealloc_cols (cols)
    else
      do ic = 1, Benkf% io% nc_v
        xv (:, IV_U , ic) = a% u    (ic,1,:,1)
        xv (:, IV_V , ic) = a% v    (ic,1,:,1)
        xv (:, IV_RH, ic) = a% rh   (ic,1,:,1)
        xv (:, IV_H , ic) = a% geof (ic,1,:,1)
        xv (:, IV_TV, ic) = a% tv   (ic,1,:,1)
      end do
    endif

    !-------------------------
    ! scale with B_ensb weight
    !-------------------------
    xv = xv * w_ens_b

    !--------------------------------------------------
    ! interpolation (coarse model grid -> observations)
    !--------------------------------------------------
    call apply_IhIv_ensb (Benkf% io, xv, y)

    !---------
    ! clean up
    !---------
    call destruct     (a)
    call destruct     (z)
    if (Benkf% n_vert /= 0) then
      call destruct   (aa)
      call destruct   (zz)
      deallocate      (aa)
      deallocate      (zz)
    endif

  end subroutine apply_B_ii_ensb

!------------------------------------------------------------------------------
  subroutine apply_Hi_ens (x, y, io)
  type (t_atm)     ,intent(inout) :: x     ! input  (ensemble model space)
  type (t_vector)  ,intent(inout) :: y     ! output (interpolation space)
  type(t_intop)    ,intent(in)    :: io    ! interpolation coefficients
  optional                        :: io    !        (default: Benkf% io)
  !----------------------------------------------------------------
  ! apply (linear) interpolation operator to (ensemble) model state
  !----------------------------------------------------------------
    real(wp) ,allocatable  :: xv(:,:,:) !  'vertical mode' space
    type (t_cols)          :: cols      ! model columns
    integer                :: ic        ! column index
    integer(i8)            :: ids
    type(t_atm) ,pointer   :: z
    target                 :: x

    !------------------
    ! VarEnkf enabled ?
    !------------------
    y = 0._wp
    if (w_ens_b == 0._wp) return

    !-----------------------------
    ! vertical wavelet transform ?
    !-----------------------------
    if (x% ub(3) == Benkf% grid% nz) then
      z => x
    else
      allocate (z)
      call construct (z, grid=x% grid)
      call v_lift_synth (x, z)
    endif

    !---------------------
    ! allocate temporaries
    !---------------------
    if (present (io)) then
      allocate (xv (Benkf% grid% nz,        io% np_v,        io% nc_v))
    else
      allocate (xv (Benkf% grid% nz, Benkf% io% np_v, Benkf% io% nc_v))
    end if

    !--------------
    ! transposition
    !--------------
    ids = COL_UV+COL_RH+COL_TV2+COL_GEO
    if (present (io)) then
      call get_cols (       io% mc, z, cols, ids)
    else
      call get_cols (Benkf% io% mc, z, cols, ids)
    end if
    do ic = 1, cols% ncol
      xv (:, IV_U , ic) = cols% col(ic)% u
      xv (:, IV_V , ic) = cols% col(ic)% v
      xv (:, IV_RH, ic) = cols% col(ic)% rh
      xv (:, IV_H , ic) = cols% col(ic)% geo
      xv (:, IV_TV, ic) = cols% col(ic)% tv
    end do

    !--------------------------------------------------
    ! interpolation (coarse model grid -> observations)
    !--------------------------------------------------
    if (present (io)) then
      call apply_IhIv_ensb (       io, xv, y)
    else
      call apply_IhIv_ensb (Benkf% io, xv, y)
    end if

    !---------
    ! clean up
    !---------
    call dealloc_cols (cols)
    if (x% ub(3) /= Benkf% grid% nz) then
      call destruct (z)
      deallocate    (z)
    endif

  end subroutine apply_Hi_ens

!------------------------------------------------------------------------------
  subroutine apply_B_mi_ensb (a_m, cbg, x, lnewpl, e_fi)
  type(t_cols)     ,intent(inout) :: a_m(:) ! analysis increment (model space)
  type(t_cols)     ,intent(in)    :: cbg(:) ! reference state
  type (t_vector)  ,intent(in)    :: x      ! input
  logical          ,intent(in)    :: lnewpl ! analysis increment on new p-levs
  type (t_vector)  ,intent(in)    :: e_fi   ! background error
  !------------------------------------------------------------------
  ! multiply a vector in interpolation space by the ensemble B matrix
  ! result is on model grid
  ! +++ work in progress +++
  !------------------------------------------------------------------

    real(wp) ,allocatable  :: xv (:,:,:) ! 'vertical mode' space
    real(wp) ,allocatable  :: logp (:,:) ! EnKF background pressure levels
    type(t_atm)  ,target   :: z          ! gradient           on coarse grid
    type(t_atm)  ,target   :: a          ! analysis increment on coarse grid
    type(t_atm)  ,pointer  :: zz         ! gradient           on coarse grid
    type(t_atm)  ,pointer  :: aa         ! analysis increment on coarse grid
    type (t_cols)          :: cols       ! model columns
    type (t_cols)          :: colp
    integer                :: ic         ! column index
    integer(i8)            :: ids        ! bit flag field for variables
    type (t_intop)         :: io         ! interpolation coefficients
    integer                :: l          ! index for optimised PSAS columns

    if (w_ens_b == 0._wp) return

    !---------------------
    ! allocate temporaries
    !---------------------
    l = Benkf% l_set
    allocate (xv   (Benkf% grid% nz, Benkf% io% np_v, Benkf% io% nc_v))
    call construct (z, template = Benkf% mean(l), alloc ='u v rh tv geof')
    call construct (a, template = Benkf% mean(1), alloc ='u v rh tv geof')

    if (Benkf% n_vert == 0) then
      aa => a
      zz => z
    else
      allocate (aa)
      allocate (zz)
      call construct (zz, template = Benkf% mean(l), nz = Benkf% j_vert,&
                             alloc = 'u v rh tv geof'                   )
      call construct (aa, template = Benkf% mean(1), nz = Benkf% j_vert,&
                             alloc = 'u v rh tv geof'                   )
    endif

    !----------------------------------------------------------
    ! adjoint interpolation (coarse model grid <- observations)
    !----------------------------------------------------------
    call apply_IhIv_ensb_t (Benkf% io, xv, x)

    !----------------------
    ! adjoint transposition
    !----------------------
    ids = COL_UV+COL_RH+COL_TV2+COL_GEO
    if (Benkf% l_set == 1) then
      call alloc_cols (cols, ncol=Benkf% io% nc_v, ke=Benkf% grid% nz,ids=ids)
      do ic = 1, cols% ncol
        cols% col(ic)% u   = xv (:, IV_U , ic)
        cols% col(ic)% v   = xv (:, IV_V , ic)
        cols% col(ic)% rh  = xv (:, IV_RH, ic)
        cols% col(ic)% geo = xv (:, IV_H , ic)
        cols% col(ic)% tv  = xv (:, IV_TV, ic)
      end do
      z = 0._wp
      call get_cols_adj (Benkf% io% mc, z, cols, ids)
      call dealloc_cols (cols)
    else
      do ic = 1, Benkf% io% nc_v
        z% u    (ic,1,:,1) = xv (:, IV_U , ic)
        z% v    (ic,1,:,1) = xv (:, IV_V , ic)
        z% rh   (ic,1,:,1) = xv (:, IV_RH, ic)
        z% geof (ic,1,:,1) = xv (:, IV_H , ic)
        z% tv   (ic,1,:,1) = xv (:, IV_TV, ic)
      end do
    endif
    deallocate (xv)

    !-------------------------------------
    ! vertical wavelet transform (adjoint)
    !-------------------------------------
    if (Benkf% n_vert > 0) call v_lift_adj (zz, z)

    !------------------------------------------
    ! localised ensemble B on coarse model grid
    !------------------------------------------
    call apply_varenkf (Benkf, aa, zz, l, 1)

    !---------------------------------------
    ! vertical wavelet transform (synthesis)
    !---------------------------------------
    if (Benkf% n_vert > 0) call v_lift_synth (aa, a)

    !---------------------------------------------
    ! set up horizontal interpolation coefficients
    ! and transposition
    !---------------------------------------------
    call stop_time ('construct (io_ensb)')
    call construct (io, cbg, Benkf% mean(1))

    !----------------------
    ! transposition
    !----------------------
    ids = COL_UV+COL_RH+COL_TV2+COL_GEO
    call get_cols (io% mc, a, cols, ids)
    allocate (xv   (Benkf% grid% nz, io% np_v, io% nc_v))
    allocate (logp (Benkf% grid% nz,           io% nc_v))

    do ic = 1, cols% ncol
      xv (:, IV_U , ic) = cols% col(ic)% u
      xv (:, IV_V , ic) = cols% col(ic)% v
      xv (:, IV_RH, ic) = cols% col(ic)% rh
      xv (:, IV_H , ic) = cols% col(ic)% geo
      xv (:, IV_TV, ic) = cols% col(ic)% tv
    end do
    call dealloc_cols (cols)

    !-------------------------------------------------
    ! get background pressure levels for interpolation
    !-------------------------------------------------
    ids = COL_P
    call get_cols (io% mc, Benkf% mean(1), colp, ids)
    do ic = 1, colp% ncol
      logp (:, ic) = log (colp% col(ic)% p)
    end do
    call dealloc_cols (colp)

    !-------------------------
    ! scale with B_ensb weight
    !-------------------------
    xv = xv * w_ens_b

    !--------------------------------------------------
    ! interpolation (coarse grid -> model grid)
    !--------------------------------------------------
    call apply_IhIv_ensb_m (xv, logp, io, a_m, cbg, lnewpl)

    !---------
    ! clean up
    !---------
    call destruct     (io)
    call destruct     (a)
    call destruct     (z)
    if (Benkf% n_vert /= 0) then
      call destruct   (aa)
      call destruct   (zz)
      deallocate      (aa)
      deallocate      (zz)
    endif

  end subroutine apply_B_mi_ensb
!------------------------------------------------------------------------------
  subroutine apply_IhIv_ensb_t (io, xv, xi)
  type(t_intop)  ,intent(in)    :: io       ! interpolation coefficients
  real(wp)       ,intent(out)   :: xv(:,:,:)! result   in 'vertical mode' space
  type(t_vector) ,intent(in)    :: xi       ! argument in 'interpolation' space
  !-----------------------------
  ! apply adjoint interpolations
  !-----------------------------
    integer  :: nc_l                    ! no. columns  at observation pts.
    integer  :: ib                      ! global 'box' index
    integer  :: ibl                     ! local  'box' index
    integer  :: ic                      ! column index (in xl)
    integer  :: icv                     ! column index (in xv)
!   integer  :: id                      ! sink index
    integer  :: il                      ! parameter index (in xl)
!   integer  :: j                       ! meridional index
    integer  :: l                       ! level index
    integer  :: i
    integer  :: k
    real(wp) :: xl (Benkf% grid% nz, &  ! columns at observation pts
                           io  % np_l)
    type(t_box)   ,pointer :: bx

    xv = 0._wp
    !----------------
    ! loop over boxes
    !----------------
    do ibl = 1, io% n_box
      bx   => io% bx(ibl)
      ib   =  bx% ib
      nc_l =  bx% nc_l
!     id   =  0

      !=====================================
      ! Iv_t: adjoint vertical interpolation
      !=====================================

      !------------------
      ! loop over columns
      !------------------
      do ic = 1, nc_l
        xl = 0._wp

        !-----------------
        ! loop over levels
        !-----------------
        do l = bx% hic(ic)% l% i + 1, bx% hic(ic)% l% i + bx% hic(ic)% l% n

          !-------------------------------------
          ! loop over observations at this level
          !-------------------------------------
          do i = bx% vic(l)% i% i + 1, &
                 bx% vic(l)% i% i + bx% vic(l)% i% n
            select case (bx% ip(i))
            case (OBS_H, OBS_HS)
              il = IL_H
            case (OBS_TV)
              il = IL_TV
            case (OBS_RH)
              il = IL_RH
            case (OBS_U)
              il = IL_U
            case (OBS_V)
              il = IL_V
            case (OBS_DUM)
            case default
              write(0,*)'apply_IhIv_ensb_t: invalid interpolation type:',&
                        bx% ip(i)
              write(0,*)'apply_IhIv_ensb_t: box,l,i=',ib,l,i
              call finish ('apply_IhIv_ensb_t','invalid interpolation type')
            end select
            select case (bx% ip(i))
!           case (OBS_TV)
!             do k = 1, covm% nwv
!               xl  (bx% vic(l)% ix+k-1 ,il) = &
!                 xl(bx% vic(l)% ix+k-1 ,il) - &
!                     xi% s(ib)% x   (i)     * &
!                    bx% vic(l)% wt(k)       * &
!                    bx% vic(l)% exzn              !+++ check exzn
!             end do
            case (OBS_DUM)
            case default
!             mr(ib)% lq(ic) = .true.
              do k = 1, bx% vic(l)% g% n
                xl (bx% vic(l)% g% i+k-1 ,il) = &
                 xl(bx% vic(l)% g% i+k-1 ,il) + &
                    bx% vic(l)% wh(k)         * &
                     xi% s(ib)% x   (i)
              end do
            end select
          end do

        end do

        !=======================================
        ! Ih_t: adjoint horizontal interpolation
        !=======================================

        !-----------------------------
        ! EnKF-B: direct interpolation
        !-----------------------------
        do k = 1, size(bx% hic(ic)% imc,1)
          icv = bx% hic(ic)% imc(k,1)
          if (icv==0) exit
          do il = 1,5
            xv   (:, il,   icv) = &
              xv (:, il,   icv) + &
              xl (:, il)        * &
              bx% hic(ic)% w (k)
          end do
        end do

        !-------------------------
        ! end of loop over columns
        !-------------------------
      end do

      !-----------------------
      ! end of loop over boxes
      !-----------------------
    end do
  end subroutine apply_IhIv_ensb_t
!------------------------------------------------------------------------------
  subroutine apply_IhIv_ensb (io, xv, xi)
  type(t_intop)  ,intent(in)    :: io       ! interpolation coefficients
  real(wp)       ,intent(in)    :: xv(:,:,:)! argument in 'vertical mode' space
  type(t_vector) ,intent(inout) :: xi       ! result   in 'interpolation' space
  !---------------------
  ! apply interpolations
  !---------------------
    integer  :: nc_l                    ! no. columns  at observation pts.
    integer  :: ib                      ! global 'box' index
    integer  :: ibl                     ! local  'box' index
    integer  :: ic                      ! column index (in xl)
    integer  :: icv                     ! column index (in xv)
!   integer  :: id                      ! sink index
    integer  :: il                      ! parameter index (in xl)
!   integer  :: j                       ! meridional index
    integer  :: l                       ! level index
    integer  :: i
    integer  :: k
    real(wp) :: xl (Benkf% grid% nz, &  ! columns at observation pts
                           io  % np_l)
    type(t_box)   ,pointer :: bx

!   xi  = 0._wp  ! add increments to previous ones from B_NMC
    !----------------
    ! loop over boxes
    !----------------
    do ibl = 1, io% n_box
      bx   => io% bx(ibl)
      ib   =  bx% ib
      nc_l =  bx% nc_l
!     id   =  0
      !------------------
      ! loop over columns
      !------------------
      do ic = 1, nc_l
        xl = 0._wp
        !=============================
        ! Ih: horizontal interpolation
        !=============================

        !-----------------------------
        ! EnKF-B: direct interpolation
        !-----------------------------
        do k = 1, size(bx% hic(ic)% imc,1)
          icv = bx% hic(ic)% imc(k,1)
          if (icv==0) exit
          do il = 1,5
            xl (:, il) = xl (:, il)       &
                       + xv (:, il,  icv) &
                       * bx% hic(ic)% w (k)
          end do
        end do

        !===========================
        ! Iv: vertical interpolation
        !===========================
        !-----------------
        ! loop over levels
        !-----------------
        do l = bx% hic(ic)% l% i + 1, bx% hic(ic)% l% i + bx% hic(ic)% l% n

          !-------------------------------------
          ! loop over observations at this level
          !-------------------------------------
          do i = bx% vic(l)% i% i + 1, &
                 bx% vic(l)% i% i + bx% vic(l)% i% n
            select case (bx% ip(i))
            case (OBS_H, OBS_HS)
              il = IL_H
            case (OBS_TV)
              il = IL_TV
            case (OBS_T)   !+++ currently required for passive t2m obs. +++!
              il = IL_TV   !+++ (first guess error estimate only)       +++!
            case (OBS_RH)
              il = IL_RH
            case (OBS_U)
              il = IL_U
            case (OBS_V)
              il = IL_V
            case (OBS_DUM)
            case default
              write(0,*)'apply_IhIv_ensb: invalid interpolation type:',&
                        bx% ip(i)
              write(0,*)'apply_IhIv_ensb: box,l,i=',ib,l,i
              call finish ('apply_IhIv_ensb','invalid interpolation type')
            end select
            select case (bx% ip(i))
!           case (OBS_TV)
!             do k = 1, covm% nwv
!               xi% s(ib)% x(i) = xi% s(ib)% x(i) - &
!                 xl(bx% vic(l)% ix+k-1 ,il)      * &
!                    bx% vic(l)% wt(k)            * &
!                    bx% vic(l)% exzn                  !+++ check exzn
!             end do
            case (OBS_DUM)
            case default
!             ml% lq(levs(l)% is) = .true.
              do k = 1, bx% vic(l)% g% n
                xi% s(ib)% x(i) = xi% s(ib)% x(i) + &
                  xl(bx% vic(l)% g% i+k-1 ,il)    * &
                  bx% vic(l)% wh(k)
              end do  ! k   (source levels)
            end select
          end do      ! i   (observations at this level)
        end do        ! l   (levels)

        !-------------------------
        ! end of loop over columns
        !-------------------------
      end do          ! ic  (columns)

      !-----------------------
      ! end of loop over boxes
      !-----------------------
    end do            ! ibl (boxes)
  end subroutine apply_IhIv_ensb
!------------------------------------------------------------------------------
  subroutine apply_IhIv_ensb_m (xw, logp, io, a_m, cbg, lnewpl)
  real(wp)     ,intent(in)    :: xw (:,:,:) ! input, vertical profiles
  real(wp)     ,intent(in)    :: logp (:,:) ! vertical coordinate (EnKF bg)
  type(t_intop),intent(inout) :: io         ! interpolation coefficients
  type(t_cols) ,intent(inout) :: a_m(:)     ! analysis increment (model space)
  type(t_cols) ,intent(in)    :: cbg(:)     ! reference state
  logical      ,intent(in)    :: lnewpl     ! analysis increment on new p-levs
  !-------------------------------------------------
  ! interpolation
  ! result is on model grid
  !-------------------------------------------------
    integer ,parameter :: nm  = 4             ! number of multi-level fields
    integer ,parameter :: ns  = 1             ! number of single-level fields
    integer            :: nvc = 0             ! number of variables/column
    integer            :: nvcf= 0             ! number of variables, full level
    integer            :: ke                  ! number of model levels
    integer            :: icv                 ! column index (in xv)
    integer            :: ic                  ! column index (in xl)
    integer            :: nc_l                ! no. columns at observation pts.
    integer            :: ib                  ! global 'box' index
!   integer            :: j                   ! meridional index
    integer            :: k                   !
    integer            :: l                   !
    integer            :: i                   !
    integer            :: il                  !
    real(wp)           :: ehs                 ! surface geopotential obs.error
    real(wp)           :: ps                  ! surface pressure
    real(wp)           :: xl (Benkf%grid% nz,&!
                              io% np_l)       ! columns at observation pts
    real(wp)           :: zl (Benkf%grid% nz) ! Benkf levels
    integer  ,allocatable :: t_int (:)        !
    real(wp) ,allocatable :: e     (:)        !
    real(wp) ,allocatable :: z     (:)
    real(wp) ,allocatable :: y     (:)
    real(wp) ,allocatable :: ph    (:)
    real(wp) ,allocatable :: lpf   (:)
    type (t_box) ,pointer :: bx               !

    call stop_time ('IhIv_ensb')
    !----------------
    ! loop over boxes
    !----------------
    do ib  = 1, size (cbg)
      bx   => io% bx(ib)
      nc_l = bx% nc_l

      !--------------------------------------
      ! allocate result variable, set to zero
      ! aready done by apply_IhIv_m !!!!!!!!!
      !--------------------------------------
!!!   call alloc_cols (a_m(ib), tmp=cbg(ib), ids=COL_T+COL_UV+COL_RH+COL_GEOH)
!!!   a_m(ib) = 0._wp

      !---------------------
      ! allocate temporaries
      !---------------------
      if (nvc == 0) then
        ke  = a_m(ib)% ke
        nvcf = nm * ke + ns
        nvc  = nvcf + ke
        allocate (t_int (nvc ))
        allocate (e     (nvc ))
        allocate (z     (nvc ))
        allocate (y     (nvc ))
        allocate (ph    (ke+1))
        allocate (lpf   (ke  ))
        !-----------------
        ! set up meta data
        !-----------------
        t_int (1        ) = OBS_H
        t_int (2:nvcf:nm) = OBS_TV
        t_int (3:nvcf:nm) = OBS_RH
        t_int (4:nvcf:nm) = OBS_U
        t_int (5:nvcf:nm) = OBS_V
        t_int (1+nvcf:  ) = OBS_H
      endif

      !------------------
      ! loop over columns
      !------------------
      do ic = 1, nc_l
        xl = 0._wp
        zl = 0._wp
        !-------------------------
        ! horizontal interpolation
        !-------------------------
        do k = 1, size(bx% hic(ic)% imc,1)
          icv = bx% hic(ic)% imc(k,1)
          if (icv==0) exit
!!!       j   = io% mc% c(icv)% idx(2)
          do il = 1,5
            xl (:, il)  = xl (:, il)       &
                        + xw (:, il,  icv) &
                        * bx% hic(ic)% w (k)
          end do
          zl   (:)      = zl   (:)      &
                        + logp (:, icv) &
                        * bx% hic(ic)% w (k)
        end do

        if (lnewpl) then
          call finish ('apply_IhIv_ensb_m','lnewpl not implemented for EnKF B')
          !---------------------------------------------------------------
          ! set up vertical interpolation coefficient for surface pressure
          !---------------------------------------------------------------
          call set_vic_ps (bx, ehs, cbg(ib)%col(ic))    !!! add GRID INFO

          !-----------------------------
          ! interpolate surface pressure
          !-----------------------------

          ps = 0
          do k = 1, bx% vic(1)% g% n
            ps = ps + xl(bx% vic(1)% g% i+k-1 ,IL_H)  * &
                         bx% vic(1)% wh(k)            * &
                         ehs
          end do
          a_m(ib)% col(ic)% s% ps = &
          a_m(ib)% col(ic)% s% ps + &
          ps * (gacc/R) * (cbg(ib)% col(ic)% s% ps / cbg(ib)% col(ic)% s% t2m)
          a_m(ib)% col(ic)% s% psr =                         &
          a_m(ib)% col(ic)% s% psr + cbg(ib)% col(ic)% s% ps &
                                   + a_m(ib)% col(ic)% s% ps
        else
          a_m(ib)% col(ic)% s% psr =                         &
          a_m(ib)% col(ic)% s% psr + cbg(ib)% col(ic)% s% ps
        endif

        !----------------------------------------------------------
        ! set up vertical interpolation coefficients for atmosphere
        !----------------------------------------------------------
        select case (a_m(ib)% levtyp)
        case (WMO3_ISOBARIC)
          !------------------------------------------------
          ! pressure levels: same for t,rh,u,v, geop.height
          !------------------------------------------------
          do k = 1, ke
            lpf (k) = log (a_m(ib)% ak(k))
            z (nm*(k-1)+2:nm*k+1) = lpf(k)     ! t,rh,u,v
            z (nvcf+k)            = lpf(k)     ! geop.height
          end do
        case default
          !------------------------------------
          ! GME/HRM hybrid pressure coordinates
          !------------------------------------
          if (a_m(ib)% vctype == VCT_P_HYB) then
            !----------------------------------------
            ! hybrid levels: full levels for t,rh,u,v
            !----------------------------------------
            do k = 1, ke + 1
              ph (k) = a_m(ib)% ak(k) + a_m(ib)% bk(k) * a_m(ib)% col(ic)% s% psr
            end do
            do k = 1, ke
              lpf (k) = log(0.5_wp * (ph (k) + ph (k+1)))
              z (nm*(k-1)+2:nm*k+1) = lpf(k)
            end do
            !----------------------------
            ! half levels for geop.height
            !----------------------------
            do k = 1, ke
              z (nvcf+k) = log (ph (k+1))
            end do
          !---------------------------++++++++++++++++++
          ! COSMO hybrid z-coordinates (and ICON so far)
          !---------------------------++++++++++++++++++
          else
            do k = 1, ke
              lpf (k) = log (cbg(ib)% col(ic)% p(k))
              z (nm*(k-1)+2:nm*k+1) = lpf(k)     ! t,rh,u,v
              z (nvcf+k)            = lpf(k)     ! geop.height
            end do
          endif
        end select
        z (1) = log (a_m(ib)% col(ic)% s% psr)
        call set_vic_atm (bx, e, z, t_int, a_m(ib)%col(ic), logp = zl)

        !------------------------------------------------
        ! interpolate atmospheric parameters (vertically)
        !------------------------------------------------
        y = 0._wp
        bx% hic(ic)% l% n = 2 * ke + 1  ! PS + 2 columns
        do l = bx% hic(ic)% l% i + 1, bx% hic(ic)% l% i + bx% hic(ic)% l% n
          do i = bx% vic(l)% i% i + 1, &
                 bx% vic(l)% i% i + bx% vic(l)% i% n
            select case (t_int (i))
            case (OBS_H, OBS_HS)
              il = IL_H
            case (OBS_TV)
              il = IL_TV
            case (OBS_RH)
              il = IL_RH
            case (OBS_U)
              il = IL_U
            case (OBS_V)
              il = IL_V
            case default
              call finish ('apply_LvWvKIhIv_m','invalid interpolation type')
            end select
             do k = 1, bx% vic(l)% g% n
               y(i) = y(i) +                         &
                      xl(bx% vic(l)% g% i+k-1, il) * &
                      bx% vic(l)% wh(k)
             end do  ! k
          end do     ! i
        end do       ! l

        !-------------
        ! store result
        !-------------
        if (.not.lnewpl) &
          a_m(ib)% col(ic)% s% ps =                    &
          a_m(ib)% col(ic)% s% ps + y (1) * (gacc/R) * &
         (cbg(ib)% col(ic)% s% ps / cbg(ib)% col(ic)% s% t2m)
        a_m(ib)%   col(ic)%    t           =             &
        a_m(ib)%   col(ic)%    t           + y (2:nvcf:nm)
        a_m(ib)%   col(ic)%    rh          =             &
        a_m(ib)%   col(ic)%    rh          + y (3:nvcf:nm)
        a_m(ib)%   col(ic)%    u           =             &
        a_m(ib)%   col(ic)%    u           + y (4:nvcf:nm)
        a_m(ib)%   col(ic)%    v           =             &
        a_m(ib)%   col(ic)%    v           + y (5:nvcf:nm)
!       a_m(ib)%   col(ic)%    geoh(1)     = 0._wp             !!! ICON ???
        a_m(ib)%   col(ic)%    geoh(2:ke+1)=                    &
        a_m(ib)%   col(ic)%    geoh(2:ke+1)+ y (1+nvcf:  ) * gacc

      end do           ! ic

      !-----------------------
      ! deallocate temporaries
      !-----------------------
      if (nvc /= 0) then
        nvc = 0
        deallocate (t_int)
        deallocate (e    )
        deallocate (z    )
        deallocate (y    )
        deallocate (ph   )
        deallocate (lpf  )
      endif

    end do             ! ib

  end subroutine apply_IhIv_ensb_m

!==============================================================================
end module mo_bg_err_ens
