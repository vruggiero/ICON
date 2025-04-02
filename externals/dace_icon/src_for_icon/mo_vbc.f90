!
!+ Routines for Variational Bias Correction (VarBC)
!
module mo_vbc
!
! Description:
!   This module holds routines for variational Bias Correction (VarBC).
!   It handles the extensions of operators H B Ht for VarBC.
!   Here we have:
!     H  == P   : calculate the bias correction from the coefficients
!                 by multiplying with the respective predictors
!     Ht == Pt  : adjoint, multiply with the respective predictors
!     B  == Bvbc: multiply with error covariance matrix for the coefficients
!                 as estimated by the online bias correction statistics
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_27        2013-11-08 Andreas Rhodin
!  implementation of Variational Bias Correction (VarBC)
! V1_28        2014/02/26 Andreas Rhodin
!  account for VarBC in preconditioner
! V1_48        2016-10-06 Robin Faulwetter
!  level-based thinning (changed interface to 'set_indx')
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==========================================================================
  !=============
  ! Modules used
  !=============
  use mo_kind,        only: wp, sp        ! working precision kind parameter
  use mo_mpi_dace,    only: dace,        &! MPI group info
                            p_sum,       &! sum over PEs
                            p_max         ! max over PEs
  use mo_exception,   only: finish,      &! abort in case of error
                            message       ! print warning
  use mo_radbiascor,  only: bccoef        ! bias correction statistics
  use mo_radbias_3dv, only: biascor_mode,&! bias correction mode
                            BC_VARBC,    &! flag value for VBC
                            precon        ! account for VarBC in preconditioner
  use mo_obs_set,     only: t_vbc         ! VBC data type
  use mo_fdbk_tables, only: OT_RAD        ! RADiance code type
  use mo_satid,       only: satname       ! derive satellite name from Id
  use mo_radbiascor,  only: m_pred        ! max.number of predictors used
  use mo_dec_matrix,  only: t_vector,    &! vector derived type
                            t_matrix,    &! matrix derived type
                            construct,   &! set up derived types
                            destruct,    &! destruct
                            BDIAG,       &! block diagonal flag value
                            FULL,        &! full   matrix block flag
                            PACKED,      &! packed matrix block flag
                            dot_product, &! dot product of vectors
                            assignment(=) ! assign value to vector
  use mo_t_obs,       only: t_obs,       &! observation derived type
                            t_spot        ! report/fov derived type
  use mo_obs_set,     only: t_obs_set     ! extended observation type
  use mo_rad,         only: chan_indx,   &! index of channel in t_rad
                            set_indx      ! search bccoef array index
  use mo_instrid,     only: rttov_instr   ! RTTOV number from instrument id
  implicit none

!------------------------------------------------------------------------------
  !================
  ! Public entities
  !================
  private
  public :: setup_vbc    ! set up variational bias correction
  public :: add_HBHvbc   ! add B_vbc in observation space for preconditioning
  public :: Hvbc_x       ! bias correction forward model
  public :: z_Hvbc       ! bias correction adjoint model
  public :: Bvbc_z       ! bias correction B application in coefficient space
  public :: z_HBvbcH     ! bias correction B application in observation space

!==============================================================================
contains
!==============================================================================
  subroutine setup_vbc (obs, lana)
  type (t_obs_set) ,intent(inout) :: obs  ! observation meta data
  logical          ,intent(in)    :: lana ! true for analysis scan
  !---------------------------------------------------------------------------
  ! Initialisation of Variational Bias Correction (VBC)
  !
  !   Data is hold in obs% vbc
  !   Information for the setup is taken from 'bccoef'
  !     (bias correction statistics files and derived coefficients)
  !   Tasks:
  !   - set up meta data: obs% vbc % meta
  !   - allocate arrays : obs% vbc % cf     ! bias correction coefficients
  !                       obs% vbc % z      ! gradient with respect to cf
  !                       obs% vbc % bc     ! bias correction per observation
  !                       obs% vbc % ic     ! auxiliary index field
  !                       obs% o(:)% bcpred ! predictor values for each FoV
  !--------------------------------------------------------------------------

    integer               :: i, j, k, n ! indices
    integer               :: ib         ! box index
    integer               :: is         ! spot index
    integer               :: n_meta     ! number of VarBC meta data entries
    integer               :: a_pred     ! max. number of active predictors
    integer               :: pred_use   ! number of predictors actively used
    type(t_spot) ,pointer :: s          ! report/fov pointer
    integer               :: ichan      ! channel number
    integer               :: instr      ! instrument number
    integer               :: l          ! channel index
    integer               :: n_pred     ! number of predictors calculated
    integer               :: satid      ! WMO satellite id
    integer               :: grid       ! interpolation grid instrument id

    !--------------------------------------------------
    ! check for Variational Bias Correction switched on
    !--------------------------------------------------
    if (iand (BC_VARBC, biascor_mode) == 0) return
    !---------------------------------------------------
    ! check if radiation bias correction data is present
    !---------------------------------------------------
    if (.not. associated (bccoef)) then
      call message ('setup_vbc','bccoef not associated')
      return
    endif
    !-----------------
    ! set up meta data
    !-----------------
    n_meta = size (bccoef)
    allocate (obs% vbc% meta (n_meta))
    do i = 1, n_meta
      obs% vbc % meta(i)% codetype = OT_RAD
      obs% vbc % meta(i)% platform = satname (bccoef(i)% i% satid)
!     obs% vbc % meta(i)% n_fov      = bccoef(i)% i% n_fov
      obs% vbc % meta(i)% n_fov      = 1 ! currently no fov dependence in VBC
      obs% vbc % meta(i)% n_chan     = bccoef(i)% i% n_chan
      obs% vbc % meta(i)% n_coeff    = bccoef(i)% n_coeff
      obs% vbc % meta(i)% n_pred     = bccoef(i)% n_pred
      obs% vbc % meta(i)% n_used     = bccoef(i)% mp
      obs% vbc % meta(i)% is         = i
      obs% vbc % meta(i)% o          = 0
      !-----------------------------------
      ! determine predictors actually used
      !-----------------------------------
      obs% vbc % meta(i)% n_tot      = (obs% vbc % meta(i)% n_fov +    &
                                        obs% vbc % meta(i)% n_used ) * &
                                        obs% vbc % meta(i)% n_chan
    end do

    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write (6,'()')
      write (6,'(a)')'  VBC Setup :'
      write (6,'()')
      write (6,'(a)')                                                 &
      '        code platform    fov   coef   pred   used   chan    tot'
      write (6,'()')
      do i = 1, n_meta
        write (6,'(4x,2i4,1x,a8,7i7)') i,                             &
                                       obs% vbc % meta(i)% codetype,  &
                                       obs% vbc % meta(i)% platform,  &
                                       obs% vbc % meta(i)% n_fov,     &
                                       obs% vbc % meta(i)% n_coeff,   &
                                       obs% vbc % meta(i)% n_pred,    &
                                       obs% vbc % meta(i)% n_used,    &
                                       obs% vbc % meta(i)% n_chan,    &
                                       obs% vbc % meta(i)% n_tot
      end do
      write (6,'()')
    end if

    !----------------
    ! allocate arrays
    !----------------
    if (lana) then
      allocate       (obs% vbc% bi)
      call construct (obs% vbc% bi, &
                      sum(obs% vbc% meta% n_tot), & ! number of elements
                      n_meta,                     & ! number of blocks
                      obs% vbc% meta% n_tot,      & ! # of elements / block
                      (/(dace% pio,i=1,n_meta)/)  ) ! processor p_io so far

      call construct (obs% vbc% z,  obs% vbc% bi, 'z_vbc', global=.true.)
      call construct (obs% vbc% cf, obs% vbc% bi, 'd_vbc', global=.true.)
      call construct (obs% vbc% ic, obs% oi,      'ic'                  )
      call construct (obs% vbc% bc, obs% oi,      'bc'                  )
      obs% vbc% z  = 0._wp
      obs% vbc% cf = 0._wp
      obs% vbc% bc = 0._wp
      obs% vbc% ic = 0
    endif
    !-----------------------------
    ! set up profile specific data
    !-----------------------------
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe)                        cycle
      n = 0
      do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_RAD) cycle
        s => obs% o(ib)% spot(is)
        !----------------------------------------
        ! set auxilliary array with channel index
        !----------------------------------------
        satid = s% hd% satid
        grid  = s% hd% grid_id
        i     = set_indx (bccoef(:)% i, satid=satid, grid=grid)
        if (lana) then
          do k = s%o%i+1, s%o%i + s%o%n
            ichan =              nint(obs% o(ib)% olev (k))
            instr = rttov_instr (int (obs% o(ib)% body (k)% lev_sig), satid)
            l = chan_indx (instr, ichan, bccoef(i)% i)
            if (l <= 0) then
              write(0,*) ' determine channel index: set, instr, chan, index =',&
                                                      i, instr, ichan, l
              call finish('setup_vbc','Failed to get channel index')
            end if
            obs% vbc% ic% s(ib)% x(k) = l
          end do
        endif
        !------------------------------------------------------------
        ! set up references to predictor values and coefficient table
        !------------------------------------------------------------
        n_pred = bccoef(i)% n_pred
        s% bcp% i = n
        s% bcp% n = n_pred
        s% ibcb   = i
        s% ibcfov = n_pred + s% phase  ! currently not used
        n = n     + n_pred
      end do
      !-----------------------------------------------
      ! allocate array components for predictor values
      !-----------------------------------------------
      if (obs% o(ib)% n_bcp == 0) then
        obs% o(ib)% n_bcp = n
        allocate (obs% o(ib)% bcpred (n))
      endif
    end do

  end subroutine setup_vbc

!==============================================================================

  subroutine Hvbc_x (bc, cf, obs)
  type(t_vector)  ,intent(inout) :: bc   ! bias correction
  type(t_vector)  ,intent(in)    :: cf   ! bias correction coefficients
  type(t_obs_set) ,intent(in)    :: obs  ! observational data
  !----------------------------------------------------------
  ! derive bias correction 'bc' for all radiance observations
  ! from bias correction coefficients 'cf'
  !----------------------------------------------------------

    integer  :: ib     ! box index
    integer  :: is     ! spot index

    !--------------------------------
    ! return if VarBC is switched off
    ! pre-set bias correction to zero
    !--------------------------------
    if (iand (BC_VARBC, biascor_mode) == 0) return
    bc = 0._wp
    !-----------------------------------
    ! loop over 'boxes' and FoVs (spots)
    !-----------------------------------
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_RAD) cycle
        !------------------------------------------
        ! derive bias correction for individual FoV
        !------------------------------------------
        call Hvbc_x_fov (bc, cf, obs, ib, is)
      end do
    end do

  end subroutine Hvbc_x

!------------------------------------------------------------------------------

  subroutine Hvbc_x_fov (bc, cf, obs, ib, is)
  type(t_vector)  ,intent(inout) :: bc   ! bias correction
  type(t_vector)  ,intent(in)    :: cf   ! bias correction coefficients
  type(t_obs_set) ,intent(in)    :: obs  ! observational data
  integer         ,intent(in)    :: ib   ! box index
  integer         ,intent(in)    :: is   ! spot index
  !------------------------------------------------------------
  ! evaluate bias correction for given FoV (box and spot index)
  !------------------------------------------------------------

    integer  :: ibcb   ! satellite/sensor index
    integer  :: ibcfov ! field of view index    ! currently not used
    integer  :: ip0    ! predictor start index in t_obs
    integer  :: istep  ! coefficients / channel, currently #of predictors+1
    integer  :: k, i   ! observation (channel) indices in t_spot
    integer  :: ic     ! channel index                 in t_bccoef
    integer  :: j      ! predictor index (used)        in t_bccoef
    integer  :: l      ! predictor index (all)         in t_spot

    !------------
    ! set indices
    !------------
    ibcb   =  obs% o(ib)% spot(is)% ibcb
    ibcfov =  obs% o(ib)% spot(is)% ibcfov   ! currently not used
    ip0    =  obs% o(ib)% spot(is)% bcp% i
    istep  =  obs% vbc%  meta(ibcb)% n_fov + &
              obs% vbc%  meta(ibcb)% n_used
    !--------------------------
    ! loop over channels in FoV
    !--------------------------
    do k = 1, obs% o(ib)% spot(is)% o% n
      i = k + obs% o(ib)% spot(is)% o% i
      ic    = obs% vbc% ic% s(ib)% x(i)
      bc% s(ib)% x(i) = 0._wp            ! preset correction to zero
      !--------------------------------------
      ! add contributions of model predictors
      !--------------------------------------
      do j = 1, bccoef(ibcb)% n_b (ic,1)
        l = bccoef(ibcb)% i_p (j,ic)
          bc% s(ib)% x(i) = bc% s(ib)% x(i) +            &
                            obs% o(ib)% bcpred (ip0+l) * &
                            cf% s(ibcb)% x((ic-1)*istep+j)
      end do
      !----------------------------------- ++++++++++++++++++++
      ! add contributions of constant term (currently only one)
      !----------------------------------- ++++++++++++++++++++
      j = obs% vbc%  meta(ibcb)% n_used + 1
      bc% s(ib)% x(i) = bc% s(ib)% x(i) + cf% s(ibcb)% x((ic-1)*istep+j)
    end do

  end subroutine Hvbc_x_fov

!==============================================================================

  subroutine z_Hvbc (z, zc, obs)
  type(t_vector)  ,intent(in)    :: z   ! gradient in observation space
  type(t_vector)  ,intent(inout) :: zc  ! gradient with respect to coeff.
  type(t_obs_set) ,intent(in)    :: obs ! observational meta data
  !------------------------------------------------------------------
  ! adjoint bias correction evaluation (in VBC context):
  ! derive gradient with respect to bias correction coefficients 'zc'
  !   from gradient with respect to observations                 'z'.
  !------------------------------------------------------------------

    integer  :: ib     ! box index
    integer  :: is     ! spot index
    integer  :: ibcb   ! satellite/sensor index

    !------------------------------------------
    ! immediate return if VarBC is switched off
    !------------------------------------------
    if (iand (BC_VARBC, biascor_mode) == 0) return
    !-----------------------------------------------------------
    ! ensure that observation error is stored for VarBC B matrix
    !-----------------------------------------------------------
    call fill_varbc_obserr (obs)
    !---------------------
    ! zero result variable
    !---------------------
    do ibcb = 1, size (zc% s)
      zc% s(ibcb)% x = 0._wp
    end do
    !-------------------------
    ! loop over field of views
    !-------------------------
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_RAD) cycle
        !---------------------------------------
        ! call adjoint H operator for single FoV
        !---------------------------------------
        call z_Hvbc_fov (z, zc, obs, ib, is)
      end do
    end do
    !---------------------------------
    ! sum final result over processors
    !---------------------------------
    do ibcb = 1, size (zc% s)
      zc% s(ibcb)% x = p_sum (zc% s(ibcb)% x)
    end do

  end subroutine z_Hvbc

!------------------------------------------------------------------------------

  subroutine z_Hvbc_fov (z, zc, obs, ib, is)
  type(t_vector)  ,intent(in)    :: z   ! gradient in observation space
  type(t_vector)  ,intent(inout) :: zc  ! gradient with respect to coeff.
  type(t_obs_set) ,intent(in)    :: obs ! observational meta data
  integer         ,intent(in)    :: ib  ! box index
  integer         ,intent(in)    :: is  ! spot index
  !-----------------------------------------------------------------
  ! adjoint bias correction evaluation for single FoV (index ib, is)
  ! derive gradient with respect to bias correction coefficients 'zc'
  !   from gradient with respect to observations                 'z'.
  !-----------------------------------------------------------------

    integer  :: ibcb   ! satellite/sensor index
    integer  :: ibcfov ! field of view index   ! currently not used
    integer  :: ip0    ! predictor start index in t_obs
    integer  :: istep  !
    integer  :: k, i   ! observation (channel) indices
    integer  :: j      ! predictor index
    integer  :: ic     ! channel index
    integer  :: l      ! predictor index
    real(wp) :: d      ! gradient in observation space

    !-----------------------------------------
    ! corresponding indices to BC coefficients
    !-----------------------------------------
    ibcb   =  obs% o(ib)% spot(is)% ibcb
    ibcfov =  obs% o(ib)% spot(is)% ibcfov   ! currently not used
    ip0    =  obs% o(ib)% spot(is)% bcp% i
    istep  =  obs% vbc%  meta(ibcb)% n_fov + &
              obs% vbc%  meta(ibcb)% n_used
    !--------------------------
    ! loop over channels in FoV
    !--------------------------
    do k = 1, obs% o(ib)% spot(is)% o% n
      i = k + obs% o(ib)% spot(is)% o% i
      ic    = obs% vbc% ic% s(ib)% x(i)
      d     = z% s(ib)% x(i)
      !--------------------------------------
      ! add contributions to model predictors
      !--------------------------------------
      do j = 1, bccoef(ibcb)% n_b (ic,1)
        l = bccoef(ibcb)% i_p (j,ic)
        zc%   s(ibcb)% x((ic-1)*istep+j) = &
          zc% s(ibcb)% x((ic-1)*istep+j) + &
          d * obs% o(ib)% bcpred (ip0+l)
      end do
      !----------------------------------- ++++++++++++++++++++
      ! add contributions to constant term (currently only one)
      !----------------------------------- ++++++++++++++++++++
      j = obs% vbc%  meta(ibcb)% n_used + 1
      zc%   s(ibcb)% x((ic-1)*istep+j) = &
        zc% s(ibcb)% x((ic-1)*istep+j) + d
    end do

  end subroutine z_Hvbc_fov
!==============================================================================

  subroutine add_HBHvbc (obs, HBH)
  type (t_obs_set) ,intent(in)    :: obs  ! observation meta data
  type (t_matrix)  ,intent(inout) :: HBH  ! HBH (B in observation space)
  !-------------------------------------------------------------------------
  ! add B_vbc in observation space for preconditioning
  ! i.e. calculates H B Ht for all pairs of radiance observations in a 'box'
  !-------------------------------------------------------------------------

    !--------
    ! indices
    !--------
    integer               :: ib         ! box index
    integer               :: isr        ! spot index rhs
    integer               :: isl        ! spot index lhs
    integer               :: ibcb       ! satellite/sensor index
    integer               :: mp         ! max. number of predictors used
    integer               :: repr       ! matrix representation (FULL, PACKED)
    character             :: tri        ! 'U'pper or 'L'ower (for PACKED repr.)
    integer               :: i1, in     ! lhs observation index bounds
    integer               :: j1, jn     ! rhs observation index bounds
    integer               :: i, j       ! observation indices
    integer               :: ic, jc     ! corresponding channel indices
    integer               :: k, n, kn   ! bookkeeping for PACKED representation
    integer               :: isl1, isln ! spot index bounds, PACKED repr.
    !--------------------
    ! temporary variables
    !--------------------
    type(t_vector)        :: d       ! adjoint variable (set to 1)
    type(t_vector)        :: bc      ! resulting covariance
    type(t_vector)        :: cf      ! intermediate: bc coefficient
    type(t_vector)        :: z       ! intermediate: gradient cf.bc
    real(wp) ,allocatable :: B (:,:) ! temporary: B_vbc,            one channel
    real(wp) ,allocatable :: cf_ (:) !            bias coefficient, one channel
    real(wp) ,allocatable ::  z_ (:) !            gradient with respect to b.c.

    !------------------------------------------
    ! immediate return if VarBC is switched off
    !------------------------------------------
    if (iand (BC_VARBC, biascor_mode) == 0) return
    if (.not. precon)                       return

    !-----------------------------------------------------------
    ! ensure that observation error is stored for VarBC B matrix
    !-----------------------------------------------------------
    call fill_varbc_obserr (obs)

    !-----------------------------
    ! allocate temporary variables
    !-----------------------------
    call construct (cf, obs% vbc% bi, global=.true.)
    call construct (z,  obs% vbc% bi, global=.true.)
    call construct (d,  obs%      oi)
    call construct (bc, obs%      oi)
    d  = 1._wp
    bc = 0._wp
    mp = maxval (bccoef(:)% mp)
    allocate (B (mp,mp), cf_(mp+1), z_(mp+1))

    !----------------
    ! loop over boxes
    !----------------
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      !----------------------------------------------
      ! check for proper matrix representation of HBH
      !----------------------------------------------
      repr = HBH% b(ib,ib)% repr ! FULL or PACKED matrix reoresentation
      tri  = HBH% b(ib,ib)% tri !  'U'pper or 'L'ower triangle stored
      n    = HBH% b(ib,ib)% n   !  number of observations in this box
      kn   = (n*(n+1))/2        !  size of packed matrix storage
      if (n == 0) cycle
      select case (repr)
      case (PACKED)
      case (FULL)
        tri = 'L'               ! arbitrary choice
      case default
        call finish ('add_HBHvbc','HBH% b /= PACKED or FULL')
      end select
      !-------------------
      ! loop over rhs FOVs
      !-------------------
      do isr = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(isr)% hd% obstype /= OT_RAD) cycle
        ibcb   =  obs% o(ib)% spot(isr)% ibcb
        z% s(ibcb)% x = 0._wp
        call z_Hvbc_fov   (d, z, obs, ib, isr)
        call Bvbc_z_instr (cf, z, obs, ibcb, B, cf_, z_)
        j1 = obs% o(ib)% spot (isr)% o% i + 1
        jn = obs% o(ib)% spot (isr)% o% i + obs% o(ib)% spot (isr)% o% n
        !-------------------------------
        ! loop over lhs FOVs
        ! (packed matrix representation)
        !-------------------------------
        select case (tri)
          case ('U')
          isl1 = 1
          isln = isr
        case ('L')
          isl1 = isr
          isln = obs% o(ib)% n_spot
        end select
        do isl = isl1, isln
          !-------------------------------------
          ! determine elements for given rhs fov
          !-------------------------------------
          if (obs% o(ib)% spot(isl)% hd% obstype /= OT_RAD) cycle
          if (obs% o(ib)% spot(isl)%     ibcb    /= ibcb  ) cycle
          call Hvbc_x_fov (bc, cf, obs, ib, isl)
          i1 = obs% o(ib)% spot (isl)% o% i + 1
          in = obs% o(ib)% spot (isl)% o% i + obs% o(ib)% spot (isl)% o% n
          !----------------------------
          ! seek corresponding channels
          !----------------------------
          i  = i1
          j  = j1
          ic = obs% vbc% ic% s(ib)% x(i)
          jc = obs% vbc% ic% s(ib)% x(j)
          do
            if (ic < jc) then
              i = i + 1
              if (i > in) exit
              ic = obs% vbc% ic% s(ib)% x(i)
              cycle
            endif
            if (jc < ic) then
              j = j + 1
              if (j > jn) exit
              jc = obs% vbc% ic% s(ib)% x(j)
              cycle
            endif
            !----------------------
            ! add to matrix element
            !----------------------
            select case (HBH% b(ib,ib)% repr)
            case (PACKED)
              select case (tri)
              case ('U')
                k = (j*(j-1))/2+i
              case ('L')
                k = kn - ((n-j+1)*(n-j+2))/2 + i-j+1
              end select
                HBH% b(ib,ib)% packed (k) = HBH% b(ib,ib)% packed (k) &
                                          + bc % s(ib)   % x    (i  )
            case (FULL)
                HBH% b(ib,ib)% full (i,j) = HBH% b(ib,ib)% full (i,j) &
                                          + bc % s(ib)   % x    (i  )
              if (i /= j) then
                HBH% b(ib,ib)% full (j,i) = HBH% b(ib,ib)% full (j,i) &
                                          + bc % s(ib)   % x    (i  )
              endif
            end select
            i = i + 1
            if (i > in) exit
            ic = obs% vbc% ic% s(ib)% x(i)
            j = j + 1
            if (j > jn) exit
            jc = obs% vbc% ic% s(ib)% x(j)
          end do
        end do
      end do
    end do

    !-------------------------------
    ! deallocate temporary variables
    !-------------------------------
    call destruct (cf)
    call destruct (z )
    call destruct (d )
    call destruct (bc)

  end subroutine add_HBHvbc

!==============================================================================

  subroutine Bvbc_z (cf, z, obs)
  type(t_vector)  ,intent(inout) :: cf  ! coefficients
  type(t_vector)  ,intent(in)    :: z   ! gradient with respect to coeff.
  type(t_obs_set) ,intent(inout) :: obs ! observational data
  !---------------------------------------------------------------------
  ! multiply a vector with the VarBC background error correlation matrix
  ! input: gradient with respect to the BC coefficients 'z'
  ! output:                             BC coefficients 'cf'
  !---------------------------------------------------------------------

    integer               :: ibcb    ! satellite/sensor index
    integer               :: mp      ! max # of active predictors
    real(wp) ,allocatable :: B (:,:) ! temporary: B  for 1 channel
    real(wp) ,allocatable :: cf_ (:) !            cf for 1 channel
    real(wp) ,allocatable ::  z_ (:) !            z  for 1 channel

    !--------------------------------
    ! return if VarBC is switched off
    !--------------------------------
    if (iand (BC_VARBC, biascor_mode) == 0) return

    !---------------------
    ! allocate temporaries
    !---------------------
    mp = maxval (bccoef(:)% mp)
    allocate (B (mp,mp), cf_(mp+1), z_(mp+1))
    !-------------------------------
    ! loop over satellites / sensors
    !-------------------------------
    do ibcb = 1, size (z% s)
      call Bvbc_z_instr (cf, z, obs, ibcb, B, cf_, z_)
    end do
    !---------------------------
    ! update VarBC cost function
    !---------------------------
    obs% vbc% J = 0.5_wp * dot_product (obs% vbc% cf, obs% vbc% z)

  end subroutine Bvbc_z

!------------------------------------------------------------------------------

  subroutine Bvbc_z_instr (cf, z, obs, ibcb, B, cf_, z_)
  type(t_vector)  ,intent(inout) :: cf      ! coefficients
  type(t_vector)  ,intent(in)    :: z       ! gradient with respect to coeff.
  type(t_obs_set) ,intent(in)    :: obs     ! observational data
  integer         ,intent(in)    :: ibcb    ! satellite/sensor index
  real(wp)        ,intent(inout) :: B (:,:) ! temporary: B  for 1 channel
  real(wp)        ,intent(inout) :: cf_ (:) !            cf for 1 channel
  real(wp)        ,intent(inout) ::  z_ (:) !            z  for 1 channel
  !---------------------------------------------------------------------
  ! multiply a vector with the VarBC background error correlation matrix
  ! application for a given set of instruments (index ibcb)
  ! input: gradient with respect to the BC coefficients 'z'
  ! output:                             BC coefficients 'cf'
  !---------------------------------------------------------------------

    integer  :: ic     ! channel index
    integer  :: istep  ! entries per channel
    integer  :: ic0    ! channel start index in vector
    integer  :: ib0    ! channel start index in B
    integer  :: nb
    integer  :: nbb
    integer  :: mp     ! max # of active predictors
    integer  :: i      ! active  predictor index

    !----------------------------------
    ! prozess for 1 satellites / sensor
    !----------------------------------
    cf% s(ibcb)% x = 0._wp
    if (bccoef(ibcb)% entries == 0._wp) return
    istep = obs% vbc%  meta(ibcb)% n_fov + &
            obs% vbc%  meta(ibcb)% n_used
    !-------------------
    ! loop over channels
    !-------------------
    do ic = 1, obs% vbc%  meta(ibcb)% n_chan
      ic0 = (ic-1)*istep
      nb  = bccoef(ibcb)% n_b(ic,1)
      nbb = bccoef(ibcb)% n_b(ic,2)
      ib0 = bccoef(ibcb)% n_b(ic,3)
      if (bccoef(ibcb)% n(ic) > 0._wp) then  ! check for BC statistics
        !------------------
        ! model predictor B
        !------------------
        b (:nb,:nb) = reshape (bccoef(ibcb)% p_b(ib0+1:ib0+nbb), (/nb,nb/)) &
                    * bccoef(ibcb)% w_vbc
        !------------------
        ! constant term rhs
        !------------------
        z_(nb+1) = z% s(ibcb)% x (ic0 + obs% vbc%  meta(ibcb)% n_used + 1)
        !-------------------------------------
        ! remove mean from model predictor rhs
        !-------------------------------------
        do i = 1, nb
          z_(i) = (z% s(ibcb)% x (ic0 + i)                 &
                - bccoef(ibcb)% p_mean  (i,ic) * z_(nb+1)) &
                * bccoef(ibcb)% p_scal  (i,ic)
        end do
        !------------------------
        ! B matrix multiplication
        !------------------------
        cf_(:nb)   = matmul (b (:nb,:nb), z_(:nb)) &
                   * bccoef(ibcb)% p_scal   (:nb,ic)
        cf_(nb+1) = z_(nb+1) * bccoef(ibcb)% w_vbc / bccoef(ibcb)% n(ic)
        !---------------------------------------------------------------
        ! scale with observation error to comply with offline B estimate
        !---------------------------------------------------------------
        cf_(:nb+1) = cf_(:nb+1) * bccoef(ibcb)% i% var (ic)
        !-----------------------------------
        ! account for predictor mean removed
        !-----------------------------------
        cf% s(ibcb)% x (ic0 + obs% vbc%  meta(ibcb)% n_used + 1) = cf_(nb+1)
        do i = 1, nb
          cf%   s(ibcb)% x (ic0 + i                                ) = cf_(i)
          cf%   s(ibcb)% x (ic0 + obs% vbc%  meta(ibcb)% n_used + 1) = &
            cf% s(ibcb)% x (ic0 + obs% vbc%  meta(ibcb)% n_used + 1) - &
            bccoef(ibcb)% p_mean  (i,ic) * cf_(i)
        end do
      endif
    end do

  end subroutine Bvbc_z_instr

!==============================================================================

  subroutine z_HBvbcH  (obs, z, bc)
  type(t_obs_set) ,intent(inout)           :: obs ! observational meta data
  type(t_vector)  ,intent(in)              :: z   ! gradient in obsv. space
  type(t_vector)  ,intent(inout) ,optional :: bc  ! bias correction
  !---------------------------------------------------
  ! application of VarBC B matrix in observation space
  ! ( apply H B H^t)
  !---------------------------------------------------

    if (iand (BC_VARBC, biascor_mode) /= 0) then
      call z_Hvbc (z, obs% vbc% z, obs)             ! H^t
      call Bvbc_z (obs% vbc% cf, obs% vbc% z,  obs) ! B
      call Hvbc_x (obs% vbc% bc, obs% vbc% cf, obs) ! H
      if (present(bc)) bc =       obs% vbc% bc
    else
      if (present(bc)) bc = 0._wp
    endif

  end subroutine z_HBvbcH

!==============================================================================

  subroutine fill_varbc_obserr (obs)
  type(t_obs_set) ,intent(in)  :: obs ! observational meta data
  !-----------------------------------------------------------------------
  ! ensure that observation error is stored for VarBC B matrix application
  !-----------------------------------------------------------------------

    integer :: ibcb
    integer :: ic
    logical :: missing
    integer :: ib
    integer :: is
    integer :: i, k

    !-----------------------------
    ! check if R is stored already
    !-----------------------------
    missing = .false.
    do ibcb = 1, size (bccoef)
      if (.not.associated(bccoef(ibcb)% i% var)) then
        missing = .true.
        allocate (bccoef(ibcb)% i% var (bccoef(ibcb)% i% n_chan))
        bccoef(ibcb)% i% var = 0._wp
      endif
    end do

    if (missing) then
      !-----------------------------------------
      ! copy R from observational meta data
      ! to bias correction coefficient meta data
      !-----------------------------------------
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
        do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype /= OT_RAD) cycle
          ibcb   =  obs% o(ib)% spot(is)% ibcb
          do k = 1, obs% o(ib)% spot(is)% o% n
            i = k + obs% o(ib)% spot(is)% o% i
            ic    = obs% vbc% ic% s(ib)% x(i)
            if (bccoef(ibcb)% i% var(ic) == 0._wp) then
              bccoef(ibcb)% i% var(ic) = obs% o(ib)% body(i)% eo ** 2
            endif
          end do
        end do
      end do
      !-----------------------
      ! communicate to all PEs
      !-----------------------
      do ibcb = 1, size (bccoef)
        bccoef(ibcb)% i% var =  p_max (bccoef(ibcb)% i% var)
      end do

    endif

  end subroutine fill_varbc_obserr

!==============================================================================

end module mo_vbc
