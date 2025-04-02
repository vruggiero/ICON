!
!+ set up background and observation error covariance matrices (PSAS)
!
MODULE mo_set_matrix
!
! Description:
!   This module holds routines to set and manipulate (background and
!   observational) error covariance matrices. For this purpose the module
!   has to know about (include the respective modules) the data structures
!   sparse matrices are represented by (module mo_dec_matrix), the data
!   structures holding information on the observation operators (module
!   mo_t_obs), the background error covariance model (mo_fg_cov) and the
!   variational quality control, i.e. the definition of non-Gaussian
!   observational errors (module mo_vqc), and the representation of
!   atmospheric fields organized in model columns (module mo_t_col).
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  fix for zero number of observations in a box
! V1_5         2009/05/25 Andreas Rhodin
!  NEC SX (USE_GET_CORR_OPT): restrict size of array rv_coi to allocate
! V1_6         2009/06/10 Harald Anlauf
!  set_Pb(get_corr_opt): fix bugs ( SX + GPSRO ); optimize revised code
! V1_9         2010/04/20 Andreas Rhodin
!  restrict max.number of attempts to increase diagonal of precond.matrix
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  technical changes
! V1_14        2011/11/08 Harald Anlauf
!  add diagnostics for vectorization
! V1_20        2012-06-18 Andreas Rhodin
!  changed comments and names of public routines in module mo_p_output
! V1_22        2013-02-13 Andreas Rhodin
!  skip setup of interpolation operators for COSMO/GPSGB operators.
! V1_23        2013-03-26 Andreas Rhodin
!  keep B in interpolation space for diagnostic printout in cofRTOVP.nc
! V1_26        2013/06/27 Andreas Rhodin
!  replace STAT_ACTIVE_0 by STAT_NOTACTIVE
! V1_27        2013-11-08 Andreas Rhodin
!  bugfix for option to keep radiance B matrix in interpolation space
! V1_28        2014/02/26 Harald Anlauf
!  Workaround for CRAYFTN optimization bug
! V1_31        2014-08-21 Andreas Rhodin
!  pass IR PC gradient (implement adjoint mode)
!  preparations for rttov specific vertical interpolation
! V1_35        2014-11-07 Andreas Rhodin
!  account for Jakobian of sink variable transformation in preconditioning
! V1_37        2014-12-23 Andreas Rhodin
!  include VarEnKF B in preconditioner
! V1_42        2015-06-08 Andreas Rhodin
!  EnVar/PSAS optimisation; option to use inflated R from 3dvar VQC in LETKF
! V1_45        2015-12-15 Harald Anlauf
!  Fix kind parameters; small optimization
! V1_47        2016-06-06 Andreas Rhodin
!  correct extraction of B_ii matrix (for 1dvar & diagnostics)
! V1_48        2016-10-06 Andreas Rhodin
!  minor cleanup
! V1_51        2017-02-24 Andreas Rhodin
!  option to account for nonlinearity of generalised humidity in preconditioner
!  namelist /psas/ fg_1col,an_1col: option to approcimate NMC-matrices
!  option to construct ensemble B preconditioner in observation space
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2004-2008  original source
! Harald Anlauf   DWD  2007       changes for SX8
! Jan Boerhout    NEC  2008       optimised 'get_corr_opt' for SX8
!==============================================================================
! Define USE_GET_CORR_OPT to prefer the vectorized version of get_corr
!#define USE_GET_CORR_OPT
!
! Define DEBUG_GET_CORR_OPT to check results from get_corr_opt with get_corr
!#define DEBUG_GET_CORR_OPT
!
!#if defined (__NEC__)
!#define USE_GET_CORR_OPT
!#endif
!==============================================================================
!
! Diagnostics for vectorization
!
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!==============================================================================
  !=============
  ! Modules used
  !=============
  !------------------------
  ! general purpose modules
  !------------------------
  use mo_kind,       only: wp, i8          ! kind parameter
  use mo_exception,  only: finish          ! abort routine
  use mo_mpi_dace,   only: dace,          &! MPI group info
                           p_send,        &! generic MPI send      routine
                           p_recv,        &! generic MPI receive   routine
                           p_bcast,       &! generic MPI broadcast routine
                           p_barrier,     &!         MPI barrier   routine
                           p_and,         &! generic MPI and       routine
                           p_sum,         &! generic MPI sum       routine
                           p_min           ! generic MPI min       routine
  use mo_time,       only: t_time          ! time+date data type
  use mo_p_output,   only: oline,         &! output line buffer
                           iol,           &! number of next line to write
                           nextline,      &! routine to increment line number
                           flush_buf,     &! routine to write buffer
                           add_line_pio    ! routine to write string on I/O PE
  use mo_fortran_units,                    &
                     only: get_unit_number,&
                           return_unit_number
  use mo_algorithms, only: jacobi,        &! calculate eigenvalues and vectors
                           index           ! sort index of an array
  use mo_matrix,     only: check_rs,      &! check eigenvalues of RS matrix
                           operator(.o.)   ! outer product
  use mo_physics,    only: rearth          ! earth radius (m)
  !----------------------------------
  ! representation of sparse matrices
  !----------------------------------
  use mo_dec_matrix, only: t_matrix,      &! covariance matrix data type
                           t_matrix_block,&!   component of t_matrix
                           t_dec_info,    &!
                           t_vector,      &! decomposed vector
                           mp,            &! matrix coefficients kind parameter
                           construct,     &! construct vectors and matrices
                           destruct,      &! destruct vectors and matrices
                           deallocate,    &! deallocate matrix block components
                           allocate_block,&! allocate matrix block components
                           operator  (*), &! matrix-vector multiply
                           operator  (-), &! difference of vectors
                           assignment(=), &! assign vectors and matrices
                           sum,           &! sum over vector elements
                           inverse,       &! inverse of matrix block
                           outer_product, &! matrix * vector
                           add_to,        &! add matrices
                           cholesky,      &! perform Cholesky decomposition
                           pack_matrix,   &! pack matrix (store only lower tr.)
                           full_matrix,   &! convert block to full represent.
                           sub_block,     &! extract submatrix
                           insert,        &! insert  submatrix
                           csr_matrix,    &! store compressed sparse row format
                           csc_matrix,    &! store compressed sparse col format
                           p_send,        &! MPI send routine
                           p_recv,        &! MPI receive routine
                           release_mem,   &! release unused memory
                           ZERO,          &! matrix representation flags
                           DIAGONAL,      &!          ''
                           MIRROR,FULL,   &!          ''
                           NOT_INIT,      &!          ''
                           CSR, CSC,      &!
                           FULL, PACKED,  &!
                           BDIAG,         &!
                           crep,          &! convert representation to char
                           srep            ! convert representation to mnemonic
  !----------------------------
  ! description of observations
  !----------------------------
  use mo_t_obs,      only: t_obs,         &! observation data type
                           t_spot,        &!   component of t_obs
                           t_coord,       &! coordinate datatype
                           rct_name,      &! get data code type name
                           TEMP,          &! TEMP/PILOT module Id
                           TOVS,          &! Radiances  module Id
                           CHR_ID,        &!
                           CHR_EVAL,      &!
                           CHR_REQ,       &!
                           CHR_NONL,      &!
                           CHR_EXP,       &!
                           CHR_COM,       &!
                           OBS_RH,        &!
                           OBS_DUM,       &!
                           OBS_CTR,       &!
                           ITY_ICOL,      &! interpolation type: column
                           ITY_ICOLS,     &!                     columns
                           ITY_MCOLS       !               model columns
  use mo_fdbk_tables,only: VN_U,          &!          wind component code
                           VN_V,          &!                         code
                           VN_Z,          &!     geopotential height code
                           OT_RAD          ! RADIANCE observation type
  use mo_obs_set,    only: t_obs_set,     &! obs. data type
                           t_obs_block,   &!
                           obs_block       !
  use mo_t_use,      only: STAT_NOTACTIVE,&! flag for inactive data
!                          STAT_ACTIVE_1, &! flag for always active data
                           STAT_DISMISS    ! flag for dismissed data
  use mo_sink,       only: bc_bd,         &! correct B (nonlinear sink var.)
                           bc_bv           ! correct B (generalised humidity)
  !----------------------------------
  ! background error covariance model
  !----------------------------------
  use mo_t_bg_err_op,only: compress_cov,  &! store covm only once
                           uncompress_cov  ! store covm on each PE
  use mo_bg_err_ens, only: apply_Hi_ens    ! interpolation (EnKF -> obs. input)
  use mo_fg_cov,     only: t_rowcol,      &! type for precalculated quantities
                           construct,     &! construct variable oftype t_rowcol
                           destruct,      &! destruct variable of type t_rowcol
                           init_fg_cov,   &! initialise dwd correlation model
                           init_spot_obs, &! set from obs.data type
                           get_corr,      &! get correlations for pair of spots
#ifdef USE_GET_CORR_OPT
                           get_corr_opt,  &! same, optimised version (NECJB)
#endif
                           L_max,         &! max. hor.length scale(unit sphere)
                           x_cut,         &! cutoff length scale
                           test_bl         ! flag to test blockdiagonal of HBHt
  use mo_cntrlvar,   only: gh_meta,       &! generalised humidity parameters
                           nl_prec         ! account for nonlinearity in precond.
  !---------------------------------
  ! atmospheric state and ensemble B
  !---------------------------------
  use mo_atm_state,  only: t_atm           ! atmospheric state derived type
  use mo_varenkf,    only: Benkf,         &! variational ensemble B matrix
!                          w_nmc_b,       &! weight of NMC  B matrix
                           w_ens_b         ! weight of EnKF B matrix
  use mo_letkf_util, only: gaspari_cohn    ! Gaspari&Cohn function
  !---------------------------------------------------
  ! observational errors (Variational Quality Control)
  !---------------------------------------------------
  use mo_vqc,        only: vqc_uncor,     &! uncorrelated observations
                           vqc_wind,      &! correlated wind components
                           vqc_corr,      &! general correlated observations
                           use_wind,      &! flag to use vqc_wind  routine
                           use_corr,      &! flag to use vqc_uncor routine
                           t_vqc,         &!
                           g_rej,         &! calculate probab. from stdev.
                           mvqc, svqc, gvqc!

#ifdef USE_GET_CORR_OPT
  use mo_fg_cov,     only: t_admin_gc,    &! get_corr_opt vector admin
                           tv_gc,         &!
                           nv_gc
#endif

  implicit none

  !================
  ! Public entities
  !================
  private
  public:: set_flags     ! set up flags in observation variable
  public:: set_Pb        ! set up the background  error covariance matrix
  public:: set_R         ! set up the observation error covariance matrix
  public:: set_H         ! set up the linearised observation operator
  public:: Pb_times_z    ! Calculate the vector matrix product Pb * z
  public:: scatter_K     ! scatter H,R over PEs
  public:: test_c        ! test routine for positive definiteness, etc.
  public:: constrain_ev  ! constrain eigenvalues
  public:: mt            ! number of tests in constrain_ev
  public:: PB_NONE,  PB_ZERO, PB_DUMMY,                  &
           PB_IDENT, PB_DIAG, PB_1DVAR, PB_BLOCKD, PB_FULL
  !=================
  ! Module variables
  !=================
  !------------------------------------
  ! valid values for parameter pb_apprx
  !------------------------------------
  integer ,parameter :: PB_NONE   = -1  ! no calculation of Pb
  integer ,parameter :: PB_ZERO   =  0  ! zero matrix
  integer ,parameter :: PB_DUMMY  =  1  ! dummy sink variables only
  integer ,parameter :: PB_IDENT  =  2  ! HPbH = Identity
  integer ,parameter :: PB_DIAG   =  3  !  Pb  = Diagonal (variance)
  integer ,parameter :: PB_1DVAR  =  4  ! no correlations between reports
  integer ,parameter :: PB_BLOCKD =  5  ! blockdiagonal matrix
  integer ,parameter :: PB_FULL   =  6  ! full matrix
  !------------
  ! Debug flags
  !------------
  logical  :: lpack    = .true.  ! packed representation for diagonal blocks
  logical  :: test_r   = .false. ! test observations for positive definiteness
  logical  :: lmirror  = .true.  ! don't calculate subdiagonal matrix blocks
  logical  :: lsparse  = .true.  ! sparse representation for offdiagonal blocks
  !-------------------------------------------
  ! number of tests in subroutine constrain_ev
  !-------------------------------------------
  integer ,parameter :: mt = 7   ! number of tests

  !======================
  ! Contained Subroutines
  !======================
contains
!==============================================================================
  subroutine set_flags (obs, time, pb_apprx)
  type (t_obs)  ,intent(inout) :: obs (:)   ! observation data type
  type (t_time) ,intent(in)    :: time      ! analysis time
  integer       ,intent(in)    :: pb_apprx  ! approximation flag
  !-----------------------------------------------------------------
  ! This routine sets flags of component 'obs% spot%' used for the
  ! distribution of observation operators in a parallel environment:
  !
  ! pe_eval:
  ! cha:
  ! n_dest:
  ! pe_dest:
  !
  !-----------------------------------------------------------------
    !----------------
    ! local variables
    !----------------
    integer                :: ib, jb       ! box index (outer loop, lhs,rhs)
    integer                :: io, jo       ! observation index (lhs,rhs)
    type (t_spot) ,pointer :: oi, oj       ! pointer to observation
    real(wp)               :: xi(3), xj(3) ! coordinates
    real(wp)               :: s

    integer                :: dest(dace% npe) ! destination PEs

    !------------------
    ! printout (header)
    !------------------
    call flush_buf
    call add_line_pio ('')
    call add_line_pio ('  setting up communication flags for observations (H)')
    call add_line_pio ('')
    call add_line_pio ('  pe  box   spot char statid     n   dest...')
    !------------------------------------
    ! initialize DWD OI correlation model
    !------------------------------------
    if (pb_apprx > PB_BLOCKD) then
      call uncompress_cov
      call init_fg_cov (time)
    endif
    !----------------------------------------------------------------
    ! loop over boxes managed on this PE (rows of correlation matrix)
    !----------------------------------------------------------------
    do ib = 1, size(obs)
      if (dace% pe /= obs(ib)% pe) cycle
      !-----------------------
      ! loop over observations
      !-----------------------
      do io = 1, obs(ib)% n_spot
        oi => obs(ib)% spot(io)
        !---------------------------------------------
        ! mark observations to be evaluated on this PE
        !---------------------------------------------
        oi% pe_eval = obs(ib)% pe
        oi% char    = ior (oi% char, CHR_EVAL + CHR_REQ)
      end do
      !-----------------------------------------------------------------------
      ! loop over boxes not managed on this PE (columns of correlation matrix)
      !-----------------------------------------------------------------------
      if (pb_apprx > PB_BLOCKD) then
        do jb=1, size(obs)
          if (dace% pe == obs(jb)% pe)   cycle
          !---------------------------------
          ! loops over pairs of observations
          !---------------------------------
          do io = 1, obs(ib)% n_spot
            oi => obs(ib)% spot(io)
            do jo = 1, obs(jb)% n_spot
              oj => obs(jb)% spot(jo)
              !--------------------------------------------------------
              ! calculate great circle distance
              !   skip if there is no correlation between observations
              !--------------------------------------------------------
              select case (oi% int_type)
              case (ITY_ICOL)
                xi = oi% col% c% x
              case (ITY_MCOLS)
                xi = nearest_x (oj% col% c% x, oi)
              case default
                call finish('set_flags','invalid int_type')
              end select
              select case (oj% int_type)
              case (ITY_ICOL)
                xj = oj% col% c% x
              case (ITY_MCOLS)
                xj = nearest_x (oi% col% c% x, oj)
              case default
                call finish('set_flags','invalid int_type')
              end select
              s = sqrt(sum( (xi - xj)**2) )

              if (s >= x_cut*l_max) cycle
              !---------------
              ! mark receivers
              !---------------
              if (iand (oj% char, CHR_NONL+CHR_EXP) == CHR_NONL+CHR_EXP ) then
                !-----------------------------------------------
                ! expensive operators evaluated on some other PE
                !   H,R must be sent to this one
                !-----------------------------------------------
                oj% char    = ior (oj% char, CHR_COM + CHR_REQ)
                oj% pe_eval = obs(jb)% pe
!               !---------
!               ! printout
!               !---------
!               call nextline
!               write(oline(iol),'(i4,i5,i7,1x,a,2i4)') &
!                dace% pe, jb, oj% id, oj% statid, -99, &
!                oj% pe_eval
              else
                !--------------------------------------------------------------
                ! cheap operators are evaluated on all PEs (including this one)
                !--------------------------------------------------------------
                oj% char    = ior (oj% char, CHR_EVAL + CHR_REQ)
                oj% pe_eval = dace% pe
              endif
              !-------------
              ! mark senders
              !-------------
              if (iand (oi% char, CHR_NONL+CHR_EXP) == CHR_NONL+CHR_EXP ) then
                !-------------------------------------------
                ! expensive operators evaluated on this PE
                !   requiring communication with some other
                !-------------------------------------------
                oi% char = ior (oi% char, CHR_COM)
                if (.not.associated (oi% pe_dest)) then
                  allocate (oi% pe_dest (dace% npe))
                  oi% pe_dest (1:dace% npe) = -1
                end if
                if (any(oi% pe_dest (:oi% n_dest) == obs(jb)% pe)) cycle
                oi% n_dest = oi% n_dest + 1
                oi% pe_dest (oi% n_dest) = obs(jb)% pe
              endif
            end do
          end do
        end do
      endif
      !-----------------------------
      ! final loop over observations
      !-----------------------------
      do io = 1, obs(ib)% n_spot
        oi => obs(ib)% spot(io)
        if (oi% n_dest > 0 .and. oi% n_dest < dace% npe) then
          dest = oi% pe_dest
          deallocate (oi% pe_dest)
          allocate   (oi% pe_dest (oi% n_dest))
          oi% pe_dest = dest (:oi% n_dest)
          !---------
          ! printout
          !---------
          call nextline
          write(oline(iol),'(i4,i5,i7,i5,1x,a,11i4)')&
            dace% pe, ib, io, oi% char, oi% statid, oi% n_dest, &
            oi% pe_dest(1:min (oi% n_dest,10))
        endif
      end do
    end do
    if (pb_apprx > PB_BLOCKD) then
      call compress_cov
    endif
    !-------------------
    ! printout (trailer)
    !-------------------
    call flush_buf

  end subroutine set_flags
!------------------------------------------------------------------------------
  function nearest_x (x0,o) result (x)
  real(wp)                   :: x(3)
  real(wp)      ,intent(in)  :: x0(3)
  type (t_spot) ,intent(in)  :: o
    real(wp)            :: s(size(o% imcol))
    integer             :: i, in(1)
    do i=1,size(o% imcol)
      s(i) = sum( (x0 - o% imcol(i)%c% x)**2)
    end do
    in = minloc (s)
    x  = o% imcol(in(1))%c% x
  end function nearest_x
!==============================================================================
  subroutine scatter_K (obs)
  type(t_obs_set),target,intent(inout) :: obs     ! observation data type
    !----------------
    ! local variables
    !----------------
    integer                       :: ib ! box         index (outer loop, rows)
    integer                       :: pe ! destination processor index
    integer                       :: is ! observation index
    type (t_spot)        ,pointer :: si ! pointer to observation
    type (t_matrix_block),pointer :: pH ! pointer to matrix block
    type (t_matrix_block)         :: H  ! temporary
    !----------------------------------------------------------------
    ! loop over boxes managed on this PE (rows of correlation matrix)
    !----------------------------------------------------------------
    do ib = 1, size(obs% o)
      pH => obs%l%H%b(ib,ib)
      call p_bcast (pH% nonzero, pH% pe)
      !-----------------------
      ! loop over observations
      !-----------------------
      do is = 1, obs% o(ib)% n_spot
        si => obs% o(ib)% spot(is)
        if (iand (si% char, CHR_COM) /= 0) then
          !-----------------
          ! matrices to send
          !-----------------
          do pe = 1, si% n_dest
            call p_send (sub_block (pH, si%o%i, si%o%n, si%i%i, si%i%n),&
                         si% pe_dest(pe))
          end do
          !--------------------
          ! matrices to receive
          !--------------------
          if (si% n_dest == 0) then
            call p_recv (H, si% pe_eval)
            if (.not.associated (pH% packed)) then
              call allocate_block (pH, H% repr, ns=pH% nonzero)
              pH% ia = 1
            endif
            call insert (pH, H, si%o%i, si%i%i)
          end if
        end if
      end do
      if (associated (pH% packed)) call release_mem (pH)
    end do
  end subroutine scatter_K
!==============================================================================
  !------------------------------------------------------
  ! setup the DWD (OI) background error correlation model
  !------------------------------------------------------
  subroutine set_Pb (P_b, obs, er, time, i_ensb,               &
                     pb_apprx, tight, single_col, verbose,     &
                     e_f, e_fi, e_fie, set_e_fi, HHt, lBii, Bii)
  type(t_matrix)         ,intent(inout):: P_b       ! matrix HBHt
  type(t_obs_set)        ,intent(in)   :: obs       ! observation data type
  real(wp)               ,intent(in)   :: er        ! earth radius
  type (t_time)          ,intent(in)   :: time      ! analysis time
  integer                ,intent(in)   :: i_ensb    ! use ensemble B matrix
  integer       ,optional,intent(in)   :: pb_apprx  ! approximation flag
  logical       ,optional,intent(in)   :: tight     ! always mirror blocks
  logical       ,optional,intent(in)   :: single_col! single column appox.
  integer       ,optional,intent(in)   :: verbose   ! 1..3
  type(t_vector),optional,intent(inout):: e_f       ! forecast error(obsv.)
  type(t_vector),optional,intent(inout):: e_fi      ! forecast error(NMC,intp.)
  type(t_vector),optional,intent(inout):: e_fie     ! forecast error(NMC+ens.)
  logical       ,optional,intent(in)   :: set_e_fi  ! set e_fi
  logical       ,optional,intent(in)   :: HHt       ! flag to apply H,Ht
  logical       ,optional,intent(in)   :: lBii      ! flag to keep Bii
  type(t_matrix),optional,intent(inout):: Bii       ! B in interp.space
!------------------------------------------------------------------------------
! The purpose of this routine is to derive an explicit representation of the
! background error covariance matrix in observation space 'Pb' (i.e. HBHt),
! currently used for preconditioning the Conjugate Gradient solver in the 3dvar
! PSAS algorithm.
!
! Mandatory parameters are observation meta data 'obs' and analysis time 'time'
! (the underlying NMC statistics have an annual cycle).
!
!
! Optional parameters 'pb_apprx' and 'tight' further specify the representation
! of the matrix 'Pb'. 'Pb' is decomposed into blocks according to the
! partitioning of observations in boxes used in the PSAS scheme.
!
! supported values for 'pb_apprx' are:
!
!   PB_IDENT  = 2  Identity, only used for testing.
!
!   PB_DIAG   = 3  Diagonal, used in the monitoring step to derive
!                  background error variances.
!
!   PB_1DVAR  = 4  No correlations between reports, used in 1D-VAR mode
!
!   PB_BLOCKD = 5  Blockdiagonal matrix (default), used in the current PSAS
!                  implementation.
!
!   PB_FULL   = 6  Full matrix, used in the former PSAS implementation.
!
! Block diagonals are returned in PACKED representation (only the lower
! triangle is stored). Off-diagonal blocks (for pb_apprx==PB_FULL) are
! returned in one of the following representations:
!
!   FULL   full representation, every matrix element is stored.
!   ZERO   all elements are zero, no element is stored.
!   CSR    compressed sparse row storage, only non-zero elements are stored
!   CSC    compressed sparse column storage, only non-zero elements are stored
!   MIRROR matrix block is not stored,
!          equal to the transposed block mirrored at the diagonal
!
!
! Optionally the background errors (sqrt of the diagonal of Pb) are returned in
! 'e_f' (HBHt, observation space) and 'e_fi' (B, interpolation space).
!
! For former implementations of dummy sink variables the background error in
! interpolation space (input to the observation operator) is 1. (and its due
! to the observation operator to handle the transformation to physical
! quantities).  For recent implementations (Infra-Red emissivity
! Principal-Component coefficients) the actual error is stored in 'e_fi'. In
! this case 'e_fi' must be present as an input parameter in order to properly
! model the variances.
!
! If 'lBii' is .true. the the matrix 'Bii' (backgrond error covariance matrix in
! interpolation space) is returned as well, (only for radiances and covariances
! within one field of view).
!
! If 'single_col' is set multi column operators (GPSRO, GPSGB) are handled
! approximately (without seperating different model columns) in the static
! NMC-B-matrix. This works only if USE_GET_CORR_OPT is not set.
!------------------------------------------------------------------------------

    !---------------------------------------------------
    ! local variables: loop indices, pointer and counter
    !---------------------------------------------------
    integer                      :: ib, jb       ! matrix block indices
    integer                      :: io, jo       ! observation type indices
    integer                      :: i, j         ! observation indices
    integer                      :: k, l         ! indices
    integer                      :: ii,jj,iii    ! spot indices
    integer                      :: iiii         ! observation index
    integer                      :: ki,kj        ! spot sizes (levels)
    integer(i8)                  :: nn0          ! nonzero elements
    integer(i8)                  :: nntot        ! total no. of elements
    integer                      :: nn0b         ! nonzero elements in box
    integer                      :: nonz         ! nonzero elems. from get_corr
    type(t_matrix_block),pointer :: b            ! pointer to block of P_b
    type(t_matrix_block),pointer :: Hi           ! pointer to block of H
    type(t_matrix_block),pointer :: Hj           ! pointer to block of H
    type(t_obs)         ,pointer :: bi, bj       ! observation box pointer
    real(wp)                     :: xi(3), xj(3) ! coordinates
    type(t_spot)        ,pointer :: oi, oj       ! observation type pointer
    integer                      :: is, js       ! spot indices
    logical                      :: lident       ! int == psas
    integer         ,allocatable :: line (:)     ! for printout
    integer                      :: pe           ! mpi processor index
    type(t_rowcol)  ,pointer     :: rc(:)        ! precalculated for row/col
    type(t_rowcol)  ,pointer     :: lhs,rhs      ! precalculated for row/col
    integer(i8)                  :: nonzero_bdpe ! elements/=0 on bdiag.   PE
    integer(i8)                  :: nonzero_odpe ! elements/=0 on offdiag. PE
    integer(i8)                  :: bsize        ! size of matrix block
    logical                      :: lHHt         ! flag to apply H,Ht
    logical                      :: llBii        ! flag to keep Bii
    integer                      :: nBii         ! number of elements in Bii
    logical                      :: ltight       ! flag to always mirror blocks
    integer                      :: pb_apr       ! approximation flag
    integer                      :: ids          ! dummy sink variable index
    logical                      :: lset_efi     ! set e_fi
    type(t_vector) ,allocatable  :: xy(:)        ! ensemble deviations
    type(t_coord)  ,pointer      :: ci           ! lhs coordinate
    type(t_coord)  ,pointer      :: cj           ! rhs coordinate
    logical                      :: lsc          ! single column approximation
    integer                      :: verb         ! verbosity level
    !-------------------------------------------------------------
    ! local variables: quantities used for the covariance matrices
    !-------------------------------------------------------------
    real(wp) :: d       ! distance in m
#if !defined (USE_GET_CORR_OPT)
    real(wp) :: s       ! secant of the unit sphere
#endif
    real(wp) :: s2      ! s^2
    real(wp) :: s2max   ! maximum value of s2 for non-vanishing correlations
    real(wp) :: dmin    ! minimum distance between spots in a box (diagnostics)
    real(wp) :: dmax    ! maximum distance between spots in a box (diagnostics)
    !----------------------------------------------
    ! temporary arrays to hold correlation matrices
    !---------------------------------------------
#if !defined (USE_GET_CORR_OPT)
    real(mp) ,allocatable :: coix (:,:)
#endif
    real(mp) ,pointer     :: coi  (:,:)
    real(mp) ,pointer     :: cop  (:,:)
    real(mp) ,allocatable :: coip (:,:)
#if defined (USE_GET_CORR_OPT)
    ! NECJB
    integer :: &
      nv_spots, mv_spots, iv, iv0, nv, &
      ni_coi, nj_coi, &
      ni_cox, nj_cox, &
      ni_spt, nj_spt, i0, j0, nv_b
#if defined (__NEC__)
    integer, parameter    :: nv_blocksize = 2048
!   restrict size of array rv_coi to allocate:
    integer, parameter    :: max_size =  512*1024*1024 ! 512M words
#else
    ! Conserve memory on scalar architectures:
    integer, parameter    :: nv_blocksize = 256
!   restrict size of array rv_coi to allocate:
    integer, parameter    :: max_size =   32*1024*1024 !  32M words
#endif
    integer,  allocatable :: iv_spots(:,:)
    real(wp), allocatable :: rv_dist (:)
    real(mp), allocatable :: rv_coi  (:,:,:,:,:)
    integer,  allocatable :: iv_nonz (:,:,:)
#endif
#if defined (DEBUG_GET_CORR_OPT)
    integer :: nonz_
#endif
    !------------------------------------
    ! initialize DWD OI correlation model
    !------------------------------------
    if (P_b% m /= P_b% n) call finish ('set_Pb','matrix is not square')
    call init_fg_cov (time)
    !------------------
    ! print header line
    !------------------
    if (dace% lpio) then
      write(6,'("----------------------------------------------------------")')
    endif
    !----------------------------------------
    ! allocate arrays, precompile information
    !----------------------------------------
    if (test_bl) lpack = .false.
    lHHt     = .true.    ;if (present(HHt))        lHHt     = HHt
    llBii    = .false.   ;if (present(lBii))       llBii    = lBii
    ltight   = .true.    ;if (present(tight))      ltight   = tight
    pb_apr   = PB_BLOCKD ;if (present(pb_apprx))   pb_apr   = pb_apprx
    lset_efi = .false.   ;if (present(set_e_fi))   lset_efi = set_e_fi
    lsc      = .false.   ;if (present(single_col)) lsc      = single_col
    verb     = 1         ;if (present(verbose))    verb     = verbose
    call uncompress_cov

    if (llBii) then
      if (.not.present(Bii)) call finish('set_Pb','Bii not present !')
      call destruct  (Bii)
      call construct (Bii, obs% ii, 'Bii', BDIAG)
    endif
    allocate (rc(P_b% m_b))
!   if (pb_apr > PB_BLOCKD) then
!     do ib=1,P_b% m_b
!       if (P_b% b (ib,ib)% pe == dace% pe) call set_rc (ib)
!     end do
!   endif

    !------------------------------------
    ! account for ensemble B contribution
    !------------------------------------
    call EnKF_precond_0 (xy, obs, i_ensb)

    !----------------
    ! loop over boxes
    !----------------
    nn0 = 0
    do jb = 1, P_b% n_b
      bj  => obs% o(jb)
      Hj  => obs% l% H% b (jb,jb)
      rhs => rc  (jb)
      do ib = 1, P_b% m_b
        b  => P_b% b(ib,jb)
        !--------------------------------------------
        ! skip if matrix block already set (mirrored)
        !                      or located on other PE
        ! else fill this block
        !--------------------------------------------
        if (b% repr >= 0)        cycle
        if (b% pe   /= dace% pe) cycle
        if (P_b% b(jb,ib)% pe == dace% pe .and. ltight .and. lmirror) then
          if ((ib > jb)) then
            P_b% b(ib,jb)% repr = MIRROR
            P_b% b(ib,jb)% nonzero = 0
            cycle
          endif
        endif
        if (P_b% b(jb,ib)% pe /= dace% pe .and. ltight .and. lmirror) then
          if ((ib > jb .and. mod(ib+jb,2)==1) .or. &
              (ib < jb .and. mod(ib+jb,2)==0)) then
                  P_b% b(ib,jb)% repr = MIRROR
                  P_b% rowwise = .false.
                  P_b% b(ib,jb)% nonzero = 0
                  cycle
          endif
        endif
        nn0b = 0
        if (ib==jb .or. pb_apr > PB_BLOCKD) then
         bi  => obs% o(ib)
         Hi  => obs% l% H% b (ib,ib)
         dmin = huge(1._wp)
         dmax = 0._wp
         !-------------------------------------------
         ! allocate Bii elements (for RADiances only)
         !-------------------------------------------
         nBii = 0
         if (llBii .and. ib==jb) then
           do i = 1, bi% n_spot
             if (bi% spot(i)% hd% obstype /= OT_RAD) cycle
             nBii = nBii + bi% spot(i)% i% n ** 2
           end do
           if (nBii > 0) &
             call allocate_block (Bii% b(ib,jb), CSC, ns = nBii)
         endif
         !-----------------------------------------
         ! approximation : identity matrix for HPbH
         !-----------------------------------------
         select case (pb_apr)
         case (PB_ZERO)
           b% repr = ZERO
           nn0b    = 0
         case (PB_DUMMY)
           call allocate_block (b, DIAGONAL)
           b% packed = 0._mp
           nn0b      = 0
!NEC$ ivdep
           do i = 1, b% m
             if (bi% t_int (i) == OBS_DUM) then
               b% packed(i) = 1._mp
               nn0b         = nn0b + 1
             endif
           end do
         case (PB_IDENT)
           call allocate_block (b, FULL)
           b% full = 0._mp
           nn0b = b% m
           do i=1,nn0b
             b% full(i,i) = 1._mp
           end do
         case default
         !-------------------------------
         ! general case, no approximation
         !-------------------------------
!         if (bi% n_lev == 0) call col2lev (bi)
!         if (bj% n_lev == 0) call col2lev (bj)
          select case (pb_apr)
          case (PB_DIAG)
            call allocate_block (b, DIAGONAL)
          case default
            call allocate_block (b, FULL)
            b% full = 0._mp
          end select
          !--------------------------------
          ! set rhs (background error info)
          !--------------------------------
          lhs => rc  (ib)
          if (.not.associated (rhs% cols)) call set_rc (jb)
          if (.not.associated (lhs% cols)) call set_rc (ib)
          !-------------------------------------------------
          ! loop over spots in the box, lower triangles only
          !-------------------------------------------------

!write (0,*)dace% pe,'ib,jb,bj% n_spot,assoc',ib,jb,bj% n_spot,associated(rhs%cols)

#if defined (USE_GET_CORR_OPT)
          nv_spots = 0
          mv_spots = bj% n_spot * bi% n_spot

          allocate( iv_spots(mv_spots,  2) )
          allocate( rv_dist (mv_spots    ) )
#endif
          s2max = (x_cut*l_max)**2
          do jo = 1, bj% n_spot
            oj => bj% spot(jo)
            js = oj% is

!if (.not.associated (rhs% cols))then
!write (0,*) dace% pe,'jo, jb, js: allocate',jo, jb, js
!if (.not.associated (rhs% cols)) call set_rc (jb)
!write (0,*) dace% pe,'jo, jb, js: allocate',jo, jb, js
!endif

            do io = bi% n_spot, 1, -1                  ! for 'L'

!if (.not.associated (rhs% cols)) call set_rc (jb)

              if (jb==ib .and. io<jo .and. lpack) exit ! representation

!if (.not.associated (rhs% cols)) call set_rc (jb)

!           do io = 1, bi% n_spot                      ! for 'U'
!             if (jb==ib .and. io>jo .and. lpack) exit
              if (pb_apr <= PB_1DVAR .and. io/=jo)            cycle
              oi => bi% spot(io)
              is = oi% is
              if (pb_apr <= PB_1DVAR .and. oi% use% state <= STAT_DISMISS) cycle
              !------------------------------
              ! finally calculate covariances
              ! (univariate so far)
              !------------------------------
              !--------------------------------------------------------
              ! calculate great circle distance
              !   skip if there is no correlation between observations
              !--------------------------------------------------------
              select case (oi% int_type)
              case (ITY_ICOL)
                xi = oi% col% c% x
              case (ITY_MCOLS)
                xi = nearest_x (oj% col% c% x, oi)
              case default
                call finish('set_Pb','invalid int_type')
              end select
              select case (oj% int_type)
              case (ITY_ICOL)
                xj = oj% col% c% x
              case (ITY_MCOLS)
                xj = nearest_x (oi% col% c% x, oj)
              case default
                call finish('set_Pb','invalid int_type')
              end select
#if !defined (USE_GET_CORR_OPT)
              s2 = sum ((xi - xj)**2)
              if (s2 >= s2max .and. i_ensb == 0) cycle
              s = sqrt (s2)
              d = 2._wp *er *asin (.5_wp * s)
              if(jb/=ib .or. io/=jo) then
                dmin=min(d,dmin)
                dmax=max(d,dmax)
              endif
#else
              s2 = (xi(1) - xj(1))**2 + &
                   (xi(2) - xj(2))**2 + &
                   (xi(3) - xj(3))**2
              if (s2 < s2max .or. i_ensb > 0) then
                nv_spots = nv_spots + 1
                iv_spots(nv_spots,1) = io
                iv_spots(nv_spots,2) = jo
                rv_dist (nv_spots  ) = s2
              end if
            end do
          end do

          nv = nv_spots
          ni_spt = 0
          nj_spt = 0
          ni_coi = 0
          nj_coi = 0
          ni_cox = 0
          nj_cox = 0
          do iv = 1, nv
            io = iv_spots(iv,1)
            jo = iv_spots(iv,2)
            d = 2._wp *er *asin (.5_wp * sqrt(rv_dist(iv)))
            if(jb/=ib .or. io/=jo) then
              dmin=min(d,dmin)
              dmax=max(d,dmax)
            endif
           !oi => bi% spot(io)
           !oj => bj% spot(jo)
#define     oi    bi% spot(io)
#define     oj    bj% spot(jo)
            !-----------------------------------
            ! Derive max. number of columns/spot
            !-----------------------------------
            ni_spt = max (ni_spt, oi% n_spt)
            nj_spt = max (nj_spt, oj% n_spt)
            !----------------------------------------------------------
            ! Derive dimensions of array rv_coi passed to get_corr_opt,
            ! taking into account that multiple columns/spot
            ! (radiooccultations!) will be split.
            !----------------------------------------------------------
            if (oi% int_type == ITY_MCOLS) then
               ni_cox = max (ni_cox, rc(ib)% spots(oi% is + 1)% n)
            else
               ni_coi = max (ni_coi, oi%i% n)
            end if
            if (oj% int_type == ITY_MCOLS) then
               nj_cox = max (nj_cox, rc(jb)% spots(oj% is + 1)% n)
            else
               nj_coi = max (nj_coi, oj%i% n)
            end if
          end do
#undef oi
#undef oj

#if defined (DEBUG_GET_CORR_OPT)
          write(0,*) dace% pe, ": set_Pb: ni_spt, nj_spt =", ni_spt, nj_spt
          write(0,*) dace% pe, ": set_Pb: ni_coi, nj_coi =", ni_coi, nj_coi
          write(0,*) dace% pe, ": set_Pb: ni_cox, nj_cox =", ni_cox, nj_cox
#endif

          ni_coi = max (ni_coi, ni_cox)
          nj_coi = max (nj_coi, nj_cox)

#if defined (DEBUG_GET_CORR_OPT)
          write(0,*) dace% pe, ": set_Pb: ni_coi, nj_coi =", ni_coi, nj_coi
#endif

!NECJB Precompute coi
!         write(0,*) dace% pe,'line,ni_coi,nj_coi,nv,ni_spt,nj_spt=', &
!                 __LINE__,ni_coi,nj_coi,nv,ni_spt,nj_spt
          ! restrict size of array rv_coi to allocate:
          nv_b = nv_blocksize
          do while (int(ni_coi*nj_coi*ni_spt*nj_spt,i8) * nv_b > max_size)
            nv_b = nv_b / 2
          end do
          nv_b = max (nv_b,1)
          allocate (rv_coi (ni_coi,nj_coi,nv_b,ni_spt,nj_spt))
          allocate (iv_nonz(nv_b,ni_spt,nj_spt))
          allocate (tv_gc  (nv_b*ni_spt*nj_spt))
          rv_coi  = 0._mp

          iv0 = 0
          do while( iv0 < nv )

            !---------
            ! loop lhs
            !---------
            do i0 = 1, ni_spt
              !---------
              ! loop rhs
              !---------
              do j0 = 1, nj_spt

                !call ftrace_region_begin( 'set_pb.1' )
                nv_gc = 0
!NEC$ ivdep
                do iv = iv0+1, min( nv, iv0+nv_b )
                  io = iv_spots(iv,1)
                  jo = iv_spots(iv,2)
                 !oi => bi% spot(io)
                 !oj => bj% spot(jo)
#define           oi    bi%spot(io)
#define           oj    bj%spot(jo)
                  !----------------
                  ! set covariances
                  !----------------
                  if( i0 <= oi% n_spt .and. &
                      j0 <= oj% n_spt         )then
                    is = oi% is
                    js = oj% is
                    i  = is+i0
                    j  = js+j0

                    nv_gc = nv_gc + 1
                    tv_gc(nv_gc)%i  = i
                    tv_gc(nv_gc)%j  = j
                    tv_gc(nv_gc)%iv = iv-iv0
                  end if
                end do
                !call ftrace_region_end  ( 'set_pb.1' )

                call get_corr_opt (lhs, rhs, ib==jb, i0, j0, iv_nonz, rv_coi )

              end do
            end do
#undef oi
#undef oj

    !call ftrace_region_begin( 'set_pb.2' )

            do iv = iv0+1, min( nv, iv0+nv_b )
              io =  iv_spots(iv,1)
              jo =  iv_spots(iv,2)
              oi => bi% spot(io)
              oj => bj% spot(jo)
              is =  oi% is
              js =  oj% is
#endif
              !-----------------------------------------------------------
              ! allocate matrices for interpolated and observational space
              !-----------------------------------------------------------
              lident = iand(oi% char, CHR_ID)/=0 &
                 .and. iand(oj% char, CHR_ID)/=0 &
                  .or. .not. lHHt
              select case (pb_apr)
              case (PB_DIAG)
                allocate (cop(oi%o% n,oj%o% n))
                cop = 0._mp
              case default
                if (lHHt) then
                  cop => b% full (oi%o% i +       1 : &
                                  oi%o% i + oi%o% n , &
                                  oj%o% i +       1 : &
                                  oj%o% i + oj%o% n )
                else
                  cop => b% full (oi%i% i +       1 : &
                                  oi%i% i + oi%i% n , &
                                  oj%i% i +       1 : &
                                  oj%i% i + oj%i% n )
                endif
              end select
              if (lident) then
                coi => cop
              else
#if defined (USE_GET_CORR_OPT)
                allocate (coi (oi%i% n, oj%i% n ))
                coi = 0._mp
#endif
                allocate (coip(oi%i% n, oj%o% n ))
                coip = 0._mp
              endif
              !----------------
              ! set covariances
              !----------------
              nonz = 0
              !---------
              ! loop lhs
              !---------
              ii = 0
              do i = is+1, is+oi% n_spt
                ki =  oi%i% n
                ci => oi% col% c
                if (oi% int_type == ITY_MCOLS) then
                  ki =  rc(ib)% spots(i)% n
                  ci => oi% imcol (i-is)% c
                endif
#if !defined (USE_GET_CORR_OPT)
                if (.not.lident) then
                  allocate (coi (ii+1:ii+ki, oj%i% n ))
                  coi = 0._mp
                endif
#endif
                !---------
                ! loop rhs
                !---------
                jj = 0
                do j = js+1, js+oj% n_spt
                  kj =  oj%i% n
                  cj => oj% col% c
                  if (oj% int_type == ITY_MCOLS) then
                    kj =  rc(jb)% spots(j)% n
                    cj => oj% imcol (j-js)% c
                  endif
#if !defined (USE_GET_CORR_OPT)

                  if (lsc) then
                    if (ii == 0 .and. jj == 0) then
                      call get_corr (lhs, rhs, i, j, .TRUE., ib==jb, nonz, &
                                     coi(ii+1:ii+ki,jj+1:jj+kj)            )
                      if (allocated (coix)) deallocate (coix)
                      if (oi% n_spt > 1 .or. oj% n_spt > 1) then
                        allocate (coix (ki,kj))
!NEC$ ivdep
                        coix = coi(ii+1:ii+ki,jj+1:jj+kj)
                      endif
                    else
!NEC$ ivdep
                      coi(ii+1:ii+ki,jj+1:jj+kj) = coix
                    endif
                  else
                    call get_corr (lhs, rhs, i, j, .TRUE., ib==jb, nonz, &
                                   coi(ii+1:ii+ki,jj+1:jj+kj)            )
                  endif
#else
#if !defined (DEBUG_GET_CORR_OPT)
                  coi(ii+1:ii+ki,jj+1:jj+kj) = &
                      rv_coi(1:ki,1:kj,iv-iv0,i-is,j-js)
                  nonz = nonz + iv_nonz(iv-iv0,i-is,j-js)
#else
                  ! Debugging aid for get_corr_opt:
                  nonz_ = nonz
                  call get_corr (lhs, rhs, i, j, .TRUE., ib==jb, nonz, &
                                 coi(ii+1:ii+ki,jj+1:jj+kj)            )
!write(0,*) "rv_coi limits:", ki, kj,iv-iv0,i-is,j-js
                  s2 = maxval (abs (   coi(ii+1:ii+ki,jj+1:jj+kj)     - &
                                    rv_coi(1:ki,1:kj,iv-iv0,i-is,j-js)) )
                  if (        s2 >  1.e-16                   .or. &
                      nonz-nonz_ /= iv_nonz(iv-iv0,i-is,j-js)     ) then
                     write(0,*) "Check get_corr vs. get_corr_opt!"
                     write(*,*) "==="
                     write(*,*) "get_corr, inner loop: i, j =", i, j
                     write(*,*) "nonz:", nonz-nonz_, iv_nonz(iv-iv0,i-is,j-js)
                     write(*,*) "delta =", s2
                     write(*,*) "   coi=",    coi(ii+1:ii+ki,jj+1:jj+kj)
                     write(*,*) "rv_coi=", rv_coi(1:ki,1:kj,iv-iv0,i-is,j-js)
                     write(*,*) "==="
                  end if
#endif /* DEBUG_GET_CORR_OPT */
#endif /* USE_GET_CORR_OPT   */
                  !----------------------------
                  ! fix dummy diagonal elements
                  !----------------------------
                  if (ib==jb .and. i==j .and. present(e_fi)) then
                    !------------------------------------------------
                    ! adjust background error (linear sink variables)
                    !------------------------------------------------
!NEC$ ivdep
                    do iii = ii+1, ii+ki
                      iiii = iii + oi%i%i
                      if (obs% o(ib)% t_int(iiii) == OBS_DUM) then
                        coi(iii,iii) = e_fi% s(ib)% x(iiii)**2
                      endif
                    end do
                    !--------------------------------------------
                    ! account for variation of Jakobian
                    ! (in nonlinear sink variable transformation)
                    !--------------------------------------------
!NEC$ ivdep
                    do ids = oi%d%i + 1, oi%d%i + oi%d%n
                      iii  = bi% sink(ids)% iobs
                      iiii = iii + oi%i%i
                      coi(iii,iii) =  bc_bd (obs% l% x% s(ib)% x(iiii), bi% sink(ids))
                    end do
                  endif
                  !-------------------------------------------
                  ! account for ensemble B matrix contribution
                  !-------------------------------------------
                  if (i_ensb == 1) then
                    call EnKF_precond_1 (coi(ii+1:ii+ki,jj+1:jj+kj), xy,    &
                                         obs, oi, oj, ci, cj, ib, jb, ii, jj)
                    nonz = size (coi)
                  endif
                  !----------------------------
                  ! fix dummy diagonal elements
                  !----------------------------
                  if (ib==jb .and. i==j .and. nl_prec) then
                    !------------------------------------------------
                    ! adjust background error (linear sink variables)
                    !------------------------------------------------
!NEC$ ivdep
                    do iii = ii+1, ii+ki
                      iiii = iii + oi%i%i
                      !----------------------------------
                      ! account for variation of Jakobian
                      ! (for generalised humidity)
                      !----------------------------------
                      if (obs% o(ib)% t_int(iiii) == OBS_RH) then
                        coi(iii,iii) =  bc_bv (obs% l% x % s(ib)% x(iiii), &
                                               obs% b% xb% s(ib)% x(iiii), &
                                               real(coi(iii,iii),wp),      &
                                               gh_meta                     )
                      endif
                    end do
                  endif
                  !---------------------------------------------------------
                  ! store interpolated diagonal background error if required
                  !---------------------------------------------------------
                  if (ib==jb .and. i==j .and. lset_efi) then
                    if (present (e_fie)) then
!NEC$ ivdep
                      do k = ii+1, ii+ki
                        e_fie% s(ib)% x (oi%i%i + k) = coi (k,k)
                      end do
                    end if
                  end if
                  !-------------------------------------
                  ! store Bii block matrices if required
                  !-------------------------------------
                  if (nBii > 0 .and. io==jo                    &
                               .and. oi% hd% obstype == OT_RAD ) then
                   call insert (Bii% b(ib,jb), coi, oi%i%i, oj%i%i, use_zero=.true.)
                  endif
                  jj = jj + kj
                end do
                !---------------------------------------------------
                ! transform interpolated to observational space, rhs
                !---------------------------------------------------
                if (.not. lident) then
                  if (iand(oj% char, CHR_ID)/=0) then
!NEC$ ivdep
                    coip(ii+1:ii+ki,:) = coi(ii+1:ii+ki,:)
                  else
!                   !--------------------
!                   ! general formulation
!                   !--------------------
!                   Hs = sub_block (Hj, oj%o%i, oj%o%n, oj%i%i, oj%i%n)
!                   do iii=ii+1, ii+ki
!                     coip(iii,:) = Hs * real(coi(iii,:),wp)
!                   end do
!                   call destruct (Hs)
FTRACE_BEGIN("set_Pb:rhs")
                    !-----------------------------------------
                    ! efficient formulation (CSC storage only)
                    !-----------------------------------------
                    do j=1,oj%i%n
                      do k = Hj%ia(oj%i%i+j), Hj%ia(oj%i%i+j+1)-1
                        l = Hj%ja(k)-oj%o%i
!DIR$ IVDEP
!NEC$ IVDEP
                        do iii=ii+1, ii+ki
                          coip(iii,l) = coip(iii,l) + Hj%packed(k) * coi(iii,j)
                        end do
                      end do
                    end do
FTRACE_END  ("set_Pb:rhs")

                  endif
#if !defined (USE_GET_CORR_OPT)
                  deallocate (coi)
#endif
                endif
                ii = ii + ki
              end do
              if (nonz/=0) then
                !---------------------------------------------------
                ! transform interpolated to observational space, lhs
                !---------------------------------------------------
                if (.not. lident) then
!                 !--------------------
!                 ! general formulation
!                 !--------------------
!                 Hs = sub_block (Hi, oi%o%i, oi%o%n, oi%i%i, oi%i%n)
!                 do j=1,oj%o% n
!                   cop (:,j) = Hs * coip(:,j)
!                 end do
!                 call destruct (Hs)
FTRACE_BEGIN("set_Pb:lhs")
                  !-----------------------------------------
                  ! efficient formulation (CSC storage only)
                  !-----------------------------------------
                  do i=1,oi%i%n
#ifdef NO_LOOP_INTERCHANGE
!NEC$ nointerchange
!DIR$ IVDEP
!NEC$ IVDEP
                    do k = Hi%ia(oi%i%i+i), Hi%ia(oi%i%i+i+1)-1
                      l = Hi%ja(k)-oi%o%i
!NEC$ loop_count_test
!DIR$ IVDEP
                      do j=1,oj%o% n
#else /* manual loop interchange for better L1 cache usage */
!NEC$ loop_count_test
!DIR$ IVDEP
                      do j=1,oj%o% n
!NEC$ nointerchange
!DIR$ IVDEP
!NEC$ IVDEP
                    do k = Hi%ia(oi%i%i+i), Hi%ia(oi%i%i+i+1)-1
                      l = Hi%ja(k)-oi%o%i
#endif
                        cop (l,j) = cop (l,j) + Hi%packed(k) * coip(i,j)
                      end do
                    end do
                  end do
FTRACE_END  ("set_Pb:lhs")

                  nn0b = nn0b + count (cop /= 0._mp)
                else
                  nn0b = nn0b + nonz
                endif
                !--------------------------------------------
                ! explicitly store block matrices if required
                !--------------------------------------------
                select case (pb_apr)
                case (PB_DIAG)
!NEC$ ivdep
                  do i = 1, oi%o% n
                    b% packed(oi%o% i+i) = cop (i,i)
                  end do
                  deallocate (cop)
                end select
              endif
              !-------------------------------
              ! add ensemble B if not yet done
              !-------------------------------
              if (i_ensb == 2) call EnKF_precond_2 (b, &
                        xy, obs, oi, oj, ib, jb)
              !-----------------------
              ! deallocate temporaries
              !-----------------------
              if (.not.lident) deallocate (coip)
#if !defined (USE_GET_CORR_OPT)
            end do
          end do
#else
              if (.not.lident) deallocate (coi)
            end do
            iv0 = iv0 + nv_b
    !call ftrace_region_end  ( 'set_pb.2' )
          end do
          deallocate( rv_coi   )
          deallocate( iv_spots )
          deallocate( rv_dist  )
          deallocate( iv_nonz  )
          deallocate( tv_gc    )
#endif
          select case (pb_apr)
          case (PB_DIAG)
            nn0b = count (b% packed /= 0._mp)
          case default
            if (test_bl .and. ib==jb) call test_c (real(b% full,wp), ib)
            nn0b = count (b% full /= 0._mp)
            if (nn0b==0) deallocate (b% full)
          end select
         end select
        endif
        !-----------------------
        ! deallocate memory asap
        !-----------------------
        if (ib /= jb) call destruct (rc(ib))
        !-----------------------
        ! return forecast errors
        !-----------------------
        if (present(e_f)) then
          if(jb==ib .and. e_f% s(ib)% pe == dace% pe) then
            select case (pb_apr)
            case (PB_DIAG)
!NEC$ ivdep
              do i = 1, bi% n_obs
                e_f% s(ib)% x (i) = sqrt (real (b% packed(i),wp))
              end do
            case default
!NEC$ ivdep
              do i = 1, bi% n_obs
                if(nn0b==0) then
                  e_f% s(ib)% x (i) = 0._wp
                else
                  e_f% s(ib)% x (i) = sqrt (real (b% full(i,i),wp))
                endif
              end do
            end select
          endif
        endif

        !----------------------------------
        ! store matrix in compressed format
        !----------------------------------
        if (nn0b==0) then
          !------------------------------------
          ! matrix block only consists of zeros
          !------------------------------------
          b% repr = ZERO
          if (lmirror .and. P_b% b(jb,ib)% pe == dace% pe) &
            P_b% b(jb,ib)% repr = ZERO
        else
          !-------------------------
          ! nonzero elements present
          !-------------------------
          if (jb == ib) then
            !----------------------------------------
            ! use packed format for diagonal elements
            !   (only store lower triangle)
            !----------------------------------------
            select case (pb_apr)
            case (PB_DIAG)
            case default
              if (lpack) call pack_matrix (b)
              nn0b = b% nonzero
            end select
          else
!            !--------------------------------------------------------
!            ! dont store subdiagonal blocks of symmetric matrix
!            !--------------------------------------------------------
!            if (lmirror) then
!              if (P_b% b(jb,ib)% pe == dace% pe) P_b% b(jb,ib)% repr = MIRROR
!            endif
            !---------------------------------------------------------------
            ! use compressed sparse row/column format for offdiagonal blocks
            !---------------------------------------------------------------
            if (lsparse) then
              if (ib>jb) then
                call csr_matrix(P_b% b(ib,jb))
              else
                call csc_matrix(P_b% b(ib,jb))
              endif
            endif
          endif
          !---------------
          ! print one line
          !---------------
          call nextline
          if (dmin == huge(1._wp)) dmin = 0._wp
          bsize =         bj% n_int
          bsize = bsize * bi% n_int
          if (bsize > 0) then
            write(oline(iol),'(3i5,2i7,i11,i15,f7.2,2f12.3)') &
              dace% pe, ib, jb, b%m, b%n, nn0b, bsize,            &
              100._wp * nn0b / bsize, dmin/1000, dmax/1000
          else
            write(oline(iol),'(3i5,2i7,i11,i15)') &
              dace% pe, ib, jb, b%m, b%n, nn0b, bsize
          endif
        endif
        nn0 = nn0 + nn0b
      end do
      if (dace% lpio) then
        write(6,*) 'set_Pb: block',jb,'/',P_b% n_b,'processed.'
      endif
      !-----------------------
      ! deallocate memory asap
      !-----------------------
!     if (Hj% pe /=dace% pe .or. pb_apr <= PB_BLOCKD) call destruct (rc(jb))
                                                      call destruct (rc(jb))
      if (Hj% pe /=dace% pe)                          call deallocate (Hj)
    end do
    !--------------
    ! print trailer
    !--------------
    if (dace% lpio) then
      write(6,'()')
      write(6,&
 '("   pe     block      m      n    nonzero          total  %         min    max distance")')
      write(6,'()')
    endif
    call flush_buf
    call nextline
    nntot = int (P_b%m,i8) * int (P_b%n,i8)
    write(oline(iol),'(i4,"       total",12x,i12,i15,f7.2)')     &
      dace% pe, nn0, nntot,                                      &
      min (100._wp,max (0._wp,100*real (nn0,wp)/max (nntot,1_i8)))
    call flush_buf
    nn0 = p_sum (nn0)
    if (dace% lpio) then
      write(6,'(" all       total",12x,i12,i15,f7.2)')             &
        nn0, nntot,                                                &
        min (100._wp,max (0._wp,100*real (nn0,wp)/max (nntot,1_i8)))
      write(6,'()')
      write(6,'()')
      write(6,*) repeat('-',79)
    endif
    if (dace% lpio .and. verb > 2) then
      write(6,'()')
      write(6,'(a)') '  Matrix Block encoding:'
      write(6,'()')
      write(6,'(4x,a1," : ",a)')crep(ZERO),  'all elements are zero'
      write(6,'(4x,a1," : ",a)')crep(FULL),  'full n x m representation'
      write(6,'(4x,a1," : ",a)')crep(PACKED),'only lower triangle is stored'
      write(6,'(4x,a1," : ",a)')crep(CSR),   'compressed sparse row storage'
      write(6,'(4x,a1," : ",a)')crep(CSC),   'compressed sparse column storage'
      write(6,'(4x,a1," : ",a)')crep(MIRROR),'not stored due to symmetry'
      write(6,'()')
      write(6,'(a)') ' PE'
    endif
    !----------------------------
    ! write block matrix encoding
    !----------------------------
    allocate (line(P_b% n_b))
    if(P_b% n_b <=256 .and. verb > 2) then
      do ib = 1, P_b% m_b
        pe = P_b% b(ib,ib)% pe
        !
        ! pack and send
        !
        if (pe == dace% pe) then
          line = P_b% b(ib,:)% repr
          if (pe /= dace% pio) call p_send (line, dace% pio, 1)
        endif
        !
        ! recv and write
        !
        if (dace% lpio) then
          if (pe /= dace% pio) call p_recv (line, pe, 1)
          write(6,'(i3,1x,256a1)') pe, crep(line)
        endif
      end do
      !
      ! write block matrix usage
      !
      if (dace% lpio) then
        write(6,'()')
        write(6,*) repeat('-',79)
        write(6,'()')
        write(6,'(a)') '  Matrix Block usage'
        write(6,'()')
        write(6,'(4x,i1," : ",a)')  0, '         empty'
        write(6,'(4x,i1," : ",a)')  1, ' <= 10 % filled'
        write(6,'(4x,   "..."  )')
        write(6,'(4x,i1," : ",a)')  9, ' <= 90 % filled'
        write(6,'(4x,i1," : ",a)') 10, ' >  90 % filled'
        write(6,'()')
        write(6,'(a)') ' PE'
      endif
      do ib = 1, P_b% m_b
        pe = P_b% b(ib,ib)% pe
        !
        ! pack and send
        !
        if (pe == dace% pe) then
          line = P_b% b(ib,:)% n * P_b% b(ib,:)% m
          where (line>0)
#if defined(__ibm__) || defined (_CRAYFTN)
            ! Work around strange bug with xlf 9.1 and bounds-checking:
            line = ceiling(float(10*P_b% b(ib,:)% nonzero) / max (line,1))
#else
            line = ceiling(float(10*P_b% b(ib,:)% nonzero)/line)
#endif
          elsewhere
            line = 0
          end where
          if (pe /= dace% pio) call p_send (line, dace% pio, 1)
        endif
        !
        ! recv and write
        !
        if (dace% lpio) then
          if (pe /= dace% pio) call p_recv (line, pe, 1)
          write(6,'(i3,1x,256i1)') pe, line
        endif
      end do
      if (dace% lpio) write(6,'()')
    endif
    deallocate (line)
    !------------------
    ! deallocate arrays
    !------------------
    call destruct (rc)
    deallocate    (rc)
    call compress_cov
    if (allocated(xy)) call destruct (xy)

    !-------------------------------------------------------------
    ! count nonzero elements in diagonal and offdiagonal locations
    !-------------------------------------------------------------
    nonzero_bdpe = 0
    nonzero_odpe = 0
    do jb = 1, P_b% m_b
      do ib = 1, P_b% m_b
        if (P_b% b(ib,jb)% nonzero < 0) cycle           ! Skip empty blocks
        if (P_b% b(ib,jb)% pe == dace% pio) then
          if (P_b% b(ib,ib)% pe == P_b% b(jb,jb)% pe) then
            nonzero_bdpe = nonzero_bdpe + P_b% b(ib,jb)% nonzero
          else
            nonzero_odpe = nonzero_odpe + P_b% b(ib,jb)% nonzero
          endif
        endif
      end do
    end do
    nonzero_bdpe = p_sum (nonzero_bdpe)
    nonzero_odpe = p_sum (nonzero_odpe)
    if (dace% lpio) then
      write(6,'()')
      write(6,'("elements on blockdiagonal: ",i10)') nonzero_bdpe
      write(6,'("elements on  off-diagonal: ",i10)') nonzero_odpe
    endif

    !-------------------------------------
    ! assign matrix elements to correct PE
    !-------------------------------------
    if (lmirror .and. ltight) then
      do ib = 1, P_b% m_b
        call p_bcast (P_b% b(ib,:)% repr, P_b% b(ib,ib)% pe)
!NEC$ ivdep
        do jb = 1, P_b% m_b
          if (jb /= ib) then
            select case (P_b% b(ib,jb)% repr)
            case (MIRROR)
              P_b% b(ib,jb)% pe = P_b% b(jb,ib)% pe
            end select
          end if
          if (P_b% b(ib,jb)% pe /= dace% pe) P_b% b(ib,jb)% repr = NOT_INIT
        end do
      end do
      P_b% rowwise = p_and (P_b% rowwise)
    endif
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine set_rc (ib)
    integer ,intent(in) :: ib                  ! box index

      type(t_obs)         ,pointer :: bi       ! observation box pointer
      type(t_spot)        ,pointer :: oi       ! observation type pointer
      integer                      :: is
      integer                      :: io
      logical                      :: le_fi
      real(wp) ,allocatable        :: efi (:)

      bi => obs% o (ib)
      call construct (rc(ib), bi% n_int, bi% n_spt)
      is = 0
      le_fi = lset_efi
      if (le_fi) le_fi = associated (e_fi% s(ib)% x)
      if(le_fi) allocate (efi (bi% n_int))
      do io = 1, bi% n_spot
        oi     => bi% spot(io)
        oi% is = is
        is     = is + oi% n_spt
        if (le_fi) then
          call init_spot_obs (rc(ib), oi, bi, &
            err= efi(oi% i%i+1:oi% i%i+oi% i%n))
        else
          call init_spot_obs (rc(ib), oi, bi)
        endif
      end do
      if (le_fi) then
        where (bi% t_int(1:bi% n_int) /= OBS_DUM) &
          e_fi% s(ib)% x = efi
      endif
    end subroutine set_rc
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine set_Pb
!------------------------------------------------------------------------------
  !--------------------------------------------------------------------
  ! Calculate the vector matrix product Pb * z .
  ! The background error covariance matrix Pb is calculated on the fly.
  ! Further optimization may be obtained by:
  !   a) accounting for the symmetry of the matrix,
  !   b) calculate the product for nonlinear operators only
  !   c) combine the routines Pb_times_z and set_Pb
  !--------------------------------------------------------------------
  subroutine Pb_times_z (Pb_z, z, obs, er, time, pb_apprx)
  type(t_vector)         ,intent(inout) :: Pb_z       ! P_b * z
  type(t_vector)         ,intent(in)    :: z          ! z
  type(t_obs)    ,target ,intent(in)    :: obs (:)    ! observation data type
  real(wp)               ,intent(in)    :: er         ! earth radius
  type (t_time)          ,intent(in)    :: time       ! analysis time
  integer                ,intent(in)    :: pb_apprx   ! 1D-Var flag

    !---------------------------------------------------
    ! local variables: loop indices, pointer and counter
    !---------------------------------------------------
    integer                      :: ib, jb       ! matrix block indices
    integer                      :: io, jo       ! observation type indices
    integer                      :: i, j         ! observation indices
    integer                      :: ii,jj        ! spot indices
    integer                      :: ki,kj        ! spot sizes (levels)
    integer                      :: nonz         ! nonzero elems. from get_corr
    type(t_obs)         ,pointer :: bi, bj       ! observation box pointer
    real(wp)                     :: xi(3), xj(3) ! coordinates
    type(t_spot)        ,pointer :: oi, oj       ! observation type pointer
    integer                      :: is, js       ! spot indices
!   logical                      :: lident       ! int == psas
    type(t_rowcol)  ,pointer     :: rc(:)        ! precalculated for row/col
    type(t_rowcol)  ,pointer     :: lhs,rhs      ! precalculated for row/col
    integer                      :: nb           ! number of boxes
    !-------------------------------------------------------------
    ! local variables: quantities used for the covariance matrices
    !-------------------------------------------------------------
    real(wp) :: d    ! distance in m
    real(wp) :: s    ! secant of the unit sphere
    real(wp) :: dmin ! minimum distance in between spots in a box (diagnostics)
    real(wp) :: dmax ! maximum distance in between spots in a box (diagnostics)
    !------------------------------------
    ! initialize DWD OI correlation model
    !------------------------------------
    Pb_z = 0._wp
    call init_fg_cov (time)
    !----------------------------------------
    ! allocate arrays, precompile information
    !----------------------------------------
    if (test_bl) lpack = .false.
    nb = size(obs)
    allocate (rc(nb))
    do ib=1,nb
      bi => obs (ib)
      call construct (rc(ib), bi% n_int, bi% n_spt)
      is = 0
      do io = 1, bi% n_spot
        oi     => bi% spot(io)
        oi% is = is
        is     = is + oi% n_spt
        call init_spot_obs (rc(ib), oi, bi)
      end do
    end do
    !----------------
    ! loop over boxes
    !----------------
    do jb = 1, nb
      do ib = 1, nb
        !--------------------------------------------
        ! skip if matrix block is located on other PE
        !--------------------------------------------
!!!     if (b% repr >= 0)            cycle
!!!     if (b% pe       /= dace% pe) cycle
        if (obs(ib)% pe /= dace% pe) cycle
        if (pb_apprx <= PB_BLOCKD .or. ib==jb) then
          bi  => obs (ib)
          bj  => obs (jb)
          lhs => rc  (ib)
          rhs => rc  (jb)
!         if (bi% n_lev == 0) call col2lev (bi)
!         if (bj% n_lev == 0) call col2lev (bj)
          !---------------------------
          ! loop over spots in the box
          !---------------------------
          dmin = huge(1._wp)
          dmax = 0._wp
          do jo = 1, bj% n_spot
            oj => bj% spot(jo)
            js = oj% is
            do io = bi% n_spot, 1, -1                  ! for 'L'
!             if (jb==ib .and. io<jo .and. lpack) exit ! representation
!           do io = 1, bi% n_spot                      ! for 'U'
!             if (jb==ib .and. io>jo .and. lpack)    exit
              if (pb_apprx <= PB_1DVAR .and. io/=jo) cycle
              oi => bi% spot(io)
              is = oi% is
              !------------------------------
              ! finally calculate covariances
              ! (univariate so far)
              !------------------------------
              !--------------------------------------------------------
              ! calculate great circle distance
              !   skip if there is no correlation between observations
              !--------------------------------------------------------
              select case (oi% int_type)
              case (ITY_ICOL)
                xi = oi% col% c% x
              case (ITY_MCOLS)
                xi = nearest_x (oj% col% c% x, oi)
              case default
                call finish('Pb_times_z','invalid int_type')
              end select
              select case (oj% int_type)
              case (ITY_ICOL)
                xj = oj% col% c% x
              case (ITY_MCOLS)
                xj = nearest_x (oi% col% c% x, oj)
              case default
                call finish('Pb_times_z','invalid int_type')
              end select
              s = sqrt(sum( (xi - xj)**2) )
              if (s >= x_cut*l_max) cycle
              d = 2._wp *er *asin (.5_wp * s)
              if(jb/=ib .or. io/=jo) then
                dmin=min(d,dmin)
                dmax=max(d,dmax)
              endif
              !----------------
              ! set covariances
              !----------------
!             lident = iand(oi% char, CHR_ID)/=0 &
!                .and. iand(oj% char, CHR_ID)/=0
              nonz = 0
              !---------
              ! loop lhs
              !---------
              ii = 0
              do i = is+1, is+oi% n_spt
                ki = oi%i% n
                if (oi% int_type == ITY_MCOLS)  ki = rc(ib)% spots(i)% n
                !---------
                ! loop rhs
                !---------
                jj = 0
                do j = js+1, js+oj% n_spt
                  kj = oj%i% n
                  if (oj% int_type == ITY_MCOLS)  kj = rc(jb)% spots(j)% n
                  call get_corr (lhs, rhs, i, j, .TRUE., ib==jb, nonz, &
                                 l=Pb_z% s(ib)% x(ii+1:ii+ki),         &
                                 r=   z% s(jb)% x(jj+1:jj+kj)          )
!!!                  call get_corr (lhs, rhs, i, j, .TRUE., ib==jb, nonz, &
!!!                                 coi(ii+1:ii+ki,jj+1:jj+kj)            )
                  jj = jj + kj
                end do
                ii = ii + ki
              end do
            end do
          end do
        endif
      end do
    end do
    !-----------------------
    ! deallocate temporaries
    !-----------------------
    call destruct (rc)
    deallocate    (rc)
  end subroutine Pb_times_z
!==============================================================================
  !----------------------------------------------
  ! setup the observation error covariance matrix
  !----------------------------------------------
  subroutine set_R (obs, y, fvqc, HBH, apprx)
  type(t_obs_set),target ,intent(inout)       :: obs    ! observation data type
  type(t_vector)     ,intent(in)              :: y      ! guess
  real(wp)           ,intent(in)              :: fvqc   ! VarQC flag
  type(t_matrix)     ,intent(in)    ,optional :: HBH    ! matrix HBH
  integer            ,intent(in)    ,optional :: apprx  ! approximation flag
    !---------------------------------------------------
    ! local variables: loop indices, pointer and counter
    !---------------------------------------------------
    integer                      :: ib          ! matrix block indices
    integer                      :: is          ! spot indices
    integer                      :: i,l, k ,j   ! observation indices
    integer                      :: nn0         ! nonzero elements
    integer                      :: nn0b        ! nonzero elements in box
    type(t_matrix_block),pointer :: b           ! matrix block pointer
    type (t_obs_block)           :: ob
!   type(t_vector)               :: j_qci       ! contribution to cost function
    type(t_vector)               :: y_o         ! model - observation
    integer                      :: lappr       ! approximation flag
    !----------------------------------------------
    ! temporary arrays to hold correlation matrices
    !----------------------------------------------
!   real(wp) ,pointer :: coi  (:,:)
    real(wp) ,pointer :: cop  (:,:)
!   real(wp) ,pointer :: coip (:,:)
    real(wp)          :: cop_tmp          ! copy of cop(i,i)
    !--------------------------------------
    ! variational quality control variables
    !--------------------------------------
    type(t_vqc)           :: tapri        ! parameters passed to vqc_corr
    logical  ,allocatable :: lvqc (:)
    integer  ,allocatable :: idx  (:)
    integer  ,allocatable :: ldx  (:)
    real(wp) ,allocatable :: R_qc (:,:)
    real(wp) ,allocatable :: wqc  (:)
    real(wp) ,allocatable :: djqc (:)
    integer  ,allocatable :: obs_stat (:)
    integer               :: ncorr
    real(wp)              :: HBH1
    real(wp)              :: HBH2(2,2)
    real(wp) ,allocatable :: HBHn(:,:)
    type (t_matrix_block) :: HBHb
    integer               :: nu, nc, nv   ! count (un)cor-,wind-observations
    integer               :: nf(0:4)      ! count VQC formulations
    integer               :: frm          ! VQC formulation
    real(wp)              :: sgm          ! VQC threshold
    !---------------------------
    ! process optional arguments
    !---------------------------
    lappr = PB_FULL; if (present(apprx)) lappr = apprx
    obs% vc% J   = 0._wp
    tapri% m_rej = mvqc
    tapri% g_rej = gvqc
!   call construct (j_qci, obs% oi)
    obs% vc% Jo = 0._wp
    call construct (y_o,   obs% oi); y_o = y - obs% obc
    if (test_bl) lpack = .false.
    HBH1 = 0._wp
    HBH2 = 0._wp
    obs% vc% y = y
    !--------------------------------
    ! loop over boxes, diagonals only
    !--------------------------------
    obs% vc% R% qual = BDIAG
    nn0 = 0
    nc  = 0
    nu  = 0
    nv  = 0
    nf  = 0
    do ib = 1, obs% vc% R% m_b
     if (obs% o(ib)% pe == dace% pe) then
      b  => obs% vc% R% b(ib,ib)
      !--------------------------------------------
      ! skip if matrix block is located on other PE
      ! else fill this block
      !--------------------------------------------
      if (b% pe /= dace% pe) cycle
      call obs_block (ob, obs, ib)
      !-------------------------------------------------------
      ! initialise all offdiagonal blocks on this PE with zero
      !-------------------------------------------------------
      obs% vc% R% b(ib,:)% repr = ZERO
      if (.not. associated (ob% o% s_vqc)) then
        allocate (ob% o% s_vqc (ob% o% n_obs))
        ob% o% s_vqc = svqc
      endif
!     if (ob% o% n_lev == 0) call col2lev (ob% o)
      call deallocate (b)
      select case (lappr)
      case (PB_DIAG)
        call allocate_block (b, DIAGONAL)
      case default
        call allocate_block (b, FULL)
        b% full = 0._mp
      end select
      nn0b = 0
      if (present(HBH)) then
        HBHb = HBH% b(ib,ib)
        call full_matrix (HBHb)
      endif
      !--------------------------------------
      ! interprete observation state variable
      !--------------------------------------
      allocate (obs_stat (ob% o% n_obs))
      obs_stat = 0
      if (ob% o% n_obs > 0) then
        where (ob% o% body% use% state <= STAT_NOTACTIVE) obs_stat = -1
!!!     where (ob% o% body% use% state >= STAT_ACTIVE_1)  obs_stat =  1
      endif
      !---------------------------
      ! loop over spots in the box
      !---------------------------
      do is = 1, ob% o% n_spot
        !-----------------------------------------------------------
        ! allocate matrices for interpolated and observational space
        !-----------------------------------------------------------
        allocate (cop (ob% o%spot(is)%o% n,ob% o%spot(is)%o% n))
        select case (lappr)
        case (PB_DIAG)
          cop = 0._wp
        case default
          cop = b% full (ob% o%spot(is)%o% i +                   1 : &
                         ob% o%spot(is)%o% i + ob% o%spot(is)%o% n , &
                         ob% o%spot(is)%o% i +                   1 : &
                         ob% o%spot(is)%o% i + ob% o%spot(is)%o% n )
        end select
        !------------------------------------------------------
        ! loop over variables in observational space
        ! no correlations between different observations so far
        !------------------------------------------------------
        !----------------------------------------------------
        ! add observational errors (diagonal matrices so far)
        !----------------------------------------------------
        if (use_corr) then
          do i=1, ob% o% spot(is)% o% n
!NEC$ ivdep
            do j=ob% R% ia (ob% o%spot(is)% o% i + i), &
                 ob% R% ia (ob% o%spot(is)% o% i + i + 1) - 1
              l = ob% R% ja (j) - ob% o%spot(is)% o% i
              cop(l,i) = cop(l,i) + ob% R% packed(j)
            end do
          end do
        else
          do i=1, ob% o% spot(is)% o% n
            do j=ob% R% ia (ob% o%spot(is)% o% i + i), &
                 ob% R% ia (ob% o%spot(is)% o% i + i + 1) - 1
              l = ob% R% ja (j) - ob% o%spot(is)% o% i
              if (i==l) cop(i,i) = cop(i,i) + ob% R% packed(j)
            end do
          end do
        endif
        !----------------------
        ! quick and dirty varqc
        !----------------------
!       if (fvqc /= 0._wp) then
        if (.TRUE.) then
          allocate (lvqc (ob% o%spot(is)%o% n))
          lvqc = .false.
          !------------------------
          ! correlated observations
          !------------------------
          if (use_corr) then
            !--------------------------------------------------
            ! identify correlated observations within this spot
            !--------------------------------------------------
            ncorr = 0
            select case (ob% o%spot(is)% hd% modtype)
            case (TOVS, TEMP)
              !------------------------------------------------
              ! geopotential height in TEMPs, not used any more
              !------------------------------------------------
              do i = 1, ob% o%spot(is)%o% n
                l = i + ob% o%spot(is)%o% i
                if (ob% R% ia (l+1) /=  ob% R% ia (l)+1) then
                  lvqc(i) = .true.
                  ncorr = ncorr + 1
                endif
              end do
            end select

            if (ncorr > 1) then
              nc    = nc    + ncorr
              nf(1) = nf(1) + ncorr  ! currently always formulation 1
              allocate (ldx(ncorr), idx(ncorr), wqc(ncorr), djqc(ncorr), &
                        R_qc(ncorr,ncorr), HBHn(ncorr,ncorr))
              HBHn = 0
              k = 0
              do  i = 1, ob% o%spot(is)%o% n
                if (lvqc(i)) then
                  l = ob% o%spot(is)%o% i + i
                  k = k + 1
                  ldx (k) = l
                  idx (k) = i
                  tapri% g_rej  = g_rej (ob% o% s_vqc(l))
                endif
              end do
              if (present (HBH)) HBHn(:,:) = HBHb% full(ldx,ldx)

              tapri% m_rej  = min (ncorr, max (ncorr/3,3))

tapri% m_rej  = min (ncorr, max (ncorr/3,3))
tapri% m_rej  = ncorr

tapri% m_rej  =  min (ncorr, mvqc)

              frm            = ob% o% f_vqc (ob% o%spot(is)%o% i+1)
              sgm            = ob% o% s_vqc (ob% o%spot(is)%o% i+1)
              call vqc_corr (y_o%         s(ib)% x(ldx   ),  &! <-  analysis - observation
                             obs% vc% Jo% s(ib)% x(ldx(1)),  &!  -> observational cost function
                             djqc,                           &!  -> gradient of cost function
                             wqc,                            &!  -> weight, probability for correct data
                             cop(idx,idx),                   &! <-  R for correct data
                             R_qc,                           &!  -> modified R, to be used in used in Newton step
                             tapri,                          &! <-  a priori probability flags
                             sgm,                            &! <-  threshold for bad data (tails)
                             .false.,                        &! <-  flag to return inverse R, not R
                             fvqc,                           &! <-  R modification flag
                             obs_stat(ldx),                  &! <-  observ.status flags
                             .false.,                        &! <-  enable printout
!                     dwda =                                 &!  -> d w / d apri
                       HBH = HBHn,                           &! <-  background error
                      form = frm                             )! <-  J_o formulation
              do  i = 1, ncorr
                obs% vc% w%  s(ib)% x(ldx(i)) = wqc (i)
                obs% vc% dJ% s(ib)% x(ldx(i)) = djqc(i)
                do k = 1, ncorr
                  cop                    (idx(i),idx(k)) = R_qc (i,k)
                end do
              end do
              deallocate (ldx, idx, wqc, djqc, R_qc, HBHn)
            endif
          endif
          !--------------------------------------
          ! correlated observations (temperature)
          !--------------------------------------
          !------------------
          ! wind observations
          !------------------
          do  i = 1, ob% o%spot(is)%o% n - 1
            if (lvqc(i)) cycle
            l = i + ob% o%spot(is)%o% i
            if (use_wind .and. ob% o% varno (l)   == VN_U  &
                         .and. ob% o% varno (l+1) == VN_V) then
              frm            = ob% o% f_vqc (l)
              nv             = nv      + 2
              nf(frm)        = nf(frm) + 2
              lvqc(i)        = .true.
              lvqc(i+1)      = .true.
              cop_tmp        = cop(i,i)
              cop(   :i-1,i) = 0._wp
              cop(i+2:   ,i) = 0._wp
              cop(i,   :i-1) = 0._wp
              cop(i,i+2:   ) = 0._wp
              if(present(HBH)) HBH2 = HBHb% full (l:l+1,l:l+1)
              call vqc_wind (y_o%         s(ib)% x(l:l+1),&! <-- analysis - obs
                             obs% vc% Jo% s(ib)% x(l    ),&! <-> cost function
                             obs% vc% dJ% s(ib)% x(l:l+1),&! --> gradient
                             obs% vc% w%  s(ib)% x(l:l+1),&! --> posterior prob.
                             cop_tmp,                     &! <-- R for corr.data
                             cop(i:i+1,i:i+1),            &! <-> (inv.) Hessian
                             ob% o% s_vqc(l),             &! <-- stdev. bad data
                             .false.,                     &! <-- return inv(R)?
                             fvqc,                        &! <-- vqc flag
                             obs_stat(l:l+1),             &! <-- status flags
                         HBH=HBH2,                        &! <-- bg.error
                        form=frm                          )! <-- formulation
            endif
          end do
          !----------------------------------------
          ! remaining observations are uncorrelated
          !----------------------------------------
          do i = 1, ob% o%spot(is)%o% n
            if (lvqc(i)) cycle
            l = i + ob% o%spot(is)%o% i
            frm            = ob% o% f_vqc(l)
            nu             = nu      + 1
            nf(frm)        = nf(frm) + 1
            lvqc(i)        = .true.
            cop_tmp        = cop(i,i)
            cop(   :i-1,i) = 0._wp
            cop(i+1:   ,i) = 0._wp
            cop(i,   :i-1) = 0._wp
            cop(i,i+1:   ) = 0._wp
            if(present(HBH)) HBH1 = HBHb% full (l,l)
            call vqc_uncor (y_o%         s(ib)% x(l), &! <-- analysis - obs
                            obs% vc% Jo% s(ib)% x(l), &! <-> cost function
                            obs% vc% dJ% s(ib)% x(l), &! --> gradient
                            obs% vc% w%  s(ib)% x(l), &! --> posterior prob.
                            cop_tmp,                  &! <-- R for corr.data
                            cop(i,i),                 &! <-- (inv.) Hessian
                            ob% o% s_vqc(l),          &! <-- stdev. bad data
                            .false.,                  &! <-- return inv(R)?
                            fvqc,                     &! <-- vqc flag
                            obs_stat(l),              &! <-- observ.status flag
                        HBH=HBH1,                     &! <-- bg.error
                       form=frm                       )! <-- formulation
          end do
          deallocate (lvqc)
        endif
        select case (lappr)
        case (PB_DIAG)
          do i = 1, ob% o%spot(is)%o% n
            b% packed(ob% o%spot(is)%o% i+i) = cop(i,i)
          end do
        case default
          b% full (ob% o%spot(is)%o% i +                   1 : &
                   ob% o%spot(is)%o% i + ob% o%spot(is)%o% n , &
                   ob% o%spot(is)%o% i +                   1 : &
                   ob% o%spot(is)%o% i + ob% o%spot(is)%o% n ) = cop
        end select
        deallocate (cop)
        obs% vc% Jo_spot% s(ib)% x(is) = &
          sum (obs% vc% Jo% s(ib)% x(ob% o%spot(is)%o% i + 1: &
                                     ob% o%spot(is)%o% i +    &
                                     ob% o%spot(is)%o% n    ) )
      end do
      if (test_r) call test_c (real(b% full,wp), ib)
      !----------------------------------
      ! store matrix in compressed format
      !----------------------------------
      select case (lappr)
      case (PB_DIAG)
      case default
        if (lpack) then
          call pack_matrix (b)
        endif
      end select
      deallocate (obs_stat)
     end if
    end do
    obs% vc% J = sum (obs% vc% Jo)
!   !----------------------------------------------
!   ! no VQC, currently handled by vqc_...-routines
!   !----------------------------------------------
!   if (fvqc == 0._wp) then
!     obs% vc% dJ = inverse_rs (obs% vc% R) * y_o
!     obs% vc% J  = 0.5_wp * sum (obs% vc% dJ * y_o)
!     obs% vc% w  = 0._wp
!   endif
!   call destruct (j_qci)
    call destruct (y_o)
    if (present(HBH)) call deallocate (HBHb)
  end subroutine set_R
!==============================================================================
  subroutine set_H (H, obs, mode)
  type(t_matrix)   ,intent(inout) :: H      ! matrix H to set
  type (t_obs_set) ,intent(in)    :: obs    ! observations
  character        ,intent(in)    :: mode   ! 'f','i'
  !------------------------------------------------------------
  ! Extract the linearised observation operator (matrix H) from
  ! the observation data structure (obs)
  !
  ! mode: 'f' forward matrix
  !       'i' inverse matrix
  !------------------------------------------------------------
    integer                    :: l     ! observation box index
    integer                    :: j     ! observation index
    type (t_spot) ,pointer     :: spt   ! pointer to meta data
    real(wp)      ,allocatable :: Hnew(:,:)     ! temporary
    integer           ,pointer :: itmp (:)      ! temporary
    real(mp)          ,pointer :: mtmp (:)      ! temporary

!   type (t_dec_info) ,pointer :: ri, ci

    select case (mode)
    case default
      call finish ('set_H','invalid mode '//mode)
    case ('f')
      H = obs% l% H
    case ('i')

    H% qual = BDIAG
!   allocate (ri, ci)
    !--------------------------
    ! loop over vector segments
    !--------------------------
    do l = 1, size(obs% o)
      !---------------------------------------
      ! handle only segments on this processor
      !---------------------------------------
      if (obs% o(l)% pe == dace% pe) then
        !-------------------------------------------
        ! allocate matrix-block with sufficient size
        !-------------------------------------------
        H% b(l,l)% repr = CSC
        H% b(l,l)% m = 0
        H% b(l,l)% n = 0
        H% b(l,l)% nonzero = 0
        H% b(l,l)% pe = dace% pe
        do j = 1, obs% o(l)% n_spot
          spt => obs% o(l)% spot(j)
          H% b(l,l)% m = H% b(l,l)% m + spt% i% n
          H% b(l,l)% n = H% b(l,l)% n + spt% o% n
          H% b(l,l)% nonzero = H% b(l,l)% nonzero + spt% o% n &
                                                  * spt% i% n
        end do
        call allocate_block (H% b(l,l), CSC, H% b(l,l)% m, H% b(l,l)% n, &
                             H% b(l,l)% nonzero)
        H% b(l,l)% ia = 1
H% b(l,l)% ja = -9
        !-----------------------------------------------
        ! loop over observations, invert block-diagonals
        !-----------------------------------------------
        do j = 1, obs% o(l)% n_spot
          spt => obs% o(l)% spot(j)
          allocate (Hnew(spt% i% n, spt% o% n))
          Hnew = inverse (sub_block (obs% l% H% b(l,l),spt%o%i,        &
                                     spt%o%n, spt%i%i, spt%i%n, FULL), &
                 max_cond_number=1.e10_wp                              )
          call insert (H% b(l,l), Hnew, spt% i% i, spt% o% i)
          deallocate (Hnew)
        end do
        !---------------------------------------------
        ! reallocate block matrix, do not waste memory
        !---------------------------------------------
        if (H% b(l,l)% nonzero /= H% b(l,l)% ia(H% b(l,l)% n+1) - 1) then
          H% b(l,l)% nonzero = H% b(l,l)% ia(H% b(l,l)% n+1) - 1
          allocate (mtmp(H% b(l,l)% nonzero))
          allocate (itmp(H% b(l,l)% nonzero))
          mtmp = H% b(l,l)% packed (:H% b(l,l)% nonzero)
          itmp = H% b(l,l)% ja     (:H% b(l,l)% nonzero)
          deallocate (H% b(l,l)% packed, H% b(l,l)% ja)
          H% b(l,l)% packed => mtmp
          H% b(l,l)% ja     => itmp
          nullify (mtmp, itmp)
        endif
      end if
    end do
!   deallocate (ri, ci)

    end select

  end subroutine set_H
!==============================================================================
  subroutine constrain_ev (A, obs, evl, fix, verbose, fdiag, scale)
  !---------------------------------------------------------------------------
  ! Determines the smallest eigenvalue of (certain submatrices of) a matrix A.
  ! The eigenvalues and eigenvectors are determined for each submatrix.
  ! Optionally the eigenvalues beyond a given bound are forced to that
  ! bound (and the matrix is reconstructed from the eigenvectors and
  ! modified eigenvalues).
  !---------------------------------------------------------------------------
  type (t_matrix) ,intent(inout)        :: A       ! matrix to test
  type (t_obs)    ,intent(in)           :: obs (:) ! info on observations
  real(wp)        ,intent(in)           :: evl(mt) ! eigenvalue bounds
  integer         ,intent(in)           :: fix     ! 1:fix eigenvalue, 2:abort
  integer         ,intent(in)           :: verbose ! 1..3
  real(wp)        ,intent(in) ,optional :: fdiag   ! apply factor on diagonal
  type (t_vector) ,intent(in) ,optional :: scale   ! scale matrix
  !----------------------------------------------------------------------------
  !
  ! Subroutine arguments:
  !
  !   evl:     Lower bound of the eigenvalues for each test.
  !            For evl<0 the respective test is not performed at all.
  !   fix:     0: do not modify the matrix
  !            1: modify small eigenvalues
  !            2: abort if eigenvalues are too small
  !            -: act only if Cholesky decompositin fails
  !   verbose: 0: be silent
  !            1: write a summary
  !            2: write one line for each defective submatrix
  !            3: write a detailed report for each defective submatrix
  !   scale:   scale the matrix by this vector:
  !            A(i,j) -> 1/sqrt(scale(i)) * A(i,j) * 1/sqrt(scale(j))
  !   fdiag:   if 'fdiag' is present and greater than 1 the expensive tests
  !            described below are not applied. Instead the factor 'fdiag'
  !            is applied to the diagonal unless the matrix is p.d.
  !
  ! Submatrices for each test:
  !
  !   1) diagonal elements.
  !   2) data of the same data type (temperature, humidity)
  !      for each individual report (station).
  !   3) all data of each individual report.
  !   4) test on symmetry (not really a test of p.d.)
  !   5) all data of the blockdiagonal matrix (as used for preconditioning).
  !   6) all data of the same report type (TEMP, SYNOP,...)
  !   7) data of same report type and data type
  !
  ! Reports written
  !
  !   Example summary (verbose>=1):
  !
  !     test             failed  hits    bound   min.ev
  !        1 (diagonal) :     0     0  0.0E+00  0.1E-01
  !        2 (rep./data):     2     2  0.0E+00 -0.7E-15
  !        3 (report)   :     1     2  0.0E+00 -0.2E-13
  !        4 (symmetry) :     0     0  0.0E+00  0.2E-13
  !        5 (block)    :     2     3  0.0E+00 -0.3E-01
  !        6 (rep.-type):     2     3  0.0E+00 -0.3E-01
  !        7 (rep-t/d-t):     1     1  0.0E+00 -0.4E-14
  !
  !     test:   number of the test.
  !     failed: number of tests failed.
  !     hits:   number of eigenvalues below the bound.
  !     bound:  eigenvalue bound.
  !     min.ev: smallest eigenvalue found.
  !
  !   Example brief report (verbose>=2):
  !
  !     test    pe block  dbkz  type   station  count  hits  e-value size file
  !        2     0     1   523  TEMP   16044        8     1 -0.2E-15   15
  !        2     0     1   523  TEMP   16044       16     1 -0.7E-15   15
  !        3     0     1   523  TEMP   16044        0     2 -0.2E-13   50
  !        5     0     1    -1                      0     2 -0.1E-13  513
  !        6     0     1    -1  TEMP                0     2 -0.1E-13  513
  !        7     0     1    -1  TEMP               16     1 -0.4E-14  141
  !        5     0     3    -1                      0     1 -0.3E-01  489
  !        6     0     3    -1  TEMP                0     1 -0.3E-01  489
  !
  !   test:    number of the test.
  !   pe:      number of the processor element (MPI parallel computation).
  !   block:   number of the blockdiagonal matrix.
  !   dbkz:    DWD Datenbankkennzahl (observation subtype).
  !   type:    observation type (TEMP, SYNOP, ..)
  !   station: station name.
  !   count:   in general ID for the data type (temperature, humidity).
  !   hits:    number of eigenvalues below the bound.
  !   e-value: smallest eigenvalue found.
  !   size:    size of the matrix tested:
  !   file:    number of file for detailed report (verbose>=3).
  !----------------------------------------------------------------------------
    target :: obs
    !---------------------------------------------------
    ! local variables: loop indices, pointer and counter
    !---------------------------------------------------
    integer ,parameter           :: mat = 10 ! max. number of attempts
                                             ! to increase diagonal
    integer                      :: ib       ! matrix block indices
    integer                      :: is       ! spot indices
    integer                      :: it       ! test index
    integer                      :: i,j,k,l  ! observation indices
    type(t_matrix_block)         :: b        ! matrix block pointer
    type(t_obs)         ,pointer :: oi       ! observation box pointer
    type(t_spot)        ,pointer :: si       ! observation operator pointer
    integer                      :: fx       ! abs (fix)
    integer                      :: info     ! error return flag of cholesky
    integer                      :: nh(mt)
    integer                      :: nt(mt)
    integer                      :: hits
    real(wp)            ,pointer :: x(:,:)
    integer                      :: n, ns, i0, i1, in
    integer                      :: nc, nspot
    integer                      :: iq,iq1,iqn,nq
    integer                      :: it1,itn,ity
    real(wp)        ,allocatable :: c (:,:)
    integer         ,allocatable :: ic(:)
    logical         ,allocatable :: lxt(:)
    logical         ,allocatable :: lxq(:)
    real(wp)        ,allocatable :: sc(:)
    logical                      :: changed
    real(wp)                     :: evm (mt)
    real(wp)                     :: evmin
    integer                      :: iu       ! unit for asyncrone output
    integer                      :: if = 0   ! number to generate file name
    real(wp)                     :: f        ! local copy of fdiag
    character(len=*),parameter   :: fh = &
"('  test    pe block  dbkz  type   station  count  hits  e-value size file')"
    character(len=*),parameter   :: fb = "(4i6,1x,a6,a10,2i6,e9.1,i5,i10)"
    real(wp) ,parameter          :: eps = 1.e-11
    !-------------
    ! write report
    !-------------
    if (all(evl<0._wp)) return
    iu = 0; if (dace% npe==1) iu = 6
    call flush_buf
    f = 0._wp; if (present(fdiag)) f = fdiag
    if = dace% pe * 1000
    if (verbose>0 .and. dace% lpio) then
!     write (6,'(a)') repeat ('-',79)
      write (6,'()')
      write (6,'(a)') '  Check matrix for small eigenvalues'
      write (6,'()')
      write (6,'(a,a )')   '    matrix = ',trim(a% name)
      write (6,'(a,l1)')   '   scaling = ',present(scale)
      write (6,'(a,i2)')   '       fix =' ,fix
      write (6,'(a,i2)')   '   verbose =' ,verbose
      write (6,'(a,f6.3)') '     fdiag =' ,f
      write (6,'()')
      if (verbose>1) write (6,fh)
    endif
    fx = abs(fix)
    !--------------------------------------
    ! loop over boxes, diagonal blocks only
    !--------------------------------------
    evm = huge(evm)
    nh  = 0
    nt  = 0
    do ib = 1, A% m_b
      !--------------------------------------------
      ! skip if matrix block is located on other PE
      ! else process this block
      !--------------------------------------------
      if (A% b(ib,ib)% pe /= dace% pe) cycle
      oi  => obs (ib)
      n   =  A% b(ib,ib)% m
      !-------------------------------------------------------------
      ! fast check on negative eigenvalues by Cholesky decomposition
      !-------------------------------------------------------------
      info    = -1
      if (fix < 0) then
        b = A% b(ib,ib)
        select case (b% repr)
        case (DIAGONAL)
          if (all (b% packed >= 0._mp)) info = 0
        case default
          call cholesky (b, info=info)
          call deallocate (b)
        end select
      endif
      if (info == 0) cycle
      !-------------------------------------------------------
      ! make a copy of the block in full representation, scale
      !-------------------------------------------------------
      allocate (sc(n))
      if (present (scale)) then
        where (scale% s(ib)% x > 0._wp)
          sc = 1._wp / sqrt(scale% s(ib)% x)
        elsewhere
          sc = 1._wp
        endwhere
      else
        sc = 1._wp
      endif
      b = A% b(ib,ib)
      if (b% repr /= FULL) call full_matrix (b)
      allocate (x(b%m,b%n))
      x = b% full
      changed = .false.
      do i = 1, n
        x (:,i) = x (:,i) * sc * sc(i)
      end do
      !-------------------------------------
      ! if (fiag>1) simply increase diagonal
      !-------------------------------------
      if (f>1._wp) then
        k = 0
        do
          k = k + 1
          if (info>0) then
            if (fx==1) changed = .true.
            if (k==1)  nt(4) = nt(4) + 1
            nh(4)            = nh(4) + 1
            do i = 1, n
              x(i,i) = x(i,i) * f
            end do
          endif
          b% full = x
          call cholesky (b, info=info)
          call allocate_block (b, FULL, dealloc=.true.)
          if (info==0) exit
          if (k > mat) call finish ('constrain_ev',                         &
            'number of attempts to increase diagonal exhausted in '//a% name)
        end do
      else
        !---------------------------
        ! loop over spots in the box
        !---------------------------
        do is = 1, oi% n_spot
          si => oi% spot(is)
          ns = si% o% n
          i1 = si% o% i + 1
          in = si% o% i + ns
          !======================================
          ! 1st test: check for diagonal elements
          !======================================
          it = 1
          if (evl(it) >= -1._wp) then
            do i = i1, in
              evm(it) = min(evm(it), x (i,i))
              if (x (i,i) < evl(it) - eps) then
                nh(it) = nh(it) + 1
                nt(it) = nt(it) + 1
                if (verbose>1) then
                  call nextline
                  write(oline(iol),fb) it,dace% pe,ib,si%hd%dbkz,      &
                    rct_name(si%hd%modtype), si%statid,i-i1+1,1,x(i,i),1
                endif
                if (fx==1) then
                  x (i,i) = evl(it)
                  changed = .true.
                endif
              endif
            end do
          endif
          !========================================================
          ! 2nd test: check for submatrix of same observed quantity
          !========================================================
          it = 2
          if (evl(it) >= -1._wp) then
            iq1 = minval (oi% varno(i1:in))
            iqn = maxval (oi% varno(i1:in))
            do i = iq1, iqn
              nq = count (oi% varno(i1:in) == i)
              if (nq == 0 ) cycle
!             if (nq <= 1 ) cycle
!             if (nq == ns) cycle
              hits = 0
              allocate (ic(nq))
              l = 0
              do k = i1, in
                if (oi% varno(k) == i) then
                  l = l + 1
                  ic(l) = k
                endif
              end do
              if (fx==1) then
                allocate (c (nq,nq))
                call check_rs (x(ic,ic), min_ev=evl(it), y=c,&
                                    hits=hits, evmin=evmin)
              else
                call check_rs (x(ic,ic), min_ev=evl(it), &
                               hits=hits, evmin=evmin)
              endif
              if (evmin < evl(it) - eps) then
                if (verbose > 2) then
                  if = if + 1
                  call test_c (x(ic,ic), if, iu)
                  call nextline
                  write(oline(iol),fb) it, dace% pe, ib, si%hd% dbkz,  &
                                       rct_name(si%hd% modtype),       &
                                       si%statid, i, hits, evmin, nq, if
                else if (verbose > 1)  then
                  call nextline
                  write(oline(iol),fb) it, dace% pe, ib, si%hd% dbkz, &
                                       rct_name(si%hd% modtype),      &
                                       si%statid, i, hits, evmin, nq
                endif
              endif
              if (fx==1) then
                if (hits>0) then
                  x(ic,ic) = c(:,:)
                  changed = .true.
                endif
                deallocate (c)
              endif
              deallocate (ic)
              if (hits > 0) then
                nh(it) = nh(it) + hits
                nt(it) = nt(it) + 1
              endif
              evm(it) = min (evm(it), evmin)
            end do
          endif
          !===========================================================
          ! 3rd test: check for submatrix of same observation operator
          !===========================================================
          it = 3
          if (evl(it) >= -1._wp) then
            hits = 0
            if (fx==1) then
              allocate (c(ns,ns))
              call check_rs (x(i1:in,i1:in), min_ev=evl(it), y=c, &
                             hits=hits, evmin=evmin)
            else
              call check_rs &
                (x(i1:in,i1:in), min_ev=evl(it), hits=hits, evmin=evmin)
            endif
            if (evmin < evl(it) - eps) then
              if (verbose > 2) then
                if = if + 1
                call test_c (x(i1:in,i1:in), if, iu)
                call nextline
                write(oline(iol),fb) it, dace% pe, ib, si%hd% dbkz,         &
                  rct_name(si%hd% modtype), si%statid, 0, hits, evmin, ns, if
              else if (verbose > 1)  then
                call nextline
                write(oline(iol),fb) it, dace% pe, ib, si%hd% dbkz,     &
                  rct_name(si%hd% modtype), si%statid, 0, hits, evmin, ns
              endif
            endif
            if (fx==1) then
              if (hits>0) then
                x(i1:in,i1:in) = c(:,:)
                changed = .true.
              endif
              deallocate (c)
            endif
            if (hits > 0) then
              nh(it) = nh(it) + hits
              nt(it) = nt(it) + 1
            endif
            evm(it) = min (evm(it), evmin)
          endif
        end do
        !=============================
        ! 4th test: check for symmetry
        !=============================
        it = 4
        if (evl(it) >= -1._wp) then
          hits = 0
          do i = 1, n
            do j = i+1, n
              evmin   = - abs (x(i,j)-x(j,i))
              evm(it) = min (evm(it), evmin)
              if (evmin < evl(it) - eps) then
                hits = hits + 1
              endif
            end do
          end do
          if (hits > 0) then
            nh(it) = nh(it) + hits
            nt(it) = nt(it) + 1
          endif
        endif
        !============================================
        ! 5th test: check for submatrix of same block
        !============================================
        it = 5
        if (evl(it) >= 0._wp) then
          hits = 0
          if (fx==1) then
            allocate (c(n,n))
            call check_rs (x, min_ev=evl(it), y=c, hits=hits, evmin=evmin)
          else
            call check_rs (x, min_ev=evl(it), hits=hits, evmin=evmin)
          endif
          if (evmin < evl(it) - eps) then
            if (verbose > 2) then
              if = if + 1
              call test_c (x, if, iu)
              call nextline
              write(oline(iol),fb) it, dace% pe,ib,-1,' ',' ', 0, hits, evmin, n, if
            else if (verbose > 1)  then
              call nextline
              write(oline(iol),fb) it, dace% pe,ib,-1,' ',' ', 0, hits, evmin, n
            endif
          endif
          if (fx==1) then
            if (hits>0) then
              x = c
              changed = .true.
            endif
            deallocate (c)
          endif
          if (hits > 0) then
            nh(it) = nh(it) + hits
            nt(it) = nt(it) + 1
          endif
          evm(it) = min (evm(it), evmin)
        endif
        !=============================================
        ! 6th test: check for all reports of same type
        !=============================================
        it = 6
        if (evl(it) >= 0._wp .or. evl(it+1) >= 0._wp) then
          !------------------------------------------------------
          ! determine report types in use, loop over report types
          !------------------------------------------------------
          nspot = oi% n_spot
          it1 = minval (oi% spot(1:nspot)% hd% modtype)
          itn = maxval (oi% spot(1:nspot)% hd% modtype)
          allocate (lxt(it1:itn)) !mask for report types
          lxt = .false.
          do is = 1,nspot
            lxt(oi% spot(is)% hd% modtype) = .true. !mask for report types
          end do
          do ity = it1, itn
            if (.not. lxt(ity)) cycle
            !------------
            ! test 6 only
            !------------
            if (evl(it) >= 0._wp)  then
              !--------------------------------------------------
              ! determine size of submatrix, allocate temporaries
              !--------------------------------------------------
              nc = sum(     oi% spot(1:nspot)% o% n,      &
                       mask=oi% spot(1:nspot)% hd% modtype == ity)
              allocate (ic(nc))
              !------------------------------
              ! set index array for submatrix
              !------------------------------
              i = 0
              do is=1,nspot
                if (oi% spot(is)% hd% modtype == ity) then
                  si => oi% spot(is)
                  ns = si% o% n
                  i0 = si% o% i
!                 ic(i+1:i+ns) = i0 + (/(k,k=1,ns)/) !+++ fails with gfortran 4.9.1
                  do k = 1, ns
                    ic(i+k) = i0 + k
                  end do
                  i=i+ns
                endif
              end do
              !---------------
              ! perform test 6
              !---------------
              hits = 0
              if (fx==1) then
                allocate (c(nc,nc))
                call check_rs(x(ic,ic),min_ev=evl(it),y=c,hits=hits,evmin=evmin)
              else
                call check_rs (x(ic,ic), min_ev=evl(it), hits=hits, evmin=evmin)
              endif
              if (hits>0) then
                if (verbose > 2) then
                  if = if + 1
                  call test_c (x(ic,ic), if, iu)
                  call nextline
                  write(oline(iol),fb) it, dace% pe, ib, -1, rct_name(ity),' ', 0,&
                                       hits, evmin, nc, if
                else if (verbose > 1)  then
                  call nextline
                  write(oline(iol),fb) it, dace% pe, ib, -1, rct_name(ity),' ', 0,&
                                       hits, evmin, nc
                endif
              endif
              if (fx==1) then
                if (hits>0) then
                  x(ic,ic) = c
                  changed = .true.
                endif
                deallocate (c)
              endif
              deallocate (ic)
              if (hits > 0) then
                nh(it) = nh(it) + hits
                nt(it) = nt(it) + 1
              endif
              evm(it) = min (evm(it), evmin)
            endif
            !========================================================
            ! 7th test: check for submatrix of same observed quantity
            !========================================================
            it = 7
            if (evl(it) >= 0._wp) then
              !----------------------------
              ! determine data types in use
              !----------------------------
              iq1 =  huge(iq1)
              iqn = -huge(iqn)
              do is=1,nspot
                if (oi% spot(is)% hd% modtype == ity) then
                  si  => oi% spot(is)
                  ns  = si% o% n
                  i1  = si% o% i + 1
                  in  = si% o% i + ns
                  iq1 = min (minval(oi% varno(i1:in)), iq1)
                  iqn = max (maxval(oi% varno(i1:in)), iqn)
                endif
              end do
              allocate (lxq(iq1:iqn))
              lxq = .false.
              do is=1,nspot
                if (oi% spot(is)% hd% modtype == ity) then
                  si  => oi% spot(is)
                  ns  = si% o% n
                  i1  = si% o% i + 1
                  in  = si% o% i + ns
                  do i= i1, in
                    lxq(oi% varno(i)) = .true.
                  end do
                endif
              end do
              !--------------------------------------------------
              ! determine size of submatrix, allocate temporaries
              !--------------------------------------------------
              do iq = iq1, iqn
                if (.not.lxq(iq)) cycle
                nc = 0
                do is=1,nspot
                  if (oi% spot(is)% hd% modtype /= ity) cycle
                  si  => oi% spot(is)
                  ns  = si% o% n
                  i1  = si% o% i + 1
                  in  = si% o% i + ns
                  nc = nc + count (oi% varno(i1:in) == iq)
                end do
                allocate (ic(nc))
                !------------------------------
                ! set index array for submatrix
                !------------------------------
                i = 0
                do is=1,nspot
                  if (oi% spot(is)% hd% modtype /= ity) cycle
                  si  => oi% spot(is)
                  ns  = si% o% n
                  i1  = si% o% i + 1
                  in  = si% o% i + ns
                  do k = i1, in
                    if (oi% varno(k) /= iq) cycle
                    i = i + 1
                    ic(i) = k
                  end do
                end do
                !---------------
                ! perform test 6
                !---------------
                hits = 0
                if (fx==1) then
                  allocate (c(nc,nc))
                  call check_rs (x(ic,ic),&
                    min_ev=evl(it),y=c,hits=hits,evmin=evmin)
                else
                  call check_rs (x(ic,ic),min_ev=evl(it),hits=hits, evmin=evmin)
                endif
                if (hits>0) then
                  if (verbose > 2) then
                    if = if + 1
                    call test_c (x(ic,ic), if, iu)
                    call nextline
                    write(oline(iol),fb) it, dace% pe, ib,-1, rct_name(ity),' ',iq,&
                                         hits, evmin, nc, if
                  else if (verbose > 1)  then
                    call nextline
                    write(oline(iol),fb) it, dace% pe, ib,-1, rct_name(ity),' ',iq,&
                                         hits, evmin, nc
                  endif
                endif
                if (fx==1) then
                  if (hits>0) then
                    x(ic,ic) = c
                    changed = .true.
                  endif
                  deallocate (c)
                endif
                deallocate (ic)
                if (hits > 0) then
                  nh(it) = nh(it) + hits
                  nt(it) = nt(it) + 1
                endif
                evm(it) = min (evm(it), evmin)
              end do ! iq
              deallocate (lxq)
            end if ! evl(7) >= 0._wp
          end do   ! type
          deallocate (lxt)
        endif
      endif
      !-----------------------------
      ! write back matrix if changed
      !-----------------------------
      if (changed) then
        sc = 1._wp / sc
        do i = 1, n
          x (:,i) = x (:,i) * sc * sc(i)
        end do
        b% full = x
        select case (A% b(ib,ib)% repr)
        case (PACKED)
          call pack_matrix (b)
          A% b(ib,ib) = b
        case (FULL)
          A% b(ib,ib) = b
        case default
          call finish ('constrain_ev(alloc)',&
            'representation not implemented: '//srep(A% b(ib,ib)% repr))
        end select
      endif
      deallocate (x)
      call deallocate (b)
      deallocate (sc)
    end do            ! boxes
    !-------------
    ! write report
    !-------------
    call flush_buf
    nh  = p_sum (nh)
    nt  = p_sum (nt)
    evm = p_min (evm)
    if (verbose > 0) then
      if (dace% lpio) then
       write (6,'()')
       write (6,'(a)')     '  test             failed  hits    bound   min.ev'
       it = 1
       if(evl(it)>=-1)&
         write (6,'(5x,a,2i6,2e9.1)')'1 (diagonal) :',nt(it),nh(it),evl(it),evm(it)
       it = 2
       if(evl(it)>=-1)&
         write (6,'(5x,a,2i6,2e9.1)')'2 (rep./data):',nt(it),nh(it),evl(it),evm(it)
       it = 3
       if(evl(it)>=-1)&
         write (6,'(5x,a,2i6,2e9.1)')'3 (report)   :',nt(it),nh(it),evl(it),evm(it)
       it = 4
       if(evl(it)>=-1)&
         write (6,'(5x,a,2i6,2e9.1)')'4 (symmetry) :',nt(it),nh(it),evl(it),evm(it)
       it = 5
       if(evl(it)>=-1)&
         write (6,'(5x,a,2i6,2e9.1)')'5 (block)    :',nt(it),nh(it),evl(it),evm(it)
       it = 6
       if(evl(it)>=-1)&
         write (6,'(5x,a,2i6,2e9.1)')'6 (rep.-type):',nt(it),nh(it),evl(it),evm(it)
       it = 7
       if(evl(it)>=-1)&
         write (6,'(5x,a,2i6,2e9.1)')'7 (rep-t/d-t):',nt(it),nh(it),evl(it),evm(it)
      endif
    endif
    !------
    ! abort
    !------
    call p_barrier
    if (fx > 1 .and. any (evm < evl - eps)) then
      call finish ('constrain_ev','eigenvalues too small in '//a% name)
    endif
  end subroutine constrain_ev
!==============================================================================
  subroutine test_c (c, ib, unit)
  !--------------------------------------------------
  ! test correlation matrix for positive definiteness
  !--------------------------------------------------
  real(wp) ,intent(in)           :: c(:,:) ! matrix to check
  integer  ,intent(in)           :: ib     ! report file number
  integer  ,intent(in) ,optional :: unit   ! unit for brief report (<0:no)

    real(wp) ,allocatable :: a  (:)
    real(wp) ,allocatable :: cc (:,:)
    real(wp) ,allocatable :: v  (:,:)
    integer  ,allocatable :: ix (:)
    real(wp) ,allocatable :: y  (:)
    real(wp) ,allocatable :: z  (:)
    integer               :: n
    integer               :: i,j
    real(wp)              :: us, mad, mid
    character(len=20)     :: file
    real(wp)              :: r
    integer               :: iu

    iu = 6 ;if(present(unit)) iu = unit
    !----------------------
    ! allocate local arrays
    !----------------------
    n = size (c,1)
    allocate (v(n,n), cc(n,n), a(n), ix(n), y(n), z(n))
    cc = c
    !-----------------
    ! normalize matrix
    !-----------------
    do j=1,n
      a(j) = cc (j,j)
    end do
    mad = maxval (a)
    mid = minval (a)
    where (a/=0._wp) a = 1._wp / sqrt(abs(a))
    do j=1,n
      cc(j,:) = cc(j,:) * a(j) * a
    end do
    !------------------
    ! test for symmetry
    !------------------
    us = 0._wp
    do j=1,n
      do i=1,j-1
        us = max (us, abs(cc(i,j)-cc(j,i)))
cc(i,j) = 0.5_wp*(cc(i,j)+cc(j,i))
cc(j,i) = cc(i,j)
      end do
    end do

    call jacobi (cc, a, v)
    ix = index (a)

    y = matmul (cc,v(:,1))
    z = matmul (cc,v(:,n))

    !---------------------------------
    ! write brief summary (to stdout?)
    !---------------------------------
    write(file,'(a12,i8.8)') 'eigenvalues_',ib
    if (iu>=0) then
      write(iu,*)
      write(iu,*) 'test_c: min, max diag.elem. =',mid,mad
      write(iu,*) 'test_c: min, max eigenvalue =',a(ix(1)),a(ix(n))
      write(iu,*) 'test_c: min, max eigenvalue =',sum(y*v(:,1)), sum(z*v(:,n))
      write(iu,*) 'test_c: eigenvalues <=> 0   =',count(a<0),count(a==0),&
                                                  count(a>0)
      write(iu,*) 'test_c: max(abs(cc(i,j)-cc(j,i))) =',us
      write(iu,*) 'test_c: filename =',file
      write(iu,*)
    endif

    !--------------------------
    ! write eigenvalues to file
    !--------------------------
    iu = get_unit_number()
    open (iu,file=file)
    write(iu,*)'# min, max diag.elem. =',mid,mad
    write(iu,*)'# min, max eigenvalue =',a(ix(1)),a(ix(n))
    write(iu,*)'# min, max eigenvalue =',sum(y*v(:,1)), sum(z*v(:,n))
    write(iu,*)'# eigenvalues <=> 0   =',count(a<0),count(a==0),count(a>0)
    write(iu,*)'# max(abs(cc(i,j)-cc(j,i))) =',us
    do i=1,n
      write(iu,*) i, a(ix(i))
    end do
    close(iu)

    !---------------------------
    ! write eigenvectors to file
    !---------------------------
    write(file,'(a12,i8.8)') 'eigenvector_',ib
    open (iu,file=file)
    if (n>10) then
      do i=1,n
        write(iu,*) i,real(v(i,1:5)), real(v(i,n-4:n))
      end do
    else
      do i=1,n
        write(iu,*) i,real(v(i,:))
      end do
    endif
    close(iu)

    !------------------------
    ! write columns ? to file
    !------------------------
    do i=1,min(10,size(v,2))
      call random_number (r)
      j = int(n*r)+1
      a = 0._wp
      a(j) = 1._wp
      v(:,i) = matmul (cc,a)
    end do
    file=''
    write(file,'(a8,i8.8)') 'columns_',ib
    open (iu,file=file)
    do i=1,n
      write(iu,*) i,v(i,1:min(10,n))
    end do
    close(iu)
    call return_unit_number(iu)

  end subroutine test_c
!==============================================================================
  subroutine EnKF_precond_0 (x, obs, i_ensb)
  !----------------------------------------------------------------------
  ! prepare preconditioning matrix (B_NMC + B_EnKF for VarEnKF)
  ! 1st step: interpolate ensemble deviations to location of observations
  !----------------------------------------------------------------------
  type(t_vector)  ,allocatable :: x(:)   ! ensemble deviations in obs. space
  type(t_obs_set) ,intent(in)  :: obs    ! observation meta data
  integer         ,intent(in)  :: i_ensb ! 0: do nothing
                                         ! 1: result in interpolation space
                                         ! 2: result on observation space
    integer        :: k  ! ensemble member index
    type(t_vector) :: z  ! temporary

    select case (i_ensb)
    case (1)
    !---------------------------------------------------
    ! result: background ensemble in interpolation space
    !         (input to observation operators)
    !---------------------------------------------------
      !-----------------------------
      ! allocate ensemble deviations
      !-----------------------------
      allocate       (x (Benkf% n_ens))
      call construct (x, obs% ii, 'Xi')
      !--------------------------------
      ! interpolate ensemble deviations
      !--------------------------------
      do k=1, Benkf% n_ens
        call apply_Hi_ens (Benkf% X (k,Benkf% l_set), x(k))
      end do
    case (2)
    !----------------------------------------------------
    ! result: background ensemble in observation space
    !         (output of (linear)  observation operators)
    !----------------------------------------------------
      allocate       (x (Benkf% n_ens))
      call construct (x, obs% oi, 'Yi')
      call construct (z, obs% ii, 'Xi')
      !--------------------------------
      ! interpolate ensemble deviations
      !--------------------------------
      do k=1, Benkf% n_ens
        call apply_Hi_ens (Benkf% X (k,Benkf% l_set), z)
        x(k) = obs% l% H * z
      end do
      call destruct (z)
    end select

  end subroutine EnKF_precond_0

!===================================================================
  subroutine EnKF_precond_1 (coi, x, obs, oi, oj, ci, cj, ib, jb, ii, jj)
  !------------------------------------------------------------
  ! prepare preconditioning matrix (B_NMC + B_EnKF for VarEnKF)
  ! in 'interpolation space'
  ! 2nd step:
  !  for each pair of reports
  !    derive B_EnKF, localize on B, sum up B_NMC + B_EnKF
  !------------------------------------------------------------
  real(mp)        ,intent(inout) :: coi (:,:) ! B matrix block, interpol. space
  type(t_vector)  ,intent(in)    :: x (:)     ! ensemble deviations
  type(t_obs_set) ,intent(in)    :: obs       ! observation meta data
  type(t_spot)    ,intent(in)    :: oi        ! lhs report meta data
  type(t_spot)    ,intent(in)    :: oj        ! rhs report meta data
  type(t_coord)   ,intent(in)    :: ci        ! lhs coordinate
  type(t_coord)   ,intent(in)    :: cj        ! rhs coordinate
  integer         ,intent(in)    :: ib, jb    ! matrix block indices
  integer         ,intent(in)    :: ii, jj    ! start indices

    integer        :: k                    ! ensemble member index
    integer        :: ni, nj               ! size of report
    integer        :: i, j, i2, j2         ! indizes
    integer        :: iii, jjj             ! indizes
    real(wp)       :: dx, lx, wx           ! x distance, length scale, weight
    real(wp)       :: w  (size(coi,1))     !   weight
    real(wp)       :: zi (size(coi,1))     ! z coordinate (log(p))
    real(wp)       :: zj (size(coi,2))     ! z coordinate (log(p))
    real(wp)       :: li, lj               ! z coordinate
    real(wp)       :: b (size(coi,1), size(coi,2)) ! b_enkf
    real(wp)       :: c (size(coi,1), size(coi,2)) ! localisation matrix
    real(wp)       :: mi(size(coi,1))              ! mask
    real(wp)       :: mj(size(coi,2))              ! mask
    real(wp)       :: m                            ! mask
    logical        :: lsym                 ! symmetric matrix block

    if (ib /= jb)               call finish ('EnKF_precond_1','ib /= jb')
    if (obs% o(ib)% pe /= dace% pe)                 &
      call finish ('EnKF_precond_1','pe /= dace% pe')

    !-----------------------------------------
    ! set horizontal localisation length scale
    !-----------------------------------------
    lx = sqrt(10._wp/3._wp) * Benkf% hloc% r_loc / rearth
    !-------------------------------
    ! set up ensemble B matrix block
    !-------------------------------
    ni   = size (coi,1)
    nj   = size (coi,2)
    iii  = oi% i% i + ii
    jjj  = oj% i% i + jj
    lsym = ib==jb .and. iii==jjj

    !-------------
    ! localisation
    !-------------
    mi = sqrt (w_ens_b)
    do i = 1,  ni
      i2 = i + iii
      select case (obs% o(ib)% t_int (i2))
      case (OBS_DUM, OBS_CTR)  ! sink variables or VarBC
        mi(i) = 0._wp
      case default
        zi(i) = Benkf% vloc% fp * log(Benkf% vloc% r_0 - obs% o(ib)% lev(i2))
      end select
    end do

    mj = sqrt (w_ens_b)
    do j = 1,  nj
      j2 = j + jjj
      select case (obs% o(jb)% t_int (j2))
      case (OBS_DUM, OBS_CTR)  ! sink variables or VarBC
        mj(j) = 0._wp
      case default
        zj(j) = Benkf% vloc% fp * log(Benkf% vloc% r_0 - obs% o(jb)% lev(j2))
      end select
    end do

    c = 1._wp
    dx = sqrt (sum ( (ci% x - cj% x ) ** 2 ))
    wx = gaspari_cohn (dx, lx)
    lj = huge(1._wp)
    do j = 1,  nj
      j2 = j + jjj
      select case (obs% o(jb)% t_int (j2))
      case (OBS_DUM, OBS_CTR)  ! sink variables or VarBC
      case default
        li = huge(1._wp)
        do i   = 1,  ni
          i2   = i + iii
          select case (obs% o(ib)% t_int (i2))
          case (OBS_DUM, OBS_CTR)   ! sink variables or VarBC
          case default
            if (lsym .and. i < j) then
              c (i,j) = c (j,i)
            else
              if   (obs% o(jb)% lev(j2) /= lj) then
                if (obs% o(ib)% lev(i2) /= li) then
                  li = obs% o(ib)% lev(i2)
                  w (i) = wx * gaspari_cohn (abs (zi(i) - zj(j)), sqrt(2._wp))
                else
                  w (i) = w (i-1)
                endif
              endif
              c (i,j) = w (i)
            endif
          end select
        end do  ! i
        lj = obs% o(jb)% lev(j2)
      end select
    end do      ! j

    !------------------------
    ! Original outer product:
    !------------------------
!   b = 0._wp
!   do k=1, Benkf% n_ens
!     b = b + (x(k)% s(ib)% x(iii+1:iii+ni) .o. x(k)% s(jb)% x(jjj+1:jjj+nj))
!   end do
    !-------------------------------------------------------
    ! Manual inlining and expansion for better optimization:
    !-------------------------------------------------------
    b = 0._wp
    do k=1, Benkf% n_ens
      do j = 1, nj
         b(:,j) = b(:,j) + x(k)% s(ib)% x(iii+1:iii+ni) * x(k)% s(jb)% x(jjj+j)
      end do
    end do
    !------------------------------------
    ! Schur product, merge with NMC B,
    ! don't touch sink variables or VarBC
    !------------------------------------
    do   j = 1,  nj
      do i = 1,  ni
        m = mi(i) * mj(j)
        coi(i,j) = (1._mp - m) * coi(i,j)        &
                 +          m  * b  (i,j) * c(i,j)
      end do
    end do

  end subroutine EnKF_precond_1

!==============================================================================

  subroutine EnKF_precond_2 (Pb, y, obs, oi, oj, ib, jb)
  !-------------------------------------------------------------
  ! prepare preconditioning matrix (B_NMC + B_EnKF for VarEnKF)
  ! in 'observation space'
  ! 2nd step:
  !  for each pair of reports
  !    derive B_EnKF, localise on HBHt, sum up B_NMC + B_EnKF
  ! Note: localisation on HBHt is faster for observation types
  ! with few output variables (compared to the number of inputs)
  ! but is inverior to localisation on B (EnKF_precond_1)
  !-------------------------------------------------------------
  type(t_matrix_block),intent(inout) :: Pb     ! block of P_b
  type(t_vector)      ,intent(in)    :: y (:)  ! ensemble deviations
  type(t_obs_set)     ,intent(in)    :: obs    ! observation meta data
  type(t_spot)        ,intent(in)    :: oi     ! lhs report meta data
  type(t_spot)        ,intent(in)    :: oj     ! rhd report meta data
  integer             ,intent(in)    :: ib, jb ! matrix block indices

    integer        :: k                    ! ensemble member index
    integer        :: ni, nj               ! size of report
    integer        :: i, j                 ! indices
    integer        :: iii, jjj             ! indices
    real(wp)       :: dx, lx, wx           ! x distance, length scale, weight
    real(wp)       :: dz, lz, wz           ! z distance, length scale, weight
    real(wp)       :: zi                   ! z coordinated (log(p))
    real(wp)       :: zj(oj%o%i+1:oj%o%i+oj%o%n)
    real(wp)       :: b
    !-------------------------------
    ! set up ensemble B matrix block
    !-------------------------------
    ni  = oi% o% n
    nj  = oj% o% n
    iii = oi% o% i
    jjj = oj% o% i
    select case (Pb% repr)
    case (DIAGONAL)
    !---------------------------------------
    ! set only diagonal elements (variances)
    ! for the first guess check
    ! without localisation
    !---------------------------------------
      do i = iii+1, iii+ni
        b = 0._wp
        do k = 1, Benkf% n_ens
          b = b + y(k)% s(ib)% x(i) ** 2
        end do
        Pb% packed (i) = (1._wp - w_ens_b) * Pb% packed (i) &
                       +          w_ens_b  * b
      end do
    case (FULL)
    !-------------------------------------
    ! set matrix block (with localisation)
    !-------------------------------------
      lx = sqrt(10._wp/3._wp) * Benkf% hloc% r_loc / rearth
      dx = sqrt (sum ( (oi% col% c% x - oj% col% c% x ) ** 2 ))
      wx = gaspari_cohn (dx, lx)
      do j = jjj+1, jjj+nj
        zj(j) = log(Benkf% vloc% r_0 - log(real (obs% o(jb)% body(j)% plev,wp)))
      end do
      do i = iii+1, iii+ni
        zi =    log(Benkf% vloc% r_0 - log(real (obs% o(ib)% body(i)% plev,wp)))
      do j = jjj+1, jjj+nj
        dz = Benkf% vloc% fp * abs (zi - zj(j))
        wz = gaspari_cohn (dz, sqrt(2._wp))
        b = 0._wp
        do k = 1, Benkf% n_ens
          b = b + y(k)% s(ib)% x(i) * y(k)% s(jb)% x(j)
        end do
        Pb% full (i,j) = (1._wp - w_ens_b) * Pb% full (i,j) &
                       +          w_ens_b  * b * wx * wz
      end do
      end do

    case default
      write(0,*)   'EnKF_precond_2:  unknown Pb% repr =',Pb% repr
      call finish ('EnKF_precond_2','unknown Pb% repr')
    end select

  end subroutine EnKF_precond_2

!==============================================================================
end module mo_set_matrix
