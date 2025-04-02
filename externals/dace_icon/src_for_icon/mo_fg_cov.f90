!
!+ explicit (PSAS) background error covariance formulation
!
MODULE mo_fg_cov
!
! Description:
!   explicit (PSAS) background error covariance formulation
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
! V1_5         2009/05/25 Harald Anlauf
!  init_spot: optimizations for SX-9
!  namelist /PSCMODEL/: emin_q allows setting global lower bound on RH error
! V1_6         2009/06/10 Harald Anlauf
!  get_corr, get_corr_opt: cleanup
! V1_8         2009/12/09 Harald Anlauf
!  optimizations for SX-9
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_15        2011/12/06 Harald Anlauf
!  formatted printing of some diagnostic output
! V1_19        2012-04-16 Harald Anlauf
!  mo_fg_cov: optimize for sxf90 rev.430+
!  geterr: hack to protect from crashes due to negative analysis errors
! V1_22        2013-02-13 Harald Anlauf
!  bugfixes for wrong book-keeping
! V1_26        2013/06/27 Harald Anlauf
!  init_spot: bugfix for interpolation weight calculation bug uncovered by ICON
!  Write diagnostics files to output directories, not the current one
! V1_27        2013-11-08 Harald Anlauf
!  geterr: limit dt(anaerr) to less than 2 weeks
!          adjust print format to catch excessive time ranges (010101)
! V1_28        2014/02/26 Robin Faulwetter
!  Introduced grtim_*_POL.
! V1_31        2014-08-21 Andreas Rhodin
!  modified rttov specific vertical interpolation (nwv_rad=3 in /observations/)
!  limit dtaerr to 72 hours, inform user about options (H.Anlauf)
!  saveguard for division/0 if background error e_h explicitly set to zero
! V1_37        2014-12-23 Andreas Rhodin
!  changes for Variational Ensemble Kalman Filter (VarEnKF)
! V1_44        2015-09-30 Harald Anlauf
!  Preparations for Jason-2
! V1_45        2015-12-15 Harald Anlauf
!  OpenMP parallelization
! V1_46        2016-02-05 Harald Anlauf
!  Fix DRIBU station id; extend SYNOP/DRIBU code for monitoring of T2m
! V1_48        2016-10-06 Andreas Rhodin
!  init_spot, get_corr: changes to handle (passive) reports
!  get_corr: correctly count number of non-zero matrix elements
! V1_50        2017-01-09 Harald Anlauf
!  get_corr, get_corr_opt: remove static variables to make these routines thread-safe
! V1_51        2017-02-24 Harald Anlauf
!  init_spot: improve vectorization
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM  2002       original explicit OI formulation
! Andreas Rhodin  DWD    2006       extend to vertical NMC formulation
! Harald Anlauf   DWD    2008,2009  optimizations for SX8,SX9
!------------------------------------------------------------------------------
!
! Diagnostics for vectorization
!
#if defined (_FTRACE) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!==============================================================================
! Enable alternative code for innermost loops of get_corr,
! which is slightly faster on NEC SX than original version.
!#define GET_CORR_NEW_INNER
!
!#if defined (__SX__)
!#define GET_CORR_NEW_INNER
!#endif
!==============================================================================
  !-------------
  ! Modules used
  !-------------
  use mo_kind,       only: wp               ! kind parameter
  use mo_exception,  only: finish           ! abort routine
  use mo_mpi_dace,   only: dace,           &! MPI group info
                           p_bcast          ! broadcast routine
  use mo_namelist,   only: position_nml,   &! routine to position nml group
                           nnml,           &! namelist fortran unit number
                           POSITIONED       ! position_nml: OK  return flag
  use mo_physics,    only: R,              &! gas constant of dry air
                           gacc,           &! gravity acceleration
                           rearth,         &! earth radius
                           e                ! 2.718...
  use mo_nr,         only: zbrent,         &! finds zeros of a function
                           bessj0, bessj1   ! Bessel functions
  use mo_time,       only: t_time,         &! time data type
                           operator(-),    &! calculate time difference
                           operator(/=),   &! compare date and time
                           imm,            &! get month
                           hours            ! convert to hours
  use mo_dwd,        only: t_ccf,          &! data type
                           nclerr,         &! climatological error coefficients
                           nfcoun,         &! climatol. 6h forecast error cf.
                           emodsm,         &! analysis error correction factor
                           ba_err,         &! background analysis error
                           lba_err,        &! background analysis error present
                           read_OI_coef,   &! initialisation routine
                         get_analysis_error ! read analysis error
  use mo_t_obs,      only: t_vic,          &! vertical interp.coef. data type
                           t_hic1,         &! horizont.interp.coef. data type
                           t_coord,        &! coordinates
                           OBS_TV, OBS_H,  &! observation type identifier
                           OBS_RH, OBS_U,  &!
                           OBS_V,  OBS_CHI,&!
                           OBS_HS, OBS_DUM,&!
                           OBS_FF, OBS_T,  &!
                           INT_H,  INT_RH, &! interpolation type identifier
                           INT_CHI,INT_TV   !
  use mo_fdbk_tables,only: OT_RAD           ! RADIANCE observation type
  use mo_atm_grid,   only: t_grid,         &! grid data type
                           to_grads,       &! write grid to grads file
                           init_ctl         ! preset .ctl file from grid
  use mo_atm_state,  only: t_atm,          &! atmospheric state data type
                           construct,      &! construct atmospheric state
                           destruct,       &! construct atmospheric state
                           to_grads         ! write atmosph. to grads file
  use mo_grads,      only: t_ctl,          &! grads  .ctl file data type
                           init_ctl,       &! preset .ctl file
                           add_var,        &! add a variable entry to ctl-file
                           write_var,      &! write data set to grads file
                           write_ctl,      &! write  .ctl file
                           destruct         ! clean up .ctl data type
  use mo_run_params, only: aux,            &! path for additional files
                           path_file        ! add path to a file name
  use mo_hum_ana,    only: ln_hum_ana_top, &! no humidity analysis above
                           e_hum_top        ! nominal humidity error above
  use mo_dec_matrix, only: mp               ! matrix coeff. kind parameter
! use coordinates,   only: cartesian        !+++ required for IBM clf V12.1
  use mo_grid_intpol,only: grid_indices     ! get indices for interpolation
  use mo_t_bg_err_op,only: covm,           &! NMC fitted covariances
                           set_vertinpc,   &! set vertical interp.coefficients
                           set_horinpc      ! set horizont.interp.coefficients
  use mo_t_obs,      only: t_obs,          &! observation data type
                           t_spot,         &!   component of t_obs
                           nwv_rad,        &! radiance vert.intp. flag
                           ITY_ICOL,       &! interpolation type: column
                           ITY_ICOLS,      &!                     columns
                           ITY_MCOLS        !               model columns
  implicit none
  !================
  ! public entities
  !================
  private
  !---------------
  ! initialization
  !---------------
  public :: init_fg_cov       ! initialisation routine for this module
  public :: read_nml_pscmodel ! read namelist /pscmodel/
  !------------------------------------------------------
  ! high level interface to calculate covariance matrices
  !------------------------------------------------------
  public :: t_rowcol    ! data type to hold precalculated information
  public :: t_dwd_spot  !   component of t_rowcol
  public :: t_dwd_col   !   component of t_rowcol
  public :: construct   ! construct type t_rowcol
  public :: zero        ! zero      type t_rowcol for reuse
  public :: destruct    ! destruct  type t_rowcol
  public :: init_spot   ! precalculate information for a spot
  public :: init_spot_obs ! set corr. data type from obs.
  public :: get_corr    ! get correlations/covariances for pair of spots
  public :: get_corr_opt! get correlations/covariances for pair of spots
  public :: IFC_CL      ! climatological error    flag (passed to init_spot)
  public :: IFC_CFC     ! climatological fc error flag (passed to init_spot)
  public :: IFC_AFC     ! actual         fc error flag (passed to init_spot)
  !----------------------------
  ! global flags and parameters
  !----------------------------
  public :: x_cut       ! cutoff length scale
  public :: L_max       ! maximum horizontal length scale (unit sphere)
  public :: Lh_max      ! maximum horizontal length scale (m)
  public :: gradsfile   ! name of plotfile (grads)
  public :: test_bl     ! flag to test blockdiagonal of HBHt for pos.def.
  public :: sec2rad     ! flag to convert secant to radiant
  public :: model       ! horizontal covariance model to use
  public :: a           ! additional parameter
  public :: L_h         ! (constant) hor. length scale from namelist
  public :: L_q         ! (constant) hor. length scale from namelist
  public :: rnu         ! e2_x / e2_v
  public :: rmu         ! constant geostrophic coupling (test)
  public :: lsvv        ! largescale contrib. to wind-wind cor.
  public :: e_h         ! (constant) height error
  public :: e_q         ! (constant) relative humidity error
  public :: e_v         ! (constant) wind error
  public :: emin_q      ! global lower bound on RH error
  !----------------------------------------------------
  ! low level routines (horizontal covariance function)
  !----------------------------------------------------
  public :: corr        ! horizontal correlation function, secant is argument
  public :: dcorr       ! derivatives of horizontal correlation function
  public :: corrh       ! horizontal correlation function, rel.humidity
  public :: set_h_corr  ! initialisation routine for horizontal correlations
  public :: L1, L1r     ! length scales
  public :: cor_h       ! as corr,  distance (m) is argument
  public :: cor_rh      ! as corrh, distance (m) is argument
  !------------------------------------------------------------------------
  ! make public for external procedures bess,d_bessj0 (workaround xlf90 bug)
  !------------------------------------------------------------------------
! public :: j1,rkn,a_dwd,rn,x_cuts,j0,na
  !----
  ! HPM
  !----
!#include "f_hpm.h"
!------------------------------------------------------------------------------
  ! Moved here from module mo_necjb_opt.f90 (J.Boerhout, NEC):
  type t_admin_gc
    integer :: i, j, iv
  end type t_admin_gc

  type(t_admin_gc), allocatable :: tv_gc(:)
  integer                       :: nv_gc

  public :: t_admin_gc, tv_gc, nv_gc
!------------------------------------------------------------------------------
  !==========================================
  ! module data type for intermediate storage
  !   of covariance matrix information
  !==========================================
  !------------------------------------
  ! hold information concerning a level
  !------------------------------------
  type t_dwd_col
    integer  :: i   ! number of consecutive levels
    real(wp) :: z   ! vertical coordinate (ln p)
    integer  :: it  ! observation type
    integer  :: itl ! observation type (.or.ed over levels)
    real(wp) :: e   ! error (std. deviation)
    real(wp) :: xv  ! transformed vertical coordinate for velocity potential
    real(wp) :: ev  ! error (std. deviation)          for wind
    real(wp) :: xq  ! transformed vertical coordinate for relative humidity
    real(wp) :: eq  ! error (std. deviation)          for relative humidity
    real(wp) :: xh  ! transformed vertical coordinate for height
    real(wp) :: eh  ! error (std. deviation)          for height
    real(wp) :: xhz ! derivative of transformed vertical coordinate for height
    real(wp) :: ehz ! derivative of error (std. deviation)          for height
    real(wp) :: et  ! error (std. deviation)          for thickness
    real(wp) :: ezn ! ehz / et
    real(wp) :: exzn! eh * xhz / et
    !-------------------------------------------------------------------------
    type(t_vic) :: vi ! vertical interpolation coefficients
  end type t_dwd_col
  !-----------------------------------
  ! hold information concerning a spot
  !-----------------------------------
  type t_dwd_spot
    !--------------------
    ! general bookkeeping
    !--------------------
    integer                  :: ih        ! spot index
    integer                  :: i         ! column start index
!   integer                  :: ic        ! matrix start index
    integer                  :: n         ! number of degrees of freedom
    !------------------------------------
    ! required for horizontal covariances
    !------------------------------------
    real(wp)                 :: x (3)     ! vector on the unit sphere
    real(wp)                 :: du(3)     ! direction vector for u
    real(wp)                 :: dv(3)     ! direction vector for v
    real(wp)                 :: xt(3)     ! vector on the transformed sphere
    real(wp)                 :: L_h       ! horizontal length scale     (T, h)
    real(wp)                 :: L_q       ! length scale for moisture
    real(wp)                 :: mu        ! geostrophic balance coefficient
    real(wp)                 :: emod      ! bg error modification factor
    !-------------------------------------------
    ! required for interpolation of coefficients
    !-------------------------------------------
    integer                  :: jlat      ! latitude index (correlation files)
    integer                  :: jlat1     ! jlat + 1
    integer                  :: ilon      ! longitude index (correlation files)
    integer                  :: ilon1     ! ilon + 1
    real(wp)                 :: wlat      ! latitudinal  weight (corr.   files)
    real(wp)                 :: wlon      ! longitudinal weight (corr.   files)
    !----------------------------------
    ! required for vertical covariances
    !----------------------------------
    real(wp)                 :: czh (7)   ! height error coefficints    (T, h)
    real(wp)                 :: cxh (7)   ! height transformation coefs (T, h)
    real(wp)                 :: cch       ! height vcorr. coefficient   (T, h)
    real(wp)                 :: czq (7)   ! humid. error coefs          (q)
    real(wp)                 :: cxq (7)   ! humid. transformation coefs (q)
    real(wp)                 :: ccq       ! humid. vcorr. coefficient   (q)
    real(wp)                 :: czv (7)   ! wind   error coefficients   (u,v)
    real(wp)                 :: cxv (7)   ! wind   transformation coefs (u,v)
    real(wp)                 :: ccv       ! wind   vcorr. coefficient   (u,v)
    type(t_dwd_col) ,pointer :: col(:)    ! column information
    !------------------------------------------------------------------------
    type(t_hic1)             :: h         ! horizontal interpolation coeff.
    real(wp)    ,allocatable :: vm(:,:,:) ! vertical modes
  end type t_dwd_spot
  !--------------------------------------------
  ! hold information concerning a row or column
  !--------------------------------------------
  type t_rowcol
    integer                   :: n                 ! no of observ. in row/colum
    integer                   :: nh                ! no of spots in row/column
    integer                   :: in                ! n set so far
    integer                   :: inh               ! largest nh set so far
    type (t_dwd_spot),pointer :: spots(:) =>NULL() ! spots
    type (t_dwd_col) ,pointer :: cols (:) =>NULL() ! column information
  end type t_rowcol

!==============================================================================
  !=================
  ! Module variables
  !=================

  !-----------------------------------
  ! data for Bessel function expansion
  !-----------------------------------
  integer  ,parameter :: mkn = 20        ! max# of coefficients, must be even
  integer  ,parameter :: na  =  8        ! number of coefficients in DWD OI
  real(wp) ,save      :: rkn   (0:mkn)   ! zeros of Bessel func. derivatives
  !--------------------------------------------------------------
  ! Bessel function expansion coefficients A_0 .. A_8
  ! cf. RM1, page 2.13
  ! note: A_0 will be modified so that corr(1)==0
  ! as specified in OI documentation and OI blockdata statement :
  !--------------------------------------------------------------
! real(wp) :: a_dwd (0:na) = (/0.15_wp,0.30_wp,0.28_wp,0.14_wp,&
!                              0.08_wp,0.05_wp,0.03_wp,0.02_wp,0.01_wp/)
  !-------------------------------
  ! as specified in file f_combes:
  !-------------------------------
  real(wp) :: a_dwd (0:na) = (/0.141326775770938E+00_wp,&! 0
                               0.282653551541877E+00_wp,&! 1
                               0.263809981439085E+00_wp,&! 2
                               0.131904990719542E+00_wp,&! 3
                               0.753742804111672E-01_wp,&! 4
                               0.471089252569796E-01_wp,&! 5
                               0.278790619670805E-01_wp,&! 6
                               0.178448608873438E-01_wp,&! 7
                               0.120975720059924E-01_wp/)! 8
  real(wp) ,save :: x_zero           ! x for corr(x)==0 for the first time
  real(wp) ,save :: x_cut            ! corr not calculated beyond, set to 0
  real(wp) ,save :: x_cutdr          ! x_cut / earth radius
  real(wp) ,save :: x_cuts           ! x_cut  scaled by L1
  real(wp) ,save :: L_max  = 4.0_wp  ! max. hor. length scale on unit sphere
  real(wp) ,save :: Lh_max = 4.e7_wp ! max. hor. length scale (m)

  !--------------------
  ! namelist /PSCMODEL/
  !--------------------
  integer,parameter:: IFC_CL = 1        !   climatological error
  integer,parameter:: IFC_CFC= 2        !   climatological fc error
  integer,parameter:: IFC_AFC= 3        !   actual         fc error
  integer          :: ifc    = IFC_AFC  ! flag for fg-error type
  character(len=8) :: model  = 'dwd-oi' ! horizontal correlation model
  logical          :: rhcexp = .true.   ! RH horiz. corr. = exp(-.5(L/L_q)**2)
  integer          :: irad   = -1       ! radian or secant is argument to c.f.
  real(wp)         :: L_h    =  0._wp   ! horizontal length scale for height
  real(wp)         :: L_q    =  0._wp   ! horizontal length scale humidity
  real(wp)         :: a_h    =  0._wp   ! vertical tuning parameter height
  real(wp)         :: a_q    =  0._wp   ! vertical tuning parameter humidity
  real(wp)         :: a_x    =  0._wp   ! vertical tuning parameter vel.pot.
  real(wp)         :: t_h    =  0._wp   ! linear transformation coef. height
  real(wp)         :: t_q    =  0._wp   ! linear transformation coef. humidity
  real(wp)         :: t_x    =  0._wp   ! linear transformation coef. vel.pot.
  real(wp)         :: e_h    = -1._wp   ! height error
  real(wp)         :: e_q    = -1._wp   ! relative humidity error
  real(wp)         :: e_v    = -1._wp   ! wind error
  real(wp)         :: rnu    = 0.1_wp   ! e2_x / e2_v
  real(wp)         :: rmu    = -9._wp   ! constant geostrophic coupling (test)
  real(wp)         :: cut    =  0._wp   ! cutoff radius (multiple of L)
  logical          :: test_bl= .false.  ! flag to test blockdiagonal of HBHt
  integer          :: fixbc  = 1        ! fix bessel coefs.: d2c/dr2|r=1 = 0
  real(wp)         :: modcvx = 0.08_wp  ! modified vertical corr. for div.wind
  real(wp)         :: a      = 0._wp    ! scale parameters for cov. models
  real(wp)         :: nw     = na       ! number of coefficients for wind corr.
  logical          :: thc    = .false.  ! transform horizontal coordinates
  real(wp)         :: lsvv   = 0._wp    ! largescale contrib. to wind-wind cor.
  real(wp)         :: scale  = 1._wp    ! scaling factor for background error
  real(wp)         :: emin_q = 0.00_wp  ! global lower bound on RH error
  real(wp)         :: grtim_NH_POL = 144._wp! time for errors to grow to
  real(wp)         :: grtim_NH     = 144._wp!   random state (hours)
  real(wp)         :: grtim_TR     =  48._wp!   n. pole, n. hemisphere, tropics
  real(wp)         :: grtim_SH     = 144._wp!   s. hemisphere
  real(wp)         :: grtim_SH_POL = 144._wp!   s. pole
  real(wp)         :: emodmin  =   0._wp! minimum value for modification factor
  real(wp)         :: ptop_v = 10._wp   ! top    pressure for polynom.exp.(hPa)
  real(wp)         :: pbot_v = 1100._wp ! bottom pressure for polynom.exp.(hPa)
  logical          :: zflip  = .false.  ! flip variances in the vertical
  logical          :: shortc = .false.  ! get_corr: shortcut single level data

  !--------------------------------------------
  ! valid correlation model parameters (model):
  !--------------------------------------------
  ! 'dwd-oi'  ! Bessel function expansion as used in DWD-OI
  ! 'c5pr_12' ! compactly supported 5. order piecewise rational functions a=1/2
  ! 'c5pr_inf'! compactly supported 5. order piecewise rational functions a=inf
  ! 'csxoar1' ! compactly supporting functions approximating SOAR (Gneiting)
  ! 'csxoar2' ! compactly supporting functions approximating TOAR (Gneiting)
  !
  ! for 'c5pr_..' cf. Gaspari+Cohn QJRMS (1999),125,723-757
  !
  ! for 'csxoar.' cf. Gneiting QJRMS (1999),125,2449-2464
  !
  namelist /PSCMODEL/ & ! namelist /PSCMODEL/ physical space covariance model
    model, rhcexp, L_h, L_q, a_h, a_q, a_x, t_h, t_q, t_x, e_h, e_q, e_v,    &
    rnu, rmu, cut, ifc, test_bl, irad, fixbc, modcvx, a, nw, thc, lsvv,      &
    scale, emin_q, grtim_NH_POL, grtim_NH, grtim_TR, grtim_SH, grtim_SH_POL, &
    emodmin, ptop_v, pbot_v, zflip, shortc

  !------------------------------------------------------------
  ! temporary variables for horizontal correlation calculations
  ! calculated for a given normalized radius
  !------------------------------------------------------------
  real(wp) :: rn = -1._wp! radius
  integer  :: ir         ! int (rn)
  real(wp) :: L1         ! second derivative d2 corr / d r2
  real(wp) :: L1r        ! L1 times earth radius
  real(wp) :: L2         ! L1^2
  real(wp) :: edL1       ! 1./L1
! real(wp) :: edr        ! 1. / earth radius
  logical  :: sec2rad    ! flag to convert secant to radiant
  real(wp) :: J0 (1:mkn) ! bessel function coefficients
  real(wp) :: J1 (1:mkn) ! bessel function coefficients
  save     :: rn, ir, L1, L2, edL1, L1r, sec2rad, J0, J1
  !-------------------------------------------------------------------------
  ! coefficients for
  ! compactly supported 5th order piecewise rational functions C0(r,a=1/2,c)
  ! cf. Gaspari+Cohn QJRMS (1999),125,723-757
  !-------------------------------------------------------------------------
  real(wp) ,parameter :: c05 = - 1._wp / 4._wp
  real(wp) ,parameter :: c04 =   1._wp / 2._wp
  real(wp) ,parameter :: c03 =   5._wp / 8._wp
  real(wp) ,parameter :: c02 = - 5._wp / 3._wp
  real(wp) ,parameter :: c00 =   1._wp
  real(wp) ,parameter :: c15 =   1._wp /12._wp
  real(wp) ,parameter :: c14 = - 1._wp / 2._wp
  real(wp) ,parameter :: c13 =   5._wp / 8._wp
  real(wp) ,parameter :: c12 =   5._wp / 3._wp
  real(wp) ,parameter :: c11 = - 5._wp
  real(wp) ,parameter :: c10 =   4._wp
  real(wp) ,parameter :: c1x = - 2._wp / 3._wp
  !-------------------------------------------------------------------------
  ! coefficients for
  ! compactly supported 5th order piecewise rational functions C0(r,a=inf,c)
  ! cf. Gaspari+Cohn QJRMS (1999),125,723-757
  !-------------------------------------------------------------------------
  real(wp) ,parameter :: f05 = -  28._wp /  33._wp
  real(wp) ,parameter :: f04 =     8._wp /  11._wp
  real(wp) ,parameter :: f03 =    20._wp /  11._wp
  real(wp) ,parameter :: f02 = -  80._wp /  33._wp
  real(wp) ,parameter :: f00 =     1._wp

  real(wp) ,parameter :: f15 =    20._wp /  33._wp
  real(wp) ,parameter :: f14 = -  16._wp /  11._wp
  real(wp) ,parameter :: f12 =   100._wp /  33._wp
  real(wp) ,parameter :: f11 = -  45._wp /  11._wp
  real(wp) ,parameter :: f10 =    51._wp /  22._wp
  real(wp) ,parameter :: f1x = -   7._wp /  44._wp

  real(wp) ,parameter :: f25 = -   4._wp /  11._wp
  real(wp) ,parameter :: f24 =    16._wp /  11._wp
  real(wp) ,parameter :: f23 = -  10._wp /  11._wp
  real(wp) ,parameter :: f22 = - 100._wp /  33._wp
  real(wp) ,parameter :: f21 =     5._wp
  real(wp) ,parameter :: f20 = -  61._wp /  22._wp
  real(wp) ,parameter :: f2x =   115._wp / 132._wp

  real(wp) ,parameter :: f35 =     4._wp /  33._wp
  real(wp) ,parameter :: f34 = -   8._wp /  11._wp
  real(wp) ,parameter :: f33 =    10._wp /  11._wp
  real(wp) ,parameter :: f32 =    80._wp /  33._wp
  real(wp) ,parameter :: f31 = -  80._wp /  11._wp
  real(wp) ,parameter :: f30 =    64._wp /  11._wp
  real(wp) ,parameter :: f3x = -  32._wp /  33._wp
  !----------------------------------------------------------
  ! coefficients for
  ! compactly supporting functions approximating SOAR or TOAR
  ! cf. Gneiting QJRMS (1999),125,2449-2464
  !----------------------------------------------------------
  real(wp), save :: gna, gnb, gnc, gne, gnf ,gng, gni, gnj, gnk
  !----------------
  ! other variables
  !----------------
  character(len=128) :: gradsfile= ''      ! name of plotfile (grads)
  logical            :: samegrid = .false. ! analysis error on 6deg global grid
  !---------------------------------
  ! checks for module initialization
  !---------------------------------
  logical :: linit_fg_cov = .false.  ! T if this module has been initialized
  logical :: lset_h_corr  = .false.  ! correlation model set? (bess.fct.coefs)
  type(t_time) ,save :: fg_time      ! date+time  first guess error is valid
  integer      ,save :: fg_month     ! month      first guess error is valid
  !-----------------------------------------------------------
  ! routines (passed as actual argument) declared as externals
  ! (required due to IBM xlf90 bug ?)
  !-----------------------------------------------------------
! real(wp) ,external :: bess, d_bessj0
  !-----------
  ! interfaces
  !-----------
  interface construct
    module procedure construct_rowcol
  end interface construct

  interface zero
    module procedure zero_rowcol
  end interface zero

  interface destruct
    module procedure destruct_rowcol
    module procedure destruct_rowcols
  end interface destruct
!==============================================================================
contains
!==============================================================================
  subroutine init_fg_cov (time)
  !----------------------------
  ! Initialize Covariance model
  !----------------------------
  type (t_time) ,intent(in) :: time ! analysis time
    if (.not. linit_fg_cov) then
      linit_fg_cov = .true.
      fg_time      = time
      fg_month     = imm (time)
      !-----------
      ! read files
      !-----------
      call read_nml_pscmodel
      call set_h_corr
      call read_OI_coef
      call get_analysis_error
      !-----------------------------------------
      ! calculate maximum value for length scale
      !-----------------------------------------
      select case (ifc)
      case (IFC_CL)
        Lh_max = maxval(nclerr% b)
      case (IFC_CFC, IFC_AFC)
        Lh_max = maxval(nfcoun% b)
      case default
        call finish('init_fg_cov','invalid flag ifc')
      end select
      L_max = Lh_max / rearth
!     edr   = 1._wp  / rearth
      !-----------------------------------------------
      ! calculate forecast error scaling factor emodsm
      !-----------------------------------------------
      if (.not. lba_err .and. ifc==IFC_AFC) call finish ('init_fg_cov',&
        'ifc==IFC_AFC but no analysis error file read')
      if (allocated(emodsm)) emodsm(:,:) = 1._wp
      if (ifc == IFC_AFC) call geterr (time)

      if (dace% lpio) print '(A,F15.9)','init_fg_cov: L_max =',L_max
    endif

    if (time /= fg_time) call finish ("init_fg_cov","time /= fg_time")

  end subroutine init_fg_cov
!------------------------------------------------------------------------------
  subroutine geterr (time)
  type (t_time) ,intent(in) :: time ! analysis time
!
!     method
!     ------
!
!        a row of first guess errors of u, v and height are calculated,
!     using:
!
!     ferr(jlon,jlev,jvar,jrow)=emodsm(jlon,jrow)*clfgrw(jlev,jvar,jrow)
!
!     where ferr   = output first guess error values
!           emodsm = smoothed error modification factors
!           clfgrw = climatological 6 hour forecast errors
!           jvar   = variable number (1=u, 2=v, 3=z)
!
!        emodsm is derived from emod by applying a 4 point smoothing
!     function at interior rows, or a 3 point smoother at polar rows.
!     the error modification factor is defined as:
!
!     emod(jlon,jrow)=anorm(jlon,jrow)+dtaerr/zgrtim*
!                     (sqrt(2.)*clnorm(jrow)-anorm(jlon,jrow))
!
!     where anorm  = normalised analysis error factor field
!           clnorm = normalised climatological error field
!           dtaerr = time difference in seconds between the analysis
!                    error file time and the current analysis
!           zgrtim = time for errors to grow to random state
!
!     by default, zgrtim = 144 hours in the extra-tropics
!                        =  48 hours in the tropics
!     with the tropics defined as the region (20,-20), and the extra-
!     tropics defined as (90,30) and (-30,-90). a smoothly mixed value
!     of zgrtim is used between the 2 regions.
!
!     the normalised analysis error factor field is calculated as:
!
!                              nevar  nelev
!                               ---    ---
!                               \      \    aerr(jlon,jlev,jvar,jrow)
!     anorm(jlon,jrow)= zrlv  *  .   *  .   -------------------------
!                               /      /     clfgrw(jlev,jvar,jrow)
!                               ---    ---
!                             jvar=1  jlev=1
!
!     where zrlv = 1./(nelev*nevar)
!           aerr = input analysis errors
!
!     the normalised climatological error factor field is defined as:
!
!                            nevar  nelev
!                             ---    ---
!                             \      \   clerrw(jlev,jvar,jrow)
!     clnorm(jrow) =  zrlv  *  .   *  .  ----------------------
!                             /      /   clfgrw(jlev,jvar,jrow)
!                             ---    ---
!                           jvar=1 jlev=1
!
!     where clerrw = climatological errors.
!
!        ferr is calculated for row merow. because the smoothing
!     function requires 3 rows of data, the analysis errors being
!     processed are for row merow+1.
!
    !----------------
    ! local variables
    !----------------
    integer ,parameter     :: nvar = 3
    integer                :: nx, ny, nz
    integer                :: nobs
    type(t_rowcol)         :: rc
    real(wp)               :: sq2
    real(wp)               :: an, cn, tn
    real(wp)               :: dtaerr, zgrtim, adlat
    real(wp)               :: smcfer
    real(wp)               :: x   (nvar)
    real(wp)               :: z   ((ba_err% grid% nz) * nvar)
    integer                :: it  ((ba_err% grid% nz) * nvar)
    real(wp)               :: ea  ((ba_err% grid% nz) * nvar)
    real(wp)               :: ef  ((ba_err% grid% nz) * nvar)
    real(wp)               :: ec  ((ba_err% grid% nz) * nvar)
    real(wp)               :: emod (ba_err% grid% nx, ba_err% grid% ny)
    integer                :: i,j,k,l
    type (t_grid) ,pointer :: g
    type (t_atm)           :: fc_err ! forecast error
    type (t_atm)           :: cl_err ! climatological error
    type (t_ctl)           :: ctl
    type (t_coord)         :: c      ! coordinates
    real(wp),    parameter :: dtmax = 3*24    ! (hours) Limit dt to 3 days
    !---------------------------------
    ! set variables to constant values
    !---------------------------------
    sq2  =  sqrt(2._wp)  ! sqrt (2)
    g    => ba_err% grid ! pointer to analysis grid
    nx   =  g% nx        ! number of longitudes
    ny   =  g% ny        ! number of latitudes
    nz   =  g% nz        ! number of levels
    nobs =  nvar * nz    ! total number of observation points
    x(:) = 0._wp         ! dummy variable
    do k=1,nz
      l = (k - 1) * nvar
!NEC$ ivdep
      z (l+1:l+nvar) = log (g% akf (k))         ! vertical coordinate
      it(l+1:l+nvar) = (/OBS_U, OBS_V, OBS_H/)  ! observation type
    end do
    dtaerr = hours (time - ba_err% time)        ! forecast time interval

    if (dace% lpio) &
       print '(A,F12.3)','geterr: forecast time interval [hours] =',dtaerr

    if(dtaerr < 0._wp) dtaerr = 3._wp
    if(dtaerr > dtmax) then
       if (dace% lpio) then
          write(*,'(A,F12.3)') 'geterr: *** dtaerr [hours] =',dtaerr
          write(*,'(A,F12.3)') '            maximum allowed:',dtmax
          write(*,'()')
          write(*,'(A)') 'If the intention is to run an analysis on a historical date,'
          write(*,'(A)') 'there exist the following options'
          write(*,'(A)') '1) Use the analysis error from a future analysis as a proxy'
          write(*,'(A)') '2) To use the climatological fc error, set:'
          write(*,'(A)') "   namelist /PSCMODEL/: ifc = 2"
          write(*,'(A)') "   and namelist /RUN/ : oldanerr_file = '/none/'"
       end if
       call finish ('geterr','bad date of previous analysis error')
    end if
    !----------------------------------
    ! allocate forecast error data type
    !----------------------------------
    call construct (rc, nobs, 1)
    !---------------------------------------------------
    ! construct climatological error data (for plotting)
    !---------------------------------------------------
    if (dace% lpio) then
      call construct (fc_err, template = ba_err)
      call construct (cl_err, template = ba_err)
    endif
    !-------------------------------------------
    ! loop over latitudes of analysis error grid
    !-------------------------------------------
    do j=1,ny
      !--------------------------------------------------------------------
      ! calculate time factor tn = dtaerr/zgrtim
      !
      !   where dtaerr = time difference between the analysis
      !                  error file time and the current analysis
      !         zgrtim = time for errors to grow to random state
      !
      !   by default, zgrtim = 144 hours in the extra-tropics
      !                      =  48 hours in the tropics
      !   with the tropics defined as the region (20,-20), and the extra-
      !   tropics defined as (90,30) and (-30,-90). a smoothly mixed value
      !   of zgrtim is used between the 2 regions.
      !-------------------------------------------------------------------
      adlat = g% dlat(j)
      if (adlat >= 60) then
        zgrtim = grtim_NH_POL
      else if ((adlat <= 60).and.(adlat >= 50))  then
        zgrtim = ((grtim_NH_POL - grtim_NH )                  &
          /(60._wp - 50._wp)) * (adlat - 50._wp) + grtim_NH
      else if ((adlat <= 50).and.(adlat >= 30)) then
        zgrtim = grtim_NH
      else if ((adlat <= 30).and.(adlat >= 20))  then
        zgrtim = ((grtim_NH - grtim_TR )                  &
          /(30._wp - 20._wp)) * (adlat - 20._wp) + grtim_TR
      else if ((adlat < 20).and.(adlat > -20))  then
        zgrtim = grtim_TR
      else if ((adlat >= -30).and.(adlat <= -20))  then
        zgrtim = ((grtim_SH - grtim_TR )                  &
          /(20._wp - 30._wp)) * (adlat + 20._wp) + grtim_TR
      else if ((adlat >= -50).and.(adlat <= -30))  then
        zgrtim = grtim_SH
      else if ((adlat >= -60).and.(adlat <= -50))  then
        zgrtim = ((grtim_SH_POL - grtim_SH )                  &
          /(50._wp - 60._wp)) * (adlat + 50._wp) + grtim_SH
      else if (adlat <= -60)  then
        zgrtim = grtim_SH_POL
      endif
      tn = min (1._wp, dtaerr / zgrtim)

if (dace% lpio) write(6,'(a,i4,4f12.6)') &
'geterr: j,adlat,zgrtim,dtaerr,tn  =',     &
         j,adlat,zgrtim,dtaerr,tn

      !--------------------------------------------
      ! loop over longitudes of analysis error grid
      !--------------------------------------------
      do i=1,nx
        c = t_coord (g% dlon(i) ,&! longitude (degree)
                     g% dlat(j) ,&! latitude  (degree)
                     x (:)      ,&! unit vector on the sphere (location)
                     x (:)      ,&! unit vector on the sphere (u-direction)
                     x (:)       )! unit vector on the sphere (v-direction)
        !----------------------------------
        ! analysis error without correction
        !----------------------------------
        rc% in        = 0
        rc% spots% ih = 0
        call init_spot (rc,         &
                        1          ,&! spot index within the box
                        c          ,&! coordinates
                        it(:)      ,&! observation type identifiers
                        z (:)      ,&! vertical coordinates (ln p)
                    lcl=IFC_CL     ,&! flag for climatological error
                    err=ec(:))       ! error std deviation
        !--------------------------------------
        ! get forecast error without correction
        !--------------------------------------
        rc% in        = 0
        rc% spots% ih = 0
        call init_spot (rc,         &
                        1          ,&! spot index within the box
                        c          ,&! coordinates
                        it(:)      ,&! observation type identifiers
                        z (:)      ,&! vertical coordinates (ln p)
                    lcl=IFC_CFC    ,&! flag for fc error without correction
                    err=ef(:))       ! error std deviation
        !--------------------------------------
        ! store statistical errors for plotting
        !--------------------------------------
        if (dace% lpio) then
          do k=1,nz
            l = (k - 1) * nvar
            fc_err% u    (i,j,k,1) = ef (l+1)
            fc_err% v    (i,j,k,1) = ef (l+2)
            fc_err% geof (i,j,k,1) = ef (l+3)
          end do
          do k=1,nz
            l = (k - 1) * nvar
            cl_err% u    (i,j,k,1) = ec (l+1)
            cl_err% v    (i,j,k,1) = ec (l+2)
            cl_err% geof (i,j,k,1) = ec (l+3)
          end do
        endif
        !-------------------
        ! get analysis error
        !-------------------
        do k=1,nz
          l = (k - 1) * nvar
          ea (l+1) = ba_err% u    (i,j,k,1)
          ea (l+2) = ba_err% v    (i,j,k,1)
          ea (l+3) = ba_err% geof (i,j,k,1)
        end do
        !--------------------------------------------------------
        ! Safeguard from bad analysis errors (hack)
        ! (not if background error e_h is set to zero explicitly)
        !--------------------------------------------------------
        if (e_h /= 0._wp) then
          if (minval (ef) <= 0) then
            write(0,*) "geterr: bad forecast error (ef):", minval (ef)
            write(*,*) "geterr: bad forecast error (ef):", minval (ef)
            ef = max (ef, 1.e-3_wp)
          end if
          if (minval (ea) <= 0) then
            write(0,*) "geterr: bad analysis error (ea):", minval (ea)
            write(*,*) "geterr: bad analysis error (ea):", minval (ea)
            ea = max (ea, 1.e-3_wp)
          end if
        end if
        !--------------------------------------------------------
        ! calculate correction factor
        ! (not if background error e_h is set to zero explicitly)
        !--------------------------------------------------------
        if (e_h /= 0._wp) then
          an = sum (ea / ef) / nobs
          cn = sum (ec / ef) / nobs
          emod (i,j) = max (emodmin, an + max (0._wp, tn * ( sq2 * cn - an)))
        else
          emod (i,j) = 1._wp
        endif
      end do
    end do
    !-------------------------
    ! smooth correction factor
    !-------------------------
    smcfer = 0.1_wp
    emodsm (:,1 ) = (1 - 3 * smcfer) * emod (:,1 ) + smcfer * emod (:,2)
    emodsm (:,ny) = (1 - 3 * smcfer) * emod (:,ny) + smcfer * emod (:,ny-1)
    emodsm (:,2:ny-1) = (1 - 4 * smcfer) * emod (:,2:ny-1)       &
                      + smcfer * (emod (:,1:ny-2) + emod (:,3:ny))
    emodsm (2:nx-1,:) = emodsm (2:nx-1,:) &
                      + smcfer * (emod (1:nx-2,:) + emod (3:nx,:))
    emodsm (1 ,:)     = emodsm (1 ,:) +  smcfer * emod (2   ,:) &
                                      +  smcfer * emod (nx  ,:)
    emodsm (nx,:)     = emodsm (nx,:) +  smcfer * emod (nx-1,:) &
                                      +  smcfer * emod (1   ,:)
    emodsm            = emodsm * scale

    if (dace% lpio) then

      write(6,'()')
      write(6,'(a,2f6.3)') 'geterr: emod   =',minval(emod),  maxval(emod)
      write(6,'(a,2f6.3)') '        emodsm =',minval(emodsm),maxval(emodsm)
      write(6,'()')
      !---------------------
      ! write GRADS plotfile
      !---------------------
      gradsfile = '^fg_err.grads' ;gradsfile = path_file (aux, gradsfile)
      call init_ctl  (ctl, gradsfile, ba_err% grid, comment =        &
                        (/'first guess and analysis error         ', &
                          'time slices: 1: climatological error   ', &
                          '             2: climat. forecast error ', &
                          '             3: previous analysis error', &
                          '             4: actual forecast error  ', &
                          '             5: actual analysis error  ', &
                          '             6: reduction: fc-ana error'/))
      call write_var (ctl, emod,'emod',comment='error modification factor')
      call write_var (ctl, emodsm, 'emodsm',&
                           comment='smoothed modification factor')
!     call to_grads  (ctl, ba_err% grid)
      call to_grads  (ctl, cl_err, t=1, comment='t: cl fc an fg an error')
      call add_var   (ctl, 't')
      call add_var   (ctl, 'rh')
      call to_grads  (ctl, fc_err, t=2)
      call to_grads  (ctl, ba_err, t=3)
      call write_ctl (ctl)
      call destruct  (ctl)
    end if
    !--------
    ! cleanup
    !--------
    call destruct (rc)
    if (dace% lpio) then
      call destruct (fc_err)
      call destruct (cl_err)
    endif
    if(nx==60 .and. ny==30 .and. g%di==6._wp .and. g%dj==6._wp) samegrid=.true.
  end subroutine geterr
!------------------------------------------------------------------------------
  subroutine read_nml_pscmodel
    integer :: ierr
    !-------------
    ! set defaults
    !-------------
    ifc          = IFC_AFC
    model        = 'dwd-oi'
    rhcexp       = .true.   ! RH horizontal corr.model is: exp(-.5(L/L_q)**2)
    irad         = -1
    L_h          =  0._wp
    L_q          =  0._wp
    a_h          =  0._wp
    a_q          =  0._wp
    a_x          =  0._wp
    e_h          = -1._wp
    e_q          = -1._wp
    e_v          = -1._wp
    t_h          =  0._wp
    t_q          =  0._wp
    t_x          =  0._wp
    rnu          = 0.1_wp
    rmu          = -9._wp   ! constant geostrophic coupling (test)
    cut          =  0._wp
    a            =  0._wp
    test_bl      = .false.  ! flag to test blockdiagonal of HBHt for pos.def.
    fixbc        = 1        ! fix bessel coefficients: 0=no fix; 1=shift; 2=scale
    nw           = na       ! number of coefficients for wind corr.
    modcvx       = 0.08_wp  ! modified vertical correlations for divergent wind
    thc          = .false.
    lsvv         = 0._wp    ! large scale contribution to wind-wind correlations
    scale        = 1._wp    ! scaling factor for the background error
    emin_q       = 0.00_wp  ! global lower bound on RH error
    grtim_NH     = 144._wp  ! time for errors to grow to
    grtim_TR     =  48._wp  !   random state (hours)
    grtim_SH     = 144._wp  !   n.,s.hemisphere, tropics
    grtim_NH_POL = -99._wp
    grtim_SH_POL = -99._wp
    emodmin      =   0._wp  ! minimum value for modification factor
    ptop_v       = 10._wp   ! top    pressure for polynom. expansion (hPa)
    pbot_v       = 1100._wp ! bottom pressure for polynom. expansion (hPa)
    zflip        = .false.  ! flip variances in the vertical
    shortc       = .false.  ! get_corr: shortcut for single level data
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('PSCMODEL', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=PSCMODEL, iostat=ierr)
        if (ierr/=0) call finish ('read_nml_pscmodel',          &
                                  'ERROR in namelist /PSCMODEL/')
#else
        read (nnml ,nml=PSCMODEL)
#endif
        if (grtim_NH_POL <= 0._wp) grtim_NH_POL = grtim_NH
        if (grtim_SH_POL <= 0._wp) grtim_SH_POL = grtim_SH
      end select
    endif
    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    call p_bcast (ifc          ,dace% pio)
    call p_bcast (model        ,dace% pio)
    call p_bcast (rhcexp       ,dace% pio)
    call p_bcast (irad         ,dace% pio)
    call p_bcast (L_h          ,dace% pio)
    call p_bcast (L_q          ,dace% pio)
    call p_bcast (a_h          ,dace% pio)
    call p_bcast (a_q          ,dace% pio)
    call p_bcast (a_x          ,dace% pio)
    call p_bcast (e_h          ,dace% pio)
    call p_bcast (e_q          ,dace% pio)
    call p_bcast (e_v          ,dace% pio)
    call p_bcast (t_h          ,dace% pio)
    call p_bcast (t_q          ,dace% pio)
    call p_bcast (t_x          ,dace% pio)
    call p_bcast (rnu          ,dace% pio)
    call p_bcast (rmu          ,dace% pio)
    call p_bcast (cut          ,dace% pio)
    call p_bcast (test_bl      ,dace% pio)
    call p_bcast (fixbc        ,dace% pio)
    call p_bcast (modcvx       ,dace% pio)
    call p_bcast (a            ,dace% pio)
    call p_bcast (nw           ,dace% pio)
    call p_bcast (thc          ,dace% pio)
    call p_bcast (lsvv         ,dace% pio)
    call p_bcast (ptop_v       ,dace% pio)
    call p_bcast (pbot_v       ,dace% pio)
    call p_bcast (scale        ,dace% pio)
    call p_bcast (emin_q       ,dace% pio)
    call p_bcast (grtim_NH_POL ,dace% pio)
    call p_bcast (grtim_NH     ,dace% pio)
    call p_bcast (grtim_TR     ,dace% pio)
    call p_bcast (grtim_SH     ,dace% pio)
    call p_bcast (grtim_SH_POL ,dace% pio)
    call p_bcast (emodmin      ,dace% pio)
    call p_bcast (zflip        ,dace% pio)
    call p_bcast (shortc       ,dace% pio)

    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') repeat ('-',79)
      write(6,'()')
      write(6,'(a)') '  namelist  /PSCMODEL/ :'
      write(6,'()')
      print '(a,i8)'   ,'    ifc          = ' ,ifc
      print '(a,a)'    ,'    model        = ' ,model
      print '(a,l8)'   ,'    rhcexp       = ' ,rhcexp
      print '(a,i8)'   ,'    irad         = ' ,irad
      print '(a,f15.6)','    L_h          = ' ,L_h
      print '(a,f15.6)','    L_q          = ' ,L_q
      print '(a,f15.6)','    a_h          = ' ,a_h
      print '(a,f15.6)','    a_q          = ' ,a_q
      print '(a,f15.6)','    a_x          = ' ,a_x
      print '(a,f15.6)','    e_h          = ' ,e_h
      print '(a,f15.6)','    e_q          = ' ,e_q
      print '(a,f15.6)','    e_v          = ' ,e_v
      print '(a,f15.6)','    t_h          = ' ,t_h
      print '(a,f15.6)','    t_q          = ' ,t_q
      print '(a,f15.6)','    t_x          = ' ,t_x
      print '(a,f15.6)','    rnu          = ' ,rnu
      print '(a,f15.6)','    rmu          = ' ,rmu
      print '(a,f15.6)','    cut          = ' ,cut
      print '(a,l8)'   ,'    test_bl      = ' ,test_bl
      print '(a,i8)'   ,'    fixbc        = ' ,fixbc
      print '(a,f15.6)','    modcvx       = ' ,modcvx
      print '(a,f15.6)','    a            = ' ,a
      print '(a,f15.6)','    nw           = ' ,nw
      print '(a,l8)'   ,'    thc          = ' ,thc
      print '(a,f15.6)','    lsvv         = ' ,lsvv
      print '(a,f15.6)','    ptop_v       = ' ,ptop_v
      print '(a,f15.6)','    pbot_v       = ' ,pbot_v
      print '(a,f15.6)','    scale        = ' ,scale
      print '(a,f15.6)','    emin_q       = ', emin_q
      print '(a,f15.6)','    grtim_NH_POL = ' ,grtim_NH_POL
      print '(a,f15.6)','    grtim_NH     = ' ,grtim_NH
      print '(a,f15.6)','    grtim_TR     = ' ,grtim_TR
      print '(a,f15.6)','    grtim_SH     = ' ,grtim_SH
      print '(a,f15.6)','    grtim_SH_POL = ' ,grtim_SH_POL
      print '(a,f15.6)','    emodmin      = ' ,emodmin
      print '(a,l8)'   ,'    zflip        = ' ,zflip
      print '(a,l8)'   ,'    shortc       = ' ,shortc
    endif

  end subroutine read_nml_pscmodel
!------------------------------------------------------------------------------
  function cor_h  (x, L)
  !---------------------------------------
  ! return horizontal correlation function
  !---------------------------------------
  real(wp) ,intent(in) :: x      ! horizontal distance (m)
  real(wp) ,intent(in) :: L      ! length scale
  real(wp)             :: cor_h  ! horizontal correlation
    real(wp) :: s
    s = 2._wp * sin (0.5_wp * x / rearth)
    cor_h = corr (s, L)
  end function cor_h
!------------------------------------------------------------------------------
  function cor_rh (x, L)
  !---------------------------------------------------------------
  ! return horizontal correlation for RH as used in DWD's OI
  !---------------------------------------------------------------
  real(wp) ,intent(in) :: x      ! horizontal distance (m)
  real(wp) ,intent(in) :: L      ! length scale
  real(wp)             :: cor_rh ! horizontal correlation
    real(wp) :: s
    s = 2._wp * sin (0.5_wp * x / rearth)
    cor_rh = corrh (s, L)
  end function cor_rh
!------------------------------------------------------------------------------
  function corrh (x,L)
  !---------------------------------------------------------------
  ! return horizontal correlation for RH as used in DWD's OI
  ! (exponential function by default)
  !---------------------------------------------------------------
  real(wp) ,intent(in) :: x     ! horizontal distance (secant on unit sphere)
  real(wp) ,intent(in) :: L     ! length scale
  real(wp)             :: corrh ! horizontal correlation

    real(wp) :: r

    if (rhcexp) then
      !-------------------------------------------
      ! convert from secant to radiant if required
      !-------------------------------------------
      if (sec2rad) then
        r = 2._wp * asin (0.5_wp * x)
      else
        r = x
      endif
      !-------------------------------
      ! scale argument by length scale
      !-------------------------------
      r = r * rearth / L
      corrh = exp (-0.5_wp * r**2)
    else
      corrh = corr (x,L)
    endif

  end function corrh
!------------------------------------------------------------------------------
  function corr (x,L)
  !---------------------------------------
  ! return horizontal correlation function
  !---------------------------------------
  real(wp) ,intent(in) :: x    ! horizontal distance (secant on unit sphere)
  real(wp) ,intent(in) :: L    ! length scale
  real(wp)             :: corr ! horizontal correlation

    real(wp) :: r
    !--------------------------------------------
    ! Initialize Covariance model if not yet done
    !--------------------------------------------
    if (.not.lset_h_corr) &
      call finish ('corr','not initialized, call set_h_corr first')
    !-------------------------------------------
    ! convert from secant to radiant if required
    !-------------------------------------------
    if (sec2rad) then
      r = 2._wp * asin (0.5_wp * x)
    else
      r = x
    endif
    !-------------------------------
    ! scale argument by length scale
    !-------------------------------
    r = r * L1r / L
    select case (model)
    !---------------------------------------------------------------
    ! horizontal correlation as used in DWD's OI
    ! (Bessel function expansion as defined in RM1, paragraph 2.3.2)
    !---------------------------------------------------------------
    case ('dwd-oi')
      corr = bess (r)
    !-------------------------------------------------------------------------
    ! compactly supported 5th order piecewise rational functions C0(r,a=1/2,c)
    ! cf. Gaspari+Cohn QJRMS (1999),125,723-757
    !-------------------------------------------------------------------------
    case ('c5pr_12')
      if (r < x_cuts) then
        ir = int (r)
        select case (ir)
        case (0 )
          corr = c00 + r *        r *( c02 + r *( c03 + r *( c04 + r*c05)))
        case (1 )
          corr = c10 + r *( c11 + r *( c12 + r *( c13 + r *( c14 + r*c15)))) &
               + c1x / r
        case (2:)
          corr = 0._wp
        case default
          call finish('corr','negative r')
        end select
      else
        corr = 0._wp
      endif
    !-------------------------------------------------------------------------
    ! compactly supported 5th order piecewise rational functions C0(r,a=inf,c)
    ! cf. Gaspari+Cohn QJRMS (1999),125,723-757
    !-------------------------------------------------------------------------
    case ('c5pr_inf')
      if (r < x_cuts) then
        ir = int (2*r)
        select case (ir)
        case (0 )
          corr = f00 + r *        r *( f02 + r *( f03 + r *( f04 + r*f05)))
        case (1 )
          corr = f10 + r *( f11 + r *( f12 + r *(       r *( f14 + r*f15)))) &
               + f1x / r
        case (2 )
          corr = f20 + r *( f21 + r *( f22 + r *( f23 + r *( f24 + r*f25)))) &
               + f2x / r
        case (3 )
          corr = f30 + r *( f31 + r *( f32 + r *( f33 + r *( f34 + r*f35)))) &
               + f3x / r
        case (4:)
          corr = 0._wp
        case default
          call finish('corr','negative r')
        end select
      else
        corr = 0._wp
      endif
    !----------------------------------------------------------
    ! compactly supporting functions approximating SOAR or TOAR
    ! cf. Gneiting QJRMS (1999),125,2449-2464
    !----------------------------------------------------------
    case ('csxoar1')
      if (r < x_cuts) then
        corr = (1._wp+(a+1)*r)*(1-r)**(a+1._wp)
      else
        corr = 0._wp
      endif
    case ('csxoar2')
      if (r < x_cuts) then
        corr = (gna + gnb * r + gnc * r*r) * (1-r)**(a+2)
      else
        corr = 0._wp
      endif

    !-----------------------------------------------------
    ! for development+debugging allow for different models
    !-----------------------------------------------------
    case ('gauss')
      corr = exp(-r**2/2._wp)
    case ('delta')
      corr = 0._wp
      if(r==0._wp) corr = 1._wp
    case ('1')
      corr = 1._wp
    case default
      call finish ('dcorr','invalid correlation model: '//model)
    end select
    rn = r
  end function corr

!------------------------------------------------------------------------------
  function bess (r)
  !---------------------------------------------------------------
  ! Bessel function expansion for DWD horizontal correlation model
  !---------------------------------------------------------------
  real(wp), intent(in) :: r    ! scaled argument, 1-> zero gradient
  real(wp)             :: bess ! result (correlation function)
    integer :: i
    if (r < x_cuts) then
      if (r /= rn) then ! skip if argument is the same as last time
        do i=1,na
          J0(i) = bessj0(r*rkn(i))
          J1(i) = bessj1(r*rkn(i))
        end do
      endif
      bess = a_dwd(0) + sum (a_dwd(1:na) * J0(1:na))
    else
      J0   = 0._wp
      J1   = 0._wp
      bess = 0._wp
    endif
  end function bess
!------------------------------------------------------------------------------

  subroutine dcorr (d1chh, d1chho1, d2chh, d1chx)
  real(wp) ,intent(out) :: d1chh   ! first  derivative of chh
  real(wp) ,intent(out) :: d2chh   ! second derivative of chh
  real(wp) ,intent(out) :: d1chho1 ! first  derivative of chh / r.
  real(wp) ,intent(out) :: d1chx   ! first  derivative, truncated expansion

    real(wp) :: sum1    ! bessel function expansion for wind correlations
    real(wp) :: sum2    ! bessel function expansion for wind correlations
    real(wp) :: r       ! argument (scaled distance, stored in rn)
    real(wp) :: w       ! weight for truncated expansion
    integer  :: i       ! index variable

    if (rn >= x_cuts) then
      d1chh   = 0._wp
      d1chho1 = 0._wp
      d2chh   = 0._wp
      d1chx   = 0._wp
    else
      select case (model)
      !------------------------------------
      ! calculate Bessel function expansion
      !------------------------------------
      case ('dwd-oi')
        sum1   =   sum (a_dwd(1:na) * rkn(1:na)    * J1(1:na))
        sum2   =   sum (a_dwd(1:na) * rkn(1:na)**2 * J0(1:na))
        d1chh  = - sum1
        if (rn /= 0._wp) then
          d1chho1 = d1chh / rn
        else
          d1chho1 = -1._wp / L2
        endif
        d2chh  = - sum2 - d1chho1
      !-----------------------------------------------------------
      ! compactly supported 5th order piecewise rational functions
      ! C0(r,a=1/2,c) cf. Gaspari+Cohn QJRMS (1999),125,723-757
      !-----------------------------------------------------------
      case ('c5pr_12')
        r = rn
        select case (ir)
        case (0 )
         d1chh =     r*( 2*c02 + r*( 3*c03 + r*(  4*c04 + r* 5*c05)))
         d2chh =         2*c02 + r*( 6*c03 + r*( 12*c04 + r*20*c05))
        case (1 )
         d1chh = c11+ r*(2*c12 + r*( 3*c13 + r*(  4*c14 + r* 5*c15)))-c1x/r**2
         d2chh =         2*c12 + r*( 6*c13 + r*( 12*c14 + r*20*c15))+2*c1x/r**3
        case (2:)
         d1chh = 0._wp
         d2chh = 0._wp
        case default
         call finish('corr','negative r')
        end select
      !-----------------------------------------------------------
      ! compactly supported 5th order piecewise rational functions
      ! C0(r,a=inf,c) cf. Gaspari+Cohn QJRMS (1999),125,723-757
      !-----------------------------------------------------------
      case ('c5pr_inf')
        r = rn
        select case (ir)
        case (0 )
         d1chh =     r*( 2*f02 + r*( 3*f03 + r*(  4*f04 + r* 5*f05)))
         d2chh =         2*f02 + r*( 6*f03 + r*( 12*f04 + r*20*f05))
        case (1 )
         d1chh = f11+ r*(2*f12 + r*(         r*(  4*f14 + r* 5*f15)))-f1x/r**2
         d2chh =         2*f12 + r*(         r*( 12*f14 + r*20*f15))+2*f1x/r**3
        case (2 )
         d1chh = f21+ r*(2*f22 + r*( 3*f23 + r*(  4*f24 + r* 5*f25)))-f2x/r**2
         d2chh =         2*f22 + r*( 6*f23 + r*( 12*f24 + r*20*f25))+2*f2x/r**3
        case (3 )
         d1chh = f31+ r*(2*f32 + r*( 3*f33 + r*(  4*f34 + r* 5*f35)))-f3x/r**2
         d2chh =         2*f32 + r*( 6*f33 + r*( 12*f34 + r*20*f35))+2*f3x/r**3
        case (4:)
         d1chh = 0._wp
         d2chh = 0._wp
        case default
         call finish('corr','negative r')
        end select
      !----------------------------------------------------------
      ! compactly supporting functions approximating SOAR or TOAR
      ! cf. Gneiting QJRMS (1999),125,2449-2464
      !----------------------------------------------------------
      case ('csxoar1')
        r = rn
        d1chh = (-2-a)*r*(a+1)*(1-r)**(a)
        d2chh = (-2-a)*  (a+1)*(1-r)**(a-1)*(1-(1+a)*r)
      case ('csxoar2')
        r = rn
        d1chh = (gne + gnf * r + gng * r*r) * (1-r)**(a+1)
        d2chh = (gni + gnj * r + gnk * r*r) * (1-r)**(a)

      case ('1')
        d1chh   = 0._wp
        d1chho1 = 0._wp
        d2chh   = 0._wp
      case default
        call finish ('dcorr','invalid correlation model: '//model)
      end select
      !-----------------------------------------------------------
      ! truncated expansion for height-wind-correlations in DWD OI
      !-----------------------------------------------------------
      select case (model)
      case ('dwd-oi')
        if (nw >= na) then
          d1chx = d1chh
        else
          sum1 = 0._wp
          do i=1,na
            w = min (nw - i, 1._wp)
            sum1 = sum1 + w * a_dwd(i) * rkn(i) * J1(i)
          end do
          d1chx = - sum1
        endif
      case default
        d1chx = d1chh
      end select
      !----------------------------------
      ! for r=0: d/dr corr is not defined
      !----------------------------------
      if (rn /= 0._wp) then
        d1chho1 = d1chh / rn
      else
        d1chho1 = -1._wp / L2
      endif
      !--------------------------------------
      ! rescale derivatives with length scale
      !--------------------------------------
      d1chh   = d1chh   * L1
      d1chho1 = d1chho1 * L2
      d2chh   = d2chh   * L2
      d1chx   = d1chx   * L1
      !-------------------------------------------
      ! convert from secant to radiant if required
      !-------------------------------------------
!      if (.not. sec2rad) then
!        cosr2   = sqrt(1._wp-(0.5_wp*rn)**2)
!        d1chh   = d1chh   * cosr2
!        d1chho1 = d1chho1 * cosr2
!        d2chh   = d2chh   * cosr2 **2
!      endif
    endif

  end subroutine dcorr
!==============================================================================
  subroutine set_h_corr
  !------------------------------------------------------
  ! initializes horizontal correlation model coefficients
  !------------------------------------------------------
    real(wp) :: x1, x2, x0 ! used to bracked zero derivatives of J_0
    real(wp) :: x,y        ! used to print simple correlation funct. plot
    real(wp) :: xold       ! used to bracket first zero of correlation fktn
    real(wp) :: d1chh      ! first  derivative of chh
    real(wp) :: d1chx      ! first  derivative, truncated expansion (h-w cor.)
    real(wp) :: d2chh      ! second derivative of chh
    real(wp) :: d1chho1    ! first  derivative of chh / rn
    integer  :: i
    integer  :: nx = 110! number of points to plot in x (simple plot of ..
    integer  :: ny = 50 ! number of points to plot in y  .. correlation funct.)
    real(wp) :: xend    ! end of plot
    real(wp) :: a_norm  ! normalization faktor so that corr(0) = 1
    real(wp) :: sumpm, suma
    real(wp) :: sump1, summ1, sump2, summ2

    sec2rad     = .false.
    x_cuts      = 100._wp
    lset_h_corr = .true.

    select case (model)
    !======================================================
    ! Bessel function expansion similar to former OI scheme
    !======================================================
    case ('dwd-oi')
      !----------------------------------------------------------------
      ! calculate locations of zero derivatives of Bessel functions J_0
      ! (the rkn are the k_n in RM1, Equ. (2.3.1)),
      !----------------------------------------------------------------
      rkn(0)=0
      x0 = 0._wp
      do i=1,mkn
        x1     = x0+1._wp
        x2     = x0+5._wp
        x0   =  zbrent (d_bessj0,x1,x2,1.e-15_wp)
        rkn(i) = x0
      end do
      !---------------------------------------------------------
      ! modify the coefficients so that d2 corr / dr2 (r=D) == 0
      !---------------------------------------------------------
      if (fixbc > 0) then
        L1r      = 1._wp
        L1       = 1._wp
        L2       = 1._wp
        a_dwd(0) = - bess (1._wp)
        rn = 1._wp
        call dcorr (d1chh, d1chho1, d2chh, d1chx)
        sumpm = 0._wp
        sump1 = 0._wp
        summ1 = 0._wp
        sump2 = 0._wp
        summ2 = 0._wp
        do i=1,na
          suma = a_dwd(i) * rkn(i)    * J1(i) &
               - a_dwd(i) * rkn(i)**2 * J0(i)
          sumpm = sumpm + suma
          if (suma > 0._wp) then
            sump1 = sump1 + suma
            sump2 = sump2 + suma /a_dwd(i)
          else
            summ1 = summ1 - suma
            summ2 = summ2 - suma /a_dwd(i)
          endif
        end do
        select case (fixbc)
        case (1)
          !------------------------------------------------------
          ! modify coeffs by adding an offset to all coefficients
          !------------------------------------------------------
          a_dwd(1:na:2) = a_dwd(1:na:2) - sumpm / (sump2+summ2)
          a_dwd(2:na:2) = a_dwd(2:na:2) + sumpm / (sump2+summ2)
        case (2)
          !------------------------------------------------------------------
          ! modify coeffs by rescaling the coeffs. with negative contribution
          !------------------------------------------------------------------
          a_dwd(2:na:2) = a_dwd(2:na:2) * sump1/summ1
        case default
          call finish ('set_h_corr',&
                       'invalid value for namelist parameter "fixbc"')
        end select
      endif
      !--------------------------------------
      ! calculate a0 and normalisation factor
      !   so that corr(r=D) == 0
      !   and     corr(r=0) == 1
      !--------------------------------------
      rn       = - 1._wp
      L1r      =   1._wp
      a_dwd(0) =   0
      a_dwd(0) = - bess (1._wp)
      !----------------------------------------
      ! calculate length scale L1 for D=1
      ! normalization factor
      !----------------------------------------
      a_norm = 1._wp / bess (0._wp)
      a_dwd  = a_norm * a_dwd
      L2   = 2._wp / sum (a_dwd(1:na) * rkn(1:na)**2)
      L1   = sqrt(L2)
      L1r  = L1
      edL1 = 1._wp/L1
      !-----------------------
      ! write coefficients k_n
      !-----------------------
      if (dace% lpio) then
        print *
        print *,'length scale L1  for D=1:', L1
        print *,'D for length scale L1=1 :', edL1
        print *,'normalization factor    :', a_norm
        print *
        print *, 'Bessel Function Expansion for horizontal correlations:'
        print *
        print *, 'i, coefficient A_i, normalized, coefficients rkn:'
        do i=0,na
          write(6,'(i2,3f20.15)') i,a_dwd(i)/a_norm,a_dwd(i),rkn(i)
        end do
        !------------------------------------
        ! simple plot of correlation function
        !------------------------------------
        print *
        print *,'Horizontal Correlation Function:'
        print *,' corr(x=0)=',bess (0._wp)
        print *,' corr(x=L)=',bess (L1   )
        print *,' corr(x=D)=',bess (1._wp)
        print *
        print *, ' x/L,    x,       C_hor(x/L):'
      endif
      !-------------------------------------------
      ! bracket first zero in correlation function
      !-------------------------------------------
      xold   = 0._wp
      x_zero = 1._wp
      do i=0,nx
        x = 1._wp/nx * i
        y = bess (x)
        if (dace% lpio) &
          write(6,'(3f7.3,1x,a)') x, x*L1, y, repeat('=',max(0,nint(ny*y)))
        if (x_zero==1._wp .and. y<0._wp)x_zero=zbrent(bess,xold,x,1.e-15_wp)
        xold = x
      end do
      xend   = 1._wp
    !=========================================================================
    ! compactly supported 5th order piecewise rational functions C0(r,a=1/2,c)
    ! cf. Gaspari+Cohn QJRMS (1999),125,723-757
    !=========================================================================
    case ('c5pr_12')
      L2 = -1._wp / (2*c02)
      xend   = 2._wp
      x_zero = 2._wp
    case ('c5pr_inf')
      L2 = -1._wp / (2*f02)
      xend   = 2._wp
      x_zero = 2._wp
    !----------------------------------------------------------
    ! compactly supporting functions approximating SOAR or TOAR
    ! cf. Gneiting QJRMS (1999),125,2449-2464
    !----------------------------------------------------------
    case ('csxoar1')
      L2     = 1._wp / ((2+a)*(a+1))
      xend   = 1._wp
      x_zero = 1._wp
    case ('csxoar2')
      gna = 1
      gnb = a+2
      gnc = ((a+2)**2-1)/3
      gne =     gnb - (a+2) * gna
      gnf = 2 * gnc - (a+3) * gnb
      gng =         - (a+4) * gnc
      gni =     gnf - (a+1) * gne
      gnj = 2 * gng - (a+2) * gnf
      gnk =         - (a+3) * gng
      L2     = -1._wp / gni
      xend   =  1._wp
      x_zero =  1._wp

    case ('1')
      L2 =  1._wp
      xend   = 5._wp
    case default
      call finish ('set_h_corr','invalid correlation model: '//model)
    end select
    !--------------------------------
    ! derive some quantities from L^2
    !--------------------------------
    L1     = sqrt(L2)
    L1r    = L1
    edL1   = 1._wp/L1
    x_cut  = xend /L1
    x_zero = x_zero/L1
    !---------------------------------------------------------------------
    ! write a file 'h_corr_xxxx' with correlation function and derivatives
    !---------------------------------------------------------------------
    if (dace% lpio) then
      open (1,file=path_file (aux, 'h_corr_'//model))
      print *
      write(1,*)'# cutoff radius x_cut    ,f :', x_cut, corr(x_cut,1._wp)
      if (cut >0._wp) write(1,*)'# cutoff radius rescaled    :', cut
      write(1,*)'# L1                        :', L1
      write(1,*)'# L1^2                      :', L2
      write(1,*)'# 1/L1                      :', edL1
      write(1,*)'#'
      write(1,*) &
      '# r/L, r, C_hor(r), d/dr C_hor, 1/r d/dr C_hor, d2/d2r C_hor:'
      write(6,*) &
      '  r/L, r, C_hor(r), d/dr C_hor, 1/r d/dr C_hor, d2/d2r C_hor:'
      do i=0,nx
        x = xend/nx * i
        y = corr(x, L1)
        call dcorr (d1chh, d1chho1, d2chh, d1chx)
        write(1,'(6f15.11)') x, x/L1, y, d1chh, d1chho1, d2chh
        if (i==0 .or. i==nx) write(6,'(6f15.11)') &
          x, x/L1, y, d1chh, d1chho1, d2chh
      end do
!     close (1)
      close (1, status="delete")    ! Don't keep that file
      print *
    endif
    !----------------------------------------------------
    ! finally switch on conversion from secant to radiant
    ! scaling from unit sphere to earth radius
    !----------------------------------------------------
    sec2rad             = (model=='dwd-oi')
    if(irad>=0) sec2rad = (irad == 1)
    L1r     = L1 * rearth
    x_cut   = xend /L1
    if (cut >0._wp .and. cut < x_cut) then
      if (dace% lpio) &
        write(6,'(/a,2f5.2/)') '# rescaling x_cut :',x_cut,cut
      x_cut   = cut
    endif
    x_cuts   = x_cut  * L1
    x_cutdr  = x_cut / rearth
if(dace% lpio) then
write(6,*) 'set_h_corr: x_cut   = ', x_cut
!write(6,*) 'set_h_corr: x_cutdr = ', x_cutdr
!write(6,*) 'set_h_corr: x_cuts  = ', x_cuts
write(6,*) 'set_h_corr: L1      = ', L1
!write(6,*) 'set_h_corr: L1^2    = ', L2
!write(6,*) 'set_h_corr: 1/L1    = ', edL1
write(6,*) 'set_h_corr: sec2rad = ', sec2rad
write(6,*) 'set_h_corr: rhcexp  = ', rhcexp
write(6,*) 'set_h_corr: fixbc   = ', fixbc
write(6,*) 'set_h_corr: na      = ', na
write(6,*) 'set_h_corr: nw      = ', nw
endif
  end subroutine set_h_corr
!==============================================================================
  function d_bessj0(x)
  !------------------
  ! Derivative of J_0
  !------------------
  real(wp) ,intent(in) ::  x
  real(wp)             ::  d_bessj0
    d_bessj0 = - bessj1 (x)
  end function d_bessj0
!==============================================================================
  subroutine construct_rowcol (rc, n, nh)
  type(t_rowcol) ,intent(out) :: rc   ! row or column to set
  integer        ,intent(in)  :: n    ! number of observations in a row
  integer        ,intent(in)  :: nh   ! number of spots in a row
  !-------------------------
  ! setup allocatable arrays
  !-------------------------
!   integer :: k
    rc% n     = n
    rc% nh    = nh
    rc% in    = 0
    rc% inh   = 0
    allocate (rc% cols  (n))
    allocate (rc% spots (nh))
!   forall (k=1:nh) rc% spots (k)% vm => NULL()
    rc% spots% ih = 0
  end subroutine construct_rowcol
!------------------------------------------------------------------------------
  subroutine zero_rowcol (rc)
  type(t_rowcol) ,intent(inout) :: rc   ! row or column to set
    rc% in        = 0
    rc% inh       = 0
    rc% spots% ih = 0
  end subroutine zero_rowcol
!------------------------------------------------------------------------------
  subroutine destruct_rowcol (rc)
  type(t_rowcol) ,intent(inout) :: rc   ! row or column to free
  !-------------------------------------
  ! deallocate module allocatable arrays
  !-------------------------------------
    integer :: k
    if (associated (rc% spots)) then
      do k=1, size(rc% spots)
        if (allocated (rc% spots(k)% vm)) deallocate (rc% spots(k)% vm)
      end do
      deallocate (rc% spots)
    endif
    if (associated (rc% cols )) deallocate (rc% cols)
  end subroutine destruct_rowcol
!------------------------------------------------------------------------------
  subroutine destruct_rowcols (rc)
  type(t_rowcol) ,intent(inout) :: rc(:)   ! row or column to free
  !-------------------------------------
  ! deallocate module allocatable arrays
  !-------------------------------------
    integer :: i
    do i = 1, size(rc)
      call destruct (rc(i))
    end do
  end subroutine destruct_rowcols
!------------------------------------------------------------------------------
  subroutine init_spot (rc, ih, coo, it, z, lcl, nwv, err)
  type(t_rowcol),intent(inout)         :: rc
  integer       ,intent(in)            :: ih     ! spot index within the box
  type(t_coord) ,intent(in)            :: coo    ! coordinates
  integer       ,intent(in)            :: it(:)  ! observation type identifiers
  real(wp)      ,intent(in)            :: z (:)  ! vertical coordinates (ln p)
  integer       ,intent(in)  ,optional :: lcl    ! clim. or fc error
  integer       ,intent(in)  ,optional :: nwv    ! vertical interpolation flag
  real(wp)      ,intent(out) ,optional :: err(:) ! error (std deviation)
  !============================================
  ! precalculate coefficients for a spot,
  ! i.e. for all observations in a model column
  !============================================
    !----------------
    ! local variables
    !----------------
    integer                   :: n       ! number of levels
    type(t_dwd_spot) ,pointer :: spt     ! spot data type
#if !defined (__NEC__)
    type(t_dwd_col)  ,pointer :: c       ! level data type
    type(t_vic)      ,pointer :: vi      ! vertical interpolation data type
#endif
    type(t_ccf)      ,pointer :: ccf     ! compressed climatology file content
    real(wp)                  :: w,w1    ! latitude weights
    real(wp)                  :: w00,w01 ! latitude/longitude weights
    real(wp)                  :: w10,w11 ! latitude/longitude weights
    real(wp)                  :: emodll  ! latitude/longitude weight
    integer                   :: j,j1    ! latitude indices
    integer                   :: k,k2,l  ! level indices
    integer                   :: nl      ! number of levels
    integer                   :: i,i1    ! longitude/polynom coefficient index
    integer                   :: ix1     ! Auxiliary index
    integer                   :: nlnp    ! number of same lnp
    real(wp)                  :: lnptop  ! ln(cutoff pressure) (  10hPa)
    real(wp)                  :: lnpbot  ! ln(cutoff pressure) (1100hPa)
    real(wp)                  :: xvz,xqz ! transformation coeff. for lnp<lnptop
    integer                   :: lf      ! forecast/climatol. errors
    real(wp)                  :: dlo     ! longitude in the range 0...360
    integer                   :: idx(9,4)! indices    from grid_indices
    integer                   :: np      ! no. points from grid_indices
    real(wp)                  :: wi (9)  ! weights    from grid_indices
    real(wp),       parameter :: zedem1 = 2*e / (e-1)
    integer,        parameter :: n3_ = max (INT_H, INT_RH, INT_CHI)
    integer                   :: lnwv    ! local vertical interpolation flag
    integer                   :: kl      ! last K level for RTTOV
#if defined (__NEC__)
    integer                   :: k2m
#endif
    !-----------------------------------
    ! Auxiliary arrays for vectorization
    !-----------------------------------
    real(wp),     allocatable :: lnp (:) ! ln (p) (p in Pa)
    real(wp),     allocatable :: lnpc(:) ! ln (p) (p in Pa) (constrained)
    real(wp),     allocatable :: lnpf(:) ! lnpc reversed
                                         ! Mask: adjacent levels with same p:
    logical                   :: samep(size(it)-1)
    integer,      allocatable :: klev(:)        ! Index of k'th p level
    integer                   :: llev(size(it)) ! Index to p level
    type(t_vic),  allocatable :: vic (:)        ! Temp. array for set_vertinpc
#if !defined (__NEC__)
    target                    :: vic            ! Temp. array for set_vertinpc
#endif
    real(wp)     ,dimension(:) &
                 ,allocatable :: eh, ev, eq, ehz, et, &! Errors
                                 xh, xv, xq, xhz       ! Derivatives of errors
    real(wp)                  :: vm(covm% nz, INT_TV)  ! Vertical modes
    !-----------------------------
    ! Vertical covariance matrices
    !-----------------------------
#if defined (__NEC__)
    ! Transpose first two dimensions of vcm
#define vcm(a,b,c) vcm_t(b,a,c)
    ! Take care of bank conflicts on NEC SX: use odd leading dimensions
    real(wp) :: vcm_t ((covm% nz/2)*2+1, (covm% nz/2)*2+1, n3_)
#else
    ! Transpose first two dimensions of temporary vcm for better vectorization
#define vcm(a,b,c) vcm_t(b,a,c)
    real(wp) :: vcm   ( covm% nz       ,  covm% nz,        n3_)
#endif

    if (.not. linit_fg_cov) call finish ('init_spot','init_fg_cov not called')
    lf  = ifc;  if (present (lcl)) lf  = lcl
    lnwv =                         -1
    if (present(nwv)) lnwv =       nwv
    if (lnwv < 1)     lnwv = covm% nwv
    !-----------------
    ! set spot pointer
    !-----------------
    n = size(it)
    spt => rc% spots(ih)
    !-------------------------------
    ! fill data type if not yet done
    !-------------------------------
    if (spt% ih == 0) then
      spt% ih      =  ih
      rc % inh     =  ih
      !-------------------
      ! set column pointer
      !-------------------
      spt% i       =  rc% in
      spt% n       =  n
      rc%  in      =  rc% in + n
      spt% col     => rc% cols (spt%i+1: spt%i+n)
      !---------------------------------
      ! store unit vectors on the sphere
      !---------------------------------
      spt% x       =  coo% x    ! (location)
      spt% du      =  coo% du   ! (u-direction)
      spt% dv      =  coo% dv   ! (v-direction)
      !---------------------------------------------------------------
      ! indices and weights for latitude dependent correlations
      ! 6degree x 6degree lat-lon grid for analysis error is hardcoded
      !---------------------------------------------------------------
      spt% wlat    =  min (30._wp, max (1._wp, (coo% dlat+87._wp)/6._wp+1))
      spt% jlat    =  int(spt% wlat)
      spt% wlat    =  spt% wlat - spt% jlat
      spt% jlat1   =  min (30, spt% jlat+1)
      dlo          =  coo% dlon
      if (dlo < 0._wp) dlo = dlo + 360._wp
      spt% wlon    =  dlo / 6._wp
      spt% ilon    =  int (spt% wlon)
      spt% wlon    =  spt% wlon - spt% ilon
      spt% ilon    =  mod (spt% ilon, 60) + 1
      spt% ilon1   =  mod (spt% ilon, 60) + 1
      !-------------------------------------------------
      ! linear interpolation of correlation coefficients
      !-------------------------------------------------
      w   = spt% wlat; w1 = 1._wp - w
      j   = spt% jlat; j1 = spt% jlat1

      ccf => nfcoun
      if (lf == IFC_CL) ccf => nclerr  ! climatological or forecast

      spt% L_h = w * ccf% b  (  j1, fg_month) + w1 * ccf% b  (  j, fg_month)
      spt% L_q = w * ccf% brh(  j1, fg_month) + w1 * ccf% brh(  j, fg_month)
      spt% mu  = w * ccf% mu (  j1, fg_month) + w1 * ccf% mu (  j, fg_month)
      if (abs (rmu) <= 1._wp) spt% mu = rmu

      spt% cch = w * ccf% czc(  j1, fg_month) + w1 * ccf% czc(  j, fg_month)
      spt% ccv = w * ccf% cuc(  j1, fg_month) + w1 * ccf% cuc(  j, fg_month)
      spt% ccq = w * ccf% crc(  j1, fg_month) + w1 * ccf% crc(  j, fg_month)

      if (covm% valid <= 1) then
        spt% czh = w * ccf% cz (:,j1, fg_month) + w1 * ccf% cz (:,j, fg_month)
        spt% cxh = w * ccf% cx (:,j1, fg_month) + w1 * ccf% cx (:,j, fg_month)

        spt% czv = w * ccf% cu (:,j1, fg_month) + w1 * ccf% cu (:,j, fg_month)
        spt% cxv = w * ccf% cy (:,j1, fg_month) + w1 * ccf% cy (:,j, fg_month)

        spt% czq = w * ccf% cr (:,j1, fg_month) + w1 * ccf% cr (:,j, fg_month)
        spt% cxq = w * ccf% ca (:,j1, fg_month) + w1 * ccf% ca (:,j, fg_month)
      end if
      !------------------------------------------------------------
      ! rescale forecast error depending on previous analysis error
      !------------------------------------------------------------
      if (lf == IFC_AFC) then
        if (samegrid) then
          i   = spt% ilon; i1 = spt% ilon1
          w00 = spt% wlon * w ; w10 = (1._wp - spt% wlon) * w
          w01 = spt% wlon * w1; w11 = (1._wp - spt% wlon) * w1
          emodll = w00 * emodsm(i1,j1) + w10 * emodsm(i ,j1) &
                 + w01 * emodsm(i1,j ) + w11 * emodsm(i ,j )
        else
          call Grid_Indices       &
                   (coo% dlon,    &! <-- geodetic longitude
                    coo% dlat,    &! <-- geodetic latitude
                    ba_err% grid, &! <-- grid data type
                    idx,          &! --> Grid point indices [Point, index]
                    wi,           &! --> Weight
                    np)            ! --> number of points returned
          emodll = 0._wp
          do i=1,np
            emodll = emodll + wi(i) * emodsm(idx(i,1),idx(i,2))
          end do
        endif
      else
        emodll = scale
      endif
      spt% emod = emodll

      if (emodll <= 0._wp) then
         ! This should never happen!
         write (0,*) "### init_spot: Lat,Lon,emodll:", &
              coo% dlat, coo% dlon, emodll
         write (0,*) "### lf, IFC_AFC, samegrid, scale=", &
              lf, IFC_AFC, samegrid, scale
         if (lf == IFC_AFC) then
            if (samegrid) then
               write (0,*) "###", i, i1, j, j1, w, w1
               write (0,*) "###", w00, w10, w01, w11
            else
               write (0,*) "#### np, wi =", np, wi(1:np)
               write (0,*) "#### idx(:,1:2) =", idx(:np,1:2)
            end if
         end if
         write (0,*) "### spt% wlat, spt% wlon =", spt% wlat, spt% wlon
         call finish ("init_spot","emodll <= 0")
      end if

      if (covm% valid <= 1) then
         spt% czh = spt% czh * emodll
         spt% czv = spt% czv * emodll
         !---------------
         ! constant error
         !---------------
         if (e_h >= 0._wp) then
            spt% czh = 0; spt% czh(1) = e_h
         endif
         if (e_q >= 0._wp) then
            spt% czq = 0; spt% czq(1) = e_q
         endif
         if (e_v >= 0._wp) then
            spt% czv = 0; spt% czv(1) = e_v
         endif
         !----------------------------------------------------
         ! constant, linear vertical coordinate transformation
         !----------------------------------------------------
         if (t_h /= 0._wp) then
            spt% cxh = 0; spt% cxh(2) = t_h
         endif
         if (t_q /= 0._wp) then
            spt% cxq = 0; spt% cxq(2) = t_q
         endif
         if (t_x /= 0._wp) then
            spt% cxv = 0; spt% cxv(2) = t_x
         endif
      end if
      !-----------------------------------------------
      ! constant vertical correlation tuning parameter
      !-----------------------------------------------
      if (a_h >  0._wp) spt% cch = a_h
      if (a_q >  0._wp) spt% ccq = a_q
      if (a_x >  0._wp) spt% ccv = a_x
      !---------------------------------
      ! for testing set:
      ! constant horizontal length scale
      !---------------------------------
      if (L_h >  0._wp) spt% L_h = L_h
      if (L_q >  0._wp) spt% L_q = L_q
      !-----------------------------------------------------
      ! transformed horizontal coordinates (deformed sphere)
      !-----------------------------------------------------
      if (thc) then
        spt% xt = spt% x * Lh_max / spt% L_h
      else
        spt% xt = spt% x
      endif
      !-----------------------------------------------------------
      ! fill column with observation type and z (ln(p)) coordinate
      !-----------------------------------------------------------
      spt% col% it = it
      spt% col% z  = z
      !------------------------------------------
      ! set horizontal interpolation coefficients
      ! for new (NMC-fitted) covariance model
      !------------------------------------------
      if (covm% valid > 0) call set_horinpc (spt% h, coo% dlat, covm)
      !----------------------------------------------
      ! search first observation on a level
      ! determine number of observations on the level
      ! 'ior' observation types on the level
      !----------------------------------------------
FTRACE_BEGIN("init_spot:levels")
      !-----------------------------------------------
      ! Partially vectorized version for NEC.
      ! Compare adjacent levels, count distinct levels
      !-----------------------------------------------
      samep(1:n-1) = (spt% col(1:n-1)% z == spt% col(2:n)% z)
      nl           = n - count (samep)
      !-------------
      ! Set defaults
      !-------------
!NEC$ ivdep
      do k = 1, n
         select case (it(k))
         case default
            spt% col(k)% itl = it(k)
         case (OBS_U, OBS_V)
            spt% col(k)% itl = ior (OBS_U, OBS_V)       ! Always join wind obs.
         end select
      end do
      spt% col(1:n)% i = 1
      !---------------------------------------
      ! Join observation types on same level.
      ! Set auxiliary arrays with indices to
      ! first obs on level (klev),
      ! mapping of obs to level number (llev).
      !---------------------------------------
      allocate (klev(nl))
      nlnp    = 1
      l       = nl
      kl      = 0
      if (n>0) then
        llev(n) = l
        spt% col(n)% i = 1
        do k = n-1,1,-1
          if (kl == 0 .and. it(k+1) == OBS_RH) kl = l
          if (samep(k)) then
            nlnp             = nlnp + 1
            spt% col(k)% i   = nlnp
            spt% col(k)% itl = ior (it(k), spt% col(k+1)% itl)
          else
            nlnp             = 1
            klev(l)          = k + 1    ! Set index to first obs on level
            l                = l - 1
          end if
          llev(k) = l                    ! Map obs # to level number
        end do
        klev(1) = 1
      endif
      if (lnwv /= 3) kl = nl
      !++++++++++++++++++++++++++++++++++++++
      ! Consistency check, may be disabled...
      !++++++++++++++++++++++++++++++++++++++
      do l = 1, nl
         if (llev(klev(l)) /= l) call finish ("init_spot","llev(klev(l)) /= l")
      end do
FTRACE_END  ("init_spot:levels")
      if (covm% valid > 0) then
        !----------------------------------------
        ! set vertical interpolation coefficients
        ! for new (NMC-fitted) covariance model
        !----------------------------------------
FTRACE_BEGIN("init_spot:vertinpc")
        allocate   (vic(nl))
        call set_vertinpc (vic, spt% col(klev(1:kl))% z, covm% logp, lnwv)
        spt% col(klev(1:nl))% vi = vic(1:nl)
FTRACE_END  ("init_spot:vertinpc")
      end if
      if (covm% valid > 1) then
        if (.not. allocated (spt% vm)) &
          allocate (spt% vm(covm% nz, INT_TV, nl))
        !---------------------------------------
        ! calculate vertical covariance matrices
        !---------------------------------------
!FTRACE_BEGIN("init_spot:vcm")
!$omp parallel do private(k,l) schedule(static)
!NEC$ shortloop
        do k=1,covm% nz
!NEC$ shortloop
        do l=1,covm% nz
          vcm (l,k,INT_H)   = spt%h%w(1) * covm% sqcvh (l,k,spt%h%ix) &
                            + spt%h%w(2) * covm% sqcvh (l,k,spt%h%ix+1)
          vcm (l,k,INT_RH)  = spt%h%w(1) * covm% sqcvq (l,k,spt%h%ix) &
                            + spt%h%w(2) * covm% sqcvq (l,k,spt%h%ix+1)
          vcm (l,k,INT_CHI) = spt%h%w(1) * covm% sqcvv (l,k,spt%h%ix) &
                            + spt%h%w(2) * covm% sqcvv (l,k,spt%h%ix+1)
        end do
        end do
!$omp end parallel do
!FTRACE_END  ("init_spot:vcm")
        !---------------------------
        ! calculate 'vertical modes'
        !---------------------------
FTRACE_BEGIN("init_spot:vmodes")   !!! Most expensive region of init_spot
#if defined (__NEC__)
!$omp parallel do private(i,k2,l,   vm) schedule(static)
#else
!$omp parallel do private(i,k2,l,vi,vm) schedule(static)
#endif
        do l = 1, nl
#if defined (__NEC__)
#define   VI    vic(l)
#else
          vi => vic(l)
#endif
          vm(:,:) = 0._wp
!NEC$ outerloop_unroll(2)
          do k2 = 1, VI% g% n
!NEC$ shortloop
           do i = 1, covm% nz
            vm(i,INT_H)   = vm(i,INT_H)   &
                          + VI% wh(k2) * vcm(VI% g% i+k2-1,i,INT_H)
            vm(i,INT_TV)  = vm(i,INT_TV)  &
                          + VI% wt(k2) * vcm(VI% g% i+k2-1,i,INT_H)
            vm(i,INT_RH)  = vm(i,INT_RH)  &
                          + VI% wh(k2) * vcm(VI% g% i+k2-1,i,INT_RH)
            vm(i,INT_CHI) = vm(i,INT_CHI) &
                          + VI% wh(k2) * vcm(VI% g% i+k2-1,i,INT_CHI)
           end do
          end do
!NEC$ shortloop
          spt% vm(:,:,l) = vm(:,:)
        end do
!$omp end parallel do
#ifdef VI
#undef VI
#endif
#ifdef vcm
#undef vcm
#endif
FTRACE_END  ("init_spot:vmodes")
      end if
FTRACE_BEGIN("init_spot:errors_calc")
      !--------------------
      ! set a few constants
      !--------------------
      if (covm% valid > 1) then
         lnptop = covm% logp(1)
         lnpbot = covm% logp(covm% nz)
      else
         lnptop = log(100._wp*ptop_v)
         lnpbot = log(100._wp*pbot_v)
      end if
      allocate (eh (nl))
      allocate (ev (nl))
      allocate (eq (nl))
      allocate (ehz(nl))
      allocate (et (nl))
      allocate (xh (nl))
      allocate (xv (nl))
      allocate (xq (nl))
      allocate (xhz(nl))
      allocate (lnp (nl))
      allocate (lnpc(nl))
      eh  (1:nl) = 0._wp
      ev  (1:nl) = 0._wp
      eq  (1:nl) = 0._wp
      ehz (1:nl) = 0._wp
      et  (1:nl) = 0._wp
      lnp (1:nl) = spt% col(klev(1:nl))% z
      lnpc(1:nl) = max (min (lnp(1:nl), lnpbot), lnptop)
      if (covm% valid > 1) then
        !----------------------------------
        ! new (NMC-fitted) covariance model
        !----------------------------------
#if defined (__NEC__)
        k2m = 0
        do l = 1, nl
           if (vic(l)% g% n > k2m) k2m = vic(l)% g% n
        end do
#endif
        if (lf == IFC_CL .and. covm% lclim > 1) then

!FTRACE_BEGIN("init_spot:ifc_cl")
#if !defined (__NEC__)  /* original scalar code */
!NEC$ ivdep
         do l = 1, nl
             vi => vic(l)
           do k2 = 1, VI% g% n
#else                   /* manual loop switch for cross-grained NEC nfort */
           do k2 = 1, k2m
!NEC$ ivdep
         do l = 1, nl
#define      VI    vic(l)
             if (k2 > VI% g% n) cycle
#endif
             ix1   = VI% g% i+k2-1
             eh (l)= eh (l)                                                 &
                   + VI% wh(k2) * (spt%h%w(1) * covm% eh_cl (ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% eh_cl (ix1,spt%h%ix+1))
             et (l)= et (l)                                                 &
                   + VI% wh(k2) * (spt%h%w(1) * covm% et_cl (ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% et_cl (ix1,spt%h%ix+1))
             ev (l)= ev (l)                                                 &
                   + VI% wh(k2) * (spt%h%w(1) * covm% ev_cl (ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% ev_cl (ix1,spt%h%ix+1))
             eq (l)= eq (l)                                                 &
                   + VI% wh(k2) * (spt%h%w(1) * covm% erh_cl(ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% erh_cl(ix1,spt%h%ix+1))
             ehz(l)= ehz(l)                                                 &
                   + VI% wt(k2) * (spt%h%w(1) * covm% eh_cl (ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% eh_cl (ix1,spt%h%ix+1))
#ifdef VI
#undef VI
#endif
           end do
          end do
!FTRACE_END  ("init_spot:ifc_cl")

        else ! lf /= IFC_CL .or. covm% lclim <= 1

!FTRACE_BEGIN("init_spot:ifc_other")
#if !defined (__NEC__)  /* original scalar code */
!NEC$ ivdep
         do l = 1, nl
             vi => vic(l)
           do k2 = 1, VI% g% n
#else                   /* manual loop switch for cross-grained NEC nfort */
           do k2 = 1, k2m
!NEC$ ivdep
         do l = 1, nl
#define      VI    vic(l)
             if (k2 > VI% g% n) cycle
#endif
             ix1   = VI% g% i+k2-1
             eh (l)= eh (l)                                              &
                   + VI% wh(k2) * (spt%h%w(1) * covm% eh (ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% eh (ix1,spt%h%ix+1))
             et (l)= et (l)                                              &
                   + VI% wh(k2) * (spt%h%w(1) * covm% et (ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% et (ix1,spt%h%ix+1))
             ev (l)= ev (l)                                              &
                   + VI% wh(k2) * (spt%h%w(1) * covm% ev (ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% ev (ix1,spt%h%ix+1))
             eq (l)= eq (l)                                              &
                   + VI% wh(k2) * (spt%h%w(1) * covm% erh(ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% erh(ix1,spt%h%ix+1))
             ehz(l)= ehz(l)                                              &
                   + VI% wt(k2) * (spt%h%w(1) * covm% eh (ix1,spt%h%ix  )&
                                  +spt%h%w(2) * covm% eh (ix1,spt%h%ix+1))
#ifdef VI
#undef VI
#endif
           end do
          end do
          eh (1:nl)   = emodll * eh (1:nl)
          et (1:nl)   = emodll * et (1:nl)
          ev (1:nl)   = emodll * ev (1:nl)
          ehz(1:nl)   = emodll * ehz(1:nl)
!FTRACE_END  ("init_spot:ifc_other")
        endif ! lf == ...
        xv (1:nl) = lnpc(1:nl)
        xq (1:nl) = lnpc(1:nl)
        xh (1:nl) = lnpc(1:nl)
        xhz(1:nl) = 1._wp
        et (1:nl) = et(1:nl) * (R / gacc) ! thickness error
      else
        !--------------------------
        ! old (OI) covariance model
        !--------------------------
        allocate (lnpf(nl))
        if (zflip) then
           lnpf(1:nl) = log(100000000._wp) - lnpc(1:nl)
        else
           lnpf(1:nl) =                      lnpc(1:nl)
        end if
        xh (1:nl)  = 0._wp
        xv (1:nl)  = 0._wp
        xq (1:nl)  = 0._wp
        xhz(1:nl)  = 0._wp
        !------------------------------------
        ! polynominal expansion for p > 10hPa
        !------------------------------------
!NEC$ unroll(7)
        do i = 7,1,-1
           eh  = eh  * lnpf +  spt% czh (i) ! e height
           ev  = ev  * lnpf +  spt% czv (i) ! e wind
           eq  = eq  * lnpf +  spt% czq (i) ! e relhum
           xh  = xh  * lnpc +  spt% cxh (i) ! x height
           xv  = xv  * lnpc +  spt% cxv (i) ! x wind
           xq  = xq  * lnpc +  spt% cxq (i) ! x relhum
        end do
        !-------------------------------
        ! calculate derivatives of
        ! height variance and coordinate
        !-------------------------------
!NEC$ unroll(6)
        do i = 7,2,-1
           xhz = xhz * lnpc + (spt% cxh (i) * (i-1))
           ehz = ehz * lnpf + (spt% czh (i) * (i-1))
        end do
        if (zflip) ehz = - ehz
        !-----------------------------------
        ! fix for p<10hPa:
        !   constant variances
        !   linear coordinate transformation
        !-----------------------------------
!NEC$ ivdep
        do l = 1, nl
           if (lnp(l) /= lnpc(l)) then
              xvz = 0._wp
              xqz = 0._wp
!NEC$ unroll(6)
              do i = 7,2,-1
                 xvz = xvz * lnpc(l) + (spt% cxv (i) * (i-1))
                 xqz = xqz * lnpc(l) + (spt% cxq (i) * (i-1))
              end do
              xh (l) = xh(l) + (lnp(l)-lnpc(l)) * xhz(l)
              xv (l) = xv(l) + (lnp(l)-lnpc(l)) * xvz
              xq (l) = xq(l) + (lnp(l)-lnpc(l)) * xqz
              ehz(l) = 0._wp
           endif
        end do
        !----------------------------------
        ! derive thickness error from
        ! height error and its derivatives.
        !----------------------------------
        et(1:nl) = sqrt (ehz(1:nl)**2 + &
                         xhz(1:nl)**2 * eh(1:nl)**2 * zedem1 / spt% cch)
      end if
      !------------------------------------------
      ! no humidity analysis in upper troposphere
      !------------------------------------------
      if (e_hum_top > 0._wp) then
!NEC$ ivdep
        do l = 1, nl
           if (spt% col(klev(l))% z < ln_hum_ana_top) eq(l) = e_hum_top
        end do
      end if
      !------------------------------
      ! lower bound on humidity error
      !------------------------------
      if (emin_q > 0._wp) then
        eq(1:nl) = max (eq(1:nl), emin_q)
      end if
      !----------------------------------------
      ! Save computed errors at levels to spots
      !----------------------------------------
!NEC$ ivdep
      do l = 1, nl
         k = klev(l)
#if defined (__NEC__)
#define  C    spt% col(k)
#else
         c => spt% col(k)
#endif
         C% eh   = eh (l)
         C% ev   = ev (l)
         C% eq   = eq (l)
         C% ehz  = ehz(l)
         C% xh   = xh (l)
         C% xv   = xv (l)
         C% xq   = xq (l)
         C% xhz  = xhz(l)
         !------------------------------
         ! rescale to temperature error.
         !------------------------------
         if (e_h /= 0._wp) then
           C% ezn  = ehz(l)          / et(l)
           C% exzn = eh (l) * xhz(l) / et(l)
         else
           C% ezn  = 0._wp
           C% exzn = 0._wp
         endif
         et(l)   = et (l) * (gacc / R)
         C% et   = et (l)
#ifdef C
#undef C
#endif
      end do
FTRACE_END  ("init_spot:errors_calc")
!FTRACE_BEGIN("init_spot:errors_set")
      !-----------------------
      ! loop over observations
      ! set variances
      !-----------------------
!NEC$ ivdep
      do k = 1, n
         l = llev(k)
         select case (spt% col(k)% it)
         case (OBS_H, OBS_HS)
            spt% col(k)% e = eh(l)
         case (OBS_TV, OBS_T)
            spt% col(k)% e = et(l)
         case (OBS_U, OBS_V, OBS_FF)
            spt% col(k)% e = ev(l)
         case (OBS_RH)
            spt% col(k)% e = eq(l)
         case (OBS_DUM, OBS_CHI)
            spt% col(k)% e = 1._wp
         case default
            spt% col(k)% e = -HUGE (1._wp)      ! Tag bad observation type
         end select
      end do
      if (any (spt% col(:n)% e == -HUGE (1._wp))) then
         do k = 1, n
            if (spt% col(k)% e < -1._wp) then
               write(0,*)
               write(0,*)'init_spot: observation type =', spt% col(k)% it
               l = llev(k)
               write(0,*) "Level:p,k,l:", exp (spt% col(k)% z),k,l
               write(0,*) "eh,et,ev,eq:", eh(l), et(l), ev(l), eq(l)
            end if
         end do
         call finish('init_spot','invalid observation type or error value')
      end if
!FTRACE_END  ("init_spot:errors_set")
    endif
    !----------------
    ! return variance
    !----------------
    if (present(err)) err(:) = spt% col(:)% e
  end subroutine init_spot
!------------------------------------------------------------------------------
  subroutine init_spot_obs (rc, spot, obs, err)
  type(t_rowcol)        ,intent(inout) :: rc     ! row/column data
  type(t_spot)          ,intent(in)    :: spot   ! spot data
  type(t_obs)           ,intent(in)    :: obs    ! observation data
  real(wp)    ,optional ,intent(inout) :: err(:) ! error (stddev)
  !----------------------------------------------------------------------
  ! Sets variable 'spot' according to the ' interpolation type' (spot%
  ! int_type) of the observation operator:
  ! ITY_ICOL:  The observation operator depends on exactly one horizontal
  !            location.
  ! ITY_ICOLS: The observation operator depends on several horizontal
  !            locations.
  ! ITY_MCOLS: The observation operator depends on several model columns.
  !----------------------------------------------------------------------
   integer                 :: m, n, m1, mn, mm, nc, mm_
   real(wp)       ,pointer :: lev(:)
!  real(wp)   ,allocatable :: p  (:)
   type (t_coord) ,pointer :: c
   integer                 :: nwv

   nwv = -1; if(spot% hd% obstype == OT_RAD) nwv = min (nwv_rad, 3)
   select case (spot% int_type)
   case (ITY_ICOL)
     n = rc% inh + 1
     call init_spot (rc,                &!
                     n,                 &! spot index within the box
                     spot% col% c,      &! coordinates
                     obs% t_int(spot%i%i+1:spot%i%i+spot%i%n) ,&!
                     obs% lev  (spot%i%i+1:spot%i%i+spot%i%n) ,&!
                 nwv=nwv,                                      &!
                 err=err)
   case (ITY_ICOLS)
     call finish &
       ('init_spot_obs','interpolation type ITY_ICOLS not implemented')
   case (ITY_MCOLS)
     if (spot% i% n /= 0) then        !  skip for verification, no 3dvar
       nc = size (spot% imcol)
       do m =  1, nc    ! loop over model columns
         c   => spot% imcol(m)% c
         mm_ =  spot% mke * spot% imcol(m)% natm + spot% imcol(m)% nsl
         mm  =  spot% i% n / nc
         m1  =  (m-1) * mm + 1           ! indices to param. ids
         mn  =  m1 +    mm - 1           !   and levels
         n = rc% inh + 1

         if(nc * mm /= spot% i% n) then
           write(0,*)  'init_spot_obs:  nc * mm /= spot% i% n'
           write(0,*)  'nc        ',nc
           write(0,*)  'mke       ',spot% mke
           write(0,*)  'natm      ',spot% imcol(m)% natm
           write(0,*)  'nsl       ',spot% imcol(m)% nsl
           write(0,*)  'mm_       ',mm_
           write(0,*)  'mm        ',mm
           write(0,*)  'spot% i% n',spot% i% n
           call finish('init_spot_obs','nc * mm /= spot% i% n')
         endif

         lev => obs% lev (spot%i%i + m1 : spot%i%i + mn)

         if (present(err)) then
           call init_spot (rc,       &!
                           n,        &! spot index within the box
                           c,        &! coordinates
                           obs% t_int(spot%i%i + m1 : spot%i%i + mn) ,&!
                           lev,      &
                       err=err(m1:mn))
         else
           call init_spot (rc,       &!
                           n,        &! spot index within the box
                           c,        &! coordinates
                           obs% t_int(spot%i%i + m1 : spot%i%i + mn) ,&!
                           lev)
         endif
       end do
     endif
   case default
     call finish ('init_spot_obs','invalid interpolation type')
   end select

  end subroutine init_spot_obs
!------------------------------------------------------------------------------
  subroutine get_corr (lhs, rhs, is, js, lcov, lsym, nonzero, c, r, l)
  type(t_rowcol)    ,intent(in)    :: lhs     ! lhs description
  type(t_rowcol)    ,intent(in)    :: rhs     ! rhs description
  integer           ,intent(in)    :: is      ! row spot index
  integer           ,intent(in)    :: js      ! column spot index
  logical           ,intent(in)    :: lcov    ! .T.:covariance, .F.:correlation
  logical           ,intent(in)    :: lsym    ! .T. for symmetric matrix
  integer ,optional ,intent(inout) :: nonzero ! number of nonzero elements
  real(mp),optional ,intent(inout) :: c (:,:) ! covariance matrix
  real(wp),optional ,intent(inout) :: l (:)   ! left hand side
  real(wp),optional ,intent(in)    :: r (:)   ! right hand side
  target                           :: c
  !========================================================================
  ! get correlation (or covariance for lcov=.t.) matrix for a pair of spots
  !========================================================================
    !----------------------------------
    ! horizontal correlation parameters
    !----------------------------------
    real(wp) :: chq     ! horizontal correlation for relative humidity
    real(wp) :: chh     ! horizontal correlation for height
    real(wp) :: clsvv   ! large scale contribution to wind-wind correlation
    real(wp) :: rlsvv   ! renormalisation due to large scale contribution
    real(wp) :: d1chh   ! first  derivative of chh
    real(wp) :: d1chx   ! first  derivative of chh, truncated for h-w corr.
    real(wp) :: d2chh   ! second derivative of chh
    real(wp) :: d1chho1 ! first  derivative of chh / 1.
    !----------------------------------------------------------------
    ! scalar product of longitudinal and transversal (l,t) directions
    ! with local (u,v) coordinate system
    !----------------------------------------------------------------
    real(wp) :: ab(3) ! vector between locations on the unit sphere
    real(wp) :: d     ! distance between locations on the unit sphere
    real(wp) :: abn   ! normalization factor
    real(wp) :: elui  ! longitudinal * u ,right hand side
    real(wp) :: elvi  ! longitudinal * v ,right hand side
    real(wp) :: etui  ! transversal  * u ,right hand side
    real(wp) :: etvi  ! transversal  * v ,right hand side
    real(wp) :: eluj  ! longitudinal * u ,left  hand side
    real(wp) :: elvj  ! longitudinal * v ,left  hand side
    real(wp) :: etuj  ! transversal  * u ,left  hand side
    real(wp) :: etvj  ! transversal  * v ,left  hand side
    !--------------------------------
    ! vertical correlation parameters
    !--------------------------------
    real(wp) :: dxh   ! difference of transformed heights for height
    real(wp) :: dxv   ! difference of transformed heights for wind
    real(wp) :: dxq   ! difference of transformed heights for humidity

    real(wp) :: dxh2  ! difference of transformed heights squared
    real(wp) :: dxv2  ! difference of transformed heights squared
    real(wp) :: dxq2  ! difference of transformed heights squared

    real(wp) :: edxh2 ! 1 / (c_h + edxh2)
    real(wp) :: eexcv ! edem1 * exp ( c_h * edxh2)

    real(wp) :: cvq   ! vertical correlation (for humidity)
    real(wp) :: cvx   ! vertical correlation (for velocity potential)
    real(wp) :: cv    ! vertical correlation (for height)
    real(wp) :: cvzl  ! first  derivative of cv (lhs)
    real(wp) :: cvzr  ! first  derivative of cv (rhs)
    real(wp) :: cvzz  ! second derivative of cv
    real(wp) :: cti   ! vertical correlation for temperature at left  hand side
    real(wp) :: ctj   ! vertical correlation for temperature at right hand side
    real(wp) :: cvtt  ! vertical correlation for temperature

    real(wp) :: c_t_vt  ! correlation temperature - transv.wind
    real(wp) :: c_h_vt  ! correlation height      - transv.wind
    real(wp) :: c_vt_t  ! correlation transv.wind - temperature
    real(wp) :: c_vt_h  ! correlation transv.wind - height
    real(wp) :: c_vt_vt ! correlation transv.wind - transv.wind
    real(wp) :: c_vl_vl ! correlation longit.wind - longit.wind

    real(wp) :: c_h     ! correlation coefficient for height
    real(wp) :: c_q     ! correlation coefficient for humidity
    real(wp) :: c_v     ! correlation coefficient for wind
    !--------------------------------------
    ! NMC fitted vertical correlation model
    !--------------------------------------
    integer               :: ill        ! level index lhs
    integer               :: ilr        ! level index rhs
    !-----------------
    ! other parameters
    !-----------------
    real(wp)                  :: r1nu      ! 1 - rnu
    real(wp)                  :: sq1nu     ! sqrt(1-rnu)
    real(wp)                  :: lh, lq
    integer                   :: il,jl
    integer                   :: i,j
    integer                   :: nonz
    integer                   :: itl_a_itl ! itl .and. itl
    integer                   :: itl_o_itl ! itl .or.  itl
    type(t_dwd_spot) ,pointer :: si, sj    ! spot data type
    type(t_dwd_col)  ,pointer :: xi, xj    ! level data type
#if !defined (GET_CORR_NEW_INNER)
    integer                   :: it, jt    ! observation types
    real(wp)                  :: cij       ! Correlation
#endif
    real(mp), pointer         :: cij_(:,:) ! Temporary covariance matrix
    real(wp) ,parameter       :: edem1 = 1 / (e-1)

    !-------------------
    ! set some constants
    !-------------------
    r1nu  = 1._wp - rnu
    sq1nu = sqrt(r1nu)
    !-----------------
    ! set spot pointer
    !-----------------
    si => lhs% spots (is)
    sj => rhs% spots (js)
    !----------------------------
    ! return for zero matrix size
    !----------------------------
    if (si% n == 0 .or. sj% n == 0) return
    !-----------------------------------------
    ! allocate and preset correlation matrix c
    !-----------------------------------------
    if (present (c)) then
       cij_ => c
    else
       allocate (cij_(si% n, sj% n))
    end if
    cij_ = 0._mp
    !-------------------------------
    ! shortcut for single level data
    !-------------------------------
    if (shortc .and. lsym .and. is==js .and. si% col(1)% i == si%n) then
      forall (i=1:min (si% n, sj% n)) cij_(i,i) = 1._mp
      nonz = sj% n
      if (lcov) then
!NEC$ ivdep
        do j = 1, sj% n
          !-----------------------------------------
          ! if requested
          ! return variances instead of correlations
          !-----------------------------------------
          cij_(j,j) = si%col(j)% e * sj%col(j)% e
        end do
      end if
      if (present (r) .and. present (l)) then
        ! (almost) DGEMV: l = l + real(C,wp)*r
!NEC$ ivdep
        do j = 1, sj% n
          l(j) = l(j) + cij_(j,j) * r(j)
        end do
      end if
    else
      !--------------------------------------
      ! distance on the surface of the planet
      !--------------------------------------
      ab   = sj% x - si% x
      d    = sqrt (sum ((si% xt - sj% xt)**2) )
      !------------------------------------------------------------
      ! set horizontal correlation functions for height
      ! and derivatives for wind correlations
      ! return from routine for large distances (zero correlations)
      !------------------------------------------------------------
      if (thc) then
        lh   = lh_max
      else
        lh   = min (si% L_h , sj% L_h)
      endif
      if (d >= x_cutdr * lh) return
      chh    = corr (d, lh)
      clsvv  = lsvv * chh
      rlsvv  = 1._wp / (1._wp + lsvv)
      call dcorr (d1chh, d1chho1, d2chh, d1chx)
      if (chh == 0._wp .and. d1chh == 0._wp) return
      !--------------------------------------------------
      ! set horizontal correlation functions for moisture
      ! default: exponential correlations (rhcexp=.true.)
      !--------------------------------------------------
      lq  = min (si% L_q , sj% L_q)
      chq = 0._wp
      if (d < x_cutdr * lq) then
        chq = corrh (d,lq)
      endif
      !----------------------------------------------
      ! set latitude dependent vertical length scales
      !----------------------------------------------
      c_h  = max (si% cch , sj% cch)
      c_q  = max (si% ccq , sj% ccq)
      c_v  = max (si% ccv , sj% ccv)
      !------------------------------------------------------
      ! transformation from local (East/North) coordinates to
      ! longitudinal/transversal coordinates
      !------------------------------------------------------
      if (d==0._wp) ab = si% du                 ! arbitrary direction for d=0
      elui =   sum (ab * si% du)                ! longit.vector in u,v system
      elvi =   sum (ab * si% dv)                !
      abn  =   1._wp / sqrt (elui**2 + elvi**2) ! normalize
      elui =   abn * elui                       !
      elvi =   abn * elvi                       !
      eluj =   sum (ab * sj% du)                ! same at left hand side
      elvj =   sum (ab * sj% dv)                !  (changed sign)
      abn  =   1._wp / sqrt (eluj**2 + elvj**2) !
      eluj = - (abn * eluj)                     !
      elvj = - (abn * elvj)                     !
      etui = - elvi                             ! rotate for tang.vector
      etvi =   elui                             !
      etuj = - elvj                             !
      etvj =   eluj

      !--------------------------
      ! loop over pairs of levels
      !--------------------------
!     call f_hpmstart(7,'get_corr: outer loop')
      il   = 1
      ill  = 0
      nonz = 0
      do
        ill =  ill + 1
        xi  => si% col (il)
        jl  =  1
        ilr =  0
        do
          ilr = ilr + 1
          xj => sj% col (jl)
          itl_a_itl = iand (xi% itl, xj % itl)
          itl_o_itl = ior  (xi% itl, xj % itl)
          !=============================
          ! NMC-fitted covariance model:
          !=============================
          if (covm% valid > 2) then
            !----------------------------------
            ! set vertical correlation function
            !----------------------------------
            !------------------
            ! height + X, T + T
            !------------------
            if (iand ( itl_o_itl, OBS_TV ) /= 0) then
              !---------------------------------------------------
              ! Optimized for SX-6 (but fast on IBM, too!)
              ! Do not touch unless you know what you are doing...
              !---------------------------------------------------
              cv   = 0
              cvzl = 0
              cvzr = 0
              cvzz = 0
!NEC$ shortloop
              do i = 1, covm% nz
                cv   = cv   + si%vm (i, INT_H,  ill) * sj%vm (i, INT_H,  ilr)
                cvzl = cvzl + si%vm (i, INT_TV, ill) * sj%vm (i, INT_H,  ilr)
                cvzr = cvzr + si%vm (i, INT_H , ill) * sj%vm (i, INT_TV, ilr)
                cvzz = cvzz + si%vm (i, INT_TV, ill) * sj%vm (i, INT_TV, ilr)
              end do
              cvzr = -cvzr
              !--------------------------------------------------------------
              ! Original (slower) version of the above code:
              ! cv   =  sum (si%vm (:, INT_H,  ill) * sj%vm (:, INT_H,  ilr))
              ! cvzl =  sum (si%vm (:, INT_TV, ill) * sj%vm (:, INT_H,  ilr))
              ! cvzr = -sum (si%vm (:, INT_H , ill) * sj%vm (:, INT_TV, ilr))
              ! cvzz =  sum (si%vm (:, INT_TV, ill) * sj%vm (:, INT_TV, ilr))
              !--------------------------------------------------------------
              cti  =  ( xi% ezn * cv + xi% exzn * cvzl )
              ctj  =  ( xj% ezn * cv - xj% exzn * cvzr )
              cvtt =  ( xi%ezn  * cv   * xj%ezn  &
                      + xi%exzn * cvzl * xj%ezn  &
                      - xi%ezn  * cvzr * xj%exzn &
                      + xi%exzn * cvzz * xj%exzn )
              c_t_vt  = sj% mu * sq1nu * d1chx * cti * rlsvv
              c_vt_t  = si% mu * sq1nu * d1chx * ctj * rlsvv
            else
            !-----------
            ! height + X
            !-----------
              cv  =  sum (si%vm (:, INT_H,  ill) * sj%vm (:, INT_H,  ilr))
            endif
            !--------
            ! RH - RH
            !--------
            if (iand ( itl_a_itl, OBS_RH) /= 0) then
              if (chq /= 0._wp) then
                cvq = sum (si%vm (:, INT_RH, ill) * sj%vm (:, INT_RH, ilr))
              else
                cvq = 0._wp
              endif
            endif
            !------------
            ! wind - wind
            !------------
            if (iand ( itl_a_itl, ior(OBS_U,OBS_V)) /= 0) then
              cvx = sum (si%vm (:, INT_CHI, ill) * sj%vm (:, INT_CHI, ilr))
              c_vt_vt = ( cv  * (r1nu * d2chh   - clsvv) &
                        + cvx *  rnu  * d1chho1        ) * rlsvv
              c_vl_vl = ( cv  * (r1nu * d1chho1 - clsvv) &
                        + cvx *  rnu  * d2chh          ) * rlsvv
            endif
            !---------
            ! X + wind
            !---------
            if (iand ( itl_o_itl, ior(OBS_U,OBS_V)) /= 0) then
              c_h_vt  = sj% mu * sq1nu * d1chx * cv * rlsvv
              c_vt_h  = si% mu * sq1nu * d1chx * cv * rlsvv
            endif
          else
          !=================================================
          ! set vertical correlation function for height + X
          !=================================================
            dxh  = xi% xh - xj% xh ! difference of vert. coord. for height
            dxh2 = dxh * dxh        ! squared differences
            edxh2 = 1._wp / (c_h + dxh2)
            eexcv = edem1 * exp ( c_h * edxh2)
            cv   =  eexcv - edem1
            !========
            ! RH - RH
            !========
            if (iand ( itl_a_itl, OBS_RH) /= 0) then
              if (chq /= 0._wp) then
                dxq  = xi% xq - xj% xq ! difference of vert.coord. for humidity
                dxq2 = dxq * dxq
                cvq  =  edem1 * exp ( c_q / (c_q + dxq2) ) - edem1
              else
                cvq  = 0._wp
              endif
            endif
            !============
            ! wind - wind
            !============
            if (iand ( itl_a_itl, ior(OBS_U,ior(OBS_V,OBS_CHI))) /= 0) then
              dxv  = xi% xv - xj% xv ! difference of vert. coord. for wind
              dxv2 = dxv * dxv
              cvx  = (edem1 * exp ( c_v / (c_v+(dxv2))) - edem1)
              if (modcvx>0._wp) &
                cvx = cvx * (1._wp-0.75_wp*dxv2) / (1._wp+modcvx*dxv2)
              c_vt_vt = ( cv  * (r1nu * d2chh   - clsvv) &
                        + cvx *  rnu  * d1chho1        ) * rlsvv
              c_vl_vl = ( cv  * (r1nu * d1chho1 - clsvv) &
                        + cvx *  rnu  * d2chh          ) * rlsvv
            endif
            !======
            ! T + T
            !======
            if (iand ( itl_o_itl, OBS_TV ) /= 0) then
              cvzl = -2._wp * ((c_h*dxh) * edxh2**2) &
                      * eexcv
              cvzr = cvzl
              cvzz = ( 2._wp - dxh2*edxh2 * ( 4._wp* c_h * edxh2 + 8._wp )) &
                     * eexcv * c_h * edxh2 * edxh2
              cti  =  ( xi% ezn * cv + xi% exzn * cvzl )
              ctj  =  ( xj% ezn * cv - xj% exzn * cvzr )
              cvtt =  ( xi%ezn  * cv   * xj%ezn  &
                      + xi%exzn * cvzl * xj%ezn  &
                      - xi%ezn  * cvzr * xj%exzn &
                      + xi%exzn * cvzz * xj%exzn )
              c_t_vt  = sj% mu * sq1nu * d1chx * cti * rlsvv
              c_vt_t  = si% mu * sq1nu * d1chx * ctj * rlsvv
            endif
            !=========
            ! X + wind
            !=========
            if (iand ( itl_o_itl, ior(OBS_U,OBS_V)) /= 0) then
              c_h_vt  = sj% mu * sq1nu * d1chx * cv * rlsvv
              c_vt_h  = si% mu * sq1nu * d1chx * cv * rlsvv
            endif
          endif
          !----------------------------------------------
          ! loop over pairs of observations at this level
          ! fill correlation matrix
          !----------------------------------------------
!         call f_hpmstart(8,'get_corr: inner loop')
#if defined (GET_CORR_NEW_INNER)
          !-----------------------------------------------------------
          ! Source code optimized for NEC SX.  I tried to make it more
          ! readable using preprocessor directives for innermost parts
          !-----------------------------------------------------------
#define BEGIN_INNER_LOOP \
          do i = il, il+xi%i-1;\
            select case (si%col(i)% it)
#define CASE(x,y) \
            case (x);\
              nonz = nonz+1;\
              cij_(i,j) = y
#define CASES(x1,x2,y) \
            case (x1,x2);\
              nonz = nonz+1;\
              cij_(i,j) = y
#define END_INNER_LOOP \
            end select;\
          end do
          !------------------------------------------------------
          ! All loops are too short to benefit from vectorization
          !------------------------------------------------------
#define VECTOR_OPTS     NEC$ NOVECTOR
          !------------------------------------------------------
          do j = jl, jl+xj%i-1
            select case (sj%col(j)% it)
            case(OBS_RH)
!VECTOR_OPTS
              BEGIN_INNER_LOOP
              CASE (OBS_RH,            cvq * chq)
              END_INNER_LOOP
            case(OBS_H, OBS_HS)
!VECTOR_OPTS
              BEGIN_INNER_LOOP
              CASES(OBS_H,OBS_HS,      cv  * chh)
              CASE (OBS_TV,          - cti * chh)
              CASE (OBS_U,        - c_vt_h * etui)
              CASE (OBS_V,        - c_vt_h * etvi)
              END_INNER_LOOP
            case(OBS_TV)
!VECTOR_OPTS
              BEGIN_INNER_LOOP
              CASES(OBS_H,OBS_HS,    - ctj * chh)
              CASE (OBS_TV,           cvtt * chh)
              CASE (OBS_U,          c_vt_t * etui)
              CASE (OBS_V,          c_vt_t * etvi)
              END_INNER_LOOP
            case(OBS_U)
!VECTOR_OPTS
              BEGIN_INNER_LOOP
              CASES(OBS_H,OBS_HS, - c_h_vt * etuj)
              CASE (OBS_TV,         c_t_vt * etuj)
              CASE (OBS_U,  etui * c_vt_vt * etuj + elui * c_vl_vl * eluj)
              CASE (OBS_V,  etvi * c_vt_vt * etuj + elvi * c_vl_vl * eluj)
              END_INNER_LOOP
            case(OBS_V)
!VECTOR_OPTS
              BEGIN_INNER_LOOP
              CASES(OBS_H,OBS_HS, - c_h_vt * etvj)
              CASE (OBS_TV,         c_t_vt * etvj)
              CASE (OBS_U,  etui * c_vt_vt * etvj + elui * c_vl_vl * elvj)
              CASE (OBS_V,  etvi * c_vt_vt * etvj + elvi * c_vl_vl * elvj)
              END_INNER_LOOP
            case(OBS_CHI)
!VECTOR_OPTS
              BEGIN_INNER_LOOP
              CASE (OBS_CHI,       cvx)
              END_INNER_LOOP
            end select
          end do
          !------------------------------------------------------
#undef BEGIN_INNER_LOOP
#undef CASE
#undef CASES
#undef END_INNER_LOOP
#undef VECTOR_OPTS
#else  /* !GET_CORR_NEW_INNER */
!NEC$ nomove
          do j = jl, jl+xj%i-1
            jt = sj%col(j)% it
! The following loop is too short to benefit from vectorization
!NEC$ nomove
!NEC$ novector
            do i = il, il+xi%i-1
              it = si%col(i)% it
!             call f_hpmstart(10,'get_corr: inner body')
              select case (it)
              case(OBS_RH)
                select case (jt)
                case(OBS_RH)
                  cij = cvq * chq
                case default
                  goto 99
                end select
              case(OBS_H, OBS_HS)
                select case (jt)
                case(OBS_H, OBS_HS)
                  cij = cv  * chh
                case(OBS_TV)
                  cij = - ctj * chh
                case(OBS_U)
                  cij = - c_h_vt * etuj
                case(OBS_V)
                  cij = - c_h_vt * etvj
                case default
                  goto 99
                end select
              case(OBS_TV)
                select case (jt)
                case(OBS_H, OBS_HS)
                  cij = - cti * chh
                case(OBS_TV)
                  cij = cvtt * chh
                case(OBS_U)
                  cij = c_t_vt * etuj
                case(OBS_V)
                  cij = c_t_vt * etvj
                case default
                  goto 99
                end select
              case(OBS_U)
                select case (jt)
                case(OBS_H, OBS_HS)
                  cij = - c_vt_h * etui
                case(OBS_TV)
                  cij = c_vt_t * etui
                case(OBS_U)
                  cij = etui * c_vt_vt * etuj &
                      + elui * c_vl_vl * eluj
                case(OBS_V)
                  cij = etui * c_vt_vt * etvj &
                      + elui * c_vl_vl * elvj
                case default
                  goto 99
                end select
              case(OBS_V)
                select case (jt)
                case(OBS_H, OBS_HS)
                  cij = - c_vt_h * etvi
                case(OBS_TV)
                  cij = c_vt_t * etvi
                case(OBS_U)
                  cij = etvi * c_vt_vt * etuj &
                      + elvi * c_vl_vl * eluj
                case(OBS_V)
                  cij = etvi * c_vt_vt * etvj &
                      + elvi * c_vl_vl * elvj
                case default
                  goto 99
                end select
              case(OBS_CHI)
                select case (jt)
                case(OBS_CHI)
                  cij = cvx
                case default
                  goto 99
                end select
              case default
                goto 99
              end select
              !------------
              ! diagnostics
              !------------
              nonz = nonz + 1
              cij_(i,j) = cij
99            continue
!             call f_hpmstop(10)
            end do
          end do
#endif /* GET_CORR_NEW_INNER */
!         call f_hpmstop(8)
          jl=jl+xj%i
          if (jl>sj%n) exit
        end do
        il=il+xi%i
        if (il>si%n) exit
      end do
      if (lsym .and. is==js) then
!NEC$ ivdep
         do i = 1, si% n
           !------------------------------
           ! force diagonal elements to 1.
           !------------------------------
           if (cij_(i,i) == 0._mp) then
             nonz      = nonz + 1
             cij_(i,i) = 1._mp
           endif
         end do
      end if
      if (lcov) then
         do j = 1, sj% n
!NEC$ ivdep
            do i = 1, si% n
               !-----------------------------------------
               ! if requested
               ! return variances instead of correlations
               !-----------------------------------------
               cij_(i,j) = cij_(i,j) * si%col(i)% e * sj%col(j)% e
            end do
         end do
      end if
      if (present (r) .and. present (l)) then
         ! (almost) DGEMV: l = l + real(C,wp)*r
         do j = 1, sj% n
!NEC$ ivdep
            do i = 1, si% n
               l(i) = l(i) + cij_(i,j) * r(j)
            end do
         end do
      end if
    endif
    if (present (c)) then
       nullify (cij_)
    else
       deallocate (cij_)
    end if
!   call f_hpmstop(7)
    if (present(nonzero)) nonzero = nonzero + nonz

!if (dace% lpio) print *, "get_corr: p_pe, nonz =", dace% pe, nonz

  end subroutine get_corr
!------------------------------------------------------------------------------
  subroutine get_corr_opt (lhs, rhs, lsym, i0, j0, nonzero, c)
  type(t_rowcol) ,intent(in)    :: lhs            ! lhs description
  type(t_rowcol) ,intent(in)    :: rhs            ! rhs description
  logical        ,intent(in)    :: lsym           ! .T. for symmetric matrix
  integer        ,intent(in)    :: i0             ! spots i index
  integer        ,intent(in)    :: j0             ! spots j index
  integer        ,intent(out)   :: nonzero(:,:,:) ! number of nonzero elements
  real(mp)       ,intent(inout) :: c(:,:,:,:,:)   ! covariance matrix
  !========================================================================
  ! get covariance matrix for a pair of spots
  !========================================================================
    !----------------------------------
    ! horizontal correlation parameters
    !----------------------------------
    real(wp) :: chq     ! horizontal correlation for relative humidity
    real(wp) :: chh     ! horizontal correlation for height
    real(wp) :: clsvv   ! large scale contribution to wind-wind correlation
    real(wp) :: rlsvv   ! renormalisation due to large scale contribution
    real(wp) :: d1chh   ! first  derivative of chh
    real(wp) :: d1chx   ! first  derivative of chh, truncated for h-w corr.
    real(wp) :: d2chh   ! second derivative of chh
    real(wp) :: d1chho1 ! first  derivative of chh / 1.
    !----------------------------------------------------------------
    ! scalar product of longitudinal and transversal (l,t) directions
    ! with local (u,v) coordinate system
    !----------------------------------------------------------------
    real(wp) :: ab(3) ! vector between locations on the unit sphere
    real(wp) :: d     ! distance between locations on the unit sphere
    real(wp) :: abn   ! normalization factor
    real(wp) :: elui  ! longitudinal * u ,right hand side
    real(wp) :: elvi  ! longitudinal * v ,right hand side
    real(wp) :: etui  ! transversal  * u ,right hand side
    real(wp) :: etvi  ! transversal  * v ,right hand side
    real(wp) :: eluj  ! longitudinal * u ,left  hand side
    real(wp) :: elvj  ! longitudinal * v ,left  hand side
    real(wp) :: etuj  ! transversal  * u ,left  hand side
    real(wp) :: etvj  ! transversal  * v ,left  hand side
    !--------------------------------
    ! vertical correlation parameters
    !--------------------------------
    real(wp) :: dxh   ! difference of transformed heights for height
    real(wp) :: dxv   ! difference of transformed heights for wind
    real(wp) :: dxq   ! difference of transformed heights for humidity

    real(wp) :: dxh2  ! difference of transformed heights squared
    real(wp) :: dxv2  ! difference of transformed heights squared
    real(wp) :: dxq2  ! difference of transformed heights squared

    real(wp) :: edxh2 ! 1 / (c_h + edxh2)
    real(wp) :: eexcv ! edem1 * exp ( c_h * edxh2)

    real(wp) :: cvq (nv_gc)  ! vertical correlation (for humidity)
    real(wp) :: cvx (nv_gc)  ! vertical correlation (for velocity potential)
    real(wp) :: cv  (nv_gc)  ! vertical correlation (for height)
    real(wp) :: cvzl(nv_gc)  ! first  derivative of cv (lhs)
    real(wp) :: cvzr(nv_gc)  ! first  derivative of cv (rhs)
    real(wp) :: cvzz(nv_gc)  ! second derivative of cv
    real(wp) :: cti    (nv_gc) ! vertical correlation for temperature at left  hand side
    real(wp) :: ctj    (nv_gc) ! vertical correlation for temperature at right hand side
    real(wp) :: cvtt   (nv_gc) ! vertical correlation for temperature

    real(wp) :: c_t_vt (nv_gc) ! correlation temperature - transv.wind
    real(wp) :: c_h_vt (nv_gc) ! correlation height      - transv.wind
    real(wp) :: c_vt_t (nv_gc) ! correlation transv.wind - temperature
    real(wp) :: c_vt_h (nv_gc) ! correlation transv.wind - height
    real(wp) :: c_vt_vt(nv_gc) ! correlation transv.wind - transv.wind
    real(wp) :: c_vl_vl(nv_gc) ! correlation longit.wind - longit.wind

    real(wp) :: c_h     ! correlation coefficient for height
    real(wp) :: c_q     ! correlation coefficient for humidity
    real(wp) :: c_v     ! correlation coefficient for wind
    !--------------------------------------
    ! NMC fitted vertical correlation model
    !--------------------------------------
    integer               :: ill        ! level index lhs
    integer               :: ilr        ! level index rhs
    !-----------------
    ! other parameters
    !-----------------
    real(wp)                  :: r1nu      ! 1 - rnu
    real(wp)                  :: sq1nu     ! sqrt(1-rnu)
    real(wp)                  :: lh, lq
    integer                   :: il,jl
    integer                   :: i,j
!   type(t_dwd_spot) ,pointer :: si, sj    ! spot data type
!   type(t_dwd_col)  ,pointer :: xi, xj    ! level data type
    integer                   :: it, jt    ! observation types
    integer                   :: itl_a_itl ! itl .and. itl
    integer                   :: itl_o_itl ! itl .or.  itl
    real(wp) ,parameter       :: edem1 = 1 / (e-1)
    !NECJB
    integer                   :: is, js         ! spot indices
    integer                   :: SI_n, SJ_n     ! degrees of freedom
    integer :: &
      iv, iiv, nv, nv_max, nll, nlr, nv_pol, ivp, ivp0, ivlen, &
      nv_o_tv, nv_o_tv_not, nv_a_rh, nv_a_uv, nv_o_uv, nil, njl, ii, jj, &
      nv_rh_rh, nv_hs_hs, nv_hs_tv, nv_hs_u,  nv_hs_v,  nv_tv_hs, &
      nv_tv_tv, nv_tv_u , nv_tv_v , nv_u_hs,  nv_u_tv,  nv_u_u  , &
      nv_u_v  , nv_v_hs , nv_v_tv , nv_v_u ,  nv_v_v ,  nv_ch_ch, &
      nv_pol2, jv, m, ixi, ixj, jvp
    real(wp) :: r
    real(wp), dimension(3,nv_gc) :: &
      rv_ab
    real(wp), dimension(nv_gc) :: &
      rv_d, rv_chh, rv_chq, rv_d1chh, rv_d1chho1, rv_d2chh, rv_d1chx, &
      rv_elui, rv_eluj, rv_elvi, rv_elvj, rv_etui, rv_etuj, rv_etvi, rv_etvj, &
      rv_c_h, rv_c_q, rv_c_v, rv_clsvv, rv_rlsvv
    integer,  dimension(nv_gc) :: &
      iv_nll, iv_nlr, iv_pol,     &
      iv_o_tv, iv_o_tv_not, iv_a_rh, iv_a_uv, iv_o_uv
    integer,  allocatable, dimension(:)   :: &  ! can be larger than nv_gc [ha]!
      iv_rh_rh, iv_hs_hs, iv_hs_tv, iv_hs_u,  iv_hs_v,  iv_tv_hs, &
      iv_tv_tv, iv_tv_u , iv_tv_v , iv_u_hs,  iv_u_tv,  iv_u_u  , &
      iv_u_v  , iv_v_hs , iv_v_tv , iv_v_u ,  iv_v_v ,  iv_ch_ch
    integer,  allocatable, dimension(:,:) :: iv_il, iv_jl
    integer,  allocatable, dimension(:,:) :: iv_pol1, iv_pol2

    integer, parameter :: MAXVL = 256        ! Maximum vector length
    real(wp), dimension(257,covm%nz,4)    :: rv_sum
!   real(wp), dimension(MAXVL+1,covm%nz,4):: rv_sum

!    real(wp), dimension(256) :: rv_cv, rv_cvzl, rv_cvzr, rv_cvzz
!    integer,  dimension(256) :: iv_is, iv_js
!!cdir vreg(rv_cv, rv_cvzl, rv_cvzr, rv_cvzz, iv_is, iv_js)

!------------------------
! set spot pointer macros
!------------------------
#define  SI  lhs%spots(is)
#define  SJ  rhs%spots(js)
#define  XI  SI%col(il)
#define  XJ  SJ%col(jl)

    c  (:,:,:,i0,j0) = 0._mp
    nonzero(:,i0,j0) = 0

    !-------------------
    ! set some constants
    !-------------------
    r1nu  = 1._wp - rnu
    sq1nu = sqrt(r1nu)

FTRACE_BEGIN("get_corr_opt:00")
    !call ftrace_region_begin( 'get_corr_opt_1' )

    if (lsym) then
      !----------------------------
      ! preset correlation matrix c
      !----------------------------
      do i = 1, min( SIZE(c,1), SIZE(c,2) )
!NEC$ ivdep
        do iiv = 1, nv_gc

          is = tv_gc(iiv)%i
          js = tv_gc(iiv)%j
          iv = tv_gc(iiv)%iv

          if (is==js .and. i <= min (SI% n, SJ% n)) then
             c(i,i,iv,i0,j0) = 1._mp
          end if
        end do
      end do
    end if

    !call ftrace_region_end  ( 'get_corr_opt_1' )
    !call ftrace_region_begin( 'get_corr_opt_2' )

    nv = 0
    do iiv = 1, nv_gc

      is = tv_gc(iiv)%i
      js = tv_gc(iiv)%j
      iv = tv_gc(iiv)%iv

      !-------------------------------
      ! shortcut for single level data
      !-------------------------------
      if (shortc .and. lsym .and. is==js .and. SI% col(1)% i == SI%n) then
        nonzero(iv,i0,j0) = SJ% n
!NEC$ ivdep
        do j = 1, SJ% n
          !-----------------------------------------
          ! return variances instead of correlations
          !-----------------------------------------
          c(j,j,iv,i0,j0) = SI%col(j)% e * SJ%col(j)% e
        end do
      else
        nv = nv + 1
        tv_gc(nv)%i  = is
        tv_gc(nv)%j  = js
        tv_gc(nv)%iv = iv
      end if
    end do

    !call ftrace_region_end  ( 'get_corr_opt_2' )
FTRACE_BEGIN("get_corr_opt:3")

    nv_max = nv
    nv = 0

    if( model == 'csxoar2' )then

      !NECJB Routines corr and dcorr inlined manually for model 'csxoar2'

!NEC$ ivdep
      do iiv = 1, nv_max

        is = tv_gc(iiv)%i
        js = tv_gc(iiv)%j
        iv = tv_gc(iiv)%iv

        !--------------------------------------
        ! distance on the surface of the planet
        !--------------------------------------
        ab(1) = SJ% x(1) - SI% x(1)
        ab(2) = SJ% x(2) - SI% x(2)
        ab(3) = SJ% x(3) - SI% x(3)
        d     = sqrt( (SI% xt(1) - SJ% xt(1))**2 + &
                      (SI% xt(2) - SJ% xt(2))**2 + &
                      (SI% xt(3) - SJ% xt(3))**2     )
        !------------------------------------------------------------
        ! set horizontal correlation functions for height
        ! and derivatives for wind correlations
        ! return from routine for large distances (zero correlations)
        !------------------------------------------------------------
        if (thc) then
          lh = lh_max
        else
          lh = min (SI% L_h , SJ% L_h)
        endif
        if (d < x_cutdr * lh) then

          !NECJB chh = corr (d, lh)
          if (sec2rad) then
            r = 2._wp * asin (0.5_wp * d)
          else
            r = d
          endif
          r = r * L1r / lh
          if (r < x_cuts) then
            chh = (gna + gnb * r + gnc * r*r) * (1-r)**(a+2)
          else
            chh = 0._wp
          endif

          !NECJB call dcorr (d1chh, d1chho1, d2chh, d1chx)
          if (r >= x_cuts) then
            d1chh   = 0._wp
            d1chho1 = 0._wp
            d2chh   = 0._wp
            d1chx   = 0._wp
          else
            d1chh = (gne + gnf * r + gng * r*r) * (1-r)**(a+1)
            d2chh = (gni + gnj * r + gnk * r*r) * (1-r)**(a)
            d1chx = d1chh
            !----------------------------------
            ! for r=0: d/dr corr is not defined
            !----------------------------------
            if (r /= 0._wp) then
              d1chho1 = d1chh / r
            else
              d1chho1 = -1._wp / L2
            endif
            !--------------------------------------
            ! rescale derivatives with length scale
            !--------------------------------------
            d1chh   = d1chh   * L1
            d1chho1 = d1chho1 * L2
            d2chh   = d2chh   * L2
            d1chx   = d1chx   * L1
          endif
          if (chh /= 0._wp .or. d1chh /= 0._wp) then
            nv = nv + 1

            tv_gc(nv)%i  = is
            tv_gc(nv)%j  = js
            tv_gc(nv)%iv = iv

            rv_ab   (1,nv) = ab(1)
            rv_ab   (2,nv) = ab(2)
            rv_ab   (3,nv) = ab(3)
            rv_d      (nv) = d
            rv_chh    (nv) = chh
            rv_d1chh  (nv) = d1chh
            rv_d1chho1(nv) = d1chho1
            rv_d2chh  (nv) = d2chh
            rv_d1chx  (nv) = d1chx
          end if
        end if
      end do

    else

      !NECJB Routines corr and dcorr called for any other model

      do iiv = 1, nv_max

        is = tv_gc(iiv)%i
        js = tv_gc(iiv)%j
        iv = tv_gc(iiv)%iv

        !--------------------------------------
        ! distance on the surface of the planet
        !--------------------------------------
        ab(1) = SJ% x(1) - SI% x(1)
        ab(2) = SJ% x(2) - SI% x(2)
        ab(3) = SJ% x(3) - SI% x(3)
        d     = sqrt( (SI% xt(1) - SJ% xt(1))**2 + &
                      (SI% xt(2) - SJ% xt(2))**2 + &
                      (SI% xt(3) - SJ% xt(3))**2     )
        !------------------------------------------------------------
        ! set horizontal correlation functions for height
        ! and derivatives for wind correlations
        ! return from routine for large distances (zero correlations)
        !------------------------------------------------------------
        if (thc) then
          lh = lh_max
        else
          lh = min (SI% L_h , SJ% L_h)
        endif
        if (d < x_cutdr * lh) then
          chh = corr (d, lh)
          call dcorr (d1chh, d1chho1, d2chh, d1chx)
          if (chh /= 0._wp .or. d1chh /= 0._wp) then
            nv = nv + 1

            tv_gc(nv)%i  = is
            tv_gc(nv)%j  = js
            tv_gc(nv)%iv = iv

            rv_ab   (1,nv) = ab(1)
            rv_ab   (2,nv) = ab(2)
            rv_ab   (3,nv) = ab(3)
            rv_d      (nv) = d
            rv_chh    (nv) = chh
            rv_d1chh  (nv) = d1chh
            rv_d1chho1(nv) = d1chho1
            rv_d2chh  (nv) = d2chh
            rv_d1chx  (nv) = d1chx
          end if
        end if
      end do
    end if
    nv_max = nv

FTRACE_END  ("get_corr_opt:3")
FTRACE_BEGIN("get_corr_opt:4")

    do iiv = 1, nv_max

      is = tv_gc(iiv)%i
      js = tv_gc(iiv)%j
!     iv = tv_gc(iiv)%iv
      d  = rv_d (iiv)

      !--------------------------------------------------
      ! set horizontal correlation functions for moisture
      ! default: exponential correlations (rhcexp=.true.)
      !--------------------------------------------------
      lq   = min (SI% L_q , SJ% L_q)
      chq = 0._wp
      if (d < x_cutdr * lq) then
        chq = corrh (d,lq)
      endif
      rv_chq(iiv) = chq
    end do

FTRACE_END  ("get_corr_opt:4")
FTRACE_BEGIN("get_corr_opt:5")

    c_h_vt  = 0._wp
    c_t_vt  = 0._wp
    c_vl_vl = 0._wp
    c_vt_h  = 0._wp
    c_vt_t  = 0._wp
    c_vt_vt = 0._wp
    cti     = 0._wp
    ctj     = 0._wp
    cv      = 0._wp
    cvq     = 0._wp
    cvtt    = 0._wp
    cvx     = 0._wp

    do iiv = 1, nv_max

      is = tv_gc(iiv)%i
      js = tv_gc(iiv)%j
!     iv = tv_gc(iiv)%iv

      ab      = rv_ab   (:,iiv)
      d       = rv_d      (iiv)
      chh     = rv_chh    (iiv)
      chq     = rv_chq    (iiv)
      d1chh   = rv_d1chh  (iiv)
      d1chho1 = rv_d1chho1(iiv)
      d2chh   = rv_d2chh  (iiv)
      d1chx   = rv_d1chx  (iiv)

      clsvv  = lsvv * chh
      rlsvv  = 1._wp / (1._wp + lsvv)
      !----------------------------------------------
      ! set latitude dependent vertical length scales
      !----------------------------------------------
      c_h  = max (SI% cch , SJ% cch)
      c_q  = max (SI% ccq , SJ% ccq)
      c_v  = max (SI% ccv , SJ% ccv)
      !------------------------------------------------------
      ! transformation from local (East/North) coordinates to
      ! longitudinal/transversal coordinates
      !------------------------------------------------------
      if (d==0._wp) ab = SI% du                 ! arbitrary direction for d=0
      elui =   sum (ab * SI% du)                ! longit.vector in u,v system
      elvi =   sum (ab * SI% dv)                !
      abn  =   1._wp / sqrt (elui**2 + elvi**2) ! normalize
      elui =   abn * elui                       !
      elvi =   abn * elvi                       !
      eluj =   sum (ab * SJ% du)                ! same at left hand side
      elvj =   sum (ab * SJ% dv)                !  (changed sign)
      abn  =   1._wp / sqrt (eluj**2 + elvj**2) !
      eluj = - (abn * eluj)                     !
      elvj = - (abn * elvj)                     !
      etui = - elvi                             ! rotate for tang.vector
      etvi =   elui                             !
      etuj = - elvj                             !
      etvj =   eluj

      rv_elui (iiv) = elui
      rv_eluj (iiv) = eluj
      rv_elvi (iiv) = elvi
      rv_elvj (iiv) = elvj
      rv_etui (iiv) = etui
      rv_etuj (iiv) = etuj
      rv_etvi (iiv) = etvi
      rv_etvj (iiv) = etvj
      rv_c_h  (iiv) = c_h
      rv_c_q  (iiv) = c_q
      rv_c_v  (iiv) = c_v
      rv_clsvv(iiv) = clsvv
      rv_rlsvv(iiv) = rlsvv

    end do

FTRACE_END  ("get_corr_opt:5")
FTRACE_BEGIN("get_corr_opt:6_1")

    nll = 0
    nlr = 0
    nil = 0
    njl = 0
    do iiv = 1, nv_max

      is = tv_gc(iiv)%i
      js = tv_gc(iiv)%j
!     iv = tv_gc(iiv)%iv
      SI_n = SI%n
      SJ_n = SJ%n

      il   = 1
      ill  = 0
      do
        ill  = ill + 1
        nil  = max( nil, XI%i )
        il   = il + XI%i
        if (il>SI_n) exit
      end do
      iv_nll(iiv) = ill
      nll = max( nll, ill )

      jl  =  1
      ilr =  0
      do
        ilr  = ilr + 1
        njl  = max( njl, XJ%i )
        jl   = jl + XJ%i
        if (jl>SJ_n) exit
      end do
      iv_nlr(iiv) = ilr
      nlr = max( nlr, ilr )

    end do

    allocate( iv_il(nv_max,nll) )
    allocate( iv_jl(nv_max,nlr) )

FTRACE_END  ("get_corr_opt:6_1")
FTRACE_BEGIN("get_corr_opt:6_2")

    do iiv = 1, nv_max

      is = tv_gc(iiv)%i
      js = tv_gc(iiv)%j

      il   = 1
      do ill = 1, iv_nll(iiv)
        iv_il(iiv,ill) = il
        il=il+SI% col (il)%i    ! XI%i
      end do

      jl  =  1
      do ilr = 1, iv_nlr(iiv)
        iv_jl(iiv,ilr) = jl
        jl=jl+SJ% col (jl)%i    ! XJ%i
      end do

    end do

FTRACE_END  ("get_corr_opt:6_2")
FTRACE_END  ("get_corr_opt:00")
FTRACE_BEGIN("get_corr_opt:99")

    !--------------------------
    ! loop over pairs of levels
    !--------------------------
    do ill = 1, nll
      do ilr = 1, nlr

FTRACE_BEGIN("get_corr_opt:7")

        nv_pol = 0
        do iiv = 1, nv_max
          if( ill > iv_nll(iiv) .or. ilr > iv_nlr(iiv) ) cycle
          nv_pol = nv_pol + 1
          iv_pol(nv_pol) = iiv
        end do

FTRACE_END  ("get_corr_opt:7")

        !=============================
        ! NMC-fitted covariance model:
        !=============================
        if (covm% valid > 2) then

FTRACE_BEGIN("get_corr_opt:8_1")

          nv_o_tv     = 0
          nv_o_tv_not = 0
          nv_a_rh     = 0
          nv_a_uv     = 0
          nv_o_uv     = 0

          do ivp = 1, nv_pol

            iiv = iv_pol(ivp)

            is = tv_gc(iiv)%i
            js = tv_gc(iiv)%j

            cv  (iiv) = 0
            cvq (iiv) = 0
            cvx (iiv) = 0
            cvzl(iiv) = 0
            cvzr(iiv) = 0
            cvzz(iiv) = 0

            il = iv_il(iiv,ill)
            jl = iv_jl(iiv,ilr)
            itl_a_itl = iand (XI% itl, XJ % itl)
            itl_o_itl = ior  (XI% itl, XJ % itl)
            if (iand ( itl_o_itl, OBS_TV ) /= 0) then
              nv_o_tv = nv_o_tv + 1
              iv_o_tv(nv_o_tv) = iiv
            else
              nv_o_tv_not = nv_o_tv_not + 1
              iv_o_tv_not(nv_o_tv_not) = iiv
            endif
            if (iand ( itl_a_itl, OBS_RH) /= 0 .and. &
                rv_chq(iiv) /= 0._wp                    ) then
              nv_a_rh = nv_a_rh + 1
              iv_a_rh(nv_a_rh) = iiv
            endif
            if (iand ( itl_a_itl, ior(OBS_U,OBS_V)) /= 0) then
              nv_a_uv = nv_a_uv + 1
              iv_a_uv(nv_a_uv) = iiv
            endif
            if (iand ( itl_o_itl, ior(OBS_U,OBS_V)) /= 0) then
              nv_o_uv = nv_o_uv + 1
              iv_o_uv(nv_o_uv) = iiv
            endif
          end do

FTRACE_END  ("get_corr_opt:8_1")

          !------------------
          ! height + X, T + T
          !------------------
#if 0
          !NECJB The following loop nest computes a number of dotproducts.
          !      Unfortunately, the vector length covm%nz is small (64),
          !      limiting the performance because of the expensive final
          !      sum collapse.
!NEC$ ivdep
          do ivp = 1, nv_o_tv
            iiv = iv_o_tv(ivp)
            is  = tv_gc(iiv)%i
            js  = tv_gc(iiv)%j
!NEC$ shortloop
            do i = 1, covm% nz
              cv  (iiv) = cv  (iiv) + SI%vm (i, INT_H,  ill) * SJ%vm (i, INT_H,  ilr)
              cvzl(iiv) = cvzl(iiv) + SI%vm (i, INT_TV, ill) * SJ%vm (i, INT_H,  ilr)
              cvzr(iiv) = cvzr(iiv) - SI%vm (i, INT_H , ill) * SJ%vm (i, INT_TV, ilr)
              cvzz(iiv) = cvzz(iiv) + SI%vm (i, INT_TV, ill) * SJ%vm (i, INT_TV, ilr)
            end do
          end do
#endif

          !NECJB The following code computes the sums of the dot products
          !      across the longer dimension of length nv_o_tv, introducing
          !      a strip  mining loop for memory efficiency, while avoiding
          !      sum-collapse  operations.

FTRACE_BEGIN("get_corr_opt:8_2")
          do ivp0 = 0, nv_o_tv-1, MAXVL
            ivlen = min( MAXVL, nv_o_tv-ivp0 )
            ! Compute the products, not the sums
!NEC$ shortloop
            do ivp = 1, ivlen
              iiv = iv_o_tv(ivp0+ivp)
              is  = tv_gc(iiv)%i
              js  = tv_gc(iiv)%j
!NEC$ shortloop
              do i = 1, covm% nz
                rv_sum(ivp,i,1) = SI%vm (i, INT_H,  ill) * SJ%vm (i, INT_H,  ilr)
                rv_sum(ivp,i,2) = SI%vm (i, INT_TV, ill) * SJ%vm (i, INT_H,  ilr)
                rv_sum(ivp,i,3) = SI%vm (i, INT_H , ill) * SJ%vm (i, INT_TV, ilr)
                rv_sum(ivp,i,4) = SI%vm (i, INT_TV, ill) * SJ%vm (i, INT_TV, ilr)
              end do
            end do
            m = covm% nz / 2
            ! Handle odd covm%nz case
            if( m*2 /= covm% nz )then
!NEC$ shortloop
              do ivp = 1, ivlen
                rv_sum(ivp,1,1) = rv_sum(ivp,1,1) + rv_sum(ivp,covm% nz,1)
                rv_sum(ivp,1,2) = rv_sum(ivp,1,2) + rv_sum(ivp,covm% nz,2)
                rv_sum(ivp,1,3) = rv_sum(ivp,1,3) + rv_sum(ivp,covm% nz,3)
                rv_sum(ivp,1,4) = rv_sum(ivp,1,4) + rv_sum(ivp,covm% nz,4)
              end do
            end if

            ! Compute the sums
            do while( m > 0 )
!NEC$ shortloop
              do i = 1, m
!NEC$ shortloop
                do ivp = 1, ivlen
                  rv_sum(ivp,i,1) = rv_sum(ivp,i,1) + rv_sum(ivp,i+m,1)
                  rv_sum(ivp,i,2) = rv_sum(ivp,i,2) + rv_sum(ivp,i+m,2)
                  rv_sum(ivp,i,3) = rv_sum(ivp,i,3) + rv_sum(ivp,i+m,3)
                  rv_sum(ivp,i,4) = rv_sum(ivp,i,4) + rv_sum(ivp,i+m,4)
                end do
              end do
              m = m / 2
            end do
            ! Store final sum
!NEC$ shortloop
            do ivp = 1, ivlen
              iiv = iv_o_tv(ivp0+ivp)
              cv  (iiv) =   rv_sum(ivp,1,1)
              cvzl(iiv) =   rv_sum(ivp,1,2)
              cvzr(iiv) = - rv_sum(ivp,1,3)
              cvzz(iiv) =   rv_sum(ivp,1,4)
            end do
          end do

FTRACE_END  ("get_corr_opt:8_2")
FTRACE_BEGIN("get_corr_opt:8_3")

!NEC$ ivdep
          do ivp = 1, nv_o_tv
            iiv = iv_o_tv(ivp)
            is = tv_gc(iiv)%i
            js = tv_gc(iiv)%j
            il = iv_il(iiv,ill)
            jl = iv_jl(iiv,ilr)
            !------------------
            ! height + X, T + T
            !------------------
            cti (iiv) =  ( XI% ezn * cv  (iiv) + XI%exzn * cvzl(iiv) )
            ctj (iiv) =  ( XJ% ezn * cv  (iiv) - XJ%exzn * cvzr(iiv) )
            cvtt(iiv) =  ( XI% ezn * cv  (iiv) * XJ%ezn  &
                         + XI%exzn * cvzl(iiv) * XJ%ezn  &
                         - XI%ezn  * cvzr(iiv) * XJ%exzn &
                         + XI%exzn * cvzz(iiv) * XJ%exzn )
            c_t_vt(iiv)  = SJ% mu * sq1nu * rv_d1chx(iiv) * cti(iiv) * rv_rlsvv(iiv)
            c_vt_t(iiv)  = SI% mu * sq1nu * rv_d1chx(iiv) * ctj(iiv) * rv_rlsvv(iiv)
          end do

FTRACE_END  ("get_corr_opt:8_3")
FTRACE_BEGIN("get_corr_opt:8_4")

          !-----------
          ! height + X
          !-----------
!NEC$ ivdep
          do ivp = 1, nv_o_tv_not
            iiv = iv_o_tv_not(ivp)
            is = tv_gc(iiv)%i
            js = tv_gc(iiv)%j
!NEC$ shortloop
            do i = 1, covm% nz
              cv(iiv) = cv(iiv) + SI%vm (i, INT_H, ill) * SJ%vm (i, INT_H, ilr)
            end do
          end do

FTRACE_END  ("get_corr_opt:8_4")
FTRACE_BEGIN("get_corr_opt:8_5")

          !--------
          ! RH - RH
          !--------
#if 0
!NEC$ ivdep
          do ivp = 1, nv_a_rh
            iiv = iv_a_rh(ivp)
            is = tv_gc(iiv)%i
            js = tv_gc(iiv)%j
!NEC$ shortloop
            do i = 1, covm% nz
              cvq(iiv) = cvq(iiv) + SI%vm (i, INT_RH, ill) * SJ%vm (i, INT_RH, ilr)
            end do
          end do
#endif
          !NECJB See comment for code section 8_2
          do ivp0 = 0, nv_a_rh-1, MAXVL
            ivlen = min( MAXVL, nv_a_rh-ivp0 )
            ! Compute the products, not the sums
!NEC$ shortloop
            do ivp = 1, ivlen
              iiv = iv_a_rh(ivp0+ivp)
              is  = tv_gc(iiv)%i
              js  = tv_gc(iiv)%j
!NEC$ shortloop
              do i = 1, covm% nz
                rv_sum(ivp,i,1) = SI%vm (i, INT_RH,  ill) * SJ%vm (i, INT_RH,  ilr)
              end do
            end do
            m = covm% nz / 2
            ! Handle odd covm%nz case
            if( m*2 /= covm% nz )then
!NEC$ shortloop
              do ivp = 1, ivlen
                rv_sum(ivp,1,1) = rv_sum(ivp,1,1) + rv_sum(ivp,covm% nz,1)
              end do
            end if

            ! Compute the sums
            do while( m > 0 )
!NEC$ shortloop
              do i = 1, m
!NEC$ shortloop
                do ivp = 1, ivlen
                  rv_sum(ivp,i,1) = rv_sum(ivp,i,1) + rv_sum(ivp,i+m,1)
                end do
              end do
              m = m / 2
            end do
            ! Store final sum
!NEC$ shortloop
            do ivp = 1, ivlen
              iiv = iv_a_rh(ivp0+ivp)
              cvq(iiv) = rv_sum(ivp,1,1)
            end do
          end do

FTRACE_END  ("get_corr_opt:8_5")
FTRACE_BEGIN("get_corr_opt:8_6")

          !------------
          ! wind - wind
          !------------
!NEC$ ivdep
          do ivp = 1, nv_a_uv
            iiv = iv_a_uv(ivp)
            is = tv_gc(iiv)%i
            js = tv_gc(iiv)%j
!NEC$ shortloop
            do i = 1, covm% nz
              cvx(iiv) = cvx(iiv) + SI%vm (i, INT_CHI, ill) * SJ%vm (i, INT_CHI, ilr)
            end do
          end do
!NEC$ ivdep
          do ivp = 1, nv_a_uv
            iiv = iv_a_uv(ivp)
            c_vt_vt(iiv) = ( cv (iiv) * (r1nu * rv_d2chh  (iiv)   - rv_clsvv(iiv)) &
                           + cvx(iiv) *  rnu  * rv_d1chho1(iiv) ) * rv_rlsvv(iiv)
            c_vl_vl(iiv) = ( cv (iiv) * (r1nu * rv_d1chho1(iiv)   - rv_clsvv(iiv)) &
                           + cvx(iiv) *  rnu  * rv_d2chh  (iiv) ) * rv_rlsvv(iiv)
          end do

FTRACE_END  ("get_corr_opt:8_6")
FTRACE_BEGIN("get_corr_opt:8_7")

          !---------
          ! X + wind
          !---------
!NEC$ ivdep
          do ivp = 1, nv_o_uv
            iiv = iv_o_uv(ivp)
            is = tv_gc(iiv)%i
            js = tv_gc(iiv)%j
            c_h_vt(iiv)  = SJ% mu * sq1nu * rv_d1chx(iiv) * cv(iiv) * rv_rlsvv(iiv)
            c_vt_h(iiv)  = SI% mu * sq1nu * rv_d1chx(iiv) * cv(iiv) * rv_rlsvv(iiv)
          end do

FTRACE_END  ("get_corr_opt:8_7")

        else

          !NECJB Code section for covm% valid <= 2
          !NECJB - Not executed in the context of the test case used;
          !        not modified because the new code is not tested.

          !call ftrace_region_begin( 'get_corr_opt_9' )

          do ivp = 1, nv_pol

            iiv = iv_pol(ivp)

            is = tv_gc(iiv)%i
            js = tv_gc(iiv)%j

            il = iv_il(iiv,ill)
            jl = iv_jl(iiv,ilr)
            itl_a_itl = iand (XI% itl, XJ % itl)
            itl_o_itl = ior  (XI% itl, XJ % itl)
            !=================================================
            ! set vertical correlation function for height + X
            !=================================================
            dxh  = XI% xh - XJ% xh ! difference of vert. coord. for height
            dxh2 = dxh * dxh        ! squared differences
            edxh2 = 1._wp / (rv_c_h(iiv) + dxh2)
            eexcv = edem1 * exp ( rv_c_h(iiv) * edxh2)
            cv(iiv)   =  eexcv - edem1
            !========
            ! RH - RH
            !========
            if (iand ( itl_a_itl, OBS_RH) /= 0) then
              if (rv_chq(iiv) /= 0._wp) then
                dxq  = XI% xq - XJ% xq ! difference of vert.coord. for humidity
                dxq2 = dxq * dxq
                cvq(iiv)  =  edem1 * exp ( rv_c_q(iiv) / (rv_c_q(iiv) + dxq2) ) - edem1
              else
                cvq(iiv)  = 0._wp
              endif
            endif
            !============
            ! wind - wind
            !============
            if (iand ( itl_a_itl, ior(OBS_U,ior(OBS_V,OBS_CHI))) /= 0) then
              dxv  = XI% xv - XJ% xv ! difference of vert. coord. for wind
              dxv2 = dxv * dxv
              cvx(iiv)  = (edem1 * exp ( rv_c_v(iiv) / (rv_c_v(iiv)+(dxv2))) - edem1)
              if (modcvx>0._wp) &
                cvx(iiv) = cvx(iiv) * (1._wp-0.75_wp*dxv2) / (1._wp+modcvx*dxv2)
              c_vt_vt(iiv) = ( cv (iiv) * (r1nu * rv_d2chh  (iiv)   - rv_clsvv(iiv)) &
                             + cvx(iiv) *  rnu  * rv_d1chho1(iiv) ) * rv_rlsvv(iiv)
              c_vl_vl(iiv) = ( cv (iiv) * (r1nu * rv_d1chho1(iiv)   - rv_clsvv(iiv)) &
                             + cvx(iiv) *  rnu  * rv_d2chh  (iiv) ) * rv_rlsvv(iiv)
            endif
            !======
            ! T + T
            !======
            if (iand ( itl_o_itl, OBS_TV ) /= 0) then
              cvzl(iiv) = -2._wp * ((rv_c_h(iiv)*dxh) * edxh2**2) &
                      * eexcv
              cvzr(iiv) = cvzl(iiv)
              cvzz(iiv) = ( 2._wp - dxh2*edxh2 * ( 4._wp* rv_c_h(iiv) * edxh2 + 8._wp )) &
                     * eexcv * rv_c_h(iiv) * edxh2 * edxh2
              cti(iiv)  =  ( XI% ezn * cv(iiv) + XI% exzn * cvzl(iiv) )
              ctj(iiv)  =  ( XJ% ezn * cv(iiv) - XJ% exzn * cvzr(iiv) )
              cvtt(iiv) =  ( XI%ezn  * cv(iiv)   * XJ%ezn  &
                      + XI%exzn * cvzl(iiv) * XJ%ezn  &
                      - XI%ezn  * cvzr(iiv) * XJ%exzn &
                      + XI%exzn * cvzz(iiv) * XJ%exzn )
              c_t_vt(iiv)  = SJ% mu * sq1nu * d1chx * cti(iiv) * rv_rlsvv(iiv)
              c_vt_t(iiv)  = SI% mu * sq1nu * d1chx * ctj(iiv) * rv_rlsvv(iiv)
            endif
            !=========
            ! X + wind
            !=========
            if (iand ( itl_o_itl, ior(OBS_U,OBS_V)) /= 0) then
              c_h_vt(iiv)  = SJ% mu * sq1nu * d1chx * cv(iiv) * rv_rlsvv(iiv)
              c_vt_h(iiv)  = SI% mu * sq1nu * d1chx * cv(iiv) * rv_rlsvv(iiv)
            endif
          end do

          !call ftrace_region_end  ( 'get_corr_opt_9' )

        end if

        !----------------------------------------------
        ! loop over pairs of observations at this level
        ! fill correlation matrix
        !----------------------------------------------
FTRACE_BEGIN("get_corr_opt:10_1a")

        allocate( iv_pol1(nv_pol        ,8) )
        allocate( iv_pol2(nv_pol*nil*njl,7) )

        !NECJB Gather data, so that we don't need indirect addressing
        !      repeatedly all nil*njl times.

        do ivp = 1, nv_pol
          iiv = iv_pol(ivp)
          is  = tv_gc(iiv)%i
          js  = tv_gc(iiv)%j
          il  = iv_il(iiv,ill)
          jl  = iv_jl(iiv,ilr)
          iv_pol1(ivp,1) = iiv
          iv_pol1(ivp,2) = tv_gc(iiv)%iv
          iv_pol1(ivp,3) = is
          iv_pol1(ivp,4) = js
          iv_pol1(ivp,5) = il
          iv_pol1(ivp,6) = jl
          iv_pol1(ivp,7) = XI%i
          iv_pol1(ivp,8) = XJ%i
        end do

        !NECJB Build index vector of relevant pairs

        nv_pol2 = 0

        do jj = 1, njl
          do ii = 1, nil
!NEC$ ivdep
            do ivp = 1, nv_pol
              ixi = iv_pol1(ivp,7)
              ixj = iv_pol1(ivp,8)
              if( jj > ixj .or. ii > ixi )cycle
              nv_pol2 = nv_pol2 + 1
              iv_pol2(nv_pol2,5) = ii
              iv_pol2(nv_pol2,6) = jj
              iv_pol2(nv_pol2,7) = ivp
            end do
          end do
        end do

        ! Allocate auxiliary index vectors with sufficient size
        m = nv_pol2
        allocate (iv_rh_rh(m), iv_hs_hs(m), iv_hs_tv(m), iv_hs_u(m), &
                  iv_hs_v (m), iv_tv_hs(m), iv_tv_tv(m), iv_tv_u(m), &
                  iv_tv_v (m), iv_u_hs (m), iv_u_tv (m), iv_u_u (m), &
                  iv_u_v  (m), iv_v_hs (m), iv_v_tv (m), iv_v_u (m), &
                  iv_v_v  (m), iv_ch_ch(m)  )

FTRACE_END  ("get_corr_opt:10_1a")
FTRACE_BEGIN("get_corr_opt:10_1b")

        nv_rh_rh = 0
        nv_hs_hs = 0
        nv_hs_tv = 0
        nv_hs_u  = 0
        nv_hs_v  = 0
        nv_tv_hs = 0
        nv_tv_tv = 0
        nv_tv_u  = 0
        nv_tv_v  = 0
        nv_u_hs  = 0
        nv_u_tv  = 0
        nv_u_u   = 0
        nv_u_v   = 0
        nv_v_hs  = 0
        nv_v_tv  = 0
        nv_v_u   = 0
        nv_v_v   = 0
        nv_ch_ch = 0

!NEC$ ivdep
        do ivp = 1, nv_pol2

          ii  = iv_pol2(ivp,5)
          jj  = iv_pol2(ivp,6)
          jvp = iv_pol2(ivp,7)
          iiv = iv_pol1(jvp,1)
          iv  = iv_pol1(jvp,2)
          is  = iv_pol1(jvp,3)
          js  = iv_pol1(jvp,4)
          il  = iv_pol1(jvp,5)
          jl  = iv_pol1(jvp,6)
          i   = il-1 + ii
          j   = jl-1 + jj
          it  = SI%col(i)% it
          jt  = SJ%col(j)% it
          iv_pol2(ivp,1) = iiv
          iv_pol2(ivp,2) = iv
          iv_pol2(ivp,3) = i
          iv_pol2(ivp,4) = j
          if( it == OBS_RH )then
            if( jt == OBS_RH )then
              nv_rh_rh = nv_rh_rh + 1
              iv_rh_rh(nv_rh_rh) = ivp
            end if
          else if (it == OBS_H .or. it == OBS_HS )then
            if (jt == OBS_H .or. jt == OBS_HS)then
              nv_hs_hs = nv_hs_hs + 1
              iv_hs_hs(nv_hs_hs) = ivp
            else if (jt == OBS_TV)then
              nv_hs_tv = nv_hs_tv + 1
              iv_hs_tv(nv_hs_tv) = ivp
            else if (jt == OBS_U)then
              nv_hs_u  = nv_hs_u  + 1
              iv_hs_u (nv_hs_u) = ivp
            else if (jt == OBS_V)then
              nv_hs_v  = nv_hs_v  + 1
              iv_hs_v (nv_hs_v) = ivp
            end if
          else if (it == OBS_TV)then
            if (jt == OBS_H .or. jt == OBS_HS)then
              nv_tv_hs = nv_tv_hs + 1
              iv_tv_hs(nv_tv_hs) = ivp
            else if (jt == OBS_TV)then
              nv_tv_tv = nv_tv_tv + 1
              iv_tv_tv(nv_tv_tv) = ivp
            else if (jt == OBS_U)then
              nv_tv_u = nv_tv_u + 1
              iv_tv_u(nv_tv_u) = ivp
            else if (jt == OBS_V)then
              nv_tv_v = nv_tv_v + 1
              iv_tv_v(nv_tv_v) = ivp
            end if
          else if (it == OBS_U)then
            if (jt == OBS_H .or. jt == OBS_HS)then
              nv_u_hs = nv_u_hs + 1
              iv_u_hs(nv_u_hs) = ivp
            else if (jt == OBS_TV)then
              nv_u_tv = nv_u_tv + 1
              iv_u_tv(nv_u_tv) = ivp
            else if (jt == OBS_U)then
              nv_u_u = nv_u_u + 1
              iv_u_u(nv_u_u) = ivp
            else if (jt == OBS_V)then
              nv_u_v = nv_u_v + 1
              iv_u_v(nv_u_v) = ivp
            end if
          else if (it == OBS_V)then
            if (jt == OBS_H .or. jt == OBS_HS)then
              nv_v_hs = nv_v_hs + 1
              iv_v_hs(nv_v_hs) = ivp
            else if (jt == OBS_TV)then
              nv_v_tv = nv_v_tv + 1
              iv_v_tv(nv_v_tv) = ivp
            else if (jt == OBS_U)then
              nv_v_u = nv_v_u + 1
              iv_v_u(nv_v_u) = ivp
            else if (jt == OBS_V)then
              nv_v_v = nv_v_v + 1
              iv_v_v(nv_v_v) = ivp
            end if
          else if (it == OBS_CHI)then
            if (jt == OBS_CHI)then
              nv_ch_ch = nv_ch_ch + 1
              iv_ch_ch(nv_ch_ch) = ivp
            end if
          end if
        end do

FTRACE_END  ("get_corr_opt:10_1b")
FTRACE_BEGIN("get_corr_opt:10_1c")

!NEC$ ivdep
        do jv = 1, nv_rh_rh
          ivp = iv_rh_rh(jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = cvq(iiv) * rv_chq(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_hs_hs
          ivp = iv_hs_hs(jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = cv(iiv)  * rv_chh(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_hs_tv
          ivp = iv_hs_tv(jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = - ctj(iiv) * rv_chh(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_hs_u
          ivp = iv_hs_u (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = - c_h_vt(iiv) * rv_etuj(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_hs_v
          ivp = iv_hs_v (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = - c_h_vt(iiv) * rv_etvj(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_tv_hs
          ivp = iv_tv_hs(jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = - cti(iiv) * rv_chh(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_tv_tv
          ivp = iv_tv_tv(jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = cvtt(iiv) * rv_chh(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_tv_u
          ivp = iv_tv_u (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = c_t_vt(iiv) * rv_etuj(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_tv_v
          ivp = iv_tv_v (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = c_t_vt(iiv) * rv_etvj(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_u_hs
          ivp = iv_u_hs (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = - c_vt_h(iiv) * rv_etui(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_u_tv
          ivp = iv_u_tv (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = c_vt_t(iiv) * rv_etui(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_u_u
          ivp = iv_u_u  (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = rv_etui(iiv) * c_vt_vt(iiv) * rv_etuj(iiv) &
                          + rv_elui(iiv) * c_vl_vl(iiv) * rv_eluj(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_u_v
          ivp = iv_u_v  (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = rv_etui(iiv) * c_vt_vt(iiv) * rv_etvj(iiv) &
                          + rv_elui(iiv) * c_vl_vl(iiv) * rv_elvj(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_v_hs
          ivp = iv_v_hs (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = - c_vt_h(iiv) * rv_etvi(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_v_tv
          ivp = iv_v_tv (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = c_vt_t(iiv) * rv_etvi(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_v_u
          ivp = iv_v_u  (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = rv_etvi(iiv) * c_vt_vt(iiv) * rv_etuj(iiv) &
                          + rv_elvi(iiv) * c_vl_vl(iiv) * rv_eluj(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_v_v
          ivp = iv_v_v  (jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = rv_etvi(iiv) * c_vt_vt(iiv) * rv_etvj(iiv) &
                          + rv_elvi(iiv) * c_vl_vl(iiv) * rv_elvj(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do
!NEC$ ivdep
        do jv = 1, nv_ch_ch
          ivp = iv_ch_ch(jv)
          iiv = iv_pol2(ivp,1)
          iv  = iv_pol2(ivp,2)
          i   = iv_pol2(ivp,3)
          j   = iv_pol2(ivp,4)
          c(i,j,iv,i0,j0) = cvx(iiv)
          nonzero(iv,i0,j0) = nonzero(iv,i0,j0) + 1
        end do

FTRACE_END  ("get_corr_opt:10_1c")

        deallocate( iv_pol1 )
        deallocate( iv_pol2 )
        deallocate( &
             iv_rh_rh, iv_hs_hs, iv_hs_tv, iv_hs_u,  iv_hs_v,  iv_tv_hs, &
             iv_tv_tv, iv_tv_u , iv_tv_v , iv_u_hs,  iv_u_tv,  iv_u_u  , &
             iv_u_v  , iv_v_hs , iv_v_tv , iv_v_u ,  iv_v_v ,  iv_ch_ch  )

      end do
    end do

FTRACE_BEGIN("get_corr_opt:11")

    do iiv = 1, nv_max

      is = tv_gc(iiv)%i
      js = tv_gc(iiv)%j
      iv = tv_gc(iiv)%iv

      do j = 1, SJ% n
!NEC$ ivdep
         do i = 1, SI% n
            !-----------------------------------------
            ! return variances instead of correlations
            !-----------------------------------------
            c(i,j,iv,i0,j0) = c(i,j,iv,i0,j0) * SI%col(i)% e * SJ%col(j)% e
         end do
      end do

    end do

FTRACE_END  ("get_corr_opt:11")

    deallocate( iv_il )
    deallocate( iv_jl )

FTRACE_END  ("get_corr_opt:99")

end subroutine get_corr_opt
!==============================================================================
end module mo_fg_cov
