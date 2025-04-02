!
!+ 1-Dimensional Multiresolution (Wavelet) Analysis
!
MODULE mo_1dmra
!
! Description:
!   1-Dimensional Multiresolution (Wavelet) Analysis
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Harald Anlauf
!  Optimisations for NEC SX
! V1_8         2009/12/09 Harald Anlauf
!  optimize for sxf90
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_19        2012-04-16 Harald Anlauf
!  optimized pwt_vec for IBM; cleanup superfluous compiler directives
! V1_43        2015-08-19 Harald Anlauf
!  pwt2: OpenMP parallelization
! V1_45        2015-12-15 Harald Anlauf
!  OpenMP parallelization
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! A. Rhodin  DWD 2005/03 original source
! H. Anlauf  DWD 2005/09 lifted (averaging) wavelets
! C.Dalelane DWD 2006/02 wavelets on the interval
! A. Rhodin  DWD 2006/06 orthogonal (Daubechies) wavelets
! H. Anlauf  DWD 2006/09 Cohen-Daubechies-Vial interval wavelets
! H. Anlauf  DWD 2007/01 Orthogonal Daubechies wavelets D14,D16,D18
! H. Anlauf  DWD 2007/02 Faster transformation routine for orthogonal wavelets
! H. Anlauf  DWD 2008/01 Vectorizing versions of some routines for NEC SX-6
!==============================================================================

  !-------------
  ! Modules used
  !-------------
  use mo_kind,      only: wp     ! working precision kind parameter
  use mo_exception, only: finish ! abort routine
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: mr_gp        ! transformation: direct          (synthesis)
  public :: mr_gp_ad     ! transformation: adjoint
  public :: gp_mr        ! transformation: inverse         (analysis)
  public :: gp_mr_ad     ! transformation: adjoint inverse (dual)
  public :: wave_1d      ! wrapper for transformations above
  public :: wave_1dt     ! wrapper for transformations on transposed x(:,:)
  public :: gp_mr_sh     ! transformation: inverse         (shifted fields)
  public :: wave_1ds     ! wrapper for transformations above
  public :: m_sh         ! returns of 1st dim of y in gp_mr_sh,wave_1ds
  public :: mr_gp_mat    ! matrix transf.: direct
  public :: gp_mr_mat    ! matrix transf.: inverse
  public :: mr_gp_ad_mat ! matrix transf.: adjoint
  public :: wave_1d_odd  ! transformation for arrays with length /= 2**n
  public :: mr_gp_odd    ! odd transformation: direct          (synthesis)
  public :: mr_gp_ad_odd ! odd transformation: adjoint
  public :: gp_mr_odd    ! odd transformation: inverse         (analysis)
  public :: gp_mr_ad_odd ! odd transformation: adjoint inverse (dual)
  public :: homg         ! set homogeneous covariance matrix
  public :: truncate     ! truncate wavelet coefficients
  public :: trunc_mat    ! truncate wavelet coefficient matrix
  public :: trunc_mat_th ! truncate wavelet coefficient matrix
  public :: basis        ! specification of wavelet basis to use
  public :: WV_NONE      !                  identity
  public :: WV_LAZY      !                  Lazy wavelet, just reorder
  public :: WV_HAAR      !                  Haar Wavelet basis
  public :: WV_LIN       !                  continuous wavelet basis
  public :: WV_CUB       !                  basis with continuous derivative
  public :: WV_LIN_LIFTED! linear B-spline (average preserving via lifting)
  public :: WV_CUB_LIFTED! cubic interpolation (average preserving via lifting)
  public :: WV_CUB_B     ! cubic B-spline wavelet, cyclic boundary condition
  public :: WV_QUA_B_INT ! quadratic B-spline, interval boundary condition
  public :: WV_QUA_B     ! quadratic B-spline, cyclic bondary condition
  public :: WV_CUB_B_INT ! cubic B-spline, interval boundary condition
  public :: WV_CUB_LIFT_INT ! cubic interpolation, average preserving, interval
  public :: WV_DAUB_4    ! Daubechies wavelet (N=2)
  public :: WV_DAUB_6    ! Daubechies wavelet (N=3)
  public :: WV_DAUB_8    ! Daubechies wavelet (N=4)
  public :: WV_DAUB_10   ! Daubechies wavelet (N=5)
  public :: WV_DAUB_12   ! Daubechies wavelet (N=6)
  public :: WV_DAUB_14   ! Daubechies wavelet (N=7)
  public :: WV_DAUB_16   ! Daubechies wavelet (N=8)
  public :: WV_DAUB_18   ! Daubechies wavelet (N=9)
  public :: WV_DAUB_20   ! Daubechies wavelet (N=10
  public :: WV_SYM_8     ! Symlet
  public :: WV_SYM_10    ! Symlet
  public :: WV_SYM_12    ! Symlet
  public :: WV_SYM_18    ! Symlet
  public :: WV_CDV_4       ! Cohen-Daubechies-Vial interval wavelet (N=2)
  public :: WV_DAUB_4_INT  ! Alias for CDV_4
  public :: WV_CDV_8       ! Cohen-Daubechies-Vial interval wavelet (N=4)
  public :: WV_DAUB_8_INT  ! Alias for CDV_8
  public :: wv_name      ! name of wavelet
  public :: test_mra     ! test routine: forward vs inverse, linear vs adjoint
  public :: ortho        ! set to .true. by 'test' for orthonormal basis
  public :: factor       ! factorize the array size
  public :: tr_name      ! name of transformation
  public :: TR_SYN, TR_ANA, TR_DUAL, TR_ADJ
  public :: ispowerof2   ! Check integer for a power of 2
  public :: ilog2        ! Integer part of max (log_2(i),0)
  public :: ortho_pwt_version   ! Testing only: 1='pwt', 2='pwt2', 3=vector
  public :: cdv_version         ! Implementation of CDV: 1=scalar, 2=vector
  public :: test_pwt_timing     ! Timing of different pwt implementations
  public :: test_cdv_timing     ! Timing of different CDV implementations
  public :: cdv_clear           ! deallocate module variables

  !==========
  ! Constants
  !==========
  !-----------------
  ! Wavelet basis id
  !-----------------
  integer ,parameter :: WV_NONE         =  1
  integer ,parameter :: WV_LAZY         =  2
  integer ,parameter :: WV_LIN          =  3
  integer ,parameter :: WV_CUB          =  4
  integer ,parameter :: WV_HAAR         =  5
  integer ,parameter :: WV_LIN_LIFTED   =  6
  integer ,parameter :: WV_CUB_LIFTED   =  7
  integer ,parameter :: WV_CUB_B        =  8
  integer ,parameter :: WV_QUA_B_INT    =  9
  integer ,parameter :: WV_QUA_B        =  10
  integer ,parameter :: WV_CUB_B_INT    =  11
  integer ,parameter :: WV_CUB_LIFT_INT =  12
  integer ,parameter :: WV_DAUB_4       =  13
  integer ,parameter :: WV_DAUB_6       =  14
  integer ,parameter :: WV_DAUB_8       =  15
  integer ,parameter :: WV_DAUB_10      =  16
  integer ,parameter :: WV_DAUB_12      =  17
  integer ,parameter :: WV_DAUB_20      =  18
  integer ,parameter :: WV_SYM_8        =  19
  integer ,parameter :: WV_SYM_10       =  20
  integer ,parameter :: WV_SYM_12       =  21
  integer ,parameter :: WV_SYM_18       =  22
  integer, parameter :: WV_CDV_4        =  23
  integer, parameter :: WV_DAUB_4_INT   =  WV_CDV_4 ! Alias
  integer, parameter :: WV_CDV_8        =  24
  integer, parameter :: WV_DAUB_8_INT   =  WV_CDV_8 ! Alias
  integer ,parameter :: WV_DAUB_14      =  25
  integer ,parameter :: WV_DAUB_16      =  26       ! D14..18 were added later
  integer ,parameter :: WV_DAUB_18      =  27
  character(len=16)  :: wv_name(1:27) = (/ &
       'NONE            ','LAZY            ','LIN             ',&
       'CUB             ','HAAR            ','LIN_LIFTED      ',&
       'CUB_LIFTED      ','CUB_B           ','QUA_B_INT       ',&
       'QUA_B           ','CUB_B_INT       ','CUB_LIFT_INT    ',&
       'DAUB_4          ','DAUB_6          ','DAUB_8          ',&
       'DAUB_10         ','DAUB_12         ','DAUB_20         ',&
       'SYM_8           ','SYM_10          ','SYM_12          ',&
       'SYM_18          ','CDV_4           ','CDV_8           ',&
       'DAUB_14         ','DAUB_16         ','DAUB_18         '/)
  ! List of Cohen-Daubechies-Vial interval wavelets
  integer, parameter :: wv_CDV_list (2) = (/ &
       WV_CDV_4, WV_CDV_8 /)
  ! List of wavelets (also) defined on an interval
  integer, parameter :: wv_interval_list (6) = (/ &
       WV_HAAR, WV_QUA_B_INT, WV_CUB_B_INT, WV_CUB_LIFT_INT, &
       wv_CDV_list(:) /)
  !----------------
  ! Transformations
  !----------------
  integer ,parameter :: TR_SYN  =  1  ! Synthesis
  integer ,parameter :: TR_ANA  = -1  ! Analysis
  integer ,parameter :: TR_DUAL =  2  ! Dual       (adjoint analysis)
  integer ,parameter :: TR_ADJ  = -2  ! Adjoint    (adjoint synthesis)
  character(len=4)   :: tr_name(-2:2) = (/'ADJ ','ANA ','XXXX','SYN ','DUAL'/)

  !=================
  ! Module variables
  !=================
  integer :: basis = WV_NONE          ! wavelet basis
  logical :: ortho = .false.
#if defined (__NEC__) || defined(__ibm__)
  integer :: ortho_pwt_version = 3    ! Orthogonal transf.: 3=pwt_vec
#else
  integer :: ortho_pwt_version = 2    ! Orthogonal transf.: 1=pwt, 2=pwt2
#endif

  !-------------------------------------------------
  ! Common constants (not well optimized by sxf90!?)
  !-------------------------------------------------
  real(wp) ,parameter :: sqrt2 = 1.41421356237309504880_wp
  real(wp) ,parameter :: edsq2 = 0.70710678118654752440_wp
  !-------------------------------
  ! coefficients used for WV_CUB_B
  !-------------------------------
! real(wp) ,parameter :: sq2       = 1.4142135623730951_wp
  real(wp) ,parameter :: sq2       = 1._wp
  real(wp) ,parameter :: sq2_64_03 = sq2 / 64 *  3
  real(wp) ,parameter :: sq2_64_05 = sq2 / 64 *  5
  real(wp) ,parameter :: sq2_64_08 = sq2 / 64 *  8
  real(wp) ,parameter :: sq2_64_12 = sq2 / 64 * 12
  real(wp) ,parameter :: sq2_64_32 = sq2 / 64 * 32
  real(wp) ,parameter :: sq2_64_40 = sq2 / 64 * 40
  real(wp) ,parameter :: sq2_64_48 = sq2 / 64 * 48

  !=========================================
  ! module variables for orthogonal wavelets
  !=========================================
  !-------------------------------------------
  ! list of wavelets handled by routine 'pwt':
  !-------------------------------------------
  integer, parameter  :: wv_ortho (13) =                           &
    (/ WV_DAUB_4,  WV_DAUB_6,  WV_DAUB_8,  WV_DAUB_10, WV_DAUB_12, &
       WV_DAUB_14, WV_DAUB_16, WV_DAUB_18, WV_DAUB_20,             &
       WV_SYM_8,   WV_SYM_10,  WV_SYM_12,  WV_SYM_18               /)
  !-----------------------------------------------------------
  ! coefficients (Daubechies, Ten Lectures on Wavelets, p.195)
  ! copied to cc, cr by routine set_coefs
  !-----------------------------------------------------------
  integer, parameter  :: mcof = 20 ! max. support
  real(wp) ,parameter :: daub4(4) =                 &
    (/0.4829629131445341_wp, 0.8365163037378077_wp, &
      0.2241438680420134_wp,-0.1294095225512603_wp /)
  real(wp) ,parameter :: daub6(6) =                 &
    (/0.3326705529500825_wp, 0.8068915093110924_wp, &
      0.4598775021184914_wp,-0.1350110200102546_wp, &
     -0.0854412738820267_wp, 0.0352262918857095_wp /)
  real(wp) ,parameter :: daub8(8) =                 &
    (/0.2303778133088964_wp, 0.7148465705529154_wp, &
      0.6308807679298587_wp,-0.0279837694168599_wp, &
     -0.1870348117190931_wp, 0.0308413818355607_wp, &
      0.0328830116668852_wp,-0.0105974017850690_wp /)
  real(wp) ,parameter :: daub10(10) =               &
    (/0.1601023979741929_wp, 0.6038292697971895_wp, &
      0.7243085284377726_wp, 0.1384281459013203_wp, &
     -0.2422948870663823_wp,-0.0322448695846381_wp, &
      0.0775714938400459_wp,-0.0062414902127983_wp, &
     -0.0125807519990820_wp, 0.0033357252854738_wp /)
  real(wp) ,parameter :: daub12(12) =               &
    (/0.1115407433501095_wp, 0.4946238903984533_wp, &
      0.7511339080210959_wp, 0.3152503517091928_wp, &
     -0.2262646939654400_wp,-0.1297668675672625_wp, &
      0.0975016055873225_wp, 0.0275228655303053_wp, &
     -0.0315820393174862_wp, 0.0005538422011614_wp, &
      0.0047772575109455_wp,-0.0010773010853085_wp /)
  real(wp), parameter :: daub14(14) =                   &
    (/0.07785205408500813_wp, 0.3965393194819123_wp,    &
      0.7291320908462274_wp,  0.4697822874051917_wp,    &
     -0.1439060039285563_wp, -0.2240361849938672_wp,    &
      0.07130921926683042_wp, 0.080612609151082_wp,     &
     -0.03802993693501439_wp,-0.016574541630667_wp,     &
      0.01255099855609955_wp, 0.0004295779729213739_wp, &
     -0.001801640704047446_wp,0.0003537137999745171_wp /)
  real(wp), parameter :: daub16(16) =                     &
    (/0.05441584224310704_wp,   0.3128715909143165_wp,    &
      0.6756307362973218_wp,    0.5853546836542239_wp,    &
     -0.01582910525637238_wp,  -0.2840155429615815_wp,    &
      0.0004724845739030209_wp, 0.1287474266204823_wp,    &
     -0.01736930100181088_wp,  -0.04408825393079791_wp,   &
      0.01398102791739956_wp,   0.00874609404740648_wp,   &
     -0.004870352993451852_wp, -0.000391740373376942_wp,  &
      0.0006754494064506183_wp,-0.0001174767841247786_wp /)
  real(wp), parameter :: daub18(18) =                      &
    (/0.03807794736388813_wp,   0.2438346746126514_wp,     &
      0.6048231236902548_wp,    0.6572880780514298_wp,     &
      0.1331973858249681_wp,   -0.2932737832793372_wp,     &
     -0.0968407832230689_wp,    0.148540749338104_wp,      &
      0.03072568147931585_wp,  -0.06763282906135907_wp,    &
      0.0002509471148277948_wp, 0.02236166212368439_wp,    &
     -0.004723204757752752_wp, -0.004281503682464633_wp,   &
      0.001847646883056686_wp,  0.0002303857635232296_wp,  &
     -0.0002519631889427889_wp, 0.00003934732031628112_wp /)
  real(wp) ,parameter :: daub20(20) =               &
    (/0.0266700579005473_wp, 0.1881768000776347_wp, &
      0.5272011889315757_wp, 0.6884590394534363_wp, &
      0.2811723436605715_wp,-0.2498464243271598_wp, &
     -0.1959462743772862_wp, 0.1273693403357541_wp, &
      0.0930573646035547_wp,-0.0713941471663501_wp, &
     -0.0294575368218399_wp, 0.0332126740593612_wp, &
      0.0036065535669870_wp,-0.0107331754833007_wp, &
      0.0013953517470688_wp, 0.0019924052951925_wp, &
     -0.0006858566949564_wp,-0.0001164668551285_wp, &
      0.0000935886703202_wp,-0.0000132642028945_wp /)
  !------------------------------------------------------
  ! symlets (Daubechies, Ten Lectures on Wavelets, p.198)
  !------------------------------------------------------
  real(wp) ,parameter :: sym8(8) =            &
    (/-0.107148901418_wp, -0.041910965125_wp, &
       0.703739068656_wp,  1.136658243408_wp, &
       0.421234534204_wp, -0.140317624179_wp, &
      -0.017824701442_wp,  0.045570345896_wp /)
  real(wp) ,parameter :: sym10(10) =          &
    (/ 0.038654795955_wp,  0.041746864422_wp, &
      -0.055344186117_wp,  0.281990696854_wp, &
       1.023052966894_wp,  0.896581648380_wp, &
       0.023478923136_wp, -0.247951362613_wp, &
      -0.029842499869_wp,  0.027632152958_wp /)
  real(wp) ,parameter :: sym12(12) =          &
    (/ 0.021784700327_wp,  0.004936612372_wp, &
      -0.166863215412_wp, -0.068323121587_wp, &
       0.694457972958_wp,  1.113892783926_wp, &
       0.477904371333_wp, -0.102724969862_wp, &
      -0.029783751299_wp,  0.063250562660_wp, &
       0.002499922093_wp, -0.011031867509_wp /)
  real(wp) ,parameter :: sym18(18) =          &
    (/ 0.001512487309_wp, -0.000669141509_wp, &
      -0.014515578553_wp,  0.012528896242_wp, &
       0.087791251554_wp, -0.025786445930_wp, &
      -0.270893783503_wp,  0.049882830959_wp, &
       0.873048407349_wp,  1.015259790832_wp, &
       0.337658923602_wp, -0.077172161097_wp, &
       0.000825140929_wp,  0.042744433602_wp, &
      -0.016303351226_wp, -0.018769396836_wp, &
       0.000876502539_wp,  0.001981193736_wp /)

  real(wp) :: cc(mcof), cr(mcof) ! coefficients copied by set_coefs
  integer  :: basis_ortho = 0    ! coefficients are copied for this basis
  integer  :: ncof               ! support of the wavelet
  integer  :: ioff, joff         ! center  of the wavelet
  !======================================================================
  ! Variables and parameters of the Cohen-Daubechies-Vial interval wavelets
  !{{{
  integer :: cdv_basis = 0              ! Current CDV basis
  logical :: cdv_precondition = .true.  ! Default: apply preconditioner to data
  logical :: cdv_keep_scaling = .false. ! Default: extended transformation
#if defined (__SX__) && !defined (__NEC__)
  integer :: cdv_version = 2            ! Implementation: 1=scalar, 2=vector
#else
  integer :: cdv_version = 1            ! Implementation: 1=scalar, 2=vector
#endif

  public  :: cdv_precondition

  ! Parameters of interval wavelet transform
  type cdv_parm_t
     integer :: N   = 0                           ! Wavelet order
     integer :: lf  = 0                           ! Filter length (grid points)
     integer :: lef = 0                           ! Edge filter length
!     real(wp), pointer :: c_fwd(:)    => NULL ()  ! Interior filter: low-pass
!     real(wp), pointer :: c_inv(:)    => NULL ()  ! Interior filter: high-pass
!     real(wp), pointer :: le_hi(:,:)  => NULL ()  ! Left-edge high-pass
!     real(wp), pointer :: le_lo(:,:)  => NULL ()  ! Left-edge low-pass
!     real(wp), pointer :: re_hi(:,:)  => NULL ()  ! Right-edge high-pass
!     real(wp), pointer :: re_lo(:,:)  => NULL ()  ! Right-edge low-pass
!     real(wp), pointer :: l_pre(:,:)  => NULL ()  ! Left-edge preconditioner
!     real(wp), pointer :: r_pre(:,:)  => NULL ()  ! Right-edge preconditioner
!     real(wp), pointer :: l_post(:,:) => NULL ()  ! Left-edge postconditioner
!     real(wp), pointer :: r_post(:,:) => NULL ()  ! Right-edge postconditioner
  end type cdv_parm_t
  type (cdv_parm_t), save :: cdv

  ! Performance optimization: avoid pointers in favor of regular allocatables
  real(wp), allocatable :: cdv_c_fwd(:)         ! Interior filter: low-pass
  real(wp), allocatable :: cdv_c_inv(:)         ! Interior filter: high-pass
  real(wp), allocatable :: cdv_le_hi(:,:)       ! Left-edge high-pass
  real(wp), allocatable :: cdv_le_lo(:,:)       ! Left-edge low-pass
  real(wp), allocatable :: cdv_re_hi(:,:)       ! Right-edge high-pass
  real(wp), allocatable :: cdv_re_lo(:,:)       ! Right-edge low-pass
  real(wp), allocatable :: cdv_l_pre(:,:)       ! Left-edge preconditioner
  real(wp), allocatable :: cdv_r_pre(:,:)       ! Right-edge preconditioner
  real(wp), allocatable :: cdv_l_post(:,:)      ! Left-edge postconditioner
  real(wp), allocatable :: cdv_r_post(:,:)      ! Right-edge postconditioner

  ! "Joining coefficients" from CDV_4 to the top-level Haar wavelet:
  real(wp), parameter :: cdv_lo42(2,4) = reshape ( (/ &
       0.3789140420125741_wp,  0.7729271227102750_wp, &
       0.5074585252828727_wp, -0.0386478571494032_wp, &
       -.1229045724050927_wp, -0.0647982858598916_wp, &
       0.2631768878482613_wp,  0.9546903026291012_wp  &
       /), shape=(/2,4/), order=(/2,1/) )

  real(wp), parameter :: cdv_hi42(2,4) = reshape ( (/ &
       0.9172342202772499_wp, -0.3279826888409800_wp, &
       -.1743693318075681_wp,  0.1438891138968081_wp, &
       0.0000000000000000_wp, -0.5392700630936935_wp, &
       0.8017600058279397_wp, -0.2576212182754090_wp  &
       /), shape=(/2,4/), order=(/2,1/) )

  ! "Joining coefficients" from CDV_8 to CDV_4:
  real(wp), parameter :: cdv_lo84(4,8) = reshape ( (/ &
       0.9060684522380685_wp, -0.0725170422067634_wp, -0.3119443371598975_wp, &
       -.1975872818331880_wp, -0.0066801499990693_wp,  0.1477351993215146_wp, &
       0.1162059305966940_wp, -0.0453552994887776_wp, &
       0.3105461732489494_wp,  0.7219558247198944_wp,  0.4608970022634533_wp, &
       0.1201638196353584_wp, -0.1563604812404894_wp, -0.2769102357548674_wp, &
       -.1049801641250508_wp,  0.2081508048731922_wp, &
       0.1268186839637496_wp,  0.0816443566048380_wp,  0.2998162152358630_wp, &
       0.4665290158481677_wp,  0.5360706243222795_wp,  0.4028109377588095_wp, &
       -.0041708949163168_wp, -0.4691123241632588_wp, &
       -.0418417503784564_wp,  0.0190198567142342_wp,  0.0618958398753194_wp, &
       0.1194393854892175_wp,  0.2178038586912327_wp,  0.3797662896412012_wp, &
       0.5183629897673542_wp,  0.7207145804038784_wp  &
       /), shape=(/4,8/), order=(/2,1/) )

  real(wp), parameter :: cdv_hi84(4,8) = reshape ( (/ &
       0.2544942503521722_wp, -0.6603426991244594_wp,  0.4089705552497403_wp, &
       0.3439926257385763_wp, -0.0167416518388335_wp, -0.3263675149741403_wp, &
       -.1983195185871304_wp,  0.2597415743341619_wp, &
       0.0000000000000000_wp, -0.1693580653517153_wp,  0.6551999014384275_wp, &
       -.6180831179653008_wp, -0.1486183582605626_wp,  0.2486443535349805_wp, &
       0.2059030497252578_wp, -0.1835666018060904_wp, &
       0.0000000000000000_wp,  0.0359466528169432_wp,  0.0000000000000000_wp, &
       -.4384081725890520_wp,  0.7863388843633551_wp, -0.3935598511689838_wp, &
       -.1259295021008718_wp,  0.1320214957688865_wp, &
       0.0000000000000000_wp, -0.0216035504572214_wp,  0.0000000000000000_wp, &
       0.1484770689534628_wp,  0.0000000000000000_wp, -0.5213317754197064_wp, &
       0.7804952430805161_wp, -0.3106898808604068_wp  &
       /), shape=(/4,8/), order=(/2,1/) )
  !}}}
  !======================================================================

contains
!==============================================================================
  subroutine wave_1ds (y, x, trans, base, keep2, shift)
  real(wp) ,intent(out)          :: y(:,:) ! result
  real(wp) ,intent(in)           :: x(:,:) ! array to transform
  integer  ,intent(in)           :: trans  ! transformation
  integer  ,intent(in) ,optional :: base   ! wavelet basis
  integer  ,intent(in) ,optional :: keep2  ! flag to keep a factor of 2
  logical  ,intent(in) ,optional :: shift  ! transform shifted fields
  !----------------------------------------------------------------------
  ! As subroutine WAVE_1D described below, but with the option (if SHIFT
  ! is passed with value .true.) to calculate all coefficients needed to
  ! represent the fields X, shifted by an arbitrary number of gridpoints,
  ! in wavelet representation.  Because the number of coefficients
  ! returned differs from the number of gridpoints, the transformation
  ! cannot be performed in place and Y is the the result argument. For
  ! details refer to the specific subroutine gp_mr_sh.
  !----------------------------------------------------------------------
    logical :: sh
    integer :: nx      ! array dimension
    integer :: n2      ! number of factors 2
    integer :: nr      ! remaining factor
    integer :: oldbase

    sh = .false.; if(present(shift)) sh=shift

    if (sh) then
      oldbase = basis
      if (present(base)) basis = base
      if (basis == WV_NONE) then
        basis = oldbase
        return
      endif
      nx = size(x,1)                   ! array size
      call factor (nx, n2, nr ,keep2)  ! factorize array size
      call gp_mr_sh (x, y, n2)
      basis = oldbase
    else
      y = x
      call wave_1d (y, trans, base, keep2)
    endif

  end subroutine wave_1ds
!------------------------------------------------------------------------------
  function m_sh (nx, keep2, shift)
  !---------------------------------------------
  ! returns the size of the first dimension of y
  ! for routines wave_1ds or gp_mr_sh
  !---------------------------------------------
  integer                       :: m_sh   ! first dimension for shifted fields
  integer ,intent(in)           :: nx     ! first dimension for regular fields
  integer ,intent(in) ,optional :: keep2  ! flag to keep a factor of 2
  logical ,intent(in) ,optional :: shift  ! transform shifted fields
    logical :: sh
    integer :: n2, nr
    sh = .true.; if (present (shift)) sh = shift
    if (sh) then
      call factor (nx, n2, nr ,keep2)
      m_sh = (n2+1)*nx
    else
      m_sh = nx
    endif
  end function m_sh
!------------------------------------------------------------------------------
  subroutine wave_1dt (x, trans, base, keep2)
  real(wp) ,intent(inout)        :: x(:,:) ! array to transform
  integer  ,intent(in)           :: trans  ! transformation
  integer  ,intent(in) ,optional :: base   ! wavelet basis
  integer  ,intent(in) ,optional :: keep2  ! flag to keep a factor of 2
    real(wp) :: xt (size(x,2)  ,size(x,1))
#if defined (__NEC__)
    integer  :: n1, n2
    real(wp) :: x_ (size(x,1)+1,size(x,2)) ! Auxiliary arrays to reduce
    real(wp) :: xt_(size(x,2)+1,size(x,1)) ! bank conflicts on transpose ()
    !----------------------------------------------------------------------
    ! NB: We introduce a significant overhead to transpose arrays here,
    !     and we do it again in the specific transform routines.
    !     Someday we could/should directly call the vector versions
    !     of the transform routines applied to the transposed arrays!
    !----------------------------------------------------------------------
    n1 = size(x,1)
    n2 = size(x,2)
    x_ (1:n1,1:n2) = x (1:n1,1:n2)
    xt (1:n2,1:n1) = transpose (x_ (1:n1,1:n2))
    call wave_1d (xt, trans, base, keep2)
    xt_(1:n2,1:n1) = xt(1:n2,1:n1)
    x(1:n1,1:n2)   = transpose (xt_(1:n2,1:n1))
#else
    xt = transpose (x)
    call wave_1d (xt, trans, base, keep2)
    x  = transpose (xt)
#endif
  end subroutine wave_1dt
!------------------------------------------------------------------------------
  subroutine wave_1d (x, trans, base, keep2)
  real(wp) ,intent(inout)        :: x(:,:) ! array to transform
  integer  ,intent(in)           :: trans  ! transformation
  integer  ,intent(in) ,optional :: base   ! wavelet basis
  integer  ,intent(in) ,optional :: keep2  ! flag to keep a factor of 2
  !-----------------------------------------------------------------------
  ! Generic wavelet transform driver routine. The row vectors of array X
  ! (first index) are transformed in place. The direction of the transform
  ! is specified by the value of TRANS. The wavelet basis specified by the
  ! optional parameter BASE overwrites the default value given by module
  ! variable BASIS. If the optional parameter KEEP2 is passed with the
  ! value .true. an additional factor of 2 is kept in the number of
  ! scaling functions on the coarsest scale.
  !-----------------------------------------------------------------------

    integer :: nx      ! array dimension
    integer :: n2      ! number of factors 2
    integer :: nr      ! remaining factor
    integer :: oldbase

    oldbase = basis
    if (present(base)) basis = base
    if (basis == WV_NONE) then
      basis = oldbase
      return
    endif

    if (any (basis == wv_CDV_list(:))) then
       call wave_cdv (x, trans, basis)
    else
       nx = size(x,1)                   ! array size
       call factor (nx, n2, nr ,keep2)  ! factorize array size

       select case (trans)
       case (TR_SYN)
          call mr_gp    (x, n2)
       case (TR_ANA)
          call gp_mr    (x, n2)
       case (TR_DUAL)
          call gp_mr_ad (x, n2)
       case (TR_ADJ)
          call mr_gp_ad (x, n2)
       case default
          write (0,*)  'wave_1d:  invalid transform =',trans
          call finish ('wave_1d','invalid transform')
       end select
    end if

    basis = oldbase

  end subroutine wave_1d
!------------------------------------------------------------------------------
  subroutine mr_gp (x, n2)
  real(wp) ,intent(inout) :: x(:,:) ! array to transform
  integer  ,intent(in)    :: n2     ! number of factors of 2
  !-------------------------------------------------------------
  ! Direct wavelet transform (wavelet -> gridpoint)
  ! Wavelet is specivied by module variable 'basis'.
  ! all columns (first index) of x are transformed independently
  !-------------------------------------------------------------

    integer           :: nr   ! remaining factor
    integer           :: l    ! loop count

    if (basis == WV_NONE) return

    nr = size(x,1) / 2**n2

    do l=1,n2                 ! loop coarse -> fine scale
      call down    (x,nr)     ! direct wavelet transform
      call scatter (x,nr)     ! reorder coefficients
      nr = nr * 2             ! remaining factor in next iteration
    end do

  end subroutine mr_gp
!------------------------------------------------------------------------------
  subroutine gp_mr (x, n2)
  real(wp) ,intent(inout) :: x(:,:) ! array to transform
  integer  ,intent(in)    :: n2     ! number of factors of 2
  !-------------------------------------------------
  ! Inverse wavelet transform (gridpoint -> wavelet)
  !-------------------------------------------------
    integer           :: nr   ! remaining factor
    integer           :: l    ! loop count

    if (basis == WV_NONE) return

    nr = size(x,1)            ! start with full length on finest scale

    do l=1,n2                 ! loop fine  ->  coarse scale
      nr = nr / 2             ! size in this iteration
      call gather (x,nr)      ! reorder coefficients
      call up     (x,nr)      ! inverse wavelet transform
    end do

  end subroutine gp_mr
!------------------------------------------------------------------------------
  subroutine gp_mr_sh (x, y, n2)
  real(wp) ,intent(in)  :: x(:,:) ! array to transform
  real(wp) ,intent(out) :: y(:,:) ! output
  integer  ,intent(in)  :: n2     ! number of factors of 2
  !-------------------------------------------------------------------------
  ! Inverse wavelet transform (gridpoint -> wavelet).
  ! All wavelet coefficients of shifted fields x will be calculated.
  ! The size of the result y (first index) is the size of x (nx) times n2+1.
  ! n2 is the number of different scales in the wavelet representation.
  ! (number of occurence of the factor 2 in nx.
  ! Ordering in y is as follows:
  !
  ! s_0..s_2*n2-1 | w_1_0..w_1_2*n2-1 | w_2_0..w_2_n2-1 |..| w_n2_0 w_n2_1
  ! 1          nx                2*nx              3*nx          (n2+1)*nx
  !
  ! with s_i   = scaling function   of the field x shifted by i gridpoints
  !      w_l_i = wavelet on scale l of the field x shifted by i gridpoints
  !
  ! for instance, transformation with the lazy wavelet gives the following
  ! result:
  !
  !   nx = 12
  !   n2 =  2
  !   x  = (1 2 3 4 5 6 7 8 9 a b c)
  !
  !   y  = (1 2 3 4 5 6 7 8 9 a b c|3 4 5 6 7 8 9 a b c 1 2|2 3 4 5 6 7 8 9 a b c 1)
  !   i  = (1                    nx                      2nx                     3nx
  !   scale 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2
  !
  !   scale=0 denotes the scaling functions. The cyclic shift for the wavelets
  !   on each scale results from the center of the wavelet basis functions is
  !   shifted with respect to the center of the scaling functions.
  !
  ! The operation count of this transform is O(nx*log(nx)).
  !-------------------------------------------------------------------------
    real(wp) ,allocatable :: z (:,:) ! temporary
    integer               :: nx      ! array dimension
    integer               :: ny      ! array dimension
    integer               :: nr      ! remaining factor
    integer               :: l       ! loop count
    integer               :: m       ! loop count
    integer               :: k, n    ! loop count, stride
    !-------------
    ! factorize nx
    !-------------
    nx = size(x,1)            ! array size
    ny = size(x,2)            ! array size

    !----------------
    ! allocate arrays
    !----------------
    allocate   (z ( 2*nx, ny))
    if (size(y,1) /= (n2+1)*nx) then
      write(0,*) 'gp_mr_sh: x ',shape(x)
      write(0,*) 'gp_mr_sh: y ',shape(y)
      write(0,*) 'gp_mr_sh: n2',n2
      write(0,*) 'gp_mr_sh: size(y,1) should be (n2+1)*size(x,1)'
      call finish ('gp_mr_sh','wrong array size')
    endif

    nr = nx                               ! start with full length,finest scale

    z(1:nx,:) = x
    do l=1,n2                             ! loop fine  ->  coarse scale
      do m = 0, nx-1, nr
        z (nx+m+1:nx+m+nr,:) = cshift (z(m+1:m+nr,:),1) ! shifted copy
      end do

      nr = nr / 2                         ! size in this iteration

      do m = 0, 2*nx-1, 2*nr              ! loop over shifted parts
        call gather (z(m+1:m+2*nr,:),nr)  ! reorder coefficients
        call up     (z(m+1:m+2*nr,:),nr)  ! inverse wavelet transform
      end do

      n = nx/nr
      do m = 0, nx-1, nr                  ! reorder parts
        k = m/nr
        y ((n2-l+1)*nx+k+1:(n2-l+1)*nx  +nx:n,:) = z (2*m+nr+1:2*(m+nr),:)
        z (            m+1:            m+nr  ,:) = z (2*m+   1:2* m+nr ,:)
      end do
    end do

!   y (1:nx,:) = z (1:nx,:)

    do m = 0, nx-1, nr                  ! reorder parts
      k = m/nr
      y (k+1:nx:n,:) = z (m+1:m+nr,:)
    end do

    deallocate (z)

  end subroutine gp_mr_sh
!------------------------------------------------------------------------------
  subroutine mr_gp_ad (x, n2)
  real(wp) ,intent(inout) :: x(:,:) ! array to transform
  integer  ,intent(in)    :: n2     ! number of factors of 2
  !--------------------------
  ! Adjoint wavelet transform
  !--------------------------
    integer           :: nr   ! remaining factor
    integer           :: l    ! loop count

    if (basis == WV_NONE) return

    nr = size(x,1)            ! start with full length on finest scale

    do l=1,n2                 ! loop fine  ->  coarse scale
      nr = nr / 2             ! size in this iteration
      call gather (x,nr)      ! reorder coefficients
      call downad (x,nr)      ! adjoint wavelet transform
    end do

  end subroutine mr_gp_ad
!------------------------------------------------------------------------------
  subroutine gp_mr_ad (x, n2)
  real(wp) ,intent(inout) :: x(:,:) ! array to transform
  integer  ,intent(in)    :: n2     ! number of factors of 2
  !---------------------------------------------------------
  ! Adjoint inverse wavelet transform (wavelet -> gridpoint)
  ! used to visualize the dual wavelet basis
  !---------------------------------------------------------

    integer           :: nr   ! remaining factor
    integer           :: l    ! loop count

    if (basis == WV_NONE) return

    nr = size(x,1) / 2**n2

    do l=1,n2                 ! loop coarse -> fine scale
      call upad    (x,nr)     ! adjoint inverse wavelet transform
      call scatter (x,nr)     ! reorder coefficients
      nr = nr * 2             ! remaining factor in next iteration
    end do

  end subroutine gp_mr_ad
!==============================================================================
  subroutine test_mra (nx, ny)
    integer, intent(in), optional :: nx, ny     ! Field size
  !-----------------------------------------
  ! test: forward vs inverse transform
  !       linear  vs adjoint transform
  ! This routine sets module variable ortho.
  !-----------------------------------------
    real(wp), dimension(:,:), allocatable :: x, y, w, g
    real(wp)         :: t      ! scalar quantity to check (rounding error size)
    integer          :: n      ! array size
    integer          :: n1, n2 ! Sizes of fictitious "horizontal" dimensions
    character(len=4) :: c      ! string to print

    if (present (nx)) then
       n1 = nx
    else
       n1 = 256
    end if
    if (present (ny)) then
       n2 = ny
    else
       n2 = n1/2        ! Default grid with aspect ratio 2:1
    end if
    allocate (x(n1,n2), y(n1,n2), w(n1,n2), g(n1,n2))

    n = size (x)
    write(6,*) 'test wavelet basis     : ',wv_name(basis),' (',basis,')'
    call random_number (x)
    call random_number (y)
    !-------------------------------
    ! test direct vs inverse routine
    ! both directions
    !-------------------------------
    w = x
    call wave_1d (w ,TR_ANA)
    g = w
    call wave_1d (g ,TR_SYN)
    t = sqrt(sum((g-x)**2)/sum(x**2))          ! norm of difference
    c = 'FAIL'; if (t<1.e-12_wp) c = 'OK'
    write(6,*) 'test forward vs inverse:',t,c  ! printout

    g = x
    call wave_1d (g ,TR_SYN)
    w = g
    call wave_1d (w ,TR_ANA)
    t = sqrt(sum((w-x)**2)/sum(x**2))
    c = 'FAIL'; if (t<1.e-12_wp) c = 'OK'
    write(6,*) 'test forward vs inverse:',t,c  ! printout

    !--------------------------------
    ! test linear  vs adjoint routine
    !--------------------------------
    g = x                                      ! x     in wavelet space
    call wave_1d (g ,TR_SYN)                   ! W (x) in grid space
    w = y                                      ! y     in grid space
    call wave_1d (w ,TR_ADJ)                   ! Wt(y) in wavelet space
    t = (sum(g*y)-sum(x*w)) / sum(x*w)         ! (W x, y) - (x, Wt y)
    c = 'FAIL'; if (abs(t)<1.e-12_wp) c = 'OK'
    write(6,*) 'test forward vs adjoint:',t,c  ! printout

    !-----------------------------
    ! test adjoint inverse routine
    !-----------------------------
    g = x                                      ! x     in wavelet space
    call wave_1d (g ,TR_DUAL)                  ! W-1t (x) in grid space
    w = y                                      ! y        in grid space
    call wave_1d (w ,TR_ANA)                   ! W-1  (y) in wavelet space
    t = (sum(g*y)-sum(x*w)) / sum(x*w)         ! (W-1t x, y) - (x, W-1 y)
    c = 'FAIL'; if (abs(t)<1.e-12_wp) c = 'OK'
    write(6,*) 'test inverse vs adjoint:',t,c  ! printout

    !--------------------
    ! test orthonormality
    ! W^-1 = Wt ?
    !--------------------
    g = x                                      ! x     in wavelet space
    call wave_1d (g ,TR_SYN)                   ! W (x) in grid space
    w = y                                      ! y      in grid space
    call wave_1d (w ,TR_ANA)                   ! W-1(y) in wavelet space
    t = (sum(g*y)-sum(x*w)) / sum(x*w)         ! (W x, y) - (x, W-1(y))
    ortho = abs(t)<1.e-12_wp                   ! set module variable 'ortho'
    c = 'FAIL'; if (ortho) c = 'OK'
    write(6,*) 'test orthonormality    :',t,c  ! printout

    !-----------------------------------------------------
    ! Report the number of vanishing moments of the filter
    ! (for periodic orthogonal and interval wavelets only)
    !-----------------------------------------------------
    if (any (basis == wv_interval_list(:))) then
       write (6,*) "Vanishing moments      :", vanishing_moment_count ()
    else if (any (basis == wv_ortho(:))) then
       write (6,*) "Vanishing moments      :", filter_vanishing_moments ()
    end if

  end subroutine test_mra
!==============================================================================
  subroutine truncate (x, n)
  !---------------------------------------------------------
  ! truncate wavelet coefficients
  ! independently for each column
  ! leave coarsest scale unchanged
  ! for the other scales keep 2 n coefficients at the center
  !---------------------------------------------------------
  real(wp)        ,intent(inout)       :: x(:,:)  ! coefficients to truncate
  integer         ,intent(in)          :: n       ! keep 2n coefficients/scale

    integer           :: nx       ! array dimension
    integer           :: n2       ! number of factors 2
    integer           :: nr       ! remaining factor
    integer           :: l        ! indices

    nx = size(x,1)                ! array size
    call factor (nx, n2, nr)      ! factorize array size

    if (n<1) return               ! do nothing for n<1

    do l=1,n2                     ! loop: coarse -> fine scale
      x(nr+1+n:2*nr-n,:) = 0._wp  ! zero coefficients
      nr = nr * 2                 ! array size in next iteration
    end do
  end subroutine truncate
!------------------------------------------------------------------------------
  subroutine truncatex (x, thresh, comment)
  !-----------------------------------------------------------------
  ! truncate wavelet coefficients
  ! independently for each column
  ! leave coarsest scale unchanged
  ! for the other scales keep coefficients >= thresh * maxval(scale)
  !-----------------------------------------------------------------
  real(wp)        ,intent(inout)       :: x(:,:)  ! coefficients to truncate
  real(wp)        ,intent(in)          :: thresh  ! threshold: 0 .. 1
  character(len=*),intent(in),optional :: comment !

    integer           :: nx                ! array dimension
    integer           :: n2                ! number of factors 2
    integer           :: nr                ! remaining factor
    integer           :: j, l              ! indices
    real(wp) ,pointer :: p (:)             ! pointer to coefs.
    real(wp)          :: xthr              ! absolute threshold value
    real(wp)          :: t (size(x,1)/2)   ! temporary
    target            :: x

    nx = size(x,1)                         ! array size
    call factor (nx, n2, nr)               ! factorize array size

    if (present(comment)) then             ! print header
      write (6,'()')
      write (6,'(a,a,a,f5.2)') '  truncating ',trim(comment),' thresh =',thresh
      write (6,'()')
      write (6,'(a)') '   l, coefs,  kept'
    endif
    if (thresh <= 0._wp) return            ! do nothing for thresh <= 0._wp
    do l=1,n2                              ! loop: coarce -> fine scale
      do j=1,size(x,2)                     ! loop: columns
        p => x(nr+1:2*nr,j)                ! pointer => coefs. of this scale
        t (:nr) = abs(p)                   ! absolute values of coefs.
        xthr    = maxval(t(:nr)) * thresh  ! absolute value of threshold
        if (present(comment)) write (6,'(i4,2i7)') l,nr, count (t(:nr)>=xthr)
        where (t(:nr) < xthr) p = 0._wp    ! set coefficients to zero
      end do
      nr = nr * 2                          ! array size of next iteration
    end do
    if (present(comment)) write (6,'()')   ! print trailer
  end subroutine truncatex
!==============================================================================
  subroutine trunc_mat_th (A, thresh, scale)
  !------------------------------------------------------------------
  ! truncate wavelet coefficients of a matrix
  ! leave coarsest scale unchanged
  ! for the other scales keep coefficients >= thresh * maxval(matrix)
  !------------------------------------------------------------------
  real(wp) ,intent(inout) :: A (:,:)          ! matrix
  real(wp) ,intent(in)    :: thresh           ! relative threshold
  logical  ,intent(in)    :: scale            ! scale with diagonal

    integer  :: n                             ! array size
    integer  :: n2                            ! number of factors 2
    integer  :: nr                            ! remaining factor
    integer  :: i,j                           ! indices
    real(wp) :: rn                            ! array size (float)
    real(wp) :: maxv                          ! max. value of coefficient
    real(wp) :: th                            ! absolute value of threshold

    if (thresh <= 0._wp) return               ! thresh <= 0 : do nothing
    n  = size(A,1)                            ! array size
    rn = n
    call factor (n, n2, nr)                   ! factorize array size
    maxv = maxval(abs(A))                     ! max. value of coefficient
    th = thresh * maxv                        ! absolute value of threshold

    do i=1,n                                  ! loop over rows
      do j=1,n                                ! loop over columns
        if (i<=nr .and. j<=nr) cycle          ! leave coarsest scale unchanged
        if (i==j)              cycle          ! leave diagonal unchanged
        if (scale) &                          ! scale with diagonal
          maxv = sqrt(abs(A(i,i)*A(j,j))) * thresh
        if (abs(A(i,j)) < th)  A(i,j) = 0._wp ! set small coefs. to zero
      end do
    end do

    write (6,*) 'trunc_mat_th:',count (A/=0._wp)/rn,'pixels out of',n,'kept,',&
      'thresh =',thresh,', maxval=',maxv

  end subroutine trunc_mat_th
!==============================================================================
  subroutine trunc_mat (A, nod)
  !------------------------------------------
  ! truncate wavelet coefficients of a matrix
  ! leave certain coefficients unchanged
  !------------------------------------------
  real(wp) ,intent(inout) :: A (:,:) ! matrix to change
  integer  ,intent(in)    :: nod     ! number of off-diagonals to keep

    integer           :: nx                 ! array dimension
    integer           :: n2, n2x, n2y       ! number of factors 2
    integer           :: nr, nrx, nry, nrm  ! remaining factor
    integer           :: m                  ! indices
    integer           :: fx,fy,ix,iy,ixy,lx,ly,iix,iiy
    integer           :: ix1,ix2,iy1,iy2
    real(wp) ,pointer :: p (:,:)
    target            :: A

    if (nod < 0._wp) return               ! nod < 0: do nothing
    nx = size(A,1)
    call factor (nx, n2, nr)

    nrx = nr
    do lx = 0, n2
      if (lx==0) then
        n2x = 1
        ix1 = 1
        ix2 = nrx
      else
        n2x = lx
        ix1 = nrx+1
        ix2 = nrx*2
      endif


      nry = nr
      do ly = 0, n2
        if (ly==0) then
          n2y = 1
          iy1 = 1
          iy2 = nry
        else
          n2y = ly
          iy1 = nry+1
          iy2 = nry*2
        endif

!       m = 0
!
!       if (lx == 0 .and. ly == 0) then
!         m = nx
!       else if (lx == ly) then
!         m = 2
!       else
!         if (lx == 0 .or. ly == 0) then
!           m = 2
!         else ! if (abs(lx-ly)<=2) then
!           m = 2
!         endif
!       endif

        m = nod

        nrm = min (nrx,nry)
        fx = nrx / nrm
        fy = nry / nrm

        if (m<nrm) then

          p => A (ix1:ix2, iy1:iy2)

          if (m==0) then
            p = 0._wp
          else
            do ix = 1, nrx
              do iy = 1,nry

                iix = (ix-1)/fx
                iiy = (iy-1)/fy
                ixy = min( mod (nrm+iix-iiy,nrm), mod (nrm+iiy-iix,nrm)) + 1
                if (ixy > m) p (ix, iy) = 0._wp

              end do
            end do

          endif
        endif

        if (ly==0) cycle
        nry = nry * 2
      end do

      if (lx==0) cycle
      nrx = nrx * 2
    end do

    write (6,*) 'trunc_mat:',count (A/=0._wp)/nx,'pixels out of',nx,'kept.'

  end subroutine trunc_mat
!==============================================================================
  subroutine gp_mr_mat (A)
  real(wp) ,intent(inout) :: A (:,:)
  !--------------------------------------------------
  ! inverse matrix wavelet transform: grid -> wavelet
  !--------------------------------------------------
    call wave_1d  (A ,TR_ANA)
    A = transpose (A)
    call wave_1d  (A ,TR_ANA)
    A = transpose (A)
  end subroutine gp_mr_mat
!------------------------------------------------------------------------------
  subroutine mr_gp_mat (A)
  real(wp) ,intent(inout) :: A (:,:)
  !-------------------------------------------------
  ! direct matrix wavelet transform: wavelet -> grid
  !-------------------------------------------------
    call wave_1d  (A ,TR_SYN)
    A = transpose (A)
    call wave_1d  (A ,TR_SYN)
    A = transpose (A)
  end subroutine mr_gp_mat
!------------------------------------------------------------------------------
  subroutine mr_gp_ad_mat (A)
  real(wp) ,intent(inout) :: A (:,:)
  !--------------------------------------------------
  ! adjoint matrix wavelet transform: wavelet -> grid
  !--------------------------------------------------
    call wave_1d  (A ,TR_ADJ)
    A = transpose (A)
    call wave_1d  (A ,TR_ADJ)
    A = transpose (A)
  end subroutine mr_gp_ad_mat
!------------------------------------------------------------------------------
  function homg (x) result (A)
  real(wp) ,intent(in) :: x (:)                ! covariance function
  real(wp)             :: A (size(x),size(x))  ! covariance matrix
  !-----------------------------------------
  ! derive a 'homogeneous' covariance matrix
  ! from the covariance function
  !-----------------------------------------
    integer :: i, n
    n = size (x)
    do i = 1, n
      A (:,i) = cshift (x,1-i)
    end do
  end function homg
!==============================================================================
  subroutine factor (n, n2, nr, keep2)
  !-------------------------
  ! factorize the array size
  !-------------------------
  integer ,intent(in)            :: n     ! array size
  integer ,intent(out)           :: n2    ! number of factors of 2
  integer ,intent(out)           :: nr    ! remaining factor
  integer ,intent(in)  ,optional :: keep2 ! flag to keep a factor of 2
    integer :: i
    n2 = 0
    nr = n
    do
      if (mod(nr,2)/=0) exit
      n2 = n2 + 1
      nr = nr / 2
    end do
    if (present (keep2)) then
      if (keep2 > 0) then
        do i = 1, keep2
          n2 = n2 - 1
          nr = nr * 2
        end do
      endif
    endif
  end subroutine factor
!------------------------------------------------------------------------------
  subroutine gather (x,n)
  !------------------------------------------------------------------------
  ! reorder coefficients in the inverse wavelet transform (grid -> wavelet)
  !------------------------------------------------------------------------
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n
    integer  :: j
    real(wp) :: t (2*n)

    if (any(basis == wv_ortho)) return

    do j=1,size(x,2)
      t (  1:  n  ) = x (1:2*n:2,j)
      t (n+1:2*n  ) = x (2:2*n:2,j)
      x (  1:2*n,j) = t
    end do
  end subroutine gather
!------------------------------------------------------------------------------
  subroutine scatter (x,n)
  !-----------------------------------------------------------------------
  ! reorder coefficients in the direct wavelet transform (wavelet -> grid)
  !-----------------------------------------------------------------------
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n
    integer  :: j
    real(wp) :: t (2*n)

    if (any(basis == wv_ortho)) return

    do j=1,size(x,2)
      t             = x (  1:2*n,j)
      x (1:2*n:2,j) = t (  1:  n  )
      x (2:2*n:2,j) = t (n+1:2*n  )
    end do
  end subroutine scatter
!------------------------------------------------------------------------------
  subroutine up (x,n)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n
  !-------------------------------------------------
  ! core of the inverse wavelet transform (analysis)
  !-------------------------------------------------
    integer  :: j
    integer  :: nn
    target   :: x
    real(wp) ,pointer :: m (:)
    real(wp) ,pointer :: d (:)
    real(wp)          :: t (n)

    if (any(basis == wv_ortho)) then
      nn = 2 * n
      if (nn < 4) return
      select case (ortho_pwt_version)
      case (1)
         do j = 1, size(x,2)
            call pwt  (x(:,j), nn, isign=+1)
         end do
      case (2)
         call pwt2    (x(:,:), nn, isign=+1)
      case (3)
         call pwt_vec (x(:,:), nn, isign=+1)    ! Vectorized version for NEC
      case default
         write (0,*) "ortho_pwt_version =", ortho_pwt_version
         call finish ("up", "unsupported pwt implementation")
      end select
      return
    end if

    do j = 1, size(x,2)
      m => x(  1:  n,j)
      d => x(n+1:2*n,j)
      select case (basis)
      case (WV_LAZY)
      case (WV_LIN)
        m = m
        d = d - 0.5_wp * (m + cshift(m,1))
      case (WV_LIN_LIFTED)
        d = d - 0.5_wp  * (m + cshift(m, 1))    ! Linear predictor
        m = m + 0.25_wp * (d + cshift(d,-1))    ! Preserve average
      case (WV_CUB)
        m = m
        d = d - (0.5625_wp * (       m    + cshift(m, 1)) - &
                 0.0625_wp * (cshift(m,2) + cshift(m,-1))   )
      case (WV_CUB_LIFTED)
        d = d - (0.5625_wp * (       m    + cshift(m, 1)) - & ! Cubic predictor
                 0.0625_wp * (cshift(m,2) + cshift(m,-1))   )
        m = m + (9*(       d    + cshift(d,-1)) &
                 - (cshift(d,1) + cshift(d,-2)) ) * 0.03125_wp
      case (WV_HAAR)
        t = edsq2 * (m + d)
        d = edsq2 * (m - d)
        m = t
      case (WV_CUB_B)
        t = sq2_64_40 *         m                    &
          + sq2_64_05 * (       d    + cshift(d,-1)) &
          - sq2_64_12 * (cshift(m,1) + cshift(m,-1)) &
          + sq2_64_03 * (cshift(d,1) + cshift(d,-2))
        d = sq2_64_48 *         d                   &
          - sq2_64_32 * (       m    + cshift(m,1)) &
          + sq2_64_08 * (cshift(d,1) + cshift(d,-1))
        m = 2._wp*t
        d = 2._wp*d
      case (WV_QUA_B_INT)
        t = (-eoshift(d,-1,m(1)) + 3._wp*(m+d) - eoshift(m,1,d(n))) * 0.25_wp
        d = (-eoshift(d,-1,m(1)) + 3._wp*(m-d) + eoshift(m,1,d(n))) * 0.25_wp
        m = t
      case (WV_QUA_B)
        t = (-cshift(d,-1) + 3._wp*(m+d) - cshift(m,1)) * 0.25_wp
        d = (-cshift(d,-1) + 3._wp*(m-d) + cshift(m,1)) * 0.25_wp
        m = t
      case (WV_CUB_B_INT)
        m = 2._wp*m
        m = m - 0.5_wp*(d+eoshift(d,-1,d(1)))
        d = d - 0.5_wp*(m+eoshift(m,1,m(n)))
        m = m + 0.375_wp*(d+eoshift(d,-1,d(1)))
      case (WV_CUB_LIFT_INT)
        d = d - (0.5625_wp * (       m    + eoshift(m, 1,m(n))) - &
                 0.0625_wp * (eoshift(m,2,m(n)) + eoshift(m,-1,m(1)))   )
        m = m + (9*(       d    + eoshift(d,-1,d(1))) &
                 - (eoshift(d,1,d(n)) + eoshift(d,-2,d(1))) ) * 0.03125_wp
      case default
        write(0,*)'mo_1dmra: invalid parameter basis =',basis
        call finish('mo_1dmra','invalid parameter basis')
      end select
    end do
  end subroutine up
!------------------------------------------------------------------------------
  subroutine upad (x,n)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n
  !----------------------------------------------
  ! core of the adjoint inverse wavelet transform
  !----------------------------------------------
    integer  :: j, nn
    target   :: x
    real(wp) ,pointer :: m (:)
    real(wp) ,pointer :: d (:)
    real(wp)          :: t (n) ! temporary
    real(wp)          :: z     ! temporary

    if (any(basis == wv_ortho)) then
      nn = 2 * n
      if (nn < 4) return
      select case (ortho_pwt_version)
      case (1)
         do j = 1, size(x,2)
            call pwt  (x(:,j), nn, isign=-1)
         end do
      case (2)
         call pwt2    (x(:,:), nn, isign=-1)
      case (3)
         call pwt_vec (x(:,:), nn, isign=-1)    ! Vectorized version for NEC
      case default
         write (0,*) "ortho_pwt_version =", ortho_pwt_version
         call finish ("upad", "unsupported pwt implementation")
      end select
      return
    end if

    do j = 1, size(x,2)
      m => x(  1:  n,j)
      d => x(n+1:2*n,j)
      select case (basis)

      case (WV_LAZY)
      case (WV_LIN)
        m = m - 0.5_wp * (d + cshift(d,-1))
      case (WV_LIN_LIFTED)
        d = d + 0.25_wp * (m + cshift(m, 1))    ! adjoint of analysis
        m = m - 0.5_wp  * (d + cshift(d,-1))
      case (WV_CUB)
        m = m - (0.5625_wp * (       d     + cshift(d,-1)) - &
                 0.0625_wp * (cshift(d,-2) + cshift(d, 1))   )
      case (WV_CUB_LIFTED)
        d = d + (9*(       m    + cshift(m, 1)) &
                 - (cshift(m,2) + cshift(m,-1)) ) * 0.03125_wp
        m = m - (0.5625_wp * (       d     + cshift(d,-1)) - &
                 0.0625_wp * (cshift(d,-2) + cshift(d, 1))   )
      case (WV_HAAR)
        t = edsq2 * (m + d)
        d = edsq2 * (m - d)
        m = t
      case (WV_CUB_B)
        t = m
        m = - sq2_64_32 * (       d     + cshift(d,-1))
        d =   sq2_64_48 *         d                     &
            + sq2_64_08 * (cshift(d,-1) + cshift(d, 1))
        m = m + sq2_64_40 *         t                     &
              - sq2_64_12 * (cshift(t,-1) + cshift(t, 1))
        d = d + sq2_64_05 * (       t     + cshift(t, 1)) &
              + sq2_64_03 * (cshift(t,-1) + cshift(t, 2))
        m = 2._wp*m
        d = 2._wp*d
      case (WV_QUA_B_INT)
        t = m
        m =     0.25_wp * eoshift(d,-1) + 0.75_wp * d; m(1) = m(1) - d(1)/4._wp
        z = d(n)
        d =   - 0.25_wp * eoshift(d, 1) - 0.75_wp * d; d(n) = d(n) + z/4._wp
        m = m - 0.25_wp * eoshift(t,-1) + 0.75_wp * t; m(1) = m(1) - t(1)/4._wp
        d = d - 0.25_wp * eoshift(t, 1) + 0.75_wp * t; d(n) = d(n) - t(n)/4._wp
      case (WV_QUA_B)
        t = m
        m =     0.25_wp * cshift(d,-1) + 0.75_wp * d
        d =   - 0.25_wp * cshift(d, 1) - 0.75_wp * d
        m = m - 0.25_wp * cshift(t,-1) + 0.75_wp * t
        d = d - 0.25_wp * cshift(t, 1) + 0.75_wp * t
      case (WV_CUB_LIFT_INT)
        d = d + (9*(       m    + eoshift(m,1)) &
                 - (eoshift(m,2) + eoshift(m,-1)) ) * 0.03125_wp
        d(1) = d(1) + 0.25_wp*m(1)
        if(n>1) d(1) = d(1) - 0.03125_wp*m(2)
        d(n) = d(n) - 0.03125_wp*m(n)
        m = m - (0.5625_wp * (       d     + eoshift(d,-1)) - &
                 0.0625_wp * (eoshift(d,-2) + eoshift(d,1))   )
        m(1) = m(1) + 0.0625_wp*d(1)
        m(n) = m(n) - 0.5_wp*d(n)
        if(n>1) m(n) = m(n) + 0.0625_wp*d(n-1)
      case (WV_CUB_B_INT)
        d = d + 3._wp*(m+eoshift(m,1))/8._wp; d(1) = d(1) + 3._wp*m(1)/8._wp
        m = m - 0.5_wp*(d+eoshift(d,-1)); m(n) = m(n) - 0.5_wp* d(n)
        d = d - 0.5_wp*(m+eoshift(m,1)); d(1) = d(1) - 0.5_wp * m(1)
        m = 2._wp*m
      case default
        call finish('upad','invalid parameter basis')
      end select
    end do
  end subroutine upad
!------------------------------------------------------------------------------
  subroutine down (x,n)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n
  !------------------------------
  ! core direct wavelet transform
  !------------------------------
    integer  :: j, nn
    target   :: x
    real(wp) ,pointer :: m (:)
    real(wp) ,pointer :: d (:)
    real(wp)          :: t (n)

    if (any(basis == wv_ortho)) then
      nn = 2 * n
      if (nn < 4) return
      select case (ortho_pwt_version)
      case (1)
         do j = 1, size(x,2)
            call pwt  (x(:,j), nn, isign=-1)
         end do
      case (2)
         call pwt2    (x(:,:), nn, isign=-1)
      case (3)
         call pwt_vec (x(:,:), nn, isign=-1)    ! Vectorized version for NEC
      case default
         write (0,*) "ortho_pwt_version =", ortho_pwt_version
         call finish ("down", "unsupported pwt implementation")
      end select
      return
    end if

    do j = 1, size(x,2)
      m => x(  1:  n,j)
      d => x(n+1:2*n,j)
      select case (basis)
      case (WV_LAZY)
      case (WV_LIN)
        d = 0.5_wp * (m + cshift(m,1)) + d
        m = m
      case (WV_LIN_LIFTED)
        m = m - 0.25_wp * (d + cshift(d,-1)) ! synthesis
        d = d + 0.5_wp  * (m + cshift(m, 1))
      case (WV_CUB)
        d = 0.5625_wp * (       m    + cshift(m, 1)) - &
            0.0625_wp * (cshift(m,2) + cshift(m,-1)) + d
        m = m
      case (WV_CUB_LIFTED)
        m = m - (9*(       d    + cshift(d,-1)) &
                 - (cshift(d,1) + cshift(d,-2)) ) * 0.03125_wp
        d = d + (0.5625_wp * (       m    + cshift(m, 1)) - &
                 0.0625_wp * (cshift(m,2) + cshift(m,-1))   )
      case (WV_HAAR)
        t = edsq2 * (m + d)
        d = edsq2 * (m - d)
        m = t
      case (WV_CUB_B)
        t = sq2_64_48 *        m &
          + sq2_64_08 *(cshift(m,-1) + cshift(m, 1)) &
          - sq2_64_05 *(       d     + cshift(d,-1)) &
          - sq2_64_03 *(cshift(d,-2) + cshift(d, 1))
        d = sq2_64_32 * (       m     + cshift(m,1)) &
          + sq2_64_40 *         d                    &
          - sq2_64_12 * (cshift(d,-1) + cshift(d,1))
        m = t
      case (WV_QUA_B_INT)
        t = (eoshift(m,-1,m(1)) - eoshift(d,-1,-d(1)) + 3._wp*(m+d)) * 0.25_wp
        d = (eoshift(m, 1,m(n)) + eoshift(d, 1,-d(n)) + 3._wp*(m-d)) * 0.25_wp
        m = t
      case (WV_QUA_B)
        t = (cshift(m,-1) - cshift(d,-1) + 3._wp*(m+d)) * 0.25_wp
        d = (cshift(m, 1) + cshift(d, 1) + 3._wp*(m-d)) * 0.25_wp
        m = t
      case (WV_CUB_B_INT)
        m = m - 0.375_wp*(d+eoshift(d,-1,d(1)))
        d = d + 0.5_wp*(m+eoshift(m,1,m(n)))
        m = m + 0.5_wp*(d+eoshift(d,-1,d(1)))
        m = 0.5_wp*m
      case (WV_CUB_LIFT_INT)
        m = m - (9*(       d    + eoshift(d,-1,d(1))) &
                 - (eoshift(d,1,d(n)) + eoshift(d,-2,d(1))) ) * 0.03125_wp
        d = d + (0.5625_wp * (       m    + eoshift(m, 1,m(n))) - &
                 0.0625_wp * (eoshift(m,2,m(n)) + eoshift(m,-1,m(1)))   )
      case default
        call finish('mo_1dmra','invalid parameter basis')
      end select

    end do
  end subroutine down
!------------------------------------------------------------------------------
  subroutine downad (x,n)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n
  !-------------------------------
  ! core adjoint wavelet transform
  !-------------------------------
    integer  :: j, nn
    target   :: x
    real(wp) ,pointer :: m (:)
    real(wp) ,pointer :: d (:)
    real(wp)          :: t (n) ! temporary
    real(wp)          :: z     ! temporary

    if (any(basis == wv_ortho)) then
      nn = 2 * n
      if (nn < 4) return
      select case (ortho_pwt_version)
      case (1)
         do j = 1, size(x,2)
            call pwt  (x(:,j), nn, isign=+1)
         end do
      case (2)
         call pwt2    (x(:,:), nn, isign=+1)
      case (3)
         call pwt_vec (x(:,:), nn, isign=+1)    ! Vectorized version for NEC
      case default
         write (0,*) "ortho_pwt_version =", ortho_pwt_version
         call finish ("downad", "unsupported pwt implementation")
      end select
      return
    end if

    do j = 1, size(x,2)
      m => x(  1:  n,j)
      d => x(n+1:2*n,j)

      select case (basis)
      case (WV_LAZY)
      case (WV_LIN)
        m = m + 0.5_wp * (d + cshift(d,-1))
      case (WV_LIN_LIFTED)
        m = m + 0.5_wp  * (d + cshift(d,-1)) ! adjoint of synthesis
        d = d - 0.25_wp * (m + cshift(m, 1))
      case (WV_CUB)
        m = m + 0.5625_wp * (       d    + cshift(d,-1)) &
              - 0.0625_wp * (cshift(d,1) + cshift(d,-2))
      case (WV_CUB_LIFTED)
        m = m + (0.5625_wp * (       d     + cshift(d,-1)) - &
                 0.0625_wp * (cshift(d,-2) + cshift(d, 1))   )
        d = d - (9*(       m    + cshift(m, 1)) &
                 - (cshift(m,2) + cshift(m,-1)) ) * 0.03125_wp
      case (WV_HAAR)
        t = edsq2 * (m + d)
        d = edsq2 * (m - d)
        m = t
      case (WV_CUB_B)
        t = m
        m = sq2_64_32 * (       d     + cshift(d,-1))
        d = sq2_64_40 *         d                    &
          - sq2_64_12 * (cshift(d, 1) + cshift(d,-1))
        m = m + sq2_64_48 *        t                     &
              + sq2_64_08 *(cshift(t, 1) + cshift(t,-1))
        d = d - sq2_64_05 *(       t     + cshift(t, 1)) &
              - sq2_64_03 *(cshift(t, 2) + cshift(t,-1))
      case (WV_QUA_B_INT)
        t = m
        m =     0.25_wp * eoshift(d,-1) + 0.75_wp * d; m(n) = m(n) + d(n) * 0.25_wp
        z = d(n)
        d =     0.25_wp * eoshift(d,-1) - 0.75_wp * d; d(n) = d(n) - z    * 0.25_wp
        m = m + 0.25_wp * eoshift(t, 1) + 0.75_wp * t; m(1) = m(1) + t(1) * 0.25_wp
        d = d - 0.25_wp * eoshift(t, 1) + 0.75_wp * t; d(1) = d(1) + t(1) * 0.25_wp
      case (WV_QUA_B)
        t = m
        m =     0.25_wp * cshift(d,-1) + 0.75_wp * d
        d =     0.25_wp * cshift(d,-1) - 0.75_wp * d
        m = m + 0.25_wp * cshift(t, 1) + 0.75_wp * t
        d = d - 0.25_wp * cshift(t, 1) + 0.75_wp * t
      case (WV_CUB_LIFT_INT)
        m = m + (0.5625_wp * (       d     + eoshift(d,-1)) - &
                 0.0625_wp * (eoshift(d,-2) + eoshift(d,1))   )
        m(1) = m(1) - 0.0625_wp*d(1)
        m(n) = m(n) + 0.5_wp*d(n)
        if(n>1) m(n) = m(n) - 0.0625_wp*d(n-1)
        d = d - (9*(       m    + eoshift(m,1)) &
                 - (eoshift(m,2) + eoshift(m,-1)) ) * 0.03125_wp
        if(n>1) d(1) = d(1) + 0.03125_wp*m(2)
        d(1) = d(1) - 0.25_wp*m(1)
        d(n) = d(n) + 0.03125_wp*m(n)
      case (WV_CUB_B_INT)
        m = m/2._wp
        d = d + (m+eoshift(m,1))/2._wp; d(1) = d(1) + 0.5_wp*m(1)
        m = m + (d+eoshift(d,-1))/2._wp; m(n) = m(n) + 0.5_wp*d(n)
        d = d - 3._wp*(m+eoshift(m,1))/8._wp; d(1) = d(1) - 0.375_wp*m(1)
      case default
        call finish('mo_1dmra','invalid parameter basis')
      end select

    end do
  end subroutine downad
!------------------------------------------------------------------------------
  subroutine normalize (A,x)
  !---------------------------------------------------
  ! scale a matrix so that the diagonal is 1
  ! A -> z A z^t  with  z = 1/sqrt(abs(diag(A)))
  ! optionally return sqrt(abs(diag) of the original A
  !---------------------------------------------------
  real(wp) ,intent(inout)           :: A (:,:) ! matrix to scale
  real(wp) ,intent(out) ,optional   :: x (:)   ! sqrt(abs(diag(A)))

    real(wp) :: z(size(A,1))
    integer  :: n, i

    n = size(A,1)
    if (size(A,2)/= n) call finish('normalize','A is not square')

    !-------------------------
    ! determine scaling factor
    !-------------------------
    do i=1,n
      if (A(i,i)/=0._wp) then
        z(i) = sqrt(abs(A(i,i)))
      else
        z(i) = 1._wp
      endif
    end do
    if (present(x)) x = z

    !------
    ! scale
    !------
    z = 1._wp / z
    call scale (A,z)

  end subroutine normalize
!------------------------------------------------------------------------------
  subroutine scale (A,x)
  !---------------------------------------------------
  ! scale a matrix
  ! A -> x A x^t  with  z = 1/sqrt(abs(diag(A)))
  !---------------------------------------------------
  real(wp) ,intent(inout) :: A (:,:) ! matrix
  real(wp) ,intent(in)    :: x (:)   ! scaling factors

    integer  :: n, i

    n = size(A,1)
    if (size(A,2)/= n) &
      call finish('scale','A is not square')

    do i=1,n
      A(:,i) = A(:,i) * x
      A(i,:) = A(i,:) * x
    end do
  end subroutine scale
!==============================================================================
  subroutine wave_1d_odd (x, trans)
  !---------------------------------------------!
  ! applies wavelet transform to arrays with    !
  ! length unequal to a power of 2              !
  ! available only for Quadratic B-spline basis !
  !---------------------------------------------!
  real(wp) ,intent(inout)        :: x(:,:) ! array to transform
  integer  ,intent(in)           :: trans  ! transformation
  integer  :: n2,nx

    nx=size(x,1)
    n2=0
    do
       if (nx<=2**n2) exit
       n2=n2+1
    end do

    basis = 9

    select case (trans)
    case (TR_SYN)
      call mr_gp_odd    (x,n2)
    case (TR_ANA)
      call gp_mr_odd    (x,n2)
    case (TR_DUAL)
      call gp_mr_ad_odd (x, n2)
    case (TR_ADJ)
      call mr_gp_ad_odd (x, n2)
    case default
      write (0,*)  'wave_1d_odd:  invalid transform =',trans
      call finish ('wave_1d_odd','invalid transform')
    end select

  end subroutine wave_1d_odd
!------------------------------------------------------------------------------
  subroutine mr_gp_odd (x,n2)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in) ::n2

    integer           :: nx,l
    integer           :: f(n2+1)

    nx = size(x,1)
    call factor_odd(nx,f)

    do l=1,n2-1
       if (f(n2-l)-f(n2-l+1)/=f(n2-l+1)) then
          call down_odd (x,f(n2-l+1))
          call scatter_odd (x,f(n2-l+1))
       else
          call down    (x,f(n2-l+1))
          call scatter (x,f(n2-l+1))
       end if
    end do

  end subroutine mr_gp_odd
!------------------------------------------------------------------------------
  subroutine gp_mr_odd (x,n2)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in) :: n2

    integer           :: nx,l
    integer           :: f(n2+1)

    nx = size(x,1)
    call factor_odd(nx,f)

    do l=1,n2-1
      if (f(l)-f(l+1)/=f(l+1)) then
         call gather_odd (x,f(l+1))
         call up_odd     (x,f(l+1))
      else
         call gather (x,f(l+1))
         call up     (x,f(l+1))
      end if
    end do

  end subroutine gp_mr_odd
!------------------------------------------------------------------------------
  subroutine factor_odd(nx,f)

    integer, intent(in) :: nx
    integer, intent(out) :: f(:)

    integer :: n2,nr,i,h1,h2

    n2=0
    do
       if (nx<=2*2**n2) exit
       n2=n2+1
    end do
    nr=nx-2**n2

    i=1
    if (nr==0) then
       do i=1,n2+1
          f(i)=nx/2**(i-1)
       end do
    else
       h1=nx
       f(1)=nx
       do i=1,n2+1
          h2=h1/2
          if (2*h2==h1) then
             f(i+1)=h2
             h1=h2
          else
             f(i+1)=h2+1
             h1=h2+1
          end if
       end do
    end if

  end subroutine factor_odd
!------------------------------------------------------------------------------
  subroutine gather_odd (x,n)

  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n
    integer  :: j
    real(wp) :: t (2*n-1)

    do j=1,size(x,2)
      t (  1:  n  ) = x (1:2*n-1:2,j)
      t (n+1:2*n-1) = x (2:2*n-1:2,j)
      x (1:2*n-1,j) = t
    end do

  end subroutine gather_odd
!------------------------------------------------------------------------------
  subroutine scatter_odd (x,n)

  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n
    integer  :: j
    real(wp) :: t (2*n-1)

    do j=1,size(x,2)
      t               = x (  1:2*n-1,j)
      x (1:2*n-1:2,j) = t (  1:  n  )
      x (2:2*n-1:2,j) = t (n+1:2*n-1)
    end do

  end subroutine scatter_odd
!------------------------------------------------------------------------------
  subroutine up_odd (x,n)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n

    integer  :: j
    target   :: x
    real(wp) ,pointer :: m (:)
    real(wp) ,pointer :: d (:)
    real(wp)          :: t (n)

    do j = 1, size(x,2)
      m => x(  1:  n,j)
      d => x(n+1:2*n-1,j)

!     QUADRATIC B-SPLINE ON INTERVAL

      t(1:n-1) = d
      t(n) = 0._wp
      d = (-eoshift(d,-1,m(1)) + 3._wp*(m(1:n-1)-d) + m(2:n))/4._wp
      m = (-eoshift(t,-1,m(1)) + 3._wp*(m+t) - eoshift(m,1,-2._wp*m(n)))/4._wp
    end do

  end subroutine up_odd
!------------------------------------------------------------------------------
  subroutine down_odd (x,n)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n

    integer  :: j
    target   :: x
    real(wp) ,pointer :: m (:)
    real(wp) ,pointer :: d (:)
    real(wp)          :: t (n),s(n)

    do j = 1, size(x,2)
      m => x(  1:  n,j)
      d => x(n+1:2*n-1,j)

!     QUADRATIC B-SPLINE ON INTERVAL

      t(1:n-1) = d
      t(n) = (-m(n-1)+d(n-1)+m(n))/7._wp
      s = (eoshift(m,-1,m(1)) - eoshift(t,-1,-t(1)) + 3._wp*(m+t))/4._wp
      t = (eoshift(m,1) + eoshift(t,1) + 3._wp*(m-t))/4._wp
      m = s
      d = t(1:n-1)

    end do

  end subroutine down_odd
!------------------------------------------------------------------------------
  subroutine mr_gp_ad_odd (x, n2)

  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in) :: n2

    integer           :: nx,l
    integer           :: f(n2+1)

    nx = size(x,1)
    call factor_odd(nx,f)

    do l=1,n2-1
      if (f(l)-f(l+1)/=f(l+1)) then
         call gather_odd (x,f(l+1))
         call downad_odd     (x,f(l+1))
      else
         call gather (x,f(l+1))
         call downad     (x,f(l+1))
      end if
    end do

  end subroutine mr_gp_ad_odd
!------------------------------------------------------------------------------
  subroutine gp_mr_ad_odd (x, n2)

  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in) ::n2

    integer           :: nx,l
    integer           :: f(n2+1)

    nx = size(x,1)
    call factor_odd(nx,f)

    do l=1,n2-1
       if (f(n2-l)-f(n2-l+1)/=f(n2-l+1)) then
          call upad_odd (x,f(n2-l+1))
          call scatter_odd (x,f(n2-l+1))
       else
          call upad    (x,f(n2-l+1))
          call scatter (x,f(n2-l+1))
       end if
    end do

  end subroutine gp_mr_ad_odd
!==============================================================================
  subroutine upad_odd (x,n)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n

    integer  :: j
    target   :: x
    real(wp) ,pointer :: m (:)
    real(wp) ,pointer :: d (:)
    real(wp)          :: t (n), s(n)

    do j = 1, size(x,2)
      m => x(  1:  n,j)
      d => x(n+1:2*n-1,j)

!     QUADRATIC B-SPLINE ON INTERVAL

      t(1:n-1) = d
      t(n) = 0._wp
      s = (-eoshift(m,-1,m(1)) + eoshift(t,-1,-t(1)) + 3._wp*(m+t))/4._wp
      s(n) = s(n) + 0.5_wp*m(n)
      d = (3._wp*m(1:n-1) - 3._wp*d - m(2:n) - eoshift(d,1))/4._wp
      m = s
    end do

  end subroutine upad_odd
!------------------------------------------------------------------------------
  subroutine downad_odd (x,n)
  real(wp) ,intent(inout) :: x(:,:)
  integer  ,intent(in)    :: n

    integer  :: j
    target   :: x
    real(wp) ,pointer :: m (:)
    real(wp) ,pointer :: d (:)
    real(wp)          :: t (n),s(n)

    do j = 1, size(x,2)
      m => x(  1:  n,j)
      d => x(n+1:2*n-1,j)

!     QUADRATIC B-SPLINE ON INTERVAL

      t(1:n-1) = d
      t(n) = (d(n-1)+3._wp*m(n))/21._wp
      s = (eoshift(t,-1,m(1)) + 3._wp*(m+t) + eoshift(m,1))/4._wp
      s(n-1) = s(n-1) - 3._wp*t(n)/4._wp
      d = (eoshift(d,-1,m(1)) + 3._wp*(m(1:n-1)-d) - m(2:n))/4._wp
      d(n-1) = d(n-1) + 3._wp*t(n)/4._wp
      m = s

   end do

  end subroutine downad_odd
!==============================================================================
! orthogonal wavelets (Daubechies + co)
! adopted from the Numerical Recipes
!------------------------------------------------------------------------------
  subroutine set_coefs
  !-----------------------------------------------------------------------
  ! Initializing routine for PWT, implementing wavelet filters as selected
  ! by the module variable BASIS. This routine is called from PWT if
  ! BASIS has changed.
  !-----------------------------------------------------------------------
    integer  :: k
    real(wp) :: s
    select case (basis)
    case (WV_DAUB_4)
      ncof = 4
      cc (1:ncof) = daub4
    case (WV_DAUB_6)
      ncof = 6
      cc (1:ncof) = daub6
    case (WV_DAUB_8)
      ncof = 8
      cc (1:ncof) = daub8
    case (WV_DAUB_10)
      ncof = 10
      cc (1:ncof) = daub10
    case (WV_DAUB_12)
      ncof = 12
      cc (1:ncof) = daub12
    case (WV_DAUB_14)
      ncof = 14
      cc (1:ncof) = daub14
    case (WV_DAUB_16)
      ncof = 16
      cc (1:ncof) = daub16
    case (WV_DAUB_18)
      ncof = 18
      cc (1:ncof) = daub18
    case (WV_DAUB_20)
      ncof = 20
      cc (1:ncof) = daub20
    case (WV_SYM_8)
      ncof = 8
      cc (1:ncof) = sym8 / sqrt(2._wp)
    case (WV_SYM_10)
      ncof = 10
      cc (1:ncof) = sym10 / sqrt(2._wp)
    case (WV_SYM_12)
      ncof = 12
      cc (1:ncof) = sym12 / sqrt(2._wp)
    case (WV_SYM_18)
      ncof = 18
      cc (1:ncof) = sym18 / sqrt(2._wp)
    case default
      call finish ('set_coefs','invalid base: '//wv_name(basis))
    end select
!   ioff = - ncof / 2 ! center the support of the wavelet
!   joff = - ncof / 2
    ioff = -        2 ! center the peak of the wavelet
    joff = - ncof + 2
    s = -1._wp
    do k=1,ncof
      cr (ncof+1-k) = s * cc (k)
      s = -s
    end do
    basis_ortho = basis
  end subroutine set_coefs
!------------------------------------------------------------------------------
  subroutine pwt (a,n,isign)
  integer ,intent(in)    :: n
  real(wp),intent(inout) :: a(n)
  integer ,intent(in)    :: isign
  !-----------------------------------------------------------------------
  ! Partial wavelet transform: Applies an arbitrary wavelet filter to data
  ! vector A(1:n) (for isign=1) or applies its transpose (for
  ! isign=-1). Used hierarchically by routines UP, DOWN and their
  ! adjoints. The filter is determined by a call to set_coefs, which
  ! initializes the module variables NCOF, CC, CR, IOFF, JOFF.
  !-----------------------------------------------------------------------
    integer  ::  i,ii,jf,jr,k,ni,nj,nh,nmod
    real(wp) :: wksp(n)
    real(wp) ::  ai,ai1

    if(basis /= basis_ortho) call set_coefs

    if (n < 4) return
    nmod = ncof * n
    nh   = n / 2
    wksp = 0._wp
    if (isign >= 0) then
      ii=1
      do i=1,n,2
        ni=i+nmod+ioff
        nj=i+nmod+joff
!NEC$ loop_count(20)
        do k=1,ncof
          jf=mod (ni+k,n)
          jr=mod (nj+k,n)
          wksp(ii)=wksp(ii)+cc(k)*a(jf+1)
          wksp(ii+nh)=wksp(ii+nh)+cr(k)*a(jr+1)
        end do
        ii=ii+1
      end do
    else
      ii=1
      do i=1,n,2
        ai=a(ii)
        ai1=a(ii+nh)
        ni=i+nmod+ioff
        nj=i+nmod+joff
!NEC$ loop_count(20)
        do k=1,ncof
          jf=mod (ni+k,n) +1
          jr=mod (nj+k,n) +1
          wksp(jf)=wksp(jf)+cc(k)*ai
          wksp(jr)=wksp(jr)+cr(k)*ai1
        end do
        ii=ii+1
      end do
    endif
    a = wksp
  end subroutine pwt
!------------------------------------------------------------------------------
  subroutine pwt_vec (a, n, isign)
    real(wp),intent(inout) :: a(:,:)
    integer ,intent(in)    :: n
    integer ,intent(in)    :: isign
    !-----------------------------------------------------------------------
    ! Vectorized version of partial wavelet transform from Numerical Recipes
    !-----------------------------------------------------------------------
    integer  :: i, ii, jf, jr, k, ni, nj, nh, nmod
    integer  :: nsz
#if defined (__NEC__)
    !-----------------------------------------------------------------------
    ! Extended first dimension of work arrays by +1 to reduce bank conflicts
    !-----------------------------------------------------------------------
    real(wp) :: a_  (n+1,size (a,2))                ! Auxiliary array
    real(wp) :: a_t (size (a,2)+1,n)                ! Transposed array
    real(wp) :: work(size (a,2)+1,n)
#else
    real(wp) :: a_t (size (a,2)  ,n)                ! Transposed array
    real(wp) :: work(size (a,2)  ,n)
#endif

    if (basis /= basis_ortho) call set_coefs

    if (n < 4) return
    nmod = ncof * n
    nh   = n / 2
    nsz  = size (a,dim=2)
#if defined (__NEC__)
    a_ (1:n,  :) = a(1:n,:)
    a_t(1:nsz,:) = transpose (a_(1:n,:))
#else
    a_t(1:nsz,:) = transpose (a(1:n,:))
#endif
    work = 0
    if (isign >= 0) then
       do ii = 1, nh
          i = 2*ii-1
          ni=i+nmod+ioff
!NEC loop_count(20)
          do k=1,ncof
             jf=mod (ni+k,n) +1
             work(1:nsz,ii   )=work(1:nsz,ii   )+cc(k)*a_t(1:nsz,jf)
          end do
       end do
       do ii = 1, nh
          i = 2*ii-1
          nj=i+nmod+joff
!NEC loop_count(20)
          do k=1,ncof
             jr=mod (nj+k,n) +1
             work(1:nsz,ii+nh)=work(1:nsz,ii+nh)+cr(k)*a_t(1:nsz,jr)
          end do
       end do
    else
       do ii = 1, nh
          i = 2*ii-1
          ni=i+nmod+ioff
          nj=i+nmod+joff
!NEC loop_count(20)
          do k=1,ncof
             jf=mod (ni+k,n) +1
             work(1:nsz,jf)=work(1:nsz,jf)+cc(k)*a_t(1:nsz,ii   )
          end do
!NEC loop_count(20)
          do k=1,ncof
             jr=mod (nj+k,n) +1
             work(1:nsz,jr)=work(1:nsz,jr)+cr(k)*a_t(1:nsz,ii+nh)
          end do
       end do
    endif
    a(1:n,:) = transpose (work(1:nsz,1:n))
  end subroutine pwt_vec
!------------------------------------------------------------------------------
  subroutine pwt2 (a, n, isign)
  real(wp),intent(inout) :: a(:,:)
  integer ,intent(in)    :: n
  integer ,intent(in)    :: isign
  !-----------------------------------------------------------------------
  ! Partial wavelet transform: Applies an arbitrary wavelet filter to the
  ! data array slice A(1:n,:) (for isign=1) or applies its transpose (for
  ! isign=-1). Used hierarchically by routines UP, DOWN and their
  ! adjoints. The filter is determined by a call to set_coefs, which
  ! initializes the module variables NCOF, CC, CR, IOFF, JOFF.
  !
  ! This implementation is compatible with 'pwt' from Numerical Recipes
  ! but empirically much faster when processing multiple long vectors.
  !-----------------------------------------------------------------------
    integer  :: i, j, k, n1, n2, nh
    ! Automatic arrays appear to be faster than allocatables...
    ! (Important note: bounds are tuned to the choices in set_coefs)
    real(wp) :: work(5-mcof:n+(mcof-2))
    !real(wp), allocatable :: work(:)

    if (n < 4) return

    if(basis /= basis_ortho) call set_coefs
    !-----------------------------------------------
    ! Allocate work array for unfolded data + filter
    !-----------------------------------------------
    n1 = min (3        + min (ioff, joff), 1)
    n2 = max (n + ncof + max (ioff, joff), n)
    !allocate (work(n1:n2))
    nh = n / 2

!$omp parallel private(i,j,k,work)
    if (isign >= 0) then
!$omp do schedule(static)
       do j = 1, size (a,2)
          !---------------------------------------------
          ! Analysis: unfold data column into work space
          ! and extend by coefficients from wrap-around
          !---------------------------------------------
          do k = n1, 0
             work(k) = a(modulo (k-1,n)+1,j)
          end do
          work(1:n) = a(1:n,j)
          do k = n+1, n2
             work(k) = a(modulo (k-1,n)+1,j)
          end do
          !----------------------------------------------------------
          ! Application of high- & low-pass filters with downsampling
          !----------------------------------------------------------
          do i = 1, nh
             a(i,j)    = dot_product (cc(1:ncof), &
                                      work(2*i+ioff+1:2*i+ioff+ncof))
             a(i+nh,j) = dot_product (cr(1:ncof), &
                                      work(2*i+joff+1:2*i+joff+ncof))
          end do
       end do
!$omp end do
    else ! isign < 0
!$omp do schedule(static)
       do j = 1, size (a,2)
          !--------------------------------------------------------------------
          ! Synthesis: scatter reconstruction of data column to extended vector
          !--------------------------------------------------------------------
          work(n1:n2) = 0
          do i = 1, nh
             !------------------------------------------------------
             ! This loop may be slow with some compilers... (*sigh*)
             !------------------------------------------------------
             work(2*i+ioff+1:2*i+ioff+ncof) = &
                  work(2*i+ioff+1:2*i+ioff+ncof) + a(i,   j) * cc(1:ncof)
             work(2*i+joff+1:2*i+joff+ncof) = &
                  work(2*i+joff+1:2*i+joff+ncof) + a(i+nh,j) * cr(1:ncof)
          end do
          !---------------------------------------------------
          ! Gather: fold extended vector back to original size
          !---------------------------------------------------
          do k = n1, 0
             i = modulo (k-1,n)+1
             work(i) = work(i) + work(k)
          end do
          do k = n+1, n2
             i = modulo (k-1,n)+1
             work(i) = work(i) + work(k)
          end do
          a(1:n,j) = work(1:n)
       end do
!$omp end do
    endif
!$omp end parallel
  end subroutine pwt2
!==============================================================================
  !======================================================================
  ! Non-periodic wavelet transform on an interval using
  ! the approach by Cohen, Daubechies, Vial, and Jawerth.
  ! We use the algorithm described in ref.[1], with the
  ! improved and corrected numerical coefficients from [2].
  !
  ! The original version of the transform does not provide
  ! filter coefficients for a continuation of the pyramidal
  ! algorithm to dyadic grids smaller than 4*N, where N is
  ! the order of the wavelet.
  ! The implementation of the algorithm below extends the
  ! transform with auxiliary ("joining") filters.  Their
  ! coefficients are determined by requiring that at any
  ! level J (where 2^J <= 2*N) the 2^J coarse scale
  ! coefficients generate all polynomials up to degree 2^J-1.
  !
  ! Note that except for the non-orthogonal preconditioning
  ! step affecting N grid points at the left and right edges
  ! the proper CD(J)V wavelet transform is orthonormal.
  !
  ! [1] A. Cohen, I. Daubechies, P. Vial,
  !     "Wavelets on the Interval and Fast Wavelet Transforms",
  !     Appl. Comp. Harm. Anal. 1 (1993) 54-81.
  ! [2] http://www.netlib.org/a/wavelets94
  !======================================================================
  subroutine wave_cdv (x, trans, basis)
    real(wp) ,intent(inout)        :: x(:,:)    ! Array to transform
    integer  ,intent(in)           :: trans     ! Transformation
    integer  ,intent(in), optional :: basis     ! Wavelet basis

    if (present (basis)) then
       if (basis /= cdv_basis) call cdv_setup (basis)
    end if
    if (cdv%N == 0) call cdv_setup ()

    select case (cdv_version)
    case (1)         ! Optimized for scalar processors
       select case (trans)
       case (TR_SYN)
          call cdv_syn    (x)
       case (TR_ANA)
          call cdv_ana    (x)
       case (TR_DUAL)
          call cdv_ana_ad (x)
       case (TR_ADJ)
          call cdv_syn_ad (x)
       case default
          write (0,*)  'wave_cdv:  invalid transform =',trans
          call finish ('wave_cdv','invalid transform')
       end select
    case (2)         ! Optimized for vector processors
       select case (trans)
       case (TR_SYN)
          call cdv_syn_vec    (x)
       case (TR_ANA)
          call cdv_ana_vec    (x)
       case (TR_DUAL)
          call cdv_ana_ad_vec (x)
       case (TR_ADJ)
          call cdv_syn_ad_vec (x)
       case default
          write (0,*)  'wave_cdv:  invalid transform =',trans
          call finish ('wave_cdv','invalid transform')
       end select
    end select
  end subroutine wave_cdv

  ! --

  subroutine cdv_setup (basis)
    integer, intent(in), optional :: basis

    if (present (basis)) then
       cdv_basis = basis
    end if

    select case (cdv_basis)
    case (WV_CDV_4)
       if (cdv%N /= 2) call cdjv4_setup ()
    case (WV_CDV_8)
       if (cdv%N /= 4) call cdjv8_setup ()
    case default
       print *, "cdv_setup: invalid interval wavelet basis:", cdv_basis
       call finish ("cdv_setup","invalid basis")
    end select
  end subroutine cdv_setup

  ! --

  subroutine cdv_clear ()
!    if (associated (cdv%c_fwd)) then
    if (allocated (cdv_c_fwd)) then
       deallocate (cdv_c_fwd, cdv_c_inv,                &
            cdv_le_hi, cdv_le_lo, cdv_re_hi, cdv_re_lo, &
            cdv_l_pre, cdv_l_post, cdv_r_pre, cdv_r_post)
    end if
    cdv%N   = 0
    cdv%lf  = 0
    cdv%lef = 0
  end subroutine cdv_clear

  ! --

  subroutine cdv_init (n)
    integer, intent(in) :: n

    ! Local variables
    integer :: lf, lef

    lf      = 2*n       ! Interior filter length
    lef     = 3*n - 1   ! Maximum edge filter length
    cdv%N   = n
    cdv%lf  = lf
    cdv%lef = lef

    allocate (cdv_c_fwd(lf), cdv_c_inv(lf), &
         cdv_le_hi(n,lef), cdv_le_lo(n,lef), &
         cdv_re_hi(n,lef), cdv_re_lo(n,lef), &
         cdv_l_pre(n,n), cdv_l_post(n,n), cdv_r_pre(n,n), cdv_r_post(n,n))
  end subroutine cdv_init

  ! --

  subroutine cdv_apply_pre (v)
    real(wp), intent(inout) :: v(:,:)
    ! Local variables
    integer :: M, N
    M = size (v, dim=1)
    N = cdv%N
    ! Apply preconditioning matrices
    v(1:N,    :) = matmul (cdv_l_pre(1:N,1:N), v(1:N,    :))
    v(M-N+1:M,:) = matmul (cdv_r_pre(1:N,1:N), v(M-N+1:M,:))
  end subroutine cdv_apply_pre

  ! --

  subroutine cdv_apply_pre_ad (v)
    real(wp), intent(inout) :: v(:,:)
    ! Local variables
    integer :: M, N
    M = size (v, dim=1)
    N = cdv%N
    ! Apply adjoint preconditioning matrices
    v(1:N,    :) = matmul (transpose (cdv_l_pre(1:N,1:N)), v(1:N,    :))
    v(M-N+1:M,:) = matmul (transpose (cdv_r_pre(1:N,1:N)), v(M-N+1:M,:))
  end subroutine cdv_apply_pre_ad

  ! --

  subroutine cdv_apply_post (v)
    real(wp), intent(inout) :: v(:,:)
    ! Local variables
    integer :: M, N
    M = size (v, dim=1)
    N = cdv%N
    ! Apply postconditioning matrices
    v(1:N,    :) = matmul (cdv_l_post(1:N,1:N), v(1:N,    :))
    v(M-N+1:M,:) = matmul (cdv_r_post(1:N,1:N), v(M-N+1:M,:))
  end subroutine cdv_apply_post

  ! --

  subroutine cdv_apply_post_ad (v)
    real(wp), intent(inout) :: v(:,:)
    ! Local variables
    integer :: M, N
    M = size (v, dim=1)
    N = cdv%N
    ! Apply adjoint postconditioning matrices
    v(1:N,    :) = matmul (transpose (cdv_l_post(1:N,1:N)), v(1:N,    :))
    v(M-N+1:M,:) = matmul (transpose (cdv_r_post(1:N,1:N)), v(M-N+1:M,:))
  end subroutine cdv_apply_post_ad

  ! --

  subroutine cdv_dyad_down (v, s, d)
    ! Perform one iteration in the cascade algorithm (analysis step)
    real(wp), intent(in)  :: v(:,:)
    real(wp), intent(out), dimension (size (v,dim=1)/2,size (v,dim=2)) :: s, d

    ! Local variables
    integer :: M, M2, N, lf, lef, k, j, P
#ifdef __NEC__
    integer :: i
#endif

    M   = size (v, dim=1)
    M2  = M/2
    P   = size (s, dim=2)
    N   = cdv%N
    lf  = cdv%lf
    lef = cdv%lef

    ! High-pass for detail coefficients:
    d(1:N,      :) = matmul (cdv_le_hi(1:N,1:lef), v(1:lef,    :)) ! Left edge
    d(M2-N+1:M2,:) = matmul (cdv_re_hi(1:N,1:lef), v(M-lef+1:M,:)) ! Right edge

    ! Low-pass for coefficients at next coarser scale:
    s(1:N,      :) = matmul (cdv_le_lo(1:N,1:lef), v(1:lef,    :)) ! Left edge
    s(M2-N+1:M2,:) = matmul (cdv_re_lo(1:N,1:lef), v(M-lef+1:M,:)) ! Right edge

#ifdef __NEC__          /* vector version */
    d(N+1:M2-N, :) = 0._wp
    s(N+1:M2-N, :) = 0._wp
!NEC$ loop_count(8)
    do i = 1, lf
       do k = N+1, M2-N
          do j = 1, P
             ! Interior filter: Convolution with downsampling
             d(k,j) = d(k,j) + cdv_c_inv(i) * v(2*k-N-1+i,j)
             ! Inner part: Low-pass filtering and decimation (downsampling)
             s(k,j) = s(k,j) + cdv_c_fwd(i) * v(2*k-N-1+i,j)
          end do
       end do
    end do
#else
!$omp parallel private(j,k)
!$omp do schedule(static)
    do j = 1, P         ! Performance note: this must remain as outer loop
       do k = N+1, M2-N
          ! Interior filter: Convolution with downsampling
          d(k,j) = dot_product (cdv_c_inv(1:lf), v(2*k-N:2*k+N-1,j))
       end do

       do k = N+1, M2-N
          ! Inner part: Low-pass filtering and decimation (downsampling)
          s(k,j) = dot_product (cdv_c_fwd(1:lf), v(2*k-N:2*k+N-1,j))
       end do
    end do
!$omp end do
!$omp end parallel
#endif
  end subroutine cdv_dyad_down

  ! --

  subroutine cdv_dyad_up (s, d, v)
    ! Perform one iteration in the cascade algorithm (synthesis step)
    real(wp), intent(out)  :: v(:,:)
    real(wp), intent(in),  dimension (size (v,dim=1)/2,size (v,dim=2)) :: s, d

    ! Local variables
    integer :: M, M2, N, lf, lef, k, j, P
#ifdef __NEC__
    integer :: i
#endif

    M2  = size (s, dim=1)
    M   = 2*M2
    P   = size (s, dim=2)
    N   = cdv%N
    lf  = cdv%lf
    lef = cdv%lef

    ! Reconstruction (inverse transform):
    v(lef+1:,   :) = 0
    ! Left edge
!### workaround for bug in sxf90 revs.381,393,400 (no problem with rev.360)
!!NEC$ novector
    v(1:lef,    :) = matmul (transpose (cdv_le_lo(1:N,1:lef)), s(1:N,:)) &
                   + matmul (transpose (cdv_le_hi(1:N,1:lef)), d(1:N,:))
    ! Right edge (beware of overlap between edge functions for M <= 6*N-3)
!### workaround for bug in sxf90 revs.381,393,400 (no problem with rev.360)
!!NEC$ novector
    v(M-lef+1:M,:) = v(M-lef+1:M,:) &
                   + matmul (transpose (cdv_re_lo(1:N,1:lef)), s(M2-N+1:M2,:))&
                   + matmul (transpose (cdv_re_hi(1:N,1:lef)), d(M2-N+1:M2,:))

    ! Interior reconstruction
#ifdef __NEC__          /* vector version */
    do k = N+1, M2-N
!NEC$ loop_count(8)
       do i = 1, lf
!NEC$ ivdep
          do j = 1, P
             v(2*k-N-1+i,j) = v(2*k-N-1+i,j)          &
                            + s(k,j) * cdv_c_fwd(i) &
                            + d(k,j) * cdv_c_inv(i)
          end do
       end do
    end do
#else
!$omp parallel do private(j,k) schedule(static)
    do j = 1, P         ! Performance note: this must remain as outer loop
       do k = N+1, M2-N
          v(2*k-N:2*k+N-1,j) = v(2*k-N:2*k+N-1,j)   &
                             + s(k,j) * cdv_c_fwd(1:lf) &
                             + d(k,j) * cdv_c_inv(1:lf)
       end do
    end do
!$omp end parallel do
#endif
  end subroutine cdv_dyad_up

  ! --
  ! Routines for transposed vectors (optimized for vector processors)
  ! --

  subroutine cdv_apply_pre_tr (v_t)
    real(wp), intent(inout) :: v_t(:,:)         ! Transposed input vector
    ! Local variables
    integer :: M, N
    M = size (v_t, dim=2)
    N = cdv%N
    ! Apply preconditioning matrices
    v_t(:,1:N    ) = matmul (v_t(:,1:N    ), transpose (cdv_l_pre(1:N,1:N)))
    v_t(:,M-N+1:M) = matmul (v_t(:,M-N+1:M), transpose (cdv_r_pre(1:N,1:N)))
  end subroutine cdv_apply_pre_tr

  ! --

  subroutine cdv_apply_pre_ad_tr (v_t)
    real(wp), intent(inout) :: v_t(:,:)         ! Transposed input vector
    ! Local variables
    integer :: M, N
    M = size (v_t, dim=2)
    N = cdv%N
    ! Apply adjoint preconditioning matrices
    v_t(:,1:N    ) = matmul (v_t(:,1:N    ), cdv_l_pre(1:N,1:N))
    v_t(:,M-N+1:M) = matmul (v_t(:,M-N+1:M), cdv_r_pre(1:N,1:N))
  end subroutine cdv_apply_pre_ad_tr

  ! --

  subroutine cdv_apply_post_tr (v_t)
    real(wp), intent(inout) :: v_t(:,:)         ! Transposed input vector
    ! Local variables
    integer  :: M, N
#if defined (__NEC__)
    integer  :: j, k
    real(wp) :: tmp(size (v_t,dim=1), cdv%N)
#endif
    M = size (v_t, dim=2)
    N = cdv%N
    !--------------------------------
    ! Apply postconditioning matrices
    !--------------------------------
#if defined (__NEC__)
    ! Force SX compiler to vectorize over first index of v_t
    tmp = 0._wp
    do j = 1, N
       do k = 1, N
          tmp(:,j) = tmp(:,j) + v_t(:,k    ) * cdv_l_post(j,k)
       end do
    end do
    v_t(:,1:N    ) = tmp
    tmp = 0._wp
    do j = 1, N
       do k = 1, N
          tmp(:,j) = tmp(:,j) + v_t(:,M-N+k) * cdv_r_post(j,k)
       end do
    end do
    v_t(:,M-N+1:M) = tmp
#else
    v_t(:,1:N    ) = matmul (v_t(:,1:N    ), transpose (cdv_l_post(1:N,1:N)))
    v_t(:,M-N+1:M) = matmul (v_t(:,M-N+1:M), transpose (cdv_r_post(1:N,1:N)))
#endif
  end subroutine cdv_apply_post_tr

  ! --

  subroutine cdv_apply_post_ad_tr (v_t)
    real(wp), intent(inout) :: v_t(:,:)         ! Transposed input vector
    ! Local variables
    integer  :: M, N
#if defined (__NEC__)
    integer  :: j, k
    real(wp) :: tmp(size (v_t,dim=1), cdv%N)
#endif
    M = size (v_t, dim=2)
    N = cdv%N
    !----------------------------------------
    ! Apply adjoint postconditioning matrices
    !----------------------------------------
#if defined (__NEC__)
    ! Force SX compiler to vectorize over first index of v_t
    tmp = 0._wp
    do j = 1, N
       do k = 1, N
          tmp(:,j) = tmp(:,j) + v_t(:,k    ) * cdv_l_post(k,j)
       end do
    end do
    v_t(:,1:N    ) = tmp
    tmp = 0._wp
    do j = 1, N
       do k = 1, N
          tmp(:,j) = tmp(:,j) + v_t(:,M-N+k) * cdv_r_post(k,j)
       end do
    end do
    v_t(:,M-N+1:M) = tmp
#else
    v_t(:,1:N    ) = matmul (v_t(:,1:N    ), cdv_l_post(1:N,1:N))
    v_t(:,M-N+1:M) = matmul (v_t(:,M-N+1:M), cdv_r_post(1:N,1:N))
#endif
  end subroutine cdv_apply_post_ad_tr

  ! --

  subroutine cdv_dyad_down_tr (v, s, d)
    ! Perform one iteration in the cascade algorithm (analysis step)
    ! (Transposed argument vectors)
    real(wp), intent(in)   :: v(:,:)
!#if defined (__NEC__)
!   ! Work around bug in optimizer of sxf90 Revs.371,372
!   real(wp), intent(out), dimension(ubound(v,dim=1),ubound(v,dim=2)/2):: s, d
!#else
    ! This appears to work with sxf90 Revs.360,381
    real(wp), intent(out), dimension (size (v,dim=1),size (v,dim=2)/2) :: s, d
!#endif

    ! Local variables
    integer :: M, M2, N, lf, lef, k, j, P

    P   = size (v, dim=1)
    M   = size (v, dim=2)
    M2  = M/2
    N   = cdv%N
    lf  = cdv%lf
    lef = cdv%lef

    ! High-pass for detail coefficients:
!   d(:,1:N      ) = matmul (v(:,1:lef    ), &
!                            transpose (cdv_le_hi(1:N,1:lef))) ! Left edge
!   d(:,M2-N+1:M2) = matmul (v(:,M-lef+1:M), &
!                            transpose (cdv_re_hi(1:N,1:lef))) ! Right edge
    ! Work around SX compiler bug with matmul:
    call dgemm ('n', 't', P, N, lef, 1._wp, &
                 v(:,1:lef    ), P, cdv_le_hi, N, 0._wp, d(:,1:N      ), P)
    call dgemm ('n', 't', P, N, lef, 1._wp, &
                 v(:,M-lef+1:M), P, cdv_re_hi, N, 0._wp, d(:,M2-N+1:M2), P)

    do k = N+1, M2-N
#ifndef __NEC__         /* optimize for SX Aurora again */
       do j = 1, P
          ! Interior filter: Convolution with downsampling
          d(j,k) = dot_product (v(j,2*k-N:2*k+N-1), cdv_c_inv(1:lf))
       end do
#else
       ! Interior filter: Convolution with downsampling
       d(:,k) = matmul (v(:,2*k-N:2*k+N-1), cdv_c_inv(1:lf)) ! Slower on SX-6!
#endif
    end do

    ! Low-pass for coefficients at next coarser scale:
!   s(:,1:N      ) = matmul (v(:,1:lef    ), &
!                            transpose (cdv_le_lo(1:N,1:lef))) ! Left edge
!   s(:,M2-N+1:M2) = matmul (v(:,M-lef+1:M), &
!                            transpose (cdv_re_lo(1:N,1:lef))) ! Right edge
    ! Work around SX compiler bug with matmul:
    call dgemm ('n', 't', P, N, lef, 1._wp, &
                 v(:,1:lef    ), P, cdv_le_lo, N, 0._wp, s(:,1:N      ), P)
    call dgemm ('n', 't', P, N, lef, 1._wp, &
                 v(:,M-lef+1:M), P, cdv_re_lo, N, 0._wp, s(:,M2-N+1:M2), P)

    do k = N+1, M2-N
#ifndef __NEC__         /* optimize for SX Aurora again */
       do j = 1, P
          ! Inner part: Low-pass filtering and decimation (downsampling)
          s(j,k) = dot_product (v(j,2*k-N:2*k+N-1), cdv_c_fwd(1:lf))
       end do
#else
       ! Inner part: Low-pass filtering and decimation (downsampling)
       s(:,k) = matmul (v(:,2*k-N:2*k+N-1), cdv_c_fwd(1:lf)) ! Slower on SX-6!
#endif
    end do

  end subroutine cdv_dyad_down_tr

  ! --

  subroutine cdv_dyad_up_tr (s, d, v)
    ! Perform one iteration in the cascade algorithm (synthesis step)
    ! (Transposed argument vectors)
    real(wp), intent(out) :: v(:,:)
    real(wp), intent(in), dimension (size (v,dim=1),size (v,dim=2)/2) :: s, d

    ! Local variables
    integer :: M, M2, N, lf, lef, k, j, P
#ifdef __NEC__
    integer :: i
#endif

    P   = size (s, dim=1)
    M2  = size (s, dim=2)
    M   = 2*M2
    N   = cdv%N
    lf  = cdv%lf
    lef = cdv%lef

    ! Reconstruction (inverse transform):
    v(:,lef+1:   ) = 0
    ! Left edge
    ! (Workaround for NEC's compiler matmul bug!)
!!CDIR NOVECTOR
!   v(:,1:lef    ) = matmul (s(:,1:N), cdv_le_lo(1:N,1:lef)) &
!                  + matmul (d(:,1:N), cdv_le_hi(1:N,1:lef))
    call dgemm ('n', 'n', P, lef, N, 1._wp, &
                 s(:,1:N      ), P, cdv_le_lo, N, 0._wp, v(:,1:lef    ), P)
    call dgemm ('n', 'n', P, lef, N, 1._wp, &
                 d(:,1:N      ), P, cdv_le_hi, N, 1._wp, v(:,1:lef    ), P)

    ! Right edge (beware of overlap between edge functions for M <= 6*N-3)
    ! (Workaround for NEC's compiler matmul bug!)
!!CDIR NOVECTOR
!   v(:,M-lef+1:M) = v(:,M-lef+1:M) &
!                  + matmul (s(:,M2-N+1:M2), cdv_re_lo(1:N,1:lef)) &
!                  + matmul (d(:,M2-N+1:M2), cdv_re_hi(1:N,1:lef))
    call dgemm ('n', 'n', P, lef, N, 1._wp, &
                 s(:,M2-N+1:M2), P, cdv_re_lo, N, 1._wp, v(:,M-lef+1:M), P)
    call dgemm ('n', 'n', P, lef, N, 1._wp, &
                 d(:,M2-N+1:M2), P, cdv_re_hi, N, 1._wp, v(:,M-lef+1:M), P)

    ! Interior reconstruction
#ifdef __NEC__          /* vector version */
    do k = N+1, M2-N
!NEC$ loop_count(8)
       do i = 1, lf
!NEC$ ivdep
          do j = 1, P
             v(j,2*k-N-1+i) = v(j,2*k-N-1+i)        &
                            + s(j,k) * cdv_c_fwd(i) &
                            + d(j,k) * cdv_c_inv(i)
          end do
       end do
    end do
#else
    do k = N+1, M2-N
       do j = 1, P
          v(j,2*k-N:2*k+N-1) = v(j,2*k-N:2*k+N-1)   &
                             + s(j,k) * cdv_c_fwd(1:lf) &
                             + d(j,k) * cdv_c_inv(1:lf)
       end do
    end do
#endif
  end subroutine cdv_dyad_up_tr

  ! --

  subroutine cdjv4_setup ()
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
    !                                                                    c
    ! Daubechies' Interval QMF's with lf = 4                             c
    !                                                                    c
    ! Coefficients by M.E. Brewster and G. Beylkin                       c
    !                                                                    c
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
    implicit none

    REAL(wp), PARAMETER :: d4h1 = 0.4829629131445341_wp
    REAL(wp), PARAMETER :: d4h2 = 0.8365163037378079_wp
    REAL(wp), PARAMETER :: d4h3 = 0.2241438680420134_wp
    REAL(wp), PARAMETER :: d4h4 =-0.1294095225512604_wp

    REAL(wp), PARAMETER :: d4g1 =  d4h4
    REAL(wp), PARAMETER :: d4g2 = -d4h3
    REAL(wp), PARAMETER :: d4g3 =  d4h2
    REAL(wp), PARAMETER :: d4g4 = -d4h1

    ! Left edge coefficients for the wavelets on the interval

    REAL(wp), PARAMETER :: d4h1l1 = 6.033325119280526E-01_wp
    REAL(wp), PARAMETER :: d4h1l2 = 6.908955318391036E-01_wp
    REAL(wp), PARAMETER :: d4h1l3 =-3.983129976982278E-01_wp

    REAL(wp), PARAMETER :: d4g1l1 =-7.965435169121828E-01_wp
    REAL(wp), PARAMETER :: d4g1l2 = 5.463927139590155E-01_wp
    REAL(wp), PARAMETER :: d4g1l3 =-2.587922483338183E-01_wp

    REAL(wp), PARAMETER :: d4h2l1 = 3.751746045244657E-02_wp
    REAL(wp), PARAMETER :: d4h2l2 = 4.573276598517686E-01_wp
    REAL(wp), PARAMETER :: d4h2l3 = 8.500881025491647E-01_wp
    REAL(wp), PARAMETER :: d4h2l4 = 2.238203569831144E-01_wp
    REAL(wp), PARAMETER :: d4h2l5 =-1.292227433543192E-01_wp

    REAL(wp), PARAMETER :: d4g2l1 =-1.003722456441388E-02_wp
    REAL(wp), PARAMETER :: d4g2l2 =-1.223510431167988E-01_wp
    REAL(wp), PARAMETER :: d4g2l3 =-2.274281116558368E-01_wp
    REAL(wp), PARAMETER :: d4g2l4 = 8.366029212236539E-01_wp
    REAL(wp), PARAMETER :: d4g2l5 =-4.830129217733038E-01_wp

    ! Right edge coefficients for the wavelets on the interval

    REAL(wp), PARAMETER :: d4h1r1 = 8.705087533498658E-01_wp
    REAL(wp), PARAMETER :: d4h1r2 = 4.348969979657030E-01_wp
    REAL(wp), PARAMETER :: d4h1r3 = 2.303890437969692E-01_wp

    REAL(wp), PARAMETER :: d4g1r1 = 2.575129194784820E-01_wp
    REAL(wp), PARAMETER :: d4g1r2 =-8.014229619903372E-01_wp
    REAL(wp), PARAMETER :: d4g1r3 = 5.398225007317715E-01_wp

    REAL(wp), PARAMETER :: d4h2r1 =-1.942334074274121E-01_wp
    REAL(wp), PARAMETER :: d4h2r2 = 1.901514184299554E-01_wp
    REAL(wp), PARAMETER :: d4h2r3 = 3.749553316456865E-01_wp
    REAL(wp), PARAMETER :: d4h2r4 = 7.675566692981142E-01_wp
    REAL(wp), PARAMETER :: d4h2r5 = 4.431490496375588E-01_wp

    REAL(wp), PARAMETER :: d4g2r1 = 3.717189665352956E-01_wp
    REAL(wp), PARAMETER :: d4g2r2 =-3.639069595708908E-01_wp
    REAL(wp), PARAMETER :: d4g2r3 =-7.175799993537222E-01_wp
    REAL(wp), PARAMETER :: d4g2r4 = 4.010695194302172E-01_wp
    REAL(wp), PARAMETER :: d4g2r5 = 2.315575950067897E-01_wp

    ! Left edge coefficients for the preconditioner

    REAL(wp), PARAMETER :: d4pc1l1 = 3.248940488989625E-01_wp
    REAL(wp), PARAMETER :: d4pc1l2 = 3.715801511588045E-02_wp
    REAL(wp), PARAMETER :: d4pc2l2 = 1.001445404981297E+00_wp

    REAL(wp), PARAMETER :: d4ipc1l1 = 3.077926491386692E+00_wp
    REAL(wp), PARAMETER :: d4ipc1l2 =-1.142045672421369E-01_wp
    REAL(wp), PARAMETER :: d4ipc2l2 = 9.985566811988880E-01_wp

    ! Right edge coefficients for the preconditioner

    REAL(wp), PARAMETER :: d4pc1r1 = 2.096292884353243E+00_wp
    REAL(wp), PARAMETER :: d4pc1r2 =-8.008132342464376E-01_wp
    REAL(wp), PARAMETER :: d4pc2r2 = 1.089843052895043E+00_wp

    REAL(wp), PARAMETER :: d4ipc1r1 = 4.770325785409152E-01_wp
    REAL(wp), PARAMETER :: d4ipc1r2 = 3.505220325509179E-01_wp
    REAL(wp), PARAMETER :: d4ipc2r2 = 9.175633109222603E-01_wp

    integer :: n, lef

    call cdv_clear ()
    n = 2
    call cdv_init (n)
    lef = cdv%lef   ! Maximum edge filter length

    ! Low- and high-pass filter coefficients
    cdv_c_fwd(:) = (/ d4h1, d4h2, d4h3, d4h4 /)
    cdv_c_inv(:) = (/ d4g1, d4g2, d4g3, d4g4 /)

    ! Preconditioner of forward transform (analysis)
    cdv_l_pre(:,:) = reshape ( &
         (/ d4pc1l1,  d4pc1l2,  0._wp, d4pc2l2 /), (/n, n/), order=(/2,1/) )
    cdv_r_pre(:,:) = reshape ( &
         (/ d4pc1r1,  d4pc1r2,  0._wp, d4pc2r2 /), (/n, n/), order=(/2,1/) )

    ! Postconditioner of inverse transform (synthesis)
    cdv_l_post(:,:) = reshape ( &
         (/ d4ipc1l1, d4ipc1l2, 0._wp, d4ipc2l2 /), (/n, n/), order=(/2,1/) )
    cdv_r_post(:,:) = reshape ( &
         (/ d4ipc1r1, d4ipc1r2, 0._wp, d4ipc2r2 /), (/n, n/), order=(/2,1/) )

    ! Left-edge low-pass
    cdv_le_lo(:,:) = reshape ( (/ &
         d4h1l1, d4h1l2, d4h1l3, 0._wp,  0._wp, &
         d4h2l1, d4h2l2, d4h2l3, d4h2l4, d4h2l5 &
         /), (/ n, lef /), order=(/2,1/) )
    ! Left-edge high-pass
    cdv_le_hi(:,:) = reshape ( (/ &
         d4g1l1, d4g1l2, d4g1l3, 0._wp,  0._wp, &
         d4g2l1, d4g2l2, d4g2l3, d4g2l4, d4g2l5   &
         /), (/ n, lef /), order=(/2,1/) )

    ! Right-edge low-pass
    cdv_re_lo(:,:) = reshape ( (/ &
         d4h1r1, d4h1r2, d4h1r3, 0._wp,  0._wp, &
         d4h2r1, d4h2r2, d4h2r3, d4h2r4, d4h2r5 &
         /), (/ n, lef /), order=(/2,1/) )
    ! Right-edge high-pass
    cdv_re_hi(:,:) = reshape ( (/ &
         d4g1r1, d4g1r2, d4g1r3, 0._wp,  0._wp, &
         d4g2r1, d4g2r2, d4g2r3, d4g2r4, d4g2r5 &
         /), (/ n, lef /), order=(/2,1/) )

    ! Flip arrays affecting right edge, to avoid flipping vectors...
    cdv_re_lo(:,:)  = cdv_re_lo(n:1:-1,lef:1:-1)
    cdv_re_hi(:,:)  = cdv_re_hi(n:1:-1,lef:1:-1)
    cdv_r_pre(:,:)  = cdv_r_pre(n:1:-1,n:1:-1)
    cdv_r_post(:,:) = cdv_r_post(n:1:-1,n:1:-1)

    !print *, "=== Setup of cdjv_4 complete ==="
  end subroutine cdjv4_setup

  ! --

  subroutine cdjv8_setup ()
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
    !                                                                    c
    ! Daubechies' Interval QMF's with lf = 8                             c
    !                                                                    c
    ! Coefficients by M.E. Brewster and G. Beylkin                       c
    !                                                                    c
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
    implicit none

    REAL(wp), PARAMETER :: d8h1 = 3.22231006040514678716159E-02_wp
    REAL(wp), PARAMETER :: d8h2 =-1.26039672620313037539160E-02_wp
    REAL(wp), PARAMETER :: d8h3 =-9.92195435766335325852080E-02_wp
    REAL(wp), PARAMETER :: d8h4 = 2.97857795605306051402901E-01_wp
    REAL(wp), PARAMETER :: d8h5 = 8.03738751805132080878805E-01_wp
    REAL(wp), PARAMETER :: d8h6 = 4.97618667632774989979605E-01_wp
    REAL(wp), PARAMETER :: d8h7 =-2.96355276460024917643691E-02_wp
    REAL(wp), PARAMETER :: d8h8 =-7.57657147895022132277461E-02_wp

    REAL(wp), PARAMETER :: d8g1 = d8h8
    REAL(wp), PARAMETER :: d8g2 =-d8h7
    REAL(wp), PARAMETER :: d8g3 = d8h6
    REAL(wp), PARAMETER :: d8g4 =-d8h5
    REAL(wp), PARAMETER :: d8g5 = d8h4
    REAL(wp), PARAMETER :: d8g6 =-d8h3
    REAL(wp), PARAMETER :: d8g7 = d8h2
    REAL(wp), PARAMETER :: d8g8 =-d8h1

    ! Left edge coefficients for the wavelets on the interval

    REAL(wp), PARAMETER :: d8h1l1 = .90975392578299541422_wp
    REAL(wp), PARAMETER :: d8h1l2 = .40416588940201723557_wp
    REAL(wp), PARAMETER :: d8h1l3 = .089040318656793899574_wp
    REAL(wp), PARAMETER :: d8h1l4 = -.011984192014966211145_wp
    REAL(wp), PARAMETER :: d8h1l5 = -.030429084139182777526_wp

    REAL(wp), PARAMETER :: d8g1l1 = .075739707618965346219_wp
    REAL(wp), PARAMETER :: d8g1l2 = -.32543917184495103305_wp
    REAL(wp), PARAMETER :: d8g1l3 = .68434908047295674723_wp
    REAL(wp), PARAMETER :: d8g1l4 = -.62004421070554711214_wp
    REAL(wp), PARAMETER :: d8g1l5 = .18858513977781955539_wp

    REAL(wp), PARAMETER :: d8h2l1 = -.27285140766949225591_wp
    REAL(wp), PARAMETER :: d8h2l2 = .50908152323149266343_wp
    REAL(wp), PARAMETER :: d8h2l3 = .62364244337236631118_wp
    REAL(wp), PARAMETER :: d8h2l4 = .46284008632544893219_wp
    REAL(wp), PARAMETER :: d8h2l5 = .24674764172786638137_wp
    REAL(wp), PARAMETER :: d8h2l6 = -.017669533291129221108_wp
    REAL(wp), PARAMETER :: d8h2l7 = -.045173645490327315413_wp

    REAL(wp), PARAMETER :: d8g2l1 = .16659595920609478902_wp
    REAL(wp), PARAMETER :: d8g2l2 = -.48478430891405802583_wp
    REAL(wp), PARAMETER :: d8g2l3 = .35646354251378829101_wp
    REAL(wp), PARAMETER :: d8g2l4 = .48398961557429770139_wp
    REAL(wp), PARAMETER :: d8g2l5 = -.60575436513594114671_wp
    REAL(wp), PARAMETER :: d8g2l6 = .034518332894524950210_wp
    REAL(wp), PARAMETER :: d8g2l7 = .088249016394632876783_wp

    REAL(wp), PARAMETER :: d8h3l1 = .12611792859809708119_wp
    REAL(wp), PARAMETER :: d8h3l2 = -.23085572680672168167_wp
    REAL(wp), PARAMETER :: d8h3l3 = -.052799235246625577463_wp
    REAL(wp), PARAMETER :: d8h3l4 = .21926517128941860346_wp
    REAL(wp), PARAMETER :: d8h3l5 = .46348072109960973065_wp
    REAL(wp), PARAMETER :: d8h3l6 = .70011971404312468766_wp
    REAL(wp), PARAMETER :: d8h3l7 = .41203257897604627778_wp
    REAL(wp), PARAMETER :: d8h3l8 = -.026222762497916067111_wp
    REAL(wp), PARAMETER :: d8h3l9 = -.067040694133818099667_wp

    REAL(wp), PARAMETER :: d8g3l1 = .20825353261585996816_wp
    REAL(wp), PARAMETER :: d8g3l2 = -.40182279321111657707_wp
    REAL(wp), PARAMETER :: d8g3l3 = -.068721487621235874441_wp
    REAL(wp), PARAMETER :: d8g3l4 = .33021351133522513671_wp
    REAL(wp), PARAMETER :: d8g3l5 = .55802129618734586368_wp
    REAL(wp), PARAMETER :: d8g3l6 = -.59949741336936802557_wp
    REAL(wp), PARAMETER :: d8g3l7 = -.069091991978684908007_wp
    REAL(wp), PARAMETER :: d8g3l8 = .027853569971929078366_wp
    REAL(wp), PARAMETER :: d8g3l9 = .071209990372730355011_wp

    REAL(wp), PARAMETER :: d8h4l1 = -.029079804272549243429_wp
    REAL(wp), PARAMETER :: d8h4l2 = .059928072294318161380_wp
    REAL(wp), PARAMETER :: d8h4l3 = .0061764277780953420663_wp
    REAL(wp), PARAMETER :: d8h4l4 = -.040210999043085892199_wp
    REAL(wp), PARAMETER :: d8h4l5 = -.039525870128122906889_wp
    REAL(wp), PARAMETER :: d8h4l6 = -.052599062573631102042_wp
    REAL(wp), PARAMETER :: d8h4l7 = .32894944796959845562_wp
    REAL(wp), PARAMETER :: d8h4l8 = .79663789672531588080_wp
    REAL(wp), PARAMETER :: d8h4l9 = .49011303364021253823_wp
    REAL(wp), PARAMETER :: d8h4l10 = -.029432877677387671779_wp
    REAL(wp), PARAMETER :: d8h4l11 = -.075247623129128381429_wp

    REAL(wp), PARAMETER :: d8g4l1 = -.065485007015422398501_wp
    REAL(wp), PARAMETER :: d8g4l2 = .13495242945354731343_wp
    REAL(wp), PARAMETER :: d8g4l3 = .013908739295079390932_wp
    REAL(wp), PARAMETER :: d8g4l4 = -.090551419457775651279_wp
    REAL(wp), PARAMETER :: d8g4l5 = -.089008573041674695765_wp
    REAL(wp), PARAMETER :: d8g4l6 = -.37334444756195190139_wp
    REAL(wp), PARAMETER :: d8g4l7 = .84046537082509037751_wp
    REAL(wp), PARAMETER :: d8g4l8 = -.31568493615246357870_wp
    REAL(wp), PARAMETER :: d8g4l9 = -.12029765089741391514_wp
    REAL(wp), PARAMETER :: d8g4l10 = .013070202799775884201_wp
    REAL(wp), PARAMETER :: d8g4l11 = .033415070904004976114_wp

    ! Right edge coefficients for the wavelets on the interval

    REAL(wp), PARAMETER :: d8h1r1 = .91547051883810854018_wp
    REAL(wp), PARAMETER :: d8h1r2 = .39191428098156657169_wp
    REAL(wp), PARAMETER :: d8h1r3 = .059477711238407501787_wp
    REAL(wp), PARAMETER :: d8h1r4 = -.025191808506643957308_wp
    REAL(wp), PARAMETER :: d8h1r5 = .064379345686262128356_wp

    REAL(wp), PARAMETER :: d8g1r1 = .19827799059709502246_wp
    REAL(wp), PARAMETER :: d8g1r2 = -.60406778678838103018_wp
    REAL(wp), PARAMETER :: d8g1r3 = .64952976030516477709_wp
    REAL(wp), PARAMETER :: d8g1r4 = -.40503095407943504406_wp
    REAL(wp), PARAMETER :: d8g1r5 = .099241947405233368263_wp

    REAL(wp), PARAMETER :: d8h2r1 = -.21916264686360397441_wp
    REAL(wp), PARAMETER :: d8h2r2 = .44880017813063045390_wp
    REAL(wp), PARAMETER :: d8h2r3 = .75400050839596823140_wp
    REAL(wp), PARAMETER :: d8h2r4 = .39377581566148591390_wp
    REAL(wp), PARAMETER :: d8h2r5 = -.15813389438931389654_wp
    REAL(wp), PARAMETER :: d8h2r6 = -.016142011904287920593_wp
    REAL(wp), PARAMETER :: d8h2r7 = .041268408805739583103_wp

    REAL(wp), PARAMETER :: d8g2r1 = -.27262735058992725801_wp
    REAL(wp), PARAMETER :: d8g2r2 = .50928674835762860473_wp
    REAL(wp), PARAMETER :: d8g2r3 = .068118566264738166854_wp
    REAL(wp), PARAMETER :: d8g2r4 = -.67353457518648896953_wp
    REAL(wp), PARAMETER :: d8g2r5 = .44993408265114559900_wp
    REAL(wp), PARAMETER :: d8g2r6 = .027190655403678115961_wp
    REAL(wp), PARAMETER :: d8g2r7 = -.069515193617030161473_wp

    REAL(wp), PARAMETER :: d8h3r1 = .012900782893843634354_wp
    REAL(wp), PARAMETER :: d8h3r2 = -.13907160056079416442_wp
    REAL(wp), PARAMETER :: d8h3r3 = .029213679496076078323_wp
    REAL(wp), PARAMETER :: d8h3r4 = .46061685367679247208_wp
    REAL(wp), PARAMETER :: d8h3r5 = .81641197419298729535_wp
    REAL(wp), PARAMETER :: d8h3r6 = .29864733458521474811_wp
    REAL(wp), PARAMETER :: d8h3r7 = -.10276635357182394708_wp
    REAL(wp), PARAMETER :: d8h3r8 = -.012574882110017058721_wp
    REAL(wp), PARAMETER :: d8h3r9 = .032148741970777129858_wp

    REAL(wp), PARAMETER :: d8g3r1 = -.0045819593397576190533_wp
    REAL(wp), PARAMETER :: d8g3r2 = -.030620323648208387712_wp
    REAL(wp), PARAMETER :: d8g3r3 = -.013887371452971361363_wp
    REAL(wp), PARAMETER :: d8g3r4 = .095048349988442514248_wp
    REAL(wp), PARAMETER :: d8g3r5 = .30158150124129352334_wp
    REAL(wp), PARAMETER :: d8g3r6 = -.80336363165658830946_wp
    REAL(wp), PARAMETER :: d8g3r7 = .49684270650616092637_wp
    REAL(wp), PARAMETER :: d8g3r8 = .029632043697507898216_wp
    REAL(wp), PARAMETER :: d8g3r9 = -.075756807782644236930_wp

    REAL(wp), PARAMETER :: d8h4r1 = -.0067756036517274764721_wp
    REAL(wp), PARAMETER :: d8h4r2 = .019132441289899007282_wp
    REAL(wp), PARAMETER :: d8h4r3 = -.017709184245907807129_wp
    REAL(wp), PARAMETER :: d8h4r4 = -.067659161740038464179_wp
    REAL(wp), PARAMETER :: d8h4r5 = -.030235884812807082946_wp
    REAL(wp), PARAMETER :: d8h4r6 = .49779208211653976270_wp
    REAL(wp), PARAMETER :: d8h4r7 = .80394959961412010659_wp
    REAL(wp), PARAMETER :: d8h4r8 = .29771110110016361257_wp
    REAL(wp), PARAMETER :: d8h4r9 = -.099108040550099571666_wp
    REAL(wp), PARAMETER :: d8h4r10 = -.012598951895850180172_wp
    REAL(wp), PARAMETER :: d8h4r11 = .032210278399291594766_wp

    REAL(wp), PARAMETER :: d8g4r1 = -.0028803051241456110658_wp
    REAL(wp), PARAMETER :: d8g4r2 = .0081331895307455420533_wp
    REAL(wp), PARAMETER :: d8g4r3 = -.0075281637990915038805_wp
    REAL(wp), PARAMETER :: d8g4r4 = -.028761869830674646117_wp
    REAL(wp), PARAMETER :: d8g4r5 = -.012853256836710187824_wp
    REAL(wp), PARAMETER :: d8g4r6 = .099201914392120836400_wp
    REAL(wp), PARAMETER :: d8g4r7 = .29778998412792446390_wp
    REAL(wp), PARAMETER :: d8g4r8 = -.80379368413398200528_wp
    REAL(wp), PARAMETER :: d8g4r9 = .49764705235425275855_wp
    REAL(wp), PARAMETER :: d8g4r10 = .029637660176292392689_wp
    REAL(wp), PARAMETER :: d8g4r11 = -.075771166782247360195_wp

    ! Left edge coefficients for the preconditioner

    REAL(wp), PARAMETER :: d8pc1l1 = 2.4899111140824833432_wp
    REAL(wp), PARAMETER :: d8pc1l2 =-2.7529885357135795429_wp
    REAL(wp), PARAMETER :: d8pc1l3 = 1.6878414467970039542_wp
    REAL(wp), PARAMETER :: d8pc1l4 =-.40222211735395905367_wp
    REAL(wp), PARAMETER :: d8pc2l2 = 1.6772105498044283777_wp
    REAL(wp), PARAMETER :: d8pc2l3 =-.70753754370482426899_wp
    REAL(wp), PARAMETER :: d8pc2l4 = .17635442878554937940_wp
    REAL(wp), PARAMETER :: d8pc3l3 = 1.1301451419680759567_wp
    REAL(wp), PARAMETER :: d8pc3l4 =-.061621215501286551203_wp
    REAL(wp), PARAMETER :: d8pc4l4 = 1.0068851564850727934_wp

    REAL(wp), PARAMETER :: d8ipc1l1 =  .40162076242167132109_wp
    REAL(wp), PARAMETER :: d8ipc1l2 =  .65922394465043988597_wp
    REAL(wp), PARAMETER :: d8ipc1l3 =-.18709674563738923634_wp
    REAL(wp), PARAMETER :: d8ipc1l4 = .033523746113502769821_wp
    REAL(wp), PARAMETER :: d8ipc2l2 = .59622806457818029458_wp
    REAL(wp), PARAMETER :: d8ipc2l3 = .37327394918930106812_wp
    REAL(wp), PARAMETER :: d8ipc2l4 =-.081584145680874628305_wp
    REAL(wp), PARAMETER :: d8ipc3l3 = .88484209935952433399_wp
    REAL(wp), PARAMETER :: d8ipc3l4 = .054152199322895070106_wp
    REAL(wp), PARAMETER :: d8ipc4l4 = .99316192473319585900_wp

    ! Right edge coefficients for the preconditioner

    REAL(wp), PARAMETER :: d8pc1r1 = .50051923113793834074_wp
    REAL(wp), PARAMETER :: d8pc1r2 = .37673863695833614763_wp
    REAL(wp), PARAMETER :: d8pc1r3 =-.00093100685477696553087_wp
    REAL(wp), PARAMETER :: d8pc1r4 =-.0073733049057012719128_wp
    REAL(wp), PARAMETER :: d8pc2r2 = .78081761658738585219_wp
    REAL(wp), PARAMETER :: d8pc2r3 = .091704627903655303768_wp
    REAL(wp), PARAMETER :: d8pc2r4 =-.018445047273518395674_wp
    REAL(wp), PARAMETER :: d8pc3r3 = 1.0023129562376633347_wp
    REAL(wp), PARAMETER :: d8pc3r4 =-.0022411542961159204922_wp
    REAL(wp), PARAMETER :: d8pc4r4 = 1.0003980780482839635_wp

    REAL(wp), PARAMETER :: d8ipc1r1 = 1.9979252300186034387_wp
    REAL(wp), PARAMETER :: d8ipc1r2 =-.96398392135616003463_wp
    REAL(wp), PARAMETER :: d8ipc1r3 = .090053578910487419586_wp
    REAL(wp), PARAMETER :: d8ipc1r4 =-.0028464600220995317358_wp
    REAL(wp), PARAMETER :: d8ipc2r2 = 1.2807088092742643304_wp
    REAL(wp), PARAMETER :: d8ipc2r3 =-.11717590207382437902_wp
    REAL(wp), PARAMETER :: d8ipc2r4 = .023350829801588023229_wp
    REAL(wp), PARAMETER :: d8ipc3r3 = .99769238118367204024_wp
    REAL(wp), PARAMETER :: d8ipc3r4 = .0022350928249024384967_wp
    REAL(wp), PARAMETER :: d8ipc4r4 = .99960208035479177484_wp

    integer :: n, lef

    call cdv_clear ()
    n = 4
    call cdv_init (n)
    lef = cdv%lef   ! Maximum edge filter length

    ! Low- and high-pass filter coefficients
    cdv_c_fwd(:) = (/ d8h1, d8h2, d8h3, d8h4, d8h5, d8h6, d8h7, d8h8 /)
    cdv_c_inv(:) = (/ d8g1, d8g2, d8g3, d8g4, d8g5, d8g6, d8g7, d8g8 /)

    ! Preconditioner of forward transform (analysis)
    cdv_l_pre(:,:) = reshape ( (/ &
         d8pc1l1, d8pc1l2, d8pc1l3, d8pc1l4, &
         0._wp,   d8pc2l2, d8pc2l3, d8pc2l4, &
         0._wp,   0._wp,   d8pc3l3, d8pc3l4, &
         0._wp,   0._wp,   0._wp,   d8pc4l4    /), (/n, n/), order=(/2,1/) )
    cdv_r_pre(:,:) = reshape ( (/ &
         d8pc1r1, d8pc1r2, d8pc1r3, d8pc1r4, &
         0._wp,   d8pc2r2, d8pc2r3, d8pc2r4, &
         0._wp,   0._wp,   d8pc3r3, d8pc3r4, &
         0._wp,   0._wp,   0._wp,   d8pc4r4    /), (/n, n/), order=(/2,1/) )

    ! Postconditioner of inverse transform (synthesis)
    cdv_l_post(:,:) = reshape ( (/ &
         d8ipc1l1, d8ipc1l2, d8ipc1l3, d8ipc1l4, &
         0._wp,    d8ipc2l2, d8ipc2l3, d8ipc2l4, &
         0._wp,    0._wp,    d8ipc3l3, d8ipc3l4, &
         0._wp,    0._wp,    0._wp,    d8ipc4l4    /), (/n,n/), order=(/2,1/) )
    cdv_r_post(:,:) = reshape ( (/ &
         d8ipc1r1, d8ipc1r2, d8ipc1r3, d8ipc1r4, &
         0._wp,    d8ipc2r2, d8ipc2r3, d8ipc2r4, &
         0._wp,    0._wp,    d8ipc3r3, d8ipc3r4, &
         0._wp,    0._wp,    0._wp,    d8ipc4r4    /), (/n,n/), order=(/2,1/) )

    ! Left-edge low-pass
    cdv_le_lo(:,:) = reshape ( (/ &
         d8h1l1, d8h1l2, d8h1l3, d8h1l4, d8h1l5, 0._wp,  0._wp,  0._wp, &
         0._wp,  0._wp,  0._wp, &
         d8h2l1, d8h2l2, d8h2l3, d8h2l4, d8h2l5, d8h2l6, d8h2l7, 0._wp, &
         0._wp,  0._wp,  0._wp, &
         d8h3l1, d8h3l2, d8h3l3, d8h3l4, d8h3l5, d8h3l6, d8h3l7, d8h3l8, &
         d8h3l9, 0._wp,  0._wp, &
         d8h4l1, d8h4l2, d8h4l3, d8h4l4, d8h4l5, d8h4l6, d8h4l7, d8h4l8, &
         d8h4l9, d8h4l10, d8h4l11 &
         /), (/ n, lef /), order=(/2,1/) )
    ! Left-edge high-pass
    cdv_le_hi(:,:) = reshape ( (/ &
         d8g1l1, d8g1l2, d8g1l3, d8g1l4, d8g1l5, 0._wp,  0._wp,  0._wp, &
         0._wp,  0._wp,  0._wp, &
         d8g2l1, d8g2l2, d8g2l3, d8g2l4, d8g2l5, d8g2l6, d8g2l7, 0._wp, &
         0._wp,  0._wp,  0._wp, &
         d8g3l1, d8g3l2, d8g3l3, d8g3l4, d8g3l5, d8g3l6, d8g3l7, d8g3l8, &
         d8g3l9, 0._wp,  0._wp, &
         d8g4l1, d8g4l2, d8g4l3, d8g4l4, d8g4l5, d8g4l6, d8g4l7, d8g4l8, &
         d8g4l9, d8g4l10, d8g4l11 &
         /), (/ n, lef /), order=(/2,1/) )

    ! Right-edge low-pass
    cdv_re_lo(:,:) = reshape ( (/ &
         d8h1r1, d8h1r2, d8h1r3, d8h1r4, d8h1r5, 0._wp,  0._wp,  0._wp, &
         0._wp,  0._wp,  0._wp, &
         d8h2r1, d8h2r2, d8h2r3, d8h2r4, d8h2r5, d8h2r6, d8h2r7, 0._wp, &
         0._wp,  0._wp,  0._wp, &
         d8h3r1, d8h3r2, d8h3r3, d8h3r4, d8h3r5, d8h3r6, d8h3r7, d8h3r8, &
         d8h3r9, 0._wp,  0._wp, &
         d8h4r1, d8h4r2, d8h4r3, d8h4r4, d8h4r5, d8h4r6, d8h4r7, d8h4r8, &
         d8h4r9, d8h4r10, d8h4r11 &
         /), (/ n, lef /), order=(/2,1/) )
    ! Right-edge high-pass
    cdv_re_hi(:,:) = reshape ( (/ &
         d8g1r1, d8g1r2, d8g1r3, d8g1r4, d8g1r5, 0._wp,  0._wp,  0._wp, &
         0._wp,  0._wp,  0._wp, &
         d8g2r1, d8g2r2, d8g2r3, d8g2r4, d8g2r5, d8g2r6, d8g2r7, 0._wp, &
         0._wp,  0._wp,  0._wp, &
         d8g3r1, d8g3r2, d8g3r3, d8g3r4, d8g3r5, d8g3r6, d8g3r7, d8g3r8, &
         d8g3r9, 0._wp,  0._wp, &
         d8g4r1, d8g4r2, d8g4r3, d8g4r4, d8g4r5, d8g4r6, d8g4r7, d8g4r8, &
         d8g4r9, d8g4r10, d8g4r11 &
         /), (/ n, lef /), order=(/2,1/) )

    ! Flip arrays affecting right edge, to avoid having to flip vectors later
    cdv_re_lo(:,:)  = cdv_re_lo(n:1:-1,lef:1:-1)
    cdv_re_hi(:,:)  = cdv_re_hi(n:1:-1,lef:1:-1)
    cdv_r_pre(:,:)  = cdv_r_pre(n:1:-1,n:1:-1)
    cdv_r_post(:,:) = cdv_r_post(n:1:-1,n:1:-1)

    !print *, "=== Setup of cdjv_8 complete ==="
  end subroutine cdjv8_setup

  ! --

  elemental function ispowerof2 (i)
    ! Returns true if integer argument i is a power of 2
    integer, intent(in) :: i
    logical             :: ispowerof2
    ispowerof2 = (iand (i, i-1) == 0)
  end function ispowerof2

  ! --

  elemental function ilog2 (i) result (l)
    ! Calculate the integer part of  max (log_2(i),0)
    integer, intent(in) :: i
    integer             :: l
    ! Local variables
    integer :: k
    l = 0               ! ilog2(0) = ilog2(1) = 0
    k = i
    do while (k > 1)    ! k >= 2
       l = l + 1
       k = k / 2
    end do
  end function ilog2

  ! --

  subroutine cdv_ana (x)
    ! Non-periodic wavelet analysis using Cohen-Daubechies-Vial wavelets
    real(wp), intent(inout) :: x(:,:)

    ! Local variables
    integer  :: M, j, n2, Mmin, nr
    real(wp), dimension (size (x,dim=1)/2,size(x,dim=2)) :: s, d
    real(wp), parameter :: sqrt_one_half = 0.70710678118654752440_wp

    if (cdv%N == 0) call cdv_setup ()

    M = size (x, dim=1)

    if (.not. ispowerof2 (M)) then
       print *, "cdv_ana: grid size is not power of 2:", M
       call finish ("cdv_ana","bad grid size")
    end if
    Mmin = 2*cdv%N
    if (M < Mmin) then
       print *, "cdv_ana: grid too small for chosen interval wavelet:", M
       call finish ("cdv_ana","grid too small")
    end if

    ! Preconditioning (except for Haar transformation that has order=1)
    if (cdv_precondition .and. cdv%N > 1) then
       call cdv_apply_pre (x)
    end if

    n2 = ilog2 (M)
    nr = M
    lvl_loop: do j = 1, n2
       nr = nr / 2
       if (nr >= Mmin) then
!### workaround for bug in sxf90 revs.393,400  (no problem with revs.360,381)
!CDIR NOIEXPAND
          call cdv_dyad_down (x(1:2*nr,:), s(1:nr,:), d(1:nr,:))
       else
          if (cdv_keep_scaling) exit lvl_loop
          select case (nr)
          case (1)      ! 1 remaining coarser-level coefficient: Haar
             d(1,:)   = sqrt_one_half * (x(1,:) - x(2,:))
             s(1,:)   = sqrt_one_half * (x(1,:) + x(2,:))
          case (2)      ! 2 remaining coarser-level coefficients: join 4->2
             d(1:2,:) = matmul (cdv_hi42, x(1:4,:))
             s(1:2,:) = matmul (cdv_lo42, x(1:4,:))
          case (4)      ! 4 remaining coarser-level coefficients: join 8->4
             d(1:4,:) = matmul (cdv_hi84, x(1:8,:))
             s(1:4,:) = matmul (cdv_lo84, x(1:8,:))
          case default
             print *, "cdv_ana: internal error trying to 'join':", nr, Mmin
             call finish ("cdv_ana","internal error trying to 'join'")
          end select
       end if
       x(1:nr,     :) = s(1:nr,:)
       x(nr+1:2*nr,:) = d(1:nr,:)
    end do lvl_loop

  end subroutine cdv_ana

  ! --

  subroutine cdv_syn (x)
    ! Non-periodic wavelet synthesis using Cohen-Daubechies-Vial wavelets
    real(wp), intent(inout) :: x(:,:)

    ! Local variables
    integer  :: M, j, n2, Mmin, nr
    real(wp), dimension (size (x,dim=1)/2,size(x,dim=2)) :: s, d
    real(wp), parameter :: sqrt_one_half = 0.70710678118654752440_wp

    if (cdv%N == 0) call cdv_setup ()

    M = size (x, dim=1)

    if (.not. ispowerof2 (M)) then
       print *, "cdv_syn: grid size is not power of 2:", M
       call finish ("cdv_syn","bad grid size")
    end if
    Mmin = 2*cdv%N
    if (M < Mmin) then
       print *, "cdv_syn: grid too small for chosen interval wavelet:", M
       call finish ("cdv_syn","grid too small")
    end if

    n2 = ilog2 (M)
    nr = 1
    lvl_loop: do j = 1, n2
       if (nr >= Mmin) then
          s(1:nr,:) = x(1:nr,     :)
          d(1:nr,:) = x(nr+1:2*nr,:)
!### workaround for bug in sxf90 revs.393,400  (no problem with revs.360,381)
!CDIR NOIEXPAND
          call cdv_dyad_up (s(1:nr,:), d(1:nr,:), x(1:2*nr,:))
       else
          if (cdv_keep_scaling) then
             nr = nr * 2
             cycle lvl_loop
          end if
          s(1:nr,:) = x(1:nr,     :)
          d(1:nr,:) = x(nr+1:2*nr,:)
          select case (nr)
          case (1)      ! 1 initial coarser-level coefficient: inverse Haar
             x(1,:)   = sqrt_one_half * (s(1,:) + d(1,:))
             x(2,:)   = sqrt_one_half * (s(1,:) - d(1,:))
          case (2)      ! 2 initial coarser-level coefficients: inverse join_2
             x(1:4,:) = matmul (transpose (cdv_lo42), s(1:2,:)) &
                      + matmul (transpose (cdv_hi42), d(1:2,:))
          case (4)      ! 4 initial coarser-level coefficients: inverse join_4
             x(1:8,:) = matmul (transpose (cdv_lo84), s(1:4,:)) &
                      + matmul (transpose (cdv_hi84), d(1:4,:))
          case default
             print *, "cdv_syn: internal error trying to 'join':", nr, Mmin
             call finish ("cdv_syn","internal error trying to 'join'")
          end select
       end if
       nr = nr * 2
    end do lvl_loop

    ! Inverse preconditioning (except for Haar transformation with order=1)
    if (cdv_precondition .and. cdv%N > 1) then
       call cdv_apply_post (x)
    end if

  end subroutine cdv_syn

  ! --

  subroutine cdv_ana_ad (x)
    ! Non-periodic adjoint wavelet analysis
    real(wp), intent(inout) :: x(:,:)

    ! Local variables
    integer  :: M, j, n2, Mmin, nr
    real(wp), dimension (size (x,dim=1)/2,size(x,dim=2)) :: s, d
    real(wp), parameter :: sqrt_one_half = 0.70710678118654752440_wp

    if (cdv%N == 0) call cdv_setup ()

    M = size (x, dim=1)

    if (.not. ispowerof2 (M)) then
       print *, "cdv_ana_ad: grid size is not power of 2:", M
       call finish ("cdv_ana_ad","bad grid size")
    end if
    Mmin = 2*cdv%N
    if (M < Mmin) then
       print *, "cdv_ana_ad: grid too small for chosen interval wavelet:", M
       call finish ("cdv_ana_ad","grid too small")
    end if

    n2 = ilog2 (M)
    nr = 1
    lvl_loop: do j = 1, n2
       if (nr >= Mmin) then
          s(1:nr,:) = x(1:nr,     :)
          d(1:nr,:) = x(nr+1:2*nr,:)
!### workaround for bug in sxf90 revs.393,400  (no problem with revs.360,381)
!CDIR NOIEXPAND
          call cdv_dyad_up (s(1:nr,:), d(1:nr,:), x(1:2*nr,:))
       else
          if (cdv_keep_scaling) then
             nr = nr * 2
             cycle lvl_loop
          end if
          s(1:nr,:) = x(1:nr,     :)
          d(1:nr,:) = x(nr+1:2*nr,:)
          select case (nr)
          case (1)      ! 1 initial coarser-level coefficient: inverse Haar
             x(1,:)   = sqrt_one_half * (s(1,:) + d(1,:))
             x(2,:)   = sqrt_one_half * (s(1,:) - d(1,:))
          case (2)      ! 2 initial coarser-level coefficients: inverse join_2
             x(1:4,:) = matmul (transpose (cdv_lo42), s(1:2,:)) &
                      + matmul (transpose (cdv_hi42), d(1:2,:))
          case (4)      ! 4 initial coarser-level coefficients: inverse join_4
             x(1:8,:) = matmul (transpose (cdv_lo84), s(1:4,:)) &
                      + matmul (transpose (cdv_hi84), d(1:4,:))
          case default
             print *, "cdv_ana_ad: internal error trying to 'join':", nr, Mmin
             call finish ("cdv_ana_ad","internal error trying to 'join'")
          end select
       end if
       nr = nr * 2
    end do lvl_loop

    ! Adjoint Preconditioning (not for Haar transformation)
    if (cdv_precondition .and. cdv%N > 1) then
       call cdv_apply_pre_ad (x)
    end if

  end subroutine cdv_ana_ad

  ! --

  subroutine cdv_syn_ad (x)
    ! Non-periodic adjoint wavelet synthesis
    real(wp), intent(inout) :: x(:,:)

    ! Local variables
    integer  :: M, j, n2, Mmin, nr
    real(wp), dimension (size (x,dim=1)/2,size(x,dim=2)) :: s, d
    real(wp), parameter :: sqrt_one_half = 0.70710678118654752440_wp

    if (cdv%N == 0) call cdv_setup ()

    M = size (x, dim=1)

    if (.not. ispowerof2 (M)) then
       print *, "cdv_syn_ad: grid size is not power of 2:", M
       call finish ("cdv_syn_ad","bad grid size")
    end if
    Mmin = 2*cdv%N
    if (M < Mmin) then
       print *, "cdv_syn_ad: grid too small for chosen interval wavelet:", M
       call finish ("cdv_syn_ad","grid too small")
    end if

    ! Adjoint Inverse Preconditioning (not for Haar transformation)
    if (cdv_precondition .and. cdv%N > 1) then
       call cdv_apply_post_ad (x)
    end if

    n2 = ilog2 (M)
    nr = M
    lvl_loop: do j = 1, n2
       nr = nr / 2
       if (nr >= Mmin) then
!### workaround for bug in sxf90 revs.393,400  (no problem with revs.360,381)
!CDIR NOIEXPAND
          call cdv_dyad_down (x(1:2*nr,:), s(1:nr,:), d(1:nr,:))
       else
          if (cdv_keep_scaling) exit lvl_loop
          select case (nr)
          case (1)      ! 1 remaining coarser-level coefficient: Haar
             d(1,:)   = sqrt_one_half * (x(1,:) - x(2,:))
             s(1,:)   = sqrt_one_half * (x(1,:) + x(2,:))
          case (2)      ! 2 remaining coarser-level coefficients: join 4->2
             d(1:2,:) = matmul (cdv_hi42, x(1:4,:))
             s(1:2,:) = matmul (cdv_lo42, x(1:4,:))
          case (4)      ! 4 remaining coarser-level coefficients: join 8->4
             d(1:4,:) = matmul (cdv_hi84, x(1:8,:))
             s(1:4,:) = matmul (cdv_lo84, x(1:8,:))
          case default
             print *, "cdv_syn_ad: internal error trying to 'join':", nr, Mmin
             call finish ("cdv_syn_ad","internal error trying to 'join'")
          end select
       end if
       x(1:nr,     :) = s(1:nr,:)
       x(nr+1:2*nr,:) = d(1:nr,:)
    end do lvl_loop

  end subroutine cdv_syn_ad

  ! --

  subroutine cdv_ana_vec (x)
    ! Non-periodic wavelet analysis using Cohen-Daubechies-Vial wavelets
    ! This variant should perform better on vector processors (e.g. SX-6)
    real(wp), intent(inout) :: x(:,:)

    ! Local variables
    integer  :: M, j, n2, Mmin, nr
    real(wp), dimension (size (x,dim=2), size (x,dim=1)  ) :: x_t
    real(wp), dimension (size (x,dim=2), size (x,dim=1)/2) :: s, d
    real(wp), parameter :: sqrt_one_half = 0.70710678118654752440_wp
#if defined (__NEC__)
    integer  :: P
    real(wp), dimension (size (x,dim=1)+1, size (x,dim=2)) :: x_
    real(wp), dimension (size (x,dim=2)+1, size (x,dim=1)) :: x_t_
#endif

    M = size (x, dim=1)

    if (.not. ispowerof2 (M)) then
       print *, "cdv_ana_vec: grid size is not power of 2:", M
       call finish ("cdv_ana_vec","bad grid size")
    end if
    Mmin = 2*cdv%N
    if (M < Mmin) then
       print *, "cdv_ana_vec: grid too small for chosen interval wavelet:", M
       call finish ("cdv_ana_vec","grid too small")
    end if

#if defined (__NEC__)
    x_(1:M,:) = x(1:M,:)
    x_t = transpose (x_(1:M,:))
#else
    x_t = transpose (x)
#endif

    ! Preconditioning (except for Haar transformation that has order=1)
    if (cdv_precondition .and. cdv%N > 1) then
       call cdv_apply_pre_tr (x_t)
    end if

    n2 = ilog2 (M)
    nr = M
    lvl_loop: do j = 1, n2
       nr = nr / 2
       if (nr >= Mmin) then
          call cdv_dyad_down_tr (x_t(:,1:2*nr), s(:,1:nr), d(:,1:nr))
       else
          if (cdv_keep_scaling) exit lvl_loop
          select case (nr)
          case (1)      ! 1 remaining coarser-level coefficient: Haar
             d(:,1)   = sqrt_one_half * (x_t(:,1) - x_t(:,2))
             s(:,1)   = sqrt_one_half * (x_t(:,1) + x_t(:,2))
          case (2)      ! 2 remaining coarser-level coefficients: join 4->2
             d(:,1:2) = matmul (x_t(:,1:4), transpose (cdv_hi42))
             s(:,1:2) = matmul (x_t(:,1:4), transpose (cdv_lo42))
          case (4)      ! 4 remaining coarser-level coefficients: join 8->4
             d(:,1:4) = matmul (x_t(:,1:8), transpose (cdv_hi84))
             s(:,1:4) = matmul (x_t(:,1:8), transpose (cdv_lo84))
          case default
             print *, "cdv_ana_vec: internal error trying to 'join':", nr, Mmin
             call finish ("cdv_ana_vec","internal error trying to 'join'")
          end select
       end if
       x_t(:,1:nr     ) = s(:,1:nr)
       x_t(:,nr+1:2*nr) = d(:,1:nr)
    end do lvl_loop

#if defined (__NEC__)
    P = size (x, dim=2)
    x_t_(1:P,:) = x_t(1:P,:)
    x = transpose (x_t_(1:P,:))
#else
    x = transpose (x_t)
#endif

  end subroutine cdv_ana_vec

  ! --

  subroutine cdv_syn_vec (x)
    ! Non-periodic wavelet synthesis using Cohen-Daubechies-Vial wavelets
    ! This variant should perform better on vector processors (e.g. SX-6)
    real(wp), intent(inout) :: x(:,:)

    ! Local variables
    integer  :: M, j, n2, Mmin, nr
    real(wp), dimension (size (x,dim=2), size (x,dim=1)  ) :: x_t
    real(wp), dimension (size (x,dim=2), size (x,dim=1)/2) :: s, d
    real(wp), parameter :: sqrt_one_half = 0.70710678118654752440_wp
#if defined (__NEC__)
    integer  :: P
    real(wp), dimension (size (x,dim=1)+1, size (x,dim=2)) :: x_
    real(wp), dimension (size (x,dim=2)+1, size (x,dim=1)) :: x_t_
#endif

    M = size (x, dim=1)

    if (.not. ispowerof2 (M)) then
       print *, "cdv_syn_vec: grid size is not power of 2:", M
       call finish ("cdv_syn_vec","bad grid size")
    end if
    Mmin = 2*cdv%N
    if (M < Mmin) then
       print *, "cdv_syn_vec: grid too small for chosen interval wavelet:", M
       call finish ("cdv_syn_vec","grid too small")
    end if

#if defined (__NEC__)
    x_(1:M,:) = x(1:M,:)
    x_t = transpose (x_(1:M,:))
#else
    x_t = transpose (x)
#endif

    n2 = ilog2 (M)
    nr = 1
    lvl_loop: do j = 1, n2
       if (nr >= Mmin) then
          s(:,1:nr) = x_t(:,1:nr     )
          d(:,1:nr) = x_t(:,nr+1:2*nr)
          call cdv_dyad_up_tr (s(:,1:nr), d(:,1:nr), x_t(:,1:2*nr))
       else
          if (cdv_keep_scaling) then
             nr = nr * 2
             cycle lvl_loop
          end if
          s(:,1:nr) = x_t(:,1:nr     )
          d(:,1:nr) = x_t(:,nr+1:2*nr)
          select case (nr)
          case (1)      ! 1 initial coarser-level coefficient: inverse Haar
             x_t(:,1)   = sqrt_one_half * (s(:,1) + d(:,1))
             x_t(:,2)   = sqrt_one_half * (s(:,1) - d(:,1))
          case (2)      ! 2 initial coarser-level coefficients: inverse join_2
             x_t(:,1:4) = matmul (s(:,1:2), cdv_lo42) &
                        + matmul (d(:,1:2), cdv_hi42)
          case (4)      ! 4 initial coarser-level coefficients: inverse join_4
             x_t(:,1:8) = matmul (s(:,1:4), cdv_lo84) &
                        + matmul (d(:,1:4), cdv_hi84)
          case default
             print *, "cdv_syn_vec: internal error trying to 'join':", nr, Mmin
             call finish ("cdv_syn_vec","internal error trying to 'join'")
          end select
       end if
       nr = nr * 2
    end do lvl_loop

    ! Inverse preconditioning (except for Haar transformation with order=1)
    if (cdv_precondition .and. cdv%N > 1) then
       call cdv_apply_post_tr (x_t)
    end if

#if defined (__NEC__)
    P = size (x, dim=2)
    x_t_(1:P,:) = x_t(1:P,:)
    x = transpose (x_t_(1:P,:))
#else
    x = transpose (x_t)
#endif

  end subroutine cdv_syn_vec

  ! --

  subroutine cdv_ana_ad_vec (x)
    ! Non-periodic adjoint wavelet analysis
    ! This variant should perform better on vector processors (e.g. SX-6)
    real(wp), intent(inout) :: x(:,:)

    ! Local variables
    integer  :: M, j, n2, Mmin, nr
    real(wp), dimension (size (x,dim=2), size (x,dim=1)  ) :: x_t
    real(wp), dimension (size (x,dim=2), size (x,dim=1)/2) :: s, d
    real(wp), parameter :: sqrt_one_half = 0.70710678118654752440_wp
#if defined (__NEC__)
    integer  :: P
    real(wp), dimension (size (x,dim=1)+1, size (x,dim=2)) :: x_
    real(wp), dimension (size (x,dim=2)+1, size (x,dim=1)) :: x_t_
#endif

    M = size (x, dim=1)

    if (.not. ispowerof2 (M)) then
       print *, "cdv_ana_ad_vec: grid size is not power of 2:", M
       call finish ("cdv_ana_ad_vec","bad grid size")
    end if
    Mmin = 2*cdv%N
    if (M < Mmin) then
       print *, "cdv_ana_ad_vec: grid too small for chosen interval wavelet:",M
       call finish ("cdv_ana_ad_vec","grid too small")
    end if

#if defined (__NEC__)
    x_(1:M,:) = x(1:M,:)
    x_t = transpose (x_(1:M,:))
#else
    x_t = transpose (x)
#endif

    n2 = ilog2 (M)
    nr = 1
    lvl_loop: do j = 1, n2
       if (nr >= Mmin) then
          s(:,1:nr) = x_t(:,1:nr     )
          d(:,1:nr) = x_t(:,nr+1:2*nr)
          call cdv_dyad_up_tr (s(:,1:nr), d(:,1:nr), x_t(:,1:2*nr))
       else
          if (cdv_keep_scaling) then
             nr = nr * 2
             cycle lvl_loop
          end if
          s(:,1:nr) = x_t(:,1:nr     )
          d(:,1:nr) = x_t(:,nr+1:2*nr)
          select case (nr)
          case (1)      ! 1 initial coarser-level coefficient: inverse Haar
             x_t(:,1)   = sqrt_one_half * (s(:,1) + d(:,1))
             x_t(:,2)   = sqrt_one_half * (s(:,1) - d(:,1))
          case (2)      ! 2 initial coarser-level coefficients: inverse join_2
             x_t(:,1:4) = matmul (s(:,1:2), cdv_lo42) &
                        + matmul (d(:,1:2), cdv_hi42)
          case (4)      ! 4 initial coarser-level coefficients: inverse join_4
             x_t(:,1:8) = matmul (s(:,1:4), cdv_lo84) &
                        + matmul (d(:,1:4), cdv_hi84)
          case default
             print *, "cdv_ana_ad: internal error trying to 'join':", nr, Mmin
             call finish ("cdv_ana_ad_vec","internal error trying to 'join'")
          end select
       end if
       nr = nr * 2
    end do lvl_loop

    ! Adjoint Preconditioning (not for Haar transformation)
    if (cdv_precondition .and. cdv%N > 1) then
       call cdv_apply_pre_ad_tr (x_t)
    end if

#if defined (__NEC__)
    P = size (x, dim=2)
    x_t_(1:P,:) = x_t(1:P,:)
    x = transpose (x_t_(1:P,:))
#else
    x = transpose (x_t)
#endif

  end subroutine cdv_ana_ad_vec

  ! --

  subroutine cdv_syn_ad_vec (x)
    ! Non-periodic adjoint wavelet synthesis
    ! This variant should perform better on vector processors (e.g. SX-6)
    real(wp), intent(inout) :: x(:,:)

    ! Local variables
    integer  :: M, j, n2, Mmin, nr
    real(wp), dimension (size (x,dim=2), size (x,dim=1)  ) :: x_t
    real(wp), dimension (size (x,dim=2), size (x,dim=1)/2) :: s, d
    real(wp), parameter :: sqrt_one_half = 0.70710678118654752440_wp
#if defined (__NEC__)
    integer  :: P
    real(wp), dimension (size (x,dim=1)+1, size (x,dim=2)) :: x_
    real(wp), dimension (size (x,dim=2)+1, size (x,dim=1)) :: x_t_
#endif

    M = size (x, dim=1)

    if (.not. ispowerof2 (M)) then
       print *, "cdv_syn_ad_vec: grid size is not power of 2:", M
       call finish ("cdv_syn_ad_vec","bad grid size")
    end if
    Mmin = 2*cdv%N
    if (M < Mmin) then
       print *, "cdv_syn_ad_vec: grid too small for chosen interval wavelet:",M
       call finish ("cdv_syn_ad_vec","grid too small")
    end if
#if defined (__NEC__)
    x_(1:M,:) = x(1:M,:)
    x_t = transpose (x_(1:M,:))
#else
    x_t = transpose (x)
#endif

    ! Adjoint Inverse Preconditioning (not for Haar transformation)
    if (cdv_precondition .and. cdv%N > 1) then
       call cdv_apply_post_ad_tr (x_t)
    end if

    n2 = ilog2 (M)
    nr = M
    lvl_loop: do j = 1, n2
       nr = nr / 2
       if (nr >= Mmin) then
          call cdv_dyad_down_tr (x_t(:,1:2*nr), s(:,1:nr), d(:,1:nr))
       else
          if (cdv_keep_scaling) exit lvl_loop
          select case (nr)
          case (1)      ! 1 remaining coarser-level coefficient: Haar
             d(:,1)   = sqrt_one_half * (x_t(:,1) - x_t(:,2))
             s(:,1)   = sqrt_one_half * (x_t(:,1) + x_t(:,2))
          case (2)      ! 2 remaining coarser-level coefficients: join 4->2
             d(:,1:2) = matmul (x_t(:,1:4), transpose (cdv_hi42))
             s(:,1:2) = matmul (x_t(:,1:4), transpose (cdv_lo42))
          case (4)      ! 4 remaining coarser-level coefficients: join 8->4
             d(:,1:4) = matmul (x_t(:,1:8), transpose (cdv_hi84))
             s(:,1:4) = matmul (x_t(:,1:8), transpose (cdv_lo84))
          case default
             print *, "cdv_syn_ad: internal error trying to 'join':", nr, Mmin
             call finish ("cdv_syn_ad_vec","internal error trying to 'join'")
          end select
       end if
       x_t(:,1:nr     ) = s(:,1:nr)
       x_t(:,nr+1:2*nr) = d(:,1:nr)
    end do lvl_loop

#if defined (__NEC__)
    P = size (x, dim=2)
    x_t_(1:P,:) = x_t(1:P,:)
    x = transpose (x_t_(1:P,:))
#else
    x = transpose (x_t)
#endif

  end subroutine cdv_syn_ad_vec

  !======================================================================

  function vanishing_moment_count () result (n)
    !------------------------------------------
    ! Determine the number of vanishing moments
    ! of a filter restricted to an interval
    !------------------------------------------
    integer :: n

    ! Local variables
    integer  :: i
    integer,  parameter :: M = 32, P = M/2
    real(wp) :: w(M,P), mom(P)
    real(wp), parameter :: tol = 1.e-14_wp

    ! Inlining bug with sxf90 revs.360,381
!CDIR NOIEXPAND
    call init_chebyshev (w)     ! Chebyshev polynomials of degrees 0...(P-1)
    call wave_1d (w, trans=TR_ANA)
    ! RMS of finest-level details
    mom(:) = sqrt (sum (w(M/2+1:M,:)**2, dim=1) / (M/2))
    ! Count the leading vanishing moments
    n = 0
    do i = 1, P
       if (mom(i) < tol) then
          n = n+1
       else
          exit
       end if
    end do
    !print '(a,i3)', "vanishing_moment_count: Number of vanishing moments:", n
  end function vanishing_moment_count

  ! --

  function filter_vanishing_moments () result (n)
    !-------------------------------------------------------------
    ! Analogous determination of the number of vanishing moments
    ! of filters on periodic intervals with filter length below 24.
    !-------------------------------------------------------------
    integer :: n

    ! Local variables
    integer  :: i
    integer,  parameter :: M = 64, P = M/2, M0 = 3*M/4, dM = 4
    real(wp) :: w(M,P), mom(P)
    real(wp), parameter :: tol = 2.e-12_wp  ! Not too tight (for symlets)!

    ! Inlining bug with sxf90 revs.360,381
!CDIR NOIEXPAND
    call init_chebyshev (w)     ! Chebyshev polynomials of degrees 0...(P-1)
    call wave_1d (w, trans=TR_ANA)
    ! Look at the RMS of some finest-level details far from both 'ends'
    mom(:) = sqrt (sum (w(M0-dM+1:M0+dM,:)**2, dim=1) / (2*dM))
    ! Count the leading vanishing moments
    n = 0
    do i = 1, P
       if (mom(i) < tol) then
          n = n+1
       else
          exit
       end if
    end do
  end function filter_vanishing_moments

  ! --

  subroutine init_chebyshev (T)
    ! Calculate leading Chebyshev polynomials over equispaced grid [-1,1]
    real(wp) :: T(:,0:)

    ! Local variables
    integer  :: M, J, i, n
    real(wp) :: x(size (T,dim=1))

    M = size (T, dim=1)         ! Number of grid points
    J = size (T, dim=2)         ! J-1 is maximum degree
    if (J == 0) return

    ! Generate equispaced grid over [-1,1]:
    x(:) = real (2* (/ (i, i=1,M) /) -M-1, wp) / (M-1)

    ! ChebyshevT(0,x) = 1
    ! ChebyshevT(1,x) = x
    ! ChebyshevT(2,x) = 2*x**2 - 1
    ! ChebyshevT(3,x) = x * (4*x**2 - 3)
    ! ChebyshevT(4,x) = 8*x**4 - 8*x**2 + 1
    ! ...
    ! Recurrence relation (n>=1): T(n+1,x) = 2*x*T(n,x) - T(n-1,x)
    T(:,0) = 1
    if (J == 1) return
    T(:,1) = x(:)
    do n = 1, J-2
       T(:,n+1) = 2*x(:)*T(:,n) - T(:,n-1)
    end do
  end subroutine init_chebyshev

  !======================================================================

  subroutine test_pwt_timing ()
    implicit none
    integer, parameter :: nx = 512      ! Dim.1 of test vector, don't change
    integer, parameter :: ny = 256      ! Dim.2 of test vector, don't change
    integer, parameter :: N_ITER = 16   ! Number of iterations
    integer, parameter :: N_PWT  = 3    ! Number of implementations

    integer :: old_basis, old_pwt_version
    integer :: tr_min, tr_max, n_tr
    integer :: i, j, k, trans
    real(wp)                   :: tmp, t0, t1
    real(wp), allocatable      :: cpu(:,:)
    real(wp), dimension(nx,ny) :: a2, b2, c2, d2    ! Arrays: 512*256

    integer, parameter :: f_trans_set(4) = &
         (/ TR_ANA, TR_SYN, TR_DUAL, TR_ADJ /)  ! "Forward" transformations
    integer, parameter :: i_trans_set(4) = &
         (/ TR_SYN, TR_ANA, TR_ADJ, TR_DUAL /)  ! "Inverse" transformations

    ! Save previous values
    old_basis       = basis
    old_pwt_version = ortho_pwt_version

    if (all (basis /= wv_ortho(:))) then
       basis=WV_DAUB_8
       print *
       print *, "test_pwt_timing: ", &
            "no suitable wavelet basis set, using default ", wv_name(basis)
    end if

    n_tr   = size   (f_trans_set)
    tr_min = minval (f_trans_set)
    tr_max = maxval (f_trans_set)

    allocate (cpu(tr_min:tr_max,N_PWT))
    cpu = 0
    call random_number (a2)

    print *
    print *, "Timings for wavelet basis: ", trim (wv_name (basis))
    print *
    do j = 1, N_PWT
       ortho_pwt_version = j
       print '(A,i2)', "ortho_pwt_version =", j
       do k = 1, n_tr
          call cpu_time (t0)
          trans = f_trans_set(k)
          print *, " 'Forward' = ", trim (tr_name(trans))
          do i = 1, N_ITER
             b2 = a2
             call wave_1d (b2, trans=trans)
          end do
          call cpu_time (t1)
          cpu(trans,j) = cpu(trans,j) + abs (t1 - t0)
          if (k == 1) then
             if (j == 1) then
                d2 = b2
             else
                tmp = sum ((b2-d2)**2) / ny
                if (tmp > 1.e-15) then
                   print '(A,g13.6)', "Check against reference FAILed:", tmp
                   call finish ("test_pwt_timing", "check program/compiler")
                end if
             end if
          end if
          call cpu_time (t0)
          trans = i_trans_set(k)
          print *, " 'Inverse' = ", trim (tr_name(trans))
          do i = 1, N_ITER
             c2 = b2
             call wave_1d (c2, trans=trans)
          end do
          call cpu_time (t1)
          cpu(trans,j) = cpu(trans,j) + abs (t1 - t0)
          tmp = sum ((c2-a2)**2) / ny
          if (tmp < 1.e-15) then
             print '(A,i0,A)', "Forward/Inverse (",k,") OK"
          else
             print '(A,i0,A,g13.6)', "Forward/Inverse (",k,") FAIL:", tmp
             call finish ("test_pwt_timing", "check program/compiler")
          end if
       end do
       print *
    end do
    print '(A,i0,A,i0,A,i0)', "Summary of CPU times: ", &
         2*N_ITER, " transformations, field size: ", nx, " * ", ny
    print *
    print '(A,1x,9A7)', "pwt_version:", &
         tr_name(f_trans_set( (/(k,k=1,n_tr)/) ))
    do j = 1, N_PWT
       print '(i6,7x,9f7.3)', j, cpu(f_trans_set( (/(k,k=1,n_tr)/) ),j)
    end do

    ! Restore previous values
    basis = old_basis
    ortho_pwt_version = old_pwt_version
  end subroutine test_pwt_timing

  !======================================================================

  subroutine test_cdv_timing ()
    implicit none
    integer, parameter :: nx = 512      ! Dim.1 of test vector, don't change
    integer, parameter :: ny = 256      ! Dim.2 of test vector, don't change
    integer, parameter :: N_ITER = 16   ! Number of iterations
    integer, parameter :: N_CDV  = 2    ! Number of implementations

    integer :: old_basis, old_cdv_version
    integer :: tr_min, tr_max, n_tr
    integer :: i, j, k, trans
    real(wp)                   :: tmp, t0, t1
    real(wp), allocatable      :: cpu(:,:)
    real(wp), dimension(ny,nx) :: a2, b2, c2    ! Arrays: 256*512

    integer, parameter :: f_trans_set(4) = &
         (/ TR_ANA, TR_SYN, TR_DUAL, TR_ADJ /)  ! "Forward" transformations
    integer, parameter :: i_trans_set(4) = &
         (/ TR_SYN, TR_ANA, TR_ADJ, TR_DUAL /)  ! "Inverse" transformations

    ! Save previous values
    old_basis       = cdv_basis
    old_cdv_version = cdv_version

    if (any (basis == wv_CDV_list)) then
       call cdv_setup (basis=basis)
    else
       call cdv_setup (basis=WV_CDV_8)
       print *
       print *, "test_cdv_timing: ", &
            "no suitable wavelet basis set, using default ", wv_name(cdv_basis)
    end if

    n_tr   = size   (f_trans_set)
    tr_min = minval (f_trans_set)
    tr_max = maxval (f_trans_set)

    allocate (cpu(tr_min:tr_max,N_CDV))
    cpu = 0
    call random_number (a2)

    print *
    print *, "Timings for wavelet basis: ", trim (wv_name (cdv_basis))
    print *
    do j = 1, N_CDV
       cdv_version = j
       print '(A,i2)', "cdv_version =", j
       do k = 1, n_tr
          call cpu_time (t0)
          trans = f_trans_set(k)
          print *, " 'Forward' = ", trim (tr_name(trans))
          do i = 1, N_ITER
             b2 = a2
             call wave_cdv (b2, trans=trans)
          end do
          call cpu_time (t1)
          cpu(trans,j) = cpu(trans,j) + abs (t1 - t0)
          call cpu_time (t0)
          trans = i_trans_set(k)
          print *, " 'Inverse' = ", trim (tr_name(trans))
          do i = 1, N_ITER
             c2 = b2
             call wave_cdv (c2, trans=trans)
          end do
          call cpu_time (t1)
          cpu(trans,j) = cpu(trans,j) + abs (t1 - t0)
          tmp = sum ((c2-a2)**2) / ny
          if (tmp < 1.e-15) then
             print '(A,i0,A)', "Forward/Inverse (",k,") OK"
          else
             print '(A,i0,A,g13.6)', "Forward/Inverse (",k,") FAIL:", tmp
             call finish ("test_cdv_timings", "check program/compiler")
          end if
       end do
       print *
    end do
    print '(A,i0,A,i0,A,i0)', "Summary of CPU times: ", &
         2*N_ITER, " transformations, field size: ", ny, " * ", nx
    print *
    print '(A,1x,9A7)', "cdv_version:", &
         tr_name(f_trans_set( (/(k,k=1,n_tr)/) ))
    do j = 1, N_CDV
       print '(i6,7x,9f7.3)', j, cpu(f_trans_set( (/(k,k=1,n_tr)/) ),j)
    end do

    ! Restore previous values
    if (any (old_basis == wv_CDV_list)) then
       call cdv_setup (basis=old_basis)
    end if
    cdv_version = old_cdv_version
  end subroutine test_cdv_timing

  !======================================================================

end module mo_1dmra
