!
!+ radiances specific observation operator derived type definition
!
MODULE mo_t_tovs
!
! Description:
!   radiances specific observation operator derived type definition
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

  !=============
  ! modules used
  !=============
  use mo_kind,        only: wp, sp, i8  ! kind parameters
  use mo_exception,   only: finish      ! abort routine
  use mo_t_obs,       only: t_obs,     &! observation derived type
                            t_spot,    &! report/fof derived type
                            new_par     ! reserve memory
  use mo_rad,         only: t_rad_set, &! specific. of instr. & channels in FOV
                            destruct,  &! t_rad_set destructor routine
                            n_set,     &
                            rad_set,   &
                            m_instr,   &
                            m_chan,    &
                            set_indx,  &
                            n_styp,    &
                            print_rad_set
  use mo_t_obs,       only: t_obs,     &
                            t_spot

  implicit none

  !================
  ! public entities
  !================
  private
  public :: t_tovs         ! observation operator specific type
  public :: t_tovs_instr   ! info on instruments in t_tovs
  public :: set_size       ! set '..._size' module variables
  public :: construct      ! t_tovs constructor routine
  public :: destruct       ! t_tovs  destructor routine
  public :: load           ! load  t_tovs
  public :: store          ! store t_tovs
  public :: get_tovs_rs    ! get rad_set corresponding to t_tovs
  public :: get_tovs_instr ! get rad_set corresponding to t_tovs
  public :: link_tovs_rs   ! link existing t_tovs with the corresponding global rad_set
  public :: add_av_cont    ! Add entry to av_cont
  public :: get_im_ch_ind
  public :: TTOVS_BASE_BIT
  public :: TTOVS_CI_BIT
  public :: TTOVS_CI_EXT_BIT
  public :: TTOVS_AV_BIT
  public :: TTOVS_CEMI_BIT
  public :: TTOVS_L2C_BIT
  public :: TTOVS_EMIS_BIT
  public :: TTOVS_FLAG_BIT
  public :: TTOVS_SINFL_BIT
  public :: TTOVS_SPEC_BIT
  public :: TTOVS_BTCS_BIT
  public :: TTOVS_TR_BIT
  public :: TTOVS_IMCL_BIT
  public :: TTOVS_IMCH_BIT
  public :: TTOVS_CLDLEV_BIT
  public :: TTOVS_BASE     ! store basic entries in t_tovs
  public :: TTOVS_CI       ! store/load t_tovs%ci
  public :: TTOVS_AV       ! store/load t_tovs%av
  public :: TTOVS_CEMI     ! store/load t_tovs%cemi
  public :: TTOVS_L2C      ! store/load t_tovs%l2c
  public :: TTOVS_EMIS     ! store/load t_tovs%emis
  public :: TTOVS_BTCS     ! store/load t_tovs%bt_cs
  public :: TTOVS_FLAG
  public :: TTOVS_SINFL
  public :: TTOVS_SPEC
  public :: TTOVS_TR
  public :: TTOVS_IMCL
  public :: TTOVS_IMCH
  public :: TTOVS_ALL
  public :: TTOVS_IM
  public :: TTOVS_CLDLEV
  public :: IMCL_FRAC
  public :: IMCL_MEAN
  public :: IMCL_STDV
  public :: IMCH_MEAN_O
  public :: IMCH_MEAN_B
  public :: IMCH_STDV
  public :: IMCH_EMIS
  public :: mx_nav
  public :: mx_nlev
  public :: mx_imch
  public :: tpp

  !-----------
  ! parameters
  !-----------
  integer            :: j

  integer, parameter :: tpp = sp              ! tovs profile precision (t_tovs%av)
!  integer, parameter :: tpp = wp              ! tovs profile precision (t_tovs%av)


  integer, parameter :: TTOVS_BASE_BIT   = 0  ! Basic info in t_tovs set
  integer, parameter :: TTOVS_CI_BIT     = 1  ! ci   array in t_tovs allocated and set
  integer, parameter :: TTOVS_CI_EXT_BIT = 2  ! ci   array in t_tovs set points to external array
  integer, parameter :: TTOVS_AV_BIT     = 3  ! av   array in t_tovs allocated and set
  integer, parameter :: TTOVS_CEMI_BIT   = 4  ! cemi array in t_tovs allocated and set
  integer, parameter :: TTOVS_L2C_BIT    = 5  ! l2c  array in t_tovs allocated and set
  integer, parameter :: TTOVS_EMIS_BIT   = 6  ! emis array in t_tovs allocated and set
  integer, parameter :: TTOVS_FLAG_BIT   = 7  ! flag array in t_tovs allocated and set
  integer, parameter :: TTOVS_SINFL_BIT  = 8  ! sinfl array in t_tovs allocated and set
  integer, parameter :: TTOVS_SPEC_BIT   = 9  ! spec array in t_tovs allocated and set
  integer, parameter :: TTOVS_BTCS_BIT   = 10 ! bt_cs array in t_tovs allocated and set
  integer, parameter :: TTOVS_TR_BIT     = 11 ! tr array in t_tovs allocated and set
  integer, parameter :: TTOVS_IMCL_BIT   = 12 ! imager info array in t_tovs allocated and set
  integer, parameter :: TTOVS_IMCH_BIT   = 13 ! imager info array in t_tovs allocated and set
  integer, parameter :: TTOVS_CLDLEV_BIT = 14 ! cldlev array in t_tovs allocated and set

  integer, parameter :: TTOVS_BASE   = 2**TTOVS_BASE_BIT    ! Basic info in t_tovs set
  integer, parameter :: TTOVS_CI     = 2**TTOVS_CI_BIT      ! ci   array in t_tovs allocated and set
  integer, parameter :: TTOVS_CI_EXT = 2**TTOVS_CI_EXT_BIT  ! ci   array in t_tovs set points to external array
  integer, parameter :: TTOVS_AV     = 2**TTOVS_AV_BIT      ! av   array in t_tovs allocated and set
  integer, parameter :: TTOVS_CEMI   = 2**TTOVS_CEMI_BIT    ! cemi array in t_tovs allocated and set
  integer, parameter :: TTOVS_L2C    = 2**TTOVS_L2C_BIT     ! l2c  array in t_tovs allocated and set
  integer, parameter :: TTOVS_EMIS   = 2**TTOVS_EMIS_BIT    ! emis array in t_tovs allocated and set
  integer, parameter :: TTOVS_FLAG   = 2**TTOVS_FLAG_BIT    ! flag array in t_tovs allocated and set
  integer, parameter :: TTOVS_SINFL  = 2**TTOVS_SINFL_BIT   ! sinfl array in t_tovs allocated and set
  integer, parameter :: TTOVS_SPEC   = 2**TTOVS_SPEC_BIT    ! spec array in t_tovs allocated and set
  integer, parameter :: TTOVS_BTCS   = 2**TTOVS_BTCS_BIT    ! btcs array in t_tovs allocated and set
  integer, parameter :: TTOVS_TR     = 2**TTOVS_TR_BIT      ! tr   array in t_tovs allocated and set
  integer, parameter :: TTOVS_IMCL   = 2**TTOVS_IMCL_BIT    ! im_cl array in t_tovs allocated and set
  integer, parameter :: TTOVS_IMCH   = 2**TTOVS_IMCH_BIT    ! im_ch array in t_tovs allocated and set
  integer, parameter :: TTOVS_CLDLEV = 2**TTOVS_CLDLEV_BIT  ! l2c  array in t_tovs allocated and set

!  The intel compiler does not accept this:
!  integer, parameter :: TTOVS_ALL  = sum( (/ (2**j,j=0,bit_size(TTOVS_BASE)-1) /) )
!  Workaround:
  integer, parameter :: TTOVS_ALL  = TTOVS_BASE + TTOVS_CI   + TTOVS_AV   + &
                                     TTOVS_CEMI + TTOVS_L2C  + TTOVS_EMIS + &
                                     TTOVS_FLAG + TTOVS_SINFL+ TTOVS_SPEC + &
                                     TTOVS_BTCS + TTOVS_TR   + TTOVS_IMCL + &
                                     TTOVS_IMCH + TTOVS_CLDLEV
  integer, parameter :: TTOVS_IM   = TTOVS_IMCL + TTOVS_IMCH

  integer, parameter :: mx_av = 15

  ! Flags for im_*_v
  integer, parameter :: IMCL_FRAC   = 1
  integer, parameter :: IMCL_MEAN   = 2
  integer, parameter :: IMCL_STDV   = 3
  integer, parameter :: mx_imcl     = 2 * 2 + 1  ! (mean,stdv) for 2 channels + cluster fraction
  integer, parameter :: IMCH_MEAN_O = 1
  integer, parameter :: IMCH_MEAN_B = 2
  integer, parameter :: IMCH_STDV   = 3
  integer, parameter :: IMCH_EMIS   = 4
  integer, parameter :: mx_imch     = 4

  !-------------------------
  ! private module variables
  !-------------------------
  integer ,save      :: tovs_size     = 0       ! size of type t_tovs
  integer ,save      ::  int_size     = 0       ! size of integer
  integer ,save      ::   sp_size     = 0       ! size of real(sp)
  integer ,save      ::   wp_size     = 0       ! size of real(sp)
  integer ,save      ::  tpp_size     = 0       ! size of real(sp)

  integer ,save      :: mx_nav        = 0       ! highest value of t_tovs%nav
  integer ,save      :: mx_nlev       = 0       ! highest valud of t_tovs%nlev

  !---------------
  ! tovs data type
  !---------------
  type t_tovs
    !--------------------------
    ! TOVS specific information
    !--------------------------
    integer           :: id_rs     =  -1      ! id of corresponding (global) rad_set
    integer           :: nchan                ! number of channels
    integer           :: init      =  0       ! bit field with TTOVS_* bit set for
                                              ! initialized entries
    integer,  pointer :: ci(:)     => NULL()  ! indices of channels in corresponding (global) t_rad_set
    integer           :: nlev      =  0       ! number of levels of bg profiles
    integer           :: scanl     = -1       ! scanline number
    ! OBSOLETE: is part of t_spot
    real(sp)          :: saza      = -999._sp ! sat.  zenith  angle [degree]
    real(sp)          :: boa       = -999._sp ! sat.  azimuth angle [degree]
    ! END OBSOLETE
    real(sp)          :: soza      = -999._sp ! solar zenith  angle [degree]
    real(sp)          :: soa       = -999._sp ! solar azimuth angle [degree]
    real(sp)          :: sasoa     = -999._sp ! sat-sol. azim.angle [degree]
    integer           :: mw_stype  = -99      ! surface type (derived from MW obs)
    integer           :: rt_stype(1:n_styp) = -99 ! surface type (derived from MW obs)
    real(sp)          :: w_styp  (1:n_styp) = 1._sp
    integer           :: ns        =  1       ! number of surface types in spot
    real(sp)          :: land_frac =  1._sp   ! land fraction (for FOV, calculated in sat_pp)
                                              ! TODO:obsolete? ,unify with w_styp
    real(sp)          :: snf       =  0._sp   ! snow fraction
    real(sp)          :: cloud_imag= -999._sp ! cloud cover of colocated imager
    real(sp)          :: cloud_top = -999._sp ! cloud top
    real(sp)          :: e_bg_ts   =  1._sp
    integer           :: nwc_flg   =  0       ! NWC_FLAG from sat_pp files
    real(sp)          :: orb_phase = -999._sp ! orbit phase
    real(sp)          :: instr_temp= -999._sp ! instrument temp.
    real(sp)          :: lwp_qin   = -999._sp ! lwp_qinzou2016 (requires all channels and fg and cannot be
                                              ! calculated everywhere, therefore stored here)
    ! level numbers
    integer           :: nld_t     =  -1      ! T dummy variable start level
    integer           :: nld_h     =  -1      ! Q dummy variable start level
    integer           :: nl_st     =  -1      ! Number of levels filled from file_strato
    ! description of av (background profiles)
    integer           :: nav       =   0
    integer           :: i_t       =  -1
    integer           :: i_q       =  -1
    integer           :: i_p       =  -1
    integer           :: i_o3      =  -1
    integer           :: i_co2     =  -1
    integer,  pointer :: flag (:)  => NULL()  ! radiance specific flags
    real(tpp),pointer :: av (:,:)  => NULL()  ! Background profile: atm.variables (lev,t:q:t_a:q_a)
    integer(i8)       :: av_cont(mx_av) =0_i8 ! contents of av, identifier from mo_t_col

    ! emissivity retrieval
    integer           :: ncemi     =  0
    real(sp), pointer :: cemi (:)  => NULL()  ! emissivity correction factor
    ! level to channel assignment
    integer           :: nl2c      =  0
    real(sp), pointer :: l2c  (:)  => NULL()  ! level to channel assignment (chan)
    ! TODO: use sp instead of wp
    real(wp), pointer :: emis (:)  => NULL()  ! emissivity (chan)
    real(sp), pointer :: sinfl(:)  => NULL()
    ! specularity
    integer           :: nspec     =  0
    real(sp), pointer :: spec (:)  => NULL()
    ! clear-sky brightness temp
    integer           :: nbtcs     =  0
    real(sp), pointer :: bt_cs (:) => NULL()  ! clear-sky TB (chan)
    ! transmission factors for biascorrection
    integer           :: ntr       =  0
    real(sp), pointer :: tr (:,:)  => NULL()  ! transmission decrease between plevels
    ! background surface variable
    real(wp)          :: ts        =  -999._wp ! skin temperature
    logical           :: l_ts_dr   =  .false.
    ! colocated imager info
    integer           :: n_im_cl   =  0       ! number of imager clusters
    integer           :: n_im_cl_v =  0       ! number of imager cluster values (2*nchan+1, see mx_im above)
    integer           :: im_cl_v(mx_imcl) = 0 ! imager cluster value "flags": 100*chan + IM_* (see IM_* above)
    real(sp), pointer :: im_cl(:,:)=> NULL()  ! imager cluster values (n_im_cl_v, n_im_cl)
    integer           :: n_im_ch   =  0       ! number imager channels
    integer           :: n_im_ch_v =  0       ! number of imager cluster values (2*nchan+1, see mx_im above)
    integer           :: im_ch_v(mx_imch)=  0 ! imager value "flags": (see IM_* above)
    real(sp), pointer :: im_ch(:,:)=> NULL()  ! imager channel values (n_im_v, n_im_ch_ch)
    integer           :: cld_flg   =  0       ! Spot cloud flag for pixel, e.g. imager cloud flag from MNW
    ! cloud_lev
    integer           :: nband     =  0
    real(sp), pointer :: cldlev(:) => NULL()  ! cloud level of band
  end type t_tovs

  !--------------------------------------------------------------
  ! auxiliary data type with information on instruments in t_tovs
  !--------------------------------------------------------------
  ! This is not part of t_tovs, since t_tovs should be as small as possible
  type t_tovs_instr
    integer           :: n_instr         = 0  ! number of instruments
    integer           :: ii(m_instr)     = 0  ! index of instrument in global t_rad_set
    integer           :: o_ch_i(m_instr) = 0  ! channel offset for instrument in t_tovs arrays
    integer           :: n_ch_i(m_instr) = 0  ! number of channels for instrument in t_tovs
  end type t_tovs_instr

  !-----------
  ! interfaces
  !-----------
  interface load
    module procedure load_tovs
  end interface load

  interface store
    module procedure store_tovs
  end interface store

  interface destruct
    module procedure destruct_tovs
  end interface destruct

  interface construct
    module procedure construct_tovs
    module procedure construct_tovs_instr
  end interface construct

!==============================================================================
contains
!==============================================================================

#if defined (_FTRACE)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif

  subroutine set_size
    !------------------------------------
    ! determine sizes of local data types
    !------------------------------------
    type (t_obs)  :: obs
    type (t_tovs) :: tovs
    tovs_size = size (transfer (tovs  ,obs% par))
    int_size  = size (transfer (1     ,obs% par))
    sp_size   = size (transfer (1._sp ,obs% par))
    wp_size   = size (transfer (1._wp ,obs% par))
    tpp_size  = size (transfer (1._tpp,obs% par))
  end subroutine set_size
!------------------------------------------------------------------------------
  subroutine load_tovs (obs, spot, tovs, rs, i_rs, ti, ci, av, cemi, l2c, emis, flag, &
                        sinfl, spec, bt_cs, tr, im_cl, im_ch, cldlev, tovs_io)
    type (t_obs),             intent(in)            :: obs       ! data of all observations
    type (t_spot),            intent(in)            :: spot      ! meta data of this observation
    type (t_tovs),   target,  intent(out), optional :: tovs      ! TOVS specific information
    type(t_rad_set), pointer, intent(out), optional :: rs
    integer,                  intent(out), optional :: i_rs
    type(t_tovs_instr),       intent(out), optional :: ti
    integer,         target,  intent(out), optional :: ci   (:)
    real(tpp),                intent(out), optional :: av (:,:)
    real(wp),                 intent(out), optional :: cemi (:)
    real(sp),                 intent(out), optional :: l2c  (:)
    real(wp),                 intent(out), optional :: emis (:)
    integer,                  intent(out), optional :: flag (:)
    real(sp),                 intent(out), optional :: sinfl(:)
    real(sp),                 intent(out), optional :: spec (:)
    real(sp),                 intent(out), optional :: bt_cs(:)
    real(sp),                 intent(out), optional :: tr (:,:)
    real(sp),                 intent(out), optional :: im_cl(:,:)
    real(sp),                 intent(out), optional :: im_ch(:,:)
    real(sp),                 intent(out), optional :: cldlev(:)
    integer,                  intent(in),  optional :: tovs_io     ! return arrays in t_tovs?

    type(t_tovs), target  :: tovs_dum
    type(t_tovs), pointer :: tovs_p
    integer  :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14
    integer  :: nchan, nlev, nl2c, ncemi, nav, nspec, nbtcs, ntr, nband
    integer  :: n_im_cl, n_im_cl_v, n_im_ch, n_im_ch_v
    integer  :: tovs_io_ ! specify, which arrays in tovs are filled
    real(sp) :: rsp  (1)  !+++ work around bug in NEC sxf90/2.0rev360-ia64
    real(tpp):: rtp  (1)
    real(wp) :: rwp  (1)
    integer  :: iwp  (1)

    FTRACE_BEGIN('load_tovs:init')
    if (tovs_size == 0) call set_size  ! if store_tovs was called on another PE

    n1    = tovs_size                  ! -> tovs

    if (spot%p%n == 0) call finish ('load_tovs','spot%p%n == 0')
    if (spot%p%n < n1) call finish ('load_tovs','spot%p%n < n1')

    if (present(tovs)) then
      tovs_p => tovs
    else
      tovs_p => tovs_dum
    end if
    if (present(tovs_io)) then
      tovs_io_ = tovs_io
    else
      tovs_io_ = TTOVS_ALL
    end if

    tovs_p  = transfer (obs% par (spot%p%i+1 : spot%p%i+n1), tovs_p)
    tovs_p%init = TTOVS_BASE

    nchan     = max(tovs_p% nchan    ,0)
    nlev      = max(tovs_p% nlev     ,0)
    nl2c      = max(tovs_p% nl2c     ,0)
    ncemi     = max(tovs_p% ncemi    ,0)
    nav       = max(tovs_p% nav      ,0)
    nspec     = max(tovs_p% nspec    ,0)
    nbtcs     = max(tovs_p% nbtcs    ,0)
    ntr       = max(tovs_p% ntr      ,0)
    n_im_cl   = max(tovs_p% n_im_cl  ,0)
    n_im_cl_v = max(tovs_p% n_im_cl_v,0)
    n_im_ch   = max(tovs_p% n_im_ch  ,0)
    n_im_ch_v = max(tovs_p% n_im_ch_v,0)
    nband     = max(tovs_p% nband    ,0)

    n2  = int_size * nchan             + n1  ! -> tovs% ci
    n3  = tpp_size * nav * nlev        + n2  ! -> tovs% av
    n4  =  sp_size * ncemi             + n3  ! -> tovs% cemi
    n5  =  sp_size * nl2c              + n4  ! -> tovs% l2c
    n6  =  wp_size * nchan             + n5  ! -> tovs% emis
    n7  = int_size * nchan             + n6  ! -> tovs% flag
    n8  =  sp_size * nchan             + n7  ! -> tovs% sinfl
    n9  =  sp_size * nspec             + n8  ! -> tovs% spec
    n10 =  sp_size * nbtcs             + n9  ! -> tovs% bt_cs
    n11 =  sp_size * ntr * nchan       + n10 ! -> tovs% tr
    n12 =  sp_size * n_im_cl*n_im_cl_v + n11 ! -> tovs% im_cl
    n13 =  sp_size * n_im_ch*n_im_ch_v + n12 ! -> tovs% im_ch
    n14 =  sp_size * nband             + n13 ! -> tovs% cldlev

    if (spot%p%n /= n14) then
      write(0,*) 'spot%hd%id,spot%p%n,n13,tovs_io_',spot%hd%id,spot%p%n,n14, tovs_io_
      call finish ('load_tovs','spot%p%n / n')
    end if
    FTRACE_END('load_tovs:init')

    FTRACE_BEGIN('load_tovs:ci')
    if (associated(tovs_p%ci)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_CI)>0) then
        allocate (tovs_p% ci  (nchan))
        tovs_p% ci    =          transfer (obs% par (spot%p%i+1+n1 : spot%p%i+n2), iwp)
        tovs_p% init  = tovs_p% init + TTOVS_CI
      end if
      if (present(ci)) then
        if (size(ci) < nchan) call finish('load_tovs','array "ci" too small')
        ci(1:nchan)   =          transfer (obs% par (spot%p%i+1+n1 : spot%p%i+n2), iwp)
        if (present(ti).and..not.btest(tovs_p%init,TTOVS_CI)) then
          tovs_p% ci => ci ! Required in get_tovs_instr/get_tovs_rs
          tovs_p% init = tovs_p% init + TTOVS_CI_EXT
        end if
      end if
    end if
    FTRACE_END('load_tovs:ci')
    FTRACE_BEGIN('load_tovs:av')
    if (associated(tovs_p%av)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_AV)>0) then
        allocate (tovs_p% av  (nlev,nav))
        tovs_p% av    = reshape (transfer (obs% par (spot%p%i+1+n2 : spot%p%i+n3), rtp),(/nlev,nav/))
        tovs_p% init  = tovs_p% init + TTOVS_AV
      end if
      if (present(av)) then
        if (size(av,1) < nlev .or. size(av,2) < nav) call finish('load_tovs','array "av" too small')
        av(1:nlev,1:nav) = real(reshape (transfer (obs% par (spot%p%i+1+n2 : spot%p%i+n3), rtp),(/nlev,nav/)), kind=wp)
      end if
    end if
    FTRACE_END('load_tovs:av')
    FTRACE_BEGIN('load_tovs:cemi')
    if (associated(tovs_p%cemi)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_CEMI)>0) then
        allocate (tovs_p% cemi(ncemi))
        tovs_p% cemi  =          transfer (obs% par (spot%p%i+1+n3 : spot%p%i+n4), rsp)
        tovs_p% init  = tovs_p% init + TTOVS_CEMI
      end if
      if (present(cemi)) then
        if (size(cemi) < ncemi) call finish('load_tovs','array "cemi" too small')
        cemi(1:ncemi) =          transfer (obs% par (spot%p%i+1+n3 : spot%p%i+n4), rsp)
      end if
    end if
    FTRACE_END('load_tovs:cemi')
    FTRACE_BEGIN('load_tovs:l2c')
    if (associated(tovs_p%l2c)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_L2C)>0) then
        allocate (tovs_p% l2c (nl2c))
        tovs_p% l2c   =          transfer (obs% par (spot%p%i+1+n4 : spot%p%i+n5), rsp)
        tovs_p% init  = tovs_p% init + TTOVS_L2C
      end if
      if (present(l2c)) then
        if (size(l2c) < nl2c) call finish('load_tovs','array "l2c" too small')
        l2c(1:nl2c)   =          transfer (obs% par (spot%p%i+1+n4 : spot%p%i+n5), rsp)
      end if
    end if
    FTRACE_END('load_tovs:l2c')
    FTRACE_BEGIN('load_tovs:emis')
    if (associated(tovs_p%emis)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_EMIS)>0) then
        allocate (tovs_p% emis(nchan))
        tovs_p% emis  =          transfer (obs% par (spot%p%i+1+n5 : spot%p%i+n6), rwp)
        tovs_p% init  = tovs_p% init + TTOVS_EMIS
      end if
      if (present(emis)) then
        if (size(emis) < nchan) call finish('load_tovs','array "emis" too small')
        emis(1:nchan) =          transfer (obs% par (spot%p%i+1+n5 : spot%p%i+n6), rwp)
      end if
    end if
    FTRACE_END('load_tovs:emis')
    FTRACE_BEGIN('load_tovs:flag')
    if (associated(tovs_p%flag)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_FLAG)>0) then
        allocate (tovs_p% flag(nchan))
        tovs_p% flag  =          transfer (obs% par (spot%p%i+1+n6 : spot%p%i+n7), iwp)
        tovs_p% init  = tovs_p% init + TTOVS_FLAG
      end if
      if (present(flag)) then
        if (size(flag) < nchan) call finish('load_tovs','array "flag" too small')
        flag(1:nchan) =          transfer (obs% par (spot%p%i+1+n6 : spot%p%i+n7), iwp)
      end if
    end if
    FTRACE_END('load_tovs:flag')
    FTRACE_BEGIN('load_tovs:sinfl')
    if (associated(tovs_p%sinfl)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_SINFL)>0) then
        allocate (tovs_p% sinfl(nchan))
        tovs_p% sinfl  =          transfer (obs% par (spot%p%i+1+n7 : spot%p%i+n8), rsp)
        tovs_p% init  = tovs_p% init + TTOVS_SINFL
      end if
      if (present(sinfl)) then
        if (size(sinfl) < nchan) call finish('load_tovs','array "sinfl" too small')
        sinfl(1:nchan) =          transfer (obs% par (spot%p%i+1+n7 : spot%p%i+n8), rsp)
      end if
    end if
    FTRACE_END('load_tovs:sinfl')
    FTRACE_BEGIN('load_tovs:spec')
    if (associated(tovs_p%spec)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_SPEC)>0) then
        allocate (tovs_p% spec(nspec))
        tovs_p% spec  =          transfer (obs% par (spot%p%i+1+n8 : spot%p%i+n9), rsp)
        tovs_p% init  = tovs_p% init + TTOVS_SPEC
      end if
      if (present(spec)) then
        if (size(spec) < nspec) call finish('load_tovs','array "spec" too small')
        spec(1:nspec) =          transfer (obs% par (spot%p%i+1+n8 : spot%p%i+n9), rsp)
      end if
    end if
    FTRACE_END('load_tovs:spec')
    FTRACE_BEGIN('load_tovs:bt_cs')
    if (associated(tovs_p%bt_cs)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_BTCS)>0) then
        allocate (tovs_p% bt_cs(nbtcs))
        tovs_p% bt_cs  =          transfer (obs% par (spot%p%i+1+n9 : spot%p%i+n10), rsp)
        tovs_p% init  = tovs_p% init + TTOVS_BTCS
      end if
      if (present(bt_cs)) then
        if (size(bt_cs) < nbtcs) call finish('load_tovs','array "bt_cs" too small')
        bt_cs(1:nbtcs) =          transfer (obs% par (spot%p%i+1+n9 : spot%p%i+n10), rsp)
      end if
    end if
    FTRACE_END('load_tovs:bt_cs')
    FTRACE_BEGIN('load_tovs:tr')
    if (associated(tovs_p%tr)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_TR)>0) then
        allocate (tovs_p% tr(nchan,ntr))
        tovs_p% tr  =  reshape (  transfer (obs% par (spot%p%i+1+n10 : spot%p%i+n11), rsp), (/nchan,ntr/))
        tovs_p% init  = tovs_p% init + TTOVS_TR
      end if
      if (present(tr)) then
        if (size(tr) < ntr*nchan) call finish('load_tovs','array "tr" too small')
        tr(1:nchan,1:ntr)= reshape(transfer(obs% par (spot%p%i+1+n10 : spot%p%i+n11), rsp), (/nchan,ntr/))
      end if
    end if
    FTRACE_END('load_tovs:tr')
    FTRACE_BEGIN('load_tovs:im_cl')
    if (associated(tovs_p%im_cl)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_IMCL)>0) then
        allocate (tovs_p% im_cl(n_im_cl_v,n_im_cl))
        tovs_p% im_cl =  reshape( transfer (obs% par (spot%p%i+1+n11 : spot%p%i+n12), rsp), (/n_im_cl_v,n_im_cl/))
        tovs_p% init  = tovs_p% init + TTOVS_IMCL
      end if
      if (present(im_cl)) then
        if (size(im_cl) < n_im_cl_v*n_im_cl) call finish('load_tovs','array "im_cl" too small')
        im_cl(1:n_im_cl_v,1:n_im_cl)= reshape(transfer(obs% par (spot%p%i+1+n11 : spot%p%i+n12), rsp), (/n_im_cl_v,n_im_cl/))
      end if
    end if
    FTRACE_END('load_tovs:im_cl')
    FTRACE_BEGIN('load_tovs:im_ch')
    if (associated(tovs_p%im_ch)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_IMCH)>0) then
        allocate (tovs_p% im_ch(n_im_ch_v,n_im_ch))
        tovs_p% im_ch =  reshape( transfer (obs% par (spot%p%i+1+n12 : spot%p%i+n13), rsp), (/n_im_ch_v,n_im_ch/))
        tovs_p% init  = tovs_p% init + TTOVS_IMCH
      end if
      if (present(im_ch)) then
        if (size(im_ch) < n_im_ch_v*n_im_ch) call finish('load_tovs','array "im_ch" too small')
        im_ch(1:n_im_ch_v,1:n_im_ch)= reshape(transfer(obs% par (spot%p%i+1+n12 : spot%p%i+n13), rsp), (/n_im_ch_v,n_im_ch/))
      end if
    end if
    FTRACE_END('load_tovs:im_*')
    FTRACE_BEGIN('load_tovs:cldlev')
    if (associated(tovs_p%cldlev)) then
      if (present(tovs).and.iand(tovs_io_,TTOVS_CLDLEV)>0) then
        allocate (tovs_p% cldlev(nband))
        tovs_p% cldlev =          transfer (obs% par (spot%p%i+1+n13 : spot%p%i+n14), rsp)
        tovs_p% init   = tovs_p% init + TTOVS_CLDLEV
      end if
      if (present(cldlev)) then
        if (size(cldlev) < nband) call finish('load_tovs','array "cldlev" too small')
        cldlev(1:nband) =          transfer (obs% par (spot%p%i+1+n13 : spot%p%i+n14), rsp)
      end if
    end if
    FTRACE_END('load_tovs:cldlev')

    FTRACE_BEGIN('load_tovs:rs/ti')
    if (present(rs) .or. present(i_rs) .or. present(ti)) then
      call get_tovs_rs(tovs_p, i_rs=i_rs, rs=rs, ti=ti)
    end if
    FTRACE_END('load_tovs:rs/ti')

  end subroutine load_tovs
!------------------------------------------------------------------------------
  subroutine store_tovs (obs, spot, tovs, ci, av, cemi, l2c, emis, flag, sinfl, spec, &
                         bt_cs, tr, im_cl, im_ch, cldlev, tovs_io)
    type (t_obs),  intent(inout)         :: obs       ! data of all observations
    type (t_spot), intent(inout)         :: spot      ! meta data of this observation
    type (t_tovs), intent(inout)         :: tovs      ! TOVS specific data type
    integer,       intent(in),  optional :: ci   (:)
    real(tpp),     intent(in),  optional :: av (:,:)
    real(wp),      intent(in),  optional :: cemi (:)
    real(sp),      intent(in),  optional :: l2c  (:)
    real(wp),      intent(in),  optional :: emis (:)
    integer,       intent(in),  optional :: flag (:)
    real(sp),      intent(in),  optional :: sinfl(:)
    real(sp),      intent(in),  optional :: spec (:)
    real(sp),      intent(in),  optional :: bt_cs(:)
    real(sp),      intent(in),  optional :: tr   (:,:)
    real(sp),      intent(in),  optional :: im_cl(:,:)
    real(sp),      intent(in),  optional :: im_ch(:,:)
    real(sp),      intent(out), optional :: cldlev(:)
    integer,       intent(in),  optional :: tovs_io  ! specify, which arrays from tovs are stored

    integer :: nlev, nchan, nl2c, ncemi, nav, nspec, nbtcs, ntr, nband
    integer :: n_im_cl, n_im_cl_v, n_im_ch, n_im_ch_v
    integer :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14
    integer :: tovs_io_
    integer :: par(1)  !+++ work around bug in NEC sxf90/2.0rev360-ia64

    nchan =     max(tovs% nchan    ,0)
    nlev  =     max(tovs% nlev     ,0)
    nl2c  =     max(tovs% nl2c     ,0)
    ncemi =     max(tovs% ncemi    ,0)
    nav   =     max(tovs% nav      ,0)
    nspec =     max(tovs% nspec    ,0)
    nbtcs =     max(tovs% nbtcs    ,0)
    ntr       = max(tovs% ntr      ,0)
    n_im_cl   = max(tovs% n_im_cl  ,0)
    n_im_cl_v = max(tovs% n_im_cl_v,0)
    n_im_ch   = max(tovs% n_im_ch  ,0)
    n_im_ch_v = max(tovs% n_im_ch_v,0)
    nband     = max(tovs% nband    ,0)

    mx_nav  = max(mx_nav,  nav)
    mx_nlev = max(mx_nlev, nlev)

    ! if (nlev      <= 0) call finish ('store_tovs','nlev<=0 (rttov not initialised)')
    if (tovs_size == 0) call set_size

    if (present(tovs_io)) then
      tovs_io_ = tovs_io
    else
      tovs_io_ = TTOVS_ALL
    end if

    n1  = tovs_size                           ! <- tovs
    n2  =  int_size * nchan             + n1  ! <- tovs% ci
    n3  =  tpp_size * nav * nlev        + n2  ! <- tovs% av
    n4  =   sp_size * ncemi             + n3  ! <- tovs% cemi
    n5  =   sp_size * nl2c              + n4  ! <- tovs% l2c
    n6  =   wp_size * nchan             + n5  ! <- tovs% emis
    n7  =  int_size * nchan             + n6  ! <- tovs% flag
    n8  =   sp_size * nchan             + n7  ! <- tovs% sinfl
    n9  =   sp_size * nspec             + n8  ! -> tovs% spec
    n10 =   sp_size * nbtcs             + n9  ! <- tovs% bt_cs
    n11 =   sp_size * ntr * nchan       + n10 ! -> tovs% tr
    n12 =   sp_size * n_im_cl*n_im_cl_v + n11 ! -> tovs% im_cl
    n13 =   sp_size * n_im_ch*n_im_ch_v + n12 ! -> tovs% im_ch
    n14 =   sp_size * nband             + n13 ! -> tovs% cldlev

    if (spot%p%n < n14 .and. tovs_io_ /= TTOVS_ALL) call finish('store_tovs',&
         'increasing size of t_tovs only possible with tovs_io == TTOVS_ALL.')

    call new_par (obs, n14, spot=spot)

    if (iand(tovs_io_,TTOVS_BASE)>0) then
      if (iand(tovs% init, TTOVS_BASE) == 0) call finish('store_tovs', &
           'refuse to store t_tovs with TTOVS_BASE not set in t_tovs%init.')
        obs% par (spot%p%i   +1 : spot%p%i+n1) = transfer(tovs                     ,par)
    end if
    if (associated(tovs% ci  )) then
      if (present(ci)) then
        obs% par (spot%p%i+n1+1 : spot%p%i+n2) = transfer(      ci  (1:nchan)      ,par)
      else if (iand(tovs_io_,TTOVS_CI)>0) then
        if (iand(tovs% init, TTOVS_CI) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%ci with TTOVS_CI not set in t_tovs%init.')
        obs% par (spot%p%i+n1+1 : spot%p%i+n2) = transfer(tovs% ci  (1:nchan)      ,par)
      end if
    end if
    if (associated(tovs% av  )) then
      if (present(av)) then
        obs% par (spot%p%i+n2+1 : spot%p%i+n3) = transfer(         av(1:nlev,1:nav),par)
      else if (iand(tovs_io_,TTOVS_AV)>0) then
        if (iand(tovs% init, TTOVS_AV) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%av with TTOVS_AV not set in t_tovs%init.')
        obs% par (spot%p%i+n2+1 : spot%p%i+n3) = transfer(tovs% av   (1:nlev,1:nav),par)
      end if
    end if
    if (associated(tovs% cemi)) then
      if (present(cemi)) then
        obs% par (spot%p%i+n3+1 : spot%p%i+n4) = transfer(      cemi(1:ncemi)      ,par)
      else if (iand(tovs_io_,TTOVS_CEMI)>0) then
        if (iand(tovs% init, TTOVS_CEMI) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%cemi with TTOVS_CEMI not set in t_tovs%init.')
        obs% par (spot%p%i+n3+1 : spot%p%i+n4) = transfer(tovs% cemi(1:ncemi)      ,par)
      end if
    end if
    if (associated(tovs% l2c )) then
      if (present(l2c)) then
        obs% par (spot%p%i+n4+1 : spot%p%i+n5) = transfer(      l2c (1:nl2c )      ,par)
      else if (iand(tovs_io_,TTOVS_L2C)>0) then
        if (iand(tovs% init, TTOVS_L2C) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%l2c with TTOVS_L2C not set in t_tovs%init.')
        obs% par (spot%p%i+n4+1 : spot%p%i+n5) = transfer(tovs% l2c (1:nl2c )      ,par)
      end if
    end if
    if (associated(tovs% emis)) then
      if (present(emis)) then
        obs% par (spot%p%i+n5+1 : spot%p%i+n6) = transfer(      emis(1:nchan)      ,par)
      else if (iand(tovs_io_,TTOVS_EMIS)>0) then
        if (iand(tovs% init, TTOVS_EMIS) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%emis with TTOVS_EMIS not set in t_tovs%init.')
        obs% par (spot%p%i+n5+1 : spot%p%i+n6) = transfer(tovs% emis(1:nchan)      ,par)
      end if
    end if
    if (associated(tovs% flag)) then
      if (present(flag)) then
        obs% par (spot%p%i+n6+1 : spot%p%i+n7) = transfer(      flag(1:nchan)      ,par)
      else if (iand(tovs_io_,TTOVS_FLAG)>0) then
        if (iand(tovs% init, TTOVS_FLAG) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%flag with TTOVS_FLAG not set in t_tovs%init.')
        obs% par (spot%p%i+n6+1 : spot%p%i+n7) = transfer(tovs% flag(1:nchan)      ,par)
      end if
    end if
    if (associated(tovs% sinfl)) then
      if (present(sinfl)) then
        obs% par (spot%p%i+n7+1 : spot%p%i+n8) = transfer(      sinfl(1:nchan)      ,par)
      else if (iand(tovs_io_,TTOVS_SINFL)>0) then
        if (iand(tovs% init, TTOVS_SINFL) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%sinfl with TTOVS_SINFL not set in t_tovs%init.')
        obs% par (spot%p%i+n7+1 : spot%p%i+n8) = transfer(tovs% sinfl(1:nchan)      ,par)
      end if
    end if
    if (associated(tovs% spec)) then
      if (present(spec)) then
        obs% par (spot%p%i+n8+1 : spot%p%i+n9) = transfer(      spec(1:nspec)      ,par)
      else if (iand(tovs_io_,TTOVS_SPEC)>0) then
        if (iand(tovs% init, TTOVS_SPEC) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%spec with TTOVS_SPEC not set in t_tovs%init.')
        obs% par (spot%p%i+n8+1 : spot%p%i+n9) = transfer(tovs% spec(1:nspec)      ,par)
      end if
    end if
    if (associated(tovs% bt_cs)) then
      if (present(bt_cs)) then
        obs% par (spot%p%i+n9+1 : spot%p%i+n10) = transfer(      bt_cs(1:nchan)      ,par)
      else if (iand(tovs_io_,TTOVS_BTCS)>0) then
        if (iand(tovs% init, TTOVS_BTCS) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%bt_cs with TTOVS_BTCS not set in t_tovs%init.')
        obs% par (spot%p%i+n9+1 : spot%p%i+n10) = transfer(tovs% bt_cs(1:nchan)      ,par)
      end if
    end if
    if (associated(tovs% tr  )) then
      if (present(tr)) then
        obs% par (spot%p%i+n10+1 : spot%p%i+n11) = transfer(      tr(1:nchan,1:ntr),par)
      else if (iand(tovs_io_,TTOVS_TR)>0) then
        if (iand(tovs% init, TTOVS_TR) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%tr with TTOVS_TR not set in t_tovs%init.')
        obs% par (spot%p%i+n10+1 : spot%p%i+n11) = transfer(tovs% tr(1:nchan,1:ntr),par)
      end if
    end if
    if (associated(tovs% im_cl)) then
      if (present(im_cl)) then
        obs% par (spot%p%i+n11+1 : spot%p%i+n12) = transfer(im_cl(1:n_im_cl_v,1:n_im_cl),par)
      else if (iand(tovs_io_,TTOVS_IMCL)>0) then
        if (iand(tovs% init, TTOVS_IMCL) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%im_cl with TTOVS_IMCL not set in t_tovs%init.')
        obs% par (spot%p%i+n11+1 : spot%p%i+n12) = transfer(tovs% im_cl(1:n_im_cl_v,1:n_im_cl),par)
      end if
    end if
    if (associated(tovs% im_ch)) then
      if (present(im_ch)) then
        obs% par (spot%p%i+n12+1 : spot%p%i+n13) = transfer(im_ch(1:n_im_ch_v,1:n_im_ch),par)
      else if (iand(tovs_io_,TTOVS_IMCH)>0) then
        if (iand(tovs% init, TTOVS_IMCH) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%im_ch with TTOVS_IMCH not set in t_tovs%init.')
        obs% par (spot%p%i+n12+1 : spot%p%i+n13) = transfer(tovs% im_ch(1:n_im_ch_v,1:n_im_ch),par)
      end if
    end if
    if (associated(tovs% cldlev)) then
      if (present(cldlev)) then
        obs% par (spot%p%i+n13+1 : spot%p%i+n14) = transfer(      cldlev(1:nband)      ,par)
      else if (iand(tovs_io_,TTOVS_CLDLEV)>0) then
        if (iand(tovs% init, TTOVS_CLDLEV) == 0) call finish('store_tovs', &
             'refuse to store t_tovs%cldlev with TTOVS_CLDLEV not set in t_tovs%init.')
        obs% par (spot%p%i+n13+1 : spot%p%i+n14) = transfer(tovs% cldlev(1:nband)      ,par)
      end if
    end if

  end subroutine store_tovs
!------------------------------------------------------------------------------
  subroutine destruct_tovs (tovs)
    type (t_tovs), intent(inout) :: tovs      ! TOVS specific data type
    if (iand(tovs%init,TTOVS_CI  )>0) then
      if (associated(tovs% ci  )) deallocate (tovs% ci  )
    end if
    if (iand(tovs%init,TTOVS_AV  )>0) then
      if (associated(tovs% av  )) deallocate (tovs% av  )
    end if
    if (iand(tovs%init,TTOVS_CEMI)>0) then
      if (associated(tovs% cemi)) deallocate (tovs% cemi)
    end if
    if (iand(tovs%init,TTOVS_L2C )>0) then
      if (associated(tovs% l2c )) deallocate (tovs% l2c )
    end if
    if (iand(tovs%init,TTOVS_EMIS)>0) then
      if (associated(tovs% emis)) deallocate (tovs% emis)
    end if
    if (iand(tovs%init,TTOVS_FLAG)>0) then
      if (associated(tovs% flag)) deallocate (tovs% flag)
    end if
    if (iand(tovs%init,TTOVS_SINFL)>0) then
      if (associated(tovs% sinfl)) deallocate (tovs% sinfl)
    end if
    if (iand(tovs%init,TTOVS_SPEC)>0) then
      if (associated(tovs% spec)) deallocate (tovs% spec)
    end if
    if (iand(tovs%init,TTOVS_BTCS)>0) then
      if (associated(tovs% bt_cs)) deallocate (tovs% bt_cs)
    end if
    if (iand(tovs%init,TTOVS_TR  )>0) then
      if (associated(tovs% tr  )) deallocate (tovs% tr  )
    end if
    if (iand(tovs%init,TTOVS_IMCL)>0) then
      if (associated(tovs% im_cl)) deallocate (tovs% im_cl)
    end if
    if (iand(tovs%init,TTOVS_IMCH)>0) then
      if (associated(tovs% im_ch)) deallocate (tovs% im_ch)
    end if
    if (iand(tovs%init,TTOVS_CLDLEV)>0) then
      if (associated(tovs% cldlev)) deallocate (tovs% cldlev)
    end if
    call construct(tovs)    ! to set tovs%init = 0
  end subroutine destruct_tovs
!------------------------------------------------------------------------------
  subroutine construct_tovs(tovs)
    type(t_tovs), intent(out) :: tovs
  end subroutine construct_tovs
!------------------------------------------------------------------------------
  subroutine construct_tovs_instr(ti)
    type(t_tovs_instr), intent(out) :: ti
  end subroutine construct_tovs_instr
!------------------------------------------------------------------------------
  subroutine get_tovs_rs(tovs, i_rs, rs, ti)
    type(t_tovs),             intent(in)            :: tovs
    integer,                  intent(out), optional :: i_rs
    type(t_rad_set), pointer, intent(out), optional :: rs
    type(t_tovs_instr),       intent(out), optional :: ti

    logical, parameter :: l_stop = .true.
    integer :: i, iset

    iset = -1
    if (tovs%id_rs > 0) then
      if (tovs%id_rs <= n_set) then
        if (tovs%id_rs == rad_set(tovs%id_rs)%id) then
          ! Usually this should be the case
          iset = tovs%id_rs
        end if
      end if
      do i = 1, n_set
        if (tovs%id_rs == rad_set(i)%id) then
          iset = i
          exit
        end if
      end do
    end if

    if (iset > 0) then
      if (present(i_rs)) i_rs   =  iset
      if (present(rs  )) rs     => rad_set(iset)
      if (present(ti  )) call get_tovs_instr(tovs, rad_set(iset), ti)
    else
      if (present(i_rs   )) i_rs    =  -1
      if (present(rs     )) rs      => null()
      if (present(ti     )) call construct(ti)
      if (l_stop) then
        write(0,*) 'get_tovs_rs failed:',tovs%id_rs,n_set,rad_set(1:n_set)%id
        call print_rad_set(rad_set(1),unit=0)
        call finish('get_tovs_rs','Did not find rad_set corresponding to t_tovs')
      end if
    end if

  end subroutine get_tovs_rs

  subroutine get_tovs_instr(tovs, rs, ti)
    type(t_tovs),       intent(in)  :: tovs
    type(t_rad_set),    intent(in)  :: rs
    type(t_tovs_instr), intent(out) :: ti

    integer :: ii, i
    logical :: l_new

    if (iand(tovs% init,TTOVS_CI+TTOVS_CI_EXT) == 0) &
         call finish('get_tovs_instr','t_tovs% ci not available.')

    ti%n_instr   = 0
    ti%o_ch_i(:) = 0
    ti%n_ch_i(:) = 0
    ti%ii(:)     = 0
    ii = 1
    do i = 1, tovs%nchan
      do while (tovs%ci(i) > rs%o_ch_i(ii) + rs%n_ch_i(ii))
        ii = ii + 1
      end do
      if (ti%n_instr == 0) then
        l_new = .true.
      else
        l_new = (ti%ii(ti%n_instr) /= ii)
      end if
      if (l_new) then
        ti%n_instr = ti%n_instr + 1
        if (ti%n_instr > m_instr) call finish('get_tovs_instr','too many instruments in t_tovs')
        ti%ii    (ti%n_instr) = ii
        ti%o_ch_i(ti%n_instr) = i-1
        ti%n_ch_i(ti%n_instr) = 1
      else
        ti%n_ch_i(ti%n_instr) = ti%n_ch_i(ti%n_instr) + 1
      end if
    end do

  end subroutine get_tovs_instr


  subroutine link_tovs_rs(sp, obs)
    type(t_spot), intent(inout) :: sp
    type(t_obs),  intent(inout) :: obs

    character(len=*),parameter :: proc = 'link_tovs_rs'
    integer,         parameter :: tovs_io = TTOVS_BASE + TTOVS_CI
    type(t_tovs)               :: tovs
    type(t_rad_set), pointer   :: rs => null()
    integer                    :: i1, in
    integer                    :: iset, ii, ic, i, j
    integer                    :: ci(m_chan)

    call load(obs, sp, tovs=tovs, ci=ci, tovs_io=tovs_io)
    if (tovs%id_rs > 0 .and. all(ci(1:tovs%nchan) > 0)) return

    iset = set_indx(rad_set, satid=int(sp%hd%satid), grid=int(sp%hd%grid_id))
    rs => rad_set(iset)
    tovs% id_rs = rs% id
    if (tovs%nchan /= sp%o%n) call finish(proc,'inconsistent channel number')
    ii = 1
    ic = 1
    do i = 1, tovs%nchan
      ci(i) = -1
      do while (obs% body(sp%o%i+i)% lev_sig /= rs%instr_wmo(ii) .and. ii <= rs%n_instr)
        ii = ii + 1
      end do
      if (ii > rs%n_instr) then
        call print_rad_set(rs, header='rs of crash in '//proc, unit=0)
        write(0,*) 'sp%hd%id=',sp%hd%id,' i=',i,' ii=',ii,' instr_wmo=', rs%instr_wmo(ii)
        write(0,*) 'lev_sig=',obs% body(sp%o%i+1:sp%o%i+sp%o%n)% lev_sig
        call finish(proc,'failed to find instrument')
      end if
      i1 = rs%o_ch_i(ii)+1
      in = rs%o_ch_i(ii)+rs%n_ch_i(ii)
      ic = max(ic, i1)
      ! Rely on correctly sorted channels
      do while (ic <= in)
        if (obs% olev(sp%o%i+i) == rs% chan(ic)) then
          ci(i) = ic
          exit
        else
          ic = ic + 1
        end if
      end do
      ! If that does not work, use slow search
      if (ci(i) <= 0) then
        do j = i1, in
          if (obs% olev(sp%o%i+i) == rs% chan(j)) then
            ci(i) = j
            exit
          end if
        end do
      end if
      if (ci(i) <= 0) then
        call print_rad_set(rs, header='rs of crash in '//proc, unit=0)
        write(0,*) 'sp%hd%id=',sp%hd%id,' i=',i,' ii=',ii,' ic=',ic,&
             ' lev_sig=',obs% body(sp%o%i+i)% lev_sig,' olev=',obs% olev(sp%o%i+i)
        call finish(proc,'channel/instrument not found')
      end if
    end do
    call store_tovs(obs, sp, tovs=tovs, ci=ci, tovs_io=tovs_io)

  end subroutine link_tovs_rs

  subroutine add_av_cont(tovs, id, i)
    type(t_tovs), intent(inout)          :: tovs
    integer(i8),  intent(in)             :: id  ! COL_* to be added
    integer,      intent(out),  optional :: i   ! index of id in t_tovs%av_cont
    integer :: j, ii

    ii = -1
    do j = 1, tovs%nav
      if (tovs%av_cont(j) == id) then
        ii = j
        exit
      end if
    end do

    if (ii < 0) then
      if (tovs%nav < mx_av) then
        tovs%nav = tovs%nav + 1
        ii = tovs%nav
        tovs%av_cont(ii) = id
      else
        do j = 1, tovs%nav
          write(0,'("tovs%av_cont(",I2,") = ",I20)') j, tovs%av_cont(j)
        end do
        write(0,'("Failed to add tovs%av_cont(",I2,") = ",I20)') tovs%nav+1, id
        call finish('add_av_cont@mo_t_tovs', 'av_cont too small. Consider to increase mx_av')
      end if
    end if
    if (present(i)) i = ii
  end subroutine add_av_cont

  subroutine get_im_ch_ind(tovs, imch, i)
    type(t_tovs), intent(inout) :: tovs
    integer,      intent(in)    :: imch
    integer,      intent(inout) :: i
    integer :: j
    i = -1
    do j = 1, tovs%n_im_ch_v
      if (tovs%im_ch_v(j) == imch) then
        i = j
        return
      elseif (tovs%im_ch_v(j) <= 0) then
        tovs%im_ch_v(j) = imch
        i = j
        return
      end if
    end do
    if (i < 0) then
      write(0,*) 'tovs%im_ch_v',tovs%im_ch_v(1:tovs%n_im_ch_v)
      write(0,*) 'imch',imch
      call finish('get_im_ch_ind', 'Failed to get entry in im_ch_v.')
    end if
  end subroutine get_im_ch_ind

end module mo_t_tovs
