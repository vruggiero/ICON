!
!+ Routines for variational quality control
!
MODULE mo_vqc
!
! Description:
!   This module holds routines for variational quality control.
!
!   For a given observation error covariance matrix and observation minus
!   analyses the cost function, its derivative and optionally the
!   inverse of its Hessian matrix (second derivative) are calculated by
!   the module procedures 'vqc_uncor', 'vqc_wind' and 'vqc_corr' for
!   uncorrelated, wind, and correlated data, respectively.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_2         2008/12/04 Andreas Rhodin
!  vqc_wind: remove unused variable e
! V1_7         2009/08/24 Andreas Rhodin
!  Print out complete namelist group /VARQC/
! V1_9         2010/04/20 Harald Anlauf
!  vqc_uncor: do not use w_min for iappr==1 for formulation 1
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_14        2011/11/08 Harald Anlauf
!  Optimization for sxf90 rev.430
! V1_15        2011/12/06 Andreas Rhodin
!  clean up handling of invalid observations in routines vqc_uncor, vqc_wind
! V1_20        2012-06-18 Andreas Rhodin
!  cleanup: remove unused variables
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2002-2008  original source
! Oliver Schmid   DWD  2005       alternative formulations
!==========================================================================
  !=============
  ! Modules used
  !=============
  !-------------------------------------------------
  ! Vector and matrix operations on 1D and 2D arrays
  !-------------------------------------------------
  use mo_matrix,    only: inverse,       &! matrix inversion (SVD)
                          inverse_rs,    &! matrix inversion (real symmetric)
                          check_rs,      &! check for eigenvalues
                          diag,          &! diagonal of matrix
                          operator (.x.),&! matrix multiplication
                          operator (.o.)  ! outer product
  !----------------
  ! kind parameters
  !----------------
  use mo_kind,      only: wp,            &! working precision
                          sp,            &! single  precision
                          i8              ! 8-byte integer kind parameter
  !---------------------------
  ! utilities to read namelist
  !---------------------------
  use mo_namelist,  only: position_nml,  &! routine to position nml group
                          nnml,          &! namelist fortran unit number
                          POSITIONED      ! position_nml: OK    return flag
  !----------------
  ! parallelization
  !----------------
  use mo_mpi_dace,  only: dace,          &! mpi communication info
                          p_bcast         ! broadcast routine
  !------------------
  ! write diagnostics
  !------------------
  use mo_p_output,  only: oline,         &! output line buffer
                          iol,           &! number of next line to write
                          nextline        ! routine to increment line number
  !-----------------
  ! error processing
  !-----------------
  use mo_exception, only: finish          ! abort routine
  implicit none
!------------------------------------------------------------------------------
  !================
  ! Public entities
  !================
  private
  public :: vqc_uncor    ! variational quality control rout., uncorrelated data
  public :: vqc_wind     ! variational quality control routine, wind data
  public :: vqc_corr     ! variational quality control routine, correlated data
  public :: t_vqc        ! Data type to pass parameters to 'vqc_corr'
  public :: read_vqc_nml ! read namelist /VARQC/
  public :: set_vqc_def  ! set default values for namelist /VARQC/
  public :: print_vqc_nml! print namelist /VARQC/
  public :: g_rej        ! function: a priory probability from rejection limit
  public :: s_rej        ! function: rejection limit from a priory probability
  public :: svqc         ! variational quality control threshold
  public :: mvqc, gvqc
  public :: w_min        ! minimum eigenvalue of the returned Hessian
  public :: nset, iset   ! max./current set of convergence parameters
  public :: iappr        ! approximation of the returned inverse Hessian
  public :: use_wind     ! flag to use vqc_wind routine
  public :: use_corr     ! flag to use vqc_corr routine
  public :: vqc_form     ! formulation of VQC: shape of J_o
!------------------------------------------------------------------------------
  !======================
  ! Data type definitions
  !======================
  !----------------------------------------------------------------------
  ! Data type 't_vqc' is passed to routine 'vcq_corr' in order to specify
  ! the limits of gross errors for a specific correlated observation.
  ! Currently only components 'm_rej', 'g_rej' are used.
  !----------------------------------------------------------------------
  type t_vqc
    integer  :: n_up    = 0
    integer  :: n_low   = 0
    integer  :: m_rej   = 0     ! Maximum number of rejected data
    real(wp) :: g_rej   = 0._wp ! A priory probability of rejected datum
    real(wp) :: g_up    = 0._wp
    real(wp) :: g_low   = 0._wp
    real(wp) :: g_diag  = 0._wp
  end type t_vqc

  !-----------------------------------------------------------------------
  ! The private data type 't_vqc_par' holds parameters to specify the
  ! strategy of routine 'vcq_corr' to compute an approximation of the cost
  ! function and its derivatives. An array 'vqc_param' of this data type
  ! defines a set of strategies used for different problem sizes.
  !
  ! Component 'dim_lim' specifies the limit on the number of correlated
  ! observations for a set of observations. The remaining parameters
  ! define the parameters used by routine 'vcq_corr' for observations
  ! belonging to this set.  These parameters are preliminary and may
  ! change in the future.
  !-----------------------------------------------------------------------
  type t_vqc_par
    integer          :: dim_lim = 0 ! limit on the number of dimensions (data)
    integer          :: n_full  = 0 ! exact up to this number of rejected data
    character(len=8) :: c_start = ''! start value for exploration
    character(len=8) :: c_final = ''! exploration strategy in final phase
    integer          :: n_lev_s = 0 ! number of levels to explore in startphase
    integer          :: n_lev_f = 0 ! number of levels to explore in finalphase
  end type t_vqc_par

  type (t_vqc_par) ,save :: empty ! set by default initialization
  !---------------------------------------------------------------------------
  ! Values for component 'c_start' are:
  !
  ! ''         ! start from top level (e.g. n_full)
  ! 'mostprob' ! start from most probable combination found
  ! 'mostpunc' ! start from most probable combination assuming no correlations
  !
  ! Values for component 'c_final' are:
  !
  ! ''         ! dont continue, all done
  ! 'descend'  ! continue until most probable state is found
  !---------------------------------------------------------------------------
!------------------------------------------------------------------------------
  !=================
  ! Namelist /VARQC/
  !=================
  !-------------------------------------------------------------------------
  ! Namelist /VARQC/ is read by subroutine 'read_vqc_nml'. (Default values
  ! are set by subroutine 'set_vqc_def'. The namelist provides parameters
  ! specifying the general behaviour of the variational quality control
  ! routine, mainly the approximation strategies followed by subroutine
  ! 'vqc_corr':
  !
  ! 'iprint'   determines the level of printout (0: none, 1:short, 2: long).
  ! 'w_min'    is the minimum eigenvalue of the returned Hessian.
  ! 'iappr'    determines the approximation of the returned inverse Hessian
  !            0: The inverse Hessian is not calculated
  !            1: an approximation is used (allways positive semidefinit)
  !            2: The exact Hessian (despite constraints on the eigenvalues
  !               is returned. Eigenvalues may be constrained as
  !               specified by parameter 'w_min'.
  !            3: As 2, but eigenvalues of R^(-1)+B^(-1) are constrained.
  ! 'vqc_par'  Approximation strategies for subroutine 'vqc_corr'
  ! 'vqc_form' Different formulations for the shape of J_o :
  !            1: J_o corresponds to a pdf which is a superposition of
  !               a Gaussian and a flat distribution. This formulation
  !               leads to a complex cost function with multiple
  !               minima and is not used any more.
  !            2: Huber Norm: J_o is quadratic for small |x| beyound
  !               the threshold and linear for larger values.
  !            3: as (2) but constant J_o for even larger |x|, not used.
  !            4: J_o has the shape of a hyperbola. This formulation
  !               is recommendet. It is very similar to the Huber Norm
  !               (2) but has a smooth transition from the quadratic
  !               to the constant part and leads to faster convergence
  !               of the variational scheme.
  !-------------------------------------------------------------------------
  integer            :: iprint     = 0        ! print flag
  real(wp)           :: w_min      = 1.e-2_wp ! minimum SV in Hessian

  integer ,parameter :: nset       = 2        ! 2 sets of different parameters
  integer            :: iset       = 1        ! set currently used
  integer            :: iappr(nset)= 2        ! approximation level for R

  logical            :: use_wind   = .true.   ! use vqc_wind
  logical            :: use_corr   = .false.  ! use vqc_corr
  integer            :: vqc_form   = 4        ! formulation of obs-costfunction
  !-------------------------------------
  ! parameters for outdated vqc_form = 1
  !-------------------------------------
  real(sp)           :: svqc       = 2._sp    ! variational quality control sigma
  integer            :: mvqc       = 3        ! max number of rejected observation
  real(wp)           :: gvqc       = 0._wp
  real(wp)           :: wrej       = 1._wp    ! weight for all rejected obsv.
  type (t_vqc_par) ,save :: vqc_par(10)       ! set of strategies

  namelist /varqc/ w_min, iappr, vqc_par, iprint, &
                   svqc, mvqc, use_wind, use_corr, wrej, vqc_form

!==============================================================================
contains
!==============================================================================
  subroutine read_vqc_nml
  !-------------------------------------
  ! subroutine to read namelist /VARQC/.
  !-------------------------------------
    integer          :: ierr  ! error flag
    integer          :: i
    !-------------
    ! set defaults
    !-------------
    call set_vqc_def

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('VARQC', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=VARQC, iostat=ierr)
        if (ierr/=0) call finish ('nml_run_flags','ERROR in namelist /VARQC/')
#else
        read (nnml ,nml=VARQC)
#endif
      end select
      !-------------------
      ! consistency checks
      !-------------------
      svqc = min (svqc, 30._sp)
      gvqc = 0._wp; if (svqc>0._sp) gvqc = 1._wp / exp (0.5_wp * svqc **2)
      do i=2,nset
        if(iappr(i)<0) iappr(i) = iappr(i-1)
      end do
!write(6,*)'read_vqc_nml: svqc,gvqc=',svqc,gvqc
    endif
    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    call p_bcast (w_min            ,dace% pio)
    call p_bcast (vqc_par% dim_lim ,dace% pio)
    call p_bcast (vqc_par% n_full  ,dace% pio)
    call p_bcast (vqc_par% c_start ,dace% pio)
    call p_bcast (vqc_par% c_final ,dace% pio)
    call p_bcast (vqc_par% n_lev_s ,dace% pio)
    call p_bcast (vqc_par% n_lev_f ,dace% pio)
    call p_bcast (svqc             ,dace% pio)
    call p_bcast (gvqc             ,dace% pio)
    call p_bcast (mvqc             ,dace% pio)
    call p_bcast (use_wind         ,dace% pio)
    call p_bcast (use_corr         ,dace% pio)
    call p_bcast (wrej             ,dace% pio)
    call p_bcast (vqc_form         ,dace% pio)
    call p_bcast (iappr            ,dace% pio)
    call p_bcast (iprint           ,dace% pio)
    !---------------
    ! print namelist
    !---------------
    call print_vqc_nml

  end subroutine read_vqc_nml
!------------------------------------------------------------------------------
  subroutine print_vqc_nml
  !-----------------------------
  ! printout of namelist /VARQC/
  !-----------------------------
    integer :: i
    if (dace% lpio) then
      write (6,'()')
      write (6,'(a)') repeat ('-',79)
      write (6,'()')
      write (6,'(a)') ' Variational Quality Control Parameters'
      write (6,'()')
      write (6,'(a)') ' Namelist /VARQC/:'
      write (6,'()')
      write (6,'(a,es8.2,a)')' w_min     = ',w_min, &
                             ' : minimum singular value in Hessian'
      write (6,'(a,2i4  ,a)')' iappr     = ',iappr, &
                             ' : approximation level for R'
      write (6,'(a,i8  ,a)') ' iprint    = ',iprint, &
                             ' : print flag'
      write (6,'(a,i8  ,a)') ' vqc_form  = ',vqc_form, &
                             ' : formulation of obs-costfunction'
      write (6,'(a,es8.2,a)')' svqc      = ',svqc, &
                             ' : variational quality control sigma'
      write (6,'(a,es8.2,a)')' gvqc      = ',gvqc
      write (6,'(a,i8  ,a)') ' mvqc      = ',mvqc, &
                             ' : max number of rejected observations'
      write (6,'(a,l8  ,a)') ' use_wind  = ',use_wind, &
                             ' : use vqc_wind'
      write (6,'(a,l8  ,a)') ' use_corr  = ',use_corr, &
                             ' : use vqc_corr'
      write (6,'(a,es8.2,a)')' wrej      = ',wrej, &
                             ' : weight for all rejected obsv.'
      write (6,'()')

      do i=1,size(vqc_par)
        if (vqc_par(i)% dim_lim == 0) exit
        write (6,'(a,i2,a)') ' vqc_par(',i,') :'
        write (6,'()')
        write (6,'(a,i8,a)') ' dim_lim = ',vqc_par(i)% dim_lim
        write (6,'(a,i8,a)') ' n_full  = ',vqc_par(i)% n_full
        write (6,'(a, a,a)') ' c_start = ',vqc_par(i)% c_start
        write (6,'(a, a,a)') ' c_final = ',vqc_par(i)% c_final
        write (6,'(a,i8,a)') ' n_lev_s = ',vqc_par(i)% n_lev_s
        write (6,'(a,i8,a)') ' n_lev_f = ',vqc_par(i)% n_lev_f
        write (6,'()')
      end do
    endif

  end subroutine print_vqc_nml
!------------------------------------------------------------------------------
  subroutine set_vqc_def
  !----------------------------------------
  ! set default values for namelist /VARQC/
  !----------------------------------------

    w_min               = 1.e-2_wp ! minimum singular value in Hessian
    iappr               = (/2,-1/) ! approximation level for inverse Hessian
    vqc_par             = empty    ! empty set of strategies
    iprint              = 0        ! print flag

    use_wind            = .true.   ! flag to use vqc_wind routine
    use_corr            = .false.  ! flag to use vqc_corr routine

    wrej                = 1._wp    ! weight for all rejected obsv.

    vqc_form = 4        ! formulation of obs-costfunction
    svqc     = 2._sp    ! variational quality control
    mvqc     = 3        ! maximum number of rejected observations

    vqc_par(1)% dim_lim = 10       ! full explorations up to 10 observations
    vqc_par(1)% n_full  = 10

    vqc_par(2)% dim_lim = 20       ! 10 rejections fully explored up to 20
    vqc_par(2)% n_full  =  5
    vqc_par(2)% c_start = 'mostprob'
    vqc_par(2)% n_lev_s = 2

    vqc_par(3)% dim_lim = 30       !  5 rejections fully explored up to 30
    vqc_par(3)% n_full  =  3
    vqc_par(3)% c_start = 'mostprob'
    vqc_par(3)% n_lev_s = 2

    vqc_par(4)% dim_lim = 50       !  3 rejections fully explored up to 50
    vqc_par(4)% n_full  =  2
    vqc_par(4)% c_start = 'mostprob'
    vqc_par(4)% n_lev_s = 1

    vqc_par(5)% dim_lim = 100      ! 2 rejections fully explored up to 100
    vqc_par(5)% n_full  =   1
    vqc_par(5)% c_start = 'mostprob'
    vqc_par(5)% n_lev_s = 1

    vqc_par(6)% dim_lim = 1000     ! 1 rejections fully explored up to 1000
    vqc_par(6)% n_full  =    1

  end subroutine set_vqc_def
!==============================================================================
  elemental function g_rej (sigma)
  !--------------------------------------------------------------------
  ! Subroutine to calculate a priory probability 'g_rej' as used in the
  ! variational quality control routines from a given rejection limit
  ! 'sigma' (standard deviation) for uncorrelated data.
  !--------------------------------------------------------------------
  real (wp) ,intent(in) :: sigma ! rejection limit (standard deviation)
  real (wp)             :: g_rej ! priory probability
    g_rej = 0._wp
    if (sigma > 0._wp) g_rej = exp (-0.5_wp * sigma**2)
  end function g_rej
!------------------------------------------------------------------------------
  elemental function s_rej (p_apr)
  !-----------------------------------------------------------------
  ! Subroutine to calculate the rejection limit (standard deviation)
  ! from the a priory probability 'g_rej'
  !-----------------------------------------------------------------
  real (wp) ,intent(in) :: p_apr ! rejection limit (standard deviation)
  real (wp)             :: s_rej ! priory probability
    s_rej = 0._wp
    if (p_apr > 0._wp) s_rej = sqrt(-2._wp * log(p_apr))
  end function s_rej
!==============================================================================
!
! For a given observation error covariance matrix and observation minus
! analyses the cost function, its derivative and optionally the
! inverse of its Hessian matrix (second derivative) are calculated by
! the module procedures 'vqc_uncor', 'vqc_wind' and 'vqc_corr' for
! uncorrelated, wind, and correlated data, respectively.
!
! The routines 'vqc_uncor', 'vqc_wind' and 'vqc_corr' calculate the cost
! function and its derivative for uncorrelated, wind, and correlated
! data, respectively. The Routines have similar arguments, in general they
! only differ in size and rank.
!
! 'x_o'   : Given analysis increment (analysis minus observation).
! 'R'     : Given observation error covariance matrix for correct data,
!           or its inverse, cf. parameter 'linv'.
! 'j_qc'  : The contribution to the cost function is added to this
!           intent(inout) variable.
! 'dj_qc' : Returned gradient of the cost function with respect to the analysis
! 'w_qc'  : Returned a posteriori probability that the datum is rejected.
! 'R_qc'  : Returned inverse Hessian. In dependence of the namelist variable
!           'iappr' this may be an approximation:
!           0: no inverse Hessian is computed.
!           1: an approximation to the inverse Hessian is computed.
!           2: the exact inverse Hessian is computed (despite lower bounds
!              on the eigenvalues of H)
! 'apri'  : Specification of the a priory probability of error. In case of
!           correlated observations this variable is a derived data type
!           and specifies both a priori probability and maximum number
!           of rejected data.
! 'linv'  : if true, not R but its inverse (the Hessian) is returned.
! 'flag'  : Modifies the behaviour of the variational quality control:
!            0.: No VQC. Cost function and gradients for correct data are
!                returned. The inverse Hessian is not calculated
!           >0.: The apriory probability of errors is at most  g_rej ( flag).
!           <0.: The apriory probability of errors is at least g_rej (-flag).
! 'lprint': Diagnostic printout is given for (lprint==.true.) only.
! 'HBH'   : Background error matrix (optional) , used to constrain eigenvalues
!           if iappr==3.
!==============================================================================
  elemental subroutine vqc_uncor (x_o, j_qc, dj_qc, w_qc, R, R_qc, apri, &
                                  linv, flag, state, dwda, HBH, form)
  !----------------------------------------------
  ! quality control routine for uncorrelated data
  !----------------------------------------------
  real(wp) ,intent(in)          :: x_o   ! analysis minus observation
  real(wp) ,intent(inout)       :: j_qc  ! cost function
  real(wp) ,intent(out)         :: dj_qc ! gradient of j_qc
  real(wp) ,intent(out)         :: w_qc  ! posteriori probab. for correct data
  real(wp) ,intent(in)          :: R     ! obs.err. cov.matrix for correct data
  real(wp) ,intent(out)         :: R_qc  ! inverse Hessian or Hessian
  real(wp) ,intent(in)          :: apri  ! a priori stdev. for error
  logical  ,intent(in)          :: linv  ! return inverse matrix instead of R
  real(wp) ,intent(in)          :: flag  ! 0: no vqc, /=0: modify criterium
  integer  ,intent(in)          :: state ! observation status flag
  real(wp) ,intent(out),optional:: dwda  ! d w / d apri
  real(wp) ,intent(in) ,optional:: HBH   ! background error
  integer  ,intent(in) ,optional:: form  ! formulation of obs-costfunction

    real(wp) :: p     ! probability
    real(wp) :: edR   ! 1 / R
    real(wp) :: g     ! a priori probability for error
    real(wp) :: e     ! exponent in cost function calculation
    real(wp) :: H     ! Hessian matrix
    real(wp) :: Bi    ! inverse background error, scaled with R
    real(wp) :: sigma ! stdev. for error, unscaled
    integer  :: f     ! formulation of obs-costfunction

    edR   = 1._wp / R                    ! recip. obs. error covariance mat.
    f     = vqc_form; if (present(form)) f = form
    sigma = apri * sqrt (R)
    if (f==2) sigma = 0.5_wp  * sigma
    if (f==3) sigma = 0.75_wp * sigma
    !------------------------------------------------
    ! no contribution to cost function for state < 0)
    !------------------------------------------------
    if (state < 0) then
      w_qc  = 0._wp                   ! posteriori probab.
      dj_qc = 0._wp                   ! gradient of cost function
!     R_qc  = edR * max(0._wp, w_min) ! Hessian
      R_qc  = edR                     ! Hessian
      if (present(dwda)) dwda = 0._wp
    !--------------------------------------------------------------
    ! no VarQC for flag==0 or state>0 (observation always accepted)
    !--------------------------------------------------------------
    else if (flag == 0 .or. state > 0 .or. f==0) then
      j_qc  = j_qc + 0.5_wp * x_o**2 * edR ! contribution to cost function
      w_qc  = 1._wp                        ! posteriori probab. (correct data)
      dj_qc = x_o * edR                    ! gradient of cost function
      R_qc  = edR                          ! Hessian
      if(present(dwda)) dwda = 0._wp
    !------------------
    ! VarQC for flag/=0
    !------------------
    else
      !----------------------------------
      ! formulation 1: pbd based approach
      !----------------------------------
      if (f == 1) then
        edR    = 1._wp / R                        ! reciproce covariance matrix
        Bi = 0._wp;
        if (present(HBH)) then
           if ( HBH /= 0._wp) Bi = 1._wp / HBH! recipr. backgrond err.
        endif
        if (flag > 0) then
          g    = min (g_rej (apri), g_rej (flag)) ! a priori probab. for error
        else
          g    = max (g_rej (apri), g_rej (-flag))
        endif
        e      = exp (min (0.5_wp * x_o**2 * edR,&
                      real(range(1._wp),wp)))
        p      = g + 1._wp / e                  ! probability
        j_qc   = j_qc - log (p / (g+1._wp))     ! contribution to cost function
        w_qc   = 1._wp - g / p                  ! probability for correct data
        w_qc   = max (w_qc, 0._wp)              ! +++ for IBM xlf -O3 option
        dj_qc  = w_qc * (x_o * edR)             ! gradient of cost function
        if (iappr(iset) == 1) then
          R_qc                = edR * w_qc      ! approximated inverse Hessian
          if (w_qc <= 0) R_qc = edR * w_min     !   saveguard for large x_o
        else                                    ! exact inverse Hessian
          H = (w_qc - g * x_o**2 * edR / (p**2 * e)) * edR
          if (iappr(iset) == 2) then
            R_qc = max(H, w_min * edR)
          else
            R_qc = max(H, w_min * edR - Bi)
          endif
        endif
        if(present(dwda)) then
          dwda = - (1._wp/e) / (g*(1._wp/e))**2 !flag not accnt.
        endif
      !----------------------------------
      ! formulation 2: quadratic + linear
      !----------------------------------
      else if (f == 2) then
        if (abs(x_o).le.sigma) then
          edR   = 1._wp / R                    ! recip. obs. error covariance mat.
          j_qc  = j_qc + 0.5_wp * x_o**2 * edR ! contribution to cost function
          w_qc  = 1._wp                        ! posteriori probab. (correct data)
          dj_qc = x_o * edR                    ! gradient of cost function
          R_qc  = edR                          ! Hessian
          if(present(dwda)) dwda = 0._wp
        else if (x_o.gt.sigma) then
          edR   = 1._wp / R                    ! recip. obs. error covariance mat.
          j_qc  = j_qc + 0.5_wp * sigma**2 * edR &
                + (x_o - sigma) * sigma * edR  ! contribution to cost function
          w_qc  = sigma / x_o                  ! posteriori probab. (correct data)
          dj_qc = sigma * edR                  ! gradient of cost function
          if (iappr(iset)==1) then
            R_qc  = edR * max(w_qc,  w_min)    ! Hessian
          else
            R_qc  = edR * max(0._wp, w_min)    ! Hessian
          endif
          if(present(dwda)) dwda = 1._wp/x_o
        else
          edR   = 1._wp / R                    ! recip. obs. error covariance mat.
          j_qc  = j_qc + 0.5_wp * sigma**2 * edR &
                - (x_o + sigma) * sigma * edR  ! contribution to cost function
          w_qc  = - sigma / x_o                ! posteriori probab. (correct data)
          dj_qc = - sigma * edR                ! gradient of cost function
          if (iappr(iset)==1) then
            R_qc  = edR * max(w_qc,  w_min)    ! Hessian
          else
            R_qc  = edR * max(0._wp, w_min)    ! Hessian
          endif
          if(present(dwda)) dwda = - 1._wp/x_o
        endif
      !--------------
      ! formulation 3
      !--------------
      else if (f == 3) then
        if (abs(x_o).le.sigma) then
          edR   = 1._wp / R                    ! recip. obs. error covariance mat.
          j_qc  = j_qc + 0.5_wp * x_o**2 * edR ! contribution to cost function
          w_qc  = 1._wp                        ! posteriori probab. (correct data)
          dj_qc = x_o * edR                    ! gradient of cost function
          R_qc  = edR                          ! Hessian
          if(present(dwda)) dwda = 0._wp
        else if (x_o.gt.sigma) then
          if (x_o.gt.2._wp*sigma) then
            edR   = 1._wp / R                  ! recip. obs. error covariance mat.
            j_qc  = j_qc + sigma**2 * edR      ! contribution to cost function
            w_qc  = 0._wp                      ! posteriori probab. (correct data)
            dj_qc = 0._wp                      ! gradient of cost function
            R_qc  = edR * max(w_qc, w_min)     ! Hessian
            if(present(dwda)) dwda = 0._wp
          else
            edR   = 1._wp / R                  ! recip. obs. error covariance mat.
            j_qc  = j_qc + sigma**2 * edR - 0.5_wp * &
                  (x_o - 2._wp *sigma)**2 * edR ! contribution to cost function
            w_qc  = - (x_o - 2._wp *sigma) / x_o
                                               ! posteriori probab. (correct data)
            dj_qc = - (x_o - 2._wp *sigma) * edR
                                               ! gradient of cost function
            if (iappr(iset)==1) then
              R_qc  = edR * max(w_qc,   w_min) ! Hessian
            else
              R_qc  = edR * max(-1._wp, w_min) ! Hessian
            endif
            if(present(dwda)) dwda = 2._wp / x_o
         endif
        else
          if (x_o.lt.-2._wp*sigma) then
            edR   = 1._wp / R                  ! recip. obs. error covariance mat.
            j_qc  = j_qc + sigma**2 * edR      ! contribution to cost function
            w_qc  = 0._wp                      ! posteriori probab. (correct data)
            dj_qc = 0._wp                      ! gradient of cost function
            R_qc  = edR * max(w_qc, w_min)     ! Hessian
            if(present(dwda)) dwda = 0._wp
          else
            edR   = 1._wp / R                  ! recip. obs. error covariance mat.
            j_qc  = j_qc + sigma**2 * edR - 0.5_wp * &
                  (x_o + 2._wp *sigma)**2 * edR ! contribution to cost function
            w_qc  = - (x_o + 2._wp *sigma) / x_o
                                               ! posteriori probab. (correct data)
            dj_qc = - (x_o + 2._wp *sigma) * edR
                                               ! gradient of cost function
            if (iappr(iset)==1) then
              R_qc  = edR * max(w_qc,   w_min) ! Hessian
            else
              R_qc  = edR * max(-1._wp, w_min) ! Hessian
            endif
            if(present(dwda)) dwda = - 2._wp / x_o
          endif
        endif
      !----------------------------------
      ! formulation 4: hyperbel
      !----------------------------------
      else if (f == 4) then
        edR   = 1._wp / R                    ! recip. obs. error covariance mat.
        j_qc  = j_qc + sqrt (0.7_wp * apri * edR) * sqrt ( x_o**2 + &
                0.7_wp * apri / edR) - 0.7_wp * apri ! contribution to cost function
        w_qc  = sqrt (0.7_wp * apri * edR) / sqrt ( x_o**2 + &
                0.7_wp * apri / edR) / edR ! posteriori probab. (correct data)
        dj_qc = sqrt (0.7_wp * apri * edR) * x_o / sqrt ( x_o**2 + &
                0.7_wp * apri / edR) ! gradient of cost function
        if (iappr(iset)==1) then           ! Hessian
          R_qc  = edR * w_qc
        else
          R_qc  = sqrt (0.7_wp * apri * edR)                           &
                * ( 1._wp - x_o**2 / ( x_o**2 + 0.7_wp * apri / edR )) &
                / sqrt ( x_o**2 + 0.7_wp * apri / edR)
        endif
        if(present(dwda)) dwda = 0._wp
      endif
    endif
    !---------------------------
    ! process optional arguments
    !---------------------------
    if(.not.linv) then
      if (  R_qc  /=0._wp) R_qc  = 1._wp / R_qc
    endif
  end subroutine vqc_uncor
!==============================================================================
  subroutine vqc_wind (x_o, j_qc, dj_qc, w_qc, R, R_qc, apri, linv, flag, &
                       state, HBH, form)
  !----------------------------------------------------------------------
  ! quality control routine for wind data.
  !
  ! The observational errors of the two wind components are assumed to be
  ! uncorrelated. Only both components may be rejected simultaniously.
  !----------------------------------------------------------------------
  real(wp)     ,intent(in)   :: x_o  (2) ! analysis minus observation
  real(wp)     ,intent(inout):: j_qc     ! cost function
  real(wp)     ,intent(out)  :: dj_qc(2) ! gradient of j_qc
  real(wp)     ,intent(out)  :: w_qc (2) ! posteriori probability for corr.data
  real(wp)     ,intent(in)   :: R        ! obs. error cov.matrix for corr. data
  real(wp)     ,intent(inout):: R_qc(2,2)! inverse Hessian
  real(wp)     ,intent(in)   :: apri     ! a priori probability for error
  logical      ,intent(in)   :: linv     ! return inverse matrix instead of R
  real(wp)     ,intent(in)   :: flag     ! 0: no vqc, /=0: modify criterium
  integer      ,intent(in)   :: state(2) ! observation status flags
  real(wp) ,intent(in) ,optional:: HBH  (2,2)! background error
  integer  ,intent(in) ,optional:: form  ! formulation of obs-costfunction

    real(wp) :: edR    ! 1 / R
    real(wp) :: g      !
    real(wp) :: edRl, edRt
    real(wp) :: Rl, Rt
    real(wp) :: x2,c,s,cc,cs,ss
    real(wp) :: Bi(2,2), tmp(2,2)
    integer  :: stat
    real(wp) :: dj_qcu ! dj_qc returned from vqc_ucor
    integer  :: f

    stat = minval (state)
    f    = vqc_form; if (present(form)) f = form
    edR  = 1._wp / R                     ! recip. obs. error covariance mat.
    !------------------------------------------------
    ! no contribution to cost function for state < 0)
    !------------------------------------------------
    if (stat < 0) then
      w_qc  = 0._wp                   ! posteriori probab.
      dj_qc = 0._wp                   ! gradient of cost function
!     R_qc  = edR * max(0._wp, w_min) ! Hessian
      R_qc (1,2) = 0._wp
      R_qc (2,1) = 0._wp
      R_qc (1,1) = edR
      R_qc (2,2) = edR
    !-------------------------------------------------------------
    ! no VarQC for flag==0 or stat>0 (observation always accepted)
    !-------------------------------------------------------------
    else if (flag == 0 .or. stat > 0 .or. f==0) then
      j_qc  = j_qc + 0.5_wp * sum(x_o**2) * edR ! contribution to cost function
      w_qc  = 1._wp                         ! posteriori probab. (correct data)
      dj_qc = x_o * edR                     ! gradient of cost function
      R_qc (1,2) = 0._wp
      R_qc (2,1) = 0._wp
      R_qc (1,1) = edR
      R_qc (2,2) = edR
    !------------------
    ! VarQC for flag/=0
    !------------------
    else
      x2    = sum(x_o**2)
      edR   = 1._wp / R                         ! reciprocal covariance matrix
      if (flag > 0) then
        g   = min (g_rej (apri), g_rej (flag))  ! a priori probab. for error
      else
        g   = max (g_rej (apri), g_rej (-flag))
      endif
      call vqc_uncor (sqrt(x2), j_qc, dj_qcu, w_qc(1), R, edRl, apri, &
                      .true., flag, stat, form=form)
      w_qc(2) = w_qc(1)
      dj_qc = w_qc * (x_o * edR)              ! gradient of cost function
      edRt  = w_qc(1) * edR                   ! Hessian,transversal direction
      Rt = max(edRt, w_min * edR)
      Rl =     edRl
      if (iappr(iset) == 1) then                ! approximated inverse Hessian
        R_qc(1,1) = Rt
        R_qc(2,2) = Rt
        R_qc(1,2) = 0._wp
        R_qc(2,1) = 0._wp
      else                                      ! exact inverse Hessian
        if (x2 == 0._wp) then                   ! zero wind case
          R_qc(1,1) = Rt
          R_qc(2,2) = Rt
          R_qc(1,2) = 0._wp
          R_qc(2,1) = 0._wp
        else                                    ! nonzero wind case:
          c  = x_o (1)                          ! Coordinate transformation
          s  = x_o (2)
          x2 = 1._wp / x2
          cc = c * c * x2
          ss = s * s * x2
          cs = c * s * x2
          R_qc(1,1) = cc * Rl + ss * Rt
          R_qc(2,2) = cc * Rt + ss * Rl
          R_qc(1,2) = cs * (Rl - Rt)
          R_qc(2,1) = R_qc (1,2)
        endif
      endif
    endif

    !----------------------------
    ! constrain vs. R^(-1)+B^(-1)
    !----------------------------
    if (iappr(iset)==3) then
      Bi   = 0._wp
      if(present(HBH)) Bi = inverse(HBH)
      tmp  = R_qc
      call check_rs ((R_qc+Bi)*R, min_ev=w_min, y=R_qc)
      R_qc = (R_qc - Bi) * edR
    endif

    !-------
    ! invert
    !-------
    if (.not.linv) then
     R_qc                    = inverse_rs (R_qc)
    endif

  end subroutine vqc_wind
!==============================================================================
  subroutine vqc_corr (x_o, j_qc, dj_qc, w_qc, R, R_qc, apri,         &
                       sgm, linv, flag, state, lprint, dwda, HBH, form)
  !--------------------------------------------
  ! quality control routine for correlated data
  !--------------------------------------------
  real(wp)   ,intent(in)          :: x_o    (:) ! observation increment
  real(wp)   ,intent(out)         :: j_qc       ! cost function
  real(wp)   ,intent(out)         :: dj_qc  (:) ! gradient of j_qc
  real(wp)   ,intent(out)         :: w_qc   (:) ! probability for correct data
  real(wp)   ,intent(in)          :: R    (:,:) ! obs. error for correct data
  real(wp)   ,intent(inout)       :: R_qc (:,:) ! obs.error used in Newton step
  type(t_vqc),intent(in)          :: apri       ! a priori probability flags
  real(wp)   ,intent(in)          :: sgm        ! threshold for tails
  logical    ,intent(in)          :: linv       ! return inverse R, not R
  real(wp)   ,intent(in)          :: flag       ! R modification flag
  integer    ,intent(in)          :: state  (:) ! observation status flags
  logical    ,intent(in)          :: lprint     ! enable printout
  real(wp)   ,intent(out),optional:: dwda       ! d w / d apri
  real(wp)   ,intent(in) ,optional:: HBH  (:,:) ! background error
  integer    ,intent(in) ,optional:: form       ! formulation of J_o
    !================
    ! local variables
    !================
    !-------------------
    ! general parameters
    !-------------------
    integer  :: n                          ! number of observations
    integer  :: minr                       ! min no of rejected obsv.
    integer  :: maxr                       ! max-min no of rejected obsv.
    integer  :: iloc(1)                    ! keeps result of maxloc function
    integer  :: ipr, ipr0                  ! local coppy of print flag
    integer  :: wall0, wall1
    real     :: cpu, cpu0, cpu1, wall
    integer  :: rate
    integer  :: f
    !----------------------------------------
    ! temporaries for univariate calculations
    !----------------------------------------
    real(wp) :: x_o_1
    real(wp) :: R_1
    real(wp) :: j_qc_1
    real(wp) :: dj_qc_1
    real(wp) :: R_qc_1
    real(wp) :: w_qc_1
    real(wp) :: fn      ! normalisation (sqrt(n))
    !--------------------------------------------------
    ! normalized covariances and observation increments
    !--------------------------------------------------
    real(wp) :: derr (size(x_o))           ! 1 / sqrt(diag(R)) (normalization)
    real(wp) :: S    (size(x_o),size(x_o)) ! normalization
    real(wp) :: C    (size(x_o),size(x_o)) ! correlation matrix
    real(wp) :: o    (size(x_o))           ! normalized analysis-observation
    real(wp) :: Bi   (size(x_o),size(x_o)) ! scaled inverse B-matrix
    !--------------------------------------------------------------------
    ! quantities valid for a specific combination (prozedure 'calc_prob')
    !--------------------------------------------------------------------
    real(wp) :: pi                         ! a posteriori probability
    real(wp) :: edRi (size(x_o),size(x_o)) ! 1 / R for combinations
    !-----------------------
    ! accumulated quantities
    !-----------------------
    real(wp) :: p                          ! accumulated probability
    real(wp) :: edR  (size(x_o),size(x_o)) ! 1 / R
    real(wp) :: djqc (size(x_o))
    real(wp) :: gsum                       ! total a priori probability
    real(wp) :: p0                         ! probability normalisation factor
    real(wp) :: wnrej (0:size(x_o))        ! probability of i rejections
    integer  :: nscan (0:size(x_o))        ! # compinations tested for i rej.
    !------------------------------------
    ! probabilities for uncorrelated data
    !------------------------------------
    real(wp) :: pu   (size(x_o))           ! probability for rejection
    logical  :: pru  (size(x_o))           ! most probable combination
    !----------------------------------------------
    ! description of combination actually evaluated
    !----------------------------------------------
    real(wp) :: proj (size(x_o))           ! projection (1.:acc, 0.:rej)
    logical  :: pri  (size(x_o))           ! initial combination to test
    logical  :: fix  (size(x_o))           ! positions kept constant
    real(wp) :: prj1 (size(x_o))           ! values of positions kept constant
    real(wp) :: g                          ! a priori probability flag
    logical  :: rest                       ! true for remaining rejections
    !------------------------------------------------
    ! index variables used for loop over combinations
    !------------------------------------------------
    integer          :: i, iu, il, i0, i1  ! index variables
    integer          :: n1, n2, n3         ! index variables
    integer          :: ii   (size(x_o))   ! index variables
    real(wp)         :: wi
    character(len=8) :: task
    integer          :: n_lev
    integer          :: m_rej              ! max number of positions to reject
    type(t_vqc_par)  :: par
    !--------------------------------------
    ! list of combinations evaluated so far
    !--------------------------------------
    integer, parameter :: m = 100000              ! maximum size of table
    integer            :: nm                      ! actual index in table
    integer            :: mm                      ! actual size of table
    integer(i8)        :: sigs (m,size(x_o)/64+1) ! bit mask for acc/rej data
    real(wp)           :: post (m)                ! posteriory probability
    integer            :: nval (m)                ! number of accepted data
    integer(i8)        :: sig    (size(x_o)/64+1) ! actual combination
    integer(i8)        :: sigalt (size(x_o)/64+1) ! last combination

    !============================================
    ! Executable statements: no VarQC for flag==0
    !============================================
    f = vqc_form; if (present(form)) f = form
                  if (flag == 0)     f = 0

    if (f == 1) then
      !==========================================
      ! VarQC for vqc_form == 1 (old formulation)
      !==========================================
      !-------
      ! timing
      !-------
      ipr  = 0;   if (lprint)          ipr  = iprint ! local copy of print flag
      ipr0 = ipr; if (.not.dace% lpio) ipr0 = 0      ! ipr0: print on IO-PE only
      if (ipr /= 0) then
        call cpu_time     (cpu0)
        call system_clock (count=wall0)
      endif
      select case (f)
      case (1)
      !-------------------------------------------------------------------
      ! J_o corresponds to a pdf which is a superposition of
      !               a Gaussian and a flat distribution. This formulation
      !               leads to a complex cost function with multiple
      !               minima and is not used any more.
      !-------------------------------------------------------------------
        !-------------------------------------
        ! set initial value of some quantities
        !-------------------------------------
        n      = size (x_o)  ! number of observations
        p      = 0._wp       ! a priori prob. for error in all data
        w_qc   = 0._wp       ! probability for error in specific datum
        dj_qc  = 0._wp       ! gradient of cost function
        p0     = 0._wp       ! probability normalisation factor
        edR    = 0._wp       ! Hessian
        wnrej  = 0._wp       ! probability for number of rejected data
        nscan  = 0           ! number of combinations tested for of rej. data
        gsum   = (1._wp + apri% g_rej) ** n ! priory prob. comb.not tested so far
        i0     = 0           ! loop indices
        i1     = 0           ! loop indices
        rest   = .false.     ! flag for final comp. of all data rejected
        mm     = 0           ! clear table of tested combinations
        sigalt = 0           ! bit mask of last best combination

        !----------------------------------------------
        ! scale observation increments and correlations
        !----------------------------------------------
        derr = sqrt (1._wp / diag (R))  ! scaling factor
        S    = (derr .o. derr)          ! scaling matrix
        C    = R   * S                  ! correlation matrix
        o    = x_o * derr               ! scaled observation increments
        Bi   = 0._wp
!       if(iappr(iset)>2 .and. present(HBH)) Bi = inverse (HBH * S) !+++IBM bug
        if(iappr(iset)>2 .and. present(HBH)) then
          Bi = inverse (HBH * S)
        endif
        !---------------------------------------------
        ! determine parameters for actual problem size
        !---------------------------------------------
        par = empty
        do i = 1, size (vqc_par)
          if (n <= vqc_par(i)% dim_lim) then
            par = vqc_par(i)
            exit
          else if (vqc_par(i)% dim_lim == 0) then
            par% dim_lim = 1000000
            par% n_full  = 1
            exit
          endif
        end do
        if (par% dim_lim == 0) &
          call finish('vqc_corr','no appropriate parameter set')
        !-----------------------------------
        ! masks for always active or passive
        !-----------------------------------
        fix                   = .false.
        where (state/=0) fix  = .true.
        prj1                  = 1._wp
        where (state<0)  prj1 = 0._wp
        m_rej = min (apri% m_rej, count(.not.fix)-1)
        !------------------------------
        ! printout of actual parameters
        !------------------------------
        if (ipr0 /= 0) then
          write (6,'()')
          write (6,'(a)') repeat('-',79)
          write (6,'()')
          write (6,'(a)') '         variational quality control'
          write (6,'()')
          write (6,'(a,i8)') ' dim     = ', n
          write (6,'(a,i8)') ' a%m_rej = ', apri% m_rej
          write (6,'(a,i8)') ' m_rej   = ', m_rej
          write (6,'(a,e14.6)') ' g_rej   = ', apri% g_rej
          write (6,'()')
          write (6,'(a,i8)') ' dim_lim = ', par% dim_lim
          write (6,'(a,i8)') ' n_full  = ', par% n_full
          write (6,'(a,a8)') ' c_start = ', par% c_start
          write (6,'(a,a8)') ' c_final = ', par% c_final
          write (6,'(a,i8)') ' n_lev_s = ', par% n_lev_s
          write (6,'(a,i8)') ' n_lev_f = ', par% n_lev_f
        endif
        !---------------------------------------------
        ! calculate probabilities of uncorrelated data
        !---------------------------------------------
        pu = apri% g_rej / (exp (-0.5_wp*o**2) + apri% g_rej)

        if (ipr0 /= 0) then
          write (6,'()')
          write (6,'(a)') '         probability for uncorrelated data'
          write (6,'()')
          write(6,'(a,2f10.6,a)')' min,max o-x:', minval(abs(o)),maxval(abs(o))
          write(6,'(a,2f10.6,a)')' min,max pi :', minval(exp (-0.5_wp*o**2)),&
                                                  maxval(exp (-0.5_wp*o**2))
          write(6,'(a,2f10.6,a)')' min,max pu :', minval(pu),maxval(pu)
          write(6,'(a,2f10.6,a)')' min,max pp :', apri% g_rej / (1+apri% g_rej)
          write(6,'()')
          write(6,'(a)')'    i iloc       o-x        pi        pu state'
        endif
        pru = .true.
        do i=1,n
          iloc = maxloc(pu)
            if (ipr0 /= 0) write(6,'(2i5,3f10.6,i2)') &
              i,iloc(1),o(iloc(1)),exp (-0.5_wp*o(iloc(1))**2),pu(iloc(1)), &
              state(iloc(1))
            if (task=='mostpunc' .and. pu(iloc(1))>0.5_wp .and. i<apri% m_rej)&
              pru(iloc(1)) = .false.
            if(fix(iloc(1))) pru(iloc(1)) = state(iloc(1)) == 1
          pu (iloc(1)) = - pu(iloc(1))
        end do
        pu = abs (pu)

        !----------------------------------
        ! set switches for full exploration
        !----------------------------------
        call set_task ('full', min (par% n_full, m_rej))
        pri = prj1 == 1._wp

        if (ipr0 /= 0) print *,'fix :',fix
        if (ipr0 /= 0) print *,'prj1:',nint(prj1)
        !===================================================
        ! loop over allowed combinations of error / no error
        !===================================================
        !----------
        ! all false
        !----------
        proj = 0._wp
        where (fix (:)) proj (:) = prj1 (:)
        maxr = count (proj==0._wp)
        !---------
        ! all true
        !---------
        proj = prj1
        minr = count (proj==0._wp)
        maxr = maxr - minr
        call calc_prob (lref=.true.)
        !---------------------------------
        ! loop over lower and upper bounds
        !---------------------------------
        do   il = 1, 1 + apri% n_low
          do iu = n, n - apri% n_up, -1
            if (il > iu) cycle
            i0 = 1
            i1 = n
            if (apri% n_low>0) i0 = 1 + apri% n_low
            if (apri% n_up >0) i1 = n - apri% n_up
!           if (apri% n_low>0) i0 = il + 1
!           if (apri% n_up >0) i1 = iu - 1

            proj = prj1
            !-----------------------
            ! loop over combinations
            !-----------------------
            call set_indices
            do
              do n3 = 1, n2
                ii(n3) = ii(n3) + 1
                do
                  if (ii(n3) > n1+1-n3)      exit
                  if (     1 > i0-1+ii(n3) ) exit
                  if (.not.fix(i0-1+ii(n3))) exit
                  ii(n3) = ii(n3) + 1
                end do
                if (ii(n3) <= n1+1-n3) exit
                if (n3==n2)            exit
                ii(n3) = 0
              end do
              do n3 = n2-1,1,-1
                if (ii(n3+1)>0) then
                  ii(n3) = max (ii(n3), ii(n3+1) + 1)
                endif
              end do
              !---------------
              ! set projection
              !---------------
              proj (    :il-1) = 0._wp
              proj (il  :iu  ) = 1._wp
              proj (iu+1:    ) = 0._wp
              do n3 = 1, n2
                if (ii(n3) > 0 .and. ii(n3) < n1+1) then
                  proj (i0-1+ii(n3)) = 0._wp
                endif
              end do
              where (.not.pri (il:iu)) proj (il:iu) = 1._wp - proj (il:iu)
              where (fix(il:iu)) proj (il:iu) = prj1(il:iu)
              i = count (proj (il:iu) == 0._wp)
              !-----------------------------------------
              ! action if all combinations are processed
              !-----------------------------------------
              if (any(ii(1:n2) > n1)) then
                if (ipr0 /= 0) print *,'best:',pri,count(pri)
                !------------------------------------
                ! determine best combination explored
                !------------------------------------
                select case (task)
                case ('mostprob')
                !-------------------------------------
                ! searching for most probable combination
                !-------------------------------------
                  if (n - count (pri) < m_rej) then
                    !-------------------------
                    ! continue with next level
                    !-------------------------
                    call best (n - count(pri)+1)
                  else
                    !--------------------------------
                    ! switch to next task if all done
                    !--------------------------------
                    call best
                    call set_task (par% c_final, par% n_lev_f)
                  endif
                case ('full')
                !-----------------
                ! full exploration
                !-----------------
                  if (n_lev >= n-1) then
                    !-----------------
                    ! exit if all done
                    !-----------------
                    exit
                  else
                    !--------------------
                    ! switch to next task
                    !--------------------
                    call best
                    call set_task (par% c_start, par% n_lev_s)
                  endif
                case ('descend')
                  call best
                end select
                !--------------------------
                ! set new start combination
                !--------------------------
                if (task == 'mostpunc') then
                  pri = pru
                  call set_task (par% c_final, par% n_lev_f)
                else
                  call sig2pri
                  if (ipr0 /= 0) &
                    print *,'best:',pri,count(pri),post(iloc(1))
                  !--------------------
                  ! exit if no progress
                  !--------------------
                  if (all(sigalt == sig)) exit
                endif
                sigalt = sig
                call set_indices
                cycle
              endif
              !------------------------------------
              ! skip combinations already processed
              !------------------------------------
              i = count(proj==0._wp) - minr
              if (i == 0   ) cycle ! all  accepted
              if (i == maxr) cycle ! all  rejected
              if (i > m_rej) cycle ! more rejected than allowed
              !--------------------------------------------------------
              ! calculate a posteriory probability for this combination
              !--------------------------------------------------------
              call calc_prob (lref=.false.)
            end do
          end do
        end do
!       !----------
!       ! all false
!       !----------
!       proj = 0._wp
!       call calc_prob
        !---------------------------
        ! all remaining combinations
        !---------------------------
        proj = 0._wp
        where (fix (:)) proj (:) = prj1 (:)
        rest = .true.
        call calc_prob (lref=.false.)
        !---------------------------
        ! calculate final quantities
        !---------------------------
        J_qc  = - log (p / p0)                 ! cost function
        dj_qc = dj_qc / p                      ! gradient of the cost function
        w_qc  = w_qc  / p                      ! probability for correct datum
        wnrej = wnrej / p
        edR   = edR / p &                      ! normalise Hessian
              + (dj_qc .o. dj_qc)
        if (iappr(iset)==3) then
          call check_rs (edR+Bi, min_ev=w_min, y=R_qc)
          R_qc  = R_qc  - Bi
        else
          call check_rs (edR, min_ev=w_min, y=R_qc)
        endif
        !--------------------
        ! denormalize, invert
        !--------------------
        dj_qc = dj_qc * derr
        R_qc  = R_qc  * S
        if (.not. linv) then
          R_qc                     = inverse_rs (R_qc)
        endif
        !---------
        ! printout
        !---------
        if (ipr0 /= 0) then
          write(6,'()')
          write(6,'(a,f12.6)') 'g_rej    :',apri% g_rej
          write(6,'(a,i5   )') 'm_rej    :',apri% m_rej
          write(6,'(a,2f10.6)') 'min,max w:',minval(w_qc),maxval(w_qc)
          write(6,'(a)') 'probability of rejected data:'
          do i = 0, n
            if (nscan(i) > 0) write(6,'(i4,f6.3,i6,i15,f9.6)') &
              i, wnrej (i), nscan(i), komb (n,i), real(nscan(i),wp)/komb (n,i)
          end do
          write(6,'()')
        endif
      case default
        write (0,*)  'vqc_corr:  VQC formulation',f,'not implemented'
        call finish ('vqc_corr','VQC formulation not implemented')
      end select
      !-------
      ! timing
      !-------
      if (ipr/=0) then
        call cpu_time     (cpu1)
        call system_clock (count=wall1)
        call system_clock (count_rate=rate)
        cpu  =  cpu1 - cpu0
        wall = (wall1-wall0)/rate
        wall = max (wall, cpu,1.e-6)
        call nextline
        write(oline(iol),'(a,i5,3f10.3)') &
          'n,cpu,wall,%:',n, cpu, wall, cpu/wall*100
      endif
    else
    !===========================================
    ! VarQC for vqc_form /= 1 (new formulations)
    !===========================================
      !--------------------------------------
      ! first calculate J_o for f==0 (no VQC)
      !--------------------------------------
      n     = size (x_o)                      ! number of observations
      R_qc  = R                               ! modified cov.matrix
      edR   = inverse (R)                     ! recip. covariance matrix
      J_qc  = 0.5_wp * sum ((x_o.x.edR) * x_o)! contribution to cost function
      w_qc  = 1._wp                           ! probability for correct data
      dj_qc = edR .x. x_o                     ! gradient of cost function
      if (present(dwda)) dwda = 0._wp         ! d w / d sigma
      if (f == 0) then
        if (linv) then
          R_qc = edR
        else
          R_qc  = R
        endif
        return
      else
        !--------------------------------------------------
        ! simple VQC for correlated data :
        ! 1) solve univariate problem for corresponding x_o
        ! 2) rescale the whole pdf
        !--------------------------------------------------
        fn     = sqrt (real(n,wp))
        x_o_1  = sqrt (2._wp * J_qc / fn)  ! x_o for R==1
        R_1    = 1._wp
        j_qc_1 = 0._wp
        !----------------------------
        ! 1) solve univariate problem
        !----------------------------
        call vqc_uncor (x_o_1,   &! <-  analysis minus observation
                        J_qc_1,  &!  -> cost function
                        dj_qc_1, &!  -> gradient of j_qc
                        w_qc_1,  &!  -> posteriori probab. for correct data
                        R_1,     &! <-  obs.err. cov.matrix for correct data
                        R_qc_1,  &!  -> inverse Hessian or Hessian
                        sgm,     &! <-  threshold for bad data (tails)
                        .TRUE.,  &! <-  return inverse matrix instead of R
                        flag,    &! <-  R modification flag
                        0,       &! <-  observation status flag
                 form = form     )! <-  formulation of obs-costfunction
        !-----------------------------------------------------------
        ! 2) rescale the pdf
        !    This corresponds to the approximation used for iset==1
        !    and not to the exact Hessian          used for iset==2.
        !    Calculations for the latter would be more complex.
        !-----------------------------------------------------------
        J_qc  = fn     * j_qc_1                 ! contribution to cost function
        w_qc  = w_qc_1                          ! probability for correct data
        dj_qc = w_qc_1 * dj_qc                  ! gradient of cost function
        if (linv) then
          R_qc = w_qc_1 * edR                   ! recip. covariance matrix
        else
          R_qc = R / w_qc_1                     ! modified cov.matrix
        endif
      endif
    endif
!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------
    subroutine set_indices
    !---------------------------------------
    ! set indices for loop over combinations
    !---------------------------------------
       n1 = i1-i0+1
       n1 = max (n1, 0)
       n2 = min (m_rej, n_lev, n1-1)
       if (n2 == 0) n1 = 0
       n2 = max (n2, 1)
       ii    =  0
       ii(1) = -1
    end subroutine set_indices
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine set_task (new_task, new_n_lev)
    !---------------------------------------------
    ! set new task and number of levels to explore
    !---------------------------------------------
    character(len=*) ,intent(in) :: new_task
    integer          ,intent(in) :: new_n_lev
       task  = new_task
       n_lev = new_n_lev
       if (ipr0 /= 0) then
         write (6,'(a,a8)') ' task    = ', task
         write (6,'(a,i8)') ' n_lev   = ', n_lev
       endif
    end subroutine set_task
!------------------------------------------------------------------------------
    subroutine calc_prob (lref)
!   use full_hessian_type_module, only: inverse
    logical ,intent(in) :: lref ! true for reference (to calculate determinant)
      !-----------------------------------------------------------------------
      ! 'calc_prob' calculates the probabilitiy of a specific combination of
      ! accepted/rejected data and accumulates the respective contributions to
      ! the cost function, its derivative and Hessian.
      !-----------------------------------------------------------------------
      real(wp)       :: det  ! determinant of inverse correlation matrix
      real(wp) ,save :: det0 ! determinant of inverse cor. (all data ok)
      integer        :: nok  ! number of correct data

      !--------------------------
      ! skip if already processed
      !--------------------------
      call proj2sig
      nm = idx ()
      if (nm /= 0) return
      !----------------------------------------------------
      ! calculate a priori probability for this combination
      !----------------------------------------------------
      nok = count (proj == 1._wp)
      i    = count (proj == 0._wp)
      wi = i
      if (rest) then
        wi = i
        i  = min  (i, m_rej + 1)
        wi = min (wi, m_rej + wrej)
      endif
      g    = apri% g_rej ** wi
!     if (rest .and. gsum>0._wp) &
!       g = gsum       ! = 1 for all remaining combinations
      gsum = gsum - g              ! remaining a priory probability
      !------------------------------------------
      ! optionally constrain a priori probability
      !------------------------------------------
      if (i /= 0) then
        if (flag > 0) then
          g   = min (g, g_rej (flag)) ! a priori probab. for error
        else
          g   = max (g, g_rej (-flag))
        endif
      endif
      !--------------------------------------------------------
      ! calculate a posteriori probability for this combination
      !--------------------------------------------------------
      edRi = inverse (C,proj=proj,det=det)     ! reciproce cov. matrix
!     if (nok == n) det0 = det
      if (lref) det0 = det
      det  = det / det0
      pi   = g * sqrt(det) *                    &!
             exp (-0.5_wp * sum ((o.x.edRi) * o))! probability of comb.
      djqc = edRi .x. o                          ! contribution to gradient
      !---------
      ! printout
      !---------
      if (ipr0 > 1) write(6,'(2i5,a,2i5,5e11.3,1x,(132l1))') &
        i0,i1,' -',mm,nok,g,sqrt(det), pi/(g*sqrt(det)),pi,gsum, proj==1._wp
      !-------------
      ! extend table
      !-------------
      call enter_table (sig,pi,i)
      !-----------------------
      ! update some quantities
      !-----------------------
      w_qc  = w_qc  + pi * proj                   ! prob. for corr. data
      p     = p     + pi                          ! probability
      p0    = p0    + g                           ! normalisation factor
      edR   = edR   + pi * (edRi-(djqc .o. djqc)) ! Hessian
      dj_qc = dj_qc + pi * djqc                   ! gradient
      wnrej(i) = wnrej(i) + pi
      nscan(i) = nscan(i) + 1
    end subroutine calc_prob
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine proj2sig
      integer :: i, i0, i1, i2
      sig = 0
      do i=1,n
        if (proj(i) == 1._wp) then
          i0 = i-1
          i2 = 1 + i0/64
          i1 = mod(i0,64)
          sig(i2) = ibset (sig(i2),i1)
        endif
      end do
    end subroutine proj2sig
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine sig2pri
      integer :: i
      do i=1,n
        pri (i) = test (sig,i-1)
      end do
    end subroutine sig2pri
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    logical function test (sig,i)
    integer(i8) ,intent(in) :: sig(:)
    integer     ,intent(in) :: i
      integer :: i1, i2
      i2 = 1 + i/64
      i1 = mod(i,64)
      test = btest (sig(i2),i1)
    end function test
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function idx ()
    integer :: idx
      integer :: i
      idx = 0
      do i=1,mm
        if (all(sigs(i,:) == sig)) idx = i
      end do
    end function idx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine best (i)
      integer ,intent(in) ,optional :: i
      if (present(i)) then
        iloc = maxloc (post(1:mm), mask = (nval(1:mm) == i))
      else
        iloc = maxloc (post(1:mm))
      endif
      sig  = sigs (iloc(1),:)
    end subroutine best
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine enter_table (sig,pi,i)
    integer(i8) ,intent(in) :: sig(:)
    real(wp)    ,intent(in) :: pi
    integer     ,intent(in) :: i
      mm = mm + 1
      sigs (mm,:) = sig
      post (mm)   = pi
      nval (mm)   = i
    end subroutine enter_table
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function komb (n,m)
    integer(i8) :: komb
    integer     :: n,m
      integer :: i
      komb = 1
      do i=m+1,n
        komb = komb * i / (i-m)
      end do
    end function komb
!------------------------------------------------------------------------------
  end subroutine vqc_corr
!==============================================================================
end module mo_vqc
