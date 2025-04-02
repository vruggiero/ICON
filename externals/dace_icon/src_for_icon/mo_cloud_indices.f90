!
!+ cloud index calculation routines for radiance observations
!
! $Id$
!
MODULE mo_cloud_indices
!
! Description:
! This module provides routines for calculating cloud indices for various satellite
! sounding instruments.
!
! Heritage: mo_surfaceemissivity.f90 and mo_clouddetection.f90
!
! Current Code Owner: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_13        2011/11/01 Detlef Pingel
!  routines regarding surface types and emissivities
! V_??         2012-09-04 Andreas Messer
!  Extend for ATMS, code cleanup, remove unneeded square roots
! V1_22        2013-02-13 Andreas Messer
!  cloud-detection: Add support for ATMS (based on AMSU-A)
! V1_28        2014/02/26 Andreas Rhodin
!  clean up to suppress cray warnings on un-initialised variables
! V2_22        2024/07/23 Robin Faulwetter
!  renamed to mo_cloud_indices.f90. Added Routines
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!  Authors:
!  Christina K"opken   DWD   2003-07-29   Initial release as surface emissivity module for amsu-a.
!  Reinhold  Hess      DWD   2007-04-26   Remove bug in surface interpretation
!  Reinhold  Hess      DWD   2008-11-12   Migration to HPC/Sun Studio Compiler with minor changes
!  Christina K"opken   DWD   2009-10-15   renaming of the module to surfaceEmissivity
!                                         for later usage for other spectral regions, too.
!                                         Adding of module to cloud detection repository.
!  Robin Faulwetter    DWD   2013-2024    Modifications in heritage modules
!                            2024-07-23   Renamed module, added Routines
!=======================================================================

  !=============
  ! Modules used
  !=============

  use mo_kind, only : wp

  implicit none
  private
  !================
  ! public entities
  !================
  public :: lwp_KaSiRu94
  public :: si_amsua
  public :: LWP_QinZou2016
  public :: TPW_QinZou2016
  public :: L_index_QinZou2016
  public :: L_index_WuLiQin2021
  public :: mw_emiss               ! determin surf. type, perform rain/cloud det. for AMSUA
  public :: inv_ind
  public :: sfchl_amsua            ! surface channels AMSU-A
  public :: sfchl_atms             ! surface channels ATMS
  public :: lwp_rain_thresh        ! rain detection threshold
  public :: ldeb
!=======================================================================
  !============================================
  ! Module parameters, variables and data types
  !============================================

  real(wp), parameter :: inv_ind = huge(0._wp)

  ! general private module parameter
  logical, save      :: init_surfemiss = .false.

  ! Thresholds
  real(wp), save     :: lwp_rain_thresh = 0.3_wp

  ! Channels and Frequencies for AMSU-A, ATMS
  integer, parameter :: sfchl_amsua(4) = (/1,2,3,15/)
  real(wp),parameter :: freq_amsua (4) = (/23.8_wp,31.4_wp,50.3_wp,89.0_wp/)
  integer, parameter :: sfchl_atms (4) = (/1,2,3,16/)
  real(wp),parameter :: freq_atms  (4) = (/23.8_wp,31.4_wp,50.3_wp,88.20_wp/)
  integer, parameter :: sfchl_mwts3(4) = (/1,2,3,18/)
  real(wp),parameter :: freq_mwts3 (4) = (/23.8_wp,31.4_wp,50.3_wp,89.0_wp/)

  ! parameters for the emissivity model ofGrody, 1988, IEEE TGARS, Vol. 26, 850-859.
  ! (Implementation similar to ECMWF, P. Bauer)
  integer, parameter :: nsfc        = 9  ! number of surface types in Grody 1988
  ! character(len=14), parameter:: s_typstr(nsfc)= (/& ! possible surface types
  !                                 'water         ',&
  !                                 'dry land      ',&
  !                                 'wet land      ',&
  !                                 'melting snow  ',&
  !                                 'dry snow      ',&
  !                                 'refrozen snow ',&
  !                                 'new ice       ',&
  !                                 '2nd year ice  ',&
  !                                 'multi-yr ice  '/)
  ! an, bn: coefficients for equation (6):
  ! used for isfc = 1,2,3,7 (water, dry land, wet land, new ice)
  real(wp), parameter :: an(nsfc) = (/ 0.061_wp,  0.950_wp, 0.282_wp,  0.950_wp,             &
                                       1.173_wp,  2.018_wp, 0.950_wp,  1.040_wp,  1.243_wp/)
  real(wp), parameter :: bn(nsfc) = (/ 0.274_wp,  0.000_wp, 0.292_wp,  0.000_wp,             &
                                      -0.230_wp, -0.844_wp, 0.000_wp, -0.107_wp, -0.310_wp/)
  ! e0, ei, w0, k: coefficients for equation (5):
  ! used for isfc = 4,5,6   (wet,dry,refrozen snow)
  ! used for isfc = 8,9     (2nd year, multiyear ice)
  real(wp), parameter :: e0(nsfc)=(/0.0_wp,0.0_wp,0.0_wp,0.76_wp,0.90_wp,0.97_wp,0.95_wp,0.93_wp,0.92_wp/)
  real(wp), parameter :: ei(nsfc)=(/0.0_wp,0.0_wp,0.0_wp,0.99_wp,0.75_wp,0.53_wp,0.95_wp,0.83_wp,0.64_wp/)
  real(wp), parameter :: w0(nsfc)=(/1.0_wp,1.0_wp,1.0_wp,9.0_wp,33.0_wp,32.0_wp,1.0_wp,31.0_wp, 31.0_wp/)
  real(wp), parameter :: k(nsfc) =(/0.0_wp,0.0_wp,0.0_wp,2.0_wp,3.0_wp,4.0_wp,0.0_wp,2.0_wp,2.0_wp/)
  ! variable for the emissivity model ofGrody, 1988, IEEE TGARS, Vol. 26, 850-859.
  ! (Implementation similar to ECMWF, P. Bauer), modelled surface emissivity:
  real(wp), target, save :: e_mod_amsua(size(freq_amsua),nsfc)
  real(wp), target, save :: e_mod_atms (size(freq_atms),nsfc)
  real(wp), target, save :: e_mod_mwts3(size(freq_atms),nsfc)

  logical :: ldeb = .false.

contains

  pure function lwp_KaSiRu94(obs) result(lwp)
    real(wp)             :: lwp
    real(wp), intent(in) :: obs(2)
    real(wp),       parameter   :: c1    =   4.299300_wp
    real(wp),       parameter   :: c2    =  -1.406920_wp
    real(wp),       parameter   :: c3    =   0.399635_wp
    real(wp),       parameter   :: t_off = 280.000000_wp  ! from KaSiRu94
!    real(wp),      parameter    :: c1 = 1.52*2.82, c2 = -1.52, c3 = 1.52*0.35, t_off = 290    ! from WeGro94
    !-----------------------------------------------------------------------------------
    ! Approximation of the liquid water path for microwave data, see Karstens et al.,
    ! Remote sensing of cloud liquid water, 1994
    ! LWP = c1 + c2 * log( t_off - obs(37 GHz V) ) + c3 * log( t_off - obs(22 GHz V) )
    !-----------------------------------------------------------------------------------
    if ( all( obs(:) < t_off ) ) then
       lwp = c1 + c2 * log(t_off - obs(2)) + c3 * log(t_off - obs(1))
    else
       lwp  = 100._wp
    end if

  end function lwp_KaSiRu94

  ! SI for AMSU-A (Zhu, Derber, Dee, Collard 2014)
  pure function si_amsua(obs) result(si)
    ! Requires channels 1, 2, and 15 from amsua
    real(kind=wp)             :: si
    real(kind=wp), intent(in) :: obs(:)

    si = -113.2_wp + (2.41_wp - 0.0049_wp * obs(1)) * obs(1) + 0.454_wp * obs(2) - obs(3)

  end function si_amsua

  pure function LWP_QinZou2016(obs, fg, l_sea) result(lwp)
    ! Qin Zhengkun and Zou Xiaolei, 2016: Development and initial assessment of a new land index
    ! for microwave humidity sounder cloud detection. J. Meteor. Res., 30(1), 012–037,
    !  doi: 10.1007/s13351-016-5076-4
    real(kind=wp) :: lwp
    real(kind=wp), intent(in) :: obs(2)
    real(kind=wp), intent(in) :: fg(2)
    logical,       intent(in) :: l_sea

    real(wp) ,parameter     :: small          =  0.01_wp

    if (l_sea) then
      if (obs(2) >= 300._wp - small) then
        lwp = 9._wp
      else
        lwp = 0.13_wp * (obs(1) - fg(1) - 33.58_wp * &
             (obs(2)-fg(2)) / (300._wp-obs(2)))
      end if
    else
      lwp = 0.85_wp * (obs(1) - fg(1)) - (obs(2) - fg(2))
    end if
    lwp = max(lwp, 0._wp)
  end function LWP_QinZou2016

  pure function TPW_QinZou2016(obs, fg, l_sea) result(tpw)
    ! See citation for LWP_QinZou2016
    real(kind=wp) :: tpw
    real(kind=wp), intent(in) :: obs(2)
    real(kind=wp), intent(in) :: fg(2)
    logical,       intent(in) :: l_sea

    real(kind=wp) :: lwp

    lwp = LWP_QinZou2016(obs, fg, l_sea)
    tpw = ( ( (obs(1)-fg(1))-7.5_wp*lwp) * 0.1_wp)**2 - lwp**2

  end function TPW_QinZou2016

  pure function L_index_QinZou2016(obs) result(L_index)
    ! See citation for LWP_QinZou2016
    real(kind=wp) :: L_index
    real(kind=wp), intent(in) :: obs(:)

    integer       :: n
    real(kind=wp) :: mu, sigma, t_norm
    real(kind=wp) :: x_num, x_denom

    L_index = inv_ind
    n = size(obs)
    if (n < 5) return

    mu    = sum(obs(1:n)) / (1.*n)
    sigma = sqrt(sum((obs(1:n)-mu)**2) / (1.*n))

    T_norm = (obs(1) - mu) / sigma

    x_num   = 2._wp * T_norm
    x_denom = (obs(2)/100._wp - 1._wp)**3

    if (abs(x_denom) > 1.E-10) then
      L_index = x_num / x_denom
    end if

    L_index = min(max(L_index, -100._wp),100._wp)

  end function L_index_QinZou2016


  pure function L_index_WuLiQin2021(obs) result(L_index)
    !> Cloud index from
    !! Wu, Z.; Li, J.; Qin, Z.
    !! Development and Evaluation of a New Method for AMSU-A Cloud
    !! Detection over Land. Remote Sens. 2021, 13, 3646.
    !! https://doi.org/10.3390/rs13183646
    !! Expects amsua channels 1-4,15 (in this order!)
    real(kind=wp) :: L_index
    real(kind=wp), intent(in) :: obs(:)

    integer       :: n
    real(kind=wp) :: mu, sigma, t_norm

    L_index = inv_ind
    n = size(obs)
    if (n /= 5) return

    mu    = sum(obs(1:n)) / (1.*n)
    sigma = sqrt(sum((obs(1:n)-mu)**2) / (1.*n))

    T_norm = (obs(3) - mu) / sigma
    L_index = 10._wp * T_norm / exp((obs(5)-200._wp)/50._wp)
    L_index = max(L_index, 0._wp)

  end function L_index_WuLiQin2021


  subroutine ini_surface_emissivity()
    ! initialization routine for the surface emissivity module
    if (.not.init_surfemiss) then
       call mw_e_grody(e_mod_amsua, freq_amsua)
       call mw_e_grody(e_mod_atms,  freq_atms)
       call mw_e_grody(e_mod_mwts3, freq_mwts3)
       init_surfemiss = .true.
    endif
  end subroutine ini_surface_emissivity

  subroutine mw_e_grody(e_mod, freq)
    real(wp), intent(out) :: e_mod(:,:)
    real(wp), intent(in)  :: freq(:)
    !   Subroutine calculating emissivities for 9 surface types using the model of
    !   Grody, 1988, IEEE TGARS, Vol. 26, 850-859.
    !   (Implementation similar to ECMWF, P. Bauer)
    !   Surface types:  1: WATER
    !                   2: DRY LAND
    !                   3: WET LAND
    !                   4: MELTING SNOW
    !                   5: DRY SNOW
    !                   6: REFROZEN SNOW
    !                   7: NEW ICE
    !                   8: 2ND YEAR ICE
    !                   9: MULTI-YEAR ICE
    !
    !   isfc=1,2,3,7 (water,dry land,wet land,new ice): equation (6)
    !   isfc=4,5,6   (wet,dry,refrozen snow)
    !          8,9   (2nd year, multiyear ice)        : equation (5)
    integer :: isfc  ! looping indices for surface

    e_mod(:,:) = 0.0

    do isfc=1,size(e_mod,2)
      if (isfc == 4 .or. isfc == 5 .or. isfc == 6 .or. &
          isfc == 8 .or. isfc == 9) then
         e_mod(:,isfc)=                                        &
              (e0(isfc)+ei(isfc)*(freq(:)/w0(isfc))** k(isfc)) &
                           /(1.0+(freq(:)/w0(isfc))** k(isfc))
      else
         e_mod(:,isfc)=an(isfc)+bn(isfc)*log10(freq(:))
      endif
    enddo
  end subroutine mw_e_grody

  subroutine mw_emiss(instr, channum,tb_bcor,zenith,frac_ls,tskin,&
                      ierr,id_sfc,e_ret, ld, lwp)
    integer, intent(in)            :: instr       ! the instrument (3=AMSU-A, 19=ATMS)
    integer, intent(in)            :: channum(:)  ! number of input frequencies
    real(wp),intent(in)            :: tb_bcor(:)  ! observed TBs per channel at obs point
    real(wp),intent(in)            :: zenith      ! zenith angle for obs point
    real(wp),intent(in)            :: frac_ls     ! land/sea fraction for obs point
    real(wp),intent(in)            :: tskin       ! skin temperature from model for obs point
    integer ,intent(out)           :: ierr        ! error code
    integer ,intent(out), optional :: id_sfc      ! surface type as determined from AMSU obs
    real(wp),intent(out), optional :: e_ret(4)    ! retrieved surface emissivity
    logical ,intent(in),  optional :: ld
    real(wp),intent(out), optional :: lwp
    ! Subroutine for determining a surface type and performing rain/cloud
    ! detection for AMSUA using window channels 1,2,3,15.
    !     Surface types:  1: water
    !                     2: dry land
    !                     3: wet land
    !                     4: melting snow
    !                     5: dry snow
    !                     6: refrozen snow
    !                     7: new ice
    !                     8: 2nd year ice
    !                     9: multi-year ice
    !                    10: clouds over water
    !                    11: rain over water or land
    !                  [ 12: coast (frac_ls used, to be done at higher level) ]
    ! Remarks : id_sfc=1: Over open water: FASTEM-2 should be used for computing
    !                     emissivities for replacing e_ret.
    !           frac_ls >0.01 (< 0.5): coast, data/emissivities should not be used
    real(wp),parameter :: pi_180 = 3.1415927/180.
    logical ,parameter :: lweng =.true.    ! switch for groody/weng land emis
    real(wp)           :: tb_win(4)        ! observed TBs in selected window channels
    real(wp)           :: ctheta           ! cosine of zenith angle
    real(wp)           :: frac_ice         ! estimated sea ice fraction for obs point
    real(wp)           :: lwp_             ! estimated liquid water path
    real(wp)           :: e_ret_(4)        ! retrieved surface emissivity
    integer            :: id_sfc_
    real(wp),pointer   :: e_mod(:,:)       ! points to current e_mod

    if (present(ld)) then
      ldeb = ld
    else
      ldeb = .false.
    end if
    if (present(lwp   )) lwp_   = inv_ind
    if (present(id_sfc)) id_sfc = 11

    ierr     = 0

    call ini_surface_emissivity()

    ! Select instrument
    select case(instr)
    case (3)
      call prepare_tb(sfchl_amsua)
      e_mod  => e_mod_amsua
    case (19)
      call prepare_tb(sfchl_atms)
      e_mod  => e_mod_atms
    case (132)
      call prepare_tb(sfchl_mwts3)
      e_mod  => e_mod_mwts3
    case default
      ierr = 1
    end select
    if (ldeb) print*,'debug_spot emiss_cloud mw_emis tb_win',tb_win,ierr

    if (ierr /= 0) return

    ctheta = cos(zenith * pi_180)

    ! reset diagnostic parameters
    frac_ice =     0.0

    ! Ocean: - set emissivity ocean
    !        - estimate sea ice fraction and sea ice type
    !          (new, second or multiyear ice)
    !        - sea ice: check for snow and snow type (wet,dry,refrozen)
    !        - emissivity e_ret updated according to ice and snow type found
    if (frac_ls >= 0.01 .and. frac_ls < 0.5) then
       id_sfc_ = 12                    ! coast
       e_ret_(:) = -9999.99
    elseif (frac_ls < 0.01) then
       id_sfc_ = 1                     ! ocean
       e_ret_(:) = e_mod(:,id_sfc_)
       call mw_seaice_ec()

       ! check snow type if enough scattering (ice and not new ice)
       if (frac_ice>0.5 .and. id_sfc_ > 7) &
         call mw_snow_ec()
    else
       !  Land: - estimate land type (wet,dry)
       !        - choice to use Ferraro emissivities for land
       !        - check for snow and snow type (wet,dry,refrozen)
       id_sfc_ = 2
       e_ret_(:) = e_mod(:,id_sfc_)
       call mw_land()
       call mw_snow()
    endif

    if (ldeb) print*,'debug_spot emiss_cloud mw_emis frac',frac_ls,frac_ice,id_sfc_
    !  Cloud / rain check:  -> surface types 10/11
    !                       - ocean: retrieve liquid water path
    !                       - land:  check for scattering
    if (frac_ls < 0.5 .and. frac_ice <= 0.5) then
       !CK IF (id_sfc_ == 1) THEN
       call mw_sea_lwp()
       if (lwp_ > 0.1) then
          id_sfc_ = 10    ! clouds over water
       endif
       if (lwp_ > lwp_rain_thresh) then
          id_sfc_ = 11    ! rain over water
       endif
       if (ldeb) print*,'debug_spot emiss_cloud mw_emis lwp',lwp_,id_sfc_
    elseif (id_sfc_ ==2 .or. id_sfc_ == 3) then       ! rain over snow free land only
       ! TODO: Check if tb_win(4) this is still valid for atms!!!
       if (tb_win(1) - tb_win(4) > 3.0) id_sfc_ = 11 ! rain over land
    endif

    if (present(lwp   )) lwp    = lwp_
    if (present(e_ret )) e_ret  = e_ret_
    if (present(id_sfc)) id_sfc = id_sfc_

  contains

    subroutine prepare_tb(sfchl)
      integer, intent(in) :: sfchl(:)
      integer :: i
      ! return -channel-number if required channel is missing
      do i=1,size(sfchl)
         if (.not.any(channum(:)==sfchl(i))) THEN
            ierr = - sfchl(i)
            return
         endif
      enddo
      ! copy needed channels into local storage
      do i=1,size(channum)
         where (sfchl(:) == channum(i)) tb_win(:) = tb_bcor(i)
      enddo
    end subroutine

    subroutine mw_seaice_ec()
      !   Subroutine to determine sea ice fraction and sea ice type, based on
      !   observed window channel MW TBs and modelled emissivities for water and
      !   ice surfaces.
      integer  ::  ifreq, isfc
      real(wp) ::  e_obs
      real(wp) ::  dtb, dtb_min
      real(wp) ::  a,b,c,d

      ! Sea ice fraction:  Sea ice fraction retrieved from local (retrieved)
      !                    emissivity vs. open water emissivity (at 23.8 ghz)
      a =  1.7340e+00 - 0.6236e+00 * ctheta
      b =  6.9906e-03 + 2.5512e-03 * ctheta
      c = -1.0666e-03
      d = -9.0883e-03

      ! TODO: Check if tb_win(4) this is still valid for atms!!!
      e_obs = a + b*tb_win(2) + c*tb_win(1) + d*tb_win(3)

      ! sea ice fraction: estimate by comparing retrieved emissivity with
      !                   expected open water and fresh sea ice emissivities
      frac_ice = (e_obs-e_mod(1,1)) / (e_mod(1,7)-e_mod(1,1))
      frac_ice = max(0.0_wp,min(1.0_wp,frac_ice))

      ! Sea ice type: compare which sea ice type (emissivity 'spectrum')
      !               fits best the observed TBs in window channels
      !               (given the estimated sea ice fraction)
      if (frac_ice > 0.5) then
         dtb_min = huge(1._wp)
         do isfc=7,9
            dtb = 0.0
            ! TODO: Check if tb_win(4) this is still valid for atms!!!
            do ifreq=1,size(e_mod,1)
               e_obs=e_mod(ifreq,1)+frac_ice*(e_mod(ifreq,isfc)-e_mod(ifreq,1))
               dtb  =dtb+(tb_win(ifreq)-tskin*e_obs)** 2
            enddo
            if (dtb < dtb_min .or. isfc == 7) then
               dtb_min = dtb
               id_sfc_  = isfc
            endif
         enddo
         e_ret_(:) = e_mod(:,1)+frac_ice*(e_mod(:,id_sfc_)-e_mod(:,1))
      endif

      return
    end subroutine mw_seaice_ec

    subroutine mw_snow_ec()
      ! Subroutine to determine snow type, based on
      ! observed window channel MW TBs and modelled emissivities for different
      ! snow types.
      integer   ::   isfc
      real (wp) ::   dtb, dtb_min

      ! Check for scattering
      ! TODO: Check if tb_win(4) this is still valid for atms!!!
      if ((tb_win (1) - tb_win (2) > 1.0 .and. tb_win (4) < 230.0) .or. &
          (tb_win (1) - tb_win (4) > 1.0 .and. tb_win (1) < 258.0)) then
         ! Snow type: compare which snow type (emissivity 'spectrum')
         !             fits best the observed TBs in window channels
         dtb_min = huge(1._wp)
         do isfc = 4, 6
            ! TODO: Check if tb_win(4) this is still valid for atms!!!
            dtb = sum((tb_win(:) - tskin * e_mod(:,isfc)) ** 2)
            if (dtb < dtb_min .or. isfc == 4) then
               dtb_min = dtb
               id_sfc_  = isfc
            endif
         enddo

         ! TODO: This check is redundant and will be always true!!
         if (frac_ice > 0.0) then ! emissivity over ocean (sea ice with snow)
            e_ret_(:) = e_mod(:,1) + frac_ice * (e_mod(:,id_sfc_) - e_mod(:,1))
         else ! emissivity over land
            e_ret_(:) = e_mod(:,id_sfc_)
         endif
      endif
    end subroutine mw_snow_ec

    subroutine mw_land()
    !-------------------------------------------------------------------------------
    ! Subroutine to determine land type, based on
    ! observed window channel MW TBs and modelled emissivities.
    !-------------------------------------------------------------------------------
      integer  :: isfc
      real(wp) :: dtb, dtb_min

      ! Land type : Compare which land type (emissivity 'spectrum')
      !             fits best the observed TBs in window channels
      dtb_min = huge(1._wp)
      do isfc=2,3
         dtb = sum((tb_win(:) - tskin * e_mod(:,isfc)) ** 2)
         if (dtb < dtb_min .or. isfc == 2) then
            dtb_min = dtb
            id_sfc_  = isfc
         endif
      enddo

      e_ret_(:) = e_mod(:,id_sfc_)

      if (lweng) then
         ! using f. weng's 23.8, 31.4, 50.3 ghz land surface emissivities
         e_ret_(1) =  7.344616e-01 + 7.651617e-04 * tb_win(1) &
                                   + 4.906263e-03 * tb_win(2) &
                                   - 4.967451e-03 * tb_win(3)
         e_ret_(2) =  6.540201e-01 - 1.602319e-03 * tb_win(1) &
                                   + 6.560924e-03 * tb_win(2) &
                                   - 3.925177e-03 * tb_win(3)
         e_ret_(3) = -4.119000e-02 - 9.091600e-03 * tb_win(1) &
                                   + 1.217200e-02 * tb_win(2) &
                                   + 4.885100e-04 * tb_win(3)
      endif
    end subroutine mw_land

    subroutine mw_snow()
      ! Subroutine to determine snow type, based on
      ! observed window channel MW TBs and modelled emissivities for different
      ! snow types.
      integer             :: isfc
      real(wp)            :: dtb, dtb_min

      ! Check for scattering
      if ((tb_win(1) - tb_win(2) > 1.0 .and. tb_win(4) < 230.0) .or. &
          (tb_win(1) - tb_win(4) > 1.0 .and. tb_win(1) < 258.0)) then
         ! Snow type: compare which snow type (emissivity 'spectrum') fits best
         ! the observed TBs in window channels
         dtb_min = huge(1._wp)
         do isfc = 4, 6
            dtb = sum((tb_win(:) - tskin * e_mod(:,isfc)) ** 2)
            if (dtb < dtb_min .or. isfc == 4) then
               dtb_min = dtb
               id_sfc_  = isfc
            endif
         enddo
         e_ret_(:) = e_mod(:,id_sfc_)
      endif
    end subroutine mw_snow

    subroutine mw_sea_lwp()
      ! Subroutine to determine liquid water path over (ice free) ocean
      real(wp) :: d0, d1, d2

      ! TODO: test parameters from ECMWF Tech. Memo.
      d0 =  8.240 - (2.622 - 1.846 * ctheta) * ctheta
      d1 =  0.754
      d2 = -2.265

      if (tb_win(1) < 285.0 .and. tb_win(2) < 285.0) then
        lwp_ = ctheta * (d0 + d1 * log(285.0 - tb_win(1)) &
               + d2 * log(285.0 - tb_win(2)))
        lwp_ = max(lwp_, 0._wp)
      else
        lwp_ = 0.0
      endif
    end subroutine mw_sea_lwp

  end subroutine mw_emiss


end module mo_cloud_indices


