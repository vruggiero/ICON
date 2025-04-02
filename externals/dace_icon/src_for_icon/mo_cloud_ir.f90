!
!+ module to flag the presence of cloud contamination in IR sounders
!
MODULE mo_cloud_ir
!
! Description:
! Flag cloud contamination in IR sounders.
! For hyperspectral sounders a rank-sorted/model difference method is provided
! (so-called MyNally-Watts scheme). This code is outdated and supersede by the
! code in the cads* subdirectory.
! For the HIRS sounder a method developped by Olaf Stiller is provided.
!
! Heritage: mo_hyperspecsounder.f90 and mo_hirs.f90
!
! Current Code Owners:
!  DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    email: robin.faulwetter@dwd.de
!  DWD, Olaf Stiller
!    phone: +49 69 8062 2910
!    email: olaf.stiller@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_13        2011/11/01 Detlef Pingel
!  module to flag the presence of cloud contamination in AIRS/IASI
! V1_22        2013-02-13 Robin Faulwetter
!  Scientific improvement in hyperspectral clouddetection.
!  Technical improvement in lev2chan.
! V1_23        2013-03-26 Robin Faulwetter
!  Implemented processing of CrIS data
! V1_31        2014-08-21 Robin Faulwetter
!  Bugfix for level2channel assignment (l2c_type=3) and IASI crossband-flagging
! V1_35        2014-11-07 Andreas Rhodin
!  new optional parameter ct_level from subroutine cloud_detect
! V1_50        2017-01-09 Robin Faulwetter
!  Restructured/unified cloud detection for radiances.
! V1_??        2018-07-26 Kristin Raykova
!   Updated the McNally&Watts cloud detection to latest version
! V2_22        2024-07-23 Robin Faulwetter
!  Renamed from mo_hyperspecsounder to mo_cloud_is, include mo_hirs
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors: (former MODULE YOMIASI):
! Andrew Collard   2006/02/01  Original Version
! Modifications:
! Andrew Collard   2006/02/01  Original export version
! Marc Schwaerz    2009/10/28  Added real datatypes.
! Kristin Raykova  2018-07-26  Updated McNally-Watts scheme
! Olaf Stiller     2021-2017   Hirs cloud detection routines
! Robin Faulwetter 2024-07-23  Added mo_hirs to this module
!
!==============================================================================
  !=============
  ! modules used
  !=============
  use mo_kind,          only: wp
  use mo_t_use,         only: STAT_PAS_REJ,       &! report is monitored + rejected
                              STAT_REJECTED        ! rejected
  use mo_rad,           only: t_rad_set,          &
                              rad_set
  use mo_namelist,      only: position_nml,       &! position namelist
                              nnml,               &! namelist Fortran unit
                              POSITIONED           ! ok code from position_nml
  use mo_exception,     only: finish               ! abort routine
  use mo_mpi_dace,      only: dace,               &! MPI group info
                              p_bcast              ! broadcast routine
  use mo_satid,         only: satid                ! derive satellite id from name


  implicit none

  !================
  ! public entities
  !================
  private
  public :: cloud_detect_hss
  public :: cloud_detect_setup_hss
  public :: cloud_detect_type
  public :: s__cloud_detect_setup
  public :: ldebug_ass
  public :: hirs_cloud_check        ! HIRS cloud check routine (called by cloud_detect)
  public :: read_nml_hirs           ! read namelist /HIRS_CLOUD_CHECK/
  !================================
  ! Module variables and data types
  !================================

  !---------------------------------------------------------
  ! Hyperspectral sounder CADS V2.X (superseeded by CADS3.1)
  !---------------------------------------------------------
  integer, parameter :: MSENSOR_AIRS    = 11
  integer, parameter :: MSENSOR_IASI    = 16
  integer, parameter :: MSENSOR_CRIS    = 27
  integer, parameter :: MSENSOR_CRISFSR = 28
  integer, parameter :: MSENSOR_IRS     = 57
  integer, parameter :: JP__Digital_Filter = 1
  integer, parameter :: jp__min_sensor_index = MSENSOR_AIRS
  integer, parameter :: jp__max_sensor_index = MSENSOR_CRISFSR

  ! Flags set by cloud_detect_hss
  ! Bits
  integer, parameter :: MNW_CLOUD   = 0
  integer, parameter :: MNW_AEROSOL = 1
  ! flag values
  integer, parameter :: CLOUDFREE = 0
  integer, parameter :: CLOUDY    = 2 ** MNW_CLOUD

  type cloud_detect_type
     integer         :: m__sensor                   =  -1     ! Unique ID for sensor
     integer         :: n__filter_method            =  -1     ! Averaging filter for cloud detection
     integer         :: n__num_bands                =  -1     ! Number of cloud detection bands
     integer,pointer :: n__gradchkinterval(:)       => null() ! Window used in gradient calculation
     integer,pointer :: n__band_size(:)             => null() ! No. of chans in each band
     integer,pointer :: n__bands(:,:)               => null() ! List of chans in bands
     integer,pointer :: n__window_width(:)          => null() ! List of filter window widths per band
     INTEGER,pointer :: n__Window_Bounds(:,:)       => null() ! (I__Max_Bands,2)
     real(wp),pointer:: r__window_grad_threshold(:) => null() ! (I__Max_Bands)

     integer,pointer :: n__bandtouse(:)       => null()  ! Which band to use for cloud height (used in cross-band)
     logical         :: l__do_quick_exit      =  .false. ! Allow quick exit
     logical         :: l__do_crossband       =  .false. ! Cross band cloud detection
     real(wp),pointer:: r__bt_threshold(:)    => null()  ! BT threshold for cloud contamination
     real(wp),pointer:: r__grad_threshold(:)  => null()  ! Gradient threshold for cloud contamination

     !-----------------------------
     ! aerosol detection parameters
     !-----------------------------
     logical         :: l__do_aerosoldetection      = .false.
     integer         :: n__num_aerosol_tests        =  -1
     integer, pointer:: n__num_aerosol_chans(:)     => null()
     integer, pointer:: n__aerosol_chans(:,:)       => null()
     integer         :: n__mean_aerosol_chans
     real(wp),pointer:: r__aerosol_tbd(:,:)         => null()

     integer         :: update_mnw                  =  0  ! MNW version
     logical         :: use_all_chans               =  .true.
  end type cloud_detect_type

  type(cloud_detect_type) :: s__cloud_detect_setup(jp__min_sensor_index:jp__max_sensor_index)
  logical :: ldebug_ass = .false.


  !---------------------
  ! HIRS cloud detection
  !---------------------
  ! derived type definition
  integer, parameter :: m_sat  = 5  ! max # of satellites with HIRS instrument
  integer, parameter :: m_chan = 9  ! number of HIRS channel parameters

  ! parameters for HIRS cloud detection
  type t_hirs_param
    integer  :: satid               = -1        ! satellite Id
    real(wp) :: thresh_dif (m_chan) = 0._wp     ! threshold
    real(wp) :: peak_shift (m_chan) = 0._wp     ! shift of peak
    real(wp) :: sig_pr     (m_chan) = 10000._wp ! additional sigma test in ..
  end type t_hirs_param                         !   .. cloud screening

  real(wp) ,parameter :: cs_fac=1.0_wp ! CS criterion: ydiff > cs_fac*variance
  real(wp) ,parameter :: cs_faB=0.8_wp ! CS criterion: fgdep_diff(i_chan+1)>cs_faB*fgdep_diff(i_chan)

  real(wp) ,parameter :: Ts_min=274._wp ! temperature threshold for
                                        ! "warm cloud over cold ocean" mode

  type(t_hirs_param) ,save :: hirs_param (m_sat) ! parameters from namelist
  integer                  :: n_param = 0        ! number of parameter sets

  target                   :: hirs_param


contains


  !---------------------------------------------------------
  ! Hyperspectral sounder CADS V2.X (superseeded by CADS3.1)
  !---------------------------------------------------------

  subroutine construct_cdt(x)
    type(cloud_detect_type), intent(out) :: x
  end subroutine construct_cdt


  subroutine cloud_detect_hss (k__sensor,                 &
                               k__nchans,                 &
                               k__chanid,                 &
                               p__modelbts,               &
                               p__obsbts,                 &
                               p__chan_level,             &
                               k__cloud_flag,             &
                               k__minlev,                 &
                               k__maxlev,                 &
                               istat,                     &
                               procid, lastcl, ct_level,  &
                               lon, lat)

  integer, intent(in)          :: k__sensor                ! Satellite sensor (AIRS/IASI/CrIS)
  integer, intent(in)          :: k__nchans                ! Number of channels
  integer, intent(in)          :: k__chanid    (k__nchans) ! Channel indicies of input channels
  real(wp),intent(in)          :: p__modelbts  (k__nchans) ! Clear brightness temperatures from Model
  real(wp),intent(in)          :: p__obsbts    (k__nchans) ! Potentially cloudy observations
  real(wp),intent(in)          :: p__chan_level(k__nchans) ! Pressure level assignment to channel (radtrcld)
  integer, intent(out)         :: k__cloud_flag(k__nchans) ! cloud flag by channel; 0=clear, 1=cloudy
  integer, intent(in)          :: k__minlev                ! level corresponding to top of
                                                           ! initial cloud search (200hPa)
  integer, intent(in)          :: k__maxlev                ! level corresponding bottom
                                                           ! of cloud search (surface, 900hPa)
  integer, intent(out)         :: istat                    ! status
  integer, intent(in) ,optional:: procid                   !
  real(wp),intent(out),optional:: lastcl                   ! level of last clear channel
  real(wp),intent(out),optional:: ct_level                 ! RTTOV-level of cloud top
  real(wp),intent(in) ,optional:: lon                      ! longitude
  real(wp),intent(in) ,optional:: lat                      ! latitude
  !----------------------------------------------------------------------------------
  ! subroutine for cloud flagging for AIRS/IASI,
  ! using digital filter or other methods
  !
  ! This software was developed within the context of
  ! the EUMETSAT Satellite Application Facility on
  ! Numerical Weather Prediction (NWP SAF), under the
  ! Cooperation Agreement dated 25 November 1998, between
  ! EUMETSAT and the Met Office, UK, by one or more partners
  ! within the NWP SAF. The partners in the NWP SAF are
  ! the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  ! Copyright 2006, EUMETSAT, All Rights Reserved.
  !
  ! Author:
  ! Phil Watts       2002/01/21
  ! Modifications:
  ! Andrew Collard   2006/02/01  Original export version
  ! Andrew Collard   2006/05/03  Allow for missing channels
  ! Andrew Collard   2006/05/04  Allow cross-band cloud detection
  !-----------------------------------------------------------------------------
!INTF
    integer              :: ist,j,i_k                ! indices
    integer              :: jband,jband2             ! indices
    logical              :: l__aerosolpresent        !

    ! Local variables - band splitting details
    integer ,pointer     :: i__bands(:,:)            ! Channel detection bands
    integer ,pointer     :: i__band_size(:)          ! no. of chans per band
    integer ,pointer     :: i__bandtouse(:)          ! Cross band definitions
    integer              :: i__num_bands             ! no. of bands
    integer              :: i__numfoundchans         ! no. of useable chans
    integer              :: i__bandnumber(k__nchans) ! Channel band indicator
    integer ,allocatable :: index(:)
    integer ,allocatable :: idchan(:)
    integer ,allocatable :: i__cloud_flag(:)
    logical              :: l__do_crossband
    real(wp),allocatable :: z__dbt(:),z__level(:)
    real(wp),allocatable :: z__dbt_smoothed(:)
    real(wp),allocatable :: grad1(:)
    real(wp),allocatable :: grad2(:)
    integer              :: start_channel
    logical              :: logtest

    ! Local variables - digital filter parameters
    integer              :: i__chan_high          ! First high channel cloud-free a-priori
    integer              :: i__chan_low           ! Lowest channel considered
    integer              :: i__firstcloudychannel ! index of lowest cloud affected channel
    integer              :: i__last_clear_channel ! Lowest cloud affected channel
    integer              :: i__method             ! Cloud detection method
    integer ,pointer     :: i__window_width(:)
    integer              :: i__WindowBounds(2)    ! Boundary of window
    integer ,pointer     :: i__Window_Chans(:)    ! Boundary of longwave window
    real(wp)             :: z__cloud_level        ! Pressure level of cloud
    integer              :: fu
    integer ,save        :: coutzeros   = 0
    integer ,save        :: coutzerosfb = 0
    integer ,save        :: coutone     = 0
    integer              :: k__imager_flag        ! Input flag from imager data

!======================================================================

    if (present (ct_level)) ct_level = -2  ! keep -2 if band 1 is not present

    if (s__cloud_detect_setup(k__sensor)%m__sensor /= k__sensor) then
      ! not initialized
      istat = 1
      return
    else
      istat = 0
    end if

    start_channel = 0
    if ( present (lastcl)) lastcl     = 1

    ! Get correct processing parameters for this sensor:
    i__method                =  s__cloud_detect_setup(k__sensor)% n__filter_method
    i__num_bands             =  s__cloud_detect_setup(k__sensor)% n__num_bands
    i__band_size             => s__cloud_detect_setup(k__sensor)% n__band_size
    i__bands                 => s__cloud_detect_setup(k__sensor)% n__bands
    i__window_width          => s__cloud_detect_setup(k__sensor)% n__window_width
    i__bandtouse             => s__cloud_detect_setup(k__sensor)% n__bandtouse
    l__do_crossband          =  s__cloud_detect_setup(k__sensor)% l__do_crossband
    i__window_chans          => null()

    if (ldebug_ass) then
       fu = 0
       write (fu,*) 'i__method           : ', i__method, jp__digital_filter
       write (fu,*) 'i__num_bands        : ', i__num_bands
       write (fu,*) 'i__band_size        : ', i__band_size
       write (fu,*) 'i__bands            : ', i__bands
       write (fu,*) 'i__window_width     : ', i__window_width
       write (fu,*) 'i__bandtouse        : ', i__bandtouse
       write (fu,*) 'l__do_crossband     : ', l__do_crossband
       write (fu,*) 'l__gradchkinterval  : ', s__cloud_detect_setup(k__sensor)% n__gradchkinterval(:)
       write (fu,*) 'l__bt_threshold     : ', s__cloud_detect_setup(k__sensor)% r__bt_threshold   (:)
       write (fu,*) 'l__grad_threshold   : ', s__cloud_detect_setup(k__sensor)% r__grad_threshold (:)
    endif ! ldebug_ass

    ! Initialise ALL channels to cloudy
    k__cloud_flag(:) = CLOUDY

    ! If using cross-band, set up an array indicating which channels correspond
    ! to which bands in k__chanid
    if (l__do_crossband) then
       i__bandnumber(:)=-1  ! Initialise
       do jband = 1, i__num_bands
          do i_k=1,k__nchans
             if (any(i__bands(:,jband) == k__chanid(i_k))) &
                  i__bandnumber(i_k)=jband
          end do
       end do
    end if

    !*1 Loop over bands
    band_loop: do jband = 1, i__num_bands         ! Loop over bands

       ! Don't bother doing the cloud detection if we're just going to use
       ! the results from another band anyway:
       if (l__do_crossband) then
          if (.not.(any(i__bandtouse(:) == jband))) cycle
       end if

       allocate (z__dbt         (i__band_size(jband))); z__dbt         (:) = 0.0_wp
       allocate (z__dbt_smoothed(i__band_size(jband))); z__dbt_smoothed(:) = 0.0_wp
       allocate (grad1          (i__band_size(jband))); grad1          (:) = 0.0_wp
       allocate (grad2          (i__band_size(jband))); grad2          (:) = 0.0_wp
       allocate (z__level (i__band_size(jband))); z__level (:) = real(k__maxlev,wp)
       allocate (i__cloud_flag  (i__band_size(jband)))
       allocate (index          (i__band_size(jband)))
       allocate (idchan         (i__band_size(jband))); idchan(:) = 1

       i__windowbounds(:) = s__cloud_detect_setup(k__sensor)% n__window_bounds(jband,:)

       !*1.1 Find channels within current band --------------------------------------
       i__numfoundchans   = 0
       if (.not.associated(i__window_chans)) allocate(i__window_chans(2))
       i__window_chans(:) = -1

       do j=1,i__band_size(jband)
          do i_k=1,k__nchans
             if (k__chanid(i_k) == i__bands(j,jband)) then
                if (p__obsbts(i_k) < 0._wp .or. p__modelbts(i_k) < 0._wp) cycle
                i__numfoundchans = i__numfoundchans + 1
                z__dbt  (i__numfoundchans) = p__obsbts    (i_k) - p__modelbts(i_k)
                z__level(i__numfoundchans) = p__chan_level(i_k)
                index (i__numfoundchans)   = i__numfoundchans
                idchan(i__numfoundchans)   = i_k
                if (s__cloud_detect_setup(k__sensor) % update_mnw >= 1) then
                   IF (k__chanid(i_k) == i__windowbounds(1)) &
                        i__window_chans(1) = i__numfoundchans
                   IF (k__chanid(i_k) == i__windowbounds(2)) &
                        i__window_chans(2) = i__numfoundchans
                end if
             endif
          enddo
       enddo
       if (i__numfoundchans == 0) then
          if (ldebug_ass) then
             write(*,*)'**cloud_detect - WARNING: Channels not found cycling band: **',jband
          endif
          deallocate(z__dbt,z__level,i__cloud_flag,index,idchan,z__dbt_smoothed,grad1,grad2)
          cycle band_loop
       endif

       i__cloud_flag(:) = CLOUDY

       if (ldebug_ass) then
         do j = 1, i__numfoundchans
           write(fu,*) 'sort_in',j,z__level(j)
         end do
       end if

       !*2 Sort according to level assignment
       call sortrx(i__numfoundchans,z__level,index)
       if (ldebug_ass) then
         do j = 1, i__numfoundchans
           write(fu,*) 'sort_out',j,index(j),z__level(index(j))
         end do
       end if

       !*2.1 Find i__chan_low - lowest channel considered based on level assignment
       j = 1
       do while (j < i__numfoundchans .and. z__level(index(j)) < real(k__maxlev, wp))
          j = j+1
       enddo

       if (j == i__numfoundchans) then
          i__chan_low = i__numfoundchans-1
       else
          i__chan_low = j
       endif
       if (s__cloud_detect_setup(k__sensor) % update_mnw >= 1) then
          if (i__chan_low <= 1) i__chan_low = 1
       end if

       !*2.1a Find i__chan_high - highest allowed channel for starting the cloud search
       j=1
       do while (j < i__numfoundchans .and. z__level(index(j)) < real(k__minlev, wp))
          j=j+1
       enddo
       i__chan_high=j

       ! define variable needed for the x-band flagging --- supposed to equal 0 for clear sky
       if (j > 1) then
          i__firstcloudychannel=index(j-1)
       else
          i__firstcloudychannel=0
       endif

       k__imager_flag = 0 ! enetering quick exit. Otherwise no clear FOVs of window channels.

       !*3  Cloud search
       if (i__method == jp__digital_filter) then
          i__cloud_flag(:) = k__cloud_flag(idchan(:))
          call cf_digital (                   &
               k__sensor,                     &
               jband,                         &
               i__numfoundchans,              &
               z__dbt(1:i__numfoundchans),    &
               index (1:i__numfoundchans),    &
               i__chan_high,                  &
               i__chan_low,                   &
               i__window_chans,               &
               i__window_width(jband),        &
               i__cloud_flag,                 &
               i__firstcloudychannel,         &
               i__last_clear_channel,         &
               k__imager_flag,                &
               z__dbt_smoothed, grad1, grad2, &
               start_channel, logtest, fu)

          ! debug: print lat/lon
          if (ldebug_ass) then
             if ((i__firstcloudychannel == 0) .and. (jband == 1)) then
                if (present(lon) .and. present(lat)) then
                   write(fu,*) 'quick exit for profile with lon/lat: ',lon,lat
                endif
                coutzerosfb = coutzerosfb + 1
             endif
          endif

          if (present(lastcl)) then
             if (jband == 1) then
                if      (i__last_clear_channel == 0) then
                   lastcl = 0
                else if (i__last_clear_channel == 1) then
                   lastcl = 1
                else
                   lastcl = z__level(i__last_clear_channel)
                endif
             endif
          endif

          if (ldebug_ass) then
             if (i__firstcloudychannel == 0) then
                coutzeros = coutzeros + 1
             else
                if (jband == 1) then
                   if (i__last_clear_channel == 1) then
                      coutone = coutone + 1
                   endif
                endif
             endif
             write (fu,*) 'after cloud detection'
             do j = 1, i__numfoundchans
                write (fu,*)  jband,j,index(j),int(z__level(index(j))),real(i__bands(index (j),jband),wp)* &
                     0.25_wp+645._wp,z__dbt(index(j)),z__dbt_smoothed(j),grad1(index(j)),grad2 (index(j))
                if (i__firstcloudychannel /= 0) then
                   write (fu,*)  i__last_clear_channel,i__firstcloudychannel, &
                        int(z__level(i__last_clear_channel)), int(z__level(i__firstcloudychannel)), &
                        logtest,i__cloud_flag (index(j))
                else
                   write (fu,*) i__last_clear_channel,i__firstcloudychannel,logtest,i__cloud_flag(index(j))
                endif
             enddo
             write (fu,*) char (10)
          endif ! ldebug_ass

          k__cloud_flag(idchan(1:i__numfoundchans)) =  i__cloud_flag(1:i__numfoundchans)

          update_mnw1: if (s__cloud_detect_setup(k__sensor) % update_mnw >= 1) then

             ! Set cloud level for cross-band:
             if (i__firstcloudychannel == 0) then ! FOV is completely clear:
                z__cloud_level = 1.e20_wp         ! Large value
             else
               if (s__cloud_detect_setup(k__sensor) % update_mnw >= 2) then
                 z__cloud_level = z__level(i__last_clear_channel)
               else
                 z__cloud_level = z__level(i__firstcloudychannel)
               end if
             end if

             ! Automatically do cross band cloud detection for all channels
             ! (whether assigned a band or not) if JBand == 1.  This can be
             ! over-ridden for the other bands.
             if (k__sensor == msensor_iasi .and. jband == 1) &
                  where(p__chan_level(:) < z__cloud_level) k__cloud_flag(:) = CLOUDFREE
          end if update_mnw1

          crossband : if (l__do_crossband) then
             old_code: if (s__cloud_detect_setup(k__sensor) % update_mnw == 0) then
                ! set cloud level for cross-band:
                if (i__firstcloudychannel == 0) then ! FOV is completely clear:
                   z__cloud_level = 1.e20_wp         ! Large value
                else
                   z__cloud_level = z__level(i__firstcloudychannel)
                end if
             end if old_code
             ! Cross Band:
             ! Loop through bands applying cloud detection to those that take their
             ! cloud detection information from the current band JBAND.
             do jband2=1,i__num_bands
                if (i__bandtouse(jband2) == jband) then
                   where(p__chan_level(:) < z__cloud_level.and.i__bandnumber == jband2) k__cloud_flag(:) = CLOUDFREE
                end if
             end do
          end if crossband

          !-------------------------------------------
          ! set cloud top level from band 1 if present
          !-------------------------------------------
          if (present(ct_level)) then
             if (jband == 1) then
                if (i__firstcloudychannel == 0) then ! FOV is completely clear:
                   ct_level = -1.                    ! clear sky
                else
                   ct_level = z__level(i__firstcloudychannel)
                end if
             end if
          end if

       else
          ! No other methods implemented as yet
          write(6,*) 'Only cf_digital (i__method=1) is currently implemented'
       endif

       ! Deallocate arrays
       deallocate (z__dbt,z__level,i__cloud_flag,index,idchan,z__dbt_smoothed,grad1,grad2)

    enddo band_loop  ! Loop over band

    if (ldebug_ass) then
       write (fu, *) 'number of cloud free profs according to mnw:          ',coutzeros
       write (fu, *) 'number of quick exits according to mnw in first band: ',coutzerosfb
       write (fu, *) 'number of cloud bad  profs according to mnw:          ',coutone
    endif

    if (s__cloud_detect_setup(k__sensor)% l__do_aerosoldetection &
         .and. associated(s__cloud_detect_setup(k__sensor)% n__num_aerosol_chans) ) then
       call aerosol_detect (  &
            k__sensor,        &! sensor
            k__nchans,        &! no. of channels
            k__chanid,        &! Channel IDs
            p__modelbts,      &! Model clear radiance estimates
            p__obsbts,        &! Observed radiances
            l__aerosolpresent) ! flag for aerosol contamination
       if (l__aerosolpresent) &
         k__cloud_flag(:) = ibset(k__cloud_flag(:), MNW_AEROSOL)
    end if

    ! Nullify pointers
    nullify(i__band_size,             &
            i__bands,                 &
            i__window_width,          &
            i__bandtouse)

  end subroutine cloud_detect_hss


  subroutine cloud_detect_setup_hss (istat, l_aerosol)
  !------------------------------------------------------------------------------
  ! Initialise cloud detection parameters for advanced infrared sounders:
  ! Default values are set and these are overriden by namelist entries
  ! The values are assigned to the cloud detections setup structure.
  !
  !
  ! This software was developed within the context of
  ! the EUMETSAT Satellite Application Facility on
  ! Numerical Weather Prediction (NWP SAF), under the
  ! Cooperation Agreement dated 25 November 1998, between
  ! EUMETSAT and the Met Office, UK, by one or more partners
  ! within the NWP SAF. The partners in the NWP SAF are
  ! the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  ! Copyright 2006, EUMETSAT, All Rights Reserved.
  !
  ! Author:
  ! Andrew Collard  ECMWF 2006/02/01
  !
  ! Modifications
  ! Andrew Collard        2006/02/06  Original code.
  ! Andrew Collard        2006/10/19  Use IASI 300 Subset Channels.
  ! Detlef Pingel         2009        pass namelist parameters via subroutine interface
  !-------------------------------------------------------------------------------
  integer,intent(out) :: istat                          ! error return code
  logical,intent(in), optional :: l_aerosol(:)
!INTF
    integer           :: j, j__sensor                   ! loop variables

    ! Cloud detection parameters
    ! (N.B. i__max_bands must be greater than 5)
    integer,parameter :: i__max_bands         =    8    ! max. number of frequency bands
    integer,parameter :: i__max_channels      = 8461    ! max. number of channels
    integer,parameter :: i__max_aerosol_chans =  600    ! orig. 200, NWPSAF 600

    type(t_rad_set), pointer :: rs => null()
    integer                  :: i, i_nml, ii, iend, istart, iset, ierr
    logical                  :: first

    !-----------------------------
    ! namelist Cloud_Detect_Coeffs
    !-----------------------------
    integer           :: m__sensor                               ! Unique ID for sensor
    integer           :: n__filter_method                        ! Filter method
    integer           :: n__num_bands                            ! Number of frequ. cloud detection bands
    integer           :: n__gradchkinterval(i__max_bands)        ! Window used in gradient calculation
    integer           :: n__band_size(i__max_bands)              ! No. of chans in each band
    integer           :: n__bands(i__max_channels,i__max_bands)  ! List of chans in bands
    integer           :: n__window_width(i__max_bands)           ! List of filter window widths per band
    integer           :: n__window_bounds(i__max_bands,2)        ! List of filter window widths per band
    real(wp)          :: r__bt_threshold(i__max_bands)           ! BT threshold for cloud contamination
    real(wp)          :: r__window_grad_threshold(i__max_bands)  ! List of filter window widths per band
    real(wp)          :: r__grad_threshold(i__max_bands)         ! Gradient threshold for cloud contamination
    logical           :: l__do_quick_exit                        ! Allow quick exit
    logical           :: l__do_crossband                         ! Cross band cloud detection
    integer           :: n__bandtouse(i__max_bands)              ! Cross band definitions: which band to use

    !---------
    ! Aerosol
    !---------
    logical           :: l__do_aerosoldetection                               ! detect aerosols
    integer           :: n__num_aerosol_tests                                 ! no. of aerosol tests
    integer           :: n__num_aerosol_chans(i__max_bands)                   ! no. of chans for each aerosol test
    integer           :: n__aerosol_chans(i__max_aerosol_chans,i__max_bands)  ! channels used for each aerosol test
    integer           :: n__mean_aerosol_chans                                ! number of channels to compute the mean
    real(wp)          :: r__aerosol_tbd(i__max_aerosol_chans,i__max_bands)

    !------------
    ! DWD options
    !------------
    integer           :: update_mnw    = 0       ! 0  : old MNW     (v1.?)
                                                 ! 1  : updated MNW (v2.3)
    logical           :: use_all_chans = .true.

    ! namelist
    NAMELIST/Cloud_Detect_Coeffs/m__sensor, n__filter_method, n__num_bands,                     &
                                 n__gradchkinterval, n__band_size, n__bands, n__window_width,   &
                                 n__window_bounds, r__window_grad_threshold,                    &
                                 r__bt_threshold, r__grad_threshold, l__do_quick_exit,          &
                                 l__do_crossband, n__bandtouse, l__do_aerosoldetection,         &
                                 n__num_aerosol_tests, n__num_aerosol_chans, n__aerosol_chans,  &
                                 n__mean_aerosol_chans, r__aerosol_tbd, update_mnw, use_all_chans

    !============================================================================

    iStat = 0
    if (i__max_bands < 5) then
       iStat = 3
       return
    endif

    !============================================================================
    ! Loop through sensors setting up cloud detection
    !============================================================================

    SensorLoop : do J__Sensor = JP__Min_Sensor_Index, JP__Max_Sensor_Index
      if (all( (/msensor_iasi, msensor_cris, msensor_crisfsr, msensor_airs /) /= j__sensor)) cycle

      !---------------------------------------
      ! Determine namelist for current sensor
      !---------------------------------------
      if (dace% lpio) then
        i_nml = -1
        first = .true.
        j     = 0
        loop_nml: do
          j = j + 1
          call position_nml ('CLOUD_DETECT_COEFFS' ,lrewind=first, status=ierr)
          if (ierr == POSITIONED) then
            m__sensor = -1
#if defined(__ibm__)
            read (nnml ,nml=CLOUD_DETECT_COEFFS ,iostat=ierr)
            if (ierr/=0) call finish ('cloud_detect_setup_hss',               &
                 'ERROR in namelist /CLOUD_DETECT_COEFFS/')
#else
            read (nnml ,nml=CLOUD_DETECT_COEFFS)
#endif
            if (m__sensor == j__sensor) then
              i_nml = j
              exit loop_nml
            end if
          else
            exit loop_nml
          endif
          first = .false.
        end do loop_nml
      end if
      call p_bcast(update_mnw,    dace% pio)
      call p_bcast(use_all_chans, dace% pio)

      n__num_bands = 0
      n__band_size = 0

      ! Set default values
      select case (j__sensor)

      !====================
      ! Set up AIRS
      !====================
      case(msensor_airs)
        n__filter_method             = jp__digital_filter
        n__num_bands                 = 5
        if (.not.use_all_chans) then
          n__band_size(:)              = 0
          n__band_size(1:n__num_bands) = (/141, 36, 53, 25, 64/)

          if (maxval(n__band_size(:)) > i__max_channels) then
            istat = 4
            return
          endif

          n__bands(:,:)= 0

          n__bands(1:n__band_size(1),1) = &
               (/ 1,   6,   7,  10,  11,  15,  16,  17,  20,  21, &
               22,  24,  27,  28,  30,  36,  39,  40,  42,  51, &
               52,  54,  55,  56,  59,  62,  63,  68,  69,  71, &
               72,  73,  74,  75,  76,  77,  78,  79,  80,  82, &
               83,  84,  86,  92,  93,  98,  99, 101, 104, 105, &
               108, 110, 111, 113, 116, 117, 123, 124, 128, 129, &
               138, 139, 144, 145, 150, 151, 156, 157, 159, 162, &
               165, 168, 169, 170, 172, 173, 174, 175, 177, 179, &
               180, 182, 185, 186, 190, 192, 193, 198, 201, 204, &
               207, 210, 213, 215, 216, 218, 221, 224, 226, 227, &
               232, 239, 248, 250, 251, 252, 253, 256, 257, 261, &
               262, 267, 272, 295, 299, 300, 305, 308, 309, 310, &
               318, 321, 325, 333, 338, 355, 362, 375, 453, 475, &
               484, 497, 528, 587, 672, 787, 791, 843, 870, 914, &
               950 /)

          n__bands(1:n__band_size(2),2) = &
               (/ 1003, 1012, 1019, 1024, 1030, 1038, 1048, 1069, 1079, 1082,  &
               1083, 1088, 1090, 1092, 1095, 1104, 1111, 1115, 1116, 1119,  &
               1120, 1123, 1130, 1138, 1142, 1178, 1199, 1206, 1221, 1237,  &
               1252, 1260, 1263, 1266, 1278, 1285 /)

          n__bands(1:n__band_size(3),3) = &
               (/ 1290, 1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &
               1466, 1471, 1477, 1479, 1488, 1500, 1519, 1520, 1538, 1545, &
               1565, 1574, 1583, 1593, 1614, 1627, 1636, 1644, 1652, 1669, &
               1674, 1681, 1694, 1708, 1717, 1723, 1740, 1748, 1751, 1756, &
               1763, 1766, 1771, 1777, 1780, 1783, 1794, 1800, 1803, 1806, &
               1812, 1826, 1843  /)

          n__bands(1:n__band_size(4),4) = &
               (/ 1852, 1865, 1866, 1867, 1868, 1869, 1872, 1873, 1875, 1876, &
               1877, 1881, 1882, 1883, 1884, 1897, 1901, 1911, 1917, 1918, &
               1921, 1923, 1924, 1928, 1937  /)

          n__bands(1:n__band_size(5),5) = &
               (/ 1938, 1939, 1941, 1946, 1947, 1948, 1958, 1971, 1973, 1988, &
               1995, 2084, 2085, 2097, 2098, 2099, 2100, 2101, 2103, 2104, &
               2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 2115, &
               2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2128, 2134, &
               2141, 2145, 2149, 2153, 2164, 2189, 2197, 2209, 2226, 2234, &
               2280, 2318, 2321, 2325, 2328, 2333, 2339, 2348, 2353, 2355, &
               2363, 2370, 2371, 2377  /)
        end if

        n__gradchkinterval(:) = 0
        n__gradchkinterval(1:n__num_bands) = (/ 5,5,5,5,5 /)
        n__window_width   (:) = 0
        n__window_width   (1:n__num_bands) = (/ 10,6,8,5,8 /)
        r__bt_threshold   (:) = 0.
        r__bt_threshold   (1:n__num_bands) = (/ 0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp/)
        r__grad_threshold (:) = 0.
        r__grad_threshold (1:n__num_bands) = (/ 0.02_wp, 0.02_wp, 0.02_wp, 0.02_wp, 0.02_wp /)

        n__window_bounds(:,:)        = 0
        n__window_bounds(1,1)        = 475
        n__window_bounds(1,2)        = 950

        r__window_grad_threshold(:)  = 0.
        r__window_grad_threshold(1)  = 0.4

        l__do_quick_exit             = .true.

        !--------------------
        ! This is cross-band:
        !--------------------
        l__do_crossband              = .false.
        n__bandtouse(:)              = 0
        n__bandtouse(1:n__num_bands) = (/ 1,2,1,4,5 /)

        !-----------------
        ! This is aerosol:
        !-----------------
        l__do_aerosoldetection                                           = .false.
        n__num_aerosol_tests                                             = 1
        n__num_aerosol_chans(:)                                          = 0
        n__num_aerosol_chans(1:n__num_aerosol_tests)                     = (/ 4 /)
        n__aerosol_chans(:,:)                                            = 0
        n__aerosol_chans(n__num_aerosol_tests,1:n__num_aerosol_chans(1)) = (/ 1260, 914, 950, 1301 /)
        n__mean_aerosol_chans                                            = 1
        r__aerosol_tbd(:,:)                                              = 0.0
        r__aerosol_tbd(n__num_aerosol_tests,1:2)                         = (/ -0.7, -2.0 /)


      !====================
      ! Set up IASI
      !====================
      case(msensor_iasi)
        n__filter_method                = jp__digital_filter
        n__num_bands                    = 5
        if (.not.use_all_chans) then
          n__band_size(:)                 = 0
          n__band_size(1:n__num_bands) = (/ 194, 15, 121, 9, 15 /)

          if (maxval(n__band_size(:)) > i__max_channels) then
            istat = 4
            return
          endif

          n__bands(:,:) = 0

          ! "IASI 354" Subset
          ! This selection was done such, that the results are backwards compatible to
          ! the DWD McNally-Watts scheme until 2019.
          n__bands(1:n__band_size(1),1) = &
               (/  16,   38,   49,   51,   55,   57,   59,   61,   63,   66, &
                   70,   72,   74,   79,   81,   83,   85,   87,   89,   92, &
                   95,   97,   99,  101,  104,  106,  109,  111,  113,  116, &
                  119,  122,  125,  128,  131,  133,  135,  138,  141,  144, &
                  146,  148,  151,  154,  157,  159,  161,  163,  165,  167, &
                  170,  173,  176,  178,  179,  180,  183,  185,  187,  189, &
                  191,  193,  195,  197,  199,  201,  203,  205,  207,  210, &
                  212,  214,  217,  219,  222,  224,  226,  228,  230,  232, &
                  234,  236,  239,  241,  242,  243,  246,  249,  252,  254, &
                  256,  258,  260,  262,  265,  267,  269,  271,  272,  273, &
                  275,  278,  280,  282,  284,  286,  288,  290,  292,  294, &
                  296,  299,  301,  303,  306,  308,  310,  312,  314,  316, &
                  318,  320,  323,  325,  327,  329,  331,  333,  335,  341, &
                  345,  347,  350,  352,  354,  356,  358,  360,  362,  364, &
                  366,  369,  371,  373,  375,  377,  379,  381,  383,  386, &
                  389,  398,  401,  404,  407,  410,  414,  416,  426,  428, &
                  432,  434,  439,  445,  457,  515,  546,  552,  559,  566, &
                  571,  573,  646,  662,  668,  756,  867,  906,  921, 1027, &
                 1046, 1121, 1133, 1191, 1194, 1271, 1786, 1805, 1884, 1991, &
                 2019, 2094, 2119, 2213 /)

          n__bands(1:n__band_size(2),2) = &
               (/1479, 1509, 1513, 1521, 1536, 1574, 1579, 1585, 1587, 1626, &
                 1639, 1643, 1652, 1658, 1671  /)

          n__bands(1:n__band_size(3),3) = &
               (/2239, 2245, 2271, 2321, 2398, 2701, 2741, 2819, 2889, 2907, &
                 2910, 2919, 2939, 2944, 2948, 2951, 2958, 2977, 2985, 2988, &
                 2991, 2993, 3002, 3008, 3014, 3027, 3029, 3036, 3047, 3049, &
                 3053, 3058, 3064, 3069, 3087, 3093, 3098, 3105, 3107, 3110, &
                 3127, 3136, 3151, 3160, 3165, 3168, 3175, 3178, 3207, 3228, &
                 3244, 3248, 3252, 3256, 3263, 3281, 3303, 3309, 3312, 3322, &
                 3339, 3375, 3378, 3411, 3438, 3440, 3442, 3444, 3446, 3448, &
                 3450, 3452, 3454, 3458, 3467, 3476, 3484, 3491, 3497, 3499, &
                 3504, 3506, 3509, 3518, 3522, 3527, 3540, 3555, 3575, 3577, &
                 3580, 3582, 3586, 3589, 3599, 3653, 3658, 3661, 3943, 4032, &
                 5130, 5368, 5371, 5379, 5381, 5383, 5397, 5399, 5401, 5403, &
                 5405, 5455, 5480, 5483, 5485, 5492, 5502, 5507, 5509, 5517, &
                 5558 /)

          n__bands(1:n__band_size(4),4) = &
               (/5992, 5994, 6003, 6350, 6463, 6601, 6962, 6980 /)

          n__bands(1:n__band_size(5),5) = &
               (/6982, 6985, 6987, 6989, 6991, 6993, 6995, 6997, 7267, 7269, &
                 7424, 7426, 7428, 7885, 8007 /)
        end if

        n__gradchkinterval(:)                 = 0
        n__gradchkinterval(1:n__num_bands)    = (/ 5,5,5,5,5 /)
        n__window_width(:)                    = 0
        n__window_width(1:n__num_bands)       = (/ 10,6,8,5,8 /)
        r__bt_threshold(:)                    = 0.
        r__grad_threshold(:)                  = 0.
        if (update_mnw /= 0) then
          r__bt_threshold(1:n__num_bands)    = (/ 0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp/)       !NWPSAF
          r__grad_threshold(1:n__num_bands)  = (/ 0.02_wp, 0.02_wp, 0.02_wp, 0.02_wp, 0.02_wp /) !NWPSAF
        else
          r__bt_threshold(1:n__num_bands)    = (/ 0.6_wp, 1.2_wp, 1.3_wp, 0.7_wp, 0.7_wp/)       !DWD
          r__grad_threshold(1:n__num_bands)  = (/ 0.03_wp, 0.05_wp, 0.1_wp, 0.05_wp, 0.05_wp /)  !DWD
        end if

        l__do_quick_exit = .true.

        n__window_bounds(:,:) = 0
        n__window_bounds(1,1) = 573
        n__window_bounds(1,2) = 2239

        r__window_grad_threshold(:) = 0.
        r__window_grad_threshold(1) = 0.4

        !--------------------
        ! This is cross-band:
        !--------------------
        l__do_crossband = .false.
        n__bandtouse(:) = 0
        !n__bandtouse(1:n__num_bands) = (/ 1,2,3,4,5 /) !NWPSAF
        n__bandtouse(1:n__num_bands) = (/ 1,2,1,4,5 /)  !DWD

        !-----------------
        ! This is aerosol:
        !-----------------
        l__do_aerosoldetection                                           = .false.
        n__num_aerosol_tests                                             = 1
        n__num_aerosol_chans(:)                                          = 0
        n__num_aerosol_chans(1:n__num_aerosol_tests)                     = (/ 4 /)
        n__aerosol_chans(:,:)                                            = 0
        n__aerosol_chans(n__num_aerosol_tests,1:n__num_aerosol_chans(1)) = (/ 1340, 2348, 1782, 2356 /)
        n__mean_aerosol_chans                                            = 11
        r__aerosol_tbd(:,:)                                              = 0.0
        r__aerosol_tbd(n__num_aerosol_tests,1:2)                         = (/ 0.2, -1.55 /)


      !====================
      ! Set up CRIS
      !====================
      case(msensor_cris,msensor_crisfsr)
        !> \todo Fill with reasonable values
        ! set up CrIS
        n__filter_method             = jp__digital_filter
        n__num_bands = 5
        if (.not.use_all_chans) then
          n__band_size(:) = 0
          n__band_size(1:n__num_bands) =(/ 152, 30, 81, 15, 42 /)

          if (maxval(n__band_size(:)) > i__max_channels) then
            istat = 4
            return
          endif

          n__bands(:,:)= 0

          ! Use the "CRIS 320" Subset   ! NWP-SAF: Use the "CRIS 300" Subset
          n__bands(1:n__band_size(1),1) = &
               (/    1,    3,    5,    8,   10,   13,   15,   17,   18,   19, &
               20,   21,   23,   25,   27,   29,   30,   32,   34,   36, &
               39,   41,   44,   46,   49,   52,   53,   54,   57,   58, &
               59,   62,   63,   65,   67,   69,   70,   71,   72,   73, &
               75,   76,   77,   79,   80,   81,   83,   85,   87,   88, &
               89,   90,   91,   92,   94,   96,   97,   99,  101,  103, &
               104,  105,  106,  107,  109,  111,  112,  114,  115,  116, &
               118,  120,  121,  123,  125,  126,  127,  128,  129,  131, &
               133,  135,  136,  137,  138,  139,  140,  141,  142,  143, &
               145,  147,  148,  149,  150,  151,  152,  153,  154,  156, &
               158,  160,  161,  163,  165,  166,  170,  172,  173,  175, &
               177,  180,  182,  185,  191,  202,  203,  204,  206,  214, &
               216,  218,  228,  233,  235,  237,  239,  246,  254,  264, &
               281,  285,  298,  305,  310,  324,  325,  326,  328,  345, &
               360,  375,  401,  428,  444,  470,  478,  492,  498,  501, &
               529,  539  /)

          n__bands(1:n__band_size(2),2) = &
               (/  562,  564,  567,  572,  578,  584,  607,  626,  639,  641, &
               643,  645,  647,  650,  653,  655,  667,  674,  676,  678, &
               680,  682,  685,  687,  689,  691,  693,  704,  707,  709  /)

          n__bands(1:n__band_size(3),3) = &
               (/  714,  716,  718,  720,  723,  726,  728,  730,  732,  734, &
               736,  738,  740,  742,  744,  746,  748,  750,  752,  754, &
               756,  758,  760,  762,  766,  768,  770,  772,  778,  780, &
               782,  784,  786,  788,  790,  792,  794,  804,  806,  809, &
               811,  812,  813,  820,  822,  824,  828,  833,  839,  843, &
               848,  851,  853,  859,  868,  872,  874,  877,  901,  915, &
               921,  951,  975,  979,  984,  993,  997, 1011, 1014, 1016, &
               1025, 1035, 1051, 1064, 1066, 1068, 1073, 1105, 1107, 1128, &
               1133  /)

          n__bands(1:n__band_size(4),4) = &
               (/ 1150, 1154, 1156, 1158, 1160, 1162, 1164, 1166, 1168, 1170, &
               1172, 1174, 1176, 1177, 1179  /)

          n__bands(1:n__band_size(5),5) = &
               (/ 1181, 1184, 1186, 1188, 1190, 1198, 1201, 1203, 1205, 1207, &
               1210, 1212, 1214, 1217, 1219, 1222, 1224, 1227, 1229, 1231, &
               1241, 1243, 1245, 1247, 1249, 1252, 1255, 1258, 1261, 1264, &
               1267, 1270, 1273, 1276, 1279, 1282, 1285, 1288, 1291, 1294, &
               1297, 1300  /)
        end if

        n__gradchkinterval(:) = 0
        n__gradchkinterval(1:n__num_bands) = (/ 5,5,5,5,5 /)

        n__window_width(:)    = 0
        n__window_width(1:n__num_bands)    = (/ 10,6,8,5,8 /)

        r__bt_threshold(:)    = 0.
        r__bt_threshold(1:n__num_bands)    = (/ 0.5, 0.5, 0.5, 0.5, 0.5/)

        r__grad_threshold(:)  = 0.
        r__grad_threshold(1:n__num_bands)  = (/ 0.02, 0.02, 0.02, 0.02, 0.02 /)

        n__window_bounds(:,:) = 0
        n__window_bounds(1,1) = 229
        n__window_bounds(1,2) = 549

        r__window_grad_threshold(:) = 0.
        r__window_grad_threshold(1) = 0.4

        l__do_quick_exit = .true.

        !--------------------
        ! This is cross-band:
        !--------------------
        l__do_crossband = .true.

        n__bandtouse(:) = 0
        n__bandtouse(1:n__num_bands) = (/ 1,1,1,1,1 /)

        !-----------------
        ! This is aerosol:
        !-----------------
        l__do_aerosoldetection                                           = .false.
        n__num_aerosol_tests                                             = 1
        n__num_aerosol_chans(:)                                          = 0
        n__num_aerosol_chans(1:n__num_aerosol_tests)                     = (/ 4 /)
        n__aerosol_chans(:,:)                                            = 0
        n__aerosol_chans(n__num_aerosol_tests,1:n__num_aerosol_chans(1)) = (/ 730, 529, 704, 732 /)
        n__mean_aerosol_chans                                            = 5
        r__aerosol_tbd(:,:)                                              = 0.0
        r__aerosol_tbd(n__num_aerosol_tests,1:2)                         = (/ 4.0, -0.5 /)

      end select

      if (present(l_aerosol)) then
        if (lbound(l_aerosol,1) <= j__sensor .and. ubound(l_aerosol,1) >= j__sensor) then
          l__do_aerosoldetection = l_aerosol(j__sensor)
        else
          call finish('cloud_detect_setup_hss','invalid shape of l_aerosol.')
        end if
      end if


      !-----------------------------------
      ! Read  namelist for current sensor
      !-----------------------------------
      if (dace% lpio .and. i_nml > 0) then
        write(6,'(1x,"Read /CLOUD_DETECT_COEFFS/ namelist number ",I2," (sensor ",I2,")")') i_nml, j__sensor
        do i = 1, i_nml
          call position_nml ('CLOUD_DETECT_COEFFS' ,lrewind=(i==1), status=ierr)
        end do
        if (ierr == POSITIONED) then
#if defined(__ibm__)
          read (nnml ,nml=CLOUD_DETECT_COEFFS ,iostat=ierr)
          if (ierr/=0) call finish ('cloud_detect_setup_hss',               &
               'ERROR in namelist /CLOUD_DETECT_COEFFS/')
#else
          read (nnml ,nml=CLOUD_DETECT_COEFFS)
#endif
        else
          if (ierr/=0) call finish ('cloud_detect_setup_hss',               &
               'Did not find namelist /CLOUD_DETECT_COEFFS/')
        endif
      end if

      m__sensor = j__sensor

      !----------------------------------------------------------------------
      call p_bcast(n__filter_method,         dace% pio)
      call p_bcast(n__num_bands,             dace% pio)
      call p_bcast(n__band_size,             dace% pio)
      call p_bcast(n__Bands,                 dace% pio)
      call p_bcast(n__gradchkinterval,       dace% pio)
      call p_bcast(n__window_width,          dace% pio)
      call p_bcast(n__window_bounds,         dace% pio)
      call p_bcast(r__window_grad_threshold, dace% pio)
      call p_bcast(r__bt_threshold,          dace% pio)
      call p_bcast(r__grad_threshold,        dace% pio)
      call p_bcast(l__do_quick_exit,         dace% pio)
      call p_bcast(l__do_crossband,          dace% pio)
      call p_bcast(n__bandtouse,             dace% pio)
      call p_bcast(l__do_aerosoldetection,   dace% pio)
      call p_bcast(n__num_aerosol_tests,     dace% pio)
      call p_bcast(n__num_aerosol_chans,     dace% pio)
      call p_bcast(n__aerosol_chans,         dace% pio)
      call p_bcast(n__mean_aerosol_chans,    dace% pio)
      call p_bcast(r__aerosol_tbd,           dace% pio)

      if (use_all_chans) then
        !----------------------------------------------------------------------
        ! Get Band information from TOVS_OBS_CHAN namelists (stored in rad_set)
        !----------------------------------------------------------------------
        do iset = 1, size(rad_set)
          rs => rad_set(iset)
          if (rs% id < 0) cycle
          if (any(rs% instr(1:rs%n_instr) == j__sensor)) then
            do ii = 1, rs%n_instr
              if (rs% instr(ii) == j__sensor) exit
            end do
            istart = rs% o_ch_i(ii) + 1
            iend   = rs% o_ch_i(ii) + rs% n_ch_i(ii)
            n__num_bands = maxval(rs%band(istart:iend))
            do j = 1, n__num_bands
              n__band_size(j) = count(rs%band(istart:iend) == j)
              n__bands(1:n__band_size(j), j) = pack(rs% chan(istart:iend), &
                   rs%band(istart:iend) == j)
            end do
          end if
        end do
      end if

      if (maxval(n__band_size(:)) > i__max_channels) then
        istat = 4
        return
      endif

      if (dace% lpio) then
         write(6,*)
         write(6,'(a)')      ' HSS cloud detection settings:'
         write(6,'(a, i6)')  '   sensor                =',m__sensor
         write(6,'(a, i6)')  '   update_mnw            =',update_mnw
         write(6,'(a, i6)')  '   filter_method         =',n__filter_method
         write(6,'(a,5x,l1)')'   use_all_chans         =',use_all_chans
         write(6,'(a, i6)')  '   num_bands             =',n__num_bands
         write(6,'(a,8i6)')  '   band_size             =',n__band_size                 &
              (1:n__num_bands)
         do j = 1, n__num_bands
           write(6,'("   bands(",I1,")              =",*(10(1x,I4),/,26x))') j, n__bands(1:n__band_size(j), j)
         end do
         write(6,'(a,8i6)')  '   grad_chk_interval     =',n__gradchkinterval           &
              (1:n__num_bands)
         write(6,'(a,8i6)')  '   window_width          =',n__window_width              &
              (1:n__num_bands)
         write(6,'(a,8i6)')  '   window_bounds 1       =',n__window_bounds             &
              (1:n__num_bands,1)
         write(6,'(a,8i6)')  '   window_bounds 2       =',n__window_bounds             &
              (1:n__num_bands,2)
         write(6,'(a,8f6.2)')'   window_grad_threshold =',r__window_grad_threshold     &
              (1:n__num_bands)
         write(6,'(a,8f6.2)')'   bt_threshold          =',r__bt_threshold              &
              (1:n__num_bands)
         write(6,'(a,8f6.2)')'   grad_threshold        =',r__grad_threshold            &
              (1:n__num_bands)
         write(6,'(a,5x,l1)')'   do_quick_exit         =',l__do_quick_exit
         write(6,'(a,5x,l1)')'   do_crossband          =',l__do_crossband
         write(6,'(a,8i6)')  '   band_to_use           =',n__BandToUse                 &
              (1:n__num_bands)
         write(6,'(a,5x,l1)')'   do_aerosol_detection  =',l__do_AerosolDetection
      end if

      !---------------------------------------------------------------
      ! Set up the s__cloud_detect_setup structure for current sensor
      !---------------------------------------------------------------
      s__cloud_detect_setup(j__sensor)%m__sensor        = m__sensor
      s__cloud_detect_setup(j__sensor)%n__filter_method = n__filter_method
      s__cloud_detect_setup(j__sensor)%n__num_bands     = n__num_bands

      allocate(s__cloud_detect_setup(j__sensor)%n__band_size(n__num_bands))
      s__cloud_detect_setup(j__sensor)%n__band_size(:)  = n__band_size(1:n__num_bands)

      allocate(s__cloud_detect_setup(j__sensor)%n__bands(maxval(n__band_size(1:n__num_bands)), n__num_bands))
      s__cloud_detect_setup(j__sensor)%n__bands(:,:)    = 0

      do j = 1, n__num_bands
         s__cloud_detect_setup(j__sensor)%n__bands(1:n__band_size(j),j) = n__bands(1:n__band_size(j),j)
      enddo

      allocate(s__cloud_detect_setup(j__sensor)%n__window_width(n__num_bands))
      s__cloud_detect_setup(j__sensor)%n__window_width(:)          = n__window_width(1:n__num_bands)

      allocate(s__cloud_detect_setup(j__sensor)%n__window_bounds(n__num_bands,2))
      s__cloud_detect_setup(j__sensor)%n__window_bounds(:,:)       = n__window_bounds(1:n__num_bands,:)

      allocate(s__cloud_detect_setup(j__sensor)%r__window_grad_threshold(n__num_bands))
      s__cloud_detect_setup(j__sensor)%r__window_grad_threshold(:) = r__window_grad_threshold(1:n__num_bands)

      allocate(s__cloud_detect_setup(j__sensor)%r__bt_threshold(n__num_bands))
      s__cloud_detect_setup(j__sensor)%r__bt_threshold(:)          = r__bt_threshold(1:n__num_bands)

      allocate(s__cloud_detect_setup(j__sensor)%r__grad_threshold(n__num_bands))
      s__cloud_detect_setup(j__sensor)%r__grad_threshold(:)        = r__grad_threshold(1:n__num_bands)

      allocate(s__cloud_detect_setup(j__sensor)%n__gradchkinterval(n__num_bands))
      s__cloud_detect_setup(j__sensor)%n__gradchkinterval(:)       = n__gradchkinterval(1:n__num_bands)

      s__cloud_detect_setup(j__sensor)%l__do_quick_exit            = l__do_quick_exit

      !------------
      ! Cross band
      !------------
      s__cloud_detect_setup(j__sensor)%l__do_crossband = l__do_crossband
      allocate(s__cloud_detect_setup(j__sensor)%n__bandtouse(n__num_bands))
      s__cloud_detect_setup(j__sensor)%n__bandtouse(:) =  n__bandtouse(1:n__num_bands)

      !---------
      ! Aerosol
      !---------
      do_aeros: if (l__do_aerosoldetection) then
         s__cloud_detect_setup(j__sensor) % l__do_aerosoldetection  =   l__do_aerosoldetection

         s__cloud_detect_setup(j__sensor) % n__num_aerosol_tests    =   n__num_aerosol_tests

         allocate( s__cloud_detect_setup(j__sensor) % n__num_aerosol_chans(n__num_aerosol_tests) )

         s__cloud_detect_setup(j__sensor) % n__num_aerosol_chans(:) =   &
              n__num_aerosol_chans(1:n__num_aerosol_tests)

         allocate(s__cloud_detect_setup(j__sensor) % n__aerosol_chans   &
              (n__num_aerosol_tests,maxval(n__num_aerosol_chans(:))))

         s__cloud_detect_setup(j__sensor) % n__aerosol_chans(:,:)   =   0
         do j = 1, n__num_aerosol_tests
            s__cloud_detect_setup(j__sensor) %                          &
                 n__aerosol_chans(j,1:n__num_aerosol_chans(j))      =   &
                 n__aerosol_chans(j,1:n__num_aerosol_chans(j))
         enddo

         s__cloud_detect_setup(j__sensor) % n__mean_aerosol_chans   =   n__mean_aerosol_chans

         allocate(s__cloud_detect_setup(j__sensor) % r__aerosol_tbd     &
              (n__num_aerosol_tests,maxval(n__num_aerosol_chans(:))))

         s__cloud_detect_setup(j__sensor) % r__aerosol_tbd(:,:)     =   0.0
         do j = 1, n__num_aerosol_tests
            s__cloud_detect_setup(j__sensor) %                          &
                 r__aerosol_tbd(j,1:n__num_aerosol_chans(j))        =   &
                 r__aerosol_tbd(j,1:n__num_aerosol_chans(j))
         enddo
      endif do_aeros

      s__cloud_detect_setup(j__sensor) % update_mnw    = update_mnw
      s__cloud_detect_setup(j__sensor) % use_all_chans = use_all_chans

   enddo sensorloop

 end subroutine cloud_detect_setup_hss


 subroutine sortrx (n, data, index)
  integer  ,intent(in)  :: n
!#if defined(__ibm__)
  ! workaround for xlf 12.1.0.3 bug, which does not properly
  ! use/allocate a temporary array for explicit size arrays:
  real(wp) ,intent(in)  :: data(:)
!#else
  !   real(wp) ,intent(in)  :: data(n)
!#endif
  integer  ,intent(out) :: index(n)
  !-----------------------------------------------------------------
  ! From Leonard J. Moss of SLAC:
  ! Here's a hybrid QuickSort I wrote a number of years ago.  It's
  ! based on suggestions in Knuth, Volume 3, and performs much better
  ! than a pure QuickSort on short or partially ordered input arrays.
  !
  !     SORTRX -- SORT, Real input, indeX output
  !
  !
  !     Input:  N     integer
  !             DATA  REAL
  !
  !     Output: INDEX integer (DIMENSION N)
  !
  ! This routine performs an in-memory sort of the first N elements of
  ! array DATA, returning into array INDEX the indices of elements of
  ! DATA arranged in ascending order.  Thus,
  !
  !    DATA(INDEX(1)) will be the smallest number in array DATA;
  !    DATA(INDEX(N)) will be the largest number in DATA.
  !
  ! The original data is not physically rearranged.  The original order
  ! of equal input values is not necessarily preserved.
  !-----------------------------------------------------------------
  ! SORTRX uses a hybrid QuickSort algorithm, based on several
  ! suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
  ! "pivot key" [my term] for dividing each subsequence is chosen to be
  ! the median of the first, last, and middle values of the subsequence;
  ! and the QuickSort is cut off when a subsequence has 9 or fewer
  ! elements, and a straight insertion sort of the entire array is done
  ! at the end.  The result is comparable to a pure insertion sort for
  ! very short arrays, and very fast for very large arrays (of order 12
  ! micro-sec/element on the 3081K for arrays of 10K elements).  It is
  ! also not subject to the poor performance of the pure QuickSort on
  ! partially ordered data.
  !
  ! Created:  15 Jul 1986  Len Moss
  !-------------------------------------------------------------------
  !
  ! The present version is a minor rewrite in Fortran 95 and takes into
  ! account the initial order through a modified final insertion sort.
  ! It performs sufficiently well on NEC SX-9.
  !
  ! Harald Anlauf  (2009).
  !-------------------------------------------------------------------
  !intf
    integer         :: lstk(31),rstk(31),istk  ! Stack
    integer         :: l,r,i,j,p,indexp,indext ! Indices
    real (kind=wp)  :: datap                   ! Pivot element
    integer,parameter:: m=9 ! QuickSort Cutoff:
                            ! Quit QuickSort-ing when a subsequence contains M or fewer
                            ! elements and finish off at end with straight insertion sort.
                            ! According to Knuth, V.3, the optimum value of M is around 9.

    ! Make initial guess for index
!CDIR ON_ADB(INDEX)
    do i=1,n
       index(i)=i
    end do

    ! If array is short, skip QuickSort and go directly to
    ! the straight insertion sort.
    if (n <= m) goto 900

    ! QuickSort:
    ! The "Qn:"s correspond roughly to steps in Algorithm Q,
    ! Knuth, V.3, PP.116-117, modified to select the median
    ! of the first, last, and middle elements as the "pivot
    ! key" (in Knuth's notation, "K").  Also modified to leave
    ! data in place and produce an INDEX array.  To simplify
    ! comments, let DATA[I]=DATA(INDEX(I)).

    ! Q1: Initialize
    istk = 0
    l    = 1
    r    = n

!CDIR ON_ADB(INDEX)
!CDIR ON_ADB(DATA)
    qsort_loop: do

       ! Q2: Sort the subsequence DATA[L]..DATA[R].
       ! At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
       ! r > R, and L <= m <= R.  (First time through, there is no
       ! DATA for l < L or r > R.)

       ! Q2.5: Select pivot key
       ! Let the pivot, P, be the midpoint of this subsequence,
       ! P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
       ! so the corresponding DATA values are in increasing order.
       ! The pivot key, DATAP, is then DATA[P].
       p      = (l+r)/2
       indexp = index(p)
       datap  = data(indexp)
       if (data(index(l)) > datap) then
          index(p) = index(l)
          index(l) = indexp
          indexp   = index(p)
          datap    = data(indexp)
       endif
       if (datap > data(index(r))) then
          if (data(index(l)) > data(index(r))) then
             index(p) = index(l)
             index(l) = index(r)
          else
             index(p) = index(r)
          endif
          index(r) = indexp
          indexp   = index(p)
          datap    = data(indexp)
       endif

       ! Now we swap values between the right and left sides and/or
       ! move DATAP until all smaller values are on the left and all
       ! larger values are on the right.  Neither the left or right
       ! side will be internally ordered yet; however, DATAP will be
       ! in its final position.
       i = l
       j = r

       do
          ! Q3: Search for datum on left >= DATAP
          ! At this point, DATA[L] <= DATAP.  We can therefore start scanning
          ! up from L, looking for a value >= DATAP (this scan is guaranteed
          ! to terminate since we initially placed DATAP near the middle of
          ! the subsequence).
          do
             i = i+1
             if (data(index(i)) >= datap) exit
          end do

          ! Q4: Search for datum on right <= DATAP
          ! At this point, DATA[R] >= DATAP.  We can therefore start scanning
          ! down from R, looking for a value <= DATAP (this scan is guaranteed
          ! to terminate since we initially placed DATAP near the middle of
          ! the subsequence).
          do
             j = j-1
             if (data(index(j)) <= datap) exit
          end do

          ! Q5: Have the two scans collided?
          if (i >= j) exit

          ! Q6: No, interchange DATA[I] <--> DATA[J] and continue
          indext   = index(i)
          index(i) = index(j)
          index(j) = indext
       end do

       ! Q7: Yes, select next subsequence to sort
       ! At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
       ! for all L <= l < I and J < r <= R.  If both subsequences are
       ! more than M elements long, push the longer one on the stack and
       ! go back to QuickSort the shorter; if only one is more than M
       ! elements long, go back and QuickSort it; otherwise, pop a
       ! subsequence off the stack and QuickSort it.
       if (r-j >= i-l .and. i-l > m) then
          istk       = istk+1
          lstk(istk) = j+1
          rstk(istk) = r
          r = i-1
       else if (i-l > r-j .and. r-j > m) then
          istk       = istk+1
          lstk(istk) = l
          rstk(istk) = i-1
          l = j+1
       else if (r-j > m) then
          l = j+1
       else if (i-l > m) then
          r = i-1
       else
          ! Q8: Pop the stack, or terminate QuickSort if empty
          if (istk < 1) exit qsort_loop
          l    = lstk(istk)
          r    = rstk(istk)
          istk = istk-1
       endif
    end do qsort_loop
900 continue

    ! Q9: Straight Insertion sort
!CDIR ON_ADB(INDEX)
!CDIR ON_ADB(DATA)
    do i=2,n
       if (     data(index(i-1)) >  data(index(i))  &
          .or. (data(index(i-1)) == data(index(i))  &
               .and. index(i-1)  >       index(i) ) ) then
          indexp=index(i)
          datap=data(indexp)
          p=i-1
          do
             index(p+1) = index(p)
             p=p-1
             if (p > 0) then
                if (     data(index(p)) >  datap   &
                   .or. (data(index(p)) == datap   &
                        .and. index(p)  >  indexp) ) cycle
             endif
             exit
          end do
          index(p+1) = indexp
       endif
    end do
    ! All done
  end subroutine sortrx


  subroutine cf_digital(k__sensor,                     &
                        k__band,                       &
                        k__numchans,                   &
                        p__dbt,                        &
                        k__index,                      &
                        k__chan_high,                  &
                        k__chan_low,                   &
                        k__chan_windows,               &
                        k__window_width,               &
                        k__cloud_flag,                 &
                        k__cloud_level,                &
                        k__last_clear_lev,             &
                        k__imager_flag,                &
                        test,test1,test2,              &
                        start_channel,logtest,fdiscr)

  integer ,intent(in)   :: k__sensor          ! satellite sensor (AIRS/IASI/CrIS)
  integer ,intent(in)   :: k__band            ! band number
  integer ,intent(in)   :: k__numchans        ! number of useable channels in band
  real(wp),intent(in)   :: p__dbt(:)          ! input ranked dBT signal
  integer ,intent(in)   :: k__index(:)        ! ranking index for dbt
  integer ,intent(in)   :: k__chan_high       ! first channel clear of high stratospheric
                                              ! model errors
  integer ,intent(in)   :: k__chan_low        ! Low channel considered for inital minimum search
  integer ,intent(in)   :: k__window_width    ! Window width for filter
  integer ,intent(in)   :: k__chan_windows(2) ! Two channels defining bounds of
                                              ! the longwave window
  integer ,intent(inout):: k__cloud_flag(:)   ! Cloud flag by channel; 0=clear, 1=cloudy
  integer ,intent(out)  :: k__cloud_level     ! Index of first cloudy channel
  integer ,intent(out)  :: k__last_clear_lev  ! Index of last clear channel
  integer, intent(in)   :: k__imager_flag     ! Input imager cloud flag

  real(wp),intent(out)  :: test(:)            !
  real(wp),intent(out)  :: test1(:)           !
  real(wp),intent(out)  :: test2(:)           !
  integer ,intent(out)  :: start_channel      !
  logical ,intent(out)  :: logtest            !
  integer ,intent(in)   :: fdiscr             !
  !-----------------------------------------------------------------------------------
  ! Cloud flagging for airs using digital filter or chi-squared methods
  !
  ! Flags the presence or otherwise of cloud contamination in AIRS
  ! channels using a rank-sorted/model difference method. methods
  ! supported are digital filter and chi-squared filter.
  !
  ! This software was developed within the context of
  ! the EUMETSAT Satellite Application Facility on
  ! Numerical Weather Prediction (NWP SAF), under the
  ! Cooperation Agreement dated 25 November 1998, between
  ! EUMETSAT and the Met Office, UK, by one or more partners
  ! within the NWP SAF. The partners in the NWP SAF are
  ! the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  ! Copyright 2006, EUMETSAT, All Rights Reserved.
  !
  ! Author:
  ! Phil Watts         2002/01/21
  !
  ! Modifications:
  ! Andrew Collard     2006/02/03  Tidy up in preparation for IASI
  ! Andrew Collard     2006/05/03  Band size is now passed in (allows
  !                                for missing channels).
  ! Andrew Collard     2006/05/04  The index of the first cloudy channel is now
  !                                returned to allow cross-band cloud detection
  !-----------------------------------------------------------------------------------
!INTF
    real(wp),allocatable:: z__dbt_smoothed(:)       ! Smoothed-ranked DBT signal
                                                    ! (extended to allow for
                                                    ! detection algorithm).
    integer             :: jch,jmin(1),jmax(1),i
    integer             :: i__buffer                ! Number of buffer channels
    integer             :: i__start_channel         ! channel to start filtering
    integer             :: i__start_channel_surf    ! 2ndary starting channel for cloud search
    integer             :: i__max_channel           ! Channel corresponding to the maximum
                                                    ! of the smoothed dBT curve
    logical             :: llcold
    logical             :: ll__window_grad_check
    logical             :: ll__startchannelchanged
    logical             :: ll__search_for_cloud_top

    ! These carry the values in S__cloud_detect_setup
    integer             :: i__gradchkinterval       ! Used in calculating the gradient
    real(wp)            :: z__bt_threshold          ! Solution contaminated threshold
    real(wp)            :: z__grad_threshold        ! Gradient threshold at which
                                                    ! to stop filter procession
    real(wp)            :: z__window_grad_threshold ! Gradient threshold for window
    integer             :: irun1                    !
    integer, save       :: ncold  = 0               !
    integer, save       :: nwarm  = 0               !
    integer, save       :: nquick = 0               !


    ! Note: Depending on the selection in the namelist, code switches between old and updated MNW version.
    mnw_version: select case(s__cloud_detect_setup(k__sensor)% update_mnw)
    case(1,2)
       !print *, 'Going into update_mnw section.'

       test (:) = 0.0_wp
       test1(:) = 0.0_wp
       test2(:) = 0.0_wp

       i__gradchkinterval       = s__cloud_detect_setup(k__sensor)% n__gradchkinterval      (k__band)
       z__bt_threshold          = s__cloud_detect_setup(k__sensor)% r__bt_threshold         (k__band)
       z__grad_threshold        = s__cloud_detect_setup(k__sensor)% r__grad_threshold       (k__band)
       z__window_grad_threshold = s__cloud_detect_setup(k__sensor)% r__window_grad_threshold(k__band)

       k__cloud_flag(:) = CLOUDY


       !*1 Smooth with boxcar (Moving Average) filter and pad averaged array
       if (ldebug_ass) then
         do i = 1, k__numchans
           write(fdiscr,*) 'raw',i,k__index(i),p__dbt(k__index(i))
         end do
       end if

       i__buffer = i__gradchkinterval
       allocate(z__dbt_smoothed(-i__buffer:k__numchans+1))
       call movinga( p__dbt(k__index(:)), &
            k__numchans        , &
            k__window_width    , &
            z__dbt_smoothed(1:k__numchans) )
       test(1:k__numchans)            = z__dbt_smoothed(1:k__numchans)
       z__dbt_smoothed(-i__buffer:0)  = z__dbt_smoothed(1)
       z__dbt_smoothed(k__numchans+1) = z__dbt_smoothed(k__numchans)

       if (ldebug_ass) then
         do i = -i__buffer, k__numchans+1
           write(fdiscr,*) 'smoothed',i,z__dbt_smoothed(i)
         end do
       end if


       !! Detection:
       !! Find i__start_channel - channel to start filtering process.
       !! Note we always start at least one channel from the end
       !! to allow gradient calculations).

       !*2 Prepare for the cloud search

       ! First define a set of key channels

       jmin = minloc(z__dbt_smoothed(k__chan_high:k__numchans))
       i__start_channel_surf = k__chan_high+jmin(1)-1

       jmin = minloc(z__dbt_smoothed(k__chan_high:k__chan_low))
       i__start_channel      = k__chan_high+jmin(1)-1


!        if (k__numchans /= k__chan_low) &
!             print *, 'Note: k__numchans is not equal to k__chan_low ! ', k__numchans, k__chan_low
!        if (k__numchans <  k__chan_low) print *, 'k__numchans lies below k__chan_low', k__numchans, k__chan_low

       ! Look for highest channel with DBT<-BT_Threshold and move I__Start_Channel
       ! there if higher than current I__Start_Channel:
       jch = i__start_channel
       startchanloop : do i=k__chan_high,k__numchans
          if (z__dbt_smoothed(i) < -z__bt_threshold .or. i == i__start_channel) then
             jch = i
             exit startchanloop
          endif
       enddo startchanloop
       i__start_channel = jch

       ! Do the same with i__start_channel_surf
       if (i__start_channel_surf /= i__start_channel) then
          jch = i__start_channel_surf
          startchanloop_surf : do i=k__chan_high,k__numchans
             if (z__dbt_smoothed(i) < -z__bt_threshold .or. i == i__start_channel_surf) then
                jch = i
                exit startchanloop_surf
             endif
          enddo startchanloop_surf
          i__start_channel_surf = jch
       end if

       ! Find the position of the equivalent maximum departure (for quick exit test)
       jmax=maxloc(z__dbt_smoothed(k__chan_high:k__numchans))
       i__max_channel = k__chan_high+jmax(1)-1

       ! Long-wave window gradient check
       ll__window_grad_check=.TRUE.
       if (all(k__chan_windows > 0)) ll__window_grad_check = &
            &    (abs(z__dbt_smoothed(k__index(k__chan_windows(1))) - &
            &    z__dbt_smoothed(k__index(k__chan_windows(2)))) &
            &    < z__window_grad_threshold)

       ! Choose scenario to be followed
       LL__Search_for_Cloud_Top=.TRUE.
       which_scen: if ( abs(z__dbt_smoothed(i__start_channel_surf)) < z__bt_threshold .and. &
            abs(z__dbt_smoothed(i__start_channel)) < z__bt_threshold .and. &
            abs(z__dbt_smoothed(i__max_channel)) < z__bt_threshold .and. &
            abs(z__dbt_smoothed(k__numchans)) < z__bt_threshold .and. &
            ll__window_grad_check .and. &
            k__imager_flag==0 .and. &
            s__cloud_detect_setup(k__sensor) % l__do_quick_exit) then

          ! Quick exit
          if (k__band == 1) nquick = nquick + 1
          if (ldebug_ass) then
             write (fdiscr, *) 'quick exit was used.'
          endif
          k__cloud_flag(k__index(1:k__numchans)) = CLOUDFREE
          k__cloud_level                         = 0 ! Special Indicator of Clear FOV
          k__last_clear_lev                      = 0
          ll__search_for_cloud_top               = .false.
          return

       elseif (abs(z__dbt_smoothed(i__start_channel)) < z__bt_threshold .and. &
            &    z__dbt_smoothed(k__numchans) > z__bt_threshold ) then ! which_scen
          ! Warm cloud start at next-to-bottom channel (allowing one channel for
          ! gradient calculations).
          llcold = .false.
          i__start_channel = k__numchans-1
       elseif (z__dbt_smoothed(i__start_channel) < -z__bt_threshold ) then
          llcold = .true.
       elseif (z__dbt_smoothed(i__start_channel) > z__bt_threshold ) then
          llcold = .false.
       else
          llcold = .true.
       endif which_scen

       search_for_cloud_top: if (ll__search_for_cloud_top) then
          ! Either cold or warm start (but not quick exit)

          jch = i__start_channel
          start_channel = jch ! old code !

          ! If the primary starting channel appears clear, start from the
          ! secondary starting channel, if it is lower. in that case also
          ! re-evaluate cold or warm start is more appropriate.

          second_check: if (i__start_channel /= i__start_channel_surf .and. &
               s__cloud_detect_setup(k__sensor)% update_mnw >= 1) then

             ll__startchannelchanged  = .false.
             if (llcold .and. ( (z__dbt_smoothed(jch-1)-z__dbt_smoothed(jch+1)) < &
                  &       z__grad_threshold .and. &
                  &       z__dbt_smoothed(jch-i__gradchkinterval)-z__dbt_smoothed(jch+1) < &
                  &       z__grad_threshold .and. &
                  &       abs(z__dbt_smoothed(jch)) < z__bt_threshold)) then
                i__start_channel = i__start_channel_surf
                ll__startchannelchanged  = .true.
             endif

             if (ll__startchannelchanged) then

                if (abs(z__dbt_smoothed(i__start_channel)) < z__bt_threshold .and. &
                     &             z__dbt_smoothed(k__numchans) > z__bt_threshold ) then
                   ! Warm cloud start at next-to-bottom channel
                   ! (allowing one channel for gradient calculations).
                   if (k__band == 1) nwarm = nwarm + 1
                   if (ldebug_ass) then
                      write (fdiscr, *) 'Warm cloud was used.', z__dbt_smoothed(i__start_channel), &
                           z__dbt_smoothed(k__numchans),      &
                           z__bt_threshold
                   endif
                   llcold = .false.
                   i__start_channel = k__numchans-1
                elseif (z__dbt_smoothed(i__start_channel) < -z__bt_threshold ) then
                   if (k__band == 1) ncold = ncold + 1
                   if (ldebug_ass) then
                      write (fdiscr, *) 'Cold cloud was used.', z__dbt_smoothed(i__start_channel), &
                           z__dbt_smoothed(k__numchans),      &
                           z__bt_threshold
                   endif
                   llcold = .true.
                elseif (z__dbt_smoothed(i__start_channel) > z__bt_threshold ) then
                   if (k__band == 1) nwarm = nwarm + 1
                   if (ldebug_ass) then
                      write (fdiscr, *) 'Warm cloud was used.', z__dbt_smoothed(i__start_channel), &
                           z__dbt_smoothed(k__numchans),      &
                           z__bt_threshold
                   endif
                   llcold = .false.
                else
                   If (k__band == 1) ncold = ncold + 1
                   if (ldebug_ass) then
                      write (fdiscr, *) 'Cold cloud was used.', z__dbt_smoothed(i__start_channel), &
                           z__dbt_smoothed(k__numchans),      &
                           z__bt_threshold
                   endif
                   llcold = .true.
                endif
                jch = i__start_channel

             endif
          endif second_check

          if (ldebug_ass) then
             write (fdiscr, *) 'ncold :', ncold
             write (fdiscr, *) 'nwarm :', nwarm
             write (fdiscr, *) 'nquick:', nquick
             write (fdiscr, *) 'start channel:', jch, k__chan_high
          endif

          logtest = llcold
          do irun1=1,k__numchans
             test1(irun1) = z__dbt_smoothed(irun1-1)-z__dbt_smoothed(irun1+1)
             test2(irun1) = z__dbt_smoothed(irun1-i__gradchkinterval)-z__dbt_smoothed(irun1+1)
          enddo

          !*3 Search for the lowest non-contaminated channel

          lowest_clear_chan: if (llcold) then
             ! If cold then progress towards higher channels whilst -ve difference is decreasing
             do while (( (z__dbt_smoothed(jch-1)                 -z__dbt_smoothed(jch+1)) > z__grad_threshold .or. &
                  (z__dbt_smoothed(jch-i__gradchkinterval)-z__dbt_smoothed(jch+1)) > z__grad_threshold .or. &
                  abs (z__dbt_smoothed(jch)) > z__bt_threshold) .and. jch > 1 )
                if (ldebug_ass) then
                   write (fdiscr, *) 'lcold: ',jch,                                       &
                        (z__dbt_smoothed(jch-1                 )-z__dbt_smoothed(jch+1)), &
                        (z__dbt_smoothed(jch-i__gradchkinterval)-z__dbt_smoothed(jch+1)), &
                        z__grad_threshold,                                                &
                        z__dbt_smoothed(jch),                                             &
                        z__bt_threshold
                endif
                jch = jch-1
             enddo
          else
             ! If warm then progress towards higher channels whilst +ve difference is  decreasing
             do while ( ((z__dbt_smoothed(jch-1)                 -z__dbt_smoothed(jch+1)) < -1.0_wp*z__grad_threshold .or. &
                  (z__dbt_smoothed(jch-i__gradchkinterval)-z__dbt_smoothed(jch+1)) < -1.0_wp*z__grad_threshold .or. &
                  abs (z__dbt_smoothed(jch)) > z__bt_threshold) .and. jch > 1 ) ! old code: jch > k__chan_high
                jch = jch-1
             enddo
          endif lowest_clear_chan
          k__cloud_level = jch

          ! Check departure at final position is not outside threshold (dwd code)
          if (abs(z__dbt_smoothed(jch)) > z__bt_threshold) then
             k__cloud_level = 1                            ! k__chan_high
             if (ldebug_ass) then
                write (fdiscr, *) 'Dirty chan high was used bad profile.', z__bt_threshold, &
                     abs (z__dbt_smoothed(jch)), jch
             endif
          endif
          if (k__cloud_level /= 1) k__cloud_flag(k__index(1:k__cloud_level-1)) = CLOUDFREE

          ! Output channel numbers for the highest cloud and lowest clear levels
          k__last_clear_lev = 1
          if (k__cloud_level>1) then ! new code
             k__last_clear_lev = k__index(k__cloud_level-1)
          else
             k__last_clear_lev = k__index(k__cloud_level)
          endif
          k__cloud_level = k__index(k__cloud_level)

       endif search_for_cloud_top

       if (allocated(z__dbt_smoothed)) deallocate(z__dbt_smoothed)

       !----------------------------------------------------------------------------------------------------------------
    case (0) ! MWN version
       !----------------------------------------------------------------------------------------------------------------

       !print *, 'Find clouds using operational code.'

       test (:) = 0.0_wp
       test1(:) = 0.0_wp
       test2(:) = 0.0_wp

       i__gradchkinterval = s__cloud_detect_setup(k__sensor)% n__gradchkinterval(k__band)
       z__bt_threshold    = s__cloud_detect_setup(k__sensor)% r__bt_threshold   (k__band)
       z__grad_threshold  = s__cloud_detect_setup(k__sensor)% r__grad_threshold (k__band)

       ! smooth with boxcar (Moving Average) filter and pad averaged array
       i__buffer = i__gradchkinterval
       allocate(z__dbt_smoothed(-i__buffer:k__numchans+1))
       call movinga( p__dbt(k__index(:)), &
            k__numchans        , &
            k__window_width    , &
            z__dbt_smoothed(1:k__numchans) )
       test(1:k__numchans)            = z__dbt_smoothed(1:k__numchans)
       z__dbt_smoothed(-i__buffer:0)  = z__dbt_smoothed(1)
       z__dbt_smoothed(k__numchans+1) = z__dbt_smoothed(k__numchans)

       ! Detection:
       ! Find i__start_channel - channel to start filtering process.
       ! Note we always start at least one channel from the end
       ! to allow gradient calculations).
       jmin             = minloc(z__dbt_smoothed(k__chan_high:k__numchans))
       i__start_channel = k__chan_high+jmin(1)-1
       if (abs(z__dbt_smoothed(i__start_channel)) < z__bt_threshold .and. &
            z__dbt_smoothed(k__numchans     )  < z__bt_threshold .and. &
            s__cloud_detect_setup(k__sensor)% l__do_quick_exit   ) then

          ! Quick exit
          if (k__band == 1) nquick = nquick + 1
          if (ldebug_ass) then
             write (fdiscr, *) 'quick exit was used.'
          endif
          k__cloud_flag(k__index(1:k__numchans)) = CLOUDFREE
          k__cloud_level                         = 0 ! Special Indicator of Clear FOV
          k__last_clear_lev                      = 0
          return

       elseif (z__dbt_smoothed(i__start_channel) > -z__bt_threshold .and. &
            z__dbt_smoothed(k__numchans     ) >  z__bt_threshold ) then

          ! warm cloud start at next-to-bottom channel
          ! (allowing one channel for gradient calculations).
          if (k__band == 1) nwarm = nwarm + 1
          if (ldebug_ass) then
             write (fdiscr, *) 'Warm cloud was used.', z__dbt_smoothed(i__start_channel), &
                  z__dbt_smoothed(k__numchans),      &
                  z__bt_threshold
          endif
          llcold           = .false.
          i__start_channel = k__numchans
       else
          ! cold cloud start
          if (k__band == 1) ncold = ncold + 1
          if (ldebug_ass) then
             write (fdiscr, *) 'Cold cloud was used.', z__dbt_smoothed(i__start_channel), &
                  z__dbt_smoothed(k__numchans),      &
                  z__bt_threshold
          endif
          llcold              = .true.
       endif

       jch           = i__start_channel
       start_channel = jch

       if (ldebug_ass) then
          write (fdiscr, *) 'ncold :', ncold
          write (fdiscr, *) 'nwarm :', nwarm
          write (fdiscr, *) 'nquick:', nquick
          write (fdiscr, *) 'start channel:', jch, k__chan_high
       endif

       logtest = llcold
       do irun1=1,k__numchans
          test1(irun1) = z__dbt_smoothed(irun1-1)-z__dbt_smoothed(irun1+1)
          test2(irun1) = z__dbt_smoothed(irun1-i__gradchkinterval)-z__dbt_smoothed(irun1+1)
       enddo

       ! if cold then progress towards higher channels whilst -ve difference is decreasing
       if (llcold) then
          do while (( (z__dbt_smoothed(jch-1)                 -z__dbt_smoothed(jch+1)) > z__grad_threshold .or. &
               (z__dbt_smoothed(jch-i__gradchkinterval)-z__dbt_smoothed(jch+1)) > z__grad_threshold .or. &
               z__dbt_smoothed(jch) < -z__bt_threshold) .and. jch > k__chan_high )
             !orig. code:   abs (z__dbt_smoothed(jch)) > z__bt_threshold) .and. jch > k__chan_high )
             if (ldebug_ass) then
                write (fdiscr, *) 'lcold: ',jch,                                       &
                     (z__dbt_smoothed(jch-1                 )-z__dbt_smoothed(jch+1)), &
                     (z__dbt_smoothed(jch-i__gradchkinterval)-z__dbt_smoothed(jch+1)), &
                     z__grad_threshold,                                                &
                     z__dbt_smoothed(jch),                                             &
                     z__bt_threshold
             endif
             jch = jch-1
          enddo
       else
          ! if warm then progress towards higher channels whilst +ve difference is  decreasing
          do while ( ((z__dbt_smoothed(jch-1)                 -z__dbt_smoothed(jch+1)) < -1.0_wp*z__grad_threshold .or. &
               (z__dbt_smoothed(jch-i__gradchkinterval)-z__dbt_smoothed(jch+1)) < -1.0_wp*z__grad_threshold .or. &
               abs (z__dbt_smoothed(jch)) > z__bt_threshold) .and. jch > k__chan_high )
             jch = jch-1
          enddo
       endif
       k__cloud_level = jch

       ! check departure at final position is not outside threshold
       if (abs(z__dbt_smoothed(jch)) > z__bt_threshold) then
          k__cloud_level = 1                            ! k__chan_high
          if (ldebug_ass) then
             write (fdiscr, *) 'dirty chan high was used bad profile.', z__bt_threshold, &
                  abs (z__dbt_smoothed(jch)), jch
          endif
       endif
       !   k__cloud_flag(k__index(k__cloud_level:k__numchans))=1
       if (k__cloud_level /= 1) k__cloud_flag(k__index(1:k__cloud_level-1)) = CLOUDFREE

       ! output channel number for first cloudy level
       k__last_clear_lev = 1
       if (k__cloud_level /= 1) then
          k__last_clear_lev = k__index(k__cloud_level-1)
       endif
       k__cloud_level = k__index(k__cloud_level)


    end select mnw_version

  end subroutine cf_digital


  subroutine aerosol_detect (k__sensor,   &       ! in
                             k__nchans,   &       ! in
                             k__chanid,   &       ! in
                             p__modelbts, &       ! in
                             p__obsbts,   &       ! in
                             k__aerosolpresent)   ! out

  integer ,intent(in)  :: k__sensor             ! sensor
  integer ,intent(in)  :: k__nchans             ! no. of channels
  integer ,intent(in)  :: k__chanid(k__nchans)  ! channel ids
  real(wp),intent(in)  :: p__modelbts(k__nchans)! model clear radiance estimates
  real(wp),intent(in)  :: p__obsbts(k__nchans)  ! observed radiances
  logical ,intent(out) :: k__aerosolpresent
  !----------------------------------------------------------------------
  ! Look for aerosol signal in AIRS/IASI radiances
  ! A unique theoretically derived aerosol signal is sought through
  ! linear regression of first guess departures.
  ! This software was developed within the context of
  ! the EUMETSAT Satellite Application Facility on
  ! Numerical Weather Prediction (NWP SAF), under the
  ! Cooperation Agreement dated 25 November 1998, between
  ! EUMETSAT and the Met Office, UK, by one or more partners
  ! within the NWP SAF. The partners in the NWP SAF are
  ! the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  ! Copyright 2006, EUMETSAT, All Rights Reserved.
  !
  ! Author:
  ! Andrew Collard  ECMWF 2006/05/17

  ! Modifications:
  ! Andrew Collard        2006/05/17  Original code.
  ! Andrew Collard        2006/10/19  Modification to present channels test
  !------------------------------------------------------------------------
!INTF
  ! Local variables

  integer               :: j, i__tb, i__k, i__test
  integer               :: i__maxchans
  integer               :: i__num_aerosol_chans
  integer               :: i__aerosol

  integer,  pointer     :: i__aerosol_chans(:)
  integer               :: i__m
  integer               :: i__mean_aerosol_chans

  real(wp)              :: z_diff_aerosol
  real(wp), pointer     :: z__tbd(:)
  real(wp), allocatable :: z_tbm(:)

  integer,  allocatable :: i__numfoundchans_1(:)
  integer               :: i__mean2_aerosol_chans

  !-----------------------------------
  ! initialise
  !-----------------------------------

  k__aerosolpresent = .false.  !0

  i__maxchans = maxval(s__cloud_detect_setup(k__sensor) % n__num_aerosol_chans(:))

  !-----------------------------------
  ! loop through tests
  !-----------------------------------

  allocate(z_tbm(i__maxchans))
  allocate(i__numfoundchans_1(i__maxchans))
  i__aerosol=0

  ae_loop : do i__test = 1, s__cloud_detect_setup(k__sensor) % n__num_aerosol_tests

     i__num_aerosol_chans   = s__cloud_detect_setup(k__sensor) % n__num_aerosol_chans(i__test)
     i__aerosol_chans => s__cloud_detect_setup(k__sensor) % n__aerosol_chans(i__test,1:i__num_aerosol_chans)
     z__tbd => s__cloud_detect_setup(k__sensor) % r__aerosol_tbd(i__test,1:i__num_aerosol_chans)
     i__mean_aerosol_chans  = s__cloud_detect_setup(k__sensor) % n__mean_aerosol_chans
     i__mean2_aerosol_chans = int(i__mean_aerosol_chans/2)+1

     do i__tb=1, i__num_aerosol_chans
        i__numfoundchans_1(i__tb)=0
        z_tbm(i__tb)=0
     enddo

     do i__k=1,k__nchans
        if (p__obsbts(i__k) <= 0.) cycle
        do i__tb=1, i__num_aerosol_chans
           do j=1, i__mean_aerosol_chans
              i__m = i__aerosol_chans(i__tb)-i__mean2_aerosol_chans+j
              if (i__m == k__chanid(i__k)) then
                 z_tbm(i__tb) = z_tbm(i__tb)+p__obsbts(i__k)
                 i__numfoundchans_1(i__tb) = i__numfoundchans_1(i__tb)+1
              endif
           enddo
        enddo
     enddo

     do i__tb=1, i__num_aerosol_chans
        if(i__numfoundchans_1(i__tb)==0) cycle ae_loop
        z_tbm(i__tb) = z_tbm(i__tb)/i__numfoundchans_1(i__tb)
     enddo

     z_diff_aerosol=z_tbm(1)-z_tbm(2)
     if (z_diff_aerosol <= z__tbd(1)) i__aerosol=i__aerosol+1
     z_diff_aerosol=z_tbm(3)-z_tbm(4)
     if (z_diff_aerosol <= z__tbd(2)) i__aerosol=i__aerosol+1

  enddo ae_loop

  if (i__aerosol==2) then
     k__aerosolpresent = .true. !1
  end if
  if (k__nchans < 6) k__aerosolpresent = .false.

  deallocate(z_tbm,i__numfoundchans_1)
  nullify(i__aerosol_chans)
  nullify(z__tbd)


  end subroutine aerosol_detect


  subroutine movinga(pv,kv,kw,pva)
  integer ,intent(in)    :: kv      ! umber of elements in pv
  real(wp),intent(in)    :: pv(kv)  ! original array to be averaged
  integer ,intent(in)    :: kw      ! length of averaging window
  real(wp),intent(inout) :: pva(kv) ! averaged array
  !-----------------------------------------------------------
  ! Moving average of REAL array:
  ! Calculate the moving average (smoothing filter) of array
  ! No error checking supplied.
  !
  ! This software was developed within the context of
  ! the EUMETSAT Satellite Application Facility on
  ! Numerical Weather Prediction (NWP SAF), under the
  ! Cooperation Agreement dated 25 November 1998, between
  ! EUMETSAT and the Met Office, UK, by one or more partners
  ! within the NWP SAF. The partners in the NWP SAF are
  ! the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  ! Copyright 2006, EUMETSAT, All Rights Reserved.
  !
  ! Author:
  ! Phil Watts         2002/01/24
  !
  ! Modifications:
  ! Andrew Collard     2001/02/06  Original export version.
  !-----------------------------------------------------------
!INTF
    integer :: inj,j,i

    pva(:)=0
    do i = 1,kv  ! loop over array elements
       inj = 0
       do j=i-kw/2,i+kw/2,1        ! loop over window
          if (j > 0 .and. j < (kv+1)) then ! if window element exists in original array
             inj = inj+1
             pva(i) = pva(i)+pv(j) ! add value
          endif
       enddo
       pva(i)=pva(i)/real(inj, wp) ! mean value
    enddo

  end subroutine movinga

  !---------------------
  ! HIRS cloud detection
  !---------------------
  subroutine hirs_cloud_check(satid, chans, o, fg, ts, state, cloudy, ierr, debug)
  integer,       intent(in)           :: satid
  integer,       intent(in)           :: chans (:)
  real(kind=wp), intent(in)           :: o     (:)
  real(kind=wp), intent(in)           :: fg    (:)
  real(kind=wp), intent(in)           :: ts
  integer,       intent(in)           :: state (:)
  integer,       intent(out)          :: cloudy(:)
  integer,       intent(out)          :: ierr
  logical,       intent(in), optional :: debug

    integer,            parameter   :: nca = 15
    integer                         :: ncr ! Number of required channels
    integer                         :: nc
    integer                         :: i, i_n, i_p, ichan
    integer                         :: ich(nca)
    integer                         :: i_param            ! namelist parameter set index
    real(wp),           allocatable :: fgdep_pr(:)
    real(wp)                        :: sigg, ydiff2, ydiff2B
    type(t_hirs_param), pointer     :: param
    integer                         :: l_condB (8) !
    integer                         :: lev_cld(15) ! for channels 9-15: corresp channel (from 1-8) for cloud sensitivity
    logical                         :: l_debug

    ierr = 0
    if (present(debug)) then
      l_debug = debug
    else
      l_debug = .false.
    end if
    cloudy = 2 ! clear

    !---------------
    ! set parameters
    !---------------
    l_condB (1:7) = 1
    l_condB (8)   = 0

    lev_cld( :)=0
    lev_cld( 2)=3
    lev_cld( 9)=7
    lev_cld(10)=8
    lev_cld(11)=7
    lev_cld(12)=5
    lev_cld(13)=8
    lev_cld(14)=7
    lev_cld(15)=7

    param => null()
    do i_param = 1, n_param
      if (hirs_param(i_param)% satid == satid) param => hirs_param(i_param)
    end do
    if (.not. associated (param)) then
      where(chans(:) <= nca) cloudy = 3
      ierr = 1
      return
    endif

    !-----------------
    ! prepare channels
    !-----------------
    nc = size(chans)
    ich = -1
    do i = 1, nc
      if (chans(i) > 0 .and. chans(i) <= nca) ich(chans(i)) = i
    end do
    ncr = 15 ! 8 is correct !!!!!!!!
    do i = 1, ncr
      if (ich(i) <= 0) then
        where(chans(:) <= nca) cloudy = 3
        ierr = -i
        return
      end if
    end do

    allocate(fgdep_pr(nc))
    !----------------------------------------
    ! compute smoothed first guess departures
    ! (needed for revised ECMWF method)
    !----------------------------------------
    do ichan = 2, 7
      i   = ich(ichan)
      i_n = ich(ichan+1)
      i_p = ich(ichan-1)
      fgdep_pr(i) = 1.0_wp / 3.0_wp * (  o (i_p) - fg(i_p) &
                                       + o (i  ) - fg(i  ) &
                                       + o (i_n) - fg(i_n) )
    end do
    i   = ich(8)
    i_n = ich(9)
    i_p = ich(7)
    fgdep_pr(i  ) = o (i) - fg(i)
    fgdep_pr(i_n) = 2 * fgdep_pr(i) + fgdep_pr(i_p)

    !----------------------------------------
    ! perform cloud check for channels 3 to 8
    !----------------------------------------
    cloudy(:)        = 2 ! clear
    cloudy(ich(3:8)) = 1 ! cloudy

    sigg = 1._wp
    if (ts < Ts_min) then
      if (fgdep_pr(ich(8)) - fgdep_pr(ich(4)) > 0.5) sigg = -1._wp
    end if
    do ichan = 3, 8
      i   = ich(ichan)
      i_n = ich(ichan+1)
      i_p = ich(ichan-1)
      ydiff2  = fgdep_pr(i  ) - fgdep_pr(i_p)
      ydiff2B = fgdep_pr(i_n) - fgdep_pr(i  )
      if ( (sigg * (param%peak_shift(ichan) - ydiff2)                   &
            > cs_fac * param% thresh_dif(ichan)) .and.                  &
           ( l_condB(ichan) == 0 .or.                                   &
             ( sigg *          (param%peak_shift(ichan+1) - ydiff2B)    &
              >sigg * cs_fab * (param%peak_shift(ichan  ) - ydiff2 ) )) &
          ) exit
      if ( abs(fgdep_pr(i)) > param% sig_pr(ichan)) exit
      cloudy(ich(i)) = 2
    end do
    !---------------------------------------------------------
    ! set cloud flag to undetermined if channels were rejected
    !---------------------------------------------------------
    do ichan = 1, 8
      i = ich(ichan)
      if (state(i) <= STAT_PAS_REJ .or. state(i) == STAT_REJECTED) &
           cloudy(ich(max(ichan-1,1):8)) = 3
    end do
    !-----------------------------------------------------
    ! relate cloud flag of channels 9..15 to channels 3..8
    !-----------------------------------------------------
    do ichan = 9, 15
      i = ich(ichan)
      if (i > 0) then
        i_n = ich(lev_cld(ichan))
        cloudy(i) = cloudy(i_n)
      end if
    end do

  end subroutine hirs_cloud_check

!==============================================================================

  subroutine read_nml_hirs
  !--------------------------------------------------
  ! reads namelist /HIRS_CLOUD_CHECK/
  ! multiple namelist groups for different satellites
  !--------------------------------------------------

    !----------------
    ! local variables
    !----------------
#if defined(__ibm__)
    integer :: ios
#endif
    integer :: ierr
    logical :: first
    integer :: sat_id  ! WMO satellite id
    integer :: i       ! channel loop index

    !----------------------------
    ! namelist /HIRS_CLOUD_CHECK/
    !----------------------------
    character(len=8)    :: sat                  ! satellite name
    real(wp)            :: thresh_dif (m_chan)
    real(wp)            :: peak_shift (m_chan)
    real(wp)            :: sig_pr     (m_chan)

    namelist /HIRS_CLOUD_CHECK/ sat, thresh_dif, peak_shift, sig_pr

    !----------------------------------------------
    ! loop over namelist groups, i.e. satellite ids
    !----------------------------------------------
    first = .true.
    do
      if (dace% lpio) then
        !-------------------
        ! set default values
        !-------------------
        sat        = ''
        thresh_dif = 0._wp
        peak_shift = 0._wp
        sig_pr     = 10000._wp
        !--------------
        ! read namelist
        !--------------
        call position_nml ('HIRS_CLOUD_CHECK' ,lrewind=first ,status=ierr)
        first=.false.
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=HIRS_CLOUD_CHECK, iostat=ios)
          if(ios/=0)call finish('read_nml_hirs',                      &
                                'ERROR in namelist /HIRS_CLOUD_CHECK/')
#else
          read (nnml ,nml=HIRS_CLOUD_CHECK)
#endif
        end select

      endif

      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit

      if (dace% lpio) then
        !------------------------
        ! check for valid entries
        !------------------------
        sat_id = satid (sat)

        if (sat_id < 1) call finish ('read_nml_hirs',          &
                                     'invalid satellite: '//sat)
        !------------------------------
        ! store entry in satellite list
        !------------------------------
        if (n_param >= m_sat) call finish('read_nml_hirs','n_param > m_sat')
        hirs_param (n_param+1)% satid      = sat_id
        hirs_param (n_param+1)% thresh_dif = thresh_dif
        hirs_param (n_param+1)% peak_shift = peak_shift
        hirs_param (n_param+1)% sig_pr     = sig_pr
        !---------
        ! printout
        !---------
        write(6,'()')
        write(6,'(a)') repeat('-',79)
        write(6,'()')
        write(6,'(a)')     ' namelist /HIRS_CLOUD_CHECK/'
        write(6,'(a,a)')   '   sat        = ',sat
        write(6,'(a,i4)')  '   satid      =' ,sat_id
        write(6,'(a)')     '   chan_param = chan thresh_dif peak_shift sig_pr'
        do i = 1, m_chan
          write(6,'(a,i4,3f11.6)') '                ',            &
                                   i, thresh_dif(i), peak_shift(i),sig_pr(i)
        end do
        write(6,'()')
      endif

      !----------------
      ! broadcast entry
      !----------------
      n_param = n_param + 1
      call p_bcast_hirs_param (hirs_param (n_param), dace% pio)
    end do

  end subroutine read_nml_hirs

!---------------------------------------------------------
! subroutine p_bcast_hirs_param (buffer, p_source, [comm])
!---------------------------------------------------------
#define DERIVED type(t_hirs_param)
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_hirs_param
#include "p_bcast.incf"
#undef  DERIVED
#undef  p_bcast_DERIVED
end module mo_cloud_ir
