MODULE CADS31_Module

!   This software was developed within the context of the EUMETSAT
!   Satellite Application Facility on Numerical Weather Prediction
!   (NWP SAF), under the Cooperation Agreement dated 7 December 2016,
!   between EUMETSAT and the Met Office, UK, by one or more partners
!   within the NWP SAF. The partners in the NWP SAF are the Met
!   Office, ECMWF, DWD and MeteoFrance.
!
!   Copyright 2020, EUMETSAT, All Rights Reserved.

!   * CADS31_Module *
!   A. Collard  ECMWF 01/02/06

!   * PURPOSE *
!   -----------
!   Sets up structures to be used in processing of advanced IR sounders.

!   * MODIFICATIONS *
!   -----------------
!   01/02/06   A.Collard   1.0   Original export version.
!   17/11/09   R.Eresmaa   1.1   Include parameters of the Quick Exit /
!                                long-wave window gradient check.
!   11/11/11   R.Eresmaa   1.2   Add processing capability for CrIS.
!   03/12/13   R.Eresmaa   2.0   Add imager-assisted cloud detection.
!   10/11/15   R.Eresmaa   2.2   Changed instrument ID naming convention.
!                                Changed aerosol detection parameters.
!   20/12/16   R.Eresmaa   2.3   Remove aerosol detection parameters.
!   05/02/19   R.Eresmaa   2.4   Explicit KIND specifications.
!   16/04/20   R.Eresmaa   3.0   Combine cloud and aerosol detection, rename.
!                                Include aerosol type recognition.
!                                Include land sensitivity parameters.
!                                Include trace gas detection. Rename.

  use mo_kind,              only: wp, i4
  
IMPLICIT NONE

public          ! All module variables here are to be exported

SAVE

INTEGER(i4), PARAMETER :: INST_ID_AIRS = 11
INTEGER(i4), PARAMETER :: INST_ID_IASI = 16
INTEGER(i4), PARAMETER :: INST_ID_CRIS = 27
INTEGER(i4), PARAMETER :: INST_ID_CRISFSR = 28
INTEGER(i4), PARAMETER :: INST_ID_IRS = 57
INTEGER(i4), PARAMETER :: INST_ID_IASING = 59
INTEGER(i4), PARAMETER :: INST_ID_IKFS2 = 94
INTEGER(i4), PARAMETER :: INST_ID_HIRAS = 97
INTEGER(i4), PARAMETER :: INST_ID_GIIRS = 98

INTEGER(i4), PARAMETER :: JP__MIN_SENSOR_INDEX = INST_ID_AIRS
INTEGER(i4), PARAMETER :: JP__MAX_SENSOR_INDEX = INST_ID_GIIRS

TYPE Aerosol_Detect_Type
INTEGER(i4) :: M__Sensor = 0_i4                         ! Unique ID for sensor
INTEGER(i4) :: N__Num_Aerosol_Tests = 0_i4              ! Number of aerosol
                                                     ! detection tests
INTEGER(i4), POINTER :: N__Num_Regression(:) => null()     ! Number of conversion
                                                     ! coefficients for AOD
INTEGER(i4), POINTER :: N__Num_Aerosol_Chans(:) => null()  ! Number of aerosol
                                                     ! detection channels
INTEGER(i4), POINTER :: N__Aerosol_Chans(:,:) => null()    ! List of aerosol
                                                     ! detection channels
INTEGER(i4)          :: N__Mean_Aerosol_Chans = 0_i4    ! Boxcar averaging window
                                                     ! width
REAL(wp), POINTER    :: R__Aerosol_TBD(:,:) => null()      ! Aerosol detection
                                                     ! thresholds
REAL(wp), POINTER    :: R__coef_AOD(:,:) => null()         ! Coefficients for
                                                     ! conversion to AOD
REAL(wp)             :: R__Rank_Thres_Coeff(3) = 0._wp   ! Coefficients to
                                                     ! restrict rejections
                                                     ! to affected channels
REAL(wp)             :: R__Unclassified_Thres = 0._wp    ! Rejection threshold for
                                                     ! unclassified aerosol
REAL(wp)             :: R__Land_Fraction_Thres = 0._wp   ! Threshold for land
                                                     ! fraction in FOV
END TYPE Aerosol_Detect_Type

TYPE Cloud_Detect_Type
INTEGER(i4)         :: M__Sensor = 0_i4                 ! Unique ID for sensor
INTEGER(i4)         :: N__Num_Bands = 0_i4              ! Number of channel bands
INTEGER(i4), POINTER :: N__GradChkInterval(:) => null()    ! Window width used in
                                                     ! gradient calculation
INTEGER(i4), POINTER :: N__Band_Size(:) => null()          ! Number of channels in
                                                     ! each band
INTEGER(i4), POINTER :: N__Bands(:,:) => null()            ! Channel lists
INTEGER(i4), POINTER :: N__Window_Width(:) => null()       ! Smoothing filter window
                                                     ! widths per band
INTEGER(i4), POINTER :: N__Window_Bounds(:,:) => null()    ! Channels in the spectral
                                                     ! window gradient check
INTEGER(i4), POINTER :: N__BandToUse(:) => null()          ! Band number assignment
                                                     ! for each channel
LOGICAL  :: L__Do_Quick_Exit = .false.                         ! On/off switch for the
                                                     ! Quick Exit scenario
LOGICAL  :: L__Do_CrossBand = .false.                          ! On/off switch for the
                                                     ! cross-band method
REAL(wp), POINTER :: R__BT_Threshold(:) => null()          ! BT threshold for cloud
                                                     ! contamination
REAL(wp), POINTER :: R__Grad_Threshold(:) => null()        ! Gradient threshold for
                                                     ! cloud contamination
REAL(wp), POINTER :: R__Window_Grad_Threshold(:) => null() ! Threshold for window
                                                     ! gradient check in QE

LOGICAL  :: L__Do_Imager_Cloud_Detection = .false.             ! On/off switch for the
                                                     ! imager cloud detection
INTEGER(i4)         :: N__Num_Imager_Chans = 0_i4       ! No. of imager channels
INTEGER(i4)         :: N__Num_Imager_Clusters = 0_i4    ! No. of clusters to be
                                                     ! expected
INTEGER(i4),POINTER :: N__Imager_Chans(:) => null()        ! List of imager channels
REAL(wp),POINTER    :: R__Stddev_Threshold(:) => null()    ! St. Dev. threshold, one
                                                     ! for each imager channel
REAL(wp)            :: R__Coverage_Threshold = 0._wp     ! Threshold for
                                                     ! fractional coverage
                                                     ! of a cluster
REAL(wp)            :: R__FG_Departure_Threshold = 0._wp ! Threshold for imager
                                                     ! FG departure

END TYPE Cloud_Detect_Type

TYPE Land_Sensitivity_Type
INTEGER(i4)         :: M__Sensor = 0_i4                 ! Unique ID for sensor
REAL(wp)            :: R__Land_Fraction_Thres = 0._wp    ! Threshold on land
                                                     ! fraction
REAL(wp)            :: R__Level_Thres            ! Threshold on normalized
                                                     ! channel height assignment
END TYPE Land_Sensitivity_Type

TYPE Trace_Gas_Detect_Type
INTEGER(i4)         :: M__Sensor = 0_i4                  ! Unique ID for sensor
INTEGER(i4)         :: N__Num_Trace_Gas_Checks = 0_i4    ! Number of trace gases
                                                      ! to be checked
INTEGER(i4),POINTER :: N__Num_Tracer_Channels(:) => null()  ! Number of gas-sensitive
                                                      ! channels
INTEGER(i4),POINTER :: N__Tracer_Channels(:,:) => null()    ! Gas-sensitive channels
INTEGER(i4),POINTER :: N__Num_Control_Channels(:) => null() ! Number of control
                                                      ! channels
INTEGER(i4),POINTER :: N__Control_Channels(:,:) => null()   ! Control channels
INTEGER(i4),POINTER :: N__Num_Flagged_Channels(:) => null() ! Number of affected
                                                      ! channels
INTEGER(i4),POINTER :: N__Flagged_Channels(:,:) => null()   ! Affected channels

REAL(wp),POINTER    :: R__D_Obs_Threshold(:) => null()      ! Observed Tb difference
                                                      ! threshold
REAL(wp),POINTER    :: R__D_Dep_Threshold(:) => null()      ! Departure difference
                                                      ! threshold
END TYPE Trace_Gas_Detect_Type


TYPE(Aerosol_Detect_Type), target :: &
&  S__CADS_Setup_Aerosol(JP__Min_Sensor_Index:JP__Max_Sensor_Index)

TYPE(Cloud_Detect_Type), target :: &
&  S__CADS_Setup_Cloud(JP__Min_Sensor_Index:JP__Max_Sensor_Index)

TYPE(Land_Sensitivity_Type), target :: &
&  S__CADS_Setup_Land(JP__Min_Sensor_Index:JP__Max_Sensor_Index)

TYPE(Trace_Gas_Detect_Type), target :: &
&  S__CADS_Setup_Trace_Gas(JP__Min_Sensor_Index:JP__Max_Sensor_Index)


END MODULE CADS31_Module
