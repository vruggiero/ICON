!+ Data module for the variables only used to compute the local information
!-------------------------------------------------------------------------------

MODULE data_obs_qc_limits

!-------------------------------------------------------------------------------
!
! Description:
!   This module contains the parameters used to derive the limits (thresholds)
!   for the quality control ('QC') of conventional observations.
!   This includes variance tables deduced from observation errors, which are
!   used as variances of the background departure (for pressure dependent
!   quality control thresholds), and which may vary from the observation error
!   variances given in module 'data_obs_cdfin'.
!
!   Note: This module belongs to a group of modules for the 'COSMO observation
!   ----  operators including quality control for conventional observations'.
!         It is shared between the COSMO model and the DACE program package !
!
! Current Code Owner (for COSMO and for DACE):
!  DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V4_28        2013/07/12 Christoph Schraff
!  Initial release, by extraction from 'data_nudge_local.f90'.
!  Previous milestones in 'data_nudge_local.f90' related to QC thresholds:
!   - 1998 : Initial release.
!   - 2000 : Radiosonde height thresholds reduced to values commonly used for
!            US-network.   Height error correlation matrix added.
!   - 2001 : Variance / threshold tables for wind adjusted to values used at
!            ECMWF.   Introduction of several additional factors to thresholds
!            for selected observations following ECMWF, and to surface pressure
!            thresholds used for the spatial consistency check.
!   - 2004-2006 : New limits / errors / factors for revised optional humidity
!                 thresholds.   Revised basic humidity QC thresholds (= new
!                 GME-OI values).   Introduction of analysis layers and standard
!                 layers for revised optional multi-level check.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! @VERSION@    @DATE@     Christoph Schraff
!  - Set parameters and variables used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - Unused variables removed (or commented out).
!  - All 'iintegers' removed and replaced by standard integers.
!  - Use of 'kind_parameters' instead of 'data_parameters'.
!
! ==============================================================================
! DACE history:
! ------------
!  1.51       2017-02-24 Andreas Rhodin : update shared modules from COSMO 5.04d
! ==============================================================================
! CAUTION: This module is used by both the DACE and COSMO model codes.       !!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
          ! 1.1  pressure levels of variance tables
PUBLIC :: nqclev , tabqclp, tabqcp , nqcalv , tabanp , nanstd ,                &
          ! 1.2  error variance tables  + 1.3  error correlations
          qcvsond, qczsond, qctsond, qcvairp, qctairp, qczcorl, nzcorl,        &
          ! 1.4  further parameters for quality control
          qctf   , qctfsu , qcvfz  , qcvairl, qcvairf, qctairf, qcvpilf,       &
          qcvdrbf, qczdrbf, qcftend, qcfctnd, qcfgood, qcfpbad, qcqinvf,       &
          qcqlapf, qcqlapl, qcqlapi, qcqfge1, qcqfge2, qcqlatl, qcqfgef

! for nudging in COSMO model code only:
PUBLIC :: qcfspat, qcftsiv, qcfbias, dtchkps

!===============================================================================

! Local Declarations:

!-------------------------------------------------------------------------------
! Global (i.e. public) Declarations:
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Section 1:  Parameters used for the threshold quality control of observations
!-------------------------------------------------------------------------------

!      1.1  Length and pressure levels used for the variance tables
!           -------------------------------------------------------

  INTEGER        , PARAMETER  :: &
    nqclev  =  15    ! number of levels in the quality control threshold tables

  REAL (KIND=wp) :: &
    tabqclp(nqclev)    ! ln(tabqcp(15))

  REAL (KIND=wp) , PARAMETER  :: &
    tabqcp (nqclev)  & ! levels of the error / threshold tables
              = (/ 100000._wp, 85000._wp, 70000._wp, 50000._wp,&
                    40000._wp, 30000._wp, 25000._wp, 20000._wp,&
                    15000._wp, 10000._wp,  7000._wp,  5000._wp,&
                     3000._wp,  2000._wp,  1000._wp /)

  INTEGER        , PARAMETER  :: &
    nqcalv  =  25    ! number of analysis layers for multi-level check plus 1

  REAL (KIND=wp) , PARAMETER  :: &
    tabanp (nqcalv)  & ! pressure limits of analysis layers for multi-lev. check
              = (/ 109900._wp,100000._wp, 97500._wp, 95000._wp,&
                    92500._wp, 90000._wp, 87500._wp, 85000._wp,&
                    82500._wp, 80000._wp, 77500._wp, 72500._wp,&
                    60000._wp, 45000._wp, 35000._wp, 27500._wp,&
                    22500._wp, 17500._wp, 12500._wp,  8500._wp,&
                     6000._wp,  4000._wp,  2500._wp,  1500._wp,&
                        1._wp  /)

  INTEGER        , PARAMETER  :: &
    nanstd (nqclev)  & ! number of analysis layers within the standard layers
                       ! [ps,925], (925,775], (775,600], (600,450], (450,p-next]
              = (/ 4, 6, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)


!      1.2  Tables of square roots of background departure variances
!      -------------------------------------------------------------

!           defined at the standard levels tabqcp(1:nqclev)

  REAL (KIND=wp) , PARAMETER  :: &
                     ! (root of) radiosonde (TEMP, PILOT) wind error variance
                     !           -----------------------------
    qcvsond (nqclev) = (/ 2.3_wp,   2.3_wp,   2.5_wp,   3.0_wp,&
                          3.5_wp,   3.7_wp,   3.5_wp,   3.5_wp,&
                          3.4_wp,   3.3_wp,   3.2_wp,   3.2_wp,&
                          3.3_wp,   3.6_wp,   4.5_wp /)           ,&
                     ! (values in O.I.: 2.0, 2.4, 2.5, 3.4, 3.6, 3.8, 3.2, 3.2
                     !                  2.4, 2.2, 2.0, 2.0, 2.5, 3.0, 4.0 )
                     !
                     ! (root of) radiosonde (TEMP, PILOT) height error variance
                     !           -------------------------------
    qczsond (nqclev) = (/ 4.3_wp,   4.4_wp,   5.2_wp,   8.4_wp,&
                          9.8_wp,  10.7_wp,  11.8_wp,  13.2_wp,&
                         15.2_wp,  18.1_wp,  19.5_wp,  22.5_wp,&
                         25.0_wp,  32.0_wp,  40.0_wp /)           ,&
                     ! ( values for European Radiosondes in O.I.: 1.5*table )
                     !
                     ! (root of) radiosonde temperature error variance
                     !           ----------------------
    qctsond (nqclev) = (/ 1.2_wp,   1.0_wp,   0.7_wp,   0.4_wp,&
                          0.4_wp,   0.5_wp,   0.5_wp,   0.6_wp,&
                          0.7_wp,   0.8_wp,   0.8_wp,   0.9_wp,&
                          0.9_wp,   1.0_wp,   1.2_wp /)           ,&
                     ! ( ECMWF: 1.7, 1.5, 1.3, 1.2, 1.2, 1.4, 1.5, 1.5,
                     !          1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.5 )
                     !
                     ! (root of) AIREP wind error variance
                     !           ----------
    qcvairp (nqclev) = (/ 2.5_wp,   2.5_wp,   3.0_wp,   3.5_wp,&
                          4.0_wp,   4.0_wp,   4.0_wp,   4.0_wp,&
                          4.0_wp,   4.0_wp,   4.0_wp,   4.0_wp,&
                          4.0_wp,   4.0_wp,   4.0_wp /)           ,&
                     ! (values in the past: 3.0 below 500 hPa, 4.0 above 300hPa)
                     !
                     ! (root of) AIREP temperature error variance
                     !           -----------------
    qctairp (nqclev) = (/ 1.2_wp,   1.0_wp,   0.7_wp,   0.5_wp,&
                          0.5_wp,   0.6_wp,   0.6_wp,   0.7_wp,&
                          0.8_wp,   0.9_wp,   1.0_wp,   1.1_wp,&
                          1.1_wp,   1.2_wp,   1.4_wp /)
                     ! ( ECMWF: 1.4, 1.3, 1.2, 1.2, 1.2, 1.3, 1.3, 1.4,
                     !          1.4, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2 )
                     !
                     ! (root of) SYNOP wind error variance (not used for quality
                     !           ----------                 control thresholds)
!   qcvsynp (nqclev) = (/ 3.6_wp,   3.6_wp,   5.8_wp,   6.8_wp,&
!                         7.8_wp,   9.8_wp,  11.0_wp,  11.8_wp,&
!                        11.8_wp,  11.8_wp,  11.8_wp,  11.8_wp,&
!                        11.8_wp,  11.8_wp,  11.8_wp /)           ,&
                     ! ( ECMWF: 3.0, 3.0, 3.0, 3.4, 3.6, 3.8, 3.2, 3.2,
                     !          2.4, 2.2, 2.0, 2.0, 2.0, 2.5, 3.0
                     !          used for sea stations only )
                     !
                     ! (root of) SYNOP height error variance (not used for qual-
                     !           ------------                 ity ctrl thresh.)
!   qczsynp (nqclev) = (/ 7.0_wp,   8.0_wp,   8.6_wp,  12.1_wp,&
!                        14.9_wp,  18.8_wp,  25.4_wp,  27.7_wp,&
!                        32.4_wp,  39.4_wp,  50.3_wp,  59.3_wp,&
!                        69.8_wp,  96.0_wp, 114.2_wp /)           ,&
                     !
!   qcvdrib =  5.4_wp  ,& ! DRIBU wind    \   (root of) error variances
!   qczdrib = 14.0_wp  ,& ! DRIBU height   >    (not used for quality
!   qczship = 14.0_wp     ! SHIP height   /      control thresholds)


!      1.3  Error correlations used for the height and thickness check
!      ---------------------------------------------------------------

!           defined at the standard levels tabqcp(1:nqclev)

  REAL (KIND=wp) :: &
    qczcorl(nqclev,nqclev)      ! radiosonde height error correlation matrix

  INTEGER        , PARAMETER  :: &
    nzcorl (nqclev,nqclev) =  & ! radiosonde height error correlation matrix
                                ! ----------------------------------- [*1000]
                             RESHAPE(                                          &
 (/1000, 716, 276,  29,   5,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
    716,1000, 733, 183,  55,   9,   3,   1,   0,   0,   0,   0,   0,   0,   0, &
    276, 733,1000, 573, 268,  77,  31,  10,   2,   0,   0,   0,   0,   0,   0, &
     29, 183, 573,1000, 851, 480, 288, 138,  49,  11,   3,   1,   0,   0,   0, &
      5,  55, 268, 851,1000, 814, 601, 364, 167,  51,  19,   7,   1,   0,   0, &
      1,   9,  77, 480, 814,1000, 935, 738, 458, 200,  93,  44,  11,   2,   0, &
      0,   3,  31, 288, 601, 935,1000, 919, 678, 361, 194, 104,  31,   7,   0, &
      0,   1,  10, 138, 364, 738, 919,1000, 895, 597, 375, 229,  84,  24,   2, &
      0,   0,   2,  49, 167, 458, 678, 895,1000, 861, 649, 460, 214,  78,   9, &
      0,   0,   0,  11,  51, 200, 361, 597, 861,1000, 929, 782, 480, 230,  42, &
      0,   0,   0,   3,  19,  93, 194, 375, 649, 929,1000, 951, 710, 412, 103, &
      0,   0,   0,   1,   7,  44, 104, 229, 460, 782, 951,1000, 878, 598, 192, &
      0,   0,   0,   0,   1,  11,  31,  84, 214, 480, 710, 878,1000, 881, 426, &
      0,   0,   0,   0,   0,   2,   7,  24,  78, 230, 412, 598, 881,1000, 725, &
      0,   0,   0,   0,   0,   0,   0,   2,   9,  42, 103, 192, 426, 725,1000/)&
                                                           , (/nqclev,nqclev/) )


!      1.4  Further parameters used for quality control (QC)
!      ------------------------------------------------

  REAL (KIND=wp) :: &
                   ! time factors (upper-air / surface-level) for QC thresholds
    qctf  (4) = (/ 0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp/),& ! upperair
    qctfsu(4) = (/ 0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp/)   ! surface-
                   ! fraction of the QC threshold at the obs. time,   !   level
                   ! which is added to the threshold for each hour
                   ! of difference between observation and model time;
                   ! indices: 1: horiz. wind, 2: 'surface' pressure,
                   !          3: temperature, 4: (relative) humidity
                   !

  REAL (KIND=wp) , PARAMETER  :: &
                   !    QC check multiples for height thresholds
    qcvfz (2) = (/ 5.48_wp , 7.75_wp /)
                   ! [  threshold  =  qcvfz(jflag) * err(obs)
                   !  = SQRT(err(f.g.)**2 + err(obs)**2) * SQRT(errlim(jflag))
                   !  ~ SQRT(3) *err(obs) * SQRT(errlim(jflag))
                   !  If the deviation exceeds the threshold resulting from
                   !    errlim(1) = 10 , then: 'probably good', check thickness,
                   !                           use as limit for thickness check
                   !    errlim(2) = 20 , then: 'probably bad', reject obs.
                   !    (cf. OI code, var. 'RMULGF') ]
                   !
!   qcvfdz         ! QC check multiple for thickness threshold
!                  ! [  errlim = 5.0625  (cf. OI code, var. 'RMULTF') ]

  REAL (KIND=wp) , PARAMETER  :: &
    qcvairl =  5.0_wp  ,& ! limit of background wind speed for rejection of
                          ! observed AIREP zero winds
                   ! additional factors used for the following QC thresholds:
    qcvairf =  0.8_wp  ,& ! (ECMWF: 0.5)  AIREP wind        thresholds
    qctairf =  1.0_wp  ,& ! (ECMWF: 1.6)  AIREP temperature thresholds
    qcvpilf =  0.9_wp  ,& ! (ECMWF: 0.8)  PILOT wind        thresholds
    qcvdrbf =  0.5_wp  ,& !               DRIBU wind        thresholds
    qczdrbf =  0.9_wp  ,& !               DRIBU pressure    thresholds
                   ! factors to derive surface pressure QC thresholds
                   ! (only) if surface pressure tendency is observed:
    qcftend =  0.6_wp  ,& ! fraction of observed tendency to be added
                          ! to pressure QC threshold
    qcfctnd =  0.8_wp  ,& ! reduction factor to the pressure threshold
                          ! as specified by namelist input
                   ! factors to derive surface pressure thresholds used for the
                   ! spatial consistency check (SSC):
    qcfgood =  0.7_wp  ,& ! factor to the main threshold to distinguish
                          ! 'good' from 'probably good' obs.
                          ! (all except 'good' obs. wiil be SSC-checked)
    qcfpbad =  1.5_wp  ,& ! factor to the main threshold to distinguish
                          ! 'probably bad' from 'bad' obs.
                          ! (all except 'bad' obs. will be used for the SSC)
    qcfspat =  0.8_wp  ,& ! reduction factor to the threshold used for
                          ! the SSC step itself
    qcftsiv =  0.8_wp  ,& ! reduction factor to the threshold used for
                          ! the IWV-SSC step itself
    qcfbias =  0.5_wp  ,& ! fraction of ABS( bias ) to be added to the
                          ! threshold used for the SSC step itself
    dtchkps =  2.0_wp  ,& ! time radius of influence for the linearly
                          ! decreasing temporal weights used in the SSC
    qcqinvf =  0.2_wp  ,& ! humidity first guess error enhancement per [K]
                          ! of observed inversion (--> inversion term)
    qcqlapf =  0.25_wp ,& ! maximum humidity first guess error enhancement
                          ! if obs has zero lapse rate ('stability term')
    qcqlapl = -0.0065_wp,&! lapse rate limit: if obs profile is more stable
                          ! then humidity first guess error is enhanced
!   qcqladz =  1.0_wp  ,& ! maximum stability factor enhancement
    qcqlapi =  0.05_wp ,& ! scaling factor of positive lapse rate for
                          ! enhancing inversion term of f.g. enhancement
    qcqfge1 =  0.10_wp ,& ! first guess humidity error at latitude > qcqlatl
    qcqfge2 =  0.15_wp ,& ! first guess humidity error at latitude < qcqlatl
    qcqlatl = 30.00_wp ,& ! latitude limit (for first guess humidity error)
!   qcqfgef (3) = (/ 0.667_wp,   1.278_wp,   1.944_wp /)
!   qcqfgef (3) = (/ 0.972_wp,   1.458_wp,   1.944_wp /)
    qcqfgef (3) = (/ 1.800_wp,   2.500_wp,   3.100_wp /)
                              ! factors for first guess humidity errors

!-------------------------------------------------------------------------------

END MODULE data_obs_qc_limits
