!
!+ provide GRIB2 table entry definitions
!
! $Id$
!
MODULE mo_t_grib_api
!
! Description:
!   provide GRIB2 table entry definitions
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_13        2011/11/01 Andreas Rhodin
!  changes for GRIB2 API
! V1_28        2014/02/26 Harald Anlauf
!  Fix GRIB2 encoding for analysis increment; clean up for ICON, IAU scheme
! V1_46        2016-02-05 Harald Anlauf
!  GRIB2: use tables 4.0, 4.7; DWD local definiton 252, encoding of mean, spread
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!------------------------------------------------------------------------------

implicit none
public

!==========================
! 0: Indicator Section (IS)
!==========================

!--------------------------------------------------------------------------
! Section 0, Code table 0: Discipline of processed data in the GRIB message
!                          number of GRIB Master Table
!--------------------------------------------------------------------------
integer, parameter ::      &!
  GRIB_0_0_METEO    =   0 ,&!
  GRIB_0_0_HYDROL   =   1 ,&!
  GRIB_0_0_LANDSURF =   2 ,&!
  GRIB_0_0_SPACE    =   3 ,&!
!                 4 -   9   ! Reserved
  GRIB_0_0_OCEAN    =  10 ,&!
!                11 - 191   ! Reserved
!               192 - 254   ! Reserved for local use
  GRIB_0_0_MISSING  = 255   ! Missing

!================================
! 1: Identification Section (IDS)
!================================

!-----------------------------------------------------------
! Section 1, Code table 0: GRIB Master Tables Version Number
!-----------------------------------------------------------
integer, parameter ::      &!
  GRIB_1_0_EXP      =   0 ,&! Experimental
  GRIB_1_0_2001     =   1 ,&! Version implemented on  7 November  2001
  GRIB_1_0_2003     =   2 ,&! Version implemented on  4 November  2003
  GRIB_1_0_2005     =   3 ,&! Version implemented on  2 November  2005
  GRIB_1_0_2007     =   4 ,&! Version implemented on  7 November  2007
  GRIB_1_0_2009     =   5 ,&! Version implemented on  4 November  2009
  GRIB_1_0_2010     =   6 ,&! Version implemented on 15 September 2010
  GRIB_1_0_PREOP    =   7 ,&! Pre-operational to be implemented next
!                  .. 254   ! Future versions
  GRIB_1_0_MISSING  = 255   ! Missing value


!----------------------------------------------------------
! Section 1, Code table 1: GRIB Local Tables Version Number
!----------------------------------------------------------
integer, parameter ::       &!
  GRIB_1_1_NOT_USED  =   0 ,&! Local tables not used
!                   .. 254   ! Number of local tables version used
  GRIB_1_1_MISSING   = 255   ! Missing


!--------------------------------------------------------
! Section 1, Code table 2: Significance of Reference Time
!--------------------------------------------------------
integer, parameter ::           &!
  GRIB_1_2_ANALYSIS      =   0 ,&! Analysis
  GRIB_1_2_FC_START      =   1 ,&! Start of forecast
  GRIB_1_2_FC_VERI       =   2 ,&! Verifying time of forecast
  GRIB_1_2_OBS           =   3 ,&! Observation time
  GRIB_1_2_MISSING       = 255   ! Missing
!                    192 - 254   ! Reserved for local use


!---------------------------------------------------
! Section 1, Code table 3: Production status of data
!---------------------------------------------------
integer, parameter ::           &!
  GRIB_1_3_OPER        =   0 ,&! Operational products
  GRIB_1_3_PRE_OPER    =   1 ,&! Operational test products
  GRIB_1_3_RESEARCH    =   2 ,&! Research products
  GRIB_1_3_RE_ANA      =   3 ,&! Re-analysis products
  GRIB_1_3_TIGGE       =   4 ,&! THORPEX Interactive Grand Global Ensemble
  GRIB_1_3_TIGGE_TEST  =   5 ,&! THORPEX Interactive Grand Global Ensemble test
!                     .. 191   ! Reserved
!                  192 - 254   ! Reserved for local use
  GRIB_1_3_MISSING     = 255   ! Missing


!--------------------------------------
! Section 1, Code table 4: Type of data
!--------------------------------------
integer, parameter ::           &!
  GRIB_1_4_ANALYSIS      =   0 ,&! Analysis products
  GRIB_1_4_FORECAST      =   1 ,&! Forecast products
  GRIB_1_4_ANA_FC        =   2 ,&! Analysis and forecast products
  GRIB_1_4_CNTR_FC       =   3 ,&! Control forecast products
  GRIB_1_4_PERT_FC       =   4 ,&! Perturbed forecast products
  GRIB_1_4_CNTR_PERT_FC  =   5 ,&! Control and perturbed forecast products
  GRIB_1_4_SAT_OBS       =   6 ,&! Processed satellite observations
  GRIB_1_4_RADAR_OBS     =   7 ,&! Processed radar observations
  GRIB_1_4_EVENT_PROB    =   8 ,&! Event Probability
!                      9 - 191 ,&! Reserved
  GRIB_1_4_MISSING       = 255 ,&! Missing
!                    192 - 254 ,&! Reserved for local use
  EDZW1_4_PERT_ANA       = 192   ! Perturbed analysis

!===========================
! 2: Local Section Use (LOC)
!===========================

!=================================
! 3: Grid definition section (GDS)
!=================================

!====================================
! 4: Product definition section (PDS)
!====================================

!------------------------------------------------------------
! Section 4, Code table 0: Product Definition Template Number
!------------------------------------------------------------
integer, parameter ::          &!
  GRIB_4_0_ANA_FC      =    0 ,&! Analysis or forecast
  GRIB_4_0_FC_ENS      =    1 ,&! Individual ensemble forecast
  GRIB_4_0_FC_ENS_ALL  =    2 ,&! Derived forecast based on all ensemble members
  GRIB_4_0_FC_ENS_REC  =    3 ,&! Derived forecast based on cluster over rectangular area
  GRIB_4_0_FC_ENS_CIR  =    4 ,&! Derived forecast based on cluster over circular area
  GRIB_4_0_FC_PROB     =    5 ,&! Probability forecasts
  GRIB_4_0_FC_PERC     =    6 ,&! Percentile forecasts
  GRIB_4_0_ANA_FC_ERR  =    7 ,&! Analysis or forecast error
  GRIB_4_0_INT         =    8 ,&! Average, accumulation, extreme, processed values in time
  GRIB_4_0_FC_PROB_INT =    9 ,&! Probability forecasts in time interval
  GRIB_4_0_FC_PERC_INT =   10 ,&! Percentile forecasts in time interval
  GRIB_4_0_FC_ENS_INT  =   11 ,&! Individual ensemble forecast in time interval
  GRIB_4_0_FC_ALL_INT  =   12 ,&! Derived forecast, all members, time interval
  GRIB_4_0_FC_REC_INT  =   13 ,&! Derived forecast, rectangular area, time interval
  GRIB_4_0_FC_CIR_INT  =   14 ,&! Derived forecast, circular area, time interval
  GRIB_4_0_AREA        =   15 ,&! Average, accumulation, extreme, processed values in area
  GRIB_4_0_RADAR       =   20 ,&! Radar product
  GRIB_4_0_SAT_DEP     =   30 ,&! Satellite product (deprecated)
  GRIB_4_0_SAT         =   31 ,&! Satellite product
  GRIB_4_0_ANA_FC_SYN  =   32 ,&! Analysis or forecast for synthetic satellite data
  GRIB_4_0_ANA_FC_CHEM =   40 ,&! Analysis or forecast for chemical constituents
  GRIB_4_0_CHEM_ENS    =   41 ,&! Individual ensemble forecast for chemical constituents
  GRIB_4_0_CHEM_INT    =   42 ,&! Average, accumulation, extreme, processed values in time
  GRIB_4_0_CHEM_ENS_INT=   43 ,&! Individual ensemble forecast, constituents, in time
  GRIB_4_0_AERO        =   44 ,&! Aerosol analysis or forecast
  GRIB_4_0_AERO_ENS    =   45 ,&! Aerosol individual ensemble forecast
  GRIB_4_0_AERO_INT    =   46 ,&! Aerosol average, accumulation, .. in time
  GRIB_4_0_AERO_ENS_INT=   47 ,&! Aerosol individual ensemble forecast in time interval
  GRIB_4_0_MULTI       =   50 ,&! multi component parameter or matrix element
  GRIB_4_0_CAT         =   51 ,&! Categorical forecasts
  GRIB_4_0_CAT_INT     =   91 ,&! Categorical forecasts in time interval
  GRIB_4_0_CCITT_IA5   =  254 ,&! CCITT IA5 character string
  GRIB_4_0_CROSS       = 1000 ,&! Cross section at a point in time
  GRIB_4_0_CROSS_INT   = 1001 ,&! Cross section processed over a range of time
  GRIB_4_0_CROSS_PROC  = 1002 ,&! Cross section processed
  GRIB_4_0_HOV         = 1100 ,&! Hovmoeller-type grid
  GRIB_4_0_HOV_PROC    = 1101 ,&! Hovmoeller-type grid with processing
  GRIB_4_0_MISSING     =65535   !

!----------------------------------------------------
! Section 4, Code table 3: Type of generating process
!----------------------------------------------------
integer, parameter ::       &!
  GRIB_4_3_ANA       =   0, &! Analysis
  GRIB_4_3_INI       =   1, &! Initialization
  GRIB_4_3_FC        =   2, &! Forecast
  GRIB_4_3_BC_FC     =   3, &! Bias corrected forecast
  GRIB_4_3_ENS_FC    =   4, &! Ensemble forecast
  GRIB_4_3_PROB_FC   =   5, &! Probability forecast
  GRIB_4_3_FC_ERR    =   6, &! Forecast error
  GRIB_4_3_ANA_ERR   =   7, &! Analysis error
  GRIB_4_3_OBS       =   8, &! Observation
  GRIB_4_3_CLIM      =   9, &! Climatological
  GRIB_4_3_PROB_W_FC =  10, &! Probability-weighted forecast
  GRIB_4_3_BC_ENS_FC =  11, &! Bias-corrected ensemble forecast
  GRIB_4_3_ANA_INC   = 201, &! Diff. analysis - first guess (DWD)
  GRIB_4_3_NUDGING   = 202, &! Nudging  (DWD)
  GRIB_4_3_NUDGCAST  = 203, &! Nudgcast (DWD)
  GRIB_4_3_DOM_FC    = 205, &! Deterministic forecast from ensemble mean (DWD)
  GRIB_4_3_MISSING   = 255   ! Missing

!---------------------------------------------------------
! Section 4, Code table 4: Indicator of unit of time range
!---------------------------------------------------------
integer, parameter ::        &!
  GRIB_4_4_MINUTE     =   0, &!
  GRIB_4_4_HOUR       =   1, &!
  GRIB_4_4_DAY        =   2, &!
  GRIB_4_4_MONTH      =   3, &!
  GRIB_4_4_YEAR       =   4, &!
  GRIB_4_4_DECADE     =   5, &!  10 years
  GRIB_4_4_NORMAL     =   6, &!  30 years
  GRIB_4_4_CENTURY    =   7, &! 100 years
  GRIB_4_4_3_HOURS    =  10, &!
  GRIB_4_4_6_HOURS    =  11, &!
  GRIB_4_4_12_HOURS   =  12, &!
  GRIB_4_4_SECOND     =  13, &!
!                   192-254, &! Reserved for local use
  GRIB_4_4_MISSING    = 254   !

!------------------------------------------
! Section 4, Code table 7: Derived forecast
!------------------------------------------
integer, parameter ::        &!
  GRIB_4_7_MEAN       =   0, &! Unweighted mean of all members
  GRIB_4_7_WMEAN      =   1, &! Weighted   mean of all members
  GRIB_4_7_STDV_C     =   2, &! Standard deviation with respect to cluster mean
  GRIB_4_7_STDV_C_N   =   3, &! Dto. normalized
  GRIB_4_7_SPREAD     =   4, &! Spread of all members
  GRIB_4_7_ANOM       =   5, &! Large anomaly index of all members
  GRIB_4_7_MEAN_C     =   6, &! Unweighted mean of the cluster members
  GRIB_4_7_QUARTILE   =   7, &! Interquartile range
  GRIB_4_7_MIN        =   8, &! Minimum of all ensemble members
  GRIB_4_7_MAX        =   9   ! Maximum of all ensemble members

!---------------------------------------------------------
! Section 4, Code table 10: Type of statistical processing
!---------------------------------------------------------
integer, parameter ::      &!
  GRIB_4_10_AVG     =   0, &! Average
  GRIB_4_10_ACCUM   =   1, &! Accumulation
  GRIB_4_10_MAX     =   2, &! Maximum
  GRIB_4_10_MIN     =   3, &! Minimum
  GRIB_4_10_DIFF    =   4, &! Difference (end of time range minus beginning)
  GRIB_4_10_RMS     =   5, &! Root mean square
  GRIB_4_10_SD      =   6, &! Standard deviation
  GRIB_4_10_COV     =   7, &! Covariance (Temporal variance)
  GRIB_4_10_DIF2    =   8, &! Difference (start of time range minus end)
  GRIB_4_10_RATIO   =   9, &! Ratio
!                 192-254, &! Reserved for local use
  GRIB_4_10_MISSING = 255   ! Missing

end module mo_t_grib_api
