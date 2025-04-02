!
!+ soil parameters relevant for ensemble SMA
!
! $Id$
!
module mo_soil
!
! Description:
!   Soil parameters relevant for ensemble SMA
!   taken from GME gmtri_V2_27 source file gme_datsoil.h
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_26        2013/06/27 Andreas Rhodin
!  initial version
! V1_28        2014/02/26 Andreas Rhodin
!  provide soil_half_level_mm  (soil level thickness in mm)
! V1_31        2014-08-21 Andreas Rhodin
!  public: n_sl, n_st, soil_half_level, soil_full_level, soil_thickness
! V1_42        2015-06-08 Andreas Rhodin
!  mo_soil: function for calculating soil moisture index
! V1_43        2015-08-19 Andreas Rhodin
!  ind_wso, wso_ind: option to normalize soil water (ice) by 0/pore volume
! V1_45        2015-12-15 Andreas Rhodin
!  make 'st' public (soil level thickness (mm))
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

!=============
! Modules used
!=============
use mo_kind,        only: wp            ! working precision kind parameter
implicit none

!================
! Public entities
!================
private
!-----------
! soil types
!-----------
public :: n_st                ! number of soil types
public :: st_name             ! names of soil types
public :: ST_ICE, ST_ROCK, ST_SAND, ST_SANDLOAM, ST_LOAM, ST_CLAYLOAM, &
          ST_CLAY, ST_PEAT, ST_SEAWATER, ST_SEAICE
!-----------------------------------------
! soil properties (depending on soil type)
!-----------------------------------------
public :: cporv               ! pore volume    (fraction of volume)
public :: cfcap               ! field capacity (fraction of volume)
public :: cpwp                ! plant wilting point  (fraction of volume)
public :: cadp                ! air dryness point  (fraction of volume)
public :: crhoc               ! soil heat capacity    J/(K*m**3)
public :: cik2                ! minimum infiltration rate  kg/(s*m**2)
!--------------------------------------
! thicknesses (depending on soil depth)
!--------------------------------------
public :: n_sl                ! number of soil layers
public :: soil_half_level     ! soil half level depth (m)
public :: soil_full_level     ! soil full level depth (m)
public :: soil_thickness      ! soil level thickness  (m)
public :: st                  ! soil level thickness  (mm)
!---------
! routines
!---------
public :: ind_wso             ! soil moisture index from soil moisture
public :: wso_ind             ! soil moisture       from soil moisture index

!=================
! Module Variables
!=================

!------------------
! Define soil types
!------------------
  integer, parameter :: n_st = 10 ! number of soil types

  integer, parameter :: &! soil type
    ST_ICE      =  1   ,&! ice
    ST_ROCK     =  2   ,&! rock
    ST_SAND     =  3   ,&! sand
    ST_SANDLOAM =  4   ,&! sandy loam
    ST_LOAM     =  5   ,&! loam
    ST_CLAYLOAM =  6   ,&! clay loam
    ST_CLAY     =  7   ,&! clay
    ST_PEAT     =  8   ,&! peat
    ST_SEAWATER =  9   ,&! sea water
    ST_SEAICE   = 10     ! sea ice

!----------------------------------------------
! list with soil type names (for printout etc.)
!----------------------------------------------
  character(len=10) ,parameter :: st_name(10) = &
    (/'ice       ',&
      'rock      ',&
      'sand      ',&
      'sandy loam',&
      'loam      ',&
      'clay loam ',&
      'clay      ',&
      'peat      ',&
      'sea water ',&
      'sea ice   '/)

!---------------------------------------
! soil properties depending on soil type
!---------------------------------------
!
! Notes:
!
! all parameters should be checked against LM version, because several
! 'constants' are supplied via Namelist in LM version
!
!    Each array element refers to the corresponding property of one
!    particular soil type, i.e.
!
!    element    soil type
!     1         ice
!     2         rock
!     3         sand
!     4         sandy loam
!     5         loam
!     6         clay loam
!     7         clay
!     8         peat
!     9         sea water
!    10         sea ice
!
  real(wp) ,parameter :: & ! pore volume    (fraction of volume)
       Cporv(10) = (/    &
       1.E-10_wp,1.E-10_wp, .364_wp, .445_wp, .455_wp,&
        .475_wp, .507_wp, .863_wp,1.E-10_wp,1.E-10_wp/)

  real(wp) ,parameter :: & ! field capacity (fraction of volume)
       Cfcap(10) = (/    &
       1.E-10_wp,1.E-10_wp, .196_wp, .260_wp, .340_wp,&
        .370_wp, .463_wp, .763_wp,1.E-10_wp,1.E-10_wp/)

  real(wp) ,parameter :: & ! plant wilting point  (fraction of volume)
       Cpwp(10) =  (/    &
        0.0_wp  , 0.0_wp  , .042_wp, .100_wp, .110_wp,&
        .185_wp, .257_wp, .265_wp, 0.0_wp  , 0.0_wp  /)

  real(wp) ,parameter :: & ! air dryness point  (fraction of volume)
       Cadp(10) =  (/    &
        0.0_wp  , 0.0_wp  , .012_wp, .030_wp, .035_wp,&
       .060_wp, .065_wp, .098_wp, 0.0_wp  , 0.0_wp  /)

  real(wp) ,parameter :: & ! soil heat capacity    J/(K*m**3)
       Crhoc(10) = (/    &
        1.92E6_wp, 2.10E6_wp, 1.28E6_wp, 1.35E6_wp, 1.42E6_wp,&
        1.50E6_wp, 1.63E6_wp, 0.58E6_wp, 4.18E6_wp, 1.92E6_wp/)

  real(wp) ,parameter :: & ! minimum infiltration rate  kg/(s*m**2)
       Cik2(10)  = (/    &
        0.0_wp   , 0.0_wp   , 0.0035_wp, 0.0023_wp, 0.0010_wp,&
        0.0006_wp, 0.0001_wp, 0.0002_wp, 0.0_wp   , 0.0_wp   /)

  !---------------------------------------------------
  ! Half-levels of soil model, hard-coded, for GRIB1
  ! The surface level is bounded by two half levels 0.
  ! These values are same as in 'mo_grib', but it is
  ! the 'thickness' in 'mm' (instead of bounds)
  !---------------------------------------------------
  integer,  parameter :: n_sl = 8 ! number of soil layers

  real(wp) ,parameter :: soil_half_level(0:9) = &
         (/ 0., 0., .01, .03, .09, .27, .81, 2.43, 7.29, 21.87 /)
  real(wp) ,parameter :: soil_full_level(0:8) = &
         (soil_half_level(0:8) + soil_half_level(1:9)) / 2
  real(wp) ,parameter :: soil_thickness (0:8) =       &
         soil_half_level(1:9) - soil_half_level(0:8)

  !----------------------------------------------------
  ! some constants for soil moisture index calculations
  !----------------------------------------------------
  integer             :: ist
  real(wp) ,parameter :: cl(n_st,3) = reshape ([cadp ,cpwp ,[(0._wp,ist=1,n_st)]],shape (cl))
  real(wp) ,parameter :: cu(n_st,3) = reshape ([cporv,cfcap,cporv               ],shape (cu))
  real(wp) ,parameter :: c (n_st,3) = cu - cl
  real(wp) ,parameter :: st(n_sl)   = soil_thickness  (1:) * 1000._wp ! in mm

!==============================================================================
contains
!==============================================================================

  elemental function ind_wso (wso, soiltyp, level, mode) result (ind)
  real(wp)             :: ind      ! soil moisture index
  real(wp) ,intent(in) :: wso      ! soil moisture
  integer  ,intent(in) :: soiltyp  ! soil type
  integer  ,intent(in) :: level    ! soil level
  integer  ,intent(in) :: mode     ! 1: use adp,porv; 2: use pwp,fcap
  !---------------------------------------
  ! soil moisture index from soil moisture
  !---------------------------------------
    if (c(soiltyp,mode) > 1.1e-10_wp) then
      ind = (wso / st(level) - cl (soiltyp, mode)) / c (soiltyp, mode)
    else
      ind = 0._wp
    endif
  end function ind_wso

!------------------------------------------------------------------------------

  elemental function wso_ind (ind, soiltyp, level, mode) result (wso)
  real(wp)             :: wso      ! soil moisture
  real(wp) ,intent(in) :: ind      ! soil moisture index
  integer  ,intent(in) :: soiltyp  ! soil type
  integer  ,intent(in) :: level    ! soil level
  integer  ,intent(in) :: mode     ! 1: use adp,porv; 2: use pwp,fcap
  !---------------------------------------
  ! soil moisture from soil moisture index
  !---------------------------------------
    wso = st(level) * (cl (soiltyp, mode) + c (soiltyp, mode) * ind)
  end function wso_ind

!==============================================================================
end module mo_soil
