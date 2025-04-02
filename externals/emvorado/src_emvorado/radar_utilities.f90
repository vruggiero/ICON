! Source module for the radar forward operator EMVORADO
!
! ---------------------------------------------------------------
! Copyright (C) 2005-2024, DWD, KIT
! Contact information: ulrich.blahak (at) dwd.de 
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

#ifdef _DACE_
#define __COSMO__
#endif


MODULE radar_utilities

!------------------------------------------------------------------------------
!
! Description:
!   This module provides some utilities which are necessary for the
!    radar forward operator EMVORADO.
!
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!
  
  USE radar_kind, ONLY : dp, sp
  
  USE radar_data, ONLY :    &
       miss_threshold, miss_value, zero_value, Z_crit_radar, dBZ_crit_radar, &
       missthr_sp, missval_sp, zeroval_sp, Z_crit_radar_sp, dBZ_crit_radar_sp

!==============================================================================

  IMPLICIT NONE

  PUBLIC

!==============================================================================

  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND( 18) !< at least 8 byte integer

  REAL(kind=dp), PARAMETER :: pi_dp     = 4.0_dp * ATAN(1.0_dp)
  REAL(kind=dp), PARAMETER :: degrad_dp = pi_dp / 180.0_dp
  REAL(kind=dp), PARAMETER :: raddeg_dp = 180.0_dp / pi_dp

  REAL(kind=sp), PARAMETER :: pi_sp     = 4.0_sp * ATAN(1.0_sp)
  REAL(kind=sp), PARAMETER :: degrad_sp = pi_sp / 180.0_sp
  REAL(kind=sp), PARAMETER :: raddeg_sp = 180.0_sp / pi_sp


  PRIVATE :: i8, &
             pi_dp, degrad_dp, raddeg_dp, &
             pi_sp, degrad_sp, raddeg_sp

#if defined(__PGI) && __PGIC__ <= 19
#  define NEED_ECL_FALLBACK
#endif

#ifdef NEED_ECL_FALLBACK
  PRIVATE :: execute_command_line
#endif

!==============================================================================

! Interface Blocks
INTERFACE geo_dist
  MODULE PROCEDURE            &
       geo_dist_dp,           &
       geo_dist_sp
END INTERFACE geo_dist

INTERFACE geo_heading
  MODULE PROCEDURE            &
       geo_heading_dp,        &
       geo_heading_sp
END INTERFACE geo_heading

INTERFACE get_range2
  MODULE PROCEDURE            &
       get_range2_dp,         &
       get_range2_sp
END INTERFACE get_range2

INTERFACE el_loc_43
  MODULE PROCEDURE            &
       el_loc_43_dp,          &
       el_loc_43_sp
END INTERFACE el_loc_43

INTERFACE el_loc_43_asin
  MODULE PROCEDURE            &
       el_loc_43_asin_dp,     &
       el_loc_43_asin_sp
END INTERFACE el_loc_43_asin

INTERFACE get_alpha3_eff_0
  MODULE PROCEDURE            &
       get_alpha3_eff_0_dp,   &
       get_alpha3_eff_0_sp
END INTERFACE get_alpha3_eff_0

INTERFACE f4_eff_horzscan
  MODULE PROCEDURE            &
       f4_eff_horzscan_dp,    &
       f4_eff_horzscan_sp
END INTERFACE f4_eff_horzscan

INTERFACE phi3_eff_horzscan
  MODULE PROCEDURE            &
       phi3_eff_horzscan_dp,    &
       phi3_eff_horzscan_sp
END INTERFACE phi3_eff_horzscan

INTERFACE theta3_eff_horzscan
  MODULE PROCEDURE            &
       theta3_eff_horzscan_dp,    &
       theta3_eff_horzscan_sp
END INTERFACE theta3_eff_horzscan

INTERFACE smth_az_horzscan
  MODULE PROCEDURE            &
       smth_az_horzscan_dp,    &
       smth_az_horzscan_sp
END INTERFACE smth_az_horzscan

INTERFACE smth_el_horzscan
  MODULE PROCEDURE            &
       smth_el_horzscan_dp,    &
       smth_el_horzscan_sp
END INTERFACE smth_el_horzscan

INTERFACE smoother
  MODULE PROCEDURE            &
       smoother_dp,           &
       smoother_sp
END INTERFACE smoother

INTERFACE refr_index_air
  MODULE PROCEDURE            &
       refr_index_air_dp,     &
       refr_index_air_sp
END INTERFACE refr_index_air

INTERFACE polar2geo
  MODULE PROCEDURE            &
       polar2geo_dp,          &
       polar2geo_sp
END INTERFACE polar2geo

INTERFACE polar2geo_old
  MODULE PROCEDURE            &
       polar2geo_old_dp,      &
       polar2geo_old_sp
END INTERFACE polar2geo_old

INTERFACE polar2geo_xy
  MODULE PROCEDURE            &
       polar2geo_xy_dp,          &
       polar2geo_xy_sp
END INTERFACE polar2geo_xy

INTERFACE phirot2phi
  MODULE PROCEDURE            &
       phirot2phi_dp,      &
       phirot2phi_sp
END INTERFACE phirot2phi

INTERFACE phi2phirot
  MODULE PROCEDURE            &
       phi2phirot_dp,      &
       phi2phirot_sp
END INTERFACE phi2phirot

INTERFACE rlarot2rla
  MODULE PROCEDURE            &
       rlarot2rla_dp,      &
       rlarot2rla_sp
END INTERFACE rlarot2rla

INTERFACE rla2rlarot
  MODULE PROCEDURE            &
       rla2rlarot_dp,      &
       rla2rlarot_sp
END INTERFACE rla2rlarot

INTERFACE init_vari
  MODULE PROCEDURE            &
       init_vari_dp_5D,       &
       init_vari_sp_5D,       &
       init_vari_int_5D,      &
       init_vari_dp_4D,       &
       init_vari_sp_4D,       &
       init_vari_int_4D,      &
       init_vari_dp_3D,       &
       init_vari_sp_3D,       &
       init_vari_int_3D,      &
       init_vari_dp_2D,       &
       init_vari_sp_2D,       &
       init_vari_int_2D,      &
       init_vari_dp_1D,       &
       init_vari_sp_1D,       &
       init_vari_int_1D
END INTERFACE

INTERFACE set_missing_and_correct0
  MODULE PROCEDURE            &
       set_missing_and_correct0_dp_3D,   &
       set_missing_and_correct0_sp_3D,   &
       set_missing_and_correct0_dp_1D,   &
       set_missing_and_correct0_sp_1D
END INTERFACE

INTERFACE dbz_to_linear
  MODULE PROCEDURE            &
       dbz_to_linear_dp_3D,   &
       dbz_to_linear_sp_3D,   &
       dbz_to_linear_dp_1D,   &
       dbz_to_linear_sp_1D,   &
       dbz_to_linear_dp_scal, &
       dbz_to_linear_sp_scal
END INTERFACE

INTERFACE linear_interpol
  MODULE PROCEDURE            &
       linear_interpol_dp,    &
       linear_interpol_sp
END INTERFACE linear_interpol

INTERFACE linear_to_dbz
  MODULE PROCEDURE            &
       linear_to_dbz_dp_3D,   &
       linear_to_dbz_sp_3D,   &
       linear_to_dbz_dp_1D,   &
       linear_to_dbz_sp_1D,   &
       linear_to_dbz_dp_scal, &
       linear_to_dbz_sp_scal
END INTERFACE

INTERFACE bubblesort
  MODULE PROCEDURE            &
       bubblesort_int,        &
       bubblesort_sp,         &
       bubblesort_dp
END INTERFACE bubblesort

CONTAINS

!==============================================================================

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the great circle distance d in m between two points
!    (lon1, lat1) and (lon2, lat2) on the sphere with radius r = r_earth + hmsl.
!
!   Input:
!     - lon1, lat1, lon2, lat2 in degrees
!     - r_earth, hmsl in m
!
!------------------------------------------------------------------------------


  ELEMENTAL FUNCTION geo_dist_dp (lon1, lat1, lon2, lat2, r_earth, hmsl) RESULT (d)

    REAL(kind=dp), INTENT(in) :: lon1, lat1, lon2, lat2, r_earth, hmsl
    REAL(kind=dp)             :: d

    REAL(kind=dp)             :: r, p1b, p2b, p1l, p2l

    p1b = lat1 * degrad_dp
    p2b = lat2 * degrad_dp
    p1l = lon1 * degrad_dp
    p2l = lon2 * degrad_dp

    r = r_earth + hmsl
    d = SIN(p2b)*SIN(p1b) + COS(p2b)*COS(p1b)*COS(p2l-p1l)
    d = r * ACOS(MIN(MAX(d, -1.0_dp), 1.0_dp))

  END FUNCTION geo_dist_dp

  ELEMENTAL FUNCTION geo_dist_sp (lon1, lat1, lon2, lat2, r_earth, hmsl) RESULT (d)

    REAL(kind=sp), INTENT(in) :: lon1, lat1, lon2, lat2, r_earth, hmsl
    REAL(kind=sp)             :: d

    REAL(kind=sp)             :: r, p1b, p2b, p1l, p2l

    p1b = lat1 * degrad_sp
    p2b = lat2 * degrad_sp
    p1l = lon1 * degrad_sp
    p2l = lon2 * degrad_sp

    r = r_earth + hmsl
    d = SIN(p2b)*SIN(p1b) + COS(p2b)*COS(p1b)*COS(p2l-p1l)
    d = r * ACOS(MIN(MAX(d, -1.0_sp), 1.0_sp))

  END FUNCTION geo_dist_sp

!==============================================================================

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the great circle direction alpha in degrees North from
!    (lon1, lat1) to (lon2, lat2) on the sphere.
!
!   Input:
!     - lon1, lat1, lon2, lat2 in degrees
!
!------------------------------------------------------------------------------


  ELEMENTAL FUNCTION geo_heading_dp (lon1, lat1, lon2, lat2) RESULT (alpha)

    REAL(kind=dp), INTENT(in) :: lon1, lat1, lon2, lat2
    REAL(kind=dp)             :: alpha

    REAL(kind=dp)             :: r, p1b, p2b, p1l, p2l, cos_zentrumswinkel, zentrumswinkel, argument

    p1b = lat1 * degrad_dp
    p2b = lat2 * degrad_dp
    p1l = lon1 * degrad_dp
    p2l = lon2 * degrad_dp

    cos_zentrumswinkel = SIN(p2b)*SIN(p1b) + COS(p2b)*COS(p1b)*COS(p2l-p1l)

    zentrumswinkel = ACOS(cos_zentrumswinkel)
    IF (ABS(zentrumswinkel) < 1e-12_dp) zentrumswinkel = 1e-12_dp

    argument = (SIN(p2b)-SIN(p1b)*cos_zentrumswinkel)/(COS(p1b)*SIN(zentrumswinkel))
    argument = MIN(MAX(argument, -1.0_dp), 1.0_dp)

    alpha = ACOS(argument)

    IF (p1l <= p2l) THEN
      alpha = alpha * raddeg_dp
    ELSE
      alpha = (2.0_dp*pi_dp - alpha) * raddeg_dp
    END IF

  END FUNCTION geo_heading_dp

  ELEMENTAL FUNCTION geo_heading_sp (lon1, lat1, lon2, lat2) RESULT (alpha)

    REAL(kind=sp), INTENT(in) :: lon1, lat1, lon2, lat2
    REAL(kind=sp)             :: alpha

    REAL(kind=sp)             :: r, p1b, p2b, p1l, p2l, cos_zentrumswinkel, zentrumswinkel, argument

    p1b = lat1 * degrad_sp
    p2b = lat2 * degrad_sp
    p1l = lon1 * degrad_sp
    p2l = lon2 * degrad_sp

    cos_zentrumswinkel = SIN(p2b)*SIN(p1b) + COS(p2b)*COS(p1b)*COS(p2l-p1l)

    zentrumswinkel = ACOS(cos_zentrumswinkel)
    IF (ABS(zentrumswinkel) < 1e-12_sp) zentrumswinkel = 1e-12_sp

    argument = (SIN(p2b)-SIN(p1b)*cos_zentrumswinkel)/(COS(p1b)*SIN(zentrumswinkel))
    argument = MIN(MAX(argument, -1.0_sp), 1.0_sp)

    alpha = ACOS(argument)

    IF (p1l <= p2l) THEN
      alpha = alpha * raddeg_sp
    ELSE
      alpha = (2.0_sp*pi_sp - alpha) * raddeg_sp
    END IF

  END FUNCTION geo_heading_sp

!==============================================================================

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the geographic longitude and latitude of a point
!    at great circle distance s (taken at MSL, not at radar height!) and
!    geographic azimut azi (degrees) relative to the point with geographic longitude
!    lon_radar and latitude lat_radar.
!
! Method:
!   Formulas from P. Snyder, Map projections used by the USGS, 1993
!
!
!   Input:
!     - s, r_earth, in meters
!     - lon_radar, lat_radar in degrees
!     - azi in degrees
!
!   Output:
!     - lon, lat in degrees
!
!------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE polar2geo_dp (lon_radar, lat_radar, r_earth, s, azi, lon, lat)

    IMPLICIT NONE

    REAL(kind=dp), INTENT(in)  :: lon_radar, lat_radar, r_earth, s, azi
    REAL(kind=dp), INTENT(out) :: lon, lat

    REAL(kind=dp)              :: x, y, sOr, denom

    IF (s > 1e-6_dp) THEN

      x = s * SIN(azi*degrad_dp)
      y = s * COS(azi*degrad_dp)

      sOr = s / r_earth

      IF (ABS(lat_radar) <= 89.99999_dp) THEN

        denom = s*COS(lat_radar*degrad_dp)*COS(sOr) - &
                y*SIN(lat_radar*degrad_dp)*SIN(sOr)
        denom = SIGN( MAX(ABS(denom),1e-12_dp), denom)
        lon = lon_radar*degrad_dp + ATAN( x*SIN(sOr) / denom )

      ELSE IF (lat_radar > 89.99999_dp) THEN
        ! North Pole:
        lon = lon_radar*degrad_dp + ATAN(-x / SIGN( MAX(ABS(y),1e-12_dp), y ) )
      ELSE
        ! South Pole:
        lon = lon_radar*degrad_dp + ATAN( x / SIGN( MAX(ABS(y),1e-12_dp), y ) )
      END IF

      lat = ASIN(  COS(sOr)*SIN(lat_radar*degrad_dp) + &
                 y*SIN(sOr)*COS(lat_radar*degrad_dp)/s)

      lon = lon * raddeg_dp
      lat = lat * raddeg_dp

    ELSE

      lon = lon_radar
      lat = lat_radar

    END IF

  END SUBROUTINE polar2geo_dp

  ELEMENTAL SUBROUTINE polar2geo_sp (lon_radar, lat_radar, r_earth, s, azi, lon, lat)

    IMPLICIT NONE

    REAL(kind=sp), INTENT(in)  :: lon_radar, lat_radar, r_earth, s, azi
    REAL(kind=sp), INTENT(out) :: lon, lat

    REAL(kind=sp)              :: x, y, sOr, denom

    IF (s > 1e-6_sp) THEN

      x = s * SIN(azi*degrad_sp)
      y = s * COS(azi*degrad_sp)

      sOr = s / r_earth

      IF (ABS(lat_radar) <= 89.999_sp) THEN

        denom = s*COS(lat_radar*degrad_sp)*COS(sOr) - &
                y*SIN(lat_radar*degrad_sp)*SIN(sOr)
        denom = SIGN( MAX(ABS(denom),1e-12_sp), denom)
        lon = lon_radar*degrad_sp + ATAN( x*SIN(sOr) / denom )

      ELSE IF (lat_radar > 89.999_sp) THEN
        ! North Pole:
        lon = lon_radar*degrad_sp + ATAN(-x / SIGN( MAX(ABS(y),1e-12_sp), y ) )
      ELSE
        ! South Pole:
        lon = lon_radar*degrad_sp + ATAN( x / SIGN( MAX(ABS(y),1e-12_sp), y ) )
      END IF

      lat = ASIN(  COS(sOr)*SIN(lat_radar*degrad_sp) + &
                 y*SIN(sOr)*COS(lat_radar*degrad_sp)/s)

      lon = lon * raddeg_sp
      lat = lat * raddeg_sp

    ELSE

      lon = lon_radar
      lat = lat_radar

    END IF

  END SUBROUTINE polar2geo_sp

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the geographic longitude and latitude of a point
!    at great circle distance s (taken at MSL, not at radar height!) and
!    geographic azimut azi (degrees) relative to the point with geographic longitude
!    lon_radar and latitude lat_radar.
!
! Method:
!   Formulas from U. Blahak, Dissertation, Appendix D
!
!   Input:
!     - s, r_earth, in meters
!     - lon_radar, lat_radar in degrees
!     - azi in degrees
!
!   Output:
!     - lon, lat in degrees
!
!------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE polar2geo_old_dp (lon_radar, lat_radar, r_earth, s, azi, lon, lat)

    IMPLICIT NONE

    REAL(kind=dp), INTENT(in)  :: lon_radar, lat_radar, r_earth, s, azi
    REAL(kind=dp), INTENT(out) :: lon, lat

    REAL(kind=dp)              :: d1, d2, d3

    d1 = s / r_earth

    ! compute geographical latitude (dependent on range, azimuth and elevation)
    ! (formula D.4 of Dissertation of Blahak)
    d2 = SIN(lat_radar*degrad_dp)*COS(d1) + COS(lat_radar*degrad_dp)*SIN(d1)*COS(azi*degrad_dp)
    ! restrict argument of asin to [-1,1]
    d2 = MIN(MAX(d2,-1.0_dp),1.0_dp)
    lat = ASIN(d2)

    ! compute geographical longitude (dependent on range, azimuth and elevation)
    ! (formula D.5 of Dissertation of Blahak)
    d3 =(COS(d1)-SIN(lat_radar*degrad_dp)*SIN(lat))/(COS(lat_radar*degrad_dp)*COS(lat))
    ! restrict argument of acos to [-1,1]
    d3 = MIN(MAX(d3,-1.0_dp),1.0_dp)
    lon = lon_radar*degrad_dp + SIGN(1.0_dp,(pi_dp - azi*degrad_dp)) * ACOS(d3)

    lon = lon * raddeg_dp
    lat = lat * raddeg_dp

  END SUBROUTINE polar2geo_old_dp

  ELEMENTAL SUBROUTINE polar2geo_old_sp (lon_radar, lat_radar, r_earth, s, azi, lon, lat)

    IMPLICIT NONE

    REAL(kind=sp), INTENT(in)  :: lon_radar, lat_radar, r_earth, s, azi
    REAL(kind=sp), INTENT(out) :: lon, lat

    REAL(kind=sp)              :: d1, d2, d3

    d1 = s / r_earth

    ! compute geographical latitude (dependent on range, azimuth and elevation)
    ! (formula D.4 of Dissertation of Blahak)
    d2 = SIN(lat_radar*degrad_sp)*COS(d1) + COS(lat_radar*degrad_sp)*SIN(d1)*COS(azi*degrad_sp)
    ! restrict argument of asin to [-1,1]
    d2 = MIN(MAX(d2,-1.0_sp),1.0_sp)
    lat = ASIN(d2)

    ! compute geographical longitude (dependent on range, azimuth and elevation)
    ! (formula D.5 of Dissertation of Blahak)
    d3 =(COS(d1)-SIN(lat_radar*degrad_sp)*SIN(lat))/(COS(lat_radar*degrad_sp)*COS(lat))
    ! restrict argument of acos to [-1,1]
    d3 = MIN(MAX(d3,-1.0_sp),1.0_sp)
    lon = lon_radar*degrad_sp + SIGN(1.0_sp,(pi_sp - azi*degrad_sp)) * ACOS(d3)

    lon = lon * raddeg_sp
    lat = lat * raddeg_sp

  END SUBROUTINE polar2geo_old_sp

!==============================================================================

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the geographic longitude and latitude of a point
!    at great circle distance s (taken at MSL, not at radar height!) and
!    geographic azimut azi relative to the point with geographic longitude
!    lon_radar and latitude lat_radar.
!   This routine is similar to polar2geo(), but instead of polar
!    coordinates s and azi it takes the kartesian coordinates x and y
!    as input, where
!
!      x = s * sin(azi)
!      y = s * cos(azi)
!
! Method:
!   Formulas from P. Snyder, Map projections used by the USGS, 1993
!
!
!   Input:
!     - x, y, r_earth in meters
!     - lon_radar, lat_radar in degrees
!
!   Output:
!     - lon, lat in degrees
!
!------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE polar2geo_xy_dp (lon_radar, lat_radar, r_earth, x, y, lon, lat)

    IMPLICIT NONE

    REAL(kind=dp), INTENT(in)  :: lon_radar, lat_radar, r_earth, x, y
    REAL(kind=dp), INTENT(out) :: lon, lat

    REAL(kind=dp)              :: s, sOr, denom

    s = SQRT(x*x + y*y)

    IF (s > 1e-6_dp) THEN

      sOr = s / r_earth

      IF (ABS(lat_radar) <= 89.99999_dp) THEN

        denom = s*COS(lat_radar*degrad_dp)*COS(sOr) - &
                y*SIN(lat_radar*degrad_dp)*SIN(sOr)
        denom = SIGN( MAX(ABS(denom),1e-12_dp), denom)
        lon = lon_radar*degrad_dp + ATAN( x*SIN(sOr) / denom )

      ELSE IF (lat_radar > 89.99999_dp) THEN
        ! North Pole: lon = lon_radar + atan (-x / y) * raddeg
        lon = lon_radar*degrad_dp + ATAN(-x / SIGN( MAX(ABS(y),1e-12_dp), y ) )
      ELSE
        ! South Pole:
        lon = lon_radar*degrad_dp + ATAN( x / SIGN( MAX(ABS(y),1e-12_dp), y ) )
      END IF

      lat = ASIN(  COS(sOr)*SIN(lat_radar*degrad_dp) + &
                 y*SIN(sOr)*COS(lat_radar*degrad_dp)/s)

      lon = lon * raddeg_dp
      lat = lat * raddeg_dp

    ELSE

      lon = lon_radar
      lat = lat_radar

    END IF

  END SUBROUTINE polar2geo_xy_dp

  ELEMENTAL SUBROUTINE polar2geo_xy_sp (lon_radar, lat_radar, r_earth, x, y, lon, lat)

    IMPLICIT NONE

    REAL(kind=sp), INTENT(in)  :: lon_radar, lat_radar, r_earth, x, y
    REAL(kind=sp), INTENT(out) :: lon, lat

    REAL(kind=sp)              :: s, sOr, denom

    s = SQRT(x*x + y*y)

    IF (s > 1e-6_sp) THEN

      sOr = s / r_earth

      IF (ABS(lat_radar) <= 89.999_sp) THEN

        denom = s*COS(lat_radar*degrad_sp)*COS(sOr) - &
                y*SIN(lat_radar*degrad_sp)*SIN(sOr)
        denom = SIGN( MAX(ABS(denom),1e-12_sp), denom)
        lon = lon_radar*degrad_sp + ATAN( x*SIN(sOr) / denom )

      ELSE IF (lat_radar > 89.999_sp) THEN
        ! North Pole: lon = lon_radar + atan (-x / y) * raddeg
        lon = lon_radar*degrad_sp + ATAN(-x / SIGN( MAX(ABS(y),1e-12_sp), y ) )
      ELSE
        ! South Pole:
        lon = lon_radar*degrad_sp + ATAN( x / SIGN( MAX(ABS(y),1e-12_sp), y ) )
      END IF

      lat = ASIN(  COS(sOr)*SIN(lat_radar*degrad_sp) + &
                 y*SIN(sOr)*COS(lat_radar*degrad_sp)/s)

      lon = lon * raddeg_sp
      lat = lat * raddeg_sp

    ELSE

      lon = lon_radar
      lat = lat_radar

    END IF

  END SUBROUTINE polar2geo_xy_sp

!==============================================================================

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the radar range along a ray assuming 4/3 earth beam
!    propagation model above the spherical earth with radius r_earth, starting at an MSL
!    height of hmsl_radar, as function of antenna elevation epsilon and
!    horizontal great circle distance from the radar station.
!
!
!   Input:
!     - r_earth, hmsl_radar in meters
!     - Elevation epsilon in degrees
!     - Great circle distance s at height hmsl_radar in meters
!
!------------------------------------------------------------------------------


  ELEMENTAL FUNCTION get_range2_sp (s, epsilon, r_earth, hmsl_radar) RESULT (range)

    REAL(kind=sp), INTENT(in) :: s, epsilon, r_earth, hmsl_radar
    REAL(kind=sp)             :: range

    REAL(kind=sp)             :: reff

    reff = 4.0_sp/3.0_sp * r_earth + hmsl_radar

    range = SIN(s/reff) * reff / COS(s/reff + epsilon*degrad_sp)

  END FUNCTION get_range2_sp

  ELEMENTAL FUNCTION get_range2_dp (s, epsilon, r_earth, hmsl_radar) RESULT (range)

    REAL(kind=dp), INTENT(in) :: s, epsilon, r_earth, hmsl_radar
    REAL(kind=dp)             :: range

    REAL(kind=dp)             :: reff

    reff = 4.0_dp/3.0_dp * r_earth + hmsl_radar

    range = SIN(s/reff) * reff / COS(s/reff + epsilon*degrad_dp)

  END FUNCTION get_range2_dp

!==============================================================================

  !-------------------------------------------------------------------------
  ! .. Local elevation angle according to the 4/3 earth model, depending on:
  !
  !    r = radar range in m
  !    s = arc distance at sea level
  !    hmsl = height of ray over MSL in m
  !    elsta = elevation of the antenna in degrees
  !    altsta = station height over MSL in m
  !    ae = earth radius in m
  !    re = 4/3 * ae
  !
  !  Formula: el_loc = atan((ae)/(ae+hmsl) * dh/ds) with
  !                        h according to Zeng et al. (2014) Eq. (14)
  !                          using r = re * sin(s/re) / cos(elsta+s/re) from the law of sines
  !
  ! ( Approxim. for r << re: eloc = elsta + s1/re = elsta + s1ore,  Notation from below )
  !
  !-------------------------------------------------------------------------

  ELEMENTAL FUNCTION el_loc_43_sp (r, elsta, altsta, r_earth) RESULT (eloc)


    IMPLICIT NONE

    !.. Input/Output parameters:

    REAL(KIND=sp), INTENT(IN) :: r, elsta, altsta, r_earth
    REAL(KIND=sp)             :: eloc

    !.. Local variables:

    REAL(KIND=sp)             :: re, h1, sinelsta, sins1ore, s1ore, &
                                 nomi, denom, coss1ore, coss1orepel, refact

    !------------------------------
    ! "Effective" 4/3 earth radius:
    re = 4.0_sp/3.0_sp * r_earth

    !------------------
    ! helper variables:

    ! calculate height of radar beam above MSL:
    h1 = SQRT((re+altsta)*(re+altsta) + r*r + 2.0_sp*(re+altsta)*r*SIN(elsta*degrad_sp)) - re

    refact = (re+altsta) / re
    sins1ore = r * COS(elsta*degrad_sp) / (re+h1)  ! s/re after Eq. (15) of Zeng et al. (2014)
    s1ore = ASIN( sins1ore )
    coss1ore = COS(s1ore)
    coss1orepel = COS(s1ore+elsta*degrad_sp)
    sinelsta = SIN(elsta*degrad_sp)

!!$    nomi = sins1ore*coss1ore + sins1ore*sins1ore*TAN(s1ore+elsta*degrad_sp) + &
!!$           refact*coss1ore*coss1orepel*sinelsta + refact*sins1ore*SIN(s1ore+elsta*degrad_sp)*sinelsta
!!$ equivalently:
    nomi = sins1ore*coss1ore + sins1ore*sins1ore*TAN(s1ore+elsta*degrad_sp) + refact*sinelsta*COS(elsta*degrad_sp)

    denom = coss1orepel * SQRT( refact*refact*coss1orepel*coss1orepel + sins1ore*sins1ore + &
            2.0_sp*refact*sins1ore*coss1orepel*sinelsta )

    eloc = raddeg_sp * ATAN( (nomi*r_earth) / (denom*(r_earth+h1)) )

  END FUNCTION el_loc_43_sp

  ELEMENTAL FUNCTION el_loc_43_dp (r, elsta, altsta, r_earth) RESULT (eloc)


    IMPLICIT NONE

    !.. Input/Output parameters:

    REAL(KIND=dp), INTENT(IN) :: r, elsta, altsta, r_earth
    REAL(KIND=dp)             :: eloc

    !.. Local variables:

    REAL(KIND=dp)             :: re, h1, sinelsta, sins1ore, s1ore, &
                                     nomi, denom, coss1ore, coss1orepel, refact

    !------------------------------
    ! "Effective" 4/3 earth radius:
    re = 4.0_dp/3.0_dp * r_earth

    !------------------
    ! helper variables:

    ! calculate height of radar beam above MSL:
    h1 = SQRT((re+altsta)*(re+altsta) + r*r + 2.0_dp*(re+altsta)*r*SIN(elsta*degrad_dp)) - re

    refact = (re+altsta) / re
    sins1ore = r * COS(elsta*degrad_dp) / (re+h1)  ! s/re after Eq. (15) of Zeng et al. (2014)
    s1ore = ASIN( sins1ore )
    coss1ore = COS(s1ore)
    coss1orepel = COS(s1ore+elsta*degrad_dp)
    sinelsta = SIN(elsta*degrad_dp)

!!$    nomi = sins1ore*coss1ore + sins1ore*sins1ore*TAN(s1ore+elsta*degrad_dp) + &
!!$           refact*coss1ore*coss1orepel*sinelsta + refact*sins1ore*SIN(s1ore+elsta*degrad_dp)*sinelsta
!!$ equivalently:
    nomi = sins1ore*coss1ore + sins1ore*sins1ore*TAN(s1ore+elsta*degrad_dp) + refact*sinelsta*COS(elsta*degrad_dp)

    denom = coss1orepel * SQRT( refact*refact*coss1orepel*coss1orepel + sins1ore*sins1ore + &
            2.0_dp*refact*sins1ore*coss1orepel*sinelsta )

    eloc = raddeg_dp * ATAN( (nomi*r_earth) / (denom*(r_earth+h1)) )

  END FUNCTION el_loc_43_dp

  !-------------------------------------------------------------------------
  ! .. Alternative formulation for the local elevation angle according
  !    to the 4/3 earth model, depending on:
  !
  !    r = radar range in m
  !    hmsl = height of ray over MSL in m
  !    elsta = elevation of the antenna in degrees
  !    altsta = station height over MSL in m
  !    ae = earth radius in m
  !
  !  Formula: el_loc = asin( dh/dr )
  !               using h(r,elsta) from Eq. (14) of Zeng et al. (2014)
  !
  ! ( Approxim. for r << re: eloc = elsta + r/re )
  !
  !-------------------------------------------------------------------------

  ELEMENTAL FUNCTION el_loc_43_asin_sp (r, elsta, altsta, r_earth) RESULT (eloc)

    IMPLICIT NONE

    !.. Input/Output parameters:

    REAL(KIND=sp), INTENT(IN) :: r, elsta, altsta, r_earth
    REAL(KIND=sp)             :: eloc

    !.. Local variables:

    REAL(KIND=sp)             :: re, nomi, denom

    !------------------------------
    ! "Effective" 4/3 earth radius:
    re = 4.0_sp/3.0_sp * r_earth

    !------------------
    ! helper variables:

    nomi = r + (re+altsta)*SIN(elsta*degrad_sp)

    denom = SQRT((re+altsta)*(re+altsta) + r*r + 2.0_sp*(re+altsta)*r*SIN(elsta*degrad_sp))

    eloc = raddeg_sp * ASIN( nomi / denom )

  END FUNCTION el_loc_43_asin_sp

  ELEMENTAL FUNCTION el_loc_43_asin_dp (r, elsta, altsta, r_earth) RESULT (eloc)

    IMPLICIT NONE

    !.. Input/Output parameters:

    REAL(KIND=dp), INTENT(IN) :: r, elsta, altsta, r_earth
    REAL(KIND=dp)             :: eloc

    !.. Local variables:

    REAL(KIND=dp)             :: re, nomi, denom

    !------------------------------
    ! "Effective" 4/3 earth radius:
    re = 4.0_dp/3.0_dp * r_earth

    !------------------
    ! helper variables:

    nomi = r + (re+altsta)*SIN(elsta*degrad_dp)

    denom = SQRT((re+altsta)*(re+altsta) + r*r + 2.0_dp*(re+altsta)*r*SIN(elsta*degrad_dp))

    eloc = raddeg_dp * ASIN( nomi / denom )

  END FUNCTION el_loc_43_asin_dp

!==============================================================================

  !==============================================================================
  !+ Module procedure for determining the half width of the
  !  effective beam weighting function
  !------------------------------------------------------------------------------

  FUNCTION get_alpha3_eff_0_dp (dalpha, phi3) RESULT(alpha3_eff_0)

    !------------------------------------------------------------------------------
    !
    ! Description: Computation of the half width parameter of the
    !              effective beam weigthing function of a scanning radar.
    !
    ! Method: Linear table interpolation; table vectors from Blahak (2008, JAOTECH).
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)   ::   dalpha, phi3
    REAL(KIND=dp)               ::   alpha3_eff_0

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    INTEGER         ::   ierr
    REAL(KIND=dp)   ::   xinter(1), finter(1)

    REAL(KIND=dp), PARAMETER :: danvec(21) = &
         (/ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, &
            1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0 /)

    REAL(KIND=dp), PARAMETER :: alvec(21)  = (/ &
         1.000, 1.015, 1.031, 1.063, 1.096,   &
         1.129, 1.180, 1.245, 1.322, 1.393,   &
         1.461, 1.554, 1.645, 1.733, 1.832,   &
         1.934, 2.433, 2.913, 3.432, 4.424,   &
         5.423 /)

    xinter(1) = dalpha / phi3

    ! .. Linear table interpolation with constant extrapolation
    !    of out-of-range values
    CALL linear_interpol(danvec, alvec, 21, &
         xinter, finter, 1, ierr)

    alpha3_eff_0 = finter(1) * phi3

  END FUNCTION get_alpha3_eff_0_dp

  FUNCTION get_alpha3_eff_0_sp (dalpha, phi3) RESULT(alpha3_eff_0)

    !------------------------------------------------------------------------------
    !
    ! Description: Computation of the half width parameter of the
    !              effective beam weigthing function of a scanning radar.
    !
    ! Method: Linear table interpolation; table vectors from Blahak (2008, JAOTECH).
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    IMPLICIT NONE

    REAL(KIND=sp), INTENT(in)   ::   dalpha, phi3
    REAL(KIND=sp)               ::   alpha3_eff_0

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    INTEGER         ::   ierr
    REAL(KIND=sp)   ::   xinter(1), finter(1)

    REAL(KIND=sp), PARAMETER :: danvec(21) = &
         (/ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, &
            1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0 /)

    REAL(KIND=sp), PARAMETER :: alvec(21)  = (/ &
         1.000, 1.015, 1.031, 1.063, 1.096,   &
         1.129, 1.180, 1.245, 1.322, 1.393,   &
         1.461, 1.554, 1.645, 1.733, 1.832,   &
         1.934, 2.433, 2.913, 3.432, 4.424,   &
         5.423 /)

    xinter(1) = dalpha / phi3

    ! .. Linear table interpolation with constant extrapolation
    !    of out-of-range values
    CALL linear_interpol(danvec, alvec, 21, &
         xinter, finter, 1, ierr)

    alpha3_eff_0 = finter(1) * phi3

  END FUNCTION get_alpha3_eff_0_sp

  !=========================================================================
  !
  ! .. Effective twoway beam function of an azimutally scanning radar
  !    after the parameterization of Blahak (2008), JAOTECH.
  !
  ! RESULT:  Value of the beam function for a specific sub-ray at
  !          azimut al and elevation el relative to the radar coordinate system.
  !
  ! INPUTS:  ista    : index of radar station
  !          az0     : azimut of beam center    [degrees]
  !          el0     : elevation of beam center [degrees]
  !          az      : azimut of the respective sub-ray [degrees]
  !          al      : elevation of the respective sub-ray [degrees]
  !
  !=========================================================================

  FUNCTION f4_eff_horzscan_dp(alpha3_eff_0,dalpha,Phi3,Theta3,az0,el0,az,el) RESULT(f4eff)

    REAL    (KIND=dp), INTENT(in)  :: alpha3_eff_0, dalpha, Phi3, Theta3, az0, el0, az, el
    REAL    (KIND=dp)              :: f4eff, phi3_eff

    phi3_eff = alpha3_eff_0 + &
         ( COS(el0*degrad_dp) - 1.0_dp ) * dalpha * &
         ( 1.0_dp - EXP(-1.5_dp*dalpha/Phi3 ) )

    f4eff = EXP(-8.0_dp*LOG(2.0_dp) * ( &
                   ( (az - az0)*COS(el*degrad_dp) / phi3_eff )**2 + &
                   ( (el - el0)                / Theta3   )**2 ) )

  END FUNCTION f4_eff_horzscan_dp

  FUNCTION f4_eff_horzscan_sp(alpha3_eff_0,dalpha,Phi3,Theta3,az0,el0,az,el) RESULT(f4eff)

    REAL    (KIND=sp), INTENT(in)  :: alpha3_eff_0, dalpha, Phi3, Theta3, az0, el0, az, el
    REAL    (KIND=sp)              :: f4eff, phi3_eff

    phi3_eff = alpha3_eff_0 + &
         ( COS(el0*degrad_sp) - 1.0_sp ) * dalpha * &
         ( 1.0_sp - EXP(-1.5_sp*dalpha/Phi3 ) )

    f4eff = EXP(-8.0_sp*LOG(2.0_sp) * ( &
                   ( (az - az0)*COS(el*degrad_sp) / phi3_eff )**2 + &
                   ( (el - el0)                / Theta3   )**2 ) )

  END FUNCTION f4_eff_horzscan_sp

  FUNCTION phi3_eff_horzscan_dp(alpha3_eff_0,dalpha,Phi3,el0) RESULT(phi3_eff)

    REAL    (KIND=dp), INTENT(in)  :: alpha3_eff_0, dalpha, Phi3, el0
    REAL    (KIND=dp)              :: phi3_eff

    phi3_eff = alpha3_eff_0 + &
         ( COS(el0*degrad_dp) - 1.0_dp ) * dalpha * &
         ( 1.0_dp - EXP(-1.5_dp*dalpha/Phi3 ) )

  END FUNCTION phi3_eff_horzscan_dp

  FUNCTION phi3_eff_horzscan_sp(alpha3_eff_0,dalpha,Phi3,el0) RESULT(phi3_eff)

    REAL    (KIND=sp), INTENT(in)  :: alpha3_eff_0, dalpha, Phi3, el0
    REAL    (KIND=sp)              :: phi3_eff

    phi3_eff = alpha3_eff_0 + &
         ( COS(el0*degrad_sp) - 1.0_sp ) * dalpha * &
         ( 1.0_sp - EXP(-1.5_sp*dalpha/Phi3 ) )

  END FUNCTION phi3_eff_horzscan_sp

  FUNCTION theta3_eff_horzscan_dp(Theta3) RESULT(theta3_eff)

    REAL    (KIND=dp), INTENT(in)  :: Theta3
    REAL    (KIND=dp)              :: theta3_eff

    theta3_eff = Theta3

  END FUNCTION theta3_eff_horzscan_dp

  FUNCTION theta3_eff_horzscan_sp(Theta3) RESULT(theta3_eff)

    REAL    (KIND=sp), INTENT(in)  :: Theta3
    REAL    (KIND=sp)              :: theta3_eff

    theta3_eff = Theta3

  END FUNCTION theta3_eff_horzscan_sp

  !=========================================================================
  !
  ! Functions for calculating the smoothing points relative to the beam center.
  ! The smoothing points are computed in a way that the azimut and elevation
  ! intervals span over the local effective beam widths times a pre-defined
  ! factor, formulated relative to the radar system.
  !
  ! The smoothing points depend on the elevation of the beam center and on
  ! the beamwidth and azimutal averaging interval of the radar.
  !
  ! See Blahak (2008), JAOTECH
  !
  ! INPUTS:  iv/ih   : vertical/horizontal index of smoothing point within
  !                    smoothing kernel, corresponding to the weights from gauleg()
  !                    (Gauss-Legendre quadrature)
  !          el0     : elevation of beam center [degrees]
  !          ista    : index of radar station
  !
  !=========================================================================

  FUNCTION smth_az_horzscan_dp(alpha3_eff_0,dalpha,Phi3,smth_interv_fact,xabscsm_h,el0) RESULT(smthpoint)

    REAL    (KIND=dp), INTENT(in)     :: alpha3_eff_0, dalpha, Phi3, smth_interv_fact, xabscsm_h, el0
    REAL    (KIND=dp)                 :: phi3_eff, smthpoint

    phi3_eff = alpha3_eff_0 + &
         ( COS(el0*degrad_dp) - 1.0_dp ) * dalpha * &
         ( 1.0_dp - EXP(-1.5_dp*dalpha/Phi3 ) )

    smthpoint = 0.5_dp * phi3_eff * smth_interv_fact * xabscsm_h

  END FUNCTION smth_az_horzscan_dp

  FUNCTION smth_az_horzscan_sp(alpha3_eff_0,dalpha,Phi3,smth_interv_fact,xabscsm_h,el0) RESULT(smthpoint)

    REAL    (KIND=sp), INTENT(in)     :: alpha3_eff_0, dalpha, Phi3, smth_interv_fact, xabscsm_h, el0
    REAL    (KIND=sp)                 :: phi3_eff, smthpoint

    phi3_eff = alpha3_eff_0 + &
         ( COS(el0*degrad_sp) - 1.0_sp ) * dalpha * &
         ( 1.0_sp - EXP(-1.5_sp*dalpha/Phi3 ) )

    smthpoint = 0.5_sp * phi3_eff * smth_interv_fact * xabscsm_h

  END FUNCTION smth_az_horzscan_sp

  FUNCTION smth_el_horzscan_dp(Theta3,smth_interv_fact,xabscsm_v) RESULT(smthpoint)

    REAL    (KIND=dp), INTENT(in)   :: Theta3, smth_interv_fact, xabscsm_v
    REAL    (KIND=dp)               :: theta3_eff, smthpoint

    theta3_eff = Theta3

    smthpoint = 0.5_dp * theta3_eff * smth_interv_fact * xabscsm_v

  END FUNCTION smth_el_horzscan_dp

  FUNCTION smth_el_horzscan_sp(Theta3,smth_interv_fact,xabscsm_v) RESULT(smthpoint)

    REAL    (KIND=sp), INTENT(in)   :: Theta3, smth_interv_fact, xabscsm_v
    REAL    (KIND=sp)               :: theta3_eff, smthpoint

    theta3_eff = Theta3

    smthpoint = 0.5_sp * theta3_eff * smth_interv_fact * xabscsm_v

  END FUNCTION smth_el_horzscan_sp

!==============================================================================

!------------------------------------------------------------------------------
!
! Generic SUBROUTINE refr_index_air (f_2d_field, ie, je, nlength, nfilt)
!
!   DP version: refr_index_air_dp
!   SP version: refr_index_air_sp
!
!------------------------------------------------------------------------------
!
! Description:
!   This function computes the refractive index of air as function of
!   temperature, total pressure and vapor pressure.
!
!------------------------------------------------------------------------------

  ELEMENTAL FUNCTION refr_index_air_dp (t, p, e) RESULT (r)

    IMPLICIT NONE

    REAL(kind=dp), INTENT(in) :: t, p, e
    REAL(kind=dp)             :: r

    r = 1.0d0 + (77.6d-2*p/t - 6.0d-2*e/t + &
                            3.75d3*e/(t*t)) * 1d-6

  END FUNCTION refr_index_air_dp

  ELEMENTAL FUNCTION refr_index_air_sp (t, p, e) RESULT (r)

    IMPLICIT NONE

    REAL(kind=sp), INTENT(in) :: t, p, e
    REAL(kind=sp)             :: r

    ! computation in DP, type cast of result to SP:
    r = 1.0d0 + (77.6d-2*DBLE(p)/DBLE(t) - 6.0d-2*DBLE(e)/DBLE(t) + &
                            3.75d3*DBLE(e)/(DBLE(t)*DBLE(t))) * 1d-6

  END FUNCTION refr_index_air_sp

!==============================================================================

!------------------------------------------------------------------------------
!
! Generic SUBROUTINE smoother (f_2d_field, ie, je, nlength, nfilt)
!
!   DP version: smoother_dp
!   SP version: smoother_sp
!
!------------------------------------------------------------------------------
!
! Description:
!   This routine smoothes an arbitrary two-dimensional field (f_2d_field) by applying
!   a 2D binomial filter with nlength symetric filter weights nfilt times.
!   The filterd field is written on f_2d_field.
!
! Method:
!   - 2D binomial filter
!   - nlength should be odd
!   - for the smoothing near the domain boundaries we assume a constant
!     extrapolation of boundary values outwards
!
!------------------------------------------------------------------------------

  SUBROUTINE smoother_dp (f_2d_field, ie, je, nlength, nfilt)

    !------------------------------------------------------------------------------
    !
    ! Parameter list:
    INTEGER       , INTENT (IN)      ::    &
         ie, je,         & ! Dimension of the field
         nlength,        & ! Filter lenght
         nfilt             ! Number of iterative filerings

    REAL (KIND=dp), INTENT (INOUT)   ::    &
         f_2d_field(:,:) ! 2-d field: unfiltered at input, filtered at output

    ! Local variables
    INTEGER                    ::    &
         i, j, m, n, nwidth, ii, jj, mm, nn, ifilt

    REAL (KIND=dp)             ::    &  ! binomial filter weights and intermediate storage
         fw_1D(nlength),              &
         tmpsum(ie,je), tmpwgt(ie,je)

    !------------------------------------------------------------------------------
    ! begin subroutine smoother

    ! 1D binomial kernel:
    fw_1D(:) = REAL( binomialcoeffs(nlength) )

    ! Loop over subsequent applications of the filter:
    DO ifilt=1, nfilt

      ! Convolution of f_2d_field with binomial kernel,
      !  separated into the two coordinate directions:
      nwidth = (nlength-1)/2
!$omp parallel workshare
      tmpsum = 0.0_dp
      tmpwgt = 0.0_dp
!$omp end parallel workshare
      DO m=-nwidth, nwidth
        mm = m + nwidth + 1
!$omp parallel do private(i,j,ii)
        DO j=1, je
          DO i=1, ie
            ii = MIN(MAX(i+m, 1), ie)
            ! Mask out grid points with missing values (miss_value):
            IF (f_2d_field(i,j) > miss_threshold .AND. f_2d_field(ii,j) > miss_threshold) THEN
              tmpsum(i,j) = tmpsum(i,j) + fw_1D(mm)*f_2d_field(ii,j)
              tmpwgt(i,j) = tmpwgt(i,j) + fw_1D(mm)
            END IF
          END DO
        END DO
!$omp end parallel do
      END DO
!$omp parallel workshare
      WHERE(tmpwgt > 0.0_dp)
        f_2d_field = tmpsum / tmpwgt
      END WHERE
!$omp end parallel workshare

!$omp parallel workshare
      tmpsum = 0.0_dp
      tmpwgt = 0.0_dp
!$omp end parallel workshare
      DO n=-nwidth, nwidth
        nn = n + nwidth + 1
!$omp parallel do private(i,j,jj)
        DO j=1, je
          jj = MIN(MAX(j+n, 1), je)
          DO i=1, ie
            ! Mask out grid points with missing values (miss_value):
            IF (f_2d_field(i,j) > miss_threshold .AND. f_2d_field(i,jj) > miss_threshold) THEN
              tmpsum(i,j) = tmpsum(i,j) + fw_1D(nn)*f_2d_field(i,jj)
              tmpwgt(i,j) = tmpwgt(i,j) + fw_1D(nn)
            END IF
          END DO
        END DO
!$omp end parallel do
      END DO
!$omp parallel workshare
      WHERE(tmpwgt > 0.0_dp)
        f_2d_field = tmpsum / tmpwgt
      END WHERE
!$omp end parallel workshare

    END DO

  END SUBROUTINE smoother_dp

  SUBROUTINE smoother_sp (f_2d_field, ie, je, nlength, nfilt)

    !------------------------------------------------------------------------------
    !
    ! Parameter list:
    INTEGER       , INTENT (IN)      ::    &
         ie, je,         & ! Dimension of the field
         nlength,        & ! Filter lenght
         nfilt             ! Number of iterative filerings

    REAL (KIND=sp), INTENT (INOUT)   ::    &
         f_2d_field(:,:) ! 2-d field: unfiltered at input, filtered at output

    ! Local variables
    INTEGER                    ::    &
         i, j, m, n, nwidth, ii, jj, mm, nn, ifilt

    REAL (KIND=sp)        ::    &  ! binomial filter weights and intermediate storage
         fw_1D(nlength), fw(nlength, nlength), &
         tmpsum(ie,je), tmpwgt(ie,je)

    !------------------------------------------------------------------------------
    ! begin subroutine smoother

    ! 1D binomial kernel:
    fw_1D(:) = REAL( binomialcoeffs(nlength) )

    ! Loop over subsequent applications of the filter:
    DO ifilt=1, nfilt

      ! Convolution of f_2d_field with binomial kernel,
      !  separated into the two coordinate directions:
      nwidth = (nlength-1)/2
!$omp parallel workshare
      tmpsum = 0.0_sp
      tmpwgt = 0.0_sp
!$omp end parallel workshare
      DO m=-nwidth, nwidth
        mm = m + nwidth + 1
!$omp parallel do private(i,j,ii)
        DO j=1, je
          DO i=1, ie
            ii = MIN(MAX(i+m, 1), ie)
            ! Mask out grid points with missing values (miss_value):
            IF (f_2d_field(i,j) > missthr_sp .AND. f_2d_field(ii,j) > missthr_sp) THEN
              tmpsum(i,j) = tmpsum(i,j) + fw_1D(mm)*f_2d_field(ii,j)
              tmpwgt(i,j) = tmpwgt(i,j) + fw_1D(mm)
            END IF
          END DO
        END DO
!$omp end parallel do
      END DO
!$omp parallel workshare
      WHERE(tmpwgt > 0.0_sp)
        f_2d_field = tmpsum / tmpwgt
      END WHERE
!$omp end parallel workshare

!$omp parallel workshare
      tmpsum = 0.0_sp
      tmpwgt = 0.0_sp
!$omp end parallel workshare
      DO n=-nwidth, nwidth
        nn = n + nwidth + 1
!$omp parallel do private(i,j,jj)
        DO j=1, je
          jj = MIN(MAX(j+n, 1), je)
          DO i=1, ie
            ! Mask out grid points with missing values (miss_value):
            IF (f_2d_field(i,j) > missthr_sp .AND. f_2d_field(i,jj) > missthr_sp) THEN
              tmpsum(i,j) = tmpsum(i,j) + fw_1D(nn)*f_2d_field(i,jj)
              tmpwgt(i,j) = tmpwgt(i,j) + fw_1D(nn)
            END IF
          END DO
        END DO
!$omp end parallel do
      END DO
!$omp parallel workshare
      WHERE(tmpwgt > 0.0_sp)
        f_2d_field = tmpsum / tmpwgt
      END WHERE
!$omp end parallel workshare

    END DO

  END SUBROUTINE smoother_sp


!------------------------------------------------------------------------------
!
! compute the n binomial coefficients of order n-1 from Pascal's Triangle:
!
!------------------------------------------------------------------------------

  FUNCTION binomialcoeffs (n) RESULT (fakt)

    INTEGER :: n, fakt(1:n)
    INTEGER :: i, j, tmp(1:n)

    ! Pascal's Triangle:
    fakt(1) = 1
    DO i=1, n-1
      tmp(:) = 1
      DO j=2,i
        tmp(j) = fakt(j-1) + fakt(j)
      END DO
      fakt(:) = tmp(:)
    END DO

  END FUNCTION binomialcoeffs


!==============================================================================

!------------------------------------------------------------------------------
!
! Routines for conversion to/from 1-dim running indices (C-style addresses) from/to
!  multi-dim. sub-indices, e.g., (i,j,k)
!
!------------------------------------------------------------------------------

  !------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE ind2sub2D( il, nm, m, n )

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Convert linear running index il to 2D field indices
    !              m, n, given the first field size
    !              nm of a field f(1:nm,1:nn).
    !
    !              It is assumed that il is defined as
    !
    !        il = m + (n-1)*nm
    !
    !              and this subroutine inverts this calculation.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)  :: il, nm
    INTEGER, INTENT(out) :: m, n

    m = MOD( (il-1)      , nm) + 1
    n =      (il-1) / nm       + 1

  END SUBROUTINE ind2sub2D

  ELEMENTAL SUBROUTINE sub2ind2D( m, n, nm, il )

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Convert indices m, n of a 2D field f(1:nm,1:nn)
    !              to a linear running index il starting at 1.
    !              nm, nn are the field dimensions.
    !
    !              il is defined as
    !
    !        il = m + (n-1)*nm
    !
    !              and the subroutine ind2sub2D() reverses the calculation.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)  :: m, n, nm
    INTEGER, INTENT(out) :: il

    il = m + (n-1)*nm

  END SUBROUTINE sub2ind2D

  ELEMENTAL SUBROUTINE sub2ind3D( m, n, o, nm, nn, il )

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Convert indices m, n, o of a 3D field f(1:nm,1:nn,1:no)
    !              to a linear running index il starting at 1.
    !              nm, nn, no are the field dimensions.
    !
    !              il is defined as
    !
    !        il = m + (n-1)*nm + (o-1)*nm*nn
    !
    !              and the subroutine ind2sub3D() reverses the calculation.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)  :: m, n, o, nm, nn
    INTEGER, INTENT(out) :: il

    il = m + (n-1)*nm + (o-1)*nm*nn

  END SUBROUTINE sub2ind3D

  ELEMENTAL SUBROUTINE ind2sub3D( il, nm, nn, m, n, o )

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Convert linear running index il to 3D field indices
    !              m, n, o, given the first 2 field sizes
    !              nm, nn of a field f(1:nm,1:nn,1:no).
    !
    !              It is assumed that il is defined as
    !
    !        il = m + (n-1)*nm + (o-1)*nm*nn
    !
    !              and this subroutine inverts this calculation.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)  :: il, nm, nn
    INTEGER, INTENT(out) :: m, n, o

    m = MOD( (il-1)              , nm) + 1
    n = MOD( (il-1) / nm         , nn) + 1
    o =      (il-1) / (nm*nn)          + 1


  END SUBROUTINE ind2sub3D

  ELEMENTAL SUBROUTINE sub2ind4D( m, n, o, p, nm, nn, no, il )

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Convert indices m, n, o, p of a
    !              4D field f(1:nm,1:nn,1:no,1:np)
    !              to a linear running index il starting at 1.
    !              nm, nn, no, np are the field dimensions.
    !
    !              il is defined as
    !
    !        il = m + (n-1)*nm + (o-1)*nm*nn + (p-1)*nm*nn*no
    !
    !              and the subroutine ind2sub4D() reverses the calculation.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)  :: m, n, o, p, nm, nn, no
    INTEGER, INTENT(out) :: il

    il = m + (n-1)*nm + (o-1)*nm*nn + (p-1)*nm*nn*no

  END SUBROUTINE sub2ind4D

  ELEMENTAL SUBROUTINE ind2sub4D( il, nm, nn, no, m, n, o, p )

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Convert linear running index il to 4D field indices
    !              m, n, o, p, given the first 3 field sizes
    !              nm, nn, no of a field f(1:nm,1:nn,1:no,1:np).
    !
    !              It is assumed that il is defined as
    !
    !        il = m + (n-1)*nm + (o-1)*nm*nn + (p-1)*nm*nn*no
    !
    !              and this subroutine inverts this calculation.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)  :: il, nm, nn, no
    INTEGER, INTENT(out) :: m, n, o, p

    m = MOD( (il-1)             , nm) + 1
    n = MOD( (il-1) / nm        , nn) + 1
    o = MOD( (il-1) / (nm*nn)   , no) + 1
    p =      (il-1) / (nm*nn*no     ) + 1


  END SUBROUTINE ind2sub4D

  ELEMENTAL SUBROUTINE sub2ind5D( m, n, o, p, q, nm, nn, no, np, il )

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Convert indices m, n, o, p, q of a
    !              5D field f(1:nm,1:nn,1:no,1:np,1:nq)
    !              to a linear running index il starting at 1.
    !              nm, nn, no, np, nq are the field dimensions.
    !
    !              il is defined as
    !
    !        il = m + (n-1)*nm + (o-1)*nm*nn + (p-1)*nm*nn*no + (q-1)*nm*nn*no*np
    !
    !              and the subroutine ind2sub5D() reverses the calculation.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)  :: m, n, o, p, q, nm, nn, no, np
    INTEGER, INTENT(out) :: il

    il = m + (n-1)*nm + (o-1)*nm*nn + (p-1)*nm*nn*no + (q-1)*nm*nn*no*np

  END SUBROUTINE sub2ind5D

  ELEMENTAL SUBROUTINE ind2sub5D( il, nm, nn, no, np, m, n, o, p, q )

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Convert linear running index il to 5D field indices
    !              m, n, o, p, q, given the first 4 field sizes
    !              nm, nn, no, np of a field f(1:nm,1:nn,1:no,1:np,1:nq).
    !
    !              It is assumed that il is defined as
    !
    !        il = m + (n-1)*nm + (o-1)*nm*nn + (p-1)*nm*nn*no + (q-1)*nm*nn*no*np
    !
    !              and this subroutine inverts this calculation.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)  :: il, nm, nn, no, np
    INTEGER, INTENT(out) :: m, n, o, p, q

    m = MOD( (il-1)              , nm) + 1
    n = MOD( (il-1) / nm         , nn) + 1
    o = MOD( (il-1) / (nm*nn)    , no) + 1
    p = MOD( (il-1) / (nm*nn*no) , np) + 1
    q =      (il-1) / (nm*nn*no*np) + 1


  END SUBROUTINE ind2sub5D

!==============================================================================

  ELEMENTAL FUNCTION round_robin(i_orig, iu, io) RESULT (i_rounded)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Compute a transformed index i_rounded from a running index i_orig
    !              (starting at 0)
    !              and two bounds iu and io (iu <= io) such that i_rounded varies
    !              periodically between iu and io.
    !
    ! This routine is useful, for example, for parallelization of computations
    ! related to a vector index i_orig over PEs in a PE ID range [iu,io].
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in) :: i_orig, iu, io
    INTEGER :: i_rounded

    i_rounded = iu + MODULO(i_orig, io-iu+1)
    
  END FUNCTION round_robin
  
!==============================================================================

  ELEMENTAL FUNCTION block_robin(i_orig, i_orig_max, iu, io) RESULT (i_rounded)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Compute a transformed index i_rounded from a running index i_orig
    !              (from 0 to i_orig_max)
    !              and two bounds iu and io (iu <= io) such that i_rounded increases
    !              blockwise from iu to io (sometimes not fully to io,
    !              if i_orig_max is not much larger than io-iu+1)
    !
    ! This routine is useful, for example, for parallelization of computations
    ! related to a vector index i_orig over PEs in a PE ID range [iu,io].
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in) :: i_orig, i_orig_max, iu, io
    INTEGER :: i_rounded

    INTEGER :: blocksize

    blocksize = (i_orig_max+1) / (io-iu+1) + 1
    i_rounded = iu + MODULO(i_orig / blocksize, io-iu+1)
    
  END FUNCTION block_robin
  
!==============================================================================

!==============================================================================
!+ Function for rotation of geographical coordinates
!------------------------------------------------------------------------------

  FUNCTION  phirot2phi_dp ( phirot, rlarot, polphi, pollam, polgam )

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This function converts phi from one rotated system to phi in another
    !   system. If the optional argument polgam is present, the other system
    !   can also be a rotated one, where polgam is the angle between the two
    !   north poles.
    !   If polgam is not present, the other system is the real geographical
    !   system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !
    ! Declarations:
    !
    !------------------------------------------------------------------------------

    ! Parameter list:
    REAL (KIND=dp),     INTENT (IN)      ::        &
         polphi,   & ! latitude of the rotated north pole
         pollam,   & ! longitude of the rotated north pole
         phirot,   & ! latitude in the rotated system
         rlarot      ! longitude in the rotated system

    REAL (KIND=dp),     INTENT (IN)      ::        &
         polgam      ! angle between the north poles of the systems

    REAL (KIND=dp)                       ::        &
         phirot2phi_dp  ! latitude in the geographical system

    ! Local variables
    REAL (KIND=dp)                       ::        &
         zsinpol, zcospol, zphis, zrlas, zarg, zgam

    REAL (KIND=dp),     PARAMETER        ::        &
         zrpi18 = 57.2957795_dp,                      &
         zpir18 = 0.0174532925_dp

    !------------------------------------------------------------------------------

    ! Begin function phirot2phi

    zsinpol     = SIN (zpir18 * polphi)
    zcospol     = COS (zpir18 * polphi)

    zphis       = zpir18 * phirot
    IF (rlarot > 180.0_dp) THEN
      zrlas = rlarot - 360.0_dp
    ELSE
      zrlas = rlarot
    ENDIF
    zrlas       = zpir18 * zrlas

    IF (polgam /= 0.0_dp) THEN
      zgam  = zpir18 * polgam
      zarg  = zsinpol*SIN (zphis) +                                           &
           zcospol*COS(zphis) * ( COS(zrlas)*COS(zgam) - SIN(zgam)*SIN(zrlas) )
    ELSE
      zarg  = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
    ENDIF

    phirot2phi_dp  = zrpi18 * ASIN (zarg)

  END FUNCTION phirot2phi_dp

  FUNCTION  phirot2phi_sp ( phirot, rlarot, polphi, pollam, polgam )

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This function converts phi from one rotated system to phi in another
    !   system. If the optional argument polgam is present, the other system
    !   can also be a rotated one, where polgam is the angle between the two
    !   north poles.
    !   If polgam is not present, the other system is the real geographical
    !   system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !
    ! Declarations:
    !
    !------------------------------------------------------------------------------

    ! Parameter list:
    REAL (KIND=sp),     INTENT (IN)      ::        &
         polphi,   & ! latitude of the rotated north pole
         pollam,   & ! longitude of the rotated north pole
         phirot,   & ! latitude in the rotated system
         rlarot      ! longitude in the rotated system

    REAL (KIND=sp),     INTENT (IN)      ::        &
         polgam      ! angle between the north poles of the systems

    REAL (KIND=sp)                       ::        &
         phirot2phi_sp  ! latitude in the geographical system

    !------------------------------------------------------------------------------

    phirot2phi_sp  = REAL( phirot2phi_dp( &
         REAL(phirot, kind=dp), &
         REAL(rlarot, kind=dp), &
         REAL(polphi, kind=dp), &
         REAL(pollam, kind=dp), &
         REAL(polgam, kind=dp)  &
         ) , kind=sp)

  END FUNCTION phirot2phi_sp

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

  FUNCTION  phi2phirot_dp ( phi, rla, polphi, pollam )

    !------------------------------------------------------------------------------
    ! Description:
    !   This routine converts phi from the real geographical system to phi
    !   in the rotated system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------
    ! Parameter list:
    REAL (KIND=dp),     INTENT (IN)      ::        &
         polphi,  & ! latitude of the rotated north pole
         pollam,  & ! longitude of the rotated north pole
         phi,     & ! latitude in the geographical system
         rla        ! longitude in the geographical system

    REAL (KIND=dp)                       ::        &
         phi2phirot_dp ! longitude in the rotated system

    ! Local variables
    REAL (KIND=dp)                           ::    &
         zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

    REAL (KIND=dp),     PARAMETER            ::    &
         zrpi18 = 57.2957795_dp,                      & !
         zpir18 = 0.0174532925_dp

    !------------------------------------------------------------------------------

    ! Begin function phi2phirot

    zsinpol  = SIN (zpir18 * polphi)
    zcospol  = COS (zpir18 * polphi)
    zlampol  =      zpir18 * pollam
    zphi     =      zpir18 * phi
    IF (rla > 180.0_dp) THEN
      zrla1  = rla - 360.0_dp
    ELSE
      zrla1  = rla
    ENDIF
    zrla     = zpir18 * zrla1

    zarg1    = SIN (zphi) * zsinpol
    zarg2    = COS (zphi) * zcospol * COS (zrla - zlampol)

    phi2phirot_dp = zrpi18 * ASIN (zarg1 + zarg2)

  END FUNCTION phi2phirot_dp

  FUNCTION  phi2phirot_sp ( phi, rla, polphi, pollam )

    !------------------------------------------------------------------------------
    ! Description:
    !   This routine converts phi from the real geographical system to phi
    !   in the rotated system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------
    ! Parameter list:
    REAL (KIND=sp),     INTENT (IN)      ::        &
         polphi,  & ! latitude of the rotated north pole
         pollam,  & ! longitude of the rotated north pole
         phi,     & ! latitude in the geographical system
         rla        ! longitude in the geographical system

    REAL (KIND=sp)                       ::        &
         phi2phirot_sp ! longitude in the rotated system

    !------------------------------------------------------------------------------

    phi2phirot_sp = REAL( phi2phirot_dp ( &
         REAL(phi   , kind=dp), &
         REAL(rla   , kind=dp), &
         REAL(polphi, kind=dp), &
         REAL(pollam, kind=dp)  &
         ), kind=sp)

  END FUNCTION phi2phirot_sp

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

  FUNCTION  rlarot2rla_dp (phirot, rlarot, polphi, pollam, polgam)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This function converts lambda from one rotated system to lambda in another
    !   system. If the optional argument polgam is present, the other system
    !   can also be a rotated one, where polgam is the angle between the two
    !   north poles.
    !   If polgam is not present, the other system is the real geographical
    !   system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    ! Modules used:    NONE
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !
    ! Declarations:
    !
    !------------------------------------------------------------------------------

    ! Parameter list:
    REAL (KIND=dp),     INTENT (IN)      ::        &
         polphi,   & ! latitude of the rotated north pole
         pollam,   & ! longitude of the rotated north pole
         phirot,   & ! latitude in the rotated system
         rlarot      ! longitude in the rotated system

    REAL (KIND=dp),     INTENT (IN)      ::        &
         polgam      ! angle between the north poles of the systems

    REAL (KIND=dp)                       ::        &
         rlarot2rla_dp  ! longitude in the geographical system

    ! Local variables
    REAL (KIND=dp)                       ::        &
         zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam

    REAL (KIND=dp),     PARAMETER        ::        &
         zrpi18 = 57.2957795_dp,                      & !
         zpir18 = 0.0174532925_dp

    !------------------------------------------------------------------------------

    ! Begin function rlarot2rla

    zsinpol = SIN (zpir18 * polphi)
    zcospol = COS (zpir18 * polphi)

    zlampol = zpir18 * pollam
    zphis   = zpir18 * phirot
    IF (rlarot > 180.0_dp) THEN
      zrlas = rlarot - 360.0_dp
    ELSE
      zrlas = rlarot
    ENDIF
    zrlas   = zpir18 * zrlas

    IF (polgam /= 0.0_dp) THEN
      zgam    = zpir18 * polgam
      zarg1   = SIN (zlampol) *                                                &
           (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
           + zcospol * SIN(zphis))                                               &
           - COS (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))

      zarg2   = COS (zlampol) *                                                &
           (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
           + zcospol * SIN(zphis))                                               &
           + SIN (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
    ELSE
      zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
           zcospol *              SIN(zphis)) -    &
           COS (zlampol) *             SIN(zrlas) * COS(zphis)
      zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
           zcospol *              SIN(zphis)) +   &
           SIN (zlampol) *             SIN(zrlas) * COS(zphis)
    ENDIF

    IF (zarg2 == 0.0_dp) zarg2 = 1.0E-20_dp

    rlarot2rla_dp = zrpi18 * ATAN2(zarg1,zarg2)

  END FUNCTION rlarot2rla_dp

  FUNCTION  rlarot2rla_sp (phirot, rlarot, polphi, pollam, polgam)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This function converts lambda from one rotated system to lambda in another
    !   system. If the optional argument polgam is present, the other system
    !   can also be a rotated one, where polgam is the angle between the two
    !   north poles.
    !   If polgam is not present, the other system is the real geographical
    !   system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    ! Modules used:    NONE
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !
    ! Declarations:
    !
    !------------------------------------------------------------------------------

    ! Parameter list:
    REAL (KIND=sp),     INTENT (IN)      ::        &
         polphi,   & ! latitude of the rotated north pole
         pollam,   & ! longitude of the rotated north pole
         phirot,   & ! latitude in the rotated system
         rlarot      ! longitude in the rotated system

    REAL (KIND=sp),     INTENT (IN)      ::        &
         polgam      ! angle between the north poles of the systems

    REAL (KIND=sp)                       ::        &
         rlarot2rla_sp  ! longitude in the geographical system

    !------------------------------------------------------------------------------

    rlarot2rla_sp = REAL( rlarot2rla_dp ( &
         REAL( phirot, kind=dp), &
         REAL( rlarot, kind=dp), &
         REAL( polphi, kind=dp), &
         REAL( pollam, kind=dp), &
         REAL( polgam, kind=dp)  &
         ) , kind=sp )

  END FUNCTION rlarot2rla_sp

!==============================================================================

!------------------------------------------------------------------------------

  FUNCTION  rla2rlarot_dp ( phi, rla, polphi, pollam, polgam )

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine converts lambda from the real geographical system to lambda
    !   in the rotated system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------
    !
    ! Parameter list:
    REAL (KIND=dp),     INTENT (IN)      ::        &
         polphi,  & ! latitude of the rotated north pole
         pollam,  & ! longitude of the rotated north pole
         phi,     & ! latitude in geographical system
         rla        ! longitude in geographical system

    REAL (KIND=dp),     INTENT (IN)      ::        &
         polgam      ! angle between the north poles of the systems

    REAL (KIND=dp)                       ::        &
         rla2rlarot_dp ! longitude in the the rotated system

    ! Local variables
    REAL (KIND=dp)                           ::    &
         zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

    REAL (KIND=dp),     PARAMETER            ::    &
         zrpi18 = 57.2957795_dp,                      & !
         zpir18 = 0.0174532925_dp

    !------------------------------------------------------------------------------

    ! Begin function rla2rlarot

    zsinpol  = SIN (zpir18 * polphi)
    zcospol  = COS (zpir18 * polphi)
    zlampol  =      zpir18 * pollam
    zphi     =      zpir18 * phi
    IF (rla > 180.0_dp) THEN
      zrla1  = rla - 360.0_dp
    ELSE
      zrla1  = rla
    ENDIF
    zrla     = zpir18 * zrla1

    zarg1    = - SIN (zrla-zlampol) * COS(zphi)
    zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)

    IF (zarg2 == 0.0_dp) zarg2 = 1.0E-20_dp

    rla2rlarot_dp = zrpi18 * ATAN2 (zarg1,zarg2)

    IF (polgam /= 0.0_dp) THEN
      rla2rlarot_dp = polgam + rla2rlarot_dp
      IF (rla2rlarot_dp > 180._dp) rla2rlarot_dp = rla2rlarot_dp-360._dp
    ENDIF

  END FUNCTION rla2rlarot_dp

  FUNCTION  rla2rlarot_sp ( phi, rla, polphi, pollam, polgam )

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine converts lambda from the real geographical system to lambda
    !   in the rotated system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------
    !
    ! Parameter list:
    REAL (KIND=sp),     INTENT (IN)      ::        &
         polphi,  & ! latitude of the rotated north pole
         pollam,  & ! longitude of the rotated north pole
         phi,     & ! latitude in geographical system
         rla        ! longitude in geographical system

    REAL (KIND=sp),     INTENT (IN)      ::        &
         polgam      ! angle between the north poles of the systems

    REAL (KIND=sp)                       ::        &
         rla2rlarot_sp ! longitude in the the rotated system

    !------------------------------------------------------------------------------

    rla2rlarot_sp = REAL( rla2rlarot_dp ( &
         REAL( phi   , kind=dp), &
         REAL( rla   , kind=dp), &
         REAL( polphi, kind=dp), &
         REAL( pollam, kind=dp), &
         REAL( polgam, kind=dp)  &
         ) , kind = sp)

  END FUNCTION rla2rlarot_sp

!==============================================================================

  SUBROUTINE get_free_funit ( funit )

    !------------------------------------------------------------------------------
    !
    ! Description: Searches for free file units for a following OPEN statement
    !
    !------------------------------------------------------------------------------
    !

    INTEGER, INTENT(out) :: funit

    INTEGER, PARAMETER   :: LUN_MIN=2001, LUN_MAX=5000
    LOGICAL              :: opened
    INTEGER              :: lun

    funit = -1
    DO lun = LUN_MIN, LUN_MAX
      INQUIRE(unit=lun, opened=opened)
      IF (.NOT. opened) THEN
        funit=lun
        EXIT
      END IF
    END DO

    IF (funit == -1) THEN
      WRITE (*,*) 'ERROR get_free_funit(): no free file unit available!'
    END IF

  END SUBROUTINE get_free_funit

!==============================================================================

!==============================================================================

!**************************************************************************************
!
! Routine for linear interpolation (1D) of the M ("anz") values "fwerte" defined at
! positions "xwerte" to the N ("anzinter") positions in "xinter". "finter" returns the
! corresponding interpolated values, and ierr contains an error flag, which
! is 0 if all "anzinter" points could be interpolated.
! For points in "xinter" which are outside the range of "xwerte",
! a constant extrapolation is applied, i.e., corresponding values
! of "finter" are equal to "fwerte(1)" or "fwerte(anz)", resp.
!
! Values in "xwerte" have to be strictly monotonically increasing, whereas
! no such constraint applies for the values in "xinter".
!
!**************************************************************************************

!===============================================================
! Routine for scalar machines:
!===============================================================

  SUBROUTINE linear_interpol_dp(xwerte, fwerte, anz, xinter, finter, anzinter, ierr, missingval)

    IMPLICIT NONE

    INTEGER,       INTENT(in)  :: anz, anzinter
    REAL(KIND=dp), INTENT(in)  :: xwerte(anz), fwerte(anz), xinter(anzinter)
    REAL(KIND=dp), INTENT(in), OPTIONAL :: missingval
    INTEGER,       INTENT(out) :: ierr
    REAL(KIND=dp), INTENT(out) :: finter(anzinter)

    INTEGER                     :: i, j, cnt
    REAL(KIND=dp)               :: xtmp, xmax, xmin
    CHARACTER(LEN=*), PARAMETER :: rname = 'linear_interpol'


    ierr = 0

    IF ( ANY(xwerte(2:anz) < xwerte(1:anz-1) ) ) THEN
      ierr = 2
      WRITE (*,*) '%%%% '//TRIM(rname)//': Error --- Xwerte must be monotonically increasing!'
      RETURN
    END IF

    xmin = MINVAL(xwerte)
    xmax = MAXVAL(xwerte)
    cnt = 0

    IF (.NOT.PRESENT(missingval)) THEN

      DO i=1, anzinter

        xtmp = xinter(i)

        IF (xtmp > xmax) THEN
          cnt = cnt + 1
          finter(i) = fwerte(anz)
        ELSEIF (xtmp <= xmin) THEN
          cnt = cnt + 1
          finter(i) = fwerte(1)
        ELSE
          intervallsuchen: DO j=1, anz-1
            IF ( xtmp > xwerte(j) .AND. xtmp <= xwerte(j+1) ) THEN
              cnt = cnt + 1
              finter(i) = fwerte(j) + (fwerte(j+1)-fwerte(j)) / &
                   (xwerte(j+1)-xwerte(j)) * (xtmp-xwerte(j))
              EXIT intervallsuchen
            END IF
          END DO intervallsuchen
        END IF

      END DO

    ELSE

      DO i=1, anzinter

        xtmp = xinter(i)

        IF (xtmp > xmax .OR. xtmp <= xmin) THEN
          cnt = cnt + 1
          finter(i) = missingval
        ELSE
          intervallsuchen2: DO j=1, anz-1
            IF ( xtmp > xwerte(j) .AND. xtmp <= xwerte(j+1) ) THEN
              cnt = cnt + 1
              finter(i) = fwerte(j) + (fwerte(j+1)-fwerte(j)) / &
                   (xwerte(j+1)-xwerte(j)) * (xtmp-xwerte(j))
              EXIT intervallsuchen2
            END IF
          END DO intervallsuchen2
        END IF

      END DO

    END IF

    IF (cnt < anzinter) ierr = 1

  END SUBROUTINE linear_interpol_dp

  SUBROUTINE linear_interpol_sp(xwerte, fwerte, anz, xinter, finter, anzinter, ierr, missingval)

    IMPLICIT NONE

    INTEGER,       INTENT(in)  :: anz, anzinter
    REAL(KIND=sp), INTENT(in)  :: xwerte(anz), fwerte(anz), xinter(anzinter)
    REAL(KIND=sp), INTENT(in), OPTIONAL :: missingval
    INTEGER,       INTENT(out) :: ierr
    REAL(KIND=sp), INTENT(out) :: finter(anzinter)

    INTEGER                     :: i, j, cnt
    REAL(KIND=sp)               :: xtmp, xmax, xmin
    CHARACTER(LEN=*), PARAMETER :: rname = 'linear_interpol'


    ierr = 0

    IF ( ANY(xwerte(2:anz) < xwerte(1:anz-1) ) ) THEN
      ierr = 2
      WRITE (*,*) '%%%% '//TRIM(rname)//': Error --- Xwerte must be monotonically increasing!'
      RETURN
    END IF

    xmin = MINVAL(xwerte)
    xmax = MAXVAL(xwerte)
    cnt = 0

    IF (.NOT.PRESENT(missingval)) THEN

      DO i=1, anzinter

        xtmp = xinter(i)

        IF (xtmp > xmax) THEN
          cnt = cnt + 1
          finter(i) = fwerte(anz)
        ELSEIF (xtmp <= xmin) THEN
          cnt = cnt + 1
          finter(i) = fwerte(1)
        ELSE
          intervallsuchen: DO j=1, anz-1
            IF ( xtmp > xwerte(j) .AND. xtmp <= xwerte(j+1) ) THEN
              cnt = cnt + 1
              finter(i) = fwerte(j) + (fwerte(j+1)-fwerte(j)) / &
                   (xwerte(j+1)-xwerte(j)) * (xtmp-xwerte(j))
              EXIT intervallsuchen
            END IF
          END DO intervallsuchen
        END IF

      END DO

    ELSE

      DO i=1, anzinter

        xtmp = xinter(i)

        IF (xtmp > xmax .OR. xtmp <= xmin) THEN
          cnt = cnt + 1
          finter(i) = missingval
        ELSE
          intervallsuchen2: DO j=1, anz-1
            IF ( xtmp > xwerte(j) .AND. xtmp <= xwerte(j+1) ) THEN
              cnt = cnt + 1
              finter(i) = fwerte(j) + (fwerte(j+1)-fwerte(j)) / &
                   (xwerte(j+1)-xwerte(j)) * (xtmp-xwerte(j))
              EXIT intervallsuchen2
            END IF
          END DO intervallsuchen2
        END IF

      END DO

    END IF

    IF (cnt < anzinter) ierr = 1

  END SUBROUTINE linear_interpol_sp

!===============================================================================

!-------------------------------------------------------------------------------

  SUBROUTINE diff_seconds ( cdate1         &
       , cdate2         &
       , isecdif, irm )

    !-------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine computes the time difference of two dates/times in seconds,
    !   given the dates/times as strings 'cdate1' and 'cdate2'
    !   in the character format 'YYYYMMDDhhmmss'.
    !   Computed is the time in seconds from 'cdate1' to 'cdate2' taking
    !   into account the correct sign of the difference.
    !   So far this is only implemented for the Gregorian calendar (itype_calendar = 0).
    !
    ! Note:
    !   If there is no century in the date (the year is a 2-digit number),
    !   then the date/time format on input should be '00YYMMDDhhmmss'!
    !
    ! Method:
    !   First checks validity of input parameters (for invalid ranges, see defintion
    !   of error status variable 'irm' below). Then straightforward computation of
    !   time difference taking into account leap years.
    !
    ! Original routine diff_minutes:   C. Schraff, 14.08.2009
    ! Modified and renamed :           U. Blahak,  10.01.2012
    !
    !-------------------------------------------------------------------------------
    ! Declarations:
    !-------------------------------------------------------------------------------
    !
    ! Subroutine arguments:
    ! ---------------------

    CHARACTER(len=14), INTENT (IN)    ::       &
         cdate1, cdate2        ! Input dates in the format 'YYYYMMDDhhmmss'

    INTEGER, INTENT (OUT)   ::       &
         isecdif     ,& ! difference of date/time 2 minus date/time 1 [seconds]
         irm(2,6)       ! error status = 0/1 : no error / error
    !       irm(1,1) = 1 : invalid year 1  ( < 0 )
    !       irm(1,2) = 1 : invalid month 1  ( < 1 .or. > 12 )
    !       irm(1,3) = 1 : invalid day 1    ( < 1 .or. > days of mm )
    !       irm(1,4) = 1 : invalid hour 1   ( < 0 .or. > 23 )
    !       irm(1,5) = 1 : invalid minute 1 ( < 0 .or. > 59 )
    !       irm(1,6) = 1 : invalid second 1 ( < 0 .or. > 59 )
    !       irm(2,1) = 1 : invalid year 2  ( < 0 )
    !       irm(2,2) = 1 : invalid month 2  ( < 1 .or. > 12 )
    !       irm(2,3) = 1 : invalid day 2    ( < 1 .or. > days of mm )
    !       irm(2,4) = 1 : invalid hour 2   ( < 0 .or. > 23 )
    !       irm(2,5) = 1 : invalid minute 2 ( < 0 .or. > 59 )
    !       irm(2,6) = 1 : invalid second 2 ( < 0 .or. > 59 )

    ! Local parameters:
    ! ----------------

    INTEGER  ::       &
         iyr_sta  , iyr_end  ,& ! year   of start date / end date  [yy] or [yyyy]
         imm_sta  , imm_end  ,& ! month  of start date / end date  [mm]
         idd_sta  , idd_end  ,& ! day    of start date / end date  [dd]
         ihh_sta  , ihh_end  ,& ! hour   of start time / end time  [hh]
         imin_sta , imin_end ,& ! minute of start time / end time  [min]
         isec_sta , isec_end    ! second of start time / end time  [sec]

    INTEGER, PARAMETER  :: &
         mdd      (12) = (/31, 28, 31, 30,  31,  30,  31,  31,  30,  31,  30,  31/),&
         mdd_offs (13) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334   &
         , 365/)

    ! Local variables:
    ! ---------------

    INTEGER          ::       &
         iyyyy_sta, iyyyy_end,& ! start / end year including century [yyyy]
         iddmax1  , iddmax2  ,& ! number of days           in start / end month
         iddjul1  , iddjul2  ,& ! Julian day (within year) of start / end date
         idddif              ,& ! difference in days between  start / end date
         isign               ,& ! = -1 if end year < start year , = 1 otherwise
         iysta    , iyend    ,& ! range of loop over years
         iys                 ,& ! loop index (over years)
         ileap               ,& ! statement function for leap year
         iy                     ! dummy variable for statement function
    !
    !------------ End of header ----------------------------------------------------

    ! Statement function for leap year:  ileap(iy) = 0:  no leap year,
    ! --------------------------------   ileap(iy) = 1:  leap year
    ileap (iy) = IABS( MOD(iy,  4) -  4) /  4 & ! every     4 years is a leapyear
         -IABS( MOD(iy,100) -100) /100 & ! but every 100 ys. is no leapyear
         +IABS( MOD(iy,400) -400) /400   ! but every 400 ys. is a leapyear

    !-------------------------------------------------------------------------------
    ! Begin Subroutine diff_seconds
    !-------------------------------------------------------------------------------

    irm     = 0
    isecdif = 0

    ! get integer date/time parts from input strings:
    ! -----------------------------------------------

    READ(cdate1, '(I4,5I2)')  iyr_sta, imm_sta, idd_sta, ihh_sta, imin_sta, isec_sta
    READ(cdate2, '(I4,5I2)')  iyr_end, imm_end, idd_end, ihh_end, imin_end, isec_end

    ! check input values
    ! ------------------

    ! check year and month values
    ! as well as hours, minutes and seconds values
    IF  ( iyr_sta  < 0 )                     irm(1,1) = 1
    IF  ( imm_sta  < 1 .OR. imm_sta  > 12 )  irm(1,2) = 1
    IF  ( idd_sta  < 1 .OR. idd_sta  > 31 )  irm(1,3) = 1
    IF  ( ihh_sta  < 0 .OR. ihh_sta  > 23 )  irm(1,4) = 1
    IF  ( imin_sta < 0 .OR. imin_sta > 59 )  irm(1,5) = 1
    IF  ( isec_sta < 0 .OR. isec_sta > 59 )  irm(1,6) = 1
    IF  ( iyr_end  < 0 )                     irm(2,1) = 1
    IF  ( imm_end  < 1 .OR. imm_end  > 12 )  irm(2,2) = 1
    IF  ( idd_end  < 1 .OR. idd_end  > 31 )  irm(2,3) = 1
    IF  ( ihh_end  < 0 .OR. ihh_end  > 23 )  irm(2,4) = 1
    IF  ( imin_end < 0 .OR. imin_end > 59 )  irm(2,5) = 1
    IF  ( isec_end < 0 .OR. isec_end > 59 )  irm(2,6) = 1

    ! if century is missing in years, then add it (to end up in range 1960 - 2059)
    IF ( MAXVAL(irm) == 0 ) THEN
      iyyyy_end  =  iyr_end
      iyyyy_sta  =  iyr_sta
      IF (iyyyy_end <  60)  iyyyy_end = 2000 + iyyyy_end
      IF (iyyyy_sta <  60)  iyyyy_sta = 2000 + iyyyy_sta
      IF (iyyyy_end < 100)  iyyyy_end = 1900 + iyyyy_end
      IF (iyyyy_sta < 100)  iyyyy_sta = 1900 + iyyyy_sta

      ! check days again, but closer:
      iddmax1  =  mdd(imm_sta)
      iddmax2  =  mdd(imm_end)
      IF (imm_sta == 2)  iddmax1  =  iddmax1 + ileap(iyyyy_sta)
      IF (imm_end == 2)  iddmax2  =  iddmax2 + ileap(iyyyy_end)
      IF ( idd_sta < 1 .OR. idd_sta > iddmax1 )  irm(1,3) = 1
      IF ( idd_end < 1 .OR. idd_end > iddmax2 )  irm(2,3) = 1
    ENDIF


    IF ( MAXVAL(irm) == 0)  THEN

      ! compute time difference in seconds
      ! ----------------------------------

      ! get Julian days (i.e. days since the beginning of the years)
      iddjul1  =  mdd_offs(imm_sta) + idd_sta
      iddjul2  =  mdd_offs(imm_end) + idd_end
      IF (imm_sta > 2)  iddjul1  =  iddjul1 + ileap(iyyyy_sta)
      IF (imm_end > 2)  iddjul2  =  iddjul2 + ileap(iyyyy_end)

      ! difference of days irrespective of years
      idddif  =  iddjul2 - iddjul1

      ! difference of days in case of different years
      IF (iyyyy_end /= iyyyy_sta) THEN
        iysta = MIN( iyyyy_sta , iyyyy_end )
        iyend = MAX( iyyyy_sta , iyyyy_end ) - 1
        isign  = 1
        IF (iyyyy_end < iyyyy_sta) THEN
          idddif = - idddif
          isign  = - 1
        ENDIF

        ! loop over years: add days from each year
        DO iys = iysta , iyend
          idddif = idddif + 365 + ileap(iys)
        ENDDO
        idddif = idddif *isign
      ENDIF

      ! difference in seconds
      isecdif  =  60 * ( 60 * ( 24*idddif + (ihh_end  - ihh_sta ) )    &
           + (imin_end - imin_sta) )  &
           + (isec_end - isec_sta)

    ELSE

      ! give detailed error message(s):

      IF ( irm(1,1) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong year    in cdate1 = '//cdate1
      END IF
      IF ( irm(1,2) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong month   in cdate1 = '//cdate1
      END IF
      IF ( irm(1,3) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong day     in cdate1 = '//cdate1
      END IF
      IF ( irm(1,4) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong hour    in cdate1 = '//cdate1
      END IF
      IF ( irm(1,5) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong minutes in cdate1 = '//cdate1
      END IF
      IF ( irm(1,6) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong seconds in cdate1 = '//cdate1
      END IF

      IF ( irm(2,1) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong year    in cdate2 = '//cdate2
      END IF
      IF ( irm(2,2) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong month   in cdate2 = '//cdate2
      END IF
      IF ( irm(2,3) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong day     in cdate2 = '//cdate2
      END IF
      IF ( irm(2,4) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong hour    in cdate2 = '//cdate2
      END IF
      IF ( irm(2,5) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong minutes in cdate2 = '//cdate2
      END IF
      IF ( irm(2,6) > 0 ) THEN
        WRITE (*,*) 'ERROR diff_seconds: wrong seconds in cdate2 = '//cdate2
      END IF

    ENDIF

    !-------------------------------------------------------------------------------
  END SUBROUTINE diff_seconds

!===============================================================================
!==============================================================================

!------------------------------------------------------------------------------

  SUBROUTINE get_utc_date (ntsteps, ystartdate, dt, itype_calendar,           &
       yactdate1, yactdate2, nactday, acthour)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine determines the actual date of this forecast step.
    !   It also works for negative steps (= negative times).
    !   Different versions of historic calendars are implemented.
    !   (Gregorian / Proleptic calendar).
    !   It can be used for very long simulations over many years.
    !
    ! Method:
    !   Using the date of the forecast-start, the number of time steps
    !   already performed and the length of the time steps, the actual
    !   date is calculated taking leap-years into consideration.
    !   The date is given in three different formats.
    !
    ! Modules used:    NONE
    !
    !------------------------------------------------------------------------------
    !
    ! Input Parameter list:
    ! ---------------------

    INTEGER         ,   INTENT(IN)   ::                           &
         itype_calendar,   & ! for specifying the calendar used
         ntsteps             ! number of actual performed time-steps

    REAL   (KIND=dp),   INTENT(IN)   ::                           &
         dt         ! time step in seconds

    CHARACTER (LEN=14), INTENT(IN)   ::                           &
         ystartdate ! start date of the forecast

    ! Output Parameter list:
    ! ----------------------

    CHARACTER (LEN=14), INTENT(OUT)  ::                           &
         yactdate1  ! actual date in the form   yyyymmddhhmmss

    CHARACTER (LEN=28), INTENT(OUT)  ::                           &
         yactdate2  ! actual date in the form   wd   dd.mm.yy  hh mm ss UTC


    INTEGER         , INTENT(OUT)    ::                           &
         nactday    ! day of the year

    REAL   (KIND=dp), INTENT(OUT)    ::                           &
         acthour    ! actual hour of the day

    ! Local variables:
    INTEGER                          ::                                   &
         month(12), monthsum(13), ileap, iweek, iy, m,                       &
         idd, imm, iyy, ihh, imi, iss,                                       &
         iday, imonth, iyear, ihour, imin, isec,                             &
         immhours, iyyhours, iyear_hours, inhour, irsec

    CHARACTER (LEN=3)                :: yweek(7)

    ! And for computing the amount of seconds of the whole forecast time,
    ! an 8-Byte INTEGER has to be used. Otherwise the computation fails after
    ! approx. 68 years!!

    INTEGER (KIND=i8) :: insec

    !------------ End of header ---------------------------------------------------

    ! Begin subroutine get_utc_date

    DATA         month  / 31 ,  28 ,  31 ,  30 ,  31 ,  30 ,       &
         31 ,  31 ,  30 ,  31 ,  30 ,  31 /
    DATA         yweek  /'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /


    ! Statementfunction: ileap(yy) = 0:  no leap year,
    !                    ileap(yy) = 1:  leap year
    ! corrected version for Gregorian / Proleptic calendar
    ! by A. Dobler, CLM Community
    ileap (iy) = IABS( MOD(iy,  4) -   4) /   4  & ! every       4 years is a leapyear
         -IABS( MOD(iy,100) - 100) / 100  & ! every     100 years is no leapyear
         +IABS( MOD(iy,400) - 400) / 400    ! but every 400 years is a leapyear

    ! Divide ystartdate in day, month, year and hour
    ! and since version 4.24: also in minutes and seconds
    ! and calculate the sums of days from the beginning of the year to the
    ! end of the months
    READ ( ystartdate, '(I4,5I2)' ) iyy, imm, idd, ihh, imi, iss

    IF     (itype_calendar == 0) THEN
      month (2)    = 28 + ileap (iyy)
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + month(m-1)
      enddo
    ELSEIF (itype_calendar == 1) THEN
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + 30
      enddo
    ELSEIF (itype_calendar == 2) THEN
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + month(m-1)
      enddo
    ENDIF

    ! Determine how many hours have passed in this year
    ! iyyhours counts to the start of the forecast (if it is a full hour)
    ! or to the full hour before start of the forecast
    iyyhours = (idd*24) + monthsum(imm)*24 + (ihh-24)

    ! insec adds the seconds between full hour before start of the forecast and the
    ! start of the forecast (if any) + the forecast length in seconds
    insec = NINT(ntsteps*dt, i8) + INT(imi, i8) * 60_i8                 &
         + INT(iss, i8) *  1_i8
    ! inhour is the same in (full) hours
    IF (insec >= 0_i8) THEN
      inhour = INT(insec/3600_i8)
    ELSE
      inhour = INT(FLOOR(REAL(insec,dp)/3600.0_dp))
    ENDIF
    iyyhours = iyyhours + inhour

    ! Take turning of the year into account
    IF     (itype_calendar == 0) THEN
      iyear_hours = 8760 + ileap(iyy)*24
    ELSEIF (itype_calendar == 1) THEN
      iyear_hours = 8640
    ELSEIF (itype_calendar == 2) THEN
      iyear_hours = 8760
    ENDIF

    IF (iyyhours < 0) THEN
      iyear    = iyy-1
      IF     (itype_calendar == 0) THEN
        iyyhours = 8760 + ileap(iyear)*24 + iyyhours
      ELSEIF (itype_calendar == 1) THEN
        iyyhours = 8640                   + iyyhours
      ELSEIF (itype_calendar == 2) THEN
        iyyhours = 8760                   + iyyhours
      ENDIF
    ELSE IF (iyyhours >= iyear_hours) THEN
      ! Take also into account if the run lasts for several years
      iyear    = iyy
      IF     (itype_calendar == 0) THEN
        iyear_hours = 8760 + ileap(iyear)*24
      ELSEIF (itype_calendar == 1) THEN
        iyear_hours = 8640
      ELSEIF (itype_calendar == 2) THEN
        iyear_hours = 8760
      ENDIF

      DO WHILE (iyyhours >= iyear_hours)
        iyyhours = iyyhours - iyear_hours
        iyear=iyear+1
        IF     (itype_calendar == 0) THEN
          iyear_hours = 8760 + ileap(iyear)*24
        ELSEIF (itype_calendar == 1) THEN
          iyear_hours = 8640
        ELSEIF (itype_calendar == 2) THEN
          iyear_hours = 8760
        ENDIF
      ENDDO
    ELSE
      iyear    =   iyy
    ENDIF

    ! calculate monthsum for actual year
    IF     (itype_calendar == 0) THEN
      month (2)    = 28 + ileap (iyear)
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + month(m-1)
      enddo
    ELSEIF (itype_calendar == 1) THEN
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + 30
      enddo
    ELSEIF (itype_calendar == 2) THEN
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + month(m-1)
      enddo
    ENDIF

    ! Determine the actual date from iyyhours
    m        = 1
    immhours = iyyhours
    DO WHILE (immhours >= 0)
      m        = m+1
      immhours = iyyhours - monthsum(m) * 24
    ENDDO
    imonth   = m-1

    immhours = iyyhours - monthsum(imonth)*24
    iday     = immhours/24 + 1
    ihour    = MOD(immhours,24)

    ! irsec are the seconds left between the last full hour and the actual time
    irsec    = INT(MODULO(insec,  3600_i8))
    imin     = INT(irsec/60)
    isec     = MODULO(irsec,  60)
    ! add the time between the last full hour and the actual time
    acthour  = REAL (ihour, dp) +  REAL (irsec, dp) / 3600.0_dp
    nactday  = monthsum(imonth) + iday + INT(acthour/24.0_dp + 0.0001)
    iweek    = MOD(monthsum(imonth) + iday + (iyear-1901) + (iyear-1901)/4, 7)+1

    ihour    = INT(acthour)

    WRITE ( yactdate1(1:4)  , '(I4.4)' ) iyear
    WRITE ( yactdate1(5:6)  , '(I2.2)' ) imonth
    WRITE ( yactdate1(7:8)  , '(I2.2)' ) iday
    WRITE ( yactdate1(9:10) , '(I2.2)' ) ihour
    WRITE ( yactdate1(11:12), '(I2.2)' ) imin
    WRITE ( yactdate1(13:14), '(I2.2)' ) isec

    IF     (itype_calendar == 0 .OR. itype_calendar == 2) THEN
      yactdate2 = yweek(iweek)//' '//yactdate1(7:8)//'.'// yactdate1(5:6)//'.' &
           //yactdate1(1:4)//'  '//yactdate1(9:10)//':' &
           //yactdate1(11:12)//':'//yactdate1(13:14)//' UTC'
    ELSEIF (itype_calendar == 1) THEN
      yactdate2 = '    '//yactdate1(7:8)//'.'// yactdate1(5:6)//'.' &
           //yactdate1(1:4)//'  '//yactdate1(9:10)//':' &
           //yactdate1(11:12)//':'//yactdate1(13:14)//' UTC'
    ENDIF

  END SUBROUTINE get_utc_date

!===============================================================================

  ! Compute new datestr YYYYMMDDhhmmss from strings inidate and forrange:
  !  FORMATS: inidate  = YYYYMMDDhh or YYYYMMDDhhmmss
  !           forrange = real number [s], might as well be negative!
  ! Might not work for very long (climatic) forecast time ranges!
  ! Only for modern calendars (not Gregorian or Proleptic).
  
  FUNCTION new_datetime(inidate,forrange) RESULT(ct)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: inidate
    REAL(KIND=dp)             :: forrange
    CHARACTER(len=14)            :: ct
    CHARACTER(len=14)            :: ds

    INTEGER :: yy ,mon ,dd ,hh ,mm ,ss
    INTEGER :: yyi,moni,ddi,hhi,mmi,ssi
    INTEGER :: i, l, ndd(12)

    ndd = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

    ! read input time from input string:
    ds(:) = '0'
    ds = TRIM(ADJUSTL(inidate))
    READ(ds,'(i4,5(i2))') yyi,moni,ddi,hhi,mmi,ssi

    IF (is_leapyear(yyi)) ndd(2) = 29

    l = NINT(forrange)  ! in s

    IF (l >= 0) THEN

      ! compute days, hours, min, sec of the forecast time:
      dd = l / 86400
      hh = (l-dd*86400) / 3600
      mm = (l-dd*86400-hh*3600) / 60
      ss = MODULO(l,60)

      ! ... and add it to the input time, taking into account
      ! different month lengths and leap years:
      mm = mm + mmi + (ss+ssi)/60
      hh = hh + hhi +  mm/60
      dd = dd + ddi +  hh/24
      ss = MODULO(ss+ssi,60)
      mm = MODULO(mm,60)
      hh = MODULO(hh,24)

      mon = moni
      yy  = yyi
      DO
        IF (dd <= ndd(mon)) EXIT
        dd = dd - ndd(mon)
        mon = mon + 1
        IF (mon > 12) THEN
          yy = yy + 1
          IF (is_leapyear(yy)) THEN
            ndd(2) = 29
          ELSE
            ndd(2) = 28
          END IF
          mon = 1
        END IF
      END DO

    ELSE

      ! compute negative days, hours, min, sec of the negative forecast time:
      dd = l / 86400
      hh = (l-dd*86400) / 3600
      mm = (l-dd*86400-hh*3600) / 60
      ss = -MODULO(-l,60)

      ! ... and add it to the input time, taking into account
      ! different month lengths and leap years:
      mm = mm + mmi +  FLOOR((ss+ssi)/60.0_dp)
      hh = hh + hhi +  FLOOR(mm/60.0_dp)
      dd = dd + ddi +  FLOOR(hh/24.0_dp)
      ss = MODULO(ss+ssi,60)
      mm = MODULO(mm,60)
      hh = MODULO(hh,24)

      mon = moni
      yy  = yyi
      DO
        IF (dd > 0) EXIT
        mon = mon - 1
        IF (mon < 1) THEN
          yy = yy - 1
          IF (is_leapyear(yy)) THEN
            ndd(2) = 29
          ELSE
            ndd(2) = 28
          END IF
          mon = 12
        END IF
        dd = dd + ndd(mon)
      END DO

    END IF
      
    ! compose the output string:
    ct(:) = ' '
    WRITE (ct,'(i4.4,5(i2.2))') yy, mon, dd, hh, mm, ss

  END FUNCTION new_datetime

!===============================================================================

  FUNCTION is_leapyear(yyyy) RESULT(flag)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: yyyy
    LOGICAL             :: flag

    flag = (MOD(yyyy, 4) == 0) .AND. ( (MOD(yyyy,100) /= 0) .OR. (MOD(yyyy,400) == 0) )

  END FUNCTION is_leapyear

!===============================================================================

  FUNCTION jul2yyyymmdd (indate) RESULT (outdate)

    ! Convert julian number of day yyyyjul to "normal" date string yyyymmdd
    ! Example: cdate(1:8) = jul2yyyymmdd('2015205')    ! cdate = '20150724'

    CHARACTER(len=*), INTENT(in)  :: indate    ! YYYYJUL
    CHARACTER(len=8)              :: outdate   ! YYYYMMDD
    INTEGER      :: i, monsum(12), year, tmp, julday, mm, dd

    IF (LEN_TRIM(indate) < 7) THEN
      WRITE (*,*) 'jul2yyyymmdd(): Wrong input date string (has to be YYYYJUL). Wrong string is: '//TRIM(indate)
      outdate = 'YYYYMMDD'
    ELSE
      READ (indate(1:4), *) year
      READ (indate(5:7), *) julday
      IF (is_leapyear(year)) THEN
        monsum = (/ 0,31,60,91,121,152,182,213,244,274,305,335 /)
      ELSE
        monsum = (/ 0,31,59,90,120,151,181,212,243,273,304,334 /)
      END IF
      i = 12
      DO WHILE (julday <= monsum(i) .AND. i > 1)
        i = i - 1
      END DO
      mm = i
      dd = MAX(julday, 1) - monsum(i)
      WRITE (outdate, '(i4.4,2i2.2)') year, mm, dd

    END IF

  END FUNCTION jul2yyyymmdd

!===============================================================================

  ! Round datestr YYYYMMDDhhmmss to the next full time interval of minutes (hour-synchronized).
  !  FORMATS: indate  = YYYYMMDDhh or YYYYMMDDhhmmss
  !           round_min = int number [min]  positive means round up, negative means round down
  ! NOTE: if indate is exactly equal to a rounding interval and it is to be rounded up (not down!),
  !       it will be rounded up to the next interval!
  FUNCTION round_datestr_min(indate,round_min) RESULT(ct)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: indate
    INTEGER                      :: round_min
    CHARACTER(len=14)            :: ct
    CHARACTER(len=14)            :: ds

    INTEGER :: yy ,mon ,dd ,hh ,mm ,ss
    INTEGER :: yyi,moni,ddi,hhi,mmi,ssi
    INTEGER :: i, l, ndd(12)

    ndd = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

    ! Make sure that round_min is smaller or equal to 60 and a divisor of 60:
    IF ( ABS(round_min) > 60 .OR. round_min == 0 .OR. MODULO(60,MAX(ABS(round_min),1)) /= 0 ) THEN
      ct(:) = '0'
      ct = TRIM(indate)
      WRITE (*,*) 'ERROR round_datestr_min: Invalid rounding interval round_min = ', round_min
      RETURN
    END IF

    ! Read input time from input string:
    ds(:) = '0'
    ds = TRIM(ADJUSTL(indate))
    READ(ds,'(i4,5(i2))') yyi,moni,ddi,hhi,mmi,ssi

    IF (is_leapyear(yyi)) ndd(2) = 29

    ! round the input minutes to full round_min interval, taking into account
    ! different month lengths and leap years:
    ss = 0
    mm = mmi - MODULO(mmi,ABS(round_min)) + MIN(SIGN(1,round_min)+1,1)*ABS(round_min)
    hh = hhi + mm/60
    dd = ddi + hh/24
    mm = MODULO(mm,60)
    hh = MODULO(hh,24)

    mon = moni
    yy  = yyi
    DO
      IF (dd <= ndd(mon)) EXIT
      dd = dd - ndd(mon)
      mon = mon + 1
      IF (mon > 12) THEN
        yy = yy + 1
        IF (is_leapyear(yy)) THEN
          ndd(2) = 29
        ELSE
          ndd(2) = 28
        END IF
        mon = 1
      END IF
    END DO

    ! Compose the output string:
    ct(:) = ' '
    WRITE (ct,'(i4.4,5(i2.2))') yy, mon, dd, hh, mm, ss

  END FUNCTION round_datestr_min

!===============================================================================

  SUBROUTINE get_filenames_in_directory(indir, strlen, fnames, nfiles, ierror, proc_id, dirlistingfile)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)  :: indir
    INTEGER, INTENT(in)           :: strlen
    CHARACTER(len=strlen), ALLOCATABLE, INTENT(out) :: fnames(:)
    INTEGER, INTENT(out) :: nfiles, ierror
    INTEGER, INTENT(in), OPTIONAL :: proc_id   ! processor identifier, if called in a parallel program
    CHARACTER(len=*), INTENT(in), OPTIONAL  :: dirlistingfile

    INTEGER :: i, funit
    CHARACTER(len=500)  :: tmpfile
    CHARACTER(len=2500) :: commandline
    CHARACTER(len=500)  :: errmsg

    ! if an explicit filename for a file containing the directory listing was given (optional),
    !  just try to read this file. Otherwise, use the local tmpfile from above and delete it after
    !  it has been read:

    IF (PRESENT(dirlistingfile)) THEN
      tmpfile(:) = ' '
      tmpfile = TRIM(dirlistingfile)
    ELSE
      tmpfile(:) = ' '
      IF (PRESENT(proc_id)) THEN
        WRITE (tmpfile, '(a,i5.5,a)') 'tmpdirlist', proc_id, '.txt'
      ELSE
        tmpfile = 'tmpdirlist.txt'
      END IF

!!$      IF (LEN_TRIM(ADJUSTL(indir)) == 0) THEN
!!$
!!$        ! Listing of regular files and links that do not point to directories:
!!$        commandline(:) = ' '
!!$        commandline = 'ls -pL1 ./ | grep -v / | wc -l > '//TRIM(tmpfile)//'; '// &
!!$                      'ls -pL1 ./ | grep -v /  >> '//TRIM(tmpfile)
!!$      ELSE
!!$
!!$        ! Listing of regular files and links that do not point to directories:
!!$        commandline(:) = ' '
!!$        commandline = 'ls -pL1 '//TRIM(ADJUSTL(indir))//'/'//' | grep -v / | wc -l > '//TRIM(tmpfile)//'; '// &
!!$                      'ls -pL1 '//TRIM(ADJUSTL(indir))//'/'//' | grep -v /  >> '//TRIM(tmpfile)
!!$      END IF

      ! Listing of regular files and links that do not point to directories:
      commandline(:) = ' '
      commandline = 'ls -pL1 '//TRIM(ADJUSTL(indir))//' | grep -v / | wc -l > '//TRIM(tmpfile)//'; '// &
                    'ls -pL1 '//TRIM(ADJUSTL(indir))//' | grep -v /  >> '//TRIM(tmpfile)
      errmsg(:) = ' '
      CALL EXECUTE_COMMAND_LINE( TRIM(commandline), wait=.TRUE., exitstat=ierror, cmdmsg=errmsg )
      IF (ierror /= 0) THEN
        WRITE (*, '(a,/,a,/,a,i7)') 'ERROR get_filenames_in_directory(): '// &
             'an "ls" via EXECUTE_COMMAND_LINE() returned the error message', &
             '"'//TRIM(errmsg)//'"', 'and error code ', ierror
      END IF
    END IF

    CALL get_free_funit ( funit )
    OPEN(funit, file=TRIM(tmpfile), form='formatted', status='old', iostat=ierror)
    IF (ierror /= 0) THEN
      WRITE (*, '(a)') 'WARNING get_filenames_in_directory(): '//TRIM(tmpfile)//' could not be opened!'
      nfiles = 0
    ELSE
      READ (funit, *, iostat=ierror) nfiles
      IF (ierror /= 0) THEN
        WRITE (*, '(a)') 'WARNING get_filenames_in_directory(): problem when reading number of files from '//TRIM(tmpfile)
        nfiles = 0
      ELSE
        ALLOCATE( fnames(nfiles) )
        fnames(:)(:) = ' '
        DO i=1, nfiles
          READ(funit, *, iostat=ierror) fnames(i)
          IF (ierror /= 0) THEN
            WRITE (*, '(a)') 'WARNING get_filenames_in_directory(): problem when reading filenames from '//TRIM(tmpfile)
            nfiles = 0
          END IF
        END DO
        CLOSE(funit)
      END IF
    END IF
    IF (.NOT.PRESENT(dirlistingfile)) THEN
      CALL EXECUTE_COMMAND_LINE( 'rm '//TRIM(tmpfile) )
    END IF

  END SUBROUTINE get_filenames_in_directory

  SUBROUTINE get_filenames_by_pattern(indir, pattern, strlen, fnames, nfiles, ierror, proc_id)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)  :: indir, pattern
    INTEGER, INTENT(in)           :: strlen
    CHARACTER(len=strlen), ALLOCATABLE, INTENT(out) :: fnames(:)
    INTEGER, INTENT(out) :: nfiles, ierror
    INTEGER, INTENT(in), OPTIONAL :: proc_id   ! processor identifier, if called in a parallel program

    INTEGER :: i, funit, ierr
    CHARACTER(len=25)  :: tmpfile1, tmpfile2
    CHARACTER(len=2500) :: commandline
    CHARACTER(len=500) :: errmsg

    tmpfile1(:) = ' '
    tmpfile2(:) = ' '
    IF (PRESENT(proc_id)) THEN
      WRITE (tmpfile1, '(a,i5.5,a)') 'tmpdirlist1', proc_id, '.txt'
      WRITE (tmpfile2, '(a,i5.5,a)') 'tmpdirlist2', proc_id, '.txt'
    ELSE
      tmpfile1 = 'tmpdirlist1.txt'
      tmpfile2 = 'tmpdirlist2.txt'
    END IF

    commandline(:) = ' '
    commandline = 'for fi in $(ls -pLd1 '//TRIM(indir)//'/'//TRIM(pattern)//' | grep -v "/$"); do '// &
         'echo $(basename $fi); done > '//TRIM(tmpfile1)//'; '// &
         'echo $(cat '//TRIM(tmpfile1)//' | wc -l) > '//TRIM(tmpfile2)//'; cat '//TRIM(tmpfile1)//' >> '//TRIM(tmpfile2)

    ! Listing of regular files and links that do not point to directories:
    CALL EXECUTE_COMMAND_LINE( TRIM(commandline), wait=.TRUE., exitstat=ierror, cmdmsg=errmsg )
    IF (ierror /= 0) THEN
      WRITE (*, '(a,/,a)') 'ERROR get_filenames_by_pattern(): '// &
           'an "ls" via EXECUTE_COMMAND_LINE() returned the error message', &
           '"'//TRIM(errmsg)//'"'
      nfiles = 0
    ELSE
      CALL get_free_funit ( funit )
      OPEN(funit, file=TRIM(tmpfile2), form='formatted', status='old', iostat=ierror)
      IF (ierror /= 0) THEN
        WRITE (*, '(a)') 'ERROR in get_filenames_by_pattern: '//TRIM(tmpfile2)//' was not found!'
        nfiles = 0
      ELSE
        READ (funit, *) nfiles
        ALLOCATE( fnames(nfiles) )
        fnames(:)(:) = ' '
        DO i=1, nfiles
          READ(funit, *) fnames(i)
        END DO
        CLOSE(funit)
      END IF
      CALL EXECUTE_COMMAND_LINE( 'rm '//TRIM(tmpfile1)//' '//TRIM(tmpfile2) )
    END IF

  END SUBROUTINE get_filenames_by_pattern

!===============================================================================

  SUBROUTINE split_string (str, sep, wordlen, words, nwords)
    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine splits a string "str" at the separator string(s) "sep"
    !   and returns the separate substrings in the allocatable string array (words(:)).
    !   Optionally, the number of words "nwords" is returned.
    !
    ! Example:  The call
    !
    !     CHARACTER(len=10), ALLOCATABLE, DIMENSION(:) :: words
    !
    !     CALL split_string ('uliga__jutta__oliverb__mama','__', LEN(words), words, nwords)
    !
    ! allocates words(4) and fills it with
    !
    !     words(1) = 'uliga     '
    !     words(2) = 'jutta     '
    !     words(3) = 'oliverb   '
    !     words(4) = 'mama      '
    !
    !  ==> the words are padded with spaces up to their length.
    !
    ! whereas
    !     CALL split_string ('__uliga__jutta__oliverb__mama__','__', LEN(words), words, nwords)
    ! returns
    !     words(1) = '          '
    !     words(2) = 'uliga     '
    !     words(3) = 'jutta     '
    !     words(4) = 'oliverb   '
    !     words(5) = 'mama      '
    !     words(6) = '          '
    !
    ! ==> "sep" at the beginning or end of "str" leads to empty words.
    !
    ! The call
    !     CALL split_string ('uliga','__', LEN(words), words, nwords)
    ! returns
    !     words(1) = 'uliga     '
    !
    ! ==> if "sep" is not in "str", str itself is returned.
    !
    ! The call with the empty string
    !     CALL split_string ('','__', LEN(words), words, nwords)
    ! returns
    !     words(1) = '          '
    !
    ! Modules used:    NONE
    !
    !------------------------------------------------------------------------------
    !
    ! Input Parameter list:
    ! ---------------------

    CHARACTER(len=*), INTENT(in)  :: str, sep
    INTEGER, INTENT(in)           :: wordlen

    ! Output Parameter list:
    ! ---------------------

    CHARACTER(len=wordlen), ALLOCATABLE, INTENT(inout) :: words(:)
    INTEGER, INTENT(out), OPTIONAL                   :: nwords

    INTEGER :: lensep, lenstr, n, pos1, pos2

    lenstr = LEN(str)
    lensep = LEN(sep)

    ! Determine how often "sep" occurs in "str":
    n = 0
    pos1 = 1
    DO
      IF (str(pos1:pos1+lensep-1) == sep) THEN
        n = n + 1
        pos1 = pos1 + lensep
      ELSE
        pos1 = pos1 + 1
      END IF
      IF (pos1 > lenstr) EXIT
    END DO

   ! Allocate the string array for the (n+1) separate words:
    IF (ALLOCATED(words)) DEALLOCATE(words)
    ALLOCATE(words(n+1))
    words(:)(:) = ' '
    IF (PRESENT(nwords)) nwords = n+1

    n = 0
    pos1 = 1
    DO
      pos2 = INDEX(str(pos1:), sep)
      IF (pos2 == 0) THEN
        n = n + 1
        words(n) = str(pos1:)
        EXIT
      END IF
      n = n + 1
      words(n) = str(pos1:pos1+pos2-2)
      pos1 = pos1 + pos2 + lensep - 1
      IF (pos1 > lenstr) EXIT
    END DO

  END SUBROUTINE split_string

!===============================================================================
! String conversion from/to lower/upper case: (NOTE: Umlauts are NOT converted!)
!===============================================================================

  FUNCTION toupper (c) RESULT (cu)

    IMPLICIT NONE

    CHARACTER(len=*) :: c
    CHARACTER(len=LEN_TRIM(c)) :: cu
    INTEGER :: i, lacode, uacode, lzcode, icode 

    lacode=IACHAR('a')
    lzcode=IACHAR('z')
    uacode=IACHAR('A')

    cu = c
    DO i = 1, LEN_TRIM(c)
      icode = IACHAR(c(i:i))
      IF (icode >= lacode .AND. icode <= lzcode) THEN
        cu(i:i) = ACHAR(icode - lacode + uacode)
      END IF
    END DO

  END FUNCTION toupper

  FUNCTION tolower (c) RESULT (cl)

    IMPLICIT NONE

    CHARACTER(len=*) :: c
    CHARACTER(len=LEN(c)) :: cl
    INTEGER :: i, lacode, uacode, uzcode, icode

    lacode=IACHAR('a')
    uacode=IACHAR('A')
    uzcode=IACHAR('Z')

    cl = c
    DO i = 1, LEN_TRIM(c)
      icode = IACHAR(c(i:i))
      IF (icode >= uacode .AND. icode <= uzcode) THEN
        cl(i:i) = ACHAR(icode - uacode + lacode)
      END IF
    END DO

  END FUNCTION tolower

!===============================================================================
  
!===============================================================================
! Routines for the overloaded model procedure init_vari(field, initval)
!===============================================================================

  SUBROUTINE init_vari_dp_5D(vari, val)

    REAL(kind=dp), INTENT(inout) :: vari(:,:,:,:,:)
    REAL(kind=dp), INTENT(in)    :: val
    INTEGER :: i, j, k, l, m

    INTEGER :: lb(5), ub(5)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k,l)
!$omp parallel do private(i,j,k,l,m)
    DO m=lb(5), ub(5)
      DO l=lb(4), ub(4)
        DO k=lb(3), ub(3)
          DO j=lb(2), ub(2)
            DO i=lb(1), ub(1)
              vari(i,j,k,l,m) = val
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_dp_5D

  SUBROUTINE init_vari_sp_5D(vari, val)

    REAL(kind=sp), INTENT(inout) :: vari(:,:,:,:,:)
    REAL(kind=sp), INTENT(in)    :: val
    INTEGER :: i, j, k, l, m

    INTEGER :: lb(5), ub(5)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k,l)
!$omp parallel do private(i,j,k,l,m)
    DO m=lb(5), ub(5)
      DO l=lb(4), ub(4)
        DO k=lb(3), ub(3)
          DO j=lb(2), ub(2)
            DO i=lb(1), ub(1)
              vari(i,j,k,l,m) = val
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_sp_5D

  SUBROUTINE init_vari_int_5D(vari, val)

    INTEGER, INTENT(inout) :: vari(:,:,:,:,:)
    INTEGER, INTENT(in)    :: val
    INTEGER :: i, j, k, l, m

    INTEGER :: lb(5), ub(5)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k,l)
!$omp parallel do  private(i,j,k,l)
    DO m=lb(5), ub(5)
      DO l=lb(4), ub(4)
        DO k=lb(3), ub(3)
          DO j=lb(2), ub(2)
            DO i=lb(1), ub(1)
              vari(i,j,k,l,m) = val
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_int_5D

  SUBROUTINE init_vari_dp_4D(vari, val)

    REAL(kind=dp), INTENT(inout) :: vari(:,:,:,:)
    REAL(kind=dp), INTENT(in)    :: val
    INTEGER :: i, j, k, l

    INTEGER :: lb(4), ub(4)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k,l)
!$omp parallel do private(i,j,k,l)
    DO l=lb(4), ub(4)
      DO k=lb(3), ub(3)
        DO j=lb(2), ub(2)
          DO i=lb(1), ub(1)
            vari(i,j,k,l) = val
          END DO
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_dp_4D

  SUBROUTINE init_vari_sp_4D(vari, val)

    REAL(kind=sp), INTENT(inout) :: vari(:,:,:,:)
    REAL(kind=sp), INTENT(in)    :: val
    INTEGER :: i, j, k, l

    INTEGER :: lb(4), ub(4)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k,l)
!$omp parallel do private(i,j,k,l)
    DO l=lb(4), ub(4)
      DO k=lb(3), ub(3)
        DO j=lb(2), ub(2)
          DO i=lb(1), ub(1)
            vari(i,j,k,l) = val
          END DO
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_sp_4D

  SUBROUTINE init_vari_int_4D(vari, val)

    INTEGER, INTENT(inout) :: vari(:,:,:,:)
    INTEGER, INTENT(in)    :: val
    INTEGER :: i, j, k, l

    INTEGER :: lb(4), ub(4)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k,l)
!$omp parallel do  private(i,j,k,l)
    DO l=lb(4), ub(4)
      DO k=lb(3), ub(3)
        DO j=lb(2), ub(2)
          DO i=lb(1), ub(1)
            vari(i,j,k,l) = val
          END DO
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_int_4D

  SUBROUTINE init_vari_dp_3D(vari, val)

    REAL(kind=dp), INTENT(inout) :: vari(:,:,:)
    REAL(kind=dp), INTENT(in)    :: val
    INTEGER :: i, j, k

    INTEGER :: lb(3), ub(3)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k)
!$omp parallel do private(i,j,k)
    DO k=lb(3), ub(3)
      DO j=lb(2), ub(2)
        DO i=lb(1), ub(1)
          vari(i,j,k) = val
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_dp_3D

  SUBROUTINE init_vari_sp_3D(vari, val)

    REAL(kind=sp), INTENT(inout) :: vari(:,:,:)
    REAL(kind=sp), INTENT(in)    :: val
    INTEGER :: i, j, k

    INTEGER :: lb(3), ub(3)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k)
!$omp parallel do private(i,j,k)
    DO k=lb(3), ub(3)
      DO j=lb(2), ub(2)
        DO i=lb(1), ub(1)
          vari(i,j,k) = val
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_sp_3D

  SUBROUTINE init_vari_int_3D(vari, val)

    INTEGER, INTENT(inout) :: vari(:,:,:)
    INTEGER, INTENT(in)    :: val
    INTEGER :: i, j, k

    INTEGER :: lb(3), ub(3)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k)
!$omp parallel do  private(i,j,k)
    DO k=lb(3), ub(3)
      DO j=lb(2), ub(2)
        DO i=lb(1), ub(1)
          vari(i,j,k) = val
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_int_3D

  SUBROUTINE init_vari_dp_2D(vari, val)

    REAL(kind=dp), INTENT(inout) :: vari(:,:)
    REAL(kind=dp), INTENT(in)    :: val
    INTEGER :: i, j

    INTEGER :: lb(2), ub(2)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(2) private(i,j)
!$omp parallel do private(i,j)
    DO j=lb(2), ub(2)
      DO i=lb(1), ub(1)
        vari(i,j) = val
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_dp_2D

  SUBROUTINE init_vari_sp_2D(vari, val)

    REAL(kind=sp), INTENT(inout) :: vari(:,:)
    REAL(kind=sp), INTENT(in)    :: val
    INTEGER :: i, j

    INTEGER :: lb(2), ub(2)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(2) private(i,j)
!$omp parallel do private(i,j)
    DO j=lb(2), ub(2)
      DO i=lb(1), ub(1)
        vari(i,j) = val
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_sp_2D

  SUBROUTINE init_vari_int_2D(vari, val)

    INTEGER, INTENT(inout) :: vari(:,:)
    INTEGER, INTENT(in)    :: val
    INTEGER :: i, j

    INTEGER :: lb(2), ub(2)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(2) private(i,j)
!$omp parallel do private(i,j)
    DO j=lb(2), ub(2)
      DO i=lb(1), ub(1)
        vari(i,j) = val
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_int_2D

  SUBROUTINE init_vari_dp_1D(vari, val)

    REAL(kind=dp), INTENT(inout) :: vari(:)
    REAL(kind=dp), INTENT(in)    :: val
    INTEGER :: i

    INTEGER :: lb(1), ub(1)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!$omp parallel do private(i)
    DO i=lb(1), ub(1)
      vari(i) = val
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_dp_1D

  SUBROUTINE init_vari_sp_1D(vari, val)

    REAL(kind=sp), INTENT(inout) :: vari(:)
    REAL(kind=sp), INTENT(in)    :: val
    INTEGER :: i

    INTEGER :: lb(1), ub(1)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!$omp parallel do private(i)
    DO i=lb(1), ub(1)
      vari(i) = val
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_sp_1D

  SUBROUTINE init_vari_int_1D(vari, val)

    INTEGER, INTENT(inout) :: vari(:)
    INTEGER, INTENT(in)    :: val
    INTEGER :: i

    INTEGER :: lb(1), ub(1)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!$omp parallel do private(i)
    DO i=lb(1), ub(1)
      vari(i) = val
    END DO
!$omp end parallel do

  END SUBROUTINE init_vari_int_1D

!===============================================================================

!===============================================================================
! Routines for the overloaded model procedure set_missing_and_correct0(field)
!===============================================================================

  SUBROUTINE set_missing_and_correct0_dp_3D(vari)

    REAL(kind=dp), INTENT(inout) :: vari(:,:,:)
    INTEGER :: i, j, k

    INTEGER :: lb(3), ub(3)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k)
!$omp parallel do private(i,j,k)
    DO k=lb(3), ub(3)
      DO j=lb(2), ub(2)
        DO i=lb(1), ub(1)
          IF (vari(i,j,k) >= miss_threshold .AND. vari(i,j,k) < dBZ_crit_radar) THEN
            vari(i,j,k) = zero_value
          END IF
          IF (vari(i,j,k) < miss_threshold) THEN
            vari(i,j,k) = miss_value
          END IF
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE set_missing_and_correct0_dp_3D

  SUBROUTINE set_missing_and_correct0_sp_3D(vari)

    REAL(kind=sp), INTENT(inout) :: vari(:,:,:)
    INTEGER :: i, j, k

    INTEGER :: lb(3), ub(3)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k)
!$omp parallel do private(i,j,k)
    DO k=lb(3), ub(3)
      DO j=lb(2), ub(2)
        DO i=lb(1), ub(1)
          IF (vari(i,j,k) >= missthr_sp .AND. vari(i,j,k) < dBZ_crit_radar_sp) THEN
            vari(i,j,k) = zeroval_sp
          END IF
          IF (vari(i,j,k) < missthr_sp) THEN
            vari(i,j,k) = missval_sp
          END IF
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE set_missing_and_correct0_sp_3D

  SUBROUTINE set_missing_and_correct0_dp_1D(vari)

    REAL(kind=dp), INTENT(inout) :: vari(:)
    INTEGER :: i

    INTEGER :: lb(1), ub(1)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!$omp parallel do private(i)
    DO i=lb(1), ub(1)
      IF (vari(i) >= miss_threshold .AND. vari(i) < dBZ_crit_radar) THEN
        vari(i) = zero_value
      END IF
      IF (vari(i) < miss_threshold) THEN
        vari(i) = miss_value
      END IF
    END DO
!$omp end parallel do

  END SUBROUTINE set_missing_and_correct0_dp_1D

  SUBROUTINE set_missing_and_correct0_sp_1D(vari)

    REAL(kind=sp), INTENT(inout) :: vari(:)
    INTEGER :: i

    INTEGER :: lb(1), ub(1)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!$omp parallel do private(i)
    DO i=lb(1), ub(1)
      IF (vari(i) >= missthr_sp .AND. vari(i) < dBZ_crit_radar_sp) THEN
        vari(i) = zeroval_sp
      END IF
      IF (vari(i) < missthr_sp) THEN
        vari(i) = missval_sp
      END IF
    END DO
!$omp end parallel do

  END SUBROUTINE set_missing_and_correct0_sp_1D

!===============================================================================

!===============================================================================
! Routines for the overloaded model procedure dbz_to_linear(field)
!===============================================================================

  SUBROUTINE dbz_to_linear_dp_3D(vari)

    REAL(kind=dp), INTENT(inout) :: vari(:,:,:)
    INTEGER :: i, j, k

    INTEGER :: lb(3), ub(3)
    REAL(kind=dp), PARAMETER :: ln10_o10 = LOG(10.0_dp)*0.1_dp

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k)
!$omp parallel do private(i,j,k)
    DO k=lb(3), ub(3)
      DO j=lb(2), ub(2)
        DO i=lb(1), ub(1)
          IF (vari(i,j,k) > dBZ_crit_radar) THEN
            vari(i,j,k) = EXP(vari(i,j,k)*ln10_o10)
          ELSE IF (vari(i,j,k) > miss_threshold) THEN
            vari(i,j,k) = 0.0_dp
          END IF
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE dbz_to_linear_dp_3D

  SUBROUTINE dbz_to_linear_sp_3D(vari)

    REAL(kind=sp), INTENT(inout) :: vari(:,:,:)
    INTEGER :: i, j, k

    INTEGER :: lb(3), ub(3)
    REAL(kind=sp), PARAMETER :: ln10_o10 = LOG(10.0_sp)*0.1_sp

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k)
!$omp parallel do private(i,j,k)
    DO k=lb(3), ub(3)
      DO j=lb(2), ub(2)
        DO i=lb(1), ub(1)
          IF (vari(i,j,k) > dBZ_crit_radar_sp) THEN
            vari(i,j,k) = EXP(vari(i,j,k)*ln10_o10)
          ELSE IF (vari(i,j,k) > missthr_sp) THEN
            vari(i,j,k) = 0.0_dp
          END IF
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE dbz_to_linear_sp_3D

  SUBROUTINE dbz_to_linear_dp_1D(vari)

    REAL(kind=dp), INTENT(inout) :: vari(:)
    INTEGER :: i

    INTEGER :: lb(1), ub(1)
    REAL(kind=dp), PARAMETER :: ln10_o10 = LOG(10.0_dp)*0.1_dp

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!$omp parallel do private(i)
    DO i=lb(1), ub(1)
      IF (vari(i) > dBZ_crit_radar) THEN
        vari(i) = EXP(vari(i)*ln10_o10)
      ELSE IF (vari(i) > miss_threshold) THEN
        vari(i) = 0.0_dp
      END IF
    END DO
!$omp end parallel do

  END SUBROUTINE dbz_to_linear_dp_1D

  SUBROUTINE dbz_to_linear_sp_1D(vari)

    REAL(kind=sp), INTENT(inout) :: vari(:)
    INTEGER :: i

    INTEGER :: lb(1), ub(1)
    REAL(kind=sp), PARAMETER :: ln10_o10 = LOG(10.0_sp)*0.1_sp

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!$omp parallel do private(i)
    DO i=lb(1), ub(1)
      IF (vari(i) > dBZ_crit_radar_sp) THEN
        vari(i) = EXP(vari(i)*ln10_o10)
      ELSE IF (vari(i) > missthr_sp) THEN
        vari(i) = 0.0_sp
     END IF
    END DO
!$omp end parallel do

  END SUBROUTINE dbz_to_linear_sp_1D

  SUBROUTINE dbz_to_linear_dp_scal(vari)

    REAL(kind=dp), INTENT(inout) :: vari

    REAL(kind=dp), PARAMETER :: ln10_o10 = LOG(10.0_dp)*0.1_dp

    IF (vari > dBZ_crit_radar) THEN
      vari = EXP(vari*ln10_o10)
    ELSE IF (vari > miss_threshold) THEN
      vari = 0.0_dp
    END IF

  END SUBROUTINE dbz_to_linear_dp_scal

  SUBROUTINE dbz_to_linear_sp_scal(vari)

    REAL(kind=sp), INTENT(inout) :: vari

    REAL(kind=sp), PARAMETER :: ln10_o10 = LOG(10.0_sp)*0.1_sp

    IF (vari > dBZ_crit_radar) THEN
      vari = EXP(vari*ln10_o10)
    ELSE IF (vari > miss_threshold) THEN
      vari = 0.0_sp
    END IF

  END SUBROUTINE dbz_to_linear_sp_scal

!===============================================================================

!===============================================================================
! Routines for the overloaded model procedure linear_to_dbz(field)
!===============================================================================

  SUBROUTINE linear_to_dbz_dp_3D(vari)

    REAL(kind=dp), INTENT(inout) :: vari(:,:,:)
    INTEGER :: i, j, k

    INTEGER :: lb(3), ub(3)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k)
!$omp parallel do private(i,j,k)
    DO k=lb(3), ub(3)
      DO j=lb(2), ub(2)
        DO i=lb(1), ub(1)
          IF (vari(i,j,k) >= Z_crit_radar) THEN
            vari(i,j,k) = 10.0_dp * LOG10(vari(i,j,k))
          ELSE IF (vari(i,j,k) >= miss_threshold) THEN
            vari(i,j,k) = zero_value
          END IF
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE linear_to_dbz_dp_3D

  SUBROUTINE linear_to_dbz_sp_3D(vari)

    REAL(kind=sp), INTENT(inout) :: vari(:,:,:)
    INTEGER :: i, j, k

    INTEGER :: lb(3), ub(3)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!!$omp parallel do collapse(3) private(i,j,k)
!$omp parallel do private(i,j,k)
    DO k=lb(3), ub(3)
      DO j=lb(2), ub(2)
        DO i=lb(1), ub(1)
          IF (vari(i,j,k) >= Z_crit_radar_sp) THEN
            vari(i,j,k) = 10.0_sp * LOG10(vari(i,j,k))
          ELSE IF (vari(i,j,k) >= missthr_sp) THEN
            vari(i,j,k) = zeroval_sp
          END IF
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE linear_to_dbz_sp_3D

  SUBROUTINE linear_to_dbz_dp_1D(vari)

    REAL(kind=dp), INTENT(inout) :: vari(:)
    INTEGER :: i

    INTEGER :: lb(1), ub(1)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!$omp parallel do private(i)
    DO i=lb(1), ub(1)
      IF (vari(i) >= Z_crit_radar) THEN
        vari(i) = 10.0_dp * LOG10(vari(i))
      ELSE IF (vari(i) >= miss_threshold) THEN
        vari(i) = zero_value
      END IF
    END DO
!$omp end parallel do

  END SUBROUTINE linear_to_dbz_dp_1D

  SUBROUTINE linear_to_dbz_sp_1D(vari)

    REAL(kind=sp), INTENT(inout) :: vari(:)
    INTEGER :: i

    INTEGER :: lb(1), ub(1)

    lb = LBOUND(vari)
    ub = UBOUND(vari)
!$omp parallel do private(i)
    DO i=lb(1), ub(1)
      IF (vari(i) >= Z_crit_radar_sp) THEN
        vari(i) = 10.0_sp * LOG10(vari(i))
      ELSE IF (vari(i) >= missthr_sp) THEN
        vari(i) = zeroval_sp
      END IF
    END DO
!$omp end parallel do

  END SUBROUTINE linear_to_dbz_sp_1D

  SUBROUTINE linear_to_dbz_dp_scal(vari)

    REAL(kind=dp), INTENT(inout) :: vari

    IF (vari >= Z_crit_radar) THEN
      vari = 10.0_dp * LOG10(vari)
    ELSE IF (vari >= miss_threshold) THEN
      vari = zero_value
    END IF

  END SUBROUTINE linear_to_dbz_dp_scal

  SUBROUTINE linear_to_dbz_sp_scal(vari)

    REAL(kind=sp), INTENT(inout) :: vari

    IF (vari >= Z_crit_radar) THEN
      vari = 10.0_sp * LOG10(vari)
    ELSE IF (vari >= miss_threshold) THEN
      vari = zero_value
    END IF

  END SUBROUTINE linear_to_dbz_sp_scal

!===============================================================================

!===============================================================================
! Routines for the overloaded model procedure bubblesort(vec)
!===============================================================================

  SUBROUTINE bubblesort_int (a)

    INTEGER, INTENT(inout) :: a(:)

    INTEGER :: i, n, tmp
    LOGICAL :: swapped
    
    n = SIZE(a)
    DO
      swapped = .FALSE.
      DO i=1, n-1
        IF (a(i) > a(i+1)) THEN
          tmp    = a(i)
          a(i)   = a(i+1)
          a(i+1) = tmp
          swapped = .TRUE.
        END IF
      END DO
      n = n-1
      IF (.NOT.swapped) EXIT
    END DO
  END SUBROUTINE bubblesort_int
  
  SUBROUTINE bubblesort_dp (a)

    REAL(kind=dp), INTENT(inout) :: a(:)

    INTEGER :: i, n, tmp
    LOGICAL :: swapped
    
    n = SIZE(a)
    DO
      swapped = .FALSE.
      DO i=1, n-1
        IF (a(i) > a(i+1)) THEN
          tmp    = a(i)
          a(i)   = a(i+1)
          a(i+1) = tmp
          swapped = .TRUE.
        END IF
      END DO
      n = n-1
      IF (.NOT.swapped) EXIT
    END DO

  END SUBROUTINE bubblesort_dp

  SUBROUTINE bubblesort_sp (a)

    REAL(kind=sp), INTENT(inout) :: a(:)

    INTEGER :: i, n, tmp
    LOGICAL :: swapped
    
    n = SIZE(a)
    DO
      swapped = .FALSE.
      DO i=1, n-1
        IF (a(i) > a(i+1)) THEN
          tmp    = a(i)
          a(i)   = a(i+1)
          a(i+1) = tmp
          swapped = .TRUE.
        END IF
      END DO
      n = n-1
      IF (.NOT.swapped) EXIT
    END DO

  END SUBROUTINE bubblesort_sp

#ifdef NEED_ECL_FALLBACK
  SUBROUTINE execute_command_line( command, wait, exitstat, cmdstat, cmdmsg )
    CHARACTER(*), INTENT(IN) :: command
    LOGICAL, INTENT(IN), OPTIONAL :: wait
    INTEGER, INTENT(OUT), OPTIONAL :: exitstat, cmdstat
    CHARACTER(*), INTENT(OUT), OPTIONAL :: cmdmsg

    INTEGER :: e, SYSTEM
    EXTERNAL :: SYSTEM

    IF (PRESENT(wait) .AND. .NOT.wait) THEN
      IF (PRESENT(exitstat)) exitstat = 1
      IF (PRESENT(cmdstat)) cmdstat = 1
      IF (PRESENT(cmdmsg)) cmdmsg = 'asynchronous execution is not supported'
      RETURN
    END IF

    e = SYSTEM(command)
    IF (PRESENT(exitstat)) exitstat = e
    IF (e /= 0) THEN
      IF (PRESENT(cmdmsg)) cmdmsg = 'unknown'
    END IF
    IF (PRESENT(cmdstat)) cmdstat = 0
  END SUBROUTINE execute_command_line
#endif

END MODULE radar_utilities
