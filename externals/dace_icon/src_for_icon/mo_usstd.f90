!
!+ Calculate US 1976 Standard Atmosphere (pressure from geopotential height)
!
MODULE mo_usstd
!
! Description:
!   Calculate US 1976 Standard Atmosphere (pressure from geopotential height).
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_22        2013-02-13 Andreas Rhodin
!  declare p_h_usstd as ELEMENTAL function
!
! Code Description:
! Language: Fortran.
! Software Standards:
!------------------------------------------------------------------------------

  !------------
  ! Moduls used
  !------------
  use mo_kind,      only: wp     ! working precision kind parameter
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: p_h_usstd ! derive pressure from geopotential height
  public :: p_r       ! ISA sea level reference pressure
  public :: t_r       ! ISA sea level reference temperature
  public :: h_p_usstd ! derive geopotential height from pressure
  public :: t_h_usstd ! derive temperature from geopotential height

  !-----------------
  ! Module constants
  !-----------------
  real(wp), parameter :: p_r = 101325._wp    ! Sea level pressure    (Pa)
  real(wp), parameter :: t_r = 288.15_wp     ! Sea level temperature (K)
  real(wp), parameter :: hc  = 34.1631947_wp ! Hydrostatic constant  (K/km)
  real(wp), parameter :: r_r = 1.225_wp      ! The sea level density (kg/m3)
                                             !   derived from quantities above
  !============================================================================
  !     L O C A L   C O N S T A N T S                                         |
  !============================================================================
  REAL(wp), PARAMETER :: GMR = 34.163195_wp  ! hydrostatic constant
  INTEGER,  PARAMETER :: NTAB=8              ! number of entries in the tables
  !============================================================================
  !     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
  !============================================================================
  REAL(wp),DIMENSION(NTAB),PARAMETER:: htab= &
       (/0.0_wp, 11.0_wp, 20.0_wp, 32.0_wp,  &
        47.0_wp, 51.0_wp, 71.0_wp, 84.852_wp/)
  REAL(wp),DIMENSION(NTAB),PARAMETER:: ttab= &
       (/288.15_wp, 216.65_wp, 216.65_wp, 228.65_wp, &
         270.65_wp, 270.65_wp, 214.65_wp, 186.946_wp/)
  REAL(wp),DIMENSION(NTAB),PARAMETER:: ptab= &
       (/1.0_wp,          2.233611E-1_wp,  5.403295E-2_wp,  8.5666784E-3_wp,&
         1.0945601E-3_wp, 6.6063531E-4_wp, 3.9046834E-5_wp, 3.68501E-6_wp/)
  REAL(wp),DIMENSION(NTAB),PARAMETER:: gtab= &
       (/-6.5_wp, 0.0_wp, 1.0_wp, 2.8_wp, 0.0_wp, -2.8_wp, -2.0_wp, 0.0_wp/)
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
  elemental function p_h_usstd (h) result (p)
  !-----------------------------------------
  ! derive pressure from geopotential height
  !-----------------------------------------
  real(wp), intent(in) :: h ! geopotential height (gpm)
  real(wp)             :: p ! pressure            (Pa)

    real(wp) :: sigma, delta, theta
    call Atmosphere ( h / 1000._wp, sigma, delta, theta)
    p = delta * p_r

  end function p_h_usstd
!------------------------------------------------------------------------------
! The 1976 Standard Atmosphere to 86 kilometers is defined. The data is taken
! from the official U.S. Government publication.
!
! The equations used are those adopted 15 October 1976 by the United States
! Committee on Extension to the Standard Atmosphere (COESA), representing 29
! U.S. scientific and engineering organizations. The values selected in 1976
! are slight modifications of those adopted in 1962. The equations and
! parameters used are documented in a book entitled U.S. Standard Atmosphere,
! 1976 published by the U.S. Government Printing Office, Washington, D.C.
!
! The Fundamental 7 layers of the Standard Atmosphere to 86 km
!
! h1 and h2 are geopotential altitude in kilometers of the lower and upper
! boundries of a layer. The gradient dT/dH is kelvins per kilometer.
!
!  h1(km)h2(km) dT/dh (K/km)
!  0     11     -6.5
!  11    20     0.0
!  20    32     1.0
!  32    47     2.8
!  47    51     0.0
!  51    71     -2.8
!  71    84.852 -2.0
!
! Note: 84.852 km geopotential=86 km geometric
!
! These data along with the sea level standard values of
! Sea level pressure = 101325 N/m2
! Sea level temperature = 288.15 K
! Hydrostatic constant = 34.1631947 kelvin/km
! define the atmosphere. The sea level density of 1.225 kg/m3 is derived from
! the fundamental quantities above.
!------------------------------------------------------------------------------
ELEMENTAL SUBROUTINE Atmosphere( h , sigma, delta, theta)
!
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
! AUTHOR - Ralph Carmichael, Public Domain Aeronautical Software
! NOTE - If alt > 86, the values returned will not be correct, but they will
!   not be too far removed from the correct values for density.
!   The reference document does not use the terms pressure and temperature
!   above 86 km.
!============================================================================
!     A R G U M E N T S                                                     |
!============================================================================
! REAL(wp),INTENT(IN)::  alt       ! geometric    altitude, km.
  REAL(wp),INTENT(IN)::  h         ! geopotential altitude, km.
  REAL(wp),INTENT(OUT):: sigma     ! density/sea-level standard density
  REAL(wp),INTENT(OUT):: delta     ! pressure/sea-level standard pressure
  REAL(wp),INTENT(OUT):: theta     ! temperature/sea-level standard temperature
!============================================================================
!     L O C A L   C O N S T A N T S                                         |
!============================================================================
! REAL(wp),PARAMETER:: REARTH = 6369.0_wp   ! radius of the Earth (km)
!============================================================================
!     L O C A L   V A R I A B L E S                                         |
!============================================================================
  INTEGER:: i,j,k            ! counters
! REAL(wp):: h               ! geopotential altitude (km)
  REAL(wp):: tgrad, tbase    ! temperature gradient and base temp of this layer
  REAL(wp):: tlocal          ! local temperature
  REAL(wp):: deltah          ! height above base of this layer
!----------------------------------------------------------------------------

! h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude

  i=1
  j=NTAB                                       ! setting up for binary search
  DO
    k=(i+j)/2                                  ! integer division
    IF (h < htab(k)) THEN
      j=k
    ELSE
      i=k
    END IF
    IF (j <= i+1) EXIT
  END DO

  tgrad=gtab(i)                                ! i will be in 1...NTAB-1
  tbase=ttab(i)
  deltah=h-htab(i)
  tlocal=tbase+tgrad*deltah
  theta=tlocal/ttab(1)                         ! temperature ratio

  IF (tgrad == 0.0_wp) THEN                    ! pressure ratio
    delta=ptab(i)*EXP(-GMR*deltah/tbase)
  ELSE
    delta=ptab(i)*(tbase/tlocal)**(GMR/tgrad)
  END IF

  sigma=delta/theta                            ! density ratio
  RETURN
END Subroutine Atmosphere
!------------------------------------------------------------------------------
  elemental function h_p_usstd (p) result (h)
    !-----------------------------------------
    ! derive geopotential height from pressure
    !-----------------------------------------
    real(wp), intent(in) :: p ! pressure            (Pa)
    real(wp)             :: h ! geopotential height (gpm)
    !==========================================================================
    !     L O C A L   V A R I A B L E S                                       |
    !==========================================================================
    INTEGER  :: i,j,k           ! indices
    real(wp) :: delta           ! pressure/sea-level standard pressure
    REAL(wp) :: tgrad, tbase    ! temperature gradient, base temp of this layer
    REAL(wp) :: tlocal          ! local temperature
    REAL(wp) :: deltah          ! height above base of this layer
    !--------------------------------------------------------------------------
    delta = p / p_r

    i = 1
    j = NTAB                    ! setting up for binary search
    do
       k=(i+j)/2
       if (delta > ptab(k)) then
          j=k
       else
          i=k
       end if
       if (j <= i+1) exit
    end do

    tbase = ttab(i)             ! i will be in 1...NTAB-1
    tgrad = gtab(i)
    if (tgrad == 0.0_wp) then
       deltah = log (ptab(i) / delta) * tbase / GMR
    else
       tlocal = tbase * exp (log (ptab(i) / delta) * tgrad / GMR)
       deltah = (tlocal - tbase) / tgrad
    end if
    h = htab(i) + deltah
    h = h * 1000._wp            ! km -> m

  end function h_p_usstd
!------------------------------------------------------------------------------
  !--------------------------------------------
  ! derive temperature from geopotential height
  !--------------------------------------------
  elemental function t_h_usstd (h) result (t)
    real(wp), intent(in) :: h   ! geopotential height (gpm)
    real(wp)             :: t   ! temperature         (K)

    real(wp) :: sigma, delta, theta

    call Atmosphere (h / 1000._wp, sigma, delta, theta)
    t = theta * t_r
  end function t_h_usstd
!------------------------------------------------------------------------------
end module mo_usstd
