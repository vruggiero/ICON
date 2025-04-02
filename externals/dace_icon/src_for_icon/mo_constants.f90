!
!+ define physical and mathematical constants
!
!------------------------------------------------------------------------------
!> define physical and mathematical constants
!>
!> Define some physical and mathematical constants.
!>
!> \todo Consistence with parameters in module mo_physics must be checked
!>
MODULE mo_constants
!
! Description:
!   Define some physical and mathematical constants.
!   Consistence with parameters in module mo_physics must be checked !!!
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
! V1_7         2009/08/24 Andreas Rhodin
!  add directives for doxygen
! V1_8         2009/12/09 Harald Anlauf
!  Add "implicit none", fix pid5
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Detlef Pingel
!  constants for IASI
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM  2002  original code
!------------------------------------------------------------------------------

  use mo_kind, ONLY: dp
  implicit none
  !--------------------------------------
  ! Public entities used by other modules
  !--------------------------------------
  private
  public :: pi       ! 3.14159265358979323846
  public :: pid5     ! pi / 5.
  Public :: omcor    ! coriolis parameter
  public :: re       ! radius of the earth
  public :: rad_con2 ! 2nd radiation constant: h*c/k:
  !-----------------------
  ! Mathematical constants
  !-----------------------
  !
  !>  pi (64 bit IEEE)
  !
  real(dp), parameter :: pi       = 3.14159265358979323846_dp
  !
  !>  2 pi
  !
  real(dp), parameter :: pi2      = 2*pi
  !
  !>  pi / 2
  !
  real(dp), parameter :: pid2     = 0.5_dp*pi
  !
  !>  pi / 3
  !
  real(dp), parameter :: pid3     = pi/3.0_dp
  !
  !>  pi / 5
  !
  real(dp), parameter :: pid5     = pi/5.0_dp
  !
  !>  180 / pi
  !
  real(dp), parameter :: pid180i  = 180.0_dp/pi

  !-------------------
  ! physical constants
  !-------------------
  !
  !>  radius of the earth [m]
  !
  real(dp), parameter :: re       = 6371229.0_dp
  !
  !>  length of a day     [s]
  !
  real(dp), parameter :: day_len  = 86164.09054_dp
  !
  !>  coriolis parameter: 2 pi / day_len
  !
  real(dp), parameter :: omcor    = 2*pi/day_len
  !
  !>  Planck constant [J s]
  !
  real(dp),parameter :: h_planck = 6.62606876e-34_dp
  !
  !>  Boltzmann constant: k [J K^-1]
  !
  real(dp),parameter :: boltz_con = 1.3806503e-23_dp
  !
  !>  Speed of light in vacuum: c [m s^-1]
  !
  real(dp),parameter :: c_vacuum  = 299792458.0_dp
  !
  !>  2nd radiation constant: h*c/k: [m K]
  !
  real(dp),parameter :: rad_con2  = h_planck * c_vacuum / boltz_con

END MODULE mo_constants
