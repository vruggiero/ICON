!
!+ Data module for some constants
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

MODULE data_constants

!-------------------------------------------------------------------------------
!
! Description:
!  Data module for some constants. Provides minimum functionality of
!  COSMO module 'data_constants' to serve the shared modules in
!  DACE.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_44        2015-09-30 Andreas Rhodin
!  New module
! V1_51        2017-02-24 Andreas Rhodin
!  changed comment lines
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!------------------------------------------------------------------------------
   use mo_kind, only: wp
   implicit none
   private

   public :: repsilon   ,&!
             rprecision ,&! precision in current floating point format
             pi         ,&! 3.1415...
             r_v        ,&! gas constant for water vapor
             r_earth    ,&! mean radius of the earth (m)
             rho_w      ,&! density of liquid water (kg/m^3)
             rho_ice    ,&! density of ice          (kg/m^3)
             K_w        ,&! dielectric constant for water
             K_ice        ! dielectric constant for ice

   real(wp), parameter :: repsilon   = 1.0e8_wp *        tiny (1.0_wp)
   real(wp), parameter :: rprecision = 10.0_wp ** (-precision (1.0_wp))
!  real(wp), parameter :: pi         =  4.0_wp  *        atan (1.0_wp)
   real(wp), parameter :: pi         = 3.141592653589793238_wp

   real(wp), parameter :: r_v     =   461.51_wp    ! gas constant for water vapor
   real(wp), parameter :: r_earth =  6371.229E3_wp ! mean radius of the earth (m)
   real(wp), parameter :: rho_w   =  1000.0_wp     ! density of liquid water (kg/m^3)
   real(wp), parameter :: rho_ice =   900.0_wp     ! density of ice          (kg/m^3)
   real(wp), parameter :: K_w     =     0.930_wp   ! dielectric constant for water
   real(wp), parameter :: K_ice   =     0.176_wp   ! dielectric constant for ice

end module data_constants

!==============================================================================
