!
!+ Defines kind parameters
!
!==============================================================================
! Define kind parameters
!
MODULE kind_parameters
!
! Description:
!   Define kind parameters.
!   This moduleprovides a common interface for modules in COSMO and DACE.
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
!
!------------------------------------------------------------------------------
   use mo_kind
   implicit none
   private
   public :: sp, dp, wp     ! real    type kind parameters
   public :: i1, i2, i4, i8 ! integer type kind parameters

end module kind_parameters

!==============================================================================
