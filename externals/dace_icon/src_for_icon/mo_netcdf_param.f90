!
!+ Interface the NetCDF include file 'netcdf.inc'
!
! $Id$
!
MODULE mo_netcdf_param
!
! Description:
!   This module provides acess to the parameters defined in the
!   file 'netcdf.inc' of the NetCDF package
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
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!------------------------------------------------------------------------------

public

include 'netcdf.inc'

end module mo_netcdf_param
