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

MODULE radar_kind

!-------------------------------------------------------------------------------
!
! Description:
!  This module declares the kind parameters for integer and real variables.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  ! Kind parameters of radar variables in case of COSMO and ICON applications:
#ifdef __COSMO__
  ! ... also for __DACE__, see above!
  USE kind_parameters, ONLY : dp, wp, sp
#endif
#ifdef __ICON__
  USE mo_kind, ONLY : dp, wp, sp
#endif


  IMPLICIT NONE

  !==============================================================================

  ! INCLUDE statements

  !==============================================================================

  PUBLIC

  ! Kind parameters of radar variables in case of other applications:
#if !(defined __COSMO__ || defined __ICON__ || defined _DACE_)
  INTEGER, PARAMETER :: dp     = KIND(1.0d0)        ! double precision (8 byte reals)
  INTEGER, PARAMETER :: wp     = dp                 ! working precision
  INTEGER, PARAMETER :: sp     = KIND(1.0)
#endif
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND( 18) !< at least 8 byte integer

#ifdef SINGLEPRECISIONRADARFWO
!!$ This is not yet fully implemented!
  INTEGER, PARAMETER :: wpfwo = sp
!!$ This will crash compilation and prevent false usage of -DSINGLEPRECISIONRADARFWO:
  -DSINGLEPRECISIONRADARFWO is not yet implemented
#else
  INTEGER, PARAMETER :: wpfwo = dp
#endif


END MODULE radar_kind

