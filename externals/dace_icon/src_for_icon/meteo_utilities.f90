!
!+ Source module for meteorological utility routines
!
! 3dvar/letkf CVS :
!
!==============================================================================

MODULE  meteo_utilities

!==============================================================================
!
! Description:
!   This module provides service utilities for meteorological calculations.
!     - no routine uses other modules, except the declarations for the
!       KIND-type parameter; the data access is by parameter list only
!     - no routine allocates dynamic memory; work space needed is
!       provided via the parameter list
!     - no derived data types are used
!
!   Routines (module procedures) currently contained:
!
!     - calrho
!       Computes the air density for time-level nx
!
!     - calrho_densities
!       Computes the air density for time-level nx
!
!     - calrho_tp_pp
!       Computes the air density for time-level nx
!
!     - calps
!       Computes the surface pressure ps from the lowest full model level.
!
!     - cloud_diag
!       computes the cloud water and - cover including subscale clouds
!
!++++++++++++++++++++++++++++ to be removed +++++++++++++++++++++++++
!     - reference_atmosphere
!       Computes the reference atmosphere and vertical coordinate
!       parameters.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     - moist_split
!       Splits the generalized relative humidity into vapor, cloud water
!       and cloud ice, according to current temperature and pressure
!       (mixed phase cloud)
!
!     - satad
!       Corrects the temperature, the specific humidity and the cloud water
!       content for condensation/evaporation.
!
!     - psat_w
!       Saturation water vapour pressure
!
!     - qsat
!       Specific humidity at saturation pressure
!
!     - tgcom
!       Computes the ground temperature.
!
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.5        1998/06/29 Guenther Doms
!  Calculation of new array hhrl in routine 'reference_atmospere'
! 1.7        1998/07/16 Guenther Doms
!  New routines for calculation of reduces surface pressure (calpmsl),
!  vertical velocity p-dot (calomega), relative humidity (calrelhum)
!  and total precipitation (calprsum).
!  Change in the calling list for calps and removal of caldpst and calrssk.
! 1.9        1998/09/16 Guenther Doms
!  New routine 'caltopdc' to calculate the top index of dry convection.
! 1.10       1998/09/29 Ulrich Schaettler
!  Introduced logical switch llm to routine reference_atmospere.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adaptations to use this module also in GME2LM
! 1.30       1999/06/24 Matthias Raschendofer
!  New routine 'synop_diag' to give an alternative calculation of 2m- and 10m-values.
!  New routine 'cloud_diag' to make the the diagnosis of cloud water and its coverage.
! 1.32       1999/08/24 Guenther Doms
!  All cal...-routines for postprocessing moved to 'pp_utilities.f90.
!  Routines 'low_winds_tem' and 'synop_diag' moved to 'near_surface.f90'.
!  Routine 'cloud_diag' adaped to software standards.
! 1.34       1999/12/10 Matthias Raschendorfer
!  Changed critical value for dq/sig in routine cloud_diag
! 2.2        2000/08/18 Matthias Raschendorfer
!  Generalization of the statistical cloud scheme in routine cloud_diag.
! 2.8        2001/07/06 Ulrich Schaettler
!  Introduced optimizations for vectorization (same as in GME2LM)
! 2.14       2002/02/15 Ulrich Schaettler
!  Modifications to SR reference_atmosphere to use the SLEVE vertical
!  coordinate
! 2.18       2002/07/16 Reinhold Schrodin
!  Eliminated variable rhde; use cf_snow instead; changed calculation of tgcom
! 3.7        2004/02/18 Ulrich Schaettler
!  Editorial changes
! 3.13       2004/12/03 Ulrich Schaettler
!  Included ak,bk in interface of SR reference_atmosphere;
!  Changes in qsat according to new version of INT2LM
!  Changed names of local statement-functions in cloud_diag because of
!    naming conflicts
! 3.16       2005/07/22 Jochen Foerstner, Ulrich Schaettler
!  Included new SR for computing calrho_densities
!  Modification in SR satad, to "save" the start values for the temperature
! 3.21       2006/12/04 Ulrich Schaettler, Jochen Foerstner
!  Put statement functions which are not used to a comment
!  Changed interface from subroutine cloud_diag
!  New subroutine calrho_tp_pp: computes air density
!  Change in cloud_diag: ql is not subject to MIN anymore
!  Change in reference_atmosphere: calculate indices klv950, klv700 hpa
! V3_23        2007/03/30 Matthias Raschendorfer
!  Renaming the used statement functions in SUBROUTINE cloud_diag into the
!  systematic names used in 'stat_funct.incf', in order to substitute the
!  explicite definitions by that include file later.
! V4_1         2007/12/04 Ulrich Schaettler
!  Editorial Changes
! V4_5         2008/09/10 Guenther Zaengl
!  Additional subroutine for new reference atmosphere,
!  Corrected calculation of base-state pressure at full model levels for old
!  reference atmosphere
! V4_6         2008/09/24 Ulrich Schaettler
!  Removing the corrected calculation of base-state pressure again
!   (must be tested further)
! V4_8         2009/02/16 Guenther Zaengl
!  Finish implementation of consistent base-state calculation for new reference
!  atmosphere; old calculation is retained for standard reference atmosphere
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann (NEC)
!  Determine MIN of a variable in SR satad, to avoid creation of index list
! V4_11        2009/11/30 Guenther Zaengl
!  New SR reference_atmosphere_BVconst with assuming a constant Brunt-Vaisala frequency
! V4_12        2010/05/11 Ulrich Schaettler
!  Changes some arguments of subroutines reference_atmosphere* to optional
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  Alignment with INT2LM
! V4_20        2011/08/31 Ulrich Schaettler
!  Replaced arguments p0, pp of subroutine cloud_diag with full pressure p_tot
! V4_24        2012/06/22 Michael Baldauf
!  Changes in subroutines for reference_atmospheres:
!    introduced new argument t0hl and its computation
!    removed computations of k-indices for special pressure levels
!  Introduced new subroutine k_index_of_pressure_levels for this task
! V4_25        2012/09/28 Ulrich Blahak, Ulrich Schaettler
!  Bugfix cloud_diag() to avoid negative values of ql: changed
!    computation of clwc(i,j,k) resp. ql by using a new formula
!    provided by M. Raschendorfer.
!  Editorial changes to adapt to INT2LM Version 1.20 (US)
! V4_27        2013/03/19 Michael Baldauf
!  Used p0ref in computations of SR reference_atmosphere_BV and introduced this
!  variable in the argument list of that SR
! V4_28        2013/07/12 Ulrich Schaettler
!  Moved subroutines for reference atmospheres to module vgrid_refatm_utils
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_3         2015-10-09 Xavier Lapillonne
!  Introduced calrho_block, calps_block (at least as intermediate solution)
! V5_4         2016-03-10 Xavier Lapillonne, Oliver Fuhrer
!  Deactivate OpenACC statements for the moment being
!  Corrected comments
!  ============================================================================
! V1_4         2009/03/26 Andreas Rhodin
!  Borrow from lm_f90 V4_6
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  update to COSMO revision V4_11; remove unused variables
! V1_44        2015-09-30 Andreas Rhodin
!  update to COSMO 5.03-beta
! V1_51        2017-02-24 Andreas Rhodin
!  update shared modules from COSMO 5.04d
! V?_??
!  explicitly set public parameters used by other modules
!  (required by ICON Compiler directives when modules are ported to ICON)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
!------------------------------------------------------------------------------

USE data_parameters , ONLY :   &
  wp,        & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE
PUBLIC :: calps, calrho
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE calrho  ( t, pp, qv, qc, qrs, p0, rho, ie, je, ke,  &
                     r_d, rvd_m_o, lacc )

!------------------------------------------------------------------------------
! Description:
!   This routine computes the air density for present time level.
!
! Method:
!   Application of the gas-law, rho = p/(r_d*tv)
!
!------------------------------------------------------------------------------
!
! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke          ! dimensions of the fields

REAL (KIND=wp),     INTENT (IN)          ::    &
  t    (ie,je,ke),   & ! temperature
  pp   (ie,je,ke),   & ! pressure perturbation
  qv   (ie,je,ke),   & ! specific humidity
  qc   (ie,je,ke),   & ! specific cloud water content
  qrs  (ie,je,ke),   & ! specific rain and snow content
  p0   (ie,je,ke),   & ! base-state pressure
  r_d            ,   & ! gas constant of dry air
  rvd_m_o              ! r_v/r_d - 1

REAL (KIND=wp),     INTENT (OUT)       ::    &
  rho  (ie,je,ke)      ! air density

LOGICAL, INTENT(IN), OPTIONAL :: lacc !this should be set to true
                                      !if the function is called within
                                      !data region where in and out varible
                                      !are present
! Local variable:
LOGICAL :: lzacc
INTEGER(KIND=iintegers) :: i,j,k

!-------------------------------------------------------------------------------
! Begin subroutine calrho

  IF (PRESENT(lacc)) THEN
     lzacc=lacc
  ELSE
     lzacc=.FALSE.
  END IF

  !NOacc data present (t,pp,qv,qc,qrs,p0,rho) if(lzacc)
  !NOacc parallel if(lzacc)
  DO k=1,ke
    !NOacc loop gang
    DO j=1,je
      !NOacc loop vector
      DO i=1,ie
        rho(i,j,k) = ( p0(i,j,k) + pp(i,j,k) ) / ( r_d*t(i,j,k)*(1.0_wp +     &
                       rvd_m_o*qv(i,j,k) - qc(i,j,k) - qrs(i,j,k)) )
      ENDDO
    ENDDO
  ENDDO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE calrho

!==============================================================================
!==============================================================================
!block version
!------------------------------------------------------------------------------

SUBROUTINE calrho_block  ( t, pp, qv, qc, qrs, p0, rho, nproma, ke,  &
                           r_d, rvd_m_o, lacc )

!------------------------------------------------------------------------------
! Description:
!   This routine computes the air density for present time level.
!
! Method:
!   Application of the gas-law, rho = p/(r_d*tv)
!
!------------------------------------------------------------------------------
!
! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  nproma, ke          ! dimensions of the fields

REAL (KIND=wp),     INTENT (IN)          ::    &
  t    (nproma,ke),   & ! temperature
  pp   (nproma,ke),   & ! pressure perturbation
  qv   (nproma,ke),   & ! specific humidity
  qc   (nproma,ke),   & ! specific cloud water content
  qrs  (nproma,ke),   & ! specific rain and snow content
  p0   (nproma,ke),   & ! base-state pressure
  r_d            ,   & ! gas constant of dry air
  rvd_m_o              ! r_v/r_d - 1

REAL (KIND=wp),     INTENT (OUT)       ::    &
  rho  (nproma,ke)      ! air density

LOGICAL, INTENT(IN), OPTIONAL :: lacc !this should be set to true
                                      !if the function is called within
                                      !data region where in and out varible
                                      !are present
! Local variable:
LOGICAL :: lzacc
INTEGER(KIND=iintegers) :: ip,k

!-------------------------------------------------------------------------------
! Begin subroutine calrho_block

  IF (PRESENT(lacc)) THEN
     lzacc=lacc
  ELSE
     lzacc=.FALSE.
  END IF

  !NOacc data present (t,pp,qv,qc,qrs,p0,rho) if(lzacc)
  !NOacc parallel if(lzacc)
  DO k=1,ke
    !NOacc loop gang vector
    DO ip=1,nproma
      rho(ip,k) = ( p0(ip,k) + pp(ip,k) ) / ( r_d*t(ip,k)*(1.0_wp +     &
           rvd_m_o*qv(ip,k) - qc(ip,k) - qrs(ip,k)) )
    ENDDO
  ENDDO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE calrho_block

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE calrho_densities( t, pp, rho_v, rho_c, rho_rs, p0, rho,        &
                             ie, je, ke, r_d, rvd_m_o )

!------------------------------------------------------------------------------
! Description:
!   This routine computes the air density for present time level.
!
! Method:
!   Application of the gas-law, but input quantitities are densities
!
!------------------------------------------------------------------------------
!
! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke            ! dimensions of the fields

REAL (KIND=wp),     INTENT (IN)          ::    &
  t     (ie,je,ke),   & ! temperature
  pp    (ie,je,ke),   & ! pressure perturbation
  rho_v (ie,je,ke),   & ! density humidity
  rho_c (ie,je,ke),   & ! density cloud water content
  rho_rs(ie,je,ke),   & ! density rain and snow content
  p0    (ie,je,ke),   & ! base-state pressure
  r_d             ,   & ! gas constant of dry air
  rvd_m_o               ! r_v/r_d - 1

REAL (KIND=wp),     INTENT (OUT)       ::    &
  rho  (ie,je,ke)       ! air density

!------------------------------------------------------------------------------
! Begin subroutine calrho_densities

  rho = ( p0 + pp ) / ( r_d * t ) - ( rvd_m_o*rho_v - rho_c - rho_rs )

END SUBROUTINE calrho_densities

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE calrho_tp_pp( tp, pp, qv, qc, qrs, t0, p0, rho,    &
  ie, je, ke, r_d, rvd_m_o )

!------------------------------------------------------------------------------
! Description:
!   This routine computes the air density for present time level.
!
! Method:
!   Application of the gas-law, rho = p/(r_d*tv)
!
!------------------------------------------------------------------------------
!
! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke           ! dimensions of the fields

REAL (KIND=wp),     INTENT (IN)          ::    &
  tp   (ie,je,ke),   & ! temperature perturbation
  pp   (ie,je,ke),   & ! pressure perturbation
  qv   (ie,je,ke),   & ! specific humidity
  qc   (ie,je,ke),   & ! specific cloud water content
  qrs  (ie,je,ke),   & ! specific rain and snow content
  t0   (ie,je,ke),   & ! base-state temperature
  p0   (ie,je,ke),   & ! base-state pressure
  r_d            ,   & ! gas constant of dry air
  rvd_m_o              ! r_v/r_d - 1

REAL (KIND=wp),     INTENT (OUT)       ::    &
  rho  (ie,je,ke)      ! air density

!------------------------------------------------------------------------------
! Begin subroutine calrho_tp_pp

  rho = ( p0 + pp ) / ( r_d*( t0 + tp ) * (1.0_wp + rvd_m_o*qv - qc - qrs) )

END SUBROUTINE calrho_tp_pp

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE calps ( ps, pp, t, qv, qc, qrs, rho0, p0, dp0,          &
                   ie, je, rvd_m_o, r_d, istart, iend, jstart, jend )

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the surface pressure ps by extrapolating the
!   nonhydrostatic pressure of the lowest full model level in a
!   hydrostatic manner.
!
!   The fields are passed twodimensional, that means the calling procedure has
!   to choose the proper level and time level.
!
! Method:
!   Formula to be described in the documentation
!
!-------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je,                     & ! dimensions of the fields
  istart, iend, jstart, jend    ! start and end-indices of the computation

REAL (KIND=wp),     INTENT (IN)          ::    &
  pp   (ie,je)    , & ! perturbation pressure
  t    (ie,je)    , & ! temperature
  qv   (ie,je)    , & ! specific water vapor content
  qc   (ie,je)    , & ! specific cloud water content
  qrs  (ie,je)    , & ! precipitation water (water loading)
  p0   (ie,je)    , & ! reference pressure at full levels
  dp0  (ie,je)    , & ! full level pressure thickness
  rho0 (ie,je)        ! reference density at full levels

REAL (KIND=wp),     INTENT (INOUT)       ::    &
  ps   (ie,je)        ! surface pressure

REAL (KIND=wp),     INTENT (IN)          ::    &
  rvd_m_o,          & ! r_v/r_d - 1
  r_d                 ! gas constant for dry air

!-------------------------------------------------------------------------------

! Begin subroutine calps

  ps (istart:iend,jstart:jend) =                                              &
       ( p0(istart:iend,jstart:jend) + pp(istart:iend,jstart:jend) )          &
   *EXP (  0.5_wp*dp0(istart:iend,jstart:jend)  /                             &
         ( t(istart:iend,jstart:jend)                                         &
           * (1.0_wp + rvd_m_o*qv(istart:iend,jstart:jend)                    &
              - qc(istart:iend,jstart:jend) - qrs(istart:iend,jstart:jend))   &
           * r_d * rho0(istart:iend,jstart:jend)) )

END SUBROUTINE calps

!==============================================================================
!==============================================================================

SUBROUTINE calps_block ( ps, pp, t, qv, qc, qrs, rho0, p0, dp0,          &
                   nproma, rvd_m_o, r_d, ipstart, ipend, lmask, lacc )
!DIR$ INLINENEVER calps_block
!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the surface pressure ps by extrapolating the
!   nonhydrostatic pressure of the lowest full model level in a
!   hydrostatic manner.
!
!   The fields are passed twodimensional, that means the calling procedure has
!   to choose the proper level and time level.
!
! Method:
!   Formula to be described in the documentation
!
!-------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  nproma,                     & ! dimensions of the fields
  ipstart, ipend                ! start and end-indices of the computation

REAL (KIND=wp),     INTENT (IN)          ::    &
  pp   (nproma)    , & ! perturbation pressure
  t    (nproma)    , & ! temperature
  qv   (nproma)    , & ! specific water vapor content
  qc   (nproma)    , & ! specific cloud water content
  qrs  (nproma)    , & ! precipitation water (water loading)
  p0   (nproma)    , & ! reference pressure at full levels
  dp0  (nproma)    , & ! full level pressure thickness
  rho0 (nproma)        ! reference density at full levels

REAL (KIND=wp),     INTENT (INOUT)       ::    &
  ps   (nproma)        ! surface pressure

REAL (KIND=wp),     INTENT (IN)          ::    &
  rvd_m_o,          & ! r_v/r_d - 1
  r_d                 ! gas constant for dry air

LOGICAL, INTENT(IN), OPTIONAL            ::    &
  lmask(nproma)

LOGICAL, INTENT(IN), OPTIONAL :: lacc !this should be set to true
                                      !if the function is called within
                                      !data region where in and out varible
                                      !are present

! Local variable:
INTEGER (KIND=iintegers) :: ip
LOGICAL :: lzmask(nproma)
LOGICAL :: lzacc

!-------------------------------------------------------------------------------

! Begin subroutine calps_block
  IF (PRESENT(lacc)) THEN
     lzacc=lacc
  ELSE
     lzacc=.FALSE.
  END IF

  !NOacc data if(lzacc)   &
  !NOacc present (ps,p0,pp,dp0,t,qv,qc,qrs,rho0) &
  !NOacc create  (lzmask) !XL_TODO make lzmask global


  IF (PRESENT(lmask)) THEN
    !NOacc kernels present(lmask)
    lzmask(:) = lmask(:)
    !NOacc end kernels
  ELSE
    !NOacc kernels
    lzmask(:) = .TRUE.
    !NOacc end kernels
  ENDIF


  !NOacc parallel if(lzacc)
  !NOacc loop gang vector
  DO ip = ipstart, ipend
      IF (lzmask(ip)) THEN
        ps (ip) = (p0(ip) + pp(ip))* EXP(0.5_wp*dp0(ip)/(t(ip)           &
                 * (1.0_wp + rvd_m_o*qv(ip) - qc(ip) - qrs(ip))            &
                 * r_d * rho0(ip)))
      ENDIF
  ENDDO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE calps_block

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE cloud_diag ( clc, clwc,                         &
                        iis, iie, ijs, ije, iks, ike,      &
                        ids, ide, jds, jde, kds, kde,      &
                        ie , je , ke , ke1,                &
                        rdv, o_m_rdv, rvd_m_o, lhocp, t0,  &
                        b1, b2w, b3, b4w, b234w, b2i, b4i, &
                        uc1, uc2, ucl, clc_diag, q_crit,   &
                        t, qv, qc, p_tot,  rcld, ps,       &
                        itype_wcld )

!------------------------------------------------------------------------------
!
! Description:
!
!     This routine calculates the area fraction of a grid box covered
!     by stratiform (non-convective) clouds.
!     If subgrid-scale condensation is required, an additional
!     saturation adjustment is done.
!
! Method:
!
!     itype_wcld = 1 :
!     The fractional cloud cover clc is determined empirically from
!     relative humidity. Also, an in-cloud water content of sugrid-scale
!     clouds is determined as a fraction of the saturation specific
!     humidity. Both quantities define the grid-volume mean cloud water
!     content.
!     itype_wcld=2:
!     A Gaussion distribution is assumed for the saturation deficit
!     dq = qt - qs where qt = qv + ql is the total water content and
!     qs is the saturation specific humidity. Using the standard deviation
!     rcld of this distribution (on input) and the conservative grid-scale
!     quantities qt and tl (liquid water temperature), a corrected liquid
!     water content is determined which contains alse the contributions from
!     subgrid-scale clouds. A corresponding cloudiness is also calculated.
!
!------------------------------------------------------------------------------

! Subroutine arguments
!----------------------

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers), INTENT (IN) :: &  ! dimensions and run indices
  ie, je, ke, ke1, itype_wcld,   & !
  iis, ijs, iks, iie, ije, ike,  & !
  ids, jds, kds, ide, jde, kde

REAL (KIND=wp),     INTENT (IN)       :: &  ! Some physical constants
  rdv, o_m_rdv, rvd_m_o,  & ! Rd/Rv, 1-Rd/Rv, Rv/Rd-1
  lhocp, t0,              & ! Lev/Cpd, freezing temperature
  b1, b3, b2w, b4w,       & ! Constants to calculate saturation vapour pressure
  b2i, b4i, b234w,        & !  "
  uc1,uc2,ucl,            & ! empirical constants to calculate cloud cover
  clc_diag,               & ! cloud cover at saturation
  q_crit                    ! critical value for normalized over-saturation

! Array arguments with intent(in):

REAL (KIND=wp),     INTENT (IN)        :: & !
  t    (ie,je,ke ),    &    ! temperature (main levels)
  qv   (ie,je,ke ),    &    ! water vapour (")
  qc   (ie,je,ke ),    &    ! cloud water  (")
  p_tot(ie,je,ke ),    &    ! full pressure (")
  rcld (ie,je,ke1),    &    ! standard deviation of saturation deficit
  ps   (ie,je)              ! surface pressure

! Array arguments with intent(out):

REAL (KIND=wp),     INTENT (OUT)        :: &
  clc (ids:ide,jds:jde,kds:kde),  & ! stratiform subgrid-scale cloud cover
  clwc(ids:ide,jds:jde,kds:kde)     ! liquid water content of ""

! Local variables and constants
! -----------------------------

INTEGER (KIND=iintegers)  :: &
  i,j,k                     ! loop indices

REAL (KIND=wp), PARAMETER :: &
  zsig_max = 1.0E-3_wp,      & ! max. standard deviation of saturation deficit
  zclwfak  = 0.005_wp,       & ! fraction of saturation specific humidity
  zuc      = 0.95_wp           ! constant for critical relative humidity

REAL (KIND=wp)            :: &
  temp, pres, ql, qt, qs,    & !
  tl, dq, gam, q, sig, uc,   & !
  zsigma, zclc1, zq_max        !

REAL (KIND=wp)            :: &
  zpsat_w, zqvap, zdqsdt,    & !statement functions and
  zpvap, zqsat, ztemp, zpres   !their formal arguments

!------------ End of header ---------------------------------------------------

! Definition of statement functions:

! saturation vapour pressure over water (zpsat_w) and over ice (zpsat_i):
  zpsat_w(ztemp) = b1 * EXP( b2w*(ztemp-b3)/(ztemp-b4w) )
! zpsat_i(ztemp) = b1 * EXP( b2i*(ztemp-b3)/(ztemp-b4i) )

! specific humidity:
  zqvap(zpvap,zpres) = rdv * zpvap / ( zpres - o_m_rdv*zpvap )

! Derivation of zqsat with respect to temperature:
  zdqsdt(ztemp,zqsat) = b234w * ( 1.0_wp + rvd_m_o*zqsat ) * zqsat &
                             / (ztemp-b4w)**2

! Begin Subroutine cloud_diag
! ---------------------------

  zq_max   = q_crit*(1.0_wp/clc_diag - 1.0_wp)

  DO k = iks, ike
    DO j = ijs, ije
      DO i = iis, iie

        ql   = qc(i,j,k)               ! cloud water content
        qt   = ql + qv(i,j,k)          ! total water content
        pres = p_tot(i,j,k)            ! pressure
        temp = t(i,j,k)                ! temperature
        tl   = temp - lhocp*ql         ! liquid water temperature
        qs   = zqvap(zpsat_w(tl),pres) ! saturation mixing ratio
        dq   = qt - qs                 ! saturation deficit
        gam  = 1.0_wp / ( 1.0_wp + lhocp*zdqsdt(tl,qs) )

        IF ( itype_wcld == 1 ) THEN

        ! Calculation of cloud cover and cloud water content
        ! using an empirical relative humidity criterion

          zsigma = pres / ps(i,j)

          ! critical relative humidity
          uc     = zuc - uc1 * zsigma * ( 1.0_wp - zsigma )  &
                   * ( 1.0_wp + uc2*(zsigma-0.5_wp) )

          ! cloud cover
          clc(i,j,k) = MAX( 0.0_wp,  &
                       MIN( 1.0_wp, clc_diag * ((qt/qs-uc)/(ucl-uc))) )**2

          ! in-cloud water content
          ql = qs * zclwfak

          ! grid-volume water content
          IF ( dq > 0.0_wp ) THEN
            zclc1 = clc_diag * ( (1.0_wp-uc)/(ucl-uc) )**2
            ql    = ql + (gam*dq-ql)*(clc(i,j,k)-zclc1)/(1.0_wp-zclc1)
          END IF
          ql = clc(i,j,k) * ql

        ELSEIF ( itype_wcld == 2 ) THEN

        ! Statistical calculation of cloud cover and cloud water content
        ! using the standard deviation of the saturation deficit

          sig = MIN ( zsig_max, rcld(i,j,k) )

          ! in case of sig=0, the method is similar to grid-scale
          ! saturation adjustment. Otherwise, a fractional cloud cover
          ! is diagnosed.
          IF ( sig <= 0.0_wp ) THEN
            clc(i,j,k) = ABS ( (SIGN(1.0_wp,dq)+1.0_wp)*0.5_wp )
            ql         = clc(i,j,k) * gam * dq
          ELSE
            q          = dq / sig
            clc(i,j,k) = MIN ( 1.0_wp, MAX ( 0.0_wp, &
                                        clc_diag * (1.0_wp+q/q_crit) ) )
            IF ( q <= - q_crit ) THEN
              ql = 0.0_wp
            ELSEIF ( q >= zq_max ) THEN
              ql = gam * dq
            ELSE
              ! problems with ql < 0:
              ! ql = gam * sig * (q+q_crit) * (q+zq_max) / (2_wp*(q_crit+zq_max))
              ! New code by Matthias (changes results)
              ql = gam * sig * zq_max *( (q + q_crit)/(zq_max + q_crit) )**2
            ENDIF
          ENDIF

        ENDIF

        clwc(i,j,k) = ql

     ENDDO
   ENDDO
 ENDDO

END SUBROUTINE cloud_diag

!==============================================================================
!==============================================================================

SUBROUTINE moist_split (t, p, grh, qvmin, qcmin, qimin, pi,                 &
                        b1, b2_w, b2_i, b3, b4_w, b4_i, Rdv, O_m_rdv,       &
                        qv, qc, qi, ie, je)

!-------------------------------------------------------------------------------
!
! Description:
!   This routine splits the generalized relative humidity into vapor, liquid
!   water and ice, accordingly to the temperature and to the pressure.
!
! Method:
!   Tests based on the temperature, allowing representation of mixed phase
!   clouds, using saturation functions above water and above ice.
!
!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

 INTEGER (KIND=iintegers), INTENT (IN)    ::  &
    ie, je

 REAL    (KIND=wp),        INTENT (IN)    ::  &
    t(ie,je),                       & ! Current temperature
    p(ie,je),                       & ! Current pressure
    grh(ie,je),                     & ! Generalized relative humidity
    qvmin, qcmin, qimin,            & ! Lowest limit values
    b1, b2_w, b2_i, b3, b4_w, b4_i, & ! physical parameters
    Rdv, O_m_rdv, pi                  ! physical parameters

 REAL    (KIND=wp),        INTENT (OUT)   ::  &
    qv(ie,je),       & ! Specific water vapor content (specific humidity)
    qc(ie,je),       & ! Specific cloud water content
    qi(ie,je)          ! Specific cloud ice content

!------------------------------------------------------------------------------

! Local scalars:
! -------------
  REAL    (KIND=wp)                      ::  &
    T1, T2, Thom, Tzero,                     & ! Threshold temperatures
    alpha(ie,je), beta(ie,je),               & ! Cloud characteristic functions
    zaq, zaqi, zbq(ie,je), zbqi(ie,je),      & ! Saturation functions
    zot1mt2, zotzmt1, ztmp

  INTEGER (KIND=iintegers) :: &
    i,j

! Definition of statement functions:
! ---------------------------------

REAL    (KIND=wp)        ::   sf_psat_w, sf_psat_i,sf_qsat, x, y, z, zi, v, w, wi

sf_psat_w (x,y,z,v,w)  = y * EXP (z*(x-v)/(x-w))              ! Saturation pressure above water

sf_psat_i (x,y,zi,v,wi)= y * EXP(zi*(x-v)/(x-wi))             ! Saturation pressure above ice

sf_qsat    (x,y,z,v)   = z * x / MAX( (y-v*x), 1.0_wp)    ! Saturation specific humidity


!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations and definitions
!------------------------------------------------------------------------------

! 1.1. Threshold temperature values

Tzero   = 271.15_wp ! threshold for heterogeneous freezing of raindrops (or +/- zero Celsius)
T1      = 248.15_wp ! threshold for het.nucleation of cloud ice (water-saturation above)
Thom    = 236.15_wp ! homogeneous ice nucleation
T2      = 236.15_wp ! threshold for hom. freezing (ice-saturation under)

zot1mt2 = 1.0_wp/(T1-T2)
zotzmt1 = 1.0_wp/(Tzero-T1)

DO i = 1, ie
  DO j = 1, je

    ! 1.2. Saturation functions

    zaq  = sf_psat_w (t(i,j), b1, b2_w, b3, b4_w)
    zaqi = sf_psat_i (t(i,j), b1, b2_i, b3, b4_i)
    zbq(i,j)  = sf_qsat   (zaq,  p(i,j), Rdv, O_m_rdv)
    zbqi(i,j) = sf_qsat   (zaqi, p(i,j), Rdv, O_m_rdv)

    ! 1.3. Alpha function (= characteristic relative humidity, function of cloud temperature)
    !      Note: Alpha varies linearly between 1.0 and psati/psatl for temperature from T2 to T1

    ztmp = MIN(1.0_wp,MAX(0.0_wp,(T1-t(i,j))*zot1mt2))
    alpha(i,j) = 1.0_wp-(1.0_wp-zaqi/zaq)*ztmp

    ! 1.4. Beta function (= characteristic ratio between ice and (ice+water) in clouds)

    beta(i,j) = 1.0_wp
    ztmp = MIN(1.0_wp,MAX(0.0_wp,(t(i,j)-T1)*zotzmt1))
    beta(i,j) = 0.5_wp*(COS(pi*ztmp)+1.0_wp)
    IF ((t(i,j) < Tzero) .AND. alpha(i,j) < 1.0_wp) THEN
      beta(i,j) = 1.0_wp
    ENDIF

  END DO
END DO

!------------------------------------------------------------------------------
! Section 2: Split accordingly to temperature and pressure
!------------------------------------------------------------------------------

  !********************************************
  !                                           *
  ! Adaptation based on an original idea      *
  !         of  Ulrike Wacker                 *
  !                                           *
  !********************************************

DO i = 1, ie
  DO j = 1, je

    ! CASE A: Saturation impossible (grh <= alpha)

    qv(i,j) = grh(i,j) * zbq(i,j)
    qc(i,j) = 0.0_wp
    qi(i,j) = 0.0_wp

    ! CASE B: Saturation possible (grh > alpha)

    ! Case B1: T >= 0 , thus no ice
    IF ((grh(i,j) > alpha(i,j)) .AND. (t(i,j) >= Tzero)) THEN
      qv(i,j) = MIN(1.0_wp, grh(i,j)) * zbq(i,j)
      qc(i,j) = MAX(0.0_wp, grh(i,j)-1.0_wp) * zbq(i,j)
      qi(i,j) = 0.0_wp
    END IF

    ! Case B2a: Thom < T < Tzero , thus mixed phase possible, saturation over water
    IF ( (grh(i,j) > alpha(i,j)) .AND. (t(i,j) < Tzero) .AND. (t(i,j) > Thom) .AND. (alpha(i,j) == 1.0_wp)) THEN
      qv(i,j) = MIN(1.0_wp, grh(i,j)) * zbq(i,j)
      qc(i,j) = MAX(0.0_wp, grh(i,j)-1.0_wp) * zbq(i,j) * (1.0_wp-beta(i,j))
      qi(i,j)=  MAX(0.0_wp, grh(i,j)-1.0_wp) * zbq(i,j) * beta(i,j)
    END IF

    ! Case B2b: Thom < T < Tzero , thus mixed phase possible, saturation over ice
    IF ( (grh(i,j) > alpha(i,j)) .AND. (t(i,j) < Tzero) .AND. (t(i,j) > Thom) .AND. (alpha(i,j) < 1.0_wp)) THEN
      qv(i,j) = MIN(1.0_wp, grh(i,j)) * zbqi(i,j)
      qc(i,j) = 0.0_wp
      qi(i,j)=  MAX(0.0_wp, grh(i,j)-1.0_wp) * zbqi(i,j)
    ENDIF

    ! Case B3: T < Thom , thus no water
    IF ( (grh(i,j) > alpha(i,j)) .AND. (t(i,j) <= Thom)) THEN
      qv(i,j) = MIN(1.0_wp, grh(i,j)) * zbqi(i,j)
      qc(i,j) = 0.0_wp
      qi(i,j)=  MAX(0.0_wp, grh(i,j)-1.0_wp) * zbqi(i,j)
    ENDIF

  END DO
END DO

!------------------------------------------------------------------------------
! Section 3: Final values check and lower limit restriction (security)
!------------------------------------------------------------------------------

 qv(:,:) = MAX(qv(:,:),qvmin)
 qc(:,:) = MAX(qc(:,:),0.0_wp)
 qi(:,:) = MAX(qi(:,:),0.0_wp)

!------------ End of the subroutine -------------------------------------------

END SUBROUTINE moist_split

!==============================================================================
!==============================================================================

SUBROUTINE SATAD ( kitera, te, qve, qce, tstart, phfe,                        &
                   zdqd  , zqdwe, zh   , ztg0  , ztgn, zdqdt0, zgqd0, zphe ,  &
                   b1, b2w, b3, b4w, b234w, rdrd, emrdrd, rddrm1, lh_v, cpdr, &
                   cp_d, idim, jdim, ilo, iup, jlo, jup )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine corrects the temperature (te), the specific humidity (qve)
!   and the cloud water content (qce) for condensation/evaporation.
!
! Method:
!   Saturation adjustment, i.e. reversible condensation/evaporation at
!   constant pressure by assuming chemical equilibrium of water and vapor
!
!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers), INTENT (IN)    ::  &
    kitera,              & !  Number of iterations in the numerical scheme
    idim, jdim,          & !  Dimension of I/O-fields
    ilo, iup, jlo, jup     !  start- and end-indices for the computation

  REAL    (KIND=wp),        INTENT (IN)    ::  &
    tstart  (idim,jdim), & ! Start temperature for iteration
    phfe    (idim,jdim)  ! Pressure (input)

  REAL    (KIND=wp),        INTENT (INOUT) ::  &
    te      (idim,jdim), & ! Temperature on input/ouput
    qve     (idim,jdim), & ! Specific humidity on input/output
    qce     (idim,jdim), & ! Specific cloud water content on input/output
    zdqd    (idim,jdim), & !
    zqdwe   (idim,jdim), & !
    zh      (idim,jdim), & !
    ztg0    (idim,jdim), & !
    ztgn    (idim,jdim), & !
    zdqdt0  (idim,jdim), & !
    zgqd0   (idim,jdim), & !
    zphe    (idim,jdim)    !

  REAL    (KIND=wp),        INTENT (IN)    ::  &
    b1, b2w, b3, b4w, b234w, rdrd, emrdrd, rddrm1, lh_v, cpdr, cp_d

! Local parameters: None
! ----------------
! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i, j,                & !  Loop indices
    nzit,                & !  Loop for iterations
    nsat,                & !  Number of saturated gridpoints
    iwrk(idim*jdim),     & !  i-index of saturated gridpoints
    jwrk(idim*jdim),     & !  j-index of saturated gridpoints
    indx                   !  loop index

  REAL    (KIND=wp)        ::  &
    zgeu  ,              & !
    zgqdu ,              & !
    zgew  ,              & !
    zqwmin,              & ! Minimum cloud water content for adjustment
    fgew  ,              & ! Name of satement function
    fgqd  ,              & ! ...
    fdqdt ,              & ! ...
    zt    ,              & ! Dummy argument for statement functions
    zge   ,              & ! ...
    zp    ,              & ! ...
    zgqd                   ! ...

  REAL    (KIND=wp)        ::  &
    minzdqd

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine satad
!-------------------------------------------------------------------------------

! STATEMENT FUNCTIONS

fgew(zt)       = b1*EXP( b2w*(zt-b3)/(zt-b4w) )
fgqd(zge,zp)   = rdrd*zge/( zp - emrdrd*zge )
fdqdt(zt,zgqd) = b234w*( 1.0_wp + rddrm1*zgqd )*zgqd/( zt-b4w )**2

  zqwmin = 1.0E-20_wp

  nsat = 0

  minzdqd= 1.0_wp

  DO j = jlo , jup
    DO i = ilo , iup

      ! "save" the start values for the temperature
      ztg0 (i,j) = tstart(i,j)

      ! correction for negative values of qv and qc
      qve (i,j) = MAX( qve(i,j), 0.0_wp )
      qce (i,j) = MAX( qce(i,j), 0.0_wp )

      ! assume first subsaturation
      zqdwe(i,j)= qve(i,j) + qce(i,j)
      te (i,j)  = te(i,j) - lh_v*qce(i,j)*cpdr
      qve(i,j)  = zqdwe(i,j)
      qce(i,j)  = 0.0_wp
      zgeu      = fgew(te(i,j))
      zgqdu     = fgqd(zgeu,phfe(i,j))
      zdqd(i,j) = zgqdu - zqdwe(i,j)
      minzdqd   = MIN(minzdqd,zdqd(i,j))

    ENDDO
  ENDDO

!NEC_CB if zdqd>=0, then for sure no points are found
  IF ( minzdqd >= 0.0_wp ) RETURN

  DO j = jlo , jup
    DO i = ilo , iup

      IF (zdqd(i,j) < 0.0_wp ) THEN
        nsat       = nsat+1
        iwrk(nsat) = i
        jwrk(nsat) = j
      ENDIF

    ENDDO
  ENDDO

  IF (nsat == 0) RETURN

! Do saturation adjustments for saturated gridpoints
! --------------------------------------------------

!NEC$ ivdep
  DO indx = 1, nsat
     i = iwrk(indx)
     j = jwrk(indx)
     zh   (i,j) = cp_d*te(i,j) + lh_v*qve(i,j)
     zphe (i,j) = phfe(i,j)
     zgew       = fgew(ztg0(i,j))
     zgqd0(i,j) = fgqd(zgew,zphe(i,j))
  ENDDO

  IF ( kitera > 1 ) THEN
    DO  nzit  = 1 , kitera-1

!NEC$ ivdep
      DO indx = 1, nsat
        i = iwrk(indx)
        j = jwrk(indx)
        zdqdt0(i,j) = fdqdt(ztg0(i,j),zgqd0(i,j))
        ztg0(i,j)   = (zh(i,j) - lh_v*(zgqd0(i,j)-zdqdt0(i,j)*ztg0(i,j)))/ &
                      ( cp_d + lh_v*zdqdt0(i,j) )
        zgew        = fgew(ztg0(i,j))
        zgqd0(i,j)  = fgqd(zgew,zphe(i,j))
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------

!NEC$ ivdep
  DO indx = 1, nsat
      i = iwrk(indx)
      j = jwrk(indx)
      zdqdt0(i,j) = fdqdt(ztg0(i,j),zgqd0(i,j))
      ztgn(i,j)   = ( zh(i,j) - lh_v*(zgqd0(i,j)-zdqdt0(i,j)*ztg0(i,j)) ) / &
                    ( cp_d + lh_v*zdqdt0(i,j) )
      zgqd0(i,j)  = zgqd0(i,j) + zdqdt0(i,j)*( ztgn(i,j)-ztg0(i,j) )
  ENDDO

! Distribute the result on gridpoints
! -----------------------------------

!NEC$ ivdep
  DO indx = 1, nsat
      i = iwrk(indx)
      j = jwrk(indx)
      te (i,j) =  ztgn(i,j)
      qve(i,j) = zgqd0(i,j)
      qce(i,j) = MAX( zqdwe(i,j) - zgqd0(i,j), zqwmin )
  ENDDO

! End of the subroutine

END SUBROUTINE satad

!*******************************************************************************

FUNCTION psat_w(tx, b1, b2_w, b3, b4_w)

!-------------------------------------------------------------------------------
!
! Description:
!   Saturation water vapour pressure (with respect to water,
!   depending on the temperature "tx")
!
!-------------------------------------------------------------------------------

  REAL (KIND=wp),     INTENT(IN) :: tx, b1, b2_w, b3, b4_w

  REAL (KIND=wp)     :: psat_w

  psat_w     = B1*EXP(B2_w*(tx - B3)/(tx - B4_w))

END FUNCTION psat_w

!==============================================================================
!==============================================================================

FUNCTION qsat(psatx, px, Rdv, O_m_rdv)

!-------------------------------------------------------------------------------
!
! Description:
!   Specific humidity at saturation pressure (depending on the
!   saturation water vapour pressure "psatx" and the air pressure "px")
!
!-------------------------------------------------------------------------------

  REAL (KIND=wp),     INTENT(IN) :: psatx, px, Rdv, O_m_rdv

  REAL (KIND=wp)     :: qsat

!US changed by Davide, because of problems with ECMWF data high up in the sky
  qsat = Rdv*psatx/MAX((px-O_m_rdv*psatx),1.0_wp)

!US old version:
! qsat = Rdv*psatx/(px-O_m_rdv*psatx)

END FUNCTION qsat

!==============================================================================
!==============================================================================

SUBROUTINE tgcom (tg, ts, tb, ws, llp, ie, je, cf_snow,                 &
                  istart, iend, jstart, jend)

!-------------------------------------------------------------------------------
!
! Description:
!   Computation of the temperature tg at the boundary layer between the ground
!   and the atmosphere. Only 2-dimensional arrays can be passed to tgcom. It
!   must be called using the desired time level.
!
! Method:
!   For grid points above water and for grid points on land that are not
!   covered with snow:   tg = ground surface temperature tb
!   For snow covered land points, tg is a function of the temperature of the
!   the snow surface ts and the ground surface temperature tb:
!       tg = ts + exp( -rhde*ws ) * (tb-ts)
!   from Version 2.18 on replaced by
!       tg = ts + ( 1. - MIN(1.,ws/cf_snow)) * (tb -ts)
!
!-------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je,                     & ! dimensions of the fields
  istart, iend, jstart, jend    ! start and end-indices of the computation

REAL (KIND=wp),     INTENT (INOUT)       ::    &
  tg (ie,je)    ! temperature at the boundary between ground and atmosphere

REAL (KIND=wp),     INTENT (IN)          ::    &
  ts (ie,je), & ! temperature of the snow surface
  tb (ie,je), & ! temperature of the ground surface
  ws (ie,je)    ! water content of snow

LOGICAL,  INTENT (IN)                    ::    &
  llp (ie,je)   ! pattern of land- and sea-points

REAL (KIND=wp),     INTENT (IN)          ::    &
  cf_snow       ! factor for the computation

!-------------------------------------------------------------------------------

! Begin subroutine tgcom

  WHERE ( (llp(istart:iend,jstart:jend) .EQV. .TRUE.) .AND.                   &
                                       (ws(istart:iend,jstart:jend) > 0.0_wp) )
      tg(istart:iend,jstart:jend) = ts(istart:iend,jstart:jend) +             &
           (1.0_wp - MIN(1.0_wp,ws(istart:iend,jstart:jend)/cf_snow)) &
             * (tb(istart:iend,jstart:jend) - ts(istart:iend,jstart:jend))
  ELSEWHERE
      tg(istart:iend,jstart:jend) =   tb(istart:iend,jstart:jend)
  END WHERE

END SUBROUTINE tgcom

!==============================================================================

END MODULE meteo_utilities
