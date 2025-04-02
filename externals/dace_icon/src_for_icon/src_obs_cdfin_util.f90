!+ Source module for utilities for the observation processing, data assimilation
!-------------------------------------------------------------------------------

#if __NEC_VERSION__ / 100 == 305    /* Work around issues with nfort 3.5.x */
!NEC$ options "-fno-reorder-logical-expression"
#endif

MODULE src_obs_cdfin_util

!-------------------------------------------------------------------------------
! Description:
!   This module provides utilities used in the (pre-)processing of observation
!   reports read from NetCDF files.
!
!   Note: This module belongs to a group of COSMO-related modules for reading
!   ----  conventional data from ('cdfin') NetCDF observation input files.
!         It is shared between the COSMO model and the DACE program package !
!
! Method:
!   It contains the following module procedures with specific tasks:
!    - obs_assign_gridpt      : assign horiz'ly an obs report to a model grid pt
!    - get_global_surface     : gather surface fields (orography, land/sea mask)
!                                 from global domain
!    - obs_assign_sort_node   : sort reports according to targeting sub-domains
!    - obs_cdf_distrib_reports: distribute reports to appropriate sub-domains
!                                 (i.e. to local nodes)
!    - obs_gather_buffer      : gather buffer arrays (used for diagn. printing)
!    - std_atmosphere         : convert variables accord. to standard atmosphere
!    - std_atm_p2z_fast       : convert p to z    accord. to standard atmosphere
!    - obs_average_layers     : average high-res profile to interpolation layers
!    - obs_rm_duplicate_levels: remove duplicate levels
!    - obs_solar_zenith_angle : compute solar zenith angle (for an array of rep)
!    - obs_td2rh              : convert dewpoint T. to model-compatible rel. hum
!    - obs_qx2rh              : convert mixing ratio to model-compatible rel hum
!    - obs_rhw2rh             : make (observed) relat. humidity model compatible
!    - obs_rh2td              : convert rel. humidity into dewpoint temperature
!    - obs_find_level         : find levels for vertical interpolation
!    - check_wind_method      : add sat frequency digit to wind comput. method
!                                 (this routine must correspond to the one in
!                                  'mo_amv.f90' of DACE)
!
!   It also contains the following elemental functions with specific tasks:
!    - f_z2p           : get (model) pressure for a specified height
!    - f_p2dp          : get approx. model layer thickness at given pressure
!    - ncfeedobs_status: get report status in feedback file format
!    \ ireplace        : paste bit pattern from one word into another word
!    - ireplace1       : paste 1 bit from one word into another (integer) word
!    - ibit1           : return 1 bit at given position from an (integer) word
!    - rmod            : MOD function for positive REALS
!    - fpvsw           : get saturation water vapour pressure from temperature
!    - fpvsi           : saturation vapour pressure over ice from temperature
!    - ftd             : dewpoint temperature from water vapour and air pressure
!    - fpv2q           : specific humidity from water vapour and air pressure
!    - fq2pv           : water vapour pressure from specific humidity + pressure
!    - ileap           : decide whether a given year is a leap year
!
!   It uses from:
!    - parallel_utilities:    - distribute_values
!                             - gather_field
!                             - global_values
!                             - scatterv_values
!                             - gatherv_values
!    - environment:           - model_abort
!    - utilities:             - phi2phirot
!                             - rla2rlarot
!
!   Data modules used:
!    - kind_parameters
!    - mo_fdbk_tables
!    - data_obs_lib_cosmo
!
!   More modules from DACE used for subroutine 'list_of_neighbours'.
!
!   Note: This module contains MPI-based communication between PE's.
!
! Current Code Owner (for COSMO and for DACE):
!  DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V4_22        2012-01-31 Christoph Schraff
!  Initial version, extracted from other modules and adapted (e.g. modified
!  routine interfaces). New routines 'obs_gather_buffer', 'std_atmosphere',
!  and 'obs_solar_zenith_angle'.
! V4_28        2013/07/12 Christoph Schraff
!  - All direct calls of mpi routine replaced by calls of new routines in module
!    'parallel_utilities' (routines 'scatterv_values', 'gatherv_values').
!  - Function 'ncfeedobs_status' moved here from module 'src_obs_cdfout_feedobs'
!  - Statement functions replaced by new elemental functions 'ireplace1',
!    'ibit1', 'rmod', 'fpvsw', 'fpvsi', 'ftd', 'fpv2q', 'fq2pv', 'ileap'.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_3         2015-10-09 Christoph Schraff
!  Conversion routine 'obs_rh2td' introduced.
! V5_4f        2017-09-01 Christoph Schraff
!  Update from DACE V1_51: Introduction of subroutine 'list_of_neighbours' to
!  make grid point assignment compatible with DACE source and applicable to
!  ICON, GME, and IFS grids (in addition to COSMO).
! V5_5a        2018-06-22 Ulrich Schaettler
!  Use kind_parameters instead of data_parameters (US)
!  Removed all iintegers (US)
!  Removed irealgrib, which is not used here
! V5_6         2019-02-27 Christoph Schraff
!  Subroutine 'check_wind_method' added to include a satellite frequency digit
!  in the entry for the satellite derived wind computation method for AMV.
! V5_6b        2019-10-16 Christoph Schraff + Andreas Rhodin
!  - Subroutine 'obs_average_layers' added for averaging high-resolution
!    profiles to pre-specified interpolation layers for superobbing.
!  - Subroutine 'obs_rm_duplicate_levels' added for removing duplicate levels
!    (i.e. levels with same pressure).
!  - Updated from DACE: modifications in 'list_of_neighbours'
!     - minor bug fixes (if n_ngb == 0),
!     - correct computation of 'nolbg' for ICON-LAM.
!  - Set subroutines and functions used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
! @VERSION@  @DATE@     Christoph Schraff
!  - Slightly more accuracy at the interpolation to superobbing levels.
!  - Bug fix at loop range for superobbing of flags.
!  - Variables related to structured grid of COSMO such as 'r_dlat' removed
!    except in code line conditioned to ifdef COSMO.
!  - Some code clean-up.
! 2021-11-03 Christoph Schraff:
!    Additional safety added in obs_average_layers (if nlint <= 1; for missing
!    obs values)
! - 24.04.2023: Optional extrapolation beyond model orography in f_z2p.
! - 24.09.2024: 'std_atm_p2z_fast' added for efficient p2z with std atmosphere.
!
! ==============================================================================
! DACE history:
! ------------
!  1.51       2017-02-24 Andreas Rhodin : update shared modules from COSMO 5.04d
! ==============================================================================
! CAUTION: This module is used by both the DACE and COSMO model codes.       !!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------
    ! for 'ncfeedobs_status':

USE mo_fdbk_tables , ONLY :   &
    ST_ACTIVE       ,& ! status: used in the assimilation
    ST_MERGED       ,& ! status: not used, merged into multi-level report
    ST_PASSIVE      ,& ! status: not used, only monitored
    ST_REJECTED     ,& ! status: not used due to suspicious quality
    ST_PAS_REJ      ,& ! status: passive and rejected
    FL_OBSTYPE      ,& ! passive report type (at obs. location)
    FL_MERGE        ,& ! (report used for) merged report (only, e.g. TEMP ABCD
                       !   or single-level aircraft used for multi-level rep.)
    FL_THIN            ! thinning

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    gather_field    ,& ! gathers the parts of a total field from all subdomains
!   gather_values   ,& ! gathers a set of values from all nod. to 1 or all nodes
    global_values   ,& ! computes global values by operating on local arrays
    scatterv_values ,& ! scatters batches of variable length to other nodes
    gatherv_values  ,& ! gathers  batches of variable length from other nodes
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY: epsy , &
    ! constants for the horizontal rotated grid and related variables
    r_pollon       ,& ! longitude of the rotated north pole (in degrees, E>0)
    r_pollat       ,& ! latitude of the rotated north pole (in degrees, N>0)
    r_polgam       ,& ! angle between the north poles of the systems
    r_dlon         ,& ! grid point distance in zonal direction (in degrees)
    r_dlat         ,& ! grid point distance in meridional direction (in degrees)
    r_startlat_tot ,& ! transformed latitude of the lower left grid point
                      ! of the total domain (in degrees, N>0)
    r_startlon_tot ,& ! transformed longitude of the lower left grid point
                      ! of the total domain (in degrees, E>0)
    r_degrad          ! factor for transforming degree to rad

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
    phi2phirot     ,& ! converts phi from the real to the rotated system
    rla2rlarot        ! converts lambda from the real to the rotated system

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
PUBLIC :: obs_assign_gridpt,  get_global_surface,  obs_assign_sort_node,       &
          obs_cdf_distrib_reports,  obs_gather_buffer,  obs_average_layers,    &
          obs_rm_duplicate_levels, std_atmosphere,  std_atm_p2z_fast,          &
          obs_solar_zenith_angle,                                              &
          obs_td2rh,  obs_qx2rh,  obs_rhw2rh, obs_rh2td,  obs_find_level,      &
          check_wind_method,  f_z2p,  f_p2dp,  ncfeedobs_status

!-------------------------------------------------------------------------------

! Private Declarations
  REAL (KIND=wp), PARAMETER ,   PRIVATE  :: &
    c1     =  1.0_wp       ! standard real constant 1.0
!   epsy   =  1.E-8_wp     ! commonly used very small value > 0

  INTEGER,        PARAMETER ,   PRIVATE  :: &
    mgp    = 16            ! maximum size of grid-point index list

!-------------------------------------------------------------------------------

! Interface Blocks

INTERFACE obs_gather_buffer
  MODULE PROCEDURE   obs_gather_buffer_int,                                    &
                     obs_gather_buffer_real
END INTERFACE

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for assigning reports to grid points
!-------------------------------------------------------------------------------

SUBROUTINE obs_assign_gridpt ( roblat, roblon, isthgh, lupair, lsat_bilin      &
                             , lvertfit, ie_tot, je_tot, hsurf_tot, fland_tot  &
                             , doromx, fdoro, imdi, nboundlines , nolbc        &
                             , lseaobs , io_tot, jo_tot, zio_tot, zjo_tot      &
                             , zsurf, lsfcob )

!-------------------------------------------------------------------------------
! Description:
!   Assign observation location horizontally to a model mass grid point.
!   Remark:  nolbc: number of grid rows at the lateral boundary,
!                   where no observations are used.
!
! Method:
!   For sea observations, for reports with missing station height, or if the
!   search radius specified is negative, then the report is assigned to model
!   (sea) grid point, which is closest in the horizontal.
!   Otherwise, the report is assigned to the grid point with the smallest
!   weighted height difference to the observation station within the positive
!   search radius, unless there is a grid point within half the mesh width
!   with a height difference of less than 40 m.
!   The horizontal search radius is  (SQRT(2) * mesh width)  for TEMP, PILOT,
!   SYNOP and DRIBU. This choice will often yield small height differences,
!   which is desirable for the use of both surface-level data and of the
!   vertical structure of sounding data, while the horizontal offset will not
!   exceed the maximum offset when the search was limited to the 4 surrounding
!   grid points.
!   (Note: Code is extended to T2m and to search at latitude > 60 deg.) 04.11.97
!
!   Note: All the required input variables are contained in the parameter list,
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments
! --------------------

  REAL    (KIND=wp)        , INTENT (IN)    ::       &
    roblat, roblon   ,& ! geographical latitude and longitude
    ! the following input variables are identical at each call
    doromx (4)       ,& ! SYNOP obs with height differences betw model orography
                        !   and station height >= 'doromx' are set passive
    fdoro  (5)          ! scaling factor to vertical distances between model
                        !   orography and station height for (z-obs > z-mod),
                        !   when assigning a grid point to the observation

  INTEGER                  , INTENT (IN)    ::       &
    ! the following input variables are identical at each call
    imdi             ,& ! missing data indicator for station height
    nboundlines      ,& ! number of overlapping boundary lines of the subdomains
    ie_tot           ,& ! number of grid pts in zonal direct.  (in total domain)
    je_tot              ! number of grid pts in meridional dir (in total domain)

  REAL    (KIND=wp)        , INTENT (IN)  ::  &
    hsurf_tot (ie_tot,je_tot) ,& ! model orography (total model domain)
    fland_tot (ie_tot,je_tot)    ! land fraction   (total model domain)

  LOGICAL                  , INTENT (IN)  ::  &
    lupair           ,& ! .true, if upper-air report
    lvertfit         ,& ! .true, if then minimum vertical distance 'dz' between
                        !           station height and model orography has given
                        !           absolute priority even if an 'appropriate
                        !           (|dz| < 40m, zdo < 0.25) grid point nearby'
                        !           is available
    lsat_bilin          ! .true, if model values are interpolated bi-linearly
                        !           instead of assigning the obs to a grid point

  INTEGER                  , INTENT (INOUT) ::       &
    isthgh           ,& ! station height in metres (may be modified for Sat obs,
                        !   i.e. if bi-linear interpolation)
    nolbc               ! number of grid rows at lateral boundaries
                        !   where obs are neglected

  LOGICAL                  , INTENT (INOUT) ::       &
    lseaobs             ! true, if obs. released from ship or dribu (IN)
                        ! true, if obs. is treated as sea obs. (OUT)
! output
  INTEGER                  , INTENT (OUT)   ::       &
    io_tot  , jo_tot    ! grid point assigned to observation.
                        ! - if io_tot = jo_tot = 0 THEN obs. report is outside
                        !     horizontal limits
                        ! - if io_tot , jo_tot < 0 THEN obs. report is outside
                        !     vertical limits or land / sea type (for synops)

  REAL    (KIND=wp)        , INTENT (OUT)   ::       &
    zio_tot, zjo_tot ,& ! exact obs. location in grid point units
    zsurf               ! height of model grid pt. to which obs. is assigned

  LOGICAL                  , INTENT (OUT)   ::       &
    lsfcob              ! false, if no surface data are to be used from profile

! Local parameters:
! ----------------

  REAL    (KIND=wp)        , PARAMETER  :: &
    rtempmx =  1.414_wp ,& ! horizontal search radius to assign a TEMP /
                               !   PILOT to a grid pt. (if (rtempmx < 0), search
                               !   is limited to the 4 neighbouring grid pts.)
    rsypmx  =  1.414_wp ,& ! as 'rtempmx', but for surface-level reports
    c0      =  0.0_wp   ,& ! standard real constant 0.0
!   c1      =  1.0_wp   ,& ! standard real constant 1.0
    c2      =  2.0_wp   ,& ! standard real constant 2.0
    c05     =  0.5_wp      ! standard real constant 0.5
!   epsy    =  1.E-8_wp    ! commonly used very small value > 0

! Local variables
! ---------------

  INTEGER                  ::       &
!   iolb    , jolb   ,& ! observation location in total grid point units
    nlandgp          ,& ! number of nearby land grid points
!   iaseek  , ieseek ,& ! grid point search range
!   jaseek  , jeseek ,& !
    i,j,k,ko         ,& ! loop indices
    n_ngb            ,& ! number of direct neighbouring grid-points
    n_ret            ,& ! total number of grid-points returned
    idx (mgp, 4)     ,& ! index list ( i/j/diamond/pe , grid-point)
    ibnd             ,& ! distance to lateral boundary (grid-point units)
    nolbg  (mgp)        ! grid-point distance to lateral boundary

  REAL    (KIND=wp)        ::       &
    zfsea            ,& ! maximum sea fraction of candidate grid points
    rseekmx          ,& ! horizontal search radius for grid pts. to assign
                        ! a observation to a gridpoint
    vfac             ,& ! modification factor for height difference
    zdo              ,& ! horizontal distance between obs and grid point
    zdomx            ,& ! maximum distance
!   zrlats           ,& ! gridpoint latitude in degrees
!   zcrlat           ,& ! cosine of grid point latitude
    fisdmn           ,& ! minimum height difference
    fisd             ,& ! height difference
    fisdps           ,& ! weighted height difference for surface pressure obs.
    fisdrh           ,& ! weighted height difference for rel. humidity observat.
    fisdtt           ,& ! weighted height difference for temperature observat.
    fisduv           ,& ! weighted height difference for wind observations
    fland (mgp)      ,& ! land fraction of grid-point
    hsurf (mgp)      ,& ! elevation of grid-point
    dist  (mgp)      ,& ! horizontal distance between obs and grid point
    w     (mgp)         ! interpolation weights

  LOGICAL                  ::       &
    linarea             ! .true, if obs inside inner model area

! CHARACTER (LEN=20) yroutine ! name of subroutine

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE obs_assign_gridpt
!-------------------------------------------------------------------------------
! yroutine = 'obs_assign_gridpt'

! Provide list of grid-points in the vicinity of a given location.
! Generic routine to be used within the DACE and COSMO sources.
! Within DACE: Generic routine to be used within the COSMO/ICON/GME/IFS model.
! Note: subr. 'list_of_neighbours' works correctly only for (rseekmx < 2.0)!

   CALL list_of_neighbours ( roblat, roblon, ie_tot, je_tot                   &
                           , hsurf_tot, fland_tot                             &
                           , zio_tot, zjo_tot, n_ngb, n_ret                   &
                           , idx, ibnd, dist, fland, hsurf, nolbg, w )
!  =======================

!  Number of lateral grid rows where obs are neglegted
   nolbc   = MAX( nolbc , INT( 1+nboundlines ) )

!  iolb    = INT( zio_tot )
!  jolb    = INT( zjo_tot )
   io_tot  = 0
   jo_tot  = 0
   ko      = 0
   zsurf   = c0
   lsfcob  = .FALSE.
   linarea = (ibnd > nolbc)
   IF (.NOT. linarea)                                                     RETURN

   lsfcob  = .TRUE.

! Re-determine 'lseaobs' (i.e. whether land or sea observation),
! to account e.g. for ships on rivers or synops/temps on tiny islands
! (input value of 'lseaobs' depends on the obs type (sea obs if ship or buoy);
!  - if a surface-level sea obs has to be assigned to a land grid point
!    then 'lseaobs' is set false to prevent its further use
!  - if an upper-air sea obs (TEMPSHIP) is assigned to a land grid point with
!    < 5 % water fraction, 'lseaobs' is set false to prevent its further use
!  - if an upper-air sea obs (TEMPSHIP) is assigned to a land grid point with
!    > 5 % water fraction, only 'lsfcob' is set false (surface level not used)
!  - Note: currently the threshold of 5 % is reduced to 1 % to allow for using
!    TEMPSHIP upper-air obs in the harbour of Hamburg
!  - land obs types may be assigned to a sea grid point but will always be used
!    (e.g. Ekofisk) )

   nlandgp = SUM (NINT   (fland(1:n_ngb)))
   zfsea   = c1 - MINVAL (fland(1:n_ngb))

!  PRINT *,'zfsea ', zfsea, zio_tot, zjo_tot, nlandgp, lseaobs, lupair
   IF (      (lseaobs) .AND. (nlandgp == n_ngb)                                 &
       .AND. (.NOT. ((lupair) .AND. (zfsea > 0.01_wp))))  lseaobs = .FALSE.

   !   if one adjacent model grid point is a land point, then the observation
   !   is considered as land observation
   IF ((.NOT. lseaobs) .AND. (nlandgp == 0)) lseaobs = .TRUE.

   !   a satellite retrieval (over land) is assigned a station height by
   !   bi-linear interpolation of the orography to the obs location
   !   (shouldnt be used any more, only for 1dvar)
   IF (lsat_bilin)  isthgh = NINT( SUM (w(:n_ngb) * hsurf (:n_ngb)))

   IF ((isthgh == imdi) .OR. (lseaobs) .OR. (lsat_bilin)) THEN
     rseekmx = c0
     vfac    = c1
   ELSEIF (lupair) THEN
     rseekmx = rtempmx
     vfac    = c1
   ELSE
     rseekmx = rsypmx
     vfac    = fdoro(5)
   ENDIF
!CS: Note: subr. 'list_of_neighbours' works correctly only for (rseekmx < 2.0)!

   !   determine the index ranges containing all candidate grid points
   IF (rseekmx < epsy) THEN
     !   if (rseekmx <= 0) then take into account only adjacent grid pts
!    jaseek = MAX( jolb - 0 , nolbc + 1 )
!    jeseek = MIN( jolb + 1 , je_tot- nolbc )
!    iaseek = MAX( iolb - 0 , nolbc + 1 )
!    ieseek = MIN( iolb + 1 , ie_tot- nolbc )
     rseekmx = 1.42_wp
!  ELSE
!    jaseek = MAX( jolb - INT(rseekmx)     , nolbc + 1 )
!    jeseek = MIN( jolb + INT(rseekmx) + 1 , je_tot- nolbc )
!    zcrlat = COS( MAX( ABS( r_startlat_tot + (jaseek-1) *r_dlat )             &
!                     , ABS( r_startlat_tot + (jeseek-1) *r_dlat ) ) *r_degrad )
!    IF (zcrlat > epsy) THEN
!      iaseek = MAX( iolb - INT(rseekmx /zcrlat)     , nolbc + 1 )
!      ieseek = MIN( iolb + INT(rseekmx /zcrlat) + 1 , ie_tot- nolbc )
!    ELSE
!      iaseek = nolbc + 1
!      ieseek = ie_tot - nolbc
!    ENDIF
   ENDIF

   IF ((lseaobs) .AND. (nlandgp == n_ngb)) THEN
     !   if radiosonde ship has to be assigned to a land g.pt. with non-zero
     !   water fraction, it is used actively as a sea obs, but its surface level
     !   is set passive
     !   ('SUM' is used here only to convert array(1) into scalar)

     k      = SUM (MINLOC (fland (:n_ngb), nolbg(:n_ngb) >= nolbc))
     io_tot = idx (k,1)
     jo_tot = idx (k,2)
     ko     = k
     lsfcob = .FALSE.

   ELSEIF ((isthgh == imdi) .OR. (lseaobs)) THEN
     !   If there is no station height, assign the observation to the grid
     !   point which is closest in the horizontal
     !   the same applies for observations over the sea, if there are model
     !   sea grid points within the search radius.

     zdomx  =  c2 + rseekmx
     DO k = 1, n_ngb
       IF (nolbg(k) <  nolbc) CYCLE
       IF (dist(k) <= zdomx) THEN
         IF ((isthgh == imdi) .OR. (NINT(fland(k)) == 0)) THEN
           zdomx  = dist(k)
           io_tot = idx (k,1)
           jo_tot = idx (k,2)
           ko     = k
           lsfcob = .TRUE.
         ELSEIF (io_tot == 0) THEN
           !    assign a sea report to a land grid point if necessary
           io_tot = idx (k,1)
           jo_tot = idx (k,2)
           ko     = k
           lsfcob = .FALSE.
         ENDIF
       ENDIF
     ENDDO
     IF (isthgh == imdi) lsfcob = .FALSE.

   ELSE

     fisdmn  = 10000.0_wp
     !   ensure that search radius is large enough to find a land grid pt.
     rseekmx = MAX( ABS(rseekmx) , SQRT( 0.5_wp ) )
     IF (nlandgp <= 3) rseekmx = MAX( rseekmx , c1 )
     IF (nlandgp <= 2) rseekmx = MAX( rseekmx , SQRT(1.25_wp) )
     IF (nlandgp <= 1) rseekmx = MAX( rseekmx , SQRT(2.0_wp ) )
     DO k = 1, n_ret
       IF (nolbg(k) <  nolbc) CYCLE
       zdo = dist(k)
       IF ((zdo <= rseekmx) .AND. (NINT(fland(k)) == 1)) THEN
         fisd = isthgh - hsurf(k)
         IF (      (zdo < 0.25_wp) .AND. (ABS(fisd) < 40.0_wp)                 &
             .AND. (.NOT. lvertfit)) THEN
           io_tot = idx (k,1)
           jo_tot = idx (k,2)
           ko     = k
           fisdmn = 0._wp
           lsfcob = .TRUE.
         ENDIF
         !   height difference 'fisd' is weighted with factor 'fdoro(5)' if
         !   the station is below model grid point; 'fisd' is always positive;
         !   for an extrapolation of 100 m and a temperature error of 12 k,
         !   the resulting pressure error is about 0.5 hPa
         fisd = ((vfac-c1)/c2 + SIGN( (vfac+c1)/c2 , fisd )) * fisd
         IF (fisd < fisdmn) THEN
           io_tot = idx (k,1)
           jo_tot = idx (k,2)
           ko     = k
           fisdmn = fisd
           lsfcob = .TRUE.
         ENDIF
       ELSEIF ((zdo <= rseekmx) .AND. (io_tot == 0)) THEN
         !   assign a land report to a sea grid point if necessary
         io_tot = idx (k,1)
         jo_tot = idx (k,2)
         ko     = k
         lsfcob = .FALSE.
       ENDIF
     ENDDO

   ENDIF

   IF (ko  > 0)    zsurf = hsurf (ko)

   IF ((.NOT. lupair) .AND. (.NOT. lsfcob)) THEN
     io_tot = - io_tot
     jo_tot = - jo_tot
   ENDIF
   IF ((io_tot <= 0) .OR. (isthgh == imdi) .OR. (lsat_bilin))             RETURN
!  IF (     (((io_tot == 0) .OR. (jo_tot == 0)) .AND. (.NOT. linarea))         &
!      .OR. (isthgh == imdi))                                             RETURN

   fisd   = isthgh - zsurf
   fisdps = ((fdoro(2)-c1)/c2 + SIGN( (fdoro(2)+c1)/c2 , fisd) ) * fisd
   fisdrh = ((fdoro(4)-c1)/c2 + SIGN( (fdoro(4)+c1)/c2 , fisd) ) * fisd
   fisdtt = ((fdoro(3)-c1)/c2 + SIGN( (fdoro(3)+c1)/c2 , fisd) ) * fisd
   fisduv = ((fdoro(1)-c1)/c2 + SIGN( (fdoro(1)+c1)/c2 , fisd) ) * fisd
   IF (      (fisdps > doromx(2)) .AND. (fisdrh > doromx(4))                   &
       .AND. (fisdtt > doromx(3)) .AND. (fisduv > doromx(1))) THEN
     IF (.NOT. lupair) THEN
       io_tot = - io_tot
       jo_tot = - jo_tot
     ELSE
       lsfcob = .FALSE.
     ENDIF
   ENDIF

!-------------------------------------------------------------------------------

RETURN
END SUBROUTINE obs_assign_gridpt

!===============================================================================

SUBROUTINE list_of_neighbours ( roblat, roblon                       &
                              , ie_tot, je_tot, hsurf_tot, fland_tot &
                              , zio_tot, zjo_tot                     &
                              , n_ngb, n_ret                         &
                              , idx, ibnd, dist, fland, hsurf        &
                              , nolbg, w                             )
!-----------------------------------------------------------------------------
! Provide list of grid-points in the vicinity of a given location.
! Generic routine to be used within the DACE and COSMO sources.
! Within DACE: Generic routine to be used within the COSMO/ICON/GME/IFS model.
!-----------------------------------------------------------------------------
! Known bugs: Subroutine works only if halo is 1 grid row.
!-----------------------------------------------------------------------------

#ifndef __COSMO__
  !
  ! access to DACE modules
  !
  USE mo_grid_intpol,     ONLY: grid_indices ! get model indices for given location
  USE mo_t_icon,          ONLY: t_grid_cells ! Derived type for grid cells
  USE data_obs_lib_cosmo, ONLY: gen_grid     ! generic grid meta date (for ICON)
  USE mo_t_obs,           ONLY: unitvector   ! unit vector on sphere from lat/lon
  USE mo_physics,         ONLY: d2r,        &! factor to transform degree to rad
                                gacc         ! gravity acceleration
  USE mo_wmo_tables, ONLY: WMO6_LATLON,     &! regular grid types
                           WMO6_ROTLL,      &!
                           WMO6_GAUSSIAN,   &!
                           DWD6_ICOSAHEDRON,&! GME  grid type
                           DWD6_ICON         ! ICON grid type
#endif

  REAL    (KIND=wp)        , INTENT (IN)    ::       &
    roblat, roblon     ! geographical latitude and longitude

  INTEGER                  , INTENT (IN)  ::  &
    ie_tot           ,& ! number of grid pts in zonal direct.  (in total domain)
    je_tot              ! number of grid pts in meridional dir (in total domain)

  REAL    (KIND=wp)        , INTENT (IN)  ::  &
    hsurf_tot (ie_tot,je_tot) ,& ! model orography (total model domain)
    fland_tot (ie_tot,je_tot)    ! land fraction   (total model domain)

  REAL    (KIND=wp)        , INTENT (OUT)   ::       &
    zio_tot, zjo_tot,& ! exact obs. location in grid point units (if applicable)
    dist  (mgp)     ,& ! distance to roblat/roblon in grid-units
    fland (mgp)     ,& ! land fraction of grid-point
    hsurf (mgp)     ,& ! elevation of grid-point
    w     (mgp)        ! interpolation weights

  INTEGER                  , INTENT (OUT)  ::  &
    n_ngb           ,& ! number of direct neighbouring grid-points
    n_ret           ,& ! number of grid-points returned
    idx (mgp, 4)    ,& ! index list ( i/j/diamond , grid-point)
    ibnd            ,& ! distance to lateral boundary (grid-point units)
    nolbg  (mgp)       ! grid-point distance to lateral boundary

! Local variables
! ---------------

  REAL    (KIND=wp)        ::       &
    rlat   , rlon    ,& ! latitude and longitude in the rotated system
    zfio   , zfjo    ,& ! factors for bilinear interpolation of sat data
    zdo              ,& ! (squared) horizontal distance between obs and grid pt.
    zrlats           ,& ! gridpoint latitude in degrees
    zcrlat           ,& ! cosine of grid point latitude
    x(3)             ,& ! unit vector on the sphere
    d_rad               ! nominal grid spacing (rad)

  INTEGER                  ::       &
    iolb    , jolb   ,& ! observation location in total grid point units
    i, j, d, k, l    ,& ! indices
    i1,j1,ie,je,ii,jj   ! bounds, index

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!---------------------------------------------------------------
! Code to be used within COSMO sources, applicable to COSMO only
!---------------------------------------------------------------
!!!#define __COSMO__
#ifdef __COSMO__

    !  calculate latitude and longitude in the rotated system
    !  and exact obs. location in grid point units
    rlon = rla2rlarot (roblat, roblon, r_pollat, r_pollon, r_polgam)
    rlat = phi2phirot (roblat, roblon, r_pollat, r_pollon)

    !  observation position in grid point units (total area)
    IF (r_dlon /= 0._wp .and. r_dlat /= 0._wp) THEN
      !  regular grid (COSMO)
      zio_tot    = 1._wp + (rlon - r_startlon_tot) /r_dlon
      zjo_tot    = 1._wp + (rlat - r_startlat_tot) /r_dlat
    ELSE
      !  irregular grid (ICON, GME) grid
      zio_tot    = 0._wp
      zjo_tot    = 0._wp
    ENDIF

    iolb = INT( zio_tot )
    jolb = INT( zjo_tot )

    !  Fill direct neighbours into list (1:n_ngb)

    ibnd  = CEILING (MIN (zio_tot, zjo_tot, ie_tot+c1-zio_tot, je_tot+c1-zjo_tot) -epsy)
    k     = 0
    DO i = MAX (iolb, 1), MIN (iolb+1, ie_tot)
      DO j = MAX (jolb, 1), MIN (jolb+1, je_tot)
        k = k + 1
        idx (k,:) = [i,j,1,-1]
        fland (k) = fland_tot (i, j)
        hsurf (k) = hsurf_tot (i, j)
        nolbg (k) = MIN (i - 1, j - 1, ie_tot - i, je_tot - j)
        zrlats    = r_startlat_tot + (j-1) * r_dlat
        zcrlat    = COS ( zrlats * r_degrad )
        zdo       = ( (i-zio_tot)*zcrlat *r_dlon/r_dlat )**2 + (j-zjo_tot)**2
        dist (k)  = SQRT( zdo )
        zfio      = ABS (zio_tot - REAL(i, wp))
        zfjo      = ABS (zjo_tot - REAL(j, wp))
        w    (k)  = (c1-zfio)*(c1-zfjo)
      ENDDO
    ENDDO
    n_ngb = k

    !  Fill halo into list (n_ngb+1:n_ret)

    DO i = MAX (iolb-1, 1), MIN (iolb+2, ie_tot)
      DO j = MAX (jolb-1, 1), MIN (jolb+2, je_tot)
        IF ( (i==iolb .OR. i==iolb+1) .AND. (j==jolb .OR. j==jolb+1) ) CYCLE
        k = k + 1
        idx (k,:) = [i,j,1,-1]
        fland (k) = fland_tot (i, j)
        hsurf (k) = hsurf_tot (i, j)
        nolbg (k) = MIN (i - 1, j - 1, ie_tot - i, je_tot - j)
        zrlats    = r_startlat_tot + (j-1) * r_dlat
        zcrlat    = COS ( zrlats * r_degrad )
        zdo       = ( (i-zio_tot)*zcrlat *r_dlon/r_dlat )**2 + (j-zjo_tot)**2
        dist (k)  = SQRT( zdo )
        w    (k)  = 0._wp
      ENDDO
    ENDDO
    n_ret = k

!----------------------------------------------------------------------
! Code to be used within DACE sources, applicable to COSMO/ICON/GME/IFS
!----------------------------------------------------------------------
#else

    TYPE(t_grid_cells) ,pointer :: c ! pointer to ICON grid cell

!   get indices of direct neighbours

    call Grid_Indices & !
           (roblon,   & ! <-- geodetic longitude
            roblat,   & ! <-- geodetic latitude
            gen_grid, & ! <-- grid data type
            idx,      & ! --> Grid point indices [Point, index]
            w,        & ! --> Weight
            n_ngb,    & ! --> number of points returned
            mode = 6  ) ! --> number of points requested

!   derive indices of the 'halo'

    d_rad = d2r * gen_grid% d_deg
    select case (gen_grid% gridtype)
    case default
      CALL model_abort (0, 99999, 'grid type not implemented', 'list_of_neighbours')
    case (DWD6_ICON)
      !-----
      ! ICON
      !-----
      zio_tot = 0._wp
      zjo_tot = 0._wp
      n_ret   = n_ngb
      c       => gen_grid% icongrid% patch% cells
      do l = 1, n_ngb
        i = idx (l,1)
        j = idx (l,2)
        do k = 1, size(c% neighbor_idx,3)
          ii = c% neighbor_idx(i,j,k)
          jj = c% neighbor_blk(i,j,k)
          if (      ii <  1              ) cycle
          if ( any (ii == idx(1:n_ngb,1) &
              .and. jj == idx(1:n_ngb,2))) cycle
          n_ret = n_ret + 1
          idx (n_ret,:) = [ii,jj,1,-1]
          w   (n_ret)   = 0._wp
          exit
        end do
      end do

    case (WMO6_LATLON, WMO6_ROTLL, WMO6_GAUSSIAN)
      !------------------
      ! COSMO / HRM / IFS
      !------------------
      if (n_ngb > 0) then
        iolb    = minval (idx (:n_ngb,1))
        jolb    = minval (idx (:n_ngb,2))
        zio_tot = iolb + sum (w(:n_ngb), idx(:n_ngb,1) == iolb+1)
        zjo_tot = jolb + sum (w(:n_ngb), idx(:n_ngb,2) == jolb+1)
        j1      = max (jolb-1, gen_grid% lbg (2))
        je      = min (jolb+2, gen_grid% ubg (2))
        i1      =      iolb-1
        ie      =      iolb+2
        if (.not. gen_grid% global) then
          i1 = max (i1, gen_grid% lbg (1))
          ie = min (ie, gen_grid% ubg (1))
        endif

        k    = n_ngb
        do ii = i1, ie
          i = ii
          if (i < gen_grid% lbg(1)) i=i +(gen_grid% ubg(1) - gen_grid% lbg(1) +1)
          if (i > gen_grid% ubg(1)) i=i -(gen_grid% ubg(1) - gen_grid% lbg(1) +1)
          do j = j1, je
            if ( (i==iolb.or.i==iolb+1) .and. (j==jolb.or.j==jolb+1) ) cycle
            k = k + 1
            idx (k,:) = [i,j,1,-1]
            w    (k)  = 0._wp
          end do
        end do
        n_ret = k
      else
        n_ret = 0
      endif

    END SELECT

!   set remaining entries

    do k = 1, n_ret
      i = idx (k,1)
      j = idx (k,2)
      d = idx (k,3)

      fland (k)  = gen_grid% lsm   (i,j,1,d)
      if (associated (gen_grid% hsurf)) then
        hsurf(k) = gen_grid% hsurf (i,j,1,d)
      else
        hsurf(k) = gen_grid% geosp (i,j,1,d) / gacc
      end if
      x         = unitvector (roblat, roblon)
      zdo       = sqrt (sum ( (x - gen_grid% xnglob (i,j,:,1)) **2))
      zdo       = zdo / d_rad
      dist (k)  = zdo

      select case (gen_grid% gridtype)
      case (DWD6_ICON)
        if (n_ngb > 0) then
          nolbg (k) =  100                      ! Default / fallback value
          if (allocated (c% c_ctrl)) then
            l = c% c_ctrl(i,j)
            if (l > 0) nolbg(k) = l
          end if
        else
          nolbg (k) = -100
        endif
        !------------------------------------------------------------------
        ! catch gridpoints with missing value for surface height (ICON-LAM)
        !------------------------------------------------------------------
        if (hsurf (k) == 9999._wp) nolbg (k) = -100
      case default
        nolbg (k) = min (i - gen_grid% lbg(1),     j - gen_grid% lbg(2),     &
                             gen_grid% ubg(1) - i,     gen_grid% ubg(2) - j  )
      end select
    end do

    if (n_ngb > 0) then
      ibnd = minval (nolbg(1:n_ngb)) + 2
    else
      ibnd = -huge (0)
    end if

#endif

END SUBROUTINE list_of_neighbours

!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for getting global surface fields
!-------------------------------------------------------------------------------

SUBROUTINE get_global_surface ( num_compute     , my_comm_id                   &
                              , ie_loc , je_loc , hsurf     , fr_land          &
                              , ie_tot , je_tot , hsurf_tot , fland_tot )

!-------------------------------------------------------------------------------
! Description:
!   Getting model orography and land/sea mask on total model domain ('global')
!   by gathering and putting together the local parts from all the sub-domains.
!
! Method:
!   Call routine 'gather_field' from module 'parallel_utilities.f90'.
!   Note:
!    - All the required input variables are contained in the parameter list,
!    - but subroutine 'model_abort' from module 'environment' is also called.
!
! Written by        :  Christoph Schraff, DWD  (original version: 20.08.04)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)  ::  &
    ie_loc         ,& ! number of grid pts in zonal direction (local sub-domain)
    je_loc         ,& ! number of grid pts in meridional dir. (local sub-domain)
    ie_tot         ,& ! number of grid pts in zonal direction (in total domain)
    je_tot         ,& ! number of grid pts in meridional dir. (in total domain)
    num_compute    ,& ! number of compute PEs
    my_comm_id        ! rank of this subdomain in the MPI communicator

  REAL    (KIND=wp)        , INTENT (IN)  ::  &
    hsurf     (ie_loc,je_loc) ,& ! model orography (local sub-domain)
    fr_land   (ie_loc,je_loc)    ! land fraction   (local sub-domain)

  REAL    (KIND=wp)        , INTENT (OUT) ::  &
    hsurf_tot (ie_tot,je_tot) ,& ! model orography (total model domain)
    fland_tot (ie_tot,je_tot)    ! land fraction   (total model domain)

! Local variables:
! ----------------
  INTEGER                  :: &
    irm                  ! error status variables
  CHARACTER (LEN=20)       :: &
    yroutine             ! name of this subroutine
  CHARACTER (LEN=30)       :: &
    yerrmsg              ! error message
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin subroutine get_global_surface
!-------------------------------------------------------------------------------

  yroutine = 'get_global_surface'

  IF (num_compute > 1) THEN

    CALL gather_field ( hsurf     , ie_loc , je_loc                            &
                      , hsurf_tot , ie_tot , je_tot , 0, irm )
!   =================

    IF (irm /= 0) WRITE( yerrmsg,'(''ERROR at gather_field'',I6)' ) irm
    IF (irm /= 0) CALL model_abort (my_comm_id, 11141, yerrmsg, yroutine, irm )

    CALL gather_field ( fr_land   , ie_loc , je_loc                            &
                      , fland_tot , ie_tot , je_tot , 0, irm )
!   =================

    IF (irm /= 0) WRITE( yerrmsg,'(''ERROR at gather_field'',I6)' ) irm
    IF (irm /= 0) CALL model_abort (my_comm_id, 11142, yerrmsg, yroutine, irm )

  ELSE
    hsurf_tot = hsurf
    fland_tot = fr_land
  ENDIF

!-------------------------------------------------------------------------------
! End subroutine get_global_surface
!-------------------------------------------------------------------------------

END SUBROUTINE get_global_surface


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for sorting reports accord. to nodes
!-------------------------------------------------------------------------------

SUBROUTINE obs_assign_sort_node ( nrep, nreproc, irprcs, iobstot, jobstot      &
                                , num_comp, i_subpos, nboundlines, my_cart_id  &
                                , irnode, irsort, iobsloc, jobsloc )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" determines for each of
!   an array of observation reports the processor of the subdomain which
!   contains the grid point to which the report is assigned.
!   This is used to distribute the reports read by processor 0 to the
!   appropriate nodes.
!   Therefore, in addition, the indices of the reports are sorted according to
!   the target processors, and optionally, the local grid point coordinates
!   of the reports are also computed.
!
! Method:
!   Reads the report indices and global grid points coordintes of the
!   observation reports and determines the output using the COSMO model
!   topology given by 'i_subpos' and 'nboundlines'.
!   Note:
!    - All the required input variables are contained in the parameter list,
!    - but subroutine 'model_abort' from module 'environment.f90' is used.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)  ::  &
    nrep               ,& ! number  of reports read from the NetCDF file
    nreproc            ,& ! number  of reports to be processed further now
    irprcs  (nreproc+1),& ! indices of reports to be processed further now
    iobstot (nrep)     ,& ! longitudinal \  global indices of grid points to
    jobstot (nrep)        ! latitudinal  /  which the reports have been assigned

  INTEGER                  , INTENT (IN)  ::  &
    num_comp           ,& ! number of compute PEs
    nboundlines        ,& ! number of overlapping boundary lines of subdomains
    my_cart_id         ,& ! rank of this subdomain in the cartesian communicator
    i_subpos (0:num_comp-1,4)  ! positions of the subdomains in the total domain
                               ! (i-, j-indices of lower left and upper right
                               ! grid pt in the order:  i_ll, j_ll, i_ur, j_ur ;
                               ! only the domain interior is considered)

  INTEGER                  , INTENT (OUT) ::  &
    irnode  (nreproc+1),& ! nodes to which the reports will be distributed
    irsort  (nreproc+1)   ! report indices sorted according to 'irnode'

  INTEGER                  , INTENT (OUT) , OPTIONAL ::  &
    iobsloc (nrep)     ,& ! longitudinal \  local indices of the assigned grid
    jobsloc (nrep)        ! latitudinal  /  points in the local sub-domain

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER                  ::  &
    irps               ,& ! index over sorted reports
    inode , irpa , irep   ! loop indices

  CHARACTER (LEN=25)       :: &
    yroutine              ! name of this subroutine
  CHARACTER (LEN=70)       :: &
    yerrmsl               ! error message

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_assign_sort_node
!-------------------------------------------------------------------------------

  yroutine = 'obs_assign_sort_node'

IF (num_comp > 1) THEN

  irps = 0
  DO  inode = 0, num_comp-1
    DO irpa = 1 , nreproc
      irep  = irprcs(irpa)
      IF (      (i_subpos(inode,1) <= iobstot(irep))                           &
          .AND. (i_subpos(inode,3) >= iobstot(irep))                           &
          .AND. (i_subpos(inode,2) <= jobstot(irep))                           &
          .AND. (i_subpos(inode,4) >= jobstot(irep))) THEN
        irps  =  irps + 1
        irnode (irps) = inode
        irsort (irps) = irep
      ENDIF
    ENDDO
  ENDDO

  IF (irps /= nreproc) THEN
    WRITE( yerrmsl,'(" irps =",I6," /= nreproc =",I6)' )  irps, nreproc
    CALL model_abort (my_cart_id, 11222, yerrmsl, yroutine)
  ENDIF

  IF ((PRESENT( iobsloc )) .AND. (PRESENT( jobsloc ))) THEN
    DO irps = 1 , nreproc
      irep  = irsort(irps)
      inode = irnode(irps)
      iobsloc(irep) = iobstot(irep) - i_subpos(inode,1) + 1 + nboundlines
      jobsloc(irep) = jobstot(irep) - i_subpos(inode,2) + 1 + nboundlines
    ENDDO
  ENDIF

ELSE     !   IF (num_comp > 1) THEN

  DO irpa = 1 , nreproc
    irnode (irpa) = my_cart_id
    irsort (irpa) = irprcs(irpa)
  ENDDO

  IF ((PRESENT( iobsloc )) .AND. (PRESENT( jobsloc ))) THEN
    DO irps = 1 , nreproc
      irep  = irsort(irps)
      iobsloc(irep) = iobstot(irep)
      jobsloc(irep) = jobstot(irep)
    ENDDO
  ENDIF

ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_assign_sort_node
!-------------------------------------------------------------------------------
END SUBROUTINE obs_assign_sort_node


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for distrib. reports to local nodes
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_distrib_reports ( itotlen, ibuftot, rbuftot, ybuftot, iylen &
                                   , nodrepn, nodleni, nodlenr, nodleny        &
                                   , nrepl  , nrecvi , nrecvr , nrecvy         &
                                   , iloclen, ibufloc, rbufloc, ybufloc        &
                                   , num_compute , my_cart_id , icomm_cart     &
                                   , imp_integers, imp_reals )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" distributes the
!   observation reports as read from NetCDF files to the different nodes,
!   such that each node obtains exactly those observations which lie in its
!   sub-domain.
!
! Method:
!   First distribute the information to each node which part (offset and length)
!   of the originating long array has to be sent to which node, and then
!   distribute the information. The distributed local information resides in
!   arrays 'ibufloc', 'rbufloc', and 'ybufloc'.
!   Note:
!    - This routine should be called only if (num_compute > 1).
!    - All the required input variables are contained in the parameter list.
!    - However, subroutines 'global_values' and 'distribute_values' from module
!      'parallel_utilities.f90' and 'model_abort' from 'environment.f90' are
!      used (by call).
!
! Initial release: Christoph Schraff, 20.12.07
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)  ::  &
    num_compute    ,& ! number of compute PEs
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_reals      ,& ! REAL      type used for MPI
    imp_integers      ! INTEGER   type used for MPI

  INTEGER                  , INTENT (IN)  ::  &
    iylen                 ,& ! length of character elements
    itotlen (3)           ,& ! length of int/real/char input buffer arrays
    iloclen (3)           ,& ! length of int/real/char output buffer arrays
    ibuftot (itotlen(1))     ! integer   input buffer array

  REAL    (KIND=wp)        , INTENT (IN)  ::  &
    rbuftot (itotlen(2))     ! real      input buffer array

  CHARACTER (LEN=iylen)    , INTENT (IN)  ::  &
    ybuftot (itotlen(3))     ! character input buffer array

  INTEGER                  , INTENT (INOUT)  ::  &
                             ! at node (my_cart_id == 0) the following variables
                             !                           are not modified
    nodrepn (num_compute) ,& ! number of            reports  \   to be
    nodleni (num_compute) ,& ! number of integer   elements   \  distributed
    nodlenr (num_compute) ,& ! number of real      elements   /  to the
    nodleny (num_compute)    ! number of character elements  /   different nodes

  INTEGER                  , INTENT (OUT) ::  &
    nrepl                 ,& ! number of            reports  \  ( to be
    nrecvi                ,& ! number of integer   elements   \   received )
    nrecvr                ,& ! number of real      elements   /  at the
    nrecvy                   ! number of character elements  /   local node

  INTEGER                  , INTENT (OUT) ::  &
    ibufloc (iloclen(1))     ! integer   output buffer array

  REAL    (KIND=wp)        , INTENT (OUT) ::  &
    rbufloc (iloclen(2))     ! real      output buffer array

  CHARACTER (LEN=iylen)    , INTENT (OUT) ::  &
    ybufloc (iloclen(3))     ! character output buffer array

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER                  ::  &
    ntoty             ,& ! total number of character elements
    ilocleny          ,& ! max. number of 1-character elements received locally
    inode, jrp, jch   ,& ! loop indices
    jerr , istat      ,& ! error indicators
    izmplcode            ! error indicator

  CHARACTER (LEN=25)       :: &
    yroutine             ! name of this subroutine
  CHARACTER (LEN=70)       :: &
    yerrmsl              ! error message

! Local arrays:
! ------------

  INTEGER                  ::  &
    nodlenyi (num_compute)    ! number of one-character elements

  INTEGER                  , ALLOCATABLE ::  &
    ibufys   (:)      ,& ! sending buffer for integers converted from characters
    ibufyr   (:)         ! receiveing buffer
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_distrib_reports
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_distrib_reports'

  IF (num_compute > 1) THEN

    ! broadcast the number of reports and (int/real/char) elements for each node

    CALL distribute_values (nodrepn, num_compute,0,imp_integers,icomm_cart,jerr)
    CALL distribute_values (nodleni, num_compute,0,imp_integers,icomm_cart,jerr)
    CALL distribute_values (nodlenr, num_compute,0,imp_integers,icomm_cart,jerr)
    CALL distribute_values (nodleny, num_compute,0,imp_integers,icomm_cart,jerr)
!   ======================

    !   get number of elements to be received locally
    nrepl   = nodrepn (my_cart_id+1)
    nrecvi  = nodleni (my_cart_id+1)
    nrecvr  = nodlenr (my_cart_id+1)
    nrecvy  = nodleny (my_cart_id+1)

    !   distribute integer elements (resulting local array is 'ibufloc')
    ibufloc (:) = 0
    IF (MAXVAL( nodleni ) > 0) THEN

      CALL scatterv_values ( ibuftot, itotlen(1), 0, nodleni                   &
                           , ibufloc, iloclen(1), my_cart_id, num_compute      &
                           , imp_integers, icomm_cart, yerrmsl, izmplcode )
!     ====================

      IF (izmplcode /= 0)                                                      &
        CALL model_abort (my_cart_id, 11162, yerrmsl, yroutine, izmplcode)
    ENDIF

    !   distribute real elements (resulting local array is 'rbufloc')
    rbufloc (:) = 0.0_wp
    IF (MAXVAL( nodlenr ) > 0) THEN

      CALL scatterv_values ( rbuftot, itotlen(2), 0, nodlenr                   &
                           , rbufloc, iloclen(2), my_cart_id, num_compute      &
                           , imp_reals, icomm_cart, yerrmsl, izmplcode )
!     ====================

      IF (izmplcode /= 0)                                                      &
        CALL model_abort (my_cart_id, 11163, yerrmsl, yroutine, izmplcode)
    ENDIF

    !   distribute character elements (resulting local array is 'ybufloc')
    ybufloc (:) = ' '
    IF (MAXVAL( nodleny ) > 0) THEN
      !   characters are converted into integer for scattering by MPI
      ntoty   = 0
      DO inode = 1 , num_compute
        nodlenyi (inode) = nodleny(inode) *iylen
        ntoty  =  ntoty  + nodleny(inode)
      ENDDO
      ilocleny = MAXVAL( nodleny )* iylen
      ALLOCATE ( ibufys (ntoty * iylen)   , STAT=istat )
      ALLOCATE ( ibufyr (ilocleny)        , STAT=istat )
      ibufys (:) = 0
      ibufyr (:) = 0
      IF (my_cart_id == 0) THEN
        DO jrp = 1 , ntoty
          DO jch = 1 , iylen
            ibufys ((jrp-1)*iylen+jch)  =  ICHAR ( ybuftot(jrp)(jch:jch) )
          ENDDO
        ENDDO
      ENDIF

      CALL scatterv_values ( ibufys , ntoty*iylen, 0, nodlenyi                 &
                           , ibufyr , ilocleny   , my_cart_id, num_compute     &
                           , imp_integers, icomm_cart, yerrmsl, izmplcode )
!     ====================

      IF (izmplcode /= 0)                                                      &
        CALL model_abort (my_cart_id, 11164, yerrmsl, yroutine, izmplcode)
      DO jrp = 1 , nrecvy
        DO jch = 1 , iylen
          ybufloc (jrp) (jch:jch)  =  CHAR ( ibufyr((jrp-1)*iylen+jch) )
        ENDDO
      ENDDO
      DEALLOCATE ( ibufys  , STAT=istat )
      DEALLOCATE ( ibufyr  , STAT=istat )
    ENDIF

! ELSE
!   nrepl  = nodrepn(1)
!   nrecvi = nodleni(1)
!   nrecvr = nodlenr(1)
!   nrecvy = nodleny(1)

  ENDIF  ! num_compute > 1

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_distrib_reports
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_distrib_reports


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for distrib. reports to local nodes
!-------------------------------------------------------------------------------

SUBROUTINE obs_gather_buffer_int ( inlen, ibufsize, ibufloc, ilentot, ibuftot  &
                                 , ireceiver, num_compute, my_cart_id          &
                                 , icomm_cart, imp_integers )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" gathers elements
!   from integer buffer arrays from all nodes in a long buffer array at node
!   'ireceiver'.
!
! Method:
!   2 cases have to be distinguished:
!   - If the number of elements ('ilenloc') in the local input buffer array is
!     already known at the call of this routine, it is given in 'inlen' >= 0 .
!   - Otherwise (i.e. if inlen <= -1), it is assumed that the whole sequence of
!     elements in the (input) buffer array 'ibufloc' is subdivided in batches.
!     In each batch, the first element contains the number of elements of the
!     batch. This makes it possible to determine the total number of elements
!     'ilenloc' in the local buffer array.
!   Gathering 'ilenloc' from all nodes, the dimension 'ilentot' of the (global)
!   target buffer array 'ibuftot' can be determined, and the local buffer arrays
!   can be gathered.
!   Note: The pointer 'ibuftot' is allocated in this routine and has to be
!         de-allocated in the calling routine afterwards.
!
! Initial release: Christoph Schraff, 05.08.09
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)    ::  &
    inlen             ,& ! local number of elements to be gathered
                         !   <= -1 : number needs to be determined
                         !   >=  0 : number is already known at call
    ibufsize          ,& ! dimension of buffer array to be gathered
    ibufloc (ibufsize),& ! (local part of) buffer array to be gathered
    ireceiver         ,& ! rank of the receiving PE
    num_compute       ,& ! number of compute PEs
    my_cart_id        ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart        ,& ! communicator for the virtual cartesian topology
    imp_integers         ! INTEGER   type used for MPI

  INTEGER                  , INTENT (OUT)   ::  &
    ilentot              ! total number of elements gathered

  INTEGER                  , POINTER        ::  &
    ibuftot (:)          ! pointer for buffer array which contains the contents
                         !   of the local buffers 'ibufloc' from all nodes

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER                  ::  &
    ilenloc           ,& ! local number of elements to be gathered
    idimtot           ,& ! dimension of 'ibuftot'  (= ilentot+1)
    implcode , istat     ! error indicator

  CHARACTER (LEN=25)       :: &
    yroutine             ! name of this subroutine
  CHARACTER (LEN=40)       :: &
    yerrmsg              ! error message

! Local arrays:
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_gather_buffer_int
!-------------------------------------------------------------------------------

  yroutine = 'obs_gather_buffer_int'

! determine the local number of elements for the current node
  IF (inlen <= -1) THEN
    ilenloc = 0
    Get_next_length:   DO
      IF (ilenloc == ibufsize)                              EXIT Get_next_length
      IF (ibufloc(ilenloc+1) == 0)                          EXIT Get_next_length
      ilenloc = ilenloc + ibufloc(ilenloc+1)
    ENDDO  Get_next_length
  ELSE
    ilenloc = inlen
  ENDIF

  IF (num_compute > 1) THEN
    !   get total number of elements to be gathered
    idimtot  =  ilenloc

    CALL global_values ( idimtot, 1,'SUM', imp_integers, icomm_cart, -1        &
                       , yerrmsg, implcode )
!   ==================

    !   for safety
    idimtot  =  idimtot + 1
    ALLOCATE ( ibuftot (idimtot)    , STAT=istat )

    !   gather the elements from all nodes and get total number of elements
    ilentot  =  0

    CALL gatherv_values ( ibufloc, ilenloc, ibufsize, num_compute              &
                        , ibuftot, ilentot, idimtot , ireceiver                &
                        , imp_integers, icomm_cart, yerrmsg, implcode )
!   ===================

    IF (implcode /= 0)                                                         &
      CALL model_abort (my_cart_id, 11112, yerrmsg, yroutine, implcode)

  ELSE

    ilentot = ilenloc
    ALLOCATE ( ibuftot (ilentot+1)    , STAT=istat )
    ibuftot (1:ilentot) = ibufloc(1:ilentot)
    ibuftot (ilentot+1) = 0

  ENDIF
! lwr_now = (ilentot >= 1)

!-------------------------------------------------------------------------------
! End Subroutine obs_gather_buffer_int
!-------------------------------------------------------------------------------
END SUBROUTINE obs_gather_buffer_int


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for distrib. reports to local nodes
!-------------------------------------------------------------------------------

SUBROUTINE obs_gather_buffer_real ( inlen, ibufsize, rbufloc, ilentot, rbuftot &
                                  , ireceiver, num_compute, my_cart_id         &
                                  , icomm_cart, imp_reals )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" gathers elements
!   from real buffer arrays from all nodes in a long buffer array at node
!   'ireceiver'.
!
! Method:
!   2 cases have to be distinguished:
!   - If the number of elements ('ilenloc') in the local input buffer array is
!     already known at the call of this routine, it is given in 'inlen' >= 0 .
!   - Otherwise (i.e. if inlen <= -1), it is assumed that the whole sequence of
!     elements in the (input) buffer array 'rbufloc' is subdivided in batches.
!     In each batch, the first element contains the number of elements of the
!     batch. This makes it possible to determine the total number of elements
!     'ilenloc' in the local buffer array.
!   Gathering 'ilenloc' from all nodes, the dimension 'ilentot' of the (global)
!   target buffer array 'rbuftot' can be determined, and the local buffer arrays
!   can be gathered.
!   Note: The pointer 'rbuftot' is allocated in this routine and has to be
!         de-allocated in the calling routine afterwards.
!
! Initial release: Christoph Schraff, 29.07.10
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)    ::  &
    inlen             ,& ! local number of elements to be gathered
                         !   <= -1 : number needs to be determined
                         !   >=  0 : number is already known at call
    ibufsize          ,& ! dimension of buffer array to be gathered
    ireceiver         ,& ! rank of the receiving PE
    num_compute       ,& ! number of compute PEs
    my_cart_id        ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart        ,& ! communicator for the virtual cartesian topology
    imp_reals            ! REAL      type used for MPI

  REAL    (KIND=wp)        , INTENT (IN)    ::  &
    rbufloc (ibufsize)   ! (local part of) buffer array to be gathered

  INTEGER                  , INTENT (OUT)   ::  &
    ilentot              ! total number of elements gathered

  REAL    (KIND=wp)        , POINTER        ::  &
    rbuftot (:)          ! pointer for buffer array which contains the contents
                         !   of the local buffers 'rbufloc' from all nodes

! Local parameters:
! ----------------

  REAL    (KIND=wp)        , PARAMETER  :: &
    c0  =  0.0_wp    ! standard real constant 0.0

! Local scalars:
! -------------

  INTEGER                  ::  &
    ilenloc           ,& ! local number of elements to be gathered
    idimtot           ,& ! dimension of 'rbuftot'  (= ilentot+1)
    implcode , istat     ! error indicator

  REAL    (KIND=wp)        ::  &
    rdimtot              ! as 'idimtot' (since imp_integers is not defined here)

  CHARACTER (LEN=25)       :: &
    yroutine             ! name of this subroutine
  CHARACTER (LEN=40)       :: &
    yerrmsg              ! error message

! Local arrays:
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_gather_buffer_real
!-------------------------------------------------------------------------------

  yroutine = 'obs_gather_buffer_real'

! determine the local number of elements for the current node
  IF (inlen <= -1) THEN
    ilenloc = 0
    Get_next_length:   DO
      IF (ilenloc == ibufsize)                              EXIT Get_next_length
      IF (NINT( rbufloc(ilenloc+1) ) == 0)                  EXIT Get_next_length
      ilenloc = ilenloc + NINT( rbufloc(ilenloc+1) )
    ENDDO  Get_next_length
  ELSE
    ilenloc = inlen
  ENDIF

  IF (num_compute > 1) THEN
    !   get total number of elements to be gathered (imp_integers: not defined!)
    rdimtot  =  REAL( ilenloc, wp )

    CALL global_values ( rdimtot, 1,'SUM', imp_reals, icomm_cart, -1           &
                       , yerrmsg, implcode )
!   ==================

    !   for safety
    idimtot  =  NINT( rdimtot ) + 1
    ALLOCATE ( rbuftot (idimtot)    , STAT=istat )

    !   gather the elements from all nodes and get total number of elements
    ilentot  =  0

    CALL gatherv_values ( rbufloc, ilenloc, ibufsize, num_compute              &
                        , rbuftot, ilentot, idimtot , ireceiver                &
                        , imp_reals, icomm_cart, yerrmsg, implcode )
!   ===================

    IF (implcode /= 0)                                                         &
      CALL model_abort (my_cart_id, 11113, yerrmsg, yroutine, implcode)

  ELSE

    ilentot = ilenloc
    ALLOCATE ( rbuftot (ilentot+1)    , STAT=istat )
    rbuftot (1:ilentot) = rbufloc(1:ilentot)
    rbuftot (ilentot+1) = c0

  ENDIF
! lwr_now = (ilentot >= 1)

!-------------------------------------------------------------------------------
! End Subroutine obs_gather_buffer_real
!-------------------------------------------------------------------------------
END SUBROUTINE obs_gather_buffer_real


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for superobbing of high-res profiles
!-------------------------------------------------------------------------------

SUBROUTINE obs_average_layers ( nlint, p_int, nvar, nflg, nlev                 &
                              , kpp, kplv, lflg0, rmdi, rmdich                 &
                              , mrepflg , repob , msupflg, supob )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" averages the observed
!   high-resolution profiles within pre-specified interpolation layers in order
!   to obtain a superobbed profile.
!   If any of the observations within the interpolation layer is passive, then
!   the superobs is also flagged as passive.
!
! Method:
!   Interpolation in log( pressure ) of the observed values to the pre-specified
!   interpolation levels, then averaging within each of the interpolation layer
!   using the values at the upper and lower bounding interpolation levels and
!   the observed values within the layer. The averaging also uses log( pressure)
!   as vertical coordinate.
!   Note:
!    - All the required input variables are contained in the parameter list,
!    - All the variables are treated independently from each other,
!      which implies that if an obs is set to missing value, the corresponding
!      obs error is not set automatically to missing value as well !!!
!    - The index 'kplv' for the pressure used for the interpolation 'kplv'
!      should be set to 'nbtplv' (not 'nbtp') in the call of this routine.
!      This pressure is the reported pressure or otherwise if height is reported
!      (e.g. wind profiler), it is pressure derived from reported height using
!      the standard atmosphere. In this way, the superobbing levels (in terms
!      of the reported level type) become independent from the model state
!      (and hence identical for all members of an ensemble).
!   Comment:
!     This routine is partly derived from the DACE routine 'average_layers'.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)    ::  &
    nlint              ,& ! number  of interpolation levels
    nvar               ,& ! number  of variables (in 'repob'  , 'supob')
                          !     (note: by setting 'nvar = mxrbdy' in the call of
                          !            this routine, the superobbing is not only
                          !            done for the observed values, but e.g.
                          !            also for the observation errors)
    nflg               ,& ! number  of variables (in 'mrepflg', 'msupflg')
    nlev               ,& ! number  of observation levels
                          ! 2 indices of pressure in 'repob', 'supob'
                          !   where pressure is either reported or otherwise:
    kpp                ,& !   - derived from height using the model state
    kplv                  !   - derived from height using the standard atmos.

  REAL    (KIND=wp)        , INTENT (IN)    ::  &
    p_int  (nlint)     ,& ! interpolation levels (= upper / lower boundaries
                          !                         of superobbing layers)
    rmdi               ,& ! commonly used missing data indicator (huge negative)
    rmdich                ! check value for missing data (= -1.E30_wp)

  LOGICAL                  , INTENT (IN)    ::  &
    lflg0        (nflg)   ! = .true. : use 'IOR'  over flags for superobbing
                          ! = .false.: use 'IAND' over flags for superobbing

  INTEGER                  , INTENT (IN)    ::  &
    mrepflg(nlev ,nflg)   ! flag words of original report body

  REAL    (KIND=wp)        , INTENT (INOUT) ::  &
    repob  (nlev ,nvar)   ! original report body (with attribute 'inout' as p in
                          !   entry 'kplv' is temporarily replaced by log(p))

  INTEGER                  , INTENT (OUT)   ::  &
    msupflg(nlint,nflg)   ! flag words of superobs

  REAL    (KIND=wp)        , INTENT (OUT)   ::  &
    supob  (nlint,nvar)   ! superobbed report body

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER                  ::  &
    klpp               ,& ! index of log(pressure) in 'repob', 'supob' (=kpp)
    klp                ,& ! index of log(pressure) in 'repob', 'supob' (=kplv)
                          !   used for the superobbing
    kob                ,& ! index of obs level closest to interpolation level
    kob2               ,& ! index of other obs level used for interpolation
    kobf1    (nlint)   ,& ! indices of obs level closest below interpol. level
    kobf2    (nlint)   ,& ! indices of obs level closest above interpol. level
    ii                 ,& ! loop index over interpolation levels
    kk                 ,& ! loop index over observation levels in int'pol. layer
    ivar               ,& ! loop index over variables (in 'repob')
    iflg                  ! loop index over variables (in 'mrepflg')

  REAL    (KIND=wp)        ::  &
    lnp_int   (nlint)  ,& ! interpolation levels
    dlnp               ,& ! difference in log( pressure )
    w1    , w2   , ww  ,& ! interpolation weights
    wsum                  ! sum of interpolation weights

  LOGICAL                  ::  &
    lval      (nlint)     ! superob to be derived for layer 'ii'
                          !   (even without an original obs within its bounds)

! CHARACTER (LEN=25)       :: &
!   yroutine              ! name of this subroutine
! CHARACTER (LEN=70)       :: &
!   yerrmsl               ! error message

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_average_layers
!-------------------------------------------------------------------------------

! yroutine = 'obs_average_layers'
  supob   (:,:) = rmdi
  msupflg (:,:) = 0
  IF (nlint < 2)                                                          RETURN

  !   it is assumed that levels with invalid pressure levels have been removed
  !   previously:
  !   temporarily replace pressure by log( pressure ), because log( pressure )
  !   (rather than pressure) is averaged (like the observations of temperature,
  !   wind, etc.) over log( pressure ) in order to obtain the pressure of the
  !   superobbed layer
  klpp              =  kpp
  klp               =  kplv
  repob   (:,klpp)  =  LOG( repob(:,kpp ) )
  repob   (:,klp )  =  LOG( repob(:,kplv) )
  lnp_int (:)       =  LOG( p_int(:) )
  DO ii = 1 , nlint
    !   get obs level 'kob' nearest to interpolation level 'ii'
    kob   =  MINLOC( ABS( lnp_int(ii) - repob(:  ,klp) ), DIM=1 )
    dlnp  =               lnp_int(ii) - repob(kob,klp)
    kob2  =  kob

    ! linear interpolation to interpolation levels
    ! --------------------------------------------
    supob (ii,:)  =  repob (kob,:)
!   IF (ABS( dlnp ) > 0.0001_wp) THEN
    IF (ABS( dlnp ) > 0.00001_wp) THEN
      !   if (dlnp > 0.) then interpolate betw. kob, kob-1=kob2
      !   if (dlnp < 0.) then interpolate betw. kob, kob+1=kob2
      kob2  =  MAX( MIN( kob - NINT( SIGN( 1._wp , dlnp ) ) , nlev ) , 1 )
      IF (      (kob2 /= kob)                                                  &
          .AND. (ABS( repob(kob,klp)-repob(kob2,klp) ) > .00001_wp)) THEN
        w1  =    (lnp_int(ii)    - repob(kob2,klp))                            &
               / (repob(kob,klp) - repob(kob2,klp))
        w2  =  1._wp - w1
        DO ivar = 1 , nvar
          IF (MIN( repob(kob,ivar), repob(kob2,ivar) ) <= rmdich) THEN
            supob (ii,ivar)  =  rmdi
          ELSE
            supob (ii,ivar)  =  w1 * repob(kob,ivar) + w2 * repob(kob2,ivar)
          ENDIF
        ENDDO
!     ELSE  ! this case is already done above: supob(ii,:) = repob(kob,:)
!       w1  =  1._wp
!       w2  =  0._wp
      ENDIF
    ENDIF
    ! for superobbing layer 'ii', original levels from kobf1(ii) to
    ! kobf2(ii+1) are used; this info is required for missing values and flags
    kobf1 (ii)  =  MIN( kob, kob2 )
    kobf2 (ii)  =  MAX( kob, kob2 )
!   WRITE(0,*)             'ZZs ', nlint, ii, kobf1(ii), kobf2(ii), kob, kob2  &
!                                , exp( supob(ii,5) ), supob(ii,3), w1,w2, dlnp
  ENDDO

  ! to account for obs level gaps in original obs profile:
  ! criterion that a superob is derived for layer 'ii':
  ! ------------------------------------------------------
  ! (bounds of superob layer 'ii' are interpolation levels 'ii', 'ii+1')
  DO ii = 2 , nlint - 2
    ! an obs level must lie within the superobbing layer
!   lval (ii) =    (kobf2(ii+1) - kobf1(ii) > 1)
    ! an obs level must lie within this or a neighbouring superobbing level
!   lval (ii) =    (kobf2(ii+1) - kobf1(ii) > 1)                               &
!              .OR.(repob(kobf1(ii  ),klp) <       lnp_int(ii-1))              &
!              .OR.(repob(kobf2(ii+1),klp) >       lnp_int(ii+2))
    ! an obs level must lie within this superobbing layer
    !                    or within half of a neighbouring superobbing level
    lval (ii) =    (kobf2(ii+1) - kobf1(ii) > 1)                               &
               .OR.(repob(kobf1(ii  ),klp) <.5_wp*(lnp_int(ii-1)+lnp_int(ii))) &
               .OR.(repob(kobf2(ii+1),klp) >.5_wp*(lnp_int(ii+2)+lnp_int(ii+1)))
  ENDDO
  !  (the top and lowest level are always derived (if possible))
  lval (nlint-1) = .TRUE.
  lval (1)       = .TRUE.

  ! derive (nlint-1) layer averages
  ! -------------------------------
  DO ii = 1 , nlint - 1
    kob   =  MAXLOC( repob(:,klp), MASK=(repob(:,klp) <= lnp_int(ii  )), DIM=1 )
    kob2  =  MINLOC( repob(:,klp), MASK=(repob(:,klp) >= lnp_int(ii+1)), DIM=1 )
!NEC$ nomove
    DO ivar = 1 , nvar
      wsum  =  0._wp
      IF (MIN( supob(ii ,ivar), supob(ii+1,ivar) ) <= rmdich) THEN
        supob (ii,ivar)  =  rmdi
      ELSEIF (kob <= kob2) THEN
        IF (MINVAL( repob(kob:kob2,ivar) ) <= rmdich) THEN
          supob (ii,ivar)  =  rmdi
        ELSE
          !   if at least 1 obs level lies within the superobbing layer
          w1  =  ABS( lnp_int(ii) - repob(kob,klp) )
          supob (ii,ivar)  =                   w1 * supob(ii ,ivar)
          supob (ii,ivar)  =  supob(ii,ivar) + w1 * repob(kob,ivar)
          wsum = wsum + w1
          DO kk = kob , kob2-1
            ww  =  ABS( repob(kk,klp) - repob(kk+1,klp) )
            supob (ii,ivar)  =  supob(ii,ivar) + ww * repob(kk  ,ivar)
            supob (ii,ivar)  =  supob(ii,ivar) + ww * repob(kk+1,ivar)
            wsum = wsum + ww
          ENDDO
          w2  =  ABS( repob(kob2,klp) - lnp_int(ii+1) )
          supob (ii,ivar)  =  supob(ii,ivar) + w2 * repob(kob2,ivar)
          supob (ii,ivar)  =  supob(ii,ivar) + w2 * supob(ii+1,ivar)
          wsum = wsum + w2
          IF (wsum >  epsy)  supob (ii,ivar)  =  supob(ii,ivar) / wsum * 0.5_wp
          IF (wsum <= epsy)  supob (ii,ivar)  =  rmdi
        ENDIF
      ELSEIF (lval(ii)) THEN
        ww  =  ABS( lnp_int(ii) - lnp_int(ii+1) )
        supob (ii,ivar)  =                   ww * supob(ii  ,ivar)
        supob (ii,ivar)  =  supob(ii,ivar) + ww * supob(ii+1,ivar)
        wsum = wsum + ww
        IF (wsum >  epsy)  supob (ii,ivar)  =  supob(ii,ivar) / wsum * 0.5_wp
        IF (wsum <= epsy)  supob (ii,ivar)  =  rmdi
      ELSE
        supob (ii,ivar)  =  rmdi
      ENDIF
      ! set missing value as soon as 1 of the original values is missing value
      ! (in which case the above algorithm for averaging does not work properly)
      DO kk = kobf1(ii) , kobf2(ii+1)
        IF (repob(kk,ivar) <= rmdich)  supob (ii,ivar) = rmdi
      ENDDO
    ENDDO
  ENDDO
  !   convert log( pressure ) (back) to pressure
  repob   (:,kplv)  =  EXP( repob(:,klp ) )
  supob   (:,kplv)  =  EXP( supob(:,klp ) )
  repob   (:,kpp )  =  EXP( repob(:,klpp) )
  supob   (:,kpp )  =  EXP( supob(:,klpp) )

  ! set flag words for superobs
  ! ---------------------------
  ! for each bit in the flag words for the superobbing layer 'ii':
  ! consider the corresponding bits from all original levels
  !   (i.e. from kobf1(ii) to kobf2(ii+1))
  ! which have been used to derive the superobs, and
  ! use the IAND or IOR condition to obtain the bits in superobs flag word
! msupflg (:,:) = 0
  DO ii = 1 , nlint - 1
    DO iflg = 1 , nflg
      msupflg (ii,iflg)  =  mrepflg(kobf1(ii),iflg)
      DO kk = kobf1(ii+1) , kobf2(ii+1)
        IF (lflg0(iflg)) THEN
          msupflg (ii,iflg)  =  IOR ( msupflg(ii,iflg) , mrepflg(kk,iflg) )
        ELSE
          msupflg (ii,iflg)  =  IAND( msupflg(ii,iflg) , mrepflg(kk,iflg) )
        ENDIF
      ENDDO
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_average_layers
!-------------------------------------------------------------------------------
END SUBROUTINE obs_average_layers


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for removing duplicate levels
!-------------------------------------------------------------------------------

SUBROUTINE obs_rm_duplicate_levels ( nlev, nvar, nflg, nbpl, rmdich            &
                                   , zmlbdy, mzmlbd , nrejlv )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" removes duplicate
!   levels in multi-level reports.
!
! Method:
!   At first remove levels with missing pressure.
!   Then remove the second of any pair of subsequent levels which are considered
!   duplicate if the difference in level pressure is less than 0.1 hPa.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)    ::  &
    nlev               ,& ! number of levels to be checked in multi-level report
    nvar               ,& ! number of 'variables' (entries) in 'zmlbdy'
    nflg               ,& ! number of 'variables' (entries) in 'mzmlbd'
    nbpl                  ! index in 'zmlbdy' for pressure (level)

  REAL    (KIND=wp)        , INTENT (IN)    ::  &
    rmdich                ! check value for missing data (= -1.E30_wp)

  REAL    (KIND=wp)        , INTENT (INOUT) ::  &
    zmlbdy (nlev ,nvar)   ! original report body

  INTEGER                  , INTENT (INOUT) ::  &
    mzmlbd (nlev ,nflg)   ! flag words of original report body

  INTEGER                  , INTENT (OUT)   ::  &
    nrejlv                ! number of rejected levels

! Local parameters: None

! Local scalars:

  INTEGER                  ::  &
    nrejlp             ,& ! number of levels with missing pressure
    nlevp              ,& ! number of levels with non-missing pressure
    klev                  ! loop index

! Local arrays: None
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_rm_duplicate_levels
!-------------------------------------------------------------------------------

  ! first remove levels with missing pressure
  nrejlp = 0
  DO klev = 1 , nlev
    IF (zmlbdy(klev,nbpl) < rmdich) THEN
      nrejlp = nrejlp + 1
    ELSEIF (nrejlp >= 1) THEN
      zmlbdy (klev-nrejlp,:) = zmlbdy(klev,:)
      mzmlbd (klev-nrejlp,:) = mzmlbd(klev,:)
    ENDIF
  ENDDO
  nlevp = nlev - nrejlp

  ! then remove duplicate levels
  nrejlv = 0
  DO klev = 2 , nlevp
    IF (zmlbdy(klev,nbpl) >= zmlbdy(klev-1,nbpl)-0.1_wp) THEN
      nrejlv = nrejlv + 1
    ELSEIF (nrejlv >= 1) THEN
      zmlbdy (klev-nrejlv,:) = zmlbdy(klev,:)
      mzmlbd (klev-nrejlv,:) = mzmlbd(klev,:)
    ENDIF
  ENDDO

  nrejlv = nrejlv + nrejlp

!-------------------------------------------------------------------------------
! End Subroutine obs_rm_duplicate_levels
!-------------------------------------------------------------------------------
END SUBROUTINE obs_rm_duplicate_levels


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for standard atmosphere information
!-------------------------------------------------------------------------------

SUBROUTINE std_atmosphere ( yvar_in, val_in, yvar_out, val_out, r_g, r_d, irm )

!-------------------------------------------------------------------------------
! Description:
!   Depending on the input variable and selected output variable, the output
!   value is determined at the level given by the input value using the 'US1976'
!   standard atmosphere.
!   Allowed combinations of input und output variables (yvar_in, yvar_out) :
!     Input     -->   yvar_in =           , Output     -->    yvar_out =
!     height [m]      'h','H','z' or 'Z'  , pressure    [Pa]  'p' or 'P'
!     height [m]      'h','H','z' or 'Z'  , temperature [K]   't' or 'T'
!     height [m]      'h','H','z' or 'Z'  , lapse rate  [K/m] 'b','B','g' or 'G'
!     pressure [Pa]   'p' or 'P'          , height      [m]   'h','H','z' or 'Z'
!     pressure [Pa]   'p' or 'P'          , temperature [K]   't' or 'T'
!     pressure [Pa]   'p' or 'P'          , lapse rate  [K/m] 'b','B','g' or 'G'
!     temperature [K] 't' or 'T'          , height      [m]   'h','H','z' or 'Z'
!     temperature [K] 't' or 'T'          , pressure    [Pa]  'p' or 'P'
!   Allowed range of input values (var_in)
!     (and of output values var_out if input is temperature) :
!                                          -800.0  <=  height
!                                              .0  <   pressure     <=  120000.0
!     Min. temperature in std. atmosphere = 186.87 <=  temperature  <   350.0
!   Note: If temperature is the input variable, the output value (height or
!         pressure) is not unique. In this case, the variable for the output
!         value ('val_out') has to contain an (input) value already at the
!         call of the routine. Then, the target level is the nearest
!         temperature level according to the input value above the level
!         given by 'val_out'.
!   Written by: C. Schraff, equations based on F77 code (dwdlib) by U. Voigt.
!
! Method:
!   Atmosphere used: 'U.S. Standard Atmosphere 1976' (NOAA, NASA)
!                    (Ref: NOAA-S/T 76-1562 Washington D.C., October 1976)
!   This atmosphere consists of a sequence of different layers.
!   Equations used:
!   - polytropic layers (dT/dz=const=-b_0): p(z) = p_0 *(1- b_0*z/T_0)^(b_h/b_0)
!   - isothermal layers (dT/dz=0)         : p(z) = p_0 *exp( (-g*z)/(r_d*T_0) )
!   '_0' indicates: at the lower boundary of the layer;
!   b  : lapse rate;   b_h: lapse rate of a homogeneous atmosphere;
!   r_d: gas constant for dry air.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments
! --------------------

  CHARACTER (LEN= 1)       , INTENT (IN)    ::       &
    yvar_in     ,& ! variable of input value
    yvar_out       ! variable of output value

  REAL    (KIND=wp)        , INTENT (IN)    ::       &
    r_g         ,& ! acceleration due to gravity
    r_d         ,& ! gas constant for dry air
    val_in         ! input value defining the vertical level

  REAL    (KIND=wp)        , INTENT (INOUT) ::       &
    val_out        ! output value (at the level given by val_in)

  INTEGER                  , INTENT (OUT)   ::       &
    irm            ! error status = 0 : no error
                   !              = 1 : invalid input variable yvar_in
                   !              = 2 : invalid output variable yvar_out
                   !              = 3 : invalid input value for val_in
                   !              = 4 : invalid input value for val_out

! Local parameters:
! ----------------

  REAL    (KIND=wp)        , PARAMETER  :: &
    c0     =  0.0_wp     ,& ! standard real constant 0.0
!   c1     =  1.0_wp     ,& ! standard real constant 1.0
!   epsy   =  1.E-8_wp   ,& ! small tolerance
    z_min  =    -800._wp ,& ! mimimum allowed value for height
    p_max  =  120000._wp ,& ! maximum allowed value for pressure
    t_max  =     350._wp    ! maximum allowed value for temperature
                                !   note: 'p_max', 't_max' should be specified
                                !         such that in the standard atmosphere,
                                !         level 'p_max' lies below level 't_max'

  INTEGER                  , PARAMETER  :: &
    nsalev = 8    ! number of levels in the standard atmosphere definition table

  REAL    (KIND=wp)        , PARAMETER  :: &
      ! definition of standard atmosphere: height z [m] / temperature T [K] /
      !                                    lapse rate b [K/m] / pressure p [hPa]
      !                                    at the lower boundary of each layer
    std_atm (4,nsalev) = RESHAPE( (/                                           &
      !       z_0       T_0        b_0           p_0
                0._wp,  288.15_wp, -.0065_wp, 101325.0000_wp,  &
            11000._wp,  216.65_wp,  c0      ,  22632.0000_wp,  &
            20000._wp,  216.65_wp,  .0010_wp,   5474.8000_wp,  &
            32000._wp,  228.65_wp,  .0028_wp,    868.0100_wp,  &
            47000._wp,  270.65_wp,  c0      ,    110.9000_wp,  &
            51000._wp,  270.65_wp, -.0028_wp,     66.9380_wp,  &
            71000._wp,  214.65_wp, -.0020_wp,      3.9564_wp,  &
            86000._wp,  186.87_wp,  c0      ,       .3733_wp/) &
                                , (/4,nsalev/) )

      !       0., 11000., 20000., 32000., 47000., 51000., 71000., 86000.,  & ! H
      !   288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.87,  & ! T
      !   -.0065,     .0,   .001,  .0028,     .0, -.0028,  -.002,     .0,  & ! b
      !  101325., 22632., 5474.8, 868.01,  110.9, 66.938, 3.9564,  .3733/) & ! p

! Local variables
! ---------------

  INTEGER                  ::       &
    ilvb        ,& ! level index for lower boundary of layer containing
                   !   value 'val_in'
    ilvbb       ,& ! level index for layer containing value 'val_ref'
    ilv  , itry    ! loop  indices

  REAL    (KIND=wp)        ::       &
    val_ref     ,& ! conditional input value (lower limit for search of target
                   !   level if input variable is temperature)
    b_hom       ,& ! lapse rate of a homogeneous atmosphere [K/m]
    dz          ,& ! height diff. between target level and level 'ilvb'
    z_0         ,& ! height      at level 'ilvb'
    t_0         ,& ! temperature at level 'ilvb'
    b_0         ,& ! lapse rate  at level 'ilvb'
    p_0            ! pressure    at level 'ilvb'

  LOGICAL                                ::       &
    lval_ok        ! result 'val_out' also meets condition from 'val_ref'

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE std_atmosphere
!-------------------------------------------------------------------------------

  irm   = 0
  b_hom = - r_g / r_d

! input variable: height  (output: pressure, temperature, or lapse rate)
! ----------------------
  IF (     (yvar_in == 'h') .OR. (yvar_in == 'H')                              &
      .OR. (yvar_in == 'z') .OR. (yvar_in == 'Z')) THEN
    IF (val_in >= z_min) THEN
      ! find (lower) level index 'ilvb' of layer which contains 'val_in'
      ilvb = 1
      DO ilv = 1 , nsalev
        IF (val_in >= std_atm(1,ilv)-epsy)  ilvb = ilv
      ENDDO
      ! compute 'val_out' by extrapolation from level 'ilvb' to height 'val_in'
      dz  = val_in - std_atm(1,ilvb)
      t_0 = std_atm(2,ilvb)
      b_0 = std_atm(3,ilvb)
      p_0 = std_atm(4,ilvb)
      IF     ((yvar_out == 'p') .OR. (yvar_out == 'P')) THEN
        IF (ABS( b_0 ) < epsy) THEN   ! isothermal atmosphere
          val_out = p_0 * EXP( -dz *r_g /(r_d *t_0) )
        ELSE                          ! polytropic atmosphere
          val_out = p_0 * (c1 + dz *b_0 /t_0) **(b_hom /b_0)
        ENDIF
      ELSEIF ((yvar_out == 't') .OR. (yvar_out == 'T')) THEN
        IF (ABS( b_0 ) < epsy) THEN   ! isothermal atmosphere
          val_out = t_0
        ELSE                          ! polytropic atmosphere
          val_out = t_0 + dz* b_0
        ENDIF
      ELSEIF (     (yvar_out == 'b') .OR. (yvar_out == 'B')                    &
              .OR. (yvar_out == 'g') .OR. (yvar_out == 'G')) THEN
        val_out = b_0
      ELSE
        irm = 2
      ENDIF
    ELSE
      irm = 3
    ENDIF

! input variable: pressure  (output: height, temperature, or lapse rate)
! ------------------------
  ELSEIF ((yvar_in == 'p') .OR. (yvar_in == 'P')) THEN
    IF ((val_in > epsy) .AND. (val_in <= p_max)) THEN
      ! find (lower) level index 'ilvb' of layer which contains 'val_in'
      ilvb = 1
      DO ilv = 1 , nsalev
        IF (val_in <= std_atm(4,ilv)+epsy)  ilvb = ilv
      ENDDO
      ! compute 'val_out' by extrapolat. from level 'ilvb' to pressure 'val_in'
      z_0 = std_atm(1,ilvb)
      t_0 = std_atm(2,ilvb)
      b_0 = std_atm(3,ilvb)
      p_0 = std_atm(4,ilvb)
      IF (     (yvar_out == 'h') .OR. (yvar_out == 'H')                        &
          .OR. (yvar_out == 'z') .OR. (yvar_out == 'Z')) THEN
        IF (ABS( b_0 ) < epsy) THEN   ! isothermal atmosphere
          val_out = z_0 - (r_d/r_g *t_0) * LOG( val_in /p_0 )
        ELSE                          ! polytropic atmosphere
          val_out = z_0 - t_0 /b_0 *(c1 - (val_in /p_0)**(b_0/b_hom))
        ENDIF
      ELSEIF ((yvar_out == 't') .OR. (yvar_out == 'T')) THEN
        IF (ABS( b_0 ) < epsy) THEN   ! isothermal atmosphere
          val_out = t_0
        ELSE                          ! polytropic atmosphere
          val_out = t_0 * (val_in /p_0) **(b_0/b_hom)
        ENDIF
      ELSEIF (     (yvar_out == 'b') .OR. (yvar_out == 'B')                    &
              .OR. (yvar_out == 'g') .OR. (yvar_out == 'G')) THEN
        val_out = b_0
      ELSE
        irm = 2
      ENDIF
    ELSE
      irm = 3
    ENDIF

! input variable: temperature  (output: height or pressure)
! ---------------------------
  ELSEIF ((yvar_in == 't') .OR. (yvar_in == 'T')) THEN

    ! find (lower) level index 'ilvbb' of layer which contains 'val_ref'
    IF ((val_in >= -MAXVAL( -std_atm(2,:) )) .AND. (val_in < t_max)) THEN
      ilvbb   = 0
      IF (     (yvar_out == 'h') .OR. (yvar_out == 'H')                        &
          .OR. (yvar_out == 'z') .OR. (yvar_out == 'Z')) THEN
        val_ref = MAX( val_out , z_min )
        DO ilv = 1 , nsalev
          IF (val_ref >= std_atm(1,ilv)-epsy)  ilvbb = ilv
        ENDDO
      ELSEIF ((yvar_out == 'p') .OR. (yvar_out == 'P')) THEN
        val_ref = val_out
        IF ((val_ref > p_max) .OR. (val_ref <= epsy))  val_ref = p_max
        DO ilv = 1 , nsalev
          IF (val_ref <= std_atm(4,ilv)+epsy)  ilvbb = ilv
        ENDDO
      ELSE
        irm = 2
      ENDIF
    ELSE
      irm = 3
    ENDIF

    IF (irm == 0) THEN
      try_loop: DO itry = 1 , 2
        ! find level index 'ilvb' ( >= ilvbb) of layer which contains 'val_in'
        ilvb = 0
        DO ilv = nsalev-1 , MAX( ilvbb,1 ) , -1
          IF (     (      (val_in >= std_atm(2,ilv  )-epsy)                    &
                    .AND. (val_in <  std_atm(2,ilv+1)+epsy))                   &
              .OR. (      (val_in <= std_atm(2,ilv  )-epsy)                    &
                    .AND. (val_in >  std_atm(2,ilv+1)+epsy)))  ilvb = ilv
        ENDDO
          ! assuming negative lapse rate below level 1
        IF ((ilvbb == 0) .AND. (val_in >= std_atm(2,1)-epsy))  ilvb = 1
        IF (ilvb == 0)  irm = 4
        IF (ilvb == 0)                                             EXIT try_loop
        ! compute 'val_out' within layer 'ilvb'
        z_0 = std_atm(1,ilvb)
        t_0 = std_atm(2,ilvb)
        b_0 = std_atm(3,ilvb)
        p_0 = std_atm(4,ilvb)
        IF (     (yvar_out == 'h') .OR. (yvar_out == 'H')                      &
            .OR. (yvar_out == 'z') .OR. (yvar_out == 'Z')) THEN
          IF (ABS( b_0 ) < epsy) THEN   ! isothermal atmosphere
            val_out = z_0
          ELSE                          ! polytropic atmosphere
            val_out = z_0 + (val_in - t_0) /b_0
          ENDIF
          ! 'val_out' is valid only if >= 'val_ref'
          lval_ok  =  (val_out >= val_ref-epsy)
        ELSEIF ((yvar_out == 'p') .OR. (yvar_out == 'P')) THEN
          IF (ABS( b_0 ) < epsy) THEN   ! isothermal atmosphere
            val_out = p_0
          ELSE                          ! polytropic atmosphere
            val_out = p_0 * (val_in /t_0) **(b_hom/b_0)
          ENDIF
          ! 'val_out' is valid only if >= 'val_ref'
          lval_ok  =  (val_out <= val_ref+epsy)
        ENDIF
        IF (lval_ok)                                               EXIT try_loop
        ! after the first try, start search once again 1 layer above 'ilvb'
        IF ((.NOT. lval_ok) .AND. (itry == 1))  ilvbb = ilvb + 1
        IF ((.NOT. lval_ok) .AND. (itry == 2))  irm   = 4
      ENDDO try_loop
    ENDIF

! other input variable: invalid
! -----------------------------
  ELSE
    irm = 1
  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine std_atmosphere
!-------------------------------------------------------------------------------
END SUBROUTINE std_atmosphere


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for standard atmosphere information
!-------------------------------------------------------------------------------

SUBROUTINE std_atm_p2z_fast ( nlev, nlvp, isrtlvp, p_in, r_g, r_d, z_miss      &
                            , z_out )

!-------------------------------------------------------------------------------
! Description:
!   Converts a pressure profile into a height profile based on the 'US1976'
!   standard atmosphere.
!   (This could also be done by subr. 'std_atmosphere', but is done here far
!    more efficiently. Note that height is computed here only for pressure
!    values  110.9 < p <= 120000.0 Pa.)
!
! Method:
!   Atmosphere used: 'U.S. Standard Atmosphere 1976' (NOAA, NASA)
!                    (Ref: NOAA-S/T 76-1562 Washington D.C., October 1976)
!   This atmosphere consists of a sequence of different layers.
!   Equations used:
!   - polytropic layers (dT/dz=const=-b_0): p(z) = p_0 *(1- b_0*z/T_0)^(b_h/b_0)
!   - isothermal layers (dT/dz=0)         : p(z) = p_0 *exp( (-g*z)/(r_d*T_0) )
!   '_0' indicates: at the lower boundary of the layer;
!   b  : lapse rate;   b_hom: lapse rate of a homogeneous atmosphere;
!   r_d: gas constant for dry air.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments
! --------------------

  INTEGER                  , INTENT (IN)    ::       &
    nlev           ,& ! total number of levels in report (assumed to be > 0)
    nlvp           ,& ! number of levels to be processed (assumed to be > 0)
    isrtlvp (nlev)    ! sorted list of indices of levels to be processed

  REAL    (KIND=wp)        , INTENT (IN)    ::       &
    r_g            ,& ! acceleration due to gravity
    r_d            ,& ! gas constant for dry air
    z_miss         ,& ! missing height value (if (p_in > p_max)
                      !                       or (p_in < 110.9Pa))
    p_in    (nlev)    ! pressure of levels to be processed

  REAL    (KIND=wp)        , INTENT (INOUT) ::       &
    z_out   (nlev)    ! standard height at pressure levels

! Local parameters:
! ----------------

  REAL    (KIND=wp)        , PARAMETER  :: &
    c0     =  0.0_wp     ,& ! standard real constant 0.0
!   c1     =  1.0_wp     ,& ! standard real constant 1.0
!   epsy   =  1.E-8_wp   ,& ! small tolerance
!   z_min  =    -800._wp ,& ! mimimum allowed value for height
    p_max  =  120000._wp    ! maximum allowed value for pressure
!   z_miss =   -9999._wp    ! missing height value (if (p_in > p_max)
                            !                       or (p_in < 110.9Pa))

  INTEGER                  , PARAMETER  :: &
    nsalev = 5    ! number of levels in the standard atmosphere definition table

  REAL    (KIND=wp)        , PARAMETER  :: &
      ! definition of standard atmosphere: height z [m] / temperature T [K] /
      !   (up to 1.109 hPa = 47km)         lapse rate b [K/m] / pressure p [hPa]
      !                                    at the lower boundary of each layer
    std_atm (4,nsalev) = RESHAPE( (/                                           &
      !       z_0       T_0        b_0           p_0
                0._wp,  288.15_wp, -.0065_wp, 101325.0000_wp,  &
            11000._wp,  216.65_wp,  c0      ,  22632.0000_wp,  &
            20000._wp,  216.65_wp,  .0010_wp,   5474.8000_wp,  &
            32000._wp,  228.65_wp,  .0028_wp,    868.0100_wp,  &
            47000._wp,  270.65_wp,  c0      ,    110.9000_wp/) &
                                , (/4,nsalev/) )

      !       0., 11000., 20000., 32000., 47000., 51000., 71000., 86000.,  & ! H
      !   288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.87,  & ! T
      !   -.0065,     .0,   .001,  .0028,     .0, -.0028,  -.002,     .0,  & ! b
      !  101325., 22632., 5474.8, 868.01,  110.9, 66.938, 3.9564,  .3733/) & ! p

! Local variables
! ---------------

  INTEGER                  ::       &
    ilvb (nlvp) ,& ! std_atm level index for lower boundary of layer containing
                   !   p_in
    klev        ,& ! loop index over levels
    ilev           ! index of vertical level to be processed

  REAL    (KIND=wp)        ::       &
    bi_hom      ,& ! inverse of lapse rate of a homogeneous atmosphere [K/m]
    z_0         ,& ! height      at level 'ilvb'
    t_0         ,& ! temperature at level 'ilvb'
    b_0         ,& ! lapse rate  at level 'ilvb'
    p_0            ! pressure    at level 'ilvb'

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------
!- Begin SUBROUTINE std_atm_p2z_fast
!-------------------------------------------------------------------------------

  z_out   = z_miss
  bi_hom  = - r_d / r_g

  DO klev = 1 , nlvp
    ilev = isrtlvp(klev)
    IF (p_in(ilev) > p_max) THEN
      ilvb (klev) = 0
    ELSEIF (p_in(ilev) > std_atm(4,2)) THEN
      ilvb (klev) = 1
    ELSEIF (p_in(ilev) > std_atm(4,3)) THEN
      ilvb (klev) = 2
    ELSEIF (p_in(ilev) > std_atm(4,4)) THEN
      ilvb (klev) = 3
    ELSEIF (p_in(ilev) > std_atm(4,5)) THEN
      ilvb (klev) = 4
    ELSE  
      ilvb (klev) = 0
    ENDIF
  ENDDO
  DO klev = 1 , nlvp
    IF (ilvb(klev) > 0) THEN
      z_0 = std_atm(1,ilvb(klev))
      t_0 = std_atm(2,ilvb(klev))
      b_0 = std_atm(3,ilvb(klev))
      p_0 = std_atm(4,ilvb(klev))
      ilev = isrtlvp(klev)
      IF (ABS( b_0 ) < epsy) THEN   ! isothermal atmosphere
        z_out (ilev) = z_0 + (bi_hom *t_0) * LOG( p_in(ilev) /p_0 )
      ELSE                          ! polytropic atmosphere
        z_out (ilev) = z_0 - t_0 /b_0 *(c1 - (p_in(ilev) /p_0)**(b_0*bi_hom))
      ENDIF
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine std_atm_p2z_fast
!-------------------------------------------------------------------------------
END SUBROUTINE std_atm_p2z_fast


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for standard atmosphere information
!-------------------------------------------------------------------------------

SUBROUTINE obs_solar_zenith_angle ( nrep, rmdich, kdate, khrmn, replat, replon &
                                  , repsolz )

!-------------------------------------------------------------------------------
! Description:
!   Compute the solar zenith angle for an array of reports (more generally:
!   for an array of 4-dim locations, consisting of date, time, latitude,
!   longitude).
!
! Method:
!   The solar zenith angle is computed using the same equations as in the
!   radiation parameterisation of the COSMO model.
!
! Initial release: Christoph Schraff, 10.02.11
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments
! --------------------

  INTEGER                  , INTENT (IN)    ::       &
    nrep             ,& ! number of reports
    kdate   (nrep)   ,& ! date of reports [yyyymmdd]
    khrmn   (nrep)      ! time of reports [hhmm]

  REAL    (KIND=wp)        , INTENT (IN)    ::       &
    replat  (nrep)   ,& ! latitude  of reports [deg]
    replon  (nrep)   ,& ! longitude of reports [deg]
    rmdich              ! check value for missing data (= -1.E30_wp)

  REAL    (KIND=wp)        , INTENT (INOUT) ::       &
    repsolz (nrep)      ! solar zenith angle of reports [deg]

! Local parameters:
! ----------------

! REAL    (KIND=wp)        , PARAMETER  :: &
!   c1  =  1.0_wp   ! standard real constant 1.0

  INTEGER                  , PARAMETER  :: &
    mdd_offs (13) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334   &
                                                                       , 365/)

! Local variables
! ---------------

  INTEGER                                ::       &
    irep             ,& ! loop  index over reports
    irps   , irpe    ,& ! start and end indices of loop over reports for which
                        !   the solar zenith angle has to be computed
    iyyyy            ,& ! year (incl. century) of observation
    imm              ,& ! month                of observation
    idd              ,& ! day                  of observation
    iddjul              ! Julian day (i.e. days since the beginning of the year)

  REAL    (KIND=wp)                      ::       &
    pi               ,& ! circle constant
    degrad           ,& ! factor for transforming degree to rad
    zhour            ,& ! hour of observation
    ztwo   , ztho    ,& ! factors in time equation (day (angle) of solar year)
    zdtzgl , zdek    ,& ! factors in time equation (day (angle) of solar year)
    zlonoon          ,& ! longitude where it is exactly noon at the obs time
    ztimrad          ,& ! longitudinal radians related to the time
    zdeksin, zdekcos ,& ! factors from time equation
    zsinlat, zcoslat ,& ! factors depending on latitude
    zcosthi, zcoszen ,& ! COS( solar zenith angle )  / above horizon
    zsolzen             ! solar zenith angle    [rad]

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE obs_solar_zenith_angle
!-------------------------------------------------------------------------------

  pi       =   4.0_wp * ATAN( c1 )
  degrad   =   pi / 180.0_wp

  !   limit loop to range of reports with missing solar zenith angle
  irps = nrep
  irpe = 1
  DO irep = 1 , nrep
    IF (repsolz(irep) <= rmdich) THEN
      irps = MIN( irps , irep )
      irpe = MAX( irpe , irep )
    ENDIF
  ENDDO

  DO irep = irps , irpe

    imm   = MOD( kdate(irep) , 10000 ) / 100
    idd   = MOD( kdate(irep) , 100 )
    iyyyy =      kdate(irep)           / 10000
    zhour =   REAL(      khrmn(irep) / 100  , wp)                               &
            + REAL( MOD( khrmn(irep) , 100 ), wp) / 60._wp
    ! get Julian day (i.e. days since the beginning of the years)
    iddjul  =  mdd_offs(imm) + idd
    IF (imm > 2)  iddjul  =  iddjul + ileap(iyyyy)

    ! use time equation as in src_radiation
    ztwo    = 0.681_wp + 0.2422_wp*(iyyyy-1949)-(iyyyy-1949)/4
    ztho    = 2.0_wp*pi*( REAL(iddjul, wp) -1.0_wp + ztwo )/365.2422_wp
    zdtzgl  = 0.000075_wp + 0.001868_wp*COS(       ztho) - 0.032077_wp*SIN(       ztho) &
                          - 0.014615_wp*COS(2.0_wp*ztho) - 0.040849_wp*SIN(2.0_wp*ztho)
    zdek    = 0.006918_wp - 0.399912_wp*COS(       ztho) + 0.070257_wp*SIN(       ztho) &
                          - 0.006758_wp*COS(2.0_wp*ztho) + 0.000907_wp*SIN(2.0_wp*ztho) &
                          - 0.002697_wp*COS(3.0_wp*ztho) + 0.001480_wp*SIN(3.0_wp*ztho)

    ! longitude-dependent part
      ! get the longitude where it is exactly noon at the obs time
    zlonoon = pi*(zhour-12._wp)/12._wp + zdtzgl
    ztimrad = zlonoon + degrad* replon(irep)

!     WRITE(0,'("SOLAR ",A ,6F7.2,F7.4,2I11,2I4)' )                            &
!           ystid (irep), zhour, replat(irep), replon(irep)                    &
!         , zdtzgl, zdek, ztimrad, degrad, kdate(irep)                         &
!         , khrmn(irep), ntotmlo, ntotml
!     WRITE(0,'("SOLAR1 ",I3,1X,A ,5I6,F7.2,2I11,2F8.4)' )  my_cart_id,        &
!           ystid (irep), ntotmlo, ntotml, irep, iddjul, iyyyy, zhour          &
!                      , kdate(irep), khrmn(irep), replat(irep), replon(irep)
!     WRITE(0,'("SOLAR2 ",I3,1X,A ,4(F9.2,1X)         )' )  my_cart_id,        &
!           ystid (irep), zdtzgl, zdek, ztimrad, degrad

    ! latitude-dependent part
    zdeksin = SIN( zdek )
    zdekcos = COS( zdek )
    zsinlat = SIN( degrad* replat(irep) )
    zcoslat = SQRT( c1 - zsinlat*zsinlat )

    ! get COS( solar zenith angle )  (above horizon)
    zcosthi = zdeksin * zsinlat + zdekcos * zcoslat * COS( ztimrad )
    zcoszen = MIN( MAX( zcosthi, -c1 ) , c1 )

    ! get solar zenith angle
    zsolzen = ACOS( zcoszen )
    repsolz (irep)  =  zsolzen / degrad
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_solar_zenith_angle
!-------------------------------------------------------------------------------
END SUBROUTINE obs_solar_zenith_angle


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for getting observed relat. humidity
!-------------------------------------------------------------------------------

SUBROUTINE obs_td2rh ( madj_hum, kdim, ztt , ztd                               &
                     , rmdi, b1, b2w, b2i, b3, b4w, b4i, t0_melt , zrhw, zrh )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" converts an array of
!   dew point observations into model compatible relative humidity.
!
! Method:
!   If the model 'knows about' (i.e. has as prognostic variables) cloud ice
!   in addition to cloud water (like in the observed nature) then observed
!   values of dew point or vapour pressure can be converted into relative
!   humidity using the standard conversion formula (which means that all values
!   are 'over water').
!   In contrast, if cloud ice is not a prognostic variable, then saturation
!   over water is required in the model in order to obtain grid-scale cloud.
!   In nature, however, saturation over ice is in principle sufficient to
!   maintain cloud (at -40 C, this corresponds to about 60 % relative humdity
!   (defined over water)). Therefore, observed saturation over ice should be
!   equivalent to saturation over water in the model. Thus, 'observed'
!   relative humidity is defined to be the ratio of observed vapour pressure
!   and saturation vapour pressure over ice.
!
!   Note: All the required input variables are contained in the parameter list,
!
! Initial release: Christoph Schraff, 20.12.07
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)     ::  &
    madj_hum     ,& ! = 0  : cloud ice exists as a (model) state variable,
                    !        and the observed dew point / vapour pressure can be
                    !        converted into relative humidity using the standard
                    !        formula
                    ! = 1  : cloud ice is not a state variable; the 'observed'
                    !        relative humidity has to be adjusted by defining it
                    !        as the ratio of observed vapour pressure and
                    !        saturation vapour pressure over ice (instead of
                    !        over water) in order to make the observation
                    !        compatible with the model (see 'Method' above)
    kdim            ! length of array

  REAL    (KIND=wp   )     , INTENT (IN)     ::  &
    ztt   (kdim) ,& ! temperature
    rmdi         ,& ! commonly used missing data indicator (large negative)
    b1           ,& ! variables for computing the saturation vapour pressure
    b2w          ,& ! over water (w) and ice (i)
    b2i          ,& !               -- " --
    b3           ,& !               -- " --
    b4w          ,& !               -- " --
    b4i          ,& !               -- " --
    t0_melt         ! melting temperature of ice

  REAL    (KIND=wp   )     , INTENT (INOUT)  ::  &
    ztd   (kdim)    ! dew point  (input: original ; output: model-compatible)

  REAL    (KIND=wp   )     , INTENT (OUT)    ::  &
    zrhw  (kdim) ,& ! relative humidity over water
    zrh   (kdim)    ! model compatible relative humidity

! Local parameters:
! ----------------

! REAL    (KIND=wp)        , PARAMETER  :: &
!   epsy   =  1.E-8_wp     ! commonly used very small value > 0

! Local scalars:
! -------------

  INTEGER                  ::  &
    kk              ! loop index

  REAL    (KIND=wp   )     ::  &
    eobs         ,& ! observed vapour pressure
    eswo         ,& ! saturation vapour pressure over water
    eseo         ,& ! saturation vapour pressure over ice
    rmdich          ! check value for missing data

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_td2rh
!-------------------------------------------------------------------------------

  rmdich  =  ABS( 0.5_wp* rmdi )

  DO kk = 1 , kdim
    IF (      (ABS( ztd(kk) ) < rmdich) .AND. (ztd(kk) > b4w+epsy)             &
        .AND. (ABS( ztt(kk) ) < rmdich) .AND. (ztt(kk) > b4w+epsy)) THEN

!         eobs : observed vapour pressure (measurement is inherently over ice)
!                 (but td as stored in data file is over water (WMO convention),
!                 which is why Magnus formula over water instead of ice is used)
      eobs       =  fpvsw( ztd(kk), b1, b2w, b3, b4w )

!         eswo : saturation vapour pressure over water at observed temperature
      eswo       =  fpvsw( ztt(kk), b1, b2w, b3, b4w )
      zrhw (kk)  =  eobs / eswo
      zrh  (kk)  =  zrhw(kk)

! convert to model compatible relative humidity
! ---------------------------------------------
      IF ((ztt(kk) <  t0_melt) .AND. (madj_hum == 1)) THEN
!           eseo : saturation vapour pressure over ice at observed temperature
        eseo      =  fpvsi( ztt(kk), b1, b2i, b3, b4i )

!           zrh  : observed relative humidity over ice   (per def. <= 100%)
        zrh (kk)  =  eobs / eseo

        IF (LOG(eobs/b1) < b2w-epsy)  ztd (kk)  =  ftd( eobs, b1, b2w, b3, b4w )
      ENDIF
    ELSE
      zrhw (kk)  =  rmdi
      zrh  (kk)  =  rmdi
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_td2rh
!-------------------------------------------------------------------------------
END SUBROUTINE obs_td2rh


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for getting observed relat. humidity
!-------------------------------------------------------------------------------

SUBROUTINE obs_qx2rh ( madj_hum, kdim, ztt, zpp, zqxw , zqvw                   &
                     , rmdi, b1, b2w, b2i, b3, b4w, b4i, rdv, t0_melt          &
                     , zqv, zrhw, zrh )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" converts an array of
!   mixing ratio observations into model compatible relative humidity.
!
! Method:
!   If the model 'knows about' (i.e. has as prognostic variables) cloud ice
!   in addition to cloud water (like in the observed nature) then observed
!   values of mixing ratio or vapour pressure can be converted into relative
!   humidity using the standard conversion formula (which means that all values
!   are 'over water').
!   In contrast, if cloud ice is not a prognostic variable, then saturation
!   over water is required in the model in order to obtain grid-scale cloud.
!   In nature, however, saturation over ice is in principle sufficient to
!   maintain cloud (at -40 C, this corresponds to about 60 % relative humdity
!   (defined over water)). Therefore, observed saturation over ice should be
!   equivalent to saturation over water in the model. Thus, 'observed'
!   relative humidity is defined to be the ratio of observed vapour pressure
!   and saturation vapour pressure over ice.
!
!   Note: All the required input variables are contained in the parameter list,
!
! Initial release: Christoph Schraff, 20.12.07
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)     ::  &
    madj_hum     ,& ! = 0  : cloud ice exists as a (model) state variable,
                    !        and the observed dew point / vapour pressure can be
                    !        converted into relative humidity using the standard
                    !        formula
                    ! = 1  : cloud ice is not a state variable; the 'observed'
                    !        relative humidity has to be adjusted by defining it
                    !        as the ratio of observed vapour pressure and
                    !        saturation vapour pressure over ice (instead of
                    !        over water) in order to make the observation
                    !        compatible with the model (see 'Method' above)
    kdim            ! length of array

  REAL    (KIND=wp   )     , INTENT (IN)     ::  &
    ztt   (kdim) ,& ! temperature
    zpp   (kdim) ,& ! pressure
    zqxw  (kdim) ,& ! mixing ratio (original value, over water)
    rmdi         ,& ! commonly used missing data indicator (large negative)
    b1           ,& ! variables for computing the saturation vapour pressure
    b2w          ,& ! over water (w) and ice (i)
    b2i          ,& !               -- " --
    b3           ,& !               -- " --
    b4w          ,& !               -- " --
    b4i          ,& !               -- " --
    rdv          ,& ! r_d / r_v
    t0_melt         ! melting temperature of ice

  REAL    (KIND=wp   )     , INTENT (INOUT)  ::  &
    zqvw  (kdim)    ! specific water vapour (over water)
                       !   (input value, but if = missing value then it is
                       !    computed from zqxw , i.e. then is it output)

  REAL    (KIND=wp   )     , INTENT (OUT)    ::  &
    zrhw  (kdim) ,& ! relative humidity over water
    zrh   (kdim) ,& ! model compatible relative humidity
    zqv   (kdim)    ! model compatible specific water vapour

! Local parameters:
! ----------------

! REAL    (KIND=wp)        , PARAMETER  :: &
!   c1     = 1.0_wp     ,& ! standard real constant 1.0
!   epsy   =  1.E-8_wp     ! commonly used very small value > 0

! Local scalars:
! -------------

  INTEGER                  ::  &
    kk              ! loop index

  REAL    (KIND=wp   )     ::  &
    eswo         ,& ! saturation vapour pressure over water
    eseo         ,& ! saturation vapour pressure over ice
    rmdich          ! check value for missing data

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_qx2rh
!-------------------------------------------------------------------------------

  rmdich  =  ABS( 0.5_wp* rmdi )

  DO kk = 1 , kdim

! compute specific humidity from mixing ratio (if specific humidity missing)
    IF ((ABS( zqvw(kk) ) > rmdich) .AND. (ABS( zqxw(kk) ) < rmdich)) THEN
      zqvw (kk)  =  zqxw(kk) / (c1 + zqxw(kk))
    ENDIF

! compute relative humidity from specific humidity
    IF (      (ABS( zqvw(kk) ) < rmdich) .AND. (ztt(kk) > b4w+epsy)            &
        .AND. (ABS( zpp (kk) ) < rmdich) .AND. (ABS( ztt(kk) ) < rmdich)) THEN

!         eswo : saturation vapour pressure over water at observed temperature
      eswo       =  fpvsw( ztt(kk), b1, b2w, b3, b4w )

!                relative humidity is over water (WMO convention)
      zrhw (kk)  =  fq2pv( zqvw(kk) , zpp(kk) , rdv ) / eswo
      zrh  (kk)  =  zrhw(kk)

! convert to model compatible relative humidity
! ---------------------------------------------
      IF ((ztt(kk) <  t0_melt) .AND. (madj_hum == 1)) THEN
!           eseo : saturation vapour pressure over ice at observed temperature
        eseo      =  fpvsi( ztt(kk), b1, b2i, b3, b4i )

!           zrh  : observed relative humidity over ice   (per def. <= 100%)
!                  (observed vapour pressure eobs = zrh *eseo = zrhw *eswo)
        zrh (kk)  =  zrhw(kk) * eswo / eseo
      ENDIF

      zqv  (kk)  =  zqvw(kk)
      IF ((c1-rdv)* zrh(kk) *eswo < zpp(kk) -epsy)                             &
        zqv  (kk)  =  fpv2q( zrh(kk) *eswo , zpp(kk) , rdv )
    ELSE
      zrhw (kk)  =  rmdi
      zrh  (kk)  =  rmdi
      zqv  (kk)  =  rmdi
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_qx2rh
!-------------------------------------------------------------------------------
END SUBROUTINE obs_qx2rh


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for getting observed relat. humidity
!-------------------------------------------------------------------------------

SUBROUTINE obs_rhw2rh ( madj_hum, kdim, ztt, zrhw                              &
                      , rmdi, b1, b2w, b2i, b3, b4w, b4i, t0_melt , zrh )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" converts an array of
!   relative humidity observations into model compatible relative humidity.
!
! Method:
!   If the model 'knows about' (i.e. has as prognostic variables) cloud ice
!   in addition to cloud water (like in the observed nature) then observed
!   values are already model compatible.
!   In contrast, if cloud ice is not a prognostic variable, then saturation
!   over water is required in the model in order to obtain grid-scale cloud.
!   In nature, however, saturation over ice is in principle sufficient to
!   maintain cloud (at -40 C, this corresponds to about 60 % relative humdity
!   (defined over water)). Therefore, observed saturation over ice should be
!   equivalent to saturation over water in the model. Thus, 'observed'
!   relative humidity is defined to be the ratio of observed vapour pressure
!   and saturation vapour pressure over ice.
!
!   Note: All the required input variables are contained in the parameter list,
!
! Initial release: Christoph Schraff, 20.12.07
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)     ::  &
    madj_hum     ,& ! = 0  : cloud ice exists as a (model) state variable,
                    !        and the observed dew point / vapour pressure can be
                    !        converted into relative humidity using the standard
                    !        formula
                    ! = 1  : cloud ice is not a state variable; the 'observed'
                    !        relative humidity has to be adjusted by defining it
                    !        as the ratio of observed vapour pressure and
                    !        saturation vapour pressure over ice (instead of
                    !        over water) in order to make the observation
                    !        compatible with the model (see 'Method' above)
    kdim            ! length of array

  REAL    (KIND=wp   )     , INTENT (IN)     ::  &
    ztt   (kdim) ,& ! temperature
    zrhw  (kdim) ,& ! original         relative humidity (over water)
    rmdi         ,& ! commonly used missing data indicator (large negative)
    b1           ,& ! variables for computing the saturation vapour pressure
    b2w          ,& ! over water (w) and ice (i)
    b2i          ,& !               -- " --
    b3           ,& !               -- " --
    b4w          ,& !               -- " --
    b4i          ,& !               -- " --
    t0_melt         ! melting temperature of ice

  REAL    (KIND=wp   )     , INTENT (OUT)    ::  &
    zrh   (kdim)    ! model-compatible relative humidity

! Local parameters:
! ----------------

! Local scalars:
! -------------

  INTEGER                  ::  &
    kk              ! loop index

  REAL    (KIND=wp   )     ::  &
    eobs         ,& ! observed vapour pressure
    eswo         ,& ! saturation vapour pressure over water
    eseo         ,& ! saturation vapour pressure over ice
    rmdich          ! check value for missing data

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_rhw2rh
!-------------------------------------------------------------------------------

  rmdich  =  ABS( 0.5_wp* rmdi )

  DO kk = 1 , kdim

! convert to model compatible relative humidity
    IF (      (madj_hum == 1)            .AND. (ztt(kk) < t0_melt-epsy)        &
        .AND. (ABS( zrhw(kk) ) < rmdich) .AND. (ztt(kk) > b4w+epsy)) THEN
!         eswo : saturation vapour pressure over water at observed temperature
      eswo       =  fpvsw( ztt(kk), b1, b2w, b3, b4w )
!         eobs : observed vapour pressure
!                 (reported relative humidity is over water by WMO convention)
      eobs       =  zrhw(kk) * eswo
!         eseo : saturation vapour pressure over ice at observed temperature
      eseo       =  fpvsi( ztt(kk), b1, b2i, b3, b4i )
!         zrh  : observed relative humidity over ice   (per def. <= 100%)
      zrh  (kk)  =  eobs / eseo

! if observed relative humidity is not model-compatible, and temperature
! is not reported, then observed relative humidity is discarded
    ELSEIF ((madj_hum == 1) .AND. (.NOT.(      (ztt(kk) >= t0_melt-epsy)       &
                                         .AND. (ABS( ztt(kk) ) < rmdich)))) THEN
      zrh  (kk)  =  rmdi

! otherwise, 'zrh' remains unchanged
    ELSE
      zrh  (kk)  =  zrhw(kk)
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_rhw2rh
!-------------------------------------------------------------------------------
END SUBROUTINE obs_rhw2rh


!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for rel. hum. to dewpoint conversion
!-------------------------------------------------------------------------------

SUBROUTINE obs_rh2td ( kdim, zrh, ztt, rmdi, b1, b2w, b3, b4w , ztd )

!-------------------------------------------------------------------------------
! Description:
!   This procedure converts an array of relative humidity data into dew point
!   temperature.
!
! Method:
!   The conversion is performed using standard formula (i.e. both relative
!   humidity and dewpoint temperature are defined 'over water' and treated
!   correspondingly).
!   No 'model compatible' adjustments are made, even if cloud ice is not a
!   prognostic model variable.
!
!   Note: All the required input variables are contained in the parameter list.
!
! Initial release: Christoph Schraff, 24.05.15
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)     ::  &
    kdim            ! length of array

  REAL    (KIND=wp   )     , INTENT (IN)     ::  &
    zrh   (kdim) ,& ! model compatible relative humidity
    ztt   (kdim) ,& ! temperature
    rmdi         ,& ! commonly used missing data indicator (large negative)
    b1           ,& ! variables for computing the saturation vapour pressure
    b2w          ,& ! over water (w) and ice (i)
    b3           ,& !               -- " --
    b4w             !               -- " --

  REAL    (KIND=wp   )     , INTENT (OUT)    ::  &
    ztd   (kdim)    ! dew point  (input: original ; output: model-compatible)

! Local parameters:
! ----------------

! REAL    (KIND=wp)        , PARAMETER  :: &
!   epsy   =  1.E-8_wp     ! commonly used very small value > 0

! Local scalars:
! -------------

  INTEGER                  ::  &
    kk              ! loop index

  REAL    (KIND=wp   )     ::  &
    eobs         ,& ! observed vapour pressure
    eswo         ,& ! saturation vapour pressure over water
    rmdich          ! check value for missing data

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_rh2td
!-------------------------------------------------------------------------------

  rmdich  =  ABS( 0.5_wp* rmdi )

  DO kk = 1 , kdim
    ztd (kk)  =  rmdi
    IF (      (ABS( zrh(kk) ) <= c1+epsy) .AND. (zrh(kk) >= -epsy)             &
        .AND. (ABS( ztt(kk) ) <  rmdich ) .AND. (ztt(kk) > b4w+epsy)) THEN
      eswo       =  fpvsw( ztt(kk), b1, b2w, b3, b4w )
      eobs       =  eswo * MAX( zrh(kk) , 0.0001_wp )
      !   (td, as e.g. stored in data files, is over water (WMO convention),
      !    thus Magnus formula over water instead of ice is used)
      IF (LOG(eobs/b1) < b2w-epsy)  ztd (kk)  =  ftd( eobs, b1, b2w, b3, b4w )
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_rh2td
!-------------------------------------------------------------------------------
END SUBROUTINE obs_rh2td


!===============================================================================

!===============================================================================
!+ Module procedure in "src_obs_cdfin_util" for getting interpolation levels
!-------------------------------------------------------------------------------

SUBROUTINE obs_find_level ( numlev, zloplev, zlopt , ilvlow, fiplow )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_util" determines the level
!   index and interpolation factor which are required to interpolate some
!   quantity given at a pre-specified set of vertical levels (e.g. standard
!   levels, given in log( pressure ) units) to a given pressure level
!   (specified in log( pressure ) units).
!   The interpolation factor allows for linear interpolation in ln(p).
!   (For instance, this can be used to interpolate pre-specified observation
!    errors given at standard levels to the pressure level of a certain
!    observation.)
!
! Method:
!   The interpolation factor allows for linear interpolation in ln(p).
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)     ::  &
    numlev              ! number of input levels (at which a quantity to be
                        !                         interpolated is pre-specified)

  REAL    (KIND=wp   )     , INTENT (IN)     ::  &
    zloplev (numlev) ,& ! log(pressure) of the input levels
    zlopt               ! log(pressure) of the target level

  INTEGER                  , INTENT (OUT)    ::  &
    ilvlow              ! nearest input level below 'zlopt'

  REAL    (KIND=wp   )     , INTENT (OUT)    ::  &
    fiplow              ! interpolation weight factor for input level 'ilvlow'

! Local parameters:
! ----------------
! Local scalars: None
! -------------
! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_find_level
!-------------------------------------------------------------------------------

! keep interpolation quantity constant below lowest input level
  IF (zlopt >= zloplev(1)-epsy) THEN
    ilvlow = 1
    fiplow = 1.0_wp
! keep inperpolation quantity constant above highest input level
  ELSEIF (zlopt <= zloplev(numlev)+epsy) THEN
    ilvlow = numlev - 1
    fiplow = 0.0_wp

! find nearest input level below target pressure
  ELSE
    ilvlow = 1
    DO WHILE ((zloplev(ilvlow+1) >= zlopt-epsy) .AND. (ilvlow < numlev-1))
      ilvlow = ilvlow + 1
    ENDDO
    IF (ABS( zlopt-zloplev(ilvlow+1) ) <  epsy) THEN
      fiplow = 0.0_wp
    ELSE
      fiplow =   (zlopt           - zloplev(ilvlow+1))                         &
               / (zloplev(ilvlow) - zloplev(ilvlow+1))
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_find_level
!-------------------------------------------------------------------------------
END SUBROUTINE obs_find_level

!===============================================================================

  SUBROUTINE check_wind_method ( iident, sfreq , istret )
  !-----------------------------------------------------------------------------
  ! This routine must correspond to the one in 'mo_amv.f90' of DACE.
  !-----------------------------------------------------------------------------
    integer   ,intent(inout) :: istret ! sat. derived wind computation method
    integer   ,intent(in)    :: iident ! satellite identifier
    real(wp)  ,intent(in)    :: sfreq  ! frequency of sat data providing AMV
    !-----------------------------------------------------------------------
    ! Code Figure 'MSDWCM' (satellite derived wind computation meth.)
    ! 0 Not used
    ! 1 Wind derived from cloud motion observed in the infrared channel
    ! 2 Wind derived from cloud motion observed in the visible channel
    ! 3 Wind derived from motion observed in the water vapour channel
    !   (only the above are used up to now)
    ! 4 Wind derived from motion observed in a combination of spectral chans
    ! 5 Wind derived from motion obs. in the water vapour ch. in clear air
    ! 6 Wind derived from motion observed in the ozon channel
    ! 7 Wind derived from cloud motion observed in water vapour channel
    !   (cloud or clear air not specified)
    !
    ! Add satellite/frequency digit (n*100) to 'istret'
    !-----------------------------------------------------------------------
    real(wp), parameter :: epsilon = 0.00001e+15_wp
    integer             :: digits

    if (istret >= 1) then
      digits = 0
      select case (iident)
      case (55:57,70)   ! Meteosat-8, Meteosat-9, Meteosat-10, Meteosat-11
        if     (abs (sfreq - 0.4721142137D15) < epsilon) then
           digits = 100                 ! VIS1: 0.6 um
        elseif (abs (sfreq - 0.399723D15)     < epsilon) then
           digits = 300                 ! VIS2 (high resolution): 0.75 um
        elseif (abs (sfreq - 0.3701142098D15) < epsilon) then
           digits = 200                 ! VIS3: 0.8 um
        elseif (abs (sfreq - 0.4796679874D14) < epsilon) then
           digits = 100                 ! WV1 : 6.2 um
        elseif (abs (sfreq - 0.4078810064D14) < epsilon) then
           digits = 200                 ! WV2 : 7.3 um
        elseif (abs (sfreq - 0.344589D14)     < epsilon) then
           digits = 100                 ! IR1 : 8.7 um
        elseif (abs (sfreq - 0.310344D14)     < epsilon) then
           digits = 200                 ! IR2 : 9.7 um
        elseif (abs (sfreq - 0.2775860013D14) < epsilon) then
           digits = 300                 ! IR3 : 10.8 um
        else
!         WRITE(0,*) '%E SATOBFREQ: channel frequency for', iident,&
!                    ' has not been recognized', sfreq
        endif
      case (254:259) ! GOES-10, GOES-11, GOES-12, GOES-13, GOES-14, GOES-15
        if     (abs (sfreq - 0.461538D15) < epsilon) then
           digits = 100                 ! VIS1 (0.65um)
        elseif (abs (sfreq - 0.405405D14) < epsilon) then
           digits = 100                 ! WV1 (7.4um)
        elseif (abs (sfreq - 0.428571D14) < epsilon) then
           digits = 200                 ! WV2 (7.0um)
        elseif (sfreq > 0.441176D14-epsilon .and. sfreq < 0.461219D14+epsilon) then
           digits = 300                 ! WV3 (6.5 till 6.8um)
        elseif (abs (sfreq - 0.280374D14) < epsilon) then
           digits = 100                 ! IR1 (10.7um)
        elseif (abs (sfreq - 0.768699D14) < epsilon) then
           digits = 200                 ! IR2 (3.9um)
        else
           WRITE(0,*) '%E SATOBFREQ: channel frequency for', iident,&
                      ' has not been recognized', sfreq
        endif
      case (270, 271) ! GOES-16, GOES-17
        if     (abs (sfreq - 0.468426D15) < epsilon) then
           digits = 100                 ! VIS1 (0.64um)
        elseif (abs (sfreq - 0.408437D14) < epsilon) then
           digits = 100                 ! WV1 (7.3um)
        elseif (abs (sfreq - 0.431356D14) < epsilon) then
           digits = 200                 ! WV2 (7.0um)
        elseif (abs (sfreq - 0.484317D14) < epsilon) then
           digits = 300                 ! WV3 (6.2um)
        elseif (abs (sfreq - 0.267672D14) < epsilon) then
           digits = 100                 ! IR1 (11.1um)
        elseif (abs (sfreq - 0.768699D14) < epsilon) then
           digits = 200                 ! IR2 (3.9um)
        else
           WRITE(0,*) '%E SATOBFREQ: channel frequency for', iident,&
                      ' has not been recognized', sfreq
        endif
      case (173, 174) ! Himawari-8/9
        if     (abs (sfreq - 0.4684257D15) < epsilon) then
           digits = 100                 ! VIS1: 0.645 um
        elseif (abs (sfreq - 0.483536D14) < epsilon) then
           digits = 100                 ! WV1 : 6.25 um
        elseif (abs (sfreq - 0.434481D14) < epsilon) then
           digits = 200                 ! WV2 : 6.95 um
        elseif (abs (sfreq - 0.410674D14) < epsilon) then
           digits = 300                 ! WV3 : 7.35 um
        elseif (abs (sfreq - 0.288261D14) < epsilon) then
           digits = 100                 ! IR1 (10.45um)
        else
           WRITE(0,*) '%E SATOBFREQ: channel frequency for', iident,&
                      ' has not been recognized', sfreq
        endif
      case default
!       WRITE(0,*) '%E SATOBFREQ: channel frequency for', iident,&
!                  ' has not been recognized', sfreq
      end select
      istret = istret + digits
    else
!     WRITE(0,*) '%E SATOBFREQ: invalid value for "stret"', istret
    endif
  END SUBROUTINE check_wind_method


!===============================================================================
!+ Module function in "src_obs_cdfin_util" for model pressure at a given height
!-------------------------------------------------------------------------------

FUNCTION f_z2p ( zzz , ke, hhl_col, p_col, t_ll, z_g, z_r_d, rmdi, lextpol )

!-------------------------------------------------------------------------------
! Description:
!   Get the model pressure for a specified height (in order to assign PILOT or
!   SATOB reports to a pressure level, if pressure itself is not reported).
!
! Remarks:
!   - if height zzz is below model surface then f_z2p is missing value 'rmdi'
!   - the parameter list contains all the required input variables.
!
! Current Code Owner:  Christoph Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
!
! Modules used:         These are declared in the module declaration section
!
!-------------------------------------------------------------------------------

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Parameter list:
! --------------
  INTEGER                  , INTENT (IN)   :: &
    ke                ! number of main model levels

  REAL    (KIND=wp)        , INTENT (IN)   :: &
    zzz            ,& ! height for which pressure is sought
    hhl_col (ke+1) ,& ! model column of height   (half levels)
    p_col   (ke)   ,& ! model column of pressure (main levels)
    t_ll           ,& ! temperature at the lowest model (main) level
    z_g            ,& ! acceleration due to gravity
    z_r_d          ,& ! gas constant for dry air
    rmdi              ! value for missing data

  LOGICAL                  , INTENT (IN)   :: &
    lextpol           ! extrapolate beyond model orography

! Local variables:
! ---------------

  INTEGER                  ::      &
    kk                ! loop index

  REAL (KIND=wp)          ::      &
    hhke           ,& ! height of lowest full model level
    zzo    , zzu   ,& ! height of upper / lower model level
    zfrac          ,& ! weight
    f_z2p             ! return value of function: pressure at height 'zzz'

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function f_z2p
!-------------------------------------------------------------------------------

  f_z2p =  rmdi
  hhke  =  0.5_wp* (hhl_col(ke+1) + hhl_col(ke))
  zzo   =  hhke

  IF (zzz >= hhke) THEN
    !  find appropriate model layer for the specified height
    DO kk = ke, 2, -1
      zzu = zzo
      zzo = 0.5_wp* (hhl_col(kk) + hhl_col(kk-1))
      IF ((zzu <= zzz) .AND. (zzo > zzz)) THEN

        !  interpolate ln(p) linear in z
        zfrac  =  ( zzz - zzu ) / ( zzo - zzu )
        f_z2p  =  EXP(  (1.0_wp-zfrac)* LOG( p_col(kk)   )                 &
                      +         zfrac * LOG( p_col(kk-1) ) )
        EXIT
      ENDIF
    ENDDO
  ELSEIF ((zzz >= hhl_col(ke+1)) .OR. (lextpol)) THEN
    !  reduce pressure from lowest model level to the specified height
    f_z2p  =  EXP( LOG(p_col(ke)) + z_g/(z_r_d* t_ll) *(hhke-zzz) )
  ENDIF

!-------------------------------------------------------------------------------
! End of module function f_z2p
!-------------------------------------------------------------------------------

END FUNCTION f_z2p

!===============================================================================
!+ Module function in "src_obs_cdfin_util" for approx. model layer thickness
!-------------------------------------------------------------------------------

FUNCTION f_p2dp ( zpp , ke, p_col, dp_col )

!-------------------------------------------------------------------------------
! Description:
!   Get the approximate model layer thickness for a specified pressure level
!   (used to adjust the scale of the vertical correlation function for
!    vertically dense upper-air (aircraft) observations in 'obs_air_correl').
!
! Remark:
!   The parameter list contains all the required input variables.
!
! Current Code Owner:  Christoph Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
!
! Modules used:         These are declared in the module declaration section
!
!-------------------------------------------------------------------------------

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Parameter list:
! --------------
  INTEGER                  , INTENT (IN)   :: &
    ke                ! number of main model levels

  REAL    (KIND=wp)        , INTENT (IN)   :: &
    zpp            ,& ! pressure for which model layer thickness is sought
    p_col   (ke)   ,& ! model column of pressure (main levels)
    dp_col  (ke)      ! model column of pressure (main levels)

! Local variables:
! ---------------

  INTEGER                  ::      &
    kdp            ,& ! index of model layer encompassing pressure zpp
    kk                ! loop index

  REAL (KIND=wp)          ::      &
    f_p2dp            ! return value of function: model layer thickness at 'zpp'

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function f_p2dp
!-------------------------------------------------------------------------------

  kdp = ke - 1

  DO kk = ke-2, 2, -1
    IF (p_col(kk) + 0.5_wp* dp_col(kk) > zpp)  kdp = kk
  ENDDO

  f_p2dp  =  dp_col(kdp)

!-------------------------------------------------------------------------------
! End of module function f_p2dp
!-------------------------------------------------------------------------------

END FUNCTION f_p2dp


!===============================================================================
!+ Module function in "src_obs_cdfin_util" for determining fdbk report status
!-------------------------------------------------------------------------------

ELEMENTAL INTEGER FUNCTION ncfeedobs_status ( modr_stat, kflags )

!-------------------------------------------------------------------------------
! Description:
!   Convert report status from ODR format (word 'nhpass') to ncfeedobs format
!   (word 'r_flags').
!
! Current Code Owner:  Christoph Schraff
!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Parameter list:
! --------------
  INTEGER                  , INTENT (IN)  :: &
    modr_stat  ,& ! status of ODR status word
    kflags        ! quality bit flag (r_flags)

! Local variables:
! ---------------
  INTEGER                                 :: &
    kflags_rp  ,& ! actual quality bit flag, which contains only the flags
                  !   related to status 'passive' or 'rejected'
    mpasflags     ! bit pattern of all flags related to status 'passive'
!   mpasflag      ! part of kflags which relate to status 'passive'
!
!------------ End of header ----------------------------------------------------

  IF (modr_stat == 0) THEN
    ncfeedobs_status = ST_ACTIVE
  ELSEIF (modr_stat == 1) THEN
    ncfeedobs_status = ST_MERGED
  ELSE
    !   at least 1 flag must be set in kflags_rp (but not flags 'NONE', 'MERGE')
    kflags_rp = ireplace1 ( kflags   , FL_MERGE  , 0, 0 )
    !   bit pattern of flags which lead to status 'passive' (not rejected)
    mpasflags = ireplace1 ( 0        , FL_OBSTYPE, 1, 0 )
    mpasflags = ireplace1 ( mpasflags, FL_THIN   , 1, 0 )
    IF (IAND( kflags_rp, mpasflags ) == 0) THEN
      !   no congruence between 'kflags_rp', 'mpasflags' --> no passive flag set
      ncfeedobs_status = ST_REJECTED
    ELSEIF (IAND( kflags_rp, NOT( mpasflags ) ) == 0) THEN
      !   no congruence between 'kflags_rp' and inverted 'mpasflags'
      !   --> no rejected flag set
      ncfeedobs_status = ST_PASSIVE
!   mpasflag = ireplace1 ( 0       , FL_OBSTYPE, kflags, FL_OBSTYPE )
!   mpasflag = ireplace1 ( mpasflag, FL_THIN   , kflags, FL_THIN )
!   IF (mpasflag == 0) THEN
!     ncfeedobs_status = ST_REJECTED
!   ELSEIF (mpasflag == ireplace1( kflags, FL_MERGE, 0, 0 )) THEN
!     ncfeedobs_status = ST_PASSIVE
    ELSE
      ncfeedobs_status = ST_PAS_REJ
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------

END FUNCTION ncfeedobs_status

!===============================================================================

! ELEMENTAL INTEGER FUNCTION ireplace   ( invar, ipos, iboc, irepl, ipsr )
  !-----------------------------------------------------------------------
! INTEGER                   , INTENT (IN)  ::  invar, ipos, iboc, irepl, ipsr
  !-----------------------------------------------------------------------------
  ! replaces 'iboc' bits starting at bit position 'ipos' of integer word 'invar'
  ! by the 'iboc' bits starting at bit position 'ipsr' from integer word 'irepl'
  !-----------------------------------------------------------------------------
  !
! ireplace = IOR( IAND( invar, NOT( ISHFT( nibits(iboc), ipos ) ) )            &
!               , ISHFT( IAND( ISHFT( irepl,-ipsr ), nibits(iboc) ), ipos ) )

! ireplace = IOR( IAND( invar, NOT( ISHFT( nibits(iboc), ipos ) ) )            &
!               , ISHFT( IBITS( irepl, ipsr, iboc ), ipos ) )
  !
! END FUNCTION ireplace

! Note: 'ireplace' could be replaced by intrinsic function 'mvbits' in F95 ff
! usage:   call mvbits (from, frompos, len, to, topos)
!
! Moves 'len' bits from positions 'frompos' through 'frompos+len-1' of 'from'
! to positions 'topos' through 'topos+len-1' of 'to'.
! The portion of argument 'to' not affected by the movement of bits is unchanged.
! The values of 'frompos+len-1' and 'topos+len-1' must be < bit_size(from).

!-------------------------------------------------------------------------------

! ELEMENTAL INTEGER FUNCTION insert  ( invar, inval, ibit )
  !--------------------------------------------------------
! INTEGER                   , INTENT (IN)  ::  invar, inval, ibit
  !----------------------------------------------------------------------------
  ! inserts bits set in 'inval' into integer word 'invar' at bit position 'ibp'
  !----------------------------------------------------------------------------
  !
! insert = IOR( invar, ISHFT( inval, ibit ) )
  !
! END FUNCTION insert

!-------------------------------------------------------------------------------

! ELEMENTAL INTEGER FUNCTION ibits  ( invar, ibp, iboc )
  !-----------------------------------------------------
! INTEGER                   , INTENT (IN)  ::  invar, ibp, iboc
  !---------------------------------------------------------------------------
  ! returns 'iboc' bits starting at bit position 'ibp' of integer word 'invar'
  !---------------------------------------------------------------------------
  !
! ibits = IAND( ISHFT( invar,-ibp ), nibits(iboc) )
  !
! END FUNCTION ibits

!-------------------------------------------------------------------------------

ELEMENTAL INTEGER FUNCTION ireplace1   ( invar, ipos, irepl, ipsr )
  !----------------------------------------------------------------
  INTEGER                   , INTENT (IN)  ::  invar, ipos, irepl, ipsr
  !-----------------------------------------------------------------------------
  ! replaces the bit at position 'ipos' of integer word 'invar'
  ! by the bit at position 'ipsr' from integer word 'irepl'
  !-----------------------------------------------------------------------------
  !
  ireplace1 = IOR( IAND( invar, NOT( ISHFT( 1, ipos ) ) )                      &
                 , ISHFT( IAND( ISHFT( irepl,-ipsr ), 1 ), ipos ) )
  !
END FUNCTION ireplace1

!-------------------------------------------------------------------------------

ELEMENTAL INTEGER FUNCTION ibit1  ( invar, ibp )
  !---------------------------------------------
  INTEGER                   , INTENT (IN)  ::  invar, ibp
  !---------------------------------------------------------------------------
  ! returns 1 bit at bit position 'ibp' of integer word 'invar'
  !---------------------------------------------------------------------------
  !
  ibit1 = IAND( ISHFT( invar,-ibp ), 1 )
  !
END FUNCTION ibit1

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION rmod  ( ztti, ziv )
  !----------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  ztti, ziv
  !---------------------------------------------------------------------------
  ! MOD function for positive REALS
  !---------------------------------------------------------------------------
  !
  rmod  =  ztti  -  ziv * REAL(INT( ztti/ziv + epsy ), wp) + epsy
  !
END FUNCTION rmod

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fpvsw  ( zt, b1, b2w, b3, b4w )
  !----------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zt, b1, b2w, b3, b4w
  !---------------------------------------------------------------------------
  ! Magnus formula for water:  input  'zt'   : temperature
  !                            output 'fpvsw': saturation water vapour pressure
  !        b1, b2w, b3, b4w :  constants of Magnus formula for water
  !---------------------------------------------------------------------------
  !
  fpvsw  =  b1 * EXP( b2w *(zt-b3) /(zt-b4w) )
  !
END FUNCTION fpvsw

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fpvsi  ( zt, b1, b2i, b3, b4i )
  !----------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zt, b1, b2i, b3, b4i
  !---------------------------------------------------------------------------
  ! Magnus formula for ice:  input  'zt'   : temperature
  !                          output 'fpvsi': saturation water vapour pressure
  !      b1, b2i, b3, b4i :  constants of Magnus formula for ice
  !---------------------------------------------------------------------------
  !
  fpvsi  =  b1 * EXP( b2i *(zt-b3) /(zt-b4i) )
  !
END FUNCTION fpvsi

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION ftd  ( zpv, b1, b2w, b3, b4w )
  !---------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zpv, b1, b2w, b3, b4w
  REAL    (KIND=wp)                        ::  zlogpv
  !---------------------------------------------------------------------------
  ! inverse of Magnus formula:  input  'zpv' :  water vapour pressure
  !                             output 'ftd' :  dewpoint temperature
  !         b1, b2w, b3, b4w :  constants of Magnus formula for water
  !---------------------------------------------------------------------------
  !
  zlogpv  =  LOG( zpv / b1 )
  ftd     =  (b3*b2w - zlogpv*b4w) / (b2w - zlogpv)
  !
END FUNCTION ftd

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fpv2q  ( zpv, zp, rdv )
  !--------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zpv, zp, rdv
  !---------------------------------------------------------------------------
  ! specific humidity from water vapour pressure 'zpv' and air pressure 'zp'
  !   (rdv = r_d / r_v
  !        = gas constant for dry air / gas constant for water vapour )
  !---------------------------------------------------------------------------
  !
  fpv2q  =  rdv * zpv / (zp - (c1-rdv)*zpv)
  !
END FUNCTION fpv2q

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fq2pv  ( zqv, zp, rdv )
  !--------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zqv, zp, rdv
  !---------------------------------------------------------------------------
  ! water vapour pressure (at T = saturation vapour press. at Td)
  ! from specific humidity 'zqv' and air pressure 'zp'
  !   (rdv = r_d / r_v
  !        = gas constant for dry air / gas constant for water vapour )
  !---------------------------------------------------------------------------
  !
  fq2pv  =  MAX( epsy , zqv ) * zp / (rdv + zqv*(c1-rdv))
  !
END FUNCTION fq2pv

!-------------------------------------------------------------------------------

ELEMENTAL INTEGER FUNCTION ileap  ( iyy )
  !--------------------------------------
  INTEGER                   , INTENT (IN)  ::  iyy    ! year [YYYY]
  !---------------------------------------------------------------------------
  ! detects leap year:  ileap(iyy) = 0:  no leap year,
  !                     ileap(iyy) = 1:  leap year
  !---------------------------------------------------------------------------
  !
  ileap =   IABS( MOD(iyy,  4) -  4) /  4   & ! every     4 years is a leapyear
          - IABS( MOD(iyy,100) -100) /100   & ! but every 100 ys. is no leapyear
          + IABS( MOD(iyy,400) -400) /400     ! but every 400 ys. is a leapyear
  !
END FUNCTION ileap

!-------------------------------------------------------------------------------

END MODULE src_obs_cdfin_util
