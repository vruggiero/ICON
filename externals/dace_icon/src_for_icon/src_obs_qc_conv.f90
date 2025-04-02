!+ Source module for COSMO quality control of conventional observations
!-------------------------------------------------------------------------------

MODULE src_obs_qc_conv

!-------------------------------------------------------------------------------
!
! Description:
!   This module 'src_obs_qc_conv' contains the routines which perform the
!   quality control of conventional observations that are also used by the
!   nudging scheme in COSMO. This includes:
!    - multi-level reports (radiosondes, profilers, aircraft)
!    - upper-air single-level reports (aircraft, satob)
!    - single-level surface reports (synop, ship, buoy)
!    - integrated water vapour (IWV) derived from ground-based GPS reports
!
!   The quality control includes a threshold quality control for each individual
!   observation (similar to a 'first guess check'), multi-level checks, and a
!   height and thickness check for (not too short) multi-level temperature data
!   designed to detect temperature biases over large vertical extents.
!   A spatial consistency check of integrated water vapour derived from
!   radiosonde humidity is not included here. This is only done inside the
!   nudging scheme of the COSMO model.
!
!   Note: This module belongs to a group of modules for the 'COSMO observation
!   ----  operators including quality control for conventional observations'.
!         It is shared between the COSMO model and the DACE program package !
!
! Method:
!   This module contains the following procedures:
!   - sing_quality_cntl     : quality control of individual single-level obs
!   - iwv_quality_cntl      : quality control of IWV 'obs'
!   - ps_quality_cntl       : quality control of (near-) surface pressure obs
!   - mult_obs_qc_fg        : quality control of multi-level obs
!   - mult_quality_cntl     : quality control of 1 level of a multi-level report
!   - mult_obs_qc_dz        : height / thickness QC for multi-level temperature
!
!   All routines have to be called by external driving routines, except for:
!   - mult_quality_cntl     : called by mult_obs_operator_qc1
!
!   Note: This module does not contain any communication between PE's.
!
!   Note on the interface of this module:
!   -------------------------------------
!   All the variables used from external data-modules are parameters,
!   Only the the subroutine (input) arguments, have to be set to their correct
!   values in the calling routines outside this module.
!
! Current Code Owner (for COSMO and for DACE):
!  DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! -------    ----       ----
! V4_28        2013/07/12 Christoph Schraff
!  Initial release, based on routines in module 'src_mult_local.f90' and
!                   'src_sing_local.f90' of COSMO V4_22.
!  Previous major milestones in the developement of 'src_mult_local.f90'
!  related to the procedures in the current module:
!   - 1998 : Initial release.
!   - 2000 : Introduction of hydrostatic height and thickness check for
!            multi-level temperature;  introduction of multi-level check.
!            Observed pressure tendency used for pressure QC thresholds.
!   - 2001 : Revised QC thresholds for aircrafts and PILOTs following ECMWF.
!            Rejection of aircraft temperature, if wind is rejected, and vice
!            versa.
!   - 2004 : New formulation of multi-level check (like in GME-OI).
!   - 2004-2006 : New option for reduced, stability-dependent QC thresholds
!                 for humdity. Revised basic humdity thresholds (new OI values).
!   - 2012 : Introduction of observation operator and QC of virtual temperature.
!  Changes with respect to V4_22:
!   - Split-up of old large routines into new routines with slim interfaces.
!   - Simulated observations instead of observation increments stored in SOR
!     array 'zsobdy', used for the feedback files.
!   - Minor bug fixes, e.g. QC threshold of passive humidity obs from upper-air
!     single-level reports, recording rejection of screen-level obs, removal of
!     QC (flags) for aircraft 'height', etc.
!   - New Option for QC of surface pressure obs against lateral boundary fields.
!   - Set GPS report to passive, if gross error occurs in QC check.
!   - Statement functions replaced by intrinsic functions.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Bug fix (array bound): if thickness is checked but not rejected then compute
!                         'zdzobs' only if (ilevb <= nlev).
!  Replaced ireals by wp (working precision) (OF)
! @VERSION@    @DATE@     Christoph Schraff + Andreas Rhodin
!  - Set subroutines used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - All 'iintegers' removed and replaced by standard integers.
!  - Use of 'kind_parameters' instead of 'data_parameters'.
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

USE data_obs_lib_cosmo , ONLY :   &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
!   c3600      ,& ! standard real constant 3600.0
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy          ! = 1.E-8_wp : commonly used very small value > 0

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

!   nupr          ! unit number of file for all the remaining information

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_record , ONLY :   &

! 1.1 ODR report header format
! ----------------------------

!   mxrhed       ,& ! header length of multi-level reports
!   nhilon       ,& ! longitude of observing station
!   nhjlat       ,& ! latitude  of observing station
!   nhalt        ,& ! station altitude [m]
!   nhtime       ,& ! time of observat. in forecast hours
!   nhzio        ,& ! latitude  of obs. station (or lowest datum) in g.pt units
!   nhzjo        ,& ! longitude of obs. station in grid pt. units
!   nhtvip       ,& ! observed multi-lev pressure interpol to lowest model level

!   mxrhdf       ,& ! header length of multi-level reports
    mxghdf       ,& ! header length of GPS reports
!   mxthdf       ,& ! header length of satellite retrieval reports
!   nhio         ,& ! (local) x-coord. of grid pt. assigned to obs
!   nhjo         ,& ! (local) y-coord. of grid pt. assigned to obs
!   nhitot       ,& ! global x-coord. of grid pt. assigned to obs
!   nhjtot       ,& ! global y-coord. of grid pt. assigned to obs
!   nhobtp       ,& ! observation type
!   nhcode       ,& ! code type
!   nhschr       ,& ! station characteristics (packed as in VOF)
    nhflag       ,& ! report flags (obs type, surf., alt., sta ID)
    nhpass       ,& ! flag for report being set to 'passive' (as in VOF)
!   nhqcfw       ,& ! status of QC and of writing to feedobs files (+ p-QC flag)
!   nhnlev       ,& ! number of obs. levels (for multi-level reports)
!   nhvqcf       ,& ! for satellite retrieval: threshold quality control flags
!   nhaexi       ,& ! for conventional report: flag for exist. of wind or temp.
!   nhuexi       ,& ! flag for existence of wind data      in multi-level report
!   nhtexi       ,& ! flag for existence of temperat. data in multi-level report
!   nhqexi       ,& ! flag for existence of humidity data  in multi-level report

! 1.3 ODR body format
! -------------------

!   maxrsl       ,& ! max. number of levels in multi-level ODR
    mxrbdy       ,& ! body length of multi-level reports
    nbtu         ,& ! u wind component [m/s]
    nbtv         ,& ! v wind component [m/s]
    nbtt         ,& ! temperature [K]
    nbtrh        ,& ! relative humidity [/]
    nbtp         ,& ! pressure [Pa]
    nbtz         ,& ! height [m]
!   nbtuer       ,& ! error of observed wind component
!   nbtter       ,& ! error of observed temperature
    nbtqer       ,& ! error of observed relative humidity
!   nbtzer       ,& ! error of observed height
    nbtlop       ,& ! LOG( pressure )
    mxrbdf       ,& ! body length of multi-level reports
    nbtflg       ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbtqcf       ,& ! threshold quality control flags      (see below: 'nb?qcf')
!   nbtlsg       ,& ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid       ,& ! level identity          (bit pattern, see below: 'nb?lid')
    mxsbdy       ,& ! body length of single-level reports
    nbsu         ,& ! u wind component [m/s]
    nbsv         ,& ! v wind component [m/s]
    nbst         ,& ! temperature [K]
    nbsrh        ,& ! relative humidity [/]
    nbsp         ,& ! pressure [Pa]
    nbsz         ,& ! height [m]
    nbsqer       ,& ! error of observed relative humidity
    nbscl        ,& ! low cloud cover
    mxsbdf       ,& ! body length of single-level reports
    nbserr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbsqcf       ,& ! threshold quality control flags (as in VOF)
    nbsflg       ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    mxgbdf       ,& ! body length of GPS reports
    nbgflg       ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbgqcf          ! threshold quality control flags      (see below: 'nb?qcf')

USE data_obs_record , ONLY :   &

! 1.4 Bit patterns for packed information in ODR body
! ---------------------------------------------------

    nvfubp       ,& ! bit pos. for main flag on wind                  nb?flg
    nvftbp       ,& ! bit pos. for main flag on temperature             "
    nvfqbp       ,& ! bit pos. for main flag on humidity                "
    nvfzbp       ,& ! bit pos. for main flag on pressure / geopot.      "
    nvfgbp       ,& ! bit pos. for main flag on IWV / ZPD               "
    nvfaoc       ,& ! no. of bits occ. by each main flag                "
    nvfbps       ,& ! bit pos. for main flags: 4: above 300 hPa         "
!   nvfboc       ,& ! no. of bits occ. for main flags                   "
    nvlidp       ,& ! level id. bit pattern                           nb?lid
    nvlido       ,& ! no. bits occ. by each indicator in level id.      "

!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru         ,& ! bit pos. for status/QC flags for horiz. wind nb?err/nb?qcf
    nvrt         ,& ! bit pos. for status/QC flags for temperature       "
    nvrq         ,& ! bit pos. for status/QC flags for humidity          "
    nvrz         ,& ! bit pos. for status/QC flags for pressure/height   "
    nvriwv       ,& ! bit pos. for status/QC flags for IWV               "
    nvrzpd       ,& ! bit pos. for status/QC flags for ZPD               "
    nvrzbc       ,& ! bit pos. for temporary flag: QC ag. LBC pressure   "

! 1.5 Further quantities related to ODR
! -------------------------------------

    ilstidp         ! character length used for printing the station ID
!   ystid           ! obs. station identity to be printed

! 4. Masking constants
! --------------------

!   nibits          ! masking constants

! end of data_obs_record

!-------------------------------------------------------------------------------

USE data_obs_qc_limits , ONLY :   &

! 2. Table for pressure dependent quality control thresholds
! ----------------------------------------------------------

    nqclev       ,& ! number of levels in the correlation scale tables
!   tabqcp       ,& ! levels of the error / threshold tables
    tabqclp      ,& ! LOG( tabqcp(15) )
    qcvsond      ,& ! (root of) radiosonde (TEMP, PILOT) wind error variance
    qczsond      ,& ! (root of) radiosonde height error variance
    qctsond      ,& ! (root of) radiosonde temperature error variance
    qcvairp      ,& ! (root of) AIREP wind error variance
    qctairp      ,& ! (root of) AIREP temperature error variance
    qctf         ,& ! temporal factor for thresholds of upper-air obs.
    qctfsu       ,& ! temporal factor for thresholds of surface-level obs.
    qcvairf      ,& ! additional factor used for airep wind thresholds
    qctairf      ,& ! additional factor used for airep temperature thresholds
    qcvpilf      ,& ! additional factor used for pilot wind thresholds
    qcvdrbf      ,& ! additional factor used for dribu wind thresholds
    qczdrbf      ,& ! additional factor used for dribu pressure thresholds
    qcvairl      ,& ! background wind speed rejection limit for zero AIREP winds
    qczcorl      ,& ! radiosonde height error correlation matrix
    qcvfz        ,& ! factor to qczsond for height (and thickness) threshold
!   qcftsiv      ,& ! reduction factor to the threshold used for IWV-SSC
    qcqinvf      ,& ! humidity f.g. error enhancement per [K] of obs'd inversion
    qcqlapf      ,& ! max. humidity f.g. error enhance. if zero lapse r. or inv.
    qcqlapl      ,& ! lapse rate limit to enhance humidity f.g. error
!   qcqladz      ,& ! maximum stability factor enhancement
    qcqlapi      ,& ! scaling factor of positive lapse rate for inversion term
    qcqfge1      ,& ! first guess humidity error at latitude > qcqlatl
    qcqfge2      ,& ! first guess humidity error at latitude < qcqlatl
    qcqlatl      ,& ! latitude limit (for first guess humidity error)
    qcqfgef      ,& ! factors for first guess humidity errors
    qcftend      ,& ! fraction of observed pressure tendency added to threshold
    qcfctnd      ,& ! reduction to namelist threshold, if tendency is observed
    qcfgood      ,& ! reduction of surface pressure threshold for 'good' obs.
    qcfpbad      ,& ! increase of pressure threshold for 'probably bad' obs.
!   qcfspat      ,& ! threshold reduction for the spatial consistency check
!   qcfbias      ,& ! fraction of bias added to spatial consistency threshold
!   dtchkps      ,& ! spatial check: radius of infl. for linear temporal weights
    nqcalv       ,& ! number of analysis layers for multi-level check plus 1
    tabanp       ,& ! pressure limits of analysis layers for multi-level check
    nanstd          ! number of analysis layers within the standard layers

! end of data_obs_qc_limits

!-------------------------------------------------------------------------------

USE mo_fdbk_tables, ONLY :   &

    OT_AIREP     ,& ! observation type for aircraft obs
    OT_DRIBU     ,& ! observation type for DRIBU drifting buoy obs
    OT_TEMP      ,& ! observation type for TEMP radiosonde
    OT_PILOT     ,& ! observation type for PILOT (incl. wind profiler, RASS,...)
    OT_SATOB     ,& ! observation type for AMV / SATOB winds
    OT_SATEM     ,& ! observation type for satellite retrievals
    FL_NO_OBS       ! no observations in report

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
          ! public subroutines
PUBLIC :: sing_quality_cntl,  ps_quality_cntl,  mult_obs_qc_fg,  mult_obs_qc_dz

! for nudging in COSMO model code only:
PUBLIC :: iwv_quality_cntl

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------------
!+ Module procedure for quality control of individual observations
!-------------------------------------------------------------------------------

SUBROUTINE sing_quality_cntl ( zobbdy, vip, ps, kb, kobtyp                     &
                             , tabdif, qcvf, qcc, ndqc, lwrqc                  &
                             , mzobbd, nlqc, zyqc )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine performs the threshold quality control for horizontal wind,
!   temperature, and (relative) humidity observations from one surface-level
!   or upper-air single-level conventional report.
!   Rejection flags are updated, including a flag to indicate that an
!   upper-air report is below model surface (pressure).
!   Documentation of rejection of observations is prepared.
!
! Method:
!   Checks of 'observation minus interpolated model value' against threshold.
!   The thresholds may depend on the pressure level of the observation point
!   (for upper-air wind and temperature) and on the difference between the
!   observation time and the current model time.
!
! Written by        :  Christoph Schraff, DWD  (original version: 13.08.98)
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
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    kb         ,& ! switch for type of obs report, and index:
                  ! - kb > 0 : index of adjacent model level below upper-air obs
                  ! - kb = 0 : upper-air obs is above top main model level
                  ! - kb =-1 : surface-level obs being QC'ed
    kobtyp     ,& ! observation type
    ndqc          ! dimension of 'zyqc'

  REAL (KIND=wp) , INTENT (IN)         ::       &
    vip  (4)   ,& ! model value inter-/extrapolated to the obs. level (pressure)
                  !   (1, 2: (u,v)-wind comp.  ;  3: temperature ;
                  !      4 : relative humidity )
    ps         ,& ! model surface pressure (at obs location)
    tabdif     ,& ! time distance between the obs. and the model state
    qcc  (4)   ,& ! vertically constant part of QC thresholds
    qcvf (4)      ! multiplication factor to the vertically varying part of the
                  !   QC thresholds (as defined in 'data_obs_qc_limits')

  LOGICAL        , INTENT (IN)          ::      &
    lwrqc         ! data with obs time close to current time, which are rejected
                  ! by the threshold quality control, are printed to a separate
                  ! file at the current timestep

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy (mxsbdy)    ! 1 level of multi-level obs body (format: see 'osgbdy')

  INTEGER        , INTENT (INOUT)      ::       &
    mzobbd (mxsbdf)    ! 1 level of multi-level obs body (format: see 'mosgbd')

  INTEGER        , INTENT (OUT)        ::       &
    nlqc               ! number of QC messages

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    zyqc   (ndqc,11)   ! information for QC control messages

! Local parameters:
! ----------------
  REAL (KIND=wp) , PARAMETER  :: &
    c100r = 0.01_wp    ! 1/100

! Local variables (scalars or automatic arrays):
! ---------------------------------------------

! REAL (KIND=wp) ::  &
!   qcc (4) = c0, 500._wp, c0, c0   ! constant part of QC thresholds

  INTEGER        ::  &
    ilva   , ilvb    ,& ! indices of correlation scale table levels that are
                        ! adjacent to the specified model level
    mflgqc (4)          ! threshold quality control bit flag

  REAL (KIND=wp) ::  &
    qctu             ,& ! QC threshold for wind vector increment
    qctt             ,& ! QC threshold for temperature increment
    qctrh            ,& ! QC threshold for relative humidity increment
    uerror           ,& ! mean wind (vector) obs error at current pressure level
    terror           ,& ! mean temperature obs error at current pressure level
    qerror           ,& ! mean rel. humidity obs error at current pressure level
    zqernor          ,& ! humidity obs error / tot error normalised by f.g error
    zlopob           ,& ! log( pressure ) at obs. point
    zf               ,& ! factor for vert. interpol. from q.c. thresholds table
    vecdif           ,& ! length of the obs. increment wind vector
    zobinc           ,& ! scalar observation increment
    zyqcioff            ! offset for indicator value for type of QC message
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine sing_quality_cntl
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Determination of thresholds
!-------------------------------------------------------------------------------

  IF (kb > 0) THEN

    ! thresholds for upper-air single-level obs
    ! -----------------------------------------

    zlopob  = LOG( zobbdy(nbsp) )

    !   vertically dependent part of thresholds
    IF ((zlopob >= tabqclp(1)) .OR. (zlopob <= tabqclp(nqclev))) THEN
      ilva = 1
      IF (zlopob <= tabqclp(nqclev)) ilva = nqclev
      ilvb = ilva
      zf   = c1
    ELSE
      ilva = 2
      DO WHILE (zlopob <= tabqclp(ilva))
        ilva = ilva + 1
      ENDDO
      ilvb = ilva - 1
      zf   = (tabqclp(ilvb) - zlopob) / (tabqclp(ilvb) - tabqclp(ilva))
    ENDIF
    uerror    = qcvairf * ((c1-zf) *qcvairp(ilvb) + zf *qcvairp(ilva))
    terror    = qctairf * ((c1-zf) *qctairp(ilvb) + zf *qctairp(ilva))
    qerror    = 0.7_wp

    !   humidity thresholds
    IF (zobbdy(nbsqer) > rmdich) THEN
      zqernor = ABS( zobbdy(nbsqer) ) / qcqfge1
      zqernor = SQRT( MIN( zqernor *zqernor + c1 , 4._wp ) )
      qerror  = MIN( zqernor * qcqfgef(3) * qcqfge1 , 0.7_wp )
    ENDIF

    !   total thresholds
    qctu   = (c1 + qctf(1) *tabdif) * (qcc(1) + qcvf(1) *uerror)
    qctt   = (c1 + qctf(3) *tabdif) * (qcc(3) + qcvf(3) *terror)
    qctrh  = (c1 + qctf(4) *tabdif) * (qcc(4) + qcvf(4) *qerror)
    IF (kobtyp == OT_DRIBU)  qctu  =  qctu * qcvdrbf

  ELSEIF (kb <= -1) THEN

    ! thresholds for surface-level obs
    ! --------------------------------

    qctu   =  (c1 + qctfsu(1) *tabdif) * qcc(1)
    qctt   =  (c1 + qctfsu(3) *tabdif) * qcc(3)
    qctrh  =  (c1 + qctfsu(4) *tabdif) * qcc(4)

  ENDIF

  !   initialize 'mflgqc', 'zqqc', 'nlqc'
  mflgqc (:) =  0
  zyqc (:,:) = c0
  nlqc       =  0
  zyqcioff   = c0
  IF (kb <= -1)  zyqcioff = 4._wp

!-------------------------------------------------------------------------------
!  Section 2: Check WIND vector obs increment against threshold
!-------------------------------------------------------------------------------

  IF ((zobbdy(nbsu ) > rmdich) .AND. (kb /= 0)) THEN
    vecdif = SQRT( (zobbdy(nbsu) - vip(1))**2 + (zobbdy(nbsv) - vip(2))**2 )
    IF (     (vecdif > qctu)                                                   &
    !   reject zero aircraft / satob wind if model wind speed exceeds 'qcvairl'
        .OR. (      ((kobtyp == OT_AIREP) .OR. (kobtyp == OT_SATOB))           &
              .AND. (vecdif > qcvairl)                                         &
              .AND. (MAX( ABS( zobbdy(nbsu) )                                  &
                        , ABS( zobbdy(nbsv) ) ) < epsy)))   mflgqc (1) = 1
    !   reject aircraft wind if threshold QC rejects temperature
    IF ((kobtyp == OT_AIREP) .AND. (zobbdy(nbst) > rmdich)) THEN
      IF (ABS( zobbdy(nbst) - vip(3) ) > qctt)    mflgqc (1) = 1
    ENDIF

    IF (mflgqc(1) == 1) THEN
      IF ((lwrqc) .OR. (      (      BTEST( mzobbd(nbserr),nvru ))             &
                        .AND. (.NOT. BTEST( mzobbd(nbsqcf),nvru )))) THEN
        nlqc = nlqc + 1
        zyqc (nlqc, 1) = zyqcioff + 1.0_wp
        zyqc (nlqc, 2) = zobbdy(nbsp) * c100r
        zyqc (nlqc, 5) = qctu
        zyqc (nlqc, 6) = zobbdy(nbsu)
        zyqc (nlqc, 7) = vip (1)
        zyqc (nlqc, 8) = zobbdy(nbsv)
        zyqc (nlqc, 9) = vip (2)
        zyqc (nlqc,10) = vecdif
      ENDIF
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!  Section 3: Check TEMPERATURE obs increment against threshold
!-------------------------------------------------------------------------------

  IF ((zobbdy(nbst ) > rmdich) .AND. (kb /= 0)) THEN
    zobinc  =  ABS( zobbdy(nbst) - vip(3) )
    !   (reject aircraft temperature if aircraft wind is rejected)
    IF ((zobinc > qctt) .OR. ((kobtyp == OT_AIREP) .AND. (mflgqc(1) == 1)))    &
      mflgqc (2) = 1

!   IF (mflgqc(2) == 1) THEN
    IF (zobinc > qctt) THEN
      IF ((lwrqc) .OR. (      (      BTEST( mzobbd(nbserr),nvrt ))             &
                        .AND. (.NOT. BTEST( mzobbd(nbsqcf),nvrt )))) THEN
        nlqc = nlqc + 1
        zyqc (nlqc, 1) = zyqcioff + 3.0_wp
        zyqc (nlqc, 2) = zobbdy(nbsp) * c100r
        zyqc (nlqc, 5) = qctt
        zyqc (nlqc, 6) = zobbdy(nbst)
        zyqc (nlqc, 7) = vip(3)
        IF (NINT( zyqcioff ) == 4)  zyqc (nlqc, 8) = zobbdy(nbst) - vip(3)
      ENDIF
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!  Section 4: Check RELATIVE HUMIDITY obs increment against threshold
!-------------------------------------------------------------------------------

  IF ((zobbdy(nbsrh) > rmdich) .AND. (kb /= 0)) THEN
    zobinc  =  ABS( zobbdy(nbsrh) - vip(4) )
    !   (reject humidity if temperature is rejected)
    IF ((zobinc > qctrh) .OR. (mflgqc(2) == 1))  mflgqc (3) = 1

    IF (mflgqc(3) == 1) THEN
!   IF (zobinc > qctrh) THEN
      IF ((lwrqc) .OR. (      (      BTEST( mzobbd(nbserr),nvrq ))             &
                        .AND. (.NOT. BTEST( mzobbd(nbsqcf),nvrq )))) THEN
        nlqc = nlqc + 1
        zyqc (nlqc, 1) = zyqcioff + 4.0_wp
        zyqc (nlqc, 2) = zobbdy(nbsp) * c100r
        zyqc (nlqc, 5) = qctrh
        zyqc (nlqc, 6) = zobbdy(nbsrh)
        zyqc (nlqc, 7) = vip(4)
        !   message for aircraft humidity also if only T-diff exceeds threshold
        IF (zobinc <= qctrh) THEN
          zyqc (nlqc, 8) = zobbdy(nbst)
          zyqc (nlqc, 9) = vip(3)
          IF (NINT( zyqcioff ) == 4)  zyqc (nlqc,10) = zobbdy(nbsrh) - vip(4)
          zyqc (nlqc, 1) = - zyqc(nlqc,1)
        ENDIF
      ENDIF
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!  Section 5: Set flag for obs pressure level being out of (model) range
!-------------------------------------------------------------------------------

  !   (also need to set flag if (kb == 0), as it can become > 0 within dtqc !)

  IF ((kb == 0) .OR. ((kb > 0) .AND. (zobbdy(nbsp) > ps)))  mflgqc (4) = 1

!-------------------------------------------------------------------------------
!  Section 6: Set threshold quality control flags in ODR
!-------------------------------------------------------------------------------

  CALL MVBITS ( mflgqc(1), 0, 1, mzobbd(nbsqcf), nvru )
  CALL MVBITS ( mflgqc(2), 0, 1, mzobbd(nbsqcf), nvrt )
  CALL MVBITS ( mflgqc(3), 0, 1, mzobbd(nbsqcf), nvrq )
  IF (kb >= 0)                                                                 &
    CALL MVBITS ( mflgqc(4), 0, 1, mzobbd(nbsqcf), nvrz )

!-------------------------------------------------------------------------------
! End of module procedure sing_quality_cntl
!-------------------------------------------------------------------------------

END SUBROUTINE sing_quality_cntl



!-------------------------------------------------------------------------------
!+ Module procedure for quality control of individual IWV observations
!-------------------------------------------------------------------------------

SUBROUTINE iwv_quality_cntl ( ziqobs, ziqmod, ziqsat, qcsiq, qcciq, qctfq      &
                            , tabdif, lwrqc , mzobbd, mzobhd , zyqc )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine performs the threshold quality control for an integrated water
!   vapour (IWV) value derived from a zenith total path delay (ZPD) observation
!   from a ground-based GPS report.
!   Rejection flags are updated.
!   Documentation of rejection of observations is prepared.
!
! Method:
!   Checks of 'observation minus interpolated model value' against threshold.
!   QC threshold, default: 0.15* IWV-saturated + 1mm.
!   ------------  The threshold is a function of the IWV of the saturated
!                 model temperature profile plus a constant value. In addition,
!                 it depends on the distance between obs and current model time.
!   The obs is also rejected, if observed IWV is < 2 kg/m2, i.e. if it is below
!   a lower limit. In this case, the gross error pre-processing flagging is
!   applied.
!
! Written by        :  Christoph Schraff, DWD  (original version: 04.04.13)
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
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  REAL (KIND=wp) , INTENT (IN)         ::       &
    ziqobs      ,& ! observed IWV (from GPS height, poss. ice-to-water adjusted)
    ziqmod      ,& ! model    IWV (from GPS antenna height)
    ziqsat      ,& ! IWV from saturated model temperature profile
    qcsiq       ,& ! IWV QC threshold, as a fraction of 'ziqsat'
    qcciq       ,& ! constant part of QC threshold for IWV
    qctfq       ,& ! temporal factor for humidity thresholds
    tabdif         ! time distance between the obs. and the model state

  LOGICAL        , INTENT (IN)          ::      &
    lwrqc          ! data with obs time close to current time which are rejected
                   ! by the threshold quality control, are printed to a separate
                   ! file at the current timestep

  INTEGER        , INTENT (INOUT)      ::       &
    mzobbd (mxgbdf) ,& ! ground-based GPS obs body   (format: see 'mogpbd')
    mzobhd (mxghdf)    ! ground-based GPS obs header (format: see 'mogphd')

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    zyqc   (11)     ! information for QC control messages

! Local parameters:
! ----------------
  REAL (KIND=wp) , PARAMETER ::  &
    thriq  =  2.0_wp   ! minimum observed IWV value for assimilation

! Local variables (scalars or automatic arrays):
! ---------------------------------------------

  INTEGER        ::  &
    mflgqc              ! threshold quality control bit flag

  REAL (KIND=wp) ::  &
    qctiq               ! time-dependent quality control threshold for IWV
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine iwv_quality_cntl
!-------------------------------------------------------------------------------

  ! determine the quality control threshold
  ! ---------------------------------------
  qctiq   = (c1 + qctfq *tabdif) * (qcsiq *ziqsat + qcciq)

  ! determine QC status (0 - 3 : good / probably good / probably bad / bad)
  ! -----------------------------------------------------------------------
  mflgqc = NINT( c05 + SIGN(c05, ABS( ziqmod -ziqobs) - qctiq) )

  ! information for control output
  ! ------------------------------
  zyqc (:) = c0
  IF ((mflgqc == 1) .OR. (ziqobs <= thriq)) THEN
    IF ((lwrqc) .OR. (      (      BTEST( mzobbd(nbgerr),nvriwv ))             &
                      .AND. (.NOT. BTEST( mzobbd(nbgqcf),nvriwv )))) THEN
      zyqc (1) = 16.0_wp
!     zyqc (2) = zobbdy(nbsp) * 0.01_wp
      zyqc (5) = qctiq
      zyqc (6) = ziqobs
      zyqc (7) = ziqmod
    ENDIF
  ENDIF

  ! set threshold quality control flag in the ODR
  ! ---------------------------------------------
  CALL MVBITS ( mflgqc, 0, 1, mzobbd(nbgqcf), nvriwv )
  CALL MVBITS ( mflgqc, 0, 1, mzobbd(nbgqcf), nvrzpd )

  ! set gross error flag in the ODR !
  ! and set no-obs report flag      !
  ! -------------------------------
  IF (ziqobs <= thriq) THEN
    mzobbd (nbgerr) = IBCLR ( mzobbd(nbgerr), nvriwv )
    mzobbd (nbgerr) = IBCLR ( mzobbd(nbgerr), nvrzpd )
    mzobbd (nbgflg) = IBSET ( mzobbd(nbgflg), nvfgbp +nvfbps(3) )
    mzobhd (nhpass) = 2
    mzobhd (nhflag) = IBSET ( mzobhd(nhflag), FL_NO_OBS )
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure iwv_quality_cntl
!-------------------------------------------------------------------------------

END SUBROUTINE iwv_quality_cntl



!-------------------------------------------------------------------------------
!+ Module procedure for quality control of individual pressure observations
!-------------------------------------------------------------------------------

SUBROUTINE ps_quality_cntl ( np, zobps, zpsvim, lsing, kobtyp, tabdif          &
                           , qccps, fqclbcp, zdpsdt, lwrqc, mxbdf              &
                           , mzobbd , mflgqc, zyqc )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine performs the threshold quality control for (near-) surface
!   pressure observation data from one surface-level or one multi-level
!   (radiosonde) report.
!   Rejection flags are updated.
!   Documentation of rejection of observations is prepared.
!
! Method:
!   Checks of 'observation minus interpolated model value' against threshold.
!   The threshold may depend on observed 3-hourly surface pressure tendency
!   on the difference between the observation time and the current model time.
!
! Written by        :  Christoph Schraff, DWD  (original version: 13.08.98)
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
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    np          ,& ! number of pressure model fields
                   !   (e.g. 2: original field, boundary field)
    kobtyp      ,& ! observation type
    mxbdf          ! dimension of 'mzobbd'

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobps  (np) ,& ! (near-) surface pressure observation value
    zpsvim (np) ,& ! model pressure interpolated to obs level ('ilvvip(1)')
    zdpsdt      ,& ! 3-hrly surface pressure tendency
    tabdif      ,& ! time distance between the obs. and the model state
    qccps       ,& ! basic QC threshold for surface pressure
    fqclbcp        ! factor to threshold used for QC checks against lateral
                   !   boundary fields (LBC)

  LOGICAL        , INTENT (IN)          ::      &
    lsing       ,& ! surface pressure obs is from single-level surface report
    lwrqc          ! data with obs time close to current time which are rejected
                   ! by the threshold quality control, are printed to a separate
                   ! file at the current timestep

  INTEGER        , INTENT (INOUT)      ::       &
    mzobbd (mxbdf) ! 1 level of multi-level obs body (format: see 'mosgbd')

  INTEGER        , INTENT (OUT)        ::       &
    mflgqc         ! threshold quality control bit flag
                   !   (0 - 3 : good / probably good / probably bad / bad)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    zyqc   (11)    ! information for QC control messages

! Local parameters:
! ----------------
  REAL (KIND=wp) , PARAMETER  :: &
    c100r = 0.01_wp    ! 1/100

! Local variables (scalars or automatic arrays):
! ---------------------------------------------

! REAL (KIND=wp) ::  &
!   qcc (4) = c0, 500._wp, c0, c0   ! constant part of QC thresholds

  LOGICAL        ::  &
    lbc         ,& ! quality control also against lateral boundary fields
    lout           ! write QC information to 'zyqc' for control output

  INTEGER        ::  &
    mflg   (np)    ! threshold quality control bit flag

  REAL (KIND=wp) ::  &
    qctps  (np)    ! QC thresholds for (surface) pressure
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine ps_quality_cntl
!-------------------------------------------------------------------------------

  lbc = (np == 2)
  zyqc (:) = c0

  ! determine the quality control threshold
  ! ---------------------------------------
  IF ((zdpsdt > rmdich) .AND. (qccps > epsy)) THEN
               qctps (1) =  qcfctnd *qccps
    IF (lbc)   qctps (2) =  fqclbcp *qctps(1)  +  qcftend *ABS( zdpsdt )
               qctps (1) =           qctps(1)  +  qcftend *ABS( zdpsdt )
  ELSE
               qctps (1) =  (c1 + qctf  (2) *tabdif) * ABS( qccps )
    IF (lsing) qctps (1) =  (c1 + qctfsu(2) *tabdif) * ABS( qccps )
    IF (kobtyp == OT_DRIBU)  qctps (1) =  qctps(1) * qczdrbf
    IF (lbc)   qctps (2) =  fqclbcp *qctps(1)
  ENDIF

  ! determine QC status (0 - 3 : good / probably good / probably bad / bad)
  ! -----------------------------------------------------------------------
  mflg (1) = 0
  IF (ABS(zobps(1) -zpsvim(1)) > qctps(1)*qcfgood)  mflg (1) = 1
  IF (ABS(zobps(1) -zpsvim(1)) > qctps(1)        )  mflg (1) = 2
  IF (ABS(zobps(1) -zpsvim(1)) > qctps(1)*qcfpbad)  mflg (1) = 3
  mflgqc = mflg(1)
  IF (lbc) THEN
    mflg (2) = 0
    IF (ABS(zobps(2) -zpsvim(2)) > qctps(2)*qcfgood)  mflg (2) = 1
    IF (ABS(zobps(2) -zpsvim(2)) > qctps(2)        )  mflg (2) = 2
    IF (ABS(zobps(2) -zpsvim(2)) > qctps(2)*qcfpbad)  mflg (2) = 3
    mflgqc = MAX( mflg(1) , mflg(2) )
  ENDIF
           mflg (1) = mflg(1) / 2
  IF (lbc) mflg (2) = mflg(2) / 2

  ! information for control output
  ! ------------------------------
  IF (mflgqc >= 2) THEN
    IF (lsing) THEN
      lout =       (.NOT. BTEST( mzobbd(nbsflg), nvfzbp+nvfbps(4) ))           &
             .AND. ((lwrqc) .OR. (      (      BTEST( mzobbd(nbserr),nvrz  ))  &
                                  .AND. (.NOT. BTEST( mzobbd(nbsqcf),nvrz  ))  &
                                  .AND. (.NOT. BTEST( mzobbd(nbsqcf),nvrzbc))))
    ELSE
      lout = lwrqc
    ENDIF
    IF (lout) THEN
      IF (mflg(1) >= 1) THEN
        zyqc (1) = 6._wp
        zyqc (2) = zobps (1)  * c100r
        zyqc (5) = qctps (1)  * c100r
        zyqc (6) = zobps (1)  * c100r
        zyqc (7) = zpsvim(1)  * c100r
      ELSEIF (lbc) THEN
        zyqc (1) = 18._wp
        zyqc (2) = zobps (2)  * c100r
        zyqc (5) = qctps (2)  * c100r
        zyqc (6) = zobps (2)  * c100r
        zyqc (7) = zpsvim(2)  * c100r
      ENDIF
      IF (.NOT. lsing)  zyqc (1) = 2._wp
    ENDIF
  ENDIF

  ! set threshold quality control flag in the ODR
  ! ---------------------------------------------
  IF (lsing) THEN
             CALL MVBITS ( mflg(1), 0, 1, mzobbd(nbsqcf), nvrz   )
    IF (lbc) CALL MVBITS ( mflg(2), 0, 1, mzobbd(nbsqcf), nvrzbc )
  ELSE
             CALL MVBITS ( mflg(1), 0, 1, mzobbd(nbtqcf), nvrz   )
    IF (lbc) CALL MVBITS ( mflg(2), 0, 1, mzobbd(nbtqcf), nvrzbc )
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure ps_quality_cntl
!-------------------------------------------------------------------------------

END SUBROUTINE ps_quality_cntl



!-------------------------------------------------------------------------------
!+ Module procedure for (threshold or first guess) quality control (QC)
!-------------------------------------------------------------------------------

SUBROUTINE mult_obs_qc_fg ( nlev, zobbdy, zt_o, kbotlev, ktoplev, kobtyp       &
                          , tabdif, mexi, oblat, qcvf, r_d, r_g                &
                          , lqc, lwrqc, lprqcm, nupr, ystid                    &
                          , mzobbd, vimtoob , nlqc, zyqc , qcc_in )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine performs the (threshold or first guess) quality control (QC)
!   of individual observations from multi-level reports, plus the so-called
!   multi-level checks.
!
! Method:
!   - quality control (QC) of individual observations :
!                             - by call of routine 'mult_quality_cntl'
!                             - stability-dependent enhancement factor to
!                               thresholds for humidity previously computed
!                               in Section 1.3
!   - multi-level QC checks : see (comments in) Section 2 (for active obs only!)
!
!
! Written by        :  Christoph Schraff, DWD  (original version: 18.09.97)
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
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    nlev       ,& ! number of vertical levels in report
    kbotlev    ,& ! number of obs levels below lowest    model level
    ktoplev    ,& ! number of obs levels above uppermost model level
    kobtyp     ,& ! observation type
    mexi   (3) ,& ! =1: active ; =0: no ; =-1: passive obs exist
                  !  (index 1: wind ; 2: temperature ; 3: humidity)
    nupr          ! file unit number for control output

  LOGICAL        , INTENT (IN)         ::       &
    lqc        ,& ! threshold quality control to be applied
    lwrqc      ,& ! data with obs time close to current time, which are rejected
                  ! by the threshold quality control, are printed to a separate
                  ! file at the current timestep
    lprqcm (2)    ! printout 1 for control at present station / node / time

! CHARACTER (LEN=ilstidp) , INTENT (IN)         ::       &
  CHARACTER (LEN=* )  , INTENT (IN)         ::       &
    ystid         ! obs. station identity to be printed

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy (nlev,mxrbdy) ,& ! multi-level obs. body (format: see 'omlbdy')
    zt_o   (nlev)        ,& ! observation-derived (dry bulb) temperature
    tabdif               ,& ! time distance between the obs. and the model state
    oblat                ,& ! latitude of observation
    qcvf   (4)           ,& ! multiplication factor to the vertically varying
                            !   part of QC thresholds (see 'data_obs_qc_limits')
    r_d                  ,& ! gas constant for dry air       (287.05)
    r_g                     ! acceleration due to gravity    (  9.80665)

  INTEGER        , INTENT (INOUT)      ::       &
    mzobbd (nlev,mxrbdf)    ! multi-level obs. body (format: see 'momlbd')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    vimtoob (nlev,5)        ! model values interpolated to obs. levels
                            !   (intent 'out' so as to set EXTRApolated values
                            !    (above + below profile) to missing after QC)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    zyqc    (nlev*4,11)     ! information for QC control messages

  INTEGER        , INTENT (OUT)        ::       &
    nlqc                    ! number of QC control messages

  REAL (KIND=wp) , INTENT (IN)    , OPTIONAL   ::       &
    qcc_in  (4)         ! vertically constant part of QC thresholds

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    c100r   =     0.01_wp     !

! Local scalars:
! -------------

  REAL (KIND=wp) ::  &
    qcc (4) = (/c0, 500._wp, c0, c0/)  ! constant part of QC thresholds
                                           !   (v, ps [Pa], T, RH)

  INTEGER        ::  &
    ivrs             ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of observed quantity: 1=(u,v); 3=T; 4=RH
    ilev             ,& ! vertical loop index over observation levels
    klev             ,& ! loop index over obs. levels below lowest model level
    ilvsfc , ilvpqm  ,& ! obs. level index of surface / humidity use limit
    ilvq             ,& ! loop index over standard pressure levels
    ilayan           ,& ! loop index over analysis layers (multi-level check)
    ilvb             ,& ! index of nearest level further below
    ilayb  , ilayt   ,& ! indices of base / top analysis layer to be rejected
    nrejstd          ,& ! number of flagged analysis layer within standard layer
    nflagw           ,& ! flag word
    nbtx             ,& ! ODR index for present variable
    nvrx             ,& ! bit pos. for status 'active' for variable
    nphpa            ,& ! pressure [hPa] of current obs. level
    ilvmlcb, ilvmlct ,& ! extent of profile rejected by multi-level check
    mandchk, mandrej ,& ! number of all / rejected mandatory level in report
    nact   , nrej       ! number of all / rejected level in (piece of) report

  REAL (KIND=wp) ::  &
    zf               ,& ! weight factor for vertical interpol. to the obs. point
    zlapse           ,& ! (approx) lapse rate
    zdtdt            ,& ! difference betw. T(ilev+1) and T extrapolated from
                        !   T(ilev) using lapse rate limit 'qcqlapl'(-0.0065K/m)
                        !   this depends on vertical distance over which the
                        !   lapse rate is computed.
    zlapdz           ,& ! stability factor enhancement depending on 'zdtdt'
    zlapsq           ,& ! stability factor: = 0 if less stable than 'qcqlapl'
                        !                   = 1 if inversion ; linear in between
    zlapiq           ,& ! scaled positive lapse rate at inversions
    zinvrs           ,& ! inversion in [K] between successive levels with q-obs
    zfqinv              ! factor by which humidity first guess error is enhanced

  LOGICAL        ::  &
    lblackls         ,& ! observation is blacklisted
    lprmlc           ,& ! printout for control of multi-level check
    lmandat             ! level is mandatory

! Local (automatic) arrays:
! -------------------------

  INTEGER        ::  &
    kqcflag (nlev,5) ,& ! quality control flag
    klevq   (nlev)   ,& ! level ind with humidity obs to be quality controlled
    kanaflg (nqcalv) ,& ! analysis layer which an observation (level) belongs to
    ilvmand (nqclev)    ! level indices of mandatory levels

  REAL (KIND=wp) ::  &
    fqerinv (nlev)      ! factor by which humidity first guess error is enhanced
                        ! in the presence of an observed inversion /stable layer

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine mult_obs_qc_fg
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Quality control: Local check for individual observations
!             --------------------------------------------------------
!             If (ABS(obs. incr.) > (threshold)) at an obs. level,
!             then set obs. to 'erroneous'
!             (in the nudging only for the subsequent 'dtqc' seconds)
!-------------------------------------------------------------------------------

! IF ((lqc) .AND. (kobtyp /= OT_SATEM)) THEN

    IF (PRESENT( qcc_in ))  qcc = qcc_in

    zyqc (:,:) = c0
    nlqc       =  0
    lprmlc     =  (lprqcm(1)) .AND. (lprqcm(2))

!-------------------------------------------------------------------------------
!  Section 1.1: Set quality control flag 'kqcflag' to be used in the
!               multi-level check (and in mult_quality_cntl):
!               - to missing (-1), if current obs is blacklisted, but active obs
!                                  exist in profile
!               - to missing (-1), if humidity above 300 hPa level
!               - to bad     ( 3), if obs error isn't defined for another reason
!               (for height obs, special rules apply, see below)
!-------------------------------------------------------------------------------

    fqerinv   (:) = c1
!   kqcflag (:,:) = 0

    DO ilev = 1 , nlev
      kqcflag (ilev,1) = -1
      kqcflag (ilev,2) = -1
      kqcflag (ilev,3) = -1
      kqcflag (ilev,4) = -1
      kqcflag (ilev,5) =  1
    ENDDO

    DO ilev = 1 , nlev
      nflagw = mzobbd(ilev,nbtflg)
      IF (zobbdy(ilev,nbtu) > rmdich) THEN
                                                     kqcflag (ilev,1) =  0
        IF (.NOT. BTEST( mzobbd(ilev,nbterr),nvru )) kqcflag (ilev,1) =  3
        lblackls  =  BTEST( nflagw, nvfubp+nvfbps(2) )
        IF ((lblackls) .AND. (mexi(1) ==  1))        kqcflag (ilev,1) = -1
        IF ((lblackls) .AND. (mexi(1) == -1))        kqcflag (ilev,1) =  0
      ENDIF
      IF (zobbdy(ilev,nbtt) > rmdich) THEN
                                                     kqcflag (ilev,3) =  0
        IF (.NOT. BTEST( mzobbd(ilev,nbterr),nvrt )) kqcflag (ilev,3) =  3
        lblackls  =  BTEST( nflagw, nvftbp+nvfbps(2) )
        IF ((lblackls) .AND. (mexi(2) ==  1))        kqcflag (ilev,3) = -1
        IF ((lblackls) .AND. (mexi(2) == -1))        kqcflag (ilev,3) =  0
      ENDIF
      IF (zobbdy(ilev,nbtrh) > rmdich) THEN
                                                     kqcflag (ilev,4) =  0
        IF (.NOT. BTEST( mzobbd(ilev,nbterr),nvrq )) kqcflag (ilev,4) =  3
        lblackls  =  BTEST( nflagw, nvfqbp+nvfbps(2) )
        IF ((lblackls) .AND. (mexi(3) ==  1))        kqcflag (ilev,4) = -1
        IF ((lblackls) .AND. (mexi(3) == -1))        kqcflag (ilev,4) =  0
        IF (      BTEST( nflagw,nvfqbp+nvfbps(4) ))  kqcflag (ilev,4) = -1
      ENDIF
      !   height obs
      IF (zobbdy(ilev,nbtz) > rmdich) THEN
                                                     kqcflag (ilev,2) =  0
        !   flag obs for multi-level check, if passive (dataset, gross,...error)
        IF (.NOT. BTEST( mzobbd(ilev,nbterr),nvrz )) kqcflag (ilev,2) =  3
        !   if only height flag is set, un-flag height obs for multi-level check
        IF (IBITS( nflagw,nvfzbp,nvfaoc ) == ISHFT( 1, nvfbps(4) ))            &
                                                     kqcflag (ilev,2) =  0
        !   do not use blacklisted height obs (even if all obs are passive)
        IF (      BTEST( nflagw,nvfzbp+nvfbps(2) ))  kqcflag (ilev,2) = -1
        !   do not use unreal (i.e. derived) z-p-obs
        IF (      BTEST( nflagw,nvfzbp+nvfbps(6) ))  kqcflag (ilev,2) = -1
        IF (      (BTEST( nflagw,nvfzbp+nvfbps(5) ))                           &
            .AND. (kobtyp == OT_AIREP))              kqcflag (ilev,2) = -1
      ENDIF
!     IF (     (zobbdy(ilev,nbtrh) > rmdich)                                   &
!         .AND.(.NOT. BTEST( mzobbd(ilev,nbterr),nvrq ))                       &
!        .AND. (lprmlc))                                                       &
!       PRINT '("kqcflag",A, F8.0,3I3)'  ,  ystid, zobbdy(ilev,nbtp)           &
!              , ibit1( nflagw,nvfqbp+nvfbps(4) )                              &
!              , ibit1( nflagw,nvfqbp+nvfbps(2) ), kqcflag(ilev,3)
    ENDDO

!-------------------------------------------------------------------------------
!  Section 1.2: Compute factor by which the humidity first guess error is
!               enhanced in the presence of an observed inversion/ stable layer:
!    - stable layer term: increases with lapse rate approaching zero
!                               and with increasing vertical distance over which
!                                                   this lapse rate is computed
!    - inversion term   : increases with increasing temperature inversion
!                               and with increasing positive lapse rate
!-------------------------------------------------------------------------------

    DO ilev = 1 , nlev
      fqerinv (ilev) = c1
    ENDDO

    klev = 0
    DO ilev = 1 , nlev
      IF (      (BTEST( mzobbd(ilev,nbterr),nvrt ))                            &
          .AND. (kqcflag(ilev,4) >= 0) .AND. (kqcflag(ilev,4) < 3)) THEN
!         .AND. (BTEST( mzobbd(ilev,nbterr),nvrq ))                            &
!         .AND. (zobbdy(ilev,nbtp  ) > 29999._wp)) THEN
        klev         = klev + 1
        klevq (klev) = ilev
      ENDIF
!     fqerinv (ilev) = c1
    ENDDO
    IF (klev > 1) THEN
      DO ilev = 1 , klev-1
        zf     = - r_g/r_d *c2                                                 &
                  /(zt_o  (klevq(ilev+1)       ) + zt_o  (klevq(ilev)       )) &
                  /(zobbdy(klevq(ilev+1),nbtlop) - zobbdy(klevq(ilev),nbtlop))
        zlapse = zf *(zt_o(klevq(ilev+1)       ) - zt_o  (klevq(ilev)       ))
        zdtdt  = MAX( c0 ,   zt_o(klevq(ilev+1))                               &
                          - (zt_o(klevq(ilev  )) + qcqlapl /zf) )
!       zlapdz = qcqladz *(c1 - c1/(c1 + c05* zdtdt))
        zlapdz =           c1 - c1/(c1 +      zdtdt)
!       zlapse = - r_g/r_d *c2                                                 &
!                 /(zt_o  (klevq(ilev+1),      ) + zt_o  (klevq(ilev),      )) &
!                 *(zt_o  (klevq(ilev+1),      ) - zt_o  (klevq(ilev),      )) &
!                 /(zobbdy(klevq(ilev+1),nbtlop) - zobbdy(klevq(ilev),nbtlop))
        zlapsq = MAX( c0, MIN( c1 , c1 - zlapse / qcqlapl ) )
        zlapiq = MAX( c0, MIN( c2 ,      zlapse / qcqlapi ) )
        zinvrs = MAX( c0, zt_o(klevq(ilev+1)) - zt_o(klevq(ilev)))
        zfqinv = c1
        IF (oblat < qcqlatl)  zfqinv = qcqfge2 / qcqfge1
        zfqinv =  zfqinv  +  zlapsq *(c1 + zlapdz) * qcqlapf                   &
                          +  zinvrs *(c1 + zlapiq) * qcqinvf
!       zfqinv =  zfqinv  +  zlapsq *qcqlapf  +  zinvrs *qcqinvf
        fqerinv (klevq(ilev  )) = MAX( fqerinv(klevq(ilev  )) , zfqinv )
        fqerinv (klevq(ilev+1)) = MAX( fqerinv(klevq(ilev+1)) , zfqinv )
        IF (lprmlc) THEN
          WRITE( nupr,'("finvq1 ",A, 2F8.0,2F7.2,3F6.3)' )                     &
                 ystid, zobbdy(klevq(ilev),nbtp), zobbdy(klevq(ilev+1),nbtp)   &
                      , zt_o  (klevq(ilev)     ), zt_o  (klevq(ilev+1)     )   &
                      , zdtdt, zlapdz, zlapiq
          WRITE( nupr,'("finvqer",A, 2F8.0,2I3,F9.5,5F6.3)' )                  &
                 ystid, zobbdy(klevq(ilev),nbtp), zobbdy(klevq(ilev+1),nbtp)   &
               , klevq(ilev), kqcflag(klevq(ilev),4), zlapse, zlapsq, zinvrs   &
                      , zfqinv, fqerinv(klevq(ilev)), fqerinv(klevq(ilev+1))
        ENDIF
      ENDDO
    ENDIF

!-------------------------------------------------------------------------------
!  Section 1.3: Perform QC for individual observations (first guess check)
!-------------------------------------------------------------------------------

    DO ilev = 1 , nlev

      CALL mult_quality_cntl ( zobbdy (ilev,:), vimtoob(ilev,1:5)              &
                             , fqerinv(ilev), kobtyp, tabdif, qcvf, qcc        &
                             , 4*nlev, lwrqc, lprqcm, nupr, ystid              &
                             , mzobbd(ilev,:), kqcflag(ilev,:), nlqc, zyqc )
!     ======================

    ENDDO
! ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Quality control: Multi-level check (for conventional(!) reports)
!             ----------------------------------
!             - If 4 subsequent mandatory-level data
!               or at least 50 % of all data are rejected in local checks,
!               then reject all data of that variable
!             - If at least 50 % of all data within 3 subsequent mandatory
!               levels are rejected, reject all data within this range
!             Note: The multi-level check is performed only for 'active' obs
!                   (in order to avoid checking active and e.g. blacklisted data
!                    together).
!                   Therefore, the flags must only be set for active obs
!-------------------------------------------------------------------------------

! IF ((lqc) .AND. (kobtyp /= OT_SATEM) .AND. (ktoplev+kbotlev < nlev)) THEN
  IF (ktoplev+kbotlev < nlev) THEN

    !   check for surface level, and for humidity use flag limit (at 300 hPa)
    ilvsfc = 0
    ilvpqm = nlev + 1
    DO ilev = nlev-ktoplev , 1+kbotlev , -1
      IF (BTEST( mzobbd(ilev,nbtlid),nvlidp(7)        ))                       &
        ilvsfc = ilev
      IF (BTEST( mzobbd(ilev,nbtflg),nvfqbp+nvfbps(4) ))                       &
        ilvpqm = ilev
    ENDDO

    DO ivar = 1 , 4

      ! check for 4 subsequent rejected mandatory-level data
      ! and 50% check of all data for short multi-level reports
      ! -------------------------------------------------------

      IF (ivar == 1)  nbtx = nbtu
      IF (ivar == 2)  nbtx = nbtz
      IF (ivar == 3)  nbtx = nbtt
      IF (ivar == 4)  nbtx = nbtrh
      ilvmlcb = nlev - ktoplev
      ilvmlct = 1 + kbotlev
      mandchk = 0
      mandrej = 0
      nact    = 0
      nrej    = 0
      DO ilev = 1+kbotlev , nlev-ktoplev
        IF (      (zobbdy(ilev,nbtx) > rmdich)                                 &
            .AND. (zobbdy(ilev,nbtp) > 9999.9_wp)) THEN
          nphpa = NINT( zobbdy(ilev,nbtp) * c100r )
          lmandat =    (nphpa == 1000) .OR. (nphpa == 850) .OR. (nphpa == 700) &
                  .OR. ((MOD( nphpa,100 ) == 0) .AND. (nphpa <= 500))          &
                  .OR. ((MOD( nphpa, 50 ) == 0) .AND. (nphpa <= 250))
          IF ((ilev == ilvsfc) .AND. (ivar /= 2))  lmandat = .FALSE.
          IF ((ilev >= ilvpqm) .AND. (ivar == 4))  lmandat = .FALSE.
          IF (lmandat) THEN
            mandchk = mandchk + 1
            ilvmand (mandchk) = ilev
          ENDIF
          IF (kqcflag(ilev,ivar) == 0) THEN
            IF ((lmandat) .AND. (mandrej < 4)) mandrej = 0
            nact = nact + 1
          ELSEIF (.NOT. (     ((ilev == ilvsfc) .AND. (ivar /= 2))             &
                         .OR. ((ilev >= ilvpqm) .AND. (ivar == 4)))) THEN
            IF (lmandat) mandrej = mandrej + 1
            nrej = nrej + 1
          ENDIF
        ENDIF
      ENDDO
      IF (     ((mandrej >= MIN( mandchk , 4 )) .AND. (mandrej > 1))           &
          .OR. ((nact <= nrej) .AND. (nrej >= 2))) THEN
        ilvmlcb = 1 + kbotlev
        ilvmlct = nlev - ktoplev
      ELSEIF (mandchk >= 4) THEN

      ! piecewise 50% check of all data for long multi-level reports
      ! ------------------------------------------------------------

        ilvmand (1)       = 1 + kbotlev
        ilvmand (mandchk) = nlev - ktoplev
        DO ilvq = 1 , mandchk - 2
          nact = 0
          nrej = 0
          DO ilev = ilvmand(ilvq) , ilvmand(ilvq+2)
            IF (kqcflag(ilev,ivar) >= 1) THEN
              IF (.NOT. (     ((ilev == ilvsfc) .AND. (ivar /= 2))             &
                         .OR. ((ilev >= ilvpqm) .AND. (ivar == 4))))           &
                nrej = nrej + 1
            ELSEIF (zobbdy(ilev,nbtx) > rmdich) THEN
              nact = nact + 1
            ENDIF
          ENDDO
          IF ((nact <= nrej) .AND. (nrej > 0)) THEN
            ilvmlcb = MIN( ilvmand(ilvq  ) , ilvmlcb )
            ilvmlct = MAX( ilvmand(ilvq+2) , ilvmlct )
          ENDIF
        ENDDO
      ENDIF

      !   determine 'kqcflag', 'kanaflg' so that old type of multi-level check
      !   keeps working  (used if (qcvf(4) < epsy))
      DO ilev = 1+kbotlev , nlev-ktoplev
        kqcflag (ilev,5) = 1
      ENDDO
      kanaflg (1) = 3

      ! -----------------------------
      ! New type of multi-level check:
      ! -----------------------------
      ! if each of 4 subsequent non-empty standard layers contains at least one
      ! observation with flag >= 1, then reject all data in all analysis layers
      ! which lie within these standard layers & contain at least 1 flagged obs
      ! -----------------------------------------------------------------------

      IF  (qcvf(4) > epsy) THEN
!     IF ((qcvf(4) > epsy) .AND. (ivar == 4)) THEN
        !  determine which analysis layers contain flagged data (kanaflg >=  1),
        !           only good data (kanaflg == 0) resp. no data (kanaflg == -1)
        DO ilayan = 1 , nqcalv - 1
          kanaflg (ilayan) = -1
          DO ilev = 1+kbotlev , nlev-ktoplev
            IF (      (zobbdy(ilev,nbtx) > rmdich)                             &
                .AND. (zobbdy(ilev,nbtp) <  tabanp(ilayan))                    &
                .AND. (zobbdy(ilev,nbtp) >= tabanp(ilayan+1))) THEN
              kanaflg (ilayan) = MAX( kanaflg(ilayan) , kqcflag(ilev,ivar) )
              !   (also determine which obs layer lies in which analysis layer)
              kqcflag (ilev,5) = ilayan
            ENDIF
          ENDDO
!         IF ((lprmlc) .AND. ((ivar == 4) .OR. (     (jo == jonl)              &
!                                               .AND.(kanaflg(ilayan) >= 0)))) &
          IF ((lprmlc) .AND. ((ivar == 4) .OR. (kanaflg(ilayan) >= 0)))        &
            WRITE( nupr,'("kanflg1  ",A ,F8.0, 3I3)' )                         &
                           ystid, tabanp(ilayan), ilayan, kanaflg(ilayan), ivar
        ENDDO

        !   multi-level check criterion: 4 consecutive standard layers
        !                                contain flags;
        !   find the interval between the lowest and the uppermost flagged
        !   analysis layers that lie within standard layers which satisfy
        !   the multi-level check criterion;
        !   determine whether there is a standard layer with only active obs
        !                                (nact == 2)
        nrej   = 0
        nact   = 0
        ilayb  = nqcalv + 1
        ilayt  = 0
        ilayan = 0
        DO ilvq = 1 , nqclev
          nrejstd = -1
          DO ilev = 1 , nanstd(ilvq)
            ilayan = ilayan + 1
            IF (kanaflg(ilayan) >= 1) THEN
              nrejstd  =  MAX( nrejstd, 0 ) + 1
              IF (nrej    == 0)  ilvb  = ilayan
              IF (nrejstd == 1)  nrej  = nrej + 1
              IF (nrej    == 4)  ilayb = MIN( ilayb, ilvb )
              IF (nrej    >= 4)  ilayt = ilayan
            ELSEIF (kanaflg(ilayan) == 0) THEN
              nrejstd  =  MAX( nrejstd, 0 )
              IF (nact    <  1)  nact  = 1
            ENDIF
          ENDDO
          IF (nrejstd == 0)  nrej  = 0
          IF ((nrej == 0) .AND. (nact == 1))  nact = 2
          IF ((nrej >= 1) .AND. (nact == 1))  nact = 0
          IF ((lprmlc) .AND. (ilvq <= 6) .AND. (ivar == 4))                    &
            WRITE( nupr,'("ilaybt ",A,5I3)') ystid,nrej,nact, ilayb,ilayt,ilayan
        ENDDO
!       IF ((lprmlc) .AND. ((ivar == 4) .OR. (jo == jonl)))                    &
        IF ((lprmlc) .AND. (ivar == 4))                                        &
          WRITE( nupr,'("ilaybt ",A,5I3)') ystid, nrej, nact, ilayb,ilayt,ilayan
        !   determine analysis layers to be rejected
        !   (if there is no active standard layer, the multi-level check
        !    criterion is assumed to apply to the whole profile)
        DO ilayan = 1 , nqcalv - 1
          IF (nact /= 0) THEN
            IF ((ilayan < ilayb) .OR. (ilayan > ilayt))                        &
              kanaflg (ilayan) = MIN( kanaflg(ilayan) , 0 )
          ENDIF
        ENDDO
        IF ((lprmlc) .AND. (ivar == 4)) THEN
          DO ilayan = 1 , nqcalv - 1
            IF (kanaflg(ilayan) >= 0) WRITE( nupr,'("kanflg2 ",A ,F8.0, 2I3)') &
                                             ystid, tabanp (ilayan), ilayan    &
                                                  , kanaflg(ilayan)
          ENDDO
        ENDIF
        !   determine interval with observations to be rejected
        ilvmlcb = nlev - ktoplev
        ilvmlct = 1 + kbotlev
        DO ilev = 1+kbotlev , nlev-ktoplev
          IF (kanaflg(kqcflag(ilev,5)) >= 1) THEN
            ilvmlct = ilev
            IF (ilvmlcb == nlev-ktoplev)  ilvmlcb = ilev
          ENDIF
        ENDDO
      ENDIF

      ! flagging of data not to be accepted by the multi-level checks
      ! -------------------------------------------------------------
      ! (multi-level check rejection of temperature implies to reject humidity,
      !  multi-level check rejection of height does not imply reject. of temp.)

      IF (ilvmlcb < ilvmlct) THEN
        IF (ivar == 1) nvrx = nvru
        IF (ivar == 2) nvrx = nvrz
        IF (ivar == 3) nvrx = nvrt
        IF (ivar == 4) nvrx = nvrq
        DO ilev = ilvmlcb , ilvmlct
          IF (kanaflg(kqcflag(ilev,5)) >= 1) THEN
            !   (by setting 'kqcflag=3' for passive obs in section 1.2:)
            !  the multi-level check is performed only for 'active' observations
            !  -----------------------------------------------------------------
            !  (in order to avoid checking active and e.g. blacklisted data
            !   together)
            !  therefore, the flags must only be set for active obs
            IF (      (      BTEST( mzobbd(ilev,nbterr),nvrx ))                &
                .AND. (.NOT. BTEST( mzobbd(ilev,nbtqcf),nvrx )))               &
              mzobbd (ilev,nbtqcf) = IBSET ( mzobbd(ilev,nbtqcf), nvrx )
            IF (      (ivar == 3)                                              &
                .AND. (      BTEST( mzobbd(ilev,nbterr),nvrq ))                &
                .AND. (.NOT. BTEST( mzobbd(ilev,nbtqcf),nvrq )))               &
              mzobbd (ilev,nbtqcf) = IBSET ( mzobbd(ilev,nbtqcf), nvrq )
          ENDIF
        ENDDO
        IF (lwrqc) THEN
          nlqc = nlqc + 1
          ivrs = ivar - MIN( ivar-1 , 1 )
          IF (ivar == 2)  ivrs = 4
          zyqc (nlqc,1) = REAL( 10 + ivrs , wp )
          zyqc (nlqc,2) = zobbdy(ilvmlcb,nbtp) * c100r
          zyqc (nlqc,5) = zobbdy(ilvmlct,nbtp) * c100r
        ENDIF
      ENDIF

    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 3: Re-set extrapolated model values to 'missing value'
!             (for the inverse observation operator, height and thickness QC,
!              and further processing in the nudging)
!-------------------------------------------------------------------------------

  !   for obs levels below the lowest main model level,
  !   set extrapolated model values to 'missing value'
  !   (so that they are not used in the nudging or the height QC)
  DO ilev = 1 , kbotlev
    DO ivar = 1 , 5
      vimtoob (ilev,ivar) = rmdi
    ENDDO
  ENDDO

  !   for obs levels above the top main model level,
  !   set extrapolated model values to 'missing value'
  DO ilev = nlev - ktoplev + 1 , nlev
    DO ivar = 1 , 5
      vimtoob (ilev,ivar) = rmdi
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! End of module procedure mult_obs_qc_fg
!-------------------------------------------------------------------------------

END SUBROUTINE mult_obs_qc_fg



!-------------------------------------------------------------------------------
!+ Module procedure for quality control of individual observations
!-------------------------------------------------------------------------------

SUBROUTINE mult_quality_cntl ( zobbdy1, vimx, fqerinv, kobtyp, tabdif          &
                             , qcvf, qcc, ndqc, lwrqc, lprqcm, nupr, ystid     &
                             , mzobbd1, kqcflag, nlqc, zyqc )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine performs the threshold quality control for one
!   observation level from one conventional multi-level report.
!   Rejection flags are updated, and simulated observation written to the
!   simulated observation record SOR, if required.
!   Rejection of observations at the observation time is documented.
!
! Method:
!   Checks of 'observation minus interpolated model value' against threshold.
!   The thresholds may depend on the pressure level of the observation point
!   (except for humidity) and on the difference between the observation time
!   and the current model time.
!
! Written by        :  Christoph Schraff, DWD  (original version: 13.08.98)
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
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    kobtyp     ,& ! observation type
    nupr       ,& ! file unit number for control output
    ndqc          ! dimension of 'zyqc'

  REAL (KIND=wp) , INTENT (IN)         ::       &
    vimx (5)   ,& ! model value inter-/extrapolated to the obs. level (pressure)
                  !   (1, 2: (u,v)-wind comp.  ;  3: temperature ;
                  !      4 : relative humidity ;  5: height)
    tabdif     ,& ! time distance between the obs. and the model state
    fqerinv    ,& ! factor by which humidity first guess error is enhanced
                  ! in the presence of an observed inversion / stable layer
    qcc  (4)   ,& ! vertically constant part of QC thresholds
    qcvf (4)      ! multiplication factor to the vertically varying part of the
                  !   QC thresholds (as defined in 'data_obs_qc_limits')

  LOGICAL        , INTENT (IN)          ::      &
    lprqcm (2) ,& ! printout for control at present station / node / time
    lwrqc         ! data with obs time close to current time, which are rejected
                  ! by the threshold quality control, are printed to a separate
                  ! file at the current timestep

! CHARACTER (LEN=ilstidp) , INTENT (IN)         ::       &
  CHARACTER (LEN=* )      , INTENT (IN)         ::       &
    ystid         ! obs. station identity to be printed

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy1 (mxrbdy)    ! 1 level of multi-level obs body (format: see 'omlbdy')

  INTEGER        , INTENT (INOUT)      ::       &
    mzobbd1 (mxrbdf) ,& ! 1 level of multi-level obs body (format: see 'momlbd')
    kqcflag (5)      ,& ! quality control flag
    nlqc                ! number of QC messages

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zyqc    (ndqc,11)   ! for QC messages

! Local parameters:
! ----------------
  REAL (KIND=wp) , PARAMETER  :: &
    c100r = 0.01_wp    ! 1/100

! Local variables (scalars or automatic arrays):
! ---------------------------------------------

! REAL (KIND=wp) ::  &
!   qcc (4) = c0, 500._wp, c0, c0   ! constant part of QC thresholds

  INTEGER        ::  &
    ilva   , ilvb    ,& ! indices of correlation scale table levels that are
                        ! adjacent to the specified model level
    mflgqc (4)          ! threshold quality control bit flag

  REAL (KIND=wp) ::  &
    qctu   , qctmu   ,& ! threshold for flag 3 / mlc  for wind vector incr.
    qctt   , qctmt   ,& ! threshold for flag 3 / mlc for temperature incr.
    qctrh  , qctmrh  ,& ! threshold for flag 3 / mlc for rel. humidity incr.
    qctz   , qctmz   ,& ! threshold for flag 3 / mlc for height increment
    uerror           ,& ! mean wind (vector) obs error at current pressure level
    terror           ,& ! mean temperature obs. error at current pressure level
    zerror           ,& ! vertical part of the height / thickness threshold
    zqermod          ,& ! first guess humidity error
    zqernor          ,& ! humidity obs error / tot error normalised by f.g error
    qerror1, qerror2 ,& ! threshold factor for flag 1, 2 to humidity increment
    qerror3          ,& ! threshold factor for flag 3 to humidity increment
    zlopob           ,& ! log( pressure ) at obs. point
    zf               ,& ! factor for vert. interpol. from q.c. thresholds table
    vecdif           ,& ! length of the obs. increment wind vector
    zobinc              ! scalar observation increment
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine mult_quality_cntl
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Determination of thresholds
!-------------------------------------------------------------------------------

  zlopob  = LOG( zobbdy1(nbtp) )

  !   vertically dependent part of thresholds
  IF ((zlopob >= tabqclp(1)) .OR. (zlopob <= tabqclp(nqclev))) THEN
    ilva = 1
    IF (zlopob <= tabqclp(nqclev)) ilva = nqclev
    ilvb = ilva
    zf   = c1
  ELSE
    ilva = 2
    DO WHILE (zlopob <= tabqclp(ilva))
      ilva = ilva + 1
    ENDDO
    ilvb = ilva - 1
    zf   = (tabqclp(ilvb) - zlopob) / (tabqclp(ilvb) - tabqclp(ilva))
  ENDIF
  IF (kobtyp == OT_AIREP) THEN
    uerror = qcvairf * ((c1-zf) *qcvairp(ilvb) + zf *qcvairp(ilva))
    terror = qctairf * ((c1-zf) *qctairp(ilvb) + zf *qctairp(ilva))
  ELSEIF (kobtyp == OT_PILOT) THEN
    uerror = qcvpilf * ((c1-zf) *qcvsond(ilvb) + zf *qcvsond(ilva))
    terror =            (c1-zf) *qctsond(ilvb) + zf *qctsond(ilva)
  ELSE
    uerror =            (c1-zf) *qcvsond(ilvb) + zf *qcvsond(ilva)
    terror =            (c1-zf) *qctsond(ilvb) + zf *qctsond(ilva)
  ENDIF
  zerror = (c1-zf) *qczsond(ilvb) + zf *qczsond(ilva)

  !   (max. RH-threshold is 70% when an obs is QC checked for the first time)
! qerror3 = 0.7_wp / (c1 + qctf(4) *wtukrsa)
  qerror3 = 0.7_wp
  qerror2 = qerror3  *  qcqfgef(2) / qcqfgef(3)
  qerror1 = qerror3  *  qcqfgef(1) / qcqfgef(3)
  IF (BTEST( mzobbd1(nbterr),nvrq )) THEN
    zqermod = qcqfge1 * fqerinv
    zqernor = ABS( zobbdy1(nbtqer) ) / zqermod
    zqernor = SQRT( MIN( zqernor *zqernor + c1 , 4._wp ) )
    qerror3 = MIN( zqernor * qcqfgef(3) * zqermod , qerror3 )
    qerror2 = MIN( zqernor * qcqfgef(2) * zqermod , qerror2 )
    qerror1 = MIN( zqernor * qcqfgef(1) * zqermod , qerror1 )
  ENDIF

  !   total thresholds
  qctu   = (c1 + qctf(1) *tabdif) * (qcc(1) + qcvf(1) *uerror)
  qctt   = (c1 + qctf(3) *tabdif) * (qcc(3) + qcvf(3) *terror)
  qctrh  = (c1 + qctf(4) *tabdif) * (qcc(4) + qcvf(4) *qerror3)
  qctz   = (c1 + qctf(3) *tabdif) * (         qcvf(2) *zerror *qcvfz(2))
  qctmu  = (c1 + qctf(1) *tabdif) * (qcc(1) + qcvf(1) *uerror)
  qctmt  = (c1 + qctf(3) *tabdif) * (qcc(3) + qcvf(3) *terror)
  qctmrh = (c1 + qctf(4) *tabdif) * (qcc(4) + qcvf(4) *qerror1)
! qctmrh = (c1 + qctf(4) *tabdif) * (qcc(4) + qcvf(4) *qerror2)
  qctmz  = (c1 + qctf(3) *tabdif) * (         qcvf(2) *zerror *qcvfz(2))

!-------------------------------------------------------------------------------
!  Section 2: Check against threshold of obs. increment of WIND vector
!-------------------------------------------------------------------------------

  mflgqc (1) = 0
  IF (zobbdy1(nbtu ) > rmdich) THEN
    vecdif     = SQRT(  (vimx(1) - zobbdy1(nbtu)) **2                          &
                      + (vimx(2) - zobbdy1(nbtv)) **2 )
    IF (     (vecdif > qctu)                                                   &
    !   reject zero aircraft wind if model wind speed exceeds 'qcvairl'
        .OR. (      (kobtyp == OT_AIREP) .AND. (vecdif > qcvairl)              &
              .AND. (MAX( ABS(zobbdy1(nbtu))                                   &
                        , ABS(zobbdy1(nbtv))) < epsy))) mflgqc (1) = 1
    !   reject aircraft wind if threshold QC rejects temperature
    IF ((kobtyp == OT_AIREP) .AND. (zobbdy1(nbtt) > rmdich)) THEN
      IF (ABS( vimx(3)-zobbdy1(nbtt) ) > qctt)    mflgqc (1) = 1
    ENDIF
    IF ((kqcflag(1) == 0) .AND. (vecdif > qctmu))  kqcflag (1) = 1
    IF ((kqcflag(1) >= 0) .AND. (mflgqc(1) == 1))  kqcflag (1) = 3
  ENDIF

  IF ((BTEST(mzobbd1(nbterr),nvru )) .AND. (mflgqc(1) == 1)) THEN
    IF ((lwrqc) .OR. (.NOT. BTEST( mzobbd1(nbtqcf),nvru ))) THEN
      nlqc = nlqc + 1
      zyqc (nlqc, 1) = 1.0_wp
      zyqc (nlqc, 2) = zobbdy1(nbtp) * c100r
      zyqc (nlqc, 5) = qctu
      zyqc (nlqc, 6) = zobbdy1(nbtu)
      zyqc (nlqc, 7) = vimx(1)
      zyqc (nlqc, 8) = zobbdy1(nbtv)
      zyqc (nlqc, 9) = vimx(2)
      zyqc (nlqc,10) = vecdif
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!  Section 3: Check against threshold of TEMPERATURE obs. increment
!-------------------------------------------------------------------------------

  mflgqc (2) = 0
  IF (zobbdy1(nbtt ) > rmdich) THEN
    zobinc     =  ABS( vimx(3) - zobbdy1(nbtt) )
    !   (reject aircraft temperature if aircraft wind is rejected)
    IF ((zobinc > qctt) .OR. ((kobtyp == OT_AIREP) .AND. (mflgqc(1) == 1)))    &
      mflgqc (2) = 1
    IF ((kqcflag(3) == 0) .AND. (zobinc > qctmt))  kqcflag (3) = 1
    IF ((kqcflag(3) >= 0) .AND. (mflgqc(2) == 1))  kqcflag (3) = 3
    IF ((kobtyp == OT_AIREP) .AND. (zobbdy1(nbtu) > c0)) THEN
      IF (kqcflag(1) >= 0)  kqcflag (1) = MAX( kqcflag(1) , kqcflag(3) )
      IF (kqcflag(3) >= 0)  kqcflag (3) = MAX( kqcflag(1) , kqcflag(3) )
    ENDIF
  ENDIF

  IF ((BTEST(mzobbd1(nbterr),nvrt )) .AND. (mflgqc(2) == 1)) THEN
    IF ((lwrqc) .OR. (.NOT. BTEST( mzobbd1(nbtqcf),nvrt ))) THEN
      nlqc = nlqc + 1
      zyqc (nlqc, 1) = 3.0_wp
      zyqc (nlqc, 2) = zobbdy1(nbtp) * c100r
      zyqc (nlqc, 5) = qctt
      zyqc (nlqc, 6) = zobbdy1(nbtt)
      zyqc (nlqc, 7) = vimx(3)
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!  Section 4: Check against threshold of RELATIVE HUMIDITY obs. increment
!-------------------------------------------------------------------------------

  mflgqc (3) = 0
  IF (zobbdy1(nbtrh) > rmdich) THEN
    zobinc     =  ABS( MIN(vimx(4),c1) - zobbdy1(nbtrh) )
    !   (reject humidity if temperature is rejected)
    IF ((zobinc > qctrh) .OR. (mflgqc(2) == 1))  mflgqc (3) = 1
    IF ((kqcflag(4) == 0) .AND. (zobinc > qctmrh))  kqcflag (4) = 1
    IF ((kqcflag(4) >= 0) .AND. (mflgqc(3) == 1))   kqcflag (4) = 3
    IF ((lprqcm(1)) .AND. ((lprqcm(2)) .OR. (mflgqc(3) == 1)))                 &
      WRITE( nupr,'("kmlqc  ",A ,F8.0, I3, 7F7.3)' )                           &
             ystid, zobbdy1(nbtp), kqcflag(4), zobbdy1(nbtrh)                  &
                  , vimx(4), qctrh, qctmrh, fqerinv, zqernor, qerror3
  ENDIF

  IF ((BTEST(mzobbd1(nbterr),nvrq )) .AND. (mflgqc(3) == 1)) THEN
    IF ((lwrqc) .OR. (      (.NOT. BTEST( mzobbd1(nbtqcf),nvrq ))              &
                      .AND. (zobinc > qctrh))) THEN
      nlqc = nlqc + 1
      zyqc (nlqc, 1) = 4.0_wp
      zyqc (nlqc, 2) = zobbdy1(nbtp) * c100r
      zyqc (nlqc, 5) = qctrh
      zyqc (nlqc, 6) = zobbdy1(nbtrh)
      zyqc (nlqc, 7) = vimx(4)
      IF (zobinc <= qctrh) THEN
        zyqc (nlqc, 8) = zobbdy1(nbtt)
        zyqc (nlqc, 9) = vimx(3)
        zyqc (nlqc, 1) = - 4.0_wp
      ENDIF
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!  Section 5: Check against threshold of HEIGHT obs. increment
!-------------------------------------------------------------------------------

  mflgqc (4) = 0
  IF (zobbdy1(nbtz ) > rmdich) THEN
    zobinc     =  ABS( vimx(5) - zobbdy1(nbtz) )
    IF (zobinc > qctz)  mflgqc (4) = 1
    IF ((kqcflag(2) == 0) .AND. (zobinc > qctmz))  kqcflag (2) = 1
    IF ((kqcflag(2) >= 0) .AND. (mflgqc(4) == 1))  kqcflag (2) = 3
  ENDIF

!-------------------------------------------------------------------------------
!  Section 6: Set threshold quality control flags in ODR
!-------------------------------------------------------------------------------

  CALL MVBITS ( mflgqc(1), 0, 1, mzobbd1(nbtqcf), nvru )
  CALL MVBITS ( mflgqc(2), 0, 1, mzobbd1(nbtqcf), nvrt )
  CALL MVBITS ( mflgqc(3), 0, 1, mzobbd1(nbtqcf), nvrq )
  CALL MVBITS ( mflgqc(4), 0, 1, mzobbd1(nbtqcf), nvrz )

! IF (lveridat) PRINT *,'QCI5 ',ystid, zobbdy1(nbtp), mzobbd1(nbtqcf), mflgqc

!-------------------------------------------------------------------------------
! End of module procedure mult_quality_cntl
!-------------------------------------------------------------------------------

END SUBROUTINE mult_quality_cntl



!-------------------------------------------------------------------------------
!+ Module procedure for height / thickness QC for multi-level temperature
!-------------------------------------------------------------------------------

SUBROUTINE mult_obs_qc_dz ( nlev, zobbdy, zdzob, tabdif, qcvfz_in, r_d, r_g    &
                          , mzobbd , zyqc , lretv_in , lrejtot, lredo_out )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine performs a quality control check of height and thickness
!   derived hydrostatically from multi-level temperature observations.
!
! Method:
!   For 'long' T-obs profiles (which exceed 3 consecutive standard levels),
!     1) first do a check of hydrostatic height:
!        - if an increment exceeds large threshold, then set lower bound of
!          the levels to be rejected
!        - if any increment exceeds small threshold, then do thickness checks
!          for whole profile
!     2) thickness check (also done for short profiles):
!        - interpolate height increments to standard levels
!        - thickness between neighbouring standard level in the low and middle
!          troposphere and between 3 standard level further above, using height
!          error correlations in the definition of the thresholds.
!   For definition of thresholds, see comments in Section 4.
!   For very short temperature profiles (vertical extent <= 50 hPa), this
!   routine is not called, so that the height / thickness check is not done.
!
! Written by        :  Christoph Schraff, DWD  (original version: 2000)
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
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    nlev                ! number of vertical levels in report

  REAL (KIND=wp) , INTENT (IN)         ::       &
    r_d              ,& ! gas constant for dry air       (287.05)
    r_g                 ! acceleration due to gravity    (  9.80665)

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zdzob  (nlev)        ,& ! profile of height observation increments
    zobbdy (nlev,mxrbdy) ,& ! multi-level obs. body (format: see 'omlbdy')
    tabdif               ,& ! time distance between the obs. and the model state
    qcvfz_in                ! multiplication factor to the vertically varying
                            !   part of the QC threshold for z
                            !   (as defined in 'data_obs_qc_limits')

  INTEGER        , INTENT (INOUT)      ::       &
    mzobbd  (nlev,mxrbdf)   ! multi-level obs. body (format: see 'momlbd')

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    zyqc    (11)            ! information for QC control messages

  LOGICAL        , INTENT (IN)    , OPTIONAL   ::       &
    lretv_in                ! if .TRUE. then obs type = satellite retrieval;
                            ! --> if 'lretv_in' does not exist then assume
                            !     conventional report (not satellite retrieval)

  LOGICAL        , INTENT (OUT)   , OPTIONAL   ::       &
    lrejtot              ,& ! reject total (temperature) profile due to (d)z-QC
    lredo_out               ! if .TRUE. then re-do the interpolation of observed
                            !   temperature and humidity to model levels later
                            !   on in the nudging
                            !  (since a part of the profile has been set passive
                            !   in the height / thickness check of this routine)

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    c100r   =   0.01_wp     !

! Local scalars:
! -------------

  INTEGER        ::  &
    ilev             ,& ! vertical loop index over observation levels
    ilva   , ilvb    ,& ! indices of obs. levels containing data of type 'ivrs'
                        ! that are adjacent to the specified model level
    ilvq             ,& ! loop index over standard pressure levels
    ilvqb  , ilvqt   ,& ! extent of profile (in standard pressure level units)
    ilevb  , ilevt   ,& ! extent of erroneous thickness (in obs. level units)
    ilveb  , ilvet   ,& ! extent of erroneous thickness (in std. level units)
    ilvactb, ilvactt ,& ! extent of active T-profile (in obs. level units)
    idlvq            ,& ! (std.) level increment for which thickness is checked
    mtypeqc             ! type (height or thickness) of quality cntl. not passed

  REAL (KIND=wp) ::  &
    zlopf            ,& ! weight factor for vertical interpol. to the obs. point
    zlperrb, zlperrt ,& ! vertical range (in LOG(p)) of erroneous thickness
    zlpertb, zlpertt ,& ! vertical range (in LOG(p)) of erroneous temperature
    qctzt            ,& ! temporal part of the height / thickness threshold
    zqqcvfz          ,& ! factor to increase thresholds above 100 hPa level
    zerror           ,& ! vertical part of the height / thickness threshold
    zcorlz           ,& ! height error correlation (between 2 obs. points)
    zcorlbb, zcorlbt ,& ! height error correlations (to base of report)
    zerrzb , zerrzt  ,& ! variances of thickness (betw: rep. base , 1 other pt)
    zvardz           ,& ! variance (mse) of thickness betw. top and base of rep.
    zthrdz           ,& ! threshold for thickness error
    zdzobs           ,& ! thickness increment over the erroneous vertical part
    zerrz            ,& ! height error threshold at lowest erroneous height obs.
    zerrb  , zerrt   ,& ! height error threshold at top and base of profile
    zfb    , zft        ! weights for interpol. from std. levels to obs. pts.

  LOGICAL        ::  &
    lredo            ,& ! if .TRUE. then re-do the interpolation of observed
                        !   temperature and humidity to model levels later on
    lretv            ,& ! if .TRUE. then obs type = satellite retrieval
!   lqciwv           ,& ! use obs for IWV spatial consistency checking
    lerrtot          ,& ! thickness between base and top of report erroneous
    lqcdz               ! perform thickness quality control

! Local (automatic) arrays:
! -------------------------

  REAL (KIND=wp) ::  &
    zlopob  (nlev)   ,& ! LOG( pressure ) of observation levels
    zdzlq   (nqclev) ,& ! profilght observation increments \   at mandatory
    zlplq   (nqclev) ,& ! LOG( pressure )                   \  pressure levels
    zerrlq  (nqclev) ,& ! height error threshold             > and at base and
    zcorlq  (nqclev)    ! height error correlation betw. 2  /  top of profile
                        !   levels with 1 level in between /

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine mult_obs_qc_dz
!-------------------------------------------------------------------------------

  lredo    = .FALSE.
  lerrtot  = .FALSE.
  zyqc (:) =  c0

  !   check for existence of 'lretv_in'
  IF (      PRESENT( lretv_in ))  lretv = lretv_in
  IF (.NOT. PRESENT( lretv_in ))  lretv = .FALSE.

! IF (lqcdz) THEN
    lqcdz   = .TRUE.

!-------------------------------------------------------------------------------
!  Section 1: Height check diagnosis
!-------------------------------------------------------------------------------

    !   auxiliary quantity
    DO ilev = 1 , nlev
      zlopob (ilev) = zobbdy(ilev,nbtlop)
    ENDDO

    zlperrb = zlopob(nlev)
    zlperrt = zlopob(1) + c1            !  + c1 is for safety
    ilevb   = nlev + 1
    ilevt   = 0

    !   determine extent of active T-profile (in obs levels)
    ilvactb = nlev
    ilvactt = 0
    DO ilev = 1 , nlev
      IF (      (      BTEST( mzobbd(ilev,nbterr),nvrt ))                      &
          .AND. (.NOT. BTEST( mzobbd(ilev,nbtqcf),nvrt ))) THEN
        IF (ilvactb == nlev) ilvactb = ilev
        ilvactt = ilev
      ENDIF
    ENDDO

! for short temperature profiles, i.e. if 3 consecutive mandatory
! levels cover the vertical extent of the multi-level report:
! one quality control check of the thickness of the whole report
! ---------------------------------------------------------------

    qctzt  =  c1 + qctf(3) *tabdif
!   qctzt  =  qctzt * qcvf(2)
    qctzt  =  qctzt * qcvfz_in

    !   determine mandatory level just below base and top of profile
    ilvqb = 1
    ilvqt = 1
    DO ilvq = 2 , nqclev - 1
      IF (tabqclp(ilvq) >= zlopob(ilvactb)-epsy) ilvqb = ilvq
      IF (tabqclp(ilvq) >  zlopob(ilvactt)+epsy) ilvqt = ilvq
    ENDDO

    !   prepare computation of height error correlation between base and top
    !   of report   (these factors are re-used further below)
    zfb  =  (tabqclp(ilvqb) - zlopob(ilvactb))                                 &
          / (tabqclp(ilvqb) - tabqclp(ilvqb+1))
    zft  =  (tabqclp(ilvqt) - zlopob(ilvactt))                                 &
          / (tabqclp(ilvqt) - tabqclp(ilvqt+1))
    IF (zlopob(ilvactb) >= tabqclp(1))       zfb = c0
    IF (zlopob(ilvactt) <= tabqclp(nqclev))  zft = c1

    !   check if 3 mandatory levels cover the vertical extent, and
    !   compute height error correlation between base and top of report
    IF (ilvqt-ilvqb <= 1) THEN
      IF (ilvqb == ilvqt) THEN
        zcorlz  = c1 - (zft - zfb) * (c1 - qczcorl(ilvqb,ilvqt+1))
      ELSE
        zcorlbb = (c1-zft) *qczcorl(ilvqb,ilvqt) + zft *qczcorl(ilvqb  ,ilvqt+1)
        zcorlbt = (c1-zft) *c1                   + zft *qczcorl(ilvqb+1,ilvqt+1)
        zcorlz  = (c1-zfb) *zcorlbb              + zfb *zcorlbt
      ENDIF

      !   compute mean square error (variance) of thickness between base and top
      !   by:
      !     var(dz(b,t),obs)  =  var(z(b))      +  var(z(t))
      !                        - 2* corr(z(b),z(t)) *SQRT( var(z(b))*var(z(t)) )
      !                       =  qczsond(b)**2  +  qczsond(t)**2
      !                        - 2* zcorlz          *( qczsond(b) * qczsond(t) )
      zerrzb = (c1-zfb) *qczsond(ilvqb) + zfb *qczsond(ilvqb+1)
      zerrzt = (c1-zft) *qczsond(ilvqt) + zft *qczsond(ilvqt+1)
      zvardz = zerrzb**2 + zerrzt**2 - 2* zcorlz *zerrzb *zerrzt

      !   compute threshold for thickness error between base and top by:
      !     thresh(dz(b,t))
      !       =  SQRT( var(dz,model)**2  +  var(dz,obs)**2 ) * errlim(z)
      !       =  var(dz,obs) * (SQRT( 1 + var(obs)**2/var(mod)**2 ) * errlim(z))
      !     where:  var(dz,model) / var(dz,obs) = x ;   x = 1 or SQRT(2)
      !             errlim(z) = threshold / variance             -------

      zthrdz = SQRT( MAX( zvardz , c0 ) ) * qcvfz(1) * qctzt

      !   prepare flagging the whole report if thickness quality control
      !   is not passed
      IF (     (ABS( zdzob(ilvactt) - zdzob(ilvactb) ) > zthrdz)               &
          .OR. (zthrdz < epsy)) THEN
        zlperrb = zlopob(ilvactb)
        zlperrt = zlopob(ilvactt)
        zdzobs  = zdzob(ilvactt) - zdzob(ilvactb)
        lerrtot = .TRUE.
        mtypeqc = 10
      ENDIF
    ENDIF

! for long temperature profiles:
! hydrostatic threshold quality control of height
! (using also obs. increments from vertical spreading)
! ----------------------------------------------------

    IF (ilvqt-ilvqb >= 2) THEN
      !   get total height QC thresholds
      lqcdz   = .FALSE.
      DO ilev = ilvactb , ilvactt
        IF (      (      BTEST( mzobbd(ilev,nbterr),nvrt ))                    &
            .AND. (.NOT. BTEST( mzobbd(ilev,nbtqcf),nvrt ))) THEN
          IF (zlopob(ilev) > tabqclp(1)) THEN
            zerror = qczsond(1)
          ELSEIF (zlopob(ilev) <= tabqclp(nqclev)) THEN
            zerror = qczsond(nqclev)
          ELSE
            ilva = 2
            DO WHILE (zlopob(ilev) <= tabqclp(ilva))
              ilva = ilva + 1
            ENDDO
            ilvb  = ilva - 1
            zlopf = (tabqclp(ilvb) -zlopob(ilev)) /(tabqclp(ilvb)-tabqclp(ilva))
            zerror = (c1-zlopf) *qczsond(ilvb) + zlopf *qczsond(ilva)
          ENDIF
          zerrt  =  zerror * qcvfz(2) * qctzt
          IF (ilev == ilvactb)  zerrb = zerrt

          !   quality control of height (i.e. for the complete profile):
          !   if rejected, determine 'zlperrb', 'ilevb'
          IF ((ilevb > nlev) .AND. (ABS( zdzob(ilev) ) > zerrt)) THEN
            ilevb = ilev
            zerrz = zerrt
          ENDIF
          IF (ABS( zdzob(ilev) ) > zerrt)  ilevt = ilev
          zerrt  =  zerror * qcvfz(1) * qctzt
          IF (ABS( zdzob(ilev) ) > zerrt)  lqcdz = .TRUE.
        ENDIF
      ENDDO
      zlperrb = zlopob(MIN( ilevb,nlev ))
      !     (levels for which zdzob > zerrt (large threshold) are rejected
      !      by the height check and not set active again)
!     IF (ilevt > 1)  zlperrt = zlopob(ilevt)
    ENDIF
    !  at this point, 'ilevb' resp. 'zlperrb' denotes the lower bound of the
    !  range of levels which are rejected by the height check
    !  (if the above line had not been commented out then
    !   [ilevb,ilevt] resp. [zlperrb,zlperrt] would have denoted
    !   the range of levels which are rejected by the height check)

!-------------------------------------------------------------------------------
!  Section 2: Thickness check diagnosis
!-------------------------------------------------------------------------------

    IF ((ilvqt-ilvqb >= 2) .AND. (lqcdz)) THEN

! for long temperature profiles with erroneous height:
! hydrostatic quality control of thickness,
! diagnosing approx. vertical extent of erroneous thickness
! ---------------------------------------------------------
! compute at mandatory levels within the profile and at base and top of profile:
! - height observation increments
! - LOG( pressure )
! - height error threshold
! - height error correlation between two levels with one level in between
!   (except level pairs 850-700 and 700-500 hPa, <-- idlvq=1 )

      zqqcvfz = qcvfz(2) / qcvfz(1)
      ilev = ilvactb
      DO ilvq = ilvqb + 1 , ilvqt
        DO WHILE (zlopob(ilev+1) > tabqclp(ilvq))
          ilev = ilev + 1
        ENDDO
        zlopf  = (zlopob(ilev) - tabqclp(ilvq)) / (zlopob(ilev) -zlopob(ilev+1))
        zdzlq (ilvq) = (c1-zlopf) *zdzob(ilev) + zlopf *zdzob(ilev+1)
        zlplq (ilvq) = tabqclp(ilvq)
        zerrlq(ilvq) = qczsond(ilvq) * qcvfz(1) * qctzt
        IF (zlplq(ilvq) < 9.21_wp) zerrlq(ilvq) = zerrlq(ilvq) * zqqcvfz
        IF (ilvq <= ilvqt-2) THEN
          zcorlq(ilvq) = qczcorl(ilvq,ilvq+2)
          IF ((ilvq == 2) .OR. (ilvq == 3))  zcorlq(ilvq) = qczcorl(ilvq,ilvq+1)
        ENDIF
      ENDDO
      zdzlq (ilvqb)   = zdzob (ilvactb)
      zlplq (ilvqb)   = zlopob(ilvactb)
      zerrlq(ilvqb)   = zerrb
      zcorlq(ilvqb)   = (c1-zfb) *qczcorl(ilvqb  ,ilvqb+2)                     &
                        +   zfb  *qczcorl(ilvqb+1,ilvqb+2)
      zdzlq (ilvqt+1) = zdzob (ilvactt)
      zlplq (ilvqt+1) = zlopob(ilvactt)
      zerrlq(ilvqt+1) = zerrt
      IF (zlplq(ilvqt+1) < 9.21_wp) zerrlq(ilvqt+1) = zerrt * zqqcvfz
      zcorlq(ilvqt-1) = (c1-zft) *qczcorl(ilvqt-1,ilvqt  )                     &
                        +   zft  *qczcorl(ilvqt-1,ilvqt+1)
      ilvet = 0

! compute threshold for thickness error between two levels (with one in between)
! [by: var(dz(b,t),obs)  =  var(z(b))      +  var(z(t))
!                         - 2* corr(z(b),z(t)) *SQRT( var(z(b))*var(z(t)) )
      DO ilvq = ilvqb , ilvqt - 1
        idlvq = 2
        IF (      ((ilvq == 2) .OR. (ilvq == 3))                               &
            .AND. (ilvq /= ilvqb) .AND. (ilvq /= ilvqt-1))  idlvq = 1
        zthrdz  =  zerrlq(ilvq) **2  +  zerrlq(ilvq+idlvq) **2                 &
                 - 2 * zcorlq(ilvq) * zerrlq(ilvq) * zerrlq(ilvq+idlvq)

        !   identify erroneous part in terms of mandatory levels and LOG(press.)
        IF (( zdzlq(ilvq+idlvq)-zdzlq(ilvq) )**2 > zthrdz) THEN
          IF (ilvet == 0) ilveb = ilvq
          ilvet = MAX( ilvet , ilvq+idlvq )
        ENDIF
      ENDDO
      IF (ilvet > 0) THEN
        lerrtot = (ilveb == ilvqb) .AND. (ilvet == ilvqt+1)
        zlperrb = MAX( zlperrb , zlplq(ilveb) )
        zlperrt = MIN( zlperrt , zlplq(ilvet) )

        !   prepare output for control
        IF ((ilveb == ilvqb) .AND. (ilvet == ilvqt+1)) THEN
          zcorlbb= (c1-zft)*qczcorl(ilvqb  ,ilvqt) +zft*qczcorl(ilvqb  ,ilvqt+1)
          zcorlbt= (c1-zft)*qczcorl(ilvqb+1,ilvqt) +zft*qczcorl(ilvqb+1,ilvqt+1)
          zcorlz = (c1-zfb) *zcorlbb              + zfb *zcorlbt
        ELSEIF (ilvet == ilvqt+1) THEN
          zcorlz = (c1-zft) *qczcorl(ilveb,ilvqt) + zft *qczcorl(ilveb,ilvqt+1)
        ELSEIF (ilveb == ilvqb  ) THEN
          zcorlz = (c1-zfb) *qczcorl(ilvqb,ilvet) + zfb *qczcorl(ilvqb+1,ilvet)
        ELSE
          zcorlz = qczcorl(ilveb,ilvet)
        ENDIF
        zthrdz  =  zerrlq(ilveb) **2  +  zerrlq(ilvet) **2                     &
                 - 2 * zcorlz * zerrlq(ilveb) * zerrlq(ilvet)
        zthrdz  = SQRT( zthrdz )
        zdzobs  = zdzlq(ilvet) - zdzlq(ilveb)
        mtypeqc = 10
      ELSE
        lerrtot = (ilevb == ilvactb)
        IF (ilevb <= nlev) THEN
          zlperrt = zlopob(ilvactt)
          zthrdz  = zerrz
          zdzobs  = zdzob(ilevb)
          mtypeqc = 9
        ENDIF
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
!  Section 3: Flagging of the part diagnosed as erroneous
!-------------------------------------------------------------------------------

! flagging of the temperature observations within the
! diagnosed erroneous vertical part of the profile, and
! setting quality control flags of temperature and height
! -------------------------------------------------------

    IF (zlperrb >= zlperrt-epsy) THEN
      zlpertt  =  zlperrt
      zlpertb  =  zlperrb
      IF (zlperrt-zlopob(ilvactt) <=  epsy)  zlpertt  =  zlopob(nlev)
      IF (zlperrb-zlopob(ilvactb) >= -epsy)  zlpertb  =  zlopob(1)
      DO ilev = 1 , nlev
        IF ((lretv) .OR. (      (zlopob(ilev) <= zlpertb+epsy)                 &
                          .AND. (zlopob(ilev) >= zlpertt-epsy))) THEN
          IF (zobbdy(ilev,nbtt ) > rmdich)                                     &
            mzobbd (ilev,nbtqcf) = IBSET ( mzobbd(ilev,nbtqcf), nvrt )
          IF (zobbdy(ilev,nbtrh) > rmdich)                                     &
            mzobbd (ilev,nbtqcf) = IBSET ( mzobbd(ilev,nbtqcf), nvrq )
          IF (zobbdy(ilev,nbtz ) > rmdich)                                     &
            mzobbd (ilev,nbtqcf) = IBSET ( mzobbd(ilev,nbtqcf), nvrz )
          !   for satellite retrievals, the complete profiles of temperature and
          !                             humidity are flagged
          !                       (Note that once they are flagged, they cannot
          !                             become unflagged later on again)
          IF (lretv)  lerrtot = .TRUE.
        ENDIF
      ENDDO

! printout message for control
!     IF ((lwrqc) .AND. (ntotqc < maxqcp)) THEN
!       ntotqc = ntotqc + 1
!       yyqc (ntotqc  ) = ystid
!       myqc (ntotqc,1) = kcdtyp
!       myqc (ntotqc,2) = mtypeqc
!       oyqc (ntotqc,1) = ztimob
!       oyqc (ntotqc,2) = EXP( zlperrb ) * c100r
!       oyqc (ntotqc,3) = zjlat
!       oyqc (ntotqc,4) = zilon
!       oyqc (ntotqc,5) = zthrdz
!       oyqc (ntotqc,6) = EXP( zlperrt ) * c100r
!       oyqc (ntotqc,7) = zdzobs
!       IF (mtypeqc == 10)                                                     &
!         oyqc (ntotqc,8) = zdzobs *r_g /(r_d *(zlperrb-zlperrt))
!     ENDIF
!     IF (lwrqc) THEN
!       nlqc = nlqc + 1
        zyqc (1) = REAL( mtypeqc , wp )
        zyqc (2) = EXP( zlperrb ) * c100r
        zyqc (5) = zthrdz
        zyqc (6) = EXP( zlperrt ) * c100r
        zyqc (7) = zdzobs
        IF (mtypeqc == 10) zyqc (8) = zdzobs *r_g /(r_d *(zlperrb-zlperrt))
!     ENDIF

! required for nudging only:
! prepare to indicate missing interpolated observed levels, and
! check whether 'mult_obs_2_modlev' needs to be called again (with iv1 == 2)
! --------------------------------------------------------------------------

      IF (lerrtot) THEN
!       DO mlev = ke , 1 , -1
!         viobtom (mlev,3) = rmdi
!         viobtom (mlev,4) = rmdi
!       ENDDO
!       mbotlv (2) = ke - mtoplv(2)
!       mbotlv (3) = ke - mtoplv(3)
      ELSE
        lredo  = .TRUE.
      ENDIF

    ENDIF
! ENDIF

  !   set 'lredo_out' and 'lrejtot' if they exist
  IF (PRESENT( lredo_out ))  lredo_out = lredo
  IF (PRESENT( lrejtot   ))  lrejtot   = lerrtot

!-------------------------------------------------------------------------------
! End of module procedure mult_obs_qc_dz
!-------------------------------------------------------------------------------

END SUBROUTINE mult_obs_qc_dz

!-------------------------------------------------------------------------------

END MODULE src_obs_qc_conv
