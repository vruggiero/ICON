!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_cdfin_print

!-------------------------------------------------------------------------------
! Description:
!   This module contains routines which write different types of information
!   onto different formatted ASCII files for diagnostic purposes of the
!   observation pre-processing.
!   Specifically, the following type of information is written:
!    - statistics on processed / active / passive / rejected reports
!    - statistics on report and data events which are the reason the reject
!      certain data
!    - rejection messages for each rejected report
!    - major parts of the reports stored in the ODR (Observation Data Record)
!    - number of reports on each sub-domain (node)
!
!   Note: This module belongs to a group of COSMO-related modules for reading
!   ----  conventional data from ('cdfin') NetCDF observation input files.
!         It is shared between the COSMO model and the DACE program package !
!
! Method:
!   This module contains the following module procedures:
!    - obs_cdf_print_statist : called by obs_org_cdf_proc
!    - obs_cdf_print_events  : called by obs_org_cdf_proc
!    - obs_cdf_print_reject  : called by obs_org_cdf_proc
!    - obs_cdf_print_odr     : called by obs_org_cdf_proc
!    - obs_print_number_rep  : called by obs_org_cdf_proc
!
!   This module also contains elemental functions, formerly statement functions:
!   - ibit1        : returns 1 bit at given bit position of given integer word
!
!   It uses from:
!    - parallel_utilities:    - global_values
!                             - gather_values
!    - environment:           - model_abort
!
!   Data modules used:
!    - kind_parameters
!    - data_obs_lib_cosmo
!    - data_obs_cdfin
!    - data_obs_record
!
! Current Code Owner (for COSMO and for DACE):
!  DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_22        2012-01-31 Christoph Schraff
!  Initial release, extracted from module 'src_obs_proc_cdf' and adapted (e.g.
!  modified routine interfaces and diagnostic arrays, 'obs_pointrs' --> 'i_cma',
!  on demand open and close of files, routines re-named and split up).
! V4_28        2013/07/12 Christoph Schraff
!  Format in yuobsdr modified (surface pressure error in Pa, not m (height)).
!  Statement functions replaced by elemental function 'ibit1'.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  Dimension of 'nevent' reduced from 3 to 2.
! V5_4d        2016-12-12 Christoph Schraff
!  Technical Adaptations
! @VERSION@    @DATE@     Christoph Schraff
!  - Set subroutines used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - (For ICON:) Number of digits used for printing the first grid point index
!    increased from 3 to 7.
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

USE data_obs_lib_cosmo, ONLY :  &

! 1. General parameters
! ---------------------

    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 2. Variables and parameters obtained by calling 'obs_cdf_interface'
! -------------------------------------------------------------------

    maxmlv     ,& ! size (level  dimension) of the  multi-level (m-l)  ODR

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    yustats    ,& ! statistics of processed reports
    yurejct    ,& ! direct reporting of rejected obs. reports
    yuobsdr    ,& ! observations stored in the observation data record
    nustat     ,& ! statistics of processed reports
    nurej      ,& ! direct reporting of rejected obs. reports
    nuodr      ,& ! observations stored in the observation data record
    nupr       ,& ! all the remaining information
    lopen_odr  ,& ! .true. if yuobsdr is open
    lopen_rej  ,& ! .true. if yurejct is open

! 6. Data type with rules for CMA obs and code types
! --------------------------------------------------

    n_cma      ,& ! number of CMA obs and code types
    t_cmatyp   ,& ! data type for information on CMA observation and code types
    cma        ,& ! array of meta data on CMA observation and code types

! 7. Functions
! ------------

    i_cma         ! function to determine the index of 'cma'
                  ! referring to a given CMA observation and code type

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

!         3.1     Format of event counters
!                 ------------------------
    mxreve     ,& ! length of report event counter array
    mxdeve     ,& ! length of data event counter array

!         3.3    Character descriptions of events and flags
!                ------------------------------------------
    crepev     ,& ! description of report events
    cdatev     ,& ! description of data events

! Section 7 :  For reporting rejection of data: Output buffer, size and formats
!-------------------------------------------------------------------------------

    outbuf     ,& ! buffer containing output for a single node
    nacout     ,& ! actual number of records stored in the output buffer
    nmxoln     ,& ! maximum length of output buffer
    istrej     ,& ! length of strings (station id) in output buffer
    nfmt1      ,& ! no pressure
    nfmt2      ,& ! excess of precipitation
    nfmt3      ,& ! no accepted data
    nfmt4      ,& ! excess of levels
    nfmt5      ,& ! several surface levels
    nfmt6      ,& ! excess of pressure tendency
    nfmt7      ,& ! excess of lapse rate
    nfmt8      ,& ! excess of wind speed shear
    nfmt9      ,& ! excess of directional shear
    nfmt10     ,& ! redundancy of surface-level report
    nfmt11     ,& ! redundancy of multi-level report
    nfmt12     ,& ! redundancy of aircraft report
    nfmt13     ,& ! redundancy of wind
    nfmt14     ,& ! redundancy of temperature
    nfmt15     ,& ! redundancy of humidity
    nfmt16     ,& ! redundancy of pressure / height
    nfmt17     ,& ! thinning of aircraft reports
    nfmt18     ,& ! exaggerated flight colocation
    nfmt19     ,& ! flight track error
    nfmt20     ,& ! message only: fog and precipitation
    nfmt21     ,& ! message only: fog and invisible sky
    nfmt22     ,& ! message only: fog and no cloud base
    nfmt23     ,& ! message only: cloud and no cloud base or fog
    nfmt24     ,& ! report (partly) blacklisted
    nfmt25     ,& ! report not on whitelist
    nfmt26        ! suspicious aircraft identity

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.     Header formats of ODR reports
!       ------------------------------------

!       1.1.1  Header formats of ODR reports: 'omlhed' and 'osghed'
!              ----------------------------------------------------
!   mxrhed     ,& ! header length of multi-level reports
!   mxshed     ,& ! header length of single-level reports
!   mxghed     ,& ! header length of GPS reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! (exact) time of observation in forecast hours
    nhsurf     ,& ! height of model grid pt. to which obs. is assigned
!   nhzio      ,& ! longitude of obs. station (or lowest datum) in grid pt. unit
!   nhzjo      ,& ! latitude  of obs. station in grid pt. units
!   nhsynt     ,& ! nominal (synoptic) time of observation in forecast hours

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mosghd'
!              ----------------------------------------------------
!   mxrhdf     ,& ! header length of multi-level reports
!   mxshdf     ,& ! header length of single-level reports
!   mxghdf     ,& ! header length of GPS reports
!   nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
!   nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot     ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
!   nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
!   nhqcfw     ,& ! threshold quality control flags for pressure, status of
!                 ! verification output
!   nhflag     ,& ! report flags (obs type, surface, altitude, station ID)
!   nhcorr     ,& ! update sequence number (station correction indicator)
!   nhcat      ,& ! data     category (from BUFR Section 1)
!   nhcats     ,& ! data sub-category (from BUFR Section 1)
!   nhkz       ,& ! DWD internal classification number (observation type)
!   nhcent     ,& ! originating centre
!   nhstid     ,& ! station identity number
!   nhdate     ,& ! absolute exact observation date [yyyymmdd]
!   nhhrmn     ,& ! absolute exact observation time [hhmm]
!   nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    nhnlev     ,& ! number of obs. levels (for multi-level reports)
!   nhvqcf     ,& ! for satellite retrieval: threshold quality control flags
!   nhaexi     ,& ! flag for exist. of wind or temperature in multi-level report
!   nhuexi     ,& ! flag for existence of wind data        in multi-level report
!   nhtexi     ,& ! flag for existence of temperature data in multi-level report
!   nhqexi     ,& ! flag for existence of humidity data    in multi-level report
!   nhrtyp     ,& ! radiosonde type    (NRARA, see WMO common code Table C2)
!   nhtrac     ,& ! tracking technique (NSASA, see WMO common code Table C7)
!   nhrad      ,& ! solar and IR radiation correction (NSR, BUFR Table 002013)
!   nhna4      ,& ! instrument type                   (NA4, BUFR Table 002003)
!   nhwce      ,& ! wind comput. enhancement (w-prof, MWCE, BUFR Table 025021)
!   nhdt       ,& ! time period of measurement (e.g. w-prof)               [s]
!   nhstyp     ,& ! surface obs: station type (buoy: MQOBL, BUFR Table 002149,
                  !                            else: NIX  , BUFR Table 002001)

!       1.1.3  Header formats of ODR reports: 'yomlhd' and 'yosghd'
!              ----------------------------------------------------

    ilstid     ,& ! character length of the station identity
    ilstidp       ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------
!       1.3.1  Body format of ODR of multi-level reports: 'omlbdy'
!              ---------------------------------------------------
!   mxrbdy     ,& ! body length of multi-level reports
    nbtu       ,& ! u wind component [m/s]
    nbtv       ,& ! v wind component [m/s]
    nbtt       ,& ! temperature [K]
    nbtrh      ,& ! relative humidity [/]
    nbtp       ,& ! pressure [Pa]
    nbtz       ,& ! height [m]
    nbtuer     ,& ! error of observed wind component
    nbtter     ,& ! error of observed temperature
    nbtqer     ,& ! error of observed rel. humidity
    nbtzer     ,& ! error of observed height
                  !  Note: nbt?er are set to the negative rms errors, if the
                  !  observations have not passed the threshold quality control

!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------
!   mxrbdf     ,& ! body length of multi-level reports
!   nbtflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
!   nbtqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
!   nbtlsg     ,& ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid     ,& ! level identity          (bit pattern, see below: 'nb?lid')

!       1.3.3  Body format of ODR of surface reports: 'osgbdy'
!              -----------------------------------------------
!   mxsbdy     ,& ! body length of single-level reports
    nbsu       ,& ! u wind component                                   [m/s]
    nbsv       ,& ! v wind component                                   [m/s]
    nbst       ,& ! temperature                                        [K]
    nbsrh      ,& ! relative humidity                                  [/]
    nbsp       ,& ! pressure                                           [Pa]
    nbsz       ,& ! height                                             [m]
    nbsuer     ,& ! error of observed wind component
    nbster     ,& ! error of observed temperature
    nbsqer     ,& ! error of observed relative humidity
    nbszer     ,& ! error of observed height
    nbsrr1,     & ! precipitation amount over 1  hour                  [mm]
    nbsrr3,     & ! precipitation amount over 3  hours                 [mm]
    nbsfg1,     & ! max. wind speed of gusts over 1 hour               [m/s]
    nbsfg3,     & ! max. wind speed of gusts over 3 hours              [m/s]
    nhkz       ,& ! DWD internal classification number (observation type)

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------
!   mxsbdf     ,& ! body length of single-level reports
!   nbsflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbserr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
!   nbsqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
!   nbslid     ,& ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                  ! else : level identity   (bit pattern, see below: 'nb?lid')

!       1.3.5  Body format of ODR of GPS reports: 'ogpbdy'
!              -------------------------------------------
!   mxgbdy     ,& ! body length of GPS reports
    nbgtze     ,& ! error in total zenith delay [mm]
    nbgzpd     ,& ! zenith path delay (total zenith delay)             [mm]
    nbgzwd     ,& ! zenith wet delay [mm]
    nbgiwv     ,& ! integrated water vapour [mm]
    nbgp       ,& ! pressure [Pa]
    nbgt       ,& ! temperature [K]
    nbgrh      ,& ! relative humidity [/]
    nbgbia     ,& ! bias correction to integrated water vapour [mm]
    nbgiwa     ,& ! adjusted (bias corrected) integrated water vapour [mm]

!       1.3.6  Body format of ODR of GPS report flags: 'mogpbd'
!              -------------------------------------------------
!   mxgbdf     ,& ! body length of GPS reports
!   nbgflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr        ! status flag word        (bit pattern, see below: 'nb?err')
!   nbgqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
!   nbglid        ! level identity          (bit pattern, see below: 'nb?lid')

USE data_obs_record, ONLY :   &

!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru       ,& ! bit pos. for status 'active' for horiz. wind      "
    nvrt       ,& ! bit pos. for status 'active' for temperature      "
    nvrq       ,& ! bit pos. for status 'active' for humidity         "
    nvrz       ,& ! bit pos. for status 'active' for pressure/height  "
    nvriwv     ,& ! bit pos. for status 'active' for IWV              "

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
!   ntotml     ,& ! tot. number of stored multi-level reports
!   ntotsg     ,& ! tot. number of stored single-level reports
!   ntotgp     ,& ! tot. number of stored GPS reports

!       2.     Observation data records (ODR)
!       -------------------------------------
    omlbdy     ,& ! body   of multi-level ODR
    omlhed     ,& ! header of multi-level ODR
    osgbdy     ,& ! body   of single-level ODR
    osghed     ,& ! header of single-level ODR
    ogpbdy     ,& ! body   of GPS ODR
    ogphed     ,& ! header of GPS ODR
    momlbd     ,& ! body   of multi-level ODR
    momlhd     ,& ! header of multi-level ODR
    mosgbd     ,& ! body   of single-level ODR
    mosghd     ,& ! header of single-level ODR
    mogpbd     ,& ! body   of GPS ODR
    mogphd     ,& ! header of GPS ODR
    yomlhd     ,& ! header of multi-level ODR
    yosghd     ,& ! header of single-level ODR
    yogphd        ! header of GPS ODR

! end of data_obs_record

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    global_values   ,& ! computes global values by operating on local arrays
    gather_values      ! gathers a set of values from all nod. to 1 or all nodes

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_gather_buffer  ! gather buffer arrays (used for diagn. printing)

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
PUBLIC :: obs_cdf_print_statist,  obs_cdf_print_events,  obs_cdf_print_reject, &
          obs_cdf_print_odr,  obs_print_number_rep

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!===============================================================================
!+ Module procedure in "src_obs_cdfin_print" for printing diagnostic statistics
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_print_statist ( iopen , num_compute , my_cart_id            &
                                 , icomm_cart , imp_integers )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfin_print" prints a summary of
!   the statistics on the processed observational reports.
!
! Method:
!   The diagnostic counters of array 'cma' are evaluated to produce a summary
!   of the reports processed and rejected.
!   Formatted printing.
!   If the program is run in parallel mode, the diagnostic arrays are summed up
!   over all nodes.
!   02.01.09
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

! Subroutine arguments: None
! --------------------

  INTEGER        , INTENT (IN)    :: &
    iopen          ,& ! = 0 : do not open or  close file with unit 'nustat'
                      ! = 1 : only   open           file with unit 'nustat'
                      ! =-1 : only            close file with unit 'nustat'
                      ! else: do     open and close file with unit 'nustat'
    num_compute    ,& ! number of compute PEs
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_integers      ! INTEGER   type used for MPI

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    inoct    (4+4*n_cma) ,& ! 1-d array of counters
    icnt_all (4)         ,& ! counters for all reports
    icnt_ob              ,& ! counter for processed reports of current obs type
    ima, ima2            ,& ! loop indices over array 'cma'
    izerror   , istat       ! error status

  CHARACTER (LEN= 5)       ::  yform10, yform11                      ! formats
  CHARACTER (LEN=42)       ::  yform72                               ! format
  CHARACTER (LEN=40)       ::  yform71                               ! format
  CHARACTER (LEN=43)       ::  yform70                               ! format
  CHARACTER (LEN=52)       ::  yform73                               ! format
  CHARACTER (LEN=49)       ::  yform18                               ! format
  CHARACTER (LEN=83)       ::  yform17                               ! format

  CHARACTER (LEN=80)       ::  yzerrmsg
  CHARACTER (LEN=75)       ::  yerrmsg    ! error message
  CHARACTER (LEN=21)       ::  yroutine   ! name of this subroutine

! Local arrays:
! ------------
!
!------------ End of header ----------------------------------------------------

  izerror  = 0
  yzerrmsg = '   '
  yroutine = 'obs_cdf_print_statist'

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_print_statist
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Summing of statistics counters over all nodes
!-------------------------------------------------------------------------------

  IF (num_compute > 1) THEN
    DO ima = 0 , n_cma
      inoct (4*ima+1)  =  cma(ima)%cnt_pr
      inoct (4*ima+2)  =  cma(ima)%cnt_ac
      inoct (4*ima+3)  =  cma(ima)%cnt_ps
      inoct (4*ima+4)  =  cma(ima)%cnt_rj
    ENDDO

    CALL global_values ( inoct, 4*n_cma+4, 'SUM', imp_integers                 &
                       , icomm_cart, 0, yzerrmsg, izerror )
!   ==================

    IF (my_cart_id == 0) THEN
      DO ima = 0 , n_cma
        cma(ima)%cnt_pr  =  inoct(4*ima+1)
        cma(ima)%cnt_ac  =  inoct(4*ima+2)
        cma(ima)%cnt_ps  =  inoct(4*ima+3)
        cma(ima)%cnt_rj  =  inoct(4*ima+4)
      ENDDO
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: For each CMA observation type, sum up over all CMA code types
!             which relate to that observation type
!             (Up to now, only CMA-code-dependent counters have been filled)
!-------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN

! sum up over all code types
    icnt_all (:) = 0
    DO ima = 0 , n_cma
      icnt_all (1)  =  icnt_all(1) + cma(ima)%cnt_pr
      icnt_all (2)  =  icnt_all(2) + cma(ima)%cnt_ac
      icnt_all (3)  =  icnt_all(3) + cma(ima)%cnt_ps
      icnt_all (4)  =  icnt_all(4) + cma(ima)%cnt_rj
    ENDDO

! sum up over all code types of each observation type separately
    DO ima = 1 , n_cma
      IF (cma(ima)%cdtyp == 0) THEN
! this element sums up a whole observation type
        DO ima2 = 1 , n_cma
          IF (cma(ima2)%obtyp == cma(ima)%obtyp) THEN
            cma(ima)%cnt_pr  =  cma(ima)%cnt_pr  +  cma(ima2)%cnt_pr
            cma(ima)%cnt_ac  =  cma(ima)%cnt_ac  +  cma(ima2)%cnt_ac
            cma(ima)%cnt_ps  =  cma(ima)%cnt_ps  +  cma(ima2)%cnt_ps
            cma(ima)%cnt_rj  =  cma(ima)%cnt_rj  +  cma(ima2)%cnt_rj
          ENDIF
        ENDDO
      ENDIF
    ENDDO

!-------------------------------------------------------------------------------
!  Section 3: Do the writing
!-------------------------------------------------------------------------------

! open file 'nustat', if required
! -------------------------------

    IF ((iopen /= 0) .AND. (iopen /= -1)) THEN
      OPEN (nustat,FILE=yustats,FORM='FORMATTED',STATUS='UNKNOWN'              &
                               ,POSITION='APPEND',IOSTAT=istat)
      IF (istat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yustats FAILED'
        CALL model_abort (my_cart_id, 7005, yerrmsg, yroutine)
      ENDIF
    ENDIF

! set formats
! -----------

    yform10 = '("0")'
    yform11 = '("1")'
    yform17 (1 :54) = '("0 ---DISTRIBUTION OF PROCESSED/ACTIVE/PASSIVE/REJECT'
    yform17 (55:83) = 'ED REPORTS FOR ASSIMILATION")'
    yform18 = '("0",42X," processed   active  passive rejected")'

    yform71 = '("0 --- observation type",I2,4X,A16,4I9)'
    yform72 = '("      code type ",I3,6X,A20,I8,3(1X,I8))'
    yform70 = '("0 --- total number of reports  ",13X,4I9)'
    yform73 = '("      reports with unknown obs/code type ",3X,4I9)'

! write table header
! ------------------

    WRITE(nustat,yform11)
    WRITE(nustat,yform10)
    WRITE(nustat,yform17)
    WRITE(nustat,yform18)
    WRITE(nustat,yform10)

! write table body
! ----------------

    WRITE(nustat,yform70) icnt_all
    WRITE(nustat,yform73) cma(0)%cnt_pr, cma(0)%cnt_ac                         &
                        , cma(0)%cnt_ps, cma(0)%cnt_rj

    icnt_ob = 0
    DO ima = 1 , n_cma
! observation types
      IF (cma(ima)%cdtyp == 0) THEN
        WRITE( nustat, yform71 )  cma(ima)%obtyp , cma(ima)%name               &
                                , cma(ima)%cnt_pr, cma(ima)%cnt_ac             &
                                , cma(ima)%cnt_ps, cma(ima)%cnt_rj
        icnt_ob = cma(ima)%cnt_pr
! code types (print only if there are processed reports with current obs type)
!     ELSEIF (cma(ima)%cnt_pr > 0) THEN
      ELSEIF (icnt_ob > 0) THEN
        WRITE( nustat, yform72 )  cma(ima)%cdtyp , cma(ima)%name               &
                                , cma(ima)%cnt_pr, cma(ima)%cnt_ac             &
                                , cma(ima)%cnt_ps, cma(ima)%cnt_rj
      ENDIF
    ENDDO

! write supplementary comments
! ----------------------------

    WRITE(nustat,'("0")' )
    WRITE(nustat,'("  --- Notes on the table above:")' )
    WRITE(nustat,'(6X,''"Rejected"/"passive" means that the whole ''           &
                    &,''report is rejected /set passive.'')' )
    WRITE(nustat,'(6X,''Partly rejected and partly passive reports are ''      &
                    &,''labeled "active".'')' )
    WRITE(nustat,'(6X,''A report can be labeled "active" even if part ''       &
                    &,''of its data is black listed.'')' )
    WRITE(nustat,'("0")' )
    WRITE(nustat,'(6X,"Reports may only be rejected or set passive for "       &
                    &,"reasons given by report")' )
    WRITE(nustat,'(6X,"  events 1 - 13 , except for events 3 and 5 "           &
                    &,"(on station altitude) for")' )
    WRITE(nustat,'(6X,"  TEMPs and PILOTs. Hence, the number of these events " &
                    &,"must equal the")' )
    WRITE(nustat,'(6X,"  number of rejected and passive reports.")' )
    WRITE(nustat,'(6X,"In the verification mode (i.e if data are written to "  &
                    &,"the VOF), reports are")' )
    WRITE(nustat,'(6X,"  rejected in case of report events 1 - 3, 8 - 13, "    &
                    &,"and event 4 if the re-")' )
    WRITE(nustat,'(6X,"  port is outside the model domain. Otherwise the "     &
                    &,"reports are set passive,")' )
    WRITE(nustat,'(6X,"  except that already passive reports are also rejected"&
                    &," if they do not")' )
    WRITE(nustat,'(6X,"  contain any data or if they are redundant and a "     &
                    &,"subset of an active")' )
    WRITE(nustat,'(6X,"  report.")' )
    WRITE(nustat,'(6X,"  Without verification, reports are always rejected "   &
                    &,"for events 1 to 13.")' )

! close file 'nustat', if required
! --------------------------------

    IF ((iopen /= 0) .AND. (iopen /= 1)) THEN
      CLOSE ( nustat )
    ENDIF

  ENDIF ! (my_cart_id == 0)

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_print_statist
!-------------------------------------------------------------------------------

END SUBROUTINE obs_cdf_print_statist


!===============================================================================
!+ Module procedure in "src_obs_cdfin_print" for printing diagnostic statistics
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_print_events ( iopen, iactio, mxeve , nevent , num_compute  &
                                , my_cart_id, icomm_cart, imp_integers )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfin_print" prints a summary of
!   report or data events. The printed report events are ordered according
!   to the order of the checks.
!
! Method:
!   The events table is inspected to produce a summary of the events.
!   Formatted printing.
!   If the program is run in parallel mode, the diagnostic arrays summed up
!   over all nodes.
!   05.01.09
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

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)    :: &
    iactio         ,& ! = 0 : print report events (from 1 to mxreve=18)
                      ! = 1 : print data events from 1 to 18
                      ! = 2 : print data events from 19 to mxdeve=37
    iopen          ,& ! = 0 : do not open or  close file with unit 'nustat'
                      ! = 1 : only   open           file with unit 'nustat'
                      ! =-1 : only            close file with unit 'nustat'
                      ! else: do     open and close file with unit 'nustat'
    mxeve          ,& ! length of first dimension of record 'nevent'
    num_compute    ,& ! number of compute PEs
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_integers      ! INTEGER   type used for MPI

  INTEGER        , INTENT (INOUT) :: &
    nevent (mxeve,n_cma)   ! record containing the event counters

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    nfrsev           ,& ! number of first event to be printed
    nlasev           ,& ! number of last  event to be printed
    nevc             ,& ! number of events to be printed
    jevt             ,& ! loop index over report events
    ima              ,& ! loop index over array 'cma'
    izerror , istat     ! error status indicator

  CHARACTER (LEN=14)       ::  yevents
  CHARACTER (LEN=62)       ::  y1900 , yform90 , yform91 , yform92   ! formats
  CHARACTER (LEN=13)       ::  yform12                               ! format
  CHARACTER (LEN= 5)       ::  yform10, yform11                      ! formats


! Local (automatic) arrays:
! -------------------------

  INTEGER        ::  &
    nevtno (mxeve)      ! event indices

  INTEGER        , ALLOCATABLE :: &
    inevent      (:)    ! 1-d array containing the event counters

  CHARACTER (LEN=80)       ::  yzerrmsg
  CHARACTER (LEN=75)       ::  yerrmsg    ! error message
  CHARACTER (LEN=20)       ::  yroutine   ! name of this subroutine
!
!------------ End of header ----------------------------------------------------

  izerror  = 0
  yzerrmsg = '   '
  yroutine = 'obs_cdf_print_events'

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_print_events
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Summing of event counters over all nodes
!-------------------------------------------------------------------------------

  IF (iactio == 0) THEN
    nfrsev  =  1
    nlasev  =  mxreve
  ELSEIF (iactio == 1) THEN
    nfrsev  =  1
    nlasev  =  18
  ELSEIF (iactio == 2) THEN
    nfrsev  =  19
    nlasev  =  mxdeve
  ENDIF

  IF (num_compute > 1) THEN
    ALLOCATE ( inevent ((nlasev-nfrsev+1)*n_cma) , STAT=istat )

    nevc  =  nlasev - nfrsev + 1
    DO ima = 1 , n_cma
      DO jevt = nfrsev , nlasev
        inevent ((ima-1)*nevc+jevt-nfrsev+1)  =  nevent(jevt,ima)
      ENDDO
    ENDDO

    CALL global_values ( inevent, n_cma*(nlasev-nfrsev+1), 'SUM', imp_integers &
                       , icomm_cart, 0, yzerrmsg, izerror )
!   ==================

    IF (my_cart_id == 0) THEN
      DO ima = 1 , n_cma
        DO jevt = nfrsev , nlasev
          nevent (jevt,ima)  =  inevent((ima-1)*nevc+jevt-nfrsev+1)
        ENDDO
      ENDDO
    ENDIF
    DEALLOCATE ( inevent , STAT=istat )
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Writing of events distribution
!-------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN

! open file 'nustat', if required
! -------------------------------

    IF ((iopen /= 0) .AND. (iopen /= -1)) THEN
      OPEN (nustat,FILE=yustats,FORM='FORMATTED',STATUS='UNKNOWN'              &
                               ,POSITION='APPEND',IOSTAT=istat)
      IF (istat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yustats FAILED'
        CALL model_abort (my_cart_id, 7005, yerrmsg, yroutine)
      ENDIF
    ENDIF

! set formats
! -----------

    yform10 = '("0")'
    yform11 = '("0")'
    yform12 = '("0 *",A20)  '
    yform90 = '(A17,I3,I3,I4,I5,I4,I3,I5,I3,I5,I3,I5,I4,I5,I3,I4,I5,I3)    '
    yform91 = '(A16,I4,I3,I3,I2,I3,I4,I5,I3,I5,I4,I3,I4,I3,I5,2I3,I5,I3)    '
    yform92 = '(A15,I4,I3,I2,I5,I3,I4,I4,I2,I4,I4,I4,I4,3I3,I5,I2,2I3,I5)   '

! write events' definition
! ------------------------

    WRITE(nustat,yform11)
    WRITE(nustat,yform10)

    IF (iactio == 0) THEN
      y1900 = yform90
      WRITE( nustat,'("1 *** REPORT EVENTS DEFINITIONS (THEIR ORDER MATCHES "  &
                    &,"THE ORDER OF THE CHECKS):")' )
      DO  jevt  = nfrsev , nlasev
        WRITE( nustat,'(8X,A)') crepev(jevt)
      ENDDO
    ELSE
      IF (iactio == 1) THEN
        y1900 = yform91
        WRITE( nustat,'("1 *** DATA EVENTS DEFINITIONS (LEVEL EVENTS APPLY "   &
                      &,"TO MULTI-LEVEL DATA ONLY,")' )
        WRITE( nustat,'(31X,"THE ORDER OF ALL EVENTS MATCHES THE ORDER OF")' )
        WRITE( nustat,'(31X,"THE CHECKS EXCEPT FOR EVENT 8")' )
      ELSE
        y1900 = yform92
        WRITE( nustat,'("1 *** DATA EVENTS DEFINITIONS (CONTINUED):")' )
      ENDIF
      DO  jevt  = nfrsev , nlasev
        WRITE( nustat,'(8X,A)') cdatev(jevt)
      ENDDO
    ENDIF

! write table header
! ------------------

    WRITE( nustat, yform11 )
    WRITE( nustat, yform10 )
    DO jevt = nfrsev , nlasev
      nevtno (jevt) = jevt
    ENDDO
    yevents = '  events      '
    WRITE( nustat, y1900   ) yevents, (nevtno(jevt), jevt=nfrsev,nlasev)
    WRITE( nustat, '("   ------ ")' )

! write table body
! ----------------

    DO ima = 1 , n_cma
! observation types
      IF (cma(ima)%cdtyp == 0) THEN
        WRITE( nustat, yform12 )  cma(ima)%name
      ELSE
        WRITE( nustat, y1900   )  cma(ima)%name                                &
                               , (nevent(jevt,ima), jevt=nfrsev,nlasev)
      ENDIF
    ENDDO

! close file 'nustat', if required
! --------------------------------

    IF ((iopen /= 0) .AND. (iopen /= 1)) THEN
      CLOSE ( nustat )
    ENDIF

  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_print_events
!-------------------------------------------------------------------------------

END SUBROUTINE obs_cdf_print_events


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_cdfin_print" for writing control output to files
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_print_reject ( lwonl , num_compute , my_cart_id             &
                                , icomm_cart , imp_integers )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfin_print" writes out the
!   information collected in buffer arrays mainly on rejection of reports
!   to file units 'nurej' or 'nuodr' for diagnostic purposes.
!
! Method:
!   If the program is run in parallel mode, the data to be printed are first
!   gathered from all PEs at the root node. The communication is designed to
!   keep the memory demands moderate.
!   09.02.98 / 20.07.98
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

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)  :: &
    num_compute    ,& ! number of compute PEs
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_integers      ! INTEGER   type used for MPI

  LOGICAL                  , INTENT (IN)  :: &
    lwonl             ! .true. for the node (sub-domain) at which file with
                      !        the unit number 'nupr' is open

! Local parameters:
! ----------------


! Local scalars:
! -------------

  INTEGER        ::  &
    istat            ,& ! error status variables
    ireceiver        ,& ! rank of the receiving PE
    onl_cart_id      ,& ! 'my_cart_id' of node where (lwonl=.true)
    ipos   , ind     ,& ! address  within a buffer
    irlen            ,& ! record length
    ilentot          ,& ! total number of elements gathered in 'allout'
    ifmt             ,& ! format identification
    klev             ,& ! level index
    npww             ,& ! extracted present weather
    nvis             ,& ! extracted visibility
    nvisfl           ,& ! flag for visibility
    npwwfl           ,& ! flag for present weather
    npsdt            ,& ! pressure tendency (Pa/3h)
    itot   , jtot    ,& ! asigned g.p.
    iol    , jol     ,& ! (local) coord. of grid pt. assigned to obs
    icdtyp           ,& ! code type
    icdtypr, icdtypa ,& !
    nreplpr          ,& ! number of replaced data
    nactpr           ,& ! number of active report
    nflaga , nflagr  ,& ! variable flags
    istalt           ,& ! station altitude
    imllev , imrlev  ,& ! number of levels of a multi-level observation
    ilid             ,& ! flag word for level information
    icl                 ! loop index over characters in a string

  REAL (KIND=wp) ::  &
    rscale           ,& ! scaling factor
    riscal           ,& ! inverse of rscale
    zwrbuf(10)       ,& ! buffer array for writing
    zrr              ,& ! precipitation amount
    rtr              ,& ! precipitation time range
    thresh           ,& ! threshold for lapse rate / wind shear
    difval           ,& ! threshold for lapse rate / wind shear
    zfidef , zfideb  ,& ! foreward and backward flight track confidences
    zpptop           ,& ! top of affected pressure range
    zoblon , zoblat  ,& ! observation longitude and latitude
    zobtim           ,& ! observation time in forecast hour units
    zio    , zjo     ,& ! (local) coord. of grid pt. assigned to obs
    zuu    , zvv     ,& ! observed wind components
    ztt    , zrh     ,& ! observed temperature and relative humidity
    zpp              ,& ! observed pressure and height
!   zuer   , zter    ,& ! observation errors for wind and temperature
!   zqer   , zzer    ,& ! observation errors for rel.hum. and height
    zppac  , zpprj   ,& ! observed active and rejected pressure
    zuuac  , zuurj   ,& ! observed active and rejected wind
    zttac  , zttrj   ,& ! observed active and rejected temperature
    zqqac  , zqqrj      ! observed active and rejected humidity

  LOGICAL        ::  &
    lstcor           ,& ! .TRUE., if second report is a station correction
    lwr_now          ,& ! .TRUE., if (more than zero) messages are written now
    lprodr              ! .TRUE.

  CHARACTER (LEN=ilstidp)  :: &
    ystid, ystidr, ystida  ! station identity

  CHARACTER (LEN=20)       :: &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=50)       :: &
    yerrmsg             ! error message

! Local (automatic) arrays:
! -------------------------

  INTEGER        , POINTER      :: &
    allout  (:)         ! buffer containing output of all nodes
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_print_reject
!-------------------------------------------------------------------------------

  yroutine    = 'obs_cdf_print_reject'
  ireceiver   = 0
  rscale      = 100._wp
  riscal      = (1._wp+epsy)/100._wp
  lprodr      = .TRUE.
  onl_cart_id = 0
  IF (lwonl)  onl_cart_id = my_cart_id

!-------------------------------------------------------------------------------
!  Section 1: Gather data records stored in buffer "outbuf" and print out
!-------------------------------------------------------------------------------

  ! note that 'allout' is allocated in obs_gather_buffer

  CALL obs_gather_buffer ( -1, nmxoln, outbuf, ilentot, allout, 0              &
                         , num_compute, my_cart_id, icomm_cart, imp_integers )
! ======================

  lwr_now = (ilentot >= 1)

  IF (my_cart_id == 0) THEN

! Open files if required
! ----------------------
    IF (lwr_now) THEN
      istat = 0
      IF  (.NOT. lopen_rej)                                                    &
        OPEN (nurej ,FILE=yurejct ,FORM='FORMATTED',STATUS='UNKNOWN'           &
                                  ,POSITION='APPEND',IOSTAT=istat)
      IF ((.NOT. lopen_odr) .AND. (istat == 0))                                &
        OPEN (nuodr ,FILE=yuobsdr ,FORM='FORMATTED',STATUS='UNKNOWN'           &
                                  ,POSITION='APPEND',IOSTAT=istat)
      IF (istat /= 0) yerrmsg = 'OPENING OF FILE yurejct or yuobsdr FAILED'
      IF (istat /= 0) CALL model_abort (my_cart_id, 7004, yerrmsg, yroutine)
      lopen_rej = .TRUE.
      lopen_odr = .TRUE.
    ENDIF

! Print out messages on rejected reports and data
! -----------------------------------------------

    ipos = 0

    Get_next_record:   DO
!   =====================

!   get record length
      irlen = allout(ipos+1)
      IF (irlen == 0)                                       EXIT Get_next_record

      ifmt  = allout(ipos+2)

      IF (    ifmt == nfmt1) THEN
!     ===========================
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+2+icl) )
        ENDDO
        WRITE(nurej,'(" SINGLE LEV REP ",a," :  NO PRESSURE ")') ystid

      ELSEIF (ifmt == nfmt2) THEN
!     ===========================
        zrr = allout(ipos+3)*riscal
        rtr = REAL (allout(ipos+4), wp)
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+4+icl) )
        ENDDO
        WRITE(nurej,'(" PRECIPITATION AMOUNT EXCEEDS LIMIT -  ",               &
            &"  DATUM REJECTED.  STID=",a," RR=",f5.1," TR=",f3.0)')           &
              ystid,zrr,rtr

      ELSEIF (ifmt == nfmt3) THEN
!     ===========================
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+2+icl) )
        ENDDO
        WRITE (nurej,      '(" SINGLE LEV REP ",a," : ",                       &
            &"NO ACCEPTED DATA IN REPORT")') ystid

      ELSEIF (ifmt == nfmt4) THEN
!     ===========================
        klev = allout(ipos+3)
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+3+icl) )
        ENDDO
        WRITE (nurej,       '(" MULTI LEV REP  ",a," : ",i5,"th LEVEL,",       &
            &" BUT ODR SIZE IS ",i5)') ystid, klev, maxmlv

      ELSEIF (ifmt == nfmt5) THEN
!     ===========================
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+2+icl) )
        ENDDO
        WRITE (nurej,       '(" MULTI LEV REP  ",a," : ",                      &
            &"SEVERAL SURFACE LEVELS")') ystid

      ELSEIF (ifmt == nfmt6) THEN
!     ===========================
        nflaga = allout(ipos+3)
        npsdt  = allout(ipos+4)
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+4+icl) )
        ENDDO
        WRITE (nurej,       '(" SINGLE LEV REP  ",a," : ",                     &
            &"PRESSURE TENDENCY:",I6," , FLAG:",I2)') ystid, npsdt, nflaga

      ELSEIF (ifmt == nfmt7) THEN
!     ===========================
        thresh = allout(ipos+3)*riscal
        difval = allout(ipos+4)*riscal
        zpp    = allout(ipos+5)*riscal
        zpptop = allout(ipos+6)*riscal
        zobtim = allout(ipos+7)*riscal
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+7+icl) )
        ENDDO
        WRITE (nurej,       '(" LAPSE RATE  ",A ,",",F5.1,": ",                &
            &"THRESHOLD , VALUE:",F5.1,F6.2," , P:",F6.0," -",F6.0)')          &
             ystid, zobtim, thresh, difval, zpp, zpptop

      ELSEIF (ifmt == nfmt8) THEN
!     ===========================
        thresh = allout(ipos+3)*riscal
        difval = allout(ipos+4)*riscal
        zpp    = allout(ipos+5)*riscal
        zpptop = allout(ipos+6)*riscal
        zobtim = allout(ipos+7)*riscal
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+7+icl) )
        ENDDO
        WRITE (nurej,       '(" WIND SPEED SHEAR  ",A ,",",F5.1,": ",          &
            &"THRESH, VALUE:",F5.1,F6.1," , P:",F6.0," -",F6.0)')              &
             ystid, zobtim, thresh, difval, zpp, zpptop

      ELSEIF (ifmt == nfmt9) THEN
!     ===========================
        thresh = allout(ipos+3)*riscal
        difval = allout(ipos+4)*riscal
        zpp    = allout(ipos+5)*riscal
        zpptop = allout(ipos+6)*riscal
        zobtim = allout(ipos+7)*riscal
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+7+icl) )
        ENDDO
        WRITE (nurej,       '(" DIRECTIONAL SHEAR  ",A ,",",F5.1,": ",         &
            &"THRESH, VALUE:",F5.0,F5.0," , P:",F6.0," -",F6.0)')              &
             ystid, zobtim, thresh, difval, zpp, zpptop

      ELSEIF (ifmt == nfmt17) THEN
!     ============================
        zpp    = allout(ipos+3)*riscal
        zobtim = allout(ipos+4)*riscal
        thresh = allout(ipos+5)*riscal
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+5+icl) )
        ENDDO
        WRITE (nurej,       '(" FLIGHT TRACK THINNING  ",A ,",",F6.2,",",F6.0  &
                            &,": (TOO) CLOSE TO OBS. TIME",F6.2)' )            &
               ystid, zobtim, zpp, thresh

      ELSEIF (ifmt == nfmt18) THEN
!     ============================
        zobtim = allout(ipos+5)*riscal
        thresh = allout(ipos+6)*riscal
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+6+icl) )
        ENDDO
        IF (allout(ipos+4) == 1) THEN
          WRITE (nurej,       '(" EXAGGERATED HORIZONTAL COLOCATION ",A        &
                              &,":",I4," REPORTS FROM",F6.2," TO",F6.2)' )     &
                 ystid, allout(ipos+3), zobtim, thresh
        ELSE
          WRITE (nurej,       '(" EXAGGERATED VERTICAL COLOCATION ",A          &
                              &,":",I4," REPORTS FROM",F6.2," TO",F6.2)' )     &
                 ystid, allout(ipos+3), zobtim, thresh
        ENDIF

      ELSEIF (ifmt == nfmt19) THEN
!     ============================
        zpp    = allout(ipos+3)*riscal
        zobtim = allout(ipos+4)*riscal
        zfidef = allout(ipos+6)*riscal
        zfideb = allout(ipos+7)*riscal
        nflaga = NINT( zfidef )
        nflagr = NINT( zfideb )
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+7+icl) )
        ENDDO
        IF (allout(ipos+5) == 1) THEN
          WRITE (nurej,       '(" FLIGHT TRACK CHECK  ",A ,",",F5.2,1X,F5.0    &
                              &,": LON SIGN, FOREWARD CONFIDENCE:",I3,I4)' )   &
                 ystid, zobtim, zpp, nflaga, nflagr
        ELSEIF (allout(ipos+5) == 2) THEN
          WRITE (nurej,       '(" FLIGHT TRACK CHECK  ",A ,",",F5.2,1X,F5.0    &
                              &,": LON SIGN, BACKWARD CONFIDENCE:",I3,I4)' )   &
                 ystid, zobtim, zpp, nflaga, nflagr
        ELSEIF (allout(ipos+5) == 3) THEN
          WRITE (nurej,       '(" FLIGHT TRACK CHECK  ",A ,",",F5.2,1X,F5.0    &
                              &,": HORIZONTAL CONFIDENCES:",2F6.1)' )          &
                 ystid, zobtim, zpp, zfidef, zfideb
        ELSEIF (allout(ipos+5) == 4) THEN
          WRITE (nurej,       '(" FLIGHT TRACK CHECK  ",A ,",",F5.2,1X,F5.0    &
                              &,": VERTICAL CONFIDENCES:",2F6.1)' )            &
                 ystid, zobtim, zpp, zfidef, zfideb
        ENDIF

      ELSEIF (ifmt == nfmt20) THEN
!     ============================
        npww   = allout(ipos+3)
        nvis   = allout(ipos+4)
        npwwfl = allout(ipos+5)
        nvisfl = allout(ipos+6)
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+6+icl) )
        ENDDO
        WRITE( nurej ,'(A,'': Fog with precip.: weather:'',I3,'' ,vis.:'',I6   &
                      &,'' ,flags:'',2I3)' )  ystid, npww, nvis, npwwfl, nvisfl

      ELSEIF (ifmt == nfmt24) THEN
!     ============================
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+11+icl) )
        ENDDO
        WRITE (nurej,       '(" STA ",A ," OBTYP",I2," BLACKLISTED (Z,V,T,Q):" &
                            &,4(I5,I4))' )                                     &
               ystid, (allout(ipos+icl), icl=3,11)

      ELSEIF (ifmt == nfmt25) THEN
!     ============================
        zobtim = allout (ipos+ 4)*riscal
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+4+icl) )
        ENDDO
        WRITE (nurej,       '(" STATION ",A ," CODE TYPE",I4                   &
                            &," NOT ON WHITELIST, ",F6.1," HRS")' )            &
               ystid, allout(ipos+ 3), zobtim

      ELSEIF (ifmt == nfmt26) THEN
!     ============================
        zobtim = allout (ipos+ 3)*riscal
        zpp    = allout (ipos+ 4)*riscal
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+4+icl) )
        ENDDO
        WRITE (nurej,       '(" STATION ",A ," : SUSPICIOUS AIRCRAFT ID",F7.1  &
                            &," HPA , ",F6.1," HRS")' )                        &
               ystid, zpp, zobtim


      ELSEIF (ifmt == nfmt10 .AND. lprodr) THEN
!     =========================================
        icdtypr= allout(ipos+3)
        zuu = rmdi
        zvv = rmdi
        ztt = rmdi
        zrh = rmdi
        zpp = rmdi
        IF (allout (ipos+ 4) /= imdi)  zuu = allout (ipos+ 4)*riscal
        IF (allout (ipos+ 5) /= imdi)  zvv = allout (ipos+ 5)*riscal
        IF (allout (ipos+ 6) /= imdi)  ztt = allout (ipos+ 6)*riscal
        IF (allout (ipos+ 7) /= imdi)  zrh = allout (ipos+ 7)*riscal
        IF (allout (ipos+ 8) /= imdi)  zpp = allout (ipos+ 8)*riscal
        icdtypa= allout (ipos+ 9)
        zoblon = allout (ipos+10)*riscal
        zoblat = allout (ipos+11)*riscal
        nactpr = allout (ipos+12)
        lstcor =(allout (ipos+13) == 1)
        nreplpr= allout (ipos+14)
        zobtim = allout (ipos+15)*riscal
        DO icl = 1 , istrej
          ystidr (icl:icl) = CHAR( allout(ipos+15       +icl) )
          ystida (icl:icl) = CHAR( allout(ipos+15+istrej+icl) )
        ENDDO
        WRITE( nuodr,'(" --> REDUNDANT: ",a  , i4                              &
                     &,"  / ACTIVE: " , 3f7.1, f5.2, f8.0, 1x, a , i4          &
                     &, 2f5.1, i2, 1x, l1, i4, f5.1)' ) ystidr, icdtypr        &
              ,zuu, zvv, ztt, zrh, zpp, ystida, icdtypa                        &
              ,zoblon, zoblat, nactpr , lstcor , nreplpr, zobtim

      ELSEIF (ifmt == nfmt11 .AND. lprodr) THEN
!     =========================================
        icdtypr= allout(ipos+3)
        itot   = allout(ipos+4)
        jtot   = allout(ipos+5)
        zoblon = allout(ipos+6)*riscal
        zoblat = allout(ipos+7)*riscal
        istalt = allout(ipos+8)
        zobtim = allout(ipos+9)*riscal
        icdtypa= allout(ipos+10)
        imllev = allout(ipos+11)
        imrlev = allout(ipos+12)
        DO icl = 1 , istrej
          ystidr (icl:icl) = CHAR( allout(ipos+12       +icl) )
          ystida (icl:icl) = CHAR( allout(ipos+12+istrej+icl) )
        ENDDO
        WRITE( nuodr,'(''--> REDUNDANT MULTI-LEVEL REPORT: '',a ,i4)' )        &
               ystidr, icdtypr
        WRITE( nuodr,'("     u      v     t    rh      p      z  ",            &
                      &"v-err t-er rh-er z-er lev")' )
        ind  = ipos + 12 + 2*istrej

        DO klev = 1, imrlev + imllev
          DO icl = 1, 10
            zwrbuf (icl) = rmdi
            IF (allout(ind+icl) /= imdi)  zwrbuf (icl) = allout(ind+icl) *riscal
          ENDDO
          ilid = allout (ind +11)
          ind  = ind  + 11
          WRITE( nuodr,'(3F7.1,F5.2, F8.0,F7.0, 2F5.1,F5.2,F5.1,I4)' )         &
                (zwrbuf(icl), icl=1,10), ilid
!                zuu, zvv, ztt, zrh, zpp, zzz, zuer, zter, zqer, zzer, ilid
          IF (klev == imrlev) THEN
            WRITE( nuodr,'(" --> ACTIVE: ",a ," G.P.=(",i7,",",i3,")  LON="    &
                         &,f6.2,"  LAT=",f6.2,"  ST.H=",i5,"  TIME(h)=",f5.1   &
                         &,"  CDTYPE=",i3," NUMLEV=",i4)')    ystida           &
                 , itot, jtot, zoblon, zoblat, istalt, zobtim, icdtypa, imllev
            WRITE( nuodr,'("     u      v     t    rh      p      z  ",        &
                          &"v-err t-er rh-er z-er lev")' )
          ENDIF
        ENDDO

      ELSEIF (ifmt == nfmt12 .AND. lprodr) THEN
!     =========================================
        zuu = rmdi
        zvv = rmdi
        ztt = rmdi
        zrh = rmdi
        zpp = rmdi
        IF (allout (ipos+ 3) /= imdi)  zuu = allout (ipos+ 3)*riscal
        IF (allout (ipos+ 4) /= imdi)  zvv = allout (ipos+ 4)*riscal
        IF (allout (ipos+ 5) /= imdi)  ztt = allout (ipos+ 5)*riscal
        IF (allout (ipos+ 6) /= imdi)  zrh = allout (ipos+ 6)*riscal
        IF (allout (ipos+ 7) /= imdi)  zpp = allout (ipos+ 7)*riscal
        icdtyp = allout (ipos+8)
        zio    = allout (ipos+9)*riscal
        zjo    = allout (ipos+10)*riscal
        iol    = allout (ipos+11)
        jol    = allout (ipos+12)
        zobtim = allout (ipos+13)*riscal
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+13+icl) )
        ENDDO
        WRITE( nuodr,'('' --> REDUNDANT AIREP: ACTIVE:'',3f7.1,f5.2,f8.0       &
                     &,1x, a , i4, 2f5.1, 2i4, f5.1)' )                        &
               zuu, zvv, ztt, zrh, zpp, ystid ,icdtyp, zio, zjo, iol,jol, zobtim

      ELSEIF (ifmt == nfmt13 .AND. lprodr) THEN
!     =========================================
        zppac = allout (ipos+ 3)*riscal
        zpprj = allout (ipos+ 4)*riscal
        zuuac = rmdi
        zuurj = rmdi
        IF (allout (ipos+ 5) /= imdi)  zuuac = allout (ipos+ 5)*riscal
        IF (allout (ipos+ 6) /= imdi)  zuurj = allout (ipos+ 6)*riscal
        nflaga = allout (ipos+ 7)
        nflagr = allout (ipos+ 8)
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+8+icl) )
        ENDDO
        WRITE( nuodr,'(1x,a,'' U: ACT/REJ: P/U/FLAG:'',2f8.0,2f6.1,2i2)')      &
               ystid ,zppac, zpprj, zuuac, zuurj, nflaga, nflagr

      ELSEIF (ifmt == nfmt14 .AND. lprodr) THEN
!     =========================================
        zppac = allout (ipos+ 3)*riscal
        zpprj = allout (ipos+ 4)*riscal
        zttac = rmdi
        zttrj = rmdi
        IF (allout (ipos+ 5) /= imdi)  zttac = allout (ipos+ 5)*riscal
        IF (allout (ipos+ 6) /= imdi)  zttrj = allout (ipos+ 6)*riscal
        nflaga = allout (ipos+ 7)
        nflagr = allout (ipos+ 8)
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+8+icl) )
        ENDDO
        WRITE( nuodr,'(1x,a,'' T: ACT/REJ: P/T/FLAG:'',2f8.0,2f6.1,2i2)')      &
               ystid ,zppac, zpprj, zttac, zttrj, nflaga, nflagr

      ELSEIF (ifmt == nfmt15 .AND. lprodr) THEN
!     =========================================
        zppac = allout (ipos+ 3)*riscal
        zpprj = allout (ipos+ 4)*riscal
        zqqac = rmdi
        zqqrj = rmdi
        IF (allout (ipos+ 5) /= imdi)  zqqac = allout (ipos+ 5)*riscal
        IF (allout (ipos+ 6) /= imdi)  zqqrj = allout (ipos+ 6)*riscal
        nflaga = allout (ipos+ 7)
        nflagr = allout (ipos+ 8)
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+8+icl) )
        ENDDO
        WRITE( nuodr,'(1x,a,'' Q: ACT/REJ: P/Q/FLAG:'',2f8.0,2f6.1,2i2)')      &
               ystid ,zppac, zpprj, zqqac, zqqrj, nflaga, nflagr

      ELSEIF (ifmt == nfmt16 .AND. lprodr) THEN
!     =========================================
        zppac = allout (ipos+ 3)*riscal
        zpprj = allout (ipos+ 4)*riscal
        zqqac = rmdi
        zqqrj = rmdi
        zttac = rmdi
        zttrj = rmdi
        IF (allout (ipos+ 5) /= imdi)  zqqac = allout (ipos+ 5)*riscal
        IF (allout (ipos+ 6) /= imdi)  zqqrj = allout (ipos+ 6)*riscal
        IF (allout (ipos+ 7) /= imdi)  zttrj = allout (ipos+ 7)*riscal
        IF (allout (ipos+ 8) /= imdi)  zttrj = allout (ipos+ 8)*riscal
        DO icl = 1 , istrej
          ystid (icl:icl) = CHAR( allout(ipos+8+icl) )
        ENDDO
        WRITE( nuodr,'(1x,a,'' Z: ACT/REJ: P/ZERR/Z:'',2f8.0,2f6.1,2f7.0)')    &
               ystid ,zppac, zpprj, zqqac, zqqrj, zttac, zttrj
      ENDIF

      ipos = ipos + irlen

! Close loop : fetch next record if necessary
! -------------------------------------------
    ENDDO Get_next_record

    IF (lopen_rej) THEN
      CLOSE (nurej)
      lopen_rej = .FALSE.
    ENDIF
  ENDIF  ! my_cart_id == 0

  DEALLOCATE ( allout    , STAT = istat )

! Print out messages on derived fog and low cloud
! -----------------------------------------------

! ! determine onl_cart_id
! IF (num_compute > 1) CALL global_values ( onl_cart_id, 1,'MAX', imp_integers &
!                                         , icomm_cart, -1, yerrmsg, implcode )
!                      ==================
! ! note that 'allout' is allocated in obs_gather_buffer
! CALL obs_gather_buffer ( -1, nmxoln, outbuf, ilentot, allout, onl_cart_id    &
!                        , num_compute, my_cart_id, icomm_cart, imp_integers )
! ======================

! IF (lwonl) THEN

!   ipos = 0

!   Get_next_message:   DO
!   ======================

!   get message length
!     irlen = allout(ipos+1)
!     IF (irlen == 0)                                      EXIT Get_next_message

!     ifmt  = allout(ipos+2)

!     IF (    ifmt == nfmt20) THEN
!     ============================
!       npww   = allout(ipos+3)
!       nvis   = allout(ipos+4)
!       npwwfl = allout(ipos+5)
!       nvisfl = allout(ipos+6)
!       DO icl = 1 , istrej
!         ystid (icl:icl) = CHAR( allout(ipos+6+icl) )
!       ENDDO
!       WRITE( nupr ,'(A,'': Fog with precip.: weather:'',I3,'' ,vis.:'',I6    &
!                    &,'' ,flags:'',2I3)' )  ystid, npww, nvis, npwwfl, nvisfl

!     ELSEIF ((ifmt >= nfmt21) .AND. (ifmt <= nfmt23)) THEN
!     =====================================================
!       nccl   = allout(ipos+3)
!       ncclfl = allout(ipos+4)
!       ncclqf = allout(ipos+5)
!       DO icl = 1 , istrej
!         ystid (icl:icl) = CHAR( allout(ipos+5+icl) )
!       ENDDO
!       IF     (ifmt == nfmt21) THEN
!         WRITE( nupr ,'(A,'': Fog, sky invisible:   low cloud'',I3            &
!                      &,'' ,flags:'',2I3)' )  ystid, nccl, ncclfl, ncclqf
!       ELSEIF (ifmt == nfmt22) THEN
!         WRITE( nupr ,'(A,'': Fog, no cloud base:   low cloud'',I3            &
!                      &,'' ,flags:'',2I3)' )  ystid, nccl, ncclfl, ncclqf
!       ELSEIF (ifmt == nfmt23) THEN
!         WRITE( nupr ,'(A,'': No cloud base or fog: low cloud'',I3            &
!                      &,'' ,flags:'',2I3)' )  ystid, nccl, ncclfl, ncclqf
!       ENDIF
!     ENDIF

!     ipos = ipos + irlen

! Close loop : fetch next message if necessary
! --------------------------------------------
!   ENDDO Get_next_message

! ENDIF  ! lwonl

! DEALLOCATE ( allout    , STAT = istat )

  DEALLOCATE ( outbuf    , STAT = istat )

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_print_reject
!-------------------------------------------------------------------------------

END SUBROUTINE obs_cdf_print_reject


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_cdfin_print" for writing control output to files
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_print_odr ( nodrtot , nodrold , num_compute , my_cart_id    &
                             , icomm_cart , imp_integers , mobname , chobtp )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfin_print" writes out a main part
!   of the ODR (Observation Data Record) to file unit 'nuodr' for diagnostic
!   purposes.
!
! Method:
!   If the program is run in parallel mode, the data to be printed are first
!   gathered from all PEs at the root node. The communication is designed to
!   keep the memory demands moderate.
!   09.02.98 / 20.07.98
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

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)  :: &
    nodrtot (3)    ,& ! total number of multi-level / single-level / GPS reports
    nodrold (3)    ,& ! number of multi-level / single-level / GPS reports
                      !   before having read new data at the current timestep
    num_compute    ,& ! number of compute PEs
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_integers   ,& ! INTEGER   type used for MPI
    mobname           ! = 0 : names of observation types derived from 'cma'
                      ! = 1 : names of observation types given via 'chobtp'

  CHARACTER (LEN=8)  , INTENT (IN)  , OPTIONAL  :: &
    chobtp (8)        ! names of observation types (if mobname=1, i.e. AOF read)

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    c100r  =  0.01_wp

! Local scalars:
! -------------

  INTEGER        ::  &
    nodrln           ,& ! dynamic length of ODR output buffer
    istat            ,& ! error status variables
    isgid            ,& ! flag for single-level observation
    imlid            ,& ! flag for multi-level observation
    igpid            ,& ! flag for GPS observation
    ireceiver        ,& ! rank of the receiving PE
    ipos   , indst   ,& ! address  within a buffer
    irlen            ,& ! record length
    ilentot          ,& ! total number of elements gathered in 'odrbufa'
    klev             ,& ! level index
    iob              ,& ! loop indices
    itot   , jtot    ,& ! asigned g.p.
    iobtyp           ,& ! observation type
    icdtyp           ,& ! code type
    istalt           ,& ! station altitude
    ihadif           ,& ! difference between station height and model height
    igporo           ,& ! model orography height
    imllev           ,& ! number of levels of a multi-level observation
    ilid             ,& ! flag word for level information
    istrlen          ,& ! length of string of station id
    icl                 ! loop index over characters in a string

  REAL (KIND=wp) ::  &
    rscale           ,& ! scaling factor
    riscal           ,& ! inverse of rscale
    zwrbuf(10)       ,& ! buffer array for writing
    zoblon , zoblat  ,& ! observation longitude and latitude
    zobtim              ! observation time in forecast hour units
!   zuer , zter      ,& ! observation errors for wind and temperatur
!   zqer , zzer      ,& ! observation errors for rel.hum. and height
!   ziwv , zzwd      ,& ! integrated water vapour / zenith wet delay
!   ztzd , ztze , zbia  ! total zenith delay / TZD error / bias correction

  LOGICAL        ::  &
    lwr_now             ! .TRUE., if (more than zero) messages are written now

  CHARACTER (LEN=ilstidp)  :: &
    ystid               ! station identity
  CHARACTER (LEN=1)        :: &
    ypass               ! indicator: single-level aircraft report set passive
  CHARACTER (LEN=5)        :: &
    yobtyp              ! name of observation type

  CHARACTER (LEN=20)       :: &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=30)       :: &
    yerrmsg             ! error message
  CHARACTER (LEN=70)       :: &
    yerrmsgl            ! error message

! Local (automatic) arrays:
! -------------------------

  INTEGER        , ALLOCATABLE  :: &
    odrbufp (:)         ! buffer containing ODR data of a single PE

  INTEGER        , POINTER      :: &
    odrbufa (:)         ! buffer containing ODR data of a all PEs
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_print_odr
!-------------------------------------------------------------------------------

  yroutine   = 'obs_cdf_print_odr'
  isgid      = 10
  igpid      = 8
  imlid      = 100
  ireceiver  = 0
  rscale     = 100._wp
  riscal     = (1._wp+epsy)/100._wp

!-------------------------------------------------------------------------------
!  Section 1:  Write local ODR data to buffer arrays
!-------------------------------------------------------------------------------

  istrlen = ilstid

! dynamic allocation of memory for the local integer buffer
! ---------------------------------------------------------
  nodrln = 0
  DO iob = nodrold(1)+1 , nodrtot(1)
    nodrln = nodrln + momlhd(iob,nhnlev)
  ENDDO
  nodrln  =   (nodrtot(1) - nodrold(1)) *(12+istrlen) + 11* nodrln             &
           +  (nodrtot(2) - nodrold(2)) *(12+istrlen  + 10)                    &
           +  (nodrtot(3) - nodrold(3)) *(12+istrlen  +  8)
  ALLOCATE ( odrbufp (nodrln+1)    , STAT=istat )
! odrbufp (:) = 0
  odrbufp (:) = imdi

! store ODR data to be printed in an intermediate integer buffer
! --------------------------------------------------------------
  ipos = 0

  DO iob = nodrold(2)+1 , nodrtot(2)
! ==================================

!   header of single-level reports (length = 12 + istrlen)
!   ------------------------------
    odrbufp(ipos+ 1) = 12 + istrlen + 10
    odrbufp(ipos+ 2) = isgid
    odrbufp(ipos+ 3) = mosghd(iob,nhitot)
    odrbufp(ipos+ 4) = mosghd(iob,nhjtot)
    odrbufp(ipos+ 5) = mosghd(iob,nhobtp)
    odrbufp(ipos+ 6) = mosghd(iob,nhcode)
    IF (mosghd(iob,nhkz) == 40005)  odrbufp(ipos+ 6) = 405
    IF (osghed(iob,nhalt) > rmdich) THEN
      odrbufp(ipos+7) = INT( osghed(iob,nhalt) )
      odrbufp(ipos+8) = INT( osghed(iob,nhalt) - osghed(iob,nhsurf) )
!   ELSE
!     odrbufp(ipos+7) = imdi
!     odrbufp(ipos+8) = imdi
    ENDIF
    odrbufp(ipos+ 9) = INT(osghed(iob,nhilon)*rscale)
    odrbufp(ipos+10) = INT(osghed(iob,nhjlat)*rscale)
    odrbufp(ipos+11) = INT(osghed(iob,nhtime)*rscale)
    odrbufp(ipos+12) = mosghd(iob,nhpass)
    ipos = ipos + 12

    DO icl = 1 , istrlen
      odrbufp(ipos+icl) = ICHAR( yosghd(iob) (icl:icl) )
    ENDDO
    ipos = ipos + istrlen

!   body of single-level reports (length = 10)
!   ----------------------------
!   DO icl = 1 , 10
!     odrbufp(ipos+icl) = imdi
!   ENDDO
    IF (osgbdy(iob,nbsu  ) > rmdich)                                           &
      odrbufp(ipos+ 1) = INT(osgbdy(iob,nbsu  )*rscale)
    IF (osgbdy(iob,nbsv  ) > rmdich)                                           &
      odrbufp(ipos+ 2) = INT(osgbdy(iob,nbsv  )*rscale)
    IF (osgbdy(iob,nbst  ) > rmdich)                                           &
      odrbufp(ipos+ 3) = INT(osgbdy(iob,nbst  )*rscale)
    IF (osgbdy(iob,nbsrh ) > rmdich)                                           &
      odrbufp(ipos+ 4) = INT(osgbdy(iob,nbsrh )*rscale)
    IF (osgbdy(iob,nbsp  ) > rmdich)                                           &
      odrbufp(ipos+ 5) = INT(osgbdy(iob,nbsp  )*rscale *c100r)
    IF (osgbdy(iob,nbsz  ) > rmdich)                                           &
      odrbufp(ipos+ 6) = INT(osgbdy(iob,nbsz  )*rscale)
!   IF (ibit1( mosgbd(iob,nbserr),nvru ) == 1)                                 &
!     odrbufp(ipos+ 7) = INT(osgbdy(iob,nbsuer)*rscale)
!   IF (ibit1( mosgbd(iob,nbserr),nvrt ) == 1)                                 &
!     odrbufp(ipos+ 8) = INT(osgbdy(iob,nbster)*rscale)
!   IF (ibit1( mosgbd(iob,nbserr),nvrq ) == 1)                                 &
!     odrbufp(ipos+ 9) = INT(osgbdy(iob,nbsqer)*rscale)
!   IF (ibit1( mosgbd(iob,nbserr),nvrz ) == 1)                                 &
!     odrbufp(ipos+10) = INT(osgbdy(iob,nbszer)*rscale)
    IF (osgbdy(iob,nbsrr1) > rmdich)                                           &
      odrbufp(ipos+ 7) = INT(osgbdy(iob,nbsrr1)*rscale)
    IF (osgbdy(iob,nbsrr3) > rmdich)                                           &
      odrbufp(ipos+ 8) = INT(osgbdy(iob,nbsrr3)*rscale)
    IF (osgbdy(iob,nbsfg1) > rmdich)                                           &
      odrbufp(ipos+ 9) = INT(osgbdy(iob,nbsfg1)*rscale)
    IF (osgbdy(iob,nbsfg3) > rmdich)                                           &
      odrbufp(ipos+10) = INT(osgbdy(iob,nbsfg3)*rscale)
    ipos = ipos + 10
  ENDDO

  Report: DO iob = nodrold(1)+1 , nodrtot(1)
! ==========================================

!   header of multi-level report (length = 12+istrlen)
!   ----------------------------
    indst = ipos + 1
    odrbufp(ipos+ 2) = imlid
    odrbufp(ipos+ 3) = momlhd(iob,nhitot)
    odrbufp(ipos+ 4) = momlhd(iob,nhjtot)
    odrbufp(ipos+ 5) = momlhd(iob,nhobtp)
    odrbufp(ipos+ 6) = momlhd(iob,nhcode)
    IF (omlhed(iob,nhalt) > rmdich) THEN
      odrbufp(ipos+7) = INT(omlhed(iob,nhalt ))
!   ELSE
!     odrbufp(ipos+7) = imdi
    ENDIF
    odrbufp(ipos+ 8) = INT(omlhed(iob,nhsurf))
    odrbufp(ipos+ 9) = INT(omlhed(iob,nhilon)*rscale)
    odrbufp(ipos+10) = INT(omlhed(iob,nhjlat)*rscale)
    odrbufp(ipos+11) = INT(omlhed(iob,nhtime)*rscale)
!   odrbufp(ipos+12) = momlhd(iob,nhnlev) * (1 - momlhd(iob,nhpass))
    odrbufp(ipos+12) = momlhd(iob,nhpass) + 10*momlhd(iob,nhnlev)
    ipos   = ipos + 12

    DO icl = 1 , istrlen
      odrbufp(ipos+icl) = ICHAR( yomlhd(iob) (icl:icl) )
    ENDDO
    ipos = ipos + istrlen

!   body of multi-level report (length=11*nlevel)
!   ---------------------------------------------

    Level: DO klev = 1,momlhd(iob,nhnlev)

!     DO icl = 1 , 10
!       odrbufp(ipos+icl) = imdi
!     ENDDO
      IF (omlbdy(iob,klev,nbtu  ) > rmdich)                                    &
        odrbufp(ipos+ 1) = INT(omlbdy(iob,klev,nbtu  )*rscale)
      IF (omlbdy(iob,klev,nbtv  ) > rmdich)                                    &
        odrbufp(ipos+ 2) = INT(omlbdy(iob,klev,nbtv  )*rscale)
      IF (omlbdy(iob,klev,nbtt  ) > rmdich)                                    &
        odrbufp(ipos+ 3) = INT(omlbdy(iob,klev,nbtt  )*rscale)
      IF (omlbdy(iob,klev,nbtrh ) > rmdich)                                    &
        odrbufp(ipos+ 4) = INT(omlbdy(iob,klev,nbtrh )*rscale)
      IF (omlbdy(iob,klev,nbtp  ) > rmdich)                                    &
        odrbufp(ipos+ 5) = INT(omlbdy(iob,klev,nbtp  )*rscale *c100r)
      IF (omlbdy(iob,klev,nbtz  ) > rmdich)                                    &
        odrbufp(ipos+ 6) = INT(omlbdy(iob,klev,nbtz  )*rscale)
      IF (ibit1( momlbd(iob,klev,nbterr),nvru ) == 1)                          &
        odrbufp(ipos+ 7) = INT(omlbdy(iob,klev,nbtuer)*rscale)
      IF (ibit1( momlbd(iob,klev,nbterr),nvrt ) == 1)                          &
        odrbufp(ipos+ 8) = INT(omlbdy(iob,klev,nbtter)*rscale)
      IF (ibit1( momlbd(iob,klev,nbterr),nvrq ) == 1)                          &
        odrbufp(ipos+ 9) = INT(omlbdy(iob,klev,nbtqer)*rscale)
      IF (ibit1( momlbd(iob,klev,nbterr),nvrz ) == 1)                          &
        odrbufp(ipos+10) = INT(omlbdy(iob,klev,nbtzer)*rscale)
      odrbufp(ipos+11) = momlbd(iob,klev,nbtlid)
      ipos   = ipos   + 11
    ENDDO Level
    odrbufp(indst) = ipos-indst+1

  ENDDO Report


  DO iob = nodrold(3)+1 , nodrtot(3)
! ==================================

!   header of GPS reports (length = 12 + istrlen)
!   ---------------------
    odrbufp(ipos+ 1) = 12 + istrlen + 10
    odrbufp(ipos+ 2) = igpid
    odrbufp(ipos+ 3) = mogphd(iob,nhitot)
    odrbufp(ipos+ 4) = mogphd(iob,nhjtot)
    odrbufp(ipos+ 5) = mogphd(iob,nhobtp)
    odrbufp(ipos+ 6) = mogphd(iob,nhcode)
    IF (ogphed(iob,nhalt) > rmdich) THEN
      odrbufp(ipos+7) = INT( ogphed(iob,nhalt) )
      odrbufp(ipos+8) = INT( ogphed(iob,nhalt) - ogphed(iob,nhsurf) )
!   ELSE
!     odrbufp(ipos+7) = imdi
!     odrbufp(ipos+8) = imdi
    ENDIF
    odrbufp(ipos+ 9) = INT(ogphed(iob,nhilon)*rscale)
    odrbufp(ipos+10) = INT(ogphed(iob,nhjlat)*rscale)
    odrbufp(ipos+11) = INT(ogphed(iob,nhtime)*rscale)
    odrbufp(ipos+12) = mogphd(iob,nhpass)
    ipos = ipos + 12

    DO icl = 1 , istrlen
      odrbufp(ipos+icl) = ICHAR( yogphd(iob) (icl:icl) )
    ENDDO
    ipos = ipos + istrlen

!   body of GPS reports (length = 8)
!   -------------------
!   DO icl = 1 , 8
!     odrbufp(ipos+icl) = imdi
!   ENDDO
    IF (ogpbdy(iob,nbgiwv) > rmdich)                                           &
      odrbufp(ipos+ 1) = INT(ogpbdy(iob,nbgiwv)*rscale)
    IF (ogpbdy(iob,nbgzwd) > rmdich)                                           &
      odrbufp(ipos+ 2) = INT(ogpbdy(iob,nbgzwd)*rscale)
    IF (ogpbdy(iob,nbgt  ) > rmdich)                                           &
      odrbufp(ipos+ 3) = INT(ogpbdy(iob,nbgt  )*rscale)
    IF (ogpbdy(iob,nbgrh ) > rmdich)                                           &
      odrbufp(ipos+ 4) = INT(ogpbdy(iob,nbgrh )*rscale)
    IF (ogpbdy(iob,nbgp  ) > rmdich)                                           &
      odrbufp(ipos+ 5) = INT(ogpbdy(iob,nbgp  )*rscale *c100r)
    IF (ogpbdy(iob,nbgzpd) > rmdich)                                           &
      odrbufp(ipos+ 6) = INT(ogpbdy(iob,nbgzpd)*rscale)
    IF (      (ogpbdy(iob,nbgtze) > rmdich)                                    &
        .AND. (ibit1( mogpbd(iob,nbgerr),nvriwv ) == 1))                       &
      odrbufp(ipos+ 7) = INT(ogpbdy(iob,nbgtze)*rscale)
    IF (ogpbdy(iob,nbgbia) > rmdich)                                           &
      odrbufp(ipos+ 8) = INT(ogpbdy(iob,nbgbia)*rscale)
    ipos = ipos + 8
  ENDDO


!-------------------------------------------------------------------------------
!  Section 2:  Gather ODR data at root node
!-------------------------------------------------------------------------------

  ! note that 'odrbufa' is allocated in obs_gather_buffer

  CALL obs_gather_buffer ( nodrln, nodrln+1, odrbufp, ilentot, odrbufa, 0      &
                         , num_compute, my_cart_id, icomm_cart, imp_integers )
! ======================

  lwr_now = (ilentot >= 1)

!-------------------------------------------------------------------------------
!  Section 3:  Print out gathered ODR data
!-------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN

    ! Open 'yuobsdr' if required
    IF ((lwr_now) .AND. (.NOT. lopen_odr)) THEN
      OPEN (nuodr ,FILE=yuobsdr ,FORM='FORMATTED',STATUS='UNKNOWN'             &
                                ,POSITION='APPEND',IOSTAT=istat)
      lopen_odr = .TRUE.
      IF (istat /= 0) yerrmsg = 'OPENING OF FILE yuobsdr FAILED'
      IF (istat /= 0) CALL model_abort (my_cart_id, 7005, yerrmsg, yroutine)
    ENDIF

    ipos = 0

    Get_next_odr:   DO
!   ==================
!     get record length
      irlen = odrbufa(ipos+1)
      IF (irlen == 0)                                          EXIT Get_next_odr
!     IF (irlen == imdi)                                       EXIT Get_next_odr
      indst = ipos + 2

      IF ((odrbufa(indst) == isgid) .OR. (odrbufa(indst) == igpid)) THEN

!       header of single-level or GPS report
!       ------------------------------------
        itot   = odrbufa(ipos+ 3)
        jtot   = odrbufa(ipos+ 4)
        iobtyp = odrbufa(ipos+ 5)
        icdtyp = odrbufa(ipos+ 6)
        istalt = odrbufa(ipos+ 7)
        ihadif = odrbufa(ipos+ 8)
        zoblon = odrbufa(ipos+ 9)*riscal
        zoblat = odrbufa(ipos+10)*riscal
        zobtim = odrbufa(ipos+11)*riscal
        ypass  = ' '
        IF (odrbufa(ipos+12) == 1) ypass  = 'M'
        IF (odrbufa(ipos+12) == 2) ypass  = 'P'
        IF (mobname == 1) THEN
          yobtyp = chobtp(iobtyp) (1:5)
        ELSE
          yobtyp = cma(i_cma( iobtyp,0 ))%name (1:5)
        ENDIF
        ipos = ipos + 12

        DO icl = 1 , ilstidp
          ystid (icl:icl) = CHAR( odrbufa(ipos+icl) )
        ENDDO
        ipos = ipos + istrlen

!       body of single-level report
!       ---------------------------
        IF (odrbufa(indst) == isgid) THEN
          DO icl = 1, 10
            zwrbuf (icl) = rmdi
            IF (odrbufa(ipos+icl) /= imdi)                                     &
              zwrbuf (icl) = odrbufa(ipos+icl) *riscal
          ENDDO
          ipos = ipos + 10

          WRITE( nuodr,'(A5,1X,A ,1X,A,F5.1,2F6.1,F5.2, F8.2,F7.0,I5,F5.1,F4.1 &
                     &,F5.2,F5.0,I5,"(",I7,",",I3,")",2(1X,F6.2),I4,I2,F6.2)') &
                 yobtyp, ystid, ypass                                          &
               ,(zwrbuf(icl), icl=1,6), istalt, (zwrbuf(icl), icl=7,10)        &
               , ihadif, itot, jtot, zoblon, zoblat, icdtyp, iobtyp, zobtim
!              , zuu ,zvv, ztt, zrh, zpp, zzz, istalt, zuer, zter, zqer, zzer  &
!              , ihadif, itot, jtot, zoblon, zoblat, icdtyp, iobtyp, zobtim

!       body of GPS report
!       ------------------
        ELSEIF (odrbufa(indst) == igpid) THEN
          DO icl = 1, 8
            zwrbuf (icl) = rmdi
            IF (odrbufa(ipos+icl) /= imdi)                                     &
              zwrbuf (icl) = odrbufa(ipos+icl) *riscal
          ENDDO
          ipos = ipos + 8

          WRITE( nuodr,'(A5,1X,A ,1X,A,F5.1,2F6.1,F5.2, F8.2,F7.1,I5,F5.1, 4X  &
                       &,F5.2,5X,I5,"(",I7,",",I3,")",2(1X,F6.2),I4,I2,F6.2)') &
                 yobtyp, ystid, ypass                                          &
               ,(zwrbuf(icl), icl=1,6), istalt, zwrbuf(7), zwrbuf(8)           &
               , ihadif, itot, jtot, zoblon, zoblat, icdtyp, iobtyp, zobtim
!              , ziwv ,zzwd, ztt, zrh, zpp, ztzd, istalt, ztze, zbia           &
!              , ihadif, itot, jtot, zoblon, zoblat, icdtyp, iobtyp, zobtim

        ENDIF

      ELSEIF (odrbufa(indst) == imlid)    THEN

!       header of multi-level report
!       ----------------------------
        itot   = odrbufa(ipos+ 3)
        jtot   = odrbufa(ipos+ 4)
        iobtyp = odrbufa(ipos+ 5)
        icdtyp = odrbufa(ipos+ 6)
        istalt = odrbufa(ipos+ 7)
        igporo = odrbufa(ipos+ 8)
        zoblon = odrbufa(ipos+ 9)*riscal
        zoblat = odrbufa(ipos+10)*riscal
        zobtim = odrbufa(ipos+11)*riscal
        ypass  = ' '
        IF (MOD( odrbufa(ipos+12), 10 ) == 1) ypass  = 'M'
        IF (MOD( odrbufa(ipos+12), 10 ) == 2) ypass  = 'P'
        IF (mobname == 1) THEN
          yobtyp = chobtp(iobtyp) (1:5)
        ELSE
          yobtyp = cma(i_cma( iobtyp,0 ))%name (1:5)
        ENDIF
        imllev = odrbufa(ipos+12) / 10
        ipos = ipos + 12

        DO icl = 1 , ilstidp
          ystid (icl:icl) = CHAR( odrbufa(ipos+icl) )
        ENDDO
        ipos = ipos + istrlen

!       print header
!       ------------
        WRITE( nuodr,'(A5,1X,A,1X,A,"  OBS/CODE TYPE:",I2,I4                   &
                                  &,"  STA/MOD HEIGHT:",I4,I5                  &
                                  &,"  LON=",F6.2," LAT=",F6.2                 &
                                  &,"  G.P.=(",I7,",",I3,")  HOUR=",F5.1)' )   &
               yobtyp, ystid ,ypass, iobtyp, icdtyp                            &
             , istalt, igporo, zoblon, zoblat, itot, jtot, zobtim

        WRITE( nuodr,'("     u      v     t    rh      p      z  "             &
                     &,"v-err t-er rh-er z-er lev")' )

!       body of multi-level report (length of one level = 11)
!       --------------------------
        DO klev = 1,imllev
!       ==================

          DO icl = 1, 10
            zwrbuf (icl) = rmdi
            IF (odrbufa(ipos+icl) /= imdi)                                     &
              zwrbuf (icl) = odrbufa(ipos+icl) *riscal
          ENDDO
          ilid = odrbufa(ipos+11)

          WRITE( nuodr,'(3F7.1,F5.2, F8.1,F7.0, 2F5.1,F5.2,F5.1,I4)' )         &
                (zwrbuf(icl), icl=1,10), ilid
!                zuu, zvv, ztt, zrh, zpp, zzz, zuer, zter, zqer, zzer, ilid

          ipos = ipos + 11
        ENDDO
      ELSE

        WRITE( yerrmsgl,'(''UNKNOWN OBS-TYPE IN ODR-BUFFER: ipos='',I6         &
                        &,'', odrbufa(indst)='',I6)' )    ipos, odrbufa(indst)
        CALL model_abort (my_cart_id, 11118, yerrmsgl, yroutine)

      ENDIF  ! single/multi-level branch
    ENDDO  Get_next_odr

    IF (lopen_odr) THEN
      CLOSE (nuodr)
      lopen_odr = .FALSE.
    ENDIF
  ENDIF ! my_cart_id==0

  DEALLOCATE ( odrbufp   , STAT = istat )
  DEALLOCATE ( odrbufa   , STAT = istat )

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_print_odr
!-------------------------------------------------------------------------------

END SUBROUTINE obs_cdf_print_odr


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_cdfin_print" for writing control output to files
!-------------------------------------------------------------------------------

SUBROUTINE obs_print_number_rep ( nodrtot , nodrold , imaxl , lwonl            &
                                , num_compute , my_cart_id , icomm_cart        &
                                , imp_integers )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfin_print" writes out statistics
!   of number of processed reports at each node to file number 'nupr'.
!
! Method:
!   If the program is run in parallel mode, the data to be printed are first
!   gathered from all PEs at the root node.
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

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)  :: &
    nodrtot (3)    ,& ! total number of multi-level / single-level / GPS reports
    nodrold (3)    ,& ! number of multi-level / single-level / GPS reports
                      !   before having read new data at the current timestep
    imaxl   (3)    ,& ! size (report dimension) of the (3 types of) ODR
    num_compute    ,& ! number of compute PEs
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_integers      ! INTEGER   type used for MPI

  LOGICAL                  , INTENT (IN)  :: &
    lwonl             ! .true. for the node (sub-domain) at which file with
                      !        the unit number 'nupr' is open

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    implcode       ,& ! for MPI error code
    ipe    , iobt     ! loop indices

  CHARACTER (LEN=20)       :: &
    yroutine          ! name of this subroutine
  CHARACTER (LEN=30)       :: &
    yerrmsg           ! error message

! Local (automatic) arrays:
! -------------------------

  INTEGER        ::  &
    numprtn (7)           ,& ! number of reports to be printed from current PE
    numprta (7,num_compute)  ! number of reports to be printed from all PEs
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_print_number_rep
!-------------------------------------------------------------------------------

  yroutine   = 'obs_print_number_rep'

  numprtn (1) = my_cart_id
  numprtn (2) = nodrold(2)
  numprtn (3) = nodrtot(2)
  numprtn (4) = nodrold(1)
  numprtn (5) = nodrtot(1)
  numprtn (6) = nodrold(3)
  numprtn (7) = nodrtot(3)

  IF (num_compute > 1) THEN
    CALL gather_values ( numprtn, numprta, 7, num_compute, imp_integers        &
                       , -1, icomm_cart, yerrmsg, implcode )
!   ==================
    IF (implcode /= 0)                                                         &
      CALL model_abort (my_cart_id, 11110, yerrmsg, yroutine, implcode)
  ELSE
    numprta (:,1) = numprtn(:)
  ENDIF

  IF (lwonl) THEN
    WRITE( nupr,'(" NUMBER OF SINGLE- AND MULTI-LEVEL AND GPS DATA TO BE"      &
                &," PRINTED:")' )
    WRITE( nupr,'(" NODE:      CART_ID  | NTOTSGO | NTOTSG | NTOTMLO | NTOTML" &
                &,                    " | NTOTGPO | NTOTGP")' )
    numprtn (:) = 0
    DO ipe = 1 , num_compute
      WRITE( nupr,'(13X,I4,7X,I5,4X,I5,5X,I5,4X,I5,5X,I5,4X,I5)' )             &
            (numprta(iobt,ipe), iobt=1,7)
      numprtn (2:7)   = numprtn(2:7) + numprta(2:7,ipe)
      numprta (2:7,1) = MAX( numprta(2:7,1) , numprta(2:7,ipe) )
    ENDDO
    WRITE( nupr,'(" MAX. LOCAL NUMBER:",3(5X,I5,4X,I5))')                      &
          (numprta(iobt,1), iobt=2,7)
    WRITE( nupr,'(" LOCAL ARRAY SIZE :",3(     14X,I5))') imaxl(2), imaxl(1)   &
                                                        , imaxl(3)
    WRITE( nupr,'(" GLOBAL NUMBER    :",3(5X,I5,4X,I5))')                      &
          (numprtn(iobt)  , iobt=2,7)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_print_number_rep
!-------------------------------------------------------------------------------

END SUBROUTINE obs_print_number_rep


!===============================================================================

ELEMENTAL INTEGER FUNCTION ibit1  ( invar, ibp )
  !---------------------------------------------
  INTEGER         , INTENT (IN)  ::  invar, ibp
  !---------------------------------------------------------------------------
  ! returns 1 bit at bit position 'ibp' of integer word 'invar'
  !---------------------------------------------------------------------------
  !
  ibit1 = IAND( ISHFT( invar,-ibp ), 1 )
  !
END FUNCTION ibit1

!-------------------------------------------------------------------------------

END MODULE src_obs_cdfin_print
