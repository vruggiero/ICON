!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_cdfin_blk

!-------------------------------------------------------------------------------
! Description:
!   This module contains routines which read and store the blacklist and
!   whitelist, and which produce a list of blacklisted vertical intervals
!   for each of a given array of reports.
!
!   Note: This module belongs to a group of COSMO-related modules for reading
!   ----  conventional data from ('cdfin') NetCDF observation input files.
!         It is shared between the COSMO model and the DACE program package !
!
! Method:
!   This module contains the following module procedures:
!    - obs_cdf_whitelist_local
!    - obs_cdf_blacklist_local
!    - obs_read_blacklist
!
!   It uses from:
!    - parallel_utilities:    - distribute_values
!    - environment:           - get_free_unit
!                             - model_abort
!
!   Data modules used:
!    - kind_parameters
!    - data_obs_lib_cosmo
!    - data_obs_cdfin
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
! V4_22        2012/01/31 Christoph Schraff
!  Initial release, extracted from module 'src_obs_proc_cdf' and adapted
!  (e.g. modified routine interfaces).
! V4_28        2013/07/12 Christoph Schraff
!  Improved comments only.
! V5_1         2014-11-28 Oliver Fuhrer
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

USE data_obs_lib_cosmo, ONLY :  &

! 1. General parameters
! ---------------------

    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 2. Variables and parameters obtained by calling 'obs_cdf_interface'
! -------------------------------------------------------------------

    ! variables related to parallelisation / domain decomposition
    num_compute    ,& ! number of compute PEs
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_integers   ,& ! INTEGER   type used for MPI
    imp_character  ,& ! CHARACTER type used for MPI

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    yucautn    ,& ! caution messages if too many obs for ODR size
                  !      incl. quality control flag for verification)
    nucautn    ,& ! caution messages if too many obs for current ODR size

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nairep        ! AIREP reports (all aircraft reports)

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

!         1.1   internal attributes of the different NetCDF input files
!               -------------------------------------------------------

    icdfinlen      ,& ! maximum length of NetCDF observation input file name

!         2.1   global blacklist and whitelists as read from file
!               -------------------------------------------------

    ilstid_blk     ,& ! assume 8-character station-IDs in Black-/Whitelist
    mxot_blk       ,& ! (max.) number of observation types in blacklist
    mx_wits        ,& ! max. number of whitelists (1 per obs code type)
    yblk_in        ,& ! file name of blacklist file (which must reside in
                      !   the same directory as the NetCDF obs input files !)
    nblack         ,& ! (total) number of blacklisted reports
    nwhite         ,& !  total  number of whitelisted stations
    n_wits         ,& ! actual  number of whitelisted stations
    kwit_obt       ,& ! observation type of whitelist
    kwit_cdt       ,& ! obs  code   type of whitelist
    kwit_frst      ,& ! index of first station of current whitelist
    kwit_last      ,& ! index of last  station of current whitelist
    iblk_frst      ,& ! index of first report of current obs type
    iblk_last      ,& ! index of last  report of current obs type
    yblk_id        ,& ! station identity of blacklisted stations
    ywit_id        ,& ! station identity of whitelisted stations
    iblk_len       ,& ! length of 'yblk_id'
    iwit_len       ,& ! length of 'ywit_id'
    iblk_pts       ,& ! number of points (wild cards) in 'yblk_id'
    iwit_pts       ,& ! number of points (wild cards) in 'ywit_id'
    iblk_obtyp     ,& ! observation type of blacklisted report
    rblk_pzlow     ,& ! lower  \  boundary (pressure) of blacklisted interval
    rblk_pzup      ,& ! upper  /                      for geopotential
    rblk_pvlow     ,& ! lower  \  boundary (pressure) of blacklisted interval
    rblk_pvup      ,& ! upper  /                      for wind
    rblk_ptlow     ,& ! lower  \  boundary (pressure) of blacklisted interval
    rblk_ptup      ,& ! upper  /                      for temperature
    rblk_pqlow     ,& ! lower  \  boundary (pressure) of blacklisted interval
    rblk_pqup      ,& ! upper  /                      for humidity

!         2.2   array of separate blacklists which are each for 1 local report
!               --------------------------------------------------------------

    maxintv        ,& ! max. number of vertical blacklist intervals per 1 stat.
    blk_loc           ! blacklists for local reports

USE data_obs_cdfin, ONLY :  &

!         7.    for reporting rejection of data: Output buffer, size and formats
!               ----------------------------------------------------------------

    outbuf     ,& ! buffer containing output for a single node
    nacout     ,& ! actual number of records stored in the output buffer
    nmxoln     ,& ! maximum length of output buffer
    istrej     ,& ! length of strings (station id) in output buffer
    nfmt24     ,& ! report (partly) blacklisted
    nfmt25        ! report not on whitelist

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    get_free_unit   ,& ! determines the next free file unit number
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
PUBLIC :: obs_cdf_whitelist_local,  obs_cdf_blacklist_local,  obs_read_blacklist

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

!===============================================================================
!+ Module procedure in "src_obs_cdfin_blk" for a blacklist for local reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_whitelist_local ( nrepl , jobtyp , jcdtyp , zobhr           &
                                   , ilstid_chk , ystidl , lblkw )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_blk" indicates for each
!   local report if it is blacklisted because it is missing on the whitelist
!   which is valid for the observation (code) type of the report.
!
! Method:
!   For each report, first determine the existence and number of the whitelist
!   which is valid for the observation (code) type of the report.
!   If such a whitelist exists then compare the station ID of the report
!   with the station ID's of the entries in that whitelist.
!   The character '.' is used as a one-character wild card in the entries.
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

  INTEGER        , INTENT (IN)  ::  &
    nrepl             ,& ! number of local reports (for each of which a list of
                         !   blacklisted vertical intervals is produced)
    jobtyp (nrepl)    ,& ! observation type of local reports
    jcdtyp (nrepl)    ,& ! obs  code   type of local reports
    ilstid_chk           ! length of (part of) station ID that can be checked

  REAL (KIND=wp) , INTENT (IN)  ::  &
    zobhr  (nrepl)       ! obs  time        of local reports

  CHARACTER (LEN=ilstid_blk) , INTENT (IN)  ::  &
    ystidl (nrepl)       ! station ID's of local reports

  LOGICAL        , INTENT (OUT) ::  &
    lblkw  (nrepl)       ! report blacklisted because station ID is missing on
                         !   the whitelist valid for its obs (code) type

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  LOGICAL        ::  &
    lwit                 ! .true. if current report is on the whitelist

  INTEGER        ::  &
    mwit   (nrepl)    ,& ! observation type of local reports
    iwit              ,& ! loop index of whitelists
    iwit0             ,& ! loop index of stations in a whitelist
    ilen              ,& ! length of output record
    irpl , jj            ! loop indices

  REAL (KIND=wp) ::  &
    rscale               ! scaling factor

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_whitelist_local
!-------------------------------------------------------------------------------

! ilstid_chk  =  MIN( ilstid , ilstid_blk )

  rscale  = 100._wp*(1._wp+epsy)

! check if the station ID (and obs type) matches entries on the whitelist
! -----------------------------------------------------------------------

! to each report, assign the whitelist valid for the same obs (code) type
! -----------------------------------------------------------------------
  mwit (:)  =  0
  DO iwit = 1 , n_wits
    DO irpl = 1 , nrepl
      IF (      (jobtyp(irpl) == kwit_obt(iwit))                               &
          .AND. (jcdtyp(irpl) == kwit_cdt(iwit)))  mwit (irpl) = iwit
    ENDDO
  ENDDO

! loop over the local reports
  DO irpl = 1 , nrepl
    lblkw (irpl) = .FALSE.

! if there is a whitelist for the report obs type,
! then indicate it if the station ID is missing on it
! ---------------------------------------------------
    IF (mwit(irpl) >= 1) THEN
      lblkw (irpl) = .TRUE.

! check the whitelist entries with the same obs type as the current report
      check_ids: DO iwit0 = kwit_frst(mwit(irpl)) , kwit_last(mwit(irpl))
! no   wild cards (".")
        IF (iwit_pts(iwit0) == 0) THEN
          IF (ystidl(irpl)(1:ilstid_chk) == ywit_id(iwit0) (1:ilstid_chk))     &
            lblkw (irpl) = .FALSE.

! only wild cards (".")
        ELSEIF (iwit_pts(iwit0) == iwit_len(iwit0)) THEN
          lblkw (irpl) = .FALSE.
! at least 1 non-".", i.e. : "A....   ", "A....B  ", "A....B. ", or ".A...B. "
! (note: to improve vectorisation, the following IF loop (which contains
!        an inner DO loop) could be run in an extra outer DO loop)
        ELSEIF (iwit_pts(iwit0) <  iwit_len(iwit0)) THEN
          lwit = .TRUE.
          DO jj = 1 , ilstid_chk
            IF (      (ywit_id(iwit0)(jj:jj) /= ystidl(irpl)(jj:jj))           &
                .AND. (ywit_id(iwit0)(jj:jj) /= '.'))      lwit = .FALSE.
          ENDDO
          IF (lwit)  lblkw (irpl) = .FALSE.
        ENDIF
        IF (.NOT. lblkw(irpl))                                    EXIT check_ids
      ENDDO check_ids

! prepare printing of control output (in file YUREJCT)
! ----------------------------------------------------

      IF (lblkw(irpl)) THEN
        ilen = 4 + 2*istrej
        IF (nacout+ilen <= nmxoln) THEN
          outbuf(nacout+ 1) = ilen
          outbuf(nacout+ 2) = nfmt25
          outbuf(nacout+ 3) = jcdtyp(irpl)
          outbuf(nacout+ 4) = NINT( zobhr(irpl)*rscale)
          DO jj = 1 , MIN( istrej , ilstid_chk )
            outbuf(nacout+4+jj) = ICHAR( ystidl(irpl) (jj:jj) )
          ENDDO
          DO jj = ilstid_chk + 1 , istrej
            outbuf(nacout+4+jj) = ICHAR( ' ' )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
      ENDIF
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_whitelist_local
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_whitelist_local


!===============================================================================
!+ Module procedure in "src_obs_cdfin_blk" for a blacklist for local reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_blacklist_local ( nrepl , jobtyp , ilstid_chk , ystidl )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_blk" produces a list of
!   blacklisted vertical intervals for each of an array of local reports.
!   These new lists allow for a more efficient use of the blacklist information
!   for multi-level reports than the original blacklist.
!
! Method:
!   Compare the station ID of each local report with the station ID's of the
!   entries in the blacklist which have the same observation type.
!   The character '.' is used as a one-character wild card in the entries.
!   If the station ID's match then the blacklisted vertical interval for any
!   variable is added to a list of blacklisted intervals which is assigned to
!   that report.
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

  INTEGER        , INTENT (IN)  ::  &
    nrepl             ,& ! number of local reports (for each of which a list of
                         !   blacklisted vertical intervals is produced)
    jobtyp (nrepl)    ,& ! observation type of local reports
    ilstid_chk           ! length of (part of) station ID that can be checked

  CHARACTER (LEN=ilstid_blk) , INTENT (IN)  ::  &
    ystidl (nrepl)       ! station ID's of local reports

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  LOGICAL                  ::  &
    lblk                 ! .true. if current report is on the blacklist

  INTEGER        ::  &
    ilen              ,& ! length of output record
    irpl , irblk , jj    ! loop indices

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_blacklist_local
!-------------------------------------------------------------------------------

! ilstid_chk  =  MIN( ilstid , ilstid_blk )

! check if the station ID (and obs type) matches entries on the blacklist
! -----------------------------------------------------------------------

! loop over the local reports (of a certain observation type)
  DO irpl = 1 , nrepl
    blk_loc(irpl)% ndim  =  0

! loop over the blacklist entries with the same obs type as the local report
    DO irblk = iblk_frst(jobtyp(irpl)) , iblk_last(jobtyp(irpl))

! no   wild cards (".")
      IF (iblk_pts(irblk) == 0) THEN
        lblk = (ystidl(irpl)(1:ilstid_chk) == yblk_id(irblk) (1:ilstid_chk))
! only wild cards (".")
      ELSEIF (iblk_pts(irblk) == iblk_len(irblk)) THEN
        lblk = .TRUE.
! at least 1 non-".", i.e. : "A....   ", "A....B  ", "A....B. ", or ".A...B. "
! (note: to improve vectorisation, the following IF loop (which contains
!        an inner DO loop) could be run in an extra outer DO loop)
      ELSEIF (iblk_pts(irblk) <  iblk_len(irblk)) THEN
        lblk = .TRUE.
        DO jj = 1 , ilstid_chk
          IF (      (yblk_id(irblk)(jj:jj) /= ystidl(irpl)(jj:jj))             &
              .AND. (yblk_id(irblk)(jj:jj) /= '.'))      lblk = .FALSE.
        ENDDO
      ENDIF

! if the station ID (and obs type) matches an entry on the blacklist, then
! complement the list of blacklisted vertical intervals for the current report
! ----------------------------------------------------------------------------

      IF (lblk) THEN
        IF ((rblk_pvlow(irblk) > epsy).AND. (blk_loc(irpl)%ndim < maxintv)) THEN
          blk_loc(irpl)% ndim                      =  blk_loc(irpl)% ndim  +  1
          blk_loc(irpl)% kvar(blk_loc(irpl)%ndim)  =  1
          blk_loc(irpl)% plow(blk_loc(irpl)%ndim)  =  rblk_pvlow(irblk)
          blk_loc(irpl)% pup (blk_loc(irpl)%ndim)  =  rblk_pvup (irblk)
        ENDIF
        IF ((rblk_pzlow(irblk) > epsy).AND. (blk_loc(irpl)%ndim < maxintv)) THEN
          blk_loc(irpl)% ndim                      =  blk_loc(irpl)% ndim  +  1
          blk_loc(irpl)% kvar(blk_loc(irpl)%ndim)  =  2
          blk_loc(irpl)% plow(blk_loc(irpl)%ndim)  =  rblk_pzlow(irblk)
          blk_loc(irpl)% pup (blk_loc(irpl)%ndim)  =  rblk_pzup (irblk)
        ENDIF
        IF ((rblk_ptlow(irblk) > epsy).AND. (blk_loc(irpl)%ndim < maxintv)) THEN
          blk_loc(irpl)% ndim                      =  blk_loc(irpl)% ndim  +  1
          blk_loc(irpl)% kvar(blk_loc(irpl)%ndim)  =  3
          blk_loc(irpl)% plow(blk_loc(irpl)%ndim)  =  rblk_ptlow(irblk)
          blk_loc(irpl)% pup (blk_loc(irpl)%ndim)  =  rblk_ptup (irblk)
        ENDIF
        IF ((rblk_pqlow(irblk) > epsy).AND. (blk_loc(irpl)%ndim < maxintv)) THEN
          blk_loc(irpl)% ndim                      =  blk_loc(irpl)% ndim  +  1
          blk_loc(irpl)% kvar(blk_loc(irpl)%ndim)  =  4
          blk_loc(irpl)% plow(blk_loc(irpl)%ndim)  =  rblk_pqlow(irblk)
          blk_loc(irpl)% pup (blk_loc(irpl)%ndim)  =  rblk_pqup (irblk)
        ENDIF
      ENDIF

! prepare printing of control output (in file YUREJCT)
! ----------------------------------------------------------------------------

      IF (lblk) THEN
        ilen = 11 + 2*istrej
        IF (      (nacout+ilen <= nmxoln)                                      &
            .AND. (     (jobtyp(irpl) /= nairep)                               &
                   .OR. (  MAX( rblk_pzlow(irblk), rblk_pvlow(irblk)           &
                              , rblk_ptlow(irblk), rblk_pqlow(irblk) )         &
                         > 9900.0_wp))) THEN
          outbuf(nacout+ 1) = ilen
          outbuf(nacout+ 2) = nfmt24
          outbuf(nacout+ 3) = jobtyp(irpl)
          outbuf(nacout+ 4) = NINT( rblk_pzlow(irblk) ) / 100
          outbuf(nacout+ 5) = NINT( rblk_pzup (irblk) ) / 100
          outbuf(nacout+ 6) = NINT( rblk_pvlow(irblk) ) / 100
          outbuf(nacout+ 7) = NINT( rblk_pvup (irblk) ) / 100
          outbuf(nacout+ 8) = NINT( rblk_ptlow(irblk) ) / 100
          outbuf(nacout+ 9) = NINT( rblk_ptup (irblk) ) / 100
          outbuf(nacout+10) = NINT( rblk_pqlow(irblk) ) / 100
          outbuf(nacout+11) = NINT( rblk_pqup (irblk) ) / 100
          DO jj = 1 , MIN( istrej , ilstid_chk )
            outbuf(nacout+11+jj) = ICHAR( ystidl(irpl) (jj:jj) )
          ENDDO
          DO jj = ilstid_chk + 1 , istrej
            outbuf(nacout+11+jj) = ICHAR( ' ' )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
      ENDIF
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_blacklist_local
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_blacklist_local


!===============================================================================
!+ Module procedure in "src_obs_cdfin_blk" for reading a blacklist / whitelist
!-------------------------------------------------------------------------------

SUBROUTINE obs_read_blacklist ( icdfdirlen , yblkdir )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_blk" reads the blacklist and
!   whitelist from a formatted ASCII-File (Filename = 'blklsttmp') and produces
!   some auxilliary lists to facilitate the further use of the blacklist.
!
! Method:
!   First, 1 node reads through the file and gets the number of reports in the
!   blacklist resp. whitelist, and the storage arrays are allocated
!   (dynamically) according to these numbers (by all nodes). Then, the nodes
!   reads and stores the black- and whitelists and broadcasts them to the other
!   nodes. Finally, the auxilliary lists are produced.
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

  INTEGER        , INTENT (IN)  ::  &
    icdfdirlen       ! max. length of name of directory where
                     !   NetCDF observation input files reside

  CHARACTER (LEN=*)        , INTENT (IN)  ::  &
    yblkdir          ! directory where the blacklist file resides

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    nublkin           ,& ! unit number of the blacklist (input) file
    ilenpath          ,& ! length of path name of blacklist file
    nblklen (4)       ,& ! number of reports in blacklist and whitelist
    nblkot  (mxot_blk),& ! number of listed reports for each observation type
    iwittyp           ,& ! format type of whitelist
    kobtyp  , kcdtyp  ,& ! observation type , observation code type
    jbl, jot,jch, iwit,& ! loop indices
    istat, nstat      ,& ! error indicators
    jerr , ierr          ! error indicators

  CHARACTER (LEN=MIN( icdfdirlen+icdfinlen, 255 ))  ::  &
    yblkfile             ! total path of blacklist (input) file
  CHARACTER (LEN=80)       :: &
    ytxt                 ! text line in blacklist file
  CHARACTER (LEN=20)       :: &
    yroutine             ! name of this subroutine
  CHARACTER (LEN=255)      :: &
    yerrmsl              ! error message

! Local arrays:
! ------------

  INTEGER        , ALLOCATABLE ::  &
    ibuf_blk    (:)   ,& ! buffer for 'distribute_values'
    iwit_obtyp  (:)   ,& ! observation type of whitelisted station
    iwit_cdtyp  (:)      ! observation code type of whitelisted station

  CHARACTER (LEN=100)      , ALLOCATABLE ::  &
    ybuf_blk    (:)      ! buffer for 'distribute_values' (must have LEN=100 !)
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_read_blacklist
!-------------------------------------------------------------------------------

! blacklist
! ---------
!!INTEGER , PARAMETER :: maxblk     = 5000
! INTEGER , PARAMETER :: ilstid_blk = 8 ! assume 8-char Stat-IDs in Blacklist
! INTEGER , PARAMETER :: mxot_blk   = 8 ! number of obs types in blacklist
! INTEGER :: nblack                ! (tot) number of blacklisted reports
! INTEGER :: nwhite                ! total number of whitelisted stations
! INTEGER :: mx_wits        = 20   ! max. number of whitelists (1 per code type)
! INTEGER :: n_wits                ! act. number of whitelists (1 per code type)
! INTEGER :: kwit_obt   (mx_wits)  ! observation type of whitelist
! INTEGER :: kwit_cdt   (mx_wits)  ! obs  code   type of whitelist
! INTEGER :: kwit_frst  (mx_wits)  ! index of first report of current whitelist
! INTEGER :: kwit_last  (mx_wits)  ! index of last  report of current whitelist
!!INTEGER :: kwit_len   (mx_wits)  ! number of whitelisted stations (per w-list)
! INTEGER :: iblk_frst (mxot_blk)  ! index of first report of current obs type
! INTEGER :: iblk_last (mxot_blk)  ! index of last  report of current obs type
! CHARACTER (LEN=ilstid_blk)                                                   &
!         :: yblk_id     (:) = ' ' ! station identity of blacklisted stations
! CHARACTER (LEN=ilstid_blk)                                                   &
!         :: ywit_id     (:) = ' ' ! station identity in all whitelists
! INTEGER :: iblk_len    (:) = 0   ! length of 'yblk_id'
! INTEGER :: iblk_pts    (:) = 0   ! number of points (wild cards) in 'yblk_id'
! INTEGER :: iblk_obtyp  (:) = 0   ! observation type of blacklisted report
! REAL    :: rblk_pzlow  (:) = 0   ! lower  \  boundary (pressure) of black-
! REAL    :: rblk_pzup   (:) = 0   ! upper  /  listed interval for geopotent.
! REAL    :: rblk_pvlow  (:) = 0   ! lower  \  boundary (pressure) of black-
! REAL    :: rblk_pvup   (:) = 0   ! upper  /  listed interval for wind
! REAL    :: rblk_ptlow  (:) = 0   ! lower  \  boundary (pressure) of black-
! REAL    :: rblk_ptup   (:) = 0   ! upper  /  listed interval for temperature
! REAL    :: rblk_pqlow  (:) = 0   ! lower  \  boundary (pressure) of black-
! REAL    :: rblk_pqup   (:) = 0   ! upper  /  listed interval for humidity
! INTEGER :: iwit_len    (:) = 0   ! length of 'ywit_id'
! INTEGER :: iwit_pts    (:) = 0   ! number of points (wild cards) in 'ywit_id'
! INTEGER :: iwit_obtyp  (:) = 0   ! observation type of whitelisted station
! INTEGER :: iwit_cdtyp  (:) = 0   ! observat. code type of whitelisted station
!
! INTEGER :: iwit_len  (:,:) = 0   ! length of 'ywit_id'
! INTEGER :: iwit_pts  (:,:) = 0   ! number of points (wild cards) in 'ywit_id'
!!INTEGER :: iwit_obtyp(:,:) = 0   ! observation type of whitelisted station
!!INTEGER :: iwit_cdtyp(:,:) = 0   ! observat. code type of whitelisted station
!-------------------------------------------------------------------------------

  yroutine = 'obs_read_blacklist'

  nblklen   (:) = 0
  iwittyp       = 0
  nblack        = 0
  nwhite        = 0
  n_wits        = 0
  kwit_frst (:) = 1
  kwit_last (:) = 1
  kwit_obt  (:) = 0
  kwit_cdt  (:) = 0

! open the blacklist file
! -----------------------

  IF (my_cart_id == 0) THEN
!   get full path of blacklist file (which must reside in same directory
!                                    as the NetCDF observation input files!)
!   yblk_in  = 'blklsttmp'
    yblkfile = yblk_in
    ilenpath = LEN_TRIM(yblkdir)
    IF (ilenpath >= 1) THEN
      IF (yblkdir(ilenpath:ilenpath) == '/')  ilenpath = ilenpath - 1
      yblkfile = yblkdir(1:ilenpath) // '/' // yblkfile(1:LEN_TRIM(yblkfile))
    ENDIF
!   open file
    PRINT *,' '
    PRINT *,' open and read BLACKLIST / WHITELIST from '                       &
           , yblkfile(1:LEN_TRIM(yblkfile))
    CALL get_free_unit (nublkin)
    OPEN ( nublkin, FILE=TRIM(yblkfile), FORM='FORMATTED', STATUS='OLD'        &
         , IOSTAT=nstat)
    IF (nstat /= 0) THEN
      WRITE( yerrmsl,'("OPENING OF BLACKLIST FILE ",A," FAILED")' )            &
             yblkfile(1:LEN_TRIM(yblkfile))
      CALL model_abort (my_cart_id, 1009, yerrmsl, yroutine)
    ENDIF
    REWIND nublkin

! get the blacklist length and number and lengths of whitelists
! -------------------------------------------------------------
!   nblack = 0
!   nwhite = -1
!   read_blklen:  DO
!     READ( nublkin,'(A)',IOSTAT=istat ) ytxt
!     IF (istat < 0)                                            EXIT read_blklen
!     IF ((ytxt(1:9) /= 'WHITELIST') .AND. (nwhite == -1)) THEN
!            nblack = nblack + 1
!     ELSE ; nwhite = nwhite + 1 ; ENDIF
!   ENDDO read_blklen

!   (the first line is a header line, therefore start counting with -1)
    nblack = -1

    read_blklen:  DO
      ytxt = '                                                            '
      READ( nublkin,'(A)',IOSTAT=istat ) ytxt
!     end of file:
      IF (istat < 0)                                            EXIT read_blklen
      IF ((ytxt(1:9) /= 'WHITELIST') .AND. (iwittyp == 0)) THEN
        nblack = nblack + 1
! new type of whitelist: new code type
      ELSEIF ((ytxt(1:9) == 'WHITELIST') .AND. (LEN_TRIM(ytxt) >= 10)) THEN
        iwittyp = 2
        n_wits  = n_wits + 1
        READ( ytxt,'(9X,1X,I2,1X,I5)' )  kwit_obt(n_wits), kwit_cdt(n_wits)
        IF (n_wits >= 2)  kwit_frst (n_wits)  =  kwit_last(n_wits-1) + 1
        kwit_last (n_wits)  =  kwit_frst(n_wits) - 1
! new type of whitelist: old code type
      ELSEIF (iwittyp == 2) THEN
        kwit_last (n_wits)  =  kwit_last(n_wits) + 1
        nwhite  =  nwhite + 1
! old type of whitelist determined
      ELSEIF (ytxt(1:9) == 'WHITELIST') THEN
        iwittyp = 1
! old type of whitelist:
      ELSE
        READ( ytxt,'(9X,I2,1X,I5)' )  kobtyp, kcdtyp
        IF (     (kobtyp /= kwit_obt(MAX(n_wits,1)))                           &
            .OR. (kcdtyp /= kwit_cdt(MAX(n_wits,1)))) THEN
          n_wits  = n_wits + 1
          kwit_obt (n_wits)  =  kobtyp
          kwit_cdt (n_wits)  =  kcdtyp
          IF (n_wits >= 2)  kwit_frst (n_wits)  =  kwit_last(n_wits-1) + 1
          kwit_last (n_wits)  =  kwit_frst(n_wits) - 1
        ENDIF
        kwit_last (n_wits)  =  kwit_last(n_wits) + 1
        nwhite  =  nwhite + 1
      ENDIF
    ENDDO read_blklen
    istat       = 0
    nblklen (1) = nblack
    nblklen (2) = nwhite
    nblklen (3) = n_wits
    nblklen (4) = iwittyp
  ENDIF

! broadcast the lengths of the blacklist and of the whitelist
! -----------------------------------------------------------

  IF (num_compute > 1) THEN

    CALL distribute_values (nblklen, 4, 0, imp_integers, icomm_cart, jerr)
!   ======================
    nblack  = nblklen(1)
    nwhite  = nblklen(2)
    n_wits  = nblklen(3)
    iwittyp = nblklen(4)
    CALL distribute_values (kwit_frst, mx_wits, 0, imp_integers,icomm_cart,jerr)
    CALL distribute_values (kwit_last, mx_wits, 0, imp_integers,icomm_cart,jerr)
    CALL distribute_values (kwit_obt , mx_wits, 0, imp_integers,icomm_cart,jerr)
    CALL distribute_values (kwit_cdt , mx_wits, 0, imp_integers,icomm_cart,jerr)
!   ======================
  ENDIF

! allocate arrays for blacklist / whitelist
! -----------------------------------------

! long-term storage
  ALLOCATE ( yblk_id    (nblack+1) , STAT=istat )
  ALLOCATE ( iblk_len   (nblack+1) , STAT=istat )
  ALLOCATE ( iblk_pts   (nblack+1) , STAT=istat )
  ALLOCATE ( iblk_obtyp (nblack+1) , STAT=istat )
  ALLOCATE ( rblk_pzlow (nblack+1) , STAT=istat )
  ALLOCATE ( rblk_pzup  (nblack+1) , STAT=istat )
  ALLOCATE ( rblk_pvlow (nblack+1) , STAT=istat )
  ALLOCATE ( rblk_pvup  (nblack+1) , STAT=istat )
  ALLOCATE ( rblk_ptlow (nblack+1) , STAT=istat )
  ALLOCATE ( rblk_ptup  (nblack+1) , STAT=istat )
  ALLOCATE ( rblk_pqlow (nblack+1) , STAT=istat )
  ALLOCATE ( rblk_pqup  (nblack+1) , STAT=istat )

  ALLOCATE ( ywit_id    (nwhite+1) , STAT=istat )
  ALLOCATE ( iwit_len   (nwhite+1) , STAT=istat )
  ALLOCATE ( iwit_pts   (nwhite+1) , STAT=istat )

! temporary buffers
  ALLOCATE ( ibuf_blk (9*nblack+1) , STAT=istat )

! read blacklisted stations with pressure intervals, read whitelist
! -----------------------------------------------------------------

  IF (my_cart_id == 0) THEN
    ALLOCATE ( iwit_obtyp (nwhite+1) , STAT=istat )
    ALLOCATE ( iwit_cdtyp (nwhite+1) , STAT=istat )

    REWIND nublkin
    IF (nblack > 0) THEN
      READ( nublkin,'(A)',IOSTAT=istat ) ytxt
      DO jbl = 1 , nblack
        READ( nublkin, '(A8,1X,I1,8(1X,I5))' )                                 &
              yblk_id(jbl), (ibuf_blk(jot), jot=jbl*9-8,jbl*9)
      ENDDO
    ENDIF
    IF (nwhite > 0) THEN
! read title line of (first) whitelist
      read_wls:  DO
        READ( nublkin,'(A)',IOSTAT=istat ) ytxt
        IF ((ytxt(1:9) == 'WHITELIST') .OR. (istat < 0))           EXIT read_wls
      ENDDO read_wls
! read whitelists
      DO iwit = 1 , n_wits
        IF (MIN(iwittyp,iwit) >= 2) READ( nublkin,'(A)',IOSTAT=istat ) ytxt
        DO jbl = kwit_frst(iwit) , kwit_last(iwit)
          READ( nublkin, '(A8,1X,I2,1X,I5)' )                                  &
                ywit_id(jbl), iwit_obtyp(jbl), iwit_cdtyp(jbl)
        ENDDO
      ENDDO
    ENDIF

! close the blacklist file
! ------------------------

    CLOSE (nublkin, IOSTAT=nstat)

! check observation (code) types in whitelists
! --------------------------------------------

    DO iwit = 1 , n_wits
      IF (     (iwit_obtyp(kwit_frst(iwit)) /= kwit_obt(iwit))                 &
          .OR. (iwit_obtyp(kwit_last(iwit)) /= kwit_obt(iwit))                 &
          .OR. (iwit_cdtyp(kwit_frst(iwit)) /= kwit_cdt(iwit))                 &
          .OR. (iwit_cdtyp(kwit_last(iwit)) /= kwit_cdt(iwit))) THEN
        OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'        &
                                   , POSITION='APPEND', IOSTAT=nstat)
        WRITE( nucautn,'("CAUTION: INCONSISTENCIES IN WHITELIST",10I4)' )      &
               n_wits, iwit, kwit_frst(iwit), kwit_last(iwit)                  &
             , iwit_obtyp(kwit_frst(iwit)), iwit_obtyp(kwit_last(iwit))        &
             , kwit_obt(iwit)             , iwit_cdtyp(kwit_frst(iwit))        &
             , iwit_cdtyp(kwit_last(iwit)), kwit_cdt(iwit)
        CLOSE (nucautn)
      ENDIF
    ENDDO
    DEALLOCATE ( iwit_obtyp , STAT=istat )
    DEALLOCATE ( iwit_cdtyp , STAT=istat )
  ENDIF

! broadcast the blacklist and the whitelist to all nodes
! and fill the long-term storage arrays
! ------------------------------------------------------

  IF ((num_compute > 1) .AND. (nblack > 0)) THEN
! arrays of characters in distribute_values must have (LEN=100)
    ALLOCATE ( ybuf_blk (nblack) , STAT=istat )
    DO jbl = 1 , nblack
      ybuf_blk (jbl)  =  yblk_id (jbl)
    ENDDO

    CALL distribute_values (ybuf_blk,   nblack, 0,imp_character,icomm_cart,ierr)
    CALL distribute_values (ibuf_blk, 9*nblack, 0,imp_integers ,icomm_cart,ierr)
!   ======================

    DO jbl = 1 , nblack
      yblk_id (jbl)  =  ybuf_blk(jbl) (1:ilstid_blk)
    ENDDO
    DEALLOCATE ( ybuf_blk , STAT=istat )
  ENDIF

  DO jbl = 1 , nblack
    iblk_obtyp (jbl)  =        ibuf_blk(jbl*9-8)
    rblk_pzlow (jbl)  =  REAL( ibuf_blk(jbl*9-7) , wp ) * 100.0_wp
    rblk_pzup  (jbl)  =  REAL( ibuf_blk(jbl*9-6) , wp ) * 100.0_wp
    rblk_pvlow (jbl)  =  REAL( ibuf_blk(jbl*9-5) , wp ) * 100.0_wp
    rblk_pvup  (jbl)  =  REAL( ibuf_blk(jbl*9-4) , wp ) * 100.0_wp
    rblk_ptlow (jbl)  =  REAL( ibuf_blk(jbl*9-3) , wp ) * 100.0_wp
    rblk_ptup  (jbl)  =  REAL( ibuf_blk(jbl*9-2) , wp ) * 100.0_wp
    rblk_pqlow (jbl)  =  REAL( ibuf_blk(jbl*9-1) , wp ) * 100.0_wp
    rblk_pqup  (jbl)  =  REAL( ibuf_blk(jbl*9  ) , wp ) * 100.0_wp
  ENDDO
  DEALLOCATE ( ibuf_blk , STAT=istat )

  IF ((num_compute > 1) .AND. (nwhite > 0)) THEN
    ALLOCATE ( ybuf_blk (nwhite) , STAT=istat )
    DO jbl = 1 , nwhite
      ybuf_blk (jbl)  =  ywit_id (jbl)
    ENDDO

!   CALL distribute_values (iwit_obtyp,nwhite, 0,imp_integers ,icomm_cart, ierr)
    CALL distribute_values (ybuf_blk  ,nwhite, 0,imp_character,icomm_cart, ierr)
!   ======================

    DO jbl = 1 , nwhite
      ywit_id (jbl)  =  ybuf_blk(jbl) (1:ilstid_blk)
    ENDDO
    DEALLOCATE ( ybuf_blk , STAT=istat )
  ENDIF

! get length of station ID's and the number of wild cards ('.')
! -------------------------------------------------------------

  DO jbl = 1 , nblack
    iblk_len (jbl) = LEN_TRIM(yblk_id(jbl))
    iblk_pts (jbl) = 0
    DO jch = 1 , ilstid_blk
      IF (yblk_id(jbl)(jch:jch) == '.')  iblk_pts (jbl) = iblk_pts(jbl) + 1
    ENDDO
  ENDDO
  DO jbl = 1 , nwhite
    iwit_len (jbl) = LEN_TRIM(ywit_id(jbl))
    iwit_pts (jbl) = 0
!   DO jch = 1 , ilstid_wls
    DO jch = 1 , ilstid_blk
      IF (ywit_id(jbl)(jch:jch) == '.')  iwit_pts (jbl) = iwit_pts(jbl) + 1
    ENDDO
  ENDDO

! get indices of first and last report for each observation type in blacklist
! ---------------------------------------------------------------------------

! get number of stations as a function of observation type
  nblkot = 0
  DO jbl = 1 , nblack
    nblkot (iblk_obtyp(jbl))  =  nblkot(iblk_obtyp(jbl)) + 1
  ENDDO
! get indices of first and last report for each observation type
  iblk_frst (1)  =  1
  iblk_last (1)  =  nblkot(1)
  DO jot = 2 , mxot_blk
    iblk_frst (jot)  =  iblk_last(jot-1) + 1
    iblk_last (jot)  =  iblk_last(jot-1) + nblkot(jot)
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_read_blacklist
!-------------------------------------------------------------------------------
END SUBROUTINE obs_read_blacklist


!-------------------------------------------------------------------------------

END MODULE src_obs_cdfin_blk
