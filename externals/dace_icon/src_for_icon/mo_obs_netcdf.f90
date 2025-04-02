!
!+ Routines to read observations from NetCDF files
!
MODULE mo_obs_netcdf
!
! Description:
!   Routines to read observations from NetCDF files generated
!   by 'BUFR2NETCDF'.
!
! Current Maintainer: DWD, Harald Anlauf, Alexander Cress
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_2         2008/12/04 Andreas Rhodin
!  changes for SATEM from IASI
! V1_4         2009/03/26 Andreas Rhodin
!  Read Scatterometer NetCDF data in parallel
! V1_6         2009/06/10 Andreas Rhodin
!  check for valid entry in MSEC (radio-occultations)
! V1_7         2009/08/24 Andreas Rhodin
!  increase dimension of array 'dbkz': 10000 -> 11000
! V1_8         2009/12/09 Andreas Rhodin
!  Add BUFR2NetCDF input routine for radio occultations
! V1_9         2010/04/20 Harald Anlauf
!  scan_obs_netcdf: write filename before crashing (GFS problem?)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  Slightly reduce verbosity; some cosmetic changes to diagnostic output
! V1_13        2011/11/01 Harald Anlauf
!  cleanup; improved diagnostics
! V1_20        2012-06-18 Harald Anlauf
!  Framework for arbitrary number of input feedback files
!  detect NetCDF observation files (for satpp)
! V1_22        2013-02-13 Harald Anlauf
!  changes for GPSRO processing
! V1_24        2013/04/26 Andreas Rhodin
!  bugfix (check if a file is not read on a specific PE could fail)
! V1_29        2014/04/02 Andreas Rhodin
!  read and process bufr2netcdf ZTD data
!  change handling for the case that station Id exceeds predefined length
! V1_31        2014-08-21 Harald Anlauf
!  guess_dbkz: handle GFZ's RO test data
!  changes for ECMWF SYNOP (BUFR2)NetCDF input
! V1_33        2014-09-19 Harald Anlauf
!  dismiss 'empty' NetCDF files from corrupted BUFR->NetCDF conversion
! V1_35        2014-11-07 Alexander Cress
!  changes for TEMP BUFR reports
! V1_37        2014-12-23 Harald Anlauf
!  new BUFR codes in default rules tables; guess_dbkz: fallbacks for SYNOP/TEMP
! V1_44        2015-09-30 Harald Anlauf
!  move type t_bufr_inv to mo_t_obs to simplify dependencies
! V1_46        2016-02-05 Harald Anlauf
!  read_obs_netcdf: fix DRIBU station ids
! V1_48        2016-10-06 Harald Anlauf
!  Add CMA as processing center for GPSRO, handle FY-3C/GNOS
! V1_51        2017-02-24 Andreas Rhodin
!  option to read GPSGB netcdf input in parallel
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Gerhard Paul   DWD  2008  original source
! Harald Anlauf  DWD  2008  optimizations and fixes for SX8
!===============================================================
!-------------
! Modules used
!-------------
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_kind,         only: i8                ! kind parameters
  use mo_exception,    only: finish            ! abort routine
  use mo_system,       only: flush             ! Flush I/O buffer
  use mo_mpi_dace,     only: dace,            &! MPI group info
                             p_bcast           ! generic broadcast routine
  use mo_run_params,   only: path_file,       &! concatenate path/filename
                             host              ! host name
  use mo_time,         only: time_yyyymmdd_hhmmss , &!
                             operator(<),     &! compare times
                             operator(<=),    &! compare times
                             zero_time,       &! default initialisation value
                             cyyyymmddhhmmss   ! derive string from time
  use mo_p_output,     only: oline,           &! output line buffer
                             iol,             &! number of next line to write
                             nextline,        &! increment line number
                             flush_buf         ! write buffer
  !---------------------------------------
  ! DACE observation data types and tables
  !---------------------------------------
  use mo_t_use,        only: t_use,           &! status variable data type
                             CHK_NOIMPL,      &! reason for dismissal
                             CHK_NONE,        &! default status
                             STAT_FORGET       !
  use mo_t_obs,        only: t_obs,           &! observation data type
                             t_head,          &! observation data type
                             derive_dbkz,     &! derive DBKZ if not present
                             source,          &! table of observation files
                             n_source,        &!     no of Report source files to read
!                            m_source,        &! max no of Report source files to read
                             bufr_inv,        &! file inventories
                             t_bufr_inv,      &! derived type definition
                             FT_NETCDF,       &! flag for NetCDF file
                             FT_MISSING,      &! flag for NetCDF file
                             FT_UNKNOWN,      &! flag for NetCDF file
                             netcdf_verb       ! verbosity of NetCDF decoding
  use mo_wigos,        only: wsi_mode,        &! WIGOS station id mode
                             wsi_prefix,      &! Prefix of "fake station id"
                             wsi_encode,      &! Derive statid,stathash from t_wsi
                             wsi_decode,      &! Derive WSI from given string
                             wsi_from_statid   ! Derive WSI from short station name?
  !-----------------------------------
  ! definition and conversion of codes
  !-----------------------------------
  use mo_obstypes,     only: t_obsid,         &! observation id table entry
                             obstype_dbkz,    &! derive obsids from dbkz
                             obstype_bufr,    &! derive obsids from BUFR type
                             dbkz_bufr         ! derive dbkz   from BUFR type
  use mo_wmo_tables,   only: WMO0_ECMWF,      &! generating center
                             WMO0_DWD,        &! DWD
                             WMO0_EUMET,      &! EUMETSAT
                             WMO0_NASA,       &! NASA (JPL)
                             WMO0_NCAR,       &! NCAR
                             WMO0_NOAA,       &! NOAA
                             WMO0_UCAR,       &! UCAR (CDAAC)
                             WMO0_CMA,        &! CMA
                             WMO0_RSMC,       &! CMA
                             WMO0_DMI,        &! DMI / ROMSAF
                             WMO0_SPIRE,      &! Spire     (commercial)
                             WMO0_GEOOPT,     &! GeoOptics (commercial)
                             WMO0_PLANET       ! PlanetiQ  (commercial)
! use mo_dwd_tables,   only:                  &! Datenbank-Kennziffer:
!                            DK_SATOB_2,      &! Satellitenbeob, SATOB, Section 2
!                            DK_SATOB_3,      &! Satellitenbeob, SATOB, Section 3
!                            DK_AMV_EUMETSAT, &!
!                            DK_AMV_GOES,     &!
!                            DK_AMV_MODIS,    &!
!                            DK_AMV_FY_X,     &!
!                            DK_AMV_MTV,      &!
!                            DK_AMV_NOAA       !
  use mo_obs_tables,   only: idb_dbk,         &! index in table rept_stat
                             rept_char,       &! observation type characteristics
                             rept_use,        &! observation type usage
                             dism_report_0     ! dismiss report
  use mo_fdbk_tables,  only: n_ot,            &! number of observation types
                             OT_SYNOP,        &! report type identifier
                             OT_SATOB,        &! report type identifier
                             OT_TEMP,         &! report type identifier
                             OT_PILOT,        &! report type identifier
                             OT_AIREP,        &! report type identifier
                             OT_DRIBU,        &! report type identifier
                             OT_PAOB,         &! report type identifier
                             OT_SCATT,        &! report type identifier
                             OT_SATEM,        &! report type identifier
                             OT_GPSRO,        &! report type identifier
                             OT_GPSGB,        &! report type identifier
                             OT_SOIL,         &! report type identifier
                             OT_WLIDAR         ! report type identifier
  !------------------------------------------------
  ! observation type specific NetCDF input routines
  !------------------------------------------------
  use mo_temp,         only: read_temp_netcdf  ! read TEMP   observations from netCDF
  use mo_synop,        only: read_synop_netcdf ! read SYNOP  observations from netCDF
  use mo_scatt,        only: read_scatt_netcdf ! read SCATT  observations from netCDF
  use mo_amv,          only: read_amv_netcdf   ! read AMV    observations from netCDF
  use mo_airep,        only: read_airep_netcdf ! read AIREP  observations from netCDF
  use mo_occ,          only: read_gpsro_netcdf ! read GPSRO  observations from netCDF
  use mo_std,          only: read_std_netcdf   ! read GPSGB  observations from netCDF
  use mo_wlidar,       only: read_wlidar_netcdf! read WLIDAR observations from netCDF
  use mo_soil_obs,     only: read_soil_netcdf  ! read SOIL   observations from netCDF
  !-----------------------------------------------------------
  ! variables from BUFR/NetCDF valid for all observation types
  !-----------------------------------------------------------
  use mo_head_netcdf,  only:                  &!
                             dimids_max,      &! Maximum number of dimensions ids
                             ncid,            &! NetCDF file id
                             edition,         &! BUFR edition
                             s1date,          &! nominal (synoptic) date
                             s1time,          &! nominal (synoptic) time
                             s1cat,           &! data category
                             s1cats,          &! data sub category
                             s1catls,         &! local data sub category
                             s1cent,          &! data centre
                             s1cents,         &! data sub centre
                             s1updat,         &! update sequence no.
                                               !    expansion in NetCDF file BUFR- data section2
                             s2ikz,           &! DWD-internal classifier
                             s2dcdate,        &! decoding date
                             s2dctime,        &! decoding time
                                               !    expansion in NetCDF file BUFR- data section4
                             mlah,mloh,       &! latitude longitude
                             mjjj,mmm,myy,    &! year month day
                             mgg,ngg,msec,    &! hour minute second
                             msec_r,          &! second ; radiooccultations in SSsss
                                               ! for exchange --------------------------------
                             stime,           &! header observation time (section1)
                             db_time,         &! data bank time
                             obs_time,        &! observation time
                             istidn,          &! WMO numeric station number combined
                             ilstidn,         &! length of station identifier in characters
                             ystidn,          &! any type of station identifier as variable
                                               ! for expansion -------------------------------
                             bdy_date,        &! observation date: CCYYMMDD
                             bdy_time,        &! observation time: HHMMSSsss; sec.123
                             mii, niii,       &! WMO block number WMO station number
                             istidb,          &! WMO buoy/platform identifier
                             istids,          &! WMO satellite identifier
                             lwsi,            &! WIGOS station identifier valid
                             wsi,             &! WIGOS station id stored representation
                             wsihash,         &! Station id hash
                             ystid,           &! any type of station identifier as array
                             lvarid            ! variable for station ID exists in NetCDF file
  use mo_head_netcdf,                         &! expansion in NetCDF file BUFR- data section1
                                               ! for netCDF variables id ---------------------
                                               ! expansion in NetCDF file BUFR- data section1
                      only:                   &!
                             !------------------------------
                             ! missing values in NetCDF file
                             !------------------------------
                             imissing,        &! NetCDF _FillValue for integer
                             rmissing,        &! NetCDF _FillValue for reals
                             ymissing,        &! NetCDF _FillValue for character
                             iascii_ymissing   ! number in ascii sorting sequence for ymissing value
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,       only: nf90_open,             &! open   NetCDF file
                          nf90_close,            &! close  NetCDF file
                          nf90_Inquire,          &!
                          nf90_Inquire_Dimension,&!
                          nf90_inq_dimid,        &!
                          nf90_inq_varid,        &!
                          nf90_get_var,          &!
                          nf90_get_att,          &!
                          nf90_strerror,         &!
                          NF90_FILL_CHAR,        &!
                          NF90_NOERR,            &!
                          NF90_NOWRITE            ! NetCDF read flag
  implicit none
  !----------------
  ! Public entities
  !----------------
  private
  public :: scan_obs_netcdf  ! scan the NetCDF file and identify the content
  public :: read_obs_netcdf  ! read variables valid for all observation types

  !-------------------------------------------------
  ! variable ID's in NetCDF file for
  ! no expansion  in NetCDF file BUFR- data section0 (total length in Bytes)
  !    expansion  in NetCDF file BUFR- data section1
  ! variable ID's in NetCDF file in    data section1  (date time category centre)
  !-------------------------------------------------
  integer              :: varid_edition  ! BUFR edition
  integer              :: varid_s1date   ! nominal (synoptic) date
  integer              :: varid_s1time   ! nominal (synoptic) time
  integer              :: varid_s1cat    ! data category
  integer              :: varid_s1cats   ! data sub category  (undefined in bufr edition 3)
  integer              :: varid_s1catls  ! local data sub category
  integer              :: varid_s1cent   ! data centre
  integer              :: varid_s1cents  ! data sub centre
  integer              :: varid_s1updat  ! update sequence no.
  !-------------------------------------------------
  !    expansion  in NetCDF file BUFR- data section2
  ! variable ID's in NetCDF file  in   data section2 (kennzahl decoding date time)
  !-------------------------------------------------
  integer              :: varid_s2ikz    ! DWD-internal classifier
  integer              :: varid_s2dcdate ! decoding date
  integer              :: varid_s2dctime ! decoding time
  !-------------------------------------------------
  ! variable ID's in NetCDF file within data section3 (number of subsets,
  !                                                    still not defined in netCDF)
  ! variable ID's in NetCDF file within data section4
  !-------------------------------------------------
  integer              :: varid_mii      ! NetCDF variable  id for  WMO-Block-Nummer
  integer              :: varid_niii     ! NetCDF variable  id for  WMO-Stationsnummer
  integer              :: varid_mabnn    ! NetCDF variable  id for  WMO buoy/platform identifier
  integer              :: varid_mi1i2    ! NetCDF variable  id for  WMO satellite identifier
  integer              :: varid_ystid    ! NetCDF variable  id for  any type of character station id
  integer              :: varid_mlah     ! NetCDF variable  id for  Latitude
  integer              :: varid_mloh     ! NetCDF variable  id for  Longitude
  integer              :: varid_mjjj     ! NetCDF variable  id for  year
  integer              :: varid_mmm      ! NetCDF variable  id for  month
  integer              :: varid_myy      ! NetCDF variable  id for  day
  integer              :: varid_mgg      ! NetCDF variable  id for  hour
  integer              :: varid_ngg      ! NetCDF variable  id for  minute
  integer              :: varid_msec     ! NetCDF variable  id for  second

!==============================================================================
contains
!==============================================================================

  subroutine scan_obs_netcdf
  !------------------------------------------------------------------
  ! scan the NetCDF file and identify the content, required to decide
  ! on which PE the input processing of this file will take place
  !------------------------------------------------------------------

    !----------------
    ! local variables
    !----------------
    logical           :: lexicdf      ! exit of NetCDF observation input file
    integer           :: j            ! loop index
    integer           :: dim_id       ! dimension variable
    integer           :: entry        ! position in source file (subset)
    !----------------------------------------------------------
    ! variables for NetCDF file concerning unlimited dimension,
    ! number of dimensions, number of attributes
    !----------------------------------------------------------
    integer           :: status       ! NetCDF status variable
    integer           :: ncdims       ! NetCDF number of dimensions defined in this NetCDF file
    character(LEN=40) :: yname        ! NetCDF dimension name
    integer           :: ncvars       ! NetCDF number of variables  defined in this NetCDF file
    integer           :: ncatts       ! NetCDF number of attributes defined in this NetCDF file
    integer           :: unlim_dimid  ! NetCDF id of unlimited dimension defined in this NetCDF file
    !---------------
    ! output formats
    !---------------
    integer                   :: len_dim      ! length of dimension
    integer                   :: len_report   ! number of reports in NetCDF file

    integer                   :: if           ! file index
    integer                   :: ir           ! observation type
    integer                   :: nr           ! number of reports in file
    integer                   :: ns           ! number of subsets in file
    character(len=192)        :: filename     ! file name
    type(t_bufr_inv) ,pointer :: bi           ! pointer to table entry
    type(t_obsid)             :: obsid        ! observation id table entry
    integer                   :: obstype      ! observation type
    integer                   :: codetype     ! CMA codetype
    logical                   :: s2miss       ! Section 2 missing
    logical                   :: scan_p       ! scan input files in parallel?
    integer                   :: n_scan       ! # pes scanning input files
    integer                   :: js, je       ! aux. variables for strip mining
    integer                   :: is, ie       ! aux. variables for strip mining
    integer                   :: pe(n_source) ! processor scanning input file

!   scan_p = .true.
    scan_p = (netcdf_verb == 0)     ! Do not parallelize scan when debugging
    if (scan_p) then
       n_scan = max (dace% npe, 1)
    else
       n_scan = 1
    end if
    pe = dace% pio

    call flush_buf
    !------------------------
    ! loop over files to read (with loop strip-mining for robustness)
    !------------------------
    do js = 1, n_source, n_scan
       je = min (js + n_scan - 1, n_source)
       do if = js, je

          if (scan_p) pe(if) = mod (if - 1, n_scan)
          !-----------------------------------------------
          ! skip if not on this pe or if not a NetCDF-file
          !-----------------------------------------------
          if (pe(if) /= dace% pe .or. source(if)% filetype /= FT_NETCDF) cycle
          !------------------------------
          ! derive full pathname/filename
          !------------------------------
          filename = path_file (source(if)% path, source(if)% file)
          !-------------------------------
          ! check existence of NetCDF file
          !-------------------------------
          inquire (file=trim(filename), exist=lexicdf)
          if (.not. lexicdf) then
             print '("NOTE: MISSING FILE ",A)', filename
             source(if)% filetype = FT_MISSING
             cycle
          end if
          !--------------------------------
          ! open NetCDF file (if it exists)
          !--------------------------------
          status = nf90_open (filename, NF90_NOWRITE, ncid)
          if (status /= NF90_NOERR) THEN
             write(6,*) 'OPENING OF '// trim(filename) // ' FAILED'
             write(6,*) '   ON HOST '// trim(host)
             call finish('scan_obs_netcdf', &
                         'open failed for NetCDF file: ' // trim (filename))
          else
             call nextline
             call nextline
             write (oline(iol),'(i3,a,a)') if ,'. netCDF file ',trim(filename)
          endif
          !-------------------------------------------------------------------
          ! inquire about number of dimensions, variables, and attributes, and
          !         about umlimited dimension
          !-------------------------------------------------------------------

          !-----------------------------------------------------------------
          ! bufrx2netcdf now creates files with unlimited dimension(reports)
          !-----------------------------------------------------------------
          status = nf90_Inquire (ncid, ncdims, ncvars, ncatts, unlim_dimid)
          if (status /= NF90_NOERR) &
               call handle_err ('scan_obs_netcdf:Inquire unlim_dimid',status)

          !-----------------------------
          ! check some netCDF parameters
          !-----------------------------
          call nextline
          call nextline
          write (oline(iol),'(a,i9)') 'number of dimensions        :', ncdims
          call nextline
          write (oline(iol),'(a,i9)') 'number of variables         :', ncvars
          call nextline
          write (oline(iol),'(a,i9)') 'unlimited dimension         :', unlim_dimid
          call nextline

          if (ncdims > dimids_max) then
             write (6,'(a,i4,a,i3,a,i2)') 'filenr',if,                     &
                  'Maximum number of dimensions in netCDF file exceeded:', &
                  ncdims,' > ', dimids_max
             call finish('scan_obs_netcdf',                                    &
                         'Maximum number of dimensions in netCDF file exceeded')
          endif

          !-----------------------------------------------------------
          ! get information about dimension for reports in netCDF file
          !-----------------------------------------------------------
          do j = 1, ncdims
             status = nf90_Inquire_Dimension(ncid, j, yname,len_dim)
             IF (status /= NF90_NOERR) then
                print *, "nf90_Inquire_Dimension failed: ncid,j,yname,len_dim =",&
                     ncid,j,yname,len_dim
                call handle_err ('scan_obs_netcdf:Inquire_Dimension',status)
             end if
             !--------------------------------------------
             ! get number of observations from unlim_dimid
             !--------------------------------------------
             if (j == unlim_dimid) then
                !------------------
                ! number of reports
                !------------------
                len_report = len_dim
                nr         = len_report
                ns         = len_report
             endif
             if (netcdf_verb > 0) then
                call nextline
                write (oline(iol),'(a,i2,a,a16,a,i6)') &
                     'nf90_Inquire_Dimension(',j,') Dimension name(o): ',&
                     trim(yname),' length(o)=',len_dim
             end if

             status = nf90_Inq_Dimid (ncid, trim(yname), dim_id)
             if (status /= NF90_NOERR) then
                print *, "nf90_Inq_Dimid failed: ncid, yname =",ncid, trim (yname)
                call handle_err ('scan_obs_netcdf:Inq_Dimid', status)
             end if
             if (netcdf_verb > 0) then
                call nextline;
                write (oline(iol),'(a,i2,a,a16,a,i6)') &
                     'nf90_Inq_Dimid        (',j,') Dimension name(i): ',&
                     trim(yname),' dim_id(o)=',dim_id
             end if
          end do
          !------------------------------------------------------------
          ! handle 'empty' files from corrupted BUFR->NetCDF conversion
          !------------------------------------------------------------
          if (len_report == 0) then
             source(if)% filetype = FT_UNKNOWN
             source(if)% obstype  = -1
             source(if)% codetype = -1
             goto 999
          end if

          !--------------------------------------------------------
          ! common station identifier (integer, character) allocate
          !--------------------------------------------------------
          allocate (s1cat  (len_report))
          allocate (s1cats (len_report))
          allocate (s1catls(len_report))
          allocate (s1cent (len_report))
          allocate (s1cents(len_report))
          allocate (s2ikz  (len_report))
          !-----------------------------------------------
          ! get variables ID's of common station variables
          !-----------------------------------------------
          status = nf90_inq_varid (ncid, 'section1_date'                   , varid_s1date )
          status = nf90_inq_varid (ncid, 'section1_time'                   , varid_s1time )
          status = nf90_inq_varid (ncid, 'section1_data_category'          , varid_s1cat  )
          status = nf90_inq_varid (ncid, 'section1_int_data_sub_category'  , varid_s1cats )
          status = nf90_inq_varid (ncid, 'section1_local_data_sub_category', varid_s1catls)
          status = nf90_inq_varid (ncid, 'section1_centre'                 , varid_s1cent )
          status = nf90_inq_varid (ncid, 'section1_subcentre'              , varid_s1cents)
!         status = nf90_inq_varid (ncid, 'section1_update_sequence_nr'     , varid_s1updat)
          status = nf90_inq_varid (ncid, 'section2_ikz'                    , varid_s2ikz  )
          s2miss = (status /= NF90_NOERR)

          status = nf90_inq_varid (ncid, 'MLAH'  ,  varid_mlah)
!         special case: ship mapping on common varid
          if (status /= NF90_NOERR)                              &
               status = nf90_inq_varid (ncid, 'MLALA', varid_mlah)
          if (status /= NF90_NOERR) then
             call error ('no horizontal coordinate defined')
          endif
          status = nf90_inq_varid (ncid, 'MLOH'  ,  varid_mloh)
!         special case: ship mapping on common varid
          if  (status /= NF90_NOERR)                             &
               status = nf90_inq_varid (ncid, 'MLOLO', varid_mloh)
          if(status /= NF90_NOERR) then
             call error ('no horizontal coordinate defined')
          endif
          status = nf90_inq_varid (ncid, 'MJJJ', varid_mjjj)
          status = nf90_inq_varid (ncid, 'MMM' , varid_mmm )
          status = nf90_inq_varid (ncid, 'MYY' , varid_myy )
          status = nf90_inq_varid (ncid, 'MGG' , varid_mgg )
          status = nf90_inq_varid (ncid, 'NGG' , varid_ngg )
          status = nf90_inq_varid (ncid, 'MSEC', varid_msec)

          !-----------------------------------------------------
          ! get obstype from dbkz if section 2 is present
          ! ECMWF BUFR have a different section 2, don't use it!
          ! (exception: Aeolus lidar data processed by ECMWF)
          !-----------------------------------------------------
          status = nf90_get_var (ncid, varid_s1cent , s1cent )
          status = nf90_get_var (ncid, varid_s1cat  , s1cat  )
          if (all (s1cent==WMO0_ECMWF) .and. all (s1cat/=23)) s2miss = .true.
          if (.not. s2miss) then
             status  = nf90_get_var (ncid, varid_s2ikz, s2ikz)
             obsid   = obstype_dbkz ( s2ikz(1) )  ! obstype
             obstype = obsid% obstype
             ir      = obstype                    ! obstype
          else
             obstype = -1
             s2ikz   = -1
          endif

          !--------------------------------
          ! Handle missing or unknown  dbkz
          !--------------------------------
          if (obstype < 0) then
             status = nf90_get_var (ncid, varid_s1cats , s1cats )
             status = nf90_get_var (ncid, varid_s1catls, s1catls)
             status = nf90_get_var (ncid, varid_s1cents, s1cents)
             if (netcdf_verb > 0) then
                call nextline
                call nextline
                oline(iol) = "scan_obs_netcdf: Section 2 is missing or dbkz unknown"
             end if
             call guess_dbkz (s1cat, s1cats, s1catls, s1cent, s1cents, s2ikz, &
                              filename)
             if (netcdf_verb > 0) then
                call nextline
                write(oline(iol),'(A,i6)') "scan_obs_netcdf: guessing DBKZ:", s2ikz(1)
             end if
             obsid   = obstype_dbkz ( s2ikz(1) )  ! obstype
             obstype = obsid% obstype
             ir      = obstype                    ! obstype
          end if
          codetype = -1
          if (obstype >= 0) codetype = obsid% codetype

          !---------------------------------------------
          ! fill table entries
          ! derived type source has got one obstype only
          !---------------------------------------------
          source(if)% obstype  = ir
          source(if)% codetype = codetype
          source(if)% reports  = nr
          source(if)% subsets  = ns

          call nextline; write (oline(iol),'(a,i10)') 'file number       =', if
          call nextline; write (oline(iol),'(a,i10)') 'Datenbank Kennzahl=', s2ikz(1)
          call nextline; write (oline(iol),'(a,i10)') 'obstype           =', obstype
          call nextline; write (oline(iol),'(a,i10)') 'codetype          =', codetype

          if (ir < 1) then
             call nextline; write(oline(iol),*) trim(filename),' : invalid obstype =',ir
             source(if)% filetype = FT_UNKNOWN
          end if

999       continue
          !------------------
          ! close netCDF file
          !------------------
          status = nf90_close (ncid)
          if (status /= NF90_NOERR) &
               call handle_err ('scan_obs_netcdf:close', status)

       end do
       !------------------------------------------------------------
       ! Redistribute information collected so far to all processors
       !------------------------------------------------------------
       do if = js, je
          call p_bcast (source(if)% filetype, pe(if))
          if (source(if)% filetype /= FT_NETCDF) cycle
          call p_bcast (source(if)% obstype , pe(if))
          call p_bcast (source(if)% codetype, pe(if))
          call p_bcast (source(if)% reports , pe(if))
          call p_bcast (source(if)% subsets , pe(if))
       end do

       do if = js, je
          entry = 0
          if (source(if)% filetype /= FT_NETCDF) cycle
          ir = source(if)% obstype
          if (ir < 1) cycle
          nr = source(if)% reports
          ns = source(if)% subsets

          bi => bufr_inv(ir)
          bi% file   (if)  = .true.
          bi% nrec         = nr + bi% nrec
          bi% nsubset      = ns + bi% nsubset
          bi% subseto(if:) = ns + bi% subseto(if:)

          source(if)% entries = entry + ns

          if (pe(if) == dace% pe) then
             if (netcdf_verb > 0) then
                call nextline
                write (oline(iol),'(a,999a1)') 'file with bit information ', &
                     merge ('*','.',bi% file)
             end if

             call nextline; write (oline(iol),'(a,i10)') 'number of reports =', nr
             call nextline; write (oline(iol),'(a,i10)') 'number of subsets =', ns
             call nextline; write (oline(iol),'(a,i10)') 'reports (obstype) =', bi% nrec
             call nextline; write (oline(iol),'(a,i10)') 'subsets (obstype) =', bi% nsubset

             if (netcdf_verb > 0) then
                call nextline; oline(iol) = 'number of subsetso'
                do is = 1, n_source, 10
                   ie = min (is + 10 - 1, n_source)
                   call nextline; write (oline(iol),'(10i10)') bi% subseto(is:ie)
                end do
                call nextline; oline(iol) = 'number of obs in source'
                do is = 1, n_source, 10
                   ie = min (is + 10 - 1, n_source)
                   call nextline; write (oline(iol),'(10i10)') source(is:ie)% entries
                end do
             end if

             call nextline; oline(iol) = '______________________________________'
             call nextline
          end if
       end do
       call flush_buf
       call flush (6)

       !------------
       ! deallocate
       !------------
       if (allocated(s1cat  )) deallocate (s1cat  )
       if (allocated(s1cats )) deallocate (s1cats )
       if (allocated(s1catls)) deallocate (s1catls)
       if (allocated(s1cent )) deallocate (s1cent )
       if (allocated(s1cents)) deallocate (s1cents)
       if (allocated(s2ikz  )) deallocate (s2ikz  )

       !------------------------------------------------------
       ! flag NetCDF files to be read from multiple processors
       ! (currently implemented for surface, GPSGB only)
       !------------------------------------------------------
       do if = js, je
          if (source(if)% filetype /= FT_NETCDF) cycle

          select case (rept_char (source(if)% obstype)% obstype)
          case (OT_SYNOP, OT_DRIBU, OT_PAOB, OT_SCATT, OT_AIREP, OT_GPSGB)
             source(if)% complete = rept_use(source(if)% obstype)% read1pe
          case default
             source(if)% complete = .true.
          end select
       end do

    end do

  contains
    !--------------------------------------------------------------------------
    subroutine error (msg)
      character(len=*), intent(in) :: msg
      write (6,'(/,a,i0,": ",a)') 'pe=', dace% pe, &
                                  msg // ': ' // trim (filename)
      call flush (6)
      call finish ('scan_obs_netcdf', msg // ': ' // trim (filename))
    end subroutine error
    !--------------------------------------------------------------------------
  end subroutine scan_obs_netcdf
!------------------------------------------------------------------------------
  subroutine read_obs_netcdf (obs)
  !------------------------------------------------------------------------
  ! 1) read variables from BUFR/NetCDF file valid for all observation types
  ! 2) branch to observation type specific routines to read the remainder
  !------------------------------------------------------------------------

    type (t_obs) ,intent (inout) :: obs
    !------------------------------
    ! read observational data files
    !------------------------------
    type (t_head)      :: head         !
    type (t_use)       :: use          ! status variable
    integer            :: bufr_type    ! BUFR message type    read
    integer            :: bufr_subtype ! BUFR message subtype read
    integer            :: centre       ! generating centre
    logical            :: lkeep        ! true to keep current observation
    integer            :: nkeep        ! number of observations to keep
    integer            :: ifile        ! file counter
    integer            :: obstype      ! observation type
    integer            :: report_subt  ! observation subtype (Datenbankkennz.)
    integer            :: report_subti ! observation subtype index
    type (t_obsid)     :: obsid        ! observation id table entry
    integer            :: i_source     ! position in source file
    !-------------------------------------
    ! statistics on bufr records processed
    !-------------------------------------
    integer            :: iobst (n_ot) ! number of reports accepted
!   integer            :: jobst (n_ot) ! number of reports read
    integer            :: jother       ! number of other reports read
    integer            :: subset0      ! last  subset in previous file
    integer            :: subset1      ! last  subset in this     file
    integer            :: subsets      ! last  subset not to read
    integer            :: subsetl      ! last  subset     to read
    integer            :: rec1         ! first record     to read
    integer            :: recs         ! last  record not to read
    integer            :: recl         ! last  record     to read
    integer            :: recn         ! n.of  records    to read

    !----------------------------------------------
    ! determine content of observational data files
    !----------------------------------------------
    character(len=192) :: filename     ! file name

    !--------------------------------------------------------------------------------------------
    ! variable for     NetCDF file concerning unlimited dimension, Maximum number of dimensions,
    !                                         number of attributes
    !--------------------------------------------------------------------------------------------
    integer            :: status       ! NetCDF status variable
    integer            :: ncdims       ! NetCDF number of dimensions defined in this NetCDF file
    character (LEN=40) :: yname        ! NetCDF dimension name
    integer            :: ncvars       ! NetCDF number of variables  defined in this NetCDF file
    integer            :: ncatts       ! NetCDF number of attributes defined in this NetCDF file
    integer            :: unlim_dimid  ! NetCDF id for of unlimited dimension defined in this NetCDF file
    integer            :: dimid_ystid  !           for  character variable

    !----------------
    ! local variables
    !----------------
    integer            :: len_report   ! number of reports in NetCDF file
                                       ! NetCDF length of unlimited dimension defined in this NetCDF file
    integer            :: len_dim      ! length of dimension
    integer            :: j            ! loop index
    integer            :: ilen         ! actual length of netCDF-dimension

    integer            :: ireport      ! loop index
    logical            :: lcharstid    ! logical for expansion of character STID
    logical            :: lwigos       ! logical for WIGOS station ids to consider

    logical            :: lpr_head     ! logical for printing head of first report
    logical            :: l_istidb     ! Allocation status of istidb
    logical            :: l_istids     ! Allocation status of istids
    logical            :: s2miss       ! Section 2 missing or incomplete
    character(10)      :: statid       ! Station ID (reduced / short)
    character(32)      :: longid       ! Station ID / expanded WSI
    integer            :: n_miss       ! Count of stations with undefined ID
    integer            :: max_miss     ! Max. no. stations to warn about
    logical            :: ltrunc       ! Station name too long?
    logical            :: longname     ! Station name encoded in long name?
    logical,allocatable:: ltsi   (:)   ! WIGOS-encoded traditional identifier
    !-----------------------------------------
    ! Characters allowed in long station names
    !-----------------------------------------
    character(*), parameter :: digits      = "0123456789"
    character(*), parameter :: uppercase   = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    character(*), parameter :: valid_chars = digits // uppercase // "-_/"


    lpr_head = .FALSE.; if (netcdf_verb > 0) lpr_head = .TRUE.
!   lpr_head = .true.
    max_miss = 10;      if (netcdf_verb > 0) max_miss = huge (0)

    call flush_buf
    call flush (6)
    !-------------
    ! set counters
    !-------------
    iobst  = 0
!   jobst  = 0
    jother = 0

    if ( lpr_head ) then
      write (6,'(a,i3,a)') 'pe=',dace% pe,'  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write (6,'(a,i3,a)') 'pe=',dace% pe,'  BEGIN: subroutine read_obs_netcdf'
      write (6,'(a,i3,a)') 'pe=',dace% pe,'  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    end if

    !------------------------
    ! loop over files to read
    !------------------------
    i_source = 0

    do ifile = 1, n_source

      !--------------------------
      ! skip if not a NetCDF-file
      !--------------------------
      if (source(ifile)% filetype /= FT_NETCDF) cycle
      if (source(ifile)% obstype  == OT_SATEM)  cycle !+++ use specific routine ++

      !------------------------------------
      ! determine processor elements to use
      !       and range of reports to read
      !-----------------------------------
      if (source(ifile)% complete) then
        if (dace% pe /= source(ifile)%pe)       cycle
      else
        obstype = source(ifile)% obstype
        if (.not. source(ifile)% used)          cycle
        subset0 = 0
        if(ifile > 1) subset0 = bufr_inv(obstype)% subseto(ifile-1)
        subsetl = bufr_inv(obstype)% subsetl
        subsets = bufr_inv(obstype)% subsets
        subset1 = bufr_inv(obstype)% subseto(ifile)
        if (subsetl<=subset0)                   cycle
        if (subsets>=subset1)                   cycle
      endif

      !------------------------------
      ! derive full pathname/filename
      !------------------------------
      filename = path_file (source(ifile)% path, source(ifile)% file)
      ! -------------------------------
      ! open NetCDF file (if it exists)
      ! -------------------------------
      status = nf90_open (filename, NF90_NOWRITE, ncid)
      if (status /= NF90_NOERR) then
        write(6,*) 'OPENING OF '// trim(filename) // ' FAILED'
        write(6,*) '   ON HOST '// trim(host)
        call finish('read_obs_netcdf',                                &
                    'open failed for NetCDF file: ' // trim (filename))
      endif

      if (dace% lpio .or. source(ifile)% complete) then
        if (netcdf_verb == 0) then
           call nextline
           write (oline(iol),'(a,i4,i4,2a)') 'pe=',dace% pe,&
                ifile ,'. netCDF file ',trim(filename)
           call nextline
           write (oline(iol),'(a,i4,2a,i0)') 'pe=',dace% pe,' read_obs_netcdf:',&
                ' obstype = ',source(ifile)%obstype
        else
           write (6,'(a,i4,i4,2a)') 'pe=',dace% pe,&
                ifile ,'. netCDF file ',trim(filename)
           write (6,'(a,i4,2a,i0)') 'pe=',dace% pe,' read_obs_netcdf:',&
                ' obstype = ',source(ifile)%obstype
           call flush (6)
        end if
      end if
      obstype = source(ifile)%obstype

      call flush (6)
      i_source = i_source + 1

      !-------------------------------------------------------------------
      ! inquire about number of dimensions, variables, and attributes, and
      !         about umlimited dimension
      !-------------------------------------------------------------------
      status = nf90_Inquire (ncid, ncdims, ncvars, ncatts, unlim_dimid)

      len_report = -1
      do j = 1, ncdims
        status = nf90_Inquire_Dimension(ncid, j, yname,len_dim)
        if (status /= NF90_NOERR) then
           print *, "nf90_Inquire_Dimension failed: pe,ncid,j,yname,len_dim =", &
                 dace% pe,ncid, j, yname,len_dim
           call handle_err ('read_obs_netcdf:Inquire_Dimension',status)
        end if
        ! get number of observations from unlim_dimid
        if ( j == unlim_dimid ) then
          len_report   = len_dim
          exit
        endif
      enddo
      if (len_report < 0) call finish('read_obs_netcdf','no unlim_dimid')

      !-----------------------------
      ! set range of reports to read
      !-----------------------------
      if (.not.source(ifile)% complete) then
        recs = max(subset0,subsets)-subset0
        recl = min(subset1,subsetl)-subset0
        rec1 = recs + 1
        recn = recl - recs
      else
        recs = 0
        rec1 = 1
        recl = len_report
        recn = len_report
      endif

      !--------------------------------------------------------
      ! common station identifier (integer, character) allocate
      !--------------------------------------------------------
      allocate (edition  (recn))
      allocate (s1date   (recn))
      allocate (s1time   (recn))
      allocate (s1cat    (recn))
      allocate (s1cats   (recn))
      allocate (s1catls  (recn))
      allocate (s1cent   (recn))
      allocate (s1cents  (recn))
      allocate (s1updat  (recn))
      allocate (s2ikz    (recn))
      allocate (s2dcdate (recn))
      allocate (s2dctime (recn))

      allocate (stime    (recn))
      allocate (db_time  (recn))

      allocate (mlah     (recn))
      allocate (mloh     (recn))
      allocate (mjjj     (recn))
      allocate (mmm      (recn))
      allocate (myy      (recn))
      allocate (mgg      (recn))
      allocate (ngg      (recn))
      allocate (msec     (recn))
      allocate (msec_r   (recn))
      allocate (bdy_date (recn))
      allocate (bdy_time (recn))

      allocate (obs_time (recn))

      !----------------------------------------------------------------------
      ! derived station identifier (integer, character) allocate / initialise
      !----------------------------------------------------------------------
      allocate (istidn  (recn))
      allocate (ystidn  (recn))
      allocate (lvarid  (recn))
      allocate (lwsi    (recn))
      allocate (ltsi    (recn))
!     due to comparison in mo_synop: default is 0
!     consistent with ident default in mo_t_obs.f90
!     istidn (:) = -1
      istidn (:) =  0
      ystidn (:) = ' '
      lvarid (:) = .FALSE.
      n_miss     = 0
      longname   = .false.

      !-----------------------------------------------
      ! get variables ID's of common station variables
      !-----------------------------------------------
      status = nf90_inq_varid (ncid, 'edition_number'                   ,  varid_edition )
      status = nf90_inq_varid (ncid, 'section1_date'                    ,  varid_s1date  )
      status = nf90_inq_varid (ncid, 'section1_time'                    ,  varid_s1time  )
      status = nf90_inq_varid (ncid, 'section1_data_category'           ,  varid_s1cat   )
      status = nf90_inq_varid (ncid, 'section1_int_data_sub_category'   ,  varid_s1cats  )
      status = nf90_inq_varid (ncid, 'section1_local_data_sub_category' ,  varid_s1catls )
      status = nf90_inq_varid (ncid, 'section1_centre'                  ,  varid_s1cent  )
      status = nf90_inq_varid (ncid, 'section1_subcentre'               ,  varid_s1cents )
      status = nf90_inq_varid (ncid, 'section1_update_sequence_nr'      ,  varid_s1updat )

      status = nf90_inq_varid (ncid, 'section2_ikz'                     ,  varid_s2ikz   )
      s2miss = (status /= NF90_NOERR)
      status = nf90_inq_varid (ncid, 'section2_decoding_date'           ,  varid_s2dcdate)
      s2miss = (status /= NF90_NOERR) .or. s2miss
      status = nf90_inq_varid (ncid, 'section2_decoding_time'           ,  varid_s2dctime)
      s2miss = (status /= NF90_NOERR) .or. s2miss

      status = nf90_inq_varid (ncid, 'MLAH'  ,  varid_mlah)
!     special case: ship mapping on common varid
      if (status /= NF90_NOERR)                                         &
      status = nf90_inq_varid (ncid, 'MLALA',  varid_mlah )
      status = nf90_inq_varid (ncid, 'MLOH' ,  varid_mloh )
!     special case: ship mapping on common varid
      if  (status /= NF90_NOERR)                                        &
      status = nf90_inq_varid (ncid, 'MLOLO',  varid_mloh )
      status = nf90_inq_varid (ncid, 'MJJJ' ,  varid_mjjj )
      status = nf90_inq_varid (ncid, 'MMM'  ,  varid_mmm  )
      status = nf90_inq_varid (ncid, 'MYY'  ,  varid_myy  )
      status = nf90_inq_varid (ncid, 'MGG'  ,  varid_mgg  )
      status = nf90_inq_varid (ncid, 'NGG'  ,  varid_ngg  )
      status = nf90_inq_varid (ncid, 'MSEC' ,  varid_msec )

      !---------------------------
      ! get section1/2 information
      !---------------------------
      status = nf90_get_var (ncid, varid_edition , edition  ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_s1date  , s1date   ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_s1time  , s1time   ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_s1cat   , s1cat    ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_s1cats  , s1cats   ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_s1catls , s1catls  ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_s1cent  , s1cent   ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_s1cents , s1cents  ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_s1updat , s1updat  ,start=(/rec1/))

      status = nf90_get_att (ncid, varid_s1date  , '_FillValue', imissing)

      !-----------------------------------------------------
      ! try to get obstype from first dbkz
      ! ECMWF BUFR have a different section 2, don't use it!
      ! (exception: Aeolus lidar data processed by ECMWF)
      !-----------------------------------------------------
      obstype = -1
      if (all(s1cent==WMO0_ECMWF) .and.  all(s1cat/=23)) s2miss = .true.
!     if (all(s1cent==WMO0_ECMWF)) s2miss = .true.
      if (.not. s2miss) then
        status  = nf90_get_var (ncid, varid_s2ikz   , s2ikz    ,start=(/rec1/))
        status  = nf90_get_var (ncid, varid_s2dcdate, s2dcdate ,start=(/rec1/))
        status  = nf90_get_var (ncid, varid_s2dctime, s2dctime ,start=(/rec1/))
        obsid   = obstype_dbkz ( s2ikz(1) )  ! obstype
        obstype = obsid% obstype
        !-------------------------------------------
        ! Handle section 2 decoding date: 00.00.2000
        !-------------------------------------------
        where (s2dcdate == 20000000)
           s2dcdate = 0 ! treat as missing
           s2dctime = 0
        end where
      else
        s2dcdate = 0 ! was: s1date, but should be treated as missing
        s2dctime = 0 ! was: s1time
      endif

      if (obstype < 0) then
        if (netcdf_verb > 0) then
          write (*,'(/,A)') "read_obs_netcdf: Section 2 is missing"
        end if
        call guess_dbkz (s1cat, s1cats, s1catls, s1cent, s1cents, s2ikz, &
                         filename)
        if (netcdf_verb > 0) then
          write(*,'(A,i6)') "read_obs_netcdf: guessing DBKZ:", s2ikz(1)
          write(*,'(A)') "read_obs_netcdf: setting decoding time to synop time"
        end if
      end if

      !-----------------
      ! read coordinates
      !-----------------
      status = nf90_get_var (ncid, varid_mlah, mlah ,start=(/rec1/))
      status = nf90_get_att (ncid, varid_mlah , '_FillValue', rmissing)
      status = nf90_get_var (ncid, varid_mloh, mloh ,start=(/rec1/))

      !------------------
      ! read date /  time
      !------------------
      status = nf90_get_var (ncid, varid_mjjj, mjjj ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_mmm , mmm  ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_myy , myy  ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_mgg , mgg  ,start=(/rec1/))
      status = nf90_get_var (ncid, varid_ngg , ngg  ,start=(/rec1/))
      !-----------------------------------
      ! check for presence of 'MSEC' entry
      !-----------------------------------
      msec   = 0
      status = nf90_inq_varid (ncid, 'MSEC' ,  varid_msec )
      if(status == NF90_NOERR) then
        if (s2ikz(1) == 1694  .or. s2ikz(1) == 1695) then
          !-------------------------------------------------
          ! float for radiooccultation, transform to integer
          !-------------------------------------------------
          status = nf90_get_var (ncid, varid_msec, msec_r ,start=(/rec1/))
          if (status == NF90_NOERR) then
            where(msec_r < 0. .or. msec_r > 60.) msec_r = 0.
            msec = nint (msec_r)
          else
            msec = 0
          endif
        else
          !---------------------------
          ! integer for other obstypes
          !---------------------------
          status = nf90_get_var (ncid, varid_msec, msec ,start=(/rec1/))
        endif
      endif

      !---------------------------------------------------
      ! derive time, data bank time(array)
      ! fix for old data +++ this is the year 2030 bug +++
      !---------------------------------------------------
      where (s1date > 20300000 .and. edition < 4) s1date = s1date - 1000000
      stime   = time_yyyymmdd_hhmmss ( s1date  , s1time  )

      where (s2dcdate /= 0)
        db_time = time_yyyymmdd_hhmmss (s2dcdate, s2dctime)
      elsewhere
        db_time = zero_time
      end where

      !------------------------------
      ! create body time date (array)
      !------------------------------
      bdy_date (:)  = -1
      bdy_time (:)  = -1
!     observation date: CCYYMMDD
      bdy_date  = mjjj * 10000 + mmm * 100 + myy

      msec = min( max ( msec, 0), 60)
!     observation time: HHMMSS   ; sec
      bdy_time  = mgg  * 10000 + ngg * 100 + msec

      !------------------------------
      ! create body time date (array)
      !------------------------------
      obs_time = time_yyyymmdd_hhmmss ( bdy_date, bdy_time)

      !-----------------------------------
      ! check for WIGOS station identifier
      !-----------------------------------
      call check_wigos (lwigos)

      !---------------------------------------------------------------
      ! Preferences / precedences for the actually used station ID:
      ! 1) numeric station ID (WMO block/station; buoy; satellite ID)
      ! 2) traditional station ID (TSI), encoded as short station name
      ! 3) WIGOS station identifier
      ! 4) long station name (as fallback and after validation only)
      !---------------------------------------------------------------

      !----------------------------------------------------
      ! get integer WMO block and station number if present
      !----------------------------------------------------
                                status = nf90_inq_varid (ncid, 'MII'  ,  varid_mii)
      if (status == NF90_NOERR) status = nf90_inq_varid (ncid, 'NIII' ,  varid_niii)
      if (status == NF90_NOERR) then
        !------------------
        ! allocate arrays
        !------------------
        allocate ( mii  (recn) )
        allocate ( niii (recn) )
        !------------------
        ! get values
        !------------------
        status = nf90_get_var (ncid, varid_mii  , mii  ,start=(/rec1/))
        status = nf90_get_var (ncid, varid_niii , niii ,start=(/rec1/))

        do ireport = 1 , recn
          if ((mii(ireport) /= imissing) .and. (niii(ireport) /= imissing)) then
            istidn(ireport) =  mii(ireport) * 1000 + niii(ireport)
            write( ystidn(ireport),'(I5.5)' )  istidn(ireport)
            if (istidn(ireport) > 999   .and. &
                istidn(ireport) < 99999       ) lvarid(ireport) = .true.
          endif
        enddo
        if (allocated(mii )) deallocate (mii)
        if (allocated(niii)) deallocate (niii)
      !---------------------------------------------------------
      ! else get integer WMO buoy/platform identifier if present
      !---------------------------------------------------------
      elseif (nf90_inq_varid (ncid,'MABNN',varid_mabnn ) == NF90_NOERR) then
        allocate ( istidb (recn) )
        status = nf90_get_var (ncid, varid_mabnn, istidb ,start=(/rec1/))
        do ireport = 1 , recn
          if (istidb(ireport) /= imissing) then
            select case (istidb(ireport))
            case (     0:  99999)
              write (ystidn(ireport),'(I5.5)') istidb(ireport)
              lvarid(ireport) = .TRUE.
            case (100000:9999999)
              write (ystidn(ireport),'(I7.7)') istidb(ireport)
              lvarid(ireport) = .TRUE.
            end select
          endif
        enddo
      !-----------------------------------------------------
      ! else get integer WMO satellite identifier if present
      ! (descriptor 0 01 007 has 10 bits)
      !-----------------------------------------------------
      elseif (nf90_inq_varid (ncid,'MI1I2',varid_mi1i2 ) == NF90_NOERR) then
        allocate ( istids (recn) )
        status = nf90_get_var (ncid, varid_mi1i2, istids ,start=(/rec1/))
        do ireport = 1 , recn
          if ((istids(ireport) /= imissing) .and. (istids(ireport) < 1023)) then
            write (ystidn(ireport),'(I4.4)') istids(ireport)
            lvarid(ireport) = .TRUE.
          endif
        enddo
      endif

      !---------------------
      ! combine integer stid
      !---------------------
      l_istidb = allocated(istidb)
      l_istids = allocated(istids)

#if 0   /* original code version */
      do ireport = 1 , recn
        if(lvarid(ireport)) then
          if (istidn (ireport) /=  0      ) then
            continue
          else if ( l_istidb ) then
            if ((istidb(ireport) /= imissing) .and. (istidb(ireport) < 9999999)) THEN
              istidn (ireport) = istidb(ireport)
            endif
          else if ( l_istids ) then
            if ((istids(ireport) /= imissing) .and. (istids(ireport) < 999)) THEN
              istidn (ireport) = istids(ireport)
            endif
          endif
        endif
      enddo
#else   /* better vectorizing version */
      if      (l_istidb) then
        do ireport = 1, recn
          if (istidn(ireport) == 0 .and. lvarid(ireport)) then
            if ((istidb(ireport) /= imissing) .and. (istidb(ireport) < 9999999)) THEN
              istidn(ireport) = istidb(ireport)
            endif
          end if
        end do
      else if (l_istids) then
        do ireport = 1, recn
          if (istidn(ireport) == 0 .and. lvarid(ireport)) then
            if ((istids(ireport) /= imissing) .and. (istids(ireport) < 1023)) THEN
              istidn(ireport) = istids(ireport)
            endif
          end if
        end do
      end if
#endif

      if (allocated(istidb )) deallocate (istidb)
      if (allocated(istids )) deallocate (istids)

      !-------------------------------------------------------
      ! check whether all stid are defined by integer section
      !-------------------------------------------------------
      lcharstid = .NOT. ALL (lvarid(1:recn))

      if (lcharstid) then
        !---------------------------
        ! get character station ID's
        !---------------------------

        ymissing    = NF90_FILL_CHAR
        dimid_ystid = -999

        IF     (nf90_inq_varid (ncid,'YDDDD' ,varid_ystid ) == NF90_NOERR) THEN
          status = nf90_inq_dimid (ncid,'YDDDD_strlen' ,dimid_ystid)
          status = nf90_get_att   (ncid, varid_ystid, '_FillValue', ymissing)
          if (status /= NF90_NOERR) ymissing = NF90_FILL_CHAR

        ELSEIF (nf90_inq_varid (ncid,'YSSOSN',varid_ystid ) == NF90_NOERR) THEN
          status = nf90_inq_dimid (ncid,'YSSOSN_strlen',dimid_ystid)
          status = nf90_get_att   (ncid, varid_ystid, '_FillValue', ymissing)
          if (status /= NF90_NOERR) ymissing = NF90_FILL_CHAR

        ELSEIF (nf90_inq_varid (ncid,'YAIRN' ,varid_ystid ) == NF90_NOERR) THEN
          status = nf90_inq_dimid (ncid,'YAIRN_strlen' ,dimid_ystid)
            status = nf90_get_att   (ncid, varid_ystid, '_FillValue', ymissing)
          if (status /= NF90_NOERR) ymissing = NF90_FILL_CHAR

        ELSEIF (nf90_inq_varid (ncid,'YCCC8' ,varid_ystid ) == NF90_NOERR) THEN
          status = nf90_inq_dimid (ncid,'YCCC8_strlen' ,dimid_ystid)
          status = nf90_get_att   (ncid, varid_ystid, '_FillValue', ymissing)
          if (status /= NF90_NOERR) ymissing = NF90_FILL_CHAR

        ELSEIF (nf90_inq_varid (ncid,'YSOSN' ,varid_ystid ) == NF90_NOERR) THEN
          status = nf90_inq_dimid (ncid,'YSOSN_strlen' ,dimid_ystid)
          status = nf90_get_att   (ncid, varid_ystid, '_FillValue', ymissing)
          if (status /= NF90_NOERR) ymissing = NF90_FILL_CHAR

        ELSEIF (nf90_inq_varid (ncid,'YXXNN' ,varid_ystid ) == NF90_NOERR) THEN
          status = nf90_inq_dimid (ncid,'YXXNN_strlen' ,dimid_ystid)
          status = nf90_get_att   (ncid, varid_ystid, '_FillValue', ymissing)
          if (status /= NF90_NOERR) ymissing = NF90_FILL_CHAR

        ELSEIF (nf90_inq_varid (ncid,'YSBPI' ,varid_ystid ) == NF90_NOERR) THEN
          status = nf90_inq_dimid (ncid,'YSBPI_strlen' ,dimid_ystid)
          status = nf90_get_att   (ncid, varid_ystid, '_FillValue', ymissing)
          if (status /= NF90_NOERR) ymissing = NF90_FILL_CHAR

        ELSEIF (nf90_inq_varid (ncid,'YATNO' ,varid_ystid ) == NF90_NOERR) THEN
          status = nf90_inq_dimid (ncid,'YATNO_strlen' ,dimid_ystid)
          status = nf90_get_att   (ncid, varid_ystid, '_FillValue', ymissing)
          if (status /= NF90_NOERR) ymissing = NF90_FILL_CHAR

        ELSEIF (nf90_inq_varid (ncid,'YLSNA' ,varid_ystid ) == NF90_NOERR) THEN
          status = nf90_inq_dimid (ncid,'YLSNA_strlen' ,dimid_ystid)
          status = nf90_get_att   (ncid, varid_ystid, '_FillValue', ymissing)
          if (status /= NF90_NOERR) ymissing = NF90_FILL_CHAR
          longname   = .true.

        ENDIF
        iascii_ymissing = iachar(ymissing)

        !----------------------------------------------------------------------
        ! WIGOS (local) identifier has higher precedence than long station name
        !----------------------------------------------------------------------
        if (lwigos .and. (dimid_ystid < 0 .or. longname)) then
           ilen = 16
           do ireport = 1, recn
              if (ltsi(ireport) .and. .not. lvarid(ireport)) then
                 ystidn(ireport) = wsi(ireport)% wigli
                 lvarid(ireport) = .true.
              end if
           end do
        end if

        if (dimid_ystid >= 1) then
          !-------------------------------------------
          !   get NetCDF string length of station ID's
          !-------------------------------------------
          status = nf90_Inquire_Dimension (ncid, dimid_ystid, len=ilen)
          if (status == NF90_NOERR) then
            if ( ilen < 1 ) then
              call finish('read_obs_netcdf',                       &
                          'string length of station ID ill defined')
            else if ( ilen > ilstidn ) then
              call nextline
              write (oline(iol),'(a,i4,a,2i6," :")') 'pe=', dace% pe,   &
                   ' WMO station id length > length of string ilstidn:',&
                   ilen, ilstidn
              call nextline
              oline(iol) = filename
              call nextline
            endif
          endif

          !-----------------------------
          !   get character station ID's
          !-----------------------------
          allocate (ystid (ilen, recn))
          ystid = ' '
          status = nf90_get_var (ncid, varid_ystid, ystid , &
                                       start=(/1   ,rec1/), &
                                       count=(/ilen,recn/))
          if (status /= NF90_NOERR) &
               CALL handle_err ('read_obs_netcdf:get_var ystid',status)
        end if

        if (dimid_ystid >= 1 .or. lwigos) then
          !--------------------------------
          ! copy to internal representation
          !--------------------------------
          do ireport = 1 , recn
            if (lvarid(ireport))                               cycle ! Keep WMO id
            if (lwsi(ireport)) then
              call wsi_encode (wsi(ireport), statid, wsihash(ireport))
              if (wsi_mode /= 1 .or. ystid(1,ireport) /= wsi_prefix) then
                ystidn(ireport) = statid
                lvarid(ireport) = .true.
                cycle
              end if
            end if
            if (iachar(ystid (1,ireport)) == iascii_ymissing)  cycle
            if (wsi_from_statid .and. ystid(1,ireport) == "0" &
                                .and. ystid(2,ireport) == "-" ) then
               longid = ""
               do j = 1, min (ilen, len (longid))
                  if (ystid(j,ireport) == ymissing) exit
                  longid(j:j) = ystid(j,ireport)
               end do
               call wsi_decode (longid, wsi(ireport), status)
               if (status == 0) then
                  !-----------------------------------------------------------
                  ! Treat WSI consistent with a WMO program as traditional ID.
                  ! This requires 18 characters for non-truncated encoding.
                  !-----------------------------------------------------------
                  if (ilen                >= 18    .and. &
                      wsi(ireport)% wigii >= 20000 .and. &
                      wsi(ireport)% wigii <  24000 .and. &
                      wsi(ireport)% wigin == 0           ) then
                     ystidn(ireport) = wsi(ireport)% wigli
                     lvarid(ireport) = .true.
                     cycle
                  end if
                  select case (wsi_mode)
                  case (1,2)
                     j = 999
                  case default !(3)
                     j = 19999
                  end select
                  if (wsi(ireport)% wigii <= j) then
                     !-------------------------------
                     ! Re-encode WSI to internal form
                     !-------------------------------
                     call wsi_encode (wsi(ireport), statid, wsihash(ireport))
                     ystidn(ireport) = statid
                     lvarid(ireport) = .true.
                     lwsi(ireport)   = .true.
                     cycle
                  else
                     wsi(ireport)% valid    = .false.
                  end if
               end if
            end if
            !-----------------------------------------------------
            ! Check for truncation of station names.  Take care of
            ! padding with "missing characters" instead of blanks.
            !-----------------------------------------------------
            ltrunc = ilen > ilstidn
            if (ltrunc) then
               ltrunc = iachar (ystid (ilstidn+1,ireport)) /= iascii_ymissing
            end if
            if (ltrunc .and. any (ystid (ilstidn+1:,ireport) /= '')) then
              !------------------------------------------
              ! skip station ids which would be truncated
              !------------------------------------------
              call nextline
              write (oline(iol),*)                                             &
                   'STATION NAME WOULD BE TRUNCATED, SKIPPED: ',ystid(:,ireport)
              cycle
            endif
            do j     = 1 , MIN( ilen, ilstidn )
              if (iachar(ystid (j,ireport)) == iascii_ymissing) then
                !--------------------------------------------
                ! in order to get all characters initialised,
                ! replace character field empty by blank
                !--------------------------------------------
                ystidn (ireport) (j:)  =  ' '
                exit
              endif
              ystidn (ireport) (j:j)  =  ystid(j,ireport)
              lvarid(ireport) = .TRUE.
            enddo
            if (longname) then
               if (verify (trim (ystidn(ireport)), valid_chars) > 0) then
                  call nextline
                  write (oline(iol),*)                                        &
                       'Station name INVALID, SKIPPED: ',trim (ystidn(ireport))
                  lvarid(ireport) = .false.
                  ystidn(ireport) = ""
               cycle
               end if
            end if
          enddo
          !----------------------------
          ! print first stid if defined
          !----------------------------
          if (lpr_head) then
          write (6,'(a,i3,a,/,a,i2,/,a,a)')              &
               'pe=', dace% pe,' character station ID:', &
               '       DB name length =', ilen,          &
               '       mod station ID = ', ystidn(1)
          if (ystid(1,1) == ymissing) then
            write (6,'(a,i3,a  )') 'pe=',dace% pe,' ori station ID is undefined'
          else
            where (ystid(1:ilen,1) == ymissing) ystid(1:ilen,1) = " "
            write (6,'(a,i3,99a)') 'pe=',dace% pe,' ori station ID = ', ystid(1:ilen,1)
          end if
          call flush (6)
          end if

          if (allocated(ystid   )) deallocate (ystid   )
        endif

        !------------------------------------
        ! check whether all stid are defined
        !------------------------------------
        lcharstid = .NOT. ALL (lvarid(1:recn))
        if (lcharstid) then
          do ireport = 1 , recn
            if (.NOT.lvarid(ireport)) then
              if (n_miss < max_miss) then
                if (n_miss == 0) call nextline
                call nextline
                write (oline(iol),'(a,i4,a,i5,i10,a)') 'pe=',dace% pe, &
                     ' read_obs_netcdf', ifile, ireport,               &
                     ' ******** station ID is undefined ********'
              end if
              n_miss = n_miss + 1
            endif
          enddo
          if (n_miss > max_miss) then
            call nextline
            write (oline(iol),'(a,i4,a,i5," :",i8,a)') 'pe=',dace% pe, &
                 ' read_obs_netcdf', ifile, n_miss - max_miss,         &
                 ' further  station ID warning(s) suppressed'
          end if
          call flush (6)
        end if
      endif ! end of get character station ID block

      !------------------------------------------------------
      ! loop over netCDF observations(station) for statistics
      !------------------------------------------------------
      do ireport = 1 ,  min(1 ,recn)

        !=======================================
        ! derive observation type specifications
        !=======================================
        report_subt  = s2ikz(ireport)
        bufr_type    = s1cat(ireport)
        bufr_subtype = s1catls(ireport)
        centre       = s1cent(ireport)
        if (report_subt >= 0) then
          !---------------------------------
          ! derive information from DWD dbkz
          !---------------------------------
          obsid      = obstype_dbkz (report_subt)
          if ( lpr_head ) then
            write (6,'(/,a,i3,a,i10,a,/,a,/,i5,i8,i9,i8,i7,i9,i6,2x,a,/)')             &
                     'pe=',dace% pe,' read_obs_netcdf ',ireport ,'. report from DWD ' ,&
                     'DWD dbkz WMO typ Loc subc centre obstype codetype used TYP LAND',&
                      obsid
!                   ----5-------8--------9--------8-----6---------9-----5--
          endif
        else
          !---------------------------------------------------------------
          ! or from bufr_type, bufr_subtype specified by generating center
          !---------------------------------------------------------------
          obsid      = obstype_bufr (bufr_type, bufr_subtype, centre)
          !------------------------
          ! optionally set DWD dbkz
          !------------------------
          if (derive_dbkz) report_subt   = obsid% dbkz
        endif
        !-------------------------------------------------------
        ! set CMA obstype, BUFR type, BUFR subtype, if not given
        !-------------------------------------------------------
        obstype                          = obsid% obstype
        if (bufr_type   <0) bufr_type    = obsid% bufrtype
        if (bufr_subtype<1 .and. obsid% centre == WMO0_ECMWF) &
                            bufr_subtype = obsid% subtype
        if (obstype < 0) cycle
        report_subti  = idb_dbk (report_subt, obstype)

        !--------------------------
        ! control of processed data
        !--------------------------
        head% obstype     = obstype
        head% dbkz        = report_subt
        head% modtype     = rept_char(obstype)% mod
        head% buf_type    = bufr_type
        head% buf_subtype = bufr_subtype
        head% codetype    = obsid% codetype
        head% idbk        = report_subti
        head% source      = ifile
        head% record      = i_source
        head% center      = s1cent(ireport)
        head% subcenter   = s1cents(ireport)
!       nsubset           = 1

        if ( lpr_head ) then
          write (6,'(a,i3,a    ,/ ,                          &
           &        6(a19, i6  ,/),                          &
           &        2(a19, a   ,/),                          &
           &        6(a19, i6  ,/) )' )                      &
           'pe=',dace% pe,' read_obs_netcdf ',               &
           ' head% obstype     ',obstype                    ,&
           ' head% dbkz        ',report_subt                ,&
           ' head% modtype     ',rept_char(obstype)% mod    ,&
           ' head% buf_type    ',bufr_type                  ,&
           ' head% buf_subtype ',bufr_subtype               ,&
           ' head% codetype    ',obsid% codetype            ,&
           ' head% time        ',cyyyymmddhhmmss (stime  (ireport:ireport)),&
           ' head% db_time     ',cyyyymmddhhmmss (db_time(ireport:ireport)),&
           ' head% idbk        ',report_subti               ,&
           ' head% source      ',ifile                      ,&
           ' head% record      ',i_source                   ,&
           ' head% center      ',s1cent(ireport)            ,&
           ' head% subcenter   ',s1cents(ireport)           ,&
           ' nsubset           ',1
        endif

        !----------------------------------------------
        ! still some observation type specific counting
        ! (this is currently dead code)
        !----------------------------------------------
!       select case (obstype)
!       case (OT_GPSRO, OT_SYNOP, OT_DRIBU, OT_PAOB, OT_SCATT, OT_SOIL)
!         jobst(obstype) = jobst(obstype) + 1
!       case (OT_SATOB)
!         select case (s2ikz(ireport))
!         case (DK_SATOB_2, DK_SATOB_3, &
!               DK_AMV_EUMETSAT, DK_AMV_GOES, DK_AMV_MODIS, &
!               DK_AMV_FY_X, DK_AMV_MTV,                    &
!               DK_AMV_NOAA)
!           jobst(obstype) = jobst(obstype) + 1
!         case default
!           jother = jother + 1
!         end select
!       case (OT_AIREP)
!         jobst(obstype) = jobst(obstype) + 1
!       case (OT_TEMP, OT_PILOT)
!         select case (s2ikz(ireport))
!         case (508, 509, 510, 511, &! PILOT      A,B,C,D geopotentielle Hoehe
!               512, 513, 514, 515, &! PILOT      A,B,C,D
!               516, 517, 518, 519, &! TEMP       A,B,C,D Mobil
!               520, 521, 522, 523, &! TEMP       A,B,C,D
!               524, 525, 526, 527, &! TEMP       A,B,C,D Mobil
!               764, 765, 766, 767, &! PILOT SHIP A,B,C,D geopotentielle Hoehe
!               768, 769, 770, 771, &! PILOT SHIP A,B,C,D
!               776, 777, 778, 779, &! TEMP  SHIP A,B,C,D
!               780, 781, 782, 783, &! TEMP  DROP A,B,C,D
!               10520, 10521      , &! TEMP  BUFR
!               10526, 10527, 10574,&! TEMP  BUFR high resolution
!               10776, 10777      , &! TEMP  SHIP BUFR
!               10782, 10783, 10785,&! TEMP  SHIP BUFR high resolution
!               10516, 10517      , &! TEMP Mobil BUFR
!               10570             , &! TEMP Mobil BUFR descent
!               10780             , &! TEMP  DROP BUFR
!               10553             , &! WINDPROFILER
!               548               , &! WINDPROFILER (u,v,(w))
!               10600)               ! VAD WIND PROFILE
!           jobst(obstype) = jobst(obstype) + 1
!         case default
!           jother = jother + 1
!         end select
!       end select
      !----------------------------------------------
      ! end of loop over netCDF observations(station)
      !----------------------------------------------
      enddo

      !--------------------------------------------
      ! call observation type specific read-routine
      !--------------------------------------------
      select case (obstype)
      case (OT_GPSRO)
        !------------------------
        ! GNSS Radio occultations
        !------------------------
        call read_gpsro_netcdf (ifile,i_source, obs, lkeep, nkeep)
        if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
      case (OT_SYNOP, OT_DRIBU, OT_PAOB)
        !------------------------------
        ! SYNOP surface data - land,sea
        !------------------------------
        call read_synop_netcdf (ifile,i_source, obs, rec1, recl, lkeep, nkeep)
        if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
      case (OT_SCATT)
        !-------------------
        ! Scatterometer data
        !-------------------
        if (rept_use(OT_SCATT) % use(CHK_NONE) > STAT_FORGET) then
          call read_scatt_netcdf (ifile,i_source, obs, rec1, recl, lkeep, nkeep)
          if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
        end if
        !---------------------------------
        ! ASCAT soil moisture observations
        ! (in same reports as wind data)
        !---------------------------------
        if (rept_use(OT_SOIL) % use(CHK_NONE) > STAT_FORGET) then
          call read_soil_netcdf (ifile, i_source, obs, head, lkeep, nkeep)
          if (lkeep) iobst(OT_SOIL) = iobst(OT_SOIL) + nkeep
        end if
      case (OT_SATOB)
        !----------------------------
        ! satellite wind observations
        !----------------------------
        call read_amv_netcdf (ifile,i_source, obs, lkeep, nkeep)
        if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
      case (OT_TEMP, OT_PILOT)
        !------------------------------------------
        ! vertical soundings (other than satellite)
        !------------------------------------------
        call read_temp_netcdf (ifile,i_source, obs, lkeep, nkeep)
        if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
      case (OT_AIREP)
        !-----------------
        ! Aircraft reports
        !-----------------
        call read_airep_netcdf (ifile,i_source, obs, rec1, recl, lkeep, nkeep)
        if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
      case (OT_GPSGB)
        !-------------------
        ! STD or ZTD reports
        !-------------------
        call read_std_netcdf ( ifile, i_source, obs, head, rec1, recl, lkeep, nkeep)
        if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
      case (OT_WLIDAR)
        !--------------------------
        ! HLOS satellite wind lidar
        !--------------------------
        call read_wlidar_netcdf (ifile,i_source, obs, head, lkeep, nkeep)
        if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
      case (OT_SOIL)
        !---------------------------------
        ! ASCAT soil moisture observations
        !---------------------------------
        call read_soil_netcdf (ifile,i_source, obs, head, lkeep, nkeep)
        if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
      case default
        call dism_report_0 (use, head, CHK_NOIMPL) !Buchfuehrung stimmt nicht
        jother = jother + 1                        !++++ lost < 0 ++++
      end select

      !------------------
      ! close netCDF file
      !------------------
      status = nf90_close (ncid)
      if (status /= NF90_NOERR) &
           call handle_err ('read_obs_netcdf:close', status)
      if (dace% lpio .or. source(ifile)% complete) then
        if (netcdf_verb == 0) then
           call nextline; oline(iol) = '____________________________________'
           call nextline
        else
           write (6,'(a,/)')           '____________________________________'
           call flush (6)
        end if
      end if


      !------------
      ! deallocate
      !------------
      if (allocated(edition )) deallocate (edition )
      if (allocated(s1date  )) deallocate (s1date  )
      if (allocated(s1time  )) deallocate (s1time  )
      if (allocated(s1cat   )) deallocate (s1cat   )
      if (allocated(s1cats  )) deallocate (s1cats  )
      if (allocated(s1catls )) deallocate (s1catls )
      if (allocated(s1cent  )) deallocate (s1cent  )
      if (allocated(s1cents )) deallocate (s1cents )
      if (allocated(s1updat )) deallocate (s1updat )
      if (allocated(s2ikz   )) deallocate (s2ikz   )
      if (allocated(s2dcdate)) deallocate (s2dcdate)
      if (allocated(s2dctime)) deallocate (s2dctime)
      if (allocated(mlah    )) deallocate (mlah    )
      if (allocated(mloh    )) deallocate (mloh    )
      if (allocated(mjjj    )) deallocate (mjjj    )
      if (allocated(mmm     )) deallocate (mmm     )
      if (allocated(myy     )) deallocate (myy     )
      if (allocated(mgg     )) deallocate (mgg     )
      if (allocated(ngg     )) deallocate (ngg     )
      if (allocated(msec    )) deallocate (msec    )
      if (allocated(msec_r  )) deallocate (msec_r  )
      if (allocated(bdy_date)) deallocate (bdy_date)
      if (allocated(bdy_time)) deallocate (bdy_time)

      if (allocated(stime   )) deallocate (stime   )
      if (allocated(db_time )) deallocate (db_time )
      if (allocated(obs_time)) deallocate (obs_time)

      if (allocated(istidn  )) deallocate (istidn  )
      if (allocated(ystidn  )) deallocate (ystidn  )
      if (allocated(lvarid  )) deallocate (lvarid  )
      if (allocated(istidb  )) deallocate (istidb  )
      if (allocated(istids  )) deallocate (istids  )
      if (allocated(ystid   )) deallocate (ystid   )

      if (allocated(ltsi    )) deallocate (ltsi    )
      if (allocated(lwsi    )) deallocate (lwsi    )
      if (allocated(wsi     )) deallocate (wsi     )
      if (allocated(wsihash )) deallocate (wsihash )
    !------------------------------
    ! end of loop over input files
    !------------------------------
    enddo

    call flush_buf
    if ( lpr_head ) then
      write (6,'(a,i3,a)')  'pe=',dace% pe,'  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write (6,'(a,i3,a)')  'pe=',dace% pe,'  ENDE : subroutine read_obs_netcdf'
      write (6,'(a,i3,a)')  'pe=',dace% pe,'  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
    endif

  contains
    !--------------------------------------------------------------------------
    subroutine check_wigos (lwigos)
      logical, intent(out) :: lwigos
      !------------------------------
      ! Check for of WIGOS station id
      !------------------------------
      integer      ,allocatable :: wigis(:)     ! WIGOS identifier series
      integer      ,allocatable :: wigii(:)     ! WIGOS issuer of identifier
      integer      ,allocatable :: wigin(:)     ! WIGOS issue number
      character(16),allocatable :: wigli(:)     ! WIGOS local identifier
      integer                   :: status       ! NetCDF status variable
      integer                   :: varid_wigis  ! variable id: identifier series
      integer                   :: varid_wigii  ! variable id: issuer of identifier
      integer                   :: varid_wigin  ! variable id: issue number
      integer                   :: varid_wigli  ! variable id: local identifier
      integer                   :: max_ii       ! max.value of issuer of identifier
      integer                   :: i            ! Loop index

      lwigos = .false.
      lwsi   = .false.
      ltsi   = .false.
      select case (wsi_mode)
      case default
         return
      case (1,2)
         max_ii =   999
      case (3)
         max_ii = 19999
      end select
      allocate (wsi(recn), wsihash(recn))
      wsihash = 0_i8

      status = nf90_inq_varid (ncid, 'MWIGIS', varid_wigis)
      if (status /= NF90_NOERR) return
      status = nf90_inq_varid (ncid, 'MWIGII', varid_wigii)
      if (status /= NF90_NOERR) return
      status = nf90_inq_varid (ncid, 'MWIGIN', varid_wigin)
      if (status /= NF90_NOERR) return
      status = nf90_inq_varid (ncid, 'YWIGLI', varid_wigli)
      if (status /= NF90_NOERR) return
      allocate (wigis(recn), wigii(recn), wigin(recn), wigli(recn))
      status = nf90_get_var (ncid, varid_wigis, wigis, start=[rec1])
      status = nf90_get_var (ncid, varid_wigii, wigii, start=[rec1])
      status = nf90_get_var (ncid, varid_wigin, wigin, start=[rec1])
      status = nf90_get_var (ncid, varid_wigli, wigli, start=[rec1])
      !-----------------------------------
      ! Technical validation of variables.
      ! We expect only identifier series 0
      ! Allow traditional IDs via WIGOS ID
      !-----------------------------------
      lwsi =       wigis == 0                           &
             .and. wigii >  0     .and. wigii <= max_ii &
             .and. wigin >= 0     .and. wigin <  65535
      ltsi =       wigis == 0                           &
             .and. wigii >= 20000 .and. wigii <  24000  &
             .and. wigin >= 0     .and. wigin <  65535
      !-------------------------------------------
      ! Local identifier must not start with blank
      !-------------------------------------------
      do i = 1, recn
         if ((lwsi(i) .or. ltsi(i)) .and. wigli(i)(1:1) == " ") then
            lwsi(i) = .false.
            ltsi(i) = .false.
         end if
      end do
      lwigos = any (lwsi) .or. any (ltsi)
      if (.not. lwigos) deallocate (wigis, wigii, wigin, wigli)
      if (.not. lwigos) return
      do i = 1, recn
         if (lwsi(i) .or. ltsi(i)) then
            wsi(i)% valid = .true.
            wsi(i)% wigis = wigis(i)
            wsi(i)% wigii = wigii(i)
            wsi(i)% wigin = wigin(i)
            wsi(i)% wigli = wigli(i)
         end if
      end do
      deallocate (wigis, wigii, wigin, wigli)
    end subroutine check_wigos
    !--------------------------------------------------------------------------
  end subroutine read_obs_netcdf
!------------------------------------------------------------------------------

  subroutine handle_err (where, status)
    character(len=*), intent(in) :: where
    integer,          intent(in) :: status
    if(status /= NF90_NOERR) then
      print *,where,' netcdf-error', TRIM(nf90_strerror(status)), &
              '*** ABORT NOW ***'
      call finish (where,' netcdf error')
    endif
  end subroutine handle_err

!------------------------------------------------------------------------------

  subroutine guess_dbkz (cat, cats, catls, cent, cents, dbkz, filename)
  !------------------------------------------------------------
  ! Guess appropriate Datenbank-Kennzahl for observational data
  ! that do not provide section 2.
  !------------------------------------------------------------
    integer,      intent(in)    :: cat  (:)  ! Data category
    integer,      intent(in)    :: cats (:)  ! Data sub-category
    integer,      intent(in)    :: catls(:)  ! Local data sub-category
    integer,      intent(in)    :: cent (:)  ! Originating center
    integer,      intent(in)    :: cents(:)  ! Originating sub-center
    integer,      intent(inout) :: dbkz (:)  ! DWD DatenBank-KennZahl
    character(*), intent(in)    :: filename  ! Observations file name
    optional                    :: filename

    integer :: kz

    if (size (cat) == 0) return

!   write(*,*) "guess_dbkz: Data category          :", cat(1)
!   write(*,*) "guess_dbkz: Data sub-category      :", cats(1)
!   write(*,*) "guess_dbkz: Local data sub-category:", catls(1)
!   write(*,*) "guess_dbkz: Generating center      :", cent(1)
!   write(*,*) "guess_dbkz: Generating sub-center  :", cents(1)

    !--------------------------------------------------------
    ! First try a lookup in the tables defined in mo_obstypes
    !--------------------------------------------------------
    kz = dbkz_bufr (cat(1), catls(1), cent(1))
    if (netcdf_verb > 1) then
       write(*,*) "guess_dbkz: dbkz_bufr returns dbkz =", kz
    end if

    if (kz < 0) then
       !----------------------------------------------------------------
       ! Last resort: if everything else fails, use hardcoded fallbacks.
       ! Use data categories as defined in WMO Common Code Table C-13
       !----------------------------------------------------------------
       select case (cat(1))
       case (0)         ! Surface data - land
          select case (cats(1))
          case (0)      ! Synoptic observations from automated land stations
             kz = 128
          case (2)      ! Synoptic observations from fixed-land stations (SYNOP)
             kz = 0
          case (10:11)  ! Aeronautical observations (METAR/SPECI)
             kz = 1
          case (14)     ! Ground-based GPS humidity observations (GPSIWV)
             kz = 94
             ! Ground-based GPS slant total delay obs. (GPSGB)
             if (catls(1) == 15) kz = 95
          end select
       case (1)         ! Surface data - sea
          select case (cats(1))
          case (0)      ! Synoptic observations (SHIP)
             kz = 256
          case (6)      ! Automated stations
             kz = 384
          case (25)     ! Buoy observation (BUOY)
             kz = 385
          end select
       case (2)         ! Vertical soundings (other than satellite)
          select case (cats(1))
          case (1)      ! Upper-wind reports from fixed-land stations (PILOT)
             kz = 508
          case (2)      ! Upper-wind reports from ships (PILOT SHIP)
             kz = 768
          case (3)      ! Upper-wind reports from mobile land stations (PILOT MOBIL)
             !???
          case (4)      ! Upper-level reports from fixed-land stations (TEMP)
             kz = 520
          case (5)      ! Upper-level reports from ships (TEMP SHIP)
             kz = 776
          case (6)      ! Upper-level reports from mobile land stations (TEMP MOBIL)
             kz = 524
          case (7)      ! Upper-level reports from dropwindsondes (TEMP DROP)
             kz = 780
          case (10)     ! Wind profiler reports
             kz = 10553
          case (11)     ! RASS temperature reports
!            kz = 556   !       (needs confirmation)
          case (14)     ! Upper-level reports from descent radiosondes,
             kz = 10574 !       fixed-land stations  (TEMP LAND  descent)
          case (15)     ! Upper-level reports from descent radiosondes,
             kz = 10785 !       ships                (TEMP SHIP  descent)
          case (16)     ! Upper-level reports from descent radiosondes,
             kz = 10570 !       mobile land stations (TEMP MOBIL descent)
          case (20)     ! ASDAR/ACARS profiles (AMDAR)
          end select
       case (3)         ! Vertical soundings (satellite)
          select case (cats(1))
!         case (1:7)    ! TOVS/AMSU/HIRS/MHS/IASI
          case (50)     ! Radio occultation sounding
             select case (cent(1))
             case (WMO0_NCAR,WMO0_UCAR,WMO0_DWD,WMO0_CMA,WMO0_RSMC,WMO0_DMI,&
                   WMO0_NASA,WMO0_NOAA,WMO0_SPIRE,WMO0_GEOOPT,WMO0_PLANET   )
                if (catls(1) ==  14) then
                   ! GPS radio occultations (e.g. CDAAC, GFZ test data)
                   kz = 1694
                end if
             case (WMO0_EUMET)
                if (catls(1) ==  14) then
                   ! GPS radio occultations (e.g. METOP-A test data)
                   kz = 1695
                end if
             end select
          case default
             if (cent(1) == WMO0_NCAR .and. cats(1) < 0 .and. catls(1) == 255) then
                ! GPS radio occultations (e.g. SAC-C, C/NOFS, old BUFR data)
                kz = 1694
             end if
          end select
       case (4)         ! Single level upper-air data (other than satellite)
          select case (cats(1))
          case (0)      ! ASDAR/ACARS (AMDAR)
            if (catls(1) ==  146 .or. &
                catls(1) ==  147      ) then
              kz = 542  ! MODES
            else
!             kz = 529  ! AMDAR
              kz = 534  ! ACARS UK MetOffice / Canada
            end if
          case (1)      ! Manual (AIREP, PIREP)
              kz = 530  ! AIREP
          case (2)      ! Mode-S
              kz = 542  !
          end select
       case (5)         ! Single level upper-air data (satellite)
       case (6)         ! Radar data
       case (7)         ! Synoptic features
       case (8)         ! Physical/chemical constituents
       case (9)         ! Dispersal and transport
       case (10)        ! Radiological data
       case (12)        ! Surface data (satellite)
         kz = 99998     ! +++ temporary for testing +++
         if      (cats(1)  ==   7) then
            kz = 1699   ! ASCAT data
         else if (cats(1)  ==   5) then
            kz = 1697   ! Quikscat
         else if (cent(1)  == 254 .and. catls(1) == 99) then
            kz = 1701   ! Altimeter (Jason-like, processed by EUMETSAT)
         else if (cent(1)  ==  99) then
            kz = 1700   ! OSCAT (processed by KNMI)
         end if
       case (21)        ! Radiances (satellite measured)
       case (22)        ! Radar (satellite) but not altimeter and scatterometer
       case (23)        ! Lidar (satellite)
         kz = 99999     ! +++ temporary for testing +++
       case (24)        ! Scatterometry (satellite)
       case (25)        ! Altimetry (satellite)
       case (26)        ! Spectrometry (satellite)
       case (30)        ! Calibration dataset (satellite)
       case (31)        ! Oceanographic data
       case (101)       ! Image data (satellite)
       end select
    end if

    if (kz < 0) then
       if (present (filename)) &
       write(0,*) "guess_dbkz: File                   : ",trim (filename)
       write(0,*) "guess_dbkz: Data category          :", cat(1)
       write(0,*) "guess_dbkz: Data sub-category      :", cats(1)
       write(0,*) "guess_dbkz: Local data sub-category:", catls(1)
       write(0,*) "guess_dbkz: Generating center      :", cent(1)
       write(0,*) "guess_dbkz: Generating sub-center  :", cents(1)
       call finish ("guess_dbkz",'Could not derive DBKZ for given data')
    end if

    dbkz = kz
  end subroutine guess_dbkz
!==============================================================================
end module mo_obs_netcdf
