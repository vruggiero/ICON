!
!+ Interface to observation operator modules from COSMO
!
MODULE mo_cosmo_obs
!
! Description:
!   Interface to observation operator modules from COSMO.
!   It organizes the reading and observation processing of
!   observation reports from NetCDF observation input files.
!
! Method:
!    This module contains the following module procedures:
!    - read_cosmo_obs     : master subroutine, calles by mo_psas.f90
!
!    - obs_cdf_prep_interface
!    - obs_cdf_prep_cma_rules
!    - obs_cdf_print_caution
!    - obs_cdf_print_eventwarn
!
!    It needs the modules
!    - src_obs_cdfin_org:   - obs_cdf_read_org
!                           - obs_cdf_interface
!                           - obs_cdf_proc_init
!                           - obs_cdfin_open_files
!                           - obs_cdf_tower_copy
!                           - obs_cdf_bcorr_rass
!                           - obs_cdf_raso_rh_bias
!                           - obs_cdf_mult_qualicheck
!    - src_obs_cdfin_print: - obs_cdf_print_statist
!                           - obs_cdf_print_events
!                           - obs_cdf_print_reject
!                           - obs_cdf_print_odr
!                           - obs_print_number_rep
!    - src_obs_cdfin_blk:   - obs_read_blacklist
!    - src_obs_proc_air:    - obs_air_org_mult
!    - src_obs_cdfin_util:  - f_p2dp
!
!    - mo_cosmo_obs_data    -  here, parameters are set which
!                             had to be new initialized,
!                             maybe update by namelist necessary ???
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_4         2009/03/26 Andreas Rhodin
!  Template for cosmo observation operator interface
! V1_5         2009/05/25 Andreas Rhodin
!  print some diagnostics
! V1_7         2009/08/24 Andreas Rhodin
!  extend example printout
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Tanja Weusthoff
!  Interface to observation operator modules from COSMO
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_22        2013-02-13 Christoph Schraff
!  adapt to COSMO multilevel observation operator
! V1_42        2015-06-08 Andreas Rhodin
!  read_cosmo_obs: fix bug (local array index adapted to 3dvar convention)
! V1_43        2015-08-19 Andreas Rhodin
!  changes for COSMO MEC
! V1_44        2015-09-30 Andreas Rhodin
!  update shared modules to COSMO 5.03-beta
! V1_47        2016-06-06 Andreas Rhodin
!  add to namelist: altopsu, thairh, mqcorr92, qcc, qccsu, qcciq, qcsiq
!  mqcorr92: change default to 0
! V1_50        2017-01-09 Andreas Rhodin
!  adapt COSMO-MEC to ICON-LAM
! V1_51        2017-02-24 Andreas Rhodin
!  compile with COSMO observation operators
! 2019-06-09 Christoph Schraff: 'lcompute_pe' added to init_par_utilities arglist.
! 2022-08-01 Christoph Schraff: SSO std. deviation collected (for 10-m wind).
! 2023-04-18 CS: 'rhtsat' + redundancy time limits 'rtmlim', 'rtmlrs(y)' etc.
!                converted into namelist parameters + added to obs_cdf_interface.
! 2023-05-16 CS: 'rtmltow', 'icdt_tws', 'itim_wp' added to NL +obs_cdf_interface.
! 2023-10-07 CS: 'icdt_rss' added to NL + obs_cdf_interface.
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!-----------------------------------------------------------------------

!
! Declarations:
!
! Modules used:
!
!-----------------------------------------------------------------------
! ++++++++++++++++++++++++++++++
! from 3dvar_environment
! ++++++++++++++++++++++++++++++

  use mo_kind,                      &
        only : wp                     ! working precision kind parameter

  use mo_exception,                 &
        only : message                ! print routine name and warning message
!       only : finish                 ! abort

  use mo_fortran_units,             &
        only : get_unit_number,     & ! reserve a unit number
               return_unit_number     ! release a unit number

  use mo_mpi_dace,                  &
        only : dace,                & ! MPI group info
               p_barrier,           & ! MPI barrier (for testing)
               p_type,              & ! neue Funktion for parallel_utilities
               p_real_wp,           &
               MPI_INTEGER,         &
               MPI_CHARACTER

  use mo_t_obs,                     &
        only : read_cdfin,          & ! read COSMO observations
               n_dace_op,           & ! length of non-empty entries in dace_op
               dace_op                ! obs types processed by dace operators

  use mo_atm_grid,                  &
        only : t_grid,              & ! derived type for grid definitions
               MO_ICON,             & ! ICON model grid
               print                  ! print content of t_grid variable

  use mo_atm_state,                 &
        only : t_atm,               & ! derived type for atmosphere
               print                  ! print content of t_atm variable

  use mo_atm_decomp,                &
        only : print_decomp           ! print decomposition


  use mo_physics,                   &
        only : gacc,                & ! gravity acceleration
               t0c,                 & ! T [K] for T=0 degree Celsius
               R,                   & ! gas constant of dry air [J/(kg*K)]
               Rd,                  & ! gas constant of water vapour
               RdRd,                & ! R/Rd
               c_p,                 & !  specific heat  [J/(kg*K)]
               d2r                    ! factor degree -> radians (pi/180)

  use mo_time,                      &
        only : t_time, print,       & ! time information
               cyyyymmddhh,         & ! character representation from 't_time'
               cyyyymmddhhmmss,     & !
               hours,               & ! real representation from 't_time'
!              days,                &
               iyyyy,               & ! integer representation from 't_time'
               imm,                 & ! months
               idd,                 & ! days
!              ihh,                 & ! hours
!              ihhmm,               & ! hours and minutes
!              imi,                 & ! minutes
               operator(-)            ! compute time differences

!#undef NO_COSMO_OBS
#ifndef NO_COSMO_OBS

! ++++++++++++++++++++++++++++++
! from "cosmo_environment"
! ++++++++++++++++++++++++++++++

  use data_obs_record                 ! use modules from COSMO


  use data_obs_lib_cosmo,                  &
      only : c0         ,& ! standard real constant 0.0
             c1         ,& ! standard real constant 1.0
             c2         ,& ! standard real constant 2.0
             c05        ,& ! standard real constant 0.5
             c3600      ,& ! standard real constant 3600.0
             rmdich     ,& ! =-1.E30_wp : check value for missing data
             epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

           ! 2b. Pointers for arrays used for reading/process. from NetCDF files
           ! -------------------------------------------------------------------

             i_subpos   ,& ! positions of the subdomains in the total domain
             i_gpscen   ,& ! array of GPS processing centres used actively
             r_av_levs  ,& ! level definition list \ for superobbing layers of
             r_av_incr  ,& ! level increment  list / high-res radiosonde reports
             r_p        ,& ! pressure (at main levels)
             r_hhl      ,& ! height   (at half levels)
             r_t_ll     ,& ! temperature at lowest model (main) level
             r_ps       ,& ! surface pressure
             r_frland   ,& ! land fraction
             r_sso_sd   ,& ! std. dev. of sub-grid scale model orography [m]

             gen_grid   ,& ! generic grid meta data
             lwonl      ,& ! .true. for the node (sub-domain) at which file with
                           !        the unit number 'nupr' is open
                           ! (i.e. where grid pt. (ionl ,jonl ) lies)

           ! 3. Tables for pressure dependent scales
           ! ---------------------------------------
             ncolev     ,& ! number of levels in the correlation scale tables
             tabcolp    ,& ! ln(tabcop(11))
             tabcop     ,& ! levels of the correlation scale tables

           ! 4. I/O device numbers for obs processing / nudging and file names
           ! -----------------------------------------------------------------

             !   file names
             yucautn    ,& ! caution messages if too many obs for ODR size
             yuprint    ,& ! all the remaining information

             !   device numbers
             nucautn    ,& ! caution messages if too many obs for ODR size
             nustat     ,& ! statistics of processed reports
             nurej      ,& ! direct reporting of rejected obs. reports
             nuodr      ,& ! observations stored in the observation data record
             nupr       ,& ! all the ramaining information

           ! ---------------------------------------------
           ! 1. CMA observation type and code type numbers
           !     --> values defined in module ###
           ! ---------------------------------------------

             nsynop     ,& ! SYNOP report
             nairep     ,& ! AIREP report (all aircraft reports)
             nsatob     ,& ! SATOB report
             ndribu     ,& ! DRIBU report
             ntemp      ,& ! TEMP  report
             npilot     ,& ! PILOT report
             nsatem     ,& ! SATEM report
             nsattv     ,& ! SATEM report
             ngps       ,& ! GPS   report
             nscatt     ,& ! SCATT report (from NetCDF only)
             nsrscd     ,& !   synop surface report
             natscd     ,& !   automatic synop surface report
             nshscd     ,& !   ship synop report
!            nabscd     ,& !   ship synop abbreviated report
!            nshred     ,& !   shred report
             natshs     ,& !   automatic ship synop report
             nmetar     ,& !   Metar
             nsytst     ,& !   test synop
             nsytmp     ,& !   surface report from TEMP
             nsytow     ,& !   surface report from tower
             naircd     ,& !   aircraft report
             ncodar     ,& !   codar report
!            ncolba     ,& !   colba report
             namdar     ,& !   amdar report
             nacar      ,& !   acar  report
             nmodes     ,& !   mode-s report
             nstbcd     ,& !   satob report
             nhrvis     ,& !   high-res. VIS wind report
             namv       ,& !   AMV   report
!            nsst       ,& !   sst report
             ndrbcd     ,& !   dribu report
!            nbathy     ,& !   bathy report
             ntesac     ,& !   tesac report (must be the same as 'nscatt')
             nldtcd     ,& !   temp land   report
             nshtcd     ,& !   temp ship   report
             nmotcd     ,& !   temp mobile report
             ntdrop     ,& !   temp drop   report
             nbtemp     ,& !   temp land   report (high-res. BUFR)
             nbship     ,& !   temp ship   report (high-res. BUFR)
             nbdrop     ,& !   temp drop   report (high-res. BUFR)
             ntdesc     ,& !   temp high-res. descending report (HA 02.04.19: as EC)
             nrocob     ,& !   rocob      report
             nrocsh     ,& !   rocob ship report
             nldpcd     ,& !   pilot land   report
             nshpcd     ,& !   pilot ship   report
             nmopcd     ,& !   pilot mobile report
             nwp_eu     ,& !   European wind profiler report
             nra_eu     ,& !   European SODAR/RASS report
             nravad     ,& !   Radar VAD wind report
             ntower     ,& !   tower profile report
             ntowic     ,& !   icos tower profile report
             npr_us     ,& !   US Wind Profiler/RASS report
             nwlidr     ,& !   ground-based wind lidar report
             nstmcd     ,& !   satem report
             nstovs     ,& !   high resolution ATOVS satellite data
!            nsmsg1     ,& !   MSG_1  satellite retrieval (nsmsg1 == 200 + idmsg1)
!            nsmsg2     ,& !   MSG_2  satellite retrieval (nsmsg1 == 200 + idmsg1)
!            nnoa15     ,& !   NOAA15 satellite retrieval (nnoaX == 200 + idnoaaX)
!            nnoa16     ,& !   NOAA16 satellite retrieval
!            nnoa17     ,& !   NOAA17 satellite retrieval
!            nnoa18     ,& !   NOAA18 satellite retrieval
             nascat     ,& !   ASCAT scatterometer report
             nqscat     ,& !   QuickScat scatterometer report
             ngpgfz     ,& !   GPS report processed by GFZ

           ! ---------------------------------------------
           ! 2. CMA obs type and code types
           !    --> values defined in module ###
           ! ---------------------------------------------

             n_cma      ,& ! number of CMA obs and code types
             t_cmatyp   ,& ! data type for info on CMA on obs + code types
             cma        ,& ! array of meta data on CMA observation + code types
             n_gpc      ,& ! number of GPS processing centres
             t_gpscen   ,& ! data type for information on GPS processing centres
             gpc        ,& ! array of meta data on GPS processing centres
             n_gpc_offs ,& ! index offset of GPS code type elements in array cma

           ! ---------------------------------------------
           ! 3. Functions
           !    --> values defined in module ###
           ! ---------------------------------------------

             i_cma      ,& ! function to determine the index of 'cma'
                           ! referring to a given CMA obs and code type
             cs         ,& ! southern latitude limit
             cn         ,& ! northern latitude limit
             cw         ,& ! western longitude limit
             ce         ,& ! eastern longitude limit


          ! ---------------------------------------------
          ! 5. Others
          ! ---------------------------------------------
             lopen_odr  ,& ! .true. if yuobsdr is open
             lopen_rej     ! .true. if yurejct is open


!------------------------------------------------------------------------------

  use data_obs_cdfin,           &

       ! 1 : NetCDF Observation Input File formats
       !     (this section is used ONLY IF obs data are read from NetCDF files)
       ! 1.1   internal attributes of the different NetCDF input files
       !      -------------------------------------------------------

      only: ntype_cdfin,      & ! number of NetCDF observation input files
            mxcdfin          ,& ! max. number of NetCDF observation input files
            icdfinlen        ,& ! maximum length of NetCDF observation input
                                ! file name
            ncdf_temp        ,& ! indicator for processing of NetCDF TEMP
            ncdf_tempship    ,& ! indicator for processing of NetCDF TEMPSHIP
            ncdf_temphirs    ,& ! indicator for proc. NetCDF TEMP high-res BUFR
            ncdf_tempdrop    ,& ! indicator for proc. NetCDF TEMP Dropsonde
            ncdf_tempdesc    ,& ! indicator for proc. NetCDF descending TEMP
            ncdf_pilot       ,& ! indicator for proc. NetCDF PILOT (z-levels)
            ncdf_pilot_p     ,& ! indicator for proc. NetCDF PILOT (p-levels)
            ncdf_amdar_ml    ,& ! indicator for proc. NetCDF AMDAR multi-level
            ncdf_amdar_vp    ,& ! indicator for proc. NetCDF AMDAR vert.profil
            ncdf_amdar       ,& ! indicator for proc. NetCDF AMDAR single-level
            ncdf_acars       ,& ! indicator for proc. NetCDF ACARS single-level
            ncdf_modes       ,& ! indicator for proc. NetCDF MODE-S KNMI format
            ncdf_modes_acr   ,& ! indicator for proc. NetCDF MODE-S ACARS fmt.
            ncdf_wprof       ,& ! indicator for proc. NetCDF wind profiler
            ncdf_rass        ,& ! indicator for proc. NetCDF RASS profiler
            ncdf_radar_vad   ,& ! indicator for proc. NetCDF radar wind prof.
            ncdf_synop       ,& ! indicator for proc. NetCDF SYNOP
            ncdf_synop_mob   ,& ! indicator for proc. NetCDF SYNOP mobile
            ncdf_ship        ,& ! indicator for proc. NetCDF SHIP
            ncdf_buoy        ,& ! indicator for proc. NetCDF BUOY
            ncdf_metar       ,& ! indicator f. proc. NetCDF METAR sfc aviation
            ncdf_gps_zenith  ,& ! indicator for proc. NetCDF GPS (ZPD / IWV)
            ncdf_ascat       ,& ! indicator for proc. NetCDF ASCAT scatterometer
            ncdf_qscat       ,& ! indicator for proc. NetCDF QuickScat scattero.
            ncdf_satob       ,& ! indicator for proc. NetCDF SATOB wind
            ncdf_acars_uk    ,& ! indicator for proc. NetCDF ACARS UK + Canada
            ncdf_acars_us    ,& ! indicator for proc. NetCDF ACARS US w. humid
            ncdf_wlidar_wp   ,& ! indicator for proc. NetCDF ground-b wind lidar
            ncdf_synop_tst   ,& ! indicator for proc. NetCDF test SYNOP
!           ycdfin           ,& ! file names of NetCDF observation input files
            n_cdfin          ,& ! number of existing NetCDF observ. input files
            icdfin              ! obs file type of NetCDF observat. input files
!           ncinid           ,& ! unit nrs of NetCDF observation input files
!           ncid             ,& ! unit number  of current NetCDF file
!           dimids              ! dimension IDs in NetCDF files


USE data_obs_cdfin, ONLY :  &

          ! 3.1    Format of event counters
          !        ------------------------
            mxreve     ,& ! length of report event counter array
            nesodr     ,& ! report number exceeding size of ODR
                          ! (==> adjust Namelist)
            nenoml     ,& ! multi-levl report not made due to ODR array size
            mxdeve     ,& ! length of data event counter array

          ! 3.2    Event counter arrays
          !        --------------------
            neventr    ,& ! counter array of report events
            neventd    ,& ! counter array of data events

          ! 5.1    Redundancy check limits
          !        -----------------------
!           rtmlim     ,& ! time limit for all reports except AIREP      [hrs]
!           rtmlrs     ,& ! time limit for radiosondes (TEMP, PILOT)     [hrs]
!           rtmlair    ,& ! time limit for reports of obs type 'AIREP'   [hrs]
                          !  (time of lowest level of multi-level ODR)

          ! 7. For reporting rejection of data: Output buffer, size and formats
          !    ----------------------------------------------------------------
            outbuf        ! buffer containing output for a single node

!-----------------------------------------------------------------------------

  use mo_cosmo_obs_data,         &

      only: read_nml_cosmo_obs  ,&
            lverpas     ,& ! on - off switch for verif. also of passive reports
            lcloud_ice  ,& ! .t.: on - off switch for cloud_ice in grid-scale
                           !      precipitation (.t. if itype_gscp > 2 in COSMO)
            icdfdirlen  ,& ! max. length of name of directory where NetCDF
                           !   observation input files reside
            mxav        ,& ! max. length of level definition list for
                           !   superobbing of high-resolution radisondes
            mxgpc       ,& ! max. number of GPS processing centres used
            mxtwex      ,& ! max. number of exceptions for tower processing
            mxbcrr      ,& ! max. number of RASS bias correction rules
            ilstidtw    ,& ! max. length of tower station ID's
            ntwex       ,& ! number of exceptions for tower processing
            nbcrr       ,& ! number of RASS bias correction rules

        ! ---------------------------------------------
        ! 0. Reading of observation reports
        ! ---------------------------------------------

            maxmlo      ,& ! max. number of multi-level reports in total domain
            maxsgo      ,& ! max. number of (surface-level and upper-air
                           !      single-level reports within the total domain
            maxgpo      ,& ! max. number of GPS reports within the total domain
!           maxtvo      ,& ! max. number of sat retrievals within total domain
            maxmlv      ,& ! max. number of obs levels in multi-level (m-l) ODR

        ! ---------------------------------------------
        ! 1. Use of stations / reports
        ! ---------------------------------------------

            qcc         ,& ! constant part of the 'quality control thresholds'
            qccsu       ,& ! surface-level data !   \ for upper-air data
            qcvf        ,& ! multiplic. factor to vertically varying part of QCT
            qcciq       ,& ! for integrated water vapour
            qcsiq       ,& ! IWV QC threshold, as a fraction of IWV of saturated
                           !       model temperature profile
            obnlat      ,& !  northern boundary of observation area
            obslat      ,& !  southern boundary of observation area
            obwlon      ,& !  western boundary of observation area
            obelon      ,& !  eastern boundary of observation area
            exnlat      ,& !  northern boundary for exclusion area
            exslat      ,& !  southern boundary for exclusion area
            exwlon      ,& !  western boundary for exclusion area
            exelon      ,& !  eastern boundary for exclusion area
            doromx      ,& ! SYNOP obs. with height differences betw.
                           ! model orography and station height larger
                           ! than 'doromx' are set passive
            altopsu     ,& ! SYNOP obs. above height 'altopsu' are set passive
            zlimv10     ,& ! additional limits for the use of 10-m wind obs:
                           !   pos/neg Laplacian of orography, surface roughness
            dhosag      ,& ! for active use of tower obs profiles:
                           ! limit in terms of height of sensor above ground [m]
                           !   below which model equivalents are computed at the
                           !         observed height of sensor above ground, and
                           !   above which model equivalents are computed at the
                           !         observed height (altitude)
            thairh      ,& ! 20.:maximum horizontal distance [km] betw the
                           ! lowest report and any single level report that
                           ! is added to a multi-level AIRCRAFT report
            av_reso     ,& ! apply superobbing if the averaged resolution of the
                           !   obs profile exceeds 'av_reso' times the model
                           !   resolution
            av_levs     ,& ! level definition list \ for superobbing layers of
            av_incr     ,& ! level increment  list / high-res radiosonde reports
            rhtsat      ,& ! relative humidity threshold above obs is set =100%
            itim_wp     ,& ! mode of correction of obs time for wind profiler
            icdt_tws    ,& ! mode for obs/code type of tower surface-level rep.
            icdt_rss    ,& ! mode for obs/code type of surface report from TEMP
            rtmlrs      ,& ! redundancy time limit for radiosondes         [hrs]
            rtmlsy      ,& ! redundancy time limit for SYNOP               [hrs]
            rtmlair     ,& ! redundancy time limit for AIREP               [hrs]
            rtmltow     ,& ! redundancy time limit for tower               [hrs]
            rtmlrsy     ,& ! redn. time limit for raso sfc. level vs Synop [hrs]
            rtmlim      ,& ! redundancy time limit for other obs           [hrs]
            lredn_repro ,& ! ensure reproducibility of redundancy check
                           !   irrespective of domain decomposition by allowing
                           !   for redundancy only between reports assigned to
                           !   the same grid point
            lsytac      ,& ! prefer TAC over BUFR Synop reports (or vice versa)
            ytwex       ,& ! tower station ID's related to exception levels
            htwex       ,& ! sensor heights a. ground [m] of exception levels
            ivtwex      ,& ! variable indicators related to exception levels:
                           !    (O: none); 1: wind; 2: T + humidity; 3: all
            ybcrr       ,& ! station ID's for RASS bias correction rules
            bcrrt       ,& ! bias correction of RASS (virtual) temperature [K]
            bcrrhl      ,& ! lower \ height limit [m] in RASS profile to apply
            bcrrhu      ,& ! upper / bias correct. 'bcrrt' for station 'ybcrr'

        ! ---------------------------------------------
        ! 2. Use of observation and code type
        ! ---------------------------------------------

            lsynop      ,&  !  if SYNOP data is used
            laircf      ,&  !  if AIREP data is used (aircraft)
            lsatob      ,&  !  if SATOB data is used
            ldribu      ,&  !  if DRIBU data is used (drifting buoy)
            ltemp       ,&  !  if TEMP  data is used
            lpilot      ,&  !  if PILOT data is used
            lsatem      ,&  !  if SATEM data is used
            lgps        ,&  !  if GPS   data is used
            lscatt      ,&  !  if SCATT data is used (scatterometer)
        !   lprodr      ,&  !  .t. for diagnostic print of obs (ODR) data
            lcd011      ,&  !  .t. if synop code  11 data is used (land synop)
            lcd014      ,&  !  .t. if synop code  14 data is used (automatic)
            lcd021      ,&  !  .t. if synop code  21 data is used (ship)
!           lcd022      ,&  !  .t. if synop code  22 data is used (ship abbrev.)
!           lcd023      ,&  !  .t. if synop code  23 data is used (shred)
            lcd024      ,&  !  .t. if synop code  24 data is used (autom. ship)
            lcd140      ,&  !  .t. if synop code 140 data is used (metar)
            lcd811      ,&  !  .t. if synop code 811 data is used (synop test)
            lcd835      ,&  !  .t. if synop code 835 data is used (surface TEMP)
            lcd839      ,&  !  .t. if synop code 839 data is used (surface tower)
            lcd041      ,&  !  .t. if airep code  41 data is used (codar)
            lcd141      ,&  !  .t. if airep code 141 data is used (airep)
!           lcd241      ,&  !  .t. if airep code 241 data is used (colba)
            lcd144      ,&  !  .t. if airep code 144 data is used (amdar)
            lcd146      ,&  !  .t. if airep code 146 data is used (mode-s)
            lcd244      ,&  !  .t. if airep code 244 data is used (acars)
            lcd088      ,&  !  .t. if satob code  88 data is used (satob)
            lcd090      ,&  !  .t. if satob code  90 data is used (amv)
!           lcd188      ,&  !  .t. if satob code 188 data is used (sst)
!           lcd063      ,&  !  .t. if dribu code  63 data is used (bathy)
            lcd064      ,&  !  .t. if dribu code  64 data is used (tesac)
            lcd165      ,&  !  .t. if dribu code 165 data is used (drift. buoy)

            lcd035      ,&  !  .t. if temp  code  35 data is used (land temp)
            lcd036      ,&  !  .t. if temp  code  36 data is used (temp ship)
            lcd037      ,&  !  .t. if temp  code  37 data is used (mobile)
            lcd135      ,&  !  .t. if temp  code 135 data is used (dropsonde)
            lcd109      ,&  !  .t. if temp  code 109 data is used (land temp hi-res)
            lcd111      ,&  !  .t. if temp  code 111 data is used (ship temp hi-res)
            lcd230      ,&  !  .t. if temp  code 230 data is used (drop temp hi-res)
            lcd231      ,&  !  .t. if temp  code 231 data is used (desc temp hi-res)
            lcd039      ,&  !  .t. if temp  code  39 data is used (rocob)
            lcd040      ,&  !  .t. if temp  code  40 data is used (rocob ship)
            lcd032      ,&  !  .t. if pilot code  32 data is used (land pilot)
            lcd033      ,&  !  .t. if pilot code  33 data is used (pilot ship)
            lcd038      ,&  !  .t. if pilot code  38 data is used (mobile)
            lcd132      ,&  !  .t. if pilot code 132 data is used (win-prof eu)
            lcd133      ,&  !  .t. if pilot code 133 data is used (sod/rass eu)
            lcd136      ,&  !  .t. if pilot code 136 data is used (pro/rass us)
            lcd137      ,&  !  .t. if pilot code 137 data is used (Radar VAD)
            lcd139      ,&  !  .t. if pilot code 139 data is used (tower)
            lcd159      ,&  !  .t. if pilot code 159 data is used (icos tower)
            lcd187      ,&  !  .t. if pilot code 187 data is used (wind lidar)
            lcd086      ,&  !  .t. if satem code  86 data is used (satem)
            lcd186      ,&  !  .t. if tovs  code 186 data is used (hi-res ATOVS)
            lcd122      ,&  !  .t. if scatt code 122 data is used (QuickScat)
            lcd123      ,&  !  .t. if scatt code 123 data is used (ASCAT)
!           lcd096      ,&  !  .t. if gps   code  96 data is used
!           mcdmsg1     ,&  !  processing / use of MSG1   code  71 data
!           mcdmsg2     ,&  !  processing / use of MSG2   code  72 data
!           mcdno15     ,&  !  processing / use of NOAA15 code 206 data
!           mcdno16     ,&  !  processing / use of NOAA16 code 207 data
!           mcdno17     ,&  !  processing / use of NOAA17 code 208 data
!           mcdno18     ,&  !  processing / use of NOAA18 code 209 data
            igpscen     ,&  ! array of GPS processing centres used actively

            ycdfdir     ,&  ! directory with NetCDF obs input + blacklist files
            nolbc       ,&  ! number of grid rows at lateral boundaries
                            !   where obs are neglected
            mqcorr92    ,&  ! switch for bias correct. for Vaisala RS92 humidity
                            !    = 0 : no correction for humidity
                            !    = 1 : correct only solar radiation bias
                            !    = 2 : correct total bias (incl. nighttime bias)

            verification_start,& ! start of time window (min.  before time_verif)
            verification_end  ,& ! end of time window (minutes before time_verif)

            ionl         ,& ! / grid point coordinates
            jonl            ! \ for standard output on nudging


  use environment,              &
        only : model_abort,     & ! aborts the program in case of errors
               comm_barrier       ! sets a synchronization point


! use utilities,                &
!       only : diff_minutes       ! compute difference in minutes between
                                  ! 2 dates/times

!------------------------------------------------------------------------------
#ifndef NOMPI

  use parallel_utilities,             &
      only : global_values        ,& ! computes global values by operating on
                                     !   local arrays
             init_par_utilities      ! initializes private variables for
                                     !   module parallel_utilities

#endif

  use src_obs_cdfin_util,             &
      only : f_p2dp                  ! get approx. model layer thickness at
                                     !   given pressure


  use src_obs_proc_air,               &
      only : obs_air_org_mult        ! production of multi-level aircraft
                                     !   reports from single-level reports

  use src_obs_cdfin_print,            &
      only : obs_cdf_print_statist,& ! prints summary of statistics on processed
                                     !   reports
             obs_cdf_print_events ,& ! prints summary of report and data events
             obs_cdf_print_reject ,& ! prints messages on report/ data rejection
             obs_cdf_print_odr    ,& ! prints a part of the ODR
             obs_print_number_rep    ! prints statistics of processed reports
                                     !   per node

  use src_obs_cdfin_org,              &
      only : obs_cdf_read_org     ,& ! organizes reading, processing +
                                     !   distribution of obs
             obs_cdf_interface    ,& ! transfers input values from DACE environ.
                                     !   to the COSMO obs operator modules
             obs_cdf_proc_init    ,& ! allocates long term storage arrays and
                                     !   opens files
             obs_cdfin_open_files ,& ! opens NetCDF observation input files
             obs_cdf_tower_copy   ,& ! applies RASS temperature bias correctionobs_cdf_bcorr_rass
             obs_cdf_bcorr_rass   ,& ! duplicates tower report for alt. proc
             obs_cdf_raso_rh_bias ,& ! bias correction of Vaisala RS92 humidity
             obs_cdf_mult_qualicheck ! model-independent multi-level gross
                                     !   error checking

  use src_obs_cdfin_blk,              &
      only : obs_read_blacklist      ! read blacklist and whitelist from file
                                     !   and store them
#endif

IMPLICIT NONE

!=============================================================================

  !----------------
  ! Public entities
  !----------------
  private
  public :: init_cosmo_obs ! initialize COSMO observation operator modules
  public :: read_cosmo_obs ! organize reading and distribution of reports
  public :: scan_cosmo_obs ! +++ empty routine so far +++

  ! variables related to parallelisation / domain decomposition
  INTEGER        ::  &
    nboundlines      ,& ! number of overlapping boundary lines of the subdomains
    my_cart_id       ,& ! rank of this subdomain in the cartesian communicator
    num_compute      ,& ! number of compute PEs
    icomm_cart       ,& ! communicator for the virtual cartesian topology
    imp_reals        ,& ! REAL      type used for MPI
    imp_integers     ,& ! INTEGER   type used for MPI
    imp_character       ! CHARACTER type used for MPI

  CHARACTER (LEN=255)     ::  &
    yerrmsg             ! error message
  CHARACTER (LEN= 20)     ::  &
    yroutine            ! name of this subroutine

  CHARACTER (LEN= 10)     ::  &
!   ydate_init       ,& ! initial time (start of forecast)
!   ydate_ref        ,& ! refer. date (here:analysis time=verif_time)
!                       ! yyyymmddhh (year, month, day, hour)
    ydate_verif         ! verification time

  CHARACTER (LEN= 14)     ::  &
    ydate_ref           ! refer. date (here:analysis time=verif_time)
                        ! yyyymmddhhmmss (year, month, day, hour, min., sec.)

  REAL (KIND=wp) ::  &
    acthr               ! actual fcst hour (with resp. to 'ydate_ref')


  INTEGER        ::  &
    ntstep           ,& ! timestep, here: actual forecast time in minutes
!   ierror           ,& ! error indicators
    verification_time   ! verification time relative to initial time


!==============================================================================

!------------------------------------------------------------------------------
! Public and Private Subroutines
!------------------------------------------------------------------------------

CONTAINS

!--------------------------------------------------------------------------
! Begin of subroutine init_cosmo_obs
!--------------------------------------------------------------------------

  subroutine init_cosmo_obs
  !----------------------------------------------
  ! initialize COSMO observation operator modules
  !----------------------------------------------
#ifndef NO_COSMO_OBS
    integer :: istat ! error return parameter

    !--------------------------
    ! Initialize some variables
    !--------------------------
    lwonl = .false.
    nupr  = -1

    ! ----------------------------------------------
    ! Read namelist parameter from mo_cosmo_obs_data
    ! ----------------------------------------------
    CALL read_nml_cosmo_obs

    !----------------------------------------------
    ! get file units and open unit 'nupr' (YUPRINT)
    ! (open YUPRINT only if ionl,jonl >= 0)
    !----------------------------------------------
!   if (dace% lpio) then
    if (ionl >= 0 .and. jonl >= 0) &
      nupr    = get_unit_number()  ! reserve a unit number
      nuodr   = get_unit_number()  ! reserve a unit number
      nurej   = get_unit_number()  ! reserve a unit number
      nustat  = get_unit_number()  ! reserve a unit number
      nucautn = get_unit_number()  ! reserve a unit number

!     print*, "nupr ......", nupr   , dace% pe
!     print*, "nuodr .....", nuodr  , dace% pe
!     print*, "nurej .....", nurej  , dace% pe
!     print*, "nustat ....", nustat , dace% pe
!     print*, "nucautn ...", nucautn, dace% pe

    istat = 0
    if (nupr >= 0)                                               &
      OPEN (nupr ,FILE=yuprint,FORM='FORMATTED',STATUS='UNKNOWN' &
                              ,POSITION='APPEND',IOSTAT=istat)
      IF (istat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yuprint FAILED'
        CALL model_abort (my_cart_id, 7005, yerrmsg, yroutine)
      ENDIF
!   endif
#endif
  end subroutine init_cosmo_obs

!--------------------------------------------------------------------------
! Begin of subroutine scan_cosmo_obs
!--------------------------------------------------------------------------

  subroutine scan_cosmo_obs
  !-----------------------------------------------------------
  ! scan the COSMO observation files and set some 3dvar-tables
  !-----------------------------------------------------------
#ifndef NO_COSMO_OBS
#endif
  end subroutine scan_cosmo_obs


!------------------------------------------------------------------------------
! Begin of subroutine read_cosmo_obs
!------------------------------------------------------------------------------

  SUBROUTINE read_cosmo_obs (grid, atm)

    !--------------------------------------------------------------------------
    ! Description:
    !    This subroutine of module "mo_cosmo_obs.f90" organizes the reading and
    !    distribution of observation reports
    !
    ! Method:
    !
    !
    !
    !
    !Current Code Owner: MCH, Tanja Weusthoff
    ! phone: +41 44 256 9 653
    ! email: tanja.weusthoff@meteoswiss.ch
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    ! Subroutine arguments:
    ! --------------------

    type (t_grid) ,pointer    ::  &
         grid                   ! grid definition

    type (t_atm)  ,intent(in) ::  &
         atm                    ! atmospheric state

#ifndef NO_COSMO_OBS

    ! Local scalars:
    !------------------

    LOGICAL        , SAVE ::  &
         lobpfrs = .TRUE.       ! variable for 'first time of obs processing'

    INTEGER        ::  &
         nodrtot  (3)       , & ! total number of multi-level / single-level /
                                ! GPS rep.
         nodrold  (3)       , & ! number of multi-level / single-level / GPS
                                ! reports before having read new data at the
                                ! current timestep (can be modified in
                                ! 'obs_cdf_del_old_rep')
         imaxl    (3)       , & ! size (report dimension) of the (4 types of)
                                ! ODR
         nexceed  (5)       , & ! nr of reports in excess of ODR array size
         nexcedml           , & ! number of rejected multi-level reports due
                                ! to ODR size
         nexceair (2)       , & ! number of multi-level reports derived from
                                ! single-level reports but in excess of array
                                ! size
         nexceold           , & ! nr of multi-level reports rejected at this
                                ! and the previous reading time due to ODR size
         min_sta (mxcdfin)  , & ! start of reading time interval [min]
         min_end (mxcdfin)  , & ! end of reading time interval [min]
!        nncdfin            , & ! if > 0 then atleast 1 NetCDF input file exists
         ntwx1              , & ! max( ntwex, 1 )
         nbcrr1             , & ! max( nbcrr, 1 )
         icdf               , & ! obs file type of NetCDF observation input file
         ilcf    ,kk        , & ! loop indices
         io      ,jo        , & ! local  horizontal coordinates of observation
         nstat, istat       , & ! error status variables
         irep                   ! loop indices

    INTEGER        , SAVE  ::  &
         imaxmll            , & ! size (report dimension) of the multi-level
                                ! (m-l) ODR
         imaxsgl            , & ! size (report dimension) of the single-level
                                ! (s-l) ODR
         imaxgpl            , & ! size (report dimension) of the
                                ! (ground-based) GPS ODR
!        imaxtvl            , & ! size (report dimension) of the satellite
                                ! retrieval ODR
         madj_hum           , & ! mode switch on adjusting observed humidity
                                ! to model
                                ! =1 : adjust observed humidity (by ratio of
                                ! saturation vapour pressure over water to the
                                ! one over ice, to be applied if cloud ice is
                                ! not a state variable
         mxgpc                  ! max. number of GPS processing centres used

    REAL (KIND=wp) ::  &
         rtmlmt           ,& ! colocation threshold for time
         zpob             ,& ! pressure at observation
         fsize  (4)          ! ratio of namelist parameter (max??o)
                             !       to ODR size (max??l)

    LOGICAL        ::  &
         lcompute_pe  ,& ! indicates whether this is a compute PE or not
         ldo_airmult  ,& ! do call 'obs_air_org_mult' because new air-obs read
         ldo_pr_out   ,& ! do call 'obs_cdf_print_odr' because new obs been read
         lcdf         ,& ! read from NetCDF files, which implies:
                         !     - 'nhflag' entry in ODR header has been set, and
                         !     - 'cma' (not 'noctps') is used for statistics
         lsvcorl      ,& ! .t. ==> diminishing of vertical correlation scales
                         !         in the presence of close observations
         lneedob (mxcdfin)    ! new observations need to be read now

    CHARACTER (LEN=20)       :: &
         yroutine       ! name of this subroutine

    ! horizontal and vertical sizes of the model fields
    INTEGER        ::  &
         ie       ,& ! number of grid pts in zonal direction (local sub-domain)
         je       ,& ! number of grid pts in meridional dir. (local sub-domain)
         ke       ,& ! number of grid pts in vertical direction (--> 'f_z2p')
         ie_tot   ,& ! number of grid pts in zonal direction (in total domain)
         je_tot   ,& ! number of grid pts in meridional dir. (in total domain)
         ie_max   ,& ! Max. of ie on all processors
         je_max      ! Max. of je on all processors


    INTEGER        , ALLOCATABLE , TARGET  ::  &
         isubpos  (:,:)    ! positions of the subdomains in the total domain
                           ! (i- and j-indices of the lower left and upper
                           ! right grid point in the order:i_ll,j_ll,i_ur,j_ur;
                           ! only the domain interior is considered, not the
                           ! boundary lines)


    ! constants for the horizontal rotated grid and related variables
    REAL (KIND=wp) ::  &
         pollon       ,& ! longitude of the rotated north pole(in degrees, E>0)
         pollat       ,& ! latitude of the rotated north pole (in degrees, N>0)
         polgam       ,& ! angle between the north poles of the systems
         dlon         ,& ! grid point distance in zonal direction (in degrees)
         dlat         ,& ! grid point distance in merid. direction (in degrees)
         startlat_tot ,& ! transformed latitude of the lower left grid point
                         ! of the total domain (in degrees, N>0)
         startlon_tot ,& ! transformed longitude of the lower left grid point
                         ! of the total domain (in degrees, E>0)
         degrad          ! factor for transforming degree to rad

    ! model fields
    REAL (KIND=wp) , ALLOCATABLE ,  TARGET  ::  &
         hhl   (:,:,:)  ,& ! height at half model levels
         t     (:,:,:)  ,& ! temperature at lowest model (main) level
         fr_land (:,:)  ,& ! land fraction
         sso_sd  (:,:)  ,& ! std. dev. of sub-grid scale model orography SSO [m]
         hsurf   (:,:)  ,& ! height of surface topography  ( m   )
         zp    (:,:,:)  ,& ! model pressure field; zp(ie,je,ke)
         zps     (:,:)     ! model surface pressure field; zps(ie,je)

    ! physical constants
    REAL (KIND=wp) ::  &
         g            ,& ! acceleration due to gravity
         tmelt        ,& ! melting temperature of ice
         r_d          ,& ! gas constant for dry air
         rdv          ,& ! r_d / r_v
         rdocp        ,& ! r_d / cp_d
         b1           ,& ! variables for computing saturation vapour pressure
         b2w          ,& ! over water (w) and ice (i)
         b2i          ,& !               -- " --
         b3           ,& !               -- " --
         b4w          ,& !               -- " --
         b4i             !               -- " --


    ! Local arrays:
    !------------------
    REAL (KIND=wp) , ALLOCATABLE  ::       &
         odp (:)        ! approx. model layer thickness at the obs level
    REAL (KIND=wp) , POINTER      ::       &
         dp0 (:,:,:)    ! model layer thickness
    INTEGER        , ALLOCATABLE  ::       &
         ievnt (:)      ! index for data event counters

    INTEGER   :: iunit

    ! dummy input parameters to obs_cdf_interface :
    INTEGER   :: krun_osse = 0       ! model run to derive obs values from file yfofin='fof'
    LOGICAL   :: llosse_fg = .false. ! f.g. check flag from 'fof' converted to 'dataset flag'
    REAL (wp) :: zfperturb = 0._wp   ! factor to obs error variances to define size of random
                                     ! perturbations added to the obs (only from yfofin='fof')
    INTEGER   :: iseed_ex  = 0       ! external seed for random number generator

    !========= End of Header ==================================================
    !==========================================================================

    yroutine = 'read_cosmo_obs'

    !-------------------------------------------------------
    ! return if COSMO observation processing is switched off
    ! derive surface pressure
    !-------------------------------------------------------
    if (.not. read_cdfin) return

    !--------------------------------------------------------------------------
    ! Section 0: Initalisation of the (MPI) environment and domain decomposion
    !            if not already done
    !--------------------------------------------------------------------------

    ! ### set nexceed to 0, otherwise it has initial value -2147483648 (WHY)
    nexceed(:) = 0


    ! ----------------------------------------
    ! 0.1 Read Information on parallel utilities from mo_atm_grid
    ! ----------------------------------------

    ALLOCATE (isubpos(0:grid% dc% nproc1 * grid% dc% nproc2-1,4))

    num_compute  = grid% dc% npe                  ! number of compute PEs
    nboundlines  = 0                              ! nr of overlapping boundary
                                                  ! lines of the subdomains
    my_cart_id   = grid% dc% myproc               ! rank of this subdomain in
                                                  ! the cartesian communicator
    icomm_cart   = grid% dc% comm                 ! communicator for virtual
                                                  ! cartesian communicator
    isubpos(:,1) = (grid% dc% isubpositions(1,:)) ! positions of the subdomains
    isubpos(:,2) = (grid% dc% isubpositions(3,:)) ! in the total domain
    isubpos(:,3) = (grid% dc% isubpositions(2,:))
    isubpos(:,4) = (grid% dc% isubpositions(4,:))
    imp_reals    = p_real_wp                      ! determines the REAL
                                                  ! type used for MPI
    imp_integers = MPI_INTEGER                    ! determines the INTEGER
                                                  ! type used for MPI
    imp_character= MPI_CHARACTER                  ! determines the CHARACTER
                                                  ! type used for MPI

    ! ----------------------------------------
    ! 0.2 Read Information on model grid from mo_atm_grid
    ! ----------------------------------------
    ie           = grid% ub (1)-grid% lb (1) + 1  ! number of grid points in x
    je           = grid% ub (2)-grid% lb (2) + 1  ! number of grid points in x
    ke           = grid% ub (3)-grid% lb (3) + 1  ! number of levels
    ie_tot       = grid% ubg(1)-grid% lbg(1) + 1  ! number of grid points in x
    je_tot       = grid% ubg(2)-grid% lbg(2) + 1  ! number of grid points in y

    ie_max       = maxval( isubpos(:,3) - isubpos(:,1) + 1)
    je_max       = maxval( isubpos(:,4) - isubpos(:,2) + 1)

    pollon       = grid% dlonr
    pollat       = grid% dlatr
    polgam       = 0                              ! angle between the north
                                                  ! poles of the system
    dlon         = grid% di                       ! grid point distance in
                                                  ! zonal direction in degree
    dlat         = grid% dj                       ! grid point distance in
                                                  ! merid. direction in degree
    startlat_tot = grid% la1
    startlon_tot = grid% lo1
    degrad       = d2r

    gen_grid     => grid                          ! grid meta data for ICON

    ! ----------------------------------------
    ! 0.3 reference time from mo_atm_state...
    ! ----------------------------------------
    ydate_verif  =  cyyyymmddhh(atm% time)         ! verification time
    !ydate_init  =  cyyyymmddhh(atm% ref_time)     ! initial time
    ydate_ref    =  cyyyymmddhhmmss(atm% time)     ! reference time

    acthr  = hours(atm% time - atm% ref_time)
    ntstep = NINT( acthr * 60._wp )

!   ! ----------------------------------------
!   ! 0.4 Read namelist parameter from mo_cosmo_obs_data
!   ! ----------------------------------------
!
!   CALL read_nml_cosmo_obs  moved to init_cosmo_obs
! ! =======================  +++++++++++++++++++++++

    IF ( ALLOCATED(zp))  DEALLOCATE ( hhl, hsurf, fr_land, sso_sd, t, zp, zps )
    ALLOCATE ( hhl     (ie,je,ke+1)                                            &
             , hsurf   (ie,je)                                                 &
             , fr_land (ie,je)                                                 &
             , sso_sd  (ie,je)                                                 &
             , t       (ie,je,ke)                                              &
             , zp      (ie,je,ke)                                              &
             , zps     (ie,je) )

    ! ----------------------------------------
    ! 0.5 Read Information on data fields (from mo_atm_grid.f90 and
    ! mo_atm_state.f90)
    ! ----------------------------------------
    if (dace% lpio) then
       write(6,'(a,i0)') '  Read information on model data fields'
    end if

    !     1. reference atmosphere
    if (associated (grid% hhl)) then
       hhl = grid% hhl(:,:,:,1)      ! geomet. height of half model levels ( m)
    else if (associated (atm% geoh)) then
       hhl = atm% geoh(:,:,:,1)/gacc ! geopot. height of half model levels ( m)
    else
       CALL model_abort (my_cart_id, 9999, "missing: hhl or geoh", yroutine)
    end if

    !     2. external parameter fields
    !        height of surface topography  ( m)
    if (associated (grid% hsurf)) then
       hsurf = grid% hsurf(grid% lb(1):grid% ub(1),grid% lb(2):grid% ub(2),1,1)
    else if (associated (grid% geosp)) then
       hsurf = grid% geosp(grid% lb(1):grid% ub(1),grid% lb(2):grid% ub(2),1,1) / gacc
    else
       CALL model_abort (my_cart_id, 9999, "missing: hsurf or geosp", yroutine)
    end if

    !        fraction of land in a grid element
    if (associated (grid% lsm)) then
       fr_land= grid% lsm  (grid% lb(1):grid% ub(1),grid% lb(2):grid% ub(2),1,1)
    else
       CALL model_abort (my_cart_id, 9999, "missing: fr_land", yroutine)
    end if

    !        std. dev. of sub-grid scale orography SSO
    if (associated (grid% sso_stdh)) then
       sso_sd = grid% sso_stdh (grid% lb(1):grid% ub(1),grid% lb(2):grid% ub(2),1,1)
    else
       if (grid% model == MO_ICON) then
          !-----------------------------------
          ! SSO is a hard requirement for ICON
          !-----------------------------------
          CALL model_abort (my_cart_id, 9999, "missing: sso", yroutine)
       else
          !------------------------------------------------------------
          ! Continue using a missing value, but issue a warning message
          !------------------------------------------------------------
          if (dace% lpio) call message (yroutine, "missing: sso")
          sso_sd = rmdich
       end if
    end if

    !     3. prognostic variables
    t      = atm% t (:,:,:,1)         ! temperature
    zp     = atm% pf(:,:,:,1)
    zps    = atm% ps(:,:,1,1)

    ! -------------------------
    ! 0.6 prepare input variables for CALL obs_cdf_interface
    ! -------------------------
    tmelt  = t0c
    r_d    = r
    rdv    = rdrd
    rdocp  = r / c_p
    g      = gacc

    b1    = 610.78_wp
    b2w   = 17.2693882_wp
    b2i   = 21.8745584_wp
    b3    = 273.16_wp
    b4w   = 35.86_wp
    b4i   = 7.66_wp

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! for testing: print out grid info and atmospheric state
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    call print (grid, verbose=.true.)

    If (dace% lpio) then
       print *
       print *,'grid% dc:'
       call print_decomp (grid% dc)
       print *
       print *,'communicator group                    :', grid% dc% comm
       print *,'Number of PEs in direction 1          :', grid% dc% nproc1
       print *,'Number of PEs in direction 2          :', grid% dc% nproc2
       print *,'nproc1 * nproc2                       :', grid% dc% npe
       print *,'Number (Rank) of own processor element:', grid% dc% myproc
       print *,'Number in direction 1                 :', grid% dc% my_num1
       print *,'Number in direction 2                 :', grid% dc% my_num2
       print *,'Number of neighbor poleward east      :', grid% dc% nbpe
       print *,'Number of neighbor antipoleward west  :', grid% dc% nbaw
       print *,'Number of neighbor poleward west      :', grid% dc% nbpw
       print *,'Number of neighbor antipoleward east  :', grid% dc% nbae
       print *,'Limits of dec. in direction 1         :', grid% dc% ilim1
       print *,'Limits of dec. in direction 2         :', grid% dc% ilim2
       print *
    endif
    call print (atm)

    !-------------------------------------------------------
    ! 0.7 Initialize the utility module parallel_utilities
    !-------------------------------------------------------
#ifndef NOMPI
    ! CS: this assumes that in the DACE environment, all PE's are compute PE's
    lcompute_pe = .true.

    CALL init_par_utilities ( ie, je, ke, ie_tot, je_tot, ke, ie_max, je_max   &
                            , grid% lb(1), grid% ub(1)                         &
                            , grid% lb(2), grid% ub(2)                         &
                            , num_compute, grid% dc% nproc1, grid% dc% nproc2  &
                            , 0, isubpos, nboundlines, icomm_cart, my_cart_id  &
                            , imp_reals, imp_integers, lcompute_pe )
  ! =======================

#endif

!   !-------------------------------------------------- +++++++++++++++++++++++
!   ! 0.8 Get file units and open unit 'nupr' (YUPRINT) moved to init_cosmo_obs
!   !-------------------------------------------------- +++++++++++++++++++++++
!
!   ! get file units and open unit 'nupr' (YUPRINT)
!   !----------------------------------------------
!   IF (lobpfrs)  THEN
!     nupr    = get_unit_number()  ! reserve a unit number
!     nuodr   = get_unit_number()  ! reserve a unit number
!     nurej   = get_unit_number()  ! reserve a unit number
!     nustat  = get_unit_number()  ! reserve a unit number
!     nucautn = get_unit_number()  ! reserve a unit number
!
!     print*, "nupr ......", nupr   , dace% pe
!     print*, "nuodr .....", nuodr  , dace% pe
!     print*, "nurej .....", nurej  , dace% pe
!     print*, "nustat ....", nustat , dace% pe
!     print*, "nucautn ...", nucautn, dace% pe
!   ENDIF
!
!   OPEN (nupr ,FILE=yuprint,FORM='FORMATTED',STATUS='UNKNOWN'                 &
!                           ,POSITION='APPEND',IOSTAT=istat)
!   IF (istat /= 0) THEN
!     yerrmsg = 'OPENING OF FILE yuprint FAILED'
!     CALL model_abort (my_cart_id, 7005, yerrmsg, yroutine)
!   ENDIF

    !------------------------
    ! 0.9 Intialise 'tabcolp'
    !------------------------
    if (lobpfrs)  then
      do kk = 1 , ncolev
        tabcolp (kk) = log( tabcop(kk) )
      enddo
    endif

    !---------------------------------------------------------------------------
    ! Section 1: Interface for the observation pre-processing:
    !            make the 'model environment' available to the COSMO Observation
    !            Processing Library COPL
    !             ('library' denotes here simply a set of modules)
    !             (this basically means: copying values from the 'model environment'
    !              into variables from data modules that are part of COPL)
    !            Remarks:
    !            - The subroutines called here are part of this module, except for
    !              'obs_cdf_interface' itself which is part of the obs proc library.
    !            - The values of all the other variables which are set outside COPL
    !              but which must be available within a subroutine of COPL are
    !              transfered directly from the 'model' environment through the
    !              argument list of the respective routine.
    !            - All the other variables which must be commonly available within
    !              COPL but can be set within the library are stored in one of the
    !              data modules that are part of COPL.
    !-------------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! 1.1 prepare the rules for the use of CMA observation and code types
    !     by filling 'cma' of type 't_cmatyp' as a function of namelist input
    ! --------------------------------------------------------------------------

    IF (lobpfrs)  CALL obs_cdf_prep_cma_rules
    !             ===========================

    if (dace% lpio)  print *,'cma rules prepared'

    ! -------------------------------------------------------------------------
    ! 1.2 prepare those input variables for 'obs_cdf_interface'
    !     which are not yet known here:
    !     - some constant variables depending on namelist input
    !     - some time-dependent model fields
    ! -------------------------------------------------------------------------

    lcdf = .true.

    CALL obs_cdf_prep_interface ( lcdf, madj_hum                               &
                                , imaxmll, imaxsgl, imaxgpl )
  ! ===========================

    if (dace% lpio)  print *,'ods_cdf_interface prepared'

    !---------------------------------------------------------------------------
    ! 1.3 make the 'model environment' available to the obs proc library
    !---------------------------------------------------------------------------

    lwonl  =         (ionl >= grid% lb(1)) .AND. (ionl <= grid% ub(1))         &
               .AND. (jonl >= grid% lb(2)) .AND. (jonl <= grid% ub(2))

!   print *,'test mode...: lwonl = ', lwonl
!   print *,'test mode...: ionl,jonl = ', ionl,jonl, dace% pe
!   print *,'test mode...: lb,ub = ', grid% lb(1),grid% lb(2), grid% ub(1), grid% ub(2)


    CALL obs_cdf_interface ( ie, je, ke, ie_tot, je_tot                        &
                           , num_compute, nboundlines, my_cart_id, icomm_cart  &
                           , imp_reals, imp_integers, imp_character            &
                           , pollon, pollat, polgam, dlon, dlat                &
                           , startlat_tot, startlon_tot, degrad                &
                           , imaxmll, imaxsgl, imaxgpl, maxmlv                 &
                           , nolbc, madj_hum, doromx, altopsu, zlimv10         &
!                          , [1.e10_wp,1.e10_wp,1.e10_wp]                      &! 1.e10 for zlimv10
                           , av_reso, rhtsat, itim_wp, icdt_tws, icdt_rss      &
                           , rtmlrs, rtmlsy, rtmlair, rtmltow                  &
                           , rtmlrsy, rtmlim, lredn_repro, lsytac, acthr       &
                           , krun_osse, llosse_fg, zfperturb, iseed_ex         &
                           , g, tmelt                                          &
                           , r_d, rdv, rdocp, b1, b2w, b2i, b3, b4w, b4i       &
                           , lverpas, lwonl, ydate_ref, mxgpc, mxav )
  ! ======================

    if (dace% lpio)  print *,'ods_cdf_interface executed'

    ! interface for arrays: set 'library 1' pointers to target arrays
    ! (target arrays must have 'target' attribute !)
    ! (this avoids the need to make copies of these arrays, by (de-)allocation)
                                   ! target arrays (must) have dimensions:
    r_p       =>  zp               ! (ie,je,ke)
    r_hhl     =>  hhl              ! (ie,je,ke+1)
    r_t_ll    =>  t (:,:,ke)       ! (ie,je)
    r_ps      =>  zps              ! (ie,je)
    r_frland  =>  fr_land          ! (ie,je)
    r_sso_sd  =>  sso_sd           ! (ie,je)
    r_av_levs =>  av_levs          ! (mxav)
    r_av_incr =>  av_incr          ! (mxav)
    i_subpos  =>  isubpos          ! (0:num_compute-1,4)
    i_gpscen  =>  igpscen          ! (mxgpc)


    !+++++++++++++++++++++++++++++++++++++++
    ! for testing print some stuff to a file
    !+++++++++++++++++++++++++++++++++++++++
    if (dace% lpio) then
       write (6, *)
!      iunit = get_unit_number()  ! reserve a unit number
!      open  (iunit, file = 'COSMO_OBS_TESTFILE')
       iunit = 6                  ! currently write to stdout
       write (iunit, *) 'written by subroutine read_cosmo_obs'
       if (associated (grid% p0)) then
         write (iunit, *) 'size grid% p0           :', size(grid% p0)
         write (iunit, *) 'maxval grid% p0         :', maxval(grid% p0)
       else
         write (iunit, *) 'grid% p0                : not associated'
       endif
       write (iunit, *) 'size atm% ps            :', size(atm% ps(:,:,1,1))
       write (iunit, *) 'maxval atm% ps          :', maxval(atm% ps(:,:,1,1))

       write (iunit, *) 'size zp                 :', size(zp)
       write (iunit, *) 'maxval zp               :', maxval(zp)

       write (iunit, *) 'size atm% pf            :', size(atm% pf(:,:,:,1))
       write (iunit, *) 'maxval atm% pf          :', maxval(atm% pf(:,:,:,1))

       write (iunit, *) 'size zps                :', size(zps)
       write (iunit, *) 'maxval zps              :', maxval(zps)

       if (associated (grid% hhl)) then
        write(iunit, *) 'size grid% hhl          :', size(grid% hhl(:,:,:,1))
        write(iunit, *) 'maxval grid% hhl        :', maxval(grid% hhl(:,:,:,1))
       end if

       write (iunit, *) 'size hhl                :', size(hhl)
       write (iunit, *) 'maxval hhl              :', maxval(hhl)

       write (iunit, *) 'size t                  :', size(t)
       write (iunit, *) 'maxval t                :', maxval(t)

       write (iunit, *) 'size fr_land            :', size(fr_land)
       write (iunit, *) 'maxval fr_land          :', maxval(fr_land)

       write (iunit, *) 'size sso_sd             :', size(sso_sd)
       write (iunit, *) 'maxval sso_sd           :', maxval(sso_sd)

       write (iunit, *) 'ie                      :', ie
       write (iunit, *) 'je                      :', je
       write (iunit, *) 'ke                      :', ke
       write (iunit, *) 'ie_tot                  :', ie_tot
       write (iunit, *) 'je_tot                  :', je_tot

       write (iunit, *) 'num_compute             :', num_compute
       write (iunit, *) 'nboundlines             :', nboundlines
       write (iunit, *) 'my_cart_id              :', my_cart_id
       write (iunit, *) 'icomm_cart              :', icomm_cart
       write (iunit, *) 'isubpos                 :', isubpos
       write (iunit, *) 'imp_reals               :', imp_reals
       write (iunit, *) 'imp_integers            :', imp_integers
       write (iunit, *) 'imp_character           :', imp_character

       write (iunit, *) 'pollon                  :', pollon
       write (iunit, *) 'pollat                  :', pollat
       write (iunit, *) 'polgam                  :', polgam
       write (iunit, *) 'dlon                    :', dlon
       write (iunit, *) 'dlat                    :', dlat
       write (iunit, *) 'startlat_tot            :', startlat_tot
       write (iunit, *) 'startlon_tot            :', startlon_tot
       write (iunit, *) 'degrad                  :', degrad

!CS: this is already written by 'read_nml_cosmo_obs'
!      write (iunit, *) 'qcc                     :', qcc
!      write (iunit, *) 'qccsu                   :', qccsu
!      write (iunit, *) 'qcvf                    :', qcvf
!!     write (iunit, *) 'qcciq                   :', qcciq
!!     write (iunit, *) 'qcsiq                   :', qcsiq
!      write (iunit, *) 'doromx                  :', doromx
!      write (iunit, *) 'altopsu                 :', altopsu
!      write (iunit, *) 'zlimv10                 :', zlimv10
!      write (iunit, *) 'dhosag                  :', dhosag
!      write (iunit, *) 'thairh                  :', thairh
!      write (iunit, *) 'av_levs                 :', r_av_levs
!      write (iunit, *) 'av_incr                 :', r_av_incr
!      write (iunit, *) 'av_reso                 :', av_reso
!      write (iunit, *) 'rhtsat                  :', rhtsat
!      write (iunit, *) 'itim_wp                 :', itim_wp
!      write (iunit, *) 'icdt_tws                :', icdt_tws
!      write (iunit, *) 'icdt_rss                :', icdt_rss
!      write (iunit, *) 'rtmlrs                  :', rtmlrs
!      write (iunit, *) 'rtmlsy                  :', rtmlsy
!      write (iunit, *) 'rtmlair                 :', rtmlair
!      write (iunit, *) 'rtmltow                 :', rtmltow
!      write (iunit, *) 'rtmlrsy                 :', rtmlrsy
!      write (iunit, *) 'rtmlim                  :', rtmlim
!      write (iunit, *) 'lredn_repro             :', lredn_repro
!      write (iunit, *) 'lsytac                  :', lsytac
!      write (iunit, *) 'imaxmll                 :', imaxmll
!      write (iunit, *) 'imaxsgl                 :', imaxsgl
!      write (iunit, *) 'imaxgpl                 :', imaxgpl
!!     write (iunit, *) 'imaxtvl                 :', imaxtvl
!      write (iunit, *) 'maxmlv                  :', maxmlv
!      write (iunit, *) 'mqcorr92                :', mqcorr92
!      write (iunit, *) 'nolbc                   :', nolbc
!      write (iunit, *) 'madj_hum                :', madj_hum
!      write (iunit, *) 'acthr                   :', acthr
!      write (iunit, *) 'ionl                    :',ionl
!      write (iunit, *) 'jonl                    :',jonl

       write (iunit, *) 'gacc                    :', gacc
       write (iunit, *) 't0c                     :', t0c
       write (iunit, *) 'r                       :', r
       write (iunit, *) 'rdrd                    :', rdrd
       write (iunit, *) 'rdocp                   :',rdocp
       write (iunit, *) 'b1                      :', b1
       write (iunit, *) 'b2w                     :',b2w
       write (iunit, *) 'b2i                     :', b2i
       write (iunit, *) 'b3                      :', b3
       write (iunit, *) 'b4w                     :',b4w
       write (iunit, *) 'b4i                     :', b4i

       write (iunit, *) 'lverpas                 :',lverpas
       write (iunit, *) 'lwonl                   :',lwonl
       write (iunit, *) 'ydate_ref               :',ydate_ref

!      close (iunit)
!      call return_unit_number (iunit)
       write (6, *)
    endif

    !---------------------------------------------------------------------------
    ! Section 2: Initialise observation processing, including opening of files,
    !            reading of blacklist, and allocation of long term storage
    !            arrays and of short term storage rejection messaging arrays
    !            ('outbuf'). (The subroutines called here are part of the
    !             obs precess. library.)
    !---------------------------------------------------------------------------

    ! allocation of arrays (ODR, 'nevent*', 'outbuf', 'rolnlv')
    !----------------------------------------------------------

    CALL obs_cdf_proc_init ( lobpfrs , n_cma)
!CS CALL obs_cdf_proc_init ( lobpfrs , n_cma , nuspecif, yuspecif )
  ! ======================

    IF (lobpfrs)  THEN

      CALL obs_read_blacklist ( icdfdirlen , ycdfdir )
    ! =======================

      if (dace% lpio)  print* , 'obs_read_blacklist executed ...'

      n_cdfin = 0

      ! check how many NetCDF observation input files exist for which obs type;
      ! get file annexes and open files (and get their file unit numbers)
      !------------------------------------------------------------------------

      IF (n_dace_op >= 1) THEN

        CALL obs_cdfin_open_files ( icdfdirlen , ycdfdir                       &
                                  , dace_op(1:n_dace_op) )
      ! =========================

      ELSE

        CALL obs_cdfin_open_files ( icdfdirlen , ycdfdir )
      ! =========================

      ENDIF

      ! write a CAUTION to YUCAUTN only if not even 1 NetCDF input file exists
      !-----------------------------------------------------------------------
      IF ((n_cdfin == 0) .AND. (my_cart_id == 0)) THEN
        OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'        &
                                   , POSITION='APPEND', IOSTAT=nstat)
        WRITE( nucautn,'("CAUTION !!!!! NO NetCDF OBSERVATION INPUT FILES IN:" &
                       &,A)' ) ycdfdir(1:LEN_TRIM(ycdfdir))
        PRINT          '("CAUTION !!!!! NO NetCDF OBSERVATION INPUT FILES IN:" &
                       &,A)' , ycdfdir(1:LEN_TRIM(ycdfdir))
        CLOSE (nucautn)
      ENDIF

      ! Set initial number of processed reports to zero
      !------------------------------------------------
      ntotml = 0
      ntotsg = 0
      ntotgp = 0
      nexceold = 0
    ENDIF

    if (dace% lpio)  print *,'ods_cdfin_open_file executed'

    !---------------------------------------------------------------------------
    ! Section 3: Determine whether and for which time interval observations
    !            have to be read now from which type of NetCDF observation
    !            input file
    !---------------------------------------------------------------------------

    do ilcf = 1 , n_cdfin
      icdf  =  icdfin(ilcf)

      !  time interval determined by verification_start and verification_end
      !  from file header of netCDF file; relative to verification_time

!!$       imoyy = iyyyy(atm% ref_time)    ! initial time
!!$       imomm = imm(atm% ref_time)
!!$       imodd = idd(atm% ref_time)
!!$       imohh = ihh(atm% ref_time)
!!$       imomin = imi(atm% ref_time)
!!$       ianyy = iyyyy(atm% time)       ! verification / analysis time
!!$       ianmm = imm(atm% time)
!!$       iandd = idd(atm% time)
!!$       ianhh = ihh(atm% time)
!!$       ianmin = imi(atm% time)
!!$
!!$       ! Time difference in minutes between observation time and initial
!!$       ! model time
!!$       CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin,              &
!!$                           ianyy, ianmm, iandd, ianhh, ianmin,              &
!!$                           verification_time, ierrf )
!!$       !================
!!$
!!$       min_sta (icdf) = verification_time + verification_start
!!$       min_end (icdf) = verification_time + verification_end
!!$       lneedob (icdf) = .TRUE.

      verification_time = 0         ! reference_time  = verification_time
!     min_sta (icdf) = verification_time + verification_start
!     min_end (icdf) = verification_time + verification_end
!     lneedob (icdf) = .TRUE.
      min_sta (ilcf) = verification_time + verification_start
      min_end (ilcf) = verification_time + verification_end
      lneedob (ilcf) = .TRUE.

      ! the period from which observations need to be read has to be extended
      ! by the temporal redundancy check limits in order to allow for correct
      ! redundancy checking across analysis time steps
      ! (otherwise 2 reports e.g. only 1 minute apart may be actively used,
      !  one each in subsequent analysis steps, even though the 2 reports should
      !  be redundant against each other and only one of them be used actively -
      !  this would often occur e.g. for Swiss wind profiler reports)
      rtmlmt  =  rtmlim
      IF (     (icdf == ncdf_temp    ) .OR. (icdf == ncdf_tempship)            &
          .OR. (icdf == ncdf_temphirs) .OR. (icdf == ncdf_tempdesc)            &
                                       .OR. (icdf == ncdf_tempdrop)            &
          .OR. (icdf == ncdf_pilot   ) .OR. (icdf == ncdf_pilot_p )) THEN
        rtmlmt         =      rtmlrs
      ELSEIF (     (icdf == ncdf_amdar)    .OR. (icdf == ncdf_acars)           &
              .OR. (icdf == ncdf_acars_uk) .OR. (icdf == ncdf_acars_us)        &
              .OR. (icdf == ncdf_modes)    .OR. (icdf == ncdf_modes_acr)       &
              .OR. (icdf == ncdf_amdar_ml) .OR. (icdf == ncdf_amdar_vp)) THEN
        rtmlmt         =      rtmlair
      ENDIF
      min_sta (ilcf) = min_sta(ilcf) - NINT( rtmlmt * 60._wp )
      min_end (ilcf) = min_end(ilcf) + NINT( rtmlmt * 60._wp )
    enddo

    if (dace% lpio) then
       print*, 'ydate_verif ', ydate_verif
       print*, 'ydate_ref ',   ydate_ref
       print*, 'verification_time ', verification_time
       print*, 'min_sta(1) ', min_sta(1)
       print*, 'min_end(1) ', min_end(1)
    endif


    !---------------------------------------------------------------------------
    ! Section 4: Read from NetCDF files, select, assign to grid points, and
    !            distribute observational reports to nodes (processors)
    !            according to the sub-domains which contain the assigned grid
    !            points; then pre-process and store the reports in the ODR
    !            arrays; finally check for redundancy of reports; also clean
    !            up the ODR arrays from reports which are not used any more
    !---------------------------------------------------------------------------

    ! get numbers of reports which have been read prior to the current timestep
    !--------------------------------------------------------------------------
    nodrold (1) = ntotml
    nodrold (2) = ntotsg
    nodrold (3) = ntotgp
    nexceed     =  0
    ldo_airmult = .FALSE.
    ldo_pr_out  = .FALSE.

    if (n_cdfin >= 1) then

      CALL obs_cdf_read_org ( lneedob(1:n_cdfin) , min_sta(1:n_cdfin) &
                            , min_end(1:n_cdfin)                      &
                            , ycdfdir, icdfdirlen                     &
                            , nodrold, nexceed(1:3)                   &
                            , ldo_airmult, ldo_pr_out                 )
    ! =====================

      if (dace% lpio)  print*, 'obs_cdf_read_org executed'

      ! bias correction depending on solar zenith angle
      ! for Vaisala RS92 radiosonde humidity
      ! -----------------------------------------------

      CALL obs_cdf_raso_rh_bias ( mqcorr92, nodrold(1), ntotml )
    ! =========================

      ! bias correction of RASS (virtual) temperature
      ! ---------------------------------------------

      nbcrr1 = MAX( nbcrr, 1 )
  
      CALL obs_cdf_bcorr_rass ( nodrold(1), ntotml, nbcrr, ilstidtw            &
                              , bcrrt (1:nbcrr1), bcrrhl(1:nbcrr1)             &
                              , bcrrhu(1:nbcrr1), ybcrr (1:nbcrr1) )
    ! =======================

    endif

    !---------------------------------------------------------------------------
    ! Section 5: Produce vertical profiles (multi-level reports) from
    !            single-level aircraft reports, and save total number of
    !            stored ODR reports
    !---------------------------------------------------------------------------

    nexceair (1) = 0
    nexceair (2) = 0

    IF (ldo_airmult) THEN
      lcdf    = .TRUE.
      lsvcorl = .TRUE.
      ALLOCATE (odp (MAX(ntotsg,1)) , STAT=istat )
      IF (associated (grid% dp0)) then
        !----------------------------------------
        ! use thickness from reference atmosphere
        !----------------------------------------
        dp0 => grid% dp0 (:,:,:,1)
      ELSE
        !-------------------------------------------
        ! reference atmosphere is not present
        ! estimate thickness from actual model state
        !-------------------------------------------
        allocate (dp0 (grid% ub(1) - grid% lb(1) + 1, &
                       grid% ub(2) - grid% lb(2) + 1, &
                       grid% ub(3) - grid% lb(3) + 1) )
        dp0 (:,:,1)      =  zp (:,:,2)    - zp (:,:,1)
        dp0 (:,:,ke)     =  zp (:,:,  ke) - zp (:,:,  ke-1)
        dp0 (:,:,2:ke-1) = (zp (:,:,3:ke) - zp (:,:,1:ke-2)) / 2._wp
      ENDIF
      DO irep = 1 , ntotsg
        io   = mosghd (irep,nhio)    ! local gridpoint indices
        jo   = mosghd (irep,nhjo)
        zpob = osgbdy (irep,nbsp)
        odp (irep)  =  f_p2dp ( zpob , ke , zp(io,jo,:) , &
                                           dp0(io,jo,:)   )
      ENDDO

      CALL obs_air_org_mult ( lsvcorl, nodrold(2), nodrold(1), thairh &
                            , lobpfrs, isubpos , nexceair , odp )
    ! =====================

      DEALLOCATE (odp , STAT=istat )
      nexceed (1)  =  nexceed(1) + nexceair(1)
    ENDIF

    !---------------------------------------------------------------------------
    ! Section 6: Model-independent quality control check for multi-level data
    !---------------------------------------------------------------------------

    IF (ntotml > nodrold(1)) THEN
      ALLOCATE ( ievnt (ntotml+1) , STAT=istat )
      DO irep = nodrold(1)+1 , ntotml
        ievnt (irep) = i_cma ( momlhd(irep,nhobtp), momlhd(irep,nhcode) )
        !               =====
      ENDDO

      CALL obs_cdf_mult_qualicheck ( nodrold(1), ntotml, maxmlv, ievnt, rdocp )
    ! ============================

      if (dace% lpio)  print*, 'obs_cdf_mult_qualicheck executed !'

      DEALLOCATE ( ievnt , STAT=istat )
    ENDIF

    !---------------------------------------------------------------------------
    ! Section 6b: Duplicating tower reports
    !             to allow for alternative observation operators
    !---------------------------------------------------------------------------

    ntwx1 = MAX( ntwex, 1 )

    CALL obs_cdf_tower_copy ( dhosag, nodrold(1), ntwex, ilstidtw              &
                            , htwex(1:ntwx1), ivtwex(1:ntwx1), ytwex(1:ntwx1)  &
                            , ntotml, nexceed(1) )
  ! =======================

    !---------------------------------------------------------------------------
    ! Section 7: Print out ODR and diagnostic output (+ de-allocate 'outbuf')
    !---------------------------------------------------------------------------

    imaxl (1)    =  imaxmll
    imaxl (2)    =  imaxsgl
    imaxl (3)    =  imaxgpl
!   imaxl (4)    =  imaxtvl
    nodrtot (1)  =  ntotml
    nodrtot (2)  =  ntotsg
    nodrtot (3)  =  ntotgp

!   print*, 'lopen_odr = ',lopen_odr, dace% pe
!   print*, 'nuodr = ',nuodr
!   print*, 'nodrtot = ', dace% pe, nodrtot(1), nodrold(1), imaxl(1),          &
!                       nexceed(1), nodrtot(2), nodrold(2), imaxl(2), lopen_odr

    IF (ldo_pr_out) THEN

      CALL obs_cdf_print_reject ( lwonl, num_compute, my_cart_id, icomm_cart   &
                                , imp_integers )
    ! =========================

      CALL obs_cdf_print_odr ( nodrtot , nodrold , num_compute , my_cart_id    &
                             , icomm_cart , imp_integers , 0 )
    ! ======================

      CALL obs_print_number_rep ( nodrtot, nodrold, imaxl, lwonl, num_compute  &
                                , my_cart_id, icomm_cart, imp_integers )
    ! =========================

    ELSE
      DEALLOCATE ( outbuf , STAT=istat )
    ENDIF

    ! PRINT *,'after  obs_cdf_print_odr ', nodrold


    !---------------------------------------------------------------------------
    ! Section 9: Write alerting ('CAUTION') messages if the ODR array size is
    !            insufficient to accommodate all observation reports;
    !            the messages include a guess by how much certain namelist
    !       parameters should be increased to render the ODR size large enough.
    !          (Due to the namelist dependency (e.g. 'maxmlo'), the called
    !          routines are part of this module rather than the cdfin-library.)
    !---------------------------------------------------------------------------

    ! write alerting messages to a specific file (yucautn), if required
    ! ----------------------------------------------------

    nexcedml     =             nexceed(1)
    nexceed (5)  =  nexceold + nexceed(1)
    nexceold     =  nexcedml

    fsize (1)    =  REAL( maxmlo ,wp ) / MAX( REAL( imaxmll ,wp ) ,c1 )
    fsize (2)    =  REAL( maxsgo ,wp ) / MAX( REAL( imaxsgl ,wp ) ,c1 )
    fsize (3)    =  REAL( maxgpo ,wp ) / MAX( REAL( imaxgpl ,wp ) ,c1 )
!   fsize (4)    =  REAL( maxtvo ,wp ) / MAX( REAL( imaxtvl ,wp ) ,c1 )

    CALL obs_cdf_print_caution ( nexceed , nexceair , imaxl , fsize )
  ! ==========================


    !---------------------------------------------------------------------------
    ! Section 10: Printing of diagnostics: Statistics on processed reports and
    !             events
    !---------------------------------------------------------------------------


    ! Print statistics on the processed observation reports (and open 'nustat')
    !--------------------------------------------------------------------------

    CALL obs_cdf_print_statist ( 1 , num_compute                               &
                                   , my_cart_id , icomm_cart , imp_integers )
  ! ==========================

    ! further summarise warning messages (to unit number 'nustat')
    !-------------------------------------------------------------

    CALL obs_cdf_print_eventwarn ( mxreve , neventr , fsize )
  ! ============================

    ! Print summary of report and data events (and finally close 'nustat')
    !---------------------------------------------------------------------

    CALL obs_cdf_print_events (  0 , 0 , mxreve , neventr , num_compute       &
                                      , my_cart_id , icomm_cart , imp_integers)
  ! =========================
    CALL obs_cdf_print_events (  0 , 1 , mxdeve , neventd , num_compute       &
                                      , my_cart_id , icomm_cart , imp_integers)
  ! =========================
    CALL obs_cdf_print_events ( -1 , 2 , mxdeve , neventd , num_compute       &
                                      , my_cart_id , icomm_cart , imp_integers)
  ! =========================

    DEALLOCATE ( neventr   , STAT = istat )
    DEALLOCATE ( neventd   , STAT = istat )

    IF (my_cart_id == 0) THEN
      IF (lopen_odr)  CLOSE (nuodr)
      IF (lopen_rej)  CLOSE (nurej)
      lopen_odr = .FALSE.
      lopen_rej = .FALSE.
    ENDIF
!   CLOSE ( nupr )

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! for testing: wait for all processors to finish; abort
!   !++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   print *, 'my_cart_id at finish ', my_cart_id
!
!   call p_barrier
!   call finish ('read_cosmo_obs','test')

    lobpfrs = .FALSE.

    !---------------------------------------------------------------------------
    ! end subroutine read_cosmo_obs
    !---------------------------------------------------------------------------

#endif

  END SUBROUTINE read_cosmo_obs

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

#ifndef NO_COSMO_OBS



!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_prep_interface ( lcdf, madj_hum                             &
                                  , imaxmll, imaxsgl, imaxgpl )

!-------------------------------------------------------------------------------
! Description:
!   This procedure of module "src_obs_proc_cdf" prepares those input variables
!   for the interface 'obs_cdf_interface' to the observation processing library
!   which are not yet known:
!    - some constant variables depending on namelist input
!    - some time-dependent model fields
!
! Method:
!   Use namelist parameters and model fields from the 'model' environment.
!   16.01.09
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------

  LOGICAL        , INTENT (IN)     ::  &
    lcdf                ! read conventional obs from NetCDF files, implies:
                        !   = .true.: this routine is always called
                        !   = .false: this routine is only called for 1DVAR

  INTEGER        , INTENT (OUT)    ::  &
!   imaxtvl          ,& ! size (report dimension) of the satellite retrieval ODR
    madj_hum            ! = 1 : adjust observed humidity (by ratio of saturation
                        !       vapour pressure over water to the one over ice,
                        !       applied if cloud ice is not a state variable

  INTEGER        , INTENT (OUT) , OPTIONAL   ::  &
    imaxmll          ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    imaxsgl          ,& ! size (report dimension) of the single-level (s-l)  ODR
    imaxgpl             ! size (report dimension) of the (ground-based) GPS  ODR

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  LOGICAL        , SAVE  ::  &
    lprepfrs = .TRUE.   ! variable for 'first time of obs. processing'

  REAL (KIND=wp) ::  &
    zfsize           ,& ! factor used to set array sizes
    zfproc              ! additional factor dependent on number of processors

! CHARACTER (LEN=20)       ::  &
!   yroutine            ! name of this subroutine
! CHARACTER (LEN=30)       ::  &
!   yerr                ! error message

! Local arrays:
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_prep_interface
!-------------------------------------------------------------------------------

! yroutine = 'obs_cdf_prep_interface'

!-------------------------------------------------------------------------------
! Section 1: Constant values, dependent on namelist parameters
!-------------------------------------------------------------------------------

  IF (lprepfrs) THEN

! determine the report dimensions of the Observation data record (ODR)
! --------------------------------------------------------------------
! (i.e. the maximum numbers (imaxmll, imaxsgl, imaxpgl, imaxtvl) of reports
!  in a subdomain as a function of the number of nodes, the temporal
!  weight functions, and the total number of reports as specified in
!  the namelist input (by: maxmlo, maxsgo, maxpgo, maxtvo))

    zfproc  = MIN( c1 , c1/10 + 9/(10 *SQRT(SQRT( c1*num_compute ))) )
!   zfproc  = MIN( c1 , c1/5  + 4/(5  *     SQRT( c1*num_compute ) ) )

!   to be adjusted:
    zfsize  = c1
!   imaxtvl =  MAX( 1 , NINT( maxtvo * zfproc * zfsize ) )

    IF ((lcdf) .AND. (PRESENT( imaxsgl )) .AND. (PRESENT( imaxmll ))           &
               .AND. (PRESENT( imaxgpl ))) THEN
!     zfsize  = MIN( (  MAX( wtukrsa, wtukara, wtuksua, tipolmx, tipmxsu )     &
!                     + MAX( wtuksue, tipmxsu ) + c1)                          &
!                   /   MAX( wtuksua+wtuksue, c2*tipmxsu, c1 )  ,  4.0_wp )
      imaxsgl = MAX( 1 , NINT( maxsgo * zfproc * zfsize ) )

!     zfsize  = MAX( wtukrsa+wtukrse, wtukara+wtukare, c2*tipolmx, c1 )
!     zfsize  = (zfsize + c1)      / zfsize
      imaxmll = MAX( 1 , NINT( maxmlo * zfproc * zfsize ) )

!     zfsize  = MAX( c1 , MIN( nstop, nudgend ) *dtdeh )
      imaxgpl = MAX( 1 , NINT( maxgpo * zfproc * zfsize ) )
      IF (lwonl) THEN
        WRITE( nupr,'(" MAXIMUM NUMBER OF REPORTS IN THE TOTAL DOMAIN:",       &
                     &"   MAXMLO=",I5,"  MAXSGO=",I5, /,                       &
                     &" MAXIMUM NUMBER OF REPORTS IN SUBDOMAINS:      ",       &
                     &"   MAXMLL=",I5,"  MAXSGL=",I5)' )                       &
               maxmlo, maxsgo, imaxmll, imaxsgl
        WRITE( nupr,'(" MAXIMUM NUMBER OF REPORTS IN THE TOTAL DOMAIN:",       &
                     &"   MAXGPO=",I5, 14X, /,                                 &
                     &" MAXIMUM NUMBER OF REPORTS IN SUBDOMAINS:      ",       &
                     &"   MAXGPL=",I5     )' )                                 &
               maxgpo,         imaxgpl
!       WRITE( nupr,'(" MAXIMUM NUMBER OF REPORTS IN THE TOTAL DOMAIN:",       &
!                    &"   MAXGPO=",I5,"  MAXTVO=",I5, /,                       &
!                    &" MAXIMUM NUMBER OF REPORTS IN SUBDOMAINS:      ",       &
!                    &"   MAXGPL=",I5,"  MAXTVL=",I5)' )                       &
!              maxgpo, maxtvo, imaxgpl, imaxtvl
      ENDIF
!   ELSEIF (lwonl) THEN
!       WRITE( nupr,'(" MAXIMUM NUMBER OF RETRIEVALS IN THE TOTAL DOMAIN:",    &
!                    &"   MAXTVO=",I5,/,                                       &
!                    &" MAXIMUM NUMBER OF RETRIEVALS IN SUBDOMAINS:      ",    &
!                    &"   MAXTVL=",I5)' )   maxtvo, imaxtvl
    ENDIF

! decide whether humidity observation values need to be adjusted
! to become compatible with the model
! (this is the case if cloud ice is not a model state variable)
! --------------------------------------------------------------------

    madj_hum  =  0
!   IF (itype_gscp <= 2)  madj_hum  =  1
    IF (.NOT. lcloud_ice)  madj_hum  =  1

  ENDIF

  lprepfrs = .FALSE.

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_prep_interface
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_prep_interface


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_prep_cma_rules

!-------------------------------------------------------------------------------
! Description:
!   This procedure of module "src_obs_proc_cdf" prepares rules for the use
!   observation and code types.
!
! Method:
!   Array 'cma' of type 't_cmatyp' is updated by the rules that depend on
!   latitute and longitude are written, using namelist parameters.
!   02.01.09
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    icma, ic, icn       ! loop indices

  LOGICAL        ::  &
    lgpc  (0:n_gpc)     ! .TRUE. if GPS processing centre shall be used actively
!   lcd097           ,&
!   lcd098           ,&
!   lcd099           ,&
!   lcd100           ,&
!   lcd101           ,&
!   lcd102           ,&
!   lcd103           ,&
!   lcd104           ,&
!   lcd105           ,&
!   lcd106           ,&
!   lcd107           ,&
!   lcd108

  CHARACTER (LEN=12)       ::  &
    ybuf                ! name of this subroutine
! CHARACTER (LEN=20)       ::  &
!   yroutine            ! name of this subroutine
! CHARACTER (LEN=30)       ::  &
!   yerr                ! error message

! Local arrays:
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_prep_cma_rules
!-------------------------------------------------------------------------------

#ifndef NO_COSMO_OBS

! yroutine = 'obs_cdf_prep_cma_rules'

  ! complete table 'cma' (for GPS code types, i.e. sub-centres)
  DO ic = 1 , n_gpc
    cma (n_gpc_offs+ic) %obtyp  =  cma (n_gpc_offs) %obtyp
    cma (n_gpc_offs+ic) %cdtyp  =  gpc(ic) %cdtyp
!   ybuf  =  cma (n_gpc_offs+ic) %name (1:12)
    ybuf  =  cma (n_gpc_offs) %name (1:12)
    cma (n_gpc_offs+ic) %name   =  ybuf(1:LEN_TRIM(ybuf)) //  ' by '           &
                                                          // gpc(ic) %name
  ENDDO

  ! determine which centres in 'gpc' shall be used actively
  lgpc (0)       =  lgps
  lgpc (1:n_gpc) = .FALSE.
  DO icn = 1 , mxgpc
    IF (igpscen(icn) >= 0) THEN
      DO ic = 1 , n_gpc
        IF (gpc(ic) %gpscen == igpscen(icn))  lgpc (ic) = .TRUE.
      ENDDO
    ENDIF
  ENDDO


  ! preparation of use of GPS processing centres: if centre is given in array
  ! 'igpscen', a corresponding logical lcdxxx is set true for active use;
  ! if centre is not in 'igpscen' then all data of this centre are set passive
! lcd108=.false.
! lcd097=.false.
! lcd098=.false.
! lcd099=.false.
! lcd100=.false.
! lcd101=.false.
! lcd102=.false.
! lcd103=.false.
! lcd104=.false.
! lcd105=.false.
! lcd106=.false.
! lcd107=.false.
! DO ic = 1 , mxgpc
!   IF (igpscen(ic) == 23)                           lcd108=.TRUE.
!   IF (igpscen(ic) == 21)                           lcd097=.TRUE.
!   IF (igpscen(ic) == 30)                           lcd098=.TRUE.
!   IF (igpscen(ic) == 24)                           lcd099=.TRUE.
!   IF (igpscen(ic) == 25)                           lcd100=.TRUE.
!   IF (igpscen(ic) == 29)                           lcd101=.TRUE.
!   IF (igpscen(ic) == 26)                           lcd102=.TRUE.
!   IF (igpscen(ic) ==  0)                           lcd103=.TRUE.
!   IF((igpscen(ic) == 32) .OR. (igpscen(ic) == 37)) lcd104=.TRUE.
!   IF (igpscen(ic) == 35)                           lcd105=.TRUE.
!   IF (igpscen(ic) == 33)                           lcd106=.TRUE.
!   IF (igpscen(ic) == 34)                           lcd107=.TRUE.
! ENDDO

!-------------------------------------------------------------------------------

  DO icma = 1 , n_cma
    cma(icma)%obslat = MAX( cma(icma)%obslat , obslat )
    cma(icma)%obnlat = MIN( cma(icma)%obnlat , obnlat )
    cma(icma)%obwlon = MAX( cma(icma)%obwlon , obwlon )
    cma(icma)%obelon = MIN( cma(icma)%obelon , obelon )
    IF (cma(icma)%obtyp == nsynop) THEN
      IF (     (.NOT. lsynop)                                                  &
          .OR. ((cma(icma)%cdtyp == nsrscd) .AND. (.NOT. lcd011))              &
          .OR. ((cma(icma)%cdtyp == natscd) .AND. (.NOT. lcd014))              &
          .OR. ((cma(icma)%cdtyp == nshscd) .AND. (.NOT. lcd021))              &
!         .OR. ((cma(icma)%cdtyp == nabscd) .AND. (.NOT. lcd022))              &
!         .OR. ((cma(icma)%cdtyp == nshred) .AND. (.NOT. lcd023))              &
          .OR. ((cma(icma)%cdtyp == natshs) .AND. (.NOT. lcd024))              &
          .OR. ((cma(icma)%cdtyp == nsytst) .AND. (.NOT. lcd811))              &
          .OR. ((cma(icma)%cdtyp == nsytmp) .AND. (.NOT. lcd835))              &
          .OR. ((cma(icma)%cdtyp == nsytow) .AND. (.NOT. lcd839))              &
          .OR.  (cma(icma)%cdtyp == nmetar)                      ) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == nairep) THEN
      IF (     (.NOT. laircf)                                                  &
          .OR. ((cma(icma)%cdtyp == ncodar) .AND. (.NOT. lcd041))              &
          .OR. ((cma(icma)%cdtyp == naircd) .AND. (.NOT. lcd141))              &
          .OR. ((cma(icma)%cdtyp == namdar) .AND. (.NOT. lcd144))              &
          .OR. ((cma(icma)%cdtyp == nacar ) .AND. (.NOT. lcd244))              &
          .OR. ((cma(icma)%cdtyp == nmodes) .AND. (.NOT. lcd146))) THEN
!         .OR. ((cma(icma)%cdtyp == ncolba) .AND. (.NOT. lcd241))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == nsatob) THEN
      IF (     (.NOT. lsatob)                                                  &
          .OR. ((cma(icma)%cdtyp == nstbcd) .AND. (.NOT. lcd088))              &
          .OR.  (cma(icma)%cdtyp == nhrvis)                                    &
          .OR. ((cma(icma)%cdtyp == namv  ) .AND. (.NOT. lcd090))) THEN
!         .OR. ((cma(icma)%cdtyp == nsst  ) .AND. (.NOT. lcd188))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == ndribu) THEN
      IF (     (.NOT. ldribu)                                                  &
!         .OR. ((cma(icma)%cdtyp == nbathy) .AND. (.NOT. lcd063))              &
          .OR. ((cma(icma)%cdtyp == ntesac) .AND. (.NOT. lcd064))              &
          .OR. ((cma(icma)%cdtyp == ndrbcd) .AND. (.NOT. lcd165))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == ntemp ) THEN
      IF (     (.NOT. ltemp )                                                  &
          .OR. ((cma(icma)%cdtyp == nldtcd) .AND. (.NOT. lcd035))              &
          .OR. ((cma(icma)%cdtyp == nshtcd) .AND. (.NOT. lcd036))              &
          .OR. ((cma(icma)%cdtyp == nmotcd) .AND. (.NOT. lcd037))              &
          .OR. ((cma(icma)%cdtyp == ntdrop) .AND. (.NOT. lcd135))              &
          .OR. ((cma(icma)%cdtyp == nbtemp) .AND. (.NOT. lcd109))              &
          .OR. ((cma(icma)%cdtyp == nbship) .AND. (.NOT. lcd111))              &
          .OR. ((cma(icma)%cdtyp == nbdrop) .AND. (.NOT. lcd230))              &
          .OR. ((cma(icma)%cdtyp == ntdesc) .AND. (.NOT. lcd231))              &
          .OR. ((cma(icma)%cdtyp == nrocob) .AND. (.NOT. lcd039))              &
          .OR. ((cma(icma)%cdtyp == nrocsh) .AND. (.NOT. lcd040))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == npilot) THEN
      IF (     (.NOT. lpilot)                                                  &
          .OR. ((cma(icma)%cdtyp == nldpcd) .AND. (.NOT. lcd032))              &
          .OR. ((cma(icma)%cdtyp == nshpcd) .AND. (.NOT. lcd033))              &
          .OR. ((cma(icma)%cdtyp == nmopcd) .AND. (.NOT. lcd038))              &
          .OR. ((cma(icma)%cdtyp == nwp_eu) .AND. (.NOT. lcd132))              &
          .OR. ((cma(icma)%cdtyp == nra_eu) .AND. (.NOT. lcd133))              &
          .OR. ((cma(icma)%cdtyp == npr_us) .AND. (.NOT. lcd136))              &
          .OR. ((cma(icma)%cdtyp == nravad) .AND. (.NOT. lcd137))              &
          .OR. ((cma(icma)%cdtyp == ntower) .AND. (.NOT. lcd139))              &
          .OR. ((cma(icma)%cdtyp == ntowic) .AND. (.NOT. lcd159))              &
          .OR. ((cma(icma)%cdtyp == nwlidr) .AND. (.NOT. lcd187))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == nsattv) THEN
      IF       (.NOT. lsatem) THEN
!     IF (     (.NOT. lsatem)                                                  &
!         .OR. ((cma(icma)%cdtyp == nsmsg1) .AND. (mcdmsg1 <= 1))              &
!         .OR. ((cma(icma)%cdtyp == nsmsg2) .AND. (mcdmsg2 <= 1))              &
!         .OR. ((cma(icma)%cdtyp == nnoa15) .AND. (mcdno15 <= 1))              &
!         .OR. ((cma(icma)%cdtyp == nnoa16) .AND. (mcdno16 <= 1))              &
!         .OR. ((cma(icma)%cdtyp == nnoa17) .AND. (mcdno17 <= 1))              &
!         .OR. ((cma(icma)%cdtyp == nnoa18) .AND. (mcdno18 <= 1))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
!     IF (         ((cma(icma)%cdtyp == nsmsg1) .AND. (mcdmsg1 >= 1))          &
!             .OR. ((cma(icma)%cdtyp == nsmsg2) .AND. (mcdmsg2 >= 1))) THEN
!       cma(icma) %attr = 1
!     ELSEIF (     ((cma(icma)%cdtyp == nnoa15) .AND. (mcdno15 >= 1))          &
!             .OR. ((cma(icma)%cdtyp == nnoa16) .AND. (mcdno16 >= 1))          &
!             .OR. ((cma(icma)%cdtyp == nnoa17) .AND. (mcdno17 >= 1))          &
!             .OR. ((cma(icma)%cdtyp == nnoa18) .AND. (mcdno18 >= 1))) THEN
!       cma(icma) %attr = 2
!     ENDIF
    ELSEIF (cma(icma)%obtyp == ngps  ) THEN
      IF ((.NOT. lgps) .OR. (.NOT. lgpc(icma-n_gpc_offs))) THEN
!     IF (     (.NOT. lgps  )                                                  &
!         .OR. ((cma(icma)%cdtyp == ngpgfz ) .AND. (.NOT. lcd108 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpasi ) .AND. (.NOT. lcd097 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpbkg ) .AND. (.NOT. lcd098 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpgop ) .AND. (.NOT. lcd099 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpieec) .AND. (.NOT. lcd100 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpsgn ) .AND. (.NOT. lcd101 ))            &
!         .OR. ((cma(icma)%cdtyp == ngplpt ) .AND. (.NOT. lcd102 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpmet ) .AND. (.NOT. lcd103 ))            &
!         .OR. ((cma(icma)%cdtyp == ngprob ) .AND. (.NOT. lcd104 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpige ) .AND. (.NOT. lcd105 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpknmi) .AND. (.NOT. lcd106 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpnga ) .AND. (.NOT. lcd107 ))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == nscatt) THEN
      IF (     (.NOT. lscatt)                                                  &
          .OR. ((cma(icma)%cdtyp == nascat) .AND. (.NOT. lcd123))              &
          .OR. ((cma(icma)%cdtyp == nqscat) .AND. (.NOT. lcd122))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSE
      cma(icma)%exslat = cs
      cma(icma)%exnlat = cn
      cma(icma)%exwlon = cw
      cma(icma)%exelon = ce
    ENDIF
!   PRINT *,'ZZ3 ', n_cma, icma, cma(icma) %obtyp, cma(icma) %cdtyp            &
!                              , cma(icma) %name , cma(icma) %exnlat
  ENDDO

!     some code types have to be added:
!       METAR, TEMP MOBIL, PILOT MOBIL, etc.
!     following code types are not used any more:
!       kcdtyp == nabscd .AND. .NOT.lcd022    .OR.                             &
!       kcdtyp == nshred .AND. .NOT.lcd023    .OR.                             &
!       kcdtyp == ncodar .AND. .NOT.lcd141    .OR.                             &
!       kcdtyp == ncolba .AND. .NOT.lcd241    .OR.                             &
!       kcdtyp == nacar  .AND. .NOT.lcd244    .OR.                             &
!       kcdtyp == nsst   .AND. .NOT.lcd188    .OR.                             &
!       kcdtyp == nbathy .AND. .NOT.lcd063    .OR.                             &
!       kcdtyp == nrocob .AND. .NOT.lcd039    .OR.                             &
!       kcdtyp == nrocsh .AND. .NOT.lcd040    .OR.                             &

#endif

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_prep_cma_rules
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_prep_cma_rules


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for printing alerting messages
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_print_caution ( nexceed , nexceair , imaxl , fsize )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_proc_cdf" writes alerting
!   ('CAUTION') messages to a specific file (yucautn) if the ODR array size
!   is insufficient to accommodate all observations;
!   the messages include an estimation by how much certain namelist parameters
!   should be increased to render the ODR size large enough.
!
! Method:
!   If the program is run in parallel mode, the maximum over all nodes of the
!   local excess is determined and used for the estimation.
!   16.01.09
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

  INTEGER        , INTENT (INOUT)  ::  &
    nexceed  (5)     ,& ! number of reports in excess of ODR array size
    nexceair (2)        ! number of multi-level reports derived from single-
                        !           level reports but in excess of array size

  INTEGER        , INTENT (IN)     ::  &
    imaxl    (3)        ! size (report dimension) of the (3 types of) ODR

  REAL (KIND=wp) , INTENT (IN)     ::  &
    fsize    (3)        ! ratio of namelist parameter (max??o)
                        !       to ODR size (max??l)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    ierr  , nstat       ! error status

  CHARACTER (LEN=30)       ::  &
    yerr                ! error message

! Local arrays:
! ------------
!
!------------ End of header ----------------------------------------------------

  yroutine = 'obs_cdf_print_caution'

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_print_caution
!-------------------------------------------------------------------------------

! collect counters
! ----------------

#ifndef NOMPI
  IF (num_compute > 1) THEN

    CALL global_values ( nexceed, 5, 'MAX',imp_integers,icomm_cart, 0,yerr,ierr)
!   ==================

    IF (ierr == 0)                                                             &

      CALL global_values (nexceair, 2,'MAX',imp_integers,icomm_cart,0,yerr,ierr)
!     ==================

    IF (ierr /= 0)  CALL model_abort (my_cart_id, 11015, yerr, yroutine)
!                   ----------------
  ENDIF
#endif

! write messages
! --------------

  IF ((my_cart_id == 0) .AND. (     (MAXVAL(nexceed(1:4)) > 0)                 &
                               .OR. (MAXVAL(nexceair(1:2)) > 0))) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (nstat /= 0) CALL model_abort (my_cart_id, 1409, yerr, yroutine)
    IF (nexceed(2) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL SINGLE-LEVEL OBS." &
                     &," BEYOND maxsgl ",I5)' ) ntstep, nexceed(2), imaxl(2)
      nexceed(2)  =  NINT( nexceed(2) *fsize(2) ) + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST" &
                     &,I6)' ) nexceed(2)
    ENDIF
    ! write alerting messages if array size is insufficient to accommodate all
    ! multi-level reports that can be derived from single-lev (aircraft) reports
    IF (nexceair(2) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",F6.3,":",I5," NEW MULTI-LEV. AIRCR"  &
                     &,"AFT REPS BEYOND ARRAY SIZE ")' ) acthr, nexceair(2)
      nexceair(2)  =  NINT( nexceair(2) *fsize(1) ) + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY AT LEAST" &
                     &,I6)' ) nexceair(2)
    ENDIF
    IF (nexceed(1) > 0) THEN
      IF (nexceair(1) > 0)                                                     &
        WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOC MULTI-LEV. AIRCR"  &
                       &,"AFT REPORTS BEYOND ARRAY SIZE ")') ntstep, nexceair(1)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL MULTI-LEVEL OBS."  &
                     &," BEYOND maxmll ",I5)' ) ntstep, nexceed(1), imaxl(1)
      nexceed(1)  =  NINT( nexceed(1) *fsize(1) ) + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY AT LEAST" &
                     &,I6)' ) nexceed(1)
      IF (nexceed(5) > nexceed(1)) THEN
        nexceed(5)  =  NINT( nexceed(5) *fsize(1) ) + 1
        WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY ABOUT " &
                       &,I6)' ) nexceed(5)
        WRITE( nucautn,'("     OR THESE AIRCRAFT REPORTS ARE ASSIMILATED AS S" &
                       &,"INGLE-LEVEL REPORTS")' )
      ENDIF
    ENDIF
    IF (nexceed(3) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL GPS (IWV) OBS."    &
                     &," BEYOND maxgpl ",I5)' ) ntstep, nexceed(3), imaxl(3)
      nexceed(3)  =  NINT( nexceed(3) *fsize(3) ) + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxgpo BY AT LEAST" &
                     &,I6)' ) nexceed(3)
    ENDIF
!   IF (nexceed(4) > 0) THEN
!     WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," SATELLITE RETRIEVALS"   &
!                    &," BEYOND maxtvl ",I5)' ) ntstep, nexceed(4), imaxl(4)
!     WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxtvo BY AT LEAST" &
!                    &,I6)' ) nexceed(4)
!   ENDIF
    CLOSE (nucautn)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_print_caution
!-------------------------------------------------------------------------------

END SUBROUTINE obs_cdf_print_caution


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for printing warning messages
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_print_eventwarn ( mxeve , nevent , fsize )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_proc_cdf" writes summarising
!   warning messages to the file with unit number 'nustat', if the size of the
!   Observation Data Record (ODR) as specified by namelist parameters has been
!   too small to accommodate all observation reports.
!
! Method:
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

  INTEGER        , INTENT (IN)         ::       &
    mxeve               ,& ! length of first dimension of record 'nevent'
    nevent (mxeve,n_cma)   ! record containing the event counters

  REAL (KIND=wp) , INTENT (IN)         ::       &
    fsize  (3)             ! ratio of namelist parameter (max??o)
                           !       to ODR size (max??l)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    nwarn  (3)       ,& ! local max of events due to ODR size limit
    ima              ,& ! loop index over array 'cma'
    izerror             ! error status

  CHARACTER (LEN=80)       ::  yzerrmsg

! Local arrays:
! ------------
!
!------------ End of header ----------------------------------------------------

  izerror  = 0
  yzerrmsg = '   '

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_print_eventwarn
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Writing of messages which warn of insufficient array size
!-------------------------------------------------------------------------------

! sum up counters for warning messages
! ------------------------------------

  nwarn (:) = 0
  DO ima = 1 , n_cma
    IF (     (cma(ima)%obtyp == nsynop) .OR. (cma(ima)%obtyp == ndribu)        &
        .OR. (cma(ima)%obtyp == nsatob) .OR. (cma(ima)%obtyp == nscatt)) THEN
      nwarn (2)  =  nwarn(2)  +  nevent(nesodr,ima)
    ELSEIF   (cma(ima)%obtyp == nairep)                                  THEN
      nwarn (2)  =  nwarn(2)  +  nevent(nesodr,ima)
      nwarn (1)  =  nwarn(1)  +  nevent(nenoml,ima)
    ELSEIF ( (cma(ima)%obtyp == ntemp ) .OR. (cma(ima)%obtyp == npilot)) THEN
      nwarn (1)  =  nwarn(1)  +  nevent(nesodr,ima)
    ELSEIF   (cma(ima)%obtyp == ngps  )                                  THEN
      nwarn (3)  =  nwarn(3)  +  nevent(nesodr,ima)
!   ELSEIF   (cma(ima)%obtyp == nsatem)                                  THEN
!     nwarn (4)  =  nwarn(4)  +  nevent(nesodr,ima)
    ENDIF
  ENDDO
#ifndef NOMPI
  IF (num_compute > 1) THEN

    CALL global_values ( nwarn, 3, 'MAX', imp_integers                         &
                       , icomm_cart, 0, yzerrmsg, izerror )
!   ===================

  ENDIF
#endif

! print warning messages
! ----------------------

  IF (my_cart_id == 0) THEN
    nwarn (1) = NINT( REAL( nwarn(1) ) * fsize(1) + 0.49_wp )
    nwarn (2) = NINT( REAL( nwarn(2) ) * fsize(2) + 0.49_wp )
    nwarn (3) = NINT( REAL( nwarn(3) ) * fsize(3) + 0.49_wp )
!   nwarn (4) = NINT( REAL( nwarn(4) ) * fsize(4) + 0.49_wp )
!   IF (MAXVAL( nwarn ) > 0) THEN
!   IF (MAX( nwarn(1),nwarn(2),nwarn(3),nwarn(4) ) > 0) THEN
    IF (MAX( nwarn(1),nwarn(2),nwarn(3) ) > 0) THEN
      WRITE( nustat,'(''1'')' )
      WRITE( nustat,'(''0'')' )
      WRITE( nustat,'(''+     !!! CAUTION !!!!! CAUTION !!!!! ''               &
                    &,''CAUTION !!!!! CAUTION !!!!! CAUTION !!!!!'')' )
      WRITE( nustat,'(''+     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~''               &
                    &,''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')' )
    ENDIF
    IF (nwarn(1) > 0)                                                          &
        WRITE( nustat,'(''0     WARNING: array size for multi-level ''         &
                    &,''observations is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXMLO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(1)
    IF (nwarn(2) > 0)                                                          &
      WRITE( nustat,'(''0     WARNING: array size for single-level ''          &
                    &,''observations is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXSGO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(2)
    IF (nwarn(3) > 0)                                                          &
        WRITE( nustat,'(''0     WARNING: array size for GPS ''                 &
                    &,''observations is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXGPO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(3)
!   IF (nwarn(4) > 0)                                                          &
!       WRITE( nustat,'(''0     WARNING: array size for SAT ''                 &
!                   &,''retrievals   is too small'', /                         &
!                   &,''      =======  to accommodate all observations'', /    &
!                   &,''       --> Increase "MAXTVO" (namelist) by'',I5        &
!                   &,'': usually ok for local obs. array'', /                 &
!                   &,''           =================                   ''      &
!                   &,'' (possibly still insufficient for'', /                 &
!                   &,''                                               ''      &
!                   &,'' the global obs increment array!)'')' )        nwarn(4)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_print_eventwarn
!-------------------------------------------------------------------------------

END SUBROUTINE obs_cdf_print_eventwarn



#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
end module mo_cosmo_obs
