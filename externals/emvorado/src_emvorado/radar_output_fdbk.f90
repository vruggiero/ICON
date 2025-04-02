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

MODULE radar_output_fdbk

  !------------------------------------------------------------------------------
  !
  ! Description: Method of the radar forward operator EMVORADO for output of
  !              feedback files for data assimilation.
  !
  ! Method:
  !   See subroutines below
  !
  !------------------------------------------------------------------------------
  !
  ! Declarations:
  !
  ! Modules used:
  !

  USE radar_kind, ONLY :  dp
  
  USE radar_data, ONLY :    &
       miss_threshold, missthr_int, miss_value, reject_threshold

  USE radar_interface, ONLY : abort_run

  USE radar_utilities, ONLY : phi2phirot, rla2rlarot

  USE radar_data,               ONLY :          &
       nradsta_max, cmaxlen,    &
       idom, ndoms_max,                         &
       my_radar_id,       & ! rank of this PE in the radar communicator (cart+radario)
       rs_meta, cart_data, dbz_meta, &
       fdbk_meta_type

  USE radar_data_namelist, ONLY :  &
       ldebug_radsim, &
       lextdbz, lweightdbz, lfall, lonline, lsode, lsmooth, &
       lvoldata_output, &
       ydirradarout, ysubdirfof, &
       itype_supobing, &
       thin_step_azi, thin_step_range, thin_step_ele, &
       itype_metric_refl_fdbk, minval_obserr_lwc

  !------------------------------------------------------------------------------

#if defined(NUDGING) && defined(NETCDF)
  USE netcdf, ONLY :  &
       NF90_NOERR, &
       nf90_put_att, &
       nf90_redef, &
       nf90_strerror, &
       nf90_global, &
       nf90_enddef

  USE mo_fdbk, ONLY : &
       t_fdbk, &
       setup_fdbk, &
       create_fdbk, &
       print_fdbk, &
       add_verification, &
       close_fdbk, &
       cleanup_fdbk, &
       open_fdbk_read, &
       open_fdbk_write, &
       read_fdbk_meta => read_meta

#ifdef __COSMO__
  USE mo_fdbk_cosmo, ONLY : &
#else
  USE mo_fdbk_emvorado, ONLY : &
#endif
       t_acc_header, t_acc_body, t_acc_radar, t_account, &
       write_report_radar_1, write_report_radar_2

  USE mo_fdbk_tables, ONLY  : &
       VN_RREFL, VN_RADVEL, VN_HEIGHT, OT_RADAR, &
       VT_FIRSTGUESS, VT_ANALYSIS, VT_FORECAST, RC_VOR, RC_ASS, &
       ST_ACTIVE, &   ! state of an obs-veri-pair: good quality of both, will be used for DA
       ST_PASSIVE, &  !       - " - : will be set to passive due to uncertain quality of obs and model equiv.
       ST_OBS_ONLY, & !       - " - : obs good quality, but no model equivalent available
       ST_REJECTED, & !       - " - : neither obs nor model equivalent have good quality
       FL_NONE, &     ! flag/check: no flag set
       FL_OBS_ERR, &  ! flag/check: observation error too large
       FL_OPERATOR    ! flag/check: observation operator not applicable
#endif

  !==============================================================================


  !================================================================================
  !================================================================================

  IMPLICIT NONE

  !================================================================================
  !================================================================================

  !==============================================================================
  ! Interface blocks for overloaded procedures:

  !==============================================================================
  ! Public and Private Subroutines

  PRIVATE

#if defined(NUDGING) && defined(NETCDF)
  PUBLIC :: obs_write_cdf_feedobs
#endif

  !==============================================================================
  ! Module variables

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

#if defined(NUDGING) && defined(NETCDF)

!!$ Caution: data fields in order (azi, range, ele), not (range, azi, ele)!

  SUBROUTINE obs_write_cdf_feedobs ( ista, itime, tot_num_observables, cvarname, time_mod, &
     rdata_obs, rquali_obs, rerr_obs, &
     rdata_mod, rlat, rlon, rheight, rpress, &
     azivec, rvec, elevec, naz, nra, nele, &
     hmax_fdbk, &
     lastpr, lrefend, lwrite_veridata,                &
     zydate_ini,                 &
     zfdbk_meta,                       &
     nexce_rep, nexce_bdy, jerr, yerr, &
     obs_is_aliased, mod_is_aliased)

  !-------------------------------------------------------------------------------
  !
  ! Description:
  !   This module procedure of module "obs_write_cdf_feedobs" writes out the
  !   radar observations, the quality flags from the radar input files,
  !   as well as the values obtained by applying the observation operator
  !   to the model values.
  !
  ! Method:
  !   If the program is run in parallel mode, this subroutine
  !   is called for one radar station, one output timestep and one obstype (vr and z)
  !   at a time and writes one feedback file per station containing all output timesteps
  !   and obstypes.
  !   This is achieved by repeatedly opening, appending, and closing
  !   the output file during each output timestep.
  !   The size of the feedback file equals the maximum possible number of data
  !   which is determined by the number of nominal ranges, azimuts and elevations of the
  !   radar station and the number of output timesteps.
  !
  ! NOTE: Only those radar bins are actually written to the file, where both observations
  !       and simulated data are not flagged as invalid or missing (miss_value).
  !       Therefore, the files may actually contain less data (and less headers, if
  !       entire records are flagged missing) than the pre-allocated
  !       filesize would allow.
  !
  !===============================================================================

  !===============================================================================
  !
  ! Local USE statements:
  !
  !===============================================================================

  ! .. values for missing data (rmdi) and checking on missing data (rmdich)
!!$  (available in a future COSMO version)
!!$  USE data_obs_lib_cosmo, ONLY : &
!!$       rmdi, rmdich
!!$  USE data_obs_record, only : &
!!$       imdi
!!$ if available from above modules, remove the local definitions of rmdi, rmdich and imdi
!!$ below !!!

  !===============================================================================

  ! Subroutine arguments:
  ! --------------------

  INTEGER , INTENT (IN)  :: &
       ista                ,& ! station nr. in the rs_meta(ista) - meta data structure
       itime               ,& ! index of obs time
       naz                 ,& ! number of azimut bins in rdata (regular 3D data set)
       nra                 ,& ! number of range bins in rdata (regular 3D data set)
       nele                ,& ! number of elevations in rdata (regular 3D data set)
       tot_num_observables    ! total number of observables in the feedback file:
                              !  relevant for size of newly created fdbk-file bodies,
                              !  so that file will be large enough to add other observables
                              !  to the file by other calls to this routine.

  CHARACTER (LEN= * )      , INTENT (IN)  :: &
       cvarname               ! Tag for type of output field: 'vr' or 'z'

  REAL    (KIND=dp)    , INTENT (IN)  :: &
       time_mod                   ,& ! Time in minutes since model start
       hmax_fdbk                  ,& ! Maximum height MSL for data in feedback file
       rdata_obs (naz, nra, nele) ,& ! Obs. radar data to be written to feedback file (3D)
       rquali_obs(naz, nra, nele) ,& ! Quality flag of radar data, directly taken from radar obs file (3D)
       rerr_obs(naz, nra, nele)   ,& ! Observation error of radar data
       rdata_mod (naz, nra, nele) ,& ! Sim. radar data to be written to feedback file (3D)
       rlat      (naz, nra, nele) ,& ! Geo. Latitude  of bin [deg] (from simulations)
       rlon      (naz, nra, nele) ,& ! Geo. Longitude of bin [deg] (from simulations)
       rheight   (naz, nra, nele) ,& ! Height of bin [m MSL] (from simulations)
       rpress    (naz, nra, nele) ,& ! Pressure at height of bin [Pa] (from simulations)
       azivec    (naz      ) ,& ! Vector of regular azimuts
       rvec      (     nra      ) ,& ! Vector of regular range bins
       elevec    (          nele)    ! Vector of regular elevations

  LOGICAL                  , INTENT (IN)  :: &
       lwrite_veridata ,& ! Write also verification data (= sim. radar data), not only obs
       lastpr      ,& ! last time, sth is written to ncfeedobs in this model run
       lrefend        ! .t.: verification reference time set to end of model run
                      !      (for NUDGE runs only; for FORECAST (LETKF) or
                      !       NUDGECAST, veri_ref_time = start of model run)

  LOGICAL        , INTENT (IN), OPTIONAL  :: &
       obs_is_aliased ,& ! flag if obs vr is aliased or not (raw data are normally aliased)
       mod_is_aliased    ! flag if simul. vr is aliased or not

!!$  CHARACTER (LEN=10)       , INTENT (IN)  :: &
  CHARACTER (LEN=14)       , INTENT (IN)  :: &
       zydate_ini      ! reference date (e.g. start of the model run) [yyyymmddhhmmss]
  !                                  (year, month, day, hour, min, sec)

  ! Structure holding info for nc file header:
  TYPE(fdbk_meta_type), INTENT(in) :: zfdbk_meta

  INTEGER , INTENT (INOUT)  :: &
       nexce_rep   ,& ! number of reports exceeding FOF size
       nexce_bdy   ,& ! number of obs     exceeding FOF size
       jerr           ! error status variable

  ! it is assumed here that LEN >= 72 !
  CHARACTER (LEN= * )      , INTENT (INOUT)  :: &
       yerr           ! error message

  ! Local parameters:
  ! ----------------

!!$ the following missing value indicators will be available in
!!$ the future from modules "data_obs_lib_cosmo" and "data_obs_record"
  INTEGER , PARAMETER  :: &
       imdi   =  2147483647    ! missing data indicator for ODR integers (2^31 -1)

  REAL    (KIND=dp)    , PARAMETER  :: &
       rmdi   = -1.E31_dp  ,& ! commonly used missing data indicator
       rmdich = -1.E30_dp     ! commonly used check value for missing data (> rmdi)
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  REAL    (KIND=dp)    , PARAMETER  :: &
       c1000r  =     0.001_dp   ,& !
       c60     =    60.0_dp     ,& !
       c0      =     0.0_dp     ,& !
       epsy    =    1e-8_dp

  REAL(kind=dp), PARAMETER :: ln10_o10 = LOG(10.0_dp)*0.1_dp  ! for dbz -> linear


  TYPE (t_fdbk), SAVE  :: fb(nradsta_max,ndoms_max)        ! feedback meta file structure for each station
!!$  !TYPE (t_account), ALLOCATABLE :: report(:)  ! reports written to feedback file
!!$  TYPE (t_acc_header) :: rep_header(naz*nele)  ! report headers written to feedback file
!!$  TYPE (t_acc_body)   :: rep_body(naz*nele,nra)! report bodys written to feedback file
!!$  INTEGER             :: rep_len(naz*nele)     ! len of report body
!!$  INTEGER             :: rep_offset(naz*nele)  ! offset of report body
  TYPE (t_acc_header) :: rep_header(nele)  ! report headers written to feedback file
  TYPE (t_acc_body)   :: rep_body(naz*nra,nele)! report bodys written to feedback file
!!$  TYPE (t_acc_radar)  :: spec_radar(naz,nele)! specific for radar rays
  TYPE (t_acc_radar)  :: spec_radar(nele*naz)! specific for radar rays
  INTEGER             :: rep_len(nele)     ! len of report body
  INTEGER             :: rep_offset(nele)  ! offset of report body
  INTEGER             :: spec_offset(nele)  ! offset of report body

  ! Local scalars:
  ! -------------

  INTEGER  :: &
       ivarno,         & ! Type of data, VN_RREFL od. VN_RADVEL
       zspec_r_flags     ! Specification flags of data:
                         !   itype_refl,lextdbz,lsmooth,lonline,lsode,lfall,obs_is_aliased,mod_is_aliased

  INTEGER ::  i, j, k, thin_i, thin_j, thin_k, ij, ik, & ! thin_k added by EB, 05/2017
       irep, icl, iradar ,& ! loop indices (reports, characters, "radar specific" dimension)
       kk             ,& ! counter or loop index over observations for feedobs file
       nfcast         ,& ! number of forecast runs to compare with observations
       nrefmin        ,& ! diff. betw. verif. ref. time and model initial time
       nzday             ! Julian day

  INTEGER ::  &
       len_body       ,& ! number of body elements in report used in feedback file
       nw_rep         ,& ! number of reports to be put into FOF
       n_body,         & ! number of bodys to be put into FOF
       n_spec            ! number of rays to be put into FOF

  REAL    (KIND=dp)    ::  &
       zsubcategory, &
       rlon_r, rlat_r, &
       zrdata_obs, zrdata_mod, zrerr_obs

  INTEGER                    ::  ihdr_pos, ibdy_pos, inw_rep, in_body, &
       iveri_exp_id       ,& ! experiment identity number of the model verif. run
       iveri_run_class    ,& ! class (main / ass.) of the model verification run
       iveri_run_type     ,& ! run type (fcst./f.g./ana.) of the verification run
       iveri_ref_date     ,& ! reference date
       iveri_ref_time     ,& ! reference time
       iveri_start        ,& ! start of verification rel. to reference time
       iveri_end          ,& ! end   of verification rel. to reference time
       iveri_forecast_time,& ! forecast time (hhhmm)
       iiveri_ens_member  ,& ! ensemble member ( -1: deterministic)
       idomain (3)        ,& ! model domain size, input for 'create_fdbk'
       imax_rep           ,& ! max. # reports in NetCDF feedobs file, type integer
       imax_radar         ,& ! max. # obs in NetCDF feedobs file, type integer
       imax_rays          ,& ! max. # rays in NetCDF feedobs file, type integer
       imax_body          ,& ! max. # obs in NetCDF feedobs file, type integer
       varid              ,& ! netcdf variable identifier, for 'add_verification'
       imiss              ,& ! missing data indicator for integers (2^31-1)
       i_supob, j_supob, ni_supob, column_ID_supob

  INTEGER                    ::  thin_naz, thin_nra, thin_nele, thin_az_step, thin_ra_step, &
                                 imaxloc(1)

  REAL       ::  &
       resolution  (2)    ,& ! model resolution (lat, lon)
       pole        (2)    ,& ! pole (lat, lon) of rotated grid
       lower_left  (2)    ,& ! lower left  corner of model domain (lat, lon)
       upper_right (2)    ,& ! upper right corner of model domain (lat, lon)
       rmich                 ! check value for missing real data   (-1.E30)

  CHARACTER (LEN=12)  :: yveri_ref_datim      ! yyyymmddhhmm reference time
  CHARACTER (LEN=12)  :: yveri_initial_date   ! yyyymmddhhmm initial time
  CHARACTER (LEN=80)  :: yveri_descript       ! description of model verif. run
  CHARACTER (LEN=64)  :: yversion(2)          ! model and version number
!!$  CHARACTER (LEN=10)  :: zyakdat1  ! actual date in the form   yyyymmddhh
!!$  CHARACTER (LEN=22)  :: zyakdat2  ! actual date, form:   wd   dd.mm.yy  hh UTC
  CHARACTER (LEN=14)  :: zyakdat1  ! actual date in the form   yyyymmddhhmmss
  CHARACTER (LEN=28)  :: zyakdat2  ! actual date, form:   wd   dd.mm.yy  hh mm ss UTC
  CHARACTER (len=cmaxlen) :: zyfilenamefof  ! path and filename of the fof-file to produce

  INTEGER , SAVE   :: &
       n_veri               ! number of verifications, always 1 in COSMO model

  LOGICAL :: &
       zobs_is_aliased ,& ! flag if obs vr is aliased or not (raw data are normally aliased)
       zmod_is_aliased    ! flag if simul. vr is aliased or not

  ! Local arrays:
  ! ------------

  INTEGER, ALLOCATABLE  :: &
       lenv(:), rind(:), kkv(:), jerrv(:), i_ind(:), j_ind(:), k_ind(:), k_report(:)

  CHARACTER (LEN=32)         :: yzroutine

  ! Parameters for itype_metric_refl_fdbk = 2 (LWC-metric instead of dBZ) and 3 ( + LWC-dependent obserr):
  REAL(kind=dp), PARAMETER   :: &
       za_lwc = 0.004_dp,    & ! factor a of the LWC-Z-Relation LWC=a*Z^b
       zb_lwc = 0.55_dp ,    & ! factor b of the LWC-Z-Relation LWC=a*Z^b
       zddbz  = 5.0_dp         ! reference Delta_dBZ for Z-dependent obserr factor

  !------------ End of header ----------------------------------------------------


  !-------------------------------------------------------------------------------
  ! Begin Subroutine obs_write_cdf_feedobs
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !  Section 0: Initialisations
  !-------------------------------------------------------------------------------


  yzroutine(:) = ' '
  yzroutine = 'obs_write_cdf_feedobs'

  IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id



  jerr       = 0
  yerr(:)    = ' '
  nexce_rep  = 0
  nexce_bdy  = 0

  IF (lwrite_veridata) THEN
    nfcast = 1   ! will be written to "n_veri" lateron
  ELSE
    nfcast = 0
  END IF

  ! difference between model initial time and verification reference time
  !   (reference time = number of minutes from ydate_ini)
  IF (      lrefend)  nrefmin  =  NINT( zfdbk_meta%hverend *c60 )
  IF (.NOT. lrefend)  nrefmin  =  0



  IF (TRIM(cvarname) == 'vr') THEN
    ivarno = VN_RADVEL
    zsubcategory = 1
  ELSE IF (TRIM(cvarname) == 'z') THEN
    ivarno = VN_RREFL
    zsubcategory = 0
  ELSE
    ivarno = -9       ! Dummy
    zsubcategory = -9 ! Dummy
    WRITE ( yerr , '(a)') &
         TRIM(yzroutine)//' ERROR: cvarname = '//TRIM(cvarname)//'is not supported!'
    jerr = 1
    RETURN
  END IF

  !-------------------------------------------------------------------------------
  !  Section 1:  Create the netCDF feedback file
  !-------------------------------------------------------------------------------

  ! .. Calculate the verification reference date in format YYYYMMDDhhmm:
  READ( zydate_ini(1:8) ,'(I8)' )  iveri_ref_date
  READ( zydate_ini(9:12),'(I4)' )  iveri_ref_time

  WRITE( yveri_ref_datim,'(I8.8,I4.4)' )  iveri_ref_date, iveri_ref_time

  iveri_end    =  NINT( zfdbk_meta%hverend *c60 )  -  nrefmin
  iveri_start  =  NINT( zfdbk_meta%hversta *c60 )  -  nrefmin
  IF (zfdbk_meta%hversta > epsy)  iveri_start  =  MAX( iveri_start , 1-nrefmin )


  !-------------------------------------------------------------------------------
  ! .. Some maximum data dimensions:

  ! .. set to 1 for no thinning:

  SELECT CASE (itype_supobing)
  CASE (0)
    thin_az_step = thin_step_azi
    thin_ra_step = thin_step_range
  CASE (1)
    thin_az_step = 1
    thin_ra_step = 1
  CASE (2)
    thin_az_step = 1
    thin_ra_step = 1
  CASE default
    WRITE(yerr,'(a,i2,a)') 'ERROR '//TRIM(yzroutine)//': Wrong itype_supobing ', &
                           itype_supobing,'! Valid values are 0, 1, or 2!'
    CALL abort_run (my_radar_id, 80072, TRIM(yerr), 'radar_process_output.f90, '//TRIM(yzroutine))
  END SELECT

  thin_naz  = (naz-1)  / thin_az_step + 1
  thin_nra  = (nra-1)  / thin_ra_step + 1
  thin_nele = rs_meta(ista)%nel_fdbk  ! (elevation thinning is already contained in rs_meta(ista)%nel_fdbk and rs_meta(ista)%ind_ele_fdbk!)

  ! .. max. number of rays per output timestep per data category (one ray = one report):
  !imax_radar = naz * nele
  ! .. max. number of rays per output timestep per data category (one ray = one elevation):
  !imax_radar = nele
  imax_radar = thin_nele
  ! .. max. possible number of reports for all output timesteps:
  imax_rep  =  imax_radar * rs_meta(ista)%nobs_times_fdbk * tot_num_observables
  ! .. max. number of body data for all timesteps (2 datasets, vr and z):
  !imax_body =  imax_rep * nra
  !imax_body =  imax_rep * thin_nra * thin_naz

  IF (itype_supobing == 0) THEN
    imax_body =  imax_rep * thin_nra * thin_naz
  ELSE
    imax_body =  imax_rep * cart_data(ista)%ni * cart_data(ista)%nj
  ENDIF


!  imax_rays = rs_meta(ista)%nobs_times_fdbk * thin_nra * thin_naz
!  imax_rays = rs_meta(ista)%nobs_times_fdbk * nele * thin_naz * tot_num_observables
  imax_rays = imax_rep * thin_naz

  !-------------------------------------------------------------------------------
  ! .. Create feedback file if necessary:

  ! Name of the fof-file to produce:
  zyfilenamefof(:) = ' '
  IF (ysubdirfof(1:1) == '/') THEN
    zyfilenamefof    = TRIM(ysubdirfof)//TRIM(rs_meta(ista)%fdbkfile)
  ELSE
    zyfilenamefof    = TRIM(ydirradarout)//TRIM(ysubdirfof)//TRIM(rs_meta(ista)%fdbkfile)
  END IF

  IF (.NOT. rs_meta(ista)%lfdbkfile_exist) THEN

    ! determine NetCDF global attributes
    ! ----------------------------------

    resolution(1)  = zfdbk_meta%dlat
    resolution(2)  = zfdbk_meta%dlon
    idomain(1)     = zfdbk_meta%ie_tot
    idomain(2)     = zfdbk_meta%je_tot
    idomain(3)     = zfdbk_meta%ke_tot
    pole(1)        = zfdbk_meta%pollat
    pole(2)        = zfdbk_meta%pollon
    lower_left(1)  = zfdbk_meta%startlat_tot                     ! / lower left corner
    lower_left(2)  = zfdbk_meta%startlon_tot                     ! \ of model domain
    upper_right(1) = zfdbk_meta%startlat_tot + (zfdbk_meta%je_tot-1) * zfdbk_meta%dlat ! / upper right corner
    upper_right(2) = zfdbk_meta%startlon_tot + (zfdbk_meta%ie_tot-1) * zfdbk_meta%dlon ! \ of model domain
    yversion(1)    = TRIM(zfdbk_meta%yglatt_source)
    IF (yversion(1)(1:5) == 'COSMO') THEN
      yversion(2)(:) = ' '
      yversion(2)    = TRIM(ADJUSTL(yversion(1)(6:LEN_TRIM(yversion(1))-5)))
      yversion(1) = yversion(1)(1:5)//REPEAT(' ',LEN_TRIM(yversion(1))-5)
    ELSE IF (yversion(1)(1:4) == 'ICON') THEN
      yversion(2)(:) = ' '
      yversion(2)    = TRIM(ADJUSTL(yversion(1)(5:LEN_TRIM(yversion(1))-4)))
      yversion(1) = yversion(1)(1:4)//REPEAT(' ',LEN_TRIM(yversion(1))-4)
    ELSE IF (yversion(1)(1:1) == '-') THEN
      yversion(2)(:) = ' '
      yversion(1) = yversion(1)(1:1)//REPEAT(' ',LEN_TRIM(yversion(1))-1)
    END IF

    ! determine NetCDF variables related to the verification data
    ! (i.e. to the current model run)
    ! -----------------------------------------------------------

    IF ((zfdbk_meta%nvers >= 0) .AND. (zfdbk_meta%nvers < 1048576)) THEN
      iveri_exp_id     =  MOD( zfdbk_meta%nvers, 16384 )
      iveri_run_class  =       zfdbk_meta%nvers/ 16384
    ELSE
      iveri_exp_id     =  imdi
      iveri_run_class  =  imdi
    ENDIF
    yveri_initial_date  =   yveri_ref_datim
    IF (lrefend) THEN
      iveri_forecast_time =   c0
      yveri_descript(:)   =   ' '
      yveri_descript      =   'nudging assimilation'
      iveri_run_type      =   VT_ANALYSIS
    ELSE
      iveri_forecast_time =   INT( zfdbk_meta%hverend + c1000r ) * 100                     &
           + MOD( NINT( zfdbk_meta%hverend *c60 ) , 60 )
      yveri_descript(:)   =   ' '
      yveri_descript      =   'deterministic forecast'
      IF (zfdbk_meta%iveri_ens_member >= 1)                                               &
           WRITE( yveri_descript,'("ensemble forecast member ",I3)' )             &
           zfdbk_meta%iveri_ens_member
      iveri_run_type      =   VT_ANALYSIS
      IF ((iveri_run_class == RC_ASS) .OR. (zfdbk_meta%hverend < 3.01_dp))            &
           iveri_run_type    =   VT_FIRSTGUESS
    ENDIF
    iiveri_ens_member = zfdbk_meta%iveri_ens_member

    IF (ldebug_radsim) WRITE(*,'(i3,2x,a)') my_radar_id, 'Creating fdbk file: '//TRIM(zyfilenamefof)

    IF ( LEN_TRIM(rs_meta(ista)% fdbkfile) > LEN(fb(ista,idom)% nc% path) ) THEN
      WRITE( yerr ,'("ERROR '//TRIM(yzroutine)//': filename of NetCDF file too long! len = ",i6", maxlen = ",i6)' )&
           LEN_TRIM(rs_meta(ista)% fdbkfile), LEN(fb(ista,idom)% nc% path)
      jerr = 35
      RETURN
    END IF

    ! create the feedback file:
    ! =========================

    !   create and initialise fdbk tables in derived type fb:
    !--------------------------------------------------------
    CALL setup_fdbk   ( fb(ista,idom)%nc )
    !    ==========

    !   change description attribute of the observation error 'e_o':
    !---------------------------------------------------------------
    DO i=1, fb(ista,idom)%nc% nvar
      IF (TRIM(fb(ista,idom)%nc %vars(i) %name) == 'e_o') THEN
        fb(ista,idom)%nc %vars(i) %longname = 'relative (not absolute!) observational error'
        EXIT
      END IF
    END DO

    !   create the fdbk-file in overwrite mode (NF_CLOBBER), see mo_t_netcdf_file.f90:
    ! --------------------------------------------------------------------------------
    CALL create_fdbk (fb(ista,idom)          ,&! feedback file meta data
    !    ===========
         TRIM(zyfilenamefof),& ! path and filename
         TRIM(yversion(1))       ,&! model
         TRIM(yversion(2))       ,&! model version
         zfdbk_meta%yglatt_institution,&! institution
         imax_rep                ,&! d_hdr
         imax_body               ,&! d_body  (later will use d_radar instead; then set d_body = 1)
         iveri_ref_date          ,&! reference date
         iveri_ref_time          ,&! reference time
         iveri_start             ,&! start of verification
         iveri_end               ,&! end   of verification
         resolution              ,&! resolution
         idomain                 ,&! COSMO domain size
         yveri_descript          ,&! comment (for history)
         yveri_ref_datim         ,&! time    (for history)
         pole       =pole        ,&! pole of rotated grid
         lower_left =lower_left  ,&! l.l. corner of domain
         upper_right=upper_right ,&! u.r. corner of domain
         n_radar    = imax_rays  ,&! d_radar = max of radar body
         opt        = 'RADAR') ! flag optional variables

    IF (nfcast >= 1)                                                           &

         CALL add_verification (fb(ista,idom)  ,&! feedback file meta data
         TRIM(yversion(1))   ,&! model
         iveri_run_type      ,&! run type
         iveri_run_class     ,&! run class
         yveri_initial_date  ,&! initial date
         iveri_forecast_time ,&! forecast time (hhhmm)
         resolution          ,&! resolution
         idomain             ,&! domain size
         yveri_descript      ,&! description
         iiveri_ens_member   ,&! ensemble member
         iveri_exp_id        ,&! experiment id
         varid )               ! variable id
    !     ====================

    rs_meta(ista)%lfdbkfile_exist = .TRUE.

    n_veri          = nfcast
    IF (ldebug_radsim) WRITE  (*,'(A)') 'Created fdbk file: '//TRIM(zyfilenamefof)

    CALL close_fdbk (fb(ista,idom))

  ENDIF   ! .not. rs_meta(ista)%lfdbkfile_exist


  ! .. (Re-)Open the existing feedback file and retrieve global attributes:
  IF (ldebug_radsim) WRITE  (*,'(A)') '(Re)Open fdbk file: '//TRIM(zyfilenamefof)

  ! .. Open for read access to read the meta data structure:
  CALL open_fdbk_read ( fb(ista,idom), TRIM(zyfilenamefof) )

  ! .. Read the more variable meta data from the file (keep the static ones from the SAVE fb(ista,idom) above):
  CALL read_fdbk_meta  ( fb(ista,idom) )
  IF (fb(ista,idom)%nc% error /= NF90_NOERR) THEN
    jerr = fb(ista,idom)%nc% error
    WRITE( yerr ,'("ERROR in read_meta, retrieving structure fb(ista,idom) ",a) ' ) TRIM(nf90_strerror(jerr))
    CALL undwech ()
    RETURN
  END IF

  ! .. Close file ...
  CALL close_fdbk ( fb(ista,idom) )

  IF (ldebug_radsim) CALL print_fdbk ( fb(ista,idom) )

  ! .. ... and open again for write access:
  CALL open_fdbk_write ( fb(ista,idom), TRIM(zyfilenamefof) )


  IF (ldebug_radsim) WRITE  (*,'(A)') 'Feedback file opened: '//TRIM(zyfilenamefof)

  !-------------------------------------------------------------------------------
  !  Section 3:  Fill 'report' which is the input for 'write_report'
  !-------------------------------------------------------------------------------

  ! determine 'nw_rep', 'n_body' (number of new reports / observations)
  ! -------------------------------------------------------------------

  nw_rep = imax_radar

!!$! It seems to be not allowed to make record allocatable, because it has a POINTER component.
!!$! However, such derived types can be automatic arrays.
!!$!   Therefore, splitted report of type(t_account) in its components
!!$!   rep_header (type t_acc_header) und rep_body (type t_acc_body) sowie
!!!$   rep_len und rep_offset!
!!$  IF (nw_rep > 0) THEN
!!$    ALLOCATE( report( nw_rep ), STAT=jerr )
!!$    IF (jerr/=0) THEN
!!$      WRITE ( yerr ,'("ERROR in memory alloc: report(nw_rep) ",I8)') nw_rep
!!$      CALL undwech ()
!!$      RETURN
!!$    ENDIF
!!$  ENDIF

  ! initialize filling of 'report'
  ! ------------------------------


#if defined (FTRACE_RADAR) && defined(_FTRACE)
CALL ftrace_region_begin(TRIM(yzroutine)//'_preparations')
#endif

  ALLOCATE ( lenv(nw_rep), k_report(nw_rep) )
  lenv(:) = 0
  k_report(:) = 0

  ik = 0

  DO thin_k = 1, rs_meta(ista)%nel_fdbk

    ! Hint: elevation thinning is already contained in rs_meta(ista)%ind_ele_fdbk(1:rs_meta(ista)%nel_fdbk)

    ! ind_ele_fdkb(thin_k) is the index into the list of actually present elevations, not necessary equal to the nominal scan strategy.
    !  Note that some radars might have from time to time less elevations than their nominal strategy (e.g., France, Denmark)
    !  in which case this prevents writing of not present obs (missing data) to fdbk.
    IF (rs_meta(ista)%ind_ele_fdbk(thin_k) <= rs_meta(ista)%nel_present .AND. &
        rs_meta(ista)%ind_ele_fdbk(thin_k) > missthr_int) THEN

      k = rs_meta(ista)%ind_ele_present(rs_meta(ista)%ind_ele_fdbk(thin_k))
 
      ! one elevation is one report:
      ik = ik + 1
      k_report(ik) = k

      DO i= 1, nra, thin_ra_step
!NEC$ ivdep
        DO j = 1, naz, thin_az_step

        ! one elevation is one report:
          irep = ik

!!$ UB 2019: changed from .AND. or .OR. to only checking the obs; .OR. was chosen by Axel
!!$     in case of radial wind plotting from feedback files, .AND. is a time bomb because
!!$     it might lead to eliminations of all values of an elevation and a crash of the LETKF!
!!$        IF (rdata_obs(j,i,k) >= miss_threshold .OR. rdata_mod(j,i,k) >= miss_threshold) THEN
!!$        IF (rdata_obs(j,i,k) >= miss_threshold .AND. rdata_mod(j,i,k) >= miss_threshold) THEN
!!$ UB 20210219: there were problems with differing n_body among different ensemble members,
!!$        caused by dealiasing the observations. The dealiasing depends on the model winds
!!$        and may lead to differing number of non-missing obs! Solution: we just write all
!!$        all of the obs to the fof files, but set their state to SC_REJECTED below.
!!$        Because superobed observations are stored sparsely among the full radar grid,
!!$        we use a new reject_value for rejected obs (which is smaller than miss_threshold but larger than miss_value)
!!$        and discriminate by checking against a new reject_threshold
!!$        (which is smaller than reject_value but larger than miss_value).
          IF (rdata_obs(j,i,k) >= reject_threshold .AND. rheight(j,i,k) <= hmax_fdbk .AND. rheight(j,i,k) >= miss_threshold) THEN
            lenv(irep) = lenv(irep) + 1
          END IF

        END DO
      END DO

    END IF
    
  END DO

  ALLOCATE ( rind(nw_rep) )
  rind(:) = 0

  rep_len(:) = 0
  i = 0
  DO irep=1,nw_rep
    IF (lenv(irep) > 0) THEN
      i = i + 1
      rep_len(i) = lenv(irep)
      rind(i) = k_report(irep)
    END IF
  END DO
  nw_rep = i

  DEALLOCATE (lenv, k_report)

  ! update total nubmers of reports in header / body:
  ! -------------------------------------------------
  fb(ista,idom)%n_hdr = fb(ista,idom)%n_hdr + nw_rep

  n_body = 0
  n_spec = 0
  rep_offset(:) = fb(ista,idom)%n_body
  spec_offset(:) = fb(ista,idom)%n_radar
  DO irep=1, nw_rep

    rep_offset(irep) = rep_offset(irep)  + n_body
    n_body = n_body + rep_len(irep)

    spec_offset(irep) = spec_offset(irep)  + n_spec
    n_spec = n_spec + thin_naz

  END DO
  fb(ista,idom)%n_body = fb(ista,idom)%n_body + n_body
  fb(ista,idom)%n_radar = fb(ista,idom)%n_radar + n_spec

  ! fill report(irep)% header
  ! -------------------------

  ! Set radar obs and mod configuration flags:
  ! Bitwise coding for itype_refl, lextdbz, lsmooth, lonline, lsode,
  !  lfall, obs_is_aliased, mod_is_aliased:

  IF (PRESENT(obs_is_aliased)) THEN
    zobs_is_aliased = obs_is_aliased
  ELSE
    zobs_is_aliased = .TRUE.   ! which is normally the case with raw radar data ...
  END IF
  IF (PRESENT(mod_is_aliased)) THEN
    zmod_is_aliased = mod_is_aliased
  ELSE
    zmod_is_aliased = .FALSE.
  END IF

  zspec_r_flags = (dbz_meta(ista)%itype_refl-1)     ! itype_refl=1-6 and reserve up to 16 for future use
  IF (lextdbz)         zspec_r_flags = IBSET(zspec_r_flags, 5)
  IF (lsmooth)         zspec_r_flags = IBSET(zspec_r_flags, 6)
  IF (lonline)         zspec_r_flags = IBSET(zspec_r_flags, 7)
  IF (lsode)           zspec_r_flags = IBSET(zspec_r_flags, 8)
  IF (lfall)           zspec_r_flags = IBSET(zspec_r_flags, 9)
  IF (zobs_is_aliased) zspec_r_flags = IBSET(zspec_r_flags,10)
  IF (zmod_is_aliased) zspec_r_flags = IBSET(zspec_r_flags,11)

  ! (recover flags with:)
  !
  ! itype_refl = MOD(zspec_r_flags, 2**5) + 1   OR   itype_refl = IBITS(zspec_r_flags,0,5) + 1
  ! lextdbz = BTEST(zspec_r_flags, 5)
  ! lsmooth = BTEST(zspec_r_flags, 6)
  ! lonline = BTEST(zspec_r_flags, 7)
  ! lsode   = BTEST(zspec_r_flags, 8)
  ! lfall   = BTEST(zspec_r_flags, 9)
  ! zobs_is_aliased = BTEST(zspec_r_flags, 10)
  ! zmod_is_aliased = BTEST(zspec_r_flags, 11)

#if defined (FTRACE_RADAR) && defined(_FTRACE)
CALL ftrace_region_end(TRIM(yzroutine)//'_preparations')
#endif


#if defined (FTRACE_RADAR) && defined(_FTRACE)
CALL ftrace_region_begin(TRIM(yzroutine)//'_header')
#endif

  DO irep=1, nw_rep

    len_body = rep_len(irep)

    rep_header(irep)% i_body        = rep_offset(irep)+1  ! offset in body
    rep_header(irep)% l_body        = len_body  ! record length in body
    rep_header(irep)% n_level       = len_body  ! record length in body

    rep_header(irep)% i_spec        = spec_offset(irep)+1  ! offset in spec
    rep_header(irep)% l_spec        = thin_naz  ! record length in spec

    rep_header(irep)% data_category = 6  ! Initially FillValue
    rep_header(irep)% sub_category  = zsubcategory  ! 0 = reflectivity, 1 = radwind
    rep_header(irep)% varno_back    = ivarno        ! VN_RREFL od. VN_RADVEL

    rep_header(irep)% center        = 78
    rep_header(irep)% sub_center    = 255
    rep_header(irep)% obstype       = OT_RADAR
    rep_header(irep)% codetype      = -1
    rep_header(irep)% lat           = rs_meta(ista)%lat
    rep_header(irep)% lon           = rs_meta(ista)%lon
    rep_header(irep)% ident         = rs_meta(ista)%station_id

    rep_header(irep)% sun_zenit     = rmdi
    rep_header(irep)% time          = NINT(time_mod) ! Obs time min relative to reftime (take time of model start as first guess value (lrefend set .false. in LETKF))

    rep_header(irep)% time_nomi     = rep_header(irep)%time

    rep_header(irep)% time_dbase    = rep_header(irep)%time ! Provisional; should be time of data base entry section2_decoding_time

    rep_header(irep)% z_station     = rs_meta(ista)%alt_msl
    rep_header(irep)% z_modsurf     = rs_meta(ista)%msl_mod  ! bilinear interpolation
    rep_header(irep)% r_state       = ST_ACTIVE

    rep_header(irep)% vnyquist      = rs_meta(ista)%ext_nyq(rind(irep),itime)  ! Nyquist velocity incl. Dual-PRF

    rep_header(irep)% r_flags       = 0          ! ???

    rep_header(irep)% spec_r_flags  = zspec_r_flags

    rep_header(irep)% r_check       = FL_NONE

    rep_header(irep)% sta_corr      = 0  ! = 0

    rep_header(irep)% mdlsfc        = 1 ! Flag bit-wise, here bit 1 set = 1 = country

    rep_header(irep)% index_x       = rs_meta(ista)%i_nearest_mod ! nearest i-gridpoint to radar station
    rep_header(irep)% index_y       = rs_meta(ista)%j_nearest_mod ! nearest j-gridpoint to radar station

    rep_header(irep)% instype       = 1      ! selbst vergeben: 1 (WMO???)

    rep_header(irep)% phase         = 1      ! Volumenscan=1, Precip-Scan=2 (selbst vergeben)

  END DO ! irep

#if defined (FTRACE_RADAR) && defined(_FTRACE)
CALL ftrace_region_end(TRIM(yzroutine)//'_header')
#endif

!!$  ! .. Extra loop for data types which are problematic for vectorization:
!!$  DO irep=1, nw_rep
!!$
!!$    rep_header(irep)% statid(4:10)   = ' '
!!$    rep_header(irep)% statid(1:1)   = rs_meta(ista)%station_name(1:1)
!!$    rep_header(irep)% statid(2:2)   = rs_meta(ista)%station_name(2:2)
!!$    rep_header(irep)% statid(3:3)   = rs_meta(ista)%station_name(3:3)
!!$
!!$  END DO


  ! Allocate report(irep)% body
  ! -----------------------

  ! .. Very inefficient on the NEC !
!!$  DO irep=1, nw_rep
!!$
!!$ !!! This is problematic, because report%body is itself a structure with
!!$ !!! allocatable pointer component "veri_data", so ALLOCATE should fail here ...
!!$ !!!   removed allocatable attribute from veri_data!!!
!!$    ALLOCATE( report(irep)% body(report(irep)% header% l_body) , STAT=jerr )
!!$    IF (jerr/=0) THEN
!!$      WRITE( yerr ,'("ERROR in memory alloc: report(irep)%body(len_body) "      &
!!$           &,2I8)' ) irep, report(irep)% header% l_body
!!$      CALL undwech ()
!!$      RETURN
!!$    ENDIF
!!$    IF (report(irep)% offset + report(irep)% header% l_body > imax_body) THEN
!!$      ! this should never happen as 'max_body' is already checked further above
!!$      WRITE( yerr , '("ERROR: MAX_BODY for feedobs file too small ",3I8,I6)' )    &
!!$           report(irep)% offset, report(irep)% header% l_body, imax_body, irep
!!$      jerr = report(irep)% offset + report(irep)% header% l_body - imax_body
!!$      CALL undwech ()
!!$      RETURN
!!$    ENDIF
!!$
!!$    IF (n_veri >= 1) THEN
!!$      DO kk = 1, report(irep)% header% l_body
!!$        ALLOCATE ( report(irep)% body(kk)% veri_data(n_veri) , STAT=jerr )
!!$        IF (jerr/=0) THEN
!!$          WRITE( yerr ,'(a,I6,I8,I3)') &
!!$               'ERROR in alloc: report(irep)%body(kk)%veri_data(n_veri) ',irep, kk, n_veri
!!$          CALL undwech ()
!!$          RETURN
!!$        ENDIF
!!$      END DO
!!$    END IF
!!$
!!$  END DO  ! irep

  ! fill report(irep)% body
  ! -----------------------

#if defined (FTRACE_RADAR) && defined(_FTRACE)
CALL ftrace_region_begin(TRIM(yzroutine)//'_body')
#endif

  ALLOCATE ( kkv(nw_rep), jerrv(nw_rep), stat=jerr)
  IF (jerr/=0) THEN
    WRITE( yerr ,'(a,I0)') &
         'ERROR in alloc: kkv(nw_rep), jerrv(nw_rep), nw_rep = ', nw_rep
    CALL undwech ()
    RETURN
  ENDIF

  ALLOCATE ( i_ind(thin_naz*thin_nra), &
             j_ind(thin_naz*thin_nra), &
             k_ind(thin_naz*thin_nra), stat=jerr)
  IF (jerr/=0) THEN
    WRITE( yerr ,'(a)') &
         'ERROR in alloc: i_ind, j_ind, k_ind'
    CALL undwech ()
    RETURN
  ENDIF


  kkv(:) = 0
!!$  spec_radar(:,:)% nbody = 0
  spec_radar(:)% nbody = 0
  jerrv(:) = 0

!$omp parallel do private(i,j,k,kk,ij,thin_i,thin_j,irep,i_ind,j_ind,k_ind), &
!$omp&            private(rlon_r,rlat_r,i_supob,j_supob,column_ID_supob,iradar,zrdata_obs,zrdata_mod,zrerr_obs)
  DO irep=1, nw_rep

    k = rind(irep)

    kk = 0

!!$    DO thin_j = 1, thin_naz
!!$
!!$      j = (thin_j-1)*thin_az_step + 1
!!$
!!$      DO thin_i = 1, thin_nra
!!$
!!$        i = (thin_i-1)*thin_ra_step + 1

    DO ij = 1, thin_naz*thin_nra

      thin_i = MODULO(ij-1,thin_nra) + 1
      i = (thin_i-1)*thin_ra_step + 1
      thin_j = (ij-1)/thin_nra + 1
      j = (thin_j-1)*thin_az_step + 1

!!$ UB: changed from .AND. or .OR. to only checking the obs; .OR. was chosen by Axel
!!$     in case of radial wind plotting from feedback files, .AND. is a time bomb because
!!$     it might lead to eliminations of all values of an elevation and a crash of the LETKF!
!!$      IF (rdata_obs(j,i,k) >= miss_threshold .OR. rdata_mod(j,i,k) >= miss_threshold) THEN
!!$      IF (rdata_obs(j,i,k) >= miss_threshold .AND. rdata_mod(j,i,k) >= miss_threshold) THEN
!!$ UB 20210219: there were problems with differing n_body among different ensemble members,
!!$        caused by dealiasing the observations. The dealiasing depends on the model winds
!!$        and may lead to differing number of non-missing obs! Solution: we just write all
!!$        all of the obs to the fof files, but set their state to SC_REJECTED below.
!!$        Because superobed observations are stored sparsely among the full radar grid,
!!$        we use a new reject_value for rejected obs (which is smaller than miss_threshold but larger than miss_value)
!!$        and discriminate by checking against a new reject_threshold
!!$        (which is smaller than reject_value but larger than miss_value).
      IF (rdata_obs(j,i,k) >= reject_threshold .AND. rheight(j,i,k) <= hmax_fdbk .AND. rheight(j,i,k) >= miss_threshold) THEN

        kk = kk + 1
        IF (kk > rep_header(irep)% l_body) THEN
          ! This should never happen!
          jerrv (irep) = irep
        ELSE
          i_ind(kk) = i
          j_ind(kk) = j
          k_ind(kk) = k
        END IF

      END IF

!!$      END DO

    END DO

    kkv(irep) = kk

    ! prepare unique vertical column ID for superobing points:
    IF (itype_supobing /= 0) THEN
      ni_supob = CEILING((cart_data(ista)%endlon - cart_data(ista)%startlon) / cart_data(ista)%dlat) + 1
    END IF

!CDIR VOVERTAKE,VOB3
    DO kk = 1, kkv(irep)

      i = i_ind(kk)
      j = j_ind(kk)
      k = k_ind(kk)

      ! compute unique vertical column ID for superobing points:
      IF (itype_supobing /= 0) THEN
        rlon_r = rla2rlarot(rlat(j,i,k), rlon(j,i,k), &
             cart_data(ista)%pollat, cart_data(ista)%pollon, cart_data(ista)%polgam)
        rlat_r = phi2phirot(rlat(j,i,k), rlon(j,i,k), &
             cart_data(ista)%pollat, cart_data(ista)%pollon)
        i_supob = NINT((rlon_r - cart_data(ista)%startlon) / cart_data(ista)%dlon) + 1
        j_supob = NINT((rlat_r - cart_data(ista)%startlat) / cart_data(ista)%dlat) + 1
        column_ID_supob = i_supob + (j_supob-1)*ni_supob
      ELSE
        column_ID_supob = -1
      END IF

      ! choose the metric / transformation of reflectivity and the obs error factor for the feedback files:
      IF (ivarno == VN_RREFL) THEN
        SELECT CASE (itype_metric_refl_fdbk)
        CASE (1)
          ! "normal" logarithmic reflectivity in dBZ:
          IF (rdata_obs(j,i,k) >= miss_threshold) THEN
            zrdata_obs = rdata_obs(j,i,k)
          ELSE
            zrdata_obs = miss_value  ! reset from reject_value to miss_value
          END IF
          zrdata_mod = rdata_mod(j,i,k)
          zrerr_obs  = rerr_obs(j,i,k)
        CASE (2)
          ! convert dBZ to an effective LWC: LWC = a*Z^b
          IF (rdata_obs(j,i,k) >= miss_threshold) THEN
            zrdata_obs = za_lwc*EXP(zb_lwc*rdata_obs(j,i,k)*ln10_o10)
          ELSE
            zrdata_obs = miss_value
          END IF
          zrerr_obs  = rerr_obs(j,i,k)
          IF (rdata_mod(j,i,k) >= miss_threshold) THEN
            zrdata_mod = za_lwc*EXP(zb_lwc*rdata_mod(j,i,k)*ln10_o10)
          ELSE
            zrdata_mod = miss_value
          END IF
        CASE (3)
          ! convert dBZ to an effective LWC: LWC = a*Z^b and define obs-dependent obserr factor:
          IF (rdata_obs(j,i,k) >= miss_threshold) THEN
            zrdata_obs = za_lwc*EXP(zb_lwc*rdata_obs(j,i,k)*ln10_o10)
            ! obs-dependent obserr factor: when multpilied by a certain Delta_dBZ (one-sided standard dev.) in the LETKF, this factor
            !  reflects a constant Delta_dBZ for larger values of dBZ_obs and a constant factor minval_obserr_lwc
            !  for smaller values: o_e = zbase_obserr + 0.5 * a * ( ( Z + 5dB)^b - (Z - 5dB)^b ) / 5dB
            zrerr_obs  = minval_obserr_lwc + 0.5_dp * za_lwc * &
                 ( EXP(zb_lwc*(rdata_obs(j,i,k)+zddbz)*ln10_o10) - &
                   EXP(zb_lwc*(rdata_obs(j,i,k)-zddbz)*ln10_o10) ) / zddbz
          ELSE
            zrdata_obs = miss_value
            zrerr_obs  = rerr_obs(j,i,k)
          END IF
          IF (rdata_mod(j,i,k) >= miss_threshold) THEN
            zrdata_mod = za_lwc*EXP(zb_lwc*rdata_mod(j,i,k)*ln10_o10)
          ELSE
            zrdata_mod = rdata_mod(j,i,k)
          END IF
        END SELECT
      ELSE
        IF (rdata_obs(j,i,k) >= miss_threshold) THEN
          zrdata_obs = rdata_obs(j,i,k)
        ELSE
          zrdata_obs = miss_value  ! reset from reject_value to miss_value
        END IF
        zrdata_mod = rdata_mod(j,i,k)
        zrerr_obs  = rerr_obs(j,i,k)
      END IF

      rep_body(kk,irep)% varno      = ivarno          ! VN_RREFL od. VN_RADVEL
      rep_body(kk,irep)% obs        = zrdata_obs
      rep_body(kk,irep)% bcor       = 0.0             ! bias corr = 0.0
      rep_body(kk,irep)% e_o        = zrerr_obs       ! Observation error (factor which is multiplied to the one from the LETKF-namelist)
      rep_body(kk,irep)% level      = rheight(j,i,k)  ! Height of radar-bins in m MSL
      rep_body(kk,irep)% dlon       = rlon(j,i,k)     ! geographic longitude, as estimated from 4/3-earth beam propagation model
      rep_body(kk,irep)% dlat       = rlat(j,i,k)     ! geographic longitude, as estimated from 4/3-earth beam propagation model
      rep_body(kk,irep)% plevel     = rpress (j,i,k)  ! corresponding pressure in Pa
      rep_body(kk,irep)% accuracy   = rep_body(kk,irep)% e_o ! Observation error 0.5 dB bzw. 0.2 m/s
      rep_body(kk,irep)% level_typ  = VN_HEIGHT
      rep_body(kk,irep)% level_sig  = column_ID_supob  ! unique ID of the vertical superobservation column
      IF (zrdata_mod >= miss_threshold .AND. zrdata_obs >= miss_threshold) THEN
        rep_body(kk,irep)% state    = ST_ACTIVE ! both model and obs have good quality, use in DA
        rep_body(kk,irep)% flags    = 0         ! 0: no check failed yet
        rep_body(kk,irep)% check    = FL_NONE   ! flag to characterize that no check in EMVORADO failed
      ELSE IF (zrdata_obs >= miss_threshold) THEN
        rep_body(kk,irep)% state    = ST_OBS_ONLY ! only obs has good quality, model value not available
        rep_body(kk,irep)% flags    = IBSET(0, FL_OPERATOR) ! set FL_OPERATOR bit to 1: observation operator not applicable
        rep_body(kk,irep)% check    = 1
      ELSE
        rep_body(kk,irep)% state    = ST_REJECTED ! neither obs nor model value has good quality, so reject obs
        rep_body(kk,irep)% flags    = IBSET(0, FL_OPERATOR) ! set FL_OPERATOR bit to 1: observation operator not applicable
        rep_body(kk,irep)% flags    = IBSET(rep_body(kk,irep)% flags, FL_OBS_ERR)  ! set FL_OBS_ERR bit to 1: obs error too large
        rep_body(kk,irep)% check    = 1
      END IF
      IF (ivarno == VN_RADVEL .AND. rs_meta(ista)%ext_nyq(k,itime) < rs_meta(ista)%vnyq_min_for_vr_active_fdbk) THEN
        ! This radial wind obs has a too small Nyquist velocity and is potentially problematic, so set to PASSIVE
        rep_body(kk,irep)% state    = ST_PASSIVE
        rep_body(kk,irep)% flags    = IBSET(0, FL_OBS_ERR) ! set FL_OBS_ERR bit to 1: obs error too large
        rep_body(kk,irep)% check    = 1
      END IF
      rep_body(kk,irep)% qual       = NINT(MIN(rquali_obs(j,i,k), 1.0_dp)) * 100
!!$            rep_body(kk,irep)% spec_index = i+(j-1)*nra+(k-1)*nra*naz  ! here: "total range-azi-ele index"
!!$            rep_body(kk,irep)% spec_index = i+(j-1)*nra                ! here: "total range-azi-index"
      rep_body(kk,irep)% spec_index = i                      ! here: range index (i)

      ! fill rep_body(kk,irep)% veri_data  (if n_veri >= 1)
      ! --------------------------------------
      IF (n_veri >= 1) THEN

        ! .. WORKS ONLY FOR N_VERI = 1 !!!

        ! Set missing values in a way that they are recognized as missing values in LETKF:
        IF (zrdata_mod >= miss_threshold) THEN
          rep_body(kk,irep)% veri_data1 = zrdata_mod
        ELSE
          rep_body(kk,irep)% veri_data1 = rmdi
        END IF

      END IF

!!$ The following is not vectorized, so that the entire loop is only partially vectorized.
!!$ However, the performance is not so bad, and up to now no better way
!!$ of doing this could be found.
      thin_j = (j-1)/thin_ra_step + 1
      iradar = thin_j + (irep-1)*thin_naz
      spec_radar(iradar)% nbody  = spec_radar(iradar)% nbody + 1

    END DO

  END DO
!$omp end parallel do

  IF ( ANY(jerrv(:) > 0) ) THEN
    ! this should never happen because of the way report(irep)%len and rep_header(irep)% l_body have been calculated above...
    imaxloc = MAXLOC(array=jerrv(:))
    WRITE( yerr ,'("ERROR: rep_header(irep)% l_body for feedobs file too small ",3I8,I6)' )    &
         rep_offset(imaxloc(1)), rep_header(imaxloc(1))% l_body, imax_body, imaxloc(1)
    jerr = rep_offset(imaxloc(1)) + rep_header(imaxloc(1))% l_body + imax_body
    CALL undwech ()
    RETURN
  END IF

!!$  DO thin_j = 1, thin_naz
!!$
!!$    j = (thin_j-1)*thin_az_step + 1
!!$
!!$!CDIR NODEP
!!$    DO irep = 1, nw_rep
!!$
!!$   !!! THIS IS THE WRONG INDEX ORDER ANYWAY. SHOULD BE spec_radar(thin_j,irep) instead!
!!$
!!$      spec_radar(irep,thin_j)% azimuth     = rs_meta(ista)%az_start + (j-1)*rs_meta(ista)%az_inc
!!$      spec_radar(irep,thin_j)% elevation   = rs_meta(ista)%el_arr(irep)
!!$      spec_radar(irep,thin_j)% range_start = rs_meta(ista)%ra_inc  ! yes, it is really the same as ra_inc!
!!$      spec_radar(irep,thin_j)% drange      = rs_meta(ista)%ra_inc
!!$      spec_radar(irep,thin_j)% nrange      = rs_meta(ista)%nra
!!$
!!$    END DO
!!$  END DO

!CDIR NODEP
!$omp parallel do private(iradar,thin_j,irep,j,k)
  DO iradar = 1, thin_naz*nw_rep

    thin_j = MOD(iradar-1,thin_naz) + 1
    irep   = (iradar-1) / thin_naz  + 1
    j = (thin_j-1)*thin_az_step + 1
    k = rind(irep)

    spec_radar(iradar)% azimuth     = rs_meta(ista)%az_start + (j-1)*rs_meta(ista)%az_inc
    spec_radar(iradar)% elevation   = rs_meta(ista)%el_arr(k)
    spec_radar(iradar)% range_start = rs_meta(ista)%ra_inc  ! yes, it is really the same as ra_inc!
    spec_radar(iradar)% drange      = rs_meta(ista)%ra_inc
    spec_radar(iradar)% nrange      = rs_meta(ista)%nra


  END DO
!$omp end parallel do

  DEALLOCATE(rind,kkv,jerrv)
  DEALLOCATE(i_ind,j_ind,k_ind)

#if defined (FTRACE_RADAR) && defined(_FTRACE)
CALL ftrace_region_end(TRIM(yzroutine)//'_body')
#endif

  !-------------------------------------------------------------------------------
  !  Section 4:  Write all obs reports (in 'report') into NetCDF feedobs file
  !-------------------------------------------------------------------------------

#if defined (FTRACE_RADAR) && defined(_FTRACE)
CALL ftrace_region_begin(TRIM(yzroutine)//'_writenetcdf')
#endif

  IF (nw_rep > 0) THEN

    ihdr_pos  =  fb(ista,idom)%n_hdr  - nw_rep + 1
    ibdy_pos  =  fb(ista,idom)%n_body - n_body + 1
    inw_rep   =  nw_rep
    in_body   =  n_body
    imiss     =  imdi
    rmich     =  rmdich

    CALL write_report_radar_1 ( fb(ista,idom), rep_header(1:nw_rep), rep_body(:,1:nw_rep), &
                            spec_radar(1:nw_rep*thin_naz), &
                            inw_rep, SIZE(rep_body,2), thin_naz, &
                            in_body, ihdr_pos, ibdy_pos, &
                            imiss, rmich, jerr, yerr, &
                            yzstatid=rs_meta(ista)%station_name(1:3))
    !   =================

    IF (jerr/=0) THEN
      yerr = 'ERROR from write_report: '//yerr
      CALL undwech ()
      RETURN
    END IF

  ENDIF  !  nw_rep

  ! print statistics
  ! ----------------------------------------------
  IF ( lastpr .AND. ldebug_radsim ) THEN
    WRITE (*,'(A,/,T10,A)')     '    NetCDF FEEDOBS FILE created: ' , TRIM(zyfilenamefof)
    WRITE (*,'(A,I10,A,I8)') '    Total number of reports in fdbk file:'        &
                           , fb(ista,idom)% n_hdr,  " ; max allowed : ", imax_rep
    WRITE (*,'(A,I10,A,I8)') '    Total number of bodys in fdbk file: '         &
                           , fb(ista,idom)% n_body, " ; max allowed : ", imax_body
  END IF


  ! set some global attributes (e.g. n_hdr and n_body) to the actual values and close FOF
  ! -------------------------------------------------------------------------------------
  jerr          = nf90_redef ( fb(ista,idom)% nc% ncid )
  IF (jerr /= nf90_noerr) THEN
    WRITE( yerr,'("ERROR in nf90_redef in '//TRIM(yzroutine)//'",a)' ) trim(nf90_strerror(jerr))
    CALL undwech ()
    RETURN
  END IF
  fb(ista,idom)% nc% error = nf90_put_att ( fb(ista,idom)%nc% ncid, NF90_GLOBAL, 'n_hdr', fb(ista,idom)% n_hdr )
  IF (fb(ista,idom)% nc% error /= nf90_noerr) THEN
    WRITE( yerr,'("ERROR in nf90_put_att in '//TRIM(yzroutine)//'",a)' ) trim(nf90_strerror(fb(ista,idom)% nc% error))
    CALL undwech ()
    RETURN
  ENDIF
  fb(ista,idom)% nc% error = nf90_put_att ( fb(ista,idom)%nc% ncid, NF90_GLOBAL, 'n_body', fb(ista,idom)% n_body )
  IF (fb(ista,idom)% nc% error /= nf90_noerr) THEN
    WRITE( yerr,'("ERROR in nf90_put_att in '//TRIM(yzroutine)//'",a)' ) trim(nf90_strerror(fb(ista,idom)% nc% error))
    CALL undwech ()
    RETURN
  ENDIF
  fb(ista,idom)% nc% error = nf90_put_att ( fb(ista,idom)%nc% ncid, NF90_GLOBAL,'n_radar', fb(ista,idom)% n_radar )
  IF (fb(ista,idom)% nc% error /= nf90_noerr) THEN
    WRITE( yerr,'("ERROR in nf90_put_att in '//TRIM(yzroutine)//'",a)' ) trim(nf90_strerror(fb(ista,idom)% nc% error))
    CALL undwech ()
    RETURN
  ENDIF
  jerr          = nf90_enddef ( fb(ista,idom)% nc% ncid )
  IF (jerr /= nf90_noerr) THEN
    WRITE( yerr,'("ERROR in nf90_enddef in '//TRIM(yzroutine)//'",a)' ) trim(nf90_strerror(jerr))
    CALL undwech ()
    RETURN
  ENDIF

#if defined (FTRACE_RADAR) && defined(_FTRACE)
CALL ftrace_region_end(TRIM(yzroutine)//'_writenetcdf')
#endif

  CALL undwech ()

  IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  !-------------------------------------------------------------------------------
  ! End of module procedure obs_write_cdf_feedobs
  !-------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE undwech
    CALL close_fdbk   ( fb(ista,idom) )
    !   ===============
    IF (lastpr) CALL cleanup_fdbk ( fb(ista,idom) )
    !   =================
  END SUBROUTINE undwech


END SUBROUTINE obs_write_cdf_feedobs

#endif     /* defined(NUDGING) && defined(NETCDF) */

END MODULE radar_output_fdbk
