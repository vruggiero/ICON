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


MODULE radar_output_utils

!------------------------------------------------------------------------------
!
! Description: Utilities of the radar forward operator EMVORADO for processing
!              of the various output methods, data and formats.
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
  
  USE radar_utilities, ONLY : get_free_funit

  USE radar_data,               ONLY :          &
       cmaxlen,    &
       radar_meta_type,         &
       idom,                          &
       my_radar_id,       & ! rank of this PE in the radar communicator (cart+radario)
       my_radar_id_dom,   & ! rank of this PE in the radar communicator (cart+radario_dom)
       my_radario_id_dom, & ! rank of this PE in the (asynchroneous) radario communicator ic
       num_radar,         & ! number of radar PEs (num_compute + num_radario)
       num_radar_dom,     & ! number of radar PEs (num_compute + num_radario_dom) per radar-active model domain
       num_radario,       & ! number of radar-IO PEs
       num_radario_dom,   & ! number of radar-IO PEs per radar-active model domain
       radario_master_dom,& ! root-PEs of the radario group for each active radar domain (in the radar_dom comm., not radario_dom-comm.!!!)
       rs_meta, dbz_meta, &
       nradsta

  USE radar_interface, ONLY : &
       abort_run,          &
       get_model_time_sec, &
       get_model_time_ddhhmmss, &
       get_obs_time_tolerance, &
       get_datetime_act, &
       get_datetime_ini, &
       it_is_time_for_radar, &
       grib2_add_modelspec_info

  USE radar_data_namelist, ONLY :  &
       ldebug_radsim, loutradwind, &
       loutdbz, loutpolstd, loutpolall, lextdbz, lout_geom, &
       lweightdbz, lfall, lonline, lsode, lsmooth, lreadmeta_from_netcdf, &
       lmds_z, lmds_vr, &
       lvoldata_output, &
       itype_obserr_vr, &
       ydirradarout

  USE radar_data_io, ONLY : radgeomoutputunit, radwindoutputunit, &
       &                    radwindobsoutputunit, radrefloutputunit, &
       &                    radreflobsoutputunit, &
       &                    zdroutputunit, zdrobsoutputunit, &
       &                    rhvoutputunit, rhvobsoutputunit, &
       &                    kdpoutputunit, kdpobsoutputunit, &
       &                    ahoutputunit, adpoutputunit, &
       &                    ldroutputunit, ldrobsoutputunit

#ifdef GRIBAPI
  USE eccodes, ONLY : GRIB_SUCCESS, codes_get_error_string
#endif

#ifdef NETCDF
  USE netcdf, ONLY :  &
       NF90_NOERR, &
       NF90_strerror
#endif
  
  !------------------------------------------------------------------------------

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

  PUBLIC ::  opendiagfiles, closediagfiles, get_fileprefix_ascii_output, &
             control_output, write_ready_radar, to_lower_or_upper, &
             check_codes_err, check_nc,  get_next_key_from_pos, &
             replace_substr_with_value

  PUBLIC :: eps_ascii

  !==============================================================================
  ! Module variables

  !.. Small epsilon-threshold:
  REAL (KIND=dp), PARAMETER :: eps_ascii = 1e-30_dp  ! |values| < eps_ascii are set to 0.0_dp in ASCII output


  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  SUBROUTINE get_filename_diagfiles (varname, rsm, filename)
    CHARACTER(len=*),      INTENT(in)    :: varname
    type(radar_meta_type), INTENT(in)    :: rsm
    CHARACTER(len=*),      INTENT(inout) :: filename

    filename(:) = ' '
    WRITE (filename, '(a,"_ID-",i6.6,"_",a)') TRIM(ydirradarout)//TRIM(varname), rsm%station_id, TRIM(rsm%scanname)

  END SUBROUTINE get_filename_diagfiles

  SUBROUTINE opendiagfiles ( openstatus, openpos )

    !------------------------------------------------------------------------------
    !
    ! Description: Opens file units for diagnostic output files
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    CHARACTER (len=*), INTENT(in) :: openstatus   ! e.g., "replace" or "old" or "unknown"
    CHARACTER (len=*), INTENT(in) :: openpos      ! e.g., "rewind" or "append"

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    ! name of control output file for radar beam geometry (height, local elevation)
    CHARACTER (LEN=cmaxlen)                   :: outfile_hr
    ! name of control output file for simul. radar radial wind
    CHARACTER (LEN=cmaxlen)                   :: outfile_vr
    ! name of control output file for obs. radar radial wind
    CHARACTER (LEN=cmaxlen)                   :: outfile_vr_obs
    ! name of control output file for simul. radar reflectivity
    CHARACTER (LEN=cmaxlen)                   :: outfile_zr
    ! name of control output file for obs. radar reflectivity
    CHARACTER (LEN=cmaxlen)                   :: outfile_zr_obs
    ! name of control output file for radar extinction
    CHARACTER (LEN=cmaxlen)                   :: outfile_ze
    ! name of control output file for radar differential extinction
    CHARACTER (LEN=cmaxlen)                   :: outfile_zed
    ! name of control output file for simul. ZDR
    CHARACTER (LEN=cmaxlen)                   :: outfile_zdr
    ! name of control output file for obs. ZDR
    CHARACTER (LEN=cmaxlen)                   :: outfile_zdr_obs
    ! name of control output file for simul. KDP
    CHARACTER (LEN=cmaxlen)                   :: outfile_kdp
    ! name of control output file for obs. KDP
    CHARACTER (LEN=cmaxlen)                   :: outfile_kdp_obs
    ! name of control output file for simul. RhoHV
    CHARACTER (LEN=cmaxlen)                   :: outfile_rhv
    ! name of control output file for obs. RhoHV
    CHARACTER (LEN=cmaxlen)                   :: outfile_rhv_obs
    ! name of control output file for simul. LDR
    CHARACTER (LEN=cmaxlen)                   :: outfile_ldr
    ! name of control output file for obs. LDR
    CHARACTER (LEN=cmaxlen)                   :: outfile_ldr_obs

    CHARACTER (LEN=32)   :: yzroutine
    INTEGER              :: ista, ista_par, ipe_rad, ios, num_io_loc

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE opendiagfiles
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'opendiagfiles'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! .. Open the file units for the control output for each data set of
    !    each radar station:
    IF (.FALSE.) THEN
!!$ Left in the code, so that in future, one might go
!!$ back to output on node 0 only for whatever reason:
      ! .. Here, a parallel output is not implemented yet, so
      !    all the files for the control output for each radar
      !    have to be opened on PE 0
      IF (my_radario_id_dom(idom) == 0) THEN
        DO ista = 1, nradsta

          IF (lout_geom) THEN
            CALL get_filename_diagfiles ('YURADGEOM', rs_meta(ista), outfile_hr)
            CALL get_free_funit(radgeomoutputunit(ista))
            OPEN(radgeomoutputunit(ista),file=TRIM(outfile_hr),iostat=ios, &
                 status=TRIM(openstatus), position=TRIM(openpos), action='write')
            CALL openerrmsg ( ios, TRIM(outfile_hr) )
          END IF

          IF (loutradwind) THEN
            CALL get_filename_diagfiles ('YURADWIND', rs_meta(ista), outfile_vr)
            CALL get_free_funit(radwindoutputunit(ista))
            OPEN(radwindoutputunit(ista),file=TRIM(outfile_vr),iostat=ios, &
                 status=TRIM(openstatus), position=TRIM(openpos), action='write')
            CALL openerrmsg ( ios, TRIM(outfile_vr) )
            IF (lreadmeta_from_netcdf) THEN
              CALL get_filename_diagfiles ('YURADWINDOBS', rs_meta(ista), outfile_vr_obs)
              CALL get_free_funit(radwindobsoutputunit(ista))
              OPEN(radwindobsoutputunit(ista),file=TRIM(outfile_vr_obs),iostat=ios, &
                   status=TRIM(openstatus), position=TRIM(openpos), action='write')
              CALL openerrmsg ( ios, TRIM(outfile_vr_obs) )
            END IF
          END IF

          IF (loutdbz .OR. (loutradwind .AND. (lmds_vr .OR. (lsmooth .AND. lweightdbz)))) THEN
            CALL get_filename_diagfiles ('YURADREFL', rs_meta(ista), outfile_zr)
            CALL get_free_funit(radrefloutputunit(ista))
            OPEN(radrefloutputunit(ista),file=TRIM(outfile_zr),iostat=ios, &
                 status=TRIM(openstatus), position=TRIM(openpos), action='write')
            CALL openerrmsg ( ios, TRIM(outfile_zr) )
          END IF

          IF (loutdbz) THEN
            IF ( (loutpolstd .OR. loutpolall) .AND. &
                 (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4) ) THEN
              CALL get_filename_diagfiles ('YUDIFFREFL', rs_meta(ista), outfile_zdr)
              CALL get_free_funit(zdroutputunit(ista))
              OPEN(zdroutputunit(ista),file=TRIM(outfile_zdr),iostat=ios, &
                   status=TRIM(openstatus), position=TRIM(openpos), action='write')
              CALL openerrmsg ( ios, TRIM(outfile_zdr) )
              IF (lreadmeta_from_netcdf) THEN
                CALL get_filename_diagfiles ('YUDIFFREFOBS', rs_meta(ista), outfile_zdr_obs)
                CALL get_free_funit(zdrobsoutputunit(ista))
                OPEN(zdrobsoutputunit(ista),file=TRIM(outfile_zdr_obs),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_zdr_obs) )
              END IF

              CALL get_filename_diagfiles ('YUCORRCOEFF', rs_meta(ista), outfile_rhv)
              CALL get_free_funit(rhvoutputunit(ista))
              OPEN(rhvoutputunit(ista),file=TRIM(outfile_rhv),iostat=ios, &
                   status=TRIM(openstatus), position=TRIM(openpos), action='write')
              CALL openerrmsg ( ios, TRIM(outfile_rhv) )
              IF (lreadmeta_from_netcdf) THEN
                CALL get_filename_diagfiles ('YUCORRCOEFFOBS', rs_meta(ista), outfile_rhv_obs)
                CALL get_free_funit(rhvobsoutputunit(ista))
                OPEN(rhvobsoutputunit(ista),file=TRIM(outfile_rhv_obs),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_rhv_obs) )
              END IF

              CALL get_filename_diagfiles ('YUSPECDIFFPHASE', rs_meta(ista), outfile_kdp)
              CALL get_free_funit(kdpoutputunit(ista))
              OPEN(kdpoutputunit(ista),file=TRIM(outfile_kdp),iostat=ios, &
                   status=TRIM(openstatus), position=TRIM(openpos), action='write')
              CALL openerrmsg ( ios, TRIM(outfile_kdp) )
              IF (lreadmeta_from_netcdf) THEN
                CALL get_filename_diagfiles ('YUSPECDIFFPHASEOBS', rs_meta(ista), outfile_kdp_obs)
                CALL get_free_funit(kdpobsoutputunit(ista))
                OPEN(kdpobsoutputunit(ista),file=TRIM(outfile_kdp_obs),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_kdp_obs) )
              END IF

              IF (loutpolall) THEN
                CALL get_filename_diagfiles ('YULINDEPOLRAT', rs_meta(ista), outfile_ldr)
                CALL get_free_funit(ldroutputunit(ista))
                OPEN(ldroutputunit(ista),file=TRIM(outfile_ldr),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_ldr) )
                IF (lreadmeta_from_netcdf) THEN
                  CALL get_filename_diagfiles ('YULINDEPOLRATOBS', rs_meta(ista), outfile_ldr_obs)
                  CALL get_free_funit(ldrobsoutputunit(ista))
                  OPEN(ldrobsoutputunit(ista),file=TRIM(outfile_ldr_obs),iostat=ios, &
                       status=TRIM(openstatus), position=TRIM(openpos), action='write')
                  CALL openerrmsg ( ios, TRIM(outfile_ldr_obs) )
                END IF
              END IF
            END IF

            IF (lextdbz .AND. &
                (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
              CALL get_filename_diagfiles ('YUEXTREFL', rs_meta(ista), outfile_ze)
              CALL get_free_funit(ahoutputunit(ista))
              OPEN(ahoutputunit(ista),file=TRIM(outfile_ze),iostat=ios, &
                   status=TRIM(openstatus), position=TRIM(openpos), action='write')
              CALL openerrmsg ( ios, TRIM(outfile_ze) )

              IF (loutpolstd .OR. loutpolall) THEN
                CALL get_filename_diagfiles ('YUDIFFEXTREFL', rs_meta(ista), outfile_zed)
                CALL get_free_funit(adpoutputunit(ista))
                OPEN(adpoutputunit(ista),file=TRIM(outfile_zed),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_zed) )
              END IF
            END IF
          END IF

          IF (lreadmeta_from_netcdf .AND. (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0))) THEN
            CALL get_filename_diagfiles ('YURADREFLOBS', rs_meta(ista), outfile_zr_obs)
            CALL get_free_funit(radreflobsoutputunit(ista))
            OPEN(radreflobsoutputunit(ista),file=TRIM(outfile_zr_obs),iostat=ios, &
                 status=TRIM(openstatus), position=TRIM(openpos), action='write')
            CALL openerrmsg ( ios, TRIM(outfile_zr_obs) )
          END IF
        END DO
      END IF
    ELSE
      IF (num_radario > 0) THEN
        ! Output on asynchoneous IO-PEs:
        num_io_loc = num_radario_dom(idom)
      ELSE
        ! Output on the compute-PEs:
        num_io_loc = num_radar_dom(idom)
      END IF
      ! .. Here, parallel output is implemented in a way that the radar
      !    station output is distributed among the compute PEs in a
      !    round-robin fashion:
      DO ista_par = 1, nradsta, num_io_loc
        DO ista = ista_par, MIN(ista_par+num_io_loc-1, nradsta)
          ipe_rad = MOD(ista-1, num_io_loc) + radario_master_dom(idom)
          IF (ipe_rad == my_radar_id_dom(idom)) THEN

            IF (lout_geom) THEN
              CALL get_filename_diagfiles ('YURADGEOM', rs_meta(ista), outfile_hr)
              CALL get_free_funit(radgeomoutputunit(ista))
              OPEN(radgeomoutputunit(ista),file=TRIM(outfile_hr),iostat=ios, &
                   status=TRIM(openstatus), position=TRIM(openpos), action='write')
              CALL openerrmsg ( ios, TRIM(outfile_hr) )
            END IF

            IF (loutradwind) THEN
              CALL get_filename_diagfiles ('YURADWIND', rs_meta(ista), outfile_vr)
              CALL get_free_funit(radwindoutputunit(ista))
              OPEN(radwindoutputunit(ista),file=TRIM(outfile_vr),iostat=ios, &
                   status=TRIM(openstatus), position=TRIM(openpos), action='write')
              CALL openerrmsg ( ios, TRIM(outfile_vr) )
              IF (lreadmeta_from_netcdf) THEN
                CALL get_filename_diagfiles ('YURADWINDOBS', rs_meta(ista), outfile_vr_obs)
                CALL get_free_funit(radwindobsoutputunit(ista))
                OPEN(radwindobsoutputunit(ista),file=TRIM(outfile_vr_obs),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_vr_obs) )
              END IF
            END IF

            IF (loutdbz .OR. (loutradwind .AND. (lmds_vr .OR. (lsmooth .AND. lweightdbz)))) THEN
              CALL get_filename_diagfiles ('YURADREFL', rs_meta(ista), outfile_zr)
              CALL get_free_funit(radrefloutputunit(ista))
              OPEN(radrefloutputunit(ista),file=TRIM(outfile_zr),iostat=ios, &
                   status=TRIM(openstatus), position=TRIM(openpos), action='write')
              CALL openerrmsg ( ios, TRIM(outfile_zr) )
            END IF

            IF (loutdbz) THEN
              IF ( (loutpolstd .OR. loutpolall) .AND. &
                   (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4) ) THEN
                CALL get_filename_diagfiles ('YUDIFFREFL', rs_meta(ista), outfile_zdr)
                CALL get_free_funit(zdroutputunit(ista))
                OPEN(zdroutputunit(ista),file=TRIM(outfile_zdr),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_zdr) )
                IF (lreadmeta_from_netcdf) THEN
                  CALL get_filename_diagfiles ('YUDIFFREFOBS', rs_meta(ista), outfile_zdr_obs)
                  CALL get_free_funit(zdrobsoutputunit(ista))
                  OPEN(zdrobsoutputunit(ista),file=TRIM(outfile_zdr_obs),iostat=ios, &
                       status=TRIM(openstatus), position=TRIM(openpos), action='write')
                  CALL openerrmsg ( ios, TRIM(outfile_zdr_obs) )
                END IF

                CALL get_filename_diagfiles ('YUCORRCOEFF', rs_meta(ista), outfile_rhv)
                CALL get_free_funit(rhvoutputunit(ista))
                OPEN(rhvoutputunit(ista),file=TRIM(outfile_rhv),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_rhv) )
                IF (lreadmeta_from_netcdf) THEN
                  CALL get_filename_diagfiles ('YUCORRCOEFFOBS', rs_meta(ista), outfile_rhv_obs)
                  CALL get_free_funit(rhvobsoutputunit(ista))
                  OPEN(rhvobsoutputunit(ista),file=TRIM(outfile_rhv_obs),iostat=ios, &
                       status=TRIM(openstatus), position=TRIM(openpos), action='write')
                  CALL openerrmsg ( ios, TRIM(outfile_rhv_obs) )
                END IF

                CALL get_filename_diagfiles ('YUSPECDIFFPHASE', rs_meta(ista), outfile_kdp)
                CALL get_free_funit(kdpoutputunit(ista))
                OPEN(kdpoutputunit(ista),file=TRIM(outfile_kdp),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_kdp) )
                IF (lreadmeta_from_netcdf) THEN
                  CALL get_filename_diagfiles ('YUSPECDIFFPHASEOBS', rs_meta(ista), outfile_kdp_obs)
                  CALL get_free_funit(kdpobsoutputunit(ista))
                  OPEN(kdpobsoutputunit(ista),file=TRIM(outfile_kdp_obs),iostat=ios, &
                       status=TRIM(openstatus), position=TRIM(openpos), action='write')
                  CALL openerrmsg ( ios, TRIM(outfile_kdp_obs) )
                END IF

                IF (loutpolall) THEN
                  CALL get_filename_diagfiles ('YULINDEPOLRAT', rs_meta(ista), outfile_ldr)
                  CALL get_free_funit(ldroutputunit(ista))
                  OPEN(ldroutputunit(ista),file=TRIM(outfile_ldr),iostat=ios, &
                       status=TRIM(openstatus), position=TRIM(openpos), action='write')
                  CALL openerrmsg ( ios, TRIM(outfile_ldr) )
                  IF (lreadmeta_from_netcdf) THEN
                    CALL get_filename_diagfiles ('YULINDEPOLRATOBS', rs_meta(ista), outfile_ldr_obs)
                    CALL get_free_funit(ldrobsoutputunit(ista))
                    OPEN(ldrobsoutputunit(ista),file=TRIM(outfile_ldr_obs),iostat=ios, &
                         status=TRIM(openstatus), position=TRIM(openpos), action='write')
                    CALL openerrmsg ( ios, TRIM(outfile_ldr_obs) )
                  END IF
                END IF
              END IF

              IF (lextdbz .AND. &
                  (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
                CALL get_filename_diagfiles ('YUEXTREFL', rs_meta(ista), outfile_ze)
                CALL get_free_funit(ahoutputunit(ista))
                OPEN(ahoutputunit(ista),file=TRIM(outfile_ze),iostat=ios, &
                     status=TRIM(openstatus), position=TRIM(openpos), action='write')
                CALL openerrmsg ( ios, TRIM(outfile_ze) )

                IF (loutpolstd .OR. loutpolall) THEN
                  CALL get_filename_diagfiles ('YUDIFFEXTREFL', rs_meta(ista), outfile_zed)
                  CALL get_free_funit(adpoutputunit(ista))
                  OPEN(adpoutputunit(ista),file=TRIM(outfile_zed),iostat=ios, &
                       status=TRIM(openstatus), position=TRIM(openpos), action='write')
                  CALL openerrmsg ( ios, TRIM(outfile_zed) )
                END IF
              END IF
            END IF

            IF (lreadmeta_from_netcdf .AND. (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0))) THEN
              CALL get_filename_diagfiles ('YURADREFLOBS', rs_meta(ista), outfile_zr_obs)
              CALL get_free_funit(radreflobsoutputunit(ista))
              OPEN(radreflobsoutputunit(ista),file=TRIM(outfile_zr_obs),iostat=ios, &
                   status=TRIM(openstatus), position=TRIM(openpos), action='write')
              CALL openerrmsg ( ios, TRIM(outfile_zr_obs) )
            END IF
          END IF
        END DO
      END DO
    END IF

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  CONTAINS

    SUBROUTINE openerrmsg ( ios, outfile )
      INTEGER, INTENT(in) :: ios
      CHARACTER (LEN=*)   :: outfile
      IF ( ios /= 0 ) THEN
        CALL abort_run (my_radar_id, 20365, &
             'ERROR: problem in ' //TRIM(yzroutine)// &
             ': error opening output diagnosis file '//TRIM(outfile) , &
             'radar_process_output.f90, '//TRIM(yzroutine))
      END IF
    END SUBROUTINE openerrmsg

  END SUBROUTINE opendiagfiles

  SUBROUTINE closediagfiles ()

    !------------------------------------------------------------------------------
    !
    ! Description: Closes file units for diagnostic output files
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32)  :: yzroutine
    INTEGER             :: ista, ista_par, ipe_rad, num_io_loc
    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE closediagfiles
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'closediagfiles'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! .. Open the file units for the control output for each data set of
    !    each radar station:
    IF (.FALSE.) THEN
!!$ Left in the code, so that in future, one might go
!!$ back to output on node 0 only for whatever reason:
      ! .. Here, a parallel output is not implemented yet, so
      !    all the files for the control output for each radar
      !    have to be opened on PE 0
      IF (my_radario_id_dom(idom) == 0) THEN
        DO ista = 1, nradsta
          IF (lout_geom) THEN
            CLOSE(radgeomoutputunit(ista))
          END IF
          IF (loutradwind) THEN
            CLOSE(radwindoutputunit(ista))
            IF (lreadmeta_from_netcdf) THEN
              CLOSE (radwindobsoutputunit(ista))
            END IF
          END IF
          IF (loutdbz .OR. (loutradwind .AND. (lmds_vr .OR. (lsmooth .AND. lweightdbz)))) THEN
            CLOSE(radrefloutputunit(ista))
            IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
              CLOSE(zdroutputunit(ista))
              CLOSE(rhvoutputunit(ista))
              CLOSE(kdpoutputunit(ista))
              IF (lreadmeta_from_netcdf) THEN
                CLOSE(zdrobsoutputunit(ista))
                CLOSE(rhvobsoutputunit(ista))
                CLOSE(kdpobsoutputunit(ista))
              END IF
              IF (loutpolall) THEN
                CLOSE(ldroutputunit(ista))
                IF (lreadmeta_from_netcdf) CLOSE(ldrobsoutputunit(ista))
              END IF
            END IF

            IF (lextdbz .AND. &
                (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
              CLOSE(ahoutputunit(ista))
              IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
                CLOSE(adpoutputunit(ista))
              END IF
            END IF
          END IF
          IF (lreadmeta_from_netcdf .AND. (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0))) THEN
            CLOSE(radreflobsoutputunit(ista))
          END IF
        END DO
      END IF
    ELSE
      IF (num_radario_dom(idom) > 0) THEN
        ! Output on asynchoneous IO-PEs:
        num_io_loc = num_radario_dom(idom)
      ELSE
        ! Output on the compute-PEs:
        num_io_loc = num_radar_dom(idom)
      END IF
      ! .. Here, parallel output is implemented in a way that the radar
      !    station output is distributed among the compute PEs in a
      !    round-robin fashion:
      DO ista_par = 1, nradsta, num_io_loc
        DO ista = ista_par, MIN(ista_par+num_io_loc-1, nradsta)
          ipe_rad = MOD(ista-1, num_io_loc) + radario_master_dom(idom)
          IF (ipe_rad == my_radar_id_dom(idom)) THEN
            IF (lout_geom) THEN
              CLOSE(radgeomoutputunit(ista))
            END IF
            IF (loutradwind) THEN
              CLOSE(radwindoutputunit(ista))
              IF (lreadmeta_from_netcdf) THEN
                CLOSE(radwindobsoutputunit(ista))
              END IF
            END IF
            IF (loutdbz .OR. (loutradwind .AND. (lmds_vr .OR. (lsmooth .AND. lweightdbz)))) THEN
              CLOSE(radrefloutputunit(ista))
              IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
                CLOSE(zdroutputunit(ista))
                CLOSE(rhvoutputunit(ista))
                CLOSE(kdpoutputunit(ista))
                IF (lreadmeta_from_netcdf) THEN
                  CLOSE(zdrobsoutputunit(ista))
                  CLOSE(rhvobsoutputunit(ista))
                  CLOSE(kdpobsoutputunit(ista))
                END IF
                IF (loutpolall) THEN
                  CLOSE(ldroutputunit(ista))
                  IF (lreadmeta_from_netcdf) CLOSE(ldrobsoutputunit(ista))
                END IF
              END IF

              IF (lextdbz .AND. &
                  (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
                CLOSE(ahoutputunit(ista))
                IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
                  CLOSE(adpoutputunit(ista))
                END IF
              END IF
            END IF
            IF (lreadmeta_from_netcdf .AND. (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0))) THEN
              CLOSE(radreflobsoutputunit(ista))
            END IF
          END IF
        END DO
      END DO
    END IF

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE closediagfiles

  !==============================================================================
  !==============================================================================

  SUBROUTINE get_fileprefix_ascii_output (varname, rsm, fileprefix)
    CHARACTER(len=*),      INTENT(in)  :: varname
    type(radar_meta_type), INTENT(in)  :: rsm
    CHARACTER(len=*),      INTENT(inout) :: fileprefix

    fileprefix(:) = ' '
    WRITE (fileprefix, '(a,"_id-",i6.6,"_",a)') TRIM(varname), rsm%station_id, TRIM(rsm%scanname)

  END SUBROUTINE get_fileprefix_ascii_output

  !==============================================================================
  !==============================================================================
  
  !------------------------------------------------------------------------------
  !
  ! Helber procedures for ostream_create_filename() and composite_create_filename()
  ! for the search and replace of <patterns> in strings with actual values
  !
  !------------------------------------------------------------------------------

  SUBROUTINE get_next_key_from_pos (fpat, pos, key, ku, ko)

    CHARACTER(len=*), INTENT(in)    :: fpat    ! string in which to search for the next key, starting at the i'th position
    INTEGER         , INTENT(inout) :: pos     ! starting address in fpat for the search for key, will be incremented in preparation for the next search
    CHARACTER(len=*), INTENT(out)   :: key     ! name of the next key from fpat
    INTEGER         , INTENT(out)   :: ku, ko  ! start and end index of the key in fpat

    INTEGER :: len_fpat

    key(:) = ' '
    len_fpat = LEN_TRIM(fpat)
    ku = INDEX(fpat(pos:len_fpat), '<')
    IF (ku > 0) THEN
      ku = ku + pos - 1
      ko = INDEX(fpat(ku:len_fpat), '>')
      IF (ko > 0) THEN
        ko = ko + ku - 1
        key(:) = ' '
        key = fpat(ku:ko)
        ! Advance new pos to slightly after the beginning of the key, so that the next search will find the next key.
        !  Note that ko would be too much, because before the next search <key> will be replaced by its actual meaning,
        !  and the string length of it could be shorter than the <key>.
        pos = ku + 1
      ELSE
        ko = -1
        pos = len_fpat
      END IF
    ELSE
      ku = -1
      ko = -1
      pos = len_fpat
    END IF

  END SUBROUTINE get_next_key_from_pos

  SUBROUTINE replace_substr_with_value (fpat, ku, ko, cval)

    CHARACTER(len=*), INTENT(inout) :: fpat    ! string in which to replace the position ku:ko with cval
    INTEGER         , INTENT(in)    :: ku, ko  ! start and end index of the characters in fpat to replace with cval
    CHARACTER(len=*), INTENT(in)    :: cval    ! replacement string

    INTEGER :: len_fpat

    len_fpat = LEN_TRIM(fpat)
    IF (ku > 1 .AND. ko < len_fpat) THEN
      fpat = fpat(1:ku-1) // TRIM(cval) // fpat(ko+1:) // REPEAT(' ',len_fpat)
    ELSE IF (ku > 1) THEN
      fpat = fpat(1:ku-1) // TRIM(cval) // REPEAT(' ',len_fpat)
    ELSE
      fpat = TRIM(cval) // fpat(ko+1:) // REPEAT(' ',len_fpat)
    END IF

  END SUBROUTINE replace_substr_with_value
     
  !=========================================================================
  !
  ! Subroutine for control output in files YURADWIND, YURADREFL, ...
  !
  ! Purpose: Output a statistics of 3D polar data field "feld_name" of the
  ! radar described by the meta data structure rad_meta.
  ! The field must be given a short name (<= 8 characters) "feldname".
  ! When searching for the max and min of the field, the search
  ! is limited to positions masked to .true. in the logical 3D field "mask_polar",
  ! which needs to have the same dimensions as "feld_polar". If
  ! all mask_polar is .false., no max / min can be found and the
  ! statistics line will contain the value in "default_missing" instead.
  !
  ! NOTE: ONLY CALL ON ONE PROCESSOR AT A TIME!
  !
  !=========================================================================

  SUBROUTINE control_output(outunit, time_mod, rad_meta, feldname, feld_polar, mask_polar, &
       default_missing, nobs)


    TYPE(radar_meta_type),   INTENT(in)  :: rad_meta
    INTEGER,                 INTENT(in)  :: outunit       , & ! unit number of the (opened) output text file
                                            nobs              ! number of valid observations
    REAL(KIND=dp),           INTENT(in)  :: time_mod          ! time in seconds since model start
    REAL(KIND=dp),           INTENT(in)  :: feld_polar(rad_meta%naz,rad_meta%nra,rad_meta%nel), default_missing
    LOGICAL,                 INTENT(in)  :: mask_polar(rad_meta%naz,rad_meta%nra,rad_meta%nel)
    CHARACTER(len=*),        INTENT(in)  :: feldname

    INTEGER                              :: i, imaxloc(3), iminloc(3)
    CHARACTER(len=30)                    :: format_out

    IF (ldebug_radsim) WRITE (*,*) 'control_output on unit ', outunit, ' on proc ', my_radar_id

    WRITE(outunit,*)
    WRITE(outunit,*)
    WRITE(outunit,'(a)') REPEAT('=', 80)
    WRITE(outunit,'("time = ",I12," seconds",10x,"field = ",a)') NINT(time_mod), TRIM(feldname)
    WRITE(outunit,'(a)') REPEAT('=', 80)
    WRITE(outunit,*)
    WRITE(outunit,'("station_id  = ",I6.6)')      rad_meta%station_id
    WRITE(outunit,'("station_name  = ",A)')      TRIM(ADJUSTL(rad_meta%station_name))
    WRITE(outunit,'("lon      = ",F12.6)')    rad_meta%lon
    WRITE(outunit,'("lat      = ",F12.6)')    rad_meta%lat
    WRITE(outunit,'("alt_agl_mod  = ",F10.2," (from namelist)")')  &
                                               rad_meta%alt_agl_mod
    WRITE(outunit,'("alt_msl  = ",F10.2," (actually used for computations)")') &
                                              rad_meta%alt_msl
    WRITE(outunit,'("mod_msl  = ",F10.2," (model orography height at radar station)")') &
                                              rad_meta%msl_mod
    WRITE(outunit,'("alt_msl_true  = ",F10.2," (from meta data of obs data file)")') &
                                              rad_meta%alt_msl_true
    IF (lsmooth) THEN
      WRITE(outunit,'("Phi3     = ",F10.2)')    rad_meta%Phi3
      WRITE(outunit,'("Theta3   = ",F10.2)')    rad_meta%Theta3
      WRITE(outunit,'("dalpha   = ",F10.2)')    rad_meta%dalpha
      WRITE(outunit,'("alpha3_eff_0 = ",F10.3)') rad_meta%alpha3_eff_0
      WRITE(outunit,'("ngpsm_h  = ",I10)')    rad_meta%ngpsm_h
      WRITE(outunit,'("ngpsm_v  = ",I10)')    rad_meta%ngpsm_v
      format_out(:) = ' '
      WRITE(format_out,'("(A,",I3.3,"(x,F8.4,"" |""))")') rad_meta%ngpsm_h
      WRITE(outunit,TRIM(format_out)) 'xabscsm_h = ', &
           (rad_meta%xabscsm_h(i),i=1,rad_meta%ngpsm_h)
      WRITE(outunit,TRIM(format_out)) 'weigsm_h = ', &
           (rad_meta%weigsm_h(i),i=1,rad_meta%ngpsm_h)
      format_out(:) = ' '
      WRITE(format_out,'("(A,",I3.3,"(x,F8.4,"" |""))")') rad_meta%ngpsm_v
      WRITE(outunit,TRIM(format_out)) 'xabscsm_v = ', &
           (rad_meta%xabscsm_v(i),i=1,rad_meta%ngpsm_v)
      WRITE(outunit,TRIM(format_out)) 'weigsm_v = ', &
           (rad_meta%weigsm_v(i),i=1,rad_meta%ngpsm_v)
    END IF
    WRITE(outunit,'("az_start = ",F10.2)')    rad_meta%az_start
    WRITE(outunit,'("daz      = ",F10.2)')    rad_meta%az_inc
    WRITE(outunit,'("naz      = ",I10)')      rad_meta%naz
    WRITE(outunit,'("dra      = ",F10.2)')    rad_meta%ra_inc
    WRITE(outunit,'("nra      = ",I10)')      rad_meta%nra
    WRITE(outunit,'("nel      = ",I10)')      rad_meta%nel
    WRITE(outunit,*)
    format_out(:) = ' '
    WRITE(format_out,'("(A,",I3.3,"(x,F6.2))")') rad_meta%nel
    WRITE(outunit,TRIM(format_out)) 'ele = ', (rad_meta%el_arr(i),i=1,rad_meta%nel)
    WRITE(outunit,'("nobs     = ",I10)') nobs

    ! A line in YURADWIND with max, min and mean field:

    imaxloc = MAXLOC(array=feld_polar, mask=mask_polar)
    iminloc = MINLOC(array=feld_polar, mask=mask_polar)
    WRITE(outunit,'(A20,3A12,5x,A20,3A12)') "MAX ","ra_max","az_max","el_max", &
         "MIN ","ra_min","az_min","el_min"
    IF (ALL(imaxloc > 0)) THEN
      WRITE(outunit,'(es20.5,3i12,5x,es20.5,3i12)')  &
           feld_polar(imaxloc(1),imaxloc(2),imaxloc(3)), &
!!$           imaxloc(1), &
!!$           imaxloc(2), &
           imaxloc(2), &
           imaxloc(1), &
           imaxloc(3), &
           feld_polar(iminloc(1),iminloc(2),iminloc(3)), &
!!$           iminloc(1), &
!!$           iminloc(2), &
           iminloc(2), &
           iminloc(1), &
           iminloc(3)
    ELSE
      WRITE(outunit,'(es20.5,3i12,5x,es20.5,3i12)')  &
           default_missing, &
           -1, &
           -1, &
           -1, &
           default_missing, &
           -1, &
           -1, &
           -1
    END IF

  END SUBROUTINE control_output


  !==============================================================================
  !==============================================================================
  !
  ! Conversion of a string to either lower or upper case characters

  FUNCTION to_lower_or_upper(string, target_case) RESULT(changed_string)

    IMPLICIT NONE

    ! in/out
    CHARACTER(LEN=*),    INTENT(IN) :: string
    CHARACTER(LEN=*),    INTENT(IN) :: target_case
    CHARACTER(LEN=LEN_TRIM(string)) :: changed_string

    ! local
    INTEGER :: jletter, lacode, uacode, lzcode, uzcode, icode

    CHARACTER(LEN=*), PARAMETER :: case_upper    = "upper"
    CHARACTER(LEN=*), PARAMETER :: case_lower    = "lower"

    !------------------------------------------------------

    lacode=IACHAR('a')
    lzcode=IACHAR('z')
    uacode=IACHAR('A')
    uzcode=IACHAR('Z')

    changed_string = string

    IF(TRIM(target_case) == case_upper) THEN
      DO jletter=1, LEN_TRIM(string)
        icode = IACHAR(string(jletter:jletter))
        IF (icode >= lacode .AND. icode <= lzcode) THEN
          changed_string(jletter:jletter) = ACHAR(icode - lacode + uacode)
        END IF
      ENDDO
    ELSEIF(TRIM(target_case) == case_lower) THEN
      DO jletter=1, LEN_TRIM(string)
        icode = IACHAR(string(jletter:jletter))
        IF (icode >= uacode .AND. icode <= uzcode) THEN
          changed_string(jletter:jletter) = ACHAR(icode - uacode + lacode)
        END IF
      ENDDO
    ELSE
      RETURN
    ENDIF

  END FUNCTION to_lower_or_upper

  !===============================================================================
  !===============================================================================

  !============================================================================
  !
  ! Routine for writing a READY file. It has to be called on one PE only!
  !
  ! The rationale of a ready file is that it does not yet exist, but
  ! for the case of "leftovers" from a crashed run or restart run,
  ! we open it with status "replace" rather than "new".
  !
  !============================================================================

  SUBROUTINE write_ready_radar ( path, file_pattern, content, ierror )

    CHARACTER(len=*), INTENT(in)  :: path, file_pattern, content
    INTEGER        ,  INTENT(out) :: ierror

    INTEGER :: funit, iret
    CHARACTER(len=cmaxlen) :: ofilename, rdy_filename, tmp_filename
    CHARACTER(len=*), PARAMETER :: yzroutine = 'emvorado::write_ready_radar'
    CHARACTER(len=*), PARAMETER :: tmp_prefix = ".."
    CHARACTER(len=14) :: model_starttime, model_validtime
    CHARACTER(len=8)  :: ddhhmmss_validtime
    CHARACTER(len=LEN(file_pattern)) :: fpat
    CHARACTER(len=20) :: key
    INTEGER           :: ku, ko, pos
    LOGICAL           :: have_time

    INTERFACE
      FUNCTION private_rename(old_filename, new_filename) RESULT(iret) BIND(C, NAME='rename')
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_char
        INTEGER(c_int) :: iret
        CHARACTER(c_char), DIMENSION(*), INTENT(in) :: old_filename
        CHARACTER(c_char), DIMENSION(*), INTENT(in) :: new_filename
      END FUNCTION private_rename
    END INTERFACE

    ofilename(:) = ' '
    ierror = 0
    
    !--------------------------------------------------------------------------
    ! .. Create filename in a way that times are organized in batches of length
    !     "content_dt" (sedonds) relative to a reference time "content_tref"
    !     (seconds since model start):

    ddhhmmss_validtime = get_model_time_ddhhmmss()
    model_starttime    = get_datetime_ini()
    model_validtime    = get_datetime_act()
    
    ofilename(:) = ' '

    IF (LEN_TRIM(file_pattern) <= 0) THEN

      ! Default name:
      !ofilename = TRIM(path)//'READY_EMVORADO_'//TRIM(model_validtime)
      ofilename = 'READY_EMVORADO_'//TRIM(model_validtime)

    ELSE
      
      !--------------------------------------------------------------------------
      ! Set the user defined file name pattern

      ! .. List of valid <keys> for filename patterns:
      !
      !       - <tmodelini>    : will result in YYYYMMDDhhmmss  (absolute datetime of model start) (optional)
      !       - <tact>         : will result in YYYYMMDDhhmmss  (absolute datetime of actual model time) \
      !       - <tvvzact>      : will result in DDhhmmss  (forecast lead time of actual model time)      / (one of both is mandatory)
      
      fpat = file_pattern ! copy original pattern

      have_time     = .FALSE.
      
      pos = 1
      DO
        ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
        CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
        IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
        ! .. and, if it is a time key, replace by its actual value:
        SELECT CASE (TRIM(key))
        CASE ('<tmodelini>')
          ! optional key, therefore no have_time = .TRUE.:
          CALL replace_substr_with_value (fpat, ku, ko, TRIM(model_starttime))
        CASE ('<tact>')
          have_time = .TRUE.
          CALL replace_substr_with_value (fpat, ku, ko, TRIM(model_validtime))
        CASE ('<tvvzact>')
          have_time = .TRUE.
          CALL replace_substr_with_value (fpat, ku, ko, TRIM(ddhhmmss_validtime))
        END SELECT
      END DO

      ! At last, check if any other keys or "<" or ">" are still remaining in the fpat. If yes, it is a wrong key:
      pos = 1
      DO
        ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
        CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
        IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
        ! .. if there is still a key left in fpat, it is a wrong key:
        WRITE (*, '(a)') &
             'ERROR '//TRIM(yzroutine)//': file_pattern "'//TRIM(file_pattern)// &
             '" contains unknown key '//TRIM(key)
        ierror = 1
        RETURN
      END DO

      ! .. if still some remaining single '<' or '>' characters are in fpat, there is an error:
      IF (INDEX(fpat,'<') > 0 .OR. INDEX(fpat,'>') > 0) THEN
        WRITE (*, '(a)') &
             'ERROR '//TRIM(yzroutine)//': file_pattern "'//TRIM(file_pattern)// &
             '" contains garbled "<" or ">" characters'
        ierror = 2
        RETURN
      END IF

      ! .. check if all mandatory keys have been specified in file_pattern:
      IF (.NOT.have_time) THEN
        WRITE (*, '(a)') &
             'ERROR '//TRIM(yzroutine)//': file_pattern "'//TRIM(file_pattern)// &
             '" does not contain a <key> for actual time (<tact>, <tvvzact>)'
        ierror = 3
        RETURN
      END IF

      ! .. at this stage, fpat contains a valid filename:
      !ofilename = TRIM(path)//TRIM(fpat)
      ofilename = TRIM(fpat)

    END IF
    
    ! Actually create ready file.
    !
    ! This procedure is carried out in two steps: First, a file with
    ! the name "tmp_prefix+rdy_filename" is created. After the file
    ! has been closed, it is then renamed into "rdy_filename" in a
    ! second step.
    ! This detour is necessary when another process polls the output
    ! directory and relies on a "complete" ready file.
    tmp_filename = TRIM(path)//tmp_prefix//TRIM(ofilename)
    rdy_filename = TRIM(path)//TRIM(ofilename)
    CALL get_free_funit( funit )

    OPEN(unit=funit, file=TRIM(tmp_filename), status='replace', form='formatted', iostat=iret)
    IF (iret == 0) THEN
      WRITE (funit,*) content
      CLOSE (funit)
    ELSE
      WRITE (*, '(a)') 'ERROR '//TRIM(yzroutine)//': could not open '//TRIM(tmp_filename)
      ierror = iret
      RETURN
    END IF

    iret = util_rename(TRIM(tmp_filename), TRIM(rdy_filename))
    IF (iret /= 0) THEN
      WRITE (*, '(a)') 'ERROR '//TRIM(yzroutine)//': could not rename '//TRIM(tmp_filename)// &
           ' to '//TRIM(rdy_filename)
      ierror = iret
      RETURN
    END IF

  CONTAINS

    FUNCTION util_rename(old_filename, new_filename) RESULT(iret)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_null_char
      INTEGER :: iret
      CHARACTER(len=*), INTENT(in) :: old_filename
      CHARACTER(len=*), INTENT(in) :: new_filename
      iret = private_rename(TRIM(old_filename)//c_null_char, TRIM(new_filename)//c_null_char)
    END FUNCTION util_rename

  END SUBROUTINE write_ready_radar
  
  !===============================================================================

#ifdef GRIBAPI
  SUBROUTINE check_codes_err (griberr, cident, caller, increrr)

    INTEGER          , INTENT(in)    :: griberr
    CHARACTER(len=*) , INTENT(in)    :: cident, caller
    INTEGER, OPTIONAL, INTENT(inout) :: increrr

    CHARACTER(len=cmaxlen)           :: griberrmsg

    IF (griberr /= GRIB_SUCCESS) THEN
      griberrmsg(:) = ' '
      CALL codes_get_error_string(griberr, griberrmsg)
      WRITE (*,'(a)') 'ERROR codes_set('//TRIM(cident)//') in '//TRIM(caller)//': '//TRIM(griberrmsg)
    END IF
    
    ! update incremented error if present
    IF (PRESENT(increrr)) increrr = increrr + griberr
    
  END SUBROUTINE check_codes_err
#endif

  !===============================================================================

#ifdef NETCDF
  FUNCTION check_nc ( istat, rinfo ) RESULT (ostat)
    INTEGER, INTENT(in)          :: istat
    CHARACTER(len=*), INTENT(in) :: rinfo
    INTEGER                      :: ostat

    ! .. return the error status, so that the below error handling might be
    !     relaxed in the future (e.g., WARNING instead of abort_run())
    !     and the calling routine can decide what to do:
    ostat = istat
    IF (istat /= NF90_NOERR) THEN
      CALL abort_run (my_radar_id, 17077, &
           'ERROR in writing NetCDF: '//TRIM(rinfo)//': '//TRIM(NF90_strerror(istat)), &
           'radar_output_methods.f90, check_nc()')
    END IF

  END FUNCTION check_nc
#endif


END MODULE radar_output_utils
