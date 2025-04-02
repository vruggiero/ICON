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

MODULE radar_output_driver

!------------------------------------------------------------------------------
!
! Description: Modules of the radar forward operator EMVORADO for processing
!              of the various output methods, data and formats:
!              volume data output, feedback file output and reflectivity composite
!              generation and output.
!              This module contains methods to collect/MPI-copy simulated radar data
!              from the compute PEs to the output PEs, for reading obs data files,
!              for producing superobservations. It uses methods from another module radar_output_methods.f90
!              (soon to be splitted from the present module!!!!)
!              for writing feedback files, for writing volume data files (ASCII, NETCDF, or BIN-format)
!              and for writing grib2 reflectivity composites on rotated lat/lon grids.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!

!!$ TODO:
!!$ - height above radar station in feedback instead of height MSL? However, station height is also given in the feedback file header ...

  USE radar_kind, ONLY : dp
  
  USE radar_interface, ONLY : &
       get_runtime_timings, &
       get_rotlatlon_domain_for_superobing

  !------------------------------------------------------------------------------


  USE radar_parallel_utilities, ONLY :  &
       global_values_radar,             &
       get_idims_all_onestation,        &
       gatherv_values_radar

  !------------------------------------------------------------------------------

  USE radar_data,               ONLY :          &
       nradsta_max, cmaxlen,                    &
       radar_meta_type, radar_data_type,        &
       ydate_ini_mod,                           &
       idom, ndoms_max,                         &
       num_compute_fwo,   & ! number of compute PEs
       num_radar,         & ! number of radar PEs (num_compute + num_radario)
       num_radario,       & ! number of radar-IO PEs
       my_radar_id,       & ! rank of this PE in the radar communicator (cart+radario)
       my_radario_id,     & ! rank of this PE in the (asynchroneous) radario communicator
       icomm_radar,       & ! communicator for the group of radar-IO PEs + compute PEs
       icomm_radario,     & ! communicator for the group of radar-IO PEs
       num_radar_dom,     & ! number of radar PEs (num_compute + num_radario_dom) per radar-active model domain
       num_radario_dom,   & ! number of radar-IO PEs per radar-active model domain
       radario_master_dom,& ! root-PEs of the radario group for each active radar domain (in the radar_dom comm., not radario_dom-comm.!!!)
       icomm_radar_dom,   & ! communicator for the group of radar-IO PEs of each domain + compute PEs
       icomm_radario_dom, & ! communicator for the group of radar-IO PEs for each domain
       lradar_pe_dom,     & ! indicates whether this is a radar PE for a certain domain or not (compute or radar-IO)
       lradario_pe_dom,   & ! indicates whether this is a radar-IO PE for a certain domain or not
       my_radar_id_dom,   & ! rank of this PE in the radar communicator (cart+radario_dom)
       my_radario_id_dom, & ! rank of this PE in the (asynchroneous) radario communicator icomm_radario_dom
       i_fwo_bubbles,     & ! Timing flag
       i_fwo_composites,  & ! Timing flag
       i_fwo_comm,        & ! Timing flag for MPI-communications
       i_fwo_out,         & ! Timing flag for output (collecting simulated data on one PE per station, sorting, ASCII-output,
                            !  reading obs data, producing feedback files)
       i_fwo_barrier,     & ! Timing flag for barrier waiting in MPI-communications (measure for load imbalance)
       raddeg, degrad,             &
       rs_meta, rs_data, cart_data, dbz_meta, &
       nradsta
  
  USE radar_data_namelist, ONLY :  &
       ldebug_radsim, loutradwind, &
       loutdbz, loutpolstd, loutpolall, lextdbz, lout_geom, &
       lweightdbz, lfall, lonline, lsode, lsmooth, lreadmeta_from_netcdf, &
       lmds_z, lmds_vr, &
       lfdbk_output, &
       supob_azi_maxsector_vr, supob_cart_resolution, &
       supob_ave_width, supob_minrange_vr, supob_minrange_z, supob_vrw, supob_rfl, &
       supob_lowthresh_z_obs, supob_lowthresh_z_sim, &
       lcomm_nonblocking_output

  USE radar_output_station,      ONLY : output_my_ista
  USE radar_output_station_smth, ONLY : output_my_ista_smth
  
  !------------------------------------------------------------------------------

#ifndef NOMPI
  USE mpi
#endif

!================================================================================
!================================================================================

  IMPLICIT NONE

#ifdef NOMPI
  include "nompi_mpif.h"
#endif

!================================================================================
!================================================================================

  PRIVATE

  PUBLIC ::  init_cart_info,                  &
             output_radar, output_radar_smth

  !==============================================================================
  ! Interface blocks for overloaded procedures:

  !==============================================================================
  ! Public and Private Subroutines

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  !------------------------------------------------------------------------------
  !+ Module procedure in radar_src for the output in case of lsmooth = .FALSE.
  !------------------------------------------------------------------------------

  SUBROUTINE output_radar(lcalc,time_mod)

    !------------------------------------------------------------------------------
    !
    ! Description: Gathers the radar data on the radar grid from each radar
    !              from all PEs to a single output PE per radar and writes
    !              the results to disc. Currently this is done as
    !              ASCII, ASCII-gzip, Fortran binary or NetCDF (internally compressed),
    !              as well as netcdf feedback files for data assimilation.
    !
    !              Has to be called on icomm_radar_dom(idom), i.e., all compute-PEs
    !              and the radar-IO-PEs allocated for the actual model domain idom.
    !              if called also on other PEs, it will hang up!
    !
    !              THIS PROCEDURE IS FOR SIMULATED RADAR DATA WITHOUT SPATIAL SMOOTHING!
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                        :: lcalc(nradsta_max)
    REAL (KIND=dp), INTENT(IN)     :: time_mod       ! seconds since model start

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) :: yzroutine
    CHARACTER (LEN=80) :: yzerrmsg

    INTEGER  :: i, ista, mpierror, ista_par, my_ista, ista_this_pe_counter, &
                ipe_rad, nobs_ista, izerror, nradsta_per_pe, num_io_loc

    INTEGER, ALLOCATABLE :: my_ista_vec(:)

    TYPE(radar_data_type)     :: rs_data_all(nradsta_max)

    INTEGER  :: ipos_all(num_radar_dom(idom)), idims_all(num_radar_dom(idom),nradsta_max), &
                message_tag, mpi_sendrequests(nradsta_max*20), &
                istatus(MPI_STATUS_SIZE)

    ! Flag for parallel output of each radar on different PEs. This
    ! flag indicates whether this PE has some radar data to output.
    LOGICAL              :: output_my_radar
    LOGICAL, ALLOCATABLE :: output_my_radar_vec(:)

    ! Flag if geometry files (hrsim, ersim) have been written to a file:

    LOGICAL, SAVE :: geom_written(nradsta_max,ndoms_max) = .FALSE., &
                     geom_all_written(ndoms_max) = .FALSE.

    ! Buffer for MPI_ALLREDUCE:
    LOGICAL :: lsendbuf(nradsta_max)


    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE output_radar
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'output_radar'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    IF (.NOT. (loutradwind .OR. loutdbz) .AND. ldebug_radsim) THEN
      WRITE (*,*) 'MESSAGE: '//TRIM(yzroutine)// &
           'No output of simulated radar data since loutradwind = loutdbz = .FALSE.'
      RETURN
    END IF

    ipos_all(:) = 0

    ! Distribute the radar data of different radars to all output PEs:
    !=================================================================

    IF (num_radario > 0) THEN
      ! Output on asynchoneous IO-PEs:
      num_io_loc = num_radario_dom(idom)
      IF (num_io_loc > 0 .AND. num_io_loc <= num_radario) THEN
        nradsta_per_pe = MAX( (nradsta-1)/num_io_loc+1 , 1)
        ALLOCATE (my_ista_vec(nradsta_per_pe))
        ALLOCATE (output_my_radar_vec(nradsta_per_pe))
        my_ista_vec = -1
        output_my_radar_vec = .FALSE.
      END IF
      ! Init of counter for no. of stations per output PE (asynchr. IO):
      ista_this_pe_counter = 0
    ELSE
      ! Output on the compute-PEs:
      num_io_loc = num_radar_dom(idom)
    END IF

#ifdef __COSMO__
    CALL get_runtime_timings (i_fwo_out)
#endif
#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_barrier)
#endif

    ! This barrier synchronizes output- and compute-PEs in case of asynchroneous output,
    !  and the  global_values_radar (geom_written(1:nradsta,idom), ...) back-propagates
    !  the info on written geometry infos from the asyn. IO PEs to all radar PEs of this domain:
    IF (num_radar > 1 .AND. num_radario > 0) THEN
      IF (lradar_pe_dom(idom)) THEN
        CALL MPI_BARRIER(icomm_radar_dom(idom), izerror)
        CALL global_values_radar (geom_written(1:nradsta,idom), nradsta, 'OR', icomm_radar_dom(idom),    &
                                  -1, yzerrmsg, mpierror)
      END IF
    END IF

#ifdef __COSMO__
    CALL get_runtime_timings (i_fwo_barrier)
#endif

    IF (num_radar > 1) THEN
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_comm)
#endif
      DO ista = 1, nradsta
        ! Prepare below calls to gatherv_values_radar(). Contains a global communication,
        !   so it is a synchro point for icomm_radar_dom(idom):
        CALL get_idims_all_onestation (rs_data(ista)%ind_intp(:,2), num_radar_dom(idom), &
             icomm_radar_dom(idom), idims_all(:,ista), yzerrmsg, mpierror)
      END DO
      ! Prepare the non-blocking comm below:
      message_tag = 0
      mpi_sendrequests(:) = MPI_REQUEST_NULL
    END IF

    ! loop over num_radario-sized chunks of radar stations:
    chunkloop: DO ista_par = 1, nradsta, num_io_loc

#ifdef __ICON__
        CALL get_runtime_timings (i_fwo_comm)
#endif

      output_my_radar = .FALSE.

      ! .. Distinction between single proc and parallel runs:
      IF (num_radar > 1) THEN

        IF (num_radario == 0) THEN
          message_tag = 0
          mpi_sendrequests(:) = MPI_REQUEST_NULL
        END IF

        istaloop: DO ista = ista_par, MIN(ista_par+num_io_loc-1, nradsta)

          IF (.NOT.lcalc(ista)) CYCLE istaloop

          ! determine the number of the PE for output of the radar with no. ista
          !  in the combined communicator icomm_radar:
          ipe_rad = MOD(ista-1, num_io_loc) + radario_master_dom(idom)

          IF (my_radar_id_dom(idom) == ipe_rad .AND. num_radario_dom(idom) > 0 .AND. .NOT. lradario_pe_dom(idom)) THEN
            WRITE (*,'(a,i4,a)') 'ERROR '//TRIM(yzroutine)//' for asynchroneous radar output: '// &
                 'output assigned to this PE (my_radar_id=', my_radar_id, ') while it is not an output PE!'
          END IF

          ! .. Gather the values for this radar on the PE with no. ipe_rad:

          IF ( (lout_geom .AND. (lonline .OR. .NOT.geom_written(ista,idom))) .OR. lreadmeta_from_netcdf ) THEN

            ! Gather height of radar bin above MSL:
            message_tag = message_tag + 1
            CALL gatherv_values_radar(rs_data(ista)%hl_loc, rs_data_all(ista)%hl_loc, ipos_all, nobs_ista, &
                 num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                 lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))

          END IF

          IF ( (lout_geom .AND. (lonline .OR. .NOT.geom_written(ista,idom))) .OR. lreadmeta_from_netcdf &
                                                                             .OR. (loutradwind .AND. lfall) ) THEN

            ! Gather local elevation angle:
            message_tag = message_tag + 1
            CALL gatherv_values_radar(rs_data(ista)%el_loc, rs_data_all(ista)%el_loc, ipos_all, nobs_ista, &
                 num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                 lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))

          END IF

          IF ( (lout_geom .OR. lreadmeta_from_netcdf) .AND. lonline ) THEN

            ! Gather arc distance from radar bin to radar station at MSL height:
            message_tag = message_tag + 1
            CALL gatherv_values_radar(rs_data(ista)%s_loc, rs_data_all(ista)%s_loc, ipos_all, nobs_ista, &
                 num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                 lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))

          END IF

          IF(loutradwind) THEN

            ! Gather radial winds:
            message_tag = message_tag + 1
            CALL gatherv_values_radar(rs_data(ista)%radwind_mod, rs_data_all(ista)%radwind_mod, ipos_all, nobs_ista, &
                 num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                 lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))

            IF(lfall) THEN
              ! Gather hydrometeor fallspeed:
              message_tag = message_tag + 1
              CALL gatherv_values_radar(rs_data(ista)%vt_mod, rs_data_all(ista)%vt_mod, ipos_all, nobs_ista, &
                   num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                   lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))

            END IF

          END IF

          IF(loutdbz .OR. (loutradwind .AND. lmds_vr)) THEN

            ! Gather radar reflectivity:
            message_tag = message_tag + 1
            CALL gatherv_values_radar(&
                 rs_data(ista)%zh_radar_mod, rs_data_all(ista)%zh_radar_mod, &
                 ipos_all, nobs_ista, num_radar_dom(idom), ipe_rad, &
                 my_radar_id_dom(idom), icomm_radar_dom(idom), &
                 yzerrmsg, mpierror, &
                 lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, &
                 isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))

            IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN

              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%zv_radar_mod, rs_data_all(ista)%zv_radar_mod, &
                   ipos_all, nobs_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))

              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%rrhv_radar_mod, rs_data_all(ista)%rrhv_radar_mod, &
                   ipos_all, nobs_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))

              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%irhv_radar_mod, rs_data_all(ista)%irhv_radar_mod, &
                   ipos_all, nobs_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))

              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%kdp_radar_mod, rs_data_all(ista)%kdp_radar_mod, &
                   ipos_all, nobs_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))

              IF (loutpolall) THEN
                message_tag = message_tag + 1
                CALL gatherv_values_radar(&
                     rs_data(ista)%zvh_radar_mod, rs_data_all(ista)%zvh_radar_mod, &
                     ipos_all, nobs_ista, num_radar_dom(idom), ipe_rad, &
                     my_radar_id_dom(idom), icomm_radar_dom(idom), &
                     yzerrmsg, mpierror, &
                     lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, &
                     isrequest=mpi_sendrequests(message_tag), &
                     idims_all_in=idims_all(:,ista))
              END IF

            END IF

            IF (lextdbz .AND. &
                (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN

              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%ah_radar_mod, rs_data_all(ista)%ah_radar_mod, &
                   ipos_all, nobs_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))

              IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
                message_tag = message_tag + 1
                CALL gatherv_values_radar(&
                     rs_data(ista)%adp_radar_mod, rs_data_all(ista)%adp_radar_mod, &
                     ipos_all, nobs_ista, num_radar_dom(idom), ipe_rad, &
                     my_radar_id_dom(idom), icomm_radar_dom(idom), &
                     yzerrmsg, mpierror, &
                     lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, &
                     isrequest=mpi_sendrequests(message_tag), &
                     idims_all_in=idims_all(:,ista))
              END IF

            END IF

          END IF

          ! hash-index of radar coordinate position (azi, range, ele):
          message_tag = message_tag + 1
          CALL gatherv_values_radar(rs_data(ista)%ind_intp(:,2), rs_data_all(ista)%radpos_all, ipos_all, nobs_ista, &
               num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
               lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
               idims_all_in=idims_all(:,ista))

          ! keep track of which radar is output on which processor:
          IF (my_radar_id_dom(idom) == ipe_rad) THEN

            IF (ldebug_radsim) THEN
              WRITE (*,'(a,i0,a,i0,a,i0)') TRIM(yzroutine)//': distribution of stations over PEs: ista = ', &
                   ista, ' ipe_rad = ', ipe_rad, ' my_radar_id_dom = ', my_radar_id_dom(idom)
            END IF

            rs_data_all(ista)%nobs = nobs_ista

            ! Book-keeping which ista shall be output on this processor:
            IF (num_radario > 0) THEN
              ! For asynchroneous output, it is the round robin collection of one or several istas,
              !  to be output later in an extra chunk loop:
              ista_this_pe_counter = ista_this_pe_counter + 1
              my_ista_vec(ista_this_pe_counter) = ista
              output_my_radar_vec(ista_this_pe_counter) = lcalc(ista)
            ELSE
              ! For synchroneous output, it is just the ista from this chunk cycle,
              !  to be output directly within this chunk loop.
              my_ista = ista
              output_my_radar = lcalc(ista)
            END IF

          END IF

        END DO istaloop

      ELSE   ! num_radar == 1

        ! Local and global radar data structure are identical
        ! No addinitional memory required:

        IF (.NOT.lcalc(ista_par)) CYCLE chunkloop

        ! Book-keeping which ista shall be output on this processor in this chunk-cylce:
        my_ista = ista_par
        output_my_radar = .TRUE.

        rs_data_all(my_ista)%nobs = SIZE(rs_data(my_ista)%ind_intp,1)

        CALL nullify_pointers (my_ista)

        rs_data_all(my_ista)%radwind_mod    => rs_data(my_ista)%radwind_mod
        rs_data_all(my_ista)%vt_mod         => rs_data(my_ista)%vt_mod
        rs_data_all(my_ista)%zh_radar_mod   => rs_data(my_ista)%zh_radar_mod
        rs_data_all(my_ista)%ah_radar_mod   => rs_data(my_ista)%ah_radar_mod
        rs_data_all(my_ista)%zv_radar_mod   => rs_data(my_ista)%zv_radar_mod
        rs_data_all(my_ista)%rrhv_radar_mod => rs_data(my_ista)%rrhv_radar_mod
        rs_data_all(my_ista)%irhv_radar_mod => rs_data(my_ista)%irhv_radar_mod
        rs_data_all(my_ista)%kdp_radar_mod  => rs_data(my_ista)%kdp_radar_mod
        rs_data_all(my_ista)%adp_radar_mod  => rs_data(my_ista)%adp_radar_mod
        rs_data_all(my_ista)%zvh_radar_mod  => rs_data(my_ista)%zvh_radar_mod
        rs_data_all(my_ista)%radpos_all     => rs_data(my_ista)%ind_intp(:,2)
        rs_data_all(my_ista)%hl_loc         => rs_data(my_ista)%hl_loc
        rs_data_all(my_ista)%el_loc         => rs_data(my_ista)%el_loc
        rs_data_all(my_ista)%s_loc          => rs_data(my_ista)%s_loc

      END IF ! num_radar > 1


      IF (num_radar > 1 .AND. num_radario == 0 .AND. lcomm_nonblocking_output) THEN
        ! Wait for the non-blocking send requests to complete:
        DO i = 1, message_tag
          CALL MPI_WAIT(mpi_sendrequests(i), istatus, mpierror)
          mpi_sendrequests(i) = MPI_REQUEST_NULL
        END DO
      END IF


#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_comm)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_out)
#endif


      IF (output_my_radar .AND. (num_radario == 0 .OR. num_radar == 1)) THEN


        ! .. Synchroneous output on the compute-PEs (in this case equal to the radar-PEs)

        ista = my_ista

        ! .. I'm the output_PE and have some data: output the data of this station:
        !     ( output_my_ista() contains no further MPI communication )
        CALL output_my_ista (time_mod, ista, &
             rs_data_all(ista)%nobs, rs_data_all(ista)%radpos_all, &
             rs_data_all(ista)%vt_mod, rs_data_all(ista)%radwind_mod, &
             rs_data_all(ista)%zh_radar_mod, rs_data_all(ista)%ah_radar_mod, &
             rs_data_all(ista)%zv_radar_mod, &
             rs_data_all(ista)%rrhv_radar_mod, rs_data_all(ista)%irhv_radar_mod, &
             rs_data_all(ista)%kdp_radar_mod, rs_data_all(ista)%adp_radar_mod, &
             rs_data_all(ista)%zvh_radar_mod, &
             rs_data_all(ista)%hl_loc, rs_data_all(ista)%el_loc, rs_data_all(ista)%s_loc, &
             geom_written(ista,idom))

        IF (num_radar > 1) THEN

          ! .. In this case, the memory has been allocated to the pointers
          !    by gatherv_values_radar(). Therefore deallocate it.
          CALL dealloc_pointers (ista)

        ELSE

          ! .. In this case, the pointers simply point to another
          !    pointers, to whom the memory has been allocated. Therefore,
          !    simply break the connection by NULLIFYING. The memory
          !    of the original pointers will be deallocated at the end
          !    of organize_radar().
          CALL nullify_pointers (ista)

        END IF

      END IF ! output_my_radar = .true.

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_out)
#endif

    END DO chunkloop ! loop over (num_io_loc packages of) stations

    !----------------------------------------------------------------

#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_comm)
#endif

    IF (num_radar > 1 .AND. num_radario > 0 .AND. lcomm_nonblocking_output) THEN
      ! Wait for the non-blocking send requests to complete:
      DO i = 1, message_tag
        CALL MPI_WAIT(mpi_sendrequests(i), istatus, mpierror)
      END DO
    END IF

    !----------------------------------------------------------------

    IF (num_radario_dom(idom) > 0 .AND. lradario_pe_dom(idom) .AND. num_radar > 1) THEN

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_out)
#endif

      ! .. Asynchroneous output on the radario-PEs (I'm one of them):

      DO i=1, ista_this_pe_counter

        IF (output_my_radar_vec(i)) THEN

          ista = my_ista_vec(i)

          ! .. I'm the output_PE and have some data: output the data of this station:
          !     ( output_my_ista() contains no further MPI communication )
          CALL output_my_ista (time_mod, ista, &
               rs_data_all(ista)%nobs, rs_data_all(ista)%radpos_all, &
               rs_data_all(ista)%vt_mod, rs_data_all(ista)%radwind_mod, &
               rs_data_all(ista)%zh_radar_mod, rs_data_all(ista)%ah_radar_mod, &
               rs_data_all(ista)%zv_radar_mod, &
               rs_data_all(ista)%rrhv_radar_mod, rs_data_all(ista)%irhv_radar_mod, &
               rs_data_all(ista)%kdp_radar_mod, rs_data_all(ista)%adp_radar_mod, &
               rs_data_all(ista)%zvh_radar_mod, &
               rs_data_all(ista)%hl_loc, rs_data_all(ista)%el_loc, rs_data_all(ista)%s_loc, &
               geom_written(ista,idom))

          ! .. In this case, the memory has been allocated to the pointers
          !    by gatherv_values_radar(). Therefore deallocate it:
          CALL dealloc_pointers (ista)

        END IF

      END DO

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_out)
#endif

    END IF

    !----------------------------------------------------------------

    IF (ALLOCATED(output_my_radar_vec)) DEALLOCATE (output_my_radar_vec)
    IF (ALLOCATED(my_ista_vec))         DEALLOCATE (my_ista_vec)


    !----------------------------------------------------------------
    ! .. A barrier is set here, so that subsequent timing results
    !     are not biased, in case num_radar > nradsta:

#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_barrier)
#endif
    IF (num_radar > 1) THEN
      IF (num_radario == 0) THEN
        CALL MPI_BARRIER(icomm_radar, izerror)
      ELSE IF (lradario_pe_dom(idom)) THEN
        CALL MPI_BARRIER(icomm_radario_dom(idom), izerror)
      END IF
    END IF

#ifdef __COSMO__
    CALL get_runtime_timings (i_fwo_barrier)
#endif
#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_comm)
#endif

    !----------------------------------------------------------------
    ! .. Distribute some infos which are needed to organize
    !    the synchroneous and asynchroneous output:

    IF (num_radar > 1) THEN
      IF (num_radario == 0) THEN

        !----------------------------------------------------------------
        ! .. Distribute the info, for which stations the feedback files
        !    have been created by any of the processors:

        IF (lreadmeta_from_netcdf .AND. lfdbk_output) THEN

          lsendbuf(:) = .FALSE.
          DO ista=1,nradsta
            lsendbuf(ista) = rs_meta(ista)%lfdbkfile_exist
          END DO
          CALL global_values_radar (lsendbuf(1:nradsta), nradsta, 'OR', icomm_radar,    &
               -1, yzerrmsg, mpierror)
          DO ista=1,nradsta
            rs_meta(ista)%lfdbkfile_exist = lsendbuf(ista)
          END DO

        END IF

        !----------------------------------------------------------------
        ! .. Distribute the info, for which stations the geometry files
        !    have been written by any of the processors:

        IF (lout_geom .AND. .NOT.lonline .AND. .NOT.geom_all_written(idom)) THEN

          CALL global_values_radar (geom_written(1:nradsta,idom), nradsta, 'OR', icomm_radar,    &
                                   -1, yzerrmsg, mpierror)
          geom_all_written(idom) = ALL(geom_written(1:nradsta,idom))

        END IF

      ELSE  ! num_radario /= 0

        IF (lradario_pe_dom(idom)) THEN

          !----------------------------------------------------------------
          ! .. Distribute the info, for which stations the feedback files
          !    have been created by any of the output processors:

          IF (lreadmeta_from_netcdf .AND. lfdbk_output) THEN

            lsendbuf(:) = .FALSE.
            DO ista=1,nradsta
              lsendbuf(ista) = rs_meta(ista)%lfdbkfile_exist
            END DO
            CALL global_values_radar (lsendbuf(1:nradsta), nradsta, 'OR', icomm_radario_dom(idom),    &
                                      -1, yzerrmsg, mpierror)
            DO ista=1,nradsta
              rs_meta(ista)%lfdbkfile_exist = lsendbuf(ista)
            END DO

          END IF

          !----------------------------------------------------------------
          ! .. Distribute the info, for which stations the geometry files
          !    have been written by any of the processors:

          IF (lout_geom .AND. .NOT.lonline .AND. .NOT.geom_all_written(idom)) THEN

            CALL global_values_radar (geom_written(1:nradsta,idom), nradsta, 'OR', icomm_radario_dom(idom),    &
                                      -1, yzerrmsg, mpierror)
            geom_all_written(idom) = ALL(geom_written(1:nradsta,idom))

          END IF

        END IF   ! lradario_pe_dom(idom)

      END IF     ! num_radario
    END IF       ! num_radar

    !----------------------------------------------------------------

#ifdef __COSMO__
    CALL get_runtime_timings (i_fwo_comm)
#endif
#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_out)
#endif

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  CONTAINS

    SUBROUTINE dealloc_pointers (ista)

      INTEGER, INTENT(in) :: ista

      IF (ASSOCIATED(rs_data_all(ista)%radwind_mod))    DEALLOCATE(rs_data_all(ista)%radwind_mod)
      IF (ASSOCIATED(rs_data_all(ista)%vt_mod))         DEALLOCATE(rs_data_all(ista)%vt_mod)
      IF (ASSOCIATED(rs_data_all(ista)%zh_radar_mod))   DEALLOCATE(rs_data_all(ista)%zh_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%ah_radar_mod))   DEALLOCATE(rs_data_all(ista)%ah_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%zv_radar_mod))   DEALLOCATE(rs_data_all(ista)%zv_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%rrhv_radar_mod)) DEALLOCATE(rs_data_all(ista)%rrhv_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%irhv_radar_mod)) DEALLOCATE(rs_data_all(ista)%irhv_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%kdp_radar_mod))  DEALLOCATE(rs_data_all(ista)%kdp_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%adp_radar_mod))  DEALLOCATE(rs_data_all(ista)%adp_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%zvh_radar_mod))  DEALLOCATE(rs_data_all(ista)%zvh_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%radpos_all))     DEALLOCATE(rs_data_all(ista)%radpos_all)
      IF (ASSOCIATED(rs_data_all(ista)%hl_loc))         DEALLOCATE(rs_data_all(ista)%hl_loc)
      IF (ASSOCIATED(rs_data_all(ista)%el_loc))         DEALLOCATE(rs_data_all(ista)%el_loc)
      IF (ASSOCIATED(rs_data_all(ista)%s_loc))          DEALLOCATE(rs_data_all(ista)%s_loc)

    END SUBROUTINE dealloc_pointers

    SUBROUTINE nullify_pointers (ista)

      INTEGER, INTENT(in) :: ista

      IF (ASSOCIATED(rs_data_all(ista)%radwind_mod))    NULLIFY(rs_data_all(ista)%radwind_mod)
      IF (ASSOCIATED(rs_data_all(ista)%vt_mod))         NULLIFY(rs_data_all(ista)%vt_mod)
      IF (ASSOCIATED(rs_data_all(ista)%zh_radar_mod))   NULLIFY(rs_data_all(ista)%zh_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%ah_radar_mod))   NULLIFY(rs_data_all(ista)%ah_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%zv_radar_mod))   NULLIFY(rs_data_all(ista)%zv_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%rrhv_radar_mod)) NULLIFY(rs_data_all(ista)%rrhv_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%irhv_radar_mod)) NULLIFY(rs_data_all(ista)%irhv_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%kdp_radar_mod))  NULLIFY(rs_data_all(ista)%kdp_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%adp_radar_mod))  NULLIFY(rs_data_all(ista)%adp_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%zvh_radar_mod))  NULLIFY(rs_data_all(ista)%zvh_radar_mod)
      IF (ASSOCIATED(rs_data_all(ista)%radpos_all))     NULLIFY(rs_data_all(ista)%radpos_all)
      IF (ASSOCIATED(rs_data_all(ista)%hl_loc))         NULLIFY(rs_data_all(ista)%hl_loc)
      IF (ASSOCIATED(rs_data_all(ista)%el_loc))         NULLIFY(rs_data_all(ista)%el_loc)
      IF (ASSOCIATED(rs_data_all(ista)%s_loc))          NULLIFY(rs_data_all(ista)%s_loc)

    END SUBROUTINE nullify_pointers

  END SUBROUTINE output_radar

  !==============================================================================
  !+ Module procedure in radar_src for the output in case of lsmooth = .TRUE.
  !------------------------------------------------------------------------------

  SUBROUTINE output_radar_smth(lcalc,time_mod)

    !------------------------------------------------------------------------------
    !
    ! Description: Gathers the radar data on the radar grid from each radar
    !              from all PEs to a single output PE per radar and writes
    !              the results to disc. Currently this is done as
    !              ASCII, ASCII-gzip, Fortran binary or NetCDF (internally compressed),
    !              as well as netcdf feedback files for data assimilation.
    !
    !              Has to be called on icomm_radar_dom(idom), i.e., all compute-PEs
    !              and the radar-IO-PEs allocated for the actual model domain idom.
    !              if called also on other PEs, it will hang up!
    !
    !              THIS PROCEDURE IS FOR SIMULATED RADAR DATA INCLUDING SPATIAL SMOOTHING!
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                        :: lcalc(nradsta_max)
    REAL (KIND=dp), INTENT(IN)     :: time_mod       ! seconds since model start


    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'output_radar_smth'
    CHARACTER (LEN=80)  :: yzerrmsg

    INTEGER :: i, ista, ista_par, my_ista, ipe_rad, nsmth_ista,  &
               mpierror, izerror, nradsta_per_pe, num_io_loc,    &
               ista_this_pe_counter

    INTEGER, ALLOCATABLE :: my_ista_vec(:)

    TYPE(radar_data_type)    :: rs_data_all(nradsta_max)

    INTEGER :: ipos_all(num_radar_dom(idom)), idims_all(num_radar_dom(idom),nradsta_max), &
                message_tag, mpi_sendrequests(nradsta_max*20), &
                istatus(MPI_STATUS_SIZE)

    ! Flag for parallel output of each radar on different PEs. This
    ! flag indicates whether this PE has some radar data to output.
    LOGICAL              :: output_my_radar
    LOGICAL, ALLOCATABLE :: output_my_radar_vec(:)

    ! Flag if geometry files (hrsim, ersim) have been written to a file:
    LOGICAL, SAVE :: geom_written(nradsta_max,ndoms_max) = .FALSE., &
                     geom_all_written(ndoms_max) = .FALSE.

    ! Buffer for MPI_ALLREDUCE:
    LOGICAL :: lsendbuf(nradsta_max)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE output_radar_smth
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    IF (.NOT. (loutradwind .OR. loutdbz)) THEN
      WRITE (*,*) 'MESSAGE: '//TRIM(yzroutine)// &
           'No output of simulated radar data since loutradwind = loutdbz = .FALSE.'
      RETURN
    END IF

    ipos_all(:) = 0

    ! Distribute the radar data of different radars to all output PEs:
    !=================================================================

    IF (num_radario > 0) THEN
      ! Output on asynchoneous IO-PEs:
      num_io_loc = num_radario_dom(idom)
      IF (num_io_loc > 0 .AND. num_io_loc <= num_radario) THEN
        nradsta_per_pe = MAX( (nradsta-1)/num_io_loc+1 , 1)
        ALLOCATE (my_ista_vec(nradsta_per_pe))
        ALLOCATE (output_my_radar_vec(nradsta_per_pe))
        my_ista_vec(:) = -1
        output_my_radar_vec = .FALSE.
      END IF
      ! Init of counter for no. of stations per output PE (asynchr. IO):
      ista_this_pe_counter = 0
    ELSE
      ! Output on the compute-PEs:
      num_io_loc = num_radar
    END IF

#ifdef __COSMO__
    CALL get_runtime_timings (i_fwo_out)
#endif
#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_barrier)
#endif

    ! This barrier synchronizes output- and compute-PEs in case of asynchroneous output,
    !  and the  global_values_radar (geom_written(1:nradsta,idom), ...) back-propagates
    !  the info on written geometry infos from the asyn. IO PEs to all radar PEs of this domain:
    IF (num_radar > 1 .AND. num_radario > 0) THEN
      IF (lradar_pe_dom(idom)) THEN
        CALL MPI_BARRIER(icomm_radar_dom(idom), izerror)
        CALL global_values_radar (geom_written(1:nradsta,idom), nradsta, 'OR', icomm_radar_dom(idom),    &
                                  -1, yzerrmsg, mpierror)
      END IF
    END IF

#ifdef __COSMO__
    CALL get_runtime_timings (i_fwo_barrier)
#endif

    IF (num_radar > 1) THEN
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_comm)
#endif
      DO ista = 1, nradsta
        ! Prepare below calls to gatherv_values_radar(). Contains a global communication,
        !   so it is a synchro point for icomm_radar_dom(idom):
        CALL get_idims_all_onestation (rs_data(ista)%ind_intp_smth(:,2), num_radar_dom(idom), &
             icomm_radar_dom(idom), idims_all(:,ista), yzerrmsg, mpierror)
      END DO
      ! Prepare the non-blocking comm below:
      message_tag = 0
      mpi_sendrequests(:) = MPI_REQUEST_NULL
    END IF

    ! loop over num_compute_fwo-sized chunks of radar stations:
    chunkloop: DO ista_par = 1, nradsta, num_io_loc

#ifdef __ICON__
        CALL get_runtime_timings (i_fwo_comm)
#endif

      output_my_radar = .FALSE.

      ! .. Distinction between single proc and parallel runs:
      IF (num_radar > 1) THEN

        IF (num_radario == 0) THEN
          message_tag = 0
          mpi_sendrequests(:) = MPI_REQUEST_NULL
        END IF

        istaloop: DO ista = ista_par, MIN(ista_par+num_io_loc-1, nradsta)

          IF (.NOT.lcalc(ista)) CYCLE istaloop

          ! determine the number of the PE for output of the radar with no. ista
          !  in the combined communicator icomm_radar:
          ipe_rad = MOD(ista-1, num_io_loc) + radario_master_dom(idom)

          IF (my_radar_id_dom(idom) == ipe_rad .AND. num_radario_dom(idom) > 0 .AND. .NOT. lradario_pe_dom(idom)) THEN
            WRITE (*,'(a,i4,a)') 'ERROR '//TRIM(yzroutine)//' for asynchroneous radar output: '// &
                 'output assigned to this PE (my_radar_id=', my_radar_id, ') while it is not an output PE!'
          END IF

          ! .. Gather the values for this radar on the PE with no. ipe_rad:
          IF ( (lout_geom .AND. (lonline .OR. .NOT.geom_written(ista,idom))) .OR. lreadmeta_from_netcdf ) THEN

            ! Gather local ray height (for testing output purposes):
            message_tag = message_tag + 1
            CALL gatherv_values_radar(rs_data(ista)%hl_loc, rs_data_all(ista)%hl_loc, ipos_all, nsmth_ista, &
                 num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                 lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))

          END IF

          IF ( (lout_geom .AND. (lonline .OR. .NOT.geom_written(ista,idom))) .OR. lreadmeta_from_netcdf &
                                                                        .OR. (loutradwind .AND. lfall)) THEN

            ! Gather local elevation angle:
            message_tag = message_tag + 1
            CALL gatherv_values_radar(rs_data(ista)%el_loc, rs_data_all(ista)%el_loc, ipos_all, nsmth_ista, &
                 num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                 lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))

          END IF

          IF ( (lout_geom .OR. lreadmeta_from_netcdf) .AND. lonline ) THEN

            ! Gather arc distance from radar bin to radar station at MSL height:
            message_tag = message_tag + 1
            CALL gatherv_values_radar(rs_data(ista)%s_loc, rs_data_all(ista)%s_loc, ipos_all, nsmth_ista, &
                 num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                 lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))

          END IF

          IF(loutradwind) THEN

            ! Gather radial winds
            message_tag = message_tag + 1
            CALL gatherv_values_radar(rs_data(ista)%radwind_mod_smth, rs_data_all(ista)%radwind_mod_smth, ipos_all, nsmth_ista, &
                 num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                 lnonblocking=lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))
            ! Deallocation right after last use to save memory:
            IF (ASSOCIATED(rs_data(ista)%radwind_mod_smth) .AND. .NOT. lcomm_nonblocking_output) &
                 DEALLOCATE(rs_data(ista)%radwind_mod_smth)

            IF(lfall) THEN
              ! Gather fall speed
              message_tag = message_tag + 1
              CALL gatherv_values_radar(rs_data(ista)%vt_mod_smth, rs_data_all(ista)%vt_mod_smth, ipos_all, nsmth_ista, &
                   num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
                   lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))
              ! Deallocation right after last use to save memory:
              IF (ASSOCIATED(rs_data(ista)%vt_mod_smth) .AND. .NOT. lcomm_nonblocking_output) &
                   DEALLOCATE(rs_data(ista)%vt_mod_smth)

            END IF

          END IF

          IF(loutdbz .OR. (loutradwind .AND. (lweightdbz .OR. lmds_vr))) THEN
            ! Gather reflectivities
            message_tag = message_tag + 1
            CALL gatherv_values_radar(&
                 rs_data(ista)%zh_radar_mod_smth, rs_data_all(ista)%zh_radar_mod_smth, &
                 ipos_all, nsmth_ista, num_radar_dom(idom), ipe_rad, &
                 my_radar_id_dom(idom), icomm_radar_dom(idom), &
                 yzerrmsg, mpierror, &
                 lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, &
                 isrequest=mpi_sendrequests(message_tag), &
                 idims_all_in=idims_all(:,ista))
            ! Deallocation right after last use to save memory:
            IF (ASSOCIATED(rs_data(ista)%zh_radar_mod_smth) .AND. &
                .NOT. lcomm_nonblocking_output) &
                 DEALLOCATE(rs_data(ista)%zh_radar_mod_smth)

            IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%zv_radar_mod_smth, rs_data_all(ista)%zv_radar_mod_smth, &
                   ipos_all, nsmth_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))
              IF (ASSOCIATED(rs_data(ista)%zv_radar_mod_smth) .AND. &
                  .NOT. lcomm_nonblocking_output) &
                   DEALLOCATE(rs_data(ista)%zv_radar_mod_smth)

              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%rrhv_radar_mod_smth, rs_data_all(ista)%rrhv_radar_mod_smth, &
                   ipos_all, nsmth_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))
              IF (ASSOCIATED(rs_data(ista)%rrhv_radar_mod_smth) .AND. &
                  .NOT. lcomm_nonblocking_output) &
                   DEALLOCATE(rs_data(ista)%rrhv_radar_mod_smth)

              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%irhv_radar_mod_smth, rs_data_all(ista)%irhv_radar_mod_smth, &
                   ipos_all, nsmth_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))
              IF (ASSOCIATED(rs_data(ista)%irhv_radar_mod_smth) .AND. &
                  .NOT. lcomm_nonblocking_output) &
                   DEALLOCATE(rs_data(ista)%irhv_radar_mod_smth)

              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%kdp_radar_mod_smth, rs_data_all(ista)%kdp_radar_mod_smth, &
                   ipos_all, nsmth_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))
              IF (ASSOCIATED(rs_data(ista)%kdp_radar_mod_smth) .AND. &
                  .NOT. lcomm_nonblocking_output) &
                   DEALLOCATE(rs_data(ista)%kdp_radar_mod_smth)

              IF (loutpolall) THEN
                message_tag = message_tag + 1
                CALL gatherv_values_radar(&
                     rs_data(ista)%zvh_radar_mod_smth, rs_data_all(ista)%zvh_radar_mod_smth, &
                     ipos_all, nsmth_ista, num_radar_dom(idom), ipe_rad, &
                     my_radar_id_dom(idom), icomm_radar_dom(idom), &
                     yzerrmsg, mpierror, &
                     lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, &
                     isrequest=mpi_sendrequests(message_tag), &
                     idims_all_in=idims_all(:,ista))
                IF (ASSOCIATED(rs_data(ista)%zvh_radar_mod_smth) .AND. &
                    .NOT. lcomm_nonblocking_output) &
                     DEALLOCATE(rs_data(ista)%zvh_radar_mod_smth)
              END IF
            END IF

            IF (lextdbz .AND. &
                (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
              ! Gather extinction coefficients
              message_tag = message_tag + 1
              CALL gatherv_values_radar(&
                   rs_data(ista)%ah_radar_mod_smth, rs_data_all(ista)%ah_radar_mod_smth, &
                   ipos_all, nsmth_ista, num_radar_dom(idom), ipe_rad, &
                   my_radar_id_dom(idom), icomm_radar_dom(idom), &
                   yzerrmsg, mpierror, &
                   lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, &
                   isrequest=mpi_sendrequests(message_tag), &
                   idims_all_in=idims_all(:,ista))
              IF (ASSOCIATED(rs_data(ista)%ah_radar_mod_smth) .AND. &
                  .NOT. lcomm_nonblocking_output) &
                   DEALLOCATE(rs_data(ista)%ah_radar_mod_smth)

              IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
                message_tag = message_tag + 1
                CALL gatherv_values_radar(&
                     rs_data(ista)%adp_radar_mod_smth, rs_data_all(ista)%adp_radar_mod_smth, &
                     ipos_all, nsmth_ista, num_radar_dom(idom), ipe_rad, &
                     my_radar_id_dom(idom), icomm_radar_dom(idom), &
                     yzerrmsg, mpierror, &
                     lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, &
                     isrequest=mpi_sendrequests(message_tag), &
                     idims_all_in=idims_all(:,ista))
                IF (ASSOCIATED(rs_data(ista)%adp_radar_mod_smth) .AND. &
                    .NOT. lcomm_nonblocking_output) &
                     DEALLOCATE(rs_data(ista)%adp_radar_mod_smth)
              END IF
            ENDIF
          ENDIF


          ! Gather linear spatial Index:
          message_tag = message_tag + 1
          CALL gatherv_values_radar(rs_data(ista)%ind_intp_smth(:,2), rs_data_all(ista)%radpos_all_smth, ipos_all, nsmth_ista, &
               num_radar_dom(idom), ipe_rad, my_radar_id_dom(idom), icomm_radar_dom(idom), yzerrmsg, mpierror, &
               lnonblocking = lcomm_nonblocking_output, message_tag=message_tag, isrequest=mpi_sendrequests(message_tag), &
               idims_all_in=idims_all(:,ista))

          ! keep track of which radar is output on which processor:
          IF (my_radar_id_dom(idom) == ipe_rad) THEN

            IF (ldebug_radsim) THEN
              WRITE (*,'(a,i0,a,i0,a,i0)') TRIM(yzroutine)//': distribution of stations over PEs: ista = ', &
                   ista, ' ipe_rad = ', ipe_rad, ' my_radar_id_dom(idom) = ', my_radar_id_dom(idom)
            END IF

            rs_data_all(ista)%nsmth = nsmth_ista

            ! Book-keeping which ista shall be output on this processor:
            IF (num_radario > 0) THEN
              ! For asynchroneous output, it is the round robin collection one or several istas,
              !  to be output later in an extra chunk loop:
              ista_this_pe_counter = ista_this_pe_counter + 1
              my_ista_vec(ista_this_pe_counter) = ista
              output_my_radar_vec(ista_this_pe_counter) = lcalc(ista)
            ELSE
              ! For synchroneous output, it is just the ista from this chunk cycle,
              !  to be output directly within this chunk loop.
              my_ista = ista
              output_my_radar = lcalc(ista)
            END IF

          END IF

        END DO istaloop

      ELSE   ! num_radar == 1

        ! local and global radar data structure are identical
        ! No addinitional memory required:

        IF (.NOT.lcalc(ista_par)) CYCLE chunkloop

        ! Book-keeping which ista shall be output on this processor in this chunk-cylce:
        my_ista = ista_par
        output_my_radar = .TRUE.

        rs_data_all(my_ista)%nsmth = SIZE(rs_data(my_ista)%ind_intp_smth,1)

        CALL nullify_pointers_smth (my_ista)

        rs_data_all(my_ista)%radpos_all_smth     => rs_data(my_ista)%ind_intp_smth(:,2)
        rs_data_all(my_ista)%vt_mod_smth         => rs_data(my_ista)%vt_mod_smth
        rs_data_all(my_ista)%radwind_mod_smth    => rs_data(my_ista)%radwind_mod_smth
        rs_data_all(my_ista)%zh_radar_mod_smth   => rs_data(my_ista)%zh_radar_mod_smth
        rs_data_all(my_ista)%ah_radar_mod_smth   => rs_data(my_ista)%ah_radar_mod_smth
        IF (dbz_meta(my_ista)%itype_refl == 1 .OR. dbz_meta(my_ista)%itype_refl > 4) THEN
          rs_data_all(my_ista)%zv_radar_mod_smth   => rs_data(my_ista)%zv_radar_mod_smth
          rs_data_all(my_ista)%rrhv_radar_mod_smth => rs_data(my_ista)%rrhv_radar_mod_smth
          rs_data_all(my_ista)%irhv_radar_mod_smth => rs_data(my_ista)%irhv_radar_mod_smth
          rs_data_all(my_ista)%kdp_radar_mod_smth  => rs_data(my_ista)%kdp_radar_mod_smth
          rs_data_all(my_ista)%adp_radar_mod_smth  => rs_data(my_ista)%adp_radar_mod_smth
          rs_data_all(my_ista)%zvh_radar_mod_smth  => rs_data(my_ista)%zvh_radar_mod_smth
        END IF
        rs_data_all(my_ista)%hl_loc              => rs_data(my_ista)%hl_loc
        rs_data_all(my_ista)%el_loc              => rs_data(my_ista)%el_loc
        rs_data_all(my_ista)%s_loc               => rs_data(my_ista)%s_loc

      END IF ! num_radar > 1


      IF (num_radar > 1 .AND. num_radario == 0 .AND. lcomm_nonblocking_output) THEN
        ! Wait for the non-blocking send requests to complete:
        DO i = 1, message_tag
          CALL MPI_WAIT(mpi_sendrequests(i), istatus, mpierror)
          mpi_sendrequests(i) = MPI_REQUEST_NULL
        END DO
        ! Deallocate the data pointers only after completion of the send-jobs:
        DO ista = ista_par, MIN(ista_par+num_io_loc-1, nradsta)
          IF (ASSOCIATED(rs_data(ista)%radwind_mod_smth))    DEALLOCATE(rs_data(ista)%radwind_mod_smth)
          IF (ASSOCIATED(rs_data(ista)%vt_mod_smth))         DEALLOCATE(rs_data(ista)%vt_mod_smth)
          IF (ASSOCIATED(rs_data(ista)%zh_radar_mod_smth))   DEALLOCATE(rs_data(ista)%zh_radar_mod_smth)
          IF (ASSOCIATED(rs_data(ista)%ah_radar_mod_smth))   DEALLOCATE(rs_data(ista)%ah_radar_mod_smth)
          IF (ASSOCIATED(rs_data(ista)%zv_radar_mod_smth))   DEALLOCATE(rs_data(ista)%zv_radar_mod_smth)
          IF (ASSOCIATED(rs_data(ista)%rrhv_radar_mod_smth)) DEALLOCATE(rs_data(ista)%rrhv_radar_mod_smth)
          IF (ASSOCIATED(rs_data(ista)%irhv_radar_mod_smth)) DEALLOCATE(rs_data(ista)%irhv_radar_mod_smth)
          IF (ASSOCIATED(rs_data(ista)%kdp_radar_mod_smth))  DEALLOCATE(rs_data(ista)%kdp_radar_mod_smth)
          IF (ASSOCIATED(rs_data(ista)%adp_radar_mod_smth))  DEALLOCATE(rs_data(ista)%adp_radar_mod_smth)
          IF (ASSOCIATED(rs_data(ista)%zvh_radar_mod_smth))  DEALLOCATE(rs_data(ista)%zvh_radar_mod_smth)
        END DO
      END IF


#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_comm)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_out)
#endif

      IF (output_my_radar .AND. (num_radario == 0 .OR. num_radar == 1)) THEN

        ! .. Synchroneous output on the compute-PEs (in this case equal to the radar-PEs)

        ista = my_ista

        ! .. I'm the output_PE and have some data: output the data of this station:
        !     ( output_my_ista_smth() contains no further MPI communication )
        CALL output_my_ista_smth (time_mod, ista, &
             rs_data_all(ista)%nsmth, rs_data_all(ista)%radpos_all_smth, &
             rs_data_all(ista)%vt_mod_smth, rs_data_all(ista)%radwind_mod_smth, &
             rs_data_all(ista)%zh_radar_mod_smth, rs_data_all(ista)%ah_radar_mod_smth, &
             rs_data_all(ista)%zv_radar_mod_smth, &
             rs_data_all(ista)%rrhv_radar_mod_smth, rs_data_all(ista)%irhv_radar_mod_smth, &
             rs_data_all(ista)%kdp_radar_mod_smth, rs_data_all(ista)%adp_radar_mod_smth, &
             rs_data_all(ista)%zvh_radar_mod_smth, &
             rs_data_all(ista)%hl_loc, rs_data_all(ista)%el_loc, rs_data_all(ista)%s_loc, &
             geom_written(ista,idom))


        IF ( num_radar > 1 ) THEN

          ! .. In this case, the memory has been allocated to the pointers
          !    itself by gatherv_values_radar(). Therefore deallocate it.
          CALL dealloc_pointers_smth (ista)

        ELSE

          ! .. In this case, the pointers simply point to another
          !    pointers, to whom the memory has been allocated. Therefore,
          !    simply break the connection by NULLIFYING. The memory
          !    of the original pointers will be deallocated at the end
          !    of organize_radar().
          CALL nullify_pointers_smth (ista)

        END IF

      END IF ! output_my_radar = .true.

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_out)
#endif

    END DO chunkloop ! loop over (num_radar packages of) stations

    !----------------------------------------------------------------

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_comm)
#endif

    IF (num_radar > 1 .AND. num_radario > 0 .AND. lcomm_nonblocking_output) THEN
      ! Wait for the non-blocking send requests to complete:
      DO i = 1, message_tag
        CALL MPI_WAIT(mpi_sendrequests(i), istatus, mpierror)
      END DO
      ! Deallocate the data pointers only after completion of the send-jobs:
      DO ista = 1, nradsta
        IF (ASSOCIATED(rs_data(ista)%radwind_mod_smth))    DEALLOCATE(rs_data(ista)%radwind_mod_smth)
        IF (ASSOCIATED(rs_data(ista)%vt_mod_smth))         DEALLOCATE(rs_data(ista)%vt_mod_smth)
        IF (ASSOCIATED(rs_data(ista)%zh_radar_mod_smth))   DEALLOCATE(rs_data(ista)%zh_radar_mod_smth)
        IF (ASSOCIATED(rs_data(ista)%ah_radar_mod_smth))   DEALLOCATE(rs_data(ista)%ah_radar_mod_smth)
        IF (ASSOCIATED(rs_data(ista)%zv_radar_mod_smth))   DEALLOCATE(rs_data(ista)%zv_radar_mod_smth)
        IF (ASSOCIATED(rs_data(ista)%rrhv_radar_mod_smth)) DEALLOCATE(rs_data(ista)%rrhv_radar_mod_smth)
        IF (ASSOCIATED(rs_data(ista)%irhv_radar_mod_smth)) DEALLOCATE(rs_data(ista)%irhv_radar_mod_smth)
        IF (ASSOCIATED(rs_data(ista)%kdp_radar_mod_smth))  DEALLOCATE(rs_data(ista)%kdp_radar_mod_smth)
        IF (ASSOCIATED(rs_data(ista)%adp_radar_mod_smth))  DEALLOCATE(rs_data(ista)%adp_radar_mod_smth)
        IF (ASSOCIATED(rs_data(ista)%zvh_radar_mod_smth))  DEALLOCATE(rs_data(ista)%zvh_radar_mod_smth)
      END DO
    END IF

    !----------------------------------------------------------------

    IF (num_radario_dom(idom) > 0 .AND. lradario_pe_dom(idom) .AND. num_radar > 1) THEN

#ifdef __ICON__
        CALL get_runtime_timings (i_fwo_out)
#endif

      ! .. Asynchroneous output on the radario-PEs (I'm one of them):

      DO i=1, ista_this_pe_counter

        IF (output_my_radar_vec(i)) THEN

          ista = my_ista_vec(i)

          ! .. I'm the output_PE and have some data: output the data of this station:
          !     ( output_my_ista() contains no further MPI communication )

          CALL output_my_ista_smth (time_mod, ista, &
               rs_data_all(ista)%nsmth, rs_data_all(ista)%radpos_all_smth, &
               rs_data_all(ista)%vt_mod_smth, rs_data_all(ista)%radwind_mod_smth, &
               rs_data_all(ista)%zh_radar_mod_smth, rs_data_all(ista)%ah_radar_mod_smth, &
               rs_data_all(ista)%zv_radar_mod_smth, &
               rs_data_all(ista)%rrhv_radar_mod_smth, rs_data_all(ista)%irhv_radar_mod_smth, &
               rs_data_all(ista)%kdp_radar_mod_smth, rs_data_all(ista)%adp_radar_mod_smth, &
               rs_data_all(ista)%zvh_radar_mod_smth, &
               rs_data_all(ista)%hl_loc, rs_data_all(ista)%el_loc, rs_data_all(ista)%s_loc, &
               geom_written(ista,idom))

          ! .. In this case, the memory has been allocated to the pointers
          !    by gatherv_values_radar(). Therefore deallocate it:
          CALL dealloc_pointers_smth (ista)

        END IF

      END DO

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_out)
#endif

    END IF


    !----------------------------------------------------------------

    IF (ALLOCATED(output_my_radar_vec)) DEALLOCATE (output_my_radar_vec)
    IF (ALLOCATED(my_ista_vec))         DEALLOCATE (my_ista_vec)


    !----------------------------------------------------------------
    ! .. A barrier is set here, so that subsequent timing results
    !     are not biased, in case num_radar > nradsta:

#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_barrier)
#endif
    IF (num_radar > 1) THEN
      IF (num_radario == 0) THEN
        CALL MPI_BARRIER(icomm_radar, izerror)
      ELSE IF (lradario_pe_dom(idom)) THEN
        CALL MPI_BARRIER(icomm_radario_dom(idom), izerror)
      END IF
    END IF

#ifdef __COSMO__
    CALL get_runtime_timings (i_fwo_barrier)
#endif
#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_comm)
#endif

    !----------------------------------------------------------------
    ! .. Distribute some infos which are needed to organize
    !    the synchroneous and asynchroneous output:

    IF (num_radar > 1) THEN
      IF (num_radario == 0) THEN

        !----------------------------------------------------------------
        ! .. Distribute the info, for which stations the feedback files
        !    have been created by any of the processors:

        IF (lreadmeta_from_netcdf .AND. lfdbk_output) THEN

          lsendbuf(:) = .FALSE.
          DO ista=1,nradsta
            lsendbuf(ista) = rs_meta(ista)%lfdbkfile_exist
          END DO
          CALL global_values_radar (lsendbuf(1:nradsta), nradsta, 'OR', icomm_radar,    &
               -1, yzerrmsg, mpierror)
          DO ista=1,nradsta
            rs_meta(ista)%lfdbkfile_exist = lsendbuf(ista)
          END DO

        END IF

        !----------------------------------------------------------------
        ! .. Distribute the info, for which stations the geometry files
        !    have been written by any of the processors:

        IF (lout_geom .AND. .NOT.lonline .AND. .NOT.geom_all_written(idom)) THEN

          CALL global_values_radar (geom_written(1:nradsta,idom), nradsta, 'OR', icomm_radar,    &
                                   -1, yzerrmsg, mpierror)
          geom_all_written(idom) = ALL(geom_written(1:nradsta,idom))

        END IF

      ELSE  ! num_radario /= 0

        IF (lradario_pe_dom(idom)) THEN

          !----------------------------------------------------------------
          ! .. Distribute the info, for which stations the feedback files
          !    have been created by any of the output processors:

          IF (lreadmeta_from_netcdf .AND. lfdbk_output) THEN

            lsendbuf(:) = .FALSE.
            DO ista=1,nradsta
              lsendbuf(ista) = rs_meta(ista)%lfdbkfile_exist
            END DO
            CALL global_values_radar (lsendbuf(1:nradsta), nradsta, 'OR', icomm_radario_dom(idom),    &
                                      -1, yzerrmsg, mpierror)
            DO ista=1,nradsta
              rs_meta(ista)%lfdbkfile_exist = lsendbuf(ista)
            END DO

          END IF

          !----------------------------------------------------------------
          ! .. Distribute the info, for which stations the geometry files
          !    have been written by any of the processors:

          IF (lout_geom .AND. .NOT.lonline .AND. .NOT.geom_all_written(idom)) THEN

            CALL global_values_radar (geom_written(1:nradsta,idom), nradsta, 'OR', icomm_radario_dom(idom),    &
                                      -1, yzerrmsg, mpierror)
            geom_all_written(idom) = ALL(geom_written(1:nradsta,idom))

          END IF

        END IF   ! lradario_pe_dom(idom)

      END IF     ! num_radario
    END IF       ! num_radar

    !----------------------------------------------------------------

#ifdef __COSMO__
    CALL get_runtime_timings (i_fwo_comm)
#endif
#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_out)
#endif

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  CONTAINS

    SUBROUTINE dealloc_pointers_smth (ista)

      INTEGER, INTENT(in) :: ista

      IF (ASSOCIATED(rs_data_all(ista)%vt_mod_smth))         DEALLOCATE(rs_data_all(ista)%vt_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%radwind_mod_smth))    DEALLOCATE(rs_data_all(ista)%radwind_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%zh_radar_mod_smth))   DEALLOCATE(rs_data_all(ista)%zh_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%ah_radar_mod_smth))   DEALLOCATE(rs_data_all(ista)%ah_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%zv_radar_mod_smth))   DEALLOCATE(rs_data_all(ista)%zv_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%rrhv_radar_mod_smth)) DEALLOCATE(rs_data_all(ista)%rrhv_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%irhv_radar_mod_smth)) DEALLOCATE(rs_data_all(ista)%irhv_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%kdp_radar_mod_smth))  DEALLOCATE(rs_data_all(ista)%kdp_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%adp_radar_mod_smth))  DEALLOCATE(rs_data_all(ista)%adp_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%zvh_radar_mod_smth))  DEALLOCATE(rs_data_all(ista)%zvh_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%radpos_all_smth))     DEALLOCATE(rs_data_all(ista)%radpos_all_smth)
      IF (ASSOCIATED(rs_data_all(ista)%hl_loc))              DEALLOCATE(rs_data_all(ista)%hl_loc)
      IF (ASSOCIATED(rs_data_all(ista)%el_loc))              DEALLOCATE(rs_data_all(ista)%el_loc)
      IF (ASSOCIATED(rs_data_all(ista)%s_loc))               DEALLOCATE(rs_data_all(ista)%s_loc)

    END SUBROUTINE dealloc_pointers_smth

    SUBROUTINE nullify_pointers_smth (ista)

      INTEGER, INTENT(in) :: ista

      IF (ASSOCIATED(rs_data_all(ista)%vt_mod_smth))         NULLIFY(rs_data_all(ista)%vt_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%radwind_mod_smth))    NULLIFY(rs_data_all(ista)%radwind_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%zh_radar_mod_smth))   NULLIFY(rs_data_all(ista)%zh_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%ah_radar_mod_smth))   NULLIFY(rs_data_all(ista)%ah_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%zv_radar_mod_smth))   NULLIFY(rs_data_all(ista)%zv_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%rrhv_radar_mod_smth)) NULLIFY(rs_data_all(ista)%rrhv_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%irhv_radar_mod_smth)) NULLIFY(rs_data_all(ista)%irhv_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%kdp_radar_mod_smth))  NULLIFY(rs_data_all(ista)%kdp_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%adp_radar_mod_smth))  NULLIFY(rs_data_all(ista)%adp_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%zvh_radar_mod_smth))  NULLIFY(rs_data_all(ista)%zvh_radar_mod_smth)
      IF (ASSOCIATED(rs_data_all(ista)%radpos_all_smth))     NULLIFY(rs_data_all(ista)%radpos_all_smth)
      IF (ASSOCIATED(rs_data_all(ista)%hl_loc))              NULLIFY(rs_data_all(ista)%hl_loc)
      IF (ASSOCIATED(rs_data_all(ista)%el_loc))              NULLIFY(rs_data_all(ista)%el_loc)
      IF (ASSOCIATED(rs_data_all(ista)%s_loc))               NULLIFY(rs_data_all(ista)%s_loc)

    END SUBROUTINE nullify_pointers_smth

  END SUBROUTINE output_radar_smth

!================================================================================
!================================================================================

  SUBROUTINE init_cart_info (idom_in)

    !------------------------------------------------------------------------------
    !
    ! Description: The superobing grid is a rotated lat/lon grid where the rotated equator should
    !  cut approximately through the superobing domain center. Normally, this superobing
    !  domain is equal to some inner part of the computational model domain, therefore
    !  the superobing grid is model specific. The setup is defined in the model specific
    !  interface module.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in) :: idom_in

    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'init_cart_info'

    INTEGER  :: ista
    REAL(KIND=dp)         :: coslat_min, &
                             pollon,   pollat,   polgam, &
                             startlon, startlat, &
                             endlon,   endlat

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations of the grid for superobservations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    CALL get_rotlatlon_domain_for_superobing ( idom_in, &
                                               pollon,   pollat,   polgam, &
                                               startlon, startlat, &
                                               endlon,   endlat)

    DO ista = 1, nradsta

      cart_data(ista)%pollon = pollon
      cart_data(ista)%pollat = pollat
      cart_data(ista)%polgam = polgam
      cart_data(ista)%startlon = startlon
      cart_data(ista)%startlat = startlat
      cart_data(ista)%endlon = endlon
      cart_data(ista)%endlat = endlat

      ! Minimum rotated latitude cosine needed for estimation of number of quasi-cartesian grid points within
      !  one radar scan radius:
      coslat_min = MIN( COS(cart_data(ista)%startlat*degrad), COS((cart_data(ista)%endlat)*degrad))

      ! cart_data(ista)%rhori = 0.025  ! horizontal influence radius for Cressman weight function in degree
      !  cart_data%rvert = 1000  ! vertical influence radius for Cressman weight function

      cart_data(ista)%resolution = supob_cart_resolution * (0.025/2800.0) ! angular resolution of quasi-cartesian coordinates [deg]
      cart_data(ista)%dlon = cart_data(ista)%resolution
      cart_data(ista)%dlat = cart_data(ista)%resolution


      cart_data(ista)%aw_max= FLOOR(0.5*supob_azi_maxsector_vr/rs_meta(ista)%az_inc)   ! aw should not be bigger than aw_max (integer)

      cart_data(ista)%width     = 0.5*supob_ave_width                  ! half constant averaging width [m]
      cart_data(ista)%width_deg = cart_data(ista)%width*(0.025/2800.0) ! half constant averaging width [deg]

      cart_data(ista)%nra_average = CEILING(cart_data(ista)%width/rs_meta(ista)%ra_inc) ! half number of range bins for averaging [-]

      cart_data(ista)%minra_vr = supob_minrange_vr   ! minimal range to radar station to do averaging for radial wind
      cart_data(ista)%minra_z  = supob_minrange_z    ! minimal range to radar station to do averaging for reflectivity

      ! The following ni_tot, nj_tot correspond to the total number of cartesian grid points.
      cart_data(ista)%ni_tot = CEILING(((cart_data(ista)%endlon - cart_data(ista)%startlon) - &
                                       2.0_dp*cart_data(ista)%width_deg) / cart_data(ista)%dlon)
      cart_data(ista)%nj_tot = CEILING(((cart_data(ista)%endlat - cart_data(ista)%startlat) - &
                                       2.0_dp*cart_data(ista)%width_deg) / cart_data(ista)%dlat)

      ! But for one station, we only need memory for the subregion which is contained within the scan radius!
      !   - total number in i-direction within a rectangle surrounding the radar domain; the +2 is for safety
      cart_data(ista)%ni = CEILING( (2.0*rs_meta(ista)%ra_inc*rs_meta(ista)%nra) / (supob_cart_resolution*coslat_min) ) + 2
      !   - total number in j-direction within a rectangle surrounding the radar domain; the +2 is for safety
      cart_data(ista)%nj = CEILING( (2.0*rs_meta(ista)%ra_inc*rs_meta(ista)%nra) / supob_cart_resolution ) + 2

    END DO

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id


  END SUBROUTINE init_cart_info

END MODULE radar_output_driver
