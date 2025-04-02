! Source module for the radar forward operator EMVORADO
!
! ---------------------------------------------------------------
! Copyright (C) 2017-2024, DWD
! Contact information: ulrich.blahak (at) dwd.de 
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE radar_mpi_init_icon

!------------------------------------------------------------------------------
!
! Description:
!   This module provides some routines which are necessary for the initialization
!   of the MPI for the radar forward operator in the ICON model.
!
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_data, ONLY : &
       dp, wp,          &
       cmaxlen,         &
       num_compute_fwo,     & ! number of compute PEs
       num_radar,       & ! total number of radar PEs (num_compute + num_radario)
       num_radar_dom,   & ! number of radar PEs (num_compute + num_radario_dom) per radar-active model domain
       num_radario,     & ! total number of radar-IO PEs
       num_radario_dom, & ! number of radar-IO PEs per radar-active model domain
       radar_master,    & ! root-PE of the radar PE group   (in the radar comm., normally = 0)
       radario_master,  & ! root-PE of the total radario group for all domains (in the radar comm., not radario-comm.!!!)
       radario_master_dom,& ! root-PEs of the radario group for each active radar domain (in the radar_dom comm., not radario_dom-comm.!!!)
       my_cart_id_fwo,      & ! rank of this PE (=subdomain) in the cartesian communicator
       my_radar_id,     & ! rank of this PE in the radar communicator (cart+radario)
       my_radar_id_dom, & ! rank of this PE in the radar communicator (cart+radario_dom)
       my_radario_id,   & ! rank of this PE in the (asynchroneous) radario communicator icomm_radario
       my_radario_id_dom,&! rank of this PE in the (asynchroneous) radario communicator icomm_radario_dom
       icomm_cart_fwo,      & ! communicator for the virtual cartesian topology
       icomm_radar,     & ! communicator for the group of radar-IO PEs (all domains) + compute PEs
       icomm_radar_dom, & ! communicator for the group of radar-IO PEs of each domain + compute PEs
       icomm_radario,   & ! communicator for the group of radar-IO PEs (all domains)
       icomm_radario_dom,&! communicator for the group of radar-IO PEs for each domain
       lcompute_pe_fwo,     & ! indicates whether this is a compute PE or not
       lradar_pe,       & ! indicates whether this is a radar PE for any domain or not (compute or radar-IO)
       lradar_pe_dom,   & ! indicates whether this is a radar PE for a certain domain or not (compute or radar-IO)
       lradario_pe,     & ! indicates whether this is a radar-IO PE for any domain or not
       lradario_pe_dom, & ! indicates whether this is a radar-IO PE for a certain domain or not
       ndoms_radar, &  ! number of model domains with active radar simulation
       list_domains_for_model, & ! list_domains_for_model(idom_radar) returns the original model domain index of the internal
                                 ! idom_radar'th EMVORADO domain. This list is created in prep_domains_radar() below.
       list_domains_for_radar    ! list_domains_for_radar(idom_model) returns the internal EMVORADO domain index of the 
                                 ! hosting model's idom_model'th domain. This list is created in prep_domains_radar() bel

#ifndef NOMPI
  USE radar_parallel_utilities, ONLY :  &
       global_values_radar
  USE mpi
#endif

  
!==============================================================================

  IMPLICIT NONE

!==============================================================================

!==============================================================================

  PRIVATE

  PUBLIC :: init_radar_mpi

!==============================================================================
!==============================================================================
!==============================================================================
!
! Intermediate variables during ICON-devel. Have to be replaced by ICON-equivalents
!  during that process!
!
!==============================================================================
!==============================================================================
!==============================================================================

  INTEGER, PARAMETER :: iint4 = kind(1)


!==============================================================================
!==============================================================================
!==============================================================================

CONTAINS

!==============================================================================


  !------------------------------------------------------------------------------

  !============================================================================
  !
  ! Has to be called from the driving model on all PEs after domain decomposition.
  ! It ONLY returns the information necessary
  !  for the driving model and sets some global variables from radar_data.f90 internally: 
  !
  !============================================================================

  SUBROUTINE init_radar_mpi (luse_radarfwo,                        & ! INPUT
       icomm_world_icon, my_world_id_icon, nproc_icon,             & ! INPUT
       icomm_cart_icon, my_cart_id_icon, num_compute_icon,         & ! INPUT
       lcompute_pe_icon,                                           & ! INPUT
       nprocio_radar_icon, radar_master_icon, radario_master_icon, & ! INPUT
       ierror, errmsg, ldebug)

    IMPLICIT NONE

    ! INPUT parameters:
    !------------------
    INTEGER, INTENT(in) :: icomm_cart_icon, my_cart_id_icon, num_compute_icon
    INTEGER, INTENT(in) :: icomm_world_icon, my_world_id_icon, &
                           nproc_icon

    INTEGER, INTENT(in) :: radar_master_icon,   & ! Start-PE of icomm_radar in icomm_world_icon
                           radario_master_icon, & ! Start-PE of icomm_radario in icomm_world_icon
                           nprocio_radar_icon     ! Number of asynchroneous radar IO PEs provided by the hosting model

    LOGICAL, INTENT(in) :: luse_radarfwo(:), lcompute_pe_icon
    INTEGER, INTENT(out):: ierror
    CHARACTER(len=*), INTENT(inout) :: errmsg
    LOGICAL, INTENT(in) :: ldebug


    ! Local parameters:
    !------------------
    INTEGER  :: fehler, idom_radar, npes_per_dom(ndoms_radar), j, ipe_u, ipe_o
    CHARACTER(len=cmaxlen)   :: yzroutine, yzerrmsg
    INTEGER  :: key_comm, color_comm


    yzroutine(:) = ' '
    yzroutine    = 'init_radar_mpi'
    errmsg(:)  = ' '
    ierror = 0
    
    !=======================================================================================
    ! Overtake MPI communicators, PEs and cartesian grid decomposition from the COSMO-model:
    !---------------------------------------------------------------------------------------

    num_compute_fwo   = num_compute_icon
    lcompute_pe_fwo   = lcompute_pe_icon
    my_cart_id_fwo    = my_cart_id_icon
    icomm_cart_fwo    = icomm_cart_icon


    !=======================================================================================
    ! Initialize the radar IO processors:
    !---------------------------------------------------------------------------------------

    npes_per_dom(:) = -1

    IF (ANY(luse_radarfwo)) THEN

      radar_master   = 0  ! root-PE of the radar PE group   (in the radar comm., normally = 0)

#ifndef NOMPI
      IF (nproc_icon > 1) THEN

        IF (nprocio_radar_icon > 0) THEN

          ! Asynchroneous mode:
          !--------------------

          IF (radario_master_icon < radar_master_icon + num_compute_fwo) THEN
            errmsg(:) = ' '
            errmsg = 'ERROR init_radar_mpi(): nprocio_radar_icon > 0,'// &
                 ' but radario_master_icon < radar_master_icon + num_compute_fwo!'
            
            ierror = 1
            RETURN
          END IF

          IF ( (radar_master_icon   <= my_world_id_icon .AND. my_world_id_icon < radar_master_icon + num_compute_fwo)  .OR. &
               (radario_master_icon <= my_world_id_icon .AND. my_world_id_icon < radario_master_icon + nprocio_radar_icon)  ) THEN
            color_comm = 0
            key_comm   = my_world_id_icon - radar_master_icon
            lradar_pe  = .TRUE.
          ELSE
            color_comm = MPI_UNDEFINED
            key_comm   = MPI_UNDEFINED
            lradar_pe  = .FALSE.
          END IF

          CALL MPI_COMM_SPLIT( icomm_world_icon, color_comm, key_comm, icomm_radar, fehler)
          IF (fehler /= 0) THEN
            errmsg(:) = ' '
            errmsg = 'ERROR init_radar_mpi(): MPI_COMM_SPLIT for icomm_radar'
            ierror = 2
            RETURN
          ENDIF

          IF (lradar_pe) THEN
            CALL MPI_COMM_RANK (icomm_radar, my_radar_id, fehler)
            CALL MPI_COMM_SIZE (icomm_radar, num_radar  , fehler)
          ELSE
            my_radar_id = MPI_UNDEFINED
            num_radar   = -1
            icomm_radar = MPI_COMM_NULL
          END IF

          IF (my_radar_id >= num_compute_fwo .AND. lradar_pe) THEN
            color_comm     = 1
            key_comm       = my_radar_id - num_compute_fwo
            lradario_pe    = .TRUE.
          ELSE
            color_comm     = MPI_UNDEFINED
            key_comm       = MPI_UNDEFINED
            lradario_pe = .FALSE.
          END IF

          CALL MPI_COMM_SPLIT( icomm_world_icon, color_comm, key_comm, icomm_radario, fehler)
          IF (fehler /= 0) THEN
            errmsg(:) = ' '
            errmsg = 'ERROR init_radar_mpi(): MPI_COMM_SPLIT for icomm_radario'
            ierror = 3
            RETURN
          ENDIF

          IF (lradario_pe) THEN
            CALL MPI_COMM_RANK (icomm_radario, my_radario_id, fehler)
            CALL MPI_COMM_SIZE (icomm_radario, num_radario  , fehler)
          ELSE
            my_radario_id  = MPI_UNDEFINED
            num_radario    = -1
            icomm_radario  = MPI_COMM_NULL
          END IF

          ! num_radario must be known on all PEs belonging to icomm_radar!
          IF (lradar_pe) THEN
            CALL global_values_radar (num_radario, 'MAX', icomm_radar, -1, yzerrmsg, fehler)

            radario_master = num_compute_fwo     ! start-PE of icomm_radario in icomm_radar

            IF (num_radar /= num_compute_fwo + num_radario) THEN
              errmsg(:) = ' '
              WRITE(errmsg, '(a,3i5)') 'ERROR init_radar_mpi(): num_radar /= num_compute_fwo + num_radario!', &
                   num_radar, num_compute_fwo, num_radario
              ierror = 4
              RETURN
            ENDIF
          END IF


          ! The same as above, but for each radar-active model domain:
          !-----------------------------------------------------------

          IF (ndoms_radar > 1) THEN

            IF (nprocio_radar_icon < ndoms_radar) THEN

              num_radar_dom(:)      = num_radar
              num_radario_dom(:)    = num_radario
              radario_master_dom(:) = radario_master
              icomm_radar_dom(:)    = icomm_radar
              icomm_radario_dom(:)  = icomm_radario
              lradar_pe_dom(:)      = lradar_pe
              lradario_pe_dom(:)    = lradario_pe
              my_radar_id_dom(:)    = my_radar_id
              my_radario_id_dom(:)  = my_radario_id

            ELSE   !  ndoms_radar < nprocio_radar_icon

              ! .. Initialisation:
              num_radar_dom(:)      = -1  ! must be negative on output-PEs which are for a different domain! 
              num_radario_dom(:)    = -1  ! must be negative on output-PEs which are for a different domain! 
              radario_master_dom(:) = -1
              icomm_radar_dom(:)    = MPI_COMM_NULL
              icomm_radario_dom(:)  = MPI_COMM_NULL
              lradar_pe_dom(:)      = .FALSE.
              lradario_pe_dom(:)    = .FALSE.
              my_radar_id_dom(:)    = MPI_UNDEFINED
              my_radario_id_dom(:)  = MPI_UNDEFINED

              npes_per_dom(1) = (nprocio_radar_icon-1) / ndoms_radar + 1
              ipe_u = radario_master_icon
              ipe_o = ipe_u + npes_per_dom(1) - 1
              DO idom_radar = 1, ndoms_radar

                j = list_domains_for_model(idom_radar)

                IF ( (radar_master_icon   <= my_world_id_icon .AND. my_world_id_icon < radar_master_icon + num_compute_fwo)  .OR. &
                     (ipe_u                <= my_world_id_icon .AND. my_world_id_icon <= ipe_o)  ) THEN
                  color_comm = 0
                  key_comm   = my_world_id_icon - radar_master_icon
                  lradar_pe_dom(j)  = .TRUE.
                ELSE
                  color_comm = MPI_UNDEFINED
                  key_comm   = MPI_UNDEFINED
                  lradar_pe_dom(j)  = .FALSE.
                END IF

                CALL MPI_COMM_SPLIT( icomm_world_icon, color_comm, key_comm, icomm_radar_dom(j), fehler)
                IF (fehler /= 0) THEN
                  errmsg(:) = ' '
                  errmsg = 'ERROR init_radar_mpi(): MPI_COMM_SPLIT for icomm_radar_dom(j)'
                  ierror = 5
                  RETURN
                ENDIF

                IF (lradar_pe_dom(j)) THEN
                  CALL MPI_COMM_RANK (icomm_radar_dom(j), my_radar_id_dom(j), fehler)
                  CALL MPI_COMM_SIZE (icomm_radar_dom(j), num_radar_dom(j)  , fehler)
                ELSE
                  my_radar_id_dom(j) = MPI_UNDEFINED
                  num_radar_dom(j)   = -1
                  icomm_radar_dom(j) = MPI_COMM_NULL
                END IF

                IF (my_radar_id_dom(j) >= num_compute_fwo .AND. lradar_pe_dom(j)) THEN
                  color_comm     = 1
                  key_comm       = my_radar_id_dom(j) - num_compute_fwo
                  lradario_pe_dom(j)    = .TRUE.
                ELSE
                  color_comm     = MPI_UNDEFINED
                  key_comm       = MPI_UNDEFINED
                  lradario_pe_dom(j) = .FALSE.
                END IF

                CALL MPI_COMM_SPLIT( icomm_world_icon, color_comm, key_comm, icomm_radario_dom(j), fehler)
                IF (fehler /= 0) THEN
                  errmsg(:) = ' '
                  errmsg = 'ERROR init_radar_mpi(): MPI_COMM_SPLIT for icomm_radario_dom(j)'
                  ierror = 6
                  RETURN
                ENDIF

                IF (lradario_pe_dom(j)) THEN
                  CALL MPI_COMM_RANK (icomm_radario_dom(j), my_radario_id_dom(j), fehler)
                  CALL MPI_COMM_SIZE (icomm_radario_dom(j), num_radario_dom(j)  , fehler)
                ELSE
                  my_radario_id_dom(j)  = MPI_UNDEFINED
                  num_radario_dom(j)    = -1
                  icomm_radario_dom(j)  = MPI_COMM_NULL
                END IF

                ! num_radario must be known on all PEs belonging to icomm_radar!
                IF (lradar_pe_dom(j)) THEN
                  CALL global_values_radar (num_radario_dom(j), 'MAX', icomm_radar_dom(j), -1, yzerrmsg, fehler)
                  radario_master_dom(j) = num_compute_fwo   ! start-PE of icomm_radario_dom(j) in icomm_radar_dom(j)
                END IF

                ! Compute npes_per_dom(i+1) for the next domain:
                IF (idom_radar < ndoms_radar) THEN
                  npes_per_dom(idom_radar+1) = (nprocio_radar_icon-(ipe_o-radario_master_icon+1)-1) / (ndoms_radar-idom_radar) + 1
                  ipe_u = ipe_o + 1
                  ipe_o = MIN(ipe_u+npes_per_dom(idom_radar+1)-1, radario_master_icon+nprocio_radar_icon)
                END IF

              END DO

            END IF

          ELSE  !  ndoms_radar > 1

            num_radar_dom(:)      = num_radar
            num_radario_dom(:)    = num_radario
            radario_master_dom(:) = radario_master
            icomm_radar_dom(:)    = icomm_radar
            icomm_radario_dom(:)  = icomm_radario
            lradar_pe_dom(:)      = lradar_pe
            lradario_pe_dom(:)    = lradario_pe
            my_radar_id_dom(:)    = my_radar_id
            my_radario_id_dom(:)  = my_radario_id

          END IF


        ELSE   !  nprocio_radar_icon > 0

          ! Synchroneous mode:
          !-------------------

          IF (radar_master_icon <= my_world_id_icon .AND. my_world_id_icon < radar_master_icon + num_compute_fwo) THEN
            color_comm = 0
            key_comm   = my_world_id_icon
            lradar_pe  = .TRUE.
          ELSE
            color_comm = MPI_UNDEFINED
            key_comm   = MPI_UNDEFINED
            lradar_pe  = .FALSE.
          END IF

          CALL MPI_COMM_SPLIT( icomm_world_icon, color_comm, key_comm, icomm_radar, fehler)
          IF (fehler /= 0) THEN
            errmsg(:) = ' '
            errmsg = 'ERROR init_radar_mpi(): MPI_COMM_SPLIT for icomm_radar'
            ierror = 7
            RETURN
          ENDIF

          IF (lradar_pe) THEN
            CALL MPI_COMM_RANK (icomm_radar, my_radar_id, fehler)
            CALL MPI_COMM_SIZE (icomm_radar, num_radar  , fehler)
          ELSE
            my_radar_id = MPI_UNDEFINED
            num_radar   = -1
            icomm_radar = MPI_COMM_NULL
          END IF

          my_radario_id  = MPI_UNDEFINED
          num_radario    = 0
          lradario_pe    = .FALSE.
          radario_master = 0     ! start-PE of icomm_radario in icomm_radar
          icomm_radario  = MPI_COMM_NULL

          num_radar_dom(:)      = num_radar
          num_radario_dom(:)    = num_radario
          radario_master_dom(:) = radario_master
          icomm_radar_dom(:)    = icomm_radar
          icomm_radario_dom(:)  = icomm_radario
          lradar_pe_dom(:)      = lradar_pe
          lradario_pe_dom(:)    = lradario_pe
          my_radar_id_dom(:)    = my_radar_id
          my_radario_id_dom(:)  = my_radario_id

        END IF

      ELSE  ! nproc_icon > 0
#endif

        my_radar_id    = 0
        my_radario_id  = 0
        num_radar      = 1
        num_radario    = 0
        icomm_radar    = icomm_world_icon ! = icomm_cart in this case
        icomm_radario  = icomm_world_icon
        lradar_pe      = .TRUE.
        lradario_pe    = .FALSE.
        radar_master   = 0     ! start-PE of icomm_radar, normally 0
        radario_master = 0     ! start-PE of icomm_radario in icomm_radar

        num_radar_dom(:)      = num_radar
        num_radario_dom(:)    = num_radario
        radario_master_dom(:) = radario_master
        icomm_radar_dom(:)    = icomm_radar
        icomm_radario_dom(:)  = icomm_radario
        lradar_pe_dom(:)      = lradar_pe
        lradario_pe_dom(:)    = lradario_pe
        my_radar_id_dom(:)    = my_radar_id
        my_radario_id_dom(:)  = my_radario_id  

#ifndef NOMPI
      END IF
#endif

      IF (ldebug) THEN
        IF (my_radar_id == 0) THEN
          WRITE(*,'(a)') 'DEBUG '//TRIM(yzroutine)//':  '                  // &
               'my_world_id_icon, lcompute_pe_fwo, my_cart_id_fwo, '       // &
               'idom_radar, ndoms_radar, npes_per_dom, '                   // &
               'lradar_pe_dom, icomm_radar_dom, my_radar_id_dom, '         // &
               'lradario_pe_dom, icomm_radario_dom, my_radario_id_dom, '   // &
               'lradario_pe, icomm_radario, my_radario_id, '               // &
               'MPI_UNDEFINED, MPI_COMM_NULL'
        END IF
        DO idom_radar = 1, ndoms_radar
          j = list_domains_for_model(idom_radar)
          WRITE(*,'(a,1x,i6,1x,L1,4x,4i12,3(4x,L1,2i12),4x,2i12)') 'DEBUG '//TRIM(yzroutine)//':', &
               my_world_id_icon, lcompute_pe_fwo, my_cart_id_fwo, &
               idom_radar, ndoms_radar, npes_per_dom(idom_radar), &
               lradar_pe_dom(j)  , icomm_radar_dom(j)  , my_radar_id_dom(j), &
               lradario_pe_dom(j), icomm_radario_dom(j), my_radario_id_dom(j), &
               lradario_pe       , icomm_radario       , my_radario_id, &
               MPI_UNDEFINED, MPI_COMM_NULL
        END DO
      END IF

    ELSE

      my_radar_id        = my_world_id_icon
      lradar_pe          = .FALSE.
      lradario_pe        = .FALSE.

    END IF   !  ANY(luse_radarfwo)

  END SUBROUTINE init_radar_mpi

END MODULE radar_mpi_init_icon

