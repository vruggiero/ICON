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


MODULE  radar_parallel_utilities

!------------------------------------------------------------------------------
!
! Description:
!   This module provides the necessary MPI routines for the radar forward
!   Operator EMRADSCOPE. It has been derived from the original COSMO-
!   module parallel_utilities.f90. The main difference is, that the
!   MPI data types for communication are no longer subroutine arguments,
!   but are hardcoded in the overloaded subroutines. The correct MPI data types
!   are private module variables and are initialized by a call to SR init_par_utilities_radar.
!   Further, only those routines needed by the
!   radar forward operator are contained.
!
!   Routines currently contained are
!
!     - distribute_values_radar
!     - gatherv_values_radar
!     - global_values_radar
!     - gather_values_radar
!
!  MPI-routines to define MPI types for some derived types:
!     - def_mpi_radar_meta_type
!     - def_mpi_polmp_type
!     - def_mpi_dbzcalc_params_type
!     - def_mpi_compmeta_type
!     - def_mpi_voldataostream_type
!
! Method:
!   Calls to the message-passing library MPI
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_kind, ONLY : dp, sp
  
  USE radar_data, ONLY :   &
       nscanstrategies_max, & ! }
       ndatakind,           & !  }
       nel_max,             & !   some maximum dimension parameters for the radar
       nel_composite_max,   & !   some maximum dimension parameters for the radar
       ngpsm_max,           & !  }  meta data derived types
       nobstimes_max,       & ! }
       noutput_fields_max,  &
       cmaxlen, cobsflen, cvarlen

!==============================================================================

!==============================================================================

#ifndef NOMPI
  USE mpi
#endif

!==============================================================================

  IMPLICIT NONE

!==============================================================================

! include statements

#ifdef NOMPI
  INCLUDE "nompi_mpif.h"
#endif

!==============================================================================

  PUBLIC

! Private Declarations

  INTEGER, PRIVATE     ::      &
    ! datatypes for message passing
    IPU_real, IPU_integer

!==============================================================================

! Interface Blocks
INTERFACE global_values_radar
  MODULE PROCEDURE                        &
    global_vectorint,                     &
    global_vectorreal,                    &
    global_vectordouble,                  &
    global_vectorlogical,                 &
    global_array2dsingle,                 &
    global_array2ddouble,                 &
    global_array3dsingle,                 &
    global_array3ddouble,                 &
    global_int,                           &
    global_real,                          &
    global_double,                        &
    global_logical
END INTERFACE

INTERFACE distribute_values_radar
  MODULE PROCEDURE                        &
    distribute_integer,                   &
    distribute_oneinteger,                &
    distribute_idouble,                   &
    distribute_onedouble,                 &
    distribute_logical,                   &
    distribute_onelogical,                &
    distribute_character,                 &
    distribute_onecharacter
END INTERFACE

INTERFACE gatherv_values_radar
  MODULE PROCEDURE                        &
    gatherv_radar_integers,               &
    gatherv_radar_reals
END INTERFACE

INTERFACE gather_values_radar
  MODULE PROCEDURE                        &
    gather_one_int
END INTERFACE

INTERFACE get_idims_all_onestation
  MODULE PROCEDURE                        &
       get_idims_all_onestation_sp,       &
       get_idims_all_onestation_dp,       &
       get_idims_all_onestation_int
END INTERFACE

!==============================================================================

CONTAINS

!==============================================================================
!+ Initializes private variables for module parallel_utilities
!------------------------------------------------------------------------------

SUBROUTINE init_par_utilities_radar                                          &
     (imp_real, imp_int)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT(IN)  ::  &
    imp_real, imp_int               ! datatypes for message passing

! Local Variables
  INTEGER              ::  &
    izstat

!==============================================================================

  ! datatypes for message passing
  IPU_real        = imp_real   ! corresponds to KIND=dp
  IPU_integer     = imp_int    ! corresponds to standard integer

END SUBROUTINE init_par_utilities_radar

!==============================================================================

!+ Defines all subroutines for the generic routine distribute_values
!------------------------------------------------------------------------------
!
! SUBROUTINE distribute_values (buffer, ibufferlen, isender,
!                               icommunicator, ierrorcode)
!
!------------------------------------------------------------------------------
!
! Description:
!  distribute_values is a generic name for several subroutines that distribute
!  values from one processor to all others. Depending on the type of the
!  first argument, the appropriate procedure is chosen.
!
! Method:
!  With the MPI_BCAST command the buffer is distributed to all other
!  processors.
!
!------------------------------------------------------------------------------

! Following are the different subroutines

!------------------------------------------------------------------------------

!+ Subroutine for array of standard integers

SUBROUTINE distribute_integer (buffer, ibufferlen, isender,     &
                               icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
INTEGER, INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER, INTENT(OUT), OPTIONAL ::                         &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                   ::                             &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

  CALL MPI_BCAST (buffer, ibufferlen, MPI_INTEGER, isender,                 &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF

END SUBROUTINE distribute_integer

!==============================================================================

!==============================================================================

!+ Subroutine for one standard integer

SUBROUTINE distribute_oneinteger(buffer, ibufferlen, isender,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
INTEGER, INTENT(INOUT)    ::                              &
  buffer                ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER, INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                   ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

  CALL MPI_BCAST (buffer, 1, MPI_INTEGER, isender, icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF

END SUBROUTINE distribute_oneinteger

!==============================================================================

!==============================================================================

!+ Subroutine for array of KIND=dp:

SUBROUTINE distribute_idouble   (buffer, ibufferlen, isender,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
REAL    (KIND=dp),         INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER, INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                   ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

  CALL MPI_BCAST (buffer, ibufferlen, MPI_DOUBLE_PRECISION, isender,   &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF

END SUBROUTINE distribute_idouble

!==============================================================================

!==============================================================================

!+ Subroutine for one KIND=dp

SUBROUTINE distribute_onedouble (buffer, ibufferlen, isender,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
REAL    (KIND=dp),      INTENT(INOUT)    ::                              &
  buffer                ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER, INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                   ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

  CALL MPI_BCAST (buffer, ibufferlen, MPI_DOUBLE_PRECISION, isender,   &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF

END SUBROUTINE distribute_onedouble

!==============================================================================

!==============================================================================

!+ Subroutine for array of default logicals

SUBROUTINE distribute_logical   (buffer, ibufferlen, isender,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
LOGICAL,                  INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER, INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                   ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

  CALL MPI_BCAST (buffer, ibufferlen, MPI_LOGICAL, isender,   &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF

END SUBROUTINE distribute_logical

!==============================================================================

!==============================================================================

!+ Subroutine for one default logical

SUBROUTINE distribute_onelogical(buffer, ibufferlen, isender,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
LOGICAL,                  INTENT(INOUT)    ::                              &
  buffer                ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER, INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                   ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

  CALL MPI_BCAST (buffer, 1, MPI_LOGICAL, isender, icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF

END SUBROUTINE distribute_onelogical

!==============================================================================

!==============================================================================

!+ Subroutine for array of characters

SUBROUTINE distribute_character (buffer, ibufferlen, isender,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
CHARACTER (LEN=100),      INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER, INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars and Arrays
INTEGER                   ::                              &
  my_comm_id, implcode, i, j

INTEGER    ::     &  ! Standard integer
  intbuf(100)

!- End of header
!------------------------------------------------------------------------------

! T3E knows no CHAR BCAST, therefore we have do to a bit here:

  CALL MPI_COMM_RANK(icommunicator, my_comm_id, implcode)

  DO i=1,ibufferlen
    IF (my_comm_id == isender) THEN
      DO j=1,100
        intbuf(j) = ICHAR ( buffer(i)(j:j) )
      ENDDO
    ENDIF

    CALL MPI_BCAST (intbuf, 100, MPI_INTEGER, isender, icommunicator, implcode)

    IF (my_comm_id /= isender ) THEN
      DO j=1,100
        buffer(i)(j:j) = CHAR (intbuf(j) )
      ENDDO
    ENDIF
  ENDDO

! and this would be the normal way
! CALL MPI_BCAST (buffer, ibufferlen, MPI_CHARACTER, isender,   &
!                 icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF

END SUBROUTINE distribute_character

!==============================================================================

!==============================================================================

!+ Subroutine for one word of characters

SUBROUTINE distribute_onecharacter (buffer, ibufferlen, isender,  &
                                    icommunicator,  ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
CHARACTER (LEN=*),        INTENT(INOUT)    ::                              &
  buffer                ! character to be broadcasted

CHARACTER (LEN=100)  :: internal_buffer   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER, INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars and Arrays
INTEGER                   ::                              &
  my_comm_id, implcode, i, j

INTEGER    ::     &  ! Standard integer
  intbuf(100)

!- End of header
!------------------------------------------------------------------------------

! T3E knows no CHAR BCAST, therefore we have do to a bit here:

  internal_buffer(:) = ' '
  internal_buffer = buffer
  CALL MPI_COMM_RANK(icommunicator, my_comm_id, implcode)

  IF (my_comm_id == isender) THEN
    DO j=1,100
      intbuf(j) = ICHAR ( internal_buffer(j:j) )
    ENDDO
  ENDIF

  CALL MPI_BCAST (intbuf, 100, MPI_INTEGER, isender, icommunicator, implcode)

  IF (my_comm_id /= isender ) THEN
    DO j=1,100
      internal_buffer(j:j) = CHAR (intbuf(j) )
    ENDDO
  ENDIF

! and this would be the normal way
! CALL MPI_BCAST (internal_buffer, 1, MPI_CHARACTER, isender,              &
!                 icommunicator, implcode)

  buffer = internal_buffer

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF

END SUBROUTINE distribute_onecharacter

!==============================================================================
!==============================================================================

! This is a hack for distribute_values() of path variables with a length > 100 characters,
!   because distribute_values() does only transfer the first 100 characters.
SUBROUTINE distribute_path_radar ( pathstring, icomm_ )

  CHARACTER(len=*), INTENT(inout):: pathstring
  INTEGER,          INTENT(in)   :: icomm_
  INTEGER                        :: plen, i, badgelen, cs, ce, ierr

  badgelen = 100
  plen = LEN(pathstring) ! NOT len_trim(), because plen needs to be the same value on all processors!

  DO i=1, plen/badgelen+1
    cs = (i-1) * badgelen + 1
    ce = MIN(cs + badgelen - 1, plen)
    CALL distribute_values_radar(pathstring(cs:ce), 100, 0, icomm_, ierr)
  END DO

END SUBROUTINE distribute_path_radar

!==============================================================================
!==============================================================================

SUBROUTINE get_idims_all_onestation_dp (rvector_loc, npes, icomm, idims_all, yerrmsg, ierror)

  REAL(kind=dp),       INTENT(in)     :: rvector_loc(:)        ! local input vector
  INTEGER,       INTENT(in)     :: npes, icomm
  INTEGER,       INTENT(out)    :: idims_all(npes)
  CHARACTER(len=*), INTENT(out) :: yerrmsg
  INTEGER,          INTENT(out) :: ierror

  INTEGER :: idim_loc, implcode
  CHARACTER(len=250) :: zyerrmsg

  yerrmsg(:)  = ' '
  zyerrmsg(:) = ' '
  ierror = 0
  idims_all(:) = -5

  ! determine dimension of input vector
  idim_loc = SIZE(rvector_loc)

  ! collect and redistribute to all PEs in icomm:
  CALL gather_values_radar ( idim_loc, idims_all, 1, npes, -1,    &
         icomm, yerrmsg, implcode )
  IF (implcode /= 0) THEN
    ierror  = 1
    yerrmsg = "get_idims_all_onestation_dp: error in gathering the dimensions of the local vectors: " // TRIM(zyerrmsg)
    RETURN
  END IF

END SUBROUTINE get_idims_all_onestation_dp

!==============================================================================

SUBROUTINE get_idims_all_onestation_sp (rvector_loc, npes, icomm, idims_all, yerrmsg, ierror)

  REAL(kind=sp),       INTENT(in)     :: rvector_loc(:)        ! local input vector
  INTEGER,       INTENT(in)     :: npes, icomm
  INTEGER,       INTENT(out)    :: idims_all(npes)
  CHARACTER(len=*), INTENT(out) :: yerrmsg
  INTEGER,          INTENT(out) :: ierror

  INTEGER :: idim_loc, implcode
  CHARACTER(len=250) :: zyerrmsg

  yerrmsg(:)  = ' '
  zyerrmsg(:) = ' '
  ierror = 0
  idims_all(:) = -5

  ! determine dimension of input vector
  idim_loc = SIZE(rvector_loc)

  ! collect and redistribute to all PEs in icomm:
  CALL gather_values_radar ( idim_loc, idims_all, 1, npes, -1,    &
         icomm, yerrmsg, implcode )
  IF (implcode /= 0) THEN
    ierror  = 1
    yerrmsg = "get_idims_all_onestation_sp: error in gathering the dimensions of the local vectors: " // TRIM(zyerrmsg)
    RETURN
  END IF

END SUBROUTINE get_idims_all_onestation_sp

!==============================================================================

SUBROUTINE get_idims_all_onestation_int (ivector_loc, npes, icomm, idims_all, yerrmsg, ierror)

  INTEGER,       INTENT(in)     :: ivector_loc(:)        ! local input vector
  INTEGER,       INTENT(in)     :: npes, icomm
  INTEGER,       INTENT(out)    :: idims_all(npes)
  CHARACTER(len=*), INTENT(out) :: yerrmsg
  INTEGER,          INTENT(out) :: ierror

  INTEGER :: idim_loc, implcode
  CHARACTER(len=250) :: zyerrmsg

  yerrmsg(:)  = ' '
  zyerrmsg(:) = ' '
  ierror = 0
  idims_all(:) = -5

  ! determine dimension of input vector
  idim_loc = SIZE(ivector_loc)

  ! collect and redistribute to all PEs in icomm:
  CALL gather_values_radar ( idim_loc, idims_all, 1, npes, -1,    &
         icomm, zyerrmsg, implcode )
  IF (implcode /= 0) THEN
    ierror  = 1
    yerrmsg = "get_idims_all_onestation_int: error in gathering the dimensions of the local vectors: " // TRIM(zyerrmsg)
    RETURN
  END IF

END SUBROUTINE get_idims_all_onestation_int

!==============================================================================
!==============================================================================
!
!+ Defines all subroutines for the generic routine gatherv_radar_values
!------------------------------------------------------------------------------
!
! SUBROUTINE gatherv_radar_values (vector_loc, vector_all, ipos_all, npes,   &
!                           ireceiver, icomm, yerrmsg, ierror)
!
!------------------------------------------------------------------------------
!
! Description:
!  This routine gathers the values of vector_loc from all nodes to
!  the PE with ID ireceiver. The local vectors may be of different (even 0) length
!  If ireceiver < 0, vector_out is send to all PEs.
!
! Method:
!  MPI_GATHERV or MPI_ALLGATHERV
!
!------------------------------------------------------------------------------

! Following are the different subroutines

!==============================================================================
!+ Subroutine for vector of KIND=dp
SUBROUTINE gatherv_radar_reals (rvector_loc, rvector_all, ipos_all, ilen_all, &
                                npes,  ireceiver, my_id, icomm, yerrmsg, ierror, &
                                lnonblocking, message_tag, isrequest, idims_all_in)

!------------------------------------------------------------------------------
!
! NOTE  rvector_all(:) is a pointer to a vector, and MUST BE INITIALIZED BEFORE INPUT
!       e.g., nullify(rvector_all)
!
! UB: Added ilen_all to the output variables, so that the case of ilen_all=0 can
!     be treated correctly by the calling program.
!     The allocation state of the output vector rvector_all will only be touched
!     on the receiving processors!
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):

INTEGER, INTENT(IN)       ::                              &
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver,            & ! PE that receives the vector in icomm
  my_id                   ! The PE-ID of the caller in icomm

REAL (KIND=dp),           INTENT(IN)       ::                              &
  rvector_loc(:)        ! local input vector

! UB>> the SX-compiler does not compile POINTER, INTENT(...) together!
!REAL (KIND=dp),     POINTER, INTENT(out)    ::                              &
REAL (KIND=dp),     POINTER    ::                                          &
  rvector_all(:)        ! global total vector; MUST BE INITIALIZED BEFORE INPUT! e.g., nullify(rvector_all)

INTEGER, INTENT(out)    ::                                &
  ipos_all(npes)      , & ! positions in output vector
  ilen_all                ! length of global vector

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER, INTENT(OUT)      ::                              &
  ierror                  ! error status variable
INTEGER, INTENT(OUT), OPTIONAL   ::                       &
  isrequest               ! MPI request for external MPI_WAIT(isrequest) in case of nonblocking comm.

INTEGER, INTENT(in), OPTIONAL ::                       &
  idims_all_in(npes)      ! array of dimensions of the input vectors of all PEs in icomm
LOGICAL, INTENT(in), OPTIONAL :: lnonblocking
INTEGER, INTENT(in), OPTIONAL :: message_tag ! message_tag = unique tag of this message among all radar stations and data types

!------------------------------------------------------------------------------

! Local variables
!------------------------------------------------------------------------------

INTEGER                   ::                              &
     implcode            , & ! error status variable for MPI routines
     idim_loc            , & ! dimension of local vector
     ipe                     ! loop variable for pe's
INTEGER                   ::                              &
     idims_all(npes), &      ! array of dimensions of all input vectors
     itag(npes),       &      ! array of communication tags for non-blocking comm.
     istatus(MPI_STATUS_SIZE), irrequest(npes)
REAL (KIND=dp),     ALLOCATABLE        ::                              &
     zrvector_loc(:)
LOGICAL :: deall_rvec, lnonblck

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_radar_reals
!------------------------------------------------------------------------------

  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  IF (PRESENT(lnonblocking)) THEN
    lnonblck = lnonblocking
    IF (.NOT.PRESENT(message_tag) .OR. .NOT.PRESENT(isrequest) .OR. .NOT.PRESENT(idims_all_in)) THEN
      ierror  = 5
      yerrmsg = 'gatherv_radar_reals: error lnonblocking=.TRUE. but missing argument(s) ' // &
                '"message_tag" or "isrequest" or "idims_all_in"'
      RETURN
    END IF
  ELSE
    lnonblck = .FALSE.
  END IF

  IF (PRESENT(idims_all_in)) THEN
    idims_all = idims_all_in
  ELSE
    CALL get_idims_all_onestation (rvector_loc, npes, icomm, idims_all, yerrmsg, implcode)
  END IF

  ! determine dimension of input vector
  idim_loc = SIZE(rvector_loc)
  IF (idim_loc /= idims_all(my_id+1)) THEN
    PRINT *, 'ULI Scheisse!!! ', my_id, idim_loc, idims_all(my_id+1)
    STOP
  END IF

  ! determine position vector
  ipos_all(1) = 0
  DO ipe = 2 , npes
     ipos_all(ipe)  = ipos_all(ipe-1) + idims_all(ipe-1)
  END DO
  ilen_all  =  ipos_all(npes) + idims_all(npes)

  IF (ilen_all /= 0) THEN

    IF (.NOT.lnonblck .OR. ireceiver < 0) THEN
      ALLOCATE(zrvector_loc(MAX(idim_loc,1)))
      zrvector_loc = 0.0_dp
      IF (idim_loc == 0) THEN
        !print *, my_id, "zero size in gatherv_radar_values(real)"
      ELSE
        zrvector_loc = rvector_loc
      END IF
    END IF

    ! MPI-Routine
    IF (ireceiver < 0) THEN
      ! allocate output vector
      IF (ASSOCIATED(rvector_all)) DEALLOCATE(rvector_all)
      ALLOCATE(rvector_all(ilen_all))
      CALL MPI_ALLGATHERV   (zrvector_loc, idim_loc, MPI_DOUBLE_PRECISION,           &
           rvector_all, idims_all, ipos_all, MPI_DOUBLE_PRECISION,           &
           icomm, implcode)
      IF (implcode /= 0) THEN
        ierror  = 2
        yerrmsg = 'gatherv_radar_reals: Error in MPI_ALLGATHERV'
      ENDIF
    ELSEIF (ireceiver < npes) THEN
      ! allocate output vector on receiver, otherwise do not change
      ! allocation state if already allocated
      deall_rvec = .FALSE.
      IF (my_id == ireceiver) THEN
        IF (ASSOCIATED(rvector_all)) DEALLOCATE(rvector_all)
        ALLOCATE(rvector_all(ilen_all))
      ELSE
        ! really necessary???
        IF (.NOT. ASSOCIATED(rvector_all)) THEN
          ALLOCATE(rvector_all(1:1))
          deall_rvec = .TRUE.
        END IF
      END IF

      IF (lnonblck) THEN

        ! Create send tags for the non-blocking communications, unique for each
        !  send PE and message:
        DO ipe=0, npes-1
          itag(ipe+1) = ipe + message_tag*npes
        END DO

        IF (my_id == ireceiver) THEN

          ! If I'm the receiver, simply copy my own data to the receive-vector, if any:
          IF (idims_all(my_id+1) > 0) THEN
            rvector_all(ipos_all(my_id+1)+1:ipos_all(my_id+1)+idims_all(my_id+1)) = &
                 rvector_loc(:)
          END IF

          ! If I'm the receiver, open receive channels from all other PEs and wait for completion:
          DO ipe=0, npes-1
            IF (idims_all(ipe+1) > 0 .AND. ipe /= ireceiver) THEN
              CALL MPI_IRECV(rvector_all(ipos_all(ipe+1)+1:ipos_all(ipe+1)+idims_all(ipe+1)), idims_all(ipe+1), &
                   MPI_DOUBLE_PRECISION, &
                   ipe, itag(ipe+1), icomm, irrequest(ipe+1), implcode)
            END IF
          END DO
          DO ipe=0, npes-1
            IF (idims_all(ipe+1) > 0 .AND. ipe /= ireceiver) THEN
              CALL MPI_WAIT(irrequest(ipe+1), istatus, implcode)
            END IF
          END DO

        END IF

        ! If I'm not the receiver, send my data to the receiver:
        IF (my_id /= ireceiver .AND. idims_all(my_id+1) > 0) THEN
          CALL MPI_ISEND(rvector_loc, idim_loc, MPI_DOUBLE_PRECISION, &
               ireceiver, itag(my_id+1), icomm, isrequest, implcode)
        ELSE
          isrequest = MPI_REQUEST_NULL
        END IF

      ELSE

        CALL MPI_GATHERV   (zrvector_loc, idim_loc, MPI_DOUBLE_PRECISION,           &
                            rvector_all, idims_all, ipos_all, MPI_DOUBLE_PRECISION,           &
                            ireceiver, icomm, implcode)
        IF (implcode /= 0) THEN
          ierror  = 2
          yerrmsg = 'gatherv_radar_reals: Error in MPI_GATHERV'
        ENDIF

      END IF

      IF (deall_rvec) THEN
        DEALLOCATE(rvector_all)
      END IF
    ELSE
      ierror  = 3
      yerrmsg = 'gatherv_radar_reals: no valid receiver: ireceiver >= number of PEs'
    ENDIF

    IF (ALLOCATED(zrvector_loc)) DEALLOCATE(zrvector_loc)

    ipos_all(:) = ipos_all(:) + 1

  ELSE

    IF (my_id == ireceiver .AND. .NOT.ASSOCIATED(rvector_all)) THEN
      ALLOCATE(rvector_all(1:1))
      rvector_all = 0.0_dp
    END IF

  END IF

END SUBROUTINE gatherv_radar_reals

!==============================================================================

!+ Subroutine for vector of integers
SUBROUTINE gatherv_radar_integers (ivector_loc, ivector_all, ipos_all, ilen_all, &
                                   npes,  ireceiver, my_id, icomm, yerrmsg, ierror, &
                                   lnonblocking, message_tag, isrequest, idims_all_in)

!------------------------------------------------------------------------------
!
! NOTE  ivector_all(:) is a pointer to a vector, and MUST BE INITIALIZED BEFORE INPUT
!       e.g., nullify(ivector_all)
!
! UB: Added ilen_all to the output variables, so that the case of ilen_all=0 can
!     be treated correctly by the calling program.
!
!------------------------------------------------------------------------------
!
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):

INTEGER, INTENT(IN)       ::                              &
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver,            & ! PE that receives the vector in icomm
  my_id                   ! The PE-ID of the caller in icomm

INTEGER,       INTENT(IN)       ::                        &
  ivector_loc(:)        ! local input vector

! UB>> the SX-compiler does not compile POINTER, INTENT(...) together!
!INTEGER, pointer, intent(out)    ::                              &
INTEGER, POINTER    ::                                    &
  ivector_all(:)        ! global total vector; MUST BE INITIALIZED BEFORE INPUT! e.g., nullify(ivector_all)

INTEGER, INTENT(out)    ::                                &
  ipos_all(npes)      , & ! positions in output vector
  ilen_all                ! length of global vector

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER, INTENT(OUT)      ::                              &
  ierror                  ! error status variable
INTEGER, INTENT(OUT), OPTIONAL   ::                       &
  isrequest               ! MPI request for external MPI_WAIT(isrequest) in case of nonblocking comm.

INTEGER, INTENT(in), OPTIONAL ::                       &
  idims_all_in(npes)      ! array of dimensions of the input vectors of all PEs in icomm
LOGICAL, INTENT(in), OPTIONAL :: lnonblocking
INTEGER, INTENT(in), OPTIONAL :: message_tag ! message_tag = unique tag of this message among all radar stations and data types

!------------------------------------------------------------------------------

! Local variables
!------------------------------------------------------------------------------

INTEGER                   ::                              &
  implcode            , & ! error status variable for MPI routines
  idim_loc            , & ! dimension of local vector
  ipe                      ! loop variable for pe's
INTEGER                   ::                              &
  idims_all(npes), &         ! array of dimensions of all input vectors
  itag(npes),       &      ! array of communication tags for non-blocking comm.
  istatus(MPI_STATUS_SIZE), irrequest(npes)
INTEGER, ALLOCATABLE      ::                              &
  zivector_loc(:)
LOGICAL :: deall_ivec, lnonblck

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_radar_integers
!------------------------------------------------------------------------------

  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  IF (PRESENT(lnonblocking)) THEN
    lnonblck = lnonblocking
    IF (.NOT.PRESENT(message_tag) .OR. .NOT.PRESENT(isrequest) .OR. .NOT.PRESENT(idims_all_in)) THEN
      ierror  = 5
      yerrmsg = 'gatherv_radar_reals: error lnonblocking=.TRUE. but missing argument(s) ' // &
                '"message_tag" or "isrequest" or "idims_all_in"'
      RETURN
    END IF
  ELSE
    lnonblck = .FALSE.
  END IF

  IF (PRESENT(idims_all_in)) THEN
    idims_all = idims_all_in
  ELSE
    CALL get_idims_all_onestation (ivector_loc, npes, icomm, idims_all, yerrmsg, implcode)
  END IF

  ! determine dimension of input vector
  idim_loc = SIZE(ivector_loc)
  IF (idim_loc /= idims_all(my_id+1)) THEN
    PRINT *, 'ULI Scheisse!!! ', my_id, idim_loc, idims_all(my_id+1)
    STOP
  END IF

  ! determine position vector
  ipos_all(1) = 0
  DO ipe = 2 , npes
    ipos_all(ipe)  = ipos_all(ipe-1) + idims_all(ipe-1)
  END DO
  ilen_all  =  ipos_all(npes) + idims_all(npes)

  IF (ilen_all /= 0) THEN

    IF (.NOT.lnonblck .OR. ireceiver < 0) THEN
      ALLOCATE(zivector_loc(MAX(idim_loc,1)))
      zivector_loc = 0.0_dp
      IF (idim_loc == 0) THEN
        !print *, my_id, "zero size in gatherv_radar_values(integer)"
      ELSE
        zivector_loc = ivector_loc
      END IF
    END IF

    ! MPI-Routine
    IF (ireceiver < 0) THEN
      IF (ASSOCIATED(ivector_all)) DEALLOCATE(ivector_all)
      ALLOCATE(ivector_all(ilen_all))
      CALL MPI_ALLGATHERV   (zivector_loc, idim_loc, MPI_INTEGER,           &
           ivector_all, idims_all, ipos_all, MPI_INTEGER,           &
           icomm, implcode)
      IF (implcode /= 0) THEN
        ierror  = 2
        yerrmsg = 'gatherv_radar_integers: Error in MPI_ALLGATHERV'
      ENDIF
    ELSEIF (ireceiver < npes) THEN
      ! allocate output vector on receiver, otherwise do not change
      ! allocation state if already allocated
      deall_ivec = .FALSE.
      IF (my_id == ireceiver) THEN
        IF (ASSOCIATED(ivector_all)) DEALLOCATE(ivector_all)
        ALLOCATE(ivector_all(ilen_all))
      ELSE
        ! really necessary???
        IF (.NOT. ASSOCIATED(ivector_all)) THEN
          ALLOCATE(ivector_all(1:1))
          deall_ivec = .TRUE.
        END IF
      END IF

      IF (lnonblck) THEN

        ! Create send tags for the non-blocking communications, unique for each
        !  send PE and message:
        DO ipe=0, npes-1
          itag(ipe+1) = ipe + message_tag*npes
        END DO

        IF (my_id == ireceiver) THEN

          ! If I'm the receiver, simply copy my own data to the receive-vector, if any:
          IF (idims_all(my_id+1) > 0) THEN
            ivector_all(ipos_all(my_id+1)+1:ipos_all(my_id+1)+idims_all(my_id+1)) = &
                 ivector_loc(:)
          END IF

          ! If I'm the receiver, open receive channels from all other PEs and wait for completion:
          DO ipe=0, npes-1
            IF (idims_all(ipe+1) > 0 .AND. ipe /= ireceiver) THEN
              CALL MPI_IRECV(ivector_all(ipos_all(ipe+1)+1:ipos_all(ipe+1)+idims_all(ipe+1)), idims_all(ipe+1), &
                   MPI_INTEGER, &
                   ipe, itag(ipe+1), icomm, irrequest(ipe+1), implcode)
            END IF
          END DO
          DO ipe=0, npes-1
            IF (idims_all(ipe+1) > 0 .AND. ipe /= ireceiver) THEN
              CALL MPI_WAIT(irrequest(ipe+1), istatus, implcode)
            END IF
          END DO

        END IF

        ! If I'm not the receiver, send my data to the receiver:
        IF (my_id /= ireceiver .AND. idims_all(my_id+1) > 0) THEN
          CALL MPI_ISEND(ivector_loc, idim_loc, MPI_INTEGER, &
               ireceiver, itag(my_id+1), icomm, isrequest, implcode)
        ELSE
          isrequest = MPI_REQUEST_NULL
        END IF

      ELSE

        CALL MPI_GATHERV   (zivector_loc, idim_loc, MPI_INTEGER,           &
             ivector_all, idims_all, ipos_all, MPI_INTEGER,           &
             ireceiver, icomm, implcode)
        IF (implcode /= 0) THEN
          ierror  = 2
          yerrmsg = 'gatherv_radar_integers: Error in MPI_GATHERV'
        ENDIF

      END IF

      IF (deall_ivec) THEN
        DEALLOCATE(ivector_all)
      END IF
    ELSE
      ierror  = 3
      yerrmsg = 'gatherv_radar_integers: no valid receiver: ireceiver >= number of PEs'
    ENDIF

    IF (ALLOCATED(zivector_loc)) DEALLOCATE(zivector_loc)

    ipos_all(:) = ipos_all(:) + 1

  ELSE

    IF (my_id == ireceiver .AND. .NOT.ASSOCIATED(ivector_all)) THEN
      ALLOCATE(ivector_all(1:1))
      ivector_all = 0
    END IF

  END IF

END SUBROUTINE gatherv_radar_integers

!==============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine collect_values
!------------------------------------------------------------------------------
!
! SUBROUTINE global_values (vector, idim, ytypop, icomm,      &
!                           ireceiver, yerrmsg, ierror)
!
!------------------------------------------------------------------------------
!
! Description:
!  This routine gathers the values of vector from all nodes and determines
!  the global values defined by the operation ytypop (MIN, MAX or SUM).
!  The values are either given to all nodes in icomm, if ireceiver is not in
!  icomm, or only to ireceiver, otherwise.
!  IF a vector of global maxima is collected, it is possible to add the
!  indices ijmax; then the global indices of the maxima are determined.
!
! Method:
!  MPI_REDUCE, MPI_ALLREDUCE
!
!------------------------------------------------------------------------------

! Following are the different subroutines

!==============================================================================

!+ Subroutine for vector of integers

SUBROUTINE global_vectorint (ivector, idim, ytypop, icomm,    &
                             ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

INTEGER, INTENT(INOUT)    ::                              &
  ivector (idim)        ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
INTEGER                   ::                              &
  ivector_out(idim),          & ! recv-buffer for MPI-routine
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_vectorint
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  ivector_out(:) = 0

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (ivector, ivector_out, idim, MPI_INTEGER, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (ivector, ivector_out, idim, MPI_INTEGER, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  ivector (:) = ivector_out(:)

END SUBROUTINE global_vectorint

!==============================================================================

!==============================================================================

!+ Subroutine for vector of KIND=sp

SUBROUTINE global_vectorreal (rvector, idim, ytypop, icomm,    &
                              ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=sp),           INTENT(INOUT)    ::                              &
  rvector (idim)        ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=sp)                             ::                              &
  rvector_out(idim)             ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_vectorreal
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  rvector_out(:) = 0.0_sp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rvector, rvector_out, idim, MPI_REAL, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rvector, rvector_out, idim, MPI_REAL, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  rvector (:) = rvector_out(:)

END SUBROUTINE global_vectorreal

!==============================================================================

!==============================================================================

!+ Subroutine for vector of KIND=dp

SUBROUTINE global_vectordouble (rvector, idim, ytypop, icomm,    &
                              ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=dp),           INTENT(INOUT)    ::                              &
  rvector (idim)        ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=dp)                             ::                              &
  rvector_out(idim)             ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_vectorreal
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  rvector_out(:) = 0.0_dp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rvector, rvector_out, idim, MPI_DOUBLE_PRECISION, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rvector, rvector_out, idim, MPI_DOUBLE_PRECISION, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  rvector (:) = rvector_out(:)

END SUBROUTINE global_vectordouble

!==============================================================================

!==============================================================================

!+ Subroutine for vector of logicals

SUBROUTINE global_vectorlogical (lvector, idim, ytypop, icomm,    &
                                 ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=*),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

LOGICAL,       INTENT(INOUT)    ::                                         &
  lvector (idim)        ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
LOGICAL                         ::                              &
  lvector_out(idim)             ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_vectorlogical
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  lvector_out(:) = .FALSE.

  SELECT CASE (TRIM(ADJUSTL(ytypop)))
  CASE ('OR','or','Or')
    nzoper = MPI_LOR
  CASE ('AND','and','And')
    nzoper = MPI_LAND
  CASE ('XOR','xor','Xor')
    nzoper = MPI_LXOR
  CASE default
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  END SELECT

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (lvector, lvector_out, idim, MPI_LOGICAL, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (lvector, lvector_out, idim, MPI_LOGICAL, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  lvector (:) = lvector_out(:)

END SUBROUTINE global_vectorlogical

!==============================================================================

!==============================================================================

!+ Subroutine for 2D array of KIND=sp

SUBROUTINE global_array2dsingle (rarray2d, idim, ytypop, icomm,    &
                                 ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  idim(2),            & ! dimension of rarray2d [ni,nj]
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=sp),           INTENT(INOUT)    ::                              &
  rarray2d (idim(1),idim(2))   ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*), INTENT(inout)           ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=sp)                             ::                              &
  rarray2d_out(idim(1),idim(2)) ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_array2dsingle
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  rarray2d_out(:,:) = 0.0_sp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rarray2d, rarray2d_out, PRODUCT(idim), MPI_DOUBLE_PRECISION, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rarray2d, rarray2d_out, PRODUCT(idim), MPI_DOUBLE_PRECISION, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  rarray2d (:,:) = rarray2d_out(:,:)

END SUBROUTINE global_array2dsingle

!==============================================================================

!==============================================================================

!+ Subroutine for 2D array of KIND=dp

SUBROUTINE global_array2ddouble (rarray2d, idim, ytypop, icomm,    &
                                 ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  idim(2),            & ! dimension of rarray2d [ni,nj]
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=dp),           INTENT(INOUT)    ::                              &
  rarray2d (idim(1),idim(2))   ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*), INTENT(inout)           ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=dp)                             ::                              &
  rarray2d_out(idim(1),idim(2)) ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_array2ddouble
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  rarray2d_out(:,:) = 0.0_dp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rarray2d, rarray2d_out, PRODUCT(idim), MPI_DOUBLE_PRECISION, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rarray2d, rarray2d_out, PRODUCT(idim), MPI_DOUBLE_PRECISION, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  rarray2d (:,:) = rarray2d_out(:,:)

END SUBROUTINE global_array2ddouble

!==============================================================================

!==============================================================================

!+ Subroutine for 3D array of KIND=sp

SUBROUTINE global_array3dsingle (rarray3d, idim, ytypop, icomm,    &
                                 ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  idim(3),            & ! dimension of rarray3d [ni,nj,nk]
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=sp),           INTENT(INOUT)    ::                              &
  rarray3d (idim(1),idim(2),idim(3))   ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*), INTENT(inout)           ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=sp)                             ::                              &
  rarray3d_out(idim(1),idim(2),idim(3)) ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_array3dsingle
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  rarray3d_out(:,:,:) = 0.0_sp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rarray3d, rarray3d_out, PRODUCT(idim), MPI_DOUBLE_PRECISION, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rarray3d, rarray3d_out, PRODUCT(idim), MPI_DOUBLE_PRECISION, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  rarray3d (:,:,:) = rarray3d_out(:,:,:)

END SUBROUTINE global_array3dsingle

!==============================================================================

!==============================================================================

!+ Subroutine for 3D array of KIND=dp

SUBROUTINE global_array3ddouble (rarray3d, idim, ytypop, icomm,    &
                                 ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  idim(3),            & ! dimension of rarray3d [ni,nj,nk]
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=dp),           INTENT(INOUT)    ::                              &
  rarray3d (idim(1),idim(2),idim(3))   ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*), INTENT(inout)           ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=dp)                             ::                              &
  rarray3d_out(idim(1),idim(2),idim(3)) ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_array3ddouble
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  rarray3d_out(:,:,:) = 0.0_dp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rarray3d, rarray3d_out, PRODUCT(idim), MPI_DOUBLE_PRECISION, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rarray3d, rarray3d_out, PRODUCT(idim), MPI_DOUBLE_PRECISION, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  rarray3d (:,:,:) = rarray3d_out(:,:,:)

END SUBROUTINE global_array3ddouble

!==============================================================================

!==============================================================================

!+ Subroutine for one integer

SUBROUTINE global_int (ivector, ytypop, icomm,    &
                       ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

INTEGER, INTENT(INOUT)    ::                              &
  ivector               ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
INTEGER                   ::                              &
  ivector_out,                & ! recv-buffer for MPI-routine
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_int
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  ivector_out = 0

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (ivector, ivector_out,    1, MPI_INTEGER, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (ivector, ivector_out,    1, MPI_INTEGER, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  ivector = ivector_out

END SUBROUTINE global_int

!==============================================================================

!+ Subroutine for one logical

SUBROUTINE global_logical (lvector, ytypop, icomm,    &
                           ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=*),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

LOGICAL, INTENT(INOUT)    ::                              &
  lvector               ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
LOGICAL                   ::                              &
  lvector_out                   ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_logical
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  lvector_out = .FALSE.

  SELECT CASE (TRIM(ADJUSTL(ytypop)))
  CASE ('OR','or','Or')
    nzoper = MPI_LOR
  CASE ('AND','and','And')
    nzoper = MPI_LAND
  CASE ('XOR','xor','Xor')
    nzoper = MPI_LXOR
  CASE default
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  END SELECT

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (lvector, lvector_out,    1, MPI_LOGICAL, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (lvector, lvector_out,    1, MPI_LOGICAL, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  lvector = lvector_out

END SUBROUTINE global_logical

!==============================================================================

!==============================================================================

!+ Subroutine for one KIND=dp

SUBROUTINE global_double (rvector, ytypop, icomm,    &
                        ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=dp),           INTENT(INOUT)    ::                          &
  rvector               ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=dp)                             ::                          &
  rvector_out                   ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_double
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  rvector_out = 0.0_dp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rvector, rvector_out,    1, MPI_DOUBLE_PRECISION, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rvector, rvector_out,    1, MPI_DOUBLE_PRECISION, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  rvector = rvector_out

END SUBROUTINE global_double

!==============================================================================

!+ Subroutine for one KIND=sp

SUBROUTINE global_real (rvector, ytypop, icomm,    &
                        ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  icomm,              & ! MPI-communicator
  ireceiver             ! ID of the node that gets the results

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=sp),           INTENT(INOUT)    ::                          &
  rvector               ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=sp)                             ::                          &
  rvector_out                   ! recv-buffer for MPI-routine

INTEGER                   ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_real
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  rvector_out = 0.0_sp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values_radar'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values_radar'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rvector, rvector_out,    1, MPI_REAL, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rvector, rvector_out,    1, MPI_REAL, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCE failed in global_values_radar'
    RETURN
  ENDIF

  rvector = rvector_out

END SUBROUTINE global_real

!==============================================================================
!==============================================================================

  !==============================================================================
  !+ Module procedure in data_radar for the generation of the derived MPI type of
  !  radar_meta_type
  !------------------------------------------------------------------------------

SUBROUTINE def_mpi_radar_meta_type (nobstimes_max_loc, mpi_radar_meta_typ)

  !------------------------------------------------------------------------------
  !
  ! Description:    Defines the derived MPI type for the radar_meta_type
  !                 Currently the array pointers need to be taken into
  !                 account in the definition of the MPI type (definition
  !                 of amount of memory), but they
  !                 have to be communicated seperately (I have not found a way
  !                 to include them in the communication of the derived type).
  !
  ! Method     :    See IBM MPI book (online and locatable with Google)
  !------------------------------------------------------------------------------
  !
  ! Subroutine / Function arguments
  ! Scalar arguments with intent(in):
  INTEGER, INTENT(IN)  :: nobstimes_max_loc
  ! Scalar arguments with intent(out):
  INTEGER, INTENT(OUT) :: mpi_radar_meta_typ
  !------------------------------------------------------------------------------
  ! Local scalars:
  INTEGER, SAVE     ::  mpi_radar_meta_typ_
  INTEGER                                              ::&
       ierr
  INTEGER (kind=MPI_ADDRESS_KIND)                      ::&
       offsets(5), lb, extent
  INTEGER                                              ::&
       mpitypes(5), blocklengths(5)
  LOGICAL, SAVE :: firstcall = .TRUE.

  !- End of header
  !==============================================================================

  !------------------------------------------------------------------------------
  !- Begin SUBROUTINE
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Section 1: Initializations
  !------------------------------------------------------------------------------

!  IF (firstcall) THEN

    ! Setup description of the 23 imp_reals variables in the structure
    ! and of the number of elevations
    ! and the number of smoothing points and weights
    ! and the number of observation times
    ! and the 3 elements in dt_obs:
    offsets(1)      = 0
    mpitypes(1)     = MPI_DOUBLE_PRECISION
    blocklengths(1) = 23 + (4+nscanstrategies_max+nobstimes_max_loc)*nel_max + 4*ngpsm_max + 4*nobstimes_max_loc + 1*3

    ! Get the size of MPI_DOUBLE_PRECISION in bytes
    CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, extent, ierr)
!!! NOTE: this only is correct if kind type wpfw0 = dp

    ! Setup description of the 21 imp_integers in the structure
    ! and the 2 vectors with dimension ndatakind
    ! and the obs_startrec,obs_endrec with dimension ndatakind*nobstimes_max_loc
    ! and the elev. indices for feedback files
    ! and the elevation index list for the composites comp_dbzsim and comp_dbzobs:
    offsets(2)      = offsets(1) + blocklengths(1)*extent
    mpitypes(2)     = MPI_INTEGER
    blocklengths(2) = 21 + nscanstrategies_max + ndatakind + 3*ndatakind*nobstimes_max_loc + 3*nel_max + nel_composite_max

    ! Get the size of MPI_INTEGER in bytes
    CALL MPI_TYPE_GET_EXTENT(MPI_INTEGER, lb, extent, ierr)

    !    ! Setup description of the 4 + ndatakind*nobstimes_max_loc logicals in the structure
    offsets(3)      = offsets(2) + blocklengths(2)*extent
    mpitypes(3)     = MPI_LOGICAL
    blocklengths(3) = 4 + ndatakind*nobstimes_max_loc

    ! Get the size of MPI_LOGICAL in bytes
    CALL MPI_TYPE_GET_EXTENT(MPI_LOGICAL, lb, extent, ierr)

    ! Setup description of the 2 strings (len = 3 and 14) in the structure
    ! and the number of observation times
    offsets(4)      = offsets(3) + blocklengths(3)*extent
    mpitypes(4)     = MPI_CHARACTER
    blocklengths(4) =   1*3                     & ! 3    = LEN(radar_meta%station_name)
         + 10                                   & ! 10   = LEN(radar_meta%scanname)
         + nobstimes_max_loc*14                 & ! 14   = LEN(radar_meta%obs_cdate)
         + ndatakind*nobstimes_max_loc*cobsflen & ! cobsflen  = LEN(radar_meta%obsfile)
         + nobstimes_max_loc*15                 & ! 15   = LEN(radar_meta%obsfile_format)
         + 1*60                                 & ! 60   = LEN(radar_meta%fdbkfile)
         + 8*cvarlen                            & ! cvarlen = LEN(radar_meta%obs_hdf5_varname_XXX)
         + cobsflen                               ! dummy buffer to compensate for memory alignment in mpi transfers

    ! Get the size of MPI_CHARACTER in bytes and set the upper boundary
    ! of the structure ( = memory block) explicitly by MPI_UB to prevent any
    ! byte-number-rounding in transfers
    CALL MPI_TYPE_GET_EXTENT(MPI_CHARACTER, lb, extent, ierr)

    ! Set the upper boundary of the memory block:
    offsets(5)      = offsets(4) + blocklengths(4)*extent
    mpitypes(5)     = MPI_UB
    blocklengths(5) = 1

    ! Define structured type and commit it
    CALL MPI_TYPE_CREATE_STRUCT(5, blocklengths, offsets, mpitypes, mpi_radar_meta_typ_, ierr)
    CALL MPI_TYPE_COMMIT(mpi_radar_meta_typ_,ierr)

!!$ PRINT *, 'ULI size mpi_radar_meta_typ: ', SUM(blocklengths * (/8, 4, 1, 1, 0/))

    firstcall = .FALSE.

!  END IF

  mpi_radar_meta_typ = mpi_radar_meta_typ_

END SUBROUTINE def_mpi_radar_meta_type

!==============================================================================

SUBROUTINE def_mpi_polmp_type (mpi_polmp_typ)

  !------------------------------------------------------------------------------
  !
  ! Description:    Defines the derived MPI type for the polMP type.
  !                 This type is defined in radar_dbzcalc_params_type.f90 and used
  !                 for shape and orientation info for Tmatrix scattering
  !                 calculations in EMVORADO.
  !
  ! Method     :    See IBM MPI book (online and locatable with Google)
  !------------------------------------------------------------------------------


  ! Subroutine / Function arguments
  !------------------------------------------------------------------------------
  ! Scalar arguments with intent(out):
  INTEGER, INTENT(OUT)                                 ::&
       mpi_polmp_typ
  !------------------------------------------------------------------------------
  ! Local scalars:
  INTEGER, SAVE    ::  mpi_polmp_typ_
  INTEGER                                              ::&
       ierr
  INTEGER (kind=MPI_ADDRESS_KIND)                      ::&
       offsets(3), lb, extent
  INTEGER                                              ::&
       mpitypes(3), blocklengths(3)
  LOGICAL, SAVE :: firstcall = .TRUE.
  !- End of header
  !==============================================================================

  !------------------------------------------------------------------------------
  !- Begin SUBROUTINE
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Section 1: Initializations
  !------------------------------------------------------------------------------

  IF (firstcall) THEN

    ! Setup description of the 8 MPI_CHARACTER word in the structure:
    offsets(1)      = 0
    mpitypes(1)     = MPI_CHARACTER
    blocklengths(1) = 8

    ! Get the size of MPI_CHARACTER in bytes
    CALL MPI_TYPE_GET_EXTENT(MPI_CHARACTER, lb, extent, ierr)

    ! Setup description of the 4 MPI_DOUBLE_PRECISION variables in the structure:
    offsets(2)      = offsets(1) + blocklengths(1)*extent
    mpitypes(2)     = MPI_DOUBLE_PRECISION
    blocklengths(2) = 4

    ! Get the size of MPI_DOUBLE_PRECISION in bytes and set the upper boundary
    ! of the structure ( = memory block) explicitly by MPI_UB to prevent any
    ! byte-number-rounding in transfers
    CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, extent, ierr)
!!! NOTE: this only is correct if kind type wpfw0 = dp

    ! Set the upper boundary of the memory block:
    offsets(3)      = offsets(2) + blocklengths(2)*extent
    mpitypes(3)     = MPI_UB
    blocklengths(3) = 1


    ! Define structured type and commit it
    CALL MPI_TYPE_CREATE_STRUCT(3, blocklengths, offsets, mpitypes, &
                                mpi_polmp_typ_, ierr)
    CALL MPI_TYPE_COMMIT(mpi_polmp_typ_,ierr)
    firstcall = .FALSE.

  END IF

  mpi_polmp_typ = mpi_polmp_typ_

END SUBROUTINE def_mpi_polmp_type

!==============================================================================

SUBROUTINE def_mpi_dbzcalc_params_type (mpi_polmp_typ, mpi_dbzcalc_params_typ)

  !------------------------------------------------------------------------------
  !
  ! Description:    Defines the derived MPI type for the dbzcalc_params type.
  !                 This type is defined in radar_data.f90 and used for
  !                 the radar reflectivity calculations in EMVORADO.
  !
  ! Method     :    See IBM MPI book (online and locatable with Google)
  !------------------------------------------------------------------------------


  ! Subroutine / Function arguments
  !------------------------------------------------------------------------------
  INTEGER, INTENT(IN)                                  ::&
       mpi_polmp_typ
  ! Scalar arguments with intent(out):
  INTEGER, INTENT(OUT)                                 ::&
       mpi_dbzcalc_params_typ
  !------------------------------------------------------------------------------
  ! Local scalars:
  INTEGER, SAVE    ::  mpi_dbzcalc_params_typ_
  INTEGER                                              ::&
       ierr
  INTEGER (kind=MPI_ADDRESS_KIND)                      ::&
       offsets(6), lb, extent
  INTEGER                                              ::&
       mpitypes(6), blocklengths(6)
  LOGICAL, SAVE :: firstcall = .TRUE.
  !- End of header
  !==============================================================================

  !------------------------------------------------------------------------------
  !- Begin SUBROUTINE
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Section 1: Initializations
  !------------------------------------------------------------------------------

  IF (firstcall) THEN

    ! Setup description of the melt.scheme params (6) per icy hymet class (4) + 3
    ! MPI_DOUBLE_PRECISION variables in the structure:
    offsets(1)      = 0
    mpitypes(1)     = MPI_DOUBLE_PRECISION
    blocklengths(1) = 4*6 + 3

    ! Get the size of MPI_DOUBLE_PRECISION in bytes
    CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, extent, ierr)
!!! NOTE: this only is correct if kind type wpfw0 = dp

    ! Setup description of the 5 MPI_INTEGER variables + 1 dummy for memory alignment in the structure:
    offsets(2)      = offsets(1) + blocklengths(1)*extent
    mpitypes(2)     = MPI_INTEGER
    blocklengths(2) = 5 + 1

    ! Get the size of MPI_INTEGER in bytes
    CALL MPI_TYPE_GET_EXTENT(MPI_INTEGER, lb, extent, ierr)


    ! Setup description of the 1 logical vector(6 elements) and two flags in the structure:
    offsets(3)      = offsets(2) + blocklengths(2)*extent
    mpitypes(3)     = MPI_LOGICAL
    blocklengths(3) = 1*6 + 2

    ! Get the size of MPI_LOGICAL in bytes
    CALL MPI_TYPE_GET_EXTENT(MPI_LOGICAL, lb, extent, ierr)

    ! Setup description of the 16 words of 12 MPI_CHARACTER each in the structure:
    offsets(4)      = offsets(3) + blocklengths(3)*extent
    mpitypes(4)     = MPI_CHARACTER
    blocklengths(4) = 16 * 12

    ! Get the size of MPI_CHARACTER in bytes
    CALL MPI_TYPE_GET_EXTENT(MPI_CHARACTER, lb, extent, ierr)

    ! Setup description of the 5 mpi_polmp_typ in the structure:
    offsets(5)      = offsets(4) + blocklengths(4)*extent
    mpitypes(5)     = mpi_polmp_typ
    blocklengths(5) = 5

    ! Get the size of mpi_polmp_typ in bytes and set the upper boundary
    ! of the structure ( = memory block) explicitly by MPI_UB to prevent any
    ! byte-number-rounding in transfers
    CALL MPI_TYPE_GET_EXTENT(mpi_polmp_typ, lb, extent, ierr)

    ! Set the upper boundary of the memory block:
    offsets(6)      = offsets(5) + blocklengths(5)*extent
    mpitypes(6)     = MPI_UB
    blocklengths(6) = 1


    ! Define structured type and commit it
    CALL MPI_TYPE_CREATE_STRUCT(6, blocklengths, offsets, mpitypes, &
                                mpi_dbzcalc_params_typ_, ierr)
    CALL MPI_TYPE_COMMIT(mpi_dbzcalc_params_typ_,ierr)
    firstcall = .FALSE.

  END IF

  mpi_dbzcalc_params_typ = mpi_dbzcalc_params_typ_

END SUBROUTINE def_mpi_dbzcalc_params_type

!==============================================================================

SUBROUTINE def_mpi_compmeta_type (mpi_compmeta_typ)

  !------------------------------------------------------------------------------
  !
  ! Description:    Defines the derived MPI type for the composite_meta_type.
  !                 This type is defined in radar_data.f90 and used for
  !                 the composite grid in EMVORADO.
  !
  ! Method     :    See IBM MPI book (online and locatable with Google)
  !------------------------------------------------------------------------------


  ! Subroutine / Function arguments
  !------------------------------------------------------------------------------
  ! Scalar arguments with intent(out):
  INTEGER, INTENT(OUT)        ::  mpi_compmeta_typ
  !------------------------------------------------------------------------------
  ! Local scalars:
  INTEGER                          ::&
       ierr
  INTEGER (kind=MPI_ADDRESS_KIND)  ::&
       offsets(5), lb, extent
  INTEGER                          ::&
       mpitypes(5), blocklengths(5)

  !- End of header
  !==============================================================================

  !------------------------------------------------------------------------------
  !- Begin SUBROUTINE
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Section 1: Initializations
  !------------------------------------------------------------------------------

  ! Setup description of the 8 MPI_DOUBLE_PRECISION variables in the structure:
  offsets(1)      = 0
  mpitypes(1)     = MPI_DOUBLE_PRECISION
  blocklengths(1) = 8

  ! Get the size of MPI_DOUBLE_PRECISION in bytes
  CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, extent, ierr)
!!! NOTE: this only is correct if kind type wpfw0 = dp


  ! Setup description of the 2 MPI_INTEGER variables in the structure:
  offsets(2)      = blocklengths(1)*extent
  mpitypes(2)     = MPI_INTEGER
  blocklengths(2) = 2

  ! Get the size of MPI_INTEGER in bytes
  CALL MPI_TYPE_GET_EXTENT(MPI_INTEGER, lb, extent, ierr)


  ! Setup description of the 0 MPI_CHARACTER each in the structure:
  offsets(3)      = offsets(2) + blocklengths(2)*extent
  mpitypes(3)     = MPI_CHARACTER
  blocklengths(3) = 0

  ! Get the size of MPI_CHARACTER in bytes
  CALL MPI_TYPE_GET_EXTENT(MPI_CHARACTER, lb, extent, ierr)


  ! Setup description of the 0 logical flag(s) in the structure:
  offsets(4)      = offsets(3) + blocklengths(3)*extent
  mpitypes(4)     = MPI_LOGICAL
  blocklengths(4) = 0

  ! Get the size of MPI_LOGICAL in bytes and set the upper boundary
  ! of the structure ( = memory block) explicitly by MPI_UB to prevent any
  ! byte-number-rounding in transfers
  CALL MPI_TYPE_GET_EXTENT(MPI_LOGICAL, lb, extent, ierr)

  ! Set the upper boundary of the memory block:
  offsets(5)      = offsets(4) + blocklengths(4)*extent
  mpitypes(5)     = MPI_UB
  blocklengths(5) = 1


  ! Define structured type and commit it
  CALL MPI_TYPE_CREATE_STRUCT(5, blocklengths, offsets, mpitypes, &
                              mpi_compmeta_typ, ierr)
  CALL MPI_TYPE_COMMIT(mpi_compmeta_typ, ierr)

END SUBROUTINE def_mpi_compmeta_type

!==============================================================================

SUBROUTINE def_mpi_voldataostream_type (mpi_voldataostream_typ)

  !------------------------------------------------------------------------------
  !
  ! Description:    Defines the derived MPI type for the type t_voldata_ostream.
  !                 This type is defined in radar_data_namelist.f90 and used for
  !                 the namelist parameters of the voldata output streams in EMVORADO.
  !
  ! Method     :    See IBM MPI book (online and locatable with Google)
  !------------------------------------------------------------------------------


  ! Subroutine / Function arguments
  !------------------------------------------------------------------------------
  ! Scalar arguments with intent(out):
  INTEGER, INTENT(OUT)        ::  mpi_voldataostream_typ
  !------------------------------------------------------------------------------
  ! Local scalars:
  INTEGER                          ::&
       ierr
  INTEGER (kind=MPI_ADDRESS_KIND)  ::&
       offsets(5), lb, extent
  INTEGER                          ::&
       mpitypes(5), blocklengths(5)

  !- End of header
  !==============================================================================

  !------------------------------------------------------------------------------
  !- Begin SUBROUTINE
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Section 1: Initializations
  !------------------------------------------------------------------------------

  ! Setup description of the 2 MPI_DOUBLE_PRECISION variables in the structure:
  offsets(1)      = 0
  mpitypes(1)     = MPI_DOUBLE_PRECISION
  blocklengths(1) = 2

  ! Get the size of MPI_DOUBLE_PRECISION in bytes
  CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, extent, ierr)
!!! NOTE: this only is correct if kind type wpfw0 = dp


  ! Setup description of the 0 MPI_INTEGER variables in the structure:
  offsets(2)      = blocklengths(1)*extent
  mpitypes(2)     = MPI_INTEGER
  blocklengths(2) = 0

  ! Get the size of MPI_INTEGER in bytes
  CALL MPI_TYPE_GET_EXTENT(MPI_INTEGER, lb, extent, ierr)


  ! Setup description of the 4*12 + noutput_fields_max*12 + 2*cmaxlen MPI_CHARACTER in the structure:
  offsets(3)      = offsets(2) + blocklengths(2)*extent
  mpitypes(3)     = MPI_CHARACTER
  blocklengths(3) = 4*12 + noutput_fields_max*12 + 2*cmaxlen

  ! Get the size of MPI_CHARACTER in bytes
  CALL MPI_TYPE_GET_EXTENT(MPI_CHARACTER, lb, extent, ierr)


  ! Setup description of the 0 logical flag(s) in the structure:
  offsets(4)      = offsets(3) + blocklengths(3)*extent
  mpitypes(4)     = MPI_LOGICAL
  blocklengths(4) = 0

  ! Get the size of MPI_LOGICAL in bytes and set the upper boundary
  ! of the structure ( = memory block) explicitly by MPI_UB to prevent any
  ! byte-number-rounding in transfers
  CALL MPI_TYPE_GET_EXTENT(MPI_LOGICAL, lb, extent, ierr)

  ! Set the upper boundary of the memory block:
  offsets(5)      = offsets(4) + blocklengths(4)*extent
  mpitypes(5)     = MPI_UB
  blocklengths(5) = 1


  ! Define structured type and commit it
  CALL MPI_TYPE_CREATE_STRUCT(5, blocklengths, offsets, mpitypes, &
                              mpi_voldataostream_typ, ierr)
  CALL MPI_TYPE_COMMIT(mpi_voldataostream_typ, ierr)

END SUBROUTINE def_mpi_voldataostream_type

!===============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine gather_values_radar
!------------------------------------------------------------------------------
!
! SUBROUTINE gather_values_radar (vector_in, vector_out, idim, npes ,   &
!                           ireceiver, icomm, yerrmsg, ierror)
!
!------------------------------------------------------------------------------
!
! Description:
!  This routine gathers the values of vector_in from all nodes to
!  the PE with ID ireceiver. If ireceiver < 0, vector_out is send to all PEs.
!
! Method:
!  MPI_GATHER or MPI_ALLGATHER
!
!------------------------------------------------------------------------------

! Following are the different subroutines

!==============================================================================

!+ Subroutine for one integer

SUBROUTINE gather_one_int  (ivector_in, ivector_out, idim, npes, &
                            ireceiver, icomm, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT(IN)       ::                              &
  idim,                 & ! dimension of ivector
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver               ! PE that receives the vector

INTEGER, INTENT(IN)       ::                              &
  ivector_in              ! subdomain field

INTEGER, INTENT(OUT)      ::                              &
  ivector_out(npes)       ! subdomain field

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER, INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                   ::                              &
  implcode                ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gather_one_int
!------------------------------------------------------------------------------

  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  ! MPI-Routine
  IF (ireceiver < 0) THEN
    CALL MPI_ALLGATHER (ivector_in,  idim, MPI_INTEGER,           &
                        ivector_out, idim, MPI_INTEGER,           &
                        icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 1
      yerrmsg = 'Error in MPI_ALLGATHER'
    ENDIF
  ELSEIF (ireceiver < npes) THEN
    CALL MPI_GATHER    (ivector_in,  idim, MPI_INTEGER,           &
                        ivector_out, idim, MPI_INTEGER,           &
                        ireceiver, icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 2
      yerrmsg = 'Error in MPI_GATHER'
    ENDIF
  ELSE
    ierror  = 3
    yerrmsg = 'no valid receiver: ireceiver >= number of PEs'
  ENDIF

END SUBROUTINE gather_one_int

!==============================================================================
!==============================================================================

END MODULE radar_parallel_utilities

!==============================================================================
