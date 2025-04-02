!+ Utility module for working with the parallel program (without I/O)
!------------------------------------------------------------------------------
#if __NEC_VERSION__ >= 30005    /* options directive supported since 3.0.5 */
!NEC$ options "-fno-outerloop-unroll"
#endif

MODULE  parallel_utilities

!------------------------------------------------------------------------------
!
! Description:
!   This module provides routines that should enable programmers, not
!   experienced in parallelization, to work with the LM. It contains
!   routines that handle tasks as the distribution or collection of
!   values.
!   This module is shared between the COSMO model and DACE program packages !
!
!   Routines currently contained are
!
!     - remark:
!     - distribute_values:
!     - gather_grib:
!     - gather_values:
!     - gatherv_radar_values
!     - scatter_grib:
!     - scatter_values:
!     - global_values:
!     - combine_subarrays:
!
!     - grid_indices:
!     - print_point:
!     - print_column:
!     - distribute_field:
!     - gather_field:
!
!     - print_lm_field:    (for GME2LM)
!     - dump_lm_field:
!     - lm_field_stats:    (for GME2LM)
!     - scatterv_values:   (for nudging/LETKF)
!     - gatherv_values:    (for nudging/LETKF)
!     - gather_all:        (for nudging)
!     - reports_all2all:   (for nudging/LETKF)
!     - exchange_profiles: (for latent heat nudging)
!
! Method:
!   Calls to the message-passing library MPI
!
!
! Current Code Owners:
!  For COSMO:                              For DACE:
!  DWD, Ulrich Schaettler                  DWD, Christoph Schraff
!  phone:  +49  69  8062 2739              phone: +49 69 8062 2725
!  fax:    +49  69  8062 3721              fax:   +49 69 8062 3721
!  email:  ulrich.schaettler@dwd.de        email: christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.8        1998/08/03 Ulrich Schaettler
!  Eliminated dependencies on grib parameters.
! 1.10       1998/09/29 Ulrich Schaettler
!  New Subroutines included; routine printfield eliminated.
! 1.12       1998/10/19 Ulrich Schaettler
!  New Subroutine for computing dot-product included.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adaptations to use this module also in GME2LM
! 1.32       1999/08/24 Guenther Doms
!  Declaration of non-used local variables deleted.
! 1.39       2000/05/03 Ulrich Schaettler
!  Check allocation of IPU_iposition that it is only done once.
! 2.8        2001/07/06 Ulrich Schaettler
!  Introduced subroutine combine_subarrays (was in src_output.f90 before)
! 2.11       2001/09/28 Ulrich Schaettler
!  Corrected a bug in the collective communications subroutines
!  (result was send to all, if nreceiver = 0 in global_values)
! 2.17       2002/05/08 Ulrich Schaettler
!  New SR scatter_grib, gather_grib for I/O communications in irealgrib-format;
!  changed combine_subarrays to a generic interface
! 3.6        2003/12/11 Ulrich Schaettler
!  Corrected a bug in Subroutine ij_local
! 3.13       2004/12/03 Ulrich Schaettler
!  Added "ireceiver" in SR gather_field, so that also other PEs can gather
!  a field (and not only PE 0).
!  Added a new SR exchange_profiles (for latent heat nudging)
! 3.18       2006/03/03 Ulrich Schaettler
!  Included gather_grib in the generic formulation of gather_values
!  Included scatter_grib in the generic formulation of scatter_values
!  Implemented generic formulation of combine_subarrays
!  Introduced idouble/isingle as KIND parameters instead of ireals/irealgrib
!  in the generic formulation of some routines
! V3_23        2007/03/30 Simone Campagna, Davide Cesari
!  In SR global_values, some vectors have to be initialized
! V4_1         2007/12/04 Ulrich Schaettler, Michael Baldauf
!  More generic implementations for subroutine distribute_values and
!  for subroutine global_values (+global_int8)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler
!  Introduced new SR dump_lm_field to dump a global field in full precision
! V4_23        2012/05/10 Burkhardt Rockel, Ulrich Schaettler
!  Modified SR distribute_field, so that another sender than PE 0 can be specified
! V4_25        2012/09/28 Carlos Osuna
!  Add new SR distribute_vartab for CHARACTER(LEN=1) TYPE for the variable table
! V4_28        2013/07/12 Christoph Schraff
!  Added new SRs 'reports_all2all' (individual info sent/received by each node),
!  'scatterv_values' (one node sends indvidual info to other nodes), and
!  'gatherv_values' (one or all nodes receive individual info from other nodes).
!  Introduced 'if(n)def NOMPI' to make this shared module compatible with 3DVAR.
! V5_1         2014-11-28 Ulrich Blahak
!  Introduced new subroutine global_vectorlogical() for global_values
!  Introduced new interface gatherv_radar_values with routines for integers and reals
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Robin Faulwetter
!  Introduced additional implementations for generic routine gatherv_values
! V5_3         2015-10-09 Hans-Juergen Panitz, Ulrich Schaettler
!  Replace KIND=8 for INTEGER variables by generic definition i8
!     from kind_parameters (HJP)
! V5_4d        2016-12-12 Michael Bender
! New subroutines to support the 3dvar MPI interface:
!     'p_abort', 'gatherv_args', 'alltoallv_args'
! V5_4f        2017-09-01 Ulrich Blahak, Christoph Schraff, Andreas Rhodin
!  Added single precision versions for SR combine_subarrays(),
!  Implemented generic version (incl. single precision) for gather_field(),
!   distribute_field(). (UB)
!  Update from DACE V1_51: Bug fixes (if NOMPI, or if my_cart_id not defined). (CS,AR)
! V5_5b        2018-10-29 Xavier Lapillonne
!  OpenACC port to GPU:
!    - SR init_par_utilities needs lcompute_pe to update device on compute PE
!    - SR combine_2D_[double|single]: added lacc
!    - FCTs i_global, j_global, ij_local: added acc-statements and
!                    modifications to avoid RETURN statement
!    - SR gather_field_[double|single]: added lacc and had to duplicate
!        MPI-call because acc-directives used do not support if, which has
!        to be written as a Fortran IF (lzacc) now.
!    - SR exchange_profiles: ported with OpenACC
! upcoming based on V5_6  Christoph Schraff
!  - Set subroutines and functions used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - (From DACE:) Set preference of MPI module over mpif.h on Cray.
!  - For DACE only: use of 'mo_mpi' replaced by 'mo_mpi_dace'.
!  - ifdef's for preference of MPI module over mpif.h adjusted to work both
!    for COSMO and in DACE.
!  - Bug fix: MPI_ALLGATHERV argument ireceiver removed in gatherv_radar_values.
!
! ==============================================================================
! DACE history:
! ------------
!  1_51       2017-02-24 Andreas Rhodin : update shared modules from COSMO 5.04d
! ==============================================================================
! CAUTION: This module is used by both the DACE and COSMO model codes.       !!!
!!!        Therefore, anybody wanting to remove or modify existing routines  !!!
!!!        in this module in the context of either of these programs must    !!!
!!!        consult the 'current code owner' of this module for the other     !!!
!!!        program, in order to allow for checking that the modification     !!!
!!!        will comply with both program packages. This must be done before  !!!
!!!        the modification is put into the Version Control System (VCS).    !!!
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
!------------------------------------------------------------------------------

USE kind_parameters , ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    sp,        & ! KIND-type parameter for real variables (single precision)
    dp,        & ! KIND-type parameter for real variables (double precision)
    i4,        & ! KIND-type parameter for integer*4 variables
    i8           ! KIND-type parameter for integer*8 variables

#ifdef NOMPI
USE environment, ONLY : model_abort   ! program abort in case of error
#endif

#ifndef NOMPI
! prefer MPI module over mpif.h if either of following compiler directives set:
! COSMO runs on CRAY  or  MPI module is defined in DACE
#if defined __CRAYXC  ||  defined HAVE_MPI_MOD
USE mpi
#endif
#endif

!==============================================================================

IMPLICIT NONE

PRIVATE
          ! routines used in COSMO and DACE
PUBLIC :: init_par_utilities,  distribute_values,  gather_values,             &
          global_values,  gather_field,  scatterv_values,  gatherv_values,    &
          gather_all,  reports_all2all,  exchange_profiles,                   &
          gatherv_radar_values

          ! routines used in COSMO only
PUBLIC :: remark,  distribute_vartab,  scatter_values,  combine_subarrays,    &
          i_global, j_global, ij_local,  distribute_field

          ! routines used in src_obs_gnss.f90 in COSMO;
          ! note that in DACE, these routines exist in 'mo_mpi_dace.f90'!
PUBLIC :: p_abort,  gatherv_args,  alltoallv_args

! unused routines: print_point, print_column, dump_lm_field

!==============================================================================

! include statements
#ifndef NOMPI
! use mpif.h only if corresponding compiler directives are not set, such that:
#ifndef __CRAYXC        /* COSMO is not compiled to run on CRAY  */
#ifndef HAVE_MPI_MOD    /* and MPI module is not defined in DACE */
INCLUDE "mpif.h"
#endif
#endif
#endif

!==============================================================================

! Private Declarations
  INTEGER                 , PRIVATE     ::      &
    ! local dimensions of a subdomain
    IPU_idim, IPU_jdim, IPU_kdim,               &

    ! dimensions of total domain
    IPU_idim_tot, IPU_jdim_tot, IPU_kdim_tot,   &

    ! maxima of ie and je, resp.
    IPU_idim_max, IPU_jdim_max,                 &

    ! for computations in the local domain without overlapping boundaries
    IPU_istart_comp, IPU_iend_comp,             &
    IPU_jstart_comp, IPU_jend_comp,             &

    ! total number of processors: numpex * numpey + numpeio
    IPU_numpe,                                  &

    ! total number of processors for computing: numpex * numpey
    IPU_numpeco,                                &

    ! number of processors in both directions of the cartesian grid
    IPU_numpex, IPU_numpey,                     &

    ! number of extra processors for doing asynchronous IO
    IPU_numpeio,                                &

    ! number of boundary lines
    IPU_nbounds,                                &

    ! communicator for the cartesian topology
    IPU_comm_cart,                              &

    ! rank of this subdomain in icommcart
    IPU_cart_id,                                &

    ! datatypes for message passing
    IPU_real, IPU_integer, IPU_single, IPU_double

  INTEGER                 , ALLOCATABLE, PRIVATE     ::      &
    ! positions of subdomains in total domain
    IPU_ipositions(:,:)

! Declare directive required for variables that are used in acc routine.
!$acc declare create(IPU_ipositions,IPU_integer,IPU_real) &
!$acc create(IPU_cart_id,IPU_comm_cart,IPU_nbounds)       &
!$acc create(IPU_numpex,IPU_numpey,IPU_numpeco,IPU_numpe) &
!$acc create(IPU_idim_tot,IPU_jdim_tot,IPU_kdim_tot)      &
!$acc create(IPU_idim,IPU_jdim,IPU_kdim)

!==============================================================================

! Interface Blocks
INTERFACE global_values
  MODULE PROCEDURE                        &
    global_vectorint,                     &
    global_vectorreal,                    &
    global_vectorlogical,                 &
    global_int,                           &
    global_int8,                          &
    global_real,                          &
    global_vectorreal_indices
END INTERFACE

INTERFACE distribute_values
  MODULE PROCEDURE                        &
    distribute_kind8,                     &
    distribute_kind4,                     &
    distribute_oneinteger,                &
    distribute_idouble,                   &
    distribute_isingle,                   &
    distribute_onedouble,                 &
    distribute_onesingle,                 &
    distribute_logical,                   &
    distribute_onelogical,                &
    distribute_character,                 &
    distribute_onecharacter
END INTERFACE

INTERFACE gather_values
  MODULE PROCEDURE                        &
    gather_integers,                      &
    gather_one_int,                       &
    gather_reals,                         &
    gather_2D_double,                     &
    gather_2D_single
END INTERFACE

INTERFACE scatter_values
  MODULE PROCEDURE                        &
    scatter_integers,                     &
    scatter_double,                       &
    scatter_single
END INTERFACE

INTERFACE combine_subarrays
  MODULE PROCEDURE                        &
    combine_2D_double,                    &
    combine_2D_single,                    &
    combine_double,                       &
    combine_single
END INTERFACE

INTERFACE scatterv_values
  MODULE PROCEDURE                        &
    scatterv_integers,                    &
    scatterv_reals
END INTERFACE

INTERFACE gatherv_values
  MODULE PROCEDURE                        &
    gatherv_integers,                     &
    gatherv_integers_2D,                  &
    gatherv_reals,                        &
    gatherv_reals_2D,                     &
    gatherv_reals_3D,                     &
    gatherv_logical_2D
END INTERFACE

INTERFACE gatherv_radar_values
  MODULE PROCEDURE                        &
    gatherv_radar_integers,               &
    gatherv_radar_reals
END INTERFACE

INTERFACE gather_field
  MODULE PROCEDURE                        &
    gather_field_double,                  &
    gather_field_single
END INTERFACE

INTERFACE distribute_field
  MODULE PROCEDURE                        &
    distribute_field_double,              &
    distribute_field_single
END INTERFACE

!==============================================================================

CONTAINS

!==============================================================================
!+ Initializes private variables for module parallel_utilities
!------------------------------------------------------------------------------

SUBROUTINE init_par_utilities                                                &
     (idim, jdim, kdim, idim_tot, jdim_tot, kdim_tot, idim_max, jdim_max,    &
      istart_comp, iend_comp, jstart_comp, jend_comp,                        &
      numpe, numpex, numpey, numpeio, ipositions, nbounds, icartcomm,        &
      icart_id, imp_real, imp_int, lcompute_pe)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER                 , INTENT(IN)  ::  &
    idim, jdim, kdim,             & ! local dimensions of a subdomain
    idim_tot, jdim_tot, kdim_tot, & ! dimensions of total domain
    idim_max, jdim_max,           & ! maxima of ie and je, resp.
    istart_comp, iend_comp,       & ! for computations in the local domains
    jstart_comp, jend_comp,       & !    without overlapping boundaries
    numpe,                        & ! total number of processors:
                                    !    numpex * numpey + numpeio
    numpex, numpey,               & ! number of processors in both directions
                                    !    of the cartesian grid
    numpeio,                      & ! number of extra processors for doing
                                    !    asynchronous IO
    nbounds,                      & ! number of boundary lines
    ipositions(0:numpex*numpey-1,4),& ! positions of subdomains in total domain
    icartcomm,                    & ! communicator for the cartesian topology
    icart_id,                     & ! rank of this subdomain in icommcart
    imp_real, imp_int               ! datatypes for message passing

  LOGICAL, INTENT(IN) ::             lcompute_pe  ! compute PE flag

! Local Variables
  INTEGER                               ::  &
    izstat

!==============================================================================

  ! dimensions of the domains
  IPU_idim        = idim
  IPU_jdim        = jdim
  IPU_kdim        = kdim
  IPU_idim_tot    = idim_tot
  IPU_jdim_tot    = jdim_tot
  IPU_kdim_tot    = kdim_tot
  IPU_idim_max    = idim_max
  IPU_jdim_max    = jdim_max

  ! for computations in the local domains without overlapping boundaries
  IPU_istart_comp = istart_comp
  IPU_iend_comp   = iend_comp
  IPU_jstart_comp = jstart_comp
  IPU_jend_comp   = jend_comp

  ! number of processors
  IPU_numpe       = numpe
  IPU_numpex      = numpex
  IPU_numpey      = numpey
  IPU_numpeco     = numpex * numpey
  IPU_numpeio     = numpeio

  ! parallel program and related variables
  IPU_nbounds     = nbounds
  IPU_comm_cart   = icartcomm
  IPU_cart_id     = icart_id

  ! datatypes for message passing
  IPU_real        = imp_real
#ifndef NOMPI
  IPU_single      = MPI_REAL4
  IPU_double      = MPI_REAL8
#else
  IPU_single      = -1
  IPU_double      = -1
#endif
  IPU_integer     = imp_int

  ! Allocate IPU_ipositions at the first call
  IF (.NOT. ALLOCATED(IPU_ipositions)) THEN
    ALLOCATE (IPU_ipositions(0:IPU_numpex*IPU_numpey-1, 4), STAT=izstat)
  ENDIF

  ! Set IPU_ipositions
  IPU_ipositions(:,:)  = ipositions(:,:)

  ! Update variable to compute PE
  ! Only the compute PE may have attached accelerator
  IF (lcompute_pe) THEN
    ! Update all "declare" module variables
    !$acc update device(IPU_ipositions,IPU_integer,IPU_real)
    !$acc update device(IPU_cart_id,IPU_comm_cart,IPU_nbounds)
    !$acc update device(IPU_numpex,IPU_numpey,IPU_numpeco,IPU_numpe)
    !$acc update device(IPU_idim_tot,IPU_jdim_tot,IPU_kdim_tot)
    !$acc update device(IPU_idim,IPU_jdim,IPU_kdim)
  END IF

END SUBROUTINE init_par_utilities

!==============================================================================

!+ Prints a remark from PE 0
!------------------------------------------------------------------------------

SUBROUTINE remark (my_id, yroutine, youtstr)

!------------------------------------------------------------------------------
!
! Description:
!   This routine writes a string and the procedure from where it is
!   called to standard output for e.g. debugging purposes.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameter list:
INTEGER                 , INTENT(IN)  ::  my_id
CHARACTER, INTENT(IN)  :: youtstr*(*), yroutine*(*)

!------------------------------------------------------------------------------

! Begin subroutine remark

IF (my_id == 0) THEN
  write(*,'(A,A,A,A)') 'Remark from Routine: ',yroutine,':  ',youtstr
ENDIF

END SUBROUTINE remark

!==============================================================================

!+ Defines all subroutines for the generic routine distribute_values
!------------------------------------------------------------------------------
!
! SUBROUTINE distribute_values (buffer, ibufferlen, isender, idatatype,
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

!+ Subroutine for array of kind=8 integers

SUBROUTINE distribute_kind8 (buffer, ibufferlen, isender, idatatype,     &
                             icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
INTEGER (KIND=i8),        INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                         &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                                    ::                             &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, ibufferlen, idatatype, isender,                 &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_kind8')
#endif

END SUBROUTINE distribute_kind8

!==============================================================================

!==============================================================================

!+ Subroutine for array of kind=4 integers

SUBROUTINE distribute_kind4 (buffer, ibufferlen, isender, idatatype,     &
                             icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
INTEGER (KIND=i4),        INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                         &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                                    ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, ibufferlen, idatatype, isender,                 &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_kind4')
#endif

END SUBROUTINE distribute_kind4

!==============================================================================

!==============================================================================

!+ Subroutine for one model integers

SUBROUTINE distribute_oneinteger(buffer, ibufferlen, isender, idatatype,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
INTEGER                 , INTENT(INOUT)    ::                              &
  buffer                ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                                    ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, 1, idatatype, isender, icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_oneinteger')
#endif

END SUBROUTINE distribute_oneinteger

!==============================================================================

!==============================================================================

!+ Subroutine for array of doubles

SUBROUTINE distribute_idouble   (buffer, ibufferlen, isender, idatatype,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
REAL    (KIND=dp),         INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                                    ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, ibufferlen, idatatype, isender,   &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_idouble')
#endif

END SUBROUTINE distribute_idouble

!==============================================================================

!==============================================================================

!+ Subroutine for one double

SUBROUTINE distribute_onedouble (buffer, ibufferlen, isender, idatatype,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
REAL    (KIND=dp),         INTENT(INOUT)    ::                              &
  buffer                ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                                    ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, ibufferlen, idatatype, isender,   &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_onedouble')
#endif

END SUBROUTINE distribute_onedouble

!==============================================================================

!==============================================================================

!+ Subroutine for array of singles

SUBROUTINE distribute_isingle   (buffer, ibufferlen, isender, idatatype,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
REAL    (KIND=sp),         INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                                    ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, ibufferlen, idatatype, isender,   &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_isingle')
#endif

END SUBROUTINE distribute_isingle

!==============================================================================

!==============================================================================

!+ Subroutine for one single

SUBROUTINE distribute_onesingle (buffer, ibufferlen, isender, idatatype,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
REAL    (KIND=sp),         INTENT(INOUT)    ::                              &
  buffer                ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                                    ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, ibufferlen, idatatype, isender,   &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_onesingle')
#endif

END SUBROUTINE distribute_onesingle

!==============================================================================

!==============================================================================

!+ Subroutine for array of default logicals

SUBROUTINE distribute_logical   (buffer, ibufferlen, isender, idatatype,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
LOGICAL,                  INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                                    ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, ibufferlen, idatatype, isender,   &
                  icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_logical')
#endif

END SUBROUTINE distribute_logical

!==============================================================================

!==============================================================================

!+ Subroutine for one default logical

SUBROUTINE distribute_onelogical(buffer, ibufferlen, isender, idatatype,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
LOGICAL,                  INTENT(INOUT)    ::                              &
  buffer                ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER                                    ::                              &
  implcode              ! local error code
!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, 1, idatatype, isender, icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_onelogical')
#endif

END SUBROUTINE distribute_onelogical

!==============================================================================

!==============================================================================

!+ Subroutine for array of characters

SUBROUTINE distribute_character (buffer, ibufferlen, isender, idatatype,     &
                                 icommunicator, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
CHARACTER (LEN=100),      INTENT(INOUT)    ::                              &
  buffer (ibufferlen)   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars and Arrays
INTEGER                                    ::                              &
  my_comm_id, implcode, i, j

INTEGER    ::     &  ! Standard integer
  intbuf(100)

!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
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
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_character')
#endif

END SUBROUTINE distribute_character

!==============================================================================

!==============================================================================

!+ Subroutine for one word of characters

SUBROUTINE distribute_onecharacter (buffer, ibufferlen, isender, idatatype,  &
                                    icommunicator,  ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
CHARACTER (LEN=*),        INTENT(INOUT)    ::                              &
  buffer                ! character to be broadcasted

CHARACTER (LEN=100)  :: internal_buffer   ! buffer to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars and Arrays
INTEGER                                    ::                              &
  my_comm_id, implcode, i, j

INTEGER    ::     &  ! Standard integer
  intbuf(100)

!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
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
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_onecharacter')
#endif

END SUBROUTINE distribute_onecharacter

!==============================================================================
!==============================================================================

!+ Subroutine for distributing array of LEN=1 characters

SUBROUTINE distribute_vartab (buffer, ibufferlen, isender, idatatype,  &
                              icommunicator,  ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  ibufferlen,         & ! length of the buffer
  isender,            & ! sending processor
  idatatype,          & ! type of buffer
  icommunicator         ! involved group of processors

! Array arguments with intent(inout):
CHARACTER (LEN=1),        INTENT(INOUT)    ::                              &
  buffer(:)              ! character to be broadcasted

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars and Arrays
INTEGER                                    ::     implcode

!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
  CALL MPI_BCAST (buffer, ibufferlen, MPI_CHARACTER, isender, icommunicator, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = implcode
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_vartab')
#endif

END SUBROUTINE distribute_vartab

!==============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine gather_values
!------------------------------------------------------------------------------
!
! SUBROUTINE gather_values (vector_in, vector_out, idim, npes , idatatype,   &
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

!+ Subroutine for vector of integers

SUBROUTINE gather_integers (ivector_in, ivector_out, idim, npes , idatatype, &
                            ireceiver, icomm, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,                 & ! dimension of ivector
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver,            & ! PE that receives the vector
  idatatype               ! data type of vector

INTEGER                 , INTENT(IN)       ::                              &
  ivector_in (idim)       ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ivector_out(idim,npes)  ! subdomain field

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode                ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gather_integers
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  ! MPI-Routine
  IF (ireceiver < 0) THEN
    CALL MPI_ALLGATHER (ivector_in,  idim, idatatype,           &
                        ivector_out, idim, idatatype,           &
                        icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 1
      yerrmsg = 'Error in MPI_ALLGATHER'
    ENDIF
  ELSEIF (ireceiver < npes) THEN
    CALL MPI_GATHER    (ivector_in,  idim, idatatype,           &
                        ivector_out, idim, idatatype,           &
                        ireceiver, icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 2
      yerrmsg = 'Error in MPI_GATHER'
    ENDIF
  ELSE
    ierror  = 3
    yerrmsg = 'no valid receiver: ireceiver >= number of PEs'
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gather_integers')
#endif

END SUBROUTINE gather_integers

!==============================================================================
!==============================================================================

!+ Subroutine for vector of reals

SUBROUTINE gather_reals    (rvector_in, rvector_out, idim, npes,  idatatype, &
                            ireceiver, icomm, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,                 & ! dimension of ivector
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver,            & ! PE that receives the vector
  idatatype               ! data type of vector

REAL (KIND=wp),           INTENT(IN)       ::                              &
  rvector_in (idim)       ! subdomain field

REAL (KIND=wp),           INTENT(OUT)      ::                              &
  rvector_out(idim,npes)  ! subdomain field

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode                ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gather_reals
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  ! MPI-Routine
  IF (ireceiver < 0) THEN
    CALL MPI_ALLGATHER (rvector_in,  idim, idatatype,           &
                        rvector_out, idim, idatatype,           &
                        icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 1
      yerrmsg = 'Error in MPI_ALLGATHER'
    ENDIF
  ELSEIF (ireceiver < npes) THEN
    CALL MPI_GATHER    (rvector_in,  idim, idatatype,           &
                        rvector_out, idim, idatatype,           &
                        ireceiver, icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 2
      yerrmsg = 'Error in MPI_GATHER'
    ENDIF
  ELSE
    ierror  = 3
    yerrmsg = 'no valid receiver: ireceiver >= number of PEs'
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gather_reals')
#endif

END SUBROUTINE gather_reals

!==============================================================================
!==============================================================================

!+ Subroutine for one integer

SUBROUTINE gather_one_int  (ivector_in, ivector_out, idim, npes, idatatype, &
                            ireceiver, icomm, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,                 & ! dimension of ivector
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver,            & ! PE that receives the vector
  idatatype               ! data type of vector

INTEGER                 , INTENT(IN)       ::                              &
  ivector_in              ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ivector_out(npes)       ! subdomain field

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode                ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gather_one_int
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  ! MPI-Routine
  IF (ireceiver < 0) THEN
    CALL MPI_ALLGATHER (ivector_in,  idim, idatatype,           &
                        ivector_out, idim, idatatype,           &
                        icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 1
      yerrmsg = 'Error in MPI_ALLGATHER'
    ENDIF
  ELSEIF (ireceiver < npes) THEN
    CALL MPI_GATHER    (ivector_in,  idim, idatatype,           &
                        ivector_out, idim, idatatype,           &
                        ireceiver, icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 2
      yerrmsg = 'Error in MPI_GATHER'
    ENDIF
  ELSE
    ierror  = 3
    yerrmsg = 'no valid receiver: ireceiver >= number of PEs'
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gather_one_int')
#endif

END SUBROUTINE gather_one_int

!==============================================================================
!==============================================================================

!+ Subroutine for 2D array of double precision reals

SUBROUTINE gather_2D_double(rvector_in, rvector_out, idim1, idim2, npes,     &
                            idatatype, ireceiver, icomm, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim1, idim2,         & ! dimensions of ivector
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver,            & ! PE that receives the vector
  idatatype               ! data type of vector

REAL (KIND=dp),           INTENT(IN)       ::                              &
  rvector_in (idim1,idim2)! subdomain field

REAL (KIND=dp),           INTENT(OUT)      ::                              &
  rvector_out(idim1,idim2,npes)  ! subdomain field

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode                ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gather_reals
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  ! MPI-Routine
  IF (ireceiver < 0) THEN
    CALL MPI_ALLGATHER (rvector_in,  idim1*idim2, idatatype,           &
                        rvector_out, idim1*idim2, idatatype,           &
                        icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 1
      yerrmsg = 'Error in MPI_ALLGATHER'
    ENDIF
  ELSEIF (ireceiver < npes) THEN
    CALL MPI_GATHER    (rvector_in,  idim1*idim2, idatatype,           &
                        rvector_out, idim1*idim2, idatatype,           &
                        ireceiver, icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 2
      yerrmsg = 'Error in MPI_GATHER'
    ENDIF
  ELSE
    ierror  = 3
    yerrmsg = 'no valid receiver: ireceiver >= number of PEs'
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gather_2D_double')
#endif

END SUBROUTINE gather_2D_double

!==============================================================================
!==============================================================================

SUBROUTINE gather_2D_single(rvector_in, rvector_out, idim1, idim2, npes,  &
                            idatatype, ireceiver, icomm, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim1, idim2,         & ! dimensions of ivector
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver,            & ! PE that receives the vector
  idatatype               ! data type of vector

REAL (KIND=sp),           INTENT(IN)       ::                              &
  rvector_in (idim1,idim2)! subdomain field

REAL (KIND=sp),           INTENT(OUT)      ::                              &
  rvector_out(idim1,idim2,npes)  ! subdomain field

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode                ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gather_reals
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  ! MPI-Routine
  IF (ireceiver < 0) THEN
    CALL MPI_ALLGATHER (rvector_in,  idim1*idim2, idatatype,           &
                        rvector_out, idim1*idim2, idatatype,           &
                        icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 1
      yerrmsg = 'Error in MPI_ALLGATHER'
    ENDIF
  ELSEIF (ireceiver < npes) THEN
    CALL MPI_GATHER    (rvector_in,  idim1*idim2, idatatype,           &
                        rvector_out, idim1*idim2, idatatype,           &
                        ireceiver, icomm, implcode)
    IF (implcode /= 0) THEN
      ierror  = 2
      yerrmsg = 'Error in MPI_GATHER'
    ENDIF
  ELSE
    ierror  = 3
    yerrmsg = 'no valid receiver: ireceiver >= number of PEs'
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gather_2D_single')
#endif

END SUBROUTINE gather_2D_single

!==============================================================================
!==============================================================================
!
!+ Defines all subroutines for the generic routine gatherv_radar_values
!------------------------------------------------------------------------------
!
! SUBROUTINE gatherv_radar_values (vector_loc, vector_all, ipos_all, npes , idatatype,   &
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
!+ Subroutine for vector of reals
SUBROUTINE gatherv_radar_reals (rvector_loc, rvector_all, ipos_all, ilen_all, &
                      npes,  idatatype, ireceiver, icomm, yerrmsg, ierror)

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
INTEGER                 , INTENT(IN)       ::                              &
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver,            & ! PE that receives the vector
  idatatype               ! data type of vector

REAL (KIND=wp),           INTENT(IN)       ::                              &
  rvector_loc(:)        ! local input vector

! UB>> the SX-compiler does not compile POINTER, INTENT(...) together!
!REAL (KIND=wp),     POINTER, INTENT(out)    ::                              &
REAL (KIND=wp),     POINTER    ::                              &
  rvector_all(:)        ! global total vector; MUST BE INITIALIZED BEFORE INPUT! e.g., nullify(rvector_all)

INTEGER                 , INTENT(out)    ::                              &
  ipos_all(npes)      , & ! positions in output vector
  ilen_all                ! length of global vector

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode            , & ! error status variable for MPI routines
  idim_loc            , & ! dimension of local vector
  ipe                     ! loop variable for pe's
INTEGER                                    ::                              &
  idims_all(npes)         ! array of dimensions of all input vectors
REAL (KIND=wp),     ALLOCATABLE            ::                              &
  zrvector_loc(:)
LOGICAL :: deall_rvec
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_reals
!------------------------------------------------------------------------------

  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  ! determine dimension of input vector
  idim_loc = SIZE(rvector_loc)

  ! gather dimensions of all input vectors in idims_all
  CALL gather_values ( idim_loc, idims_all, 1, npes, IPU_integer, -1,    &
       icomm, yerrmsg, implcode )
  IF (implcode /= 0) THEN
    ierror  = 1
    yerrmsg = "gatherv_radar_reals: error in gathering the dimensions of the local vectors: " // yerrmsg
    RETURN
  ENDIF

  ! determine position vector
  ipos_all(1) = 0
  DO ipe = 2 , npes
     ipos_all(ipe)  = ipos_all(ipe-1) + idims_all(ipe-1)
  END DO
  ilen_all  =  ipos_all(npes) + idims_all(npes)

  IF (ilen_all /= 0) THEN

    ALLOCATE(zrvector_loc(MAX(idim_loc,1)))
    zrvector_loc = 0.0_wp

    IF (idim_loc == 0) THEN
      !print *, IPU_cart_id, "zero size in gatherv_radar_values(real)"
    ELSE
      zrvector_loc = rvector_loc
    END IF

    ! MPI-Routine
    IF (ireceiver < 0) THEN
      ! allocate output vector
      IF (ASSOCIATED(rvector_all)) DEALLOCATE(rvector_all)
      ALLOCATE(rvector_all(ilen_all))
#ifndef NOMPI
      CALL MPI_ALLGATHERV   (zrvector_loc, idim_loc, idatatype,           &
           rvector_all, idims_all, ipos_all, idatatype,           &
           icomm, implcode)
#else
      CALL model_abort (-1, -1, 'MPI not available', 'gatherv_radar_reals')
#endif
      IF (implcode /= 0) THEN
        ierror  = 2
        yerrmsg = 'gatherv_radar_reals: Error in MPI_ALLGATHERV'
      ENDIF
    ELSEIF (ireceiver < npes) THEN
      ! allocate output vector on receiver, otherwise do not change
      ! allocation state if already allocated
      deall_rvec = .FALSE.
      IF (IPU_cart_id == ireceiver) THEN
        IF (ASSOCIATED(rvector_all)) DEALLOCATE(rvector_all)
        ALLOCATE(rvector_all(ilen_all))
      ELSE
        ! really necessary???
        IF (.NOT. ASSOCIATED(rvector_all)) THEN
          ALLOCATE(rvector_all(1:1))
          deall_rvec = .TRUE.
        END IF
      END IF
#ifndef NOMPI
      CALL MPI_GATHERV   (zrvector_loc, idim_loc, idatatype,           &
           rvector_all, idims_all, ipos_all, idatatype,           &
           ireceiver, icomm, implcode)
#else
      CALL model_abort (-1, -1, 'MPI not available', 'gatherv_radar_reals')
#endif
      IF (implcode /= 0) THEN
        ierror  = 2
        yerrmsg = 'gatherv_radar_reals: Error in MPI_GATHERV'
      ENDIF
      IF (deall_rvec) THEN
        DEALLOCATE(rvector_all)
      END IF
    ELSE
      ierror  = 3
      yerrmsg = 'gatherv_radar_reals: no valid receiver: ireceiver >= number of PEs'
    ENDIF

    DEALLOCATE(zrvector_loc)

    ipos_all(:) = ipos_all(:) + 1

  ELSE

    IF (IPU_cart_id == ireceiver .AND. .NOT.ASSOCIATED(rvector_all)) THEN
      ALLOCATE(rvector_all(1:1))
      rvector_all = 0.0_wp
    END IF

  END IF

END SUBROUTINE gatherv_radar_reals

!==============================================================================

!+ Subroutine for vector of integers
SUBROUTINE gatherv_radar_integers (ivector_loc, ivector_all, ipos_all,      &
              ilen_all, npes,  idatatype, ireceiver, icomm, yerrmsg, ierror)

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
INTEGER                 , INTENT(IN)       ::                              &
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  ireceiver,            & ! PE that receives the vector
  idatatype               ! data type of vector

INTEGER                 ,       INTENT(IN)       ::                              &
  ivector_loc(:)        ! local input vector

! UB>> the SX-compiler does not compile POINTER, INTENT(...) together!
!INTEGER                 , pointer, intent(out)    ::                              &
INTEGER                 , POINTER    ::                              &
  ivector_all(:)        ! global total vector; MUST BE INITIALIZED BEFORE INPUT! e.g., nullify(ivector_all)

INTEGER                 , INTENT(out)    ::                              &
  ipos_all(npes)      , & ! positions in output vector
  ilen_all                ! length of global vector

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode            , & ! error status variable for MPI routines
  idim_loc            , & ! dimension of local vector
  ipe                      ! loop variable for pe's
INTEGER                                    ::                              &
  idims_all(npes)         ! array of dimensions of all input vectors
INTEGER                 , ALLOCATABLE            ::                              &
  zivector_loc(:)
LOGICAL :: deall_ivec
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_reals
!------------------------------------------------------------------------------

  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  ! determine dimension of input vector
  idim_loc = SIZE(ivector_loc)

  ! gather dimensions of all input vectors in idims_all
  CALL gather_values ( idim_loc, idims_all, 1, npes, IPU_integer, -1,    &
       icomm, yerrmsg, implcode )
  IF (implcode /= 0) THEN
    ierror  = 1
    yerrmsg = "gatherv_radar_integers: error in gathering the dimensions of the local vectors: " // yerrmsg
    RETURN
  ENDIF

  ! determine position vector
  ipos_all(1) = 0
  DO ipe = 2 , npes
    ipos_all(ipe)  = ipos_all(ipe-1) + idims_all(ipe-1)
  END DO
  ilen_all  =  ipos_all(npes) + idims_all(npes)

  IF (ilen_all /= 0) THEN

    ALLOCATE(zivector_loc(MAX(idim_loc,1)))
    zivector_loc = 0.0_wp

    IF (idim_loc == 0) THEN
      !print *, IPU_cart_id, "zero size in gatherv_radar_values(integer)"
    ELSE
      zivector_loc = ivector_loc
    END IF

    ! MPI-Routine
    IF (ireceiver < 0) THEN
      IF (ASSOCIATED(ivector_all)) DEALLOCATE(ivector_all)
      ALLOCATE(ivector_all(ilen_all))
#ifndef NOMPI
      CALL MPI_ALLGATHERV   (zivector_loc, idim_loc, idatatype,           &
           ivector_all, idims_all, ipos_all, idatatype,           &
           icomm, implcode)
#else
      CALL model_abort (-1, -1, 'MPI not available', 'gatherv_radar_integers')
#endif
      IF (implcode /= 0) THEN
        ierror  = 2
        yerrmsg = 'gatherv_radar_integers: Error in MPI_ALLGATHERV'
      ENDIF
    ELSEIF (ireceiver < npes) THEN
      ! allocate output vector on receiver, otherwise do not change
      ! allocation state if already allocated
      deall_ivec = .FALSE.
      IF (IPU_cart_id == ireceiver) THEN
        IF (ASSOCIATED(ivector_all)) DEALLOCATE(ivector_all)
        ALLOCATE(ivector_all(ilen_all))
      ELSE
        ! really necessary???
        IF (.NOT. ASSOCIATED(ivector_all)) THEN
          ALLOCATE(ivector_all(1:1))
          deall_ivec = .TRUE.
        END IF
      END IF
#ifndef NOMPI
      CALL MPI_GATHERV   (zivector_loc, idim_loc, idatatype,           &
           ivector_all, idims_all, ipos_all, idatatype,           &
           ireceiver, icomm, implcode)
#else
      CALL model_abort (-1, -1, 'MPI not available', 'gatherv_radar_integers')
#endif
      IF (implcode /= 0) THEN
        ierror  = 2
        yerrmsg = 'gatherv_radar_integers: Error in MPI_GATHERV'
      ENDIF
      IF (deall_ivec) THEN
        DEALLOCATE(ivector_all)
      END IF
    ELSE
      ierror  = 3
      yerrmsg = 'gatherv_radar_integers: no valid receiver: ireceiver >= number of PEs'
    ENDIF

    DEALLOCATE(zivector_loc)

    ipos_all(:) = ipos_all(:) + 1

  ELSE

    IF (IPU_cart_id == ireceiver .AND. .NOT.ASSOCIATED(ivector_all)) THEN
      ALLOCATE(ivector_all(1:1))
      ivector_all = 0
    END IF

  END IF

END SUBROUTINE gatherv_radar_integers

!==============================================================================
!==============================================================================

!+ Defines all subroutines for the generic routine scatter_values
!------------------------------------------------------------------------------
!
! SUBROUTINE scatter_values (vector_in, vector_out, idim, npes, idatatype,   &
!                           isender, icomm, yerrmsg, ierror)
!
!------------------------------------------------------------------------------
!
! Description:
!  This routine scatters the values of vector_in from the PE with ID isender
!  to all nodes.
!
! Method:
!  MPI_SCATTER
!
!==============================================================================

! Following are the different subroutines

!==============================================================================

!+ Subroutine for vector of integers

SUBROUTINE scatter_integers (ivector_in, ivector_out, idim, npes, idatatype, &
                            isender, icomm, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,                 & ! dimension of ivector
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  isender,              & ! PE that receives the vector
  idatatype               ! data type of vector

INTEGER                 , INTENT(IN)       ::                              &
  ivector_in (idim,npes)  ! global field

INTEGER                 , INTENT(OUT)      ::                              &
  ivector_out(idim)       ! subdomain field

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode                ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine scatter_integers
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  IF (npes > 1) THEN
    IF ( (0 <= isender) .AND. (isender < npes) ) THEN
      CALL MPI_SCATTER (ivector_in,  idim, idatatype,           &
                        ivector_out, idim, idatatype,           &
                        isender, icomm, implcode)
      IF (implcode /= 0) THEN
        ierror  = 1
        yerrmsg = 'Error in MPI_SCATTER'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid sender'
    ENDIF
  ELSE
    ivector_out(:) = ivector_in(:,1)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'scatter_integers')
#endif

END SUBROUTINE scatter_integers

!==============================================================================
!==============================================================================

!+ Subroutine for vector of reals

SUBROUTINE scatter_double   (rvector_in, rvector_out, idim, npes, idatatype, &
                            isender, icomm, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,                 & ! dimension of ivector
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  isender,              & ! PE that receives the vector
  idatatype               ! data type of vector

REAL (KIND=dp),           INTENT(IN)       ::                              &
  rvector_in (idim,npes)  ! global field

REAL (KIND=dp),           INTENT(OUT)      ::                              &
  rvector_out(idim)       ! subdomain field

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode                ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine scatter_reals
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  IF (npes > 1) THEN
    IF ( (0 <= isender) .AND. (isender < npes) ) THEN
      CALL MPI_SCATTER (rvector_in,  idim, idatatype,           &
                        rvector_out, idim, idatatype,           &
                        isender, icomm, implcode)
      IF (implcode /= 0) THEN
        ierror  = 1
        yerrmsg = 'Error in MPI_SCATTER'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid sender'
    ENDIF
  ELSE
    rvector_out(:) = rvector_in(:,1)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'scatter_double')
#endif

END SUBROUTINE scatter_double

!==============================================================================
!==============================================================================

!+ Subroutine for scattering an input field

SUBROUTINE scatter_single (rvector_in, rvector_out, idim, npes, idatatype, &
                           isender, icomm, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,                 & ! dimension of ivector
  npes,                 & ! number of PEs
  icomm,                & ! MPI-communicator
  isender,              & ! PE that receives the vector
  idatatype               ! data type of vector

REAL (KIND=sp),           INTENT(IN)       ::                              &
  rvector_in (idim,npes)  ! global field

REAL (KIND=sp),           INTENT(OUT)      ::                              &
  rvector_out(idim)       ! subdomain field

CHARACTER (LEN= *)      , INTENT(OUT)      ::                              &
  yerrmsg                 ! string for error messages

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                  ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                              &
  implcode                ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine scatter_grib
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode = 0
  ierror   = 0
  yerrmsg  = '   '

  IF (npes > 1) THEN
    IF ( (0 <= isender) .AND. (isender < npes) ) THEN
      CALL MPI_SCATTER (rvector_in,  idim, idatatype,           &
                        rvector_out, idim, idatatype,           &
                        isender, icomm, implcode)
      IF (implcode /= 0) THEN
        ierror  = 1
        yerrmsg = 'Error in MPI_SCATTER'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid sender'
    ENDIF
  ELSE
    rvector_out(:) = rvector_in(:,1)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'scatter_single')
#endif

END SUBROUTINE scatter_single

!==============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine collect_values
!------------------------------------------------------------------------------
!
! SUBROUTINE global_values (vector, idim, ytypop, idatatype, icomm,      &
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

SUBROUTINE global_vectorint (ivector, idim, ytypop, idatatype, icomm,    &
                             ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver,          & ! ID of the node that gets the results
  idatatype             ! data type of vector

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

INTEGER                 , INTENT(INOUT)    ::                              &
  ivector (idim)        ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

#ifndef NOMPI
! Local Variables
INTEGER                                    ::                              &
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
    yerrmsg = 'no valid operation type in global_values'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (ivector, ivector_out, idim, idatatype, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (ivector, ivector_out, idim, idatatype, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCED failed in global_values'
    RETURN
  ENDIF

  ivector (:) = ivector_out(:)
#else
  CALL model_abort (-1, -1, 'MPI not available', 'global_vectorint')
#endif

END SUBROUTINE global_vectorint

!==============================================================================

!==============================================================================

!+ Subroutine for vector of reals

SUBROUTINE global_vectorreal (rvector, idim, ytypop, idatatype, icomm,    &
                              ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver,          & ! ID of the node that gets the results
  idatatype             ! data type of vector

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=wp),           INTENT(INOUT)    ::                              &
  rvector (idim)        ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

#ifndef NOMPI
! Local Variables
REAL (KIND=wp)                             ::                              &
  rvector_out(idim)             ! recv-buffer for MPI-routine

INTEGER                                    ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_vectorreal
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  rvector_out(:) = 0.0_wp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rvector, rvector_out, idim, idatatype, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rvector, rvector_out, idim, idatatype, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCED failed in global_values'
    RETURN
  ENDIF

  rvector (:) = rvector_out(:)
#else
  CALL model_abort (-1, -1, 'MPI not available', 'global_vectorreal')
#endif

END SUBROUTINE global_vectorreal

!==============================================================================

!==============================================================================

!+ Subroutine for vector of logicals

SUBROUTINE global_vectorlogical (lvector, idim, ytypop, idatatype, icomm,    &
                                 ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver,          & ! ID of the node that gets the results
  idatatype             ! data type of vector

CHARACTER (LEN=*),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

LOGICAL,       INTENT(INOUT)    ::                              &
  lvector (idim)        ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
LOGICAL                         ::                              &
  lvector_out(idim)             ! recv-buffer for MPI-routine

INTEGER                                    ::                              &
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

#ifndef NOMPI
  SELECT CASE (TRIM(ADJUSTL(ytypop)))
  CASE ('OR','or','Or')
    nzoper = MPI_LOR
  CASE ('AND','and','And')
    nzoper = MPI_LAND
  CASE ('XOR','xor','Xor')
    nzoper = MPI_LXOR
  CASE default
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values'
    RETURN
  END SELECT

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (lvector, lvector_out, idim, idatatype, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (lvector, lvector_out, idim, idatatype, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCED failed in global_values'
    RETURN
  ENDIF
#else   /* NOMPI */
  lvector_out(1:idim) = lvector(1:idim)     ! Trivial reduction for nzsize=1
#endif

  lvector (:) = lvector_out(:)

END SUBROUTINE global_vectorlogical

!==============================================================================

!==============================================================================

!+ Subroutine for one integer

SUBROUTINE global_int (ivector, idim, ytypop, idatatype, icomm,    &
                       ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver,          & ! ID of the node that gets the results
  idatatype             ! data type of vector

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

INTEGER                 , INTENT(INOUT)    ::                              &
  ivector               ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
INTEGER                                    ::                              &
  ivector_out,                & ! recv-buffer for MPI-routine
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_int
!------------------------------------------------------------------------------

#ifndef NOMPI
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
    yerrmsg = 'no valid operation type in global_values'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (ivector, ivector_out,    1, idatatype, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (ivector, ivector_out,    1, idatatype, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCED failed in global_values'
    RETURN
  ENDIF

  ivector = ivector_out
#else
  CALL model_abort (-1, -1, 'MPI not available', 'global_int')
#endif

END SUBROUTINE global_int

!==============================================================================

!==============================================================================

!+ Subroutine for one 8 Byte-integer

SUBROUTINE global_int8 (ivector, idim, ytypop, idatatype, icomm,    &
                       ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver,          & ! ID of the node that gets the results
  idatatype             ! data type of vector

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

INTEGER (KIND=i8),        INTENT(INOUT)    ::                              &
  ivector               ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
INTEGER                                    ::                              &
  ivector_out,                & ! recv-buffer for MPI-routine
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_int
!------------------------------------------------------------------------------

#ifndef NOMPI
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
    yerrmsg = 'no valid operation type in global_values'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (ivector, ivector_out,    1, idatatype, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (ivector, ivector_out,    1, idatatype, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCED failed in global_values'
    RETURN
  ENDIF

  ivector = ivector_out
#else
  CALL model_abort (-1, -1, 'MPI not available', 'global_int8')
#endif

END SUBROUTINE global_int8

!==============================================================================

!==============================================================================

!+ Subroutine for one real

SUBROUTINE global_real (rvector, idim, ytypop, idatatype, icomm,    &
                        ireceiver, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver,          & ! ID of the node that gets the results
  idatatype             ! data type of vector

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=wp),           INTENT(INOUT)    ::                              &
  rvector               ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                ! error status variable

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=wp)                             ::                              &
  rvector_out                   ! recv-buffer for MPI-routine

INTEGER                                    ::                              &
  nzsize,                     & ! size of communicator icomm
  nzoper                        ! for type of operation

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_real
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  ierror   = 0
  rvector_out = 0.0_wp

  IF     (ytypop == 'SUM') THEN
    nzoper = MPI_SUM
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values'
    RETURN
  ENDIF

  ! Get size of communicator icomm
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLREDUCE (rvector, rvector_out,    1, idatatype, nzoper,   &
                        icomm, ierror)
  ELSE
    CALL MPI_REDUCE    (rvector, rvector_out,    1, idatatype, nzoper,   &
                        ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_REDUCED failed in global_values'
    RETURN
  ENDIF

  rvector = rvector_out
#else
  CALL model_abort (-1, -1, 'MPI not available', 'global_real')
#endif

END SUBROUTINE global_real

!==============================================================================

!==============================================================================

!+ Subroutine for vector of reals with indices

SUBROUTINE global_vectorreal_indices (rvector, idim, ytypop, idatatype, icomm,&
                                      ireceiver, ijmax, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim,               & ! dimension of ivector
  icomm,              & ! MPI-communicator
  ireceiver,          & ! ID of the node that gets the results
  idatatype             ! data type of vector

CHARACTER (LEN=3),        INTENT(IN)       ::                              &
  ytypop                ! type of operation: MIN, MAX or SUM

REAL (KIND=wp),           INTENT(INOUT)    ::                              &
  rvector (idim)        ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                ! error status variable

INTEGER                 , INTENT(INOUT)    ::                              &
  ijmax(2, idim)        ! for local / global indices

CHARACTER (LEN=*)                          ::                              &
  yerrmsg               ! for error message

!------------------------------------------------------------------------------

! Local Variables
REAL (KIND=wp)                             ::                              &
  rvector_out(idim,IPU_numpe),& ! for gathering all MAX-values
  zgmax

INTEGER                                    ::                              &
  indices(2,idim,IPU_numpe),  & ! for gathering all indices
  nzsize,                     & ! size of communicator icomm
  myz_id,                     & ! ID of this PE in communicator icomm
  nzoper, izindex, l, n

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine global_vectorreal_indices
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  ierror   = 0
  rvector_out(:,:) = 0.0_wp

  IF     (ytypop == 'SUM') THEN
    ierror  = 1
    yerrmsg = 'cannot determine indices for summation'
    RETURN
  ELSEIF (ytypop == 'MAX') THEN
    nzoper = MPI_MAX
  ELSEIF (ytypop == 'MIN') THEN
    nzoper = MPI_MIN
  ELSE
    ierror  = 1
    yerrmsg = 'no valid operation type in global_values'
    RETURN
  ENDIF

  ! Get size of communicator icomm and rank of this PE
  CALL MPI_COMM_SIZE (icomm, nzsize, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_SIZE failed in global_values'
    RETURN
  ENDIF

  CALL MPI_COMM_RANK (icomm, myz_id, ierror)
  IF (ierror /= 0) THEN
    ierror  = 2
    yerrmsg = 'MPI_COMM_RANK failed in global_values'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize) THEN
    CALL MPI_ALLGATHER (rvector,     idim, IPU_real,                     &
                        rvector_out, idim, IPU_real, icomm, ierror)

    CALL MPI_ALLGATHER (ijmax,   2*idim, IPU_integer,                    &
                        indices, 2*idim, IPU_integer, icomm, ierror)
  ELSE
    CALL MPI_GATHER    (rvector,     idim, IPU_real,                     &
                        rvector_out, idim, IPU_real, ireceiver, icomm, ierror)

    CALL MPI_GATHER    (ijmax,   2*idim, IPU_integer,                    &
                        indices, 2*idim, IPU_integer, ireceiver, icomm, ierror)
  ENDIF
  IF (ierror /= 0) THEN
    ierror  = 3
    yerrmsg = 'MPI_GATHER failed in global_values'
    RETURN
  ENDIF

  IF (ireceiver < 0 .OR. ireceiver >= nzsize .OR. ireceiver == myz_id) THEN
    ! Determine the global maximas and the corresponding indices in
    ! the root process
    DO l = 1, idim
      ! Initialize variables
      zgmax   = rvector_out(l,1)
      izindex = 1
      DO n = 2,nzsize
        ! Determine the maxima and in which processor it is
        IF (rvector_out(l,n) > zgmax) THEN
          zgmax   = rvector_out(l,n)
          izindex = n
        ENDIF
      ENDDO

      ! Put Maximum into rvector and determine global indices
      rvector(l) = zgmax
      ijmax(1,l) = IPU_ipositions(izindex-1,1) + indices(1,l,izindex)   &
                                               - IPU_nbounds-1
      ijmax(2,l) = IPU_ipositions(izindex-1,2) + indices(2,l,izindex)   &
                                               - IPU_nbounds-1
    ENDDO
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'global_vectorreal_indices')
#endif

END SUBROUTINE global_vectorreal_indices

!==============================================================================
!+ Defines all subroutines for the generic routine combine_subarrays
!------------------------------------------------------------------------------
!
! SUBROUTINE combine_subarrays (procarray, dslocal)
!
! Description:
!  combine_subarrays puts the different subparts of a total field together
!  in a result field dslocal, which is processed by routines from the
!  grib library.
!    combine_[double,single] return a one-dimensional field dslocal
!    combine_2D_[double,single] return a two-dimensional field dslocal
!
!------------------------------------------------------------------------------
!+ put parts from subdomains together in a total field for double precision
!------------------------------------------------------------------------------

SUBROUTINE combine_2D_double (procarray, dslocal, lacc)

!------------------------------------------------------------------------------

! Parameterlist
REAL (KIND=dp),     INTENT(IN)  ::                              &
        procarray(IPU_idim_max,IPU_jdim_max,IPU_numpeco)
REAL (KIND=dp),     INTENT(OUT) ::                              &
        dslocal(IPU_idim_tot,IPU_jdim_tot)

LOGICAL, INTENT(IN), OPTIONAL ::                                &
        lacc                          ! if true assumes input and
                                      ! output on GPU

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER                               :: i
INTEGER                               :: iz_is, iz_ie, iz_js, iz_je
INTEGER                               :: jz_is, jz_ie, jz_js, jz_je
LOGICAL                               :: lzacc  ! logical for OpenACC
!
!- End of header
!==============================================================================

! Initialize GPU flag
IF (PRESENT(lacc)) THEN
  lzacc = lacc
ELSE
  lzacc = .FALSE.
ENDIF

DO i=0, IPU_numpeco-1
   IF ( i < IPU_numpey ) THEN
      iz_is = IPU_ipositions(i,1) - IPU_nbounds
   ELSE
      iz_is = IPU_ipositions(i,1)
   ENDIF
   IF ( mod(i,IPU_numpey)==0 ) THEN
      iz_js = IPU_ipositions(i,2) - IPU_nbounds
   ELSE
      iz_js = IPU_ipositions(i,2)
   ENDIF
   IF (( i >= IPU_numpeco-IPU_numpey ).AND.(i <= IPU_numpeco-1)) THEN
      iz_ie = IPU_ipositions(i,3) + IPU_nbounds
   ELSE
      iz_ie = IPU_ipositions(i,3)
   ENDIF
   IF ( mod(i+1,IPU_numpey)==0 ) THEN
      iz_je = IPU_ipositions(i,4) + IPU_nbounds
   ELSE
      iz_je = IPU_ipositions(i,4)
   ENDIF

   jz_is = iz_is - IPU_ipositions(i,1) + IPU_nbounds + 1
   jz_ie = iz_ie - IPU_ipositions(i,1) + IPU_nbounds + 1
   jz_js = iz_js - IPU_ipositions(i,2) + IPU_nbounds + 1
   jz_je = iz_je - IPU_ipositions(i,2) + IPU_nbounds + 1

   !$acc kernels present(dslocal,procarray) if (lzacc)
   dslocal(iz_is:iz_ie,iz_js:iz_je) = procarray(jz_is:jz_ie,jz_js:jz_je,i+1)
   !$acc end kernels

ENDDO

END SUBROUTINE combine_2D_double

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+ put parts from subdomains together in a total field for single precision
!------------------------------------------------------------------------------

SUBROUTINE combine_2D_single (procarray, dslocal, lacc)

!------------------------------------------------------------------------------

! Parameterlist
REAL (KIND=sp),     INTENT(IN)  ::                              &
        procarray(IPU_idim_max,IPU_jdim_max,IPU_numpeco)
REAL (KIND=sp),     INTENT(OUT) ::                              &
        dslocal(IPU_idim_tot,IPU_jdim_tot)

LOGICAL, INTENT(IN), OPTIONAL ::                                &
        lacc                          ! if true assumes input and
                                      ! output on GPU

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER                               :: i
INTEGER                               :: iz_is, iz_ie, iz_js, iz_je
INTEGER                               :: jz_is, jz_ie, jz_js, jz_je
LOGICAL                               :: lzacc  ! logical for OpenACC
!
!- End of header
!==============================================================================

! Initialize GPU flag
IF (PRESENT(lacc)) THEN
  lzacc = lacc
ELSE
  lzacc = .FALSE.
ENDIF

DO i=0, IPU_numpeco-1
   IF ( i < IPU_numpey ) THEN
      iz_is = IPU_ipositions(i,1) - IPU_nbounds
   ELSE
      iz_is = IPU_ipositions(i,1)
   ENDIF
   IF ( mod(i,IPU_numpey)==0 ) THEN
      iz_js = IPU_ipositions(i,2) - IPU_nbounds
   ELSE
      iz_js = IPU_ipositions(i,2)
   ENDIF
   IF (( i >= IPU_numpeco-IPU_numpey ).AND.(i <= IPU_numpeco-1)) THEN
      iz_ie = IPU_ipositions(i,3) + IPU_nbounds
   ELSE
      iz_ie = IPU_ipositions(i,3)
   ENDIF
   IF ( mod(i+1,IPU_numpey)==0 ) THEN
      iz_je = IPU_ipositions(i,4) + IPU_nbounds
   ELSE
      iz_je = IPU_ipositions(i,4)
   ENDIF

   jz_is = iz_is - IPU_ipositions(i,1) + IPU_nbounds + 1
   jz_ie = iz_ie - IPU_ipositions(i,1) + IPU_nbounds + 1
   jz_js = iz_js - IPU_ipositions(i,2) + IPU_nbounds + 1
   jz_je = iz_je - IPU_ipositions(i,2) + IPU_nbounds + 1

   !$acc kernels present(dslocal,procarray) if (lzacc)
   dslocal(iz_is:iz_ie,iz_js:iz_je) = procarray(jz_is:jz_ie,jz_js:jz_je,i+1)
   !$acc end kernels

ENDDO

END SUBROUTINE combine_2D_single

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+ put parts from subdomains together in a total field for double precision
!------------------------------------------------------------------------------

SUBROUTINE combine_double (procarray, dslocal)

!------------------------------------------------------------------------------

! Parameterlist
REAL (KIND=dp),    INTENT(IN)  ::                              &
        procarray(IPU_idim_max,IPU_jdim_max,IPU_numpeco)
REAL (KIND=dp),    INTENT(OUT) ::                              &
        dslocal(IPU_idim_tot*IPU_jdim_tot)

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER                               :: i, ids, jds, index_ds, ipa, jpa
INTEGER                               :: iz_is, iz_ie, iz_js, iz_je
INTEGER                               :: jz_is, jz_ie, jz_js, jz_je
!
!- End of header
!==============================================================================

DO i=0, IPU_numpeco-1
   IF ( i < IPU_numpey ) THEN
      iz_is = IPU_ipositions(i,1) - IPU_nbounds
   ELSE
      iz_is = IPU_ipositions(i,1)
   ENDIF
   IF ( mod(i,IPU_numpey)==0 ) THEN
      iz_js = IPU_ipositions(i,2) - IPU_nbounds
   ELSE
      iz_js = IPU_ipositions(i,2)
   ENDIF
   IF (( i >= IPU_numpeco-IPU_numpey ).AND.(i <= IPU_numpeco-1)) THEN
      iz_ie = IPU_ipositions(i,3) + IPU_nbounds
   ELSE
      iz_ie = IPU_ipositions(i,3)
   ENDIF
   IF ( mod(i+1,IPU_numpey)==0 ) THEN
      iz_je = IPU_ipositions(i,4) + IPU_nbounds
   ELSE
      iz_je = IPU_ipositions(i,4)
   ENDIF

   jz_is = iz_is - IPU_ipositions(i,1) + IPU_nbounds + 1
   jz_ie = iz_ie - IPU_ipositions(i,1) + IPU_nbounds + 1
   jz_js = iz_js - IPU_ipositions(i,2) + IPU_nbounds + 1
   jz_je = iz_je - IPU_ipositions(i,2) + IPU_nbounds + 1

   jds = iz_js - 1
   DO jpa = jz_js, jz_je
     jds = jds + 1
     ids = iz_is - 1
     DO ipa = jz_is, jz_ie
       ids = ids + 1
       index_ds = (jds-1) * IPU_idim_tot + ids
       dslocal (index_ds) = procarray(ipa,jpa,i+1)
     ENDDO
   ENDDO

!  dslocal(iz_is:iz_ie,iz_js:iz_je) = procarray(jz_is:jz_ie,jz_js:jz_je,i+1)

ENDDO

END SUBROUTINE combine_double

!==============================================================================
!------------------------------------------------------------------------------
!+ put parts from subdomains together in a total field for single precision
!------------------------------------------------------------------------------

SUBROUTINE combine_single (procarray, dslocal)

!------------------------------------------------------------------------------

! Parameterlist
REAL (KIND=sp),    INTENT(IN)  ::                              &
        procarray(IPU_idim_max,IPU_jdim_max,IPU_numpeco)
REAL (KIND=sp),    INTENT(OUT) ::                              &
        dslocal(IPU_idim_tot*IPU_jdim_tot)

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER                               :: i, ids, jds, index_ds, ipa, jpa
INTEGER                               :: iz_is, iz_ie, iz_js, iz_je
INTEGER                               :: jz_is, jz_ie, jz_js, jz_je
!
!- End of header
!==============================================================================

DO i=0, IPU_numpeco-1
   IF ( i < IPU_numpey ) THEN
      iz_is = IPU_ipositions(i,1) - IPU_nbounds
   ELSE
      iz_is = IPU_ipositions(i,1)
   ENDIF
   IF ( mod(i,IPU_numpey)==0 ) THEN
      iz_js = IPU_ipositions(i,2) - IPU_nbounds
   ELSE
      iz_js = IPU_ipositions(i,2)
   ENDIF
   IF (( i >= IPU_numpeco-IPU_numpey ).AND.(i <= IPU_numpeco-1)) THEN
      iz_ie = IPU_ipositions(i,3) + IPU_nbounds
   ELSE
      iz_ie = IPU_ipositions(i,3)
   ENDIF
   IF ( mod(i+1,IPU_numpey)==0 ) THEN
      iz_je = IPU_ipositions(i,4) + IPU_nbounds
   ELSE
      iz_je = IPU_ipositions(i,4)
   ENDIF

   jz_is = iz_is - IPU_ipositions(i,1) + IPU_nbounds + 1
   jz_ie = iz_ie - IPU_ipositions(i,1) + IPU_nbounds + 1
   jz_js = iz_js - IPU_ipositions(i,2) + IPU_nbounds + 1
   jz_je = iz_je - IPU_ipositions(i,2) + IPU_nbounds + 1

   jds = iz_js - 1
   DO jpa = jz_js, jz_je
     jds = jds + 1
     ids = iz_is - 1
     DO ipa = jz_is, jz_ie
       ids = ids + 1
       index_ds = (jds-1) * IPU_idim_tot + ids
       dslocal (index_ds) = procarray(ipa,jpa,i+1)
     ENDDO
   ENDDO

!  dslocal(iz_is:iz_ie,iz_js:iz_je) = procarray(jz_is:jz_ie,jz_js:jz_je,i+1)

ENDDO

END SUBROUTINE combine_single

!==============================================================================

!==============================================================================

FUNCTION i_global (i_loc)
!$acc routine seq

!------------------------------------------------------------------------------
!
! Description:
!  Function i_global determines the global i-index that belongs to the local
!  index i_loc and the subdomain my_cart_id. If i_loc is not in the range
!  1..ie, an error message is printed and the program is aborted.
!  Note that, if the grid point belonging to i_loc is located at the boundary
!  of a subdomain, it also belongs to (several) other subdomains.
!
! Method:
!  Using isubpos from Module data_parallel.
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                 , INTENT (IN) ::  &
    i_loc     ! local i-index of this subdomain

  INTEGER                               ::  &
    i_global

!------------ End of header ---------------------------------------------------

  ! Test for wrong values
  IF (i_loc < 1 .OR. i_loc > IPU_idim) THEN
#ifndef _OPENACC
    !print not supported in acc routine
    PRINT *, ' *** Error in Function i_global *** '
    PRINT *, ' *** Index ',i_loc,' is not in the range 1..',IPU_idim
#endif
    i_global = -9999
  ELSE
    ! Compute global index
    i_global = IPU_ipositions(IPU_cart_id, 1) - IPU_nbounds - 1 + i_loc
  ENDIF

END FUNCTION i_global

!==============================================================================

!==============================================================================

FUNCTION j_global (j_loc)
!$acc routine seq

!------------------------------------------------------------------------------
!
! Description:
!  Function j_global determines the global i-index that belongs to the local
!  index j_loc and the subdomain my_cart_id. If j_loc is not in the range
!  1..ie, an error message is printed and the program is aborted.
!  Note that, if the grid point belonging to j_loc is located at the boundary
!  of a subdomain, it also belongs to (several) other subdomains.
!
! Method:
!  Using isubpos from Module data_parallel.
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                 , INTENT (IN) ::  &
    j_loc     ! local j-index of this subdomain


  INTEGER                               ::  &
    j_global

!------------ End of header ---------------------------------------------------

  ! Test for wrong values
  IF (j_loc < 1 .OR. j_loc > IPU_jdim) THEN
#ifndef _OPENACC
    !print not supported in acc routine
    PRINT *, ' *** Error in Function j_global *** '
    PRINT *, ' *** Index ',j_loc,' is not in the range 1..',IPU_jdim
#endif
    j_global = -9999
  ELSE
    ! Compute global index
    j_global = IPU_ipositions(IPU_cart_id, 2) - IPU_nbounds - 1 + j_loc
  ENDIF

END FUNCTION j_global

!==============================================================================

!==============================================================================

SUBROUTINE ij_local (i_glob, j_glob, i_loc, j_loc, isubdomain, ierror)
!$acc routine seq

!------------------------------------------------------------------------------
!
! Description:
!  Given the global i- and j-indices of a point, subroutine ij_local determines
!  the subdomain, to which this grid point belongs and the local i- and
!  j-indices in this subdomain. If the grid point (i_glob,j_glob) belongs to
!  several subdomains, only the unique subdomain where it belongs to the
!  interior is considered.
!
! Method:
!  Using isubpos from Module data_parallel.
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                 , INTENT (IN) ::  &
    i_glob,  & ! local i-index of this subdomain
    j_glob     ! local j-index of this subdomain


  INTEGER                 , INTENT(OUT) ::  &
    i_loc,   & ! global i-index
    j_loc,   & ! global j-index
    isubdomain ! number of the processor in the cartesian communicator

! Local variables:
! ----------------

  INTEGER                               ::  &
    ierror, i, j, ix, jy, ilow, ihig, jlow, jhig

!------------ End of header ---------------------------------------------------

  ! Test for wrong values
  ierror = 0
  IF (i_glob < 1 .OR. i_glob > IPU_idim_tot) THEN
#ifndef _OPENACC
    !print not supported in acc routine
    PRINT *, ' *** Error in Routine ij_local *** '
    PRINT *, ' *** Index ',i_glob,' is not in the range 1..',IPU_idim_tot
#endif
    ierror = 1
  ENDIF
  IF (j_glob < 1 .OR. j_glob > IPU_jdim_tot) THEN
#ifndef _OPENACC
    PRINT *, ' *** Error in Routine ij_local *** '
    PRINT *, ' *** Index ',j_glob,' is not in the range 1..',IPU_jdim_tot
#endif
    ierror = 1
  ENDIF

  ! Only execute this part if no error
  ! Return statement not supported in acc routine
  IF (ierror == 0) THEN

    ! Search from 0..IPU_numpex-1 to determine the subdomain in x-direction to
    ! which i_glob belongs
    DO i = 0, IPU_numpex-1
      IF (i == 0) THEN
        ilow = IPU_ipositions(i*IPU_numpey,1)-IPU_nbounds
      ELSE
        ilow = IPU_ipositions(i*IPU_numpey,1)
      ENDIF
      IF (i == IPU_numpex-1) THEN
        ihig = IPU_ipositions(i*IPU_numpey,3)+IPU_nbounds
      ELSE
        ihig = IPU_ipositions(i*IPU_numpey,3)
      ENDIF
      IF ( (ilow <= i_glob) .AND. (i_glob <= ihig) ) THEN
        ix = i
        EXIT
      ENDIF
    ENDDO

    ! Search from 0..nprocy-1 to determine the subdomain in y-direction to
    ! which j_glob belongs
    DO j = 0, IPU_numpey-1
      IF (j == 0) THEN
        jlow = IPU_ipositions(j,2)-IPU_nbounds
      ELSE
        jlow = IPU_ipositions(j,2)
      ENDIF
      IF (j == IPU_numpey-1) THEN
        jhig = IPU_ipositions(j,4)+IPU_nbounds
      ELSE
        jhig = IPU_ipositions(j,4)
      ENDIF
      IF ( (jlow <= j_glob) .AND. (j_glob <= jhig) ) THEN
        jy = j
        EXIT
      ENDIF
    ENDDO

    ! Set the output variables
    isubdomain = ix * IPU_numpey + jy
    i_loc      = 1 + IPU_nbounds + (i_glob - IPU_ipositions(isubdomain,1))
    j_loc      = 1 + IPU_nbounds + (j_glob - IPU_ipositions(isubdomain,2))
  ELSE
    ! Set to unrealistic values in case of error
    isubdomain = -9999
    i_loc      = -9999
    j_loc      = -9999
  ENDIF

END SUBROUTINE ij_local

!==============================================================================

!+ To print a point of a real 2D field in the parallel program
!------------------------------------------------------------------------------

SUBROUTINE print_point (i_glob, j_glob, rfield, idim, jdim, ytext,     &
                             iout, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  The subroutine print_point prints the point (i_glob,j_glob) from a
!  two-dimensional model field. i_glob and j_glob have to be the global
!  indices of the grid point, i.e. the indices regarding the total domain.
!  For additional output information a message up to 60 characters can
!  be printed by using the character variable ytext.
!  With the variable iout it is possible to direct the output either to
!  standard output (iout=0) or to a file (30 <= iout <= 49) with unit number
!  iout. This file has to be opened by the user before calling this subroutine.
!
! Method:
!  With routine ij_local it is determined which processor owns the grid point
!  with the given global indices (i_glob,j_glob).
!  If iout=0, the ytext and the grid point value is printed to standard output,
!  otherwise the information is sent to PE 0 to write it to the file iout.
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  i_glob, j_glob,     & ! global indices of the point to be printed
  idim , jdim,        & ! dimensions of rfield
  iout                  ! unit number of file where to write to

! Array arguments with intent(inout):
REAL    (KIND=wp),        INTENT(IN)       ::                              &
  rfield(idim,jdim)     ! field from which a point is printed

! Character arguments with intent(in):
CHARACTER (LEN=*), INTENT(IN)              ::                              &
  ytext  ! for additional output information

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                ! error status variable

!------------------------------------------------------------------------------

#ifndef NOMPI
! Local Variables
INTEGER                                    ::                              &
  i_loc, j_loc,       & ! local indices of the point to be printed
  isub,               & ! subdomain that contains the point to be printed
  ilentext,           & ! lenght of ytext
  iztag,              & ! tag of message
  izrequest(MPI_STATUS_SIZE), & ! for MPI_RECV
  implcode              ! local error code

CHARACTER (LEN=60)                         ::                              &
  yline                 ! storing ytext and value of field in character format

REAL    (KIND=wp)                          ::                              &
  rdata                 ! value of the field at grid point i_glob, j_glob

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations and Checks
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  iztag    = 8001
  yline    = ytext
  ilentext = LEN (ytext)

  ! Check the input values
  IF (idim /= IPU_idim .OR. jdim /= IPU_jdim) THEN
    PRINT *, ' WARNING: Wrong input values for print_point:'
    PRINT *, '      idim = ',idim,' and ie = ',IPU_idim
    PRINT *, '      jdim = ',jdim,' and je = ',IPU_jdim
    ierror = 1
  ENDIF

  IF (iout /= 0 .AND. (iout < 30 .OR. iout > 49)) THEN
    PRINT *, ' WARNING: Wrong input values for print_point: iout = ', iout
    ierror = 2
  ENDIF

  IF (i_glob < 1 .OR. i_glob > IPU_idim_tot) THEN
    PRINT *, ' WARNING: Wrong input values for print_point:'
    PRINT *, '      i_glob = ',i_glob,' and ie = ',IPU_idim_tot
    ierror = 3
  ENDIF

  IF (j_glob < 1 .OR. j_glob > IPU_jdim_tot) THEN
    PRINT *, '      j_glob = ',j_glob,' and je = ',IPU_jdim_tot
    ierror = 4
  ENDIF

  IF (ierror /= 0) RETURN

!------------------------------------------------------------------------------
! Section 2: Determine the subdomain with i_glob, j_glob
!------------------------------------------------------------------------------

  CALL ij_local (i_glob, j_glob, i_loc, j_loc, isub, ierror)
  IF (ierror /= 0) RETURN
  IF (isub == IPU_cart_id) THEN
    rdata = rfield (i_loc,j_loc)
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Print text and value; if necessary, send data to PE 0
!------------------------------------------------------------------------------

  IF (iout == 0) THEN
    IF (isub == IPU_cart_id) THEN
      ! subdomain which contains the grid point, just prints it
      WRITE (*,*) yline(1:ilentext), rdata
    ENDIF
  ELSE
    IF (isub /= 0) THEN
      ! Send data to PE 0
      IF (isub == IPU_cart_id) THEN
        CALL MPI_SEND (rdata, 1, IPU_real,    0, iztag, IPU_comm_cart,    &
                       implcode)
      ELSEIF (IPU_cart_id == 0) THEN
        CALL MPI_RECV (rdata, 1, IPU_real, isub, iztag, IPU_comm_cart,    &
                       izrequest, implcode)
        WRITE (iout,*) yline(1:ilentext), rdata
      ENDIF
    ELSE
      ! PE 0 has the data
      WRITE (iout,*) yline(1:ilentext), rdata
    ENDIF
  ENDIF

  CALL MPI_BARRIER (IPU_comm_cart, implcode)
#else
  CALL model_abort (-1, -1, 'MPI not available', 'print_point')
#endif

END SUBROUTINE print_point

!==============================================================================

!==============================================================================

!+ Defines all subroutines for the generic routine print_column
!------------------------------------------------------------------------------

SUBROUTINE print_column (i_glob, j_glob, rfield, idim, jdim, kdim, ytext, &
                              iout, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  The subroutine print_column prints the column over the point (i_glob,j_glob)
!  from a three-dimensional model field. i_glob and j_glob have to be the
!  global indices of the grid point, i.e. the indices regarding the total
!  domain. For additional output information a message up to 60 characters can
!  be printed by using the character variable ytext.
!  With the variable iout it is possible to direct the output either to
!  standard output (iout=0) or to a file (30 <= iout <= 49) with unit number
!  iout. This file has to be opened by the user before calling this subroutine.
!
! Method:
!  With routine ij_local it is determined which processor owns the grid point
!  with the given global indices (i_glob,j_glob).
!  If iout=0, the ytext and the grid point value is printed to standard output,
!  otherwise the information is sent to PE 0 to write it to the file iout.
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  i_glob, j_glob,     & ! global indices of the point to be printed
  idim , jdim, kdim,  & ! dimensions of rfield
  iout                  ! unit number of file where to write to

! Array arguments with intent(inout):
REAL    (KIND=wp),        INTENT(IN)       ::                              &
  rfield(idim,jdim,kdim)! field from which a point is printed

! Character arguments with intent(in):
CHARACTER (LEN=*), INTENT(IN)              ::                              &
  ytext  ! for additional output information

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                ! error status variable

!------------------------------------------------------------------------------

#ifndef NOMPI
! Local Variables
INTEGER                                    ::                              &
  i_loc, j_loc,       & ! local indices of the point to be printed
  isub,               & ! subdomain that contains the point to be printed
  ilentext,           & ! lenght of ytext
  iztag,              & ! tag of message
  izrequest(MPI_STATUS_SIZE), & ! for MPI_RECV
  implcode, k           ! local error code

CHARACTER (LEN=60)                         ::                              &
  yline                 ! storing ytext and value of field in character format

REAL    (KIND=wp)                          ::                              &
  rdata(kdim)           ! value of the field at grid point i_glob, j_glob

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine print_real_column
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations and Checks
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0
  iztag    = 8001
  yline    = ytext
  ilentext = LEN (ytext)

  ! Check the input values
  IF (idim /= IPU_idim .OR. jdim /= IPU_jdim .OR. kdim /= IPU_kdim) THEN
    PRINT *, ' WARNING: Wrong input values for print_point:'
    PRINT *, '      idim = ',idim,' and ie = ',IPU_idim
    PRINT *, '      jdim = ',jdim,' and je = ',IPU_jdim
    PRINT *, '      kdim = ',kdim,' and ke = ',IPU_kdim
    ierror = 1
  ENDIF

  IF (iout /= 0 .AND. (iout < 30 .OR. iout > 49)) THEN
    PRINT *, ' WARNING: Wrong input values for print_point: iout = ', iout
    ierror = 2
  ENDIF

  IF (i_glob < 1 .OR. i_glob > IPU_idim_tot) THEN
    PRINT *, ' WARNING: Wrong input values for print_point:'
    PRINT *, '      i_glob = ',i_glob,' and ie = ',IPU_idim_tot
    ierror = 3
  ENDIF

  IF (j_glob < 1 .OR. j_glob > IPU_jdim_tot) THEN
    PRINT *, '      j_glob = ',j_glob,' and je = ',IPU_jdim_tot
    ierror = 4
  ENDIF

  IF (ierror /= 0) RETURN

!------------------------------------------------------------------------------
! Section 2: Determine the subdomain with i_glob, j_glob
!------------------------------------------------------------------------------

  CALL ij_local (i_glob, j_glob, i_loc, j_loc, isub, ierror)
  IF (ierror /= 0) RETURN
  IF (isub == IPU_cart_id) THEN
    DO k = 1, kdim
      rdata(k) = rfield (i_loc,j_loc,k)
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Print text and value; if necessary, send data to PE 0
!------------------------------------------------------------------------------

  IF (iout == 0) THEN
    IF (isub == IPU_cart_id) THEN
      ! subdomain which contains the grid point, just prints it
      WRITE (*,*) yline(1:ilentext)
      DO k = 1, kdim
        WRITE (*,*) '  k = ',k,':  ',rdata(k)
      ENDDO
    ENDIF
  ELSE
    IF (isub /= 0) THEN
      ! Send data to PE 0
      IF (isub == IPU_cart_id) THEN
        CALL MPI_SEND (rdata, kdim, IPU_real,    0, iztag, IPU_comm_cart,    &
                       implcode)
      ELSEIF (IPU_cart_id == 0) THEN
        CALL MPI_RECV (rdata, kdim, IPU_real, isub, iztag, IPU_comm_cart,    &
                       izrequest, implcode)
      ENDIF
    ENDIF

    ! PE 0 now has the data
    IF (IPU_cart_id == 0) THEN
      WRITE (iout,*) yline(1:ilentext)
      DO k = 1, kdim
        WRITE (iout,*) '  k = ',k,':  ',rdata(k)
      ENDDO
    ENDIF
  ENDIF

  CALL MPI_BARRIER (IPU_comm_cart, implcode)
#else
  CALL model_abort (-1, -1, 'MPI not available', 'print_column')
#endif

END SUBROUTINE print_column

!==============================================================================

!==============================================================================

!+ Defines all subroutines for the generic routine distribute_field
!------------------------------------------------------------------------------
!
!SUBROUTINE distribute_field (rfield_tot, idim_tot, jdim_tot,              &
!                             rfield_loc, idim_loc, jdim_loc, isender, ierror)
!
!------------------------------------------------------------------------------
!
! Description:
!  This routine takes a two-dimensional total domain field, decomposes it into
!  the parts for the subdomains and distributes these parts to the
!  corresponding processors in the cartesian communicator group. The routine
!  has to be called by all processors in this group. The total field has to
!  be in processor isender.
!
! Method:
!  The subparts of the total domain field are copied into a special field
!  procarray, which is distributed using MPI_SCATTER.
!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! double precision version:
!------------------------------------------------------------------------------

SUBROUTINE distribute_field_double (rfield_tot, idim_tot, jdim_tot,        &
                             rfield_loc, idim_loc, jdim_loc, isender, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim_tot, jdim_tot, & ! dimensions of the total field
  idim_loc, jdim_loc, & ! dimensions of the subdomain
  isender               ! PE that contains the whole field and sends the parts

REAL    (KIND=dp),        INTENT(IN)       ::                              &
  rfield_tot(idim_tot,jdim_tot) ! total domain field that has to be distributed

! Scalar arguments with intent(out):
REAL    (KIND=dp),        INTENT(OUT)      ::                              &
  rfield_loc(idim_loc,jdim_loc) ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                        ! error status variable

!------------------------------------------------------------------------------

! Local Variables
INTEGER                                    ::                              &
  ie_p, je_p,             & ! to adjust dimensions of the different subdomains
  izlo, izup, jzlo, jzup, & ! lower and upper i- and j-dimensions
  implcode, n               ! local error code

REAL    (KIND=dp)                          ::                              &
  rzprocarray(IPU_idim_max,IPU_jdim_max,0:IPU_numpex*IPU_numpey-1),        &
                        ! holds all subdomains of total field
  rzsubarray (IPU_idim_max,IPU_jdim_max)
                        ! subdomain for this processor

!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
!------------------------------------------------------------------------------
! Section 1: Initializations and Checks
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0

  ! Check the input values
  IF (idim_tot /= IPU_idim_tot .OR. jdim_tot /= IPU_jdim_tot) THEN
    PRINT *, ' WARNING: Wrong input values for distribute_field:'
    PRINT *, '      idim_tot = ',idim_tot,' and ie_tot = ',IPU_idim_tot
    PRINT *, '      jdim_tot = ',jdim_tot,' and je_tot = ',IPU_jdim_tot
    ierror = 1
  ENDIF

  IF (idim_loc /= IPU_idim .OR. jdim_loc /= IPU_jdim) THEN
    PRINT *, ' WARNING: Wrong input values for distribute_field:'
    PRINT *, '      idim_loc = ',idim_loc,' and ie = ',IPU_idim
    PRINT *, '      jdim_loc = ',jdim_loc,' and je = ',IPU_jdim
    ierror = 1
  ENDIF

  IF (ierror /= 0) RETURN

  ! This subroutine has to be called by all routines in the cartesian grid
  CALL MPI_BARRIER (IPU_comm_cart, implcode)

!------------------------------------------------------------------------------
! Section 2: Copy total field into parts of procarray in processor 0
!------------------------------------------------------------------------------

  IF (IPU_cart_id == isender) THEN
    DO n = 0, IPU_numpex*IPU_numpey-1
      izlo = IPU_ipositions(n,1)
      izup = IPU_ipositions(n,3)
      jzlo = IPU_ipositions(n,2)
      jzup = IPU_ipositions(n,4)


      ie_p = izup - izlo + 1 + 2*IPU_nbounds
      je_p = jzup - jzlo + 1 + 2*IPU_nbounds
      rzprocarray(1:ie_p,1:je_p,n) =                       &
         rfield_tot(izlo-IPU_nbounds : izup+IPU_nbounds,   &
                    jzlo-IPU_nbounds : jzup+IPU_nbounds)
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Distribute procarray
!------------------------------------------------------------------------------

  IF (IPU_numpex*IPU_numpey > 1) THEN
    CALL MPI_SCATTER ( rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_double,   &
                       rzsubarray , IPU_idim_max*IPU_jdim_max, IPU_double,   &
                       isender, IPU_comm_cart, implcode)
    IF (implcode /= 0) THEN
      ierror = 2
      RETURN
    ENDIF
  ELSE
    rzsubarray(:,:) = rzprocarray(:,:,0)
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Copy received subpart into field_loc in all processors
!------------------------------------------------------------------------------

  rfield_loc(1:IPU_idim,1:IPU_jdim) = rzsubarray(1:IPU_idim,1:IPU_jdim)
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_field')
#endif

END SUBROUTINE distribute_field_double

!------------------------------------------------------------------------------
! single precision version:
!------------------------------------------------------------------------------

SUBROUTINE distribute_field_single (rfield_tot, idim_tot, jdim_tot,        &
                             rfield_loc, idim_loc, jdim_loc, isender, ierror)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim_tot, jdim_tot, & ! dimensions of the total field
  idim_loc, jdim_loc, & ! dimensions of the subdomain
  isender               ! PE that contains the whole field and sends the parts

REAL    (KIND=sp),        INTENT(IN)       ::                              &
  rfield_tot(idim_tot,jdim_tot) ! total domain field that has to be distributed

! Scalar arguments with intent(out):
REAL    (KIND=sp),        INTENT(OUT)      ::                              &
  rfield_loc(idim_loc,jdim_loc) ! subdomain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                        ! error status variable

!------------------------------------------------------------------------------

! Local Variables
INTEGER                                    ::                              &
  ie_p, je_p,             & ! to adjust dimensions of the different subdomains
  izlo, izup, jzlo, jzup, & ! lower and upper i- and j-dimensions
  implcode, n               ! local error code

REAL    (KIND=sp)                          ::                              &
  rzprocarray(IPU_idim_max,IPU_jdim_max,0:IPU_numpex*IPU_numpey-1),        &
                        ! holds all subdomains of total field
  rzsubarray (IPU_idim_max,IPU_jdim_max)
                        ! subdomain for this processor

!- End of header
!------------------------------------------------------------------------------

#ifndef NOMPI
!------------------------------------------------------------------------------
! Section 1: Initializations and Checks
!------------------------------------------------------------------------------

  ! Initializations
  ierror   = 0

  ! Check the input values
  IF (idim_tot /= IPU_idim_tot .OR. jdim_tot /= IPU_jdim_tot) THEN
    PRINT *, ' WARNING: Wrong input values for distribute_field:'
    PRINT *, '      idim_tot = ',idim_tot,' and ie_tot = ',IPU_idim_tot
    PRINT *, '      jdim_tot = ',jdim_tot,' and je_tot = ',IPU_jdim_tot
    ierror = 1
  ENDIF

  IF (idim_loc /= IPU_idim .OR. jdim_loc /= IPU_jdim) THEN
    PRINT *, ' WARNING: Wrong input values for distribute_field:'
    PRINT *, '      idim_loc = ',idim_loc,' and ie = ',IPU_idim
    PRINT *, '      jdim_loc = ',jdim_loc,' and je = ',IPU_jdim
    ierror = 1
  ENDIF

  IF (ierror /= 0) RETURN

  ! This subroutine has to be called by all routines in the cartesian grid
  CALL MPI_BARRIER (IPU_comm_cart, implcode)

!------------------------------------------------------------------------------
! Section 2: Copy total field into parts of procarray in processor 0
!------------------------------------------------------------------------------

  IF (IPU_cart_id == isender) THEN
    DO n = 0, IPU_numpex*IPU_numpey-1
      izlo = IPU_ipositions(n,1)
      izup = IPU_ipositions(n,3)
      jzlo = IPU_ipositions(n,2)
      jzup = IPU_ipositions(n,4)


      ie_p = izup - izlo + 1 + 2*IPU_nbounds
      je_p = jzup - jzlo + 1 + 2*IPU_nbounds
      rzprocarray(1:ie_p,1:je_p,n) =                       &
         rfield_tot(izlo-IPU_nbounds : izup+IPU_nbounds,   &
                    jzlo-IPU_nbounds : jzup+IPU_nbounds)
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Distribute procarray
!------------------------------------------------------------------------------

  IF (IPU_numpex*IPU_numpey > 1) THEN
    CALL MPI_SCATTER ( rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_single,   &
                       rzsubarray , IPU_idim_max*IPU_jdim_max, IPU_single,   &
                       isender, IPU_comm_cart, implcode)
    IF (implcode /= 0) THEN
      ierror = 2
      RETURN
    ENDIF
  ELSE
    rzsubarray(:,:) = rzprocarray(:,:,0)
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Copy received subpart into field_loc in all processors
!------------------------------------------------------------------------------

  rfield_loc(1:IPU_idim,1:IPU_jdim) = rzsubarray(1:IPU_idim,1:IPU_jdim)
#else
  CALL model_abort (-1, -1, 'MPI not available', 'distribute_field')
#endif

END SUBROUTINE distribute_field_single

!==============================================================================
!==============================================================================
!+ dumps  LM fields for the total domain
!------------------------------------------------------------------------------

SUBROUTINE dump_lm_field                                                 &
   (rfield_loc, yname, ydimension, ilev, factor, bias,                   &
    istart, iend, jstart, jend, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  This routine gathers the local parts of a total field from all subdomains
!  in the cartesian communicator to processor 0 and prints the part of the
!  total field specified by (istart:iend,jstart:jend) scaled with factor and
!  bias. It has to be called by all processors in this group.
!
!  Caution:
!  The routines specified for the generic routine print_lm_field allocate
!  memory of their own. This has not been avoided because this is regarded
!  to be a routine used only for debug purposes, so that performance and
!  memory usage are not very important.
!
! Method:
!  Call to gather_field
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::     &
  ilev,               & ! level of field_loc in the atmosphere
  istart, iend,       & ! start and end-indices for the part of the
  jstart, jend          ! field that has to be printed

REAL (KIND=wp),           INTENT(IN)       ::     &
  rfield_loc(IPU_idim,IPU_jdim) ! subdomain field

REAL (KIND=wp),           INTENT(IN)       ::     &
  factor, bias          ! for scaling the output according to GRIB

CHARACTER (LEN=*),        INTENT(IN)       ::     &
  yname,              & ! name of the variable
  ydimension            ! dimension (unit) of the variable

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT)      ::     &
  ierror                        ! error status variable

CHARACTER (LEN=*)                          ::     &
  yerrmsg               ! for error string

!------------------------------------------------------------------------------

! Local Variables
INTEGER                                    ::                              &
  istat, implcode, n, i1, i2, i, j, nout    ! local error code

REAL (KIND=wp)                             ::                              &
  rzfield_tot(IPU_idim_tot,IPU_jdim_tot)     ! total field

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations and Checks
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  istat    = 0
  ierror   = 0
  yerrmsg  = '   '

  ! This subroutine has to be called by all routines in icomm
  IF (IPU_numpeco > 1) THEN
    CALL MPI_BARRIER (IPU_comm_cart, implcode)
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Gather the field
!------------------------------------------------------------------------------

  IF (IPU_numpeco > 1) THEN
    CALL gather_field (rfield_loc,  IPU_idim, IPU_jdim,                       &
                       rzfield_tot, IPU_idim_tot, IPU_jdim_tot, 0, ierror)
  ELSE
    rzfield_tot(:,:) = rfield_loc(:,:)
  ENDIF

  IF (ierror /= 0) RETURN

!------------------------------------------------------------------------------
! Section 3: Print the field
!------------------------------------------------------------------------------

  IF (IPU_cart_id == 0) THEN
    nout = 47
    OPEN  (nout, FILE=yname, FORM='FORMATTED', STATUS='UNKNOWN', IOSTAT=istat)
    rzfield_tot(:,:) = (rzfield_tot(:,:) + bias) * factor
    WRITE (nout,'(A,A,A,I3,A,A,A,F8.3,A,F8.1)')                           &
      '  Field: ', yname, '  Level: ', ilev, '  Dimension: ', ydimension, &
      '  Factor: ', 1.0_wp/factor, '  Bias: ', bias
    DO j = jstart, jend
      DO i = istart, iend
        WRITE (nout,'(I3,I3,F40.25)') i,j, rzfield_tot(i,j)
      ENDDO
    ENDDO
    CLOSE (nout, STATUS='KEEP')
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'dump_lm_field')
#endif

END SUBROUTINE dump_lm_field

!==============================================================================

!==============================================================================
!+ Defines all subroutines for the generic routine gather_field
!------------------------------------------------------------------------------
!
!SUBROUTINE gather_field (rfield_loc, idim_loc, jdim_loc,                  &
!                         rfield_tot, idim_tot, jdim_tot, ireceiver, ierror)
!
!------------------------------------------------------------------------------
!
! Description:
!  This routine gathers the parts of a total field from all subdomains in the
!  cartesian communicator group to processor ireceiver. If ireceiver < 0, it
!  will be gathered to all PEs in this group. It has to be called by all
!  processors in this group.
!
! Method:
!  The subparts of the total domain field are copied into a special field
!  subarray, which is gathered using MPI_GATHER.
!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! double precision version:
!------------------------------------------------------------------------------

SUBROUTINE gather_field_double (rfield_loc, idim_loc, jdim_loc,            &
                                rfield_tot, idim_tot, jdim_tot,            &
                                ireceiver, ierror, lacc)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim_tot, jdim_tot, & ! dimensions of the total field
  idim_loc, jdim_loc, & ! dimensions of the subdomain
  ireceiver             ! PE that should get the field. If ireceiver < 0
                        ! the field is gathered to all PEs

REAL    (KIND=dp),        INTENT(IN)       ::                              &
  rfield_loc(idim_loc,jdim_loc) ! subdomain field

! Scalar arguments with intent(out):
REAL    (KIND=dp),        INTENT(OUT)      ::                              &
  rfield_tot(idim_tot,jdim_tot) ! total domain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                        ! error status variable

LOGICAL, INTENT(IN), OPTIONAL              ::                              &
  lacc                          ! if true assumes input and
                                ! output on GPU

!------------------------------------------------------------------------------

! Local Variables
INTEGER                                    ::                              &
  izlo, izup, jzlo, jzup,                 & ! indices for local subdomains
  izlo_tot, izup_tot, jzlo_tot, jzup_tot, & ! indices for total domain
  implcode, n                               ! local error code

REAL    (KIND=dp)                          ::                              &
  rzprocarray(IPU_idim_max,IPU_jdim_max,0:IPU_numpeco-1), &
                                    ! holds all subdomains of total field
  rzsubarray (IPU_idim_max,IPU_jdim_max)    ! subdomain for this processor

LOGICAL :: lzacc      ! local logical for gpu execution

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations and Checks
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  ierror   = 0

  ! Initialize GPU flag
  IF (PRESENT(lacc)) THEN
    lzacc = lacc
  ELSE
    lzacc = .FALSE.
  ENDIF

  ! Check the input values
  IF (idim_tot /= IPU_idim_tot .OR. jdim_tot /= IPU_jdim_tot) THEN
    PRINT *, ' WARNING: Wrong input values for gather_field:'
    PRINT *, '      idim_tot = ',idim_tot,' and ie_tot = ',IPU_idim_tot
    PRINT *, '      jdim_tot = ',jdim_tot,' and je_tot = ',IPU_jdim_tot
    ierror = 1
  ENDIF

  IF (idim_loc /= IPU_idim .OR. jdim_loc /= IPU_jdim) THEN
    PRINT *, ' WARNING: Wrong input values for gather_field:'
    PRINT *, '      idim_loc = ',idim_loc,' and ie = ',IPU_idim
    PRINT *, '      jdim_loc = ',jdim_loc,' and je = ',IPU_jdim
    ierror = 1
  ENDIF

  IF (ierror /= 0) RETURN

  !$acc enter data create (rzprocarray,rzsubarray)

!------------------------------------------------------------------------------
! Section 2: Copy subpart into izsubarray
!------------------------------------------------------------------------------

  !$acc kernels present(rfield_loc,rzsubarray) if (lzacc)
  rzsubarray(1:IPU_idim,1:IPU_jdim) = rfield_loc(1:IPU_idim,1:IPU_jdim)
  !$acc end kernels

!------------------------------------------------------------------------------
! Section 3: Gather the subdomains
!------------------------------------------------------------------------------

  IF (IPU_numpeco > 1) THEN
    IF (ireceiver < 0) THEN

      !XL: We need to duplicate the call here because if statement is not
      !    supported with OpenACC host_data directive
      IF (lzacc) THEN
#ifdef NO_MPI_HOST_DATA
        !$acc update host(rzsubarray)
#else
        !$acc host_data use_device(rzsubarray,rzprocarray)
#endif
        CALL MPI_ALLGATHER (rzsubarray,  IPU_idim_max*IPU_jdim_max, IPU_double,&
                          rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_double,&
                          IPU_comm_cart,implcode)
#ifdef NO_MPI_HOST_DATA
        !$acc update device(rzprocarray)
#else
        !$acc end host_data
#endif
      ELSE
        CALL MPI_ALLGATHER (rzsubarray,  IPU_idim_max*IPU_jdim_max, IPU_double,&
                          rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_double,&
                          IPU_comm_cart,implcode)
      ENDIF

    ELSEIF (ireceiver < IPU_numpeco) THEN

      !XL: We need to duplicate the call here because if statement is not
      !    supported with OpenACC host_data directive
      IF (lzacc) THEN
#ifdef NO_MPI_HOST_DATA
        !$acc update host(rzsubarray)
#else
        !$acc host_data use_device(rzsubarray,rzprocarray)
#endif
        CALL MPI_GATHER ( rzsubarray , IPU_idim_max*IPU_jdim_max, IPU_double,  &
             rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_double,  &
             ireceiver  , IPU_comm_cart, implcode)
#ifdef NO_MPI_HOST_DATA
        !$acc update device(rzprocarray)
#else
        !$acc end host_data
#endif
      ELSE
        CALL MPI_GATHER ( rzsubarray , IPU_idim_max*IPU_jdim_max, IPU_double,  &
             rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_double,  &
             ireceiver  , IPU_comm_cart, implcode)
      ENDIF

      IF (implcode /= 0) THEN
        ierror = 2
        RETURN
      ENDIF

    ELSE
      ierror  = 3
      RETURN
    ENDIF
  ELSE
    !$acc kernels present(rzprocarray,rzsubarray) if(lzacc)
    rzprocarray(:,:,0) = rzsubarray(:,:)
    !$acc end kernels
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Combine the different parts to the total field
!------------------------------------------------------------------------------

  IF (ireceiver < 0) THEN
    CALL combine_subarrays (rzprocarray, rfield_tot, lzacc)
  ELSEIF (ireceiver < IPU_numpeco) THEN
    IF (IPU_cart_id == ireceiver) THEN
      CALL combine_subarrays (rzprocarray, rfield_tot, lzacc)
    ENDIF
  ENDIF

  !$acc exit data delete(rzprocarray,rzsubarray)

#else
  CALL model_abort (-1, -1, 'MPI not available', 'gather_field')
#endif

END SUBROUTINE gather_field_double

!------------------------------------------------------------------------------
! single precision version:
!------------------------------------------------------------------------------

SUBROUTINE gather_field_single (rfield_loc, idim_loc, jdim_loc,            &
                                rfield_tot, idim_tot, jdim_tot,            &
                                ireceiver, ierror, lacc)

!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  idim_tot, jdim_tot, & ! dimensions of the total field
  idim_loc, jdim_loc, & ! dimensions of the subdomain
  ireceiver             ! PE that should get the field. If ireceiver < 0
                        ! the field is gathered to all PEs

REAL    (KIND=sp),        INTENT(IN)       ::                              &
  rfield_loc(idim_loc,jdim_loc) ! subdomain field

! Scalar arguments with intent(out):
REAL    (KIND=sp),        INTENT(OUT)      ::                              &
  rfield_tot(idim_tot,jdim_tot) ! total domain field

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                        ! error status variable

LOGICAL, INTENT(IN), OPTIONAL              ::                              &
  lacc                          ! if true assumes input and
                                ! output on GPU

!------------------------------------------------------------------------------

! Local Variables
INTEGER                                    ::                              &
  izlo, izup, jzlo, jzup,                 & ! indices for local subdomains
  izlo_tot, izup_tot, jzlo_tot, jzup_tot, & ! indices for total domain
  implcode, n                               ! local error code

REAL    (KIND=sp)                          ::                              &
  rzprocarray(IPU_idim_max,IPU_jdim_max,0:IPU_numpeco-1), &
                                    ! holds all subdomains of total field
  rzsubarray (IPU_idim_max,IPU_jdim_max)    ! subdomain for this processor

LOGICAL :: lzacc      ! local logical for gpu execution

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations and Checks
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  ierror   = 0

  ! Initialize GPU flag
  IF (PRESENT(lacc)) THEN
    lzacc = lacc
  ELSE
    lzacc = .FALSE.
  ENDIF

  ! Check the input values
  IF (idim_tot /= IPU_idim_tot .OR. jdim_tot /= IPU_jdim_tot) THEN
    PRINT *, ' WARNING: Wrong input values for gather_field:'
    PRINT *, '      idim_tot = ',idim_tot,' and ie_tot = ',IPU_idim_tot
    PRINT *, '      jdim_tot = ',jdim_tot,' and je_tot = ',IPU_jdim_tot
    ierror = 1
  ENDIF

  IF (idim_loc /= IPU_idim .OR. jdim_loc /= IPU_jdim) THEN
    PRINT *, ' WARNING: Wrong input values for gather_field:'
    PRINT *, '      idim_loc = ',idim_loc,' and ie = ',IPU_idim
    PRINT *, '      jdim_loc = ',jdim_loc,' and je = ',IPU_jdim
    ierror = 1
  ENDIF

  IF (ierror /= 0) RETURN

  !$acc enter data create (rzprocarray,rzsubarray)

!------------------------------------------------------------------------------
! Section 2: Copy subpart into izsubarray
!------------------------------------------------------------------------------

  !$acc kernels present(rfield_loc,rzsubarray) if (lzacc)
  rzsubarray(1:IPU_idim,1:IPU_jdim) = rfield_loc(1:IPU_idim,1:IPU_jdim)
  !$acc end kernels

!------------------------------------------------------------------------------
! Section 3: Gather the subdomains
!------------------------------------------------------------------------------

  IF (IPU_numpeco > 1) THEN
    IF (ireceiver < 0) THEN

      !XL: We need to duplicate the call here because if statement is not
      !    supported with OpenACC host_data directive
      IF (lzacc) THEN
#ifdef NO_MPI_HOST_DATA
        !$acc update host(rzsubarray)
#else
        !$acc host_data use_device(rzsubarray,rzprocarray)
#endif
        CALL MPI_ALLGATHER (rzsubarray,  IPU_idim_max*IPU_jdim_max, IPU_single,&
                            rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_single,&
                            IPU_comm_cart,implcode)
#ifdef NO_MPI_HOST_DATA
        !$acc update device(rzprocarray)
#else
        !$acc end host_data
#endif
      ELSE
        CALL MPI_ALLGATHER (rzsubarray,  IPU_idim_max*IPU_jdim_max, IPU_single,&
                            rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_single,&
                            IPU_comm_cart,implcode)
      ENDIF

    ELSEIF (ireceiver < IPU_numpeco) THEN

      !XL: We need to duplicate the call here because if statement is not
      !    supported with OpenACC host_data directive
      IF (lzacc) THEN
#ifdef NO_MPI_HOST_DATA
        !$acc update host(rzsubarray)
#else
        !$acc host_data use_device(rzsubarray,rzprocarray)
#endif
        CALL MPI_GATHER ( rzsubarray , IPU_idim_max*IPU_jdim_max, IPU_single,  &
                          rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_single,  &
                          ireceiver  , IPU_comm_cart, implcode)
#ifdef NO_MPI_HOST_DATA
        !$acc update device(rzprocarray)
#else
        !$acc end host_data
#endif
      ELSE
        CALL MPI_GATHER ( rzsubarray , IPU_idim_max*IPU_jdim_max, IPU_single,  &
                          rzprocarray, IPU_idim_max*IPU_jdim_max, IPU_single,  &
                          ireceiver  , IPU_comm_cart, implcode)
      ENDIF

      IF (implcode /= 0) THEN
        ierror = 2
        RETURN
      ENDIF
    ELSE
      ierror  = 3
      RETURN
    ENDIF
  ELSE
    !$acc kernels present(rzprocarray,rzsubarray) if(lzacc)
    rzprocarray(:,:,0) = rzsubarray(:,:)
    !$acc end kernels
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Combine the different parts to the total field
!------------------------------------------------------------------------------

  IF (ireceiver < 0) THEN
    CALL combine_subarrays (rzprocarray, rfield_tot, lzacc)
  ELSEIF (ireceiver < IPU_numpeco) THEN
    IF (IPU_cart_id == ireceiver) THEN
      CALL combine_subarrays (rzprocarray, rfield_tot, lzacc)
    ENDIF
  ENDIF

  !$acc exit data delete(rzprocarray,rzsubarray)

#else
  CALL model_abort (-1, -1, 'MPI not available', 'gather_field')
#endif

END SUBROUTINE gather_field_single

!==============================================================================
!==============================================================================
!+ prints LM fields for the total domain
!------------------------------------------------------------------------------

SUBROUTINE print_lm_field                                                &
   (rfield_loc, yname, ydimension, ilev, factor, bias, nout,             &
    istart, iend, jstart, jend, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  This routine gathers the local parts of a total field from all subdomains
!  in the cartesian communicator to processor 0 and prints the part of the
!  total field specified by (istart:iend,jstart:jend) scaled with factor and
!  bias. It has to be called by all processors in this group.
!
!  Caution:
!  The routines specified for the generic routine print_lm_field allocate
!  memory of their own. This has not been avoided because this is regarded
!  to be a routine used only for debug purposes, so that performance and
!  memory usage are not very important.
!
! Method:
!  Call to gather_field
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::     &
  ilev,               & ! level of field_loc in the atmosphere
  nout,               & ! unit number of output file
  istart, iend,       & ! start and end-indices for the part of the
  jstart, jend          ! field that has to be printed

REAL (KIND=wp),           INTENT(IN)       ::     &
  rfield_loc(IPU_idim,IPU_jdim) ! subdomain field

REAL (KIND=wp),           INTENT(IN)       ::     &
  factor, bias          ! for scaling the output according to GRIB

CHARACTER (LEN=*),        INTENT(IN)       ::     &
  yname,              & ! name of the variable
  ydimension            ! dimension (unit) of the variable

! Scalar arguments with intent(out):
INTEGER                 , INTENT(OUT)      ::     &
  ierror                        ! error status variable

CHARACTER (LEN=*)                          ::     &
  yerrmsg               ! for error string

!------------------------------------------------------------------------------

! Local Variables
INTEGER                                    ::                              &
  implcode, n, i1, i2, i, j     ! local error code

REAL (KIND=wp)                             ::                              &
  rzfield_tot(IPU_idim_tot,IPU_jdim_tot)     ! total field

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations and Checks
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  ierror   = 0
  yerrmsg  = '   '

  ! This subroutine has to be called by all routines in icomm
  IF (IPU_numpeco > 1) THEN
    CALL MPI_BARRIER (IPU_comm_cart, implcode)
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Gather the field
!------------------------------------------------------------------------------

  CALL gather_field (rfield_loc,  IPU_idim, IPU_jdim,                       &
                     rzfield_tot, IPU_idim_tot, IPU_jdim_tot, 0, ierror)
  IF (ierror /= 0) RETURN

!------------------------------------------------------------------------------
! Section 3: Print the field
!------------------------------------------------------------------------------

  IF (IPU_cart_id == 0) THEN
    rzfield_tot(:,:) = (rzfield_tot(:,:) + bias) * factor
    i1 = istart
    i2 = MIN(istart+24, iend)
    DO WHILE (i1 <= iend)
      WRITE (nout,'(A,A,A,I3,A,A,A,F8.3,A,F8.1)')                           &
        '  Field: ', yname, '  Level: ', ilev, '  Dimension: ', ydimension, &
        '  Factor: ', 1.0_wp/factor, '  Bias: ', bias
      WRITE (nout,'(A,25I5)') '   I=', (i,i=i1,i2)
      WRITE (nout,'(A)')      '     '
      DO j = jend, jstart, - 1
        WRITE (nout,'(A,I4,25I5)') '  J=',j,(NINT(rzfield_tot(i,j)),i = i1,i2)
      ENDDO
      i1 = i1 + 25
      i2 = MIN (i2 + 25,iend)
    ENDDO
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'print_lm_field')
#endif

END SUBROUTINE print_lm_field

!==============================================================================
!==============================================================================
!+ Computes min, max and meanvalue of a distributed 2D field
!------------------------------------------------------------------------------

SUBROUTINE lm_field_stats (field_loc, ystring, factor, nout, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This routines determines the minimum, the maximum and the meanvalue
!   of a two-dimensional field. For minimum and maximum also the locations
!   are computed. If the field is distributed in the parallel program,
!   every PE computes the values of its subpart which are gathered to all PEs.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                              &
  nout                       ! unit number for the output

REAL    (KIND=wp),        INTENT(IN)       ::                              &
  field_loc(IPU_idim,IPU_jdim), & ! subdomain field
  factor                          ! for scaling the output

INTEGER                 , INTENT(OUT)      ::                              &
  ierror                     ! error status variable

CHARACTER (LEN=*),        INTENT(IN)       ::                              &
  ystring                    ! string for output

CHARACTER (LEN=*),        INTENT(OUT)      ::                              &
  yerrmsg                    ! error message

!------------------------------------------------------------------------------

! Local Variables
INTEGER                                    ::                              &
  ijzmin(2),                   & ! for getting imin and jmin
  ijzmax(2),                   & ! for getting imax and jmax
  izminproc,                   & ! processor that has minimum
  izmaxproc,                   & ! processor that has maximum
  izmplcode, n                   ! local error code

REAL    (KIND=wp)                          ::                              &
  zmin, zmax, zmean          ! min, max and meanvalue of the total 2D field

REAL    (KIND=wp)                          ::                              &
  rzpevalues (7,0:IPU_numpeco-1),&! holds the values of all PEs
  rzvalues   (7)            ! values for this PE

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations and Checks
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  ierror   = 0
  yerrmsg  = '   '
  zmin     = 0.0_wp
  zmax     = 0.0_wp
  zmean    = 0.0_wp

!------------------------------------------------------------------------------
! Section 2: Determine values of this subfield
!------------------------------------------------------------------------------

  zmin     = MINVAL (field_loc(:,:))
  ijzmin   = MINLOC (field_loc(:,:))
  zmax     = MAXVAL (field_loc(:,:))
  ijzmax   = MAXLOC (field_loc(:,:))
  zmean    = SUM    (field_loc(IPU_istart_comp:IPU_iend_comp,      &
                               IPU_jstart_comp:IPU_jend_comp))

!------------------------------------------------------------------------------
! Section 3: Gather the values to PE 0
!------------------------------------------------------------------------------

  IF (IPU_numpeco > 1) THEN
    rzvalues (1) = zmin
    rzvalues (2) = REAL (ijzmin(1), wp)
    rzvalues (3) = REAL (ijzmin(2), wp)
    rzvalues (4) = zmax
    rzvalues (5) = REAL (ijzmax(1), wp)
    rzvalues (6) = REAL (ijzmax(2), wp)
    rzvalues (7) = zmean

    CALL MPI_GATHER (rzvalues, 7, IPU_real, rzpevalues, 7, IPU_real,  &
                     0, IPU_comm_cart, izmplcode)
    IF (izmplcode /= 0) THEN
      ierror  = 1
      yerrmsg = 'Error in MPI_GATHER'
      RETURN
    ENDIF

    IF (IPU_cart_id == 0) THEN
      zmin      = rzpevalues (1,0)
      izminproc = 0
      ijzmin(1) = INT (rzpevalues (2,0))
      ijzmin(2) = INT (rzpevalues (3,0))
      zmax      = rzpevalues (4,0)
      izmaxproc = 0
      ijzmax(1) = INT (rzpevalues (5,0))
      ijzmax(2) = INT (rzpevalues (6,0))
      zmean     = rzpevalues (7,0)

      DO n = 1, IPU_numpeco-1
        zmean = zmean + rzpevalues(7,n)
        IF (rzpevalues(1,n) < zmin) THEN
          zmin      = rzpevalues (1,n)
          izminproc = n
          ijzmin(1) = INT (rzpevalues (2,n))
          ijzmin(2) = INT (rzpevalues (3,n))
        ENDIF
        IF (rzpevalues(4,n) > zmax) THEN
          zmax      = rzpevalues (4,n)
          izmaxproc = n
          ijzmax(1) = INT (rzpevalues (5,n))
          ijzmax(2) = INT (rzpevalues (6,n))
        ENDIF
      ENDDO

      ! Compute grid indices of maximum in the total domain
      ijzmax(1) = IPU_ipositions(izmaxproc,1) - IPU_nbounds - 1 + ijzmax(1)
      ijzmax(2) = IPU_ipositions(izmaxproc,2) - IPU_nbounds - 1 + ijzmax(2)

      ! Compute grid indices of minimum in the total domain
      ijzmin(1) = IPU_ipositions(izminproc,1) - IPU_nbounds - 1 + ijzmin(1)
      ijzmin(2) = IPU_ipositions(izminproc,2) - IPU_nbounds - 1 + ijzmin(2)
    ENDIF
  ENDIF

  ! Compute the meanvalue of the total field
  zmean = zmean / REAL (IPU_idim_tot * IPU_jdim_tot, wp)

!------------------------------------------------------------------------------
! Section 4: Output of variables
!------------------------------------------------------------------------------

  IF (IPU_cart_id == 0) THEN
    WRITE (nout,'(A)') ystring(1:LEN_TRIM(ystring))
    WRITE (nout,'(A,F15.5)')                                                &
           '     Mean    = ', zmean*factor
    WRITE (nout,'(A,F15.5,A,I5,A,I5)')                                      &
           '     Maximum = ', zmax*factor, ' at GP: i = ',ijzmax(1),        &
                                                  ' j = ',ijzmax(2)
    WRITE (nout,'(A,F15.5,A,I5,A,I5)')                                      &
           '     Minimum = ', zmin*factor, ' at GP: i = ',ijzmin(1),        &
                                                  ' j = ',ijzmin(2)
    WRITE (nout,'(A)') '   '
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'lm_field_stats')
#endif

END SUBROUTINE lm_field_stats

!==============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine scatterv_values
!------------------------------------------------------------------------------
!
! SUBROUTINE scatterv_values ( vector_in , idimtot, isender, ivec_len        &
!                            , vector_out, idimloc, my_id  , npes            &
!                            , idatatype, icomm, yerrmsg, ierror )
!
!------------------------------------------------------------------------------
!
! Description:
!  This routine scatters batches of variable length of vector_in from the PE
!  with ID isender individually to the other nodes (so that each node receives
!  one of these batches).
!
! Method:
!  MPI_SCATTERV
!
! Written by:  DWD, Christoph Schraff  (19.04.13)
!==============================================================================

! Following are the different subroutines

!==============================================================================
!+ Subroutine for vector of integers
!------------------------------------------------------------------------------

SUBROUTINE scatterv_integers ( ivector_in , idimtot, isender, ivec_len       &
                             , ivector_out, idimloc, my_id  , npes           &
                             , idatatype, icomm, yerrmsg, ierror )

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                                 &
  idimtot             ,& ! dimension of total vector containing all batches
                         !   (condition: idimtot >= SUM( ivec_len ))
  idimloc             ,& ! dimension of receiving local vector
                         !   (condition: idimloc >= MAXVAL( ivec_len ))
  npes                ,& ! number of PEs
  ivec_len   (npes)   ,& ! length of batches of ivector_in sent to different PEs
  isender             ,& ! PE that scatters (sends) the vector
  my_id               ,& ! local PE in the communicator
  icomm               ,& ! MPI-communicator
  idatatype              ! data type of vector

INTEGER                 , INTENT(IN)       ::                                 &
  ivector_in (idimtot)   ! total vector containing all batches to be scattered

INTEGER                 , INTENT(OUT)      ::                                 &
  ivector_out(idimloc)   ! vector for receiving local batch

CHARACTER (LEN= *)      , INTENT(OUT)      ::                                 &
  yerrmsg                ! string for error messages

INTEGER                 , INTENT(OUT)      ::                                 &
  ierror                 ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                                 &
  ivec_offset(npes)   ,& ! offset of batches of rvector_in sent to different PEs
  nrecv               ,& ! length of batch received by local node
  inode               ,& ! loop index over PEs
  implcode               ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine scatterv_integers
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode    = 0
  ierror      = 0
  yerrmsg     = '   '
  ivector_out = 0

  IF (npes > 1) THEN
    IF ( (0 <= isender) .AND. (isender < npes) ) THEN
      ivec_offset (1) = 0
      DO inode = 2 , npes
        ivec_offset (inode)  =  ivec_offset(inode-1) + ivec_len(inode-1)
      ENDDO
      nrecv = ivec_len(my_id+1)
      CALL MPI_SCATTERV ( ivector_in,  ivec_len, ivec_offset, idatatype       &
                        , ivector_out, nrecv                , idatatype       &
                        , isender, icomm, implcode )
      IF (implcode /= 0) THEN
        ierror  = 1
        yerrmsg = 'Error in MPI_SCATTERV'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid sender'
    ENDIF
  ELSE
    ivector_out(:) = ivector_in(:)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'scatterv_integers')
#endif

END SUBROUTINE scatterv_integers

!==============================================================================
!+ Subroutine for vector of reals
!------------------------------------------------------------------------------

SUBROUTINE scatterv_reals    ( rvector_in , idimtot, isender, ivec_len       &
                             , rvector_out, idimloc, my_id  , npes           &
                             , idatatype, icomm, yerrmsg, ierror )

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                                 &
  idimtot             ,& ! dimension of total vector containing all batches
                         !   (condition: idimtot >= SUM( ivec_len ),
                         !               must be identical on all PEs)
  idimloc             ,& ! dimension of receiving local vector
                         !   (condition: idimloc >= MAXVAL( ivec_len ),
                         !               must be identical on all PEs)
  npes                ,& ! number of PEs
  ivec_len   (npes)   ,& ! length of batches of rvector_in sent to different PEs
  isender             ,& ! PE that scatters (sends) the vector
  my_id               ,& ! local PE in the communicator
  icomm               ,& ! MPI-communicator
  idatatype              ! data type of vector

REAL (KIND=wp),           INTENT(IN)       ::                                 &
  rvector_in (idimtot)   ! total vector containing all batches to be scattered

REAL (KIND=wp),           INTENT(OUT)      ::                                 &
  rvector_out(idimloc)   ! vector for receiving local batch

CHARACTER (LEN= *)      , INTENT(OUT)      ::                                 &
  yerrmsg                ! string for error messages

INTEGER                 , INTENT(OUT)      ::                                 &
  ierror                 ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                                 &
  ivec_offset(npes)   ,& ! offset of batches of rvector_in sent to different PEs
  nrecv               ,& ! length of batch received by local node
  inode               ,& ! loop index over PEs
  implcode               ! error status variable for MPI routines

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine scatterv_reals
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode    = 0
  ierror      = 0
  yerrmsg     = '   '
  rvector_out = 0._wp

  IF (npes > 1) THEN
    IF ( (0 <= isender) .AND. (isender < npes) ) THEN
      ivec_offset (1) = 0
      DO inode = 2 , npes
        ivec_offset (inode)  =  ivec_offset(inode-1) + ivec_len(inode-1)
      ENDDO
      nrecv = ivec_len(my_id+1)
      CALL MPI_SCATTERV ( rvector_in,  ivec_len, ivec_offset, idatatype       &
                        , rvector_out, nrecv                , idatatype       &
                        , isender, icomm, implcode )
      IF (implcode /= 0) THEN
        ierror  = 1
        yerrmsg = 'Error in MPI_SCATTERV'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid sender'
    ENDIF
  ELSE
    rvector_out(:) = rvector_in(:)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'scatterv_reals')
#endif

END SUBROUTINE scatterv_reals

!==============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine gatherv_values
!------------------------------------------------------------------------------
!
! SUBROUTINE gatherv_values ( vector_in , ilenloc, idimloc, my_id  , npes    &
!                           , vector_out, ilentot, idimtot, ireceiver        &
!                           , idatatype, icomm, yerrmsg, ierror )
!
!------------------------------------------------------------------------------
!
! Description:
!  This routine gathers a batch of variable length from each node at the PE
!  with ID ireceiver (or at all PEs, if ireceiver is < 0) and stores them in
!  vector_out (according to the order of PE IDs).
!
! Method:
!  MPI_GATHERV and MPI_ALLGATHERV
!
! Written by:  DWD, Christoph Schraff  (19.04.13)
!==============================================================================

! Following are the different subroutines

!==============================================================================
!+ Subroutine for vector of integers
!------------------------------------------------------------------------------

SUBROUTINE gatherv_integers ( ivector_in , ilenloc, idimloc, npes             &
                            , ivector_out, ilentot, idimtot, ireceiver        &
                            , idatatype, icomm, yerrmsg, ierror, init)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                                 &
  idimloc             ,& ! dimension of sending local vector
                         !   (condition: idimloc >= MAXVAL( ilenloc ),
                         !               must be identical on all PEs)
  idimtot             ,& ! dimension of total vector containing all batches
                         !   (condition: idimtot >= SUM( ilenloc ),
                         !               must be identical on all PEs)
  npes                ,& ! number of PEs
  ilenloc             ,& ! length of local batch sent to receiving PEs
                         !   (may be different on different PEs)
  ireceiver           ,& ! PE that gathers (receives) the vector
  icomm               ,& ! MPI-communicator
  idatatype              ! data type of vector

INTEGER                 , INTENT(IN)       ::                                 &
  ivector_in (idimloc)   ! local vector with local batch to be sent to receiver

LOGICAL                 , INTENT(IN), OPTIONAL ::                             &
  init                   ! whether output array shall be initialized

INTEGER                 , INTENT(OUT)      ::                                 &
  ilentot                ! total number of elements gathered (= SUM( ilenloc ))

INTEGER                 , INTENT(INOUT)    ::                                 &
  ivector_out(idimtot)   ! total vector for receiving all gathered elements

CHARACTER (LEN= *)      , INTENT(OUT)      ::                                 &
  yerrmsg                ! string for error messages

INTEGER                 , INTENT(OUT)      ::                                 &
  ierror                 ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                                 &
  iall_len   (npes)   ,& ! length in rvector_out of batches received from PEs
  iall_offset(npes)   ,& ! offset in rvector_out of batches received from PEs
  inode               ,& ! loop index over PEs
  implcode               ! error status variable for MPI routines

LOGICAL                                    ::                                 &
  ini                    ! whether output array shall be initialized


!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_integers
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode    = 0
  ierror      = 0
  yerrmsg     = '   '
  ilentot     = 0
  ini         = .true.
  if (present(init)) ini = init
  if (ini) ivector_out = 0

  IF (npes > 1) THEN
    IF (ireceiver < npes) THEN
      !   gather the number of elements for each node

      CALL gather_values ( ilenloc, iall_len, 1, npes, IPU_integer, -1        &
                         , icomm, yerrmsg, implcode )
!     ==================

      IF (implcode == 0) THEN
        !   get offset of the first element from each node in 'ivector_out'
        iall_offset (1) = 0
        DO inode = 2 , npes
          iall_offset (inode)  =  iall_offset(inode-1) + iall_len(inode-1)
        ENDDO
        ilentot  =  iall_offset(npes) + iall_len(npes)

        !   perform the gathering
        IF ((ireceiver >= 0) .AND. (ilentot > 0)) THEN

          CALL MPI_GATHERV ( ivector_in , ilenloc              , idatatype    &
                           , ivector_out, iall_len, iall_offset, idatatype    &
                           , ireceiver, icomm, implcode )
!         ================
        ELSEIF (ilentot > 0) THEN

          CALL MPI_ALLGATHERV ( ivector_in , ilenloc              , idatatype &
                              , ivector_out, iall_len, iall_offset, idatatype &
                              , icomm, implcode )
!         ===================
        ENDIF
        IF (implcode /= 0) THEN
          ierror  = 1
          yerrmsg = 'Error in MPI_GATHERV'
        ENDIF
      ELSE
        ierror  = 3
        yerrmsg = 'Error in gatherv_integers'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid receiver'
    ENDIF
  ELSE
    ivector_out(:) = ivector_in(:)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gatherv_integers')
#endif

END SUBROUTINE gatherv_integers

!==============================================================================
!+ Subroutine for 2D-vector of integers
!------------------------------------------------------------------------------

SUBROUTINE gatherv_integers_2D( ivector_in , ilenloc, idimloc, npes           &
                              , ivector_out, ilentot, idimtot, idim1          &
                              , ireceiver  , idatatype, icomm, yerrmsg, ierror&
                              , init)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                                 &
  idimloc             ,& ! length of 2nd dimension of sending local vector
                         !   (condition: idimloc >= MAXVAL( ilenloc ),
                         !               must be identical on all PEs)
  idimtot             ,& ! length of 2nd dimension of total vector containing all batches
                         !   (condition: idimtot >= SUM( ilenloc ),
                         !               must be identical on all PEs)
  npes                ,& ! number of PEs
  ilenloc             ,& ! length of local batch sent to receiving PEs
                         !   (may be different on different PEs)
  ireceiver           ,& ! PE that gathers (receives) the vector
  icomm               ,& ! MPI-communicator
  idatatype           ,& ! data type of vector
  idim1                  ! length of first dimension
                         !    (condition: identical on all PEs and for sending and receiving)

INTEGER                 , INTENT(IN)       ::                                 &
  ivector_in (idim1,idimloc)   ! local vector with local batch to be sent to receiver

LOGICAL                 , INTENT(IN), OPTIONAL ::                             &
  init                   ! whether output array shall be initialized

INTEGER                 , INTENT(OUT)      ::                                 &
  ilentot                ! total number of elements gathered (= SUM( ilenloc ))

INTEGER                 , INTENT(INOUT)    ::                                 &
  ivector_out(idim1,idimtot)   ! total vector for receiving all gathered elements

CHARACTER (LEN= *)      , INTENT(OUT)      ::                                 &
  yerrmsg                ! string for error messages

INTEGER                 , INTENT(OUT)      ::                                 &
  ierror                 ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                                 &
  iall_len   (npes)   ,& ! length in rvector_out of batches received from PEs
  iall_offset(npes)   ,& ! offset in rvector_out of batches received from PEs
  inode               ,& ! loop index over PEs
  implcode               ! error status variable for MPI routines

LOGICAL                                    ::                                 &
  ini                    ! whether output array shall be initialized

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_integers
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode    = 0
  ierror      = 0
  yerrmsg     = '   '
  ilentot     = 0
  ini         = .true.
  if (present(init)) ini = init
  if (ini) ivector_out = 0

  IF (npes > 1) THEN
    IF (ireceiver < npes) THEN
      !   gather the number of elements for each node

      CALL gather_values ( ilenloc, iall_len, 1, npes, IPU_integer, -1        &
                         , icomm, yerrmsg, implcode )
!     ==================

      IF (implcode == 0) THEN
        !   get offset of the first element from each node in 'ivector_out'
        iall_offset (1) = 0
        DO inode = 2 , npes
          iall_offset (inode)  =  iall_offset(inode-1) + iall_len(inode-1)
        ENDDO
        ilentot  =  iall_offset(npes) + iall_len(npes)

        !   perform the gathering
        IF ((ireceiver >= 0) .AND. (ilentot > 0)) THEN

          CALL MPI_GATHERV ( ivector_in , ilenloc*idim1 , idatatype           &
                           , ivector_out, iall_len*idim1, iall_offset*idim1   &
                           , idatatype  , ireceiver    , icomm, implcode )
!         ================
        ELSEIF (ilentot > 0) THEN

          CALL MPI_ALLGATHERV ( ivector_in , ilenloc*idim1 , idatatype        &
                              , ivector_out, iall_len*idim1, iall_offset*idim1&
                              , idatatype  , icomm, implcode )
!         ===================
        ENDIF
        IF (implcode /= 0) THEN
          ierror  = 1
          yerrmsg = 'Error in MPI_GATHERV'
        ENDIF
      ELSE
        ierror  = 3
        yerrmsg = 'Error in gatherv_integers'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid receiver'
    ENDIF
  ELSE
    ivector_out(:,:) = ivector_in(:,:)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gatherv_integers')
#endif

END SUBROUTINE gatherv_integers_2D

!==============================================================================
!+ Subroutine for vector of reals
!------------------------------------------------------------------------------

SUBROUTINE gatherv_reals ( rvector_in , ilenloc, idimloc, npes                &
                         , rvector_out, ilentot, idimtot, ireceiver           &
                         , idatatype, icomm, yerrmsg, ierror, init)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                                 &
  idimloc             ,& ! dimension of sending local vector
                         !   (condition: idimloc >= MAXVAL( ilenloc ),
                         !               must be identical on all PEs)
  idimtot             ,& ! dimension of total vector containing all batches
                         !   (condition: idimtot >= SUM( ilenloc ),
                         !               must be identical on all PEs)
  npes                ,& ! number of PEs
  ilenloc             ,& ! length of local batch sent to receiving PEs
                         !   (may be different on different PEs)
  ireceiver           ,& ! PE that gathers (receives) the vector
  icomm               ,& ! MPI-communicator
  idatatype              ! data type of vector

REAL (KIND=wp),           INTENT(IN)       ::                                 &
  rvector_in (idimloc)   ! local vector with local batch to be sent to receiver

LOGICAL                 , INTENT(IN), OPTIONAL ::                             &
  init                   ! whether output array shall be initialized

INTEGER                 , INTENT(OUT)      ::                                 &
  ilentot                ! total number of elements gathered (= SUM( ilenloc ))

REAL (KIND=wp),           INTENT(INOUT)    ::                                 &
  rvector_out(idimtot)   ! total vector for receiving all gathered elements

CHARACTER (LEN= *)      , INTENT(OUT)      ::                                 &
  yerrmsg                ! string for error messages

INTEGER                 , INTENT(OUT)      ::                                 &
  ierror                 ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                                 &
  iall_len   (npes)   ,& ! length in rvector_out of batches received from PEs
  iall_offset(npes)   ,& ! offset in rvector_out of batches received from PEs
  inode               ,& ! loop index over PEs
  implcode               ! error status variable for MPI routines

LOGICAL                                    ::                                 &
  ini                    ! whether output array shall be initialized

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_reals
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode    = 0
  ierror      = 0
  yerrmsg     = '   '
  ilentot     = 0
  ini         = .true.
  if (present(init)) ini = init
  if (ini) rvector_out = 0

  IF (npes > 1) THEN
    IF (ireceiver < npes) THEN
      !   gather the number of elements for each node

      CALL gather_values ( ilenloc, iall_len, 1, npes, IPU_integer, -1        &
                         , icomm, yerrmsg, implcode )
!     ==================

      IF (implcode == 0) THEN
        !   get offset of the first element from each node in 'rvector_out'
        iall_offset (1) = 0
        DO inode = 2 , npes
          iall_offset (inode)  =  iall_offset(inode-1) + iall_len(inode-1)
        ENDDO
        ilentot  =  iall_offset(npes) + iall_len(npes)

        !   perform the gathering
        IF ((ireceiver >= 0) .AND. (ilentot > 0)) THEN

          CALL MPI_GATHERV ( rvector_in , ilenloc              , idatatype    &
                           , rvector_out, iall_len, iall_offset, idatatype    &
                           , ireceiver, icomm, implcode )
!         ================
        ELSEIF (ilentot > 0) THEN

          CALL MPI_ALLGATHERV ( rvector_in , ilenloc              , idatatype &
                              , rvector_out, iall_len, iall_offset, idatatype &
                              , icomm, implcode )
!         ===================
        ENDIF
        IF (implcode /= 0) THEN
          ierror  = 1
          yerrmsg = 'Error in MPI_GATHERV'
        ENDIF
      ELSE
        ierror  = 3
        yerrmsg = 'Error in gatherv_reals'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid receiver'
    ENDIF
  ELSE
    rvector_out(:) = rvector_in(:)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gatherv_reals')
#endif

END SUBROUTINE gatherv_reals


!==============================================================================
!+ Subroutine for vector of reals
!------------------------------------------------------------------------------

SUBROUTINE gatherv_reals_2D ( rvector_in , ilenloc, idimloc, npes         &
                            , rvector_out, ilentot, idimtot, idim1        &
                            , ireceiver, idatatype, icomm, yerrmsg, ierror&
                            , init)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                                 &
  idimloc             ,& ! dimension of sending local vector
                         !   (condition: idimloc >= MAXVAL( ilenloc ),
                         !               must be identical on all PEs)
  idimtot             ,& ! dimension of total vector containing all batches
                         !   (condition: idimtot >= SUM( ilenloc ),
                         !               must be identical on all PEs)
  npes                ,& ! number of PEs
  ilenloc             ,& ! length of local batch sent to receiving PEs
                         !   (may be different on different PEs)
  ireceiver           ,& ! PE that gathers (receives) the vector
  icomm               ,& ! MPI-communicator
  idatatype           ,& ! data type of vector
  idim1                  ! length of first dimension
                         !    (condition: identical on all PEs and for sending and receiving)

REAL (KIND=wp),           INTENT(IN)       ::                                 &
  rvector_in (idim1,idimloc)   ! local vector with local batch to be sent to receiver

LOGICAL                 , INTENT(IN), OPTIONAL ::                             &
  init                   ! whether output array shall be initialized

INTEGER                 , INTENT(OUT)      ::                                 &
  ilentot                ! total number of elements gathered (= SUM( ilenloc ))

REAL (KIND=wp),           INTENT(INOUT)    ::                                 &
  rvector_out(idim1,idimtot)   ! total vector for receiving all gathered elements

CHARACTER (LEN= *)      , INTENT(OUT)      ::                                 &
  yerrmsg                ! string for error messages

INTEGER                 , INTENT(OUT)      ::                                 &
  ierror                 ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                                 &
  iall_len   (npes)   ,& ! length in rvector_out of batches received from PEs
  iall_offset(npes)   ,& ! offset in rvector_out of batches received from PEs
  inode               ,& ! loop index over PEs
  implcode               ! error status variable for MPI routines

LOGICAL                                    ::                                 &
  ini                    ! whether output array shall be initialized

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_reals
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode    = 0
  ierror      = 0
  yerrmsg     = '   '
  ilentot     = 0
  ini         = .true.
  if (present(init)) ini = init
  if (ini) rvector_out = 0

  IF (npes > 1) THEN
    IF (ireceiver < npes) THEN
      !   gather the number of elements for each node

      CALL gather_values ( ilenloc, iall_len, 1, npes, IPU_integer, -1        &
                         , icomm, yerrmsg, implcode )
!     ==================

      IF (implcode == 0) THEN
        !   get offset of the first element from each node in 'rvector_out'
        iall_offset (1) = 0
        DO inode = 2 , npes
          iall_offset (inode)  =  iall_offset(inode-1) + iall_len(inode-1)
        ENDDO
        ilentot  =  iall_offset(npes) + iall_len(npes)

        !   perform the gathering
        IF ((ireceiver >= 0) .AND. (ilentot > 0)) THEN

          CALL MPI_GATHERV ( rvector_in , ilenloc*idim1 , idatatype           &
                           , rvector_out, iall_len*idim1, iall_offset*idim1   &
                           , idatatype  , ireceiver    , icomm, implcode )
!         ================
        ELSEIF (ilentot > 0) THEN

          CALL MPI_ALLGATHERV ( rvector_in , ilenloc*idim1 , idatatype         &
                              , rvector_out, iall_len*idim1, iall_offset*idim1 &
                              , idatatype  , icomm, implcode )
!         ===================
        ENDIF
        IF (implcode /= 0) THEN
          ierror  = 1
          yerrmsg = 'Error in MPI_GATHERV'
        ENDIF
      ELSE
        ierror  = 3
        yerrmsg = 'Error in gatherv_reals'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid receiver'
    ENDIF
  ELSE
    rvector_out(:,:) = rvector_in(:,:)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gatherv_reals')
#endif

END SUBROUTINE gatherv_reals_2D


!==============================================================================
!+ Subroutine for vector of reals
!------------------------------------------------------------------------------

SUBROUTINE gatherv_reals_3D ( rvector_in , ilenloc, idimloc, npes         &
                            , rvector_out, ilentot, idimtot, idim1, idim2 &
                            , ireceiver, idatatype, icomm, yerrmsg, ierror&
                            , init)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                                 &
  idimloc             ,& ! dimension of sending local vector
                         !   (condition: idimloc >= MAXVAL( ilenloc ),
                         !               must be identical on all PEs)
  idimtot             ,& ! dimension of total vector containing all batches
                         !   (condition: idimtot >= SUM( ilenloc ),
                         !               must be identical on all PEs)
  npes                ,& ! number of PEs
  ilenloc             ,& ! length of local batch sent to receiving PEs
                         !   (may be different on different PEs)
  ireceiver           ,& ! PE that gathers (receives) the vector
  icomm               ,& ! MPI-communicator
  idatatype           ,& ! data type of vector
  idim1               ,& ! length of first dimension
  idim2                  ! length of first dimension
                         !    (condition: identical on all PEs and for sending and receiving)

REAL (KIND=wp),           INTENT(IN)       ::                                 &
  rvector_in (idim1,idim2,idimloc)   ! local vector with local batch to be sent to receiver

LOGICAL                 , INTENT(IN), OPTIONAL ::                             &
  init                   ! whether output array shall be initialized

INTEGER                 , INTENT(OUT)      ::                                 &
  ilentot                ! total number of elements gathered (= SUM( ilenloc ))

REAL (KIND=wp),           INTENT(INOUT)    ::                                 &
  rvector_out(idim1,idim2,idimtot)   ! total vector for receiving all gathered elements

CHARACTER (LEN= *)      , INTENT(OUT)      ::                                 &
  yerrmsg                ! string for error messages

INTEGER                 , INTENT(OUT)      ::                                 &
  ierror                 ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                                 &
  iall_len   (npes)   ,& ! length in rvector_out of batches received from PEs
  iall_offset(npes)   ,& ! offset in rvector_out of batches received from PEs
  inode               ,& ! loop index over PEs
  implcode            ,& ! error status variable for MPI routines
  idim                   !

LOGICAL                                    ::                                 &
  ini                    ! whether output array shall be initialized

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_reals
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode    = 0
  ierror      = 0
  yerrmsg     = '   '
  ilentot     = 0
  ini         = .true.
  if (present(init)) ini = init
  if (ini) rvector_out = 0
  idim        = idim1*idim2

  IF (npes > 1) THEN
    IF (ireceiver < npes) THEN
      !   gather the number of elements for each node

      CALL gather_values ( ilenloc, iall_len, 1, npes, IPU_integer, -1        &
                         , icomm, yerrmsg, implcode )
!     ==================

      IF (implcode == 0) THEN
        !   get offset of the first element from each node in 'rvector_out'
        iall_offset (1) = 0
        DO inode = 2 , npes
          iall_offset (inode)  =  iall_offset(inode-1) + iall_len(inode-1)
        ENDDO
        ilentot  =  iall_offset(npes) + iall_len(npes)

        !   perform the gathering
        IF ((ireceiver >= 0) .AND. (ilentot > 0)) THEN

          CALL MPI_GATHERV ( rvector_in , ilenloc*idim  , idatatype           &
                           , rvector_out, iall_len*idim , iall_offset*idim    &
                           , idatatype  , ireceiver    , icomm, implcode )
!         ================
        ELSEIF (ilentot > 0) THEN

          CALL MPI_ALLGATHERV ( rvector_in , ilenloc*idim  , idatatype         &
                              , rvector_out, iall_len*idim , iall_offset*idim  &
                              , idatatype  , icomm, implcode )
!         ===================
        ENDIF
        IF (implcode /= 0) THEN
          ierror  = 1
          yerrmsg = 'Error in MPI_GATHERV'
        ENDIF
      ELSE
        ierror  = 3
        yerrmsg = 'Error in gatherv_reals'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid receiver'
    ENDIF
  ELSE
    rvector_out(:,:,:) = rvector_in(:,:,:)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gatherv_reals')
#endif

END SUBROUTINE gatherv_reals_3D
!==============================================================================
!+ Subroutine for 2D-vector of integers
!------------------------------------------------------------------------------

SUBROUTINE gatherv_logical_2D(  lvector_in , ilenloc, idimloc, npes           &
                              , lvector_out, ilentot, idimtot, idim1          &
                              , ireceiver  , idatatype, icomm, yerrmsg, ierror&
                              , init)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER                 , INTENT(IN)       ::                                 &
  idimloc             ,& ! length of 2nd dimension of sending local vector
                         !   (condition: idimloc >= MAXVAL( ilenloc ),
                         !               must be identical on all PEs)
  idimtot             ,& ! length of 2nd dimension of total vector containing all batches
                         !   (condition: idimtot >= SUM( ilenloc ),
                         !               must be identical on all PEs)
  npes                ,& ! number of PEs
  ilenloc             ,& ! length of local batch sent to receiving PEs
                         !   (may be different on different PEs)
  ireceiver           ,& ! PE that gathers (receives) the vector
  icomm               ,& ! MPI-communicator
  idatatype           ,& ! data type of vector
  idim1                  ! length of first dimension
                         !    (condition: identical on all PEs and for sending and receiving)

LOGICAL                 , INTENT(IN)       ::                                 &
  lvector_in (idim1,idimloc)   ! local vector with local batch to be sent to receiver

LOGICAL                 , INTENT(IN), OPTIONAL ::                             &
  init                   ! whether output array shall be initialized

INTEGER                 , INTENT(OUT)      ::                                 &
  ilentot                ! total number of elements gathered (= SUM( ilenloc ))

LOGICAL                 , INTENT(INOUT)    ::                                 &
  lvector_out(idim1,idimtot)   ! total vector for receiving all gathered elements

CHARACTER (LEN= *)      , INTENT(OUT)      ::                                 &
  yerrmsg                ! string for error messages

INTEGER                 , INTENT(OUT)      ::                                 &
  ierror                 ! error status variable

!------------------------------------------------------------------------------

! Local variables
INTEGER                                    ::                                 &
  iall_len   (npes)   ,& ! length in rvector_out of batches received from PEs
  iall_offset(npes)   ,& ! offset in rvector_out of batches received from PEs
  inode               ,& ! loop index over PEs
  implcode               ! error status variable for MPI routines

LOGICAL                                    ::                                 &
  ini                    ! whether output array shall be initialized

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine gatherv_integers
!------------------------------------------------------------------------------

#ifndef NOMPI
  ! Initializations
  implcode    = 0
  ierror      = 0
  yerrmsg     = '   '
  ilentot     = 0
  ini         = .true.
  if (present(init)) ini = init
  if (ini) lvector_out = .FALSE.

  IF (npes > 1) THEN
    IF (ireceiver < npes) THEN
      !   gather the number of elements for each node

      CALL gather_values ( ilenloc, iall_len, 1, npes, IPU_integer, -1        &
                         , icomm, yerrmsg, implcode )
!     ==================

      IF (implcode == 0) THEN
        !   get offset of the first element from each node in 'lvector_out'
        iall_offset (1) = 0
        DO inode = 2 , npes
          iall_offset (inode)  =  iall_offset(inode-1) + iall_len(inode-1)
        ENDDO
        ilentot  =  iall_offset(npes) + iall_len(npes)

        !   perform the gathering
        IF ((ireceiver >= 0) .AND. (ilentot > 0)) THEN

          CALL MPI_GATHERV ( lvector_in , ilenloc*idim1 , idatatype           &
                           , lvector_out, iall_len*idim1, iall_offset*idim1   &
                           , idatatype  , ireceiver    , icomm, implcode )
!         ================
        ELSEIF (ilentot > 0) THEN

          CALL MPI_ALLGATHERV ( lvector_in , ilenloc*idim1 , idatatype        &
                              , lvector_out, iall_len*idim1, iall_offset*idim1&
                              , idatatype  , icomm, implcode )
!         ===================
        ENDIF
        IF (implcode /= 0) THEN
          ierror  = 1
          yerrmsg = 'Error in MPI_GATHERV'
        ENDIF
      ELSE
        ierror  = 3
        yerrmsg = 'Error in gatherv_integers'
      ENDIF
    ELSE
      ierror  = 2
      yerrmsg = 'no valid receiver'
    ENDIF
  ELSE
    lvector_out(:,:) = lvector_in(:,:)
  ENDIF
#else
  CALL model_abort (-1, -1, 'MPI not available', 'gatherv_integers')
#endif

END SUBROUTINE gatherv_logical_2D

!==============================================================================
!==============================================================================
!+ Gathering 2-dim. arrays from all nodes, giving the collection to all nodes
!------------------------------------------------------------------------------

SUBROUTINE gather_all ( nlenint, nlenreal, nlenchar, maxcol, istrlen         &
                      , ivect  , rvect   , yvect   , ncol  , nsumcol         &
                      , yerrmsg, ierror )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine gathers (collects) from all nodes the first 'ncol' columns
!   of the local 2-dimensional arrays 'ivect', 'rvect', resp. 'yvect' of types
!   integer, real, resp. character. The length of the columns is 'nlenint',
!   'nlenreal', resp. 'nlenchar'. The resulting (collected) arrays are written
!   back to the arrays 'ivect', 'rvect', resp. 'yvect', and number of collected
!   columns are given by (the updated) 'ncol' (and by 'nsumcol').
!   The order of the placement of the local columns in the resulting arrays is
!   according to the rank of the node from which the columns have been gathered.
!   If the sum of the local number of columns 'nsumcol' exceeds the specified
!   maximum number of columns 'maxcol' for the receiving arrays, only the first
!   'maxcol' columns are written to these arrays, and the output value of 'ncol'
!   equals 'maxcol'.
!   (Note: If 'nlenint', 'nlenreal', or 'nlenchar' are specified to be zero,
!          the first dimension of arrays 'ivect', 'rvect', resp. 'yvect' has
!          nevertheless length 1 rather than zero.)
!
! Method:
!   Use of MPI-routines MPI_ALLGATHER to gather the local numbers of columns and
!   of MPI_ALLGATHERV to gather the local columns.
!   Characters are converted to integers, gathered as integers, and reconverted
!   to characters.
!
!------------------------------------------------------------------------------

! Local parameters:
! ----------------

  INTEGER                 , PARAMETER           ::       &
    i1  =  1

! Subroutine arguments:
! --------------------

  INTEGER                 , INTENT (IN)         ::       &
    nlenint                     ,& ! number of elements in each integer   column
    nlenreal                    ,& ! number of elements in each real      column
    nlenchar                    ,& ! number of elements in each character column
    maxcol                      ,& ! maximum number of columns
    istrlen                        ! string length of each character element

  INTEGER                 , INTENT (INOUT)      ::       &
    ncol                           ! number of    columns   to be /been gathered

  INTEGER                 , INTENT (INOUT)      ::       &
    ivect(MAX(nlenint,i1) ,maxcol) ! contains the int. data to be /been gathered

  REAL    (KIND=wp)       , INTENT (INOUT)      ::       &
    rvect(MAX(nlenreal,i1),maxcol) ! contains the real data to be /been gathered

  CHARACTER (LEN=istrlen) , INTENT (INOUT)      ::       &
    yvect(MAX(nlenchar,i1),maxcol) ! contains the char data to be /been gathered

  INTEGER                 , INTENT (OUT)        ::       &
    nsumcol                        ! sum of local number of columns

  CHARACTER (LEN= *)      , INTENT (OUT)        ::       &
    yerrmsg

  INTEGER                 , INTENT (OUT)        ::       &
    ierror

! Local scalars:
! -------------

  INTEGER                    ::       &
    izmplcode         ,& ! for MPI error code
    icomm             ,& ! MPI communicator
    my_id             ,& ! ID of current node in 'icomm'
    iscount           ,& ! number of columns at local node
    iscdim , ircdim   ,& ! dimensions of sending / receiving counters of columns
    iisdim            ,& ! number of integer   elements sent to each node
    irsdim            ,& ! number of real      elements sent to each node
    iysdim            ,& ! number of character elements sent to each node
    ncolloc           ,& ! number of columns at local node only
    icol , irow , icl ,& ! loop indices for conversion CHARACTER <--> INTEGER
    icc , nzerr

  INTEGER                    ::       &
    imp_char_gather      ! determines the correct CHARACTER type used in the
                         ! model for MPI

  CHARACTER (LEN=25) yroutine

! Local (automatic) arrays:
! -------------------------

  INTEGER                    ::       &
    ircount  (IPU_numpeco)  ,& ! number of columns
    iirdim   (IPU_numpeco)  ,& ! number of integer   elements received from each node
    irrdim   (IPU_numpeco)  ,& ! number of real      elements received from each node
    iyrdim   (IPU_numpeco)  ,& ! number of character elements received from each node
    iidispl  (IPU_numpeco)  ,& ! displacement of incoming integer data within 'ivect'
    irdispl  (IPU_numpeco)  ,& ! displacement of incoming real    data within 'rvect'
    iydispl  (IPU_numpeco)     ! displacement of incoming charac. data within 'yvect'

  INTEGER                    ::       &
    isbuf    (MAX(nlenint ,i1)          ,maxcol)    ! for intermediate storage

  REAL    (KIND=wp)          ::       &
    rsbuf    (MAX(nlenreal,i1)          ,maxcol)    ! for intermediate storage

! Normal MPI version
! CHARACTER (LEN=istrlen)    ::       &
!   ysbuf    (MAX(nlenchar,i1)          ,maxcol) ,& ! for intermediate storage
!   yrbuf    (MAX(nlenchar,i1)          ,maxcol)    ! for intermediate storage
! T3E version (knows no CHAR BCAST)
  INTEGER                    ::       &
    ysbuf    (MAX(nlenchar,i1), istrlen ,maxcol) ,& ! for intermediate storage
    yrbuf    (MAX(nlenchar,i1), istrlen ,maxcol)    ! for intermediate storage

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE gather_all
!------------------------------------------------------------------------------

  ierror  = 0
  yerrmsg = '   '

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

#ifndef NOMPI

! Normal MPI version
  ! Determine the type CHARACTERs for MPI and other variables.
  ! If the KIND-type parameters in data_parameters are changed, the
  ! variables here have to be changed accordingly.
  imp_char_gather = IPU_integer

  yroutine  = 'gather_all'
  icomm     =  IPU_comm_cart
  my_id     =  IPU_cart_id

  ncolloc = ncol
  IF (ncolloc > 0) THEN
    IF (nlenint  > 0) isbuf (1:nlenint ,1:ncolloc) = ivect(1:nlenint ,1:ncolloc)
    IF (nlenreal > 0) rsbuf (1:nlenreal,1:ncolloc) = rvect(1:nlenreal,1:ncolloc)
! Normal MPI version
!   IF (nlenchar > 0) ysbuf (1:nlenchar,1:ncolloc) = yvect(1:nlenchar,1:ncolloc)
! T3E version (knows no CHAR BCAST)
    IF ((nlenchar > 0) .AND. (istrlen > 0)) THEN
      DO    icol = 1 , ncolloc
!NEC$ noouterloop_unroll
        DO  irow = 1 , nlenchar
          DO icl = 1 , istrlen
            ysbuf (irow,icl,icol) = ICHAR( yvect(irow,icol) (icl:icl) )
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- Section 2: Gathering the number of columns from & for all nodes
!------------------------------------------------------------------------------

  iscount = ncolloc
  iscdim  = 1
  ircdim  = 1

  CALL MPI_ALLGATHER ( iscount, iscdim, IPU_integer                           &
                     , ircount, ircdim, IPU_integer, icomm, izmplcode )
! ==================

  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_ALLGATHER failed in gather_all'
  ENDIF

!------------------------------------------------------------------------------
!- Section 3: Gathering the columns from & for all nodes
!------------------------------------------------------------------------------

  ncol    = 0
  nsumcol = 0
  DO icc = 1, IPU_numpeco
    iidispl (icc) = nlenint  * ncol
    irdispl (icc) = nlenreal * ncol
    iydispl (icc) = nlenchar * ncol * istrlen
    nsumcol       = nsumcol  + ircount(icc)
! avoid writing beyond array bounds of receiving (global) arrays
    IF   (nsumcol > maxcol)    ircount(icc) = maxcol - ncol
    ncol          = ncol     + ircount(icc)
    iirdim  (icc) = nlenint  * ircount(icc)
    irrdim  (icc) = nlenreal * ircount(icc)
    iyrdim  (icc) = nlenchar * ircount(icc) * istrlen
  ENDDO
  ncolloc = ircount(my_id+1)
  iisdim  = nlenint  * ncolloc
  irsdim  = nlenreal * ncolloc
  iysdim  = nlenchar * ncolloc * istrlen

  IF (ncol > 0) THEN
    IF (nlenint > 0) THEN

      CALL MPI_ALLGATHERV ( isbuf, iisdim         , IPU_integer               &
                          , ivect, iirdim, iidispl, IPU_integer               &
                          , icomm, izmplcode )
!     ===================

      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_ALLGATHERV failed in gather_all'
      ENDIF
    ENDIF

    IF (nlenreal > 0) THEN

      CALL MPI_ALLGATHERV ( rsbuf, irsdim         , IPU_real                  &
                          , rvect, irrdim, irdispl, IPU_real                  &
                          , icomm, izmplcode )
!     ===================

      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_ALLGATHERV failed in gather_all'
      ENDIF
    ENDIF

    IF ((nlenchar > 0) .AND. (istrlen > 0)) THEN

      CALL MPI_ALLGATHERV ( ysbuf, iysdim         , imp_char_gather            &
                          , yrbuf, iyrdim, iydispl, imp_char_gather            &
                          , icomm, izmplcode )
!     ===================

      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_ALLGATHERV failed in gather_all'
      ENDIF
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- Section 4: Copy INTEGER buffers to CHARACTER and LOGICAL result fields
!------------------------------------------------------------------------------

  IF (ncol > 0) THEN
! Normal MPI version
!   IF (nlenchar > 0) yvect (1:nlenchar,1:ncol) = yrbuf (1:nlenchar,1:ncol)
! T3E version (knows no CHAR BCAST)
    IF ((nlenchar > 0) .AND. (istrlen > 0)) THEN
      DO    icol = 1 , ncol
        DO  irow = 1 , nlenchar
          DO icl = 1 , istrlen
            yvect (irow,icol) (icl:icl) = CHAR( yrbuf(irow,icl,icol) )
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

#else
  CALL model_abort (-1, -1, 'MPI not available', 'gather_all')
#endif

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE gather_all

!===============================================================================
!===============================================================================
!+ Individual (amount of) information sent and received by each individual node
!-------------------------------------------------------------------------------

SUBROUTINE reports_all2all ( iscount , ircount , npe                           &
                           , maxsrep, isendbuf, rsendbuf, nelemi , nelemr      &
                           , maxrrep, irecvbuf, rrecvbuf, yerrmsg, ierror )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine exchanges 'reports' between nodes (processors) such that each
!   node sends and receives other reports compared to other nodes. The number
!   of reports may also be different for each pair of communicating nodes and
!   for the two directions of sending.
!   A 'report' consists of a fixed number of integer elements and another
!   fixed number of real elements.
!
! Method:
!   Use of MPI-routine MPI_ALLTOALLV.
!   The sending and receiving buffers are 2-dimensional arrays.
!   The reports to be sent and amount of reports to be sent and received by
!   each individual node are input variables and have to be specified prior to
!   the call of this routine.
!
! Written by:  DWD, Christoph Schraff  (18.04.13)
!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                 , INTENT (IN)         ::       &
    npe             ,& ! number of nodes
    iscount  (npe)  ,& ! number of reports sent  by local node to other nodes
    ircount  (npe)  ,& ! number reports received by local node from other nodes
    nelemi          ,& ! number of integer elements in each report
    nelemr          ,& ! number of real    elements in each report
    maxsrep         ,& ! maximum number of reports in  sending  buffers
    maxrrep            ! maximum number of reports in receiving buffers

  REAL    (KIND=wp)       , INTENT (IN)         ::       &
    rsendbuf (nelemr, maxsrep)     ! contains the real data sent by local node

  INTEGER                 , INTENT (IN)         ::       &
    isendbuf (nelemi, maxsrep)     ! contains integer data sent by local node

  REAL    (KIND=wp)       , INTENT (OUT)        ::       &
    rrecvbuf (nelemr, maxrrep)     ! real data to be received by local node

  INTEGER                 , INTENT (OUT)        ::       &
    irecvbuf (nelemi, maxrrep)     ! integer data to be received by local node

  CHARACTER (LEN= *)      , INTENT (OUT)        ::       &
    yerrmsg            ! error message (LEN must be >= 40)

  INTEGER                 , INTENT (OUT)        ::       &
    ierror             ! error return code

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER                    ::       &
    izmplcode       ,& ! for MPI error code
    icomm           ,& ! MPI communicator
    my_id           ,& ! ID of current node in 'icomm'
    inode           ,& ! node number
    jerr , ierr        ! error codes

  CHARACTER (LEN=25) yroutine
  CHARACTER (LEN=50) yerr      ! error message

! Local (automatic) arrays:
! -------------------------

  INTEGER                    ::       &
    isrcount (npe)  ,& ! number       of real    elements   sent   by local node
    isicount (npe)  ,& ! number       of integer elements   sent   by local node
    irrcount (npe)  ,& ! number       of real    elements received by local node
    iricount (npe)  ,& ! number       of integer elements received by local node
    isrdispl (npe)  ,& ! displacement of real    elements   sent   by local node
    isidispl (npe)  ,& ! displacement of integer elements   sent   by local node
    irrdispl (npe)  ,& ! displacement of real    elements received by local node
    iridispl (npe)     ! displacement of integer elements received by local node
!
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE reports_all2all
!-------------------------------------------------------------------------------

  ierror   = 0
  yerrmsg  = '   '
  irecvbuf = 0
  rrecvbuf = 0._wp

!-------------------------------------------------------------------------------
!- Section 1: Initializations
!-------------------------------------------------------------------------------

#ifndef NOMPI

! Normal MPI version
  yroutine  = 'reports_all2all'
  icomm     =  IPU_comm_cart
  my_id     =  IPU_cart_id

! Scatter and gather the multi-level reports
! ------------------------------------------

  DO inode = 1, npe
    IF (inode >= 2) THEN
      isrdispl (inode)  =  isrdispl(inode-1) + isrcount(inode-1)
      isidispl (inode)  =  isidispl(inode-1) + isicount(inode-1)
      irrdispl (inode)  =  irrdispl(inode-1) + irrcount(inode-1)
      iridispl (inode)  =  iridispl(inode-1) + iricount(inode-1)
    ELSE
      isrdispl (inode)  =  0
      isidispl (inode)  =  0
      irrdispl (inode)  =  0
      iridispl (inode)  =  0
    ENDIF
    isrcount (inode)  =  iscount(inode) * nelemr
    isicount (inode)  =  iscount(inode) * nelemi
    irrcount (inode)  =  ircount(inode) * nelemr
    iricount (inode)  =  ircount(inode) * nelemi
  ENDDO

  CALL MPI_ALLTOALLV ( rsendbuf , isrcount , isrdispl , IPU_real               &
                     , rrecvbuf , irrcount , irrdispl , IPU_real               &
                     , icomm    , izmplcode )
! ==================

  jerr = ABS( izmplcode )

  CALL global_values ( jerr, 1, 'MAX', IPU_integer, icomm, -1, yerr, ierr )
! ==================
  IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
    WRITE( yerrmsg,'("* MPI_ALLTOALLV REAL *",3I6)' ) jerr, ierr, izmplcode
  ELSE

    CALL MPI_ALLTOALLV ( isendbuf , isicount , isidispl , IPU_integer          &
                       , irecvbuf , iricount , iridispl , IPU_integer          &
                       , icomm    , izmplcode )
  ! ==================

    jerr = ABS( izmplcode )

    CALL global_values( jerr, 1, 'MAX', IPU_integer, icomm, -1, yerr, ierr )
  ! ==================
    IF ((jerr /= 0) .OR. (ierr /= 0))                                          &
      WRITE( yerrmsg,'("* MPI_ALLTOALLV INT *",3I6)' ) jerr, ierr, izmplcode
  ENDIF
  ierror  =  izmplcode

#else
  CALL model_abort (-1, -1, 'MPI not available', 'reports_all2all')
#endif

!-------------------------------------------------------------------------------
!- End of the Subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE reports_all2all

!==============================================================================
!==============================================================================

SUBROUTINE exchange_profiles (nlev, nprof, i_list, j_list, node_list,    &
                              field_3d, ktop, kbot, profiles)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine of the module 'parallel_utilities' exchanges vertical
!   profiles (columns) of the model field 'field_3d' between different
!   processors. Any number of profiles from other nodes may be requested,
!   they can be located arbitrarily at the different nodes.
!   The number of profiles 'nprof' and lists of their local coordinates
!   at the nodes where they are to be found and fetched from have to be
!   specified in the lists 'i_list,j_list,node_list'.
!   The exchange may be limited to the levels between ktop and kbot to
!   restrict the necessary data exchange/communication.
!
! Method:
!   The 'order lists' of ij-coordinates requested by each processor are
!   exchanged between all processors using MPI_ALLTOALLV.
!   The data are sent and received only between those processors where a data
!   exchange is necessary (i.e. where the order list has not been empty).
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
!------------------------------------------------------------------------------
! Scalar arguments with intent(in):
!---------------------------------

  INTEGER                 , INTENT(IN)             ::       &
    nprof             ,& ! number of profiles requested from other nodes
    kbot              ,& ! index for lowest model layer in transfered profile
    ktop                 ! index for uppest model layer in transfered profile


! Array arguments with intent(in):
!---------------------------------
  INTEGER                 , INTENT(IN)             ::       &
    i_list(nprof)     ,& ! list of i indeces of requested profiles
                         ! (coord. at node concerned)
    j_list(nprof)     ,& ! list of j indeces of requested profiles
                         ! (coord. at node concerned)
    node_list(nprof)     ! list of nodes at which the profiles can be found

  REAL (KIND=wp),     INTENT(IN)                   ::       &
    field_3d(IPU_idim,IPU_jdim,IPU_kdim)   ! model field of which profiles
                                           ! are to be exchanged

! Array arguments with intent(out):
!---------------------------------
  REAL (KIND=wp),     INTENT(OUT)                  ::       &
    profiles(ktop:kbot,nprof)   ! fetched profiles of the requested model field

! Local parameters, scalars, arrays :
!------------------------------------------------------------------------------
! Local scalars:
!---------------

  INTEGER                                ::       &
    nlev             ,& ! number of levels in requested profile transfer
    nfetch           ,& ! number of profiles to be fetched from other nodes
    len,ip           ,& ! loop indeces
    iproc            ,& ! loop index over the nodes/processors
    inode            ,& ! node number
    istat,izlocstat  ,& ! error flags for data array allocation
    implcode         ,& ! error code due to MPI communication
    ind,i,j,k,index     ! indeces

! Local arrays:
!---------------

! arrays for distribution of 'order lists' for the profiles
   INTEGER                               ::  &
     ij_order(2*nprof)        ,& ! buffer array for sending index lists
     num_from_node(nprof)     ,& ! place of each profile in the list from a PE
     icount_order(IPU_numpeco),& ! array for lengths of sending buffer sections
     icount_distr(IPU_numpeco),& ! array for lengths of recving buffer sections
     idispl_order(IPU_numpeco),& ! array for starts of sending buffer sections
     idispl_distr(IPU_numpeco)   ! array for receiving recving buffer sections

   INTEGER                 , ALLOCATABLE   ::  &
     ij_distr(:)                 ! buffer for receiving index lists of all PEs

! arrays for exchange of the profile data
   REAL  (KIND=wp),     ALLOCATABLE                 ::  &
     rrecvbuf(:)                   ,& ! buffer array for receiving profiles
     rsendbuf(:)                      ! buffer array for sending profiles

!------------------------------------------------------------------------------

#ifndef NOMPI
!------------------------------------------------------------------------------
! Section 1 : Exchange of the 'order lists' of profiles needed from other nodes
!             Allocate space for send and receive buffers
!------------------------------------------------------------------------------

! Note OpenACC: When running on accelerator, section 1 is computed
! on the CPU. This assumes  i_list, j_list and node_list are up to date on the CPU

! sort info on needed profiles according to node-number into sendbuffer ij_order
  icount_order(:) = 0
  ij_order(:) = 0
  len = 1

  DO iproc = 1, IPU_numpeco
     inode = iproc-1
     nfetch = 0
     DO ip = 1, nprof
        IF (node_list(ip) == inode ) THEN
            nfetch = nfetch + 1
            ij_order(len)   = i_list(ip)
            ij_order(len+1) = j_list(ip)
            len = len + 2
!           retain the info on the 'send/receive list place' of this profile
            num_from_node(ip) = nfetch
        ENDIF
     ENDDO
     icount_order(iproc) = 2 * nfetch
  ENDDO
  idispl_order(:) = 0
  DO iproc = 2, IPU_numpeco
     idispl_order(iproc) = idispl_order(iproc-1) + icount_order(iproc-1)
  ENDDO

! distribute the information on number of ordered profiles,
  CALL MPI_ALLTOALL ( icount_order, 1, IPU_integer,    &
                      icount_distr, 1, IPU_integer, IPU_comm_cart, implcode )

  if (implcode /= 0) then
     PRINT *,'MPI_ALLTOALL: errorcode: ',implcode
  endif

  idispl_distr(1) = 0
  DO iproc = 2, IPU_numpeco
     idispl_distr(iproc) = idispl_distr(iproc-1) + icount_distr(iproc-1)
  ENDDO

! allocate send/receive buffers
  istat = 0
  ALLOCATE ( ij_distr( SUM(icount_distr) ), STAT=izlocstat )
  istat = istat + izlocstat
  ALLOCATE (rsendbuf(INT(nlev*SUM(icount_distr)/2 + 0.01_wp)),  &
                                                              STAT=izlocstat)
  istat = istat + izlocstat
  ALLOCATE (rrecvbuf(INT(nlev*SUM(icount_order)/2 + 0.01_wp)),  &
                                                              STAT=izlocstat)
  istat = istat + izlocstat

  if (istat /= 0) then
     PRINT *,'ERROR in ALLOCATION of buffer needed for exchange of profiles'
  endif

! distribute all the information of needed profiles between all nodes
  CALL MPI_ALLTOALLV ( ij_order, icount_order, idispl_order, IPU_integer     &
                     , ij_distr, icount_distr, idispl_distr, IPU_integer     &
                     , IPU_comm_cart , implcode )
  if (implcode /= 0) then
     PRINT *,'MPI_ALLTOALLV(ij): Fehlercode: ',implcode
  endif

!------------------------------------------------------------------------------
! Section 2 : store the profiles requested by other processes in the
!             sending buffer
!------------------------------------------------------------------------------

  !Start  execution on accelerator
  !$acc data                                                         &
  !$acc present(field_3d,profiles)                                   &
  !$acc copyin(ij_distr,idispl_distr,icount_distr)                   &
  !$acc create(rsendbuf)

  DO iproc = 1, IPU_numpeco
    IF(icount_distr(iproc) /= 0) THEN
      ! fill all profiles to send into the node iproc into sending buffer
      !$acc parallel loop gang vector private(ind,i,j,index)
      DO ip = 1, icount_distr(iproc), 2
        ind = idispl_distr(iproc) + ip
        i = ij_distr(ind)
        j = ij_distr(ind+1)
        DO k = ktop, kbot
          index = INT(k - ktop + 1 + (idispl_distr(iproc) + ip - 1) &
                                          / 2 * nlev + 0.01_wp)
          rsendbuf(index) = field_3d(i,j,k)
        ENDDO
      ENDDO
      !$acc end parallel loop
    ENDIF
  ENDDO

  ! rsendbuf used in section 3 on CPU
  !$acc update host(rsendbuf)

!------------------------------------------------------------------------------
! Section 3 : Send profiles needed from this node to the nodes requesting them
!------------------------------------------------------------------------------

  icount_order = INT(icount_order * nlev / 2 + 0.01_wp)
  idispl_order = INT(idispl_order * nlev / 2 + 0.01_wp)
  icount_distr = INT(icount_distr * nlev / 2 + 0.01_wp)
  idispl_distr = INT(idispl_distr * nlev / 2 + 0.01_wp)

  CALL MPI_ALLTOALLV ( rsendbuf, icount_distr, idispl_distr, IPU_real     &
                     , rrecvbuf, icount_order, idispl_order, IPU_real     &
                     , IPU_comm_cart , implcode )

  if (implcode /= 0) then
     PRINT *,'MPI_ALLTOALLV(prof): Fehlercode: ',implcode
  endif

  ! Update required since profiles was initialized on GPU
  ! The sequential loop below is executed on CPU
  !$acc update host(profiles)
  DO iproc = 1, IPU_numpeco
     inode = iproc-1
     IF(icount_order(iproc) /= 0) THEN
! store the data received from this node into output array (at correct position)
        nfetch = 0
        DO ip = 1, nprof
           IF (node_list(ip) == inode) THEN
              nfetch = nfetch + 1
              DO k = ktop, kbot
                 index = k - ktop + 1 + idispl_order(iproc) + &
                                     (nfetch - 1) * nlev
                 profiles(k,ip) = rrecvbuf(index)
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !Copy back profiles to accelerator
  !$acc update device(profiles)

  !$acc end data

  DEALLOCATE ( ij_distr, STAT=izlocstat )
  DEALLOCATE ( rsendbuf, STAT=izlocstat )
  DEALLOCATE ( rrecvbuf, STAT=izlocstat )

  RETURN
#else
  CALL model_abort (-1, -1, 'MPI not available', 'exchange_profiles')
#endif

!------------------------------------------------------------------------------
! End of subroutine
!------------------------------------------------------------------------------

END SUBROUTINE exchange_profiles


!==============================================================================
! The subroutines defined below have been copied from the 3dvar module
! 'mo_mpi_dace' (3dvar/basic/mo_mpi_dace.f90) and slightly modified in order to
! run within the COSMO environment. These routines are required by the include
! files
! p_gather_derived.incf
! p_alltoall_derived.incf
! which support the MPI exchange of derived types.
!
!==============================================================================

  SUBROUTINE p_abort
  !-------------------------------------------------------------------
  ! this routine should be used instead of abort, util_abort() or STOP
  ! in all routines for proper clean up of all PEs
  !-------------------------------------------------------------------

!    EXTERNAL util_exit
#ifndef NOMPI
    INTEGER :: p_error                     ! MPI error number
    INTEGER :: nerr = 0    ! write error message to error out
#endif

#ifndef NOMPI
    CALL MPI_ABORT (MPI_COMM_WORLD, 1, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ABORT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
#else
!   CALL util_exit(1)
!    CALL abort('p_abort called')
    CALL model_abort (-1, -1, 'MPI_ABORT called', 'p_abort')
#endif

  END SUBROUTINE p_abort

!------------------------------------------------------------------------------
  SUBROUTINE gatherv_args (bufsize, sbuflen, rbuflen, groupsize, root, comm, &
                           recvcounts, p_sendcount, p_recvcounts, p_rdispls)

#ifdef __COSMO__
    USE data_parallel, ONLY:  p_pe => my_cart_id
#else
    USE mo_mpi_dace,   ONLY:  p_pe => d_pe
#endif

  !-------------------------------------------++++++++++++++++
  ! prepare arguments for MPI_GATHERV routines NOT YET TESTED
  !-------------------------------------------++++++++++++++++
  INTEGER ,INTENT(in)                   :: bufsize
  INTEGER ,INTENT(in)                   :: sbuflen
  INTEGER ,INTENT(in)                   :: rbuflen
  INTEGER ,INTENT(in)                   :: groupsize
  INTEGER ,INTENT(in)                   :: root
  INTEGER ,INTENT(in)                   :: comm
  INTEGER ,INTENT(in) ,OPTIONAL         :: recvcounts   (:)
  INTEGER ,INTENT(out)                  :: p_sendcount
  INTEGER ,POINTER                      :: p_recvcounts (:)
  INTEGER ,POINTER                      :: p_rdispls    (:)

  INTEGER :: p_error, p_int, i
  INTEGER :: nerr = 0    ! write error message to error out

#ifdef NOMPI
    !call abort ('gatherv_args: should never be called ifdef NOMPI')
  CALL model_abort (-1, -1, 'MPI_GATHER should never be called ifdef NOMPI', &
                            'gatherv_args')
#else

    ALLOCATE (p_recvcounts (groupsize))
    ALLOCATE (p_rdispls    (groupsize))

    p_sendcount = bufsize * sbuflen
    IF (PRESENT (recvcounts)) THEN
      p_recvcounts = recvcounts
    ELSE
      p_int = MPI_INTEGER                          ! redefine p_int for COSMO
      CALL MPI_GATHER (sbuflen,      1, p_int, &
                       p_recvcounts, 1, p_int, &
                       root, comm, p_error)

!MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE,
!           ROOT, COMM, IERROR)

    ENDIF
    if (p_pe == root) then
       IF (SUM (p_recvcounts) /= rbuflen) THEN
!          WRITE (nerr,'(a)') ' gatherv_args: sum (recvcounts) /= rbuflen'
          WRITE (nerr, *) ' gatherv_args: sum (recvcounts) /= rbuflen', &
                          rbuflen,sum (p_recvcounts), p_recvcounts
          CALL p_abort
       ENDIF
       p_recvcounts = bufsize * p_recvcounts

       p_rdispls (1) = 0
       DO i = 2, groupsize
          p_rdispls (i) = p_rdispls (i-1) + p_recvcounts (i-1)
       END DO
    end if

#endif

  END SUBROUTINE gatherv_args
!------------------------------------------------------------------------------

  SUBROUTINE alltoallv_args (bufsize, sbuflen, rbuflen, groupsize, comm,     &
                             sendcounts, recvcounts,                         &
                             p_sendcounts, p_sdispls, p_recvcounts, p_rdispls)

#ifdef __COSMO__
    USE data_parallel, ONLY:  p_pe => my_cart_id
#else
    USE mo_mpi_dace,   ONLY:  p_pe => d_pe
#endif

  !---------------------------------------------
  ! prepare arguments for MPI_ALLTOALLV routines
  !---------------------------------------------
  INTEGER ,INTENT(in)                   :: bufsize
  INTEGER ,INTENT(in)                   :: sbuflen
  INTEGER ,INTENT(in)                   :: rbuflen
  INTEGER ,INTENT(in)                   :: groupsize
  INTEGER ,INTENT(in)                   :: comm
  INTEGER ,INTENT(in)                   :: sendcounts   (:)
  INTEGER ,INTENT(in) ,OPTIONAL         :: recvcounts   (:)
  INTEGER ,POINTER                      :: p_sendcounts (:)
  INTEGER ,POINTER                      :: p_sdispls    (:)
  INTEGER ,POINTER                      :: p_recvcounts (:)
  INTEGER ,POINTER                      :: p_rdispls    (:)

#ifdef NOMPI
!    call abort ('alltoallv_args: should never be called ifdef NOMPI')
    CALL model_abort (-1, -1,                                            &
                      'MPI_ALLTOALL should never be called ifdef NOMPI', &
                      'gatherv_args')
#else

    INTEGER :: p_error, p_int, i
    INTEGER :: nerr = 0    ! write error message to error out

    IF (SUM (sendcounts) /= sbuflen) THEN
!      WRITE (nerr,*) p_pe,' alltoallv_args:  sum (sendcounts) /= sbuflen',&
!                                             sum (sendcounts) ,  sbuflen
      WRITE (nerr,*) p_pe,                                              &
                     ' alltoallv_args:  sum (sendcounts) /= sbuflen',   &
                     sum (sendcounts) ,  sbuflen

      CALL p_abort
    ENDIF

    ALLOCATE (p_sendcounts (groupsize))
    ALLOCATE (p_recvcounts (groupsize))
    ALLOCATE (p_sdispls    (groupsize))
    ALLOCATE (p_rdispls    (groupsize))

    p_sendcounts = bufsize * sendcounts
    IF (PRESENT (recvcounts)) THEN
      p_recvcounts = recvcounts
    ELSE
      p_int = MPI_INTEGER                          ! redefine p_int for COSMO
      CALL MPI_ALLTOALL (  sendcounts, 1, p_int, &
                         p_recvcounts, 1, p_int, &
                         comm, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_ALLTOALL failed in alltoallv_args.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      ENDIF
    ENDIF
    IF (SUM (p_recvcounts) /= rbuflen) THEN
      WRITE (nerr,'(a)') ' alltoallv_args: sum (recvcounts) /= rbuflen'

!+++++++++++++++++++++++++++++++++
!write some additional diagnostics
!+++++++++++++++++++++++++++++++++
!WRITE (nerr,*) '    p_pe               =',p_pe

WRITE (nerr,*) '    bufsize            =',bufsize
WRITE (nerr,*) '    sbuflen            =',sbuflen
WRITE (nerr,*) '    rbuflen            =',rbuflen
WRITE (nerr,*) '    groupsize          =',groupsize
WRITE (nerr,*) '    comm               =',comm
WRITE (nerr,*) '    size(sendcounts)   =',size(sendcounts)
WRITE (nerr,*) '    sendcounts         =',sendcounts

WRITE (nerr,*) '    p_sendcounts       =',p_sendcounts
WRITE (nerr,*) '    p_recvcounts       =',p_recvcounts

WRITE (nerr,*) '    PRESENT(recvcounts)=',PRESENT(recvcounts)
if(PRESENT(recvcounts)) then
WRITE (nerr,*) '    size(recvcounts)   =',size(recvcounts)
WRITE (nerr,*) '    recvcounts         =',recvcounts
endif

      CALL p_abort
    ENDIF
    p_recvcounts = bufsize * p_recvcounts

    p_sdispls (1) = 0
    p_rdispls (1) = 0
    DO i = 2, groupsize
      p_sdispls (i) = p_sdispls (i-1) + p_sendcounts (i-1)
      p_rdispls (i) = p_rdispls (i-1) + p_recvcounts (i-1)
    END DO

#endif

  END SUBROUTINE alltoallv_args

!==============================================================================
!==============================================================================

END MODULE parallel_utilities

!==============================================================================
