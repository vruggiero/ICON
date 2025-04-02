!
!+ parallelised linear algebra tool box
!
!------------------------------------------------------------------------------
#if defined (__NEC__)
! This compiler directive must appear before the MODULE statement!
!OPTION! -pvctl matmul
!
! Note: compile this module explicitly with the option -Npi,
!       otherwise sxf90 rev.360 grinds a long time using lots of memory (32G!)
#endif
!------------------------------------------------------------------------------
!
MODULE mo_dec_matrix
!
! Description:
!   This module holds the declaration of data types 't_matrix'
!   for matrices and 't_vector' for vectors in block decomposed form
!   to be stored on different processors in a distributed memory
!   environment. Arithmetic operations +, -, *, /, \dots may be
!   performed on these data types. Based on these data types and operators
!   preconditioned conjugate gradient algorithms are implemented in module
!   'mo_solver'.
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_2         2008/12/04 Harald Anlauf
!  add_matrix_to_matrix: set kind of some real constants to mp
! V1_4         2009/03/26 Andreas Rhodin
!  fix for NEC SX (loopcount); changes for LETKF
! V1_5         2009/05/25 Harald Anlauf
!  Disable ftrace regions
! V1_6         2009/06/10 Andreas Rhodin
!  new subroutine local; fixes for SX9
! V1_7         2009/08/24 Andreas Rhodin
!  do not maintain linked list for trace of memory usage
! V1_8         2009/12/09 Harald Anlauf
!  pack_block, unpack_block, diag_packed: vectorize for SX-9
! V1_9         2010/04/20 Andreas Rhodin
!  new specific routine: vector_equal_real (operator (==))
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  extend functionality of subroutine print_matrix[_block]
! V1_13        2011/11/01 Andreas Rhodin
!  new routines:         release_ivector_mem, assign_ivector_to_i4array,
!                        any(boolean_vector), scalar_minus_vector,
!                        realloc_matrix_block
!  assign_vector:        construct destination if not done
!  pack_matrix:          convert to FULL matrix representation first
!  insert_full_block:    increase memory if not sufficient
!  inverse_matrix_block: provide oprional parameter max_cond_number
! V1_15        2011/12/06 Andreas Rhodin
!  implement inverse_rs for CSC and CSR representation
! V1_19        2012-04-16 Andreas Rhodin
!  define specific 'destruct' for arrays of derived type
! V1_20        2012-06-18 Andreas Rhodin
!  new specific routine assign_scalar_to_vector
! V1_22        2013-02-13 Robin Faulwetter
!  new subroutine insert_full_block_2; adaptation to 50/51 levels for RTTOV
! V1_23        2013-03-26 Andreas Rhodin
!  new vector/vector and vector/real operators: all,<,<=,>,>=
!  single + double precision version for insert_full_block
! V1_26        2013/06/27 Andreas Rhodin
!  new specific routine gather_re (gather real vector for ensemble)
!  new subroutine assign_vectorS (fixes bug in LETKF revision 9063)
! V1_27        2013-11-08 Andreas Rhodin
!  new specific assignment operator: t_ivector = integer
! V1_28        2014/02/26 Andreas Rhodin
!  add_matrix_to_matrix: allow more combinations of matrix representations
! V1_37        2014-12-23 Harald Anlauf
!  modify_vcov: implement stdv_ratio format version 2
! V1_42        2015-06-08 Harald Anlauf
!  OpenMP parallelization
! V1_45        2015-12-15 Harald Anlauf
!  OpenMP optimization
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD  2001-2008  original code
! Oliver Schmid   DWD        2003-2004  .., LU decomposition
! Harald Anlauf   DWD        2008       optimizations for SX8
!============================================================================
!
! Diagnostics for vectorization
!
#if defined (_FTRACE) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!============================================================================
#include "tr15581.incf"
!============================================================================
!
! Maintain linked list for trace of memory usage ?
!
! #define LINKED_LIST
!============================================================================
  !-------------
  ! Modules used
  !-------------
  use mo_kind,       only: wp,sp,dp,i4,i8   ! precision kind parameters
  use mo_exception,  only: finish           ! abort routine
  use mo_allocate,   only: call_level,     &! level a routine is called at
                           enter_function, &! to be called at start of routines
                           leave_function   ! to be called at end of routines
  use mo_mpi_dace,   only: dace,           &! MPI grop info
                           p_bcast,        &! overloaded MPI broadcast routine
                           p_send,         &! overloaded MPI send routine
                           p_recv,         &! overloaded MPI receive routine
                           p_sum,          &! overloaded MPI sum routine
                           p_or,           &! overloaded MPI or  routine
                           p_and,          &! overloaded MPI and routine
                           p_max,          &! overloaded MPI maximum
                           p_min,          &! overloaded MPI minimum
#ifndef NO_MPI3
                           p_ibcast,       &! generic non-blocking MPI bcast
                           p_irequest,     &! index of next non-blocking request
                           p_mrequest,     &! maximum outstanding requests
#endif
                           p_waitall,      &! wait for MPI requests to complete
                           crc              ! checksum calculation (debugging)
  use mo_dialog,     only: ask_menue,      &! interactive user interface
                           ask              ! interactive user interface
  use mo_dace_string,only: char3            ! convert: integer -> char(len=3)
  use mo_matrix,     only: inverse,        &! invert a matrix
                           inverse_rs,     &! invert a real symmetric matrix
                           sqrt_rs,        &! square root of r.s. matrix
                           pow_rs           ! A ** x
  use slatec_module, only: sort             ! sort
  use mo_p_output,   only: oline,          &! output line buffer
                           iol,            &! index of next line to write
                           nextline,       &! routine to increment line number
                           flush_buf,      &! routine to write buffer
                           add_line_pio     ! routine to write string on I/O PE
  use mo_random,     only: random_state_t, &! random generator state type
                           destruct,       &! deallocate random generator state
                           random_number,  &! uniform random numbers
                           random_gauss     ! normal distribution
  use mo_random_seed,only: construct_stream ! new random generator, derive seed
  use mo_sparskit,   only: csrjad,         &! Convert CSR -> JAD
                           jadcsr,         &!         JAD -> CSR
                           csrcsc,         &!        CSR <-> CSC
                           amuxj            ! JAD matrix times vector

  implicit none

  !-----------------------------------
  ! External BLAS/LAPACK routines used
  !-----------------------------------
! real(wp) ,external :: ddot   ! BLAS dot_product
  real(wp) ,external :: dnrm2  ! BLAS 2-norm
  external           :: dpptrf ! LAPACK Cholesky factorisation
  external           :: dpptrs ! LAPACK (Cholesky fac.) inv. matrix multiply
  external           :: dspmv  ! BLAS   (packed) matrix vector multiplication
  external           :: dgemv  ! BLAS   matrix vector multiplication
  external           :: dgetrf ! LAPACK LU factorisation
  external           :: dgetrs ! LAPACK (LU fac.) inv. matrix multiply
!
!EOX
!
! !PUBLIC TYPES:
!
  !
  ! The following public entities are provided by module MO_DEC_MATRIX:
  !
  ! Data types T_MATRIX and T_VECTOR hold matrices vectors in block
  ! decomposed form to be stored on different processors in a distributed
  ! memory environment. For a given vector or matrix information on the
  ! distribution is stored in the data type T_DEC_INFO.
  !
!BOX
  private
!EOX
  public :: t_matrix   ! block decomposed matrix
  public :: t_vector   ! decomposed vector
  public :: t_ivector  ! decomposed integer-vector
  public :: t_bvector  ! decomposed boolean-vector
  public :: t_dec_info ! decomposition meta data
  public :: t_perm     ! meta data for CSR(PERM), JAD representation
  !
  ! The blocks or segments are stored in array components of the following
  ! types:
  !
  public :: t_matrix_block  ! matrix blocks
  public :: t_vector_segm   ! vector segments
  public :: t_ivector_segm  ! vector segments
  public :: t_bvector_segm  ! vector segments
  public :: t_info_block    ! decomposition meta data
  public :: mp, ip          ! kind parameter coefficients and indices
!
! !PUBLIC DATA MEMBERS:
!
  !
  ! Matrices (or matrix blocks) may have the following qualities
  ! (component %qual):
  !
  public :: GENERAL  ! general representation, no restrictions
  public :: BDIAG    ! block diagonal
  public :: BLOWER   ! strict lower block triangular matrix
  public :: BUPPER   ! strict upper block triangular matrix
  public :: SYMMETRIC! symmetric                     matrix
  public :: UPPER    ! strict upper       triangular matrix
  public :: LOWER    ! strict lower       triangular matrix
  !
  ! Matrix blocks may be stored in different representations (%repr):
  !
  public :: NOT_INIT ! matrix block not initialized (located on other PE)
  public :: ZERO     ! all elements are zero, no storage required
  public :: IDENT    ! identity matrix, no storage required
  public :: DIAGONAL ! diagonal matrix
  public :: FULL     ! full n x m representation
  public :: PACKED   ! only upper triangle is stored
  public :: CHOLES   ! Cholesky decomposed (inverse matrices)
  public :: LU       ! LU decomposed (inverse matrices)
  public :: CSR      ! compressed sparse row representation
  public :: CSC      ! compressed sparse column representation
  public :: JAD      ! compressed sparse JAgged Diagonal representation
  public :: MIRROR   ! this i,j block is not stored on this pe, use j,i instead
  public :: NEST     ! nested matrix
  public :: SVD      ! singular value decomposition
  public :: srep     ! string array to convert representation to mnemonic
!
! !PUBLIC MEMBER FUNCTIONS:
!
  !
  ! For initialization and deinitialization the routines CONSTRUCT and
  ! DESTRUCT are defined on variables of data types T_MATRIX,
  ! T_VECTOR and T_DEC_INFO. DELETE_STORAGE is used to deallocate
  ! pointer components of return arguments from functions. Subroutine
  ! DEC_MATRIX_MEM diagnoses the memory used.
  !
  public :: construct
  public :: construct_unity
  public :: destruct
  public :: deallocate
  public :: reallocate
  public :: delete_storage
  public :: dec_matrix_mem
  !
  public :: convert                 ! Convert between representations
  public :: set_perm                ! Set permutation information for CSR/CSC
  public :: destruct_perm           ! Clear permutation information
  !
  ! MPI communication routines
  ! The following routines send/receive/broadcast matrix blocks.
  !
  public :: p_send
  public :: p_recv
  public :: p_bcast
  !
  ! The following routines gather information on one or all PE's
  !
  public :: global             ! store vector content on all processors
  public :: local              ! store vector content on 1 processor only
  public :: gather             ! gather  vector content on 1 processor
  public :: scatter            ! scatter vector content on 1 processor
  public :: release_mem        ! release temporary allocated memory
  !
  ! Arithmetic operations are defined on variables of data types
  ! T_MATRIX, T_VECTOR and T_DEC_INFO:
  !
  public :: operator   (+)     ! vector + vector
  public :: operator   (-)     ! vector - vector
  public :: operator   (*)     ! matrix * vector; scalar * vector
  public :: operator   (/)     ! 'elemental' vector / vector
  public :: operator   (**)    ! vector ** integer, matrix ** real
  public :: operator   (==)    ! vector == vector, real
  public :: operator   (/=)    ! vector /= vector, real
  public :: operator   (< )    ! vector <  vector, real
  public :: operator   (<=)    ! vector <= vector, real
  public :: operator   (> )    ! vector >  vector, real
  public :: operator   (>=)    ! vector >= vector, real
  public :: assignment (=)     ! vector = vector,scalar; matrix = matrix
  public :: outer_product      ! outer product matrix*vector / vector*matrix
  public :: assign             ! matrix = matrix with optional parameters
  public :: insert             ! insert submatrix in block
  public :: insert_full_block_2! overwrite submatrix in block
  public :: reorder            ! reorder the elements of a matrix or vector
  public :: add_to             ! matrix = matrix + matrix
  public :: diag               ! diagonal (matrix)
  public :: dot_product        ! dot product
  public :: norm2              ! 2-norm of a vector
  public :: abs                ! 'elemental' absolute value
  public :: max                ! 'elemental' maximum value
  public :: min                ! 'elemental' minimum value
  public :: log                ! 'elemental' logarithm
  public :: amax               ! maximum of absolute value
  public :: sum                ! vector element sum
  public :: count              ! count true elements of a vector
  public :: any                ! any element  of boolean vector
  public :: all                ! all elements of boolean vector
  public :: sqrt               ! vector, matrix sqrt
  public :: pow                ! matrix A ** x
  public :: update             ! A = A + u v^t
  public :: maxval             ! maximum value of elements
  public :: minval             ! minimum value of elements
  public :: size               ! size of vector
  public :: pack_matrix        ! store matrix in packed format
  public :: full_matrix        ! store matrix in full   format
  public :: csr_matrix         ! store matrix in compressed sparse row format
  public :: csc_matrix         ! store matrix in compressed sparse row format
  public :: cholesky           ! calculate Cholesky factorisation
  public :: LU_solv            ! calculate LU factorisation
  public :: scale              ! Diagonal*Blockdiagonal*Digonalmatrix
  public :: random_number      ! set vector to random field
  public :: random_gauss       ! set vector to random field
  public :: print              ! print routine
  public :: explore            ! interactively show matrix or vector
  public :: allocate_block     ! allocate components of a matrix block
  public :: get_row            ! return one row of the matrix
  public :: get_row_sparse     ! return one row (sparse representation)
  public :: crep               ! character used for visualisation
  public :: inverse            ! inverse matrix
  public :: inverse_rs_matrix  ! inverse real symmetric matrix
  public :: sub_block          ! return submatrix
  public :: random_seed        ! set random seed
  public :: print_crc          ! print vector/matrix checksums (for debugging)
!EOP
!==============================================================================
! set pointer component memory to an invalid value before deallocation.
  logical, parameter :: setinvalid = .true.
!
! To save memory store matrix coefficients with reduced precision
  integer, parameter :: mp = sp  ! wp or sp
  integer, parameter :: ip = i4  ! i4 or i2
  integer, parameter :: mi = huge(0_ip)
!==============================================================================
  !
  ! For a given vector or matrix detailed information on the distribution
  ! over processors is stored in the data type T_DEC_INFO.
  !
  ! For each block, T_DEC_BLOCK specifies the number of vector components
  ! N and the processor element PE the vector elements and matrix diagonal
  ! blocks shall be stored on.
  !
  type t_info_block
    integer              :: n       =  0      ! number of elements
    integer              :: pe      = -1      ! processor hosting this block
    type(random_state_t) :: seed              ! random generator state
  end type t_info_block
  !
  ! T_DEC_INFO holds the total number of elements N, the number of
  ! blocks N_B as well as the array B of type T_DEC_INFO:
  !
  type t_dec_info
    integer                      :: n     =  0      ! number of elements
    integer                      :: n_b   =  0      ! number of blocks
    type (t_info_block) ,pointer :: b (:) => NULL() ! info blocks
  end type t_dec_info
  !
  ! T_PERM holds the permutation vector and metadata for the JAD and
  ! for the CSRPERM/CSCPERM representation of a CSR/CSC matrix.
  !
  type t_perm
    integer          :: jdiag      =  -1      ! Number of jagged diagonals
    integer          :: ngroup     =  -1      ! Number of groups (CSRPERM)
    integer          :: m          =   0      ! Number of rows
    integer          :: maxnz      =   0      ! Maximum no. of nonzeros/row
    integer _POINTER :: iperm  (:) => NULL () ! Permutation vector (1:m)
    integer _POINTER :: ia     (:) => NULL () ! Index to jagged d. (1:jdiag+1)
    integer _POINTER :: nzgroup(:) => NULL () ! Nonzeros in group  (1:ngroup)
    integer _POINTER :: xgroup (:) => NULL () ! Group index        (1:ngroup+1)
  end type t_perm
!==============================================================================
!BOP
!
! !DATATYPE: t_matrix
!
! !DESCRIPTION:
!
  ! The data type 't_matrix' holds a symmetric 'n'x'n'
  ! matrix decomposed in 'n_b' times 'n_b' blocks 'b(:,:)'
  ! of data type 't_matrix'. More specific information on the
  ! decompostition may be retrieved by the reference 'info'. Special
  ! characteristics of the matrix (e.g. blog diagonal matrices) are
  ! indicated by the component 'qual' with allowed values given by the
  ! predefined constants 'GENERAL' or 'BDIAG'. The latter value
  ! indicates that only diagonal blocks 'b'(i,i) are present.
  !
  ! Each matrix block is stored in a variable of type 't_matrix_block'
  ! defining the number of rows 'm' and columns 'n' within this
  ! block as well as the number of nonzero elements 'nonzero'. (For
  ! nonzero==0 no further storage is required).  Matrix blocks are
  ! distributed over processor elements in a parallel environment.  The
  ! index of the processor element hosting this block is specified in {\bf
  ! pe'. Currently the parallisation strategy is to hold all elements of a
  ! row on the same processor element.  Depending on the value of the
  ! representation flag 'repr' the matrix elements are stored in
  ! different representations:
  !
  ! \begin{description}
  ! \item'ZERO}:\\
  ! All elements of the block are zero, no storage is required.
  ! \item'IDENT}:\\
  ! The matrix block is the identity matrix, no storage is required.
  ! \item'DIAGONAL}:\\
  ! The matrix block is a diagonal matrix. The diagonal is stored in
  ! component 'packed(1:n)}.
  ! \item'FULL}:\\
  ! The block is stored in full representation in component {\bf
  ! full(1:n,1:m)}.
  ! \item'PACKED}:\\
  ! The block is symmetric, only the upper triangle is stored in component
  ! 'packed}
  ! \item'CHOLES}:\\
  ! The block represents the inverse of a matrix block, stored in Cholesky
  ! decomposed form in component 'packed(:)}
  ! \item'LU}:\\
  ! The block represents the inverse of a matrix block, stored in LU
  ! decomposed form in component 'full(:,:)}
  ! \item'CSR, CSC}:\\
  ! If the number of nonzero elements is small, the matrix may be stored
  ! in a compressed sparse row or column representation using {\bf
  ! packed(:)}, 'ia(:)} and 'ja(:)}.
  ! \item'MIRROR}:\\
  ! Because of symmetry, this block 'b(i,j)} is not stored, the
  ! transpose of 'b(j,i)} shall be used instead.
  ! \end{description}
  !
! !DEFINITION:
  !
  ! The matrix may have the following characteristics (% qual):
  !
  integer ,parameter :: GENERAL  = 0 ! general representation, no restrictions
  integer ,parameter :: BDIAG    = 1 ! block diagonal representation
  integer ,parameter :: BUPPER   = 2 ! strict upper block triangular matrix
  integer ,parameter :: BLOWER   = 3 ! strict lower block triangular matrix
  integer ,parameter :: SYMMETRIC= 4 ! symmetric                     matrix
  integer ,parameter :: UPPER    = 5 ! strict upper       triangular matrix
  integer ,parameter :: LOWER    = 6 ! strict lower       triangular matrix
  !
  ! Each block may be represented in one of the following forms (% repr)
  !
  integer  ,parameter :: NOT_INIT= -1 ! not initialised
  integer  ,parameter :: ZERO    =  0 ! all elements are zero, not stored
  integer  ,parameter :: IDENT   =  1 ! identity matrix, no storage required
  integer  ,parameter :: DIAGONAL=  2 ! diagonal matrix
  integer  ,parameter :: FULL    =  3 ! full n x m representation
  integer  ,parameter :: PACKED  =  4 ! only upper or lower triangle is stored
  integer  ,parameter :: CHOLES  =  5 ! Cholesky decomposed (inverse matrices)
  integer  ,parameter :: CSR     =  6 ! compressed sparse row representation
  integer  ,parameter :: CSC     =  7 ! compressed sparse column representation
  integer  ,parameter :: MIRROR  =  9 ! this i,j block is not stored, use j,i
  integer  ,parameter :: NEST    = 10 ! nested matrices
  integer  ,parameter :: SVD     = 11 ! singular value decomposition
  integer  ,parameter :: LU      = 12 ! LU decomposition
  integer  ,parameter :: JAD     = 13 ! JAgged Diagonal representation
  !
  type t_matrix_block
    integer             :: m          =  0       ! number of rows
    integer             :: n          =  0       ! number of columns
    integer             :: nonzero    = -1       ! number of nonzero elements
    integer             :: pe         = -1       ! processor hosting this block
    integer             :: repr       = -1       ! representation (FULL,..)
    character           :: tri        = ''       ! triangle ('U' or 'L') stored
    integer             :: alloc_l    = 0
    integer             :: qual       = GENERAL  ! SYMMETRIC UPPER LOWER
    real(mp)   _POINTER :: full (:,:) => NULL()  ! coefficients   FULL   repr.
    real(mp)   _POINTER :: packed (:) => NULL()  ! coefficients   PACKED repr.
    integer    _POINTER :: ia     (:) => NULL()  ! row indices    CSR    repr.
    integer(ip)_POINTER :: ja     (:) => NULL()  ! column indices CSR    repr.
    type(t_matrix) ,pointer :: nest   => NULL()  ! nested matrices
    type(t_perm)   ,pointer :: perm   => NULL()  ! JAD/CSRPERM/CSCPERM metadata
  end type t_matrix_block
  !
  type t_matrix
    type(t_dec_info)     ,pointer :: rinfo   => NULL()  ! decomp. info (rows)
    type(t_dec_info)     ,pointer :: cinfo   => NULL()  ! decomp. info (cols)
    character(len=8)              :: name    = '???? '  ! name
    integer                       :: m       =  0       ! number of rows
    integer                       :: n       =  0       ! number of cols
    integer                       :: m_b     =  0       ! number of row blocks
    integer                       :: n_b     =  0       ! number of column bl.
    integer                       :: qual    =  GENERAL ! qualities
    integer                       :: alloc_l =  0       ! allocation level
    logical                       :: rowwise = .true.   ! rowwise storage
    type(t_matrix_block) ,pointer :: b (:,:) => NULL()  ! matrix blocks
    type(t_matrix)       ,pointer :: next    => NULL()  ! linked list for
    type(t_matrix)       ,pointer :: prev    => NULL()  !   memory usage diagn.
  end type t_matrix
!EOP
  character,parameter :: crep (-1:13) = &
    (/'@',' ','1','\\'(1:1),'+','x','I','-','|','?','.','n','s','%','j'/)
  character(len=8) ,parameter :: srep (-1:13) = &
    (/'NOT_INIT','ZERO    ','IDENT   ','DIAGONAL','FULL    ','PACKED  ',&
      'CHOLES  ','CSR     ','CSC     ','??????  ','MIRROR  ','NEST    ',&
      'SVD     ','LU      ','JAD     '/)
!==============================================================================
  !
  ! Vectors are decomposed into segments corresponding to the respective
  ! blocks of decomposed matrices. The N coefficients of a segment are
  ! stored on processor element PE within the component X(:). More
  ! specific information on the decompostition may be retrieved by the
  ! reference INFO to data type T_DEC_INFO.
  !
  type t_vector_segm
    integer           :: pe   = -1      ! processor hosting this block
    integer           :: n    =  0      ! number of elements
    real(wp) _POINTER :: x(:) => NULL() ! coefficients
  end type t_vector_segm
  !
  ! The segments of the decomposed vector are kept in the component S(:)
  ! of the data type T_VECTOR. N and specifies the total number of
  ! elements and N_S the number of segments. ALLOC_L indicates if the
  ! vector is a local variable within a subroutine or function. This
  ! information is required for the deallocation of pointer componebts in
  ! function return arguments. Arrays segments may be either distributed
  ! over processors ore stored rdundantly on each processor. The latter
  ! case is indicated by GLOBAL=.false.
  !
  type t_vector
    type (t_dec_info)    ,pointer :: info    => NULL()  ! decomposition info
    character(len=8)              :: name    = '???? '  ! name
    integer                       :: n       =  0       ! number of elements
    integer                       :: n_s     =  0       ! number of segments
    integer                       :: alloc_l =  0       ! allocation level
    logical                       :: global  = .false.  ! global allocation
    type (t_vector_segm) ,pointer :: s (:)   => NULL()  ! vector blocks
    type (t_vector)      ,pointer :: next    => NULL()  ! linked list for
    type (t_vector)      ,pointer :: prev    => NULL()  !   memory usage diagn.
  end type t_vector
!==============================================================================
  !
  ! Same as type t_vector but for an integer vector
  !
  type t_ivector_segm
    integer          :: pe   = -1      ! processor hosting this block
    integer          :: n    =  0      ! number of elements
    integer _POINTER :: x(:) => NULL() ! coefficients
  end type t_ivector_segm
  !
  type t_ivector
    type (t_dec_info)     ,pointer :: info    => NULL()  ! decomposition info
    character(len=8)               :: name    = '???? '  ! name
    integer                        :: n       =  0       ! number of elements
    integer                        :: n_s     =  0       ! number of segments
    integer                        :: alloc_l =  0       ! allocation level
    logical                        :: global  = .false.  ! global allocation
    type (t_ivector_segm) ,pointer :: s (:)   => NULL()  ! vector blocks
    type (t_ivector)      ,pointer :: next    => NULL()  ! linked list for
    type (t_ivector)      ,pointer :: prev    => NULL()  !   memory usage diagn.
  end type t_ivector
!==============================================================================
  !
  ! Same as type t_vector but for a boolean vector
  !
  type t_bvector_segm
    integer          :: pe   = -1      ! processor hosting this block
    integer          :: n    =  0      ! number of elements
    logical _POINTER :: x(:) => NULL() ! coefficients
  end type t_bvector_segm
  !
  type t_bvector
    type (t_dec_info)     ,pointer :: info    => NULL()  ! decomposition info
    character(len=8)               :: name    = '???? '  ! name
    integer                        :: n       =  0       ! number of elements
    integer                        :: n_s     =  0       ! number of segments
    integer                        :: alloc_l =  0       ! allocation level
    logical                        :: global  = .false.  ! global allocation
    type (t_bvector_segm) ,pointer :: s (:)   => NULL()  ! vector blocks
    type (t_bvector)      ,pointer :: next    => NULL()  ! linked list for
    type (t_bvector)      ,pointer :: prev    => NULL()  !   memory usage diagn.
  end type t_bvector
!==============================================================================
  !
  ! Private entities for accounting of allocated memory
  !
  type (t_matrix)  ,pointer :: first_m => NULL()
  type (t_matrix)  ,pointer :: last_m  => NULL()
  type (t_vector)  ,pointer :: first_v => NULL()
  type (t_vector)  ,pointer :: last_v  => NULL()
  type (t_ivector) ,pointer :: first_i => NULL()
  type (t_ivector) ,pointer :: last_i  => NULL()
  type (t_bvector) ,pointer :: first_b => NULL()
  type (t_bvector) ,pointer :: last_b  => NULL()
!==============================================================================
  !
  ! Interfaces for the public routines and operators are defined below.
  !
  interface construct
    module procedure construct_matrix
    module procedure construct_vector
    module procedure construct_vectors
    module procedure construct_ivector
    module procedure construct_bvector
    module procedure construct_dec_info
    module procedure construct_info_block
    module procedure construct_matrix_block
  end interface construct

  interface construct_unity
    module procedure construct_ident
  end interface construct_unity

  interface delete_storage
    module procedure delete_dec_vector
    module procedure delete_dec_integer
    module procedure delete_dec_boolean
    module procedure delete_dec_matrix
    module procedure delete_matrix_block
  end interface delete_storage

  interface deallocate
    module procedure dealloc_matrix_block
    module procedure dealloc_perm
  end interface deallocate

  interface reallocate
    module procedure realloc_matrix_block
  end interface reallocate

  interface destruct
    module procedure destruct_matrix
    module procedure destruct_matrix_block
    module procedure destruct_vector
    module procedure destruct_vectors
    module procedure destruct_vectors2
    module procedure destruct_ivector
    module procedure destruct_bvector
    module procedure destruct_dec_info
    module procedure destruct_perm_
  end interface destruct

  interface set_invalid
    module procedure set_invalid_matrix_block
  end interface set_invalid

  interface assignment (=)
    module procedure assign_vector
    module procedure assign_vectors
    module procedure assign_ivector
    module procedure assign_bvector
    module procedure assign_vector_to_dparray
    module procedure assign_vector_to_sparray
    module procedure assign_ivector_to_i4array
    module procedure assign_scalar_to_vector
    module procedure assign_scalar_to_vectors
    module procedure assign_scalar_to_bvector
    module procedure assign_int_to_ivector
    module procedure assign_matrix_to_matrix
    module procedure assign_matrix_block
    module procedure assign_matrix_block_array
  end interface assignment (=)

  interface assign
    module procedure assign_matrix_optional
  end interface assign

  interface global
    module procedure global_0 ! scalar  version
    module procedure global_1 ! vector  version
    module procedure global_b ! logical version
  end interface global

  interface local
    module procedure local_0 ! scalar version
    module procedure local_1 ! vector version
  end interface local

  interface insert
    module procedure insert_full_block
    module procedure insert_full_dp
    module procedure insert_sparse_block
  end interface insert

  interface reorder
    module procedure reorder_vector
    module procedure reorder_matrix
  end interface

  interface add_to
    module procedure add_matrix_to_matrix
    module procedure add_diagonal_to_matrix
  end interface add_to

  interface random_gauss
    module procedure random_gauss_vector
  end interface random_gauss

  interface operator (+)
    module procedure add_vectors
  end interface operator (+)

  interface operator (-)
    module procedure subtract_vectors
    module procedure vector_minus_scalar
    module procedure scalar_minus_vector
    module procedure minus_vector
  end interface operator (-)

  interface operator (*)
    module procedure matrix_times_vector
    module procedure vector_times_matrix
    module procedure matrix_times_vector_b
    module procedure vector_times_matrix_b
    module procedure scalar_times_vector
    module procedure vector_times_vector
    module procedure matrix_times_matrix
    module procedure scalar_times_matrix
  end interface operator (*)

  interface operator (/)
    module procedure scalar_over_vector
    module procedure vector_over_vector
    module procedure vector_over_scalar
    module procedure vector_over_int
  end interface operator (/)

  interface operator (**)
    module procedure vector_exp_integer
    module procedure pow_rs_matrix_op
  end interface operator (**)

  interface operator (==)
    module procedure vector_equal_vector
    module procedure vector_equal_real
  end interface operator (==)

  interface operator (/=)
    module procedure vector_ne_vector
    module procedure vector_ne_real
  end interface operator (/=)

  interface operator (<)
    module procedure vector_lt_vector
    module procedure vector_lt_real
  end interface operator (<)

  interface operator (<=)
    module procedure vector_le_vector
    module procedure vector_le_real
  end interface operator (<=)

  interface operator (>)
    module procedure vector_gt_vector
    module procedure vector_gt_real
  end interface operator (>)

  interface operator (>=)
    module procedure vector_ge_vector
    module procedure vector_ge_real
  end interface operator (>=)

  interface pow
    module procedure pow_rs_matrix
  end interface pow

  interface update
    module procedure update_matrix_block
  end interface

  interface abs
    module procedure abs_vector
  end interface abs

  interface log
    module procedure log_vector
  end interface log

  interface max
    module procedure max_vector_real
    module procedure max_vector_vector
  end interface max

  interface min
    module procedure min_vector_real
    module procedure min_vector_vector
  end interface min

  interface pack_matrix
    module procedure pack_dec_matrix
    module procedure pack_block
  end interface pack_matrix

  interface full_matrix
    module procedure unpack_block
  end interface full_matrix

  interface csr_matrix
    module procedure csr_block
  end interface csr_matrix

  interface csc_matrix
    module procedure csc_block
  end interface csc_matrix

  interface set_perm
    module procedure set_perm_block         ! Set permutation information
  end interface

  interface destruct_perm
    module procedure destruct_perm_block
  end interface

  interface convert
    module procedure convert_matrix_block   ! Convert between representations
  end interface

  interface random_number
    module procedure random_vector
  end interface random_number

  interface dot_product
    module procedure dot_vector
  end interface dot_product

  interface sum
    module procedure sum_vector
  end interface sum

  interface count
    module procedure count_vector
  end interface count

  interface any
    module procedure any_vector
  end interface any

  interface all
    module procedure all_vector
  end interface all

  interface sqrt
    module procedure sqrt_vector
    module procedure sqrt_rs_matrix
    module procedure sqrt_rs_matrix_block
  end interface sqrt

  interface maxval
    module procedure maxval_vector
  end interface maxval

  interface minval
    module procedure minval_vector
  end interface minval

  interface size
    module procedure size_of_vector
  end interface size

  interface print
    module procedure print_matrix_block
    module procedure print_vector
    module procedure print_matrix
  end interface print

  interface explore
    module procedure explore_vector
    module procedure explore_matrix
  end interface

  interface diag
    module procedure diag_matrix
    module procedure diag_matrix_block
  end interface diag

  interface inverse
    module procedure inverse_matrix_block
    module procedure inverse_matrix
  end interface inverse

  interface cholesky
    module procedure cholesky_block
    module procedure cholesky
  end interface cholesky

  interface scale
    module procedure scale_matrix
    module procedure scale_block
  end interface

!+++++++++++++++++++++++++++++++++++++++++
! interface 'inverse_rs' did not work with
! NAG compiler (segmentation fault) :
!+++++++++++++++++++++++++++++++++++++++++
  interface inverse_rs_matrix
    module procedure inverse_rs_matrix_block
    module procedure inverse_rs_matrix
  end interface inverse_rs_matrix

  interface p_send
    module procedure send_matrix_block
    module procedure send_perm
  end interface p_send

  interface p_recv
    module procedure recv_matrix_block
    module procedure recv_perm
  end interface p_recv

  interface p_bcast
    module procedure bcast_matrix_block
    module procedure bcast_perm
  end interface p_bcast

  interface gather
    module procedure gather_r   ! real
    module procedure gather_re  ! real vector
    module procedure gather_b   ! boolean
    module procedure gather_i   ! integer
  end interface gather

  interface scatter
    module procedure scatter_r  ! real
    module procedure scatter_re ! real vector
  end interface scatter

  interface release_mem
    module procedure release_vector_mem
    module procedure release_vector_mem_1
    module procedure release_ivector_mem
    module procedure release_bvector_mem
    module procedure release_mblock_mem
  end interface release_mem

  interface random_seed
    module procedure put_seed
  end interface random_seed

  interface outer_product
    module procedure outer_matrix_vector
    module procedure outer_vector_matrix
  end interface outer_product

  interface print_crc
    module procedure print_crc_vector
    module procedure print_crc_matrix
  end interface print_crc
  
!==============================================================================
contains
!==============================================================================
  !
  ! Initialization (CONSTRUCT) routines to initialize variables of derived
  ! data types before their first use are defined below.
  !
!------------------------------------------------------------------------------
  subroutine construct_dec_info (z, n, n_b, nb, pe)
  type (t_dec_info) ,intent(out) :: z
  integer           ,intent(in)  :: n      ! number of elements
  integer           ,intent(in)  :: n_b    ! number of blocks
  integer ,optional ,intent(in)  :: nb (:) ! number of elements per block
  integer ,optional ,intent(in)  :: pe (:) ! processor element  of  block
    integer :: i
    z% n   = n
    z% n_b = n_b
    allocate (z% b (n_b))
    if (present(nb) .and. present(pe)) then
      do i=1,n_b
        call construct (z% b(i), nb(i), pe(i))
      end do
    endif
  end subroutine construct_dec_info
!------------------------------------------------------------------------------
  subroutine construct_info_block (z, n, pe)
  type (t_info_block) ,intent(out) :: z
  integer             ,intent(in)  :: n   ! number of elements in this block
  integer             ,intent(in)  :: pe  ! processor
    z% n   = n
    z% pe  = pe
  end subroutine construct_info_block
!------------------------------------------------------------------------------
  subroutine construct_ident (a, ri)
  type (t_matrix)   ,intent(out),target   :: a
  type (t_dec_info) ,pointer              :: ri
    type (t_matrix_block) ,pointer :: b
    integer                        :: i, j
    a% rinfo   => ri
    a% cinfo   => ri
    a% m       =  a% rinfo% n
    a% n       =  a% cinfo% n
    a% m_b     =  a% rinfo% n_b
    a% n_b     =  a% cinfo% n_b
    a% qual    =  BDIAG
    a% name    =  'I'
    a% alloc_l =  call_level
    allocate (a% b (a% m_b, a% n_b))
    do j = 1, a% n_b
      do i = 1, a% m_b
        b => a% b (i, j)
        b% m  = a% rinfo% b(i) % n
        b% n  = a% cinfo% b(j) % n
        b% pe = a% rinfo% b(i) % pe          ! associate rows with PEs
        if (b% pe == dace% pe) then
          if (i==j) then
            call allocate_block (b, IDENT, call_l= a% alloc_l)
          else
            call allocate_block (b, ZERO,  call_l= a% alloc_l)
          endif
        endif
      end do
    end do
#if defined (LINKED_LIST)
    if (associated (last_m)) then
      last_m% next => a
      a%      prev => last_m
      last_m       => a
    else
      first_m => a
      last_m  => a
    endif
#endif
  end subroutine construct_ident
!------------------------------------------------------------------------------
  !
  ! Initialization routine for a block decomposed matrix. Dimension and
  ! decomposition information is taken from the INFO variable. Arrays on
  ! the block level are not allocated yet.
  !
  subroutine construct_matrix (a, ri, name, qual, ci)
  type (t_matrix)   ,intent(out),target   :: a
  type (t_dec_info) ,pointer              :: ri
  character(len=*)  ,intent(in) ,optional :: name
  integer           ,intent(in) ,optional :: qual
  type (t_dec_info) ,pointer    ,optional :: ci
    type (t_matrix_block) ,pointer :: b
    integer                        :: i, j
    a% rinfo   => ri
    a% cinfo   => ri;        if (present(ci))    a% cinfo => ci
    a% m       =  a% rinfo% n
    a% n       =  a% cinfo% n
    a% m_b     =  a% rinfo% n_b
    a% n_b     =  a% cinfo% n_b
    a% qual    =  GENERAL;   if (present(qual))  a% qual  = qual
    a% name    =  'matrix?'; if (present(name))  a% name  = name
    a% alloc_l =  call_level
    allocate (a% b (a% m_b, a% n_b))
    do j = 1, a% n_b
      do i = 1, a% m_b
        b => a% b (i, j)
        b% m  = a% rinfo% b(i) % n
        b% n  = a% cinfo% b(j) % n
        b% pe = a% rinfo% b(i) % pe          ! associate rows with PEs
        if (b% pe == dace% pe) then
          if ((a% qual==BDIAG  .and. i/=j) .or.   &
              (a% qual==BUPPER .and. i<=j) .or.   &
              (a% qual==BLOWER .and. i>=j)      ) &
            call allocate_block (b, ZERO, dealloc=.false.)
        endif
      end do
    end do
#if defined (LINKED_LIST)
    if (associated (last_m)) then
      last_m% next => a
      a%      prev => last_m
      last_m       => a
    else
      first_m => a
      last_m  => a
    endif
#endif
  end subroutine construct_matrix
!------------------------------------------------------------------------------
  subroutine realloc_matrix_block (b, m, n, ns, call_l)
  type (t_matrix_block) ,intent(inout) :: b       ! matrix block to init.
  integer     ,optional ,intent(in)    :: m       ! number of rows
  integer     ,optional ,intent(in)    :: n       ! number of columns
  integer     ,optional ,intent(in)    :: ns      ! size of b% packed (CSR,CSC)
  integer     ,optional ,intent(in)    :: call_l  ! value for allocation level
  !---------------------------------------------
  ! re-allocate matrix block with different size
  !---------------------------------------------

    real(mp) ,pointer :: packed(:)
    real(mp) ,pointer :: full  (:,:)
    integer  ,pointer :: ia    (:)
    integer  ,pointer :: ja    (:)
    integer           :: level
    integer           :: repr
    integer           :: l,k

    !----------------------------
    ! keep old pointer components
    !----------------------------
    packed => b% packed ;nullify (b% packed)
    full   => b% full   ;nullify (b% full)
    ia     => b% ia     ;nullify (b% ia)
    ja     => b% ja     ;nullify (b% ja)
    if (associated(b% nest)) call finish('reallocate','associated (nest)')
    if (associated(b% perm)) call finish('reallocate','associated (perm)')

    !---------
    ! allocate
    !---------
    level = b% alloc_l ;if (present(call_l)) level = call_l
    repr  = b% repr
    call allocate_block (b, repr, m=m, n=n, ns=ns, call_l=level)

    !------------------------
    ! copy pointer components
    !------------------------
    if (associated (packed)) then
      l = min (size(packed), size(b% packed))
      b% packed (1:l)  = packed (1:l)
      b% packed (l+1:) = 0._mp
    endif
    if (associated (ia)) then
      l = min (size(ia), size(b% ia))
      b% ia (1:l)  = ia (1:l)
      b% ia (l+1:) = 0
    endif
    if (associated (ja)) then
      l = min (size(ja), size(b% ja))
      b% ja (1:l)  = ja (1:l)
      b% ja (l+1:) = 0
    endif
    if (associated (full)) then
      l = min (size(full,1), size(b% full,1))
      k = min (size(full,2), size(b% full,2))
      b% full (1:l ,1:k ) = full (1:l,1:k)
      b% full (l+1:, :  ) = 0._mp
      b% full (1:l ,k+1:) = 0._mp
    endif

  end subroutine realloc_matrix_block
!------------------------------------------------------------------------------
  !
  ! Initialisation routine for a matrix block. Pointer components are
  ! allocated according to the information in the components M, N and the
  ! argument repr
  !
  subroutine allocate_block (b, repr, m, n, ns, perm, dealloc, call_l)
  type (t_matrix_block) ,intent(inout) :: b       ! matrix block to init.
  integer               ,intent(in)    :: repr    ! representation
  integer     ,optional ,intent(in)    :: m       ! number of rows
  integer     ,optional ,intent(in)    :: n       ! number of columns
  integer     ,optional ,intent(in)    :: ns      ! size of b% packed (CSR,CSC)
  logical     ,optional ,intent(in)    :: perm    ! allocate perm  (default=F)
  logical     ,optional ,intent(in)    :: dealloc ! first dealloc. (default=T)
  integer     ,optional ,intent(in)    :: call_l  ! value for allocation level

    logical :: deal
    logical :: notsquare
    logical :: alloc_perm
    !-------------------------------------------
    ! deallocate pointer components if requested
    !-------------------------------------------
    deal = .true.; if(present(dealloc)) deal = dealloc
    if(deal) call deallocate (b)
    !-----------------------------------------------
    ! Allocate derived type for permutation metadata
    !-----------------------------------------------
    alloc_perm = .false.; if(present(perm)) alloc_perm = perm
    if (alloc_perm) allocate (b% perm)

    if (present(m)) b% m = m
    if (present(n)) b% n = n
    b% alloc_l = call_level ;if (present(call_l)) b% alloc_l = call_l
    b% pe      = dace% pe
    b% repr    = repr
    nullify (b% packed)
    nullify (b% full)
    nullify (b% ia)
    nullify (b% ja)
    nullify (b% nest)
    notsquare = b% m /= b% n
    if (b% m > huge(0_ip) .or. b% n > huge(0_ip)) &
      call finish ('allocate_block','ia,ja representation insufficient')
    select case (repr)
    case (IDENT,DIAGONAL,PACKED)
      if(notsquare) call finish('allocate_block',&
                                'matrix is not square: '//srep(repr))
    end select
    select case (repr)
    case (NOT_INIT)
      b% nonzero = 0
    case (ZERO)
      b% nonzero = 0
    case (IDENT)
    case (DIAGONAL)
      allocate (b% packed (b%m))
      b% packed = 0._mp
    case (FULL)
      allocate (b% full (b%m,b%n))
      b% full = 0._mp
    case (PACKED)
      allocate (b% packed (b%m*(b%m+1)/2))
      b% packed = 0._mp
    case (NEST)
      allocate (b% nest)
    case (CSR)
      if (.not.present(ns)) call finish ('allocate_block','ns not present')
      b% nonzero = ns
      allocate (b% packed (ns))
      allocate (b% ja (ns))
      allocate (b% ia (b%m+1))
      b% ia(:) = 1
      if (huge(0_ip)<b% n) call finish('allocate_block',            &
                                       'kind parameter ip too small')
    case (CSC)
      if (.not.present(ns)) call finish ('allocate_block','ns not present')
      b% nonzero = ns
      allocate (b% packed (ns))
      allocate (b% ja (ns))
      allocate (b% ia (b%n+1))
      b% ia(:) = 1
      if (huge(0_ip)<b% m) call finish('allocate_block',            &
                                       'kind parameter ip too small')
    case (JAD)
      if (.not.present(ns)) call finish ('allocate_block','ns not present')
      b% nonzero = ns
      allocate (b% packed (ns))
      allocate (b% ja (ns))
      if (huge(0_ip)<b% m) call finish('allocate_block',            &
                                       'kind parameter ip too small')
      ! Always allocate derived type for metadata
      if (.not. associated (b% perm)) allocate (b% perm)
    case default
      call finish ('allocate_block','unknown flag '//srep(repr))
    end select
  end subroutine allocate_block
!------------------------------------------------------------------------------
  subroutine construct_matrix_block (b)
  type (t_matrix_block) ,intent(out) :: b
  end subroutine construct_matrix_block
!------------------------------------------------------------------------------
  ! Initialization routine for the permutation metadata
  subroutine allocate_perm (perm, m, jdiag, ngroup, dealloc)
    type (t_perm)     ,intent(inout) :: perm    ! metadata
    integer ,optional ,intent(in)    :: m       ! number of rows
    integer ,optional ,intent(in)    :: jdiag   ! number of jagged diagonals
    integer ,optional ,intent(in)    :: ngroup  ! number of groups (CSRPERM)
    logical ,optional ,intent(in)    :: dealloc ! first dealloc. (default=T)
    !-------------------------------------------
    ! deallocate pointer components if requested
    !-------------------------------------------
    logical :: deal
    deal = .true.; if (present (dealloc)) deal = dealloc
    if(deal) then
       call deallocate (perm)
    else
       nullify (perm% iperm)
       nullify (perm% ia)
       nullify (perm% nzgroup)
       nullify (perm% xgroup)
    end if

    if (present (m)     ) perm% m      = m
    if (present (jdiag) ) perm% jdiag  = jdiag
    if (present (ngroup)) perm% ngroup = ngroup
!write (0,*) "pe=", dace% pe," allocate_perm:", perm% m, perm% jdiag, perm% ngroup
    if (perm% m == 0) return
    if (perm% m <  0) call finish ("allocate_perm","m < 0")
    if (perm% jdiag < 0 .and. perm% ngroup < 0) return

    allocate (perm% iperm (perm% m))
    perm% iperm = 0

    if (perm% jdiag >= 0) then
       allocate (perm% ia(perm% jdiag+1))
       perm% ia = 1
    end if

    if (perm% ngroup > 0) then
       allocate (perm% nzgroup(perm% ngroup))
       allocate (perm% xgroup (perm% ngroup+1))
       perm% nzgroup = 0
       perm% xgroup  = 1
    end if
  end subroutine allocate_perm
!------------------------------------------------------------------------------
  !
  ! Initialization routine for a distributed vector. Dimension and
  ! decomposition information is taken from the INFO variable. Arrays are
  ! allocated.
  !
  subroutine construct_vector (x, info, name, global)
  type (t_vector)   ,intent(out),target   :: x
  type (t_dec_info) ,pointer              :: info
  character(len=*)  ,intent(in) ,optional :: name
  logical           ,intent(in) ,optional :: global
    type (t_vector_segm) ,pointer :: s
    integer                       :: i
    !
    ! set components
    !
    x% info    => info
    x% n       =  info% n
    x% n_s     =  info% n_b
    x% alloc_l =  call_level
    x% global  =  .false.; if (present(global)) x% global = global
    x% name =  'vector?' ; if (present(name))   x% name   = name
    if (dace% npe == 1) x% global = .true.
    !
    ! allocate and set pointer components
    !
    allocate (x% s (x% n_s))
    do i = 1, x% n_s
      s => x% s(i)
      s% n  = info% b(i)% n
      s% pe = info% b(i)% pe
      if (s% pe == dace% pe .or. x% global) allocate (s% x (s% n))
    end do
    !
    ! maintain linked list
    !
#if defined (LINKED_LIST)
    if (x% alloc_l == 0) then
      if (associated (last_v)) then
        last_v% next => x
        x%      prev => last_v
        last_v       => x
      else
        first_v => x
        last_v  => x
      endif
    endif
#endif
  end subroutine construct_vector
!------------------------------------------------------------------------------
  subroutine construct_vectors (x, info, name, global)
  type (t_vector)   ,intent(out)          :: x (:)
  type (t_dec_info) ,pointer              :: info
  character(len=*)  ,intent(in) ,optional :: name
  logical           ,intent(in) ,optional :: global
  !-------------------------------
  ! Initialize an array of vectors
  !-------------------------------
    integer :: i
    do i=1,size(x)
      call construct (x(i), info, name, global)
    end do
  end subroutine construct_vectors
!------------------------------------------------------------------------------
  !
  ! Initialization routine for a distributed ivector. Dimension and
  ! decomposition information is taken from the INFO variable. Arrays are
  ! allocated.
  !
  subroutine construct_ivector (x, info, name, global)
  type (t_ivector)  ,intent(out),target   :: x
  type (t_dec_info) ,pointer              :: info
  character(len=*)  ,intent(in) ,optional :: name
  logical           ,intent(in) ,optional :: global
    type (t_ivector_segm) ,pointer :: s
    integer                       :: i
    !
    ! set components
    !
    x% info    => info
    x% n       =  info% n
    x% n_s     =  info% n_b
    x% alloc_l =  call_level
    x% global  =  .false.; if (present(global)) x% global = global
    x% name =  'vector?' ; if (present(name))   x% name   = name
    if (dace% npe == 1) x% global = .true.
    !
    ! allocate and set pointer components
    !
    allocate (x% s (x% n_s))
    do i = 1, x% n_s
      s => x% s(i)
      s% n  = info% b(i)% n
      s% pe = info% b(i)% pe
      if (s% pe == dace% pe .or. x% global) allocate (s% x (s% n))
    end do
    !
    ! maintain linked list
    !
#if defined (LINKED_LIST)
    if (x% alloc_l == 0) then
      if (associated (last_i)) then
        last_i% next => x
        x%      prev => last_i
        last_i       => x
      else
        first_i => x
        last_i  => x
      endif
    endif
#endif
  end subroutine construct_ivector
!------------------------------------------------------------------------------
  !
  ! Initialization routine for a distributed bvector. Dimension and
  ! decomposition information is taken from the INFO variable. Arrays are
  ! allocated.
  !
  subroutine construct_bvector (x, info, name, global)
  type (t_bvector)  ,intent(out),target   :: x
  type (t_dec_info) ,pointer              :: info
  character(len=*)  ,intent(in) ,optional :: name
  logical           ,intent(in) ,optional :: global
    type (t_bvector_segm) ,pointer :: s
    integer                       :: i
    !
    ! set components
    !
    x% info    => info
    x% n       =  info% n
    x% n_s     =  info% n_b
    x% alloc_l =  call_level
    x% global  =  .false.; if (present(global)) x% global = global
    x% name =  'vector?' ; if (present(name))   x% name   = name
    if (dace% npe == 1) x% global = .true.
    !
    ! allocate and set pointer components
    !
    allocate (x% s (x% n_s))
    do i = 1, x% n_s
      s => x% s(i)
      s% n  = info% b(i)% n
      s% pe = info% b(i)% pe
      if (s% pe == dace% pe .or. x% global) allocate (s% x (s% n))
    end do
    !
    ! maintain linked list
    !
#if defined (LINKED_LIST)
    if (x% alloc_l == 0) then
      if (associated (last_i)) then
        last_b% next => x
        x%      prev => last_b
        last_b       => x
      else
        first_b => x
        last_b  => x
      endif
    endif
#endif
  end subroutine construct_bvector
!==============================================================================
  !
  ! The deinitialisation routines (DESTRUCT) deallocates memory and
  ! resets flags to default values.
  !
!------------------------------------------------------------------------------
  subroutine destruct_dec_info (z)
  type (t_dec_info) ,intent(inout) :: z
    type (t_dec_info) :: empty_info
    integer :: i
    if (associated (z% b)) then
      do i=1,size(z% b)
        call destruct (z% b(i)% seed)
      end do
      deallocate (z% b)
    endif
    z = empty_info
  end subroutine destruct_dec_info
!------------------------------------------------------------------------------
  recursive subroutine destruct_matrix (a)
  type (t_matrix) ,intent(inout) :: a
    integer             :: i, j
    if (associated (a% b)) then
      do j = 1, size(a% b, 2)
        do i = 1, size(a% b, 1)
          call deallocate (a% b(i,j))
        end do
      end do
      deallocate (a% b)
    endif
    a% rinfo => NULL()
    a% cinfo => NULL()
    a% m     =  0
    a% n     =  0
    a% m_b   =  0
    a% n_b   =  0
    a% qual  =  -1
    !
    ! maintain linked list
    !
    if (associated (a% next)) then
      a% next% prev => a% prev
    else
      last_m        => a% prev
    endif
    if (associated (a% prev)) then
      a% prev% next => a% next
    else
      first_m       => a% next
    endif
    nullify (a% next)
    nullify (a% prev)
  end subroutine destruct_matrix
!------------------------------------------------------------------------------
  subroutine delete_matrix_block (x)
  type (t_matrix_block) ,intent(in) :: x
    target                        :: x
    type(t_matrix_block) ,pointer :: p
    if (x% alloc_l == call_level) then
      p => x
      call deallocate (p)
    endif
  end subroutine delete_matrix_block
!------------------------------------------------------------------------------
  recursive subroutine dealloc_matrix_block (b)
  type (t_matrix_block) ,intent (inout) :: b
    if (setinvalid) call set_invalid (b)
    if (associated (b% full  )) deallocate (b% full)
    if (associated (b% packed)) deallocate (b% packed)
    if (associated (b% ia    )) deallocate (b% ia)
    if (associated (b% ja    )) deallocate (b% ja)
    if (associated (b% nest  )) then
      call destruct (b% nest)
      deallocate    (b% nest)
    endif
    if (associated (B% perm)) then
       call deallocate (b% perm)
       deallocate (B% perm)
    end if
  end subroutine dealloc_matrix_block
!------------------------------------------------------------------------------
  subroutine destruct_perm_block (B)
    type (t_matrix_block), intent(inout) :: B
    if (associated (B% perm)) then
       call destruct (B% perm)
       deallocate (B% perm)
    end if
    nullify (B% perm)
  end subroutine destruct_perm_block
!------------------------------------------------------------------------------
  subroutine destruct_perm_ (perm)
    type (t_perm), intent(inout) :: perm
    perm% m      =  0
    perm% maxnz  =  0
    perm% jdiag  = -1
    perm% ngroup = -1
    call deallocate (perm)
  end subroutine destruct_perm_
!------------------------------------------------------------------------------
  subroutine dealloc_perm (perm)
    type (t_perm), intent(inout) :: perm
    if (associated (perm% iperm))   deallocate (perm% iperm)
    if (associated (perm% ia   ))   deallocate (perm% ia)
    if (associated (perm% nzgroup)) deallocate (perm% nzgroup)
    if (associated (perm% xgroup))  deallocate (perm% xgroup)
    nullify (perm% iperm, perm% ia, perm% xgroup, perm% nzgroup)
  end subroutine dealloc_perm
!------------------------------------------------------------------------------
  subroutine destruct_matrix_block (b)
  type (t_matrix_block) ,intent (inout) :: b
    call dealloc_matrix_block (b)
    call construct_matrix_block (b)
  end subroutine destruct_matrix_block
!------------------------------------------------------------------------------
  subroutine set_invalid_matrix_block (b)
  type (t_matrix_block) ,intent (inout) :: b
    if (associated (b% full  )) b% full   = -huge(1._mp)
    if (associated (b% packed)) b% packed = -huge(1._mp)
    if (associated (b% ia    )) b% ia     = -huge(1)
    if (associated (b% ja    )) b% ja     = -huge(1_ip)
    if (associated (b% nest  )) then
    endif
  end subroutine set_invalid_matrix_block
!------------------------------------------------------------------------------
  subroutine destruct_vector (x)
  type (t_vector) ,intent(inout) :: x
    type (t_vector_segm) ,pointer :: s
    integer                       :: i
    !
    ! deallocate pointer component
    !
    if (associated (x% s)) then
      do i = 1, size(x% s)
        s => x% s(i)
        if (associated (s% x)) deallocate (s% x)
      end do
      deallocate (x% s)
    endif
    !
    ! maintain linked list
    !
    if (x% alloc_l == 0) then
      if (associated (x% next)) then
        x% next% prev => x% prev
      else
        last_v        => x% prev
      endif
      if (associated (x% prev)) then
        x% prev% next => x% next
      else
        first_v       => x% next
      endif
    endif
    !
    ! reset components
    !
    x% n       =  0
    x% n_s     =  0
    x% alloc_l =  0
  end subroutine destruct_vector
!------------------------------------------------------------------------------
  subroutine destruct_vectors (x)
  type (t_vector) ,intent(inout) :: x (:)
  !---------------------------------------------
  ! Deallocate components of an array of vectors
  !---------------------------------------------
    integer :: i
    do i=1,size(x)
      call destruct (x(i))
    end do
  end subroutine destruct_vectors
!------------------------------------------------------------------------------
  subroutine destruct_vectors2 (x)
  type (t_vector) ,intent(inout) :: x (:,:)
  !---------------------------------------------
  ! Deallocate components of an array of vectors
  !---------------------------------------------
    integer :: i,j
    do j=1,size(x,2)
      do i=1,size(x,1)
        call destruct (x(i,j))
      end do
    end do
  end subroutine destruct_vectors2
!------------------------------------------------------------------------------
  subroutine destruct_ivector (x)
  type (t_ivector) ,intent(inout) :: x
    type (t_ivector_segm) ,pointer :: s
    integer                        :: i
    !
    ! deallocate pointer component
    !
    if (associated (x% s)) then
      do i = 1, size(x% s)
        s => x% s(i)
        if (associated (s% x)) deallocate (s% x)
      end do
      deallocate (x% s)
    endif
    !
    ! maintain linked list
    !
    if (x% alloc_l == 0) then
      if (associated (x% next)) then
        x% next% prev => x% prev
      else
        last_i        => x% prev
      endif
      if (associated (x% prev)) then
        x% prev% next => x% next
      else
        first_i       => x% next
      endif
    endif
    !
    ! reset components
    !
    x% n       =  0
    x% n_s     =  0
    x% alloc_l =  0
  end subroutine destruct_ivector
!------------------------------------------------------------------------------
  subroutine destruct_bvector (x)
  type (t_bvector) ,intent(inout) :: x
    type (t_bvector_segm) ,pointer :: s
    integer                        :: i
    !
    ! deallocate pointer component
    !
    if (associated (x% s)) then
      do i = 1, size(x% s)
        s => x% s(i)
        if (associated (s% x)) deallocate (s% x)
      end do
      deallocate (x% s)
    endif
    !
    ! maintain linked list
    !
    if (x% alloc_l == 0) then
      if (associated (x% next)) then
        x% next% prev => x% prev
      else
        last_b        => x% prev
      endif
      if (associated (x% prev)) then
        x% prev% next => x% next
      else
        first_b       => x% next
      endif
    endif
    !
    ! reset components
    !
    x% n       =  0
    x% n_s     =  0
    x% alloc_l =  0
  end subroutine destruct_bvector
!------------------------------------------------------------------------------
  !
  ! The delete routine is required to deallocate components of function
  ! return arguments.
  !
  subroutine delete_dec_vector (x)
#if defined (__GFORTRAN__) && (__GNUC__ * 100 + __GNUC_MINOR__) >= 407
  ! Workaround for GNU Fortran 4.7 (prerelease)
  type (t_vector)             :: x
#else
  type (t_vector) ,intent(in) :: x
#endif
    type (t_vector_segm) ,pointer :: s
    type (t_vector_segm) ,pointer :: z(:)
    integer                       :: i
    if (x% alloc_l == call_level) then
      if (associated (x% s)) then
        do i = 1, size(x% s)
          s => x% s(i)
          if (associated (s% x)) deallocate (s% x)
        end do
        z=> x% s
        deallocate (z)
      endif
      !
      ! maintain linked list
      !
!      if (associated (x% next)) then
!        x% next% prev => x% prev
!      else
!        last_v        => x% prev
!      endif
!      if (associated (x% prev)) then
!        x% prev% next => x% next
!      else
!        first_v       => x% next
!      endif
    endif
  end subroutine delete_dec_vector
!------------------------------------------------------------------------------
  !
  ! The delete routine is required to deallocate components of function
  ! return arguments.
  !
  subroutine delete_dec_integer (x)
  type (t_ivector) ,intent(in) :: x
    type (t_ivector_segm) ,pointer :: s
    type (t_ivector_segm) ,pointer :: z(:)
    integer                       :: i
    if (x% alloc_l == call_level) then
      if (associated (x% s)) then
        do i = 1, size(x% s)
          s => x% s(i)
          if (associated (s% x)) deallocate (s% x)
        end do
        z=> x% s
        deallocate (z)
      endif
      !
      ! maintain linked list
      !
!      if (associated (x% next)) then
!        x% next% prev => x% prev
!      else
!        last_v        => x% prev
!      endif
!      if (associated (x% prev)) then
!        x% prev% next => x% next
!      else
!        first_v       => x% next
!      endif
    endif
  end subroutine delete_dec_integer
!------------------------------------------------------------------------------
  !
  ! The delete routine is required to deallocate components of function
  ! return arguments.
  !
  subroutine delete_dec_boolean (x)
  type (t_bvector) ,intent(in) :: x
    type (t_bvector_segm) ,pointer :: s
    type (t_bvector_segm) ,pointer :: z(:)
    integer                       :: i
    if (x% alloc_l == call_level) then
      if (associated (x% s)) then
        do i = 1, size(x% s)
          s => x% s(i)
          if (associated (s% x)) deallocate (s% x)
        end do
        z=> x% s
        deallocate (z)
      endif
      !
      ! maintain linked list
      !
!      if (associated (x% next)) then
!        x% next% prev => x% prev
!      else
!        last_v        => x% prev
!      endif
!      if (associated (x% prev)) then
!        x% prev% next => x% next
!      else
!        first_v       => x% next
!      endif
    endif
  end subroutine delete_dec_boolean
!------------------------------------------------------------------------------
  !
  ! The delete routine is required to deallocate components of function
  ! return arguments.
  !
  subroutine delete_dec_matrix (x)
  type (t_matrix) ,intent(in) :: x
    target                   :: x
    type (t_matrix) ,pointer :: p
    if (x% alloc_l == call_level) then
      p => x
      call destruct (p)
    endif
  end subroutine delete_dec_matrix
!==============================================================================
  !
  ! The following routines perform the high level organization of the
  ! arithmetic operations +,-,* =, i.e. and transformations between
  ! matrix representations. They call the required communication
  ! routines and arithmetic routines on the block level.
  !
  ! Currently single processor mode is assumed.
  !
  ! For FULL and PACKED representation, blas matrix-vector multiplication
  ! routines are called. For CSR storage, hand written routines are
  ! currently used.
  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: global
!
! !DESCRIPTION:

  !
  ! store vector content on all processors
  !

! !INTERFACE:
!
  subroutine global_0 (x, tmp)
  type (t_vector)   ,intent(inout) :: x
  logical ,optional ,intent(in)    :: tmp
!
!EOP
    type (t_vector_segm) ,pointer :: s
    integer                       :: i
!
    if (.not. x% global) then
      x% global = .true.
      if(present(tmp)) x% global = .not. tmp
      do i = 1, x% n_s
        s => x% s(i)
        if (.not.associated(s% x)) allocate (s% x (s% n))
        if (s% n == 0)             cycle

#ifdef NO_MPI3                          /* fallback for MPI 2.x */
        call p_bcast  (s% x, s% pe)
#else                                   /* we do have MPI >=3.0 */
        call p_ibcast (s% x, s% pe)
        if (p_irequest > p_mrequest) call p_waitall ()
#endif

      end do
      call p_waitall ()
    endif
  end subroutine global_0
!------------------------------------------------------------------------------
  subroutine global_1 (x, tmp)
  type (t_vector)   ,intent(inout) :: x(:)
  logical ,optional ,intent(in)    :: tmp
    integer :: i
    do i=1,size(x)
      call global (x(i))
    end do
  end subroutine global_1
!------------------------------------------------------------------------------
  subroutine global_b (x, tmp)
  type (t_bvector)  ,intent(inout) :: x
  logical ,optional ,intent(in)    :: tmp

    type (t_bvector_segm) ,pointer :: s
    integer                        :: i

    if (.not. x% global) then
      x% global = .true.
      if(present(tmp)) x% global = .not. tmp
      do i = 1, x% n_s
        s => x% s(i)
        if (.not.associated(s% x)) allocate (s% x (s% n))
        if (s% n == 0)             cycle
        call p_bcast (s% x, s% pe)
      end do
    endif
  end subroutine global_b
!------------------------------------------------------------------------------
  subroutine gather_r (x, dest, isgm)
  type (t_vector) ,intent(in)           :: x    ! vector
  integer         ,intent(in) ,optional :: dest ! PE to gather segments on
  integer         ,intent(in) ,optional :: isgm ! only this segment
  !------------------------------------------
  ! gather real vector content on 1 processor
  !------------------------------------------
    type (t_vector_segm) ,pointer :: s
    integer                       :: i
    integer                       :: is
    is = 0; if(present(isgm)) is = isgm
    if (.not. x% global) then
      !-------------------
      ! loop over segments
      !-------------------
      do i = 1, x% n_s
        if (present(isgm).and.is/=i) cycle
        s => x% s(i)
        !-----------------
        ! gather on one PE
        !-----------------
        if (present (dest)) then
          if (s% pe /= dest) then   ! no action if segment is on destination PE
            if (dace% pe == dest) then  ! gather on destination PE
              if (associated(s% x)) then
                if (size(s% x) /= s% n) deallocate (s% x)
              endif
              if (.not.associated(s% x)) allocate (s% x (s% n))
              call p_recv (s% x, s% pe, p_tag=1)
            else if (s% pe == dace% pe) then
              call p_send (s% x, dest, p_tag=1)
            endif
          endif
        !------------------
        ! gather on all PEs
        !------------------
        else
          if (.not.associated(s% x)) allocate (s% x (s% n))
          if (s% n == 0)             cycle
          call p_bcast (s% x, s% pe)
        endif
      end do
    endif
  end subroutine gather_r
!------------------------------------------------------------------------------
  subroutine scatter_r (x, src)
  type (t_vector) ,intent(in) :: x    ! vector
  integer         ,intent(in) :: src  ! PE to scatter segments from
  !----------------------------
  ! scatter real vector content
  !----------------------------
    type (t_vector_segm) ,pointer :: s
    integer                       :: i
    !-------------------
    ! loop over segments
    !-------------------
    do i = 1, x% n_s
      s => x% s(i)
      if (.not. x% global) then
        !------------------
        ! scatter on one PE
        !------------------
        if (s% pe /= src) then   ! no action if segment is on source PE
          if (dace% pe == src) then  ! send from source PE
            call p_send (s% x, s% pe, p_tag=1)
          else if (s% pe == dace% pe) then
            if (associated (s% x)) then
              if (size(s% x) /= s% n) deallocate (s% x)
            endif
            if (.not.associated(s% x)) allocate (s% x (s% n))
            call p_recv (s% x, src, p_tag=1)
          endif
        endif
      else
        !-------------------
        ! scatter to all PEs
        !-------------------
        if (associated (s% x)) then
          if (size(s% x) /= s% n) deallocate (s% x)
        endif
        if (.not.associated (s% x)) allocate (s% x (s% n))
        if (s% n == 0)              cycle
        call p_bcast (s% x, src)
      endif
    end do
  end subroutine scatter_r
!------------------------------------------------------------------------------
  subroutine gather_re (x, dest, isgm)
  type (t_vector) ,intent(in)           :: x(:) ! vector
  integer         ,intent(in) ,optional :: dest ! PE to gather segments on
  integer         ,intent(in) ,optional :: isgm ! only this segment
  !------------------------------------------
  ! gather real vector content on 1 processor
  !------------------------------------------
    integer                  :: i

    !-------------------
    ! loop over elements
    !-------------------
    do i = 1, size(x)
      call gather (x(i), dest, isgm)
    end do

  end subroutine gather_re
!------------------------------------------------------------------------------
  subroutine scatter_re (x, src)
  type (t_vector) ,intent(in)           :: x(:) ! vector
  integer         ,intent(in) ,optional :: src  ! PE to scatter segments from
  !------------------------------------------
  ! scatter real vector content on 1 processor
  !------------------------------------------
    integer                  :: i

    !-------------------
    ! loop over elements
    !-------------------
    do i = 1, size(x)
      call scatter (x(i), src)
    end do

  end subroutine scatter_re
!------------------------------------------------------------------------------
  subroutine gather_i (x, dest, isgm)
  type (t_ivector) ,intent(in)           :: x    ! vector
  integer          ,intent(in) ,optional :: dest ! PE to gather segments on
  integer          ,intent(in) ,optional :: isgm ! only this segment
  !---------------------------------------------
  ! gather integer vector content on 1 processor
  !---------------------------------------------
    type (t_ivector_segm) ,pointer :: s
    integer                        :: i
    integer                        :: is
    is = 0; if(present(isgm)) is = isgm
    if (.not. x% global) then
      !-------------------
      ! loop over segments
      !-------------------
      do i = 1, x% n_s
        if (present(isgm).and.is/=i) cycle
        s => x% s(i)
        !-----------------
        ! gather on one PE
        !-----------------
        if (present (dest)) then
          if (s% pe /= dest) then   ! no action if segment is on destination PE
            if (dace% pe == dest) then  ! gather on destination PE
              if (associated(s% x)) then
                if (size(s% x) /= s% n) deallocate (s% x)
              endif
              if (.not.associated(s% x)) allocate (s% x (s% n))
              call p_recv (s% x, s% pe, p_tag=1)
            else if (s% pe == dace% pe) then
              call p_send (s% x, dest, p_tag=1)
            endif
          endif
        !------------------
        ! gather on all PEs
        !------------------
        else
          if (.not.associated(s% x)) allocate (s% x (s% n))
          if (s% n == 0)             cycle
          call p_bcast (s% x, s% pe)
        endif
      end do
    endif
  end subroutine gather_i
!------------------------------------------------------------------------------
  subroutine gather_b (x, dest, isgm)
  type (t_bvector) ,intent(in)           :: x    ! vector
  integer          ,intent(in) ,optional :: dest ! PE to gather segments on
  integer          ,intent(in) ,optional :: isgm ! only this segment
  !---------------------------------------------
  ! gather boolean vector content on 1 processor
  !---------------------------------------------
    type (t_bvector_segm) ,pointer :: s
    integer                        :: i
    integer                        :: is
    is = 0; if(present(isgm)) is = isgm
    if (.not. x% global) then
      !-------------------
      ! loop over segments
      !-------------------
      do i = 1, x% n_s
        if (present(isgm).and.is/=i) cycle
        s => x% s(i)
        !-----------------
        ! gather on one PE
        !-----------------
        if (present (dest)) then
          if (s% pe /= dest) then   ! no action if segment is on destination PE
            if (dace% pe == dest) then  ! gather on destination PE
              if (associated(s% x)) then
                if (size(s% x) /= s% n) deallocate (s% x)
              endif
              if (.not.associated(s% x)) allocate (s% x (s% n))
              call p_recv (s% x, s% pe, p_tag=1)
            else if (s% pe == dace% pe) then
              call p_send (s% x, dest, p_tag=1)
            endif
          endif
        !------------------
        ! gather on all PEs
        !------------------
        else
          if (.not.associated(s% x)) allocate (s% x (s% n))
          if (s% n == 0)             cycle
          call p_bcast (s% x, s% pe)
        endif
      end do
    endif
  end subroutine gather_b
!------------------------------------------------------------------------------
  subroutine release_vector_mem_1 (x)
  type (t_vector)   ,intent(in) :: x(:)
    integer :: i
    do i=1,size(x)
      call release_vector_mem (x(i))
    end do
  end subroutine release_vector_mem_1
!------------------------------------------------------------------------------
  subroutine release_vector_mem (x)
  type (t_vector) ,intent(in)     :: x
    type (t_vector_segm) ,pointer :: s
    integer                       :: i
    if (.not. x% global) then
      do i = 1, x% n_s
        s => x% s(i)
        if (s% pe /= dace% pe .and. associated(s% x)) deallocate (s% x)
      end do
    end if
  end subroutine release_vector_mem
!------------------------------------------------------------------------------
  subroutine release_ivector_mem (x)
  type (t_ivector) ,intent(in)     :: x
    type (t_ivector_segm) ,pointer :: s
    integer                        :: i
    if (.not. x% global) then
      do i = 1, x% n_s
        s => x% s(i)
        if (s% pe /= dace% pe .and. associated(s% x)) deallocate (s% x)
      end do
    end if
  end subroutine release_ivector_mem
!------------------------------------------------------------------------------
  subroutine local_1 (x)
  type (t_vector)   ,intent(inout) :: x(:)
    integer :: i
    do i=1,size(x)
      call local_0 (x(i))
    end do
  end subroutine local_1
!------------------------------------------------------------------------------
  subroutine local_0 (x)
  type (t_vector) ,intent(inout)  :: x
    type (t_vector_segm) ,pointer :: s
    integer                       :: i
    x% global = .false.
    do i = 1, x% n_s
      s => x% s(i)
      if (s% pe /= dace% pe .and. associated(s% x)) deallocate (s% x)
    end do
  end subroutine local_0
!------------------------------------------------------------------------------
  subroutine release_bvector_mem (x)
  type (t_bvector) ,intent(in)     :: x
    type (t_bvector_segm) ,pointer :: s
    integer                        :: i
    if (.not. x% global) then
      do i = 1, x% n_s
        s => x% s(i)
        if (s% pe /= dace% pe .and. associated(s% x)) deallocate (s% x)
      end do
    end if
  end subroutine release_bvector_mem
!------------------------------------------------------------------------------
  subroutine release_mblock_mem (a)
  type (t_matrix_block) ,intent(inout) :: a
  !--------------------------------------------
  ! release unused memory within a matrix block
  !--------------------------------------------
    integer(ip) ,pointer :: itmp (:)
    real   (mp) ,pointer ::  tmp (:)
    integer              :: n
    !-------------------------------------
    ! deallocate unused pointer components
    !-------------------------------------
    select case (a% repr)
    case (ZERO,IDENT,MIRROR,NOT_INIT)
      if (associated(a% full))   deallocate (a% full)
      if (associated(a% packed)) deallocate (a% packed)
      if (associated(a% ia))     deallocate (a% ia)
      if (associated(a% ja))     deallocate (a% ja)
      if (associated(a% nest))   deallocate (a% nest)
    case (DIAGONAL,PACKED,CHOLES)
      if (associated(a% full))   deallocate (a% full)
      if (associated(a% ia))     deallocate (a% ia)
      if (associated(a% ja))     deallocate (a% ja)
      if (associated(a% nest))   deallocate (a% nest)
    case (FULL)
      if (associated(a% packed)) deallocate (a% packed)
      if (associated(a% ia))     deallocate (a% ia)
      if (associated(a% ja))     deallocate (a% ja)
      if (associated(a% nest))   deallocate (a% nest)
    case (LU)
      if (associated(a% packed)) deallocate (a% packed)
      if (associated(a% ja))     deallocate (a% ja)
      if (associated(a% nest))   deallocate (a% nest)
    case (NEST)
      if (associated(a% full))   deallocate (a% full)
      if (associated(a% packed)) deallocate (a% packed)
      if (associated(a% ia))     deallocate (a% ia)
      if (associated(a% ja))     deallocate (a% ja)
    case (CSR,CSC)
      if (associated(a% full))   deallocate (a% full)
      if (associated(a% nest))   deallocate (a% nest)
    case default
      call finish ('release_mblock_mem','invalid repr '//srep(a% repr))
    end select

    !--------------------------------------
    ! shrink components to appropriate size
    !--------------------------------------
    select case (a% repr)
    case (CSR,CSC)
      if (a% repr == CSR) n = a% m + 1
      if (a% repr == CSC) n = a% n + 1

      if (size(a% packed) < a% nonzero) &
        call finish('release_mblock_mem','packed is too small')
      if (size(a% ja    ) < a% nonzero) &
        call finish('release_mblock_mem','ja is too small')
      if (size(a% ia    ) < n) &
        call finish('release_mblock_mem','ia is too small')

      if (size(a% packed) > a% nonzero) then
        allocate (tmp (a% nonzero))
        tmp = a% packed(:a% nonzero)
        deallocate (a% packed)
        a% packed => tmp
      endif

      if (size(a% ja) > a% nonzero) then
        allocate (itmp (a% nonzero))
        itmp = a% ja(:a% nonzero)
        deallocate (a% ja)
        a% ja => itmp
      endif

      if (size(a% ia) > n) then
        allocate (itmp (n))
        itmp = a% ia(:n)
        deallocate (a% ia)
        a% ia => itmp
      endif
    end select

  end subroutine release_mblock_mem
!------------------------------------------------------------------------------
  function matrix_times_vector (a, x) result (y)
  !----------------
  ! matrix * vector
  !----------------
  type (t_vector)                     :: y
  type (t_matrix) ,intent(in)         :: a
  type (t_vector) ,intent(in) ,target :: x
    integer                  :: i, j
    type (t_vector) ,pointer :: z
    call enter_function
    call construct (y, a% rinfo, name='A*x', global=.not.a%rowwise)
    !---------------------------------------------------
    ! gather right hand side X if not globally allocated
    !---------------------------------------------------
    if (x% global .or. a% qual == BDIAG) then
      z => x
    else
      allocate (z)
      call construct (z, x% info, global = .true.)
      z = x
    endif
    !-----------------------
    ! matrix-vector multiply
    !-----------------------
! OpenMP disabled here since some LAPACK/BLAS libraries may not be thread-safe
!!$omp parallel do private(i,j) schedule(dynamic)
    do i = 1, a% m_b
      if (y% s(i)% pe == dace% pe .or. .not. a% rowwise) then
        if (a% rowwise .and. any(a% b(i,:)% pe /= dace% pe)) then
          write(0,*)dace% pe,i,':',a% b(i,:)% pe
          call finish ('matrix_times_vector','pe mismatch: '//a% name)
        endif
        y% s(i)% x = 0._wp
        do j = 1, a% n_b
          if (a% qual == BDIAG  .and. i/=j) cycle
          if (a% qual == BUPPER .and. i<=j) cycle
          if (a% qual == BLOWER .and. i>=j) cycle
          if (a% b(i,j)% pe /= dace% pe)        cycle
          select case (a% b(i,j)% repr)
          case (ZERO)
          case default
            call matrix_times_vector_block(y% s(i)% x, a% b(i,j), &
                                           z% s(j)% x, .false.)
          case (MIRROR)
            call matrix_times_vector_block(y% s(i)% x, a% b(j,i), &
                                           z% s(j)% x, .true.)
          end select
        end do
      end if
    end do
!!$omp end parallel do
    !--------
    ! cleanup
    !--------
    if (.not. (x% global .or. a% qual == BDIAG)) then
      call destruct (z)
      deallocate (z)
    endif
    if (.not.a%rowwise) then
      call p_sum_vector (y)
      y% global = .false.
      call release_mem (y)
    endif
    call delete_storage (a)
    call delete_storage (x)
    call leave_function
  end function matrix_times_vector
!------------------------------------------------------------------------------
  function vector_times_matrix (x, a) result (y)
  !-----------------------------------------------------------
  ! vector * matrix (multiply vector by the transposed matrix)
  !-----------------------------------------------------------
  type (t_vector)                     :: y
  type (t_vector) ,intent(in) ,target :: x
  type (t_matrix) ,intent(in)         :: a
    integer                  :: i, j
    type (t_vector) ,pointer :: z
    !-----------------
    ! check dimensions
    !-----------------
    if (x% n /= a% m) then
      write(0,*)   'vector_times_matrix',x%name,x%n,a%name,a%m
      call finish ('vector_times_matrix','dimensions do not fit')
    endif
    if (x% n_s /= a% m_b) then
      write(0,*)   'vector_times_matrix',x%name,x%n_s,a%name,a%m_b
      call finish ('vector_times_matrix','number of blocks does not fit')
    endif
    !------------------------
    ! construct result vector
    !------------------------
    call enter_function
    call construct (y, a% cinfo, name='x*A', global=.true.)
    !---------------------------------------------------
    ! gather left hand side X if not globally allocated
    !---------------------------------------------------
    if (x% global .or. a% qual == BDIAG .or. a% rowwise) then
      z => x
    else
      allocate (z)
      call construct (z, x% info, global = .true.)
      z = x
    endif
    !-----------------------
    ! matrix-vector multiply
    !-----------------------
!$omp parallel do private(i,j) schedule(dynamic)
    do j = 1, a% n_b
      y% s(j)% x = 0._wp
      do i = 1, a% m_b
        if (associated (z% s(i)% x)) then
          if (a% qual == BDIAG  .and. i/=j) cycle
          if (a% qual == BUPPER .and. i<=j) cycle
          if (a% qual == BLOWER .and. i>=j) cycle
          select case (a% b(i,j)% repr)
          case (ZERO)
          case default
            call matrix_times_vector_block(y% s(j)% x, a% b(i,j), &
                                           z% s(i)% x, .true.)
          case (MIRROR)
            call matrix_times_vector_block(y% s(j)% x, a% b(j,i), &
                                           z% s(i)% x, .false.)
          end select
        end if
      end do
    end do
!$omp end parallel do
    !--------
    ! cleanup
    !--------
    if (.not. (x% global .or. a% qual == BDIAG .or. a% rowwise)) then
      call destruct (z)
      deallocate (z)
    endif
    call p_sum_vector (y)
    y% global = .false.
    call release_mem (y)
    call delete_storage (a)
    call delete_storage (x)
    call leave_function
  end function vector_times_matrix
!------------------------------------------------------------------------------
  function matrix_times_array (a, x) result (y)
  !---------------
  ! matrix * array
  !---------------
  type (t_matrix) ,intent(in) :: a
  real(wp)        ,intent(in) :: x (:)
  real(wp)                    :: y (a% m)
    integer :: ib,jb,i,j,n,m
    !-------
    ! checks
    !-------
    if (size(x) /= a% n)  then
      write (0,*) 'matrix_times_array: size(x) /= a% n)'
      write (0,*) '  matrix:',a%name,', m,n,m_b,n_b=', a%m,a%n,a%m_b,a%n_b
      write (0,*) '  size(x) =',size(x)
      call finish ('matrix_times_array','size(x) /= a% n')
    endif
    if (a% qual /= BDIAG) then
      write (0,*) 'matrix_times_array (',a%name,'): qual=',a% qual
      write (0,*) 'matrix_times_array (',a%name,'): non BDIAG not implemented'
      call finish('matrix_times_array',            'non BDIAG not implemented')
    endif
    !-----------------------
    ! matrix-vector multiply
    !-----------------------
    y = 0._wp
    i = 0
    j = 0
    do ib= 1, a% m_b
      jb = ib
      m = a% b(ib,jb)% m
      n = a% b(ib,jb)% n
      select case (a% b(ib,jb)% repr)
      case (ZERO)
      case default
        call matrix_times_vector_block (y(i+1:i+m), a% b(ib,jb),&
                                        x(j+1:j+n), .false.     )
      case (MIRROR)
        call matrix_times_vector_block (y(i+1:i+m), a% b(jb,ib),&
                                        x(j+1:j+n), .true.)
      end select
      i = i + m
      j = j + n
    end do
    call enter_function
    call delete_storage (a)
    call leave_function
  end function matrix_times_array
!------------------------------------------------------------------------------
  function array_times_matrix (x, a) result (y)
  !-----------------------------------------------------
  ! array * matrix (multyply array by transposed matrix)
  !-----------------------------------------------------
  real(wp)        ,intent(in) :: x (:)
  type (t_matrix) ,intent(in) :: a
  real(wp)                    :: y (a% n)
    integer :: ib,jb,i,j,n,m
    !-------
    ! checks
    !-------
    if (size(x) /= a% m)  call finish ('array_times_matrix','size(x) /= a% m')
    if (a% qual /= BDIAG) then
      write (0,*) 'array_times_matrix (',a%name,'): qual=',a% qual
      write (0,*) 'array_times_matrix (',a%name,'): non BDIAG not implemented'
      call finish('array_times_matrix',            'non BDIAG not implemented')
    endif
    !-----------------------
    ! matrix-vector multiply
    !-----------------------
    y = 0._wp
    i = 0
    j = 0
    do ib= 1, a% m_b
      jb = ib
      m = a% b(ib,jb)% m
      n = a% b(ib,jb)% n
      select case (a% b(ib,jb)% repr)
      case (ZERO)
      case default
        call matrix_times_vector_block (y(j+1:j+n), a% b(ib,jb),&
                                        x(i+1:i+m), .true.     )
      case (MIRROR)
        call matrix_times_vector_block (y(j+1:j+n), a% b(jb,ib),&
                                        x(i+1:i+m), .false.)
      end select
      i = i + m
      j = j + n
    end do
    call enter_function
    call delete_storage (a)
    call leave_function
  end function array_times_matrix
!------------------------------------------------------------------------------
  recursive function scalar_times_matrix (a, X) result (y)
  !--------------
  ! real * matrix
  !--------------
  type (t_matrix)             :: Y
  real(wp)        ,intent(in) :: a
  type (t_matrix) ,intent(in) :: X
    integer                        :: i, j
    type (t_matrix_block) ,pointer :: xb, yb
    call enter_function
    call construct (y, ri=X% rinfo, ci=X% cinfo, qual=X%qual, name='a*X')
    do i = 1, X% m_b
      do j = 1, X% n_b
        xb => X% b(i,j)
        yb => Y% b(i,j)
        yb% repr = xb% repr
        yb% pe   = xb% pe
        yb% m    = xb% m
        yb% n    = xb% n
        yb% tri  = xb% tri
        if (yb% pe == dace% pe) then
          call allocate_block (yb, xb% repr, xb% m, xb% n, dealloc=.false.)
          if (associated(xb% full))   yb% full   = a * xb% full
          if (associated(xb% packed)) yb% packed = a * xb% packed
          if (associated(xb% nest))   yb% nest   = a * xb% nest
        endif
      end do
    end do
    call delete_storage (X)
    call leave_function
  end function scalar_times_matrix
!------------------------------------------------------------------------------
  recursive function matrix_times_matrix (a, b) result (y)
  !----------------
  ! matrix * matrix
  !----------------
  type (t_matrix)             :: y
  type (t_matrix) ,intent(in) :: a
  type (t_matrix) ,intent(in) :: b
  !-------------------------------------------------
  ! matrix multiplication: Yik = Aij * Bjk .
  ! currently we allow blockdiagonal matrices B only
  ! (so there is no need for communication)
  ! therefore: Yik = Aii * Bik
  !-------------------------------------------------
    integer                        :: i,j
    integer                        :: m, n
    integer                        :: ns
    type (t_matrix_block)          :: ab
    type (t_matrix_block)          :: bb
    type (t_matrix_block) ,pointer :: aii
    type (t_matrix_block) ,pointer :: bij
    type (t_matrix_block) ,pointer :: yij
    !-------------------
    ! consistency checks
    !-------------------
    if (a% qual /= BDIAG) call finish('matrix_times_matrix','A is not b.diag.')
    if (a% n    /= b% m)  call finish('matrix_times_matrix','a% n /= b% m')
    if (a% n_b  /= b% m_b)call finish('matrix_times_matrix','a% n_b /= b% m_b')
    !---------------
    ! initialisation
    !---------------
    call enter_function
    call construct (y, ri=a% rinfo, ci=b% cinfo, qual=b%qual, name='A*B')
    !------------------------
    ! loop over matrix blocks
    !------------------------
    do  i= 1, a% m_b
    aii => a% b(i,i)
    do j= 1, a% n_b
      bij => b% b(i,j)
      yij => y% b(i,j)
      m = aii% m
      n = bij% n
      if (aii% pe == dace% pe .and. bij% pe == dace% pe) then
        !--------------------------------------------
        ! both matrix elements are located on this PE
        !--------------------------------------------
        yij% alloc_l = y% alloc_l
        !-------------------------------
        ! treat some simple cases first:
        ! both matrix elements are ZERO
        !-------------------------------
        if (aii% repr == ZERO .or. bij% repr == ZERO) then
          yij% repr = ZERO
        !---------------------------
        ! one factor is the identity
        !---------------------------
        else if (aii% repr == IDENT) then
          yij = bij
        else if (bij% repr == IDENT) then
          yij = aii
        !---------------------------------------------
        ! matrix elements have the same representation
        !---------------------------------------------
        else if (aii% repr == bij% repr) then
          !-----------------------
          ! check for correct size
          !-----------------------
          if (aii% n /= bij% m) then
            write(6,'(i3,1x,a/2i5,3(1x,a8),3i5,l2,4i5)')               &
              dace% pe,'matrix*matrix: i,j, reps, pes, a:mn,b:mn',i ,j,&
              srep(yij% repr), srep(aii% repr),                        &
              srep(bij% repr), yij% pe,                                &
              aii% pe, bij% pe, dace% pe==bij% pe,                     &
              aii% m, aii% n, bij% m, bij% n
            call finish('matrix_times_matrix','aii% n /= bij% m')
          endif

          select case (aii% repr)
          case (ZERO)
          case (IDENT)
          case (DIAGONAL)
            call allocate_block(yij, DIAGONAL, m, n)
            yij% packed = aii% packed * bij% packed
          case (FULL)
            call allocate_block(yij, FULL, m, n)
            yij% full = matmul (aii% full, bij% full)
          case (NEST)
            yij% nest = aii% nest * bij% nest
          case (PACKED)
            call allocate_block(yij, FULL, m, n)
            call allocate_block(ab, aii% repr, aii% m,aii% n)
            call allocate_block(bb, bij% repr, bij% m,bij% n)
            ab% packed = aii% packed
            bb% packed = bij% packed
            ab% tri    = aii% tri
            bb% tri    = bij% tri
            call unpack_block (ab)
            call unpack_block (bb)
            yij% full = matmul (ab% full, bb% full)
            call deallocate (ab)
            call deallocate (bb)
          case default
            write(6,'(i3,1x,a/2i5,3(1x,a8),3i5,l2)')        &
              dace% pe,'matrix*matrix: i,j, reps, pes',i ,j,&
              srep(yij% repr), srep(aii% repr),             &
              srep(bij% repr), yij% pe,                     &
              aii% pe, bij% pe, dace% pe==bij% pe
            call finish('matrix_times_matrix',              &
                        'not implemented: '//srep(aii% repr))
          end select
        !---------------------------------------------
        ! operation on matrices in full representation
        !---------------------------------------------
        else if ((aii% repr == FULL   &
             .or. aii% repr == PACKED &
             .or. aii% repr == CSR    &
             .or. aii% repr == CSC)   &
            .and.(bij% repr == FULL   &
             .or. bij% repr == PACKED &
             .or. bij% repr == CSR    &
             .or. bij% repr == CSC)   ) then
          !-----------------------
          ! check for correct size
          !-----------------------
          if (aii% n /= bij% m) then
            write(6,'(i3,1x,a/2i5,3(1x,a8),3i5,l2,4i5)')               &
              dace% pe,'matrix*matrix: i,j, reps, pes, a:mn,b:mn',i ,j,&
              srep(yij% repr), srep(aii% repr),                        &
              srep(bij% repr), yij% pe,                                &
              aii% pe, bij% pe, dace% pe==bij% pe,                     &
              aii% m, aii% n, bij% m, bij% n
            call finish('matrix_times_matrix','aii% n /= bij% m')
          endif

          call allocate_block(yij, FULL, m, n, call_l=y%alloc_l)
          !-----------------------------
          ! get full representation of A
          !-----------------------------
          if (aii% repr == FULL) then
            ab% full => aii% FULL
          else
            ns = size (aii% packed)
            call allocate_block(ab, aii% repr, m,aii% n, ns = ns)
            ab% packed = aii% packed
            ab% tri    = aii% tri
            if (associated (aii% ia)) ab% ia = aii% ia
            if (associated (aii% ja)) ab% ja = aii% ja
            call unpack_block (ab)
          endif
          !-----------------------------
          ! get full representation of B
          !-----------------------------
          if (bij% repr == FULL) then
            bb% full => bij% FULL
          else
            ns = size (bij% packed)
            call allocate_block(bb, bij% repr, bij% m, n, ns=ns)
            bb% packed = bij% packed
            bb% tri    = bij% tri
            if (associated (bij% ia)) bb% ia = bij% ia
            if (associated (bij% ja)) bb% ja = bij% ja
            call unpack_block (bb)
          endif
          !-------
          ! matmul
          !-------
          yij% full = matmul (ab% full, bb% full)
          !-----------------------
          ! deallocate temporaries
          !-----------------------
          if (aii% repr /= FULL) call deallocate (ab)
          if (bij% repr /= FULL) call deallocate (bb)
        else
          write(0,'(i3,1x,a/2i5,3(1x,a8),3i5,l2)')         &
            dace% pe,'matrix*matrix: i,j, reps, pes',i ,j, &
            srep(yij% repr), srep(aii% repr),              &
            srep(bij% repr), yij% pe,                      &
            aii% pe, bij% pe, dace% pe==bij% pe
          call finish('matrix_times_matrix','a% repr /= b% repr')
        endif
      elseif (aii% pe /= dace% pe .and. bij% pe /= dace% pe) then
        !---------------------------------------------------
        ! none of the matrix elements is  located on this PE
        !---------------------------------------------------
      else
        call finish('matrix_times_matrix',&
                    'processor mismatch: '//a%name//b%name)
      endif
    end do
    end do
    !--------
    ! cleanup
    !--------
    call delete_storage (a)
    call delete_storage (b)
    call leave_function
  end function matrix_times_matrix
!------------------------------------------------------------------------------
  function scalar_times_vector (a, x) result (y)
  !--------------
  ! real * vector
  !--------------
  type (t_vector)             :: y
  real(wp)        ,intent(in) :: a
  type (t_vector) ,intent(in) :: x
    integer :: i
    call enter_function
    call construct (y, x% info, 's*v', x% global)
    do i = 1, x% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = a * x% s(i)% x
    end do
    call delete_storage (x)
    call leave_function
  end function scalar_times_vector
!------------------------------------------------------------------------------
  subroutine assign_vectors (y, x)
  !----------------------------
  ! vector(:) = vector(:)
  ! no check for overlap so far
  !----------------------------
  type (t_vector) ,intent(inout) :: y (:)
  type (t_vector) ,intent(in)    :: x (:)
    integer :: i
    if (size(x) /= size(y)) call finish ('assign_vectors','size mismatch')
    do i = 1, size (y)
      y(i) = x(i)
    end do
  end subroutine assign_vectors
!------------------------------------------------------------------------------
  subroutine assign_vector (y, x)
  !----------------
  ! vector = vector
  !----------------
  type (t_vector) ,intent(inout) :: y
  type (t_vector) ,intent(in)    :: x
    integer                   :: i
    type(t_dec_info) ,pointer :: info
    !----------------------------------
    ! construct destination if not done
    !----------------------------------
    if (.not.associated (y% info)) y% info => x% info
    info => y% info
    if (.not.associated (y% s   )) call construct (y, info)
    !-------------
    ! copy content
    !-------------
    call enter_function
    do i = 1, y% n_s
      if (associated(y% s(i)% x) .and. associated(x% s(i)% x)) &
        y% s(i)% x = x% s(i)% x
      if (y% global .and..not. x% global) &
        call p_bcast (y% s(i)% x, y% s(i)% pe)
    end do
    call delete_storage (x)
    call leave_function
  end subroutine assign_vector
!------------------------------------------------------------------------------
  subroutine assign_bvector (y, x)
  !----------------
  ! vector = vector
  !----------------
  type (t_bvector) ,intent(inout) :: y
  type (t_bvector) ,intent(in)    :: x
    integer :: i
    call enter_function
    do i = 1, y% n_s
      if (associated(y% s(i)% x) .and. associated(x% s(i)% x)) &
        y% s(i)% x = x% s(i)% x
      if (y% global .and..not. x% global) &
        call p_bcast (y% s(i)% x, y% s(i)% pe)
    end do
    call delete_storage (x)
    call leave_function
  end subroutine assign_bvector
!------------------------------------------------------------------------------
  subroutine assign_scalar_to_bvector (y, x)
  !------------------
  ! bvector = logical
  !------------------
  type (t_bvector) ,intent(inout) :: y
  logical          ,intent(in)    :: x
    integer :: i
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = x
    end do
  end subroutine assign_scalar_to_bvector
!------------------------------------------------------------------------------
  subroutine assign_ivector (y, x)
  !----------------
  ! vector = vector
  !----------------
  type (t_ivector) ,intent(inout) :: y
  type (t_ivector) ,intent(in)    :: x
    integer :: i
    call enter_function
    do i = 1, y% n_s
      if (associated(y% s(i)% x) .and. associated(x% s(i)% x)) &
        y% s(i)% x = x% s(i)% x
      if (y% global .and..not. x% global) &
        call p_bcast (y% s(i)% x, y% s(i)% pe)
    end do
    call delete_storage (x)
    call leave_function
  end subroutine assign_ivector
!------------------------------------------------------------------------------
  subroutine assign_vector_to_dparray (a, x)
  real(dp)        ,intent(out) :: a (:)
  type (t_vector) ,intent(in)  :: x
    integer                       :: i, i0, n
    type (t_vector_segm) ,pointer :: s
    call enter_function
    i0 = 0
    do i = 1, x% n_s
      s => x% s(i)
      if (.not. associated (s% x)) &
        call finish('assign_vector_to_array','segment is not associated')
      n = s% n
      a (i0+1:i0+n) = s% x
      i0 = i0 + n
    end do
    call delete_storage (x)
    call leave_function
  end  subroutine assign_vector_to_dparray
!------------------------------------------------------------------------------
  subroutine assign_ivector_to_i4array (a, x)
  integer(i4)      ,intent(out) :: a (:)
  type (t_ivector) ,intent(in)  :: x
    integer                        :: i, i0, n
    type (t_ivector_segm) ,pointer :: s
    call enter_function
    i0 = 0
    do i = 1, x% n_s
      s => x% s(i)
      if (.not. associated (s% x)) &
        call finish('assign_vector_to_array','segment is not associated')
      n = s% n
      a (i0+1:i0+n) = s% x
      i0 = i0 + n
    end do
    call delete_storage (x)
    call leave_function
  end  subroutine assign_ivector_to_i4array
!------------------------------------------------------------------------------
  subroutine assign_vector_to_sparray (a, x)
  real(sp)        ,intent(out) :: a (:)
  type (t_vector) ,intent(in)  :: x
    integer                       :: i, i0, n
    type (t_vector_segm) ,pointer :: s
    call enter_function
    i0 = 0
    do i = 1, x% n_s
      s => x% s(i)
      if (.not. associated (s% x)) &
        call finish('assign_vector_to_array','segment is not associated')
      n = s% n
      a (i0+1:i0+n) = s% x
      i0 = i0 + n
    end do
    call delete_storage (x)
    call leave_function
  end  subroutine assign_vector_to_sparray
!------------------------------------------------------------------------------
  subroutine assign_scalar_to_vectors (y, x)
  !-----------------
  ! vector(:) = real
  !-----------------
  type (t_vector) ,intent(inout) :: y(:)
  real(wp)        ,intent(in)    :: x
    integer :: i, k
    do k = 1, size(y)
      do i = 1, y(k)% n_s
        if (associated (y(k)% s(i)% x)) y(k)% s(i)% x = x
      end do
    end do
  end subroutine assign_scalar_to_vectors
!------------------------------------------------------------------------------
  subroutine assign_scalar_to_vector (y, x)
  !--------------
  ! vector = real
  !--------------
  type (t_vector) ,intent(inout) :: y
  real(wp)        ,intent(in)    :: x
    integer :: i
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = x
    end do
  end subroutine assign_scalar_to_vector
!------------------------------------------------------------------------------
  subroutine assign_int_to_ivector (y, x)
  !------------------
  ! ivector = integer
  !------------------
  type (t_ivector) ,intent(inout) :: y
  integer          ,intent(in)    :: x
    integer :: i
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = x
    end do
  end subroutine assign_int_to_ivector
!------------------------------------------------------------------------------
  function diag_matrix(x) result (y)
  type (t_vector)                       :: y
  type (t_matrix) ,intent(in)           :: x
    integer                             :: i
    if (x% m /= x% n) call finish ('diag_matrix','matrix is not square')
    call enter_function
    call construct (y, x% rinfo)
    do i = 1, y% n_s

!     if (y% s(i)% pe  == dace% pe) y% s(i)% x = diag(x% b(i,i)) !+++ IBM xlf bug
      if (y% s(i)% pe  == dace% pe) then
        y% s(i)% x = diag(x% b(i,i))
      endif
    enddo
    call delete_storage (x)
    call leave_function
  end function diag_matrix
!------------------------------------------------------------------------------
  function inverse_matrix (x) result (y)
  type (t_matrix)                       :: y
  type (t_matrix) ,intent(in)           :: x
    integer :: i, j
    if (x% qual /= BDIAG) call finish ('inverse_matrix','matrix must be BDIAG')
    if (x% m_b  /= x%n_b) call finish ('inverse_matrix','x% m_b /= x% n_b')
    call enter_function
    call construct (y, ri=x% cinfo, ci=x% rinfo, name='inv(X)',qual=BDIAG)
    do   i = 1, x% m_b
      do j = 1, x% n_B
        if (y% b(i,j)% pe == dace% pe) then
          if (i==j) then
            y% b(i,j) = inverse (x% b(i,j))
          else
            call allocate_block (y% b(i,j), ZERO)
          endif
        endif
      end do
    enddo
    call delete_storage (x)
    call leave_function
  end function inverse_matrix
!------------------------------------------------------------------------------
  function inverse_rs_matrix (x, max_ev, min_ev) result (y)
  type (t_matrix)                       :: y
  type (t_matrix) ,intent(in)           :: x
  real(wp)        ,intent(in), optional :: min_ev
  real(wp)        ,intent(in), optional :: max_ev
    integer :: i
    if (x% qual /= BDIAG) call finish ('inverse_rs_matrix','matrix not BDIAG')
    if (x% m_b  /= x%n_b) call finish ('inverse_rs_matrix','x% m_b /= x% n_b')
    call enter_function
    call construct (y, ri=x% rinfo, ci=x% cinfo, name='inv(X)',qual=BDIAG)
    do i = 1, x% m_b
      if (y% b(i,i)% pe == dace% pe) then
        y% b(i,i) = inverse_rs_matrix_block (x% b(i,i), max_ev, min_ev)
      endif
    enddo
    call delete_storage (x)
    call leave_function
  end function inverse_rs_matrix
!------------------------------------------------------------------------------
  function sqrt_rs_matrix (x, min_ev) result (y)
  type (t_matrix)                       :: y
  type (t_matrix) ,intent(in)           :: x
  real(wp)        ,intent(in), optional :: min_ev
    integer :: i
    if (x% qual /= BDIAG) call finish ('sqrt_rs_matrix','matrix not BDIAG')
    if (x% m_b  /= x%n_b) call finish ('sqrt_rs_matrix','x% m_b /= x% n_b')
    call enter_function
    call construct (y, ri=x% rinfo, ci=x% cinfo, name='inv(X)',qual=BDIAG)
    do i = 1, x% m_b
      if (y% b(i,i)% pe == dace% pe) then
        y% b(i,i) = sqrt_rs_matrix_block (x% b(i,i), min_ev)
      endif
    enddo
    call delete_storage (x)
    call leave_function
  end function sqrt_rs_matrix
!------------------------------------------------------------------------------
  function pow_rs_matrix_op (x, z) result (y)
  type (t_matrix)                       :: y
  type (t_matrix) ,intent(in)           :: x
  real(wp)        ,intent(in)           :: z
    integer :: i
    if (x% qual /= BDIAG) call finish ('pow_rs_matrix','matrix not BDIAG')
    if (x% m_b  /= x%n_b) call finish ('pow_rs_matrix','x% m_b /= x% n_b')
    call enter_function
    call construct (y, ri=x% rinfo, ci=x% cinfo, name='inv(X)',qual=BDIAG)
    do i = 1, x% m_b
      if (y% b(i,i)% pe == dace% pe) then
        y% b(i,i) = pow_rs_matrix_block (x% b(i,i), z)
      endif
    enddo
    call delete_storage (x)
    call leave_function
  end function pow_rs_matrix_op
!------------------------------------------------------------------------------
  function pow_rs_matrix (x, z, min_ev) result (y)
  type (t_matrix)                       :: y
  type (t_matrix) ,intent(in)           :: x
  real(wp)        ,intent(in)           :: z
  real(wp)        ,intent(in), optional :: min_ev
    integer :: i
    if (x% qual /= BDIAG) call finish ('pow_rs_matrix','matrix not BDIAG')
    if (x% m_b  /= x%n_b) call finish ('pow_rs_matrix','x% m_b /= x% n_b')
    call enter_function
    call construct (y, ri=x% rinfo, ci=x% cinfo, name='inv(X)',qual=BDIAG)
    do i = 1, x% m_b
      if (y% b(i,i)% pe == dace% pe) then
        y% b(i,i) = pow_rs_matrix_block (x% b(i,i), z, min_ev)
      endif
    enddo
    call delete_storage (x)
    call leave_function
  end function pow_rs_matrix
!------------------------------------------------------------------------------
  function sqrt_vector (x) result (y)
  !--------------
  ! sqrt (vector)
  !--------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x
    integer :: i
#if defined (__NEC__)
    integer :: l
#endif
    call enter_function
    call construct (y, x% info, global=x% global)
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) then
#if defined (__NEC__)
        do l=1,size(y% s(i)% x)
          y% s(i)% x(l) = sqrt (x% s(i)% x(l))
        enddo
#else
        ! ****  96 Loop count is greater than that assumed by the compiler :
        ! loop-count=7371  (from bg_err_rndm)
        y% s(i)% x = sqrt (x% s(i)% x)
#endif
      endif
    end do
    call delete_storage (x)
    call leave_function
  end function sqrt_vector
!------------------------------------------------------------------------------
  function scalar_over_vector (a, x) result (y)
  !--------------
  ! real / vector
  !--------------
  type (t_vector)             :: y
  real(wp)        ,intent(in) :: a
  type (t_vector) ,intent(in) :: x
    integer :: i
    call enter_function
    call construct (y, x% info, global=x% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = a / x% s(i)% x
    end do
    call delete_storage (x)
    call leave_function
  end function scalar_over_vector
!------------------------------------------------------------------------------
  function vector_over_scalar (x, a) result (y)
  !--------------
  ! vector / real
  !--------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x
  real(wp)        ,intent(in) :: a
    real(wp) :: eda
    integer  :: i
    call enter_function
    call construct (y, x% info, global=x% global)
    eda = 1._wp / a
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x% s(i)% x * eda
    end do
    call delete_storage (x)
    call leave_function
  end function vector_over_scalar
!------------------------------------------------------------------------------
  function vector_over_int (x, a) result (y)
  !-----------------
  ! vector / integer
  !-----------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x
  integer         ,intent(in) :: a
    real(wp) :: eda
    integer  :: i
    call enter_function
    call construct (y, x% info, global=x% global)
    eda = 1._wp / a
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x% s(i)% x * eda
    end do
    call delete_storage (x)
    call leave_function
  end function vector_over_int
!------------------------------------------------------------------------------
  subroutine assign_matrix_to_matrix (y, x)
  !------------------------------------------------------------------
  ! matrix = matrix (without optional arguments, overloaded with (=))
  !------------------------------------------------------------------
  type (t_matrix) ,intent(inout)        :: y
  type (t_matrix) ,intent(in)           :: x
    call assign_matrix_optional (y, x)
  end subroutine assign_matrix_to_matrix
!------------------------------------------------------------------------------
  subroutine assign_matrix_optional (y, x, qual)
  !---------------------------------------------------
  ! matrix = matrix assignment with optional arguments
  !---------------------------------------------------
  type (t_matrix) ,intent(inout)        :: y
  type (t_matrix) ,intent(in)           :: x
  integer         ,intent(in) ,optional :: qual
    integer                        :: i, j, alloc_l
    type (t_matrix_block) ,pointer :: yb,xb
    alloc_l = y% alloc_l
    call destruct (y)
    y% rinfo   => x% rinfo
    y% cinfo   => x% cinfo
    y% m       =  x% m
    y% n       =  x% n
    y% m_b     =  x% m_b
    y% n_b     =  x% n_b
    y% qual    =  x% qual; if (present (qual)) y% qual = qual
    y% alloc_l = alloc_l
    allocate (y% b (y% m_b, y% n_b))
    do i = 1, y% m_b
      do j = 1, y% n_b
        yb => y% b (i,j)
        xb => x% b (i,j)
        if (y% qual == BDIAG) then
          yb% pe = x% b (i,i)% pe
          if (yb% pe == dace% pe) then
            if (j/=i) then
              call allocate_block (yb, ZERO, call_l=y% alloc_l)
            else
              call assign_matrix_block (yb, xb)
            endif
          else
            yb% repr = -1
          endif
        else
          yb% pe = xb% pe
          if (.not. x% rowwise) call finish ('assign_matrix_optional',&
                                             '.not. x% rowwise')
          if (dace% pe == xb% pe) then
            if ((y% qual == BLOWER .and. j>=i) .or. &
                (y% qual == BUPPER .and. j<=i)      ) then
              call allocate_block (yb, ZERO, call_l=y% alloc_l)
            else
              call assign_matrix_block (yb, xb)
            endif
          else
            yb% repr = -1
          endif
        endif
      end do
    end do
    call enter_function
    call delete_storage (x)
    call leave_function
  end subroutine assign_matrix_optional
!------------------------------------------------------------------------------
  subroutine assign_matrix_block_array (y, x)
    real(wp)              ,intent(out) :: y (:,:)
    type (t_matrix_block) ,intent(in)  :: x
      type (t_matrix_block) :: z
      if (size(y,1)/=x%m .or. size(y,2)/=x%n) then
        write(0,*) 'assign_matrix_block_array: shape(y) =',shape(y)
        write(0,*) 'assign_matrix_block_array: shape(x) =',x%m,x%n
        call finish ('assign_matrix_block_array','shape(y)/=shape(x)')
      endif
      call enter_function
      z = x
      call full_matrix (z)
      y = z% full
      call deallocate (z)
      call delete_storage (x)
      call leave_function
  end subroutine assign_matrix_block_array
!------------------------------------------------------------------------------
  subroutine assign_matrix_block (y, x)
    type (t_matrix_block) ,intent(inout) :: y
    type (t_matrix_block) ,intent(in)    :: x
      !---------------
      ! deallocate lhs
      !---------------
      call deallocate (y)
      call enter_function
      !----------------------
      ! copy meta information
      !----------------------
      y% m       = x% m
      y% n       = x% n
      y% nonzero = x% nonzero
      y% pe      = x% pe
      y% repr    = x% repr
      y% tri     = x% tri
      y% qual    = x% qual
      !------------------------------
      ! allocate arrays, copy content
      !------------------------------
      if (y% pe == dace% pe) then
        if (associated (x% full)) then
          allocate (y% full(y% m, y% n))
          y% full = x% full
        endif
        if (associated (x% packed)) then
          allocate (y% packed (size (x% packed)))
          y% packed = x% packed
          y% tri    = x% tri
        endif
        if (associated (x% ia)) then
          allocate (y% ia (size (x% ia)))
          y% ia = x% ia
        endif
        if (associated (x% ja)) then
          allocate (y% ja (size (x% ja)))
          y% ja = x% ja
        endif
        if (associated (x% nest)) then
          allocate (y% nest)
          y% nest = x% nest
        endif
      endif
      call delete_storage (x)
      call leave_function
  end subroutine assign_matrix_block
!------------------------------------------------------------------------------
  subroutine insert_full_dp (y, x, oi, oj, use_zero)
  type (t_matrix_block) ,intent(inout) :: y
  real (dp)             ,intent(in)    :: x (:,:)
  integer               ,intent(in)    :: oi   ! row    offset
  integer               ,intent(in)    :: oj   ! column offset
  logical               ,intent(in), OPTIONAL :: use_zero
    call insert_full_block (y, real(x,sp), oi, oj, use_zero)
  end subroutine insert_full_dp
!------------------------------------------------------------------------------
  subroutine insert_full_block (y, x, oi, oj, use_zero)
  type (t_matrix_block) ,intent(inout) :: y
  real (sp)             ,intent(in)    :: x (:,:)
  integer               ,intent(in)    :: oi   ! row    offset
  integer               ,intent(in)    :: oj   ! column offset
  logical               ,intent(in), OPTIONAL :: use_zero
  !---------------------------------------
  ! insert a submatrix into a matrix block
  !---------------------------------------
    integer :: ni      ! row    size
    integer :: nj      ! column size
    integer :: ns      ! number of nonzero elements
    integer :: i, j, k ! indices

    logical :: l_zero

    l_zero = .false.
    if (present(use_zero)) l_zero = use_zero

    call enter_function
      ni = size(x,1)
      nj = size(x,2)

      select case (y% repr)
      case default
        call finish ('insert_full_block','invalid repr '//srep(y% repr))

      case (CSC)
        k = y% ia (oj + 1)
        !----------------------------------
        ! increase memory if not sufficient
        !----------------------------------
        ns = count(x/=0._sp)
        if (size(y% packed)-k+1 < ns) then
          call reallocate (y, ns=k-1+ns)
        endif
        !---------------------------
        ! copy non-zero coefficients
        !---------------------------
        do j = 1, nj                     ! columns
          do i = 1, ni                   ! rows
            if (x(i,j) /= 0._sp .or. l_zero) then
              y% packed (k) = x(i,j)     ! coefficient
              y% ja     (k) = oi + i     ! row index
              k = k + 1
            endif
          end do                         ! rows
          y% ia (oj +j+1) = k            ! column index
        end do                           ! columns
        y% ia (oj +j+2:)  = k            ! zero remaining entries
        y% nonzero        = k-1          ! count non-zero entries
      end select

    call leave_function

  end subroutine insert_full_block
!------------------------------------------------------------------------------
  subroutine insert_full_block_2 (y, x, oi, oj)
  type (t_matrix_block) ,intent(inout) :: y
  real (wp)             ,intent(in)    :: x (:,:)
  integer               ,intent(in)    :: oi   ! row    offset
  integer               ,intent(in)    :: oj   ! column offset
  !----------------------------------------------------------------------------
  ! insert a submatrix into a matrix block which is already allocated
  !
  ! this routine must only be called, if insert_full_block was called before,
  ! with the same oi and oj, the same size of x and only nonzero elements in x.
  !
  ! this routine is called by the vectorised mode TSK_K for radiances to
  ! fill the H matrix into the place previously allocated.
  !----------------------------------------------------------------------------
    integer :: ni      ! row    size
    integer :: nj      ! column size
    integer :: ns      ! number of nonzero elements
    integer :: i, j, k ! indices

    call enter_function
      ni = size(x,1)
      nj = size(x,2)

      select case (y% repr)
      case default
        call finish ('insert_full_block_2','invalid repr '//srep(y% repr))

      case (CSC)
        k = y% ia (oj + 1)
        !----------------------------------
        ! programming error if memory is not sufficient
        !----------------------------------
        ns = size(x)
        if (size(y% packed)-k+1 < ns) then
          call finish ('insert_full_block_2','invalid size')
        endif
        !---------------------------
        ! copy non-zero coefficients
        !---------------------------
        do j = 1, nj                     ! columns
          do i = 1, ni                   ! rows
            y% packed (k) = x(i,j)       ! coefficient
            k = k + 1
          end do                         ! rows
        end do                           ! columns
      end select

    call leave_function

  end subroutine insert_full_block_2

!------------------------------------------------------------------------------
  subroutine insert_sparse_block (y, x, oi, oj)
  type (t_matrix_block) ,intent(inout) :: y
  type (t_matrix_block) ,intent(in)    :: x
  integer               ,intent(in)    :: oi   ! row    offset
  integer               ,intent(in)    :: oj   ! column offset
  !---------------------------------------
  ! insert a submatrix into a matrix block
  !---------------------------------------
    integer :: jx, jy, kx, ky ! indices

    call enter_function

      if (y% repr /= x% repr) &
        call finish ('insert_sparse_block','representations do not match')

      select case (y% repr)
      case default
        call finish ('insert_sparse_block','invalid repr '//srep(y% repr))

      case (CSC)
        do jx = 1, x% n                            ! columns
          jy = jx + oj
          ky = y% ia (jy)
          do kx = x% ia (jx), x% ia (jx+1)-1       ! rows
            if (x% packed(kx) /= 0._mp) then
              y% packed (ky) = x% packed(kx)       ! coefficient
              y% ja     (ky) = x% ja    (kx) + oi  ! row index
              ky = ky + 1
            endif
          end do                                   ! rows
          y% ia (jy+1 ) = ky                       ! column index
        end do                                     ! columns
        y% ia (jy+2:) = ky                         ! zero remaining entries
        y% nonzero    = ky-1                       ! count non-zero entries
      end select

    call delete_storage (x)
    call leave_function

  end subroutine insert_sparse_block
!------------------------------------------------------------------------------
  function sub_block (x, oi, ni, oj, nj, repr) result (y)
  type (t_matrix_block)                       :: y
  type (t_matrix_block) ,intent(in)           :: x
  integer               ,intent(in)           :: oi   ! row    offset
  integer               ,intent(in)           :: ni   ! row    size
  integer               ,intent(in)           :: oj   ! column offset
  integer               ,intent(in)           :: nj   ! column size
  integer               ,intent(in) ,optional :: repr ! representation
  !----------------------------------------
  ! extract a submatrix from a matrix block
  !----------------------------------------
    integer :: jy, jx, ix, ky, kx, rep
    !-----------------
    ! meta information
    !-----------------
    call enter_function
    rep = x% repr; if(present(repr)) rep = repr
    call allocate_block (y, rep, ni, nj, ns=ni*nj)
    y% pe      = x% pe
    !-------------
    ! extract body
    !-------------
    if (rep == x% repr) then
      select case (x% repr)
      case (ZERO)
        y% nonzero = 0
      case (FULL)
        y% full    = x% full(oi+1:oi+ni, oj+1:oj+nj)
        y% nonzero = count (y% full /= 0._mp)
      case (CSC)
        ky = 1
        y% ia(1) = ky
        do jy = 1, nj
          jx = jy + oj
          do kx = x% ia (jx), x% ia (jx+1)-1
            if (x% ja    (kx) >  oi+ni) exit
            if (x% ja    (kx) <= oi   ) cycle
            if (x% packed(kx) == 0._mp) cycle
            y% packed (ky) = x% packed (kx)
            y% ja     (ky) = x% ja     (kx) - oi
            ky = ky + 1
          end do
          y% ia(jy+1) = ky
        end do
        y% nonzero = ky-1
      case (PACKED)
        y% tri = x% tri
        y% nonzero = 0
        select case (x%tri)
        case ('U')
          ky = 0
          kx = 0
          do jx = 1, x%n
            if (jx >  oj + nj) exit
            if (jx <= oj) then
              kx = kx + jx
            else
              do ix = 1, jx
                kx=kx+1
                if (ix > oi .and. ix <= oi+ni) then
                  ky = ky + 1
                  y% packed(ky) = x% packed(kx)
                  if (x% packed(kx) /= 0._mp) y% nonzero = y% nonzero + 1
                endif
              end do
            end if
          end do
        case ('L')
          ky = 0
          kx = 0
          do jx = 1, x%n
            if (jx >  oj + nj) exit
            if (jx <= oj) then
              kx = kx + x%n - jx + 1
            else
              do ix = jx, x%n
                kx=kx+1
                if (ix > oi .and. ix <= oi+ni) then
                  ky = ky + 1
                  y% packed(ky) = x% packed(kx)
                  if (x% packed(kx) /= 0._mp) y% nonzero = y% nonzero + 1
                endif
              end do
            end if
          end do
        case default
          call finish ('sub_block','invalid triangle: '//x%tri)
        end select
      case default
        call finish('sub_block',                                       &
                    'not implemented: '//srep(rep)//'<-'//srep(x% repr))
      end select
    else if (rep == FULL) then
      select case (x% repr)
      case (CSC)
        y% full    = 0._mp
        y% nonzero = 0
        do jy = 1, nj
          jx = jy + oj
          do kx = x% ia (jx), x% ia (jx+1)-1
            if (x% ja    (kx) >  oi+ni) exit
            if (x% ja    (kx) <= oi   ) cycle
            if (x% packed(kx) == 0._mp) cycle
            y% full (x% ja(kx)-oi, jy) = x% packed (kx)
            y% nonzero = y% nonzero + 1
          end do
        end do
        y% nonzero = count (y% full /= 0._mp)
      case (PACKED)
        y% tri = x% tri
        y% nonzero = 0
        select case (x%tri)
        case ('U')
          kx = 0
          do jx = 1, x%n
            if (jx >  oj + nj) exit
            if (jx <= oj) then
              kx = kx + jx
            else
              do ix = 1, jx
                kx=kx+1
                if (ix > oi .and. ix <= oi+ni) then
                  y% full (ix-oi, jx-oj) = x% full (ix,jx)
                  if (x% full (ix,jx) /= 0._mp) y% nonzero = y% nonzero + 1
                endif
              end do
            end if
          end do
        case ('L')
          kx = 0
          do jx = 1, x%n
            if (jx >  oj + nj) exit
            if (jx <= oj) then
              kx = kx + x%n - jx + 1
            else
              do ix = jx, x%n
                kx=kx+1
                if (ix > oi .and. ix <= oi+ni) then
                  y% full (ix-oi, jx-oj) = x% full (ix,jx)
                  if (x% full (ix,jx) /= 0._mp) y% nonzero = y% nonzero + 1
                endif
              end do
            end if
          end do
        case default
          call finish ('sub_block','invalid triangle: '//x%tri)
        end select
      case default
        call finish('sub_block',                                       &
                    'not implemented: '//srep(rep)//'<-'//srep(x% repr))
      end select
    else
      call finish ('sub_block',                                         &
                   'invalid representation '//srep(rep)//'<-'//srep(rep))
    endif
    call release_mem    (y)
    call delete_storage (x)
    call leave_function
  end function sub_block
!------------------------------------------------------------------------------
  function inverse_matrix_block (x, max_cond_number) result (y)
    type (t_matrix_block)                       :: y
    type (t_matrix_block) ,intent(in)           :: x
    real(wp)              ,intent(in) ,optional :: max_cond_number
      integer               :: i, yrepr
      real(wp) ,allocatable :: xt (:,:)
      call enter_function
      yrepr = x% repr
      if (x% repr == CSC .or. x% repr == CSR) yrepr = FULL
      call allocate_block (y, yrepr, x%n, x%m)
      select case (x% repr)
      case (ZERO)
      case (IDENT)
      case (DIAGONAL)
        allocate (y% packed (size (x% packed)))
        do i=1,y% n
          if (x% packed(i) /= 0._mp) then
            y% packed(i) = 1._mp / x% packed(i)
          else
            y% packed(i) = 0._mp
          endif
        end do
      case (FULL)
        y% full = inverse (real(x% full,wp), max_cond_number=max_cond_number)
      case (PACKED)
        y% packed = x% packed
        y% tri    = x% tri
        call unpack_block (y)
        y% full = inverse (real(y% full,wp))
        call pack_block (y, x% tri)
      case (CSC, CSR)
        allocate (xt(x%m,x%n))
        xt = x
        y% full = inverse (xt)
        deallocate (xt)
        select case (x% repr)
        case (CSC)
          call csc_matrix (y)
        case (CSR)
          call csr_matrix (y)
        end select
      case default
        call finish ('inverse_matrix_block','not implemented: '//srep(x% repr))
      end select
      call delete_storage(x)
      call leave_function
  end function inverse_matrix_block
!------------------------------------------------------------------------------
  function inverse_rs_matrix_block (x, max_ev, min_ev) result (y)
    type (t_matrix_block)                       :: y
    type (t_matrix_block) ,intent(in)           :: x
    real(wp)              ,intent(in), optional :: min_ev
    real(wp)              ,intent(in), optional :: max_ev
      integer               :: i, yrepr
      real(wp) ,allocatable :: xt (:,:)
      call enter_function

      yrepr = x% repr
      if (x% repr == CSC .or. x% repr == CSR) yrepr = FULL
      call allocate_block (y, yrepr, x%n, x%m)

      select case (x% repr)
      case (ZERO)
      case (IDENT)
      case (DIAGONAL)
        allocate (y% packed (size (x% packed)))
        do i=1,y% n
          if (x% packed(i) /= 0._mp) then
            y% packed(i) = 1._mp / x% packed(i)
          else
            y% packed(i) = 0._mp
          endif
        end do
      case (FULL)
        y% full = inverse_rs (real(x% full,wp), max_ev, min_ev)
      case (PACKED)
        y% packed = x% packed
        y% tri    = x% tri
        call unpack_block (y)
        y% full = inverse_rs (real(y% full,wp), max_ev, min_ev)
        call pack_block (y, x% tri)
      case (CSC, CSR)
        allocate (xt(x%m,x%n))
        xt = x
        y% full = inverse_rs (xt, max_ev, min_ev)
        deallocate (xt)
        select case (x% repr)
        case (CSC)
          call csc_matrix (y)
        case (CSR)
          call csr_matrix (y)
        end select
      case default
        call finish ('inverse_rs_matrix_block',&
                     'not implemented: '//srep(x% repr))
      end select
      call delete_storage(x)
      call leave_function
  end function inverse_rs_matrix_block
!------------------------------------------------------------------------------
  function sqrt_rs_matrix_block (x, min_ev) result (y)
    type (t_matrix_block)                       :: y
    type (t_matrix_block) ,intent(in)           :: x
    real(wp)              ,intent(in), optional :: min_ev
      call enter_function
      call allocate_block (y, x% repr, x%n, x%m)
      select case (x% repr)
      case (ZERO)
      case (IDENT)
      case (DIAGONAL)
        allocate (y% packed (size (x% packed)))
        y% packed = sqrt (x% packed)
      case (FULL)
        y% full = sqrt_rs (real(x% full,wp), min_ev)
      case (PACKED)
        y% packed = x% packed
        y% tri    = x% tri
        call unpack_block (y)
        y% full = sqrt_rs (real(y% full,wp), min_ev)
        call pack_block (y, x% tri)
      case default
        call finish ('sqrt_rs_matrix_block',&
                     'not implemented: '//srep(x% repr))
      end select
      call delete_storage(x)
      call leave_function
  end function sqrt_rs_matrix_block
!------------------------------------------------------------------------------
  function pow_rs_matrix_block (x, z, min_ev) result (y)
    type (t_matrix_block)                       :: y
    type (t_matrix_block) ,intent(in)           :: x
    real(wp)              ,intent(in)           :: z
    real(wp)              ,intent(in), optional :: min_ev
      integer :: i, yrepr
      real(wp), allocatable :: xt (:,:)
      call enter_function

      yrepr = x% repr
      if (x% repr == CSC .or. x% repr == CSR) yrepr = FULL
      call allocate_block (y, yrepr, x%n, x%m)

      select case (x% repr)
      case (ZERO)
      case (IDENT)
      case (DIAGONAL)
        allocate (y% packed (size (x% packed)))
        do i=1,y% n
          if (x% packed(i) /= 0._mp) then
            y% packed(i) = x% packed(i) ** z
          else
            y% packed(i) = 0._mp
          endif
        end do
      case (FULL)
        y% full = pow_rs (real(x% full,wp), z, min_ev)
      case (PACKED)
        y% packed = x% packed
        y% tri    = x% tri
        call unpack_block (y)
        y% full = pow_rs (real(y% full,wp), z, min_ev)
        call pack_block (y, x% tri)
      case (CSC, CSR)
        allocate (xt(x%m,x%n))
        xt = x
        y% full = pow_rs (real (xt,wp), z, min_ev)
        deallocate (xt)
        select case (x% repr)
        case (CSC)
          call csc_matrix (y)
        case (CSR)
          call csr_matrix (y)
        end select
      case default
        call finish ('pow_rs_matrix_block',&
                     'not implemented: '//srep(x% repr))
      end select
      call delete_storage(x)
      call leave_function
  end function pow_rs_matrix_block
!------------------------------------------------------------------------------
  subroutine add_diagonal_to_matrix (y, x)
  type (t_matrix) ,intent(inout) :: y  ! matrix to add diagonal to
  type (t_vector) ,intent(in)    :: x  ! diagonal
  !--------------------------------
  ! add diagonal (vector) to matrix
  !--------------------------------
    integer :: l
    call enter_function
    do l = 1, y% n_b
      call add_diagonal_to_block (y% b(l,l), x% s(l)% x)
    end do
    call delete_storage (x)
    call leave_function
  end subroutine add_diagonal_to_matrix
!------------------------------------------------------------------------------
  subroutine add_diagonal_to_block (y, x)
  type (t_matrix_block) ,intent(inout) :: y     ! matrix block
  real(wp)              ,intent(in)    :: x (:) ! diagonal
  !-----------------------------------------------------------
  ! add diagonal to a matrix block (different representations)
  !-----------------------------------------------------------
    integer :: i, j, k

    if (y% pe == dace% pe) then
      select case (y% repr)
      case (ZERO)
        call allocate_block (y, DIAGONAL, call_l = y% alloc_l)
        do i = 1, y% m
          y% packed (i) = x(i)
        end do
      case (IDENT)
        call allocate_block (y, DIAGONAL, call_l = y% alloc_l)
        do i = 1, y% m
          y% packed (i) = 1._wp + x(i)
        end do
      case (FULL)
          do i = 1, y% m
            y% full(i,i) =  y% full(i,i)+x(i)
          end do
      case (PACKED)
        k=0
        select case (y% tri)
        case ('U')
          do j = 1, y% n
            do i = 1, j
              k=k+1
              if (i.eq.j) y% packed(k) =  y% packed(k)+x(i)
            end do
          end do
        case ('L')
          do j = 1, y% n
            do i = j, y% n
              k=k+1
              if (i.eq.j) y% packed(k) =  y% packed(k)+x(i)
            end do
          end do
        case default
          call finish ('add_diagonal_to_block','invalid triangle: '//y% tri)
        end select
      case default
        call finish ('add_diagonal_to_block',&
                     'not implemented for representation '//srep(y% repr))
      end select
    endif

  end subroutine add_diagonal_to_block
!------------------------------------------------------------------------------
  subroutine outer_matrix_vector (y, z)
  type (t_matrix) ,intent(inout) :: y
  type (t_vector) ,intent(in)    :: z
    type (t_matrix_block) ,pointer :: b
    type (t_vector_segm)  ,pointer :: zs
    integer                        :: i, j, k, l
    call enter_function
    do l = 1, y% n_b
      b  => y% b(l,l)
      zs => z% s(l)
      if (b% pe == dace% pe) then
        select case (b% repr)
        case (ZERO)
        case (DIAGONAL)
          do i = 1, b% m
            b% packed(i) = b% packed(i) * zs% x(i)
          end do
        case (FULL)
          do i = 1, b% m
            b% full(:,i) =  b% full(:,i) * zs% x(i)
          end do
        case (PACKED)
          select case (b%tri)
          case ('U')
            k = 0
            do j = 1, b%n
              do i = 1, j
                k=k+1
!               b% full(i,j) = b% packed(k)
                b% packed(k) = b% packed(k) * zs% x(j)
              end do
            end do
          case ('L')
            k = 0
            do j = 1, b%n
              do i = j, b%n
                k=k+1
!               b% full (i,j)= b% packed(k)
                b% packed(k) = b% packed(k) * zs% x(j)
              end do
            end do
          case default
            call finish ('unpack_block','invalid triangle: '//b%tri)
          end select
        case default
          call finish('outer_matrix_vector',                            &
                      'representation not implemented: '//srep(b% repr))
        end select
      endif
    end do
    call delete_storage (z)
    call leave_function
  end subroutine outer_matrix_vector
!------------------------------------------------------------------------------
  subroutine outer_vector_matrix (x, y)
  type (t_vector) ,intent(in)    :: x
  type (t_matrix) ,intent(inout) :: y
    type (t_matrix_block) ,pointer :: b
    type (t_vector_segm)  ,pointer :: xs
    integer                        :: i, j, k, l
    call enter_function
    do l = 1, y% n_b
      b  => y% b(l,l)
      xs => x% s(l)
      if (b% pe == dace% pe) then
        select case (b% repr)
        case (ZERO)
        case (DIAGONAL)
          do i = 1, b% m
            b% packed(i) = b% packed(i) * xs% x(i)
          end do
        case (FULL)
          do i = 1, b% m
            b% full(i,:) =  b% full(i,:) * xs% x(i)
          end do
        case (PACKED)
          select case (b%tri)
          case ('U')
            k = 0
            do j = 1, b%n
              do i = 1, j
                k=k+1
!               b% full (i,j)= b% packed(k)
                b% packed(k) = b% packed(k) * xs% x(i)
              end do
            end do
          case ('L')
            k = 0
            do j = 1, b%n
              do i = j, b%n
                k=k+1
!               b% full (i,j)= b% packed(k)
                b% packed(k) = b% packed(k) * xs% x(i)
              end do
            end do
          case default
            call finish ('unpack_block','invalid triangle: '//b%tri)
          end select
        case default
          call finish('outer_vector_matrix',                            &
                      'representation not implemented: '//srep(b% repr))
        end select
      endif
    end do
    call delete_storage (x)
    call leave_function
  end subroutine outer_vector_matrix
!------------------------------------------------------------------------------
  subroutine add_matrix_to_matrix (y, x)
  !--------------------------------------------
  ! add matrix X to matrix Y element by element
  !--------------------------------------------
  type (t_matrix) ,intent(inout) :: y
  type (t_matrix) ,intent(in)    :: x
    integer                        :: i, j, l, i1, j1
    type (t_matrix_block) ,pointer :: yb, xb
    real(wp) ,allocatable :: diag(:)
    !--------------------------------------
    ! check for consistent matrix qualities
    !--------------------------------------
    if ((y% qual == BDIAG  .and. x% qual /= BDIAG)   .or. &
        (y% qual == BUPPER .and. x% qual == GENERAL) .or. &
        (y% qual == BUPPER .and. x% qual == BLOWER)  .or. &
        (y% qual == BLOWER .and. x% qual == GENERAL) .or. &
        (y% qual == BLOWER .and. x% qual == BUPPER)) then
      call finish ('add_matrix_to_matrix','matrix quality mismatch')
    endif
    call enter_function
    do i = 1, y% m_b
      do j = 1, y% n_b
        yb => y% b (i,j)
        xb => x% b (i,j)
        !---------------------------------------------------
        ! if right hand side == ZERO or undefined do nothing
        !---------------------------------------------------
        if ((xb % repr == ZERO) .or. &
            (xb % repr == NOT_INIT)) then
        else if (xb% pe == dace% pe .and. yb% pe == dace% pe) then
          !------------------------------------
          ! if left hand side == ZERO just copy
          !------------------------------------
          if (yb% repr == ZERO) then
            call assign_matrix_block (yb, x% b(i,j))
          !---------------------------------------------
          ! else add (y,x must have same representation)
          !---------------------------------------------
          else if (yb% repr == xb% repr) then
            select case (yb% repr)
            case (ZERO)
            case (FULL)
              yb% full = yb% full + xb% full
            case (PACKED, DIAGONAL)
              if (yb% tri /= xb% tri) &
                call finish ('add_matrix_to_matrix',&
                             'triangles do not match')
#if defined (__NEC__)
              do l=1,size(xb% packed)
                yb% packed(l) = yb% packed(l) + xb% packed(l)
              end do
#else
              ! ****  96 Loop count is greater than that assumed by the compiler :
              !       loop-count=13695
              yb% packed = yb% packed + xb% packed
#endif
            case default
              call finish ('add_matrix_to_matrix',&
                           'not implemented for representation')
            end select
          !--------------------------
          ! else if rhs is IDENT: add
          !--------------------------
          else if (xb% repr == IDENT) then
            select case (yb% repr)
!           case (-1)
!             yb% repr = IDENT
            case (DIAGONAL)
              do l = 1, yb%n
                yb% packed(l) = yb% packed(l) + 1._mp
              end do
            case (FULL)
              do l = 1, yb%n
                yb% full(l,l) = yb% full(l,l) + 1._mp
              end do
            case (PACKED)
              l=0
              select case (yb% tri)
              case ('U')
                do j1 = 1, yb% n
                  do i1 = 1, j1
                  l=l+1
                  if (i1.eq.j1) yb% packed(l) =  yb% packed(l) + 1._mp
                  end do
                end do
              case ('L')
                do j1 = 1, yb% n
                  do i1 = j1, yb% n
                   l=l+1
                   if (i1.eq.j1) yb% packed(l) =  yb% packed(l) + 1._mp
                   end do
                end do
              case default
                call finish ('scale','invalid triangle: '//yb% tri)
              end select
            case default
              call finish ('add_unity_to_matrix',&
                           'not implemented for representation')
            end select
          elseif (yb% repr == diagonal) then
            !--------------------------------------------------------------
            ! lhs and rhs don't have same representation: some simple cases
            !--------------------------------------------------------------
            allocate (diag (yb% m))
            diag = yb% packed
            call assign_matrix_block   (yb, x% b(i,j))
            call add_diagonal_to_block (yb, diag)
            deallocate (diag)
          else
            !------------------------------------------------
            ! else lhs and rhs don't have same representation
            !------------------------------------------------
            write (0,*) 'add_matrix_to_matrix, i,j=', i, j,    &
              '; lhs: pe,repr=', y%name, yb%pe, srep(yb%repr), &
              '; rhs: pe,repr=', x%name, xb%pe, srep(xb%repr)
            call finish ("add_matrix_to_matrix",&
                         "lhs and rhs don't have same representation")
          endif
        elseif (xb% pe /= dace% pe .and. yb% pe /= dace% pe) then
        else
          write(0,*)  'add_matrix_to_matrix: i,j,y%pe,x%pe=',i,j,yb%pe,xb%pe
          call finish('add_matrix_to_matrix','processor mismatch: '&
                      //y%name//x%name)
        endif
      end do
    end do
    call delete_storage (x)
    call leave_function
  end subroutine add_matrix_to_matrix
!------------------------------------------------------------------------------
  subroutine pack_dec_matrix (x, tri)
  !-----------------------------------------
  ! transform: FULL to PACKED representation
  !-----------------------------------------
  type (t_matrix)     ,intent(inout) :: x
  character ,optional ,intent(in)    :: tri
    integer                          :: i
    if (x% m   /= x% n  ) call finish('pack_dec_matrix','matrix is not square')
    if (x% m_b /= x% n_b) call finish('pack_dec_matrix','matrix is not square')
    do i = 1, x% m_b
      call pack_block (x% b(i,i), tri)
    end do
  end subroutine pack_dec_matrix
!------------------------------------------------------------------------------
  subroutine cholesky (x, tri, info)
  !--------------------------------------------------
  ! perform Cholesky factorization of diagonal blocks
  !--------------------------------------------------
  type (t_matrix)     ,intent(inout) :: x
  character ,optional ,intent(in)    :: tri
  integer   ,optional ,intent(out)   :: info
    integer                          :: i
    integer                          :: linfo(x% m_b) ! 'info' for blocks

    if (x% m   /= x% n  ) call finish('cholesky','matrix is not square')
    if (x% m_b /= x% n_b) call finish('cholesky','matrix is not square')
    if (present(info)) info = 0
    linfo = 0
! OpenMP disabled here since some LAPACK/BLAS libraries may not be thread-safe
!!$omp parallel do private(i) shared(linfo) schedule(dynamic)
    do i = 1, x% m_b
      if (x% b(i,i)% pe == dace% pe) call cholesky_block (x% b(i,i), tri, &
                                                      info=linfo(i))
#ifndef _OPENMP
      if(present(info)) then
        if (linfo(i) /= 0) exit
      endif
#endif
    end do
!!$omp end parallel do
    if (present(info)) info = p_max (maxval (abs (linfo)))
  end subroutine cholesky
!------------------------------------------------------------------------------
  subroutine LU_solv (x)
  !--------------------------------------------------
  ! perform LU factorization of diagonal blocks
  !--------------------------------------------------
  type (t_matrix)     ,intent(inout) :: x
    integer                          :: i
    if (x% m   /= x% n  ) call finish('LU','matrix is not square')
    if (x% m_b /= x% n_b) call finish('LU','matrix is not square')
    do i = 1, x% m_b
      if (x% b(i,i)% pe == dace% pe) call LU_block (x% b(i,i))
    end do
  end subroutine LU_solv
!------------------------------------------------------------------------------
  subroutine scale_matrix (x, vl, vr)
  type (t_matrix) ,intent(inout) :: x  ! matrix
  type (t_vector) ,intent(in)    :: vl ! vector to multiply from the left
  type (t_vector) ,intent(in)    :: vr ! vector to multiply from the right
    type (t_matrix_block) ,pointer :: xb
    type (t_vector_segm)  ,pointer :: sl, sr
    integer                        :: l, m
    if (x% m   /= vl% n  ) call finish('scale','x% m /= vl% n')
    if (x% n   /= vr% n  ) call finish('scale','x% n /= vl% n')
    if (x% m_b /= vl% n_s) call finish('scale','x% m_b /= vl% n_s')
    if (x% n_b /= vr% n_s) call finish('scale','x% n_b /= vr% n_s')
    call enter_function
    do l = 1, x% m_b
      do m = 1, x% n_b
        xb =>  x% b(l,m)
        if (xb% pe   /= dace% pe)     cycle
        if (xb% repr == ZERO)     cycle
        if (xb% repr == MIRROR)   cycle
        if (xb% repr == NOT_INIT) cycle
        sl => vl% s(l)
        sr => vr% s(m)
        select case (xb% repr)
        case (IDENT)
          call allocate_block (xb, DIAGONAL, call_l = x% alloc_l)
          xb% packed (:) = sl% x(:)*sr% x(:)
        case default
          call scale_block (xb, sl% x, sr% x)
        end select
      end do
    end do
    call delete_storage (vl)
    call delete_storage (vr)
    call leave_function
  end subroutine scale_matrix
!------------------------------------------------------------------------------
  subroutine scale_block (xb, l, r)
  type(t_matrix_block) ,intent(inout) :: xb    ! matrix
  real(wp)   ,optional ,intent(in)    :: l(:)  ! vector to multiply (left)
  real(wp)   ,optional ,intent(in)    :: r(:)  ! vector to multiply (right)

    integer :: i, j, k
    integer :: lr

    lr = 0
    if (present(l)) lr=lr+1
    if (present(r)) lr=lr+2
    if (lr==0) return

    if (xb% pe == dace% pe) then
      select case (lr)
      case (3)
        select case (xb% repr)
        case (IDENT)
          call allocate_block (xb, DIAGONAL)
          xb% packed (:) = l(:)*r(:)
        case (FULL)
          do j = 1, xb% n
            do i = 1, xb% m
              xb% full(i,j) =  xb% full(i,j)*l(i)*r(j)
            end do
          end do
        case (PACKED)
          k=0
          select case (xb% tri)
          case ('U')
            do j = 1, xb% n
              do i = 1, j
                k=k+1
                xb% packed(k) =  xb% packed(k)*l(i)*r(j)
              end do
            end do
          case ('L')
            do j = 1, xb% n
              do i = j, xb% n
                k=k+1
                xb% packed(k) =  xb% packed(k)*l(i)*r(j)
              end do
            end do
          case default
            call finish ('scale','invalid triangle: '//xb% tri)
          end select
        case (CSR)
          do i = 1, xb% m
            do k = xb% ia (i), xb% ia (i+1) - 1
              j = xb% ja (k)
              xb% packed (k) = xb% packed (k) * l(i) * r(j)
            end do
          end do
        case (CSC)
          do j = 1, xb% m
            do k = xb% ia (j), xb% ia (j+1) - 1
              i = xb% ja (k)
              xb% packed (k) = xb% packed (k) * l(i) * r(j)
            end do
          end do
        case default
          call finish ('scale',&
                       'not implemented for representation '//srep(xb% repr))
        end select

      case (1)
        select case (xb% repr)
        case (IDENT)
          call allocate_block (xb, DIAGONAL)
          xb% packed (:) = l(:)
        case (FULL)
          do j = 1, xb% n
            do i = 1, xb% m
              xb% full(i,j) =  xb% full(i,j)*l(i)
            end do
          end do
        case (PACKED)
          k=0
          select case (xb% tri)
          case ('U')
            do j = 1, xb% n
              do i = 1, j
                k=k+1
                xb% packed(k) =  xb% packed(k)*l(i)
              end do
            end do
          case ('L')
            do j = 1, xb% n
              do i = j, xb% n
                k=k+1
                xb% packed(k) =  xb% packed(k)*l(i)
              end do
            end do
          case default
            call finish ('scale','invalid triangle: '//xb% tri)
          end select
        case (CSR)
          do i = 1, xb% m
            do k = xb% ia (i), xb% ia (i+1) - 1
              j = xb% ja (k)
              xb% packed (k) = xb% packed (k) * l(i)
            end do
          end do
        case (CSC)
          do j = 1, xb% m
            do k = xb% ia (j), xb% ia (j+1) - 1
              i = xb% ja (k)
              xb% packed (k) = xb% packed (k) * l(i)
            end do
          end do
        case default
          call finish ('scale',&
                       'not implemented for representation '//srep(xb% repr))
        end select

      case (2)
        select case (xb% repr)
        case (IDENT)
          call allocate_block (xb, DIAGONAL)
          xb% packed (:) = r(:)
        case (FULL)
          do j = 1, xb% n
            do i = 1, xb% m
              xb% full(i,j) =  xb% full(i,j)*r(j)
            end do
          end do
        case (PACKED)
          k=0
          select case (xb% tri)
          case ('U')
            do j = 1, xb% n
              do i = 1, j
                k=k+1
                xb% packed(k) =  xb% packed(k)*r(j)
              end do
            end do
          case ('L')
            do j = 1, xb% n
              do i = j, xb% n
                k=k+1
                xb% packed(k) =  xb% packed(k)*r(j)
              end do
            end do
          case default
            call finish ('scale','invalid triangle: '//xb% tri)
          end select
        case (CSR)
          do i = 1, xb% m
            do k = xb% ia (i), xb% ia (i+1) - 1
              j = xb% ja (k)
              xb% packed (k) = xb% packed (k) * r(j)
            end do
          end do
        case (CSC)
          do j = 1, xb% m
            do k = xb% ia (j), xb% ia (j+1) - 1
              i = xb% ja (k)
              xb% packed (k) = xb% packed (k) * r(j)
            end do
          end do
        case default
          call finish ('scale',&
                       'not implemented for representation '//srep(xb% repr))
        end select
      end select
    endif
  end subroutine scale_block
!------------------------------------------------------------------------------
  function vector_times_vector (x1, x2) result (y)
  !----------------
  ! vector * vector
  !----------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'v*v', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x * x2% s(i)% x
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
#if defined (__GFORTRAN__) && (__GNUC__ * 100 + __GNUC_MINOR__) >= 407
    ! Workaround for GNU Fortran 4.7 (prerelease)
    if (call_level < -1) print *, "vector_times_vector: barrier"
#endif
  end function vector_times_vector
!------------------------------------------------------------------------------
  function vector_over_vector (x1, x2) result (y)
  !----------------
  ! vector / vector
  !----------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'v/v', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x / x2% s(i)% x
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
#if defined (__GFORTRAN__) && (__GNUC__ * 100 + __GNUC_MINOR__) >= 407
    ! Workaround for GNU Fortran 4.7 (prerelease)
    if (call_level < -1) print *, "vector_over_vector: barrier"
#endif
  end function vector_over_vector
!------------------------------------------------------------------------------
  function add_vectors (x1, x2) result (y)
  !----------------
  ! vector + vector
  !----------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
#if defined (__NEC__)
    integer :: l
#endif
    call enter_function
    call construct (y, x1% info, 'v+v', x1% global .and. x2% global)
    do i = 1, y% n_s
#if defined (__NEC__)
      if (associated(y% s(i)% x)) then
        do l=1,size(y% s(i)% x)
          y% s(i)% x(l) = x1% s(i)% x(l) + x2% s(i)% x(l)
        enddo
      endif
#else
      ! ****  96 Loop count is greater than that assumed by the compiler
      if (associated (y% s(i)% x)) y% s(i)% x = x1% s(i)% x + x2% s(i)% x
#endif
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
#if defined (__GFORTRAN__) && (__GNUC__ * 100 + __GNUC_MINOR__) >= 407
    ! Workaround for GNU Fortran 4.7 (prerelease)
    if (call_level < -1) print *, "add_vectors: barrier"
#endif
  end function add_vectors
!------------------------------------------------------------------------------
  function subtract_vectors (x1, x2) result (y)
  !----------------
  ! vector - vector
  !----------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
#if defined (__NEC__)
    integer :: l
#endif
    call enter_function
    call construct (y, x1% info, 'v-v', x1% global .and. x2% global)
    do i = 1, y% n_s
#if defined (__NEC__)
      if (associated(y% s(i)% x)) then
        do l=1,size(y% s(i)% x)
          y% s(i)% x(l) = x1% s(i)% x(l) - x2% s(i)% x(l)
        enddo
      endif
#else
      ! ****  96 Loop count is greater than that assumed by the compiler
      if (associated (y% s(i)% x)) y% s(i)% x = x1% s(i)% x - x2% s(i)% x
#endif
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
#if defined (__GFORTRAN__) && (__GNUC__ * 100 + __GNUC_MINOR__) >= 407
    ! Workaround for GNU Fortran 4.7 (prerelease)
    if (call_level < -1) print *, "subtract_vectors: barrier"
#endif
  end function subtract_vectors
!------------------------------------------------------------------------------
  function minus_vector (x) result (y)
  !---------
  ! - vector
  !---------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x
    integer :: i
    call enter_function
    call construct (y, x% info, '-v', x% global)
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = - x% s(i)% x
    end do
    call delete_storage (x)
    call leave_function
  end function minus_vector
!------------------------------------------------------------------------------
  function abs_vector (x) result (y)
  !-------------
  ! abs (vector)
  !-------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x
    integer :: i
    call enter_function
    call construct (y, x% info, '|v|', x% global)
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = abs (x% s(i)% x)
    end do
    call delete_storage (x)
    call leave_function
  end function abs_vector

!------------------------------------------------------------------------------
  function log_vector (x) result (y)
  !-------------
  ! log (vector)
  !-------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x
    integer :: i
#if defined (__NEC__)
    integer :: l
#endif
    call enter_function
    call construct (y, x% info, '|v|', x% global)
    do i = 1, y% n_s
#if defined (__NEC__)
      if (associated (y% s(i)% x)) then
        do l=1,size(y% s(i)% x)
          y% s(i)% x(l) = log (x% s(i)% x(l))
        enddo
      endif
#else
      ! ****  96 Loop count is greater than that assumed by the compiler :
      if (associated (y% s(i)% x)) y% s(i)% x = log (x% s(i)% x)
#endif
    end do
    call delete_storage (x)
    call leave_function
  end function log_vector
!------------------------------------------------------------------------------
  function max_vector_real (x1, x2) result (y)
  !-------------------
  ! max (vector, real)
  !-------------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x1
  real(wp)        ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'max(v,s)', x1% global)
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = max (x1% s(i)% x, x2)
    end do
    call delete_storage (x1)
    call leave_function
  end function max_vector_real
!------------------------------------------------------------------------------
  function max_vector_vector (x1, x2) result (y)
  !---------------------
  ! max (vector, vector)
  !---------------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'max(v,v)', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = max(x1% s(i)% x, x2% s(i)% x)
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
  end function max_vector_vector
!------------------------------------------------------------------------------
  function min_vector_real (x1, x2) result (y)
  !-------------------
  ! min (vector, real)
  !-------------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x1
  real(wp)        ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'min(v,r)', x1% global)
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = min (x1% s(i)% x, x2)
    end do
    call delete_storage (x1)
    call leave_function
  end function min_vector_real
!------------------------------------------------------------------------------
  function min_vector_vector (x1, x2) result (y)
  !---------------------
  ! min (vector, vector)
  !---------------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'min(v,v)', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = min(x1% s(i)% x, x2% s(i)% x)
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
  end function min_vector_vector
!------------------------------------------------------------------------------
  function vector_minus_scalar (x1, x2) result (y)
  !--------------
  ! vector - real
  !--------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x1
  real(wp)        ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'v-s', x1% global)
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = x1% s(i)% x - x2
    end do
    call delete_storage (x1)
    call leave_function
  end function vector_minus_scalar
!------------------------------------------------------------------------------
  function scalar_minus_vector (x1, x2) result (y)
  !--------------
  ! real - vector
  !--------------
  type (t_vector)             :: y
  real(wp)        ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x2% info, 's-v', x2% global)
    do i = 1, y% n_s
      if (associated (y% s(i)% x)) y% s(i)% x = x1 - x2% s(i)% x
    end do
    call delete_storage (x2)
    call leave_function
  end function scalar_minus_vector
!------------------------------------------------------------------------------
  function vector_exp_integer (x, i) result (y)
  !--------------
  ! vector - real
  !--------------
  type (t_vector)             :: y
  type (t_vector) ,intent(in) :: x
  integer         ,intent(in) :: i
    integer :: j
#if defined (__NEC__)
    integer :: l
#endif
    call enter_function
    call construct (y, x% info, 'v**i', x% global)
    do j = 1, y% n_s
      if (associated (y% s(j)% x)) then
#if defined (__NEC__)
        do l=1,size(y% s(j)% x)
          y% s(j)% x(l) = x% s(j)% x(l) ** i
        enddo
#else
        ! ****  96 Loop count is greater than that assumed by the compiler :
        !          loop-count=7371   (from bg_err_rndm)
        y% s(j)% x = x% s(j)% x ** i
#endif
      endif
    end do
    call delete_storage (x)
    call leave_function
  end function vector_exp_integer
!-----------------------------------------------------------------------------
  function vector_equal_vector (x1, x2) result (y)
  !-------------------
  ! (vector == vector)
  !-------------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1==x2', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x == x2% s(i)% x
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
  end function vector_equal_vector
!-----------------------------------------------------------------------------
  function vector_equal_real (x1, x2) result (y)
  !-----------------
  ! (vector /= real)
  !-----------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  real(wp)        ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1==x2', x1% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x == x2
    end do
    call delete_storage (x1)
    call leave_function
  end function vector_equal_real
!-----------------------------------------------------------------------------
  function vector_ne_vector (x1, x2) result (y)
  !-------------------
  ! (vector /= vector)
  !-------------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1/=x2', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x /= x2% s(i)% x
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
  end function vector_ne_vector
!-----------------------------------------------------------------------------
  function vector_ne_real (x1, x2) result (y)
  !-----------------
  ! (vector /= real)
  !-----------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  real(wp)        ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1/=x2', x1% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x /= x2
    end do
    call delete_storage (x1)
    call leave_function
  end function vector_ne_real
!-----------------------------------------------------------------------------
  function vector_lt_vector (x1, x2) result (y)
  !-------------------
  ! (vector < vector)
  !-------------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1<x2', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x < x2% s(i)% x
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
  end function vector_lt_vector
!-----------------------------------------------------------------------------
  function vector_lt_real (x1, x2) result (y)
  !-----------------
  ! (vector < real)
  !-----------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  real(wp)        ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1<x2', x1% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x < x2
    end do
    call delete_storage (x1)
    call leave_function
  end function vector_lt_real
!-----------------------------------------------------------------------------
  function vector_le_vector (x1, x2) result (y)
  !-------------------
  ! (vector <= vector)
  !-------------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1<=x2', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x <= x2% s(i)% x
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
  end function vector_le_vector
!-----------------------------------------------------------------------------
  function vector_le_real (x1, x2) result (y)
  !-----------------
  ! (vector <= real)
  !-----------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  real(wp)        ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1<=x2', x1% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x <= x2
    end do
    call delete_storage (x1)
    call leave_function
  end function vector_le_real
!-----------------------------------------------------------------------------
  function vector_gt_vector (x1, x2) result (y)
  !-------------------
  ! (vector > vector)
  !-------------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1>x2', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x > x2% s(i)% x
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
  end function vector_gt_vector
!-----------------------------------------------------------------------------
  function vector_gt_real (x1, x2) result (y)
  !-----------------
  ! (vector > real)
  !-----------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  real(wp)        ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1>x2', x1% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x > x2
    end do
    call delete_storage (x1)
    call leave_function
  end function vector_gt_real
!-----------------------------------------------------------------------------
  function vector_ge_vector (x1, x2) result (y)
  !-------------------
  ! (vector >= vector)
  !-------------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1>=x2', x1% global .and. x2% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x >= x2% s(i)% x
    end do
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
  end function vector_ge_vector
!-----------------------------------------------------------------------------
  function vector_ge_real (x1, x2) result (y)
  !-----------------
  ! (vector >= real)
  !-----------------
  type (t_bvector)            :: y
  type (t_vector) ,intent(in) :: x1
  real(wp)        ,intent(in) :: x2
    integer :: i
    call enter_function
    call construct (y, x1% info, 'x1>=x2', x1% global)
    do i = 1, y% n_s
      if (associated(y% s(i)% x)) y% s(i)% x = x1% s(i)% x >= x2
    end do
    call delete_storage (x1)
    call leave_function
  end function vector_ge_real
!------------------------------------------------------------------------------
  function dot_vector (x1, x2) result (y)
  !-----------------------------
  ! dot_product (vector, vector)
  !-----------------------------
  real(wp)                    :: y
  type (t_vector) ,intent(in) :: x1
  type (t_vector) ,intent(in) :: x2
    integer  :: i
    real(wp) :: z (x1% n_s)
    call enter_function
    z = 0._wp
    if (x1% global .and. x2% global) then
      !---------------------------
      ! globally allocated fields:
      !---------------------------
      do i = 1, x1% n_s
!        z(i) = ddot (x1% s(i)% n ,x1% s(i)% x ,1 ,x2% s(i)% x, 1)
         z(i) = dot_product       (x1% s(i)% x,    x2% s(i)% x)
      end do
    else
      !--------------------------
      ! locally allocated fields:
      !--------------------------
      do i = 1, x1% n_s
        if (dace% pe == x1% s(i)% pe) &
!         z(i) = ddot (x1% s(i)% n ,x1% s(i)% x ,1 ,x2% s(i)% x, 1)
          z(i) = dot_product       (x1% s(i)% x,    x2% s(i)% x)
      end do
      z = p_sum (z)
    endif
    y = sum (z)
    call delete_storage (x1)
    call delete_storage (x2)
    call leave_function
  end function dot_vector
!------------------------------------------------------------------------------
  function sum_vector (x, m) result (y)
  !---------------------------
  ! sum elements of the vector
  !---------------------------
  real(wp)                               :: y
  type (t_vector)  ,intent(in)           :: x
  type (t_bvector) ,intent(in) ,optional :: m ! not called mask due to
    integer  :: i                             ! ifort compiler bug
    real(wp) :: z (x% n_s)
    call enter_function
    y = 0._wp
    if (x% global) then
      !---------------------------
      ! globally allocated fields:
      !---------------------------
      if (present (m)) then
        call gather (m)
        do i = 1, x% n_s
          y = y + sum (x% s(i)% x, mask=m% s(i)% x)
        end do
        call release_mem (m)
      else
        do i = 1, x% n_s
          y = y + sum (x% s(i)% x)
        end do
      endif
    else
      !--------------------------
      ! locally allocated fields:
      !--------------------------
      z = 0._wp
      do i = 1, x% n_s
        if (dace% pe == x% s(i)% pe) then
          if (present (m)) then
            z (i) = sum (x% s(i)% x, mask=m% s(i)% x)
          else
            z (i) = sum (x% s(i)% x)
          endif
        endif
      end do
      z = p_sum (z)
      do i = 1, x% n_s
        y = y + z (i)
      end do
    endif
    call delete_storage (x)
    call leave_function
  end function sum_vector
!------------------------------------------------------------------------------
  function count_vector (x) result (y)
  !----------------------------------
  ! count true elements of the vector
  !----------------------------------
  integer                      :: y
  type (t_bvector) ,intent(in) :: x
    integer  :: i
    integer  :: z
    call enter_function
    if (x% global) then
      !---------------------------
      ! globally allocated fields:
      !---------------------------
      y = 0
      do i = 1, x% n_s
        y = y + count (x% s(i)% x)
      end do
    else
      !--------------------------
      ! locally allocated fields:
      !--------------------------
      z = 0
      do i = 1, x% n_s
        if (dace% pe == x% s(i)% pe) z = z + count (x% s(i)% x)
      end do
      y = p_sum (z)
    endif
    call delete_storage (x)
    call leave_function
  end function count_vector
!------------------------------------------------------------------------------
  function any_vector (x) result (y)
  !----------------------------------
  ! count true elements of the vector
  !----------------------------------
  logical                      :: y
  type (t_bvector) ,intent(in) :: x
    integer  :: i
    logical  :: z
    call enter_function
    if (x% global) then
      !---------------------------
      ! globally allocated fields:
      !---------------------------
      y = .false.
      do i = 1, x% n_s
        y = y .or. any (x% s(i)% x)
        if (y) exit
      end do
    else
      !--------------------------
      ! locally allocated fields:
      !--------------------------
      z = .false.
      do i = 1, x% n_s
        if (dace% pe == x% s(i)% pe) z = z .or. any (x% s(i)% x)
        if (z) exit
      end do
      y = p_or (z)
    endif
    call delete_storage (x)
    call leave_function
  end function any_vector
!------------------------------------------------------------------------------
  function all_vector (x) result (y)
  !----------------------------------
  ! count true elements of the vector
  !----------------------------------
  logical                      :: y
  type (t_bvector) ,intent(in) :: x
    integer  :: i
    logical  :: z
    call enter_function
    if (x% global) then
      !---------------------------
      ! globally allocated fields:
      !---------------------------
      y = .true.
      do i = 1, x% n_s
        y = y .and. all (x% s(i)% x)
      end do
    else
      !--------------------------
      ! locally allocated fields:
      !--------------------------
      z = .true.
      do i = 1, x% n_s
        if (dace% pe == x% s(i)% pe) z = z .and. all (x% s(i)% x)
      end do
      y = p_and (z)
    endif
    call delete_storage (x)
    call leave_function
  end function all_vector
!------------------------------------------------------------------------------
  function maxval_vector (x) result (y)
  !---------------------------
  ! maxval elements of the vector
  !---------------------------
  real(wp)                    :: y
  type (t_vector) ,intent(in) :: x
    integer  :: i
    call enter_function
    y = - huge (y)
    if (x% global) then
      !---------------------------
      ! globally allocated fields:
      !---------------------------
      do i = 1, x% n_s
        y = max (y, maxval (x% s(i)% x))
      end do
    else
      !--------------------------
      ! locally allocated fields:
      !--------------------------
      do i = 1, x% n_s
        if (dace% pe == x% s(i)% pe) y = max (y, maxval (x% s(i)% x))
      end do
      y = p_max (y)
    endif
    call delete_storage (x)
    call leave_function
  end function maxval_vector
!------------------------------------------------------------------------------
  function minval_vector (x) result (y)
  !---------------------------
  ! minval elements of the vector
  !---------------------------
  real(wp)                    :: y
  type (t_vector) ,intent(in) :: x
    integer  :: i
    call enter_function
    y = huge (y)
    if (x% global) then
      !---------------------------
      ! globally allocated fields:
      !---------------------------
      do i = 1, x% n_s
        y = min (y, minval (x% s(i)% x))
      end do
    else
      !--------------------------
      ! locally allocated fields:
      !--------------------------
      do i = 1, x% n_s
        if (dace% pe == x% s(i)% pe) y = min (y, minval (x% s(i)% x))
      end do
      y = p_min (y)
    endif
    call delete_storage (x)
    call leave_function
  end function minval_vector
!------------------------------------------------------------------------------
  function size_of_vector (x)
  !---------------------------------
  ! number of elements of the vector
  !---------------------------------
  integer                     :: size_of_vector
  type (t_vector) ,intent(in) :: x
    call enter_function
    size_of_vector = x% n
    call delete_storage (x)
    call leave_function
  end function size_of_vector
!------------------------------------------------------------------------------
  function norm2 (x) result (y)
  !-------------------
  ! 2-norm of a vector
  !-------------------
  real(wp)                    :: y
  type (t_vector) ,intent(in) :: x
    integer :: i
    real(wp) :: z (x% n_s)
    call enter_function
    y = 0._wp
    if (x% global) then
      !---------------------------
      ! globally allocated fields:
      !---------------------------
      do i = 1, x% n_s
        y = y + dnrm2 (x% s(i)% n ,x% s(i)% x ,1) ** 2
      end do
    else
      !--------------------------
      ! locally allocated fields:
      !--------------------------
      z = 0._wp
      do i = 1, x% n_s
        if (dace% pe == x% s(i)% pe) z (i) = dnrm2 (x% s(i)% n ,x% s(i)% x ,1) ** 2
      end do
      z = p_sum (z)
      do i = 1, x% n_s
        y = y + z (i)
      end do
    endif
    y = sqrt (y)
    call delete_storage (x)
    call leave_function
  end function norm2
!------------------------------------------------------------------------------
  function amax (x) result (y)
  !--------------------------
  ! maximum of absolute value
  !--------------------------
  real(wp)                    :: y
  type (t_vector) ,intent(in) :: x
    integer :: i
    call enter_function
    y = 0._wp
    do i = 1, x% n_s
      if(associated(x% s(i)% x)) y = max (y, maxval(abs(x% s(i)% x)))
    end do
    if (.not.x% global) y = p_max (y)
    call delete_storage (x)
    call leave_function
  end function amax
!------------------------------------------------------------------------------
  subroutine random_gauss_vector (x)
  !--------------------------------------------------
  ! set vector to random field (for testing purposes)
  !--------------------------------------------------
  type (t_vector) ,intent(inout) :: x
    integer :: i
    type (t_vector_segm) ,pointer :: s
    do i = 1, x% n_s
      s => x% s(i)
      if (s% pe == dace% pe) call random_gauss (s% x, x% info% b(i)% seed)
      if (x% global) call p_bcast (s% x, s% pe)
    end do
  end subroutine random_gauss_vector
!------------------------------------------------------------------------------
  subroutine random_vector (x)
  !--------------------------------------------------
  ! set vector to random field (for testing purposes)
  !--------------------------------------------------
  type (t_vector) ,intent(inout) :: x
    integer :: i
    type (t_vector_segm) ,pointer :: s
    do i = 1, x% n_s
      s => x% s(i)
      if (dace% lpio) then
        !----------------------------------------
        ! on dace% pe==p_io: call random_number, send
        !----------------------------------------
        if (.not. associated (s% x)) allocate (s% x(s% n))
        call random_number (s% x)
        if (x% global) then
          call p_bcast (s% x, dace% pio)
        else if (dace% pe /= s% pe) then
          call p_send  (s% x, s% pe, p_tag=i)
        endif
        if (.not.(x% global.or.dace% pe==s% pe)) deallocate (s% x)
      else
        !-----------------------
        ! on dace% pe/=p_io: receive
        !-----------------------
        if (x% global) then
          call p_bcast (s% x, dace% pio)
        else if (.not.dace% lpio) then
          if (s% pe == dace% pe) call p_recv  (s% x, dace% pio, p_tag=i)
        endif
      endif
    end do
  end subroutine random_vector
!------------------------------------------------------------------------------
  subroutine print_vector (x, comment, unit, pe, verbose)
  !-------------
  ! print vector
  !-------------
  type (t_vector)  ,intent(in)           :: x       ! vector to print
  character(len=*) ,intent(in) ,optional :: comment !
  integer          ,intent(in) ,optional :: unit    ! output file unit
  integer          ,intent(in) ,optional :: pe      ! pe to use
  integer          ,intent(in) ,optional :: verbose ! 0: metadata only
                                                    ! 1: segments
                                                    ! 3: elements
    integer                       :: is, i
    integer                       :: lunit, lpe, lverb
    type (t_vector_segm) ,pointer :: s
    real(wp)                      :: mi, ma

    lunit = 6    ;if(present(unit   )) lunit = unit
    lpe   = dace% pio ;if(present(pe     )) lpe   = pe
    lverb = 0    ;if(present(verbose)) lverb = verbose

    call enter_function
    mi = minval(x)
    ma = maxval(x)
    if (dace% pe == lpe) then
      write(lunit,'()')
      write(lunit,'(a)') 'type t_vector:'
      if (present(comment)) write(lunit,'(a,a )') '  ',comment
      write(lunit,'(a,a )') '  name    = ',x% name
      write(lunit,'(a,i5)') '  n       =' ,x% n
      write(lunit,'(a,i5)') '  n_s     =' ,x% n_s
      write(lunit,'(a,i5)') '  alloc_l =' ,x% alloc_l
      write(lunit,'(a,l5)') '  global  =' ,x% global
      write(lunit,*)         ' min,max =' ,mi,ma
      write(lunit,'(a)')    '   is    i       x'
    endif
    if (lverb > 0) then
      do is = 1, x% n_s
        s => x% s(is)
        if (dace% pe == lpe) then
          !------------------------------
          ! on dace% pe==lpe receve and print
          !------------------------------
          if (.not. associated (s% x)) allocate (s% x(s% n))
          if (.not.(x% global.or.dace% pe==s% pe)) call p_recv (s%x, s%pe,p_tag=is)
          if (lverb > 1) then
            do i = 1, s% n
              write (lunit,'(2i5,f15.7)') is, i, s% x(i)
            end do
          else
            write (lunit,'(2i5,2f15.7)') is, s% n, minval(s% x),maxval(s% x)
          endif
          if (.not.(x% global.or.dace% pe==s% pe)) deallocate (s% x)
        else
          !--------------------
          ! on p_pe/=lpe: send
          !--------------------
          if (.not. x% global .and. dace% pe /= lpe) then
            if (s% pe == dace% pe) call p_send  (s% x, lpe, p_tag=is)
          endif
        endif
      end do
    endif
    call delete_storage (x)
    call leave_function
  end subroutine print_vector
!------------------------------------------------------------------------------
  subroutine print_matrix (x, comment, unit, pe, verbose)
  !-------------
  ! print matrix
  !-------------
  type (t_matrix)  ,intent(in)           :: x       ! matrix to print
  character(len=*) ,intent(in) ,optional :: comment !
  integer          ,intent(in) ,optional :: unit    ! output file unit
  integer          ,intent(in) ,optional :: pe      ! pe to use
  integer          ,intent(in) ,optional :: verbose ! 0: metadata only
                                                    ! 1: blocks
                                                    ! 3: elements
    integer                       :: lunit, lpe, lverb
    integer                       :: i,j

    lunit = 6    ;if(present(unit   )) lunit = unit
    lpe   = dace% pio ;if(present(pe     )) lpe   = pe
    lverb = 0    ;if(present(verbose)) lverb = verbose

    call enter_function
    if (dace% pe == lpe) then
      write(lunit,'()')
      write(lunit,'(a,i5)') '  p  e    = ',lpe
      write(lunit,'(a)') 'type t_matrix:'
      if (present(comment)) write(lunit,'(a,a )') '  ',comment
      write(lunit,'(a,a )') '  name    = ',x% name
      write(lunit,'(a,i5)') '  m       =' ,x% m
      write(lunit,'(a,i5)') '  n       =' ,x% n
      write(lunit,'(a,i5)') '  m_b     =' ,x% m_b
      write(lunit,'(a,i5)') '  n_b     =' ,x% n_b
      write(lunit,'(a,i5)') '  alloc_l =' ,x% alloc_l
      write(lunit,'(a,a )') '  associated:'
      write(lunit,'(a,l1)') '  rinfo   :' ,associated (x% rinfo)
      write(lunit,'(a,l1)') '  cinfo   :' ,associated (x% cinfo)
      write(lunit,'(a,l1)') '  b       :' ,associated (x% b)
      if (associated (x% b))&
      write(lunit,'(a,2i5)')'  shape(b)=' ,shape(x% b)
      write(lunit,'(a,l1)') '  next    :' ,associated (x% next)
      write(lunit,'(a,l1)') '  prev    :' ,associated (x% prev)
      if (verbose >= 1) then
        do   j = 1, x% n_b   ! number of column blocks
          do i = 1, x% m_b   ! number of row    blocks
            write (lunit,*)
            write (lunit,'(a,i5)') ' j_b        = ',j
            write (lunit,'(a,i5)') ' i_b        = ',i
            call print (x% b(i,j), unit=lunit, pe=lpe, verbose=lverb)
          end do
        end do
      endif
    endif
    call delete_storage (x)
    call leave_function
  end subroutine print_matrix
!------------------------------------------------------------------------------
  subroutine explore_vector (x)
  type (t_vector)  ,intent(in) :: x ! vector to show
    integer   :: verbose, pe
    character :: c
    verbose = 0
    pe      = dace% pio
    call enter_function
    do
      call print (x, unit=0, verbose=verbose)
      c = ask_menue ('explore vector',      &
                     'v verbosity',         &
                     'e exit explore',      &
                     'x quit program'       )
      select case (c)
      case ('v')
        call ask (verbose, '0:metadata 1:sections 2:elements')
      case ('p')
        do
          call ask (pe, 'pe (processor element) to print')
          if (pe>0 .and. pe<dace% npe) exit
          write(0,*)'pe must be between',0,'and',dace% npe-1
        end do
      case ('e')
        exit
      case ('x')
        call finish ('explore_vector','exit on user request')
      end select
    end do
    call delete_storage (x)
    call leave_function
  end subroutine explore_vector
!------------------------------------------------------------------------------
  recursive subroutine explore_matrix (x)
  type (t_matrix)  ,intent(in) :: x ! matrix to show
    integer   :: verbose, pe, i, j, ij(2)
    character :: c
    verbose = 0
    pe      = dace% pio
    call enter_function
    do
      call print (x, unit=0, verbose=verbose)
      c = ask_menue ('explore matrix',      &
                     'v verbosity',         &
                     'p set p_io',          &
                     'b write block information',&
                     'e exit explore',      &
                     'x quit program'       )
      select case (c)
      case ('v')
        call ask (verbose, '0:metadata 1:sections 2:elements')
      case ('e')
        exit
      case ('x')
        call finish ('explore_matrix','exit on user request')
      case ('b')
        do
          call ask (ij,'block indices')
          i = ij(1); j=ij(2)
          if (i>0 .and. i<=x%m_b .and. j>0 .and. j<=x%n_b) exit
          write(0,*)'i must be between',1,'and',x%m_b
          write(0,*)'j must be between',1,'and',x%n_b
        end do
        call print (x% b(i,j), pe=pe)
        if (x% b(i,j) %repr == NEST .and. associated(x% b(i,j) %nest)) then
          call explore (x% b(i,j) %nest)
        end if
      end select
    end do
    call delete_storage (x)
    call leave_function
  end subroutine explore_matrix
!------------------------------------------------------------------------------
  subroutine reorder_vector (y, x, segm, segm_idx)
  type(t_vector)  ,intent(inout) :: y
  type(t_vector)  ,intent(in)    :: x
  type(t_ivector) ,intent(in)    :: segm
  type(t_ivector) ,intent(in)    :: segm_idx
  !-----------------------------------------------------------------------
  ! Reorders the elements according to the indices given in the
  ! integer-vectors SEGM and SEGM_IDX. The element order of these vectors
  ! is the same as that of the output vextor Y.  SEGM denotes the segment
  ! number and SEGM_IDX the index within the segment in the input array X.
  !-----------------------------------------------------------------------
    call finish('reorder_vector','not yet implemented')
  end subroutine reorder_vector
!------------------------------------------------------------------------------
  subroutine reorder_matrix (y, x, box, box_idx, qual)
  type(t_matrix)  ,intent(inout)        :: y
  type(t_matrix)  ,intent(in)           :: x
  type(t_ivector) ,intent(in)           :: box
  type(t_ivector) ,intent(in)           :: box_idx
  integer         ,intent(in) ,optional :: qual
    call finish('reorder_matrix','not yet implemented')
  end subroutine reorder_matrix
!==============================================================================
  !
  ! The following routines perform the arithmetic operations +,-,* =, on
  ! the block level:
  !
!------------------------------------------------------------------------------
  function matrix_times_vector_b (B, x) result (y)
  !----------------------------
  ! matrix * vector, overloaded
  !----------------------------
  type (t_matrix_block) ,intent(in)    :: B
  real(wp)              ,intent(in)    :: x (:)
  real(wp)                             :: y (b% m)
    y = 0._wp
    call matrix_times_vector_block (y, B, x, t=.false.)
  end function matrix_times_vector_b
!------------------------------------------------------------------------------
  function vector_times_matrix_b (x, B) result (y)
  !----------------------------
  ! vector * matrix, overloaded
  !----------------------------
  real(wp)              ,intent(in)    :: x (:)
  type (t_matrix_block) ,intent(in)    :: B
  real(wp)                             :: y (b% n)
    y = 0._wp
    call matrix_times_vector_block (y, B, x, t=.true.)
  end function vector_times_matrix_b
!------------------------------------------------------------------------------
  subroutine update_matrix_block (A, u, v)
  type (t_matrix_block) ,intent(inout) :: A
  real(wp)              ,intent(in)    :: u (:)
  real(wp)              ,intent(in)    :: v (:)
    select case (A% repr)
    case (CSR)
      call update_csr       (A% packed, A% ja, A% ia, u, v, A%m)
    case (CSC)
      call update_csc       (A% packed, A% ja, A% ia, u, v, A%n)
    case default
      call finish('update_matrix_block',&
                  'invalid matrix representation '//srep(A% repr))
    end select
  end subroutine update_matrix_block
!------------------------------------------------------------------------------
  recursive subroutine matrix_times_vector_block (y, b, x, t)
  !--------------------
  ! y = B   * x + y  or
  ! y = B^t * x + y
  !--------------------
  real(wp)              ,intent(inout) :: y (:)
  type (t_matrix_block) ,intent(in)    :: b
  real(wp)              ,intent(in)    :: x (:)
  logical     ,optional ,intent(in)    :: t ! transpose flag
    integer               :: info
    logical               :: tr
    real(wp) ,allocatable :: tmp(:)
    if (b%m==0 .or. b%n==0) return
    tr = .false.; if (present(t)) tr = t
    select case (b% repr)
    case (ZERO)
      return
    case (IDENT)
      y = y + x
    case (DIAGONAL)
      y = y + b% packed * x
    case (FULL)
FTRACE_BEGIN("matrix_times_vector_block:dgemv")
      if (tr) then
        call DGEMV ('T',b%m, b%n, 1._wp, real(b% full,wp), b%m, x, 1, 1._wp, y, 1)
      else
        call DGEMV ('N',b%m, b%n, 1._wp, real(b% full,wp), b%m, x, 1, 1._wp, y, 1)
      endif
FTRACE_END  ("matrix_times_vector_block:dgemv")
    case (PACKED)
      if (tr) call finish ('matrix_times_vector_block','transposed PACKED')
      call DSPMV (b%tri, b%m, 1._wp, real(b% packed,wp), x, 1, 1._wp, y, 1)
    case (CSR)
FTRACE_BEGIN("matrix_times_vector_block:csr")
      if (tr) then
        if (b% qual == SYMMETRIC) then
          if (associated (b% perm)) then
            call csrperm_times_vector (b% packed, b% ja, b% ia, x, y, b% perm)
          else
            call csr_times_vector (b% packed, b% ja, b% ia, x, y, b%m)
          end if
        else    ! not symmetric
          call csc_times_vector (b% packed, b% ja, b% ia, x, y, b%m)
        end if
      else  ! not tr
        if (associated (b% perm)) then
          call csrperm_times_vector (b% packed, b% ja, b% ia, x, y, b% perm)
        else
          call csr_times_vector (b% packed, b% ja, b% ia, x, y, b%m)
        end if
      endif ! not tr
FTRACE_END  ("matrix_times_vector_block:csr")
    case (CSC)
FTRACE_BEGIN("matrix_times_vector_block:csc")
      if (tr) then
        if (associated (b% perm)) then
          call csrperm_times_vector (b% packed, b% ja, b% ia, x, y, b% perm)
        else
          call csr_times_vector (b% packed, b% ja, b% ia, x, y, b%n)
        end if
      else  ! not tr
        if (b% qual == SYMMETRIC) then
          if (associated (b% perm)) then
            call csrperm_times_vector (b% packed, b% ja, b% ia, x, y, b% perm)
          else
            call csr_times_vector (b% packed, b% ja, b% ia, x, y, b%m)
          end if
        else    ! not symmetric
          call csc_times_vector (b% packed, b% ja, b% ia, x, y, b%n)
        end if
      endif ! not tr
FTRACE_END  ("matrix_times_vector_block:csc")
    case (CHOLES)
      if (tr) call finish ('matrix_times_vector_block','transposed CHOLES')
      allocate (tmp(b% m))
      tmp = x
FTRACE_BEGIN("matrix_times_vector_block:dpptrs")
      call dpptrs (b%tri, b%m, 1, real(b% packed,wp), tmp, b%m, info)
FTRACE_END  ("matrix_times_vector_block:dpptrs")
      if (info/=0) call finish ('matrix_times_vector_block','info/=0')
      y = y + tmp
      deallocate (tmp)
    case (LU)
      if (tr) call finish ('matrix_times_vector_block','transposed LU')
      allocate (tmp(b% m))
      tmp = x
      call dgetrs ('N', b%m, 1, real(b% full,wp), b%m, b% ia, tmp, b%m, info)
      if (info/=0) call finish ('matrix_times_vector_block','info/=0')
      y = y + tmp
      deallocate (tmp)
    case (JAD)
      if (.not. associated (b% perm)) &
           call finish ("matrix_times_vector_block","JAD but missing metadata")
      if (tr) then
        if (b% qual /= SYMMETRIC) &
             call finish ("matrix_times_vector_block","JAD and transpose")
        ! b is symmetric, fallthrough...
      end if
      if (b%perm% jdiag <= 0 .or. size (b%perm% ia) < 1 .or. &
           .not. associated (b%perm% ia)) then
         print *, "jdiag =", b%perm% jdiag
         print *, "size (b%perm% ia) =", size (b%perm% ia)
         print *, "associated (b%perm% ia) =", associated (b%perm% ia)
         call finish ("matrix_times_vector_block","metadate for amuxj")
      end if
      allocate (tmp(b% m))
!FTRACE_BEGIN("matrix_times_vector_block:amuxj")
      call amuxj (b% m, x, tmp, b%perm% jdiag, b% packed, b% ja, b%perm% ia)
      y(b%perm% iperm) = y(b%perm% iperm) + tmp
!FTRACE_END  ("matrix_times_vector_block:amuxj")
      deallocate (tmp)
    case (NEST)
      if (tr) then
        y = y + array_times_matrix (x, b% nest)
      else
        y = y + matrix_times_array (b% nest, x)
      endif
    case default
      call finish('matrix_times_vector_block',&
                  'invalid matrix representation '//srep(b% repr))
    end select
  end subroutine matrix_times_vector_block
!------------------------------------------------------------------------------
  !
  ! Store the upper or lower triangle of symmetric matrix b%full(:,:) in packed
  ! form in b%packed(:)
  !
  ! upper triangle (tri='U'):
  !
  !   b% packed(i + (j-1)*j/2) = b% full(i,j) for 1<=i<=j
  !
  ! lower triangle (tri='L'):
  !
  !   b% packed(i + (j-1)*(2n-j)/2) = b% full(i,j) for j<=i<=n
  !
  ! The packed storage scheme is illustrated by the following example
  ! when N = 4 and tri = 'U':
  !
  ! Two-dimensional storage of the symmetric matrix B%FULL:
  !
  !   b11 b12 b13 b14
  !       b22 b23 b24
  !           b33 b34     (bij = bji)
  !               b44
  !
  ! Packed storage of the upper triangle of B%PACKED:
  !
  ! B%PACKED = [ b11, b12, b22, b13, b23, b33, b14, b24, b34, b44 ]
  !
  subroutine pack_block (b, tri)
  type (t_matrix_block) ,intent(inout) :: b
  character ,optional   ,intent(in)    :: tri
    integer :: n, i, j, k
    integer :: nz           ! Non-zero mirror elements (needed for SX-9)
    integer(kind=i8) :: n_

    !--------------------------------------------
    ! convert to FULL matrix representation first
    !--------------------------------------------
    select case (b% repr)
    case (FULL)                                    ! no prior unpacking
    case (ZERO, IDENT, DIAGONAL, NOT_INIT, PACKED) ! no action at all
      return
    case default                                   ! prior unpacking
      call unpack_block (b)
    end select

    !---------------------------------
    ! convert to PACKED representation
    !---------------------------------
    if (b% repr == FULL) then
      n_ = int(b%n, kind=i8)
      n_ = n_ * (n_+1)/2
      if (n_ > huge(n)) call finish ('pack_block','matrix too large')
      n = int(n_)
      allocate (b% packed(n))
      k = 0
      b% nonzero = 0
      b% tri     = 'L'; if (present(tri)) b% tri     = tri
      select case (b%tri)
      case ('U')
        nz = 0
        do j = 1, b%n
          do i = 1, j
            k=k+1
            b% packed(k) = b% full (i,j)
            if (b% packed(k)/=0._mp) then
              b% nonzero = b% nonzero + 1
!             if(i/=j) b% nonzero = b% nonzero + 1
              if(i/=j) nz = nz + 1
            endif
          end do
        end do
        b% nonzero = b% nonzero + nz
      case ('L')
        nz = 0
        do j = 1, b%n
          do i = j, b%n
            k=k+1
            b% packed(k) = b% full (i,j)
            if (b% packed(k)/=0._mp) then
              b% nonzero = b% nonzero + 1
!             if(i/=j) b% nonzero = b% nonzero + 1
              if(i/=j) nz = nz + 1
            endif
          end do
        end do
        b% nonzero = b% nonzero + nz
      case default
        call finish ('pack_block','invalid triangle: '//b%tri)
      end select
      deallocate (b% full)
      b% repr = PACKED
    else
      call finish ('pack_block','invalid representation: '//srep(b% repr))
    endif
  end subroutine pack_block
!------------------------------------------------------------------------------
  !
  ! Restore the matrix stored in packed form into full form
  !
  ! upper triangle (tri='U'):
  !
  !   b% packed(i + (j-1)*j/2) = b% full(i,j) for 1<=i<=j
  !
  ! lower triangle (tri='L'):
  !
  !   b% packed(i + (j-1)*(2n-j)/2) = b% full(i,j) for j<=i<=n
  !
  subroutine unpack_block (b)
  type (t_matrix_block) ,intent(inout) :: b
    integer :: n, i, j, k
    integer :: nz           ! Non-zero mirror elements (needed for SX-9)
    !----------------------------------
    ! no action for FULL representation
    !----------------------------------
    if (b% repr == FULL) return
    !--------------------------------
    ! unpack for other representation
    !--------------------------------
    allocate (b% full(b%m,b%n))
    if (b% repr == PACKED) then
      n = b%n * (b%n+1) / 2
      if (size(b%packed)/=n) then
        write(0,*)'unpack_block: m,n,n*(n+1)/2,size=',b%m,b%n,n,size(b%packed)
        call finish ('unpack_block','n*(n+1)/2 /= size')
      endif
      k = 0
      b% nonzero = 0
      select case (b%tri)
      case ('U')
        nz = 0
        do j = 1, b%n
!NEC$ ivdep
          do i = 1, j
            k=k+1
            b% full (i,j) = b% packed(k)
            if(i/=j) b% full (j,i) = b% packed(k)
            if (b% packed(k)/=0._mp) then
              b% nonzero = b% nonzero + 1
!             if(i/=j) b% nonzero = b% nonzero + 1
              if(i/=j) nz = nz + 1
            endif
          end do
        end do
        b% nonzero = b% nonzero + nz
      case ('L')
        nz = 0
        do j = 1, b%n
!NEC$ ivdep
          do i = j, b%n
            k=k+1
            b% full (i,j) = b% packed(k)
            if(i/=j) b% full (j,i) = b% packed(k)
            if (b% packed(k)/=0._mp) then
              b% nonzero = b% nonzero + 1
!             if(i/=j) b% nonzero = b% nonzero + 1
              if(i/=j) nz = nz + 1
            endif
          end do
        end do
        b% nonzero = b% nonzero + nz
      case default
        call finish ('unpack_block','invalid triangle: '//b%tri)
      end select
      deallocate (b% packed)
      b% tri  = ''
    else if (b% repr == ZERO) then
      b% full = 0._mp
    else if (b% repr == IDENT) then
      b% full = 0._mp
      do i = 1, b%n
        b% full(i,i) = 1._mp
      enddo
    else if (b% repr == DIAGONAL) then
      b% full = 0._mp
      do i = 1, b%n
        b% full(i,i) = b% packed(i)
      enddo
      deallocate (b% packed)
    else if (b% repr == CSR) then
      b% full = 0._mp
      do i=1,b%m
        k=b% ia (i+1) - 1
        do  j=b% ia (i), k
          b% full(i,b% ja (j)) = b% packed (j)
        end do
      end do
      deallocate (b% packed, b% ia, b% ja)
    else if (b% repr == CSC) then
      b% full = 0._mp
      do i=1,b%n
        k=b% ia (i+1) - 1
        do  j=b% ia (i), k
          b% full(b% ja (j),i) = b% packed (j)
        end do
      end do
      deallocate (b% packed, b% ia, b% ja)
    else
      call finish ('unpack_block','invalid representation: '//srep(b% repr))
    endif
    b% repr = FULL
  end subroutine unpack_block
!------------------------------------------------------------------------------
  function diag_packed (b) result (d)
  type (t_matrix_block) ,intent(in) :: b
  real(wp)                          :: d(b%n)
  !---------------------------------------------
  ! return the diagonal of a packed matrix block
  !---------------------------------------------
    integer :: j, k, n2
    select case (b%tri)
    case ('U')
      k = 0
      do j = 1, b%n
!       k = k + j               ! Recurrence relation for scalar machines
        k = j*(j+1)/2           ! Vectorized index calculation
        d(j) = b% packed(k)
      end do
    case ('L')
      k = 1
      n2 = 2*b%n
      do j = 1, b%n
        k = ((n2+3-j)*j - n2)/2 ! Vectorized index calculation
        d(j) = b% packed(k)
!       k = k + b%n-j+1         ! Recurrence relation for scalar machines
      end do
    end select
  end function diag_packed
!------------------------------------------------------------------------------
  !
  ! Store the upper triangle of symmetric matrix b%full(:,:) in
  ! row compressed storage format in b%packed(:)
  !
  !      The nonzero entries of the coefficient matrix m are stored
  ! row-by-row in the array a.  To identify the individual nonzero
  ! entries in each row, we need to know in which column each entry
  ! lies.  The column indices which correspond to the nonzero entries
  ! of m are stored in the array ja;  i.e., if  a(k) = m(i,j),  then
  ! ja(k) = j.  In addition, we need to know where each row starts and
  ! how long it is.  The index positions in ja and a where the rows of
  ! m begin are stored in the array ia;  i.e., if m(i,j) is the first
  ! nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
  ! ia(i) = k.  Moreover, the index in ja and a of the first location
  ! following the last element in the last row is stored in ia(n+1).
  ! thus, the number of entries in the i-th row is given by
  ! ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
  ! consecutively in
  !         a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
  ! and the corresponding column indices are stored consecutively in
  !         ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
  ! for example, the 5 by 5 matrix
  !             ( 1. 0. 2. 0. 0.)
  !             ( 0. 3. 0. 0. 0.)
  !         m = ( 0. 4. 5. 6. 0.)
  !             ( 0. 0. 0. 7. 0.)
  !             ( 0. 0. 0. 8. 9.)
  ! would be stored as
  !            \ 1  2  3  4  5  6  7  8  9
  !         ---+--------------------------
  !         ia \ 1  3  4  7  8 10
  !         ja \ 1  3  2  2  3  4  4  4  5
  !          a \ 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
  !
  subroutine csr_block (b)
  type (t_matrix_block) ,intent(inout) :: b
    integer :: nm, n0, i, j, k
    if (b% repr == FULL) then
      if (b% nonzero <0) b% nonzero = count (b% full /= 0._mp)
      n0 = b% nonzero
      nm = b% n * b% m
      !
      ! compress if we gain memory
      !
      if (2*n0+min(b%n,b%m)+1<nm) then
        allocate (b% packed (n0))
        allocate (b% ja     (n0))
        allocate (b% ia  (b%m+1))
        k = 1
        do i=1,b%m
          b% ia (i) = k
          do j=1,b%n
            if (b% full(i,j) /= 0._mp) then
              b% packed (k) = b% full(i,j)
              b% ja     (k) = j
              k = k + 1
            endif
          end do
        end do
        b% ia (b%m+1) = k
        deallocate (b% full)
        b% repr = CSR
      endif
    endif
  end subroutine csr_block
!------------------------------------------------------------------------------
  subroutine CSC_block (b)
  type (t_matrix_block) ,intent(inout) :: b
    integer :: nm, n0, i, j, k
    if (b% repr == FULL) then
      if (b% nonzero <0) b% nonzero = count (b% full /= 0._mp)
      n0 = b% nonzero
      nm = b% n * b% m
      !
      ! compress if we gain memory
      !
      if (2*n0+min(b%n,b%m)+1<nm) then
        allocate (b% packed (n0))
        allocate (b% ja     (n0))
        allocate (b% ia  (b%n+1))
        k = 1
        do j=1,b%n
          b% ia (j) = k
          do i=1,b%m
            if (b% full(i,j) /= 0._mp) then
              b% packed (k) = b% full(i,j)
              b% ja     (k) = I
              k = k + 1
            endif
          end do
        end do
        b% ia (b%n+1) = k
        deallocate (b% full)
        b% repr = CSC
      endif
    endif
  end subroutine csc_block
!------------------------------------------------------------------------------
  !
  ! Perform Cholesky factorization
  !
  subroutine cholesky_block (b, tri, info)
  type (t_matrix_block) ,intent(inout) :: b    ! matrix to invert
  character ,optional   ,intent(in)    :: tri  ! upper or lower triangle
  integer   ,optional   ,intent(out)   :: info ! /=0 for error
    integer               :: inf
    real(wp) ,allocatable :: p (:)
    !--------------------------------------------------------
    ! check for correct representation, keep zero size matrix
    !--------------------------------------------------------
    if (present(info)) info = 0
    if (b% n /= b% m)      call finish ('cholesky_block',&
      'm /= n : '//srep(b% repr)//char3(b%m)//' x '//char3(b%n))
    if (b% n == 0) return
    if (b% repr == FULL) call pack_block (b, tri)
    if (b% repr /= PACKED) call finish ('cholesky_block',&
      'repr /= PACKED : '//srep(b% repr)//char3(b%m)//' x '//char3(b%n))
    !-------------------------------
    ! perform Cholesky factorization
    !-------------------------------
    allocate (p(size(b%packed)))
    p = b%packed
    call dpptrf( b%tri, b%m, p, inf )
    if(inf < 0) then
      write(0,*)'cholesky_block: the',-inf,'-th argument had an illegal value.'
      call finish ('cholesky_block','info < 0')
    else if(inf > 0) then
      if (present (info)) then
        info = inf
      else
        write(0,*)'the leading minor of order ',info,'is not positive definite'
        write(0,*)'and the factorization could not be completed.'
        call finish ('cholesky_block','info > 0')
      endif
    endif
    b% repr = CHOLES
    b% packed = p
    deallocate (p)
  end subroutine cholesky_block
!------------------------------------------------------------------------------
  !
  ! Perform LU factorization
  !
  subroutine LU_block (b)
  type (t_matrix_block) ,intent(inout) :: b
    integer  :: info
    real(wp) :: f (b%m,b%n)
    select case (b% repr)
    case (IDENT)
      return
    case (DIAGONAL)
      b% packed = 1._mp / b% packed
      return
    case (PACKED)
      call unpack_block (b)
    end select
    if (b% repr /= FULL) call finish ('LU_block','repr/= FULL')
    if (associated (b% ia)) deallocate(b% ia)
    allocate(b% ia(b% n))
    f = b% full
    call dgetrf(b% m,b% n,f,b% n, b% ia, info)
    if (info/=0) then
      write(0,*) 'LU_block: =',info
      if(info<0) write(0,*) 'the',-info,'-th argument had an illegal value.'
      if(info>0) then
        write(0,*)'U(',info,',',info,') is exactly zero. The factorization'
        write(0,*)'has been completed, but the factor U is exactly singular,'
        write(0,*)'and division by zero will occur if it is used to solve'
        write(0,*)' a system of equations'
      endif
      call finish ('LU','info/=0')
    endif
    b% repr = LU
    b% full = f
  end subroutine LU_block
!------------------------------------------------------------------------------
  subroutine print_matrix_block (b, unit, pe, verbose)
  type (t_matrix_block), intent(in)           :: b
  integer,               intent(in) ,optional :: unit    ! output file unit
  integer,               intent(in) ,optional :: pe
  integer,               intent(in) ,optional :: verbose ! 1: blocks
                                                         ! 3: elements

    integer :: i, j, k, i1, in, lpe, lunit, lverb
    lpe   = dace% pio  ;if (present(pe))      lpe   = pe
    lunit = 6     ;if (present(unit))    lunit = unit
    lverb = 3     ;if (present(verbose)) lverb = verbose
    if (lpe == dace% pe) then

      write (lunit,*)
      write (lunit,*) 'matrix_block:'
      write (lunit,'(a,i8)') ' m          = ',b% m
      write (lunit,'(a,i8)') ' n          = ',b% n
      write (lunit,'(a,i8)') ' nonzero    = ',b% nonzero
      write (lunit,'(a,i8)') ' pe         = ',b% pe
      write (lunit,'(a,i8)') ' repr       = ',b% repr
      write (lunit,'(a,a)')  ' repr       = ',srep(b% repr)
      write (lunit,'(a,a)')  ' tri        = ',b% tri
      write (lunit,'(a,i8)') ' alloc_l    = ',b% alloc_l
      write (lunit,'(a)')    ' associated :'
      write (lunit,'(a,l1)') ' full       : ',associated(b% full)
      write (lunit,'(a,l1)') ' packed     : ',associated(b% packed)
      write (lunit,'(a,l1)') ' ia         : ',associated(b% ia)
      write (lunit,'(a,l1)') ' ja         : ',associated(b% ja)

      select case (b% repr)
      case (FULL)
        write (lunit,'(a,2i5)')' shape      = ',shape (b% full)
        write (lunit,*)         'min,maxval = ',minval(b% full), maxval(b% full)
        if (lverb==3) then
          do   j = 1, b% n   ! number of columns
            do i = 1, b% m   ! number of rows
              write (lunit,'(a,5x,2i5,f10.3)') '   i,j,b_ij = ',i,j,b% full(i,j)
            end do
          end do
        endif
      case (CSR)
        write (lunit,'(a)')    ' repr       = CSR (compressed sparse row)'
        write (lunit,'(a,2i8)')' shape      = ',shape (b% packed)
        write (lunit,*)         'min,maxval = ',minval(b% packed), &
                                                maxval(b% packed)
        if (lverb==3) then
          do i=1,b% m
            do k=b%ia(i),b%ia(i+1)
              write (lunit,'(a,3i5,f10.3)')' k,i,j,b_ij = ',k,i,b%ja(k),b%packed(k)
            end do
          end do
        endif
      case (CSC)
        write (lunit,'(a)')    ' repr       = CSC (compressed sparse column)'
        write (lunit,'(a,2i8)')' shape      = ',shape (b% packed)
        write (lunit,*)         'min,maxval = ',minval(b% packed), &
                                                maxval(b% packed)
        if (lverb==3) then
          do i=1,b% m
            do k=b%ia(i),b%ia(i+1)
              write (lunit,'(a,3i5,f10.3)')' k,i,j,b_ij = ',k,b%ja(k),i,b%packed(k)
            end do
          end do
        endif
      case (PACKED)
        write (lunit,'(a)')    ' repr       = PACKED ('//b%tri//' triangle only)'
        write (lunit,'(a,2i5)')' shape      = ',shape (b% packed)
        write (lunit,*)         'min,maxval = ',minval(b% packed), maxval(b% packed)
        if (lverb==3) then
          k = 0
          do j = 1, b%n
            select case (b%tri)
            case ('U')
              i1 = 1; in = j
            case ('L')
              i1 = j; in = b%n
            end select
            do i = i1, in
              k=k+1
              if (b% packed(k)/=0._mp) then
                write (lunit,'(a,3i5,f10.3)') ' k,i,j,b_ij = ',k,i,j,b%packed(k)
              endif
            end do
          end do
        endif
      case default
      end select
    endif
  end subroutine print_matrix_block
!------------------------------------------------------------------------------
  function diag_matrix_block (b) result (diag)
  type(t_matrix_block), intent(in) :: b
  real(wp)                         :: diag (b% n)

    integer :: i, k
    if (b% n /= b% m) call finish ('diag_matrix_block','matrix is not square')
    select case (b% repr)
    case (ZERO)
      diag = 0._wp
    case (IDENT)
      diag = 1._wp
    case (DIAGONAL)
      diag = b% packed
    case (PACKED)
      diag = diag_packed (b)
    case (FULL)
      do i = 1, b% n
        diag(i) = b% full(i,i)
      end do
    case (CSR, CSC)
      do i=1, b% n
        diag(i) = 0._wp
        do k=b% ia(i), b% ia(i+1)-1
          if (b% ja(k) == i) diag(i) = b% packed(k)
        end do
      end do
    case default
      if (lbound(srep,1) <= b% repr .and. b% repr <= ubound(srep,1)) &
        call finish ('diag_matrix_block','not implemented for '//srep(b% repr))
      call finish   ('diag_matrix_block','unknown representation')
    end select

  end function diag_matrix_block
!==============================================================================
  !
  ! Low level routines for matrix - vector multiplication
  ! (compressed sparse row/column storage scheme)
  !
  subroutine csr_times_vector (a, ja, ia, x, y, m)
  real(mp)   ,intent(in)    ::  a (:) ! coefficients   of matrix A
  integer(ip),intent(in)    :: ja (:) ! column indices of matrix A
  integer    ,intent(in)    :: ia (:) ! indices to a,ja for row indices
  real(wp)   ,intent(in)    ::  x (:) ! right hand side
  real(wp)   ,intent(inout) ::  y (:) ! left hand side
  integer    ,intent(in)    ::  m     ! number of rows
#ifdef HAVE_F2008_CONTIGUOUS
  contiguous :: a, ja, ia
#endif
    integer :: i
    integer :: k1, k2
!$omp parallel do private(i,k1,k2) firstprivate(m) schedule(auto)
    do i=1,m
      k1 = ia(i)
      k2 = ia(i+1)-1
      y(i) = y(i) + sum (a(k1:k2) * x(ja(k1:k2)))
    end do
!$omp end parallel do
  end subroutine csr_times_vector
!------------------------------------------------------------------------------
  subroutine csc_times_vector (a, ja, ia, x, y, n)
  real(mp)   ,intent(in)    ::  a (:) ! coefficients of matrix A
  integer(ip),intent(in)    :: ja (:) ! row indices  of matrix A
  integer    ,intent(in)    :: ia (:) ! indices to a,ia for column indices
  real(wp)   ,intent(in)    ::  x (:) ! right hand side
  real(wp)   ,intent(inout) ::  y (:) ! left hand side
  integer    ,intent(in)    ::  n     ! number of columns
#ifdef HAVE_F2008_CONTIGUOUS
  contiguous :: a, ja, ia
#endif
    integer :: i, j, k
    do j=1,n                    ! Outer loop j: columns of A
!NEC$ loop_count_test
!NEC$ ivdep
!DIR$ IVDEP
      do k = ia(j), ia(j+1)-1   ! Inner loop i: rows of (sparse) A
        i = ja(k)               ! (the i's are distinct for different k's)
        y(i) = y(i) + a(k) * x(j)
      end do
    end do
  end subroutine csc_times_vector
!------------------------------------------------------------------------------
  subroutine update_csr (a, ja, ia, u, v, m)
  real(mp)   ,intent(inout) ::  a (:) ! coefficients   of matrix A
  integer(ip),intent(in)    :: ja (:) ! column indices of matrix A
  integer    ,intent(in)    :: ia (:) ! indices to a,ja for row indices
  real(wp)   ,intent(in)    ::  u (:) ! right hand side
  real(wp)   ,intent(in)    ::  v (:) ! left hand side
  integer    ,intent(in)    ::  m     ! number of rows
    integer :: i
    integer :: k1, k2
    do i=1,m
      k1 = ia(i)
      k2 = ia(i+1)-1
      a(k1:k2) = a(k1:k2) + u(i) * v(ja(k1:k2))
    end do
  end subroutine update_csr
!------------------------------------------------------------------------------
  subroutine update_csc (a, ja, ia, u, v, n)
  real(mp)   ,intent(inout) ::  a (:) ! coefficients   of matrix A
  integer(ip),intent(in)    :: ja (:) ! column indices of matrix A
  integer    ,intent(in)    :: ia (:) ! indices to a,ja for row indices
  real(wp)   ,intent(in)    ::  u (:) ! right hand side
  real(wp)   ,intent(in)    ::  v (:) ! left hand side
  integer    ,intent(in)    ::  n     ! number of columns
    integer :: i, j, k
    do j=1,n
!NEC$ ivdep
      do k = ia(j), ia(j+1)-1
        i = ja(k)
        a(k) = a(k) + u(i) * v(j)
      end do
    end do
  end subroutine update_csc
!==============================================================================
  subroutine get_row (row, c, i, pe)
  !=============================
  ! return one row of the matrix
  !=============================
  real(wp)          ,intent(out) :: row(:) ! ith row of c
  type (t_matrix)   ,intent(in)  :: c      ! matrix
  integer           ,intent(in)  :: i      ! row index
  integer ,optional ,intent(in)  :: pe     ! return row on processor PE only

    integer :: ib, jb ! indices of matrix blocks
    integer :: pe_row ! index of processor which holds the requested row
    integer :: i0, in ! bounds for rows stored in current matrix block
    integer :: ii     ! requested row index within current matrix block
    integer :: k      ! offset for columns stored in current matrix block
    integer :: j      ! column index in current matrix block
    integer :: l      ! index to coefficients in sparse representation
    type (t_matrix_block) ,pointer :: b

    !--------------------------------
    ! loop over rows of matrix blocks
    !--------------------------------
    i0 = 0
    do ib = 1, c% m_b
      in = i0 + c% rinfo% b(ib)% n
      if (i>i0 .and. i<=in) then
        !-------------------------------------
        ! row is located within this block row
        !-------------------------------------
        pe_row = c% rinfo% b(ib)% pe
        if (pe_row == dace% pe) then
          !---------------------------------------
          ! block row is located on this processor
          !---------------------------------------
          ii = i - i0
          row = 0._wp
          !------------------------
          ! loop over block columns
          !------------------------
          k   = 0
          do jb = 1, c% n_b
            b => c% b(ib,jb)
            if (b% repr == MIRROR) then
              !------------------------------------------------
              ! extract from block b(jb,ib) instead of b(ib,jb)
              !------------------------------------------------
              b => c% b(jb,ib)
              select case (b% repr)
              case (ZERO)      ! all matrix elements are zero
              case (CSR)       ! compressed sparse row representation
                do j = 1, b% m
                  do l= b% ia(j), b% ia(j+1)-1
                    if (b% ja (l) == ii) row (k+j) = b% packed(l)
                  end do
                end do
              case (CSC)       ! compressed sparse column representation
                do l = b% ia(ii), b% ia(ii+1)-1
                  row (k+b% ja(l)) = b% packed(l)
                end do
              case (FULL)      ! full n x m representation
                row(k+1 : k+b% m) = b% full (:,ii)
              case default
                call finish ('get_row','invalid matrix block representation')
              end select
              k = k + b% m
            else
              !----------------------------
              ! extract from block b(ib,jb)
              !----------------------------
              select case (b% repr)
              case (ZERO)        ! all matrix elements are zero
              case (IDENT)       ! identity matrix
                row (k+ii) = 1._wp
              case (DIAGONAL)    ! diagonal matrix
                row (k+ii) = b% packed(ii)
              case (PACKED)
                select case (b% tri)
                case ('U')       ! only upper triangle is stored
                  do j=1,ii
                    row (k+j) = b% packed(j +(ii-1)*ii/2)
                  end do
                  do j=ii+1,b% n
                    row (k+j) = b% packed(ii+( j-1)* j/2)
                  end do
                case ('L')       ! only lower triangle is stored
                  do j=1,ii
                    row (k+j) = b% packed(ii + (j-1)*(2*b%n-j)/2)
                  end do
                  do j=ii+1,b% n
                    row (k+j) = b% packed(j + (ii-1)*(2*b%n-ii)/2)
                  end do
                end select
              case (CSR)         ! compressed sparse row representation
                do l = b% ia(ii), b% ia(ii+1)-1
                  row (k+b% ja(l)) = b% packed(l)
                end do
              case (CSC)         ! compressed sparse column representation
                do j = 1, b% n
                  do l= b% ia(j), b% ia(j+1)-1
                    if (b% ja (l) == ii) row (k+j) = b% packed(l)
                  end do
                end do
              case (FULL)        ! full n x m representation
                row(k+1 : k+b% n) = b% full (ii,:)
              case default
                call finish ('get_row','invalid matrix block representation')
              end select
              k = k + b% n
            end if
          end do
        endif
        exit
      endif
      i0 = in
    end do
    !-------------------------
    ! send, receive, broadcast
    !-------------------------
    if (.not. present(pe)) then
      call p_bcast (row, pe_row)
    else if (pe /= pe_row) then
      if (dace% pe == pe_row) call p_send (row, pe    , 1)
      if (dace% pe == pe)     call p_recv (row, pe_row, 1)
    endif
  end subroutine get_row
!------------------------------------------------------------------------------
  subroutine get_row_sparse (row, idx, n, c, i, pe, d_first)
  !=============================
  ! return one row of the matrix
  !=============================
  real(wp)          ,intent(out) :: row(:)  ! ith row of c (nonzero elements)
  integer           ,intent(out) :: idx(:)  ! column index
  integer           ,intent(out) :: n       ! number of nonzero elements
  type (t_matrix)   ,intent(in)  :: c       ! matrix
  integer           ,intent(in)  :: i       ! row index
  integer ,optional ,intent(in)  :: pe      ! return row on processor PE only
  logical ,optional ,intent(in)  :: d_first ! put diagonal in first place

    integer :: ib, jb ! indices of matrix blocks
    integer :: pe_row ! index of processor which holds the requested row
    integer :: i0, in ! bounds for rows stored in current matrix block
    integer :: ii     ! requested row index within current matrix block
    integer :: k      ! offset for columns stored in current matrix block
    integer :: j      ! column index in current matrix block
    integer :: l      ! index to coefficients in sparse representation
    integer :: id     ! diagonal index
    logical :: d1     ! put diagonal in first place
    integer :: ier    ! 'sort' error return flag
    integer ,allocatable           :: iperm (:) ! index vector for permutation
    type (t_matrix_block) ,pointer :: b         ! pointer to matrix block

    !--------------------------------
    ! loop over rows of matrix blocks
    !--------------------------------
    d1 = .false.; if (present(d_first)) d1 = d_first
    id = 0;       if (d1)               id = i
    i0 = 0
    n  = 0
    do ib = 1, c% m_b
      in = i0 + c% rinfo% b(ib)% n
      if (i>i0 .and. i<=in) then
        !-------------------------------------
        ! row is located within this block row
        !-------------------------------------
        pe_row = c% rinfo% b(ib)% pe
        if (pe_row == dace% pe) then
          !---------------------------------------
          ! block row is located on this processor
          !---------------------------------------
          ii = i - i0
          !------------------------
          ! loop over block columns
          !------------------------
          k   = 0
          do jb = 1, c% n_b
            b => c% b(ib,jb)
            if (b% repr == MIRROR) then
              !------------------------------------------------
              ! extract from block b(jb,ib) instead of b(ib,jb)
              !------------------------------------------------
              b => c% b(jb,ib)
              select case (b% repr)
              case (ZERO)      ! all matrix elements are zero
              case (CSR)       ! compressed sparse row representation
                do j = 1, b% m
                  do l= b% ia(j), b% ia(j+1)-1
                    if (b% ja (l) == ii .and. b% packed(l)/=0._mp) then
                      n = n + 1
                      idx (n) = k+j
                      row (n) = b% packed(l)
                      if (idx(n)==id) idx (n) = 0
                    endif
                  end do
                end do
              case (CSC)       ! compressed sparse column representation
                do l = b% ia(ii), b% ia(ii+1)-1
                  if (b% packed(l)/=0._mp) then
                    n = n + 1
                    idx (n) = k+b% ja(l)
                    row (n) = b% packed(l)
                    if (idx(n)==id) idx (n) = 0
                  endif
                end do
              case (FULL)      ! full n x m representation
                do l = 1,b% m
                  if (b% full (l,ii) /= 0._mp) then
                    n = n + 1
                    idx (n) = k+l
                    row (n) = b% full (l,ii)
                    if (idx(n)==id) idx (n) = 0
                  endif
                end do
              case default
                call finish ('get_row_sparse',&
                             'invalid matrix block representation')
              end select
              k = k + b% m
            else
              !----------------------------
              ! extract from block b(ib,jb)
              !----------------------------
              select case (b% repr)
              case (ZERO)        ! all matrix elements are zero
              case (IDENT)       ! identity matrix
                n = n + 1
                idx (n) = k+ii
                row (n) = 1._wp
                if (idx(n)==id) idx (n) = 0
              case (DIAGONAL)    ! diagonal matrix
                if (b% packed(ii) /= 0._mp) then
                  n = n + 1
                  idx (n) = k+ii
                  row (n) = b% packed(ii)
                  if (idx(n)==id) idx (n) = 0
                endif
              case (PACKED)
                select case (b% tri)
                case ('U')       ! only upper triangle is stored
                  do j=1,ii
                    if (b% packed(j +(ii-1)*ii/2) /= 0._mp) then
                      n = n + 1
                      idx (n) = k+j
                      row (n) = b% packed(j +(ii-1)*ii/2)
                      if (idx(n)==id) idx (n) = 0
                    endif
                  end do
                  do j=ii+1,b% n
                    if (b% packed(ii+( j-1)* j/2) /= 0._mp) then
                      n = n + 1
                      idx (n) = k+j
                      row (n) = b% packed(ii+( j-1)* j/2)
                      if (idx(n)==id) idx (n) = 0
                    endif
                  end do
                case ('L')       ! only lower triangle is stored
                  do j=1,ii
                    if (b% packed(ii + (j-1)*(2*b%n-j)/2) /= 0._mp) then
                      n = n + 1
                      idx (n) = k+j
                      row (n) = b% packed(ii + (j-1)*(2*b%n-j)/2)
                      if (idx(n)==id) idx (n) = 0
                    endif
                  end do
                  do j=ii+1,b% n
                    if (b% packed(j + (ii-1)*(2*b%n-ii)/2) /= 0._mp) then
                      n = n + 1
                      idx (n) = k+j
                      row (n) = b% packed(j + (ii-1)*(2*b%n-ii)/2)
                      if (idx(n)==id) idx (n) = 0
                    endif
                  end do
                end select
              case (CSR)         ! compressed sparse row representation
                do l = b% ia(ii), b% ia(ii+1)-1
                  if (b% packed(l) /= 0._mp) then
                    n = n + 1
                    idx (n) = k+b% ja(l)
                    row (n) = b% packed(l)
                    if (idx(n)==id) idx (n) = 0
                  endif
                end do
              case (CSC)         ! compressed sparse column representation
                do j = 1, b% n
                  do l= b% ia(j), b% ia(j+1)-1
                    if (b% ja (l) == ii .and. b% packed(l) /= 0._mp) then
                      n = n + 1
                      idx (n) = k+j
                      row (n) = b% packed(l)
                      if (idx(n)==id) idx (n) = 0
                    endif
                  end do
                end do
              case (FULL)        ! full n x m representation
                do l = 1,b% n
                  if (b% full (ii,l) /= 0._mp) then
                    n = n + 1
                    idx (n) = k+l
                    row (n) = b% full (ii,l)
                    if (idx(n)==id) idx (n) = 0
                  endif
                end do
              case default
                call finish ('get_row_sparse',&
                             'invalid matrix block representation')
              end select
              k = k + b% n
            end if
          end do
        endif
        exit
      endif
      i0 = in
    end do
    !-----
    ! sort
    !-----
    allocate (iperm(n))
    call sort (idx(1:n), iperm, 1, ier)
    idx(1:n) = idx(iperm)
    row(1:n) = row(iperm)
    deallocate (iperm)
    !------------------------------
    ! put diagonal into first place
    !------------------------------
    if (d1) then
      if (idx(1) /= 0) then
        idx (2:n+1) = idx (1:n)
        row (2:n+1) = row (1:n)
        row (1)     = 0._wp
      endif
      idx(1) = id
    endif
    !-------------------------
    ! send, receive, broadcast
    !-------------------------
    if (.not. present(pe)) then
      call  p_bcast (      n , pe_row)
      call  p_bcast (row(1:n), pe_row)
      call  p_bcast (idx(1:n), pe_row)
    else if (pe /= pe_row) then
      if (dace% pe == pe_row) then
        call p_send (      n , pe    , 1)
        call p_send (row(1:n), pe    , 1)
        call p_send (idx(1:n), pe    , 1)
      endif
      if (dace% pe == pe) then
        call p_recv (      n , pe_row, 1)
        call p_recv (row(1:n), pe_row, 1)
        call p_recv (idx(1:n), pe_row, 1)
      endif
    endif
  end subroutine get_row_sparse
!==============================================================================
! routine for memory usage diagnostics
!------------------------------------------------------------------------------
  subroutine dec_matrix_mem (bytes_used, count, tmu, comment)
  integer(i8)      ,intent(out) :: bytes_used
  integer          ,intent(in)  :: count
  integer          ,intent(in)  :: tmu
  character(len=*) ,intent(in)  :: comment

    type(t_matrix)       ,pointer :: m
    type(t_matrix_block) ,pointer :: b
    type(t_vector)       ,pointer :: v
    type(t_vector_segm)  ,pointer :: s
    integer                       :: elements, sum_elements
    integer                       :: bytes,    sum_bytes
    integer                       :: i, j
    character(len=*), parameter   :: fmt = '(a8,1x,a8,1x,i4    ,1x,i10,1x,i10)'
    character(len=*), parameter   :: fm2 = '(a8,1x,a8,1x," all",1x,i10,1x,i10)'
    !
    ! write header
    !
    call add_line_pio ('')
    call add_line_pio ('memory usage')
    call add_line_pio ('type     name       pe   elements      bytes')
    !                   12345678 12345678 1234 1234567890 1234567890
    !
    ! determine memory usage of matrices
    !
    m => first_m
    sum_elements = 0
    sum_bytes    = 0
    do
      elements = 0
      bytes    = 0
      if (.not.associated (m)) exit
      do   i = 1, m% m_b
        do j = 1, m% n_b
          b => m% b(i,j)
          if (associated (b%full)) then
            elements = elements + size (b%full)
            bytes    = bytes    + size (b%full) * size(transfer(1._mp,(/' '/)))
          endif
          if (associated (b%packed)) then
            elements = elements + size (b%packed)
            bytes    = bytes    + size (b%packed)*size(transfer(1._mp,(/' '/)))
          endif
          if (associated (b%ia)) then
            bytes    = bytes    + size (b%ia) * size(transfer(1,(/' '/)))
          endif
          if (associated (b%ja)) then
            bytes    = bytes    + size (b%ja) * size(transfer(1,(/' '/)))
          endif
        end do
      end do
      !
      ! write report for each matrix
      !
      if (tmu >= 3) then
        call nextline
        write (oline(iol),fmt) 't_matrix', m% name, dace% pe, elements, bytes
      endif
      sum_elements = sum_elements + elements
      sum_bytes    = sum_bytes    + bytes
      elements = p_sum (elements)
      bytes    = p_sum (bytes)
      if (tmu >= 3 .and. dace% npe > 1) then
        call flush_buf
        if (dace% lpio) then
          call nextline
          write (oline(iol),fm2) 't_matrix',m% name, elements, bytes
        endif
      endif
      m => m% next
    end do
    !
    ! write summary for this PE
    !
    if (tmu >= 3) then
      call nextline
      write (oline(iol),fmt)'t_matrix','--all---',dace% pe, sum_elements, sum_bytes
    endif
    !
    ! write summary for all PEs
    !
    bytes_used   = sum_bytes
    sum_elements = p_sum (sum_elements)
    sum_bytes    = p_sum (sum_bytes)
    if (tmu >= 3 .and. dace% npe > 1) then
      call flush_buf
      if (dace% lpio) then
        call nextline
        write (oline(iol),fm2) 't_matrix','--all---',sum_elements, sum_bytes
      endif
    endif
    call flush_buf
    !
    ! determine memory usage of vectors
    !
    v => first_v
    sum_elements = 0
    sum_bytes    = 0
    do
      elements = 0
      bytes    = 0
      if (.not.associated (v)) exit
      do   i = 1, v% n_s
        s => v% s(i)
        if (associated (s% x)) then
          elements = elements + size (s% x)
          bytes    = bytes    + size (s% x) * size(transfer(1._wp,(/' '/)))
        endif
      end do
      !
      ! write report for each vector
      !
!     if (tmu >= 3) then
!       call nextline
!       write (oline(iol),fmt) 't_vector', v% name, dace% pe, elements, bytes
!     endif
      sum_elements = sum_elements + elements
      sum_bytes    = sum_bytes    + bytes
      elements = p_sum (elements)
      bytes    = p_sum (bytes)
!     if (dace% npe > 1) then
      if (tmu >= 3) then
        call flush_buf
        if (dace% lpio) then
          call nextline
          write (oline(iol),fm2) 't_vector',v% name, elements, bytes
        endif
      endif
!     endif
      v => v% next
    end do
    !
    ! write summary for the PE
    !
    if (tmu >= 3) then
      call nextline
      write (oline(iol),fmt)'t_vector','--all---',dace% pe, sum_elements, sum_bytes
    endif
    !
    ! write summary for all PEs
    !
    bytes_used   = bytes_used + sum_bytes
    sum_elements = p_sum (sum_elements)
    sum_bytes    = p_sum (sum_bytes)
    if (tmu >= 3 .and. dace% npe > 1) then
    call flush_buf
      if (dace% lpio) then
        call nextline
        write (oline(iol),fm2) 't_vector','--all---',sum_elements, sum_bytes
      endif
    endif
    call add_line_pio('')
    call flush_buf

  end subroutine dec_matrix_mem
!==============================================================================
! Communication routines for matrix blocks
! (send receive bcast)
!------------------------------------------------------------------------------
  subroutine send_matrix_block (mb, dest, tag, comm)
  type (t_matrix_block) ,intent(in) :: mb
  integer               ,intent(in) :: dest
  integer   ,optional   ,intent(in) :: tag
  integer   ,optional   ,intent(in) :: comm

    integer :: ltag, lcom
    integer :: count = 0

    call enter_function
    ltag = 1          ;if(present(tag))  ltag = tag
    lcom = dace% comm ;if(present(comm)) lcom = comm

    if (count == 0) count = size(transfer(mb,(/' '/)))

    call p_send_derivedtype (mb, count, dest, ltag, lcom)

    if (associated (mb% full))   call p_send (mb% full,  dest, ltag, comm=lcom)
    if (associated (mb% packed)) call p_send (mb% packed,dest, ltag, comm=lcom)
    if (associated (mb% ia))     call p_send (mb% ia,    dest, ltag, comm=lcom)
    if (associated (mb% ja))     call p_send (mb% ja,    dest, ltag, comm=lcom)
    if (associated (mb% nest))   call finish ('send_matrix_block','repr=NEST')

!write (0,*) "send_matrix_block:", dace% pe, count, mb% repr, associated (mb%perm)
    if (associated (mb% perm))   call p_send (mb% perm,  dest, ltag, comm=lcom)

    call delete_storage (mb)
    call leave_function

  end subroutine send_matrix_block
!------------------------------------------------------------------------------
  subroutine recv_matrix_block (mb, source, tag, comm)
  type (t_matrix_block) ,intent(inout) :: mb
  integer               ,intent(in)    :: source
  integer   ,optional   ,intent(in)    :: tag
  integer   ,optional   ,intent(in)    :: comm

    integer :: ltag, lcom
    integer :: count = 0

    ltag = 1          ;if(present(tag))  ltag = tag
    lcom = dace% comm ;if(present(comm)) lcom = comm

    if (count == 0) count = size(transfer(mb,(/' '/)))

    call deallocate (mb)
    call p_recv_derivedtype (mb, count, source, ltag, lcom)
    call allocate_block (mb, mb% repr, ns=size(mb% packed), &
                         perm=associated (mb% perm), dealloc=.false.)

    if(associated(mb% full))   call p_recv (mb% full,  source, ltag, comm=lcom)
    if(associated(mb% packed)) call p_recv (mb% packed,source, ltag, comm=lcom)
    if(associated(mb% ia))     call p_recv (mb% ia,    source, ltag, comm=lcom)
    if(associated(mb% ja))     call p_recv (mb% ja,    source, ltag, comm=lcom)
    if(associated(mb% nest))   call finish ('recv_matrix_block','repr=NEST')

!write (0,*) "recv_matrix_block:", dace% pe, count, mb% repr, associated (mb%perm)
    if (associated (mb% perm)) call p_recv (mb% perm,  source, ltag, comm=lcom)

  end subroutine recv_matrix_block
!------------------------------------------------------------------------------
  subroutine bcast_matrix_block (mb, source, comm)
    type (t_matrix_block) ,intent(inout) :: mb
    integer               ,intent(in)    :: source
    integer   ,optional   ,intent(in)    :: comm

    integer :: lcom
    integer :: count = 0

    lcom = dace% comm ;if(present(comm)) lcom = comm

    if (count == 0) count = size(transfer(mb,(/' '/)))

    if (dace% pe /= source) call deallocate (mb)
    call p_bcast_derivedtype (mb, count, source, lcom)
    if (dace% pe /= source) call allocate_block (mb, mb% repr, ns=mb% nonzero, &
                                                 perm=associated (mb% perm),   &
                                                 dealloc=.false.)

    if (associated (mb% full))   call p_bcast (mb% full,   source, lcom)
    if (associated (mb% packed)) call p_bcast (mb% packed, source, lcom)
    if (associated (mb% ia))     call p_bcast (mb% ia,     source, lcom)
    if (associated (mb% ja))     call p_bcast (mb% ja,     source, lcom)
    if (associated (mb% nest))   call finish ('bcast_matrix_block','repr=NEST')

!write (0,*) "bcast_matrix_block:", dace% pe, count, mb% repr, associated (mb%perm)
    if (associated (mb% perm))   call p_bcast (mb% perm,   source, lcom)

  end subroutine bcast_matrix_block
!------------------------------------------------------------------------------
  subroutine p_sum_vector (x)
  type (t_vector) ,intent(inout) :: x
    integer :: i
    do i=1, x% n_s
      if (.not.associated(x% s(i)% x)) &
        call finish ('p_sum_vector','vector is not globally allocated')
      x% s(i)% x = p_sum (x% s(i)% x)
    end do
  end subroutine p_sum_vector
!------------------------------------------------------------------------------
  subroutine put_seed (info)
  type (t_dec_info) ,intent(inout) :: info   ! decomposition info
    integer                        :: i
    do i= 1,info% n_b
      call construct_stream (info% b(i)% seed, use=info% b(i)% pe==dace% pe)
    end do
  end subroutine put_seed
!==============================================================================
  subroutine set_perm_block (B)
    type (t_matrix_block), intent(inout) :: B

    call destruct_perm (B)

    select case (B% repr)
    case (CSR)
       allocate (B% perm)
       call set_perm_csr (B% m, B% ia, B% perm)
    case (CSC)
       call finish ("set_perm_block", "conversion from CSC not implemented")
    case default
       write (0,*) "set_perm_block: B% repr =", B% repr
       call finish ("set_perm_block", "conversion not possible")
    end select

  end subroutine set_perm_block
  !----------------------------------------------------------------------------
  subroutine set_perm_csr (m, ia, perm)
    !------------------------------
    ! Setup permutation information
    !------------------------------
    integer, intent(in) :: m            ! Number of rows in matrix
    integer, intent(in) :: ia(:)        ! Row indices in CSR representation
    type(t_perm)        :: perm         ! Permutation info

    integer :: i, ipos, istart, nz, maxnz, igroup

    integer, allocatable :: nz_in_row(:)
    integer, allocatable :: rows_in_bucket(:)
    integer, allocatable :: ipnz(:)

    if (size (ia) /= m+1) then
       write (0,*) "m, size (ia) =", m, size (ia)
       call finish ("set_perm_csr","inconsistent metadata")
    end if
    !----------------------------------------------------------------------
    ! Count nonzeros in each row, determine maximum of nonzeros in any row.
    !----------------------------------------------------------------------
    allocate(nz_in_row(m))
    nz_in_row(1:m) = ia(2:m+1)-ia(1:m)
    maxnz = maxval (nz_in_row(1:m))
    perm% m     = m
    perm% maxnz = maxnz
    !-----------------------------------------------------------------
    ! Count rows with same number of nonzeros.  Handle empty rows too.
    !-----------------------------------------------------------------
    allocate(rows_in_bucket(0:maxnz))
    rows_in_bucket(:) = 0
!NEC$ novector
    do i = 1, m
       nz = nz_in_row(i)
       rows_in_bucket(nz) = rows_in_bucket(nz)+1
    enddo
    !----------------------------------------------------
    ! Determine grouping information:
    ! nzgroup: count of nonzeros per row in current group
    ! xgroup:  index to start of group
    !----------------------------------------------------
    perm% ngroup = count (rows_in_bucket(:) > 0)
    allocate(perm% nzgroup(perm% ngroup))
    allocate(perm% xgroup (perm% ngroup+1))
    perm% xgroup = 0
    istart = 1
    igroup = 1
!NEC$ novector
    do i = 0, maxnz
       if (rows_in_bucket(i) > 0) then
          perm% nzgroup(igroup) = i
          perm% xgroup (igroup) = istart
          igroup = igroup + 1
          istart = istart + rows_in_bucket(i)
       endif
    enddo
    perm% xgroup(igroup) = istart

    if (igroup /= perm% ngroup+1 .or.istart /= m+1) then
       print *, igroup, perm% ngroup, istart, m
       call finish ("set_perm_csr","inconsistent grouping information")
    end if
    !----------------------------
    ! Setup of permutation vector
    !----------------------------
    allocate(perm% iperm  (m))
    allocate(ipnz(0:maxnz))
    ipnz(0) = 1
!NEC$ novector
    do i = 1, maxnz
       ipnz(i) = ipnz(i-1) + rows_in_bucket(i-1)
    enddo
!NEC$ novector
    do i = 1, m
       nz = nz_in_row(i)
       ipos = ipnz(nz)
       perm% iperm(ipos) = i
       ipnz(nz) = ipnz(nz)+1
    enddo
    deallocate (ipnz)

    deallocate (rows_in_bucket, nz_in_row)
  end subroutine set_perm_csr
!==============================================================================
  subroutine csrperm_times_vector (a, ja, ia, x, y, perm)
    !------------------------------------------------------------
    ! Sparse matrix multiply (CSR) using a given row permutation.
    ! Fortran95 version by Luis Kornblueh.
    !
    ! [1] Eduardo F. D'Azevedo, Mark R. Fahey, Richard T. Mills,
    !     "Vectorized sparse matrix multiply for compressed
    !      row storage format", ICCS 2005,
    !     Lecture Notes in Computer Science, 3514 (2005) 99-106.
    !------------------------------------------------------------
    real(mp)    ,intent(in)    ::  a (:) ! coefficients   of matrix A
    integer(ip) ,intent(in)    :: ja (:) ! column indices of matrix A
    integer     ,intent(in)    :: ia (:) ! indices to a,ja for row indices
    real(wp)    ,intent(in)    ::  x (:) ! right hand side
    real(wp)    ,intent(inout) ::  y (:) ! left hand side
    type(t_perm),intent(in)    :: perm   ! Permutation info

    integer :: i, j, istart, iend, jstart, jend, igroup, nz
    integer :: ipos, isize, iold

    integer, parameter :: NB = 512      ! Block size
    integer  :: ix(NB)          ! Renamed from ip to avoid collision with kind
    real(wp) :: yp(NB)          ! Permuted result

    if (perm% ngroup < 0) call finish ("csrperm_times_vector","perm not set!")
    if (perm% maxnz + 1       <  perm% ngroup .or. &
         size (perm% iperm  ) /= perm% m      .or. &
         size (perm% nzgroup) /= perm% ngroup .or. &
         size (perm% xgroup)  /= perm% ngroup+1) then
       print *, "csrperm_times_vector:"
       print *, perm% m, perm% maxnz, perm% ngroup
       print *, size (perm% iperm), size (perm% nzgroup), size (perm% xgroup)
       print *, "xgroup:", perm% xgroup(perm% ngroup+1), maxval (perm% xgroup)
       call finish ("csrperm_times_vector","bad metadata")
    end if

    do igroup = 1, perm% ngroup
       jstart = perm% xgroup (igroup)
       jend   = perm% xgroup (igroup+1) - 1
       nz     = perm% nzgroup(igroup)
       !--------------------------------------------------------------------
       ! Handle special cases where the number of nonzeros per row is 0 or 1
       !--------------------------------------------------------------------
       if (nz == 0) then
          y(perm% iperm(jstart:jend)) = 0.0_wp
       else if (nz == 1) then
          do i = jstart, jend
             iold = perm% iperm(i)
             ipos = ia(iold)
             y(iold) = a(ipos)*x(ja(ipos))
          enddo
       else
          !----------------------------------------------------------
          ! Work through current group in chunks of NB rows at a time
          !----------------------------------------------------------
          do istart = jstart, jend, NB
             iend  = min (jend, istart+(NB-1))
             isize = iend-istart+1
             !-------------------------------------------------------------
             ! ix(i) points to the beginning of the ith row of a chunk in a
             !-------------------------------------------------------------
             ix(1:isize) = ia(perm% iperm(istart:iend))
             yp(1:isize) = 0.0_wp
             if (nz > isize) then
                !----------------------------------------------------------
                ! The number of nonzeros per row exceeds the number of rows
                ! in the chunk.  We should vectorize along nz, i.e. one row
                ! at a time, as in the usual CSR case.
                !----------------------------------------------------------
                do i = 1, isize
                   do j = 1, nz
                      ipos = ix(i)+(j-1)
                      yp(i) = yp(i)+a(ipos)*x(ja(ipos))
                   enddo
                enddo
             else
                !----------------------------------------------------------
                ! There are enough rows in the chunk to vectorize across
                ! the rows, i.e. to operate with "columns" of the chunk.
                !----------------------------------------------------------
                do j = 1, nz
                   do i = 1, isize
                      ipos = ix(i)+(j-1)
                      yp(i) = yp(i)+a(ipos)*x(ja(ipos))
                   enddo
                enddo
             endif
             y(perm% iperm(istart:iend)) = yp(1:isize)
          enddo
       endif
    enddo

  end subroutine csrperm_times_vector
  !============================================================================
  subroutine send_perm (perm, dest, tag, comm)
    type(t_perm),      intent(in) :: perm
    integer,           intent(in) :: dest
    integer, optional, intent(in) :: tag
    integer, optional, intent(in) :: comm

    integer :: ltag, lcom
    integer :: count = 0

    if (count == 0) count = size (transfer (perm, (/' '/)))

    ltag = 2          ;if (present (tag))  ltag = tag
    lcom = dace% comm ;if (present (comm)) lcom = comm

    call p_send_derivedtype2 (perm, count, dest, ltag, lcom)

    if (associated (perm% iperm)  ) &
         call p_send (perm% iperm,   dest, ltag, comm=lcom)
    if (associated (perm% ia)     ) &
         call p_send (perm% ia,      dest, ltag, comm=lcom)
    if (associated (perm% nzgroup)) &
         call p_send (perm% nzgroup, dest, ltag, comm=lcom)
    if (associated (perm% xgroup) ) &
         call p_send (perm% xgroup,  dest, ltag, comm=lcom)

  end subroutine send_perm
  !----------------------------------------------------------------------------
  subroutine recv_perm (perm, source, tag, comm)
    type(t_perm),      intent(inout) :: perm
    integer,           intent(in)    :: source
    integer, optional, intent(in)    :: tag
    integer, optional, intent(in)    :: comm

    integer :: ltag, lcom
    integer :: count = 0

    if (count == 0) count = size (transfer (perm, (/' '/)))

    ltag = 2          ;if (present (tag))  ltag = tag
    lcom = dace% comm ;if (present (comm)) lcom = comm

    call deallocate (perm)
    call p_recv_derivedtype2 (perm, count, source, ltag, lcom)
    call allocate_perm (perm, dealloc=.false.)

    if (associated (perm% iperm)  ) then
       call p_recv (perm% iperm,   source, ltag, comm=lcom)
    end if
    if (associated (perm% ia)     ) then
       call p_recv (perm% ia,      source, ltag, comm=lcom)
    end if
    if (associated (perm% nzgroup)) then
       call p_recv (perm% nzgroup, source, ltag, comm=lcom)
    end if
    if (associated (perm% xgroup) ) then
       call p_recv (perm% xgroup,  source, ltag, comm=lcom)
    end if

  end subroutine recv_perm
  !----------------------------------------------------------------------------
  subroutine bcast_perm (perm, source, comm)
    type(t_perm),      intent(inout) :: perm
    integer,           intent(in)    :: source
    integer, optional, intent(in)    :: comm

    integer :: lcom
    integer :: count = 0

    if (count == 0) count = size (transfer (perm, (/' '/)))

    lcom = dace% comm ;if (present (comm)) lcom = comm

    if (dace% pe /= source) call deallocate (perm)
    call p_bcast_derivedtype2 (perm, count, source, lcom)
    if (dace% pe /= source) call allocate_perm (perm, dealloc=.false.)

    if (associated (perm% iperm)  ) call p_bcast (perm% iperm,   source, lcom)
    if (associated (perm% ia)     ) call p_bcast (perm% ia,      source, lcom)
    if (associated (perm% nzgroup)) call p_bcast (perm% nzgroup, source, lcom)
    if (associated (perm% xgroup) ) call p_bcast (perm% xgroup,  source, lcom)

  end subroutine bcast_perm
!==============================================================================
  subroutine convert_matrix_block (b, repr)
    type(t_matrix_block), intent(inout) :: b
    integer,              intent(in)    :: repr
    !------------------------------------------------------------
    ! Convert a matrix from one sparse representation to another.
    ! Currently supported representations: CSR, CSC, JAD
    ! Restrictions: B must be square.
    !------------------------------------------------------------
    integer               :: m, nz, jdiag
    real(mp), allocatable :: a_tmp(:)
    integer,  allocatable :: ia_tmp(:), ja_tmp(:)

    if (repr < 0 .or. repr > ubound (srep,1)) then
       write (0,*) "Conversion out of range: repr =", repr
       call finish ("convert_matrix_block","cannot convert")
    end if
    if (b% repr < 0 .or. b% repr > ubound (srep,1)) then
       write (0,*) "Bad representation of input matrix: repr =", b% repr
       call finish ("convert_matrix_block","bad matrix")
    end if
    if (b% m /= b% n) then
       print *, "Non-square matrix:", b% m, "*", b% n
       call finish ("convert_matrix_block","unsupported: non-square matrix")
    end if

    if (repr == b% repr) return
    !--------------------------------------------------------
    ! Convert from input representation to (intermediate) CSR
    !--------------------------------------------------------
    m  = b% m
    nz = b% nonzero
    select case (b% repr)
    case (CSR)
       ! We're already set
    case (CSC)
       allocate (a_tmp(nz), ja_tmp(nz), ia_tmp(m+1))
       call csrcsc (m,1,1, b% packed, b% ja, b% ia, a_tmp, ja_tmp, ia_tmp)
       b% packed = a_tmp
       b% ja     = ja_tmp
       b% ia     = ia_tmp
       deallocate (a_tmp, ja_tmp, ia_tmp)
    case (JAD)
       allocate (b% ia(m+1))
       allocate (a_tmp(nz), ja_tmp(nz))
       call jadcsr (m, b%perm% jdiag, b% packed, b% ja, b%perm% ia, &
                       b%perm% iperm, a_tmp, ja_tmp, b% ia)
       b% packed = a_tmp
       b% ja     = ja_tmp
       deallocate (a_tmp, ja_tmp)
    case default
       write (0,*) "Input: srep = ", trim (srep (b% repr))
       call finish ("convert_matrix_block", &
            "cannot convert from" // trim (srep (b% repr)))
    end select
    b% repr = CSR               ! Intermediate representation is plain CSR
    call destruct_perm (b)
    !-------------------------------------------------------
    ! Convert from intermediate CSR to output representation
    !-------------------------------------------------------
    select case (repr)
    case (CSR)
       ! We're already set
    case (CSC)
       allocate (a_tmp(nz), ja_tmp(nz), ia_tmp(m+1))
       call csrcsc (m,1,1, b% packed, b% ja, b% ia, a_tmp, ja_tmp, ia_tmp)
       b% packed = a_tmp
       b% ja     = ja_tmp
       b% ia     = ia_tmp
       deallocate (a_tmp, ja_tmp, ia_tmp)
       b% repr = CSC
    case (JAD)
       allocate (b% perm)
       b%perm% m = m
       allocate (b%perm% iperm(m))          ! Row permutation
       allocate (a_tmp(nz), ja_tmp(nz), ia_tmp(m+1))
       call csrjad (m, b% packed, b% ja, b% ia, &
                       jdiag, b%perm% iperm, a_tmp, ja_tmp, ia_tmp)
       b%perm% jdiag = jdiag
       deallocate (b% ia)
       nullify (b% ia)
       allocate (b%perm% ia(jdiag+1))
       b% packed  = a_tmp
       b% ja      = ja_tmp
       b%perm% ia = ia_tmp(1:jdiag+1)
       deallocate (a_tmp, ja_tmp, ia_tmp)
       b% repr = JAD
!write(0,*) "pe=", int(dace% pe,1), &
!     "convert_matrix_block[JAD]:", size (b%perm% iperm),size (b%perm% ia)

    case default
       write (0,*) "Required conversion: srep = ", trim (srep (repr))
       call finish ("convert_matrix_block", &
            "cannot convert to" // trim (srep (repr)))
    end select
  end subroutine convert_matrix_block
!==============================================================================
  subroutine print_crc_vector(v, hint)
    !------------------------------------------------------------
    ! print crc checksum of vector, useful for debugging
    !------------------------------------------------------------
    type(t_vector),   intent(in)           :: v
    character(len=*), intent(in), optional :: hint
    character(len=120) :: msg = ''
    integer            :: i
    msg = 'crc_vector '//trim(v%name)
    if (present(hint)) msg=trim(msg)//' '//trim(hint)
    do i = 1, v%n_s
      if (dace%pe == v%s(i)%pe) then
        write(*,*) trim(msg),i,crc(v%s(i)%x(1:v%s(i)%n))
      end if
    end do
  end subroutine print_crc_vector
!==============================================================================
  subroutine print_crc_matrix(m, hint)
    !------------------------------------------------------------
    ! print crc checksum of matrix, useful for debugging
    !------------------------------------------------------------
    type(t_matrix),   intent(in), target   :: m
    character(len=*), intent(in), optional :: hint
    type(t_matrix_block), pointer :: b
    character(len=120) :: msg = ''
    integer            :: i,j
    msg = 'crc_matrix'
    if (present(hint)) msg=trim(msg)//' '//trim(hint)
    do i = 1, size(m%b,1)
      do j = 1, size(m%b,2)
        b => m%b(i,j)
        if (b%pe == dace%pe) then
          if (associated(b%full)) then
            write(*,*) trim(msg),' full',shape(b%full),crc(reshape(b%full,(/size(b%full)/)))
          elseif (associated(b%packed)) then
            write(*,*) trim(msg),' packed',shape(b%packed),crc(b%packed)
          end if
        end if
      end do
    end do
  end subroutine print_crc_matrix

end module mo_dec_matrix
